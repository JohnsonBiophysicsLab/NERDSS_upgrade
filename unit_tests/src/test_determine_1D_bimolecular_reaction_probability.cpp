/*! \file test_determine_1D_bimolecular_reaction_probability.cpp
 *
 * ### Unit test for src/reactions/determine_1D_bimolecular_reaction_probability.cpp
 *
 * Function under test:
 * \code
 *   void determine_1D_bimolecular_reaction_probability(
 *       int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
 *       BiMolData& biMolData, const Parameters& params,
 *       std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
 *       const std::vector<ForwardRxn>& forwardRxns,
 *       const std::vector<BackRxn>& backRxns);
 * \endcode
 *
 * The routine evaluates whether two interfaces that live on a one-dimensional
 * fiber are close enough to react during the current time step and, if so,
 * records the encounter (crossbase / mycrossint / crossrxn / probvec / ncross)
 * and computes the (possibly reweighted) association probability.
 *
 * The behaviour that is verified here:
 *   1. Encounter bookkeeping when the pair is inside both the 1D (x-axis) and
 *      the 3D "same fiber" cutoffs.
 *   2. No bookkeeping at all when the pair is beyond the 1D cutoff.
 *   3. No bookkeeping when the x-separation is small but the molecules are on
 *      *different* fibers (large 3D separation).
 *   4. A zero reaction rate leaves the stored probability at exactly 0.
 *   5. A molecule that just dissociated has its probability evaluation skipped
 *      (encounter is still recorded, probability stays 0).
 *   6. Promoter (3D -> 1D binding) pairs use a smaller RMax because sigma_x = 0.
 *   7. Separations smaller than the binding radius are clamped up to the
 *      binding radius before the probability is evaluated.
 *   8. Reweighting history: a previous separation >= RMax resets the norm to
 *      unity (probability identical to a fresh pair), while a previous
 *      separation inside RMax applies a finite, positive reweighting factor.
 *
 * All console output goes to stderr so the reader can follow exactly which
 * scenario is being exercised and which criteria are used to pass.
 */

#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local test fixture helpers (unique "d1dbrp_" prefix to avoid collisions with
// the other translation units that get linked into the same test binary).
// -----------------------------------------------------------------------------
namespace {

//! Binding radius shared by every scenario below (nm).
constexpr double kD1dbrpBindRadius = 1.0;
//! Total diffusion constant used for the pair (nm^2/us).
constexpr double kD1dbrpDtot = 1.0;
//! Time step used for the pair (us).
constexpr double kD1dbrpTimeStep = 1.0;

/*! \brief Everything the function under test needs, bundled for convenience. */
struct D1dbrpScenario {
    Parameters params;                    //!< simulation parameters (time step)
    std::vector<Molecule> moleculeList;   //!< exactly two molecules
    std::vector<Complex> complexList;     //!< one complex per molecule
    std::vector<ForwardRxn> forwardRxns;  //!< a single 1D forward reaction
    std::vector<BackRxn> backRxns;        //!< unused by the function, kept empty
    BiMolData biMolData;                  //!< pair data handed to the function
};

/*! \brief Build a minimal molecule with a single interface at \p coord.
 *
 * The interface coordinate is what the function actually measures, so the COM
 * is simply placed at the same point to keep the geometry trivial.
 */
Molecule d1dbrp_make_molecule(int index, int comIndex, const Coord& coord)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.comCoord = coord;

    Molecule::Iface iface;
    iface.coord = coord;
    iface.index = 0;
    iface.relIndex = 0;
    iface.molTypeIndex = 0;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a one-member complex that owns molecule \p molIndex. */
Complex d1dbrp_make_complex(int index, int molIndex, const Coord& coord)
{
    Complex com;
    com.index = index;
    com.comCoord = coord;
    com.memberList.clear();
    com.memberList.push_back(molIndex);
    com.D = Coord(kD1dbrpDtot / 2.0, kD1dbrpDtot / 2.0, kD1dbrpDtot / 2.0);
    com.Dr = Coord(0.0, 0.0, 0.0);
    com.ncross = 0;
    return com;
}

/*! \brief Assemble a two-molecule scenario separated by (dx, dy, dz).
 *
 * \param[in] dx        x-axis separation of the two reacting interfaces (the
 *                      coordinate that matters for the 1D cutoff).
 * \param[in] dy,dz     transverse separation, used only by the 3D "same fiber"
 *                      test inside the function.
 * \param[in] rate      macroscopic forward rate of the single reaction.
 */
D1dbrpScenario d1dbrp_make_scenario(double dx, double dy, double dz, double rate)
{
    D1dbrpScenario s;

    // --- parameters -----------------------------------------------------
    s.params.timeStep = kD1dbrpTimeStep;

    // --- molecules: molecule 0 at the origin, molecule 1 offset ---------
    s.moleculeList.push_back(d1dbrp_make_molecule(0, 0, Coord(0.0, 0.0, 0.0)));
    s.moleculeList.push_back(d1dbrp_make_molecule(1, 1, Coord(dx, dy, dz)));

    // --- one complex per molecule so that ncross can be tracked apart ---
    s.complexList.push_back(d1dbrp_make_complex(0, 0, Coord(0.0, 0.0, 0.0)));
    s.complexList.push_back(d1dbrp_make_complex(1, 1, Coord(dx, dy, dz)));

    // --- a single 1D forward reaction -----------------------------------
    ForwardRxn rxn;
    rxn.bindRadius = kD1dbrpBindRadius;
    // The 3D->1D conversion area; kact = rate / area3Dto1D inside the function.
    rxn.area3Dto1D = 4.0 * M_PI * kD1dbrpBindRadius * kD1dbrpBindRadius;
    rxn.isSymmetric = false;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    RxnBase::RateState rateState;
    rateState.rate = rate;
    rxn.rateList.push_back(rateState);
    s.forwardRxns.push_back(rxn);

    // --- pair data ------------------------------------------------------
    s.biMolData = BiMolData(/*pro1Index*/ 0, /*pro2Index*/ 1,
        /*com1Index*/ 0, /*com2Index*/ 1,
        /*relIface1*/ 0, /*relIface2*/ 0,
        /*absIface1*/ 0, /*absIface2*/ 0,
        /*Dtot*/ kD1dbrpDtot,
        /*magMol1*/ 1.0, /*magMol2*/ 1.0);

    return s;
}

//! Analytic 1D cutoff used by the function for a normal (non-promoter) pair.
double d1dbrp_rmax_normal()
{
    return 4.0 * std::sqrt(2.0 * kD1dbrpDtot * kD1dbrpTimeStep) + kD1dbrpBindRadius;
}

//! Analytic 1D cutoff used by the function for a promoter (3D->1D) pair.
double d1dbrp_rmax_promoter()
{
    return 4.0 * std::sqrt(2.0 * kD1dbrpDtot * kD1dbrpTimeStep);
}

//! Analytic 3D "same fiber" cutoff used by the function.
double d1dbrp_rmax_3d()
{
    return 3.0 * std::sqrt(6.0 * kD1dbrpDtot * kD1dbrpTimeStep) + kD1dbrpBindRadius;
}

/*! \brief Convenience wrapper that calls the function under test. */
void d1dbrp_invoke(D1dbrpScenario& s)
{
    determine_1D_bimolecular_reaction_probability(
        /*simItr*/ 1, /*rxnIndex*/ 0, /*rateIndex*/ 0,
        /*isStateChangeBackRxn*/ false, s.biMolData, s.params,
        s.moleculeList, s.complexList, s.forwardRxns, s.backRxns);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: pair inside both cutoffs -> full encounter bookkeeping + probability.
// -----------------------------------------------------------------------------
void test_d1dbrp_within_rmax_records_encounter()
{
    std::cerr << "\n[TEST] test_d1dbrp_within_rmax_records_encounter\n"
              << "  Source file: determine_1D_bimolecular_reaction_probability.cpp\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    dx = 1.5 nm (inside RMax1D = " << d1dbrp_rmax_normal()
              << "), dy = dz = 0 (same fiber, RMax3D = " << d1dbrp_rmax_3d() << ").\n"
              << "  Criteria:    crossbase/mycrossint/crossrxn/probvec each gain one\n"
              << "               entry on BOTH molecules, both complexes' ncross grow\n"
              << "               by 1, the stored probability is positive and <= 1,\n"
              << "               and the reweighting history is appended on the\n"
              << "               lower-index molecule (proA == 0).\n";

    D1dbrpScenario s = d1dbrp_make_scenario(/*dx*/ 1.5, /*dy*/ 0.0, /*dz*/ 0.0, /*rate*/ 10.0);
    d1dbrp_invoke(s);

    // --- encounter lists on molecule 0 ---------------------------------
    EXPECT_EQ(s.moleculeList[0].crossbase.size(), 1u)
        << "molecule 0 should have recorded exactly one crossing partner";
    EXPECT_EQ(s.moleculeList[1].crossbase.size(), 1u)
        << "molecule 1 should have recorded exactly one crossing partner";
    if (!s.moleculeList[0].crossbase.empty())
        EXPECT_EQ(s.moleculeList[0].crossbase[0], 1) << "partner of molecule 0 must be molecule 1";
    if (!s.moleculeList[1].crossbase.empty())
        EXPECT_EQ(s.moleculeList[1].crossbase[0], 0) << "partner of molecule 1 must be molecule 0";

    // --- interface indices --------------------------------------------
    EXPECT_EQ(s.moleculeList[0].mycrossint.size(), 1u) << "one crossing interface on molecule 0";
    EXPECT_EQ(s.moleculeList[1].mycrossint.size(), 1u) << "one crossing interface on molecule 1";
    if (!s.moleculeList[0].mycrossint.empty())
        EXPECT_EQ(s.moleculeList[0].mycrossint[0], 0) << "relIface1 should be stored";
    if (!s.moleculeList[1].mycrossint.empty())
        EXPECT_EQ(s.moleculeList[1].mycrossint[0], 0) << "relIface2 should be stored";

    // --- reaction identifiers -----------------------------------------
    EXPECT_EQ(s.moleculeList[0].crossrxn.size(), 1u) << "one crossrxn record on molecule 0";
    if (!s.moleculeList[0].crossrxn.empty()) {
        EXPECT_EQ(s.moleculeList[0].crossrxn[0][0], 0) << "rxnIndex should be 0";
        EXPECT_EQ(s.moleculeList[0].crossrxn[0][1], 0) << "rateIndex should be 0";
        EXPECT_EQ(s.moleculeList[0].crossrxn[0][2], 0) << "isStateChangeBackRxn should be false(0)";
    }

    // --- ncross on each parent complex --------------------------------
    EXPECT_EQ(s.complexList[0].ncross, 1) << "complex 0 ncross should be incremented once";
    EXPECT_EQ(s.complexList[1].ncross, 1) << "complex 1 ncross should be incremented once";

    // --- probability --------------------------------------------------
    ASSERT_EQ(s.moleculeList[0].probvec.size(), 1u)
        << "molecule 0 should hold exactly one probability entry";
    ASSERT_EQ(s.moleculeList[1].probvec.size(), 1u)
        << "molecule 1 should hold exactly one probability entry";
    const double prob = s.moleculeList[0].probvec.back();
    std::cerr << "  Computed association probability = " << prob << '\n';
    EXPECT_TRUE(std::isfinite(prob)) << "probability must be a finite number";
    EXPECT_GT(prob, 0.0) << "a contact-range pair with a positive rate must have p > 0";
    EXPECT_LE(prob, 1.000001) << "probability must not exceed unity";
    EXPECT_DOUBLE_EQ(s.moleculeList[1].probvec.back(), prob)
        << "both partners must store the same probability";

    // --- reweighting history is stored on proA (the lower index) -------
    ASSERT_EQ(s.moleculeList[0].currprevsep.size(), 1u)
        << "proA (molecule 0) should store the current separation";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].currprevsep.back(), 1.5)
        << "stored separation should be the x-axis distance (1.5 nm)";
    EXPECT_EQ(s.moleculeList[0].currlist.size(), 1u) << "currlist gains the partner index";
    if (!s.moleculeList[0].currlist.empty())
        EXPECT_EQ(s.moleculeList[0].currlist.back(), 1) << "partner should be molecule 1";
    EXPECT_EQ(s.moleculeList[0].currmyface.size(), 1u) << "currmyface gains one entry";
    EXPECT_EQ(s.moleculeList[0].currpface.size(), 1u) << "currpface gains one entry";
    ASSERT_EQ(s.moleculeList[0].currprevnorm.size(), 1u) << "currprevnorm gains one entry";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].currprevnorm.back(), 1.0)
        << "with no previous encounter the norm must stay at unity";
    ASSERT_EQ(s.moleculeList[0].currps_prev.size(), 1u) << "currps_prev gains one entry";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].currps_prev.back(), 1.0 - prob)
        << "currps_prev must equal 1 - p";

    // Molecule 1 is proB: it must not accumulate reweighting history.
    EXPECT_TRUE(s.moleculeList[1].currprevsep.empty())
        << "proB should not store reweighting history";
}

// -----------------------------------------------------------------------------
// Test 2: pair beyond the 1D cutoff -> nothing at all is recorded.
// -----------------------------------------------------------------------------
void test_d1dbrp_outside_rmax_no_encounter()
{
    std::cerr << "\n[TEST] test_d1dbrp_outside_rmax_no_encounter\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    dx = 40 nm, far beyond RMax1D = " << d1dbrp_rmax_normal()
              << " and RMax3D = " << d1dbrp_rmax_3d() << ".\n"
              << "  Criteria:    all encounter containers stay empty and ncross stays 0.\n";

    D1dbrpScenario s = d1dbrp_make_scenario(/*dx*/ 40.0, /*dy*/ 0.0, /*dz*/ 0.0, /*rate*/ 10.0);
    d1dbrp_invoke(s);

    EXPECT_TRUE(s.moleculeList[0].crossbase.empty()) << "no crossing should be recorded (mol 0)";
    EXPECT_TRUE(s.moleculeList[1].crossbase.empty()) << "no crossing should be recorded (mol 1)";
    EXPECT_TRUE(s.moleculeList[0].probvec.empty()) << "no probability entry should be created";
    EXPECT_TRUE(s.moleculeList[1].probvec.empty()) << "no probability entry should be created";
    EXPECT_TRUE(s.moleculeList[0].currlist.empty()) << "no reweighting history should be created";
    EXPECT_EQ(s.complexList[0].ncross, 0) << "complex 0 ncross must remain 0";
    EXPECT_EQ(s.complexList[1].ncross, 0) << "complex 1 ncross must remain 0";
}

// -----------------------------------------------------------------------------
// Test 3: small x-separation but large 3D separation -> different fibers, so the
//         pair must be rejected even though the 1D test alone would pass.
// -----------------------------------------------------------------------------
void test_d1dbrp_different_fiber_rejected()
{
    std::cerr << "\n[TEST] test_d1dbrp_different_fiber_rejected\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    dx = 1.0 nm (inside RMax1D = " << d1dbrp_rmax_normal()
              << ") but dy = 50 nm so the 3D distance (" << std::sqrt(1.0 + 2500.0)
              << ") exceeds RMax3D = " << d1dbrp_rmax_3d() << ".\n"
              << "  Criteria:    the pair is NOT treated as an encounter (different fiber).\n";

    D1dbrpScenario s = d1dbrp_make_scenario(/*dx*/ 1.0, /*dy*/ 50.0, /*dz*/ 0.0, /*rate*/ 10.0);
    d1dbrp_invoke(s);

    EXPECT_TRUE(s.moleculeList[0].crossbase.empty())
        << "molecules on different fibers must not be paired";
    EXPECT_TRUE(s.moleculeList[0].probvec.empty())
        << "no probability entry for molecules on different fibers";
    EXPECT_EQ(s.complexList[0].ncross, 0) << "ncross must not be incremented";
    EXPECT_EQ(s.complexList[1].ncross, 0) << "ncross must not be incremented";
}

// -----------------------------------------------------------------------------
// Test 4: zero reaction rate -> encounter recorded, probability left at zero.
// -----------------------------------------------------------------------------
void test_d1dbrp_zero_rate_keeps_probability_zero()
{
    std::cerr << "\n[TEST] test_d1dbrp_zero_rate_keeps_probability_zero\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    dx = 1.5 nm (inside RMax) but the reaction rate is 0.\n"
              << "  Criteria:    the encounter is still recorded, but probvec stays at 0\n"
              << "               and no reweighting history is stored.\n";

    D1dbrpScenario s = d1dbrp_make_scenario(/*dx*/ 1.5, /*dy*/ 0.0, /*dz*/ 0.0, /*rate*/ 0.0);
    d1dbrp_invoke(s);

    ASSERT_EQ(s.moleculeList[0].probvec.size(), 1u)
        << "the encounter itself should still push a probability slot";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].probvec.back(), 0.0)
        << "a zero rate must leave the probability at exactly 0";
    EXPECT_DOUBLE_EQ(s.moleculeList[1].probvec.back(), 0.0)
        << "a zero rate must leave the partner probability at exactly 0";
    EXPECT_TRUE(s.moleculeList[0].currprevsep.empty())
        << "no reweighting history should be created when the rate is 0";
    EXPECT_EQ(s.complexList[0].ncross, 1) << "the encounter still increments ncross";
}

// -----------------------------------------------------------------------------
// Test 5: a molecule that just dissociated -> probability evaluation skipped.
// -----------------------------------------------------------------------------
void test_d1dbrp_dissociated_skips_probability()
{
    std::cerr << "\n[TEST] test_d1dbrp_dissociated_skips_probability\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    dx = 1.5 nm (inside RMax) but molecule 0 has\n"
              << "               isDissociated == true.\n"
              << "  Criteria:    the encounter is recorded (overlap avoidance) but the\n"
              << "               probability stays 0 and no history is stored.\n";

    D1dbrpScenario s = d1dbrp_make_scenario(/*dx*/ 1.5, /*dy*/ 0.0, /*dz*/ 0.0, /*rate*/ 10.0);
    s.moleculeList[0].isDissociated = true;
    d1dbrp_invoke(s);

    EXPECT_EQ(s.moleculeList[0].crossbase.size(), 1u)
        << "the crossing is still needed so overlap can be resolved";
    ASSERT_EQ(s.moleculeList[0].probvec.size(), 1u) << "a probability slot is still pushed";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].probvec.back(), 0.0)
        << "a just-dissociated molecule must not be given a reaction probability";
    EXPECT_TRUE(s.moleculeList[0].currlist.empty())
        << "no reweighting bookkeeping for a just-dissociated molecule";
}

// -----------------------------------------------------------------------------
// Test 6: promoter pairs (3D -> 1D binding, sigma_x = 0) use a smaller RMax.
// -----------------------------------------------------------------------------
void test_d1dbrp_promoter_uses_smaller_rmax()
{
    std::cerr << "\n[TEST] test_d1dbrp_promoter_uses_smaller_rmax\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    dx = 6.0 nm, which lies between the promoter cutoff ("
              << d1dbrp_rmax_promoter() << ") and the normal cutoff ("
              << d1dbrp_rmax_normal() << ").\n"
              << "  Criteria:    a non-promoter pair reacts (encounter recorded) while an\n"
              << "               otherwise identical promoter/non-promoter pair does not,\n"
              << "               because sigma_x is dropped from RMax for 3D->1D binding.\n";

    // (a) Neither molecule is a promoter -> the larger cutoff applies.
    D1dbrpScenario plain = d1dbrp_make_scenario(/*dx*/ 6.0, 0.0, 0.0, /*rate*/ 10.0);
    d1dbrp_invoke(plain);
    EXPECT_EQ(plain.moleculeList[0].crossbase.size(), 1u)
        << "dx = 6.0 is inside the normal RMax, so the encounter must be recorded";

    // (b) Molecule 0 is a promoter -> the smaller cutoff applies.
    D1dbrpScenario promoter = d1dbrp_make_scenario(/*dx*/ 6.0, 0.0, 0.0, /*rate*/ 10.0);
    promoter.moleculeList[0].isPromoter = true;
    d1dbrp_invoke(promoter);
    EXPECT_TRUE(promoter.moleculeList[0].crossbase.empty())
        << "dx = 6.0 is outside the promoter RMax, so no encounter should be recorded";
    EXPECT_EQ(promoter.complexList[0].ncross, 0)
        << "no ncross increment for a promoter pair beyond its RMax";

    // (c) The reverse promoter orientation must behave identically.
    D1dbrpScenario promoter2 = d1dbrp_make_scenario(/*dx*/ 6.0, 0.0, 0.0, /*rate*/ 10.0);
    promoter2.moleculeList[1].isPromoter = true;
    d1dbrp_invoke(promoter2);
    EXPECT_TRUE(promoter2.moleculeList[1].crossbase.empty())
        << "promoter check must be symmetric in the two reactants";
}

// -----------------------------------------------------------------------------
// Test 7: separations below the binding radius are clamped up to it (unless the
//         pair binds from 3D to 1D, where sigma_x = 0 and no clamp is applied).
// -----------------------------------------------------------------------------
void test_d1dbrp_clamps_separation_to_bind_radius()
{
    std::cerr << "\n[TEST] test_d1dbrp_clamps_separation_to_bind_radius\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    dx = 0.4 nm, below the binding radius of "
              << kD1dbrpBindRadius << " nm.\n"
              << "  Criteria:    for a normal pair the stored separation is clamped to\n"
              << "               the binding radius; for a promoter pair (sigma_x = 0)\n"
              << "               the raw separation is kept.\n";

    // (a) Normal pair -> R1 is clamped up to bindRadius before being stored.
    D1dbrpScenario plain = d1dbrp_make_scenario(/*dx*/ 0.4, 0.0, 0.0, /*rate*/ 10.0);
    d1dbrp_invoke(plain);
    ASSERT_EQ(plain.moleculeList[0].currprevsep.size(), 1u)
        << "an inside-RMax pair with a positive rate must store its separation";
    EXPECT_DOUBLE_EQ(plain.moleculeList[0].currprevsep.back(), kD1dbrpBindRadius)
        << "separation below sigma must be clamped up to sigma";
    std::cerr << "  Normal pair: stored sep = " << plain.moleculeList[0].currprevsep.back()
              << ", probability = " << plain.moleculeList[0].probvec.back() << '\n';

    // (b) Promoter pair -> no clamp, the raw x-distance is used.
    D1dbrpScenario promoter = d1dbrp_make_scenario(/*dx*/ 0.4, 0.0, 0.0, /*rate*/ 10.0);
    promoter.moleculeList[0].isPromoter = true;
    d1dbrp_invoke(promoter);
    ASSERT_EQ(promoter.moleculeList[0].currprevsep.size(), 1u)
        << "the promoter pair is well inside RMax and must be evaluated";
    EXPECT_DOUBLE_EQ(promoter.moleculeList[0].currprevsep.back(), 0.4)
        << "3D->1D binding uses sigma_x = 0, so the raw separation is kept";
    std::cerr << "  Promoter pair: stored sep = " << promoter.moleculeList[0].currprevsep.back()
              << ", probability = " << promoter.moleculeList[0].probvec.back() << '\n';
    EXPECT_TRUE(std::isfinite(promoter.moleculeList[0].probvec.back()))
        << "the promoter probability must be finite";
}

// -----------------------------------------------------------------------------
// Test 8: reweighting history with a previous separation >= RMax resets the
//         normalization to unity, i.e. the probability equals a fresh pair's.
// -----------------------------------------------------------------------------
void test_d1dbrp_far_previous_separation_resets_norm()
{
    std::cerr << "\n[TEST] test_d1dbrp_far_previous_separation_resets_norm\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    the pair has a stored history whose previous separation\n"
              << "               (100 nm) is >= RMax, which forces currnorm back to 1.\n"
              << "  Criteria:    the probability equals the probability of an identical\n"
              << "               pair with no history, and currprevnorm == 1.\n";

    // Baseline: identical geometry, no reweighting history at all.
    D1dbrpScenario baseline = d1dbrp_make_scenario(/*dx*/ 1.5, 0.0, 0.0, /*rate*/ 10.0);
    d1dbrp_invoke(baseline);
    ASSERT_EQ(baseline.moleculeList[0].probvec.size(), 1u) << "baseline must be evaluated";
    const double baseProb = baseline.moleculeList[0].probvec.back();

    // Same geometry, but with a "far" previous encounter recorded on proA.
    D1dbrpScenario withHistory = d1dbrp_make_scenario(/*dx*/ 1.5, 0.0, 0.0, /*rate*/ 10.0);
    withHistory.moleculeList[0].prevlist.push_back(1);    // partner index
    withHistory.moleculeList[0].prevmyface.push_back(0);  // my interface
    withHistory.moleculeList[0].prevpface.push_back(0);   // partner interface
    withHistory.moleculeList[0].prevsep.push_back(100.0); // >= RMax -> reset
    withHistory.moleculeList[0].prevnorm.push_back(0.25); // should be ignored
    withHistory.moleculeList[0].ps_prev.push_back(0.5);   // should be ignored
    d1dbrp_invoke(withHistory);

    ASSERT_EQ(withHistory.moleculeList[0].probvec.size(), 1u)
        << "the pair with history must still be evaluated";
    std::cerr << "  baseline p = " << baseProb
              << ", with far history p = " << withHistory.moleculeList[0].probvec.back() << '\n';
    EXPECT_DOUBLE_EQ(withHistory.moleculeList[0].probvec.back(), baseProb)
        << "a previous separation beyond RMax must restart reweighting (norm = 1)";
    ASSERT_EQ(withHistory.moleculeList[0].currprevnorm.size(), 1u)
        << "the new norm must be recorded";
    EXPECT_DOUBLE_EQ(withHistory.moleculeList[0].currprevnorm.back(), 1.0)
        << "currnorm must be reset to unity";
}

// -----------------------------------------------------------------------------
// Test 9: reweighting history with a previous separation inside RMax applies a
//         finite, positive reweighting factor to the raw probability.
// -----------------------------------------------------------------------------
void test_d1dbrp_close_previous_separation_applies_norm()
{
    std::cerr << "\n[TEST] test_d1dbrp_close_previous_separation_applies_norm\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    the pair has a stored history whose previous separation\n"
              << "               (1.6 nm) is inside RMax, with prevnorm = 1 and\n"
              << "               ps_prev = 0.5, so pirr/pfree reweighting is applied.\n"
              << "  Criteria:    currprevnorm is finite and > 0, the stored probability\n"
              << "               equals rawProbability * currnorm (checked through the\n"
              << "               currps_prev == 1 - p relation) and stays <= 1.\n";

    D1dbrpScenario s = d1dbrp_make_scenario(/*dx*/ 1.5, 0.0, 0.0, /*rate*/ 10.0);
    s.moleculeList[0].prevlist.push_back(1);
    s.moleculeList[0].prevmyface.push_back(0);
    s.moleculeList[0].prevpface.push_back(0);
    s.moleculeList[0].prevsep.push_back(1.6);  // inside RMax -> reweight
    s.moleculeList[0].prevnorm.push_back(1.0);
    s.moleculeList[0].ps_prev.push_back(0.5);
    d1dbrp_invoke(s);

    ASSERT_EQ(s.moleculeList[0].probvec.size(), 1u) << "the pair must be evaluated";
    ASSERT_EQ(s.moleculeList[0].currprevnorm.size(), 1u) << "a new norm must be recorded";
    const double norm = s.moleculeList[0].currprevnorm.back();
    const double prob = s.moleculeList[0].probvec.back();
    std::cerr << "  reweighted norm = " << norm << ", probability = " << prob << '\n';

    EXPECT_TRUE(std::isfinite(norm)) << "the reweighting factor must be finite";
    EXPECT_GT(norm, 0.0) << "the reweighting factor must be positive";
    EXPECT_TRUE(std::isfinite(prob)) << "the reweighted probability must be finite";
    EXPECT_GE(prob, 0.0) << "the reweighted probability must be non-negative";
    EXPECT_LE(prob, 1.000001) << "the reweighted probability must not exceed unity";
    EXPECT_DOUBLE_EQ(s.moleculeList[1].probvec.back(), prob)
        << "both partners store the same reweighted probability";
    ASSERT_EQ(s.moleculeList[0].currps_prev.size(), 1u) << "currps_prev must be recorded";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].currps_prev.back(), 1.0 - prob)
        << "currps_prev must be 1 - p for the next step's reweighting";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].currprevsep.back(), 1.5)
        << "the current separation (1.5 nm) must be stored for the next step";
}

// -----------------------------------------------------------------------------
// Test 10: the reweighting history is matched on (partner, myface, pface), so a
//          history entry for a different partner must be ignored.
// -----------------------------------------------------------------------------
void test_d1dbrp_history_for_other_partner_ignored()
{
    std::cerr << "\n[TEST] test_d1dbrp_history_for_other_partner_ignored\n"
              << "  Function:    determine_1D_bimolecular_reaction_probability\n"
              << "  Scenario:    proA carries a history entry that refers to a different\n"
              << "               partner index (7), so it must not be applied.\n"
              << "  Criteria:    currprevnorm == 1 and the probability matches the\n"
              << "               no-history baseline.\n";

    // Baseline without any history.
    D1dbrpScenario baseline = d1dbrp_make_scenario(/*dx*/ 1.5, 0.0, 0.0, /*rate*/ 10.0);
    d1dbrp_invoke(baseline);
    ASSERT_EQ(baseline.moleculeList[0].probvec.size(), 1u) << "baseline must be evaluated";
    const double baseProb = baseline.moleculeList[0].probvec.back();

    // Same geometry, history entry pointing at an unrelated partner.
    D1dbrpScenario s = d1dbrp_make_scenario(/*dx*/ 1.5, 0.0, 0.0, /*rate*/ 10.0);
    s.moleculeList[0].prevlist.push_back(7);   // NOT the current partner
    s.moleculeList[0].prevmyface.push_back(0);
    s.moleculeList[0].prevpface.push_back(0);
    s.moleculeList[0].prevsep.push_back(1.6);
    s.moleculeList[0].prevnorm.push_back(0.1); // would strongly change p if used
    s.moleculeList[0].ps_prev.push_back(0.5);
    d1dbrp_invoke(s);

    ASSERT_EQ(s.moleculeList[0].currprevnorm.size(), 1u) << "a norm must still be recorded";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].currprevnorm.back(), 1.0)
        << "history for a different partner must not be applied";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].probvec.back(), baseProb)
        << "probability must match the no-history baseline";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each scenario runs in its own TEST so that a failure in
// one does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(Determine1DBimolecularReactionProbability, WithinRmaxRecordsEncounter)
{
    test_d1dbrp_within_rmax_records_encounter();
}
TEST(Determine1DBimolecularReactionProbability, OutsideRmaxNoEncounter)
{
    test_d1dbrp_outside_rmax_no_encounter();
}
TEST(Determine1DBimolecularReactionProbability, DifferentFiberRejected)
{
    test_d1dbrp_different_fiber_rejected();
}
TEST(Determine1DBimolecularReactionProbability, ZeroRateKeepsProbabilityZero)
{
    test_d1dbrp_zero_rate_keeps_probability_zero();
}
TEST(Determine1DBimolecularReactionProbability, DissociatedSkipsProbability)
{
    test_d1dbrp_dissociated_skips_probability();
}
TEST(Determine1DBimolecularReactionProbability, PromoterUsesSmallerRmax)
{
    test_d1dbrp_promoter_uses_smaller_rmax();
}
TEST(Determine1DBimolecularReactionProbability, ClampsSeparationToBindRadius)
{
    test_d1dbrp_clamps_separation_to_bind_radius();
}
TEST(Determine1DBimolecularReactionProbability, FarPreviousSeparationResetsNorm)
{
    test_d1dbrp_far_previous_separation_resets_norm();
}
TEST(Determine1DBimolecularReactionProbability, ClosePreviousSeparationAppliesNorm)
{
    test_d1dbrp_close_previous_separation_applies_norm();
}
TEST(Determine1DBimolecularReactionProbability, HistoryForOtherPartnerIgnored)
{
    test_d1dbrp_history_for_other_partner_ignored();
}