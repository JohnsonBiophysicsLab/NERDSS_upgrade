/*! \file test_determine_2D_bimolecular_reaction_probability.cpp
 *
 * ### Unit tests for src/reactions/determine_2D_bimolecular_reaction_probability.cpp
 *
 * Function under test:
 *
 * \code
 * void determine_2D_bimolecular_reaction_probability(int simItr, int rxnIndex, int rateIndex,
 *      bool isStateChangeBackRxn, unsigned& DDTableIndex, double* tableIDs, BiMolData& biMolData,
 *      const Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
 *      const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
 *      Membrane& membraneObject, std::vector<gsl_matrix*>& normMatrices,
 *      std::vector<gsl_matrix*>& survMatrices, std::vector<gsl_matrix*>& pirMatrices)
 * \endcode
 *
 * The routine does the following, and each step is checked below:
 *
 *   1. Adds the rotational contribution of both complexes to biMolData.Dtot,
 *      using Dr = 2*|com-iface|*(1 - cos(sqrt(2*Dr.z*dt))) and
 *      Dtot += (Dr1 + Dr2)/(4*dt).
 *   2. Truncates Dtot to a small number of significant figures (so that only a
 *      limited number of unique 2D lookup tables ever get built), and clamps
 *      values below 1e-50 to exactly zero.
 *   3. Computes RMax = 3.5*sqrt(4*Dtot*dt) + bindRadius and asks get_distance()
 *      whether the two interfaces are within RMax.  If they are, a placeholder
 *      probability of 0 is pushed onto both molecules' probvec.
 *   4. If (and only if) the pair is within RMax, the reaction rate is > 0 and
 *      neither molecule just dissociated, it
 *        - looks for an existing 2D lookup table matching (ktemp, Dtot),
 *          where ktemp = rate/length3Dto2D (doubled for symmetric reactions),
 *        - builds a new table (surv/norm/pir matrices) if none matches,
 *          recording (ktemp, Dtot) in tableIDs and incrementing DDTableIndex,
 *        - evaluates the (re-weighted) reaction probability via get_prevSurv()
 *          and DDpirr_pfree_ratio_ps(),
 *        - stores the probability in both molecules' probvec and records the
 *          re-weighting bookkeeping on the *lower indexed* molecule (proA).
 *
 * NOTE: the "full path" test really does build a 2D lookup table by numerical
 * integration, so it is the slowest test in this file.  All the other tests are
 * arranged to avoid table construction so they run instantly.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gsl/gsl_matrix.h>
#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (all prefixed d2dbp_ so they cannot collide with other tests).
// -----------------------------------------------------------------------------
namespace {

//! Everything the function under test needs, bundled for convenience.
struct D2dbpSystem {
    Parameters params;
    Membrane membrane;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
};

/*! \brief Build a minimal single-interface Molecule.
 *
 * The interface coordinate is what get_distance() measures, so it is supplied
 * separately from the center of mass.
 */
Molecule d2dbp_make_molecule(int index, int comIndex, const Coord& com, const Coord& ifaceCrd)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isDissociated = false;
    mol.isEmpty = false;

    // One free (unbound) interface, absolute index == molecule index so the two
    // molecules have distinguishable interface indices.
    Molecule::Iface iface;
    iface.coord = ifaceCrd;
    iface.index = index;
    iface.relIndex = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);
    mol.freelist.push_back(0);

    return mol;
}

/*! \brief Build a Complex holding exactly one member molecule.
 *
 * Only Dr.z is read by the function under test (2D rotational diffusion about
 * the membrane normal), but D is filled in for realism.
 */
Complex d2dbp_make_complex(int index, int memberMol, const Coord& com, double drZ)
{
    Complex targCom;
    targCom.index = index;
    targCom.comCoord = com;
    targCom.mass = 1.0;
    targCom.radius = 1.0;

    // Membrane-bound complex: no motion along z.
    targCom.D.x = 0.25;
    targCom.D.y = 0.25;
    targCom.D.z = 0.0;
    targCom.Dr.x = drZ;
    targCom.Dr.y = drZ;
    targCom.Dr.z = drZ;

    targCom.memberList.push_back(memberMol);
    targCom.OnSurface = true;
    targCom.ncross = 0;
    return targCom;
}

/*! \brief Build a complete two-molecule / two-complex 2D system.
 *
 * \param[in] sepDistance x-separation between the two reacting interfaces (nm)
 * \param[in] rate        macroscopic 3D rate stored in the ForwardRxn
 * \param[in] timeStep    simulation time step (us)
 * \param[in] drZ         rotational diffusion constant of both complexes
 */
D2dbpSystem d2dbp_build_system(double sepDistance, double rate, double timeStep, double drZ)
{
    D2dbpSystem sys;

    // --- parameters -------------------------------------------------------
    sys.params.timeStep = timeStep;
    sys.params.max2DRxns = 4; // keep the tableIDs array tiny for the test
    sys.params.numMolTypes = 1;

    // --- boundary ---------------------------------------------------------
    sys.membrane.isSphere = false;
    sys.membrane.isBox = true;
    sys.membrane.waterBox.x = 500.0;
    sys.membrane.waterBox.y = 500.0;
    sys.membrane.waterBox.z = 500.0;
    sys.membrane.waterBox.volume = 500.0 * 500.0 * 500.0;

    // --- molecules sitting on the lower membrane leaflet ------------------
    const double zSurface = -sys.membrane.waterBox.z / 2.0;
    sys.moleculeList.push_back(
        d2dbp_make_molecule(0, 0, Coord(0.0, 0.0, zSurface), Coord(0.0, 0.0, zSurface)));
    sys.moleculeList.push_back(
        d2dbp_make_molecule(1, 1, Coord(sepDistance, 0.0, zSurface), Coord(sepDistance, 0.0, zSurface)));

    // --- one complex per molecule ----------------------------------------
    sys.complexList.push_back(d2dbp_make_complex(0, 0, Coord(0.0, 0.0, zSurface), drZ));
    sys.complexList.push_back(d2dbp_make_complex(1, 1, Coord(sepDistance, 0.0, zSurface), drZ));

    // --- one bimolecular forward reaction --------------------------------
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.bindRadius = 1.0;
    rxn.length3Dto2D = 10.0; // ktemp = rate / length3Dto2D
    rxn.isOnMem = true;
    rxn.isSymmetric = false;
    rxn.isReversible = false;
    rxn.rateList.emplace_back(rate, std::vector<std::vector<RxnIface>>{});
    rxn.reactantListNew.emplace_back("m0iface", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("m1iface", 0, 1, 0, '\0', false);
    rxn.productListNew.emplace_back("m0iface", 0, 2, 0, '\0', true);
    rxn.productListNew.emplace_back("m1iface", 0, 3, 0, '\0', true);
    rxn.intReactantList = { 0, 1 };
    rxn.intProductList = { 2, 3 };
    sys.forwardRxns.push_back(rxn);

    return sys;
}

/*! \brief Build the BiMolData describing "interface 0 of mol 0 meets interface 0 of mol 1". */
BiMolData d2dbp_make_bimoldata(double dtot, double magMol1, double magMol2, bool swapOrder = false)
{
    if (!swapOrder)
        return BiMolData(0, 1, 0, 1, 0, 0, 0, 1, dtot, magMol1, magMol2);
    // pro1Index > pro2Index: the routine must then treat molecule 0 as "proA"
    return BiMolData(1, 0, 1, 0, 0, 0, 1, 0, dtot, magMol1, magMol2);
}

/*! \brief Clear all per-timestep bookkeeping so the routine can be called again. */
void d2dbp_reset_step_state(std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList)
{
    for (auto& mol : moleculeList) {
        mol.probvec.clear();
        mol.crossbase.clear();
        mol.mycrossint.clear();
        mol.crossrxn.clear();
        mol.currlist.clear();
        mol.currmyface.clear();
        mol.currpface.clear();
        mol.currprevnorm.clear();
        mol.currps_prev.clear();
        mol.currprevsep.clear();
        mol.prevlist.clear();
        mol.prevmyface.clear();
        mol.prevpface.clear();
        mol.prevnorm.clear();
        mol.ps_prev.clear();
        mol.prevsep.clear();
    }
    for (auto& com : complexList)
        com.ncross = 0;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: rotational contribution to Dtot.
//
// Both complexes get Dr.z = 0.01 and |com-iface| = 1, dt = 1.  Then
//   Dr = 2*1*(1 - cos(sqrt(2*0.01*1))) = 0.0199667
//   Dtot = 0.5 + (Dr + Dr)/(4*1) = 0.50998...
// which after the 2-significant-figure truncation becomes exactly 0.51.
//
// The molecules are placed far apart so no lookup table is ever built.
// -----------------------------------------------------------------------------
void test_d2dbp_rotational_contribution_to_dtot()
{
    std::cerr << "\n[TEST] test_d2dbp_rotational_contribution_to_dtot\n"
              << "  Source file: determine_2D_bimolecular_reaction_probability.cpp\n"
              << "  Function:    determine_2D_bimolecular_reaction_probability\n"
              << "  Scenario:    rotational diffusion of both complexes must be folded\n"
              << "               into biMolData.Dtot, then rounded to 2 sig figs.\n"
              << "  Pass:        Dtot(0.5) -> 0.51, no probabilities generated.\n";

    // Far apart (200 nm) so get_distance() returns false -> no 2D table needed.
    D2dbpSystem sys = d2dbp_build_system(/*sep*/ 200.0, /*rate*/ 100.0, /*dt*/ 1.0, /*drZ*/ 0.01);

    BiMolData biMolData = d2dbp_make_bimoldata(/*Dtot*/ 0.5, /*magMol1*/ 1.0, /*magMol2*/ 1.0);

    unsigned DDTableIndex { 0 };
    std::vector<double> tableIDs(2 * sys.params.max2DRxns, 0.0);
    std::vector<gsl_matrix*> normMatrices, survMatrices, pirMatrices;

    std::cerr << "  Dtot before call = " << biMolData.Dtot << '\n';
    determine_2D_bimolecular_reaction_probability(/*simItr*/ 1, /*rxnIndex*/ 0, /*rateIndex*/ 0,
        /*isStateChangeBackRxn*/ false, DDTableIndex, tableIDs.data(), biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.membrane,
        normMatrices, survMatrices, pirMatrices);
    std::cerr << "  Dtot after  call = " << biMolData.Dtot << " (expected 0.51)\n";

    // Rotation must strictly increase Dtot, and the result must be the rounded value.
    EXPECT_GT(biMolData.Dtot, 0.5) << "rotational diffusion should add to Dtot";
    EXPECT_NEAR(biMolData.Dtot, 0.51, 1e-12)
        << "0.5 + rotation (0.00998) should be truncated to 0.51";

    // No pair within RMax -> nothing recorded anywhere.
    EXPECT_EQ(sys.moleculeList[0].probvec.size(), 0u) << "no probability for a far pair";
    EXPECT_EQ(sys.moleculeList[1].probvec.size(), 0u) << "no probability for a far pair";
    EXPECT_EQ(DDTableIndex, 0u) << "no 2D lookup table should be created";
    EXPECT_TRUE(survMatrices.empty()) << "surv matrices must stay empty";
}

// -----------------------------------------------------------------------------
// Test 2: the significant-figure truncation of Dtot over every branch, plus the
//         clamp of tiny values to exactly zero.
//
// Rotational diffusion is switched off (Dr.z = 0) so Dtot only changes through
// the truncation logic, and the molecules stay far apart so nothing else runs.
// -----------------------------------------------------------------------------
void test_d2dbp_dtot_significant_figure_truncation()
{
    std::cerr << "\n[TEST] test_d2dbp_dtot_significant_figure_truncation\n"
              << "  Source file: determine_2D_bimolecular_reaction_probability.cpp\n"
              << "  Function:    determine_2D_bimolecular_reaction_probability\n"
              << "  Scenario:    Dtot is quantised so only a few unique 2D tables exist.\n"
              << "  Pass:        each input maps to the documented quantised value.\n";

    // {input Dtot, expected quantised Dtot}
    const std::vector<std::pair<double, double>> cases {
        { 0.5, 0.5 }, // >0.1  -> 2 decimal places
        { 12.3, 12.3 }, // >10   -> same rule (d*100 rounded)
        { 0.123456, 0.12 }, // >0.1  -> 0.12
        { 0.0567, 0.06 }, // <0.1  -> 0.06
        { 0.00123, 0.001 }, // <0.01 -> 0.001
        { 0.000456, 0.0005 }, // <0.001 -> 0.0005
        { 0.0000234, 0.00002 }, // <0.0001 -> 0.00002
        { 1.0e-60, 0.0 }, // below 1e-50 -> clamped to exactly 0
    };

    for (const auto& oneCase : cases) {
        // 200 nm apart; even for Dtot = 12.3 (RMax ~ 25.6 nm) this is outside RMax.
        D2dbpSystem sys = d2dbp_build_system(/*sep*/ 200.0, /*rate*/ 100.0, /*dt*/ 1.0, /*drZ*/ 0.0);
        BiMolData biMolData = d2dbp_make_bimoldata(oneCase.first, 1.0, 1.0);

        unsigned DDTableIndex { 0 };
        std::vector<double> tableIDs(2 * sys.params.max2DRxns, 0.0);
        std::vector<gsl_matrix*> normMatrices, survMatrices, pirMatrices;

        determine_2D_bimolecular_reaction_probability(1, 0, 0, false, DDTableIndex, tableIDs.data(),
            biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
            sys.membrane, normMatrices, survMatrices, pirMatrices);

        std::cerr << "  -> Dtot " << oneCase.first << " quantised to " << biMolData.Dtot
                  << " (expected " << oneCase.second << ")\n";
        EXPECT_NEAR(biMolData.Dtot, oneCase.second, 1e-12)
            << "Dtot " << oneCase.first << " should quantise to " << oneCase.second;

        // Nothing should have been created for a far pair.
        EXPECT_EQ(DDTableIndex, 0u) << "no table for an out-of-range pair";
    }
}

// -----------------------------------------------------------------------------
// Test 3: pair outside RMax -> the routine is a complete no-op with respect to
//         probabilities, tables and re-weighting bookkeeping.
// -----------------------------------------------------------------------------
void test_d2dbp_outside_rmax_is_noop()
{
    std::cerr << "\n[TEST] test_d2dbp_outside_rmax_is_noop\n"
              << "  Source file: determine_2D_bimolecular_reaction_probability.cpp\n"
              << "  Function:    determine_2D_bimolecular_reaction_probability\n"
              << "  Scenario:    interfaces separated by 200 nm, RMax ~ 5.9 nm.\n"
              << "  Pass:        probvec/currlist untouched, DDTableIndex stays 0.\n";

    D2dbpSystem sys = d2dbp_build_system(/*sep*/ 200.0, /*rate*/ 100.0, /*dt*/ 1.0, /*drZ*/ 0.0);
    BiMolData biMolData = d2dbp_make_bimoldata(0.5, 1.0, 1.0);

    unsigned DDTableIndex { 0 };
    std::vector<double> tableIDs(2 * sys.params.max2DRxns, 0.0);
    std::vector<gsl_matrix*> normMatrices, survMatrices, pirMatrices;

    const double rMax = 3.5 * std::sqrt(4.0 * 0.5 * sys.params.timeStep) + sys.forwardRxns[0].bindRadius;
    std::cerr << "  RMax = " << rMax << " nm, actual separation = 200 nm\n";

    determine_2D_bimolecular_reaction_probability(5, 0, 0, false, DDTableIndex, tableIDs.data(),
        biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
        sys.membrane, normMatrices, survMatrices, pirMatrices);

    EXPECT_EQ(sys.moleculeList[0].probvec.size(), 0u) << "mol 0 probvec must stay empty";
    EXPECT_EQ(sys.moleculeList[1].probvec.size(), 0u) << "mol 1 probvec must stay empty";
    EXPECT_EQ(sys.moleculeList[0].currlist.size(), 0u) << "no re-weighting record expected";
    EXPECT_EQ(DDTableIndex, 0u) << "no lookup table should be built";
    EXPECT_EQ(tableIDs[0], 0.0) << "tableIDs must not be written";
    EXPECT_TRUE(normMatrices.empty()) << "norm matrices must stay empty";
    EXPECT_TRUE(pirMatrices.empty()) << "pir matrices must stay empty";
}

// -----------------------------------------------------------------------------
// Test 4: pair inside RMax but the reaction rate is zero.
//
// get_distance() succeeds so a placeholder 0 probability is pushed for each
// molecule, but the (rate > 0) guard prevents any table construction and any
// re-weighting bookkeeping.
// -----------------------------------------------------------------------------
void test_d2dbp_zero_rate_pushes_zero_probability()
{
    std::cerr << "\n[TEST] test_d2dbp_zero_rate_pushes_zero_probability\n"
              << "  Source file: determine_2D_bimolecular_reaction_probability.cpp\n"
              << "  Function:    determine_2D_bimolecular_reaction_probability\n"
              << "  Scenario:    interfaces 1.2 nm apart (inside RMax) but rate == 0.\n"
              << "  Pass:        one placeholder prob of 0 per molecule, no table.\n";

    D2dbpSystem sys = d2dbp_build_system(/*sep*/ 1.2, /*rate*/ 0.0, /*dt*/ 1.0, /*drZ*/ 0.0);
    BiMolData biMolData = d2dbp_make_bimoldata(0.5, 1.0, 1.0);

    unsigned DDTableIndex { 0 };
    std::vector<double> tableIDs(2 * sys.params.max2DRxns, 0.0);
    std::vector<gsl_matrix*> normMatrices, survMatrices, pirMatrices;

    determine_2D_bimolecular_reaction_probability(7, 0, 0, false, DDTableIndex, tableIDs.data(),
        biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
        sys.membrane, normMatrices, survMatrices, pirMatrices);

    // The placeholder pushes prove that get_distance() reported "within RMax".
    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u) << "mol 0 should get one placeholder entry";
    ASSERT_EQ(sys.moleculeList[1].probvec.size(), 1u) << "mol 1 should get one placeholder entry";
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].probvec.back(), 0.0) << "zero rate -> zero probability";
    EXPECT_DOUBLE_EQ(sys.moleculeList[1].probvec.back(), 0.0) << "zero rate -> zero probability";

    EXPECT_EQ(DDTableIndex, 0u) << "zero rate must not trigger 2D table construction";
    EXPECT_EQ(sys.moleculeList[0].currlist.size(), 0u) << "no re-weighting record for zero rate";
    std::cerr << "  probvec = " << sys.moleculeList[0].probvec.back()
              << ", DDTableIndex = " << DDTableIndex << '\n';
}

// -----------------------------------------------------------------------------
// Test 5: a just-dissociated molecule must not be allowed to re-bind this step.
//
// The pair is inside RMax so the placeholder zeros are pushed, but the
// isDissociated guard skips the whole probability evaluation.
// -----------------------------------------------------------------------------
void test_d2dbp_dissociated_molecule_skips_probability()
{
    std::cerr << "\n[TEST] test_d2dbp_dissociated_molecule_skips_probability\n"
              << "  Source file: determine_2D_bimolecular_reaction_probability.cpp\n"
              << "  Function:    determine_2D_bimolecular_reaction_probability\n"
              << "  Scenario:    inside RMax, rate > 0, but mol 0 just dissociated.\n"
              << "  Pass:        probability stays 0 and no table is created.\n";

    D2dbpSystem sys = d2dbp_build_system(/*sep*/ 1.2, /*rate*/ 100.0, /*dt*/ 1.0, /*drZ*/ 0.0);
    sys.moleculeList[0].isDissociated = true; // just came apart this timestep

    BiMolData biMolData = d2dbp_make_bimoldata(0.5, 1.0, 1.0);

    unsigned DDTableIndex { 0 };
    std::vector<double> tableIDs(2 * sys.params.max2DRxns, 0.0);
    std::vector<gsl_matrix*> normMatrices, survMatrices, pirMatrices;

    determine_2D_bimolecular_reaction_probability(11, 0, 0, false, DDTableIndex, tableIDs.data(),
        biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
        sys.membrane, normMatrices, survMatrices, pirMatrices);

    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u) << "placeholder still pushed (within RMax)";
    ASSERT_EQ(sys.moleculeList[1].probvec.size(), 1u) << "placeholder still pushed (within RMax)";
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].probvec.back(), 0.0)
        << "dissociated molecule must have zero binding probability";
    EXPECT_DOUBLE_EQ(sys.moleculeList[1].probvec.back(), 0.0)
        << "partner of a dissociated molecule must have zero probability";
    EXPECT_EQ(DDTableIndex, 0u) << "no table for a dissociated pair";
    EXPECT_EQ(sys.moleculeList[0].currlist.size(), 0u) << "no re-weighting record";
}

// -----------------------------------------------------------------------------
// Test 6 (the heavy one): the complete reactive path.
//
// A single 2D lookup table is built by numerical integration and then reused by
// several follow-up calls so the cost is paid only once.  Sub-cases:
//
//   A: fresh call            -> table built, DDTableIndex 1, prob == get_prevSurv()
//   B: identical repeat call -> table reused (DDTableIndex still 1), same prob
//   C: prevsep >= RMax       -> re-weighting restarted (currnorm = 1), same prob
//   D: prevsep <  RMax       -> DDpirr_pfree_ratio_ps() re-weighting path
//   E: isSymmetric == true   -> ktemp doubled; rate halved so it hits table A
//   F: swapped pro indices   -> bookkeeping still stored on the lower index
// -----------------------------------------------------------------------------
void test_d2dbp_full_reactive_path_and_table_reuse()
{
    std::cerr << "\n[TEST] test_d2dbp_full_reactive_path_and_table_reuse\n"
              << "  Source file: determine_2D_bimolecular_reaction_probability.cpp\n"
              << "  Function:    determine_2D_bimolecular_reaction_probability\n"
              << "  Scenario:    interfaces 1.2 nm apart, bindRadius 1 nm, rate > 0.\n"
              << "  Pass:        a 2D table is built once, reused afterwards, and the\n"
              << "               reaction probability matches get_prevSurv().\n"
              << "  NOTE:        this builds real 2D matrices -> slowest test here.\n";

    const double timeStep = 0.1;
    const double dtotInput = 0.5;
    const double rate = 100.0; // ktemp = rate / length3Dto2D = 10 nm^2/us
    const double sep = 1.2; // interface separation (bindRadius is 1.0)

    D2dbpSystem sys = d2dbp_build_system(sep, rate, timeStep, /*drZ*/ 0.0);
    const double bindRadius = sys.forwardRxns[0].bindRadius;
    const double ktemp = rate / sys.forwardRxns[0].length3Dto2D;
    const double rMax = 3.5 * std::sqrt(4.0 * dtotInput * timeStep) + bindRadius;

    std::cerr << "  bindRadius = " << bindRadius << ", Dtot = " << dtotInput
              << ", dt = " << timeStep << ", RMax = " << rMax << ", ktemp = " << ktemp << '\n';

    unsigned DDTableIndex { 0 };
    std::vector<double> tableIDs(2 * sys.params.max2DRxns, 0.0);
    std::vector<gsl_matrix*> normMatrices, survMatrices, pirMatrices;

    // ---------------- sub-case A: first (table building) call ----------------
    std::cerr << "  [A] first call: expect a brand new 2D lookup table...\n";
    {
        BiMolData biMolData = d2dbp_make_bimoldata(dtotInput, 1.0, 1.0);
        d2dbp_reset_step_state(sys.moleculeList, sys.complexList);

        determine_2D_bimolecular_reaction_probability(100, 0, 0, false, DDTableIndex, tableIDs.data(),
            biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
            sys.membrane, normMatrices, survMatrices, pirMatrices);

        // Dtot should be unchanged (no rotation, and 0.5 is already quantised).
        EXPECT_NEAR(biMolData.Dtot, dtotInput, 1e-12) << "Dtot should stay 0.5";

        // Exactly one table registered, with the right identifiers.
        EXPECT_EQ(DDTableIndex, 1u) << "one unique 2D reaction table expected";
        EXPECT_NEAR(tableIDs[0], ktemp, 1e-10) << "first tableID slot holds ktemp";
        EXPECT_NEAR(tableIDs[sys.params.max2DRxns], dtotInput, 1e-10)
            << "second tableID dimension holds Dtot";

        // Matrices allocated with the sizes size_lookup() dictates.
        ASSERT_EQ(survMatrices.size(), 1u) << "one surv matrix";
        ASSERT_EQ(normMatrices.size(), 1u) << "one norm matrix";
        ASSERT_EQ(pirMatrices.size(), 1u) << "one pir matrix";
        const size_t veclen = size_lookup(bindRadius, dtotInput, sys.params, rMax);
        std::cerr << "      size_lookup() veclen = " << veclen
                  << ", survMatrix is " << survMatrices[0]->size1 << "x" << survMatrices[0]->size2
                  << ", pirMatrix is " << pirMatrices[0]->size1 << "x" << pirMatrices[0]->size2 << '\n';
        EXPECT_EQ(survMatrices[0]->size1, 2u) << "surv matrix has 2 rows";
        EXPECT_EQ(survMatrices[0]->size2, veclen) << "surv matrix width == size_lookup()";
        EXPECT_EQ(normMatrices[0]->size1, 2u) << "norm matrix has 2 rows";
        EXPECT_EQ(normMatrices[0]->size2, veclen) << "norm matrix width == size_lookup()";
        EXPECT_EQ(pirMatrices[0]->size1, veclen) << "pir matrix is square, veclen rows";
        EXPECT_EQ(pirMatrices[0]->size2, veclen) << "pir matrix is square, veclen cols";

        // A real probability was produced for both partners.
        ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u) << "mol 0 gets one probability";
        ASSERT_EQ(sys.moleculeList[1].probvec.size(), 1u) << "mol 1 gets one probability";
        const double prob = sys.moleculeList[0].probvec.back();
        std::cerr << "      reaction probability = " << prob << '\n';
        EXPECT_GT(prob, 0.0) << "a reactive pair should have a positive probability";
        EXPECT_LE(prob, 1.000001) << "probability must not exceed unity";
        EXPECT_DOUBLE_EQ(sys.moleculeList[1].probvec.back(), prob)
            << "both partners must record the same probability";

        // With no re-weighting history, currnorm == 1, so the probability must
        // be exactly the survival-table lookup.
        const double expected = get_prevSurv(survMatrices[0], dtotInput, timeStep,
            sys.moleculeList[0].currprevsep.back(), bindRadius);
        EXPECT_NEAR(prob, expected, 1e-12)
            << "unweighted probability must equal get_prevSurv() of the stored separation";

        // Re-weighting bookkeeping is stored on the lower-indexed molecule (proA = 0).
        ASSERT_EQ(sys.moleculeList[0].currlist.size(), 1u) << "proA records the encounter";
        EXPECT_EQ(sys.moleculeList[0].currlist.back(), 1) << "partner index is 1";
        EXPECT_EQ(sys.moleculeList[0].currmyface.back(), 0) << "own interface index is 0";
        EXPECT_EQ(sys.moleculeList[0].currpface.back(), 0) << "partner interface index is 0";
        EXPECT_NEAR(sys.moleculeList[0].currprevsep.back(), sep, 1e-9)
            << "stored separation should be the actual interface distance";
        EXPECT_NEAR(sys.moleculeList[0].currprevnorm.back(), 1.0, 1e-12)
            << "no history -> normalisation of 1";
        EXPECT_NEAR(sys.moleculeList[0].currps_prev.back(), 1.0 - prob, 1e-12)
            << "stored survival must be 1 - probability";
        EXPECT_EQ(sys.moleculeList[1].currlist.size(), 0u)
            << "the higher-indexed molecule stores nothing";
    }

    // Remember the baseline probability for comparison in later sub-cases.
    double baselineProb = 0.0;
    {
        BiMolData biMolData = d2dbp_make_bimoldata(dtotInput, 1.0, 1.0);
        d2dbp_reset_step_state(sys.moleculeList, sys.complexList);

        // ---------------- sub-case B: identical call reuses the table --------
        std::cerr << "  [B] repeat call: expect the existing table to be reused...\n";
        determine_2D_bimolecular_reaction_probability(101, 0, 0, false, DDTableIndex, tableIDs.data(),
            biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
            sys.membrane, normMatrices, survMatrices, pirMatrices);

        EXPECT_EQ(DDTableIndex, 1u) << "identical (ktemp, Dtot) must not allocate a new table";
        EXPECT_EQ(survMatrices.size(), 1u) << "matrix vectors must not grow";
        ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u);
        baselineProb = sys.moleculeList[0].probvec.back();
        std::cerr << "      probability = " << baselineProb << '\n';
        EXPECT_GT(baselineProb, 0.0) << "reused table must still produce a probability";
    }

    // ---------------- sub-case C: previous separation outside RMax ------------
    std::cerr << "  [C] history with prevsep >= RMax: re-weighting must restart (currnorm = 1)...\n";
    {
        BiMolData biMolData = d2dbp_make_bimoldata(dtotInput, 1.0, 1.0);
        d2dbp_reset_step_state(sys.moleculeList, sys.complexList);

        // Encounter history: same partner/faces, but the pair was previously in 3D.
        sys.moleculeList[0].prevlist.push_back(1);
        sys.moleculeList[0].prevmyface.push_back(0);
        sys.moleculeList[0].prevpface.push_back(0);
        sys.moleculeList[0].prevnorm.push_back(0.25); // must be discarded
        sys.moleculeList[0].ps_prev.push_back(0.5);
        sys.moleculeList[0].prevsep.push_back(100.0); // >= RMax

        determine_2D_bimolecular_reaction_probability(102, 0, 0, false, DDTableIndex, tableIDs.data(),
            biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
            sys.membrane, normMatrices, survMatrices, pirMatrices);

        ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u);
        std::cerr << "      probability = " << sys.moleculeList[0].probvec.back()
                  << " (baseline " << baselineProb << ")\n";
        EXPECT_NEAR(sys.moleculeList[0].currprevnorm.back(), 1.0, 1e-12)
            << "prevsep >= RMax must reset the normalisation to 1";
        EXPECT_NEAR(sys.moleculeList[0].probvec.back(), baselineProb, 1e-12)
            << "with currnorm reset the probability equals the unweighted baseline";
        EXPECT_EQ(DDTableIndex, 1u) << "still only one table";
    }

    // ---------------- sub-case D: previous separation inside RMax -------------
    std::cerr << "  [D] history with prevsep < RMax: DDpirr_pfree_ratio_ps() re-weighting...\n";
    {
        BiMolData biMolData = d2dbp_make_bimoldata(dtotInput, 1.0, 1.0);
        d2dbp_reset_step_state(sys.moleculeList, sys.complexList);

        sys.moleculeList[0].prevlist.push_back(1);
        sys.moleculeList[0].prevmyface.push_back(0);
        sys.moleculeList[0].prevpface.push_back(0);
        sys.moleculeList[0].prevnorm.push_back(1.0);
        sys.moleculeList[0].ps_prev.push_back(0.9);
        sys.moleculeList[0].prevsep.push_back(1.4); // inside RMax, outside sigma

        determine_2D_bimolecular_reaction_probability(103, 0, 0, false, DDTableIndex, tableIDs.data(),
            biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
            sys.membrane, normMatrices, survMatrices, pirMatrices);

        ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u);
        ASSERT_EQ(sys.moleculeList[0].currprevnorm.size(), 1u);
        const double prob = sys.moleculeList[0].probvec.back();
        const double norm = sys.moleculeList[0].currprevnorm.back();
        std::cerr << "      re-weighted probability = " << prob << ", currnorm = " << norm << '\n';

        // The re-weighted values must remain physical and self-consistent.
        EXPECT_TRUE(std::isfinite(prob)) << "re-weighted probability must be finite";
        EXPECT_GE(prob, 0.0) << "probability cannot be negative";
        EXPECT_TRUE(std::isfinite(norm)) << "normalisation must be finite";
        EXPECT_GT(norm, 0.0) << "normalisation must be positive";
        EXPECT_NEAR(sys.moleculeList[0].currps_prev.back(), 1.0 - prob, 1e-12)
            << "stored survival must always be 1 - probability";
        EXPECT_DOUBLE_EQ(sys.moleculeList[1].probvec.back(), prob)
            << "both partners share the re-weighted probability";
        EXPECT_EQ(DDTableIndex, 1u) << "no new table for a re-weighted call";
    }

    // ---------------- sub-case E: symmetric reaction doubles ktemp -----------
    std::cerr << "  [E] symmetric reaction with half the rate must map onto the SAME table...\n";
    {
        // rate/2 with isSymmetric == true gives ktemp*2 == the tabulated ktemp,
        // so no new table may be allocated.
        sys.forwardRxns[0].isSymmetric = true;
        sys.forwardRxns[0].rateList[0].rate = rate / 2.0;

        BiMolData biMolData = d2dbp_make_bimoldata(dtotInput, 1.0, 1.0);
        d2dbp_reset_step_state(sys.moleculeList, sys.complexList);

        determine_2D_bimolecular_reaction_probability(104, 0, 0, false, DDTableIndex, tableIDs.data(),
            biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
            sys.membrane, normMatrices, survMatrices, pirMatrices);

        EXPECT_EQ(DDTableIndex, 1u)
            << "symmetric doubling of ktemp should match the existing table (no new table)";
        ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u);
        std::cerr << "      probability = " << sys.moleculeList[0].probvec.back() << '\n';
        EXPECT_NEAR(sys.moleculeList[0].probvec.back(), baselineProb, 1e-12)
            << "same effective ktemp/Dtot -> same probability as the baseline";

        // restore the reaction for the final sub-case
        sys.forwardRxns[0].isSymmetric = false;
        sys.forwardRxns[0].rateList[0].rate = rate;
    }

    // ---------------- sub-case F: swapped molecule order ---------------------
    std::cerr << "  [F] swapped pro1/pro2: bookkeeping must land on the lower molecule index...\n";
    {
        BiMolData biMolData = d2dbp_make_bimoldata(dtotInput, 1.0, 1.0, /*swapOrder*/ true);
        d2dbp_reset_step_state(sys.moleculeList, sys.complexList);

        determine_2D_bimolecular_reaction_probability(105, 0, 0, false, DDTableIndex, tableIDs.data(),
            biMolData, sys.params, sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns,
            sys.membrane, normMatrices, survMatrices, pirMatrices);

        EXPECT_EQ(DDTableIndex, 1u) << "still one table";
        ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u) << "mol 0 still gets a probability";
        ASSERT_EQ(sys.moleculeList[1].probvec.size(), 1u) << "mol 1 still gets a probability";
        EXPECT_NEAR(sys.moleculeList[0].probvec.back(), baselineProb, 1e-12)
            << "probability is independent of the argument order";

        // proA is always min(pro1Index, pro2Index) == 0.
        ASSERT_EQ(sys.moleculeList[0].currlist.size(), 1u)
            << "the lower-indexed molecule must own the re-weighting record";
        EXPECT_EQ(sys.moleculeList[0].currlist.back(), 1) << "record points at molecule 1";
        EXPECT_EQ(sys.moleculeList[1].currlist.size(), 0u)
            << "the higher-indexed molecule must not own a record";
        std::cerr << "      currlist owner is molecule 0 as required\n";
    }

    // ---------------- clean up the GSL allocations ---------------------------
    std::cerr << "  Releasing GSL matrices...\n";
    for (auto* m : survMatrices)
        if (m != nullptr)
            gsl_matrix_free(m);
    for (auto* m : normMatrices)
        if (m != nullptr)
            gsl_matrix_free(m);
    for (auto* m : pirMatrices)
        if (m != nullptr)
            gsl_matrix_free(m);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so that
// a failure in one does not stop the others from running.
// -----------------------------------------------------------------------------
TEST(Determine2DBimolecularReactionProbability, RotationalContributionToDtot)
{
    test_d2dbp_rotational_contribution_to_dtot();
}

TEST(Determine2DBimolecularReactionProbability, DtotSignificantFigureTruncation)
{
    test_d2dbp_dtot_significant_figure_truncation();
}

TEST(Determine2DBimolecularReactionProbability, OutsideRMaxIsNoop)
{
    test_d2dbp_outside_rmax_is_noop();
}

TEST(Determine2DBimolecularReactionProbability, ZeroRatePushesZeroProbability)
{
    test_d2dbp_zero_rate_pushes_zero_probability();
}

TEST(Determine2DBimolecularReactionProbability, DissociatedMoleculeSkipsProbability)
{
    test_d2dbp_dissociated_molecule_skips_probability();
}

TEST(Determine2DBimolecularReactionProbability, FullReactivePathAndTableReuse)
{
    test_d2dbp_full_reactive_path_and_table_reuse();
}