/*! \file test_determine_3D_bimolecular_reaction_probability.cpp
 *
 * ### Unit test for src/reactions/determine_3D_bimolecular_reaction_probability.cpp
 *
 * Function under test:
 *
 *     void determine_3D_bimolecular_reaction_probability(
 *         int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
 *         BiMolData& biMolData, const Parameters& params,
 *         std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
 *         const std::vector<ForwardRxn>& forwardRxns,
 *         const std::vector<BackRxn>& backRxns);
 *
 * What the routine does (and therefore what we verify here):
 *   1. Adds a rotational contribution to `biMolData.Dtot` for each of the two
 *      complexes.  Which formula is used depends on whether the complex has
 *      `D.z == 0` (membrane bound, 2D rotation, /(4*dt)) or `D.z != 0`
 *      (solution, 3D rotation, /(6*dt)).
 *   2. Computes Rmax = 3*sqrt(6*Dtot*dt) + bindRadius and asks get_distance()
 *      whether the two reacting interfaces are inside Rmax.
 *   3. If inside Rmax, a placeholder 0 probability is appended to the probvec
 *      of BOTH molecules.
 *   4. If neither molecule just dissociated and the reaction rate is > 0, the
 *      association probability is evaluated (passocF), optionally reweighted by
 *      the pair history stored in prev* vectors, written into both probvecs and
 *      recorded in the curr* bookkeeping vectors of the lower-indexed molecule.
 *   5. kact is doubled when either complex is OnSurface (and not on a fiber),
 *      and doubled again for a symmetric reaction.
 *
 * All console output is verbose so the reader can follow exactly which scenario
 * is being exercised and what the pass criterion is.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (file-internal, uniquely prefixed d3bp_ to avoid collisions).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Bundle of everything the function under test needs. */
struct D3bpFixture {
    Parameters params {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {}; // unused by the routine but required by the signature
};

/*! \brief Create a single-interface molecule.
 *
 * \param[in] index        index of this molecule inside moleculeList
 * \param[in] comIndex     index of the parent complex inside complexList
 * \param[in] comCoord     center of mass coordinate
 * \param[in] ifaceCoord   absolute coordinate of the (single) interface
 * \param[in] absIfaceIndex absolute state index given to the interface
 */
Molecule d3bp_make_molecule(int index, int comIndex, const Coord& comCoord,
    const Coord& ifaceCoord, int absIfaceIndex)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isEmpty = false;
    mol.isDissociated = false;
    mol.comCoord = comCoord;

    // One interface, coincident with whatever coordinate the caller asked for.
    Molecule::Iface iface;
    iface.coord = ifaceCoord;
    iface.index = absIfaceIndex;
    iface.relIndex = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // The single interface is free, so it can react.
    mol.freelist.clear();
    mol.freelist.push_back(0);

    return mol;
}

/*! \brief Create a complex owning exactly one member molecule.
 *
 * \param[in] index    index inside complexList
 * \param[in] memberId index of the member molecule inside moleculeList
 * \param[in] comCoord center of mass coordinate
 * \param[in] transD   translational diffusion constant used for x, y and z
 * \param[in] rotD     rotational diffusion constant used for x, y and z
 */
Complex d3bp_make_complex(int index, int memberId, const Coord& comCoord,
    double transD, double rotD)
{
    Complex com;
    com.index = index;
    com.comCoord = comCoord;
    com.radius = 1.0;
    com.mass = 1.0;
    com.D = Coord(transD, transD, transD);
    com.Dr = Coord(rotD, rotD, rotD);
    com.isEmpty = false;
    com.OnSurface = false;
    com.onFiber = false;
    com.memberList.clear();
    com.memberList.push_back(memberId);
    return com;
}

/*! \brief Build a two-molecule / two-complex system with one bimolecular reaction.
 *
 * Molecule 0 (complex 0) sits at the origin; molecule 1 (complex 1) is displaced
 * along +x by `ifaceSeparation`, so the interface-interface distance is exactly
 * `ifaceSeparation`.
 *
 * \param[in] ifaceSeparation distance (nm) between the two reacting interfaces
 * \param[in] rate            micro association rate for rateList[0]
 */
D3bpFixture d3bp_make_fixture(double ifaceSeparation, double rate)
{
    D3bpFixture sys;

    // Small time step (microseconds), as used throughout NERDSS.
    sys.params.timeStep = 0.1;
    sys.params.numMolTypes = 1;

    // Two molecules, each with one interface, separated along the x axis.
    sys.moleculeList.push_back(
        d3bp_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0), 0));
    sys.moleculeList.push_back(d3bp_make_molecule(
        1, 1, Coord(ifaceSeparation, 0.0, 0.0), Coord(ifaceSeparation, 0.0, 0.0), 1));

    // Two independent complexes, both diffusing in solution (D.z != 0).
    sys.complexList.push_back(d3bp_make_complex(0, 0, Coord(0.0, 0.0, 0.0), 10.0, 0.01));
    sys.complexList.push_back(
        d3bp_make_complex(1, 1, Coord(ifaceSeparation, 0.0, 0.0), 10.0, 0.01));

    // One forward (bimolecular) reaction between the two interfaces.
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = 1.0;
    rxn.isSymmetric = false;
    rxn.isOnMem = false;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.bindRadSameCom = 1.1;
    rxn.rateList.emplace_back(rate, std::vector<std::vector<RxnIface>> {});
    sys.forwardRxns.push_back(rxn);

    return sys;
}

/*! \brief Convenience constructor for the BiMolData describing our pair.
 *
 * pro1 = 0 (complex 0, rel iface 0, abs iface 0),
 * pro2 = 1 (complex 1, rel iface 0, abs iface 1).
 */
BiMolData d3bp_make_bimoldata(double dtot, double magMol1, double magMol2)
{
    return BiMolData(/*pro1*/ 0, /*pro2*/ 1, /*com1*/ 0, /*com2*/ 1,
        /*relIface1*/ 0, /*relIface2*/ 0, /*absIface1*/ 0, /*absIface2*/ 1,
        dtot, magMol1, magMol2);
}

/*! \brief Reproduce the rotational Dtot contribution of one complex.
 *
 * Mirrors the arithmetic in the routine so we can validate it independently.
 */
double d3bp_expected_rot_contribution(double transDz, double rotDz, double mag, double dt)
{
    if (std::abs(transDz - 0.0) < 1E-10) {
        // Membrane-bound complex: 2D rotation, normalised by 4*dt.
        const double cf = std::cos(std::sqrt(2.0 * rotDz * dt));
        return (2.0 * mag * (1.0 - cf)) / (4.0 * dt);
    }
    // Solution complex: 3D rotation, normalised by 6*dt.
    const double cf = std::cos(std::sqrt(4.0 * rotDz * dt));
    return (2.0 * mag * (1.0 - cf)) / (6.0 * dt);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: rotational contributions are added to Dtot with the correct branch.
// -----------------------------------------------------------------------------
void test_d3bp_rotational_dtot_contribution()
{
    std::cerr << "\n[TEST] test_d3bp_rotational_dtot_contribution\n"
              << "  Source file:   determine_3D_bimolecular_reaction_probability.cpp\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      complex 0 is membrane bound (D.z == 0 -> /4dt branch),\n"
              << "                 complex 1 diffuses in solution (D.z != 0 -> /6dt branch).\n"
              << "                 Molecules are far apart so only Dtot bookkeeping runs.\n"
              << "  Pass criteria: Dtot equals the hand-computed rotational sum.\n";

    // Interfaces 100 nm apart: guaranteed outside Rmax (~11 nm here).
    D3bpFixture sys = d3bp_make_fixture(100.0, 100.0);

    // Complex 0 pinned to the membrane: no z translation.
    sys.complexList[0].D = Coord(1.0, 1.0, 0.0);
    sys.complexList[0].Dr = Coord(0.01, 0.01, 0.01);
    // Complex 1 free in solution.
    sys.complexList[1].D = Coord(10.0, 10.0, 10.0);
    sys.complexList[1].Dr = Coord(0.02, 0.02, 0.02);

    const double startDtot = 20.0;
    const double magMol1 = 4.0;
    const double magMol2 = 9.0;
    BiMolData biMolData = d3bp_make_bimoldata(startDtot, magMol1, magMol2);

    const double expected = startDtot
        + d3bp_expected_rot_contribution(0.0, 0.01, magMol1, sys.params.timeStep)
        + d3bp_expected_rot_contribution(10.0, 0.02, magMol2, sys.params.timeStep);

    std::cerr << "  Calling determine_3D_bimolecular_reaction_probability...\n";
    determine_3D_bimolecular_reaction_probability(/*simItr*/ 1, /*rxnIndex*/ 0,
        /*rateIndex*/ 0, /*isStateChangeBackRxn*/ false, biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns);

    std::cerr << "  Dtot after call = " << biMolData.Dtot
              << ", expected = " << expected << '\n';

    EXPECT_NEAR(biMolData.Dtot, expected, 1e-10)
        << "Dtot must include both rotational contributions using the correct branches";
    EXPECT_GT(biMolData.Dtot, startDtot)
        << "Rotational diffusion should strictly increase Dtot for nonzero Dr";
}

// -----------------------------------------------------------------------------
// Test 2: zero rotational diffusion leaves Dtot untouched.
// -----------------------------------------------------------------------------
void test_d3bp_zero_rotation_leaves_dtot()
{
    std::cerr << "\n[TEST] test_d3bp_zero_rotation_leaves_dtot\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      both complexes have Dr.z == 0, molecules far apart.\n"
              << "  Pass criteria: cos(0) == 1 so no rotational term is added; Dtot is\n"
              << "                 exactly the input value.\n";

    D3bpFixture sys = d3bp_make_fixture(100.0, 100.0);
    sys.complexList[0].Dr = Coord(0.0, 0.0, 0.0);
    sys.complexList[1].Dr = Coord(0.0, 0.0, 0.0);

    BiMolData biMolData = d3bp_make_bimoldata(15.0, 4.0, 9.0);

    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns);

    std::cerr << "  Dtot after call = " << biMolData.Dtot << " (expected 15)\n";
    EXPECT_DOUBLE_EQ(biMolData.Dtot, 15.0)
        << "Dtot must be unchanged when both rotational diffusion constants are zero";
}

// -----------------------------------------------------------------------------
// Test 3: pair outside Rmax -> nothing is recorded at all.
// -----------------------------------------------------------------------------
void test_d3bp_outside_rmax_records_nothing()
{
    std::cerr << "\n[TEST] test_d3bp_outside_rmax_records_nothing\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      interfaces separated by 100 nm while Rmax is ~11 nm.\n"
              << "  Pass criteria: no probability appended to either probvec and no\n"
              << "                 reweighting bookkeeping written.\n";

    D3bpFixture sys = d3bp_make_fixture(100.0, 100.0);
    BiMolData biMolData = d3bp_make_bimoldata(20.0, 4.0, 4.0);

    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns);

    std::cerr << "  probvec sizes: mol0 = " << sys.moleculeList[0].probvec.size()
              << ", mol1 = " << sys.moleculeList[1].probvec.size() << '\n';

    EXPECT_EQ(sys.moleculeList[0].probvec.size(), 0u)
        << "molecule 0 probvec should stay empty when outside Rmax";
    EXPECT_EQ(sys.moleculeList[1].probvec.size(), 0u)
        << "molecule 1 probvec should stay empty when outside Rmax";
    EXPECT_EQ(sys.moleculeList[0].currlist.size(), 0u)
        << "no reweighting partner should be recorded when outside Rmax";
    EXPECT_EQ(sys.moleculeList[0].currprevsep.size(), 0u)
        << "no separation should be recorded when outside Rmax";
}

// -----------------------------------------------------------------------------
// Test 4: pair inside Rmax with a positive rate -> probability + bookkeeping.
// -----------------------------------------------------------------------------
void test_d3bp_inside_rmax_sets_probability()
{
    std::cerr << "\n[TEST] test_d3bp_inside_rmax_sets_probability\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      interfaces 2 nm apart (bindRadius 1 nm), rate = 100.\n"
              << "  Pass criteria: one probability appended to each probvec, both equal,\n"
              << "                 0 < p <= 1, and the curr* bookkeeping of the lower\n"
              << "                 indexed molecule holds the separation/partner/faces\n"
              << "                 with currps_prev == 1 - p.\n";

    D3bpFixture sys = d3bp_make_fixture(2.0, 100.0);
    BiMolData biMolData = d3bp_make_bimoldata(20.0, 4.0, 4.0);

    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns);

    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u)
        << "one probability must be appended for molecule 0";
    ASSERT_EQ(sys.moleculeList[1].probvec.size(), 1u)
        << "one probability must be appended for molecule 1";

    const double prob = sys.moleculeList[0].probvec.back();
    std::cerr << "  association probability = " << prob << '\n';

    EXPECT_GT(prob, 0.0) << "probability should be strictly positive inside Rmax";
    EXPECT_LE(prob, 1.0) << "probability must not exceed unity for this time step";
    EXPECT_DOUBLE_EQ(sys.moleculeList[1].probvec.back(), prob)
        << "both partners must be assigned the identical pair probability";

    // Reweighting bookkeeping is stored on the lower-indexed molecule (proA == 0).
    ASSERT_EQ(sys.moleculeList[0].currlist.size(), 1u)
        << "partner index must be recorded on the lower-indexed molecule";
    EXPECT_EQ(sys.moleculeList[0].currlist.back(), 1)
        << "recorded partner should be molecule 1";
    ASSERT_EQ(sys.moleculeList[0].currmyface.size(), 1u);
    EXPECT_EQ(sys.moleculeList[0].currmyface.back(), 0) << "own interface index is 0";
    ASSERT_EQ(sys.moleculeList[0].currpface.size(), 1u);
    EXPECT_EQ(sys.moleculeList[0].currpface.back(), 0) << "partner interface index is 0";
    ASSERT_EQ(sys.moleculeList[0].currprevsep.size(), 1u);
    EXPECT_NEAR(sys.moleculeList[0].currprevsep.back(), 2.0, 1e-9)
        << "recorded separation should be the interface-interface distance";
    ASSERT_EQ(sys.moleculeList[0].currprevnorm.size(), 1u);
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].currprevnorm.back(), 1.0)
        << "with no previous history the reweighting norm must be 1";
    ASSERT_EQ(sys.moleculeList[0].currps_prev.size(), 1u);
    EXPECT_NEAR(sys.moleculeList[0].currps_prev.back(), 1.0 - prob, 1e-12)
        << "stored survival probability must be 1 - p";

    // The partner molecule keeps no curr* history (only proA does).
    EXPECT_EQ(sys.moleculeList[1].currlist.size(), 0u)
        << "higher-indexed molecule should not store the pair history";
}

// -----------------------------------------------------------------------------
// Test 5: zero rate -> placeholder zero probability, no bookkeeping.
// -----------------------------------------------------------------------------
void test_d3bp_zero_rate_gives_zero_probability()
{
    std::cerr << "\n[TEST] test_d3bp_zero_rate_gives_zero_probability\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      interfaces 2 nm apart but rateList[0].rate == 0.\n"
              << "  Pass criteria: the placeholder 0 is appended to both probvecs and\n"
              << "                 remains 0; no reweighting bookkeeping written.\n";

    D3bpFixture sys = d3bp_make_fixture(2.0, 0.0);
    BiMolData biMolData = d3bp_make_bimoldata(20.0, 4.0, 4.0);

    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns);

    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u)
        << "placeholder probability is still appended when inside Rmax";
    ASSERT_EQ(sys.moleculeList[1].probvec.size(), 1u);

    std::cerr << "  probvec[0] = " << sys.moleculeList[0].probvec.back() << '\n';
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].probvec.back(), 0.0)
        << "a zero rate must leave the probability at 0";
    EXPECT_DOUBLE_EQ(sys.moleculeList[1].probvec.back(), 0.0)
        << "a zero rate must leave the partner probability at 0";
    EXPECT_EQ(sys.moleculeList[0].currlist.size(), 0u)
        << "no history recorded for a zero-rate reaction";
}

// -----------------------------------------------------------------------------
// Test 6: a just-dissociated molecule blocks the probability evaluation.
// -----------------------------------------------------------------------------
void test_d3bp_dissociated_molecule_skips_evaluation()
{
    std::cerr << "\n[TEST] test_d3bp_dissociated_molecule_skips_evaluation\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      interfaces 2 nm apart, positive rate, but molecule 1\n"
              << "                 has isDissociated == true.\n"
              << "  Pass criteria: placeholder 0 is appended to both probvecs and stays 0,\n"
              << "                 and no reweighting bookkeeping is written.\n";

    D3bpFixture sys = d3bp_make_fixture(2.0, 100.0);
    sys.moleculeList[1].isDissociated = true; // just came apart this time step
    BiMolData biMolData = d3bp_make_bimoldata(20.0, 4.0, 4.0);

    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns);

    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u)
        << "placeholder probability is appended because the pair is inside Rmax";
    std::cerr << "  probvec[0] = " << sys.moleculeList[0].probvec.back() << '\n';

    EXPECT_DOUBLE_EQ(sys.moleculeList[0].probvec.back(), 0.0)
        << "dissociated partner must prevent rebinding this step";
    EXPECT_DOUBLE_EQ(sys.moleculeList[1].probvec.back(), 0.0)
        << "dissociated partner must prevent rebinding this step";
    EXPECT_EQ(sys.moleculeList[0].currlist.size(), 0u)
        << "no history should be recorded when association is skipped";
}

// -----------------------------------------------------------------------------
// Test 7: OnSurface complexes double kact -> larger probability.
// -----------------------------------------------------------------------------
void test_d3bp_onsurface_doubles_rate()
{
    std::cerr << "\n[TEST] test_d3bp_onsurface_doubles_rate\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      identical geometry/diffusion run twice, once with\n"
              << "                 complex 0 flagged OnSurface (onFiber == false).\n"
              << "  Pass criteria: the OnSurface run yields a strictly larger probability\n"
              << "                 because kact is multiplied by 2.\n";

    // Baseline: no complex on the surface.
    D3bpFixture base = d3bp_make_fixture(2.0, 100.0);
    BiMolData baseData = d3bp_make_bimoldata(20.0, 4.0, 4.0);
    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, baseData, base.params,
        base.moleculeList, base.complexList, base.forwardRxns, base.backRxns);
    ASSERT_EQ(base.moleculeList[0].probvec.size(), 1u);
    const double baseProb = base.moleculeList[0].probvec.back();

    // Surface case: only the OnSurface flag changes so Dtot is identical.
    D3bpFixture surf = d3bp_make_fixture(2.0, 100.0);
    surf.complexList[0].OnSurface = true;
    surf.complexList[0].onFiber = false;
    BiMolData surfData = d3bp_make_bimoldata(20.0, 4.0, 4.0);
    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, surfData, surf.params,
        surf.moleculeList, surf.complexList, surf.forwardRxns, surf.backRxns);
    ASSERT_EQ(surf.moleculeList[0].probvec.size(), 1u);
    const double surfProb = surf.moleculeList[0].probvec.back();

    std::cerr << "  baseline probability = " << baseProb
              << ", OnSurface probability = " << surfProb << '\n';

    EXPECT_GT(surfProb, baseProb)
        << "doubling kact for a surface-bound complex must raise the probability";
    EXPECT_LE(surfProb, 1.0) << "probability must remain physical (<= 1)";
}

// -----------------------------------------------------------------------------
// Test 8: symmetric reaction doubles kact -> larger probability.
// -----------------------------------------------------------------------------
void test_d3bp_symmetric_reaction_doubles_rate()
{
    std::cerr << "\n[TEST] test_d3bp_symmetric_reaction_doubles_rate\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      identical setup run twice, once with\n"
              << "                 forwardRxns[0].isSymmetric == true (A(a)+A(a) case).\n"
              << "  Pass criteria: the symmetric run yields a strictly larger probability.\n";

    D3bpFixture asym = d3bp_make_fixture(2.0, 100.0);
    BiMolData asymData = d3bp_make_bimoldata(20.0, 4.0, 4.0);
    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, asymData, asym.params,
        asym.moleculeList, asym.complexList, asym.forwardRxns, asym.backRxns);
    ASSERT_EQ(asym.moleculeList[0].probvec.size(), 1u);
    const double asymProb = asym.moleculeList[0].probvec.back();

    D3bpFixture sym = d3bp_make_fixture(2.0, 100.0);
    sym.forwardRxns[0].isSymmetric = true;
    BiMolData symData = d3bp_make_bimoldata(20.0, 4.0, 4.0);
    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, symData, sym.params,
        sym.moleculeList, sym.complexList, sym.forwardRxns, sym.backRxns);
    ASSERT_EQ(sym.moleculeList[0].probvec.size(), 1u);
    const double symProb = sym.moleculeList[0].probvec.back();

    std::cerr << "  asymmetric probability = " << asymProb
              << ", symmetric probability = " << symProb << '\n';

    EXPECT_GT(symProb, asymProb)
        << "a symmetric reaction doubles kact and therefore the probability";
    EXPECT_LE(symProb, 1.0) << "probability must remain physical (<= 1)";
}

// -----------------------------------------------------------------------------
// Test 9: previous-encounter history triggers reweighting.
// -----------------------------------------------------------------------------
void test_d3bp_previous_history_reweights_probability()
{
    std::cerr << "\n[TEST] test_d3bp_previous_history_reweights_probability\n"
              << "  Function:      determine_3D_bimolecular_reaction_probability\n"
              << "  Scenario:      molecule 0 already carries a prev* record for the pair\n"
              << "                 (partner 1, own face 0, partner face 0), so the routine\n"
              << "                 must apply pirr_pfree_ratio_psF reweighting.\n"
              << "  Pass criteria: currprevnorm is finite/positive and the stored\n"
              << "                 probability equals (unweighted p) * currprevnorm; the\n"
              << "                 stored survival value is 1 - p.\n";

    // First get the unweighted probability for exactly the same geometry.
    D3bpFixture plain = d3bp_make_fixture(2.0, 100.0);
    BiMolData plainData = d3bp_make_bimoldata(20.0, 4.0, 4.0);
    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, plainData, plain.params,
        plain.moleculeList, plain.complexList, plain.forwardRxns, plain.backRxns);
    ASSERT_EQ(plain.moleculeList[0].probvec.size(), 1u);
    const double unweighted = plain.moleculeList[0].probvec.back();

    // Now the same system, but with a stored previous encounter.
    D3bpFixture hist = d3bp_make_fixture(2.0, 100.0);
    hist.moleculeList[0].prevlist.push_back(1); // partner molecule index
    hist.moleculeList[0].prevmyface.push_back(0); // own interface
    hist.moleculeList[0].prevpface.push_back(0); // partner interface
    hist.moleculeList[0].prevnorm.push_back(1.0); // previous normalisation
    hist.moleculeList[0].ps_prev.push_back(0.95); // previous survival prob.
    hist.moleculeList[0].prevsep.push_back(3.0); // previous separation (nm)

    BiMolData histData = d3bp_make_bimoldata(20.0, 4.0, 4.0);
    determine_3D_bimolecular_reaction_probability(1, 0, 0, false, histData, hist.params,
        hist.moleculeList, hist.complexList, hist.forwardRxns, hist.backRxns);

    ASSERT_EQ(hist.moleculeList[0].probvec.size(), 1u);
    ASSERT_EQ(hist.moleculeList[0].currprevnorm.size(), 1u);

    const double weighted = hist.moleculeList[0].probvec.back();
    const double norm = hist.moleculeList[0].currprevnorm.back();

    std::cerr << "  unweighted p = " << unweighted << ", reweight norm = " << norm
              << ", reweighted p = " << weighted << '\n';

    EXPECT_TRUE(std::isfinite(norm)) << "reweighting norm must be a finite number";
    EXPECT_GT(norm, 0.0) << "reweighting norm must be positive";
    EXPECT_NEAR(weighted, unweighted * norm, 1e-10)
        << "stored probability must be the raw probability times the reweighting norm";
    EXPECT_NEAR(hist.moleculeList[0].currps_prev.back(), 1.0 - weighted, 1e-12)
        << "stored survival probability must be 1 - reweighted p";
    EXPECT_DOUBLE_EQ(hist.moleculeList[1].probvec.back(), weighted)
        << "both partners must share the reweighted probability";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each scenario runs independently so a failure in one does
// not stop the remaining checks.
// -----------------------------------------------------------------------------
TEST(Determine3DBimolecularReactionProbability, RotationalDtotContribution)
{
    test_d3bp_rotational_dtot_contribution();
}
TEST(Determine3DBimolecularReactionProbability, ZeroRotationLeavesDtot)
{
    test_d3bp_zero_rotation_leaves_dtot();
}
TEST(Determine3DBimolecularReactionProbability, OutsideRmaxRecordsNothing)
{
    test_d3bp_outside_rmax_records_nothing();
}
TEST(Determine3DBimolecularReactionProbability, InsideRmaxSetsProbability)
{
    test_d3bp_inside_rmax_sets_probability();
}
TEST(Determine3DBimolecularReactionProbability, ZeroRateGivesZeroProbability)
{
    test_d3bp_zero_rate_gives_zero_probability();
}
TEST(Determine3DBimolecularReactionProbability, DissociatedMoleculeSkipsEvaluation)
{
    test_d3bp_dissociated_molecule_skips_evaluation();
}
TEST(Determine3DBimolecularReactionProbability, OnSurfaceDoublesRate)
{
    test_d3bp_onsurface_doubles_rate();
}
TEST(Determine3DBimolecularReactionProbability, SymmetricReactionDoublesRate)
{
    test_d3bp_symmetric_reaction_doubles_rate();
}
TEST(Determine3DBimolecularReactionProbability, PreviousHistoryReweightsProbability)
{
    test_d3bp_previous_history_reweights_probability();
}