/*! \file test_associate_box.cpp
 *
 * ### Unit test for src/reactions/associate_box.cpp
 *
 * Function under test:
 *
 *     void associate_box(long long int iter, int ifaceIndex1, int ifaceIndex2,
 *                        Molecule& reactMol1, Molecule& reactMol2,
 *                        Complex& reactCom1, Complex& reactCom2,
 *                        const Parameters& params, ForwardRxn& currRxn,
 *                        std::vector<Molecule>& moleculeList,
 *                        std::vector<MolTemplate>& molTemplateList,
 *                        std::map<std::string,int>& observablesList,
 *                        copyCounters& counterArrays,
 *                        std::vector<Complex>& complexList,
 *                        Membrane& membraneObject,
 *                        const std::vector<ForwardRxn>& forwardRxns,
 *                        const std::vector<BackRxn>& backRxns,
 *                        std::ofstream& assocDissocFile);
 *
 * associate_box() performs a physical association event inside a rectangular
 * (box) simulation volume.  Its behaviour splits into three distinct paths that
 * this file exercises:
 *
 *   1. LOOP CLOSURE  - both reacting molecules already belong to the *same*
 *                      complex.  No geometry is touched, counterArrays.nLoops
 *                      is incremented, and only the bond bookkeeping runs.
 *   2. SUCCESSFUL    - two different complexes.  The complexes are translated
 *      ASSOCIATION     towards each other until the two reacting interfaces are
 *                      exactly `bindRadius` apart, all sanity checks pass, and
 *                      the second complex is merged into (and absorbed by) the
 *                      first.
 *   3. CANCELLED     - a sanity check (here: steric overlap of the centres of
 *      ASSOCIATION     mass, controlled by Parameters::overlapSepLimit and
 *                      MolTemplate::checkOverlap) fires, the temporary
 *                      coordinates are discarded and the system is left exactly
 *                      as it was found.
 *
 * All the test systems below use *point* molecules (MolTemplate::isPoint ==
 * true, the single interface sits on the centre of mass) and a ForwardRxn whose
 * association angles are left at their default quiet_NaN() values.  Both of
 * these choices make associate_box() skip every rotation, so the resulting
 * coordinates are completely deterministic and can be predicted analytically:
 *
 *     sigma        = iface1 - iface2
 *     displaceFrac = (|sigma| - bindRadius) / |sigma|
 *     transVec_i   = +/- sigma * (D_i / (D_1 + D_2)) * displaceFrac
 */

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"

namespace {

// Tolerance used when comparing coordinates.  associate_box() only performs
// additions/multiplications on our inputs, so agreement is essentially exact.
constexpr double kAssocBoxTol = 1e-9;

// Absolute interface-state indices used by the test reaction.
constexpr int kReactantAbsIndex = 0; // state of a *free* interface "a"
constexpr int kProductAbsIndex = 2;  // state of a *bound* interface "a"
constexpr int kNumSpeciesStates = 4; // size of the copy-number bookkeeping arrays

/*! \brief Container holding a complete, self consistent mini simulation.
 *
 * Everything associate_box() dereferences must be fully initialised, otherwise
 * the routine walks off the end of a vector.  This struct keeps all of those
 * pieces together so each test can build a fresh, independent system.
 */
struct AssocBoxSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    copyCounters counterArrays {};
    ForwardRxn rxn {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::map<std::string, int> observablesList {};
    std::ofstream assocDissocFile {}; //!< deliberately left closed -> no file output
};

/*! \brief Initialise the global GSL random number generator (deterministic seed).
 *
 * associate_box() only draws a random number on the 1D/fiber path, which these
 * tests do not exercise, but the generator must nevertheless be valid because
 * other translation units share it.
 */
void assocbox_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}

/*! \brief Euclidean distance between two coordinates (test-local helper). */
double assocbox_distance(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Build the single MolTemplate ("A") shared by every molecule.
 *
 * The template describes a point molecule with exactly one interface, no
 * overlap checking, no transition counting and no destruction reactions, which
 * keeps every optional branch of associate_box() switched off.
 */
MolTemplate assocbox_make_template()
{
    MolTemplate temp {};
    temp.molName = "A";
    temp.molTypeIndex = 0;
    temp.copies = 3;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord { 1.0, 1.0, 1.0 };   // equal diffusion => each complex moves half of sigma
    temp.Dr = Coord { 0.01, 0.01, 0.01 };
    temp.isPoint = true;      // -> associate_box skips all rotations
    temp.isRod = false;
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPromoter = false;
    temp.checkOverlap = false; // -> check_for_structure_overlap is a no-op
    temp.countTransition = false;
    temp.canDestroy = false;
    temp.excludeVolumeBound = false;
    temp.insideCompartment = false;
    temp.outsideCompartment = false;
    temp.crossesCompartment = false;

    // A single interface named "a" with a single (unbound) state.
    Interface iface {};
    iface.name = "a";
    iface.index = 0;
    iface.iCoord = Coord { 0.0, 0.0, 0.0 }; // point molecule: interface on the COM
    iface.stateList.emplace_back(std::string("a"), kReactantAbsIndex);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build one point Molecule with a single free interface at \p com. */
Molecule assocbox_make_molecule(int index, int comIndex, const Coord& com)
{
    Molecule mol {};
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = 0;
    mol.myComIndex = comIndex;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.isEmpty = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    Molecule::Iface iface {};
    iface.coord = com; // point molecule: interface coincides with the COM
    iface.index = kReactantAbsIndex;
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.stateIden = '\0';
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    mol.freelist.push_back(0); // interface 0 is currently free

    return mol;
}

/*! \brief Assemble a complete test system.
 *
 * \param[in] sameComplex if true the two reactants share complex 0 (loop
 *                        closure); otherwise each reactant owns its own complex.
 *
 * Molecules 0 and 1 are the reactants (placed at x = -5 and x = +5).  Molecule 2
 * is an uninvolved bystander far away; it is used to verify the "zero out the
 * reaction probabilities of my neighbours" bookkeeping at the end of
 * associate_box().
 */
void assocbox_build_system(AssocBoxSystem& sys, bool sameComplex)
{
    assocbox_init_rng();

    /* ---------------- parameters ---------------- */
    sys.params.numMolTypes = 1;
    sys.params.numTotalSpecies = kNumSpeciesStates;
    sys.params.timeStep = 1.0;
    sys.params.overlapSepLimit = 0.1;   // far smaller than bindRadius -> no overlap cancel
    sys.params.scaleMaxDisplace = 100.0; // very permissive displacement check
    sys.params.nItr = 1;
    Parameters::dt = sys.params.timeStep;

    /* ---------------- boundary (a 100 nm cube) ---------------- */
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    sys.membrane.isBox = true;
    sys.membrane.isSphere = false;
    sys.membrane.hasCompartment = false;
    sys.membrane.implicitLipid = false;
    sys.membrane.TwoD = false;

    /* ---------------- templates ---------------- */
    MolTemplate::numMolTypes = 1;              // MUST be set before building Complexes
    MolTemplate::numEachMolType = std::vector<int> { 3 };
    sys.molTemplateList.clear();
    sys.molTemplateList.push_back(assocbox_make_template());

    /* ---------------- molecules ---------------- */
    sys.moleculeList.clear();
    sys.moleculeList.push_back(assocbox_make_molecule(0, 0, Coord { -5.0, 0.0, 0.0 }));
    sys.moleculeList.push_back(
        assocbox_make_molecule(1, sameComplex ? 0 : 1, Coord { 5.0, 0.0, 0.0 }));
    sys.moleculeList.push_back(assocbox_make_molecule(2, 2, Coord { 30.0, 30.0, 30.0 }));

    /* ---------------- complexes ---------------- */
    // Complex(mol, temp) sets index = mol.index, memberList = {mol.index} and
    // sizes numEachMol / lastNumberUpdateItrEachMol from MolTemplate::numMolTypes.
    sys.complexList.clear();
    sys.complexList.emplace_back(sys.moleculeList[0], sys.molTemplateList[0]);
    sys.complexList.emplace_back(sys.moleculeList[1], sys.molTemplateList[0]);
    sys.complexList.emplace_back(sys.moleculeList[2], sys.molTemplateList[0]);
    for (int i = 0; i < 3; ++i) {
        sys.complexList[i].index = i;
        sys.complexList[i].isEmpty = false;
        sys.complexList[i].OnSurface = false; // solution (3D) complexes
        sys.complexList[i].onFiber = false;
        sys.complexList[i].ncross = 0;
        sys.complexList[i].trajStatus = TrajStatus::none;
    }

    if (sameComplex) {
        // Fold molecule 1 into complex 0 and retire complex 1.
        sys.complexList[0].memberList.push_back(1);
        sys.complexList[0].numEachMol[0] = 2;
        sys.complexList[0].mass = 2.0;
        sys.complexList[1].memberList.clear();
        sys.complexList[1].isEmpty = true;
    }

    /* ---------------- copy counters ---------------- */
    sys.counterArrays.copyNumSpecies.assign(kNumSpeciesStates, 10);
    sys.counterArrays.nBoundPairs.assign(16, 0); // generously sized on purpose
    sys.counterArrays.proPairlist.assign(16, 0);
    sys.counterArrays.singleDouble.assign(kNumSpeciesStates, 0);
    sys.counterArrays.implicitDouble.assign(kNumSpeciesStates, false);
    sys.counterArrays.canDissociate.assign(kNumSpeciesStates, false);
    sys.counterArrays.canDissociate[kProductAbsIndex] = true; // exercise bindPairList
    sys.counterArrays.bindPairList.assign(kNumSpeciesStates, std::vector<int> {});
    sys.counterArrays.bindPairListIL2D.assign(kNumSpeciesStates, std::vector<int> {});
    sys.counterArrays.bindPairListIL3D.assign(kNumSpeciesStates, std::vector<int> {});
    sys.counterArrays.nLoops = 0;
    sys.counterArrays.nAssocSuccess = 0;
    sys.counterArrays.nCancelOverlapPartner = 0;
    sys.counterArrays.nCancelOverlapSystem = 0;
    sys.counterArrays.nCancelSpanBox = 0;
    sys.counterArrays.nCancelDisplace3D = 0;
    sys.counterArrays.nCancelDisplace2D = 0;
    sys.counterArrays.nCancelDisplace3Dto2D = 0;

    // track_association_events() writes into these event histograms.
    init_association_events(sys.counterArrays);
    const auto minSize = static_cast<size_t>(sys.counterArrays.eventArraySize);
    if (sys.counterArrays.events3D.size() < minSize)
        sys.counterArrays.events3D.resize(minSize, 0);
    if (sys.counterArrays.events2D.size() < minSize)
        sys.counterArrays.events2D.resize(minSize, 0);
    if (sys.counterArrays.events3Dto2D.size() < minSize)
        sys.counterArrays.events3Dto2D.resize(minSize, 0);

    /* ---------------- the reaction A(a) + A(a) <-> A(a!1).A(a!1) ---------------- */
    sys.rxn = ForwardRxn {};
    sys.rxn.rxnType = ReactionType::bimolecular;
    sys.rxn.bindRadius = 2.0;
    sys.rxn.isReversible = true;
    sys.rxn.conjBackRxnIndex = 0;
    sys.rxn.absRxnIndex = 0;
    sys.rxn.relRxnIndex = 0;
    sys.rxn.isOnMem = false;
    sys.rxn.isObserved = false;
    sys.rxn.isSymmetric = true;
    // NOTE: assocAngles are left at their default quiet_NaN() so every
    //       theta/phi/omega rotation inside associate_box() is skipped.
    sys.rxn.reactantListNew.emplace_back("a", 0, kReactantAbsIndex, 0, '\0', false);
    sys.rxn.reactantListNew.emplace_back("a", 0, kReactantAbsIndex, 0, '\0', false);
    sys.rxn.productListNew.emplace_back("a", 0, kProductAbsIndex, 0, '\0', true);
    sys.rxn.rateList.emplace_back();
    sys.rxn.rateList.back().rate = 1.0;

    sys.forwardRxns.clear();
    sys.forwardRxns.push_back(sys.rxn);
    sys.backRxns.clear();
}

} // namespace

// -----------------------------------------------------------------------------
// TEST 1: loop closure (both reactants already in the same complex).
// -----------------------------------------------------------------------------
void test_assocbox_loop_closure()
{
    std::cerr << "\n[TEST] test_assocbox_loop_closure\n"
              << "  Source file:   src/reactions/associate_box.cpp\n"
              << "  Function:      associate_box (reactCom1.index == reactCom2.index)\n"
              << "  Scenario:      molecules 0 and 1 are already members of complex 0,\n"
              << "                 so binding them closes a loop.\n"
              << "  Pass criteria: nLoops is incremented, NO coordinates change, both\n"
              << "                 molecules are flagged 'propagated', and the normal\n"
              << "                 bond bookkeeping (isBound/partner/bndlist/freelist,\n"
              << "                 copy numbers) is still performed.\n";

    AssocBoxSystem sys;
    assocbox_build_system(sys, /*sameComplex=*/true);

    const Coord com0Before = sys.moleculeList[0].comCoord;
    const Coord com1Before = sys.moleculeList[1].comCoord;

    std::cerr << "  Calling associate_box(iter=1, iface1=0, iface2=0)...\n";
    associate_box(1, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0],
        sys.complexList[0], sys.params, sys.rxn, sys.moleculeList, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, sys.assocDissocFile);

    // ---- loop counter -------------------------------------------------------
    EXPECT_EQ(sys.counterArrays.nLoops, 1) << "Closing a loop must increment counterArrays.nLoops";
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 0)
        << "The loop-closure branch never reaches the nAssocSuccess counter";

    // ---- geometry untouched -------------------------------------------------
    std::cerr << "  Checking that the loop-closure branch leaves geometry alone\n";
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, com0Before.x, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[0].comCoord.y, com0Before.y, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[0].comCoord.z, com0Before.z, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, com1Before.x, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[1].comCoord.y, com1Before.y, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[1].comCoord.z, com1Before.z, kAssocBoxTol);

    // ---- trajectory status --------------------------------------------------
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "Members of the looping complex must be marked as propagated";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "Members of the looping complex must be marked as propagated";

    // ---- bond bookkeeping still happens -------------------------------------
    std::cerr << "  Checking the interface/bond bookkeeping shared by both branches\n";
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "iface 0 of mol 0 should be bound";
    EXPECT_TRUE(sys.moleculeList[1].interfaceList[0].isBound) << "iface 0 of mol 1 should be bound";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, 1);
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].interaction.partnerIndex, 0);
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIfaceIndex, 0);
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].interaction.partnerIfaceIndex, 0);
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kProductAbsIndex)
        << "Bound interface must adopt the product state index";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].index, kProductAbsIndex);
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.conjBackRxn, sys.rxn.conjBackRxnIndex)
        << "Reversible reaction must record the conjugate back-reaction index";

    EXPECT_EQ(sys.moleculeList[0].bndlist.size(), 1u);
    EXPECT_EQ(sys.moleculeList[1].bndlist.size(), 1u);
    EXPECT_EQ(sys.moleculeList[0].bndpartner.size(), 1u);
    EXPECT_EQ(sys.moleculeList[1].bndpartner.size(), 1u);
    EXPECT_TRUE(sys.moleculeList[0].freelist.empty()) << "The only free interface got bound";
    EXPECT_TRUE(sys.moleculeList[1].freelist.empty()) << "The only free interface got bound";

    // ---- copy numbers -------------------------------------------------------
    std::cerr << "  copyNumSpecies[reactant]=" << sys.counterArrays.copyNumSpecies[kReactantAbsIndex]
              << " copyNumSpecies[product]=" << sys.counterArrays.copyNumSpecies[kProductAbsIndex] << '\n';
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kReactantAbsIndex], 8)
        << "Both reactant interfaces (same abs index) are decremented: 10 - 2 = 8";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kProductAbsIndex], 11)
        << "The product state is incremented once: 10 + 1 = 11";
}

// -----------------------------------------------------------------------------
// TEST 2: successful association of two separate complexes - geometry.
// -----------------------------------------------------------------------------
void test_assocbox_two_complex_geometry()
{
    std::cerr << "\n[TEST] test_assocbox_two_complex_geometry\n"
              << "  Source file:   src/reactions/associate_box.cpp\n"
              << "  Function:      associate_box (two distinct complexes, 3D)\n"
              << "  Scenario:      point molecules at x=-5 and x=+5, bindRadius = 2 nm,\n"
              << "                 equal diffusion constants, all association angles NaN\n"
              << "                 so no rotation is applied.\n"
              << "  Pass criteria: each complex is translated by half of the separation\n"
              << "                 deficit, i.e. mol0 -> x=-1, mol1 -> x=+1, giving an\n"
              << "                 interface-interface distance of exactly bindRadius.\n";

    AssocBoxSystem sys;
    assocbox_build_system(sys, /*sameComplex=*/false);

    const double bindRadius = sys.rxn.bindRadius;
    const double sepBefore = assocbox_distance(
        sys.moleculeList[0].interfaceList[0].coord, sys.moleculeList[1].interfaceList[0].coord);
    std::cerr << "  Interface separation before = " << sepBefore << " nm (bindRadius = " << bindRadius << ")\n";

    std::cerr << "  Calling associate_box(iter=1, iface1=0, iface2=0)...\n";
    associate_box(1, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0],
        sys.complexList[1], sys.params, sys.rxn, sys.moleculeList, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, sys.assocDissocFile);

    const double sepAfter = assocbox_distance(
        sys.moleculeList[0].interfaceList[0].coord, sys.moleculeList[1].interfaceList[0].coord);
    std::cerr << "  Interface separation after  = " << sepAfter << " nm\n";
    std::cerr << "  mol0 COM = (" << sys.moleculeList[0].comCoord.x << ", "
              << sys.moleculeList[0].comCoord.y << ", " << sys.moleculeList[0].comCoord.z << ")\n";
    std::cerr << "  mol1 COM = (" << sys.moleculeList[1].comCoord.x << ", "
              << sys.moleculeList[1].comCoord.y << ", " << sys.moleculeList[1].comCoord.z << ")\n";

    // The association must have succeeded (otherwise nothing moved).
    ASSERT_EQ(sys.counterArrays.nAssocSuccess, 1)
        << "Association was unexpectedly cancelled - the geometry checks below are meaningless";

    EXPECT_NEAR(sepAfter, bindRadius, 1e-9)
        << "After association the two reacting interfaces must sit exactly at sigma";

    // Symmetric diffusion constants -> each complex covers half of the deficit.
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, -1.0, kAssocBoxTol) << "mol0 should move -5 -> -1";
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, 1.0, kAssocBoxTol) << "mol1 should move +5 -> +1";
    // sigma has no y/z component, so those coordinates must be untouched.
    EXPECT_NEAR(sys.moleculeList[0].comCoord.y, 0.0, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[0].comCoord.z, 0.0, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[1].comCoord.y, 0.0, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[1].comCoord.z, 0.0, kAssocBoxTol);

    // Point molecules: the interface coordinate follows the centre of mass.
    EXPECT_NEAR(sys.moleculeList[0].interfaceList[0].coord.x, sys.moleculeList[0].comCoord.x, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[1].interfaceList[0].coord.x, sys.moleculeList[1].comCoord.x, kAssocBoxTol);

    // The uninvolved bystander must not have moved at all.
    EXPECT_NEAR(sys.moleculeList[2].comCoord.x, 30.0, kAssocBoxTol) << "bystander molecule must not move";
    EXPECT_NEAR(sys.moleculeList[2].comCoord.y, 30.0, kAssocBoxTol);
    EXPECT_NEAR(sys.moleculeList[2].comCoord.z, 30.0, kAssocBoxTol);

    // Temporary association coordinates must be released again.
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty()) << "tmpICoords must be cleared after association";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty()) << "tmpICoords must be cleared after association";
}

// -----------------------------------------------------------------------------
// TEST 3: successful association - complex merging, counters and neighbour lists.
// -----------------------------------------------------------------------------
void test_assocbox_two_complex_bookkeeping()
{
    std::cerr << "\n[TEST] test_assocbox_two_complex_bookkeeping\n"
              << "  Source file:   src/reactions/associate_box.cpp\n"
              << "  Function:      associate_box (two distinct complexes, 3D)\n"
              << "  Scenario:      same system as the geometry test, but now we inspect\n"
              << "                 the bookkeeping: complex merge, copy numbers, bound\n"
              << "                 pair counters, bindPairList and neighbour probvec.\n"
              << "  Pass criteria: complex 1 is destroyed, its member joins complex 0,\n"
              << "                 nAssocSuccess == 1, copy numbers shift by -2/+1, the\n"
              << "                 dissociable product is registered in bindPairList, and\n"
              << "                 the reacting molecules stop avoiding their neighbours.\n";

    AssocBoxSystem sys;
    assocbox_build_system(sys, /*sameComplex=*/false);

    // Set up a mutual "encounter" between molecule 0 and the bystander molecule 2
    // so that we can verify the probability zeroing at the end of associate_box().
    sys.moleculeList[0].crossbase.push_back(2);
    sys.moleculeList[0].mycrossint.push_back(0);
    sys.moleculeList[0].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    sys.moleculeList[0].probvec.push_back(0.5);

    sys.moleculeList[2].crossbase.push_back(0);
    sys.moleculeList[2].mycrossint.push_back(0);
    sys.moleculeList[2].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    sys.moleculeList[2].probvec.push_back(0.5);

    std::cerr << "  Calling associate_box(iter=7, iface1=0, iface2=0)...\n";
    associate_box(7, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0],
        sys.complexList[1], sys.params, sys.rxn, sys.moleculeList, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, sys.assocDissocFile);

    // ---- success counter ----------------------------------------------------
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "A successful association must be counted";
    EXPECT_EQ(sys.counterArrays.nLoops, 0) << "This is not a loop closure";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapPartner, 0) << "No overlap cancellation expected";
    EXPECT_EQ(sys.counterArrays.nCancelSpanBox, 0) << "The dimer is far smaller than the box";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapSystem, 0) << "No system overlap expected";
    EXPECT_EQ(sys.counterArrays.nCancelDisplace3D, 0) << "The 2 nm displacement is well below the limit";

    // ---- complex merge ------------------------------------------------------
    std::cerr << "  complexList[0].memberList.size() = " << sys.complexList[0].memberList.size()
              << ", complexList[1].isEmpty = " << std::boolalpha << sys.complexList[1].isEmpty << '\n';
    EXPECT_EQ(sys.complexList[0].memberList.size(), 2u)
        << "The surviving complex must contain both molecules";
    EXPECT_TRUE(sys.complexList[1].isEmpty) << "The absorbed complex must be flagged as empty";
    EXPECT_TRUE(sys.complexList[1].memberList.empty()) << "The absorbed complex must have no members";
    EXPECT_EQ(sys.moleculeList[1].myComIndex, sys.complexList[0].index)
        << "Molecule 1 must now point at the surviving complex";
    EXPECT_EQ(sys.moleculeList[0].myComIndex, sys.complexList[0].index);
    EXPECT_EQ(sys.complexList[0].numEachMol[0], 2)
        << "update_properties() must recount the molecules of type 0";
    EXPECT_NEAR(sys.complexList[0].mass, 2.0, 1e-9) << "The merged complex has the summed mass";

    // ---- trajectory status --------------------------------------------------
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated));
    EXPECT_EQ(static_cast<int>(sys.moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated));

    // ---- copy numbers / bound pairs ----------------------------------------
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kReactantAbsIndex], 8)
        << "Two free interfaces of the same species are consumed";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kProductAbsIndex], 11)
        << "One bound (product) species is created";
    const int boundPairSum = std::accumulate(
        sys.counterArrays.nBoundPairs.begin(), sys.counterArrays.nBoundPairs.end(), 0);
    std::cerr << "  Sum over nBoundPairs = " << boundPairSum << '\n';
    EXPECT_GE(boundPairSum, 1) << "update_Nboundpairs() must record the new A-A bond";

    // ---- dissociation pool --------------------------------------------------
    ASSERT_EQ(sys.counterArrays.bindPairList.size(), static_cast<size_t>(kNumSpeciesStates));
    ASSERT_EQ(sys.counterArrays.bindPairList[kProductAbsIndex].size(), 1u)
        << "canDissociate[product] is true, so the new bond must join the dissociation pool";
    EXPECT_EQ(sys.counterArrays.bindPairList[kProductAbsIndex][0], sys.moleculeList[0].index)
        << "The molecule owning reactantListNew[0].relIfaceIndex is stored";

    // ---- neighbour / crossing bookkeeping ----------------------------------
    std::cerr << "  Checking that neighbour reaction probabilities were zeroed\n";
    ASSERT_EQ(sys.moleculeList[2].probvec.size(), 1u);
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[0], 0.0)
        << "Neighbours of a reacted molecule must have their probability set to zero";
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty()) << "reactMol1.crossbase is cleared";
    EXPECT_TRUE(sys.moleculeList[1].crossbase.empty()) << "reactMol2.crossbase is cleared";
    EXPECT_EQ(sys.complexList[0].ncross, -1) << "The product complex is marked with ncross = -1";

    // ---- interface state ----------------------------------------------------
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound);
    EXPECT_TRUE(sys.moleculeList[1].interfaceList[0].isBound);
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, 1);
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].interaction.partnerIndex, 0);
    EXPECT_TRUE(sys.moleculeList[0].freelist.empty());
    EXPECT_TRUE(sys.moleculeList[1].freelist.empty());
}

// -----------------------------------------------------------------------------
// TEST 4: cancelled association because of a steric (COM-COM) overlap.
// -----------------------------------------------------------------------------
void test_assocbox_cancelled_by_overlap()
{
    std::cerr << "\n[TEST] test_assocbox_cancelled_by_overlap\n"
              << "  Source file:   src/reactions/associate_box.cpp\n"
              << "  Function:      associate_box -> check_for_structure_overlap\n"
              << "  Scenario:      MolTemplate::checkOverlap is switched on and\n"
              << "                 Parameters::overlapSepLimit (5 nm) is made larger than\n"
              << "                 the binding radius (2 nm), so the two point molecules\n"
              << "                 necessarily overlap once pushed to sigma.\n"
              << "  Pass criteria: association is cancelled -> nCancelOverlapPartner == 1,\n"
              << "                 nAssocSuccess == 0, coordinates/complexes/copy numbers\n"
              << "                 are all restored to their pre-call values and the\n"
              << "                 temporary association coordinates are released.\n";

    AssocBoxSystem sys;
    assocbox_build_system(sys, /*sameComplex=*/false);

    // Turn the overlap check on and make it impossible to satisfy.
    sys.molTemplateList[0].checkOverlap = true;
    sys.params.overlapSepLimit = 5.0; // > bindRadius (2.0)

    const Coord com0Before = sys.moleculeList[0].comCoord;
    const Coord com1Before = sys.moleculeList[1].comCoord;

    std::cerr << "  Calling associate_box(iter=3, iface1=0, iface2=0)...\n";
    associate_box(3, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0],
        sys.complexList[1], sys.params, sys.rxn, sys.moleculeList, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, sys.assocDissocFile);

    std::cerr << "  nCancelOverlapPartner = " << sys.counterArrays.nCancelOverlapPartner
              << ", nAssocSuccess = " << sys.counterArrays.nAssocSuccess << '\n';

    // ---- cancellation counters ---------------------------------------------
    EXPECT_EQ(sys.counterArrays.nCancelOverlapPartner, 1)
        << "A partner-overlap cancellation must be recorded";
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 0) << "No successful association may be recorded";

    // ---- system fully restored ---------------------------------------------
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, com0Before.x, kAssocBoxTol)
        << "Cancelled association must leave the real coordinates untouched";
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, com1Before.x, kAssocBoxTol)
        << "Cancelled association must leave the real coordinates untouched";
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty())
        << "Temporary association coordinates must be cleared on cancellation";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty())
        << "Temporary association coordinates must be cleared on cancellation";

    // ---- no bond formed -----------------------------------------------------
    EXPECT_FALSE(sys.moleculeList[0].interfaceList[0].isBound) << "No bond may be formed";
    EXPECT_FALSE(sys.moleculeList[1].interfaceList[0].isBound) << "No bond may be formed";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kReactantAbsIndex)
        << "The interface must keep its free-state index";
    EXPECT_TRUE(sys.moleculeList[0].bndlist.empty());
    EXPECT_TRUE(sys.moleculeList[1].bndpartner.empty());
    EXPECT_EQ(sys.moleculeList[0].freelist.size(), 1u) << "The interface is still free";
    EXPECT_EQ(sys.moleculeList[1].freelist.size(), 1u) << "The interface is still free";

    // ---- complexes untouched ------------------------------------------------
    EXPECT_EQ(sys.complexList[0].memberList.size(), 1u) << "Complexes must not be merged";
    EXPECT_EQ(sys.complexList[1].memberList.size(), 1u) << "Complexes must not be merged";
    EXPECT_FALSE(sys.complexList[1].isEmpty) << "The second complex must survive";
    EXPECT_EQ(sys.moleculeList[1].myComIndex, 1) << "Molecule 1 keeps its own complex";

    // ---- copy numbers untouched --------------------------------------------
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kReactantAbsIndex], 10)
        << "Copy numbers may not change on a cancelled association";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kProductAbsIndex], 10)
        << "Copy numbers may not change on a cancelled association";
    EXPECT_TRUE(sys.counterArrays.bindPairList[kProductAbsIndex].empty())
        << "No entry may be added to the dissociation pool";
}

// -----------------------------------------------------------------------------
// TEST 5: observable tracking for an observed reaction product.
// -----------------------------------------------------------------------------
void test_assocbox_observable_tracking()
{
    std::cerr << "\n[TEST] test_assocbox_observable_tracking\n"
              << "  Source file:   src/reactions/associate_box.cpp\n"
              << "  Function:      associate_box (observable increment tail)\n"
              << "  Scenario:      ForwardRxn::isObserved is true with label \"AA\" and the\n"
              << "                 observables map already carries that label.  A loop\n"
              << "                 closure is used so the result is fully deterministic.\n"
              << "  Pass criteria: observablesList[\"AA\"] is incremented by one; an unknown\n"
              << "                 label leaves the map untouched.\n";

    // --- case A: the label exists in the map --------------------------------
    {
        AssocBoxSystem sys;
        assocbox_build_system(sys, /*sameComplex=*/true);
        sys.rxn.isObserved = true;
        sys.rxn.observeLabel = "AA";
        sys.observablesList["AA"] = 4;

        std::cerr << "  [A] observable \"AA\" starts at " << sys.observablesList["AA"] << '\n';
        associate_box(1, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0],
            sys.complexList[0], sys.params, sys.rxn, sys.moleculeList, sys.molTemplateList,
            sys.observablesList, sys.counterArrays, sys.complexList, sys.membrane,
            sys.forwardRxns, sys.backRxns, sys.assocDissocFile);
        std::cerr << "  [A] observable \"AA\" ends at   " << sys.observablesList["AA"] << '\n';

        EXPECT_EQ(sys.observablesList["AA"], 5)
            << "An observed reaction must increment its observable counter";
        EXPECT_EQ(sys.observablesList.size(), 1u) << "No new observable entry may be created";
    }

    // --- case B: the label is NOT in the map --------------------------------
    {
        AssocBoxSystem sys;
        assocbox_build_system(sys, /*sameComplex=*/true);
        sys.rxn.isObserved = true;
        sys.rxn.observeLabel = "NOT_TRACKED";
        sys.observablesList["AA"] = 4;

        std::cerr << "  [B] reaction observes an untracked label; map must be unchanged\n";
        associate_box(1, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0],
            sys.complexList[0], sys.params, sys.rxn, sys.moleculeList, sys.molTemplateList,
            sys.observablesList, sys.counterArrays, sys.complexList, sys.membrane,
            sys.forwardRxns, sys.backRxns, sys.assocDissocFile);

        EXPECT_EQ(sys.observablesList.size(), 1u)
            << "An unknown observable label must not insert a new map entry";
        EXPECT_EQ(sys.observablesList["AA"], 4)
            << "An unrelated observable must not be modified";
    }
}

// -----------------------------------------------------------------------------
// TEST 6: an irreversible reaction must not store a conjugate back reaction.
// -----------------------------------------------------------------------------
void test_assocbox_irreversible_no_back_rxn()
{
    std::cerr << "\n[TEST] test_assocbox_irreversible_no_back_rxn\n"
              << "  Source file:   src/reactions/associate_box.cpp\n"
              << "  Function:      associate_box (conjBackRxn assignment)\n"
              << "  Scenario:      the same loop-closure system but with\n"
              << "                 ForwardRxn::isReversible == false.\n"
              << "  Pass criteria: the bond is still formed, but Interaction::conjBackRxn\n"
              << "                 keeps its default value of -1 for both partners.\n";

    AssocBoxSystem sys;
    assocbox_build_system(sys, /*sameComplex=*/true);
    sys.rxn.isReversible = false;
    sys.rxn.conjBackRxnIndex = -1;

    std::cerr << "  Calling associate_box with an irreversible reaction...\n";
    associate_box(1, 0, 0, sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0],
        sys.complexList[0], sys.params, sys.rxn, sys.moleculeList, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, sys.assocDissocFile);

    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "The bond is still formed";
    EXPECT_TRUE(sys.moleculeList[1].interfaceList[0].isBound) << "The bond is still formed";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.conjBackRxn, -1)
        << "Irreversible reactions leave conjBackRxn at its default (-1)";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].interaction.conjBackRxn, -1)
        << "Irreversible reactions leave conjBackRxn at its default (-1)";
    EXPECT_EQ(sys.moleculeList[0].bndRxnList.size(), 1u)
        << "The forward reaction index is recorded in bndRxnList regardless";
    EXPECT_EQ(sys.moleculeList[0].bndRxnList[0], sys.rxn.relRxnIndex);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario is run in its own TEST so that a failure
// in one path does not stop the remaining paths from executing.
// -----------------------------------------------------------------------------
TEST(AssociateBoxTest, LoopClosure) { test_assocbox_loop_closure(); }
TEST(AssociateBoxTest, TwoComplexGeometry) { test_assocbox_two_complex_geometry(); }
TEST(AssociateBoxTest, TwoComplexBookkeeping) { test_assocbox_two_complex_bookkeeping(); }
TEST(AssociateBoxTest, CancelledByOverlap) { test_assocbox_cancelled_by_overlap(); }
TEST(AssociateBoxTest, ObservableTracking) { test_assocbox_observable_tracking(); }
TEST(AssociateBoxTest, IrreversibleNoBackRxn) { test_assocbox_irreversible_no_back_rxn(); }