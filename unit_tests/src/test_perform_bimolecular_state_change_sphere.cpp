/*! \file test_perform_bimolecular_state_change_sphere.cpp
 *
 * ### Unit test for src/reactions/perform_bimolecular_state_change_sphere.cpp
 *
 * The function under test is
 *
 *     void perform_bimolecular_state_change_sphere(
 *              int stateChangeIface, int facilitatorIface, std::array<int,3>& rxnItr,
 *              Molecule& stateChangeMol, Molecule& facilitatorMol,
 *              Complex& stateChangeCom, Complex& facilitatorCom,
 *              copyCounters&, const Parameters&,
 *              std::vector<ForwardRxn>&, std::vector<BackRxn>&,
 *              std::vector<Molecule>&, std::vector<Complex>&,
 *              std::vector<MolTemplate>&, std::map<std::string,int>&, Membrane&)
 *
 * It performs a *bimolecular state change* inside a spherical simulation
 * volume: a "facilitator" molecule is brought into contact with a
 * "state change" molecule, the pair is oriented, checked for overlap, and -
 * if nothing cancels the event - the state of one interface is flipped and
 * the species/observable counters are updated.
 *
 * The routine touches a great deal of machinery (association geometry,
 * boundary reflection, overlap checks, counters, observables).  To keep the
 * test deterministic and to keep the focus on *this* file, every molecule
 * type used here is declared `isPoint = true`, which makes the routine skip
 * all the theta/phi/omega rotations (that branch is exercised by the
 * association tests).  What remains and is verified here is:
 *
 *   1. the "move to sigma" displacement of the two complexes,
 *   2. the forcing of the pair centre-of-mass back to its starting value,
 *   3. the actual interface state change (forward and back reaction),
 *   4. the species copy-number bookkeeping,
 *   5. the observable increment/decrement,
 *   6. the zeroing of the reaction probabilities of crossing partners,
 *   7. the bookkeeping performed when the event is cancelled by overlap,
 *   8. the 2D (both complexes on the sphere surface) code path.
 *
 * NOTE on the geometry assertions: in this source file the displacement
 * vector built from -sigma (`transVec1`) is applied to `stateChangeCom`
 * while the vector built from +sigma (`transVec2`) is applied to
 * `facilitatorCom`, where `sigma = reactIface1(facilitator) -
 * reactIface2(stateChange)`.  The two complexes therefore move *apart*
 * rather than to contact.  The tests below assert the behaviour the code
 * actually has, and the expected numbers are derived from that code.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Small constants / helpers used by all of the tests in this file.
// ---------------------------------------------------------------------------

const double kBscPi = 3.14159265358979323846;

//! Radius of the spherical boundary used for the 3D (solution) tests.
const double kBscSolutionSphereR = 500.0;
//! Radius of the spherical boundary used for the 2D (membrane) test.
const double kBscMembraneSphereR = 100.0;

/*! \brief Euclidean distance between two coordinates. */
double bsc_dist(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Magnitude of a coordinate treated as a vector from the origin. */
double bsc_mag(const Coord& a) { return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }

/*! \brief Make sure the global GSL generator exists.
 *
 * Some of the boundary-condition helpers called from the routine under test
 * may resample a trajectory, which uses the global `r`.  `gtest_main.cpp`
 * defines `r` but leaves it null, so allocate it here if nobody else has.
 */
void bsc_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}

// ---------------------------------------------------------------------------
// Builders for the minimal system the routine needs.
// ---------------------------------------------------------------------------

/*! \brief Build a "point" molecule template (all interfaces sit on the COM).
 *
 * Point templates make perform_bimolecular_state_change_sphere() skip every
 * orientation rotation, which keeps the geometry analytically predictable.
 */
MolTemplate bsc_make_template(int typeIndex, const std::string& name, const std::string& ifaceName,
    const std::vector<std::pair<char, int>>& states, const Coord& D, bool checkOverlap)
{
    MolTemplate mt;
    mt.molTypeIndex = typeIndex;
    mt.molName = name;
    mt.copies = 1;
    mt.mass = 1.0;
    mt.radius = 1.0;
    mt.isPoint = true; // <- skips theta/phi/omega rotations
    mt.isRod = false;
    mt.isLipid = false;
    mt.isImplicitLipid = false;
    mt.isPromoter = false;
    mt.checkOverlap = checkOverlap;
    mt.comCoord = Coord(0.0, 0.0, 0.0);
    mt.D = D;
    mt.Dr = Coord(0.0, 0.0, 0.0);

    Interface iface;
    iface.index = 0;
    iface.name = ifaceName;
    iface.iCoord = Coord(0.0, 0.0, 0.0); // point molecule: interface == COM
    for (const auto& s : states)
        iface.stateList.emplace_back(s.first, s.second); // State(char iden, int index)
    iface.set_ifaceAndStateNames();
    mt.interfaceList.push_back(iface);

    return mt;
}

/*! \brief Build a single-interface molecule sitting at \p pos. */
Molecule bsc_make_molecule(int index, int comIndex, int molTypeIndex, const Coord& pos, int absIfaceIndex,
    char stateIden, int stateIndex)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = pos;

    Molecule::Iface iface;
    iface.coord = pos; // point molecule
    iface.index = absIfaceIndex;
    iface.relIndex = 0;
    iface.stateIden = stateIden;
    iface.stateIndex = stateIndex;
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    mol.freelist.push_back(0);
    return mol;
}

/*! \brief Build a complex holding exactly one member molecule. */
Complex bsc_make_complex(int index, int memberMol, int memberType, const Coord& pos, const Coord& D, bool onSurface)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.comCoord = pos;
    com.memberList.push_back(memberMol);
    com.numEachMol = std::vector<int>(2, 0);
    com.numEachMol[memberType] = 1;
    com.lastNumberUpdateItrEachMol.resize(2);
    com.D = D;
    com.Dr = Coord(0.0, 0.0, 0.0);
    com.mass = 1.0;
    com.radius = 1.0;
    com.isEmpty = false;
    com.OnSurface = onSurface;
    com.onFiber = false;
    com.linksToSurface = 0;
    com.ncross = 2; // non-trivial value, the routine must set it to -1
    com.trajStatus = TrajStatus::none;
    return com;
}

/*! \brief Convenience wrapper around the RxnIface constructor. */
RxnIface bsc_make_rxn_iface(const std::string& name, int molTypeIndex, int absIndex, int relIndex, char state)
{
    return RxnIface(name, molTypeIndex, absIndex, relIndex, state, false);
}

/*! \struct BscSphereSystem
 * \brief All of the containers the routine under test needs, in one place.
 */
struct BscSphereSystem {
    Parameters params;
    Membrane membrane;
    std::vector<MolTemplate> molTemplateList;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::map<std::string, int> observablesList;
    copyCounters counterArrays;
};

/*! \brief Assemble the whole test system.
 *
 * Layout (indices):
 *   molecule 0 / complex 0 : molecule "A", type 0, the *state change* molecule
 *   molecule 1 / complex 1 : molecule "B", type 1, the *facilitator*
 *   molecule 2 / complex 2 : a bystander that "crosses" both reactants; it is
 *                            used to check that the routine zeroes the
 *                            reaction probabilities of crossing partners.
 *
 * \param[in] onMembrane   place both complexes on the sphere surface with
 *                         D.z == 0 so the 2D (membrane) branch is taken.
 * \param[in] checkOverlap flag both templates for overlap checking (used to
 *                         force the association to be cancelled).
 * \param[in] startState   starting state identity of A's interface ('U' or 'P').
 */
BscSphereSystem bsc_build_system(bool onMembrane, bool checkOverlap, char startState)
{
    BscSphereSystem sys;

    // ---------------- parameters ----------------
    sys.params.timeStep = 1.0;
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 8;
    sys.params.rMaxLimit = 50.0;
    sys.params.scaleMaxDisplace = 1.0e6; // never cancel because of displacement
    // A huge overlap limit makes check_for_structure_overlap() cancel the event.
    sys.params.overlapSepLimit = checkOverlap ? 100.0 : 0.1;

    // ---------------- spherical boundary ----------------
    const double sphereR = onMembrane ? kBscMembraneSphereR : kBscSolutionSphereR;
    sys.membrane.isSphere = true;
    sys.membrane.isBox = false;
    sys.membrane.hasCompartment = false;
    sys.membrane.implicitLipid = false;
    sys.membrane.sphereR = sphereR;
    sys.membrane.sphereVol = (4.0 / 3.0) * kBscPi * sphereR * sphereR * sphereR;
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double>{ 2 * sphereR, 2 * sphereR, 2 * sphereR });
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";

    // ---------------- diffusion constants ----------------
    // In the 2D case both complexes must have D.z == 0 so that DzSum < 1e-14.
    const Coord dState = onMembrane ? Coord(0.5, 0.5, 0.0) : Coord(1.0, 1.0, 1.0);
    const Coord dFacil = onMembrane ? Coord(0.5, 0.5, 0.0) : Coord(3.0, 3.0, 3.0);

    // ---------------- molecule templates ----------------
    // Type 0 ("A") carries two states: 'U' (absolute index 0) and 'P' (index 1).
    sys.molTemplateList.push_back(
        bsc_make_template(0, "A", "a", { { 'U', 0 }, { 'P', 1 } }, dState, checkOverlap));
    // Type 1 ("B") is the facilitator, single stateless interface (index 2).
    sys.molTemplateList.push_back(
        bsc_make_template(1, "B", "b", { { '\0', 2 } }, dFacil, checkOverlap));

    // ---------------- coordinates ----------------
    Coord posA;
    Coord posB;
    Coord posBystander;
    if (onMembrane) {
        // Both reactants sit on the sphere, separated by a 0.1 rad geodesic
        // (arc length = 10 nm for R = 100 nm).
        posA = Coord(sphereR, 0.0, 0.0);
        posB = Coord(sphereR * std::cos(0.1), sphereR * std::sin(0.1), 0.0);
        posBystander = Coord(0.0, 0.0, sphereR);
    } else {
        posA = Coord(0.0, 0.0, 0.0);
        posB = Coord(10.0, 0.0, 0.0);
        posBystander = Coord(0.0, 0.0, 40.0);
    }

    // ---------------- molecules ----------------
    const int startAbsIndex = (startState == 'U') ? 0 : 1;
    const int startRelIndex = (startState == 'U') ? 0 : 1;
    sys.moleculeList.push_back(bsc_make_molecule(0, 0, 0, posA, startAbsIndex, startState, startRelIndex));
    sys.moleculeList.push_back(bsc_make_molecule(1, 1, 1, posB, 2, '\0', 0));
    sys.moleculeList.push_back(bsc_make_molecule(2, 2, 1, posBystander, 2, '\0', 0));

    // Cross lists: the bystander (index 2) has both reactants in its crossbase,
    // and both reactants list the bystander.  The routine must set the
    // bystander's probabilities for those two partners to zero.
    sys.moleculeList[0].crossbase.push_back(2);
    sys.moleculeList[0].mycrossint.push_back(0);
    sys.moleculeList[0].crossrxn.push_back(std::array<int, 3>{ 0, 0, 0 });
    sys.moleculeList[0].probvec.push_back(0.25);

    sys.moleculeList[1].crossbase.push_back(2);
    sys.moleculeList[1].mycrossint.push_back(0);
    sys.moleculeList[1].crossrxn.push_back(std::array<int, 3>{ 0, 0, 0 });
    sys.moleculeList[1].probvec.push_back(0.75);

    sys.moleculeList[2].crossbase.push_back(0);
    sys.moleculeList[2].crossbase.push_back(1);
    sys.moleculeList[2].mycrossint.push_back(0);
    sys.moleculeList[2].mycrossint.push_back(0);
    sys.moleculeList[2].crossrxn.push_back(std::array<int, 3>{ 0, 0, 0 });
    sys.moleculeList[2].crossrxn.push_back(std::array<int, 3>{ 0, 0, 0 });
    sys.moleculeList[2].probvec.push_back(0.25);
    sys.moleculeList[2].probvec.push_back(0.75);

    // ---------------- complexes ----------------
    sys.complexList.push_back(bsc_make_complex(0, 0, 0, posA, dState, onMembrane));
    sys.complexList.push_back(bsc_make_complex(1, 1, 1, posB, dFacil, onMembrane));
    sys.complexList.push_back(bsc_make_complex(2, 2, 1, posBystander, dFacil, onMembrane));

    // ---------------- reactions ----------------
    // Forward: B(b) + A(a~U) -> B(b) + A(a~P)
    ForwardRxn fr;
    fr.rxnType = ReactionType::biMolStateChange;
    fr.absRxnIndex = 0;
    fr.relRxnIndex = 0;
    fr.hasStateChange = true;
    fr.isReversible = true;
    fr.conjBackRxnIndex = 0;
    fr.bindRadius = 1.0;
    fr.bindRadius2D = 1.0;
    fr.isObserved = true;
    fr.observeLabel = "A_P";
    fr.reactantListNew.push_back(bsc_make_rxn_iface("b", 1, 2, 0, '\0'));
    fr.reactantListNew.push_back(bsc_make_rxn_iface("a", 0, 0, 0, 'U'));
    fr.productListNew.push_back(bsc_make_rxn_iface("b", 1, 2, 0, '\0'));
    fr.productListNew.push_back(bsc_make_rxn_iface("a", 0, 1, 0, 'P')); // <- used as newState
    fr.stateChangeIface = std::pair<RxnIface, RxnIface>(fr.reactantListNew[1], fr.productListNew[1]);
    fr.rateList.emplace_back();
    fr.rateList.back().rate = 1.0;
    sys.forwardRxns.push_back(fr);

    // Back: B(b) + A(a~P) -> B(b) + A(a~U)
    BackRxn br;
    br.rxnType = ReactionType::biMolStateChange;
    br.absRxnIndex = 1;
    br.relRxnIndex = 0;
    br.hasStateChange = true;
    br.conjForwardRxnIndex = 0;
    br.isObserved = true;
    br.observeLabel = "A_U";
    br.reactantListNew = fr.productListNew;
    br.productListNew.push_back(bsc_make_rxn_iface("b", 1, 2, 0, '\0'));
    br.productListNew.push_back(bsc_make_rxn_iface("a", 0, 0, 0, 'U')); // <- used as newState
    br.stateChangeIface = std::pair<RxnIface, RxnIface>(fr.productListNew[1], fr.reactantListNew[1]);
    br.rateList.emplace_back();
    br.rateList.back().rate = 1.0;
    sys.backRxns.push_back(br);

    // ---------------- counters / observables ----------------
    sys.counterArrays.copyNumSpecies = std::vector<int>(8, 50);
    sys.observablesList["A_P"] = 5;
    sys.observablesList["A_U"] = 7;

    return sys;
}

/*! \brief Invoke the routine under test on molecules 0 (state change) and 1 (facilitator). */
void bsc_run(BscSphereSystem& sys, std::array<int, 3>& rxnItr)
{
    perform_bimolecular_state_change_sphere(/*stateChangeIface=*/0, /*facilitatorIface=*/0, rxnItr,
        sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0], sys.complexList[1], sys.counterArrays,
        sys.params, sys.forwardRxns, sys.backRxns, sys.moleculeList, sys.complexList, sys.molTemplateList,
        sys.observablesList, sys.membrane);
}

} // namespace

// ---------------------------------------------------------------------------
// Test 1: forward reaction in solution (both complexes 3D).
// ---------------------------------------------------------------------------
void test_bsc_sphere_forward_in_solution()
{
    std::cerr << "\n[TEST] test_bsc_sphere_forward_in_solution\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_sphere.cpp\n"
              << "  Function:    perform_bimolecular_state_change_sphere (3D branch)\n"
              << "  Scenario:    A(a~U) at (0,0,0) with D=(1,1,1) and facilitator B at\n"
              << "               (10,0,0) with D=(3,3,3) inside a sphere of radius 500.\n"
              << "  Criteria:    A's interface state flips U->P, copy numbers and the\n"
              << "               observable are updated, both complexes are displaced by\n"
              << "               the amounts the source computes and the pair centre of\n"
              << "               mass is restored to its starting value.\n";

    bsc_init_rng();
    BscSphereSystem sys = bsc_build_system(/*onMembrane=*/false, /*checkOverlap=*/false, /*startState=*/'U');

    // Book-keeping of the pre-reaction state so we can compare afterwards.
    const Coord startA = sys.moleculeList[0].comCoord;
    const Coord startB = sys.moleculeList[1].comCoord;
    const double startPairComX = 0.5 * (startA.x + startB.x); // equal masses -> plain average
    const int startCopyU = sys.counterArrays.copyNumSpecies[0];
    const int startCopyP = sys.counterArrays.copyNumSpecies[1];
    const int startObs = sys.observablesList["A_P"];

    std::cerr << "  Initial separation = " << bsc_dist(startA, startB) << " nm, bindRadius = "
              << sys.forwardRxns[0].bindRadius << " nm\n";

    std::array<int, 3> rxnItr { 0, 0, 0 }; // {rxnIndex, rateIndex, isBackRxn=false}
    bsc_run(sys, rxnItr);

    // ---- 1. the state change itself ----
    std::cerr << "  -> Checking the interface state changed from 'U' to 'P'\n";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'P')
        << "state identity should be taken from forwardRxns[0].productListNew[1]";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 1)
        << "absolute interface index should become the product's absIfaceIndex";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, 1)
        << "relative state index should be 1 ('P' is the second state in the template)";

    // ---- 2. species copy numbers ----
    // NOTE: the source indexes copyNumSpecies with the *relative* state index,
    // i.e. 0 for 'U' and 1 for 'P' in this template.
    std::cerr << "  -> Checking copyNumSpecies[0] decremented and copyNumSpecies[1] incremented\n";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], startCopyU - 1) << "old state count should drop by one";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], startCopyP + 1) << "new state count should rise by one";

    // ---- 3. observable ----
    std::cerr << "  -> Checking the forward observable \"A_P\" was incremented\n";
    EXPECT_EQ(sys.observablesList["A_P"], startObs + 1) << "forward reactions increment their observable";
    EXPECT_EQ(sys.observablesList["A_U"], 7) << "the back-reaction observable must be untouched";

    // ---- 4. geometry ----
    // sigma = iface(B) - iface(A) = (10,0,0); displaceFrac = (10-1)/10 = 0.9.
    // transVec1 = -sigma * (D_A/DSum) * 0.9 = (-2.25,0,0)  -> applied to A's complex
    // transVec2 = +sigma * (D_B/DSum) * 0.9 = (+6.75,0,0)  -> applied to B's complex
    // The pair COM then moves from x=5 to x=7.25 and is forced back by -2.25,
    // giving A at x=-4.5 and B at x=14.5.
    std::cerr << "  -> Checking the displacement produced by the 'move to sigma' step\n";
    std::cerr << "     final A.x = " << sys.moleculeList[0].comCoord.x
              << ", final B.x = " << sys.moleculeList[1].comCoord.x << '\n';
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, -4.5, 1e-9) << "state-change molecule moved by -2.25 twice";
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, 14.5, 1e-9) << "facilitator moved by +6.75 then -2.25";
    EXPECT_NEAR(sys.moleculeList[0].comCoord.y, 0.0, 1e-12) << "no motion perpendicular to sigma";
    EXPECT_NEAR(sys.moleculeList[0].comCoord.z, 0.0, 1e-12) << "no motion perpendicular to sigma";

    std::cerr << "  -> Checking the pair centre of mass was restored\n";
    const double endPairComX = 0.5 * (sys.moleculeList[0].comCoord.x + sys.moleculeList[1].comCoord.x);
    EXPECT_NEAR(endPairComX, startPairComX, 1e-9)
        << "the routine translates the pair back onto its original centre of mass";

    // Interfaces of point molecules must follow their centre of mass exactly.
    EXPECT_NEAR(sys.moleculeList[0].interfaceList[0].coord.x, sys.moleculeList[0].comCoord.x, 1e-12)
        << "point molecule: interface must stay on the COM";
    EXPECT_NEAR(sys.moleculeList[1].interfaceList[0].coord.x, sys.moleculeList[1].comCoord.x, 1e-12)
        << "point molecule: interface must stay on the COM";

    // ---- 5. temporary association coordinates cleared ----
    std::cerr << "  -> Checking temporary association coordinates were cleared\n";
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty()) << "tmpICoords must be cleared on success";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty()) << "tmpICoords must be cleared on success";

    // ---- 6. trajectory status / cross lists ----
    std::cerr << "  -> Checking trajStatus, ncross and the crossing lists\n";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "the state-change molecule must be marked propagated";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "the facilitator must be marked propagated";
    EXPECT_EQ(sys.complexList[0].ncross, -1) << "state-change complex ncross is set to -1";
    EXPECT_EQ(sys.complexList[1].ncross, -1) << "facilitator complex ncross is set to -1";
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty()) << "state-change molecule crossbase is cleared";
    EXPECT_TRUE(sys.moleculeList[1].crossbase.empty()) << "facilitator crossbase is cleared";

    // The bystander listed both reactants, so both of its probabilities go to 0.
    ASSERT_EQ(sys.moleculeList[2].probvec.size(), 2u) << "bystander should still hold two crossing entries";
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[0], 0.0)
        << "probability of reacting with the state-change molecule must be zeroed";
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[1], 0.0)
        << "probability of reacting with the facilitator must be zeroed";

    // ---- 7. complex properties refreshed ----
    std::cerr << "  -> Checking Complex::update_properties() was applied\n";
    EXPECT_NEAR(sys.complexList[0].comCoord.x, sys.moleculeList[0].comCoord.x, 1e-9)
        << "single-member complex COM should track its molecule";
    EXPECT_NEAR(sys.complexList[1].comCoord.x, sys.moleculeList[1].comCoord.x, 1e-9)
        << "single-member complex COM should track its molecule";
}

// ---------------------------------------------------------------------------
// Test 2: back reaction (rxnItr[2] != 0) in solution.
// ---------------------------------------------------------------------------
void test_bsc_sphere_back_reaction()
{
    std::cerr << "\n[TEST] test_bsc_sphere_back_reaction\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_sphere.cpp\n"
              << "  Function:    perform_bimolecular_state_change_sphere (back-reaction path)\n"
              << "  Scenario:    A starts in state 'P' and rxnItr = {0,0,1}, so the routine\n"
              << "               reads backRxns[0].productListNew[1] and resolves the binding\n"
              << "               radius through backRxns[0].conjForwardRxnIndex.\n"
              << "  Criteria:    state flips P->U, copy numbers swap the other way and the\n"
              << "               back-reaction observable is DECREMENTED.\n";

    bsc_init_rng();
    BscSphereSystem sys = bsc_build_system(/*onMembrane=*/false, /*checkOverlap=*/false, /*startState=*/'P');

    const int startCopyU = sys.counterArrays.copyNumSpecies[0];
    const int startCopyP = sys.counterArrays.copyNumSpecies[1];
    const int startObsU = sys.observablesList["A_U"];
    const int startObsP = sys.observablesList["A_P"];

    std::array<int, 3> rxnItr { 0, 0, 1 }; // isStateChangeBackRxn == true
    bsc_run(sys, rxnItr);

    std::cerr << "  -> Checking the interface state changed back from 'P' to 'U'\n";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'U')
        << "back reaction product state should be 'U'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 0) << "absolute index should be that of state 'U'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, 0) << "relative state index should be 0";

    std::cerr << "  -> Checking the species counters moved from 'P' back to 'U'\n";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], startCopyP - 1) << "'P' count should drop by one";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], startCopyU + 1) << "'U' count should rise by one";

    std::cerr << "  -> Checking the observable of the back reaction was decremented\n";
    EXPECT_EQ(sys.observablesList["A_U"], startObsU - 1) << "back reactions decrement their observable";
    EXPECT_EQ(sys.observablesList["A_P"], startObsP) << "the forward observable must be untouched";

    // The geometry is identical to the forward case since the conjugate forward
    // reaction supplies the same binding radius.
    std::cerr << "  -> Checking that the geometry matches the forward case (same bindRadius)\n";
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, -4.5, 1e-9) << "same displacement as the forward reaction";
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, 14.5, 1e-9) << "same displacement as the forward reaction";
}

// ---------------------------------------------------------------------------
// Test 3: the event is cancelled because the structures overlap.
// ---------------------------------------------------------------------------
void test_bsc_sphere_cancelled_by_overlap()
{
    std::cerr << "\n[TEST] test_bsc_sphere_cancelled_by_overlap\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_sphere.cpp\n"
              << "  Function:    perform_bimolecular_state_change_sphere (cancel path)\n"
              << "  Scenario:    both templates set checkOverlap=true and\n"
              << "               params.overlapSepLimit=100 nm, which is far larger than any\n"
              << "               separation the pair can reach, so check_for_structure_overlap()\n"
              << "               must cancel the event.\n"
              << "  Criteria:    nothing is committed: state, coordinates, counters,\n"
              << "               observables, ncross and the crossing lists are unchanged,\n"
              << "               and the temporary coordinates are released.\n";

    bsc_init_rng();
    BscSphereSystem sys = bsc_build_system(/*onMembrane=*/false, /*checkOverlap=*/true, /*startState=*/'U');

    const Coord startA = sys.moleculeList[0].comCoord;
    const Coord startB = sys.moleculeList[1].comCoord;
    const int startCopyU = sys.counterArrays.copyNumSpecies[0];
    const int startCopyP = sys.counterArrays.copyNumSpecies[1];
    const int startObs = sys.observablesList["A_P"];

    std::array<int, 3> rxnItr { 0, 0, 0 };
    bsc_run(sys, rxnItr);

    std::cerr << "  -> Checking the state was NOT changed\n";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'U') << "state must remain 'U' after cancelling";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, 0) << "relative state index must remain 0";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 0) << "absolute index must remain 0";

    std::cerr << "  -> Checking the real coordinates were not modified\n";
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, startA.x, 1e-12) << "state-change molecule must not move";
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, startB.x, 1e-12) << "facilitator must not move";

    std::cerr << "  -> Checking the temporary association coordinates were released\n";
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty()) << "tmpICoords must be cleared when cancelling";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty()) << "tmpICoords must be cleared when cancelling";

    std::cerr << "  -> Checking counters, observables and cross lists were untouched\n";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], startCopyU) << "no copy-number bookkeeping on cancel";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], startCopyP) << "no copy-number bookkeeping on cancel";
    EXPECT_EQ(sys.observablesList["A_P"], startObs) << "no observable update on cancel";
    EXPECT_EQ(sys.complexList[0].ncross, 2) << "ncross must keep its incoming value on cancel";
    EXPECT_EQ(sys.complexList[1].ncross, 2) << "ncross must keep its incoming value on cancel";
    EXPECT_EQ(sys.moleculeList[0].crossbase.size(), 1u) << "crossbase must not be cleared on cancel";
    EXPECT_EQ(sys.moleculeList[1].crossbase.size(), 1u) << "crossbase must not be cleared on cancel";
    ASSERT_EQ(sys.moleculeList[2].probvec.size(), 2u);
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[0], 0.25) << "partner probabilities untouched on cancel";
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[1], 0.75) << "partner probabilities untouched on cancel";

    std::cerr << "  -> Checking trajStatus is still 'none'\n";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "no propagation is recorded when the event is cancelled";
}

// ---------------------------------------------------------------------------
// Test 4: both complexes on the sphere surface (the 2D / membrane branch).
// ---------------------------------------------------------------------------
void test_bsc_sphere_on_membrane_branch()
{
    std::cerr << "\n[TEST] test_bsc_sphere_on_membrane_branch\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_sphere.cpp\n"
              << "  Function:    perform_bimolecular_state_change_sphere (2D sphere branch)\n"
              << "  Scenario:    both complexes have D.z == 0 and sit on a sphere of radius\n"
              << "               100 nm, separated by a 0.1 rad geodesic (10 nm arc).  This\n"
              << "               drives the geodesic code path: calc_bindRadius2D(),\n"
              << "               get_geodesic_distance(), find_position_after_association()\n"
              << "               and Complex::update_association_coords_sphere().\n"
              << "  Criteria:    the reaction completes (state flips, counters/observable\n"
              << "               updated), all coordinates stay finite and both molecules\n"
              << "               remain on (or very near) the spherical surface.\n";

    bsc_init_rng();
    BscSphereSystem sys = bsc_build_system(/*onMembrane=*/true, /*checkOverlap=*/false, /*startState=*/'U');

    const int startCopyU = sys.counterArrays.copyNumSpecies[0];
    const int startCopyP = sys.counterArrays.copyNumSpecies[1];
    const int startObs = sys.observablesList["A_P"];

    std::cerr << "  Initial |A| = " << bsc_mag(sys.moleculeList[0].comCoord)
              << ", |B| = " << bsc_mag(sys.moleculeList[1].comCoord)
              << ", sphere radius = " << sys.membrane.sphereR << '\n';

    std::array<int, 3> rxnItr { 0, 0, 0 };
    bsc_run(sys, rxnItr);

    // ---- the state change and the bookkeeping still happen on the membrane ----
    std::cerr << "  -> Checking the state change was applied on the membrane path\n";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'P') << "state should flip to 'P'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, 1) << "relative state index should be 1";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], startCopyU - 1) << "'U' count should drop";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], startCopyP + 1) << "'P' count should rise";
    EXPECT_EQ(sys.observablesList["A_P"], startObs + 1) << "observable should be incremented";

    // ---- the geodesic maths must not produce NaN/inf ----
    std::cerr << "  -> Checking all coordinates are finite after the geodesic move\n";
    for (int m = 0; m < 2; ++m) {
        const Coord& c = sys.moleculeList[m].comCoord;
        EXPECT_TRUE(std::isfinite(c.x)) << "molecule " << m << " x coordinate must be finite";
        EXPECT_TRUE(std::isfinite(c.y)) << "molecule " << m << " y coordinate must be finite";
        EXPECT_TRUE(std::isfinite(c.z)) << "molecule " << m << " z coordinate must be finite";
    }

    // ---- the surface-bound molecules must remain on the sphere ----
    // translate_on_sphere()/rotate_on_sphere() preserve the radius; a 1 nm
    // tolerance also absorbs any small inward push from the boundary check.
    const double magA = bsc_mag(sys.moleculeList[0].comCoord);
    const double magB = bsc_mag(sys.moleculeList[1].comCoord);
    std::cerr << "  Final |A| = " << magA << ", |B| = " << magB << " (sphere radius "
              << sys.membrane.sphereR << ")\n";
    EXPECT_NEAR(magA, sys.membrane.sphereR, 1.0) << "state-change molecule should stay on the sphere";
    EXPECT_NEAR(magB, sys.membrane.sphereR, 1.0) << "facilitator should stay on the sphere";

    // ---- the usual post-conditions ----
    std::cerr << "  -> Checking temporary coordinates, trajStatus and cross lists\n";
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty()) << "tmpICoords must be cleared on success";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty()) << "tmpICoords must be cleared on success";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "the state-change molecule must be marked propagated";
    EXPECT_EQ(sys.complexList[0].ncross, -1) << "state-change complex ncross is set to -1";
    EXPECT_EQ(sys.complexList[1].ncross, -1) << "facilitator complex ncross is set to -1";
    ASSERT_EQ(sys.moleculeList[2].probvec.size(), 2u);
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[0], 0.0) << "crossing probabilities must be zeroed";
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[1], 0.0) << "crossing probabilities must be zeroed";
}

// ---------------------------------------------------------------------------
// Test 5: an unobserved reaction and a missing observable label must not
//         disturb the observables map.
// ---------------------------------------------------------------------------
void test_bsc_sphere_unobserved_reaction()
{
    std::cerr << "\n[TEST] test_bsc_sphere_unobserved_reaction\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_sphere.cpp\n"
              << "  Function:    perform_bimolecular_state_change_sphere (observable handling)\n"
              << "  Scenario:    (a) the forward reaction has isObserved = false, and\n"
              << "               (b) a second run uses a label that is not in the map.\n"
              << "  Criteria:    the state change still happens and the observables map is\n"
              << "               left with exactly the entries (and values) it started with.\n";

    bsc_init_rng();

    // (a) isObserved == false --------------------------------------------------
    {
        BscSphereSystem sys = bsc_build_system(false, false, 'U');
        sys.forwardRxns[0].isObserved = false;
        const std::size_t startSize = sys.observablesList.size();

        std::array<int, 3> rxnItr { 0, 0, 0 };
        bsc_run(sys, rxnItr);

        std::cerr << "  -> (a) Checking an unobserved reaction leaves the map alone\n";
        EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'P')
            << "the state change must still be performed";
        EXPECT_EQ(sys.observablesList["A_P"], 5) << "an unobserved reaction must not touch observables";
        EXPECT_EQ(sys.observablesList.size(), startSize) << "no new observable entries may appear";
    }

    // (b) observed, but the label is not registered ----------------------------
    {
        BscSphereSystem sys = bsc_build_system(false, false, 'U');
        sys.forwardRxns[0].observeLabel = "not_registered";
        const std::size_t startSize = sys.observablesList.size();

        std::array<int, 3> rxnItr { 0, 0, 0 };
        bsc_run(sys, rxnItr);

        std::cerr << "  -> (b) Checking an unknown observable label is handled gracefully\n";
        EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'P')
            << "the state change must still be performed";
        EXPECT_EQ(sys.observablesList.size(), startSize)
            << "an unknown label must not silently insert a new observable";
        EXPECT_EQ(sys.observablesList["A_P"], 5) << "existing observables must be unchanged";
    }
}

// ---------------------------------------------------------------------------
// GoogleTest wrappers - each named helper is run inside its own TEST so all of
// them execute even if one of them reports a failure.
// ---------------------------------------------------------------------------
TEST(PerformBimolecularStateChangeSphere, ForwardInSolution) { test_bsc_sphere_forward_in_solution(); }
TEST(PerformBimolecularStateChangeSphere, BackReaction) { test_bsc_sphere_back_reaction(); }
TEST(PerformBimolecularStateChangeSphere, CancelledByOverlap) { test_bsc_sphere_cancelled_by_overlap(); }
TEST(PerformBimolecularStateChangeSphere, OnMembraneBranch) { test_bsc_sphere_on_membrane_branch(); }
TEST(PerformBimolecularStateChangeSphere, UnobservedReaction) { test_bsc_sphere_unobserved_reaction(); }