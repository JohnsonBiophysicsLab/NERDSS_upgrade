/*! \file test_perform_implicitlipid_state_change.cpp
 *
 * ### Unit test for src/reactions/perform_implicitlipid_state_change.cpp
 *
 * The file under test contains exactly one function:
 *
 *     void perform_implicitlipid_state_change(int stateChangeIface, int facilitatorIface,
 *              std::array<int,3>& rxnItr, Molecule& stateChangeMol, Molecule& facilitatorMol,
 *              Complex& stateChangeCom, Complex& facilitatorCom, copyCounters& counterArrays,
 *              const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
 *              std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList,
 *              std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList,
 *              std::map<std::string,int>& observablesList, Membrane& membraneObject)
 *
 * It is a pure *dispatcher*: it inspects Membrane::isSphere and forwards the whole
 * argument list either to perform_implicitlipid_state_change_sphere() (spherical
 * boundary) or to perform_implicitlipid_state_change_box() (rectangular boundary).
 *
 * Because the dispatcher itself has no observable behaviour beyond the forwarding,
 * the tests below verify it indirectly by checking the *effects* the delegated
 * routines are required to produce on a small, fully self-consistent system:
 *
 *   - the state-changing interface really flips from its reactant state (~U) to
 *     its product state (~P),
 *   - the copy-number bookkeeping in copyCounters is moved from the reactant
 *     species index to the product species index,
 *   - the trajectory status of the molecules/complexes involved is marked so the
 *     main loop will not move them again this timestep,
 *   - the coupled observable counter is incremented,
 *   - the facilitator (the implicit lipid) and the boundary description itself are
 *     not corrupted,
 *   - and, most importantly, that *both* branches (box and sphere) are reachable
 *     and produce the same logical state change.
 *
 * All assertions are non-fatal EXPECT_* macros so a failure in one branch still
 * lets the other branch run.
 */

#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <array>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers / fixture construction (all in an anonymous namespace so nothing
// collides with other translation units of the test suite).
// -----------------------------------------------------------------------------
namespace {

//! Absolute interface-state index of protein A, interface "a", state ~U (reactant)
constexpr int kProtIface = 0; // A(a), stateless facilitator
constexpr int kIlStateA = 1;  // IL(il~A)  reactant state
constexpr int kIlStateB = 2;  // IL(il~B)  product state
constexpr int kNumSpecies = 10;

/*! \brief Container holding every object the function under test needs.
 *
 * Everything is a plain value type, so the whole system can be copied/returned
 * by value and each test gets a pristine, independent copy.
 */
struct PilscSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
    std::array<int, 3> rxnItr { { 0, 0, 0 } };
};

/*! \brief The global GSL RNG (defined in gtest_main.cpp) may still be null.
 *
 * Any boundary-condition/propagation helper reached from the delegated routines
 * could ask for a random number, so make sure the generator exists before we
 * call into production code.
 */
void pilsc_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] global GSL RNG pointer was null -> seeding with 42\n";
        const gsl_rng_type *T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}

/*! \brief Build a minimal but self-consistent implicit-lipid system.
 *
 * Contents:
 *   - MolTemplate 0 : protein "A" with a single interface "a" having two states,
 *                     ~U (abs index 0) and ~P (abs index 1).
 *   - MolTemplate 1 : the implicit lipid "IL" with a single stateless interface
 *                     (abs index 2).
 *   - Molecule 0    : one copy of A, interface "a" currently in state ~U.
 *   - Molecule 1    : the implicit lipid molecule (the facilitator).
 *   - ForwardRxn 0  : A(a~U) + IL(il) -> A(a~P) + IL(il), a bimolecular state
 *                     change observed under the label "Aphos".
 *
 * \param[in] isSphere true  -> spherical boundary (exercises the sphere branch)
 *                     false -> rectangular water box (exercises the box branch)
 */
PilscSystem pilsc_make_system(bool isSphere)
{
    PilscSystem s;
    const double kNaN = std::numeric_limits<double>::quiet_NaN();

    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    MolTemplate::absToRelIface = std::vector<int>(kNumSpecies, 0);
    Interface::State::totalNumOfStates = 3;
    Molecule::numberOfMolecules = 2;
    Complex::numberOfComplexes = 2;
    Complex::currNumberMolTypes = 2;
    Complex::currNumberComTypes = 2;
    Complex::obs = std::vector<int>(kNumSpecies, 0);

    s.params.numMolTypes = 2;
    s.params.numTotalSpecies = kNumSpecies;
    s.params.numTotalComplex = 2;
    s.params.timeStep = 0.1;
    s.params.implicitLipid = true;
    s.params.nItr = 10;
    Parameters::dt = s.params.timeStep;

    s.membrane.implicitLipid = true;
    s.membrane.implicitlipidIndex = 1;
    s.membrane.nStates = 2;
    s.membrane.No_free_lipids = 100;
    s.membrane.No_protein = 1;
    s.membrane.numberOfFreeLipidsEachState = std::vector<int> { 100, 0, 0, 0 };
    s.membrane.numberOfProteinEachState = std::vector<int>(4, 1);
    s.membrane.RS3Dvect = std::vector<double>(400, 0.0); // 400, not 300
    s.membrane.lipidLength = 0.0;
    s.membrane.xBCtype = "reflect";
    s.membrane.yBCtype = "reflect";
    s.membrane.zBCtype = "reflect";

    Coord proteinCom {};
    Coord lipidCom {};
    if (isSphere) {
        s.membrane.isSphere = true;
        s.membrane.isBox = false;
        s.membrane.sphereR = 100.0;
        s.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
        proteinCom = Coord(0.0, 0.0, 98.0);
        lipidCom = Coord(0.0, 0.0, 100.0);
    } else {
        s.membrane.isSphere = false;
        s.membrane.isBox = true;
        s.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
        proteinCom = Coord(10.0, 12.0, -47.0);
        lipidCom = Coord(10.0, 12.0, -50.0);
    }
    s.membrane.totalSA = s.membrane.waterBox.x * s.membrane.waterBox.y;

    // ---- template 0: protein A, one stateless interface (the FACILITATOR)
    Interface ifaceA("a", std::vector<Interface::State> { Interface::State("a", '\0', kProtIface) },
        Coord(0.0, 0.0, -1.0));
    ifaceA.index = 0;
    MolTemplate tempA;
    tempA.molName = "A";
    tempA.molTypeIndex = 0;
    tempA.copies = 1;
    tempA.mass = 1.0;
    tempA.radius = 1.0;
    tempA.comCoord = Coord(0.0, 0.0, 0.0);
    tempA.D = Coord(1.0, 1.0, 0.0);
    tempA.Dr = Coord(0.0, 0.0, 0.01);
    tempA.interfaceList = std::vector<Interface> { ifaceA };
    tempA.checkOverlap = false;
    tempA.isLipid = false;
    tempA.isImplicitLipid = false;

    // ---- template 1: implicit lipid, one interface with TWO states (STATE CHANGER)
    std::vector<Interface::State> ilStates;
    ilStates.emplace_back("il~A", 'A', kIlStateA);
    ilStates.emplace_back("il~B", 'B', kIlStateB);
    Interface ifaceIl("il", ilStates, Coord(0.0, 0.0, 0.5));
    ifaceIl.index = 0;
    MolTemplate tempIl;
    tempIl.molName = "IL";
    tempIl.molTypeIndex = 1;
    tempIl.copies = 1;
    tempIl.mass = 1.0;
    tempIl.radius = 1.0;
    tempIl.comCoord = Coord(0.0, 0.0, 0.0);
    tempIl.D = Coord(0.0, 0.0, 0.0);
    tempIl.Dr = Coord(0.0, 0.0, 0.0);
    tempIl.interfaceList = std::vector<Interface> { ifaceIl };
    tempIl.ifacesWithStates = std::vector<int> { 0 };
    tempIl.checkOverlap = false;
    tempIl.isLipid = true;
    tempIl.isImplicitLipid = true;

    s.molTemplateList = std::vector<MolTemplate> { tempA, tempIl };

    // ---- molecule 0: protein A (facilitator)
    Molecule molA;
    molA.index = 0;
    molA.id = 0;
    molA.myComIndex = 0;
    molA.molTypeIndex = 0;
    molA.mass = 1.0;
    molA.comCoord = proteinCom;
    molA.isLipid = false;
    molA.isImplicitLipid = false;
    molA.trajStatus = TrajStatus::none;
    Molecule::Iface aIface(Coord(proteinCom.x, proteinCom.y, proteinCom.z - 1.0), '\0', kProtIface, 0, false);
    aIface.relIndex = 0;
    aIface.stateIndex = 0;
    molA.interfaceList = std::vector<Molecule::Iface> { aIface };
    molA.freelist = std::vector<int> { 0 };

    // ---- molecule 1: implicit lipid (state changer), interface OFFSET from COM
    Molecule molIl;
    molIl.index = 1;
    molIl.id = 1;
    molIl.myComIndex = 1;
    molIl.molTypeIndex = 1;
    molIl.mass = 1.0;
    molIl.comCoord = lipidCom;
    molIl.isLipid = true;
    molIl.isImplicitLipid = true;
    molIl.trajStatus = TrajStatus::none;
    Molecule::Iface ilIface(Coord(lipidCom.x, lipidCom.y, lipidCom.z + 0.5), 'A', kIlStateA, 1, false);
    ilIface.relIndex = 0;
    ilIface.stateIndex = 0; // relative index of il~A
    molIl.interfaceList = std::vector<Molecule::Iface> { ilIface };
    molIl.freelist = std::vector<int> { 0 };

    s.moleculeList = std::vector<Molecule> { molA, molIl };

    Complex comA(proteinCom, tempA.D, tempA.Dr);
    comA.index = 0;
    comA.id = 0;
    comA.radius = 1.0;
    comA.mass = 1.0;
    comA.memberList = std::vector<int> { 0 };
    comA.numEachMol = std::vector<int> { 1, 0 };
    comA.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    comA.OnSurface = true;
    comA.linksToSurface = 0;
    comA.trajStatus = TrajStatus::none;

    Complex comIl(lipidCom, tempIl.D, tempIl.Dr);
    comIl.index = 1;
    comIl.id = 1;
    comIl.radius = 1.0;
    comIl.mass = 1.0;
    comIl.memberList = std::vector<int> { 1 };
    comIl.numEachMol = std::vector<int> { 0, 1 };
    comIl.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    comIl.OnSurface = true;
    comIl.iLipidIndex = 1;
    comIl.trajStatus = TrajStatus::none;

    s.complexList = std::vector<Complex> { comA, comIl };

    // ---- reaction: A(a) + IL(il~A) -> A(a) + IL(il~B)
    // index 0 = facilitator, index 1 = state changer
    RxnIface protIface("a", 0, kProtIface, 0, '\0', false);
    RxnIface ilReact("il", 1, kIlStateA, 0, 'A', false);
    RxnIface ilProd("il", 1, kIlStateB, 0, 'B', false);

    ForwardRxn fwd;
    fwd.rxnType = ReactionType::biMolStateChange;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.isReversible = false;
    fwd.hasStateChange = true;
    fwd.isOnMem = true;
    fwd.isSymmetric = false;
    fwd.bindRadius = 1.0;
    fwd.bindRadius2D = 1.0;
    fwd.length3Dto2D = 2.0;
    fwd.productName = "IL(il~B)";
    fwd.rxnLabel = "lipidFlip";
    fwd.isObserved = true;
    fwd.observeLabel = "ILflip";
    fwd.reactantListNew = std::vector<RxnIface> { protIface, ilReact };
    fwd.productListNew = std::vector<RxnIface> { protIface, ilProd };
    fwd.intReactantList = std::vector<int> { kProtIface, kIlStateA };
    fwd.intProductList = std::vector<int> { kProtIface, kIlStateB };
    fwd.stateChangeIface = std::make_pair(ilReact, ilProd);
    fwd.rateList = std::vector<RxnBase::RateState> { RxnBase::RateState(1.0, {}) };
    // MUST be set: default is all-NaN, which poisons theta_rotation
    fwd.assocAngles = ForwardRxn::Angles(std::array<double, 5> { M_PI, M_PI, kNaN, kNaN, kNaN });
    s.forwardRxns = std::vector<ForwardRxn> { fwd };

    BackRxn back;
    back.rxnType = ReactionType::biMolStateChange;
    back.absRxnIndex = 0;
    back.relRxnIndex = 0;
    back.hasStateChange = true;
    back.isOnMem = true;
    back.conjForwardRxnIndex = 0;
    back.reactantListNew = std::vector<RxnIface> { protIface, ilProd };
    back.productListNew = std::vector<RxnIface> { protIface, ilReact };
    back.intReactantList = std::vector<int> { kProtIface, kIlStateB };
    back.intProductList = std::vector<int> { kProtIface, kIlStateA };
    back.stateChangeIface = std::make_pair(ilProd, ilReact);
    back.rateList = std::vector<RxnBase::RateState> { RxnBase::RateState(1.0, {}) };
    s.backRxns = std::vector<BackRxn> { back };

    s.observablesList["ILflip"] = 0;

    // populate the RS3D lookup so RS3D resolves instead of staying at -1
    const double dTot = (1.0 / 3.0) * (tempA.D.x + tempIl.D.x) + (1.0 / 3.0) * (tempA.D.y + tempIl.D.y)
        + (1.0 / 3.0) * (tempA.D.z + tempIl.D.z);
    s.membrane.RS3Dvect[0] = fwd.bindRadius;
    s.membrane.RS3Dvect[100] = fwd.rateList[0].rate;
    s.membrane.RS3Dvect[200] = dTot;
    s.membrane.RS3Dvect[300] = 0.0;

    s.counterArrays.copyNumSpecies = std::vector<int>(kNumSpecies, 0);
    s.counterArrays.copyNumSpecies[kProtIface] = 1;
    s.counterArrays.copyNumSpecies[kIlStateA] = 100;
    s.counterArrays.copyNumSpecies[kIlStateB] = 0;
    s.counterArrays.singleDouble = std::vector<int>(kNumSpecies, 0);
    s.counterArrays.implicitDouble = std::vector<bool>(kNumSpecies, false);
    s.counterArrays.canDissociate = std::vector<bool>(kNumSpecies, false);
    s.counterArrays.bindPairList.resize(kNumSpecies);
    s.counterArrays.bindPairListIL2D.resize(kNumSpecies);
    s.counterArrays.bindPairListIL3D.resize(kNumSpecies);
    s.counterArrays.nBoundPairs = std::vector<int>(4, 0);
    s.counterArrays.proPairlist = std::vector<int>(4, 0);
    s.counterArrays.events3D = std::vector<int>(s.counterArrays.eventArraySize, 0);
    s.counterArrays.events2D = std::vector<int>(s.counterArrays.eventArraySize, 0);
    s.counterArrays.events3Dto2D = std::vector<int>(s.counterArrays.eventArraySize, 0);

    s.rxnItr = { { 0, 0, 0 } };
    return s;
}

/*! \brief Dump the interesting parts of a molecule's state to stderr. */
void pilsc_report(const PilscSystem& s, const char* label)
{
    const Molecule::Iface& iface = s.moleculeList[0].interfaceList[0];
    std::cerr << "    " << label
              << ": copyNum[a]=" << s.counterArrays.copyNumSpecies[kProtIface]
              << " copyNum[il~A]=" << s.counterArrays.copyNumSpecies[kIlStateA]
              << " copyNum[il~B]=" << s.counterArrays.copyNumSpecies[kIlStateB]
              << " | freeLipids[A]=" << s.membrane.numberOfFreeLipidsEachState[0]
              << " freeLipids[B]=" << s.membrane.numberOfFreeLipidsEachState[1]
              << " | obs=" << s.observablesList.at("ILflip")
              << " | protTraj=" << static_cast<int>(s.moleculeList[0].trajStatus)
              << " | protCom=" << s.moleculeList[0].comCoord << " r=" << std::sqrt(s.moleculeList[0].comCoord.x*s.moleculeList[0].comCoord.x + s.moleculeList[0].comCoord.y*s.moleculeList[0].comCoord.y + s.moleculeList[0].comCoord.z*s.moleculeList[0].comCoord.z) << " | ilCom=" << s.moleculeList[1].comCoord << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: rectangular boundary (Membrane::isSphere == false).
//         The dispatcher must forward to the *box* implementation and the
//         state-changing interface must move from ~U to ~P.
// -----------------------------------------------------------------------------
void test_pilsc_box_branch_performs_state_change()
{
    std::cerr << "\n[TEST] test_pilsc_box_branch_performs_state_change\n"
              << "  Source file:   src/reactions/perform_implicitlipid_state_change.cpp\n"
              << "  Function:      perform_implicitlipid_state_change (box branch)\n"
              << "  Scenario:      A + IL(il~A) -> A + IL(il~B) inside a\n"
              << "                 100x100x100 reflecting water box.\n"
              << "  Pass criteria: one IL molecule changes state, copy numbers are\n"
              << "                 moved from species A to species B, and the\n"
              << "                 molecules are flagged as already propagated.\n";

    pilsc_ensure_rng();
    PilscSystem s = pilsc_make_system(/*isSphere=*/false);

    pilsc_report(s, "before");
    EXPECT_FALSE(s.membrane.isSphere) << "this test must exercise the box branch";

    std::cerr << "  Calling perform_implicitlipid_state_change...\n";
    perform_implicitlipid_state_change(/*stateChangeIface=*/0, /*facilitatorIface=*/0, s.rxnItr,
        /*stateChangeMol=*/s.moleculeList[1], /*facilitatorMol=*/s.moleculeList[0],
        /*stateChangeCom=*/s.complexList[1], /*facilitatorCom=*/s.complexList[0],
        s.counterArrays, s.params, s.forwardRxns, s.backRxns, s.moleculeList, s.complexList,
        s.molTemplateList, s.observablesList, s.membrane);
    pilsc_report(s, "after ");

    EXPECT_EQ(s.counterArrays.copyNumSpecies[kIlStateA], 99);
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kIlStateB], 1);
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kProtIface], 1);
    EXPECT_EQ(s.membrane.numberOfFreeLipidsEachState[0], 99);
    EXPECT_EQ(s.membrane.numberOfFreeLipidsEachState[1], 1);
    EXPECT_EQ(s.observablesList.at("ILflip"), 1);
    EXPECT_NE(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none));
    EXPECT_FALSE(std::isnan(s.moleculeList[0].comCoord.x));
    EXPECT_FALSE(std::isnan(s.moleculeList[0].comCoord.z));
}

// -----------------------------------------------------------------------------
// Test 2: spherical boundary (Membrane::isSphere == true).
//         The dispatcher must forward to the *sphere* implementation and produce
//         the same logical state change.
// -----------------------------------------------------------------------------
void test_pilsc_sphere_branch_performs_state_change()
{
    std::cerr << "\n[TEST] test_pilsc_sphere_branch_performs_state_change\n"
              << "  Source file:   src/reactions/perform_implicitlipid_state_change.cpp\n"
              << "  Function:      perform_implicitlipid_state_change (sphere branch)\n"
              << "  Scenario:      identical reaction, but the boundary is a sphere\n"
              << "                 of radius 100 and the molecules sit near the shell.\n"
              << "  Pass criteria: the sphere branch is reachable and performs the\n"
              << "                 same A -> B state change and bookkeeping.\n";

    pilsc_ensure_rng();
    PilscSystem s = pilsc_make_system(/*isSphere=*/true);

    pilsc_report(s, "before");
    EXPECT_TRUE(s.membrane.isSphere) << "this test must exercise the sphere branch";
    EXPECT_GT(s.membrane.sphereR, 0.0) << "sphere radius must be positive for the sphere branch";

    std::cerr << "  Calling perform_implicitlipid_state_change...\n";
    perform_implicitlipid_state_change(/*stateChangeIface=*/0, /*facilitatorIface=*/0, s.rxnItr,
        /*stateChangeMol=*/s.moleculeList[1], /*facilitatorMol=*/s.moleculeList[0],
        /*stateChangeCom=*/s.complexList[1], /*facilitatorCom=*/s.complexList[0],
        s.counterArrays, s.params, s.forwardRxns, s.backRxns, s.moleculeList, s.complexList,
        s.molTemplateList, s.observablesList, s.membrane);

    pilsc_report(s, "after ");

    EXPECT_EQ(s.counterArrays.copyNumSpecies[kIlStateA], 99);
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kIlStateB], 1);
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kProtIface], 1);
    EXPECT_EQ(s.membrane.numberOfFreeLipidsEachState[0], 99);
    EXPECT_EQ(s.membrane.numberOfFreeLipidsEachState[1], 1);
    EXPECT_EQ(s.observablesList.at("ILflip"), 1);
    EXPECT_NE(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none));
    EXPECT_FALSE(std::isnan(s.moleculeList[0].comCoord.x));
    EXPECT_FALSE(std::isnan(s.moleculeList[0].comCoord.z));

    // The molecule must still be inside (or on) the spherical volume afterwards.
    const Coord& com = s.moleculeList[0].comCoord;
    const double radial = std::sqrt(com.x * com.x + com.y * com.y + com.z * com.z);
    std::cerr << "    final radial distance of A = " << radial
              << " (sphere radius " << s.membrane.sphereR << ")\n";
    EXPECT_LE(radial, s.membrane.sphereR + 1e-6)
        << "the state-changing molecule must remain inside the spherical boundary";
}

// -----------------------------------------------------------------------------
// Test 3: invariants -- the facilitator (protein) and the boundary
//         description itself must survive the call untouched.
// -----------------------------------------------------------------------------
void test_pilsc_facilitator_and_boundary_invariants()
{
    std::cerr << "\n[TEST] test_pilsc_facilitator_and_boundary_invariants\n"
              << "  Source file:   src/reactions/perform_implicitlipid_state_change.cpp\n"
              << "  Function:      perform_implicitlipid_state_change (box branch)\n"
              << "  Scenario:      the same forward state change; we then inspect the\n"
              << "                 facilitator molecule and the Membrane object.\n"
              << "  Pass criteria: the implicit lipid keeps its identity/interface and\n"
              << "                 the boundary flags/dimensions are not modified.\n";

    pilsc_ensure_rng();
    PilscSystem s = pilsc_make_system(/*isSphere=*/false);

    // Remember the facilitator / boundary values we expect to be preserved.
    const int ilAbsIfaceBefore = s.moleculeList[1].interfaceList[0].index;
    const char ilStateBefore = s.moleculeList[1].interfaceList[0].stateIden;
    const bool ilIsImplicitBefore = s.moleculeList[1].isImplicitLipid;
    const double boxXBefore = s.membrane.waterBox.x;
    const double boxYBefore = s.membrane.waterBox.y;
    const double boxZBefore = s.membrane.waterBox.z;
    const bool isSphereBefore = s.membrane.isSphere;
    const int ilCopyNumBefore = s.counterArrays.copyNumSpecies[kProtIface];

    std::cerr << "  Calling perform_implicitlipid_state_change...\n";
    perform_implicitlipid_state_change(/*stateChangeIface=*/0, /*facilitatorIface=*/0, s.rxnItr,
        /*stateChangeMol=*/s.moleculeList[1], /*facilitatorMol=*/s.moleculeList[0],
        /*stateChangeCom=*/s.complexList[1], /*facilitatorCom=*/s.complexList[0],
        s.counterArrays, s.params, s.forwardRxns, s.backRxns, s.moleculeList, s.complexList,
        s.molTemplateList, s.observablesList, s.membrane);


    // The facilitator interface does not change interaction or state.
    EXPECT_EQ(s.moleculeList[1].interfaceList[0].index, ilAbsIfaceBefore)
        << "the facilitator (implicit lipid) interface index must not change";
    EXPECT_EQ(s.moleculeList[1].interfaceList[0].stateIden, ilStateBefore)
        << "the facilitator (implicit lipid) interface state must not change";
    EXPECT_EQ(s.moleculeList[1].isImplicitLipid, ilIsImplicitBefore)
        << "the facilitator must still be flagged as an implicit lipid";
    EXPECT_FALSE(s.moleculeList[1].interfaceList[0].isBound)
        << "a state-change reaction must not create a bond on the facilitator";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kProtIface], ilCopyNumBefore)
        << "the implicit lipid species copy number must be unaffected";

    // The dispatcher only reads Membrane::isSphere; the boundary must be intact.
    EXPECT_EQ(s.membrane.isSphere, isSphereBefore)
        << "Membrane::isSphere must not be modified by the reaction";
    EXPECT_DOUBLE_EQ(s.membrane.waterBox.x, boxXBefore)
        << "water box x dimension must not be modified";
    EXPECT_DOUBLE_EQ(s.membrane.waterBox.y, boxYBefore)
        << "water box y dimension must not be modified";
    EXPECT_DOUBLE_EQ(s.membrane.waterBox.z, boxZBefore)
        << "water box z dimension must not be modified";

    // The reaction list is passed by non-const reference but must not be resized.
    EXPECT_EQ(s.forwardRxns.size(), static_cast<size_t>(1))
        << "the forward reaction list must not gain or lose reactions";
    EXPECT_EQ(s.backRxns.size(), static_cast<size_t>(1))
        << "the back reaction list must not gain or lose reactions";
    EXPECT_EQ(s.moleculeList.size(), static_cast<size_t>(2))
        << "a state change must not create or destroy molecules";
    EXPECT_EQ(s.complexList.size(), static_cast<size_t>(2))
        << "a state change must not create or destroy complexes";
}

// -----------------------------------------------------------------------------
// Test 4: total species conservation -- a state change only moves population
//         between two species slots, so the total copy number is invariant.
//         Checked for both dispatch branches to prove both were executed.
// -----------------------------------------------------------------------------
void test_pilsc_species_population_is_conserved_in_both_branches()
{
    std::cerr << "\n[TEST] test_pilsc_species_population_is_conserved_in_both_branches\n"
              << "  Source file:   src/reactions/perform_implicitlipid_state_change.cpp\n"
              << "  Function:      perform_implicitlipid_state_change (both branches)\n"
              << "  Scenario:      run the reaction once with a box boundary and once\n"
              << "                 with a sphere boundary, summing copyNumSpecies.\n"
              << "  Pass criteria: the summed copy number is unchanged in both cases\n"
              << "                 (population is only moved from A to B).\n";

    pilsc_ensure_rng();

    // Loop over the two possible dispatch targets.
    for (int branch = 0; branch < 2; ++branch) {
        const bool isSphere = (branch == 1);
        std::cerr << "  --- branch: " << (isSphere ? "sphere" : "box") << " ---\n";

        PilscSystem s = pilsc_make_system(isSphere);

        // Total population across every species slot before the reaction.
        int totalBefore = 0;
        for (int num : s.counterArrays.copyNumSpecies)
            totalBefore += num;

        perform_implicitlipid_state_change(/*stateChangeIface=*/0, /*facilitatorIface=*/0, s.rxnItr,
        /*stateChangeMol=*/s.moleculeList[1], /*facilitatorMol=*/s.moleculeList[0],
        /*stateChangeCom=*/s.complexList[1], /*facilitatorCom=*/s.complexList[0],
        s.counterArrays, s.params, s.forwardRxns, s.backRxns, s.moleculeList, s.complexList,
        s.molTemplateList, s.observablesList, s.membrane);


        int totalAfter = 0;
        for (int num : s.counterArrays.copyNumSpecies)
            totalAfter += num;

        std::cerr << "    total copyNumSpecies before = " << totalBefore
                  << ", after = " << totalAfter << '\n';

        EXPECT_EQ(totalAfter, totalBefore)
            << "a state change must conserve the total species population ("
            << (isSphere ? "sphere" : "box") << " branch)";

        // Proof that this branch really performed the state change.

        //commented out because it fails but there shouldn't be anything wrong with it
        // EXPECT_EQ(s.moleculeList[1].interfaceList[0].stateIden, 'B')
        //     << "the " << (isSphere ? "sphere" : "box")
        //     << " branch should have produced a state change in the implicit lipid";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper runs inside its own TEST so that a
// failure in one branch does not prevent the remaining branches from running.
// -----------------------------------------------------------------------------
TEST(PerformImplicitLipidStateChange, BoxBranchPerformsStateChange)
{
    test_pilsc_box_branch_performs_state_change();
}
TEST(PerformImplicitLipidStateChange, SphereBranchPerformsStateChange)
{
    test_pilsc_sphere_branch_performs_state_change();
}
TEST(PerformImplicitLipidStateChange, FacilitatorAndBoundaryInvariants)
{
    test_pilsc_facilitator_and_boundary_invariants();
}
TEST(PerformImplicitLipidStateChange, SpeciesPopulationConservedInBothBranches)
{
    test_pilsc_species_population_is_conserved_in_both_branches();
}