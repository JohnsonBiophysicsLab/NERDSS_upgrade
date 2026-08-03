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
constexpr int kPilscStateU = 0;
//! Absolute interface-state index of protein A, interface "a", state ~P (product)
constexpr int kPilscStateP = 1;
//! Absolute interface index of the implicit lipid's binding interface (facilitator)
constexpr int kPilscIlIface = 2;
//! Number of species slots reserved in the copy-number arrays (generously large)
constexpr int kPilscNumSpecies = 10;

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
    //! {reaction index, rate index, isStateChangeBackRxn}
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
        std::cerr << "  [setup] global GSL RNG pointer was null -> seeding with srand_gsl(1)\n";
        srand_gsl(1);
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

    /* ---------------------------------------------------------------- statics */
    // Several NERDSS classes carry static bookkeeping that the production code
    // reads; give it consistent values for this two-molecule system.
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    // absolute interface index -> relative interface index (all molecules here
    // only have one interface, so every absolute index maps to relative 0).
    MolTemplate::absToRelIface = std::vector<int>(kPilscNumSpecies, 0);
    Interface::State::totalNumOfStates = 3;
    Molecule::numberOfMolecules = 2;
    Complex::numberOfComplexes = 2;
    Complex::currNumberMolTypes = 2;
    Complex::currNumberComTypes = 2;
    Complex::obs = std::vector<int>(kPilscNumSpecies, 0);

    /* ------------------------------------------------------------- parameters */
    s.params.numMolTypes = 2;
    s.params.numTotalSpecies = kPilscNumSpecies;
    s.params.numTotalComplex = 2;
    s.params.timeStep = 0.1;
    s.params.implicitLipid = true;
    s.params.nItr = 10;
    Parameters::dt = s.params.timeStep;

    /* ---------------------------------------------------------------- membrane */
    s.membrane.implicitLipid = true;
    s.membrane.implicitlipidIndex = 1; // molecule index of the implicit lipid
    s.membrane.nStates = 1;
    s.membrane.No_free_lipids = 100;
    s.membrane.No_protein = 1;
    s.membrane.numberOfFreeLipidsEachState = std::vector<int>(4, 100);
    s.membrane.numberOfProteinEachState = std::vector<int>(4, 1);
    // The RS3D lookup table is indexed in some implicit-lipid helpers; provide a
    // safely sized table of zeros so no out-of-range access can occur.
    s.membrane.RS3Dvect = std::vector<double>(300, 0.0);
    s.membrane.lipidLength = 0.0;
    s.membrane.xBCtype = "reflect";
    s.membrane.yBCtype = "reflect";
    s.membrane.zBCtype = "reflect";

    Coord proteinCom {};
    Coord lipidCom {};

    if (isSphere) {
        // Spherical boundary of radius 100; put both molecules at a generic
        // (non-degenerate) point so no spherical-coordinate singularity is hit.
        s.membrane.isSphere = true;
        s.membrane.isBox = false;
        s.membrane.sphereR = 100.0;
        s.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
        proteinCom = Coord(30.0, 40.0, -80.0);   // |r| ~= 94.3 -> inside the sphere
        lipidCom = Coord(31.8, 42.4, -84.8);     // essentially on the shell
    } else {
        // Rectangular 100 x 100 x 100 water box, lipid on the bottom face.
        s.membrane.isSphere = false;
        s.membrane.isBox = true;
        s.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
        proteinCom = Coord(10.0, 12.0, -49.0);
        lipidCom = Coord(10.0, 12.0, -50.0);
    }
    s.membrane.totalSA = s.membrane.waterBox.x * s.membrane.waterBox.y;

    /* ----------------------------------------------------- molecule templates */
    // Protein A: one interface with two states (~U reactant, ~P product).
    std::vector<Interface::State> aStates;
    aStates.emplace_back("a~U", 'U', kPilscStateU);
    aStates.emplace_back("a~P", 'P', kPilscStateP);
    Interface ifaceA("a", aStates, Coord(0.0, 0.0, -1.0));
    ifaceA.index = 0;

    MolTemplate tempA;
    tempA.molName = "A";
    tempA.molTypeIndex = 0;
    tempA.copies = 1;
    tempA.mass = 1.0;
    tempA.radius = 1.0;
    tempA.comCoord = Coord(0.0, 0.0, 0.0);
    tempA.D = Coord(1.0, 1.0, 0.0);   // membrane bound -> no z translation
    tempA.Dr = Coord(0.0, 0.0, 0.01); // membrane bound -> rotation about z only
    tempA.interfaceList = std::vector<Interface> { ifaceA };
    tempA.ifacesWithStates = std::vector<int> { 0 };
    tempA.checkOverlap = false;
    tempA.isLipid = false;
    tempA.isImplicitLipid = false;
    tempA.isPoint = false;
    tempA.isRod = false;

    // The implicit lipid template: a stateless single interface, no diffusion.
    std::vector<Interface::State> ilStates;
    ilStates.emplace_back("il", '\0', kPilscIlIface);
    Interface ifaceIl("il", ilStates, Coord(0.0, 0.0, 0.0));
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
    tempIl.isLipid = true;
    tempIl.isImplicitLipid = true;

    s.molTemplateList = std::vector<MolTemplate> { tempA, tempIl };

    /* ------------------------------------------------------------- molecules */
    // Molecule 0: protein A, interface "a" in the reactant state ~U.
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
    Molecule::Iface molAIface(Coord(proteinCom.x, proteinCom.y, proteinCom.z - 1.0), 'U', kPilscStateU, 0, false);
    molAIface.relIndex = 0;
    molAIface.stateIndex = 0; // relative index of ~U inside stateList
    molA.interfaceList = std::vector<Molecule::Iface> { molAIface };
    molA.freelist = std::vector<int> { 0 };

    // Molecule 1: the implicit lipid (facilitator of the state change).
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
    Molecule::Iface molIlIface(lipidCom, '\0', kPilscIlIface, 1, false);
    molIlIface.relIndex = 0;
    molIlIface.stateIndex = 0;
    molIl.interfaceList = std::vector<Molecule::Iface> { molIlIface };
    molIl.freelist = std::vector<int> { 0 };

    s.moleculeList = std::vector<Molecule> { molA, molIl };

    /* -------------------------------------------------------------- complexes */
    Complex comA(proteinCom, tempA.D, tempA.Dr);
    comA.index = 0;
    comA.id = 0;
    comA.radius = 1.0;
    comA.mass = 1.0;
    comA.memberList = std::vector<int> { 0 };
    comA.numEachMol = std::vector<int> { 1, 0 };
    comA.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    comA.OnSurface = true;      // implicit-lipid reactions happen at the surface
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

    /* --------------------------------------------------------------- reaction */
    // A(a~U) + IL(il) <-> A(a~P) + IL(il) : bimolecular state change.
    RxnIface reactA("a", 0, kPilscStateU, 0, 'U', false);
    RxnIface prodA("a", 0, kPilscStateP, 0, 'P', false);
    RxnIface ilFacilitator("il", 1, kPilscIlIface, 0, '\0', false);

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
    fwd.productName = "A(a~P)";
    fwd.rxnLabel = "phosphorylation";
    fwd.isObserved = true;
    fwd.observeLabel = "Aphos";
    fwd.reactantListNew = std::vector<RxnIface> { reactA, ilFacilitator };
    fwd.productListNew = std::vector<RxnIface> { prodA, ilFacilitator };
    fwd.intReactantList = std::vector<int> { kPilscStateU, kPilscIlIface };
    fwd.intProductList = std::vector<int> { kPilscStateP, kPilscIlIface };
    fwd.stateChangeIface = std::make_pair(reactA, prodA);
    fwd.rateList = std::vector<RxnBase::RateState> { RxnBase::RateState(1.0, {}) };
    s.forwardRxns = std::vector<ForwardRxn> { fwd };

    // A matching back reaction so any lookup into backRxns is well defined.
    BackRxn back;
    back.rxnType = ReactionType::biMolStateChange;
    back.absRxnIndex = 0;
    back.relRxnIndex = 0;
    back.hasStateChange = true;
    back.isOnMem = true;
    back.conjForwardRxnIndex = 0;
    back.reactantListNew = std::vector<RxnIface> { prodA, ilFacilitator };
    back.productListNew = std::vector<RxnIface> { reactA, ilFacilitator };
    back.intReactantList = std::vector<int> { kPilscStateP, kPilscIlIface };
    back.intProductList = std::vector<int> { kPilscStateU, kPilscIlIface };
    back.stateChangeIface = std::make_pair(prodA, reactA);
    back.rateList = std::vector<RxnBase::RateState> { RxnBase::RateState(1.0, {}) };
    s.backRxns = std::vector<BackRxn> { back };

    /* ------------------------------------------------------------ observables */
    s.observablesList["Aphos"] = 0;

    /* --------------------------------------------------------- copy counters */
    s.counterArrays.copyNumSpecies = std::vector<int>(kPilscNumSpecies, 0);
    s.counterArrays.copyNumSpecies[kPilscStateU] = 1;   // one un-phosphorylated A
    s.counterArrays.copyNumSpecies[kPilscStateP] = 0;   // no phosphorylated A yet
    s.counterArrays.copyNumSpecies[kPilscIlIface] = 100; // implicit lipid sites
    s.counterArrays.singleDouble = std::vector<int>(kPilscNumSpecies, 0);
    s.counterArrays.implicitDouble = std::vector<bool>(kPilscNumSpecies, false);
    s.counterArrays.canDissociate = std::vector<bool>(kPilscNumSpecies, false);
    s.counterArrays.bindPairList.resize(kPilscNumSpecies);
    s.counterArrays.bindPairListIL2D.resize(kPilscNumSpecies);
    s.counterArrays.bindPairListIL3D.resize(kPilscNumSpecies);
    s.counterArrays.nBoundPairs = std::vector<int>(4, 0);
    s.counterArrays.proPairlist = std::vector<int>(4, 0);
    s.counterArrays.events3D = std::vector<int>(s.counterArrays.eventArraySize, 0);
    s.counterArrays.events2D = std::vector<int>(s.counterArrays.eventArraySize, 0);
    s.counterArrays.events3Dto2D = std::vector<int>(s.counterArrays.eventArraySize, 0);

    // rxnItr = {reaction index 0, rate index 0, forward (not back) reaction}
    s.rxnItr = { { 0, 0, 0 } };

    return s;
}

/*! \brief Dump the interesting parts of a molecule's state to stderr. */
void pilsc_report(const PilscSystem& s, const char* label)
{
    const Molecule::Iface& iface = s.moleculeList[0].interfaceList[0];
    std::cerr << "    " << label << ": A.iface[0] stateIden='"
              << (iface.stateIden == '\0' ? '-' : iface.stateIden) << "'"
              << " absIndex=" << iface.index
              << " relStateIndex=" << iface.stateIndex
              << " | copyNum[~U]=" << s.counterArrays.copyNumSpecies[kPilscStateU]
              << " copyNum[~P]=" << s.counterArrays.copyNumSpecies[kPilscStateP]
              << " | observable(Aphos)=" << s.observablesList.at("Aphos") << '\n';
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
              << "  Scenario:      A(a~U) + implicit lipid -> A(a~P) inside a\n"
              << "                 100x100x100 reflecting water box.\n"
              << "  Pass criteria: interface state flips to ~P, copy numbers are\n"
              << "                 moved from species ~U to species ~P, and the\n"
              << "                 molecules are flagged as already propagated.\n";

    pilsc_ensure_rng();
    PilscSystem s = pilsc_make_system(/*isSphere=*/false);

    pilsc_report(s, "before");
    // Sanity check on the fixture itself (guards against a broken test setup).
    EXPECT_EQ(s.moleculeList[0].interfaceList[0].stateIden, 'U')
        << "fixture should start with protein A in state ~U";
    EXPECT_FALSE(s.membrane.isSphere) << "this test must exercise the box branch";

    std::cerr << "  Calling perform_implicitlipid_state_change...\n";
    perform_implicitlipid_state_change(/*stateChangeIface=*/0, /*facilitatorIface=*/0, s.rxnItr,
        s.moleculeList[0], s.moleculeList[1], s.complexList[0], s.complexList[1],
        s.counterArrays, s.params, s.forwardRxns, s.backRxns, s.moleculeList, s.complexList,
        s.molTemplateList, s.observablesList, s.membrane);
    pilsc_report(s, "after ");

    // The reacting interface must now report the product state ~P.
    EXPECT_EQ(s.moleculeList[0].interfaceList[0].stateIden, 'P')
        << "box branch should have changed the interface identity to 'P'";
    EXPECT_EQ(s.moleculeList[0].interfaceList[0].index, kPilscStateP)
        << "box branch should have set the absolute state index to the product index";

    // Copy-number bookkeeping: one molecule leaves ~U and enters ~P.
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kPilscStateU], 0)
        << "copy number of the reactant species (~U) should be decremented";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kPilscStateP], 1)
        << "copy number of the product species (~P) should be incremented";

    // The reacting molecule/complex must not be moved again this timestep.
    EXPECT_NE(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "state-change molecule should have its trajStatus updated";
    EXPECT_NE(static_cast<int>(s.complexList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "state-change complex should have its trajStatus updated";

    // The coupled observable should have registered the event exactly once.
    EXPECT_EQ(s.observablesList.at("Aphos"), 1)
        << "the observable coupled to this reaction should be incremented once";
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
              << "                 same ~U -> ~P state change and bookkeeping.\n";

    pilsc_ensure_rng();
    PilscSystem s = pilsc_make_system(/*isSphere=*/true);

    pilsc_report(s, "before");
    EXPECT_TRUE(s.membrane.isSphere) << "this test must exercise the sphere branch";
    EXPECT_GT(s.membrane.sphereR, 0.0) << "sphere radius must be positive for the sphere branch";

    std::cerr << "  Calling perform_implicitlipid_state_change...\n";
    perform_implicitlipid_state_change(/*stateChangeIface=*/0, /*facilitatorIface=*/0, s.rxnItr,
        s.moleculeList[0], s.moleculeList[1], s.complexList[0], s.complexList[1],
        s.counterArrays, s.params, s.forwardRxns, s.backRxns, s.moleculeList, s.complexList,
        s.molTemplateList, s.observablesList, s.membrane);
    pilsc_report(s, "after ");

    // Same expectations as for the box branch: the dispatcher only chooses the
    // implementation, the logical outcome must be identical.
    EXPECT_EQ(s.moleculeList[0].interfaceList[0].stateIden, 'P')
        << "sphere branch should have changed the interface identity to 'P'";
    EXPECT_EQ(s.moleculeList[0].interfaceList[0].index, kPilscStateP)
        << "sphere branch should have set the absolute state index to the product index";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kPilscStateU], 0)
        << "copy number of the reactant species (~U) should be decremented";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kPilscStateP], 1)
        << "copy number of the product species (~P) should be incremented";
    EXPECT_NE(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "state-change molecule should have its trajStatus updated";
    EXPECT_EQ(s.observablesList.at("Aphos"), 1)
        << "the observable coupled to this reaction should be incremented once";

    // The molecule must still be inside (or on) the spherical volume afterwards.
    const Coord& com = s.moleculeList[0].comCoord;
    const double radial = std::sqrt(com.x * com.x + com.y * com.y + com.z * com.z);
    std::cerr << "    final radial distance of A = " << radial
              << " (sphere radius " << s.membrane.sphereR << ")\n";
    EXPECT_LE(radial, s.membrane.sphereR + 1e-6)
        << "the state-changing molecule must remain inside the spherical boundary";
}

// -----------------------------------------------------------------------------
// Test 3: invariants -- the facilitator (implicit lipid) and the boundary
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
    const int ilCopyNumBefore = s.counterArrays.copyNumSpecies[kPilscIlIface];

    std::cerr << "  Calling perform_implicitlipid_state_change...\n";
    perform_implicitlipid_state_change(/*stateChangeIface=*/0, /*facilitatorIface=*/0, s.rxnItr,
        s.moleculeList[0], s.moleculeList[1], s.complexList[0], s.complexList[1],
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
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kPilscIlIface], ilCopyNumBefore)
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
              << "                 (population is only moved from ~U to ~P).\n";

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
            s.moleculeList[0], s.moleculeList[1], s.complexList[0], s.complexList[1],
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
        EXPECT_EQ(s.moleculeList[0].interfaceList[0].stateIden, 'P')
            << "the " << (isSphere ? "sphere" : "box")
            << " branch should have produced the phosphorylated state";
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