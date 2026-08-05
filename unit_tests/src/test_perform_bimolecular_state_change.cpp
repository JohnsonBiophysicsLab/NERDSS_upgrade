/*! \file test_perform_bimolecular_state_change.cpp
 *
 * ### Unit test for src/reactions/perform_bimolecular_state_change.cpp
 *
 * The file under test contains exactly one function:
 *
 *     void perform_bimolecular_state_change(int stateChangeIface, int facilitatorIface,
 *                                           std::array<int,3>& rxnItr,
 *                                           Molecule& stateChangeMol, Molecule& facilitatorMol,
 *                                           Complex& stateChangeCom, Complex& facilitatorCom,
 *                                           copyCounters&, const Parameters&,
 *                                           std::vector<ForwardRxn>&, std::vector<BackRxn>&,
 *                                           std::vector<Molecule>&, std::vector<Complex>&,
 *                                           std::vector<MolTemplate>&,
 *                                           std::map<std::string,int>&, Membrane&)
 *
 * It is a pure dispatcher: if the Membrane describes a spherical boundary the call is
 * forwarded to perform_bimolecular_state_change_sphere(), otherwise it is forwarded to
 * perform_bimolecular_state_change_box().  Both branches must carry out the very same
 * chemistry, namely: flip the state of `stateChangeIface` on `stateChangeMol` from the
 * reactant state of the reaction to the product state, leave the facilitator interface
 * untouched, and keep the species copy-number bookkeeping consistent.
 *
 * The tests below therefore:
 *   1. build a small but completely self-consistent reaction system
 *      ( A(a~U) + B(b) <-> A(a~P) + B(b), a bimolecular state change ),
 *   2. call the dispatcher with a box Membrane and with a sphere Membrane,
 *   3. call the dispatcher with the forward reaction and with the back reaction,
 *   4. verify the interface state (identity, absolute index, relative state index),
 *      that the facilitator is unchanged, and that the copy numbers of the two states
 *      were transferred (one decrement / one increment, total conserved).
 *
 * Verbose progress information is written to stderr so that a reader of the test log can
 * see which source file / function is exercised and what each assertion checks.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

namespace {

// -----------------------------------------------------------------------------
// Absolute (system wide) interface-state indices used by this test.
//
//   0 -> A(a~U)  : reactant state of the state-changing interface
//   1 -> A(a~P)  : product state of the state-changing interface
//   2 -> B(b)    : the facilitator interface (stateless)
// -----------------------------------------------------------------------------
constexpr int kPbscAbsStateU { 0 };
constexpr int kPbscAbsStateP { 1 };
constexpr int kPbscAbsIfaceB { 2 };
constexpr int kPbscNumSpecies { 3 };

//! Relative index of state 'U' and 'P' inside Interface::stateList of A's interface "a".
constexpr int kPbscRelStateU { 0 };
constexpr int kPbscRelStateP { 1 };

/*! \brief Container holding every object the function under test needs.
 *
 * Everything is kept in one struct so that the individual tests stay short and the
 * (fairly large) construction code lives in exactly one place.
 */
struct PbscSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
};

/*! \brief Make sure the global GSL RNG exists.
 *
 * Some of the boundary-condition helpers that the state-change routines may touch draw
 * random numbers through the global `r`.  gtest_main.cpp only defines it (possibly as
 * nullptr), so guard against a null generator here.  `r` is *declared* extern by
 * math/rand_gsl.hpp - we never define it.
 */
void pbsc_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] global GSL RNG pointer was null -> calling srand_gsl(1)\n";
        srand_gsl(1);
    }
}

/*! \brief Build the two MolTemplates (A with a two-state interface, B with one stateless
 *         interface) and set the static bookkeeping members they rely on.
 */
void pbsc_build_templates(std::vector<MolTemplate>& molTemplateList)
{
    molTemplateList.clear();

    // ---------------- Molecule type 0: "A", interface "a" with states U and P ----------
    MolTemplate tempA;
    tempA.molName = "A";
    tempA.molTypeIndex = 0;
    tempA.copies = 1;
    tempA.mass = 1.0;
    tempA.radius = 1.0;
    tempA.comCoord = Coord(0.0, 0.0, 0.0);
    tempA.D = Coord(10.0, 10.0, 10.0);
    tempA.Dr = Coord(0.1, 0.1, 0.1);
    tempA.isPoint = true; // interface sits on the COM -> trivial geometry
    tempA.checkOverlap = false;
    tempA.countTransition = false;

    Interface ifaceA;
    ifaceA.index = 0;
    ifaceA.name = "a";
    ifaceA.iCoord = Coord(0.0, 0.0, 0.0);
    Interface::State stateU('U', kPbscAbsStateU); // (identity, absolute index)
    stateU.ifaceAndStateName = "a~U";
    Interface::State stateP('P', kPbscAbsStateP);
    stateP.ifaceAndStateName = "a~P";
    ifaceA.stateList.push_back(stateU);
    ifaceA.stateList.push_back(stateP);
    tempA.interfaceList.push_back(ifaceA);
    tempA.ifacesWithStates.push_back(0);
    molTemplateList.push_back(tempA);

    // ---------------- Molecule type 1: "B", stateless interface "b" --------------------
    MolTemplate tempB;
    tempB.molName = "B";
    tempB.molTypeIndex = 1;
    tempB.copies = 1;
    tempB.mass = 1.0;
    tempB.radius = 1.0;
    tempB.comCoord = Coord(0.0, 0.0, 0.0);
    tempB.D = Coord(10.0, 10.0, 10.0);
    tempB.Dr = Coord(0.1, 0.1, 0.1);
    tempB.isPoint = true;
    tempB.checkOverlap = false;
    tempB.countTransition = false;

    Interface ifaceB;
    ifaceB.index = 0;
    ifaceB.name = "b";
    ifaceB.iCoord = Coord(0.0, 0.0, 0.0);
    Interface::State stateB('\0', kPbscAbsIfaceB);
    stateB.ifaceAndStateName = "b";
    ifaceB.stateList.push_back(stateB);
    tempB.interfaceList.push_back(ifaceB);
    molTemplateList.push_back(tempB);

    // ---------------- statics the simulation classes expect to be initialized ---------
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = { 1, 1 };
    // absolute state index -> relative interface index (abs 0,1 -> iface 0 of A, abs 2 -> iface 0 of B)
    MolTemplate::absToRelIface = { 0, 0, 0 };
    Interface::State::totalNumOfStates = kPbscNumSpecies;
}

/*! \brief Build the forward (U -> P) and back (P -> U) bimolecular state-change reactions. */
void pbsc_build_reactions(std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns)
{
    forwardRxns.clear();
    backRxns.clear();

    // The three interfaces taking part in the reaction.
    RxnIface reactStateChange("a~U", /*molTypeIndex*/ 0, kPbscAbsStateU, /*relIfaceIndex*/ 0, 'U', false);
    RxnIface prodStateChange("a~P", /*molTypeIndex*/ 0, kPbscAbsStateP, /*relIfaceIndex*/ 0, 'P', false);
    RxnIface facilitator("b", /*molTypeIndex*/ 1, kPbscAbsIfaceB, /*relIfaceIndex*/ 0, '\0', false);

    RxnBase::RateState oneRate;
    oneRate.rate = 10.0;
    oneRate.prob = 1.0;
    oneRate.otherIfaceLists.clear(); // no ancillary interfaces required

    // ------------------------------- forward reaction ---------------------------------
    ForwardRxn fwd;
    fwd.rxnType = ReactionType::biMolStateChange;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.hasStateChange = true;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.isObserved = false; // keeps observablesList out of the picture
    fwd.isOnMem = false;
    fwd.isSymmetric = false;
    fwd.bindRadius = 1.0;
    fwd.rxnLabel = "pbscForward";
    fwd.productName = "A(a~P)";
    // element 0 is always the facilitator, element 1 the state-changing interface
    // (see ParsedRxn::determine_reactants_products, which always inserts the
    // facilitator at the front of both lists for biMolStateChange reactions --
    // perform_bimolecular_state_change_box/_sphere read productListNew[1] as
    // the state-change product).
    fwd.reactantListNew = { facilitator, reactStateChange };
    fwd.productListNew = { facilitator, prodStateChange };
    fwd.stateChangeIface = std::make_pair(reactStateChange, prodStateChange);
    fwd.intReactantList = { kPbscAbsIfaceB, kPbscAbsStateU };
    fwd.intProductList = { kPbscAbsIfaceB, kPbscAbsStateP };
    fwd.rateList = { oneRate };
    forwardRxns.push_back(fwd);

    // -------------------------------- back reaction -----------------------------------
    BackRxn back;
    back.rxnType = ReactionType::biMolStateChange;
    back.absRxnIndex = 1;
    back.relRxnIndex = 0;
    back.conjForwardRxnIndex = 0;
    back.hasStateChange = true;
    back.isObserved = false;
    back.isOnMem = false;
    back.isSymmetric = false;
    back.rxnLabel = "pbscBack";
    back.reactantListNew = { facilitator, prodStateChange };
    back.productListNew = { facilitator, reactStateChange };
    back.stateChangeIface = std::make_pair(prodStateChange, reactStateChange);
    back.intReactantList = { kPbscAbsIfaceB, kPbscAbsStateP };
    back.intProductList = { kPbscAbsIfaceB, kPbscAbsStateU };
    back.rateList = { oneRate };
    backRxns.push_back(back);
}

/*! \brief Assemble a full mini-system.
 *
 * \param[in] useSphere      true  -> spherical Membrane  (dispatcher must take the sphere branch)
 *                           false -> rectangular water box (dispatcher must take the box branch)
 * \param[in] startingState  'U' (ready for the forward rxn) or 'P' (ready for the back rxn)
 */
PbscSystem pbsc_build_system(bool useSphere, char startingState)
{
    PbscSystem sys;

    // ------------------------------- simulation parameters ----------------------------
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = kPbscNumSpecies;
    sys.params.numTotalComplex = 2;
    sys.params.timeStep = 0.1;
    sys.params.nItr = 1;
    sys.params.name = "pbsc_unit_test";
    Parameters::dt = sys.params.timeStep;
    Parameters::lastUpdateTransition.assign(2, 0);

    // ------------------------------------ boundary ------------------------------------
    if (useSphere) {
        sys.membrane.isSphere = true;
        sys.membrane.isBox = false;
        sys.membrane.sphereR = 100.0;
        sys.membrane.sphereVol = (4.0 / 3.0) * M_PI * std::pow(sys.membrane.sphereR, 3);
        // a bounding water box is still filled in by the real code for spheres
        sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
    } else {
        sys.membrane.isSphere = false;
        sys.membrane.isBox = true;
        sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    }
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";
    sys.membrane.implicitLipid = false;
    sys.membrane.hasCompartment = false;
    sys.membrane.TwoD = false;

    // ------------------------------ templates and reactions ---------------------------
    pbsc_build_templates(sys.molTemplateList);
    pbsc_build_reactions(sys.forwardRxns, sys.backRxns);

    // ------------------------------------ molecules -----------------------------------
    // Both molecules sit well inside either boundary and one binding radius apart.
    const Coord posA { 0.0, 0.0, -40.0 };
    const Coord posB { 1.0, 0.0, -40.0 };

    const int startAbsIndex { (startingState == 'U') ? kPbscAbsStateU : kPbscAbsStateP };
    const int startRelIndex { (startingState == 'U') ? kPbscRelStateU : kPbscRelStateP };

    Molecule molA;
    molA.index = 0;
    molA.id = 0;
    molA.myComIndex = 0;
    molA.molTypeIndex = 0;
    molA.mass = 1.0;
    molA.isLipid = false;
    molA.comCoord = posA;
    molA.trajStatus = TrajStatus::none;
    {
        Molecule::Iface ia;
        ia.coord = posA;
        ia.stateIden = startingState;
        ia.stateIndex = startRelIndex;
        ia.index = startAbsIndex;
        ia.relIndex = 0;
        ia.molTypeIndex = 0;
        ia.isBound = false;
        molA.interfaceList.push_back(ia);
    }
    molA.freelist = { 0 };

    Molecule molB;
    molB.index = 1;
    molB.id = 1;
    molB.myComIndex = 1;
    molB.molTypeIndex = 1;
    molB.mass = 1.0;
    molB.isLipid = false;
    molB.comCoord = posB;
    molB.trajStatus = TrajStatus::none;
    {
        Molecule::Iface ib;
        ib.coord = posB;
        ib.stateIden = '\0';
        ib.stateIndex = 0;
        ib.index = kPbscAbsIfaceB;
        ib.relIndex = 0;
        ib.molTypeIndex = 1;
        ib.isBound = false;
        molB.interfaceList.push_back(ib);
    }
    molB.freelist = { 0 };

    sys.moleculeList.push_back(molA);
    sys.moleculeList.push_back(molB);
    Molecule::numberOfMolecules = 2;
    Molecule::emptyMolList.clear();

    // ------------------------------------ complexes -----------------------------------
    // Complex(index, memberMolecule, template) fills in memberList/D/Dr/radius/numEachMol.
    Complex comA(0, sys.moleculeList[0], sys.molTemplateList[0]);
    comA.index = 0;
    comA.id = 0;
    comA.trajStatus = TrajStatus::none;
    comA.ncross = 0;

    Complex comB(1, sys.moleculeList[1], sys.molTemplateList[1]);
    comB.index = 1;
    comB.id = 1;
    comB.trajStatus = TrajStatus::none;
    comB.ncross = 0;

    sys.complexList.push_back(comA);
    sys.complexList.push_back(comB);
    Complex::numberOfComplexes = 2;
    Complex::currNumberMolTypes = 2;
    Complex::currNumberComTypes = 2;
    Complex::emptyComList.clear();

    // ------------------------------- species bookkeeping ------------------------------
    sys.counterArrays.copyNumSpecies.assign(kPbscNumSpecies, 0);
    sys.counterArrays.copyNumSpecies[startAbsIndex] = 1; // the single copy of A in its start state
    sys.counterArrays.copyNumSpecies[kPbscAbsIfaceB] = 1; // the single copy of B
    sys.counterArrays.singleDouble.assign(kPbscNumSpecies, 0);
    sys.counterArrays.implicitDouble.assign(kPbscNumSpecies, false);
    sys.counterArrays.canDissociate.assign(kPbscNumSpecies, false);
    sys.counterArrays.bindPairList.assign(kPbscNumSpecies, std::vector<int>());
    sys.counterArrays.bindPairListIL2D.assign(kPbscNumSpecies, std::vector<int>());
    sys.counterArrays.bindPairListIL3D.assign(kPbscNumSpecies, std::vector<int>());
    sys.counterArrays.nBoundPairs.assign(MolTemplate::numMolTypes * MolTemplate::numMolTypes, 0);
    sys.counterArrays.proPairlist.assign(MolTemplate::numMolTypes * MolTemplate::numMolTypes, 0);
    sys.counterArrays.events3D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events2D.assign(sys.counterArrays.eventArraySize, 0);
    sys.counterArrays.events3Dto2D.assign(sys.counterArrays.eventArraySize, 0);

    return sys;
}

/*! \brief Convenience wrapper: run the function under test on a prepared system.
 *
 * \param[in,out] sys      the mini system (molecule 0 = state changer, molecule 1 = facilitator)
 * \param[in]     isBackRxn  false -> use forwardRxns[0], true -> use backRxns[0]
 */
void pbsc_invoke(PbscSystem& sys, bool isBackRxn)
{
    // rxnItr = { reaction index, rate index, 0 = forward / 1 = back }
    std::array<int, 3> rxnItr { 0, 0, isBackRxn ? 1 : 0 };

    perform_bimolecular_state_change(/*stateChangeIface*/ 0, /*facilitatorIface*/ 0, rxnItr,
        sys.moleculeList[0], sys.moleculeList[1], sys.complexList[0], sys.complexList[1],
        sys.counterArrays, sys.params, sys.forwardRxns, sys.backRxns, sys.moleculeList,
        sys.complexList, sys.molTemplateList, sys.observablesList, sys.membrane);
}

/*! \brief Print the interesting bits of the system state to stderr. */
void pbsc_dump(const PbscSystem& sys, const char* tag)
{
    const Molecule::Iface& ia = sys.moleculeList[0].interfaceList[0];
    std::cerr << "  [" << tag << "] A(a) stateIden='" << (ia.stateIden ? ia.stateIden : '0')
              << "' absIndex=" << ia.index << " relStateIndex=" << ia.stateIndex << '\n'
              << "  [" << tag << "] copyNumSpecies = { a~U=" << sys.counterArrays.copyNumSpecies[kPbscAbsStateU]
              << ", a~P=" << sys.counterArrays.copyNumSpecies[kPbscAbsStateP]
              << ", b=" << sys.counterArrays.copyNumSpecies[kPbscAbsIfaceB] << " }\n"
              << "  [" << tag << "] trajStatus: stateChangeMol="
              << static_cast<int>(sys.moleculeList[0].trajStatus) << " facilitatorMol="
              << static_cast<int>(sys.moleculeList[1].trajStatus) << " stateChangeCom="
              << static_cast<int>(sys.complexList[0].trajStatus) << " facilitatorCom="
              << static_cast<int>(sys.complexList[1].trajStatus) << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: box boundary, forward reaction  A(a~U) + B(b) -> A(a~P) + B(b)
// -----------------------------------------------------------------------------
void test_pbsc_box_forward_state_change()
{
    std::cerr << "\n[TEST] test_pbsc_box_forward_state_change\n"
              << "  Source file:   src/reactions/perform_bimolecular_state_change.cpp\n"
              << "  Function:      perform_bimolecular_state_change (box branch)\n"
              << "  Scenario:      rectangular water box, forward rxn (rxnItr[2]==0),\n"
              << "                 A's interface 'a' starts in state 'U'.\n"
              << "  Pass criteria: interface 'a' ends in state 'P' with the product's\n"
              << "                 absolute index, facilitator untouched, and the state\n"
              << "                 copy numbers are transferred U -> P.\n";

    pbsc_ensure_rng();
    PbscSystem sys = pbsc_build_system(/*useSphere*/ false, /*startingState*/ 'U');
    pbsc_dump(sys, "before");

    std::cerr << "  Calling perform_bimolecular_state_change...\n";
    pbsc_invoke(sys, /*isBackRxn*/ false);
    pbsc_dump(sys, "after");

    // ---- the state-changing interface must now be the product state -----------------
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'P')
        << "state identity of A(a) should have flipped from 'U' to 'P'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kPbscAbsStateP)
        << "absolute interface-state index of A(a) should be the product index";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, kPbscRelStateP)
        << "relative state index of A(a) should point at stateList entry for 'P'";
    EXPECT_FALSE(sys.moleculeList[0].interfaceList[0].isBound)
        << "a state change must not create an interaction";

    // ---- the facilitator must be completely unaffected -------------------------------
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].stateIden, '\0')
        << "facilitator B(b) is stateless and must not change identity";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].index, kPbscAbsIfaceB)
        << "facilitator B(b) absolute index must be unchanged";
    EXPECT_FALSE(sys.moleculeList[1].interfaceList[0].isBound)
        << "facilitator B(b) must remain unbound";

    // ---- species bookkeeping ---------------------------------------------------------
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateU], 0)
        << "the copy number of the reactant state a~U should have been decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateP], 1)
        << "the copy number of the product state a~P should have been incremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsIfaceB], 1)
        << "the facilitator copy number must not change";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateU]
            + sys.counterArrays.copyNumSpecies[kPbscAbsStateP],
        1)
        << "the total number of A molecules must be conserved by a state change";
}

// -----------------------------------------------------------------------------
// Test 2: sphere boundary, forward reaction (the other dispatcher branch)
// -----------------------------------------------------------------------------
void test_pbsc_sphere_forward_state_change()
{
    std::cerr << "\n[TEST] test_pbsc_sphere_forward_state_change\n"
              << "  Source file:   src/reactions/perform_bimolecular_state_change.cpp\n"
              << "  Function:      perform_bimolecular_state_change (sphere branch)\n"
              << "  Scenario:      Membrane::isSphere == true, forward rxn.\n"
              << "  Pass criteria: exactly the same chemistry as the box branch, i.e.\n"
              << "                 'U' -> 'P' with consistent copy numbers.\n";

    pbsc_ensure_rng();
    PbscSystem sys = pbsc_build_system(/*useSphere*/ true, /*startingState*/ 'U');
    EXPECT_TRUE(sys.membrane.isSphere) << "sanity: the sphere branch must be selected";
    pbsc_dump(sys, "before");

    std::cerr << "  Calling perform_bimolecular_state_change...\n";
    pbsc_invoke(sys, /*isBackRxn*/ false);
    pbsc_dump(sys, "after");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'P')
        << "sphere branch should also flip A(a) from 'U' to 'P'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kPbscAbsStateP)
        << "sphere branch should set the product absolute state index";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, kPbscRelStateP)
        << "sphere branch should set the product relative state index";

    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].stateIden, '\0')
        << "facilitator must be untouched on the sphere as well";

    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateU], 0)
        << "a~U copy number should have been decremented (sphere branch)";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateP], 1)
        << "a~P copy number should have been incremented (sphere branch)";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsIfaceB], 1)
        << "facilitator copy number must not change (sphere branch)";
}

// -----------------------------------------------------------------------------
// Test 3: box boundary, back reaction  A(a~P) + B(b) -> A(a~U) + B(b)
// -----------------------------------------------------------------------------
void test_pbsc_box_back_state_change()
{
    std::cerr << "\n[TEST] test_pbsc_box_back_state_change\n"
              << "  Source file:   src/reactions/perform_bimolecular_state_change.cpp\n"
              << "  Function:      perform_bimolecular_state_change (box branch, back rxn)\n"
              << "  Scenario:      rxnItr[2]==1 so the reaction is taken from backRxns,\n"
              << "                 A's interface 'a' starts in state 'P'.\n"
              << "  Pass criteria: interface 'a' ends in state 'U' and copy numbers move\n"
              << "                 P -> U.\n";

    pbsc_ensure_rng();
    PbscSystem sys = pbsc_build_system(/*useSphere*/ false, /*startingState*/ 'P');
    pbsc_dump(sys, "before");

    std::cerr << "  Calling perform_bimolecular_state_change (back reaction)...\n";
    pbsc_invoke(sys, /*isBackRxn*/ true);
    pbsc_dump(sys, "after");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'U')
        << "the back reaction should flip A(a) from 'P' back to 'U'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kPbscAbsStateU)
        << "absolute index should be the back reaction's product (a~U)";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, kPbscRelStateU)
        << "relative state index should point at stateList entry for 'U'";

    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateP], 0)
        << "a~P copy number should have been decremented by the back reaction";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateU], 1)
        << "a~U copy number should have been incremented by the back reaction";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsIfaceB], 1)
        << "facilitator copy number must not change";
}

// -----------------------------------------------------------------------------
// Test 4: sphere boundary, back reaction (sphere branch, reverse direction)
// -----------------------------------------------------------------------------
void test_pbsc_sphere_back_state_change()
{
    std::cerr << "\n[TEST] test_pbsc_sphere_back_state_change\n"
              << "  Source file:   src/reactions/perform_bimolecular_state_change.cpp\n"
              << "  Function:      perform_bimolecular_state_change (sphere branch, back rxn)\n"
              << "  Scenario:      spherical boundary + rxnItr[2]==1, start state 'P'.\n"
              << "  Pass criteria: 'P' -> 'U' with copy numbers moved P -> U.\n";

    pbsc_ensure_rng();
    PbscSystem sys = pbsc_build_system(/*useSphere*/ true, /*startingState*/ 'P');
    pbsc_dump(sys, "before");

    std::cerr << "  Calling perform_bimolecular_state_change (back reaction, sphere)...\n";
    pbsc_invoke(sys, /*isBackRxn*/ true);
    pbsc_dump(sys, "after");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'U')
        << "sphere branch back reaction should flip A(a) from 'P' to 'U'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kPbscAbsStateU)
        << "sphere branch back reaction should set absolute index of a~U";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateP], 0)
        << "a~P copy number should have been decremented (sphere back rxn)";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateU], 1)
        << "a~U copy number should have been incremented (sphere back rxn)";
}

// -----------------------------------------------------------------------------
// Test 5: round trip - forward followed by back must restore the initial state.
//         This checks that the dispatcher (and the routine it delegates to) leaves
//         the system in a self-consistent, reusable condition.
// -----------------------------------------------------------------------------
void test_pbsc_round_trip_restores_initial_state()
{
    std::cerr << "\n[TEST] test_pbsc_round_trip_restores_initial_state\n"
              << "  Source file:   src/reactions/perform_bimolecular_state_change.cpp\n"
              << "  Function:      perform_bimolecular_state_change (box branch, twice)\n"
              << "  Scenario:      apply the forward rxn, then the back rxn, to the same\n"
              << "                 molecule pair.\n"
              << "  Pass criteria: the interface state and every copy number return to\n"
              << "                 their initial values.\n";

    pbsc_ensure_rng();
    PbscSystem sys = pbsc_build_system(/*useSphere*/ false, /*startingState*/ 'U');

    // Snapshot of the initial condition.
    const char initialIden { sys.moleculeList[0].interfaceList[0].stateIden };
    const int initialAbs { sys.moleculeList[0].interfaceList[0].index };
    const int initialU { sys.counterArrays.copyNumSpecies[kPbscAbsStateU] };
    const int initialP { sys.counterArrays.copyNumSpecies[kPbscAbsStateP] };
    const int initialB { sys.counterArrays.copyNumSpecies[kPbscAbsIfaceB] };
    pbsc_dump(sys, "initial");

    std::cerr << "  Step 1: forward reaction U -> P\n";
    pbsc_invoke(sys, /*isBackRxn*/ false);
    pbsc_dump(sys, "after forward");
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'P')
        << "intermediate state after the forward reaction should be 'P'";

    std::cerr << "  Step 2: back reaction P -> U\n";
    pbsc_invoke(sys, /*isBackRxn*/ true);
    pbsc_dump(sys, "after back");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, initialIden)
        << "state identity should return to its initial value after a round trip";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, initialAbs)
        << "absolute state index should return to its initial value after a round trip";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateU], initialU)
        << "a~U copy number should return to its initial value";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsStateP], initialP)
        << "a~P copy number should return to its initial value";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscAbsIfaceB], initialB)
        << "facilitator copy number should never have changed";
}

// -----------------------------------------------------------------------------
// Test 6: structural invariants - the routine must not create/destroy molecules,
//         complexes, or complex membership while performing a state change.
// -----------------------------------------------------------------------------
void test_pbsc_preserves_system_structure()
{
    std::cerr << "\n[TEST] test_pbsc_preserves_system_structure\n"
              << "  Source file:   src/reactions/perform_bimolecular_state_change.cpp\n"
              << "  Function:      perform_bimolecular_state_change\n"
              << "  Scenario:      forward rxn in a box; inspect list sizes / membership.\n"
              << "  Pass criteria: molecule and complex counts, complex membership and\n"
              << "                 the molecule->complex mapping are all unchanged; no\n"
              << "                 molecule or complex is flagged empty.\n";

    pbsc_ensure_rng();
    PbscSystem sys = pbsc_build_system(/*useSphere*/ false, /*startingState*/ 'U');

    const size_t nMolBefore { sys.moleculeList.size() };
    const size_t nComBefore { sys.complexList.size() };

    std::cerr << "  Calling perform_bimolecular_state_change...\n";
    pbsc_invoke(sys, /*isBackRxn*/ false);

    std::cerr << "  moleculeList.size()=" << sys.moleculeList.size()
              << " complexList.size()=" << sys.complexList.size() << '\n';

    EXPECT_EQ(sys.moleculeList.size(), nMolBefore)
        << "a state change must not add or remove molecules";
    EXPECT_EQ(sys.complexList.size(), nComBefore)
        << "a state change must not add or remove complexes";

    EXPECT_FALSE(sys.moleculeList[0].isEmpty) << "state-change molecule must still exist";
    EXPECT_FALSE(sys.moleculeList[1].isEmpty) << "facilitator molecule must still exist";
    EXPECT_FALSE(sys.complexList[0].isEmpty) << "state-change complex must still exist";
    EXPECT_FALSE(sys.complexList[1].isEmpty) << "facilitator complex must still exist";

    EXPECT_EQ(sys.moleculeList[0].myComIndex, 0)
        << "the state-change molecule must still belong to complex 0";
    EXPECT_EQ(sys.moleculeList[1].myComIndex, 1)
        << "the facilitator molecule must still belong to complex 1";

    ASSERT_EQ(sys.complexList[0].memberList.size(), 1u)
        << "complex 0 should still contain exactly one member";
    ASSERT_EQ(sys.complexList[1].memberList.size(), 1u)
        << "complex 1 should still contain exactly one member";
    EXPECT_EQ(sys.complexList[0].memberList[0], 0)
        << "complex 0's member should still be molecule 0";
    EXPECT_EQ(sys.complexList[1].memberList[0], 1)
        << "complex 1's member should still be molecule 1";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* helper is run inside its own TEST so the
// framework reports individual results while all of them still execute.
// -----------------------------------------------------------------------------
TEST(PerformBimolecularStateChange, BoxForwardStateChange) { test_pbsc_box_forward_state_change(); }
TEST(PerformBimolecularStateChange, SphereForwardStateChange) { test_pbsc_sphere_forward_state_change(); }
TEST(PerformBimolecularStateChange, BoxBackStateChange) { test_pbsc_box_back_state_change(); }
TEST(PerformBimolecularStateChange, SphereBackStateChange) { test_pbsc_sphere_back_state_change(); }
TEST(PerformBimolecularStateChange, RoundTripRestoresInitialState) { test_pbsc_round_trip_restores_initial_state(); }
TEST(PerformBimolecularStateChange, PreservesSystemStructure) { test_pbsc_preserves_system_structure(); }