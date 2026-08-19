/*! \file test_check_for_unimolstatechange_reactions.cpp
 *
 * ### Unit test for src/reactions/check_for_unimolstatechange_reactions.cpp
 *
 * Function under test:
 * \code
 * void check_for_unimolstatechange_reactions(unsigned simItr, Parameters&,
 *          std::vector<Molecule>&, std::vector<Complex>&, SimulVolume&,
 *          const std::vector<ForwardRxn>&, const std::vector<BackRxn>&,
 *          const std::vector<CreateDestructRxn>&, std::vector<MolTemplate>&,
 *          std::map<std::string,int>&, copyCounters&, Membrane&)
 * \endcode
 *
 * The routine loops over every Molecule in the system and, for every interface
 * that owns a state list, tries to fire a unimolecular state-change reaction.
 * Whether the reaction actually fires is decided by
 *      prob = 1 - exp(-rate * timeStep * 1e-6)   compared against rand_gsl().
 *
 * That gives us two fully deterministic regimes that we can test:
 *   * rate == 0                -> prob == 0 -> the reaction can never fire.
 *   * rate == 1e12 (huge)      -> prob == 1 -> the reaction always fires.
 *
 * With those two regimes we exercise:
 *   - the early-out guards (isEmpty, isGhosted, params.hasUniMolStateChange,
 *     empty ifacesWithStates, TrajStatus != none, bound interface),
 *   - the explicit-molecule state change path (state flip, copyNumSpecies
 *     bookkeeping, TrajStatus propagation, observable counting),
 *   - the implicit-lipid path (copyNumSpecies plus
 *     Membrane::numberOfFreeLipidsEachState bookkeeping and the temporary
 *     interface-index shifting that must be undone again).
 *
 * Everything is printed to stderr so the reader can follow which behaviour is
 * being checked and what the pass criterion is.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

namespace {

// -----------------------------------------------------------------------------
// Constants describing the tiny reaction network used by every test:
//
//   molecule "A" has a single interface "a" that can be in state 'A' or 'B'
//   forward reaction (index 0): a~A -> a~B   (absolute state 0 -> 1)
//   back    reaction (index 0): a~B -> a~A   (absolute state 1 -> 0)
// -----------------------------------------------------------------------------
constexpr int kCfuscStateAAbs = 0; //!< absolute interface-state index of a~A
constexpr int kCfuscStateBAbs = 1; //!< absolute interface-state index of a~B

//! A rate so large that prob == 1.0 for a 1 us time step (reaction always fires)
constexpr double kCfuscHugeRate = 1.0e12;

/*! \brief Container holding an entire (minimal) simulation system.
 *
 * Keeping everything in one struct means each test can build a pristine system
 * with one call and then poke at exactly the field it wants to vary.
 */
struct CfuscSystem {
    Parameters params;
    Membrane membrane;
    SimulVolume simulVolume;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<MolTemplate> molTemplateList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<CreateDestructRxn> createDestructRxns;
    std::map<std::string, int> observablesList;
    copyCounters counterArrays;
};

/*! \brief Make sure the global GSL RNG used by rand_gsl() is allocated.
 *
 * gtest_main.cpp only defines `gsl_rng* r = nullptr;`, and the function under
 * test always evaluates rand_gsl(), so the generator has to exist before the
 * first call or we would dereference a null pointer.
 */
void cfusc_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] global GSL rng was null -> calling srand_gsl(12345)\n";
         const gsl_rng_type *T;
         T = gsl_rng_default;
         r = gsl_rng_alloc(T);
         gsl_rng_set(r, 42);
    }
}

/*! \brief Build a RateState with the requested rate and no ancillary interfaces.
 *
 * `otherIfaceLists` gets exactly one *empty* inner list, which is how the
 * parser represents "this reaction requires no ancillary interfaces". That
 * keeps hasIntangibles() happy without giving it an empty outer vector.
 */
RxnBase::RateState cfusc_make_rate_state(double rate)
{
    RxnBase::RateState rateState;
    rateState.rate = rate;
    rateState.prob = 0.0;
    rateState.otherIfaceLists.emplace_back(); // one empty list == no requirements
    return rateState;
}

/*! \brief Build the complete two-state test system.
 *
 * \param[in] forwardRate rate of the a~A -> a~B reaction (1/s)
 * \param[in] backRate    rate of the a~B -> a~A reaction (1/s)
 *
 * The returned system contains one molecule (state 'A'), one complex holding
 * it, one MolTemplate with a two-state interface, one forward and one back
 * state-change reaction, and a copyCounters whose copyNumSpecies is {1, 0}.
 */
CfuscSystem cfusc_build_system(double forwardRate, double backRate = 0.0)
{
    CfuscSystem sys;

    // ---- Parameters: 1 us time step, state changes enabled ------------------
    sys.params.timeStep = 1.0;
    sys.params.hasUniMolStateChange = true;

    // ---- Membrane: plain box, no implicit lipid by default -------------------
    sys.membrane.isSphere = false;
    sys.membrane.waterBox.x = 100.0;
    sys.membrane.waterBox.y = 100.0;
    sys.membrane.waterBox.z = 100.0;
    sys.membrane.nStates = 0;
    sys.membrane.numberOfFreeLipidsEachState = { 0, 0 };

    // ---- The two interface states ------------------------------------------
    Interface::State stateA("a~A", 'A', kCfuscStateAAbs);
    Interface::State stateB("a~B", 'B', kCfuscStateBAbs);
    // Both states know about the (forward 0, back 0) state-change reaction pair.
    stateA.stateChangeRxns.push_back(std::make_pair(0, 0));
    stateB.stateChangeRxns.push_back(std::make_pair(0, 0));
    stateA.myForwardRxns.push_back(0u);
    stateB.myForwardRxns.push_back(0u);

    // ---- MolTemplate with the single, two-state interface -------------------
    Interface iface("a", std::vector<Interface::State> { stateA, stateB }, Coord(0.0, 0.0, 0.0));
    iface.index = 0;

    MolTemplate molTemplate;
    molTemplate.molName = "A";
    molTemplate.molTypeIndex = 0;
    molTemplate.copies = 1;
    molTemplate.interfaceList.push_back(iface);
    molTemplate.ifacesWithStates = { 0 }; // relative index 0 owns states
    sys.molTemplateList.push_back(molTemplate);

    // ---- The molecule, currently in state 'A' -------------------------------
    Molecule mol;
    mol.index = 0;
    mol.molTypeIndex = 0;
    mol.myComIndex = 0;
    mol.comCoord = Coord(0.0, 0.0, 0.0);
    mol.trajStatus = TrajStatus::none;
    mol.isEmpty = false;
    mol.isGhosted = false;
    mol.isImplicitLipid = false;

    Molecule::Iface molIface;
    molIface.coord = Coord(0.0, 0.0, 0.0);
    molIface.stateIden = 'A';
    molIface.stateIndex = 0;             // relative index within stateList
    molIface.index = kCfuscStateAAbs;    // absolute state index
    molIface.relIndex = 0;               // index in Molecule::interfaceList
    molIface.molTypeIndex = 0;
    molIface.isBound = false;
    mol.interfaceList.push_back(molIface);
    sys.moleculeList.push_back(mol);

    // ---- The complex owning that molecule -----------------------------------
    Complex com;
    com.index = 0;
    com.comCoord = Coord(0.0, 0.0, 0.0);
    com.memberList.push_back(0);
    com.trajStatus = TrajStatus::none;
    sys.complexList.push_back(com);

    // ---- Reactions ----------------------------------------------------------
    // The reactant/product lists are duplicated (index 0 and 1) exactly the way
    // the parser fills them for single-reactant reactions, so any helper that
    // peeks at element 1 stays inside the vector.
    RxnIface reactantA("a", 0, kCfuscStateAAbs, 0, 'A', false);
    RxnIface productB("a", 0, kCfuscStateBAbs, 0, 'B', false);

    ForwardRxn fwd;
    fwd.rxnType = ReactionType::uniMolStateChange;
    fwd.hasStateChange = true;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.rxnLabel = "AtoB";
    fwd.reactantListNew = { reactantA, reactantA };
    fwd.productListNew = { productB, productB };
    fwd.intReactantList = { kCfuscStateAAbs };
    fwd.intProductList = { kCfuscStateBAbs };
    fwd.stateChangeIface = std::make_pair(reactantA, productB);
    fwd.rateList.push_back(cfusc_make_rate_state(forwardRate));
    sys.forwardRxns.push_back(fwd);

    BackRxn back;
    back.rxnType = ReactionType::uniMolStateChange;
    back.hasStateChange = true;
    back.conjForwardRxnIndex = 0;
    back.absRxnIndex = 1;
    back.relRxnIndex = 0;
    back.rxnLabel = "BtoA";
    back.reactantListNew = { productB, productB };
    back.productListNew = { reactantA, reactantA };
    back.intReactantList = { kCfuscStateBAbs };
    back.intProductList = { kCfuscStateAAbs };
    back.stateChangeIface = std::make_pair(productB, reactantA);
    back.rateList.push_back(cfusc_make_rate_state(backRate));
    sys.backRxns.push_back(back);

    // ---- Species counters: one copy of a~A, zero copies of a~B --------------
    sys.counterArrays.copyNumSpecies = { 1, 0 };

    return sys;
}

/*! \brief Thin wrapper that forwards a CfuscSystem to the function under test. */
void cfusc_run(CfuscSystem& sys, unsigned simItr = 0)
{
    cfusc_ensure_rng();
    check_for_unimolstatechange_reactions(simItr, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.membrane);
}

/*! \brief Dump the interesting parts of the system state to stderr. */
void cfusc_report(const CfuscSystem& sys, const char* tag)
{
    if (sys.moleculeList.empty()) {
        std::cerr << "  [" << tag << "] moleculeList is empty\n";
    } else {
        const auto& mi = sys.moleculeList[0].interfaceList[0];
        std::cerr << "  [" << tag << "] iface stateIden='" << mi.stateIden
                  << "' absIndex=" << mi.index << " relStateIndex=" << mi.stateIndex
                  << " relIndex=" << mi.relIndex
                  << " molTrajStatus=" << static_cast<int>(sys.moleculeList[0].trajStatus) << '\n';
    }
    std::cerr << "  [" << tag << "] copyNumSpecies = {" << sys.counterArrays.copyNumSpecies[0]
              << ", " << sys.counterArrays.copyNumSpecies[1] << "}\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a huge forward rate must flip the interface state and update all the
//         associated bookkeeping.
// -----------------------------------------------------------------------------
void test_cfusc_state_change_fires()
{
    std::cerr << "\n[TEST] test_cfusc_state_change_fires\n"
              << "  Source file:   check_for_unimolstatechange_reactions.cpp\n"
              << "  Scenario:      one free molecule in state 'A', forward rate 1e12 /s\n"
              << "                 => prob == 1.0, so the reaction must fire.\n"
              << "  Pass criteria: iface state becomes 'B' (abs index 1, rel index 1),\n"
              << "                 copyNumSpecies moves 1 -> 0 and 0 -> 1, and both the\n"
              << "                 molecule and its complex are marked 'propagated'.\n";

    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    cfusc_report(sys, "before");
    cfusc_run(sys);
    cfusc_report(sys, "after ");

    // The state itself must have flipped via Iface::change_state().
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'B')
        << "interface state identity should have changed from 'A' to 'B'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kCfuscStateBAbs)
        << "absolute interface-state index should be that of the product a~B";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, 1)
        << "relative state index should point at stateList[1] (state 'B')";

    // Species bookkeeping: one copy leaves a~A and one arrives at a~B.
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 0)
        << "copy number of the reactant state a~A should be decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 1)
        << "copy number of the product state a~B should be incremented";

    // Both molecule and parent complex are flagged so nothing moves them again.
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "molecule trajStatus should be set to propagated after a state change";
    EXPECT_EQ(static_cast<int>(sys.complexList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "parent complex trajStatus should be set to propagated after a state change";
}

// -----------------------------------------------------------------------------
// Test 2: a zero rate gives prob == 0, so nothing may ever change.
// -----------------------------------------------------------------------------
void test_cfusc_zero_rate_does_nothing()
{
    std::cerr << "\n[TEST] test_cfusc_zero_rate_does_nothing\n"
              << "  Scenario:      identical system but the forward rate is 0 /s\n"
              << "                 => prob == 0, which can never exceed rand_gsl().\n"
              << "  Pass criteria: state, copy numbers and trajStatus are untouched.\n";

    CfuscSystem sys = cfusc_build_system(0.0);
    cfusc_run(sys);
    cfusc_report(sys, "after ");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'A')
        << "state must stay 'A' when the reaction probability is zero";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kCfuscStateAAbs)
        << "absolute state index must stay unchanged";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 1)
        << "copy number of a~A must stay at 1";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 0)
        << "copy number of a~B must stay at 0";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must stay 'none' when no reaction happens";
}

// -----------------------------------------------------------------------------
// Test 3: the global params.hasUniMolStateChange switch disables the whole body.
// -----------------------------------------------------------------------------
void test_cfusc_flag_disables_everything()
{
    std::cerr << "\n[TEST] test_cfusc_flag_disables_everything\n"
              << "  Scenario:      huge rate, but params.hasUniMolStateChange == false.\n"
              << "  Pass criteria: the molecule is left completely alone.\n";

    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    sys.params.hasUniMolStateChange = false; // master switch off
    cfusc_run(sys);
    cfusc_report(sys, "after ");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'A')
        << "no state change is allowed when hasUniMolStateChange is false";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 1)
        << "copy numbers must be untouched when the feature is disabled";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 0)
        << "copy numbers must be untouched when the feature is disabled";
}

// -----------------------------------------------------------------------------
// Test 4: destroyed (isEmpty) and ghosted (isGhosted) molecules are skipped.
// -----------------------------------------------------------------------------
void test_cfusc_skips_empty_and_ghosted()
{
    std::cerr << "\n[TEST] test_cfusc_skips_empty_and_ghosted\n"
              << "  Scenario:      huge rate, but the molecule is flagged isEmpty\n"
              << "                 (first pass) and isGhosted (second pass).\n"
              << "  Pass criteria: both flags cause an early 'continue', so nothing\n"
              << "                 about the molecule or the counters changes.\n";

    // --- isEmpty -------------------------------------------------------------
    CfuscSystem emptySys = cfusc_build_system(kCfuscHugeRate);
    emptySys.moleculeList[0].isEmpty = true;
    cfusc_run(emptySys);
    cfusc_report(emptySys, "isEmpty");
    EXPECT_EQ(emptySys.moleculeList[0].interfaceList[0].stateIden, 'A')
        << "an isEmpty molecule must never react";
    EXPECT_EQ(emptySys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 0)
        << "an isEmpty molecule must not touch the species counters";

    // --- isGhosted -----------------------------------------------------------
    CfuscSystem ghostSys = cfusc_build_system(kCfuscHugeRate);
    ghostSys.moleculeList[0].isGhosted = true;
    cfusc_run(ghostSys);
    cfusc_report(ghostSys, "ghosted");
    EXPECT_EQ(ghostSys.moleculeList[0].interfaceList[0].stateIden, 'A')
        << "a ghosted (MPI copy) molecule must never react";
    EXPECT_EQ(ghostSys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 0)
        << "a ghosted molecule must not touch the species counters";
}

// -----------------------------------------------------------------------------
// Test 5: only molecules whose trajStatus is 'none' may change state.
// -----------------------------------------------------------------------------
void test_cfusc_requires_trajstatus_none()
{
    std::cerr << "\n[TEST] test_cfusc_requires_trajstatus_none\n"
              << "  Scenario:      huge rate, but the molecule was already propagated\n"
              << "                 earlier in the time step.\n"
              << "  Pass criteria: the explicit-molecule branch is skipped entirely.\n";

    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    sys.moleculeList[0].trajStatus = TrajStatus::propagated;
    cfusc_run(sys);
    cfusc_report(sys, "after ");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'A')
        << "an already propagated molecule must not change state";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 1)
        << "counters must be unchanged for an already propagated molecule";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 0)
        << "counters must be unchanged for an already propagated molecule";
}

// -----------------------------------------------------------------------------
// Test 6: a MolTemplate with no state-bearing interfaces short-circuits.
// -----------------------------------------------------------------------------
void test_cfusc_no_state_ifaces()
{
    std::cerr << "\n[TEST] test_cfusc_no_state_ifaces\n"
              << "  Scenario:      huge rate, but MolTemplate::ifacesWithStates is empty.\n"
              << "  Pass criteria: nothing is examined, nothing changes.\n";

    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    sys.molTemplateList[0].ifacesWithStates.clear();
    cfusc_run(sys);
    cfusc_report(sys, "after ");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'A')
        << "with no state-bearing interfaces there is nothing to change";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 0)
        << "counters must not be touched when ifacesWithStates is empty";
}

// -----------------------------------------------------------------------------
// Test 7: bound interfaces fall into the (currently unimplemented) TODO branch.
// -----------------------------------------------------------------------------
void test_cfusc_bound_interface_is_ignored()
{
    std::cerr << "\n[TEST] test_cfusc_bound_interface_is_ignored\n"
              << "  Scenario:      huge rate, but the state-bearing interface is bound.\n"
              << "  Pass criteria: the bound branch is a TODO, so no state change and\n"
              << "                 no counter update may occur.\n";

    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    sys.moleculeList[0].interfaceList[0].isBound = true;
    sys.moleculeList[0].interfaceList[0].interaction.partnerIndex = 1;
    sys.moleculeList[0].interfaceList[0].interaction.partnerIfaceIndex = 0;
    cfusc_run(sys);
    cfusc_report(sys, "after ");

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'A')
        << "bound interfaces are not handled yet and must not change state";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 1)
        << "counters must be unchanged for a bound interface";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 0)
        << "counters must be unchanged for a bound interface";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must remain 'none' for a bound interface";
}

// -----------------------------------------------------------------------------
// Test 8: when the reaction is observed the matching observable is incremented,
//         and an unknown observable label is tolerated.
// -----------------------------------------------------------------------------
void test_cfusc_observable_bookkeeping()
{
    std::cerr << "\n[TEST] test_cfusc_observable_bookkeeping\n"
              << "  Scenario:      firing reaction flagged isObserved with label 'AtoB'.\n"
              << "  Pass criteria: observablesList[\"AtoB\"] goes 0 -> 1; an unknown\n"
              << "                 label leaves the map size alone (no insert, no crash).\n";

    // --- known label ---------------------------------------------------------
    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    sys.forwardRxns[0].isObserved = true;
    sys.forwardRxns[0].observeLabel = "AtoB";
    sys.observablesList["AtoB"] = 0;
    cfusc_run(sys);
    std::cerr << "  [after ] observablesList[\"AtoB\"] = " << sys.observablesList["AtoB"] << '\n';
    EXPECT_EQ(sys.observablesList["AtoB"], 1)
        << "the observable of the fired reaction should have been incremented once";

    // --- unknown label -------------------------------------------------------
    CfuscSystem sys2 = cfusc_build_system(kCfuscHugeRate);
    sys2.forwardRxns[0].isObserved = true;
    sys2.forwardRxns[0].observeLabel = "notDefined";
    cfusc_run(sys2);
    std::cerr << "  [after ] observablesList.size() = " << sys2.observablesList.size() << '\n';
    EXPECT_EQ(sys2.observablesList.size(), 0u)
        << "an undefined observable label must not be inserted into the map";
    // The state change itself must still have happened.
    EXPECT_EQ(sys2.moleculeList[0].interfaceList[0].stateIden, 'B')
        << "an undefined observable label must not prevent the state change";
}

// -----------------------------------------------------------------------------
// Test 9: implicit-lipid branch with a huge rate: species counters *and* the
//         per-state free-lipid counters are updated, but the interface state of
//         the implicit-lipid pseudo-molecule itself is left alone.
// -----------------------------------------------------------------------------
void test_cfusc_implicit_lipid_state_change()
{
    std::cerr << "\n[TEST] test_cfusc_implicit_lipid_state_change\n"
              << "  Scenario:      implicit-lipid molecule, nStates == 1, one free lipid\n"
              << "                 in state 'A', huge forward rate.\n"
              << "  Pass criteria: copyNumSpecies moves 1 -> 0 / 0 -> 1 and\n"
              << "                 numberOfFreeLipidsEachState moves 1 -> 0 / 0 -> 1,\n"
              << "                 while the pseudo-molecule's own iface is unchanged\n"
              << "                 (the IL branch never calls change_state()).\n";

    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    sys.moleculeList[0].isImplicitLipid = true;
    sys.molTemplateList[0].isImplicitLipid = true;
    sys.membrane.implicitLipid = true;
    sys.membrane.nStates = 1;                                  // only state 'A' is looped over
    sys.membrane.numberOfFreeLipidsEachState = { 1, 0 };       // one free a~A lipid
    sys.membrane.implicitlipidIndex = 0;

    std::cerr << "  [before] freeLipids = {" << sys.membrane.numberOfFreeLipidsEachState[0] << ", "
              << sys.membrane.numberOfFreeLipidsEachState[1] << "}\n";
    cfusc_run(sys);
    cfusc_report(sys, "after ");
    std::cerr << "  [after ] freeLipids = {" << sys.membrane.numberOfFreeLipidsEachState[0] << ", "
              << sys.membrane.numberOfFreeLipidsEachState[1] << "}\n";

    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 0)
        << "implicit-lipid state change should decrement the reactant species count";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 1)
        << "implicit-lipid state change should increment the product species count";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], 0)
        << "the free-lipid count of the reactant state should be decremented";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[1], 1)
        << "the free-lipid count of the product state should be incremented";

    // The implicit-lipid branch only bookkeeps; it must not rewrite the iface.
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'A')
        << "the implicit-lipid pseudo-molecule's own state identity is not rewritten";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kCfuscStateAAbs)
        << "the temporary index shift must be undone, leaving the original abs index";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].relIndex, 0)
        << "the temporary relIndex shift must be undone";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, 0)
        << "the temporary stateIndex shift must be undone";
}

// -----------------------------------------------------------------------------
// Test 10: implicit-lipid branch with zero rates. Nothing may change, and the
//          temporary index shifting done per lipid state has to be reverted.
// -----------------------------------------------------------------------------
void test_cfusc_implicit_lipid_zero_rate_restores_indices()
{
    std::cerr << "\n[TEST] test_cfusc_implicit_lipid_zero_rate_restores_indices\n"
              << "  Scenario:      implicit-lipid molecule, nStates == 2 with 2 and 1 free\n"
              << "                 lipids, both reaction rates 0 /s.\n"
              << "  Pass criteria: counters and free-lipid counts are untouched and the\n"
              << "                 interface's index/relIndex/stateIndex are restored\n"
              << "                 after the temporary per-state shifting.\n";

    CfuscSystem sys = cfusc_build_system(0.0, 0.0);
    sys.moleculeList[0].isImplicitLipid = true;
    sys.molTemplateList[0].isImplicitLipid = true;
    sys.membrane.implicitLipid = true;
    sys.membrane.nStates = 2;
    sys.membrane.numberOfFreeLipidsEachState = { 2, 1 };
    sys.membrane.implicitlipidIndex = 0;
    sys.counterArrays.copyNumSpecies = { 2, 1 };

    cfusc_run(sys);
    cfusc_report(sys, "after ");
    std::cerr << "  [after ] freeLipids = {" << sys.membrane.numberOfFreeLipidsEachState[0] << ", "
              << sys.membrane.numberOfFreeLipidsEachState[1] << "}\n";

    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 2)
        << "zero-rate implicit-lipid reactions must not change species counts";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 1)
        << "zero-rate implicit-lipid reactions must not change species counts";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], 2)
        << "zero-rate implicit-lipid reactions must not change free-lipid counts";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[1], 1)
        << "zero-rate implicit-lipid reactions must not change free-lipid counts";

    // The routine shifts these three fields by the lipid state index and then
    // subtracts the same amount again; verify the round trip is exact.
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kCfuscStateAAbs)
        << "absolute interface index must be restored after the loop";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].relIndex, 0)
        << "relIndex must be restored after the loop";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, 0)
        << "stateIndex must be restored after the loop";
}

// -----------------------------------------------------------------------------
// Test 11: an empty molecule list is handled gracefully (no iteration at all).
// -----------------------------------------------------------------------------
void test_cfusc_empty_molecule_list()
{
    std::cerr << "\n[TEST] test_cfusc_empty_molecule_list\n"
              << "  Scenario:      no molecules at all in the system.\n"
              << "  Pass criteria: the call returns without touching the counters.\n";

    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    sys.moleculeList.clear();
    sys.complexList.clear();
    cfusc_run(sys);
    cfusc_report(sys, "after ");

    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 1)
        << "an empty molecule list must leave the counters untouched";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 0)
        << "an empty molecule list must leave the counters untouched";
}

// -----------------------------------------------------------------------------
// Test 12: repeated calls after a successful state change. Because the molecule
//          is now 'propagated', a second call in the same step must be a no-op.
// -----------------------------------------------------------------------------
void test_cfusc_second_call_is_noop()
{
    std::cerr << "\n[TEST] test_cfusc_second_call_is_noop\n"
              << "  Scenario:      call the routine twice in a row with a huge rate.\n"
              << "  Pass criteria: the first call flips 'A'->'B' and marks the molecule\n"
              << "                 propagated; the second call therefore does nothing,\n"
              << "                 so copy numbers stay at {0, 1}.\n";

    CfuscSystem sys = cfusc_build_system(kCfuscHugeRate);
    cfusc_run(sys, 0);
    cfusc_report(sys, "call#1");
    ASSERT_FALSE(sys.moleculeList.empty()); // guard the indexing below

    // Sanity: the first call really did fire.
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'B')
        << "first call should have performed the state change";

    cfusc_run(sys, 1);
    cfusc_report(sys, "call#2");

    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateAAbs], 0)
        << "second call must not double count the reactant species";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCfuscStateBAbs], 1)
        << "second call must not double count the product species";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'B')
        << "state should still be 'B' after the (no-op) second call";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper runs inside its own TEST so
// every scenario is reported separately and a failure in one does not stop the
// others from running.
// -----------------------------------------------------------------------------
TEST(CheckForUnimolStateChangeReactions, StateChangeFires) { test_cfusc_state_change_fires(); }
TEST(CheckForUnimolStateChangeReactions, ZeroRateDoesNothing) { test_cfusc_zero_rate_does_nothing(); }
TEST(CheckForUnimolStateChangeReactions, FlagDisablesEverything) { test_cfusc_flag_disables_everything(); }
TEST(CheckForUnimolStateChangeReactions, SkipsEmptyAndGhosted) { test_cfusc_skips_empty_and_ghosted(); }
TEST(CheckForUnimolStateChangeReactions, RequiresTrajStatusNone) { test_cfusc_requires_trajstatus_none(); }
TEST(CheckForUnimolStateChangeReactions, NoStateIfaces) { test_cfusc_no_state_ifaces(); }
TEST(CheckForUnimolStateChangeReactions, BoundInterfaceIsIgnored) { test_cfusc_bound_interface_is_ignored(); }
TEST(CheckForUnimolStateChangeReactions, ObservableBookkeeping) { test_cfusc_observable_bookkeeping(); }
TEST(CheckForUnimolStateChangeReactions, ImplicitLipidStateChange) { test_cfusc_implicit_lipid_state_change(); }
TEST(CheckForUnimolStateChangeReactions, ImplicitLipidZeroRateRestoresIndices)
{
    test_cfusc_implicit_lipid_zero_rate_restores_indices();
}
TEST(CheckForUnimolStateChangeReactions, EmptyMoleculeList) { test_cfusc_empty_molecule_list(); }
TEST(CheckForUnimolStateChangeReactions, SecondCallIsNoop) { test_cfusc_second_call_is_noop(); }