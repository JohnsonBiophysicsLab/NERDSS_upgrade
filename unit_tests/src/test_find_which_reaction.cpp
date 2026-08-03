/*! \file test_find_which_reaction.cpp
 *
 * ### Unit test for src/reactions/find_which_reaction.cpp
 *
 * Function under test:
 *
 *     void find_which_reaction(int ifaceIndex1, int ifaceIndex2, int& rxnIndex,
 *                              int& rateIndex, bool& isStateChangeBackRxn,
 *                              const Interface::State& currState,
 *                              const Molecule& reactMol1, const Molecule& reactMol2,
 *                              const std::vector<ForwardRxn>& forwardRxns,
 *                              const std::vector<BackRxn>& backRxns,
 *                              const std::vector<MolTemplate>& molTemplateList)
 *
 * The routine walks the list of ForwardRxns attached to the state of the first
 * reacting interface (`currState.myForwardRxns`) and tries to find a reaction
 * whose `reactantListNew` contains BOTH of the reacting interfaces.  When such
 * a reaction is found, it selects a rate state (via hasIntangibles()) and
 * returns the relative reaction index and the rate index through the reference
 * output parameters.  Special handling exists for
 *
 *   - implicit lipid partners (reactant index of an implicit lipid is forced to 1)
 *   - bimolecular state-change reactions, where the *products* are searched and
 *     the conjugate BackRxn rate list is used (isStateChangeBackRxn == true)
 *
 * The tests below build minimal, hand-made reaction/molecule structures so that
 * each of those code paths can be exercised in isolation.  All reaction rate
 * states are created with two EMPTY ancillary ("other") interface lists, which
 * means no ancillary interfaces are required for the reaction to occur, so the
 * rate state always matches.
 */

#include "reactions/shared_reaction_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// `totMatches` is a global counter defined in unit_tests/src/gtest_main.cpp and
// declared extern in classes/class_Rxns.hpp; find_which_reaction() increments it
// every time it resolves a bimolecular reaction.

namespace {

// -----------------------------------------------------------------------------
// Helper: build a RxnIface (one reactant/product entry of a reaction).
// -----------------------------------------------------------------------------
RxnIface fwr_make_rxn_iface(int molTypeIndex, int absIfaceIndex, int relIfaceIndex)
{
    RxnIface oneIface;
    oneIface.ifaceName = "iface" + std::to_string(absIfaceIndex);
    oneIface.molTypeIndex = molTypeIndex;
    oneIface.absIfaceIndex = absIfaceIndex;
    oneIface.relIfaceIndex = relIfaceIndex;
    oneIface.requiresState = '\0'; // no particular state required
    oneIface.requiresInteraction = false; // does not need to be bound
    return oneIface;
}

// -----------------------------------------------------------------------------
// Helper: build a RateState with two EMPTY ancillary interface lists.
//
// find_which_reaction() always dereferences otherIfaceLists[0] and [1] when
// scoring a match, so both lists must exist.  Leaving them empty means the
// reaction has no ancillary requirements and therefore always applies.
// -----------------------------------------------------------------------------
RxnBase::RateState fwr_make_rate_state(double rate)
{
    RxnBase::RateState oneRate;
    oneRate.rate = rate;
    oneRate.otherIfaceLists.resize(2); // one (empty) list per reactant
    return oneRate;
}

// -----------------------------------------------------------------------------
// Helper: build a molecule with a set of interfaces whose absolute (state)
// indices are given by `ifaceAbsIndices`.
// -----------------------------------------------------------------------------
Molecule fwr_make_molecule(int molTypeIndex, const std::vector<int>& ifaceAbsIndices,
    bool isImplicitLipid = false)
{
    Molecule mol;
    mol.molTypeIndex = molTypeIndex;
    mol.isImplicitLipid = isImplicitLipid;
    mol.index = 0;
    mol.myComIndex = 0;

    for (std::size_t i = 0; i < ifaceAbsIndices.size(); ++i) {
        Molecule::Iface oneIface;
        oneIface.index = ifaceAbsIndices[i]; // absolute (state) index
        oneIface.relIndex = static_cast<int>(i);
        oneIface.stateIndex = 0;
        oneIface.molTypeIndex = molTypeIndex;
        oneIface.isBound = false;
        mol.interfaceList.push_back(oneIface);
    }
    return mol;
}

// -----------------------------------------------------------------------------
// Helper: two molecule templates, neither of which is an implicit lipid.
// The reaction reactant entries reference these by molTypeIndex, so the list
// must be large enough for every molTypeIndex used in a test.
// -----------------------------------------------------------------------------
std::vector<MolTemplate> fwr_make_templates()
{
    std::vector<MolTemplate> molTemplateList(2);
    molTemplateList[0].molTypeIndex = 0;
    molTemplateList[0].molName = "A";
    molTemplateList[0].isImplicitLipid = false;
    molTemplateList[1].molTypeIndex = 1;
    molTemplateList[1].molName = "B";
    molTemplateList[1].isImplicitLipid = false;
    return molTemplateList;
}

// -----------------------------------------------------------------------------
// Helper: a plain bimolecular ForwardRxn between two absolute interface indices.
// -----------------------------------------------------------------------------
ForwardRxn fwr_make_bimolecular_rxn(int relRxnIndex, int absIface1, int absIface2, int numRates)
{
    ForwardRxn oneRxn;
    oneRxn.rxnType = ReactionType::bimolecular;
    oneRxn.relRxnIndex = relRxnIndex;
    oneRxn.absRxnIndex = relRxnIndex;
    oneRxn.conjBackRxnIndex = -1; // no conjugate back reaction (default)
    oneRxn.reactantListNew.push_back(fwr_make_rxn_iface(0, absIface1, 0));
    oneRxn.reactantListNew.push_back(fwr_make_rxn_iface(1, absIface2, 0));
    oneRxn.intReactantList = { absIface1, absIface2 };
    for (int i = 0; i < numRates; ++i)
        oneRxn.rateList.push_back(fwr_make_rate_state(1.0 + i));
    return oneRxn;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: The interface state lists no forward reactions at all.
//         => the output parameters must be left completely untouched.
// -----------------------------------------------------------------------------
void test_fwr_no_forward_reactions()
{
    std::cerr << "\n[TEST] test_fwr_no_forward_reactions\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      currState.myForwardRxns is empty, so the main\n"
              << "                 loop body never executes.\n"
              << "  Pass criteria: rxnIndex/rateIndex/isStateChangeBackRxn keep\n"
              << "                 their sentinel input values.\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };
    std::vector<ForwardRxn> forwardRxns; // deliberately empty
    std::vector<BackRxn> backRxns; // deliberately empty

    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 1 }) };

    Interface::State currState; // myForwardRxns left empty

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    const unsigned long matchesBefore { totMatches };
    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex
              << " isStateChangeBackRxn=" << std::boolalpha << isStateChangeBackRxn << '\n';

    EXPECT_EQ(rxnIndex, -1) << "rxnIndex must not be assigned when no reactions are listed";
    EXPECT_EQ(rateIndex, -1) << "rateIndex must not be assigned when no reactions are listed";
    EXPECT_FALSE(isStateChangeBackRxn) << "isStateChangeBackRxn must remain false";
    EXPECT_EQ(totMatches, matchesBefore) << "totMatches must not be incremented";
}

// -----------------------------------------------------------------------------
// Test 2: A forward reaction exists but neither of the reacting interfaces
//         appears in its reactantListNew.
//         => nothing is assigned (reactIndex1 and reactIndex2 remain -1).
// -----------------------------------------------------------------------------
void test_fwr_no_matching_reactants()
{
    std::cerr << "\n[TEST] test_fwr_no_matching_reactants\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      the single listed reaction wants absolute\n"
              << "                 interfaces 100 and 101, the molecules carry 0 and 1.\n"
              << "  Pass criteria: no reaction is selected (outputs untouched).\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };
    std::vector<ForwardRxn> forwardRxns { fwr_make_bimolecular_rxn(3, 100, 101, 1) };
    std::vector<BackRxn> backRxns;

    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 1 }) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0); // points at forwardRxns[0]

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    const unsigned long matchesBefore { totMatches };
    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex << '\n';

    EXPECT_EQ(rxnIndex, -1) << "rxnIndex must stay -1 when the reactants do not match";
    EXPECT_EQ(rateIndex, -1) << "rateIndex must stay -1 when the reactants do not match";
    EXPECT_FALSE(isStateChangeBackRxn) << "isStateChangeBackRxn must remain false";
    EXPECT_EQ(totMatches, matchesBefore) << "totMatches must not be incremented";
}

// -----------------------------------------------------------------------------
// Test 3: Only the FIRST reacting interface is found in the reactant list.
//         => both reactant indices are required, so nothing is assigned.
// -----------------------------------------------------------------------------
void test_fwr_only_one_reactant_matches()
{
    std::cerr << "\n[TEST] test_fwr_only_one_reactant_matches\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      reaction reactants are (0, 101); molecule 1 carries\n"
              << "                 interface 0 (match) and molecule 2 carries 1 (no match).\n"
              << "  Pass criteria: partial matches are rejected, outputs untouched.\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };
    std::vector<ForwardRxn> forwardRxns { fwr_make_bimolecular_rxn(4, 0, 101, 1) };
    std::vector<BackRxn> backRxns;

    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 1 }) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0);

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex << '\n';

    EXPECT_EQ(rxnIndex, -1) << "a reaction with only one matching reactant must be skipped";
    EXPECT_EQ(rateIndex, -1) << "rateIndex must stay -1 for a partial match";
    EXPECT_FALSE(isStateChangeBackRxn) << "isStateChangeBackRxn must remain false";
}

// -----------------------------------------------------------------------------
// Test 4: A well formed bimolecular reaction with exactly one rate state that
//         has no ancillary interface requirements.
//         => the relative reaction index and rate index 0 are returned and the
//            global totMatches counter is incremented by one.
// -----------------------------------------------------------------------------
void test_fwr_single_bimolecular_match()
{
    std::cerr << "\n[TEST] test_fwr_single_bimolecular_match\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      both reacting interfaces (abs 0 and abs 1) are found\n"
              << "                 in reactantListNew; the reaction has one rate state\n"
              << "                 with no required ancillary interfaces.\n"
              << "  Pass criteria: rxnIndex == relRxnIndex (7), rateIndex == 0,\n"
              << "                 isStateChangeBackRxn == false, totMatches += 1.\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };
    std::vector<ForwardRxn> forwardRxns { fwr_make_bimolecular_rxn(7, 0, 1, 1) };
    std::vector<BackRxn> backRxns;

    // Molecule 1 reacts through its interface 0 (absolute index 0),
    // molecule 2 reacts through its interface 0 (absolute index 1).
    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 1 }) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0);

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    const unsigned long matchesBefore { totMatches };
    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex
              << " totMatches delta=" << (totMatches - matchesBefore) << '\n';

    EXPECT_EQ(rxnIndex, 7) << "the relative index of the matching reaction should be returned";
    EXPECT_EQ(rateIndex, 0) << "the only rate state (index 0) should be selected";
    EXPECT_FALSE(isStateChangeBackRxn) << "a plain bimolecular reaction is not a state change";
    EXPECT_EQ(totMatches - matchesBefore, 1ul) << "totMatches should be incremented exactly once";
}

// -----------------------------------------------------------------------------
// Test 5: The same reaction, but with TWO equally valid rate states.
//         When several rate states match, the code picks the one requiring the
//         most ancillary interfaces; with all counts equal to zero the first
//         entry in the match list (rate index 0) wins.
// -----------------------------------------------------------------------------
void test_fwr_multiple_rate_states()
{
    std::cerr << "\n[TEST] test_fwr_multiple_rate_states\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      reaction has two rate states, both with zero\n"
              << "                 ancillary interfaces, so both match.\n"
              << "  Pass criteria: the 'best fit' fallback selects rate index 0 and\n"
              << "                 rxnIndex is still the matching reaction (11).\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };
    std::vector<ForwardRxn> forwardRxns { fwr_make_bimolecular_rxn(11, 0, 1, 2) };
    std::vector<BackRxn> backRxns;

    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 1 }) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0);

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    const unsigned long matchesBefore { totMatches };
    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex
              << " totMatches delta=" << (totMatches - matchesBefore) << '\n';

    EXPECT_EQ(rxnIndex, 11) << "the matching reaction's relative index should be returned";
    EXPECT_EQ(rateIndex, 0) << "with equal ancillary counts the first matching rate is used";
    EXPECT_FALSE(isStateChangeBackRxn) << "isStateChangeBackRxn must remain false";
    EXPECT_EQ(totMatches - matchesBefore, 1ul) << "one resolved reaction => one match counted";
}

// -----------------------------------------------------------------------------
// Test 6: Several forward reactions listed on the state; only the second one
//         actually contains both reacting interfaces.
//         => the loop must keep scanning and return the second reaction.
// -----------------------------------------------------------------------------
void test_fwr_scans_all_listed_reactions()
{
    std::cerr << "\n[TEST] test_fwr_scans_all_listed_reactions\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      state lists two reactions; the first does not match\n"
              << "                 the reacting interfaces, the second one does.\n"
              << "  Pass criteria: rxnIndex == relRxnIndex of the SECOND reaction (21).\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };
    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(fwr_make_bimolecular_rxn(20, 50, 51, 1)); // does not match
    forwardRxns.push_back(fwr_make_bimolecular_rxn(21, 0, 1, 1)); // matches

    std::vector<BackRxn> backRxns;

    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 1 }) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0);
    currState.myForwardRxns.push_back(1);

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex << '\n';

    EXPECT_EQ(rxnIndex, 21) << "the second (matching) reaction should be selected";
    EXPECT_EQ(rateIndex, 0) << "its only rate state should be selected";
    EXPECT_FALSE(isStateChangeBackRxn) << "isStateChangeBackRxn must remain false";
}

// -----------------------------------------------------------------------------
// Test 7: The second reactant is an implicit lipid.
//
// The interface index of an implicit lipid is never compared: instead, the
// reaction's reactant entry is identified through
// molTemplateList[...].isImplicitLipid and the reactant index is forced to 1.
// -----------------------------------------------------------------------------
void test_fwr_implicit_lipid_partner()
{
    std::cerr << "\n[TEST] test_fwr_implicit_lipid_partner\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      molecule 2 is an implicit lipid and the second\n"
              << "                 reactant of the reaction belongs to an implicit\n"
              << "                 lipid MolTemplate (its abs iface index differs).\n"
              << "  Pass criteria: reaction is matched anyway; rxnIndex == 31,\n"
              << "                 rateIndex == 0.\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };
    molTemplateList[1].isImplicitLipid = true; // template of the lipid partner

    // Reactant 0 -> molecule 1 interface (abs 0); reactant 1 -> implicit lipid
    // (abs 55, which intentionally does not equal any real interface index).
    std::vector<ForwardRxn> forwardRxns { fwr_make_bimolecular_rxn(31, 0, 55, 1) };
    std::vector<BackRxn> backRxns;

    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 55 }, /*isImplicitLipid=*/true) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0);

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    const unsigned long matchesBefore { totMatches };
    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex
              << " totMatches delta=" << (totMatches - matchesBefore) << '\n';

    EXPECT_EQ(rxnIndex, 31) << "an implicit-lipid reaction should still be resolved";
    EXPECT_EQ(rateIndex, 0) << "the single rate state should be selected";
    EXPECT_FALSE(isStateChangeBackRxn) << "this is not a state change reaction";
    EXPECT_EQ(totMatches - matchesBefore, 1ul) << "one resolved reaction => one match counted";
}

// -----------------------------------------------------------------------------
// Test 8: Bimolecular STATE CHANGE reaction resolved through the conjugate
//         back reaction.
//
// Setup: the forward reaction's reactantListNew does NOT contain the reacting
// interfaces (so reactIndex1/reactIndex2 stay -1), but its productListNew does.
// Because rxnType == biMolStateChange and conjBackRxnIndex > 0, the code scans
// the products, then evaluates the rate states of backRxns[conjBackRxnIndex].
// A single matching back rate results in isStateChangeBackRxn == true and
// rxnIndex == relRxnIndex (rateIndex is deliberately left untouched by the
// implementation on this path).
// -----------------------------------------------------------------------------
void test_fwr_bimolecular_state_change_back_rxn()
{
    std::cerr << "\n[TEST] test_fwr_bimolecular_state_change_back_rxn\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      biMolStateChange reaction with conjBackRxnIndex=1;\n"
              << "                 the reacting interfaces appear in productListNew,\n"
              << "                 not in reactantListNew.\n"
              << "  Pass criteria: isStateChangeBackRxn == true and rxnIndex == 41.\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };

    // The forward (state change) reaction: reactants are unrelated absolute
    // indices (100/101) so they cannot be matched, while the products carry the
    // absolute indices actually held by the two molecules (0 and 1).
    ForwardRxn stateChangeRxn;
    stateChangeRxn.rxnType = ReactionType::biMolStateChange;
    stateChangeRxn.relRxnIndex = 41;
    stateChangeRxn.absRxnIndex = 41;
    stateChangeRxn.conjBackRxnIndex = 1; // must be > 0 to enter the branch
    stateChangeRxn.reactantListNew.push_back(fwr_make_rxn_iface(0, 100, 0));
    stateChangeRxn.reactantListNew.push_back(fwr_make_rxn_iface(1, 101, 0));
    stateChangeRxn.productListNew.push_back(fwr_make_rxn_iface(0, 0, 0));
    stateChangeRxn.productListNew.push_back(fwr_make_rxn_iface(1, 1, 0));
    stateChangeRxn.intReactantList = { 100, 101 };
    stateChangeRxn.rateList.push_back(fwr_make_rate_state(1.0)); // needed for
    // the pointer arithmetic inside the state-change branch

    std::vector<ForwardRxn> forwardRxns { stateChangeRxn };

    // backRxns[1] is the conjugate reaction whose rate list is scanned; give it
    // exactly one rate state with no ancillary interface requirements.
    std::vector<BackRxn> backRxns(2);
    backRxns[1].rxnType = ReactionType::biMolStateChange;
    backRxns[1].relRxnIndex = 1;
    backRxns[1].rateList.push_back(fwr_make_rate_state(2.0));

    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 1 }) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0);

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex
              << " isStateChangeBackRxn=" << std::boolalpha << isStateChangeBackRxn << '\n';

    EXPECT_TRUE(isStateChangeBackRxn)
        << "a bimolecular state change matched through its products must set the flag";
    EXPECT_EQ(rxnIndex, 41) << "the relative index of the state change reaction should be returned";
}

// -----------------------------------------------------------------------------
// Test 9: A biMolStateChange reaction whose products also fail to match must be
//         skipped entirely (no output assignment, no crash).
// -----------------------------------------------------------------------------
void test_fwr_state_change_no_product_match()
{
    std::cerr << "\n[TEST] test_fwr_state_change_no_product_match\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      biMolStateChange reaction where neither the\n"
              << "                 reactants nor the products contain the reacting\n"
              << "                 interfaces.\n"
              << "  Pass criteria: outputs untouched, isStateChangeBackRxn == false.\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };

    ForwardRxn stateChangeRxn;
    stateChangeRxn.rxnType = ReactionType::biMolStateChange;
    stateChangeRxn.relRxnIndex = 51;
    stateChangeRxn.conjBackRxnIndex = 1;
    stateChangeRxn.reactantListNew.push_back(fwr_make_rxn_iface(0, 100, 0));
    stateChangeRxn.reactantListNew.push_back(fwr_make_rxn_iface(1, 101, 0));
    stateChangeRxn.productListNew.push_back(fwr_make_rxn_iface(0, 102, 0));
    stateChangeRxn.productListNew.push_back(fwr_make_rxn_iface(1, 103, 0));
    stateChangeRxn.intReactantList = { 100, 101 };
    stateChangeRxn.rateList.push_back(fwr_make_rate_state(1.0));

    std::vector<ForwardRxn> forwardRxns { stateChangeRxn };

    std::vector<BackRxn> backRxns(2);
    backRxns[1].rateList.push_back(fwr_make_rate_state(2.0));

    Molecule mol1 { fwr_make_molecule(0, { 0 }) };
    Molecule mol2 { fwr_make_molecule(1, { 1 }) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0);

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };

    find_which_reaction(0, 0, rxnIndex, rateIndex, isStateChangeBackRxn, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);

    std::cerr << "  Result: rxnIndex=" << rxnIndex << " rateIndex=" << rateIndex
              << " isStateChangeBackRxn=" << std::boolalpha << isStateChangeBackRxn << '\n';

    EXPECT_EQ(rxnIndex, -1) << "an unmatched state change reaction must not be selected";
    EXPECT_EQ(rateIndex, -1) << "rateIndex must stay at its sentinel value";
    EXPECT_FALSE(isStateChangeBackRxn) << "isStateChangeBackRxn must remain false";
}

// -----------------------------------------------------------------------------
// Test 10: Non-zero interface indices.  Verify that ifaceIndex1/ifaceIndex2 are
//          really used to look up the reacting interface inside each molecule's
//          interfaceList (i.e. the second interface of each molecule reacts).
// -----------------------------------------------------------------------------
void test_fwr_uses_supplied_interface_indices()
{
    std::cerr << "\n[TEST] test_fwr_uses_supplied_interface_indices\n"
              << "  Source file:   src/reactions/find_which_reaction.cpp\n"
              << "  Function:      find_which_reaction()\n"
              << "  Scenario:      each molecule has two interfaces; the reaction is\n"
              << "                 defined for the SECOND interface of each molecule\n"
              << "                 (abs 5 and abs 6).\n"
              << "  Pass criteria: passing ifaceIndex1=1, ifaceIndex2=1 resolves the\n"
              << "                 reaction (rxnIndex 61) while ifaceIndex 0,0 does not.\n";

    std::vector<MolTemplate> molTemplateList { fwr_make_templates() };
    std::vector<ForwardRxn> forwardRxns { fwr_make_bimolecular_rxn(61, 5, 6, 1) };
    std::vector<BackRxn> backRxns;

    // interface 0 -> abs 3/4 (not in the reaction), interface 1 -> abs 5/6.
    Molecule mol1 { fwr_make_molecule(0, { 3, 5 }) };
    Molecule mol2 { fwr_make_molecule(1, { 4, 6 }) };

    Interface::State currState;
    currState.myForwardRxns.push_back(0);

    // First call with the WRONG interface indices: nothing must be found.
    int rxnIndexWrong { -1 };
    int rateIndexWrong { -1 };
    bool flagWrong { false };
    find_which_reaction(0, 0, rxnIndexWrong, rateIndexWrong, flagWrong, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);
    std::cerr << "  With ifaceIndex (0,0): rxnIndex=" << rxnIndexWrong << '\n';
    EXPECT_EQ(rxnIndexWrong, -1) << "interfaces 0/0 are not reactants of this reaction";

    // Second call with the correct interface indices: the reaction is found.
    int rxnIndexRight { -1 };
    int rateIndexRight { -1 };
    bool flagRight { false };
    find_which_reaction(1, 1, rxnIndexRight, rateIndexRight, flagRight, currState, mol1, mol2,
        forwardRxns, backRxns, molTemplateList);
    std::cerr << "  With ifaceIndex (1,1): rxnIndex=" << rxnIndexRight
              << " rateIndex=" << rateIndexRight << '\n';
    EXPECT_EQ(rxnIndexRight, 61) << "interfaces 1/1 must resolve to the defined reaction";
    EXPECT_EQ(rateIndexRight, 0) << "the single rate state should be selected";
    EXPECT_FALSE(flagRight) << "isStateChangeBackRxn must remain false";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: one TEST per scenario so a failure in one does not stop
// the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(FindWhichReaction, NoForwardReactions) { test_fwr_no_forward_reactions(); }
TEST(FindWhichReaction, NoMatchingReactants) { test_fwr_no_matching_reactants(); }
TEST(FindWhichReaction, OnlyOneReactantMatches) { test_fwr_only_one_reactant_matches(); }
TEST(FindWhichReaction, SingleBimolecularMatch) { test_fwr_single_bimolecular_match(); }
TEST(FindWhichReaction, MultipleRateStates) { test_fwr_multiple_rate_states(); }
TEST(FindWhichReaction, ScansAllListedReactions) { test_fwr_scans_all_listed_reactions(); }
TEST(FindWhichReaction, ImplicitLipidPartner) { test_fwr_implicit_lipid_partner(); }
TEST(FindWhichReaction, BimolecularStateChangeBackRxn) { test_fwr_bimolecular_state_change_back_rxn(); }
TEST(FindWhichReaction, StateChangeNoProductMatch) { test_fwr_state_change_no_product_match(); }
TEST(FindWhichReaction, UsesSuppliedInterfaceIndices) { test_fwr_uses_supplied_interface_indices(); }