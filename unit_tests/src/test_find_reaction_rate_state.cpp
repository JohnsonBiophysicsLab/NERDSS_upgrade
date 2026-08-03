/*! \file test_find_reaction_rate_state.cpp
 *
 * ### Unit test for src/reactions/find_reaction_rate_state.cpp
 *
 * Function under test:
 *
 *     int find_reaction_rate_state(int simItr, int relIfaceIndex1, int relIfaceIndex2,
 *                                  const Molecule& reactMol1, const Molecule& reactMol2,
 *                                  const BackRxn& backRxn,
 *                                  const std::vector<MolTemplate>& molTemplateList);
 *
 * What the function does:
 *   1. It walks BackRxn::reactantListNew and figures out which entry of that list
 *      corresponds to reactMol1/relIfaceIndex1 (reactIndex1) and which corresponds to
 *      reactMol2/relIfaceIndex2 (reactIndex2).  If either cannot be located the function
 *      prints an error and calls exit(1) -- every test below therefore always supplies a
 *      reactant list that *does* contain both reacting interfaces, so we never hit that
 *      fatal branch.
 *   2. It then loops over BackRxn::rateList and asks hasIntangibles() whether the two
 *      molecules carry the ancillary (required) interfaces of each rate state.
 *        - exactly one matching rate state  -> return that rate state's index
 *        - no matching rate state           -> return -1
 *        - more than one matching rate state-> return the index of the match with the
 *                                              largest number of ancillary interfaces
 *                                              (ties resolved in favour of the first one)
 *
 * The tests below exercise each of those outcomes and print, verbosely, which source
 * file/function is being probed and what the pass criterion is.
 */

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "classes/class_Rxns.hpp"
#include "reactions/shared_reaction_functions.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Small builders used by all tests.  All names are prefixed with "frrs_"
// (find_reaction_rate_state) so they cannot collide with other test files.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Create one Molecule::Iface with fully consistent bookkeeping fields.
 *
 * All identity fields (absolute index, relative index, parent molecule type, state
 * character and bound status) are set explicitly so that a matching RxnIface can be
 * created with exactly the same values -- this makes interface comparison unambiguous.
 */
Molecule::Iface frrs_make_iface(int absIndex, int relIndex, int molTypeIndex, char state, bool isBound)
{
    Molecule::Iface iface;
    iface.coord = Coord { 0.0, 0.0, 0.0 };
    iface.stateIden = state;
    iface.stateIndex = 0;
    iface.index = absIndex; // absolute interface-state index
    iface.relIndex = relIndex; // index within Molecule::interfaceList
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = isBound;
    return iface;
}

/*! \brief Create a Molecule of a given type holding the supplied interface list. */
Molecule frrs_make_molecule(int molTypeIndex, int molIndex, const std::vector<Molecule::Iface>& ifaceList)
{
    Molecule mol;
    mol.molTypeIndex = molTypeIndex;
    mol.index = molIndex;
    mol.myComIndex = molIndex;
    mol.comCoord = Coord { 0.0, 0.0, 0.0 };
    mol.interfaceList = ifaceList;
    return mol;
}

/*! \brief Create a reaction-side interface descriptor (RxnIface).
 *
 * Every identity field is filled in so the descriptor can be made to either match or
 * deliberately mismatch a Molecule::Iface built by frrs_make_iface().
 */
RxnIface frrs_make_rxn_iface(
    const std::string& name, int molTypeIndex, int absIndex, int relIndex, char state, bool requiresInteraction)
{
    RxnIface rxnIface;
    rxnIface.ifaceName = name;
    rxnIface.molTypeIndex = molTypeIndex;
    rxnIface.absIfaceIndex = absIndex;
    rxnIface.relIfaceIndex = relIndex;
    rxnIface.requiresState = state;
    rxnIface.requiresInteraction = requiresInteraction;
    return rxnIface;
}

/*! \brief Create a RxnBase::RateState with two ancillary-interface lists.
 *
 * find_reaction_rate_state() always indexes otherIfaceLists[0] and [1], so both lists
 * must exist even when they are empty.
 */
RxnBase::RateState frrs_make_rate_state(
    double rate, const std::vector<RxnIface>& ancillary1, const std::vector<RxnIface>& ancillary2)
{
    std::vector<std::vector<RxnIface>> otherIfaceLists { ancillary1, ancillary2 };
    RxnBase::RateState rateState(rate, otherIfaceLists);
    return rateState;
}

/*! \brief Assemble a minimal but self-consistent BackRxn. */
BackRxn frrs_make_back_rxn(const std::vector<RxnIface>& reactants, const std::vector<RxnBase::RateState>& rates)
{
    BackRxn backRxn;
    backRxn.rxnType = ReactionType::bimolecular;
    backRxn.absRxnIndex = 0;
    backRxn.relRxnIndex = 0;
    backRxn.reactantListNew = reactants;
    backRxn.rateList = rates;
    return backRxn;
}

/*! \brief Two dummy MolTemplates.
 *
 * These are only consumed by the fatal error branch of the function (which we never
 * trigger), but they are built properly anyway so the argument is always valid.
 */
std::vector<MolTemplate> frrs_make_templates()
{
    std::vector<MolTemplate> molTemplateList(2);
    for (int i = 0; i < 2; ++i) {
        molTemplateList[i].molTypeIndex = i;
        molTemplateList[i].molName = (i == 0) ? "A" : "B";
        molTemplateList[i].radius = 1.0;

        Interface iface;
        iface.index = 0;
        iface.name = (i == 0) ? "a1" : "b1";
        Interface::State state;
        state.iden = 'U';
        state.index = 10 * (i + 1);
        state.ifaceAndStateName = iface.name + "~U";
        iface.stateList.push_back(state);
        molTemplateList[i].interfaceList.push_back(iface);
    }
    return molTemplateList;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: exactly one rate state, no ancillary interfaces required.
// -----------------------------------------------------------------------------
void test_frrs_single_rate_state_no_ancillary()
{
    std::cerr << "\n[TEST] test_frrs_single_rate_state_no_ancillary\n"
              << "  Source file:   src/reactions/find_reaction_rate_state.cpp\n"
              << "  Function:      find_reaction_rate_state(BackRxn overload)\n"
              << "  Scenario:      one rate state whose ancillary interface lists are both\n"
              << "                 empty, so the rate state trivially applies.\n"
              << "  Pass criteria: the function returns rate-state index 0.\n";

    const std::vector<MolTemplate> molTemplateList = frrs_make_templates();

    // Molecule 1 (type 0) with a single reacting interface (relIndex 0, absIndex 10).
    Molecule reactMol1 = frrs_make_molecule(0, 0, { frrs_make_iface(10, 0, 0, 'U', true) });
    // Molecule 2 (type 1) with a single reacting interface (relIndex 0, absIndex 20).
    Molecule reactMol2 = frrs_make_molecule(1, 1, { frrs_make_iface(20, 0, 1, 'U', true) });

    // The dissociation reaction lists both reacting interfaces, mol1 first.
    std::vector<RxnIface> reactants { frrs_make_rxn_iface("a1", 0, 10, 0, 'U', true),
        frrs_make_rxn_iface("b1", 1, 20, 0, 'U', true) };

    // One rate state, no ancillary requirements at all.
    std::vector<RxnBase::RateState> rates { frrs_make_rate_state(1.0, {}, {}) };

    BackRxn backRxn = frrs_make_back_rxn(reactants, rates);

    std::cerr << "  Calling find_reaction_rate_state(simItr=0, relIface1=0, relIface2=0, ...)\n";
    const int result = find_reaction_rate_state(0, 0, 0, reactMol1, reactMol2, backRxn, molTemplateList);
    std::cerr << "  Returned rate-state index = " << result << " (expected 0)\n";

    EXPECT_EQ(result, 0) << "With a single unconditional rate state, index 0 must be returned";
}

// -----------------------------------------------------------------------------
// Test 2: the only rate state demands an ancillary interface the molecules lack.
// -----------------------------------------------------------------------------
void test_frrs_no_matching_rate_state()
{
    std::cerr << "\n[TEST] test_frrs_no_matching_rate_state\n"
              << "  Source file:   src/reactions/find_reaction_rate_state.cpp\n"
              << "  Function:      find_reaction_rate_state(BackRxn overload)\n"
              << "  Scenario:      the single rate state requires an ancillary interface\n"
              << "                 (absIndex 99, state 'Z') that neither molecule owns.\n"
              << "  Pass criteria: the function reports 'no match' by returning -1.\n";

    const std::vector<MolTemplate> molTemplateList = frrs_make_templates();

    Molecule reactMol1 = frrs_make_molecule(0, 0, { frrs_make_iface(10, 0, 0, 'U', true) });
    Molecule reactMol2 = frrs_make_molecule(1, 1, { frrs_make_iface(20, 0, 1, 'U', true) });

    std::vector<RxnIface> reactants { frrs_make_rxn_iface("a1", 0, 10, 0, 'U', true),
        frrs_make_rxn_iface("b1", 1, 20, 0, 'U', true) };

    // Ancillary requirement on the mol1 side that cannot possibly be satisfied.
    std::vector<RxnIface> impossible { frrs_make_rxn_iface("aX", 0, 99, 7, 'Z', true) };
    std::vector<RxnBase::RateState> rates { frrs_make_rate_state(1.0, impossible, {}) };

    BackRxn backRxn = frrs_make_back_rxn(reactants, rates);

    std::cerr << "  Calling find_reaction_rate_state with an unsatisfiable ancillary list\n";
    const int result = find_reaction_rate_state(0, 0, 0, reactMol1, reactMol2, backRxn, molTemplateList);
    std::cerr << "  Returned rate-state index = " << result << " (expected -1)\n";

    EXPECT_EQ(result, -1) << "When no rate state's ancillary interfaces are present, -1 is expected";
}

// -----------------------------------------------------------------------------
// Test 3: several matching rate states -> the most specific one (most ancillary
//         interfaces) must win.
// -----------------------------------------------------------------------------
void test_frrs_prefers_most_ancillary_interfaces()
{
    std::cerr << "\n[TEST] test_frrs_prefers_most_ancillary_interfaces\n"
              << "  Source file:   src/reactions/find_reaction_rate_state.cpp\n"
              << "  Function:      find_reaction_rate_state(BackRxn overload)\n"
              << "  Scenario:      rate state 0 is unconditional (0 ancillary ifaces) and\n"
              << "                 rate state 1 requires one ancillary interface that mol1\n"
              << "                 actually carries -- both therefore match.\n"
              << "  Pass criteria: the more specific rate state (index 1) is chosen.\n";

    const std::vector<MolTemplate> molTemplateList = frrs_make_templates();

    // mol1 owns the reacting interface (rel 0 / abs 10) AND an extra interface
    // (rel 1 / abs 11, state 'U', unbound) that will be used as the ancillary one.
    Molecule reactMol1 = frrs_make_molecule(
        0, 0, { frrs_make_iface(10, 0, 0, 'U', true), frrs_make_iface(11, 1, 0, 'U', false) });
    Molecule reactMol2 = frrs_make_molecule(1, 1, { frrs_make_iface(20, 0, 1, 'U', true) });

    std::vector<RxnIface> reactants { frrs_make_rxn_iface("a1", 0, 10, 0, 'U', true),
        frrs_make_rxn_iface("b1", 1, 20, 0, 'U', true) };

    // Ancillary descriptor whose every field mirrors mol1's second interface.
    std::vector<RxnIface> ancillaryOnMol1 { frrs_make_rxn_iface("a2", 0, 11, 1, 'U', false) };

    std::vector<RxnBase::RateState> rates {
        frrs_make_rate_state(1.0, {}, {}), // index 0: generic, 0 ancillary interfaces
        frrs_make_rate_state(5.0, ancillaryOnMol1, {}) // index 1: specific, 1 ancillary interface
    };

    BackRxn backRxn = frrs_make_back_rxn(reactants, rates);

    std::cerr << "  Calling find_reaction_rate_state with 2 candidate rate states\n";
    const int result = find_reaction_rate_state(0, 0, 0, reactMol1, reactMol2, backRxn, molTemplateList);
    std::cerr << "  Returned rate-state index = " << result << " (expected 1, the specific state)\n";

    EXPECT_EQ(result, 1) << "Among multiple matches the state with the most ancillary interfaces must win";
}

// -----------------------------------------------------------------------------
// Test 4: several matching rate states, all equally unspecific -> first one wins.
// -----------------------------------------------------------------------------
void test_frrs_tie_returns_first_match()
{
    std::cerr << "\n[TEST] test_frrs_tie_returns_first_match\n"
              << "  Source file:   src/reactions/find_reaction_rate_state.cpp\n"
              << "  Function:      find_reaction_rate_state(BackRxn overload)\n"
              << "  Scenario:      two rate states, both with empty ancillary lists, so both\n"
              << "                 match and both have the same 'specificity' of zero.\n"
              << "  Pass criteria: the tie is broken deterministically in favour of index 0.\n";

    const std::vector<MolTemplate> molTemplateList = frrs_make_templates();

    Molecule reactMol1 = frrs_make_molecule(0, 0, { frrs_make_iface(10, 0, 0, 'U', true) });
    Molecule reactMol2 = frrs_make_molecule(1, 1, { frrs_make_iface(20, 0, 1, 'U', true) });

    std::vector<RxnIface> reactants { frrs_make_rxn_iface("a1", 0, 10, 0, 'U', true),
        frrs_make_rxn_iface("b1", 1, 20, 0, 'U', true) };

    std::vector<RxnBase::RateState> rates { frrs_make_rate_state(1.0, {}, {}), frrs_make_rate_state(2.0, {}, {}) };

    BackRxn backRxn = frrs_make_back_rxn(reactants, rates);

    std::cerr << "  Calling find_reaction_rate_state with 2 equally generic rate states\n";
    const int result = find_reaction_rate_state(0, 0, 0, reactMol1, reactMol2, backRxn, molTemplateList);
    std::cerr << "  Returned rate-state index = " << result << " (expected 0)\n";

    EXPECT_EQ(result, 0) << "A tie in ancillary-interface count must resolve to the first match";
}

// -----------------------------------------------------------------------------
// Test 5: the reactant list stores the molecules in the opposite order.  The
//         reactant-index bookkeeping inside the function must cope with this.
// -----------------------------------------------------------------------------
void test_frrs_reversed_reactant_order()
{
    std::cerr << "\n[TEST] test_frrs_reversed_reactant_order\n"
              << "  Source file:   src/reactions/find_reaction_rate_state.cpp\n"
              << "  Function:      find_reaction_rate_state(BackRxn overload)\n"
              << "  Scenario:      BackRxn::reactantListNew lists molecule 2 first, so the\n"
              << "                 internally derived reactIndex1/reactIndex2 are swapped.\n"
              << "  Pass criteria: the reactants are still located (no exit) and the single\n"
              << "                 unconditional rate state (index 0) is returned.\n";

    const std::vector<MolTemplate> molTemplateList = frrs_make_templates();

    Molecule reactMol1 = frrs_make_molecule(0, 0, { frrs_make_iface(10, 0, 0, 'U', true) });
    Molecule reactMol2 = frrs_make_molecule(1, 1, { frrs_make_iface(20, 0, 1, 'U', true) });

    // Deliberately reversed: entry 0 belongs to reactMol2, entry 1 to reactMol1.
    std::vector<RxnIface> reactants { frrs_make_rxn_iface("b1", 1, 20, 0, 'U', true),
        frrs_make_rxn_iface("a1", 0, 10, 0, 'U', true) };

    std::vector<RxnBase::RateState> rates { frrs_make_rate_state(3.0, {}, {}) };

    BackRxn backRxn = frrs_make_back_rxn(reactants, rates);

    std::cerr << "  Calling find_reaction_rate_state with a reversed reactant list\n";
    const int result = find_reaction_rate_state(0, 0, 0, reactMol1, reactMol2, backRxn, molTemplateList);
    std::cerr << "  Returned rate-state index = " << result << " (expected 0)\n";

    EXPECT_EQ(result, 0) << "Reactant ordering inside the reaction must not change the returned index";
}

// -----------------------------------------------------------------------------
// Test 6: symmetric reaction (identical molecule types and identical relative
//         interface indices for both reactants).
// -----------------------------------------------------------------------------
void test_frrs_symmetric_reaction()
{
    std::cerr << "\n[TEST] test_frrs_symmetric_reaction\n"
              << "  Source file:   src/reactions/find_reaction_rate_state.cpp\n"
              << "  Function:      find_reaction_rate_state(BackRxn overload)\n"
              << "  Scenario:      a homodimer dissociation: both reactants are type 0 and\n"
              << "                 use relative interface 0, so the two reactant-list entries\n"
              << "                 are indistinguishable and must be assigned one each.\n"
              << "  Pass criteria: both reactant indices are found (no exit(1)) and rate\n"
              << "                 state 0 is returned.\n";

    const std::vector<MolTemplate> molTemplateList = frrs_make_templates();

    // Two copies of the same molecule type bound through the same interface.
    Molecule reactMol1 = frrs_make_molecule(0, 0, { frrs_make_iface(10, 0, 0, 'U', true) });
    Molecule reactMol2 = frrs_make_molecule(0, 1, { frrs_make_iface(10, 0, 0, 'U', true) });

    std::vector<RxnIface> reactants { frrs_make_rxn_iface("a1", 0, 10, 0, 'U', true),
        frrs_make_rxn_iface("a1", 0, 10, 0, 'U', true) };

    std::vector<RxnBase::RateState> rates { frrs_make_rate_state(4.0, {}, {}) };

    BackRxn backRxn = frrs_make_back_rxn(reactants, rates);

    std::cerr << "  Calling find_reaction_rate_state for a symmetric (homodimer) back reaction\n";
    const int result = find_reaction_rate_state(0, 0, 0, reactMol1, reactMol2, backRxn, molTemplateList);
    std::cerr << "  Returned rate-state index = " << result << " (expected 0)\n";

    EXPECT_EQ(result, 0) << "A symmetric back reaction with one rate state must return index 0";
}

// -----------------------------------------------------------------------------
// Test 7: a rate state that requires ancillary interfaces on *both* reactants.
// -----------------------------------------------------------------------------
void test_frrs_ancillary_on_both_reactants()
{
    std::cerr << "\n[TEST] test_frrs_ancillary_on_both_reactants\n"
              << "  Source file:   src/reactions/find_reaction_rate_state.cpp\n"
              << "  Function:      find_reaction_rate_state(BackRxn overload)\n"
              << "  Scenario:      rate state 0 requires one ancillary interface on mol1 and\n"
              << "                 one on mol2; mol1 has its requirement, mol2 does NOT.\n"
              << "  Pass criteria: the partial match is rejected, so -1 is returned.\n";

    const std::vector<MolTemplate> molTemplateList = frrs_make_templates();

    // mol1 has the ancillary interface (abs 11); mol2 only has its reacting interface.
    Molecule reactMol1 = frrs_make_molecule(
        0, 0, { frrs_make_iface(10, 0, 0, 'U', true), frrs_make_iface(11, 1, 0, 'U', false) });
    Molecule reactMol2 = frrs_make_molecule(1, 1, { frrs_make_iface(20, 0, 1, 'U', true) });

    std::vector<RxnIface> reactants { frrs_make_rxn_iface("a1", 0, 10, 0, 'U', true),
        frrs_make_rxn_iface("b1", 1, 20, 0, 'U', true) };

    std::vector<RxnIface> ancillaryOnMol1 { frrs_make_rxn_iface("a2", 0, 11, 1, 'U', false) };
    // mol2 is asked for an interface (abs 77) it does not have at all.
    std::vector<RxnIface> ancillaryOnMol2 { frrs_make_rxn_iface("b2", 1, 77, 3, 'P', true) };

    std::vector<RxnBase::RateState> rates { frrs_make_rate_state(9.0, ancillaryOnMol1, ancillaryOnMol2) };

    BackRxn backRxn = frrs_make_back_rxn(reactants, rates);

    std::cerr << "  Calling find_reaction_rate_state where only half the ancillary set is present\n";
    const int result = find_reaction_rate_state(0, 0, 0, reactMol1, reactMol2, backRxn, molTemplateList);
    std::cerr << "  Returned rate-state index = " << result << " (expected -1)\n";

    EXPECT_EQ(result, -1) << "All ancillary interfaces of a rate state must be present for it to match";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario runs in its own TEST so that a failure in
// one does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(FindReactionRateState, SingleRateStateNoAncillary) { test_frrs_single_rate_state_no_ancillary(); }
TEST(FindReactionRateState, NoMatchingRateState) { test_frrs_no_matching_rate_state(); }
TEST(FindReactionRateState, PrefersMostAncillaryInterfaces) { test_frrs_prefers_most_ancillary_interfaces(); }
TEST(FindReactionRateState, TieReturnsFirstMatch) { test_frrs_tie_returns_first_match(); }
TEST(FindReactionRateState, ReversedReactantOrder) { test_frrs_reversed_reactant_order(); }
TEST(FindReactionRateState, SymmetricReaction) { test_frrs_symmetric_reaction(); }
TEST(FindReactionRateState, AncillaryOnBothReactants) { test_frrs_ancillary_on_both_reactants(); }