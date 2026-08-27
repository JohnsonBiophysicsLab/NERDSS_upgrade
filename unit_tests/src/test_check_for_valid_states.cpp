/*! \file test_check_for_valid_states.cpp
 *
 * ### Unit test for src/parser/check_for_valid_states.cpp
 *
 * Function under test:
 *
 *     void check_for_valid_states(size_t parsedMolIndex,
 *                                 ParsedMol& targMol,
 *                                 ParsedRxn& parsedRxn,
 *                                 const std::vector<MolTemplate>& molTemplateList);
 *
 * Behaviour that is being verified here (read directly out of the
 * implementation, not guessed):
 *
 *   1. The MolTemplate whose `molName` matches `targMol.molName` is located and
 *      `targMol.molTypeIndex` is overwritten with that template's
 *      `molTypeIndex` (NOT the position of the template in the list).
 *
 *   2. For every interface of the parsed molecule the matching template
 *      Interface is located by name. Then:
 *        a. If the template interface has exactly ONE state, nothing at all
 *           happens - regardless of whether the parsed interface declared a
 *           state or not (both `if` branches test `stateList.size() != 1`).
 *        b. If the template interface has MORE than one state and the parsed
 *           interface declared no state ('\0'), the interface is recorded in
 *           `parsedRxn.noStateList` under the key `parsedMolIndex` so the
 *           reaction can later be split into one reaction per state.
 *        c. If the template interface has more than one state AND the parsed
 *           interface declared a state, the state is validated. For a
 *           `ReactionType::destruction` reaction the function only reports the
 *           match; for any other reaction type it delegates to
 *           check_for_state_change().
 *
 * Paths deliberately NOT exercised (they terminate the process and would kill
 * the whole gtest binary):
 *   - unknown molecule name  -> `exit(1)`
 *   - unknown interface name -> `exit(120)`
 *   - unknown state          -> `exit(1)`
 *
 * The "valid state on a non-destruction reaction" path is likewise only
 * exercised through the destruction branch, because the other branch calls
 * check_for_state_change(), a separate translation unit with its own
 * termination behaviour that is out of scope for this test.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

/*! \brief Build a MolTemplate with a chosen name, type index and interfaces.
 *
 * \param[in] molName      Name used by check_for_valid_states() to find the template.
 * \param[in] molTypeIndex The index the function is expected to copy into targMol.
 * \param[in] ifaceSpecs   One entry per interface: {interface name, string of
 *                         state identifiers}. An EMPTY state string means the
 *                         interface gets exactly one (unnamed) state, i.e. a
 *                         "stateless" interface.
 */
MolTemplate cfvs_make_template(const std::string& molName, int molTypeIndex,
    const std::vector<std::pair<std::string, std::string>>& ifaceSpecs)
{
    MolTemplate temp {};
    temp.molName = molName;
    temp.molTypeIndex = molTypeIndex;

    int absIndex { 0 }; // running "absolute" state index
    int relIndex { 0 }; // running relative interface index

    for (const auto& spec : ifaceSpecs) {
        Interface iface {};
        iface.name = spec.first;
        iface.index = relIndex++;
        iface.iCoord = Coord { 0.0, 0.0, 0.0 };

        if (spec.second.empty()) {
            // Single, state-less state: iden stays '\0'.
            Interface::State st {};
            st.iden = '\0';
            st.index = absIndex++;
            st.ifaceAndStateName = spec.first;
            iface.stateList.push_back(st);
        } else {
            // One State object per requested identifier character.
            for (char stateIden : spec.second) {
                Interface::State st {};
                st.iden = stateIden;
                st.index = absIndex++;
                st.ifaceAndStateName = spec.first + "~" + stateIden;
                iface.stateList.push_back(st);
            }
        }
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Build a ParsedMol (a molecule as read out of a reaction line).
 *
 * \param[in] molName      Molecule name, must match a template for the test to
 *                         avoid the exit(1) path.
 * \param[in] ifaces       One entry per interface: {interface name, declared
 *                         state}. Use '\0' for "no state declared".
 * \param[in] specieIndex  Which reactant species this molecule belongs to.
 */
ParsedMol cfvs_make_parsed_mol(const std::string& molName,
    const std::vector<std::pair<std::string, char>>& ifaces, int specieIndex)
{
    ParsedMol mol {};
    mol.molName = molName;
    mol.specieIndex = specieIndex;
    // molTypeIndex intentionally left at its default (-1) so the tests can show
    // that check_for_valid_states() is the thing that fills it in.
    for (const auto& oneIface : ifaces) {
        mol.interfaceList.emplace_back(oneIface.first, oneIface.second, /*isBound*/ false,
            /*bondIndex*/ -1, Involvement::possible, specieIndex);
    }
    return mol;
}

/*! \brief Build a minimal ParsedRxn whose type avoids check_for_state_change().
 *
 * ReactionType::destruction is used so the "valid explicit state" branch stops
 * after printing instead of delegating to another translation unit.
 */
ParsedRxn cfvs_make_destruction_rxn()
{
    ParsedRxn rxn {};
    rxn.rxnType = ReactionType::destruction;
    return rxn;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the molecule's type index is taken from the matching MolTemplate.
// -----------------------------------------------------------------------------
void test_cfvs_assigns_mol_type_index()
{
    std::cerr << "\n[TEST] test_cfvs_assigns_mol_type_index\n"
              << "  Source file: src/parser/check_for_valid_states.cpp\n"
              << "  Function:    check_for_valid_states()\n"
              << "  Scenario:    two templates whose molTypeIndex values differ from\n"
              << "               their position in molTemplateList; the parsed molecule\n"
              << "               matches the second one by name.\n"
              << "  Pass:        targMol.molTypeIndex == the template's molTypeIndex (3),\n"
              << "               not the list position (1).\n";

    // Template list: "alpha" is at position 0 but carries molTypeIndex 7,
    //                "beta"  is at position 1 but carries molTypeIndex 3.
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(cfvs_make_template("alpha", 7, { { "a1", "" } }));
    molTemplateList.push_back(cfvs_make_template("beta", 3, { { "b1", "" } }));

    // Parsed reactant "beta" with its single state-less interface.
    ParsedMol targMol = cfvs_make_parsed_mol("beta", { { "b1", '\0' } }, 0);
    ParsedRxn parsedRxn = cfvs_make_destruction_rxn();

    std::cerr << "  Before call: targMol.molTypeIndex = " << targMol.molTypeIndex << '\n';
    check_for_valid_states(0, targMol, parsedRxn, molTemplateList);
    std::cerr << "  After  call: targMol.molTypeIndex = " << targMol.molTypeIndex << '\n';

    EXPECT_EQ(targMol.molTypeIndex, 3)
        << "molTypeIndex must be copied from the matching MolTemplate, not the list index";

    // A state-less interface never produces a noStateList entry.
    EXPECT_TRUE(parsedRxn.noStateList.empty())
        << "A single-state interface must not be queued for reaction splitting";
}

// -----------------------------------------------------------------------------
// Test 2: interfaces that have exactly one state are always ignored.
// -----------------------------------------------------------------------------
void test_cfvs_single_state_interfaces_are_ignored()
{
    std::cerr << "\n[TEST] test_cfvs_single_state_interfaces_are_ignored\n"
              << "  Source file: src/parser/check_for_valid_states.cpp\n"
              << "  Function:    check_for_valid_states()\n"
              << "  Scenario:    a molecule with two interfaces, each having exactly one\n"
              << "               state in the template. One parsed interface declares no\n"
              << "               state, the other declares a state anyway.\n"
              << "  Pass:        noStateList stays empty in both cases, because both\n"
              << "               branches in the implementation require stateList.size() != 1.\n";

    std::vector<MolTemplate> molTemplateList;
    // Both "s1" and "s2" have a single, state-less state.
    molTemplateList.push_back(cfvs_make_template("solo", 0, { { "s1", "" }, { "s2", "" } }));

    // s1 declares no state, s2 declares 'Q' which does not even exist. Because
    // the template interface has only one state the function never validates it.
    ParsedMol targMol = cfvs_make_parsed_mol("solo", { { "s1", '\0' }, { "s2", 'Q' } }, 0);
    ParsedRxn parsedRxn = cfvs_make_destruction_rxn();

    check_for_valid_states(0, targMol, parsedRxn, molTemplateList);

    std::cerr << "  noStateList size after call = " << parsedRxn.noStateList.size() << '\n';
    EXPECT_EQ(parsedRxn.noStateList.size(), 0u)
        << "Single-state interfaces must never be added to noStateList";

    EXPECT_EQ(targMol.molTypeIndex, 0) << "molTypeIndex should still be assigned";

    // The parsed interfaces themselves must be untouched by this function.
    EXPECT_EQ(targMol.interfaceList[0].state, '\0')
        << "check_for_valid_states must not modify the declared state";
    EXPECT_EQ(targMol.interfaceList[1].state, 'Q')
        << "check_for_valid_states must not modify the declared state";
}

// -----------------------------------------------------------------------------
// Test 3: a multi-state interface with no declared state is queued for splitting.
// -----------------------------------------------------------------------------
void test_cfvs_multi_state_without_state_is_recorded()
{
    std::cerr << "\n[TEST] test_cfvs_multi_state_without_state_is_recorded\n"
              << "  Source file: src/parser/check_for_valid_states.cpp\n"
              << "  Function:    check_for_valid_states()\n"
              << "  Scenario:    template interface 'p' has states U and P; the parsed\n"
              << "               interface declares no state.\n"
              << "  Pass:        exactly one entry is inserted into parsedRxn.noStateList\n"
              << "               under key parsedMolIndex, and the stored IfaceInfo is a\n"
              << "               copy of the parsed interface.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(cfvs_make_template("kinase", 2, { { "p", "UP" } }));

    ParsedMol targMol = cfvs_make_parsed_mol("kinase", { { "p", '\0' } }, 0);
    ParsedRxn parsedRxn = cfvs_make_destruction_rxn();

    const size_t parsedMolIndex { 0 };
    check_for_valid_states(parsedMolIndex, targMol, parsedRxn, molTemplateList);

    std::cerr << "  noStateList size after call = " << parsedRxn.noStateList.size() << '\n';
    ASSERT_EQ(parsedRxn.noStateList.size(), 1u)
        << "A stateless declaration on a multi-state interface must be recorded once";

    // The key of the multimap entry is the parsed molecule index handed in.
    auto itr = parsedRxn.noStateList.begin();
    EXPECT_EQ(itr->first, static_cast<int>(parsedMolIndex))
        << "noStateList must be keyed by the parsedMolIndex argument";
    EXPECT_EQ(itr->second.ifaceName, std::string { "p" })
        << "The stored IfaceInfo must be a copy of the offending interface";
    EXPECT_EQ(itr->second.state, '\0')
        << "The stored IfaceInfo must still carry the empty state identifier";

    EXPECT_EQ(targMol.molTypeIndex, 2) << "molTypeIndex should be assigned from the template";
}

// -----------------------------------------------------------------------------
// Test 4: several stateless declarations on the same molecule all get recorded.
// -----------------------------------------------------------------------------
void test_cfvs_multiple_stateless_interfaces_all_recorded()
{
    std::cerr << "\n[TEST] test_cfvs_multiple_stateless_interfaces_all_recorded\n"
              << "  Source file: src/parser/check_for_valid_states.cpp\n"
              << "  Function:    check_for_valid_states()\n"
              << "  Scenario:    a molecule with three interfaces: two multi-state ones\n"
              << "               without a declared state and one single-state one.\n"
              << "  Pass:        noStateList holds exactly two entries, both under the\n"
              << "               same key, and the single-state interface is absent.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(cfvs_make_template("multi", 5,
        { { "m1", "AB" }, { "m2", "XYZ" }, { "plain", "" } }));

    ParsedMol targMol = cfvs_make_parsed_mol(
        "multi", { { "m1", '\0' }, { "m2", '\0' }, { "plain", '\0' } }, 0);
    ParsedRxn parsedRxn = cfvs_make_destruction_rxn();

    const size_t parsedMolIndex { 1 };
    check_for_valid_states(parsedMolIndex, targMol, parsedRxn, molTemplateList);

    std::cerr << "  noStateList size after call = " << parsedRxn.noStateList.size() << '\n';
    EXPECT_EQ(parsedRxn.noStateList.size(), 2u)
        << "Only the two multi-state interfaces should be queued";
    EXPECT_EQ(parsedRxn.noStateList.count(static_cast<int>(parsedMolIndex)), 2u)
        << "Both entries must share the parsedMolIndex key";

    // Verify the exact interface names that were recorded, and that 'plain'
    // (a single-state interface) was skipped.
    bool foundM1 { false };
    bool foundM2 { false };
    bool foundPlain { false };
    for (const auto& entry : parsedRxn.noStateList) {
        if (entry.second.ifaceName == "m1")
            foundM1 = true;
        if (entry.second.ifaceName == "m2")
            foundM2 = true;
        if (entry.second.ifaceName == "plain")
            foundPlain = true;
    }
    EXPECT_TRUE(foundM1) << "Interface m1 (states A,B) should be queued";
    EXPECT_TRUE(foundM2) << "Interface m2 (states X,Y,Z) should be queued";
    EXPECT_FALSE(foundPlain) << "Interface plain has one state and must not be queued";
}

// -----------------------------------------------------------------------------
// Test 5: an explicitly declared, valid state on a destruction reaction is
//         accepted and does NOT create a noStateList entry.
// -----------------------------------------------------------------------------
void test_cfvs_valid_declared_state_destruction()
{
    std::cerr << "\n[TEST] test_cfvs_valid_declared_state_destruction\n"
              << "  Source file: src/parser/check_for_valid_states.cpp\n"
              << "  Function:    check_for_valid_states()\n"
              << "  Scenario:    template interface 'p' has states U and P, the parsed\n"
              << "               interface declares 'P', and the reaction type is\n"
              << "               ReactionType::destruction.\n"
              << "  Pass:        the state is accepted (no termination), nothing is added\n"
              << "               to noStateList, and check_for_state_change() is skipped\n"
              << "               because the reaction is a destruction reaction.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(cfvs_make_template("kinase", 4, { { "p", "UP" } }));

    ParsedMol targMol = cfvs_make_parsed_mol("kinase", { { "p", 'P' } }, 0);
    ParsedRxn parsedRxn = cfvs_make_destruction_rxn();

    // Sanity: confirm the reaction really is a destruction reaction, otherwise
    // the implementation would call into check_for_state_change().
    ASSERT_EQ(parsedRxn.rxnType, ReactionType::destruction)
        << "Test set-up error: this test relies on the destruction branch";

    check_for_valid_states(0, targMol, parsedRxn, molTemplateList);

    std::cerr << "  noStateList size after call = " << parsedRxn.noStateList.size() << '\n';
    EXPECT_TRUE(parsedRxn.noStateList.empty())
        << "An explicitly declared valid state must not be queued for splitting";
    EXPECT_EQ(targMol.molTypeIndex, 4) << "molTypeIndex should be assigned from the template";
    EXPECT_EQ(targMol.interfaceList[0].state, 'P')
        << "The declared state must be left untouched";

    // A destruction reaction must not have been flagged as a state change here.
    EXPECT_FALSE(parsedRxn.hasStateChange)
        << "check_for_state_change() must be skipped for destruction reactions";
}

// -----------------------------------------------------------------------------
// Test 6: entries from separate molecules accumulate under distinct keys.
// -----------------------------------------------------------------------------
void test_cfvs_entries_accumulate_across_calls()
{
    std::cerr << "\n[TEST] test_cfvs_entries_accumulate_across_calls\n"
              << "  Source file: src/parser/check_for_valid_states.cpp\n"
              << "  Function:    check_for_valid_states()\n"
              << "  Scenario:    the function is called twice on the same ParsedRxn with\n"
              << "               different parsedMolIndex values (reactant 0 and 1).\n"
              << "  Pass:        noStateList ends up with one entry per call, filed under\n"
              << "               the respective parsedMolIndex key.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(cfvs_make_template("first", 0, { { "f", "AB" } }));
    molTemplateList.push_back(cfvs_make_template("second", 1, { { "s", "CD" } }));

    ParsedMol mol0 = cfvs_make_parsed_mol("first", { { "f", '\0' } }, 0);
    ParsedMol mol1 = cfvs_make_parsed_mol("second", { { "s", '\0' } }, 1);
    ParsedRxn parsedRxn = cfvs_make_destruction_rxn();

    check_for_valid_states(0, mol0, parsedRxn, molTemplateList);
    std::cerr << "  after first call, noStateList size = " << parsedRxn.noStateList.size() << '\n';
    EXPECT_EQ(parsedRxn.noStateList.size(), 1u) << "First call should add exactly one entry";

    check_for_valid_states(1, mol1, parsedRxn, molTemplateList);
    std::cerr << "  after second call, noStateList size = " << parsedRxn.noStateList.size() << '\n';
    EXPECT_EQ(parsedRxn.noStateList.size(), 2u)
        << "Second call should append rather than replace";

    // Each parsed molecule index owns exactly one entry.
    EXPECT_EQ(parsedRxn.noStateList.count(0), 1u) << "Key 0 should hold the 'first' interface";
    EXPECT_EQ(parsedRxn.noStateList.count(1), 1u) << "Key 1 should hold the 'second' interface";

    // Confirm the right interface landed under the right key.
    auto itr0 = parsedRxn.noStateList.find(0);
    ASSERT_NE(itr0, parsedRxn.noStateList.end()) << "Key 0 must exist";
    EXPECT_EQ(itr0->second.ifaceName, std::string { "f" })
        << "Interface 'f' belongs to parsed molecule 0";

    auto itr1 = parsedRxn.noStateList.find(1);
    ASSERT_NE(itr1, parsedRxn.noStateList.end()) << "Key 1 must exist";
    EXPECT_EQ(itr1->second.ifaceName, std::string { "s" })
        << "Interface 's' belongs to parsed molecule 1";

    // Both molecules should have picked up their template type index.
    EXPECT_EQ(mol0.molTypeIndex, 0) << "mol0 should be typed as template 'first'";
    EXPECT_EQ(mol1.molTypeIndex, 1) << "mol1 should be typed as template 'second'";
}

// -----------------------------------------------------------------------------
// Test 7: a molecule with no interfaces at all only gets its type index set.
// -----------------------------------------------------------------------------
void test_cfvs_empty_interface_list()
{
    std::cerr << "\n[TEST] test_cfvs_empty_interface_list\n"
              << "  Source file: src/parser/check_for_valid_states.cpp\n"
              << "  Function:    check_for_valid_states()\n"
              << "  Scenario:    the parsed molecule declares no interfaces at all, so the\n"
              << "               interface loop body never executes.\n"
              << "  Pass:        molTypeIndex is still assigned and noStateList stays empty.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(cfvs_make_template("bare", 9, { { "onlyIface", "AB" } }));

    ParsedMol targMol = cfvs_make_parsed_mol("bare", {}, 0);
    ParsedRxn parsedRxn = cfvs_make_destruction_rxn();

    check_for_valid_states(0, targMol, parsedRxn, molTemplateList);

    EXPECT_EQ(targMol.molTypeIndex, 9)
        << "The type index lookup happens before the interface loop";
    EXPECT_TRUE(parsedRxn.noStateList.empty())
        << "No interfaces means nothing can be queued for splitting";
    EXPECT_TRUE(targMol.interfaceList.empty())
        << "The function must not invent interfaces";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_cfvs_* helper is invoked from its own
// TEST so a failure in one does not stop the others from running.
// -----------------------------------------------------------------------------
TEST(CheckForValidStates, AssignsMolTypeIndex) { test_cfvs_assigns_mol_type_index(); }
TEST(CheckForValidStates, SingleStateInterfacesAreIgnored) { test_cfvs_single_state_interfaces_are_ignored(); }
TEST(CheckForValidStates, MultiStateWithoutStateIsRecorded) { test_cfvs_multi_state_without_state_is_recorded(); }
TEST(CheckForValidStates, MultipleStatelessInterfacesAllRecorded) { test_cfvs_multiple_stateless_interfaces_all_recorded(); }
TEST(CheckForValidStates, ValidDeclaredStateDestruction) { test_cfvs_valid_declared_state_destruction(); }
TEST(CheckForValidStates, EntriesAccumulateAcrossCalls) { test_cfvs_entries_accumulate_across_calls(); }
TEST(CheckForValidStates, EmptyInterfaceList) { test_cfvs_empty_interface_list(); }