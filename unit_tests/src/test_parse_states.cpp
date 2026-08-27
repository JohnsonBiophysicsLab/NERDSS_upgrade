/*! \file test_parse_states.cpp
 *
 * ### Unit test for ../src/parser/parse_states.cpp
 *
 * Function under test:
 *
 *     void parse_states(std::string& line, MolTemplate& molTemplate)
 *
 * Behaviour of the function (read directly from the implementation):
 *   1. Everything before the *first* `~` in `line` is taken as the interface
 *      name.  That prefix (and the `~`) is erased from `line`.
 *   2. Every remaining `~`-delimited token is collected as a state name; the
 *      final token (which has no trailing `~`) is added after the loop.
 *      NOTE: the trailing token is *not* erased from `line`, so on return
 *      `line` holds the text of the last state.
 *   3. The interface whose `name` matches the parsed interface name is located
 *      in `molTemplate.interfaceList`.  For every state token a
 *      `Interface::State(char iden, int index)` is **appended** to that
 *      interface's `stateList` with
 *          iden  = std::toupper(token[0])
 *          index = -1                    (a placeholder, per the TODO)
 *
 * IMPORTANT: if the interface name cannot be found the function calls
 * `exit(1)`.  That path is therefore deliberately *not* exercised here - doing
 * so would terminate the whole gtest binary.  Every test below hands
 * parse_states an interface name that really exists on the MolTemplate.
 */

#include "classes/class_MolTemplate.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// The prototype lives in "parser/parser_functions.hpp", but that header pulls in
// a large slice of the simulation (and split.cpp).  Only the single free
// function is required here, so it is declared directly to keep this
// translation unit small and free of link-time surprises.  The signature must
// match parser_functions.hpp exactly.
// -----------------------------------------------------------------------------
void parse_states(std::string& line, MolTemplate& molTemplate);

namespace {

/*! \brief Build a minimal MolTemplate holding a set of *stateless* interfaces.
 *
 * Each interface is given a name, a relative index and an empty stateList,
 * which is exactly the situation parse_states expects while a .mol file is
 * being read.
 *
 * \param[in] molName    name given to the molecule (used only for messages)
 * \param[in] ifaceNames list of interface names to create, in order
 * \return a fully initialised MolTemplate
 */
MolTemplate ps_make_template(const std::string& molName,
    const std::vector<std::string>& ifaceNames)
{
    MolTemplate molTemplate {};
    molTemplate.molName = molName;
    molTemplate.molTypeIndex = 0;

    for (unsigned i { 0 }; i < ifaceNames.size(); ++i) {
        // Interface(name, iCoord) leaves stateList empty, which is what we want.
        Interface oneIface { ifaceNames[i], Coord { 0.0, 0.0, 0.0 } };
        oneIface.index = static_cast<int>(i);
        molTemplate.interfaceList.push_back(oneIface);
    }
    return molTemplate;
}

/*! \brief Convenience printer so the console log shows the parsed states. */
void ps_dump_states(const Interface& iface)
{
    std::cerr << "    interface \"" << iface.name << "\" now has "
              << iface.stateList.size() << " state(s):";
    for (const auto& state : iface.stateList)
        std::cerr << " '" << state.iden << "'(index=" << state.index << ")";
    std::cerr << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the ordinary case - "iface~state1~state2" produces two states whose
//         identities are the upper-cased first characters of the tokens.
// -----------------------------------------------------------------------------
void test_parse_states_two_states()
{
    std::cerr << "\n[TEST] test_parse_states_two_states\n"
              << "  Source file: src/parser/parse_states.cpp\n"
              << "  Function:    parse_states()\n"
              << "  Input line:  \"a~U~P\" on a molecule with interfaces {a, b}\n"
              << "  Pass criteria: interface 'a' gains exactly 2 states with\n"
              << "                 idens 'U' and 'P', both with index == -1.\n";

    MolTemplate molTemplate { ps_make_template("pip2", { "a", "b" }) };
    std::string line { "a~U~P" };

    parse_states(line, molTemplate);
    ps_dump_states(molTemplate.interfaceList[0]);

    // Two state tokens were supplied, so two states must have been appended.
    EXPECT_EQ(molTemplate.interfaceList[0].stateList.size(), 2u)
        << "interface 'a' should have exactly two states";

    if (molTemplate.interfaceList[0].stateList.size() == 2u) {
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[0].iden, 'U')
            << "first state identity should be 'U'";
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[1].iden, 'P')
            << "second state identity should be 'P'";
        // The implementation uses -1 as a placeholder absolute index.
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[0].index, -1)
            << "state index is a placeholder and must be -1";
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[1].index, -1)
            << "state index is a placeholder and must be -1";
        // State(char,int) does not populate ifaceAndStateName.
        EXPECT_TRUE(molTemplate.interfaceList[0].stateList[0].ifaceAndStateName.empty())
            << "ifaceAndStateName is not set by parse_states";
    }

    // The other interface must be untouched.
    EXPECT_TRUE(molTemplate.interfaceList[1].stateList.empty())
        << "interface 'b' should not have received any states";
}

// -----------------------------------------------------------------------------
// Test 2: `line` is passed by reference and is consumed by the parser.  After
//         the call it holds only the *last* state token, because that token is
//         never erased inside the while-loop.
// -----------------------------------------------------------------------------
void test_parse_states_consumes_input_line()
{
    std::cerr << "\n[TEST] test_parse_states_consumes_input_line\n"
              << "  Source file: src/parser/parse_states.cpp\n"
              << "  Function:    parse_states()\n"
              << "  Scenario:    verify the side effect on the input string.\n"
              << "  Pass criteria: after parsing \"a~U~P\" the string equals \"P\"\n"
              << "                 (the trailing token is never erased).\n";

    MolTemplate molTemplate { ps_make_template("pip2", { "a" }) };
    std::string line { "a~U~P" };

    parse_states(line, molTemplate);

    std::cerr << "    line after the call = \"" << line << "\"\n";
    EXPECT_EQ(line, std::string { "P" })
        << "parse_states should leave the final state token in the input line";
}

// -----------------------------------------------------------------------------
// Test 3: multi-character, lower-case state names.  Only the first character is
//         used and it is upper-cased via std::toupper.
// -----------------------------------------------------------------------------
void test_parse_states_multichar_lowercase()
{
    std::cerr << "\n[TEST] test_parse_states_multichar_lowercase\n"
              << "  Source file: src/parser/parse_states.cpp\n"
              << "  Function:    parse_states()\n"
              << "  Input line:  \"site~unbound~phos\"\n"
              << "  Pass criteria: state identities are 'U' and 'P' (first char,\n"
              << "                 upper-cased) rather than whole words.\n";

    MolTemplate molTemplate { ps_make_template("kinase", { "site" }) };
    std::string line { "site~unbound~phos" };

    parse_states(line, molTemplate);
    ps_dump_states(molTemplate.interfaceList[0]);

    EXPECT_EQ(molTemplate.interfaceList[0].stateList.size(), 2u)
        << "two state tokens were provided";
    if (molTemplate.interfaceList[0].stateList.size() == 2u) {
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[0].iden, 'U')
            << "'unbound' must be reduced to its upper-cased first char 'U'";
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[1].iden, 'P')
            << "'phos' must be reduced to its upper-cased first char 'P'";
    }
}

// -----------------------------------------------------------------------------
// Test 4: a single state, i.e. "iface~state".  The while-loop body never runs;
//         the single token is added by the emplace_back after the loop.
// -----------------------------------------------------------------------------
void test_parse_states_single_state()
{
    std::cerr << "\n[TEST] test_parse_states_single_state\n"
              << "  Source file: src/parser/parse_states.cpp\n"
              << "  Function:    parse_states()\n"
              << "  Input line:  \"a~U\"\n"
              << "  Pass criteria: exactly one state with iden 'U' is appended and\n"
              << "                 the input line is reduced to \"U\".\n";

    MolTemplate molTemplate { ps_make_template("pip2", { "a" }) };
    std::string line { "a~U" };

    parse_states(line, molTemplate);
    ps_dump_states(molTemplate.interfaceList[0]);

    EXPECT_EQ(molTemplate.interfaceList[0].stateList.size(), 1u)
        << "one state token means one appended State";
    if (!molTemplate.interfaceList[0].stateList.empty()) {
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[0].iden, 'U')
            << "the single state identity should be 'U'";
    }
    EXPECT_EQ(line, std::string { "U" })
        << "the trailing token remains in the input line";
}

// -----------------------------------------------------------------------------
// Test 5: a line with no `~` at all.  find('~') returns npos, so the whole line
//         becomes the interface name; erase(0, npos + 1) wraps to erase(0, 0)
//         and therefore does nothing, and the (unchanged) line is then also
//         used as the single state token.  So "b" registers state 'B' on
//         interface "b" - an odd but well-defined behaviour worth pinning down.
// -----------------------------------------------------------------------------
void test_parse_states_no_tilde()
{
    std::cerr << "\n[TEST] test_parse_states_no_tilde\n"
              << "  Source file: src/parser/parse_states.cpp\n"
              << "  Function:    parse_states()\n"
              << "  Input line:  \"b\" (no '~' separator at all)\n"
              << "  Pass criteria: the whole line is used both as the interface\n"
              << "                 name and as the single state token, giving one\n"
              << "                 state with iden 'B'; the line is unchanged.\n";

    MolTemplate molTemplate { ps_make_template("pip2", { "a", "b" }) };
    std::string line { "b" };

    parse_states(line, molTemplate);
    ps_dump_states(molTemplate.interfaceList[1]);

    // Interface "a" is not the target and must remain stateless.
    EXPECT_TRUE(molTemplate.interfaceList[0].stateList.empty())
        << "interface 'a' should not be modified";

    EXPECT_EQ(molTemplate.interfaceList[1].stateList.size(), 1u)
        << "interface 'b' should have picked up exactly one state";
    if (!molTemplate.interfaceList[1].stateList.empty()) {
        EXPECT_EQ(molTemplate.interfaceList[1].stateList[0].iden, 'B')
            << "state identity comes from upper-casing the line's first char";
    }

    EXPECT_EQ(line, std::string { "b" })
        << "with no '~' present the input line is left untouched";
}

// -----------------------------------------------------------------------------
// Test 6: the target interface is chosen by *name*, not by position, and states
//         are appended to whatever is already present in stateList.
// -----------------------------------------------------------------------------
void test_parse_states_appends_to_correct_interface()
{
    std::cerr << "\n[TEST] test_parse_states_appends_to_correct_interface\n"
              << "  Source file: src/parser/parse_states.cpp\n"
              << "  Function:    parse_states()\n"
              << "  Scenario:    molecule with interfaces {a, b, c}; interface 'b'\n"
              << "               already carries one state before parsing.\n"
              << "  Pass criteria: only 'b' is modified and the new states are\n"
              << "                 appended after the pre-existing one.\n";

    MolTemplate molTemplate { ps_make_template("ap2", { "a", "b", "c" }) };

    // Pre-load interface 'b' with one state to prove the parser appends.
    molTemplate.interfaceList[1].stateList.emplace_back('X', 7);

    std::string line { "b~U~P" };
    parse_states(line, molTemplate);
    ps_dump_states(molTemplate.interfaceList[1]);

    // Neighbouring interfaces must be untouched.
    EXPECT_TRUE(molTemplate.interfaceList[0].stateList.empty())
        << "interface 'a' should be untouched";
    EXPECT_TRUE(molTemplate.interfaceList[2].stateList.empty())
        << "interface 'c' should be untouched";

    // 1 pre-existing + 2 parsed = 3 states, in that order.
    EXPECT_EQ(molTemplate.interfaceList[1].stateList.size(), 3u)
        << "parse_states appends rather than replaces the state list";
    if (molTemplate.interfaceList[1].stateList.size() == 3u) {
        EXPECT_EQ(molTemplate.interfaceList[1].stateList[0].iden, 'X')
            << "the pre-existing state must still be first";
        EXPECT_EQ(molTemplate.interfaceList[1].stateList[0].index, 7)
            << "the pre-existing state's index must not be altered";
        EXPECT_EQ(molTemplate.interfaceList[1].stateList[1].iden, 'U')
            << "first parsed state should follow the pre-existing one";
        EXPECT_EQ(molTemplate.interfaceList[1].stateList[2].iden, 'P')
            << "second parsed state should be last";
    }
}

// -----------------------------------------------------------------------------
// Test 7: non-alphabetic state tokens.  std::toupper leaves digits and symbols
//         unchanged, so they are stored verbatim as the state identity.
// -----------------------------------------------------------------------------
void test_parse_states_nonalpha_tokens()
{
    std::cerr << "\n[TEST] test_parse_states_nonalpha_tokens\n"
              << "  Source file: src/parser/parse_states.cpp\n"
              << "  Function:    parse_states()\n"
              << "  Input line:  \"a~1~2\"\n"
              << "  Pass criteria: digits pass through std::toupper unchanged, so\n"
              << "                 the identities are '1' and '2'.\n";

    MolTemplate molTemplate { ps_make_template("pip2", { "a" }) };
    std::string line { "a~1~2" };

    parse_states(line, molTemplate);
    ps_dump_states(molTemplate.interfaceList[0]);

    EXPECT_EQ(molTemplate.interfaceList[0].stateList.size(), 2u)
        << "two tokens should give two states";
    if (molTemplate.interfaceList[0].stateList.size() == 2u) {
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[0].iden, '1')
            << "digit '1' is unaffected by toupper";
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[1].iden, '2')
            << "digit '2' is unaffected by toupper";
    }
}

// -----------------------------------------------------------------------------
// Test 8: three or more states, to confirm the while-loop iterates correctly
//         and preserves ordering.
// -----------------------------------------------------------------------------
void test_parse_states_three_states_order()
{
    std::cerr << "\n[TEST] test_parse_states_three_states_order\n"
              << "  Source file: src/parser/parse_states.cpp\n"
              << "  Function:    parse_states()\n"
              << "  Input line:  \"a~U~P~D\"\n"
              << "  Pass criteria: three states appear in the same order as the\n"
              << "                 tokens: 'U', 'P', 'D'.\n";

    MolTemplate molTemplate { ps_make_template("pip2", { "a" }) };
    std::string line { "a~U~P~D" };

    parse_states(line, molTemplate);
    ps_dump_states(molTemplate.interfaceList[0]);

    EXPECT_EQ(molTemplate.interfaceList[0].stateList.size(), 3u)
        << "three tokens should give three states";
    if (molTemplate.interfaceList[0].stateList.size() == 3u) {
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[0].iden, 'U')
            << "state order must match token order (1st)";
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[1].iden, 'P')
            << "state order must match token order (2nd)";
        EXPECT_EQ(molTemplate.interfaceList[0].stateList[2].iden, 'D')
            << "state order must match token order (3rd)";
    }

    // The last token still sits in the input string.
    EXPECT_EQ(line, std::string { "D" })
        << "the final token is left in the input line";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so that
// a failure in one scenario does not prevent the remaining scenarios running.
// -----------------------------------------------------------------------------
TEST(ParseStates, TwoStates) { test_parse_states_two_states(); }
TEST(ParseStates, ConsumesInputLine) { test_parse_states_consumes_input_line(); }
TEST(ParseStates, MulticharLowercase) { test_parse_states_multichar_lowercase(); }
TEST(ParseStates, SingleState) { test_parse_states_single_state(); }
TEST(ParseStates, NoTilde) { test_parse_states_no_tilde(); }
TEST(ParseStates, AppendsToCorrectInterface) { test_parse_states_appends_to_correct_interface(); }
TEST(ParseStates, NonAlphaTokens) { test_parse_states_nonalpha_tokens(); }
TEST(ParseStates, ThreeStatesOrder) { test_parse_states_three_states_order(); }