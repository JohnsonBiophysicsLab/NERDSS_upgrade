/*! \file test_areSameExceptState.cpp
 *
 * ### Unit test for ../src/parser/areSameExceptState.cpp
 *
 * The single function under test is:
 *
 * \code
 * bool areSameExceptState(const ParsedMol::IfaceInfo& iface1,
 *                         const ParsedMol::IfaceInfo& iface2);
 * \endcode
 *
 * It returns true only when **all four** of the following hold:
 *   1. `ifaceName`     of both interfaces are equal (std::string comparison),
 *   2. `isBound`       of both interfaces are equal,
 *   3. `state`         of the two interfaces are **different**,
 *   4. `speciesIndex`  of both interfaces are equal.
 *
 * Every other data member of ParsedMol::IfaceInfo (molTypeIndex, absIndex,
 * relIndex, bondIndex, indexFound, ifaceAndStateName, ifaceRxnStatus) is
 * deliberately ignored by the function; the tests below verify that as well.
 *
 * Verbose progress messages are written to stderr so that the reader can follow
 * exactly which source file, which function, and which criterion is being
 * exercised by each assertion.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with "assxs_" so they cannot collide with helpers
// belonging to other translation units of the combined gtest binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully-initialised ParsedMol::IfaceInfo.
 *
 * The default constructor of IfaceInfo leaves every member at its in-class
 * initialiser, which is already well defined, so we simply start from a
 * default-constructed object and then assign the four members the function
 * under test actually inspects.  All remaining members are left at their
 * defaults unless a test explicitly overrides them.
 *
 * \param[in] name          value for IfaceInfo::ifaceName
 * \param[in] state         value for IfaceInfo::state ('\0' means "no state")
 * \param[in] isBound       value for IfaceInfo::isBound
 * \param[in] speciesIndex  value for IfaceInfo::speciesIndex
 * \return A populated IfaceInfo object.
 */
ParsedMol::IfaceInfo assxs_make_iface(const std::string& name, char state, bool isBound, int speciesIndex)
{
    ParsedMol::IfaceInfo iface {};
    iface.ifaceName = name;
    iface.state = state;
    iface.isBound = isBound;
    iface.speciesIndex = speciesIndex;

    // The following members are NOT read by areSameExceptState(); they are set
    // here only so the object is never left partially uninitialised.
    iface.molTypeIndex = 0;
    iface.absIndex = 0;
    iface.relIndex = 0;
    iface.bondIndex = -1;
    iface.indexFound = false;
    iface.ifaceRxnStatus = Involvement::none;
    iface.ifaceAndStateName = (state != '\0') ? (name + "~" + state) : name;

    return iface;
}

/*! \brief Small convenience printer used by the verbose console output. */
void assxs_describe(const char* label, const ParsedMol::IfaceInfo& iface)
{
    std::cerr << "      " << label << ": name=\"" << iface.ifaceName << "\", state='"
              << (iface.state == '\0' ? '0' : iface.state) << "', isBound=" << std::boolalpha << iface.isBound
              << ", speciesIndex=" << iface.speciesIndex << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: The canonical "true" case -- everything matches except the state.
// -----------------------------------------------------------------------------
void test_assxs_true_when_only_state_differs()
{
    std::cerr << "\n[TEST] test_assxs_true_when_only_state_differs\n"
              << "  Source file: src/parser/areSameExceptState.cpp\n"
              << "  Function:    areSameExceptState()\n"
              << "  Scenario:    two interfaces with identical name / isBound /\n"
              << "               speciesIndex but different state characters.\n"
              << "  Pass:        the function must return true.\n";

    // "a" on species 0, unbound, in state U versus the same interface in state P.
    ParsedMol::IfaceInfo unphos = assxs_make_iface("a", 'U', false, 0);
    ParsedMol::IfaceInfo phos = assxs_make_iface("a", 'P', false, 0);
    assxs_describe("iface1", unphos);
    assxs_describe("iface2", phos);

    EXPECT_TRUE(areSameExceptState(unphos, phos))
        << "Interfaces that differ only by state must be reported as 'same except state'";

    // The comparison is completely symmetric, so swapping arguments must agree.
    std::cerr << "  -> Verifying the comparison is symmetric (swap the arguments)\n";
    EXPECT_TRUE(areSameExceptState(phos, unphos))
        << "areSameExceptState must give the same answer with swapped arguments";
}

// -----------------------------------------------------------------------------
// Test 2: Identical objects -- the state check (state1 != state2) fails.
// -----------------------------------------------------------------------------
void test_assxs_false_when_states_are_identical()
{
    std::cerr << "\n[TEST] test_assxs_false_when_states_are_identical\n"
              << "  Function:    areSameExceptState()\n"
              << "  Scenario:    two interfaces that are byte-for-byte identical,\n"
              << "               including the state character.\n"
              << "  Pass:        the function must return false, because the\n"
              << "               predicate explicitly requires state1 != state2.\n";

    ParsedMol::IfaceInfo iface1 = assxs_make_iface("a", 'U', false, 0);
    ParsedMol::IfaceInfo iface2 = assxs_make_iface("a", 'U', false, 0);
    assxs_describe("iface1", iface1);
    assxs_describe("iface2", iface2);

    EXPECT_FALSE(areSameExceptState(iface1, iface2))
        << "Identical interfaces are NOT 'same except state' -- the states must differ";

    // A self-comparison is the degenerate version of the same situation.
    std::cerr << "  -> Verifying comparing an interface with itself returns false\n";
    EXPECT_FALSE(areSameExceptState(iface1, iface1))
        << "Comparing an interface with itself must return false (states are equal)";
}

// -----------------------------------------------------------------------------
// Test 3: Different interface names must break the match.
// -----------------------------------------------------------------------------
void test_assxs_false_when_names_differ()
{
    std::cerr << "\n[TEST] test_assxs_false_when_names_differ\n"
              << "  Function:    areSameExceptState()\n"
              << "  Scenario:    the interface names differ (\"a\" vs \"b\"), everything\n"
              << "               else is a valid 'same except state' configuration.\n"
              << "  Pass:        the function must return false.\n";

    ParsedMol::IfaceInfo ifaceA = assxs_make_iface("a", 'U', false, 0);
    ParsedMol::IfaceInfo ifaceB = assxs_make_iface("b", 'P', false, 0);
    assxs_describe("iface1", ifaceA);
    assxs_describe("iface2", ifaceB);

    EXPECT_FALSE(areSameExceptState(ifaceA, ifaceB))
        << "Interfaces with different names can never be 'same except state'";

    // The string comparison is case sensitive: "A" is not "a".
    std::cerr << "  -> Verifying the name comparison is case sensitive\n";
    ParsedMol::IfaceInfo ifaceUpper = assxs_make_iface("A", 'P', false, 0);
    EXPECT_FALSE(areSameExceptState(ifaceA, ifaceUpper))
        << "std::string comparison is case sensitive, so \"a\" != \"A\"";
}

// -----------------------------------------------------------------------------
// Test 4: Different binding status must break the match.
// -----------------------------------------------------------------------------
void test_assxs_false_when_isBound_differs()
{
    std::cerr << "\n[TEST] test_assxs_false_when_isBound_differs\n"
              << "  Function:    areSameExceptState()\n"
              << "  Scenario:    same name / speciesIndex and differing states, but\n"
              << "               one interface is bound and the other is free.\n"
              << "  Pass:        the function must return false.\n";

    ParsedMol::IfaceInfo freeIface = assxs_make_iface("a", 'U', false, 0);
    ParsedMol::IfaceInfo boundIface = assxs_make_iface("a", 'P', true, 0);
    assxs_describe("iface1", freeIface);
    assxs_describe("iface2", boundIface);

    EXPECT_FALSE(areSameExceptState(freeIface, boundIface))
        << "A bound and an unbound interface must not be reported as 'same except state'";

    // Both bound with different states is a legitimate match again.
    std::cerr << "  -> Control: both bound, differing states => true\n";
    ParsedMol::IfaceInfo boundOther = assxs_make_iface("a", 'U', true, 0);
    EXPECT_TRUE(areSameExceptState(boundIface, boundOther))
        << "Two bound interfaces differing only by state must return true";
}

// -----------------------------------------------------------------------------
// Test 5: Different speciesIndex must break the match.
// -----------------------------------------------------------------------------
void test_assxs_false_when_speciesIndex_differs()
{
    std::cerr << "\n[TEST] test_assxs_false_when_speciesIndex_differs\n"
              << "  Function:    areSameExceptState()\n"
              << "  Scenario:    identical name / isBound and differing states, but\n"
              << "               the interfaces live on different reaction species.\n"
              << "  Pass:        the function must return false.\n";

    ParsedMol::IfaceInfo species0 = assxs_make_iface("a", 'U', false, 0);
    ParsedMol::IfaceInfo species1 = assxs_make_iface("a", 'P', false, 1);
    assxs_describe("iface1", species0);
    assxs_describe("iface2", species1);

    EXPECT_FALSE(areSameExceptState(species0, species1))
        << "Interfaces on different species (speciesIndex 0 vs 1) must return false";

    // Default-constructed speciesIndex is -1; two such interfaces still match
    // each other because -1 == -1.
    std::cerr << "  -> Control: both interfaces keep the default speciesIndex (-1)\n";
    ParsedMol::IfaceInfo defaultA {};
    ParsedMol::IfaceInfo defaultB {};
    defaultA.ifaceName = "a";
    defaultB.ifaceName = "a";
    defaultA.state = 'U';
    defaultB.state = 'P';
    // isBound defaults to false and speciesIndex defaults to -1 for both.
    EXPECT_TRUE(areSameExceptState(defaultA, defaultB))
        << "Two default speciesIndex values (-1) compare equal, so the call must return true";
}

// -----------------------------------------------------------------------------
// Test 6: The "no state" character '\0' participates like any other value.
// -----------------------------------------------------------------------------
void test_assxs_handles_null_state_character()
{
    std::cerr << "\n[TEST] test_assxs_handles_null_state_character\n"
              << "  Function:    areSameExceptState()\n"
              << "  Scenario:    one interface has no declared state ('\\0') while the\n"
              << "               other declares an explicit state.\n"
              << "  Pass:        '\\0' != 'P' so the call returns true; two '\\0'\n"
              << "               states compare equal so the call returns false.\n";

    ParsedMol::IfaceInfo noState = assxs_make_iface("a", '\0', false, 0);
    ParsedMol::IfaceInfo hasState = assxs_make_iface("a", 'P', false, 0);
    assxs_describe("iface1", noState);
    assxs_describe("iface2", hasState);

    EXPECT_TRUE(areSameExceptState(noState, hasState))
        << "'\\0' and 'P' are different characters, so the interfaces differ only by state";

    std::cerr << "  -> Two stateless interfaces have equal states => false\n";
    ParsedMol::IfaceInfo noStateAgain = assxs_make_iface("a", '\0', false, 0);
    EXPECT_FALSE(areSameExceptState(noState, noStateAgain))
        << "Two interfaces that both lack a state have equal states, so the call must return false";
}

// -----------------------------------------------------------------------------
// Test 7: Members the function does not look at must not influence the result.
// -----------------------------------------------------------------------------
void test_assxs_ignores_unrelated_members()
{
    std::cerr << "\n[TEST] test_assxs_ignores_unrelated_members\n"
              << "  Function:    areSameExceptState()\n"
              << "  Scenario:    the two interfaces differ wildly in molTypeIndex,\n"
              << "               absIndex, relIndex, bondIndex, indexFound,\n"
              << "               ifaceAndStateName and ifaceRxnStatus, but agree on\n"
              << "               name / isBound / speciesIndex and differ in state.\n"
              << "  Pass:        the function only inspects the four relevant members,\n"
              << "               so it must still return true.\n";

    ParsedMol::IfaceInfo iface1 = assxs_make_iface("a", 'U', true, 3);
    iface1.molTypeIndex = 7;
    iface1.absIndex = 11;
    iface1.relIndex = 2;
    iface1.bondIndex = 1;
    iface1.indexFound = true;
    iface1.ifaceRxnStatus = Involvement::interactionChange;
    iface1.ifaceAndStateName = "completely~different";

    ParsedMol::IfaceInfo iface2 = assxs_make_iface("a", 'P', true, 3);
    iface2.molTypeIndex = 99;
    iface2.absIndex = -5;
    iface2.relIndex = 42;
    iface2.bondIndex = -1;
    iface2.indexFound = false;
    iface2.ifaceRxnStatus = Involvement::ancillary;
    iface2.ifaceAndStateName = "another~name";

    assxs_describe("iface1", iface1);
    assxs_describe("iface2", iface2);
    std::cerr << "      (molTypeIndex, absIndex, relIndex, bondIndex, indexFound,\n"
              << "       ifaceRxnStatus and ifaceAndStateName all differ)\n";

    EXPECT_TRUE(areSameExceptState(iface1, iface2))
        << "Only ifaceName, isBound, state and speciesIndex are considered by the predicate";
}

// -----------------------------------------------------------------------------
// Test 8: Objects created through the real IfaceInfo constructors.
// -----------------------------------------------------------------------------
void test_assxs_with_constructor_built_ifaces()
{
    std::cerr << "\n[TEST] test_assxs_with_constructor_built_ifaces\n"
              << "  Function:    areSameExceptState()\n"
              << "  Scenario:    build the interfaces with the production constructor\n"
              << "               IfaceInfo(name, state, isBound, Involvement, speciesIndex)\n"
              << "               instead of assigning members by hand.\n"
              << "  Pass:        same expectations as the hand-built objects.\n";

    // Constructor #1: (ifaceName, state, isBound, ifaceRxnStatus, speciesIndex).
    // It also fills ifaceAndStateName as "name~state", which the predicate ignores.
    ParsedMol::IfaceInfo ctorUnphos { "a", 'U', false, Involvement::possible, 0 };
    ParsedMol::IfaceInfo ctorPhos { "a", 'P', false, Involvement::stateChange, 0 };
    assxs_describe("iface1", ctorUnphos);
    assxs_describe("iface2", ctorPhos);
    std::cerr << "      derived ifaceAndStateName: \"" << ctorUnphos.ifaceAndStateName << "\" vs \""
              << ctorPhos.ifaceAndStateName << "\"\n";

    EXPECT_TRUE(areSameExceptState(ctorUnphos, ctorPhos))
        << "Constructor-built interfaces differing only by state must return true";

    // Same constructor, but now the species indices differ -> false.
    std::cerr << "  -> Same construction but with speciesIndex 0 vs 1\n";
    ParsedMol::IfaceInfo ctorOtherSpecies { "a", 'P', false, Involvement::stateChange, 1 };
    EXPECT_FALSE(areSameExceptState(ctorUnphos, ctorOtherSpecies))
        << "Constructor-built interfaces on different species must return false";

    // Same constructor, but with different interface names -> false.
    std::cerr << "  -> Same construction but with names \"a\" vs \"b\"\n";
    ParsedMol::IfaceInfo ctorOtherName { "b", 'P', false, Involvement::stateChange, 0 };
    EXPECT_FALSE(areSameExceptState(ctorUnphos, ctorOtherName))
        << "Constructor-built interfaces with different names must return false";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is run inside its own TEST so that a
// failure in one scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(AreSameExceptStateTest, TrueWhenOnlyStateDiffers) { test_assxs_true_when_only_state_differs(); }
TEST(AreSameExceptStateTest, FalseWhenStatesAreIdentical) { test_assxs_false_when_states_are_identical(); }
TEST(AreSameExceptStateTest, FalseWhenNamesDiffer) { test_assxs_false_when_names_differ(); }
TEST(AreSameExceptStateTest, FalseWhenIsBoundDiffers) { test_assxs_false_when_isBound_differs(); }
TEST(AreSameExceptStateTest, FalseWhenSpeciesIndexDiffers) { test_assxs_false_when_speciesIndex_differs(); }
TEST(AreSameExceptStateTest, HandlesNullStateCharacter) { test_assxs_handles_null_state_character(); }
TEST(AreSameExceptStateTest, IgnoresUnrelatedMembers) { test_assxs_ignores_unrelated_members(); }
TEST(AreSameExceptStateTest, WithConstructorBuiltIfaces) { test_assxs_with_constructor_built_ifaces(); }