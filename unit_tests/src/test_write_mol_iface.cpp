/*! \file test_write_mol_iface.cpp
 *
 * ### Unit test for ../src/parser/write_mol_iface.cpp
 *
 * The file under test contains exactly one free function:
 *
 * \code
 *     std::string write_mol_iface(std::string mol, std::string iface);
 * \endcode
 *
 * The function builds the BNGL-style "molecule(interface)" string used
 * throughout the parser for display and for species naming, i.e. it simply
 * concatenates:
 *
 *     mol + "(" + iface + ")"
 *
 * There is no validation, no trimming, and no case conversion performed inside
 * the function, so the tests below verify:
 *   - the exact concatenation format for ordinary input,
 *   - that empty arguments still produce the surrounding parentheses,
 *   - that whitespace, tildes/state characters, bond markers, digits and other
 *     "special" characters are passed through verbatim,
 *   - that the arguments are taken by value and the caller's strings are not
 *     modified,
 *   - that the function is deterministic (same input -> same output), and
 *   - that the resulting string length is exactly len(mol)+len(iface)+2.
 *
 * Verbose progress output is written to stderr so the reader can follow which
 * source file / function is under test and what criteria each check uses.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>

// -----------------------------------------------------------------------------
// Test 1: Ordinary molecule and interface names.
//
// Pass criteria: the output is exactly "<mol>(<iface>)".
// -----------------------------------------------------------------------------
void test_wmi_basic_formatting()
{
    std::cerr << "\n[TEST] test_wmi_basic_formatting\n"
              << "  Source file: src/parser/write_mol_iface.cpp\n"
              << "  Function:    write_mol_iface(std::string, std::string)\n"
              << "  Scenario:    ordinary lowercase molecule/interface names.\n"
              << "  Criteria:    result must equal mol + \"(\" + iface + \")\".\n";

    // A typical NERDSS molecule/interface pair, e.g. clathrin with leg 1.
    const std::string result = write_mol_iface("clathrin", "leg1");
    std::cerr << "  write_mol_iface(\"clathrin\", \"leg1\") -> \"" << result << "\"\n";
    EXPECT_EQ(result, "clathrin(leg1)")
        << "Basic concatenation must produce clathrin(leg1)";

    // A single-character pair, the smallest non-empty case.
    const std::string shortResult = write_mol_iface("a", "b");
    std::cerr << "  write_mol_iface(\"a\", \"b\") -> \"" << shortResult << "\"\n";
    EXPECT_EQ(shortResult, "a(b)")
        << "Single character names must produce a(b)";

    // Another realistic pair to double-check no hidden separators are inserted.
    const std::string apResult = write_mol_iface("ap2", "beta");
    std::cerr << "  write_mol_iface(\"ap2\", \"beta\") -> \"" << apResult << "\"\n";
    EXPECT_EQ(apResult, "ap2(beta)")
        << "Names containing digits must be concatenated verbatim";
}

// -----------------------------------------------------------------------------
// Test 2: Empty arguments.
//
// The function performs no validation, so an empty molecule name and/or an
// empty interface name still produce the surrounding parentheses.
//
// Pass criteria: "()" for two empty strings, "mol()" and "(iface)" otherwise.
// -----------------------------------------------------------------------------
void test_wmi_empty_arguments()
{
    std::cerr << "\n[TEST] test_wmi_empty_arguments\n"
              << "  Source file: src/parser/write_mol_iface.cpp\n"
              << "  Function:    write_mol_iface(std::string, std::string)\n"
              << "  Scenario:    empty molecule name, empty interface name, both.\n"
              << "  Criteria:    parentheses are always emitted; no validation.\n";

    // Both arguments empty -> just the parentheses.
    const std::string bothEmpty = write_mol_iface("", "");
    std::cerr << "  write_mol_iface(\"\", \"\") -> \"" << bothEmpty << "\"\n";
    EXPECT_EQ(bothEmpty, "()")
        << "Two empty strings must still produce the parentheses only";

    // Empty interface -> "mol()".
    const std::string emptyIface = write_mol_iface("clathrin", "");
    std::cerr << "  write_mol_iface(\"clathrin\", \"\") -> \"" << emptyIface << "\"\n";
    EXPECT_EQ(emptyIface, "clathrin()")
        << "An empty interface name must yield mol()";

    // Empty molecule -> "(iface)".
    const std::string emptyMol = write_mol_iface("", "leg1");
    std::cerr << "  write_mol_iface(\"\", \"leg1\") -> \"" << emptyMol << "\"\n";
    EXPECT_EQ(emptyMol, "(leg1)")
        << "An empty molecule name must yield (iface)";
}

// -----------------------------------------------------------------------------
// Test 3: Special characters are passed through verbatim.
//
// The parser uses '~' for interface states and '!' for bond indices; the
// function must not strip, escape or reorder any of these.
//
// Pass criteria: every character of the inputs survives in order.
// -----------------------------------------------------------------------------
void test_wmi_special_characters_preserved()
{
    std::cerr << "\n[TEST] test_wmi_special_characters_preserved\n"
              << "  Source file: src/parser/write_mol_iface.cpp\n"
              << "  Function:    write_mol_iface(std::string, std::string)\n"
              << "  Scenario:    state ('~'), bond ('!'), nested parens, spaces.\n"
              << "  Criteria:    all characters are copied through unchanged.\n";

    // Interface with a state identifier, as written in BNGL (iface~State).
    const std::string stateResult = write_mol_iface("pip2", "head~P");
    std::cerr << "  write_mol_iface(\"pip2\", \"head~P\") -> \"" << stateResult << "\"\n";
    EXPECT_EQ(stateResult, "pip2(head~P)")
        << "The '~' state separator must be preserved";

    // Interface with an explicit bond index (iface!1).
    const std::string bondResult = write_mol_iface("ap2", "beta!1");
    std::cerr << "  write_mol_iface(\"ap2\", \"beta!1\") -> \"" << bondResult << "\"\n";
    EXPECT_EQ(bondResult, "ap2(beta!1)")
        << "The '!' bond marker and index must be preserved";

    // Whitespace is not trimmed by the function.
    const std::string spaceResult = write_mol_iface(" mol ", " iface ");
    std::cerr << "  write_mol_iface(\" mol \", \" iface \") -> \"" << spaceResult << "\"\n";
    EXPECT_EQ(spaceResult, " mol ( iface )")
        << "Leading/trailing whitespace must NOT be trimmed";

    // Characters that already look like the delimiters are not escaped.
    const std::string nestedResult = write_mol_iface("a(b", "c)d");
    std::cerr << "  write_mol_iface(\"a(b\", \"c)d\") -> \"" << nestedResult << "\"\n";
    EXPECT_EQ(nestedResult, "a(b(c)d)")
        << "Parentheses inside the arguments must be copied verbatim";
}

// -----------------------------------------------------------------------------
// Test 4: Arguments are taken by value, so the caller's strings are untouched.
//
// Pass criteria: after the call the original std::string objects still hold
// their original contents.
// -----------------------------------------------------------------------------
void test_wmi_arguments_not_modified()
{
    std::cerr << "\n[TEST] test_wmi_arguments_not_modified\n"
              << "  Source file: src/parser/write_mol_iface.cpp\n"
              << "  Function:    write_mol_iface(std::string, std::string)\n"
              << "  Scenario:    call with named lvalue strings.\n"
              << "  Criteria:    the caller's strings are unchanged (pass by value).\n";

    std::string mol { "clathrin" };
    std::string iface { "leg1" };

    // Keep copies of the originals so we can compare afterwards.
    const std::string molBefore { mol };
    const std::string ifaceBefore { iface };

    const std::string result = write_mol_iface(mol, iface);
    std::cerr << "  Result:      \"" << result << "\"\n";
    std::cerr << "  mol after:   \"" << mol << "\" (was \"" << molBefore << "\")\n";
    std::cerr << "  iface after: \"" << iface << "\" (was \"" << ifaceBefore << "\")\n";

    EXPECT_EQ(result, "clathrin(leg1)")
        << "The returned string must still be correctly formatted";
    EXPECT_EQ(mol, molBefore)
        << "The molecule argument must not be modified by the call";
    EXPECT_EQ(iface, ifaceBefore)
        << "The interface argument must not be modified by the call";
}

// -----------------------------------------------------------------------------
// Test 5: Determinism and length invariant.
//
// The function is pure, so repeated calls with identical input must produce
// identical output, and the length is always len(mol) + len(iface) + 2 for the
// two parentheses.
// -----------------------------------------------------------------------------
void test_wmi_deterministic_and_length()
{
    std::cerr << "\n[TEST] test_wmi_deterministic_and_length\n"
              << "  Source file: src/parser/write_mol_iface.cpp\n"
              << "  Function:    write_mol_iface(std::string, std::string)\n"
              << "  Scenario:    repeat calls and check the output size formula.\n"
              << "  Criteria:    identical output every call; size == n1+n2+2.\n";

    const std::string mol { "someLongMoleculeName" };
    const std::string iface { "someLongInterfaceName" };

    // Call the function twice with exactly the same input.
    const std::string first = write_mol_iface(mol, iface);
    const std::string second = write_mol_iface(mol, iface);

    std::cerr << "  First call:  \"" << first << "\"\n";
    std::cerr << "  Second call: \"" << second << "\"\n";
    EXPECT_EQ(first, second)
        << "The function is pure; repeated calls must give identical results";

    // The only characters added are the two parentheses.
    const size_t expectedLength = mol.size() + iface.size() + 2;
    std::cerr << "  Length: " << first.size() << " (expected " << expectedLength << ")\n";
    EXPECT_EQ(first.size(), expectedLength)
        << "Output length must be len(mol) + len(iface) + 2";

    // Sanity-check the exact positions of the delimiters.
    EXPECT_EQ(first.substr(0, mol.size()), mol)
        << "The output must start with the molecule name";
    EXPECT_EQ(first[mol.size()], '(')
        << "An opening parenthesis must directly follow the molecule name";
    EXPECT_EQ(first.back(), ')')
        << "The output must end with a closing parenthesis";
    EXPECT_EQ(first.substr(mol.size() + 1, iface.size()), iface)
        << "The interface name must sit between the parentheses";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* helper is invoked from its own TEST so
// that the framework reports individual results and every case still executes
// even if an earlier one records a failure (all checks use non-fatal EXPECT_*).
// -----------------------------------------------------------------------------
TEST(WriteMolIfaceTest, BasicFormatting) { test_wmi_basic_formatting(); }
TEST(WriteMolIfaceTest, EmptyArguments) { test_wmi_empty_arguments(); }
TEST(WriteMolIfaceTest, SpecialCharactersPreserved) { test_wmi_special_characters_preserved(); }
TEST(WriteMolIfaceTest, ArgumentsNotModified) { test_wmi_arguments_not_modified(); }
TEST(WriteMolIfaceTest, DeterministicAndLength) { test_wmi_deterministic_and_length(); }