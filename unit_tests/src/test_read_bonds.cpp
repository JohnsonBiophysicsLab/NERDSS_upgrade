/*! \file test_read_bonds.cpp
 *
 * ### Unit test for ../src/parser/read_bonds.cpp
 *
 * Exercises the single function defined in that translation unit:
 *
 *     void read_bonds(int numBonds, std::ifstream& molFile, MolTemplate& molTemplate)
 *
 * `read_bonds()` reads `numBonds` lines of the form "<atom1> <atom2>" from an
 * already-open molecule file and converts each pair into a PSF-style bond index
 * pair stored in `molTemplate.bondList`.  The index convention implemented by
 * the source is:
 *
 *   - the string "com"          -> index 0
 *   - interfaceList[i].name     -> index i + 1
 *
 * Both atom names are lower-cased before matching, anything trailing on the
 * line is discarded (ignore up to '\n'), and the resulting pair is sorted
 * ascending before being appended.  If **either** atom of **any** line cannot
 * be resolved, the function prints a warning, wipes the *entire* bondList
 * (including bonds added before the call) and returns immediately.
 *
 * The tests below build tiny temporary molecule files on disk (the function
 * signature demands a std::ifstream, so a stringstream cannot be substituted)
 * and check each of those documented behaviours.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers.  All names carry an "rb_" (read_bonds) prefix so they cannot
// collide with helpers defined by other test translation units in the suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Write `content` to a file called `fileName` and return that name.
 *
 * read_bonds() insists on a std::ifstream, so the test data must live in a
 * real file on disk for the duration of the test.
 */
std::string rb_write_temp_file(const std::string& fileName, const std::string& content)
{
    std::ofstream out(fileName);
    out << content;
    out.close();
    std::cerr << "    (wrote temp file \"" << fileName << "\" containing:)\n";
    std::cerr << "      ---8<---\n";
    std::cerr << content;
    std::cerr << "      --->8---\n";
    return fileName;
}

/*! \brief Build a fully-initialised MolTemplate with the requested interfaces.
 *
 * Only `molName` (used for the warning message) and `interfaceList` (used for
 * the name lookup) matter to read_bonds(), but everything touched by the
 * function is filled in explicitly so the object is never read uninitialised.
 */
MolTemplate rb_make_template(const std::string& molName, const std::vector<std::string>& ifaceNames)
{
    MolTemplate molTemplate;
    molTemplate.molName = molName;
    molTemplate.comCoord = Coord { 0.0, 0.0, 0.0 };
    molTemplate.molTypeIndex = 0;
    molTemplate.bondList.clear();

    for (unsigned i { 0 }; i < ifaceNames.size(); ++i) {
        Interface iface;
        iface.name = ifaceNames[i];
        iface.index = static_cast<int>(i);
        iface.iCoord = Coord { static_cast<double>(i) + 1.0, 0.0, 0.0 };
        molTemplate.interfaceList.push_back(iface);
    }
    return molTemplate;
}

/*! \brief Pretty-print a bondList so failures are easy to read on the console. */
void rb_dump_bonds(const MolTemplate& molTemplate)
{
    std::cerr << "    resulting bondList (size " << molTemplate.bondList.size() << "):";
    for (const auto& bond : molTemplate.bondList)
        std::cerr << " [" << bond[0] << "," << bond[1] << "]";
    std::cerr << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: The common case -- every bond connects the COM to an interface.
//         "com" must map to 0 and interface i must map to i + 1.
// -----------------------------------------------------------------------------
void test_rb_com_to_interface_bonds()
{
    std::cerr << "\n[TEST] test_rb_com_to_interface_bonds\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario:      three 'com <iface>' bond lines are parsed.\n"
              << "  Pass criteria: bondList == {{0,1},{0,2},{0,3}} in file order.\n";

    // Molecule with three interfaces -> expected indices 1, 2 and 3.
    MolTemplate molTemplate = rb_make_template("testmol", { "a", "b", "c" });

    const std::string fileName = rb_write_temp_file(
        "test_read_bonds_com_iface.tmp", "com a\ncom b\ncom c\n");

    std::ifstream molFile(fileName);
    ASSERT_TRUE(molFile.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  Calling read_bonds(3, ...)\n";
    read_bonds(3, molFile, molTemplate);
    molFile.close();
    rb_dump_bonds(molTemplate);

    // Three valid lines -> three bonds appended, none discarded.
    EXPECT_EQ(molTemplate.bondList.size(), 3u)
        << "three valid bond lines should produce three bonds";

    if (molTemplate.bondList.size() == 3u) {
        // "com" -> 0, interface index i -> i + 1, and each pair is sorted.
        EXPECT_EQ(molTemplate.bondList[0][0], 0) << "bond 0 should start at the COM (index 0)";
        EXPECT_EQ(molTemplate.bondList[0][1], 1) << "interface 'a' should map to index 1";
        EXPECT_EQ(molTemplate.bondList[1][0], 0) << "bond 1 should start at the COM (index 0)";
        EXPECT_EQ(molTemplate.bondList[1][1], 2) << "interface 'b' should map to index 2";
        EXPECT_EQ(molTemplate.bondList[2][0], 0) << "bond 2 should start at the COM (index 0)";
        EXPECT_EQ(molTemplate.bondList[2][1], 3) << "interface 'c' should map to index 3";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 2: A bond between two interfaces (no COM involved), and the guarantee
//         that the stored pair is sorted ascending regardless of file order.
// -----------------------------------------------------------------------------
void test_rb_interface_to_interface_and_sorting()
{
    std::cerr << "\n[TEST] test_rb_interface_to_interface_and_sorting\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario:      'a b' then the reversed pair 'c a'.\n"
              << "  Pass criteria: both stored pairs are sorted ascending,\n"
              << "                 i.e. {1,2} and {1,3}.\n";

    MolTemplate molTemplate = rb_make_template("testmol", { "a", "b", "c" });

    // The second line lists the higher index first; read_bonds must sort it.
    const std::string fileName = rb_write_temp_file(
        "test_read_bonds_iface_sort.tmp", "a b\nc a\n");

    std::ifstream molFile(fileName);
    ASSERT_TRUE(molFile.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  Calling read_bonds(2, ...)\n";
    read_bonds(2, molFile, molTemplate);
    molFile.close();
    rb_dump_bonds(molTemplate);

    EXPECT_EQ(molTemplate.bondList.size(), 2u)
        << "two valid interface-interface lines should produce two bonds";

    if (molTemplate.bondList.size() == 2u) {
        // "a b" -> {1,2}, already ascending.
        EXPECT_EQ(molTemplate.bondList[0][0], 1) << "'a' should map to index 1";
        EXPECT_EQ(molTemplate.bondList[0][1], 2) << "'b' should map to index 2";
        // "c a" -> {3,1} before the std::sort, {1,3} after it.
        EXPECT_EQ(molTemplate.bondList[1][0], 1) << "reversed pair must be sorted: low index first";
        EXPECT_EQ(molTemplate.bondList[1][1], 3) << "reversed pair must be sorted: high index second";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 3: Atom names are case-insensitive because read_bonds() lower-cases
//         both tokens before comparing them.
// -----------------------------------------------------------------------------
void test_rb_case_insensitive_names()
{
    std::cerr << "\n[TEST] test_rb_case_insensitive_names\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario:      tokens written as 'COM' / 'A' / 'B'.\n"
              << "  Pass criteria: they still resolve to indices 0/1/2.\n";

    MolTemplate molTemplate = rb_make_template("testmol", { "a", "b" });

    // Mixed / upper case on purpose - the source applies ::towlower to both.
    const std::string fileName = rb_write_temp_file(
        "test_read_bonds_case.tmp", "COM A\nB A\n");

    std::ifstream molFile(fileName);
    ASSERT_TRUE(molFile.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  Calling read_bonds(2, ...)\n";
    read_bonds(2, molFile, molTemplate);
    molFile.close();
    rb_dump_bonds(molTemplate);

    EXPECT_EQ(molTemplate.bondList.size(), 2u)
        << "upper-case atom names must still be recognised";

    if (molTemplate.bondList.size() == 2u) {
        EXPECT_EQ(molTemplate.bondList[0][0], 0) << "'COM' should be lower-cased and map to 0";
        EXPECT_EQ(molTemplate.bondList[0][1], 1) << "'A' should be lower-cased and map to 1";
        EXPECT_EQ(molTemplate.bondList[1][0], 1) << "'B A' should sort to {1,2}: first element 1";
        EXPECT_EQ(molTemplate.bondList[1][1], 2) << "'B A' should sort to {1,2}: second element 2";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 4: Anything after the two atom names on a line is ignored, because the
//         function skips to the end of line after reading the second token.
// -----------------------------------------------------------------------------
void test_rb_trailing_text_ignored()
{
    std::cerr << "\n[TEST] test_rb_trailing_text_ignored\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario:      first line has trailing junk after 'com a'.\n"
              << "  Pass criteria: junk is skipped and the *next* line still\n"
              << "                 parses as a normal bond ('a b' -> {1,2}).\n";

    MolTemplate molTemplate = rb_make_template("testmol", { "a", "b" });

    // "# comment here" must be swallowed by the ignore-to-newline call.
    const std::string fileName = rb_write_temp_file(
        "test_read_bonds_trailing.tmp", "com a  # comment here\na b\n");

    std::ifstream molFile(fileName);
    ASSERT_TRUE(molFile.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  Calling read_bonds(2, ...)\n";
    read_bonds(2, molFile, molTemplate);
    molFile.close();
    rb_dump_bonds(molTemplate);

    EXPECT_EQ(molTemplate.bondList.size(), 2u)
        << "trailing text must not break parsing of the following line";

    if (molTemplate.bondList.size() == 2u) {
        EXPECT_EQ(molTemplate.bondList[0][0], 0) << "line 1 should still yield the COM index 0";
        EXPECT_EQ(molTemplate.bondList[0][1], 1) << "line 1 should still yield interface index 1";
        EXPECT_EQ(molTemplate.bondList[1][0], 1) << "line 2 should yield {1,2}: first element";
        EXPECT_EQ(molTemplate.bondList[1][1], 2) << "line 2 should yield {1,2}: second element";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 5: An unknown atom name aborts the whole bond block: the bondList is
//         cleared (including bonds parsed earlier in the same call) and the
//         function returns without touching later lines.
// -----------------------------------------------------------------------------
void test_rb_unknown_atom_clears_and_aborts()
{
    std::cerr << "\n[TEST] test_rb_unknown_atom_clears_and_aborts\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario:      line 2 names an interface that does not exist,\n"
              << "                 line 3 would otherwise be a valid bond.\n"
              << "  Pass criteria: bondList ends up EMPTY (the good line 1 bond\n"
              << "                 is discarded and line 3 is never read).\n";

    MolTemplate molTemplate = rb_make_template("testmol", { "a", "b" });

    // 'zzz' is not an interface of this template and is not "com".
    const std::string fileName = rb_write_temp_file(
        "test_read_bonds_unknown.tmp", "com a\ncom zzz\ncom b\n");

    std::ifstream molFile(fileName);
    ASSERT_TRUE(molFile.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  Calling read_bonds(3, ...) -- a warning message is expected below\n";
    read_bonds(3, molFile, molTemplate);
    molFile.close();
    rb_dump_bonds(molTemplate);

    EXPECT_TRUE(molTemplate.bondList.empty())
        << "an unresolvable atom name must wipe the entire bondList";
    EXPECT_EQ(molTemplate.bondList.size(), 0u)
        << "no bonds at all should survive an invalid bond line";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 6: Bonds already present in the template before the call are appended
//         to on success, but destroyed on failure.  Both halves are checked.
// -----------------------------------------------------------------------------
void test_rb_preexisting_bonds()
{
    std::cerr << "\n[TEST] test_rb_preexisting_bonds\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario (a):  template already holds one bond, file is valid.\n"
              << "  Pass criteria: the old bond survives and the new one is appended.\n"
              << "  Scenario (b):  same starting state, file is invalid.\n"
              << "  Pass criteria: the pre-existing bond is cleared as well.\n";

    // --- Scenario (a): valid input appends to what is already there. --------
    MolTemplate goodTemplate = rb_make_template("testmol", { "a", "b" });
    goodTemplate.bondList.push_back(std::array<int, 2> { { 7, 8 } }); // sentinel bond

    const std::string goodFile = rb_write_temp_file(
        "test_read_bonds_append.tmp", "com b\n");

    std::ifstream goodStream(goodFile);
    ASSERT_TRUE(goodStream.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  (a) Calling read_bonds(1, ...) on a template with 1 existing bond\n";
    read_bonds(1, goodStream, goodTemplate);
    goodStream.close();
    rb_dump_bonds(goodTemplate);

    EXPECT_EQ(goodTemplate.bondList.size(), 2u)
        << "the new bond should be appended, not replace the existing list";
    if (goodTemplate.bondList.size() == 2u) {
        EXPECT_EQ(goodTemplate.bondList[0][0], 7) << "pre-existing bond must be untouched";
        EXPECT_EQ(goodTemplate.bondList[0][1], 8) << "pre-existing bond must be untouched";
        EXPECT_EQ(goodTemplate.bondList[1][0], 0) << "appended bond should be {0,2}: first element";
        EXPECT_EQ(goodTemplate.bondList[1][1], 2) << "appended bond should be {0,2}: second element";
    }
    std::remove(goodFile.c_str());

    // --- Scenario (b): invalid input clears everything, old bonds included. --
    MolTemplate badTemplate = rb_make_template("testmol", { "a", "b" });
    badTemplate.bondList.push_back(std::array<int, 2> { { 7, 8 } }); // same sentinel

    const std::string badFile = rb_write_temp_file(
        "test_read_bonds_clear.tmp", "nosuchiface b\n");

    std::ifstream badStream(badFile);
    ASSERT_TRUE(badStream.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  (b) Calling read_bonds(1, ...) with an invalid atom name\n";
    read_bonds(1, badStream, badTemplate);
    badStream.close();
    rb_dump_bonds(badTemplate);

    EXPECT_TRUE(badTemplate.bondList.empty())
        << "an invalid bond line clears the list including pre-existing bonds";
    std::remove(badFile.c_str());
}

// -----------------------------------------------------------------------------
// Test 7: numBonds == 0 is a no-op: nothing is read from the stream and the
//         bondList is left exactly as it was.
// -----------------------------------------------------------------------------
void test_rb_zero_bonds_is_noop()
{
    std::cerr << "\n[TEST] test_rb_zero_bonds_is_noop\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario:      numBonds == 0 with a non-empty file.\n"
              << "  Pass criteria: bondList unchanged AND the stream is still\n"
              << "                 positioned at the very first token.\n";

    MolTemplate molTemplate = rb_make_template("testmol", { "a" });

    const std::string fileName = rb_write_temp_file(
        "test_read_bonds_zero.tmp", "com a\n");

    std::ifstream molFile(fileName);
    ASSERT_TRUE(molFile.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  Calling read_bonds(0, ...)\n";
    read_bonds(0, molFile, molTemplate);
    rb_dump_bonds(molTemplate);

    EXPECT_TRUE(molTemplate.bondList.empty())
        << "requesting zero bonds must not add anything to the bondList";

    // Nothing should have been consumed, so the first token is still available.
    std::string firstToken;
    molFile >> firstToken;
    std::cerr << "    first token still available in the stream: \"" << firstToken << "\"\n";
    EXPECT_EQ(firstToken, std::string { "com" })
        << "read_bonds(0, ...) must not consume any input";

    molFile.close();
    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 8: Only `numBonds` lines are consumed -- the rest of the file stays
//         available to the caller (this is how parse_molFile() continues).
// -----------------------------------------------------------------------------
void test_rb_stops_after_numbonds_lines()
{
    std::cerr << "\n[TEST] test_rb_stops_after_numbonds_lines\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario:      file has 3 bond lines plus a trailing keyword,\n"
              << "                 but only 2 bonds are requested.\n"
              << "  Pass criteria: 2 bonds parsed and the 3rd line is still the\n"
              << "                 next thing readable from the stream.\n";

    MolTemplate molTemplate = rb_make_template("testmol", { "a", "b" });

    const std::string fileName = rb_write_temp_file(
        "test_read_bonds_partial.tmp", "com a\ncom b\nendkeyword\n");

    std::ifstream molFile(fileName);
    ASSERT_TRUE(molFile.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  Calling read_bonds(2, ...)\n";
    read_bonds(2, molFile, molTemplate);
    rb_dump_bonds(molTemplate);

    EXPECT_EQ(molTemplate.bondList.size(), 2u)
        << "exactly numBonds lines should be turned into bonds";

    // The parser must have stopped exactly at the start of line 3.
    std::string nextToken;
    molFile >> nextToken;
    std::cerr << "    next token left in the stream: \"" << nextToken << "\"\n";
    EXPECT_EQ(nextToken, std::string { "endkeyword" })
        << "read_bonds should leave the stream positioned after the last bond line";

    molFile.close();
    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 9: Asking for more bonds than the file contains.  The extraction fails,
//         the tokens stay empty, neither matches an atom name, so the function
//         takes its "invalid atom" branch and clears the list.
// -----------------------------------------------------------------------------
void test_rb_more_bonds_than_lines()
{
    std::cerr << "\n[TEST] test_rb_more_bonds_than_lines\n"
              << "  Source file:   src/parser/read_bonds.cpp\n"
              << "  Function:      read_bonds()\n"
              << "  Scenario:      numBonds == 2 but the file holds only 1 line.\n"
              << "  Pass criteria: the failed extraction is treated like an\n"
              << "                 invalid atom, so bondList ends up EMPTY.\n";

    MolTemplate molTemplate = rb_make_template("testmol", { "a" });

    const std::string fileName = rb_write_temp_file(
        "test_read_bonds_short.tmp", "com a\n");

    std::ifstream molFile(fileName);
    ASSERT_TRUE(molFile.is_open()) << "temporary molecule file could not be opened";

    std::cerr << "  Calling read_bonds(2, ...) -- a warning message is expected below\n";
    read_bonds(2, molFile, molTemplate);
    molFile.close();
    rb_dump_bonds(molTemplate);

    EXPECT_TRUE(molTemplate.bondList.empty())
        << "running out of input mid-block must clear the bondList";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so that a
// failure in one scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(ReadBonds, ComToInterfaceBonds) { test_rb_com_to_interface_bonds(); }
TEST(ReadBonds, InterfaceToInterfaceAndSorting) { test_rb_interface_to_interface_and_sorting(); }
TEST(ReadBonds, CaseInsensitiveNames) { test_rb_case_insensitive_names(); }
TEST(ReadBonds, TrailingTextIgnored) { test_rb_trailing_text_ignored(); }
TEST(ReadBonds, UnknownAtomClearsAndAborts) { test_rb_unknown_atom_clears_and_aborts(); }
TEST(ReadBonds, PreexistingBonds) { test_rb_preexisting_bonds(); }
TEST(ReadBonds, ZeroBondsIsNoop) { test_rb_zero_bonds_is_noop(); }
TEST(ReadBonds, StopsAfterNumBondsLines) { test_rb_stops_after_numbonds_lines(); }
TEST(ReadBonds, MoreBondsThanLines) { test_rb_more_bonds_than_lines(); }