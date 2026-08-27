/*! \file test_read_internal_coordinates.cpp
 *
 * ### Unit test for ../src/parser/read_internal_coordinates.cpp
 *
 * Function under test:
 *
 *     void read_internal_coordinates(std::ifstream& molFile,
 *                                    MolTemplate& molTemplate)
 *
 * Behaviour that the implementation actually has (and which the assertions
 * below check):
 *
 *   * It repeatedly extracts four whitespace separated tokens of the form
 *     `<name> <x> <y> <z>` from the stream.
 *   * If `<name>` lowercases to exactly "com" the values are stored in
 *     `molTemplate.comCoord`; otherwise a new `Interface` is appended to
 *     `molTemplate.interfaceList` using the ORIGINAL (not lowercased) name.
 *   * After each successfully parsed record the rest of the line is discarded
 *     (this is how trailing comments are skipped).
 *   * The loop stops at the first line that cannot be parsed as
 *     `string double double double`.  The stream is then `clear()`ed and
 *     rewound (`seekg`) to the beginning of that offending line, so the caller
 *     can keep parsing keywords.
 *   * Existing entries in `interfaceList` are kept -- new interfaces are
 *     appended, never assigned.
 *
 * Every test writes a small temporary file, opens it with an std::ifstream,
 * calls the function, and inspects both the resulting MolTemplate and the
 * stream position.  Verbose progress information is written to stderr.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

// -----------------------------------------------------------------------------
// Small local helpers (unique "ric_" prefix so they cannot collide with other
// translation units in the combined gtest binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Write `contents` to a uniquely named temporary file.
 *
 * \param[in] tag       Short tag appended to the file name so each test has its
 *                      own private file.
 * \param[in] contents  Exact text to write.
 * \return The name of the file that was created.
 */
std::string ric_write_temp_file(const std::string& tag, const std::string& contents)
{
    const std::string fileName = "test_read_internal_coordinates_" + tag + ".tmp";
    std::ofstream out(fileName.c_str());
    out << contents;
    out.close();
    return fileName;
}

/*! \brief Delete a temporary file created by ric_write_temp_file(). */
void ric_remove_temp_file(const std::string& fileName)
{
    std::remove(fileName.c_str());
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a plain coordinate block containing a COM line and two interfaces.
// -----------------------------------------------------------------------------
void test_ric_com_and_interfaces()
{
    std::cerr << "\n[TEST] test_ric_com_and_interfaces\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    block with 'COM' plus two interface lines.\n"
              << "  Pass:        comCoord holds the COM values, interfaceList\n"
              << "               holds the two interfaces in file order with\n"
              << "               their names and coordinates preserved.\n";

    // The trailing newline keeps the final ignore() from hitting EOF.
    const std::string fileName = ric_write_temp_file("basic",
        "COM 0.0 0.0 0.0\n"
        "a1 1.0 2.0 3.0\n"
        "a2 -4.5 5.25 -6.125\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    MolTemplate molTemplate;
    std::cerr << "  Calling read_internal_coordinates()... (function echoes the\n"
              << "  parsed coordinates to stdout)\n";
    read_internal_coordinates(molFile, molTemplate);

    // COM should have been captured into comCoord, not into interfaceList.
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.x, 0.0) << "COM x should be 0.0";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.y, 0.0) << "COM y should be 0.0";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.z, 0.0) << "COM z should be 0.0";

    // Exactly two interfaces should have been appended (COM is excluded).
    ASSERT_EQ(molTemplate.interfaceList.size(), 2u)
        << "Two non-COM lines should produce exactly two interfaces";

    std::cerr << "  Parsed " << molTemplate.interfaceList.size() << " interfaces.\n";

    // First interface: name and coordinates.
    EXPECT_EQ(molTemplate.interfaceList[0].name, std::string("a1"))
        << "First interface should be named 'a1'";
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.x, 1.0);
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.y, 2.0);
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.z, 3.0);

    // Second interface: negative and fractional values must survive intact.
    EXPECT_EQ(molTemplate.interfaceList[1].name, std::string("a2"))
        << "Second interface should be named 'a2'";
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[1].iCoord.x, -4.5);
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[1].iCoord.y, 5.25);
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[1].iCoord.z, -6.125);

    // The Interface(name, coord) constructor leaves index at its default -1;
    // this function does not assign relative indices.
    EXPECT_EQ(molTemplate.interfaceList[0].index, -1)
        << "read_internal_coordinates does not set Interface::index";
    EXPECT_EQ(molTemplate.interfaceList[1].index, -1)
        << "read_internal_coordinates does not set Interface::index";

    // The function calls clear() before returning, so the stream must be usable.
    EXPECT_FALSE(molFile.bad()) << "Stream should not be in a bad state";
    EXPECT_FALSE(molFile.fail()) << "clear() should have removed the failbit";

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// Test 2: the "COM" keyword is matched case-insensitively, but interface names
//         keep their original capitalisation.
// -----------------------------------------------------------------------------
void test_ric_com_keyword_is_case_insensitive()
{
    std::cerr << "\n[TEST] test_ric_com_keyword_is_case_insensitive\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    'CoM' as the centre-of-mass keyword and a mixed\n"
              << "               case interface name 'SiteA'.\n"
              << "  Pass:        'CoM' fills comCoord; 'SiteA' is stored as an\n"
              << "               interface with its original spelling.\n";

    const std::string fileName = ric_write_temp_file("case",
        "CoM 7.0 8.0 9.0\n"
        "SiteA 1.0 0.0 0.0\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    MolTemplate molTemplate;
    read_internal_coordinates(molFile, molTemplate);

    // 'CoM' lowercases to "com" and therefore sets the centre of mass.
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.x, 7.0) << "'CoM' should set comCoord.x";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.y, 8.0) << "'CoM' should set comCoord.y";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.z, 9.0) << "'CoM' should set comCoord.z";

    // Only the non-COM record becomes an interface.
    ASSERT_EQ(molTemplate.interfaceList.size(), 1u)
        << "Only 'SiteA' should have become an interface";
    EXPECT_EQ(molTemplate.interfaceList[0].name, std::string("SiteA"))
        << "The interface name must keep its original capitalisation";

    std::cerr << "  comCoord = (" << molTemplate.comCoord.x << ", "
              << molTemplate.comCoord.y << ", " << molTemplate.comCoord.z
              << "), interface name = '" << molTemplate.interfaceList[0].name << "'\n";

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// Test 3: the comparison against "com" is an exact string compare, so a name
//         that merely starts with "com" is treated as a normal interface.
// -----------------------------------------------------------------------------
void test_ric_name_starting_with_com_is_an_interface()
{
    std::cerr << "\n[TEST] test_ric_name_starting_with_com_is_an_interface\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    an interface literally named 'combo'.\n"
              << "  Pass:        'combo' is NOT mistaken for the COM keyword; it\n"
              << "               becomes a regular interface and comCoord stays\n"
              << "               at its default (0,0,0).\n";

    const std::string fileName = ric_write_temp_file("combo",
        "combo 1.0 1.0 1.0\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    MolTemplate molTemplate;
    read_internal_coordinates(molFile, molTemplate);

    // comCoord is untouched because no exact "com" record was seen.
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.x, 0.0) << "comCoord should stay default";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.y, 0.0) << "comCoord should stay default";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.z, 0.0) << "comCoord should stay default";

    ASSERT_EQ(molTemplate.interfaceList.size(), 1u)
        << "'combo' should have been stored as an interface";
    EXPECT_EQ(molTemplate.interfaceList[0].name, std::string("combo"));
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.x, 1.0);

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// Test 4: everything after the fourth token on a line is discarded, which is how
//         trailing comments are handled.
// -----------------------------------------------------------------------------
void test_ric_trailing_comments_are_ignored()
{
    std::cerr << "\n[TEST] test_ric_trailing_comments_are_ignored\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    each coordinate line carries a '#' comment.\n"
              << "  Pass:        the comments are consumed by ignore() and do not\n"
              << "               break parsing of the following lines.\n";

    const std::string fileName = ric_write_temp_file("comments",
        "COM 0.0 0.0 0.0 # centre of mass\n"
        "a1 1.0 0.0 0.0 # first binding site\n"
        "a2 0.0 1.0 0.0 # second binding site\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    MolTemplate molTemplate;
    read_internal_coordinates(molFile, molTemplate);

    // Both commented interface lines must still be parsed correctly.
    ASSERT_EQ(molTemplate.interfaceList.size(), 2u)
        << "Comments must not stop the parse loop";
    EXPECT_EQ(molTemplate.interfaceList[0].name, std::string("a1"));
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.x, 1.0);
    EXPECT_EQ(molTemplate.interfaceList[1].name, std::string("a2"));
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[1].iCoord.y, 1.0);

    std::cerr << "  Both commented lines parsed; interfaceList size = "
              << molTemplate.interfaceList.size() << '\n';

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// Test 5: parsing stops at the first non-coordinate line and the stream is
//         rewound so the caller can read that line itself.
// -----------------------------------------------------------------------------
void test_ric_rewinds_to_first_non_coordinate_line()
{
    std::cerr << "\n[TEST] test_ric_rewinds_to_first_non_coordinate_line\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    the coordinate block is followed by a keyword\n"
              << "               line ('bonds = 1') and a bond definition.\n"
              << "  Pass:        only the coordinate lines are consumed, and after\n"
              << "               the call the very next getline() returns the\n"
              << "               keyword line verbatim.\n";

    const std::string fileName = ric_write_temp_file("rewind",
        "COM 0.0 0.0 0.0\n"
        "a1 1.0 0.0 0.0\n"
        "bonds = 1\n"
        "com a1\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    MolTemplate molTemplate;
    read_internal_coordinates(molFile, molTemplate);

    // Only the single coordinate line before the keyword became an interface.
    ASSERT_EQ(molTemplate.interfaceList.size(), 1u)
        << "Parsing must stop at the 'bonds = 1' line";
    EXPECT_EQ(molTemplate.interfaceList[0].name, std::string("a1"));

    // The stream must be usable (clear() was called) and positioned at the
    // beginning of the offending line.
    EXPECT_FALSE(molFile.fail())
        << "The failbit must be cleared before the function returns";

    std::string nextLine;
    const bool gotLine = static_cast<bool>(std::getline(molFile, nextLine));
    EXPECT_TRUE(gotLine) << "There should still be content left in the stream";
    EXPECT_EQ(nextLine, std::string("bonds = 1"))
        << "The stream should be rewound to the start of the keyword line";

    std::cerr << "  Next line after the call = \"" << nextLine << "\"\n";

    // And the line after that should still be intact as well.
    std::string bondLine;
    EXPECT_TRUE(static_cast<bool>(std::getline(molFile, bondLine)))
        << "The bond definition line should still be readable";
    EXPECT_EQ(bondLine, std::string("com a1"));

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// Test 6: an "empty" coordinate block (the very first line is not a coordinate)
//         leaves the MolTemplate untouched and the stream at position zero.
// -----------------------------------------------------------------------------
void test_ric_empty_block_leaves_stream_untouched()
{
    std::cerr << "\n[TEST] test_ric_empty_block_leaves_stream_untouched\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    the stream starts with a keyword line, so nothing\n"
              << "               can be parsed as a coordinate record.\n"
              << "  Pass:        comCoord stays (0,0,0), interfaceList stays empty,\n"
              << "               and the first line is still readable afterwards.\n";

    const std::string fileName = ric_write_temp_file("empty",
        "bonds = 0\n"
        "state = a1~P\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    MolTemplate molTemplate;
    read_internal_coordinates(molFile, molTemplate);

    // Nothing should have been stored.
    EXPECT_TRUE(molTemplate.interfaceList.empty())
        << "No coordinate lines means no interfaces";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.x, 0.0) << "comCoord should stay default";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.y, 0.0) << "comCoord should stay default";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.z, 0.0) << "comCoord should stay default";

    // The stream must have been rewound to the position it had on entry (0).
    EXPECT_FALSE(molFile.fail()) << "The failbit must have been cleared";
    std::string firstLine;
    EXPECT_TRUE(static_cast<bool>(std::getline(molFile, firstLine)))
        << "The first line should still be available";
    EXPECT_EQ(firstLine, std::string("bonds = 0"))
        << "The stream should be back at the start of the first line";

    std::cerr << "  First line still readable as \"" << firstLine << "\"\n";

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// Test 7: new interfaces are appended, existing entries are preserved.
// -----------------------------------------------------------------------------
void test_ric_appends_to_existing_interface_list()
{
    std::cerr << "\n[TEST] test_ric_appends_to_existing_interface_list\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    the MolTemplate already owns one interface before\n"
              << "               the call (emplace_back is used, not assignment).\n"
              << "  Pass:        the pre-existing interface stays at index 0 and\n"
              << "               the newly parsed one is appended at index 1.\n";

    const std::string fileName = ric_write_temp_file("append",
        "COM 0.0 0.0 0.0\n"
        "new1 2.0 0.0 0.0\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    // Pre-load a template with one interface so we can check appending.
    MolTemplate molTemplate;
    molTemplate.interfaceList.emplace_back(std::string("preexisting"),
                                           Coord { -1.0, -1.0, -1.0 });
    ASSERT_EQ(molTemplate.interfaceList.size(), 1u)
        << "Sanity check on the pre-loaded interface";

    read_internal_coordinates(molFile, molTemplate);

    // One interface should have been appended -> two in total.
    ASSERT_EQ(molTemplate.interfaceList.size(), 2u)
        << "The new interface must be appended, not overwrite the old one";

    // Old entry untouched.
    EXPECT_EQ(molTemplate.interfaceList[0].name, std::string("preexisting"))
        << "The pre-existing interface must remain at index 0";
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.x, -1.0);

    // New entry appended at the back.
    EXPECT_EQ(molTemplate.interfaceList[1].name, std::string("new1"))
        << "The parsed interface must be appended at index 1";
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[1].iCoord.x, 2.0);

    std::cerr << "  interfaceList now = [" << molTemplate.interfaceList[0].name
              << ", " << molTemplate.interfaceList[1].name << "]\n";

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// Test 8: a later COM record overwrites an earlier one (comCoord is assigned,
//         not accumulated), and scientific notation parses correctly.
// -----------------------------------------------------------------------------
void test_ric_repeated_com_and_scientific_notation()
{
    std::cerr << "\n[TEST] test_ric_repeated_com_and_scientific_notation\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    two COM records plus an interface written in\n"
              << "               scientific notation.\n"
              << "  Pass:        comCoord holds the LAST COM record, and the\n"
              << "               exponent notation is converted correctly by the\n"
              << "               stream's double extraction.\n";

    const std::string fileName = ric_write_temp_file("overwrite",
        "COM 1.0 1.0 1.0\n"
        "s1 1.5e2 -2.5e-1 0.0\n"
        "COM 3.0 4.0 5.0\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    MolTemplate molTemplate;
    read_internal_coordinates(molFile, molTemplate);

    // The second COM assignment wins.
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.x, 3.0)
        << "The last COM record should overwrite the first";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.y, 4.0)
        << "The last COM record should overwrite the first";
    EXPECT_DOUBLE_EQ(molTemplate.comCoord.z, 5.0)
        << "The last COM record should overwrite the first";

    // Scientific notation must round-trip through operator>>(double).
    ASSERT_EQ(molTemplate.interfaceList.size(), 1u)
        << "Only 's1' is a non-COM record";
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.x, 150.0)
        << "1.5e2 should parse as 150.0";
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.y, -0.25)
        << "-2.5e-1 should parse as -0.25";
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[0].iCoord.z, 0.0);

    std::cerr << "  Final comCoord = (" << molTemplate.comCoord.x << ", "
              << molTemplate.comCoord.y << ", " << molTemplate.comCoord.z << ")\n";

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// Test 9: a file that consists solely of coordinate lines -- the loop runs to
//         end-of-file, and the stream is left at EOF but not in a failed state.
// -----------------------------------------------------------------------------
void test_ric_consumes_whole_file()
{
    std::cerr << "\n[TEST] test_ric_consumes_whole_file\n"
              << "  Source file: src/parser/read_internal_coordinates.cpp\n"
              << "  Function:    read_internal_coordinates()\n"
              << "  Scenario:    the file contains nothing but coordinate lines.\n"
              << "  Pass:        every record is parsed, the failbit is cleared,\n"
              << "               and no further line can be read.\n";

    const std::string fileName = ric_write_temp_file("wholefile",
        "COM 0.0 0.0 0.0\n"
        "a1 1.0 0.0 0.0\n"
        "a2 0.0 1.0 0.0\n"
        "a3 0.0 0.0 1.0\n");

    std::ifstream molFile(fileName.c_str());
    ASSERT_TRUE(molFile.is_open()) << "Could not open the temporary mol file";

    MolTemplate molTemplate;
    read_internal_coordinates(molFile, molTemplate);

    // All three non-COM records must be present, in order.
    ASSERT_EQ(molTemplate.interfaceList.size(), 3u)
        << "All three interface lines should be parsed";
    EXPECT_EQ(molTemplate.interfaceList[0].name, std::string("a1"));
    EXPECT_EQ(molTemplate.interfaceList[1].name, std::string("a2"));
    EXPECT_EQ(molTemplate.interfaceList[2].name, std::string("a3"));
    EXPECT_DOUBLE_EQ(molTemplate.interfaceList[2].iCoord.z, 1.0);

    // clear() was called, so failbit/eofbit are reset before the seek.
    EXPECT_FALSE(molFile.fail())
        << "The stream must not be left in a failed state";

    // Nothing is left to read: the stream was rewound to the end of the file.
    std::string leftOver;
    EXPECT_FALSE(static_cast<bool>(std::getline(molFile, leftOver)))
        << "The whole file should have been consumed";

    std::cerr << "  All " << molTemplate.interfaceList.size()
              << " interfaces parsed; nothing left in the stream.\n";

    molFile.close();
    ric_remove_temp_file(fileName);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each named helper runs inside its own TEST so that a
// failure in one does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(ReadInternalCoordinates, ComAndInterfaces) { test_ric_com_and_interfaces(); }
TEST(ReadInternalCoordinates, ComKeywordIsCaseInsensitive) { test_ric_com_keyword_is_case_insensitive(); }
TEST(ReadInternalCoordinates, NameStartingWithComIsAnInterface) { test_ric_name_starting_with_com_is_an_interface(); }
TEST(ReadInternalCoordinates, TrailingCommentsAreIgnored) { test_ric_trailing_comments_are_ignored(); }
TEST(ReadInternalCoordinates, RewindsToFirstNonCoordinateLine) { test_ric_rewinds_to_first_non_coordinate_line(); }
TEST(ReadInternalCoordinates, EmptyBlockLeavesStreamUntouched) { test_ric_empty_block_leaves_stream_untouched(); }
TEST(ReadInternalCoordinates, AppendsToExistingInterfaceList) { test_ric_appends_to_existing_interface_list(); }
TEST(ReadInternalCoordinates, RepeatedComAndScientificNotation) { test_ric_repeated_com_and_scientific_notation(); }
TEST(ReadInternalCoordinates, ConsumesWholeFile) { test_ric_consumes_whole_file(); }