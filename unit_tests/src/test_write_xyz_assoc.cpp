/*! \file test_write_xyz_assoc.cpp
 *
 * ### Unit test for src/io/write_xyz_assoc.cpp
 *
 * Function under test:
 *
 *     void write_xyz_assoc(std::string filename,
 *                          const Complex& reactCom1,
 *                          const Complex& reactCom2,
 *                          const std::vector<Molecule>& moleculeList);
 *
 * The routine is a debugging helper that dumps the *temporary* association
 * coordinates (`Molecule::tmpComCoord` and `Molecule::tmpICoords`) of the two
 * complexes involved in an association event to an XYZ-style text file:
 *
 *     <totAtoms>
 *     mol output final
 *     A <com coord of complex-1 member>
 *     A <each tmp interface coord of that member>
 *     ...
 *     B <com coord of complex-2 member>
 *     B <each tmp interface coord of that member>
 *     ...
 *
 * Notable implementation detail that the tests below deliberately pin down:
 * the atom counter on the first line is computed from
 * `Molecule::interfaceList.size()`, while the actual coordinate lines are
 * emitted from `Molecule::tmpICoords`.  When those two containers have
 * different sizes (which happens if the temporary association coordinates were
 * never populated) the header count and the number of written lines disagree.
 *
 * Because the exact text formatting of `operator<<(std::ostream&, const Coord&)`
 * is an implementation detail, the tests parse every floating point number out
 * of each line instead of doing brittle string comparisons.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Helper: build a minimal Molecule carrying temporary association coordinates.
//
// \param[in] com            value assigned to Molecule::tmpComCoord
// \param[in] tmpIfaceCoords values assigned to Molecule::tmpICoords (these are
//                           the coordinates that actually get written out)
// \param[in] numIfaceSlots  how many entries Molecule::interfaceList gets; the
//                           function under test uses *this* size for the atom
//                           count on the first line of the file
// -----------------------------------------------------------------------------
Molecule wxa_make_molecule(const Coord& com, const std::vector<Coord>& tmpIfaceCoords,
    std::size_t numIfaceSlots)
{
    Molecule mol;
    mol.tmpComCoord = com;
    mol.tmpICoords = tmpIfaceCoords;
    // interfaceList only needs the correct *size* for this test; contents are
    // irrelevant because write_xyz_assoc never dereferences them.
    mol.interfaceList.resize(numIfaceSlots);
    return mol;
}

// -----------------------------------------------------------------------------
// Helper: build a Complex whose memberList points at the given molecule indices.
// -----------------------------------------------------------------------------
Complex wxa_make_complex(const std::vector<int>& members)
{
    Complex com;
    com.memberList = members;
    return com;
}

// -----------------------------------------------------------------------------
// Helper: read a whole text file into a vector of lines.
// Returns an empty vector if the file could not be opened.
// -----------------------------------------------------------------------------
std::vector<std::string> wxa_read_lines(const std::string& filename)
{
    std::vector<std::string> lines;
    std::ifstream in(filename);
    if (!in.is_open())
        return lines;

    std::string line;
    while (std::getline(in, line))
        lines.push_back(line);
    return lines;
}

// -----------------------------------------------------------------------------
// Helper: pull every floating point number out of a line, treating anything
// that is not part of a number (spaces, commas, brackets, the leading atom
// label, ...) as a separator. This keeps the tests independent of the exact
// Coord streaming format.
// -----------------------------------------------------------------------------
std::vector<double> wxa_extract_doubles(const std::string& line)
{
    std::vector<double> values;
    std::size_t i = 0;
    while (i < line.size()) {
        const char c = line[i];
        const bool startsNumber = std::isdigit(static_cast<unsigned char>(c))
            || ((c == '-' || c == '+' || c == '.') && (i + 1 < line.size())
                && std::isdigit(static_cast<unsigned char>(line[i + 1])));
        if (startsNumber) {
            std::size_t consumed = 0;
            try {
                const double v = std::stod(line.substr(i), &consumed);
                values.push_back(v);
                i += (consumed > 0) ? consumed : 1;
            } catch (...) {
                ++i;
            }
        } else {
            ++i;
        }
    }
    return values;
}

// -----------------------------------------------------------------------------
// Helper: count how many lines begin with the given atom label character.
// -----------------------------------------------------------------------------
int wxa_count_label(const std::vector<std::string>& lines, char label)
{
    int count = 0;
    for (std::size_t i = 2; i < lines.size(); ++i) { // skip the two header lines
        if (!lines[i].empty() && lines[i][0] == label)
            ++count;
    }
    return count;
}

// -----------------------------------------------------------------------------
// Helper: compare the three numbers found on a coordinate line with a Coord.
// Emits verbose output and non-fatal expectations only.
// -----------------------------------------------------------------------------
void wxa_expect_line_matches_coord(const std::string& line, const Coord& expected,
    const std::string& what)
{
    const std::vector<double> nums = wxa_extract_doubles(line);
    std::cerr << "      " << what << " line = \"" << line << "\" -> parsed "
              << nums.size() << " number(s)\n";

    EXPECT_EQ(nums.size(), 3u)
        << what << ": expected exactly three numbers (x y z) on the line";
    if (nums.size() == 3u) {
        EXPECT_NEAR(nums[0], expected.x, 1e-4) << what << ": x mismatch";
        EXPECT_NEAR(nums[1], expected.y, 1e-4) << what << ": y mismatch";
        EXPECT_NEAR(nums[2], expected.z, 1e-4) << what << ": z mismatch";
    }
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Nominal case - both complexes contain molecules whose tmpICoords and
//         interfaceList sizes agree.
//
// Pass criteria:
//   * the file exists and has 2 header lines + one line per written atom
//   * line 1 == total atom count (members + their interfaces)
//   * line 2 == "mol output final"
//   * complex 1 atoms are labelled 'A', complex 2 atoms labelled 'B'
//   * the coordinates written are the *temporary* association coordinates
// -----------------------------------------------------------------------------
void test_wxa_nominal_two_complexes()
{
    std::cerr << "\n[TEST] test_wxa_nominal_two_complexes\n"
              << "  Source file: src/io/write_xyz_assoc.cpp\n"
              << "  Function:    write_xyz_assoc()\n"
              << "  Scenario:    complex 1 has one molecule with 2 interfaces,\n"
              << "               complex 2 has one molecule with 1 interface.\n"
              << "  Criteria:    header count, header text, A/B labels, coords.\n";

    const std::string filename { "test_wxa_nominal.xyz" };

    // Molecule 0 -> complex 1 ('A'), molecule 1 -> complex 2 ('B').
    const Coord com0 { 1.5, -2.25, 3.0 };
    const Coord if0a { 2.5, -2.25, 3.0 };
    const Coord if0b { 1.5, -1.25, 3.0 };
    const Coord com1 { -10.0, 0.5, 7.75 };
    const Coord if1a { -9.0, 0.5, 7.75 };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wxa_make_molecule(com0, { if0a, if0b }, 2));
    moleculeList.push_back(wxa_make_molecule(com1, { if1a }, 1));

    Complex reactCom1 = wxa_make_complex({ 0 });
    Complex reactCom2 = wxa_make_complex({ 1 });

    std::cerr << "  Calling write_xyz_assoc(\"" << filename << "\", ...)\n";
    write_xyz_assoc(filename, reactCom1, reactCom2, moleculeList);

    const std::vector<std::string> lines = wxa_read_lines(filename);
    std::cerr << "  File contains " << lines.size() << " line(s)\n";

    // Expected: 2 members + 2 interfaces + 1 interface = 5 atoms.
    const int expectedAtoms = 5;
    const std::size_t expectedLines = 2u + static_cast<std::size_t>(expectedAtoms);

    EXPECT_FALSE(lines.empty()) << "write_xyz_assoc must create a readable file";
    EXPECT_EQ(lines.size(), expectedLines)
        << "file should hold 2 header lines plus one line per written atom";

    if (lines.size() >= 2u) {
        // First line: atom count.
        int reportedAtoms = -1;
        std::istringstream(lines[0]) >> reportedAtoms;
        std::cerr << "    header atom count = " << reportedAtoms
                  << " (expected " << expectedAtoms << ")\n";
        EXPECT_EQ(reportedAtoms, expectedAtoms)
            << "first line must be members(2) + interfaces(3) = 5";

        // Second line: fixed comment string.
        std::cerr << "    header comment    = \"" << lines[1] << "\"\n";
        EXPECT_EQ(lines[1], std::string("mol output final"))
            << "second line must be the literal comment 'mol output final'";
    }

    // Label bookkeeping: 3 'A' lines (com + 2 ifaces), 2 'B' lines (com + 1 iface).
    const int numA = wxa_count_label(lines, 'A');
    const int numB = wxa_count_label(lines, 'B');
    std::cerr << "    'A' lines = " << numA << " (expect 3), 'B' lines = " << numB
              << " (expect 2)\n";
    EXPECT_EQ(numA, 3) << "complex 1 should emit its COM plus its two interfaces";
    EXPECT_EQ(numB, 2) << "complex 2 should emit its COM plus its one interface";

    // Verify the actual coordinate values, in order.
    if (lines.size() == expectedLines) {
        wxa_expect_line_matches_coord(lines[2], com0, "A com");
        wxa_expect_line_matches_coord(lines[3], if0a, "A iface 0");
        wxa_expect_line_matches_coord(lines[4], if0b, "A iface 1");
        wxa_expect_line_matches_coord(lines[5], com1, "B com");
        wxa_expect_line_matches_coord(lines[6], if1a, "B iface 0");

        // Ordering guarantee: all 'A' lines must precede all 'B' lines.
        EXPECT_EQ(lines[2][0], 'A');
        EXPECT_EQ(lines[3][0], 'A');
        EXPECT_EQ(lines[4][0], 'A');
        EXPECT_EQ(lines[5][0], 'B');
        EXPECT_EQ(lines[6][0], 'B');
    }

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 2: Multiple molecules per complex.
//
// Pass criteria: the atom count and per-label line counts scale with the number
// of member molecules, and the members are written in memberList order.
// -----------------------------------------------------------------------------
void test_wxa_multiple_members_per_complex()
{
    std::cerr << "\n[TEST] test_wxa_multiple_members_per_complex\n"
              << "  Source file: src/io/write_xyz_assoc.cpp\n"
              << "  Function:    write_xyz_assoc()\n"
              << "  Scenario:    complex 1 owns molecules {0,1}, complex 2 owns {2}.\n"
              << "  Criteria:    atom count = 3 members + 4 interfaces = 7, and the\n"
              << "               member molecules appear in memberList order.\n";

    const std::string filename { "test_wxa_multi.xyz" };

    const Coord c0 { 0.0, 0.0, 0.0 };
    const Coord c1 { 1.0, 1.0, 1.0 };
    const Coord c2 { -5.0, -5.0, -5.0 };

    std::vector<Molecule> moleculeList;
    // Molecule 0: 1 interface, molecule 1: 2 interfaces, molecule 2: 1 interface.
    moleculeList.push_back(wxa_make_molecule(c0, { Coord { 0.5, 0.0, 0.0 } }, 1));
    moleculeList.push_back(
        wxa_make_molecule(c1, { Coord { 1.5, 1.0, 1.0 }, Coord { 1.0, 1.5, 1.0 } }, 2));
    moleculeList.push_back(wxa_make_molecule(c2, { Coord { -4.5, -5.0, -5.0 } }, 1));

    Complex reactCom1 = wxa_make_complex({ 0, 1 });
    Complex reactCom2 = wxa_make_complex({ 2 });

    std::cerr << "  Calling write_xyz_assoc(\"" << filename << "\", ...)\n";
    write_xyz_assoc(filename, reactCom1, reactCom2, moleculeList);

    const std::vector<std::string> lines = wxa_read_lines(filename);
    std::cerr << "  File contains " << lines.size() << " line(s)\n";

    // 3 member COMs + 4 interfaces = 7 atoms -> 9 lines total.
    EXPECT_EQ(lines.size(), 9u) << "2 header lines + 7 atom lines expected";

    if (!lines.empty()) {
        int reportedAtoms = -1;
        std::istringstream(lines[0]) >> reportedAtoms;
        std::cerr << "    header atom count = " << reportedAtoms << " (expect 7)\n";
        EXPECT_EQ(reportedAtoms, 7) << "3 COMs + 4 interfaces = 7 atoms";
    }

    const int numA = wxa_count_label(lines, 'A');
    const int numB = wxa_count_label(lines, 'B');
    std::cerr << "    'A' lines = " << numA << " (expect 5), 'B' lines = " << numB
              << " (expect 2)\n";
    EXPECT_EQ(numA, 5) << "complex 1: 2 COMs + 3 interfaces";
    EXPECT_EQ(numB, 2) << "complex 2: 1 COM + 1 interface";

    // Member ordering: molecule 0's COM must be written before molecule 1's COM.
    if (lines.size() == 9u) {
        wxa_expect_line_matches_coord(lines[2], c0, "A member 0 com");
        wxa_expect_line_matches_coord(lines[4], c1, "A member 1 com");
        wxa_expect_line_matches_coord(lines[7], c2, "B member 0 com");
    }

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 3: Both complexes empty.
//
// Pass criteria: the file is still created, the atom count is 0, and only the
// two header lines are present.
// -----------------------------------------------------------------------------
void test_wxa_empty_complexes()
{
    std::cerr << "\n[TEST] test_wxa_empty_complexes\n"
              << "  Source file: src/io/write_xyz_assoc.cpp\n"
              << "  Function:    write_xyz_assoc()\n"
              << "  Scenario:    both complexes have empty memberLists.\n"
              << "  Criteria:    file exists, atom count is 0, only headers written.\n";

    const std::string filename { "test_wxa_empty.xyz" };

    std::vector<Molecule> moleculeList; // no molecules at all
    Complex reactCom1 = wxa_make_complex({});
    Complex reactCom2 = wxa_make_complex({});

    std::cerr << "  Calling write_xyz_assoc(\"" << filename << "\", ...)\n";
    write_xyz_assoc(filename, reactCom1, reactCom2, moleculeList);

    const std::vector<std::string> lines = wxa_read_lines(filename);
    std::cerr << "  File contains " << lines.size() << " line(s)\n";

    EXPECT_EQ(lines.size(), 2u) << "an empty association should write only the headers";

    if (!lines.empty()) {
        int reportedAtoms = -1;
        std::istringstream(lines[0]) >> reportedAtoms;
        std::cerr << "    header atom count = " << reportedAtoms << " (expect 0)\n";
        EXPECT_EQ(reportedAtoms, 0) << "no members -> zero atoms";
    }
    if (lines.size() >= 2u) {
        EXPECT_EQ(lines[1], std::string("mol output final"))
            << "comment line must still be written";
    }

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 4: Documented quirk - the header count comes from interfaceList while the
//         emitted coordinate lines come from tmpICoords.
//
// Here molecule 0 declares 3 interfaces but only has 1 temporary interface
// coordinate.  The header therefore reports 1 COM + 3 interfaces = 4 atoms while
// only 2 coordinate lines are actually written.
//
// Pass criteria: header count == 4, written atom lines == 2 (i.e. the mismatch
// is reproduced, proving the count is driven by interfaceList).
// -----------------------------------------------------------------------------
void test_wxa_count_uses_interfaceList_not_tmpICoords()
{
    std::cerr << "\n[TEST] test_wxa_count_uses_interfaceList_not_tmpICoords\n"
              << "  Source file: src/io/write_xyz_assoc.cpp\n"
              << "  Function:    write_xyz_assoc()\n"
              << "  Scenario:    a molecule has 3 interfaceList entries but only 1\n"
              << "               tmpICoords entry (temp coords not fully populated).\n"
              << "  Criteria:    header count is driven by interfaceList (=4) while\n"
              << "               the number of coordinate lines follows tmpICoords (=2).\n";

    const std::string filename { "test_wxa_mismatch.xyz" };

    std::vector<Molecule> moleculeList;
    // 3 declared interfaces, but only one temporary interface coordinate.
    moleculeList.push_back(
        wxa_make_molecule(Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 0.0, 0.0 } }, 3));

    Complex reactCom1 = wxa_make_complex({ 0 });
    Complex reactCom2 = wxa_make_complex({}); // nothing on the 'B' side

    std::cerr << "  Calling write_xyz_assoc(\"" << filename << "\", ...)\n";
    write_xyz_assoc(filename, reactCom1, reactCom2, moleculeList);

    const std::vector<std::string> lines = wxa_read_lines(filename);
    std::cerr << "  File contains " << lines.size() << " line(s)\n";

    // Header count: 1 member + 3 interfaceList entries = 4.
    if (!lines.empty()) {
        int reportedAtoms = -1;
        std::istringstream(lines[0]) >> reportedAtoms;
        std::cerr << "    header atom count = " << reportedAtoms
                  << " (expect 4, from interfaceList)\n";
        EXPECT_EQ(reportedAtoms, 4)
            << "the count is computed from interfaceList.size(), not tmpICoords.size()";
    }

    // Written lines: 1 COM + 1 tmp interface coordinate = 2 atom lines.
    const int writtenAtomLines = static_cast<int>(lines.size()) - 2;
    std::cerr << "    written atom lines = " << writtenAtomLines
              << " (expect 2, from tmpICoords)\n";
    EXPECT_EQ(writtenAtomLines, 2)
        << "only tmpComCoord plus the single tmpICoords entry can be written";
    EXPECT_EQ(wxa_count_label(lines, 'A'), 2) << "both written lines belong to complex 1";
    EXPECT_EQ(wxa_count_label(lines, 'B'), 0) << "complex 2 is empty, so no 'B' lines";

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 5: Re-writing the same filename truncates the previous contents.
//
// Pass criteria: after a second, smaller dump the file must contain only the new
// data (std::ofstream opens in truncating mode by default).
// -----------------------------------------------------------------------------
void test_wxa_overwrites_existing_file()
{
    std::cerr << "\n[TEST] test_wxa_overwrites_existing_file\n"
              << "  Source file: src/io/write_xyz_assoc.cpp\n"
              << "  Function:    write_xyz_assoc()\n"
              << "  Scenario:    call the writer twice on the same filename, the\n"
              << "               second call dumping fewer atoms.\n"
              << "  Criteria:    file is truncated, not appended.\n";

    const std::string filename { "test_wxa_overwrite.xyz" };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(
        wxa_make_molecule(Coord { 0.0, 0.0, 0.0 },
            { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } }, 2));
    moleculeList.push_back(
        wxa_make_molecule(Coord { 5.0, 5.0, 5.0 }, { Coord { 6.0, 5.0, 5.0 } }, 1));

    // First dump: both complexes populated -> 2 COMs + 3 interfaces = 5 atoms.
    Complex big1 = wxa_make_complex({ 0 });
    Complex big2 = wxa_make_complex({ 1 });
    std::cerr << "  First call: dumping 5 atoms\n";
    write_xyz_assoc(filename, big1, big2, moleculeList);

    const std::size_t firstLineCount = wxa_read_lines(filename).size();
    std::cerr << "    line count after first call  = " << firstLineCount << "\n";
    EXPECT_EQ(firstLineCount, 7u) << "2 headers + 5 atoms expected on the first dump";

    // Second dump: only molecule 1 in complex 2 -> 1 COM + 1 interface = 2 atoms.
    Complex small1 = wxa_make_complex({});
    Complex small2 = wxa_make_complex({ 1 });
    std::cerr << "  Second call: dumping 2 atoms to the same filename\n";
    write_xyz_assoc(filename, small1, small2, moleculeList);

    const std::vector<std::string> lines = wxa_read_lines(filename);
    std::cerr << "    line count after second call = " << lines.size() << "\n";
    EXPECT_EQ(lines.size(), 4u)
        << "the file must be truncated: 2 headers + 2 atoms only";

    if (!lines.empty()) {
        int reportedAtoms = -1;
        std::istringstream(lines[0]) >> reportedAtoms;
        std::cerr << "    header atom count = " << reportedAtoms << " (expect 2)\n";
        EXPECT_EQ(reportedAtoms, 2) << "the new header must describe the new contents";
    }
    EXPECT_EQ(wxa_count_label(lines, 'A'), 0) << "complex 1 is empty on the second dump";
    EXPECT_EQ(wxa_count_label(lines, 'B'), 2) << "complex 2 contributes COM + 1 interface";

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 6: The routine reports the *temporary* association coordinates, not the
//         permanent ones.  A molecule whose comCoord differs strongly from its
//         tmpComCoord proves which member is used.
// -----------------------------------------------------------------------------
void test_wxa_uses_temporary_coordinates()
{
    std::cerr << "\n[TEST] test_wxa_uses_temporary_coordinates\n"
              << "  Source file: src/io/write_xyz_assoc.cpp\n"
              << "  Function:    write_xyz_assoc()\n"
              << "  Scenario:    comCoord and tmpComCoord differ; likewise the\n"
              << "               permanent interface coord differs from tmpICoords.\n"
              << "  Criteria:    the file contains the tmp* values.\n";

    const std::string filename { "test_wxa_tmpcrds.xyz" };

    const Coord permanentCom { 100.0, 100.0, 100.0 };
    const Coord tmpCom { -1.0, -2.0, -3.0 };
    const Coord tmpIface { -1.0, -2.0, -4.0 };

    Molecule mol = wxa_make_molecule(tmpCom, { tmpIface }, 1);
    // Deliberately set the permanent coordinates to something very different.
    mol.comCoord = permanentCom;
    mol.interfaceList[0].coord = permanentCom;

    std::vector<Molecule> moleculeList { mol };
    Complex reactCom1 = wxa_make_complex({ 0 });
    Complex reactCom2 = wxa_make_complex({});

    std::cerr << "  Calling write_xyz_assoc(\"" << filename << "\", ...)\n";
    write_xyz_assoc(filename, reactCom1, reactCom2, moleculeList);

    const std::vector<std::string> lines = wxa_read_lines(filename);
    std::cerr << "  File contains " << lines.size() << " line(s)\n";
    EXPECT_EQ(lines.size(), 4u) << "2 headers + COM + 1 interface expected";

    if (lines.size() == 4u) {
        wxa_expect_line_matches_coord(lines[2], tmpCom, "tmpComCoord");
        wxa_expect_line_matches_coord(lines[3], tmpIface, "tmpICoords[0]");

        // Guard against the permanent coordinates leaking into the output.
        EXPECT_EQ(lines[2].find("100"), std::string::npos)
            << "permanent comCoord (100,100,100) must not appear in the dump";
        EXPECT_EQ(lines[3].find("100"), std::string::npos)
            << "permanent interface coord must not appear in the dump";
    }

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: every named helper runs inside its own TEST so that a
// failure in one scenario does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(WriteXyzAssoc, NominalTwoComplexes) { test_wxa_nominal_two_complexes(); }
TEST(WriteXyzAssoc, MultipleMembersPerComplex) { test_wxa_multiple_members_per_complex(); }
TEST(WriteXyzAssoc, EmptyComplexes) { test_wxa_empty_complexes(); }
TEST(WriteXyzAssoc, CountUsesInterfaceList) { test_wxa_count_uses_interfaceList_not_tmpICoords(); }
TEST(WriteXyzAssoc, OverwritesExistingFile) { test_wxa_overwrites_existing_file(); }
TEST(WriteXyzAssoc, UsesTemporaryCoordinates) { test_wxa_uses_temporary_coordinates(); }