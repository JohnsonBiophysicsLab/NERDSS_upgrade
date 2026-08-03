/*! \file test_write_xyz.cpp
 *
 * ### Unit test for ../src/io/write_xyz.cpp
 *
 * Function under test:
 *
 *     void write_xyz(std::string filename, const Parameters& params,
 *                    const std::vector<Molecule>& moleculeList,
 *                    const std::vector<MolTemplate>& molTemplateList)
 *
 * Behaviour that is verified here (taken directly from the implementation):
 *
 *   1. Line 1 of the file is `params.numTotalUnits`.
 *   2. Line 2 of the file is the literal string "mol output final".
 *   3. For every non-implicit-lipid Molecule the molecule name is written
 *      followed by the center-of-mass coordinate, then one line per
 *      interface (same name, interface coordinate).
 *   4. Molecules whose MolTemplate has `isImplicitLipid == true` are skipped
 *      entirely (neither COM nor interfaces are written).
 *   5. If fewer coordinate lines were written than `params.numTotalUnits`,
 *      the remainder is padded with lines of the form
 *      `EMTY 150.000000 150.000000 150.000000`.
 *   6. If more lines than `numTotalUnits` were already written, no padding
 *      occurs (the while loop condition is immediately false).
 *
 * Every test writes to a temporary file, reads it back, and then deletes it.
 * Verbose progress information is printed to stderr so the reader can follow
 * exactly which behaviour is being probed and what the pass criteria are.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (file-local so they cannot collide with other test files).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal MolTemplate with a name and implicit-lipid flag.
 *
 * write_xyz() only ever reads `molName` and `isImplicitLipid` from a
 * MolTemplate, so nothing else needs to be populated.
 */
MolTemplate wxyz_make_template(const std::string& name, int typeIndex, bool isImplicitLipid = false)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.isImplicitLipid = isImplicitLipid;
    return temp;
}

/*! \brief Build a Molecule of a given type with a COM and a list of interfaces.
 *
 * \param[in] typeIndex index into molTemplateList (used for the printed name).
 * \param[in] com       center of mass coordinate written on the first line.
 * \param[in] ifaceCrds one coordinate per interface, each printed on its own line.
 */
Molecule wxyz_make_molecule(int typeIndex, const Coord& com, const std::vector<Coord>& ifaceCrds)
{
    Molecule mol;
    mol.molTypeIndex = typeIndex;
    mol.comCoord = com;
    mol.interfaceList.clear();
    for (const auto& crd : ifaceCrds) {
        Molecule::Iface iface;
        iface.coord = crd;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Read a whole text file into a vector of lines.
 *
 * Empty trailing lines produced by the final std::endl are not included.
 */
std::vector<std::string> wxyz_read_lines(const std::string& filename)
{
    std::vector<std::string> lines;
    std::ifstream in(filename);
    std::string line;
    while (std::getline(in, line)) {
        // Skip a completely empty trailing line, if any.
        if (line.empty() && in.eof())
            break;
        lines.push_back(line);
    }
    return lines;
}

/*! \brief Return the first whitespace-delimited token of a line.
 *
 * write_xyz() pads the molecule name with std::setw(4), so leading blanks are
 * expected and must be stripped before comparing names.
 */
std::string wxyz_first_token(const std::string& line)
{
    std::istringstream iss(line);
    std::string token;
    iss >> token;
    return token;
}

/*! \brief Extract every floating point number found in a string.
 *
 * The exact separator used by `operator<<(std::ostream&, const Coord&)` is an
 * implementation detail (space, tab, comma, ...), so instead of assuming a
 * layout we simply scan the line for numeric literals.  This keeps the test
 * robust while still verifying the numeric content.
 */
std::vector<double> wxyz_extract_doubles(const std::string& line)
{
    std::vector<double> vals;
    const char* p = line.c_str();
    while (*p != '\0') {
        const bool startsNumber = std::isdigit(static_cast<unsigned char>(*p))
            || ((*p == '-' || *p == '+' || *p == '.')
                && std::isdigit(static_cast<unsigned char>(*(p + 1))));
        if (startsNumber) {
            char* end = nullptr;
            const double v = std::strtod(p, &end);
            if (end == p) {
                ++p; // Defensive: could not convert, move on.
                continue;
            }
            vals.push_back(v);
            p = end;
        } else {
            ++p;
        }
    }
    return vals;
}

/*! \brief Count how many lines start with the "EMTY" padding token. */
std::size_t wxyz_count_empty_lines(const std::vector<std::string>& lines)
{
    std::size_t count { 0 };
    for (const auto& line : lines) {
        if (wxyz_first_token(line) == "EMTY")
            ++count;
    }
    return count;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: header + normal molecule/interface output.
// -----------------------------------------------------------------------------
void test_wxyz_basic_output()
{
    std::cerr << "\n[TEST] test_wxyz_basic_output\n"
              << "  Source file:   src/io/write_xyz.cpp\n"
              << "  Function:      write_xyz()\n"
              << "  Scenario:      two molecules (3 + 2 coordinate lines) and\n"
              << "                 numTotalUnits set exactly to 5.\n"
              << "  Pass criteria: line 1 == numTotalUnits, line 2 == the\n"
              << "                 'mol output final' banner, then one COM line\n"
              << "                 and one line per interface, in order, with\n"
              << "                 the correct molecule names and coordinates.\n";

    // Two molecule types: "AAA" (2 interfaces) and "BBB" (1 interface).
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wxyz_make_template("AAA", 0));
    molTemplateList.push_back(wxyz_make_template("BBB", 1));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wxyz_make_molecule(0, Coord { 1.0, 2.0, 3.0 },
        { Coord { 1.5, 2.5, 3.5 }, Coord { -1.5, -2.5, -3.5 } }));
    moleculeList.push_back(wxyz_make_molecule(1, Coord { 10.0, 20.0, 30.0 },
        { Coord { 11.0, 21.0, 31.0 } }));

    // 1 COM + 2 ifaces + 1 COM + 1 iface == 5 coordinate lines expected.
    Parameters params;
    params.numTotalUnits = 5;

    const std::string filename { "test_wxyz_basic_output.xyz" };
    std::cerr << "  Writing file '" << filename << "'...\n";
    write_xyz(filename, params, moleculeList, molTemplateList);

    const std::vector<std::string> lines = wxyz_read_lines(filename);
    std::cerr << "  Read back " << lines.size() << " lines.\n";

    // Two header lines + five coordinate lines.
    EXPECT_EQ(lines.size(), static_cast<std::size_t>(7))
        << "Expected 2 header lines plus 5 coordinate lines";

    if (lines.size() >= 2) {
        // Header line 1: the total number of units.
        EXPECT_EQ(wxyz_first_token(lines[0]), std::string("5"))
            << "First line must be params.numTotalUnits";
        // Header line 2: fixed banner text.
        EXPECT_EQ(lines[1], std::string("mol output final"))
            << "Second line must be the literal banner 'mol output final'";
    }

    // Expected sequence of names and coordinates.
    const std::vector<std::string> expectNames { "AAA", "AAA", "AAA", "BBB", "BBB" };
    const std::vector<Coord> expectCrds {
        Coord { 1.0, 2.0, 3.0 },
        Coord { 1.5, 2.5, 3.5 },
        Coord { -1.5, -2.5, -3.5 },
        Coord { 10.0, 20.0, 30.0 },
        Coord { 11.0, 21.0, 31.0 },
    };

    for (std::size_t i = 0; i < expectNames.size(); ++i) {
        const std::size_t lineIdx = i + 2; // skip the two header lines
        if (lineIdx >= lines.size()) {
            ADD_FAILURE() << "Missing coordinate line " << i;
            continue;
        }
        std::cerr << "    line " << lineIdx << ": \"" << lines[lineIdx] << "\"\n";

        // Name check (leading setw(4) padding is stripped by first_token).
        EXPECT_EQ(wxyz_first_token(lines[lineIdx]), expectNames[i])
            << "Molecule name mismatch on coordinate line " << i;

        // Coordinate check: exactly three numbers with the expected values.
        const std::vector<double> vals = wxyz_extract_doubles(lines[lineIdx]);
        EXPECT_EQ(vals.size(), static_cast<std::size_t>(3))
            << "Each coordinate line should contain exactly 3 numbers";
        if (vals.size() == 3) {
            EXPECT_NEAR(vals[0], expectCrds[i].x, 1e-4) << "x mismatch on line " << i;
            EXPECT_NEAR(vals[1], expectCrds[i].y, 1e-4) << "y mismatch on line " << i;
            EXPECT_NEAR(vals[2], expectCrds[i].z, 1e-4) << "z mismatch on line " << i;
        }
    }

    // No padding should have been needed.
    EXPECT_EQ(wxyz_count_empty_lines(lines), static_cast<std::size_t>(0))
        << "No EMTY padding expected when numWritten == numTotalUnits";

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 2: molecules belonging to an implicit-lipid template are skipped.
// -----------------------------------------------------------------------------
void test_wxyz_implicit_lipid_skipped()
{
    std::cerr << "\n[TEST] test_wxyz_implicit_lipid_skipped\n"
              << "  Source file:   src/io/write_xyz.cpp\n"
              << "  Function:      write_xyz()\n"
              << "  Scenario:      molTemplateList[1].isImplicitLipid == true and\n"
              << "                 one molecule of that type is in the list.\n"
              << "  Pass criteria: the implicit lipid name never appears in the\n"
              << "                 file and only the explicit molecule's lines\n"
              << "                 (COM + 1 interface) are written.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wxyz_make_template("PRO", 0, /*isImplicitLipid=*/false));
    molTemplateList.push_back(wxyz_make_template("ILIP", 1, /*isImplicitLipid=*/true));

    std::vector<Molecule> moleculeList;
    // Explicit molecule: contributes 2 lines (COM + 1 interface).
    moleculeList.push_back(wxyz_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 0.0, 0.0 } }));
    // Implicit lipid: should contribute nothing at all.
    moleculeList.push_back(wxyz_make_molecule(1, Coord { 99.0, 99.0, 99.0 }, { Coord { 98.0, 98.0, 98.0 } }));

    Parameters params;
    params.numTotalUnits = 2; // exactly the two explicit lines -> no padding

    const std::string filename { "test_wxyz_implicit_lipid.xyz" };
    std::cerr << "  Writing file '" << filename << "'...\n";
    write_xyz(filename, params, moleculeList, molTemplateList);

    const std::vector<std::string> lines = wxyz_read_lines(filename);
    std::cerr << "  Read back " << lines.size() << " lines.\n";

    // 2 headers + 2 coordinate lines from the explicit molecule only.
    EXPECT_EQ(lines.size(), static_cast<std::size_t>(4))
        << "Implicit lipid molecule must not produce any output lines";

    // Verify the implicit lipid name is nowhere in the file.
    bool foundImplicitName { false };
    for (const auto& line : lines) {
        if (line.find("ILIP") != std::string::npos)
            foundImplicitName = true;
    }
    EXPECT_FALSE(foundImplicitName) << "Implicit lipid name 'ILIP' should never be written";

    // The two written lines should belong to the explicit molecule.
    if (lines.size() >= 4) {
        EXPECT_EQ(wxyz_first_token(lines[2]), std::string("PRO")) << "First body line should be the explicit molecule COM";
        EXPECT_EQ(wxyz_first_token(lines[3]), std::string("PRO")) << "Second body line should be the explicit interface";
    }

    // And no padding, since numWritten == numTotalUnits == 2.
    EXPECT_EQ(wxyz_count_empty_lines(lines), static_cast<std::size_t>(0))
        << "No EMTY padding expected here";

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 3: padding with EMTY lines when fewer units were written than declared.
// -----------------------------------------------------------------------------
void test_wxyz_pads_with_empty_lines()
{
    std::cerr << "\n[TEST] test_wxyz_pads_with_empty_lines\n"
              << "  Source file:   src/io/write_xyz.cpp\n"
              << "  Function:      write_xyz()\n"
              << "  Scenario:      only 2 coordinate lines are produced but\n"
              << "                 numTotalUnits == 6.\n"
              << "  Pass criteria: 4 trailing 'EMTY' lines with the coordinate\n"
              << "                 (150, 150, 150) are appended.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wxyz_make_template("MOL", 0));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wxyz_make_molecule(0, Coord { 5.0, 5.0, 5.0 }, { Coord { 6.0, 5.0, 5.0 } }));

    Parameters params;
    params.numTotalUnits = 6; // 2 real lines -> 4 EMTY lines expected

    const std::string filename { "test_wxyz_padding.xyz" };
    std::cerr << "  Writing file '" << filename << "'...\n";
    write_xyz(filename, params, moleculeList, molTemplateList);

    const std::vector<std::string> lines = wxyz_read_lines(filename);
    std::cerr << "  Read back " << lines.size() << " lines.\n";

    // 2 headers + 2 real lines + 4 padded lines == 8.
    EXPECT_EQ(lines.size(), static_cast<std::size_t>(8))
        << "File should contain 2 headers, 2 molecule lines and 4 padded lines";

    // Exactly four padding lines.
    EXPECT_EQ(wxyz_count_empty_lines(lines), static_cast<std::size_t>(4))
        << "Exactly numTotalUnits - numWritten == 4 EMTY lines expected";

    // Each padding line must carry the sentinel coordinate 150,150,150.
    for (std::size_t i = 4; i < lines.size(); ++i) {
        std::cerr << "    padding line " << i << ": \"" << lines[i] << "\"\n";
        EXPECT_EQ(wxyz_first_token(lines[i]), std::string("EMTY"))
            << "Padding lines must be labelled EMTY";
        const std::vector<double> vals = wxyz_extract_doubles(lines[i]);
        EXPECT_EQ(vals.size(), static_cast<std::size_t>(3))
            << "Padding line should contain 3 coordinate values";
        if (vals.size() == 3) {
            EXPECT_NEAR(vals[0], 150.0, 1e-4) << "Padding x should be 150";
            EXPECT_NEAR(vals[1], 150.0, 1e-4) << "Padding y should be 150";
            EXPECT_NEAR(vals[2], 150.0, 1e-4) << "Padding z should be 150";
        }
    }

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 4: no padding when more lines were written than numTotalUnits declares.
// -----------------------------------------------------------------------------
void test_wxyz_no_padding_when_over_count()
{
    std::cerr << "\n[TEST] test_wxyz_no_padding_when_over_count\n"
              << "  Source file:   src/io/write_xyz.cpp\n"
              << "  Function:      write_xyz()\n"
              << "  Scenario:      3 coordinate lines are produced while\n"
              << "                 numTotalUnits == 1 (deliberately too small).\n"
              << "  Pass criteria: the padding loop never runs, so no 'EMTY'\n"
              << "                 line is present and all molecule data is kept.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wxyz_make_template("BIG", 0));

    std::vector<Molecule> moleculeList;
    // 1 COM + 2 interfaces == 3 coordinate lines.
    moleculeList.push_back(wxyz_make_molecule(0, Coord { 0.0, 0.0, 0.0 },
        { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } }));

    Parameters params;
    params.numTotalUnits = 1;

    const std::string filename { "test_wxyz_over_count.xyz" };
    std::cerr << "  Writing file '" << filename << "'...\n";
    write_xyz(filename, params, moleculeList, molTemplateList);

    const std::vector<std::string> lines = wxyz_read_lines(filename);
    std::cerr << "  Read back " << lines.size() << " lines.\n";

    // 2 headers + 3 coordinate lines, no padding.
    EXPECT_EQ(lines.size(), static_cast<std::size_t>(5))
        << "All molecule lines must be written even if numTotalUnits is too small";
    EXPECT_EQ(wxyz_count_empty_lines(lines), static_cast<std::size_t>(0))
        << "Padding loop must not execute when numWritten >= numTotalUnits";

    // The header still reports the (incorrect) user-supplied count verbatim.
    if (!lines.empty()) {
        EXPECT_EQ(wxyz_first_token(lines[0]), std::string("1"))
            << "Header must echo params.numTotalUnits unmodified";
    }

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 5: empty molecule list -> file consists of headers plus pure padding.
// -----------------------------------------------------------------------------
void test_wxyz_empty_molecule_list()
{
    std::cerr << "\n[TEST] test_wxyz_empty_molecule_list\n"
              << "  Source file:   src/io/write_xyz.cpp\n"
              << "  Function:      write_xyz()\n"
              << "  Scenario:      moleculeList is empty and numTotalUnits == 3.\n"
              << "  Pass criteria: the file still opens successfully and holds\n"
              << "                 2 header lines plus 3 'EMTY' padding lines.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wxyz_make_template("NONE", 0));

    std::vector<Molecule> moleculeList; // deliberately empty

    Parameters params;
    params.numTotalUnits = 3;

    const std::string filename { "test_wxyz_empty_list.xyz" };
    std::cerr << "  Writing file '" << filename << "'...\n";
    write_xyz(filename, params, moleculeList, molTemplateList);

    // Verify the file was actually created (the function opens it in ctor).
    std::ifstream check(filename);
    EXPECT_TRUE(check.good()) << "write_xyz should create the output file";
    check.close();

    const std::vector<std::string> lines = wxyz_read_lines(filename);
    std::cerr << "  Read back " << lines.size() << " lines.\n";

    EXPECT_EQ(lines.size(), static_cast<std::size_t>(5))
        << "Expected 2 header lines and 3 padded lines";
    EXPECT_EQ(wxyz_count_empty_lines(lines), static_cast<std::size_t>(3))
        << "All three coordinate slots should be filled with EMTY placeholders";

    if (lines.size() >= 2) {
        EXPECT_EQ(wxyz_first_token(lines[0]), std::string("3"))
            << "Header must report numTotalUnits == 3";
        EXPECT_EQ(lines[1], std::string("mol output final"))
            << "Banner line must be present even with no molecules";
    }

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// Test 6: a molecule with no interfaces contributes only its COM line, and
//         calling write_xyz twice on the same file truncates (does not append).
// -----------------------------------------------------------------------------
void test_wxyz_no_interfaces_and_overwrite()
{
    std::cerr << "\n[TEST] test_wxyz_no_interfaces_and_overwrite\n"
              << "  Source file:   src/io/write_xyz.cpp\n"
              << "  Function:      write_xyz()\n"
              << "  Scenario:      (a) a molecule with an empty interfaceList,\n"
              << "                 (b) the same filename written twice.\n"
              << "  Pass criteria: (a) exactly one coordinate line per such\n"
              << "                 molecule, (b) the second write replaces the\n"
              << "                 first (std::ofstream truncates by default).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wxyz_make_template("PNT", 0));

    // A single "point" molecule: COM only, no interfaces.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wxyz_make_molecule(0, Coord { -7.25, 8.5, 0.0 }, {}));

    Parameters params;
    params.numTotalUnits = 1;

    const std::string filename { "test_wxyz_overwrite.xyz" };

    // ---- First write: two molecules, so the file is comparatively long.
    std::vector<Molecule> longerList = moleculeList;
    longerList.push_back(wxyz_make_molecule(0, Coord { 1.0, 1.0, 1.0 },
        { Coord { 2.0, 1.0, 1.0 }, Coord { 1.0, 2.0, 1.0 } }));
    Parameters longerParams;
    longerParams.numTotalUnits = 4;

    std::cerr << "  First write (4 coordinate lines)...\n";
    write_xyz(filename, longerParams, longerList, molTemplateList);
    const std::vector<std::string> firstLines = wxyz_read_lines(filename);
    std::cerr << "  First write produced " << firstLines.size() << " lines.\n";
    EXPECT_EQ(firstLines.size(), static_cast<std::size_t>(6))
        << "First write should contain 2 headers + 4 coordinate lines";

    // ---- Second write: only the point molecule, file must be truncated.
    std::cerr << "  Second write (1 coordinate line) to the same filename...\n";
    write_xyz(filename, params, moleculeList, molTemplateList);
    const std::vector<std::string> secondLines = wxyz_read_lines(filename);
    std::cerr << "  Second write produced " << secondLines.size() << " lines.\n";

    // 2 headers + 1 COM line (no interfaces, no padding needed).
    EXPECT_EQ(secondLines.size(), static_cast<std::size_t>(3))
        << "Second write must truncate the file, not append to it";
    EXPECT_EQ(wxyz_count_empty_lines(secondLines), static_cast<std::size_t>(0))
        << "No padding expected: 1 line written for numTotalUnits == 1";

    if (secondLines.size() >= 3) {
        EXPECT_EQ(wxyz_first_token(secondLines[2]), std::string("PNT"))
            << "The single body line should be the point molecule's COM";
        const std::vector<double> vals = wxyz_extract_doubles(secondLines[2]);
        EXPECT_EQ(vals.size(), static_cast<std::size_t>(3))
            << "COM line should carry three coordinate values";
        if (vals.size() == 3) {
            EXPECT_NEAR(vals[0], -7.25, 1e-4) << "COM x should be -7.25";
            EXPECT_NEAR(vals[1], 8.5, 1e-4) << "COM y should be 8.5";
            EXPECT_NEAR(vals[2], 0.0, 1e-4) << "COM z should be 0.0";
        }
    }

    std::remove(filename.c_str());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each helper is run inside its own TEST so that a failure
// in one scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(WriteXyzTest, BasicOutput) { test_wxyz_basic_output(); }
TEST(WriteXyzTest, ImplicitLipidSkipped) { test_wxyz_implicit_lipid_skipped(); }
TEST(WriteXyzTest, PadsWithEmptyLines) { test_wxyz_pads_with_empty_lines(); }
TEST(WriteXyzTest, NoPaddingWhenOverCount) { test_wxyz_no_padding_when_over_count(); }
TEST(WriteXyzTest, EmptyMoleculeList) { test_wxyz_empty_molecule_list(); }
TEST(WriteXyzTest, NoInterfacesAndOverwrite) { test_wxyz_no_interfaces_and_overwrite(); }