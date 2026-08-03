/*! \file test_write_traj.cpp
 *
 * ### Unit test for ../src/io/write_traj.cpp
 *
 * This test exercises the single function defined in that file:
 *
 *     void write_traj(long long int iter, std::ofstream& trajFile,
 *                     const Parameters& params,
 *                     const std::vector<Molecule>& moleculeList,
 *                     const std::vector<MolTemplate>& molTemplateList,
 *                     const Membrane& membraneObject)
 *
 * write_traj() appends one "frame" to an XYZ-style trajectory file:
 *
 *     <params.numTotalUnits>
 *     iteration: <iter>
 *     NAME  <center of mass coordinate>
 *     NAME  <interface 0 coordinate>
 *     NAME  <interface 1 coordinate>
 *     ...            (repeated for every non-empty, non-implicit molecule)
 *
 * Where NAME is the molecule template name truncated to 4 characters and
 * right-aligned in a field of width 4 (std::setw(4)).
 *
 * Because the function only produces text, every test writes to a temporary
 * file and then reads the file back to verify:
 *   - the two header lines,
 *   - the number and order of the coordinate lines,
 *   - that empty / implicit-lipid molecules are skipped,
 *   - that the molecule name is truncated/padded to exactly 4 characters,
 *   - that the coordinate text matches what Coord's operator<< produces
 *     when the stream is in std::fixed mode,
 *   - that data is flushed to disk before the stream is closed,
 *   - that repeated calls append additional frames.
 *
 * Verbose console output describes each step so a reader can follow along.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Small local helpers (prefixed wtraj_ so they cannot collide with other tests)
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a MolTemplate that only carries the fields write_traj() uses.
 *
 * write_traj() only reads MolTemplate::molName, so nothing else is needed.
 *
 * \param[in] name      Name of the molecule type (may be longer than 4 chars).
 * \param[in] typeIndex Index this template will occupy in molTemplateList.
 */
MolTemplate wtraj_make_template(const std::string& name, int typeIndex)
{
    MolTemplate molTemp;
    molTemp.molName = name;
    molTemp.molTypeIndex = typeIndex;
    return molTemp;
}

/*! \brief Build a Molecule with a COM coordinate and a set of interfaces.
 *
 * \param[in] typeIndex   Index into molTemplateList (selects the printed name).
 * \param[in] com         Center of mass coordinate.
 * \param[in] ifaceCoords Coordinates of each interface to attach.
 */
Molecule wtraj_make_molecule(int typeIndex, const Coord& com,
    const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.molTypeIndex = typeIndex;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isImplicitLipid = false;

    mol.interfaceList.clear();
    for (const auto& crd : ifaceCoords) {
        Molecule::Iface iface;
        iface.coord = crd;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Read every line of a text file into a vector of strings. */
std::vector<std::string> wtraj_read_lines(const std::string& fileName)
{
    std::vector<std::string> lines;
    std::ifstream inFile(fileName);
    std::string oneLine;
    while (std::getline(inFile, oneLine))
        lines.push_back(oneLine);
    return lines;
}

/*! \brief Reproduce exactly how write_traj() renders a Coord.
 *
 * write_traj() puts the stream into std::fixed and then uses Coord's
 * operator<<, so we do the same here to build the expected substring without
 * hard-coding any particular formatting of Coord.
 */
std::string wtraj_expected_coord_text(const Coord& crd)
{
    std::ostringstream oss;
    oss << std::fixed << crd;
    return oss.str();
}

/*! \brief Render a molecule-name field the way write_traj() does.
 *
 * The name is first truncated to 4 characters (substr(0,4)) and then written
 * with std::setw(4), i.e. right-aligned/padded to width 4.
 */
std::string wtraj_expected_name_field(const std::string& molName)
{
    std::ostringstream oss;
    oss << std::setw(4) << molName.substr(0, 4);
    return oss.str();
}

/*! \brief A tiny Membrane; write_traj() does not currently use it, but the
 *         signature requires one. */
Membrane wtraj_make_membrane()
{
    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    return membraneObject;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Header lines + total number of written lines for a simple system.
// -----------------------------------------------------------------------------
void test_wtraj_header_and_line_count()
{
    std::cerr << "\n[TEST] test_wtraj_header_and_line_count\n"
              << "  Source file:   src/io/write_traj.cpp\n"
              << "  Function:      write_traj()\n"
              << "  Scenario:      2 molecules, 2 interfaces each, iteration 7.\n"
              << "  Pass criteria: line 1 == params.numTotalUnits, line 2 ==\n"
              << "                 \"iteration: 7\", and 6 coordinate lines follow\n"
              << "                 (1 COM + 2 interfaces per molecule).\n";

    const std::string fileName { "test_wtraj_header_and_line_count.xyz" };

    // --- Set up the fake system -------------------------------------------
    Parameters params;
    params.numTotalUnits = 6; // 2 molecules * (1 COM + 2 ifaces)

    std::vector<MolTemplate> molTemplateList { wtraj_make_template("AAAA", 0) };

    std::vector<Molecule> moleculeList {
        wtraj_make_molecule(0, Coord { 0.0, 0.0, 0.0 },
            { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } }),
        wtraj_make_molecule(0, Coord { 5.0, 5.0, 5.0 },
            { Coord { 6.0, 5.0, 5.0 }, Coord { 5.0, 6.0, 5.0 } })
    };

    Membrane membraneObject { wtraj_make_membrane() };

    // --- Call the function under test -------------------------------------
    std::cerr << "  Writing one frame to \"" << fileName << "\"...\n";
    {
        std::ofstream trajFile(fileName);
        ASSERT_TRUE(trajFile.is_open()) << "Could not open temporary trajectory file";
        write_traj(7, trajFile, params, moleculeList, molTemplateList, membraneObject);
    }

    // --- Verify -----------------------------------------------------------
    std::vector<std::string> lines = wtraj_read_lines(fileName);
    std::cerr << "  File contains " << lines.size() << " lines.\n";

    // Expect 2 header lines + 6 coordinate lines.
    EXPECT_EQ(lines.size(), static_cast<size_t>(8))
        << "Expected 2 header lines plus 6 coordinate lines";

    if (lines.size() >= 2) {
        std::cerr << "  Line 1: \"" << lines[0] << "\"\n";
        std::cerr << "  Line 2: \"" << lines[1] << "\"\n";
        // First line is the total number of units.
        EXPECT_EQ(lines[0], std::string("6"))
            << "First line must be params.numTotalUnits";
        // Second line is the iteration marker.
        EXPECT_EQ(lines[1], std::string("iteration: 7"))
            << "Second line must read \"iteration: <iter>\"";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 2: Empty and implicit-lipid molecules are skipped entirely.
// -----------------------------------------------------------------------------
void test_wtraj_skips_empty_and_implicit_molecules()
{
    std::cerr << "\n[TEST] test_wtraj_skips_empty_and_implicit_molecules\n"
              << "  Source file:   src/io/write_traj.cpp\n"
              << "  Function:      write_traj()\n"
              << "  Scenario:      3 molecules: one normal, one isEmpty, one\n"
              << "                 isImplicitLipid.\n"
              << "  Pass criteria: only the normal molecule contributes lines\n"
              << "                 (1 COM + 1 interface = 2 coordinate lines).\n";

    const std::string fileName { "test_wtraj_skips_empty_and_implicit.xyz" };

    Parameters params;
    params.numTotalUnits = 2;

    std::vector<MolTemplate> molTemplateList {
        wtraj_make_template("GOOD", 0),
        wtraj_make_template("SKIP", 1)
    };

    // Molecule 0: a normal molecule that must be written.
    Molecule good = wtraj_make_molecule(0, Coord { 1.0, 2.0, 3.0 },
        { Coord { 1.5, 2.0, 3.0 } });

    // Molecule 1: flagged as destroyed/empty -> must be skipped.
    Molecule emptyMol = wtraj_make_molecule(1, Coord { 10.0, 10.0, 10.0 },
        { Coord { 11.0, 10.0, 10.0 } });
    emptyMol.isEmpty = true;

    // Molecule 2: an implicit lipid -> must be skipped.
    Molecule implicitMol = wtraj_make_molecule(1, Coord { -10.0, -10.0, -10.0 },
        { Coord { -11.0, -10.0, -10.0 } });
    implicitMol.isImplicitLipid = true;

    std::vector<Molecule> moleculeList { good, emptyMol, implicitMol };
    Membrane membraneObject { wtraj_make_membrane() };

    std::cerr << "  Writing frame with 1 valid + 2 skipped molecules...\n";
    {
        std::ofstream trajFile(fileName);
        ASSERT_TRUE(trajFile.is_open()) << "Could not open temporary trajectory file";
        write_traj(0, trajFile, params, moleculeList, molTemplateList, membraneObject);
    }

    std::vector<std::string> lines = wtraj_read_lines(fileName);
    std::cerr << "  File contains " << lines.size() << " lines (expect 4).\n";

    // 2 header lines + 2 coordinate lines from the single valid molecule.
    EXPECT_EQ(lines.size(), static_cast<size_t>(4))
        << "Empty and implicit-lipid molecules must not be written";

    // No line may mention the skipped template name.
    for (size_t i = 0; i < lines.size(); ++i) {
        EXPECT_EQ(lines[i].find("SKIP"), std::string::npos)
            << "Line " << i << " (\"" << lines[i]
            << "\") must not reference a skipped molecule";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 3: Molecule names are truncated to 4 characters and padded to width 4.
// -----------------------------------------------------------------------------
void test_wtraj_name_truncation_and_padding()
{
    std::cerr << "\n[TEST] test_wtraj_name_truncation_and_padding\n"
              << "  Source file:   src/io/write_traj.cpp\n"
              << "  Function:      write_traj()\n"
              << "  Scenario:      one template with a long name (\"LONGNAME\")\n"
              << "                 and one with a short name (\"AB\").\n"
              << "  Pass criteria: printed field is substr(0,4) of the name,\n"
              << "                 right-aligned in a width-4 column, followed\n"
              << "                 by a single space.\n";

    const std::string fileName { "test_wtraj_name_truncation.xyz" };

    Parameters params;
    params.numTotalUnits = 2;

    std::vector<MolTemplate> molTemplateList {
        wtraj_make_template("LONGNAME", 0),
        wtraj_make_template("AB", 1)
    };

    // Each molecule has only a COM (no interfaces) to keep the file tiny.
    std::vector<Molecule> moleculeList {
        wtraj_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, {}),
        wtraj_make_molecule(1, Coord { 1.0, 1.0, 1.0 }, {})
    };

    Membrane membraneObject { wtraj_make_membrane() };

    std::cerr << "  Writing frame with long and short molecule names...\n";
    {
        std::ofstream trajFile(fileName);
        ASSERT_TRUE(trajFile.is_open()) << "Could not open temporary trajectory file";
        write_traj(1, trajFile, params, moleculeList, molTemplateList, membraneObject);
    }

    std::vector<std::string> lines = wtraj_read_lines(fileName);
    ASSERT_GE(lines.size(), static_cast<size_t>(4))
        << "Expected 2 header lines plus 2 coordinate lines";

    const std::string expectedLong = wtraj_expected_name_field("LONGNAME"); // "LONG"
    const std::string expectedShort = wtraj_expected_name_field("AB"); // "  AB"

    std::cerr << "  Coordinate line 1: \"" << lines[2] << "\"\n"
              << "    expected name field: \"" << expectedLong << "\"\n";
    std::cerr << "  Coordinate line 2: \"" << lines[3] << "\"\n"
              << "    expected name field: \"" << expectedShort << "\"\n";

    // The first 4 characters of each coordinate line are the name field.
    EXPECT_EQ(lines[2].substr(0, 4), expectedLong)
        << "\"LONGNAME\" should be truncated to its first 4 characters";
    EXPECT_EQ(lines[3].substr(0, 4), expectedShort)
        << "\"AB\" should be right-aligned in a width-4 field";

    // A single space separates the name field from the coordinate text.
    EXPECT_EQ(lines[2][4], ' ') << "A space must follow the name field";
    EXPECT_EQ(lines[3][4], ' ') << "A space must follow the name field";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 4: COM line comes first, then interfaces in order, with correct text.
// -----------------------------------------------------------------------------
void test_wtraj_coordinate_order_and_values()
{
    std::cerr << "\n[TEST] test_wtraj_coordinate_order_and_values\n"
              << "  Source file:   src/io/write_traj.cpp\n"
              << "  Function:      write_traj()\n"
              << "  Scenario:      one molecule with 3 interfaces at distinct\n"
              << "                 coordinates.\n"
              << "  Pass criteria: the COM line is written first, then each\n"
              << "                 interface in list order, and each line ends\n"
              << "                 with the std::fixed rendering of its Coord.\n";

    const std::string fileName { "test_wtraj_coordinate_order.xyz" };

    Parameters params;
    params.numTotalUnits = 4;

    std::vector<MolTemplate> molTemplateList { wtraj_make_template("PRO1", 0) };

    const Coord comCrd { 1.25, -2.5, 3.75 };
    const Coord ifaceA { 10.5, 0.0, 0.0 };
    const Coord ifaceB { 0.0, 20.25, 0.0 };
    const Coord ifaceC { 0.0, 0.0, -30.125 };

    std::vector<Molecule> moleculeList {
        wtraj_make_molecule(0, comCrd, { ifaceA, ifaceB, ifaceC })
    };

    Membrane membraneObject { wtraj_make_membrane() };

    std::cerr << "  Writing frame with one molecule and three interfaces...\n";
    {
        std::ofstream trajFile(fileName);
        ASSERT_TRUE(trajFile.is_open()) << "Could not open temporary trajectory file";
        write_traj(42, trajFile, params, moleculeList, molTemplateList, membraneObject);
    }

    std::vector<std::string> lines = wtraj_read_lines(fileName);
    ASSERT_GE(lines.size(), static_cast<size_t>(6))
        << "Expected 2 header lines plus 4 coordinate lines";

    // Build the expected coordinate text for each line, in write order.
    const std::vector<Coord> expectedOrder { comCrd, ifaceA, ifaceB, ifaceC };

    for (size_t i = 0; i < expectedOrder.size(); ++i) {
        const std::string expectedText = wtraj_expected_coord_text(expectedOrder[i]);
        const std::string& actualLine = lines[i + 2]; // skip the 2 header lines
        std::cerr << "  Coordinate line " << i << ": \"" << actualLine << "\"\n"
                  << "    expected to contain: \"" << expectedText << "\"\n";

        // Each line must carry the correct molecule name...
        EXPECT_EQ(actualLine.substr(0, 4), wtraj_expected_name_field("PRO1"))
            << "Line " << i << " should be labelled with the template name";
        // ...and the correct coordinate text, in the correct order.
        EXPECT_NE(actualLine.find(expectedText), std::string::npos)
            << "Line " << i << " should contain the coordinate rendering \""
            << expectedText << '"';
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 5: Data is flushed to disk before the stream is closed.
// -----------------------------------------------------------------------------
void test_wtraj_output_is_flushed()
{
    std::cerr << "\n[TEST] test_wtraj_output_is_flushed\n"
              << "  Source file:   src/io/write_traj.cpp\n"
              << "  Function:      write_traj()\n"
              << "  Scenario:      read the file back while the ofstream is\n"
              << "                 still open.\n"
              << "  Pass criteria: the frame is already visible on disk because\n"
              << "                 write_traj() calls std::endl / std::flush.\n";

    const std::string fileName { "test_wtraj_flush.xyz" };

    Parameters params;
    params.numTotalUnits = 2;

    std::vector<MolTemplate> molTemplateList { wtraj_make_template("FLSH", 0) };
    std::vector<Molecule> moleculeList {
        wtraj_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 1.0, 1.0 } })
    };
    Membrane membraneObject { wtraj_make_membrane() };

    // Deliberately keep the stream open while reading the file back.
    std::ofstream trajFile(fileName);
    ASSERT_TRUE(trajFile.is_open()) << "Could not open temporary trajectory file";
    std::cerr << "  Writing frame (stream deliberately left open)...\n";
    write_traj(3, trajFile, params, moleculeList, molTemplateList, membraneObject);

    std::vector<std::string> lines = wtraj_read_lines(fileName);
    std::cerr << "  Lines visible on disk before close: " << lines.size() << '\n';

    // 2 header lines + 2 coordinate lines should already be present.
    EXPECT_EQ(lines.size(), static_cast<size_t>(4))
        << "write_traj() must flush its output as it goes";
    if (!lines.empty())
        EXPECT_EQ(lines[0], std::string("2")) << "Header should already be flushed";

    trajFile.close();
    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 6: Repeated calls append additional frames rather than overwriting.
// -----------------------------------------------------------------------------
void test_wtraj_appends_multiple_frames()
{
    std::cerr << "\n[TEST] test_wtraj_appends_multiple_frames\n"
              << "  Source file:   src/io/write_traj.cpp\n"
              << "  Function:      write_traj()\n"
              << "  Scenario:      call write_traj() twice on the same stream\n"
              << "                 with iterations 0 and 100.\n"
              << "  Pass criteria: the file holds two complete frames and both\n"
              << "                 iteration markers appear in order.\n";

    const std::string fileName { "test_wtraj_multi_frame.xyz" };

    Parameters params;
    params.numTotalUnits = 2;

    std::vector<MolTemplate> molTemplateList { wtraj_make_template("MULT", 0) };
    std::vector<Molecule> moleculeList {
        wtraj_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, { Coord { 2.0, 0.0, 0.0 } })
    };
    Membrane membraneObject { wtraj_make_membrane() };

    std::cerr << "  Writing frame for iteration 0 then iteration 100...\n";
    {
        std::ofstream trajFile(fileName);
        ASSERT_TRUE(trajFile.is_open()) << "Could not open temporary trajectory file";
        write_traj(0, trajFile, params, moleculeList, molTemplateList, membraneObject);
        write_traj(100, trajFile, params, moleculeList, molTemplateList, membraneObject);
    }

    std::vector<std::string> lines = wtraj_read_lines(fileName);
    std::cerr << "  File contains " << lines.size() << " lines (expect 8).\n";

    // Each frame is 2 header lines + 2 coordinate lines => 8 total.
    EXPECT_EQ(lines.size(), static_cast<size_t>(8))
        << "Two frames should have been appended";

    if (lines.size() >= 6) {
        std::cerr << "  Frame 1 iteration line: \"" << lines[1] << "\"\n";
        std::cerr << "  Frame 2 iteration line: \"" << lines[5] << "\"\n";
        EXPECT_EQ(lines[1], std::string("iteration: 0"))
            << "First frame must report iteration 0";
        EXPECT_EQ(lines[5], std::string("iteration: 100"))
            << "Second frame must report iteration 100";
        // Both frames repeat the unit-count header.
        EXPECT_EQ(lines[4], std::string("2"))
            << "Second frame must start with the unit count again";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 7: An all-empty system still produces the two header lines.
// -----------------------------------------------------------------------------
void test_wtraj_no_writable_molecules()
{
    std::cerr << "\n[TEST] test_wtraj_no_writable_molecules\n"
              << "  Source file:   src/io/write_traj.cpp\n"
              << "  Function:      write_traj()\n"
              << "  Scenario:      moleculeList is empty.\n"
              << "  Pass criteria: only the two header lines are written; no\n"
              << "                 coordinate lines and no crash.\n";

    const std::string fileName { "test_wtraj_no_molecules.xyz" };

    Parameters params;
    params.numTotalUnits = 0;

    std::vector<MolTemplate> molTemplateList { wtraj_make_template("NONE", 0) };
    std::vector<Molecule> moleculeList {}; // no molecules at all
    Membrane membraneObject { wtraj_make_membrane() };

    std::cerr << "  Writing frame for an empty molecule list...\n";
    {
        std::ofstream trajFile(fileName);
        ASSERT_TRUE(trajFile.is_open()) << "Could not open temporary trajectory file";
        write_traj(11, trajFile, params, moleculeList, molTemplateList, membraneObject);
    }

    std::vector<std::string> lines = wtraj_read_lines(fileName);
    std::cerr << "  File contains " << lines.size() << " lines (expect 2).\n";

    EXPECT_EQ(lines.size(), static_cast<size_t>(2))
        << "Only the two header lines should be written";
    if (lines.size() >= 2) {
        EXPECT_EQ(lines[0], std::string("0")) << "Unit count header should be 0";
        EXPECT_EQ(lines[1], std::string("iteration: 11"))
            << "Iteration header should still be written";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each helper is run inside its own TEST so that a failure
// in one scenario does not stop the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(WriteTraj, HeaderAndLineCount) { test_wtraj_header_and_line_count(); }
TEST(WriteTraj, SkipsEmptyAndImplicitMolecules) { test_wtraj_skips_empty_and_implicit_molecules(); }
TEST(WriteTraj, NameTruncationAndPadding) { test_wtraj_name_truncation_and_padding(); }
TEST(WriteTraj, CoordinateOrderAndValues) { test_wtraj_coordinate_order_and_values(); }
TEST(WriteTraj, OutputIsFlushed) { test_wtraj_output_is_flushed(); }
TEST(WriteTraj, AppendsMultipleFrames) { test_wtraj_appends_multiple_frames(); }
TEST(WriteTraj, NoWritableMolecules) { test_wtraj_no_writable_molecules(); }