/*! \file test_write_psf.cpp
 *
 * ### Unit test for ../src/io/write_psf.cpp
 *
 * The function under test is:
 *
 *     void write_psf(const Parameters& params,
 *                    const std::vector<Molecule>& moleculeList,
 *                    const std::vector<MolTemplate>& molTemplateList)
 *
 * `write_psf()` writes a CHARMM-style PSF topology file describing every
 * molecule (center of mass + interfaces) in the system to the fixed path
 * "DATA/system.psf".  It has no return value, so the only way to test it is to
 * call it with a hand-built system and then parse the file it produced.
 *
 * The tests below check, for a set of small synthetic systems:
 *   - the fixed PSF header block (title, remark lines, static counters),
 *   - the "!NATOM" count and the number of atom records actually written
 *     (1 COM record + 1 record per interface, per molecule),
 *   - that molecules whose MolTemplate is an implicit lipid are skipped,
 *   - the "EMY" (empty) padding records emitted up to params.numTotalUnits,
 *   - the per-interface segment/colour naming rule
 *     (interfaces 0-2 -> "N", 3-5 -> "LPA", 6+ -> "O"; the COM -> "O"),
 *   - the "!NBOND" count and the number of bond indices written.
 *
 * All assertions are non-fatal (EXPECT_*) so every test runs to completion.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <sys/stat.h>

#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// Path constants matching the hard-coded output location inside write_psf().
const char* kWpsfDir = "DATA";
const char* kWpsfFile = "DATA/system.psf";

/*! \brief Make sure the DATA/ directory exists so the ofstream can open.
 *
 * write_psf() opens "DATA/system.psf" without creating the directory, so the
 * test harness has to do it.  Failure (e.g. because it already exists) is
 * intentionally ignored.
 */
void wpsf_ensure_data_dir()
{
    mkdir(kWpsfDir, 0755);
}

/*! \brief Delete any previously written PSF so each test starts clean. */
void wpsf_remove_output()
{
    std::remove(kWpsfFile);
}

/*! \brief Read the whole PSF file back as a vector of lines.
 *
 * Returns an empty vector if the file could not be opened.
 */
std::vector<std::string> wpsf_read_lines()
{
    std::vector<std::string> lines;
    std::ifstream in(kWpsfFile);
    if (!in.is_open())
        return lines;
    std::string line;
    while (std::getline(in, line))
        lines.push_back(line);
    return lines;
}

/*! \brief Split a line on whitespace into tokens (PSF records are padded). */
std::vector<std::string> wpsf_tokens(const std::string& line)
{
    std::vector<std::string> toks;
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok)
        toks.push_back(tok);
    return toks;
}

/*! \brief Find the index of the first line containing \p needle, or -1. */
int wpsf_find_line(const std::vector<std::string>& lines, const std::string& needle)
{
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].find(needle) != std::string::npos)
            return static_cast<int>(i);
    }
    return -1;
}

/*! \brief Parse the leading integer of a line such as "      6 !NATOM".
 *
 * Returns -1 if the line has no leading integer.
 */
int wpsf_leading_int(const std::string& line)
{
    std::istringstream iss(line);
    int value = -1;
    if (iss >> value)
        return value;
    return -1;
}

/*! \brief Build a MolTemplate with \p numIfaces named interfaces.
 *
 * \param typeIndex  value stored in molTypeIndex (also the molTemplateList slot)
 * \param name       molecule name; write_psf() prints only name.substr(0,3)
 * \param numIfaces  number of interfaces on the template
 * \param bonds      optional COM/interface bond list (drives the !NBOND count)
 * \param isImplicit marks the template as an implicit lipid (skipped by writer)
 */
MolTemplate wpsf_make_template(int typeIndex, const std::string& name, int numIfaces,
    const std::vector<std::array<int, 2>>& bonds, bool isImplicit)
{
    MolTemplate temp;
    temp.molTypeIndex = typeIndex;
    temp.molName = name;
    temp.isImplicitLipid = isImplicit;
    temp.bondList = bonds;

    // Interface names are written with substr(0,3), so keep them 3 chars long
    // and unique: "i00", "i01", ...
    temp.interfaceList.clear();
    for (int i = 0; i < numIfaces; ++i) {
        Interface iface;
        std::ostringstream oss;
        oss << 'i' << (i / 10) << (i % 10);
        iface.name = oss.str();
        iface.index = i;
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Build a Molecule matching a template's interface count.
 *
 * write_psf() iterates over mol.interfaceList but indexes
 * molTemplateList[...].interfaceList, so the two sizes must agree.
 */
Molecule wpsf_make_molecule(int index, int typeIndex, int numIfaces, double mass)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = typeIndex;
    mol.mass = mass;
    mol.myComIndex = index;
    mol.comCoord = Coord(0.0, 0.0, 0.0);

    mol.interfaceList.clear();
    for (int i = 0; i < numIfaces; ++i) {
        Molecule::Iface iface;
        iface.relIndex = i;
        iface.molTypeIndex = typeIndex;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: The fixed header block, including the two static counters that
//         write_psf() reports in its REMARKS lines.
// -----------------------------------------------------------------------------
void test_wpsf_header_block()
{
    std::cerr << "\n[TEST] test_wpsf_header_block\n"
              << "  Source file:   src/io/write_psf.cpp\n"
              << "  Function:      write_psf()\n"
              << "  Scenario:      one 2-interface molecule; inspect the PSF header.\n"
              << "  Pass criteria: line 0 == \"PSF CMAP CHEQ\", a \"2 !NTITLE\" line,\n"
              << "                 and REMARKS lines echoing Molecule::numberOfMolecules\n"
              << "                 and Complex::numberOfComplexes.\n";

    wpsf_ensure_data_dir();
    wpsf_remove_output();

    // Static counters are echoed verbatim into the REMARKS lines.
    Molecule::numberOfMolecules = 7;
    Complex::numberOfComplexes = 3;

    std::vector<MolTemplate> molTemplateList { wpsf_make_template(0, "PROT", 2, {}, false) };
    std::vector<Molecule> moleculeList { wpsf_make_molecule(0, 0, 2, 10.0) };

    Parameters params;
    params.numTotalUnits = 3; // 1 COM + 2 interfaces -> no EMY padding

    std::cerr << "  Calling write_psf()...\n";
    write_psf(params, moleculeList, molTemplateList);

    std::vector<std::string> lines = wpsf_read_lines();
    std::cerr << "  Read back " << lines.size() << " lines from " << kWpsfFile << '\n';

    // The file must exist and have at least the 8-line header block.
    ASSERT_GE(lines.size(), 8u) << "PSF header requires at least 8 lines; was DATA/ writable?";

    // Line 0: magic header string.
    EXPECT_EQ(lines[0], "PSF CMAP CHEQ") << "First line must be the PSF magic string";

    // Line 1: blank (the header is written as "PSF CMAP CHEQ\n\n").
    EXPECT_TRUE(lines[1].empty()) << "Second line should be blank";

    // Line 2: right-padded "2" followed by " !NTITLE".
    EXPECT_NE(lines[2].find("!NTITLE"), std::string::npos) << "Third line should declare !NTITLE";
    EXPECT_EQ(wpsf_leading_int(lines[2]), 2) << "NTITLE count is hard-coded to 2";

    // Lines 3-5: the three REMARKS lines.
    EXPECT_EQ(lines[3], "REMARKS PSF for entire system") << "First REMARKS line is a fixed string";
    EXPECT_NE(lines[4].find("total molecules: 7"), std::string::npos)
        << "REMARKS should report Molecule::numberOfMolecules (7)";
    EXPECT_NE(lines[5].find("total complexes: 3"), std::string::npos)
        << "REMARKS should report Complex::numberOfComplexes (3)";

    // Line 6 blank, line 7 the !NATOM declaration.
    EXPECT_TRUE(lines[6].empty()) << "A blank line precedes the !NATOM record";
    EXPECT_NE(lines[7].find("!NATOM"), std::string::npos) << "Line 7 should declare !NATOM";
}

// -----------------------------------------------------------------------------
// Test 2: !NATOM count and the atom records themselves.
// -----------------------------------------------------------------------------
void test_wpsf_atom_count_and_records()
{
    std::cerr << "\n[TEST] test_wpsf_atom_count_and_records\n"
              << "  Source file:   src/io/write_psf.cpp\n"
              << "  Function:      write_psf()\n"
              << "  Scenario:      two molecules with 2 and 3 interfaces.\n"
              << "  Pass criteria: !NATOM == (2+1)+(3+1) == 7, exactly 7 atom records\n"
              << "                 are written, they are numbered 1..7, and each\n"
              << "                 carries the molecule name (first 3 chars) and mass.\n";

    wpsf_ensure_data_dir();
    wpsf_remove_output();

    Molecule::numberOfMolecules = 2;
    Complex::numberOfComplexes = 2;

    std::vector<MolTemplate> molTemplateList {
        wpsf_make_template(0, "ALPHA", 2, {}, false), // -> "ALP"
        wpsf_make_template(1, "BETA", 3, {}, false) //  -> "BET"
    };
    std::vector<Molecule> moleculeList {
        wpsf_make_molecule(0, 0, 2, 11.0),
        wpsf_make_molecule(1, 1, 3, 22.0)
    };

    const int expectedNAtom = (2 + 1) + (3 + 1); // 7

    Parameters params;
    params.numTotalUnits = static_cast<unsigned>(expectedNAtom); // no EMY padding

    std::cerr << "  Calling write_psf()...\n";
    write_psf(params, moleculeList, molTemplateList);

    std::vector<std::string> lines = wpsf_read_lines();
    ASSERT_FALSE(lines.empty()) << "No PSF output was produced";

    const int natomLine = wpsf_find_line(lines, "!NATOM");
    ASSERT_GE(natomLine, 0) << "Could not locate the !NATOM record";
    std::cerr << "  !NATOM line reads: \"" << lines[natomLine] << "\"\n";
    EXPECT_EQ(wpsf_leading_int(lines[natomLine]), expectedNAtom)
        << "!NATOM must equal sum over molecules of (interfaces + 1)";

    // Atom records run from the line after !NATOM up to the blank lines that
    // precede the !NBOND record.
    const int nbondLine = wpsf_find_line(lines, "!NBOND");
    ASSERT_GT(nbondLine, natomLine) << "Could not locate the !NBOND record after !NATOM";

    std::vector<std::string> atomLines;
    for (int i = natomLine + 1; i < nbondLine; ++i) {
        if (!wpsf_tokens(lines[i]).empty())
            atomLines.push_back(lines[i]);
    }
    std::cerr << "  Counted " << atomLines.size() << " non-blank atom records\n";
    EXPECT_EQ(static_cast<int>(atomLines.size()), expectedNAtom)
        << "Number of atom records must match the !NATOM count";

    // Every record must be sequentially numbered starting at 1 and hold 9 fields.
    for (size_t i = 0; i < atomLines.size(); ++i) {
        std::vector<std::string> toks = wpsf_tokens(atomLines[i]);
        EXPECT_EQ(toks.size(), 9u) << "Atom record " << (i + 1) << " should have 9 fields";
        if (toks.size() == 9u) {
            EXPECT_EQ(wpsf_leading_int(atomLines[i]), static_cast<int>(i + 1))
                << "Atom records must be numbered consecutively from 1";
        }
    }

    // Record 1 is the COM of the first molecule ("ALP" / "COM"), records 2-3 are
    // its interfaces; record 4 is the COM of the second molecule ("BET").
    if (atomLines.size() == static_cast<size_t>(expectedNAtom)) {
        std::vector<std::string> first = wpsf_tokens(atomLines[0]);
        EXPECT_EQ(first[1], "ALP") << "Molecule name is truncated to 3 characters";
        EXPECT_EQ(first[3], "COM") << "First record of a molecule is its center of mass";
        EXPECT_EQ(first[7], "11") << "COM record should carry the molecule mass (11)";

        std::vector<std::string> fourth = wpsf_tokens(atomLines[3]);
        EXPECT_EQ(fourth[1], "BET") << "Fourth record belongs to the second molecule";
        EXPECT_EQ(fourth[3], "COM") << "Fourth record is the second molecule's COM";
        EXPECT_EQ(fourth[7], "22") << "Second molecule's mass (22) should be written";
    }
}

// -----------------------------------------------------------------------------
// Test 3: Implicit-lipid molecules must be completely skipped.
// -----------------------------------------------------------------------------
void test_wpsf_implicit_lipid_skipped()
{
    std::cerr << "\n[TEST] test_wpsf_implicit_lipid_skipped\n"
              << "  Source file:   src/io/write_psf.cpp\n"
              << "  Function:      write_psf()\n"
              << "  Scenario:      one explicit molecule (2 ifaces) plus one implicit\n"
              << "                 lipid molecule (4 ifaces).\n"
              << "  Pass criteria: !NATOM == 3 (only the explicit molecule counts) and\n"
              << "                 the implicit lipid's name never appears in the file.\n";

    wpsf_ensure_data_dir();
    wpsf_remove_output();

    Molecule::numberOfMolecules = 2;
    Complex::numberOfComplexes = 2;

    std::vector<MolTemplate> molTemplateList {
        wpsf_make_template(0, "REAL", 2, {}, false), // explicit -> "REA"
        wpsf_make_template(1, "LIPID", 4, {}, true) // implicit -> "LIP"
    };
    std::vector<Molecule> moleculeList {
        wpsf_make_molecule(0, 0, 2, 5.0),
        wpsf_make_molecule(1, 1, 4, 6.0)
    };

    Parameters params;
    params.numTotalUnits = 3; // only the explicit molecule's 3 records

    std::cerr << "  Calling write_psf()...\n";
    write_psf(params, moleculeList, molTemplateList);

    std::vector<std::string> lines = wpsf_read_lines();
    ASSERT_FALSE(lines.empty()) << "No PSF output was produced";

    const int natomLine = wpsf_find_line(lines, "!NATOM");
    ASSERT_GE(natomLine, 0) << "Could not locate the !NATOM record";
    EXPECT_EQ(wpsf_leading_int(lines[natomLine]), 3)
        << "Implicit lipids must not contribute to the !NATOM count";

    // The truncated implicit-lipid name must not be present anywhere.
    EXPECT_EQ(wpsf_find_line(lines, "LIP"), -1)
        << "No record should be written for an implicit lipid molecule";
    EXPECT_GE(wpsf_find_line(lines, "REA"), 0)
        << "The explicit molecule should still be written";
}

// -----------------------------------------------------------------------------
// Test 4: "EMY" padding records are emitted up to params.numTotalUnits.
// -----------------------------------------------------------------------------
void test_wpsf_empty_padding_records()
{
    std::cerr << "\n[TEST] test_wpsf_empty_padding_records\n"
              << "  Source file:   src/io/write_psf.cpp\n"
              << "  Function:      write_psf()\n"
              << "  Scenario:      one molecule producing 3 records but\n"
              << "                 params.numTotalUnits == 6.\n"
              << "  Pass criteria: exactly 3 additional \"EMY\" placeholder records are\n"
              << "                 written, numbered 4, 5 and 6.\n";

    wpsf_ensure_data_dir();
    wpsf_remove_output();

    Molecule::numberOfMolecules = 1;
    Complex::numberOfComplexes = 1;

    std::vector<MolTemplate> molTemplateList { wpsf_make_template(0, "PROT", 2, {}, false) };
    std::vector<Molecule> moleculeList { wpsf_make_molecule(0, 0, 2, 4.0) };

    Parameters params;
    params.numTotalUnits = 6; // 3 real records + 3 padding records

    std::cerr << "  Calling write_psf()...\n";
    write_psf(params, moleculeList, molTemplateList);

    std::vector<std::string> lines = wpsf_read_lines();
    ASSERT_FALSE(lines.empty()) << "No PSF output was produced";

    // Collect all EMY placeholder records.
    std::vector<std::string> emyLines;
    for (const std::string& line : lines) {
        std::vector<std::string> toks = wpsf_tokens(line);
        if (toks.size() > 1 && toks[1] == "EMY")
            emyLines.push_back(line);
    }
    std::cerr << "  Found " << emyLines.size() << " EMY placeholder records\n";
    EXPECT_EQ(emyLines.size(), 3u)
        << "numTotalUnits(6) - writtenRecords(3) == 3 EMY records expected";

    // They continue the running atom index (4, 5, 6).
    for (size_t i = 0; i < emyLines.size(); ++i) {
        EXPECT_EQ(wpsf_leading_int(emyLines[i]), static_cast<int>(4 + i))
            << "EMY records must continue the atom numbering";
    }

    // Sanity: !NATOM still only counts the real molecule records.
    const int natomLine = wpsf_find_line(lines, "!NATOM");
    ASSERT_GE(natomLine, 0) << "Could not locate the !NATOM record";
    EXPECT_EQ(wpsf_leading_int(lines[natomLine]), 3)
        << "!NATOM counts real molecule records only, not EMY padding";
}

// -----------------------------------------------------------------------------
// Test 5: Interface segment/"colour" naming rule.
// -----------------------------------------------------------------------------
void test_wpsf_interface_colour_rule()
{
    std::cerr << "\n[TEST] test_wpsf_interface_colour_rule\n"
              << "  Source file:   src/io/write_psf.cpp\n"
              << "  Function:      write_psf()\n"
              << "  Scenario:      a single molecule with 8 interfaces.\n"
              << "  Pass criteria: COM record uses \"O\"; interfaces 0-2 use \"N\",\n"
              << "                 interfaces 3-5 use \"LPA\" and interfaces 6+ use \"O\";\n"
              << "                 interface names are the template names.\n";

    wpsf_ensure_data_dir();
    wpsf_remove_output();

    Molecule::numberOfMolecules = 1;
    Complex::numberOfComplexes = 1;

    const int numIfaces = 8;
    std::vector<MolTemplate> molTemplateList {
        wpsf_make_template(0, "MULTI", numIfaces, {}, false)
    };
    std::vector<Molecule> moleculeList { wpsf_make_molecule(0, 0, numIfaces, 9.0) };

    Parameters params;
    params.numTotalUnits = static_cast<unsigned>(numIfaces + 1); // no EMY records

    std::cerr << "  Calling write_psf()...\n";
    write_psf(params, moleculeList, molTemplateList);

    std::vector<std::string> lines = wpsf_read_lines();
    ASSERT_FALSE(lines.empty()) << "No PSF output was produced";

    const int natomLine = wpsf_find_line(lines, "!NATOM");
    const int nbondLine = wpsf_find_line(lines, "!NBOND");
    ASSERT_GE(natomLine, 0) << "Could not locate the !NATOM record";
    ASSERT_GT(nbondLine, natomLine) << "Could not locate the !NBOND record";

    std::vector<std::vector<std::string>> records;
    for (int i = natomLine + 1; i < nbondLine; ++i) {
        std::vector<std::string> toks = wpsf_tokens(lines[i]);
        if (!toks.empty())
            records.push_back(toks);
    }

    ASSERT_EQ(records.size(), static_cast<size_t>(numIfaces + 1))
        << "Expected 1 COM record plus " << numIfaces << " interface records";

    // Record 0 is the center of mass: name "COM", colour field "O".
    EXPECT_EQ(records[0][3], "COM") << "First record must be the COM";
    EXPECT_EQ(records[0][4], "O") << "COM record uses the \"O\" segment label";

    // Records 1..8 are the interfaces; verify the 3-way colour rule.
    for (int i = 0; i < numIfaces; ++i) {
        const std::vector<std::string>& rec = records[i + 1];
        ASSERT_EQ(rec.size(), 9u) << "Interface record " << i << " should have 9 fields";

        // Expected interface name is the template name ("i00", "i01", ...).
        std::ostringstream expectedName;
        expectedName << 'i' << (i / 10) << (i % 10);
        EXPECT_EQ(rec[3], expectedName.str())
            << "Interface record " << i << " should print the template interface name";

        const std::string expectedColour = (i < 3) ? "N" : ((i < 6) ? "LPA" : "O");
        EXPECT_EQ(rec[4], expectedColour)
            << "Interface " << i << " expected segment label \"" << expectedColour << '"';
        std::cerr << "    iface " << i << " -> name=" << rec[3] << " label=" << rec[4] << '\n';
    }
}

// -----------------------------------------------------------------------------
// Test 6: !NBOND count and number of bond indices written.
// -----------------------------------------------------------------------------
void test_wpsf_bond_section()
{
    std::cerr << "\n[TEST] test_wpsf_bond_section\n"
              << "  Source file:   src/io/write_psf.cpp\n"
              << "  Function:      write_psf()\n"
              << "  Scenario:      two molecules of a template holding 2 bonds each.\n"
              << "  Pass criteria: !NBOND == 4 and exactly 8 integers (2 per bond)\n"
              << "                 follow the !NBOND record.\n";

    wpsf_ensure_data_dir();
    wpsf_remove_output();

    Molecule::numberOfMolecules = 2;
    Complex::numberOfComplexes = 1;

    // COM(0) -> iface1(1) and COM(0) -> iface2(2).
    std::vector<std::array<int, 2>> bonds { { { 0, 1 } }, { { 0, 2 } } };
    std::vector<MolTemplate> molTemplateList { wpsf_make_template(0, "BOND", 2, bonds, false) };
    std::vector<Molecule> moleculeList {
        wpsf_make_molecule(0, 0, 2, 1.0),
        wpsf_make_molecule(1, 0, 2, 1.0)
    };

    Parameters params;
    params.numTotalUnits = 6; // 2 molecules x (2 ifaces + COM)

    std::cerr << "  Calling write_psf()...\n";
    write_psf(params, moleculeList, molTemplateList);

    std::vector<std::string> lines = wpsf_read_lines();
    ASSERT_FALSE(lines.empty()) << "No PSF output was produced";

    const int nbondLine = wpsf_find_line(lines, "!NBOND");
    ASSERT_GE(nbondLine, 0) << "Could not locate the !NBOND record";
    std::cerr << "  !NBOND line reads: \"" << lines[nbondLine] << "\"\n";
    EXPECT_EQ(wpsf_leading_int(lines[nbondLine]), 4)
        << "!NBOND must be the summed bondList size over all written molecules";

    // Everything after the !NBOND line is a whitespace-separated list of atom
    // indices, two per bond.
    int indexCount = 0;
    for (size_t i = nbondLine + 1; i < lines.size(); ++i) {
        indexCount += static_cast<int>(wpsf_tokens(lines[i]).size());
    }
    std::cerr << "  Counted " << indexCount << " bond indices after the !NBOND record\n";
    EXPECT_EQ(indexCount, 8) << "4 bonds x 2 atom indices == 8 integers expected";
}

// -----------------------------------------------------------------------------
// Test 7: A system without any bond definitions produces an empty bond section.
// -----------------------------------------------------------------------------
void test_wpsf_no_bonds_section()
{
    std::cerr << "\n[TEST] test_wpsf_no_bonds_section\n"
              << "  Source file:   src/io/write_psf.cpp\n"
              << "  Function:      write_psf()\n"
              << "  Scenario:      one molecule whose template declares no bonds.\n"
              << "  Pass criteria: !NBOND == 0 and no bond indices follow it.\n";

    wpsf_ensure_data_dir();
    wpsf_remove_output();

    Molecule::numberOfMolecules = 1;
    Complex::numberOfComplexes = 1;

    std::vector<MolTemplate> molTemplateList { wpsf_make_template(0, "NOBOND", 2, {}, false) };
    std::vector<Molecule> moleculeList { wpsf_make_molecule(0, 0, 2, 3.0) };

    Parameters params;
    params.numTotalUnits = 3;

    std::cerr << "  Calling write_psf()...\n";
    write_psf(params, moleculeList, molTemplateList);

    std::vector<std::string> lines = wpsf_read_lines();
    ASSERT_FALSE(lines.empty()) << "No PSF output was produced";

    const int nbondLine = wpsf_find_line(lines, "!NBOND");
    ASSERT_GE(nbondLine, 0) << "Could not locate the !NBOND record";
    EXPECT_EQ(wpsf_leading_int(lines[nbondLine]), 0)
        << "A template with an empty bondList must yield !NBOND == 0";

    // No numeric bond entries may follow.
    int indexCount = 0;
    for (size_t i = nbondLine + 1; i < lines.size(); ++i)
        indexCount += static_cast<int>(wpsf_tokens(lines[i]).size());
    std::cerr << "  Counted " << indexCount << " tokens after the !NBOND record\n";
    EXPECT_EQ(indexCount, 0) << "No bond indices should be written when there are no bonds";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named helper runs inside its own TEST so the
// framework reports individual results while still executing all of them.
// -----------------------------------------------------------------------------
TEST(WritePsf, HeaderBlock) { test_wpsf_header_block(); }
TEST(WritePsf, AtomCountAndRecords) { test_wpsf_atom_count_and_records(); }
TEST(WritePsf, ImplicitLipidSkipped) { test_wpsf_implicit_lipid_skipped(); }
TEST(WritePsf, EmptyPaddingRecords) { test_wpsf_empty_padding_records(); }
TEST(WritePsf, InterfaceColourRule) { test_wpsf_interface_colour_rule(); }
TEST(WritePsf, BondSection) { test_wpsf_bond_section(); }
TEST(WritePsf, NoBondsSection) { test_wpsf_no_bonds_section(); }