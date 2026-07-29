/*! \file test_write_crds.cpp
 *
 * ### Unit test for ../src/io/write_crds.cpp
 *
 * The file under test provides the overloaded error-dump routine:
 *
 *     void write_crds(const std::vector<Complex>& Complexlist,
 *                     const std::vector<Molecule>& bases)
 *
 * Behaviour of the function (from the source):
 *   - It loops over every Complex in `Complexlist` using the *positional*
 *     index `i` (NOT Complex::index).
 *   - For every molecule index `memMol` inside `Complexlist[i].memberList`
 *     it opens the file  "out/c<i>_p<memMol>_error.dat".
 *   - It writes one header line containing
 *         bases[memMol].molTypeIndex  <space>  bases[memMol].myComIndex
 *   - It then delegates to Molecule::write_crd_file(out) which appends the
 *     coordinate block for that molecule.
 *
 * Therefore the tests below verify:
 *   1. The naming scheme of the produced files (positional complex index and
 *      the molecule index taken from memberList).
 *   2. The contents of the header line (molTypeIndex then myComIndex).
 *   3. That the coordinate block is actually appended (file has > 1 line).
 *   4. Degenerate inputs (empty complex list / empty memberList) create no
 *      files and do not crash.
 *
 * NOTE: the function writes into a relative directory called "out"; the test
 * creates that directory first, otherwise the std::ofstream would silently
 * fail and nothing would be written.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (all prefixed with wc_ so they cannot collide with other tests)
// -----------------------------------------------------------------------------
namespace {

/*! \brief Make sure the relative "out" directory exists.
 *
 * write_crds() opens files under "out/"; if that directory is missing the
 * ofstream construction fails silently and no output is produced. Creating it
 * here makes the test deterministic regardless of where the binary is run.
 */
void wc_make_out_dir()
{
    // 0755 permissions; ignore the error if it already exists.
    mkdir("out", 0755);
}

/*! \brief Build the exact filename write_crds() is expected to produce.
 *
 * \param[in] complexPos positional index of the Complex within Complexlist
 * \param[in] molIndex   molecule index as stored in Complex::memberList
 */
std::string wc_expected_name(int complexPos, int molIndex)
{
    return "out/c" + std::to_string(complexPos) + "_p" + std::to_string(molIndex)
        + "_error.dat";
}

/*! \brief Returns true if the given path can be opened for reading. */
bool wc_file_exists(const std::string& path)
{
    std::ifstream in(path);
    return in.good();
}

/*! \brief Reads the whole file into a string (empty string if unreadable). */
std::string wc_read_file(const std::string& path)
{
    std::ifstream in(path);
    if (!in.good())
        return std::string {};
    std::ostringstream oss;
    oss << in.rdbuf();
    return oss.str();
}

/*! \brief Counts the number of newline-terminated lines in a file. */
int wc_count_lines(const std::string& path)
{
    std::ifstream in(path);
    if (!in.good())
        return 0;
    int count { 0 };
    std::string line;
    while (std::getline(in, line))
        ++count;
    return count;
}

/*! \brief Deletes a file if it exists (best effort, ignores failures). */
void wc_remove_file(const std::string& path) { std::remove(path.c_str()); }

/*! \brief Construct a minimal, self-consistent Molecule for the dump.
 *
 * The molecule gets a center of mass plus a single interface so that
 * Molecule::write_crd_file() has something to print.
 *
 * \param[in] molTypeIndex value written as the first field of the header line
 * \param[in] myComIndex   value written as the second field of the header line
 * \param[in] com          center-of-mass coordinate
 */
Molecule wc_make_molecule(int molTypeIndex, int myComIndex, const Coord& com)
{
    Molecule mol;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = myComIndex;
    mol.comCoord = com;
    mol.mass = 1.0;
    mol.isEmpty = false;

    // One interface located slightly off the COM.
    Molecule::Iface iface;
    iface.coord = Coord { com.x + 1.0, com.y, com.z };
    iface.index = 0;
    iface.relIndex = 0;
    iface.molTypeIndex = molTypeIndex;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Construct a Complex owning the supplied molecule indices.
 *
 * \param[in] index       the Complex::index field (deliberately allowed to
 *                        differ from the positional index to prove the
 *                        function uses the loop counter for the file name)
 * \param[in] memberList  molecule indices belonging to this complex
 */
Complex wc_make_complex(int index, const std::vector<int>& memberList)
{
    Complex com;
    com.index = index;
    com.memberList = memberList;
    com.comCoord = Coord { 0.0, 0.0, 0.0 };
    com.radius = 1.0;
    com.mass = 1.0;
    com.isEmpty = false;
    return com;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a single complex with a single molecule.
//         Verifies the file name, the header line values, and that the
//         coordinate block from write_crd_file() was appended.
// -----------------------------------------------------------------------------
void test_wc_single_complex_single_molecule()
{
    std::cerr << "\n[TEST] test_wc_single_complex_single_molecule\n"
              << "  Source file:   src/io/write_crds.cpp\n"
              << "  Function:      write_crds(const std::vector<Complex>&,\n"
              << "                            const std::vector<Molecule>&)\n"
              << "  Scenario:      one Complex owning one Molecule.\n"
              << "  Pass criteria: file out/c0_p0_error.dat exists, its first\n"
              << "                 line is \"<molTypeIndex> <myComIndex>\", and\n"
              << "                 additional coordinate lines follow.\n";

    wc_make_out_dir();

    // Molecule 0 belongs to complex 0 and is of molecule type 3.
    std::vector<Molecule> bases;
    bases.push_back(wc_make_molecule(/*molTypeIndex=*/3, /*myComIndex=*/0,
        Coord { 1.5, -2.5, 3.5 }));

    // A single complex whose memberList holds molecule index 0.
    std::vector<Complex> complexList;
    complexList.push_back(wc_make_complex(/*index=*/0, { 0 }));

    const std::string path = wc_expected_name(0, 0);

    // Start from a clean slate so a stale file cannot mask a failure.
    wc_remove_file(path);
    ASSERT_FALSE(wc_file_exists(path))
        << "Pre-condition: " << path << " should not exist before the call";

    std::cerr << "  Calling write_crds()...\n";
    write_crds(complexList, bases);

    // 1) The file must have been created with the documented name.
    EXPECT_TRUE(wc_file_exists(path))
        << "write_crds() should have created " << path;

    // 2) The header line must be "molTypeIndex myComIndex".
    std::ifstream in(path);
    std::string firstLine;
    if (in.good())
        std::getline(in, firstLine);
    std::cerr << "  First line read back: \"" << firstLine << "\"\n";

    std::istringstream headerStream(firstLine);
    int readMolTypeIndex { -999 };
    int readMyComIndex { -999 };
    headerStream >> readMolTypeIndex >> readMyComIndex;

    EXPECT_EQ(readMolTypeIndex, 3)
        << "First header field should be Molecule::molTypeIndex (3)";
    EXPECT_EQ(readMyComIndex, 0)
        << "Second header field should be Molecule::myComIndex (0)";

    // 3) write_crd_file() must have appended at least one more line.
    const int lineCount = wc_count_lines(path);
    std::cerr << "  Total lines in file: " << lineCount << '\n';
    EXPECT_GT(lineCount, 1)
        << "The coordinate block written by Molecule::write_crd_file() should "
           "add lines after the header";

    // The file should not be empty in any case.
    EXPECT_FALSE(wc_read_file(path).empty()) << "Dump file should not be empty";

    // Clean up the artifact created by this test.
    wc_remove_file(path);
}

// -----------------------------------------------------------------------------
// Test 2: several complexes, several members each.
//         Verifies that every (complex position, molecule index) pair yields a
//         distinct file and that the header of each file matches its molecule.
// -----------------------------------------------------------------------------
void test_wc_multiple_complexes_and_members()
{
    std::cerr << "\n[TEST] test_wc_multiple_complexes_and_members\n"
              << "  Source file:   src/io/write_crds.cpp\n"
              << "  Function:      write_crds\n"
              << "  Scenario:      two Complexes; the first owns molecules 0 and\n"
              << "                 2, the second owns molecule 1.\n"
              << "  Pass criteria: one file per (complex position, member) pair\n"
              << "                 with the correct header for that molecule.\n";

    wc_make_out_dir();

    // Three molecules with distinct type/complex indices so the headers differ.
    std::vector<Molecule> bases;
    bases.push_back(wc_make_molecule(/*molTypeIndex=*/0, /*myComIndex=*/0,
        Coord { 0.0, 0.0, 0.0 })); // molecule 0 -> complex 0
    bases.push_back(wc_make_molecule(/*molTypeIndex=*/1, /*myComIndex=*/1,
        Coord { 10.0, 0.0, 0.0 })); // molecule 1 -> complex 1
    bases.push_back(wc_make_molecule(/*molTypeIndex=*/2, /*myComIndex=*/0,
        Coord { 0.0, 10.0, 0.0 })); // molecule 2 -> complex 0

    // Complex at position 0 holds molecules {0, 2}; position 1 holds {1}.
    std::vector<Complex> complexList;
    complexList.push_back(wc_make_complex(/*index=*/0, { 0, 2 }));
    complexList.push_back(wc_make_complex(/*index=*/1, { 1 }));

    // Expected files: complex position paired with each member index.
    const std::string p00 = wc_expected_name(0, 0);
    const std::string p02 = wc_expected_name(0, 2);
    const std::string p11 = wc_expected_name(1, 1);

    wc_remove_file(p00);
    wc_remove_file(p02);
    wc_remove_file(p11);

    std::cerr << "  Calling write_crds()...\n";
    write_crds(complexList, bases);

    // Every expected file must exist.
    EXPECT_TRUE(wc_file_exists(p00)) << "Missing expected file " << p00;
    EXPECT_TRUE(wc_file_exists(p02)) << "Missing expected file " << p02;
    EXPECT_TRUE(wc_file_exists(p11)) << "Missing expected file " << p11;

    // Small lambda to check the header of one dump file.
    auto checkHeader = [](const std::string& path, int expTypeIndex,
                           int expComIndex) {
        std::ifstream in(path);
        std::string line;
        if (in.good())
            std::getline(in, line);
        std::cerr << "  " << path << " header: \"" << line << "\"\n";
        std::istringstream iss(line);
        int t { -999 };
        int c { -999 };
        iss >> t >> c;
        EXPECT_EQ(t, expTypeIndex)
            << path << " first field should be molTypeIndex " << expTypeIndex;
        EXPECT_EQ(c, expComIndex)
            << path << " second field should be myComIndex " << expComIndex;
    };

    // Headers must reflect the *molecule* the file belongs to.
    checkHeader(p00, 0, 0);
    checkHeader(p02, 2, 0);
    checkHeader(p11, 1, 1);

    // No file should be produced for a pairing that does not exist, e.g.
    // molecule 1 was never a member of the complex at position 0.
    const std::string bogus = wc_expected_name(0, 1);
    std::cerr << "  Checking that unrelated pairing " << bogus
              << " was NOT created\n";
    EXPECT_FALSE(wc_file_exists(bogus))
        << bogus << " should not exist: molecule 1 is not in complex 0";

    wc_remove_file(p00);
    wc_remove_file(p02);
    wc_remove_file(p11);
}

// -----------------------------------------------------------------------------
// Test 3: the loop uses the positional index, not Complex::index.
//         A single complex placed at position 0 but carrying index = 7 must
//         still produce "c0_..." and never "c7_...".
// -----------------------------------------------------------------------------
void test_wc_uses_positional_complex_index()
{
    std::cerr << "\n[TEST] test_wc_uses_positional_complex_index\n"
              << "  Source file:   src/io/write_crds.cpp\n"
              << "  Function:      write_crds\n"
              << "  Scenario:      Complex::index (7) differs from its position\n"
              << "                 in the vector (0).\n"
              << "  Pass criteria: out/c0_p0_error.dat is written and\n"
              << "                 out/c7_p0_error.dat is NOT written.\n";

    wc_make_out_dir();

    std::vector<Molecule> bases;
    bases.push_back(wc_make_molecule(/*molTypeIndex=*/5, /*myComIndex=*/7,
        Coord { -1.0, -1.0, -1.0 }));

    // Position 0 in the vector, but Complex::index deliberately set to 7.
    std::vector<Complex> complexList;
    complexList.push_back(wc_make_complex(/*index=*/7, { 0 }));

    const std::string positional = wc_expected_name(0, 0); // expected
    const std::string byIndex = wc_expected_name(7, 0);    // must not appear

    wc_remove_file(positional);
    wc_remove_file(byIndex);

    std::cerr << "  Calling write_crds()...\n";
    write_crds(complexList, bases);

    EXPECT_TRUE(wc_file_exists(positional))
        << "File name should be built from the loop position (expected "
        << positional << ")";
    EXPECT_FALSE(wc_file_exists(byIndex))
        << "File name must not be built from Complex::index (" << byIndex
        << " should not exist)";

    // The header should still carry the molecule's own bookkeeping indices.
    std::ifstream in(positional);
    std::string line;
    if (in.good())
        std::getline(in, line);
    std::cerr << "  Header of " << positional << ": \"" << line << "\"\n";
    std::istringstream iss(line);
    int t { -999 };
    int c { -999 };
    iss >> t >> c;
    EXPECT_EQ(t, 5) << "molTypeIndex should be 5";
    EXPECT_EQ(c, 7) << "myComIndex should be 7 (the molecule's own value)";

    wc_remove_file(positional);
    wc_remove_file(byIndex);
}

// -----------------------------------------------------------------------------
// Test 4: degenerate inputs.
//         a) An empty complex list must do nothing and must not crash.
//         b) A complex with an empty memberList must produce no file.
// -----------------------------------------------------------------------------
void test_wc_degenerate_inputs()
{
    std::cerr << "\n[TEST] test_wc_degenerate_inputs\n"
              << "  Source file:   src/io/write_crds.cpp\n"
              << "  Function:      write_crds\n"
              << "  Scenario a:    empty Complexlist and empty bases list.\n"
              << "  Scenario b:    one Complex with an empty memberList.\n"
              << "  Pass criteria: no files are created and the call returns\n"
              << "                 normally (no crash).\n";

    wc_make_out_dir();

    // ---- Scenario a: nothing at all to write. -------------------------------
    std::vector<Complex> emptyComplexList;
    std::vector<Molecule> emptyBases;

    const std::string neverWritten = wc_expected_name(0, 0);
    wc_remove_file(neverWritten);

    std::cerr << "  Calling write_crds() with empty vectors...\n";
    write_crds(emptyComplexList, emptyBases);

    EXPECT_FALSE(wc_file_exists(neverWritten))
        << "No file should be produced from an empty complex list";

    // ---- Scenario b: a complex that owns no molecules. ----------------------
    std::vector<Molecule> bases;
    bases.push_back(wc_make_molecule(/*molTypeIndex=*/0, /*myComIndex=*/0,
        Coord { 0.0, 0.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(wc_make_complex(/*index=*/0, {})); // empty memberList

    std::cerr << "  Calling write_crds() with an empty memberList...\n";
    write_crds(complexList, bases);

    EXPECT_FALSE(wc_file_exists(neverWritten))
        << "A complex with no members should not produce any dump file";

    std::cerr << "  Degenerate-input calls completed without crashing.\n";

    wc_remove_file(neverWritten);
}

// -----------------------------------------------------------------------------
// Test 5: repeated invocation truncates (rather than appends to) the dump file.
//         std::ofstream opens in truncate mode by default, so calling
//         write_crds() twice must leave the same line count as calling it once.
// -----------------------------------------------------------------------------
void test_wc_overwrites_existing_file()
{
    std::cerr << "\n[TEST] test_wc_overwrites_existing_file\n"
              << "  Source file:   src/io/write_crds.cpp\n"
              << "  Function:      write_crds\n"
              << "  Scenario:      write_crds() is called twice on the same data.\n"
              << "  Pass criteria: the dump file is truncated each time, so the\n"
              << "                 line count after two calls equals the count\n"
              << "                 after one call.\n";

    wc_make_out_dir();

    std::vector<Molecule> bases;
    bases.push_back(wc_make_molecule(/*molTypeIndex=*/1, /*myComIndex=*/0,
        Coord { 2.0, 2.0, 2.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(wc_make_complex(/*index=*/0, { 0 }));

    const std::string path = wc_expected_name(0, 0);
    wc_remove_file(path);

    std::cerr << "  First call to write_crds()...\n";
    write_crds(complexList, bases);
    const int linesAfterFirst = wc_count_lines(path);
    std::cerr << "  Lines after first call:  " << linesAfterFirst << '\n';
    EXPECT_GT(linesAfterFirst, 0) << "First call should produce content";

    std::cerr << "  Second call to write_crds()...\n";
    write_crds(complexList, bases);
    const int linesAfterSecond = wc_count_lines(path);
    std::cerr << "  Lines after second call: " << linesAfterSecond << '\n';

    EXPECT_EQ(linesAfterSecond, linesAfterFirst)
        << "The file should be overwritten (truncated), not appended to";

    wc_remove_file(path);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each scenario runs in its own TEST so a failure in one
// does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(WriteCrdsTest, SingleComplexSingleMolecule) { test_wc_single_complex_single_molecule(); }
TEST(WriteCrdsTest, MultipleComplexesAndMembers) { test_wc_multiple_complexes_and_members(); }
TEST(WriteCrdsTest, UsesPositionalComplexIndex) { test_wc_uses_positional_complex_index(); }
TEST(WriteCrdsTest, DegenerateInputs) { test_wc_degenerate_inputs(); }
TEST(WriteCrdsTest, OverwritesExistingFile) { test_wc_overwrites_existing_file(); }