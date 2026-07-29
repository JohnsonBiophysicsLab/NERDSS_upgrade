/*! \file test_write_complex_crds.cpp
 *
 * ### Unit test for ../src/io/write_complex_crds.cpp
 *
 * This test exercises the single function defined in that file:
 *
 *     void write_complex_crds(std::string name,
 *                             const Complex& complex1,
 *                             const Complex& complex2,
 *                             std::vector<Molecule>& moleculeList)
 *
 * Behaviour under test
 * --------------------
 * For every member Molecule index `mp` in `complex1.memberList` and in
 * `complex2.memberList` the function opens
 *
 *     out/c<complex.index>_p<mp>_<name>.dat
 *
 * writes a single header line containing
 *
 *     <moleculeList[mp].molTypeIndex> <moleculeList[mp].myComIndex>
 *
 * and then delegates to `Molecule::write_crd_file()` to dump the coordinates.
 *
 * Consequently the observable behaviour we can verify is:
 *   1. one file per member molecule of *each* complex, with the exact expected
 *      file name (index of the complex, index of the molecule, and the passed
 *      `name` tag),
 *   2. the header line holds the molTypeIndex / myComIndex of the molecule
 *      referenced by the member list (i.e. the *molecule list index* is used
 *      to look up the data, not the loop counter),
 *   3. coordinate data is appended after the header (file is longer than the
 *      header alone),
 *   4. a complex with an empty member list produces no files,
 *   5. different `name` arguments produce different files.
 *
 * NOTE: the function writes into the relative directory "out/". If that
 * directory does not exist the std::ofstream silently fails and nothing is
 * written, so every test first makes sure the directory exists and reports
 * loudly if it could not be created.
 *
 * Verbose progress information is printed to stderr so a reader of the test
 * log can see which source file / function is exercised and what each
 * assertion checks.
 */

#include <sys/stat.h>
#include <sys/types.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "io/io.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Small local helpers (anonymous namespace so they cannot collide with helpers
// from other translation units in the combined test binary).
// -----------------------------------------------------------------------------
namespace {

//! Directory the function under test hard-codes as its output location.
const char* kWccOutDir = "out";

/*! \brief Make sure the "out" directory exists (create it if needed).
 *  \return true if the directory exists (or was created) after the call.
 */
bool wcc_ensure_out_dir()
{
    struct stat info {};
    if (stat(kWccOutDir, &info) == 0) {
        // Path exists -- make sure it is really a directory.
        return (info.st_mode & S_IFDIR) != 0;
    }
    // Try to create it with permissive permissions; ignore EEXIST races.
    if (mkdir(kWccOutDir, 0777) == 0)
        return true;
    return stat(kWccOutDir, &info) == 0 && (info.st_mode & S_IFDIR) != 0;
}

/*! \brief Reproduce the exact file name the function under test builds. */
std::string wcc_expected_file(int complexIndex, int molIndex, const std::string& name)
{
    return std::string(kWccOutDir) + "/c" + std::to_string(complexIndex) + "_p"
        + std::to_string(molIndex) + "_" + name + ".dat";
}

/*! \brief Does the given path exist and is it readable? */
bool wcc_file_exists(const std::string& path)
{
    std::ifstream in(path);
    return in.good();
}

/*! \brief Read the whole file into a string ("" if the file cannot be opened). */
std::string wcc_read_file(const std::string& path)
{
    std::ifstream in(path);
    if (!in.good())
        return std::string {};
    std::string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    return contents;
}

/*! \brief Read the first line of a file ("" if it cannot be opened). */
std::string wcc_read_first_line(const std::string& path)
{
    std::ifstream in(path);
    std::string line;
    if (in.good())
        std::getline(in, line);
    return line;
}

/*! \brief Delete a file, ignoring failures (used for test setup/teardown). */
void wcc_remove_file(const std::string& path) { std::remove(path.c_str()); }

/*! \brief Build a minimal, self-consistent Molecule.
 *
 * \param[in] index         index of this Molecule inside moleculeList
 * \param[in] molTypeIndex  MolTemplate index written on the header line
 * \param[in] comIndex      parent complex index written on the header line
 * \param[in] com           center-of-mass coordinate
 * \param[in] numIfaces     how many interfaces to create (offset from the COM)
 */
Molecule wcc_make_molecule(int index, int molTypeIndex, int comIndex, const Coord& com, int numIfaces)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = comIndex;
    mol.comCoord = com;
    mol.mass = 1.0;

    mol.interfaceList.clear();
    for (int i = 0; i < numIfaces; ++i) {
        Molecule::Iface iface;
        // Offset each interface slightly so the coordinates are distinguishable.
        iface.coord = Coord { com.x + 1.0 + i, com.y + 2.0 + i, com.z + 3.0 + i };
        iface.relIndex = i;
        iface.index = i;
        iface.molTypeIndex = molTypeIndex;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Build a Complex with a specified index and member molecule list. */
Complex wcc_make_complex(int index, const std::vector<int>& members)
{
    Complex com;
    com.index = index;
    com.memberList = members;
    com.comCoord = Coord { 0.0, 0.0, 0.0 };
    com.radius = 1.0;
    com.mass = static_cast<double>(members.size());
    return com;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: both complexes have members -> one correctly named file per member.
// -----------------------------------------------------------------------------
void test_wcc_creates_one_file_per_member()
{
    std::cerr << "\n[TEST] test_wcc_creates_one_file_per_member\n"
              << "  Source file:   src/io/write_complex_crds.cpp\n"
              << "  Function:      write_complex_crds()\n"
              << "  Scenario:      complex 0 owns molecules {0,1}, complex 1 owns {2}.\n"
              << "  Pass criteria: exactly the files out/c<ci>_p<mi>_<name>.dat exist.\n";

    // The output directory must exist or the ofstreams silently do nothing.
    const bool haveDir = wcc_ensure_out_dir();
    EXPECT_TRUE(haveDir) << "Could not create/find the 'out' directory required by "
                            "write_complex_crds(); file checks below will fail.";
    if (!haveDir)
        return; // nothing further can be verified

    // Three molecules; complex 0 owns molecules 0 and 1, complex 1 owns molecule 2.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wcc_make_molecule(0, /*molTypeIndex*/ 0, /*comIndex*/ 0, Coord { 0.0, 0.0, 0.0 }, 2));
    moleculeList.push_back(wcc_make_molecule(1, /*molTypeIndex*/ 1, /*comIndex*/ 0, Coord { 5.0, 0.0, 0.0 }, 1));
    moleculeList.push_back(wcc_make_molecule(2, /*molTypeIndex*/ 2, /*comIndex*/ 1, Coord { -5.0, 1.0, 2.0 }, 3));

    Complex complex1 = wcc_make_complex(/*index*/ 0, { 0, 1 });
    Complex complex2 = wcc_make_complex(/*index*/ 1, { 2 });

    const std::string tag = "unittest_basic";

    // Names we expect to be produced.
    const std::string f0 = wcc_expected_file(0, 0, tag);
    const std::string f1 = wcc_expected_file(0, 1, tag);
    const std::string f2 = wcc_expected_file(1, 2, tag);

    // Start from a clean slate so we know the files come from this call.
    wcc_remove_file(f0);
    wcc_remove_file(f1);
    wcc_remove_file(f2);

    std::cerr << "  Calling write_complex_crds(\"" << tag << "\", ...)\n";
    write_complex_crds(tag, complex1, complex2, moleculeList);

    // Each member molecule of each complex must have produced a file.
    std::cerr << "  Checking that " << f0 << " was created\n";
    EXPECT_TRUE(wcc_file_exists(f0)) << "Missing file for complex 0 / molecule 0: " << f0;
    std::cerr << "  Checking that " << f1 << " was created\n";
    EXPECT_TRUE(wcc_file_exists(f1)) << "Missing file for complex 0 / molecule 1: " << f1;
    std::cerr << "  Checking that " << f2 << " was created\n";
    EXPECT_TRUE(wcc_file_exists(f2)) << "Missing file for complex 1 / molecule 2: " << f2;

    // Files that would correspond to wrong (complex, molecule) pairings must NOT
    // exist -- this verifies the naming uses the molecule list index, not the
    // loop position, and that each complex uses its own index.
    const std::string wrongPairing = wcc_expected_file(1, 0, tag); // complex 1 does not own mol 0
    std::cerr << "  Checking that the bogus pairing " << wrongPairing << " was NOT created\n";
    EXPECT_FALSE(wcc_file_exists(wrongPairing))
        << "File named with a (complex, molecule) pairing that does not exist was created: "
        << wrongPairing;

    // Teardown: keep the working tree tidy for subsequent tests.
    wcc_remove_file(f0);
    wcc_remove_file(f1);
    wcc_remove_file(f2);
}

// -----------------------------------------------------------------------------
// Test 2: header line holds molTypeIndex and myComIndex of the referenced mol.
// -----------------------------------------------------------------------------
void test_wcc_header_line_contents()
{
    std::cerr << "\n[TEST] test_wcc_header_line_contents\n"
              << "  Source file:   src/io/write_complex_crds.cpp\n"
              << "  Function:      write_complex_crds()\n"
              << "  Scenario:      molecule 1 has molTypeIndex=7, myComIndex=3.\n"
              << "  Pass criteria: first line of its file is exactly \"7 3\".\n";

    const bool haveDir = wcc_ensure_out_dir();
    EXPECT_TRUE(haveDir) << "Could not create/find the 'out' directory.";
    if (!haveDir)
        return;

    std::vector<Molecule> moleculeList;
    // Molecule 0 is deliberately different so we can be sure the correct entry
    // in moleculeList is consulted.
    moleculeList.push_back(wcc_make_molecule(0, /*molTypeIndex*/ 99, /*comIndex*/ 99, Coord { 0.0, 0.0, 0.0 }, 1));
    moleculeList.push_back(wcc_make_molecule(1, /*molTypeIndex*/ 7, /*comIndex*/ 3, Coord { 1.5, -2.5, 3.5 }, 2));

    // Only molecule 1 is a member; complex2 is intentionally empty here.
    Complex complex1 = wcc_make_complex(/*index*/ 4, { 1 });
    Complex complex2 = wcc_make_complex(/*index*/ 5, {});

    const std::string tag = "unittest_header";
    const std::string path = wcc_expected_file(4, 1, tag);
    wcc_remove_file(path);

    std::cerr << "  Calling write_complex_crds(\"" << tag << "\", ...)\n";
    write_complex_crds(tag, complex1, complex2, moleculeList);

    ASSERT_TRUE(true); // keep going even if the file is missing (checked below)
    EXPECT_TRUE(wcc_file_exists(path)) << "Expected output file was not created: " << path;

    const std::string firstLine = wcc_read_first_line(path);
    std::cerr << "  First line read back: \"" << firstLine << "\"\n";
    EXPECT_EQ(firstLine, std::string("7 3"))
        << "Header line should be '<molTypeIndex> <myComIndex>' of molecule 1";

    wcc_remove_file(path);
}

// -----------------------------------------------------------------------------
// Test 3: coordinate data is written after the header (file longer than header).
// -----------------------------------------------------------------------------
void test_wcc_writes_coordinate_body()
{
    std::cerr << "\n[TEST] test_wcc_writes_coordinate_body\n"
              << "  Source file:   src/io/write_complex_crds.cpp\n"
              << "  Function:      write_complex_crds() -> Molecule::write_crd_file()\n"
              << "  Scenario:      a molecule with 3 interfaces is written out.\n"
              << "  Pass criteria: the file contains more than just the header line.\n";

    const bool haveDir = wcc_ensure_out_dir();
    EXPECT_TRUE(haveDir) << "Could not create/find the 'out' directory.";
    if (!haveDir)
        return;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(
        wcc_make_molecule(0, /*molTypeIndex*/ 0, /*comIndex*/ 0, Coord { 10.0, 20.0, 30.0 }, /*numIfaces*/ 3));

    Complex complex1 = wcc_make_complex(/*index*/ 0, { 0 });
    Complex complex2 = wcc_make_complex(/*index*/ 1, {}); // no members

    const std::string tag = "unittest_body";
    const std::string path = wcc_expected_file(0, 0, tag);
    wcc_remove_file(path);

    std::cerr << "  Calling write_complex_crds(\"" << tag << "\", ...)\n";
    write_complex_crds(tag, complex1, complex2, moleculeList);

    EXPECT_TRUE(wcc_file_exists(path)) << "Expected output file was not created: " << path;

    const std::string contents = wcc_read_file(path);
    // The header the function itself writes is "0 0\n" (4 characters). Anything
    // beyond that must come from Molecule::write_crd_file().
    const std::string header = "0 0\n";
    std::cerr << "  File length = " << contents.size() << " bytes (header alone would be "
              << header.size() << " bytes)\n";
    EXPECT_GT(contents.size(), header.size())
        << "Nothing appears to have been written after the header line; "
           "Molecule::write_crd_file() output is missing";

    // Sanity: the contents should still begin with the header.
    EXPECT_EQ(contents.compare(0, header.size(), header), 0)
        << "File should begin with the '<molTypeIndex> <myComIndex>' header line";

    wcc_remove_file(path);
}

// -----------------------------------------------------------------------------
// Test 4: two empty member lists -> no files created at all.
// -----------------------------------------------------------------------------
void test_wcc_empty_member_lists_write_nothing()
{
    std::cerr << "\n[TEST] test_wcc_empty_member_lists_write_nothing\n"
              << "  Source file:   src/io/write_complex_crds.cpp\n"
              << "  Function:      write_complex_crds()\n"
              << "  Scenario:      both complexes have empty memberLists.\n"
              << "  Pass criteria: no output files are produced and no crash occurs.\n";

    const bool haveDir = wcc_ensure_out_dir();
    EXPECT_TRUE(haveDir) << "Could not create/find the 'out' directory.";
    if (!haveDir)
        return;

    // A molecule exists but is not referenced by any complex.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wcc_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, 1));

    Complex complex1 = wcc_make_complex(/*index*/ 11, {});
    Complex complex2 = wcc_make_complex(/*index*/ 12, {});

    const std::string tag = "unittest_empty";
    const std::string bogus1 = wcc_expected_file(11, 0, tag);
    const std::string bogus2 = wcc_expected_file(12, 0, tag);
    wcc_remove_file(bogus1);
    wcc_remove_file(bogus2);

    std::cerr << "  Calling write_complex_crds(\"" << tag << "\", ...) with empty member lists\n";
    write_complex_crds(tag, complex1, complex2, moleculeList);

    std::cerr << "  Verifying no file " << bogus1 << " was produced\n";
    EXPECT_FALSE(wcc_file_exists(bogus1))
        << "No file should be written for a complex with an empty memberList";
    std::cerr << "  Verifying no file " << bogus2 << " was produced\n";
    EXPECT_FALSE(wcc_file_exists(bogus2))
        << "No file should be written for a complex with an empty memberList";
}

// -----------------------------------------------------------------------------
// Test 5: the `name` argument is part of the file name, so two calls with
//         different tags produce two distinct files for the same molecule.
// -----------------------------------------------------------------------------
void test_wcc_name_tag_distinguishes_files()
{
    std::cerr << "\n[TEST] test_wcc_name_tag_distinguishes_files\n"
              << "  Source file:   src/io/write_complex_crds.cpp\n"
              << "  Function:      write_complex_crds()\n"
              << "  Scenario:      the same complex is dumped twice with different tags.\n"
              << "  Pass criteria: both tag-specific files exist independently.\n";

    const bool haveDir = wcc_ensure_out_dir();
    EXPECT_TRUE(haveDir) << "Could not create/find the 'out' directory.";
    if (!haveDir)
        return;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wcc_make_molecule(0, /*molTypeIndex*/ 3, /*comIndex*/ 2, Coord { 1.0, 1.0, 1.0 }, 1));

    Complex complex1 = wcc_make_complex(/*index*/ 2, { 0 });
    Complex complex2 = wcc_make_complex(/*index*/ 9, {});

    const std::string tagBefore = "unittest_before";
    const std::string tagAfter = "unittest_after";
    const std::string pathBefore = wcc_expected_file(2, 0, tagBefore);
    const std::string pathAfter = wcc_expected_file(2, 0, tagAfter);
    wcc_remove_file(pathBefore);
    wcc_remove_file(pathAfter);

    std::cerr << "  Calling write_complex_crds with tag \"" << tagBefore << "\"\n";
    write_complex_crds(tagBefore, complex1, complex2, moleculeList);
    std::cerr << "  Calling write_complex_crds with tag \"" << tagAfter << "\"\n";
    write_complex_crds(tagAfter, complex1, complex2, moleculeList);

    EXPECT_TRUE(wcc_file_exists(pathBefore)) << "Missing tag-specific file: " << pathBefore;
    EXPECT_TRUE(wcc_file_exists(pathAfter)) << "Missing tag-specific file: " << pathAfter;

    // Both dumps describe the same molecule, so their contents must match.
    const std::string contentsBefore = wcc_read_file(pathBefore);
    const std::string contentsAfter = wcc_read_file(pathAfter);
    std::cerr << "  Comparing the two dumps (" << contentsBefore.size() << " vs "
              << contentsAfter.size() << " bytes)\n";
    EXPECT_EQ(contentsBefore, contentsAfter)
        << "Dumping the same unchanged molecule twice should give identical contents";

    wcc_remove_file(pathBefore);
    wcc_remove_file(pathAfter);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* function runs inside its own TEST so
// individual results are reported while all of them still execute.
// -----------------------------------------------------------------------------
TEST(WriteComplexCrds, CreatesOneFilePerMember) { test_wcc_creates_one_file_per_member(); }
TEST(WriteComplexCrds, HeaderLineContents) { test_wcc_header_line_contents(); }
TEST(WriteComplexCrds, WritesCoordinateBody) { test_wcc_writes_coordinate_body(); }
TEST(WriteComplexCrds, EmptyMemberListsWriteNothing) { test_wcc_empty_member_lists_write_nothing(); }
TEST(WriteComplexCrds, NameTagDistinguishesFiles) { test_wcc_name_tag_distinguishes_files(); }