/*! \file test_write_timestep_information.cpp
 *
 * ### Unit test for src/io/write_timestep_information.cpp
 *
 * The function under test is:
 *
 *     void write_timestep_information(long long int simItr,
 *                                     std::ofstream& outFile,
 *                                     std::ofstream& molecTypesFile,
 *                                     std::ofstream& textTimeStatFile,
 *                                     const Parameters& params,
 *                                     std::vector<std::vector<int>>& molecTypesList,
 *                                     std::vector<Molecule>& moleculeList,
 *                                     std::vector<Complex>& complexList,
 *                                     std::vector<MolTemplate>& molTemplateList)
 *
 * What the routine does (and therefore what we verify):
 *   1. Resets the static counter Complex::currNumberComTypes and re-counts how
 *      many *distinct* complex "types" exist in complexList.  A type is
 *      identified by a product-of-(value+1) hash stored in adjmat[i][0].
 *   2. Grows the static counter Complex::currNumberMolTypes and appends one line
 *      of the form "P<molTypeIndex>:<count>" to molecTypesFile for every complex
 *      type that has never been seen before (across all calls).
 *   3. Writes one line to outFile: "<simItr>\t<simTime>\t<count of type 0>\t...".
 *   4. Writes the composition of every complex type to textTimeStatFile.
 *   5. Records new complex types in molecTypesList (caller supplied storage).
 *
 * NOTE on a known quirk of the implementation: while looping over the members of
 * complex i, the code inspects `moleculeList[i].molTypeIndex` (the complex index)
 * instead of the member molecule index.  The tests below are written against the
 * *actual* behaviour, so a single complex is always reported as containing
 * memberList.size() copies of moleculeList[complexIndex].molTypeIndex.  This is
 * called out explicitly in the console output so the reader is not surprised.
 *
 * All output is captured by writing to temporary files, closing them, and then
 * reading the bytes back for comparison.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// Temporary file names used to capture the three output streams.  Unique prefix
// keeps them from colliding with any other test in the suite.
const char* kWtsiOutFile = "wtsi_timestep_out.tmp";
const char* kWtsiMolecTypesFile = "wtsi_molec_types.tmp";
const char* kWtsiTextStatFile = "wtsi_text_stat.tmp";

/*! \brief Slurps an entire text file into a std::string. */
std::string wtsi_read_file(const char* path)
{
    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

/*! \brief Builds the minimal Molecule the routine actually reads (molTypeIndex). */
Molecule wtsi_make_molecule(int index, int molTypeIndex)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = index; // one molecule per complex in most of the tests
    return mol;
}

/*! \brief Builds the minimal Complex the routine actually reads (memberList). */
Complex wtsi_make_complex(int index, const std::vector<int>& memberList)
{
    Complex com;
    com.index = index;
    com.memberList = memberList;
    return com;
}

/*! \brief Creates numMolTypes default MolTemplates (only .size() is used). */
std::vector<MolTemplate> wtsi_make_mol_templates(int numMolTypes)
{
    std::vector<MolTemplate> templates;
    for (int i = 0; i < numMolTypes; ++i) {
        MolTemplate temp;
        temp.molTypeIndex = i;
        temp.molName = "P" + std::to_string(i);
        templates.push_back(temp);
    }
    return templates;
}

/*! \brief Allocates the caller-owned molecTypesList scratch storage.
 *
 * The routine indexes molecTypesList[Complex::currNumberMolTypes][k] with
 * k < 2*numUniqueMols+2, so we allocate the same width the routine uses for its
 * internal adjacency matrix (2 + 2*numMolTypes) and plenty of rows.
 */
std::vector<std::vector<int>> wtsi_make_molec_types_list(int numMolTypes, int rows)
{
    const int width = 2 + 2 * numMolTypes;
    return std::vector<std::vector<int>>(rows, std::vector<int>(width, 0));
}

/*! \brief Captured contents of the three output files after one call. */
struct WtsiOutputs {
    std::string timestepLine; //!< contents of outFile
    std::string molecTypes; //!< contents of molecTypesFile
    std::string textStat; //!< contents of textTimeStatFile
};

/*! \brief Opens the three temp files, invokes the function, and reads them back.
 *
 * The files are truncated on every call, so each invocation reports only what
 * that invocation wrote.
 */
WtsiOutputs wtsi_run(long long int simItr, const Parameters& params,
    std::vector<std::vector<int>>& molecTypesList, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList)
{
    std::ofstream outFile(kWtsiOutFile, std::ios::trunc);
    std::ofstream molecTypesFile(kWtsiMolecTypesFile, std::ios::trunc);
    std::ofstream textTimeStatFile(kWtsiTextStatFile, std::ios::trunc);

    write_timestep_information(simItr, outFile, molecTypesFile, textTimeStatFile, params,
        molecTypesList, moleculeList, complexList, molTemplateList);

    // Flush/close so the bytes are visible to the reader below.
    outFile.close();
    molecTypesFile.close();
    textTimeStatFile.close();

    WtsiOutputs captured;
    captured.timestepLine = wtsi_read_file(kWtsiOutFile);
    captured.molecTypes = wtsi_read_file(kWtsiMolecTypesFile);
    captured.textStat = wtsi_read_file(kWtsiTextStatFile);
    return captured;
}

/*! \brief Removes the temporary files created by wtsi_run(). */
void wtsi_cleanup_files()
{
    std::remove(kWtsiOutFile);
    std::remove(kWtsiMolecTypesFile);
    std::remove(kWtsiTextStatFile);
}

/*! \brief RAII guard for the static Complex counters the routine depends on.
 *
 * write_timestep_information() both reads and writes Complex::numberOfComplexes,
 * Complex::currNumberComTypes and Complex::currNumberMolTypes.  This guard puts
 * them into a well-defined state for the test and restores the previous values
 * afterwards so no other test in the suite is disturbed.
 */
class WtsiStaticGuard {
public:
    explicit WtsiStaticGuard(int numComplexes)
        : savedNumberOfComplexes_(Complex::numberOfComplexes)
        , savedComTypes_(Complex::currNumberComTypes)
        , savedMolTypes_(Complex::currNumberMolTypes)
    {
        Complex::numberOfComplexes = numComplexes; // must match complexList.size()
        Complex::currNumberComTypes = 0;
        Complex::currNumberMolTypes = 0;
    }

    ~WtsiStaticGuard()
    {
        Complex::numberOfComplexes = savedNumberOfComplexes_;
        Complex::currNumberComTypes = savedComTypes_;
        Complex::currNumberMolTypes = savedMolTypes_;
    }

private:
    int savedNumberOfComplexes_;
    int savedComTypes_;
    int savedMolTypes_;
};

/*! \brief Dumps the captured output with escaped tabs/newlines for readability. */
std::string wtsi_escape(const std::string& raw)
{
    std::string out;
    for (char c : raw) {
        if (c == '\t')
            out += "\\t";
        else if (c == '\n')
            out += "\\n";
        else
            out += c;
    }
    return out;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: two complexes with *different* molecule types.
//         Expect two distinct complex types, one molecTypes line each, and a
//         timestep line reporting a population of 1 for each type.
// -----------------------------------------------------------------------------
void test_wtsi_two_distinct_complex_types()
{
    std::cerr << "\n[TEST] test_wtsi_two_distinct_complex_types\n"
              << "  Source file: src/io/write_timestep_information.cpp\n"
              << "  Function:    write_timestep_information()\n"
              << "  Scenario:    2 complexes, each with 1 member; complex 0 maps to\n"
              << "               molecule type 0 and complex 1 to molecule type 1.\n"
              << "  Pass:        currNumberComTypes == 2, currNumberMolTypes == 2,\n"
              << "               molecTypes file lists 'P0:1' and 'P1:1', and the\n"
              << "               timestep line reports one complex of each type.\n";

    // Simulation parameters: timeStep 0 keeps the printed simulation time simple.
    Parameters params;
    params.timeStep = 0.1;

    std::vector<MolTemplate> molTemplateList = wtsi_make_mol_templates(2);
    std::vector<Molecule> moleculeList { wtsi_make_molecule(0, 0), wtsi_make_molecule(1, 1) };
    std::vector<Complex> complexList { wtsi_make_complex(0, { 0 }), wtsi_make_complex(1, { 1 }) };
    std::vector<std::vector<int>> molecTypesList = wtsi_make_molec_types_list(2, 16);

    WtsiStaticGuard guard(static_cast<int>(complexList.size()));

    // simItr == 0 so the printed simulation time is exactly "0".
    WtsiOutputs captured = wtsi_run(0, params, molecTypesList, moleculeList, complexList, molTemplateList);

    std::cerr << "  outFile        = \"" << wtsi_escape(captured.timestepLine) << "\"\n"
              << "  molecTypesFile = \"" << wtsi_escape(captured.molecTypes) << "\"\n"
              << "  textStatFile   = \"" << wtsi_escape(captured.textStat) << "\"\n"
              << "  currNumberComTypes = " << Complex::currNumberComTypes
              << ", currNumberMolTypes = " << Complex::currNumberMolTypes << "\n";

    // Two structurally different complexes -> two complex types.
    EXPECT_EQ(Complex::currNumberComTypes, 2)
        << "Two complexes with different molecule types must yield 2 complex types";
    // Both types are brand new, so both are recorded as new molecule types.
    EXPECT_EQ(Complex::currNumberMolTypes, 2)
        << "Both freshly seen complex types must be appended to the molecule type list";

    // Timestep line: iteration, simulation time, then one population per type.
    EXPECT_EQ(captured.timestepLine, "0\t0\t1\t1\t\n")
        << "Timestep line should be '<itr>\\t<time>\\t<pop type0>\\t<pop type1>\\t'";

    // Species file: one "P<type>:<count>" line per newly observed complex type.
    EXPECT_EQ(captured.molecTypes, "P0:1\nP1:1\n")
        << "molecTypesFile should list one composition line per new complex type";

    // Human readable composition file, tab separated, terminated by std::endl.
    EXPECT_EQ(captured.textStat, "P0:1\tP1:1\t\n")
        << "textTimeStatFile should list the composition of every complex type";

    // molecTypesList must have been populated with the complex-type hashes.
    EXPECT_NE(molecTypesList[0][0], 0) << "Row 0 of molecTypesList should hold a complex-type id";
    EXPECT_NE(molecTypesList[1][0], 0) << "Row 1 of molecTypesList should hold a complex-type id";
    EXPECT_NE(molecTypesList[0][0], molecTypesList[1][0])
        << "The two recorded complex types must have different ids";

    wtsi_cleanup_files();
}

// -----------------------------------------------------------------------------
// Test 2: two *identical* complexes.
//         Expect a single complex type whose population is reported as 2.
// -----------------------------------------------------------------------------
void test_wtsi_duplicate_complexes_are_merged()
{
    std::cerr << "\n[TEST] test_wtsi_duplicate_complexes_are_merged\n"
              << "  Source file: src/io/write_timestep_information.cpp\n"
              << "  Function:    write_timestep_information()\n"
              << "  Scenario:    2 complexes that both map to molecule type 0.\n"
              << "  Pass:        currNumberComTypes == 1, currNumberMolTypes == 1,\n"
              << "               and the single reported population equals 2.\n";

    Parameters params;
    params.timeStep = 0.1;

    std::vector<MolTemplate> molTemplateList = wtsi_make_mol_templates(2);
    // Both molecules share molTypeIndex 0 => both complexes hash identically.
    std::vector<Molecule> moleculeList { wtsi_make_molecule(0, 0), wtsi_make_molecule(1, 0) };
    std::vector<Complex> complexList { wtsi_make_complex(0, { 0 }), wtsi_make_complex(1, { 1 }) };
    std::vector<std::vector<int>> molecTypesList = wtsi_make_molec_types_list(2, 16);

    WtsiStaticGuard guard(static_cast<int>(complexList.size()));

    WtsiOutputs captured = wtsi_run(0, params, molecTypesList, moleculeList, complexList, molTemplateList);

    std::cerr << "  outFile        = \"" << wtsi_escape(captured.timestepLine) << "\"\n"
              << "  molecTypesFile = \"" << wtsi_escape(captured.molecTypes) << "\"\n"
              << "  textStatFile   = \"" << wtsi_escape(captured.textStat) << "\"\n"
              << "  currNumberComTypes = " << Complex::currNumberComTypes
              << ", currNumberMolTypes = " << Complex::currNumberMolTypes << "\n";

    EXPECT_EQ(Complex::currNumberComTypes, 1)
        << "Identical complexes must collapse into a single complex type";
    EXPECT_EQ(Complex::currNumberMolTypes, 1)
        << "Only one new complex type should be appended to the molecule type list";

    // Only one column, and its population is 2 (both complexes counted).
    EXPECT_EQ(captured.timestepLine, "0\t0\t2\t\n")
        << "The single complex type should have a population of 2";
    EXPECT_EQ(captured.molecTypes, "P0:1\n")
        << "Only the one distinct complex composition should be written";
    EXPECT_EQ(captured.textStat, "P0:1\t\n")
        << "Only the one distinct complex composition should be written";

    wtsi_cleanup_files();
}

// -----------------------------------------------------------------------------
// Test 3: a single complex holding three members.
//         Expect the composition to report a count of 3 for the molecule type.
// -----------------------------------------------------------------------------
void test_wtsi_multi_member_complex_counts_members()
{
    std::cerr << "\n[TEST] test_wtsi_multi_member_complex_counts_members\n"
              << "  Source file: src/io/write_timestep_information.cpp\n"
              << "  Function:    write_timestep_information()\n"
              << "  Scenario:    1 complex whose memberList has 3 entries.\n"
              << "  Note:        the implementation reads moleculeList[complexIndex]\n"
              << "               for every member, so all 3 members are attributed to\n"
              << "               molecule type of moleculeList[0].\n"
              << "  Pass:        composition strings report 'P0:3' and the population\n"
              << "               of the single complex type is 1.\n";

    Parameters params;
    params.timeStep = 0.1;

    std::vector<MolTemplate> molTemplateList = wtsi_make_mol_templates(2);
    // Three molecules, all of type 0; the routine only ever looks at index 0 here.
    std::vector<Molecule> moleculeList { wtsi_make_molecule(0, 0), wtsi_make_molecule(1, 0),
        wtsi_make_molecule(2, 0) };
    std::vector<Complex> complexList { wtsi_make_complex(0, { 0, 1, 2 }) };
    std::vector<std::vector<int>> molecTypesList = wtsi_make_molec_types_list(2, 16);

    WtsiStaticGuard guard(static_cast<int>(complexList.size()));

    WtsiOutputs captured = wtsi_run(0, params, molecTypesList, moleculeList, complexList, molTemplateList);

    std::cerr << "  outFile        = \"" << wtsi_escape(captured.timestepLine) << "\"\n"
              << "  molecTypesFile = \"" << wtsi_escape(captured.molecTypes) << "\"\n"
              << "  textStatFile   = \"" << wtsi_escape(captured.textStat) << "\"\n"
              << "  currNumberComTypes = " << Complex::currNumberComTypes
              << ", currNumberMolTypes = " << Complex::currNumberMolTypes << "\n";

    EXPECT_EQ(Complex::currNumberComTypes, 1) << "A single complex yields exactly one complex type";
    EXPECT_EQ(Complex::currNumberMolTypes, 1) << "That one complex type is new and must be recorded";

    // The three members are collapsed into "P0:3".
    EXPECT_EQ(captured.molecTypes, "P0:3\n")
        << "A complex with three members of type 0 should be written as 'P0:3'";
    EXPECT_EQ(captured.textStat, "P0:3\t\n")
        << "The human readable composition should also read 'P0:3'";
    EXPECT_EQ(captured.timestepLine, "0\t0\t1\t\n")
        << "There is only one complex of this type in the system";

    // The recorded row should say: <hash>, 1 unique mol type, type index 0, count 3.
    EXPECT_EQ(molecTypesList[0][1], 1) << "Row should record exactly one unique molecule type";
    EXPECT_EQ(molecTypesList[0][2], 0) << "Row should record molecule type index 0";
    EXPECT_EQ(molecTypesList[0][3], 3) << "Row should record a member count of 3";

    wtsi_cleanup_files();
}

// -----------------------------------------------------------------------------
// Test 4: non-zero iteration / timestep -> the printed simulation time must be
//         simItr * timeStep * 1e-6 (microseconds -> seconds).
// -----------------------------------------------------------------------------
void test_wtsi_writes_iteration_and_time()
{
    std::cerr << "\n[TEST] test_wtsi_writes_iteration_and_time\n"
              << "  Source file: src/io/write_timestep_information.cpp\n"
              << "  Function:    write_timestep_information()\n"
              << "  Scenario:    simItr = 1000000, timeStep = 1.0 us.\n"
              << "  Pass:        the timestep line begins with the iteration number\n"
              << "               followed by simItr*timeStep*1e-6 == 1 second.\n";

    Parameters params;
    params.timeStep = 1.0; // microseconds

    std::vector<MolTemplate> molTemplateList = wtsi_make_mol_templates(1);
    std::vector<Molecule> moleculeList { wtsi_make_molecule(0, 0) };
    std::vector<Complex> complexList { wtsi_make_complex(0, { 0 }) };
    std::vector<std::vector<int>> molecTypesList = wtsi_make_molec_types_list(1, 16);

    WtsiStaticGuard guard(static_cast<int>(complexList.size()));

    const long long int simItr = 1000000; // 1e6 * 1 us = 1 second
    WtsiOutputs captured = wtsi_run(simItr, params, molecTypesList, moleculeList, complexList, molTemplateList);

    std::cerr << "  outFile = \"" << wtsi_escape(captured.timestepLine) << "\"\n";

    // Iteration number must be the first field.
    EXPECT_EQ(captured.timestepLine.rfind("1000000\t", 0), 0u)
        << "The line must start with the iteration number followed by a tab";
    // Simulation time (seconds) must be the second field: 1e6 * 1.0 * 1e-6 = 1.
    EXPECT_EQ(captured.timestepLine, "1000000\t1\t1\t\n")
        << "Expected '<itr>\\t<time in s>\\t<population>\\t'";

    wtsi_cleanup_files();
}

// -----------------------------------------------------------------------------
// Test 5: calling the routine twice.
//         Complex::currNumberComTypes is recomputed from scratch each call while
//         Complex::currNumberMolTypes accumulates only genuinely new types, so a
//         repeated call must not append anything to molecTypesFile.
// -----------------------------------------------------------------------------
void test_wtsi_repeated_call_does_not_duplicate_types()
{
    std::cerr << "\n[TEST] test_wtsi_repeated_call_does_not_duplicate_types\n"
              << "  Source file: src/io/write_timestep_information.cpp\n"
              << "  Function:    write_timestep_information()\n"
              << "  Scenario:    the same system is written twice in a row.\n"
              << "  Pass:        currNumberComTypes is 2 after both calls, the static\n"
              << "               molecule-type counter stays at 2, and the second call\n"
              << "               writes nothing to molecTypesFile.\n";

    Parameters params;
    params.timeStep = 0.1;

    std::vector<MolTemplate> molTemplateList = wtsi_make_mol_templates(2);
    std::vector<Molecule> moleculeList { wtsi_make_molecule(0, 0), wtsi_make_molecule(1, 1) };
    std::vector<Complex> complexList { wtsi_make_complex(0, { 0 }), wtsi_make_complex(1, { 1 }) };
    // The SAME molecTypesList is reused across both calls, which is what allows the
    // routine to recognize already-known complex types.
    std::vector<std::vector<int>> molecTypesList = wtsi_make_molec_types_list(2, 16);

    WtsiStaticGuard guard(static_cast<int>(complexList.size()));

    // ---- first call: everything is new ------------------------------------
    WtsiOutputs first = wtsi_run(0, params, molecTypesList, moleculeList, complexList, molTemplateList);
    std::cerr << "  after call 1: comTypes = " << Complex::currNumberComTypes
              << ", molTypes = " << Complex::currNumberMolTypes
              << ", molecTypesFile = \"" << wtsi_escape(first.molecTypes) << "\"\n";
    EXPECT_EQ(Complex::currNumberComTypes, 2) << "First call should find two complex types";
    EXPECT_EQ(Complex::currNumberMolTypes, 2) << "First call should record two new molecule types";
    EXPECT_EQ(first.molecTypes, "P0:1\nP1:1\n") << "First call writes both compositions";

    // ---- second call: nothing is new -------------------------------------
    WtsiOutputs second = wtsi_run(1, params, molecTypesList, moleculeList, complexList, molTemplateList);
    std::cerr << "  after call 2: comTypes = " << Complex::currNumberComTypes
              << ", molTypes = " << Complex::currNumberMolTypes
              << ", molecTypesFile = \"" << wtsi_escape(second.molecTypes) << "\"\n";

    EXPECT_EQ(Complex::currNumberComTypes, 2)
        << "Complex type count is recomputed each call and must stay at 2";
    EXPECT_EQ(Complex::currNumberMolTypes, 2)
        << "No new complex types appeared, so the static molecule-type counter must not grow";
    EXPECT_TRUE(second.molecTypes.empty())
        << "Second call must not re-write already known compositions, got: \""
        << wtsi_escape(second.molecTypes) << "\"";

    // Populations are unchanged; only the iteration/time fields differ.
    EXPECT_EQ(second.timestepLine, "1\t1e-07\t1\t1\t\n")
        << "Second call should report iteration 1 with the same populations";

    wtsi_cleanup_files();
}

// -----------------------------------------------------------------------------
// Test 6: robustness - a complex with an empty memberList.
//         The routine must not crash and must still emit a (degenerate) record.
// -----------------------------------------------------------------------------
void test_wtsi_empty_member_list_is_handled()
{
    std::cerr << "\n[TEST] test_wtsi_empty_member_list_is_handled\n"
              << "  Source file: src/io/write_timestep_information.cpp\n"
              << "  Function:    write_timestep_information()\n"
              << "  Scenario:    a single complex with no members at all.\n"
              << "  Pass:        no crash; one complex type with zero unique molecule\n"
              << "               types, population 1, and an empty composition string.\n";

    Parameters params;
    params.timeStep = 0.1;

    std::vector<MolTemplate> molTemplateList = wtsi_make_mol_templates(2);
    std::vector<Molecule> moleculeList { wtsi_make_molecule(0, 0) };
    std::vector<Complex> complexList { wtsi_make_complex(0, {}) }; // no members
    std::vector<std::vector<int>> molecTypesList = wtsi_make_molec_types_list(2, 16);

    WtsiStaticGuard guard(static_cast<int>(complexList.size()));

    WtsiOutputs captured = wtsi_run(0, params, molecTypesList, moleculeList, complexList, molTemplateList);

    std::cerr << "  outFile        = \"" << wtsi_escape(captured.timestepLine) << "\"\n"
              << "  molecTypesFile = \"" << wtsi_escape(captured.molecTypes) << "\"\n"
              << "  textStatFile   = \"" << wtsi_escape(captured.textStat) << "\"\n"
              << "  currNumberComTypes = " << Complex::currNumberComTypes
              << ", currNumberMolTypes = " << Complex::currNumberMolTypes << "\n";

    EXPECT_EQ(Complex::currNumberComTypes, 1)
        << "Even an empty complex counts as exactly one complex type";
    EXPECT_EQ(Complex::currNumberMolTypes, 1)
        << "The empty complex type is new and must be recorded once";
    EXPECT_EQ(molecTypesList[0][1], 0)
        << "The recorded row should report zero unique molecule types";

    // Population line still lists the single (empty) complex.
    EXPECT_EQ(captured.timestepLine, "0\t0\t1\t\n")
        << "One complex of the (empty) type should be reported";
    // No "P<type>:<count>" tokens, just the terminating newline / tab.
    EXPECT_EQ(captured.molecTypes, "\n")
        << "An empty complex has no composition tokens, only the line break";
    EXPECT_EQ(captured.textStat, "\t\n")
        << "An empty complex contributes only the field separator";

    wtsi_cleanup_files();
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs in its own TEST so that a failure
// in one scenario does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(WriteTimestepInformation, TwoDistinctComplexTypes) { test_wtsi_two_distinct_complex_types(); }
TEST(WriteTimestepInformation, DuplicateComplexesAreMerged) { test_wtsi_duplicate_complexes_are_merged(); }
TEST(WriteTimestepInformation, MultiMemberComplexCountsMembers) { test_wtsi_multi_member_complex_counts_members(); }
TEST(WriteTimestepInformation, WritesIterationAndTime) { test_wtsi_writes_iteration_and_time(); }
TEST(WriteTimestepInformation, RepeatedCallDoesNotDuplicateTypes) { test_wtsi_repeated_call_does_not_duplicate_types(); }
TEST(WriteTimestepInformation, EmptyMemberListIsHandled) { test_wtsi_empty_member_list_is_handled(); }