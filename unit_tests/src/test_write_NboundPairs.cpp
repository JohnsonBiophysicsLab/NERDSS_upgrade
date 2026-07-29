/*! \file test_write_NboundPairs.cpp
 *
 * ### Unit test for src/io/write_NboundPairs.cpp
 *
 * The single function under test is:
 *
 *     void write_NboundPairs(copyCounters& counterArrays, std::ofstream& outfile,
 *                           int it, const Parameters& params,
 *                           std::vector<Molecule>& moleculeList)
 *
 * The function appends exactly one tab-delimited line to `outfile`:
 *
 *   column 0                      : simulation time, computed as
 *                                   (it - params.itrRestartFrom) * params.timeStep * 1E-6
 *                                   + params.timeRestartFrom
 *   columns 1 .. N                 : counterArrays.nBoundPairs[ proPairlist[i] ] for each
 *                                   entry i of counterArrays.proPairlist (N = proPairlist.size())
 *   final 8 columns (in order)     : nLoops, nCancelOverlapPartner, nCancelOverlapSystem,
 *                                   nCancelSpanBox, nCancelDisplace2D, nCancelDisplace3D,
 *                                   nCancelDisplace3Dto2D, nAssocSuccess
 *
 * Since the function has no return value, every test writes to a temporary file,
 * reads the file back, tokenizes the line(s) and asserts on the values found.
 * All checks use non-fatal EXPECT_* macros so a failure never aborts the suite.
 */

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_copyCounters.hpp"
#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Helper: split one line of the output file into whitespace-delimited tokens.
// -----------------------------------------------------------------------------
std::vector<std::string> wnbp_tokenize(const std::string& line)
{
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok)
        tokens.push_back(tok);
    return tokens;
}

// -----------------------------------------------------------------------------
// Helper: read every non-empty line of a file and return them.
// -----------------------------------------------------------------------------
std::vector<std::string> wnbp_read_lines(const std::string& fileName)
{
    std::vector<std::string> lines;
    std::ifstream in(fileName);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty())
            lines.push_back(line);
    }
    return lines;
}

// -----------------------------------------------------------------------------
// Helper: build a copyCounters object with fully deterministic contents.
//
// nBoundPairs is filled with the pattern value = 100 + index so that any
// mis-indexing by the function under test is immediately visible.
// -----------------------------------------------------------------------------
copyCounters wnbp_make_counters(const std::vector<int>& proPairlist, size_t nBoundPairsSize)
{
    copyCounters counters;
    counters.proPairlist = proPairlist;

    counters.nBoundPairs.clear();
    for (size_t i = 0; i < nBoundPairsSize; ++i)
        counters.nBoundPairs.push_back(static_cast<int>(100 + i));

    // The eight trailing diagnostic counters - all distinct so ordering is testable.
    counters.nLoops = 11;
    counters.nCancelOverlapPartner = 22;
    counters.nCancelOverlapSystem = 33;
    counters.nCancelSpanBox = 44;
    counters.nCancelDisplace2D = 55;
    counters.nCancelDisplace3D = 66;
    counters.nCancelDisplace3Dto2D = 77;
    counters.nAssocSuccess = 88;

    return counters;
}

// -----------------------------------------------------------------------------
// Helper: default Parameters used by most tests (no restart offsets).
// -----------------------------------------------------------------------------
Parameters wnbp_make_params()
{
    Parameters params;
    params.numMolTypes = 2;
    params.timeStep = 1.0; // microseconds
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;
    params.implicitLipid = false;
    return params;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Basic single-line output.
//
// Criteria: exactly one line is written, the token count equals
// 1 (time) + proPairlist.size() + 8 (diagnostic counters).
// -----------------------------------------------------------------------------
void test_wnbp_basic_line_layout()
{
    std::cerr << "\n[TEST] test_wnbp_basic_line_layout\n"
              << "  Source file:   src/io/write_NboundPairs.cpp\n"
              << "  Function:      write_NboundPairs\n"
              << "  Scenario:      one call with a 3-entry proPairlist.\n"
              << "  Pass criteria: one output line containing 1 + 3 + 8 = 12 tokens.\n";

    const std::string fileName = "test_wnbp_basic.dat";

    // proPairlist selects nBoundPairs entries 0, 1 and 2.
    copyCounters counters = wnbp_make_counters({ 0, 1, 2 }, 4);
    Parameters params = wnbp_make_params();
    std::vector<Molecule> moleculeList; // unused by the function, but required

    std::ofstream outfile(fileName);
    ASSERT_TRUE(outfile.is_open()) << "Could not open temporary output file";

    std::cerr << "  Calling write_NboundPairs with it = 0 ...\n";
    write_NboundPairs(counters, outfile, 0, params, moleculeList);
    outfile.close();

    std::vector<std::string> lines = wnbp_read_lines(fileName);
    EXPECT_EQ(lines.size(), 1u) << "write_NboundPairs should emit exactly one line per call";

    if (!lines.empty()) {
        std::vector<std::string> tokens = wnbp_tokenize(lines[0]);
        std::cerr << "  Line written: \"" << lines[0] << "\" (" << tokens.size() << " tokens)\n";
        EXPECT_EQ(tokens.size(), 1u + counters.proPairlist.size() + 8u)
            << "Expected time + one column per proPairlist entry + 8 diagnostic counters";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 2: The bound-pair columns must be looked up THROUGH proPairlist.
//
// Criteria: with proPairlist = {2, 0, 3} the written columns must be
// nBoundPairs[2], nBoundPairs[0], nBoundPairs[3] == 102, 100, 103.
// -----------------------------------------------------------------------------
void test_wnbp_bound_pair_indexing()
{
    std::cerr << "\n[TEST] test_wnbp_bound_pair_indexing\n"
              << "  Source file:   src/io/write_NboundPairs.cpp\n"
              << "  Function:      write_NboundPairs\n"
              << "  Scenario:      proPairlist holds out-of-order indices {2,0,3}.\n"
              << "  Pass criteria: written columns are nBoundPairs[2], [0], [3]\n"
              << "                 i.e. 102, 100, 103 for the seeded pattern.\n";

    const std::string fileName = "test_wnbp_indexing.dat";

    copyCounters counters = wnbp_make_counters({ 2, 0, 3 }, 5);
    Parameters params = wnbp_make_params();
    std::vector<Molecule> moleculeList;

    std::ofstream outfile(fileName);
    ASSERT_TRUE(outfile.is_open()) << "Could not open temporary output file";
    write_NboundPairs(counters, outfile, 0, params, moleculeList);
    outfile.close();

    std::vector<std::string> lines = wnbp_read_lines(fileName);
    ASSERT_EQ(lines.size(), 1u) << "Expected a single output line";

    std::vector<std::string> tokens = wnbp_tokenize(lines[0]);
    std::cerr << "  Line written: \"" << lines[0] << "\"\n";
    ASSERT_GE(tokens.size(), 4u) << "Line is too short to hold the bound-pair columns";

    // Token 0 is the time; tokens 1..3 are the bound-pair columns.
    EXPECT_EQ(tokens[1], "102") << "First bound-pair column should be nBoundPairs[2]";
    EXPECT_EQ(tokens[2], "100") << "Second bound-pair column should be nBoundPairs[0]";
    EXPECT_EQ(tokens[3], "103") << "Third bound-pair column should be nBoundPairs[3]";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 3: The trailing eight diagnostic counters must appear in a fixed order.
//
// Criteria: the last eight tokens equal nLoops, nCancelOverlapPartner,
// nCancelOverlapSystem, nCancelSpanBox, nCancelDisplace2D, nCancelDisplace3D,
// nCancelDisplace3Dto2D, nAssocSuccess.
// -----------------------------------------------------------------------------
void test_wnbp_diagnostic_counter_order()
{
    std::cerr << "\n[TEST] test_wnbp_diagnostic_counter_order\n"
              << "  Source file:   src/io/write_NboundPairs.cpp\n"
              << "  Function:      write_NboundPairs\n"
              << "  Scenario:      distinct values are placed in each diagnostic counter.\n"
              << "  Pass criteria: last 8 columns are 11,22,33,44,55,66,77,88 in order.\n";

    const std::string fileName = "test_wnbp_counters.dat";

    copyCounters counters = wnbp_make_counters({ 1 }, 3);
    Parameters params = wnbp_make_params();
    std::vector<Molecule> moleculeList;

    std::ofstream outfile(fileName);
    ASSERT_TRUE(outfile.is_open()) << "Could not open temporary output file";
    write_NboundPairs(counters, outfile, 0, params, moleculeList);
    outfile.close();

    std::vector<std::string> lines = wnbp_read_lines(fileName);
    ASSERT_EQ(lines.size(), 1u) << "Expected a single output line";

    std::vector<std::string> tokens = wnbp_tokenize(lines[0]);
    std::cerr << "  Line written: \"" << lines[0] << "\"\n";
    ASSERT_EQ(tokens.size(), 1u + 1u + 8u)
        << "Expected time + 1 bound-pair column + 8 diagnostic counters";

    // The expected trailing sequence, in the exact order the function writes them.
    const std::vector<std::string> expected {
        "11", // nLoops
        "22", // nCancelOverlapPartner
        "33", // nCancelOverlapSystem
        "44", // nCancelSpanBox
        "55", // nCancelDisplace2D
        "66", // nCancelDisplace3D
        "77", // nCancelDisplace3Dto2D
        "88" // nAssocSuccess
    };

    const size_t offset = tokens.size() - expected.size();
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(tokens[offset + i], expected[i])
            << "Diagnostic counter column " << i << " is out of order or wrong";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 4: The time column with no restart offsets.
//
// Criteria: time == it * timeStep * 1E-6 (microseconds -> seconds).
// -----------------------------------------------------------------------------
void test_wnbp_time_column_no_restart()
{
    std::cerr << "\n[TEST] test_wnbp_time_column_no_restart\n"
              << "  Source file:   src/io/write_NboundPairs.cpp\n"
              << "  Function:      write_NboundPairs\n"
              << "  Scenario:      itrRestartFrom = 0, timeRestartFrom = 0,\n"
              << "                 timeStep = 2 us, it = 1000.\n"
              << "  Pass criteria: time column == 1000 * 2 * 1e-6 = 0.002 s.\n";

    const std::string fileName = "test_wnbp_time_norestart.dat";

    copyCounters counters = wnbp_make_counters({ 0 }, 2);
    Parameters params = wnbp_make_params();
    params.timeStep = 2.0; // microseconds

    std::vector<Molecule> moleculeList;

    const int it = 1000;
    std::ofstream outfile(fileName);
    ASSERT_TRUE(outfile.is_open()) << "Could not open temporary output file";
    write_NboundPairs(counters, outfile, it, params, moleculeList);
    outfile.close();

    std::vector<std::string> lines = wnbp_read_lines(fileName);
    ASSERT_EQ(lines.size(), 1u) << "Expected a single output line";

    std::vector<std::string> tokens = wnbp_tokenize(lines[0]);
    ASSERT_FALSE(tokens.empty()) << "Output line should not be empty";

    const double writtenTime = std::stod(tokens[0]);
    const double expectedTime = it * params.timeStep * 1E-6;
    std::cerr << "  Written time = " << writtenTime << " s, expected = " << expectedTime << " s\n";
    EXPECT_NEAR(writtenTime, expectedTime, 1e-12)
        << "Time column should be it * timeStep * 1E-6 when no restart offsets are used";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 5: The time column when restarting from a previous simulation.
//
// Criteria: time == (it - itrRestartFrom) * timeStep * 1E-6 + timeRestartFrom.
// -----------------------------------------------------------------------------
void test_wnbp_time_column_with_restart()
{
    std::cerr << "\n[TEST] test_wnbp_time_column_with_restart\n"
              << "  Source file:   src/io/write_NboundPairs.cpp\n"
              << "  Function:      write_NboundPairs\n"
              << "  Scenario:      itrRestartFrom = 500, timeRestartFrom = 0.5 s,\n"
              << "                 timeStep = 1 us, it = 1500.\n"
              << "  Pass criteria: time column == (1500-500)*1e-6 + 0.5 = 0.501 s.\n";

    const std::string fileName = "test_wnbp_time_restart.dat";

    copyCounters counters = wnbp_make_counters({ 0 }, 2);
    Parameters params = wnbp_make_params();
    params.timeStep = 1.0;
    params.itrRestartFrom = 500;
    params.timeRestartFrom = 0.5;

    std::vector<Molecule> moleculeList;

    const int it = 1500;
    std::ofstream outfile(fileName);
    ASSERT_TRUE(outfile.is_open()) << "Could not open temporary output file";
    write_NboundPairs(counters, outfile, it, params, moleculeList);
    outfile.close();

    std::vector<std::string> lines = wnbp_read_lines(fileName);
    ASSERT_EQ(lines.size(), 1u) << "Expected a single output line";

    std::vector<std::string> tokens = wnbp_tokenize(lines[0]);
    ASSERT_FALSE(tokens.empty()) << "Output line should not be empty";

    const double writtenTime = std::stod(tokens[0]);
    const double expectedTime
        = (it - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom;
    std::cerr << "  Written time = " << writtenTime << " s, expected = " << expectedTime << " s\n";
    EXPECT_NEAR(writtenTime, expectedTime, 1e-9)
        << "Time column must include the restart iteration and time offsets";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 6: Empty proPairlist (degenerate but legal input).
//
// Criteria: the line still contains the time column plus the 8 diagnostic
// counters, i.e. exactly 9 tokens, and no crash occurs.
// -----------------------------------------------------------------------------
void test_wnbp_empty_pair_list()
{
    std::cerr << "\n[TEST] test_wnbp_empty_pair_list\n"
              << "  Source file:   src/io/write_NboundPairs.cpp\n"
              << "  Function:      write_NboundPairs\n"
              << "  Scenario:      proPairlist is empty (no reactive pairs).\n"
              << "  Pass criteria: line has exactly 1 + 0 + 8 = 9 tokens.\n";

    const std::string fileName = "test_wnbp_empty.dat";

    copyCounters counters = wnbp_make_counters({}, 3);
    Parameters params = wnbp_make_params();
    std::vector<Molecule> moleculeList;

    std::ofstream outfile(fileName);
    ASSERT_TRUE(outfile.is_open()) << "Could not open temporary output file";
    write_NboundPairs(counters, outfile, 10, params, moleculeList);
    outfile.close();

    std::vector<std::string> lines = wnbp_read_lines(fileName);
    EXPECT_EQ(lines.size(), 1u) << "A line should still be written for an empty pair list";

    if (!lines.empty()) {
        std::vector<std::string> tokens = wnbp_tokenize(lines[0]);
        std::cerr << "  Line written: \"" << lines[0] << "\" (" << tokens.size() << " tokens)\n";
        EXPECT_EQ(tokens.size(), 9u)
            << "Expected only the time column and the 8 diagnostic counters";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 7: Repeated calls append one line each, with increasing time stamps.
//
// Criteria: three calls produce three lines; time stamps are strictly
// increasing and each line has an identical token count.
// -----------------------------------------------------------------------------
void test_wnbp_multiple_calls_append()
{
    std::cerr << "\n[TEST] test_wnbp_multiple_calls_append\n"
              << "  Source file:   src/io/write_NboundPairs.cpp\n"
              << "  Function:      write_NboundPairs\n"
              << "  Scenario:      three successive calls at it = 0, 100, 200.\n"
              << "  Pass criteria: three lines, equal token counts, strictly\n"
              << "                 increasing time stamps.\n";

    const std::string fileName = "test_wnbp_multi.dat";

    copyCounters counters = wnbp_make_counters({ 0, 1 }, 3);
    Parameters params = wnbp_make_params();
    std::vector<Molecule> moleculeList;

    std::ofstream outfile(fileName);
    ASSERT_TRUE(outfile.is_open()) << "Could not open temporary output file";

    const std::vector<int> iterations { 0, 100, 200 };
    for (int it : iterations) {
        std::cerr << "  Calling write_NboundPairs with it = " << it << " ...\n";
        write_NboundPairs(counters, outfile, it, params, moleculeList);
    }
    outfile.close();

    std::vector<std::string> lines = wnbp_read_lines(fileName);
    EXPECT_EQ(lines.size(), iterations.size())
        << "Each call should append exactly one line to the stream";

    double previousTime = -1.0;
    for (size_t i = 0; i < lines.size(); ++i) {
        std::vector<std::string> tokens = wnbp_tokenize(lines[i]);
        EXPECT_EQ(tokens.size(), 1u + counters.proPairlist.size() + 8u)
            << "Line " << i << " has an unexpected number of columns";

        if (!tokens.empty()) {
            const double currentTime = std::stod(tokens[0]);
            std::cerr << "    line " << i << " time = " << currentTime << " s\n";
            EXPECT_GT(currentTime, previousTime)
                << "Time stamps should increase with the iteration counter";
            previousTime = currentTime;
        }
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 8: Zeroed counters are written as zeros (not skipped).
//
// Criteria: with every counter set to zero the line still contains the full
// set of columns, all equal to "0" (apart from the time column).
// -----------------------------------------------------------------------------
void test_wnbp_all_zero_counters()
{
    std::cerr << "\n[TEST] test_wnbp_all_zero_counters\n"
              << "  Source file:   src/io/write_NboundPairs.cpp\n"
              << "  Function:      write_NboundPairs\n"
              << "  Scenario:      every counter (bound pairs and diagnostics) is 0.\n"
              << "  Pass criteria: all non-time columns are printed as \"0\".\n";

    const std::string fileName = "test_wnbp_zeros.dat";

    // Build counters manually so everything is zero.
    copyCounters counters;
    counters.proPairlist = { 0, 1 };
    counters.nBoundPairs = { 0, 0 };
    counters.nLoops = 0;
    counters.nCancelOverlapPartner = 0;
    counters.nCancelOverlapSystem = 0;
    counters.nCancelSpanBox = 0;
    counters.nCancelDisplace2D = 0;
    counters.nCancelDisplace3D = 0;
    counters.nCancelDisplace3Dto2D = 0;
    counters.nAssocSuccess = 0;

    Parameters params = wnbp_make_params();
    std::vector<Molecule> moleculeList;

    std::ofstream outfile(fileName);
    ASSERT_TRUE(outfile.is_open()) << "Could not open temporary output file";
    write_NboundPairs(counters, outfile, 0, params, moleculeList);
    outfile.close();

    std::vector<std::string> lines = wnbp_read_lines(fileName);
    ASSERT_EQ(lines.size(), 1u) << "Expected a single output line";

    std::vector<std::string> tokens = wnbp_tokenize(lines[0]);
    std::cerr << "  Line written: \"" << lines[0] << "\"\n";
    EXPECT_EQ(tokens.size(), 1u + counters.proPairlist.size() + 8u)
        << "Zero-valued columns must still be written";

    // Every column after the time stamp should literally be "0".
    for (size_t i = 1; i < tokens.size(); ++i) {
        EXPECT_EQ(tokens[i], "0") << "Column " << i << " should be printed as 0";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named helper is invoked from its own TEST so that
// individual failures are reported separately while all tests still execute.
// -----------------------------------------------------------------------------
TEST(WriteNboundPairs, BasicLineLayout) { test_wnbp_basic_line_layout(); }
TEST(WriteNboundPairs, BoundPairIndexing) { test_wnbp_bound_pair_indexing(); }
TEST(WriteNboundPairs, DiagnosticCounterOrder) { test_wnbp_diagnostic_counter_order(); }
TEST(WriteNboundPairs, TimeColumnNoRestart) { test_wnbp_time_column_no_restart(); }
TEST(WriteNboundPairs, TimeColumnWithRestart) { test_wnbp_time_column_with_restart(); }
TEST(WriteNboundPairs, EmptyPairList) { test_wnbp_empty_pair_list(); }
TEST(WriteNboundPairs, MultipleCallsAppend) { test_wnbp_multiple_calls_append(); }
TEST(WriteNboundPairs, AllZeroCounters) { test_wnbp_all_zero_counters(); }