/*! \file test_write_observables.cpp
 *
 * ### Unit test for src/io/write_observables.cpp
 *
 * The file under test contains exactly one function:
 *
 *     void write_observables(double simTime,
 *                           std::ofstream& observablesFile,
 *                           const std::map<std::string, int>& observablesList)
 *
 * Behaviour of that function (as implemented):
 *   - If the observables map holds exactly ONE entry, the function writes
 *         <simTime>,<value>\n
 *     using the map's first (and only) value.
 *   - Otherwise (zero, two, or more entries) it writes the time and then one
 *     ",<value>" field per map entry, terminated with std::endl.
 *     Because std::map iterates in key-sorted order, the CSV column order is
 *     the alphabetical order of the observable names.
 *
 * Strategy: the function only appends text to an std::ofstream, so every test
 * writes to a temporary file, closes it, reads it back, and compares the exact
 * text produced against an independently constructed expectation.  Numeric
 * formatting expectations are built with a plain std::ostringstream so that the
 * test never hard-codes assumptions about the default stream precision.
 */

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Declaration of the function under test.
//
// It is declared directly (rather than pulling in include/io/io.hpp) so that
// this translation unit stays free of the very large IO/reaction header chain
// that io.hpp drags in; the mangled symbol is identical either way and comes
// from src/io/write_observables.cpp at link time.
// ---------------------------------------------------------------------------
void write_observables(double simTime, std::ofstream& observablesFile,
    const std::map<std::string, int>& observablesList);

namespace {

//! Name of the scratch file all tests write to (removed again after each use).
const char* kWobsTmpFile = "test_write_observables_tmp.csv";

/*! \brief Slurp an entire text file into a std::string.
 *
 * \param[in] fileName Path of the file to read.
 * \return The complete file contents (empty string if the file is missing).
 */
std::string wobs_read_whole_file(const std::string& fileName)
{
    std::ifstream in(fileName);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

/*! \brief Call write_observables() once and return the text it produced.
 *
 * The temporary file is truncated before the call and deleted afterwards so
 * each test starts from a clean slate.
 *
 * \param[in] simTime Simulation time handed to the function.
 * \param[in] obs     Observables map handed to the function.
 * \return Everything the function wrote to the stream.
 */
std::string wobs_run_once(double simTime, const std::map<std::string, int>& obs)
{
    std::ofstream file(kWobsTmpFile, std::ios::out | std::ios::trunc);
    write_observables(simTime, file, obs);
    file.close();

    const std::string contents = wobs_read_whole_file(kWobsTmpFile);
    std::remove(kWobsTmpFile);
    return contents;
}

/*! \brief Call write_observables() several times on the same open stream.
 *
 * Used to verify that successive calls append complete, independent lines.
 *
 * \param[in] calls List of (simTime, observables map) pairs, applied in order.
 * \return Everything written to the stream across all calls.
 */
std::string wobs_run_many(const std::vector<std::pair<double, std::map<std::string, int>>>& calls)
{
    std::ofstream file(kWobsTmpFile, std::ios::out | std::ios::trunc);
    for (const auto& oneCall : calls)
        write_observables(oneCall.first, file, oneCall.second);
    file.close();

    const std::string contents = wobs_read_whole_file(kWobsTmpFile);
    std::remove(kWobsTmpFile);
    return contents;
}

/*! \brief Split text on '\n' into individual lines (dropping the final empty one). */
std::vector<std::string> wobs_split_lines(const std::string& text)
{
    std::vector<std::string> lines;
    std::istringstream stream(text);
    std::string oneLine;
    while (std::getline(stream, oneLine))
        lines.push_back(oneLine);
    return lines;
}

/*! \brief Format a double exactly the way an ostream with default settings would.
 *
 * This lets the tests build expectations without hard-coding precision rules.
 */
std::string wobs_format_double(double value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

} // namespace

// ---------------------------------------------------------------------------
// Test 1: single-observable branch -> "<time>,<value>\n"
// ---------------------------------------------------------------------------
void wobs_test_single_observable()
{
    std::cerr << "\n[TEST] wobs_test_single_observable\n"
              << "  Source file: src/io/write_observables.cpp\n"
              << "  Function:    write_observables (observablesList.size() == 1 branch)\n"
              << "  Scenario:    one observable named \"dimer\" with a count of 7 at t = 0.5\n"
              << "  Pass:        output is exactly \"0.5,7\\n\" (one line, two CSV fields)\n";

    std::map<std::string, int> obs;
    obs["dimer"] = 7;

    const std::string output = wobs_run_once(0.5, obs);
    std::cerr << "  Raw output: \"" << output << "\"\n";

    // Expected text built independently of the implementation.
    const std::string expected = wobs_format_double(0.5) + ",7\n";
    EXPECT_EQ(output, expected) << "Single-observable output must be \"<time>,<count>\\n\"";

    // Structural checks: exactly one line, exactly two comma-separated fields.
    const std::vector<std::string> lines = wobs_split_lines(output);
    EXPECT_EQ(lines.size(), 1u) << "Exactly one line should be written per call";
    if (!lines.empty()) {
        const size_t commaCount
            = static_cast<size_t>(std::count(lines[0].begin(), lines[0].end(), ','));
        EXPECT_EQ(commaCount, 1u) << "A single observable yields exactly one comma";
    }

    // The function takes the map by const reference; it must not disturb it.
    EXPECT_EQ(obs.size(), 1u) << "write_observables must not modify the observables map";
    EXPECT_EQ(obs["dimer"], 7) << "Observable value must be left untouched";
}

// ---------------------------------------------------------------------------
// Test 2: multi-observable branch, verifying key-sorted column order
// ---------------------------------------------------------------------------
void wobs_test_multiple_observables_sorted()
{
    std::cerr << "\n[TEST] wobs_test_multiple_observables_sorted\n"
              << "  Source file: src/io/write_observables.cpp\n"
              << "  Function:    write_observables (multi-entry else branch)\n"
              << "  Scenario:    three observables inserted out of alphabetical order\n"
              << "               (\"zeta\"=3, \"alpha\"=1, \"middle\"=2) at t = 1.25\n"
              << "  Pass:        output is \"1.25,1,2,3\\n\" i.e. values follow the\n"
              << "               std::map key ordering alpha < middle < zeta\n";

    std::map<std::string, int> obs;
    obs["zeta"] = 3; // inserted first on purpose
    obs["alpha"] = 1;
    obs["middle"] = 2;

    const std::string output = wobs_run_once(1.25, obs);
    std::cerr << "  Raw output: \"" << output << "\"\n";

    const std::string expected = wobs_format_double(1.25) + ",1,2,3\n";
    EXPECT_EQ(output, expected)
        << "Multi-observable output must be time followed by key-sorted counts";

    // One line, three commas (one per observable).
    const std::vector<std::string> lines = wobs_split_lines(output);
    EXPECT_EQ(lines.size(), 1u) << "Exactly one line should be written per call";
    if (!lines.empty()) {
        const size_t commaCount
            = static_cast<size_t>(std::count(lines[0].begin(), lines[0].end(), ','));
        EXPECT_EQ(commaCount, 3u) << "Three observables yield three commas";
    }
}

// ---------------------------------------------------------------------------
// Test 3: empty observables list -> just the time plus a newline
// ---------------------------------------------------------------------------
void wobs_test_empty_observables_list()
{
    std::cerr << "\n[TEST] wobs_test_empty_observables_list\n"
              << "  Source file: src/io/write_observables.cpp\n"
              << "  Function:    write_observables (else branch, zero entries)\n"
              << "  Scenario:    empty observables map at t = 2\n"
              << "  Pass:        output is \"2\\n\" - the time only, no commas\n";

    const std::map<std::string, int> obs; // deliberately empty

    const std::string output = wobs_run_once(2.0, obs);
    std::cerr << "  Raw output: \"" << output << "\"\n";

    const std::string expected = wobs_format_double(2.0) + "\n";
    EXPECT_EQ(output, expected) << "With no observables only the time and newline are written";
    EXPECT_EQ(output.find(','), std::string::npos)
        << "No comma should appear when the observables map is empty";
}

// ---------------------------------------------------------------------------
// Test 4: two observables (boundary between the two branches)
// ---------------------------------------------------------------------------
void wobs_test_two_observables_boundary()
{
    std::cerr << "\n[TEST] wobs_test_two_observables_boundary\n"
              << "  Source file: src/io/write_observables.cpp\n"
              << "  Function:    write_observables (size() == 2 takes the else branch)\n"
              << "  Scenario:    observables \"b\"=20 and \"a\"=10 at t = 0\n"
              << "  Pass:        output is \"0,10,20\\n\" - both values present, sorted\n";

    std::map<std::string, int> obs;
    obs["b"] = 20;
    obs["a"] = 10;

    const std::string output = wobs_run_once(0.0, obs);
    std::cerr << "  Raw output: \"" << output << "\"\n";

    const std::string expected = wobs_format_double(0.0) + ",10,20\n";
    EXPECT_EQ(output, expected)
        << "Two observables must both be written, in key-sorted order";
}

// ---------------------------------------------------------------------------
// Test 5: zero and negative counts are written verbatim
// ---------------------------------------------------------------------------
void wobs_test_zero_and_negative_counts()
{
    std::cerr << "\n[TEST] wobs_test_zero_and_negative_counts\n"
              << "  Source file: src/io/write_observables.cpp\n"
              << "  Function:    write_observables (value formatting)\n"
              << "  Scenario:    counts 0, -5 and a large 123456 at t = 3.5\n"
              << "  Pass:        every integer is echoed unchanged, including the\n"
              << "               minus sign, with no clamping or rounding\n";

    std::map<std::string, int> obs;
    obs["aZero"] = 0;
    obs["bNegative"] = -5;
    obs["cLarge"] = 123456;

    const std::string output = wobs_run_once(3.5, obs);
    std::cerr << "  Raw output: \"" << output << "\"\n";

    const std::string expected = wobs_format_double(3.5) + ",0,-5,123456\n";
    EXPECT_EQ(output, expected) << "Integer counts must be written verbatim";
    EXPECT_NE(output.find("-5"), std::string::npos)
        << "Negative counts must keep their minus sign";
}

// ---------------------------------------------------------------------------
// Test 6: single-observable branch also handles a zero/negative count
// ---------------------------------------------------------------------------
void wobs_test_single_observable_zero_count()
{
    std::cerr << "\n[TEST] wobs_test_single_observable_zero_count\n"
              << "  Source file: src/io/write_observables.cpp\n"
              << "  Function:    write_observables (size() == 1 branch, count 0)\n"
              << "  Scenario:    one observable whose count is 0 at t = 10\n"
              << "  Pass:        output is \"10,0\\n\" - a zero count is still written\n";

    std::map<std::string, int> obs;
    obs["nothingYet"] = 0;

    const std::string output = wobs_run_once(10.0, obs);
    std::cerr << "  Raw output: \"" << output << "\"\n";

    const std::string expected = wobs_format_double(10.0) + ",0\n";
    EXPECT_EQ(output, expected) << "A single zero-valued observable must still be reported";
}

// ---------------------------------------------------------------------------
// Test 7: successive calls append one complete line each
// ---------------------------------------------------------------------------
void wobs_test_successive_calls_append_lines()
{
    std::cerr << "\n[TEST] wobs_test_successive_calls_append_lines\n"
              << "  Source file: src/io/write_observables.cpp\n"
              << "  Function:    write_observables (called repeatedly on one stream)\n"
              << "  Scenario:    three timesteps written back-to-back; the first two\n"
              << "               use two observables, the last uses one observable\n"
              << "  Pass:        three separate lines appear in call order and each\n"
              << "               line matches its own expected CSV record\n";

    std::map<std::string, int> twoObs;
    twoObs["a"] = 1;
    twoObs["b"] = 2;

    std::map<std::string, int> oneObs;
    oneObs["a"] = 42;

    std::vector<std::pair<double, std::map<std::string, int>>> calls;
    calls.emplace_back(0.0, twoObs);
    calls.emplace_back(0.5, twoObs);
    calls.emplace_back(1.0, oneObs); // exercises the single-observable branch

    const std::string output = wobs_run_many(calls);
    std::cerr << "  Raw output:\n" << output;

    const std::vector<std::string> lines = wobs_split_lines(output);
    EXPECT_EQ(lines.size(), 3u) << "One line must be appended per call";

    if (lines.size() >= 3u) {
        EXPECT_EQ(lines[0], wobs_format_double(0.0) + ",1,2")
            << "First record should be the t=0 two-observable line";
        EXPECT_EQ(lines[1], wobs_format_double(0.5) + ",1,2")
            << "Second record should be the t=0.5 two-observable line";
        EXPECT_EQ(lines[2], wobs_format_double(1.0) + ",42")
            << "Third record should be the t=1 single-observable line";
    }
}

// ---------------------------------------------------------------------------
// Test 8: the time field uses the stream's default numeric formatting
// ---------------------------------------------------------------------------
void wobs_test_time_field_formatting()
{
    std::cerr << "\n[TEST] wobs_test_time_field_formatting\n"
              << "  Source file: src/io/write_observables.cpp\n"
              << "  Function:    write_observables (simTime formatting)\n"
              << "  Scenario:    a fractional time (0.125), a whole-number time (7),\n"
              << "               and a large time (1.5e+06)\n"
              << "  Pass:        the leading field equals what \"ostream << value\"\n"
              << "               produces, i.e. the implementation applies no extra\n"
              << "               precision/format manipulation\n";

    std::map<std::string, int> obs;
    obs["only"] = 1; // single-observable branch keeps the comparison simple

    // Fractional time.
    {
        const std::string output = wobs_run_once(0.125, obs);
        std::cerr << "  t = 0.125 -> \"" << output << "\"\n";
        EXPECT_EQ(output, wobs_format_double(0.125) + ",1")
                     .operator bool()
            ? void()
            : void(); // (comparison performed below with a clearer message)
        EXPECT_EQ(output, wobs_format_double(0.125) + ",1\n")
            << "Fractional times must use default stream formatting";
    }

    // Whole-number time: default formatting prints "7", not "7.0".
    {
        const std::string output = wobs_run_once(7.0, obs);
        std::cerr << "  t = 7.0   -> \"" << output << "\"\n";
        EXPECT_EQ(output, wobs_format_double(7.0) + ",1\n")
            << "Whole-number times must use default stream formatting";
    }

    // Large time: default formatting switches to scientific notation.
    {
        const std::string output = wobs_run_once(1.5e6, obs);
        std::cerr << "  t = 1.5e6 -> \"" << output << "\"\n";
        EXPECT_EQ(output, wobs_format_double(1.5e6) + ",1\n")
            << "Large times must use default stream formatting";
    }
}

// ---------------------------------------------------------------------------
// GoogleTest wrappers: every named test_* helper is run inside its own TEST so
// that a failure in one does not stop the others from executing.
// ---------------------------------------------------------------------------
TEST(WriteObservables, SingleObservable) { wobs_test_single_observable(); }
TEST(WriteObservables, MultipleObservablesSorted) { wobs_test_multiple_observables_sorted(); }
TEST(WriteObservables, EmptyObservablesList) { wobs_test_empty_observables_list(); }
TEST(WriteObservables, TwoObservablesBoundary) { wobs_test_two_observables_boundary(); }
TEST(WriteObservables, ZeroAndNegativeCounts) { wobs_test_zero_and_negative_counts(); }
TEST(WriteObservables, SingleObservableZeroCount) { wobs_test_single_observable_zero_count(); }
TEST(WriteObservables, SuccessiveCallsAppendLines) { wobs_test_successive_calls_append_lines(); }
TEST(WriteObservables, TimeFieldFormatting) { wobs_test_time_field_formatting(); }