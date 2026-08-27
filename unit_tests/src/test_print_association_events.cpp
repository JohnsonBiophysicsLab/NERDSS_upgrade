/*! \file test_print_association_events.cpp
 *
 * ### Unit test for ../src/reactions/print_association_events.cpp
 *
 * This test exercises the single function defined in that translation unit:
 *
 *     void print_association_events(copyCounters& counterArrays,
 *                                   std::ofstream& outfile,
 *                                   int it,
 *                                   Parameters params)
 *
 * The routine writes a human readable summary of the association events that
 * have been binned into three histograms held on `copyCounters`:
 *
 *     counterArrays.events3D      -> printed first,  prefix "3D"
 *     counterArrays.events3Dto2D  -> printed second, prefix "3Dto2D"
 *     counterArrays.events2D      -> printed last,   prefix "2D"
 *
 * The binning/labelling rules taken directly from the implementation are:
 *
 *   * a header line is always written:
 *         "time (s): " << (it - params.itrRestartFrom) * params.timeStep * 1E-6
 *                          + params.timeRestartFrom
 *   * a bin is only written when its counter is strictly greater than zero
 *   * bin 0                        -> "<prefix> dimers"
 *   * bins 1 .. maxSingles-1 (1..9)-> "<prefix> n=<i>"
 *   * bins maxSingles .. size-2    -> "<prefix> n=<n> to <n+9>" with
 *                                     n = (i - maxSingles + 1) * 10
 *   * the final bin (size-1)       -> "<prefix> n><n>"  (open ended bucket)
 *   * the number of bins scanned is `counterArrays.eventArraySize`, NOT the
 *     size of the underlying vectors.
 *
 * The tests below drive the function with a real std::ofstream backed by a
 * temporary file, read the file back, and assert on the exact text produced.
 *
 * NOTE: the function indexes events3D/events2D/events3Dto2D from 0 up to
 * eventArraySize-1 without any bounds check, so every helper below sizes those
 * vectors to at least eventArraySize before calling into the code under test.
 */

#include "classes/class_Parameters.hpp"
#include "classes/class_copyCounters.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>

namespace {

//! Temporary file used as the destination std::ofstream for every test.
const char* kPaeTmpFile = "test_print_association_events_output.txt";

/*! \brief Build a fully initialised copyCounters object.
 *
 * `print_association_events` reads exactly `eventArraySize` elements out of the
 * three event vectors, so all three are sized accordingly and zero filled.
 *
 * \param[in] arraySize number of histogram bins to allocate / advertise.
 * \return a copyCounters ready to hand to the function under test.
 */
copyCounters pae_make_counters(int arraySize)
{
    copyCounters counters;
    counters.eventArraySize = arraySize;
    counters.events3D.assign(arraySize, 0);
    counters.events2D.assign(arraySize, 0);
    counters.events3Dto2D.assign(arraySize, 0);
    return counters;
}

/*! \brief Build a Parameters object with only the fields the function reads. */
Parameters pae_make_params(double timeStep, long long itrRestartFrom, double timeRestartFrom)
{
    Parameters params;
    params.timeStep = timeStep;
    params.itrRestartFrom = itrRestartFrom;
    params.timeRestartFrom = timeRestartFrom;
    return params;
}

/*! \brief Run the function under test and return everything it wrote.
 *
 * Opens the temporary file truncated, calls print_association_events, closes
 * the stream (so all buffered output is flushed), slurps the file back into a
 * string and finally deletes the file.
 */
std::string pae_run(copyCounters& counters, int it, const Parameters& params)
{
    {
        std::ofstream outFile(kPaeTmpFile, std::ios::out | std::ios::trunc);
        // Fail loudly (but non-fatally) if the stream could not be opened.
        EXPECT_TRUE(outFile.is_open()) << "Could not open temporary output file " << kPaeTmpFile;
        print_association_events(counters, outFile, it, params);
    } // ofstream destructor closes and flushes here

    std::ifstream inFile(kPaeTmpFile);
    std::stringstream buffer;
    buffer << inFile.rdbuf();
    inFile.close();
    std::remove(kPaeTmpFile);
    return buffer.str();
}

/*! \brief Split a blob of text into individual lines (dropping the newlines). */
std::vector<std::string> pae_split_lines(const std::string& text)
{
    std::vector<std::string> lines;
    std::istringstream iss(text);
    std::string line;
    while (std::getline(iss, line))
        lines.push_back(line);
    return lines;
}

/*! \brief Reproduce the exact header line the implementation emits.
 *
 * The value goes through operator<< on an ostream with default formatting
 * (6 significant digits), so we build the expected string the same way rather
 * than guessing at a printf format.
 */
std::string pae_expected_time_line(int it, const Parameters& params)
{
    std::ostringstream oss;
    oss << "time (s): "
        << (it - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom;
    return oss.str();
}

/*! \brief Index of a line inside the output, or -1 if it is absent. */
int pae_line_index(const std::vector<std::string>& lines, const std::string& target)
{
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i] == target)
            return static_cast<int>(i);
    }
    return -1;
}

/*! \brief Convenience predicate built on top of pae_line_index(). */
bool pae_has_line(const std::vector<std::string>& lines, const std::string& target)
{
    return pae_line_index(lines, target) >= 0;
}

/*! \brief Dump the captured output to stderr so failures are easy to diagnose. */
void pae_echo(const std::vector<std::string>& lines)
{
    std::cerr << "    captured output (" << lines.size() << " line(s)):\n";
    for (const auto& line : lines)
        std::cerr << "      | " << line << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: with every counter at zero only the time header line is produced.
// -----------------------------------------------------------------------------
void test_pae_header_only_when_no_events()
{
    std::cerr << "\n[TEST] test_pae_header_only_when_no_events\n"
              << "  Source file:   src/reactions/print_association_events.cpp\n"
              << "  Function:      print_association_events\n"
              << "  Scenario:      all three histograms are entirely zero.\n"
              << "  Pass criteria: exactly one line is written, the time header.\n";

    copyCounters counters = pae_make_counters(20);
    Parameters params = pae_make_params(/*timeStep*/ 1.0, /*itrRestartFrom*/ 0, /*timeRestartFrom*/ 0.0);
    const int it = 0;

    const std::vector<std::string> lines = pae_split_lines(pae_run(counters, it, params));
    pae_echo(lines);

    // Nothing but the header should have been written.
    EXPECT_EQ(lines.size(), 1u) << "Only the time header should be written when no events exist";
    if (!lines.empty()) {
        EXPECT_EQ(lines[0], pae_expected_time_line(it, params))
            << "Header line text must match the documented format";
    }
}

// -----------------------------------------------------------------------------
// Test 2: the header time value obeys
//         (it - itrRestartFrom) * timeStep * 1E-6 + timeRestartFrom
// -----------------------------------------------------------------------------
void test_pae_time_header_value()
{
    std::cerr << "\n[TEST] test_pae_time_header_value\n"
              << "  Source file:   src/reactions/print_association_events.cpp\n"
              << "  Function:      print_association_events\n"
              << "  Scenario:      restart offsets are non-zero, so the printed\n"
              << "                 time must be (it - itrRestartFrom)*timeStep*1e-6\n"
              << "                 + timeRestartFrom (timeStep is in microseconds).\n"
              << "  Pass criteria: the numeric value parsed back from the header\n"
              << "                 matches the hand computed value.\n";

    copyCounters counters = pae_make_counters(20);
    // 0.5 us timestep, simulation restarted at iteration 200 / t = 1.5e-3 s.
    Parameters params = pae_make_params(0.5, 200, 1.5e-3);
    const int it = 1000;

    // Hand computed: (1000-200) * 0.5 us = 400 us = 4.0e-4 s, plus 1.5e-3 s.
    const double expectedSeconds = 4.0e-4 + 1.5e-3;

    const std::vector<std::string> lines = pae_split_lines(pae_run(counters, it, params));
    pae_echo(lines);

    ASSERT_FALSE(lines.empty()) << "No output was produced at all";

    // The whole line should match what an ostream would have produced.
    EXPECT_EQ(lines[0], pae_expected_time_line(it, params))
        << "Header line should match ostream default formatting of the time";

    // Parse the numeric tail of the header and compare against the analytic value.
    const std::string prefix = "time (s): ";
    ASSERT_EQ(lines[0].compare(0, prefix.size(), prefix), 0)
        << "Header line should begin with \"time (s): \"";
    const double parsed = std::stod(lines[0].substr(prefix.size()));
    std::cerr << "    parsed time = " << parsed << " s, expected " << expectedSeconds << " s\n";
    EXPECT_NEAR(parsed, expectedSeconds, 1e-12)
        << "Header time must equal (it - itrRestartFrom)*timeStep*1e-6 + timeRestartFrom";
}

// -----------------------------------------------------------------------------
// Test 3: labelling rules for the 3D histogram across every bin category.
// -----------------------------------------------------------------------------
void test_pae_3d_bin_labels()
{
    std::cerr << "\n[TEST] test_pae_3d_bin_labels\n"
              << "  Source file:   src/reactions/print_association_events.cpp\n"
              << "  Function:      print_association_events (events3D loop)\n"
              << "  Scenario:      populate bin 0 (dimers), bins 1 and 9 (singles),\n"
              << "                 bins 10 and 18 (decade ranges) and bin 19 (the\n"
              << "                 open-ended final bucket) of a 20 bin histogram.\n"
              << "  Pass criteria: each bin prints with its documented label and\n"
              << "                 its stored count.\n";

    const int arraySize = 20;
    copyCounters counters = pae_make_counters(arraySize);
    counters.events3D[0] = 5;   // dimerisation bucket
    counters.events3D[1] = 7;   // single n
    counters.events3D[9] = 2;   // last single n (maxSingles - 1)
    counters.events3D[10] = 3;  // first decade bucket: n = 10 -> "10 to 19"
    counters.events3D[18] = 6;  // a = 9  -> n = 90 -> "90 to 99"
    counters.events3D[19] = 4;  // final bucket: a = 10 -> n = 100 -> "n>100"

    Parameters params = pae_make_params(1.0, 0, 0.0);
    const std::vector<std::string> lines = pae_split_lines(pae_run(counters, 0, params));
    pae_echo(lines);

    // Header + six populated bins.
    EXPECT_EQ(lines.size(), 7u) << "Header plus one line per populated 3D bin expected";

    EXPECT_TRUE(pae_has_line(lines, "3D dimers: 5")) << "Bin 0 must be labelled \"3D dimers\"";
    EXPECT_TRUE(pae_has_line(lines, "3D n=1: 7")) << "Bin 1 must be labelled \"3D n=1\"";
    EXPECT_TRUE(pae_has_line(lines, "3D n=9: 2")) << "Bin 9 must be labelled \"3D n=9\"";
    EXPECT_TRUE(pae_has_line(lines, "3D n=10 to 19: 3"))
        << "Bin 10 must be the first decade bucket, 10 through 19";
    EXPECT_TRUE(pae_has_line(lines, "3D n=90 to 99: 6"))
        << "Bin 18 must be the decade bucket 90 through 99";
    EXPECT_TRUE(pae_has_line(lines, "3D n>100: 4"))
        << "The final bin must be the open-ended \"n>100\" bucket";
}

// -----------------------------------------------------------------------------
// Test 4: bins whose counter is not strictly positive are skipped entirely.
// -----------------------------------------------------------------------------
void test_pae_skips_nonpositive_bins()
{
    std::cerr << "\n[TEST] test_pae_skips_nonpositive_bins\n"
              << "  Source file:   src/reactions/print_association_events.cpp\n"
              << "  Function:      print_association_events\n"
              << "  Scenario:      one bin holds a positive count, one holds zero\n"
              << "                 and one holds a negative count.\n"
              << "  Pass criteria: only the strictly positive bin is written; the\n"
              << "                 guard in the code is `> 0`, not `!= 0`.\n";

    copyCounters counters = pae_make_counters(20);
    counters.events3D[2] = 11; // printed
    counters.events3D[3] = 0;  // skipped (zero)
    counters.events3D[4] = -8; // skipped (negative, because the test is > 0)

    Parameters params = pae_make_params(1.0, 0, 0.0);
    const std::vector<std::string> lines = pae_split_lines(pae_run(counters, 0, params));
    pae_echo(lines);

    EXPECT_EQ(lines.size(), 2u) << "Only the header and the single positive bin should appear";
    EXPECT_TRUE(pae_has_line(lines, "3D n=2: 11")) << "The positive bin must be printed";
    EXPECT_FALSE(pae_has_line(lines, "3D n=3: 0")) << "A zero counter must not be printed";
    EXPECT_FALSE(pae_has_line(lines, "3D n=4: -8")) << "A negative counter must not be printed";
}

// -----------------------------------------------------------------------------
// Test 5: the three histograms are emitted in the order 3D, 3Dto2D, 2D and each
//         uses its own text prefix.
// -----------------------------------------------------------------------------
void test_pae_channel_prefixes_and_ordering()
{
    std::cerr << "\n[TEST] test_pae_channel_prefixes_and_ordering\n"
              << "  Source file:   src/reactions/print_association_events.cpp\n"
              << "  Function:      print_association_events\n"
              << "  Scenario:      populate the same bins in all three histograms.\n"
              << "  Pass criteria: labels use the \"3D\", \"3Dto2D\" and \"2D\"\n"
              << "                 prefixes and appear in exactly that order.\n";

    copyCounters counters = pae_make_counters(20);
    // Bin 0 (dimers) and bin 5 (single n) populated in every channel.
    counters.events3D[0] = 1;
    counters.events3D[5] = 2;
    counters.events3Dto2D[0] = 3;
    counters.events3Dto2D[5] = 4;
    counters.events2D[0] = 5;
    counters.events2D[5] = 6;

    Parameters params = pae_make_params(1.0, 0, 0.0);
    const std::vector<std::string> lines = pae_split_lines(pae_run(counters, 0, params));
    pae_echo(lines);

    // Header plus two lines per channel.
    EXPECT_EQ(lines.size(), 7u) << "Header plus two populated bins in each of three channels";

    // Every expected label, with the right prefix, must be present.
    const int i3Ddimer = pae_line_index(lines, "3D dimers: 1");
    const int i3Dn5 = pae_line_index(lines, "3D n=5: 2");
    const int iMixDimer = pae_line_index(lines, "3Dto2D dimers: 3");
    const int iMixN5 = pae_line_index(lines, "3Dto2D n=5: 4");
    const int i2Ddimer = pae_line_index(lines, "2D dimers: 5");
    const int i2Dn5 = pae_line_index(lines, "2D n=5: 6");

    EXPECT_GE(i3Ddimer, 0) << "\"3D dimers\" line missing";
    EXPECT_GE(i3Dn5, 0) << "\"3D n=5\" line missing";
    EXPECT_GE(iMixDimer, 0) << "\"3Dto2D dimers\" line missing";
    EXPECT_GE(iMixN5, 0) << "\"3Dto2D n=5\" line missing";
    EXPECT_GE(i2Ddimer, 0) << "\"2D dimers\" line missing";
    EXPECT_GE(i2Dn5, 0) << "\"2D n=5\" line missing";

    // Ordering: the whole 3D block precedes the 3Dto2D block, which precedes 2D.
    if (i3Dn5 >= 0 && iMixDimer >= 0)
        EXPECT_LT(i3Dn5, iMixDimer) << "All 3D lines must be written before any 3Dto2D line";
    if (iMixN5 >= 0 && i2Ddimer >= 0)
        EXPECT_LT(iMixN5, i2Ddimer) << "All 3Dto2D lines must be written before any 2D line";
    // And within a channel, bins come out in increasing index order.
    if (i3Ddimer >= 0 && i3Dn5 >= 0)
        EXPECT_LT(i3Ddimer, i3Dn5) << "Bins within a channel must be written in index order";
}

// -----------------------------------------------------------------------------
// Test 6: eventArraySize (not the vector size) controls both how many bins are
//         scanned and where the open-ended final bucket lands.
// -----------------------------------------------------------------------------
void test_pae_respects_eventArraySize()
{
    std::cerr << "\n[TEST] test_pae_respects_eventArraySize\n"
              << "  Source file:   src/reactions/print_association_events.cpp\n"
              << "  Function:      print_association_events\n"
              << "  Scenario:      the event vectors are longer (30) than the\n"
              << "                 advertised eventArraySize (12).\n"
              << "  Pass criteria: bins at or beyond index 12 are ignored, and the\n"
              << "                 final scanned bin (index 11) becomes the\n"
              << "                 open-ended bucket \"n>20\" because\n"
              << "                 a = 11 - 10 + 1 = 2 and n = 2*10 = 20.\n";

    // Deliberately allocate more storage than eventArraySize advertises.
    copyCounters counters = pae_make_counters(30);
    counters.eventArraySize = 12;

    counters.events3D[10] = 21; // decade bucket: a = 1 -> n = 10 -> "10 to 19"
    counters.events3D[11] = 22; // final scanned bin -> "n>20"
    counters.events3D[12] = 99; // beyond eventArraySize: must be ignored
    counters.events3D[25] = 98; // far beyond: must be ignored

    Parameters params = pae_make_params(1.0, 0, 0.0);
    const std::vector<std::string> lines = pae_split_lines(pae_run(counters, 0, params));
    pae_echo(lines);

    EXPECT_EQ(lines.size(), 3u) << "Header plus only the two in-range populated bins";
    EXPECT_TRUE(pae_has_line(lines, "3D n=10 to 19: 21"))
        << "Bin 10 should still be the 10-to-19 decade bucket";
    EXPECT_TRUE(pae_has_line(lines, "3D n>20: 22"))
        << "Bin eventArraySize-1 becomes the open-ended bucket, here n>20";

    // Nothing from beyond the advertised array size should leak into the output.
    const std::string all = std::accumulate(
        lines.begin(), lines.end(), std::string {},
        [](const std::string& acc, const std::string& s) { return acc + s + "\n"; });
    EXPECT_EQ(all.find("99"), std::string::npos)
        << "Counter stored past eventArraySize must not be printed";
    EXPECT_EQ(all.find("98"), std::string::npos)
        << "Counter stored far past eventArraySize must not be printed";
}

// -----------------------------------------------------------------------------
// Test 7: repeated calls append to the same stream, each preceded by a fresh
//         header line reflecting the new iteration number.
// -----------------------------------------------------------------------------
void test_pae_repeated_calls_append()
{
    std::cerr << "\n[TEST] test_pae_repeated_calls_append\n"
              << "  Source file:   src/reactions/print_association_events.cpp\n"
              << "  Function:      print_association_events\n"
              << "  Scenario:      the function is called twice on one open stream\n"
              << "                 with different iteration numbers.\n"
              << "  Pass criteria: two header lines are present with the correct\n"
              << "                 (different) times, and each is followed by the\n"
              << "                 event lines valid at that call.\n";

    copyCounters counters = pae_make_counters(20);
    Parameters params = pae_make_params(/*timeStep*/ 2.0, /*itrRestartFrom*/ 0, /*timeRestartFrom*/ 0.0);

    counters.events2D[0] = 1; // one dimerisation on the membrane so far

    {
        std::ofstream outFile(kPaeTmpFile, std::ios::out | std::ios::trunc);
        EXPECT_TRUE(outFile.is_open()) << "Could not open temporary output file";

        // First snapshot at iteration 100.
        print_association_events(counters, outFile, 100, params);

        // Simulation continues: another dimer forms and a larger complex binds.
        counters.events2D[0] = 2;
        counters.events2D[3] = 1;

        // Second snapshot at iteration 200.
        print_association_events(counters, outFile, 200, params);
    }

    std::ifstream inFile(kPaeTmpFile);
    std::stringstream buffer;
    buffer << inFile.rdbuf();
    inFile.close();
    std::remove(kPaeTmpFile);

    const std::vector<std::string> lines = pae_split_lines(buffer.str());
    pae_echo(lines);

    // Expected: header, "2D dimers: 1", header, "2D dimers: 2", "2D n=3: 1".
    EXPECT_EQ(lines.size(), 5u) << "Two headers plus one then two event lines expected";

    const int firstHeader = pae_line_index(lines, pae_expected_time_line(100, params));
    const int secondHeader = pae_line_index(lines, pae_expected_time_line(200, params));
    EXPECT_GE(firstHeader, 0) << "Header for iteration 100 missing";
    EXPECT_GE(secondHeader, 0) << "Header for iteration 200 missing";
    if (firstHeader >= 0 && secondHeader >= 0)
        EXPECT_LT(firstHeader, secondHeader) << "Snapshots must appear in call order";

    EXPECT_TRUE(pae_has_line(lines, "2D dimers: 1")) << "First snapshot should report 1 dimer";
    EXPECT_TRUE(pae_has_line(lines, "2D dimers: 2")) << "Second snapshot should report 2 dimers";
    EXPECT_TRUE(pae_has_line(lines, "2D n=3: 1"))
        << "Second snapshot should report the newly populated n=3 bin";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is invoked from its own TEST so that a
// failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(PrintAssociationEvents, HeaderOnlyWhenNoEvents) { test_pae_header_only_when_no_events(); }
TEST(PrintAssociationEvents, TimeHeaderValue) { test_pae_time_header_value(); }
TEST(PrintAssociationEvents, ThreeDBinLabels) { test_pae_3d_bin_labels(); }
TEST(PrintAssociationEvents, SkipsNonPositiveBins) { test_pae_skips_nonpositive_bins(); }
TEST(PrintAssociationEvents, ChannelPrefixesAndOrdering) { test_pae_channel_prefixes_and_ordering(); }
TEST(PrintAssociationEvents, RespectsEventArraySize) { test_pae_respects_eventArraySize(); }
TEST(PrintAssociationEvents, RepeatedCallsAppend) { test_pae_repeated_calls_append(); }