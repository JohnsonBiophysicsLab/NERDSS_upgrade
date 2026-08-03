/*! \file test_write_all_species.cpp
 *
 * ### Unit test for ../src/io/write_all_species.cpp
 *
 * The file under test defines a single function:
 *
 * \code
 *   void write_all_species(double simTime, std::ofstream& speciesFile,
 *                         const copyCounters& counterArray);
 * \endcode
 *
 * The function writes one comma separated line to the supplied output file
 * stream:
 *
 *   <simTime>,<copyNumSpecies[0]>,<copyNumSpecies[1]>,...\n
 *
 * i.e. the simulation time first (with no leading comma), then every element of
 * copyCounters::copyNumSpecies each prefixed by a comma, then std::endl.
 *
 * Because the function only takes an std::ofstream (not a generic std::ostream)
 * these tests write to a temporary file on disk and read the contents back to
 * verify the exact bytes produced.
 *
 * NOTE: we deliberately declare the prototype of the function under test here
 * instead of including "io/io.hpp".  io.hpp transitively pulls in
 * class_SimulVolume.hpp which `#include`s split.cpp; including that from more
 * than one translation unit of the combined gtest binary can produce duplicate
 * symbols.  The declaration below matches the definition exactly, so the test
 * still links against the real implementation.
 */

#include "classes/class_copyCounters.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Prototype of the function under test (defined in src/io/write_all_species.cpp).
void write_all_species(double simTime, std::ofstream& speciesFile,
    const copyCounters& counterArray);

namespace {

//! Name of the scratch file used by every test in this file.
const char* kWasTmpFileName = "test_write_all_species_output.tmp";

/*! \brief Read the whole scratch file back into a std::string.
 *
 * \return The exact contents of the file (empty string if it cannot be opened).
 */
std::string was_read_tmp_file()
{
    std::ifstream in(kWasTmpFileName);
    if (!in.is_open()) {
        std::cerr << "  !! could not reopen temporary file for reading\n";
        return std::string {};
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

/*! \brief Call write_all_species() once and return what landed on disk.
 *
 * The file is truncated first, so the returned string contains exactly the one
 * line written by the function under test.
 */
std::string was_write_one_line(double simTime, const copyCounters& counters)
{
    std::ofstream out(kWasTmpFileName, std::ios::out | std::ios::trunc);
    EXPECT_TRUE(out.is_open()) << "temporary output file must be writable";
    if (out.is_open()) {
        write_all_species(simTime, out, counters);
        out.close();
    }
    return was_read_tmp_file();
}

/*! \brief Build the string we expect write_all_species() to produce.
 *
 * The expectation is built with an std::ostringstream that has the same default
 * formatting flags as the std::ofstream used by the function under test, so the
 * floating point rendering of simTime is guaranteed to match without hardcoding
 * iostream default precision behaviour into the test.
 */
std::string was_expected_line(double simTime, const std::vector<int>& species)
{
    std::ostringstream oss;
    oss << simTime;
    for (int elem : species)
        oss << ',' << elem;
    oss << '\n'; // std::endl -> newline (+ flush, which leaves no extra bytes)
    return oss.str();
}

//! Remove the scratch file so tests do not leave litter behind.
void was_cleanup_tmp_file() { std::remove(kWasTmpFileName); }

} // namespace

// -----------------------------------------------------------------------------
// Test 1: A counter object with no species at all.
//         Expected: only the time value and a newline; no commas whatsoever.
// -----------------------------------------------------------------------------
void test_was_empty_species_list()
{
    std::cerr << "\n[TEST] test_was_empty_species_list\n"
              << "  Source file:   src/io/write_all_species.cpp\n"
              << "  Function:      write_all_species()\n"
              << "  Scenario:      copyCounters::copyNumSpecies is empty.\n"
              << "  Pass criteria: file holds just \"<simTime>\\n\" with no comma.\n";

    copyCounters counters; // copyNumSpecies default-constructed (empty)
    const double simTime = 0.0;

    const std::string actual = was_write_one_line(simTime, counters);
    const std::string expected = was_expected_line(simTime, counters.copyNumSpecies);

    std::cerr << "  Written line: \"" << actual << "\" (expected \"" << expected << "\")\n";

    EXPECT_EQ(actual, expected) << "empty species list must yield only the time value";
    EXPECT_EQ(actual.find(','), std::string::npos)
        << "no comma may appear when there are no species";
    EXPECT_FALSE(actual.empty()) << "a line must still be written";
    if (!actual.empty()) {
        EXPECT_EQ(actual.back(), '\n') << "the record must be newline terminated";
    }

    was_cleanup_tmp_file();
}

// -----------------------------------------------------------------------------
// Test 2: A typical populated species list.
//         Expected: time, then each count in order, comma separated.
// -----------------------------------------------------------------------------
void test_was_multiple_species()
{
    std::cerr << "\n[TEST] test_was_multiple_species\n"
              << "  Source file:   src/io/write_all_species.cpp\n"
              << "  Function:      write_all_species()\n"
              << "  Scenario:      five species counts written at simTime = 1.5.\n"
              << "  Pass criteria: output equals \"1.5,10,0,3,7,42\\n\" (order kept,\n"
              << "                 no trailing comma).\n";

    copyCounters counters;
    counters.copyNumSpecies = { 10, 0, 3, 7, 42 };
    const double simTime = 1.5;

    const std::string actual = was_write_one_line(simTime, counters);
    const std::string expected = was_expected_line(simTime, counters.copyNumSpecies);

    std::cerr << "  Written line: \"" << actual << "\"\n";
    std::cerr << "  Expected    : \"" << expected << "\"\n";

    EXPECT_EQ(actual, expected) << "species counts must be emitted in order";
    EXPECT_EQ(actual, std::string("1.5,10,0,3,7,42\n"))
        << "explicit byte-for-byte check of the expected CSV record";

    // There must be exactly one comma per species entry.
    const size_t commaCount = std::count(actual.begin(), actual.end(), ',');
    std::cerr << "  Comma count = " << commaCount << " (expected "
              << counters.copyNumSpecies.size() << ")\n";
    EXPECT_EQ(commaCount, counters.copyNumSpecies.size())
        << "one leading comma per species, none trailing";

    was_cleanup_tmp_file();
}

// -----------------------------------------------------------------------------
// Test 3: A single species -- checks the boundary between "no comma" and
//         "comma per element" logic.
// -----------------------------------------------------------------------------
void test_was_single_species()
{
    std::cerr << "\n[TEST] test_was_single_species\n"
              << "  Source file:   src/io/write_all_species.cpp\n"
              << "  Function:      write_all_species()\n"
              << "  Scenario:      exactly one species count.\n"
              << "  Pass criteria: output is \"<simTime>,<count>\\n\".\n";

    copyCounters counters;
    counters.copyNumSpecies = { 99 };
    const double simTime = 2.0;

    const std::string actual = was_write_one_line(simTime, counters);

    std::cerr << "  Written line: \"" << actual << "\"\n";

    EXPECT_EQ(actual, std::string("2,99\n"))
        << "one species produces exactly one comma separated value";

    was_cleanup_tmp_file();
}

// -----------------------------------------------------------------------------
// Test 4: Negative and large magnitude counts must be written verbatim (the
//         function performs no filtering or clamping).
// -----------------------------------------------------------------------------
void test_was_negative_and_large_counts()
{
    std::cerr << "\n[TEST] test_was_negative_and_large_counts\n"
              << "  Source file:   src/io/write_all_species.cpp\n"
              << "  Function:      write_all_species()\n"
              << "  Scenario:      counts include a negative value and INT_MAX-ish value.\n"
              << "  Pass criteria: values are streamed verbatim, unmodified.\n";

    copyCounters counters;
    counters.copyNumSpecies = { -5, 0, 2147483647 };
    const double simTime = 12.0;

    const std::string actual = was_write_one_line(simTime, counters);
    const std::string expected = was_expected_line(simTime, counters.copyNumSpecies);

    std::cerr << "  Written line: \"" << actual << "\"\n";

    EXPECT_EQ(actual, expected) << "values must not be altered by the writer";
    EXPECT_NE(actual.find("-5"), std::string::npos) << "negative count must appear";
    EXPECT_NE(actual.find("2147483647"), std::string::npos) << "large count must appear";

    was_cleanup_tmp_file();
}

// -----------------------------------------------------------------------------
// Test 5: Formatting of the simulation time.  We compare against an
//         ostringstream using the same default flags, and additionally spot
//         check a few representative time values.
// -----------------------------------------------------------------------------
void test_was_time_formatting()
{
    std::cerr << "\n[TEST] test_was_time_formatting\n"
              << "  Source file:   src/io/write_all_species.cpp\n"
              << "  Function:      write_all_species()\n"
              << "  Scenario:      several simTime values (integral, fractional, tiny).\n"
              << "  Pass criteria: simTime is streamed with default ostream formatting\n"
              << "                 and is the first field of the record.\n";

    copyCounters counters;
    counters.copyNumSpecies = { 1, 2 };

    const std::vector<double> times { 0.0, 1.0, 0.125, 1.0e-6, 123456.75 };

    for (double simTime : times) {
        const std::string actual = was_write_one_line(simTime, counters);
        const std::string expected = was_expected_line(simTime, counters.copyNumSpecies);

        std::cerr << "  simTime = " << simTime << " -> \"" << actual << "\"\n";

        EXPECT_EQ(actual, expected)
            << "time field must match default ostream formatting for " << simTime;

        // The time must be the very first field (nothing before the first comma
        // other than the time itself).
        const size_t firstComma = actual.find(',');
        EXPECT_NE(firstComma, std::string::npos) << "there should be species after the time";
        if (firstComma != std::string::npos) {
            std::ostringstream timeOnly;
            timeOnly << simTime;
            EXPECT_EQ(actual.substr(0, firstComma), timeOnly.str())
                << "first CSV field must be the simulation time";
        }
    }

    was_cleanup_tmp_file();
}

// -----------------------------------------------------------------------------
// Test 6: Repeated calls on the same open stream must append one line per call
//         (the function writes exactly one record and does not rewind).
// -----------------------------------------------------------------------------
void test_was_multiple_calls_append_rows()
{
    std::cerr << "\n[TEST] test_was_multiple_calls_append_rows\n"
              << "  Source file:   src/io/write_all_species.cpp\n"
              << "  Function:      write_all_species()\n"
              << "  Scenario:      three consecutive calls on one open ofstream,\n"
              << "                 with the counts changing between calls.\n"
              << "  Pass criteria: file contains three lines in call order.\n";

    copyCounters counters;
    counters.copyNumSpecies = { 1, 1 };

    {
        std::ofstream out(kWasTmpFileName, std::ios::out | std::ios::trunc);
        EXPECT_TRUE(out.is_open()) << "temporary output file must be writable";

        write_all_species(0.0, out, counters);

        // Change the counts to make sure the second row reflects new state.
        counters.copyNumSpecies = { 2, 3 };
        write_all_species(0.5, out, counters);

        // Shrink the species list for the third row.
        counters.copyNumSpecies = { 4 };
        write_all_species(1.0, out, counters);

        out.close();
    }

    const std::string contents = was_read_tmp_file();
    std::cerr << "  File contents:\n" << contents;

    // Split into lines for verification.
    std::vector<std::string> lines;
    std::istringstream iss(contents);
    std::string line;
    while (std::getline(iss, line))
        lines.push_back(line);

    std::cerr << "  Line count = " << lines.size() << " (expected 3)\n";
    EXPECT_EQ(lines.size(), static_cast<size_t>(3))
        << "one record must be written per call";

    if (lines.size() >= 3) {
        EXPECT_EQ(lines[0], std::string("0,1,1")) << "first record";
        EXPECT_EQ(lines[1], std::string("0.5,2,3")) << "second record";
        EXPECT_EQ(lines[2], std::string("1,4")) << "third record (shorter species list)";
    }

    was_cleanup_tmp_file();
}

// -----------------------------------------------------------------------------
// Test 7: Only copyNumSpecies is consulted.  Filling in every other member of
//         copyCounters must not change the output at all.
// -----------------------------------------------------------------------------
void test_was_ignores_other_counters()
{
    std::cerr << "\n[TEST] test_was_ignores_other_counters\n"
              << "  Source file:   src/io/write_all_species.cpp\n"
              << "  Function:      write_all_species()\n"
              << "  Scenario:      two copyCounters with identical copyNumSpecies but\n"
              << "                 wildly different values for every other member.\n"
              << "  Pass criteria: both produce byte-identical output.\n";

    // Reference object: only copyNumSpecies populated.
    copyCounters plain;
    plain.copyNumSpecies = { 5, 6, 7 };

    // Noisy object: same species counts, everything else filled in.
    copyCounters noisy;
    noisy.copyNumSpecies = { 5, 6, 7 };
    noisy.nBoundPairs = { 11, 12, 13 };
    noisy.proPairlist = { 1, 2 };
    noisy.nLoops = 4;
    noisy.nCancelOverlapPartner = 5;
    noisy.nCancelOverlapSystem = 6;
    noisy.nCancelDisplace2D = 7;
    noisy.nCancelDisplace3D = 8;
    noisy.nCancelDisplace3Dto2D = 9;
    noisy.nCancelSpanBox = 10;
    noisy.nAssocSuccess = 11;
    noisy.singleDouble = { 0, 1, 0 };
    noisy.implicitDouble = { true, false, true };
    noisy.canDissociate = { false, true, false };
    noisy.bindPairList = { { 1, 2 }, { 3 } };
    noisy.events3D = { 1, 2, 3 };
    noisy.events2D = { 4, 5 };
    noisy.events3Dto2D = { 6 };

    const double simTime = 3.25;

    const std::string plainOut = was_write_one_line(simTime, plain);
    const std::string noisyOut = was_write_one_line(simTime, noisy);

    std::cerr << "  plain  -> \"" << plainOut << "\"\n";
    std::cerr << "  noisy  -> \"" << noisyOut << "\"\n";

    EXPECT_EQ(plainOut, noisyOut)
        << "only copyNumSpecies may influence the written record";
    EXPECT_EQ(noisyOut, std::string("3.25,5,6,7\n"))
        << "explicit expected record for the noisy counter object";

    was_cleanup_tmp_file();
}

// -----------------------------------------------------------------------------
// Test 8: The stream must remain healthy (good()) after the call so callers can
//         keep appending, and the record must be flushed (std::endl) so the
//         bytes are visible to a separate reader even before close().
// -----------------------------------------------------------------------------
void test_was_stream_state_and_flush()
{
    std::cerr << "\n[TEST] test_was_stream_state_and_flush\n"
              << "  Source file:   src/io/write_all_species.cpp\n"
              << "  Function:      write_all_species()\n"
              << "  Scenario:      inspect the stream and the file while the stream is\n"
              << "                 still open (std::endl should flush).\n"
              << "  Pass criteria: stream still good() and the bytes are already on disk.\n";

    copyCounters counters;
    counters.copyNumSpecies = { 8, 9 };
    const double simTime = 7.0;

    {
        std::ofstream out(kWasTmpFileName, std::ios::out | std::ios::trunc);
        EXPECT_TRUE(out.is_open()) << "temporary output file must be writable";

        write_all_species(simTime, out, counters);

        // std::endl flushes, so the data should be readable before close().
        EXPECT_TRUE(out.good()) << "stream must remain usable after the write";

        const std::string flushed = was_read_tmp_file();
        std::cerr << "  Contents before close(): \"" << flushed << "\"\n";
        EXPECT_EQ(flushed, std::string("7,8,9\n"))
            << "std::endl should have flushed the record to disk";

        out.close();
    }

    // And of course the same content is there after closing.
    const std::string finalContents = was_read_tmp_file();
    std::cerr << "  Contents after close():  \"" << finalContents << "\"\n";
    EXPECT_EQ(finalContents, std::string("7,8,9\n")) << "content unchanged by close()";

    was_cleanup_tmp_file();
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* helper runs inside its own TEST so the
// framework reports individual results, and all EXPECT_* (non-fatal) assertions
// let every test complete even if some checks fail.
// -----------------------------------------------------------------------------
TEST(WriteAllSpecies, EmptySpeciesList) { test_was_empty_species_list(); }
TEST(WriteAllSpecies, MultipleSpecies) { test_was_multiple_species(); }
TEST(WriteAllSpecies, SingleSpecies) { test_was_single_species(); }
TEST(WriteAllSpecies, NegativeAndLargeCounts) { test_was_negative_and_large_counts(); }
TEST(WriteAllSpecies, TimeFormatting) { test_was_time_formatting(); }
TEST(WriteAllSpecies, MultipleCallsAppendRows) { test_was_multiple_calls_append_rows(); }
TEST(WriteAllSpecies, IgnoresOtherCounters) { test_was_ignores_other_counters(); }
TEST(WriteAllSpecies, StreamStateAndFlush) { test_was_stream_state_and_flush(); }