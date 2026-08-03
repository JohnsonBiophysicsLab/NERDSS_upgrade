/*! \file test_track_association_events.cpp
 *
 * ### Unit test for src/reactions/track_association_events.cpp
 *
 * Function under test:
 *
 *     void track_association_events(Complex& reactCom1, Complex& reactCom2,
 *                                   bool transitionToSurface, bool isOnMembrane,
 *                                   copyCounters& counterArrays)
 *
 * The routine is a pure book-keeping function. Given the two complexes that
 * are about to associate, it works out
 *
 *   1. the size (number of member molecules) of the *smaller* of the two
 *      complexes,
 *   2. maps that size onto a bin index in the event histogram, and
 *   3. increments exactly one bin in exactly one of the three histograms
 *      (events3D, events2D, events3Dto2D) depending on the dimensionality
 *      flags passed in.
 *
 * The binning rules implemented by the source are:
 *
 *   - 1 + 1                        -> index 0   (a "dimerization" event)
 *   - smaller size n, 1 <= n <= 9  -> index n
 *   - smaller size n, n >= 10      -> index 9 + (n / 10)   (integer division)
 *   - the index is clamped to (eventArraySize - 1)
 *
 * The dimensionality selection rules are:
 *
 *   - transitionToSurface == true             -> events3Dto2D  (takes priority)
 *   - transitionToSurface == false && isOnMembrane == true -> events2D
 *   - otherwise                               -> events3D
 *
 * Every test below prints what it is doing, which bin it expects, and why.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <numeric>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (uniquely prefixed with "tae_" so they cannot collide with
// helpers defined by other translation units in the combined test binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a Complex whose memberList contains \p numMembers dummy entries.
 *
 * track_association_events() only ever reads memberList.size(), so the actual
 * molecule indices stored are irrelevant; we fill them with 0..numMembers-1 for
 * readability while debugging.
 *
 * \param[in] numMembers Number of member molecules the complex should report.
 * \return A Complex with the requested member count.
 */
Complex tae_make_complex(int numMembers)
{
    Complex com;
    com.memberList.clear();
    for (int i = 0; i < numMembers; ++i)
        com.memberList.push_back(i);
    return com;
}

/*! \brief Create a zeroed copyCounters object with event histograms allocated.
 *
 * copyCounters defaults eventArraySize to 20 but leaves the three histogram
 * vectors empty, so the caller (normally init_association_events()) must size
 * them. We do that here explicitly so this test does not depend on any other
 * source file.
 *
 * \param[in] arraySize Number of bins in each of the three histograms.
 * \return A copyCounters with events3D/events2D/events3Dto2D all zeroed.
 */
copyCounters tae_make_counters(int arraySize = 20)
{
    copyCounters counterArrays;
    counterArrays.eventArraySize = arraySize;
    counterArrays.events3D.assign(arraySize, 0);
    counterArrays.events2D.assign(arraySize, 0);
    counterArrays.events3Dto2D.assign(arraySize, 0);
    return counterArrays;
}

/*! \brief Sum every bin of a histogram; used to prove only one bin changed. */
int tae_sum(const std::vector<int>& hist)
{
    return std::accumulate(hist.begin(), hist.end(), 0);
}

/*! \brief Verify that exactly one bin of one histogram was incremented.
 *
 * \param[in] counterArrays   The counters after the call under test.
 * \param[in] expectHist      Pointer to the histogram expected to be hit.
 * \param[in] expectIndex     The bin expected to hold the single increment.
 * \param[in] otherA          One of the two histograms expected to stay zero.
 * \param[in] otherB          The other histogram expected to stay zero.
 * \param[in] label           Human readable description for the log output.
 */
void tae_expect_single_hit(const std::vector<int>& expectHist, int expectIndex,
    const std::vector<int>& otherA, const std::vector<int>& otherB,
    const std::string& label)
{
    std::cerr << "    checking that " << label << " bin " << expectIndex
              << " holds the only increment\n";

    // The expected bin must contain exactly one event.
    EXPECT_EQ(expectHist.at(expectIndex), 1)
        << label << " bin " << expectIndex << " should have been incremented once";

    // No other bin of the same histogram may have been touched.
    EXPECT_EQ(tae_sum(expectHist), 1)
        << label << " should contain exactly one event in total";

    // The other two histograms must be completely untouched.
    EXPECT_EQ(tae_sum(otherA), 0) << "a non-selected histogram was modified";
    EXPECT_EQ(tae_sum(otherB), 0) << "a non-selected histogram was modified";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: 1 + 1 association in 3D is recorded as a dimerization in bin 0.
// -----------------------------------------------------------------------------
void test_tae_dimerization_in_3d()
{
    std::cerr << "\n[TEST] test_tae_dimerization_in_3d\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    monomer (1 member) + monomer (1 member), in solution.\n"
              << "  Criteria:    events3D[0] == 1 (bin 0 is reserved for 1+1),\n"
              << "               events2D and events3Dto2D untouched.\n";

    copyCounters counterArrays = tae_make_counters();
    Complex com1 = tae_make_complex(1);
    Complex com2 = tae_make_complex(1);

    // transitionToSurface = false, isOnMembrane = false -> pure 3D event.
    track_association_events(com1, com2, false, false, counterArrays);

    tae_expect_single_hit(counterArrays.events3D, 0, counterArrays.events2D,
        counterArrays.events3Dto2D, "events3D");

    // The function must not disturb the complexes it inspects.
    EXPECT_EQ(com1.memberList.size(), 1u) << "reactCom1 memberList must not change";
    EXPECT_EQ(com2.memberList.size(), 1u) << "reactCom2 memberList must not change";
}

// -----------------------------------------------------------------------------
// Test 2: monomer joining a larger complex is binned by the *smaller* size (1).
// -----------------------------------------------------------------------------
void test_tae_monomer_joins_larger_complex()
{
    std::cerr << "\n[TEST] test_tae_monomer_joins_larger_complex\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    1 member + 5 members (and the mirrored 5 + 1).\n"
              << "  Criteria:    index == smaller size == 1, so events3D[1] == 1,\n"
              << "               and the result is independent of argument order.\n";

    // Case A: small complex passed first.
    {
        copyCounters counterArrays = tae_make_counters();
        Complex small = tae_make_complex(1);
        Complex large = tae_make_complex(5);
        std::cerr << "    case A: (1, 5)\n";
        track_association_events(small, large, false, false, counterArrays);
        tae_expect_single_hit(counterArrays.events3D, 1, counterArrays.events2D,
            counterArrays.events3Dto2D, "events3D");
    }

    // Case B: arguments swapped -- the smaller size still selects the bin.
    {
        copyCounters counterArrays = tae_make_counters();
        Complex large = tae_make_complex(5);
        Complex small = tae_make_complex(1);
        std::cerr << "    case B: (5, 1) -- expect the same bin as case A\n";
        track_association_events(large, small, false, false, counterArrays);
        tae_expect_single_hit(counterArrays.events3D, 1, counterArrays.events2D,
            counterArrays.events3Dto2D, "events3D");
    }
}

// -----------------------------------------------------------------------------
// Test 3: for smaller sizes 2..9 the bin index equals the smaller size.
// -----------------------------------------------------------------------------
void test_tae_small_sizes_map_directly()
{
    std::cerr << "\n[TEST] test_tae_small_sizes_map_directly\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    smaller complex of size n for n = 2..9, larger = n + 7.\n"
              << "  Criteria:    events3D[n] == 1 for each n (no bundling below 10).\n";

    for (int n = 2; n <= 9; ++n) {
        copyCounters counterArrays = tae_make_counters();
        Complex small = tae_make_complex(n);
        Complex large = tae_make_complex(n + 7);

        std::cerr << "    smaller size " << n << " + larger size " << (n + 7)
                  << " -> expecting bin " << n << '\n';
        track_association_events(small, large, false, false, counterArrays);

        EXPECT_EQ(counterArrays.events3D.at(n), 1)
            << "smaller size " << n << " should land in bin " << n;
        EXPECT_EQ(tae_sum(counterArrays.events3D), 1)
            << "only one bin should be incremented for smaller size " << n;
    }
}

// -----------------------------------------------------------------------------
// Test 4: sizes of 10 and above are bundled with a spacing of 10.
// -----------------------------------------------------------------------------
void test_tae_large_sizes_are_bundled()
{
    std::cerr << "\n[TEST] test_tae_large_sizes_are_bundled\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    smaller sizes >= 10 are grouped: index = 9 + n/10.\n"
              << "  Criteria:    10..19 -> bin 10, 20..29 -> bin 11, 35 -> bin 12.\n";

    // Each pair is {smaller complex size, expected histogram bin}.
    const std::vector<std::pair<int, int>> cases {
        { 10, 10 }, // 9 + 10/10 = 10
        { 15, 10 }, // 9 + 15/10 = 10 (same bundle as 10)
        { 19, 10 }, // 9 + 19/10 = 10 (top of the bundle)
        { 20, 11 }, // 9 + 20/10 = 11
        { 29, 11 }, // 9 + 29/10 = 11
        { 35, 12 }, // 9 + 35/10 = 12
    };

    for (const auto& oneCase : cases) {
        const int smallerSize = oneCase.first;
        const int expectedBin = oneCase.second;

        copyCounters counterArrays = tae_make_counters();
        Complex small = tae_make_complex(smallerSize);
        // Larger partner is deliberately much bigger so it never wins the min().
        Complex large = tae_make_complex(smallerSize + 100);

        std::cerr << "    smaller size " << smallerSize << " -> expecting bin "
                  << expectedBin << '\n';
        track_association_events(small, large, false, false, counterArrays);

        EXPECT_EQ(counterArrays.events3D.at(expectedBin), 1)
            << "smaller size " << smallerSize << " should be bundled into bin "
            << expectedBin;
        EXPECT_EQ(tae_sum(counterArrays.events3D), 1)
            << "exactly one bin should be incremented for smaller size "
            << smallerSize;
    }
}

// -----------------------------------------------------------------------------
// Test 5: enormous complexes are clamped into the final bin (no overflow).
// -----------------------------------------------------------------------------
void test_tae_index_is_clamped_to_last_bin()
{
    std::cerr << "\n[TEST] test_tae_index_is_clamped_to_last_bin\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    smaller size 110 with a 20-bin histogram would give\n"
              << "               index 9 + 11 = 20, which is out of range.\n"
              << "  Criteria:    the index is clamped to eventArraySize - 1 == 19,\n"
              << "               so events3D[19] == 1 and nothing overflows.\n";

    copyCounters counterArrays = tae_make_counters(20);
    Complex small = tae_make_complex(110);
    Complex large = tae_make_complex(500);

    track_association_events(small, large, false, false, counterArrays);

    tae_expect_single_hit(counterArrays.events3D, 19, counterArrays.events2D,
        counterArrays.events3Dto2D, "events3D");

    // The histogram must not have been resized by the call.
    EXPECT_EQ(static_cast<int>(counterArrays.events3D.size()), 20)
        << "events3D should still contain exactly eventArraySize bins";
}

// -----------------------------------------------------------------------------
// Test 6: the clamp respects a non-default eventArraySize.
// -----------------------------------------------------------------------------
void test_tae_clamp_respects_custom_array_size()
{
    std::cerr << "\n[TEST] test_tae_clamp_respects_custom_array_size\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    eventArraySize shrunk to 12; smaller size 50 gives\n"
              << "               raw index 9 + 5 = 14, which exceeds the array.\n"
              << "  Criteria:    the event is recorded in the last bin, index 11.\n";

    copyCounters counterArrays = tae_make_counters(12);
    Complex small = tae_make_complex(50);
    Complex large = tae_make_complex(80);

    track_association_events(small, large, false, false, counterArrays);

    tae_expect_single_hit(counterArrays.events3D, 11, counterArrays.events2D,
        counterArrays.events3Dto2D, "events3D");
}

// -----------------------------------------------------------------------------
// Test 7: the dimensionality flags select the correct histogram.
// -----------------------------------------------------------------------------
void test_tae_dimensionality_selection()
{
    std::cerr << "\n[TEST] test_tae_dimensionality_selection\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    same 1+1 event recorded with each flag combination.\n"
              << "  Criteria:    (false,false) -> events3D,\n"
              << "               (false,true)  -> events2D,\n"
              << "               (true,false)  -> events3Dto2D,\n"
              << "               (true,true)   -> events3Dto2D (transition wins).\n";

    // 3D: neither flag set.
    {
        copyCounters counterArrays = tae_make_counters();
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(1);
        std::cerr << "    transitionToSurface=false, isOnMembrane=false\n";
        track_association_events(a, b, false, false, counterArrays);
        tae_expect_single_hit(counterArrays.events3D, 0, counterArrays.events2D,
            counterArrays.events3Dto2D, "events3D");
    }

    // 2D: already on the membrane, not transitioning.
    {
        copyCounters counterArrays = tae_make_counters();
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(1);
        std::cerr << "    transitionToSurface=false, isOnMembrane=true\n";
        track_association_events(a, b, false, true, counterArrays);
        tae_expect_single_hit(counterArrays.events2D, 0, counterArrays.events3D,
            counterArrays.events3Dto2D, "events2D");
    }

    // 3D -> 2D: transitioning onto the surface.
    {
        copyCounters counterArrays = tae_make_counters();
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(1);
        std::cerr << "    transitionToSurface=true, isOnMembrane=false\n";
        track_association_events(a, b, true, false, counterArrays);
        tae_expect_single_hit(counterArrays.events3Dto2D, 0, counterArrays.events3D,
            counterArrays.events2D, "events3Dto2D");
    }

    // Both flags set: transitionToSurface is checked first, so it wins.
    {
        copyCounters counterArrays = tae_make_counters();
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(1);
        std::cerr << "    transitionToSurface=true, isOnMembrane=true"
                  << " (transition must take priority)\n";
        track_association_events(a, b, true, true, counterArrays);
        tae_expect_single_hit(counterArrays.events3Dto2D, 0, counterArrays.events3D,
            counterArrays.events2D, "events3Dto2D");
    }
}

// -----------------------------------------------------------------------------
// Test 8: repeated calls accumulate rather than overwrite.
// -----------------------------------------------------------------------------
void test_tae_repeated_calls_accumulate()
{
    std::cerr << "\n[TEST] test_tae_repeated_calls_accumulate\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    three 1+1 3D events, two 1+4 2D events and one\n"
              << "               3D->2D event, all pushed into the same counters.\n"
              << "  Criteria:    events3D[0]==3, events2D[1]==2, events3Dto2D[0]==1\n"
              << "               and the grand total of all bins is 6.\n";

    copyCounters counterArrays = tae_make_counters();

    // Three dimerizations in solution -> events3D[0] should reach 3.
    for (int i = 0; i < 3; ++i) {
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(1);
        track_association_events(a, b, false, false, counterArrays);
    }

    // Two "monomer adds to tetramer" events on the membrane -> events2D[1] == 2.
    for (int i = 0; i < 2; ++i) {
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(4);
        track_association_events(a, b, false, true, counterArrays);
    }

    // One dimerization that lands on the surface -> events3Dto2D[0] == 1.
    {
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(1);
        track_association_events(a, b, true, false, counterArrays);
    }

    std::cerr << "    events3D[0]=" << counterArrays.events3D[0]
              << ", events2D[1]=" << counterArrays.events2D[1]
              << ", events3Dto2D[0]=" << counterArrays.events3Dto2D[0] << '\n';

    EXPECT_EQ(counterArrays.events3D.at(0), 3)
        << "three 3D dimerizations should accumulate in events3D[0]";
    EXPECT_EQ(counterArrays.events2D.at(1), 2)
        << "two 2D monomer-addition events should accumulate in events2D[1]";
    EXPECT_EQ(counterArrays.events3Dto2D.at(0), 1)
        << "one 3D->2D dimerization should land in events3Dto2D[0]";

    // Nothing may leak into other bins: total across all three histograms is 6.
    const int total = tae_sum(counterArrays.events3D) + tae_sum(counterArrays.events2D)
        + tae_sum(counterArrays.events3Dto2D);
    EXPECT_EQ(total, 6) << "exactly six events were recorded in total";
}

// -----------------------------------------------------------------------------
// Test 9: a 1 + 1 event is *not* confused with the "smaller size 1" bin.
// -----------------------------------------------------------------------------
void test_tae_dimerization_is_distinct_from_monomer_addition()
{
    std::cerr << "\n[TEST] test_tae_dimerization_is_distinct_from_monomer_addition\n"
              << "  Source file: src/reactions/track_association_events.cpp\n"
              << "  Function:    track_association_events\n"
              << "  Scenario:    one 1+1 event and one 1+2 event in the same counters.\n"
              << "  Criteria:    the dimerization goes to bin 0 while the monomer\n"
              << "               addition goes to bin 1 -- they must not share a bin.\n";

    copyCounters counterArrays = tae_make_counters();

    // 1 + 1 -> bin 0 (dimerization special case).
    {
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(1);
        track_association_events(a, b, false, false, counterArrays);
    }
    // 1 + 2 -> bin 1 (smaller size is 1, but this is not a dimerization).
    {
        Complex a = tae_make_complex(1);
        Complex b = tae_make_complex(2);
        track_association_events(a, b, false, false, counterArrays);
    }

    std::cerr << "    events3D[0]=" << counterArrays.events3D[0]
              << " (expect 1), events3D[1]=" << counterArrays.events3D[1]
              << " (expect 1)\n";

    EXPECT_EQ(counterArrays.events3D.at(0), 1)
        << "the 1+1 event belongs in the dedicated dimerization bin 0";
    EXPECT_EQ(counterArrays.events3D.at(1), 1)
        << "the 1+2 event belongs in bin 1 (smaller size == 1)";
    EXPECT_EQ(tae_sum(counterArrays.events3D), 2)
        << "only the two expected bins should have been incremented";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is executed inside its own TEST so a
// failure in one scenario does not stop the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(TrackAssociationEvents, DimerizationIn3D) { test_tae_dimerization_in_3d(); }
TEST(TrackAssociationEvents, MonomerJoinsLargerComplex) { test_tae_monomer_joins_larger_complex(); }
TEST(TrackAssociationEvents, SmallSizesMapDirectly) { test_tae_small_sizes_map_directly(); }
TEST(TrackAssociationEvents, LargeSizesAreBundled) { test_tae_large_sizes_are_bundled(); }
TEST(TrackAssociationEvents, IndexIsClampedToLastBin) { test_tae_index_is_clamped_to_last_bin(); }
TEST(TrackAssociationEvents, ClampRespectsCustomArraySize) { test_tae_clamp_respects_custom_array_size(); }
TEST(TrackAssociationEvents, DimensionalitySelection) { test_tae_dimensionality_selection(); }
TEST(TrackAssociationEvents, RepeatedCallsAccumulate) { test_tae_repeated_calls_accumulate(); }
TEST(TrackAssociationEvents, DimerizationDistinctFromMonomerAddition)
{
    test_tae_dimerization_is_distinct_from_monomer_addition();
}