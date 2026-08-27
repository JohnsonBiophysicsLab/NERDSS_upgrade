/*! \file test_init_association_events.cpp
 *
 * ### Unit test for ../src/reactions/init_association_events.cpp
 *
 * This file exercises the single function defined in that source file:
 *
 *     void init_association_events(copyCounters& counterArrays)
 *
 * Behaviour of the function (read directly from the implementation):
 *   - It reads `counterArrays.eventArraySize` (call it `arraySize`).
 *   - It then loops `arraySize` times and, on each iteration, *push_back*s a
 *     single 0 onto each of the three histogram vectors:
 *         counterArrays.events3D
 *         counterArrays.events3Dto2D
 *         counterArrays.events2D
 *   - It does NOT clear the vectors first, so calling it twice appends a
 *     second block of zeros rather than re-initializing.
 *   - The local variables `spacing` and `maxSingles` are unused by the loop,
 *     so they cannot be observed from outside; nothing else in copyCounters is
 *     touched.
 *
 * The tests below therefore verify: resulting sizes, that every element is
 * zero, that unrelated copyCounters members are untouched, the append (not
 * reset) semantics on a repeated call, and the degenerate zero/negative
 * `eventArraySize` cases (where the loop body never executes).
 */

#include "classes/class_copyCounters.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Small local helper: report whether every element of a vector<int> is zero and
// print a verbose message describing what is being checked.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Check that every entry of a histogram vector is exactly 0.
 *
 * \param[in] vec   The histogram vector produced by init_association_events().
 * \param[in] name  Human readable name used for the console/assertion message.
 */
void iae_expect_all_zero(const std::vector<int>& vec, const char* name)
{
    std::cerr << "    checking that all " << vec.size() << " entries of " << name
              << " are initialized to 0\n";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        // Use EXPECT (non-fatal) so that a single bad element does not abort
        // the remainder of the suite.
        EXPECT_EQ(vec[i], 0) << name << "[" << i << "] should have been initialized to 0";
    }
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Default construction -- eventArraySize defaults to 20 in the
//         copyCounters class definition, so all three histograms must end up
//         with exactly 20 zeroed bins.
// -----------------------------------------------------------------------------
void test_iae_default_array_size()
{
    std::cerr << "\n[TEST] test_iae_default_array_size\n"
              << "  Source file:   src/reactions/init_association_events.cpp\n"
              << "  Function:      init_association_events(copyCounters&)\n"
              << "  Scenario:      freshly default-constructed copyCounters, whose\n"
              << "                 eventArraySize member defaults to 20.\n"
              << "  Pass criteria: events3D / events3Dto2D / events2D each have 20\n"
              << "                 elements and every element equals 0.\n";

    copyCounters counterArrays; // eventArraySize defaults to 20

    // Sanity check on the precondition we rely on for the expected size.
    std::cerr << "  Precondition: counterArrays.eventArraySize = "
              << counterArrays.eventArraySize << '\n';
    EXPECT_EQ(counterArrays.eventArraySize, 20)
        << "copyCounters::eventArraySize is expected to default to 20";

    // The three histograms should start out empty before the call.
    EXPECT_TRUE(counterArrays.events3D.empty()) << "events3D should start empty";
    EXPECT_TRUE(counterArrays.events3Dto2D.empty()) << "events3Dto2D should start empty";
    EXPECT_TRUE(counterArrays.events2D.empty()) << "events2D should start empty";

    std::cerr << "  Calling init_association_events...\n";
    init_association_events(counterArrays);

    const std::size_t expected = static_cast<std::size_t>(counterArrays.eventArraySize);
    std::cerr << "  After the call: events3D.size()=" << counterArrays.events3D.size()
              << ", events3Dto2D.size()=" << counterArrays.events3Dto2D.size()
              << ", events2D.size()=" << counterArrays.events2D.size()
              << " (expected " << expected << " for each)\n";

    EXPECT_EQ(counterArrays.events3D.size(), expected)
        << "events3D should contain eventArraySize bins";
    EXPECT_EQ(counterArrays.events3Dto2D.size(), expected)
        << "events3Dto2D should contain eventArraySize bins";
    EXPECT_EQ(counterArrays.events2D.size(), expected)
        << "events2D should contain eventArraySize bins";

    iae_expect_all_zero(counterArrays.events3D, "events3D");
    iae_expect_all_zero(counterArrays.events3Dto2D, "events3Dto2D");
    iae_expect_all_zero(counterArrays.events2D, "events2D");
}

// -----------------------------------------------------------------------------
// Test 2: A user-supplied, non-default eventArraySize must be honoured exactly
//         -- the number of bins is driven solely by that member.
// -----------------------------------------------------------------------------
void test_iae_custom_array_size()
{
    std::cerr << "\n[TEST] test_iae_custom_array_size\n"
              << "  Source file:   src/reactions/init_association_events.cpp\n"
              << "  Function:      init_association_events(copyCounters&)\n"
              << "  Scenario:      eventArraySize is overridden to 7 before the call.\n"
              << "  Pass criteria: each histogram has exactly 7 zeroed elements.\n";

    copyCounters counterArrays;
    counterArrays.eventArraySize = 7; // deliberately different from the default 20

    std::cerr << "  Calling init_association_events with eventArraySize = "
              << counterArrays.eventArraySize << " ...\n";
    init_association_events(counterArrays);

    std::cerr << "  After the call: events3D.size()=" << counterArrays.events3D.size()
              << ", events3Dto2D.size()=" << counterArrays.events3Dto2D.size()
              << ", events2D.size()=" << counterArrays.events2D.size() << '\n';

    EXPECT_EQ(counterArrays.events3D.size(), 7u)
        << "events3D size must follow the custom eventArraySize";
    EXPECT_EQ(counterArrays.events3Dto2D.size(), 7u)
        << "events3Dto2D size must follow the custom eventArraySize";
    EXPECT_EQ(counterArrays.events2D.size(), 7u)
        << "events2D size must follow the custom eventArraySize";

    iae_expect_all_zero(counterArrays.events3D, "events3D");
    iae_expect_all_zero(counterArrays.events3Dto2D, "events3Dto2D");
    iae_expect_all_zero(counterArrays.events2D, "events2D");

    // The function should not modify eventArraySize itself.
    EXPECT_EQ(counterArrays.eventArraySize, 7)
        << "init_association_events must not alter eventArraySize";
}

// -----------------------------------------------------------------------------
// Test 3: The three histograms must all have identical length after the call,
//         because the same loop pushes one element onto each per iteration.
// -----------------------------------------------------------------------------
void test_iae_all_three_histograms_match()
{
    std::cerr << "\n[TEST] test_iae_all_three_histograms_match\n"
              << "  Source file:   src/reactions/init_association_events.cpp\n"
              << "  Function:      init_association_events(copyCounters&)\n"
              << "  Scenario:      one shared loop fills 3D, 3D->2D and 2D histograms.\n"
              << "  Pass criteria: all three vectors end up with the same size.\n";

    copyCounters counterArrays;
    counterArrays.eventArraySize = 13;

    std::cerr << "  Calling init_association_events...\n";
    init_association_events(counterArrays);

    std::cerr << "  Sizes: 3D=" << counterArrays.events3D.size()
              << ", 3Dto2D=" << counterArrays.events3Dto2D.size()
              << ", 2D=" << counterArrays.events2D.size() << '\n';

    EXPECT_EQ(counterArrays.events3D.size(), counterArrays.events2D.size())
        << "events3D and events2D must be filled in lockstep";
    EXPECT_EQ(counterArrays.events3D.size(), counterArrays.events3Dto2D.size())
        << "events3D and events3Dto2D must be filled in lockstep";
}

// -----------------------------------------------------------------------------
// Test 4: Unrelated copyCounters members must be left completely untouched --
//         the function only appends to the three event histograms.
// -----------------------------------------------------------------------------
void test_iae_leaves_other_counters_untouched()
{
    std::cerr << "\n[TEST] test_iae_leaves_other_counters_untouched\n"
              << "  Source file:   src/reactions/init_association_events.cpp\n"
              << "  Function:      init_association_events(copyCounters&)\n"
              << "  Scenario:      pre-populate unrelated counters/vectors, then call.\n"
              << "  Pass criteria: those members retain their original values.\n";

    copyCounters counterArrays;

    // Seed a selection of unrelated members with recognizable values.
    counterArrays.nLoops = 42;
    counterArrays.nAssocSuccess = 7;
    counterArrays.nCancelSpanBox = 3;
    counterArrays.nCancelOverlapPartner = 5;
    counterArrays.nCancelOverlapSystem = 6;
    counterArrays.nCancelDisplace2D = 1;
    counterArrays.nCancelDisplace3D = 2;
    counterArrays.nCancelDisplace3Dto2D = 4;
    counterArrays.nBoundPairs = std::vector<int>{ 11, 22, 33 };
    counterArrays.copyNumSpecies = std::vector<int>{ 100, 200 };
    counterArrays.proPairlist = std::vector<int>{ 9 };

    std::cerr << "  Calling init_association_events...\n";
    init_association_events(counterArrays);

    // Scalar counters must be unchanged.
    EXPECT_EQ(counterArrays.nLoops, 42) << "nLoops must not be modified";
    EXPECT_EQ(counterArrays.nAssocSuccess, 7) << "nAssocSuccess must not be modified";
    EXPECT_EQ(counterArrays.nCancelSpanBox, 3) << "nCancelSpanBox must not be modified";
    EXPECT_EQ(counterArrays.nCancelOverlapPartner, 5)
        << "nCancelOverlapPartner must not be modified";
    EXPECT_EQ(counterArrays.nCancelOverlapSystem, 6)
        << "nCancelOverlapSystem must not be modified";
    EXPECT_EQ(counterArrays.nCancelDisplace2D, 1) << "nCancelDisplace2D must not be modified";
    EXPECT_EQ(counterArrays.nCancelDisplace3D, 2) << "nCancelDisplace3D must not be modified";
    EXPECT_EQ(counterArrays.nCancelDisplace3Dto2D, 4)
        << "nCancelDisplace3Dto2D must not be modified";

    // Vector members other than the three histograms must be unchanged.
    ASSERT_EQ(counterArrays.nBoundPairs.size(), 3u)
        << "nBoundPairs should keep its original size";
    EXPECT_EQ(counterArrays.nBoundPairs[0], 11) << "nBoundPairs[0] must not be modified";
    EXPECT_EQ(counterArrays.nBoundPairs[1], 22) << "nBoundPairs[1] must not be modified";
    EXPECT_EQ(counterArrays.nBoundPairs[2], 33) << "nBoundPairs[2] must not be modified";

    ASSERT_EQ(counterArrays.copyNumSpecies.size(), 2u)
        << "copyNumSpecies should keep its original size";
    EXPECT_EQ(counterArrays.copyNumSpecies[0], 100) << "copyNumSpecies[0] must not be modified";
    EXPECT_EQ(counterArrays.copyNumSpecies[1], 200) << "copyNumSpecies[1] must not be modified";

    ASSERT_EQ(counterArrays.proPairlist.size(), 1u)
        << "proPairlist should keep its original size";
    EXPECT_EQ(counterArrays.proPairlist[0], 9) << "proPairlist[0] must not be modified";

    std::cerr << "  Unrelated counters verified intact; histograms now sized "
              << counterArrays.events3D.size() << '\n';
}

// -----------------------------------------------------------------------------
// Test 5: The implementation uses push_back and never clears the vectors, so a
//         second call APPENDS another block of zeros rather than resetting.
//         This documents the actual (append) semantics of the code.
// -----------------------------------------------------------------------------
void test_iae_second_call_appends()
{
    std::cerr << "\n[TEST] test_iae_second_call_appends\n"
              << "  Source file:   src/reactions/init_association_events.cpp\n"
              << "  Function:      init_association_events(copyCounters&)\n"
              << "  Scenario:      the function is invoked twice on the same object.\n"
              << "  Pass criteria: because the code uses push_back without clearing,\n"
              << "                 the histograms double in length (2 * eventArraySize)\n"
              << "                 and all elements are still 0.\n";

    copyCounters counterArrays;
    counterArrays.eventArraySize = 5;

    std::cerr << "  First call...\n";
    init_association_events(counterArrays);
    EXPECT_EQ(counterArrays.events3D.size(), 5u) << "first call should add 5 bins";

    std::cerr << "  Second call (should append, not reset)...\n";
    init_association_events(counterArrays);

    std::cerr << "  Sizes after two calls: 3D=" << counterArrays.events3D.size()
              << ", 3Dto2D=" << counterArrays.events3Dto2D.size()
              << ", 2D=" << counterArrays.events2D.size() << " (expected 10 each)\n";

    EXPECT_EQ(counterArrays.events3D.size(), 10u)
        << "second call appends another eventArraySize entries to events3D";
    EXPECT_EQ(counterArrays.events3Dto2D.size(), 10u)
        << "second call appends another eventArraySize entries to events3Dto2D";
    EXPECT_EQ(counterArrays.events2D.size(), 10u)
        << "second call appends another eventArraySize entries to events2D";

    iae_expect_all_zero(counterArrays.events3D, "events3D");
    iae_expect_all_zero(counterArrays.events3Dto2D, "events3Dto2D");
    iae_expect_all_zero(counterArrays.events2D, "events2D");
}

// -----------------------------------------------------------------------------
// Test 6: A pre-existing (non-empty) histogram is preserved and the new zeroed
//         bins are appended after it -- again a direct consequence of push_back.
// -----------------------------------------------------------------------------
void test_iae_preserves_preexisting_entries()
{
    std::cerr << "\n[TEST] test_iae_preserves_preexisting_entries\n"
              << "  Source file:   src/reactions/init_association_events.cpp\n"
              << "  Function:      init_association_events(copyCounters&)\n"
              << "  Scenario:      events3D already holds a non-zero sentinel value.\n"
              << "  Pass criteria: the sentinel survives at index 0 and the new bins\n"
              << "                 (indices 1..eventArraySize) are all 0.\n";

    copyCounters counterArrays;
    counterArrays.eventArraySize = 4;

    // Put a recognizable sentinel in place before calling.
    counterArrays.events3D.push_back(99);

    std::cerr << "  Calling init_association_events with a pre-seeded events3D...\n";
    init_association_events(counterArrays);

    std::cerr << "  events3D.size()=" << counterArrays.events3D.size()
              << " (expected 1 sentinel + 4 new bins = 5)\n";

    ASSERT_EQ(counterArrays.events3D.size(), 5u)
        << "the 4 new bins should be appended after the single pre-existing entry";
    EXPECT_EQ(counterArrays.events3D[0], 99)
        << "the pre-existing sentinel value must be preserved at index 0";
    for (std::size_t i = 1; i < counterArrays.events3D.size(); ++i) {
        EXPECT_EQ(counterArrays.events3D[i], 0)
            << "newly appended bin events3D[" << i << "] should be 0";
    }

    // The other two histograms had no sentinel, so they get exactly 4 zeros.
    EXPECT_EQ(counterArrays.events3Dto2D.size(), 4u)
        << "events3Dto2D receives exactly eventArraySize new bins";
    EXPECT_EQ(counterArrays.events2D.size(), 4u)
        << "events2D receives exactly eventArraySize new bins";
    iae_expect_all_zero(counterArrays.events3Dto2D, "events3Dto2D");
    iae_expect_all_zero(counterArrays.events2D, "events2D");
}

// -----------------------------------------------------------------------------
// Test 7: Degenerate sizes. With eventArraySize == 0 (or negative) the `for`
//         loop condition fails immediately, so nothing is appended at all.
//         Neither case reaches any exit()/abort() path in the code.
// -----------------------------------------------------------------------------
void test_iae_zero_and_negative_array_size()
{
    std::cerr << "\n[TEST] test_iae_zero_and_negative_array_size\n"
              << "  Source file:   src/reactions/init_association_events.cpp\n"
              << "  Function:      init_association_events(copyCounters&)\n"
              << "  Scenario:      eventArraySize set to 0, then to a negative value.\n"
              << "  Pass criteria: the loop body never runs, so all three histograms\n"
              << "                 remain empty in both cases.\n";

    // --- eventArraySize == 0 -------------------------------------------------
    {
        copyCounters counterArrays;
        counterArrays.eventArraySize = 0;

        std::cerr << "  Calling init_association_events with eventArraySize = 0 ...\n";
        init_association_events(counterArrays);

        std::cerr << "    resulting sizes: 3D=" << counterArrays.events3D.size()
                  << ", 3Dto2D=" << counterArrays.events3Dto2D.size()
                  << ", 2D=" << counterArrays.events2D.size() << '\n';

        EXPECT_TRUE(counterArrays.events3D.empty())
            << "no bins should be created when eventArraySize is 0";
        EXPECT_TRUE(counterArrays.events3Dto2D.empty())
            << "no bins should be created when eventArraySize is 0";
        EXPECT_TRUE(counterArrays.events2D.empty())
            << "no bins should be created when eventArraySize is 0";
    }

    // --- eventArraySize < 0 --------------------------------------------------
    {
        copyCounters counterArrays;
        counterArrays.eventArraySize = -3;

        std::cerr << "  Calling init_association_events with eventArraySize = -3 ...\n";
        init_association_events(counterArrays);

        std::cerr << "    resulting sizes: 3D=" << counterArrays.events3D.size()
                  << ", 3Dto2D=" << counterArrays.events3Dto2D.size()
                  << ", 2D=" << counterArrays.events2D.size() << '\n';

        EXPECT_TRUE(counterArrays.events3D.empty())
            << "a negative eventArraySize must not create any bins";
        EXPECT_TRUE(counterArrays.events3Dto2D.empty())
            << "a negative eventArraySize must not create any bins";
        EXPECT_TRUE(counterArrays.events2D.empty())
            << "a negative eventArraySize must not create any bins";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper is invoked from its own TEST so
// the framework reports the results separately while still running every case
// (all assertions are non-fatal EXPECT_* except for the size preconditions that
// guard subsequent indexing).
// -----------------------------------------------------------------------------
TEST(InitAssociationEvents, DefaultArraySize) { test_iae_default_array_size(); }
TEST(InitAssociationEvents, CustomArraySize) { test_iae_custom_array_size(); }
TEST(InitAssociationEvents, AllThreeHistogramsMatch) { test_iae_all_three_histograms_match(); }
TEST(InitAssociationEvents, LeavesOtherCountersUntouched) { test_iae_leaves_other_counters_untouched(); }
TEST(InitAssociationEvents, SecondCallAppends) { test_iae_second_call_appends(); }
TEST(InitAssociationEvents, PreservesPreexistingEntries) { test_iae_preserves_preexisting_entries(); }
TEST(InitAssociationEvents, ZeroAndNegativeArraySize) { test_iae_zero_and_negative_array_size(); }