/*! \file test_update_Nboundpairs.cpp
 *
 * ### Unit test for ../src/reactions/update_Nboundpairs.cpp
 *
 * Function under test:
 *
 *     void update_Nboundpairs(int ptype1, int ptype2, int chg,
 *                             const Parameters& params,
 *                             copyCounters& counterArrays)
 *
 * Behaviour (read directly out of the implementation):
 *
 *   * `Npro = params.numMolTypes`.
 *   * If `ptype2 == -1` the event is a "bind to surface" event and the bin
 *     index used is `Npro*Npro + ptype1`, i.e. the surface bins live in a
 *     block that sits immediately after the Npro x Npro protein-pair block.
 *   * Otherwise the bin index is `ptype1*Npro + ptype2`, but if
 *     `ptype1 > ptype2` the indices are swapped so that the index becomes
 *     `ptype2*Npro + ptype1`.  This makes the storage *upper triangular*:
 *     the pair (a,b) and the pair (b,a) always land in the same bin.
 *   * The bin is updated with `+= chg`, so `chg` may be positive (bond
 *     formed) or negative (bond broken), and repeated calls accumulate.
 *
 * The function performs **no** bounds checking, so every test below sizes
 * `counterArrays.nBoundPairs` to `Npro*Npro + Npro` entries before calling.
 *
 * Verbose progress information is written to stderr so a reader of the test
 * log can follow exactly which behaviour is being probed and what the pass
 * criterion is.
 */

#include "classes/class_Parameters.hpp"
#include "classes/class_copyCounters.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Declaration of the function under test.
//
// The public declaration lives in reactions/shared_reaction_functions.hpp, but
// that header transitively pulls in a large amount of unrelated machinery
// (SimulVolume, GSL matrices, split.cpp ...).  Since we only need this one
// symbol, we declare it directly here with exactly the signature used by the
// implementation.
// -----------------------------------------------------------------------------
void update_Nboundpairs(int ptype1, int ptype2, int chg, const Parameters& params,
    copyCounters& counterArrays);

namespace {

/*! \brief Build a Parameters object that only carries the field the function
 *         under test actually reads (numMolTypes).
 *
 * \param[in] numMolTypes Number of distinct molecule types in the "system".
 * \return A fully-default Parameters with numMolTypes assigned.
 */
Parameters nbp_make_params(int numMolTypes)
{
    Parameters params; // every other member has an in-class initializer
    params.numMolTypes = numMolTypes;
    return params;
}

/*! \brief Build a copyCounters whose nBoundPairs vector is large enough for
 *         both the Npro x Npro protein-pair block and the Npro surface bins.
 *
 * update_Nboundpairs indexes nBoundPairs without any bounds check, so the
 * vector must be at least Npro*Npro + Npro entries long.
 *
 * \param[in] numMolTypes Number of distinct molecule types.
 * \return A copyCounters with an all-zero nBoundPairs array.
 */
copyCounters nbp_make_counters(int numMolTypes)
{
    copyCounters counterArrays;
    counterArrays.nBoundPairs.assign(
        static_cast<size_t>(numMolTypes) * numMolTypes + numMolTypes, 0);
    return counterArrays;
}

/*! \brief Print the whole nBoundPairs array to stderr for diagnostics. */
void nbp_dump_counters(const copyCounters& counterArrays, const char* label)
{
    std::cerr << "    " << label << " nBoundPairs = [";
    for (size_t i = 0; i < counterArrays.nBoundPairs.size(); ++i) {
        if (i != 0)
            std::cerr << ", ";
        std::cerr << counterArrays.nBoundPairs[i];
    }
    std::cerr << "]\n";
}

/*! \brief Count how many bins in nBoundPairs are non-zero. */
int nbp_count_nonzero(const copyCounters& counterArrays)
{
    int nonZero = 0;
    for (int val : counterArrays.nBoundPairs) {
        if (val != 0)
            ++nonZero;
    }
    return nonZero;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: A simple, ordered protein pair (ptype1 < ptype2) lands in the bin
//         ptype1*Npro + ptype2, and only that bin changes.
// -----------------------------------------------------------------------------
void nbp_test_ordered_pair_indexing()
{
    std::cerr << "\n[TEST] nbp_test_ordered_pair_indexing\n"
              << "  Source file:   src/reactions/update_Nboundpairs.cpp\n"
              << "  Function:      update_Nboundpairs()\n"
              << "  Scenario:      ptype1 (=0) < ptype2 (=2), Npro = 3, chg = +1.\n"
              << "  Pass criteria: bin index 0*3+2 = 2 becomes 1 and every other\n"
              << "                 bin remains exactly 0.\n";

    const int Npro = 3;
    Parameters params = nbp_make_params(Npro);
    copyCounters counterArrays = nbp_make_counters(Npro);

    std::cerr << "  Calling update_Nboundpairs(0, 2, +1, ...)\n";
    update_Nboundpairs(0, 2, 1, params, counterArrays);
    nbp_dump_counters(counterArrays, "after:");

    // ptype1 == 0, ptype2 == 2  ->  ind = 0*3 + 2 = 2
    EXPECT_EQ(counterArrays.nBoundPairs[2], 1)
        << "Pair (0,2) with Npro=3 must be stored in bin 0*3+2 = 2";

    // Nothing else may have been touched.
    EXPECT_EQ(nbp_count_nonzero(counterArrays), 1)
        << "Exactly one bin should be non-zero after a single update";
}

// -----------------------------------------------------------------------------
// Test 2: The mapping is symmetric.  Calling with (a,b) and with (b,a) must
//         update the *same* bin, because the implementation swaps the indices
//         whenever ptype1 > ptype2.
// -----------------------------------------------------------------------------
void nbp_test_symmetric_pair_mapping()
{
    std::cerr << "\n[TEST] nbp_test_symmetric_pair_mapping\n"
              << "  Source file:   src/reactions/update_Nboundpairs.cpp\n"
              << "  Function:      update_Nboundpairs()\n"
              << "  Scenario:      call once with (1,3) then once with (3,1), Npro = 4.\n"
              << "  Pass criteria: both calls accumulate into the single bin\n"
              << "                 min*Npro+max = 1*4+3 = 7 (value 2), and the\n"
              << "                 mirror bin 3*4+1 = 13 is never used (stays 0).\n";

    const int Npro = 4;
    Parameters params = nbp_make_params(Npro);
    copyCounters counterArrays = nbp_make_counters(Npro);

    std::cerr << "  Calling update_Nboundpairs(1, 3, +1, ...)\n";
    update_Nboundpairs(1, 3, 1, params, counterArrays);
    nbp_dump_counters(counterArrays, "after first call: ");

    // Snapshot after the "ordered" call so we can compare against the swapped one.
    const int afterOrdered = counterArrays.nBoundPairs[1 * Npro + 3];
    EXPECT_EQ(afterOrdered, 1) << "Ordered call (1,3) must fill bin 1*4+3 = 7";

    std::cerr << "  Calling update_Nboundpairs(3, 1, +1, ...) (indices reversed)\n";
    update_Nboundpairs(3, 1, 1, params, counterArrays);
    nbp_dump_counters(counterArrays, "after second call:");

    // The reversed call must accumulate into the very same bin.
    EXPECT_EQ(counterArrays.nBoundPairs[1 * Npro + 3], 2)
        << "Reversed call (3,1) must accumulate into the same upper-triangular bin 7";

    // The "lower triangular" mirror slot must never be written to.
    EXPECT_EQ(counterArrays.nBoundPairs[3 * Npro + 1], 0)
        << "Bin 3*4+1 = 13 is the mirror slot and must remain untouched";

    // Only the single shared bin should be non-zero.
    EXPECT_EQ(nbp_count_nonzero(counterArrays), 1)
        << "Both calls must land in exactly one shared bin";
}

// -----------------------------------------------------------------------------
// Test 3: A self-pair (ptype1 == ptype2) has no swap applied and lands on the
//         diagonal bin i*Npro + i.
// -----------------------------------------------------------------------------
void nbp_test_self_pair_diagonal()
{
    std::cerr << "\n[TEST] nbp_test_self_pair_diagonal\n"
              << "  Source file:   src/reactions/update_Nboundpairs.cpp\n"
              << "  Function:      update_Nboundpairs()\n"
              << "  Scenario:      homodimer event with ptype1 == ptype2 == 2, Npro = 3.\n"
              << "  Pass criteria: diagonal bin 2*3+2 = 8 receives the change; no\n"
              << "                 index swap occurs because the swap only fires\n"
              << "                 when ptype1 > ptype2 (strictly).\n";

    const int Npro = 3;
    Parameters params = nbp_make_params(Npro);
    copyCounters counterArrays = nbp_make_counters(Npro);

    std::cerr << "  Calling update_Nboundpairs(2, 2, +5, ...)\n";
    update_Nboundpairs(2, 2, 5, params, counterArrays);
    nbp_dump_counters(counterArrays, "after:");

    EXPECT_EQ(counterArrays.nBoundPairs[2 * Npro + 2], 5)
        << "Homodimer pair (2,2) must be stored in the diagonal bin 2*3+2 = 8";
    EXPECT_EQ(nbp_count_nonzero(counterArrays), 1)
        << "Only the diagonal bin should be non-zero";
}

// -----------------------------------------------------------------------------
// Test 4: Surface binding.  When ptype2 == -1 the implementation switches to the
//         surface block that starts at offset Npro*Npro.
// -----------------------------------------------------------------------------
void nbp_test_surface_binding_index()
{
    std::cerr << "\n[TEST] nbp_test_surface_binding_index\n"
              << "  Source file:   src/reactions/update_Nboundpairs.cpp\n"
              << "  Function:      update_Nboundpairs()\n"
              << "  Scenario:      ptype2 == -1 flags a bind-to-surface event.\n"
              << "                 Npro = 3, ptype1 sweeps 0, 1 and 2.\n"
              << "  Pass criteria: each ptype1 updates bin Npro*Npro + ptype1,\n"
              << "                 i.e. bins 9, 10 and 11, and nothing in the\n"
              << "                 protein-pair block (bins 0..8) changes.\n";

    const int Npro = 3;
    Parameters params = nbp_make_params(Npro);
    copyCounters counterArrays = nbp_make_counters(Npro);

    // Give each surface bin a distinguishable value so we can tell them apart.
    std::cerr << "  Calling update_Nboundpairs(0, -1, +1, ...)\n";
    update_Nboundpairs(0, -1, 1, params, counterArrays);
    std::cerr << "  Calling update_Nboundpairs(1, -1, +2, ...)\n";
    update_Nboundpairs(1, -1, 2, params, counterArrays);
    std::cerr << "  Calling update_Nboundpairs(2, -1, +3, ...)\n";
    update_Nboundpairs(2, -1, 3, params, counterArrays);
    nbp_dump_counters(counterArrays, "after:");

    EXPECT_EQ(counterArrays.nBoundPairs[Npro * Npro + 0], 1)
        << "Surface bin for molecule type 0 is Npro*Npro+0 = 9";
    EXPECT_EQ(counterArrays.nBoundPairs[Npro * Npro + 1], 2)
        << "Surface bin for molecule type 1 is Npro*Npro+1 = 10";
    EXPECT_EQ(counterArrays.nBoundPairs[Npro * Npro + 2], 3)
        << "Surface bin for molecule type 2 is Npro*Npro+2 = 11";

    // The whole protein-protein block must be untouched by surface events.
    for (int i = 0; i < Npro * Npro; ++i) {
        EXPECT_EQ(counterArrays.nBoundPairs[i], 0)
            << "Protein-pair bin " << i << " must not be modified by a surface event";
    }
}

// -----------------------------------------------------------------------------
// Test 5: Accumulation semantics.  The function uses `+=`, so positive and
//         negative changes accumulate and can be driven back to zero.
// -----------------------------------------------------------------------------
void nbp_test_accumulation_and_decrement()
{
    std::cerr << "\n[TEST] nbp_test_accumulation_and_decrement\n"
              << "  Source file:   src/reactions/update_Nboundpairs.cpp\n"
              << "  Function:      update_Nboundpairs()\n"
              << "  Scenario:      repeated association (+1) events on pair (0,1)\n"
              << "                 followed by dissociation (-1) events, Npro = 2.\n"
              << "  Pass criteria: the bin tracks the running sum exactly, and a\n"
              << "                 matching number of +1/-1 calls returns it to 0.\n";

    const int Npro = 2;
    Parameters params = nbp_make_params(Npro);
    copyCounters counterArrays = nbp_make_counters(Npro);

    const int bin = 0 * Npro + 1; // pair (0,1) -> bin 1

    // Three association events.
    for (int i = 0; i < 3; ++i) {
        std::cerr << "  Association event " << (i + 1)
                  << ": update_Nboundpairs(0, 1, +1, ...)\n";
        update_Nboundpairs(0, 1, 1, params, counterArrays);
        EXPECT_EQ(counterArrays.nBoundPairs[bin], i + 1)
            << "Bin " << bin << " should hold the running count of bound pairs";
    }

    // Two dissociation events (note the reversed argument order to also confirm
    // the symmetric mapping applies to negative changes too).
    std::cerr << "  Dissociation event 1: update_Nboundpairs(1, 0, -1, ...)\n";
    update_Nboundpairs(1, 0, -1, params, counterArrays);
    EXPECT_EQ(counterArrays.nBoundPairs[bin], 2)
        << "A -1 change delivered with reversed indices must decrement the same bin";

    std::cerr << "  Dissociation event 2: update_Nboundpairs(1, 0, -2, ...)\n";
    update_Nboundpairs(1, 0, -2, params, counterArrays);
    nbp_dump_counters(counterArrays, "after:");

    EXPECT_EQ(counterArrays.nBoundPairs[bin], 0)
        << "Equal numbers of formed and broken bonds must return the bin to zero";
    EXPECT_EQ(nbp_count_nonzero(counterArrays), 0)
        << "No other bin should have been disturbed during accumulation";
}

// -----------------------------------------------------------------------------
// Test 6: A change of zero is a genuine no-op, and pre-existing bin contents
//         are preserved rather than overwritten.
// -----------------------------------------------------------------------------
void nbp_test_zero_change_is_noop()
{
    std::cerr << "\n[TEST] nbp_test_zero_change_is_noop\n"
              << "  Source file:   src/reactions/update_Nboundpairs.cpp\n"
              << "  Function:      update_Nboundpairs()\n"
              << "  Scenario:      pre-seed every bin with a known value, then call\n"
              << "                 with chg = 0, Npro = 2.\n"
              << "  Pass criteria: the whole nBoundPairs array is bit-for-bit\n"
              << "                 unchanged (the function adds, it does not assign).\n";

    const int Npro = 2;
    Parameters params = nbp_make_params(Npro);
    copyCounters counterArrays = nbp_make_counters(Npro);

    // Seed the array with recognisable, non-zero values.
    for (size_t i = 0; i < counterArrays.nBoundPairs.size(); ++i)
        counterArrays.nBoundPairs[i] = static_cast<int>(i) + 100;

    const std::vector<int> before = counterArrays.nBoundPairs;
    nbp_dump_counters(counterArrays, "before:");

    std::cerr << "  Calling update_Nboundpairs(0, 1, 0, ...) (chg == 0)\n";
    update_Nboundpairs(0, 1, 0, params, counterArrays);
    std::cerr << "  Calling update_Nboundpairs(1, -1, 0, ...) (surface, chg == 0)\n";
    update_Nboundpairs(1, -1, 0, params, counterArrays);
    nbp_dump_counters(counterArrays, "after: ");

    ASSERT_EQ(counterArrays.nBoundPairs.size(), before.size())
        << "update_Nboundpairs must not resize the nBoundPairs vector";
    for (size_t i = 0; i < before.size(); ++i) {
        EXPECT_EQ(counterArrays.nBoundPairs[i], before[i])
            << "Bin " << i << " must be unchanged when chg == 0";
    }
}

// -----------------------------------------------------------------------------
// Test 7: Distinct pairs must not collide.  Fill in every upper-triangular pair
//         for a small Npro and confirm each one occupies its own bin.
// -----------------------------------------------------------------------------
void nbp_test_distinct_pairs_do_not_collide()
{
    std::cerr << "\n[TEST] nbp_test_distinct_pairs_do_not_collide\n"
              << "  Source file:   src/reactions/update_Nboundpairs.cpp\n"
              << "  Function:      update_Nboundpairs()\n"
              << "  Scenario:      Npro = 3; register one bond for every unique\n"
              << "                 unordered pair (0,0) (0,1) (0,2) (1,1) (1,2) (2,2).\n"
              << "  Pass criteria: six distinct upper-triangular bins each hold 1,\n"
              << "                 the three lower-triangular mirrors hold 0, and\n"
              << "                 the three surface bins hold 0.\n";

    const int Npro = 3;
    Parameters params = nbp_make_params(Npro);
    copyCounters counterArrays = nbp_make_counters(Npro);

    // Register one bond for every unique unordered pair.
    for (int a = 0; a < Npro; ++a) {
        for (int b = a; b < Npro; ++b) {
            std::cerr << "  Calling update_Nboundpairs(" << a << ", " << b << ", +1, ...)"
                      << " -> expected bin " << (a * Npro + b) << '\n';
            update_Nboundpairs(a, b, 1, params, counterArrays);
        }
    }
    nbp_dump_counters(counterArrays, "after:");

    // Every upper-triangular bin (including the diagonal) should hold exactly 1.
    for (int a = 0; a < Npro; ++a) {
        for (int b = a; b < Npro; ++b) {
            EXPECT_EQ(counterArrays.nBoundPairs[a * Npro + b], 1)
                << "Unique pair (" << a << ',' << b << ") should occupy its own bin "
                << (a * Npro + b);
        }
    }

    // Every strictly lower-triangular bin must remain unused.
    for (int a = 1; a < Npro; ++a) {
        for (int b = 0; b < a; ++b) {
            EXPECT_EQ(counterArrays.nBoundPairs[a * Npro + b], 0)
                << "Lower-triangular mirror bin " << (a * Npro + b)
                << " must never be written to";
        }
    }

    // Surface bins are untouched by protein-protein events.
    for (int i = 0; i < Npro; ++i) {
        EXPECT_EQ(counterArrays.nBoundPairs[Npro * Npro + i], 0)
            << "Surface bin " << (Npro * Npro + i)
            << " must not be modified by protein-pair events";
    }

    // Total bookkeeping check: 6 unique pairs for Npro == 3.
    EXPECT_EQ(nbp_count_nonzero(counterArrays), 6)
        << "Npro=3 has 3*(3+1)/2 = 6 unique unordered pairs";
}

// -----------------------------------------------------------------------------
// Test 8: The bin layout depends on params.numMolTypes.  The same (ptype1,
//         ptype2) arguments must map to different indices for different Npro.
// -----------------------------------------------------------------------------
void nbp_test_index_depends_on_numMolTypes()
{
    std::cerr << "\n[TEST] nbp_test_index_depends_on_numMolTypes\n"
              << "  Source file:   src/reactions/update_Nboundpairs.cpp\n"
              << "  Function:      update_Nboundpairs()\n"
              << "  Scenario:      the same pair (1,2) is registered twice, once in\n"
              << "                 a system with Npro = 3 and once with Npro = 5.\n"
              << "  Pass criteria: the bin used is 1*3+2 = 5 in the first case and\n"
              << "                 1*5+2 = 7 in the second, showing the stride is\n"
              << "                 taken from params.numMolTypes.\n";

    // --- Npro = 3 -----------------------------------------------------------
    {
        const int Npro = 3;
        Parameters params = nbp_make_params(Npro);
        copyCounters counterArrays = nbp_make_counters(Npro);

        std::cerr << "  [Npro=3] Calling update_Nboundpairs(1, 2, +1, ...)\n";
        update_Nboundpairs(1, 2, 1, params, counterArrays);
        nbp_dump_counters(counterArrays, "[Npro=3] after:");

        EXPECT_EQ(counterArrays.nBoundPairs[1 * 3 + 2], 1)
            << "With Npro=3 the pair (1,2) belongs in bin 5";
        EXPECT_EQ(nbp_count_nonzero(counterArrays), 1)
            << "Only one bin should be non-zero for Npro=3";
    }

    // --- Npro = 5 -----------------------------------------------------------
    {
        const int Npro = 5;
        Parameters params = nbp_make_params(Npro);
        copyCounters counterArrays = nbp_make_counters(Npro);

        std::cerr << "  [Npro=5] Calling update_Nboundpairs(1, 2, +1, ...)\n";
        update_Nboundpairs(1, 2, 1, params, counterArrays);
        nbp_dump_counters(counterArrays, "[Npro=5] after:");

        EXPECT_EQ(counterArrays.nBoundPairs[1 * 5 + 2], 1)
            << "With Npro=5 the same pair (1,2) belongs in bin 7";
        EXPECT_EQ(counterArrays.nBoundPairs[1 * 3 + 2], 0)
            << "Bin 5 (the Npro=3 answer) must be empty when Npro=5";
        EXPECT_EQ(nbp_count_nonzero(counterArrays), 1)
            << "Only one bin should be non-zero for Npro=5";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* helper is invoked from its own TEST
// so that a failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(UpdateNboundPairs, OrderedPairIndexing) { nbp_test_ordered_pair_indexing(); }
TEST(UpdateNboundPairs, SymmetricPairMapping) { nbp_test_symmetric_pair_mapping(); }
TEST(UpdateNboundPairs, SelfPairDiagonal) { nbp_test_self_pair_diagonal(); }
TEST(UpdateNboundPairs, SurfaceBindingIndex) { nbp_test_surface_binding_index(); }
TEST(UpdateNboundPairs, AccumulationAndDecrement) { nbp_test_accumulation_and_decrement(); }
TEST(UpdateNboundPairs, ZeroChangeIsNoop) { nbp_test_zero_change_is_noop(); }
TEST(UpdateNboundPairs, DistinctPairsDoNotCollide) { nbp_test_distinct_pairs_do_not_collide(); }
TEST(UpdateNboundPairs, IndexDependsOnNumMolTypes) { nbp_test_index_depends_on_numMolTypes(); }