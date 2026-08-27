/*! \file test_measure_overlap_protein_interfaces.cpp
 *
 * ### Unit test for src/reactions/measure_overlap_protein_interfaces.cpp
 *
 * Function under test:
 *
 *     void measure_overlap_protein_interfaces(Molecule base1, Molecule baseTmp,
 *                                             bool& flagCancel)
 *
 * Semantics taken directly from the implementation:
 *
 *   - `base1` represents a molecule that already sits in the system, so its
 *     *real* coordinates (`base1.interfaceList[i].coord`) are used.
 *   - `baseTmp` represents the molecule that is currently associating, so its
 *     *temporary* association coordinates (`baseTmp.tmpICoords[m]`) are used.
 *   - The comparison is done over every interface pair (i, m).
 *   - The hard-coded squared binding radius is `bindrad2 = 1.0`, so the test
 *     is a strict `d2 < 1.0` (i.e. separation strictly less than 1.0 nm).
 *   - `flagCancel` is only ever set to `true`; it is never cleared. A caller
 *     that passes in `true` gets `true` back regardless of geometry.
 *   - NOTE (implementation detail that the tests pin down): the inner loop is
 *     bounded by `baseTmp.interfaceList.size()` but indexes into
 *     `baseTmp.tmpICoords`. The two containers must therefore be the same size
 *     for the call to be safe. Every molecule built below keeps them in sync,
 *     except for one deliberate case where `interfaceList` is *empty* (which
 *     simply means the inner loop never runs -- still safe, no out-of-bounds).
 *   - Both molecules are taken **by value**, so the caller's objects can never
 *     be modified by this routine.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "reactions/association/association.hpp"

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace + unique "mopi_" prefix so that nothing
// collides with the other translation units linked into the test binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully-initialized Molecule whose *real* interface
 *         coordinates are the ones supplied.
 *
 * The temporary association coordinates are also filled in (mirroring the real
 * coordinates) so that `interfaceList.size() == tmpICoords.size()`, which is
 * what the function under test requires for a safe inner loop.
 *
 * \param[in] ifaceCoords real coordinates for each interface
 * \return a Molecule ready to be handed to measure_overlap_protein_interfaces
 */
Molecule mopi_make_molecule(const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.index = 0;
    mol.myComIndex = 0;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.comCoord = Coord { 0.0, 0.0, 0.0 };

    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface;
        iface.coord = ifaceCoords[i];
        iface.relIndex = static_cast<int>(i);
        iface.index = static_cast<int>(i);
        iface.molTypeIndex = 0;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
    }

    // Mirror the real coordinates into the temporary association coordinates,
    // which is exactly what Molecule::set_tmp_association_coords() does.
    mol.set_tmp_association_coords();

    return mol;
}

/*! \brief Overwrite the temporary association coordinates of a molecule.
 *
 * The real interface coordinates are deliberately left alone so the tests can
 * prove that only `tmpICoords` is consulted for the associating partner.
 *
 * \param[in,out] mol       molecule to modify
 * \param[in]     tmpCoords new temporary coordinates (must match interface count)
 */
void mopi_set_tmp_coords(Molecule& mol, const std::vector<Coord>& tmpCoords)
{
    mol.tmpICoords.clear();
    for (const auto& crd : tmpCoords)
        mol.tmpICoords.push_back(crd);
}

/*! \brief Convenience: Euclidean distance between two coordinates. */
double mopi_distance(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: well-separated interfaces must NOT cancel the association.
// -----------------------------------------------------------------------------
void test_mopi_no_overlap_leaves_flag_false()
{
    std::cerr << "\n[TEST] test_mopi_no_overlap_leaves_flag_false\n"
              << "  Source file: src/reactions/measure_overlap_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_protein_interfaces()\n"
              << "  Scenario:    one interface each, separated by 5.0 nm which is far\n"
              << "               beyond the hard-coded binding radius of 1.0 nm.\n"
              << "  Pass rule:   flagCancel must still be false on return.\n";

    // System molecule with a single interface at the origin.
    Molecule base1 = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });

    // Associating molecule: real coords irrelevant, tmpICoords put 5 nm away.
    Molecule baseTmp = mopi_make_molecule({ Coord { 5.0, 0.0, 0.0 } });
    mopi_set_tmp_coords(baseTmp, { Coord { 5.0, 0.0, 0.0 } });

    std::cerr << "  Interface separation = "
              << mopi_distance(base1.interfaceList[0].coord, baseTmp.tmpICoords[0])
              << " nm (cutoff 1.0 nm)\n";

    bool flagCancel { false };
    measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

    EXPECT_FALSE(flagCancel)
        << "Interfaces 5 nm apart are outside the 1 nm cutoff; association must not be cancelled";

    std::cerr << "  Result: flagCancel = " << std::boolalpha << flagCancel << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: overlapping interfaces must cancel the association.
// -----------------------------------------------------------------------------
void test_mopi_overlap_sets_flag_true()
{
    std::cerr << "\n[TEST] test_mopi_overlap_sets_flag_true\n"
              << "  Source file: src/reactions/measure_overlap_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_protein_interfaces()\n"
              << "  Scenario:    one interface each, separated by only 0.5 nm which is\n"
              << "               inside the hard-coded 1.0 nm binding radius.\n"
              << "  Pass rule:   flagCancel must be flipped to true.\n";

    Molecule base1 = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });

    Molecule baseTmp = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });
    mopi_set_tmp_coords(baseTmp, { Coord { 0.5, 0.0, 0.0 } });

    std::cerr << "  Interface separation = "
              << mopi_distance(base1.interfaceList[0].coord, baseTmp.tmpICoords[0])
              << " nm (cutoff 1.0 nm)\n";

    bool flagCancel { false };
    measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

    EXPECT_TRUE(flagCancel)
        << "Interfaces 0.5 nm apart overlap within the 1 nm cutoff; association must be cancelled";

    std::cerr << "  Result: flagCancel = " << std::boolalpha << flagCancel << '\n';
}

// -----------------------------------------------------------------------------
// Test 3: the comparison is strict (`d2 < bindrad2`), so a separation of
//         exactly 1.0 nm does NOT cancel, while just under 1.0 nm does.
// -----------------------------------------------------------------------------
void test_mopi_cutoff_is_strictly_less_than()
{
    std::cerr << "\n[TEST] test_mopi_cutoff_is_strictly_less_than\n"
              << "  Source file: src/reactions/measure_overlap_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_protein_interfaces()\n"
              << "  Scenario:    probe the exact cutoff. The code compares d2 < 1.0, so\n"
              << "               a separation of exactly 1.0 nm gives d2 == 1.0 and must\n"
              << "               NOT cancel, whereas 0.999999 nm must cancel.\n"
              << "  Pass rule:   false at exactly 1.0 nm, true just inside it.\n";

    Molecule base1 = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });

    // --- Exactly on the cutoff: dx = 1.0 gives d2 = 1.0 exactly in IEEE-754. ---
    {
        Molecule baseTmp = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });
        mopi_set_tmp_coords(baseTmp, { Coord { 1.0, 0.0, 0.0 } });

        std::cerr << "  Sub-case A: separation exactly 1.0 nm -> d2 = 1.0\n";

        bool flagCancel { false };
        measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

        EXPECT_FALSE(flagCancel)
            << "d2 == bindrad2 (1.0) fails the strict `<` test, so it must not cancel";
        std::cerr << "            flagCancel = " << std::boolalpha << flagCancel << '\n';
    }

    // --- Just inside the cutoff. ---
    {
        Molecule baseTmp = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });
        mopi_set_tmp_coords(baseTmp, { Coord { 0.999999, 0.0, 0.0 } });

        std::cerr << "  Sub-case B: separation 0.999999 nm -> d2 slightly below 1.0\n";

        bool flagCancel { false };
        measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

        EXPECT_TRUE(flagCancel)
            << "d2 just below bindrad2 satisfies the strict `<` test, so it must cancel";
        std::cerr << "            flagCancel = " << std::boolalpha << flagCancel << '\n';
    }
}

// -----------------------------------------------------------------------------
// Test 4: the flag is write-only-true; a pre-set true flag survives even when
//         nothing overlaps (the routine never resets it).
// -----------------------------------------------------------------------------
void test_mopi_flag_is_never_cleared()
{
    std::cerr << "\n[TEST] test_mopi_flag_is_never_cleared\n"
              << "  Source file: src/reactions/measure_overlap_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_protein_interfaces()\n"
              << "  Scenario:    caller already decided to cancel (flagCancel = true) and\n"
              << "               the geometry is far from overlapping.\n"
              << "  Pass rule:   flagCancel must remain true; the routine only ever sets\n"
              << "               it to true and never writes false.\n";

    Molecule base1 = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });
    Molecule baseTmp = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });
    mopi_set_tmp_coords(baseTmp, { Coord { 100.0, 100.0, 100.0 } });

    bool flagCancel { true };
    measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

    EXPECT_TRUE(flagCancel)
        << "An already-true cancellation flag must be preserved (the routine never clears it)";

    std::cerr << "  Result: flagCancel = " << std::boolalpha << flagCancel
              << " (input was true, geometry was non-overlapping)\n";
}

// -----------------------------------------------------------------------------
// Test 5: with several interfaces on each molecule a single overlapping pair is
//         enough to cancel; a configuration where every pair is far apart is not.
// -----------------------------------------------------------------------------
void test_mopi_multiple_interfaces_pairwise_scan()
{
    std::cerr << "\n[TEST] test_mopi_multiple_interfaces_pairwise_scan\n"
              << "  Source file: src/reactions/measure_overlap_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_protein_interfaces()\n"
              << "  Scenario:    3 interfaces on base1 x 3 interfaces on baseTmp. The\n"
              << "               routine scans every one of the 9 pairs.\n"
              << "  Pass rule:   false when all 9 pairs are far apart; true as soon as a\n"
              << "               single pair (here base1[2] vs baseTmp[0]) overlaps.\n";

    // base1 has three widely spread interfaces.
    Molecule base1 = mopi_make_molecule({
        Coord { 0.0, 0.0, 0.0 },
        Coord { 10.0, 0.0, 0.0 },
        Coord { 20.0, 0.0, 0.0 },
    });

    // --- Sub-case A: nothing overlaps. ---
    {
        Molecule baseTmp = mopi_make_molecule({
            Coord { 0.0, 0.0, 0.0 },
            Coord { 0.0, 0.0, 0.0 },
            Coord { 0.0, 0.0, 0.0 },
        });
        mopi_set_tmp_coords(baseTmp, {
                                         Coord { 0.0, 50.0, 0.0 },
                                         Coord { 10.0, 50.0, 0.0 },
                                         Coord { 20.0, 50.0, 0.0 },
                                     });

        std::cerr << "  Sub-case A: all 9 interface pairs are >= 50 nm apart\n";

        bool flagCancel { false };
        measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

        EXPECT_FALSE(flagCancel)
            << "No pair among the 9 combinations is within 1 nm, so nothing should cancel";
        std::cerr << "            flagCancel = " << std::boolalpha << flagCancel << '\n';
    }

    // --- Sub-case B: exactly one pair overlaps (base1[2] vs baseTmp[0]). ---
    {
        Molecule baseTmp = mopi_make_molecule({
            Coord { 0.0, 0.0, 0.0 },
            Coord { 0.0, 0.0, 0.0 },
            Coord { 0.0, 0.0, 0.0 },
        });
        mopi_set_tmp_coords(baseTmp, {
                                         Coord { 20.3, 0.0, 0.0 }, // 0.3 nm from base1[2] -> overlap
                                         Coord { 10.0, 50.0, 0.0 },
                                         Coord { 20.0, 50.0, 0.0 },
                                     });

        std::cerr << "  Sub-case B: base1[2]=(20,0,0) vs baseTmp.tmp[0]=(20.3,0,0), "
                  << "separation "
                  << mopi_distance(base1.interfaceList[2].coord, baseTmp.tmpICoords[0])
                  << " nm\n";

        bool flagCancel { false };
        measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

        EXPECT_TRUE(flagCancel)
            << "A single overlapping interface pair out of nine must cancel the association";
        std::cerr << "            flagCancel = " << std::boolalpha << flagCancel << '\n';
    }
}

// -----------------------------------------------------------------------------
// Test 6: the routine reads base1's REAL coordinates and baseTmp's TEMPORARY
//         coordinates -- and nothing else.
// -----------------------------------------------------------------------------
void test_mopi_uses_real_coords_for_base1_and_tmp_for_basetmp()
{
    std::cerr << "\n[TEST] test_mopi_uses_real_coords_for_base1_and_tmp_for_basetmp\n"
              << "  Source file: src/reactions/measure_overlap_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_protein_interfaces()\n"
              << "  Scenario:    craft molecules where the 'wrong' coordinate arrays would\n"
              << "               give the opposite answer, proving which arrays are read.\n"
              << "  Pass rule:   only base1.interfaceList[].coord and baseTmp.tmpICoords[]\n"
              << "               influence the outcome.\n";

    // --- Sub-case A: baseTmp's real coords overlap but its tmp coords do not. ---
    {
        Molecule base1 = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });

        // Real interface coordinate sits right on top of base1 ...
        Molecule baseTmp = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });
        // ... but the temporary (association) coordinate is far away.
        mopi_set_tmp_coords(baseTmp, { Coord { 30.0, 0.0, 0.0 } });

        std::cerr << "  Sub-case A: baseTmp real coord overlaps, tmp coord is 30 nm away\n";

        bool flagCancel { false };
        measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

        EXPECT_FALSE(flagCancel)
            << "baseTmp.interfaceList coords must be ignored; only tmpICoords are compared";
        std::cerr << "            flagCancel = " << std::boolalpha << flagCancel << '\n';
    }

    // --- Sub-case B: base1's tmp coords overlap but its real coords do not. ---
    {
        // base1 real interface is far from the associating partner ...
        Molecule base1 = mopi_make_molecule({ Coord { 30.0, 0.0, 0.0 } });
        // ... while its (unused) temporary coordinate sits on the partner.
        mopi_set_tmp_coords(base1, { Coord { 0.0, 0.0, 0.0 } });

        Molecule baseTmp = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });
        mopi_set_tmp_coords(baseTmp, { Coord { 0.0, 0.0, 0.0 } });

        std::cerr << "  Sub-case B: base1 tmp coord overlaps, real coord is 30 nm away\n";

        bool flagCancel { false };
        measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

        EXPECT_FALSE(flagCancel)
            << "base1.tmpICoords must be ignored; only its real interfaceList coords are compared";
        std::cerr << "            flagCancel = " << std::boolalpha << flagCancel << '\n';
    }
}

// -----------------------------------------------------------------------------
// Test 7: degenerate interface lists. The loop bounds come from
//         base1.interfaceList.size() and baseTmp.interfaceList.size(), so an
//         empty list on either side means zero comparisons are performed.
// -----------------------------------------------------------------------------
void test_mopi_empty_interface_lists_do_nothing()
{
    std::cerr << "\n[TEST] test_mopi_empty_interface_lists_do_nothing\n"
              << "  Source file: src/reactions/measure_overlap_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_protein_interfaces()\n"
              << "  Scenario:    (a) base1 has no interfaces, (b) baseTmp has no\n"
              << "               interfaces even though tmpICoords is populated and would\n"
              << "               overlap. The inner loop bound is interfaceList.size(),\n"
              << "               so no comparison happens in either case.\n"
              << "  Pass rule:   flagCancel stays false in both cases.\n";

    // --- Sub-case A: base1 has no interfaces at all -> outer loop never runs. ---
    {
        Molecule base1 = mopi_make_molecule({}); // empty interfaceList / tmpICoords
        Molecule baseTmp = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });
        mopi_set_tmp_coords(baseTmp, { Coord { 0.0, 0.0, 0.0 } });

        std::cerr << "  Sub-case A: base1 has 0 interfaces (outer loop body never executes)\n";

        bool flagCancel { false };
        measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

        EXPECT_FALSE(flagCancel)
            << "With no interfaces on base1 there is nothing to compare, so no cancellation";
        std::cerr << "            flagCancel = " << std::boolalpha << flagCancel << '\n';
    }

    // --- Sub-case B: baseTmp.interfaceList empty -> inner loop never runs even
    //     though tmpICoords holds a perfectly overlapping coordinate. This also
    //     documents that the inner loop bound is interfaceList.size(), NOT
    //     tmpICoords.size(). ---
    {
        Molecule base1 = mopi_make_molecule({ Coord { 0.0, 0.0, 0.0 } });

        Molecule baseTmp = mopi_make_molecule({}); // no interfaces
        // Deliberately give it a temporary coordinate right on top of base1.
        mopi_set_tmp_coords(baseTmp, { Coord { 0.0, 0.0, 0.0 } });

        std::cerr << "  Sub-case B: baseTmp has 0 interfaces but 1 overlapping tmpICoord\n"
                  << "              (inner loop is bounded by interfaceList.size(), so it\n"
                  << "               never executes and the overlap is never seen)\n";

        bool flagCancel { false };
        measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

        EXPECT_FALSE(flagCancel)
            << "Inner loop is bounded by baseTmp.interfaceList.size(); an empty list means "
               "the overlapping tmpICoord is never examined";
        std::cerr << "            flagCancel = " << std::boolalpha << flagCancel << '\n';
    }
}

// -----------------------------------------------------------------------------
// Test 8: both molecules are taken by value, so the caller's copies must be
//         completely untouched by the call.
// -----------------------------------------------------------------------------
void test_mopi_arguments_are_passed_by_value()
{
    std::cerr << "\n[TEST] test_mopi_arguments_are_passed_by_value\n"
              << "  Source file: src/reactions/measure_overlap_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_protein_interfaces()\n"
              << "  Scenario:    call with overlapping geometry (which triggers the cancel\n"
              << "               branch) and then re-inspect the caller's molecules.\n"
              << "  Pass rule:   coordinates and container sizes are unchanged, because the\n"
              << "               Molecules are parameters passed by value.\n";

    Molecule base1 = mopi_make_molecule({ Coord { 1.0, 2.0, 3.0 } });
    Molecule baseTmp = mopi_make_molecule({ Coord { 4.0, 5.0, 6.0 } });
    mopi_set_tmp_coords(baseTmp, { Coord { 1.2, 2.0, 3.0 } }); // 0.2 nm -> overlap

    // Snapshot the state we expect to survive the call.
    const Coord base1Before = base1.interfaceList[0].coord;
    const Coord baseTmpRealBefore = baseTmp.interfaceList[0].coord;
    const Coord baseTmpTmpBefore = baseTmp.tmpICoords[0];
    const std::size_t base1SizeBefore = base1.interfaceList.size();
    const std::size_t baseTmpSizeBefore = baseTmp.tmpICoords.size();

    bool flagCancel { false };
    measure_overlap_protein_interfaces(base1, baseTmp, flagCancel);

    // The geometry does overlap, so the flag should be set...
    EXPECT_TRUE(flagCancel)
        << "Sanity check: the chosen geometry (0.2 nm separation) must trigger cancellation";

    // ...but neither molecule may have been modified.
    EXPECT_TRUE(base1.interfaceList[0].coord == base1Before)
        << "base1's real interface coordinate must be unchanged (passed by value)";
    EXPECT_TRUE(baseTmp.interfaceList[0].coord == baseTmpRealBefore)
        << "baseTmp's real interface coordinate must be unchanged (passed by value)";
    EXPECT_TRUE(baseTmp.tmpICoords[0] == baseTmpTmpBefore)
        << "baseTmp's temporary coordinate must be unchanged (passed by value)";
    EXPECT_EQ(base1.interfaceList.size(), base1SizeBefore)
        << "base1's interface count must be unchanged";
    EXPECT_EQ(baseTmp.tmpICoords.size(), baseTmpSizeBefore)
        << "baseTmp's temporary coordinate count must be unchanged";

    std::cerr << "  Caller state after the call: base1 iface = " << base1.interfaceList[0].coord
              << ", baseTmp tmp iface = " << baseTmp.tmpICoords[0] << '\n';
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper is invoked from its own TEST so
// that a failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(MeasureOverlapProteinInterfaces, NoOverlapLeavesFlagFalse)
{
    test_mopi_no_overlap_leaves_flag_false();
}

TEST(MeasureOverlapProteinInterfaces, OverlapSetsFlagTrue)
{
    test_mopi_overlap_sets_flag_true();
}

TEST(MeasureOverlapProteinInterfaces, CutoffIsStrictlyLessThan)
{
    test_mopi_cutoff_is_strictly_less_than();
}

TEST(MeasureOverlapProteinInterfaces, FlagIsNeverCleared)
{
    test_mopi_flag_is_never_cleared();
}

TEST(MeasureOverlapProteinInterfaces, MultipleInterfacesPairwiseScan)
{
    test_mopi_multiple_interfaces_pairwise_scan();
}

TEST(MeasureOverlapProteinInterfaces, UsesRealCoordsForBase1AndTmpForBaseTmp)
{
    test_mopi_uses_real_coords_for_base1_and_tmp_for_basetmp();
}

TEST(MeasureOverlapProteinInterfaces, EmptyInterfaceListsDoNothing)
{
    test_mopi_empty_interface_lists_do_nothing();
}

TEST(MeasureOverlapProteinInterfaces, ArgumentsArePassedByValue)
{
    test_mopi_arguments_are_passed_by_value();
}