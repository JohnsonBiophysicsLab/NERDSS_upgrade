/*! \file test_determine_parent_complex.cpp
 *
 * ### Unit test for src/reactions/determine_parent_complex.cpp
 *
 * Function under test:
 *
 *     bool determine_parent_complex(int pro1Index, int pro2Index, int newComIndex,
 *                                  std::vector<Molecule>& moleculeList,
 *                                  std::vector<Complex>& complexList)
 *
 * The routine is called right after a bond between `pro1Index` and `pro2Index`
 * has been broken. It has to decide whether the parent Complex actually falls
 * apart into two pieces, or whether the two molecules are still connected
 * through some other path (a closed loop / a doubly bound pair). Behaviour:
 *
 *   * returns **true**  -> the complex stays together. Both dissociating
 *                          molecules keep the original complex index and the
 *                          member lists are left untouched.
 *   * returns **false** -> the complex is split. `complexList[c1]` keeps the
 *                          molecules connected to pro1, `complexList[newComIndex]`
 *                          receives the molecules connected to pro2, member
 *                          lists are sorted, `Complex::index` of the new
 *                          complex is set, and every molecule's `myComIndex`
 *                          is updated.
 *
 * Connectivity is expressed purely through `Molecule::bndpartner`, where the
 * value -1 marks a slot whose bond has been broken. The tests below therefore
 * build tiny synthetic topologies by hand and check both return value and all
 * side effects.
 *
 * NOTE: only `classes/class_Molecule_Complex.hpp` is included; the function
 * itself is declared locally (identical signature to the one in
 * include/reactions/unimolecular/unimolecular_reactions.hpp) so that this
 * translation unit does not have to drag in the whole reaction/SimulVolume
 * header chain.
 */

#include "classes/class_Molecule_Complex.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// Declaration of the function under test (see note above).
bool determine_parent_complex(int pro1Index, int pro2Index, int newComIndex,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList);

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace -> no symbol collisions with other tests)
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal Molecule for connectivity testing.
 *
 * \param[in] index    the molecule's own index inside moleculeList
 * \param[in] comIndex the index of the Complex this molecule belongs to
 * \param[in] partners the bndpartner list (-1 marks a broken/empty slot)
 */
Molecule dpc_make_molecule(int index, int comIndex, const std::vector<int>& partners)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.bndpartner = partners;
    mol.comCoord = Coord { 0.0, 0.0, 0.0 };
    return mol;
}

/*! \brief Build a minimal Complex holding the given member molecule indices. */
Complex dpc_make_complex(int index, const std::vector<int>& members)
{
    Complex com;
    com.index = index;
    com.memberList = members;
    return com;
}

/*! \brief Print and compare a member list against the expected content.
 *
 * Uses non-fatal EXPECT_* assertions so that all checks in a test run even if
 * one of them fails.
 */
void dpc_check_member_list(const std::vector<int>& actual, const std::vector<int>& expected,
    const std::string& label)
{
    std::cerr << "    " << label << " = [";
    for (size_t i = 0; i < actual.size(); ++i)
        std::cerr << (i ? " " : "") << actual[i];
    std::cerr << "]  expected = [";
    for (size_t i = 0; i < expected.size(); ++i)
        std::cerr << (i ? " " : "") << expected[i];
    std::cerr << "]\n";

    EXPECT_EQ(actual.size(), expected.size()) << label << ": wrong number of members";
    if (actual.size() == expected.size()) {
        for (size_t i = 0; i < actual.size(); ++i)
            EXPECT_EQ(actual[i], expected[i]) << label << ": member " << i << " differs";
    }
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: A simple dimer (0-1) whose only bond was just broken must split.
// -----------------------------------------------------------------------------
void test_dpc_simple_dimer_splits()
{
    std::cerr << "\n[TEST] test_dpc_simple_dimer_splits\n"
              << "  Source file:   src/reactions/determine_parent_complex.cpp\n"
              << "  Function:      determine_parent_complex()\n"
              << "  Topology:      dimer 0-1, the single bond has been broken\n"
              << "  Pass criteria: returns false (complex splits), complex 0 keeps\n"
              << "                 molecule 0, complex 1 receives molecule 1, and\n"
              << "                 both myComIndex values are updated.\n";

    // Two molecules, no remaining bonds (empty bndpartner lists).
    std::vector<Molecule> moleculeList {
        dpc_make_molecule(0, 0, {}),
        dpc_make_molecule(1, 0, {})
    };

    // Complex 0 currently owns both molecules; complex 1 is the empty target.
    std::vector<Complex> complexList {
        dpc_make_complex(0, { 0, 1 }),
        dpc_make_complex(1, {})
    };

    const bool keptTogether = determine_parent_complex(0, 1, 1, moleculeList, complexList);
    std::cerr << "  Return value = " << std::boolalpha << keptTogether << "\n";

    // The two molecules are no longer connected -> the complex must be split.
    EXPECT_FALSE(keptTogether) << "A broken dimer must be reported as split (false)";

    dpc_check_member_list(complexList[0].memberList, { 0 }, "complexList[0].memberList");
    dpc_check_member_list(complexList[1].memberList, { 1 }, "complexList[1].memberList");

    EXPECT_EQ(complexList[1].index, 1) << "The new complex must have its index set to newComIndex";
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "Molecule 0 must stay in the parent complex";
    EXPECT_EQ(moleculeList[1].myComIndex, 1) << "Molecule 1 must move to the new complex";
}

// -----------------------------------------------------------------------------
// Test 2: -1 entries in bndpartner are sentinels for broken bonds and must not
//         be mistaken for a doubly bound pair.
// -----------------------------------------------------------------------------
void test_dpc_sentinel_partners_are_ignored()
{
    std::cerr << "\n[TEST] test_dpc_sentinel_partners_are_ignored\n"
              << "  Source file:   src/reactions/determine_parent_complex.cpp\n"
              << "  Function:      determine_parent_complex()\n"
              << "  Topology:      dimer 0-1 where both molecules keep two -1 slots\n"
              << "  Pass criteria: the duplicate-partner test skips -1 entries, so the\n"
              << "                 complex still splits (returns false).\n";

    // Both molecules have two empty (-1) interface slots. If the -1 guard in the
    // 'doubly bound' test were missing, these would look like duplicate partners.
    std::vector<Molecule> moleculeList {
        dpc_make_molecule(0, 0, { -1, -1 }),
        dpc_make_molecule(1, 0, { -1, -1 })
    };
    std::vector<Complex> complexList {
        dpc_make_complex(0, { 0, 1 }),
        dpc_make_complex(1, {})
    };

    const bool keptTogether = determine_parent_complex(0, 1, 1, moleculeList, complexList);
    std::cerr << "  Return value = " << std::boolalpha << keptTogether << "\n";

    EXPECT_FALSE(keptTogether) << "-1 sentinels must not trigger the doubly-bound branch";
    dpc_check_member_list(complexList[0].memberList, { 0 }, "complexList[0].memberList");
    dpc_check_member_list(complexList[1].memberList, { 1 }, "complexList[1].memberList");
}

// -----------------------------------------------------------------------------
// Test 3: pro1/pro2 are asymmetric - the FIRST index keeps the parent complex.
// -----------------------------------------------------------------------------
void test_dpc_reversed_indices_split()
{
    std::cerr << "\n[TEST] test_dpc_reversed_indices_split\n"
              << "  Source file:   src/reactions/determine_parent_complex.cpp\n"
              << "  Function:      determine_parent_complex()\n"
              << "  Topology:      dimer 0-1, called with pro1=1 and pro2=0\n"
              << "  Pass criteria: molecule 1 (pro1) keeps complex 0, molecule 0 (pro2)\n"
              << "                 moves into the new complex 1.\n";

    std::vector<Molecule> moleculeList {
        dpc_make_molecule(0, 0, { -1 }),
        dpc_make_molecule(1, 0, { -1 })
    };
    std::vector<Complex> complexList {
        dpc_make_complex(0, { 0, 1 }),
        dpc_make_complex(1, {})
    };

    // Note the swapped order of the two dissociating molecules.
    const bool keptTogether = determine_parent_complex(1, 0, 1, moleculeList, complexList);
    std::cerr << "  Return value = " << std::boolalpha << keptTogether << "\n";

    EXPECT_FALSE(keptTogether) << "The broken dimer must still be split";
    dpc_check_member_list(complexList[0].memberList, { 1 }, "complexList[0].memberList");
    dpc_check_member_list(complexList[1].memberList, { 0 }, "complexList[1].memberList");
    EXPECT_EQ(moleculeList[1].myComIndex, 0) << "pro1 (molecule 1) must keep the parent complex";
    EXPECT_EQ(moleculeList[0].myComIndex, 1) << "pro2 (molecule 0) must move to the new complex";
}

// -----------------------------------------------------------------------------
// Test 4: linear trimer 0-1-2, break 0-1: molecule 2 has to follow molecule 1.
// -----------------------------------------------------------------------------
void test_dpc_chain_trimer_splits()
{
    std::cerr << "\n[TEST] test_dpc_chain_trimer_splits\n"
              << "  Source file:   src/reactions/determine_parent_complex.cpp\n"
              << "  Function:      determine_parent_complex()\n"
              << "  Topology:      chain 0-1-2, the 0-1 bond has been broken\n"
              << "  Pass criteria: returns false; complex 0 = {0}, complex 1 = {1,2}\n"
              << "                 because molecule 2 is directly bound to pro2 (=1).\n";

    // 0: bond to 1 broken            -> {-1}
    // 1: bond to 0 broken, bond to 2 -> {-1, 2}
    // 2: bond to 1                   -> {1}
    std::vector<Molecule> moleculeList {
        dpc_make_molecule(0, 0, { -1 }),
        dpc_make_molecule(1, 0, { -1, 2 }),
        dpc_make_molecule(2, 0, { 1 })
    };
    std::vector<Complex> complexList {
        dpc_make_complex(0, { 0, 1, 2 }),
        dpc_make_complex(1, {})
    };

    const bool keptTogether = determine_parent_complex(0, 1, 1, moleculeList, complexList);
    std::cerr << "  Return value = " << std::boolalpha << keptTogether << "\n";

    EXPECT_FALSE(keptTogether) << "An open chain must split when the terminal bond breaks";
    dpc_check_member_list(complexList[0].memberList, { 0 }, "complexList[0].memberList");
    dpc_check_member_list(complexList[1].memberList, { 1, 2 }, "complexList[1].memberList");
    EXPECT_EQ(moleculeList[2].myComIndex, 1) << "Molecule 2 must follow molecule 1 into complex 1";
}

// -----------------------------------------------------------------------------
// Test 5: longer chain 0-1-2-3 - molecule 3 is only reachable through molecule 2
//         and therefore exercises the deeper "recheck" search of the routine.
// -----------------------------------------------------------------------------
void test_dpc_deep_search_assigns_distant_member()
{
    std::cerr << "\n[TEST] test_dpc_deep_search_assigns_distant_member\n"
              << "  Source file:   src/reactions/determine_parent_complex.cpp\n"
              << "  Function:      determine_parent_complex()\n"
              << "  Topology:      chain 0-1-2-3, the 0-1 bond has been broken\n"
              << "  Pass criteria: returns false; molecule 3, which is bound to neither\n"
              << "                 dissociating molecule, is found by the recursive\n"
              << "                 'recheck' search and placed in complex 1 = {1,2,3}.\n";

    std::vector<Molecule> moleculeList {
        dpc_make_molecule(0, 0, { -1 }), // 0: only bond (to 1) is broken
        dpc_make_molecule(1, 0, { -1, 2 }), // 1: broken bond to 0, bond to 2
        dpc_make_molecule(2, 0, { 1, 3 }), // 2: bonds to 1 and 3
        dpc_make_molecule(3, 0, { 2 }) // 3: bond to 2 only -> needs deep search
    };
    std::vector<Complex> complexList {
        dpc_make_complex(0, { 0, 1, 2, 3 }),
        dpc_make_complex(1, {})
    };

    const bool keptTogether = determine_parent_complex(0, 1, 1, moleculeList, complexList);
    std::cerr << "  Return value = " << std::boolalpha << keptTogether << "\n";

    EXPECT_FALSE(keptTogether) << "The open chain must be reported as split";
    dpc_check_member_list(complexList[0].memberList, { 0 }, "complexList[0].memberList");
    dpc_check_member_list(complexList[1].memberList, { 1, 2, 3 }, "complexList[1].memberList");
    EXPECT_EQ(moleculeList[3].myComIndex, 1)
        << "Indirectly connected molecule 3 must end up in the new complex";
}

// -----------------------------------------------------------------------------
// Test 6: a third molecule bound to BOTH dissociating molecules keeps everything
//         inside one complex (first early-exit branch of the routine).
// -----------------------------------------------------------------------------
void test_dpc_shared_partner_keeps_one_complex()
{
    std::cerr << "\n[TEST] test_dpc_shared_partner_keeps_one_complex\n"
              << "  Source file:   src/reactions/determine_parent_complex.cpp\n"
              << "  Function:      determine_parent_complex()\n"
              << "  Topology:      molecule 2 is bound to both 0 and 1 (triangle),\n"
              << "                 the 0-1 bond has been broken\n"
              << "  Pass criteria: returns true, both dissociating molecules keep\n"
              << "                 complex 0 and no member list is modified.\n";

    std::vector<Molecule> moleculeList {
        dpc_make_molecule(0, 0, { -1, 2 }), // broken bond to 1, still bound to 2
        dpc_make_molecule(1, 0, { -1, 2 }), // broken bond to 0, still bound to 2
        dpc_make_molecule(2, 0, { 0, 1 }) // bound to both dissociating molecules
    };
    std::vector<Complex> complexList {
        dpc_make_complex(0, { 0, 1, 2 }),
        dpc_make_complex(1, {})
    };

    const bool keptTogether = determine_parent_complex(0, 1, 1, moleculeList, complexList);
    std::cerr << "  Return value = " << std::boolalpha << keptTogether << "\n";

    EXPECT_TRUE(keptTogether) << "A shared binding partner must keep the complex intact";
    dpc_check_member_list(complexList[0].memberList, { 0, 1, 2 }, "complexList[0].memberList");
    EXPECT_TRUE(complexList[1].memberList.empty())
        << "The new complex must stay empty when nothing dissociates";
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "Molecule 0 must remain in complex 0";
    EXPECT_EQ(moleculeList[1].myComIndex, 0)
        << "Molecule 1 must be reset back to complex 0 (it is pre-set to newComIndex)";
}

// -----------------------------------------------------------------------------
// Test 7: a doubly bound pair somewhere in the complex (the same partner listed
//         twice in bndpartner) also forces the complex to stay together.
// -----------------------------------------------------------------------------
void test_dpc_doubly_bound_pair_keeps_one_complex()
{
    std::cerr << "\n[TEST] test_dpc_doubly_bound_pair_keeps_one_complex\n"
              << "  Source file:   src/reactions/determine_parent_complex.cpp\n"
              << "  Function:      determine_parent_complex()\n"
              << "  Topology:      chain 0-1-2=3 where 2 and 3 are bound twice,\n"
              << "                 the 0-1 bond has been broken\n"
              << "  Pass criteria: the duplicate entry in bndpartner is detected, the\n"
              << "                 function returns true and nothing is split.\n";

    std::vector<Molecule> moleculeList {
        dpc_make_molecule(0, 0, { -1 }),
        dpc_make_molecule(1, 0, { -1, 2 }),
        dpc_make_molecule(2, 0, { 1, 3, 3 }), // molecule 3 appears twice
        dpc_make_molecule(3, 0, { 2, 2 }) // molecule 2 appears twice
    };
    std::vector<Complex> complexList {
        dpc_make_complex(0, { 0, 1, 2, 3 }),
        dpc_make_complex(1, {})
    };

    const bool keptTogether = determine_parent_complex(0, 1, 1, moleculeList, complexList);
    std::cerr << "  Return value = " << std::boolalpha << keptTogether << "\n";

    EXPECT_TRUE(keptTogether) << "A doubly bound pair must keep the complex as one unit";
    dpc_check_member_list(complexList[0].memberList, { 0, 1, 2, 3 }, "complexList[0].memberList");
    EXPECT_TRUE(complexList[1].memberList.empty()) << "The new complex must stay empty";
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "Molecule 0 must remain in complex 0";
    EXPECT_EQ(moleculeList[1].myComIndex, 0) << "Molecule 1 must be reset back to complex 0";
}

// -----------------------------------------------------------------------------
// Test 8: a closed ring 0-1-2-3-0. Breaking one bond leaves the ring connected,
//         which is only discovered by the shared-partner check near the end of
//         the routine.
// -----------------------------------------------------------------------------
void test_dpc_ring_stays_one_complex()
{
    std::cerr << "\n[TEST] test_dpc_ring_stays_one_complex\n"
              << "  Source file:   src/reactions/determine_parent_complex.cpp\n"
              << "  Function:      determine_parent_complex()\n"
              << "  Topology:      ring 0-1-2-3-0, the 0-1 bond has been broken\n"
              << "  Pass criteria: the remaining path 1-2-3-0 is detected as a closed\n"
              << "                 loop, so the function returns true and the member\n"
              << "                 lists are unchanged.\n";

    std::vector<Molecule> moleculeList {
        dpc_make_molecule(0, 0, { -1, 3 }), // broken bond to 1, bond to 3
        dpc_make_molecule(1, 0, { -1, 2 }), // broken bond to 0, bond to 2
        dpc_make_molecule(2, 0, { 1, 3 }),
        dpc_make_molecule(3, 0, { 2, 0 })
    };
    std::vector<Complex> complexList {
        dpc_make_complex(0, { 0, 1, 2, 3 }),
        dpc_make_complex(1, {})
    };

    const bool keptTogether = determine_parent_complex(0, 1, 1, moleculeList, complexList);
    std::cerr << "  Return value = " << std::boolalpha << keptTogether << "\n";

    EXPECT_TRUE(keptTogether) << "Breaking one bond of a ring must not split the complex";
    dpc_check_member_list(complexList[0].memberList, { 0, 1, 2, 3 }, "complexList[0].memberList");
    EXPECT_TRUE(complexList[1].memberList.empty()) << "The new complex must stay empty for a ring";
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "Molecule 0 must remain in complex 0";
    EXPECT_EQ(moleculeList[1].myComIndex, 0) << "Molecule 1 must be reset back to complex 0";
    // The other ring members must not have been touched either.
    EXPECT_EQ(moleculeList[2].myComIndex, 0) << "Molecule 2 must remain in complex 0";
    EXPECT_EQ(moleculeList[3].myComIndex, 0) << "Molecule 3 must remain in complex 0";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each scenario gets its own TEST so that the framework
// reports individual results while all scenarios are still executed.
// -----------------------------------------------------------------------------
TEST(DetermineParentComplex, SimpleDimerSplits) { test_dpc_simple_dimer_splits(); }
TEST(DetermineParentComplex, SentinelPartnersAreIgnored) { test_dpc_sentinel_partners_are_ignored(); }
TEST(DetermineParentComplex, ReversedIndicesSplit) { test_dpc_reversed_indices_split(); }
TEST(DetermineParentComplex, ChainTrimerSplits) { test_dpc_chain_trimer_splits(); }
TEST(DetermineParentComplex, DeepSearchAssignsDistantMember) { test_dpc_deep_search_assigns_distant_member(); }
TEST(DetermineParentComplex, SharedPartnerKeepsOneComplex) { test_dpc_shared_partner_keeps_one_complex(); }
TEST(DetermineParentComplex, DoublyBoundPairKeepsOneComplex) { test_dpc_doubly_bound_pair_keeps_one_complex(); }
TEST(DetermineParentComplex, RingStaysOneComplex) { test_dpc_ring_stays_one_complex(); }