/*! \file test_determine_parent_complex_IL.cpp
 *
 * ### Unit test for src/reactions/determine_parent_complex_IL.cpp
 *
 * Function under test:
 *
 *     bool determine_parent_complex_IL(int pro1Index, int pro2Index,
 *                                      int newComIndex,
 *                                      std::vector<Molecule>& moleculeList,
 *                                      std::vector<Complex>& complexList,
 *                                      int ILindexMol);
 *
 * The routine is called after a bond between two molecules has been broken in a
 * system that contains an *implicit lipid* (IL).  Starting at `pro1Index` it
 * walks the bound-partner graph (`Molecule::bndpartner`), deliberately skipping
 * the implicit lipid molecule `ILindexMol` so that links through the membrane
 * do not count as physical connectivity.
 *
 *   * If `pro2Index` is reachable, the two molecules are still part of the same
 *     physical complex -> the function returns **true** and leaves both
 *     molecules assigned to the original complex index.
 *
 *   * If `pro2Index` is *not* reachable, the complex must be split -> the
 *     function returns **false**, keeps the connected component of `pro1Index`
 *     in the original complex, moves every remaining original member into
 *     `complexList[newComIndex]`, stamps that complex's `index`, and repairs
 *     `Molecule::myComIndex` for every member of both complexes.
 *
 * Each test below prints what it is doing and what the pass criterion is, and
 * uses non-fatal EXPECT_* assertions so that every test runs to completion.
 */

#include "reactions/unimolecular/unimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace so names cannot collide with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal Molecule for connectivity tests.
 *
 * Only the fields the function under test actually reads/writes are populated:
 * the molecule's own index, its parent complex index, and its bound partner
 * list.
 *
 * \param[in] index    index of this molecule in moleculeList
 * \param[in] comIndex index of the complex this molecule currently belongs to
 * \param[in] partners indices of the molecules this molecule is bound to
 */
Molecule dpcil_make_molecule(int index, int comIndex, const std::vector<int>& partners)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.bndpartner = partners;
    return mol;
}

/*! \brief Build a Complex holding the given member molecule indices. */
Complex dpcil_make_complex(int index, const std::vector<int>& members)
{
    Complex com;
    com.index = index;
    com.memberList = members;
    return com;
}

/*! \brief Convenience predicate: is `val` present in the vector? */
bool dpcil_contains(const std::vector<int>& vec, int val)
{
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

/*! \brief Print a member list to stderr for human inspection. */
void dpcil_print_list(const char* label, const std::vector<int>& vec)
{
    std::cerr << "    " << label << " = {";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i != 0)
            std::cerr << ", ";
        std::cerr << vec[i];
    }
    std::cerr << "}\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the two molecules remain connected through a third molecule.
//         Topology (after the bond 0-2 was broken):  0 - 1 - 2
// -----------------------------------------------------------------------------
void test_dpcil_still_connected_returns_true()
{
    std::cerr << "\n[TEST] test_dpcil_still_connected_returns_true\n"
              << "  Source file: src/reactions/determine_parent_complex_IL.cpp\n"
              << "  Function:    determine_parent_complex_IL\n"
              << "  Scenario:    molecules 0 and 2 are still linked via molecule 1\n"
              << "               (0-1-2 chain), so the complex must NOT be split.\n"
              << "  Pass:        returns true, both molecules keep the original\n"
              << "               complex index (0), and no complex is modified.\n";

    // moleculeList: 0<->1, 1<->2.  Molecule 3 is the implicit lipid (unused here).
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dpcil_make_molecule(0, 0, { 1 }));
    moleculeList.push_back(dpcil_make_molecule(1, 0, { 0, 2 }));
    moleculeList.push_back(dpcil_make_molecule(2, 0, { 1 }));
    moleculeList.push_back(dpcil_make_molecule(3, 1, {})); // implicit lipid

    // Complex 0 holds all three real molecules; complex 1 is the (empty) target.
    std::vector<Complex> complexList;
    complexList.push_back(dpcil_make_complex(0, { 0, 1, 2 }));
    complexList.push_back(dpcil_make_complex(1, {}));

    const int newComIndex = 1;
    const int ILindexMol = 3;

    std::cerr << "  Calling determine_parent_complex_IL(pro1=0, pro2=2, newCom=1, IL=3)\n";
    const bool stillConnected
        = determine_parent_complex_IL(0, 2, newComIndex, moleculeList, complexList, ILindexMol);

    // The molecules are connected, so the routine must report "true".
    EXPECT_TRUE(stillConnected)
        << "0 and 2 are linked via 1, so the complex should not be split";

    // Both reacting molecules must be restored to the original complex index.
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "pro1 must stay in complex 0";
    EXPECT_EQ(moleculeList[2].myComIndex, 0)
        << "pro2 must be put back into complex 0 (it was speculatively set to newComIndex)";

    // The original complex membership is untouched, and the new complex stays empty.
    dpcil_print_list("complexList[0].memberList", complexList[0].memberList);
    dpcil_print_list("complexList[1].memberList", complexList[1].memberList);
    EXPECT_EQ(complexList[0].memberList.size(), 3u)
        << "original complex should still hold all three molecules";
    EXPECT_TRUE(complexList[1].memberList.empty())
        << "new complex should remain empty when nothing was split off";
}

// -----------------------------------------------------------------------------
// Test 2: the two molecules are no longer connected -> the complex must split.
//         Topology after breaking: {0-1}  and  {2-3}
// -----------------------------------------------------------------------------
void test_dpcil_disconnected_splits_complex()
{
    std::cerr << "\n[TEST] test_dpcil_disconnected_splits_complex\n"
              << "  Source file: src/reactions/determine_parent_complex_IL.cpp\n"
              << "  Function:    determine_parent_complex_IL\n"
              << "  Scenario:    complex 0 originally held {0,1,2,3} but the bond\n"
              << "               between the 0-1 pair and the 2-3 pair is gone.\n"
              << "  Pass:        returns false; complex 0 keeps {0,1}; complex 1\n"
              << "               receives {2,3}; complex 1 index is stamped; and\n"
              << "               every molecule's myComIndex is repaired.\n";

    // Two disjoint dimers, all four molecules initially in complex 0.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dpcil_make_molecule(0, 0, { 1 }));
    moleculeList.push_back(dpcil_make_molecule(1, 0, { 0 }));
    moleculeList.push_back(dpcil_make_molecule(2, 0, { 3 }));
    moleculeList.push_back(dpcil_make_molecule(3, 0, { 2 }));
    moleculeList.push_back(dpcil_make_molecule(4, 1, {})); // implicit lipid

    std::vector<Complex> complexList;
    complexList.push_back(dpcil_make_complex(0, { 0, 1, 2, 3 }));
    complexList.push_back(dpcil_make_complex(-1, {})); // index intentionally wrong

    const int newComIndex = 1;
    const int ILindexMol = 4;

    std::cerr << "  Calling determine_parent_complex_IL(pro1=0, pro2=2, newCom=1, IL=4)\n";
    const bool stillConnected
        = determine_parent_complex_IL(0, 2, newComIndex, moleculeList, complexList, ILindexMol);

    // Molecule 2 cannot be reached from molecule 0 -> split required.
    EXPECT_FALSE(stillConnected)
        << "0 and 2 are in disjoint components, so the complex must be split";

    dpcil_print_list("complexList[0].memberList", complexList[0].memberList);
    dpcil_print_list("complexList[1].memberList", complexList[1].memberList);

    // Complex 0 keeps exactly the connected component containing molecule 0.
    EXPECT_EQ(complexList[0].memberList.size(), 2u)
        << "old complex should keep exactly two members";
    EXPECT_TRUE(dpcil_contains(complexList[0].memberList, 0))
        << "molecule 0 must stay in the old complex";
    EXPECT_TRUE(dpcil_contains(complexList[0].memberList, 1))
        << "molecule 1 (bound to 0) must stay in the old complex";

    // Complex 1 receives the remaining original members.
    EXPECT_EQ(complexList[1].memberList.size(), 2u)
        << "new complex should receive exactly two members";
    EXPECT_TRUE(dpcil_contains(complexList[1].memberList, 2))
        << "molecule 2 must move to the new complex";
    EXPECT_TRUE(dpcil_contains(complexList[1].memberList, 3))
        << "molecule 3 (bound to 2) must move to the new complex";

    // The new complex must know its own index.
    EXPECT_EQ(complexList[1].index, newComIndex)
        << "new complex index must be set to newComIndex";

    // Every molecule's parent-complex back-pointer must be consistent.
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "molecule 0 -> complex 0";
    EXPECT_EQ(moleculeList[1].myComIndex, 0) << "molecule 1 -> complex 0";
    EXPECT_EQ(moleculeList[2].myComIndex, newComIndex) << "molecule 2 -> complex 1";
    EXPECT_EQ(moleculeList[3].myComIndex, newComIndex) << "molecule 3 -> complex 1";
}

// -----------------------------------------------------------------------------
// Test 3: a path that runs *through the implicit lipid* must be ignored.
//         Topology: 0 - IL(2) - 1, with ILindexMol == 2.
// -----------------------------------------------------------------------------
void test_dpcil_implicit_lipid_link_is_ignored()
{
    std::cerr << "\n[TEST] test_dpcil_implicit_lipid_link_is_ignored\n"
              << "  Source file: src/reactions/determine_parent_complex_IL.cpp\n"
              << "  Function:    determine_parent_complex_IL\n"
              << "  Scenario:    molecules 0 and 1 are only linked through the\n"
              << "               implicit lipid (molecule 2), which the search\n"
              << "               must skip.\n"
              << "  Pass:        returns false (treated as disconnected) and the\n"
              << "               complex is split into {0} and {1}.\n";

    // Both proteins are bound to the implicit lipid (index 2) and to nothing else.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dpcil_make_molecule(0, 0, { 2 }));
    moleculeList.push_back(dpcil_make_molecule(1, 0, { 2 }));
    moleculeList.push_back(dpcil_make_molecule(2, 2, { 0, 1 })); // the implicit lipid

    std::vector<Complex> complexList;
    complexList.push_back(dpcil_make_complex(0, { 0, 1 }));
    complexList.push_back(dpcil_make_complex(1, {}));
    complexList.push_back(dpcil_make_complex(2, { 2 })); // IL's own complex

    const int newComIndex = 1;
    const int ILindexMol = 2;

    std::cerr << "  Calling determine_parent_complex_IL(pro1=0, pro2=1, newCom=1, IL=2)\n";
    const bool stillConnected
        = determine_parent_complex_IL(0, 1, newComIndex, moleculeList, complexList, ILindexMol);

    // The only path is through the IL, which is skipped -> disconnected.
    EXPECT_FALSE(stillConnected)
        << "a connection only via the implicit lipid must not count as connectivity";

    dpcil_print_list("complexList[0].memberList", complexList[0].memberList);
    dpcil_print_list("complexList[1].memberList", complexList[1].memberList);

    EXPECT_EQ(complexList[0].memberList.size(), 1u)
        << "old complex should keep only molecule 0";
    EXPECT_TRUE(dpcil_contains(complexList[0].memberList, 0))
        << "molecule 0 stays in the old complex";
    EXPECT_EQ(complexList[1].memberList.size(), 1u)
        << "new complex should receive only molecule 1";
    EXPECT_TRUE(dpcil_contains(complexList[1].memberList, 1))
        << "molecule 1 moves to the new complex";

    // Back-pointers repaired, and the IL itself must be untouched.
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "molecule 0 -> complex 0";
    EXPECT_EQ(moleculeList[1].myComIndex, newComIndex) << "molecule 1 -> complex 1";
    EXPECT_EQ(moleculeList[2].myComIndex, 2)
        << "the implicit lipid must not be reassigned by this routine";
}

// -----------------------------------------------------------------------------
// Test 4: control for test 3 -- if the intermediate molecule is *not* declared
//         to be the implicit lipid, the same topology is reported as connected.
// -----------------------------------------------------------------------------
void test_dpcil_non_il_intermediate_keeps_connection()
{
    std::cerr << "\n[TEST] test_dpcil_non_il_intermediate_keeps_connection\n"
              << "  Source file: src/reactions/determine_parent_complex_IL.cpp\n"
              << "  Function:    determine_parent_complex_IL\n"
              << "  Scenario:    identical topology to the previous test (0 - 2 - 1)\n"
              << "               but ILindexMol = -1, so molecule 2 is a normal\n"
              << "               molecule and provides a real path.\n"
              << "  Pass:        returns true and no complex is split.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dpcil_make_molecule(0, 0, { 2 }));
    moleculeList.push_back(dpcil_make_molecule(1, 0, { 2 }));
    moleculeList.push_back(dpcil_make_molecule(2, 0, { 0, 1 })); // ordinary bridge

    std::vector<Complex> complexList;
    complexList.push_back(dpcil_make_complex(0, { 0, 1, 2 }));
    complexList.push_back(dpcil_make_complex(1, {}));

    const int newComIndex = 1;
    const int ILindexMol = -1; // no implicit lipid in the system

    std::cerr << "  Calling determine_parent_complex_IL(pro1=0, pro2=1, newCom=1, IL=-1)\n";
    const bool stillConnected
        = determine_parent_complex_IL(0, 1, newComIndex, moleculeList, complexList, ILindexMol);

    EXPECT_TRUE(stillConnected)
        << "with no IL to skip, molecule 2 bridges 0 and 1 -> still one complex";

    // No split happened: memberships unchanged, new complex still empty.
    dpcil_print_list("complexList[0].memberList", complexList[0].memberList);
    dpcil_print_list("complexList[1].memberList", complexList[1].memberList);
    EXPECT_EQ(complexList[0].memberList.size(), 3u)
        << "original complex membership should be unchanged";
    EXPECT_TRUE(complexList[1].memberList.empty())
        << "new complex should stay empty";

    // Both reacting molecules point back at the original complex.
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "pro1 -> complex 0";
    EXPECT_EQ(moleculeList[1].myComIndex, 0) << "pro2 restored to complex 0";
}

// -----------------------------------------------------------------------------
// Test 5: a molecule listed in the original complex but not reachable from
//         either reactant (e.g. an already-detached fragment) is placed into the
//         new complex, because the routine assigns *everything* not in the
//         connected component of pro1 to newComIndex.
// -----------------------------------------------------------------------------
void test_dpcil_unreachable_member_goes_to_new_complex()
{
    std::cerr << "\n[TEST] test_dpcil_unreachable_member_goes_to_new_complex\n"
              << "  Source file: src/reactions/determine_parent_complex_IL.cpp\n"
              << "  Function:    determine_parent_complex_IL\n"
              << "  Scenario:    complex 0 lists {0,1,2}; molecule 0 is isolated,\n"
              << "               molecules 1 and 2 are bound to each other.\n"
              << "  Pass:        returns false; complex 0 keeps only {0}; both 1\n"
              << "               and 2 (everything outside pro1's component) end\n"
              << "               up in the new complex.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dpcil_make_molecule(0, 0, {})); // isolated reactant
    moleculeList.push_back(dpcil_make_molecule(1, 0, { 2 }));
    moleculeList.push_back(dpcil_make_molecule(2, 0, { 1 }));

    std::vector<Complex> complexList;
    complexList.push_back(dpcil_make_complex(0, { 0, 1, 2 }));
    complexList.push_back(dpcil_make_complex(0, {})); // wrong index on purpose

    const int newComIndex = 1;
    const int ILindexMol = -1;

    std::cerr << "  Calling determine_parent_complex_IL(pro1=0, pro2=1, newCom=1, IL=-1)\n";
    const bool stillConnected
        = determine_parent_complex_IL(0, 1, newComIndex, moleculeList, complexList, ILindexMol);

    EXPECT_FALSE(stillConnected)
        << "molecule 0 has no partners, so it cannot reach molecule 1";

    dpcil_print_list("complexList[0].memberList", complexList[0].memberList);
    dpcil_print_list("complexList[1].memberList", complexList[1].memberList);

    // pro1's component is just {0}.
    EXPECT_EQ(complexList[0].memberList.size(), 1u)
        << "old complex should contain only the isolated molecule 0";
    EXPECT_TRUE(dpcil_contains(complexList[0].memberList, 0))
        << "molecule 0 stays in the old complex";

    // Everything else from the original member list moves over.
    EXPECT_EQ(complexList[1].memberList.size(), 2u)
        << "new complex should absorb the remaining two original members";
    EXPECT_TRUE(dpcil_contains(complexList[1].memberList, 1))
        << "molecule 1 moves to the new complex";
    EXPECT_TRUE(dpcil_contains(complexList[1].memberList, 2))
        << "molecule 2 moves to the new complex";

    EXPECT_EQ(complexList[1].index, newComIndex)
        << "new complex index must be corrected to newComIndex";

    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "molecule 0 -> complex 0";
    EXPECT_EQ(moleculeList[1].myComIndex, newComIndex) << "molecule 1 -> complex 1";
    EXPECT_EQ(moleculeList[2].myComIndex, newComIndex) << "molecule 2 -> complex 1";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each helper is run inside its own TEST so that all of
// them execute (and report) independently even if some assertions fail.
// -----------------------------------------------------------------------------
TEST(DetermineParentComplexIL, StillConnectedReturnsTrue) { test_dpcil_still_connected_returns_true(); }
TEST(DetermineParentComplexIL, DisconnectedSplitsComplex) { test_dpcil_disconnected_splits_complex(); }
TEST(DetermineParentComplexIL, ImplicitLipidLinkIsIgnored) { test_dpcil_implicit_lipid_link_is_ignored(); }
TEST(DetermineParentComplexIL, NonILIntermediateKeepsConnection) { test_dpcil_non_il_intermediate_keeps_connection(); }
TEST(DetermineParentComplexIL, UnreachableMemberGoesToNewComplex) { test_dpcil_unreachable_member_goes_to_new_complex(); }