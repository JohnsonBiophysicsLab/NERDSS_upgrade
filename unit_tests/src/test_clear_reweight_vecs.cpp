/*! \file test_clear_reweight_vecs.cpp
 *
 * ### Unit test for ../src/trajectory_functions/clear_reweight_vecs.cpp
 *
 * The function under test is:
 *
 *     void clear_reweight_vecs(Molecule& oneMol)
 *
 * Its behaviour (with the "write reweight factor" diagnostic block commented
 * out in the source) is a pure book-keeping operation on the six reweighting
 * vectors carried by every Molecule:
 *
 *     prevlist  <- currlist          currlist      cleared
 *     prevmyface<- currmyface        currmyface    cleared
 *     prevpface <- currpface         currpface     cleared
 *     prevnorm  <- currprevnorm      currprevnorm  cleared
 *     ps_prev   <- currps_prev       currps_prev   cleared
 *     prevsep   <- currprevsep       currprevsep   cleared
 *
 * Note the slightly asymmetric member names: the "current" partners of
 * `prevnorm`, `ps_prev` and `prevsep` are `currprevnorm`, `currps_prev` and
 * `currprevsep` respectively. The tests below verify the mapping explicitly so
 * that a future rename/mis-wiring would be caught.
 *
 * Pass criteria used throughout:
 *   - after the call, each prev* vector equals the pre-call value of its
 *     matching curr* vector (element-by-element, and same size);
 *   - after the call, every curr* vector is empty;
 *   - nothing else on the Molecule is modified.
 *
 * Verbose progress messages are written to stderr so a reader of the test log
 * can see exactly which source file / function is under test and what each
 * assertion is checking.
 */

#include "trajectory_functions/trajectory_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers. All names are prefixed with `crv_` (clear_reweight_vecs) so
// they cannot collide with helpers defined by other test translation units that
// are linked into the same gtest binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Pretty-print an int vector to stderr for the test log. */
void crv_dump_int_vec(const std::string& name, const std::vector<int>& v)
{
    std::cerr << "      " << name << " (size " << v.size() << ") = [";
    for (std::size_t i = 0; i < v.size(); ++i) {
        std::cerr << v[i];
        if (i + 1 < v.size())
            std::cerr << ", ";
    }
    std::cerr << "]\n";
}

/*! \brief Pretty-print a double vector to stderr for the test log. */
void crv_dump_dbl_vec(const std::string& name, const std::vector<double>& v)
{
    std::cerr << "      " << name << " (size " << v.size() << ") = [";
    for (std::size_t i = 0; i < v.size(); ++i) {
        std::cerr << v[i];
        if (i + 1 < v.size())
            std::cerr << ", ";
    }
    std::cerr << "]\n";
}

/*!
 * \brief Compare two int vectors element-by-element with non-fatal assertions.
 *
 * EXPECT_* (rather than ASSERT_*) is used everywhere so that a failure in one
 * comparison does not abort the remaining checks or the rest of the suite.
 */
void crv_expect_int_vec_eq(const std::string& what, const std::vector<int>& actual,
    const std::vector<int>& expected)
{
    EXPECT_EQ(actual.size(), expected.size())
        << what << ": vector sizes should match after the copy";

    // Only walk the overlapping region so a size mismatch cannot read past the
    // end of either vector.
    const std::size_t n = std::min(actual.size(), expected.size());
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_EQ(actual[i], expected[i])
            << what << ": element " << i << " should have been copied verbatim";
    }
}

/*! \brief Compare two double vectors element-by-element (exact bit equality).
 *
 * The function under test performs a plain vector assignment, so the values are
 * copied bit-for-bit; EXPECT_DOUBLE_EQ is therefore an appropriate criterion
 * (no arithmetic is applied that could introduce rounding).
 */
void crv_expect_dbl_vec_eq(const std::string& what, const std::vector<double>& actual,
    const std::vector<double>& expected)
{
    EXPECT_EQ(actual.size(), expected.size())
        << what << ": vector sizes should match after the copy";

    const std::size_t n = std::min(actual.size(), expected.size());
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_DOUBLE_EQ(actual[i], expected[i])
            << what << ": element " << i << " should have been copied verbatim";
    }
}

/*!
 * \brief Build a fully initialised Molecule whose curr* reweighting vectors
 *        hold easily recognisable, mutually distinct values.
 *
 * Every field the function touches is populated, and a handful of unrelated
 * fields are set as well so we can later prove they were left alone.
 *
 * The three integer vectors get disjoint value ranges (100s, 200s, 300s) and
 * the three double vectors get disjoint magnitudes (0.x, 10.x, 100.x) so that a
 * cross-wired assignment inside the function (e.g. prevnorm <- currps_prev)
 * would be immediately visible.
 */
Molecule crv_make_molecule_with_curr_data()
{
    Molecule mol;

    // ---- identity / bookkeeping fields (must survive the call untouched) ----
    mol.index = 7;
    mol.id = 7;
    mol.myComIndex = 3;
    mol.molTypeIndex = 1;
    mol.mySubVolIndex = 11;
    mol.mass = 2.5;
    mol.isLipid = false;
    mol.isEmpty = false;
    mol.isImplicitLipid = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = Coord{ 1.25, -2.5, 3.75 };

    // ---- "current" reweighting data: the inputs to the function ----
    mol.currlist = { 101, 102, 103 };     // partner molecule indices
    mol.currmyface = { 201, 202, 203 };   // this molecule's reacting interfaces
    mol.currpface = { 301, 302, 303 };    // partner's reacting interfaces
    mol.currprevnorm = { 0.1, 0.2, 0.3 }; // normalisation factors
    mol.currps_prev = { 10.1, 10.2, 10.3 };   // survival probabilities
    mol.currprevsep = { 100.1, 100.2, 100.3 }; // separations

    // ---- "previous" reweighting data: stale values that must be overwritten ----
    mol.prevlist = { -1, -2 };
    mol.prevmyface = { -3, -4 };
    mol.prevpface = { -5, -6 };
    mol.prevnorm = { -0.5, -0.6 };
    mol.ps_prev = { -0.7, -0.8 };
    mol.prevsep = { -0.9, -1.0 };

    return mol;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: The six curr* vectors are copied into their prev* counterparts, and
//         the curr* vectors are emptied.
// -----------------------------------------------------------------------------
void crv_test_copies_curr_into_prev_and_clears()
{
    std::cerr << "\n[TEST] crv_test_copies_curr_into_prev_and_clears\n"
              << "  Source file:   src/trajectory_functions/clear_reweight_vecs.cpp\n"
              << "  Function:      clear_reweight_vecs(Molecule&)\n"
              << "  Scenario:      a molecule carrying three current reweighting\n"
              << "                 entries and two stale previous entries.\n"
              << "  Pass criteria: each prev* vector equals the pre-call curr*\n"
              << "                 vector, and every curr* vector ends empty.\n";

    Molecule mol = crv_make_molecule_with_curr_data();

    // Snapshot the "current" data before the call so we can compare afterwards.
    const std::vector<int> expectedList = mol.currlist;
    const std::vector<int> expectedMyFace = mol.currmyface;
    const std::vector<int> expectedPFace = mol.currpface;
    const std::vector<double> expectedNorm = mol.currprevnorm;
    const std::vector<double> expectedPsPrev = mol.currps_prev;
    const std::vector<double> expectedSep = mol.currprevsep;

    std::cerr << "  Input state (before the call):\n";
    crv_dump_int_vec("currlist    ", mol.currlist);
    crv_dump_int_vec("currmyface  ", mol.currmyface);
    crv_dump_int_vec("currpface   ", mol.currpface);
    crv_dump_dbl_vec("currprevnorm", mol.currprevnorm);
    crv_dump_dbl_vec("currps_prev ", mol.currps_prev);
    crv_dump_dbl_vec("currprevsep ", mol.currprevsep);

    std::cerr << "  Calling clear_reweight_vecs...\n";
    clear_reweight_vecs(mol);

    std::cerr << "  Output state (after the call):\n";
    crv_dump_int_vec("prevlist    ", mol.prevlist);
    crv_dump_int_vec("prevmyface  ", mol.prevmyface);
    crv_dump_int_vec("prevpface   ", mol.prevpface);
    crv_dump_dbl_vec("prevnorm    ", mol.prevnorm);
    crv_dump_dbl_vec("ps_prev     ", mol.ps_prev);
    crv_dump_dbl_vec("prevsep     ", mol.prevsep);

    // --- the six copies: prev* must now hold the old curr* contents ---
    crv_expect_int_vec_eq("prevlist <- currlist", mol.prevlist, expectedList);
    crv_expect_int_vec_eq("prevmyface <- currmyface", mol.prevmyface, expectedMyFace);
    crv_expect_int_vec_eq("prevpface <- currpface", mol.prevpface, expectedPFace);
    crv_expect_dbl_vec_eq("prevnorm <- currprevnorm", mol.prevnorm, expectedNorm);
    crv_expect_dbl_vec_eq("ps_prev <- currps_prev", mol.ps_prev, expectedPsPrev);
    crv_expect_dbl_vec_eq("prevsep <- currprevsep", mol.prevsep, expectedSep);

    // --- the six clears: every curr* vector must be empty ---
    EXPECT_TRUE(mol.currlist.empty()) << "currlist should be cleared";
    EXPECT_TRUE(mol.currmyface.empty()) << "currmyface should be cleared";
    EXPECT_TRUE(mol.currpface.empty()) << "currpface should be cleared";
    EXPECT_TRUE(mol.currprevnorm.empty()) << "currprevnorm should be cleared";
    EXPECT_TRUE(mol.currps_prev.empty()) << "currps_prev should be cleared";
    EXPECT_TRUE(mol.currprevsep.empty()) << "currprevsep should be cleared";

    // Report the sizes explicitly as well, which makes a partial clear obvious.
    EXPECT_EQ(mol.currlist.size(), 0u) << "currlist size should be 0";
    EXPECT_EQ(mol.currmyface.size(), 0u) << "currmyface size should be 0";
    EXPECT_EQ(mol.currpface.size(), 0u) << "currpface size should be 0";
    EXPECT_EQ(mol.currprevnorm.size(), 0u) << "currprevnorm size should be 0";
    EXPECT_EQ(mol.currps_prev.size(), 0u) << "currps_prev size should be 0";
    EXPECT_EQ(mol.currprevsep.size(), 0u) << "currprevsep size should be 0";
}

// -----------------------------------------------------------------------------
// Test 2: The copy must be a true overwrite, not an append. The stale prev*
//         values present before the call must be gone afterwards.
// -----------------------------------------------------------------------------
void crv_test_overwrites_stale_prev_values()
{
    std::cerr << "\n[TEST] crv_test_overwrites_stale_prev_values\n"
              << "  Source file:   src/trajectory_functions/clear_reweight_vecs.cpp\n"
              << "  Function:      clear_reweight_vecs(Molecule&)\n"
              << "  Scenario:      prev* already contains 2 stale entries while\n"
              << "                 curr* contains 3 fresh entries.\n"
              << "  Pass criteria: prev* ends with exactly 3 entries (assignment,\n"
              << "                 not append) and none of the sentinel negative\n"
              << "                 stale values survive.\n";

    Molecule mol = crv_make_molecule_with_curr_data();

    // Deliberately different sizes: 2 stale entries vs 3 fresh entries.
    ASSERT_EQ(mol.prevlist.size(), 2u) << "test fixture should start with 2 stale entries";
    ASSERT_EQ(mol.currlist.size(), 3u) << "test fixture should start with 3 fresh entries";

    std::cerr << "  Stale previous data before the call:\n";
    crv_dump_int_vec("prevlist    ", mol.prevlist);
    crv_dump_dbl_vec("prevnorm    ", mol.prevnorm);

    std::cerr << "  Calling clear_reweight_vecs...\n";
    clear_reweight_vecs(mol);

    // Assignment semantics: prev* must now have the size of the old curr*, i.e.
    // 3, not 2 (unchanged) and not 5 (appended).
    EXPECT_EQ(mol.prevlist.size(), 3u)
        << "prevlist should be replaced (size 3), not appended to (size 5)";
    EXPECT_EQ(mol.prevmyface.size(), 3u) << "prevmyface should be replaced, size 3";
    EXPECT_EQ(mol.prevpface.size(), 3u) << "prevpface should be replaced, size 3";
    EXPECT_EQ(mol.prevnorm.size(), 3u) << "prevnorm should be replaced, size 3";
    EXPECT_EQ(mol.ps_prev.size(), 3u) << "ps_prev should be replaced, size 3";
    EXPECT_EQ(mol.prevsep.size(), 3u) << "prevsep should be replaced, size 3";

    // The stale values were all negative sentinels; none of the fresh values are
    // negative, so no element of any prev* vector may be negative now.
    for (std::size_t i = 0; i < mol.prevlist.size(); ++i)
        EXPECT_GT(mol.prevlist[i], 0) << "stale sentinel survived in prevlist[" << i << "]";
    for (std::size_t i = 0; i < mol.prevmyface.size(); ++i)
        EXPECT_GT(mol.prevmyface[i], 0) << "stale sentinel survived in prevmyface[" << i << "]";
    for (std::size_t i = 0; i < mol.prevpface.size(); ++i)
        EXPECT_GT(mol.prevpface[i], 0) << "stale sentinel survived in prevpface[" << i << "]";
    for (std::size_t i = 0; i < mol.prevnorm.size(); ++i)
        EXPECT_GT(mol.prevnorm[i], 0.0) << "stale sentinel survived in prevnorm[" << i << "]";
    for (std::size_t i = 0; i < mol.ps_prev.size(); ++i)
        EXPECT_GT(mol.ps_prev[i], 0.0) << "stale sentinel survived in ps_prev[" << i << "]";
    for (std::size_t i = 0; i < mol.prevsep.size(); ++i)
        EXPECT_GT(mol.prevsep[i], 0.0) << "stale sentinel survived in prevsep[" << i << "]";

    std::cerr << "  Previous data after the call:\n";
    crv_dump_int_vec("prevlist    ", mol.prevlist);
    crv_dump_dbl_vec("prevnorm    ", mol.prevnorm);
}

// -----------------------------------------------------------------------------
// Test 3: Verify that each prev* member is fed from the correct curr* member.
//         The curr/prev name pairing is not uniform in the Molecule struct
//         (currprevnorm -> prevnorm, currps_prev -> ps_prev,
//          currprevsep -> prevsep), so a mis-wiring is easy to introduce.
// -----------------------------------------------------------------------------
void crv_test_pairs_are_not_crosswired()
{
    std::cerr << "\n[TEST] crv_test_pairs_are_not_crosswired\n"
              << "  Source file:   src/trajectory_functions/clear_reweight_vecs.cpp\n"
              << "  Function:      clear_reweight_vecs(Molecule&)\n"
              << "  Scenario:      each curr* vector holds a distinct magnitude\n"
              << "                 band, so a swapped assignment is detectable.\n"
              << "  Pass criteria: prevnorm values lie in [0,1), ps_prev values in\n"
              << "                 [10,11), prevsep values in [100,101); and the\n"
              << "                 integer vectors keep their 100/200/300 bands.\n";

    Molecule mol = crv_make_molecule_with_curr_data();

    std::cerr << "  Calling clear_reweight_vecs...\n";
    clear_reweight_vecs(mol);

    // Integer vectors: 100-band, 200-band, 300-band respectively.
    for (std::size_t i = 0; i < mol.prevlist.size(); ++i) {
        EXPECT_GE(mol.prevlist[i], 100);
        EXPECT_LT(mol.prevlist[i], 200)
            << "prevlist[" << i << "] should come from currlist (100-band)";
    }
    for (std::size_t i = 0; i < mol.prevmyface.size(); ++i) {
        EXPECT_GE(mol.prevmyface[i], 200);
        EXPECT_LT(mol.prevmyface[i], 300)
            << "prevmyface[" << i << "] should come from currmyface (200-band)";
    }
    for (std::size_t i = 0; i < mol.prevpface.size(); ++i) {
        EXPECT_GE(mol.prevpface[i], 300);
        EXPECT_LT(mol.prevpface[i], 400)
            << "prevpface[" << i << "] should come from currpface (300-band)";
    }

    // Double vectors: sub-unit band, 10-band, 100-band respectively.
    for (std::size_t i = 0; i < mol.prevnorm.size(); ++i) {
        EXPECT_GE(mol.prevnorm[i], 0.0);
        EXPECT_LT(mol.prevnorm[i], 1.0)
            << "prevnorm[" << i << "] should come from currprevnorm (0.x band)";
    }
    for (std::size_t i = 0; i < mol.ps_prev.size(); ++i) {
        EXPECT_GE(mol.ps_prev[i], 10.0);
        EXPECT_LT(mol.ps_prev[i], 11.0)
            << "ps_prev[" << i << "] should come from currps_prev (10.x band)";
    }
    for (std::size_t i = 0; i < mol.prevsep.size(); ++i) {
        EXPECT_GE(mol.prevsep[i], 100.0);
        EXPECT_LT(mol.prevsep[i], 101.0)
            << "prevsep[" << i << "] should come from currprevsep (100.x band)";
    }

    std::cerr << "  All six prev* vectors fell inside their expected value bands.\n";
}

// -----------------------------------------------------------------------------
// Test 4: Calling the function on a molecule whose curr* vectors are already
//         empty must wipe the prev* vectors (they get assigned an empty vector).
//         This is the "molecule saw no encounters this timestep" case.
// -----------------------------------------------------------------------------
void crv_test_empty_curr_wipes_prev()
{
    std::cerr << "\n[TEST] crv_test_empty_curr_wipes_prev\n"
              << "  Source file:   src/trajectory_functions/clear_reweight_vecs.cpp\n"
              << "  Function:      clear_reweight_vecs(Molecule&)\n"
              << "  Scenario:      molecule had encounters last step (prev* full)\n"
              << "                 but none this step (curr* empty).\n"
              << "  Pass criteria: all prev* and all curr* vectors end empty.\n";

    Molecule mol = crv_make_molecule_with_curr_data();

    // Simulate "no encounters this timestep" by emptying every current vector.
    mol.currlist.clear();
    mol.currmyface.clear();
    mol.currpface.clear();
    mol.currprevnorm.clear();
    mol.currps_prev.clear();
    mol.currprevsep.clear();

    // Sanity check the fixture: prev* is still populated going in.
    ASSERT_FALSE(mol.prevlist.empty()) << "fixture should still hold stale prev data";

    std::cerr << "  Calling clear_reweight_vecs with empty curr* vectors...\n";
    clear_reweight_vecs(mol);

    // Assigning an empty vector must leave prev* empty too.
    EXPECT_TRUE(mol.prevlist.empty()) << "prevlist should be emptied by an empty currlist";
    EXPECT_TRUE(mol.prevmyface.empty()) << "prevmyface should be emptied";
    EXPECT_TRUE(mol.prevpface.empty()) << "prevpface should be emptied";
    EXPECT_TRUE(mol.prevnorm.empty()) << "prevnorm should be emptied";
    EXPECT_TRUE(mol.ps_prev.empty()) << "ps_prev should be emptied";
    EXPECT_TRUE(mol.prevsep.empty()) << "prevsep should be emptied";

    // And the curr* vectors obviously remain empty.
    EXPECT_TRUE(mol.currlist.empty()) << "currlist should remain empty";
    EXPECT_TRUE(mol.currmyface.empty()) << "currmyface should remain empty";
    EXPECT_TRUE(mol.currpface.empty()) << "currpface should remain empty";
    EXPECT_TRUE(mol.currprevnorm.empty()) << "currprevnorm should remain empty";
    EXPECT_TRUE(mol.currps_prev.empty()) << "currps_prev should remain empty";
    EXPECT_TRUE(mol.currprevsep.empty()) << "currprevsep should remain empty";

    std::cerr << "  All twelve reweighting vectors are empty, as expected.\n";
}

// -----------------------------------------------------------------------------
// Test 5: Calling the function twice in a row (as happens on consecutive
//         timesteps with no new encounters) must leave everything empty. The
//         second call is effectively idempotent on an already-drained molecule.
// -----------------------------------------------------------------------------
void crv_test_two_consecutive_calls_drain_everything()
{
    std::cerr << "\n[TEST] crv_test_two_consecutive_calls_drain_everything\n"
              << "  Source file:   src/trajectory_functions/clear_reweight_vecs.cpp\n"
              << "  Function:      clear_reweight_vecs(Molecule&)\n"
              << "  Scenario:      two timesteps in a row with no refill of curr*.\n"
              << "  Pass criteria: after call #1 prev* holds the data; after call\n"
              << "                 #2 both prev* and curr* are empty; a third call\n"
              << "                 changes nothing.\n";

    Molecule mol = crv_make_molecule_with_curr_data();
    const std::size_t nCurr = mol.currlist.size();

    // --- first call: data migrates curr* -> prev* ---
    std::cerr << "  Call #1 (curr* populated)...\n";
    clear_reweight_vecs(mol);
    EXPECT_EQ(mol.prevlist.size(), nCurr)
        << "after the first call prev* should hold the migrated data";
    EXPECT_TRUE(mol.currlist.empty()) << "after the first call curr* should be empty";

    // --- second call: nothing left to migrate, so prev* is wiped ---
    std::cerr << "  Call #2 (curr* now empty)...\n";
    clear_reweight_vecs(mol);
    EXPECT_TRUE(mol.prevlist.empty()) << "second call should wipe prevlist";
    EXPECT_TRUE(mol.prevmyface.empty()) << "second call should wipe prevmyface";
    EXPECT_TRUE(mol.prevpface.empty()) << "second call should wipe prevpface";
    EXPECT_TRUE(mol.prevnorm.empty()) << "second call should wipe prevnorm";
    EXPECT_TRUE(mol.ps_prev.empty()) << "second call should wipe ps_prev";
    EXPECT_TRUE(mol.prevsep.empty()) << "second call should wipe prevsep";

    // --- third call: fully drained molecule, must remain drained ---
    std::cerr << "  Call #3 (everything already empty; expect a no-op)...\n";
    clear_reweight_vecs(mol);
    EXPECT_TRUE(mol.prevlist.empty()) << "third call must not resurrect prevlist";
    EXPECT_TRUE(mol.currlist.empty()) << "third call must not resurrect currlist";
    EXPECT_TRUE(mol.prevnorm.empty()) << "third call must not resurrect prevnorm";
    EXPECT_TRUE(mol.currprevnorm.empty()) << "third call must not resurrect currprevnorm";

    std::cerr << "  Repeated calls on a drained molecule are stable.\n";
}

// -----------------------------------------------------------------------------
// Test 6: The function must touch only the twelve reweighting vectors. All the
//         other Molecule state (identity, coordinates, flags, interface list,
//         and the unrelated association vectors) must be preserved exactly.
// -----------------------------------------------------------------------------
void crv_test_leaves_unrelated_fields_untouched()
{
    std::cerr << "\n[TEST] crv_test_leaves_unrelated_fields_untouched\n"
              << "  Source file:   src/trajectory_functions/clear_reweight_vecs.cpp\n"
              << "  Function:      clear_reweight_vecs(Molecule&)\n"
              << "  Scenario:      a molecule with populated identity fields,\n"
              << "                 coordinates, an interface, and non-reweighting\n"
              << "                 association vectors (freelist/bndlist/probvec).\n"
              << "  Pass criteria: every one of those fields is bit-identical\n"
              << "                 after the call.\n";

    Molecule mol = crv_make_molecule_with_curr_data();

    // Populate a few fields that clear_reweight_vecs must never modify.
    Molecule::Iface iface;
    iface.coord = Coord{ 1.25, -2.5, 4.75 };
    iface.relIndex = 0;
    iface.index = 5;
    iface.molTypeIndex = mol.molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    mol.freelist = { 0 };
    mol.bndlist = { 4, 5 };
    mol.bndpartner = { 6 };
    mol.probvec = { 0.25, 0.5 };
    mol.crossbase = { 9 };
    mol.mycrossint = { 8 };

    // Snapshot everything we intend to verify.
    const int expIndex = mol.index;
    const int expId = mol.id;
    const int expComIndex = mol.myComIndex;
    const int expMolTypeIndex = mol.molTypeIndex;
    const int expSubVolIndex = mol.mySubVolIndex;
    const double expMass = mol.mass;
    const bool expIsLipid = mol.isLipid;
    const bool expIsEmpty = mol.isEmpty;
    const Coord expCom = mol.comCoord;
    const Coord expIfaceCoord = mol.interfaceList[0].coord;
    const std::vector<int> expFreeList = mol.freelist;
    const std::vector<int> expBndList = mol.bndlist;
    const std::vector<int> expBndPartner = mol.bndpartner;
    const std::vector<double> expProbVec = mol.probvec;
    const std::vector<int> expCrossBase = mol.crossbase;
    const std::vector<int> expMyCrossInt = mol.mycrossint;
    const std::size_t expIfaceCount = mol.interfaceList.size();

    std::cerr << "  Calling clear_reweight_vecs...\n";
    clear_reweight_vecs(mol);

    // --- scalar identity fields ---
    EXPECT_EQ(mol.index, expIndex) << "index must not be modified";
    EXPECT_EQ(mol.id, expId) << "id must not be modified";
    EXPECT_EQ(mol.myComIndex, expComIndex) << "myComIndex must not be modified";
    EXPECT_EQ(mol.molTypeIndex, expMolTypeIndex) << "molTypeIndex must not be modified";
    EXPECT_EQ(mol.mySubVolIndex, expSubVolIndex) << "mySubVolIndex must not be modified";
    EXPECT_DOUBLE_EQ(mol.mass, expMass) << "mass must not be modified";
    EXPECT_EQ(mol.isLipid, expIsLipid) << "isLipid must not be modified";
    EXPECT_EQ(mol.isEmpty, expIsEmpty) << "isEmpty must not be modified";

    // --- geometry ---
    EXPECT_DOUBLE_EQ(mol.comCoord.x, expCom.x) << "comCoord.x must not be modified";
    EXPECT_DOUBLE_EQ(mol.comCoord.y, expCom.y) << "comCoord.y must not be modified";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, expCom.z) << "comCoord.z must not be modified";

    EXPECT_EQ(mol.interfaceList.size(), expIfaceCount)
        << "the interface list must not gain or lose entries";
    if (mol.interfaceList.size() == expIfaceCount && expIfaceCount > 0) {
        EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.x, expIfaceCoord.x)
            << "interface coord.x must not be modified";
        EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.y, expIfaceCoord.y)
            << "interface coord.y must not be modified";
        EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.z, expIfaceCoord.z)
            << "interface coord.z must not be modified";
    }

    // --- other association vectors, which share the Molecule struct but are
    //     not part of the reweighting book-keeping ---
    crv_expect_int_vec_eq("freelist", mol.freelist, expFreeList);
    crv_expect_int_vec_eq("bndlist", mol.bndlist, expBndList);
    crv_expect_int_vec_eq("bndpartner", mol.bndpartner, expBndPartner);
    crv_expect_dbl_vec_eq("probvec", mol.probvec, expProbVec);
    crv_expect_int_vec_eq("crossbase", mol.crossbase, expCrossBase);
    crv_expect_int_vec_eq("mycrossint", mol.mycrossint, expMyCrossInt);

    std::cerr << "  All unrelated Molecule fields survived the call unchanged.\n";
}

// -----------------------------------------------------------------------------
// Test 7: Operating on one molecule out of a list must not disturb its
//         neighbours. clear_reweight_vecs takes a reference, so this guards
//         against accidental aliasing in future refactors.
// -----------------------------------------------------------------------------
void crv_test_only_target_molecule_is_affected()
{
    std::cerr << "\n[TEST] crv_test_only_target_molecule_is_affected\n"
              << "  Source file:   src/trajectory_functions/clear_reweight_vecs.cpp\n"
              << "  Function:      clear_reweight_vecs(Molecule&)\n"
              << "  Scenario:      a three-molecule list; only element 1 is passed.\n"
              << "  Pass criteria: element 1 is drained/migrated while elements 0\n"
              << "                 and 2 keep both their curr* and prev* data.\n";

    // Build a small list of identical molecules and give each a unique index.
    std::vector<Molecule> moleculeList;
    for (int i = 0; i < 3; ++i) {
        Molecule m = crv_make_molecule_with_curr_data();
        m.index = i;
        m.id = i;
        moleculeList.push_back(m);
    }

    const std::size_t nCurr = moleculeList[1].currlist.size();

    std::cerr << "  Calling clear_reweight_vecs on moleculeList[1] only...\n";
    clear_reweight_vecs(moleculeList[1]);

    // --- the targeted molecule migrated its data ---
    EXPECT_TRUE(moleculeList[1].currlist.empty())
        << "target molecule currlist should be cleared";
    EXPECT_EQ(moleculeList[1].prevlist.size(), nCurr)
        << "target molecule prevlist should hold the migrated data";

    // --- the untouched neighbours kept their original curr* data ---
    EXPECT_EQ(moleculeList[0].currlist.size(), nCurr)
        << "moleculeList[0] currlist should be untouched";
    EXPECT_EQ(moleculeList[2].currlist.size(), nCurr)
        << "moleculeList[2] currlist should be untouched";
    EXPECT_EQ(moleculeList[0].currprevnorm.size(), nCurr)
        << "moleculeList[0] currprevnorm should be untouched";
    EXPECT_EQ(moleculeList[2].currprevsep.size(), nCurr)
        << "moleculeList[2] currprevsep should be untouched";

    // --- and their stale prev* data (size 2 sentinels) is still stale ---
    EXPECT_EQ(moleculeList[0].prevlist.size(), 2u)
        << "moleculeList[0] prevlist should still hold its 2 stale entries";
    EXPECT_EQ(moleculeList[2].prevlist.size(), 2u)
        << "moleculeList[2] prevlist should still hold its 2 stale entries";

    std::cerr << "  Neighbouring molecules were correctly left alone.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper runs inside its own TEST so a
// failure in one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(ClearReweightVecs, CopiesCurrIntoPrevAndClears) { crv_test_copies_curr_into_prev_and_clears(); }
TEST(ClearReweightVecs, OverwritesStalePrevValues) { crv_test_overwrites_stale_prev_values(); }
TEST(ClearReweightVecs, PairsAreNotCrosswired) { crv_test_pairs_are_not_crosswired(); }
TEST(ClearReweightVecs, EmptyCurrWipesPrev) { crv_test_empty_curr_wipes_prev(); }
TEST(ClearReweightVecs, TwoConsecutiveCallsDrainEverything) { crv_test_two_consecutive_calls_drain_everything(); }
TEST(ClearReweightVecs, LeavesUnrelatedFieldsUntouched) { crv_test_leaves_unrelated_fields_untouched(); }
TEST(ClearReweightVecs, OnlyTargetMoleculeIsAffected) { crv_test_only_target_molecule_is_affected(); }