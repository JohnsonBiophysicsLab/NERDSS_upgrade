/*! \file test_debug_molecule_complex_missmatch.cpp
 *
 * ### Unit tests for src/debug/debug_molecule_complex_missmatch.cpp
 *
 * That translation unit provides two functions:
 *
 *   1. void debug_molecule_complex_missmatch(MpiContext&, vector<Molecule>&,
 *                                            vector<Complex>&, string)
 *      -- a consistency checker whose *entire* body is wrapped in `if (DEBUG)`.
 *         `DEBUG` is `#define DEBUG false` in include/macro.hpp, so in the build
 *         used by this test suite the function is compiled as a pure no-op.
 *         When DEBUG is enabled the function calls error(), which terminates the
 *         process, so the only behaviour that can be safely asserted here is the
 *         no-op behaviour: it must return, must not throw, and must not modify
 *         the Molecule/Complex containers even when they are inconsistent.
 *
 *   2. void check_complex_index(MpiContext&, vector<Molecule>&,
 *                               vector<Complex>&)
 *      -- always active. It walks moleculeList and throws std::runtime_error
 *         (after printing the offending molecule) whenever a non-empty molecule
 *         has myComIndex outside the range [0, complexList.size()).
 *
 * Pass criteria used below:
 *   * debug_molecule_complex_missmatch leaves every field it inspects untouched
 *     and never throws (with DEBUG == false).
 *   * check_complex_index throws std::runtime_error exactly for live molecules
 *     with an out-of-range complex index, and returns quietly otherwise
 *     (including for molecules flagged isEmpty, which it skips).
 */

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_SimulVolume.hpp"
#include "debug/debug.hpp"
#include "macro.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers. All names are prefixed with `dmcm_` (debug_molecule_complex_
// missmatch) so they cannot collide with helpers from other test files.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal but *fully initialized* Molecule.
 *
 * Both functions under test read `index`, `id`, `myComIndex` and `isEmpty`, and
 * check_complex_index additionally calls Molecule::print() on failure, which
 * touches comCoord and the bookkeeping vectors. Everything those paths read is
 * initialized here so nothing is left dangling.
 *
 * \param[in] index    Index of the molecule inside moleculeList.
 * \param[in] comIndex Value stored in Molecule::myComIndex.
 * \param[in] isEmpty  Whether the molecule is flagged as destroyed/void.
 */
Molecule dmcm_make_molecule(int index, int comIndex, bool isEmpty = false)
{
    Molecule mol;
    mol.index = index;
    mol.id = index; // unique id, used only in error strings
    mol.myComIndex = comIndex;
    mol.complexId = comIndex;
    mol.molTypeIndex = 0;
    mol.mySubVolIndex = 0;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isEmpty = isEmpty;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = Coord { 0.0, 0.0, 0.0 };

    // One interface coincident with the centre of mass keeps geometry trivial.
    Molecule::Iface iface;
    iface.coord = Coord { 0.0, 0.0, 0.0 };
    iface.relIndex = 0;
    iface.index = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // Association bookkeeping lists are read by Molecule::print(); leave empty
    // but well-formed.
    mol.freelist.clear();
    mol.bndlist.clear();
    mol.bndpartner.clear();

    return mol;
}

/*! \brief Build a minimal Complex owning the given member molecule indices. */
Complex dmcm_make_complex(int index, const std::vector<int>& members, bool isEmpty = false)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.isEmpty = isEmpty;
    com.memberList = members;
    com.numEachMol = std::vector<int>(1, static_cast<int>(members.size()));
    com.lastNumberUpdateItrEachMol = std::vector<long long int>(1, 0);
    com.comCoord = Coord { 0.0, 0.0, 0.0 };
    com.mass = static_cast<double>(members.size());
    com.radius = 1.0;
    com.D = Coord { 1.0, 1.0, 1.0 };
    com.Dr = Coord { 0.01, 0.01, 0.01 };
    com.trajStatus = TrajStatus::none;
    return com;
}

/*! \brief Populate the MpiContext fields dereferenced by Molecule::print().
 *
 * Molecule::print() (invoked by check_complex_index just before it throws)
 * dereferences mpiContext.membraneObject and mpiContext.simulVolume and divides
 * by simulVolume.subCellSize.x, so both are given non-degenerate values here.
 */
void dmcm_setup_context(MpiContext& ctx, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, Membrane& membraneObject, SimulVolume& simulVolume)
{
    membraneObject.waterBox.x = 100.0;
    membraneObject.waterBox.y = 100.0;
    membraneObject.waterBox.z = 100.0;

    simulVolume.subCellSize = Coord { 10.0, 10.0, 10.0 };
    simulVolume.numSubCells.x = 10;
    simulVolume.numSubCells.y = 10;
    simulVolume.numSubCells.z = 10;
    simulVolume.numSubCells.tot = 1000;

    ctx.rank = 0;
    ctx.nprocs = 1;
    ctx.xOffset = 0;
    ctx.membraneObject = &membraneObject;
    ctx.simulVolume = &simulVolume;
    ctx.moleculeList = &moleculeList;
    ctx.complexList = &complexList;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: debug_molecule_complex_missmatch() on perfectly consistent data.
// -----------------------------------------------------------------------------
void test_dmcm_consistent_data_is_noop()
{
    std::cerr << "\n[TEST] test_dmcm_consistent_data_is_noop\n"
              << "  Source file:   debug/debug_molecule_complex_missmatch.cpp\n"
              << "  Function:      debug_molecule_complex_missmatch()\n"
              << "  Scenario:      2 complexes, 3 molecules, all cross references\n"
              << "                 (mol.myComIndex <-> com.memberList) consistent.\n"
              << "  Pass criteria: the call returns without throwing and does not\n"
              << "                 alter any molecule or complex field.\n";

    // Molecule 0 is treated as the implicit lipid slot by the checker
    // (`if (!mol.index) continue;`), molecules 1 and 2 are real members.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dmcm_make_molecule(0, 0)); // "implicit lipid" slot
    moleculeList.push_back(dmcm_make_molecule(1, 0));
    moleculeList.push_back(dmcm_make_molecule(2, 1));

    std::vector<Complex> complexList;
    complexList.push_back(dmcm_make_complex(0, std::vector<int> { 0, 1 }));
    complexList.push_back(dmcm_make_complex(1, std::vector<int> { 2 }));

    Membrane membraneObject;
    SimulVolume simulVolume;
    MpiContext mpiContext;
    dmcm_setup_context(mpiContext, moleculeList, complexList, membraneObject, simulVolume);

    std::cerr << "  Calling debug_molecule_complex_missmatch(\"consistent data\")...\n";
    EXPECT_NO_THROW(debug_molecule_complex_missmatch(
        mpiContext, moleculeList, complexList, std::string(" consistent data")))
        << "Consistent data must never raise an exception";

    // Nothing must have been added, removed or renumbered.
    EXPECT_EQ(moleculeList.size(), 3u) << "moleculeList size must be preserved";
    EXPECT_EQ(complexList.size(), 2u) << "complexList size must be preserved";
    EXPECT_EQ(moleculeList[1].myComIndex, 0) << "molecule 1 must still point at complex 0";
    EXPECT_EQ(moleculeList[2].myComIndex, 1) << "molecule 2 must still point at complex 1";
    EXPECT_EQ(complexList[0].memberList.size(), 2u) << "complex 0 member count must be preserved";
    EXPECT_EQ(complexList[1].memberList.size(), 1u) << "complex 1 member count must be preserved";
    EXPECT_EQ(complexList[1].memberList[0], 2) << "complex 1 must still own molecule 2";

    std::cerr << "  Data unchanged after the call, as expected.\n";
}

// -----------------------------------------------------------------------------
// Test 2: debug_molecule_complex_missmatch() with deliberately broken data.
//
// NOTE: with DEBUG enabled this input would make the function call error(),
// which terminates the process. macro.hpp defines DEBUG as false, so the body is
// compiled out and the function must behave as a no-op. The inconsistent input
// is therefore only exercised when DEBUG is false; the `if (DEBUG)` guard below
// keeps this test harmless should the macro ever be flipped on.
// -----------------------------------------------------------------------------
void test_dmcm_inconsistent_data_is_noop_when_debug_disabled()
{
    std::cerr << "\n[TEST] test_dmcm_inconsistent_data_is_noop_when_debug_disabled\n"
              << "  Source file:   debug/debug_molecule_complex_missmatch.cpp\n"
              << "  Function:      debug_molecule_complex_missmatch()\n"
              << "  Scenario:      molecule points at a complex that does not list it,\n"
              << "                 and a complex lists a molecule index that is out of\n"
              << "                 range for moleculeList.\n"
              << "  Pass criteria: with DEBUG == false the call is a no-op: it returns,\n"
              << "                 does not throw, and leaves the broken data as-is.\n";

    if (DEBUG) {
        std::cerr << "  DEBUG is enabled in macro.hpp -> the checker would call error()\n"
                  << "  and terminate the process. Skipping the inconsistent-data case.\n";
        SUCCEED() << "Skipped because DEBUG is enabled";
        return;
    }

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dmcm_make_molecule(0, 0)); // implicit-lipid slot
    // Molecule 1 claims complex 0, but complex 0's memberList will not contain it.
    moleculeList.push_back(dmcm_make_molecule(1, 0));

    std::vector<Complex> complexList;
    // Complex 0 lists molecule index 99, which is far outside moleculeList.
    complexList.push_back(dmcm_make_complex(0, std::vector<int> { 0, 99 }));

    Membrane membraneObject;
    SimulVolume simulVolume;
    MpiContext mpiContext;
    dmcm_setup_context(mpiContext, moleculeList, complexList, membraneObject, simulVolume);

    std::cerr << "  Calling debug_molecule_complex_missmatch(\"broken data\")...\n";
    EXPECT_NO_THROW(debug_molecule_complex_missmatch(
        mpiContext, moleculeList, complexList, std::string(" deliberately broken data")))
        << "With DEBUG disabled the checker must be an inert no-op";

    // The function is a checker, never a fixer: the broken values must survive.
    EXPECT_EQ(moleculeList.size(), 2u) << "moleculeList must not be resized";
    EXPECT_EQ(complexList.size(), 1u) << "complexList must not be resized";
    EXPECT_EQ(moleculeList[1].myComIndex, 0) << "molecule 1's (wrong) complex index must be untouched";
    ASSERT_EQ(complexList[0].memberList.size(), 2u) << "memberList must not be edited";
    EXPECT_EQ(complexList[0].memberList[1], 99) << "the bogus member index must be untouched";

    std::cerr << "  Broken data preserved verbatim, as expected for a no-op.\n";
}

// -----------------------------------------------------------------------------
// Test 3: debug_molecule_complex_missmatch() with empty containers.
// -----------------------------------------------------------------------------
void test_dmcm_empty_containers()
{
    std::cerr << "\n[TEST] test_dmcm_empty_containers\n"
              << "  Source file:   debug/debug_molecule_complex_missmatch.cpp\n"
              << "  Function:      debug_molecule_complex_missmatch()\n"
              << "  Scenario:      empty moleculeList and empty complexList.\n"
              << "  Pass criteria: returns cleanly, containers stay empty.\n";

    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;

    Membrane membraneObject;
    SimulVolume simulVolume;
    MpiContext mpiContext;
    dmcm_setup_context(mpiContext, moleculeList, complexList, membraneObject, simulVolume);

    EXPECT_NO_THROW(debug_molecule_complex_missmatch(
        mpiContext, moleculeList, complexList, std::string(" empty containers")))
        << "Empty containers must be handled gracefully";

    EXPECT_TRUE(moleculeList.empty()) << "moleculeList must remain empty";
    EXPECT_TRUE(complexList.empty()) << "complexList must remain empty";

    std::cerr << "  Empty containers handled without incident.\n";
}

// -----------------------------------------------------------------------------
// Test 4: check_complex_index() accepts fully valid indices.
// -----------------------------------------------------------------------------
void test_dmcm_check_complex_index_valid()
{
    std::cerr << "\n[TEST] test_dmcm_check_complex_index_valid\n"
              << "  Source file:   debug/debug_molecule_complex_missmatch.cpp\n"
              << "  Function:      check_complex_index()\n"
              << "  Scenario:      every live molecule has 0 <= myComIndex < complexList.size().\n"
              << "  Pass criteria: no exception is thrown and no molecule is modified.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dmcm_make_molecule(0, 0)); // lowest valid index
    moleculeList.push_back(dmcm_make_molecule(1, 1));
    moleculeList.push_back(dmcm_make_molecule(2, 2)); // highest valid index

    std::vector<Complex> complexList;
    complexList.push_back(dmcm_make_complex(0, std::vector<int> { 0 }));
    complexList.push_back(dmcm_make_complex(1, std::vector<int> { 1 }));
    complexList.push_back(dmcm_make_complex(2, std::vector<int> { 2 }));

    Membrane membraneObject;
    SimulVolume simulVolume;
    MpiContext mpiContext;
    dmcm_setup_context(mpiContext, moleculeList, complexList, membraneObject, simulVolume);

    std::cerr << "  Calling check_complex_index() with 3 valid molecules...\n";
    EXPECT_NO_THROW(check_complex_index(mpiContext, moleculeList, complexList))
        << "Valid complex indices must not raise std::runtime_error";

    // The routine is read-only for valid input.
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "molecule 0's complex index must be unchanged";
    EXPECT_EQ(moleculeList[1].myComIndex, 1) << "molecule 1's complex index must be unchanged";
    EXPECT_EQ(moleculeList[2].myComIndex, 2) << "molecule 2's complex index must be unchanged";
    EXPECT_EQ(moleculeList.size(), 3u) << "moleculeList must not be resized";

    std::cerr << "  All indices accepted; molecules untouched.\n";
}

// -----------------------------------------------------------------------------
// Test 5: check_complex_index() ignores molecules flagged isEmpty.
// -----------------------------------------------------------------------------
void test_dmcm_check_complex_index_skips_empty_molecules()
{
    std::cerr << "\n[TEST] test_dmcm_check_complex_index_skips_empty_molecules\n"
              << "  Source file:   debug/debug_molecule_complex_missmatch.cpp\n"
              << "  Function:      check_complex_index()\n"
              << "  Scenario:      a destroyed molecule (isEmpty == true) carries the\n"
              << "                 sentinel myComIndex == -1, which is out of range.\n"
              << "  Pass criteria: `if (mol.isEmpty) continue;` makes the function skip\n"
              << "                 it, so no exception is thrown.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dmcm_make_molecule(0, 0)); // live, valid
    // Destroyed molecule: Molecule::destroy() sets myComIndex to -1 and isEmpty
    // to true. That index is illegal but must be tolerated.
    moleculeList.push_back(dmcm_make_molecule(1, -1, /*isEmpty=*/true));

    std::vector<Complex> complexList;
    complexList.push_back(dmcm_make_complex(0, std::vector<int> { 0 }));

    Membrane membraneObject;
    SimulVolume simulVolume;
    MpiContext mpiContext;
    dmcm_setup_context(mpiContext, moleculeList, complexList, membraneObject, simulVolume);

    std::cerr << "  Calling check_complex_index() with one destroyed molecule...\n";
    EXPECT_NO_THROW(check_complex_index(mpiContext, moleculeList, complexList))
        << "Molecules flagged isEmpty must be skipped, not reported";

    EXPECT_TRUE(moleculeList[1].isEmpty) << "the destroyed molecule must stay flagged empty";
    EXPECT_EQ(moleculeList[1].myComIndex, -1) << "the sentinel index must be left alone";

    std::cerr << "  Destroyed molecule correctly skipped.\n";
}

// -----------------------------------------------------------------------------
// Test 6: check_complex_index() on an empty moleculeList.
// -----------------------------------------------------------------------------
void test_dmcm_check_complex_index_empty_list()
{
    std::cerr << "\n[TEST] test_dmcm_check_complex_index_empty_list\n"
              << "  Source file:   debug/debug_molecule_complex_missmatch.cpp\n"
              << "  Function:      check_complex_index()\n"
              << "  Scenario:      empty moleculeList (loop body never executes).\n"
              << "  Pass criteria: returns without throwing.\n";

    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    complexList.push_back(dmcm_make_complex(0, std::vector<int> {}));

    Membrane membraneObject;
    SimulVolume simulVolume;
    MpiContext mpiContext;
    dmcm_setup_context(mpiContext, moleculeList, complexList, membraneObject, simulVolume);

    EXPECT_NO_THROW(check_complex_index(mpiContext, moleculeList, complexList))
        << "An empty moleculeList must be a trivial success";

    std::cerr << "  Empty moleculeList handled without incident.\n";
}

// -----------------------------------------------------------------------------
// Test 7: check_complex_index() throws for a negative complex index.
//
// The failure path prints the offending molecule via Molecule::print() before
// throwing; the MpiContext supplied by dmcm_setup_context() carries the Membrane
// and SimulVolume pointers that printing dereferences.
// -----------------------------------------------------------------------------
void test_dmcm_check_complex_index_negative_throws()
{
    std::cerr << "\n[TEST] test_dmcm_check_complex_index_negative_throws\n"
              << "  Source file:   debug/debug_molecule_complex_missmatch.cpp\n"
              << "  Function:      check_complex_index()\n"
              << "  Scenario:      a live molecule carries myComIndex == -1.\n"
              << "  Pass criteria: std::runtime_error is thrown (after the molecule is\n"
              << "                 dumped to stdout).\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dmcm_make_molecule(0, 0)); // valid
    moleculeList.push_back(dmcm_make_molecule(1, -1)); // live but index < 0

    std::vector<Complex> complexList;
    complexList.push_back(dmcm_make_complex(0, std::vector<int> { 0 }));

    Membrane membraneObject;
    SimulVolume simulVolume;
    MpiContext mpiContext;
    dmcm_setup_context(mpiContext, moleculeList, complexList, membraneObject, simulVolume);

    std::cerr << "  Calling check_complex_index(); a molecule dump on stdout is expected...\n";
    EXPECT_THROW(check_complex_index(mpiContext, moleculeList, complexList), std::runtime_error)
        << "A negative myComIndex on a live molecule must raise std::runtime_error";

    std::cerr << "  std::runtime_error raised as expected.\n";
}

// -----------------------------------------------------------------------------
// Test 8: check_complex_index() throws when myComIndex == complexList.size().
//         This pins down the exclusive upper bound of the accepted range.
// -----------------------------------------------------------------------------
void test_dmcm_check_complex_index_out_of_range_throws()
{
    std::cerr << "\n[TEST] test_dmcm_check_complex_index_out_of_range_throws\n"
              << "  Source file:   debug/debug_molecule_complex_missmatch.cpp\n"
              << "  Function:      check_complex_index()\n"
              << "  Scenario:      complexList holds 2 complexes and a live molecule\n"
              << "                 carries myComIndex == 2 (one past the end).\n"
              << "  Pass criteria: std::runtime_error is thrown, confirming the upper\n"
              << "                 bound is exclusive.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dmcm_make_molecule(0, 0)); // valid
    moleculeList.push_back(dmcm_make_molecule(1, 1)); // valid (last legal index)
    moleculeList.push_back(dmcm_make_molecule(2, 2)); // one past the end -> invalid

    std::vector<Complex> complexList;
    complexList.push_back(dmcm_make_complex(0, std::vector<int> { 0 }));
    complexList.push_back(dmcm_make_complex(1, std::vector<int> { 1 }));

    Membrane membraneObject;
    SimulVolume simulVolume;
    MpiContext mpiContext;
    dmcm_setup_context(mpiContext, moleculeList, complexList, membraneObject, simulVolume);

    std::cerr << "  complexList.size() == " << complexList.size()
              << "; offending molecule myComIndex == " << moleculeList[2].myComIndex << '\n';
    EXPECT_THROW(check_complex_index(mpiContext, moleculeList, complexList), std::runtime_error)
        << "myComIndex == complexList.size() must be rejected";

    std::cerr << "  std::runtime_error raised as expected.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper is invoked from its own TEST so
// a failure in one case does not prevent the remaining cases from running.
// -----------------------------------------------------------------------------
TEST(DebugMoleculeComplexMissmatch, ConsistentDataIsNoop) { test_dmcm_consistent_data_is_noop(); }
TEST(DebugMoleculeComplexMissmatch, InconsistentDataIsNoopWhenDebugDisabled)
{
    test_dmcm_inconsistent_data_is_noop_when_debug_disabled();
}
TEST(DebugMoleculeComplexMissmatch, EmptyContainers) { test_dmcm_empty_containers(); }
TEST(CheckComplexIndex, ValidIndices) { test_dmcm_check_complex_index_valid(); }
TEST(CheckComplexIndex, SkipsEmptyMolecules) { test_dmcm_check_complex_index_skips_empty_molecules(); }
TEST(CheckComplexIndex, EmptyMoleculeList) { test_dmcm_check_complex_index_empty_list(); }
TEST(CheckComplexIndex, NegativeIndexThrows) { test_dmcm_check_complex_index_negative_throws(); }
TEST(CheckComplexIndex, OutOfRangeIndexThrows) { test_dmcm_check_complex_index_out_of_range_throws(); }