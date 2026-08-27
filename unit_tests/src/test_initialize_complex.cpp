/*! \file test_initialize_complex.cpp
 *
 * ### Unit test for ../src/system_setup/initialize_complex.cpp
 *
 * The single function under test is:
 *
 *     Complex initialize_complex(const Molecule& mol, const MolTemplate& molTemp)
 *
 * It builds a brand new Complex that wraps exactly one Molecule.  The
 * behaviour that must hold, as read directly from the implementation, is:
 *
 *   - comCoord, D, Dr, radius and mass are copied from the Molecule / MolTemplate
 *   - memberList contains exactly one element: mol.index
 *   - isEmpty is explicitly set to false
 *   - index    == the value of the static Complex::numberOfComplexes *before* the call
 *   - id       == the value of the static Complex::maxID           *before* the call
 *   - Complex::numberOfComplexes and Complex::maxID are both incremented by one
 *   - numEachMol is sized to MolTemplate::numMolTypes, zeroed, then the slot
 *     belonging to molTemp.molTypeIndex is incremented to 1
 *   - lastNumberUpdateItrEachMol is resized to MolTemplate::numMolTypes
 *
 * IMPORTANT: the function indexes numEachMol with molTemp.molTypeIndex without
 * any bounds check, so every test below sets MolTemplate::numMolTypes large
 * enough before calling.  The static counters are saved and restored by each
 * test so that this translation unit does not perturb the rest of the suite.
 */

#include "system_setup/system_setup.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Small helpers, uniquely prefixed so they cannot clash with other test files.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully initialized Molecule for use as the seed of a Complex.
 *
 * initialize_complex() only reads mol.comCoord and mol.index, but we fill in
 * the rest of the commonly used fields so the object is never half-initialized.
 */
Molecule ic_make_molecule(int index, int molTypeIndex, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = molTypeIndex;
    mol.comCoord = com;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.myComIndex = -1;
    mol.interfaceList.clear();
    return mol;
}

/*! \brief Build a MolTemplate carrying distinctive diffusion / size values.
 *
 * Distinctive numbers make it obvious in the assertions which member of the
 * template each Complex field was copied from.
 */
MolTemplate ic_make_moltemplate(int molTypeIndex, const std::string& name)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.mass = 2.5;
    temp.radius = 3.75;
    temp.D = Coord(11.0, 12.0, 13.0);
    temp.Dr = Coord(0.11, 0.12, 0.13);
    temp.copies = 1;
    temp.interfaceList.clear();
    return temp;
}

/*! \brief RAII-style snapshot/restore of every static this function mutates.
 *
 * All generated tests share one binary, so leaving Complex::numberOfComplexes,
 * Complex::maxID or MolTemplate::numMolTypes modified would leak into other
 * translation units.  Constructing this struct snapshots them; destroying it
 * puts them back exactly as they were.
 */
struct IcStaticsGuard {
    int savedNumberOfComplexes;
    int savedMaxID;
    unsigned savedNumMolTypes;

    IcStaticsGuard()
        : savedNumberOfComplexes(Complex::numberOfComplexes)
        , savedMaxID(Complex::maxID)
        , savedNumMolTypes(MolTemplate::numMolTypes)
    {
    }

    ~IcStaticsGuard()
    {
        Complex::numberOfComplexes = savedNumberOfComplexes;
        Complex::maxID = savedMaxID;
        MolTemplate::numMolTypes = savedNumMolTypes;
    }
};

} // namespace

// -----------------------------------------------------------------------------
// Test 1: every value-carrying field is copied straight out of the inputs.
// -----------------------------------------------------------------------------
void test_ic_copies_geometry_and_diffusion()
{
    std::cerr << "\n[TEST] test_ic_copies_geometry_and_diffusion\n"
              << "  Source file:   src/system_setup/initialize_complex.cpp\n"
              << "  Function:      initialize_complex()\n"
              << "  Scenario:      one molecule at (1.5, -2.5, 3.25) built from a\n"
              << "                 template with D=(11,12,13), Dr=(0.11,0.12,0.13),\n"
              << "                 radius=3.75, mass=2.5.\n"
              << "  Pass criteria: comCoord/D/Dr/radius/mass are copied verbatim,\n"
              << "                 memberList == {mol.index}, isEmpty == false.\n";

    IcStaticsGuard guard; // restores the statics when this test returns

    // Two molecule types exist in this scenario; our molecule is type 0.
    MolTemplate::numMolTypes = 2;

    const Coord com { 1.5, -2.5, 3.25 };
    Molecule mol = ic_make_molecule(/*index=*/7, /*molTypeIndex=*/0, com);
    MolTemplate temp = ic_make_moltemplate(/*molTypeIndex=*/0, "A");

    std::cerr << "  Calling initialize_complex()...\n";
    Complex newCom = initialize_complex(mol, temp);

    // --- center of mass is taken from the Molecule, not the template ---------
    EXPECT_DOUBLE_EQ(newCom.comCoord.x, com.x) << "comCoord.x must come from the Molecule";
    EXPECT_DOUBLE_EQ(newCom.comCoord.y, com.y) << "comCoord.y must come from the Molecule";
    EXPECT_DOUBLE_EQ(newCom.comCoord.z, com.z) << "comCoord.z must come from the Molecule";

    // --- translational diffusion comes from the template ---------------------
    EXPECT_DOUBLE_EQ(newCom.D.x, temp.D.x) << "D.x must be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(newCom.D.y, temp.D.y) << "D.y must be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(newCom.D.z, temp.D.z) << "D.z must be copied from the MolTemplate";

    // --- rotational diffusion comes from the template ------------------------
    EXPECT_DOUBLE_EQ(newCom.Dr.x, temp.Dr.x) << "Dr.x must be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(newCom.Dr.y, temp.Dr.y) << "Dr.y must be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(newCom.Dr.z, temp.Dr.z) << "Dr.z must be copied from the MolTemplate";

    // --- scalar properties ---------------------------------------------------
    EXPECT_DOUBLE_EQ(newCom.radius, temp.radius) << "radius must be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(newCom.mass, temp.mass)
        << "mass must be copied from the MolTemplate (NOT from the Molecule)";

    // The Molecule mass is deliberately different (1.0 vs 2.5); confirm that the
    // template value, not the molecule value, is what landed in the Complex.
    EXPECT_NE(newCom.mass, mol.mass)
        << "the implementation uses molTemp.mass, so it should differ from mol.mass here";

    // --- membership ----------------------------------------------------------
    ASSERT_EQ(newCom.memberList.size(), static_cast<size_t>(1))
        << "a freshly created Complex holds exactly one member Molecule";
    EXPECT_EQ(newCom.memberList[0], mol.index)
        << "the single member must be the index of the seeding Molecule";

    // --- the complex is live, not a recycled empty slot ----------------------
    EXPECT_FALSE(newCom.isEmpty) << "initialize_complex explicitly sets isEmpty = false";

    std::cerr << "  Resulting Complex: index=" << newCom.index << ", id=" << newCom.id
              << ", radius=" << newCom.radius << ", mass=" << newCom.mass << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: index/id are taken from the statics and both counters are bumped.
// -----------------------------------------------------------------------------
void test_ic_assigns_index_and_id_from_statics()
{
    std::cerr << "\n[TEST] test_ic_assigns_index_and_id_from_statics\n"
              << "  Source file:   src/system_setup/initialize_complex.cpp\n"
              << "  Function:      initialize_complex()\n"
              << "  Scenario:      the static counters are seeded to known values\n"
              << "                 before a single call.\n"
              << "  Pass criteria: index == pre-call numberOfComplexes,\n"
              << "                 id    == pre-call maxID, and both statics grow by 1.\n";

    IcStaticsGuard guard;

    MolTemplate::numMolTypes = 1;

    // Seed the statics with easily recognizable values.
    Complex::numberOfComplexes = 42;
    Complex::maxID = 100;

    const int expectedIndex = Complex::numberOfComplexes;
    const int expectedId = Complex::maxID;

    Molecule mol = ic_make_molecule(/*index=*/0, /*molTypeIndex=*/0, Coord { 0.0, 0.0, 0.0 });
    MolTemplate temp = ic_make_moltemplate(/*molTypeIndex=*/0, "A");

    std::cerr << "  Before call: numberOfComplexes=" << Complex::numberOfComplexes
              << ", maxID=" << Complex::maxID << '\n';

    Complex newCom = initialize_complex(mol, temp);

    std::cerr << "  After  call: numberOfComplexes=" << Complex::numberOfComplexes
              << ", maxID=" << Complex::maxID << " (new index=" << newCom.index
              << ", new id=" << newCom.id << ")\n";

    // index is read before the counter is incremented.
    EXPECT_EQ(newCom.index, expectedIndex)
        << "index should equal Complex::numberOfComplexes as it was before the call";

    // id uses a post-increment, so it also gets the pre-call value.
    EXPECT_EQ(newCom.id, expectedId)
        << "id should equal Complex::maxID as it was before the call (post-increment)";

    // Both global counters must advance by exactly one.
    EXPECT_EQ(Complex::numberOfComplexes, expectedIndex + 1)
        << "Complex::numberOfComplexes must be incremented once per created Complex";
    EXPECT_EQ(Complex::maxID, expectedId + 1)
        << "Complex::maxID must be incremented once per created Complex";
}

// -----------------------------------------------------------------------------
// Test 3: successive calls hand out consecutive, unique indices and ids.
// -----------------------------------------------------------------------------
void test_ic_successive_calls_are_unique()
{
    std::cerr << "\n[TEST] test_ic_successive_calls_are_unique\n"
              << "  Source file:   src/system_setup/initialize_complex.cpp\n"
              << "  Function:      initialize_complex()\n"
              << "  Scenario:      three complexes are created back to back.\n"
              << "  Pass criteria: index and id increase by exactly one each time\n"
              << "                 and never repeat.\n";

    IcStaticsGuard guard;

    MolTemplate::numMolTypes = 2;

    Complex::numberOfComplexes = 0;
    Complex::maxID = 0;

    MolTemplate tempA = ic_make_moltemplate(/*molTypeIndex=*/0, "A");
    MolTemplate tempB = ic_make_moltemplate(/*molTypeIndex=*/1, "B");

    Molecule mol0 = ic_make_molecule(0, 0, Coord { 1.0, 0.0, 0.0 });
    Molecule mol1 = ic_make_molecule(1, 1, Coord { 0.0, 1.0, 0.0 });
    Molecule mol2 = ic_make_molecule(2, 0, Coord { 0.0, 0.0, 1.0 });

    Complex c0 = initialize_complex(mol0, tempA);
    Complex c1 = initialize_complex(mol1, tempB);
    Complex c2 = initialize_complex(mol2, tempA);

    std::cerr << "  Created indices: " << c0.index << ", " << c1.index << ", " << c2.index << '\n';
    std::cerr << "  Created ids:     " << c0.id << ", " << c1.id << ", " << c2.id << '\n';

    // Consecutive indices, starting from the seeded zero.
    EXPECT_EQ(c0.index, 0) << "first complex should receive index 0";
    EXPECT_EQ(c1.index, 1) << "second complex should receive index 1";
    EXPECT_EQ(c2.index, 2) << "third complex should receive index 2";

    // Consecutive ids as well.
    EXPECT_EQ(c0.id, 0) << "first complex should receive id 0";
    EXPECT_EQ(c1.id, 1) << "second complex should receive id 1";
    EXPECT_EQ(c2.id, 2) << "third complex should receive id 2";

    // The global count must reflect all three creations.
    EXPECT_EQ(Complex::numberOfComplexes, 3)
        << "three calls must leave Complex::numberOfComplexes at 3";
    EXPECT_EQ(Complex::maxID, 3) << "three calls must leave Complex::maxID at 3";

    // Each complex still owns exactly its own seeding molecule.
    ASSERT_EQ(c0.memberList.size(), static_cast<size_t>(1));
    ASSERT_EQ(c1.memberList.size(), static_cast<size_t>(1));
    ASSERT_EQ(c2.memberList.size(), static_cast<size_t>(1));
    EXPECT_EQ(c0.memberList[0], mol0.index) << "complex 0 should own molecule 0";
    EXPECT_EQ(c1.memberList[0], mol1.index) << "complex 1 should own molecule 1";
    EXPECT_EQ(c2.memberList[0], mol2.index) << "complex 2 should own molecule 2";
}

// -----------------------------------------------------------------------------
// Test 4: numEachMol bookkeeping -- correct size, correct single tally.
// -----------------------------------------------------------------------------
void test_ic_numEachMol_bookkeeping()
{
    std::cerr << "\n[TEST] test_ic_numEachMol_bookkeeping\n"
              << "  Source file:   src/system_setup/initialize_complex.cpp\n"
              << "  Function:      initialize_complex()\n"
              << "  Scenario:      4 molecule types exist; the seed molecule is type 2.\n"
              << "  Pass criteria: numEachMol has 4 entries, all zero except slot 2\n"
              << "                 which holds 1.\n";

    IcStaticsGuard guard;

    // Four distinct molecule types in this hypothetical system.
    MolTemplate::numMolTypes = 4;

    const int targetType = 2;
    Molecule mol = ic_make_molecule(/*index=*/5, targetType, Coord { 0.0, 0.0, 0.0 });
    MolTemplate temp = ic_make_moltemplate(targetType, "C");

    Complex newCom = initialize_complex(mol, temp);

    // The vector is constructed as std::vector<int>(MolTemplate::numMolTypes),
    // i.e. value-initialized to zero, then a single slot is incremented.
    ASSERT_EQ(newCom.numEachMol.size(), static_cast<size_t>(MolTemplate::numMolTypes))
        << "numEachMol must be sized to MolTemplate::numMolTypes";

    std::cerr << "  numEachMol =";
    for (int count : newCom.numEachMol)
        std::cerr << ' ' << count;
    std::cerr << '\n';

    for (int typeItr = 0; typeItr < static_cast<int>(newCom.numEachMol.size()); ++typeItr) {
        if (typeItr == targetType) {
            EXPECT_EQ(newCom.numEachMol[typeItr], 1)
                << "slot " << typeItr << " (the seed molecule's type) must be exactly 1";
        } else {
            EXPECT_EQ(newCom.numEachMol[typeItr], 0)
                << "slot " << typeItr << " holds no molecules and must be 0";
        }
    }

    // The per-molecule-type "last update iteration" tracker is resized (and so
    // value-initialized to zero) alongside numEachMol.
    ASSERT_EQ(newCom.lastNumberUpdateItrEachMol.size(),
        static_cast<size_t>(MolTemplate::numMolTypes))
        << "lastNumberUpdateItrEachMol must be resized to MolTemplate::numMolTypes";
    for (size_t itr = 0; itr < newCom.lastNumberUpdateItrEachMol.size(); ++itr) {
        EXPECT_EQ(newCom.lastNumberUpdateItrEachMol[itr], 0LL)
            << "resize() value-initializes lastNumberUpdateItrEachMol entry " << itr << " to 0";
    }
}

// -----------------------------------------------------------------------------
// Test 5: fields the function does NOT touch keep their class defaults.
// -----------------------------------------------------------------------------
void test_ic_untouched_fields_keep_defaults()
{
    std::cerr << "\n[TEST] test_ic_untouched_fields_keep_defaults\n"
              << "  Source file:   src/system_setup/initialize_complex.cpp\n"
              << "  Function:      initialize_complex()\n"
              << "  Scenario:      a lipid-flagged molecule is used as the seed.\n"
              << "  Pass criteria: initialize_complex does not inspect isLipid, so\n"
              << "                 OnSurface/onFiber/linksToSurface/ncross/trajStatus\n"
              << "                 all keep their in-class default values.\n";

    IcStaticsGuard guard;

    MolTemplate::numMolTypes = 1;

    // Deliberately flag the molecule as a lipid.  Complex's own constructor
    // would set OnSurface for such a molecule, but initialize_complex() never
    // looks at that flag -- this test pins that documented difference down.
    Molecule mol = ic_make_molecule(/*index=*/0, /*molTypeIndex=*/0, Coord { 0.0, 0.0, -5.0 });
    mol.isLipid = true;
    mol.isPromoter = true;

    MolTemplate temp = ic_make_moltemplate(/*molTypeIndex=*/0, "Lipid");
    temp.isLipid = true;
    temp.isPromoter = true;

    Complex newCom = initialize_complex(mol, temp);

    std::cerr << "  OnSurface=" << std::boolalpha << newCom.OnSurface
              << ", onFiber=" << newCom.onFiber
              << ", linksToSurface=" << newCom.linksToSurface
              << ", ncross=" << newCom.ncross << '\n';

    EXPECT_FALSE(newCom.OnSurface)
        << "initialize_complex never sets OnSurface, so it keeps its default false";
    EXPECT_FALSE(newCom.onFiber)
        << "initialize_complex never sets onFiber, so it keeps its default false";
    EXPECT_FALSE(newCom.tmpOnSurface)
        << "tmpOnSurface is untouched and keeps its default false";
    EXPECT_EQ(newCom.linksToSurface, 0) << "linksToSurface is untouched and defaults to 0";
    EXPECT_EQ(newCom.iLipidIndex, 0) << "iLipidIndex is untouched and defaults to 0";
    EXPECT_EQ(newCom.ncross, 0) << "ncross is untouched and defaults to 0";
    EXPECT_EQ(newCom.trajStatus, TrajStatus::none)
        << "trajStatus is untouched and defaults to TrajStatus::none";

    // The propagation vectors are default-constructed to the origin.
    EXPECT_DOUBLE_EQ(newCom.trajTrans.x, 0.0) << "trajTrans starts at the origin";
    EXPECT_DOUBLE_EQ(newCom.trajTrans.y, 0.0) << "trajTrans starts at the origin";
    EXPECT_DOUBLE_EQ(newCom.trajTrans.z, 0.0) << "trajTrans starts at the origin";
    EXPECT_DOUBLE_EQ(newCom.trajRot.x, 0.0) << "trajRot starts at the origin";
    EXPECT_DOUBLE_EQ(newCom.trajRot.y, 0.0) << "trajRot starts at the origin";
    EXPECT_DOUBLE_EQ(newCom.trajRot.z, 0.0) << "trajRot starts at the origin";
}

// -----------------------------------------------------------------------------
// Test 6: the returned Complex is an independent copy of its inputs.
// -----------------------------------------------------------------------------
void test_ic_returned_complex_is_independent()
{
    std::cerr << "\n[TEST] test_ic_returned_complex_is_independent\n"
              << "  Source file:   src/system_setup/initialize_complex.cpp\n"
              << "  Function:      initialize_complex()\n"
              << "  Scenario:      after construction, the seeding Molecule and\n"
              << "                 MolTemplate are mutated.\n"
              << "  Pass criteria: the Complex retains the values captured at the\n"
              << "                 moment of the call (everything is copied by value).\n";

    IcStaticsGuard guard;

    MolTemplate::numMolTypes = 1;

    Molecule mol = ic_make_molecule(/*index=*/3, /*molTypeIndex=*/0, Coord { 9.0, 8.0, 7.0 });
    MolTemplate temp = ic_make_moltemplate(/*molTypeIndex=*/0, "A");

    Complex newCom = initialize_complex(mol, temp);

    // Snapshot what the Complex received.
    const double capturedX = newCom.comCoord.x;
    const double capturedRadius = newCom.radius;
    const double capturedDx = newCom.D.x;
    const int capturedMember = newCom.memberList[0];

    std::cerr << "  Mutating the source Molecule and MolTemplate after the call...\n";
    mol.comCoord = Coord { -1.0, -1.0, -1.0 };
    mol.index = 999;
    temp.radius = 1000.0;
    temp.D = Coord(-1.0, -1.0, -1.0);

    // Nothing about the already-built Complex may change.
    EXPECT_DOUBLE_EQ(newCom.comCoord.x, capturedX)
        << "comCoord was copied by value and must not follow the Molecule";
    EXPECT_DOUBLE_EQ(newCom.comCoord.x, 9.0) << "comCoord.x should still be the original 9.0";
    EXPECT_DOUBLE_EQ(newCom.radius, capturedRadius)
        << "radius was copied by value and must not follow the MolTemplate";
    EXPECT_DOUBLE_EQ(newCom.D.x, capturedDx)
        << "D was copied by value and must not follow the MolTemplate";
    EXPECT_EQ(newCom.memberList[0], capturedMember)
        << "memberList holds a copy of the index, not a live reference";
    EXPECT_EQ(newCom.memberList[0], 3) << "memberList should still hold the original index 3";

    std::cerr << "  Complex survived source mutation with comCoord.x=" << newCom.comCoord.x
              << ", radius=" << newCom.radius << ", member=" << newCom.memberList[0] << '\n';
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so that a
// failure in one does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(InitializeComplexTest, CopiesGeometryAndDiffusion) { test_ic_copies_geometry_and_diffusion(); }
TEST(InitializeComplexTest, AssignsIndexAndIdFromStatics) { test_ic_assigns_index_and_id_from_statics(); }
TEST(InitializeComplexTest, SuccessiveCallsAreUnique) { test_ic_successive_calls_are_unique(); }
TEST(InitializeComplexTest, NumEachMolBookkeeping) { test_ic_numEachMol_bookkeeping(); }
TEST(InitializeComplexTest, UntouchedFieldsKeepDefaults) { test_ic_untouched_fields_keep_defaults(); }
TEST(InitializeComplexTest, ReturnedComplexIsIndependent) { test_ic_returned_complex_is_independent(); }