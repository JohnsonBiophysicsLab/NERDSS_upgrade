// =============================================================================
// Unit tests for ../src/classes/class_Molecule_Complex.cpp
//
// This test file exercises the free functions, operators, and member
// functions defined in class_Molecule_Complex.cpp.  Because this source file
// implements the behaviour of the Molecule and Complex classes (plus a couple
// of free helper functions and operator overloads), the tests below build up
// small Molecule / Complex objects and verify the observable behaviour of the
// functions that do not require the full simulation environment.
//
// The tests use gtest and prefer non-fatal EXPECT_* assertions so that every
// test runs to completion, even if an individual check fails.  Verbose output
// is written to stderr describing which function (in which file) is being
// tested and what each test is doing.
//
// NOTE: There is intentionally NO main() here; this file is meant to be linked
//       into a larger test suite.
// =============================================================================

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vector.hpp"

#include <gtest/gtest.h>

#include <array>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Forward declaration of the free helper function implemented in the source
// file under test.  It is not declared in a header we necessarily include, so
// we declare it here to be able to call it directly.
// -----------------------------------------------------------------------------
bool skipLine(std::string line);

// -----------------------------------------------------------------------------
// Test the free function skipLine().
//
// skipLine() returns true when a line should be skipped, i.e. when it is
// empty or starts with a '#' comment character.  Otherwise it returns false.
// -----------------------------------------------------------------------------
void test_mc_skipLine() {
  std::cerr << "\n[TEST] class_Molecule_Complex.cpp :: skipLine()\n";
  std::cerr << "       Verifying that empty lines and comment lines are "
               "skipped, while normal lines are not.\n";

  // An empty string should be skipped.
  std::cerr << "       - Checking empty string is skipped.\n";
  EXPECT_TRUE(skipLine("")) << "Empty string should be skipped.";

  // A string starting with '#' is a comment and should be skipped.
  std::cerr << "       - Checking comment line (leading '#') is skipped.\n";
  EXPECT_TRUE(skipLine("# this is a comment"))
      << "Comment line should be skipped.";

  // A normal, non-empty line without a leading '#' should NOT be skipped.
  std::cerr << "       - Checking normal content line is NOT skipped.\n";
  EXPECT_FALSE(skipLine("molecule A"))
      << "Normal content line should not be skipped.";

  // A line with a '#' NOT at the start should not be skipped.
  std::cerr << "       - Checking line with trailing '#' is NOT skipped.\n";
  EXPECT_FALSE(skipLine("value 3 # trailing comment"))
      << "Line where '#' is not the first character should not be skipped.";
}

// -----------------------------------------------------------------------------
// Test the operator<< overload for std::array<double, 3>.
//
// The operator writes the three elements separated by single spaces.
// -----------------------------------------------------------------------------
void test_mc_array_ostream_operator() {
  std::cerr << "\n[TEST] class_Molecule_Complex.cpp :: "
               "operator<<(ostream&, std::array<double,3>)\n";
  std::cerr << "       Verifying that a 3-element array is streamed as "
               "space-separated values.\n";

  std::array<double, 3> arr{1.0, 2.0, 3.0};
  std::ostringstream oss;
  oss << arr;

  std::cerr << "       - Streamed representation: \"" << oss.str() << "\"\n";
  // The elements should appear separated by single spaces.
  EXPECT_EQ(oss.str(), std::string("1 2 3"))
      << "Array should be streamed as '1 2 3'.";
}

// -----------------------------------------------------------------------------
// Test Molecule::Iface::change_state().
//
// change_state() sets the relative state index, absolute index, and state
// identifier char of an interface.
// -----------------------------------------------------------------------------
void test_mc_iface_change_state() {
  std::cerr << "\n[TEST] class_Molecule_Complex.cpp :: "
               "Molecule::Iface::change_state()\n";
  std::cerr << "       Verifying that change_state updates stateIndex, index, "
               "and stateIden.\n";

  Molecule::Iface iface;
  std::cerr << "       - Calling change_state(2, 42, 'A').\n";
  iface.change_state(2, 42, 'A');

  EXPECT_EQ(iface.stateIndex, 2) << "stateIndex should be updated to 2.";
  EXPECT_EQ(iface.index, 42) << "index should be updated to 42.";
  EXPECT_EQ(iface.stateIden, 'A') << "stateIden should be updated to 'A'.";
}

// -----------------------------------------------------------------------------
// Test Molecule::Interaction::clear().
//
// clear() resets all partner and reaction bookkeeping fields to -1.
// -----------------------------------------------------------------------------
void test_mc_interaction_clear() {
  std::cerr << "\n[TEST] class_Molecule_Complex.cpp :: "
               "Molecule::Interaction::clear()\n";
  std::cerr << "       Verifying that clear() resets all interaction fields to "
               "-1.\n";

  Molecule::Interaction interaction;
  // Give the fields non-default values so we can confirm they are reset.
  interaction.partnerIfaceIndex = 7;
  interaction.partnerIndex = 8;
  interaction.partnerId = 9;
  interaction.conjBackRxn = 10;

  std::cerr << "       - Calling clear() after setting non-default values.\n";
  interaction.clear();

  EXPECT_EQ(interaction.partnerIfaceIndex, -1)
      << "partnerIfaceIndex should be reset to -1.";
  EXPECT_EQ(interaction.partnerIndex, -1)
      << "partnerIndex should be reset to -1.";
  EXPECT_EQ(interaction.partnerId, -1) << "partnerId should be reset to -1.";
  EXPECT_EQ(interaction.conjBackRxn, -1)
      << "conjBackRxn should be reset to -1.";
}

// -----------------------------------------------------------------------------
// Helper: build a minimal Molecule with a given index / type / com coord.
// This keeps the equality tests readable.
// -----------------------------------------------------------------------------
static Molecule make_simple_molecule(int myComIndex, int molTypeIndex,
                                      int index, double x, double y, double z) {
  Molecule mol;
  mol.myComIndex = myComIndex;
  mol.molTypeIndex = molTypeIndex;
  mol.index = index;
  mol.comCoord.x = x;
  mol.comCoord.y = y;
  mol.comCoord.z = z;
  return mol;
}

// -----------------------------------------------------------------------------
// Test Molecule::operator== and Molecule::operator!=.
//
// Two Molecules compare equal when their myComIndex, molTypeIndex, index, and
// comCoord all match.
// -----------------------------------------------------------------------------
void test_mc_molecule_equality() {
  std::cerr << "\n[TEST] class_Molecule_Complex.cpp :: "
               "Molecule::operator== / operator!=\n";
  std::cerr << "       Verifying equality is based on complex index, type "
               "index, index, and comCoord.\n";

  Molecule a = make_simple_molecule(0, 1, 5, 1.0, 2.0, 3.0);
  Molecule b = make_simple_molecule(0, 1, 5, 1.0, 2.0, 3.0);
  Molecule c = make_simple_molecule(0, 1, 6, 1.0, 2.0, 3.0);  // different index

  std::cerr << "       - Comparing two identical molecules (expect equal).\n";
  EXPECT_TRUE(a == b) << "Identical molecules should be equal.";
  EXPECT_FALSE(a != b) << "Identical molecules should not be unequal.";

  std::cerr << "       - Comparing molecules with different index (expect "
               "unequal).\n";
  EXPECT_TRUE(a != c) << "Molecules with different index should be unequal.";
  EXPECT_FALSE(a == c) << "Molecules with different index should not be equal.";
}

// -----------------------------------------------------------------------------
// Test Molecule::set_tmp_association_coords() and
// Molecule::clear_tmp_association_coords().
//
// set_tmp_association_coords copies comCoord into tmpComCoord and each
// interface coordinate into tmpICoords.  clear_tmp_association_coords zeros
// tmpComCoord and empties tmpICoords.
// -----------------------------------------------------------------------------
void test_mc_tmp_association_coords() {
  std::cerr << "\n[TEST] class_Molecule_Complex.cpp :: "
               "Molecule::set_tmp_association_coords() / "
               "clear_tmp_association_coords()\n";
  std::cerr << "       Verifying temporary association coordinates are set and "
               "cleared correctly.\n";

  // Build a molecule with a com coordinate and two interfaces.
  Molecule mol = make_simple_molecule(0, 0, 0, 4.0, 5.0, 6.0);
  Molecule::Iface iface1;
  iface1.coord.x = 1.0;
  iface1.coord.y = 1.0;
  iface1.coord.z = 1.0;
  Molecule::Iface iface2;
  iface2.coord.x = 2.0;
  iface2.coord.y = 2.0;
  iface2.coord.z = 2.0;
  mol.interfaceList.push_back(iface1);
  mol.interfaceList.push_back(iface2);

  std::cerr << "       - Calling set_tmp_association_coords().\n";
  mol.set_tmp_association_coords();

  // tmpComCoord should mirror comCoord.
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.x, 4.0) << "tmpComCoord.x should match comCoord.x.";
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.y, 5.0) << "tmpComCoord.y should match comCoord.y.";
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.z, 6.0) << "tmpComCoord.z should match comCoord.z.";

  // tmpICoords should contain both interface coordinates.
  ASSERT_EQ(mol.tmpICoords.size(), 2u)
      << "tmpICoords should hold both interface coordinates.";
  EXPECT_DOUBLE_EQ(mol.tmpICoords[0].x, 1.0);
  EXPECT_DOUBLE_EQ(mol.tmpICoords[1].x, 2.0);

  std::cerr << "       - Calling clear_tmp_association_coords().\n";
  mol.clear_tmp_association_coords();

  // After clearing, tmpICoords should be empty and tmpComCoord zeroed.
  EXPECT_TRUE(mol.tmpICoords.empty())
      << "tmpICoords should be empty after clearing.";
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.x, 0.0) << "tmpComCoord.x should be zeroed.";
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.y, 0.0) << "tmpComCoord.y should be zeroed.";
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.z, 0.0) << "tmpComCoord.z should be zeroed.";
}

// -----------------------------------------------------------------------------
// Test Molecule::update_association_coords().
//
// When tmpICoords is empty, update_association_coords fills tmpComCoord from
// (vec + comCoord) and pushes (vec + iface.coord) for each interface.  When
// tmpICoords already contains coords, it shifts the existing tmp coordinates
// by the vector instead.
// -----------------------------------------------------------------------------
void test_mc_update_association_coords() {
  std::cerr << "\n[TEST] class_Molecule_Complex.cpp :: "
               "Molecule::update_association_coords()\n";
  std::cerr << "       Verifying association coordinates are shifted by the "
               "supplied vector (empty and non-empty tmp cases).\n";

  Molecule mol = make_simple_molecule(0, 0, 0, 10.0, 10.0, 10.0);
  Molecule::Iface iface;
  iface.coord.x = 1.0;
  iface.coord.y = 1.0;
  iface.coord.z = 1.0;
  mol.interfaceList.push_back(iface);

  // Build a translation vector.
  Vector vec;
  vec.x = 2.0;
  vec.y = 3.0;
  vec.z = 4.0;

  std::cerr << "       - First call with empty tmpICoords (fills from "
               "comCoord/interfaceList).\n";
  mol.update_association_coords(vec);

  // tmpComCoord = vec + comCoord.
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.x, 12.0);
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.y, 13.0);
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.z, 14.0);
  ASSERT_EQ(mol.tmpICoords.size(), 1u);
  // tmpICoords[0] = vec + iface.coord.
  EXPECT_DOUBLE_EQ(mol.tmpICoords[0].x, 3.0);
  EXPECT_DOUBLE_EQ(mol.tmpICoords[0].y, 4.0);
  EXPECT_DOUBLE_EQ(mol.tmpICoords[0].z, 5.0);

  std::cerr << "       - Second call with non-empty tmpICoords (shifts "
               "existing tmp coords).\n";
  mol.update_association_coords(vec);

  // tmpComCoord shifted again: 12+2, 13+3, 14+4.
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.x, 14.0);
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.y, 16.0);
  EXPECT_DOUBLE_EQ(mol.tmpComCoord.z, 18.0);
  // tmpICoords shifted again: 3+2, 4+3, 5+4.
  EXPECT_DOUBLE_EQ(mol.tmpICoords[0].x, 5.0);
  EXPECT_DOUBLE_EQ(mol.tmpICoords[0].y, 7.0);
  EXPECT_DOUBLE_EQ(mol.tmpICoords[0].z, 9.0);
}

// -----------------------------------------------------------------------------
// Test Complex::translate().
//
// translate() shifts every member molecule's comCoord and every interface
// coordinate by the supplied vector.
// -----------------------------------------------------------------------------
void test_mc_complex_translate() {
  std::cerr << "\n[TEST] class_Molecule_Complex.cpp :: Complex::translate()\n";
  std::cerr << "       Verifying that every member molecule and its interfaces "
               "are shifted by the vector.\n";

  // Build two member molecules with interfaces.
  std::vector<Molecule> moleculeList;
  Molecule mol0 = make_simple_molecule(0, 0, 0, 0.0, 0.0, 0.0);
  Molecule::Iface iface0;
  iface0.coord.x = 1.0;
  iface0.coord.y = 1.0;
  iface0.coord.z = 1.0;
  mol0.interfaceList.push_back(iface0);

  Molecule mol1 = make_simple_molecule(0, 0, 1, 5.0, 5.0, 5.0);
  Molecule::Iface iface1;
  iface1.coord.x = 6.0;
  iface1.coord.y = 6.0;
  iface1.coord.z = 6.0;
  mol1.interfaceList.push_back(iface1);

  moleculeList.push_back(mol0);
  moleculeList.push_back(mol1);

  // Build a complex whose member list references molecules 0 and 1.
  Complex complexObj;
  complexObj.memberList.push_back(0);
  complexObj.memberList.push_back(1);

  // Translate by (1, 2, 3).
  Vector transVec;
  transVec.x = 1.0;
  transVec.y = 2.0;
  transVec.z = 3.0;

  std::cerr << "       - Translating complex by vector (1, 2, 3).\n";
  complexObj.translate(transVec, moleculeList);

  // Molecule 0's com should now be (1, 2, 3), its interface (2, 3, 4).
  EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 1.0);
  EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, 2.0);
  EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 3.0);
  EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, 2.0);
  EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.y, 3.0);
  EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.z, 4.0);

  // Molecule 1's com should now be (6, 7, 8), its interface (7, 8, 9).
  EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 6.0);
  EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.y, 7.0);
  EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.z, 8.0);
  EXPECT_DOUBLE_EQ(moleculeList[1].interfaceList[0].coord.x, 7.0);
  EXPECT_DOUBLE_EQ(moleculeList[1].interfaceList[0].coord.y, 8.0);
  EXPECT_DOUBLE_EQ(moleculeList[1].interfaceList[0].coord.z, 9.0);
}

// =============================================================================
// gtest wrappers.  Each TEST simply calls the corresponding test_* function so
// that the assertions above are executed and reported by the gtest framework.
// =============================================================================

TEST(MoleculeComplexTest, SkipLine) { test_mc_skipLine(); }
TEST(MoleculeComplexTest, ArrayOstreamOperator) {
  test_mc_array_ostream_operator();
}
TEST(MoleculeComplexTest, IfaceChangeState) { test_mc_iface_change_state(); }
TEST(MoleculeComplexTest, InteractionClear) { test_mc_interaction_clear(); }
TEST(MoleculeComplexTest, MoleculeEquality) { test_mc_molecule_equality(); }
TEST(MoleculeComplexTest, TmpAssociationCoords) {
  test_mc_tmp_association_coords();
}
TEST(MoleculeComplexTest, UpdateAssociationCoords) {
  test_mc_update_association_coords();
}
TEST(MoleculeComplexTest, ComplexTranslate) { test_mc_complex_translate(); }