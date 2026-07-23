#include "boundary_conditions/reflect_functions.hpp"

#include <array>
#include <iostream>
#include <vector>

#include "gtest/gtest.h"

// -----------------------------------------------------------------------------
// Unit tests for:
//   ../src/boundary_conditions/reflect_traj_tmp_crds_box.cpp
//   void reflect_traj_tmp_crds_box(const Parameters&, std::vector<Molecule>&,
//                                  Complex&, std::array<double,3>&,
//                                  const Membrane&, double)
//
// The function inspects the *tmp* coordinates of the molecules that belong to
// a target complex, together with a proposed translation vector `traj`, and
// modifies `traj` so that the complex is reflected back inside the simulation
// water box when a proposed move would have carried (part of) it outside.
//
// It does NOT modify the molecule / complex coordinates directly; it only
// mutates the passed-in `traj` array. Therefore all assertions below inspect
// the resulting values of `traj`.
// -----------------------------------------------------------------------------

namespace {

// Convenient constant for the (cubic) test box we use throughout.
constexpr double kBoxSide = 100.0;   // Box spans [-50, +50] in each axis.
constexpr double kHalfBox = kBoxSide / 2.0;

// -----------------------------------------------------------------------------
// Helper: build a single-molecule complex centered at `com` with one interface
// located at `com + ifaceOffset`.
//
// The helper wires up all the fields that reflect_traj_tmp_crds_box() actually
// reads:
//   * Molecule::tmpComCoord      - molecule temporary center of mass
//   * Molecule::tmpICoords       - temporary interface coordinates
//   * Molecule::interfaceList    - only its .size() is used
//   * Complex::tmpComCoord       - complex temporary center of mass
//   * Complex::memberList        - indices into moleculeList
//   * Complex::radius            - coarse bounding radius
//   * Complex::OnSurface / tmpOnSurface
// -----------------------------------------------------------------------------
void BuildSingleMoleculeComplex(std::vector<Molecule>& moleculeList, Complex& com,
                                const Coord& comCoord, const Coord& ifaceOffset,
                                double radius) {
    // Create one molecule and give it a temporary COM and one interface.
    Molecule mol {};
    mol.tmpComCoord = comCoord;

    // Interface is expressed as an absolute temporary coordinate.
    Coord iface { comCoord.x + ifaceOffset.x, comCoord.y + ifaceOffset.y,
                  comCoord.z + ifaceOffset.z };
    mol.tmpICoords.clear();
    mol.tmpICoords.push_back(iface);

    // interfaceList only needs the same size as tmpICoords (only .size() used).
    mol.interfaceList.resize(mol.tmpICoords.size());

    moleculeList.clear();
    moleculeList.push_back(mol);

    // Configure the complex to contain molecule index 0.
    com.memberList.clear();
    com.memberList.push_back(0);
    com.tmpComCoord = comCoord;
    com.radius = radius;
    com.index = 0;
    com.OnSurface = false;
    com.tmpOnSurface = false;
}

// Helper to create a cubic membrane/water-box of side kBoxSide.
Membrane MakeCubicMembrane() {
    Membrane mem {};
    mem.waterBox.x = kBoxSide;
    mem.waterBox.y = kBoxSide;
    mem.waterBox.z = kBoxSide;
    return mem;
}

// -----------------------------------------------------------------------------
// Test 1: A complex safely inside the box must leave `traj` unchanged.
// -----------------------------------------------------------------------------
void test_reflectTmpBox_InsideBoxUnchanged() {
    std::cerr << "\n[TEST] reflect_traj_tmp_crds_box: complex fully inside box\n"
              << "  Source: reflect_traj_tmp_crds_box.cpp\n"
              << "  Criteria: traj must be returned unmodified.\n";

    Parameters params {};
    Membrane mem = MakeCubicMembrane();
    std::vector<Molecule> moleculeList;
    Complex com {};

    // Small complex at the origin, well within [-50, +50].
    BuildSingleMoleculeComplex(moleculeList, com, Coord { 0.0, 0.0, 0.0 },
                               Coord { 5.0, 0.0, 0.0 }, /*radius=*/10.0);

    // Proposed small move that keeps everything comfortably inside the box.
    std::array<double, 3> traj { 5.0, -3.0, 2.0 };
    const std::array<double, 3> expected = traj;

    reflect_traj_tmp_crds_box(params, moleculeList, com, traj, mem, /*RS3Dinput=*/0.0);

    std::cerr << "  traj after = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
    EXPECT_DOUBLE_EQ(traj[0], expected[0]) << "traj[0] should be unchanged inside box";
    EXPECT_DOUBLE_EQ(traj[1], expected[1]) << "traj[1] should be unchanged inside box";
    EXPECT_DOUBLE_EQ(traj[2], expected[2]) << "traj[2] should be unchanged inside box";
}

// -----------------------------------------------------------------------------
// Test 2: A move that carries the complex past +X must be reflected back.
//
// Setup: box half-width = 50, COM at origin, interface at +5 in X, traj[0]=60.
//   farthest +X point (posWall) = 0 + 60 + 5 = 65
//   posWall exceeds +50 by 15  -> traj[0] -= 2*(65-50) = 30
//   expected traj[0] = 60 - 30 = 30
// -----------------------------------------------------------------------------
void test_reflectTmpBox_OutsidePosX() {
    std::cerr << "\n[TEST] reflect_traj_tmp_crds_box: reflect off +X wall\n"
              << "  Source: reflect_traj_tmp_crds_box.cpp\n"
              << "  Criteria: traj[0] reflected by 2*(posWall-posX); Y,Z unchanged.\n";

    Parameters params {};
    Membrane mem = MakeCubicMembrane();
    std::vector<Molecule> moleculeList;
    Complex com {};

    BuildSingleMoleculeComplex(moleculeList, com, Coord { 0.0, 0.0, 0.0 },
                               Coord { 5.0, 0.0, 0.0 }, /*radius=*/10.0);

    std::array<double, 3> traj { 60.0, 0.0, 0.0 };

    reflect_traj_tmp_crds_box(params, moleculeList, com, traj, mem, /*RS3Dinput=*/0.0);

    std::cerr << "  traj after = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
    // posWall = 65, posX = 50 -> traj[0] = 60 - 2*(65-50) = 30
    EXPECT_DOUBLE_EQ(traj[0], 30.0) << "traj[0] should be reflected back inside +X wall";
    EXPECT_DOUBLE_EQ(traj[1], 0.0) << "traj[1] should be untouched";
    EXPECT_DOUBLE_EQ(traj[2], 0.0) << "traj[2] should be untouched";
}

// -----------------------------------------------------------------------------
// Test 3: A move past -X must be reflected back the other way.
//
// Setup: COM at origin, interface at +5 in X, traj[0] = -60.
//   negWall (farthest -X) = COM point 0 + (-60) + 0 = -60
//   negWall below -50 by -10 -> traj[0] -= 2*(negWall-negX) = 2*(-60-(-50)) = -20
//   expected traj[0] = -60 - (-20) = -40
// -----------------------------------------------------------------------------
void test_reflectTmpBox_OutsideNegX() {
    std::cerr << "\n[TEST] reflect_traj_tmp_crds_box: reflect off -X wall\n"
              << "  Source: reflect_traj_tmp_crds_box.cpp\n"
              << "  Criteria: traj[0] reflected by 2*(negWall-negX).\n";

    Parameters params {};
    Membrane mem = MakeCubicMembrane();
    std::vector<Molecule> moleculeList;
    Complex com {};

    BuildSingleMoleculeComplex(moleculeList, com, Coord { 0.0, 0.0, 0.0 },
                               Coord { 5.0, 0.0, 0.0 }, /*radius=*/10.0);

    std::array<double, 3> traj { -60.0, 0.0, 0.0 };

    reflect_traj_tmp_crds_box(params, moleculeList, com, traj, mem, /*RS3Dinput=*/0.0);

    std::cerr << "  traj after = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
    EXPECT_DOUBLE_EQ(traj[0], -40.0) << "traj[0] should be reflected back inside -X wall";
    EXPECT_DOUBLE_EQ(traj[1], 0.0) << "traj[1] should be untouched";
    EXPECT_DOUBLE_EQ(traj[2], 0.0) << "traj[2] should be untouched";
}

// -----------------------------------------------------------------------------
// Test 4: A move past +Y must be reflected on the Y axis only.
//
// Setup: COM at origin, interface at +5 in Y, traj[1] = 60.
//   posWall = 0 + 60 + 5 = 65 -> traj[1] = 60 - 2*(65-50) = 30
// -----------------------------------------------------------------------------
void test_reflectTmpBox_OutsidePosY() {
    std::cerr << "\n[TEST] reflect_traj_tmp_crds_box: reflect off +Y wall\n"
              << "  Source: reflect_traj_tmp_crds_box.cpp\n"
              << "  Criteria: traj[1] reflected; X,Z unchanged.\n";

    Parameters params {};
    Membrane mem = MakeCubicMembrane();
    std::vector<Molecule> moleculeList;
    Complex com {};

    BuildSingleMoleculeComplex(moleculeList, com, Coord { 0.0, 0.0, 0.0 },
                               Coord { 0.0, 5.0, 0.0 }, /*radius=*/10.0);

    std::array<double, 3> traj { 0.0, 60.0, 0.0 };

    reflect_traj_tmp_crds_box(params, moleculeList, com, traj, mem, /*RS3Dinput=*/0.0);

    std::cerr << "  traj after = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
    EXPECT_DOUBLE_EQ(traj[1], 30.0) << "traj[1] should be reflected back inside +Y wall";
    EXPECT_DOUBLE_EQ(traj[0], 0.0) << "traj[0] should be untouched";
    EXPECT_DOUBLE_EQ(traj[2], 0.0) << "traj[2] should be untouched";
}

// -----------------------------------------------------------------------------
// Test 5: A move past +Z must be reflected on the Z axis only.
//
// Setup: COM at origin, interface at +5 in Z, traj[2] = 60.
//   posWall = 0 + 60 + 5 = 65 -> traj[2] = 60 - 2*(65-50) = 30
// (RS3Dinput = 0 and OnSurface=false so negZ = -50, posZ = +50.)
// -----------------------------------------------------------------------------
void test_reflectTmpBox_OutsidePosZ() {
    std::cerr << "\n[TEST] reflect_traj_tmp_crds_box: reflect off +Z wall\n"
              << "  Source: reflect_traj_tmp_crds_box.cpp\n"
              << "  Criteria: traj[2] reflected; X,Y unchanged.\n";

    Parameters params {};
    Membrane mem = MakeCubicMembrane();
    std::vector<Molecule> moleculeList;
    Complex com {};

    BuildSingleMoleculeComplex(moleculeList, com, Coord { 0.0, 0.0, 0.0 },
                               Coord { 0.0, 0.0, 5.0 }, /*radius=*/10.0);

    std::array<double, 3> traj { 0.0, 0.0, 60.0 };

    reflect_traj_tmp_crds_box(params, moleculeList, com, traj, mem, /*RS3Dinput=*/0.0);

    std::cerr << "  traj after = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
    EXPECT_DOUBLE_EQ(traj[2], 30.0) << "traj[2] should be reflected back inside +Z wall";
    EXPECT_DOUBLE_EQ(traj[0], 0.0) << "traj[0] should be untouched";
    EXPECT_DOUBLE_EQ(traj[1], 0.0) << "traj[1] should be untouched";
}

// -----------------------------------------------------------------------------
// Test 6: Verify that RS3Dinput raises the -Z wall when the complex is not on
//         the surface. With RS3Dinput=5, negZ = -50 + 5 = -45, so a move that
//         would sit at -46 (inside the geometric box but below negZ) reflects.
//
// Setup: COM at origin, interface at -5 in Z, traj[2] = -41.
//   negWall = 0 + (-41) + (-5) = -46
//   negWall < negZ (-45) by -1 -> traj[2] -= 2*(negWall - negZ) = 2*(-46-(-45)) = -2
//   expected traj[2] = -41 - (-2) = -39
// -----------------------------------------------------------------------------
void test_reflectTmpBox_RS3DRaisesNegZWall() {
    std::cerr << "\n[TEST] reflect_traj_tmp_crds_box: RS3Dinput raises -Z wall\n"
              << "  Source: reflect_traj_tmp_crds_box.cpp\n"
              << "  Criteria: negZ = -halfBox + RS3D; reflection uses raised wall.\n";

    Parameters params {};
    Membrane mem = MakeCubicMembrane();
    std::vector<Molecule> moleculeList;
    Complex com {};

    BuildSingleMoleculeComplex(moleculeList, com, Coord { 0.0, 0.0, 0.0 },
                               Coord { 0.0, 0.0, -5.0 }, /*radius=*/10.0);
    // Ensure it is treated as a volume complex (not on surface) so RS3D applies.
    com.OnSurface = false;
    com.tmpOnSurface = false;

    std::array<double, 3> traj { 0.0, 0.0, -41.0 };

    reflect_traj_tmp_crds_box(params, moleculeList, com, traj, mem, /*RS3Dinput=*/5.0);

    std::cerr << "  traj after = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
    EXPECT_DOUBLE_EQ(traj[2], -39.0) << "traj[2] should reflect off the RS3D-raised -Z wall";
    EXPECT_DOUBLE_EQ(traj[0], 0.0) << "traj[0] should be untouched";
    EXPECT_DOUBLE_EQ(traj[1], 0.0) << "traj[1] should be untouched";
}

// -----------------------------------------------------------------------------
// Test 7: When the complex is on the surface, RS3D is forced to 0, so negZ
//         collapses to -halfBox even if RS3Dinput is non-zero. The same move
//         from Test 6 should now stay inside and NOT be reflected.
// -----------------------------------------------------------------------------
void test_reflectTmpBox_OnSurfaceIgnoresRS3D() {
    std::cerr << "\n[TEST] reflect_traj_tmp_crds_box: OnSurface forces RS3D=0\n"
              << "  Source: reflect_traj_tmp_crds_box.cpp\n"
              << "  Criteria: negZ = -halfBox; move at -46 stays inside, no reflection.\n";

    Parameters params {};
    Membrane mem = MakeCubicMembrane();
    std::vector<Molecule> moleculeList;
    Complex com {};

    BuildSingleMoleculeComplex(moleculeList, com, Coord { 0.0, 0.0, 0.0 },
                               Coord { 0.0, 0.0, -5.0 }, /*radius=*/10.0);
    // Marking the complex on-surface should zero out RS3D internally.
    com.OnSurface = true;
    com.tmpOnSurface = true;

    std::array<double, 3> traj { 0.0, 0.0, -41.0 };
    const std::array<double, 3> expected = traj;  // -46 is still >= -50, no reflection

    reflect_traj_tmp_crds_box(params, moleculeList, com, traj, mem, /*RS3Dinput=*/5.0);

    std::cerr << "  traj after = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
    EXPECT_DOUBLE_EQ(traj[2], expected[2]) << "traj[2] should NOT reflect when on surface";
    EXPECT_DOUBLE_EQ(traj[0], expected[0]) << "traj[0] should be untouched";
    EXPECT_DOUBLE_EQ(traj[1], expected[1]) << "traj[1] should be untouched";
}

}  // namespace

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each TEST simply invokes one of the named helper
// functions above so that failures are individually reported while every test
// still runs (all assertions are non-fatal EXPECT_* macros).
// -----------------------------------------------------------------------------
TEST(ReflectTrajTmpCrdsBoxTest, InsideBoxUnchanged) {
    test_reflectTmpBox_InsideBoxUnchanged();
}

TEST(ReflectTrajTmpCrdsBoxTest, OutsidePosX) {
    test_reflectTmpBox_OutsidePosX();
}

TEST(ReflectTrajTmpCrdsBoxTest, OutsideNegX) {
    test_reflectTmpBox_OutsideNegX();
}

TEST(ReflectTrajTmpCrdsBoxTest, OutsidePosY) {
    test_reflectTmpBox_OutsidePosY();
}

TEST(ReflectTrajTmpCrdsBoxTest, OutsidePosZ) {
    test_reflectTmpBox_OutsidePosZ();
}

TEST(ReflectTrajTmpCrdsBoxTest, RS3DRaisesNegZWall) {
    test_reflectTmpBox_RS3DRaisesNegZWall();
}

TEST(ReflectTrajTmpCrdsBoxTest, OnSurfaceIgnoresRS3D) {
    test_reflectTmpBox_OnSurfaceIgnoresRS3D();
}