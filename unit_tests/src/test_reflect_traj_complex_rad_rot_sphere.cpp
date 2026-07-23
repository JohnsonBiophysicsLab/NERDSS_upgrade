// Unit test for reflect_traj_complex_rad_rot_sphere.cpp
//
// This test exercises the free function
//    reflect_traj_complex_rad_rot_sphere(...)
// which is responsible for keeping a molecular "Complex" that lives *inside*
// a spherical boundary from crossing (spanning) that boundary. When a proposed
// translation ("trajTrans") would push the complex (or one of its interfaces)
// outside the sphere, the function reflects/moves the translation vector back
// so that the complex remains inside.
//
// Because the routine relies on several NERDSS data structures (Parameters,
// Molecule, Complex, Membrane, Coord, Vector, ...) we build up minimal but
// valid instances of these objects and then verify the observable behavior.

#include "boundary_conditions/reflect_functions.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

namespace {

// -----------------------------------------------------------------------------
// Small helper that prints a banner describing which source file / function is
// currently under test.  Keeps the console output verbose and self-describing.
// -----------------------------------------------------------------------------
void rtcrrs_print_banner(const std::string& what) {
    std::cerr << "\n========================================================\n";
    std::cerr << "[TEST] File   : reflect_traj_complex_rad_rot_sphere.cpp\n";
    std::cerr << "[TEST] Func   : reflect_traj_complex_rad_rot_sphere()\n";
    std::cerr << "[TEST] Case   : " << what << "\n";
    std::cerr << "========================================================\n";
}

// -----------------------------------------------------------------------------
// Helper that builds a single Molecule with a given center-of-mass coordinate
// and (by default) an empty interface list.  An empty interface list keeps the
// inner loop in the tested routine well-defined but simple.
// -----------------------------------------------------------------------------
Molecule rtcrrs_make_molecule(double x, double y, double z) {
    Molecule mol;
    mol.comCoord.x = x;
    mol.comCoord.y = y;
    mol.comCoord.z = z;
    // Leave mol.interfaceList empty so the interface loop is a no-op.
    return mol;
}

// -----------------------------------------------------------------------------
// Helper to build a Complex that references molecule index 0 in moleculeList.
// The trajRot is left at zero so the Euler rotation matrix is the identity,
// which makes the geometry easy to reason about analytically.
// -----------------------------------------------------------------------------
Complex rtcrrs_make_complex(double comX, double comY, double comZ,
                            double transX, double transY, double transZ,
                            double complexRadius, bool onSurface) {
    Complex com;
    com.index = 0;

    // Center of mass of the complex.
    com.comCoord.x = comX;
    com.comCoord.y = comY;
    com.comCoord.z = comZ;

    // Proposed translation for this timestep.
    com.trajTrans.x = transX;
    com.trajTrans.y = transY;
    com.trajTrans.z = transZ;

    // Zero rotation -> identity rotation matrix.
    com.trajRot.x = 0.0;
    com.trajRot.y = 0.0;
    com.trajRot.z = 0.0;

    com.radius = complexRadius;
    com.OnSurface = onSurface;

    // This complex is made up of a single member molecule (index 0).
    com.memberList.clear();
    com.memberList.push_back(0);

    return com;
}

// -----------------------------------------------------------------------------
// TEST 1: Complex flagged as OnSurface.
//
// According to the source, when targCom.OnSurface == true the routine takes the
// (empty) "on surface" branch and performs NO modification of trajTrans.
// We verify that the translation vector is returned completely unchanged.
// -----------------------------------------------------------------------------
void test_rtcrrs_on_surface_leaves_traj_unchanged() {
    rtcrrs_print_banner("OnSurface complex should NOT modify trajTrans");

    Parameters params;               // Default parameters are sufficient here.
    Membrane membraneObject;         // Default membrane object.
    const double radius = 100.0;     // Sphere radius passed to the routine.
    const double RS3Dinput = 5.0;    // Reflecting surface (ignored when OnSurface).

    // Build a molecule and a complex that is marked as being on the surface.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rtcrrs_make_molecule(0.0, 0.0, 0.0));

    // A translation that (in the "inside" case) would clearly escape the sphere.
    Complex targCom = rtcrrs_make_complex(0.0, 0.0, 0.0,
                                          500.0, 0.0, 0.0,
                                          /*complexRadius=*/10.0,
                                          /*onSurface=*/true);

    // Record the input translation so we can compare afterward.
    const double origX = targCom.trajTrans.x;
    const double origY = targCom.trajTrans.y;
    const double origZ = targCom.trajTrans.z;

    std::cerr << "[INFO] Calling routine with OnSurface = true, trajTrans = ("
              << origX << ", " << origY << ", " << origZ << ")\n";

    reflect_traj_complex_rad_rot_sphere(params, moleculeList, targCom,
                                        membraneObject, radius, RS3Dinput);

    // Criteria: because OnSurface complexes skip the reflection logic, the
    // translation vector must be exactly what we passed in.
    std::cerr << "[CHECK] trajTrans should be unchanged for OnSurface complex.\n";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, origX)
        << "trajTrans.x changed for an OnSurface complex.";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, origY)
        << "trajTrans.y changed for an OnSurface complex.";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, origZ)
        << "trajTrans.z changed for an OnSurface complex.";
}

// -----------------------------------------------------------------------------
// TEST 2: Complex safely inside the sphere.
//
// When the complex (comCoord + trajTrans magnitude + radius) stays within
// sphereR, the routine's outer distance test fails and trajTrans should remain
// unchanged.
// -----------------------------------------------------------------------------
void test_rtcrrs_inside_sphere_no_change() {
    rtcrrs_print_banner("Complex well inside sphere should NOT modify trajTrans");

    Parameters params;
    Membrane membraneObject;
    const double radius = 100.0;     // Sphere radius.
    const double RS3Dinput = 0.0;    // No reflecting-surface offset.

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rtcrrs_make_molecule(0.0, 0.0, 0.0));

    // Small translation, small radius: comCoord(0,0,0)+trajTrans(1,0,0) => r=1,
    // plus radius 2 => 3 << sphereR 100, so it clearly fits.
    Complex targCom = rtcrrs_make_complex(0.0, 0.0, 0.0,
                                          1.0, 0.0, 0.0,
                                          /*complexRadius=*/2.0,
                                          /*onSurface=*/false);

    const double origX = targCom.trajTrans.x;
    const double origY = targCom.trajTrans.y;
    const double origZ = targCom.trajTrans.z;

    std::cerr << "[INFO] Calling routine, complex fits comfortably inside sphere.\n";

    reflect_traj_complex_rad_rot_sphere(params, moleculeList, targCom,
                                        membraneObject, radius, RS3Dinput);

    // Criteria: the fast-path distance check should have found it fits, so
    // nothing about the translation should change.
    std::cerr << "[CHECK] trajTrans should be unchanged when complex fits inside.\n";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, origX)
        << "trajTrans.x changed even though complex fit inside the sphere.";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, origY)
        << "trajTrans.y changed even though complex fit inside the sphere.";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, origZ)
        << "trajTrans.z changed even though complex fit inside the sphere.";
}

// -----------------------------------------------------------------------------
// TEST 3: Complex whose proposed translation would escape the sphere.
//
// Here we choose a translation that pushes the single member molecule far
// outside the sphere. The routine should detect this ("outside == true"),
// reflect trajTrans back toward the interior, and set the recheck flag (which
// invokes reflect_traj_check_span_sphere).
//
// We analytically expect (with identity rotation, molecule at origin,
// comCoord at origin):
//    curr = trajTrans = (200,0,0), rtmp = 200, targR becomes 200
//    lamda = -2*(200 - sphereR)/200 with sphereR = radius = 100
//          = -2*(100)/200 = -1
//    trajTrans_new = trajTrans + lamda*targcrds = (200,0,0) + (-1)*(200,0,0)
//                  = (0,0,0)
// After the internal reflect_traj_check_span_sphere call the complex should be
// well within the sphere, so we mainly assert the translation was changed and
// that the resulting position no longer lies far outside the sphere.
// -----------------------------------------------------------------------------
void test_rtcrrs_reflect_when_outside() {
    rtcrrs_print_banner("Complex escaping sphere should be reflected inward");

    Parameters params;
    Membrane membraneObject;
    const double radius = 100.0;     // Sphere radius.
    const double RS3Dinput = 0.0;    // No reflecting-surface offset for clarity.

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rtcrrs_make_molecule(0.0, 0.0, 0.0));

    // Large translation pushes the complex far outside the 100-radius sphere.
    Complex targCom = rtcrrs_make_complex(0.0, 0.0, 0.0,
                                          200.0, 0.0, 0.0,
                                          /*complexRadius=*/10.0,
                                          /*onSurface=*/false);

    const double origX = targCom.trajTrans.x;
    const double origY = targCom.trajTrans.y;
    const double origZ = targCom.trajTrans.z;

    std::cerr << "[INFO] Calling routine with escaping trajTrans = ("
              << origX << ", " << origY << ", " << origZ << ")\n";

    reflect_traj_complex_rad_rot_sphere(params, moleculeList, targCom,
                                        membraneObject, radius, RS3Dinput);

    std::cerr << "[INFO] Post-call trajTrans = ("
              << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", "
              << targCom.trajTrans.z << ")\n";

    // Criteria 1: The translation must have been modified (it was escaping).
    std::cerr << "[CHECK] trajTrans must change because complex was escaping.\n";
    const bool changed =
        (std::fabs(targCom.trajTrans.x - origX) > 1e-9) ||
        (std::fabs(targCom.trajTrans.y - origY) > 1e-9) ||
        (std::fabs(targCom.trajTrans.z - origZ) > 1e-9);
    EXPECT_TRUE(changed)
        << "trajTrans was NOT modified even though the complex escaped the sphere.";

    // Criteria 2: The resulting COM position should be brought back inside (or
    // at least not far outside) the sphere.  We compute the resulting radial
    // distance of the complex COM after the reflection.
    const double newX = targCom.comCoord.x + targCom.trajTrans.x;
    const double newY = targCom.comCoord.y + targCom.trajTrans.y;
    const double newZ = targCom.comCoord.z + targCom.trajTrans.z;
    const double newR = std::sqrt(newX * newX + newY * newY + newZ * newZ);

    std::cerr << "[CHECK] Reflected COM radial distance (" << newR
              << ") should be <= sphere radius (" << radius << ").\n";
    EXPECT_LE(newR, radius + 1e-6)
        << "Complex COM is still outside the sphere after reflection (r=" << newR
        << ", sphereR=" << radius << ").";
}

}  // namespace

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each TEST simply invokes the corresponding test_*
// helper so that all assertions run even if some fail (no fatal ASSERT used).
// -----------------------------------------------------------------------------
TEST(ReflectTrajComplexRadRotSphereTest, OnSurfaceLeavesTrajUnchanged) {
    test_rtcrrs_on_surface_leaves_traj_unchanged();
}

TEST(ReflectTrajComplexRadRotSphereTest, InsideSphereNoChange) {
    test_rtcrrs_inside_sphere_no_change();
}

TEST(ReflectTrajComplexRadRotSphereTest, ReflectWhenOutside) {
    test_rtcrrs_reflect_when_outside();
}