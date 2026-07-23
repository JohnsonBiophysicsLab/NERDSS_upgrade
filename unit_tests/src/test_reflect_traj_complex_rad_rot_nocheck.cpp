#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

#include <gtest/gtest.h>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Unit tests for:
//   ../src/boundary_conditions/reflect_traj_complex_rad_rot_nocheck.cpp
//
// The function under test:
//   void reflect_traj_complex_rad_rot_nocheck(const Parameters& params,
//                                             Complex& targCom,
//                                             std::vector<Molecule>& moleculeList,
//                                             const Membrane& membraneObject,
//                                             double RS3Dinput)
//
// This is a *dispatch* function.  It inspects membraneObject.isSphere and then
// forwards the call to one of two implementations:
//    - reflect_traj_complex_rad_rot_nocheck_sphere(...)  when isSphere == true
//    - reflect_traj_complex_rad_rot_nocheck_box(...)      when isSphere == false
//
// Because the actual reflection math lives in the sub-functions, these tests
// focus on:
//    1. Verifying that the dispatch routine runs without crashing on both
//       branches (sphere / box).
//    2. Verifying that, for a box geometry, a complex whose trial translation
//       would push it beyond a wall gets its trajectory adjusted (reflected)
//       back inward.
// ---------------------------------------------------------------------------

namespace {

// Helper that builds a very simple, single-molecule Complex sitting at the
// origin so that we can exercise the reflection routine in a controlled way.
//
// The returned moleculeList is filled by reference; the returned Complex holds
// its member index (0) and its own bookkeeping coordinates.
Complex reflect_dispatch_MakeSimpleComplex(std::vector<Molecule>& moleculeList)
{
    moleculeList.clear();

    // Build a single molecule located at the origin with one interface.
    Molecule mol {};
    mol.comCoord = Coord { 0.0, 0.0, 0.0 };

    // Give the molecule a single interface offset slightly from the COM so the
    // rotation matrix has something non-trivial to act upon.
    Molecule::Iface iface {};
    iface.coord = Coord { 1.0, 0.0, 0.0 };
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    moleculeList.push_back(mol);

    // Build the complex that owns molecule index 0.
    Complex targCom {};
    targCom.comCoord = Coord { 0.0, 0.0, 0.0 };
    targCom.radius = 1.0;
    targCom.OnSurface = false;
    targCom.memberList.clear();
    targCom.memberList.push_back(0);

    // Start with no trial translation; individual tests set this as needed.
    targCom.trajTrans = Vector { 0.0, 0.0, 0.0 };

    return targCom;
}

}  // namespace

// ---------------------------------------------------------------------------
// Test 1: Box branch executes and reflects an out-of-bounds complex.
// ---------------------------------------------------------------------------
void test_reflect_dispatch_box_branch_reflects()
{
    std::cerr << "\n[TEST] reflect_traj_complex_rad_rot_nocheck (BOX branch)\n"
              << "  File under test: reflect_traj_complex_rad_rot_nocheck.cpp\n"
              << "  Criteria: with isSphere=false the box implementation is\n"
              << "            invoked and a complex pushed past the +X wall has\n"
              << "            its trajTrans.x reduced (reflected inward).\n";

    // Build the molecule list and complex.
    std::vector<Molecule> moleculeList;
    Complex targCom = reflect_dispatch_MakeSimpleComplex(moleculeList);

    // Configure a cubic water box of side length 20 (walls at +/-10).
    Membrane membraneObject {};
    membraneObject.isSphere = false;
    membraneObject.waterBox.x = 20.0;
    membraneObject.waterBox.y = 20.0;
    membraneObject.waterBox.z = 20.0;

    // Give the complex a large +X translation that would drive it well past
    // the +X wall (10.0).  We expect the routine to reflect it back inward.
    targCom.trajTrans = Vector { 50.0, 0.0, 0.0 };
    const double originalTransX = targCom.trajTrans.x;

    Parameters params {};
    const double RS3Dinput = 0.0;

    std::cerr << "  Initial trajTrans.x = " << originalTransX
              << " (would place complex outside +X wall at 10.0)\n";

    // Call the dispatch function; must not throw / crash.
    reflect_traj_complex_rad_rot_nocheck(params, targCom, moleculeList,
                                         membraneObject, RS3Dinput);

    std::cerr << "  Reflected trajTrans.x = " << targCom.trajTrans.x << "\n";

    // The reflected translation should be strictly smaller than the original
    // (it was pushed back toward the interior of the box).
    EXPECT_LT(targCom.trajTrans.x, originalTransX)
        << "Box branch should reduce trajTrans.x for an out-of-bounds complex.";
}

// ---------------------------------------------------------------------------
// Test 2: Box branch leaves an in-bounds complex essentially unchanged.
// ---------------------------------------------------------------------------
void test_reflect_dispatch_box_branch_inside_unchanged()
{
    std::cerr << "\n[TEST] reflect_traj_complex_rad_rot_nocheck (BOX branch, inside)\n"
              << "  Criteria: a complex whose translation keeps it inside all\n"
              << "            walls should have its trajTrans left unchanged.\n";

    std::vector<Molecule> moleculeList;
    Complex targCom = reflect_dispatch_MakeSimpleComplex(moleculeList);

    Membrane membraneObject {};
    membraneObject.isSphere = false;
    membraneObject.waterBox.x = 20.0;
    membraneObject.waterBox.y = 20.0;
    membraneObject.waterBox.z = 20.0;

    // A tiny translation keeps everything well inside the box.
    targCom.trajTrans = Vector { 0.5, 0.5, 0.5 };
    const Coord original = targCom.trajTrans;

    Parameters params {};
    const double RS3Dinput = 0.0;

    reflect_traj_complex_rad_rot_nocheck(params, targCom, moleculeList,
                                         membraneObject, RS3Dinput);

    std::cerr << "  trajTrans after call: (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n";

    // In-bounds complex should be untouched (allow a tiny numerical tolerance).
    EXPECT_NEAR(targCom.trajTrans.x, original.x, 1e-9)
        << "trajTrans.x should be unchanged for in-bounds complex.";
    EXPECT_NEAR(targCom.trajTrans.y, original.y, 1e-9)
        << "trajTrans.y should be unchanged for in-bounds complex.";
    EXPECT_NEAR(targCom.trajTrans.z, original.z, 1e-9)
        << "trajTrans.z should be unchanged for in-bounds complex.";
}

// ---------------------------------------------------------------------------
// Test 3: Sphere branch executes without crashing.
// ---------------------------------------------------------------------------
void test_reflect_dispatch_sphere_branch_runs()
{
    std::cerr << "\n[TEST] reflect_traj_complex_rad_rot_nocheck (SPHERE branch)\n"
              << "  Criteria: with isSphere=true the sphere implementation is\n"
              << "            invoked and the call completes without error.\n";

    std::vector<Molecule> moleculeList;
    Complex targCom = reflect_dispatch_MakeSimpleComplex(moleculeList);

    // Configure a spherical membrane geometry.
    Membrane membraneObject {};
    membraneObject.isSphere = true;
    membraneObject.sphereR = 10.0;   // sphere radius

    // Small translation so the sphere routine has valid geometry to work with.
    targCom.trajTrans = Vector { 1.0, 0.0, 0.0 };

    Parameters params {};
    const double RS3Dinput = 0.0;

    std::cerr << "  Invoking sphere branch with sphereR = "
              << membraneObject.sphereR << "\n";

    // The main assertion here is simply that the dispatch reaches the sphere
    // implementation and returns without throwing.
    EXPECT_NO_THROW({
        reflect_traj_complex_rad_rot_nocheck(params, targCom, moleculeList,
                                             membraneObject, RS3Dinput);
    }) << "Sphere branch should execute without throwing.";

    std::cerr << "  Sphere branch completed. trajTrans = ("
              << targCom.trajTrans.x << ", " << targCom.trajTrans.y << ", "
              << targCom.trajTrans.z << ")\n";
}

// ---------------------------------------------------------------------------
// GoogleTest wrappers: each TEST forwards to the corresponding test_* helper.
// Using EXPECT_* (non-fatal) assertions ensures all tests run to completion.
// ---------------------------------------------------------------------------
TEST(ReflectTrajComplexRadRotNoCheckTest, BoxBranchReflects)
{
    test_reflect_dispatch_box_branch_reflects();
}

TEST(ReflectTrajComplexRadRotNoCheckTest, BoxBranchInsideUnchanged)
{
    test_reflect_dispatch_box_branch_inside_unchanged();
}

TEST(ReflectTrajComplexRadRotNoCheckTest, SphereBranchRuns)
{
    test_reflect_dispatch_sphere_branch_runs();
}