#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

#include <array>
#include <vector>

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Unit tests for: ../src/boundary_conditions/reflect_traj_tmp_crds.cpp
//
// The function under test is:
//     void reflect_traj_tmp_crds(const Parameters& params,
//                                std::vector<Molecule>& moleculeList,
//                                Complex& targCom,
//                                std::array<double, 3>& traj,
//                                const Membrane& membraneObject,
//                                double RS3Dinput,
//                                bool isInsideCompartment)
//
// This routine is a *dispatcher*.  Based on the geometry flags in the passed
// Membrane object and the compartment flag it decides which of the lower-level
// reflection routines to call:
//     - reflect_traj_tmp_crds_sphere()       (box is a sphere OR inside compartment)
//     - reflect_traj_tmp_crds_box()          (box is a rectangular water box)
//     - reflect_traj_tmp_crds_compartment()  (protein enforces compartment BC)
//
// Because the routine only ever *reads* the input state and *modifies* the
// passed in `traj` array (it does not modify molecule trajectory vectors),
// these tests verify that:
//   1. Calling the function does not crash for the various geometry paths.
//   2. When the complex is well inside the boundary the trajectory is left
//      unchanged (no reflection needed).
//   3. When the complex would leave the box, the trajectory is modified so that
//      the complex is reflected back inside.
// -----------------------------------------------------------------------------

namespace {

// Helper: build a minimal, self-consistent single-molecule Complex placed at a
// requested center-of-mass position, so we can drive the dispatcher along its
// various code paths.  All temporary coordinates are made consistent with the
// COM as required by the routine's contract.
void ReflectTrajTmpCrds_BuildSingleMoleculeComplex(
    std::vector<Molecule>& moleculeList, Complex& targCom, const Coord& comPos)
{
    moleculeList.clear();

    // Create one molecule whose tmp COM matches the complex tmp COM.
    Molecule mol {};
    mol.index = 0;
    mol.tmpComCoord = comPos;
    // Give the molecule a single interface at the COM (keeps geometry simple).
    mol.interfaceList.resize(1);
    mol.tmpICoords.clear();
    mol.tmpICoords.push_back(comPos);
    // By default do not enforce compartment boundary conditions.
    mol.enforceCompartmentBC = false;
    moleculeList.push_back(mol);

    // Configure the complex to contain the single molecule above.
    targCom.index = 0;
    targCom.memberList.clear();
    targCom.memberList.push_back(0);
    targCom.tmpComCoord = comPos;
    targCom.radius = 1.0;       // small radius so it fits well inside the box
    targCom.OnSurface = false;  // not restricted to a surface
    targCom.D.z = 1.0;          // free to move in z
}

// Helper: build a rectangular-box Membrane (not a sphere) with a symmetric
// water box of the requested edge length.
Membrane ReflectTrajTmpCrds_BuildBoxMembrane(double edge)
{
    Membrane membrane {};
    membrane.isSphere = false;
    membrane.waterBox.x = edge;
    membrane.waterBox.y = edge;
    membrane.waterBox.z = edge;
    membrane.sphereR = 0.0;
    membrane.compartmentR = 0.0;
    return membrane;
}

// Helper: build a spherical Membrane with the requested sphere radius.
Membrane ReflectTrajTmpCrds_BuildSphereMembrane(double sphereR, double compartmentR)
{
    Membrane membrane {};
    membrane.isSphere = true;
    membrane.sphereR = sphereR;
    membrane.compartmentR = compartmentR;
    membrane.waterBox.x = 2.0 * sphereR;
    membrane.waterBox.y = 2.0 * sphereR;
    membrane.waterBox.z = 2.0 * sphereR;
    return membrane;
}

// -----------------------------------------------------------------------------
// Test 1: Box geometry, complex well inside the box.
// Criteria: The trajectory must remain unchanged since no reflection is needed.
// -----------------------------------------------------------------------------
void test_reflect_traj_tmp_crds_box_inside_no_change()
{
    std::cerr << "[RUN] test_reflect_traj_tmp_crds_box_inside_no_change\n";
    std::cerr << "  Testing reflect_traj_tmp_crds() (box path, inside box) "
                 "in reflect_traj_tmp_crds.cpp\n";

    Parameters params {};
    std::vector<Molecule> moleculeList;
    Complex targCom {};

    // Place the complex at the box origin (deep inside a 100 unit box).
    ReflectTrajTmpCrds_BuildSingleMoleculeComplex(moleculeList, targCom, Coord { 0.0, 0.0, 0.0 });
    Membrane membrane = ReflectTrajTmpCrds_BuildBoxMembrane(100.0);

    // A small trajectory step that keeps the complex safely inside the box.
    std::array<double, 3> traj { 1.0, 2.0, 3.0 };
    const std::array<double, 3> expected = traj;  // expect no modification

    std::cerr << "  Input traj = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";

    // isInsideCompartment == false so the box branch is taken.
    reflect_traj_tmp_crds(params, moleculeList, targCom, traj, membrane, 0.0, false);

    std::cerr << "  Output traj = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";

    // A complex deep inside the box should not be reflected.
    EXPECT_NEAR(traj[0], expected[0], 1e-9) << "traj[0] changed unexpectedly (inside box)";
    EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] changed unexpectedly (inside box)";
    EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] changed unexpectedly (inside box)";
}

// -----------------------------------------------------------------------------
// Test 2: Box geometry, complex pushed past the +X wall.
// Criteria: The x-component of the trajectory must be reduced (reflected back)
//           so that it does not carry the complex out of the box.
// -----------------------------------------------------------------------------
void test_reflect_traj_tmp_crds_box_reflect_posX()
{
    std::cerr << "[RUN] test_reflect_traj_tmp_crds_box_reflect_posX\n";
    std::cerr << "  Testing reflect_traj_tmp_crds() (box path, reflection off +X wall) "
                 "in reflect_traj_tmp_crds.cpp\n";

    Parameters params {};
    std::vector<Molecule> moleculeList;
    Complex targCom {};

    // Box edge = 20 -> walls at +/-10.  Place the complex near the +X wall.
    ReflectTrajTmpCrds_BuildSingleMoleculeComplex(moleculeList, targCom, Coord { 9.0, 0.0, 0.0 });
    Membrane membrane = ReflectTrajTmpCrds_BuildBoxMembrane(20.0);

    // A large +X step that would carry the complex out through the +X wall.
    std::array<double, 3> traj { 5.0, 0.0, 0.0 };

    std::cerr << "  Complex COM x = " << targCom.tmpComCoord.x
              << ", +X wall = " << membrane.waterBox.x / 2.0 << "\n";
    std::cerr << "  Input traj[0] = " << traj[0] << " (would end at x="
              << targCom.tmpComCoord.x + traj[0] << ")\n";

    reflect_traj_tmp_crds(params, moleculeList, targCom, traj, membrane, 0.0, false);

    std::cerr << "  Output traj[0] = " << traj[0] << " (now ends at x="
              << targCom.tmpComCoord.x + traj[0] << ")\n";

    // After reflection the resulting x position must be back inside the box.
    const double finalX = targCom.tmpComCoord.x + traj[0];
    EXPECT_LE(finalX, membrane.waterBox.x / 2.0 + 1e-9)
        << "Complex still outside +X wall after reflection";
    // The trajectory must have actually been reduced (reflected).
    EXPECT_LT(traj[0], 5.0) << "traj[0] was not reduced despite leaving the box";
}

// -----------------------------------------------------------------------------
// Test 3: Sphere geometry, complex near the center.
// Criteria: The dispatcher must take the sphere path without crashing and, for
//           a complex well inside the sphere, leave the trajectory unchanged.
// -----------------------------------------------------------------------------
void test_reflect_traj_tmp_crds_sphere_inside_no_change()
{
    std::cerr << "[RUN] test_reflect_traj_tmp_crds_sphere_inside_no_change\n";
    std::cerr << "  Testing reflect_traj_tmp_crds() (sphere path, inside sphere) "
                 "in reflect_traj_tmp_crds.cpp\n";

    Parameters params {};
    std::vector<Molecule> moleculeList;
    Complex targCom {};

    // Sphere radius = 100, complex at center.
    ReflectTrajTmpCrds_BuildSingleMoleculeComplex(moleculeList, targCom, Coord { 0.0, 0.0, 0.0 });
    Membrane membrane = ReflectTrajTmpCrds_BuildSphereMembrane(100.0, 50.0);

    // Small step that stays well inside the sphere.
    std::array<double, 3> traj { 1.0, 1.0, 1.0 };
    const std::array<double, 3> expected = traj;

    std::cerr << "  Sphere radius = " << membrane.sphereR << "\n";
    std::cerr << "  Input traj = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";

    reflect_traj_tmp_crds(params, moleculeList, targCom, traj, membrane, 0.0, false);

    std::cerr << "  Output traj = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";

    // Deep inside a big sphere: no reflection expected.
    EXPECT_NEAR(traj[0], expected[0], 1e-9) << "traj[0] changed unexpectedly (inside sphere)";
    EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] changed unexpectedly (inside sphere)";
    EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] changed unexpectedly (inside sphere)";
}

// -----------------------------------------------------------------------------
// Test 4: isInsideCompartment == true path.
// Criteria: The dispatcher must call reflect_traj_tmp_crds_sphere using the
//           compartment radius.  For a complex near the center of a large
//           compartment the trajectory should not change.
// -----------------------------------------------------------------------------
void test_reflect_traj_tmp_crds_inside_compartment()
{
    std::cerr << "[RUN] test_reflect_traj_tmp_crds_inside_compartment\n";
    std::cerr << "  Testing reflect_traj_tmp_crds() (inside-compartment path) "
                 "in reflect_traj_tmp_crds.cpp\n";

    Parameters params {};
    std::vector<Molecule> moleculeList;
    Complex targCom {};

    // Compartment radius = 100, complex at center.
    ReflectTrajTmpCrds_BuildSingleMoleculeComplex(moleculeList, targCom, Coord { 0.0, 0.0, 0.0 });
    Membrane membrane = ReflectTrajTmpCrds_BuildSphereMembrane(200.0, 100.0);

    std::array<double, 3> traj { 0.5, 0.5, 0.5 };
    const std::array<double, 3> expected = traj;

    std::cerr << "  Compartment radius = " << membrane.compartmentR << "\n";
    std::cerr << "  Input traj = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";

    // isInsideCompartment == true selects the compartment sphere branch.
    reflect_traj_tmp_crds(params, moleculeList, targCom, traj, membrane, 0.0, true);

    std::cerr << "  Output traj = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";

    // Deep inside a big compartment: no reflection expected.
    EXPECT_NEAR(traj[0], expected[0], 1e-9) << "traj[0] changed unexpectedly (inside compartment)";
    EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] changed unexpectedly (inside compartment)";
    EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] changed unexpectedly (inside compartment)";
}

}  // namespace

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each TEST simply invokes the corresponding test_*
// helper so that all assertions are executed and reported, and non-fatal
// EXPECT_* failures do not abort the remaining tests.
// -----------------------------------------------------------------------------

TEST(ReflectTrajTmpCrdsTest, BoxInsideNoChange)
{
    test_reflect_traj_tmp_crds_box_inside_no_change();
}

TEST(ReflectTrajTmpCrdsTest, BoxReflectPosX)
{
    test_reflect_traj_tmp_crds_box_reflect_posX();
}

TEST(ReflectTrajTmpCrdsTest, SphereInsideNoChange)
{
    test_reflect_traj_tmp_crds_sphere_inside_no_change();
}

TEST(ReflectTrajTmpCrdsTest, InsideCompartment)
{
    test_reflect_traj_tmp_crds_inside_compartment();
}