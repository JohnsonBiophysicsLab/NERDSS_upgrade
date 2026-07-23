/*! \file test_reflect_traj_check_span.cpp
 *
 * ### Unit test for ../src/boundary_conditions/reflect_traj_check_span.cpp
 *
 * This test exercises the function:
 *     void reflect_traj_check_span(const Parameters&, Complex&,
 *                                  std::vector<Molecule>&, const Membrane&, double)
 *
 * The function is a simple dispatcher: depending on whether the Membrane
 * describes a sphere or a box, it forwards the call to either
 * reflect_traj_check_span_sphere() or reflect_traj_check_span_box().
 *
 * Because the underlying routines only reflect a complex back inside the
 * simulation volume when it would otherwise exit, we can verify behaviour
 * by:
 *   1. Placing a complex safely inside the volume with a small trajectory
 *      (should be left essentially untouched), and
 *   2. Placing a complex so a trajectory step would push it outside the
 *      volume (should have its trajectory adjusted so it stays inside).
 *
 * Verbose console output is emitted so the reader can follow which source
 * file / function is being tested and what each assertion checks.
 */

#include <array>
#include <iostream>
#include <vector>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Helper: build a minimal single-interface molecule at a given coordinate.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Construct a very small molecule that is part of a complex.
 *
 * The molecule has a single interface located at the same point as its
 * center of mass, which keeps geometry simple for the boundary check.
 *
 * \param[in] com The center-of-mass coordinate to assign.
 * \return A Molecule ready to be placed into a moleculeList.
 */
Molecule rtcs_make_molecule(const Coord& com)
{
    Molecule mol;
    mol.comCoord = com;

    // Give the molecule a single interface coincident with the COM.
    Molecule::Iface iface;
    iface.coord = com;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // The molecule belongs to complex index 0.
    mol.myComIndex = 0;
    return mol;
}

/*! \brief Build a Complex that owns a single member molecule (index 0). */
Complex rtcs_make_complex(const Coord& com)
{
    Complex targCom;
    targCom.comCoord = com;

    // Diffusion constants (used only when the move must be resampled).
    targCom.D.x = 1.0;
    targCom.D.y = 1.0;
    targCom.D.z = 1.0;
    targCom.Dr.x = 0.01;
    targCom.Dr.y = 0.01;
    targCom.Dr.z = 0.01;

    // No initial rotation, no initial translation.
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot.x = 0.0;
    targCom.trajRot.y = 0.0;
    targCom.trajRot.z = 0.0;

    // The single member molecule lives at moleculeList index 0.
    targCom.memberList.clear();
    targCom.memberList.push_back(0);

    // Not restricted to a 2D surface.
    targCom.OnSurface = false;

    return targCom;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Box boundary, complex well inside -> trajectory left essentially
//         unchanged (no reflection required).
// -----------------------------------------------------------------------------
void test_rtcs_box_inside_no_reflection()
{
    std::cerr << "\n[TEST] test_rtcs_box_inside_no_reflection\n"
              << "  Source file:   reflect_traj_check_span.cpp\n"
              << "  Function:      reflect_traj_check_span (box branch)\n"
              << "  Scenario:      complex is well inside the box with a tiny\n"
              << "                 trajectory that keeps it inside.\n"
              << "  Pass criteria: trajTrans is left unchanged (no reflection).\n";

    // Parameters: a small time step is all that is needed.
    Parameters params;
    params.timeStep = 1.0;

    // A cubic water box centered on the origin, side length 100.
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.waterBox.x = 100.0;
    membraneObject.waterBox.y = 100.0;
    membraneObject.waterBox.z = 100.0;

    // Place the complex at the origin - deep inside the box.
    Complex targCom = rtcs_make_complex(Coord{ 0.0, 0.0, 0.0 });
    std::vector<Molecule> moleculeList{ rtcs_make_molecule(Coord{ 0.0, 0.0, 0.0 }) };

    // Give it a tiny displacement that will not reach any wall.
    targCom.trajTrans.x = 1.0;
    targCom.trajTrans.y = -1.0;
    targCom.trajTrans.z = 0.5;

    const double origX = targCom.trajTrans.x;
    const double origY = targCom.trajTrans.y;
    const double origZ = targCom.trajTrans.z;

    const double RS3Dinput = 0.0;

    std::cerr << "  Calling reflect_traj_check_span...\n";
    reflect_traj_check_span(params, targCom, moleculeList, membraneObject, RS3Dinput);

    // Because the complex never approaches a wall, the trajectory must be
    // untouched.
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, origX)
        << "trajTrans.x should be unchanged for an interior complex";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, origY)
        << "trajTrans.y should be unchanged for an interior complex";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, origZ)
        << "trajTrans.z should be unchanged for an interior complex";

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: Box boundary, complex would exit +X wall -> trajectory reflected so
//         the complex stays inside the box.
// -----------------------------------------------------------------------------
void test_rtcs_box_reflects_off_wall()
{
    std::cerr << "\n[TEST] test_rtcs_box_reflects_off_wall\n"
              << "  Source file:   reflect_traj_check_span.cpp\n"
              << "  Function:      reflect_traj_check_span (box branch)\n"
              << "  Scenario:      a step would push the complex past +X wall.\n"
              << "  Pass criteria: after the call the final +X coordinate is\n"
              << "                 no longer outside the box wall.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.waterBox.x = 100.0;
    membraneObject.waterBox.y = 100.0;
    membraneObject.waterBox.z = 100.0;

    // Positive X wall of the box.
    const double posX = membraneObject.waterBox.x / 2.0; // 50

    // Place the complex near the +X wall.
    Complex targCom = rtcs_make_complex(Coord{ 45.0, 0.0, 0.0 });
    std::vector<Molecule> moleculeList{ rtcs_make_molecule(Coord{ 45.0, 0.0, 0.0 }) };

    // A translation that would push it out past +X (45 + 20 = 65 > 50).
    targCom.trajTrans.x = 20.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;

    const double RS3Dinput = 0.0;

    std::cerr << "  Calling reflect_traj_check_span...\n";
    reflect_traj_check_span(params, targCom, moleculeList, membraneObject, RS3Dinput);

    // Compute the final X position after the (possibly reflected) trajectory.
    const double finalX = targCom.comCoord.x + targCom.trajTrans.x;

    std::cerr << "  Final X position after reflection = " << finalX
              << " (wall at " << posX << ")\n";

    // The reflected position should now be inside (or at) the wall, within a
    // small numerical tolerance.
    EXPECT_LE(finalX, posX + 1e-6)
        << "Complex should be reflected back inside the +X wall";
}

// -----------------------------------------------------------------------------
// Test 3: Sphere boundary, complex well inside -> trajectory left unchanged.
// -----------------------------------------------------------------------------
void test_rtcs_sphere_inside_no_reflection()
{
    std::cerr << "\n[TEST] test_rtcs_sphere_inside_no_reflection\n"
              << "  Source file:   reflect_traj_check_span.cpp\n"
              << "  Function:      reflect_traj_check_span (sphere branch)\n"
              << "  Scenario:      complex is deep inside a spherical volume.\n"
              << "  Pass criteria: trajTrans is left unchanged (no reflection).\n";

    Parameters params;
    params.timeStep = 1.0;

    // A spherical boundary of radius 100.
    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Place complex near the center of the sphere.
    Complex targCom = rtcs_make_complex(Coord{ 0.0, 0.0, 0.0 });
    std::vector<Molecule> moleculeList{ rtcs_make_molecule(Coord{ 0.0, 0.0, 0.0 }) };

    // Small displacement that stays well inside the sphere.
    targCom.trajTrans.x = 1.0;
    targCom.trajTrans.y = 1.0;
    targCom.trajTrans.z = 1.0;

    const double origX = targCom.trajTrans.x;
    const double origY = targCom.trajTrans.y;
    const double origZ = targCom.trajTrans.z;

    const double RS3Dinput = 0.0;

    std::cerr << "  Calling reflect_traj_check_span...\n";
    reflect_traj_check_span(params, targCom, moleculeList, membraneObject, RS3Dinput);

    // Interior complexes should not be adjusted.
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, origX)
        << "trajTrans.x should be unchanged for an interior complex (sphere)";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, origY)
        << "trajTrans.y should be unchanged for an interior complex (sphere)";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, origZ)
        << "trajTrans.z should be unchanged for an interior complex (sphere)";

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 4: Sphere boundary, complex would exit the sphere -> trajectory adjusted
//         so the final radial distance is not larger than the sphere radius.
// -----------------------------------------------------------------------------
void test_rtcs_sphere_reflects_off_wall()
{
    std::cerr << "\n[TEST] test_rtcs_sphere_reflects_off_wall\n"
              << "  Source file:   reflect_traj_check_span.cpp\n"
              << "  Function:      reflect_traj_check_span (sphere branch)\n"
              << "  Scenario:      a step would push the complex outside the sphere.\n"
              << "  Pass criteria: the final radial distance is <= sphere radius.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Place the complex near the spherical shell.
    Complex targCom = rtcs_make_complex(Coord{ 95.0, 0.0, 0.0 });
    std::vector<Molecule> moleculeList{ rtcs_make_molecule(Coord{ 95.0, 0.0, 0.0 }) };

    // A displacement that would push it outside the sphere (95 + 20 = 115 > 100).
    targCom.trajTrans.x = 20.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;

    const double RS3Dinput = 0.0;

    std::cerr << "  Calling reflect_traj_check_span...\n";
    reflect_traj_check_span(params, targCom, moleculeList, membraneObject, RS3Dinput);

    // Compute the resulting radial distance from the sphere center.
    const double fx = targCom.comCoord.x + targCom.trajTrans.x;
    const double fy = targCom.comCoord.y + targCom.trajTrans.y;
    const double fz = targCom.comCoord.z + targCom.trajTrans.z;
    const double radial = std::sqrt(fx * fx + fy * fy + fz * fz);

    std::cerr << "  Final radial distance after reflection = " << radial
              << " (sphere radius " << membraneObject.sphereR << ")\n";

    // The reflected position should be inside (or on) the spherical shell,
    // within a small tolerance.
    EXPECT_LE(radial, membraneObject.sphereR + 1e-6)
        << "Complex should be reflected back inside the spherical boundary";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* function is executed inside a TEST so
// the framework reports individual results while still running all of them.
// -----------------------------------------------------------------------------
TEST(ReflectTrajCheckSpan, BoxInsideNoReflection) { test_rtcs_box_inside_no_reflection(); }
TEST(ReflectTrajCheckSpan, BoxReflectsOffWall) { test_rtcs_box_reflects_off_wall(); }
TEST(ReflectTrajCheckSpan, SphereInsideNoReflection) { test_rtcs_sphere_inside_no_reflection(); }
TEST(ReflectTrajCheckSpan, SphereReflectsOffWall) { test_rtcs_sphere_reflects_off_wall(); }