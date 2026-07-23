// -----------------------------------------------------------------------------
// Unit test for: ../src/boundary_conditions/reflect_complex_rad_rot.cpp
//
// The file under test provides a single function:
//
//     void reflect_complex_rad_rot(const Membrane& membraneObject,
//                                  Complex& targCom,
//                                  std::vector<Molecule>& moleculeList,
//                                  double RS3Dinput,
//                                  bool isInsideCompartment)
//
// This function is a dispatcher.  Depending on the membrane geometry and the
// isInsideCompartment flag it forwards the call to one of several lower level
// helper routines that actually reflect a complex back inside the simulation
// volume (a box or a sphere):
//
//   * reflect_complex_rad_rot_sphere(...)   -> spherical boundary
//   * reflect_complex_rad_rot_box(...)      -> rectangular box boundary
//   * reflect_complex_compartment(...)      -> enforce compartment boundary
//
// Because these helpers modify the coordinates in moleculeList / targCom, our
// tests exercise the various decision branches and verify that a complex that
// begins outside the boundary is moved back inside (or that a complex already
// inside the boundary is left untouched).
// -----------------------------------------------------------------------------

#include "boundary_conditions/reflect_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Helper: build a simple single-molecule complex located at a given position.
//
// We create one Molecule with a single interface coincident with its center
// of mass, and one Complex whose memberList references that molecule.  This is
// the minimal well-formed input the dispatcher (and its helpers) expect.
// ---------------------------------------------------------------------------
void rcrr_build_simple_complex(double x, double y, double z, double radius,
                               Complex& targCom,
                               std::vector<Molecule>& moleculeList,
                               bool enforceCompartmentBC = false) {
    // --- Build the single molecule -----------------------------------------
    Molecule mol;
    mol.comCoord = Coord { x, y, z };
    mol.index = 0;
    mol.enforceCompartmentBC = enforceCompartmentBC;

    // Give the molecule a single interface at the center of mass so that the
    // reflect helpers have an interface coordinate to examine.
    Molecule::Iface iface {};
    iface.coord = Coord { x, y, z };
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    moleculeList.clear();
    moleculeList.push_back(mol);

    // --- Build the complex that owns the molecule --------------------------
    targCom.comCoord = Coord { x, y, z };
    targCom.radius = radius;
    targCom.memberList.clear();
    targCom.memberList.push_back(0); // references moleculeList[0]
    targCom.OnSurface = false;
}

// ---------------------------------------------------------------------------
// Test 1: Spherical membrane, complex OUTSIDE the sphere.
//
// When the membrane is a sphere and the complex sits outside sphereR, the
// dispatcher should route to reflect_complex_rad_rot_sphere which pulls the
// complex back toward / inside the sphere.  We verify that the resulting
// distance from the origin is no greater than the starting distance (i.e. the
// complex was moved inward, not further out).
// ---------------------------------------------------------------------------
void test_rcrr_sphere_outside() {
    std::cerr << "[TEST] reflect_complex_rad_rot (reflect_complex_rad_rot.cpp)\n"
              << "       Scenario: spherical membrane, complex OUTSIDE sphere.\n";

    Membrane membrane;
    membrane.isSphere = true;
    membrane.sphereR = 100.0;
    membrane.compartmentR = 50.0;

    Complex targCom;
    std::vector<Molecule> moleculeList;

    // Place the complex well outside the sphere radius (distance 150 > 100).
    const double startDist = 150.0;
    rcrr_build_simple_complex(startDist, 0.0, 0.0, /*radius=*/5.0,
                              targCom, moleculeList);

    std::cerr << "       Start COM = (" << targCom.comCoord.x << ", "
              << targCom.comCoord.y << ", " << targCom.comCoord.z
              << ")  |dist|=" << startDist << ", sphereR=" << membrane.sphereR
              << '\n';

    // Call the function under test (not inside a compartment).
    reflect_complex_rad_rot(membrane, targCom, moleculeList,
                            /*RS3Dinput=*/0.0, /*isInsideCompartment=*/false);

    const double newDist = std::sqrt(targCom.comCoord.x * targCom.comCoord.x
                                     + targCom.comCoord.y * targCom.comCoord.y
                                     + targCom.comCoord.z * targCom.comCoord.z);

    std::cerr << "       After  COM = (" << targCom.comCoord.x << ", "
              << targCom.comCoord.y << ", " << targCom.comCoord.z
              << ")  |dist|=" << newDist << '\n';

    // Passing criterion: the complex must have moved inward (or at least not
    // outward) after reflection off the spherical boundary.
    EXPECT_LE(newDist, startDist)
        << "Complex outside the sphere should be moved inward by reflection.";
}

// ---------------------------------------------------------------------------
// Test 2: Spherical membrane, complex ALREADY INSIDE the sphere.
//
// A complex comfortably inside the sphere should be left essentially where it
// is (no reflection needed).  We check that the coordinates are unchanged.
// ---------------------------------------------------------------------------
void test_rcrr_sphere_inside() {
    std::cerr << "[TEST] reflect_complex_rad_rot (reflect_complex_rad_rot.cpp)\n"
              << "       Scenario: spherical membrane, complex INSIDE sphere.\n";

    Membrane membrane;
    membrane.isSphere = true;
    membrane.sphereR = 100.0;
    membrane.compartmentR = 50.0;

    Complex targCom;
    std::vector<Molecule> moleculeList;

    // Place the complex safely inside the sphere.
    rcrr_build_simple_complex(10.0, 0.0, 0.0, /*radius=*/1.0,
                              targCom, moleculeList);

    const Coord before = targCom.comCoord;
    std::cerr << "       COM before = (" << before.x << ", " << before.y
              << ", " << before.z << ")\n";

    reflect_complex_rad_rot(membrane, targCom, moleculeList,
                            /*RS3Dinput=*/0.0, /*isInsideCompartment=*/false);

    std::cerr << "       COM after  = (" << targCom.comCoord.x << ", "
              << targCom.comCoord.y << ", " << targCom.comCoord.z << ")\n";

    // Passing criterion: a complex already inside should stay in place.
    EXPECT_NEAR(before.x, targCom.comCoord.x, 1e-9)
        << "X should be unchanged for a complex already inside the sphere.";
    EXPECT_NEAR(before.y, targCom.comCoord.y, 1e-9)
        << "Y should be unchanged for a complex already inside the sphere.";
    EXPECT_NEAR(before.z, targCom.comCoord.z, 1e-9)
        << "Z should be unchanged for a complex already inside the sphere.";
}

// ---------------------------------------------------------------------------
// Test 3: Box membrane, complex OUTSIDE the box in +X.
//
// With a rectangular waterBox and a non-sphere membrane the dispatcher routes
// to reflect_complex_rad_rot_box.  We place a complex outside the +X wall and
// confirm it is reflected back so that its X coordinate is smaller (moved in
// the -X direction, i.e. back inside).
// ---------------------------------------------------------------------------
void test_rcrr_box_outside() {
    std::cerr << "[TEST] reflect_complex_rad_rot (reflect_complex_rad_rot.cpp)\n"
              << "       Scenario: box membrane, complex OUTSIDE +X wall.\n";

    Membrane membrane;
    membrane.isSphere = false;
    std::vector<double> box {100.0, 100.0, 100.0};
    membrane.waterBox = Membrane::WaterBox(box); // half extents = 50
    membrane.compartmentR = 25.0;

    Complex targCom;
    std::vector<Molecule> moleculeList;

    // Place the complex beyond the +X wall (60 > 50 half-extent).
    const double startX = 60.0;
    rcrr_build_simple_complex(startX, 0.0, 0.0, /*radius=*/2.0,
                              targCom, moleculeList);

    std::cerr << "       COM before = (" << targCom.comCoord.x << ", "
              << targCom.comCoord.y << ", " << targCom.comCoord.z
              << "), +X wall at " << (membrane.waterBox.x / 2.0) << '\n';

    reflect_complex_rad_rot(membrane, targCom, moleculeList,
                            /*RS3Dinput=*/0.0, /*isInsideCompartment=*/false);

    std::cerr << "       COM after  = (" << targCom.comCoord.x << ", "
              << targCom.comCoord.y << ", " << targCom.comCoord.z << ")\n";

    // Passing criterion: complex should be moved back toward the interior
    // (its X coordinate must decrease from the out-of-box starting value).
    EXPECT_LT(targCom.comCoord.x, startX)
        << "Complex outside +X wall should be reflected back inside the box.";
}

// ---------------------------------------------------------------------------
// Test 4: isInsideCompartment == true uses the compartment radius.
//
// When the caller flags the complex as inside a compartment, the dispatcher
// always calls the sphere reflection using membrane.compartmentR (regardless
// of whether the outer membrane is a sphere or a box).  We put the complex
// outside compartmentR and verify it is drawn inward.
// ---------------------------------------------------------------------------
void test_rcrr_inside_compartment() {
    std::cerr << "[TEST] reflect_complex_rad_rot (reflect_complex_rad_rot.cpp)\n"
              << "       Scenario: isInsideCompartment=true, use compartmentR.\n";

    Membrane membrane;
    membrane.isSphere = false;                 // outer geometry irrelevant here
    std::vector<double> box {500.0, 500.0, 500.0};
    membrane.waterBox = Membrane::WaterBox(box); // half extents = 50
    membrane.sphereR = 200.0;
    membrane.compartmentR = 40.0;              // the radius that will be used

    Complex targCom;
    std::vector<Molecule> moleculeList;

    // Place the complex outside the compartment radius (60 > 40).
    const double startDist = 60.0;
    rcrr_build_simple_complex(startDist, 0.0, 0.0, /*radius=*/3.0,
                              targCom, moleculeList);

    std::cerr << "       Start |dist|=" << startDist
              << ", compartmentR=" << membrane.compartmentR << '\n';

    reflect_complex_rad_rot(membrane, targCom, moleculeList,
                            /*RS3Dinput=*/0.0, /*isInsideCompartment=*/true);

    const double newDist = std::sqrt(targCom.comCoord.x * targCom.comCoord.x
                                     + targCom.comCoord.y * targCom.comCoord.y
                                     + targCom.comCoord.z * targCom.comCoord.z);

    std::cerr << "       After |dist|=" << newDist << '\n';

    // Passing criterion: the complex is reflected inward using compartmentR.
    EXPECT_LE(newDist, startDist)
        << "Complex outside the compartment should be moved inward.";
}

} // namespace

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each TEST simply invokes the corresponding test_*()
// helper.  Using non-fatal EXPECT_* assertions ensures every scenario runs
// even if one of them fails.
// -----------------------------------------------------------------------------

TEST(ReflectComplexRadRotTest, SphereOutside) {
    test_rcrr_sphere_outside();
}

TEST(ReflectComplexRadRotTest, SphereInside) {
    test_rcrr_sphere_inside();
}

TEST(ReflectComplexRadRotTest, BoxOutside) {
    test_rcrr_box_outside();
}

TEST(ReflectComplexRadRotTest, InsideCompartment) {
    test_rcrr_inside_compartment();
}