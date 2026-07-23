#include <vector>
#include <iostream>

#include "gtest/gtest.h"

// Include the headers that provide the types and function declarations used by
// the source file under test.  reflect_functions.hpp is expected to declare
// reflect_traj_complex_rad_rot() as well as the helper reflect_* functions.
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

// ---------------------------------------------------------------------------
// Unit tests for: ../src/boundary_conditions/reflect_traj_complex_rad_rot.cpp
//
// The function under test is:
//
//     void reflect_traj_complex_rad_rot(const Parameters& params,
//                                        std::vector<Molecule>& moleculeList,
//                                        Complex& targCom,
//                                        const Membrane& membraneObject,
//                                        double RS3Dinput,
//                                        bool isInsideCompartment);
//
// This function is a dispatcher that decides, based on the geometry of the
// simulation box (sphere vs. box), the compartment enforcement flag, and the
// isInsideCompartment argument, which reflection routine to call in order to
// keep a molecular complex inside its boundaries.
//
// Because the function delegates to helper functions and does not return a
// value, these tests exercise the different control-flow branches and verify
// that the routine can be called safely and that the state it manipulates
// (the complex trajectory translation) remains well-defined (finite).
// ---------------------------------------------------------------------------

namespace {

// ---------------------------------------------------------------------------
// Helper: build a minimal single-molecule complex located inside a cubic box.
// Returns a fully populated (moleculeList, complex) pair through references.
// ---------------------------------------------------------------------------
void reflectTraj_buildSimpleSystem(std::vector<Molecule>& moleculeList,
                                   Complex& targCom,
                                   bool enforceCompartmentBC) {
    // Create a single molecule with one interface.
    Molecule mol;
    mol.index = 0;
    mol.enforceCompartmentBC = enforceCompartmentBC;

    // Place the molecule center-of-mass at the origin.
    mol.comCoord = Coord{0.0, 0.0, 0.0};

    // Give the molecule a single interface coincident with the COM so that
    // rotation/translation math stays trivially bounded.
    Molecule::Iface iface;
    iface.coord = Coord{0.0, 0.0, 0.0};
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    moleculeList.clear();
    moleculeList.push_back(mol);

    // Build the complex referencing that single molecule.
    targCom.index = 0;
    targCom.memberList.clear();
    targCom.memberList.push_back(0);
    targCom.comCoord = Coord{0.0, 0.0, 0.0};
    targCom.trajTrans = Vector {0.0, 0.0, 0.0};
    targCom.radius = 1.0;
}

// ---------------------------------------------------------------------------
// Helper: build a Membrane object configured as a cubic box.
// ---------------------------------------------------------------------------
Membrane reflectTraj_makeBoxMembrane() {
    Membrane membrane;
    membrane.isSphere = false;
    membrane.waterBox.x = 100.0;
    membrane.waterBox.y = 100.0;
    membrane.waterBox.z = 100.0;
    return membrane;
}

// ---------------------------------------------------------------------------
// Helper: build a Membrane object configured as a sphere.
// ---------------------------------------------------------------------------
Membrane reflectTraj_makeSphereMembrane() {
    Membrane membrane;
    membrane.isSphere = true;
    membrane.sphereR = 50.0;
    membrane.compartmentR = 25.0;
    return membrane;
}

// ---------------------------------------------------------------------------
// Test 1: Box geometry, molecule fully inside, no compartment enforcement.
//
// Criteria: the routine must run without crashing and leave the trajectory
// translation finite.  Because the molecule starts well inside the box, we
// do not expect the translation to be modified.
// ---------------------------------------------------------------------------
void test_reflectTraj_BoxInsideNoCompartment() {
    std::cerr << "[RUN] test_reflectTraj_BoxInsideNoCompartment\n"
              << "  Source file : reflect_traj_complex_rad_rot.cpp\n"
              << "  Function    : reflect_traj_complex_rad_rot()\n"
              << "  Scenario    : cubic box, molecule inside, "
              << "no compartment BC, isInsideCompartment=false\n";

    Parameters params;                 // Default-constructed parameters.
    std::vector<Molecule> moleculeList;
    Complex targCom;
    reflectTraj_buildSimpleSystem(moleculeList, targCom,
                                  /*enforceCompartmentBC=*/false);

    Membrane membrane = reflectTraj_makeBoxMembrane();
    const double RS3Dinput = 0.0;

    // Call the function under test; should not throw or crash.
    reflect_traj_complex_rad_rot(params, moleculeList, targCom, membrane,
                                 RS3Dinput, /*isInsideCompartment=*/false);

    // The resulting translation should be a finite (well-defined) value.
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.x))
        << "trajTrans.x is not finite after reflection.";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.y))
        << "trajTrans.y is not finite after reflection.";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.z))
        << "trajTrans.z is not finite after reflection.";

    std::cerr << "  Result: trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n"
              << "[DONE] test_reflectTraj_BoxInsideNoCompartment\n\n";
}

// ---------------------------------------------------------------------------
// Test 2: Box geometry, compartment enforcement ON.
//
// This exercises the extra branch that calls
// reflect_traj_complex_compartment() when the member molecule's
// enforceCompartmentBC flag is set.
//
// Criteria: run without crashing, keep translation finite.
// ---------------------------------------------------------------------------
void test_reflectTraj_BoxWithCompartmentBC() {
    std::cerr << "[RUN] test_reflectTraj_BoxWithCompartmentBC\n"
              << "  Source file : reflect_traj_complex_rad_rot.cpp\n"
              << "  Function    : reflect_traj_complex_rad_rot()\n"
              << "  Scenario    : cubic box, compartment BC enforced\n";

    Parameters params;
    std::vector<Molecule> moleculeList;
    Complex targCom;
    reflectTraj_buildSimpleSystem(moleculeList, targCom,
                                  /*enforceCompartmentBC=*/true);

    Membrane membrane = reflectTraj_makeBoxMembrane();
    membrane.compartmentR = 25.0;  // Needed by the compartment routine.
    const double RS3Dinput = 0.0;

    reflect_traj_complex_rad_rot(params, moleculeList, targCom, membrane,
                                 RS3Dinput, /*isInsideCompartment=*/false);

    // Verify the translation remains finite (no NaN/inf produced).
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.x))
        << "trajTrans.x is not finite after compartment reflection.";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.y))
        << "trajTrans.y is not finite after compartment reflection.";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.z))
        << "trajTrans.z is not finite after compartment reflection.";

    std::cerr << "  Result: trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n"
              << "[DONE] test_reflectTraj_BoxWithCompartmentBC\n\n";
}

// ---------------------------------------------------------------------------
// Test 3: Sphere geometry (isSphere == true), outside compartment.
//
// This exercises the branch that calls
// reflect_traj_complex_rad_rot_sphere() using membrane.sphereR.
//
// Criteria: run without crashing, keep translation finite.
// ---------------------------------------------------------------------------
void test_reflectTraj_SphereOutsideCompartment() {
    std::cerr << "[RUN] test_reflectTraj_SphereOutsideCompartment\n"
              << "  Source file : reflect_traj_complex_rad_rot.cpp\n"
              << "  Function    : reflect_traj_complex_rad_rot()\n"
              << "  Scenario    : spherical boundary, "
              << "isInsideCompartment=false\n";

    Parameters params;
    std::vector<Molecule> moleculeList;
    Complex targCom;
    reflectTraj_buildSimpleSystem(moleculeList, targCom,
                                  /*enforceCompartmentBC=*/false);

    Membrane membrane = reflectTraj_makeSphereMembrane();
    const double RS3Dinput = 0.0;

    reflect_traj_complex_rad_rot(params, moleculeList, targCom, membrane,
                                 RS3Dinput, /*isInsideCompartment=*/false);

    EXPECT_TRUE(std::isfinite(targCom.trajTrans.x))
        << "trajTrans.x is not finite (sphere, outside compartment).";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.y))
        << "trajTrans.y is not finite (sphere, outside compartment).";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.z))
        << "trajTrans.z is not finite (sphere, outside compartment).";

    std::cerr << "  Result: trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n"
              << "[DONE] test_reflectTraj_SphereOutsideCompartment\n\n";
}

// ---------------------------------------------------------------------------
// Test 4: isInsideCompartment == true.
//
// This exercises the else-branch which always calls
// reflect_traj_complex_rad_rot_sphere() with membrane.compartmentR.
//
// Criteria: run without crashing, keep translation finite.
// ---------------------------------------------------------------------------
void test_reflectTraj_InsideCompartment() {
    std::cerr << "[RUN] test_reflectTraj_InsideCompartment\n"
              << "  Source file : reflect_traj_complex_rad_rot.cpp\n"
              << "  Function    : reflect_traj_complex_rad_rot()\n"
              << "  Scenario    : isInsideCompartment=true, "
              << "reflect against compartmentR\n";

    Parameters params;
    std::vector<Molecule> moleculeList;
    Complex targCom;
    reflectTraj_buildSimpleSystem(moleculeList, targCom,
                                  /*enforceCompartmentBC=*/false);

    Membrane membrane = reflectTraj_makeSphereMembrane();
    const double RS3Dinput = 0.0;

    reflect_traj_complex_rad_rot(params, moleculeList, targCom, membrane,
                                 RS3Dinput, /*isInsideCompartment=*/true);

    EXPECT_TRUE(std::isfinite(targCom.trajTrans.x))
        << "trajTrans.x is not finite (inside compartment).";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.y))
        << "trajTrans.y is not finite (inside compartment).";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.z))
        << "trajTrans.z is not finite (inside compartment).";

    std::cerr << "  Result: trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n"
              << "[DONE] test_reflectTraj_InsideCompartment\n\n";
}

}  // namespace

// ---------------------------------------------------------------------------
// GoogleTest wrappers.  Each TEST simply forwards to the corresponding
// test_* helper defined above.  Grouping the logic in helpers keeps the
// verbose console output separate from the framework registration.
// ---------------------------------------------------------------------------

// Verifies dispatch for a cubic box with no compartment enforcement.
TEST(ReflectTrajComplexRadRotTest, BoxInsideNoCompartment) {
    test_reflectTraj_BoxInsideNoCompartment();
}

// Verifies dispatch for a cubic box with compartment enforcement.
TEST(ReflectTrajComplexRadRotTest, BoxWithCompartmentBC) {
    test_reflectTraj_BoxWithCompartmentBC();
}

// Verifies dispatch for a spherical boundary, outside the compartment.
TEST(ReflectTrajComplexRadRotTest, SphereOutsideCompartment) {
    test_reflectTraj_SphereOutsideCompartment();
}

// Verifies dispatch for the isInsideCompartment==true branch.
TEST(ReflectTrajComplexRadRotTest, InsideCompartment) {
    test_reflectTraj_InsideCompartment();
}