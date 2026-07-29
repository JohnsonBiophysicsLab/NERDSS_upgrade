/*! \file test_reflect_traj_complex_compartment.cpp
 *
 * ### Unit test for src/boundary_conditions/reflect_traj_complex_compartment.cpp
 *
 * Function under test:
 *
 *     void reflect_traj_complex_compartment(const Parameters& params,
 *                                           std::vector<Molecule>& moleculeList,
 *                                           Complex& targCom,
 *                                           const Membrane& membraneObject,
 *                                           double RS3Dinput)
 *
 * Behaviour of the routine (as implemented):
 *
 *   1. The "compartment" is a sphere centred on the origin whose radius is
 *      `membraneObject.compartmentR + RS3D`, where RS3D == RS3Dinput unless the
 *      complex sits on the membrane surface (`targCom.OnSurface == true`), in
 *      which case RS3D is forced to zero.
 *   2. A guard is applied first: the routine only does work when the
 *      *translated* complex centre of mass plus the complex bounding radius lies
 *      completely inside the compartment sphere
 *      (`|comCoord + trajTrans| + radius < sphereR`).
 *   3. If that guard passes, every member molecule centre of mass and every
 *      interface is propagated (translation + Euler rotation from
 *      `targCom.trajRot`) and the *smallest* radial distance to the origin is
 *      recorded. If that smallest distance is inside the compartment the
 *      translation vector `trajTrans` is augmented radially outwards by
 *      `lamda * targcrds`, with `lamda = -2 (targR - sphereR)/targR`, which
 *      mirrors the offending point through the compartment shell so that its new
 *      radial distance becomes exactly `2*sphereR - targR`.
 *   4. Only `targCom.trajTrans` is modified — no coordinates are committed.
 *
 * Each test below prints what is being exercised and what the pass criteria are.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (unique rtcc_ prefix so they cannot collide with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal Molecule with a COM and an arbitrary set of interfaces.
 *
 * \param[in] com          Centre of mass of the molecule.
 * \param[in] ifaceCoords  Absolute coordinates for each interface.
 * \return A Molecule owned by complex index 0.
 */
Molecule rtcc_make_molecule(const Coord& com, const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.comCoord = com;
    mol.interfaceList.clear();
    for (const auto& oneIfaceCrd : ifaceCoords) {
        Molecule::Iface iface;
        iface.coord = oneIfaceCrd;
        mol.interfaceList.push_back(iface);
    }
    // Every molecule created here belongs to the single complex we build below.
    mol.myComIndex = 0;
    return mol;
}

/*! \brief Build a Complex holding exactly one member molecule (moleculeList[0]).
 *
 * Translation and rotation trajectories start out at zero so each test can set
 * them explicitly.
 *
 * \param[in] com     Complex centre of mass.
 * \param[in] radius  Bounding radius of the complex (used by the routine's guard).
 */
Complex rtcc_make_complex(const Coord& com, double radius)
{
    Complex targCom;
    targCom.comCoord = com;
    targCom.radius = radius;

    // Diffusion constants are not used by this routine (the resample branch is
    // commented out in the source) but are set to sane values anyway.
    targCom.D = Coord(1.0, 1.0, 1.0);
    targCom.Dr = Coord(0.01, 0.01, 0.01);

    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot.x = 0.0;
    targCom.trajRot.y = 0.0;
    targCom.trajRot.z = 0.0;

    targCom.memberList.clear();
    targCom.memberList.push_back(0);

    targCom.OnSurface = false;
    return targCom;
}

/*! \brief Radial distance from the origin of a point after the complex move.
 *
 * Mirrors exactly what the routine does: rotate the (point - COM) vector with
 * the Euler matrix built from `trajRot`, then add COM + trajTrans.
 */
double rtcc_radial_after_move(const Complex& targCom, const std::array<double, 9>& M, const Coord& pt)
{
    const double vx = pt.x - targCom.comCoord.x;
    const double vy = pt.y - targCom.comCoord.y;
    const double vz = pt.z - targCom.comCoord.z;

    const double rx = M[0] * vx + M[1] * vy + M[2] * vz;
    const double ry = M[3] * vx + M[4] * vy + M[5] * vz;
    const double rz = M[6] * vx + M[7] * vy + M[8] * vz;

    const double fx = targCom.comCoord.x + targCom.trajTrans.x + rx;
    const double fy = targCom.comCoord.y + targCom.trajTrans.y + ry;
    const double fz = targCom.comCoord.z + targCom.trajTrans.z + rz;

    return std::sqrt(fx * fx + fy * fy + fz * fz);
}

/*! \brief Smallest radial distance to the origin over all molecule COMs and
 *         interfaces of the complex, after applying the current trajectory.
 *
 * This is exactly the quantity the routine calls `targR`.
 */
double rtcc_min_radial_distance(const Complex& targCom, const std::vector<Molecule>& moleculeList)
{
    std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);

    double minR = std::numeric_limits<double>::max();
    for (auto memMol : targCom.memberList) {
        const Molecule& mol = moleculeList[memMol];
        minR = std::min(minR, rtcc_radial_after_move(targCom, M, mol.comCoord));
        for (const auto& iface : mol.interfaceList)
            minR = std::min(minR, rtcc_radial_after_move(targCom, M, iface.coord));
    }
    return minR;
}

/*! \brief Radial distance of the (untranslated, unrotated) complex COM. */
double rtcc_magnitude(const Coord& c) { return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z); }

} // namespace

// -----------------------------------------------------------------------------
// Test 1: The guard `|com + trajTrans| + radius < sphereR` fails.
//         Nothing at all should be touched.
// -----------------------------------------------------------------------------
void test_rtcc_guard_blocks_reflection()
{
    std::cerr << "\n[TEST] test_rtcc_guard_blocks_reflection\n"
              << "  Source file:   reflect_traj_complex_compartment.cpp\n"
              << "  Function:      reflect_traj_complex_compartment\n"
              << "  Scenario:      the complex bounding sphere is NOT fully enclosed\n"
              << "                 by the compartment (|com|+radius >= compartmentR),\n"
              << "                 so the routine's entry guard must reject it.\n"
              << "  Pass criteria: trajTrans, trajRot and molecule coords unchanged.\n";

    Parameters params;
    params.timeStep = 1.0;

    // Compartment sphere of radius 10 centred on the origin.
    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 10.0;

    // |com| = 9, radius = 2  ->  9 + 2 = 11 which is NOT < 10.
    Complex targCom = rtcc_make_complex(Coord { 9.0, 0.0, 0.0 }, 2.0);
    std::vector<Molecule> moleculeList { rtcc_make_molecule(Coord { 9.0, 0.0, 0.0 },
        { Coord { 9.0, 0.0, 0.0 } }) };

    std::cerr << "  |com| = " << rtcc_magnitude(targCom.comCoord)
              << ", complex radius = " << targCom.radius
              << ", compartment radius = " << membraneObject.compartmentR << '\n';

    reflect_traj_complex_compartment(params, moleculeList, targCom, membraneObject, 0.0);

    // The trajectory must not have been modified at all.
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 0.0) << "trajTrans.x must stay 0 when the guard rejects";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "trajTrans.y must stay 0 when the guard rejects";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0) << "trajTrans.z must stay 0 when the guard rejects";

    // The routine never rotates or commits coordinates.
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, 0.0) << "trajRot must never be modified";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, 0.0) << "trajRot must never be modified";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, 0.0) << "trajRot must never be modified";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 9.0) << "molecule coordinates must not be committed";

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: Molecule COM sits inside the compartment -> mirrored outwards.
//         The arithmetic is fully deterministic, so exact values are checked.
// -----------------------------------------------------------------------------
void test_rtcc_reflects_molecule_out_of_compartment()
{
    std::cerr << "\n[TEST] test_rtcc_reflects_molecule_out_of_compartment\n"
              << "  Source file:   reflect_traj_complex_compartment.cpp\n"
              << "  Function:      reflect_traj_complex_compartment\n"
              << "  Scenario:      single small molecule at radius 2 inside a\n"
              << "                 compartment of radius 10, zero rotation.\n"
              << "  Pass criteria: lamda = -2(2-10)/2 = 8, so trajTrans becomes\n"
              << "                 8*(2,0,0) = (16,0,0) and the final radial\n"
              << "                 distance is 2*10 - 2 = 18 (mirrored shell).\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 10.0;

    // Complex COM at (2,0,0), tiny bounding radius so the guard passes (3 < 10).
    Complex targCom = rtcc_make_complex(Coord { 2.0, 0.0, 0.0 }, 1.0);
    std::vector<Molecule> moleculeList { rtcc_make_molecule(Coord { 2.0, 0.0, 0.0 },
        { Coord { 2.0, 0.0, 0.0 } }) };

    reflect_traj_complex_compartment(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";

    // Exact expected translation: purely radial, along +x.
    EXPECT_NEAR(targCom.trajTrans.x, 16.0, 1e-12) << "trajTrans.x should be lamda*targcrds.x = 8*2 = 16";
    EXPECT_NEAR(targCom.trajTrans.y, 0.0, 1e-12) << "no y displacement is expected";
    EXPECT_NEAR(targCom.trajTrans.z, 0.0, 1e-12) << "no z displacement is expected";

    // The offending point is mirrored through the compartment shell.
    const double finalRadial = rtcc_min_radial_distance(targCom, moleculeList);
    std::cerr << "  Final minimum radial distance = " << finalRadial
              << " (expected 2*10 - 2 = 18, must be >= compartment radius)\n";
    EXPECT_NEAR(finalRadial, 18.0, 1e-9) << "reflected point should sit at 2*sphereR - targR";
    EXPECT_GE(finalRadial, membraneObject.compartmentR)
        << "after reflection nothing may remain inside the compartment";
}

// -----------------------------------------------------------------------------
// Test 3: A pre-existing trajTrans must be *added to*, not overwritten.
// -----------------------------------------------------------------------------
void test_rtcc_accumulates_existing_translation()
{
    std::cerr << "\n[TEST] test_rtcc_accumulates_existing_translation\n"
              << "  Source file:   reflect_traj_complex_compartment.cpp\n"
              << "  Function:      reflect_traj_complex_compartment\n"
              << "  Scenario:      COM (0,0,3) with an existing step (0,0,1) inside a\n"
              << "                 compartment of radius 10; targR = 4.\n"
              << "  Pass criteria: lamda = -2(4-10)/4 = 3 so trajTrans becomes\n"
              << "                 (0,0,1) + 3*(0,0,4) = (0,0,13) and the final\n"
              << "                 radial distance is 2*10 - 4 = 16.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 10.0;

    Complex targCom = rtcc_make_complex(Coord { 0.0, 0.0, 3.0 }, 1.0);
    std::vector<Molecule> moleculeList { rtcc_make_molecule(Coord { 0.0, 0.0, 3.0 },
        { Coord { 0.0, 0.0, 3.0 } }) };

    // Pre-existing propagation step along +z.
    targCom.trajTrans.z = 1.0;

    reflect_traj_complex_compartment(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";

    EXPECT_NEAR(targCom.trajTrans.z, 13.0, 1e-12)
        << "the radial correction must be added onto the incoming translation";
    EXPECT_NEAR(targCom.trajTrans.x, 0.0, 1e-12) << "no x displacement is expected";
    EXPECT_NEAR(targCom.trajTrans.y, 0.0, 1e-12) << "no y displacement is expected";

    const double finalRadial = rtcc_min_radial_distance(targCom, moleculeList);
    std::cerr << "  Final minimum radial distance = " << finalRadial << " (expected 16)\n";
    EXPECT_NEAR(finalRadial, 16.0, 1e-9) << "reflected point should sit at 2*sphereR - targR";
}

// -----------------------------------------------------------------------------
// Test 4: An interface deeper inside the compartment than the molecule COM must
//         be the point that drives the reflection (targR is a minimum).
// -----------------------------------------------------------------------------
void test_rtcc_interface_drives_reflection()
{
    std::cerr << "\n[TEST] test_rtcc_interface_drives_reflection\n"
              << "  Source file:   reflect_traj_complex_compartment.cpp\n"
              << "  Function:      reflect_traj_complex_compartment\n"
              << "  Scenario:      molecule COM at radius 5 but one interface at\n"
              << "                 radius 2; compartment radius 10.\n"
              << "  Pass criteria: the interface (smallest radius) sets targR = 2,\n"
              << "                 giving trajTrans = (16,0,0) and a mirrored\n"
              << "                 interface radius of 2*10 - 2 = 18.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 10.0;

    // Guard: |com| + radius = 5 + 1 = 6 < 10, so the routine proceeds.
    Complex targCom = rtcc_make_complex(Coord { 5.0, 0.0, 0.0 }, 1.0);
    std::vector<Molecule> moleculeList { rtcc_make_molecule(Coord { 5.0, 0.0, 0.0 },
        { Coord { 2.0, 0.0, 0.0 } }) };

    const double beforeMin = rtcc_min_radial_distance(targCom, moleculeList);
    std::cerr << "  Minimum radial distance before the call = " << beforeMin
              << " (should be the interface at 2)\n";
    EXPECT_NEAR(beforeMin, 2.0, 1e-12) << "the interface must be the innermost point";

    reflect_traj_complex_compartment(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";

    EXPECT_NEAR(targCom.trajTrans.x, 16.0, 1e-12)
        << "reflection must use the interface radius (2), not the molecule COM radius (5)";

    const double afterMin = rtcc_min_radial_distance(targCom, moleculeList);
    std::cerr << "  Minimum radial distance after the call = " << afterMin << " (expected 18)\n";
    EXPECT_NEAR(afterMin, 18.0, 1e-9) << "innermost point must be mirrored to 2*sphereR - targR";
    EXPECT_GE(afterMin, membraneObject.compartmentR) << "nothing may remain inside the compartment";
}

// -----------------------------------------------------------------------------
// Test 5: If every propagated point is already outside the compartment shell no
//         reflection is performed (inside stays false).
// -----------------------------------------------------------------------------
void test_rtcc_no_reflection_when_all_points_outside()
{
    std::cerr << "\n[TEST] test_rtcc_no_reflection_when_all_points_outside\n"
              << "  Source file:   reflect_traj_complex_compartment.cpp\n"
              << "  Function:      reflect_traj_complex_compartment\n"
              << "  Scenario:      complex COM is inside the guard region but the only\n"
              << "                 member molecule (and its interface) sit far outside\n"
              << "                 the compartment shell.\n"
              << "  Pass criteria: trajTrans is left unchanged (no point inside).\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 10.0;

    // Guard passes: |com| + radius = 1 + 0.1 < 10.
    Complex targCom = rtcc_make_complex(Coord { 1.0, 0.0, 0.0 }, 0.1);
    // ... but the member molecule is displaced far beyond the compartment shell.
    std::vector<Molecule> moleculeList { rtcc_make_molecule(Coord { 20.0, 0.0, 0.0 },
        { Coord { 21.0, 0.0, 0.0 } }) };

    const double beforeMin = rtcc_min_radial_distance(targCom, moleculeList);
    std::cerr << "  Minimum radial distance before the call = " << beforeMin
              << " (compartment radius " << membraneObject.compartmentR << ")\n";
    EXPECT_GT(beforeMin, membraneObject.compartmentR)
        << "test precondition: all points must start outside the compartment";

    reflect_traj_complex_compartment(params, moleculeList, targCom, membraneObject, 0.0);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 0.0) << "trajTrans.x must be untouched";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "trajTrans.y must be untouched";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0) << "trajTrans.z must be untouched";

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 6: RS3Dinput enlarges the effective compartment radius, but only for
//         complexes that are NOT on the membrane surface.
// -----------------------------------------------------------------------------
void test_rtcc_onsurface_ignores_RS3D()
{
    std::cerr << "\n[TEST] test_rtcc_onsurface_ignores_RS3D\n"
              << "  Source file:   reflect_traj_complex_compartment.cpp\n"
              << "  Function:      reflect_traj_complex_compartment\n"
              << "  Scenario:      identical geometry (point at radius 2, compartmentR\n"
              << "                 = 10, RS3Dinput = 5) run twice: once with\n"
              << "                 OnSurface = false (sphereR = 15) and once with\n"
              << "                 OnSurface = true (RS3D forced to 0 -> sphereR = 10).\n"
              << "  Pass criteria: trajTrans.x = 26 (final radius 28) when off the\n"
              << "                 surface, and 16 (final radius 18) when on it.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 10.0;

    const double RS3Dinput = 5.0;

    // ---- Case A: complex in solution, RS3D is honoured -> sphereR = 15. -----
    Complex comOff = rtcc_make_complex(Coord { 2.0, 0.0, 0.0 }, 1.0);
    comOff.OnSurface = false;
    std::vector<Molecule> molListOff { rtcc_make_molecule(Coord { 2.0, 0.0, 0.0 },
        { Coord { 2.0, 0.0, 0.0 } }) };

    reflect_traj_complex_compartment(params, molListOff, comOff, membraneObject, RS3Dinput);

    std::cerr << "  OnSurface=false: trajTrans.x = " << comOff.trajTrans.x
              << " (expected 26, mirrored radius 28)\n";
    EXPECT_NEAR(comOff.trajTrans.x, 26.0, 1e-12)
        << "RS3Dinput should enlarge the compartment radius to 15 for a solution complex";
    EXPECT_NEAR(rtcc_min_radial_distance(comOff, molListOff), 28.0, 1e-9)
        << "final radius should be 2*15 - 2 = 28";

    // ---- Case B: complex on the surface, RS3D forced to 0 -> sphereR = 10. --
    Complex comOn = rtcc_make_complex(Coord { 2.0, 0.0, 0.0 }, 1.0);
    comOn.OnSurface = true;
    std::vector<Molecule> molListOn { rtcc_make_molecule(Coord { 2.0, 0.0, 0.0 },
        { Coord { 2.0, 0.0, 0.0 } }) };

    reflect_traj_complex_compartment(params, molListOn, comOn, membraneObject, RS3Dinput);

    std::cerr << "  OnSurface=true : trajTrans.x = " << comOn.trajTrans.x
              << " (expected 16, mirrored radius 18)\n";
    EXPECT_NEAR(comOn.trajTrans.x, 16.0, 1e-12)
        << "RS3Dinput must be ignored for a complex flagged OnSurface";
    EXPECT_NEAR(rtcc_min_radial_distance(comOn, molListOn), 18.0, 1e-9)
        << "final radius should be 2*10 - 2 = 18";

    // The two branches must be measurably different, proving RS3D was applied.
    EXPECT_GT(comOff.trajTrans.x, comOn.trajTrans.x)
        << "the solution complex must be pushed further than the surface complex";
}

// -----------------------------------------------------------------------------
// Test 7: The rotation stored in trajRot must be taken into account when
//         locating the innermost point.
// -----------------------------------------------------------------------------
void test_rtcc_accounts_for_rotation()
{
    std::cerr << "\n[TEST] test_rtcc_accounts_for_rotation\n"
              << "  Source file:   reflect_traj_complex_compartment.cpp\n"
              << "  Function:      reflect_traj_complex_compartment\n"
              << "  Scenario:      complex COM at (3,0,0) with a member molecule offset\n"
              << "                 by (1,0,0) and a 90 degree rotation about z, so the\n"
              << "                 rotated molecule sits at radius sqrt(10).\n"
              << "  Pass criteria: the innermost rotated point starts inside the\n"
              << "                 compartment and ends up mirrored at exactly\n"
              << "                 2*compartmentR - targR (>= compartmentR).\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 10.0;

    // Guard: |com| + radius = 3 + 1.5 = 4.5 < 10, so the routine proceeds.
    Complex targCom = rtcc_make_complex(Coord { 3.0, 0.0, 0.0 }, 1.5);
    std::vector<Molecule> moleculeList { rtcc_make_molecule(Coord { 4.0, 0.0, 0.0 },
        { Coord { 4.0, 0.0, 0.0 } }) };

    // Rotate a quarter turn about z; the (1,0,0) offset maps into the y axis.
    targCom.trajRot.z = M_PI / 2.0;

    // Reference value computed with the very same rotation matrix helper, so the
    // expectation does not depend on the Euler convention used internally.
    const double beforeMin = rtcc_min_radial_distance(targCom, moleculeList);
    std::cerr << "  Minimum rotated radial distance before the call = " << beforeMin << '\n';
    EXPECT_LT(beforeMin, membraneObject.compartmentR)
        << "test precondition: the rotated point must start inside the compartment";

    reflect_traj_complex_compartment(params, moleculeList, targCom, membraneObject, 0.0);

    // trajRot itself must never be modified by this routine.
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, M_PI / 2.0) << "trajRot must be left alone";

    const double afterMin = rtcc_min_radial_distance(targCom, moleculeList);
    const double expected = 2.0 * membraneObject.compartmentR - beforeMin;
    std::cerr << "  Minimum rotated radial distance after the call  = " << afterMin
              << " (expected " << expected << ")\n";

    EXPECT_NEAR(afterMin, expected, 1e-9)
        << "the rotated innermost point should be mirrored to 2*sphereR - targR";
    EXPECT_GE(afterMin, membraneObject.compartmentR)
        << "after reflection no rotated point may remain inside the compartment";

    // Some radial displacement must have been generated (non-trivial reflection).
    const double transMag = std::sqrt(targCom.trajTrans.x * targCom.trajTrans.x
        + targCom.trajTrans.y * targCom.trajTrans.y + targCom.trajTrans.z * targCom.trajTrans.z);
    std::cerr << "  |trajTrans| after the call = " << transMag << '\n';
    EXPECT_GT(transMag, 0.0) << "a reflection must produce a non-zero translation";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario runs independently so that a failure in
// one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(ReflectTrajComplexCompartment, GuardBlocksReflection) { test_rtcc_guard_blocks_reflection(); }
TEST(ReflectTrajComplexCompartment, ReflectsMoleculeOutOfCompartment) { test_rtcc_reflects_molecule_out_of_compartment(); }
TEST(ReflectTrajComplexCompartment, AccumulatesExistingTranslation) { test_rtcc_accumulates_existing_translation(); }
TEST(ReflectTrajComplexCompartment, InterfaceDrivesReflection) { test_rtcc_interface_drives_reflection(); }
TEST(ReflectTrajComplexCompartment, NoReflectionWhenAllPointsOutside) { test_rtcc_no_reflection_when_all_points_outside(); }
TEST(ReflectTrajComplexCompartment, OnSurfaceIgnoresRS3D) { test_rtcc_onsurface_ignores_RS3D(); }
TEST(ReflectTrajComplexCompartment, AccountsForRotation) { test_rtcc_accounts_for_rotation(); }