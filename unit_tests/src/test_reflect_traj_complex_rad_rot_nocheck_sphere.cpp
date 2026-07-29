/*! \file test_reflect_traj_complex_rad_rot_nocheck_sphere.cpp
 *
 * ### Unit test for src/boundary_conditions/reflect_traj_complex_rad_rot_nocheck_sphere.cpp
 *
 * Function under test:
 *
 *     void reflect_traj_complex_rad_rot_nocheck_sphere(const Parameters& params,
 *                                                      Complex& targCom,
 *                                                      std::vector<Molecule>& moleculeList,
 *                                                      const Membrane& membraneObject,
 *                                                      double RS3Dinput)
 *
 * Behaviour that is verified here (straight from the source):
 *
 *   1. If `targCom.OnSurface == true` the reflecting surface offset RS3D is
 *      forced to zero and the routine does *nothing at all* (a complex living
 *      on the sphere surface only moves in theta/phi, so it can never leave).
 *   2. Otherwise the effective radius is `sphereR = membraneObject.sphereR - RS3D`.
 *      If `|comCoord + trajTrans| + targCom.radius <= sphereR` the routine
 *      returns immediately (a cheap bounding-sphere early-out) and the
 *      trajectory is untouched.
 *   3. If that cheap test fails, every member molecule COM and every member
 *      interface is rotated by the Euler matrix built from `targCom.trajRot`
 *      and translated by `trajTrans`; the point furthest from the sphere centre
 *      is located.  A correction
 *          lamda  = -2 * (rtmp - sphereR) / rtmp
 *          dtrans = lamda * targcrds
 *      is then *added* to `trajTrans` so that the furthest point is mirrored
 *      back inside the sphere.  Note that the correction can be exactly zero
 *      if no real member point actually pokes outside the sphere (rtmp keeps
 *      its initial value of sphereR).
 *   4. Only `targCom.trajTrans` is modified — molecule coordinates and the
 *      complex COM are left alone.
 *
 * Every test prints what file/function is exercised, what the scenario is and
 * what the pass criterion is, so failures are self-describing.
 */

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (file-static, unique names so they cannot clash with the rest
// of the test suite).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal molecule: one interface, coincident with the COM.
 *
 * \param[in] com   center-of-mass coordinate of the molecule
 * \param[in] iface coordinate of its single interface
 * \return a Molecule that belongs to complex index 0
 */
Molecule rtcrrncs_make_molecule(const Coord& com, const Coord& iface)
{
    Molecule mol;
    mol.comCoord = com;

    Molecule::Iface oneIface;
    oneIface.coord = iface;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(oneIface);

    mol.myComIndex = 0; // single complex in all of these tests
    return mol;
}

/*! \brief Build a Complex whose memberList references moleculeList entries 0..n-1.
 *
 * \param[in] com      complex center of mass
 * \param[in] radius   bounding radius used by the cheap early-out test
 * \param[in] nMembers number of member molecules (indices 0..nMembers-1)
 */
Complex rtcrrncs_make_complex(const Coord& com, double radius, int nMembers)
{
    Complex targCom;
    targCom.comCoord = com;
    targCom.radius = radius;

    // Diffusion constants: non-zero so nothing is interpreted as membrane-bound
    // purely from the diffusion values (the routine keys off OnSurface anyway).
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Coord { 0.01, 0.01, 0.01 };

    // No proposed motion yet; individual tests set these.
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot = Coord { 0.0, 0.0, 0.0 };

    targCom.memberList.clear();
    for (int i = 0; i < nMembers; ++i)
        targCom.memberList.push_back(i);

    targCom.OnSurface = false;
    return targCom;
}

/*! \brief Recompute, exactly the way the routine does, the largest radial
 *         distance reached by any member COM/interface after the (possibly
 *         corrected) trajectory is applied.
 *
 * This lets the tests assert on the physically meaningful outcome ("nothing
 * pokes out of the sphere") instead of on internal bookkeeping.
 */
double rtcrrncs_max_member_radius(const Complex& targCom, const std::vector<Molecule>& moleculeList)
{
    std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);
    double maxR = 0.0;

    for (auto memMol : targCom.memberList) {
        // Member center of mass.
        Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
        Vector rotComVec { matrix_rotate(comVec, M) };
        Coord curr { targCom.comCoord + targCom.trajTrans + rotComVec };
        double currR = curr.get_magnitude();
        if (currR > maxR)
            maxR = currR;

        // Every interface of that member.
        for (const auto& iface : moleculeList[memMol].interfaceList) {
            Vector ifaceVec { iface.coord - targCom.comCoord };
            Vector rotIfaceVec { matrix_rotate(ifaceVec, M) };
            Coord currI { targCom.comCoord + targCom.trajTrans + rotIfaceVec };
            double currIR = currI.get_magnitude();
            if (currIR > maxR)
                maxR = currIR;
        }
    }
    return maxR;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: complex flagged OnSurface -> routine must be a complete no-op.
// -----------------------------------------------------------------------------
void test_rtcrrncs_on_surface_is_noop()
{
    std::cerr << "\n[TEST] test_rtcrrncs_on_surface_is_noop\n"
              << "  Source file:   reflect_traj_complex_rad_rot_nocheck_sphere.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_sphere\n"
              << "  Scenario:      targCom.OnSurface == true and the proposed step\n"
              << "                 would take the complex far outside the sphere.\n"
              << "  Pass criteria: trajTrans, trajRot, comCoord and the molecule\n"
              << "                 coordinates are all left untouched.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Complex sitting essentially on the shell, marked as surface-bound.
    Complex targCom = rtcrrncs_make_complex(Coord { 100.0, 0.0, 0.0 }, 1.0, 1);
    targCom.OnSurface = true;
    targCom.trajTrans.x = 50.0; // absurdly large: would leave the sphere
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;

    std::vector<Molecule> moleculeList {
        rtcrrncs_make_molecule(Coord { 100.0, 0.0, 0.0 }, Coord { 100.0, 0.0, 0.0 })
    };

    // A non-zero RS3D is supplied on purpose: for OnSurface complexes the
    // routine overrides it with zero, but either way nothing should happen.
    const double RS3Dinput = 10.0;

    std::cerr << "  Calling reflect_traj_complex_rad_rot_nocheck_sphere (RS3Dinput = "
              << RS3Dinput << ")...\n";
    reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, RS3Dinput);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 50.0) << "OnSurface complex: trajTrans.x must not change";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "OnSurface complex: trajTrans.y must not change";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0) << "OnSurface complex: trajTrans.z must not change";

    EXPECT_DOUBLE_EQ(targCom.comCoord.x, 100.0) << "COM of the complex must never be edited here";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 100.0) << "Member molecule COM must never be edited here";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, 100.0)
        << "Member interface coordinate must never be edited here";

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: interior complex, cheap bounding-sphere test succeeds -> early return.
// -----------------------------------------------------------------------------
void test_rtcrrncs_interior_early_return()
{
    std::cerr << "\n[TEST] test_rtcrrncs_interior_early_return\n"
              << "  Source file:   reflect_traj_complex_rad_rot_nocheck_sphere.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_sphere\n"
              << "  Scenario:      |comCoord + trajTrans| + radius <= sphereR, i.e. the\n"
              << "                 whole bounding sphere of the complex stays inside.\n"
              << "  Pass criteria: trajTrans returned bit-for-bit unchanged.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Complex near the origin, bounding radius 5, small step.
    Complex targCom = rtcrrncs_make_complex(Coord { 0.0, 0.0, 0.0 }, 5.0, 1);
    targCom.trajTrans.x = 1.0;
    targCom.trajTrans.y = -1.0;
    targCom.trajTrans.z = 0.5;

    std::vector<Molecule> moleculeList {
        rtcrrncs_make_molecule(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 })
    };

    const double origX = targCom.trajTrans.x;
    const double origY = targCom.trajTrans.y;
    const double origZ = targCom.trajTrans.z;

    std::cerr << "  Calling reflect_traj_complex_rad_rot_nocheck_sphere (RS3Dinput = 0)...\n";
    reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, 0.0);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, origX) << "Interior complex: trajTrans.x must be unchanged";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, origY) << "Interior complex: trajTrans.y must be unchanged";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, origZ) << "Interior complex: trajTrans.z must be unchanged";

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 3: complex would exit through the shell -> exact reflection arithmetic.
// -----------------------------------------------------------------------------
void test_rtcrrncs_reflects_point_complex()
{
    std::cerr << "\n[TEST] test_rtcrrncs_reflects_point_complex\n"
              << "  Source file:   reflect_traj_complex_rad_rot_nocheck_sphere.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_sphere\n"
              << "  Scenario:      single point-like molecule at r = 90 with a +20 x step,\n"
              << "                 sphere radius 100, so the step lands at r = 110.\n"
              << "  Pass criteria: lamda = -2*(110-100)/110 gives dtrans = (-20,0,0), so\n"
              << "                 trajTrans becomes (0,0,0) and the mirrored point sits at\n"
              << "                 r = 2*R - rtmp = 90 <= sphereR.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // radius = 0 so the cheap test is decided purely by the COM position.
    Complex targCom = rtcrrncs_make_complex(Coord { 90.0, 0.0, 0.0 }, 0.0, 1);
    targCom.trajTrans.x = 20.0;

    std::vector<Molecule> moleculeList {
        rtcrrncs_make_molecule(Coord { 90.0, 0.0, 0.0 }, Coord { 90.0, 0.0, 0.0 })
    };

    std::cerr << "  Calling reflect_traj_complex_rad_rot_nocheck_sphere (RS3Dinput = 0)...\n";
    reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, 0.0);

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";

    // Expected exact correction: 20 + (-2*(110-100)/110)*110 = 20 - 20 = 0.
    EXPECT_NEAR(targCom.trajTrans.x, 0.0, 1e-9)
        << "trajTrans.x should be reduced to 0 by the mirror correction";
    EXPECT_NEAR(targCom.trajTrans.y, 0.0, 1e-9) << "y component of the correction must stay 0";
    EXPECT_NEAR(targCom.trajTrans.z, 0.0, 1e-9) << "z component of the correction must stay 0";

    // Physical criterion: nothing sticks out of the sphere anymore.
    const double maxR = rtcrrncs_max_member_radius(targCom, moleculeList);
    std::cerr << "  Largest member radius after reflection = " << maxR
              << " (sphere radius " << membraneObject.sphereR << ")\n";
    EXPECT_LE(maxR, membraneObject.sphereR + 1e-6)
        << "After reflection every member point must lie inside the sphere";

    // The routine must not have moved anything but the trajectory.
    EXPECT_DOUBLE_EQ(targCom.comCoord.x, 90.0) << "Complex COM must not be modified";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 90.0) << "Molecule COM must not be modified";
}

// -----------------------------------------------------------------------------
// Test 4: cheap test fails (bounding radius is huge) but no real member point
//         is outside -> rtmp keeps its initial value sphereR, lamda == 0.
// -----------------------------------------------------------------------------
void test_rtcrrncs_conservative_bound_no_correction()
{
    std::cerr << "\n[TEST] test_rtcrrncs_conservative_bound_no_correction\n"
              << "  Source file:   reflect_traj_complex_rad_rot_nocheck_sphere.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_sphere\n"
              << "  Scenario:      targCom.radius is deliberately over-estimated (99.9) so\n"
              << "                 the cheap bounding test fails, but the only member point\n"
              << "                 sits at r = 1 and therefore never exceeds sphereR.\n"
              << "  Pass criteria: rtmp stays at sphereR -> lamda = 0 -> trajTrans unchanged.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    Complex targCom = rtcrrncs_make_complex(Coord { 0.0, 0.0, 0.0 }, 99.9, 1);
    targCom.trajTrans.x = 1.0; // |0 + 1| + 99.9 = 100.9 > 100 -> no early return

    std::vector<Molecule> moleculeList {
        rtcrrncs_make_molecule(Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 })
    };

    std::cerr << "  Calling reflect_traj_complex_rad_rot_nocheck_sphere (RS3Dinput = 0)...\n";
    reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, 0.0);

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";

    EXPECT_NEAR(targCom.trajTrans.x, 1.0, 1e-12)
        << "No member point is outside, so the correction must be exactly zero";
    EXPECT_NEAR(targCom.trajTrans.y, 0.0, 1e-12) << "y trajectory must stay zero";
    EXPECT_NEAR(targCom.trajTrans.z, 0.0, 1e-12) << "z trajectory must stay zero";
}

// -----------------------------------------------------------------------------
// Test 5: RS3Dinput shrinks the effective reflecting radius for a complex that
//         is *not* on the surface.
// -----------------------------------------------------------------------------
void test_rtcrrncs_rs3d_shrinks_radius()
{
    std::cerr << "\n[TEST] test_rtcrrncs_rs3d_shrinks_radius\n"
              << "  Source file:   reflect_traj_complex_rad_rot_nocheck_sphere.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_sphere\n"
              << "  Scenario:      identical geometry run twice: the final radius is 95 with\n"
              << "                 membrane sphereR = 100.\n"
              << "                   (a) RS3Dinput = 0  -> effective radius 100, 95 is inside\n"
              << "                       -> no correction expected;\n"
              << "                   (b) RS3Dinput = 10 -> effective radius 90, 95 is outside\n"
              << "                       -> correction expected, final r <= 90.\n"
              << "  Pass criteria: (a) trajTrans unchanged; (b) trajTrans reduced and the\n"
              << "                 furthest member point is within sphereR - RS3D.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // ---- (a) RS3Dinput = 0: nothing should happen ----------------------------
    {
        Complex targCom = rtcrrncs_make_complex(Coord { 85.0, 0.0, 0.0 }, 0.0, 1);
        targCom.trajTrans.x = 10.0; // lands at r = 95 < 100

        std::vector<Molecule> moleculeList {
            rtcrrncs_make_molecule(Coord { 85.0, 0.0, 0.0 }, Coord { 85.0, 0.0, 0.0 })
        };

        std::cerr << "  (a) Calling with RS3Dinput = 0...\n";
        reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, 0.0);

        std::cerr << "      trajTrans.x = " << targCom.trajTrans.x << " (expected 10)\n";
        EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 10.0)
            << "With RS3D = 0 the step of 10 keeps r = 95 < 100 and must be kept";
    }

    // ---- (b) RS3Dinput = 10: reflection required -----------------------------
    {
        const double RS3Dinput = 10.0;
        const double effectiveR = membraneObject.sphereR - RS3Dinput; // 90

        Complex targCom = rtcrrncs_make_complex(Coord { 85.0, 0.0, 0.0 }, 0.0, 1);
        targCom.trajTrans.x = 10.0; // would land at r = 95 > 90

        std::vector<Molecule> moleculeList {
            rtcrrncs_make_molecule(Coord { 85.0, 0.0, 0.0 }, Coord { 85.0, 0.0, 0.0 })
        };

        std::cerr << "  (b) Calling with RS3Dinput = " << RS3Dinput
                  << " (effective radius " << effectiveR << ")...\n";
        reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, RS3Dinput);

        std::cerr << "      trajTrans.x = " << targCom.trajTrans.x << " (expected ~0)\n";
        // lamda = -2*(95-90)/95 ; dtrans = lamda*95 = -10 -> trajTrans.x = 0
        EXPECT_NEAR(targCom.trajTrans.x, 0.0, 1e-9)
            << "With RS3D = 10 the outward step must be mirrored away";
        EXPECT_LT(targCom.trajTrans.x, 10.0)
            << "The corrected trajectory must be smaller than the proposed one";

        const double maxR = rtcrrncs_max_member_radius(targCom, moleculeList);
        std::cerr << "      Largest member radius = " << maxR
                  << " (must be <= " << effectiveR << ")\n";
        EXPECT_LE(maxR, effectiveR + 1e-6)
            << "After reflection the complex must be inside the RS3D-shrunk sphere";
    }
}

// -----------------------------------------------------------------------------
// Test 6: member offsets and the Euler rotation matrix are taken into account.
// -----------------------------------------------------------------------------
void test_rtcrrncs_uses_rotated_member_coords()
{
    std::cerr << "\n[TEST] test_rtcrrncs_uses_rotated_member_coords\n"
              << "  Source file:   reflect_traj_complex_rad_rot_nocheck_sphere.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_sphere\n"
              << "  Scenario:      the complex COM stays inside the sphere but one member\n"
              << "                 molecule, offset from the COM and rotated by a non-zero\n"
              << "                 trajRot, pokes outside after the proposed step.\n"
              << "  Pass criteria: trajTrans is corrected (not equal to the proposed value)\n"
              << "                 and, recomputing the rotated member positions, the\n"
              << "                 furthest point ends up inside the sphere.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // COM at r = 0 but the (single) member sits 95 nm away along +x.
    // The bounding radius is set to 95 so the cheap test fails once we step.
    Complex targCom = rtcrrncs_make_complex(Coord { 0.0, 0.0, 0.0 }, 95.0, 1);
    targCom.trajTrans.x = 10.0; // furthest member ends up near r = 105
    targCom.trajRot = Coord { 0.05, -0.03, 0.02 }; // small but non-trivial rotation

    std::vector<Molecule> moleculeList {
        rtcrrncs_make_molecule(Coord { 95.0, 0.0, 0.0 }, Coord { 95.0, 0.0, 0.0 })
    };

    // Sanity: before the call the rotated member really is outside the sphere.
    const double before = rtcrrncs_max_member_radius(targCom, moleculeList);
    std::cerr << "  Largest member radius BEFORE the call = " << before << "\n";
    EXPECT_GT(before, membraneObject.sphereR)
        << "Test setup problem: the member should start outside the sphere";

    std::cerr << "  Calling reflect_traj_complex_rad_rot_nocheck_sphere (RS3Dinput = 0)...\n";
    reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, 0.0);

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";

    // The correction must actually have fired.
    EXPECT_LT(targCom.trajTrans.x, 10.0)
        << "The x trajectory must be pulled back because a rotated member left the sphere";

    const double after = rtcrrncs_max_member_radius(targCom, moleculeList);
    std::cerr << "  Largest member radius AFTER  the call = " << after
              << " (sphere radius " << membraneObject.sphereR << ")\n";
    EXPECT_LE(after, membraneObject.sphereR + 1e-6)
        << "After reflection the rotated member must be back inside the sphere";

    // trajRot itself is only read, never written.
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, 0.05) << "trajRot.x must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, -0.03) << "trajRot.y must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, 0.02) << "trajRot.z must not be modified";
}

// -----------------------------------------------------------------------------
// Test 7: multi-member complex — the *furthest* point drives the correction.
// -----------------------------------------------------------------------------
void test_rtcrrncs_multi_member_uses_furthest_point()
{
    std::cerr << "\n[TEST] test_rtcrrncs_multi_member_uses_furthest_point\n"
              << "  Source file:   reflect_traj_complex_rad_rot_nocheck_sphere.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_sphere\n"
              << "  Scenario:      two member molecules at r = 80 and r = 92 (plus an\n"
              << "                 interface even further out at r = 94); the +x step of 10\n"
              << "                 pushes the furthest interface to r = 104.\n"
              << "  Pass criteria: the correction is computed from the FURTHEST point\n"
              << "                 (r = 104), i.e. dtrans.x = -2*(104-100)/104*104 = -8,\n"
              << "                 leaving trajTrans.x ~= 2 and the furthest point at r = 96.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // COM chosen at the origin, zero rotation so the arithmetic is exact.
    Complex targCom = rtcrrncs_make_complex(Coord { 0.0, 0.0, 0.0 }, 94.0, 2);
    targCom.trajTrans.x = 10.0;
    targCom.trajRot = Coord { 0.0, 0.0, 0.0 }; // identity rotation matrix

    std::vector<Molecule> moleculeList;
    // Member 0: modest distance, should NOT drive the correction.
    moleculeList.push_back(rtcrrncs_make_molecule(Coord { 80.0, 0.0, 0.0 }, Coord { 80.0, 0.0, 0.0 }));
    // Member 1: COM at 92, interface even further out at 94 -> the winner.
    moleculeList.push_back(rtcrrncs_make_molecule(Coord { 92.0, 0.0, 0.0 }, Coord { 94.0, 0.0, 0.0 }));

    std::cerr << "  Calling reflect_traj_complex_rad_rot_nocheck_sphere (RS3Dinput = 0)...\n";
    reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, 0.0);

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";

    // rtmp = 94 + 10 = 104 ; lamda = -2*4/104 ; dtrans.x = lamda*104 = -8.
    EXPECT_NEAR(targCom.trajTrans.x, 2.0, 1e-9)
        << "Correction must be derived from the furthest interface (r = 104), giving 10 - 8 = 2";

    const double maxR = rtcrrncs_max_member_radius(targCom, moleculeList);
    std::cerr << "  Largest member radius after reflection = " << maxR
              << " (expected 96, must be <= " << membraneObject.sphereR << ")\n";
    EXPECT_NEAR(maxR, 96.0, 1e-9) << "Mirrored furthest point should sit at 2*R - rtmp = 96";
    EXPECT_LE(maxR, membraneObject.sphereR + 1e-6)
        << "Every member point must be inside the sphere after the correction";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each scenario is its own TEST so a failure in one does
// not stop the others from running (all assertions are non-fatal EXPECT_*).
// -----------------------------------------------------------------------------
TEST(ReflectTrajComplexRadRotNocheckSphere, OnSurfaceIsNoop) { test_rtcrrncs_on_surface_is_noop(); }
TEST(ReflectTrajComplexRadRotNocheckSphere, InteriorEarlyReturn) { test_rtcrrncs_interior_early_return(); }
TEST(ReflectTrajComplexRadRotNocheckSphere, ReflectsPointComplex) { test_rtcrrncs_reflects_point_complex(); }
TEST(ReflectTrajComplexRadRotNocheckSphere, ConservativeBoundNoCorrection) { test_rtcrrncs_conservative_bound_no_correction(); }
TEST(ReflectTrajComplexRadRotNocheckSphere, RS3DShrinksRadius) { test_rtcrrncs_rs3d_shrinks_radius(); }
TEST(ReflectTrajComplexRadRotNocheckSphere, UsesRotatedMemberCoords) { test_rtcrrncs_uses_rotated_member_coords(); }
TEST(ReflectTrajComplexRadRotNocheckSphere, MultiMemberUsesFurthestPoint) { test_rtcrrncs_multi_member_uses_furthest_point(); }