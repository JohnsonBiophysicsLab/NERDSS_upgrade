/*! \file test_reflect_traj_tmp_crds_sphere.cpp
 *
 * ### Unit test for ../src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp
 *
 * Function under test:
 * \code
 *   void reflect_traj_tmp_crds_sphere(const Parameters& params,
 *                                    std::vector<Molecule>& moleculeList,
 *                                    Complex& targCom,
 *                                    std::array<double, 3>& traj,
 *                                    const Membrane& membraneObject,
 *                                    double radius,
 *                                    double RS3Dinput);
 * \endcode
 *
 * The routine evaluates whether a *tentatively* moved complex (association uses
 * the tmp* coordinate set) would poke outside a spherical boundary of radius
 *   sphereR = radius - RS3D
 * where RS3D == RS3Dinput unless the complex is (or will be) on the surface, in
 * which case RS3D == 0.
 *
 * Algorithm summary (used to derive the expected values in every test below):
 *   1. A cheap pre-check: if |tmpComCoord + traj| + targCom.radius <= sphereR the
 *      complex cannot possibly stick out and the function returns immediately,
 *      leaving `traj` untouched.
 *   2. Otherwise every member molecule COM and every member interface is
 *      displaced by `traj` (the rotation matrix is the identity here) and the
 *      largest radial excursion `dr = |curr| - sphereR` is recorded together
 *      with the offending coordinate `targcrds`.
 *   3. If any point is outside (dr > 0) the trajectory is *replaced* by
 *         traj = lamda * targcrds,  lamda = -2 * dr / |targcrds|
 *      i.e. traj = -2 * dr * (unit vector along targcrds).
 *
 * Note that the function never modifies the complex or molecule coordinates --
 * only the `traj` array that is accumulated by the caller.
 *
 * Every test prints what it is doing, which source file/function it exercises,
 * and the criteria used to decide pass/fail.
 */

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary_conditions/reflect_functions.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (prefixed rttcs_ = Reflect Traj Tmp Crds Sphere) so that they
// cannot collide with helpers in other translation units of the test suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal molecule whose tmp coordinates are explicitly set.
 *
 * \param[in] tmpCom  temporary center-of-mass coordinate (used by the routine)
 * \param[in] tmpIfaces temporary interface coordinates (absolute, not relative)
 * \param[in] comIndex index of the parent complex
 *
 * The regular (non-tmp) coordinates are set to the same values so the molecule
 * is self-consistent, but the routine under test only reads the tmp* fields.
 */
Molecule rttcs_make_molecule(const Coord& tmpCom, const std::vector<Coord>& tmpIfaces, int comIndex = 0)
{
    Molecule mol;
    mol.comCoord = tmpCom;
    mol.tmpComCoord = tmpCom;
    mol.myComIndex = comIndex;

    // interfaceList and tmpICoords must have identical sizes: the routine loops
    // over interfaceList.size() but indexes tmpICoords.
    mol.interfaceList.clear();
    mol.tmpICoords.clear();
    for (const auto& oneIface : tmpIfaces) {
        Molecule::Iface iface;
        iface.coord = oneIface;
        mol.interfaceList.push_back(iface);
        mol.tmpICoords.push_back(oneIface);
    }
    return mol;
}

/*! \brief Build a complex holding the given member molecule indices.
 *
 * \param[in] tmpCom      temporary center of mass of the complex
 * \param[in] comRadius   bounding radius used by the cheap pre-check
 * \param[in] memberList  indices into moleculeList
 */
Complex rttcs_make_complex(const Coord& tmpCom, double comRadius, const std::vector<int>& memberList)
{
    Complex targCom;
    targCom.comCoord = tmpCom;
    targCom.tmpComCoord = tmpCom;
    targCom.radius = comRadius;

    // Diffusion constants are not used by this routine, but give sane values.
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Coord { 0.01, 0.01, 0.01 };

    targCom.memberList = memberList;
    targCom.OnSurface = false;
    targCom.tmpOnSurface = false;
    return targCom;
}

/*! \brief Magnitude of a plain xyz triple (used for expected-value math). */
double rttcs_magnitude(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z); }

/*! \brief Pretty-print a trajectory triple to stderr. */
void rttcs_print_traj(const char* label, const std::array<double, 3>& traj)
{
    std::cerr << "    " << label << " = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: complex deep inside the sphere -> the cheap pre-check short-circuits
//         and `traj` must be returned bit-for-bit unchanged.
// -----------------------------------------------------------------------------
void test_rttcs_deep_inside_leaves_traj_unchanged()
{
    std::cerr << "\n[TEST] test_rttcs_deep_inside_leaves_traj_unchanged\n"
              << "  Source file:   src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_sphere()\n"
              << "  Scenario:      complex sits at the origin of a R=100 sphere with a\n"
              << "                 1 nm step; |com+traj| + comRadius is far below R.\n"
              << "  Pass criteria: traj is returned completely unmodified.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject; // membraneObject itself is unused by this routine
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // One molecule whose COM and single interface both sit at the origin.
    std::vector<Molecule> moleculeList { rttcs_make_molecule(Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 0.0, 0.0 } }) };
    Complex targCom = rttcs_make_complex(Coord { 0.0, 0.0, 0.0 }, /*comRadius*/ 2.0, { 0 });

    std::array<double, 3> traj { 1.0, -1.0, 0.5 };
    const std::array<double, 3> orig = traj;

    std::cerr << "  Calling reflect_traj_tmp_crds_sphere(radius=100, RS3Dinput=0)...\n";
    reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 0.0);
    rttcs_print_traj("traj after call", traj);

    // Nothing should have been touched.
    EXPECT_DOUBLE_EQ(traj[0], orig[0]) << "traj[0] must be unchanged for an interior complex";
    EXPECT_DOUBLE_EQ(traj[1], orig[1]) << "traj[1] must be unchanged for an interior complex";
    EXPECT_DOUBLE_EQ(traj[2], orig[2]) << "traj[2] must be unchanged for an interior complex";
}

// -----------------------------------------------------------------------------
// Test 2: the cheap pre-check fires (large bounding radius) but the detailed
//         per-molecule / per-interface scan finds nothing outside, so `traj`
//         must still be unchanged.
// -----------------------------------------------------------------------------
void test_rttcs_precheck_triggers_but_no_reflection()
{
    std::cerr << "\n[TEST] test_rttcs_precheck_triggers_but_no_reflection\n"
              << "  Source file:   src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_sphere()\n"
              << "  Scenario:      |com| (95) + comRadius (20) > R (100) so the coarse\n"
              << "                 test flags the complex, but the only real member point\n"
              << "                 is at r = 95 < 100 (dr stays 0).\n"
              << "  Pass criteria: traj is left unchanged (no false-positive reflection).\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Molecule COM and interface coincide with the complex COM at r = 95.
    const Coord com { 95.0, 0.0, 0.0 };
    std::vector<Molecule> moleculeList { rttcs_make_molecule(com, { com }) };
    // Deliberately inflated bounding radius so that the pre-check succeeds.
    Complex targCom = rttcs_make_complex(com, /*comRadius*/ 20.0, { 0 });

    std::array<double, 3> traj { 0.0, 0.0, 0.0 };

    std::cerr << "  Calling reflect_traj_tmp_crds_sphere(radius=100, RS3Dinput=0)...\n";
    reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 0.0);
    rttcs_print_traj("traj after call", traj);

    EXPECT_DOUBLE_EQ(traj[0], 0.0) << "traj[0] must stay 0: no member point is outside the sphere";
    EXPECT_DOUBLE_EQ(traj[1], 0.0) << "traj[1] must stay 0: no member point is outside the sphere";
    EXPECT_DOUBLE_EQ(traj[2], 0.0) << "traj[2] must stay 0: no member point is outside the sphere";
}

// -----------------------------------------------------------------------------
// Test 3: the molecule COM would leave the sphere -> traj is replaced by
//         -2 * dr * unit(targcrds), and the resulting position lands back
//         inside the sphere.
// -----------------------------------------------------------------------------
void test_rttcs_com_outside_is_reflected()
{
    std::cerr << "\n[TEST] test_rttcs_com_outside_is_reflected\n"
              << "  Source file:   src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_sphere()\n"
              << "  Scenario:      com = (95,0,0), traj = (20,0,0) -> candidate point at\n"
              << "                 r = 115 with R = 100 so dr = 15.\n"
              << "  Pass criteria: traj == -2*dr*unit(targcrds) == (-30,0,0) and the\n"
              << "                 reflected position |com+traj| = 65 <= R.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    const Coord com { 95.0, 0.0, 0.0 };
    std::vector<Molecule> moleculeList { rttcs_make_molecule(com, { com }) };
    Complex targCom = rttcs_make_complex(com, /*comRadius*/ 1.0, { 0 });

    std::array<double, 3> traj { 20.0, 0.0, 0.0 };

    // Analytic expectation: targcrds = (115,0,0), dr = 15, lamda = -2*15/115
    const double sphereR = 100.0;
    const double targMag = rttcs_magnitude(115.0, 0.0, 0.0);
    const double dr = targMag - sphereR;
    const double lamda = -2.0 * dr / targMag;
    const double expectedX = lamda * 115.0;

    std::cerr << "  Expected traj[0] = " << expectedX << " (analytic -2*dr = -30)\n";
    std::cerr << "  Calling reflect_traj_tmp_crds_sphere(radius=100, RS3Dinput=0)...\n";
    reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 0.0);
    rttcs_print_traj("traj after call", traj);

    EXPECT_NEAR(traj[0], expectedX, 1e-9) << "traj[0] must equal lamda*targcrds.x";
    EXPECT_NEAR(traj[0], -30.0, 1e-9) << "closed-form value for this geometry is -30";
    EXPECT_NEAR(traj[1], 0.0, 1e-12) << "no y displacement is expected for a pure +x excursion";
    EXPECT_NEAR(traj[2], 0.0, 1e-12) << "no z displacement is expected for a pure +x excursion";

    // The reflected end point must be back inside the sphere.
    const double finalR = rttcs_magnitude(com.x + traj[0], com.y + traj[1], com.z + traj[2]);
    std::cerr << "    reflected |com+traj| = " << finalR << " (sphere R = " << sphereR << ")\n";
    EXPECT_LE(finalR, sphereR + 1e-9) << "reflected position must lie inside the spherical boundary";
}

// -----------------------------------------------------------------------------
// Test 4: only an *interface* pokes out (the COM stays inside) -> the interface
//         excursion must drive the reflection.
// -----------------------------------------------------------------------------
void test_rttcs_interface_outside_is_reflected()
{
    std::cerr << "\n[TEST] test_rttcs_interface_outside_is_reflected\n"
              << "  Source file:   src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_sphere()\n"
              << "  Scenario:      com+traj at r = 95 (inside) but an interface offset by\n"
              << "                 +9 nm reaches r = 104 with R = 100 so dr = 4.\n"
              << "  Pass criteria: traj == (-8,0,0) == -2*dr along +x.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Complex COM at r = 90, interface 9 nm further out along +x.
    const Coord com { 90.0, 0.0, 0.0 };
    const Coord iface { 99.0, 0.0, 0.0 };
    std::vector<Molecule> moleculeList { rttcs_make_molecule(com, { iface }) };
    // comRadius = 9 makes the coarse pre-check succeed: 95 + 9 > 100.
    Complex targCom = rttcs_make_complex(com, /*comRadius*/ 9.0, { 0 });

    std::array<double, 3> traj { 5.0, 0.0, 0.0 };

    // Analytic expectation: interface candidate = (104,0,0), dr = 4 -> traj = -8 x̂
    std::cerr << "  Expected traj = (-8, 0, 0)\n";
    std::cerr << "  Calling reflect_traj_tmp_crds_sphere(radius=100, RS3Dinput=0)...\n";
    reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 0.0);
    rttcs_print_traj("traj after call", traj);

    EXPECT_NEAR(traj[0], -8.0, 1e-9) << "interface excursion of 4 nm must yield traj[0] = -8";
    EXPECT_NEAR(traj[1], 0.0, 1e-12) << "no y component expected";
    EXPECT_NEAR(traj[2], 0.0, 1e-12) << "no z component expected";
}

// -----------------------------------------------------------------------------
// Test 5: several members -> the *largest* radial violation must be the one
//         that determines the reflection vector.
// -----------------------------------------------------------------------------
void test_rttcs_largest_violation_wins()
{
    std::cerr << "\n[TEST] test_rttcs_largest_violation_wins\n"
              << "  Source file:   src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_sphere()\n"
              << "  Scenario:      three members: r = 80 (inside), r = 110 (dr = 10) and\n"
              << "                 r = 130 (dr = 30) with R = 100.\n"
              << "  Pass criteria: traj follows the r = 130 point along -z, i.e. (0,0,-60),\n"
              << "                 proving the maximum dr is selected.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Members: inside, moderately outside (+y), badly outside (+z).
    std::vector<Molecule> moleculeList {
        rttcs_make_molecule(Coord { 80.0, 0.0, 0.0 }, { Coord { 80.0, 0.0, 0.0 } }),
        rttcs_make_molecule(Coord { 0.0, 110.0, 0.0 }, { Coord { 0.0, 110.0, 0.0 } }),
        rttcs_make_molecule(Coord { 0.0, 0.0, 130.0 }, { Coord { 0.0, 0.0, 130.0 } }),
    };
    // Complex COM at the origin with a large bounding radius so the pre-check fires.
    Complex targCom = rttcs_make_complex(Coord { 0.0, 0.0, 0.0 }, /*comRadius*/ 140.0, { 0, 1, 2 });

    std::array<double, 3> traj { 0.0, 0.0, 0.0 };

    // Analytic expectation: worst point (0,0,130), dr = 30 -> traj = -60 ẑ
    std::cerr << "  Expected traj = (0, 0, -60)\n";
    std::cerr << "  Calling reflect_traj_tmp_crds_sphere(radius=100, RS3Dinput=0)...\n";
    reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 0.0);
    rttcs_print_traj("traj after call", traj);

    EXPECT_NEAR(traj[0], 0.0, 1e-12) << "worst violation is purely along +z, so x must be 0";
    EXPECT_NEAR(traj[1], 0.0, 1e-12) << "the smaller +y violation must be superseded";
    EXPECT_NEAR(traj[2], -60.0, 1e-9) << "largest dr (30) must give traj[2] = -60";
}

// -----------------------------------------------------------------------------
// Test 6: RS3Dinput shrinks the effective boundary (sphereR = radius - RS3D)
//         only when the complex is NOT on the surface.
// -----------------------------------------------------------------------------
void test_rttcs_rs3d_shrinks_effective_radius()
{
    std::cerr << "\n[TEST] test_rttcs_rs3d_shrinks_effective_radius\n"
              << "  Source file:   src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_sphere()\n"
              << "  Scenario A:    OnSurface = false, radius = 100, RS3Dinput = 10 so the\n"
              << "                 effective boundary is 90; a point at r = 95 is outside.\n"
              << "  Scenario B:    OnSurface = true  -> RS3D forced to 0, boundary stays\n"
              << "                 100 and the same r = 95 point is inside.\n"
              << "  Pass criteria: A reflects to traj = (-10,0,0); B leaves traj at zero.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    const Coord com { 95.0, 0.0, 0.0 };

    // ---- Scenario A: 3D complex, RS3D applies -------------------------------
    {
        std::vector<Molecule> moleculeList { rttcs_make_molecule(com, { com }) };
        Complex targCom = rttcs_make_complex(com, /*comRadius*/ 1.0, { 0 });
        targCom.OnSurface = false;
        targCom.tmpOnSurface = false;

        std::array<double, 3> traj { 0.0, 0.0, 0.0 };

        std::cerr << "  [A] Calling with radius=100, RS3Dinput=10, OnSurface=false...\n";
        reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 10.0);
        rttcs_print_traj("[A] traj after call", traj);

        // sphereR = 90, targcrds = (95,0,0), dr = 5, traj = -2*5 = -10 along x.
        EXPECT_NEAR(traj[0], -10.0, 1e-9) << "RS3D=10 must shrink the boundary to 90 and reflect by -10";
        EXPECT_NEAR(traj[1], 0.0, 1e-12) << "no y component expected";
        EXPECT_NEAR(traj[2], 0.0, 1e-12) << "no z component expected";
    }

    // ---- Scenario B: complex on the surface, RS3D ignored -------------------
    {
        std::vector<Molecule> moleculeList { rttcs_make_molecule(com, { com }) };
        Complex targCom = rttcs_make_complex(com, /*comRadius*/ 1.0, { 0 });
        targCom.OnSurface = true; // forces RS3D = 0 inside the routine

        std::array<double, 3> traj { 0.0, 0.0, 0.0 };

        std::cerr << "  [B] Calling with radius=100, RS3Dinput=10, OnSurface=true...\n";
        reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 10.0);
        rttcs_print_traj("[B] traj after call", traj);

        // sphereR stays 100 and 95 + 1 < 100, so the pre-check exits immediately.
        EXPECT_DOUBLE_EQ(traj[0], 0.0) << "surface complex ignores RS3D: no reflection expected";
        EXPECT_DOUBLE_EQ(traj[1], 0.0) << "surface complex ignores RS3D: no reflection expected";
        EXPECT_DOUBLE_EQ(traj[2], 0.0) << "surface complex ignores RS3D: no reflection expected";
    }

    // ---- Scenario C: tmpOnSurface alone must also disable RS3D --------------
    {
        std::vector<Molecule> moleculeList { rttcs_make_molecule(com, { com }) };
        Complex targCom = rttcs_make_complex(com, /*comRadius*/ 1.0, { 0 });
        targCom.OnSurface = false;
        targCom.tmpOnSurface = true; // the "will be on surface" flag

        std::array<double, 3> traj { 0.0, 0.0, 0.0 };

        std::cerr << "  [C] Calling with radius=100, RS3Dinput=10, tmpOnSurface=true...\n";
        reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 10.0);
        rttcs_print_traj("[C] traj after call", traj);

        EXPECT_DOUBLE_EQ(traj[0], 0.0) << "tmpOnSurface must also force RS3D to 0";
        EXPECT_DOUBLE_EQ(traj[1], 0.0) << "tmpOnSurface must also force RS3D to 0";
        EXPECT_DOUBLE_EQ(traj[2], 0.0) << "tmpOnSurface must also force RS3D to 0";
    }
}

// -----------------------------------------------------------------------------
// Test 7: the routine must be side-effect free with respect to coordinates --
//         only `traj` may change.
// -----------------------------------------------------------------------------
void test_rttcs_coordinates_are_not_modified()
{
    std::cerr << "\n[TEST] test_rttcs_coordinates_are_not_modified\n"
              << "  Source file:   src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_sphere()\n"
              << "  Scenario:      a reflecting configuration (com would exit the sphere).\n"
              << "  Pass criteria: targCom.tmpComCoord, molecule tmpComCoord and tmpICoords\n"
              << "                 are byte-identical before and after the call.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    const Coord com { 96.0, 0.0, 0.0 };
    const Coord iface { 98.0, 1.0, 0.0 };
    std::vector<Molecule> moleculeList { rttcs_make_molecule(com, { iface }) };
    Complex targCom = rttcs_make_complex(com, /*comRadius*/ 5.0, { 0 });

    // Snapshot the coordinates we expect to be preserved.
    const Coord comBefore = targCom.tmpComCoord;
    const Coord molComBefore = moleculeList[0].tmpComCoord;
    const Coord ifaceBefore = moleculeList[0].tmpICoords[0];

    std::array<double, 3> traj { 10.0, 0.0, 0.0 }; // pushes the COM to r = 106

    std::cerr << "  Calling reflect_traj_tmp_crds_sphere(radius=100, RS3Dinput=0)...\n";
    reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 0.0);
    rttcs_print_traj("traj after call", traj);

    // The reflection itself should have happened (sanity check on the setup).
    EXPECT_LT(traj[0], 0.0) << "setup sanity: the +x excursion should produce a negative traj[0]";

    // Coordinates must be untouched.
    EXPECT_TRUE(targCom.tmpComCoord == comBefore) << "complex tmpComCoord must not be modified";
    EXPECT_TRUE(moleculeList[0].tmpComCoord == molComBefore) << "molecule tmpComCoord must not be modified";
    EXPECT_TRUE(moleculeList[0].tmpICoords[0] == ifaceBefore) << "molecule tmpICoords must not be modified";

    std::cerr << "    complex tmpComCoord still = " << targCom.tmpComCoord << "\n"
              << "    molecule tmpComCoord still = " << moleculeList[0].tmpComCoord << "\n"
              << "    molecule tmpICoords[0] still = " << moleculeList[0].tmpICoords[0] << "\n";
}

// -----------------------------------------------------------------------------
// Test 8: a molecule with no interfaces (interfaceList empty) must be handled
//         gracefully -- the COM check alone drives the decision.
// -----------------------------------------------------------------------------
void test_rttcs_handles_molecule_without_interfaces()
{
    std::cerr << "\n[TEST] test_rttcs_handles_molecule_without_interfaces\n"
              << "  Source file:   src/boundary_conditions/reflect_traj_tmp_crds_sphere.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_sphere()\n"
              << "  Scenario:      single member molecule with zero interfaces placed so\n"
              << "                 that its COM exits the sphere (r = 120, R = 100).\n"
              << "  Pass criteria: no crash, and traj = -2*dr along +x = (-40,0,0).\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Empty interface list -> the interface loop body never executes.
    std::vector<Molecule> moleculeList { rttcs_make_molecule(Coord { 100.0, 0.0, 0.0 }, {}) };
    Complex targCom = rttcs_make_complex(Coord { 100.0, 0.0, 0.0 }, /*comRadius*/ 1.0, { 0 });

    std::array<double, 3> traj { 20.0, 0.0, 0.0 }; // candidate COM at r = 120, dr = 20

    std::cerr << "  Expected traj = (-40, 0, 0)\n";
    std::cerr << "  Calling reflect_traj_tmp_crds_sphere(radius=100, RS3Dinput=0)...\n";
    reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, 100.0, 0.0);
    rttcs_print_traj("traj after call", traj);

    EXPECT_NEAR(traj[0], -40.0, 1e-9) << "dr = 20 must give traj[0] = -40 even without interfaces";
    EXPECT_NEAR(traj[1], 0.0, 1e-12) << "no y component expected";
    EXPECT_NEAR(traj[2], 0.0, 1e-12) << "no z component expected";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each named helper runs inside its own TEST so that a
// failure in one scenario does not stop the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(ReflectTrajTmpCrdsSphere, DeepInsideLeavesTrajUnchanged) { test_rttcs_deep_inside_leaves_traj_unchanged(); }
TEST(ReflectTrajTmpCrdsSphere, PrecheckTriggersButNoReflection) { test_rttcs_precheck_triggers_but_no_reflection(); }
TEST(ReflectTrajTmpCrdsSphere, ComOutsideIsReflected) { test_rttcs_com_outside_is_reflected(); }
TEST(ReflectTrajTmpCrdsSphere, InterfaceOutsideIsReflected) { test_rttcs_interface_outside_is_reflected(); }
TEST(ReflectTrajTmpCrdsSphere, LargestViolationWins) { test_rttcs_largest_violation_wins(); }
TEST(ReflectTrajTmpCrdsSphere, RS3DShrinksEffectiveRadius) { test_rttcs_rs3d_shrinks_effective_radius(); }
TEST(ReflectTrajTmpCrdsSphere, CoordinatesAreNotModified) { test_rttcs_coordinates_are_not_modified(); }
TEST(ReflectTrajTmpCrdsSphere, HandlesMoleculeWithoutInterfaces) { test_rttcs_handles_molecule_without_interfaces(); }