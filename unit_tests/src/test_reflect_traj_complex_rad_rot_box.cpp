/*! \file test_reflect_traj_complex_rad_rot_box.cpp
 *
 * ### Unit test for src/boundary_conditions/reflect_traj_complex_rad_rot_box.cpp
 *
 * Function under test:
 *
 *     void reflect_traj_complex_rad_rot_box(const Parameters& params,
 *                                          std::vector<Molecule>& moleculeList,
 *                                          Complex& targCom,
 *                                          const Membrane& membraneObject,
 *                                          double RS3Dinput)
 *
 * What the function does (and therefore what we verify):
 *
 *   * It builds an Euler rotation matrix from `targCom.trajRot` and uses it to
 *     find where every member molecule COM and every member interface *would*
 *     end up after the proposed translation + rotation.
 *   * The six walls of the rectangular water box are at +/- waterBox.{x,y,z}/2,
 *     except the -Z wall which is raised by the implicit-lipid reflecting
 *     surface RS3D (RS3D is forced to 0 when the complex is on the surface).
 *   * For each axis, if anything pokes through a wall by a distance d, the
 *     translation along that axis is corrected by -2*d (a specular reflection).
 *   * If the correction could push the complex out of the *opposite* wall (or
 *     the complex sticks out of both walls at once), the routine delegates to
 *     reflect_traj_check_span_box() which resamples the move.
 *
 * Every test prints what is being exercised and what the pass criterion is.
 */

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Small helpers used by every test in this file. All names are prefixed with
// "rtcrrb_" (Reflect Traj Complex Rad Rot Box) to avoid collisions with the
// other translation units linked into the same gtest binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal molecule whose single interface sits on its COM.
 *
 * \param[in] com          Center-of-mass coordinate for the molecule.
 * \param[in] parentComIdx Index of the parent complex (always 0 here).
 */
Molecule rtcrrb_make_molecule(const Coord& com, int parentComIdx = 0)
{
    Molecule mol;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = com; // interface coincident with the COM keeps geometry simple
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    mol.myComIndex = parentComIdx;
    return mol;
}

/*! \brief Build a complex with a given COM, bounding radius and member list.
 *
 * Diffusion constants are set small-but-nonzero by default; individual tests
 * override them when they want deterministic behaviour out of the resampling
 * code path.
 */
Complex rtcrrb_make_complex(const Coord& com, double radius, const std::vector<int>& members)
{
    Complex targCom;
    targCom.comCoord = com;
    targCom.radius = radius;

    targCom.D.x = 1.0;
    targCom.D.y = 1.0;
    targCom.D.z = 1.0;
    targCom.Dr.x = 0.01;
    targCom.Dr.y = 0.01;
    targCom.Dr.z = 0.01;

    // No proposed motion until the test sets it.
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot.x = 0.0;
    targCom.trajRot.y = 0.0;
    targCom.trajRot.z = 0.0;

    targCom.memberList = members;
    targCom.OnSurface = false;
    targCom.index = 0;

    return targCom;
}

/*! \brief Build a cubic (or rectangular) water box membrane object. */
Membrane rtcrrb_make_box(double x, double y, double z)
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.waterBox.x = x;
    membraneObject.waterBox.y = y;
    membraneObject.waterBox.z = z;
    membraneObject.waterBox.volume = x * y * z;
    return membraneObject;
}

/*! \brief Default parameters (only the timestep matters to the routine). */
Parameters rtcrrb_make_params()
{
    Parameters params;
    params.timeStep = 1.0;
    return params;
}

/*! \brief Reproduce the routine's own geometry: final position of a point.
 *
 * The function evaluates `com + trajTrans + M * (point - com)` for every
 * member COM/interface, so we use exactly the same expression to validate the
 * post-condition "everything is inside the box".
 */
Coord rtcrrb_final_position(const Complex& targCom, const Coord& point)
{
    std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);
    Vector rel { point - targCom.comCoord };

    Coord out;
    out.x = targCom.comCoord.x + targCom.trajTrans.x
        + (M[0] * rel.x + M[1] * rel.y + M[2] * rel.z);
    out.y = targCom.comCoord.y + targCom.trajTrans.y
        + (M[3] * rel.x + M[4] * rel.y + M[5] * rel.z);
    out.z = targCom.comCoord.z + targCom.trajTrans.z
        + (M[6] * rel.x + M[7] * rel.y + M[8] * rel.z);
    return out;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a complex safely inside the box must not be touched at all.
// -----------------------------------------------------------------------------
void test_rtcrrb_interior_no_change()
{
    std::cerr << "\n[TEST] test_rtcrrb_interior_no_change\n"
              << "  Source file:   reflect_traj_complex_rad_rot_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_box\n"
              << "  Scenario:      100x100x100 box, complex at the origin with a\n"
              << "                 1 nm step; nothing can reach a wall.\n"
              << "  Pass criteria: trajTrans and trajRot are bit-for-bit unchanged.\n";

    Parameters params = rtcrrb_make_params();
    Membrane membraneObject = rtcrrb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrb_make_molecule(Coord { 0.0, 0.0, 0.0 }) };
    Complex targCom = rtcrrb_make_complex(Coord { 0.0, 0.0, 0.0 }, 2.0, { 0 });

    targCom.trajTrans.x = 1.0;
    targCom.trajTrans.y = -1.0;
    targCom.trajTrans.z = 0.5;

    reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Resulting trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 1.0) << "interior complex: x step must not change";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, -1.0) << "interior complex: y step must not change";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.5) << "interior complex: z step must not change";

    // Rotation is never modified by this routine unless resampling happens.
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, 0.0) << "rotation must be untouched";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, 0.0) << "rotation must be untouched";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, 0.0) << "rotation must be untouched";
}

// -----------------------------------------------------------------------------
// Test 2: exceeding the +X wall by 15 nm must subtract 2*15 = 30 from trajTrans.x
//         and must leave Y and Z alone.
// -----------------------------------------------------------------------------
void test_rtcrrb_positive_x_reflection()
{
    std::cerr << "\n[TEST] test_rtcrrb_positive_x_reflection\n"
              << "  Source file:   reflect_traj_complex_rad_rot_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_box (+X wall)\n"
              << "  Scenario:      molecule at x=45 in a 100 nm box takes a +20 nm\n"
              << "                 step, i.e. it overshoots the wall (x=50) by 15 nm.\n"
              << "  Pass criteria: trajTrans.x == 20 - 2*15 == -10, y/z unchanged,\n"
              << "                 and the final position sits inside the wall.\n";

    Parameters params = rtcrrb_make_params();
    Membrane membraneObject = rtcrrb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrb_make_molecule(Coord { 45.0, 0.0, 0.0 }) };
    Complex targCom = rtcrrb_make_complex(Coord { 45.0, 0.0, 0.0 }, 2.0, { 0 });

    targCom.trajTrans.x = 20.0;
    targCom.trajTrans.y = 3.0;
    targCom.trajTrans.z = -4.0;

    reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Resulting trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n";

    EXPECT_NEAR(targCom.trajTrans.x, -10.0, 1e-12)
        << "specular reflection off the +X wall should give -10 nm";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 3.0) << "y step must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, -4.0) << "z step must not be modified";

    Coord fin = rtcrrb_final_position(targCom, moleculeList[0].comCoord);
    std::cerr << "  Final molecule x = " << fin.x << " (wall at 50)\n";
    EXPECT_LE(fin.x, 50.0 + 1e-6) << "reflected molecule must end up inside the +X wall";
}

// -----------------------------------------------------------------------------
// Test 3: exceeding the -Y wall must add 2*|d| to trajTrans.y.
// -----------------------------------------------------------------------------
void test_rtcrrb_negative_y_reflection()
{
    std::cerr << "\n[TEST] test_rtcrrb_negative_y_reflection\n"
              << "  Source file:   reflect_traj_complex_rad_rot_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_box (-Y wall)\n"
              << "  Scenario:      molecule at y=-46 takes a -8 nm step, overshooting\n"
              << "                 the -Y wall (y=-50) by 4 nm.\n"
              << "  Pass criteria: trajTrans.y == -8 - 2*(-4) == 0, x/z unchanged.\n";

    Parameters params = rtcrrb_make_params();
    Membrane membraneObject = rtcrrb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrb_make_molecule(Coord { 0.0, -46.0, 0.0 }) };
    Complex targCom = rtcrrb_make_complex(Coord { 0.0, -46.0, 0.0 }, 2.0, { 0 });

    targCom.trajTrans.x = 1.0;
    targCom.trajTrans.y = -8.0;
    targCom.trajTrans.z = 2.0;

    reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Resulting trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n";

    EXPECT_NEAR(targCom.trajTrans.y, 0.0, 1e-12)
        << "specular reflection off the -Y wall should cancel the -8 nm step";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 1.0) << "x step must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 2.0) << "z step must not be modified";

    Coord fin = rtcrrb_final_position(targCom, moleculeList[0].comCoord);
    std::cerr << "  Final molecule y = " << fin.y << " (wall at -50)\n";
    EXPECT_GE(fin.y, -50.0 - 1e-6) << "reflected molecule must end up inside the -Y wall";
}

// -----------------------------------------------------------------------------
// Test 4: the -Z wall is raised by RS3D for a complex that is *not* on the
//         membrane surface.
// -----------------------------------------------------------------------------
void test_rtcrrb_z_wall_uses_rs3d_when_not_on_surface()
{
    std::cerr << "\n[TEST] test_rtcrrb_z_wall_uses_rs3d_when_not_on_surface\n"
              << "  Source file:   reflect_traj_complex_rad_rot_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_box (-Z wall, RS3D)\n"
              << "  Scenario:      RS3D = 5 raises the -Z wall from -50 to -45;\n"
              << "                 molecule at z=-40 takes a -15 nm step (z=-55),\n"
              << "                 i.e. 10 nm past the raised wall.\n"
              << "  Pass criteria: trajTrans.z == -15 - 2*(-10) == +5.\n";

    Parameters params = rtcrrb_make_params();
    Membrane membraneObject = rtcrrb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrb_make_molecule(Coord { 0.0, 0.0, -40.0 }) };
    Complex targCom = rtcrrb_make_complex(Coord { 0.0, 0.0, -40.0 }, 2.0, { 0 });
    targCom.OnSurface = false; // -> RS3D input is honoured

    targCom.trajTrans.z = -15.0;

    const double RS3Dinput = 5.0;
    reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, RS3Dinput);

    std::cerr << "  Resulting trajTrans.z = " << targCom.trajTrans.z << "\n";

    EXPECT_NEAR(targCom.trajTrans.z, 5.0, 1e-12)
        << "reflection should use the RS3D-shifted -Z wall at -45";

    Coord fin = rtcrrb_final_position(targCom, moleculeList[0].comCoord);
    std::cerr << "  Final molecule z = " << fin.z << " (effective wall at -45)\n";
    EXPECT_GE(fin.z, -45.0 - 1e-6) << "molecule must not be below the reflecting surface";
}

// -----------------------------------------------------------------------------
// Test 5: for a complex bound to the surface RS3D is forced to zero, so the
//         same geometry reflects off the true box wall instead.
// -----------------------------------------------------------------------------
void test_rtcrrb_z_wall_ignores_rs3d_when_on_surface()
{
    std::cerr << "\n[TEST] test_rtcrrb_z_wall_ignores_rs3d_when_on_surface\n"
              << "  Source file:   reflect_traj_complex_rad_rot_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_box (-Z wall, OnSurface)\n"
              << "  Scenario:      identical to the previous test but the complex is\n"
              << "                 flagged OnSurface, so RS3D is internally set to 0\n"
              << "                 and the wall stays at -50 (overshoot 5 nm).\n"
              << "  Pass criteria: trajTrans.z == -15 - 2*(-5) == -5, i.e. different\n"
              << "                 from the not-on-surface case (+5).\n";

    Parameters params = rtcrrb_make_params();
    Membrane membraneObject = rtcrrb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrb_make_molecule(Coord { 0.0, 0.0, -40.0 }) };
    Complex targCom = rtcrrb_make_complex(Coord { 0.0, 0.0, -40.0 }, 2.0, { 0 });
    targCom.OnSurface = true; // -> RS3D forced to 0 inside the routine

    targCom.trajTrans.z = -15.0;

    const double RS3Dinput = 5.0; // deliberately non-zero; must be ignored
    reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, RS3Dinput);

    std::cerr << "  Resulting trajTrans.z = " << targCom.trajTrans.z << "\n";

    EXPECT_NEAR(targCom.trajTrans.z, -5.0, 1e-12)
        << "an on-surface complex must reflect off the real -Z wall at -50";

    Coord fin = rtcrrb_final_position(targCom, moleculeList[0].comCoord);
    std::cerr << "  Final molecule z = " << fin.z << " (wall at -50)\n";
    EXPECT_GE(fin.z, -50.0 - 1e-6) << "molecule must stay inside the box bottom";
}

// -----------------------------------------------------------------------------
// Test 6: the routine measures member molecules/interfaces, not just the
//         complex COM. Here the COM stays inside but a member sticks out.
// -----------------------------------------------------------------------------
void test_rtcrrb_member_offset_drives_reflection()
{
    std::cerr << "\n[TEST] test_rtcrrb_member_offset_drives_reflection\n"
              << "  Source file:   reflect_traj_complex_rad_rot_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_box (member scan)\n"
              << "  Scenario:      complex COM at x=40 with a member 8 nm further out\n"
              << "                 (x=48) takes a +5 nm step. The COM lands at 45\n"
              << "                 (inside) but the member lands at 53 (3 nm outside).\n"
              << "  Pass criteria: trajTrans.x == 5 - 2*3 == -1 (driven by the member).\n";

    Parameters params = rtcrrb_make_params();
    Membrane membraneObject = rtcrrb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrb_make_molecule(Coord { 48.0, 0.0, 0.0 }) };
    Complex targCom = rtcrrb_make_complex(Coord { 40.0, 0.0, 0.0 }, 8.0, { 0 });

    targCom.trajTrans.x = 5.0;

    reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Resulting trajTrans.x = " << targCom.trajTrans.x << "\n";

    EXPECT_NEAR(targCom.trajTrans.x, -1.0, 1e-12)
        << "the outermost member (not the COM) sets the reflection distance";

    Coord fin = rtcrrb_final_position(targCom, moleculeList[0].comCoord);
    std::cerr << "  Final member x = " << fin.x << " (wall at 50)\n";
    EXPECT_LE(fin.x, 50.0 + 1e-6) << "the offending member must be brought back inside";
}

// -----------------------------------------------------------------------------
// Test 7: the proposed rotation is taken into account. A 180 deg rotation about
//         z folds the protruding member back inside, so no reflection is needed.
// -----------------------------------------------------------------------------
void test_rtcrrb_rotation_is_applied()
{
    std::cerr << "\n[TEST] test_rtcrrb_rotation_is_applied\n"
              << "  Source file:   reflect_traj_complex_rad_rot_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_box (rotation matrix)\n"
              << "  Scenario:      complex COM at x=45 with a member offset +8 nm\n"
              << "                 (x=53, nominally outside the wall) but with a\n"
              << "                 trajRot of pi about z, which maps the offset to -8.\n"
              << "  Pass criteria: no reflection at all -> trajTrans stays (0,0,0),\n"
              << "                 proving the rotation matrix is used in the test.\n";

    Parameters params = rtcrrb_make_params();
    Membrane membraneObject = rtcrrb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrb_make_molecule(Coord { 53.0, 0.0, 0.0 }) };
    Complex targCom = rtcrrb_make_complex(Coord { 45.0, 0.0, 0.0 }, 8.0, { 0 });

    // No translation, but a half turn about the z axis.
    targCom.trajRot.z = M_PI;

    reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Resulting trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n";

    EXPECT_NEAR(targCom.trajTrans.x, 0.0, 1e-12)
        << "rotated member lands at x=37, inside the box, so no x reflection";
    EXPECT_NEAR(targCom.trajTrans.y, 0.0, 1e-12) << "no y reflection expected";
    EXPECT_NEAR(targCom.trajTrans.z, 0.0, 1e-12) << "no z reflection expected";

    Coord fin = rtcrrb_final_position(targCom, moleculeList[0].comCoord);
    std::cerr << "  Final (rotated) member x = " << fin.x << " (wall at 50)\n";
    EXPECT_LE(fin.x, 50.0 + 1e-6) << "the rotated member must be inside the +X wall";
}

// -----------------------------------------------------------------------------
// Test 8: the "recheck" branch. When the reflection would push the complex out
//         of the opposite wall the routine delegates to
//         reflect_traj_check_span_box(), which resamples the move.
//
//         To keep this deterministic the diffusion constants are set to zero:
//         the resampled displacement/rotation is then zero, which leaves the
//         (28 nm wide) complex comfortably inside the 30 nm box.
// -----------------------------------------------------------------------------
void test_rtcrrb_recheck_span_branch()
{
    std::cerr << "\n[TEST] test_rtcrrb_recheck_span_branch\n"
              << "  Source file:   reflect_traj_complex_rad_rot_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_box (recheck branch)\n"
              << "  Scenario:      narrow 30 nm box holding a 28 nm wide two-molecule\n"
              << "                 complex; a +4 nm step pokes out of the +X wall and\n"
              << "                 the 2*d correction would poke out of the -X wall,\n"
              << "                 which triggers reflect_traj_check_span_box().\n"
              << "                 D and Dr are zero so the resampled move is zero.\n"
              << "  Pass criteria: the call returns and every member ends up inside\n"
              << "                 both X walls (-15 <= x <= 15).\n";

    // reflect_traj_check_span_box() resamples with the GSL RNG; make sure the
    // generator exists (gtest_main only declares/defines the pointer as null).
    if (r == nullptr) {
        std::cerr << "  (GSL RNG was not initialized yet -> seeding it here)\n";
        const gsl_rng_type *T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
    std::cerr << "  GSL RNG initialized.\n";
    std::flush(std::cerr);
    Parameters params = rtcrrb_make_params();
    Membrane membraneObject = rtcrrb_make_box(30.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList {
        rtcrrb_make_molecule(Coord { -14.0, 0.0, 0.0 }),
        rtcrrb_make_molecule(Coord { 14.0, 0.0, 0.0 })
    };
    Complex targCom = rtcrrb_make_complex(Coord { 0.0, 0.0, 0.0 }, 14.0, { 0, 1 });

    // Zero diffusion -> deterministic (zero) resampling inside the span check.
    targCom.D.x = 0.0;
    targCom.D.y = 0.0;
    targCom.D.z = 0.0;
    targCom.Dr.x = 0.0;
    targCom.Dr.y = 0.0;
    targCom.Dr.z = 0.0;

    targCom.trajTrans.x = 4.0;

    reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, 0.0);

    std::cerr << "  Resulting trajTrans = (" << targCom.trajTrans.x << ", "
              << targCom.trajTrans.y << ", " << targCom.trajTrans.z << ")\n"
              << "  Resulting trajRot   = (" << targCom.trajRot.x << ", "
              << targCom.trajRot.y << ", " << targCom.trajRot.z << ")\n";

    // The routine must not leave NaN/inf behind.
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.x)) << "trajTrans.x must be finite";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.y)) << "trajTrans.y must be finite";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.z)) << "trajTrans.z must be finite";

    // Every member (and its interface) must be inside both X walls afterwards.
    const double posWall = membraneObject.waterBox.x / 2.0;
    for (std::size_t i = 0; i < moleculeList.size(); ++i) {
        Coord fin = rtcrrb_final_position(targCom, moleculeList[i].comCoord);
        std::cerr << "  Member " << i << " final x = " << fin.x
                  << " (walls at +/-" << posWall << ")\n";
        EXPECT_LE(fin.x, posWall + 1e-6) << "member " << i << " must be inside the +X wall";
        EXPECT_GE(fin.x, -posWall - 1e-6) << "member " << i << " must be inside the -X wall";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper runs inside its own TEST so a
// failure in one does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(ReflectTrajComplexRadRotBox, InteriorNoChange) { test_rtcrrb_interior_no_change(); }
TEST(ReflectTrajComplexRadRotBox, PositiveXReflection) { test_rtcrrb_positive_x_reflection(); }
TEST(ReflectTrajComplexRadRotBox, NegativeYReflection) { test_rtcrrb_negative_y_reflection(); }
TEST(ReflectTrajComplexRadRotBox, ZWallUsesRS3D) { test_rtcrrb_z_wall_uses_rs3d_when_not_on_surface(); }
TEST(ReflectTrajComplexRadRotBox, ZWallIgnoresRS3DOnSurface) { test_rtcrrb_z_wall_ignores_rs3d_when_on_surface(); }
TEST(ReflectTrajComplexRadRotBox, MemberOffsetDrivesReflection) { test_rtcrrb_member_offset_drives_reflection(); }
TEST(ReflectTrajComplexRadRotBox, RotationIsApplied) { test_rtcrrb_rotation_is_applied(); }
TEST(ReflectTrajComplexRadRotBox, RecheckSpanBranch) { test_rtcrrb_recheck_span_branch(); }