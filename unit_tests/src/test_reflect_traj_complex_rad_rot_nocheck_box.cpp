/*! \file test_reflect_traj_complex_rad_rot_nocheck_box.cpp
 *
 * ### Unit tests for src/boundary_conditions/reflect_traj_complex_rad_rot_nocheck_box.cpp
 *
 * Function under test:
 * \code
 *   void reflect_traj_complex_rad_rot_nocheck_box(const Parameters& params,
 *                                                Complex& targCom,
 *                                                std::vector<Molecule>& moleculeList,
 *                                                const Membrane& membraneObject,
 *                                                double RS3Dinput);
 * \endcode
 *
 * ### What the function does
 * It looks at where a Complex *would* end up after applying its current
 * translational (`trajTrans`) and rotational (`trajRot`) trajectory, and, if any
 * molecule center-of-mass or interface would leave the rectangular water box, it
 * shortens (reflects) `trajTrans` by twice the overshoot in that dimension:
 *
 *     trajTrans.d -= 2 * (farthestPoint_d - wall_d)
 *
 * Notes that drive the tests below:
 *   - Only `trajTrans` is modified. Coordinates (`comCoord`, interface coords)
 *     are never touched, and there is no re-checking after the reflection
 *     ("nocheck").
 *   - A dimension is only examined at all if the sphere of radius
 *     `targCom.radius` around the *translated* complex COM could poke out of
 *     the box (the `canBeOutsideX/Y/Z` short-circuit).
 *   - The negative-Z wall is offset by RS3D: `negZ = -z/2 + RS3D`, where
 *     RS3D == 0 whenever `targCom.OnSurface == true`, otherwise RS3Dinput.
 *   - The rotation matrix built from `trajRot` is applied to the
 *     (member coordinate - complex COM) vectors when locating the extreme
 *     points, so rotation can by itself cause a reflection.
 *
 * All assertions are non-fatal EXPECT_* so every scenario runs even on failure,
 * and each test prints what it is doing and its pass criteria to stderr.
 */

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Small helpers used to build minimal Molecule / Complex / Membrane objects.
// The rtcrrnb_ prefix keeps these unique inside the combined test binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a molecule with a given COM and an explicit list of interfaces.
 *
 * \param[in] com        Center-of-mass coordinate of the molecule.
 * \param[in] ifaceCrds  Absolute interface coordinates (may be empty).
 * \param[in] comIndex   Index of the parent complex (always 0 in these tests).
 */
Molecule rtcrrnb_make_molecule(const Coord& com, const std::vector<Coord>& ifaceCrds, int comIndex = 0)
{
    Molecule mol;
    mol.comCoord = com;
    mol.myComIndex = comIndex;
    mol.interfaceList.clear();
    for (const auto& crd : ifaceCrds) {
        Molecule::Iface iface;
        iface.coord = crd;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Build a complex with a COM, bounding radius, and member indices.
 *
 * Trajectories start at zero; individual tests set them as needed.
 */
Complex rtcrrnb_make_complex(const Coord& com, double radius, const std::vector<int>& members)
{
    Complex targCom;
    targCom.comCoord = com;
    targCom.radius = radius;

    // Diffusion constants are unused by this function, but set them to sane
    // non-zero values so the object resembles a real simulation complex.
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Coord { 0.01, 0.01, 0.01 };

    // No trajectory to begin with (identity rotation, zero translation).
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot = Coord { 0.0, 0.0, 0.0 };

    targCom.memberList = members;
    targCom.OnSurface = false;
    return targCom;
}

/*! \brief Build a rectangular (non-spherical) water box membrane. */
Membrane rtcrrnb_make_box(double x, double y, double z)
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

/*! \brief A Parameters object with just a timestep (nothing else is read). */
Parameters rtcrrnb_make_params()
{
    Parameters params;
    params.timeStep = 1.0;
    return params;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: A complex safely inside the box keeps its trajectory untouched.
// -----------------------------------------------------------------------------
void test_rtcrrnb_interior_no_change()
{
    std::cerr << "\n[TEST] test_rtcrrnb_interior_no_change\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      100^3 box, complex at origin (radius 5) with a\n"
              << "                 tiny translation that cannot reach any wall.\n"
              << "  Pass criteria: trajTrans is bit-for-bit unchanged because all\n"
              << "                 canBeOutsideX/Y/Z tests are false.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    // One molecule sitting on the complex COM, with a coincident interface.
    std::vector<Molecule> moleculeList { rtcrrnb_make_molecule(Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 0.0, 0.0 } }) };
    Complex targCom = rtcrrnb_make_complex(Coord { 0.0, 0.0, 0.0 }, 5.0, { 0 });

    targCom.trajTrans.x = 1.0;
    targCom.trajTrans.y = -2.0;
    targCom.trajTrans.z = 0.5;

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    std::cerr << "  Resulting trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y << ", "
              << targCom.trajTrans.z << ")\n";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 1.0) << "interior complex: trajTrans.x must not change";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, -2.0) << "interior complex: trajTrans.y must not change";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.5) << "interior complex: trajTrans.z must not change";
}

// -----------------------------------------------------------------------------
// Test 2: Overshoot of the +X wall is reflected by exactly 2 * overshoot.
// -----------------------------------------------------------------------------
void test_rtcrrnb_reflect_pos_x()
{
    std::cerr << "\n[TEST] test_rtcrrnb_reflect_pos_x\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      100^3 box (walls at +/-50), complex COM at x=45,\n"
              << "                 radius 10, trajTrans.x = +20 => would land at 65.\n"
              << "  Pass criteria: overshoot = 15, so trajTrans.x = 20 - 2*15 = -10,\n"
              << "                 and Y/Z trajectories are untouched.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrnb_make_molecule(Coord { 45.0, 0.0, 0.0 }, { Coord { 45.0, 0.0, 0.0 } }) };
    Complex targCom = rtcrrnb_make_complex(Coord { 45.0, 0.0, 0.0 }, 10.0, { 0 });
    targCom.trajTrans.x = 20.0;

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    const double finalX = targCom.comCoord.x + targCom.trajTrans.x;
    std::cerr << "  trajTrans.x = " << targCom.trajTrans.x << " => final x = " << finalX << " (wall at +50)\n";

    EXPECT_NEAR(targCom.trajTrans.x, -10.0, 1e-12) << "trajTrans.x should be reflected to -10 (20 - 2*15)";
    EXPECT_LE(finalX, 50.0 + 1e-9) << "reflected complex must end up inside the +X wall";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "Y trajectory must be untouched";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0) << "Z trajectory must be untouched";
}

// -----------------------------------------------------------------------------
// Test 3: Overshoot of the -X wall is reflected in the opposite direction.
// -----------------------------------------------------------------------------
void test_rtcrrnb_reflect_neg_x()
{
    std::cerr << "\n[TEST] test_rtcrrnb_reflect_neg_x\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      complex COM at x=-45, radius 10, trajTrans.x = -20\n"
              << "                 => would land at -65 (wall at -50).\n"
              << "  Pass criteria: overshoot = -15, so trajTrans.x = -20 - 2*(-15) = +10.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList {
        rtcrrnb_make_molecule(Coord { -45.0, 0.0, 0.0 }, { Coord { -45.0, 0.0, 0.0 } })
    };
    Complex targCom = rtcrrnb_make_complex(Coord { -45.0, 0.0, 0.0 }, 10.0, { 0 });
    targCom.trajTrans.x = -20.0;

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    const double finalX = targCom.comCoord.x + targCom.trajTrans.x;
    std::cerr << "  trajTrans.x = " << targCom.trajTrans.x << " => final x = " << finalX << " (wall at -50)\n";

    EXPECT_NEAR(targCom.trajTrans.x, 10.0, 1e-12) << "trajTrans.x should be reflected to +10";
    EXPECT_GE(finalX, -50.0 - 1e-9) << "reflected complex must end up inside the -X wall";
}

// -----------------------------------------------------------------------------
// Test 4: The Y branch works on a non-cubic box (independent wall positions).
// -----------------------------------------------------------------------------
void test_rtcrrnb_reflect_y_noncubic_box()
{
    std::cerr << "\n[TEST] test_rtcrrnb_reflect_y_noncubic_box\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      box 100 x 60 x 100 (Y walls at +/-30), complex COM\n"
              << "                 at y=25, radius 10, trajTrans.y = +10 => y=35.\n"
              << "  Pass criteria: overshoot = 5 => trajTrans.y = 10 - 2*5 = 0, and the\n"
              << "                 X/Z dimensions are never even examined.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 60.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrnb_make_molecule(Coord { 0.0, 25.0, 0.0 }, { Coord { 0.0, 25.0, 0.0 } }) };
    Complex targCom = rtcrrnb_make_complex(Coord { 0.0, 25.0, 0.0 }, 10.0, { 0 });
    targCom.trajTrans.x = 3.0; // safe in X (walls at +/-50)
    targCom.trajTrans.y = 10.0; // violates +Y wall
    targCom.trajTrans.z = -4.0; // safe in Z

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    const double finalY = targCom.comCoord.y + targCom.trajTrans.y;
    std::cerr << "  trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y << ", "
              << targCom.trajTrans.z << ") => final y = " << finalY << " (wall at +30)\n";

    EXPECT_NEAR(targCom.trajTrans.y, 0.0, 1e-12) << "trajTrans.y should be reflected to 0 (10 - 2*5)";
    EXPECT_LE(finalY, 30.0 + 1e-9) << "reflected complex must end up inside the +Y wall";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 3.0) << "X trajectory should not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, -4.0) << "Z trajectory should not be modified";
}

// -----------------------------------------------------------------------------
// Test 5: The -Z wall is raised by RS3D when the complex is NOT on the surface.
// -----------------------------------------------------------------------------
void test_rtcrrnb_negz_uses_rs3d_offset()
{
    std::cerr << "\n[TEST] test_rtcrrnb_negz_uses_rs3d_offset\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      100^3 box, OnSurface = false, RS3Dinput = 5 so the\n"
              << "                 effective -Z wall is at -45. Complex COM z=-40,\n"
              << "                 radius 10, trajTrans.z = -10 => z = -50.\n"
              << "  Pass criteria: overshoot = -5 => trajTrans.z = -10 - 2*(-5) = 0,\n"
              << "                 i.e. the molecule stops at the reflecting surface.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList {
        rtcrrnb_make_molecule(Coord { 0.0, 0.0, -40.0 }, { Coord { 0.0, 0.0, -40.0 } })
    };
    Complex targCom = rtcrrnb_make_complex(Coord { 0.0, 0.0, -40.0 }, 10.0, { 0 });
    targCom.OnSurface = false; // => RS3D = RS3Dinput
    targCom.trajTrans.z = -10.0;

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 5.0);

    const double finalZ = targCom.comCoord.z + targCom.trajTrans.z;
    std::cerr << "  trajTrans.z = " << targCom.trajTrans.z << " => final z = " << finalZ
              << " (effective -Z wall at -45)\n";

    EXPECT_NEAR(targCom.trajTrans.z, 0.0, 1e-12) << "trajTrans.z should be reflected to 0 given negZ = -45";
    EXPECT_GE(finalZ, -45.0 - 1e-9) << "final z must be at or above the RS3D-shifted wall";
}

// -----------------------------------------------------------------------------
// Test 6: When the complex is on the surface, RS3D is forced to zero.
// -----------------------------------------------------------------------------
void test_rtcrrnb_onsurface_ignores_rs3d()
{
    std::cerr << "\n[TEST] test_rtcrrnb_onsurface_ignores_rs3d\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      identical geometry to the previous test, but\n"
              << "                 OnSurface = true so RS3D is forced to 0 and the\n"
              << "                 -Z wall sits at the true box edge, -50.\n"
              << "  Pass criteria: z lands exactly on -50, no overshoot, so\n"
              << "                 trajTrans.z stays at -10.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList {
        rtcrrnb_make_molecule(Coord { 0.0, 0.0, -40.0 }, { Coord { 0.0, 0.0, -40.0 } })
    };
    Complex targCom = rtcrrnb_make_complex(Coord { 0.0, 0.0, -40.0 }, 10.0, { 0 });
    targCom.OnSurface = true; // => RS3D = 0 regardless of RS3Dinput
    targCom.trajTrans.z = -10.0;

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 5.0);

    std::cerr << "  trajTrans.z = " << targCom.trajTrans.z << " (expected unchanged -10)\n";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, -10.0)
        << "OnSurface complexes must ignore RS3Dinput; -50 is exactly on the wall, so no reflection";
}

// -----------------------------------------------------------------------------
// Test 7: An interface, not the COM, can be the point that triggers reflection.
// -----------------------------------------------------------------------------
void test_rtcrrnb_interface_drives_reflection()
{
    std::cerr << "\n[TEST] test_rtcrrnb_interface_drives_reflection\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      COM at x=40 (inside), radius 12 so the dimension is\n"
              << "                 examined, and an interface 12 nm out at x=52 which is\n"
              << "                 past the +50 wall. No translation applied.\n"
              << "  Pass criteria: overshoot of the interface = 2 => trajTrans.x = -4,\n"
              << "                 pulling the interface back to exactly x = 48.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    // Molecule COM at 40, interface sticking out to 52 (offset +12 from complex COM).
    std::vector<Molecule> moleculeList { rtcrrnb_make_molecule(Coord { 40.0, 0.0, 0.0 }, { Coord { 52.0, 0.0, 0.0 } }) };
    Complex targCom = rtcrrnb_make_complex(Coord { 40.0, 0.0, 0.0 }, 12.0, { 0 });

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    const double finalIfaceX = 52.0 + targCom.trajTrans.x;
    std::cerr << "  trajTrans.x = " << targCom.trajTrans.x << " => final interface x = " << finalIfaceX << "\n";

    EXPECT_NEAR(targCom.trajTrans.x, -4.0, 1e-12) << "trajTrans.x should be -2*(52-50) = -4";
    EXPECT_LE(finalIfaceX, 50.0 + 1e-9) << "the offending interface must be brought back inside the box";
}

// -----------------------------------------------------------------------------
// Test 8: The farthest member of a multi-molecule complex sets the reflection.
// -----------------------------------------------------------------------------
void test_rtcrrnb_farthest_member_wins()
{
    std::cerr << "\n[TEST] test_rtcrrnb_farthest_member_wins\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      two-molecule complex, COM at origin, radius 60 so all\n"
              << "                 dimensions get examined. Members at x=20 (inside) and\n"
              << "                 x=56 (outside the +50 wall).\n"
              << "  Pass criteria: the x=56 member sets posWallX, so trajTrans.x = -12,\n"
              << "                 while Y and Z remain zero (nothing sticks out there).\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList {
        rtcrrnb_make_molecule(Coord { 20.0, 0.0, 0.0 }, { Coord { 20.0, 0.0, 0.0 } }, 0),
        rtcrrnb_make_molecule(Coord { 56.0, 0.0, 0.0 }, { Coord { 56.0, 0.0, 0.0 } }, 0)
    };
    Complex targCom = rtcrrnb_make_complex(Coord { 0.0, 0.0, 0.0 }, 60.0, { 0, 1 });

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    std::cerr << "  trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y << ", "
              << targCom.trajTrans.z << ")\n";

    EXPECT_NEAR(targCom.trajTrans.x, -12.0, 1e-12) << "reflection must be driven by the farthest member (-2*(56-50))";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "no Y overshoot exists, so trajTrans.y stays 0";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0) << "no Z overshoot exists, so trajTrans.z stays 0";
}

// -----------------------------------------------------------------------------
// Test 9: Rotation alone (no translation) can push a member out and cause a
//         reflection, proving trajRot is folded into the extreme-point search.
// -----------------------------------------------------------------------------
void test_rtcrrnb_rotation_causes_reflection()
{
    std::cerr << "\n[TEST] test_rtcrrnb_rotation_causes_reflection\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      COM at x=-44 with an interface offset of +10 (x=-34),\n"
              << "                 zero translation but trajRot = (0,0,pi). The pi rotation\n"
              << "                 flips the +10 offset to -10, putting the interface at -54,\n"
              << "                 outside the -50 wall.\n"
              << "  Pass criteria: trajTrans.x becomes +8 (= -2*(-54 - -50)), so the\n"
              << "                 rotated interface ends up exactly on the wall.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    // Sanity note printed for the reader: with trajRot = (0,0,pi) the x-row of the
    // Euler matrix maps a pure +x offset onto -x (cos(pi) = -1).
    std::array<double, 9> M = create_euler_rotation_matrix(Coord { 0.0, 0.0, M_PI });
    std::cerr << "  Euler matrix x-row for (0,0,pi) = [" << M[0] << ", " << M[1] << ", " << M[2] << "]\n";

    std::vector<Molecule> moleculeList {
        rtcrrnb_make_molecule(Coord { -44.0, 0.0, 0.0 }, { Coord { -34.0, 0.0, 0.0 } })
    };
    Complex targCom = rtcrrnb_make_complex(Coord { -44.0, 0.0, 0.0 }, 12.0, { 0 });
    targCom.trajRot = Coord { 0.0, 0.0, M_PI };
    // No translation at all -- any change must come from the rotation.

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    // Rotated interface position (offset +10 -> -10 relative to COM), plus the
    // reflection that the function just computed.
    const double rotatedIfaceX = targCom.comCoord.x + targCom.trajTrans.x + (M[0] * 10.0);
    std::cerr << "  trajTrans.x = " << targCom.trajTrans.x << " => rotated interface x = " << rotatedIfaceX
              << " (wall at -50)\n";

    EXPECT_NEAR(targCom.trajTrans.x, 8.0, 1e-9) << "rotation-induced overshoot of 4 nm must give trajTrans.x = +8";
    EXPECT_GE(rotatedIfaceX, -50.0 - 1e-9) << "the rotated interface must be brought back inside the -X wall";
}

// -----------------------------------------------------------------------------
// Test 10: The function must not modify any coordinates, only trajTrans.
// -----------------------------------------------------------------------------
void test_rtcrrnb_coordinates_untouched()
{
    std::cerr << "\n[TEST] test_rtcrrnb_coordinates_untouched\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      a strongly out-of-box move is reflected; afterwards we\n"
              << "                 re-inspect the complex COM, the molecule COM, the\n"
              << "                 interface coordinate and trajRot.\n"
              << "  Pass criteria: only trajTrans changed; every coordinate and the\n"
              << "                 rotation vector are identical to their input values.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList { rtcrrnb_make_molecule(Coord { 48.0, 0.0, 0.0 }, { Coord { 49.0, 1.0, 2.0 } }) };
    Complex targCom = rtcrrnb_make_complex(Coord { 48.0, 0.0, 0.0 }, 10.0, { 0 });
    targCom.trajRot = Coord { 0.0, 0.0, 0.0 };
    targCom.trajTrans.x = 30.0; // grossly out of the box

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    std::cerr << "  trajTrans.x after reflection = " << targCom.trajTrans.x << "\n"
              << "  complex COM  = (" << targCom.comCoord.x << ", " << targCom.comCoord.y << ", " << targCom.comCoord.z
              << ")\n"
              << "  molecule COM = (" << moleculeList[0].comCoord.x << ", " << moleculeList[0].comCoord.y << ", "
              << moleculeList[0].comCoord.z << ")\n";

    // The trajectory must have been shortened (proof the reflection happened).
    EXPECT_LT(targCom.trajTrans.x, 30.0) << "an out-of-box move must be shortened by the reflection";

    // Coordinates are read-only for this routine.
    EXPECT_DOUBLE_EQ(targCom.comCoord.x, 48.0) << "complex COM x must not be modified";
    EXPECT_DOUBLE_EQ(targCom.comCoord.y, 0.0) << "complex COM y must not be modified";
    EXPECT_DOUBLE_EQ(targCom.comCoord.z, 0.0) << "complex COM z must not be modified";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 48.0) << "molecule COM x must not be modified";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, 49.0) << "interface x must not be modified";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.y, 1.0) << "interface y must not be modified";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.z, 2.0) << "interface z must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, 0.0) << "trajRot must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, 0.0) << "trajRot must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, 0.0) << "trajRot must not be modified";
}

// -----------------------------------------------------------------------------
// Test 11: A complex thin in a dimension but with a big radius: the radius-based
//          pre-check gates whether a dimension is even inspected.
// -----------------------------------------------------------------------------
void test_rtcrrnb_radius_gate_skips_dimension()
{
    std::cerr << "\n[TEST] test_rtcrrnb_radius_gate_skips_dimension\n"
              << "  Source:        reflect_traj_complex_rad_rot_nocheck_box.cpp\n"
              << "  Function:      reflect_traj_complex_rad_rot_nocheck_box\n"
              << "  Scenario:      an interface actually sits outside the +X wall (x=60),\n"
              << "                 but the complex radius is tiny (0.1) and its COM is at\n"
              << "                 the origin, so canBeOutsideX evaluates to false.\n"
              << "  Pass criteria: no reflection happens -- this documents the 'nocheck'\n"
              << "                 short-circuit behaviour of the routine.\n";

    Parameters params = rtcrrnb_make_params();
    Membrane membraneObject = rtcrrnb_make_box(100.0, 100.0, 100.0);

    // Deliberately inconsistent geometry: radius does not enclose the interface.
    std::vector<Molecule> moleculeList { rtcrrnb_make_molecule(Coord { 0.0, 0.0, 0.0 }, { Coord { 60.0, 0.0, 0.0 } }) };
    Complex targCom = rtcrrnb_make_complex(Coord { 0.0, 0.0, 0.0 }, 0.1, { 0 });
    targCom.trajTrans.x = 1.0;

    reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, 0.0);

    std::cerr << "  trajTrans.x = " << targCom.trajTrans.x << " (expected unchanged 1.0)\n";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 1.0)
        << "with a radius too small to reach the wall the X dimension is skipped entirely";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each scenario is its own reported test case, and all of
// them run even if earlier ones fail (only non-fatal EXPECT_* used above).
// -----------------------------------------------------------------------------
TEST(ReflectTrajComplexRadRotNocheckBox, InteriorNoChange) { test_rtcrrnb_interior_no_change(); }
TEST(ReflectTrajComplexRadRotNocheckBox, ReflectPosX) { test_rtcrrnb_reflect_pos_x(); }
TEST(ReflectTrajComplexRadRotNocheckBox, ReflectNegX) { test_rtcrrnb_reflect_neg_x(); }
TEST(ReflectTrajComplexRadRotNocheckBox, ReflectYNonCubicBox) { test_rtcrrnb_reflect_y_noncubic_box(); }
TEST(ReflectTrajComplexRadRotNocheckBox, NegZUsesRS3DOffset) { test_rtcrrnb_negz_uses_rs3d_offset(); }
TEST(ReflectTrajComplexRadRotNocheckBox, OnSurfaceIgnoresRS3D) { test_rtcrrnb_onsurface_ignores_rs3d(); }
TEST(ReflectTrajComplexRadRotNocheckBox, InterfaceDrivesReflection) { test_rtcrrnb_interface_drives_reflection(); }
TEST(ReflectTrajComplexRadRotNocheckBox, FarthestMemberWins) { test_rtcrrnb_farthest_member_wins(); }
TEST(ReflectTrajComplexRadRotNocheckBox, RotationCausesReflection) { test_rtcrrnb_rotation_causes_reflection(); }
TEST(ReflectTrajComplexRadRotNocheckBox, CoordinatesUntouched) { test_rtcrrnb_coordinates_untouched(); }
TEST(ReflectTrajComplexRadRotNocheckBox, RadiusGateSkipsDimension) { test_rtcrrnb_radius_gate_skips_dimension(); }