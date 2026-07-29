/*! \file test_reflect_traj_check_span_box.cpp
 *
 * ### Unit test for src/boundary_conditions/reflect_traj_check_span_box.cpp
 *
 * Function under test:
 *
 *     void reflect_traj_check_span_box(const Parameters& params,
 *                                      Complex& targCom,
 *                                      std::vector<Molecule>& moleculeList,
 *                                      const Membrane& membraneObject,
 *                                      double RS3Dinput)
 *
 * What the function does (and therefore what we verify):
 *   - It walks over every member Molecule of `targCom`, and over every
 *     interface of every member Molecule, and rotates/translates those points
 *     by the complex's proposed move (`trajRot`, `trajTrans`).
 *   - From those points it finds the extreme (farthest) coordinate in each of
 *     the +X/-X, +Y/-Y, +Z/-Z directions.
 *   - If an extreme point would end up outside a wall of the rectangular water
 *     box, the corresponding component of `targCom.trajTrans` is *reflected*:
 *         trajTrans -= 2 * (overshoot)
 *   - The -Z wall is offset by RS3D (the 3D reflecting surface), but only when
 *     the complex is NOT already on the surface (`targCom.OnSurface == false`).
 *   - If the complex sticks out of *both* walls of a dimension (i.e. it spans
 *     the box), the move is declared failed, a brand new random move is drawn
 *     (GaussV/GSL RNG), and the check is repeated up to a bounded number of
 *     iterations.
 *   - The function only ever modifies `targCom.trajTrans` / `targCom.trajRot`;
 *     it must not move the stored coordinates of the complex or its molecules.
 *
 * Everything below uses a zero rotation vector so the internal Euler rotation
 * matrix is exactly the identity, which makes the reflected translation values
 * exactly predictable.
 */

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/rand_gsl.hpp"

#include <gtest/gtest.h>

// The GSL random number generator is defined in unit_tests/src/gtest_main.cpp.
// It is declared extern by math/rand_gsl.hpp; we only need to make sure it is
// allocated before we exercise the code path that resamples a move.

namespace {

// -----------------------------------------------------------------------------
// Helper: make sure the GSL RNG is usable.
//
// The "move failed" branch of the function under test calls GaussV(), which
// dereferences the global gsl_rng pointer `r`.  gtest_main.cpp defines it as
// nullptr, so if no other test has seeded it yet we do it here.
// -----------------------------------------------------------------------------
void rtcsb_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  (initialising GSL RNG with srand_gsl(1) for resampling path)\n";
        const gsl_rng_type *T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 1);
    }
}

/*! \brief Build a minimal Molecule belonging to complex 0.
 *
 * \param[in] com          center-of-mass coordinate of the molecule
 * \param[in] ifaceOffset  offset (added to com) of the single interface
 *
 * A single interface keeps the geometry easy to reason about: the extreme
 * point of the complex is either the molecule COM or this one interface.
 */
Molecule rtcsb_make_molecule(const Coord& com, const Coord& ifaceOffset = Coord { 0.0, 0.0, 0.0 })
{
    Molecule mol;
    mol.comCoord = com;
    mol.myComIndex = 0;

    Molecule::Iface iface;
    iface.coord = Coord { com.x + ifaceOffset.x, com.y + ifaceOffset.y, com.z + ifaceOffset.z };
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a Complex centered at `com` owning `nMembers` molecules.
 *
 * The member indices are 0 .. nMembers-1 so they line up with a moleculeList
 * built in the same order.  Diffusion constants are non-zero so that, should
 * the move need to be resampled, GaussV() produces a meaningful displacement.
 */
Complex rtcsb_make_complex(const Coord& com, int nMembers)
{
    Complex targCom;
    targCom.index = 0;
    targCom.comCoord = com;

    // Translational / rotational diffusion constants (only used on resample).
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Coord { 0.01, 0.01, 0.01 };

    // Start with no proposed move at all; each test sets what it needs.
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot = Coord { 0.0, 0.0, 0.0 }; // identity rotation matrix

    targCom.memberList.clear();
    for (int i = 0; i < nMembers; ++i)
        targCom.memberList.push_back(i);

    targCom.OnSurface = false;

    return targCom;
}

/*! \brief A cubic water box of side `side` centered on the origin. */
Membrane rtcsb_make_box(double side)
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.waterBox.x = side;
    membraneObject.waterBox.y = side;
    membraneObject.waterBox.z = side;
    membraneObject.waterBox.volume = side * side * side;
    return membraneObject;
}

/*! \brief Simple Parameters object; only timeStep is consulted. */
Parameters rtcsb_make_params()
{
    Parameters params;
    params.timeStep = 1.0;
    return params;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: complex deep inside the box.
// Pass criteria: trajTrans and trajRot untouched, stored coordinates untouched.
// -----------------------------------------------------------------------------
void test_rtcsb_interior_move_untouched()
{
    std::cerr << "\n[TEST] test_rtcsb_interior_move_untouched\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    single molecule at the origin of a 100 nm box with a\n"
              << "               1 nm-scale step; nothing comes near a wall.\n"
              << "  Criteria:    trajTrans/trajRot unchanged, coordinates unchanged.\n";

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);

    Complex targCom = rtcsb_make_complex(Coord { 0.0, 0.0, 0.0 }, 1);
    std::vector<Molecule> moleculeList { rtcsb_make_molecule(Coord { 0.0, 0.0, 0.0 }) };

    // A small, definitely-safe displacement.
    targCom.trajTrans.x = 1.0;
    targCom.trajTrans.y = -2.0;
    targCom.trajTrans.z = 0.5;

    std::cerr << "  Proposed trajTrans = (1, -2, 0.5)\n";

    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, 0.0);

    // No reflection should have been applied.
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 1.0) << "trajTrans.x must be untouched inside the box";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, -2.0) << "trajTrans.y must be untouched inside the box";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.5) << "trajTrans.z must be untouched inside the box";

    // The rotation vector must also be left alone (no resampling happened).
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, 0.0) << "trajRot.x must be untouched";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, 0.0) << "trajRot.y must be untouched";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, 0.0) << "trajRot.z must be untouched";

    // The function must not commit the move: stored coordinates are unchanged.
    EXPECT_DOUBLE_EQ(targCom.comCoord.x, 0.0) << "complex COM must not be moved by this function";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 0.0) << "molecule COM must not be moved by this function";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, 0.0)
        << "interface coordinate must not be moved by this function";

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y << ", "
              << targCom.trajTrans.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: reflection off the +X wall.
//
// com.x = 45, step +20 -> farthest point 65, wall at +50, overshoot 15.
// Expected reflected step: 20 - 2*15 = -10, final x = 35.
// -----------------------------------------------------------------------------
void test_rtcsb_reflect_positive_x_wall()
{
    std::cerr << "\n[TEST] test_rtcsb_reflect_positive_x_wall\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    COM at x=45 in a 100 nm box, proposed step +20 nm.\n"
              << "  Criteria:    trajTrans.x becomes 20 - 2*(65-50) = -10, and the\n"
              << "               resulting position (35) is inside the +X wall.\n";

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);
    const double posX = membraneObject.waterBox.x / 2.0;

    Complex targCom = rtcsb_make_complex(Coord { 45.0, 0.0, 0.0 }, 1);
    std::vector<Molecule> moleculeList { rtcsb_make_molecule(Coord { 45.0, 0.0, 0.0 }) };

    targCom.trajTrans.x = 20.0;

    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, 0.0);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, -10.0) << "trajTrans.x should be reflected to -10";
    // Untouched dimensions must remain zero.
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "trajTrans.y should stay 0 (no Y violation)";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0) << "trajTrans.z should stay 0 (no Z violation)";

    const double finalX = targCom.comCoord.x + targCom.trajTrans.x;
    std::cerr << "  Final X = " << finalX << " (wall at " << posX << ")\n";
    EXPECT_LE(finalX, posX + 1e-9) << "reflected complex must end up inside the +X wall";
}

// -----------------------------------------------------------------------------
// Test 3: reflection off the -Y wall.
//
// com.y = -45, step -20 -> farthest point -65, wall at -50, overshoot -15.
// Expected reflected step: -20 - 2*(-15) = +10, final y = -35.
// -----------------------------------------------------------------------------
void test_rtcsb_reflect_negative_y_wall()
{
    std::cerr << "\n[TEST] test_rtcsb_reflect_negative_y_wall\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    COM at y=-45 in a 100 nm box, proposed step -20 nm.\n"
              << "  Criteria:    trajTrans.y becomes -20 - 2*(-65+50) = +10 and the\n"
              << "               result (-35) is inside the -Y wall.\n";

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);
    const double negY = -membraneObject.waterBox.y / 2.0;

    Complex targCom = rtcsb_make_complex(Coord { 0.0, -45.0, 0.0 }, 1);
    std::vector<Molecule> moleculeList { rtcsb_make_molecule(Coord { 0.0, -45.0, 0.0 }) };

    targCom.trajTrans.y = -20.0;

    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, 0.0);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 10.0) << "trajTrans.y should be reflected to +10";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 0.0) << "trajTrans.x should stay 0 (no X violation)";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0) << "trajTrans.z should stay 0 (no Z violation)";

    const double finalY = targCom.comCoord.y + targCom.trajTrans.y;
    std::cerr << "  Final Y = " << finalY << " (wall at " << negY << ")\n";
    EXPECT_GE(finalY, negY - 1e-9) << "reflected complex must end up inside the -Y wall";
}

// -----------------------------------------------------------------------------
// Test 4: the interface coordinates - not just the molecule COM - define the
//         extreme point used for the reflection.
//
// COM at x=40 with its interface 9 nm further out in +X; step +10.
//   using the COM only      -> overshoot 0   -> no reflection
//   using the interface too -> point 59, overshoot 9 -> step 10-18 = -8
// -----------------------------------------------------------------------------
void test_rtcsb_interface_defines_extreme_point()
{
    std::cerr << "\n[TEST] test_rtcsb_interface_defines_extreme_point\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    molecule COM at x=40 (safe) but its interface at\n"
              << "               x=49; proposed step +10 nm pushes the interface to 59.\n"
              << "  Criteria:    trajTrans.x = 10 - 2*(59-50) = -8, proving interfaces\n"
              << "               are included in the span check.\n";

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);
    const double posX = membraneObject.waterBox.x / 2.0;

    Complex targCom = rtcsb_make_complex(Coord { 40.0, 0.0, 0.0 }, 1);
    std::vector<Molecule> moleculeList {
        rtcsb_make_molecule(Coord { 40.0, 0.0, 0.0 }, Coord { 9.0, 0.0, 0.0 })
    };

    targCom.trajTrans.x = 10.0;

    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, 0.0);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, -8.0)
        << "the interface, not the COM, should set the overshoot (expected -8)";

    const double finalIfaceX = moleculeList[0].interfaceList[0].coord.x + targCom.trajTrans.x;
    std::cerr << "  Final interface X = " << finalIfaceX << " (wall at " << posX << ")\n";
    EXPECT_LE(finalIfaceX, posX + 1e-9) << "the outermost interface must end up inside the +X wall";
}

// -----------------------------------------------------------------------------
// Test 5: multi-molecule complex - the reflection is driven by whichever member
//         sticks out farthest.
//
// Members at x = 30 and x = 48, step +10 -> farthest 58, overshoot 8,
// expected step 10 - 16 = -6.
// -----------------------------------------------------------------------------
void test_rtcsb_multi_member_uses_farthest_member()
{
    std::cerr << "\n[TEST] test_rtcsb_multi_member_uses_farthest_member\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    two-molecule complex (x=30 and x=48), step +10 nm.\n"
              << "  Criteria:    reflection uses the farthest member (48+10=58), so\n"
              << "               trajTrans.x = 10 - 2*(58-50) = -6, and both members\n"
              << "               finish inside the +X wall.\n";

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);
    const double posX = membraneObject.waterBox.x / 2.0;

    // Complex COM placed between the two molecules; its exact value only
    // matters through the (molecule - complex) vectors, which are exact here.
    Complex targCom = rtcsb_make_complex(Coord { 39.0, 0.0, 0.0 }, 2);
    std::vector<Molecule> moleculeList {
        rtcsb_make_molecule(Coord { 30.0, 0.0, 0.0 }),
        rtcsb_make_molecule(Coord { 48.0, 0.0, 0.0 })
    };
    moleculeList[1].myComIndex = 0;

    targCom.trajTrans.x = 10.0;

    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, 0.0);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, -6.0)
        << "reflection must be based on the farthest member (expected -6)";

    for (std::size_t i = 0; i < moleculeList.size(); ++i) {
        const double finalX = moleculeList[i].comCoord.x + targCom.trajTrans.x;
        std::cerr << "  member " << i << " final X = " << finalX << '\n';
        EXPECT_LE(finalX, posX + 1e-9) << "member " << i << " must end up inside the +X wall";
    }
}

// -----------------------------------------------------------------------------
// Test 6: the -Z wall is raised by RS3Dinput when the complex is in solution.
//
// Box z in [-50, 50]; RS3Dinput = 5 -> effective floor at -45.
// COM at z=-40, step -10 -> point -50, overshoot -5 -> step -10 + 10 = 0.
// -----------------------------------------------------------------------------
void test_rtcsb_rs3d_raises_negative_z_wall()
{
    std::cerr << "\n[TEST] test_rtcsb_rs3d_raises_negative_z_wall\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    OnSurface == false, RS3Dinput = 5, so the reflecting\n"
              << "               floor sits at z = -45 instead of -50. COM z=-40 with\n"
              << "               a -10 nm step would reach -50.\n"
              << "  Criteria:    trajTrans.z = -10 - 2*(-50 + 45) = 0, i.e. the move is\n"
              << "               reflected off the raised surface.\n";

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);
    const double RS3Dinput = 5.0;
    const double negZeff = -membraneObject.waterBox.z / 2.0 + RS3Dinput;

    Complex targCom = rtcsb_make_complex(Coord { 0.0, 0.0, -40.0 }, 1);
    targCom.OnSurface = false; // in solution -> RS3D is applied
    std::vector<Molecule> moleculeList { rtcsb_make_molecule(Coord { 0.0, 0.0, -40.0 }) };

    targCom.trajTrans.z = -10.0;

    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, RS3Dinput);

    EXPECT_NEAR(targCom.trajTrans.z, 0.0, 1e-12)
        << "trajTrans.z should be reflected to 0 off the RS3D-raised floor";

    const double finalZ = targCom.comCoord.z + targCom.trajTrans.z;
    std::cerr << "  Final Z = " << finalZ << " (effective floor at " << negZeff << ")\n";
    EXPECT_GE(finalZ, negZeff - 1e-9) << "complex must stay above the RS3D reflecting surface";
}

// -----------------------------------------------------------------------------
// Test 7: for a complex already on the surface, RS3D is ignored (RS3D = 0), so
//         the -Z wall is the real box wall at -50.
//
// COM at z=-48 with no step: outside the raised floor (-45) but inside the real
// wall (-50), so the move must NOT be modified.
// -----------------------------------------------------------------------------
void test_rtcsb_onsurface_ignores_rs3d()
{
    std::cerr << "\n[TEST] test_rtcsb_onsurface_ignores_rs3d\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    OnSurface == true with RS3Dinput = 5. COM at z=-48 is\n"
              << "               below the raised floor (-45) but above the true wall.\n"
              << "  Criteria:    RS3D must be ignored, so trajTrans stays exactly 0.\n";

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);

    Complex targCom = rtcsb_make_complex(Coord { 0.0, 0.0, -48.0 }, 1);
    targCom.OnSurface = true; // on the membrane -> RS3D forced to 0
    std::vector<Molecule> moleculeList { rtcsb_make_molecule(Coord { 0.0, 0.0, -48.0 }) };

    // No proposed motion at all.
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;

    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, 5.0);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0)
        << "an on-surface complex must not be reflected off the RS3D offset";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 0.0) << "trajTrans.x must stay 0";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "trajTrans.y must stay 0";

    std::cerr << "  Result trajTrans.z = " << targCom.trajTrans.z << " (unchanged, as required)\n";
}

// -----------------------------------------------------------------------------
// Test 8: a complex that spans the whole box in X triggers the "move failed"
//         path: the trajectory is resampled from the RNG (which needs the GSL
//         generator), and the routine must still terminate with finite values.
// -----------------------------------------------------------------------------
void test_rtcsb_spanning_complex_resamples_and_terminates()
{
    std::cerr << "\n[TEST] test_rtcsb_spanning_complex_resamples_and_terminates\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    a two-molecule complex reaching x=-60 and x=+60 in a\n"
              << "               100 nm box sticks out of BOTH X walls, which forces\n"
              << "               the resampling branch (GaussV + nocheck reflection).\n"
              << "  Criteria:    the call terminates (bounded iteration count) and all\n"
              << "               trajectory components remain finite numbers.\n";

    rtcsb_ensure_rng();

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);

    Complex targCom = rtcsb_make_complex(Coord { 0.0, 0.0, 0.0 }, 2);
    std::vector<Molecule> moleculeList {
        rtcsb_make_molecule(Coord { -60.0, 0.0, 0.0 }),
        rtcsb_make_molecule(Coord { 60.0, 0.0, 0.0 })
    };
    moleculeList[1].myComIndex = 0;

    // Any step at all: the complex is already outside both X walls.
    targCom.trajTrans.x = 1.0;

    std::cerr << "  Calling reflect_traj_check_span_box (expect internal resampling)...\n";
    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, 0.0);
    std::cerr << "  Returned without hanging.\n";

    // The resampled move must be a usable (finite) number in every component.
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.x)) << "trajTrans.x must remain finite after resampling";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.y)) << "trajTrans.y must remain finite after resampling";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.z)) << "trajTrans.z must remain finite after resampling";
    EXPECT_TRUE(std::isfinite(targCom.trajRot.x)) << "trajRot.x must remain finite after resampling";
    EXPECT_TRUE(std::isfinite(targCom.trajRot.y)) << "trajRot.y must remain finite after resampling";
    EXPECT_TRUE(std::isfinite(targCom.trajRot.z)) << "trajRot.z must remain finite after resampling";

    // Even in the failure case the stored coordinates must not be committed.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, -60.0) << "molecule 0 COM must not be modified";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 60.0) << "molecule 1 COM must not be modified";

    std::cerr << "  Result trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y << ", "
              << targCom.trajTrans.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 9: simultaneous violations in X, Y and Z are each reflected independently.
//
// COM at (45, -45, 45); step (+20, -20, +20):
//   X: 65 -> 20 - 2*15 = -10   (final  35)
//   Y: -65 -> -20 + 2*15 = +10 (final -35)
//   Z: 65 -> 20 - 2*15 = -10   (final  35)
// -----------------------------------------------------------------------------
void test_rtcsb_reflects_all_three_dimensions()
{
    std::cerr << "\n[TEST] test_rtcsb_reflects_all_three_dimensions\n"
              << "  Source file: reflect_traj_check_span_box.cpp\n"
              << "  Function:    reflect_traj_check_span_box\n"
              << "  Scenario:    a corner molecule at (45,-45,45) stepping (+20,-20,+20)\n"
              << "               violates three different walls at once.\n"
              << "  Criteria:    each component is reflected independently to\n"
              << "               (-10, +10, -10) and the final point is inside the box.\n";

    Parameters params = rtcsb_make_params();
    Membrane membraneObject = rtcsb_make_box(100.0);
    const double half = membraneObject.waterBox.x / 2.0;

    Complex targCom = rtcsb_make_complex(Coord { 45.0, -45.0, 45.0 }, 1);
    std::vector<Molecule> moleculeList { rtcsb_make_molecule(Coord { 45.0, -45.0, 45.0 }) };

    targCom.trajTrans.x = 20.0;
    targCom.trajTrans.y = -20.0;
    targCom.trajTrans.z = 20.0;

    reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, 0.0);

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, -10.0) << "X component should reflect to -10";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 10.0) << "Y component should reflect to +10";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, -10.0) << "Z component should reflect to -10";

    const double fx = targCom.comCoord.x + targCom.trajTrans.x;
    const double fy = targCom.comCoord.y + targCom.trajTrans.y;
    const double fz = targCom.comCoord.z + targCom.trajTrans.z;
    std::cerr << "  Final position = (" << fx << ", " << fy << ", " << fz << "), box half-width " << half << '\n';

    EXPECT_LE(fx, half + 1e-9) << "final X must be inside +X wall";
    EXPECT_GE(fy, -half - 1e-9) << "final Y must be inside -Y wall";
    EXPECT_LE(fz, half + 1e-9) << "final Z must be inside +Z wall";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario runs in its own TEST so that a failure in
// one does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(ReflectTrajCheckSpanBox, InteriorMoveUntouched) { test_rtcsb_interior_move_untouched(); }
TEST(ReflectTrajCheckSpanBox, ReflectPositiveXWall) { test_rtcsb_reflect_positive_x_wall(); }
TEST(ReflectTrajCheckSpanBox, ReflectNegativeYWall) { test_rtcsb_reflect_negative_y_wall(); }
TEST(ReflectTrajCheckSpanBox, InterfaceDefinesExtremePoint) { test_rtcsb_interface_defines_extreme_point(); }
TEST(ReflectTrajCheckSpanBox, MultiMemberUsesFarthestMember) { test_rtcsb_multi_member_uses_farthest_member(); }
TEST(ReflectTrajCheckSpanBox, RS3DRaisesNegativeZWall) { test_rtcsb_rs3d_raises_negative_z_wall(); }
TEST(ReflectTrajCheckSpanBox, OnSurfaceIgnoresRS3D) { test_rtcsb_onsurface_ignores_rs3d(); }
TEST(ReflectTrajCheckSpanBox, SpanningComplexResamplesAndTerminates)
{
    test_rtcsb_spanning_complex_resamples_and_terminates();
}
TEST(ReflectTrajCheckSpanBox, ReflectsAllThreeDimensions) { test_rtcsb_reflects_all_three_dimensions(); }