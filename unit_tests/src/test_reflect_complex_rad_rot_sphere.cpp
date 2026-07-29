/*! \file test_reflect_complex_rad_rot_sphere.cpp
 *
 * ### Unit test for src/boundary_conditions/reflect_complex_rad_rot_sphere.cpp
 *
 * Function under test:
 *
 *     void reflect_complex_rad_rot_sphere(const Membrane& membraneObject,
 *                                         Complex& targCom,
 *                                         std::vector<Molecule>& moleculeList,
 *                                         double radius,
 *                                         double RS3Dinput)
 *
 * What the routine does (and therefore what we verify here):
 *
 *   1. It selects the effective spherical boundary:
 *        - targCom.OnSurface == true  -> sphereR = membraneObject.sphereR
 *        - targCom.OnSurface == false -> sphereR = radius - RS3Dinput
 *   2. A cheap early-out test is performed: only if
 *      (|targCom.comCoord| + targCom.radius) > sphereR can any part of the
 *      complex possibly stick out of the sphere.
 *   3. If a member molecule COM or interface really is outside the sphere, the
 *      whole complex (COM, member COMs and every interface) is rigidly
 *      translated inward along the radial direction of the *farthest* offending
 *      point, by twice the amount it overshot ("mirror" reflection).  The check
 *      is repeated until nothing sticks out any more.
 *   4. Finally, if the complex contains a lipid whose COM is not sitting on the
 *      membrane sphere (|com| != membraneObject.sphereR, tolerance 1e-4), all
 *      *member molecules* are translated so that lipid lands exactly on the
 *      membrane.  NOTE: in this last step the Complex::comCoord itself is
 *      intentionally (in the source as written) NOT updated - the test records
 *      that documented behaviour.
 *
 * Every test prints what file/function is exercised, what scenario is built and
 * what the pass criterion is, so the console log is self-describing.
 */

#include "boundary_conditions/reflect_functions.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Small local helpers (prefixed rcrrs_ = Reflect Complex Rad Rot Sphere) so that
// they cannot collide with helpers from any other test translation unit.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal Molecule belonging to complex index 0.
 *
 * \param[in] com      Center-of-mass coordinate of the molecule.
 * \param[in] ifaces   Absolute coordinates for each interface.
 * \param[in] isLipid  Flag marking the molecule as a (membrane bound) lipid.
 */
Molecule rcrrs_make_molecule(const Coord& com, const std::vector<Coord>& ifaces, bool isLipid = false)
{
    Molecule mol;
    mol.comCoord = com;
    mol.isLipid = isLipid;
    mol.myComIndex = 0;

    mol.interfaceList.clear();
    for (const auto& oneCrd : ifaces) {
        Molecule::Iface iface;
        iface.coord = oneCrd;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Build a Complex owning the supplied member molecule indices. */
Complex rcrrs_make_complex(
    const Coord& com, double complexRadius, const std::vector<int>& members, bool onSurface)
{
    Complex targCom;
    targCom.index = 0;
    targCom.comCoord = com;
    targCom.radius = complexRadius;
    targCom.memberList = members;
    targCom.OnSurface = onSurface;

    // No propagation is involved in this routine, but keep the fields sane.
    targCom.D = Coord(1.0, 1.0, 1.0);
    targCom.Dr = Coord(0.01, 0.01, 0.01);
    return targCom;
}

/*! \brief Radial distance of a Coord from the sphere center (non-mutating). */
double rcrrs_mag(const Coord& c) { return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z); }

/*! \brief Largest radial distance over all member COMs and all interfaces. */
double rcrrs_max_radius(const Complex& targCom, const std::vector<Molecule>& moleculeList)
{
    double maxR = 0.0;
    for (auto memMol : targCom.memberList) {
        maxR = std::max(maxR, rcrrs_mag(moleculeList[memMol].comCoord));
        for (const auto& iface : moleculeList[memMol].interfaceList)
            maxR = std::max(maxR, rcrrs_mag(iface.coord));
    }
    return maxR;
}

/*! \brief Pretty-print a Coord to stderr. */
void rcrrs_print_coord(const char* label, const Coord& c)
{
    std::cerr << "    " << label << " = (" << c.x << ", " << c.y << ", " << c.z
              << ")  |r| = " << rcrrs_mag(c) << "\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a complex comfortably inside the sphere must not be touched at all.
//         The early-out test (|com| + radius > sphereR) is false here, so the
//         function should return without modifying anything.
// -----------------------------------------------------------------------------
void test_rcrrs_inside_leaves_everything_unchanged()
{
    std::cerr << "\n[TEST] test_rcrrs_inside_leaves_everything_unchanged\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere\n"
              << "  Scenario:      complex near the sphere center, nothing can\n"
              << "                 possibly stick out (|com| + radius < sphereR).\n"
              << "  Pass criteria: complex COM, molecule COM and interface coords\n"
              << "                 are all bit-for-bit unchanged.\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // OnSurface == true -> the routine uses membraneObject.sphereR (100).
    std::vector<Molecule> moleculeList { rcrrs_make_molecule(Coord(0.0, 0.0, 0.0), { Coord(1.0, 0.0, 0.0) }) };
    Complex targCom = rcrrs_make_complex(Coord(0.0, 0.0, 0.0), 2.0, { 0 }, true);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (radius=100, RS3D=0)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 100.0, 0.0);

    rcrrs_print_coord("complex COM ", targCom.comCoord);
    rcrrs_print_coord("molecule COM", moleculeList[0].comCoord);
    rcrrs_print_coord("interface   ", moleculeList[0].interfaceList[0].coord);

    EXPECT_DOUBLE_EQ(targCom.comCoord.x, 0.0) << "interior complex COM.x must not move";
    EXPECT_DOUBLE_EQ(targCom.comCoord.y, 0.0) << "interior complex COM.y must not move";
    EXPECT_DOUBLE_EQ(targCom.comCoord.z, 0.0) << "interior complex COM.z must not move";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 0.0) << "interior molecule COM must not move";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, 1.0)
        << "interior interface coordinate must not move";
}

// -----------------------------------------------------------------------------
// Test 2: early-out condition is true (bounding sphere crosses the boundary) but
//         no actual point is outside -> still no modification.
// -----------------------------------------------------------------------------
void test_rcrrs_bounding_overlap_but_no_point_outside()
{
    std::cerr << "\n[TEST] test_rcrrs_bounding_overlap_but_no_point_outside\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere\n"
              << "  Scenario:      |com| + complex radius > sphereR (so the cheap\n"
              << "                 test triggers) but every COM/interface is still\n"
              << "                 inside the sphere.\n"
              << "  Pass criteria: nothing is translated (the 'outside' flag stays\n"
              << "                 false and the while-loop never runs).\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 10.0;

    std::vector<Molecule> moleculeList { rcrrs_make_molecule(Coord(9.0, 0.0, 0.0), { Coord(9.5, 0.0, 0.0) }) };
    // radius 5 => 9 + 5 = 14 > 10 triggers the check, but 9 and 9.5 are inside.
    Complex targCom = rcrrs_make_complex(Coord(9.0, 0.0, 0.0), 5.0, { 0 }, true);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (sphereR=10)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 10.0, 0.0);

    rcrrs_print_coord("complex COM ", targCom.comCoord);
    rcrrs_print_coord("molecule COM", moleculeList[0].comCoord);
    rcrrs_print_coord("interface   ", moleculeList[0].interfaceList[0].coord);

    EXPECT_DOUBLE_EQ(targCom.comCoord.x, 9.0) << "no point is outside -> complex COM must stay";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 9.0) << "no point is outside -> molecule COM must stay";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, 9.5)
        << "no point is outside -> interface must stay";
}

// -----------------------------------------------------------------------------
// Test 3: OnSurface == false -> the boundary is (radius - RS3Dinput).  A single
//         molecule with its interface hanging outside is mirrored back inside.
//         With RS3Dinput = 0 and radius = 100 the arithmetic is exact:
//           farthest point = interface at r = 106
//           lamda = -2*(106-100)/106, dtrans = lamda * (106,0,0) = (-12,0,0)
// -----------------------------------------------------------------------------
void test_rcrrs_reflects_using_radius_argument()
{
    std::cerr << "\n[TEST] test_rcrrs_reflects_using_radius_argument\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere\n"
              << "  Scenario:      OnSurface = false, radius = 100, RS3Dinput = 0,\n"
              << "                 so sphereR = 100.  The interface sits at r = 106.\n"
              << "  Pass criteria: whole complex is rigidly shifted by (-12,0,0)\n"
              << "                 (twice the 6 nm overshoot) and everything ends\n"
              << "                 up inside r <= 100.\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    // Deliberately different from the 'radius' argument to prove it is unused here.
    membraneObject.sphereR = 1000.0;

    std::vector<Molecule> moleculeList { rcrrs_make_molecule(Coord(105.0, 0.0, 0.0), { Coord(106.0, 0.0, 0.0) }) };
    Complex targCom = rcrrs_make_complex(Coord(105.0, 0.0, 0.0), 1.0, { 0 }, false);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (radius=100, RS3D=0)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 100.0, 0.0);

    rcrrs_print_coord("complex COM ", targCom.comCoord);
    rcrrs_print_coord("molecule COM", moleculeList[0].comCoord);
    rcrrs_print_coord("interface   ", moleculeList[0].interfaceList[0].coord);

    // Exact expected shift of -12 nm along x.
    EXPECT_NEAR(targCom.comCoord.x, 93.0, 1e-9) << "complex COM should be mirrored to x = 93";
    EXPECT_NEAR(moleculeList[0].comCoord.x, 93.0, 1e-9) << "molecule COM should be mirrored to x = 93";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.x, 94.0, 1e-9)
        << "interface should be mirrored to x = 94 (rigid-body shift)";

    // Transverse components must not be touched (motion is purely radial here).
    EXPECT_NEAR(targCom.comCoord.y, 0.0, 1e-12) << "no y motion expected for a purely radial reflection";
    EXPECT_NEAR(targCom.comCoord.z, 0.0, 1e-12) << "no z motion expected for a purely radial reflection";

    // And the global invariant: nothing sticks out any more.
    const double maxR = rcrrs_max_radius(targCom, moleculeList);
    std::cerr << "    max radial distance after reflection = " << maxR << " (boundary 100)\n";
    EXPECT_LE(maxR, 100.0 + 1e-9) << "after reflection no point may lie outside the sphere";
}

// -----------------------------------------------------------------------------
// Test 4: RS3Dinput shrinks the effective boundary when OnSurface == false.
//         radius = 100, RS3Dinput = 20 -> sphereR = 80.  A molecule at r = 90
//         overshoots by 10 and must be mirrored to r = 70.
// -----------------------------------------------------------------------------
void test_rcrrs_rs3dinput_shrinks_boundary()
{
    std::cerr << "\n[TEST] test_rcrrs_rs3dinput_shrinks_boundary\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere\n"
              << "  Scenario:      OnSurface = false, radius = 100, RS3Dinput = 20,\n"
              << "                 so the reflecting surface is at r = 80 and the\n"
              << "                 molecule sits at r = 90.\n"
              << "  Pass criteria: molecule/complex mirrored to r = 70 and the\n"
              << "                 maximum radius is <= 80.\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0; // unused on this code path (no lipids present)

    std::vector<Molecule> moleculeList { rcrrs_make_molecule(Coord(90.0, 0.0, 0.0), { Coord(90.0, 0.0, 0.0) }) };
    Complex targCom = rcrrs_make_complex(Coord(90.0, 0.0, 0.0), 1.0, { 0 }, false);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (radius=100, RS3D=20)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 100.0, 20.0);

    rcrrs_print_coord("complex COM ", targCom.comCoord);
    rcrrs_print_coord("molecule COM", moleculeList[0].comCoord);

    EXPECT_NEAR(targCom.comCoord.x, 70.0, 1e-9)
        << "with sphereR = radius - RS3Dinput = 80, r=90 must mirror to r=70";
    EXPECT_NEAR(moleculeList[0].comCoord.x, 70.0, 1e-9) << "member molecule must follow the complex";

    const double maxR = rcrrs_max_radius(targCom, moleculeList);
    std::cerr << "    max radial distance after reflection = " << maxR << " (boundary 80)\n";
    EXPECT_LE(maxR, 80.0 + 1e-9) << "no point may lie outside the reduced boundary";
}

// -----------------------------------------------------------------------------
// Test 5: OnSurface == true -> the membrane radius is used and the 'radius'
//         argument is ignored.  We pass a huge 'radius' (1000): if it were used
//         no reflection would happen at all, so a reflection proves the branch.
// -----------------------------------------------------------------------------
void test_rcrrs_onsurface_uses_membrane_radius()
{
    std::cerr << "\n[TEST] test_rcrrs_onsurface_uses_membrane_radius\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere\n"
              << "  Scenario:      OnSurface = true, membraneObject.sphereR = 50 but\n"
              << "                 the 'radius' argument is 1000.  Molecule at r=55.\n"
              << "  Pass criteria: reflection happens against r = 50 (result r = 45),\n"
              << "                 which can only occur if membraneObject.sphereR is\n"
              << "                 the boundary that was used.\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 50.0;

    std::vector<Molecule> moleculeList { rcrrs_make_molecule(Coord(55.0, 0.0, 0.0), { Coord(55.0, 0.0, 0.0) }) };
    Complex targCom = rcrrs_make_complex(Coord(55.0, 0.0, 0.0), 1.0, { 0 }, true);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (radius=1000, RS3D=0)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 1000.0, 0.0);

    rcrrs_print_coord("complex COM ", targCom.comCoord);
    rcrrs_print_coord("molecule COM", moleculeList[0].comCoord);

    EXPECT_NEAR(targCom.comCoord.x, 45.0, 1e-9)
        << "OnSurface complexes must reflect off membraneObject.sphereR (50), not 'radius'";
    EXPECT_NEAR(moleculeList[0].comCoord.x, 45.0, 1e-9) << "member molecule must follow the complex";

    const double maxR = rcrrs_max_radius(targCom, moleculeList);
    std::cerr << "    max radial distance after reflection = " << maxR << " (boundary 50)\n";
    EXPECT_LE(maxR, 50.0 + 1e-9) << "no point may lie outside the membrane sphere";
}

// -----------------------------------------------------------------------------
// Test 6: repeated reflection.  The first mirror move pulls the COM inside but
//         throws an interface out on another side, so a second pass of the
//         while-loop is required.  We only assert the final invariant (nothing
//         outside) plus rigidity of the motion, which is robust to the exact
//         iteration arithmetic.
// -----------------------------------------------------------------------------
void test_rcrrs_multiple_reflections_converge()
{
    std::cerr << "\n[TEST] test_rcrrs_multiple_reflections_converge\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere\n"
              << "  Scenario:      sphereR = 10.  Molecule COM at (10.5,0,0) and an\n"
              << "                 interface at (0,10.4,0): reflecting the COM inward\n"
              << "                 pushes the interface out, forcing a second pass.\n"
              << "  Pass criteria: the loop converges (no exit(1)), the final maximum\n"
              << "                 radius is <= 10, and COM/interface separation is\n"
              << "                 preserved (rigid-body translation only).\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 10.0;

    const Coord molCom0(10.5, 0.0, 0.0);
    const Coord iface0(0.0, 10.4, 0.0);

    std::vector<Molecule> moleculeList { rcrrs_make_molecule(molCom0, { iface0 }) };
    Complex targCom = rcrrs_make_complex(molCom0, 1.0, { 0 }, true);

    // Internal geometry (interface relative to molecule COM) before the call.
    const Coord relBefore(iface0.x - molCom0.x, iface0.y - molCom0.y, iface0.z - molCom0.z);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (sphereR=10)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 10.0, 0.0);

    rcrrs_print_coord("complex COM ", targCom.comCoord);
    rcrrs_print_coord("molecule COM", moleculeList[0].comCoord);
    rcrrs_print_coord("interface   ", moleculeList[0].interfaceList[0].coord);

    // 1) Everything must now be inside the sphere.
    const double maxR = rcrrs_max_radius(targCom, moleculeList);
    std::cerr << "    max radial distance after reflections = " << maxR << " (boundary 10)\n";
    EXPECT_LE(maxR, 10.0 + 1e-9) << "iterated reflection must place every point inside the sphere";

    // 2) The complex actually moved (it started outside).
    EXPECT_TRUE(std::fabs(targCom.comCoord.x - molCom0.x) > 1e-9
        || std::fabs(targCom.comCoord.y - molCom0.y) > 1e-9)
        << "an out-of-sphere complex must be displaced";

    // 3) Rigidity: the interface-to-COM vector is unchanged.
    const Coord relAfter(moleculeList[0].interfaceList[0].coord.x - moleculeList[0].comCoord.x,
        moleculeList[0].interfaceList[0].coord.y - moleculeList[0].comCoord.y,
        moleculeList[0].interfaceList[0].coord.z - moleculeList[0].comCoord.z);
    EXPECT_NEAR(relAfter.x, relBefore.x, 1e-9) << "reflection must be a pure translation (x)";
    EXPECT_NEAR(relAfter.y, relBefore.y, 1e-9) << "reflection must be a pure translation (y)";
    EXPECT_NEAR(relAfter.z, relBefore.z, 1e-9) << "reflection must be a pure translation (z)";

    // 4) Complex COM and member COM received identical shifts.
    EXPECT_NEAR(targCom.comCoord.x, moleculeList[0].comCoord.x, 1e-9)
        << "complex COM and single-member COM must stay coincident";
    EXPECT_NEAR(targCom.comCoord.y, moleculeList[0].comCoord.y, 1e-9)
        << "complex COM and single-member COM must stay coincident";
}

// -----------------------------------------------------------------------------
// Test 7: lipid re-anchoring.  A complex fully inside the sphere (so no
//         reflection occurs) contains a lipid at r = 90 while the membrane is at
//         r = 100.  All member molecules must be shifted by (+10,0,0) so the
//         lipid lands exactly on the membrane.  Documented quirk: the Complex
//         comCoord is NOT updated by this step.
// -----------------------------------------------------------------------------
void test_rcrrs_lipid_snapped_back_onto_membrane()
{
    std::cerr << "\n[TEST] test_rcrrs_lipid_snapped_back_onto_membrane\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere (lipid fix-up block)\n"
              << "  Scenario:      membrane sphereR = 100, complex is inside so no\n"
              << "                 reflection triggers, but a member lipid sits at\n"
              << "                 r = 90 instead of on the membrane.\n"
              << "  Pass criteria: every member molecule (and interface) is shifted\n"
              << "                 by (+10,0,0) so the lipid COM lands at r = 100;\n"
              << "                 Complex::comCoord is left unchanged (as coded).\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    // Index 0: the lipid (off the membrane by 10 nm).  Index 1: a partner protein.
    std::vector<Molecule> moleculeList {
        rcrrs_make_molecule(Coord(90.0, 0.0, 0.0), { Coord(90.0, 0.0, 1.0) }, /*isLipid=*/true),
        rcrrs_make_molecule(Coord(89.0, 0.0, 0.0), { Coord(88.0, 0.0, 0.0) }, /*isLipid=*/false)
    };
    // radius 1 -> 90 + 1 = 91 < 100, therefore the reflection block is skipped.
    Complex targCom = rcrrs_make_complex(Coord(90.0, 0.0, 0.0), 1.0, { 0, 1 }, true);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (sphereR=100)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 100.0, 0.0);

    rcrrs_print_coord("lipid COM   ", moleculeList[0].comCoord);
    rcrrs_print_coord("lipid iface ", moleculeList[0].interfaceList[0].coord);
    rcrrs_print_coord("protein COM ", moleculeList[1].comCoord);
    rcrrs_print_coord("protein ifc ", moleculeList[1].interfaceList[0].coord);
    rcrrs_print_coord("complex COM ", targCom.comCoord);

    // Lipid COM must land exactly on the membrane sphere.
    EXPECT_NEAR(rcrrs_mag(moleculeList[0].comCoord), membraneObject.sphereR, 1e-9)
        << "lipid COM must be placed exactly on the membrane sphere";
    EXPECT_NEAR(moleculeList[0].comCoord.x, 100.0, 1e-9) << "lipid COM should move 90 -> 100 along x";

    // The rest of the complex is translated by the same (+10,0,0) vector.
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.x, 100.0, 1e-9)
        << "lipid interface must be translated with its molecule";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.z, 1.0, 1e-9)
        << "translation is purely along x, so the interface z is unchanged";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 99.0, 1e-9)
        << "partner molecule COM must be translated by the same (+10,0,0)";
    EXPECT_NEAR(moleculeList[1].interfaceList[0].coord.x, 98.0, 1e-9)
        << "partner interface must be translated by the same (+10,0,0)";

    // Documented behaviour of the source: the Complex COM is not touched here.
    EXPECT_NEAR(targCom.comCoord.x, 90.0, 1e-9)
        << "as implemented, the lipid fix-up does not update Complex::comCoord";
}

// -----------------------------------------------------------------------------
// Test 8: a lipid that is already on the membrane (within the 1e-4 tolerance)
//         must not trigger any translation.
// -----------------------------------------------------------------------------
void test_rcrrs_lipid_on_membrane_not_moved()
{
    std::cerr << "\n[TEST] test_rcrrs_lipid_on_membrane_not_moved\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere (lipid fix-up block)\n"
              << "  Scenario:      lipid is off the membrane by only 5e-5 nm, which is\n"
              << "                 below the 1e-4 tolerance used by the routine.\n"
              << "  Pass criteria: no molecule is translated at all.\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    const double lipidR = 100.00005; // within 1e-4 of the membrane radius
    std::vector<Molecule> moleculeList {
        rcrrs_make_molecule(Coord(lipidR, 0.0, 0.0), { Coord(lipidR, 0.0, 0.0) }, /*isLipid=*/true)
    };
    // Complex COM at the origin with a tiny radius -> the reflection early-out is
    // false (0 + 1 < 100), so only the lipid block can act.
    Complex targCom = rcrrs_make_complex(Coord(0.0, 0.0, 0.0), 1.0, { 0 }, true);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (sphereR=100)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 100.0, 0.0);

    rcrrs_print_coord("lipid COM   ", moleculeList[0].comCoord);
    rcrrs_print_coord("complex COM ", targCom.comCoord);

    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, lipidR)
        << "a lipid within 1e-4 of the membrane must be left alone";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, lipidR)
        << "its interface must be left alone as well";
    EXPECT_DOUBLE_EQ(targCom.comCoord.x, 0.0) << "complex COM must be untouched";
}

// -----------------------------------------------------------------------------
// Test 9: a complex outside the sphere in a non axis-aligned direction is
//         mirrored back along its own radial direction.  This checks that the
//         reflection is genuinely radial (direction preserved, magnitude fixed).
// -----------------------------------------------------------------------------
void test_rcrrs_reflection_is_radial_offaxis()
{
    std::cerr << "\n[TEST] test_rcrrs_reflection_is_radial_offaxis\n"
              << "  Source file:   reflect_complex_rad_rot_sphere.cpp\n"
              << "  Function:      reflect_complex_rad_rot_sphere\n"
              << "  Scenario:      single molecule at (6,8,0) => r = 10 with the\n"
              << "                 boundary at r = 8 (OnSurface, membrane sphereR=8).\n"
              << "  Pass criteria: the point is mirrored to r = 2*8 - 10 = 6 while\n"
              << "                 keeping the same radial direction (3/5, 4/5, 0).\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 8.0;

    const Coord start(6.0, 8.0, 0.0); // r = 10
    std::vector<Molecule> moleculeList { rcrrs_make_molecule(start, { start }) };
    Complex targCom = rcrrs_make_complex(start, 0.5, { 0 }, true);

    std::cerr << "  Calling reflect_complex_rad_rot_sphere (sphereR=8)...\n";
    reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, 8.0, 0.0);

    rcrrs_print_coord("molecule COM", moleculeList[0].comCoord);

    // Mirror of r = 10 about the shell at 8 is r = 6, direction unchanged.
    EXPECT_NEAR(rcrrs_mag(moleculeList[0].comCoord), 6.0, 1e-9)
        << "radial distance must become 2*sphereR - r = 6";
    EXPECT_NEAR(moleculeList[0].comCoord.x, 3.6, 1e-9) << "x should be 6 * (6/10) = 3.6";
    EXPECT_NEAR(moleculeList[0].comCoord.y, 4.8, 1e-9) << "y should be 6 * (8/10) = 4.8";
    EXPECT_NEAR(moleculeList[0].comCoord.z, 0.0, 1e-12) << "z should remain zero";

    const double maxR = rcrrs_max_radius(targCom, moleculeList);
    std::cerr << "    max radial distance after reflection = " << maxR << " (boundary 8)\n";
    EXPECT_LE(maxR, 8.0 + 1e-9) << "reflected complex must lie inside the sphere";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario runs in its own TEST so a failure in one
// does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(ReflectComplexRadRotSphere, InsideLeavesEverythingUnchanged)
{
    test_rcrrs_inside_leaves_everything_unchanged();
}
TEST(ReflectComplexRadRotSphere, BoundingOverlapButNoPointOutside)
{
    test_rcrrs_bounding_overlap_but_no_point_outside();
}
TEST(ReflectComplexRadRotSphere, ReflectsUsingRadiusArgument)
{
    test_rcrrs_reflects_using_radius_argument();
}
TEST(ReflectComplexRadRotSphere, RS3DInputShrinksBoundary) { test_rcrrs_rs3dinput_shrinks_boundary(); }
TEST(ReflectComplexRadRotSphere, OnSurfaceUsesMembraneRadius)
{
    test_rcrrs_onsurface_uses_membrane_radius();
}
TEST(ReflectComplexRadRotSphere, MultipleReflectionsConverge)
{
    test_rcrrs_multiple_reflections_converge();
}
TEST(ReflectComplexRadRotSphere, LipidSnappedBackOntoMembrane)
{
    test_rcrrs_lipid_snapped_back_onto_membrane();
}
TEST(ReflectComplexRadRotSphere, LipidOnMembraneNotMoved) { test_rcrrs_lipid_on_membrane_not_moved(); }
TEST(ReflectComplexRadRotSphere, ReflectionIsRadialOffAxis) { test_rcrrs_reflection_is_radial_offaxis(); }