/*! \file test_reflect_complex_compartment.cpp
 *
 * ### Unit test for src/boundary_conditions/reflect_complex_compartment.cpp
 *
 * Function under test:
 *
 *     void reflect_complex_compartment(const Membrane& membraneObject,
 *                                      Complex&              targCom,
 *                                      std::vector<Molecule>& moleculeList,
 *                                      double                 RS3Dinput)
 *
 * ### What the function does (and therefore what we assert on)
 *
 * The routine pushes a Complex *out* of a spherical compartment that is
 * centred on the origin.  The effective compartment radius used is
 *
 *     sphereR = membraneObject.compartmentR                (if targCom.D.z < 1e-8)
 *     sphereR = membraneObject.compartmentR + RS3Dinput    (otherwise)
 *
 * Steps performed by the routine:
 *   1. `canBeInsideX` is true only when `|comCoord| + targCom.radius < sphereR`,
 *      i.e. the whole bounding sphere of the complex sits inside the
 *      compartment.  If this test fails the routine returns without touching
 *      anything.
 *   2. If it is true, every member-molecule centre of mass and every one of
 *      their interfaces is scanned and the point with the *smallest* distance
 *      to the origin (that is still smaller than sphereR) is remembered.
 *      Call that distance `rtmp` and that point `targcrds`.
 *   3. A single rigid translation
 *          dtrans = (-2 * (rtmp - sphereR) / rtmp) * targcrds
 *      is applied to the Complex COM, all member COMs and all interfaces.
 *      The effect is that the selected closest point ends up at radius
 *          2 * sphereR - rtmp
 *      along its original direction (i.e. it is mirrored through the
 *      compartment surface), and the whole complex moves rigidly with it.
 *   4. If no interior point is found (e.g. an empty member list) the
 *      translation degenerates to exactly zero and nothing moves.
 *
 * Every test below prints what it is doing and what the pass criteria are.
 */

#include "boundary_conditions/reflect_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (file-local, uniquely prefixed to avoid collisions in the suite)
// -----------------------------------------------------------------------------
namespace {

//! Tolerance used for all floating point comparisons in this file.
constexpr double kRccTol = 1e-9;

/*! \brief Build a Molecule with a given COM and an arbitrary list of interfaces.
 *
 * \param[in] com          Centre of mass coordinate of the molecule.
 * \param[in] ifaceCoords  Absolute coordinates for each interface.
 * \param[in] molIndex     Index this molecule will have inside moleculeList.
 * \return A Molecule ready to be pushed into a moleculeList.
 */
Molecule rcc_make_molecule(const Coord& com, const std::vector<Coord>& ifaceCoords, int molIndex)
{
    Molecule mol;
    mol.comCoord = com;
    mol.index = molIndex;
    mol.myComIndex = 0; // every molecule created here belongs to complex 0

    mol.interfaceList.clear();
    for (const auto& ic : ifaceCoords) {
        Molecule::Iface iface(ic);
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Build a Complex with a given COM, bounding radius, D.z and members.
 *
 * D.z is the switch used by the function under test to decide whether the
 * complex is membrane bound (D.z ~ 0, RS3Dinput ignored) or free in 3D
 * (D.z > 0, RS3Dinput added to the compartment radius).
 */
Complex rcc_make_complex(const Coord& com, double radius, double Dz, const std::vector<int>& members)
{
    Complex targCom;
    targCom.index = 0;
    targCom.comCoord = com;
    targCom.radius = radius;

    // Only D.z matters for the branch selection, but fill all of them in.
    targCom.D.x = (Dz < 1e-8) ? 1.0 : 1.0;
    targCom.D.y = (Dz < 1e-8) ? 1.0 : 1.0;
    targCom.D.z = Dz;

    targCom.memberList = members;
    return targCom;
}

/*! \brief Convenience: radial distance of a Coord from the origin. */
double rcc_radius(const Coord& c) { return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z); }

/*! \brief Build a Membrane describing a spherical compartment of radius R. */
Membrane rcc_make_membrane(double compartmentR)
{
    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = compartmentR;
    return membraneObject;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: complex whose bounding sphere is NOT fully inside the compartment.
//         `canBeInsideX` is false -> absolutely nothing must change.
// -----------------------------------------------------------------------------
void test_rcc_outside_compartment_no_change()
{
    std::cerr << "\n[TEST] test_rcc_outside_compartment_no_change\n"
              << "  Source file:   reflect_complex_compartment.cpp\n"
              << "  Function:      reflect_complex_compartment\n"
              << "  Scenario:      complex COM well outside the compartment\n"
              << "                 (|COM| + radius >= compartmentR).\n"
              << "  Pass criteria: no coordinate (complex COM, molecule COM or\n"
              << "                 interface) is modified.\n";

    const double compartmentR = 50.0;
    Membrane membraneObject = rcc_make_membrane(compartmentR);

    // Complex sits at radius 100, far outside the 50 nm compartment.
    std::vector<Molecule> moleculeList { rcc_make_molecule(Coord(100.0, 0.0, 0.0),
        { Coord(103.0, 0.0, 0.0) }, 0) };
    Complex targCom = rcc_make_complex(Coord(100.0, 0.0, 0.0), 5.0, /*Dz=*/1.0, { 0 });

    std::cerr << "  Before: complex COM radius = " << rcc_radius(targCom.comCoord)
              << ", compartmentR = " << compartmentR << '\n';

    reflect_complex_compartment(membraneObject, targCom, moleculeList, /*RS3Dinput=*/0.0);

    std::cerr << "  After:  complex COM radius = " << rcc_radius(targCom.comCoord) << '\n';

    // Complex COM must be untouched.
    EXPECT_NEAR(targCom.comCoord.x, 100.0, kRccTol) << "complex COM x must not move";
    EXPECT_NEAR(targCom.comCoord.y, 0.0, kRccTol) << "complex COM y must not move";
    EXPECT_NEAR(targCom.comCoord.z, 0.0, kRccTol) << "complex COM z must not move";

    // Molecule COM and its interface must be untouched.
    EXPECT_NEAR(moleculeList[0].comCoord.x, 100.0, kRccTol) << "molecule COM x must not move";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.x, 103.0, kRccTol)
        << "interface x must not move";
}

// -----------------------------------------------------------------------------
// Test 2: complex fully inside the compartment is mirrored through the surface.
//         Verifies the exact displacement formula and rigid-body behaviour.
// -----------------------------------------------------------------------------
void test_rcc_inside_is_pushed_outward()
{
    std::cerr << "\n[TEST] test_rcc_inside_is_pushed_outward\n"
              << "  Source file:   reflect_complex_compartment.cpp\n"
              << "  Function:      reflect_complex_compartment\n"
              << "  Scenario:      one molecule (COM + 1 interface) entirely inside a\n"
              << "                 compartment of radius 70 (D.z > 0, RS3Dinput = 20\n"
              << "                 -> sphereR = 50 + 20 = 70).\n"
              << "  Pass criteria: the closest point (radius 10) is mirrored to\n"
              << "                 2*sphereR - 10 = 130, and every other coordinate\n"
              << "                 is shifted by the SAME translation vector.\n";

    const double compartmentR = 50.0;
    const double RS3Dinput = 20.0;
    const double sphereR = compartmentR + RS3Dinput; // D.z > 0 -> offset applies
    Membrane membraneObject = rcc_make_membrane(compartmentR);

    // Closest point to the origin is the molecule COM at radius 10.
    std::vector<Molecule> moleculeList { rcc_make_molecule(Coord(10.0, 0.0, 0.0),
        { Coord(12.0, 0.0, 0.0) }, 0) };
    Complex targCom = rcc_make_complex(Coord(11.0, 0.0, 0.0), /*radius=*/2.0, /*Dz=*/1.0, { 0 });

    // Expected translation, computed independently of the implementation.
    const double rtmp = 10.0;
    const double lamda = -2.0 * (rtmp - sphereR) / rtmp; // = 2*(70-10)/10 = 12
    const Coord expectedDtrans(lamda * 10.0, 0.0, 0.0); // = (120, 0, 0)

    std::cerr << "  Expected lambda = " << lamda << ", expected dtrans.x = "
              << expectedDtrans.x << '\n';

    reflect_complex_compartment(membraneObject, targCom, moleculeList, RS3Dinput);

    std::cerr << "  After: molecule COM radius = " << rcc_radius(moleculeList[0].comCoord)
              << " (expected " << 2.0 * sphereR - rtmp << ")\n";

    // The mirrored closest point must sit at 2*sphereR - rtmp.
    EXPECT_NEAR(rcc_radius(moleculeList[0].comCoord), 2.0 * sphereR - rtmp, 1e-8)
        << "closest point should be mirrored through the compartment surface";

    // Rigid translation: every coordinate moved by exactly the same vector.
    EXPECT_NEAR(targCom.comCoord.x, 11.0 + expectedDtrans.x, 1e-8)
        << "complex COM must be shifted by dtrans";
    EXPECT_NEAR(moleculeList[0].comCoord.x, 10.0 + expectedDtrans.x, 1e-8)
        << "molecule COM must be shifted by dtrans";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.x, 12.0 + expectedDtrans.x, 1e-8)
        << "interface must be shifted by dtrans";

    // Perpendicular components must be untouched (translation was purely +x).
    EXPECT_NEAR(targCom.comCoord.y, 0.0, kRccTol) << "no y motion expected";
    EXPECT_NEAR(targCom.comCoord.z, 0.0, kRccTol) << "no z motion expected";
}

// -----------------------------------------------------------------------------
// Test 3: the point used for the reflection is the one CLOSEST to the origin,
//         no matter whether it is a molecule COM or an interface.
// -----------------------------------------------------------------------------
void test_rcc_selects_closest_point()
{
    std::cerr << "\n[TEST] test_rcc_selects_closest_point\n"
              << "  Source file:   reflect_complex_compartment.cpp\n"
              << "  Function:      reflect_complex_compartment\n"
              << "  Scenario:      two molecules with several interfaces; the minimum\n"
              << "                 radius (5, an interface of molecule 1) determines\n"
              << "                 the reflection.  compartmentR = 100, RS3Dinput = 0.\n"
              << "  Pass criteria: dtrans = 2*(100-5)/5 * (0,5,0) = (0,190,0) is\n"
              << "                 applied to every coordinate and the closest point\n"
              << "                 ends up at radius 2*100-5 = 195.\n";

    const double compartmentR = 100.0;
    Membrane membraneObject = rcc_make_membrane(compartmentR);

    // Radii present: 20 (mol0 COM), 30 (mol0 iface), 10 (mol1 COM), 5 (mol1 iface).
    std::vector<Molecule> moleculeList {
        rcc_make_molecule(Coord(20.0, 0.0, 0.0), { Coord(30.0, 0.0, 0.0) }, 0),
        rcc_make_molecule(Coord(0.0, 10.0, 0.0), { Coord(0.0, 5.0, 0.0) }, 1)
    };
    Complex targCom = rcc_make_complex(Coord(10.0, 0.0, 0.0), /*radius=*/5.0, /*Dz=*/1.0, { 0, 1 });

    const double rtmp = 5.0;
    const double lamda = -2.0 * (rtmp - compartmentR) / rtmp; // = 38
    const double dy = lamda * 5.0; // = 190

    std::cerr << "  Expected dtrans = (0, " << dy << ", 0)\n";

    reflect_complex_compartment(membraneObject, targCom, moleculeList, /*RS3Dinput=*/0.0);

    std::cerr << "  After: closest-point radius = "
              << rcc_radius(moleculeList[1].interfaceList[0].coord)
              << " (expected " << 2.0 * compartmentR - rtmp << ")\n";

    // The selected (minimum radius) point is mirrored.
    EXPECT_NEAR(moleculeList[1].interfaceList[0].coord.y, 5.0 + dy, 1e-8)
        << "the minimum-radius interface should be mirrored to 2*sphereR - rtmp";
    EXPECT_NEAR(rcc_radius(moleculeList[1].interfaceList[0].coord), 2.0 * compartmentR - rtmp, 1e-8)
        << "mirrored radius must be 2*sphereR - rtmp";

    // Everything else translated by the same (0, dy, 0) vector.
    EXPECT_NEAR(targCom.comCoord.x, 10.0, 1e-8) << "complex COM x unchanged (pure +y shift)";
    EXPECT_NEAR(targCom.comCoord.y, 0.0 + dy, 1e-8) << "complex COM y shifted by dtrans";
    EXPECT_NEAR(moleculeList[0].comCoord.x, 20.0, 1e-8) << "mol0 COM x unchanged";
    EXPECT_NEAR(moleculeList[0].comCoord.y, 0.0 + dy, 1e-8) << "mol0 COM y shifted by dtrans";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.y, 0.0 + dy, 1e-8)
        << "mol0 interface y shifted by dtrans";
    EXPECT_NEAR(moleculeList[1].comCoord.y, 10.0 + dy, 1e-8) << "mol1 COM y shifted by dtrans";
}

// -----------------------------------------------------------------------------
// Test 4: membrane-bound complex (D.z < 1e-8) must IGNORE RS3Dinput.
// -----------------------------------------------------------------------------
void test_rcc_membrane_bound_ignores_rs3d()
{
    std::cerr << "\n[TEST] test_rcc_membrane_bound_ignores_rs3d\n"
              << "  Source file:   reflect_complex_compartment.cpp\n"
              << "  Function:      reflect_complex_compartment\n"
              << "  Scenario:      D.z = 0 (membrane bound).  The same geometry is\n"
              << "                 reflected twice, once with RS3Dinput = 0 and once\n"
              << "                 with RS3Dinput = 20.\n"
              << "  Pass criteria: both runs give the identical result computed from\n"
              << "                 sphereR = compartmentR (= 50), i.e. RS3Dinput is\n"
              << "                 ignored for membrane-bound complexes.\n";

    const double compartmentR = 50.0;
    Membrane membraneObject = rcc_make_membrane(compartmentR);

    // Expected displacement using sphereR == compartmentR and rtmp == 10.
    const double rtmp = 10.0;
    const double lamdaExpected = -2.0 * (rtmp - compartmentR) / rtmp; // = 8
    const double expectedX = 10.0 + lamdaExpected * 10.0; // 10 + 80 = 90

    // --- run A: RS3Dinput = 0 -------------------------------------------------
    std::vector<Molecule> molsA { rcc_make_molecule(Coord(10.0, 0.0, 0.0), {}, 0) };
    Complex comA = rcc_make_complex(Coord(10.0, 0.0, 0.0), /*radius=*/1.0, /*Dz=*/0.0, { 0 });
    reflect_complex_compartment(membraneObject, comA, molsA, /*RS3Dinput=*/0.0);

    // --- run B: RS3Dinput = 20 (must be ignored) ------------------------------
    std::vector<Molecule> molsB { rcc_make_molecule(Coord(10.0, 0.0, 0.0), {}, 0) };
    Complex comB = rcc_make_complex(Coord(10.0, 0.0, 0.0), /*radius=*/1.0, /*Dz=*/0.0, { 0 });
    reflect_complex_compartment(membraneObject, comB, molsB, /*RS3Dinput=*/20.0);

    std::cerr << "  run A (RS3D=0)  new molecule COM x = " << molsA[0].comCoord.x << '\n'
              << "  run B (RS3D=20) new molecule COM x = " << molsB[0].comCoord.x << '\n'
              << "  expected (sphereR = compartmentR) x = " << expectedX << '\n';

    EXPECT_NEAR(molsA[0].comCoord.x, expectedX, 1e-8)
        << "membrane-bound reflection must use compartmentR only";
    EXPECT_NEAR(molsB[0].comCoord.x, expectedX, 1e-8)
        << "RS3Dinput must be ignored when D.z < 1e-8";
    EXPECT_NEAR(molsA[0].comCoord.x, molsB[0].comCoord.x, 1e-8)
        << "both RS3Dinput values must give identical results for D.z = 0";
}

// -----------------------------------------------------------------------------
// Test 5: 3D complex (D.z > 0) must ADD RS3Dinput to the compartment radius.
// -----------------------------------------------------------------------------
void test_rcc_three_d_uses_rs3d_offset()
{
    std::cerr << "\n[TEST] test_rcc_three_d_uses_rs3d_offset\n"
              << "  Source file:   reflect_complex_compartment.cpp\n"
              << "  Function:      reflect_complex_compartment\n"
              << "  Scenario A:    D.z = 1, compartmentR = 50, RS3Dinput = 20 so\n"
              << "                 sphereR = 70; complex at radius 55 with bounding\n"
              << "                 radius 2 => 57 < 70 -> reflection happens.\n"
              << "  Scenario B:    same geometry but RS3Dinput = 0 so sphereR = 50 and\n"
              << "                 55 + 2 > 50 -> nothing happens.\n"
              << "  Pass criteria: A moves the molecule to 2*70-55 = 85, B leaves it\n"
              << "                 exactly where it was.\n";

    const double compartmentR = 50.0;
    Membrane membraneObject = rcc_make_membrane(compartmentR);

    // --- scenario A: offset makes the complex "inside" ------------------------
    std::vector<Molecule> molsA { rcc_make_molecule(Coord(55.0, 0.0, 0.0), {}, 0) };
    Complex comA = rcc_make_complex(Coord(55.0, 0.0, 0.0), /*radius=*/2.0, /*Dz=*/1.0, { 0 });
    reflect_complex_compartment(membraneObject, comA, molsA, /*RS3Dinput=*/20.0);

    const double sphereRA = compartmentR + 20.0; // 70
    const double expectedA = 2.0 * sphereRA - 55.0; // 85

    std::cerr << "  scenario A new molecule COM x = " << molsA[0].comCoord.x
              << " (expected " << expectedA << ")\n";
    EXPECT_NEAR(molsA[0].comCoord.x, expectedA, 1e-8)
        << "with RS3Dinput = 20 the point should mirror through r = 70";

    // --- scenario B: no offset -> complex is not inside, nothing changes ------
    std::vector<Molecule> molsB { rcc_make_molecule(Coord(55.0, 0.0, 0.0), {}, 0) };
    Complex comB = rcc_make_complex(Coord(55.0, 0.0, 0.0), /*radius=*/2.0, /*Dz=*/1.0, { 0 });
    reflect_complex_compartment(membraneObject, comB, molsB, /*RS3Dinput=*/0.0);

    std::cerr << "  scenario B new molecule COM x = " << molsB[0].comCoord.x
              << " (expected 55, unchanged)\n";
    EXPECT_NEAR(molsB[0].comCoord.x, 55.0, kRccTol)
        << "without the RS3D offset the complex is outside and must not move";
    EXPECT_NEAR(comB.comCoord.x, 55.0, kRccTol)
        << "complex COM must not move when canBeInsideX is false";
}

// -----------------------------------------------------------------------------
// Test 6: exactly on the boundary (|COM| + radius == sphereR) -> strict '<'
//         comparison fails, so nothing must move.
// -----------------------------------------------------------------------------
void test_rcc_exact_boundary_no_change()
{
    std::cerr << "\n[TEST] test_rcc_exact_boundary_no_change\n"
              << "  Source file:   reflect_complex_compartment.cpp\n"
              << "  Function:      reflect_complex_compartment\n"
              << "  Scenario:      |COM| + radius == sphereR exactly (40 + 10 == 50).\n"
              << "  Pass criteria: the strict '<' test fails, so no coordinate moves.\n";

    const double compartmentR = 50.0;
    Membrane membraneObject = rcc_make_membrane(compartmentR);

    std::vector<Molecule> moleculeList { rcc_make_molecule(Coord(40.0, 0.0, 0.0),
        { Coord(45.0, 0.0, 0.0) }, 0) };
    Complex targCom = rcc_make_complex(Coord(40.0, 0.0, 0.0), /*radius=*/10.0, /*Dz=*/1.0, { 0 });

    reflect_complex_compartment(membraneObject, targCom, moleculeList, /*RS3Dinput=*/0.0);

    std::cerr << "  After: complex COM x = " << targCom.comCoord.x
              << ", molecule COM x = " << moleculeList[0].comCoord.x << " (both expected 40)\n";

    EXPECT_NEAR(targCom.comCoord.x, 40.0, kRccTol) << "complex COM must stay at 40";
    EXPECT_NEAR(moleculeList[0].comCoord.x, 40.0, kRccTol) << "molecule COM must stay at 40";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.x, 45.0, kRccTol)
        << "interface must stay at 45";
}

// -----------------------------------------------------------------------------
// Test 7: complex flagged as "inside" but with no member molecules.  No interior
//         point can be found, so the computed translation degenerates to zero.
// -----------------------------------------------------------------------------
void test_rcc_empty_member_list_no_change()
{
    std::cerr << "\n[TEST] test_rcc_empty_member_list_no_change\n"
              << "  Source file:   reflect_complex_compartment.cpp\n"
              << "  Function:      reflect_complex_compartment\n"
              << "  Scenario:      canBeInsideX is true but memberList is empty, so no\n"
              << "                 interior point can be found and the fallback point\n"
              << "                 sits exactly on the surface -> lambda = 0.\n"
              << "  Pass criteria: the complex COM is unchanged (zero translation) and\n"
              << "                 the (empty) molecule list is untouched.\n";

    const double compartmentR = 50.0;
    Membrane membraneObject = rcc_make_membrane(compartmentR);

    std::vector<Molecule> moleculeList; // deliberately empty
    Complex targCom = rcc_make_complex(Coord(5.0, 6.0, 7.0), /*radius=*/1.0, /*Dz=*/1.0, {});

    reflect_complex_compartment(membraneObject, targCom, moleculeList, /*RS3Dinput=*/0.0);

    std::cerr << "  After: complex COM = (" << targCom.comCoord.x << ", "
              << targCom.comCoord.y << ", " << targCom.comCoord.z
              << ") expected (5, 6, 7)\n";

    EXPECT_NEAR(targCom.comCoord.x, 5.0, kRccTol) << "COM x must not move (lambda == 0)";
    EXPECT_NEAR(targCom.comCoord.y, 6.0, kRccTol) << "COM y must not move (lambda == 0)";
    EXPECT_NEAR(targCom.comCoord.z, 7.0, kRccTol) << "COM z must not move (lambda == 0)";
    EXPECT_TRUE(moleculeList.empty()) << "molecule list must remain empty";
}

// -----------------------------------------------------------------------------
// Test 8: the reflection preserves internal geometry (all pairwise separations
//         inside the complex are conserved by a rigid translation).
// -----------------------------------------------------------------------------
void test_rcc_preserves_internal_geometry()
{
    std::cerr << "\n[TEST] test_rcc_preserves_internal_geometry\n"
              << "  Source file:   reflect_complex_compartment.cpp\n"
              << "  Function:      reflect_complex_compartment\n"
              << "  Scenario:      two molecules with interfaces are reflected out of a\n"
              << "                 compartment of radius 100 (RS3Dinput = 0, D.z > 0).\n"
              << "  Pass criteria: the COM-COM distance and every COM-interface\n"
              << "                 distance are unchanged (pure rigid translation).\n";

    const double compartmentR = 100.0;
    Membrane membraneObject = rcc_make_membrane(compartmentR);

    std::vector<Molecule> moleculeList {
        rcc_make_molecule(Coord(12.0, 3.0, -4.0), { Coord(14.0, 3.0, -4.0), Coord(12.0, 5.0, -4.0) }, 0),
        rcc_make_molecule(Coord(-8.0, 9.0, 2.0), { Coord(-8.0, 9.0, 5.0) }, 1)
    };
    Complex targCom = rcc_make_complex(Coord(2.0, 6.0, -1.0), /*radius=*/20.0, /*Dz=*/2.0, { 0, 1 });

    // Record the internal distances before the call.
    auto dist = [](const Coord& a, const Coord& b) {
        return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)
            + (a.z - b.z) * (a.z - b.z));
    };
    const double d_com_com = dist(moleculeList[0].comCoord, moleculeList[1].comCoord);
    const double d_m0_i0 = dist(moleculeList[0].comCoord, moleculeList[0].interfaceList[0].coord);
    const double d_m0_i1 = dist(moleculeList[0].comCoord, moleculeList[0].interfaceList[1].coord);
    const double d_m1_i0 = dist(moleculeList[1].comCoord, moleculeList[1].interfaceList[0].coord);

    std::cerr << "  Before: d(mol0,mol1) = " << d_com_com
              << ", d(mol0,iface0) = " << d_m0_i0 << '\n';

    reflect_complex_compartment(membraneObject, targCom, moleculeList, /*RS3Dinput=*/0.0);

    const double d_com_com_after = dist(moleculeList[0].comCoord, moleculeList[1].comCoord);
    const double d_m0_i0_after = dist(moleculeList[0].comCoord, moleculeList[0].interfaceList[0].coord);
    const double d_m0_i1_after = dist(moleculeList[0].comCoord, moleculeList[0].interfaceList[1].coord);
    const double d_m1_i0_after = dist(moleculeList[1].comCoord, moleculeList[1].interfaceList[0].coord);

    std::cerr << "  After:  d(mol0,mol1) = " << d_com_com_after
              << ", d(mol0,iface0) = " << d_m0_i0_after << '\n';

    EXPECT_NEAR(d_com_com_after, d_com_com, 1e-8) << "COM-COM distance must be conserved";
    EXPECT_NEAR(d_m0_i0_after, d_m0_i0, 1e-8) << "mol0-iface0 distance must be conserved";
    EXPECT_NEAR(d_m0_i1_after, d_m0_i1, 1e-8) << "mol0-iface1 distance must be conserved";
    EXPECT_NEAR(d_m1_i0_after, d_m1_i0, 1e-8) << "mol1-iface0 distance must be conserved";

    // Sanity check: the complex actually moved (it started well inside).
    const bool moved = std::abs(targCom.comCoord.x - 2.0) > 1e-8
        || std::abs(targCom.comCoord.y - 6.0) > 1e-8
        || std::abs(targCom.comCoord.z + 1.0) > 1e-8;
    EXPECT_TRUE(moved) << "an interior complex should have been displaced by the reflection";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named helper is executed inside its own TEST so that
// a failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(ReflectComplexCompartment, OutsideCompartmentNoChange) { test_rcc_outside_compartment_no_change(); }
TEST(ReflectComplexCompartment, InsideIsPushedOutward) { test_rcc_inside_is_pushed_outward(); }
TEST(ReflectComplexCompartment, SelectsClosestPoint) { test_rcc_selects_closest_point(); }
TEST(ReflectComplexCompartment, MembraneBoundIgnoresRS3D) { test_rcc_membrane_bound_ignores_rs3d(); }
TEST(ReflectComplexCompartment, ThreeDUsesRS3DOffset) { test_rcc_three_d_uses_rs3d_offset(); }
TEST(ReflectComplexCompartment, ExactBoundaryNoChange) { test_rcc_exact_boundary_no_change(); }
TEST(ReflectComplexCompartment, EmptyMemberListNoChange) { test_rcc_empty_member_list_no_change(); }
TEST(ReflectComplexCompartment, PreservesInternalGeometry) { test_rcc_preserves_internal_geometry(); }