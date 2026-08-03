/*! \file test_transform.cpp
 *
 * ### Unit test for src/reactions/transform.cpp
 *
 * Function under test:
 *
 *     void transform(Coord& reactIface, Molecule& reactMol1,
 *                    Molecule& reactMol2, const Vector& axis)
 *
 * `transform()` is a helper used by the association machinery
 * (calculate_phi() / calculate_omega()).  It rigidly rotates the *temporary*
 * association coordinates (`Molecule::tmpComCoord` and `Molecule::tmpICoords`)
 * of two molecules about the point `reactIface`, so that the supplied `axis`
 * ends up lying on the z-axis.  The rotation quaternion is built as
 *
 *     rotAxis = zhat x axis        (normalised)
 *     theta   = angle(axis, zhat)
 *     Q       = ( cos(-theta/2), sin(-theta/2) * rotAxis )
 *
 * A rotation of -theta about (zhat x axis) maps `axis` onto **+z**, which is the
 * behaviour verified below.
 *
 * Important detail exercised by these tests: `Vector::dot_theta()` uses the
 * *stored* `magnitude` member, therefore every axis handed to `transform()`
 * must have had `calc_magnitude()` called on it first.  The helper
 * `tfm_axis()` below does that.
 *
 * The real coordinates (`comCoord`, `interfaceList`) are *not* touched by
 * `transform()` - only the tmp association coordinates are - and that is also
 * asserted here.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (file-internal so they cannot collide with other test files).
// -----------------------------------------------------------------------------
namespace {

//! Tolerance used for all floating point geometry comparisons.
constexpr double kTfmTol = 1e-9;

/*! \brief Build a Vector and pre-compute its magnitude.
 *
 * `transform()` calls `Vector::dot_theta()` on the axis, which relies on the
 * cached `magnitude` field, so it must be populated before the call.
 */
Vector tfm_axis(double x, double y, double z)
{
    Vector v(x, y, z);
    v.calc_magnitude();
    return v;
}

/*! \brief Build a minimal Molecule holding temporary association coordinates.
 *
 * `transform()` only reads/writes `tmpComCoord` and `tmpICoords`, so those are
 * the fields that matter.  `comCoord` is also filled in so we can prove the
 * real (non-temporary) coordinates are left alone.
 */
Molecule tfm_molecule(const Coord& com, const std::vector<Coord>& ifaces)
{
    Molecule mol;
    mol.comCoord = com; // "real" coordinate - must not be modified
    mol.tmpComCoord = com; // temporary coordinate - will be rotated
    mol.tmpICoords = ifaces; // temporary interface coordinates - rotated
    return mol;
}

//! Euclidean distance between two Coords.
double tfm_distance(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

//! Pretty-print a coordinate to stderr for the verbose test log.
void tfm_print(const char* label, const Coord& c)
{
    std::cerr << "      " << label << " = (" << c.x << ", " << c.y << ", " << c.z << ")\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: axis already on +z  ->  early return, absolutely nothing changes.
// -----------------------------------------------------------------------------
void tfm_test_no_rotation_when_axis_is_plus_z()
{
    std::cerr << "\n[TEST] tfm_test_no_rotation_when_axis_is_plus_z\n"
              << "  Source file:   src/reactions/transform.cpp\n"
              << "  Function:      transform()\n"
              << "  Scenario:      axis == (0,0,1); angle to z-axis is exactly 0,\n"
              << "                 so the function must return before rotating.\n"
              << "  Pass criteria: every tmp coordinate is bitwise unchanged.\n";

    Coord reactIface { 0.0, 0.0, 0.0 };
    Molecule mol1 = tfm_molecule(Coord { 1.0, 2.0, 3.0 }, { Coord { 4.0, 5.0, 6.0 } });
    Molecule mol2 = tfm_molecule(Coord { -1.0, -2.0, -3.0 }, { Coord { -4.0, -5.0, -6.0 } });

    Vector axis = tfm_axis(0.0, 0.0, 1.0);

    std::cerr << "  Calling transform() with axis = (0,0,1)...\n";
    transform(reactIface, mol1, mol2, axis);

    tfm_print("mol1.tmpComCoord", mol1.tmpComCoord);
    tfm_print("mol2.tmpComCoord", mol2.tmpComCoord);

    // Nothing at all should have moved.
    EXPECT_DOUBLE_EQ(mol1.tmpComCoord.x, 1.0) << "mol1 COM x must be untouched";
    EXPECT_DOUBLE_EQ(mol1.tmpComCoord.y, 2.0) << "mol1 COM y must be untouched";
    EXPECT_DOUBLE_EQ(mol1.tmpComCoord.z, 3.0) << "mol1 COM z must be untouched";
    EXPECT_DOUBLE_EQ(mol1.tmpICoords[0].x, 4.0) << "mol1 iface x must be untouched";
    EXPECT_DOUBLE_EQ(mol1.tmpICoords[0].y, 5.0) << "mol1 iface y must be untouched";
    EXPECT_DOUBLE_EQ(mol1.tmpICoords[0].z, 6.0) << "mol1 iface z must be untouched";

    EXPECT_DOUBLE_EQ(mol2.tmpComCoord.x, -1.0) << "mol2 COM x must be untouched";
    EXPECT_DOUBLE_EQ(mol2.tmpComCoord.y, -2.0) << "mol2 COM y must be untouched";
    EXPECT_DOUBLE_EQ(mol2.tmpComCoord.z, -3.0) << "mol2 COM z must be untouched";
    EXPECT_DOUBLE_EQ(mol2.tmpICoords[0].x, -4.0) << "mol2 iface x must be untouched";
    EXPECT_DOUBLE_EQ(mol2.tmpICoords[0].y, -5.0) << "mol2 iface y must be untouched";
    EXPECT_DOUBLE_EQ(mol2.tmpICoords[0].z, -6.0) << "mol2 iface z must be untouched";
}

// -----------------------------------------------------------------------------
// Test 2: axis anti-parallel to z (angle == M_PI) -> early return as well.
// -----------------------------------------------------------------------------
void tfm_test_no_rotation_when_axis_is_minus_z()
{
    std::cerr << "\n[TEST] tfm_test_no_rotation_when_axis_is_minus_z\n"
              << "  Source file:   src/reactions/transform.cpp\n"
              << "  Function:      transform()\n"
              << "  Scenario:      axis == (0,0,-5); angle to z-axis is exactly M_PI,\n"
              << "                 which is the second early-return condition.\n"
              << "  Pass criteria: every tmp coordinate is unchanged.\n";

    Coord reactIface { 2.0, 2.0, 2.0 };
    Molecule mol1 = tfm_molecule(Coord { 3.0, 0.0, 0.0 }, { Coord { 0.0, 3.0, 0.0 } });
    Molecule mol2 = tfm_molecule(Coord { 0.0, 0.0, 3.0 }, {});

    // Non-unit length on purpose: the ratio dot/(|a||b|) is still exactly -1.
    Vector axis = tfm_axis(0.0, 0.0, -5.0);

    std::cerr << "  Calling transform() with axis = (0,0,-5)...\n";
    transform(reactIface, mol1, mol2, axis);

    tfm_print("mol1.tmpComCoord", mol1.tmpComCoord);
    tfm_print("mol1.tmpICoords[0]", mol1.tmpICoords[0]);
    tfm_print("mol2.tmpComCoord", mol2.tmpComCoord);

    EXPECT_DOUBLE_EQ(mol1.tmpComCoord.x, 3.0) << "anti-parallel axis must not rotate mol1";
    EXPECT_DOUBLE_EQ(mol1.tmpComCoord.y, 0.0) << "anti-parallel axis must not rotate mol1";
    EXPECT_DOUBLE_EQ(mol1.tmpComCoord.z, 0.0) << "anti-parallel axis must not rotate mol1";
    EXPECT_DOUBLE_EQ(mol1.tmpICoords[0].y, 3.0) << "anti-parallel axis must not rotate ifaces";
    EXPECT_DOUBLE_EQ(mol2.tmpComCoord.z, 3.0) << "anti-parallel axis must not rotate mol2";
}

// -----------------------------------------------------------------------------
// Test 3: axis == +x, reactIface at the origin.  Fully deterministic rotation.
//
// rotAxis = zhat x xhat = +yhat, theta = pi/2, rotation angle = -pi/2 about y,
// which maps (x,y,z) -> (-z, y, x).
// -----------------------------------------------------------------------------
void tfm_test_x_axis_explicit_rotation()
{
    std::cerr << "\n[TEST] tfm_test_x_axis_explicit_rotation\n"
              << "  Source file:   src/reactions/transform.cpp\n"
              << "  Function:      transform()\n"
              << "  Scenario:      axis == (1,0,0), reactIface at the origin.\n"
              << "                 The quaternion is a -90 deg rotation about +y,\n"
              << "                 i.e. (x,y,z) -> (-z, y, x).\n"
              << "  Pass criteria: every rotated tmp coordinate matches the\n"
              << "                 analytically expected value (tol 1e-9), and\n"
              << "                 the axis direction ends up on +z.\n";

    Coord reactIface { 0.0, 0.0, 0.0 };

    // mol1 sits along the axis (so it should end up on the z axis).
    Molecule mol1 = tfm_molecule(Coord { 1.0, 0.0, 0.0 },
        { Coord { 2.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } });
    // mol2 sits along z (so it should end up along -x).
    Molecule mol2 = tfm_molecule(Coord { 0.0, 0.0, 1.0 },
        { Coord { 0.0, 0.0, 2.0 }, Coord { 0.0, -1.0, 0.0 } });

    Vector axis = tfm_axis(1.0, 0.0, 0.0);

    std::cerr << "  Calling transform() with axis = (1,0,0)...\n";
    transform(reactIface, mol1, mol2, axis);

    tfm_print("mol1.tmpComCoord   (expect  0, 0, 1)", mol1.tmpComCoord);
    tfm_print("mol1.tmpICoords[0] (expect  0, 0, 2)", mol1.tmpICoords[0]);
    tfm_print("mol1.tmpICoords[1] (expect  0, 1, 0)", mol1.tmpICoords[1]);
    tfm_print("mol2.tmpComCoord   (expect -1, 0, 0)", mol2.tmpComCoord);
    tfm_print("mol2.tmpICoords[0] (expect -2, 0, 0)", mol2.tmpICoords[0]);
    tfm_print("mol2.tmpICoords[1] (expect  0,-1, 0)", mol2.tmpICoords[1]);

    // Molecule 1: the point on the axis must land on +z.
    EXPECT_NEAR(mol1.tmpComCoord.x, 0.0, kTfmTol) << "mol1 COM x after -90deg about y";
    EXPECT_NEAR(mol1.tmpComCoord.y, 0.0, kTfmTol) << "mol1 COM y after -90deg about y";
    EXPECT_NEAR(mol1.tmpComCoord.z, 1.0, kTfmTol) << "axis direction must map to +z";

    EXPECT_NEAR(mol1.tmpICoords[0].x, 0.0, kTfmTol) << "mol1 iface0 x";
    EXPECT_NEAR(mol1.tmpICoords[0].y, 0.0, kTfmTol) << "mol1 iface0 y";
    EXPECT_NEAR(mol1.tmpICoords[0].z, 2.0, kTfmTol) << "mol1 iface0 z";

    // A point on the rotation axis (y) is invariant.
    EXPECT_NEAR(mol1.tmpICoords[1].x, 0.0, kTfmTol) << "point on rotation axis: x unchanged";
    EXPECT_NEAR(mol1.tmpICoords[1].y, 1.0, kTfmTol) << "point on rotation axis: y unchanged";
    EXPECT_NEAR(mol1.tmpICoords[1].z, 0.0, kTfmTol) << "point on rotation axis: z unchanged";

    // Molecule 2 is rotated by exactly the same quaternion.
    EXPECT_NEAR(mol2.tmpComCoord.x, -1.0, kTfmTol) << "mol2 COM x after -90deg about y";
    EXPECT_NEAR(mol2.tmpComCoord.y, 0.0, kTfmTol) << "mol2 COM y after -90deg about y";
    EXPECT_NEAR(mol2.tmpComCoord.z, 0.0, kTfmTol) << "mol2 COM z after -90deg about y";

    EXPECT_NEAR(mol2.tmpICoords[0].x, -2.0, kTfmTol) << "mol2 iface0 x";
    EXPECT_NEAR(mol2.tmpICoords[0].y, 0.0, kTfmTol) << "mol2 iface0 y";
    EXPECT_NEAR(mol2.tmpICoords[0].z, 0.0, kTfmTol) << "mol2 iface0 z";

    EXPECT_NEAR(mol2.tmpICoords[1].x, 0.0, kTfmTol) << "mol2 iface1 x (on rotation axis)";
    EXPECT_NEAR(mol2.tmpICoords[1].y, -1.0, kTfmTol) << "mol2 iface1 y (on rotation axis)";
    EXPECT_NEAR(mol2.tmpICoords[1].z, 0.0, kTfmTol) << "mol2 iface1 z (on rotation axis)";
}

// -----------------------------------------------------------------------------
// Test 4: arbitrary axis, non-origin reactIface -> the axis direction must end
//         up on the z-axis (this is the whole purpose of the function).
// -----------------------------------------------------------------------------
void tfm_test_arbitrary_axis_maps_onto_z_axis()
{
    std::cerr << "\n[TEST] tfm_test_arbitrary_axis_maps_onto_z_axis\n"
              << "  Source file:   src/reactions/transform.cpp\n"
              << "  Function:      transform()\n"
              << "  Scenario:      axis == (1,1,1) (non-unit length), reactIface\n"
              << "                 away from the origin, and a molecule placed\n"
              << "                 exactly along the axis from reactIface.\n"
              << "  Pass criteria: after the call, the molecule-to-reactIface\n"
              << "                 vector has x == y == 0 and z == +|vector|,\n"
              << "                 i.e. the axis now lies on +z.\n";

    Coord reactIface { 1.0, -2.0, 0.5 };
    const Coord origIface = reactIface;

    // Displacement of 2 * axis from the interface: this is "on the axis".
    const Coord onAxis { reactIface.x + 2.0, reactIface.y + 2.0, reactIface.z + 2.0 };
    const double axisLen = std::sqrt(3.0) * 2.0; // |(2,2,2)|

    Molecule mol1 = tfm_molecule(onAxis, { onAxis });
    Molecule mol2 = tfm_molecule(Coord { reactIface.x, reactIface.y, reactIface.z + 1.0 }, {});

    Vector axis = tfm_axis(1.0, 1.0, 1.0);

    std::cerr << "  Calling transform() with axis = (1,1,1)...\n";
    transform(reactIface, mol1, mol2, axis);

    // Vector from the (unchanged) rotation origin to the rotated COM.
    const double rx = mol1.tmpComCoord.x - reactIface.x;
    const double ry = mol1.tmpComCoord.y - reactIface.y;
    const double rz = mol1.tmpComCoord.z - reactIface.z;
    std::cerr << "      rotated on-axis vector = (" << rx << ", " << ry << ", " << rz
              << "), expected (0, 0, " << axisLen << ")\n";

    EXPECT_NEAR(rx, 0.0, 1e-8) << "the on-axis vector must have no x component left";
    EXPECT_NEAR(ry, 0.0, 1e-8) << "the on-axis vector must have no y component left";
    EXPECT_NEAR(rz, axisLen, 1e-8) << "the on-axis vector must lie along +z with its length preserved";

    // The same must hold for the interface coordinate that was on the axis.
    const double ix = mol1.tmpICoords[0].x - reactIface.x;
    const double iy = mol1.tmpICoords[0].y - reactIface.y;
    const double iz = mol1.tmpICoords[0].z - reactIface.z;
    EXPECT_NEAR(ix, 0.0, 1e-8) << "on-axis interface must have no x component left";
    EXPECT_NEAR(iy, 0.0, 1e-8) << "on-axis interface must have no y component left";
    EXPECT_NEAR(iz, axisLen, 1e-8) << "on-axis interface must lie along +z";

    // reactIface is the rotation origin and must never be modified.
    EXPECT_DOUBLE_EQ(reactIface.x, origIface.x) << "reactIface.x must not be modified";
    EXPECT_DOUBLE_EQ(reactIface.y, origIface.y) << "reactIface.y must not be modified";
    EXPECT_DOUBLE_EQ(reactIface.z, origIface.z) << "reactIface.z must not be modified";
}

// -----------------------------------------------------------------------------
// Test 5: the transformation is a rigid body rotation about reactIface, applied
//         identically to both molecules -> all internal distances survive.
// -----------------------------------------------------------------------------
void tfm_test_preserves_rigid_geometry()
{
    std::cerr << "\n[TEST] tfm_test_preserves_rigid_geometry\n"
              << "  Source file:   src/reactions/transform.cpp\n"
              << "  Function:      transform()\n"
              << "  Scenario:      an arbitrary axis (2,-1,3) and a non-origin\n"
              << "                 reactIface, with two multi-interface molecules.\n"
              << "  Pass criteria: (a) every distance to reactIface is preserved,\n"
              << "                 (b) the mol1-mol2 COM distance is preserved,\n"
              << "                 (c) intra-molecular distances are preserved.\n";

    Coord reactIface { -3.0, 4.0, 1.0 };

    Molecule mol1 = tfm_molecule(Coord { 0.0, 0.0, 0.0 },
        { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 2.0, -1.0 } });
    Molecule mol2 = tfm_molecule(Coord { -6.0, 7.0, 4.0 },
        { Coord { -5.0, 7.0, 4.0 } });

    // Capture the pre-rotation geometry.
    const double d1Com = tfm_distance(mol1.tmpComCoord, reactIface);
    const double d1If0 = tfm_distance(mol1.tmpICoords[0], reactIface);
    const double d1If1 = tfm_distance(mol1.tmpICoords[1], reactIface);
    const double d2Com = tfm_distance(mol2.tmpComCoord, reactIface);
    const double d2If0 = tfm_distance(mol2.tmpICoords[0], reactIface);
    const double dCross = tfm_distance(mol1.tmpComCoord, mol2.tmpComCoord);
    const double dIntra1 = tfm_distance(mol1.tmpComCoord, mol1.tmpICoords[1]);

    Vector axis = tfm_axis(2.0, -1.0, 3.0);

    std::cerr << "  Calling transform() with axis = (2,-1,3)...\n";
    transform(reactIface, mol1, mol2, axis);

    tfm_print("mol1.tmpComCoord", mol1.tmpComCoord);
    tfm_print("mol2.tmpComCoord", mol2.tmpComCoord);

    // (a) rotation about reactIface -> radial distances unchanged.
    EXPECT_NEAR(tfm_distance(mol1.tmpComCoord, reactIface), d1Com, 1e-8)
        << "mol1 COM distance to reactIface must be conserved";
    EXPECT_NEAR(tfm_distance(mol1.tmpICoords[0], reactIface), d1If0, 1e-8)
        << "mol1 iface0 distance to reactIface must be conserved";
    EXPECT_NEAR(tfm_distance(mol1.tmpICoords[1], reactIface), d1If1, 1e-8)
        << "mol1 iface1 distance to reactIface must be conserved";
    EXPECT_NEAR(tfm_distance(mol2.tmpComCoord, reactIface), d2Com, 1e-8)
        << "mol2 COM distance to reactIface must be conserved";
    EXPECT_NEAR(tfm_distance(mol2.tmpICoords[0], reactIface), d2If0, 1e-8)
        << "mol2 iface0 distance to reactIface must be conserved";

    // (b) both molecules get the *same* quaternion -> their separation is fixed.
    EXPECT_NEAR(tfm_distance(mol1.tmpComCoord, mol2.tmpComCoord), dCross, 1e-8)
        << "mol1-mol2 COM separation must be conserved (same rotation applied)";

    // (c) internal rigidity of a molecule.
    EXPECT_NEAR(tfm_distance(mol1.tmpComCoord, mol1.tmpICoords[1]), dIntra1, 1e-8)
        << "mol1 COM-to-interface distance must be conserved";

    // Sanity: something actually moved (this axis is not on z).
    const bool moved = std::fabs(mol1.tmpComCoord.x - 0.0) > kTfmTol
        || std::fabs(mol1.tmpComCoord.y - 0.0) > kTfmTol
        || std::fabs(mol1.tmpComCoord.z - 0.0) > kTfmTol;
    EXPECT_TRUE(moved) << "a non-z axis must actually produce a rotation";
}

// -----------------------------------------------------------------------------
// Test 6: the *real* coordinates must not be touched, only the tmp ones.
// -----------------------------------------------------------------------------
void tfm_test_real_coords_untouched()
{
    std::cerr << "\n[TEST] tfm_test_real_coords_untouched\n"
              << "  Source file:   src/reactions/transform.cpp\n"
              << "  Function:      transform()\n"
              << "  Scenario:      a rotating axis (0,1,0) is applied while the\n"
              << "                 molecules also carry real comCoord values.\n"
              << "  Pass criteria: Molecule::comCoord is bitwise unchanged while\n"
              << "                 Molecule::tmpComCoord has moved.\n";

    Coord reactIface { 0.0, 0.0, 0.0 };
    Molecule mol1 = tfm_molecule(Coord { 0.0, 0.0, 5.0 }, { Coord { 0.0, 0.0, 6.0 } });
    Molecule mol2 = tfm_molecule(Coord { 0.0, 0.0, -5.0 }, {});

    // Remember the real coordinates.
    const Coord realCom1 = mol1.comCoord;
    const Coord realCom2 = mol2.comCoord;

    Vector axis = tfm_axis(0.0, 1.0, 0.0); // 90 deg away from z -> real rotation

    std::cerr << "  Calling transform() with axis = (0,1,0)...\n";
    transform(reactIface, mol1, mol2, axis);

    tfm_print("mol1.comCoord    (real, must be 0,0,5)", mol1.comCoord);
    tfm_print("mol1.tmpComCoord (temporary, rotated) ", mol1.tmpComCoord);

    // Real coordinates untouched.
    EXPECT_DOUBLE_EQ(mol1.comCoord.x, realCom1.x) << "real comCoord.x must not change";
    EXPECT_DOUBLE_EQ(mol1.comCoord.y, realCom1.y) << "real comCoord.y must not change";
    EXPECT_DOUBLE_EQ(mol1.comCoord.z, realCom1.z) << "real comCoord.z must not change";
    EXPECT_DOUBLE_EQ(mol2.comCoord.z, realCom2.z) << "real comCoord.z must not change (mol2)";

    // Temporary coordinates did move: (0,0,5) rotated -90 deg about (zhat x yhat)
    // = -xhat leaves the z-axis, so the tmp z component can no longer be 5.
    EXPECT_GT(std::fabs(mol1.tmpComCoord.z - 5.0), 1e-6)
        << "tmpComCoord should have been rotated away from (0,0,5)";
    // ... but its distance from the rotation origin is still 5.
    EXPECT_NEAR(tfm_distance(mol1.tmpComCoord, reactIface), 5.0, 1e-8)
        << "rotation must preserve the radius of the tmp coordinate";
}

// -----------------------------------------------------------------------------
// Test 7: robustness - molecules with no temporary interface coordinates.
// -----------------------------------------------------------------------------
void tfm_test_handles_empty_interface_lists()
{
    std::cerr << "\n[TEST] tfm_test_handles_empty_interface_lists\n"
              << "  Source file:   src/reactions/transform.cpp\n"
              << "  Function:      transform()\n"
              << "  Scenario:      both molecules have empty tmpICoords vectors.\n"
              << "  Pass criteria: no crash; the COMs are still rotated rigidly\n"
              << "                 and the interface lists stay empty.\n";

    Coord reactIface { 0.0, 0.0, 0.0 };
    Molecule mol1 = tfm_molecule(Coord { 3.0, 0.0, 0.0 }, {});
    Molecule mol2 = tfm_molecule(Coord { 0.0, 4.0, 0.0 }, {});

    Vector axis = tfm_axis(1.0, 0.0, 0.0); // -90 deg about +y

    std::cerr << "  Calling transform() with empty tmpICoords...\n";
    transform(reactIface, mol1, mol2, axis);

    tfm_print("mol1.tmpComCoord (expect 0,0,3)", mol1.tmpComCoord);
    tfm_print("mol2.tmpComCoord (expect 0,4,0)", mol2.tmpComCoord);

    EXPECT_TRUE(mol1.tmpICoords.empty()) << "mol1 interface list must remain empty";
    EXPECT_TRUE(mol2.tmpICoords.empty()) << "mol2 interface list must remain empty";

    // (x,y,z) -> (-z, y, x):  (3,0,0) -> (0,0,3) and (0,4,0) -> (0,4,0).
    EXPECT_NEAR(mol1.tmpComCoord.x, 0.0, kTfmTol) << "mol1 COM x";
    EXPECT_NEAR(mol1.tmpComCoord.y, 0.0, kTfmTol) << "mol1 COM y";
    EXPECT_NEAR(mol1.tmpComCoord.z, 3.0, kTfmTol) << "mol1 COM z (axis mapped to +z)";
    EXPECT_NEAR(mol2.tmpComCoord.x, 0.0, kTfmTol) << "mol2 COM x (on rotation axis)";
    EXPECT_NEAR(mol2.tmpComCoord.y, 4.0, kTfmTol) << "mol2 COM y (on rotation axis)";
    EXPECT_NEAR(mol2.tmpComCoord.z, 0.0, kTfmTol) << "mol2 COM z (on rotation axis)";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - one TEST per scenario so every case runs and reports
// independently even if an earlier one fails (all assertions are non-fatal).
// -----------------------------------------------------------------------------
TEST(TransformAssociation, NoRotationWhenAxisIsPlusZ) { tfm_test_no_rotation_when_axis_is_plus_z(); }
TEST(TransformAssociation, NoRotationWhenAxisIsMinusZ) { tfm_test_no_rotation_when_axis_is_minus_z(); }
TEST(TransformAssociation, XAxisExplicitRotation) { tfm_test_x_axis_explicit_rotation(); }
TEST(TransformAssociation, ArbitraryAxisMapsOntoZAxis) { tfm_test_arbitrary_axis_maps_onto_z_axis(); }
TEST(TransformAssociation, PreservesRigidGeometry) { tfm_test_preserves_rigid_geometry(); }
TEST(TransformAssociation, RealCoordsUntouched) { tfm_test_real_coords_untouched(); }
TEST(TransformAssociation, HandlesEmptyInterfaceLists) { tfm_test_handles_empty_interface_lists(); }