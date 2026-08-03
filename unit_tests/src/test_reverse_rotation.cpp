/*! \file test_reverse_rotation.cpp
 *
 * ### Unit test for src/reactions/reverse_rotation.cpp
 *
 * Function under test:
 * \code
 *   void reverse_rotation(Coord& reactIface1, Molecule& reactMol1, Molecule& reactMol2,
 *                         Complex& reactCom1, Complex& reactCom2,
 *                         Quat rotQuatPos, Quat rotQuatNeg,
 *                         std::vector<Molecule>& moleculeList);
 * \endcode
 *
 * The routine is used during association to *undo* a rotation that was applied
 * with `rotate()`.  It does exactly three things:
 *
 *   1. inverts both rotation quaternions (`Quat::inverse()`),
 *   2. rotates the first complex (`reactCom1`) about `reactIface1` with the
 *      inverted "positive" quaternion,
 *   3. rotates the second complex (`reactCom2`) about `reactIface1` with the
 *      inverted "negative" quaternion.
 *
 * `rotate()` operates on the *temporary* association coordinates
 * (`Molecule::tmpComCoord` / `Molecule::tmpICoords`), so all of the checks below
 * inspect those, and additionally verify that the permanent coordinates
 * (`comCoord` / `interfaceList[i].coord`) are left alone.
 *
 * Because the sign convention of `Quat::rotate()` is an implementation detail,
 * the "expected" coordinates are computed here with the very same `Quat` class
 * (rotating a point about the same origin with the inverted quaternion).  This
 * keeps the assertions convention independent while still being an exact,
 * analytic check that `reverse_rotation()` applies the *inverse* of each
 * quaternion, to the *correct* complex.
 */

#include "classes/class_Quat.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

//! Numerical tolerance used for all floating point comparisons in this file.
constexpr double kRrTol = 1e-9;

/*! \brief Build a unit rotation quaternion from an axis and an angle (radians).
 *
 * Q = ( cos(angle/2), sin(angle/2) * axisHat )
 */
Quat rr_make_unit_quat(double angle, double ax, double ay, double az)
{
    const double norm = std::sqrt(ax * ax + ay * ay + az * az);
    ax /= norm;
    ay /= norm;
    az /= norm;
    const double half = 0.5 * angle;
    return Quat(std::cos(half), std::sin(half) * ax, std::sin(half) * ay, std::sin(half) * az);
}

/*! \brief Create a minimal Molecule with a COM, a list of interfaces, and the
 *         temporary association coordinates already initialised.
 *
 * \param[in] comIndex index of the parent Complex (Molecule::myComIndex)
 * \param[in] com      center-of-mass coordinate
 * \param[in] ifaces   interface coordinates (absolute, not relative)
 */
Molecule rr_make_molecule(int comIndex, const Coord& com, const std::vector<Coord>& ifaces)
{
    Molecule mol;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.comCoord = com;

    for (const auto& oneIface : ifaces) {
        Molecule::Iface iface;
        iface.coord = oneIface;
        mol.interfaceList.push_back(iface);
    }

    // Association always works on the temporary coordinates; mirror the real
    // ones into tmpComCoord / tmpICoords the same way associate() does.
    mol.set_tmp_association_coords();
    return mol;
}

/*! \brief Create a Complex owning the supplied moleculeList indices. */
Complex rr_make_complex(const Coord& com, const std::vector<int>& members)
{
    Complex targCom;
    targCom.comCoord = com;
    targCom.tmpComCoord = com;
    targCom.memberList = members;

    // Diffusion constants are irrelevant here but are given sane values.
    targCom.D.x = 1.0;
    targCom.D.y = 1.0;
    targCom.D.z = 1.0;
    targCom.Dr.x = 0.01;
    targCom.Dr.y = 0.01;
    targCom.Dr.z = 0.01;
    return targCom;
}

/*! \brief Rotate a single point about an origin with the given quaternion.
 *
 * Uses the project's own Quat::rotate(), so the expected values follow the same
 * sign/handedness convention as the code under test.
 */
Coord rr_rotate_point(const Coord& point, const Coord& origin, Quat q)
{
    Coord diff = point - origin; // vector from the rotation origin to the point
    Vector vec(diff);
    q.rotate(vec);
    return Coord { origin.x + vec.x, origin.y + vec.y, origin.z + vec.z };
}

/*! \brief EXPECT that two Coords agree component-wise (non-fatal). */
void rr_expect_coord_near(const Coord& actual, const Coord& expected, const std::string& what)
{
    EXPECT_NEAR(actual.x, expected.x, kRrTol) << what << ": x component mismatch";
    EXPECT_NEAR(actual.y, expected.y, kRrTol) << what << ": y component mismatch";
    EXPECT_NEAR(actual.z, expected.z, kRrTol) << what << ": z component mismatch";
}

/*! \brief Convenience printer used for the verbose console output. */
void rr_print_coord(const std::string& label, const Coord& c)
{
    std::cerr << "      " << label << " = (" << c.x << ", " << c.y << ", " << c.z << ")\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: identity quaternions must leave every coordinate untouched.
//         The inverse of the identity quaternion is the identity, so a
//         "reverse rotation" by the identity is a no-op.
// -----------------------------------------------------------------------------
void test_rr_identity_quaternions_are_a_noop()
{
    std::cerr << "\n[TEST] test_rr_identity_quaternions_are_a_noop\n"
              << "  Source file:   src/reactions/reverse_rotation.cpp\n"
              << "  Function:      reverse_rotation()\n"
              << "  Scenario:      both quaternions are the identity (1,0,0,0).\n"
              << "  Pass criteria: all temporary coordinates of both complexes are\n"
              << "                 numerically unchanged.\n";

    // Two molecules, one per complex, each with a single interface.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rr_make_molecule(0, Coord { 1.0, 2.0, 3.0 }, { Coord { 1.0, 2.0, 4.0 } }));
    moleculeList.push_back(rr_make_molecule(1, Coord { -3.0, 0.5, 2.0 }, { Coord { -3.0, 0.5, 1.0 } }));

    Complex reactCom1 = rr_make_complex(moleculeList[0].comCoord, { 0 });
    Complex reactCom2 = rr_make_complex(moleculeList[1].comCoord, { 1 });

    // Remember the starting temporary coordinates.
    const Coord com0Before = moleculeList[0].tmpComCoord;
    const Coord iface0Before = moleculeList[0].tmpICoords[0];
    const Coord com1Before = moleculeList[1].tmpComCoord;
    const Coord iface1Before = moleculeList[1].tmpICoords[0];

    Coord reactIface1 { 0.0, 0.0, 0.0 }; // rotation origin
    Quat identityPos(1.0, 0.0, 0.0, 0.0);
    Quat identityNeg(1.0, 0.0, 0.0, 0.0);

    std::cerr << "    Calling reverse_rotation() with identity quaternions...\n";
    reverse_rotation(reactIface1, moleculeList[0], moleculeList[1], reactCom1, reactCom2,
        identityPos, identityNeg, moleculeList);

    rr_expect_coord_near(moleculeList[0].tmpComCoord, com0Before, "complex 1 molecule tmpComCoord");
    rr_expect_coord_near(moleculeList[0].tmpICoords[0], iface0Before, "complex 1 molecule tmpICoords[0]");
    rr_expect_coord_near(moleculeList[1].tmpComCoord, com1Before, "complex 2 molecule tmpComCoord");
    rr_expect_coord_near(moleculeList[1].tmpICoords[0], iface1Before, "complex 2 molecule tmpICoords[0]");

    rr_print_coord("complex 1 tmpComCoord after", moleculeList[0].tmpComCoord);
    rr_print_coord("complex 2 tmpComCoord after", moleculeList[1].tmpComCoord);
}

// -----------------------------------------------------------------------------
// Test 2: a single non-trivial quaternion is applied *inverted* and about the
//         supplied origin (reactIface1), and only the temporary coordinates are
//         touched.
// -----------------------------------------------------------------------------
void test_rr_applies_inverse_about_given_origin()
{
    std::cerr << "\n[TEST] test_rr_applies_inverse_about_given_origin\n"
              << "  Source file:   src/reactions/reverse_rotation.cpp\n"
              << "  Function:      reverse_rotation()\n"
              << "  Scenario:      rotQuatPos = 40 deg about z, rotQuatNeg = identity,\n"
              << "                 rotation origin is a non-zero point.\n"
              << "  Pass criteria: complex 1's temporary coordinates equal the result\n"
              << "                 of rotating them with rotQuatPos.inverse() about\n"
              << "                 that origin; complex 2 is untouched; permanent\n"
              << "                 coordinates are untouched; the origin is untouched.\n";

    // Complex 1 owns two molecules; the first one has two interfaces so that we
    // also verify that every interface of every member is rotated.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rr_make_molecule(0, Coord { 4.0, 0.0, 1.0 },
        { Coord { 5.0, 0.0, 1.0 }, Coord { 4.0, 1.0, 1.0 } }));
    moleculeList.push_back(rr_make_molecule(0, Coord { 3.0, -2.0, 1.0 }, { Coord { 3.0, -3.0, 1.0 } }));
    moleculeList.push_back(rr_make_molecule(1, Coord { -6.0, 2.0, -1.0 }, { Coord { -6.0, 2.0, 0.0 } }));

    Complex reactCom1 = rr_make_complex(Coord { 3.5, -1.0, 1.0 }, { 0, 1 });
    Complex reactCom2 = rr_make_complex(Coord { -6.0, 2.0, -1.0 }, { 2 });

    // Rotation origin: deliberately not the coordinate system origin.
    Coord reactIface1 { 1.0, -2.0, 0.5 };
    const Coord originBefore = reactIface1;

    Quat rotQuatPos = rr_make_unit_quat(40.0 * M_PI / 180.0, 0.0, 0.0, 1.0);
    Quat rotQuatNeg(1.0, 0.0, 0.0, 0.0); // identity -> complex 2 must not move

    // Compute the expected coordinates using the *inverse* quaternion.
    Quat invPos = rotQuatPos.inverse();
    const Coord expMol0Com = rr_rotate_point(moleculeList[0].tmpComCoord, reactIface1, invPos);
    const Coord expMol0If0 = rr_rotate_point(moleculeList[0].tmpICoords[0], reactIface1, invPos);
    const Coord expMol0If1 = rr_rotate_point(moleculeList[0].tmpICoords[1], reactIface1, invPos);
    const Coord expMol1Com = rr_rotate_point(moleculeList[1].tmpComCoord, reactIface1, invPos);
    const Coord expMol1If0 = rr_rotate_point(moleculeList[1].tmpICoords[0], reactIface1, invPos);

    // Snapshot of things that must not change.
    const Coord mol0RealCom = moleculeList[0].comCoord;
    const Coord mol0RealIf0 = moleculeList[0].interfaceList[0].coord;
    const Coord mol2TmpCom = moleculeList[2].tmpComCoord;
    const Coord mol2TmpIf0 = moleculeList[2].tmpICoords[0];

    std::cerr << "    Calling reverse_rotation()...\n";
    reverse_rotation(reactIface1, moleculeList[0], moleculeList[2], reactCom1, reactCom2,
        rotQuatPos, rotQuatNeg, moleculeList);

    // (a) every member molecule/interface of complex 1 got the inverse rotation
    rr_expect_coord_near(moleculeList[0].tmpComCoord, expMol0Com, "member 0 tmpComCoord");
    rr_expect_coord_near(moleculeList[0].tmpICoords[0], expMol0If0, "member 0 tmpICoords[0]");
    rr_expect_coord_near(moleculeList[0].tmpICoords[1], expMol0If1, "member 0 tmpICoords[1]");
    rr_expect_coord_near(moleculeList[1].tmpComCoord, expMol1Com, "member 1 tmpComCoord");
    rr_expect_coord_near(moleculeList[1].tmpICoords[0], expMol1If0, "member 1 tmpICoords[0]");

    // (b) a real rotation actually happened (guards against a silent no-op)
    EXPECT_GT(std::fabs(moleculeList[0].tmpComCoord.y - mol0RealCom.y), 1e-6)
        << "a 40 degree rotation about z should move the member molecule";

    // (c) complex 2 was rotated with the identity, i.e. not moved at all
    rr_expect_coord_near(moleculeList[2].tmpComCoord, mol2TmpCom, "complex 2 tmpComCoord (identity)");
    rr_expect_coord_near(moleculeList[2].tmpICoords[0], mol2TmpIf0, "complex 2 tmpICoords[0] (identity)");

    // (d) permanent (non-temporary) coordinates must be untouched
    rr_expect_coord_near(moleculeList[0].comCoord, mol0RealCom, "member 0 permanent comCoord");
    rr_expect_coord_near(moleculeList[0].interfaceList[0].coord, mol0RealIf0,
        "member 0 permanent interfaceList[0].coord");

    // (e) the rotation origin itself must be untouched
    rr_expect_coord_near(reactIface1, originBefore, "rotation origin reactIface1");

    rr_print_coord("member 0 tmpComCoord expected", expMol0Com);
    rr_print_coord("member 0 tmpComCoord actual  ", moleculeList[0].tmpComCoord);
}

// -----------------------------------------------------------------------------
// Test 3: the two quaternions must not be swapped: rotQuatPos^-1 belongs to
//         reactCom1 and rotQuatNeg^-1 belongs to reactCom2.
// -----------------------------------------------------------------------------
void test_rr_quaternions_map_to_correct_complexes()
{
    std::cerr << "\n[TEST] test_rr_quaternions_map_to_correct_complexes\n"
              << "  Source file:   src/reactions/reverse_rotation.cpp\n"
              << "  Function:      reverse_rotation()\n"
              << "  Scenario:      rotQuatPos = 30 deg about z, rotQuatNeg = -50 deg about x\n"
              << "                 (deliberately different axes and angles).\n"
              << "  Pass criteria: complex 1 follows rotQuatPos.inverse(), complex 2\n"
              << "                 follows rotQuatNeg.inverse() (not swapped).\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rr_make_molecule(0, Coord { 2.0, 1.0, 0.0 }, { Coord { 3.0, 1.0, 0.0 } }));
    moleculeList.push_back(rr_make_molecule(1, Coord { 0.0, 3.0, 2.0 }, { Coord { 0.0, 4.0, 2.0 } }));

    Complex reactCom1 = rr_make_complex(moleculeList[0].comCoord, { 0 });
    Complex reactCom2 = rr_make_complex(moleculeList[1].comCoord, { 1 });

    Coord reactIface1 { 0.0, 0.0, 0.0 };

    Quat rotQuatPos = rr_make_unit_quat(30.0 * M_PI / 180.0, 0.0, 0.0, 1.0);
    Quat rotQuatNeg = rr_make_unit_quat(-50.0 * M_PI / 180.0, 1.0, 0.0, 0.0);

    Quat invPos = rotQuatPos.inverse();
    Quat invNeg = rotQuatNeg.inverse();

    // Expected results, and the "wrong" (swapped) results used as a contrast.
    const Coord expCom1 = rr_rotate_point(moleculeList[0].tmpComCoord, reactIface1, invPos);
    const Coord expCom2 = rr_rotate_point(moleculeList[1].tmpComCoord, reactIface1, invNeg);
    const Coord swappedCom1 = rr_rotate_point(moleculeList[0].tmpComCoord, reactIface1, invNeg);
    const Coord swappedCom2 = rr_rotate_point(moleculeList[1].tmpComCoord, reactIface1, invPos);

    std::cerr << "    Calling reverse_rotation()...\n";
    reverse_rotation(reactIface1, moleculeList[0], moleculeList[1], reactCom1, reactCom2,
        rotQuatPos, rotQuatNeg, moleculeList);

    rr_expect_coord_near(moleculeList[0].tmpComCoord, expCom1,
        "complex 1 must use rotQuatPos.inverse()");
    rr_expect_coord_near(moleculeList[1].tmpComCoord, expCom2,
        "complex 2 must use rotQuatNeg.inverse()");

    // Sanity: the two candidate results really are distinguishable, so the
    // assertions above are meaningful.
    const double contrast1 = std::fabs(expCom1.y - swappedCom1.y) + std::fabs(expCom1.z - swappedCom1.z);
    const double contrast2 = std::fabs(expCom2.y - swappedCom2.y) + std::fabs(expCom2.z - swappedCom2.z);
    EXPECT_GT(contrast1, 1e-6) << "test setup should distinguish pos/neg quaternions for complex 1";
    EXPECT_GT(contrast2, 1e-6) << "test setup should distinguish pos/neg quaternions for complex 2";

    rr_print_coord("complex 1 tmpComCoord actual", moleculeList[0].tmpComCoord);
    rr_print_coord("complex 2 tmpComCoord actual", moleculeList[1].tmpComCoord);
}

// -----------------------------------------------------------------------------
// Test 4: round trip.  rotate() followed by reverse_rotation() with the same
//         quaternions must restore the original temporary coordinates.  This is
//         exactly how association uses the function.
// -----------------------------------------------------------------------------
void test_rr_round_trip_restores_coordinates()
{
    std::cerr << "\n[TEST] test_rr_round_trip_restores_coordinates\n"
              << "  Source file:   src/reactions/reverse_rotation.cpp\n"
              << "  Function:      reverse_rotation() (paired with rotate())\n"
              << "  Scenario:      forward rotate() both complexes, then undo with\n"
              << "                 reverse_rotation() using the same quaternions.\n"
              << "  Pass criteria: (a) the forward rotation actually moved things,\n"
              << "                 (b) the coordinates are restored to their original\n"
              << "                     values to within 1e-9.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rr_make_molecule(0, Coord { 5.0, -1.0, 2.0 },
        { Coord { 6.0, -1.0, 2.0 }, Coord { 5.0, -1.0, 3.0 } }));
    moleculeList.push_back(rr_make_molecule(1, Coord { -4.0, 3.0, -2.0 }, { Coord { -4.0, 2.0, -2.0 } }));

    Complex reactCom1 = rr_make_complex(moleculeList[0].comCoord, { 0 });
    Complex reactCom2 = rr_make_complex(moleculeList[1].comCoord, { 1 });

    // Originals, before any rotation at all.
    const Coord origCom1 = moleculeList[0].tmpComCoord;
    const Coord origIf1a = moleculeList[0].tmpICoords[0];
    const Coord origIf1b = moleculeList[0].tmpICoords[1];
    const Coord origCom2 = moleculeList[1].tmpComCoord;
    const Coord origIf2 = moleculeList[1].tmpICoords[0];

    Coord reactIface1 { 0.5, 0.5, 0.5 };

    Quat rotQuatPos = rr_make_unit_quat(70.0 * M_PI / 180.0, 1.0, 2.0, -1.0);
    Quat rotQuatNeg = rr_make_unit_quat(-70.0 * M_PI / 180.0, 1.0, 2.0, -1.0);

    // ---- forward rotation, exactly as the association routines do it -------
    std::cerr << "    Applying forward rotate() to both complexes...\n";
    rotate(reactIface1, rotQuatPos, reactCom1, moleculeList);
    rotate(reactIface1, rotQuatNeg, reactCom2, moleculeList);

    const double moved1 = std::fabs(moleculeList[0].tmpComCoord.x - origCom1.x)
        + std::fabs(moleculeList[0].tmpComCoord.y - origCom1.y)
        + std::fabs(moleculeList[0].tmpComCoord.z - origCom1.z);
    const double moved2 = std::fabs(moleculeList[1].tmpComCoord.x - origCom2.x)
        + std::fabs(moleculeList[1].tmpComCoord.y - origCom2.y)
        + std::fabs(moleculeList[1].tmpComCoord.z - origCom2.z);
    EXPECT_GT(moved1, 1e-6) << "forward rotation should have displaced complex 1";
    EXPECT_GT(moved2, 1e-6) << "forward rotation should have displaced complex 2";
    std::cerr << "      forward displacement (complex 1) = " << moved1 << "\n"
              << "      forward displacement (complex 2) = " << moved2 << "\n";

    // ---- undo it ----------------------------------------------------------
    std::cerr << "    Calling reverse_rotation() to undo the forward rotation...\n";
    reverse_rotation(reactIface1, moleculeList[0], moleculeList[1], reactCom1, reactCom2,
        rotQuatPos, rotQuatNeg, moleculeList);

    rr_expect_coord_near(moleculeList[0].tmpComCoord, origCom1, "restored complex 1 tmpComCoord");
    rr_expect_coord_near(moleculeList[0].tmpICoords[0], origIf1a, "restored complex 1 tmpICoords[0]");
    rr_expect_coord_near(moleculeList[0].tmpICoords[1], origIf1b, "restored complex 1 tmpICoords[1]");
    rr_expect_coord_near(moleculeList[1].tmpComCoord, origCom2, "restored complex 2 tmpComCoord");
    rr_expect_coord_near(moleculeList[1].tmpICoords[0], origIf2, "restored complex 2 tmpICoords[0]");

    rr_print_coord("complex 1 tmpComCoord original", origCom1);
    rr_print_coord("complex 1 tmpComCoord restored", moleculeList[0].tmpComCoord);
}

// -----------------------------------------------------------------------------
// Test 5: the quaternion arguments are taken *by value*, so inverting them
//         inside reverse_rotation() must not be visible to the caller.
// -----------------------------------------------------------------------------
void test_rr_quaternion_arguments_are_by_value()
{
    std::cerr << "\n[TEST] test_rr_quaternion_arguments_are_by_value\n"
              << "  Source file:   src/reactions/reverse_rotation.cpp\n"
              << "  Function:      reverse_rotation()\n"
              << "  Scenario:      inspect the caller's quaternions after the call.\n"
              << "  Pass criteria: rotQuatPos/rotQuatNeg still hold their original\n"
              << "                 w,x,y,z (they are passed by value and inverted\n"
              << "                 only inside the function).\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rr_make_molecule(0, Coord { 1.0, 0.0, 0.0 }, { Coord { 2.0, 0.0, 0.0 } }));
    moleculeList.push_back(rr_make_molecule(1, Coord { -1.0, 0.0, 0.0 }, { Coord { -2.0, 0.0, 0.0 } }));

    Complex reactCom1 = rr_make_complex(moleculeList[0].comCoord, { 0 });
    Complex reactCom2 = rr_make_complex(moleculeList[1].comCoord, { 1 });

    Coord reactIface1 { 0.0, 0.0, 0.0 };

    Quat rotQuatPos = rr_make_unit_quat(25.0 * M_PI / 180.0, 0.0, 1.0, 0.0);
    Quat rotQuatNeg = rr_make_unit_quat(-25.0 * M_PI / 180.0, 0.0, 1.0, 0.0);

    // Copies of the caller's quaternions before the call.
    const Quat posBefore = rotQuatPos;
    const Quat negBefore = rotQuatNeg;

    std::cerr << "    Calling reverse_rotation()...\n";
    reverse_rotation(reactIface1, moleculeList[0], moleculeList[1], reactCom1, reactCom2,
        rotQuatPos, rotQuatNeg, moleculeList);

    EXPECT_NEAR(rotQuatPos.w, posBefore.w, kRrTol) << "caller's rotQuatPos.w must be unchanged";
    EXPECT_NEAR(rotQuatPos.x, posBefore.x, kRrTol) << "caller's rotQuatPos.x must be unchanged";
    EXPECT_NEAR(rotQuatPos.y, posBefore.y, kRrTol) << "caller's rotQuatPos.y must be unchanged";
    EXPECT_NEAR(rotQuatPos.z, posBefore.z, kRrTol) << "caller's rotQuatPos.z must be unchanged";

    EXPECT_NEAR(rotQuatNeg.w, negBefore.w, kRrTol) << "caller's rotQuatNeg.w must be unchanged";
    EXPECT_NEAR(rotQuatNeg.x, negBefore.x, kRrTol) << "caller's rotQuatNeg.x must be unchanged";
    EXPECT_NEAR(rotQuatNeg.y, negBefore.y, kRrTol) << "caller's rotQuatNeg.y must be unchanged";
    EXPECT_NEAR(rotQuatNeg.z, negBefore.z, kRrTol) << "caller's rotQuatNeg.z must be unchanged";

    std::cerr << "      rotQuatPos after call = (" << rotQuatPos.w << ", " << rotQuatPos.x << ", "
              << rotQuatPos.y << ", " << rotQuatPos.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 6: rotations are rigid.  Undoing a rotation must preserve every
//         interface-to-COM distance and the distance to the rotation origin.
// -----------------------------------------------------------------------------
void test_rr_preserves_rigid_geometry()
{
    std::cerr << "\n[TEST] test_rr_preserves_rigid_geometry\n"
              << "  Source file:   src/reactions/reverse_rotation.cpp\n"
              << "  Function:      reverse_rotation()\n"
              << "  Scenario:      apply a reverse rotation to a two-interface molecule.\n"
              << "  Pass criteria: |iface - COM| and |COM - origin| are conserved,\n"
              << "                 i.e. the operation is a rigid-body rotation.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rr_make_molecule(0, Coord { 3.0, 4.0, 0.0 },
        { Coord { 4.0, 4.0, 0.0 }, Coord { 3.0, 4.0, 2.0 } }));
    moleculeList.push_back(rr_make_molecule(1, Coord { -3.0, -4.0, 0.0 }, { Coord { -3.0, -5.0, 0.0 } }));

    Complex reactCom1 = rr_make_complex(moleculeList[0].comCoord, { 0 });
    Complex reactCom2 = rr_make_complex(moleculeList[1].comCoord, { 1 });

    Coord reactIface1 { 0.0, 0.0, 0.0 };

    // Distances before the call.
    auto dist = [](const Coord& a, const Coord& b) {
        const double dx = a.x - b.x;
        const double dy = a.y - b.y;
        const double dz = a.z - b.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };
    const double d0Before = dist(moleculeList[0].tmpICoords[0], moleculeList[0].tmpComCoord);
    const double d1Before = dist(moleculeList[0].tmpICoords[1], moleculeList[0].tmpComCoord);
    const double rBefore = dist(moleculeList[0].tmpComCoord, reactIface1);

    Quat rotQuatPos = rr_make_unit_quat(115.0 * M_PI / 180.0, -1.0, 0.5, 2.0);
    Quat rotQuatNeg = rr_make_unit_quat(15.0 * M_PI / 180.0, 0.0, 0.0, 1.0);

    std::cerr << "    Calling reverse_rotation()...\n";
    reverse_rotation(reactIface1, moleculeList[0], moleculeList[1], reactCom1, reactCom2,
        rotQuatPos, rotQuatNeg, moleculeList);

    const double d0After = dist(moleculeList[0].tmpICoords[0], moleculeList[0].tmpComCoord);
    const double d1After = dist(moleculeList[0].tmpICoords[1], moleculeList[0].tmpComCoord);
    const double rAfter = dist(moleculeList[0].tmpComCoord, reactIface1);

    EXPECT_NEAR(d0After, d0Before, 1e-9) << "interface 0 to COM distance must be conserved";
    EXPECT_NEAR(d1After, d1Before, 1e-9) << "interface 1 to COM distance must be conserved";
    EXPECT_NEAR(rAfter, rBefore, 1e-9) << "COM to rotation origin distance must be conserved";

    std::cerr << "      |iface0-COM|: before = " << d0Before << ", after = " << d0After << "\n"
              << "      |iface1-COM|: before = " << d1Before << ", after = " << d1After << "\n"
              << "      |COM-origin|: before = " << rBefore << ", after = " << rAfter << "\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is run inside its own TEST so that a
// failure in one does not stop the others (all assertions are non-fatal).
// -----------------------------------------------------------------------------
TEST(ReverseRotation, IdentityQuaternionsAreNoop) { test_rr_identity_quaternions_are_a_noop(); }
TEST(ReverseRotation, AppliesInverseAboutGivenOrigin) { test_rr_applies_inverse_about_given_origin(); }
TEST(ReverseRotation, QuaternionsMapToCorrectComplexes) { test_rr_quaternions_map_to_correct_complexes(); }
TEST(ReverseRotation, RoundTripRestoresCoordinates) { test_rr_round_trip_restores_coordinates(); }
TEST(ReverseRotation, QuaternionArgumentsAreByValue) { test_rr_quaternion_arguments_are_by_value(); }
TEST(ReverseRotation, PreservesRigidGeometry) { test_rr_preserves_rigid_geometry(); }