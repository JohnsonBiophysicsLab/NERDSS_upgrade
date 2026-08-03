/*! \file test_rotate.cpp
 *
 * ### Unit test for src/reactions/rotate.cpp
 *
 * Function under test:
 *
 *     void rotate(Coord& rotOrigin, Quat& rotQuat, Complex& targCom,
 *                 std::vector<Molecule>& moleculeList)
 *
 * The function performs a rigid-body quaternion rotation of an entire Complex
 * about an arbitrary origin, operating **only** on the temporary association
 * coordinates (`Molecule::tmpComCoord` and `Molecule::tmpICoords`) of the
 * molecules listed in `Complex::memberList`.
 *
 * The tests below verify the mathematical and bookkeeping properties that must
 * hold regardless of the internal quaternion sign/handedness convention:
 *
 *   1. An identity quaternion must leave every coordinate untouched.
 *   2. A 180 degree rotation about z maps (x,y,z) -> (-x,-y,z) relative to the
 *      rotation origin (this result is independent of rotation handedness).
 *   3. A point sitting exactly on the rotation origin is a fixed point.
 *   4. Rotation is rigid: all distances to the rotation origin and all pairwise
 *      distances between rotated points are preserved.
 *   5. Only molecules present in Complex::memberList are modified.
 *   6. The permanent coordinates (comCoord / interfaceList) are never touched;
 *      only the tmp* association coordinates change.
 *   7. Applying a 90 degree rotation twice is equivalent to one 180 degree
 *      rotation (composition sanity check, also handedness independent).
 *
 * Verbose progress information is written to stderr so a reader can follow
 * exactly what is being exercised and what the pass criteria are.
 */

#include "classes/class_Coord.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Quat.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (kept in an anonymous namespace so they cannot collide with
// helpers from other translation units in the combined test binary).
// -----------------------------------------------------------------------------
namespace {

//! Tolerance used for all floating point comparisons in this file.
constexpr double kRotateFnTol = 1e-10;

/*! \brief Build a Molecule whose permanent and temporary coordinates coincide.
 *
 * \param[in] com      Center-of-mass coordinate.
 * \param[in] ifaces   Interface coordinates (absolute, not relative).
 * \param[in] comIndex Index of the parent Complex.
 * \return A Molecule ready to be inserted into a moleculeList.
 */
Molecule rotatefn_make_molecule(const Coord& com, const std::vector<Coord>& ifaces,
    int comIndex)
{
    Molecule mol;
    mol.comCoord = com;
    mol.tmpComCoord = com; // rotate() operates on the tmp coords only

    mol.interfaceList.clear();
    mol.tmpICoords.clear();
    for (const auto& oneIface : ifaces) {
        Molecule::Iface iface;
        iface.coord = oneIface; // permanent copy (must remain untouched)
        mol.interfaceList.push_back(iface);
        mol.tmpICoords.push_back(oneIface); // temporary copy (gets rotated)
    }

    mol.myComIndex = comIndex;
    return mol;
}

/*! \brief Build a Complex owning the given member molecule indices. */
Complex rotatefn_make_complex(const std::vector<int>& members)
{
    Complex targCom;
    targCom.index = 0;
    targCom.memberList = members;
    // Diffusion constants are irrelevant for rotate(), but set sane values.
    targCom.D = Coord(1.0, 1.0, 1.0);
    targCom.Dr = Coord(0.01, 0.01, 0.01);
    return targCom;
}

/*! \brief Euclidean distance between two Coords. */
double rotatefn_dist(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Convenience: compare two Coords component-wise with a tolerance. */
void rotatefn_expect_coord_near(const Coord& actual, const Coord& expected,
    const char* what)
{
    EXPECT_NEAR(actual.x, expected.x, kRotateFnTol) << what << " : x component";
    EXPECT_NEAR(actual.y, expected.y, kRotateFnTol) << what << " : y component";
    EXPECT_NEAR(actual.z, expected.z, kRotateFnTol) << what << " : z component";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Identity quaternion => nothing moves.
// -----------------------------------------------------------------------------
void test_rotatefn_identity_quaternion()
{
    std::cerr << "\n[TEST] test_rotatefn_identity_quaternion\n"
              << "  Source file:   src/reactions/rotate.cpp\n"
              << "  Function:      rotate(Coord&, Quat&, Complex&, vector<Molecule>&)\n"
              << "  Scenario:      rotate the complex with the identity quaternion\n"
              << "                 Q = (1, 0, 0, 0) about the origin.\n"
              << "  Pass criteria: every tmpComCoord and tmpICoord is unchanged.\n";

    // One molecule with two interfaces at arbitrary positions.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rotatefn_make_molecule(
        Coord(3.0, -4.0, 5.0), { Coord(4.0, -4.0, 5.0), Coord(3.0, -3.0, 5.5) }, 0));

    Complex targCom = rotatefn_make_complex({ 0 });

    Coord rotOrigin(0.0, 0.0, 0.0);
    Quat identity(1.0, 0.0, 0.0, 0.0); // w = 1 => zero rotation angle

    std::cerr << "  Calling rotate()...\n";
    rotate(rotOrigin, identity, targCom, moleculeList);

    // The identity rotation must be a no-op.
    rotatefn_expect_coord_near(moleculeList[0].tmpComCoord, Coord(3.0, -4.0, 5.0),
        "tmpComCoord after identity rotation");
    rotatefn_expect_coord_near(moleculeList[0].tmpICoords[0], Coord(4.0, -4.0, 5.0),
        "tmpICoords[0] after identity rotation");
    rotatefn_expect_coord_near(moleculeList[0].tmpICoords[1], Coord(3.0, -3.0, 5.5),
        "tmpICoords[1] after identity rotation");

    std::cerr << "  Result tmpComCoord = (" << moleculeList[0].tmpComCoord.x << ", "
              << moleculeList[0].tmpComCoord.y << ", "
              << moleculeList[0].tmpComCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: 180 degree rotation about the z axis, rotation origin at (0,0,0).
//         (x, y, z) -> (-x, -y, z). This expectation is independent of the
//         handedness convention used inside Quat::rotate().
// -----------------------------------------------------------------------------
void test_rotatefn_180_degrees_about_z()
{
    std::cerr << "\n[TEST] test_rotatefn_180_degrees_about_z\n"
              << "  Source file:   src/reactions/rotate.cpp\n"
              << "  Function:      rotate()\n"
              << "  Scenario:      Q = (0, 0, 0, 1) i.e. a pi rotation about z,\n"
              << "                 rotation origin at (0,0,0).\n"
              << "  Pass criteria: (x,y,z) -> (-x,-y,z) for the COM and every\n"
              << "                 interface of every member molecule.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rotatefn_make_molecule(
        Coord(2.0, 1.0, 7.0), { Coord(3.0, 1.0, 7.0), Coord(2.0, -2.0, 7.0) }, 0));

    Complex targCom = rotatefn_make_complex({ 0 });

    Coord rotOrigin(0.0, 0.0, 0.0);
    // Unit quaternion with w = cos(pi/2) = 0 and axis component sin(pi/2) = 1.
    Quat rotZ180(0.0, 0.0, 0.0, 1.0);

    std::cerr << "  Calling rotate()...\n";
    rotate(rotOrigin, rotZ180, targCom, moleculeList);

    rotatefn_expect_coord_near(moleculeList[0].tmpComCoord, Coord(-2.0, -1.0, 7.0),
        "tmpComCoord after pi rotation about z");
    rotatefn_expect_coord_near(moleculeList[0].tmpICoords[0], Coord(-3.0, -1.0, 7.0),
        "tmpICoords[0] after pi rotation about z");
    rotatefn_expect_coord_near(moleculeList[0].tmpICoords[1], Coord(-2.0, 2.0, 7.0),
        "tmpICoords[1] after pi rotation about z");

    std::cerr << "  Result tmpComCoord = (" << moleculeList[0].tmpComCoord.x << ", "
              << moleculeList[0].tmpComCoord.y << ", "
              << moleculeList[0].tmpComCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 3: The rotation origin itself is a fixed point, and rotation about a
//         non-zero origin still gives the analytically expected pi result.
// -----------------------------------------------------------------------------
void test_rotatefn_nonzero_origin_is_fixed_point()
{
    std::cerr << "\n[TEST] test_rotatefn_nonzero_origin_is_fixed_point\n"
              << "  Source file:   src/reactions/rotate.cpp\n"
              << "  Function:      rotate()\n"
              << "  Scenario:      pi rotation about z with rotation origin at\n"
              << "                 (10, -5, 3). Molecule 0's COM sits exactly on\n"
              << "                 that origin; molecule 1 is offset from it.\n"
              << "  Pass criteria: the point on the origin does not move, and the\n"
              << "                 offset point is mirrored through the origin in\n"
              << "                 x and y while keeping its z value.\n";

    const Coord origin(10.0, -5.0, 3.0);

    std::vector<Molecule> moleculeList;
    // Molecule 0: COM exactly on the rotation origin, interface offset in x.
    moleculeList.push_back(rotatefn_make_molecule(
        origin, { Coord(origin.x + 2.0, origin.y, origin.z) }, 0));
    // Molecule 1: offset in x and y, and displaced in z.
    moleculeList.push_back(rotatefn_make_molecule(
        Coord(origin.x + 1.0, origin.y + 4.0, origin.z + 6.0), {}, 0));

    Complex targCom = rotatefn_make_complex({ 0, 1 });

    Coord rotOrigin = origin;
    Quat rotZ180(0.0, 0.0, 0.0, 1.0);

    std::cerr << "  Calling rotate()...\n";
    rotate(rotOrigin, rotZ180, targCom, moleculeList);

    // The point coincident with the rotation origin must be a fixed point.
    rotatefn_expect_coord_near(moleculeList[0].tmpComCoord, origin,
        "COM located on the rotation origin must not move");

    // Interface offset +2 in x becomes offset -2 in x.
    rotatefn_expect_coord_near(moleculeList[0].tmpICoords[0],
        Coord(origin.x - 2.0, origin.y, origin.z),
        "interface mirrored through the rotation origin");

    // Offset (+1, +4, +6) becomes (-1, -4, +6) relative to the origin.
    rotatefn_expect_coord_near(moleculeList[1].tmpComCoord,
        Coord(origin.x - 1.0, origin.y - 4.0, origin.z + 6.0),
        "second molecule COM mirrored in x and y, z preserved");

    std::cerr << "  Fixed point check OK; offset molecule now at ("
              << moleculeList[1].tmpComCoord.x << ", "
              << moleculeList[1].tmpComCoord.y << ", "
              << moleculeList[1].tmpComCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 4: Rigid-body property. A 90 degree rotation about an arbitrary unit
//         axis must preserve (a) each point's distance to the rotation origin
//         and (b) all pairwise distances between the rotated points.
// -----------------------------------------------------------------------------
void test_rotatefn_preserves_rigid_body_geometry()
{
    std::cerr << "\n[TEST] test_rotatefn_preserves_rigid_body_geometry\n"
              << "  Source file:   src/reactions/rotate.cpp\n"
              << "  Function:      rotate()\n"
              << "  Scenario:      90 degree rotation about the unit axis\n"
              << "                 (1,1,1)/sqrt(3) with origin (1,2,3).\n"
              << "  Pass criteria: distances to the rotation origin and all\n"
              << "                 pairwise distances are unchanged (rigid body).\n";

    const Coord origin(1.0, 2.0, 3.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rotatefn_make_molecule(
        Coord(4.0, 2.0, 3.0), { Coord(5.0, 2.0, 3.0), Coord(4.0, 4.0, 3.0) }, 0));
    moleculeList.push_back(rotatefn_make_molecule(
        Coord(1.0, 2.0, 9.0), { Coord(1.5, 2.5, 9.5) }, 0));

    Complex targCom = rotatefn_make_complex({ 0, 1 });

    // Record all pre-rotation distances.
    const double d0Origin = rotatefn_dist(moleculeList[0].tmpComCoord, origin);
    const double d1Origin = rotatefn_dist(moleculeList[1].tmpComCoord, origin);
    const double dIface0Origin = rotatefn_dist(moleculeList[0].tmpICoords[0], origin);
    const double dPair = rotatefn_dist(moleculeList[0].tmpComCoord,
        moleculeList[1].tmpComCoord);
    const double dComToIface = rotatefn_dist(moleculeList[0].tmpComCoord,
        moleculeList[0].tmpICoords[1]);

    Coord rotOrigin = origin;
    // 90 degree rotation: w = cos(45deg), axis*sin(45deg) with unit axis.
    const double halfAng = M_PI / 4.0;
    const double axisComp = std::sin(halfAng) / std::sqrt(3.0);
    Quat rotQuat(std::cos(halfAng), axisComp, axisComp, axisComp);

    std::cerr << "  Calling rotate()...\n";
    rotate(rotOrigin, rotQuat, targCom, moleculeList);

    // (a) Distances to the rotation origin must be invariant.
    EXPECT_NEAR(rotatefn_dist(moleculeList[0].tmpComCoord, origin), d0Origin, 1e-9)
        << "molecule 0 COM distance to rotation origin must be preserved";
    EXPECT_NEAR(rotatefn_dist(moleculeList[1].tmpComCoord, origin), d1Origin, 1e-9)
        << "molecule 1 COM distance to rotation origin must be preserved";
    EXPECT_NEAR(rotatefn_dist(moleculeList[0].tmpICoords[0], origin), dIface0Origin, 1e-9)
        << "interface distance to rotation origin must be preserved";

    // (b) Internal geometry (pairwise distances) must be invariant.
    EXPECT_NEAR(rotatefn_dist(moleculeList[0].tmpComCoord,
                    moleculeList[1].tmpComCoord),
        dPair, 1e-9)
        << "COM-COM distance must be preserved by a rigid rotation";
    EXPECT_NEAR(rotatefn_dist(moleculeList[0].tmpComCoord,
                    moleculeList[0].tmpICoords[1]),
        dComToIface, 1e-9)
        << "COM-interface distance must be preserved by a rigid rotation";

    // A 90 degree rotation about a non-degenerate axis must actually move the
    // points (guards against a silently broken quaternion).
    EXPECT_GT(rotatefn_dist(moleculeList[0].tmpComCoord, Coord(4.0, 2.0, 3.0)), 1e-6)
        << "a 90 degree rotation should visibly displace an off-axis point";

    std::cerr << "  Distances preserved; molecule 0 COM moved to ("
              << moleculeList[0].tmpComCoord.x << ", "
              << moleculeList[0].tmpComCoord.y << ", "
              << moleculeList[0].tmpComCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 5: Only molecules referenced by Complex::memberList are rotated.
// -----------------------------------------------------------------------------
void test_rotatefn_only_member_molecules_move()
{
    std::cerr << "\n[TEST] test_rotatefn_only_member_molecules_move\n"
              << "  Source file:   src/reactions/rotate.cpp\n"
              << "  Function:      rotate()\n"
              << "  Scenario:      moleculeList holds 3 molecules but the complex\n"
              << "                 memberList only contains indices {0, 2}.\n"
              << "  Pass criteria: molecules 0 and 2 are rotated, molecule 1 keeps\n"
              << "                 its original temporary coordinates.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rotatefn_make_molecule(Coord(1.0, 0.0, 0.0),
        { Coord(2.0, 0.0, 0.0) }, 0)); // member
    moleculeList.push_back(rotatefn_make_molecule(Coord(0.0, 5.0, 0.0),
        { Coord(0.0, 6.0, 0.0) }, 1)); // NOT a member
    moleculeList.push_back(rotatefn_make_molecule(Coord(0.0, 0.0, 4.0),
        { Coord(1.0, 0.0, 4.0) }, 0)); // member

    Complex targCom = rotatefn_make_complex({ 0, 2 });

    Coord rotOrigin(0.0, 0.0, 0.0);
    Quat rotZ180(0.0, 0.0, 0.0, 1.0);

    std::cerr << "  Calling rotate()...\n";
    rotate(rotOrigin, rotZ180, targCom, moleculeList);

    // Members: mirrored in x and y.
    rotatefn_expect_coord_near(moleculeList[0].tmpComCoord, Coord(-1.0, 0.0, 0.0),
        "member molecule 0 COM should be rotated");
    rotatefn_expect_coord_near(moleculeList[0].tmpICoords[0], Coord(-2.0, 0.0, 0.0),
        "member molecule 0 interface should be rotated");
    rotatefn_expect_coord_near(moleculeList[2].tmpComCoord, Coord(0.0, 0.0, 4.0),
        "member molecule 2 COM lies on the z axis so only x/y flip (both zero)");
    rotatefn_expect_coord_near(moleculeList[2].tmpICoords[0], Coord(-1.0, 0.0, 4.0),
        "member molecule 2 interface should be rotated");

    // Non-member must be completely untouched.
    rotatefn_expect_coord_near(moleculeList[1].tmpComCoord, Coord(0.0, 5.0, 0.0),
        "non-member molecule COM must NOT be rotated");
    rotatefn_expect_coord_near(moleculeList[1].tmpICoords[0], Coord(0.0, 6.0, 0.0),
        "non-member molecule interface must NOT be rotated");

    std::cerr << "  Non-member molecule remains at ("
              << moleculeList[1].tmpComCoord.x << ", "
              << moleculeList[1].tmpComCoord.y << ", "
              << moleculeList[1].tmpComCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 6: The permanent (non-temporary) coordinates must not be modified.
// -----------------------------------------------------------------------------
void test_rotatefn_permanent_coords_untouched()
{
    std::cerr << "\n[TEST] test_rotatefn_permanent_coords_untouched\n"
              << "  Source file:   src/reactions/rotate.cpp\n"
              << "  Function:      rotate()\n"
              << "  Scenario:      rotate a complex by pi about z and inspect the\n"
              << "                 permanent Molecule::comCoord / interfaceList.\n"
              << "  Pass criteria: permanent coordinates are unchanged while the\n"
              << "                 tmp* coordinates have changed.\n";

    const Coord origCom(2.0, 3.0, -1.0);
    const Coord origIface(2.0, 4.0, -1.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rotatefn_make_molecule(origCom, { origIface }, 0));

    Complex targCom = rotatefn_make_complex({ 0 });

    Coord rotOrigin(0.0, 0.0, 0.0);
    Quat rotZ180(0.0, 0.0, 0.0, 1.0);

    std::cerr << "  Calling rotate()...\n";
    rotate(rotOrigin, rotZ180, targCom, moleculeList);

    // Permanent copies must be pristine.
    rotatefn_expect_coord_near(moleculeList[0].comCoord, origCom,
        "permanent comCoord must not be modified by rotate()");
    rotatefn_expect_coord_near(moleculeList[0].interfaceList[0].coord, origIface,
        "permanent interfaceList coord must not be modified by rotate()");

    // Temporary copies must have moved.
    rotatefn_expect_coord_near(moleculeList[0].tmpComCoord,
        Coord(-origCom.x, -origCom.y, origCom.z),
        "tmpComCoord should hold the rotated position");
    rotatefn_expect_coord_near(moleculeList[0].tmpICoords[0],
        Coord(-origIface.x, -origIface.y, origIface.z),
        "tmpICoords should hold the rotated interface position");

    std::cerr << "  Permanent comCoord still (" << moleculeList[0].comCoord.x << ", "
              << moleculeList[0].comCoord.y << ", " << moleculeList[0].comCoord.z
              << "), tmpComCoord now (" << moleculeList[0].tmpComCoord.x << ", "
              << moleculeList[0].tmpComCoord.y << ", "
              << moleculeList[0].tmpComCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 7: Composition check. Two successive 90 degree rotations about z must
//         equal a single 180 degree rotation about z, whatever the handedness.
// -----------------------------------------------------------------------------
void test_rotatefn_two_90_equals_one_180()
{
    std::cerr << "\n[TEST] test_rotatefn_two_90_equals_one_180\n"
              << "  Source file:   src/reactions/rotate.cpp\n"
              << "  Function:      rotate()\n"
              << "  Scenario:      apply a 90 degree z rotation twice to one\n"
              << "                 complex, and a single 180 degree z rotation to\n"
              << "                 an identical copy.\n"
              << "  Pass criteria: both complexes end up at the same coordinates.\n";

    const Coord startCom(3.0, 1.0, -2.0);
    const Coord startIface(4.0, 1.5, -2.0);
    Coord rotOrigin(0.0, 0.0, 0.0);

    // --- Path A: two 90 degree rotations -------------------------------------
    std::vector<Molecule> listA;
    listA.push_back(rotatefn_make_molecule(startCom, { startIface }, 0));
    Complex comA = rotatefn_make_complex({ 0 });

    const double halfAng = M_PI / 4.0; // half of 90 degrees
    Quat rotZ90(std::cos(halfAng), 0.0, 0.0, std::sin(halfAng));

    std::cerr << "  Calling rotate() twice with the 90 degree quaternion...\n";
    rotate(rotOrigin, rotZ90, comA, listA);
    rotate(rotOrigin, rotZ90, comA, listA);

    // --- Path B: one 180 degree rotation ------------------------------------
    std::vector<Molecule> listB;
    listB.push_back(rotatefn_make_molecule(startCom, { startIface }, 0));
    Complex comB = rotatefn_make_complex({ 0 });

    Quat rotZ180(0.0, 0.0, 0.0, 1.0);

    std::cerr << "  Calling rotate() once with the 180 degree quaternion...\n";
    rotate(rotOrigin, rotZ180, comB, listB);

    // The two paths must agree (within numerical tolerance).
    EXPECT_NEAR(listA[0].tmpComCoord.x, listB[0].tmpComCoord.x, 1e-9)
        << "2 x 90 deg vs 1 x 180 deg: COM x must match";
    EXPECT_NEAR(listA[0].tmpComCoord.y, listB[0].tmpComCoord.y, 1e-9)
        << "2 x 90 deg vs 1 x 180 deg: COM y must match";
    EXPECT_NEAR(listA[0].tmpComCoord.z, listB[0].tmpComCoord.z, 1e-9)
        << "2 x 90 deg vs 1 x 180 deg: COM z must match";

    EXPECT_NEAR(listA[0].tmpICoords[0].x, listB[0].tmpICoords[0].x, 1e-9)
        << "2 x 90 deg vs 1 x 180 deg: interface x must match";
    EXPECT_NEAR(listA[0].tmpICoords[0].y, listB[0].tmpICoords[0].y, 1e-9)
        << "2 x 90 deg vs 1 x 180 deg: interface y must match";
    EXPECT_NEAR(listA[0].tmpICoords[0].z, listB[0].tmpICoords[0].z, 1e-9)
        << "2 x 90 deg vs 1 x 180 deg: interface z must match";

    std::cerr << "  Path A COM = (" << listA[0].tmpComCoord.x << ", "
              << listA[0].tmpComCoord.y << ", " << listA[0].tmpComCoord.z << ")\n"
              << "  Path B COM = (" << listB[0].tmpComCoord.x << ", "
              << listB[0].tmpComCoord.y << ", " << listB[0].tmpComCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 8: An empty memberList must be handled gracefully (no crash, no change).
// -----------------------------------------------------------------------------
void test_rotatefn_empty_member_list()
{
    std::cerr << "\n[TEST] test_rotatefn_empty_member_list\n"
              << "  Source file:   src/reactions/rotate.cpp\n"
              << "  Function:      rotate()\n"
              << "  Scenario:      a Complex with an empty memberList.\n"
              << "  Pass criteria: rotate() returns without touching any molecule.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(rotatefn_make_molecule(Coord(1.0, 2.0, 3.0),
        { Coord(1.0, 2.0, 4.0) }, 0));

    Complex targCom = rotatefn_make_complex({}); // deliberately empty

    Coord rotOrigin(0.0, 0.0, 0.0);
    Quat rotZ180(0.0, 0.0, 0.0, 1.0);

    std::cerr << "  Calling rotate() on an empty complex...\n";
    rotate(rotOrigin, rotZ180, targCom, moleculeList);

    // Nothing should have been modified since no member was listed.
    rotatefn_expect_coord_near(moleculeList[0].tmpComCoord, Coord(1.0, 2.0, 3.0),
        "molecule not listed in an empty memberList must not move");
    rotatefn_expect_coord_near(moleculeList[0].tmpICoords[0], Coord(1.0, 2.0, 4.0),
        "interface of an unlisted molecule must not move");

    std::cerr << "  No coordinates were modified, as expected.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper runs inside its own TEST so that a
// failure in one scenario does not prevent the remaining scenarios from running
// (all assertions above are non-fatal EXPECT_* checks).
// -----------------------------------------------------------------------------
TEST(RotateAssociation, IdentityQuaternion) { test_rotatefn_identity_quaternion(); }
TEST(RotateAssociation, PiRotationAboutZ) { test_rotatefn_180_degrees_about_z(); }
TEST(RotateAssociation, NonZeroOriginFixedPoint) { test_rotatefn_nonzero_origin_is_fixed_point(); }
TEST(RotateAssociation, PreservesRigidBodyGeometry) { test_rotatefn_preserves_rigid_body_geometry(); }
TEST(RotateAssociation, OnlyMemberMoleculesMove) { test_rotatefn_only_member_molecules_move(); }
TEST(RotateAssociation, PermanentCoordsUntouched) { test_rotatefn_permanent_coords_untouched(); }
TEST(RotateAssociation, TwoNinetyEqualsOneEighty) { test_rotatefn_two_90_equals_one_180(); }
TEST(RotateAssociation, EmptyMemberList) { test_rotatefn_empty_member_list(); }