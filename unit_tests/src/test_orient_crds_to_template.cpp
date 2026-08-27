/*! \file test_orient_crds_to_template.cpp
 *
 * ### Unit test for src/reactions/orient_crds_to_template.cpp
 *
 * Function under test:
 *
 *     Quat orient_crds_to_template(const MolTemplate& oneTemplate, Molecule& targMol)
 *
 * What the function does (read from the implementation):
 *
 *   1. It builds a first quaternion (`firstRot`) that rotates the vector
 *      (targMol.tmpICoords[0] - targMol.tmpComCoord) onto the template vector
 *      (oneTemplate.interfaceList[0].iCoord - oneTemplate.comCoord).  The sign of
 *      the rotation angle is verified/flipped by actually performing the rotation
 *      and comparing against the template coordinate.
 *   2. Every entry of targMol.tmpICoords is *replaced in place* by the rotated,
 *      **center-of-mass relative** vector (note: this is a side effect the test
 *      checks explicitly).
 *   3. If the molecule has more than one interface it searches the template for
 *      the first interface that is *not* collinear with interface 0.  If none is
 *      found and the template is not a rod, `firstRot` is returned immediately.
 *      Otherwise a second rotation about the (already aligned) interface-0 axis
 *      is computed and applied so that the remaining interfaces line up too.
 *   4. The composite quaternion `secondRot * firstRot` (or just `firstRot` when
 *      the molecule has a single tmp interface) is returned.
 *
 * Notes / deliberate omissions:
 *   - The second-rotation branch subtracts `tmpComCoord` a second time from the
 *     already COM-relative tmpICoords, so every test below puts the molecule's
 *     centre of mass (and the template's) at the origin, which is the regime in
 *     which the routine is used by association().
 *   - A perfect 180 degree misalignment of interface 0 is a mathematical
 *     degeneracy for this algorithm (the cross product is the zero vector, so no
 *     rotation axis exists).  That degenerate input is intentionally NOT tested
 *     here because the routine cannot recover from it.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers.  Everything is prefixed with "octt_" (orient_crds_to_template)
// so that the names cannot collide with the rest of the unit-test suite.
// -----------------------------------------------------------------------------
namespace {

//! Tolerance used when comparing rotated coordinates (pure trig round-off only).
constexpr double kOcttCrdTol = 1e-9;

/*! \brief Rotate a coordinate about an arbitrary axis by a given angle.
 *
 * Used only to *construct* test inputs, so that the expected answer of the
 * function under test is known analytically (it must undo this rotation).
 */
Coord octt_rotate_coord(const Coord& crd, const Coord& axis, double angle)
{
    Vector ax { axis };
    ax.normalize();
    const double sa { std::sin(angle / 2.0) };
    Quat rotQuat { std::cos(angle / 2.0), sa * ax.x, sa * ax.y, sa * ax.z };
    rotQuat = rotQuat.unit();

    Vector vec { crd };
    rotQuat.rotate(vec);
    return Coord { vec.x, vec.y, vec.z };
}

/*! \brief Apply a quaternion to a coordinate (treated as a vector from origin). */
Coord octt_apply_quat(Quat rotQuat, const Coord& crd)
{
    Vector vec { crd };
    rotQuat.rotate(vec);
    return Coord { vec.x, vec.y, vec.z };
}

/*! \brief Build a MolTemplate with its COM at the origin and the given interfaces. */
MolTemplate octt_make_template(const std::vector<Coord>& ifaceCoords, bool isRod = false)
{
    MolTemplate oneTemplate {};
    oneTemplate.molName = "octtMol";
    oneTemplate.molTypeIndex = 0;
    oneTemplate.comCoord = Coord { 0.0, 0.0, 0.0 };
    oneTemplate.isRod = isRod;
    oneTemplate.radius = 1.0;
    oneTemplate.mass = 1.0;

    for (unsigned ifaceItr { 0 }; ifaceItr < ifaceCoords.size(); ++ifaceItr) {
        Interface iface {};
        iface.index = static_cast<int>(ifaceItr);
        iface.name = "i" + std::to_string(ifaceItr);
        iface.iCoord = ifaceCoords[ifaceItr];
        // one (default) state per interface, exactly as the parser would create
        iface.stateList.emplace_back(iface.name, static_cast<int>(ifaceItr));
        oneTemplate.interfaceList.emplace_back(iface);
    }
    return oneTemplate;
}

/*! \brief Build a Molecule with its COM at the origin and the given interfaces.
 *
 * set_tmp_association_coords() is called so that tmpComCoord/tmpICoords are
 * populated exactly the way association() populates them before calling the
 * function under test.
 */
Molecule octt_make_molecule(const std::vector<Coord>& ifaceCoords)
{
    Molecule targMol {};
    targMol.index = 0;
    targMol.myComIndex = 0;
    targMol.molTypeIndex = 0;
    targMol.mass = 1.0;
    targMol.comCoord = Coord { 0.0, 0.0, 0.0 };

    for (unsigned ifaceItr { 0 }; ifaceItr < ifaceCoords.size(); ++ifaceItr) {
        Molecule::Iface iface {};
        iface.coord = ifaceCoords[ifaceItr];
        iface.relIndex = static_cast<int>(ifaceItr);
        iface.index = static_cast<int>(ifaceItr);
        iface.molTypeIndex = 0;
        targMol.interfaceList.emplace_back(iface);
    }

    targMol.set_tmp_association_coords(); // tmpComCoord = comCoord, tmpICoords = iface coords
    return targMol;
}

/*! \brief EXPECT that two coordinates agree component-wise within tol. */
void octt_expect_coord_near(const Coord& got, const Coord& expected, const std::string& what)
{
    EXPECT_NEAR(got.x, expected.x, kOcttCrdTol) << what << " -- x component";
    EXPECT_NEAR(got.y, expected.y, kOcttCrdTol) << what << " -- y component";
    EXPECT_NEAR(got.z, expected.z, kOcttCrdTol) << what << " -- z component";
}

/*! \brief Magnitude of a coordinate treated as a vector from the origin. */
double octt_mag(const Coord& crd)
{
    return std::sqrt(crd.x * crd.x + crd.y * crd.y + crd.z * crd.z);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: single-interface molecule.
//         The routine only has to perform the first rotation, and the sign of
//         the rotation angle must be corrected internally (the naive +angle
//         rotation sends +y to -x, not to +x).
// -----------------------------------------------------------------------------
void test_octt_single_interface_alignment()
{
    std::cerr << "\n[TEST] test_octt_single_interface_alignment\n"
              << "  Source file:   src/reactions/orient_crds_to_template.cpp\n"
              << "  Function:      orient_crds_to_template()\n"
              << "  Scenario:      template interface is at (1,0,0); the molecule's\n"
              << "                 tmp interface is at (0,1,0) (a 90 deg mismatch).\n"
              << "  Pass criteria: tmpICoords[0] becomes (1,0,0), the returned quat\n"
              << "                 is a unit quaternion, and applying it to the\n"
              << "                 original interface vector also yields (1,0,0).\n";

    MolTemplate oneTemplate { octt_make_template({ Coord { 1.0, 0.0, 0.0 } }) };
    Molecule targMol { octt_make_molecule({ Coord { 0.0, 1.0, 0.0 } }) };

    // keep a copy of the input, the function rewrites tmpICoords in place
    const Coord origIface { targMol.tmpICoords[0] };

    Quat rotQuat { orient_crds_to_template(oneTemplate, targMol) };
    std::cerr << "  Returned quaternion: " << rotQuat << '\n';
    std::cerr << "  tmpICoords[0] after call: " << targMol.tmpICoords[0] << '\n';

    // 1. the tmp coordinate has been rotated onto the template coordinate
    octt_expect_coord_near(targMol.tmpICoords[0], Coord { 1.0, 0.0, 0.0 },
        "rotated tmp interface 0 should equal the template interface");

    // 2. the returned quaternion is a rotation (unit) quaternion
    EXPECT_NEAR(rotQuat.mag(), 1.0, 1e-12)
        << "the returned quaternion must have unit magnitude";

    // 3. the returned quaternion reproduces the same rotation on the raw input
    octt_expect_coord_near(octt_apply_quat(rotQuat, origIface), Coord { 1.0, 0.0, 0.0 },
        "returned quat applied to the original interface vector");

    // 4. the function must not disturb the molecule's real (non-tmp) coordinates
    EXPECT_DOUBLE_EQ(targMol.comCoord.x, 0.0) << "comCoord must be untouched";
    EXPECT_DOUBLE_EQ(targMol.comCoord.y, 0.0) << "comCoord must be untouched";
    EXPECT_DOUBLE_EQ(targMol.comCoord.z, 0.0) << "comCoord must be untouched";
    EXPECT_DOUBLE_EQ(targMol.interfaceList[0].coord.y, 1.0)
        << "interfaceList coordinates must be untouched";
    EXPECT_DOUBLE_EQ(targMol.tmpComCoord.x, 0.0) << "tmpComCoord must be untouched";
}

// -----------------------------------------------------------------------------
// Test 2: two non-collinear interfaces, molecule is an exactly rotated copy of
//         the template.  Both rotations are exercised and the routine must
//         recover the inverse of the applied rotation.
// -----------------------------------------------------------------------------
void test_octt_two_interfaces_recovers_known_rotation()
{
    std::cerr << "\n[TEST] test_octt_two_interfaces_recovers_known_rotation\n"
              << "  Source file:   src/reactions/orient_crds_to_template.cpp\n"
              << "  Function:      orient_crds_to_template()\n"
              << "  Scenario:      the molecule's tmp coordinates are the template\n"
              << "                 coordinates rotated by 0.7 rad about (1,2,3).\n"
              << "  Pass criteria: after the call every tmpICoord equals the matching\n"
              << "                 template iCoord (i.e. the rotation was undone), and\n"
              << "                 the returned quat performs the same mapping.\n";

    const std::vector<Coord> tmplCrds { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 2.0, 0.0 } };
    const Coord axis { 1.0, 2.0, 3.0 };
    const double angle { 0.7 };

    std::vector<Coord> molCrds;
    for (const auto& crd : tmplCrds)
        molCrds.push_back(octt_rotate_coord(crd, axis, angle));

    MolTemplate oneTemplate { octt_make_template(tmplCrds) };
    Molecule targMol { octt_make_molecule(molCrds) };

    const std::vector<Coord> origCrds { targMol.tmpICoords };

    Quat rotQuat { orient_crds_to_template(oneTemplate, targMol) };
    std::cerr << "  Returned quaternion: " << rotQuat << '\n';

    for (unsigned i { 0 }; i < tmplCrds.size(); ++i) {
        std::cerr << "  tmpICoords[" << i << "] after call: " << targMol.tmpICoords[i]
                  << "  (expected " << tmplCrds[i] << ")\n";
        octt_expect_coord_near(targMol.tmpICoords[i], tmplCrds[i],
            "tmp interface " + std::to_string(i) + " should be oriented to the template");

        // the returned composite quaternion must perform the very same mapping
        octt_expect_coord_near(octt_apply_quat(rotQuat, origCrds[i]), tmplCrds[i],
            "returned quat applied to original interface " + std::to_string(i));
    }

    EXPECT_NEAR(rotQuat.mag(), 1.0, 1e-12)
        << "the composite quaternion must still be a unit quaternion";
}

// -----------------------------------------------------------------------------
// Test 3: molecule already in template orientation -> identity quaternion and
//         unchanged coordinates.
// -----------------------------------------------------------------------------
void test_octt_identity_when_already_oriented()
{
    std::cerr << "\n[TEST] test_octt_identity_when_already_oriented\n"
              << "  Source file:   src/reactions/orient_crds_to_template.cpp\n"
              << "  Function:      orient_crds_to_template()\n"
              << "  Scenario:      the molecule tmp coordinates are identical to the\n"
              << "                 template coordinates (nothing to do).\n"
              << "  Pass criteria: returned quaternion is the identity (w=+/-1, xyz=0)\n"
              << "                 and the tmp coordinates are unchanged.\n";

    const std::vector<Coord> crds { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } };
    MolTemplate oneTemplate { octt_make_template(crds) };
    Molecule targMol { octt_make_molecule(crds) };

    Quat rotQuat { orient_crds_to_template(oneTemplate, targMol) };
    std::cerr << "  Returned quaternion: " << rotQuat << '\n';

    // identity rotation: the vector part is zero and |w| == 1
    EXPECT_NEAR(std::abs(rotQuat.w), 1.0, 1e-12) << "identity quaternion should have |w| == 1";
    EXPECT_NEAR(rotQuat.x, 0.0, 1e-12) << "identity quaternion should have x == 0";
    EXPECT_NEAR(rotQuat.y, 0.0, 1e-12) << "identity quaternion should have y == 0";
    EXPECT_NEAR(rotQuat.z, 0.0, 1e-12) << "identity quaternion should have z == 0";

    for (unsigned i { 0 }; i < crds.size(); ++i) {
        octt_expect_coord_near(targMol.tmpICoords[i], crds[i],
            "tmp interface " + std::to_string(i) + " must be unchanged");
    }
}

// -----------------------------------------------------------------------------
// Test 4: every template interface is collinear with interface 0 and the
//         template is NOT flagged as a rod -> the routine returns after the
//         first rotation only.  For a collinear (1D) molecule that single
//         rotation is nevertheless sufficient to align all interfaces.
// -----------------------------------------------------------------------------
void test_octt_all_collinear_returns_first_rotation()
{
    std::cerr << "\n[TEST] test_octt_all_collinear_returns_first_rotation\n"
              << "  Source file:   src/reactions/orient_crds_to_template.cpp\n"
              << "  Function:      orient_crds_to_template()\n"
              << "  Scenario:      template interfaces (0,0,1) and (0,0,2) are\n"
              << "                 collinear, isRod == false, so the early-return\n"
              << "                 'allcollinear' branch is taken.\n"
              << "  Pass criteria: both tmp interfaces are aligned with the template\n"
              << "                 by the first rotation alone.\n";

    const std::vector<Coord> tmplCrds { Coord { 0.0, 0.0, 1.0 }, Coord { 0.0, 0.0, 2.0 } };
    const Coord axis { 0.0, 1.0, 0.0 }; // rotate about y: +z -> +x for angle = pi/2
    const double angle { M_PI / 2.0 };

    std::vector<Coord> molCrds;
    for (const auto& crd : tmplCrds)
        molCrds.push_back(octt_rotate_coord(crd, axis, angle));

    MolTemplate oneTemplate { octt_make_template(tmplCrds, /*isRod=*/false) };
    Molecule targMol { octt_make_molecule(molCrds) };

    const std::vector<Coord> origCrds { targMol.tmpICoords };
    std::cerr << "  Molecule tmp interface 0 before call: " << origCrds[0] << '\n';

    Quat rotQuat { orient_crds_to_template(oneTemplate, targMol) };
    std::cerr << "  Returned quaternion (should be firstRot only): " << rotQuat << '\n';

    for (unsigned i { 0 }; i < tmplCrds.size(); ++i) {
        octt_expect_coord_near(targMol.tmpICoords[i], tmplCrds[i],
            "collinear tmp interface " + std::to_string(i) + " aligned by firstRot");
        octt_expect_coord_near(octt_apply_quat(rotQuat, origCrds[i]), tmplCrds[i],
            "returned quat applied to original collinear interface " + std::to_string(i));
    }

    EXPECT_NEAR(rotQuat.mag(), 1.0, 1e-12) << "firstRot must be a unit quaternion";
}

// -----------------------------------------------------------------------------
// Test 5: same collinear geometry but with isRod == true, which forces the
//         routine into the second-rotation branch even though the projections
//         degenerate to zero-length vectors (the second rotation must come out
//         as the identity so the result is unchanged).
//         NOTE: this path legitimately prints a "vectors with magnitude 0"
//         warning from Vector::dot_theta() to stdout.
// -----------------------------------------------------------------------------
void test_octt_rod_collinear_second_rotation_is_identity()
{
    std::cerr << "\n[TEST] test_octt_rod_collinear_second_rotation_is_identity\n"
              << "  Source file:   src/reactions/orient_crds_to_template.cpp\n"
              << "  Function:      orient_crds_to_template()\n"
              << "  Scenario:      collinear interfaces with isRod == true, so the\n"
              << "                 second-rotation branch runs on degenerate (zero)\n"
              << "                 projections.\n"
              << "  Pass criteria: the second rotation contributes nothing, so both\n"
              << "                 tmp interfaces still land on the template.\n"
              << "                 (a dot_theta magnitude warning on stdout is expected)\n";

    const std::vector<Coord> tmplCrds { Coord { 0.0, 0.0, 1.0 }, Coord { 0.0, 0.0, 2.0 } };
    const Coord axis { 1.0, 0.0, 0.0 };
    const double angle { 0.9 };

    std::vector<Coord> molCrds;
    for (const auto& crd : tmplCrds)
        molCrds.push_back(octt_rotate_coord(crd, axis, angle));

    MolTemplate oneTemplate { octt_make_template(tmplCrds, /*isRod=*/true) };
    Molecule targMol { octt_make_molecule(molCrds) };

    const std::vector<Coord> origCrds { targMol.tmpICoords };

    Quat rotQuat { orient_crds_to_template(oneTemplate, targMol) };
    std::cerr << "  Returned quaternion: " << rotQuat << '\n';

    for (unsigned i { 0 }; i < tmplCrds.size(); ++i) {
        octt_expect_coord_near(targMol.tmpICoords[i], tmplCrds[i],
            "rod tmp interface " + std::to_string(i) + " should still be aligned");
        octt_expect_coord_near(octt_apply_quat(rotQuat, origCrds[i]), tmplCrds[i],
            "returned quat applied to original rod interface " + std::to_string(i));
    }

    EXPECT_NEAR(rotQuat.mag(), 1.0, 1e-12)
        << "secondRot * firstRot must still be a unit quaternion";
}

// -----------------------------------------------------------------------------
// Test 6: three interfaces where interface 1 is collinear with interface 0.
//         This exercises the search loop that skips collinear interfaces and
//         settles on interface 2 for the second rotation.
// -----------------------------------------------------------------------------
void test_octt_skips_collinear_interface_in_search()
{
    std::cerr << "\n[TEST] test_octt_skips_collinear_interface_in_search\n"
              << "  Source file:   src/reactions/orient_crds_to_template.cpp\n"
              << "  Function:      orient_crds_to_template()\n"
              << "  Scenario:      template interfaces (1,0,0), (2,0,0) [collinear with\n"
              << "                 the first] and (0,0,1); the search loop must skip\n"
              << "                 index 1 and use index 2 for the second rotation.\n"
              << "  Pass criteria: all three tmp interfaces are oriented onto the\n"
              << "                 template after a known rotation is applied.\n";

    const std::vector<Coord> tmplCrds { Coord { 1.0, 0.0, 0.0 }, Coord { 2.0, 0.0, 0.0 },
        Coord { 0.0, 0.0, 1.0 } };
    const Coord axis { -1.0, 3.0, 2.0 };
    const double angle { 1.3 };

    std::vector<Coord> molCrds;
    for (const auto& crd : tmplCrds)
        molCrds.push_back(octt_rotate_coord(crd, axis, angle));

    MolTemplate oneTemplate { octt_make_template(tmplCrds) };
    Molecule targMol { octt_make_molecule(molCrds) };

    const std::vector<Coord> origCrds { targMol.tmpICoords };

    Quat rotQuat { orient_crds_to_template(oneTemplate, targMol) };
    std::cerr << "  Returned quaternion: " << rotQuat << '\n';

    for (unsigned i { 0 }; i < tmplCrds.size(); ++i) {
        std::cerr << "  tmpICoords[" << i << "] after call: " << targMol.tmpICoords[i]
                  << "  (expected " << tmplCrds[i] << ")\n";
        octt_expect_coord_near(targMol.tmpICoords[i], tmplCrds[i],
            "tmp interface " + std::to_string(i) + " should be oriented to the template");
        octt_expect_coord_near(octt_apply_quat(rotQuat, origCrds[i]), tmplCrds[i],
            "returned quat applied to original interface " + std::to_string(i));
    }
}

// -----------------------------------------------------------------------------
// Test 7: rigid-body invariants.  Independently of which branch is taken, the
//         routine may only rotate: interface lengths and the angles between
//         interfaces must be preserved, and the tmp COM must not move.
// -----------------------------------------------------------------------------
void test_octt_preserves_rigid_body_geometry()
{
    std::cerr << "\n[TEST] test_octt_preserves_rigid_body_geometry\n"
              << "  Source file:   src/reactions/orient_crds_to_template.cpp\n"
              << "  Function:      orient_crds_to_template()\n"
              << "  Scenario:      a three-interface molecule in an arbitrary\n"
              << "                 orientation is oriented to its template.\n"
              << "  Pass criteria: |iface| and the pairwise dot products between the\n"
              << "                 interface vectors are unchanged (pure rotation), and\n"
              << "                 tmpComCoord is untouched.\n";

    const std::vector<Coord> tmplCrds { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 2.0, 0.0 },
        Coord { 0.0, 0.0, 3.0 } };
    const Coord axis { 2.0, -1.0, 0.5 };
    const double angle { -0.45 };

    std::vector<Coord> molCrds;
    for (const auto& crd : tmplCrds)
        molCrds.push_back(octt_rotate_coord(crd, axis, angle));

    MolTemplate oneTemplate { octt_make_template(tmplCrds) };
    Molecule targMol { octt_make_molecule(molCrds) };

    const std::vector<Coord> origCrds { targMol.tmpICoords };
    const Coord origTmpCom { targMol.tmpComCoord };

    orient_crds_to_template(oneTemplate, targMol);

    // number of interfaces is unchanged
    ASSERT_EQ(targMol.tmpICoords.size(), origCrds.size())
        << "the routine must not add or remove tmp interfaces";

    // lengths conserved
    for (unsigned i { 0 }; i < origCrds.size(); ++i) {
        const double before { octt_mag(origCrds[i]) };
        const double after { octt_mag(targMol.tmpICoords[i]) };
        std::cerr << "  |iface " << i << "|  before = " << before << ", after = " << after << '\n';
        EXPECT_NEAR(after, before, kOcttCrdTol)
            << "interface " << i << " magnitude must be conserved by a pure rotation";
    }

    // pairwise dot products conserved (i.e. internal angles unchanged)
    for (unsigned i { 0 }; i < origCrds.size(); ++i) {
        for (unsigned j { i + 1 }; j < origCrds.size(); ++j) {
            const double dotBefore { origCrds[i].x * origCrds[j].x + origCrds[i].y * origCrds[j].y
                + origCrds[i].z * origCrds[j].z };
            const double dotAfter { targMol.tmpICoords[i].x * targMol.tmpICoords[j].x
                + targMol.tmpICoords[i].y * targMol.tmpICoords[j].y
                + targMol.tmpICoords[i].z * targMol.tmpICoords[j].z };
            EXPECT_NEAR(dotAfter, dotBefore, kOcttCrdTol)
                << "dot product of interfaces " << i << " and " << j << " must be conserved";
        }
    }

    // the centre of mass is never touched by this routine
    EXPECT_DOUBLE_EQ(targMol.tmpComCoord.x, origTmpCom.x) << "tmpComCoord.x must be unchanged";
    EXPECT_DOUBLE_EQ(targMol.tmpComCoord.y, origTmpCom.y) << "tmpComCoord.y must be unchanged";
    EXPECT_DOUBLE_EQ(targMol.tmpComCoord.z, origTmpCom.z) << "tmpComCoord.z must be unchanged";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario lives in its own TEST so that a failure in
// one does not stop the others from running.
// -----------------------------------------------------------------------------
TEST(OrientCrdsToTemplate, SingleInterfaceAlignment) { test_octt_single_interface_alignment(); }
TEST(OrientCrdsToTemplate, TwoInterfacesRecoversKnownRotation) { test_octt_two_interfaces_recovers_known_rotation(); }
TEST(OrientCrdsToTemplate, IdentityWhenAlreadyOriented) { test_octt_identity_when_already_oriented(); }
TEST(OrientCrdsToTemplate, AllCollinearReturnsFirstRotation) { test_octt_all_collinear_returns_first_rotation(); }
TEST(OrientCrdsToTemplate, RodCollinearSecondRotationIsIdentity) { test_octt_rod_collinear_second_rotation_is_identity(); }
TEST(OrientCrdsToTemplate, SkipsCollinearInterfaceInSearch) { test_octt_skips_collinear_interface_in_search(); }
TEST(OrientCrdsToTemplate, PreservesRigidBodyGeometry) { test_octt_preserves_rigid_body_geometry(); }