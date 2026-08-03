/*! \file test_determine_normal.cpp
 *
 * ### Unit test for src/reactions/determine_normal.cpp
 *
 * Function under test:
 *
 *     Vector determine_normal(Vector normal, const MolTemplate& molTemplate,
 *                             Molecule oneMol);
 *
 * What the function does (from the source):
 *   1. If the Molecule has no interfaces (it is a "point"), the zero vector
 *      {0,0,0} is returned immediately.
 *   2. Otherwise the Molecule's temporary association coordinates
 *      (tmpComCoord / tmpICoords) are re-centered on the COM.
 *   3. orient_crds_to_template() produces the quaternion that rotates the
 *      Molecule's real coordinates onto its MolTemplate (internal) coordinates.
 *   4. That quaternion is normalized, INVERTED, and used to rotate the supplied
 *      `normal` (which is expressed in the template/internal frame) into the
 *      Molecule's real frame.
 *   5. A sanity check reverses the rotation on every interface; if any
 *      interface does not come back to where it started the zero vector is
 *      returned.
 *   6. The resulting normal is normalized and returned.
 *
 * Testing strategy (all output is echoed to stderr so a reader can follow
 * exactly which behaviour is being probed and what the pass criterion is):
 *   - point molecule            -> zero vector
 *   - already-aligned molecule  -> normal returned unchanged (identity rotation)
 *   - un-normalized input       -> output has unit magnitude
 *   - molecule rotated 90 deg about z, normal along z -> z axis is invariant
 *   - molecule rotated 90 deg about x, normal along z -> normal lands on y axis
 *   - non-zero COM offset       -> handled by the re-centering step
 *   - the Molecule argument is taken BY VALUE, so the caller's copy must be
 *     untouched
 *   - single-interface (rod) molecule -> result must be either unit length or
 *     the documented zero-vector failure value
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace so nothing collides with the rest of the
// generated test suite).
// -----------------------------------------------------------------------------
namespace {

//! Tolerance used for all floating point comparisons in this file.
constexpr double kDnTol = 1e-6;

/*! \brief Build a MolTemplate whose internal interface coordinates are given.
 *
 * \param[in] ifaceCoords internal (COM-relative) interface coordinates
 * \param[in] isRod       flag telling the orientation code that the molecule is
 *                        strictly one dimensional
 */
MolTemplate dn_make_template(const std::vector<Coord>& ifaceCoords, bool isRod)
{
    MolTemplate temp;
    temp.molName = "dnTestMol";
    temp.molTypeIndex = 0;
    temp.comCoord = Coord(0.0, 0.0, 0.0); // templates are always centered
    temp.isRod = isRod;
    temp.isPoint = ifaceCoords.empty();
    temp.radius = 1.0;
    temp.interfaceList.clear();

    for (std::size_t i { 0 }; i < ifaceCoords.size(); ++i) {
        Interface iface;
        iface.index = static_cast<int>(i);
        iface.name = "iface" + std::to_string(i);
        iface.iCoord = ifaceCoords[i];
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Build a Molecule with absolute interface coordinates.
 *
 * determine_normal() only looks at interfaceList (for the emptiness check) and
 * at the temporary association coordinates (tmpComCoord / tmpICoords), so both
 * sets are filled in consistently here - exactly as
 * Molecule::set_tmp_association_coords() would do.
 *
 * \param[in] com         absolute center-of-mass coordinate
 * \param[in] ifaceCoords absolute interface coordinates
 */
Molecule dn_make_molecule(const Coord& com, const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.index = 0;
    mol.myComIndex = 0;
    mol.molTypeIndex = 0;
    mol.comCoord = com;

    mol.interfaceList.clear();
    for (std::size_t i { 0 }; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface;
        iface.coord = ifaceCoords[i];
        iface.index = static_cast<int>(i);
        iface.relIndex = static_cast<int>(i);
        iface.molTypeIndex = 0;
        mol.interfaceList.push_back(iface);
    }

    // Temporary association coordinates -- these are what determine_normal uses.
    mol.tmpComCoord = com;
    mol.tmpICoords.clear();
    for (const auto& crd : ifaceCoords)
        mol.tmpICoords.push_back(crd);

    return mol;
}

//! Rotate a coordinate by +90 degrees about the z axis: (x,y,z) -> (-y,x,z).
Coord dn_rot90_z(const Coord& c) { return Coord(-c.y, c.x, c.z); }

//! Rotate a coordinate by +90 degrees about the x axis: (x,y,z) -> (x,-z,y).
Coord dn_rot90_x(const Coord& c) { return Coord(c.x, -c.z, c.y); }

//! Pretty-print a Vector to stderr.
void dn_report(const std::string& label, Vector v)
{
    std::cerr << "    " << label << " = (" << v.x << ", " << v.y << ", " << v.z
              << "), |v| = " << v.get_magnitude() << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a Molecule with no interfaces is a "point" and has no normal.
// -----------------------------------------------------------------------------
void test_dn_point_molecule_returns_zero()
{
    std::cerr << "\n[TEST] test_dn_point_molecule_returns_zero\n"
              << "  Source file:   src/reactions/determine_normal.cpp\n"
              << "  Function:      determine_normal()\n"
              << "  Scenario:      Molecule::interfaceList is empty (a point molecule).\n"
              << "  Pass criteria: the early-return branch fires and (0,0,0) comes back.\n";

    // Empty template and empty molecule => the first branch of the function.
    MolTemplate molTemplate { dn_make_template({}, /*isRod=*/false) };
    Molecule oneMol { dn_make_molecule(Coord(0.0, 0.0, 0.0), {}) };

    Vector normal(0.0, 0.0, 1.0);
    Vector result { determine_normal(normal, molTemplate, oneMol) };
    dn_report("returned normal", result);

    EXPECT_DOUBLE_EQ(result.x, 0.0) << "point molecule must return x == 0";
    EXPECT_DOUBLE_EQ(result.y, 0.0) << "point molecule must return y == 0";
    EXPECT_DOUBLE_EQ(result.z, 0.0) << "point molecule must return z == 0";
}

// -----------------------------------------------------------------------------
// Test 2: a Molecule already in its template orientation => identity rotation,
//         so the supplied normal must come straight back out.
// -----------------------------------------------------------------------------
void test_dn_aligned_molecule_returns_input_normal()
{
    std::cerr << "\n[TEST] test_dn_aligned_molecule_returns_input_normal\n"
              << "  Source file:   src/reactions/determine_normal.cpp\n"
              << "  Function:      determine_normal()\n"
              << "  Scenario:      molecule coordinates == template coordinates,\n"
              << "                 COM at the origin, input normal = (0,0,1).\n"
              << "  Pass criteria: rotation is the identity, so (0,0,1) is returned.\n";

    // Two non-collinear interfaces so the orientation is fully determined.
    std::vector<Coord> ifaces { Coord(1.0, 0.0, 0.0), Coord(0.0, 1.0, 0.0) };
    MolTemplate molTemplate { dn_make_template(ifaces, /*isRod=*/false) };
    Molecule oneMol { dn_make_molecule(Coord(0.0, 0.0, 0.0), ifaces) };

    Vector normal(0.0, 0.0, 1.0);
    Vector result { determine_normal(normal, molTemplate, oneMol) };
    dn_report("returned normal", result);

    EXPECT_NEAR(result.x, 0.0, kDnTol) << "identity rotation must leave x == 0";
    EXPECT_NEAR(result.y, 0.0, kDnTol) << "identity rotation must leave y == 0";
    EXPECT_NEAR(result.z, 1.0, kDnTol) << "identity rotation must leave z == 1";
    EXPECT_NEAR(result.get_magnitude(), 1.0, kDnTol) << "returned normal must be a unit vector";
}

// -----------------------------------------------------------------------------
// Test 3: the function ends with normal.normalize(), so an over-long input
//         normal must come back with unit magnitude.
// -----------------------------------------------------------------------------
void test_dn_output_is_normalized()
{
    std::cerr << "\n[TEST] test_dn_output_is_normalized\n"
              << "  Source file:   src/reactions/determine_normal.cpp\n"
              << "  Function:      determine_normal()\n"
              << "  Scenario:      aligned molecule, input normal = (0,0,5)\n"
              << "                 (deliberately NOT unit length).\n"
              << "  Pass criteria: output is (0,0,1) -- normalize() was applied.\n";

    std::vector<Coord> ifaces { Coord(1.0, 0.0, 0.0), Coord(0.0, 1.0, 0.0) };
    MolTemplate molTemplate { dn_make_template(ifaces, /*isRod=*/false) };
    Molecule oneMol { dn_make_molecule(Coord(0.0, 0.0, 0.0), ifaces) };

    Vector normal(0.0, 0.0, 5.0);
    Vector result { determine_normal(normal, molTemplate, oneMol) };
    dn_report("returned normal", result);

    EXPECT_NEAR(result.get_magnitude(), 1.0, kDnTol)
        << "an un-normalized input normal must be normalized on output";
    EXPECT_NEAR(result.z, 1.0, kDnTol) << "direction must be preserved (+z)";
    EXPECT_NEAR(result.x, 0.0, kDnTol) << "x component must stay 0";
    EXPECT_NEAR(result.y, 0.0, kDnTol) << "y component must stay 0";
}

// -----------------------------------------------------------------------------
// Test 4: molecule rotated 90 degrees about z. The z axis is invariant under
//         such a rotation (and under its inverse), so a z-normal is unchanged.
//         This makes the test insensitive to the rotation-sign convention.
// -----------------------------------------------------------------------------
void test_dn_rotation_about_z_leaves_z_normal_invariant()
{
    std::cerr << "\n[TEST] test_dn_rotation_about_z_leaves_z_normal_invariant\n"
              << "  Source file:   src/reactions/determine_normal.cpp\n"
              << "  Function:      determine_normal()\n"
              << "  Scenario:      molecule coordinates are the template rotated\n"
              << "                 +90 deg about z; internal normal = (0,0,1).\n"
              << "  Pass criteria: the z axis is the rotation axis, therefore the\n"
              << "                 returned normal is still (0,0,1).\n";

    std::vector<Coord> templateIfaces { Coord(1.0, 0.0, 0.0), Coord(0.0, 1.0, 0.0) };
    std::vector<Coord> molIfaces { dn_rot90_z(templateIfaces[0]), dn_rot90_z(templateIfaces[1]) };

    std::cerr << "    template ifaces: (1,0,0) (0,1,0)\n"
              << "    molecule ifaces: (" << molIfaces[0].x << ',' << molIfaces[0].y << ','
              << molIfaces[0].z << ") (" << molIfaces[1].x << ',' << molIfaces[1].y << ','
              << molIfaces[1].z << ")\n";

    MolTemplate molTemplate { dn_make_template(templateIfaces, /*isRod=*/false) };
    Molecule oneMol { dn_make_molecule(Coord(0.0, 0.0, 0.0), molIfaces) };

    Vector normal(0.0, 0.0, 1.0);
    Vector result { determine_normal(normal, molTemplate, oneMol) };
    dn_report("returned normal", result);

    EXPECT_NEAR(result.get_magnitude(), 1.0, kDnTol) << "result must be a unit vector";
    EXPECT_NEAR(result.z, 1.0, kDnTol) << "rotation axis component must be untouched";
    EXPECT_NEAR(result.x, 0.0, kDnTol) << "x must remain 0 for a pure z rotation";
    EXPECT_NEAR(result.y, 0.0, kDnTol) << "y must remain 0 for a pure z rotation";
}

// -----------------------------------------------------------------------------
// Test 5: molecule rotated 90 degrees about x. Now the z-normal of the template
//         frame must be carried onto the +/- y axis of the molecule frame.
//         Only the magnitude of each component is asserted, because the overall
//         sign depends on the quaternion convention inside
//         orient_crds_to_template(); the actual sign is printed for the reader.
// -----------------------------------------------------------------------------
void test_dn_rotation_about_x_moves_normal_to_y_axis()
{
    std::cerr << "\n[TEST] test_dn_rotation_about_x_moves_normal_to_y_axis\n"
              << "  Source file:   src/reactions/determine_normal.cpp\n"
              << "  Function:      determine_normal()\n"
              << "  Scenario:      molecule coordinates are the template rotated\n"
              << "                 +90 deg about x; internal normal = (0,0,1).\n"
              << "  Pass criteria: the returned normal lies along the y axis\n"
              << "                 (|y| == 1, x == z == 0). The sign is reported\n"
              << "                 but not asserted (convention dependent).\n";

    std::vector<Coord> templateIfaces { Coord(1.0, 0.0, 0.0), Coord(0.0, 1.0, 0.0) };
    std::vector<Coord> molIfaces { dn_rot90_x(templateIfaces[0]), dn_rot90_x(templateIfaces[1]) };

    std::cerr << "    molecule ifaces: (" << molIfaces[0].x << ',' << molIfaces[0].y << ','
              << molIfaces[0].z << ") (" << molIfaces[1].x << ',' << molIfaces[1].y << ','
              << molIfaces[1].z << ")\n";

    MolTemplate molTemplate { dn_make_template(templateIfaces, /*isRod=*/false) };
    Molecule oneMol { dn_make_molecule(Coord(0.0, 0.0, 0.0), molIfaces) };

    Vector normal(0.0, 0.0, 1.0);
    Vector result { determine_normal(normal, molTemplate, oneMol) };
    dn_report("returned normal", result);
    std::cerr << "    (sign of the y component reported only: "
              << (result.y >= 0.0 ? "+y" : "-y") << ")\n";

    EXPECT_NEAR(result.get_magnitude(), 1.0, kDnTol) << "result must be a unit vector";
    EXPECT_NEAR(std::abs(result.y), 1.0, kDnTol)
        << "the internal z-normal must be carried onto the molecule's y axis";
    EXPECT_NEAR(result.x, 0.0, kDnTol) << "x must remain 0 for a pure x rotation";
    EXPECT_NEAR(result.z, 0.0, kDnTol) << "z must vanish after a 90 deg x rotation";
}

// -----------------------------------------------------------------------------
// Test 6: the first thing the function does is subtract tmpComCoord from every
//         interface, so a translated molecule must give the same answer as the
//         centered one.
// -----------------------------------------------------------------------------
void test_dn_handles_nonzero_com_offset()
{
    std::cerr << "\n[TEST] test_dn_handles_nonzero_com_offset\n"
              << "  Source file:   src/reactions/determine_normal.cpp\n"
              << "  Function:      determine_normal()\n"
              << "  Scenario:      same aligned geometry as the identity test, but\n"
              << "                 the whole molecule is translated to (10,-5,3).\n"
              << "  Pass criteria: re-centering makes the answer identical to the\n"
              << "                 centered case, i.e. (0,0,1).\n";

    const Coord com(10.0, -5.0, 3.0);
    std::vector<Coord> templateIfaces { Coord(1.0, 0.0, 0.0), Coord(0.0, 1.0, 0.0) };
    // Absolute molecule interface coordinates = COM + internal coordinates.
    std::vector<Coord> molIfaces { com + templateIfaces[0], com + templateIfaces[1] };

    MolTemplate molTemplate { dn_make_template(templateIfaces, /*isRod=*/false) };
    Molecule oneMol { dn_make_molecule(com, molIfaces) };

    Vector normal(0.0, 0.0, 1.0);
    Vector result { determine_normal(normal, molTemplate, oneMol) };
    dn_report("returned normal", result);

    EXPECT_NEAR(result.x, 0.0, kDnTol) << "translation must not affect the normal (x)";
    EXPECT_NEAR(result.y, 0.0, kDnTol) << "translation must not affect the normal (y)";
    EXPECT_NEAR(result.z, 1.0, kDnTol) << "translation must not affect the normal (z)";
    EXPECT_NEAR(result.get_magnitude(), 1.0, kDnTol) << "result must be a unit vector";
}

// -----------------------------------------------------------------------------
// Test 7: `oneMol` is a by-value parameter, and the function mutates it
//         heavily (re-centering + rotations). The caller's Molecule must be
//         completely unaffected.
// -----------------------------------------------------------------------------
void test_dn_does_not_modify_caller_molecule()
{
    std::cerr << "\n[TEST] test_dn_does_not_modify_caller_molecule\n"
              << "  Source file:   src/reactions/determine_normal.cpp\n"
              << "  Function:      determine_normal()\n"
              << "  Scenario:      the Molecule is passed BY VALUE while the body\n"
              << "                 re-centers and rotates tmpComCoord/tmpICoords.\n"
              << "  Pass criteria: the caller's Molecule still holds its original\n"
              << "                 temporary association coordinates.\n";

    const Coord com(2.0, 3.0, 4.0);
    std::vector<Coord> templateIfaces { Coord(1.0, 0.0, 0.0), Coord(0.0, 1.0, 0.0) };
    std::vector<Coord> molIfaces { com + templateIfaces[0], com + templateIfaces[1] };

    MolTemplate molTemplate { dn_make_template(templateIfaces, /*isRod=*/false) };
    Molecule oneMol { dn_make_molecule(com, molIfaces) };

    // Snapshot before the call.
    const Coord comBefore { oneMol.tmpComCoord };
    const std::vector<Coord> ifacesBefore { oneMol.tmpICoords };

    Vector normal(0.0, 0.0, 1.0);
    Vector result { determine_normal(normal, molTemplate, oneMol) };
    dn_report("returned normal", result);

    // The temporary COM must be untouched.
    EXPECT_DOUBLE_EQ(oneMol.tmpComCoord.x, comBefore.x) << "caller's tmpComCoord.x changed";
    EXPECT_DOUBLE_EQ(oneMol.tmpComCoord.y, comBefore.y) << "caller's tmpComCoord.y changed";
    EXPECT_DOUBLE_EQ(oneMol.tmpComCoord.z, comBefore.z) << "caller's tmpComCoord.z changed";

    // Every temporary interface coordinate must be untouched.
    ASSERT_EQ(oneMol.tmpICoords.size(), ifacesBefore.size());
    for (std::size_t i { 0 }; i < ifacesBefore.size(); ++i) {
        EXPECT_DOUBLE_EQ(oneMol.tmpICoords[i].x, ifacesBefore[i].x)
            << "caller's tmpICoords[" << i << "].x changed";
        EXPECT_DOUBLE_EQ(oneMol.tmpICoords[i].y, ifacesBefore[i].y)
            << "caller's tmpICoords[" << i << "].y changed";
        EXPECT_DOUBLE_EQ(oneMol.tmpICoords[i].z, ifacesBefore[i].z)
            << "caller's tmpICoords[" << i << "].z changed";
    }
}

// -----------------------------------------------------------------------------
// Test 8: single-interface (rod) molecule. Only one interface-to-COM vector is
//         available, so the orientation is under-determined. The documented
//         contract is that the function either returns a unit normal or the
//         zero vector (the failure value), and that it must not crash.
// -----------------------------------------------------------------------------
void test_dn_single_interface_rod_molecule()
{
    std::cerr << "\n[TEST] test_dn_single_interface_rod_molecule\n"
              << "  Source file:   src/reactions/determine_normal.cpp\n"
              << "  Function:      determine_normal()\n"
              << "  Scenario:      rod molecule with a single interface that is\n"
              << "                 already in its template orientation.\n"
              << "  Pass criteria: the call completes and returns either a unit\n"
              << "                 vector (success) or exactly (0,0,0) -- the\n"
              << "                 documented failure value.\n";

    std::vector<Coord> ifaces { Coord(0.0, 0.0, 1.0) };
    MolTemplate molTemplate { dn_make_template(ifaces, /*isRod=*/true) };
    Molecule oneMol { dn_make_molecule(Coord(0.0, 0.0, 0.0), ifaces) };

    Vector normal(1.0, 0.0, 0.0);
    Vector result { determine_normal(normal, molTemplate, oneMol) };
    dn_report("returned normal", result);

    const double mag { result.get_magnitude() };
    const bool isZeroVector { mag < kDnTol };
    const bool isUnitVector { std::abs(mag - 1.0) < kDnTol };

    std::cerr << "    classified as: " << (isZeroVector ? "zero vector (failure path)"
                                                        : (isUnitVector ? "unit vector (success path)"
                                                                        : "NEITHER (unexpected)"))
              << '\n';

    EXPECT_TRUE(isZeroVector || isUnitVector)
        << "result must be either a unit normal or the (0,0,0) failure value, got |v| = " << mag;
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each scenario gets its own test case so that all of
// them execute even if one reports a failure (no fatal assertions are used
// except the size guard, which protects a subsequent loop).
// -----------------------------------------------------------------------------
TEST(DetermineNormal, PointMoleculeReturnsZero) { test_dn_point_molecule_returns_zero(); }
TEST(DetermineNormal, AlignedMoleculeReturnsInputNormal) { test_dn_aligned_molecule_returns_input_normal(); }
TEST(DetermineNormal, OutputIsNormalized) { test_dn_output_is_normalized(); }
TEST(DetermineNormal, RotationAboutZLeavesZNormalInvariant) { test_dn_rotation_about_z_leaves_z_normal_invariant(); }
TEST(DetermineNormal, RotationAboutXMovesNormalToYAxis) { test_dn_rotation_about_x_moves_normal_to_y_axis(); }
TEST(DetermineNormal, HandlesNonzeroComOffset) { test_dn_handles_nonzero_com_offset(); }
TEST(DetermineNormal, DoesNotModifyCallerMolecule) { test_dn_does_not_modify_caller_molecule(); }
TEST(DetermineNormal, SingleInterfaceRodMolecule) { test_dn_single_interface_rod_molecule(); }