/*! \file test_theta_rotation.cpp
 *
 * ### Unit test for src/reactions/theta_rotation.cpp
 *
 * The function under test is
 *
 * \code
 * void theta_rotation(Coord& reactIface1, Coord& reactIface2,
 *                     Molecule& reactMol1, Molecule& reactMol2,
 *                     double targAngle,
 *                     Complex& reactCom1, Complex& reactCom2,
 *                     std::vector<Molecule>& moleculeList);
 * \endcode
 *
 * It rotates the two associating Complexes (using their *temporary*
 * association coordinates, i.e. Molecule::tmpComCoord and
 * Molecule::tmpICoords) about the reacting interface `reactIface1` so that
 * the angle
 *
 *      theta = angle( sigma , v1 )
 *      sigma = reactIface1 - reactIface2
 *      v1    = reactIface1 - reactMol1.tmpComCoord
 *
 * becomes `targAngle`.
 *
 * The tests below check:
 *   1. The early-out for point particles (|v1| == 0).
 *   2. The early-out when theta is already at the requested value.
 *   3. That an actual rotation drives theta to the target value.
 *   4. That the rotation is a *rigid body* rotation about reactIface1, i.e.
 *      every distance to reactIface1 is conserved and reactIface1 itself
 *      never moves. This invariant holds no matter how the total rotation
 *      is split between the two complexes by determine_rotation_angles().
 *   5. That every member Molecule of a multi-molecule Complex is moved.
 *
 * Verbose output is written to stderr so a reader can follow exactly which
 * source file/function is exercised and what the pass criteria are.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers. Everything is prefixed with `thetarot_` so that the symbols
// cannot collide with other tests linked into the same binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully initialised Molecule with temporary association coords.
 *
 * theta_rotation() (and the rotate() helper it delegates to) works exclusively
 * on tmpComCoord / tmpICoords, so those must be populated. The "real"
 * coordinates are filled in as well so that the object is in a consistent
 * state.
 *
 * \param[in] index    index of this Molecule in moleculeList
 * \param[in] comIndex index of the parent Complex
 * \param[in] com      centre of mass coordinate
 * \param[in] ifaces   interface coordinates (absolute, not relative)
 */
Molecule thetarot_make_molecule(int index, int comIndex, const Coord& com,
                                const std::vector<Coord>& ifaces)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;

    mol.comCoord = com;
    mol.tmpComCoord = com;

    for (const auto& oneIfaceCrd : ifaces) {
        Molecule::Iface iface;
        iface.coord = oneIfaceCrd;
        iface.index = 0;
        iface.relIndex = static_cast<int>(mol.interfaceList.size());
        iface.molTypeIndex = 0;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);

        // the temporary (association) copy of the same coordinate
        mol.tmpICoords.push_back(oneIfaceCrd);
    }
    return mol;
}

/*! \brief Build a Complex owning the given member Molecule indices.
 *
 * Non-zero, isotropic translational and rotational diffusion constants are
 * used so that determine_rotation_angles() takes its generic (solution)
 * branch and splits the rotation between the two complexes.
 */
Complex thetarot_make_complex(int index, const std::vector<int>& members, const Coord& com)
{
    Complex targCom;
    targCom.index = index;
    targCom.id = index;
    targCom.comCoord = com;
    targCom.tmpComCoord = com;
    targCom.memberList = members;
    targCom.mass = static_cast<double>(members.size());
    targCom.radius = 1.0;
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Coord { 0.1, 0.1, 0.1 };
    targCom.isEmpty = false;
    targCom.OnSurface = false;
    return targCom;
}

/*! \brief Euclidean distance between two coordinates. */
double thetarot_dist(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Independently recompute the theta angle the function targets.
 *
 * theta = angle between sigma = (iface1 - iface2) and v1 = (iface1 - com1).
 * Computed here from first principles (not via Vector::dot_theta) so the
 * test does not simply mirror the implementation.
 */
double thetarot_theta(const Coord& iface1, const Coord& iface2, const Coord& com1)
{
    const double vx = iface1.x - com1.x;
    const double vy = iface1.y - com1.y;
    const double vz = iface1.z - com1.z;
    const double sx = iface1.x - iface2.x;
    const double sy = iface1.y - iface2.y;
    const double sz = iface1.z - iface2.z;

    const double dot = vx * sx + vy * sy + vz * sz;
    const double magV = std::sqrt(vx * vx + vy * vy + vz * vz);
    const double magS = std::sqrt(sx * sx + sy * sy + sz * sz);

    double cosTheta = dot / (magV * magS);
    if (cosTheta > 1.0)
        cosTheta = 1.0;
    if (cosTheta < -1.0)
        cosTheta = -1.0;
    return std::acos(cosTheta);
}

/*! \brief A copy of every temporary coordinate in the system. */
struct ThetaRotSnapshot {
    std::vector<Coord> comCoords;
    std::vector<std::vector<Coord>> ifaceCoords;
};

/*! \brief Take a snapshot of all tmpComCoord / tmpICoords values. */
ThetaRotSnapshot thetarot_snapshot(const std::vector<Molecule>& moleculeList)
{
    ThetaRotSnapshot snap;
    for (const auto& mol : moleculeList) {
        snap.comCoords.push_back(mol.tmpComCoord);
        snap.ifaceCoords.push_back(mol.tmpICoords);
    }
    return snap;
}

/*! \brief Assert that no temporary coordinate changed at all (bit-for-bit). */
void thetarot_expect_unchanged(const ThetaRotSnapshot& before,
                               const std::vector<Molecule>& moleculeList)
{
    for (std::size_t molItr = 0; molItr < moleculeList.size(); ++molItr) {
        const Molecule& mol = moleculeList[molItr];
        EXPECT_DOUBLE_EQ(mol.tmpComCoord.x, before.comCoords[molItr].x)
            << "tmpComCoord.x of molecule " << molItr << " should not move";
        EXPECT_DOUBLE_EQ(mol.tmpComCoord.y, before.comCoords[molItr].y)
            << "tmpComCoord.y of molecule " << molItr << " should not move";
        EXPECT_DOUBLE_EQ(mol.tmpComCoord.z, before.comCoords[molItr].z)
            << "tmpComCoord.z of molecule " << molItr << " should not move";

        for (std::size_t ifaceItr = 0; ifaceItr < mol.tmpICoords.size(); ++ifaceItr) {
            EXPECT_DOUBLE_EQ(mol.tmpICoords[ifaceItr].x, before.ifaceCoords[molItr][ifaceItr].x)
                << "tmpICoords[" << ifaceItr << "].x of molecule " << molItr << " should not move";
            EXPECT_DOUBLE_EQ(mol.tmpICoords[ifaceItr].y, before.ifaceCoords[molItr][ifaceItr].y)
                << "tmpICoords[" << ifaceItr << "].y of molecule " << molItr << " should not move";
            EXPECT_DOUBLE_EQ(mol.tmpICoords[ifaceItr].z, before.ifaceCoords[molItr][ifaceItr].z)
                << "tmpICoords[" << ifaceItr << "].z of molecule " << molItr << " should not move";
        }
    }
}

/*! \brief Verify the rigid-body invariant of a rotation about `origin`.
 *
 * Every centre of mass and every interface of every molecule must keep its
 * distance to the rotation origin, and molecules inside one complex must keep
 * their mutual distances.
 */
void thetarot_expect_rigid_about_origin(const ThetaRotSnapshot& before,
                                        const std::vector<Molecule>& moleculeList,
                                        const Coord& originBefore, const Coord& originAfter)
{
    for (std::size_t molItr = 0; molItr < moleculeList.size(); ++molItr) {
        const double dBefore = thetarot_dist(before.comCoords[molItr], originBefore);
        const double dAfter = thetarot_dist(moleculeList[molItr].tmpComCoord, originAfter);
        EXPECT_NEAR(dAfter, dBefore, 1e-9)
            << "COM of molecule " << molItr << " changed its distance to the rotation origin";

        for (std::size_t ifaceItr = 0; ifaceItr < moleculeList[molItr].tmpICoords.size(); ++ifaceItr) {
            const double diBefore = thetarot_dist(before.ifaceCoords[molItr][ifaceItr], originBefore);
            const double diAfter
                = thetarot_dist(moleculeList[molItr].tmpICoords[ifaceItr], originAfter);
            EXPECT_NEAR(diAfter, diBefore, 1e-9)
                << "Interface " << ifaceItr << " of molecule " << molItr
                << " changed its distance to the rotation origin";
        }
    }
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: point particle (interface sits on the centre of mass).
// -----------------------------------------------------------------------------
void test_thetarot_point_particle_no_rotation()
{
    std::cerr << "\n[TEST] test_thetarot_point_particle_no_rotation\n"
              << "  Source file:   src/reactions/theta_rotation.cpp\n"
              << "  Function:      theta_rotation()\n"
              << "  Scenario:      reactMol1's reacting interface coincides with its\n"
              << "                 centre of mass, so |v1| == 0 (a POINT particle).\n"
              << "  Pass criteria: the routine returns immediately and no temporary\n"
              << "                 coordinate in the system is modified.\n";

    // Molecule 0 is a point: interface == COM.
    std::vector<Molecule> moleculeList {
        thetarot_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 0.0, 0.0 } }),
        thetarot_make_molecule(1, 1, Coord { 1.0, -3.0, 0.0 }, { Coord { 1.0, -2.0, 0.0 } })
    };

    Complex com1 = thetarot_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 });
    Complex com2 = thetarot_make_complex(1, { 1 }, Coord { 1.0, -3.0, 0.0 });

    const ThetaRotSnapshot before = thetarot_snapshot(moleculeList);

    const double targAngle = M_PI / 2.0;
    std::cerr << "  Requested target theta = " << targAngle << " rad\n";
    std::cerr << "  Calling theta_rotation()...\n";

    theta_rotation(moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0],
        moleculeList[1], targAngle, com1, com2, moleculeList);

    // Nothing at all may have moved.
    thetarot_expect_unchanged(before, moleculeList);
    std::cerr << "  Verified: all tmpComCoord / tmpICoords are untouched.\n";
}

// -----------------------------------------------------------------------------
// Test 2: theta already equals the requested target -> nothing happens.
// -----------------------------------------------------------------------------
void test_thetarot_already_at_target_angle()
{
    std::cerr << "\n[TEST] test_thetarot_already_at_target_angle\n"
              << "  Source file:   src/reactions/theta_rotation.cpp\n"
              << "  Function:      theta_rotation()\n"
              << "  Scenario:      v1 = +x, sigma = +y so theta = pi/2 exactly, and\n"
              << "                 the target angle is also pi/2.\n"
              << "  Pass criteria: |target - current| < 1e-8 so the routine performs\n"
              << "                 no rotation; every coordinate is unchanged.\n";

    // mol0 COM (0,0,0), iface (1,0,0)  -> v1    = (1, 0,0)
    // mol1 iface (1,-2,0)              -> sigma = (0, 2,0)   => theta = pi/2
    std::vector<Molecule> moleculeList {
        thetarot_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 0.0, 0.0 } }),
        thetarot_make_molecule(1, 1, Coord { 1.0, -3.0, 0.0 }, { Coord { 1.0, -2.0, 0.0 } })
    };

    Complex com1 = thetarot_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 });
    Complex com2 = thetarot_make_complex(1, { 1 }, Coord { 1.0, -3.0, 0.0 });

    const double thetaBefore = thetarot_theta(
        moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0].tmpComCoord);
    std::cerr << "  Initial theta = " << thetaBefore << " rad (expected pi/2 = "
              << M_PI / 2.0 << ")\n";
    EXPECT_NEAR(thetaBefore, M_PI / 2.0, 1e-12) << "test setup should start at pi/2";

    const ThetaRotSnapshot before = thetarot_snapshot(moleculeList);

    std::cerr << "  Calling theta_rotation() with target = pi/2...\n";
    theta_rotation(moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0],
        moleculeList[1], M_PI / 2.0, com1, com2, moleculeList);

    thetarot_expect_unchanged(before, moleculeList);
    std::cerr << "  Verified: no coordinate moved because theta already matched.\n";
}

// -----------------------------------------------------------------------------
// Test 3: already anti-parallel (theta == pi) with target pi -> nothing happens.
// -----------------------------------------------------------------------------
void test_thetarot_already_antiparallel()
{
    std::cerr << "\n[TEST] test_thetarot_already_antiparallel\n"
              << "  Source file:   src/reactions/theta_rotation.cpp\n"
              << "  Function:      theta_rotation()\n"
              << "  Scenario:      sigma and v1 are exactly anti-parallel (theta = pi)\n"
              << "                 and the requested target angle is pi as well.\n"
              << "  Pass criteria: no rotation is performed (early return path).\n";

    // mol0 COM (0,0,0), iface (1,0,0) -> v1    = ( 1,0,0)
    // mol1 iface (2,0,0)              -> sigma = (-1,0,0)  => theta = pi
    std::vector<Molecule> moleculeList {
        thetarot_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 0.0, 0.0 } }),
        thetarot_make_molecule(1, 1, Coord { 3.0, 0.0, 0.0 }, { Coord { 2.0, 0.0, 0.0 } })
    };

    Complex com1 = thetarot_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 });
    Complex com2 = thetarot_make_complex(1, { 1 }, Coord { 3.0, 0.0, 0.0 });

    const double thetaBefore = thetarot_theta(
        moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0].tmpComCoord);
    std::cerr << "  Initial theta = " << thetaBefore << " rad (expected pi = " << M_PI << ")\n";
    EXPECT_NEAR(thetaBefore, M_PI, 1e-12) << "test setup should start anti-parallel";

    const ThetaRotSnapshot before = thetarot_snapshot(moleculeList);

    std::cerr << "  Calling theta_rotation() with target = pi...\n";
    theta_rotation(moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0],
        moleculeList[1], M_PI, com1, com2, moleculeList);

    thetarot_expect_unchanged(before, moleculeList);
    std::cerr << "  Verified: anti-parallel configuration was left alone.\n";
}

// -----------------------------------------------------------------------------
// Test 4: an actual rotation from pi/2 to pi/3.
// -----------------------------------------------------------------------------
void test_thetarot_rotates_to_target()
{
    std::cerr << "\n[TEST] test_thetarot_rotates_to_target\n"
              << "  Source file:   src/reactions/theta_rotation.cpp\n"
              << "  Function:      theta_rotation()\n"
              << "  Scenario:      theta starts at pi/2, target is pi/3. Both complexes\n"
              << "                 have identical, non-zero diffusion constants so the\n"
              << "                 rotation is shared between them.\n"
              << "  Pass criteria: (a) the recomputed theta equals pi/3 to 1e-6,\n"
              << "                 (b) reactIface1 (the rotation origin) never moves,\n"
              << "                 (c) every distance to reactIface1 is conserved.\n";

    std::vector<Molecule> moleculeList {
        thetarot_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 0.0, 0.0 } }),
        thetarot_make_molecule(1, 1, Coord { 1.0, -3.0, 0.0 }, { Coord { 1.0, -2.0, 0.0 } })
    };

    Complex com1 = thetarot_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 });
    Complex com2 = thetarot_make_complex(1, { 1 }, Coord { 1.0, -3.0, 0.0 });

    const Coord originBefore = moleculeList[0].tmpICoords[0]; // reactIface1
    const ThetaRotSnapshot before = thetarot_snapshot(moleculeList);

    const double thetaBefore = thetarot_theta(
        moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0].tmpComCoord);
    const double targAngle = M_PI / 3.0;
    std::cerr << "  Initial theta = " << thetaBefore << " rad, target = " << targAngle << " rad\n";

    std::cerr << "  Calling theta_rotation()...\n";
    theta_rotation(moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0],
        moleculeList[1], targAngle, com1, com2, moleculeList);

    const Coord originAfter = moleculeList[0].tmpICoords[0];
    const double thetaAfter = thetarot_theta(
        moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0].tmpComCoord);
    std::cerr << "  Resulting theta = " << thetaAfter << " rad\n";

    // (a) the geometric goal of the routine
    EXPECT_NEAR(thetaAfter, targAngle, 1e-6)
        << "theta_rotation should drive theta onto the requested target angle";

    // (b) the reacting interface of molecule 1 is the rotation origin: it must not move
    EXPECT_NEAR(originAfter.x, originBefore.x, 1e-12) << "reactIface1.x must stay put";
    EXPECT_NEAR(originAfter.y, originBefore.y, 1e-12) << "reactIface1.y must stay put";
    EXPECT_NEAR(originAfter.z, originBefore.z, 1e-12) << "reactIface1.z must stay put";

    // (c) rigid body: distances to the rotation origin conserved everywhere
    thetarot_expect_rigid_about_origin(before, moleculeList, originBefore, originAfter);

    // Something really did move (otherwise (a)-(c) would be trivially satisfied)
    const double comMoved = thetarot_dist(moleculeList[0].tmpComCoord, before.comCoords[0]);
    std::cerr << "  Molecule 0 COM displacement = " << comMoved << " nm\n";
    EXPECT_GT(comMoved, 1e-6) << "the first complex should actually have been rotated";
}

// -----------------------------------------------------------------------------
// Test 5: rotate from pi/2 all the way to pi (target is a parallel/limit angle).
// -----------------------------------------------------------------------------
void test_thetarot_rotates_to_pi()
{
    std::cerr << "\n[TEST] test_thetarot_rotates_to_pi\n"
              << "  Source file:   src/reactions/theta_rotation.cpp\n"
              << "  Function:      theta_rotation()\n"
              << "  Scenario:      theta starts at pi/2 and the target is exactly pi,\n"
              << "                 exercising the areSameAngle(targAngle, M_PI) branch\n"
              << "                 taken after the rotation.\n"
              << "  Pass criteria: the recomputed theta equals pi to 1e-5 (acos is\n"
              << "                 numerically flat near -1, hence the looser bound).\n";

    std::vector<Molecule> moleculeList {
        thetarot_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 0.0, 0.0 } }),
        thetarot_make_molecule(1, 1, Coord { 1.0, -3.0, 0.0 }, { Coord { 1.0, -2.0, 0.0 } })
    };

    Complex com1 = thetarot_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 });
    Complex com2 = thetarot_make_complex(1, { 1 }, Coord { 1.0, -3.0, 0.0 });

    const Coord originBefore = moleculeList[0].tmpICoords[0];
    const ThetaRotSnapshot before = thetarot_snapshot(moleculeList);

    std::cerr << "  Initial theta = "
              << thetarot_theta(moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0],
                     moleculeList[0].tmpComCoord)
              << " rad, target = " << M_PI << " rad\n";

    std::cerr << "  Calling theta_rotation()...\n";
    theta_rotation(moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0],
        moleculeList[1], M_PI, com1, com2, moleculeList);

    const Coord originAfter = moleculeList[0].tmpICoords[0];
    const double thetaAfter = thetarot_theta(
        moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0].tmpComCoord);
    std::cerr << "  Resulting theta = " << thetaAfter << " rad\n";

    EXPECT_NEAR(thetaAfter, M_PI, 1e-5)
        << "theta_rotation should be able to reach the limiting angle pi";

    // The rotation origin must still be fixed and the motion must be rigid.
    EXPECT_NEAR(originAfter.x, originBefore.x, 1e-12) << "reactIface1.x must stay put";
    EXPECT_NEAR(originAfter.y, originBefore.y, 1e-12) << "reactIface1.y must stay put";
    EXPECT_NEAR(originAfter.z, originBefore.z, 1e-12) << "reactIface1.z must stay put";
    thetarot_expect_rigid_about_origin(before, moleculeList, originBefore, originAfter);
}

// -----------------------------------------------------------------------------
// Test 6: the second complex holds two molecules -> both must be rotated.
// -----------------------------------------------------------------------------
void test_thetarot_multi_molecule_complex_moves_all_members()
{
    std::cerr << "\n[TEST] test_thetarot_multi_molecule_complex_moves_all_members\n"
              << "  Source file:   src/reactions/theta_rotation.cpp\n"
              << "  Function:      theta_rotation() (via rotate())\n"
              << "  Scenario:      complex 1 owns molecules 1 and 2; only molecule 1\n"
              << "                 carries the reacting interface.\n"
              << "  Pass criteria: molecule 2 is moved as well, the internal distance\n"
              << "                 between molecules 1 and 2 is conserved, and the\n"
              << "                 final theta reaches the target.\n";

    std::vector<Molecule> moleculeList {
        thetarot_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 0.0, 0.0 } }),
        thetarot_make_molecule(1, 1, Coord { 1.0, -3.0, 0.0 }, { Coord { 1.0, -2.0, 0.0 } }),
        thetarot_make_molecule(2, 1, Coord { 1.0, -5.0, 0.0 }, { Coord { 1.0, -4.0, 0.0 } })
    };

    Complex com1 = thetarot_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 });
    Complex com2 = thetarot_make_complex(1, { 1, 2 }, Coord { 1.0, -4.0, 0.0 });

    const Coord originBefore = moleculeList[0].tmpICoords[0];
    const ThetaRotSnapshot before = thetarot_snapshot(moleculeList);
    const double internalDistBefore
        = thetarot_dist(moleculeList[1].tmpComCoord, moleculeList[2].tmpComCoord);

    const double targAngle = M_PI / 3.0;
    std::cerr << "  Initial theta = "
              << thetarot_theta(moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0],
                     moleculeList[0].tmpComCoord)
              << " rad, target = " << targAngle << " rad\n";

    std::cerr << "  Calling theta_rotation()...\n";
    theta_rotation(moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0],
        moleculeList[1], targAngle, com1, com2, moleculeList);

    const Coord originAfter = moleculeList[0].tmpICoords[0];
    const double thetaAfter = thetarot_theta(
        moleculeList[0].tmpICoords[0], moleculeList[1].tmpICoords[0], moleculeList[0].tmpComCoord);
    std::cerr << "  Resulting theta = " << thetaAfter << " rad\n";
    EXPECT_NEAR(thetaAfter, targAngle, 1e-6)
        << "theta must reach the target even with a multi-molecule complex";

    // Molecule 2 belongs to the rotated complex, so it must have been displaced.
    const double mol2Moved = thetarot_dist(moleculeList[2].tmpComCoord, before.comCoords[2]);
    std::cerr << "  Molecule 2 (non-reacting member) displacement = " << mol2Moved << " nm\n";
    EXPECT_GT(mol2Moved, 1e-6) << "every member of the complex must be rotated";

    // The complex must move as a rigid body: internal distance conserved.
    const double internalDistAfter
        = thetarot_dist(moleculeList[1].tmpComCoord, moleculeList[2].tmpComCoord);
    std::cerr << "  Internal COM-COM distance before/after = " << internalDistBefore << " / "
              << internalDistAfter << " nm\n";
    EXPECT_NEAR(internalDistAfter, internalDistBefore, 1e-9)
        << "the internal geometry of the complex must be preserved";

    // And all distances to the rotation origin are conserved for every molecule.
    thetarot_expect_rigid_about_origin(before, moleculeList, originBefore, originAfter);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is invoked from its own TEST so that
// a failure in one scenario does not prevent the remaining ones from running.
// -----------------------------------------------------------------------------
TEST(ThetaRotation, PointParticleNoRotation) { test_thetarot_point_particle_no_rotation(); }
TEST(ThetaRotation, AlreadyAtTargetAngle) { test_thetarot_already_at_target_angle(); }
TEST(ThetaRotation, AlreadyAntiparallel) { test_thetarot_already_antiparallel(); }
TEST(ThetaRotation, RotatesToTarget) { test_thetarot_rotates_to_target(); }
TEST(ThetaRotation, RotatesToPi) { test_thetarot_rotates_to_pi(); }
TEST(ThetaRotation, MultiMoleculeComplexMovesAllMembers)
{
    test_thetarot_multi_molecule_complex_moves_all_members();
}