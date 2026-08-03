/*! \file test_check_bases.cpp
 *
 * ### Unit test for src/reactions/check_bases.cpp
 *
 * The single function under test is
 *
 *     void check_bases(bool& cancelAssoc,
 *                      const Coord& reactIface1, const Coord& reactIface2,
 *                      int ifaceIndex1, int ifaceIndex2,
 *                      const Molecule& reactMol1, const Molecule& reactMol2,
 *                      const Complex& reactCom1, const Complex& reactCom2,
 *                      const ForwardRxn& currRxn,
 *                      const std::vector<Molecule>& moleculeList,
 *                      const std::vector<MolTemplate>& molTemplateList);
 *
 * check_bases() is the post-association sanity check.  It never "un-cancels" an
 * association; it only ever sets `cancelAssoc = true` and returns early.  The
 * checks it performs, in order, are:
 *
 *   1. theta1  -- compares the target theta1 of the reaction against the angle
 *                 between sigma and (iface1 - tmpCOM1).  NOTE: as written in the
 *                 source the comparison is inverted with respect to the other
 *                 angle checks: a *match* within 1e-6 sets cancelAssoc.  The
 *                 tests below document that implemented behaviour.
 *   2. theta2  -- cancels when the target theta2 matches neither the measured
 *                 angle nor (2*pi - measured angle), compared at roundv()
 *                 precision (4 decimal places).
 *   3. phi1 / phi2 / omega -- skipped when the target angle is NaN, otherwise
 *                 recomputed (calculate_phi()/calculate_omega()) and compared.
 *   4. conservedMags()  -- COM->interface vector lengths must be unchanged
 *                 between the real coordinates and the temporary (association)
 *                 coordinates, for every member of both complexes.
 *   5. conservedRigid() -- internal rigidity of both complexes must be intact.
 *
 * The tests build a small, fully specified two-molecule scene so that every
 * measured angle is known analytically, then flip one input at a time and check
 * how `cancelAssoc` responds.  All assertions are non-fatal (EXPECT_*) so the
 * whole file always runs to completion.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers.  Everything lives in an anonymous namespace and/or carries the
// "cb_" (check_bases) prefix so nothing collides with the rest of the suite.
// -----------------------------------------------------------------------------
namespace {

//! Convenience: NaN, which is what ForwardRxn::assocAngles defaults to and what
//! makes check_bases() skip the phi/omega branches entirely.
const double kCbNaN = std::numeric_limits<double>::quiet_NaN();

/*! \brief Angle (radians) that check_bases() will measure for theta1 in the
 * scene built by cb_build_scene(): sigma = (-1,0,0), v1 = (1,-1,0) => 3*pi/4. */
const double kCbTheta1Measured = 3.0 * M_PI / 4.0;

/*! \brief Angle (radians) that check_bases() will measure for theta2 in the
 * scene built by cb_build_scene(): -sigma = (1,0,0), v2 = (0,-1,0) => pi/2. */
const double kCbTheta2Measured = M_PI / 2.0;

/*! \brief A tiny container holding everything check_bases() needs. */
struct CbScene {
    std::vector<Molecule> moleculeList; //!< molecule 0 and molecule 1
    std::vector<MolTemplate> molTemplateList; //!< one template per molecule type
    Complex com1; //!< parent complex of molecule 0
    Complex com2; //!< parent complex of molecule 1
    Coord iface1; //!< coordinate passed as reactIface1
    Coord iface2; //!< coordinate passed as reactIface2
};

/*! \brief Build a molecule whose temporary (association) coordinates initially
 * mirror its real coordinates -- i.e. a configuration that trivially satisfies
 * conservedMags() and conservedRigid().
 *
 * \param[in] index        index of this molecule inside moleculeList
 * \param[in] molTypeIndex index of the matching MolTemplate
 * \param[in] com          center of mass coordinate
 * \param[in] ifaceCoords  absolute interface coordinates
 */
Molecule cb_make_molecule(int index, int molTypeIndex, const Coord& com,
    const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = index; // one molecule per complex in these tests
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.tmpComCoord = com; // association coords start identical

    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface {};
        iface.coord = ifaceCoords[i];
        iface.index = static_cast<int>(i);
        iface.relIndex = static_cast<int>(i);
        iface.molTypeIndex = molTypeIndex;
        mol.interfaceList.push_back(iface);
        mol.tmpICoords.push_back(ifaceCoords[i]);
    }
    return mol;
}

/*! \brief Build a single-member Complex. */
Complex cb_make_complex(int index, int memberMolIndex, const Coord& com)
{
    Complex complex;
    complex.index = index;
    complex.comCoord = com;
    complex.tmpComCoord = com;
    complex.mass = 1.0;
    complex.radius = 2.0;
    complex.memberList.clear();
    complex.memberList.push_back(memberMolIndex);
    complex.D = Coord { 10.0, 10.0, 10.0 };
    complex.Dr = Coord { 0.01, 0.01, 0.01 };
    complex.OnSurface = false;
    return complex;
}

/*! \brief Build a MolTemplate whose internal coordinates match the molecule's
 * COM->interface vectors.  Needed because calculate_phi()/calculate_omega() can
 * re-orient a molecule onto its template.
 *
 * \param[in] molTypeIndex          template index
 * \param[in] name                  template (molecule) name
 * \param[in] internalIfaceCoords   interface coords relative to a COM at origin
 */
MolTemplate cb_make_template(int molTypeIndex, const std::string& name,
    const std::vector<Coord>& internalIfaceCoords)
{
    MolTemplate molTemplate;
    molTemplate.molTypeIndex = molTypeIndex;
    molTemplate.molName = name;
    molTemplate.comCoord = Coord { 0.0, 0.0, 0.0 };
    molTemplate.mass = 1.0;
    molTemplate.D = Coord { 10.0, 10.0, 10.0 };
    molTemplate.Dr = Coord { 0.01, 0.01, 0.01 };
    molTemplate.isRod = false;
    molTemplate.isPoint = false;
    molTemplate.isLipid = false;
    molTemplate.isImplicitLipid = false;
    molTemplate.checkOverlap = false;
    molTemplate.copies = 1;

    double maxRadius = 0.0;
    for (std::size_t i = 0; i < internalIfaceCoords.size(); ++i) {
        Interface iface {};
        iface.index = static_cast<int>(i);
        iface.name = "i" + std::to_string(i);
        iface.iCoord = internalIfaceCoords[i];

        // Give the interface a single (default) state so the template is
        // well-formed for any code that inspects stateList.
        Interface::State state {};
        state.index = static_cast<int>(i);
        state.iden = '\0';
        state.ifaceAndStateName = iface.name;
        iface.stateList.push_back(state);

        molTemplate.interfaceList.push_back(iface);

        Coord c = internalIfaceCoords[i];
        const double mag = std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
        if (mag > maxRadius)
            maxRadius = mag;
    }
    molTemplate.radius = (maxRadius > 0.0) ? maxRadius : 0.0001;
    return molTemplate;
}

/*! \brief Build the reference scene used by every test.
 *
 * Geometry (all coordinates absolute, association/tmp coords == real coords):
 *
 *   molecule 0 : COM (0,1,0)  interface (1,0,0)   -> v1 = ( 1,-1,0), |v1| = sqrt(2)
 *   molecule 1 : COM (2,1,0)  interface (2,0,0)   -> v2 = ( 0,-1,0), |v2| = 1
 *   sigma      = iface1 - iface2 = (-1,0,0), |sigma| = 1 (== bindRadius)
 *
 * Therefore the angles check_bases() will *measure* are
 *   theta1 = acos(-1/sqrt(2)) = 3*pi/4
 *   theta2 = acos(0)          = pi/2
 */
CbScene cb_build_scene()
{
    CbScene scene;

    const Coord com1 { 0.0, 1.0, 0.0 };
    const Coord iface1 { 1.0, 0.0, 0.0 };
    const Coord com2 { 2.0, 1.0, 0.0 };
    const Coord iface2 { 2.0, 0.0, 0.0 };

    scene.moleculeList.push_back(cb_make_molecule(0, 0, com1, { iface1 }));
    scene.moleculeList.push_back(cb_make_molecule(1, 1, com2, { iface2 }));

    scene.molTemplateList.push_back(
        cb_make_template(0, "A", { Coord { iface1.x - com1.x, iface1.y - com1.y, iface1.z - com1.z } }));
    scene.molTemplateList.push_back(
        cb_make_template(1, "B", { Coord { iface2.x - com2.x, iface2.y - com2.y, iface2.z - com2.z } }));

    scene.com1 = cb_make_complex(0, 0, com1);
    scene.com2 = cb_make_complex(1, 1, com2);

    scene.iface1 = iface1;
    scene.iface2 = iface2;

    return scene;
}

/*! \brief Build a ForwardRxn with the requested target association angles.
 *
 * Pass kCbNaN for any angle that should be skipped by check_bases().
 */
ForwardRxn cb_make_rxn(double theta1, double theta2, double phi1, double phi2, double omega)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.bindRadius = 1.0; // matches |sigma| of the reference scene
    rxn.bindRadius2D = 1.0;
    rxn.isReversible = false;
    rxn.norm1 = Vector(0.0, 0.0, 1.0);
    rxn.norm2 = Vector(0.0, 0.0, 1.0);
    rxn.assocAngles = ForwardRxn::Angles(theta1, theta2, phi1, phi2, omega);

    // Populate the reactant list so the reaction object is well formed.
    rxn.reactantListNew.push_back(RxnIface("i0", 0, 0, 0, '\0', false));
    rxn.reactantListNew.push_back(RxnIface("i0", 1, 1, 0, '\0', false));

    return rxn;
}

/*! \brief Rotate a coordinate about the z-axis (used for the rigid-body test). */
Coord cb_rotate_z(const Coord& c, double ang)
{
    return Coord { c.x * std::cos(ang) - c.y * std::sin(ang),
        c.x * std::sin(ang) + c.y * std::cos(ang),
        c.z };
}

/*! \brief Thin wrapper that runs check_bases() on a scene and reports the flag. */
bool cb_run_check(CbScene& scene, const ForwardRxn& rxn, bool initialFlag = false)
{
    bool cancelAssoc { initialFlag };
    check_bases(cancelAssoc, scene.iface1, scene.iface2, /*ifaceIndex1=*/0, /*ifaceIndex2=*/0,
        scene.moleculeList[0], scene.moleculeList[1], scene.com1, scene.com2, rxn,
        scene.moleculeList, scene.molTemplateList);
    std::cerr << "    check_bases() returned cancelAssoc = " << std::boolalpha << cancelAssoc << '\n';
    return cancelAssoc;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a fully self-consistent configuration is accepted.
//
// theta1 is deliberately given a value that does NOT match the measured angle
// (because the source cancels on a match), theta2 is given the measured value,
// and phi1/phi2/omega are NaN so they are skipped.  The temporary coordinates
// equal the real coordinates, so conservedMags()/conservedRigid() must pass.
// -----------------------------------------------------------------------------
void test_check_bases_accepts_consistent_geometry()
{
    std::cerr << "\n[TEST] test_check_bases_accepts_consistent_geometry\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases()\n"
              << "  Scenario:      theta1 target (0.5) differs from measured "
              << kCbTheta1Measured << ", theta2 target == measured "
              << kCbTheta2Measured << ", phi1/phi2/omega are NaN,\n"
              << "                 tmp coords == real coords.\n"
              << "  Pass criteria: cancelAssoc stays false (association accepted).\n";

    CbScene scene = cb_build_scene();
    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, /*theta2=*/kCbTheta2Measured, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_FALSE(cancelAssoc)
        << "A geometrically consistent association with NaN dihedrals should not be cancelled";
}

// -----------------------------------------------------------------------------
// Test 2: the theta1 branch.
//
// The implementation cancels when |theta1_target - theta1_measured| < 1e-6, so
// handing it the measured value must set the flag.  This test pins that
// (inverted-looking) behaviour so a future change is noticed.
// -----------------------------------------------------------------------------
void test_check_bases_theta1_branch()
{
    std::cerr << "\n[TEST] test_check_bases_theta1_branch\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- theta1 comparison\n"
              << "  Scenario:      theta1 target set to the measured value "
              << kCbTheta1Measured << " (3*pi/4).\n"
              << "  Pass criteria: cancelAssoc becomes true, because the source\n"
              << "                 cancels when |target - measured| < 1e-6.\n";

    CbScene scene = cb_build_scene();
    ForwardRxn rxn = cb_make_rxn(kCbTheta1Measured, kCbTheta2Measured, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_TRUE(cancelAssoc)
        << "check_bases() sets cancelAssoc when theta1 is within 1e-6 of the measured angle";
}

// -----------------------------------------------------------------------------
// Test 3: the theta2 branch rejects a mismatching target angle.
// -----------------------------------------------------------------------------
void test_check_bases_theta2_mismatch_cancels()
{
    std::cerr << "\n[TEST] test_check_bases_theta2_mismatch_cancels\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- theta2 comparison\n"
              << "  Scenario:      theta2 target = 1.0 rad while the measured angle is "
              << kCbTheta2Measured << ".\n"
              << "  Pass criteria: cancelAssoc becomes true (theta2 matches neither\n"
              << "                 the measured angle nor 2*pi - measured).\n";

    CbScene scene = cb_build_scene();
    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, /*theta2=*/1.0, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_TRUE(cancelAssoc) << "A theta2 target of 1.0 rad must be rejected for a pi/2 geometry";
}

// -----------------------------------------------------------------------------
// Test 4: the theta2 branch also accepts the (2*pi - measured) alternative.
//
// 2*pi - pi/2 == 3*pi/2, so a target of 3*pi/2 satisfies the second clause of
// the theta2 test and must NOT cancel.
// -----------------------------------------------------------------------------
void test_check_bases_theta2_supplementary_angle_accepted()
{
    std::cerr << "\n[TEST] test_check_bases_theta2_supplementary_angle_accepted\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- theta2 comparison (2*pi - measured)\n"
              << "  Scenario:      theta2 target = 3*pi/2 = " << (3.0 * M_PI / 2.0)
              << ", measured = " << kCbTheta2Measured << ", and 2*pi - measured = "
              << (2.0 * M_PI - kCbTheta2Measured) << ".\n"
              << "  Pass criteria: cancelAssoc stays false (second clause satisfied).\n";

    CbScene scene = cb_build_scene();
    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, /*theta2=*/3.0 * M_PI / 2.0, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_FALSE(cancelAssoc)
        << "theta2 == 2*pi - measured must be accepted by check_bases()";
}

// -----------------------------------------------------------------------------
// Test 5: NaN dihedrals are skipped, finite bogus dihedrals are rejected.
//
// phi1 is given an impossible value (99 rad), which can never equal whatever
// calculate_phi() reports, so the phi1 branch must fire.
// -----------------------------------------------------------------------------
void test_check_bases_phi1_mismatch_cancels()
{
    std::cerr << "\n[TEST] test_check_bases_phi1_mismatch_cancels\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- phi1 branch (calculate_phi)\n"
              << "  Scenario:      phi1 target = 99.0 rad, an angle no geometry can produce.\n"
              << "  Pass criteria: cancelAssoc becomes true (phi1 branch is entered\n"
              << "                 because phi1 is not NaN, and the comparison fails).\n";

    CbScene scene = cb_build_scene();
    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, kCbTheta2Measured, /*phi1=*/99.0, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_TRUE(cancelAssoc) << "A non-NaN, unattainable phi1 target must cancel the association";
}

// -----------------------------------------------------------------------------
// Test 6: same idea for phi2 (phi1 left NaN so the phi2 branch is reached).
// -----------------------------------------------------------------------------
void test_check_bases_phi2_mismatch_cancels()
{
    std::cerr << "\n[TEST] test_check_bases_phi2_mismatch_cancels\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- phi2 branch (calculate_phi)\n"
              << "  Scenario:      phi1 = NaN (skipped), phi2 target = 99.0 rad.\n"
              << "  Pass criteria: cancelAssoc becomes true.\n";

    CbScene scene = cb_build_scene();
    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, kCbTheta2Measured, kCbNaN, /*phi2=*/99.0, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_TRUE(cancelAssoc) << "A non-NaN, unattainable phi2 target must cancel the association";
}

// -----------------------------------------------------------------------------
// Test 7: the omega branch.  Both thetas are far from pi in this reaction, so
// calculate_omega() uses the COM->interface vectors (not the norms).
// -----------------------------------------------------------------------------
void test_check_bases_omega_mismatch_cancels()
{
    std::cerr << "\n[TEST] test_check_bases_omega_mismatch_cancels\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- omega branch (calculate_omega)\n"
              << "  Scenario:      phi1/phi2 = NaN (skipped), omega target = 99.0 rad.\n"
              << "  Pass criteria: cancelAssoc becomes true.\n";

    CbScene scene = cb_build_scene();
    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, kCbTheta2Measured, kCbNaN, kCbNaN, /*omega=*/99.0);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_TRUE(cancelAssoc) << "A non-NaN, unattainable omega target must cancel the association";
}

// -----------------------------------------------------------------------------
// Test 8: conservedMags() failure for the first complex.
//
// The temporary interface coordinate of molecule 0 is moved so that its
// COM->interface length changes from sqrt(2) to 5, while the coordinates handed
// to check_bases() for the angle tests are untouched.  Only the magnitude
// conservation check can therefore be the reason for cancelling.
// -----------------------------------------------------------------------------
void test_check_bases_broken_magnitude_complex1_cancels()
{
    std::cerr << "\n[TEST] test_check_bases_broken_magnitude_complex1_cancels\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- conservedMags(reactCom1, ...)\n"
              << "  Scenario:      molecule 0 tmp interface moved from (1,0,0) to (5,1,0),\n"
              << "                 changing its COM->interface length from sqrt(2) to 5.\n"
              << "  Pass criteria: cancelAssoc becomes true (magnitude not conserved).\n";

    CbScene scene = cb_build_scene();
    // Break only the temporary (association) interface coordinate.
    scene.moleculeList[0].tmpICoords[0] = Coord { 5.0, 1.0, 0.0 };

    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, kCbTheta2Measured, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_TRUE(cancelAssoc)
        << "check_bases() must cancel when a COM->interface length of complex 1 changed";
}

// -----------------------------------------------------------------------------
// Test 9: conservedMags() failure for the second complex.
// -----------------------------------------------------------------------------
void test_check_bases_broken_magnitude_complex2_cancels()
{
    std::cerr << "\n[TEST] test_check_bases_broken_magnitude_complex2_cancels\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- conservedMags(reactCom2, ...)\n"
              << "  Scenario:      molecule 1 tmp interface moved from (2,0,0) to (2,5,0),\n"
              << "                 changing its COM->interface length from 1 to 4.\n"
              << "  Pass criteria: cancelAssoc becomes true (magnitude not conserved).\n";

    CbScene scene = cb_build_scene();
    scene.moleculeList[1].tmpICoords[0] = Coord { 2.0, 5.0, 0.0 };

    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, kCbTheta2Measured, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_TRUE(cancelAssoc)
        << "check_bases() must cancel when a COM->interface length of complex 2 changed";
}

// -----------------------------------------------------------------------------
// Test 10: a pure rigid-body motion of the association coordinates is accepted.
//
// Every temporary coordinate (COMs and interfaces) is the real coordinate
// rotated 30 degrees about z.  All lengths and internal angles are preserved, so
// conservedMags()/conservedRigid() must both succeed and the measured thetas are
// unchanged.
// -----------------------------------------------------------------------------
void test_check_bases_rigid_rotation_accepted()
{
    const double ang = M_PI / 6.0; // 30 degrees

    std::cerr << "\n[TEST] test_check_bases_rigid_rotation_accepted\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- conservedMags()/conservedRigid()\n"
              << "  Scenario:      all tmp coordinates are the real coordinates rotated "
              << (ang * 180.0 / M_PI) << " deg about z.\n"
              << "  Pass criteria: cancelAssoc stays false; a rigid-body move conserves\n"
              << "                 both vector magnitudes and internal rigidity.\n";

    CbScene scene = cb_build_scene();

    // Rotate the temporary coordinates of both molecules.
    for (auto& mol : scene.moleculeList) {
        mol.tmpComCoord = cb_rotate_z(mol.comCoord, ang);
        for (std::size_t i = 0; i < mol.interfaceList.size(); ++i)
            mol.tmpICoords[i] = cb_rotate_z(mol.interfaceList[i].coord, ang);
    }
    // The reacting interface coordinates handed to check_bases() must be the
    // rotated ones so the geometry stays self-consistent.
    scene.iface1 = scene.moleculeList[0].tmpICoords[0];
    scene.iface2 = scene.moleculeList[1].tmpICoords[0];
    scene.com1.tmpComCoord = scene.moleculeList[0].tmpComCoord;
    scene.com2.tmpComCoord = scene.moleculeList[1].tmpComCoord;

    // Rotation preserves relative geometry, so the measured thetas are the same
    // as in the unrotated scene.
    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, kCbTheta2Measured, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_FALSE(cancelAssoc)
        << "A rigid-body rotation of the association coordinates must be accepted";
}

// -----------------------------------------------------------------------------
// Test 11: check_bases() never clears an already-set cancel flag.
//
// The scene/reaction pair used here is the "accepted" configuration from test 1,
// but the caller arrives with cancelAssoc already true (e.g. an earlier overlap
// check failed).  The flag must survive.
// -----------------------------------------------------------------------------
void test_check_bases_never_clears_cancel_flag()
{
    std::cerr << "\n[TEST] test_check_bases_never_clears_cancel_flag\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases()\n"
              << "  Scenario:      an otherwise acceptable association is checked while\n"
              << "                 cancelAssoc is already true on entry.\n"
              << "  Pass criteria: cancelAssoc is still true on return -- the function\n"
              << "                 only ever sets the flag, it never resets it.\n";

    CbScene scene = cb_build_scene();
    ForwardRxn rxn = cb_make_rxn(/*theta1=*/0.5, kCbTheta2Measured, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn, /*initialFlag=*/true);
    EXPECT_TRUE(cancelAssoc) << "check_bases() must not clear a cancel flag set by a previous check";
}

// -----------------------------------------------------------------------------
// Test 12: the theta comparisons are performed against the *temporary*
// (association) center of mass, not the real one.
//
// Here the real COMs are left alone but molecule 0's tmpComCoord is moved to the
// interface's mirror point so that v1 = iface1 - tmpCOM1 becomes (1,0,0), i.e.
// exactly anti-parallel to sigma => measured theta1 = pi.  Feeding theta1 = pi
// then triggers the (match-cancels) theta1 branch, proving the tmp COM is the
// quantity actually used.
// -----------------------------------------------------------------------------
void test_check_bases_uses_tmp_com_coords()
{
    std::cerr << "\n[TEST] test_check_bases_uses_tmp_com_coords\n"
              << "  Source file:   src/reactions/check_bases.cpp\n"
              << "  Function:      check_bases() -- v1/v2 built from tmpComCoord\n"
              << "  Scenario:      molecule 0 tmpComCoord moved to (0,0,0) so that\n"
              << "                 v1 = (1,0,0) is anti-parallel to sigma = (-1,0,0)\n"
              << "                 and the measured theta1 becomes pi.\n"
              << "  Pass criteria: passing theta1 = pi trips the theta1 branch\n"
              << "                 (cancelAssoc true), which can only happen if the\n"
              << "                 temporary COM was used for the measurement.\n";

    CbScene scene = cb_build_scene();

    // Move only the temporary COM; keep the tmp interface where it is so the
    // COM->interface length changes -- irrelevant here because the theta1 branch
    // returns before the conservation checks are reached.
    scene.moleculeList[0].tmpComCoord = Coord { 0.0, 0.0, 0.0 };

    ForwardRxn rxn = cb_make_rxn(/*theta1=*/M_PI, kCbTheta2Measured, kCbNaN, kCbNaN, kCbNaN);

    const bool cancelAssoc = cb_run_check(scene, rxn);
    EXPECT_TRUE(cancelAssoc)
        << "theta1 must be measured from tmpComCoord; with tmpCOM1 at the origin it is pi";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* helper is invoked from its own TEST so
// every scenario is reported individually and a failure in one does not stop the
// others.
// -----------------------------------------------------------------------------
TEST(CheckBasesTest, AcceptsConsistentGeometry) { test_check_bases_accepts_consistent_geometry(); }
TEST(CheckBasesTest, Theta1Branch) { test_check_bases_theta1_branch(); }
TEST(CheckBasesTest, Theta2MismatchCancels) { test_check_bases_theta2_mismatch_cancels(); }
TEST(CheckBasesTest, Theta2SupplementaryAngleAccepted) { test_check_bases_theta2_supplementary_angle_accepted(); }
TEST(CheckBasesTest, Phi1MismatchCancels) { test_check_bases_phi1_mismatch_cancels(); }
TEST(CheckBasesTest, Phi2MismatchCancels) { test_check_bases_phi2_mismatch_cancels(); }
TEST(CheckBasesTest, OmegaMismatchCancels) { test_check_bases_omega_mismatch_cancels(); }
TEST(CheckBasesTest, BrokenMagnitudeComplex1Cancels) { test_check_bases_broken_magnitude_complex1_cancels(); }
TEST(CheckBasesTest, BrokenMagnitudeComplex2Cancels) { test_check_bases_broken_magnitude_complex2_cancels(); }
TEST(CheckBasesTest, RigidRotationAccepted) { test_check_bases_rigid_rotation_accepted(); }
TEST(CheckBasesTest, NeverClearsCancelFlag) { test_check_bases_never_clears_cancel_flag(); }
TEST(CheckBasesTest, UsesTmpComCoords) { test_check_bases_uses_tmp_com_coords(); }