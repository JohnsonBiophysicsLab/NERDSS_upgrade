/*! \file test_omega_rotation.cpp
 *
 * ### Unit test for src/reactions/omega_rotation.cpp
 *
 * Function under test:
 *
 *     void omega_rotation(Coord& reactIface1, Coord& reactIface2, int ifaceIndex2,
 *                         Molecule& reactMol1, Molecule& reactMol2,
 *                         Complex& reactCom1, Complex& reactCom2, double targOmega,
 *                         const ForwardRxn& currRxn,
 *                         std::vector<Molecule>& moleculeList,
 *                         const std::vector<MolTemplate>& molTemplateList)
 *
 * `omega_rotation()` is one of the four rigid-body rotations performed during
 * association.  It rotates the two associating complexes about the "sigma" axis
 * (the vector joining the two reacting interfaces) until the
 * (COM1)-(iface1)-(sigma)-(iface2)-(COM2) dihedral angle - omega - matches the
 * value requested by the reaction.
 *
 * All of the geometry manipulated here lives in the *temporary* association
 * coordinates (`Molecule::tmpComCoord` / `Molecule::tmpICoords`), so the test
 * builds a minimal but fully self-consistent association scenario:
 *
 *   complex 0 : molecule 0 (reactant 1) + molecule 2 (an innocent bystander)
 *   complex 1 : molecule 1 (reactant 2)
 *
 *   reacting interface of molecule 0 : (0,0,0)
 *   reacting interface of molecule 1 : (1,0,0)      -> sigma is the x-axis, |sigma| = 1
 *
 *   COM of molecule 0 : (0,-1,0)  -> (COM-iface) projects onto -y
 *   COM of molecule 1 : (1,0,-1)  -> (COM-iface) projects onto -z
 *
 * The two projections are perpendicular, so the initial |omega| is pi/2.  From
 * there we ask for a variety of target omegas and check the outcome.
 *
 * Things that are verified:
 *   1. The scenario really starts at |omega| = pi/2 (sanity of the fixture).
 *   2. Asking for the angle the system already has performs *no* rotation.
 *   3. Asking for omega = pi converges to +/-pi.
 *   4. Asking for omega = 0 converges to 0.
 *   5. Asking for a generic angle (pi/4) converges to that magnitude.
 *   6. The rotation is rigid: intra-complex distances and every
 *      |COM - interface| magnitude are conserved.
 *   7. The two reacting interfaces lie *on* the rotation axis, so they must not
 *      move at all, and the sigma separation must be conserved.
 *   8. The rotation actually moves molecules (i.e. something happened).
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "reactions/association/association.hpp"

namespace {

// -----------------------------------------------------------------------------
// Small geometry helpers (prefixed to stay unique inside the whole test suite)
// -----------------------------------------------------------------------------

/*! \brief Euclidean distance between two Coords. */
double omegarot_dist(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Container holding a complete, self-consistent association scenario. */
struct OmegaRotScenario {
    std::vector<Molecule> moleculeList; //!< molecules 0,1 react; molecule 2 rides along
    std::vector<Complex> complexList; //!< complex 0 = {mol 0, mol 2}, complex 1 = {mol 1}
    std::vector<MolTemplate> molTemplateList; //!< a single template shared by all molecules
    ForwardRxn rxn; //!< the association reaction being performed
};

/*! \brief Build a molecule with the given COM and interface coordinates.
 *
 * The temporary association coordinates are initialised from the "real"
 * coordinates, because omega_rotation() (like the rest of association) operates
 * exclusively on the temporary coordinates.
 */
Molecule omegarot_make_molecule(int index, int comIndex, const Coord& com,
    const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0; // everybody uses molTemplateList[0]
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.comCoord = com;

    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface;
        iface.coord = ifaceCoords[i];
        iface.index = static_cast<int>(i); // absolute index (irrelevant here)
        iface.relIndex = static_cast<int>(i); // relative index inside the molecule
        iface.stateIndex = 0;
        iface.molTypeIndex = 0;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
    }

    // Copies comCoord -> tmpComCoord and interfaceList -> tmpICoords.
    mol.set_tmp_association_coords();
    return mol;
}

/*! \brief Build a complex owning the listed member molecules. */
Complex omegarot_make_complex(int index, const Coord& com, const std::vector<int>& members)
{
    Complex targCom;
    targCom.index = index;
    targCom.comCoord = com;
    targCom.tmpComCoord = com;
    targCom.memberList = members;
    targCom.mass = static_cast<double>(members.size());
    targCom.radius = 3.0;

    // Both complexes get identical, non-zero diffusion constants so that
    // determine_rotation_angles() splits the required rotation between them.
    targCom.D = Coord(10.0, 10.0, 10.0);
    targCom.Dr = Coord(0.1, 0.1, 0.1);
    targCom.OnSurface = false;
    return targCom;
}

/*! \brief Build the single MolTemplate shared by every molecule of the scenario.
 *
 * Internal (template) coordinates: COM at the origin, interface 0 at (0,1,0) and
 * interface 1 at (1,0,0).  Every molecule created below is an exact rigid
 * rotation of this template, so any code that orients a molecule onto its
 * template (e.g. determine_normal()) stays well behaved.
 */
MolTemplate omegarot_make_template()
{
    MolTemplate molTemplate;
    molTemplate.molName = "testMol";
    molTemplate.molTypeIndex = 0;
    molTemplate.comCoord = Coord(0.0, 0.0, 0.0);
    molTemplate.mass = 1.0;
    molTemplate.radius = 1.0;
    molTemplate.D = Coord(10.0, 10.0, 10.0);
    molTemplate.Dr = Coord(0.1, 0.1, 0.1);
    molTemplate.isRod = false; // NOT a rod: omega uses the COM-iface vectors
    molTemplate.isPoint = false; // NOT a point: COM and ifaces are distinct
    molTemplate.isLipid = false;
    molTemplate.isImplicitLipid = false;
    molTemplate.checkOverlap = false;
    molTemplate.copies = 3;

    const std::vector<Coord> internalCrds { Coord(0.0, 1.0, 0.0), Coord(1.0, 0.0, 0.0) };
    for (std::size_t i = 0; i < internalCrds.size(); ++i) {
        Interface iface(std::string("i") + std::to_string(i), internalCrds[i]);
        iface.index = static_cast<int>(i);

        // One (unnamed) state per interface, created without invoking any
        // side-effect-carrying constructor.
        Interface::State state;
        state.index = static_cast<int>(i);
        state.iden = '\0';
        iface.stateList.push_back(state);

        molTemplate.interfaceList.push_back(iface);
    }
    return molTemplate;
}

/*! \brief Assemble the whole two-complex association scenario described above. */
OmegaRotScenario omegarot_build_scenario()
{
    OmegaRotScenario scenario;

    scenario.molTemplateList.push_back(omegarot_make_template());

    // Molecule 0 (reactant 1): identity orientation of the template.
    scenario.moleculeList.push_back(omegarot_make_molecule(0, 0, Coord(0.0, -1.0, 0.0),
        { Coord(0.0, 0.0, 0.0), Coord(1.0, -1.0, 0.0) }));

    // Molecule 1 (reactant 2): template rotated +90 deg about x, so its
    // (COM - iface) vector points along -z instead of -y.
    scenario.moleculeList.push_back(omegarot_make_molecule(1, 1, Coord(1.0, 0.0, -1.0),
        { Coord(1.0, 0.0, 0.0), Coord(2.0, 0.0, -1.0) }));

    // Molecule 2: a passive member of complex 0 used to verify rigid motion.
    scenario.moleculeList.push_back(omegarot_make_molecule(2, 0, Coord(0.0, -3.0, 0.0),
        { Coord(0.0, -2.0, 0.0), Coord(1.0, -3.0, 0.0) }));

    scenario.complexList.push_back(omegarot_make_complex(0, Coord(0.0, -2.0, 0.0), { 0, 2 }));
    scenario.complexList.push_back(omegarot_make_complex(1, Coord(1.0, 0.0, -1.0), { 1 }));

    // The reaction itself.  theta1/theta2 are pi/2 (not 0 or pi), which keeps
    // calculate_omega() on the "use the COM-iface vectors" branch.
    scenario.rxn.rxnType = ReactionType::bimolecular;
    scenario.rxn.bindRadius = 1.0;
    scenario.rxn.norm1 = Vector(0.0, 0.0, 1.0);
    scenario.rxn.norm2 = Vector(0.0, 0.0, 1.0);
    scenario.rxn.assocAngles = ForwardRxn::Angles(M_PI / 2.0, M_PI / 2.0, M_PI / 2.0, M_PI / 2.0, M_PI);
    scenario.rxn.reactantListNew.emplace_back("i0", 0, 0, 0, '\0', false);
    scenario.rxn.reactantListNew.emplace_back("i0", 0, 0, 0, '\0', false);

    return scenario;
}

/*! \brief Compute the current omega of the scenario exactly the way
 *         omega_rotation() does internally (sigma = iface1 - iface2, normalised).
 */
double omegarot_measure_omega(const OmegaRotScenario& scenario)
{
    const Coord& iface1 = scenario.moleculeList[0].tmpICoords[0];
    const Coord& iface2 = scenario.moleculeList[1].tmpICoords[0];
    Vector sigma { iface1 - iface2 };
    sigma.normalize();
    return calculate_omega(iface1, 0, sigma, scenario.rxn, scenario.moleculeList[0],
        scenario.moleculeList[1], scenario.molTemplateList);
}

/*! \brief Invoke the function under test with the requested target angle. */
void omegarot_run(OmegaRotScenario& scenario, double targOmega)
{
    omega_rotation(scenario.moleculeList[0].tmpICoords[0], // reactIface1
        scenario.moleculeList[1].tmpICoords[0], // reactIface2
        0, // ifaceIndex2
        scenario.moleculeList[0], scenario.moleculeList[1],
        scenario.complexList[0], scenario.complexList[1],
        targOmega, scenario.rxn, scenario.moleculeList, scenario.molTemplateList);
}

/*! \brief Flatten every temporary coordinate (COM + interfaces) into one list. */
std::vector<Coord> omegarot_snapshot(const std::vector<Molecule>& moleculeList)
{
    std::vector<Coord> snapshot;
    for (const auto& mol : moleculeList) {
        snapshot.push_back(mol.tmpComCoord);
        for (const auto& iface : mol.tmpICoords)
            snapshot.push_back(iface);
    }
    return snapshot;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the fixture itself - the scenario must start at |omega| = pi/2.
// -----------------------------------------------------------------------------
void test_omegarot_initial_geometry_sanity()
{
    std::cerr << "\n[TEST] test_omegarot_initial_geometry_sanity\n"
              << "  Source file:   src/reactions/omega_rotation.cpp\n"
              << "  Helper tested: calculate_omega() on the freshly built scenario\n"
              << "  Scenario:      sigma along +x with |sigma| = 1, the two\n"
              << "                 (COM - iface) vectors pointing along -y and -z.\n"
              << "  Pass criteria: |omega| == pi/2 and |sigma| == bindRadius (1.0).\n";

    OmegaRotScenario scenario = omegarot_build_scenario();

    const double sigmaLen = omegarot_dist(scenario.moleculeList[0].tmpICoords[0],
        scenario.moleculeList[1].tmpICoords[0]);
    std::cerr << "  |sigma| = " << sigmaLen << " (expected 1.0)\n";
    EXPECT_NEAR(sigmaLen, 1.0, 1e-12)
        << "The two reacting interfaces should be one binding radius apart";

    const double omega = omegarot_measure_omega(scenario);
    std::cerr << "  initial omega = " << omega << " rad (|omega| expected "
              << M_PI / 2.0 << ")\n";
    EXPECT_NEAR(std::abs(omega), M_PI / 2.0, 1e-6)
        << "Perpendicular COM-iface projections must give |omega| = pi/2";
}

// -----------------------------------------------------------------------------
// Test 2: requesting the angle the system already has must be a no-op.
// -----------------------------------------------------------------------------
void test_omegarot_no_rotation_when_target_matches()
{
    std::cerr << "\n[TEST] test_omegarot_no_rotation_when_target_matches\n"
              << "  Function:      omega_rotation()\n"
              << "  Scenario:      targOmega is set to the *current* omega.\n"
              << "  Pass criteria: the early-return branch fires, so every\n"
              << "                 temporary coordinate is bit-for-bit unchanged.\n";

    OmegaRotScenario scenario = omegarot_build_scenario();

    const double currOmega = omegarot_measure_omega(scenario);
    std::cerr << "  current omega = " << currOmega << " -> requesting the same value\n";

    const std::vector<Coord> before = omegarot_snapshot(scenario.moleculeList);
    omegarot_run(scenario, currOmega);
    const std::vector<Coord> after = omegarot_snapshot(scenario.moleculeList);

    ASSERT_EQ(before.size(), after.size()) << "Snapshot sizes must match";
    for (std::size_t i = 0; i < before.size(); ++i) {
        EXPECT_DOUBLE_EQ(after[i].x, before[i].x) << "coord " << i << " x changed";
        EXPECT_DOUBLE_EQ(after[i].y, before[i].y) << "coord " << i << " y changed";
        EXPECT_DOUBLE_EQ(after[i].z, before[i].z) << "coord " << i << " z changed";
    }
    std::cerr << "  All " << before.size() << " temporary coordinates were left untouched.\n";
}

// -----------------------------------------------------------------------------
// Test 3: converge onto omega = pi (the "anti-parallel" special case).
// -----------------------------------------------------------------------------
void test_omegarot_rotates_to_pi()
{
    std::cerr << "\n[TEST] test_omegarot_rotates_to_pi\n"
              << "  Function:      omega_rotation()\n"
              << "  Scenario:      start at |omega| = pi/2, request omega = pi.\n"
              << "  Pass criteria: recomputed |omega| == pi (either sign is\n"
              << "                 accepted, since +pi and -pi are the same dihedral).\n";

    OmegaRotScenario scenario = omegarot_build_scenario();
    std::cerr << "  omega before = " << omegarot_measure_omega(scenario) << '\n';

    omegarot_run(scenario, M_PI);

    const double finalOmega = omegarot_measure_omega(scenario);
    std::cerr << "  omega after  = " << finalOmega << " (target " << M_PI << ")\n";
    EXPECT_NEAR(std::abs(finalOmega), M_PI, 1e-6)
        << "omega_rotation() should have driven omega to +/-pi";
}

// -----------------------------------------------------------------------------
// Test 4: converge onto omega = 0 (the "parallel" special case).
// -----------------------------------------------------------------------------
void test_omegarot_rotates_to_zero()
{
    std::cerr << "\n[TEST] test_omegarot_rotates_to_zero\n"
              << "  Function:      omega_rotation()\n"
              << "  Scenario:      start at |omega| = pi/2, request omega = 0.\n"
              << "  Pass criteria: recomputed |omega| == 0.\n";

    OmegaRotScenario scenario = omegarot_build_scenario();
    std::cerr << "  omega before = " << omegarot_measure_omega(scenario) << '\n';

    omegarot_run(scenario, 0.0);

    const double finalOmega = omegarot_measure_omega(scenario);
    std::cerr << "  omega after  = " << finalOmega << " (target 0)\n";
    EXPECT_NEAR(std::abs(finalOmega), 0.0, 1e-6)
        << "omega_rotation() should have driven omega to 0";
}

// -----------------------------------------------------------------------------
// Test 5: converge onto a generic (non special-cased) target angle.
// -----------------------------------------------------------------------------
void test_omegarot_rotates_to_generic_angle()
{
    std::cerr << "\n[TEST] test_omegarot_rotates_to_generic_angle\n"
              << "  Function:      omega_rotation()\n"
              << "  Scenario:      start at |omega| = pi/2, request omega = pi/4.\n"
              << "                 This exercises the generic branch, including the\n"
              << "                 'rotated the wrong way -> reverse and retry' logic.\n"
              << "  Pass criteria: recomputed |omega| == pi/4.\n";

    OmegaRotScenario scenario = omegarot_build_scenario();
    const double target = M_PI / 4.0;
    std::cerr << "  omega before = " << omegarot_measure_omega(scenario) << '\n';

    omegarot_run(scenario, target);

    const double finalOmega = omegarot_measure_omega(scenario);
    std::cerr << "  omega after  = " << finalOmega << " (target " << target << ")\n";
    EXPECT_NEAR(std::abs(finalOmega), target, 1e-6)
        << "omega_rotation() should have driven |omega| to the requested pi/4";
}

// -----------------------------------------------------------------------------
// Test 6: the motion has to be rigid for both complexes.
// -----------------------------------------------------------------------------
void test_omegarot_preserves_rigid_geometry()
{
    std::cerr << "\n[TEST] test_omegarot_preserves_rigid_geometry\n"
              << "  Function:      omega_rotation()\n"
              << "  Scenario:      rotate to omega = pi and re-measure internal\n"
              << "                 distances of both complexes.\n"
              << "  Pass criteria: every |COM - interface| magnitude and the\n"
              << "                 intra-complex COM-COM distance are conserved.\n";

    OmegaRotScenario scenario = omegarot_build_scenario();

    // Record every |COM - iface| magnitude and the complex-0 internal distance.
    std::vector<double> magsBefore;
    for (const auto& mol : scenario.moleculeList) {
        for (const auto& iface : mol.tmpICoords)
            magsBefore.push_back(omegarot_dist(mol.tmpComCoord, iface));
    }
    const double intraBefore = omegarot_dist(scenario.moleculeList[0].tmpComCoord,
        scenario.moleculeList[2].tmpComCoord);

    omegarot_run(scenario, M_PI);

    std::vector<double> magsAfter;
    for (const auto& mol : scenario.moleculeList) {
        for (const auto& iface : mol.tmpICoords)
            magsAfter.push_back(omegarot_dist(mol.tmpComCoord, iface));
    }
    const double intraAfter = omegarot_dist(scenario.moleculeList[0].tmpComCoord,
        scenario.moleculeList[2].tmpComCoord);

    ASSERT_EQ(magsBefore.size(), magsAfter.size())
        << "The number of interfaces must not change";
    for (std::size_t i = 0; i < magsBefore.size(); ++i) {
        std::cerr << "  |COM-iface| #" << i << ": " << magsBefore[i] << " -> "
                  << magsAfter[i] << '\n';
        EXPECT_NEAR(magsAfter[i], magsBefore[i], 1e-9)
            << "COM-interface magnitude #" << i << " was not conserved";
    }

    std::cerr << "  complex-0 internal COM-COM distance: " << intraBefore << " -> "
              << intraAfter << '\n';
    EXPECT_NEAR(intraAfter, intraBefore, 1e-9)
        << "Molecules inside the same complex must move rigidly together";
}

// -----------------------------------------------------------------------------
// Test 7: the reacting interfaces sit on the rotation axis and must not move.
// -----------------------------------------------------------------------------
void test_omegarot_reacting_interfaces_are_fixed()
{
    std::cerr << "\n[TEST] test_omegarot_reacting_interfaces_are_fixed\n"
              << "  Function:      omega_rotation()\n"
              << "  Scenario:      rotate to omega = pi.  The rotation axis is the\n"
              << "                 line through both reacting interfaces, so they\n"
              << "                 are invariant under the rotation.\n"
              << "  Pass criteria: both reacting interface coordinates, and the\n"
              << "                 sigma separation between them, are unchanged.\n";

    OmegaRotScenario scenario = omegarot_build_scenario();

    const Coord iface1Before = scenario.moleculeList[0].tmpICoords[0];
    const Coord iface2Before = scenario.moleculeList[1].tmpICoords[0];
    const double sigmaBefore = omegarot_dist(iface1Before, iface2Before);

    omegarot_run(scenario, M_PI);

    const Coord iface1After = scenario.moleculeList[0].tmpICoords[0];
    const Coord iface2After = scenario.moleculeList[1].tmpICoords[0];
    const double sigmaAfter = omegarot_dist(iface1After, iface2After);

    std::cerr << "  iface1: " << iface1Before << " -> " << iface1After << '\n';
    std::cerr << "  iface2: " << iface2Before << " -> " << iface2After << '\n';
    std::cerr << "  |sigma|: " << sigmaBefore << " -> " << sigmaAfter << '\n';

    EXPECT_NEAR(iface1After.x, iface1Before.x, 1e-9) << "reacting iface 1 x moved";
    EXPECT_NEAR(iface1After.y, iface1Before.y, 1e-9) << "reacting iface 1 y moved";
    EXPECT_NEAR(iface1After.z, iface1Before.z, 1e-9) << "reacting iface 1 z moved";

    EXPECT_NEAR(iface2After.x, iface2Before.x, 1e-9) << "reacting iface 2 x moved";
    EXPECT_NEAR(iface2After.y, iface2Before.y, 1e-9) << "reacting iface 2 y moved";
    EXPECT_NEAR(iface2After.z, iface2Before.z, 1e-9) << "reacting iface 2 z moved";

    EXPECT_NEAR(sigmaAfter, sigmaBefore, 1e-9)
        << "The binding separation sigma must be conserved by an omega rotation";
}

// -----------------------------------------------------------------------------
// Test 8: something actually moved (guards against a silent no-op).
// -----------------------------------------------------------------------------
void test_omegarot_actually_moves_molecules()
{
    std::cerr << "\n[TEST] test_omegarot_actually_moves_molecules\n"
              << "  Function:      omega_rotation()\n"
              << "  Scenario:      rotate from |omega| = pi/2 to omega = pi.\n"
              << "  Pass criteria: at least one temporary COM coordinate changed,\n"
              << "                 proving the rotation really was applied.\n";

    OmegaRotScenario scenario = omegarot_build_scenario();

    const Coord com0Before = scenario.moleculeList[0].tmpComCoord;
    const Coord com1Before = scenario.moleculeList[1].tmpComCoord;
    const Coord com2Before = scenario.moleculeList[2].tmpComCoord;

    omegarot_run(scenario, M_PI);

    const double moved0 = omegarot_dist(com0Before, scenario.moleculeList[0].tmpComCoord);
    const double moved1 = omegarot_dist(com1Before, scenario.moleculeList[1].tmpComCoord);
    const double moved2 = omegarot_dist(com2Before, scenario.moleculeList[2].tmpComCoord);

    std::cerr << "  displacement of molecule 0 COM (complex 0): " << moved0 << '\n';
    std::cerr << "  displacement of molecule 1 COM (complex 1): " << moved1 << '\n';
    std::cerr << "  displacement of molecule 2 COM (complex 0): " << moved2 << '\n';

    const bool anythingMoved = (moved0 > 1e-9) || (moved1 > 1e-9) || (moved2 > 1e-9);
    EXPECT_TRUE(anythingMoved)
        << "omega_rotation() was asked for a different angle but nothing moved";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario is reported separately, and every
// assertion is non-fatal so that all of the tests run even on failure.
// -----------------------------------------------------------------------------
TEST(OmegaRotation, InitialGeometrySanity) { test_omegarot_initial_geometry_sanity(); }
TEST(OmegaRotation, NoRotationWhenTargetMatches) { test_omegarot_no_rotation_when_target_matches(); }
TEST(OmegaRotation, RotatesToPi) { test_omegarot_rotates_to_pi(); }
TEST(OmegaRotation, RotatesToZero) { test_omegarot_rotates_to_zero(); }
TEST(OmegaRotation, RotatesToGenericAngle) { test_omegarot_rotates_to_generic_angle(); }
TEST(OmegaRotation, PreservesRigidGeometry) { test_omegarot_preserves_rigid_geometry(); }
TEST(OmegaRotation, ReactingInterfacesAreFixed) { test_omegarot_reacting_interfaces_are_fixed(); }
TEST(OmegaRotation, ActuallyMovesMolecules) { test_omegarot_actually_moves_molecules(); }