/*! \file test_calculate_omega.cpp
 *
 * ### Unit test for src/reactions/calculate_omega.cpp
 *
 * Function under test:
 *
 * \code
 * double calculate_omega(Coord reactIface1, int reactIface2, Vector& sigma,
 *                        const ForwardRxn& currRxn, Molecule reactMol1,
 *                        Molecule reactMol2,
 *                        const std::vector<MolTemplate>& molTemplateList);
 * \endcode
 *
 * `calculate_omega()` computes the association dihedral angle \f$\omega\f$, i.e. the
 * angle
 *
 *     (COM1 - iface1) --- sigma --- (COM2 - iface2)
 *
 * measured *around* the sigma axis.  Internally it
 *   1. calls transform() to re-align both (temporary/association) coordinate sets so
 *      that sigma points purely along the z-axis,
 *   2. builds the two "arm" vectors
 *        - from the COM-to-interface vectors (generic case), or
 *        - from the molecule normals (when theta1 or theta2 == PI),
 *   3. orthographically projects both arms onto the xy-plane, and
 *   4. returns the angle between the projections, with the sign flipped when the
 *      cross product points along +z (guarded so that exactly 0 and exactly PI are
 *      never sign flipped).
 *
 * ### Test strategy
 *
 * The *magnitude* of omega is a rigid-body invariant: it is the angle between the
 * components of the two arm vectors that are perpendicular to sigma (the classic
 * dihedral angle).  That value can be computed independently of the internal
 * frame-alignment machinery, so every assertion below compares
 * `std::abs(calculate_omega(...))` against such an independently computed
 * reference.  This makes the tests immune to the internal sign convention of
 * transform() (whether sigma is mapped onto +z or -z), while still strongly
 * constraining the geometry.
 *
 * The *sign* behaviour is tested relatively: mirroring the second molecule's arm
 * must flip the sign of omega while preserving its magnitude.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

namespace {

// -----------------------------------------------------------------------------
// Small geometry helpers (file-local, prefixed to avoid collisions).
// -----------------------------------------------------------------------------

/*! \brief Euclidean magnitude of a Coord treated as a vector from the origin. */
double comega_mag(const Coord& c)
{
    return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
}

/*! \brief Dot product of two Coords treated as vectors. */
double comega_dot(const Coord& a, const Coord& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/*! \brief Independently computed reference dihedral angle (always in [0, PI]).
 *
 * Removes the component of each arm vector that is parallel to the sigma axis and
 * returns the angle between the two residual (perpendicular) vectors.  This is,
 * by definition, |omega|.
 *
 * \param[in] sigma The sigma axis (need not be normalized).
 * \param[in] v1    Arm vector of molecule 1 (iface1 - com1).
 * \param[in] v2    Arm vector of molecule 2 (iface2 - com2).
 */
double comega_reference_dihedral(const Coord& sigma, const Coord& v1, const Coord& v2)
{
    const double sMag = comega_mag(sigma);
    const Coord sHat { sigma.x / sMag, sigma.y / sMag, sigma.z / sMag };

    // Project out the sigma-parallel components.
    const double d1 = comega_dot(v1, sHat);
    const double d2 = comega_dot(v2, sHat);
    const Coord p1 { v1.x - d1 * sHat.x, v1.y - d1 * sHat.y, v1.z - d1 * sHat.z };
    const Coord p2 { v2.x - d2 * sHat.x, v2.y - d2 * sHat.y, v2.z - d2 * sHat.z };

    // Angle between the perpendicular residuals, clamped for numerical safety.
    double cosAng = comega_dot(p1, p2) / (comega_mag(p1) * comega_mag(p2));
    cosAng = std::max(-1.0, std::min(1.0, cosAng));
    return std::acos(cosAng);
}

/*! \brief Rigidly rotate a point: first about z by \p az, then about x by \p ax. */
Coord comega_rotate(const Coord& p, double az, double ax)
{
    const double x1 = p.x * std::cos(az) - p.y * std::sin(az);
    const double y1 = p.x * std::sin(az) + p.y * std::cos(az);
    const double z1 = p.z;

    const double x2 = x1;
    const double y2 = y1 * std::cos(ax) - z1 * std::sin(ax);
    const double z2 = y1 * std::sin(ax) + z1 * std::cos(ax);
    return Coord { x2, y2, z2 };
}

// -----------------------------------------------------------------------------
// Builders for the minimal NERDSS objects the function needs.
// -----------------------------------------------------------------------------

/*! \brief Build a single-interface Molecule with valid association (tmp) coordinates.
 *
 * calculate_omega() works exclusively on the temporary association coordinates
 * (tmpComCoord / tmpICoords), so set_tmp_association_coords() must be called.
 */
Molecule comega_make_molecule(int molTypeIndex, const Coord& com, const Coord& ifaceCoord)
{
    Molecule mol;
    mol.molTypeIndex = molTypeIndex;
    mol.index = molTypeIndex;
    mol.comCoord = com;

    mol.interfaceList.clear();
    mol.interfaceList.emplace_back(ifaceCoord); // Molecule::Iface(const Coord&)
    mol.interfaceList[0].index = 0;
    mol.interfaceList[0].relIndex = 0;
    mol.interfaceList[0].molTypeIndex = molTypeIndex;

    // Copy comCoord/interfaceList into tmpComCoord/tmpICoords.
    mol.set_tmp_association_coords();
    return mol;
}

/*! \brief Two trivial MolTemplates, only needed by the theta == PI (normal) branch. */
std::vector<MolTemplate> comega_make_templates()
{
    std::vector<MolTemplate> templates;
    for (int i = 0; i < 2; ++i) {
        MolTemplate temp;
        temp.molTypeIndex = i;
        temp.molName = (i == 0) ? "comegaA" : "comegaB";
        temp.comCoord = Coord { 0.0, 0.0, 0.0 };
        temp.radius = 1.0;
        temp.isPoint = false;
        temp.isRod = false;
        temp.D = Coord { 1.0, 1.0, 1.0 };
        temp.Dr = Coord { 0.01, 0.01, 0.01 };

        // One interface, deliberately NOT parallel to the interface vectors used in
        // the tests, so the internal template-orientation quaternion is well defined.
        Interface iface("i0", Coord { 1.0, 0.0, 0.0 });
        iface.index = 0;
        temp.interfaceList.push_back(iface);

        templates.push_back(temp);
    }
    return templates;
}

/*! \brief Build a ForwardRxn carrying the requested theta angles and unit normals. */
ForwardRxn comega_make_rxn(double theta1, double theta2)
{
    ForwardRxn rxn;
    rxn.bindRadius = 1.0;
    rxn.assocAngles = ForwardRxn::Angles(theta1, theta2, M_PI / 2.0, M_PI / 2.0, 0.0);
    rxn.norm1 = Vector(0.0, 0.0, 1.0);
    rxn.norm1.calc_magnitude();
    rxn.norm2 = Vector(0.0, 0.0, 1.0);
    rxn.norm2.calc_magnitude();
    return rxn;
}

/*! \brief Full one-shot evaluation of calculate_omega() for a given geometry.
 *
 * sigma is taken as (iface2 - iface1), which is the physical sigma vector for a
 * binding event, and is kept in that fixed orientation for every test so that
 * relative sign comparisons remain meaningful.
 */
double comega_evaluate(const Coord& com1, const Coord& iface1, const Coord& com2, const Coord& iface2,
    double theta1 = M_PI / 2.0, double theta2 = M_PI / 2.0)
{
    Molecule mol1 = comega_make_molecule(0, com1, iface1);
    Molecule mol2 = comega_make_molecule(1, com2, iface2);

    Vector sigma(iface2.x - iface1.x, iface2.y - iface1.y, iface2.z - iface1.z);
    sigma.calc_magnitude(); // dot_theta()/normalize() rely on the cached magnitude

    ForwardRxn rxn = comega_make_rxn(theta1, theta2);
    const std::vector<MolTemplate> templates = comega_make_templates();

    // reactIface2 == 0: molecule 2 owns exactly one interface.
    return calculate_omega(mol1.tmpICoords[0], 0, sigma, rxn, mol1, mol2, templates);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: orthogonal arms -> |omega| must be PI/2.
// -----------------------------------------------------------------------------
void test_comega_perpendicular_arms()
{
    std::cerr << "\n[TEST] test_comega_perpendicular_arms\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      sigma along +x, arm1 along +y, arm2 along +z.\n"
              << "  Pass criteria: |omega| == PI/2 (the two arms are perpendicular\n"
              << "                 in the plane normal to sigma).\n";

    // arm1 = iface1 - com1 = (0,1,0); arm2 = iface2 - com2 = (0,0,1)
    const Coord com1 { 0.0, 0.0, 0.0 };
    const Coord iface1 { 0.0, 1.0, 0.0 };
    const Coord iface2 { 4.0, 1.0, 0.0 }; // sigma = (4,0,0) -> +x
    const Coord com2 { 4.0, 1.0, -1.0 };

    const double omega = comega_evaluate(com1, iface1, com2, iface2);
    std::cerr << "  Returned omega = " << omega << " rad, |omega| = " << std::abs(omega)
              << " (expected " << M_PI / 2.0 << ")\n";

    EXPECT_NEAR(std::abs(omega), M_PI / 2.0, 1e-8)
        << "Perpendicular arms about sigma must give |omega| = PI/2";
    EXPECT_LE(std::abs(omega), M_PI + 1e-9) << "omega must always lie within [-PI, PI]";
}

// -----------------------------------------------------------------------------
// Test 2: anti-parallel arms -> omega must be exactly +PI (sign flip suppressed).
// -----------------------------------------------------------------------------
void test_comega_antiparallel_arms()
{
    std::cerr << "\n[TEST] test_comega_antiparallel_arms\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      arm1 = +y, arm2 = -y, sigma along +x.\n"
              << "  Pass criteria: omega == +PI. The implementation explicitly does\n"
              << "                 NOT flip the sign when |omega| is PI, so the\n"
              << "                 result must be positive PI.\n";

    const Coord com1 { 0.0, 0.0, 0.0 };
    const Coord iface1 { 0.0, 1.0, 0.0 }; // arm1 = (0,1,0)
    const Coord iface2 { 4.0, 1.0, 0.0 };
    const Coord com2 { 4.0, 2.0, 0.0 }; // arm2 = (0,-1,0)

    const double omega = comega_evaluate(com1, iface1, com2, iface2);
    std::cerr << "  Returned omega = " << omega << " rad (expected +" << M_PI << ")\n";

    EXPECT_NEAR(omega, M_PI, 1e-8) << "Anti-parallel arms must return exactly +PI (no sign flip)";
    EXPECT_GT(omega, 0.0) << "PI must be returned with a positive sign";
}

// -----------------------------------------------------------------------------
// Test 3: parallel arms -> omega must be (numerically) zero.
// -----------------------------------------------------------------------------
void test_comega_parallel_arms()
{
    std::cerr << "\n[TEST] test_comega_parallel_arms\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      arm1 = arm2 = +y, sigma along +x (eclipsed).\n"
              << "  Pass criteria: omega == 0 (and no spurious sign flip).\n";

    const Coord com1 { 0.0, 0.0, 0.0 };
    const Coord iface1 { 0.0, 1.0, 0.0 }; // arm1 = (0,1,0)
    const Coord iface2 { 4.0, 1.0, 0.0 };
    const Coord com2 { 4.0, 0.0, 0.0 }; // arm2 = (0,1,0)

    const double omega = comega_evaluate(com1, iface1, com2, iface2);
    std::cerr << "  Returned omega = " << omega << " rad (expected 0)\n";

    EXPECT_NEAR(omega, 0.0, 1e-6) << "Eclipsed (parallel) arms must give omega = 0";
}

// -----------------------------------------------------------------------------
// Test 4: sigma-parallel components of the arms must be projected out.
// -----------------------------------------------------------------------------
void test_comega_sigma_component_projected_out()
{
    std::cerr << "\n[TEST] test_comega_sigma_component_projected_out\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      arm1 = (2,1,0), arm2 = (-3,0,1), sigma along +x.\n"
              << "                 Both arms carry a large component along sigma.\n"
              << "  Pass criteria: |omega| == PI/2, i.e. only the perpendicular\n"
              << "                 components (+y and +z) determine the dihedral.\n";

    const Coord com1 { 0.0, 0.0, 0.0 };
    const Coord iface1 { 2.0, 1.0, 0.0 }; // arm1 = (2,1,0)
    const Coord iface2 { 6.0, 1.0, 0.0 }; // sigma = (4,0,0) -> +x
    const Coord com2 { 9.0, 1.0, -1.0 }; // arm2 = (-3,0,1)

    const double omega = comega_evaluate(com1, iface1, com2, iface2);
    const double reference = comega_reference_dihedral(
        Coord { 4.0, 0.0, 0.0 }, Coord { 2.0, 1.0, 0.0 }, Coord { -3.0, 0.0, 1.0 });

    std::cerr << "  Returned |omega| = " << std::abs(omega) << ", reference dihedral = " << reference
              << " (expected " << M_PI / 2.0 << ")\n";

    EXPECT_NEAR(reference, M_PI / 2.0, 1e-12) << "Sanity check of the reference computation";
    EXPECT_NEAR(std::abs(omega), reference, 1e-8)
        << "Components of the arms along sigma must not influence omega";
}

// -----------------------------------------------------------------------------
// Test 5: mirroring one arm must flip the sign of omega but keep its magnitude.
// -----------------------------------------------------------------------------
void test_comega_sign_flips_with_chirality()
{
    std::cerr << "\n[TEST] test_comega_sign_flips_with_chirality\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      identical geometry except arm2 = +z vs arm2 = -z.\n"
              << "  Pass criteria: equal magnitudes, opposite signs (the sign of\n"
              << "                 omega tracks the handedness of the dihedral).\n";

    const Coord com1 { 0.0, 0.0, 0.0 };
    const Coord iface1 { 0.0, 1.0, 0.0 }; // arm1 = (0,1,0)
    const Coord iface2 { 4.0, 1.0, 0.0 }; // sigma = +x

    // arm2 = +z  and  arm2 = -z : mirror images through the (sigma, arm1) plane.
    const double omegaPlus = comega_evaluate(com1, iface1, Coord { 4.0, 1.0, -1.0 }, iface2);
    const double omegaMinus = comega_evaluate(com1, iface1, Coord { 4.0, 1.0, 1.0 }, iface2);

    std::cerr << "  omega(arm2=+z) = " << omegaPlus << ", omega(arm2=-z) = " << omegaMinus << '\n';

    EXPECT_NEAR(std::abs(omegaPlus), std::abs(omegaMinus), 1e-8)
        << "Mirrored geometries must have the same |omega|";
    EXPECT_NEAR(omegaPlus, -omegaMinus, 1e-8)
        << "Mirrored geometries must have opposite omega signs";
    EXPECT_GT(std::abs(omegaPlus), 1e-6) << "The test geometry must not be degenerate";
}

// -----------------------------------------------------------------------------
// Test 6: sweep the dihedral through a full turn and compare with the reference.
// -----------------------------------------------------------------------------
void test_comega_full_rotation_sweep()
{
    std::cerr << "\n[TEST] test_comega_full_rotation_sweep\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      arm2 is rotated about sigma in 30 degree steps\n"
              << "                 through a full 360 degree turn.\n"
              << "  Pass criteria: |omega| equals the folded rotation angle\n"
              << "                 min(phi, 2*PI - phi) at every step, and always\n"
              << "                 stays inside [0, PI].\n";

    const Coord com1 { 0.0, 0.0, 0.0 };
    const Coord iface1 { 0.0, 1.0, 0.0 }; // arm1 = (0,1,0)
    const Coord iface2 { 4.0, 1.0, 0.0 }; // sigma = (4,0,0) -> +x

    for (int step = 0; step < 12; ++step) {
        const double phi = step * (M_PI / 6.0); // 0, 30, 60 ... 330 degrees
        // arm2 = (0, cos(phi), sin(phi)) -> dihedral with arm1 (+y) is phi.
        const Coord arm2 { 0.0, std::cos(phi), std::sin(phi) };
        const Coord com2 { iface2.x - arm2.x, iface2.y - arm2.y, iface2.z - arm2.z };

        const double omega = comega_evaluate(com1, iface1, com2, iface2);
        const double expected = std::min(phi, 2.0 * M_PI - phi); // folded into [0, PI]

        std::cerr << "    phi = " << (step * 30) << " deg -> omega = " << omega
                  << ", |omega| expected " << expected << '\n';

        EXPECT_NEAR(std::abs(omega), expected, 1e-8)
            << "|omega| must match the imposed dihedral at phi = " << (step * 30) << " degrees";
        EXPECT_LE(std::abs(omega), M_PI + 1e-9)
            << "omega must remain inside [-PI, PI] at phi = " << (step * 30) << " degrees";
    }
}

// -----------------------------------------------------------------------------
// Test 7: oblique geometry checked against the independent reference dihedral.
// -----------------------------------------------------------------------------
void test_comega_oblique_geometry_matches_reference()
{
    std::cerr << "\n[TEST] test_comega_oblique_geometry_matches_reference\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      an arbitrary, non-axis-aligned arrangement of the\n"
              << "                 two molecules and sigma.\n"
              << "  Pass criteria: |omega| equals the dihedral angle computed with an\n"
              << "                 independent projection formula.\n";

    const Coord com1 { 0.5, -0.3, 0.2 };
    const Coord iface1 { 1.0, 0.4, -0.1 }; // arm1 = (0.5, 0.7, -0.3)
    const Coord iface2 { 3.0, 1.4, -0.6 }; // sigma = (2.0, 1.0, -0.5)
    const Coord com2 { 2.2, 2.0, 0.9 }; // arm2 = (0.8, -0.6, -1.5)

    const Coord sigma { iface2.x - iface1.x, iface2.y - iface1.y, iface2.z - iface1.z };
    const Coord arm1 { iface1.x - com1.x, iface1.y - com1.y, iface1.z - com1.z };
    const Coord arm2 { iface2.x - com2.x, iface2.y - com2.y, iface2.z - com2.z };

    const double omega = comega_evaluate(com1, iface1, com2, iface2);
    const double reference = comega_reference_dihedral(sigma, arm1, arm2);

    std::cerr << "  Returned |omega| = " << std::abs(omega) << ", reference dihedral = " << reference
              << '\n';

    EXPECT_NEAR(std::abs(omega), reference, 1e-8)
        << "|omega| must equal the independently computed dihedral for oblique geometry";
}

// -----------------------------------------------------------------------------
// Test 8: |omega| must be invariant under a rigid rotation of the whole system.
// -----------------------------------------------------------------------------
void test_comega_rigid_rotation_invariance()
{
    std::cerr << "\n[TEST] test_comega_rigid_rotation_invariance\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      the same oblique geometry, then rotated rigidly\n"
              << "                 (30 deg about z followed by 40 deg about x).\n"
              << "  Pass criteria: |omega| is unchanged, proving the internal\n"
              << "                 sigma-to-z alignment works for any orientation.\n";

    const Coord com1 { 0.5, -0.3, 0.2 };
    const Coord iface1 { 1.0, 0.4, -0.1 };
    const Coord iface2 { 3.0, 1.4, -0.6 };
    const Coord com2 { 2.2, 2.0, 0.9 };

    const double omegaOriginal = comega_evaluate(com1, iface1, com2, iface2);

    // Apply one rigid rotation to every point; sigma is re-derived from the
    // rotated interfaces inside comega_evaluate(), so it rotates as well.
    const double az = 30.0 * M_PI / 180.0;
    const double ax = 40.0 * M_PI / 180.0;
    const double omegaRotated = comega_evaluate(comega_rotate(com1, az, ax), comega_rotate(iface1, az, ax),
        comega_rotate(com2, az, ax), comega_rotate(iface2, az, ax));

    std::cerr << "  |omega| original = " << std::abs(omegaOriginal)
              << ", |omega| rotated = " << std::abs(omegaRotated) << '\n';

    EXPECT_NEAR(std::abs(omegaRotated), std::abs(omegaOriginal), 1e-8)
        << "|omega| must be invariant under a rigid rotation of the whole system";
}

// -----------------------------------------------------------------------------
// Test 9: the theta == PI branch (omega derived from the molecule normals).
// -----------------------------------------------------------------------------
void test_comega_normal_branch_when_theta_is_pi()
{
    std::cerr << "\n[TEST] test_comega_normal_branch_when_theta_is_pi\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega() (theta == PI code path)\n"
              << "  Scenario:      theta1 = PI, so the COM-interface vectors are\n"
              << "                 parallel to sigma and determine_normal() is used\n"
              << "                 to build the two arm vectors instead.\n"
              << "  Pass criteria: the call completes and returns a finite angle in\n"
              << "                 [-PI, PI] (the normal-based value itself depends on\n"
              << "                 the MolTemplate frames, so only the contract of the\n"
              << "                 return value is asserted).\n";

    // Geometry in which both arms are (anti)parallel to sigma, the situation for
    // which theta == PI is the physically meaningful description.
    const Coord com1 { -1.0, 0.0, 0.0 };
    const Coord iface1 { 0.0, 0.0, 0.0 }; // arm1 = (1,0,0)  || sigma
    const Coord iface2 { 4.0, 0.0, 0.0 }; // sigma = (4,0,0)
    const Coord com2 { 5.0, 0.0, 0.0 }; // arm2 = (-1,0,0) || sigma

    const double omega = comega_evaluate(com1, iface1, com2, iface2, M_PI, M_PI / 2.0);
    std::cerr << "  Returned omega (normal branch) = " << omega << " rad\n";

    const bool isFinite = std::isfinite(omega);
    EXPECT_TRUE(isFinite) << "The theta == PI branch must return a finite angle";
    if (isFinite) {
        EXPECT_LE(std::abs(omega), M_PI + 1e-9)
            << "The normal-based omega must still lie within [-PI, PI]";
    }
}

// -----------------------------------------------------------------------------
// Test 10: the caller's molecules / interface coordinate must not be modified.
// -----------------------------------------------------------------------------
void test_comega_inputs_are_not_modified()
{
    std::cerr << "\n[TEST] test_comega_inputs_are_not_modified\n"
              << "  Source file:   src/reactions/calculate_omega.cpp\n"
              << "  Function:      calculate_omega()\n"
              << "  Scenario:      calculate_omega() takes reactIface1, reactMol1 and\n"
              << "                 reactMol2 BY VALUE, even though it internally calls\n"
              << "                 transform() which rewrites the tmp coordinates.\n"
              << "  Pass criteria: the caller's Coord and both Molecules retain their\n"
              << "                 original tmp coordinates after the call.\n";

    Molecule mol1 = comega_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 });
    Molecule mol2 = comega_make_molecule(1, Coord { 4.0, 1.0, -1.0 }, Coord { 4.0, 1.0, 0.0 });

    // Remember the pre-call state.
    const Coord com1Before = mol1.tmpComCoord;
    const Coord iface1Before = mol1.tmpICoords[0];
    const Coord com2Before = mol2.tmpComCoord;
    const Coord iface2Before = mol2.tmpICoords[0];

    Coord reactIface1 = mol1.tmpICoords[0];
    const Coord reactIface1Before = reactIface1;

    Vector sigma(iface2Before.x - iface1Before.x, iface2Before.y - iface1Before.y,
        iface2Before.z - iface1Before.z);
    sigma.calc_magnitude();

    ForwardRxn rxn = comega_make_rxn(M_PI / 2.0, M_PI / 2.0);
    const std::vector<MolTemplate> templates = comega_make_templates();

    const double omega = calculate_omega(reactIface1, 0, sigma, rxn, mol1, mol2, templates);
    std::cerr << "  Returned omega = " << omega << " rad; now verifying inputs are untouched.\n";

    // The local Coord handed in by value must be unchanged.
    EXPECT_TRUE(reactIface1 == reactIface1Before)
        << "reactIface1 is passed by value and must not be modified in the caller";

    // Both molecules' temporary association coordinates must be unchanged.
    EXPECT_TRUE(mol1.tmpComCoord == com1Before) << "mol1.tmpComCoord must be unchanged";
    EXPECT_TRUE(mol1.tmpICoords[0] == iface1Before) << "mol1.tmpICoords[0] must be unchanged";
    EXPECT_TRUE(mol2.tmpComCoord == com2Before) << "mol2.tmpComCoord must be unchanged";
    EXPECT_TRUE(mol2.tmpICoords[0] == iface2Before) << "mol2.tmpICoords[0] must be unchanged";

    // Sanity: the returned angle is still the expected PI/2 for this geometry.
    EXPECT_NEAR(std::abs(omega), M_PI / 2.0, 1e-8)
        << "The by-value call must still produce the correct |omega|";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: one TEST per named test_* function so every scenario is
// reported (and executed) independently, using only non-fatal expectations.
// -----------------------------------------------------------------------------
TEST(CalculateOmegaTest, PerpendicularArms) { test_comega_perpendicular_arms(); }
TEST(CalculateOmegaTest, AntiparallelArms) { test_comega_antiparallel_arms(); }
TEST(CalculateOmegaTest, ParallelArms) { test_comega_parallel_arms(); }
TEST(CalculateOmegaTest, SigmaComponentProjectedOut) { test_comega_sigma_component_projected_out(); }
TEST(CalculateOmegaTest, SignFlipsWithChirality) { test_comega_sign_flips_with_chirality(); }
TEST(CalculateOmegaTest, FullRotationSweep) { test_comega_full_rotation_sweep(); }
TEST(CalculateOmegaTest, ObliqueGeometryMatchesReference) { test_comega_oblique_geometry_matches_reference(); }
TEST(CalculateOmegaTest, RigidRotationInvariance) { test_comega_rigid_rotation_invariance(); }
TEST(CalculateOmegaTest, NormalBranchWhenThetaIsPi) { test_comega_normal_branch_when_theta_is_pi(); }
TEST(CalculateOmegaTest, InputsAreNotModified) { test_comega_inputs_are_not_modified(); }