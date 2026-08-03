/*! \file test_calculate_phi.cpp
 *
 * ### Unit test for src/reactions/calculate_phi.cpp
 *
 * Function under test:
 * \code
 *   double calculate_phi(Coord reactIface1, int ifaceIndex2,
 *                        Molecule reactMol1, Molecule reactMol2,
 *                        const Vector& normal, Vector axis,
 *                        const ForwardRxn& currRxn,
 *                        const std::vector<MolTemplate>& molTemplateList);
 * \endcode
 *
 * `calculate_phi` returns the (sigma)-(interface)-(normal) dihedral angle used by
 * association.  Its algorithm is:
 *
 *   1. `transform()` rigidly rotates the temporary coordinates of both reactants so
 *      that the supplied `axis` (the COM->interface vector of reactant 1) lies on the
 *      z axis.
 *   2. sigma is formed as (reacting iface of mol 1) - (reacting iface of mol 2).
 *   3. `determine_normal()` maps the *internal frame* normal of reactant 1 into the
 *      current (transformed) frame.
 *   4. Both vectors are orthographically projected onto the xy plane (their z
 *      components are dropped) and the angle between the projections is returned,
 *      with a sign convention based on the z component of their cross product.
 *
 * ## How the expected values below are derived (no dependence on the internal
 *    rotation that step 1 happens to pick)
 *
 * Let `R` be the rotation applied by `transform()`; by construction `R * axis = z`.
 * Because a rotation preserves dot products, the following two situations give
 * analytically known answers no matter what `R` is:
 *
 *   - If the real-space normal lies **inside** the plane spanned by (axis, sigma),
 *     then after the projection the two projected vectors are (anti)parallel, so
 *     |phi| == 0 (parallel to sigma) or |phi| == PI (antiparallel to sigma).
 *   - If the real-space normal is **perpendicular** to both `axis` and `sigma`,
 *     then the projected vectors are orthogonal, so |phi| == PI/2.
 *
 * The test geometry places molecule 1 exactly in its MolTemplate orientation, so
 * `determine_normal()` must return the supplied internal-frame normal rotated by the
 * same `R`, which is what makes the reasoning above valid.
 *
 * Geometry used by every test below (pre-transform / lab frame):
 *
 *      mol1 COM        = (0, 0, 0)      mol1 iface0 (reacting) = (1, 0, 0)
 *      mol1 iface1     = (0, 1, 0)      => axis  = iface0 - COM = (1, 0, 0)
 *      mol2 COM        = (1, 2, 0)      mol2 iface0 (reacting) = (1, 1, 0)
 *      sigma           = iface1 - iface2 = (1,0,0) - (1,1,0) = (0, -1, 0)
 *
 *   => plane(axis, sigma) is the xy plane.
 *      normal (0,-1,0) is parallel      to sigma   -> phi ~   0
 *      normal (0, 1,0) is antiparallel  to sigma   -> |phi| ~ PI
 *      normal (0, 0,1) is perpendicular to both    -> |phi| ~ PI/2
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

//! Numerical tolerance used for all angle comparisons (rotations introduce round-off).
constexpr double kCphiTol = 1e-6;

/*! \brief Build a Vector and make sure its cached magnitude is valid.
 *
 * Several Vector member functions (dot_theta, normalize) rely on the cached
 * magnitude, so it is always computed explicitly here.
 */
Vector cphi_vec(double x, double y, double z)
{
    Vector v(x, y, z);
    v.calc_magnitude();
    return v;
}

/*! \brief Create a minimal but self-consistent MolTemplate.
 *
 * `determine_normal()` compares a Molecule's current interface coordinates against
 * the template's internal coordinates in order to recover the molecule's
 * orientation, so the interface list is the only part of the template that matters
 * for this test.  States/reactions are irrelevant to the geometry and are omitted.
 */
MolTemplate cphi_make_template(int typeIndex, const std::string& name, const std::vector<Coord>& ifaceCoords)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.comCoord = Coord { 0.0, 0.0, 0.0 }; // templates are always centred on the origin
    temp.isRod = false;
    temp.isPoint = false;
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.mass = 1.0;
    temp.D = Coord { 10.0, 10.0, 10.0 };
    temp.Dr = Coord { 0.1, 0.1, 0.1 };

    double maxRadius { 0.0 };
    for (std::size_t i { 0 }; i < ifaceCoords.size(); ++i) {
        Interface iface("i" + std::to_string(i), ifaceCoords[i]);
        iface.index = static_cast<int>(i);
        temp.interfaceList.push_back(iface);

        Coord c { ifaceCoords[i] };
        const double mag { std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z) };
        if (mag > maxRadius)
            maxRadius = mag;
    }
    temp.radius = maxRadius;
    return temp;
}

/*! \brief Create a Molecule whose real *and* temporary coordinates are identical.
 *
 * Association works on the tmp* coordinate set, but both sets are filled in so the
 * molecule is in a fully consistent state before the call.
 */
Molecule cphi_make_molecule(int typeIndex, int molIndex, const Coord& com, const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.molTypeIndex = typeIndex;
    mol.index = molIndex;
    mol.myComIndex = molIndex;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.comCoord = com;
    mol.tmpComCoord = com;

    for (std::size_t i { 0 }; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface(ifaceCoords[i]);
        iface.index = static_cast<int>(i);
        iface.relIndex = static_cast<int>(i);
        iface.molTypeIndex = typeIndex;
        mol.interfaceList.push_back(iface);
        mol.tmpICoords.push_back(ifaceCoords[i]);
    }
    return mol;
}

/*! \brief Everything `calculate_phi()` needs for one call. */
struct CphiConfig {
    std::vector<MolTemplate> molTemplateList {};
    Molecule mol1 {};
    Molecule mol2 {};
    Coord reactIface1 {};
    Vector axis {};
    int ifaceIndex2 { 0 };
    ForwardRxn rxn {};
};

/*! \brief Build the standard geometry documented at the top of this file.
 *
 * \param[in] decoyFirstPartnerIface when true, molecule 2 gets a far away *decoy*
 *            interface at relative index 0 and the real reacting interface is moved
 *            to relative index 1.  Used to prove `ifaceIndex2` is honoured.
 */
CphiConfig cphi_standard_config(bool decoyFirstPartnerIface = false)
{
    CphiConfig cfg;

    // Template 0 ("A") == molecule 1; template 1 ("B") == molecule 2.
    cfg.molTemplateList.push_back(
        cphi_make_template(0, "A", { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } }));
    cfg.molTemplateList.push_back(
        cphi_make_template(1, "B", { Coord { 0.0, -1.0, 0.0 }, Coord { -1.0, 0.0, 0.0 } }));

    // Molecule 1 sits exactly in its template orientation (identity rotation).
    cfg.mol1 = cphi_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 },
        { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 1.0, 0.0 } });

    if (!decoyFirstPartnerIface) {
        cfg.mol2 = cphi_make_molecule(1, 1, Coord { 1.0, 2.0, 0.0 },
            { Coord { 1.0, 1.0, 0.0 }, Coord { 0.0, 2.0, 0.0 } });
        cfg.ifaceIndex2 = 0;
    } else {
        cfg.mol2 = cphi_make_molecule(1, 1, Coord { 1.0, 2.0, 0.0 },
            { Coord { 9.0, 9.0, 9.0 }, Coord { 1.0, 1.0, 0.0 } });
        cfg.ifaceIndex2 = 1; // the *real* reacting interface
    }

    // Reacting interface of molecule 1 and the rotation axis (COM -> iface vector).
    cfg.reactIface1 = Coord { 1.0, 0.0, 0.0 };
    cfg.axis = cphi_vec(1.0, 0.0, 0.0);
    cfg.axis.normalize();

    // currRxn is not read by calculate_phi, but a realistic reaction is supplied.
    cfg.rxn.bindRadius = 1.0;
    cfg.rxn.norm1 = cphi_vec(0.0, 0.0, 1.0);
    cfg.rxn.norm2 = cphi_vec(0.0, 0.0, 1.0);
    cfg.rxn.assocAngles = ForwardRxn::Angles(M_PI / 2.0, M_PI / 2.0, M_PI / 2.0, M_PI / 2.0, M_PI);
    cfg.rxn.isSymmetric = false;
    cfg.rxn.isOnMem = false;

    return cfg;
}

/*! \brief Convenience wrapper: build the standard config and evaluate phi. */
double cphi_phi_for_normal(const Vector& internalNormal, bool decoyFirstPartnerIface = false)
{
    CphiConfig cfg { cphi_standard_config(decoyFirstPartnerIface) };
    return calculate_phi(cfg.reactIface1, cfg.ifaceIndex2, cfg.mol1, cfg.mol2, internalNormal, cfg.axis, cfg.rxn,
        cfg.molTemplateList);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: normal parallel to sigma  ->  phi == 0
// -----------------------------------------------------------------------------
void test_cphi_normal_parallel_to_sigma_is_zero()
{
    std::cerr << "\n[TEST] test_cphi_normal_parallel_to_sigma_is_zero\n"
              << "  Source file:   src/reactions/calculate_phi.cpp\n"
              << "  Function:      calculate_phi()\n"
              << "  Scenario:      internal normal (0,-1,0) is parallel to\n"
              << "                 sigma = (0,-1,0) and lies in plane(axis, sigma).\n"
              << "  Pass criteria: |phi| ~ 0 (projections are parallel).\n";

    const Vector normal { cphi_vec(0.0, -1.0, 0.0) };
    const double phi { cphi_phi_for_normal(normal) };

    std::cerr << "  calculate_phi returned phi = " << phi << " rad\n";

    // A NaN here would mean one of the projected vectors collapsed to zero length.
    EXPECT_FALSE(std::isnan(phi)) << "phi must be a finite number, not NaN";
    EXPECT_NEAR(std::abs(phi), 0.0, kCphiTol)
        << "A normal parallel to sigma must give a dihedral of 0 radians";
}

// -----------------------------------------------------------------------------
// Test 2: normal antiparallel to sigma  ->  |phi| == PI
// -----------------------------------------------------------------------------
void test_cphi_normal_antiparallel_to_sigma_is_pi()
{
    std::cerr << "\n[TEST] test_cphi_normal_antiparallel_to_sigma_is_pi\n"
              << "  Source file:   src/reactions/calculate_phi.cpp\n"
              << "  Function:      calculate_phi()\n"
              << "  Scenario:      internal normal (0,1,0) is antiparallel to\n"
              << "                 sigma = (0,-1,0) (still coplanar with the axis).\n"
              << "  Pass criteria: |phi| ~ PI.  The implementation deliberately does\n"
              << "                 NOT flip the sign at exactly PI, so the sign is\n"
              << "                 only reported, not asserted (+PI == -PI).\n";

    const Vector normal { cphi_vec(0.0, 1.0, 0.0) };
    const double phi { cphi_phi_for_normal(normal) };

    std::cerr << "  calculate_phi returned phi = " << phi << " rad (|phi| = " << std::abs(phi) << ")\n";

    EXPECT_FALSE(std::isnan(phi)) << "phi must be a finite number, not NaN";
    EXPECT_NEAR(std::abs(phi), M_PI, kCphiTol)
        << "A normal antiparallel to sigma must give a dihedral of +/- PI radians";
}

// -----------------------------------------------------------------------------
// Test 3: normal perpendicular to both axis and sigma  ->  |phi| == PI/2
// -----------------------------------------------------------------------------
void test_cphi_normal_perpendicular_is_half_pi()
{
    std::cerr << "\n[TEST] test_cphi_normal_perpendicular_is_half_pi\n"
              << "  Source file:   src/reactions/calculate_phi.cpp\n"
              << "  Function:      calculate_phi()\n"
              << "  Scenario:      internal normal (0,0,1) is perpendicular to both\n"
              << "                 axis = (1,0,0) and sigma = (0,-1,0).\n"
              << "  Pass criteria: |phi| ~ PI/2 (projections are orthogonal).\n";

    const Vector normal { cphi_vec(0.0, 0.0, 1.0) };
    const double phi { cphi_phi_for_normal(normal) };

    std::cerr << "  calculate_phi returned phi = " << phi << " rad (expected |phi| = " << M_PI / 2.0 << ")\n";

    EXPECT_FALSE(std::isnan(phi)) << "phi must be a finite number, not NaN";
    EXPECT_NEAR(std::abs(phi), M_PI / 2.0, kCphiTol)
        << "A normal perpendicular to axis and sigma must give |phi| = PI/2";
}

// -----------------------------------------------------------------------------
// Test 4: the sign convention.  Reversing the normal must flip the sign of phi
//         while leaving its magnitude alone (this exercises the
//         "if (test.z > 0) phi = -phi" branch).
// -----------------------------------------------------------------------------
void test_cphi_sign_flips_with_reversed_normal()
{
    std::cerr << "\n[TEST] test_cphi_sign_flips_with_reversed_normal\n"
              << "  Source file:   src/reactions/calculate_phi.cpp\n"
              << "  Function:      calculate_phi() (sign convention branch)\n"
              << "  Scenario:      evaluate phi for normals (0,0,+1) and (0,0,-1).\n"
              << "  Pass criteria: both magnitudes ~ PI/2, the two results have\n"
              << "                 opposite signs and sum to ~0.\n";

    const double phiPlus { cphi_phi_for_normal(cphi_vec(0.0, 0.0, 1.0)) };
    const double phiMinus { cphi_phi_for_normal(cphi_vec(0.0, 0.0, -1.0)) };

    std::cerr << "  phi(+z normal) = " << phiPlus << " rad\n"
              << "  phi(-z normal) = " << phiMinus << " rad\n"
              << "  sum            = " << (phiPlus + phiMinus) << "\n";

    EXPECT_FALSE(std::isnan(phiPlus)) << "phi for +z normal must not be NaN";
    EXPECT_FALSE(std::isnan(phiMinus)) << "phi for -z normal must not be NaN";

    EXPECT_NEAR(std::abs(phiPlus), M_PI / 2.0, kCphiTol) << "|phi| for +z normal should be PI/2";
    EXPECT_NEAR(std::abs(phiMinus), M_PI / 2.0, kCphiTol) << "|phi| for -z normal should be PI/2";

    // Opposite signs: their sum cancels and their product is negative.
    EXPECT_NEAR(phiPlus + phiMinus, 0.0, kCphiTol)
        << "Reversing the normal must flip the sign of the dihedral";
    EXPECT_LT(phiPlus * phiMinus, 0.0)
        << "phi values for opposite normals must have opposite signs";
}

// -----------------------------------------------------------------------------
// Test 5: the returned angle is always finite and inside [-PI, PI] for an
//         arbitrary (non-degenerate) normal.
// -----------------------------------------------------------------------------
void test_cphi_result_is_finite_and_in_range()
{
    std::cerr << "\n[TEST] test_cphi_result_is_finite_and_in_range\n"
              << "  Source file:   src/reactions/calculate_phi.cpp\n"
              << "  Function:      calculate_phi()\n"
              << "  Scenario:      an arbitrary skew normal (0.3,-0.7,0.5).\n"
              << "  Pass criteria: phi is finite and -PI <= phi <= PI, and the\n"
              << "                 function is deterministic (two identical calls\n"
              << "                 return bit-identical results).\n";

    const Vector normal { cphi_vec(0.3, -0.7, 0.5) };
    const double phi1 { cphi_phi_for_normal(normal) };
    const double phi2 { cphi_phi_for_normal(normal) };

    std::cerr << "  first call  phi = " << phi1 << " rad\n"
              << "  second call phi = " << phi2 << " rad\n";

    EXPECT_TRUE(std::isfinite(phi1)) << "phi must be finite (no NaN/Inf)";
    EXPECT_LE(phi1, M_PI + kCphiTol) << "phi must not exceed +PI";
    EXPECT_GE(phi1, -M_PI - kCphiTol) << "phi must not be below -PI";

    // calculate_phi uses no random numbers, so repeated calls must agree exactly.
    EXPECT_DOUBLE_EQ(phi1, phi2) << "calculate_phi must be deterministic for identical input";
}

// -----------------------------------------------------------------------------
// Test 6: `ifaceIndex2` really selects which interface of reactant 2 forms sigma.
// -----------------------------------------------------------------------------
void test_cphi_honours_partner_interface_index()
{
    std::cerr << "\n[TEST] test_cphi_honours_partner_interface_index\n"
              << "  Source file:   src/reactions/calculate_phi.cpp\n"
              << "  Function:      calculate_phi() (sigma construction)\n"
              << "  Scenario:      the same reacting partner interface (1,1,0) is\n"
              << "                 placed at relative index 0 in one run and at\n"
              << "                 relative index 1 (behind a far away decoy\n"
              << "                 interface at (9,9,9)) in the other run.\n"
              << "  Pass criteria: both runs return the same phi, proving the decoy\n"
              << "                 interface is ignored and ifaceIndex2 is used.\n";

    const Vector normal { cphi_vec(0.0, 0.0, 1.0) };
    const double phiIndex0 { cphi_phi_for_normal(normal, /*decoyFirstPartnerIface=*/false) };
    const double phiIndex1 { cphi_phi_for_normal(normal, /*decoyFirstPartnerIface=*/true) };

    std::cerr << "  phi (reacting iface at index 0) = " << phiIndex0 << " rad\n"
              << "  phi (reacting iface at index 1) = " << phiIndex1 << " rad\n";

    EXPECT_FALSE(std::isnan(phiIndex0)) << "phi must not be NaN when ifaceIndex2 == 0";
    EXPECT_FALSE(std::isnan(phiIndex1)) << "phi must not be NaN when ifaceIndex2 == 1";
    EXPECT_NEAR(phiIndex0, phiIndex1, kCphiTol)
        << "sigma must be built from the interface named by ifaceIndex2 only";
}

// -----------------------------------------------------------------------------
// Test 7: calculate_phi takes its Coord and Molecule arguments *by value*, so the
//         caller's objects must be untouched even though transform() rewrites the
//         temporary coordinates internally.
// -----------------------------------------------------------------------------
void test_cphi_does_not_modify_caller_arguments()
{
    std::cerr << "\n[TEST] test_cphi_does_not_modify_caller_arguments\n"
              << "  Source file:   src/reactions/calculate_phi.cpp\n"
              << "  Function:      calculate_phi() (pass-by-value contract)\n"
              << "  Scenario:      remember the caller's interface coordinate and the\n"
              << "                 tmp coordinates of both molecules, call the\n"
              << "                 function, then compare.\n"
              << "  Pass criteria: reactIface1, both tmpComCoords and every entry of\n"
              << "                 both tmpICoords lists are unchanged.\n";

    CphiConfig cfg { cphi_standard_config() };

    // Snapshot everything transform() would touch if the arguments were references.
    const Coord ifaceBefore { cfg.reactIface1 };
    const Coord mol1ComBefore { cfg.mol1.tmpComCoord };
    const Coord mol2ComBefore { cfg.mol2.tmpComCoord };
    const std::vector<Coord> mol1IfacesBefore { cfg.mol1.tmpICoords };
    const std::vector<Coord> mol2IfacesBefore { cfg.mol2.tmpICoords };

    const Vector normal { cphi_vec(0.0, 0.0, 1.0) };
    const double phi { calculate_phi(cfg.reactIface1, cfg.ifaceIndex2, cfg.mol1, cfg.mol2, normal, cfg.axis, cfg.rxn,
        cfg.molTemplateList) };

    std::cerr << "  calculate_phi returned phi = " << phi << " rad\n"
              << "  checking that the caller's copies were not mutated...\n";

    EXPECT_TRUE(cfg.reactIface1 == ifaceBefore)
        << "reactIface1 is a by-value parameter and must not change in the caller";
    EXPECT_TRUE(cfg.mol1.tmpComCoord == mol1ComBefore)
        << "reactMol1.tmpComCoord must not change in the caller";
    EXPECT_TRUE(cfg.mol2.tmpComCoord == mol2ComBefore)
        << "reactMol2.tmpComCoord must not change in the caller";

    ASSERT_EQ(cfg.mol1.tmpICoords.size(), mol1IfacesBefore.size())
        << "reactMol1 interface count must not change";
    for (std::size_t i { 0 }; i < mol1IfacesBefore.size(); ++i) {
        EXPECT_TRUE(cfg.mol1.tmpICoords[i] == mol1IfacesBefore[i])
            << "reactMol1.tmpICoords[" << i << "] must not change in the caller";
    }

    ASSERT_EQ(cfg.mol2.tmpICoords.size(), mol2IfacesBefore.size())
        << "reactMol2 interface count must not change";
    for (std::size_t i { 0 }; i < mol2IfacesBefore.size(); ++i) {
        EXPECT_TRUE(cfg.mol2.tmpICoords[i] == mol2IfacesBefore[i])
            << "reactMol2.tmpICoords[" << i << "] must not change in the caller";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so that a
// failure in one scenario does not stop the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(CalculatePhi, NormalParallelToSigmaIsZero) { test_cphi_normal_parallel_to_sigma_is_zero(); }
TEST(CalculatePhi, NormalAntiparallelToSigmaIsPi) { test_cphi_normal_antiparallel_to_sigma_is_pi(); }
TEST(CalculatePhi, NormalPerpendicularIsHalfPi) { test_cphi_normal_perpendicular_is_half_pi(); }
TEST(CalculatePhi, SignFlipsWithReversedNormal) { test_cphi_sign_flips_with_reversed_normal(); }
TEST(CalculatePhi, ResultIsFiniteAndInRange) { test_cphi_result_is_finite_and_in_range(); }
TEST(CalculatePhi, HonoursPartnerInterfaceIndex) { test_cphi_honours_partner_interface_index(); }
TEST(CalculatePhi, DoesNotModifyCallerArguments) { test_cphi_does_not_modify_caller_arguments(); }