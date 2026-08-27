/*! \file test_calc_one_angular_displacement.cpp
 *
 * ### Unit test for src/reactions/calc_one_angular_displacement.cpp
 *
 * Function under test:
 *
 *     double calc_one_angular_displacement(int ifaceIndex1,
 *                                          Molecule& reactMol1,
 *                                          Complex& reactCom1);
 *
 * The routine measures how far the (interface - centre-of-mass) vector of a
 * molecule has rotated between its *original* coordinates
 * (`interfaceList[i].coord - comCoord`) and its *temporary association*
 * coordinates (`tmpICoords[i] - tmpComCoord`).
 *
 * Implementation details that the assertions below are written against:
 *
 *   - v1 = tmpICoords[i]           - tmpComCoord      (final / temporary vector)
 *   - v2 = interfaceList[i].coord  - comCoord         (original vector)
 *   - If |v1| < 1E-12 the function returns exactly 0.0 (a "point" particle,
 *     the interface sits on the COM, so there is nothing to rotate).
 *   - Otherwise the angle is `v2.dot_theta(v1)`, i.e. acos of the normalised
 *     dot product, which is always in [0, pi].
 *   - The sign is then *flipped* (theta -> -theta) when the z component of
 *     the normalised cross product v2 x v1 is strictly positive AND the angle
 *     is neither ~0 nor ~pi.  Note that only the **z** component is inspected,
 *     so a rotation about the x or y axis is never sign-flipped.
 *   - The Complex argument is completely unused by the implementation.
 *
 * The test only builds fully initialised Molecule/Complex objects; nothing here
 * touches a code path that calls exit()/abort().
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with "coad_" so they cannot clash with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully initialised Molecule with a set of interfaces.
 *
 * \param[in] com            centre of mass in the "original" coordinate set
 * \param[in] ifaceCoords    absolute interface coordinates, original set
 * \param[in] tmpCom         centre of mass in the "temporary" (association) set
 * \param[in] tmpIfaceCoords absolute interface coordinates, temporary set
 *
 * Both coordinate sets are stored, because calc_one_angular_displacement()
 * compares one against the other.
 */
Molecule coad_make_molecule(const Coord& com, const std::vector<Coord>& ifaceCoords,
    const Coord& tmpCom, const std::vector<Coord>& tmpIfaceCoords)
{
    Molecule mol;

    // Bookkeeping fields - filled in so nothing is left uninitialised.
    mol.index = 0;
    mol.myComIndex = 0;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isEmpty = false;

    // "Real" coordinates.
    mol.comCoord = com;
    for (size_t i = 0; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface;
        iface.coord = ifaceCoords[i];
        iface.relIndex = static_cast<int>(i);
        iface.index = static_cast<int>(i);
        iface.molTypeIndex = 0;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
    }

    // Temporary association coordinates.
    mol.tmpComCoord = tmpCom;
    mol.tmpICoords = tmpIfaceCoords;

    return mol;
}

/*! \brief Build a minimal, fully initialised Complex owning molecule index 0.
 *
 * The function under test never reads this object, but it must be a valid
 * reference, so we still populate the fields a Complex normally carries.
 */
Complex coad_make_complex(const Coord& com)
{
    Complex targCom;
    targCom.comCoord = com;
    targCom.index = 0;
    targCom.mass = 1.0;
    targCom.radius = 1.0;
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Coord { 0.01, 0.01, 0.01 };
    targCom.memberList.push_back(0);
    targCom.numEachMol.push_back(1);
    targCom.isEmpty = false;
    targCom.OnSurface = false;
    return targCom;
}

//! Tolerance used for all angular comparisons (radians).
constexpr double coadTol = 1e-12;

} // namespace

// -----------------------------------------------------------------------------
// Test 1: point particle -> the temporary interface sits exactly on the
//         temporary COM, so |v1| == 0 and the function must short-circuit to 0.
// -----------------------------------------------------------------------------
void test_coad_point_particle_returns_zero()
{
    std::cerr << "\n[TEST] test_coad_point_particle_returns_zero\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Function:    calc_one_angular_displacement\n"
              << "  Scenario:    tmpICoords[0] == tmpComCoord, so the temporary\n"
              << "               interface-to-COM vector has zero length.\n"
              << "  Criteria:    the early-return branch fires and yields exactly 0.0.\n";

    // Original vector is (1,0,0); temporary vector is the zero vector.
    Molecule mol = coad_make_molecule(
        /*com*/ Coord { 0.0, 0.0, 0.0 },
        /*ifaceCoords*/ { Coord { 1.0, 0.0, 0.0 } },
        /*tmpCom*/ Coord { 5.0, 5.0, 5.0 },
        /*tmpIfaceCoords*/ { Coord { 5.0, 5.0, 5.0 } });
    Complex com = coad_make_complex(Coord { 0.0, 0.0, 0.0 });

    const double theta = calc_one_angular_displacement(0, mol, com);
    std::cerr << "  Returned angle = " << theta << " rad\n";

    EXPECT_DOUBLE_EQ(theta, 0.0)
        << "A zero-length temporary vector must return exactly 0.0";
}

// -----------------------------------------------------------------------------
// Test 2: no rotation at all -> both vectors identical, angle must be 0.
//         The cross product is the zero vector so no sign flip can occur.
// -----------------------------------------------------------------------------
void test_coad_identical_vectors_return_zero()
{
    std::cerr << "\n[TEST] test_coad_identical_vectors_return_zero\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    original and temporary interface-to-COM vectors\n"
              << "               are identical, (1,0,0), but the two COMs differ\n"
              << "               (a pure translation, no rotation).\n"
              << "  Criteria:    the returned angle is 0 within " << coadTol << " rad.\n";

    Molecule mol = coad_make_molecule(
        /*com*/ Coord { 0.0, 0.0, 0.0 },
        /*ifaceCoords*/ { Coord { 1.0, 0.0, 0.0 } },
        /*tmpCom*/ Coord { -20.0, 3.0, 8.0 },
        /*tmpIfaceCoords*/ { Coord { -19.0, 3.0, 8.0 } }); // tmpCom + (1,0,0)
    Complex com = coad_make_complex(Coord { 0.0, 0.0, 0.0 });

    const double theta = calc_one_angular_displacement(0, mol, com);
    std::cerr << "  Returned angle = " << theta << " rad (expected 0)\n";

    EXPECT_NEAR(theta, 0.0, coadTol)
        << "Pure translation must not register any angular displacement";
}

// -----------------------------------------------------------------------------
// Test 3: +90 degrees about +z.  v2 = (1,0,0) -> v1 = (0,1,0).
//         v2 x v1 = (0,0,1), z > 0, so the sign is flipped to -pi/2.
// -----------------------------------------------------------------------------
void test_coad_ninety_degrees_about_z_is_flipped()
{
    std::cerr << "\n[TEST] test_coad_ninety_degrees_about_z_is_flipped\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    original vector (1,0,0) rotates to (0,1,0),\n"
              << "               i.e. +90 deg about the +z axis.\n"
              << "  Criteria:    cross-product z component is +1 (> 0), so the\n"
              << "               implementation flips the sign: expect -pi/2.\n";

    Molecule mol = coad_make_molecule(
        /*com*/ Coord { 0.0, 0.0, 0.0 },
        /*ifaceCoords*/ { Coord { 1.0, 0.0, 0.0 } },
        /*tmpCom*/ Coord { 0.0, 0.0, 0.0 },
        /*tmpIfaceCoords*/ { Coord { 0.0, 1.0, 0.0 } });
    Complex com = coad_make_complex(Coord { 0.0, 0.0, 0.0 });

    const double theta = calc_one_angular_displacement(0, mol, com);
    std::cerr << "  Returned angle = " << theta << " rad (expected " << -M_PI / 2.0 << ")\n";

    EXPECT_NEAR(theta, -M_PI / 2.0, coadTol)
        << "Rotation producing a +z cross product must be reported as negative";
}

// -----------------------------------------------------------------------------
// Test 4: -90 degrees about z.  v2 = (0,1,0) -> v1 = (1,0,0).
//         v2 x v1 = (0,0,-1), z < 0, so no flip: result stays +pi/2.
// -----------------------------------------------------------------------------
void test_coad_minus_ninety_degrees_about_z_is_positive()
{
    std::cerr << "\n[TEST] test_coad_minus_ninety_degrees_about_z_is_positive\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    original vector (0,1,0) rotates to (1,0,0),\n"
              << "               i.e. -90 deg about the +z axis.\n"
              << "  Criteria:    cross-product z component is -1 (not > 0), so no\n"
              << "               sign flip occurs: expect +pi/2.\n";

    Molecule mol = coad_make_molecule(
        /*com*/ Coord { 0.0, 0.0, 0.0 },
        /*ifaceCoords*/ { Coord { 0.0, 1.0, 0.0 } },
        /*tmpCom*/ Coord { 0.0, 0.0, 0.0 },
        /*tmpIfaceCoords*/ { Coord { 1.0, 0.0, 0.0 } });
    Complex com = coad_make_complex(Coord { 0.0, 0.0, 0.0 });

    const double theta = calc_one_angular_displacement(0, mol, com);
    std::cerr << "  Returned angle = " << theta << " rad (expected " << M_PI / 2.0 << ")\n";

    EXPECT_NEAR(theta, M_PI / 2.0, coadTol)
        << "Rotation producing a -z cross product must keep the positive sign";
}

// -----------------------------------------------------------------------------
// Test 5: rotation about the x axis.  v2 = (0,1,0) -> v1 = (0,0,1).
//         v2 x v1 = (1,0,0), whose z component is 0, so the sign is never
//         flipped regardless of the rotation direction.
// -----------------------------------------------------------------------------
void test_coad_rotation_about_x_never_flips_sign()
{
    std::cerr << "\n[TEST] test_coad_rotation_about_x_never_flips_sign\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    90 deg rotation about the x axis, in both\n"
              << "               directions: (0,1,0)->(0,0,1) and (0,0,1)->(0,1,0).\n"
              << "  Criteria:    the cross product lies along +/-x so its z component\n"
              << "               is 0; both directions must return +pi/2.\n";

    Complex com = coad_make_complex(Coord { 0.0, 0.0, 0.0 });

    // Direction 1: y -> z, cross = (+1, 0, 0)
    Molecule molA = coad_make_molecule(
        Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 1.0, 0.0 } },
        Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 0.0, 1.0 } });
    const double thetaA = calc_one_angular_displacement(0, molA, com);
    std::cerr << "  (0,1,0)->(0,0,1) angle = " << thetaA << " rad\n";
    EXPECT_NEAR(thetaA, M_PI / 2.0, coadTol)
        << "Rotation about +x must be reported as +pi/2 (z component of cross is 0)";

    // Direction 2: z -> y, cross = (-1, 0, 0)
    Molecule molB = coad_make_molecule(
        Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 0.0, 1.0 } },
        Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 1.0, 0.0 } });
    const double thetaB = calc_one_angular_displacement(0, molB, com);
    std::cerr << "  (0,0,1)->(0,1,0) angle = " << thetaB << " rad\n";
    EXPECT_NEAR(thetaB, M_PI / 2.0, coadTol)
        << "Rotation about -x must also be reported as +pi/2 (no sign information)";
}

// -----------------------------------------------------------------------------
// Test 6: exactly anti-parallel vectors -> the angle is pi and, because
//         (pi - |theta|) is not greater than 1E-12, the sign flip is suppressed
//         even though the cross product is degenerate.
// -----------------------------------------------------------------------------
void test_coad_antiparallel_returns_pi()
{
    std::cerr << "\n[TEST] test_coad_antiparallel_returns_pi\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    original vector (1,0,0) flips to (-1,0,0).\n"
              << "  Criteria:    dot_theta clamps the cosine to -1 so acos gives pi;\n"
              << "               the pi-guard prevents any sign flip. Expect +pi.\n";

    Molecule mol = coad_make_molecule(
        Coord { 0.0, 0.0, 0.0 }, { Coord { 2.0, 0.0, 0.0 } },
        Coord { 0.0, 0.0, 0.0 }, { Coord { -3.0, 0.0, 0.0 } });
    Complex com = coad_make_complex(Coord { 0.0, 0.0, 0.0 });

    const double theta = calc_one_angular_displacement(0, mol, com);
    std::cerr << "  Returned angle = " << theta << " rad (expected " << M_PI << ")\n";

    EXPECT_NEAR(theta, M_PI, coadTol)
        << "Anti-parallel vectors must give an angular displacement of pi";
}

// -----------------------------------------------------------------------------
// Test 7: the result depends only on the *directions* of the two vectors, not
//         on their lengths or on where the two centres of mass happen to be.
//         45 degrees about +z -> flipped to -pi/4.
// -----------------------------------------------------------------------------
void test_coad_is_translation_and_scale_invariant()
{
    std::cerr << "\n[TEST] test_coad_is_translation_and_scale_invariant\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    45 deg rotation about +z, with the two coordinate\n"
              << "               sets translated far apart and the temporary vector\n"
              << "               scaled up by a factor of 7.\n"
              << "  Criteria:    only the relative direction matters, so the answer\n"
              << "               is -pi/4 (flipped, because cross.z > 0).\n";

    // Original vector (1,0,0) anchored at an arbitrary COM.
    const Coord com { 10.0, -3.0, 7.0 };
    const Coord iface { com.x + 1.0, com.y + 0.0, com.z + 0.0 };

    // Temporary vector 7*(1,1,0) anchored at a completely different COM.
    const Coord tmpCom { -5.0, 2.0, 0.0 };
    const Coord tmpIface { tmpCom.x + 7.0, tmpCom.y + 7.0, tmpCom.z + 0.0 };

    Molecule mol = coad_make_molecule(com, { iface }, tmpCom, { tmpIface });
    Complex cmplx = coad_make_complex(com);

    const double theta = calc_one_angular_displacement(0, mol, cmplx);
    std::cerr << "  Returned angle = " << theta << " rad (expected " << -M_PI / 4.0 << ")\n";

    EXPECT_NEAR(theta, -M_PI / 4.0, coadTol)
        << "Angle must be independent of translation and of vector magnitudes";
}

// -----------------------------------------------------------------------------
// Test 8: the ifaceIndex1 argument really selects the requested interface.
//         Three interfaces are supplied; only interface 1 is rotated.
// -----------------------------------------------------------------------------
void test_coad_uses_the_requested_interface_index()
{
    std::cerr << "\n[TEST] test_coad_uses_the_requested_interface_index\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    a molecule with three interfaces where only\n"
              << "               interface index 1 has been rotated (by 90 deg\n"
              << "               about +z); the other two are unmoved.\n"
              << "  Criteria:    index 0 and 2 return 0, index 1 returns -pi/2.\n";

    const Coord com { 0.0, 0.0, 0.0 };
    std::vector<Coord> ifaces {
        Coord { 1.0, 0.0, 0.0 }, // interface 0: unchanged
        Coord { 1.0, 0.0, 0.0 }, // interface 1: will be rotated
        Coord { 0.0, 0.0, 2.0 } // interface 2: unchanged
    };
    std::vector<Coord> tmpIfaces {
        Coord { 1.0, 0.0, 0.0 }, // identical to original
        Coord { 0.0, 1.0, 0.0 }, // rotated +90 deg about z
        Coord { 0.0, 0.0, 2.0 } // identical to original
    };

    Molecule mol = coad_make_molecule(com, ifaces, com, tmpIfaces);
    Complex cmplx = coad_make_complex(com);

    const double theta0 = calc_one_angular_displacement(0, mol, cmplx);
    const double theta1 = calc_one_angular_displacement(1, mol, cmplx);
    const double theta2 = calc_one_angular_displacement(2, mol, cmplx);

    std::cerr << "  angle(iface 0) = " << theta0 << " rad (expected 0)\n";
    std::cerr << "  angle(iface 1) = " << theta1 << " rad (expected " << -M_PI / 2.0 << ")\n";
    std::cerr << "  angle(iface 2) = " << theta2 << " rad (expected 0)\n";

    EXPECT_NEAR(theta0, 0.0, coadTol) << "Interface 0 was not rotated";
    EXPECT_NEAR(theta1, -M_PI / 2.0, coadTol) << "Interface 1 was rotated by 90 deg about +z";
    EXPECT_NEAR(theta2, 0.0, coadTol) << "Interface 2 was not rotated";
}

// -----------------------------------------------------------------------------
// Test 9: the Complex argument is read-only as far as this routine is concerned.
//         We snapshot a few fields and verify none of them change.
// -----------------------------------------------------------------------------
void test_coad_does_not_modify_inputs()
{
    std::cerr << "\n[TEST] test_coad_does_not_modify_inputs\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    call the routine and then compare the molecule's\n"
              << "               coordinates and the complex's state to snapshots.\n"
              << "  Criteria:    the function is a pure measurement - nothing in the\n"
              << "               Molecule or Complex may be altered.\n";

    Molecule mol = coad_make_molecule(
        Coord { 1.0, 2.0, 3.0 }, { Coord { 2.0, 2.0, 3.0 } },
        Coord { -1.0, -1.0, -1.0 }, { Coord { -1.0, 0.0, -1.0 } });
    Complex cmplx = coad_make_complex(Coord { 1.0, 2.0, 3.0 });

    // Snapshots taken before the call.
    const Coord comBefore = mol.comCoord;
    const Coord ifaceBefore = mol.interfaceList[0].coord;
    const Coord tmpComBefore = mol.tmpComCoord;
    const Coord tmpIfaceBefore = mol.tmpICoords[0];
    const Coord complexComBefore = cmplx.comCoord;
    const double complexRadiusBefore = cmplx.radius;
    const size_t complexMembersBefore = cmplx.memberList.size();

    const double theta = calc_one_angular_displacement(0, mol, cmplx);
    std::cerr << "  Returned angle = " << theta << " rad\n";

    EXPECT_TRUE(mol.comCoord == comBefore) << "Molecule COM must not change";
    EXPECT_TRUE(mol.interfaceList[0].coord == ifaceBefore)
        << "Molecule interface coordinate must not change";
    EXPECT_TRUE(mol.tmpComCoord == tmpComBefore) << "Temporary COM must not change";
    EXPECT_TRUE(mol.tmpICoords[0] == tmpIfaceBefore)
        << "Temporary interface coordinate must not change";
    EXPECT_TRUE(cmplx.comCoord == complexComBefore) << "Complex COM must not change";
    EXPECT_DOUBLE_EQ(cmplx.radius, complexRadiusBefore) << "Complex radius must not change";
    EXPECT_EQ(cmplx.memberList.size(), complexMembersBefore)
        << "Complex member list must not change";
}

// -----------------------------------------------------------------------------
// Test 10: degenerate original vector.  |v1| is fine but |v2| is zero, so
//          Vector::dot_theta() takes its "magnitude too small" branch, prints a
//          warning and returns 0.0; the cross product is also zero so no flip.
// -----------------------------------------------------------------------------
void test_coad_zero_original_vector_returns_zero()
{
    std::cerr << "\n[TEST] test_coad_zero_original_vector_returns_zero\n"
              << "  Source file: calc_one_angular_displacement.cpp\n"
              << "  Scenario:    the ORIGINAL interface sits exactly on the original\n"
              << "               COM (|v2| == 0) while the temporary vector is finite.\n"
              << "  Criteria:    Vector::dot_theta detects the zero-length vector,\n"
              << "               emits a warning and returns 0.0, so the routine\n"
              << "               reports 0.0 as well.\n"
              << "  (A 'WARNING: ... magnitude 0' line from Vector::dot_theta on\n"
              << "   stdout is expected here and is not a failure.)\n";

    Molecule mol = coad_make_molecule(
        /*com*/ Coord { 4.0, 4.0, 4.0 },
        /*ifaceCoords*/ { Coord { 4.0, 4.0, 4.0 } }, // coincides with the COM
        /*tmpCom*/ Coord { 0.0, 0.0, 0.0 },
        /*tmpIfaceCoords*/ { Coord { 1.0, 0.0, 0.0 } }); // finite temporary vector
    Complex cmplx = coad_make_complex(Coord { 4.0, 4.0, 4.0 });

    const double theta = calc_one_angular_displacement(0, mol, cmplx);
    std::cerr << "  Returned angle = " << theta << " rad (expected 0)\n";

    EXPECT_NEAR(theta, 0.0, coadTol)
        << "A zero-length original vector must yield a zero angular displacement";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each helper runs inside its own TEST so that a failure
// in one scenario does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(CalcOneAngularDisplacement, PointParticleReturnsZero) { test_coad_point_particle_returns_zero(); }
TEST(CalcOneAngularDisplacement, IdenticalVectorsReturnZero) { test_coad_identical_vectors_return_zero(); }
TEST(CalcOneAngularDisplacement, NinetyDegreesAboutZIsFlipped) { test_coad_ninety_degrees_about_z_is_flipped(); }
TEST(CalcOneAngularDisplacement, MinusNinetyDegreesAboutZIsPositive) { test_coad_minus_ninety_degrees_about_z_is_positive(); }
TEST(CalcOneAngularDisplacement, RotationAboutXNeverFlipsSign) { test_coad_rotation_about_x_never_flips_sign(); }
TEST(CalcOneAngularDisplacement, AntiparallelReturnsPi) { test_coad_antiparallel_returns_pi(); }
TEST(CalcOneAngularDisplacement, TranslationAndScaleInvariant) { test_coad_is_translation_and_scale_invariant(); }
TEST(CalcOneAngularDisplacement, UsesRequestedInterfaceIndex) { test_coad_uses_the_requested_interface_index(); }
TEST(CalcOneAngularDisplacement, DoesNotModifyInputs) { test_coad_does_not_modify_inputs(); }
TEST(CalcOneAngularDisplacement, ZeroOriginalVectorReturnsZero) { test_coad_zero_original_vector_returns_zero(); }