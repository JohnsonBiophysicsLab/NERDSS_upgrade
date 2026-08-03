/*! \file test_areSameAngle.cpp
 *
 * ### Unit test for src/reactions/areSameAngle.cpp
 *
 * This file exercises the single free function defined in that source file:
 *
 *     bool areSameAngle(double ang1, double ang2)
 *
 * The implementation is a tolerance-based floating point comparison:
 *
 *     return std::abs(ang1 - ang2) < 1E-4;
 *
 * So the pass criteria for every assertion below is simply:
 *   - two angles whose absolute difference is strictly less than 1E-4 must
 *     compare as "the same angle" (function returns true), and
 *   - two angles whose absolute difference is >= 1E-4 must compare as
 *     "different angles" (function returns false).
 *
 * Note that the function does NOT compare magnitudes (the commented-out
 * legacy behaviour), so an angle and its negation are only "the same" when
 * the angle is (nearly) zero. There is also no 2*pi wrap-around handling,
 * so pi and -pi are reported as different. Both of those properties are
 * pinned down explicitly by the tests below so a future behaviour change is
 * caught.
 *
 * All output is written to stderr so the reader can follow exactly which
 * function in which source file is being exercised and with what inputs.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

// -----------------------------------------------------------------------------
// Small local helper: run one comparison and report it verbosely.
//
// Using EXPECT_EQ (non-fatal) keeps every case running even if one fails.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Call areSameAngle() and check the result against an expectation.
 *
 * \param[in] ang1     first angle passed to areSameAngle()
 * \param[in] ang2     second angle passed to areSameAngle()
 * \param[in] expected the value areSameAngle() is required to return
 * \param[in] why      human readable description of the pass criterion
 */
void asa_check(double ang1, double ang2, bool expected, const char* why)
{
    const bool actual = areSameAngle(ang1, ang2);
    std::cerr << "    areSameAngle(" << ang1 << ", " << ang2 << ") -> "
              << std::boolalpha << actual
              << "  | expected " << expected
              << "  | |diff| = " << std::abs(ang1 - ang2)
              << "  | " << why << '\n';

    EXPECT_EQ(actual, expected)
        << "areSameAngle(" << ang1 << ", " << ang2 << ") returned " << actual
        << " but " << expected << " was required (" << why << ")";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: exactly identical inputs must always be reported as the same angle.
// -----------------------------------------------------------------------------
void asa_test_identical_angles()
{
    std::cerr << "\n[TEST] asa_test_identical_angles\n"
              << "  Source file:   src/reactions/areSameAngle.cpp\n"
              << "  Function:      areSameAngle(double, double)\n"
              << "  Scenario:      both arguments are bit-for-bit identical.\n"
              << "  Pass criteria: |ang1 - ang2| == 0 < 1E-4, so true.\n";

    asa_check(0.0, 0.0, true, "zero compared with zero");
    asa_check(1.0, 1.0, true, "positive value compared with itself");
    asa_check(-1.0, -1.0, true, "negative value compared with itself");
    asa_check(M_PI, M_PI, true, "pi compared with itself");
    asa_check(-M_PI_2, -M_PI_2, true, "-pi/2 compared with itself");
    asa_check(1234.5678, 1234.5678, true, "large value compared with itself");
}

// -----------------------------------------------------------------------------
// Test 2: differences strictly smaller than the 1E-4 tolerance are "same".
// -----------------------------------------------------------------------------
void asa_test_within_tolerance()
{
    std::cerr << "\n[TEST] asa_test_within_tolerance\n"
              << "  Source file:   src/reactions/areSameAngle.cpp\n"
              << "  Function:      areSameAngle(double, double)\n"
              << "  Scenario:      angles differ by less than the 1E-4 tolerance.\n"
              << "  Pass criteria: function returns true for every pair.\n";

    // Difference of 1E-5, an order of magnitude below the tolerance.
    asa_check(1.0, 1.0 + 1e-5, true, "difference 1e-5 is below tolerance");
    asa_check(1.0 + 1e-5, 1.0, true, "same pair with arguments swapped (symmetry)");

    // Difference just barely below the tolerance.
    asa_check(0.0, 9.99e-5, true, "difference 9.99e-5 is just below tolerance");
    asa_check(0.0, -9.99e-5, true, "negative difference of same magnitude");

    // Tolerance is absolute, so it also applies around large angles.
    asa_check(M_PI, M_PI - 1e-6, true, "near-pi pair differing by 1e-6");
    asa_check(-2.0, -2.0 - 5e-5, true, "negative pair differing by 5e-5");
}

// -----------------------------------------------------------------------------
// Test 3: differences at or above the 1E-4 tolerance are "different".
// -----------------------------------------------------------------------------
void asa_test_outside_tolerance()
{
    std::cerr << "\n[TEST] asa_test_outside_tolerance\n"
              << "  Source file:   src/reactions/areSameAngle.cpp\n"
              << "  Function:      areSameAngle(double, double)\n"
              << "  Scenario:      angles differ by 1E-4 or more.\n"
              << "  Pass criteria: function returns false for every pair\n"
              << "                 (the comparison is strict '<', so exactly\n"
              << "                  1E-4 must NOT be considered the same).\n";

    // Comfortably outside the tolerance.
    asa_check(0.0, 1e-3, false, "difference 1e-3 exceeds tolerance");
    asa_check(1.0, 1.5, false, "difference 0.5 exceeds tolerance");
    asa_check(-1.0, 1.0, false, "difference 2.0 exceeds tolerance");
    asa_check(M_PI, M_PI_2, false, "pi vs pi/2 are clearly different");

    // Just above the tolerance: 1.01e-4.
    asa_check(0.0, 1.01e-4, false, "difference 1.01e-4 is just above tolerance");
    asa_check(0.0, -1.01e-4, false, "negative difference just above tolerance");
}

// -----------------------------------------------------------------------------
// Test 4: the boundary itself. The implementation uses a strict '<', so a
//         difference of exactly 1E-4 must be reported as different.
// -----------------------------------------------------------------------------
void asa_test_tolerance_boundary()
{
    std::cerr << "\n[TEST] asa_test_tolerance_boundary\n"
              << "  Source file:   src/reactions/areSameAngle.cpp\n"
              << "  Function:      areSameAngle(double, double)\n"
              << "  Scenario:      difference is exactly the tolerance 1E-4.\n"
              << "  Pass criteria: strict '<' comparison means false.\n";

    // 0.0 and 1e-4 are both exactly representable operations here, so the
    // subtraction yields exactly 1e-4 and the strict comparison fails.
    asa_check(0.0, 1e-4, false, "exactly 1e-4 apart -> not within strict tolerance");
    asa_check(1e-4, 0.0, false, "same boundary case with arguments swapped");
    asa_check(0.0, -1e-4, false, "exactly -1e-4 apart -> also outside");
}

// -----------------------------------------------------------------------------
// Test 5: the function is symmetric in its arguments and reflexive.
//         These are properties association code relies on when checking
//         association angles from either reactant's point of view.
// -----------------------------------------------------------------------------
void asa_test_symmetry_and_reflexivity()
{
    std::cerr << "\n[TEST] asa_test_symmetry_and_reflexivity\n"
              << "  Source file:   src/reactions/areSameAngle.cpp\n"
              << "  Function:      areSameAngle(double, double)\n"
              << "  Scenario:      sweep a set of angle pairs and confirm\n"
              << "                 areSameAngle(a,b) == areSameAngle(b,a) and\n"
              << "                 areSameAngle(a,a) == true.\n"
              << "  Pass criteria: symmetry holds for all pairs and every\n"
              << "                 angle equals itself.\n";

    const double angles[] = { -M_PI, -1.5, -1e-5, 0.0, 1e-5, 0.5, M_PI_2, M_PI, 10.0 };
    const int numAngles = static_cast<int>(sizeof(angles) / sizeof(angles[0]));

    for (int i = 0; i < numAngles; ++i) {
        // Reflexivity: an angle is always the same as itself.
        const bool self = areSameAngle(angles[i], angles[i]);
        std::cerr << "    reflexive: areSameAngle(" << angles[i] << ", "
                  << angles[i] << ") -> " << std::boolalpha << self << '\n';
        EXPECT_TRUE(self) << "angle " << angles[i] << " must equal itself";

        // Symmetry: order of arguments must not matter.
        for (int j = 0; j < numAngles; ++j) {
            const bool forward = areSameAngle(angles[i], angles[j]);
            const bool reverse = areSameAngle(angles[j], angles[i]);
            EXPECT_EQ(forward, reverse)
                << "symmetry violated for (" << angles[i] << ", " << angles[j]
                << "): forward=" << forward << " reverse=" << reverse;
        }
    }
    std::cerr << "    symmetry checked for all " << numAngles << "x" << numAngles
              << " angle combinations.\n";
}

// -----------------------------------------------------------------------------
// Test 6: document the *absence* of magnitude comparison and of 2*pi
//         wrap-around handling. The commented-out legacy implementation
//         compared |ang1| == |ang2|; the current one does not. These
//         assertions lock in the current, documented behaviour.
// -----------------------------------------------------------------------------
void asa_test_sign_and_wraparound_behavior()
{
    std::cerr << "\n[TEST] asa_test_sign_and_wraparound_behavior\n"
              << "  Source file:   src/reactions/areSameAngle.cpp\n"
              << "  Function:      areSameAngle(double, double)\n"
              << "  Scenario:      opposite-sign angles and 2*pi-equivalent\n"
              << "                 angles.\n"
              << "  Pass criteria: the function compares raw values, so\n"
              << "                 +x vs -x is false (unless x ~ 0) and\n"
              << "                 pi vs -pi / 0 vs 2*pi are also false.\n";

    // Opposite signs: NOT the same angle under the current implementation.
    asa_check(M_PI_2, -M_PI_2, false, "pi/2 vs -pi/2 (no magnitude compare)");
    asa_check(1.0, -1.0, false, "1 vs -1 (no magnitude compare)");
    asa_check(M_PI, -M_PI, false, "pi vs -pi (no 2*pi wrap-around)");

    // A tiny angle and its negation ARE the same, because their difference
    // (2e-5) is still inside the tolerance. This is a tolerance effect, not
    // a magnitude comparison.
    asa_check(1e-5, -1e-5, true, "+/-1e-5 differ by 2e-5, inside tolerance");

    // No modular arithmetic: 0 and 2*pi are treated as different.
    asa_check(0.0, 2.0 * M_PI, false, "0 vs 2*pi (no modular reduction)");
}

// -----------------------------------------------------------------------------
// Test 7: behaviour with non-finite inputs. NaN comparisons are always false
//         (any arithmetic with NaN yields NaN, and NaN < 1E-4 is false), and
//         infinities differ from finite values.
// -----------------------------------------------------------------------------
void asa_test_non_finite_inputs()
{
    std::cerr << "\n[TEST] asa_test_non_finite_inputs\n"
              << "  Source file:   src/reactions/areSameAngle.cpp\n"
              << "  Function:      areSameAngle(double, double)\n"
              << "  Scenario:      NaN and infinite inputs (the association\n"
              << "                 angle defaults are quiet NaN).\n"
              << "  Pass criteria: any comparison involving NaN is false;\n"
              << "                 infinity vs a finite angle is false.\n";

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    // NaN never compares "same", not even with itself, because
    // std::abs(NaN - NaN) = NaN and (NaN < 1E-4) is false.
    asa_check(nan, nan, false, "NaN vs NaN -> comparison with NaN is false");
    asa_check(nan, 0.0, false, "NaN vs 0 -> false");
    asa_check(0.0, nan, false, "0 vs NaN -> false (symmetric)");

    // Infinity differs from any finite angle by infinity.
    asa_check(inf, 0.0, false, "+inf vs 0 -> false");
    asa_check(-inf, 0.0, false, "-inf vs 0 -> false");

    // inf - inf is NaN, so even identical infinities are reported different.
    asa_check(inf, inf, false, "+inf vs +inf -> inf-inf is NaN -> false");
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper is invoked from its own TEST
// so the framework reports individual results while every case still runs
// (all assertions above are non-fatal EXPECT_* style).
// -----------------------------------------------------------------------------
TEST(AreSameAngleTest, IdenticalAngles) { asa_test_identical_angles(); }
TEST(AreSameAngleTest, WithinTolerance) { asa_test_within_tolerance(); }
TEST(AreSameAngleTest, OutsideTolerance) { asa_test_outside_tolerance(); }
TEST(AreSameAngleTest, ToleranceBoundary) { asa_test_tolerance_boundary(); }
TEST(AreSameAngleTest, SymmetryAndReflexivity) { asa_test_symmetry_and_reflexivity(); }
TEST(AreSameAngleTest, SignAndWraparoundBehavior) { asa_test_sign_and_wraparound_behavior(); }
TEST(AreSameAngleTest, NonFiniteInputs) { asa_test_non_finite_inputs(); }