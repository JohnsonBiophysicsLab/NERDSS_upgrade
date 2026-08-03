/*! \file test_areParallel.cpp
 *
 * ### Unit test for ../src/reactions/areParallel.cpp
 *
 * The file under test contains exactly one free function:
 *
 *     bool areParallel(const double& angle)
 *     {
 *         return angle == M_PI || angle == 0;
 *     }
 *
 * Semantics: the function reports whether two vectors, whose separation angle
 * (in radians) is given by `angle`, are parallel or anti-parallel.  It does so
 * with an *exact* floating point comparison against M_PI and 0 — there is no
 * tolerance.  The tests below therefore verify:
 *
 *   1. The two "true" inputs: exactly 0 and exactly M_PI.
 *   2. Common "false" inputs: pi/2, pi/4, -M_PI, 2*M_PI, etc.
 *   3. The strictness of the comparison: values a few ULPs away from 0 or M_PI
 *      must report false (this documents the behaviour, tolerance-free).
 *   4. Signed-zero behaviour: -0.0 == 0.0 in IEEE-754, so -0.0 is parallel.
 *   5. Special values: NaN and +/-infinity must report false.
 *   6. That values produced by typical trig round-trips (acos(1), acos(-1))
 *      still compare exactly equal, since these are the realistic call sites.
 *
 * Every check uses EXPECT_* (non-fatal) so all assertions execute even when
 * an earlier one fails, and every step prints to stderr what is being tested.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <limits>

// -----------------------------------------------------------------------------
// Small helper that runs areParallel() and echoes the input/result pair so the
// console log is self-documenting.  Kept in an anonymous namespace so it cannot
// collide with symbols from other translation units in the test suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Call areParallel() and print the input and the returned value.
 *
 * \param[in] label  Human readable description of the input.
 * \param[in] angle  The angle (radians) handed to areParallel().
 * \return The value returned by areParallel().
 */
bool ap_probe(const char* label, double angle)
{
    const bool result = areParallel(angle);
    std::cerr << "    areParallel(" << label << " = " << angle
              << ") -> " << (result ? "true" : "false") << '\n';
    return result;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the two inputs that must return true.
//
// Criteria: exactly 0.0 (vectors parallel) and exactly M_PI (vectors
// anti-parallel) both report "parallel".
// -----------------------------------------------------------------------------
void test_areparallel_true_cases()
{
    std::cerr << "\n[TEST] test_areparallel_true_cases\n"
              << "  Source file: src/reactions/areParallel.cpp\n"
              << "  Function:    areParallel(const double&)\n"
              << "  Scenario:    the only two angles the function accepts.\n"
              << "  Pass:        areParallel(0.0) and areParallel(M_PI) are true.\n";

    // Angle of 0 radians: the vectors point in the same direction.
    EXPECT_TRUE(ap_probe("0.0", 0.0))
        << "An angle of exactly 0 radians must be reported as parallel";

    // Angle of pi radians: the vectors point in opposite directions, which the
    // association code also treats as "parallel" (co-linear).
    EXPECT_TRUE(ap_probe("M_PI", M_PI))
        << "An angle of exactly M_PI radians must be reported as parallel";
}

// -----------------------------------------------------------------------------
// Test 2: representative angles that must return false.
//
// Criteria: any angle that is neither exactly 0 nor exactly M_PI is not
// parallel.  We include a negative pi, since the implementation does NOT take
// the absolute value, and multiples of pi, which are also not special-cased.
// -----------------------------------------------------------------------------
void test_areparallel_false_cases()
{
    std::cerr << "\n[TEST] test_areparallel_false_cases\n"
              << "  Source file: src/reactions/areParallel.cpp\n"
              << "  Function:    areParallel(const double&)\n"
              << "  Scenario:    ordinary non-parallel angles, plus -M_PI and\n"
              << "               2*M_PI which the implementation does not fold.\n"
              << "  Pass:        every one of these returns false.\n";

    // Perpendicular vectors.
    EXPECT_FALSE(ap_probe("M_PI/2", M_PI / 2.0))
        << "pi/2 (perpendicular) must not be reported as parallel";

    // 45 degrees.
    EXPECT_FALSE(ap_probe("M_PI/4", M_PI / 4.0))
        << "pi/4 must not be reported as parallel";

    // 135 degrees.
    EXPECT_FALSE(ap_probe("3*M_PI/4", 3.0 * M_PI / 4.0))
        << "3*pi/4 must not be reported as parallel";

    // Negative pi: mathematically anti-parallel, but the comparison is exact
    // and unsigned-agnostic only for zero, so this must be false.
    EXPECT_FALSE(ap_probe("-M_PI", -M_PI))
        << "-M_PI is not compared against, so it must return false";

    // A full turn: again not special-cased.
    EXPECT_FALSE(ap_probe("2*M_PI", 2.0 * M_PI))
        << "2*M_PI is not compared against, so it must return false";

    // An arbitrary small non-zero angle.
    EXPECT_FALSE(ap_probe("0.1", 0.1))
        << "0.1 rad must not be reported as parallel";

    // An arbitrary large angle.
    EXPECT_FALSE(ap_probe("10.0", 10.0))
        << "10.0 rad must not be reported as parallel";
}

// -----------------------------------------------------------------------------
// Test 3: the comparison is exact — no tolerance band around 0 or M_PI.
//
// Criteria: values one ULP (and one small epsilon) away from 0 or M_PI must
// return false.  This documents the tolerance-free behaviour of the function.
// -----------------------------------------------------------------------------
void test_areparallel_exactness()
{
    std::cerr << "\n[TEST] test_areparallel_exactness\n"
              << "  Source file: src/reactions/areParallel.cpp\n"
              << "  Function:    areParallel(const double&)\n"
              << "  Scenario:    angles a hair away from 0 and M_PI.\n"
              << "  Pass:        all return false (comparison is exact, no tolerance).\n";

    // One ULP above zero (smallest positive denormal) - not equal to 0.
    const double tinyPositive = std::numeric_limits<double>::denorm_min();
    EXPECT_FALSE(ap_probe("denorm_min()", tinyPositive))
        << "The smallest positive subnormal is not exactly 0, so not parallel";

    // Machine epsilon above zero.
    const double epsAboveZero = std::numeric_limits<double>::epsilon();
    EXPECT_FALSE(ap_probe("epsilon()", epsAboveZero))
        << "epsilon is not exactly 0, so not parallel";

    // One ULP below and one ULP above M_PI via std::nextafter.
    const double justBelowPi = std::nextafter(M_PI, 0.0);
    const double justAbovePi = std::nextafter(M_PI, 4.0);
    EXPECT_FALSE(ap_probe("nextafter(M_PI, 0)", justBelowPi))
        << "One ULP below M_PI is not exactly M_PI, so not parallel";
    EXPECT_FALSE(ap_probe("nextafter(M_PI, 4)", justAbovePi))
        << "One ULP above M_PI is not exactly M_PI, so not parallel";

    // Sanity check that the perturbed values really differ from M_PI, otherwise
    // the assertions above would be vacuous on an exotic platform.
    EXPECT_NE(justBelowPi, M_PI) << "nextafter should have produced a different value";
    EXPECT_NE(justAbovePi, M_PI) << "nextafter should have produced a different value";
}

// -----------------------------------------------------------------------------
// Test 4: signed zero.
//
// Criteria: IEEE-754 states -0.0 == 0.0, therefore areParallel(-0.0) must be
// true even though the bit pattern differs from +0.0.
// -----------------------------------------------------------------------------
void test_areparallel_signed_zero()
{
    std::cerr << "\n[TEST] test_areparallel_signed_zero\n"
              << "  Source file: src/reactions/areParallel.cpp\n"
              << "  Function:    areParallel(const double&)\n"
              << "  Scenario:    negative zero input.\n"
              << "  Pass:        -0.0 compares equal to 0 and is reported parallel.\n";

    const double negZero = -0.0;

    // Confirm we really have a negative zero (bit-wise different from +0.0).
    EXPECT_TRUE(std::signbit(negZero))
        << "The test input should carry a negative sign bit";

    // ... and yet it must compare equal to 0 and therefore be "parallel".
    EXPECT_TRUE(ap_probe("-0.0", negZero))
        << "-0.0 == 0.0 in IEEE-754, so areParallel(-0.0) must be true";
}

// -----------------------------------------------------------------------------
// Test 5: special IEEE values.
//
// Criteria: NaN never compares equal to anything (including itself), and the
// infinities are not equal to 0 or M_PI, so all must return false.
// -----------------------------------------------------------------------------
void test_areparallel_special_values()
{
    std::cerr << "\n[TEST] test_areparallel_special_values\n"
              << "  Source file: src/reactions/areParallel.cpp\n"
              << "  Function:    areParallel(const double&)\n"
              << "  Scenario:    NaN, +infinity and -infinity inputs.\n"
              << "  Pass:        all three return false (no equality can hold).\n";

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double posInf = std::numeric_limits<double>::infinity();
    const double negInf = -std::numeric_limits<double>::infinity();

    // NaN compares unequal to every value, so the function must return false.
    EXPECT_FALSE(ap_probe("quiet_NaN", nan))
        << "NaN cannot equal 0 or M_PI, so areParallel(NaN) must be false";

    // Infinities are finite-value-unequal, so also false.
    EXPECT_FALSE(ap_probe("+infinity", posInf))
        << "+inf is neither 0 nor M_PI, so areParallel(+inf) must be false";
    EXPECT_FALSE(ap_probe("-infinity", negInf))
        << "-inf is neither 0 nor M_PI, so areParallel(-inf) must be false";
}

// -----------------------------------------------------------------------------
// Test 6: realistic call sites.
//
// In the association code the angle handed to areParallel() usually comes from
// an acos() of a normalized dot product.  Verify that the canonical parallel /
// anti-parallel dot products (+1 and -1) round-trip through acos() to values
// that areParallel() still recognises, and that an in-between dot product does
// not.
// -----------------------------------------------------------------------------
void test_areparallel_trig_roundtrip()
{
    std::cerr << "\n[TEST] test_areparallel_trig_roundtrip\n"
              << "  Source file: src/reactions/areParallel.cpp\n"
              << "  Function:    areParallel(const double&)\n"
              << "  Scenario:    angles obtained from acos() of dot products,\n"
              << "               the way association code actually calls it.\n"
              << "  Pass:        acos(1) -> parallel, acos(-1) -> parallel,\n"
              << "               acos(0) and acos(0.5) -> not parallel.\n";

    // Dot product of +1 for unit vectors => angle 0 exactly.
    const double angleFromDotPlusOne = std::acos(1.0);
    EXPECT_DOUBLE_EQ(angleFromDotPlusOne, 0.0)
        << "acos(1.0) should yield exactly 0";
    EXPECT_TRUE(ap_probe("acos(1.0)", angleFromDotPlusOne))
        << "A dot product of +1 corresponds to parallel vectors";

    // Dot product of -1 for unit vectors => angle M_PI exactly.
    const double angleFromDotMinusOne = std::acos(-1.0);
    EXPECT_DOUBLE_EQ(angleFromDotMinusOne, M_PI)
        << "acos(-1.0) should yield exactly M_PI";
    EXPECT_TRUE(ap_probe("acos(-1.0)", angleFromDotMinusOne))
        << "A dot product of -1 corresponds to anti-parallel (still 'parallel')";

    // Orthogonal vectors: dot product 0 => pi/2, definitely not parallel.
    EXPECT_FALSE(ap_probe("acos(0.0)", std::acos(0.0)))
        << "A dot product of 0 (orthogonal) must not be parallel";

    // 60 degrees: dot product 0.5, not parallel.
    EXPECT_FALSE(ap_probe("acos(0.5)", std::acos(0.5)))
        << "A dot product of 0.5 (60 degrees) must not be parallel";
}

// -----------------------------------------------------------------------------
// Test 7: determinism / purity.
//
// areParallel() has no state; verify repeated calls with the same input always
// yield the same answer, and that calling it with other values in between does
// not perturb the result.
// -----------------------------------------------------------------------------
void test_areparallel_determinism()
{
    std::cerr << "\n[TEST] test_areparallel_determinism\n"
              << "  Source file: src/reactions/areParallel.cpp\n"
              << "  Function:    areParallel(const double&)\n"
              << "  Scenario:    repeated interleaved calls.\n"
              << "  Pass:        results are stable across repetitions (pure function).\n";

    // Interleave "true" and "false" inputs several times; results must be stable.
    for (int i = 0; i < 3; ++i) {
        std::cerr << "  Iteration " << i << ":\n";
        EXPECT_TRUE(areParallel(0.0)) << "areParallel(0.0) must always be true";
        EXPECT_FALSE(areParallel(M_PI / 3.0)) << "areParallel(pi/3) must always be false";
        EXPECT_TRUE(areParallel(M_PI)) << "areParallel(M_PI) must always be true";
        EXPECT_FALSE(areParallel(-1.0)) << "areParallel(-1.0) must always be false";
    }
    std::cerr << "    3 interleaved iterations completed with stable results\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: one TEST per named test_* function so the framework
// reports individual results while still executing every scenario.
// -----------------------------------------------------------------------------
TEST(AreParallelTest, TrueCases) { test_areparallel_true_cases(); }
TEST(AreParallelTest, FalseCases) { test_areparallel_false_cases(); }
TEST(AreParallelTest, Exactness) { test_areparallel_exactness(); }
TEST(AreParallelTest, SignedZero) { test_areparallel_signed_zero(); }
TEST(AreParallelTest, SpecialValues) { test_areparallel_special_values(); }
TEST(AreParallelTest, TrigRoundTrip) { test_areparallel_trig_roundtrip(); }
TEST(AreParallelTest, Determinism) { test_areparallel_determinism(); }