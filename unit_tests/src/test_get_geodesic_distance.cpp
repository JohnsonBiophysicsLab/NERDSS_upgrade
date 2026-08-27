/*! \file test_get_geodesic_distance.cpp
 *
 * ### Unit test for src/reactions/get_geodesic_distance.cpp
 *
 * The single function under test is:
 *
 *     double get_geodesic_distance(Coord intFace1, Coord intFace2)
 *
 * Behaviour, read directly from the implementation:
 *   - r1 = |intFace1|, r2 = |intFace2|                (Coord::get_magnitude)
 *   - dotProduct = x1*x2 + y1*y2 + z1*z2
 *   - theta = acos( dotProduct / (r1 * r2) )          (angle between the vectors)
 *   - meanR = (r1 + r2) / 2                           (arithmetic mean radius)
 *   - returns meanR * theta                           (arc length on a sphere of
 *                                                      radius meanR)
 *
 * So the return value is the great-circle (geodesic) arc length between the two
 * points, evaluated on a sphere whose radius is the average of the two input
 * radii.  Note the arguments are taken **by value**, so the caller's Coord
 * objects can never be modified.
 *
 * Every test below prints what is being exercised and what the pass criterion
 * is, and uses non-fatal EXPECT_* assertions so a single failure does not stop
 * the rest of the suite.
 */

#include "classes/class_Coord.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <limits>

namespace {

//! Tolerance used for all floating point comparisons in this file.
//  The computation is a handful of multiplications plus one acos(), so results
//  are accurate to well within 1e-12 of the analytic answer.
constexpr double kGgdTol = 1e-12;

/*! \brief Convenience helper that prints the two inputs and the returned value.
 *
 * Keeps the console output uniform across all of the tests below.
 */
double ggd_call_and_report(const char* label, const Coord& a, const Coord& b)
{
    Coord aCopy { a };
    Coord bCopy { b };
    const double dist = get_geodesic_distance(aCopy, bCopy);
    std::cerr << "    " << label << ": p1 = (" << a.x << ", " << a.y << ", " << a.z
              << "), p2 = (" << b.x << ", " << b.y << ", " << b.z
              << ")  ->  geodesic distance = " << dist << '\n';
    return dist;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Coincident / parallel points must give a zero arc length.
// -----------------------------------------------------------------------------
void test_ggd_zero_angle_gives_zero_distance()
{
    std::cerr << "\n[TEST] test_ggd_zero_angle_gives_zero_distance\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    two points that lie along exactly the same ray\n"
              << "               from the origin, so the enclosed angle is 0.\n"
              << "  Criterion:   theta == 0 so meanR * theta == 0 regardless of\n"
              << "               how far apart the two radii are.\n";

    // Identical points: dot = r^2, ratio is exactly 1, acos(1) == 0.
    // Values are chosen so dotProduct / (r1*r2) is exactly representable and
    // acos() cannot be handed a value slightly greater than 1.
    const double d1 = ggd_call_and_report("identical points", Coord { 3.0, 0.0, 0.0 },
        Coord { 3.0, 0.0, 0.0 });
    EXPECT_NEAR(d1, 0.0, kGgdTol) << "Identical interfaces must have zero geodesic separation";

    // Co-linear, same direction, different magnitudes: dot = 5, r1*r2 = 5.
    // theta is still 0, so the mean radius (3) is multiplied by zero.
    const double d2 = ggd_call_and_report("co-linear, different radii", Coord { 1.0, 0.0, 0.0 },
        Coord { 5.0, 0.0, 0.0 });
    EXPECT_NEAR(d2, 0.0, kGgdTol)
        << "Points on the same ray have angle 0, hence zero arc length even with unequal radii";
}

// -----------------------------------------------------------------------------
// Test 2: Orthogonal points on a common sphere -> quarter of a great circle.
// -----------------------------------------------------------------------------
void test_ggd_orthogonal_points_on_common_sphere()
{
    std::cerr << "\n[TEST] test_ggd_orthogonal_points_on_common_sphere\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    two perpendicular points at the same radius R.\n"
              << "  Criterion:   theta == pi/2 and meanR == R, so the answer is\n"
              << "               exactly R * pi / 2 (a quarter great circle).\n";

    const double R = 10.0;

    // (R,0,0) and (0,R,0): dot = 0 so acos(0) = pi/2.
    const double expectedXY = R * (M_PI / 2.0);
    const double dXY = ggd_call_and_report("x-axis vs y-axis", Coord { R, 0.0, 0.0 },
        Coord { 0.0, R, 0.0 });
    EXPECT_NEAR(dXY, expectedXY, kGgdTol)
        << "Perpendicular points at radius " << R << " should be R*pi/2 = " << expectedXY << " apart";

    // Same check in a different plane, to make sure all three components of the
    // dot product are being used.
    const double dYZ = ggd_call_and_report("y-axis vs z-axis", Coord { 0.0, R, 0.0 },
        Coord { 0.0, 0.0, R });
    EXPECT_NEAR(dYZ, expectedXY, kGgdTol)
        << "The y/z pair must give the same quarter-circle arc as the x/y pair";
}

// -----------------------------------------------------------------------------
// Test 3: Antipodal points -> half of a great circle.
// -----------------------------------------------------------------------------
void test_ggd_antipodal_points()
{
    std::cerr << "\n[TEST] test_ggd_antipodal_points\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    diametrically opposed points at radius R.\n"
              << "  Criterion:   dot/(r1*r2) == -1 so theta == pi, giving R*pi,\n"
              << "               the largest value the function can return for a\n"
              << "               given radius.\n";

    const double R = 4.0;
    const double expected = R * M_PI;

    const double dist = ggd_call_and_report("antipodal along x", Coord { R, 0.0, 0.0 },
        Coord { -R, 0.0, 0.0 });
    EXPECT_NEAR(dist, expected, kGgdTol)
        << "Antipodal points at radius " << R << " should be half a great circle, R*pi = " << expected;

    // The maximum possible arc for this radius; anything larger would mean acos
    // returned something greater than pi, which cannot happen.
    EXPECT_LE(dist, expected + kGgdTol) << "Arc length can never exceed meanR * pi";
}

// -----------------------------------------------------------------------------
// Test 4: A 45 degree separation, verifying the general acos() path.
// -----------------------------------------------------------------------------
void test_ggd_forty_five_degree_separation()
{
    std::cerr << "\n[TEST] test_ggd_forty_five_degree_separation\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    p1 = (1,0,0), p2 = (1,1,0).\n"
              << "               dot = 1, r1 = 1, r2 = sqrt(2), so the cosine is\n"
              << "               1/sqrt(2) and theta = pi/4.\n"
              << "  Criterion:   result == ((1 + sqrt(2))/2) * (pi/4).\n";

    const Coord p1 { 1.0, 0.0, 0.0 };
    const Coord p2 { 1.0, 1.0, 0.0 };

    const double r1 = 1.0;
    const double r2 = std::sqrt(2.0);
    const double meanR = (r1 + r2) / 2.0;
    const double expected = meanR * (M_PI / 4.0);

    const double dist = ggd_call_and_report("45 degree pair", p1, p2);
    std::cerr << "    expected = meanR(" << meanR << ") * pi/4 = " << expected << '\n';

    EXPECT_NEAR(dist, expected, kGgdTol)
        << "General angle case must use acos(dot/(r1*r2)) and the arithmetic mean radius";
}

// -----------------------------------------------------------------------------
// Test 5: The mean radius really is the arithmetic mean of the two magnitudes.
// -----------------------------------------------------------------------------
void test_ggd_uses_arithmetic_mean_radius()
{
    std::cerr << "\n[TEST] test_ggd_uses_arithmetic_mean_radius\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    perpendicular points with *different* radii,\n"
              << "               p1 = (1,0,0) (r=1) and p2 = (0,3,0) (r=3).\n"
              << "  Criterion:   theta is still pi/2, and the radius used is the\n"
              << "               arithmetic mean (1+3)/2 = 2, so the answer is\n"
              << "               2 * pi/2 = pi.  This distinguishes the mean from,\n"
              << "               e.g., the min, max or geometric mean radius.\n";

    const double expected = 2.0 * (M_PI / 2.0); // == pi
    const double dist = ggd_call_and_report("r=1 vs r=3, perpendicular", Coord { 1.0, 0.0, 0.0 },
        Coord { 0.0, 3.0, 0.0 });

    std::cerr << "    expected = pi = " << expected << '\n';
    EXPECT_NEAR(dist, expected, kGgdTol)
        << "Mismatched radii must be combined as (r1+r2)/2 before scaling the angle";

    // Sanity checks against the alternative radius definitions, so a regression
    // that swaps the mean for min/max/geometric mean is caught explicitly.
    const double minRadiusAnswer = 1.0 * (M_PI / 2.0);
    const double maxRadiusAnswer = 3.0 * (M_PI / 2.0);
    const double geoMeanAnswer = std::sqrt(3.0) * (M_PI / 2.0);
    EXPECT_GT(std::abs(dist - minRadiusAnswer), 1e-6) << "Result must not be based on the minimum radius";
    EXPECT_GT(std::abs(dist - maxRadiusAnswer), 1e-6) << "Result must not be based on the maximum radius";
    EXPECT_GT(std::abs(dist - geoMeanAnswer), 1e-6) << "Result must not be based on the geometric mean radius";
}

// -----------------------------------------------------------------------------
// Test 6: The function is symmetric in its two arguments.
// -----------------------------------------------------------------------------
void test_ggd_symmetry_in_arguments()
{
    std::cerr << "\n[TEST] test_ggd_symmetry_in_arguments\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    call the function with the two interfaces in\n"
              << "               both orders.\n"
              << "  Criterion:   the dot product and the mean radius are both\n"
              << "               symmetric, so f(a,b) must equal f(b,a) exactly.\n";

    const Coord a { 2.0, -1.0, 0.5 };
    const Coord b { -0.75, 3.0, 1.25 };

    const double forward = ggd_call_and_report("forward order", a, b);
    const double reverse = ggd_call_and_report("reversed order", b, a);

    EXPECT_DOUBLE_EQ(forward, reverse)
        << "get_geodesic_distance must be symmetric under exchange of its arguments";
}

// -----------------------------------------------------------------------------
// Test 7: Scaling both points scales the distance linearly.
// -----------------------------------------------------------------------------
void test_ggd_scales_linearly_with_radius()
{
    std::cerr << "\n[TEST] test_ggd_scales_linearly_with_radius\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    multiply both interface coordinates by k = 7.\n"
              << "  Criterion:   the angle is scale invariant while meanR scales\n"
              << "               by k, so the returned arc length must scale by k.\n";

    const Coord a { 1.0, 2.0, -3.0 };
    const Coord b { -2.0, 0.5, 1.0 };
    const double k = 7.0;

    const double base = ggd_call_and_report("unscaled pair", a, b);
    const double scaled = ggd_call_and_report("scaled by k=7",
        Coord { k * a.x, k * a.y, k * a.z }, Coord { k * b.x, k * b.y, k * b.z });

    std::cerr << "    k * base = " << (k * base) << '\n';
    // A relative tolerance is used here because the magnitudes involved are
    // larger than in the analytic cases above.
    EXPECT_NEAR(scaled, k * base, 1e-10)
        << "Uniformly scaling both interfaces by k must scale the geodesic distance by k";
}

// -----------------------------------------------------------------------------
// Test 8: Arguments are taken by value, so the caller's Coords are untouched.
// -----------------------------------------------------------------------------
void test_ggd_does_not_modify_inputs()
{
    std::cerr << "\n[TEST] test_ggd_does_not_modify_inputs\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    the signature takes both Coords by value.\n"
              << "  Criterion:   the caller's Coord objects still hold their\n"
              << "               original x/y/z after the call.\n";

    Coord p1 { 1.5, -2.5, 3.5 };
    Coord p2 { -1.0, 4.0, 0.25 };

    const double dist = get_geodesic_distance(p1, p2);
    std::cerr << "    returned distance = " << dist << '\n';

    // Coord::operator== compares components after roundv() (4 decimal places),
    // so compare the raw doubles directly here to be strict about "unmodified".
    EXPECT_DOUBLE_EQ(p1.x, 1.5) << "p1.x must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(p1.y, -2.5) << "p1.y must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(p1.z, 3.5) << "p1.z must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(p2.x, -1.0) << "p2.x must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(p2.y, 4.0) << "p2.y must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(p2.z, 0.25) << "p2.z must be unchanged (pass by value)";
}

// -----------------------------------------------------------------------------
// Test 9: Results are always finite, non-negative, and bounded by meanR * pi
//         for a spread of arbitrary inputs.
// -----------------------------------------------------------------------------
void test_ggd_result_is_bounded_and_non_negative()
{
    std::cerr << "\n[TEST] test_ggd_result_is_bounded_and_non_negative\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    a handful of arbitrary non-zero interface pairs.\n"
              << "  Criterion:   because acos() returns a value in [0, pi], every\n"
              << "               result must be finite, >= 0, and <= meanR * pi.\n";

    const Coord pairs[][2] = {
        { Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 2.0 } },
        { Coord { 5.0, 5.0, 5.0 }, Coord { -5.0, 1.0, 0.0 } },
        { Coord { 0.1, 0.2, 0.3 }, Coord { 0.3, -0.2, 0.1 } },
        { Coord { 100.0, 0.0, 0.0 }, Coord { 99.0, 1.0, -1.0 } },
    };

    for (const auto& pair : pairs) {
        const Coord& a = pair[0];
        const Coord& b = pair[1];

        Coord aCopy { a };
        Coord bCopy { b };
        const double rA = aCopy.get_magnitude();
        const double rB = bCopy.get_magnitude();
        const double meanR = (rA + rB) / 2.0;
        const double upperBound = meanR * M_PI;

        const double dist = ggd_call_and_report("bounded check", a, b);
        std::cerr << "      meanR = " << meanR << ", upper bound = " << upperBound << '\n';

        EXPECT_TRUE(std::isfinite(dist)) << "Distance must be finite for non-zero interface vectors";
        EXPECT_GE(dist, -kGgdTol) << "Arc length can never be negative";
        EXPECT_LE(dist, upperBound + 1e-9) << "Arc length can never exceed meanR * pi";
    }
}

// -----------------------------------------------------------------------------
// Test 10: Documented degenerate input.  A zero-length interface vector makes
//          r1*r2 == 0, so the implementation performs 0.0/0.0 and hands NaN to
//          acos().  The function does not guard against this and does not exit,
//          so it is safe to exercise here; we simply document the behaviour.
// -----------------------------------------------------------------------------
void test_ggd_zero_vector_yields_nan()
{
    std::cerr << "\n[TEST] test_ggd_zero_vector_yields_nan\n"
              << "  Source file: src/reactions/get_geodesic_distance.cpp\n"
              << "  Function:    get_geodesic_distance()\n"
              << "  Scenario:    one interface sits exactly at the origin, so its\n"
              << "               magnitude is 0 and the cosine is computed as\n"
              << "               0.0 / 0.0.\n"
              << "  Criterion:   the function has no guard for this case, so the\n"
              << "               IEEE result is NaN.  The call is still safe (no\n"
              << "               exit()/abort()), so we document the behaviour\n"
              << "               rather than treating it as an error.\n";

    const double dist = ggd_call_and_report("origin vs unit x", Coord { 0.0, 0.0, 0.0 },
        Coord { 1.0, 0.0, 0.0 });

    EXPECT_TRUE(std::isnan(dist))
        << "A zero-magnitude interface makes the cosine 0/0, which propagates NaN through acos()";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so that a
// failure in one case does not prevent the remaining cases from running.
// -----------------------------------------------------------------------------
TEST(GetGeodesicDistance, ZeroAngleGivesZeroDistance) { test_ggd_zero_angle_gives_zero_distance(); }
TEST(GetGeodesicDistance, OrthogonalPointsOnCommonSphere) { test_ggd_orthogonal_points_on_common_sphere(); }
TEST(GetGeodesicDistance, AntipodalPoints) { test_ggd_antipodal_points(); }
TEST(GetGeodesicDistance, FortyFiveDegreeSeparation) { test_ggd_forty_five_degree_separation(); }
TEST(GetGeodesicDistance, UsesArithmeticMeanRadius) { test_ggd_uses_arithmetic_mean_radius(); }
TEST(GetGeodesicDistance, SymmetryInArguments) { test_ggd_symmetry_in_arguments(); }
TEST(GetGeodesicDistance, ScalesLinearlyWithRadius) { test_ggd_scales_linearly_with_radius(); }
TEST(GetGeodesicDistance, DoesNotModifyInputs) { test_ggd_does_not_modify_inputs(); }
TEST(GetGeodesicDistance, ResultIsBoundedAndNonNegative) { test_ggd_result_is_bounded_and_non_negative(); }
TEST(GetGeodesicDistance, ZeroVectorYieldsNaN) { test_ggd_zero_vector_yields_nan(); }