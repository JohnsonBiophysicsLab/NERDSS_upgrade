/*! \file test_requiresSignFlip.cpp
 *
 * ### Unit test for src/reactions/requiresSignFlip.cpp
 *
 * The single function under test is
 *
 *     bool requiresSignFlip(Vector axis, Vector v1, Vector v2)
 *
 * Its job is to decide, for two vectors `v1` and `v2` that are looked at down a
 * given rotation `axis`, whether the sign of the angle between them has to be
 * flipped.  Internally it
 *
 *   1. builds the rotation that lines `axis` up with a cartesian axis,
 *      normally the z axis (`u = z x axis`, `theta = angle(z, axis)`),
 *   2. falls back to the x axis when `z x axis` is the zero vector, i.e. when
 *      `axis` is parallel to +z or -z,
 *   3. rotates `v1`, `v2` (and `axis`) with that rotation, possibly reversing
 *      it if the first attempt turned `axis` the wrong way,
 *   4. orthographically projects `v1` and `v2` onto the plane normal to the
 *      chosen cartesian axis and returns whether the cross product of the two
 *      projections points along the positive cartesian axis.
 *
 * Behaviour that the tests below pin down (derived by hand from the source and
 * from Vector/Quat in class_Vector.cpp / class_Quat.cpp):
 *
 *   * For an `axis` that is NOT parallel to z the net rotation applied to
 *     `v1`/`v2` maps `axis` onto +z, therefore the return value is exactly
 *     `((v1 x v2) . axis) > 0`.
 *   * For an `axis` that IS parallel to +z or -z the x-axis fallback is used and
 *     the resulting sign convention is INVERTED, i.e. the return value is
 *     `((v1 x v2) . axis) < 0`.  These tests document that (surprising) branch
 *     rather than "fixing" it.
 *   * `Vector::dot_theta()` needs the cached `magnitude` member, so callers must
 *     call `calc_magnitude()` on `axis` before calling the function.  A separate
 *     test documents what happens if they forget.
 *   * `axis`, `v1` and `v2` are taken BY VALUE, so the caller's vectors are not
 *     modified.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Small helpers used by all of the tests below.  They carry the `rsf_` prefix so
// they cannot collide with anything else linked into the test binary.
// -----------------------------------------------------------------------------

/*! \brief Build a Vector and initialise its cached magnitude.
 *
 * requiresSignFlip() calls Vector::dot_theta() on the axis, and dot_theta()
 * uses the cached `magnitude` member (it does NOT recompute it), so any axis
 * handed to the function must have had calc_magnitude() called on it.
 */
Vector rsf_vec(double x, double y, double z)
{
    Vector v(x, y, z);
    v.calc_magnitude();
    return v;
}

/*! \brief Plain (non-normalised) cross product, used as the independent
 *         reference implementation.
 *
 * Note that Vector::cross() normalises its result; the SIGN of every component
 * is identical either way, which is all these tests care about.
 */
Vector rsf_raw_cross(const Vector& a, const Vector& b)
{
    Vector c(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    c.calc_magnitude();
    return c;
}

/*! \brief Scalar triple product (v1 x v2) . axis -- the right-hand-rule sign. */
double rsf_triple(const Vector& v1, const Vector& v2, const Vector& axis)
{
    const Vector c = rsf_raw_cross(v1, v2);
    return (c.x * axis.x) + (c.y * axis.y) + (c.z * axis.z);
}

/*! \brief Pretty printer for console output. */
std::string rsf_str(const Vector& v)
{
    std::ostringstream oss;
    oss << '(' << v.x << ", " << v.y << ", " << v.z << ')';
    return oss.str();
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: general axis (not parallel to z), pair built so the answer is known.
//
// For any axis and any v1 not parallel to it, choosing v2 = axis x v1 gives
//     (v1 x v2) . axis = |v1|^2 |axis|^2 - (v1.axis)^2 > 0,
// so the function must return true, and false with the arguments swapped.
// -----------------------------------------------------------------------------
void test_rsf_general_axis_right_hand_rule()
{
    std::cerr << "\n[requiresSignFlip.cpp] Testing requiresSignFlip() -- general axis, constructed pair\n"
              << "  Setup:         v2 = axis x v1, so (v1 x v2).axis is guaranteed positive.\n"
              << "  Pass criteria: requiresSignFlip(axis, v1, v2) == true and\n"
              << "                 requiresSignFlip(axis, v2, v1) == false.\n";

    // {axis, v1}.  Every axis here is comfortably away from +/-z so the generic
    // (z-axis) branch of the function is exercised.
    struct Case {
        double ax, ay, az;
        double vx, vy, vz;
    };
    const std::vector<Case> cases {
        { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, //  +x axis
        { -1.0, 0.0, 0.0, 0.0, 0.0, 1.0 }, //  -x axis
        { 0.0, 1.0, 0.0, 1.0, 0.0, 0.0 }, //  +y axis
        { 0.0, -1.0, 0.0, 0.0, 0.0, 1.0 }, //  -y axis
        { 1.0, 1.0, 0.0, 0.0, 0.0, 1.0 }, //  diagonal in the xy plane
        { 1.0, 2.0, 3.0, 1.0, 0.0, 0.0 }, //  generic oblique axis
        { 2.0, -1.0, 1.0, 0.0, 1.0, 0.0 }, //  generic oblique axis
        { -3.0, 1.0, -1.0, 1.0, 1.0, 1.0 }, //  generic oblique axis, negative z part
        { 2.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, //  NON-unit axis (magnitude must not matter)
    };

    for (const auto& c : cases) {
        const Vector axis = rsf_vec(c.ax, c.ay, c.az);
        const Vector v1 = rsf_vec(c.vx, c.vy, c.vz);
        const Vector v2 = rsf_raw_cross(axis, v1); // by construction (v1 x v2).axis > 0

        const double triple = rsf_triple(v1, v2, axis);

        std::cerr << "  -> axis " << rsf_str(axis) << ", v1 " << rsf_str(v1) << ", v2 = axis x v1 "
                  << rsf_str(v2) << ", (v1 x v2).axis = " << triple << '\n';

        // Sanity check on the test setup itself, not on the code under test.
        EXPECT_GT(triple, 0.0) << "test setup error: the constructed pair should have a positive triple product";

        const bool forward = requiresSignFlip(axis, v1, v2);
        const bool swapped = requiresSignFlip(axis, v2, v1);
        std::cerr << "     requiresSignFlip(axis, v1, v2) = " << std::boolalpha << forward
                  << ", requiresSignFlip(axis, v2, v1) = " << swapped << '\n';

        EXPECT_TRUE(forward) << "positive triple product (v1 x v2).axis must give true for axis " << rsf_str(axis);
        EXPECT_FALSE(swapped) << "swapping the two vectors must flip the answer for axis " << rsf_str(axis);
    }
}

// -----------------------------------------------------------------------------
// Test 2: general axis, arbitrary (non-orthogonally-constructed) vector pairs.
//
// Verifies the full rule: result == ((v1 x v2) . axis > 0).
// -----------------------------------------------------------------------------
void test_rsf_general_axis_matches_triple_product()
{
    std::cerr << "\n[requiresSignFlip.cpp] Testing requiresSignFlip() -- general axis, arbitrary pairs\n"
              << "  Pass criteria: result == ((v1 x v2) . axis > 0) for every non-degenerate case.\n";

    struct Case {
        double ax, ay, az;
        double v1x, v1y, v1z;
        double v2x, v2y, v2z;
    };
    const std::vector<Case> cases {
        { 1.0, 2.0, 3.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, // triple = +3
        { 1.0, 2.0, 3.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0 }, // triple = -3
        { 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0 }, // triple = +1
        { 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0 }, // triple = -1
        { 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 }, // triple = +1
        { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, // triple = +1
        { 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0 }, // triple = -1
    };

    for (const auto& c : cases) {
        const Vector axis = rsf_vec(c.ax, c.ay, c.az);
        const Vector v1 = rsf_vec(c.v1x, c.v1y, c.v1z);
        const Vector v2 = rsf_vec(c.v2x, c.v2y, c.v2z);

        const double triple = rsf_triple(v1, v2, axis);
        const bool expected = (triple > 0.0);
        const bool actual = requiresSignFlip(axis, v1, v2);

        std::cerr << "  -> axis " << rsf_str(axis) << ", v1 " << rsf_str(v1) << ", v2 " << rsf_str(v2)
                  << " : triple = " << triple << ", expected " << std::boolalpha << expected << ", got "
                  << actual << '\n';

        // Guard against accidentally adding a degenerate case to the table.
        EXPECT_NE(triple, 0.0) << "test setup error: degenerate (co-planar) case in the table";
        EXPECT_EQ(actual, expected)
            << "for a non-z axis the result must be the right-hand-rule sign of (v1 x v2).axis";
    }
}

// -----------------------------------------------------------------------------
// Test 3: axis parallel to +z or -z -> the x-axis fallback branch.
//
// Here `z x axis` is the zero vector, so the code switches to `x x axis`.  In
// that branch the returned sign is the OPPOSITE of the right-hand-rule sign,
// i.e. result == ((v1 x v2) . axis < 0).  This test documents that behaviour.
// -----------------------------------------------------------------------------
void test_rsf_z_aligned_axis_uses_x_axis_fallback()
{
    std::cerr << "\n[requiresSignFlip.cpp] Testing requiresSignFlip() -- axis parallel to z (x-axis fallback)\n"
              << "  Note:          z x axis == 0 here, so the routine falls back to the x axis and the\n"
              << "                 returned sign convention is INVERTED with respect to the general case.\n"
              << "  Pass criteria: result == ((v1 x v2) . axis < 0).\n";

    struct Case {
        double ax, ay, az;
        double v1x, v1y, v1z;
        double v2x, v2y, v2z;
    };
    const std::vector<Case> cases {
        { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, // +z, x then y
        { 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0 }, // +z, y then x
        { 0.0, 0.0, 5.0, 1.0, 1.0, 0.0, -1.0, 1.0, 0.0 }, // +z, non-unit axis
        { 0.0, 0.0, 1.0, 1.0, 0.0, 5.0, 0.0, 1.0, -2.0 }, // +z, vectors with z parts
        { 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, // -z, x then y
        { 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0 }, // -z, y then x
    };

    for (const auto& c : cases) {
        const Vector axis = rsf_vec(c.ax, c.ay, c.az);
        const Vector v1 = rsf_vec(c.v1x, c.v1y, c.v1z);
        const Vector v2 = rsf_vec(c.v2x, c.v2y, c.v2z);

        const double triple = rsf_triple(v1, v2, axis);
        const bool expected = (triple < 0.0); // inverted convention in this branch
        const bool actual = requiresSignFlip(axis, v1, v2);

        std::cerr << "  -> axis " << rsf_str(axis) << ", v1 " << rsf_str(v1) << ", v2 " << rsf_str(v2)
                  << " : triple = " << triple << ", expected " << std::boolalpha << expected << ", got "
                  << actual << '\n';

        EXPECT_NE(triple, 0.0) << "test setup error: degenerate (co-planar) case in the table";
        EXPECT_EQ(actual, expected)
            << "for an axis parallel to z the x-axis fallback inverts the sign convention";
    }
}

// -----------------------------------------------------------------------------
// Test 4: swapping the two vectors always flips the answer (antisymmetry).
// -----------------------------------------------------------------------------
void test_rsf_swapping_vectors_flips_result()
{
    std::cerr << "\n[requiresSignFlip.cpp] Testing antisymmetry of requiresSignFlip()\n"
              << "  Pass criteria: requiresSignFlip(axis, v1, v2) != requiresSignFlip(axis, v2, v1)\n"
              << "                 for every non-degenerate (non-coplanar-with-axis) pair.\n";

    struct Case {
        double ax, ay, az;
        double v1x, v1y, v1z;
        double v2x, v2y, v2z;
    };
    const std::vector<Case> cases {
        { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 }, // simple x axis
        { 1.0, 2.0, 3.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, // oblique axis
        { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, // z axis (fallback branch)
        { 0.0, 0.0, -1.0, 1.0, 1.0, 0.0, -1.0, 1.0, 0.0 }, // -z axis (fallback branch)
    };

    for (const auto& c : cases) {
        const Vector axis = rsf_vec(c.ax, c.ay, c.az);
        const Vector v1 = rsf_vec(c.v1x, c.v1y, c.v1z);
        const Vector v2 = rsf_vec(c.v2x, c.v2y, c.v2z);

        const bool forward = requiresSignFlip(axis, v1, v2);
        const bool swapped = requiresSignFlip(axis, v2, v1);

        std::cerr << "  -> axis " << rsf_str(axis) << " : forward = " << std::boolalpha << forward
                  << ", swapped = " << swapped << '\n';

        EXPECT_NE(forward, swapped) << "the two orderings must give opposite answers for axis " << rsf_str(axis);
    }
}

// -----------------------------------------------------------------------------
// Test 5: parallel (or identical) vectors have a zero projected cross product,
//         so the strict `> 0` comparison reports false in either order.
//
// Using v2 == v1 and v2 == 2*v1 keeps the arithmetic bit-exact: identical (or
// exactly power-of-two scaled) inputs go through identical rotations, so the
// cross product components cancel exactly to 0.0.
// -----------------------------------------------------------------------------
void test_rsf_parallel_vectors_return_false()
{
    std::cerr << "\n[requiresSignFlip.cpp] Testing requiresSignFlip() with parallel vectors\n"
              << "  Pass criteria: the projected cross product is exactly zero, and because the code\n"
              << "                 tests `> 0` the function returns false in both orderings.\n";

    struct Case {
        double ax, ay, az;
        double vx, vy, vz;
        double scale; // v2 = scale * v1
    };
    const std::vector<Case> cases {
        { 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 }, // identical vectors, x axis
        { 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0 }, // parallel vectors, x axis
        { 1.0, 2.0, 3.0, 1.0, 0.0, 0.0, 4.0 }, // parallel vectors, oblique axis
        { 0.0, 0.0, 1.0, 1.0, 2.0, 0.0, 2.0 }, // parallel vectors, z axis (fallback branch)
    };

    for (const auto& c : cases) {
        const Vector axis = rsf_vec(c.ax, c.ay, c.az);
        const Vector v1 = rsf_vec(c.vx, c.vy, c.vz);
        const Vector v2 = rsf_vec(c.vx * c.scale, c.vy * c.scale, c.vz * c.scale);

        const bool forward = requiresSignFlip(axis, v1, v2);
        const bool swapped = requiresSignFlip(axis, v2, v1);

        std::cerr << "  -> axis " << rsf_str(axis) << ", v1 " << rsf_str(v1) << ", v2 " << rsf_str(v2)
                  << " : forward = " << std::boolalpha << forward << ", swapped = " << swapped << '\n';

        EXPECT_FALSE(forward) << "parallel vectors give a zero cross product, which is not > 0";
        EXPECT_FALSE(swapped) << "parallel vectors give a zero cross product in either order";
    }
}

// -----------------------------------------------------------------------------
// Test 6: components along the axis are projected away and must not change the
//         answer.
// -----------------------------------------------------------------------------
void test_rsf_component_along_axis_is_ignored()
{
    std::cerr << "\n[requiresSignFlip.cpp] Testing that components parallel to the axis are ignored\n"
              << "  Pass criteria: adding any multiple of `axis` to v1 and/or v2 leaves the result\n"
              << "                 unchanged, because the routine projects onto the plane normal to axis.\n";

    const Vector axis = rsf_vec(1.0, 0.0, 0.0);
    const Vector v1 = rsf_vec(0.0, 1.0, 0.0);
    const Vector v2 = rsf_vec(0.0, 0.0, 1.0);

    // Same perpendicular parts, but now with large components along +/- axis.
    const Vector v1Shifted = rsf_vec(3.0, 1.0, 0.0);
    const Vector v2Shifted = rsf_vec(-7.0, 0.0, 1.0);

    const bool base = requiresSignFlip(axis, v1, v2);
    const bool shiftedFirst = requiresSignFlip(axis, v1Shifted, v2);
    const bool shiftedBoth = requiresSignFlip(axis, v1Shifted, v2Shifted);

    std::cerr << "  -> axis " << rsf_str(axis) << '\n'
              << "     base   : v1 " << rsf_str(v1) << ", v2 " << rsf_str(v2) << " -> " << std::boolalpha
              << base << '\n'
              << "     shift1 : v1 " << rsf_str(v1Shifted) << ", v2 " << rsf_str(v2) << " -> " << shiftedFirst
              << '\n'
              << "     shift2 : v1 " << rsf_str(v1Shifted) << ", v2 " << rsf_str(v2Shifted) << " -> "
              << shiftedBoth << '\n';

    // (y x z).x = +1 > 0, so the baseline answer is true for this configuration.
    EXPECT_TRUE(base) << "(v1 x v2).axis = +1 for the baseline configuration, so the result must be true";
    EXPECT_EQ(shiftedFirst, base) << "adding a multiple of axis to v1 must not change the answer";
    EXPECT_EQ(shiftedBoth, base) << "adding a multiple of axis to both vectors must not change the answer";
}

// -----------------------------------------------------------------------------
// Test 7: the three arguments are taken by value, so the caller's vectors must
//         come back untouched even though the function rotates its copies.
// -----------------------------------------------------------------------------
void test_rsf_arguments_are_passed_by_value()
{
    std::cerr << "\n[requiresSignFlip.cpp] Testing that requiresSignFlip() does not modify its arguments\n"
              << "  Pass criteria: axis, v1 and v2 (components and cached magnitude) are unchanged after\n"
              << "                 the call, because the signature takes them by value.\n";

    Vector axis = rsf_vec(1.0, 2.0, 3.0);
    Vector v1 = rsf_vec(1.0, 0.0, 0.0);
    Vector v2 = rsf_vec(0.0, 1.0, 0.0);

    // Remember every field we care about.
    const double ax = axis.x, ay = axis.y, az = axis.z, am = axis.magnitude;
    const double b1x = v1.x, b1y = v1.y, b1z = v1.z, b1m = v1.magnitude;
    const double b2x = v2.x, b2y = v2.y, b2z = v2.z, b2m = v2.magnitude;

    const bool result = requiresSignFlip(axis, v1, v2);
    std::cerr << "  -> call returned " << std::boolalpha << result << "; now checking the inputs...\n";

    EXPECT_DOUBLE_EQ(axis.x, ax) << "axis.x must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(axis.y, ay) << "axis.y must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(axis.z, az) << "axis.z must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(axis.magnitude, am) << "axis.magnitude must be unchanged (pass by value)";

    EXPECT_DOUBLE_EQ(v1.x, b1x) << "v1.x must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(v1.y, b1y) << "v1.y must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(v1.z, b1z) << "v1.z must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(v1.magnitude, b1m) << "v1.magnitude must be unchanged (pass by value)";

    EXPECT_DOUBLE_EQ(v2.x, b2x) << "v2.x must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(v2.y, b2y) << "v2.y must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(v2.z, b2z) << "v2.z must be unchanged (pass by value)";
    EXPECT_DOUBLE_EQ(v2.magnitude, b2m) << "v2.magnitude must be unchanged (pass by value)";
}

// -----------------------------------------------------------------------------
// Test 8: the axis must have its cached magnitude initialised.
//
// Vector::dot_theta() reads the cached `magnitude` member; when it is still the
// default 0 it prints a warning and returns 0.0.  In that case the internal
// rotation collapses to the identity, the projection is taken in the raw xy
// plane instead of the plane normal to `axis`, and the answer is generally
// wrong.  This test documents that precondition rather than the "right" answer.
// -----------------------------------------------------------------------------
void test_rsf_axis_magnitude_precondition()
{
    std::cerr << "\n[requiresSignFlip.cpp] Testing the axis-magnitude precondition\n"
              << "  Note:          Vector::dot_theta() uses the cached magnitude; with magnitude == 0 it\n"
              << "                 prints a warning and returns 0, so the internal rotation becomes the\n"
              << "                 identity and the xy-plane projection is used instead.\n"
              << "  Pass criteria: uninitialised magnitude -> false (projection of v2 onto xy is zero),\n"
              << "                 initialised magnitude   -> true  ((v1 x v2).axis = +1 > 0).\n"
              << "  (One or more 'magnitude 0' warnings from Vector::dot_theta are expected below.)\n";

    Vector axis(1.0, 0.0, 0.0); // magnitude deliberately left at its default 0
    const Vector v1 = rsf_vec(0.0, 1.0, 0.0);
    const Vector v2 = rsf_vec(0.0, 0.0, 1.0);

    const bool uninitialised = requiresSignFlip(axis, v1, v2);
    std::cerr << "  -> with axis.magnitude == 0 the call returned " << std::boolalpha << uninitialised << '\n';
    EXPECT_FALSE(uninitialised)
        << "with no rotation applied the xy projection of v2 = (0,0,1) is the zero vector, giving false";

    axis.calc_magnitude(); // fix the precondition
    const bool initialised = requiresSignFlip(axis, v1, v2);
    std::cerr << "  -> with axis.magnitude == " << axis.magnitude << " the call returned " << initialised
              << '\n';
    EXPECT_TRUE(initialised) << "(y x z).x = +1 > 0, so the correctly initialised axis must give true";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* function runs inside its own TEST so
// that a failure in one does not stop the remaining checks from running.
// -----------------------------------------------------------------------------
TEST(RequiresSignFlipTest, GeneralAxisRightHandRule) { test_rsf_general_axis_right_hand_rule(); }
TEST(RequiresSignFlipTest, GeneralAxisMatchesTripleProduct) { test_rsf_general_axis_matches_triple_product(); }
TEST(RequiresSignFlipTest, ZAlignedAxisUsesXAxisFallback) { test_rsf_z_aligned_axis_uses_x_axis_fallback(); }
TEST(RequiresSignFlipTest, SwappingVectorsFlipsResult) { test_rsf_swapping_vectors_flips_result(); }
TEST(RequiresSignFlipTest, ParallelVectorsReturnFalse) { test_rsf_parallel_vectors_return_false(); }
TEST(RequiresSignFlipTest, ComponentAlongAxisIsIgnored) { test_rsf_component_along_axis_is_ignored(); }
TEST(RequiresSignFlipTest, ArgumentsArePassedByValue) { test_rsf_arguments_are_passed_by_value(); }
TEST(RequiresSignFlipTest, AxisMagnitudePrecondition) { test_rsf_axis_magnitude_precondition(); }