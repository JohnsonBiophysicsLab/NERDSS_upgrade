/*! \file test_create_arbitrary_vector.cpp
 *
 * ### Unit test for src/reactions/create_arbitrary_vector.cpp
 *
 * Function under test:
 * \code
 *     Vector create_arbitrary_vector(Vector& vec);
 * \endcode
 *
 * Behaviour of the function (from the source):
 *   1. It builds the two cartesian reference axes x=(1,0,0) and y=(0,1,0).
 *   2. It **normalizes the input vector in place** (note: `vec` is taken by
 *      non-const reference, so the caller's vector is modified).
 *   3. It returns the cross product of the (now normalized) input vector with
 *      the x-axis, *unless* the input vector is parallel (angle == 0) or
 *      antiparallel (angle == pi) to the x-axis, in which case the cross
 *      product is taken with the y-axis instead.
 *
 * The purpose of the "unless" branch is to guarantee that the returned vector
 * is never the degenerate zero vector: crossing a vector with an axis it is
 * co-linear with would produce (0,0,0).
 *
 * So the invariants that every test below checks are:
 *   - the input vector has unit length after the call (side effect),
 *   - the returned vector is orthogonal to the (normalized) input vector,
 *   - the returned vector is non-degenerate (non-zero magnitude).
 *
 * Verbose progress information is written to std::cerr so the reader can
 * follow which function is being exercised and what each assertion checks.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with cav_ to avoid collisions with other test files).
// -----------------------------------------------------------------------------
namespace {

//! Numerical tolerance used for all floating point comparisons.
constexpr double kCavTol = 1e-12;

/*! \brief Compute the Euclidean length of a Vector without mutating it.
 *
 * Coord::get_magnitude() is a non-const member function, so we compute the
 * magnitude by hand here to keep the tested objects untouched.
 */
double cav_magnitude(const Vector& v)
{
    return std::sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

/*! \brief Plain dot product of two Vectors (no mutation, no magnitude fields). */
double cav_dot(const Vector& a, const Vector& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

/*! \brief Pretty-print a Vector to std::cerr for the verbose log. */
void cav_print(const char* label, const Vector& v)
{
    std::cerr << "    " << label << " = (" << v.x << ", " << v.y << ", " << v.z
              << "), |v| = " << cav_magnitude(v) << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: A general vector (not co-linear with the x-axis).
//         Expected: cross product with the x-axis is returned.
// -----------------------------------------------------------------------------
void test_cav_general_vector_uses_x_axis()
{
    std::cerr << "\n[TEST] test_cav_general_vector_uses_x_axis\n"
              << "  Source file:   src/reactions/create_arbitrary_vector.cpp\n"
              << "  Function:      create_arbitrary_vector(Vector&)\n"
              << "  Scenario:      input is the z-axis (0,0,1), i.e. 90 degrees\n"
              << "                 away from the x-axis so the x-axis branch is taken.\n"
              << "  Pass criteria: result == cross((0,0,1),(1,0,0)) == (0,1,0),\n"
              << "                 result is orthogonal to the input and non-zero.\n";

    Vector vec(0.0, 0.0, 1.0);
    cav_print("input  vec", vec);

    Vector result = create_arbitrary_vector(vec);
    cav_print("output res", result);
    cav_print("vec after ", vec);

    // cross((0,0,1),(1,0,0)) = (0*0 - 1*0, 1*1 - 0*0, 0*0 - 0*1) = (0,1,0)
    EXPECT_NEAR(result.x, 0.0, kCavTol) << "x component of the cross product should be 0";
    EXPECT_NEAR(result.y, 1.0, kCavTol) << "y component of the cross product should be 1";
    EXPECT_NEAR(result.z, 0.0, kCavTol) << "z component of the cross product should be 0";

    // The returned vector must be orthogonal to the input.
    EXPECT_NEAR(cav_dot(result, vec), 0.0, kCavTol)
        << "returned vector must be orthogonal to the input vector";

    // ...and must never be the degenerate zero vector.
    EXPECT_GT(cav_magnitude(result), kCavTol)
        << "returned vector must be non-degenerate (non-zero length)";
}

// -----------------------------------------------------------------------------
// Test 2: The function normalizes its (by-reference) argument in place.
// -----------------------------------------------------------------------------
void test_cav_normalizes_input_in_place()
{
    std::cerr << "\n[TEST] test_cav_normalizes_input_in_place\n"
              << "  Source file:   src/reactions/create_arbitrary_vector.cpp\n"
              << "  Function:      create_arbitrary_vector(Vector&)\n"
              << "  Scenario:      pass the non-unit vector (3,4,0) (length 5).\n"
              << "  Pass criteria: after the call vec == (0.6,0.8,0) with unit length,\n"
              << "                 documenting the in-place normalization side effect.\n";

    Vector vec(3.0, 4.0, 0.0);
    std::cerr << "  Input length before call = " << cav_magnitude(vec)
              << " (expected 5)\n";
    EXPECT_NEAR(cav_magnitude(vec), 5.0, kCavTol)
        << "sanity check: the input vector should start with length 5";

    Vector result = create_arbitrary_vector(vec);
    cav_print("vec after ", vec);
    cav_print("output res", result);

    // The caller's vector is normalized by the function.
    EXPECT_NEAR(vec.x, 0.6, kCavTol) << "normalized x should be 3/5 = 0.6";
    EXPECT_NEAR(vec.y, 0.8, kCavTol) << "normalized y should be 4/5 = 0.8";
    EXPECT_NEAR(vec.z, 0.0, kCavTol) << "normalized z should be 0";
    EXPECT_NEAR(cav_magnitude(vec), 1.0, kCavTol)
        << "the input vector must have unit length after the call";

    // (3,4,0) is not co-linear with x, so cross with the x-axis is used:
    // cross((0.6,0.8,0),(1,0,0)) = (0.8*0-0*0, 0*1-0.6*0, 0.6*0-0.8*1) = (0,0,-0.8)
    EXPECT_NEAR(result.x, 0.0, kCavTol) << "expected cross product x component 0";
    EXPECT_NEAR(result.y, 0.0, kCavTol) << "expected cross product y component 0";
    EXPECT_NEAR(result.z, -0.8, kCavTol) << "expected cross product z component -0.8";

    // Orthogonality invariant.
    EXPECT_NEAR(cav_dot(result, vec), 0.0, kCavTol)
        << "returned vector must be orthogonal to the normalized input";
}

// -----------------------------------------------------------------------------
// Test 3: Input parallel to the x-axis -> must fall back to the y-axis so the
//         result is not the zero vector.
// -----------------------------------------------------------------------------
void test_cav_parallel_to_x_axis_falls_back_to_y()
{
    std::cerr << "\n[TEST] test_cav_parallel_to_x_axis_falls_back_to_y\n"
              << "  Source file:   src/reactions/create_arbitrary_vector.cpp\n"
              << "  Function:      create_arbitrary_vector(Vector&)\n"
              << "  Scenario:      input (2,0,0) is parallel to the x-axis (angle 0),\n"
              << "                 which would make cross(vec, x_axis) == (0,0,0).\n"
              << "  Pass criteria: the y-axis branch is used instead, giving\n"
              << "                 cross((1,0,0),(0,1,0)) == (0,0,1) (non-degenerate).\n";

    Vector vec(2.0, 0.0, 0.0);
    cav_print("input  vec", vec);

    Vector result = create_arbitrary_vector(vec);
    cav_print("vec after ", vec);
    cav_print("output res", result);

    // Input is normalized to the pure x unit vector.
    EXPECT_NEAR(vec.x, 1.0, kCavTol) << "normalized input should be the x unit vector";
    EXPECT_NEAR(vec.y, 0.0, kCavTol) << "normalized input y should be 0";
    EXPECT_NEAR(vec.z, 0.0, kCavTol) << "normalized input z should be 0";

    // The critical property: the result must NOT be the zero vector.
    EXPECT_GT(cav_magnitude(result), kCavTol)
        << "the y-axis fallback must avoid returning a degenerate zero vector";

    // cross((1,0,0),(0,1,0)) = (0,0,1)
    EXPECT_NEAR(result.x, 0.0, kCavTol) << "expected result x component 0";
    EXPECT_NEAR(result.y, 0.0, kCavTol) << "expected result y component 0";
    EXPECT_NEAR(result.z, 1.0, kCavTol) << "expected result z component 1";

    // Still orthogonal to the input.
    EXPECT_NEAR(cav_dot(result, vec), 0.0, kCavTol)
        << "returned vector must be orthogonal to the input vector";
}

// -----------------------------------------------------------------------------
// Test 4: Input antiparallel to the x-axis (angle == pi) -> also uses the
//         y-axis branch.
// -----------------------------------------------------------------------------
void test_cav_antiparallel_to_x_axis_falls_back_to_y()
{
    std::cerr << "\n[TEST] test_cav_antiparallel_to_x_axis_falls_back_to_y\n"
              << "  Source file:   src/reactions/create_arbitrary_vector.cpp\n"
              << "  Function:      create_arbitrary_vector(Vector&)\n"
              << "  Scenario:      input (-3,0,0) is antiparallel to the x-axis\n"
              << "                 (angle == pi), the second degenerate case.\n"
              << "  Pass criteria: result == cross((-1,0,0),(0,1,0)) == (0,0,-1),\n"
              << "                 i.e. non-zero and orthogonal to the input.\n";

    Vector vec(-3.0, 0.0, 0.0);
    cav_print("input  vec", vec);

    Vector result = create_arbitrary_vector(vec);
    cav_print("vec after ", vec);
    cav_print("output res", result);

    EXPECT_NEAR(vec.x, -1.0, kCavTol) << "normalized input should be the -x unit vector";
    EXPECT_NEAR(cav_magnitude(vec), 1.0, kCavTol) << "input must be unit length after call";

    // cross((-1,0,0),(0,1,0)) = (0,0,-1)
    EXPECT_NEAR(result.x, 0.0, kCavTol) << "expected result x component 0";
    EXPECT_NEAR(result.y, 0.0, kCavTol) << "expected result y component 0";
    EXPECT_NEAR(result.z, -1.0, kCavTol) << "expected result z component -1";

    EXPECT_GT(cav_magnitude(result), kCavTol)
        << "antiparallel case must also avoid a degenerate zero vector";
    EXPECT_NEAR(cav_dot(result, vec), 0.0, kCavTol)
        << "returned vector must be orthogonal to the input vector";
}

// -----------------------------------------------------------------------------
// Test 5: Sweep a batch of arbitrary vectors and verify the general contract:
//         input normalized, output orthogonal, output non-degenerate, and the
//         output length equal to |sin(angle to the chosen axis)|.
// -----------------------------------------------------------------------------
void test_cav_orthogonality_contract_for_many_vectors()
{
    std::cerr << "\n[TEST] test_cav_orthogonality_contract_for_many_vectors\n"
              << "  Source file:   src/reactions/create_arbitrary_vector.cpp\n"
              << "  Function:      create_arbitrary_vector(Vector&)\n"
              << "  Scenario:      loop over a set of arbitrary (non-degenerate)\n"
              << "                 vectors of different lengths and orientations.\n"
              << "  Pass criteria: for every input -> |vec| == 1 after the call,\n"
              << "                 dot(result, vec) == 0, |result| > 0, and\n"
              << "                 |result| == |sin(theta)| w.r.t. the chosen axis.\n";

    // A handful of arbitrary directions: none of them is co-linear with x, so
    // all of them should take the cross-with-x-axis branch.
    const std::vector<Vector> inputs {
        Vector(1.0, 1.0, 1.0),
        Vector(-2.0, 5.0, 0.5),
        Vector(0.0, 1.0, 0.0), // the y-axis itself (90 deg from x)
        Vector(0.1, -0.2, 0.3),
        Vector(7.0, 0.0, 7.0),
        Vector(-1.0, -1.0, 4.0),
    };

    for (std::size_t i = 0; i < inputs.size(); ++i) {
        Vector vec = inputs[i]; // working copy, will be normalized in place
        std::cerr << "  --- case " << i << " ---\n";
        cav_print("input  vec", vec);

        Vector result = create_arbitrary_vector(vec);
        cav_print("vec after ", vec);
        cav_print("output res", result);

        // 1. Side effect: the input becomes a unit vector.
        EXPECT_NEAR(cav_magnitude(vec), 1.0, 1e-10)
            << "case " << i << ": input vector must be normalized in place";

        // 2. Orthogonality: the cross product is perpendicular to its operands.
        EXPECT_NEAR(cav_dot(result, vec), 0.0, 1e-10)
            << "case " << i << ": result must be orthogonal to the input";

        // 3. Non-degeneracy: the whole point of the function.
        EXPECT_GT(cav_magnitude(result), 1e-10)
            << "case " << i << ": result must not be the zero vector";

        // 4. Because none of these inputs is co-linear with x, the result is
        //    also perpendicular to the x-axis, and its length is |sin(theta)|
        //    where theta is the angle between the (unit) input and the x-axis.
        const Vector xAxis(1.0, 0.0, 0.0);
        EXPECT_NEAR(cav_dot(result, xAxis), 0.0, 1e-10)
            << "case " << i << ": result must be orthogonal to the x-axis too";

        const double cosTheta = cav_dot(vec, xAxis); // vec is a unit vector here
        const double expectedLen = std::sqrt(std::max(0.0, 1.0 - (cosTheta * cosTheta)));
        std::cerr << "    cos(theta) w.r.t x-axis = " << cosTheta
                  << ", expected |result| = " << expectedLen << '\n';
        EXPECT_NEAR(cav_magnitude(result), expectedLen, 1e-10)
            << "case " << i << ": |a x b| should equal |a||b|sin(theta) with |a|=|b|=1";
    }
}

// -----------------------------------------------------------------------------
// Test 6: A vector that is *almost* parallel to the x-axis still takes the
//         x-axis branch (the guard is an exact 0/pi comparison), producing a
//         very short but still non-zero vector. This documents a known
//         numerical limitation of the implementation.
// -----------------------------------------------------------------------------
void test_cav_nearly_parallel_vector_still_nonzero()
{
    std::cerr << "\n[TEST] test_cav_nearly_parallel_vector_still_nonzero\n"
              << "  Source file:   src/reactions/create_arbitrary_vector.cpp\n"
              << "  Function:      create_arbitrary_vector(Vector&)\n"
              << "  Scenario:      input (1, 1e-6, 0) is nearly (but not exactly)\n"
              << "                 parallel to the x-axis, so the exact 0/pi guard\n"
              << "                 does not trigger and the x-axis branch is used.\n"
              << "  Pass criteria: the returned vector is still non-zero and\n"
              << "                 orthogonal to the input, although very short.\n";

    Vector vec(1.0, 1.0e-6, 0.0);
    cav_print("input  vec", vec);

    Vector result = create_arbitrary_vector(vec);
    cav_print("vec after ", vec);
    cav_print("output res", result);

    // Input normalized as usual.
    EXPECT_NEAR(cav_magnitude(vec), 1.0, 1e-12)
        << "input must still be normalized in the near-degenerate case";

    // Non-zero, though tiny (~1e-6): documents the numerical sensitivity.
    EXPECT_GT(cav_magnitude(result), 0.0)
        << "result should still be non-zero for a nearly-parallel input";
    EXPECT_LT(cav_magnitude(result), 1e-3)
        << "result is expected to be very short for a nearly-parallel input";

    // Orthogonality still holds (up to the small magnitude involved).
    EXPECT_NEAR(cav_dot(result, vec), 0.0, 1e-12)
        << "result must remain orthogonal to the input";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper runs inside its own TEST so
// that a failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(CreateArbitraryVector, GeneralVectorUsesXAxis) { test_cav_general_vector_uses_x_axis(); }
TEST(CreateArbitraryVector, NormalizesInputInPlace) { test_cav_normalizes_input_in_place(); }
TEST(CreateArbitraryVector, ParallelToXAxisFallsBackToY) { test_cav_parallel_to_x_axis_falls_back_to_y(); }
TEST(CreateArbitraryVector, AntiparallelToXAxisFallsBackToY) { test_cav_antiparallel_to_x_axis_falls_back_to_y(); }
TEST(CreateArbitraryVector, OrthogonalityContractForManyVectors) { test_cav_orthogonality_contract_for_many_vectors(); }
TEST(CreateArbitraryVector, NearlyParallelVectorStillNonzero) { test_cav_nearly_parallel_vector_still_nonzero(); }