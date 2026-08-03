/*! \file test_angleSignIsCorrect.cpp
 *
 * ### Unit test for src/reactions/angleSignIsCorrect.cpp
 *
 * Function under test:
 *
 *     bool angleSignIsCorrect(const Vector& vec1, const Vector& vec2)
 *
 * Behaviour of the implementation (used to derive the pass criteria):
 *
 *   * If BOTH vectors have |z| > 1e-12, the vectors are projected onto the
 *     x-z plane (their y components are discarded) and the function returns
 *     true when the y component of (proj(vec1) x proj(vec2)) is negative.
 *     For y-zeroed vectors the cross product y component is
 *         y = z1 * x2 - x1 * z2
 *
 *   * Otherwise (at least one vector has an essentially zero z component) the
 *     vectors are projected onto the x-y plane (their z components are
 *     discarded) and the function returns true when the z component of
 *     (proj(vec1) x proj(vec2)) is negative.  For z-zeroed vectors,
 *         z = x1 * y2 - y1 * x2
 *
 *   * Degenerate (parallel / anti-parallel) pairs give a zero cross product
 *     component, which is NOT < 0, so the function must return false.
 *
 * The tests below cover both branches, the y/z component that must be
 * ignored in each branch, the 1e-12 tolerance that selects the branch, and
 * the degenerate zero-cross-product cases.  Verbose progress information is
 * written to stderr so a reader of the test log can see exactly which case is
 * being exercised.
 */

#include "reactions/association/association.hpp"
#include "classes/class_Vector.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

namespace {

/*! \brief Convenience factory so the tests read like the math they check. */
Vector asic_make_vector(double x, double y, double z)
{
    return Vector(x, y, z);
}

/*! \brief Pretty-print a Vector to stderr (Vector's operator<< needs a
 *         non-const reference, so we format the components manually). */
void asic_print_vector(const char* label, const Vector& vec)
{
    std::cerr << "    " << label << " = (" << vec.x << ", " << vec.y << ", "
              << vec.z << ")\n";
}

/*! \brief Independent re-implementation of the expected result.
 *
 * This mirrors the documented contract of angleSignIsCorrect() and is used to
 * cross-check the production function over a small sweep of inputs.
 */
bool asic_expected_result(const Vector& v1, const Vector& v2)
{
    if (std::abs(v1.z) > 1E-12 && std::abs(v2.z) > 1E-12) {
        // x-z plane projection: cross product y component.
        return (v1.z * v2.x - v1.x * v2.z) < 0.0;
    }
    // x-y plane projection: cross product z component.
    return (v1.x * v2.y - v1.y * v2.x) < 0.0;
}

} // namespace

// -----------------------------------------------------------------------------
// Branch 1: both vectors have a significant z component -> x-z projection.
// -----------------------------------------------------------------------------
void asic_test_xz_branch_sign()
{
    std::cerr << "\n[TEST] asic_test_xz_branch_sign\n"
              << "  Source file: src/reactions/angleSignIsCorrect.cpp\n"
              << "  Function:    angleSignIsCorrect (x-z projection branch)\n"
              << "  Criteria:    returns true iff (z1*x2 - x1*z2) < 0\n";

    // vec1 = (1,0,1), vec2 = (0,0,1):  z1*x2 - x1*z2 = 1*0 - 1*1 = -1 < 0 -> true
    Vector v1 = asic_make_vector(1.0, 0.0, 1.0);
    Vector v2 = asic_make_vector(0.0, 0.0, 1.0);
    asic_print_vector("vec1", v1);
    asic_print_vector("vec2", v2);
    std::cerr << "  Expected: true (cross-y = -1)\n";
    EXPECT_TRUE(angleSignIsCorrect(v1, v2))
        << "Negative x-z cross product y component must give true";

    // Swapping the arguments flips the sign of the cross product -> false.
    std::cerr << "  Swapping the arguments should flip the result to false\n";
    EXPECT_FALSE(angleSignIsCorrect(v2, v1))
        << "Positive x-z cross product y component must give false";
}

// -----------------------------------------------------------------------------
// Branch 1 detail: in the x-z branch the y components must be discarded.
// -----------------------------------------------------------------------------
void asic_test_xz_branch_ignores_y()
{
    std::cerr << "\n[TEST] asic_test_xz_branch_ignores_y\n"
              << "  Source file: src/reactions/angleSignIsCorrect.cpp\n"
              << "  Function:    angleSignIsCorrect (x-z projection branch)\n"
              << "  Criteria:    y components must not influence the result\n";

    // Same x and z components as the previous test, but with large, opposite
    // y components added.  The result must be unchanged (true).
    Vector plain1 = asic_make_vector(1.0, 0.0, 1.0);
    Vector plain2 = asic_make_vector(0.0, 0.0, 1.0);
    Vector noisy1 = asic_make_vector(1.0, 25.0, 1.0);
    Vector noisy2 = asic_make_vector(0.0, -13.0, 1.0);

    asic_print_vector("plain1", plain1);
    asic_print_vector("noisy1", noisy1);
    asic_print_vector("plain2", plain2);
    asic_print_vector("noisy2", noisy2);

    const bool plainResult = angleSignIsCorrect(plain1, plain2);
    const bool noisyResult = angleSignIsCorrect(noisy1, noisy2);
    std::cerr << "  plain result = " << std::boolalpha << plainResult
              << ", noisy result = " << noisyResult << "\n";

    EXPECT_TRUE(plainResult) << "Baseline x-z case should be true";
    EXPECT_EQ(plainResult, noisyResult)
        << "Adding y components must not change the x-z branch result";
}

// -----------------------------------------------------------------------------
// Branch 2: at least one z component is (numerically) zero -> x-y projection.
// -----------------------------------------------------------------------------
void asic_test_xy_branch_sign()
{
    std::cerr << "\n[TEST] asic_test_xy_branch_sign\n"
              << "  Source file: src/reactions/angleSignIsCorrect.cpp\n"
              << "  Function:    angleSignIsCorrect (x-y projection branch)\n"
              << "  Criteria:    returns true iff (x1*y2 - y1*x2) < 0\n";

    // x-axis then y-axis: x1*y2 - y1*x2 = 1*1 - 0*0 = +1 -> false
    Vector xAxis = asic_make_vector(1.0, 0.0, 0.0);
    Vector yAxis = asic_make_vector(0.0, 1.0, 0.0);
    asic_print_vector("vec1 (x-axis)", xAxis);
    asic_print_vector("vec2 (y-axis)", yAxis);
    std::cerr << "  Expected: false (cross-z = +1)\n";
    EXPECT_FALSE(angleSignIsCorrect(xAxis, yAxis))
        << "Positive x-y cross product z component must give false";

    // y-axis then x-axis: 0*0 - 1*1 = -1 -> true
    std::cerr << "  Reversed order expected: true (cross-z = -1)\n";
    EXPECT_TRUE(angleSignIsCorrect(yAxis, xAxis))
        << "Negative x-y cross product z component must give true";
}

// -----------------------------------------------------------------------------
// Branch 2 detail: only one vector needs a vanishing z to force the x-y branch.
// -----------------------------------------------------------------------------
void asic_test_xy_branch_selected_when_one_z_is_zero()
{
    std::cerr << "\n[TEST] asic_test_xy_branch_selected_when_one_z_is_zero\n"
              << "  Source file: src/reactions/angleSignIsCorrect.cpp\n"
              << "  Function:    angleSignIsCorrect (branch selection)\n"
              << "  Criteria:    a single zero z component forces the x-y branch\n";

    // If the x-z branch were taken for these vectors we would evaluate
    //   z1*x2 - x1*z2 = 1*1 - 1*0 = +1  -> false
    // The x-y branch instead evaluates
    //   x1*y2 - y1*x2 = 1*1 - 1*1 = 0   -> false as well, so choose values
    // that distinguish the two branches clearly:
    //   vec1 = (1, 1, 1), vec2 = (1, 0, 0)
    //   x-z branch: z1*x2 - x1*z2 = 1*1 - 1*0 = +1 -> false
    //   x-y branch: x1*y2 - y1*x2 = 1*0 - 1*1 = -1 -> true  (this is expected)
    Vector v1 = asic_make_vector(1.0, 1.0, 1.0);
    Vector v2 = asic_make_vector(1.0, 0.0, 0.0); // z == 0 -> x-y branch
    asic_print_vector("vec1", v1);
    asic_print_vector("vec2 (z == 0)", v2);
    std::cerr << "  Expected: true (x-y branch, cross-z = -1)\n";
    EXPECT_TRUE(angleSignIsCorrect(v1, v2))
        << "A zero z component in either vector must select the x-y branch";
}

// -----------------------------------------------------------------------------
// The 1e-12 tolerance decides which branch is taken.
// -----------------------------------------------------------------------------
void asic_test_z_tolerance_boundary()
{
    std::cerr << "\n[TEST] asic_test_z_tolerance_boundary\n"
              << "  Source file: src/reactions/angleSignIsCorrect.cpp\n"
              << "  Function:    angleSignIsCorrect (1e-12 z tolerance)\n"
              << "  Criteria:    |z| <= 1e-12 counts as zero -> x-y branch,\n"
              << "               |z| >  1e-12 -> x-z branch\n";

    // Case A: z magnitudes below the tolerance -> x-y branch is used.
    //   x-y branch: x1*y2 - y1*x2 = 1*1 - 0*0 = +1 -> false
    Vector tiny1 = asic_make_vector(1.0, 0.0, 1e-13);
    Vector tiny2 = asic_make_vector(0.0, 1.0, 1e-13);
    asic_print_vector("tiny1 (|z| = 1e-13)", tiny1);
    asic_print_vector("tiny2 (|z| = 1e-13)", tiny2);
    std::cerr << "  Expected: false (x-y branch, cross-z = +1)\n";
    EXPECT_FALSE(angleSignIsCorrect(tiny1, tiny2))
        << "|z| below 1e-12 must be treated as zero (x-y branch)";

    // Case B: identical x/y layout but z magnitudes above the tolerance ->
    // x-z branch is used instead.
    //   x-z branch: z1*x2 - x1*z2 = 1e-6*0 - 1*1e-6 = -1e-6 -> true
    Vector big1 = asic_make_vector(1.0, 0.0, 1e-6);
    Vector big2 = asic_make_vector(0.0, 1.0, 1e-6);
    asic_print_vector("big1 (|z| = 1e-6)", big1);
    asic_print_vector("big2 (|z| = 1e-6)", big2);
    std::cerr << "  Expected: true (x-z branch, cross-y = -1e-6)\n";
    EXPECT_TRUE(angleSignIsCorrect(big1, big2))
        << "|z| above 1e-12 must select the x-z branch and flip the result";
}

// -----------------------------------------------------------------------------
// Degenerate cases: a zero cross product component is not negative -> false.
// -----------------------------------------------------------------------------
void asic_test_degenerate_pairs_return_false()
{
    std::cerr << "\n[TEST] asic_test_degenerate_pairs_return_false\n"
              << "  Source file: src/reactions/angleSignIsCorrect.cpp\n"
              << "  Function:    angleSignIsCorrect (degenerate input)\n"
              << "  Criteria:    zero cross product component -> false\n";

    // Parallel vectors in the x-z branch: 1*2 - 1*2 = 0 -> false
    Vector par1 = asic_make_vector(1.0, 0.0, 1.0);
    Vector par2 = asic_make_vector(2.0, 0.0, 2.0);
    asic_print_vector("parallel vec1", par1);
    asic_print_vector("parallel vec2", par2);
    EXPECT_FALSE(angleSignIsCorrect(par1, par2))
        << "Parallel vectors give a zero cross product -> false";

    // Anti-parallel vectors in the x-y branch: 1*0 - 0*(-1) = 0 -> false
    Vector anti1 = asic_make_vector(1.0, 0.0, 0.0);
    Vector anti2 = asic_make_vector(-1.0, 0.0, 0.0);
    asic_print_vector("anti-parallel vec1", anti1);
    asic_print_vector("anti-parallel vec2", anti2);
    EXPECT_FALSE(angleSignIsCorrect(anti1, anti2))
        << "Anti-parallel vectors give a zero cross product -> false";

    // Identical vectors always give a zero cross product -> false
    Vector same = asic_make_vector(3.0, -2.0, 4.0);
    asic_print_vector("identical vectors", same);
    EXPECT_FALSE(angleSignIsCorrect(same, same))
        << "A vector crossed with itself has a zero cross product -> false";

    // Zero vector: both components vanish, and z == 0 selects the x-y branch.
    Vector zero = asic_make_vector(0.0, 0.0, 0.0);
    asic_print_vector("zero vector", zero);
    EXPECT_FALSE(angleSignIsCorrect(zero, same))
        << "The zero vector cannot produce a negative cross product component";
    EXPECT_FALSE(angleSignIsCorrect(same, zero))
        << "The zero vector cannot produce a negative cross product component";
}

// -----------------------------------------------------------------------------
// Cross-check the production function against an independent re-implementation
// over a deterministic sweep of inputs (both branches are exercised).
// -----------------------------------------------------------------------------
void asic_test_matches_reference_formula()
{
    std::cerr << "\n[TEST] asic_test_matches_reference_formula\n"
              << "  Source file: src/reactions/angleSignIsCorrect.cpp\n"
              << "  Function:    angleSignIsCorrect (sweep vs. reference)\n"
              << "  Criteria:    result equals the documented sign formula for\n"
              << "               every sampled vector pair\n";

    // A deterministic set of component values covering positive, negative and
    // (numerically) zero z values so both branches are visited.
    const double comps[] = { -2.0, -0.5, 0.0, 0.75, 3.0 };
    const int n = static_cast<int>(sizeof(comps) / sizeof(comps[0]));

    int checked = 0;
    int mismatches = 0;
    for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) {
            for (int c = 0; c < n; ++c) {
                // vec1 sweeps all components; vec2 is derived from the same
                // pool with rotated indices to keep the loop count modest.
                Vector v1 = asic_make_vector(comps[a], comps[b], comps[c]);
                Vector v2 = asic_make_vector(comps[c], comps[a], comps[b]);

                const bool expected = asic_expected_result(v1, v2);
                const bool actual = angleSignIsCorrect(v1, v2);
                ++checked;
                if (expected != actual) {
                    ++mismatches;
                    std::cerr << "    MISMATCH for vec1 = (" << v1.x << ", "
                              << v1.y << ", " << v1.z << "), vec2 = (" << v2.x
                              << ", " << v2.y << ", " << v2.z
                              << "): expected " << std::boolalpha << expected
                              << " got " << actual << "\n";
                }
                // Non-fatal per-sample check so the whole sweep always runs.
                EXPECT_EQ(expected, actual)
                    << "angleSignIsCorrect disagrees with the reference sign "
                       "formula for vec1 = (" << v1.x << ", " << v1.y << ", "
                    << v1.z << "), vec2 = (" << v2.x << ", " << v2.y << ", "
                    << v2.z << ")";
            }
        }
    }

    std::cerr << "  Sampled " << checked << " vector pairs, " << mismatches
              << " mismatch(es)\n";
    EXPECT_EQ(mismatches, 0) << "No sampled pair should disagree with the reference";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each helper is run inside its own TEST so a failure in
// one scenario still lets the remaining scenarios execute.
// -----------------------------------------------------------------------------
TEST(AngleSignIsCorrect, XZBranchSign) { asic_test_xz_branch_sign(); }
TEST(AngleSignIsCorrect, XZBranchIgnoresY) { asic_test_xz_branch_ignores_y(); }
TEST(AngleSignIsCorrect, XYBranchSign) { asic_test_xy_branch_sign(); }
TEST(AngleSignIsCorrect, XYBranchSelectedWhenOneZIsZero)
{
    asic_test_xy_branch_selected_when_one_z_is_zero();
}
TEST(AngleSignIsCorrect, ZToleranceBoundary) { asic_test_z_tolerance_boundary(); }
TEST(AngleSignIsCorrect, DegeneratePairsReturnFalse) { asic_test_degenerate_pairs_return_false(); }
TEST(AngleSignIsCorrect, MatchesReferenceFormula) { asic_test_matches_reference_formula(); }