/*! \file test_get_prevNorm.cpp
 *
 * ### Unit test for src/reactions/get_prevNorm.cpp
 *
 * Function under test:
 *
 *     double get_prevNorm(gsl_matrix* normMatrix, double RStepSize,
 *                         double r0, double bindRadius)
 *
 * Behaviour, as read directly from the implementation:
 *
 *   1. The lookup index is `floor((r0 - bindRadius) / RStepSize)`, clamped up
 *      to 0 when it would be negative (i.e. when r0 < bindRadius).
 *   2. Row 0 of `normMatrix` holds the separation (r) values of the table,
 *      row 1 holds the corresponding normalization values.
 *   3. The return value is the linear interpolation
 *
 *          (val[i] * (r[i+1] - r0) + val[i+1] * (r0 - r[i])) / RStepSize
 *
 *      Note that the divisor is the *caller supplied* RStepSize, not
 *      (r[i+1] - r[i]) taken from the table. The function therefore trusts the
 *      caller that row 0 of the matrix is a uniform grid of spacing RStepSize
 *      starting at bindRadius; the tests below exercise both the "consistent"
 *      case and one deliberately inconsistent case to pin down that behaviour.
 *
 * IMPORTANT: the function reads column `index + 1`. If r0 sits at (or beyond)
 * the very last column of the table, `index + 1` is out of range and
 * gsl_matrix_get() invokes the GSL error handler, which aborts the process.
 * Every test below therefore keeps r0 strictly below the final grid point.
 * That failure mode is intentionally *not* exercised, because aborting would
 * kill the whole gtest binary rather than fail a single case.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_matrix.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <iostream>
#include <vector>

namespace {

/*! \brief Allocate a 2 x N gsl_matrix table for get_prevNorm().
 *
 * Row 0 receives the separation values, row 1 receives the normalization
 * values. The caller owns the returned matrix and must gsl_matrix_free() it.
 *
 * \param[in] rVals    separation grid (goes into row 0)
 * \param[in] normVals tabulated norm values (goes into row 1)
 */
gsl_matrix* gpn_build_table(const std::vector<double>& rVals,
    const std::vector<double>& normVals)
{
    // Both rows must be the same length or the table is malformed.
    EXPECT_EQ(rVals.size(), normVals.size())
        << "Test helper misuse: r and norm vectors must have equal length";

    gsl_matrix* table = gsl_matrix_alloc(2, rVals.size());
    for (std::size_t i = 0; i < rVals.size(); ++i) {
        gsl_matrix_set(table, 0, i, rVals[i]);
        gsl_matrix_set(table, 1, i, normVals[i]);
    }
    return table;
}

/*! \brief Reference implementation of the linear function used for the
 *         "consistent table" tests: f(r) = 2r + 1.
 *
 * Interpolating a linear function on a uniform grid must reproduce the
 * function exactly, which gives us an independent expected value.
 */
double gpn_linear_ref(double r) { return 2.0 * r + 1.0; }

} // namespace

// -----------------------------------------------------------------------------
// Test 1: r0 sitting exactly on a grid point returns the tabulated value.
// -----------------------------------------------------------------------------
void test_gpn_exact_grid_points()
{
    std::cerr << "\n[TEST] test_gpn_exact_grid_points\n"
              << "  Source file:   src/reactions/get_prevNorm.cpp\n"
              << "  Function:      get_prevNorm()\n"
              << "  Scenario:      r0 lands exactly on a tabulated grid point.\n"
              << "  Pass criteria: the interpolation collapses to the stored\n"
              << "                 row-1 value for that column.\n";

    // Uniform grid: bindRadius = 1.0, step = 0.25 (a power of two, so every
    // grid coordinate and every division below is exact in binary floating
    // point -- this lets us assert with EXPECT_DOUBLE_EQ).
    const double bindRadius = 1.0;
    const double RStepSize = 0.25;
    const int nCols = 9; // r = 1.00 .. 3.00

    std::vector<double> rVals;
    std::vector<double> normVals;
    for (int i = 0; i < nCols; ++i) {
        const double r = bindRadius + i * RStepSize;
        rVals.push_back(r);
        normVals.push_back(gpn_linear_ref(r));
    }

    gsl_matrix* table = gpn_build_table(rVals, normVals);

    // Probe several interior grid points. The final column (index nCols-1) is
    // deliberately skipped: index+1 would be out of range there.
    for (int i = 0; i < nCols - 1; ++i) {
        const double r0 = rVals[i];
        const double got = get_prevNorm(table, RStepSize, r0, bindRadius);
        std::cerr << "    r0 = " << r0 << " -> get_prevNorm = " << got
                  << " (tabulated value " << normVals[i] << ")\n";
        EXPECT_DOUBLE_EQ(got, normVals[i])
            << "At an exact grid point the result must equal the stored value";
    }

    gsl_matrix_free(table);
}

// -----------------------------------------------------------------------------
// Test 2: interpolation between grid points reproduces a linear table exactly.
// -----------------------------------------------------------------------------
void test_gpn_linear_interpolation()
{
    std::cerr << "\n[TEST] test_gpn_linear_interpolation\n"
              << "  Source file:   src/reactions/get_prevNorm.cpp\n"
              << "  Function:      get_prevNorm()\n"
              << "  Scenario:      r0 falls between two grid points of a table\n"
              << "                 that tabulates the linear function 2r + 1.\n"
              << "  Pass criteria: linear interpolation of a linear function\n"
              << "                 reproduces the function value at r0.\n";

    const double bindRadius = 1.0;
    const double RStepSize = 0.25;
    const int nCols = 9;

    std::vector<double> rVals;
    std::vector<double> normVals;
    for (int i = 0; i < nCols; ++i) {
        const double r = bindRadius + i * RStepSize;
        rVals.push_back(r);
        normVals.push_back(gpn_linear_ref(r));
    }

    gsl_matrix* table = gpn_build_table(rVals, normVals);

    // Midpoints and off-centre points, all strictly inside the table so that
    // index + 1 stays in range.
    const std::vector<double> probes { 1.125, 1.625, 2.0625, 2.4, 2.74 };

    for (double r0 : probes) {
        const double got = get_prevNorm(table, RStepSize, r0, bindRadius);
        const double want = gpn_linear_ref(r0);
        std::cerr << "    r0 = " << r0 << " -> get_prevNorm = " << got
                  << ", expected 2*r0+1 = " << want << '\n';
        EXPECT_NEAR(got, want, 1e-12)
            << "Interpolating a linear table must reproduce the linear function";
    }

    gsl_matrix_free(table);
}

// -----------------------------------------------------------------------------
// Test 3: r0 below bindRadius clamps the index to 0 (linear extrapolation).
// -----------------------------------------------------------------------------
void test_gpn_below_bind_radius_clamps_index()
{
    std::cerr << "\n[TEST] test_gpn_below_bind_radius_clamps_index\n"
              << "  Source file:   src/reactions/get_prevNorm.cpp\n"
              << "  Function:      get_prevNorm()\n"
              << "  Scenario:      r0 < bindRadius, so floor(...) is negative\n"
              << "                 and the code clamps the index to 0.\n"
              << "  Pass criteria: the value is computed from columns 0 and 1,\n"
              << "                 i.e. a straight-line extrapolation below the\n"
              << "                 first grid point (no clamping of the value).\n";

    const double bindRadius = 1.0;
    const double RStepSize = 0.25;
    const int nCols = 9;

    std::vector<double> rVals;
    std::vector<double> normVals;
    for (int i = 0; i < nCols; ++i) {
        const double r = bindRadius + i * RStepSize;
        rVals.push_back(r);
        normVals.push_back(gpn_linear_ref(r));
    }

    gsl_matrix* table = gpn_build_table(rVals, normVals);

    // r0 = 0.875 gives (0.875 - 1.0)/0.25 = -0.5, floor = -1, clamped to 0.
    // Columns 0 and 1 are then used, so for a linear table the result is the
    // straight-line extrapolation, which equals 2*r0 + 1 = 2.75.
    const double r0Below = 0.875;
    const double gotBelow = get_prevNorm(table, RStepSize, r0Below, bindRadius);
    std::cerr << "    r0 = " << r0Below << " (below bindRadius " << bindRadius
              << ") -> get_prevNorm = " << gotBelow
              << ", expected extrapolation " << gpn_linear_ref(r0Below) << '\n';
    EXPECT_NEAR(gotBelow, gpn_linear_ref(r0Below), 1e-12)
        << "Negative indices are clamped to 0 and the value extrapolates linearly";

    // r0 exactly at bindRadius: index = 0 exactly, result is the first value.
    const double gotAt = get_prevNorm(table, RStepSize, bindRadius, bindRadius);
    std::cerr << "    r0 = bindRadius = " << bindRadius
              << " -> get_prevNorm = " << gotAt
              << ", expected first tabulated value " << normVals[0] << '\n';
    EXPECT_DOUBLE_EQ(gotAt, normVals[0])
        << "r0 == bindRadius should return the very first tabulated value";

    // A deeply negative r0 also clamps to index 0 (no crash, no clamping of the
    // returned value -- the extrapolation simply continues the first segment).
    const double r0Far = 0.0;
    const double gotFar = get_prevNorm(table, RStepSize, r0Far, bindRadius);
    std::cerr << "    r0 = " << r0Far << " -> get_prevNorm = " << gotFar
              << ", expected extrapolation " << gpn_linear_ref(r0Far) << '\n';
    EXPECT_NEAR(gotFar, gpn_linear_ref(r0Far), 1e-12)
        << "Any r0 below bindRadius uses columns 0 and 1 of the table";

    gsl_matrix_free(table);
}

// -----------------------------------------------------------------------------
// Test 4: with a non-linear table, the result matches the hand-computed
//         two-point interpolation formula (not the true function value).
// -----------------------------------------------------------------------------
void test_gpn_nonlinear_matches_formula()
{
    std::cerr << "\n[TEST] test_gpn_nonlinear_matches_formula\n"
              << "  Source file:   src/reactions/get_prevNorm.cpp\n"
              << "  Function:      get_prevNorm()\n"
              << "  Scenario:      table holds arbitrary (non-linear) values.\n"
              << "  Pass criteria: the result equals the hand evaluated formula\n"
              << "                 (v_i*(r_{i+1}-r0) + v_{i+1}*(r0-r_i))/step,\n"
              << "                 i.e. plain linear interpolation, and does NOT\n"
              << "                 attempt any higher-order reconstruction.\n";

    const double bindRadius = 2.0;
    const double RStepSize = 0.5;

    // Uniform r grid, arbitrary norm values.
    const std::vector<double> rVals { 2.0, 2.5, 3.0, 3.5, 4.0 };
    const std::vector<double> normVals { 10.0, 4.0, 7.0, -1.0, 0.5 };

    gsl_matrix* table = gpn_build_table(rVals, normVals);

    // r0 = 3.2 -> (3.2 - 2.0)/0.5 = 2.4 -> floor = 2, so columns 2 and 3.
    const double r0 = 3.2;
    const std::size_t idx = 2;
    const double expected = (normVals[idx] * (rVals[idx + 1] - r0)
                                + normVals[idx + 1] * (r0 - rVals[idx]))
        / RStepSize;

    const double got = get_prevNorm(table, RStepSize, r0, bindRadius);
    std::cerr << "    r0 = " << r0 << " -> index " << idx
              << ", get_prevNorm = " << got
              << ", hand-computed expected = " << expected << '\n';
    EXPECT_NEAR(got, expected, 1e-12)
        << "Result must equal the documented two-point interpolation formula";

    // Sanity check on the bracketing property: because the interpolation weight
    // is between 0 and 1 here, the answer must lie between the two neighbours.
    const double lo = std::min(normVals[idx], normVals[idx + 1]);
    const double hi = std::max(normVals[idx], normVals[idx + 1]);
    std::cerr << "    checking result is bracketed by neighbours [" << lo << ", "
              << hi << "]\n";
    EXPECT_GE(got, lo - 1e-12) << "Interior interpolation cannot undershoot both neighbours";
    EXPECT_LE(got, hi + 1e-12) << "Interior interpolation cannot overshoot both neighbours";

    gsl_matrix_free(table);
}

// -----------------------------------------------------------------------------
// Test 5: the r values used in the formula come from row 0 of the matrix, while
//         the index and the divisor come from bindRadius / RStepSize.
// -----------------------------------------------------------------------------
void test_gpn_uses_row0_positions_and_step_divisor()
{
    std::cerr << "\n[TEST] test_gpn_uses_row0_positions_and_step_divisor\n"
              << "  Source file:   src/reactions/get_prevNorm.cpp\n"
              << "  Function:      get_prevNorm()\n"
              << "  Scenario:      the table's row-0 grid deliberately does NOT\n"
              << "                 match bindRadius/RStepSize.\n"
              << "  Pass criteria: the function still uses row 0 for r_i/r_{i+1}\n"
              << "                 and RStepSize as the divisor, producing the\n"
              << "                 (physically meaningless but well defined)\n"
              << "                 value we compute by hand. This documents the\n"
              << "                 implicit precondition on the caller.\n";

    // Index selection uses bindRadius = 1.0 and step = 0.25 ...
    const double bindRadius = 1.0;
    const double RStepSize = 0.25;

    // ... but row 0 is a completely different (coarse) grid.
    const std::vector<double> rVals { 0.0, 10.0, 20.0, 30.0 };
    const std::vector<double> normVals { 1.0, 2.0, 3.0, 4.0 };

    gsl_matrix* table = gpn_build_table(rVals, normVals);

    // r0 = 1.6 -> (1.6 - 1.0)/0.25 = 2.4 -> floor = 2, so columns 2 and 3.
    const double r0 = 1.6;
    const std::size_t idx = 2;
    const double expected = (normVals[idx] * (rVals[idx + 1] - r0)
                                + normVals[idx + 1] * (r0 - rVals[idx]))
        / RStepSize;

    const double got = get_prevNorm(table, RStepSize, r0, bindRadius);
    std::cerr << "    r0 = " << r0 << " -> index " << idx
              << " (from bindRadius/RStepSize), r_i = " << rVals[idx]
              << ", r_i+1 = " << rVals[idx + 1] << '\n';
    std::cerr << "    get_prevNorm = " << got << ", hand-computed expected = "
              << expected << '\n';
    EXPECT_NEAR(got, expected, 1e-9)
        << "Positions come from row 0; the divisor is the supplied RStepSize";

    gsl_matrix_free(table);
}

// -----------------------------------------------------------------------------
// Test 6: monotone table -> monotone interpolated output as r0 increases.
// -----------------------------------------------------------------------------
void test_gpn_monotonic_table_gives_monotonic_output()
{
    std::cerr << "\n[TEST] test_gpn_monotonic_table_gives_monotonic_output\n"
              << "  Source file:   src/reactions/get_prevNorm.cpp\n"
              << "  Function:      get_prevNorm()\n"
              << "  Scenario:      table values decrease monotonically with r.\n"
              << "  Pass criteria: successive calls with increasing r0 return\n"
              << "                 non-increasing values (piecewise-linear\n"
              << "                 interpolation preserves monotonicity).\n";

    const double bindRadius = 1.0;
    const double RStepSize = 0.25;
    const int nCols = 9;

    std::vector<double> rVals;
    std::vector<double> normVals;
    for (int i = 0; i < nCols; ++i) {
        const double r = bindRadius + i * RStepSize;
        rVals.push_back(r);
        // Strictly decreasing sequence: 1.0, 0.9, 0.8, ...
        normVals.push_back(1.0 - 0.1 * i);
    }

    gsl_matrix* table = gpn_build_table(rVals, normVals);

    // Walk r0 across the table, staying strictly below the last grid point.
    double previous = get_prevNorm(table, RStepSize, bindRadius, bindRadius);
    std::cerr << "    starting value at r0 = " << bindRadius << " is " << previous << '\n';

    for (double r0 = bindRadius + 0.05; r0 < rVals[nCols - 1] - 1e-6; r0 += 0.05) {
        const double got = get_prevNorm(table, RStepSize, r0, bindRadius);
        EXPECT_LE(got, previous + 1e-12)
            << "Interpolated norm should not increase for a decreasing table (r0 = "
            << r0 << ')';
        previous = got;
    }
    std::cerr << "    final value is " << previous
              << " (monotonic decrease confirmed across the table)\n";

    gsl_matrix_free(table);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* routine runs inside its own TEST so
// that a failure in one does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(GetPrevNorm, ExactGridPoints) { test_gpn_exact_grid_points(); }
TEST(GetPrevNorm, LinearInterpolation) { test_gpn_linear_interpolation(); }
TEST(GetPrevNorm, BelowBindRadiusClampsIndex) { test_gpn_below_bind_radius_clamps_index(); }
TEST(GetPrevNorm, NonLinearMatchesFormula) { test_gpn_nonlinear_matches_formula(); }
TEST(GetPrevNorm, UsesRow0PositionsAndStepDivisor) { test_gpn_uses_row0_positions_and_step_divisor(); }
TEST(GetPrevNorm, MonotonicTableGivesMonotonicOutput) { test_gpn_monotonic_table_gives_monotonic_output(); }