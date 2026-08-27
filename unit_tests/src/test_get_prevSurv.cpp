/*! \file test_get_prevSurv.cpp
 *
 * ### Unit test for ../src/reactions/get_prevSurv.cpp
 *
 * Function under test:
 *
 *     double get_prevSurv(const gsl_matrix* survMatrix, double Dtot,
 *                         double deltaT, double r0, double bindRadius)
 *
 * The routine performs a table look-up followed by a linear interpolation on a
 * pre-computed 2D "survival" matrix.  The layout of that matrix is:
 *
 *     row 0 : the radial grid values  r_i
 *     row 1 : the survival value at each grid point,  S(r_i)
 *
 * The grid spacing used by the lookup is fixed by the implementation as
 *
 *     RstepSize = sqrt(Dtot * deltaT) / 50
 *
 * and the column index is
 *
 *     index = floor((r0 - bindRadius) / RstepSize),   clamped so index >= 0
 *
 * The returned value is
 *
 *     (S_i * (r_{i+1} - r0) + S_{i+1} * (r0 - r_i)) / RstepSize
 *
 * which, when the table is built on the same RstepSize grid (r_{i+1} - r_i ==
 * RstepSize), is exactly the linear interpolation of S at r0.
 *
 * The tests below build their own survival matrices (so the expected answers
 * can be derived analytically) and check:
 *   - a constant table returns the constant for any r0 (partition-of-unity),
 *   - a linear table is reproduced exactly (interpolation is exact for lines),
 *   - grid-node values are returned at the nodes,
 *   - r0 below bindRadius clamps the index to 0 and extrapolates from the
 *     first segment (instead of reading out of bounds / negative indices),
 *   - only Dtot*deltaT matters (through RstepSize), not the two separately,
 *   - rows beyond row 1 of the matrix are ignored.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_matrix.h>
#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with gps_ == "get_prevSurv" to avoid collisions with
// other translation units in the combined test binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Allocate and fill a survival matrix in the layout get_prevSurv expects.
 *
 * \param[in] rVals    radial grid points, written to row 0.
 * \param[in] survVals survival values at those grid points, written to row 1.
 * \param[in] nRows    total number of rows to allocate (>= 2).  Any rows past
 *                     row 1 are filled with an obvious sentinel so we can prove
 *                     the function never reads them.
 * \return newly allocated gsl_matrix; caller must gsl_matrix_free() it.
 */
gsl_matrix* gps_make_matrix(const std::vector<double>& rVals,
                            const std::vector<double>& survVals,
                            size_t nRows = 2)
{
    const size_t nCols = rVals.size();
    gsl_matrix* m = gsl_matrix_alloc(nRows, nCols);
    gsl_matrix_set_all(m, -12345.0); // sentinel for the unused rows

    for (size_t i = 0; i < nCols; ++i) {
        gsl_matrix_set(m, 0, i, rVals[i]);   // row 0: r grid
        gsl_matrix_set(m, 1, i, survVals[i]); // row 1: survival values
    }
    return m;
}

/*! \brief Build the r grid the implementation implicitly assumes:
 *         r_i = bindRadius + i * RstepSize.
 */
std::vector<double> gps_make_r_grid(double bindRadius, double rStepSize, size_t nCols)
{
    std::vector<double> rVals(nCols);
    for (size_t i = 0; i < nCols; ++i)
        rVals[i] = bindRadius + static_cast<double>(i) * rStepSize;
    return rVals;
}

/*! \brief The RstepSize the implementation uses internally. */
double gps_step_size(double Dtot, double deltaT) { return std::sqrt(Dtot * deltaT) / 50.0; }

/*! \brief Re-statement of the interpolation formula, used only where the test
 *         needs to say *which* table segment should have been selected.
 */
double gps_segment_value(double val01, double val02, double r01, double r02,
                         double r0, double rStepSize)
{
    return (val01 * (r02 - r0) + val02 * (r0 - r01)) / rStepSize;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: A constant survival table must return that constant for every r0.
//
// This is a "partition of unity" check: the two interpolation weights
// (r02 - r0) and (r0 - r01) always sum to RstepSize on the assumed grid, so a
// flat table has to come back flat regardless of where r0 falls.
// -----------------------------------------------------------------------------
void test_gps_constant_table_returns_constant()
{
    std::cerr << "\n[TEST] test_gps_constant_table_returns_constant\n"
              << "  Source file: src/reactions/get_prevSurv.cpp\n"
              << "  Function:    get_prevSurv\n"
              << "  Scenario:    survival table is a constant 0.75 everywhere.\n"
              << "  Criteria:    interpolated value == 0.75 for many r0 values.\n";

    const double Dtot = 1.0;
    const double deltaT = 1.0;
    const double bindRadius = 1.0;
    const double rStep = gps_step_size(Dtot, deltaT); // sqrt(1*1)/50 = 0.02
    const size_t nCols = 200;

    std::cerr << "  RstepSize = sqrt(Dtot*deltaT)/50 = " << rStep << '\n';

    std::vector<double> rVals = gps_make_r_grid(bindRadius, rStep, nCols);
    std::vector<double> survVals(nCols, 0.75);
    gsl_matrix* survMatrix = gps_make_matrix(rVals, survVals);

    // Sample r0 across many grid cells, including non-node positions.
    for (int k = 0; k < 20; ++k) {
        const double r0 = bindRadius + (k * 3 + 0.37) * rStep;
        const double result = get_prevSurv(survMatrix, Dtot, deltaT, r0, bindRadius);
        EXPECT_NEAR(result, 0.75, 1e-12)
            << "Constant table must return the constant; failed at r0 = " << r0;
    }
    std::cerr << "  Sampled 20 positions, all returned 0.75 within 1e-12.\n";

    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// Test 2: A linear survival table is reproduced exactly.
//
// Linear interpolation reproduces linear data exactly, so for
// S(r) = a + b*(r - bindRadius) the function must return a + b*(r0 - bindRadius)
// for any r0 inside the table.
// -----------------------------------------------------------------------------
void test_gps_linear_table_is_exact()
{
    std::cerr << "\n[TEST] test_gps_linear_table_is_exact\n"
              << "  Source file: src/reactions/get_prevSurv.cpp\n"
              << "  Function:    get_prevSurv\n"
              << "  Scenario:    S(r) = a + b*(r - bindRadius), a=0.2, b=1.5.\n"
              << "  Criteria:    lookup returns the analytic line value at r0.\n";

    const double Dtot = 4.0;
    const double deltaT = 0.25;
    const double bindRadius = 2.0;
    const double rStep = gps_step_size(Dtot, deltaT); // sqrt(1)/50 = 0.02
    const size_t nCols = 300;

    const double a = 0.2;
    const double b = 1.5;

    std::vector<double> rVals = gps_make_r_grid(bindRadius, rStep, nCols);
    std::vector<double> survVals(nCols);
    for (size_t i = 0; i < nCols; ++i)
        survVals[i] = a + b * (rVals[i] - bindRadius);

    gsl_matrix* survMatrix = gps_make_matrix(rVals, survVals);

    // Probe several offsets: exactly on nodes, a quarter way, mid-cell, etc.
    const double offsets[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 5.5, 17.3, 100.9 };
    for (double off : offsets) {
        const double r0 = bindRadius + off * rStep;
        const double expected = a + b * (r0 - bindRadius);
        const double result = get_prevSurv(survMatrix, Dtot, deltaT, r0, bindRadius);
        std::cerr << "    r0 - bindRadius = " << (r0 - bindRadius)
                  << "  expected " << expected << "  got " << result << '\n';
        EXPECT_NEAR(result, expected, 1e-10)
            << "Linear table should be interpolated exactly at r0 = " << r0;
    }

    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// Test 3: Landing on a grid node returns that node's tabulated value.
//
// Note the interpolation is continuous, so it does not matter whether floating
// point rounding puts the node into segment [k-1,k] or [k,k+1] -- both give the
// same answer, S(r_k).
// -----------------------------------------------------------------------------
void test_gps_returns_node_value_at_grid_points()
{
    std::cerr << "\n[TEST] test_gps_returns_node_value_at_grid_points\n"
              << "  Source file: src/reactions/get_prevSurv.cpp\n"
              << "  Function:    get_prevSurv\n"
              << "  Scenario:    r0 sits exactly on tabulated grid points of a\n"
              << "               non-linear (exponential-like) table.\n"
              << "  Criteria:    the tabulated value at that node is returned.\n";

    const double Dtot = 1.0;
    const double deltaT = 1.0;
    const double bindRadius = 1.0;
    const double rStep = gps_step_size(Dtot, deltaT); // 0.02
    const size_t nCols = 150;

    std::vector<double> rVals = gps_make_r_grid(bindRadius, rStep, nCols);
    std::vector<double> survVals(nCols);
    for (size_t i = 0; i < nCols; ++i)
        survVals[i] = std::exp(-(rVals[i] - bindRadius)); // smooth, non-linear

    gsl_matrix* survMatrix = gps_make_matrix(rVals, survVals);

    // Check several nodes.  Tolerance is loose enough to absorb the round-off in
    // reconstructing r0 = bindRadius + k*rStep, but tight enough to fail if the
    // wrong table entry were used (neighbouring entries differ by ~2e-2).
    for (size_t k : { size_t { 0 }, size_t { 1 }, size_t { 7 }, size_t { 42 }, size_t { 100 } }) {
        const double r0 = rVals[k];
        const double result = get_prevSurv(survMatrix, Dtot, deltaT, r0, bindRadius);
        std::cerr << "    node " << k << " (r0 = " << r0 << ") expected "
                  << survVals[k] << "  got " << result << '\n';
        EXPECT_NEAR(result, survVals[k], 1e-9)
            << "Value at grid node " << k << " should be the tabulated value";
    }

    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// Test 4: r0 smaller than bindRadius clamps the index to 0.
//
// floor((r0 - bindRadius)/RstepSize) is negative in that case; the code clamps
// it to 0 so the first table segment [r_0, r_1] is used and the formula
// extrapolates.  We verify the returned number equals the extrapolation from
// that first segment (and NOT, e.g., simply the value at r_0).
// -----------------------------------------------------------------------------
void test_gps_clamps_negative_index_to_zero()
{
    std::cerr << "\n[TEST] test_gps_clamps_negative_index_to_zero\n"
              << "  Source file: src/reactions/get_prevSurv.cpp\n"
              << "  Function:    get_prevSurv\n"
              << "  Scenario:    r0 < bindRadius so floor() would give a negative\n"
              << "               column index.\n"
              << "  Criteria:    index is clamped to 0 and the first table segment\n"
              << "               is used for the (extrapolating) interpolation.\n";

    const double Dtot = 1.0;
    const double deltaT = 1.0;
    const double bindRadius = 1.0;
    const double rStep = gps_step_size(Dtot, deltaT); // 0.02
    const size_t nCols = 50;

    std::vector<double> rVals = gps_make_r_grid(bindRadius, rStep, nCols);
    std::vector<double> survVals(nCols);
    for (size_t i = 0; i < nCols; ++i)
        survVals[i] = std::exp(-(rVals[i] - bindRadius)); // non-linear on purpose

    gsl_matrix* survMatrix = gps_make_matrix(rVals, survVals);

    // r0 half a step BELOW the binding radius -> raw index would be -1.
    const double r0 = bindRadius - 0.5 * rStep;

    // Expected: the clamped segment is columns 0 and 1.
    const double expected = gps_segment_value(survVals[0], survVals[1],
                                              rVals[0], rVals[1], r0, rStep);

    const double result = get_prevSurv(survMatrix, Dtot, deltaT, r0, bindRadius);
    std::cerr << "    r0 = " << r0 << " (bindRadius = " << bindRadius << ")\n"
              << "    expected (segment 0-1 extrapolation) = " << expected
              << "  got " << result << '\n';

    EXPECT_NEAR(result, expected, 1e-12)
        << "Negative index must clamp to 0 and use the first table segment";

    // Sanity: for this decaying table the extrapolated value must exceed S(r_0),
    // which is what distinguishes "clamped and extrapolated" from "returned S(r_0)".
    EXPECT_GT(result, survVals[0])
        << "Extrapolation below bindRadius should exceed S at the first node";

    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// Test 5: Only the product Dtot*deltaT matters, because RstepSize depends on
// sqrt(Dtot*deltaT).  Two different (Dtot, deltaT) pairs with the same product
// must give bit-comparable answers on the same table.
// -----------------------------------------------------------------------------
void test_gps_depends_only_on_Dtot_times_deltaT()
{
    std::cerr << "\n[TEST] test_gps_depends_only_on_Dtot_times_deltaT\n"
              << "  Source file: src/reactions/get_prevSurv.cpp\n"
              << "  Function:    get_prevSurv\n"
              << "  Scenario:    (Dtot=4, deltaT=1) vs (Dtot=1, deltaT=4) --\n"
              << "               same product, hence same RstepSize.\n"
              << "  Criteria:    both calls return the same value; a pair with a\n"
              << "               different product returns something different.\n";

    const double bindRadius = 0.5;

    const double DtotA = 4.0, deltaTA = 1.0; // product 4 -> RstepSize = 2/50
    const double DtotB = 1.0, deltaTB = 4.0; // product 4 -> same RstepSize
    const double DtotC = 1.0, deltaTC = 1.0; // product 1 -> different RstepSize

    const double rStepAB = gps_step_size(DtotA, deltaTA);
    std::cerr << "  RstepSize(A) = " << rStepAB
              << ", RstepSize(B) = " << gps_step_size(DtotB, deltaTB)
              << ", RstepSize(C) = " << gps_step_size(DtotC, deltaTC) << '\n';

    // Build the table on the A/B grid so A and B are "correct" look-ups.
    const size_t nCols = 200;
    std::vector<double> rVals = gps_make_r_grid(bindRadius, rStepAB, nCols);
    std::vector<double> survVals(nCols);
    for (size_t i = 0; i < nCols; ++i)
        survVals[i] = 1.0 / (1.0 + (rVals[i] - bindRadius)); // monotone, non-linear

    gsl_matrix* survMatrix = gps_make_matrix(rVals, survVals);

    const double r0 = bindRadius + 3.4 * rStepAB;

    const double resA = get_prevSurv(survMatrix, DtotA, deltaTA, r0, bindRadius);
    const double resB = get_prevSurv(survMatrix, DtotB, deltaTB, r0, bindRadius);
    const double resC = get_prevSurv(survMatrix, DtotC, deltaTC, r0, bindRadius);

    std::cerr << "    result(A) = " << resA << ", result(B) = " << resB
              << ", result(C) = " << resC << '\n';

    EXPECT_DOUBLE_EQ(resA, resB)
        << "Equal Dtot*deltaT products must give identical results";
    EXPECT_NE(resA, resC)
        << "A different Dtot*deltaT product changes RstepSize and the result";

    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// Test 6: Only rows 0 and 1 of the matrix are consulted.
//
// We build a 5-row matrix whose extra rows are filled with an obviously bogus
// sentinel value and confirm the answer is unchanged relative to a 2-row matrix
// carrying the same data.
// -----------------------------------------------------------------------------
void test_gps_ignores_rows_beyond_first_two()
{
    std::cerr << "\n[TEST] test_gps_ignores_rows_beyond_first_two\n"
              << "  Source file: src/reactions/get_prevSurv.cpp\n"
              << "  Function:    get_prevSurv\n"
              << "  Scenario:    identical data stored in a 2-row and a 5-row\n"
              << "               matrix (extra rows filled with -12345).\n"
              << "  Criteria:    both matrices produce exactly the same result.\n";

    const double Dtot = 2.0;
    const double deltaT = 0.5;
    const double bindRadius = 1.25;
    const double rStep = gps_step_size(Dtot, deltaT); // sqrt(1)/50 = 0.02
    const size_t nCols = 120;

    std::vector<double> rVals = gps_make_r_grid(bindRadius, rStep, nCols);
    std::vector<double> survVals(nCols);
    for (size_t i = 0; i < nCols; ++i)
        survVals[i] = 0.9 - 0.3 * (rVals[i] - bindRadius);

    gsl_matrix* twoRow = gps_make_matrix(rVals, survVals, 2);
    gsl_matrix* fiveRow = gps_make_matrix(rVals, survVals, 5);

    const double r0 = bindRadius + 11.75 * rStep;
    const double resTwo = get_prevSurv(twoRow, Dtot, deltaT, r0, bindRadius);
    const double resFive = get_prevSurv(fiveRow, Dtot, deltaT, r0, bindRadius);

    std::cerr << "    2-row result = " << resTwo
              << ", 5-row result = " << resFive << '\n';

    EXPECT_DOUBLE_EQ(resTwo, resFive)
        << "Rows beyond row 1 must not influence the look-up";

    // And the value should still be finite / physically sensible for this table.
    EXPECT_TRUE(std::isfinite(resTwo)) << "Result must be a finite number";

    gsl_matrix_free(twoRow);
    gsl_matrix_free(fiveRow);
}

// -----------------------------------------------------------------------------
// Test 7: Monotone table -> monotone look-up.
//
// A decreasing survival table must produce a non-increasing sequence of
// interpolated values as r0 increases.  This guards against index/weight sign
// mistakes that would otherwise still pass the "constant table" check.
// -----------------------------------------------------------------------------
void test_gps_monotone_table_gives_monotone_output()
{
    std::cerr << "\n[TEST] test_gps_monotone_table_gives_monotone_output\n"
              << "  Source file: src/reactions/get_prevSurv.cpp\n"
              << "  Function:    get_prevSurv\n"
              << "  Scenario:    strictly decreasing survival table, r0 swept\n"
              << "               upward across many cells.\n"
              << "  Criteria:    successive look-ups are non-increasing.\n";

    const double Dtot = 1.0;
    const double deltaT = 1.0;
    const double bindRadius = 1.0;
    const double rStep = gps_step_size(Dtot, deltaT);
    const size_t nCols = 400;

    std::vector<double> rVals = gps_make_r_grid(bindRadius, rStep, nCols);
    std::vector<double> survVals(nCols);
    for (size_t i = 0; i < nCols; ++i)
        survVals[i] = std::exp(-2.0 * (rVals[i] - bindRadius)); // strictly decreasing

    gsl_matrix* survMatrix = gps_make_matrix(rVals, survVals);

    double previous = get_prevSurv(survMatrix, Dtot, deltaT, bindRadius, bindRadius);
    int nChecked = 0;
    for (int k = 1; k < 60; ++k) {
        const double r0 = bindRadius + (k * 0.5) * rStep; // half-cell steps
        const double current = get_prevSurv(survMatrix, Dtot, deltaT, r0, bindRadius);
        EXPECT_LE(current, previous + 1e-12)
            << "Look-up should not increase for a decreasing table (step " << k << ')';
        previous = current;
        ++nChecked;
    }
    std::cerr << "  Verified monotonicity across " << nChecked << " half-cell steps;"
              << " final value = " << previous << '\n';

    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so that
// a failure in one case does not stop the remaining cases from running.
// -----------------------------------------------------------------------------
TEST(GetPrevSurvTest, ConstantTableReturnsConstant) { test_gps_constant_table_returns_constant(); }
TEST(GetPrevSurvTest, LinearTableIsExact) { test_gps_linear_table_is_exact(); }
TEST(GetPrevSurvTest, ReturnsNodeValueAtGridPoints) { test_gps_returns_node_value_at_grid_points(); }
TEST(GetPrevSurvTest, ClampsNegativeIndexToZero) { test_gps_clamps_negative_index_to_zero(); }
TEST(GetPrevSurvTest, DependsOnlyOnDtotTimesDeltaT) { test_gps_depends_only_on_Dtot_times_deltaT(); }
TEST(GetPrevSurvTest, IgnoresRowsBeyondFirstTwo) { test_gps_ignores_rows_beyond_first_two(); }
TEST(GetPrevSurvTest, MonotoneTableGivesMonotoneOutput) { test_gps_monotone_table_gives_monotone_output(); }