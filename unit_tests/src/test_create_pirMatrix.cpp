/*! \file test_create_pirMatrix.cpp
 *
 * ### Unit test for src/reactions/create_pirMatrix.cpp
 *
 * The single function under test is
 *
 *     void create_pirMatrix(gsl_matrix*& pirMatrix, double bindRadius, double Dtot,
 *                           double kr, double comRMax, double RStepSize,
 *                           const Parameters& params);
 *
 * It walks a 2D grid of separations,
 *
 *     r  (column index, "itr2") from bindRadius to comRMax + RStepSize
 *     r0 (row    index, "itr1") from bindRadius to comRMax + RStepSize
 *
 * and for every (r0, r) pair numerically integrates pir_function() and stores the
 * result with gsl_matrix_set(pirMatrix, itr1, itr2, result).
 *
 * Since the physical value of the integral is not something we can reproduce
 * analytically here, the tests verify the *contract* of the routine:
 *
 *   1. It writes exactly the expected rectangular block of the supplied matrix
 *      (dimension determined by bindRadius / comRMax / RStepSize) and touches
 *      nothing outside of it.  We detect this by pre-filling the whole matrix
 *      with a sentinel value.
 *   2. Every value it writes is a finite (non-NaN, non-inf) number.
 *   3. The routine is deterministic - two identical calls produce bit-identical
 *      matrices.
 *   4. Growing comRMax grows the filled block.
 *   5. The routine tolerates different rate / diffusion inputs and still returns
 *      only finite numbers.
 *
 * All grid sizes are kept deliberately small (large RStepSize) because each grid
 * point requires a numerical quadrature, which is expensive.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace {

//! Sentinel written into every matrix element before the call under test.
//! Any element still holding this value afterwards was *not* touched.
constexpr double kCpmSentinel = -12345.0;

/*! \brief Reproduces, exactly, the loop bound arithmetic used inside
 *         create_pirMatrix() so we know how many rows/columns will be written.
 *
 * The accumulation order (start at bindRadius, repeatedly += RStepSize) is
 * duplicated on purpose so floating point round-off matches the source.
 */
size_t cpm_expected_dim(double bindRadius, double comRMax, double RStepSize)
{
    size_t count { 0 };
    double index { bindRadius };
    while (index <= comRMax + RStepSize) {
        ++count;
        index += RStepSize;
    }
    return count;
}

/*! \brief Builds the minimal Parameters object the routine reads (timeStep only). */
Parameters cpm_make_params(double timeStep)
{
    Parameters params;
    params.timeStep = timeStep; // microseconds; only field used by create_pirMatrix
    return params;
}

/*! \brief Allocates a square gsl_matrix and fills it with the sentinel value. */
gsl_matrix* cpm_make_sentinel_matrix(size_t dim)
{
    gsl_matrix* m = gsl_matrix_alloc(dim, dim);
    gsl_matrix_set_all(m, kCpmSentinel);
    return m;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the routine fills exactly the expected block and leaves the padding
//         rows/columns untouched.
// -----------------------------------------------------------------------------
void test_cpm_fills_expected_block()
{
    std::cerr << "\n[TEST] test_cpm_fills_expected_block\n"
              << "  Source file: src/reactions/create_pirMatrix.cpp\n"
              << "  Function:    create_pirMatrix()\n"
              << "  Scenario:    pre-fill the matrix with a sentinel, run the routine,\n"
              << "               then confirm which elements were overwritten.\n"
              << "  Pass criteria: every element of the expected dim x dim block holds\n"
              << "               a finite number (not the sentinel), and the extra\n"
              << "               padding row/column still holds the sentinel.\n";

    // Small, cheap grid: bindRadius = 1, comRMax = 3, RStepSize = 0.5
    const double bindRadius { 1.0 };
    const double Dtot { 1.0 };
    const double kr { 10.0 };
    const double comRMax { 3.0 };
    const double RStepSize { 0.5 };
    const Parameters params { cpm_make_params(0.1) };

    // Expected number of grid points along each axis.
    const size_t dim { cpm_expected_dim(bindRadius, comRMax, RStepSize) };
    // Two extra rows/columns of padding so we can prove nothing spills over.
    const size_t alloc { dim + 2 };

    std::cerr << "  Expected filled dimension = " << dim
              << ", allocated matrix = " << alloc << " x " << alloc << '\n';
    ASSERT_GT(dim, 1u) << "Test setup error: grid must have more than one point";

    gsl_matrix* pirMatrix = cpm_make_sentinel_matrix(alloc);

    std::cerr << "  Calling create_pirMatrix(bindRadius=" << bindRadius
              << ", Dtot=" << Dtot << ", kr=" << kr << ", comRMax=" << comRMax
              << ", RStepSize=" << RStepSize << ", timeStep=" << params.timeStep << ")\n";
    create_pirMatrix(pirMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);

    // ---- every element in the expected block must have been written -----------
    size_t untouchedInside { 0 };
    size_t nonFiniteInside { 0 };
    for (size_t row = 0; row < dim; ++row) {
        for (size_t col = 0; col < dim; ++col) {
            const double val = gsl_matrix_get(pirMatrix, row, col);
            if (val == kCpmSentinel)
                ++untouchedInside;
            if (!std::isfinite(val))
                ++nonFiniteInside;
        }
    }
    EXPECT_EQ(untouchedInside, 0u)
        << "Every (r0, r) grid point inside the expected block should be written";
    EXPECT_EQ(nonFiniteInside, 0u)
        << "Every written pir value must be a finite number (no NaN/inf)";

    // ---- the padding must be untouched ---------------------------------------
    size_t touchedPadding { 0 };
    for (size_t idx = 0; idx < alloc; ++idx) {
        // Last allocated row and last allocated column are pure padding.
        if (gsl_matrix_get(pirMatrix, alloc - 1, idx) != kCpmSentinel)
            ++touchedPadding;
        if (gsl_matrix_get(pirMatrix, idx, alloc - 1) != kCpmSentinel)
            ++touchedPadding;
    }
    EXPECT_EQ(touchedPadding, 0u)
        << "create_pirMatrix must not write outside the dim x dim grid block";

    std::cerr << "  Sample values: pir[0][0] = " << gsl_matrix_get(pirMatrix, 0, 0)
              << ", pir[" << dim - 1 << "][" << dim - 1 << "] = "
              << gsl_matrix_get(pirMatrix, dim - 1, dim - 1) << '\n';

    gsl_matrix_free(pirMatrix);
}

// -----------------------------------------------------------------------------
// Test 2: the written values are finite and not all identical (the table really
//         does vary with r and r0).
// -----------------------------------------------------------------------------
void test_cpm_values_vary_across_grid()
{
    std::cerr << "\n[TEST] test_cpm_values_vary_across_grid\n"
              << "  Source file: src/reactions/create_pirMatrix.cpp\n"
              << "  Function:    create_pirMatrix()\n"
              << "  Scenario:    inspect the populated table for variation.\n"
              << "  Pass criteria: all values finite AND at least two elements differ,\n"
              << "               proving the (r, r0) dependence is actually evaluated.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 1.0 };
    const double kr { 10.0 };
    const double comRMax { 2.5 };
    const double RStepSize { 0.5 };
    const Parameters params { cpm_make_params(0.1) };

    const size_t dim { cpm_expected_dim(bindRadius, comRMax, RStepSize) };
    gsl_matrix* pirMatrix = cpm_make_sentinel_matrix(dim + 1);

    std::cerr << "  Grid dimension = " << dim << "; calling create_pirMatrix...\n";
    create_pirMatrix(pirMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);

    // Collect statistics over the filled block.
    double minVal { gsl_matrix_get(pirMatrix, 0, 0) };
    double maxVal { minVal };
    bool allFinite { true };
    for (size_t row = 0; row < dim; ++row) {
        for (size_t col = 0; col < dim; ++col) {
            const double val = gsl_matrix_get(pirMatrix, row, col);
            allFinite = allFinite && std::isfinite(val);
            if (val < minVal)
                minVal = val;
            if (val > maxVal)
                maxVal = val;
        }
    }

    std::cerr << "  min = " << minVal << ", max = " << maxVal << '\n';
    EXPECT_TRUE(allFinite) << "All table entries must be finite";
    EXPECT_NE(minVal, maxVal)
        << "The table should not be constant - pir depends on both r and r0";

    gsl_matrix_free(pirMatrix);
}

// -----------------------------------------------------------------------------
// Test 3: determinism - identical inputs produce identical tables.
// -----------------------------------------------------------------------------
void test_cpm_is_deterministic()
{
    std::cerr << "\n[TEST] test_cpm_is_deterministic\n"
              << "  Source file: src/reactions/create_pirMatrix.cpp\n"
              << "  Function:    create_pirMatrix()\n"
              << "  Scenario:    call the routine twice with identical arguments.\n"
              << "  Pass criteria: every element of the two tables is bit-identical\n"
              << "               (the quadrature uses no random input).\n";

    const double bindRadius { 1.0 };
    const double Dtot { 2.0 };
    const double kr { 5.0 };
    const double comRMax { 2.5 };
    const double RStepSize { 0.5 };
    const Parameters params { cpm_make_params(0.1) };

    const size_t dim { cpm_expected_dim(bindRadius, comRMax, RStepSize) };
    gsl_matrix* first = cpm_make_sentinel_matrix(dim + 1);
    gsl_matrix* second = cpm_make_sentinel_matrix(dim + 1);

    std::cerr << "  Filling two " << dim << " x " << dim << " tables...\n";
    create_pirMatrix(first, bindRadius, Dtot, kr, comRMax, RStepSize, params);
    create_pirMatrix(second, bindRadius, Dtot, kr, comRMax, RStepSize, params);

    size_t mismatches { 0 };
    for (size_t row = 0; row < dim; ++row) {
        for (size_t col = 0; col < dim; ++col) {
            if (gsl_matrix_get(first, row, col) != gsl_matrix_get(second, row, col))
                ++mismatches;
        }
    }
    std::cerr << "  Mismatching elements = " << mismatches << '\n';
    EXPECT_EQ(mismatches, 0u) << "create_pirMatrix should be deterministic";

    gsl_matrix_free(first);
    gsl_matrix_free(second);
}

// -----------------------------------------------------------------------------
// Test 4: a larger comRMax must produce a larger filled block.
// -----------------------------------------------------------------------------
void test_cpm_grid_grows_with_rmax()
{
    std::cerr << "\n[TEST] test_cpm_grid_grows_with_rmax\n"
              << "  Source file: src/reactions/create_pirMatrix.cpp\n"
              << "  Function:    create_pirMatrix()\n"
              << "  Scenario:    run with comRMax = 2.0 and again with comRMax = 3.0,\n"
              << "               using a sentinel to count how many rows were filled.\n"
              << "  Pass criteria: the larger comRMax fills strictly more rows/columns,\n"
              << "               and the count matches the loop bound arithmetic.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 1.0 };
    const double kr { 10.0 };
    const double RStepSize { 0.5 };
    const Parameters params { cpm_make_params(0.1) };

    const double smallRMax { 2.0 };
    const double largeRMax { 3.0 };

    const size_t smallDim { cpm_expected_dim(bindRadius, smallRMax, RStepSize) };
    const size_t largeDim { cpm_expected_dim(bindRadius, largeRMax, RStepSize) };
    std::cerr << "  Predicted dims: small = " << smallDim << ", large = " << largeDim << '\n';
    EXPECT_GT(largeDim, smallDim)
        << "A larger comRMax must yield more grid points";

    // Allocate both matrices at the *large* size so we can count filled cells
    // without the risk of an out-of-range write.
    const size_t alloc { largeDim + 2 };

    gsl_matrix* smallMat = cpm_make_sentinel_matrix(alloc);
    gsl_matrix* largeMat = cpm_make_sentinel_matrix(alloc);

    std::cerr << "  Calling create_pirMatrix twice (comRMax = " << smallRMax
              << " then " << largeRMax << ")...\n";
    create_pirMatrix(smallMat, bindRadius, Dtot, kr, smallRMax, RStepSize, params);
    create_pirMatrix(largeMat, bindRadius, Dtot, kr, largeRMax, RStepSize, params);

    // Count the number of leading columns whose element (0, col) was written.
    auto countFilledColumns = [alloc](const gsl_matrix* m) -> size_t {
        size_t filled { 0 };
        for (size_t col = 0; col < alloc; ++col) {
            if (gsl_matrix_get(m, 0, col) != kCpmSentinel)
                ++filled;
        }
        return filled;
    };

    const size_t smallFilled { countFilledColumns(smallMat) };
    const size_t largeFilled { countFilledColumns(largeMat) };
    std::cerr << "  Observed filled columns: small = " << smallFilled
              << ", large = " << largeFilled << '\n';

    EXPECT_EQ(smallFilled, smallDim)
        << "Filled column count should match the predicted grid size (small comRMax)";
    EXPECT_EQ(largeFilled, largeDim)
        << "Filled column count should match the predicted grid size (large comRMax)";
    EXPECT_GT(largeFilled, smallFilled)
        << "Growing comRMax must grow the populated block";

    gsl_matrix_free(smallMat);
    gsl_matrix_free(largeMat);
}

// -----------------------------------------------------------------------------
// Test 5: robustness against different physical inputs (rate and diffusion
//         constant) - the table must still contain only finite numbers.
// -----------------------------------------------------------------------------
void test_cpm_robust_to_rate_and_diffusion()
{
    std::cerr << "\n[TEST] test_cpm_robust_to_rate_and_diffusion\n"
              << "  Source file: src/reactions/create_pirMatrix.cpp\n"
              << "  Function:    create_pirMatrix()\n"
              << "  Scenario:    sweep a few (kr, Dtot, bindRadius) combinations.\n"
              << "  Pass criteria: for each combination the whole grid is written and\n"
              << "               contains only finite values (quadrature never blows up).\n";

    struct Case {
        double bindRadius;
        double Dtot;
        double kr;
        const char* label;
    };

    // A slow rate, a fast rate, and a larger binding radius / diffusion constant.
    const std::vector<Case> cases {
        { 1.0, 1.0, 0.5, "slow association rate" },
        { 1.0, 1.0, 100.0, "fast association rate" },
        { 2.0, 4.0, 10.0, "larger sigma and Dtot" },
    };

    const double RStepSize { 0.75 };
    const Parameters params { cpm_make_params(0.1) };

    for (const auto& oneCase : cases) {
        // Keep the grid small: comRMax = bindRadius + 2.
        const double comRMax { oneCase.bindRadius + 2.0 };
        const size_t dim { cpm_expected_dim(oneCase.bindRadius, comRMax, RStepSize) };
        gsl_matrix* pirMatrix = cpm_make_sentinel_matrix(dim + 1);

        std::cerr << "  Case '" << oneCase.label << "': bindRadius="
                  << oneCase.bindRadius << ", Dtot=" << oneCase.Dtot
                  << ", kr=" << oneCase.kr << ", comRMax=" << comRMax
                  << ", dim=" << dim << '\n';

        create_pirMatrix(pirMatrix, oneCase.bindRadius, oneCase.Dtot, oneCase.kr,
            comRMax, RStepSize, params);

        size_t bad { 0 };
        for (size_t row = 0; row < dim; ++row) {
            for (size_t col = 0; col < dim; ++col) {
                const double val = gsl_matrix_get(pirMatrix, row, col);
                if (val == kCpmSentinel || !std::isfinite(val))
                    ++bad;
            }
        }
        EXPECT_EQ(bad, 0u) << "Case '" << oneCase.label
                           << "': every grid point must be written with a finite value";
        std::cerr << "    -> bad elements = " << bad
                  << ", pir[0][0] = " << gsl_matrix_get(pirMatrix, 0, 0) << '\n';

        gsl_matrix_free(pirMatrix);
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each helper runs inside its own TEST so that a failure
// in one scenario does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(CreatePirMatrix, FillsExpectedBlock) { test_cpm_fills_expected_block(); }
TEST(CreatePirMatrix, ValuesVaryAcrossGrid) { test_cpm_values_vary_across_grid(); }
TEST(CreatePirMatrix, IsDeterministic) { test_cpm_is_deterministic(); }
TEST(CreatePirMatrix, GridGrowsWithRmax) { test_cpm_grid_grows_with_rmax(); }
TEST(CreatePirMatrix, RobustToRateAndDiffusion) { test_cpm_robust_to_rate_and_diffusion(); }