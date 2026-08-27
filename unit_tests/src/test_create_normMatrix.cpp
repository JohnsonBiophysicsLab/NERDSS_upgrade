/*! \file test_create_normMatrix.cpp
 *
 * ### Unit test for src/reactions/create_normMatrix.cpp
 *
 * The single function under test is:
 *
 *     void create_normMatrix(gsl_matrix*& normMatrix, double bindRadius,
 *                            double Dtot, double kr, double comRMax,
 *                            double RStepSize, const Parameters& params)
 *
 * Behaviour of the routine (read directly from the implementation):
 *
 *   - It does **not** allocate `normMatrix`; the caller must hand in a matrix
 *     that already has at least 2 rows and enough columns for the whole
 *     r0-grid.  Writing past the end would make GSL abort the process, so the
 *     tests below mirror the loop arithmetic exactly and then add a safety
 *     margin of extra columns.
 *   - It walks an r0 grid starting exactly at `bindRadius` and stepping by
 *     `RStepSize` while `RIndex <= comRMax + RStepSize` (note the "+ one extra
 *     step" in the loop guard - the grid deliberately runs one step past
 *     comRMax).
 *   - For grid point `itr` it stores the r0 value in row 0 and the value of
 *     \f$\int_{bindRadius}^{\infty} norm\_function(x)\,dx\f$ in row 1, using
 *     gsl_integration_qagiu with epsAbs = epsRel = 1e-6.
 *   - The integrand parameters are (a = bindRadius, D = Dtot, k = kr,
 *     t = params.timeStep, r0 = current grid point).  `rho` and `r` are left
 *     at their default zero values by IntegrandParams.
 *
 * The tests therefore verify the *grid layout*, the *stored integral values*
 * (by repeating the integration independently with the same integrand), the
 * fact that untouched columns are left alone, and determinism.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace {

//! Sentinel value pre-written into every matrix element.  Any element still
//! holding this value after the call was *not* touched by create_normMatrix().
constexpr double kCnmSentinel = -98765.4321;

/*!
 * \brief Reproduce, bit for bit, the loop arithmetic of create_normMatrix().
 *
 * create_normMatrix() accumulates `RIndex += RStepSize` in a while loop guarded
 * by `RIndex <= comRMax + RStepSize`.  Re-using the identical accumulation here
 * means our expected grid suffers exactly the same floating point drift as the
 * production code, so the column count and the individual r0 values match
 * exactly rather than approximately.
 *
 * \return number of grid points (== number of columns the routine will write).
 */
int cnm_expected_grid(double bindRadius, double comRMax, double RStepSize,
    std::vector<double>& expectedR0)
{
    expectedR0.clear();
    double RIndex { bindRadius };
    while (RIndex <= comRMax + RStepSize) {
        expectedR0.push_back(RIndex);
        RIndex += RStepSize;
    }
    return static_cast<int>(expectedR0.size());
}

/*!
 * \brief Independently evaluate the same semi-infinite integral the routine does.
 *
 * Uses the very same integrand (norm_function), the same lower bound
 * (bindRadius) and the same absolute/relative tolerances (1e-6).  The GSL
 * subdivision limit differs (a much smaller workspace is used here to keep the
 * test lightweight), which is harmless as long as the integral converges in far
 * fewer subdivisions than the limit - in that case the returned value is
 * identical.
 */
double cnm_reference_integral(double bindRadius, double Dtot, double kr, double timeStep, double r0)
{
    IntegrandParams intParams; // all members are value-initialised (rho = 0, r = 0)
    intParams.a = bindRadius;
    intParams.D = Dtot;
    intParams.k = kr;
    intParams.t = timeStep;
    intParams.r0 = r0;

    gsl_function F;
    F.function = &norm_function;
    F.params = reinterpret_cast<void*>(&intParams);

    const std::size_t limit { 10000 };
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);

    double result { 0.0 };
    double error { 0.0 };
    gsl_integration_qagiu(&F, bindRadius, 1e-6, 1e-6, limit, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

/*!
 * \brief Allocate a 2 x (n + slack) matrix and stamp every element with the
 *        sentinel so we can later see which columns were written.
 */
gsl_matrix* cnm_make_matrix(int numGridPoints, int slack)
{
    gsl_matrix* m = gsl_matrix_alloc(2, static_cast<std::size_t>(numGridPoints + slack));
    gsl_matrix_set_all(m, kCnmSentinel);
    return m;
}

/*!
 * \brief Count how many columns of row 0 differ from the sentinel.
 */
int cnm_count_written_columns(const gsl_matrix* m)
{
    int written { 0 };
    for (std::size_t j = 0; j < m->size2; ++j) {
        if (gsl_matrix_get(m, 0, j) != kCnmSentinel)
            written = static_cast<int>(j) + 1;
    }
    return written;
}

/*! RAII helper: switch the GSL error handler off for the duration of a test.
 *
 * By default GSL aborts the whole process when e.g. an integration fails to
 * reach the requested tolerance.  That would take down the entire gtest binary,
 * so the handler is disabled while the routine runs and restored afterwards.
 */
class CnmErrorHandlerGuard {
public:
    CnmErrorHandlerGuard() { oldHandler_ = gsl_set_error_handler_off(); }
    ~CnmErrorHandlerGuard() { gsl_set_error_handler(oldHandler_); }

private:
    gsl_error_handler_t* oldHandler_ { nullptr };
};

} // namespace

// -----------------------------------------------------------------------------
// Test 1: grid layout - row 0 must hold exactly bindRadius, bindRadius+step, ...
//         and the number of written columns must follow the loop guard
//         "RIndex <= comRMax + RStepSize" (i.e. one step past comRMax).
// -----------------------------------------------------------------------------
void test_cnm_grid_layout()
{
    std::cerr << "\n[TEST] test_cnm_grid_layout\n"
              << "  Source file:   src/reactions/create_normMatrix.cpp\n"
              << "  Function:      create_normMatrix()\n"
              << "  Checking:      row 0 stores the r0 grid, starting at bindRadius\n"
              << "                 and stepping by RStepSize, running one step past\n"
              << "                 comRMax because of the '<= comRMax + RStepSize' guard.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 1.0 };
    const double kr { 10.0 };
    const double comRMax { 3.0 };
    const double RStepSize { 0.5 };

    Parameters params;
    params.timeStep = 0.1; // only timeStep is consumed by the routine

    // Work out the expected grid using the identical accumulation.
    std::vector<double> expectedR0;
    const int expectedCols = cnm_expected_grid(bindRadius, comRMax, RStepSize, expectedR0);
    std::cerr << "  Expected number of grid points: " << expectedCols << '\n';

    // Allocate with slack so a mis-sized loop cannot walk off the matrix.
    gsl_matrix* normMatrix = cnm_make_matrix(expectedCols, 4);

    {
        CnmErrorHandlerGuard guard;
        std::cerr << "  Calling create_normMatrix(bindRadius=" << bindRadius
                  << ", Dtot=" << Dtot << ", kr=" << kr << ", comRMax=" << comRMax
                  << ", RStepSize=" << RStepSize << ", timeStep=" << params.timeStep << ")\n";
        create_normMatrix(normMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);
    }

    // --- number of columns actually written -------------------------------
    const int written = cnm_count_written_columns(normMatrix);
    std::cerr << "  Columns written by the routine: " << written << '\n';
    EXPECT_EQ(written, expectedCols)
        << "create_normMatrix should write exactly one column per grid point";

    // --- the r0 values stored in row 0 ------------------------------------
    for (int j = 0; j < expectedCols && j < written; ++j) {
        const double got = gsl_matrix_get(normMatrix, 0, static_cast<std::size_t>(j));
        std::cerr << "    col " << j << ": r0 stored = " << got
                  << " (expected " << expectedR0[j] << ")\n";
        EXPECT_DOUBLE_EQ(got, expectedR0[j])
            << "Row 0, column " << j << " should hold the r0 grid value";
    }

    // The very first grid point must be the binding radius itself.
    if (written > 0) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(normMatrix, 0, 0), bindRadius)
            << "The grid must start exactly at bindRadius";
    }

    // The grid must be strictly increasing.
    for (int j = 1; j < written; ++j) {
        EXPECT_GT(gsl_matrix_get(normMatrix, 0, static_cast<std::size_t>(j)),
            gsl_matrix_get(normMatrix, 0, static_cast<std::size_t>(j - 1)))
            << "r0 grid must increase monotonically (column " << j << ')';
    }

    // The last grid point must be beyond comRMax (the loop runs one step past).
    if (written > 0) {
        const double last = gsl_matrix_get(normMatrix, 0, static_cast<std::size_t>(written - 1));
        std::cerr << "  Last grid point = " << last << " (comRMax = " << comRMax << ")\n";
        EXPECT_GT(last, comRMax)
            << "The loop guard 'RIndex <= comRMax + RStepSize' means the grid extends past comRMax";
        EXPECT_LE(last, comRMax + RStepSize + 1e-12)
            << "The grid must not extend more than one step past comRMax";
    }

    gsl_matrix_free(normMatrix);
}

// -----------------------------------------------------------------------------
// Test 2: row 1 must equal the semi-infinite integral of norm_function taken
//         from bindRadius to infinity, for the r0 sitting in row 0.
// -----------------------------------------------------------------------------
void test_cnm_matches_direct_integration()
{
    std::cerr << "\n[TEST] test_cnm_matches_direct_integration\n"
              << "  Source file:   src/reactions/create_normMatrix.cpp\n"
              << "  Function:      create_normMatrix()\n"
              << "  Checking:      row 1 holds qagiu(norm_function, bindRadius -> inf)\n"
              << "                 evaluated with a = bindRadius, D = Dtot, k = kr,\n"
              << "                 t = params.timeStep and r0 = the row-0 grid value.\n"
              << "  Pass criteria: values agree with an independent integration of the\n"
              << "                 same integrand to within a 1e-6 relative tolerance.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 1.0 };
    const double kr { 10.0 };
    const double comRMax { 2.5 };
    const double RStepSize { 0.5 };

    Parameters params;
    params.timeStep = 0.1;

    std::vector<double> expectedR0;
    const int expectedCols = cnm_expected_grid(bindRadius, comRMax, RStepSize, expectedR0);
    gsl_matrix* normMatrix = cnm_make_matrix(expectedCols, 4);

    {
        CnmErrorHandlerGuard guard;
        std::cerr << "  Calling create_normMatrix with " << expectedCols << " grid points...\n";
        create_normMatrix(normMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);
    }

    const int written = cnm_count_written_columns(normMatrix);
    ASSERT_GT(written, 0) << "The routine must write at least one grid point";

    for (int j = 0; j < written; ++j) {
        const double r0 = gsl_matrix_get(normMatrix, 0, static_cast<std::size_t>(j));
        const double stored = gsl_matrix_get(normMatrix, 1, static_cast<std::size_t>(j));

        double reference { 0.0 };
        {
            CnmErrorHandlerGuard guard;
            reference = cnm_reference_integral(bindRadius, Dtot, kr, params.timeStep, r0);
        }

        std::cerr << "    col " << j << ": r0 = " << r0
                  << ", stored norm = " << stored
                  << ", independent integral = " << reference << '\n';

        // The stored value must be a usable finite number.
        EXPECT_TRUE(std::isfinite(stored))
            << "Row 1, column " << j << " must be finite (not NaN/Inf)";

        // Compare on a relative scale so both tiny and large integrals work.
        const double scale = std::max(1.0, std::fabs(reference));
        EXPECT_NEAR(stored, reference, 1e-6 * scale)
            << "Row 1, column " << j << " should reproduce the norm_function integral";
    }

    gsl_matrix_free(normMatrix);
}

// -----------------------------------------------------------------------------
// Test 3: columns beyond the r0 grid must be left completely untouched, and the
//         routine must be deterministic (no RNG, no hidden state).
// -----------------------------------------------------------------------------
void test_cnm_untouched_columns_and_determinism()
{
    std::cerr << "\n[TEST] test_cnm_untouched_columns_and_determinism\n"
              << "  Source file:   src/reactions/create_normMatrix.cpp\n"
              << "  Function:      create_normMatrix()\n"
              << "  Checking:      (a) columns past the grid still hold the sentinel,\n"
              << "                 (b) two identical calls produce identical matrices.\n";

    const double bindRadius { 1.5 };
    const double Dtot { 0.5 };
    const double kr { 5.0 };
    const double comRMax { 3.0 };
    const double RStepSize { 0.5 };

    Parameters params;
    params.timeStep = 0.2;

    std::vector<double> expectedR0;
    const int expectedCols = cnm_expected_grid(bindRadius, comRMax, RStepSize, expectedR0);
    const int slack = 5;

    gsl_matrix* first = cnm_make_matrix(expectedCols, slack);
    gsl_matrix* second = cnm_make_matrix(expectedCols, slack);

    {
        CnmErrorHandlerGuard guard;
        std::cerr << "  First call to create_normMatrix (" << expectedCols << " grid points)...\n";
        create_normMatrix(first, bindRadius, Dtot, kr, comRMax, RStepSize, params);
        std::cerr << "  Second call with identical arguments...\n";
        create_normMatrix(second, bindRadius, Dtot, kr, comRMax, RStepSize, params);
    }

    // --- (a) trailing columns must still hold the sentinel -----------------
    const int written = cnm_count_written_columns(first);
    std::cerr << "  Columns written = " << written << ", matrix has "
              << first->size2 << " columns; checking the trailing "
              << (static_cast<int>(first->size2) - written) << " are untouched.\n";
    for (std::size_t j = static_cast<std::size_t>(written); j < first->size2; ++j) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(first, 0, j), kCnmSentinel)
            << "Row 0, column " << j << " lies past the grid and must not be written";
        EXPECT_DOUBLE_EQ(gsl_matrix_get(first, 1, j), kCnmSentinel)
            << "Row 1, column " << j << " lies past the grid and must not be written";
    }

    // --- (b) determinism ---------------------------------------------------
    for (std::size_t j = 0; j < first->size2; ++j) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(first, 0, j), gsl_matrix_get(second, 0, j))
            << "Row 0, column " << j << " differs between two identical calls";
        EXPECT_DOUBLE_EQ(gsl_matrix_get(first, 1, j), gsl_matrix_get(second, 1, j))
            << "Row 1, column " << j << " differs between two identical calls";
    }
    std::cerr << "  Both calls produced identical matrices (deterministic).\n";

    gsl_matrix_free(first);
    gsl_matrix_free(second);
}

// -----------------------------------------------------------------------------
// Test 4: the grid start is driven purely by bindRadius, and changing the
//         time step changes the stored integrals (the integrand depends on t).
// -----------------------------------------------------------------------------
void test_cnm_parameter_sensitivity()
{
    std::cerr << "\n[TEST] test_cnm_parameter_sensitivity\n"
              << "  Source file:   src/reactions/create_normMatrix.cpp\n"
              << "  Function:      create_normMatrix()\n"
              << "  Checking:      (a) the r0 grid begins at whatever bindRadius is given,\n"
              << "                 (b) params.timeStep really is forwarded to the integrand,\n"
              << "                     so a different time step yields different norm values.\n";

    const double Dtot { 1.0 };
    const double kr { 10.0 };
    const double comRMax { 2.0 };
    const double RStepSize { 0.5 };

    // --- (a) a shifted binding radius shifts the whole grid ---------------
    const double bindRadiusA { 1.0 };
    const double bindRadiusB { 2.0 };

    Parameters params;
    params.timeStep = 0.1;

    std::vector<double> gridA;
    std::vector<double> gridB;
    const int colsA = cnm_expected_grid(bindRadiusA, comRMax, RStepSize, gridA);
    const int colsB = cnm_expected_grid(bindRadiusB, comRMax, RStepSize, gridB);

    gsl_matrix* matA = cnm_make_matrix(colsA, 4);
    gsl_matrix* matB = cnm_make_matrix(colsB, 4);

    {
        CnmErrorHandlerGuard guard;
        std::cerr << "  Call with bindRadius = " << bindRadiusA << " ...\n";
        create_normMatrix(matA, bindRadiusA, Dtot, kr, comRMax, RStepSize, params);
        std::cerr << "  Call with bindRadius = " << bindRadiusB << " ...\n";
        create_normMatrix(matB, bindRadiusB, Dtot, kr, comRMax, RStepSize, params);
    }

    ASSERT_GT(cnm_count_written_columns(matA), 0) << "grid A must be non-empty";
    ASSERT_GT(cnm_count_written_columns(matB), 0) << "grid B must be non-empty";

    std::cerr << "    grid A starts at " << gsl_matrix_get(matA, 0, 0)
              << ", grid B starts at " << gsl_matrix_get(matB, 0, 0) << '\n';
    EXPECT_DOUBLE_EQ(gsl_matrix_get(matA, 0, 0), bindRadiusA)
        << "The first grid point must equal the supplied bindRadius";
    EXPECT_DOUBLE_EQ(gsl_matrix_get(matB, 0, 0), bindRadiusB)
        << "The first grid point must equal the supplied bindRadius";

    // A larger bindRadius reaches comRMax sooner, so it needs no more columns.
    std::cerr << "    grid A columns = " << colsA << ", grid B columns = " << colsB << '\n';
    EXPECT_LE(colsB, colsA)
        << "A larger bindRadius must not produce more grid points for the same comRMax";

    gsl_matrix_free(matB);

    // --- (b) a different time step must alter the integrals ----------------
    Parameters paramsLongStep;
    paramsLongStep.timeStep = 10.0 * params.timeStep;

    gsl_matrix* matLong = cnm_make_matrix(colsA, 4);
    {
        CnmErrorHandlerGuard guard;
        std::cerr << "  Call with timeStep = " << paramsLongStep.timeStep
                  << " (10x the original) ...\n";
        create_normMatrix(matLong, bindRadiusA, Dtot, kr, comRMax, RStepSize, paramsLongStep);
    }

    // The r0 grid does not depend on the time step at all.
    const int writtenA = cnm_count_written_columns(matA);
    for (int j = 0; j < writtenA; ++j) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(matA, 0, static_cast<std::size_t>(j)),
            gsl_matrix_get(matLong, 0, static_cast<std::size_t>(j)))
            << "Row 0 (the r0 grid) must be independent of params.timeStep";
    }

    // At least one integral value must change, proving timeStep is forwarded.
    bool anyDifferent { false };
    for (int j = 0; j < writtenA; ++j) {
        const double shortStep = gsl_matrix_get(matA, 1, static_cast<std::size_t>(j));
        const double longStep = gsl_matrix_get(matLong, 1, static_cast<std::size_t>(j));
        std::cerr << "    col " << j << ": norm(t=" << params.timeStep << ") = " << shortStep
                  << ", norm(t=" << paramsLongStep.timeStep << ") = " << longStep << '\n';
        if (shortStep != longStep)
            anyDifferent = true;
    }
    EXPECT_TRUE(anyDifferent)
        << "params.timeStep is passed into the integrand, so the stored norms must change";

    gsl_matrix_free(matA);
    gsl_matrix_free(matLong);
}

// -----------------------------------------------------------------------------
// Test 5: degenerate grid - comRMax equal to bindRadius still yields the two
//         grid points implied by the "<= comRMax + RStepSize" guard.
// -----------------------------------------------------------------------------
void test_cnm_minimal_grid()
{
    std::cerr << "\n[TEST] test_cnm_minimal_grid\n"
              << "  Source file:   src/reactions/create_normMatrix.cpp\n"
              << "  Function:      create_normMatrix()\n"
              << "  Checking:      with comRMax == bindRadius the loop still runs the\n"
              << "                 extra step allowed by '<= comRMax + RStepSize'.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 1.0 };
    const double kr { 10.0 };
    const double comRMax { 1.0 }; // identical to bindRadius
    const double RStepSize { 0.5 };

    Parameters params;
    params.timeStep = 0.1;

    std::vector<double> expectedR0;
    const int expectedCols = cnm_expected_grid(bindRadius, comRMax, RStepSize, expectedR0);
    std::cerr << "  Expected grid points for the degenerate case: " << expectedCols << '\n';
    // bindRadius (<= 1.5) and bindRadius+step == 1.5 (<= 1.5) both qualify.
    EXPECT_GE(expectedCols, 2)
        << "The loop guard allows one full step past comRMax, so >= 2 points are expected";

    gsl_matrix* normMatrix = cnm_make_matrix(expectedCols, 4);
    {
        CnmErrorHandlerGuard guard;
        std::cerr << "  Calling create_normMatrix on the degenerate grid...\n";
        create_normMatrix(normMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);
    }

    const int written = cnm_count_written_columns(normMatrix);
    std::cerr << "  Columns written = " << written << '\n';
    EXPECT_EQ(written, expectedCols)
        << "Even for comRMax == bindRadius the grid follows the loop guard exactly";

    for (int j = 0; j < written; ++j) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(normMatrix, 0, static_cast<std::size_t>(j)), expectedR0[j])
            << "Row 0, column " << j << " should hold the expected r0 value";
        EXPECT_TRUE(std::isfinite(gsl_matrix_get(normMatrix, 1, static_cast<std::size_t>(j))))
            << "Row 1, column " << j << " should hold a finite integral value";
    }

    gsl_matrix_free(normMatrix);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each named helper is executed inside its own TEST so a
// failure in one does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(CreateNormMatrix, GridLayout) { test_cnm_grid_layout(); }
TEST(CreateNormMatrix, MatchesDirectIntegration) { test_cnm_matches_direct_integration(); }
TEST(CreateNormMatrix, UntouchedColumnsAndDeterminism) { test_cnm_untouched_columns_and_determinism(); }
TEST(CreateNormMatrix, ParameterSensitivity) { test_cnm_parameter_sensitivity(); }
TEST(CreateNormMatrix, MinimalGrid) { test_cnm_minimal_grid(); }