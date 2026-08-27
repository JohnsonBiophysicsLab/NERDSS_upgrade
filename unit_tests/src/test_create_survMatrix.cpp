/*! \file test_create_survMatrix.cpp
 *
 * ### Unit test for ../src/reactions/create_survMatrix.cpp
 *
 * The single function under test is
 *
 *     void create_survMatrix(gsl_matrix*& survMatrix, double bindRadius, double Dtot,
 *                            double kr, double comRMax, double RStepSize,
 *                            const Parameters& params);
 *
 * Behaviour of the function, read directly from the implementation:
 *
 *  - It walks a local variable `RIndex` starting at `bindRadius` and stepping by
 *    `RStepSize` for as long as `RIndex <= comRMax + RStepSize`.
 *  - For every step it writes the current separation into row 0 of `survMatrix`
 *    (column `ctr`) and the value returned by `integrator(...)` into row 1.
 *  - The `kr < 1.0/0.0` test is `kr < +infinity`, which is true for **every**
 *    finite rate constant, so the "1.0 - result" branch is only reachable when
 *    `kr` is exactly `+inf` (or NaN).  We deliberately do not feed an infinite
 *    rate constant to the adaptive GSL quadrature - that is not a supported
 *    input and the integration is not guaranteed to terminate - so that branch
 *    is documented here rather than executed.
 *  - The caller owns the matrix: create_survMatrix never allocates or resizes
 *    it, so every test below allocates a matrix that is big enough and pre-fills
 *    it with a sentinel value so we can see exactly which cells were touched.
 *
 * Because the actual survival probabilities come out of an adaptive numerical
 * integration, the assertions here concentrate on the parts of the contract that
 * are exactly defined by the code (which cells are written, what the separation
 * grid is, determinism) plus finiteness/"was actually written" checks on the
 * integrated values.  Diagnostic values are printed to stderr for eyeballing.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_matrix.h>
#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <iostream>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace so the names cannot collide with the rest
// of the test suite).
// -----------------------------------------------------------------------------
namespace {

//! Value used to pre-fill the matrix so untouched cells are recognisable.
constexpr double kCsmSentinel = -999.0;

/*! \brief Reproduce, bit for bit, the loop bound used inside create_survMatrix.
 *
 * The production loop accumulates `RIndex += RStepSize` in double precision, so
 * the only reliable way to predict how many columns get filled is to perform
 * exactly the same accumulation here.
 */
size_t csm_expected_columns(double bindRadius, double comRMax, double RStepSize)
{
    size_t ctr { 0 };
    double RIndex { bindRadius };
    while (RIndex <= comRMax + RStepSize) {
        ++ctr;
        RIndex += RStepSize;
    }
    return ctr;
}

/*! \brief Allocate a 2 x cols matrix and fill every cell with the sentinel. */
gsl_matrix* csm_alloc_sentinel_matrix(size_t cols)
{
    gsl_matrix* mat = gsl_matrix_alloc(2, cols);
    gsl_matrix_set_all(mat, kCsmSentinel);
    return mat;
}

/*! \brief Count the leading columns of row 0 that were overwritten. */
size_t csm_count_written_columns(const gsl_matrix* mat)
{
    size_t written { 0 };
    for (size_t col = 0; col < mat->size2; ++col) {
        if (gsl_matrix_get(mat, 0, col) == kCsmSentinel)
            break;
        ++written;
    }
    return written;
}

/*! \brief Minimal Parameters object: create_survMatrix only reads timeStep. */
Parameters csm_make_params(double timeStep)
{
    Parameters params;
    params.timeStep = timeStep; // microseconds
    return params;
}

/*! \brief Dump a filled matrix to stderr so a human can inspect the table. */
void csm_dump(const gsl_matrix* mat, size_t nWritten, const char* label)
{
    std::cerr << "    " << label << " (" << nWritten << " columns):\n";
    for (size_t col = 0; col < nWritten; ++col) {
        std::cerr << "      r0 = " << gsl_matrix_get(mat, 0, col)
                  << "   surv = " << gsl_matrix_get(mat, 1, col) << '\n';
    }
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the separation grid stored in row 0 must be exactly
//         bindRadius, bindRadius + RStepSize, bindRadius + 2*RStepSize, ...
// -----------------------------------------------------------------------------
void test_csm_row0_separation_grid()
{
    std::cerr << "\n[TEST] test_csm_row0_separation_grid\n"
              << "  Source file: src/reactions/create_survMatrix.cpp\n"
              << "  Function:    create_survMatrix()\n"
              << "  Scenario:    fill a survival table for a 2D reaction with\n"
              << "               sigma = 1 nm, Dtot = 0.1 nm^2/us, ka = 10 nm^2/us,\n"
              << "               dt = 0.1 us, Rmax = 1.6 nm, step = 0.05 nm.\n"
              << "  Pass:        row 0 holds the accumulated separation grid,\n"
              << "               starting exactly at the binding radius.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 0.1 };
    const double kr { 10.0 };
    const double comRMax { 1.6 };
    const double RStepSize { 0.05 };
    const Parameters params { csm_make_params(0.1) };

    const size_t expectedCols { csm_expected_columns(bindRadius, comRMax, RStepSize) };
    std::cerr << "  Expected number of filled columns: " << expectedCols << '\n';

    // Allocate a matrix with slack so we can also detect over-writes.
    gsl_matrix* survMatrix = csm_alloc_sentinel_matrix(expectedCols + 8);

    create_survMatrix(survMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);

    const size_t written { csm_count_written_columns(survMatrix) };
    EXPECT_EQ(written, expectedCols)
        << "create_survMatrix should fill exactly the number of columns implied by its loop";

    // Re-create the same accumulation to compare bit-for-bit.
    double RIndex { bindRadius };
    for (size_t col = 0; col < written; ++col) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(survMatrix, 0, col), RIndex)
            << "row 0, column " << col << " should hold the accumulated separation";
        RIndex += RStepSize;
    }

    // The very first entry is, by construction, the binding radius itself.
    if (written > 0) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(survMatrix, 0, 0), bindRadius)
            << "the table must start at contact (r0 == bindRadius)";
    }

    csm_dump(survMatrix, written, "survMatrix");
    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// Test 2: every integrated value written into row 1 must be a real number and
//         must actually have been written (i.e. differ from the sentinel).
// -----------------------------------------------------------------------------
void test_csm_row1_values_are_finite()
{
    std::cerr << "\n[TEST] test_csm_row1_values_are_finite\n"
              << "  Source file: src/reactions/create_survMatrix.cpp\n"
              << "  Function:    create_survMatrix() -> integrator()/survival_function()\n"
              << "  Pass:        every row-1 entry is written (!= sentinel) and finite\n"
              << "               (no NaN / no infinity coming out of the quadrature).\n";

    const double bindRadius { 1.0 };
    const double Dtot { 0.1 };
    const double kr { 10.0 };
    const double comRMax { 1.4 };
    const double RStepSize { 0.05 };
    const Parameters params { csm_make_params(0.1) };

    const size_t expectedCols { csm_expected_columns(bindRadius, comRMax, RStepSize) };
    gsl_matrix* survMatrix = csm_alloc_sentinel_matrix(expectedCols + 4);

    create_survMatrix(survMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);

    const size_t written { csm_count_written_columns(survMatrix) };
    EXPECT_GT(written, 0u) << "at least one column must be produced for a valid Rmax";

    for (size_t col = 0; col < written; ++col) {
        const double value { gsl_matrix_get(survMatrix, 1, col) };
        EXPECT_NE(value, kCsmSentinel)
            << "row 1, column " << col << " was never written by create_survMatrix";
        EXPECT_TRUE(std::isfinite(value))
            << "row 1, column " << col << " is not finite (value = " << value << ')';
    }

    csm_dump(survMatrix, written, "survMatrix");
    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// Test 3: the function must not touch any cell beyond the last column it fills.
// -----------------------------------------------------------------------------
void test_csm_does_not_write_past_last_column()
{
    std::cerr << "\n[TEST] test_csm_does_not_write_past_last_column\n"
              << "  Source file: src/reactions/create_survMatrix.cpp\n"
              << "  Function:    create_survMatrix()\n"
              << "  Scenario:    matrix is allocated with 6 spare columns.\n"
              << "  Pass:        the spare columns of BOTH rows still hold the sentinel.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 0.1 };
    const double kr { 10.0 };
    const double comRMax { 1.3 };
    const double RStepSize { 0.05 };
    const Parameters params { csm_make_params(0.1) };

    const size_t expectedCols { csm_expected_columns(bindRadius, comRMax, RStepSize) };
    const size_t spare { 6 };
    gsl_matrix* survMatrix = csm_alloc_sentinel_matrix(expectedCols + spare);

    create_survMatrix(survMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);

    std::cerr << "  Filled columns: " << expectedCols
              << ", checking columns " << expectedCols << " .. "
              << (expectedCols + spare - 1) << " are untouched\n";

    for (size_t col = expectedCols; col < expectedCols + spare; ++col) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(survMatrix, 0, col), kCsmSentinel)
            << "row 0, column " << col << " should not have been written";
        EXPECT_DOUBLE_EQ(gsl_matrix_get(survMatrix, 1, col), kCsmSentinel)
            << "row 1, column " << col << " should not have been written";
    }

    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// Test 4: the number of columns produced scales the way the loop says it does -
//         halving the step size roughly doubles the table, and increasing Rmax
//         lengthens the table.
// -----------------------------------------------------------------------------
void test_csm_column_count_follows_step_and_rmax()
{
    std::cerr << "\n[TEST] test_csm_column_count_follows_step_and_rmax\n"
              << "  Source file: src/reactions/create_survMatrix.cpp\n"
              << "  Function:    create_survMatrix()\n"
              << "  Pass:        column count matches the loop bound for two different\n"
              << "               step sizes, a finer step gives more columns, and a\n"
              << "               larger comRMax gives more columns.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 0.1 };
    const double kr { 10.0 };
    const Parameters params { csm_make_params(0.1) };

    // (a) coarse step
    const double coarseStep { 0.10 };
    const double comRMax { 1.4 };
    const size_t coarseExpected { csm_expected_columns(bindRadius, comRMax, coarseStep) };
    gsl_matrix* coarse = csm_alloc_sentinel_matrix(coarseExpected + 4);
    create_survMatrix(coarse, bindRadius, Dtot, kr, comRMax, coarseStep, params);
    const size_t coarseWritten { csm_count_written_columns(coarse) };
    std::cerr << "  step = " << coarseStep << " -> " << coarseWritten
              << " columns (expected " << coarseExpected << ")\n";
    EXPECT_EQ(coarseWritten, coarseExpected);

    // (b) fine step over the same range
    const double fineStep { 0.05 };
    const size_t fineExpected { csm_expected_columns(bindRadius, comRMax, fineStep) };
    gsl_matrix* fine = csm_alloc_sentinel_matrix(fineExpected + 4);
    create_survMatrix(fine, bindRadius, Dtot, kr, comRMax, fineStep, params);
    const size_t fineWritten { csm_count_written_columns(fine) };
    std::cerr << "  step = " << fineStep << " -> " << fineWritten
              << " columns (expected " << fineExpected << ")\n";
    EXPECT_EQ(fineWritten, fineExpected);
    EXPECT_GT(fineWritten, coarseWritten)
        << "a finer RStepSize must produce a longer table over the same range";

    // (c) larger Rmax with the coarse step
    const double biggerRMax { 2.0 };
    const size_t biggerExpected { csm_expected_columns(bindRadius, biggerRMax, coarseStep) };
    gsl_matrix* bigger = csm_alloc_sentinel_matrix(biggerExpected + 4);
    create_survMatrix(bigger, bindRadius, Dtot, kr, biggerRMax, coarseStep, params);
    const size_t biggerWritten { csm_count_written_columns(bigger) };
    std::cerr << "  comRMax = " << biggerRMax << " -> " << biggerWritten
              << " columns (expected " << biggerExpected << ")\n";
    EXPECT_EQ(biggerWritten, biggerExpected);
    EXPECT_GT(biggerWritten, coarseWritten)
        << "a larger comRMax must produce a longer table for the same step size";

    gsl_matrix_free(coarse);
    gsl_matrix_free(fine);
    gsl_matrix_free(bigger);
}

// -----------------------------------------------------------------------------
// Test 5: the routine is deterministic - identical inputs give an identical
//         table (it uses no random numbers and keeps no hidden state).
// -----------------------------------------------------------------------------
void test_csm_is_deterministic()
{
    std::cerr << "\n[TEST] test_csm_is_deterministic\n"
              << "  Source file: src/reactions/create_survMatrix.cpp\n"
              << "  Function:    create_survMatrix()\n"
              << "  Pass:        two calls with identical arguments produce bit-identical\n"
              << "               matrices (no RNG, no carried-over state).\n";

    const double bindRadius { 1.0 };
    const double Dtot { 0.1 };
    const double kr { 10.0 };
    const double comRMax { 1.4 };
    const double RStepSize { 0.05 };
    const Parameters params { csm_make_params(0.1) };

    const size_t expectedCols { csm_expected_columns(bindRadius, comRMax, RStepSize) };
    gsl_matrix* first = csm_alloc_sentinel_matrix(expectedCols + 2);
    gsl_matrix* second = csm_alloc_sentinel_matrix(expectedCols + 2);

    create_survMatrix(first, bindRadius, Dtot, kr, comRMax, RStepSize, params);
    create_survMatrix(second, bindRadius, Dtot, kr, comRMax, RStepSize, params);

    const size_t written { csm_count_written_columns(first) };
    EXPECT_EQ(written, csm_count_written_columns(second))
        << "both calls should fill the same number of columns";

    for (size_t col = 0; col < written; ++col) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(first, 0, col), gsl_matrix_get(second, 0, col))
            << "row 0 differs between two identical calls at column " << col;
        EXPECT_DOUBLE_EQ(gsl_matrix_get(first, 1, col), gsl_matrix_get(second, 1, col))
            << "row 1 differs between two identical calls at column " << col;
    }

    std::cerr << "  Compared " << written << " columns; tables are identical.\n";

    gsl_matrix_free(first);
    gsl_matrix_free(second);
}

// -----------------------------------------------------------------------------
// Test 6: the tabulated quantity really depends on the separation r0 and on the
//         intrinsic rate constant - i.e. the integrator is being driven with the
//         per-column parameters and is not returning a constant.
// -----------------------------------------------------------------------------
void test_csm_values_depend_on_inputs()
{
    std::cerr << "\n[TEST] test_csm_values_depend_on_inputs\n"
              << "  Source file: src/reactions/create_survMatrix.cpp\n"
              << "  Function:    create_survMatrix() -> integrator()/survival_function()\n"
              << "  Pass:        (a) the value at contact differs from the value at the\n"
              << "                   largest separation (r0 is fed to the integrand), and\n"
              << "               (b) changing the intrinsic rate ka changes the table.\n";

    const double bindRadius { 1.0 };
    const double Dtot { 0.1 };
    const double comRMax { 1.4 };
    const double RStepSize { 0.05 };
    const Parameters params { csm_make_params(0.1) };

    const size_t expectedCols { csm_expected_columns(bindRadius, comRMax, RStepSize) };

    // (a) small rate constant
    const double krSmall { 1.0 };
    gsl_matrix* small = csm_alloc_sentinel_matrix(expectedCols + 2);
    create_survMatrix(small, bindRadius, Dtot, krSmall, comRMax, RStepSize, params);
    const size_t written { csm_count_written_columns(small) };

    if (written >= 2) {
        const double atContact { gsl_matrix_get(small, 1, 0) };
        const double atFarEnd { gsl_matrix_get(small, 1, written - 1) };
        std::cerr << "  ka = " << krSmall << ": value at r0 = " << bindRadius << " is "
                  << atContact << ", value at r0 = "
                  << gsl_matrix_get(small, 0, written - 1) << " is " << atFarEnd << '\n';
        EXPECT_NE(atContact, atFarEnd)
            << "the tabulated value must vary with the initial separation r0";
    } else {
        ADD_FAILURE() << "expected at least two columns to compare separations";
    }

    // (b) much larger rate constant over the identical grid
    const double krLarge { 1000.0 };
    gsl_matrix* large = csm_alloc_sentinel_matrix(expectedCols + 2);
    create_survMatrix(large, bindRadius, Dtot, krLarge, comRMax, RStepSize, params);

    if (written > 0) {
        const double smallContact { gsl_matrix_get(small, 1, 0) };
        const double largeContact { gsl_matrix_get(large, 1, 0) };
        std::cerr << "  contact value with ka = " << krSmall << ": " << smallContact
                  << ", with ka = " << krLarge << ": " << largeContact << '\n';
        EXPECT_NE(smallContact, largeContact)
            << "changing the intrinsic association rate must change the tabulated value";

        // Both tables must still describe the same separation grid.
        for (size_t col = 0; col < written; ++col) {
            EXPECT_DOUBLE_EQ(gsl_matrix_get(small, 0, col), gsl_matrix_get(large, 0, col))
                << "the separation grid must not depend on the rate constant (column "
                << col << ')';
        }
    }

    gsl_matrix_free(small);
    gsl_matrix_free(large);
}

// -----------------------------------------------------------------------------
// Test 7: degenerate range - when comRMax is smaller than the binding radius by
//         more than one step, the loop body never runs and nothing is written.
//         (With comRMax + RStepSize < bindRadius the while-condition is false on
//         the very first evaluation.)
// -----------------------------------------------------------------------------
void test_csm_empty_range_writes_nothing()
{
    std::cerr << "\n[TEST] test_csm_empty_range_writes_nothing\n"
              << "  Source file: src/reactions/create_survMatrix.cpp\n"
              << "  Function:    create_survMatrix()\n"
              << "  Scenario:    comRMax + RStepSize < bindRadius, so the while-loop\n"
              << "               condition fails immediately.\n"
              << "  Pass:        the matrix is left completely untouched.\n";

    const double bindRadius { 5.0 };
    const double Dtot { 0.1 };
    const double kr { 10.0 };
    const double comRMax { 1.0 };  // far below the binding radius
    const double RStepSize { 0.05 };
    const Parameters params { csm_make_params(0.1) };

    const size_t expectedCols { csm_expected_columns(bindRadius, comRMax, RStepSize) };
    std::cerr << "  Expected number of filled columns: " << expectedCols << '\n';
    EXPECT_EQ(expectedCols, 0u) << "helper and production loop must agree that no step is taken";

    gsl_matrix* survMatrix = csm_alloc_sentinel_matrix(4);
    create_survMatrix(survMatrix, bindRadius, Dtot, kr, comRMax, RStepSize, params);

    for (size_t col = 0; col < survMatrix->size2; ++col) {
        EXPECT_DOUBLE_EQ(gsl_matrix_get(survMatrix, 0, col), kCsmSentinel)
            << "row 0, column " << col << " should be untouched for an empty range";
        EXPECT_DOUBLE_EQ(gsl_matrix_get(survMatrix, 1, col), kCsmSentinel)
            << "row 1, column " << col << " should be untouched for an empty range";
    }

    gsl_matrix_free(survMatrix);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each named helper is run inside its own TEST so a single
// failure never prevents the remaining cases from executing.
// -----------------------------------------------------------------------------
TEST(CreateSurvMatrix, Row0SeparationGrid) { test_csm_row0_separation_grid(); }
TEST(CreateSurvMatrix, Row1ValuesAreFinite) { test_csm_row1_values_are_finite(); }
TEST(CreateSurvMatrix, DoesNotWritePastLastColumn) { test_csm_does_not_write_past_last_column(); }
TEST(CreateSurvMatrix, ColumnCountFollowsStepAndRmax) { test_csm_column_count_follows_step_and_rmax(); }
TEST(CreateSurvMatrix, IsDeterministic) { test_csm_is_deterministic(); }
TEST(CreateSurvMatrix, ValuesDependOnInputs) { test_csm_values_depend_on_inputs(); }
TEST(CreateSurvMatrix, EmptyRangeWritesNothing) { test_csm_empty_range_writes_nothing(); }