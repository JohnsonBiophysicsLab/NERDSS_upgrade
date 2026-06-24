// Unit test for calc_pirr() in:
//   reactions/bimolecular/2D_reaction_table_functions.cpp
//
// calc_pirr() performs a 2D lookup / interpolation over a pair of GSL
// matrices (a "pir" probability matrix and a "surv" coordinate matrix).
// It maps an arbitrary coordinate pair (r, r0) to a probability value
// using indices computed from RStepSize and a, then either returns the
// raw matrix entry (when both indices coincide) or an inverse-distance
// weighted interpolation of the four surrounding grid points.
//
// We test:
//   1. The "same index" branch (indexr0 == indexr) returns the direct entry.
//   2. Negative indices are clamped to 0.
//   3. The interpolation branch produces the expected inverse-distance
//      weighted average.
//
// Build requires GSL and gtest.

#include <gtest/gtest.h>
#include <gsl/gsl_matrix.h>
#include <cmath>
#include <cstdlib>
#include <iostream>

// Forward declaration of the function under test (defined in the .cpp).
double calc_pirr(gsl_matrix* pirMatrix, gsl_matrix* survMatrix, double RStepSize,
                 double r, double r0, double a);

// ---------------------------------------------------------------------------
// Helper assertion functions: print verbose diagnostics to stderr and abort
// the whole program with exit(1) if a check fails.
// ---------------------------------------------------------------------------

// Compare two doubles within a small tolerance.
static void require_close(double actual, double expected, const char* label)
{
    const double tol = 1e-9;
    std::cerr << "    [CHECK] " << label
              << " : actual=" << actual
              << " expected=" << expected
              << " (tol=" << tol << ")\n";
    if (std::fabs(actual - expected) > tol) {
        std::cerr << "    [FAIL ] " << label
                  << " : difference " << std::fabs(actual - expected)
                  << " exceeds tolerance " << tol << "\n";
        std::exit(1);
    }
    std::cerr << "    [ OK  ] " << label << "\n";
}

// Verify a boolean condition.
static void require_true(bool condition, const char* label)
{
    std::cerr << "    [CHECK] " << label << " : "
              << (condition ? "true" : "false") << "\n";
    if (!condition) {
        std::cerr << "    [FAIL ] " << label << "\n";
        std::exit(1);
    }
    std::cerr << "    [ OK  ] " << label << "\n";
}

// ---------------------------------------------------------------------------
// Small utilities to build the GSL matrices used as lookup tables.
// ---------------------------------------------------------------------------

// Allocate a rows x cols matrix and zero-fill it.
static gsl_matrix* make_matrix(size_t rows, size_t cols)
{
    gsl_matrix* m = gsl_matrix_alloc(rows, cols);
    gsl_matrix_set_zero(m);
    return m;
}

// ---------------------------------------------------------------------------
// Test 1: the "same index" branch.
//
// If r and r0 map to the same grid index, calc_pirr should return exactly the
// matrix entry pirMatrix[index][index] without doing any interpolation.
// ---------------------------------------------------------------------------
void test_same_index_returns_direct_entry()
{
    std::cerr << "\n=== test_same_index_returns_direct_entry ===\n";
    std::cerr << "Testing calc_pirr() branch where indexr0 == indexr.\n";

    const double RStepSize = 1.0;
    const double a = 0.0;

    // Build a 3x3 pir matrix with recognizable entries.
    gsl_matrix* pir = make_matrix(3, 3);
    gsl_matrix_set(pir, 0, 0, 0.11);
    gsl_matrix_set(pir, 1, 1, 0.55); // entry we expect to be returned
    gsl_matrix_set(pir, 2, 2, 0.99);

    // surv matrix is unused on this branch but must exist; one row suffices.
    gsl_matrix* surv = make_matrix(1, 3);

    // Choose r and r0 such that both map to index 1:
    //   index = floor((value - a) / RStepSize) = floor(1.5) = 1
    const double r = 1.5;
    const double r0 = 1.5;

    std::cerr << "  Inputs: r=" << r << " r0=" << r0
              << " RStepSize=" << RStepSize << " a=" << a << "\n";

    double result = calc_pirr(pir, surv, RStepSize, r, r0, a);
    require_close(result, 0.55, "same-index returns pir[1][1]");

    gsl_matrix_free(pir);
    gsl_matrix_free(surv);
}

// ---------------------------------------------------------------------------
// Test 2: negative index clamping.
//
// When (r - a) or (r0 - a) is negative, the computed index is negative and
// must be clamped to 0. With both clamped to 0 the same-index branch is taken
// and pir[0][0] is returned.
// ---------------------------------------------------------------------------
void test_negative_indices_are_clamped()
{
    std::cerr << "\n=== test_negative_indices_are_clamped ===\n";
    std::cerr << "Testing calc_pirr() clamps negative indices to 0.\n";

    const double RStepSize = 1.0;
    const double a = 5.0; // makes (r - a) negative for small r

    gsl_matrix* pir = make_matrix(3, 3);
    gsl_matrix_set(pir, 0, 0, 0.42); // expected return after clamping

    gsl_matrix* surv = make_matrix(1, 3);

    // r and r0 below 'a' -> indices floor(negative) -> clamped to 0.
    const double r = 1.0;
    const double r0 = 1.0;

    std::cerr << "  Inputs: r=" << r << " r0=" << r0
              << " RStepSize=" << RStepSize << " a=" << a << "\n";

    double result = calc_pirr(pir, surv, RStepSize, r, r0, a);
    require_close(result, 0.42, "clamped negative indices return pir[0][0]");

    gsl_matrix_free(pir);
    gsl_matrix_free(surv);
}

// ---------------------------------------------------------------------------
// Test 3: inverse-distance interpolation branch.
//
// When indexr0 != indexr, calc_pirr interpolates between four grid points
// A, B, C, D using inverse-distance weighting. We construct a controlled
// scenario and compute the expected value with the same formula the function
// uses, then confirm they match.
// ---------------------------------------------------------------------------
void test_interpolation_branch()
{
    std::cerr << "\n=== test_interpolation_branch ===\n";
    std::cerr << "Testing calc_pirr() inverse-distance interpolation branch.\n";

    const double RStepSize = 1.0;
    const double a = 0.0;

    // r0 -> index 0, r -> index 1  (so indexr0 != indexr triggers interp).
    const double r0 = 0.5; // floor(0.5) = 0
    const double r = 1.5;  // floor(1.5) = 1
    const int indexr0 = 0;
    const int indexr = 1;

    // pir matrix values at the four corners used by the algorithm.
    gsl_matrix* pir = make_matrix(3, 3);
    gsl_matrix_set(pir, indexr0,     indexr,     0.10); // valA
    gsl_matrix_set(pir, indexr0,     indexr + 1, 0.20); // valB
    gsl_matrix_set(pir, indexr0 + 1, indexr + 1, 0.30); // valC
    gsl_matrix_set(pir, indexr0 + 1, indexr,     0.40); // valD

    // surv matrix supplies the grid coordinates (row 0 only).
    gsl_matrix* surv = make_matrix(1, 3);
    gsl_matrix_set(surv, 0, indexr0,     0.0); // used for *r0 coords A,B
    gsl_matrix_set(surv, 0, indexr,      1.0); // used for *r coords A,D
    gsl_matrix_set(surv, 0, indexr + 1,  2.0); // used for *r coords B,C and r0 coords C,D

    // --- Manually reproduce the function's expected computation ---
    double Ar0 = gsl_matrix_get(surv, 0, indexr0);
    double Ar  = gsl_matrix_get(surv, 0, indexr);
    double valA = gsl_matrix_get(pir, indexr0, indexr);
    double distA = std::sqrt(std::pow(Ar0 - r0, 2) + std::pow(Ar - r, 2));

    double Br0 = gsl_matrix_get(surv, 0, indexr0);
    double Br  = gsl_matrix_get(surv, 0, indexr + 1);
    double valB = gsl_matrix_get(pir, indexr0, indexr + 1);
    double distB = std::sqrt(std::pow(Br0 - r0, 2) + std::pow(Br - r, 2));

    double Cr0 = gsl_matrix_get(surv, 0, indexr0 + 1);
    double Cr  = gsl_matrix_get(surv, 0, indexr + 1);
    double valC = gsl_matrix_get(pir, indexr0 + 1, indexr + 1);
    double distC = std::sqrt(std::pow(Cr0 - r0, 2) + std::pow(Cr - r, 2));

    double Dr0 = gsl_matrix_get(surv, 0, indexr0 + 1);
    double Dr  = gsl_matrix_get(surv, 0, indexr);
    double valD = gsl_matrix_get(pir, indexr0 + 1, indexr);
    double distD = std::sqrt(std::pow(Dr0 - r0, 2) + std::pow(Dr - r, 2));

    double expected =
        ((valA / distA) + (valB / distB) + (valC / distC) + (valD / distD))
        / ((1 / distA) + (1 / distB) + (1 / distC) + (1 / distD));

    std::cerr << "  Inputs: r=" << r << " r0=" << r0
              << " RStepSize=" << RStepSize << " a=" << a << "\n";
    std::cerr << "  Independently computed expected interpolation="
              << expected << "\n";

    // Sanity: distances must be positive so we never divide by zero.
    require_true(distA > 0 && distB > 0 && distC > 0 && distD > 0,
                 "all corner distances are positive (no div-by-zero)");

    double result = calc_pirr(pir, surv, RStepSize, r, r0, a);
    require_close(result, expected, "interpolated probability matches formula");

    gsl_matrix_free(pir);
    gsl_matrix_free(surv);
}

// ---------------------------------------------------------------------------
// GoogleTest wrappers so the suite integrates with gtest while still using
// our verbose require_* helpers.
// ---------------------------------------------------------------------------
TEST(CalcPirrTest, SameIndexReturnsDirectEntry) { test_same_index_returns_direct_entry(); }
TEST(CalcPirrTest, NegativeIndicesAreClamped)   { test_negative_indices_are_clamped(); }
TEST(CalcPirrTest, InterpolationBranch)         { test_interpolation_branch(); }

// ---------------------------------------------------------------------------
// main: announce the file/function under test, then run all gtest cases.
// A test passes if every require_close/require_true succeeds