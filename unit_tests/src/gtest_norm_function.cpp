// Unit test for norm_function from 2D_reaction_table_functions.cpp
// This test verifies the behavior of norm_function, which computes the
// integrand used in 2D reaction table calculations. The function depends on
// the modified Bessel function of the first kind (scaled), I0.

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_sf_bessel.h>
#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

// ---------------------------------------------------------------------------
// Helper functions for assertions. These print a clear message to stderr and
// exit with a non-zero code if the test fails, matching the required style.
// ---------------------------------------------------------------------------

// Check that 'actual' is close to 'expected' within a small tolerance.
void require_close(double actual, double expected, const std::string& label)
{
    const double tol = 1e-9;
    if (std::abs(actual - expected) > tol) {
        std::cerr << "FAILED [" << label << "]: expected " << expected
                  << " but got " << actual << " (diff "
                  << std::abs(actual - expected) << ")\n";
        std::exit(1);
    }
}

// Check that a boolean condition holds true.
void require_true(bool condition, const std::string& label)
{
    if (!condition) {
        std::cerr << "FAILED [" << label << "]: condition was false\n";
        std::exit(1);
    }
}

// ---------------------------------------------------------------------------
// Reference implementation of norm_function, used to compute the expected
// value independently so we can compare against the function under test.
// ---------------------------------------------------------------------------
double reference_norm(double x, const IntegrandParams& params)
{
    // temp is the argument used by the scaled Bessel function I0
    double temp = x * params.r0 / (2.0 * params.D * params.t);
    // temp2 is the leading prefactor x / (2 D t)
    double temp2 = x / (2.0 * params.D * params.t);
    // temp3 is the exponential term
    double temp3 = std::exp(temp - (params.r0 * params.r0 + (x * x))
                                       / (4.0 * params.t * params.D));
    // multiply prefactor, exponential, and the scaled Bessel function
    return temp2 * temp3 * gsl_sf_bessel_I0_scaled(temp);
}

// ---------------------------------------------------------------------------
// Test: norm_function returns the same value as our reference computation
// for a typical set of parameters.
// ---------------------------------------------------------------------------
void test_matches_reference()
{
    // Set up a realistic set of integrand parameters.
    IntegrandParams params;
    params.r0 = 1.0; // initial separation
    params.D = 0.5;  // diffusion coefficient
    params.t = 2.0;  // time

    // Test the function at several different x values.
    for (double x : { 0.1, 0.5, 1.0, 2.0, 5.0 }) {
        double actual = norm_function(x, &params);
        double expected = reference_norm(x, params);
        require_close(actual, expected,
            "norm_function matches reference at x=" + std::to_string(x));
    }
}

// ---------------------------------------------------------------------------
// Test: at x = 0 the prefactor (x / (2 D t)) is zero, so the integrand must
// be exactly zero regardless of the other parameters.
// ---------------------------------------------------------------------------
void test_zero_at_x_zero()
{
    IntegrandParams params;
    params.r0 = 1.0;
    params.D = 0.5;
    params.t = 2.0;

    double actual = norm_function(0.0, &params);
    require_close(actual, 0.0, "norm_function is zero at x=0");
}

// ---------------------------------------------------------------------------
// Test: the integrand should always be non-negative for positive parameters,
// since each factor (prefactor, exponential, scaled Bessel I0) is positive.
// ---------------------------------------------------------------------------
void test_non_negative()
{
    IntegrandParams params;
    params.r0 = 2.0;
    params.D = 1.0;
    params.t = 1.0;

    for (double x : { 0.1, 0.5, 1.0, 3.0, 10.0 }) {
        double value = norm_function(x, &params);
        require_true(value >= 0.0,
            "norm_function non-negative at x=" + std::to_string(x));
    }
}

// ---------------------------------------------------------------------------
// Test: when r0 = 0, the Bessel argument temp becomes 0 and I0_scaled(0) = 1.
// This lets us verify a simplified closed-form expression.
// ---------------------------------------------------------------------------
void test_r0_zero_simplification()
{
    IntegrandParams params;
    params.r0 = 0.0; // initial separation zero
    params.D = 0.5;
    params.t = 2.0;

    double x = 1.5;
    // With r0 = 0: temp = 0, I0_scaled(0) = 1, so the result simplifies to
    // (x / (2 D t)) * exp(-x^2 / (4 t D)).
    double expected = (x / (2.0 * params.D * params.t))
                      * std::exp(-(x * x) / (4.0 * params.t * params.D));
    double actual = norm_function(x, &params);
    require_close(actual, expected,
        "norm_function simplifies correctly when r0=0");
}

// ---------------------------------------------------------------------------
// GoogleTest wrappers so the test functions are picked up by the test runner.
// ---------------------------------------------------------------------------
TEST(NormFunctionTest, MatchesReference)     { test_matches_reference(); }
TEST(NormFunctionTest, ZeroAtXZero)          { test_zero_at_x_zero(); }
TEST(NormFunctionTest, NonNegative)          { test_non_negative(); }
TEST(NormFunctionTest, R0ZeroSimplification) { test_r0_zero_simplification(); }

