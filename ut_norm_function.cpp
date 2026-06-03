#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include <gsl/gsl_sf_bessel.h>

double norm_function(double x, void* p);

static void require_close(double actual, double expected, const char* label)
{
    const double tol = 1e-9;
    double diff = std::fabs(actual - expected);
    double scale = std::max(1.0, std::fabs(expected));
    if (diff > tol * scale)
    {
        std::cerr << "FAILED: " << label << " : actual=" << actual
                  << " expected=" << expected << " diff=" << diff << std::endl;
        std::exit(1);
    }
}

static void require_true(bool condition, const char* label)
{
    if (!condition)
    {
        std::cerr << "FAILED: " << label << std::endl;
        std::exit(1);
    }
}

// Reference implementation matching the function under test.
static double reference_norm_function(double x, double r0, double D, double t)
{
    double temp = x * r0 / (2.0 * D * t);
    double temp2 = x / (2.0 * D * t);
    double temp3 = std::exp(temp - (r0 * r0 + (x * x)) / (4.0 * t * D));
    return temp2 * temp3 * gsl_sf_bessel_I0_scaled(temp);
}

static IntegrandParams make_params(double r0, double D, double t)
{
    IntegrandParams params {};
    params.r0 = r0;
    params.D = D;
    params.t = t;
    return params;
}

void test_norm_function_basic()
{
    IntegrandParams params = make_params(1.0, 1.0, 1.0);
    double x = 2.0;
    double actual = norm_function(x, &params);
    double expected = reference_norm_function(x, params.r0, params.D, params.t);
    require_close(actual, expected, "norm_function basic value");
}

void test_norm_function_multiple_points()
{
    IntegrandParams params = make_params(2.0, 0.5, 1.5);
    double xs[] = { 0.5, 1.0, 2.5, 5.0, 10.0 };
    for (double x : xs)
    {
        double actual = norm_function(x, &params);
        double expected = reference_norm_function(x, params.r0, params.D, params.t);
        require_close(actual, expected, "norm_function multiple points");
    }
}

void test_norm_function_zero_at_x_zero()
{
    IntegrandParams params = make_params(1.0, 1.0, 1.0);
    double actual = norm_function(0.0, &params);
    // x=0 makes temp2 = 0, so result must be exactly 0.
    require_close(actual, 0.0, "norm_function at x=0");
}

void test_norm_function_positive()
{
    IntegrandParams params = make_params(1.0, 2.0, 0.75);
    double xs[] = { 0.1, 1.0, 3.0, 7.0 };
    for (double x : xs)
    {
        double actual = norm_function(x, &params);
        require_true(actual > 0.0, "norm_function positive for positive x");
    }
}

void test_norm_function_finite()
{
    IntegrandParams params = make_params(3.0, 1.0, 2.0);
    double xs[] = { 0.0, 0.5, 5.0, 50.0 };
    for (double x : xs)
    {
        double actual = norm_function(x, &params);
        require_true(std::isfinite(actual), "norm_function returns finite value");
    }
}

int main()
{
    test_norm_function_basic();
    test_norm_function_multiple_points();
    test_norm_function_zero_at_x_zero();
    test_norm_function_positive();
    test_norm_function_finite();

    std::cout << "All tests passed." << std::endl;
    return 0;
}
