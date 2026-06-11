#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_sf_bessel.h>
#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

bool verbose_output_enabled()
{
    const char* value = std::getenv("NERDSS_TEST_VERBOSE");
    return value != nullptr && std::string(value) != "0";
}

void log_detail(const std::string& label, double actual, double expected)
{
    if (!verbose_output_enabled()) {
        return;
    }

    std::cout << label << ": actual=" << actual
              << " expected=" << expected
              << " diff=" << std::abs(actual - expected) << '\n';
}

double reference_norm(double x, const IntegrandParams& params)
{
    const double temp = x * params.r0 / (2.0 * params.D * params.t);
    const double prefactor = x / (2.0 * params.D * params.t);
    const double exponential = std::exp(
        temp - (params.r0 * params.r0 + (x * x)) / (4.0 * params.t * params.D));

    return prefactor * exponential * gsl_sf_bessel_I0_scaled(temp);
}

IntegrandParams make_params(double r0, double diffusion, double time)
{
    IntegrandParams params {};
    params.r0 = r0;
    params.D = diffusion;
    params.t = time;
    return params;
}

void expect_close(double actual, double expected, const std::string& label)
{
    log_detail(label, actual, expected);
    EXPECT_NEAR(actual, expected, 1.0e-9) << label;
}

} // namespace

TEST(NormFunctionTest, MatchesReferenceAcrossSamplePoints)
{
    const IntegrandParams params = make_params(1.0, 0.5, 2.0);

    for (double x : { 0.1, 0.5, 1.0, 2.0, 5.0 }) {
        expect_close(norm_function(x, const_cast<IntegrandParams*>(&params)),
            reference_norm(x, params),
            "norm_function at x=" + std::to_string(x));
    }
}

TEST(NormFunctionTest, IsZeroAtOrigin)
{
    IntegrandParams params = make_params(1.0, 0.5, 2.0);

    expect_close(norm_function(0.0, &params), 0.0, "norm_function at x=0");
}

TEST(NormFunctionTest, IsNonNegativeForPositiveParameters)
{
    IntegrandParams params = make_params(2.0, 1.0, 1.0);

    for (double x : { 0.1, 0.5, 1.0, 3.0, 10.0 }) {
        const double value = norm_function(x, &params);
        if (verbose_output_enabled()) {
            std::cout << "norm_function at x=" << x << ": value=" << value << '\n';
        }
        EXPECT_GE(value, 0.0) << "x=" << x;
    }
}

TEST(NormFunctionTest, SimplifiesWhenInitialSeparationIsZero)
{
    IntegrandParams params = make_params(0.0, 0.5, 2.0);
    const double x = 1.5;
    const double expected = (x / (2.0 * params.D * params.t))
        * std::exp(-(x * x) / (4.0 * params.t * params.D));

    expect_close(norm_function(x, &params), expected, "r0=0 simplification");
}
