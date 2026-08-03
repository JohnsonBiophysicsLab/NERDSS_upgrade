/*! \file test_integrator.cpp
 *
 * ### Unit test for src/reactions/integrator.cpp
 *
 * The function under test is:
 *
 *     double integrator(gsl_function F, IntegrandParams params,
 *                       gsl_integration_workspace* w, double r0, double bindrad,
 *                       double Dtot, double kr, double deltat, char* funcID,
 *                       double (*f)(double, void*));
 *
 * It is a thin (but not trivial) wrapper around the GSL adaptive quadrature
 * routines used when the 2D reaction lookup tables are generated:
 *
 *   1. It first tries `gsl_integration_qagiu` on the semi-infinite domain
 *      [0, inf) with absolute/relative tolerances of 1e-7.
 *   2. If that fails it retries with a weaker (1e-6) criterion.
 *   3. If that also fails it probes the raw integrand pointer `f` to find a
 *      finite upper bound `umax` and falls back to `gsl_integration_qag`.
 *
 * Since the numerical answer is the only observable output, the tests below
 * feed the integrator analytically-known, well-behaved (rapidly decaying)
 * integrands and check the returned value against the exact result.  A
 * "probe counter" integrand is also used to verify that the fallback branch is
 * *not* entered when the primary quadrature succeeds.
 *
 * Notes on safety:
 *   - `integrator()` passes a subinterval limit of 1e6 to GSL, so the workspace
 *     handed to it must be allocated with at least that many subintervals,
 *     otherwise GSL raises a fatal error.
 *   - The GSL error handler is switched off for the duration of each test so
 *     that a non-convergent integration returns a status code (which the code
 *     under test is written to handle) instead of aborting the whole test
 *     binary.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

// -----------------------------------------------------------------------------
// Local helpers: analytic integrands and small utilities.
// All names are prefixed with `integ_` so they cannot collide with symbols
// coming from other translation units in the combined test binary.
// -----------------------------------------------------------------------------
namespace {

//! Number of subintervals in the workspace. Must be >= the 1e6 limit that
//! integrator() hands to the GSL routines.
constexpr size_t kIntegWorkspaceSize = 1000000;

//! Counter incremented every time the "raw" integrand pointer (the `f`
//! argument of integrator()) is invoked. Only the fallback branch calls it.
int integ_fallbackProbeCalls = 0;

/*! \brief f(x) = exp(-x); integral over [0, inf) == 1. */
double integ_exp_decay(double x, void* /*p*/) { return std::exp(-x); }

/*! \brief f(x) = exp(-x^2); integral over [0, inf) == sqrt(pi)/2. */
double integ_gaussian(double x, void* /*p*/) { return std::exp(-x * x); }

/*! \brief f(x) = x*exp(-x); integral over [0, inf) == 1 (Gamma(2)). */
double integ_x_exp_decay(double x, void* /*p*/) { return x * std::exp(-x); }

/*! \brief f(x) = exp(-x)*cos(x); integral over [0, inf) == 1/2.
 *  An oscillatory-but-decaying integrand, similar in spirit to the Bessel
 *  based integrands actually used by the 2D tables. */
double integ_damped_cos(double x, void* /*p*/) { return std::exp(-x) * std::cos(x); }

/*! \brief f(x) = 0 everywhere; integral == 0 exactly. */
double integ_zero(double /*x*/, void* /*p*/) { return 0.0; }

/*! \brief f(u) = u*exp(-D*t*u^2) using the IntegrandParams handed through
 *  gsl_function::params; integral over [0, inf) == 1/(2*D*t).
 *
 *  This verifies that integrator() forwards the user parameter block untouched
 *  to the GSL quadrature machinery. */
double integ_param_scaled(double u, void* p)
{
    const IntegrandParams* ip = static_cast<const IntegrandParams*>(p);
    return u * std::exp(-ip->D * ip->t * u * u);
}

/*! \brief Decaying probe used as the `f` argument; counts its invocations.
 *
 *  Because it decays like exp(-x), the `while (|f(umax)| > 1e-10)` loop inside
 *  the fallback branch of integrator() would still terminate if it were ever
 *  reached. */
double integ_fallback_probe(double x, void* /*p*/)
{
    ++integ_fallbackProbeCalls;
    return std::exp(-x);
}

/*! \brief RAII-ish helper bundling the GSL workspace and error-handler state.
 *
 *  Constructing it allocates a large workspace and disables the fatal GSL
 *  error handler; destroying it restores the handler and frees the workspace.
 */
struct IntegTestEnv {
    gsl_integration_workspace* w { nullptr };
    gsl_error_handler_t* oldHandler { nullptr };

    IntegTestEnv()
    {
        oldHandler = gsl_set_error_handler_off();
        w = gsl_integration_workspace_alloc(kIntegWorkspaceSize);
    }
    ~IntegTestEnv()
    {
        if (w != nullptr)
            gsl_integration_workspace_free(w);
        gsl_set_error_handler(oldHandler);
    }
};

/*! \brief Builds a default-ish IntegrandParams block.
 *
 *  The values mimic a plausible 2D-table call (sigma, diffusion constant,
 *  association rate, time step), even though the analytic integrands used here
 *  mostly ignore them. */
IntegrandParams integ_make_params(double a, double D, double k, double t)
{
    IntegrandParams p(a, D, k, t);
    p.r0 = 1.5 * a;
    p.r = 1.5 * a;
    p.rho = 0.01;
    return p;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: exp(-x) over [0, inf) must integrate to exactly 1.
// -----------------------------------------------------------------------------
void test_integrator_exponential_decay()
{
    std::cerr << "\n[TEST] test_integrator_exponential_decay\n"
              << "  Source file:   src/reactions/integrator.cpp\n"
              << "  Function:      integrator()\n"
              << "  Integrand:     f(x) = exp(-x) on [0, inf)\n"
              << "  Pass criteria: returned value == 1.0 within 1e-6\n";

    IntegTestEnv env;
    ASSERT_NE(env.w, nullptr) << "GSL workspace allocation failed";

    // Build the gsl_function that GSL will evaluate.
    IntegrandParams params { integ_make_params(1.0, 1.0, 10.0, 0.1) };
    gsl_function F;
    F.function = &integ_exp_decay;
    F.params = &params;

    char funcID[] = "exp_decay";

    std::cerr << "  Calling integrator()...\n";
    const double result = integrator(F, params, env.w,
        /*r0=*/params.r0, /*bindrad=*/params.a, /*Dtot=*/params.D,
        /*kr=*/params.k, /*deltat=*/params.t, funcID, &integ_exp_decay);

    std::cerr << "  Result = " << result << " (expected 1)\n";
    EXPECT_NEAR(result, 1.0, 1e-6) << "integral of exp(-x) on [0,inf) must be 1";
    EXPECT_TRUE(std::isfinite(result)) << "result must be a finite number";
}

// -----------------------------------------------------------------------------
// Test 2: exp(-x^2) over [0, inf) must integrate to sqrt(pi)/2.
// -----------------------------------------------------------------------------
void test_integrator_gaussian_tail()
{
    std::cerr << "\n[TEST] test_integrator_gaussian_tail\n"
              << "  Source file:   src/reactions/integrator.cpp\n"
              << "  Function:      integrator()\n"
              << "  Integrand:     f(x) = exp(-x^2) on [0, inf)\n"
              << "  Pass criteria: returned value == sqrt(pi)/2 within 1e-6\n";

    IntegTestEnv env;
    ASSERT_NE(env.w, nullptr) << "GSL workspace allocation failed";

    IntegrandParams params { integ_make_params(1.0, 5.0, 1.0, 1.0) };
    gsl_function F;
    F.function = &integ_gaussian;
    F.params = &params;

    char funcID[] = "gaussian";

    const double expected = 0.5 * std::sqrt(M_PI);

    std::cerr << "  Calling integrator()...\n";
    const double result = integrator(F, params, env.w, params.r0, params.a,
        params.D, params.k, params.t, funcID, &integ_gaussian);

    std::cerr << "  Result = " << result << " (expected " << expected << ")\n";
    EXPECT_NEAR(result, expected, 1e-6)
        << "integral of exp(-x^2) on [0,inf) must be sqrt(pi)/2";
}

// -----------------------------------------------------------------------------
// Test 3: x*exp(-x) over [0, inf) must integrate to 1 (Gamma(2) = 1!).
//         This integrand rises before decaying, exercising the adaptive
//         subdivision more than a pure exponential.
// -----------------------------------------------------------------------------
void test_integrator_peaked_integrand()
{
    std::cerr << "\n[TEST] test_integrator_peaked_integrand\n"
              << "  Source file:   src/reactions/integrator.cpp\n"
              << "  Function:      integrator()\n"
              << "  Integrand:     f(x) = x*exp(-x) on [0, inf) (peaked at x=1)\n"
              << "  Pass criteria: returned value == 1.0 within 1e-6\n";

    IntegTestEnv env;
    ASSERT_NE(env.w, nullptr) << "GSL workspace allocation failed";

    IntegrandParams params { integ_make_params(2.0, 1.0, 100.0, 0.01) };
    gsl_function F;
    F.function = &integ_x_exp_decay;
    F.params = &params;

    char funcID[] = "x_exp_decay";

    std::cerr << "  Calling integrator()...\n";
    const double result = integrator(F, params, env.w, params.r0, params.a,
        params.D, params.k, params.t, funcID, &integ_x_exp_decay);

    std::cerr << "  Result = " << result << " (expected 1)\n";
    EXPECT_NEAR(result, 1.0, 1e-6) << "integral of x*exp(-x) on [0,inf) must be 1";
}

// -----------------------------------------------------------------------------
// Test 4: exp(-x)*cos(x) over [0, inf) must integrate to 1/2. Oscillatory
//         integrands are the ones that historically caused the fallback branch
//         to be needed, so this checks a mildly oscillatory case still works.
// -----------------------------------------------------------------------------
void test_integrator_oscillatory_integrand()
{
    std::cerr << "\n[TEST] test_integrator_oscillatory_integrand\n"
              << "  Source file:   src/reactions/integrator.cpp\n"
              << "  Function:      integrator()\n"
              << "  Integrand:     f(x) = exp(-x)*cos(x) on [0, inf)\n"
              << "  Pass criteria: returned value == 0.5 within 1e-6\n";

    IntegTestEnv env;
    ASSERT_NE(env.w, nullptr) << "GSL workspace allocation failed";

    IntegrandParams params { integ_make_params(1.0, 1.0, 1.0, 1.0) };
    gsl_function F;
    F.function = &integ_damped_cos;
    F.params = &params;

    char funcID[] = "damped_cos";

    std::cerr << "  Calling integrator()...\n";
    const double result = integrator(F, params, env.w, params.r0, params.a,
        params.D, params.k, params.t, funcID, &integ_damped_cos);

    std::cerr << "  Result = " << result << " (expected 0.5)\n";
    EXPECT_NEAR(result, 0.5, 1e-6)
        << "integral of exp(-x)*cos(x) on [0,inf) must be 1/2";
}

// -----------------------------------------------------------------------------
// Test 5: A zero integrand must give exactly zero. This also confirms the
//         function returns a value (rather than, e.g., leaving `result`
//         uninitialised) for the degenerate case.
// -----------------------------------------------------------------------------
void test_integrator_zero_integrand()
{
    std::cerr << "\n[TEST] test_integrator_zero_integrand\n"
              << "  Source file:   src/reactions/integrator.cpp\n"
              << "  Function:      integrator()\n"
              << "  Integrand:     f(x) = 0 on [0, inf)\n"
              << "  Pass criteria: returned value == 0.0 exactly\n";

    IntegTestEnv env;
    ASSERT_NE(env.w, nullptr) << "GSL workspace allocation failed";

    IntegrandParams params { integ_make_params(1.0, 1.0, 1.0, 1.0) };
    gsl_function F;
    F.function = &integ_zero;
    F.params = &params;

    char funcID[] = "zero";

    std::cerr << "  Calling integrator()...\n";
    const double result = integrator(F, params, env.w, params.r0, params.a,
        params.D, params.k, params.t, funcID, &integ_zero);

    std::cerr << "  Result = " << result << " (expected 0)\n";
    EXPECT_DOUBLE_EQ(result, 0.0) << "integral of the zero function must be 0";
}

// -----------------------------------------------------------------------------
// Test 6: The IntegrandParams block attached to gsl_function::params must be
//         forwarded unmodified, so an integrand whose decay rate is taken from
//         the parameters produces the parameter-dependent analytic answer:
//
//              int_0^inf u*exp(-D*t*u^2) du = 1/(2*D*t)
//
//         Two different (D, t) pairs are checked with the same workspace.
// -----------------------------------------------------------------------------
void test_integrator_uses_params_block()
{
    std::cerr << "\n[TEST] test_integrator_uses_params_block\n"
              << "  Source file:   src/reactions/integrator.cpp\n"
              << "  Function:      integrator()\n"
              << "  Integrand:     f(u) = u*exp(-D*t*u^2), D and t read from\n"
              << "                 the IntegrandParams block\n"
              << "  Pass criteria: returned value == 1/(2*D*t) within 1e-6 for\n"
              << "                 two different parameter sets\n";

    IntegTestEnv env;
    ASSERT_NE(env.w, nullptr) << "GSL workspace allocation failed";

    char funcID[] = "param_scaled";

    // --- First parameter set: D = 1, t = 2  => expected 0.25 -----------------
    {
        IntegrandParams params { integ_make_params(1.0, 1.0, 10.0, 2.0) };
        gsl_function F;
        F.function = &integ_param_scaled;
        F.params = &params;

        const double expected = 1.0 / (2.0 * params.D * params.t);
        std::cerr << "  Calling integrator() with D=" << params.D
                  << ", t=" << params.t << " (expected " << expected << ")\n";
        const double result = integrator(F, params, env.w, params.r0, params.a,
            params.D, params.k, params.t, funcID, &integ_param_scaled);
        std::cerr << "  Result = " << result << '\n';
        EXPECT_NEAR(result, expected, 1e-6)
            << "parameters D/t must be honoured by the integration";
    }

    // --- Second parameter set: D = 4, t = 0.5 => expected 0.25 as well, but
    //     via completely different numbers, so a hard-coded answer would fail
    //     the first case. Use D = 4, t = 0.25 => expected 0.5 -----------------
    {
        IntegrandParams params { integ_make_params(1.0, 4.0, 10.0, 0.25) };
        gsl_function F;
        F.function = &integ_param_scaled;
        F.params = &params;

        const double expected = 1.0 / (2.0 * params.D * params.t);
        std::cerr << "  Calling integrator() with D=" << params.D
                  << ", t=" << params.t << " (expected " << expected << ")\n";
        const double result = integrator(F, params, env.w, params.r0, params.a,
            params.D, params.k, params.t, funcID, &integ_param_scaled);
        std::cerr << "  Result = " << result << '\n';
        EXPECT_NEAR(result, expected, 1e-6)
            << "a second parameter set must give the matching analytic answer";
    }
}

// -----------------------------------------------------------------------------
// Test 7: The raw integrand pointer `f` is only used inside the truncation
//         fallback. For a smooth, rapidly decaying integrand the primary
//         `gsl_integration_qagiu` call succeeds, so `f` must never be invoked.
// -----------------------------------------------------------------------------
void test_integrator_fallback_not_triggered()
{
    std::cerr << "\n[TEST] test_integrator_fallback_not_triggered\n"
              << "  Source file:   src/reactions/integrator.cpp\n"
              << "  Function:      integrator()\n"
              << "  Scenario:      well-behaved integrand exp(-x); the `f`\n"
              << "                 pointer is a counting probe.\n"
              << "  Pass criteria: probe call count stays 0 (fallback branch\n"
              << "                 never entered) and result == 1.0\n";

    IntegTestEnv env;
    ASSERT_NE(env.w, nullptr) << "GSL workspace allocation failed";

    integ_fallbackProbeCalls = 0;

    IntegrandParams params { integ_make_params(1.0, 1.0, 1.0, 1.0) };
    gsl_function F;
    F.function = &integ_exp_decay; // GSL integrates this
    F.params = &params;

    char funcID[] = "probe";

    std::cerr << "  Calling integrator() with probe as the raw `f` pointer...\n";
    const double result = integrator(F, params, env.w, params.r0, params.a,
        params.D, params.k, params.t, funcID, &integ_fallback_probe);

    std::cerr << "  Result = " << result
              << ", probe invocations = " << integ_fallbackProbeCalls << '\n';
    EXPECT_NEAR(result, 1.0, 1e-6) << "primary quadrature should give 1.0";
    EXPECT_EQ(integ_fallbackProbeCalls, 0)
        << "the truncation fallback must not run for a convergent integral";
}

// -----------------------------------------------------------------------------
// Test 8: Repeated calls with the same workspace must be reproducible and
//         independent of one another (the workspace is reset internally by GSL
//         on every call). This mirrors how the 2D table generators reuse a
//         single workspace for thousands of integrals.
// -----------------------------------------------------------------------------
void test_integrator_workspace_reuse()
{
    std::cerr << "\n[TEST] test_integrator_workspace_reuse\n"
              << "  Source file:   src/reactions/integrator.cpp\n"
              << "  Function:      integrator()\n"
              << "  Scenario:      three integrals evaluated back-to-back with\n"
              << "                 the same gsl_integration_workspace.\n"
              << "  Pass criteria: each integral matches its analytic value and\n"
              << "                 repeating the first one reproduces it bitwise\n";

    IntegTestEnv env;
    ASSERT_NE(env.w, nullptr) << "GSL workspace allocation failed";

    IntegrandParams params { integ_make_params(1.0, 1.0, 1.0, 1.0) };
    char funcID[] = "reuse";

    // First integral: exp(-x) -> 1
    gsl_function F1;
    F1.function = &integ_exp_decay;
    F1.params = &params;
    const double first = integrator(F1, params, env.w, params.r0, params.a,
        params.D, params.k, params.t, funcID, &integ_exp_decay);
    std::cerr << "  Call 1 (exp(-x))      = " << first << " (expected 1)\n";
    EXPECT_NEAR(first, 1.0, 1e-6) << "first integral must be 1";

    // Second integral, same workspace: exp(-x)*cos(x) -> 0.5
    gsl_function F2;
    F2.function = &integ_damped_cos;
    F2.params = &params;
    const double second = integrator(F2, params, env.w, params.r0, params.a,
        params.D, params.k, params.t, funcID, &integ_damped_cos);
    std::cerr << "  Call 2 (exp(-x)cos x) = " << second << " (expected 0.5)\n";
    EXPECT_NEAR(second, 0.5, 1e-6)
        << "second integral must be unaffected by workspace reuse";

    // Third call repeats the first integral: must reproduce it exactly.
    const double repeat = integrator(F1, params, env.w, params.r0, params.a,
        params.D, params.k, params.t, funcID, &integ_exp_decay);
    std::cerr << "  Call 3 (repeat of 1)  = " << repeat << '\n';
    EXPECT_DOUBLE_EQ(repeat, first)
        << "repeating an integral with the reused workspace must be deterministic";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper runs inside its own TEST so
// that a failure in one does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(IntegratorTest, ExponentialDecay) { test_integrator_exponential_decay(); }
TEST(IntegratorTest, GaussianTail) { test_integrator_gaussian_tail(); }
TEST(IntegratorTest, PeakedIntegrand) { test_integrator_peaked_integrand(); }
TEST(IntegratorTest, OscillatoryIntegrand) { test_integrator_oscillatory_integrand(); }
TEST(IntegratorTest, ZeroIntegrand) { test_integrator_zero_integrand(); }
TEST(IntegratorTest, UsesParamsBlock) { test_integrator_uses_params_block(); }
TEST(IntegratorTest, FallbackNotTriggered) { test_integrator_fallback_not_triggered(); }
TEST(IntegratorTest, WorkspaceReuse) { test_integrator_workspace_reuse(); }