/*! \file test_pir_function.cpp
 *
 * ### Unit test for src/reactions/pir_function.cpp
 *
 * The file under test provides a single free function:
 *
 *     double pir_function(double x, void* p);
 *
 * It is used as a GSL integrand while building the 2D reaction lookup
 * tables (see create_pirMatrix()).  The integrand is
 *
 *     f(x) = x * exp(-D t x^2) * P(r) * T(r0) / (2 pi)
 *
 * where P and T are built from the cylindrical Bessel functions
 * j0/j1/y0/y1 evaluated at the binding radius `a`, and where the
 * radiation (finite k) or absorbing (k == infinity) boundary condition
 * selects which pair of coefficients (alp, bet) is used.
 *
 * Because the function is a pure mathematical expression of the values in
 * the IntegrandParams struct, the strategy here is:
 *
 *   1. Re-implement the documented formula independently in the test file
 *      and require exact agreement (within a tight tolerance).
 *   2. Verify the two code branches (finite k vs. k == infinity) are both
 *      reachable and produce the expected, different results.
 *   3. Verify analytic invariants that must hold regardless of branch:
 *        - symmetry under swapping r and r0 (P and T just trade places),
 *        - the exp(-D t x^2) time-decay factor,
 *        - f(x) == x * P * T / (2 pi) when t == 0,
 *        - k -> very large converges to the absorbing-boundary branch,
 *        - the input parameter struct is not modified (read-only access).
 *
 * All assertions are non-fatal EXPECT_* so every case runs even on failure,
 * and each step prints what is being tested to stderr.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {

//------------------------------------------------------------------------------
// Independent reference implementation of the integrand.
//
// This mirrors the documented mathematics of pir_function() but is written
// from scratch here so that a regression in the production code (a swapped
// sign, a dropped 2*pi, a wrong Bessel order, ...) shows up as a mismatch.
//------------------------------------------------------------------------------
double pirfunc_reference(double x, const IntegrandParams& p)
{
    const double h = 2.0 * M_PI * p.a * p.D;

    double alp = 0.0;
    double bet = 0.0;
    if (p.k < std::numeric_limits<double>::infinity()) {
        // Radiation (partially absorbing) boundary condition.
        alp = h * x * j1(x * p.a) + p.k * j0(x * p.a);
        bet = h * x * y1(x * p.a) + p.k * y0(x * p.a);
    } else {
        // Fully absorbing boundary condition.
        alp = j0(x * p.a);
        bet = y0(x * p.a);
    }

    const double tet = std::sqrt(alp * alp + bet * bet);

    const double P = (j0(x * p.r) * bet - y0(x * p.r) * alp) / tet;
    const double T = (j0(x * p.r0) * bet - y0(x * p.r0) * alp) / tet;

    return x * std::exp(-p.D * p.t * x * x) * P * T / (2.0 * M_PI);
}

/*! \brief Build a physically reasonable set of integrand parameters.
 *
 * Units follow the rest of NERDSS: lengths in nm, D in nm^2/us, t in us.
 */
IntegrandParams pirfunc_make_params(double k)
{
    IntegrandParams params;
    params.a = 1.0; // binding radius (sigma), nm
    params.D = 0.5; // total diffusion constant, nm^2/us
    params.k = k; // association rate constant
    params.r0 = 2.0; // initial separation, nm
    params.r = 3.0; // current separation, nm
    params.t = 1.0; // elapsed time, us
    params.rho = 0.01; // lipid density (unused by pir_function)
    return params;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: finite k (radiation boundary) must match the reference formula for a
//         sweep of integration variables x.
// -----------------------------------------------------------------------------
void pirfunc_test_finite_k_matches_reference()
{
    std::cerr << "\n[TEST] pirfunc_test_finite_k_matches_reference\n"
              << "  Source file:   src/reactions/pir_function.cpp\n"
              << "  Function:      pir_function(double x, void* p)\n"
              << "  Scenario:      finite k -> radiation boundary branch.\n"
              << "  Pass criteria: value matches an independent reference\n"
              << "                 implementation of the integrand to 1e-12.\n";

    IntegrandParams params = pirfunc_make_params(/*k=*/10.0);

    // A sweep of x values spanning the region where the integrand is
    // non-negligible.  x is strictly positive because y0(0) diverges.
    const std::vector<double> xValues { 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0 };

    for (double x : xValues) {
        const double expected = pirfunc_reference(x, params);
        const double actual = pir_function(x, static_cast<void*>(&params));

        std::cerr << "    x = " << x << " -> pir_function = " << actual
                  << ", reference = " << expected << '\n';

        // Both implementations use identical elementary operations, so exact
        // agreement is expected; a tiny absolute tolerance guards rounding.
        EXPECT_NEAR(actual, expected, 1e-12)
            << "pir_function disagreed with the reference formula at x = " << x;
        EXPECT_TRUE(std::isfinite(actual))
            << "pir_function returned a non-finite value at x = " << x;
    }
}

// -----------------------------------------------------------------------------
// Test 2: k == infinity selects the absorbing-boundary branch, which must also
//         match the reference and must differ from the finite-k result.
// -----------------------------------------------------------------------------
void pirfunc_test_absorbing_branch()
{
    std::cerr << "\n[TEST] pirfunc_test_absorbing_branch\n"
              << "  Source file:   src/reactions/pir_function.cpp\n"
              << "  Function:      pir_function (else branch, k == infinity)\n"
              << "  Scenario:      k set to +infinity so 'k < 1.0/0.0' is false.\n"
              << "  Pass criteria: matches the absorbing-boundary reference and\n"
              << "                 differs from the finite-k evaluation.\n";

    IntegrandParams absorbing = pirfunc_make_params(std::numeric_limits<double>::infinity());
    IntegrandParams finiteK = pirfunc_make_params(/*k=*/10.0);

    const std::vector<double> xValues { 0.05, 0.5, 1.5, 4.0 };

    for (double x : xValues) {
        const double expected = pirfunc_reference(x, absorbing);
        const double actual = pir_function(x, static_cast<void*>(&absorbing));
        const double finiteVal = pir_function(x, static_cast<void*>(&finiteK));

        std::cerr << "    x = " << x << " -> absorbing = " << actual
                  << " (reference " << expected << "), finite k = " << finiteVal << '\n';

        EXPECT_NEAR(actual, expected, 1e-12)
            << "Absorbing branch disagreed with reference at x = " << x;
        EXPECT_TRUE(std::isfinite(actual))
            << "Absorbing branch produced a non-finite value at x = " << x;

        // The two boundary conditions are physically distinct, so the results
        // should not coincide (only a pathological x could make them equal).
        EXPECT_NE(actual, finiteVal)
            << "Absorbing and finite-k branches should differ at x = " << x;
    }
}

// -----------------------------------------------------------------------------
// Test 3: the integrand is symmetric under swapping r and r0, since P and T
//         are the same functional form evaluated at those two radii.
// -----------------------------------------------------------------------------
void pirfunc_test_symmetry_r_r0()
{
    std::cerr << "\n[TEST] pirfunc_test_symmetry_r_r0\n"
              << "  Source file:   src/reactions/pir_function.cpp\n"
              << "  Function:      pir_function\n"
              << "  Scenario:      swap the r and r0 fields of IntegrandParams.\n"
              << "  Pass criteria: f(r, r0) == f(r0, r) because the product P*T\n"
              << "                 is symmetric in the two radii.\n";

    IntegrandParams forward = pirfunc_make_params(/*k=*/5.0);
    IntegrandParams swapped = forward;
    swapped.r = forward.r0;
    swapped.r0 = forward.r;

    const std::vector<double> xValues { 0.2, 0.9, 2.5 };

    for (double x : xValues) {
        const double a = pir_function(x, static_cast<void*>(&forward));
        const double b = pir_function(x, static_cast<void*>(&swapped));

        std::cerr << "    x = " << x << " -> f(r=" << forward.r << ",r0=" << forward.r0
                  << ") = " << a << ", f(r=" << swapped.r << ",r0=" << swapped.r0
                  << ") = " << b << '\n';

        EXPECT_NEAR(a, b, 1e-12)
            << "Integrand should be symmetric in r <-> r0 at x = " << x;
    }
}

// -----------------------------------------------------------------------------
// Test 4: the only time dependence is the exp(-D t x^2) prefactor, so the ratio
//         of two evaluations at different times is analytically known.
// -----------------------------------------------------------------------------
void pirfunc_test_time_decay_factor()
{
    std::cerr << "\n[TEST] pirfunc_test_time_decay_factor\n"
              << "  Source file:   src/reactions/pir_function.cpp\n"
              << "  Function:      pir_function\n"
              << "  Scenario:      evaluate at two different times t1 < t2.\n"
              << "  Pass criteria: f(t2)/f(t1) == exp(-D*(t2-t1)*x^2).\n";

    IntegrandParams p1 = pirfunc_make_params(/*k=*/8.0);
    p1.t = 0.5;
    IntegrandParams p2 = p1;
    p2.t = 2.0;

    const std::vector<double> xValues { 0.3, 1.0, 2.0 };

    for (double x : xValues) {
        const double f1 = pir_function(x, static_cast<void*>(&p1));
        const double f2 = pir_function(x, static_cast<void*>(&p2));
        const double expectedRatio = std::exp(-p1.D * (p2.t - p1.t) * x * x);

        std::cerr << "    x = " << x << " -> f(t=0.5) = " << f1 << ", f(t=2.0) = " << f2
                  << ", ratio = " << (f2 / f1) << " (expected " << expectedRatio << ")\n";

        // Guard against dividing by a value that happens to be ~0.
        if (std::fabs(f1) > 1e-30) {
            EXPECT_NEAR(f2 / f1, expectedRatio, 1e-10)
                << "Time-decay prefactor is incorrect at x = " << x;
        } else {
            std::cerr << "      (skipping ratio check: f(t=0.5) is essentially zero)\n";
        }

        // The exponential factor is < 1 for t2 > t1, so magnitude must shrink.
        EXPECT_LE(std::fabs(f2), std::fabs(f1) + 1e-15)
            << "Magnitude should not grow with increasing time at x = " << x;
    }
}

// -----------------------------------------------------------------------------
// Test 5: at t == 0 the exponential is unity, so f reduces to x*P*T/(2 pi).
// -----------------------------------------------------------------------------
void pirfunc_test_zero_time_reduces_to_xPT()
{
    std::cerr << "\n[TEST] pirfunc_test_zero_time_reduces_to_xPT\n"
              << "  Source file:   src/reactions/pir_function.cpp\n"
              << "  Function:      pir_function\n"
              << "  Scenario:      t = 0 removes the exponential damping.\n"
              << "  Pass criteria: f(x) == x * P(r) * T(r0) / (2 pi).\n";

    IntegrandParams params = pirfunc_make_params(/*k=*/3.0);
    params.t = 0.0;

    const std::vector<double> xValues { 0.25, 1.0, 3.0 };

    for (double x : xValues) {
        // Recompute P and T locally to build the expected value explicitly.
        const double h = 2.0 * M_PI * params.a * params.D;
        const double alp = h * x * j1(x * params.a) + params.k * j0(x * params.a);
        const double bet = h * x * y1(x * params.a) + params.k * y0(x * params.a);
        const double tet = std::sqrt(alp * alp + bet * bet);
        const double P = (j0(x * params.r) * bet - y0(x * params.r) * alp) / tet;
        const double T = (j0(x * params.r0) * bet - y0(x * params.r0) * alp) / tet;
        const double expected = x * P * T / (2.0 * M_PI);

        const double actual = pir_function(x, static_cast<void*>(&params));

        std::cerr << "    x = " << x << " -> f = " << actual << ", x*P*T/(2pi) = " << expected
                  << '\n';

        EXPECT_NEAR(actual, expected, 1e-12)
            << "At t = 0 the integrand should reduce to x*P*T/(2 pi) at x = " << x;
    }
}

// -----------------------------------------------------------------------------
// Test 6: as k grows very large, the finite-k branch must converge to the
//         absorbing-boundary result (the k factor cancels in P/tet and T/tet).
// -----------------------------------------------------------------------------
void pirfunc_test_large_k_approaches_absorbing()
{
    std::cerr << "\n[TEST] pirfunc_test_large_k_approaches_absorbing\n"
              << "  Source file:   src/reactions/pir_function.cpp\n"
              << "  Function:      pir_function\n"
              << "  Scenario:      compare huge-but-finite k with k == infinity.\n"
              << "  Pass criteria: the finite-k value converges to the absorbing\n"
              << "                 branch value (relative difference < 1e-6).\n";

    IntegrandParams absorbing = pirfunc_make_params(std::numeric_limits<double>::infinity());
    IntegrandParams hugeK = pirfunc_make_params(/*k=*/1.0e12);

    const std::vector<double> xValues { 0.1, 0.7, 2.0 };

    for (double x : xValues) {
        const double fAbs = pir_function(x, static_cast<void*>(&absorbing));
        const double fBig = pir_function(x, static_cast<void*>(&hugeK));

        // Use a scale-aware tolerance so tiny values are not over-constrained.
        const double scale = std::max(1.0e-12, std::fabs(fAbs));

        std::cerr << "    x = " << x << " -> absorbing = " << fAbs << ", k=1e12 = " << fBig
                  << ", |diff|/scale = " << (std::fabs(fBig - fAbs) / scale) << '\n';

        EXPECT_NEAR(fBig, fAbs, 1e-6 * scale)
            << "Large finite k should converge to the absorbing branch at x = " << x;
    }
}

// -----------------------------------------------------------------------------
// Test 7: pir_function must treat its parameter block as read-only.  The GSL
//         integrator reuses the same struct across many evaluations, so any
//         mutation would silently corrupt the lookup tables.
// -----------------------------------------------------------------------------
void pirfunc_test_params_are_not_modified()
{
    std::cerr << "\n[TEST] pirfunc_test_params_are_not_modified\n"
              << "  Source file:   src/reactions/pir_function.cpp\n"
              << "  Function:      pir_function\n"
              << "  Scenario:      call the integrand and re-inspect the struct.\n"
              << "  Pass criteria: every IntegrandParams field is unchanged.\n";

    IntegrandParams params = pirfunc_make_params(/*k=*/7.5);

    // Snapshot the inputs before the call.
    const double a0 = params.a;
    const double D0 = params.D;
    const double k0 = params.k;
    const double r00 = params.r0;
    const double rr0 = params.r;
    const double t0 = params.t;
    const double rho0 = params.rho;

    const double value = pir_function(1.25, static_cast<void*>(&params));
    std::cerr << "    pir_function(1.25, params) = " << value << '\n';

    EXPECT_DOUBLE_EQ(params.a, a0) << "field 'a' (binding radius) was modified";
    EXPECT_DOUBLE_EQ(params.D, D0) << "field 'D' (diffusion constant) was modified";
    EXPECT_DOUBLE_EQ(params.k, k0) << "field 'k' (rate) was modified";
    EXPECT_DOUBLE_EQ(params.r0, r00) << "field 'r0' was modified";
    EXPECT_DOUBLE_EQ(params.r, rr0) << "field 'r' was modified";
    EXPECT_DOUBLE_EQ(params.t, t0) << "field 't' was modified";
    EXPECT_DOUBLE_EQ(params.rho, rho0) << "field 'rho' was modified";
}

// -----------------------------------------------------------------------------
// Test 8: robustness sweep - the integrand should stay finite over a broad
//         range of x for both a reflecting-like (k == 0) and a moderate k, and
//         it should decay toward zero for large x because of exp(-D t x^2).
// -----------------------------------------------------------------------------
void pirfunc_test_finiteness_and_large_x_decay()
{
    std::cerr << "\n[TEST] pirfunc_test_finiteness_and_large_x_decay\n"
              << "  Source file:   src/reactions/pir_function.cpp\n"
              << "  Function:      pir_function\n"
              << "  Scenario:      sweep x for k = 0 and k = 25, then probe large x.\n"
              << "  Pass criteria: all values finite; large-x value is ~0.\n";

    const std::vector<double> kValues { 0.0, 25.0 };
    const std::vector<double> xValues { 1e-4, 1e-2, 0.1, 1.0, 3.0, 8.0, 20.0 };

    for (double k : kValues) {
        IntegrandParams params = pirfunc_make_params(k);
        std::cerr << "  -> k = " << k << '\n';

        for (double x : xValues) {
            const double f = pir_function(x, static_cast<void*>(&params));
            std::cerr << "     x = " << x << " -> f = " << f << '\n';

            EXPECT_TRUE(std::isfinite(f))
                << "Integrand must be finite for k = " << k << ", x = " << x;

            // Cross-check against the reference for every sampled point.
            EXPECT_NEAR(f, pirfunc_reference(x, params), 1e-12)
                << "Reference mismatch for k = " << k << ", x = " << x;
        }

        // With D*t = 0.5 the Gaussian factor at x = 40 is exp(-800), i.e. zero
        // to double precision, so the integrand must vanish there.
        const double fLarge = pir_function(40.0, static_cast<void*>(&params));
        std::cerr << "     x = 40 -> f = " << fLarge << " (expected ~0)\n";
        EXPECT_NEAR(fLarge, 0.0, 1e-30)
            << "Integrand should be numerically zero at large x for k = " << k;
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named helper is run inside its own TEST so that a
// failure in one scenario does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(PirFunction, FiniteKMatchesReference) { pirfunc_test_finite_k_matches_reference(); }
TEST(PirFunction, AbsorbingBranch) { pirfunc_test_absorbing_branch(); }
TEST(PirFunction, SymmetryInRAndR0) { pirfunc_test_symmetry_r_r0(); }
TEST(PirFunction, TimeDecayFactor) { pirfunc_test_time_decay_factor(); }
TEST(PirFunction, ZeroTimeReducesToXPT) { pirfunc_test_zero_time_reduces_to_xPT(); }
TEST(PirFunction, LargeKApproachesAbsorbing) { pirfunc_test_large_k_approaches_absorbing(); }
TEST(PirFunction, ParamsAreNotModified) { pirfunc_test_params_are_not_modified(); }
TEST(PirFunction, FinitenessAndLargeXDecay) { pirfunc_test_finiteness_and_large_x_decay(); }