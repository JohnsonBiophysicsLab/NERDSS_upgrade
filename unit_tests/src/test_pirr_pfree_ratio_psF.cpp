/*! \file test_pirr_pfree_ratio_psF.cpp
 *
 * ### Unit test for src/reactions/pirr_pfree_ratio_psF.cpp
 *
 * The single function under test is
 *
 *     double pirr_pfree_ratio_psF(double rCurr, double r0, double tCurr,
 *                                 double Dtot, double bindrad, double alpha,
 *                                 double ps_prev, double rtol);
 *
 * It computes the ratio between
 *   * `pirr`  : the (normalized) irreversible Green's function for a pair of
 *               molecules that started at separation `r0` and is evaluated at
 *               separation `rCurr` after a time `tCurr`, for a partially
 *               absorbing boundary characterized by `alpha`, and
 *   * `pfree` : the free (no-reaction) propagator with no diffusion inside the
 *               binding radius, renormalized so that it integrates to one from
 *               sigma to infinity,
 * scaled by the previously accumulated survival probability `ps_prev`:
 *
 *     ratio = pirr / (pfree * ps_prev)
 *
 * with the special case that whenever `|pirr - pfree*ps_prev| < rtol` the
 * function short-circuits and returns exactly 1.0.
 *
 * The implementation additionally contains two numerically distinct code paths
 * for the `exp(e1) * erfc(ef1)` term:
 *   1. the direct `exp(e1)*erfc(ef1)` evaluation, and
 *   2. an overflow safe path via the Faddeeva function w(z) which is taken when
 *      `exp(e1)` overflows to +inf.
 * Both branches are exercised below.
 *
 * Strategy of the tests:
 *   * An *independent* reference implementation is written in this file which
 *     evaluates the same closed-form expressions but uses the scaled
 *     complementary error function `erfcx(x) = exp(x^2) erfc(x)` instead of
 *     either `exp()*erfc()` or the Faddeeva function.  Because
 *         exp(2ab + b^2) * erfc(a+b) == exp(-a^2) * erfcx(a+b)
 *     this reference is mathematically identical but numerically stable for all
 *     inputs, which lets us check the production code in *both* of its branches.
 *   * Structural properties (rtol clamping, 1/ps_prev scaling, monotonicity in
 *     alpha, positivity/finiteness) are checked separately so a failure points
 *     at a specific piece of behaviour.
 */

#include "math/Faddeeva.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

// -----------------------------------------------------------------------------
// Independent reference implementation of the quantity under test.
//
// Everything here mirrors the closed-form math of the production routine, but
// the potentially overflowing product exp(e1)*erfc(ef1) is replaced by the
// algebraically identical, always-finite form exp(-a^2)*erfcx(a+b).
// -----------------------------------------------------------------------------
double pirrpf_reference_ratio(double rCurr, double r0, double tCurr, double Dtot,
    double bindrad, double alpha, double ps_prev, double rtol)
{
    const double fDt = 4.0 * Dtot * tCurr; // 4 D t
    const double sq_fDt = std::sqrt(fDt);
    const double sqrt_t = std::sqrt(tCurr);

    const double f1 = 1.0 / std::sqrt(4.0 * M_PI * tCurr);
    const double f2 = 1.0 / (4.0 * M_PI * r0 * std::sqrt(Dtot));

    const double dist = rCurr - r0; // "direct" image distance
    const double sep = rCurr + r0 - 2.0 * bindrad; // reflected image distance

    // Sum of the direct and the mirror Gaussian.
    const double term1 = f1 * (std::exp(-dist * dist / fDt) + std::exp(-sep * sep / fDt));

    // Absorbing (radiation boundary) correction, evaluated in a scaled fashion:
    //   exp(2ab + b^2) erfc(a+b) == exp(-a^2) erfcx(a+b)
    const double a = sep / sq_fDt;
    const double b = sqrt_t * alpha;
    const double term2 = std::exp(-a * a) * Faddeeva::erfcx(a + b);

    const double pirr = (term1 - alpha * term2) * f2 / rCurr;

    // Normalization of pfree: integral of the free propagator from sigma to inf.
    const double cof = f1 * f2;
    const double c1 = 4.0 * M_PI * cof;
    const double ndist = bindrad - r0;
    const double nadist = bindrad + r0;
    const double sq_P = std::sqrt(M_PI);
    const double n1 = -0.5 * fDt * std::exp(-ndist * ndist / fDt)
        - 0.5 * sq_fDt * sq_P * r0 * std::erf(-ndist / sq_fDt);
    const double n2 = 0.5 * fDt * std::exp(-nadist * nadist / fDt)
        + 0.5 * sq_fDt * sq_P * r0 * std::erf(nadist / sq_fDt);
    const double pnorm = 1.0 - c1 * (n1 + n2);

    const double adist = rCurr + r0;
    const double pfree = cof / rCurr
        * (std::exp(-dist * dist / fDt) - std::exp(-adist * adist / fDt)) / pnorm;

    // The production code clamps to unity when the two densities are within rtol.
    if (std::abs(pirr - pfree * ps_prev) < rtol)
        return 1.0;

    return pirr / (pfree * ps_prev);
}

/*! \brief Helper: is the exponent used by the production code going to overflow?
 *
 * The production code takes the Faddeeva branch exactly when exp(e1) == inf,
 * with e1 = 2*a*b + alpha^2 * t.  Reproducing e1 here lets the test report which
 * branch it believes it is exercising.
 */
double pirrpf_exponent_e1(double rCurr, double r0, double tCurr, double Dtot,
    double bindrad, double alpha)
{
    const double fDt = 4.0 * Dtot * tCurr;
    const double a = (rCurr + r0 - 2.0 * bindrad) / std::sqrt(fDt);
    const double b = std::sqrt(tCurr) * alpha;
    return 2.0 * a * b + alpha * alpha * tCurr;
}

/*! \brief Relative comparison with a small absolute floor, reported verbosely. */
void pirrpf_expect_relclose(double actual, double expected, double relTol, const char* what)
{
    const double tol = relTol * std::fabs(expected) + 1e-14;
    std::cerr << "      " << what << ": actual = " << std::setprecision(12) << actual
              << ", reference = " << expected << ", tol = " << tol << '\n';
    EXPECT_NEAR(actual, expected, tol) << what;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the returned ratio matches the independent closed-form reference for
//         a set of physically reasonable parameter combinations (direct
//         exp()*erfc() branch of the implementation).
// -----------------------------------------------------------------------------
void test_pirrpf_matches_analytic_reference()
{
    std::cerr << "\n[TEST] test_pirrpf_matches_analytic_reference\n"
              << "  Source file:   src/reactions/pirr_pfree_ratio_psF.cpp\n"
              << "  Function:      pirr_pfree_ratio_psF()\n"
              << "  Scenario:      several typical (r, r0, t, D, sigma, alpha, ps) sets.\n"
              << "  Pass criteria: value agrees with the erfcx-based reference to 1e-9\n"
              << "                 relative error.\n";

    // Each row: rCurr, r0, tCurr, Dtot, bindrad, alpha, ps_prev
    const double cases[][7] = {
        { 1.5, 1.2, 0.5, 20.0, 1.0, 5.0, 0.90 }, // "generic" 3D case
        { 1.0, 1.0, 0.1, 10.0, 1.0, 1.0, 1.00 }, // starts and ends at contact
        { 3.0, 1.5, 1.0, 5.0, 1.0, 0.5, 0.75 }, // well separated, slow binding
        { 2.0, 2.0, 0.01, 100.0, 1.0, 20.0, 0.99 }, // fast diffusion, short step
        { 1.2, 5.0, 0.25, 15.0, 1.0, 2.0, 0.50 }, // moving inward from far away
    };

    for (const auto& c : cases) {
        const double rCurr = c[0], r0 = c[1], t = c[2], D = c[3];
        const double sigma = c[4], alpha = c[5], ps = c[6];
        const double rtol = 0.0; // disable the "close enough -> 1.0" clamp

        std::cerr << "  -> r=" << rCurr << " r0=" << r0 << " t=" << t << " D=" << D
                  << " sigma=" << sigma << " alpha=" << alpha << " ps_prev=" << ps
                  << " (e1 = " << pirrpf_exponent_e1(rCurr, r0, t, D, sigma, alpha) << ")\n";

        const double actual = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, rtol);
        const double expected = pirrpf_reference_ratio(rCurr, r0, t, D, sigma, alpha, ps, rtol);

        pirrpf_expect_relclose(actual, expected, 1e-9, "pirr/pfree ratio");

        // A physically meaningful ratio must be finite and positive here.
        EXPECT_TRUE(std::isfinite(actual)) << "ratio must be finite for well posed input";
        EXPECT_GT(actual, 0.0) << "ratio must be positive for well posed input";
    }
}

// -----------------------------------------------------------------------------
// Test 2: when rtol is large the two densities are considered indistinguishable
//         and the function must return exactly 1.0.
// -----------------------------------------------------------------------------
void test_pirrpf_rtol_clamps_to_unity()
{
    std::cerr << "\n[TEST] test_pirrpf_rtol_clamps_to_unity\n"
              << "  Function:      pirr_pfree_ratio_psF()\n"
              << "  Scenario:      rtol is huge, so |pirr - pfree*ps| < rtol always.\n"
              << "  Pass criteria: returned value is exactly 1.0 (short-circuit path).\n";

    const double rCurr = 1.5, r0 = 1.2, t = 0.5, D = 20.0, sigma = 1.0;
    const double alpha = 5.0, ps = 0.9;

    // With rtol = 1e6 the difference of two densities (order 1e-3) is inside it.
    const double clamped = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, 1.0e6);
    std::cerr << "  -> with rtol = 1e6 the function returned " << clamped << '\n';
    EXPECT_DOUBLE_EQ(clamped, 1.0) << "large rtol must force the ratio to unity";

    // Sanity check: with rtol = 0 the same input does *not* give exactly 1.
    const double unclamped = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, 0.0);
    std::cerr << "  -> with rtol = 0   the function returned " << unclamped << '\n';
    EXPECT_NE(unclamped, 1.0) << "with rtol = 0 the true (non-unity) ratio should be returned";
}

// -----------------------------------------------------------------------------
// Test 3: outside the clamping region the ratio must scale as 1/ps_prev, since
//         ps_prev only appears in the denominator.
// -----------------------------------------------------------------------------
void test_pirrpf_ps_prev_inverse_scaling()
{
    std::cerr << "\n[TEST] test_pirrpf_ps_prev_inverse_scaling\n"
              << "  Function:      pirr_pfree_ratio_psF()\n"
              << "  Scenario:      identical geometry, ps_prev = 1.0 vs 0.5 vs 0.25.\n"
              << "  Pass criteria: ratio(ps) * ps is constant (ratio ~ 1/ps_prev).\n";

    const double rCurr = 1.8, r0 = 1.3, t = 0.4, D = 25.0, sigma = 1.0;
    const double alpha = 3.0, rtol = 0.0;

    const double ratio1 = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, 1.00, rtol);
    const double ratio2 = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, 0.50, rtol);
    const double ratio4 = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, 0.25, rtol);

    std::cerr << "  -> ratio(ps=1.00) = " << ratio1 << '\n'
              << "  -> ratio(ps=0.50) = " << ratio2 << " (expect ~2x)\n"
              << "  -> ratio(ps=0.25) = " << ratio4 << " (expect ~4x)\n";

    pirrpf_expect_relclose(ratio2, 2.0 * ratio1, 1e-12, "ratio at ps=0.5 vs 2*ratio at ps=1");
    pirrpf_expect_relclose(ratio4, 4.0 * ratio1, 1e-12, "ratio at ps=0.25 vs 4*ratio at ps=1");
}

// -----------------------------------------------------------------------------
// Test 4: monotonicity in alpha.  alpha is proportional to the intrinsic
//         association rate; a more absorbing boundary removes probability from
//         pirr, so the ratio must decrease as alpha increases.
// -----------------------------------------------------------------------------
void test_pirrpf_monotonic_in_alpha()
{
    std::cerr << "\n[TEST] test_pirrpf_monotonic_in_alpha\n"
              << "  Function:      pirr_pfree_ratio_psF()\n"
              << "  Scenario:      alpha = 0, 0.1, 1, 10, 100 with everything else fixed.\n"
              << "  Pass criteria: the ratio decreases strictly with increasing alpha,\n"
              << "                 stays finite and positive, and alpha = 0 (purely\n"
              << "                 reflecting) gives a ratio above unity.\n";

    const double rCurr = 1.6, r0 = 1.4, t = 0.3, D = 18.0, sigma = 1.0;
    const double ps = 1.0, rtol = 0.0;

    const double alphas[] = { 0.0, 0.1, 1.0, 10.0, 100.0 };
    double previous = std::numeric_limits<double>::infinity();

    for (double alpha : alphas) {
        const double ratio = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, rtol);
        std::cerr << "  -> alpha = " << std::setw(6) << alpha << "  ratio = " << ratio << '\n';

        EXPECT_TRUE(std::isfinite(ratio)) << "ratio must be finite for alpha = " << alpha;
        EXPECT_GT(ratio, 0.0) << "ratio must remain positive for alpha = " << alpha;
        EXPECT_LT(ratio, previous) << "ratio must decrease as alpha (absorption) grows";

        // Cross-check against the independent reference at every alpha, which
        // also spans both numerical branches as alpha grows.
        pirrpf_expect_relclose(ratio,
            pirrpf_reference_ratio(rCurr, r0, t, D, sigma, alpha, ps, rtol), 1e-7,
            "ratio vs reference");

        previous = ratio;
    }

    // alpha == 0 -> reflecting boundary, pirr keeps all of its probability while
    // pfree is normalized with an absorbing-like image term, so ratio > 1.
    const double reflecting = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, 0.0, 1.0, rtol);
    std::cerr << "  -> reflecting limit (alpha = 0) ratio = " << reflecting
              << " (must exceed 1)\n";
    EXPECT_GT(reflecting, 1.0) << "no absorption must leave more density than the free case";
}

// -----------------------------------------------------------------------------
// Test 5: the Faddeeva (overflow) branch.  Parameters are chosen so that
//         e1 = 2ab + alpha^2 t exceeds ~709 and exp(e1) becomes +inf, forcing
//         the implementation into the complex-w() code path.
// -----------------------------------------------------------------------------
void test_pirrpf_faddeeva_overflow_branch()
{
    std::cerr << "\n[TEST] test_pirrpf_faddeeva_overflow_branch\n"
              << "  Function:      pirr_pfree_ratio_psF() -- std::isinf(exp(e1)) branch\n"
              << "  Scenario:      very large alpha so that exp(e1) overflows to +inf.\n"
              << "  Pass criteria: the result stays finite/positive and matches the\n"
              << "                 overflow-free erfcx reference.\n";

    const double rCurr = 1.5, r0 = 1.2, t = 0.5, D = 20.0, sigma = 1.0;
    const double ps = 1.0, rtol = 0.0;
    const double alpha = 40.0; // b = alpha*sqrt(t) ~ 28.3 -> b^2 ~ 800

    const double e1 = pirrpf_exponent_e1(rCurr, r0, t, D, sigma, alpha);
    std::cerr << "  -> e1 = " << e1 << ", exp(e1) = " << std::exp(e1)
              << (std::isinf(std::exp(e1)) ? "  (overflow: Faddeeva branch taken)\n"
                                           : "  (NO overflow: direct branch taken!)\n");
    // Make sure the test is really probing the intended branch.
    EXPECT_TRUE(std::isinf(std::exp(e1)))
        << "test parameters no longer trigger the exp() overflow branch";

    const double actual = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, rtol);
    const double expected = pirrpf_reference_ratio(rCurr, r0, t, D, sigma, alpha, ps, rtol);

    EXPECT_TRUE(std::isfinite(actual)) << "the Faddeeva branch must not produce inf/NaN";
    EXPECT_GT(actual, 0.0) << "the Faddeeva branch must not produce a negative density ratio";
    // The direct and mirror Gaussians nearly cancel here, so allow a looser
    // relative tolerance than in the well-conditioned cases.
    pirrpf_expect_relclose(actual, expected, 1e-5, "overflow-branch ratio vs reference");
}

// -----------------------------------------------------------------------------
// Test 6: sweep the current separation over a grid and make sure the function is
//         numerically well behaved (finite, positive) everywhere outside sigma.
// -----------------------------------------------------------------------------
void test_pirrpf_finite_over_separation_grid()
{
    std::cerr << "\n[TEST] test_pirrpf_finite_over_separation_grid\n"
              << "  Function:      pirr_pfree_ratio_psF()\n"
              << "  Scenario:      r0 fixed at contact, rCurr swept from sigma to sigma+3.\n"
              << "  Pass criteria: every value is finite, positive and equal to the\n"
              << "                 reference implementation.\n";

    const double sigma = 1.0;
    const double r0 = sigma; // start exactly at contact
    const double t = 0.2, D = 12.0, alpha = 4.0, ps = 0.95, rtol = 0.0;

    for (int i = 0; i <= 6; ++i) {
        const double rCurr = sigma + 0.5 * static_cast<double>(i);
        const double ratio = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, rtol);
        const double expected = pirrpf_reference_ratio(rCurr, r0, t, D, sigma, alpha, ps, rtol);

        std::cerr << "  -> rCurr = " << std::setw(5) << rCurr << "  ratio = " << ratio << '\n';

        EXPECT_FALSE(std::isnan(ratio)) << "ratio is NaN at rCurr = " << rCurr;
        EXPECT_FALSE(std::isinf(ratio)) << "ratio is infinite at rCurr = " << rCurr;
        EXPECT_GT(ratio, 0.0) << "ratio must be positive at rCurr = " << rCurr;
        pirrpf_expect_relclose(ratio, expected, 1e-8, "grid point ratio vs reference");
    }
}

// -----------------------------------------------------------------------------
// Test 7: determinism / purity -- the routine has no state, so repeated calls
//         with identical arguments must return bit-identical values.
// -----------------------------------------------------------------------------
void test_pirrpf_is_deterministic()
{
    std::cerr << "\n[TEST] test_pirrpf_is_deterministic\n"
              << "  Function:      pirr_pfree_ratio_psF()\n"
              << "  Scenario:      the same arguments are passed three times.\n"
              << "  Pass criteria: all three return values are bit-identical.\n";

    const double rCurr = 2.2, r0 = 1.1, t = 0.35, D = 30.0, sigma = 1.0;
    const double alpha = 7.5, ps = 0.8, rtol = 1e-14;

    const double first = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, rtol);
    const double second = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, rtol);
    const double third = pirr_pfree_ratio_psF(rCurr, r0, t, D, sigma, alpha, ps, rtol);

    std::cerr << "  -> values: " << std::setprecision(17) << first << ", " << second << ", "
              << third << '\n';

    EXPECT_DOUBLE_EQ(first, second) << "the function must be free of hidden state";
    EXPECT_DOUBLE_EQ(second, third) << "the function must be free of hidden state";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- one per named test_* routine so that a failure in one
// scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(PirrPfreeRatioPsF, MatchesAnalyticReference) { test_pirrpf_matches_analytic_reference(); }
TEST(PirrPfreeRatioPsF, RtolClampsToUnity) { test_pirrpf_rtol_clamps_to_unity(); }
TEST(PirrPfreeRatioPsF, PsPrevInverseScaling) { test_pirrpf_ps_prev_inverse_scaling(); }
TEST(PirrPfreeRatioPsF, MonotonicInAlpha) { test_pirrpf_monotonic_in_alpha(); }
TEST(PirrPfreeRatioPsF, FaddeevaOverflowBranch) { test_pirrpf_faddeeva_overflow_branch(); }
TEST(PirrPfreeRatioPsF, FiniteOverSeparationGrid) { test_pirrpf_finite_over_separation_grid(); }
TEST(PirrPfreeRatioPsF, IsDeterministic) { test_pirrpf_is_deterministic(); }