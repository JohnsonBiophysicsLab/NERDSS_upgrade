/*! \file test_passocF_1D.cpp
 *
 * ### Unit test for ../src/reactions/passocF_1D.cpp
 *
 * Function under test:
 *
 *     double passocF_1D(double r0, double tCurr, double Dtot,
 *                       double bindRadius, double ka)
 *
 * This is the 1D association ("survival"/reaction) probability from the
 * Smoluchowski reaction-diffusion model.  Its implementation has two
 * distinct branches:
 *
 *   1. Dtot == 0  -- a purely static pair.  The routine short-circuits and
 *      returns 1.0 when the separation is not more than 1e-6 beyond the
 *      binding radius, and 0.0 otherwise.
 *
 *   2. Dtot > 0   -- the analytic expression
 *
 *          sqrtfDt    = sqrt(4 * Dtot * tCurr)
 *          sq_scale_t = ka * sqrt(tCurr / Dtot)
 *          sep        = (r0 - bindRadius) / sqrtfDt
 *          p          = erfc(sep) - exp(-sep^2) * erfcx(sep + sq_scale_t)
 *
 *      where erfcx() is the scaled complementary error function supplied by
 *      the vendored Faddeeva package.
 *
 * The tests below verify both branches, the exact algebraic identities that
 * follow from the formula (ka == 0 must give exactly zero probability), an
 * independent re-derivation of one value using erfcx(y) = exp(y^2)*erfc(y),
 * the physically expected monotonic trends, the fact that the result only
 * depends on the two dimensionless groups (sep, sq_scale_t), and the
 * ka -> infinity limit.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "math/Faddeeva.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

// -----------------------------------------------------------------------------
// Test 1: the Dtot == 0 short-circuit branch.
//
// With no diffusion the molecules cannot move, so the routine only asks
// whether they are already (essentially) in contact.  The threshold used in
// the source is a strict "r0 - bindRadius > 1e-6" test, so a separation of
// exactly 1e-6 still counts as "in contact" and returns 1.0.
// -----------------------------------------------------------------------------
void passocf1d_test_zero_diffusion()
{
    std::cerr << "\n[TEST] passocf1d_test_zero_diffusion\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (Dtot == 0 branch)\n"
              << "  Criteria:    returns 1.0 when (r0 - bindRadius) <= 1e-6,\n"
              << "               0.0 when the pair is farther apart than that.\n";

    const double bindRadius = 1.0;
    const double tCurr = 0.5;
    const double ka = 10.0;
    const double Dtot = 0.0;

    // Exactly at contact -> difference is 0, which is not > 1e-6 -> 1.0
    double pContact = passocF_1D(bindRadius, tCurr, Dtot, bindRadius, ka);
    std::cerr << "  r0 == bindRadius            -> p = " << pContact << '\n';
    EXPECT_DOUBLE_EQ(pContact, 1.0)
        << "A static pair sitting exactly at sigma must associate with probability 1";

    // Inside the binding radius -> negative difference -> 1.0
    double pInside = passocF_1D(bindRadius - 0.25, tCurr, Dtot, bindRadius, ka);
    std::cerr << "  r0 <  bindRadius            -> p = " << pInside << '\n';
    EXPECT_DOUBLE_EQ(pInside, 1.0)
        << "A static pair closer than sigma must associate with probability 1";

    // Just inside the 1e-6 tolerance -> still 1.0
    double pTol = passocF_1D(bindRadius + 1.0e-7, tCurr, Dtot, bindRadius, ka);
    std::cerr << "  r0 == bindRadius + 1e-7     -> p = " << pTol << '\n';
    EXPECT_DOUBLE_EQ(pTol, 1.0)
        << "Separation of 1e-7 is within the hard-coded 1e-6 contact tolerance";

    // Clearly beyond the tolerance -> 0.0
    double pFar = passocF_1D(bindRadius + 1.0e-3, tCurr, Dtot, bindRadius, ka);
    std::cerr << "  r0 == bindRadius + 1e-3     -> p = " << pFar << '\n';
    EXPECT_DOUBLE_EQ(pFar, 0.0)
        << "A static pair separated by more than 1e-6 can never associate";
}

// -----------------------------------------------------------------------------
// Test 2: a zero intrinsic rate constant must give exactly zero probability.
//
// With ka == 0 we have sq_scale_t == 0 so the expression collapses to
//     erfc(sep) - exp(-sep^2) * erfcx(sep)
//   = erfc(sep) - exp(-sep^2) * exp(sep^2) * erfc(sep)
//   = 0
// analytically.  Numerically it should be zero to round-off.
// -----------------------------------------------------------------------------
void passocf1d_test_zero_rate_gives_zero()
{
    std::cerr << "\n[TEST] passocf1d_test_zero_rate_gives_zero\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (diffusive branch, ka = 0)\n"
              << "  Criteria:    p is analytically identically zero, so the\n"
              << "               computed value must be zero to round-off.\n";

    const double bindRadius = 1.0;
    const double Dtot = 2.0;
    const double tCurr = 0.75;
    const double ka = 0.0;

    // Sample a handful of starting separations at and beyond contact.
    const std::vector<double> separations{ 0.0, 0.1, 0.5, 2.0, 5.0 };
    for (double ds : separations) {
        double p = passocF_1D(bindRadius + ds, tCurr, Dtot, bindRadius, ka);
        std::cerr << "  r0 - sigma = " << ds << " -> p = " << p << '\n';
        EXPECT_NEAR(p, 0.0, 1.0e-12)
            << "With ka = 0 the two erfc terms cancel exactly; p must vanish";
    }
}

// -----------------------------------------------------------------------------
// Test 3: at contact (r0 == bindRadius) the formula reduces to a closed form.
//
// sep == 0, so p = erfc(0) - exp(0) * erfcx(sq_scale_t) = 1 - erfcx(ka*sqrt(t/D)).
// We recompute that closed form directly from Faddeeva::erfcx and compare.
// -----------------------------------------------------------------------------
void passocf1d_test_at_contact_matches_closed_form()
{
    std::cerr << "\n[TEST] passocf1d_test_at_contact_matches_closed_form\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (diffusive branch, r0 == sigma)\n"
              << "  Criteria:    p == 1 - erfcx(ka*sqrt(t/D)) when sep == 0.\n";

    const double bindRadius = 1.5;
    const double Dtot = 3.0;
    const double tCurr = 0.25;

    const std::vector<double> rates{ 0.5, 2.0, 25.0 };
    for (double ka : rates) {
        const double sqScaleT = ka * std::sqrt(tCurr / Dtot);
        const double expected = 1.0 - Faddeeva::erfcx(sqScaleT);

        double p = passocF_1D(bindRadius, tCurr, Dtot, bindRadius, ka);
        std::cerr << "  ka = " << ka << ", sq_scale_t = " << sqScaleT
                  << " -> p = " << p << ", expected = " << expected << '\n';

        EXPECT_NEAR(p, expected, 1.0e-14)
            << "At contact the formula must reduce to 1 - erfcx(ka*sqrt(t/D))";
        // A probability at contact must be physically sensible.
        EXPECT_GE(p, 0.0) << "Probability at contact must be non-negative";
        EXPECT_LE(p, 1.0) << "Probability at contact must not exceed unity";
    }
}

// -----------------------------------------------------------------------------
// Test 4: independent re-derivation using the identity
//         erfcx(y) = exp(y^2) * erfc(y).
//
// For moderate arguments (y ~ 1) that identity is numerically well behaved, so
// we can rebuild the expected answer purely from std::erfc/std::exp and check
// that the routine (which uses Faddeeva::erfcx) agrees.
//
//   expected = erfc(sep) - exp(y^2 - sep^2) * erfc(y),  y = sep + sq_scale_t
// -----------------------------------------------------------------------------
void passocf1d_test_matches_erfcx_identity()
{
    std::cerr << "\n[TEST] passocf1d_test_matches_erfcx_identity\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (diffusive branch)\n"
              << "  Criteria:    matches an independent evaluation built from\n"
              << "               std::erfc using erfcx(y) = exp(y^2)*erfc(y).\n";

    // Chosen so that sqrt(4*D*t) = 2, sep = 0.25 and sq_scale_t = 0.75,
    // giving y = 1.0 where the identity is numerically harmless.
    const double Dtot = 1.0;
    const double tCurr = 1.0;
    const double bindRadius = 1.0;
    const double r0 = bindRadius + 0.5; // sep = 0.5 / 2 = 0.25
    const double ka = 0.75;             // sq_scale_t = 0.75 * sqrt(1/1)

    const double sqrtfDt = std::sqrt(4.0 * Dtot * tCurr);
    const double sep = (r0 - bindRadius) / sqrtfDt;
    const double sqScaleT = ka * std::sqrt(tCurr / Dtot);
    const double y = sep + sqScaleT;

    const double expected
        = std::erfc(sep) - std::exp(y * y - sep * sep) * std::erfc(y);

    double p = passocF_1D(r0, tCurr, Dtot, bindRadius, ka);
    std::cerr << "  sep = " << sep << ", sq_scale_t = " << sqScaleT
              << ", y = " << y << '\n';
    std::cerr << "  p = " << p << ", independent expected = " << expected << '\n';

    EXPECT_NEAR(p, expected, 1.0e-10)
        << "Faddeeva-based evaluation must agree with the std::erfc identity";
}

// -----------------------------------------------------------------------------
// Test 5: monotonic behaviour in the starting separation r0.
//
// Physically, starting farther apart can only lower the chance of reacting in
// the given time window, so p must decrease as r0 grows and must tend to zero
// for very large separations.
// -----------------------------------------------------------------------------
void passocf1d_test_monotonic_in_r0()
{
    std::cerr << "\n[TEST] passocf1d_test_monotonic_in_r0\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (diffusive branch)\n"
              << "  Criteria:    p strictly decreases with increasing r0 and\n"
              << "               decays to ~0 at very large separations.\n";

    const double bindRadius = 1.0;
    const double Dtot = 1.0;
    const double tCurr = 1.0;
    const double ka = 5.0;

    double prev = passocF_1D(bindRadius, tCurr, Dtot, bindRadius, ka);
    std::cerr << "  r0 - sigma = 0 -> p = " << prev << '\n';

    for (int i = 1; i <= 6; ++i) {
        const double ds = 0.5 * i;
        double curr = passocF_1D(bindRadius + ds, tCurr, Dtot, bindRadius, ka);
        std::cerr << "  r0 - sigma = " << ds << " -> p = " << curr << '\n';

        EXPECT_LT(curr, prev)
            << "Association probability must fall off with starting separation";
        EXPECT_GE(curr, 0.0) << "Probability must remain non-negative";
        prev = curr;
    }

    // Far away, essentially nothing can reach the binding radius in one step.
    double pFar = passocF_1D(bindRadius + 100.0, tCurr, Dtot, bindRadius, ka);
    std::cerr << "  r0 - sigma = 100 -> p = " << pFar << '\n';
    EXPECT_NEAR(pFar, 0.0, 1.0e-9)
        << "At 100 * sqrt(4 D t)/2 separations the probability must vanish";
}

// -----------------------------------------------------------------------------
// Test 6: monotonic behaviour in the elapsed time and in the intrinsic rate.
//
// Longer times and larger intrinsic rates can only increase the probability of
// having reacted.
// -----------------------------------------------------------------------------
void passocf1d_test_monotonic_in_time_and_rate()
{
    std::cerr << "\n[TEST] passocf1d_test_monotonic_in_time_and_rate\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (diffusive branch)\n"
              << "  Criteria:    p increases with tCurr (fixed ka) and with ka\n"
              << "               (fixed tCurr).\n";

    const double bindRadius = 1.0;
    const double r0 = 1.2;
    const double Dtot = 2.0;

    // ---- increasing time at fixed rate ----
    {
        const double ka = 4.0;
        double prev = passocF_1D(r0, 0.01, Dtot, bindRadius, ka);
        std::cerr << "  t = 0.01 -> p = " << prev << '\n';
        const std::vector<double> times{ 0.05, 0.2, 1.0, 5.0, 50.0 };
        for (double t : times) {
            double curr = passocF_1D(r0, t, Dtot, bindRadius, ka);
            std::cerr << "  t = " << t << " -> p = " << curr << '\n';
            EXPECT_GT(curr, prev)
                << "Probability of having reacted must grow with elapsed time";
            prev = curr;
        }
        // Long-time limit in 1D: essentially certain to react.
        double pLong = passocF_1D(r0, 1.0e8, Dtot, bindRadius, ka);
        std::cerr << "  t = 1e8 -> p = " << pLong << '\n';
        EXPECT_GT(pLong, 0.99)
            << "In 1D the pair is recurrent; long times should approach unity";
        EXPECT_LE(pLong, 1.0) << "Probability must never exceed unity";
    }

    // ---- increasing rate at fixed time ----
    {
        const double tCurr = 0.5;
        double prev = passocF_1D(r0, tCurr, Dtot, bindRadius, 0.0);
        std::cerr << "  ka = 0 -> p = " << prev << '\n';
        const std::vector<double> rates{ 0.1, 1.0, 10.0, 100.0, 1000.0 };
        for (double ka : rates) {
            double curr = passocF_1D(r0, tCurr, Dtot, bindRadius, ka);
            std::cerr << "  ka = " << ka << " -> p = " << curr << '\n';
            EXPECT_GT(curr, prev)
                << "A larger intrinsic rate must give a larger reaction probability";
            prev = curr;
        }
    }
}

// -----------------------------------------------------------------------------
// Test 7: the returned value is a valid probability over a parameter grid.
//
// For every physical input (r0 >= bindRadius, positive D, t, ka) the result
// must be finite and lie in [0, 1].
// -----------------------------------------------------------------------------
void passocf1d_test_bounds_over_grid()
{
    std::cerr << "\n[TEST] passocf1d_test_bounds_over_grid\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (diffusive branch)\n"
              << "  Criteria:    finite result in [0,1] for a grid of physical\n"
              << "               inputs (r0 >= sigma).\n";

    const double bindRadius = 1.0;
    const std::vector<double> diffusions{ 0.001, 0.1, 1.0, 100.0 };
    const std::vector<double> times{ 1.0e-6, 1.0e-3, 1.0, 1.0e3 };
    const std::vector<double> rates{ 0.0, 0.5, 50.0, 5000.0 };
    const std::vector<double> gaps{ 0.0, 0.01, 1.0, 25.0 };

    int checked = 0;
    for (double D : diffusions) {
        for (double t : times) {
            for (double ka : rates) {
                for (double gap : gaps) {
                    double p = passocF_1D(bindRadius + gap, t, D, bindRadius, ka);
                    ++checked;

                    EXPECT_TRUE(std::isfinite(p))
                        << "Non-finite probability for D=" << D << " t=" << t
                        << " ka=" << ka << " gap=" << gap;
                    EXPECT_GE(p, -1.0e-12)
                        << "Negative probability for D=" << D << " t=" << t
                        << " ka=" << ka << " gap=" << gap;
                    EXPECT_LE(p, 1.0 + 1.0e-12)
                        << "Probability above unity for D=" << D << " t=" << t
                        << " ka=" << ka << " gap=" << gap;
                }
            }
        }
    }
    std::cerr << "  Checked " << checked << " parameter combinations; all must be\n"
              << "  finite and within [0,1].\n";
}

// -----------------------------------------------------------------------------
// Test 8: the result depends only on the two dimensionless groups
//         sep = (r0 - sigma)/sqrt(4*D*t)  and  s = ka*sqrt(t/D).
//
// Two very different physical parameter sets that share the same (sep, s) must
// therefore produce the same number.
// -----------------------------------------------------------------------------
void passocf1d_test_scaling_invariance()
{
    std::cerr << "\n[TEST] passocf1d_test_scaling_invariance\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (diffusive branch)\n"
              << "  Criteria:    identical (sep, ka*sqrt(t/D)) pairs give\n"
              << "               identical probabilities.\n";

    const double bindRadius = 1.0;

    // Reference set: D = 1, t = 1 -> sqrt(4Dt) = 2, gap = 0.6 -> sep = 0.30,
    //                ka = 3 -> s = 3.
    const double pRef = passocF_1D(bindRadius + 0.6, 1.0, 1.0, bindRadius, 3.0);

    // Rescaled set: D = 4, t = 1 -> sqrt(4Dt) = 4, so gap must be 1.2 to keep
    //               sep = 0.30; s = ka*sqrt(1/4) = ka/2 -> ka = 6.
    const double pScaled = passocF_1D(bindRadius + 1.2, 1.0, 4.0, bindRadius, 6.0);

    // Another rescaling in time: D = 1, t = 4 -> sqrt(4Dt) = 4, gap = 1.2;
    //               s = ka*sqrt(4/1) = 2*ka -> ka = 1.5.
    const double pTimeScaled
        = passocF_1D(bindRadius + 1.2, 4.0, 1.0, bindRadius, 1.5);

    std::cerr << "  reference   p = " << pRef << '\n';
    std::cerr << "  D-rescaled  p = " << pScaled << '\n';
    std::cerr << "  t-rescaled  p = " << pTimeScaled << '\n';

    EXPECT_NEAR(pScaled, pRef, 1.0e-12)
        << "Rescaling the diffusion constant while holding sep and s fixed "
           "must not change the probability";
    EXPECT_NEAR(pTimeScaled, pRef, 1.0e-12)
        << "Rescaling the time while holding sep and s fixed must not change "
           "the probability";
}

// -----------------------------------------------------------------------------
// Test 9: the ka -> infinity (diffusion limited) limit.
//
// As sq_scale_t grows, erfcx(sep + s) ~ 1/(sqrt(pi)*(sep+s)) -> 0, so the
// expression tends to erfc(sep): the classic absorbing-boundary result.
// This also exercises Faddeeva::erfcx in its large-argument regime, where a
// naive exp(y^2)*erfc(y) would overflow / underflow.
// -----------------------------------------------------------------------------
void passocf1d_test_large_ka_limit()
{
    std::cerr << "\n[TEST] passocf1d_test_large_ka_limit\n"
              << "  Source file: src/reactions/passocF_1D.cpp\n"
              << "  Function:    passocF_1D (diffusive branch, ka -> infinity)\n"
              << "  Criteria:    p -> erfc(sep), the absorbing-boundary result,\n"
              << "               and stays finite for huge ka.\n";

    const double bindRadius = 1.0;
    const double Dtot = 1.0;
    const double tCurr = 1.0;
    const double r0 = bindRadius + 0.5;      // sep = 0.5 / sqrt(4) = 0.25
    const double sep = (r0 - bindRadius) / std::sqrt(4.0 * Dtot * tCurr);
    const double absorbingLimit = std::erfc(sep);

    double pBig = passocF_1D(r0, tCurr, Dtot, bindRadius, 1.0e8);
    double pHuge = passocF_1D(r0, tCurr, Dtot, bindRadius, 1.0e30);

    std::cerr << "  sep = " << sep << ", erfc(sep) = " << absorbingLimit << '\n';
    std::cerr << "  ka = 1e8  -> p = " << pBig << '\n';
    std::cerr << "  ka = 1e30 -> p = " << pHuge << '\n';

    EXPECT_TRUE(std::isfinite(pBig)) << "Large ka must not produce NaN/Inf";
    EXPECT_TRUE(std::isfinite(pHuge)) << "Very large ka must not produce NaN/Inf";
    EXPECT_NEAR(pBig, absorbingLimit, 1.0e-6)
        << "For ka = 1e8 the answer should already be at the erfc(sep) limit";
    EXPECT_NEAR(pHuge, absorbingLimit, 1.0e-12)
        << "For effectively infinite ka the answer must equal erfc(sep)";
    // The finite-rate answer can never exceed the diffusion-limited answer.
    double pModerate = passocF_1D(r0, tCurr, Dtot, bindRadius, 1.0);
    std::cerr << "  ka = 1    -> p = " << pModerate << '\n';
    EXPECT_LT(pModerate, absorbingLimit)
        << "A finite intrinsic rate must give less reaction than the "
           "diffusion-limited (absorbing) case";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is run inside its own TEST so that a
// failure in one scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(PassocF1D, ZeroDiffusionBranch) { passocf1d_test_zero_diffusion(); }
TEST(PassocF1D, ZeroRateGivesZero) { passocf1d_test_zero_rate_gives_zero(); }
TEST(PassocF1D, AtContactMatchesClosedForm) { passocf1d_test_at_contact_matches_closed_form(); }
TEST(PassocF1D, MatchesErfcxIdentity) { passocf1d_test_matches_erfcx_identity(); }
TEST(PassocF1D, MonotonicInR0) { passocf1d_test_monotonic_in_r0(); }
TEST(PassocF1D, MonotonicInTimeAndRate) { passocf1d_test_monotonic_in_time_and_rate(); }
TEST(PassocF1D, BoundsOverParameterGrid) { passocf1d_test_bounds_over_grid(); }
TEST(PassocF1D, ScalingInvariance) { passocf1d_test_scaling_invariance(); }
TEST(PassocF1D, LargeKaLimit) { passocf1d_test_large_ka_limit(); }