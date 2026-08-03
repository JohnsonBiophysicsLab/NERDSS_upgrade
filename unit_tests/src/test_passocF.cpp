/*! \file test_passocF.cpp
 *
 * ### Unit test for ../src/reactions/passocF.cpp
 *
 * Function under test:
 *
 *     double passocF(double r0, double tCurr, double Dtot,
 *                    double bindRadius, double alpha, double cof)
 *
 * `passocF` evaluates the analytic 3D association probability of a pair of
 * molecules that start at separation `r0` after a time `tCurr`, for a
 * partially absorbing (radiation) boundary condition at `bindRadius`:
 *
 *     sep  = (r0 - bindRadius) / sqrt(4 * Dtot * tCurr)
 *     ef1  = sep + alpha * sqrt(tCurr)
 *     e1   = 2 * sep * sqrt(tCurr) * alpha + alpha^2 * tCurr   ( = ef1^2 - sep^2 )
 *     f1   = cof * bindRadius / r0
 *
 *     p    = [ erfc(sep) - exp(e1) * erfc(ef1) ] * f1
 *          = [ erfc(sep) - exp(-sep^2) * erfcx(ef1) ] * f1
 *
 * The source contains two numerical branches which must give the *same*
 * mathematical answer:
 *   1. the direct branch, used when exp(e1) does not overflow, and
 *   2. the scaled branch, used when exp(e1) overflows to +inf, in which case
 *      the scaled complementary error function erfcx(ef1) is obtained from the
 *      Faddeeva function, Re{w(i*ef1)} == erfcx(ef1).
 *
 * The tests below therefore:
 *   - compare passocF against an independent reference built from
 *     std::erfc() and Faddeeva::erfcx() (numerically stable for both ranges),
 *   - verify the analytically required limits (alpha == 0 -> p == 0),
 *   - verify physical monotonicity (p grows with time, shrinks with r0),
 *   - verify the bound 0 <= p <= cof*bindRadius/r0,
 *   - verify exact linearity in the prefactor `cof`, and
 *   - explicitly exercise the Faddeeva (overflow) branch and check that it is
 *     continuous with the direct branch across the overflow threshold.
 *
 * Verbose progress information is written to stderr so the reader can follow
 * exactly which behaviour is being probed.
 */

#include "math/Faddeeva.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Independent reference implementation of the analytic expression.
//
// It uses the *scaled* complementary error function everywhere
// (erfcx(x) = exp(x^2) * erfc(x)), so it never overflows and is valid for both
// the small-alpha and the huge-alpha regimes.  This makes it a genuine
// cross-check of the branch logic inside passocF().
// -----------------------------------------------------------------------------
double passocf_reference(double r0, double tCurr, double Dtot, double bindRadius,
    double alpha, double cof)
{
    const double sqrtfDt = std::sqrt(4.0 * Dtot * tCurr);
    const double sep = (r0 - bindRadius) / sqrtfDt;
    const double ef1 = sep + alpha * std::sqrt(tCurr);
    const double f1 = cof * bindRadius / r0;

    // erfc(sep) - exp(-sep^2) * erfcx(ef1)
    return (std::erfc(sep) - std::exp(-sep * sep) * Faddeeva::erfcx(ef1)) * f1;
}

// -----------------------------------------------------------------------------
// Helper that reports which internal branch of passocF() a parameter set hits.
// The source switches to the Faddeeva code path when exp(e1) is infinite.
// -----------------------------------------------------------------------------
bool passocf_hits_faddeeva_branch(double r0, double tCurr, double Dtot,
    double bindRadius, double alpha)
{
    const double sqrtfDt = std::sqrt(4.0 * Dtot * tCurr);
    const double sep = (r0 - bindRadius) / sqrtfDt;
    const double e1 = 2.0 * sep * std::sqrt(tCurr) * alpha + alpha * alpha * tCurr;
    return std::isinf(std::exp(e1));
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: alpha == 0 means no intrinsic reactivity, so the probability must be
//         identically zero (term1 and term2 become bit-for-bit identical).
// -----------------------------------------------------------------------------
void test_passocf_zero_alpha_gives_zero_probability()
{
    std::cerr << "\n[TEST] test_passocf_zero_alpha_gives_zero_probability\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF()\n"
              << "  Scenario:      alpha = 0 (no intrinsic association rate).\n"
              << "  Pass criteria: p == 0 exactly, because erfc(sep) and\n"
              << "                 exp(0)*erfc(sep) cancel identically.\n";

    // A handful of otherwise unrelated geometries / times.
    const double bindRadius = 1.0;
    const double cof = 1.0;
    const std::vector<double> r0List { 1.0, 1.5, 5.0 };
    const std::vector<double> tList { 1e-6, 1e-3, 1.0 };
    const double Dtot = 10.0;

    for (double r0 : r0List) {
        for (double t : tList) {
            const double p = passocF(r0, t, Dtot, bindRadius, /*alpha=*/0.0, cof);
            std::cerr << "    r0 = " << r0 << ", tCurr = " << t
                      << " -> p = " << p << '\n';
            EXPECT_DOUBLE_EQ(p, 0.0)
                << "With alpha = 0 the association probability must vanish "
                   "(r0 = " << r0 << ", tCurr = " << t << ")";
        }
    }
}

// -----------------------------------------------------------------------------
// Test 2: Molecules starting exactly at contact (r0 == bindRadius, sep == 0).
//         The analytic result collapses to (1 - erfcx(alpha*sqrt(t))) * cof.
// -----------------------------------------------------------------------------
void test_passocf_at_contact_matches_closed_form()
{
    std::cerr << "\n[TEST] test_passocf_at_contact_matches_closed_form\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF()\n"
              << "  Scenario:      r0 == bindRadius so sep == 0.\n"
              << "  Pass criteria: p == (1 - erfcx(alpha*sqrt(tCurr))) * cof,\n"
              << "                 within 1e-12 relative/absolute tolerance.\n";

    const double bindRadius = 1.0;
    const double r0 = bindRadius; // start at contact -> sep = 0
    const double Dtot = 10.0;
    const double cof = 1.0;

    const std::vector<double> alphaList { 0.5, 5.0, 50.0 };
    const std::vector<double> tList { 1e-8, 1e-6, 1e-3 };

    for (double alpha : alphaList) {
        for (double t : tList) {
            const double x = alpha * std::sqrt(t);
            const double expected = (1.0 - Faddeeva::erfcx(x)) * cof;
            const double p = passocF(r0, t, Dtot, bindRadius, alpha, cof);

            std::cerr << "    alpha = " << alpha << ", tCurr = " << t
                      << " -> p = " << p << " (expected " << expected << ")\n";
            EXPECT_NEAR(p, expected, 1e-12)
                << "Contact-start probability must match the closed form "
                   "(alpha = " << alpha << ", tCurr = " << t << ")";
        }
    }
}

// -----------------------------------------------------------------------------
// Test 3: General parameter sweep versus the independent reference formula.
//         This is the core correctness check of the direct (non-overflow) path.
// -----------------------------------------------------------------------------
void test_passocf_matches_reference_formula()
{
    std::cerr << "\n[TEST] test_passocf_matches_reference_formula\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF()\n"
              << "  Scenario:      sweep over r0, tCurr, Dtot, alpha, cof.\n"
              << "  Pass criteria: |p - reference| <= 1e-12 where reference =\n"
              << "                 [erfc(sep) - exp(-sep^2)*erfcx(ef1)] * f1.\n";

    const double bindRadius = 1.0;
    const std::vector<double> r0List { 1.0, 1.2, 2.0, 10.0 };
    const std::vector<double> tList { 1e-7, 1e-4, 1e-1 };
    const std::vector<double> dList { 0.5, 10.0, 100.0 };
    const std::vector<double> alphaList { 0.1, 2.0, 30.0 };
    const std::vector<double> cofList { 1.0, 0.25 };

    int checked = 0;
    for (double r0 : r0List) {
        for (double t : tList) {
            for (double D : dList) {
                for (double alpha : alphaList) {
                    for (double cof : cofList) {
                        const double p = passocF(r0, t, D, bindRadius, alpha, cof);
                        const double ref
                            = passocf_reference(r0, t, D, bindRadius, alpha, cof);

                        EXPECT_NEAR(p, ref, 1e-12)
                            << "Mismatch with reference for r0 = " << r0
                            << ", tCurr = " << t << ", Dtot = " << D
                            << ", alpha = " << alpha << ", cof = " << cof;
                        ++checked;
                    }
                }
            }
        }
    }
    std::cerr << "    Compared " << checked
              << " parameter combinations against the reference formula.\n";
}

// -----------------------------------------------------------------------------
// Test 4: The probability must never be negative and never exceed the
//         prefactor cof*bindRadius/r0 (the t -> infinity absorbing limit).
// -----------------------------------------------------------------------------
void test_passocf_bounds_are_respected()
{
    std::cerr << "\n[TEST] test_passocf_bounds_are_respected\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF()\n"
              << "  Scenario:      sweep over times / separations / rates.\n"
              << "  Pass criteria: 0 <= p <= cof*bindRadius/r0 and p is finite.\n";

    const double bindRadius = 2.0;
    const double cof = 1.0;
    const double Dtot = 25.0;

    const std::vector<double> r0List { 2.0, 2.5, 4.0, 20.0 };
    const std::vector<double> tList { 1e-9, 1e-5, 1e-2, 1.0 };
    const std::vector<double> alphaList { 0.01, 1.0, 100.0, 1.0e5 };

    for (double r0 : r0List) {
        const double upper = cof * bindRadius / r0;
        for (double t : tList) {
            for (double alpha : alphaList) {
                const double p = passocF(r0, t, Dtot, bindRadius, alpha, cof);

                EXPECT_TRUE(std::isfinite(p))
                    << "Probability must be finite (r0 = " << r0
                    << ", tCurr = " << t << ", alpha = " << alpha << ")";
                EXPECT_GE(p, 0.0)
                    << "Probability must be non-negative (r0 = " << r0
                    << ", tCurr = " << t << ", alpha = " << alpha << ")";
                EXPECT_LE(p, upper + 1e-12)
                    << "Probability must not exceed cof*bindRadius/r0 = " << upper
                    << " (r0 = " << r0 << ", tCurr = " << t
                    << ", alpha = " << alpha << ")";
            }
        }
    }
    std::cerr << "    All sampled probabilities stayed inside [0, cof*sigma/r0].\n";
}

// -----------------------------------------------------------------------------
// Test 5: Physical monotonicity in time - waiting longer can only increase the
//         chance that the pair has associated.
// -----------------------------------------------------------------------------
void test_passocf_increases_with_time()
{
    std::cerr << "\n[TEST] test_passocf_increases_with_time\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF()\n"
              << "  Scenario:      fixed geometry, increasing tCurr.\n"
              << "  Pass criteria: p(t_{i+1}) > p(t_i) for the sampled times.\n";

    const double bindRadius = 1.0;
    const double r0 = 1.3;
    const double Dtot = 10.0;
    const double alpha = 5.0;
    const double cof = 1.0;

    const std::vector<double> tList { 1e-7, 1e-6, 1e-5, 1e-4, 1e-3 };
    double prev = -1.0;
    for (double t : tList) {
        const double p = passocF(r0, t, Dtot, bindRadius, alpha, cof);
        std::cerr << "    tCurr = " << t << " -> p = " << p << '\n';
        EXPECT_GT(p, prev) << "Association probability should grow with time "
                              "(tCurr = " << t << ")";
        prev = p;
    }
}

// -----------------------------------------------------------------------------
// Test 6: Physical monotonicity in the starting separation - the farther apart
//         the pair starts, the smaller the association probability.
// -----------------------------------------------------------------------------
void test_passocf_decreases_with_separation()
{
    std::cerr << "\n[TEST] test_passocf_decreases_with_separation\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF()\n"
              << "  Scenario:      fixed time, increasing initial separation r0.\n"
              << "  Pass criteria: p(r0_{i+1}) < p(r0_i).\n";

    const double bindRadius = 1.0;
    const double tCurr = 1e-4;
    const double Dtot = 10.0;
    const double alpha = 5.0;
    const double cof = 1.0;

    const std::vector<double> r0List { 1.0, 1.05, 1.2, 1.5, 3.0 };
    double prev = std::numeric_limits<double>::infinity();
    for (double r0 : r0List) {
        const double p = passocF(r0, tCurr, Dtot, bindRadius, alpha, cof);
        std::cerr << "    r0 = " << r0 << " -> p = " << p << '\n';
        EXPECT_LT(p, prev) << "Association probability should fall off with the "
                              "initial separation (r0 = " << r0 << ")";
        prev = p;
    }

    // Very large separations must effectively give zero probability.
    const double pFar = passocF(1000.0, tCurr, Dtot, bindRadius, alpha, cof);
    std::cerr << "    r0 = 1000 (far field) -> p = " << pFar << '\n';
    EXPECT_NEAR(pFar, 0.0, 1e-12)
        << "A pair starting very far apart should have essentially zero "
           "association probability";
}

// -----------------------------------------------------------------------------
// Test 7: The `cof` argument is a pure multiplicative prefactor.
// -----------------------------------------------------------------------------
void test_passocf_is_linear_in_cof()
{
    std::cerr << "\n[TEST] test_passocf_is_linear_in_cof\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF()\n"
              << "  Scenario:      identical inputs, cof = 1 vs cof = 2 vs 0.\n"
              << "  Pass criteria: p(2*cof) == 2*p(cof) and p(cof = 0) == 0.\n";

    const double bindRadius = 1.0;
    const double r0 = 1.4;
    const double tCurr = 1e-5;
    const double Dtot = 8.0;
    const double alpha = 3.0;

    const double p1 = passocF(r0, tCurr, Dtot, bindRadius, alpha, 1.0);
    const double p2 = passocF(r0, tCurr, Dtot, bindRadius, alpha, 2.0);
    const double p0 = passocF(r0, tCurr, Dtot, bindRadius, alpha, 0.0);

    std::cerr << "    p(cof=1) = " << p1 << ", p(cof=2) = " << p2
              << ", p(cof=0) = " << p0 << '\n';

    EXPECT_NEAR(p2, 2.0 * p1, 1e-15 * std::fabs(2.0 * p1) + 1e-18)
        << "Doubling cof must double the returned probability";
    EXPECT_DOUBLE_EQ(p0, 0.0) << "cof = 0 must zero out the probability";
}

// -----------------------------------------------------------------------------
// Test 8: Explicitly exercise the Faddeeva (exp overflow) branch.
//         For a huge alpha the boundary is fully absorbing and the result must
//         approach erfc(sep) * cof * bindRadius / r0 from below.
// -----------------------------------------------------------------------------
void test_passocf_faddeeva_overflow_branch()
{
    std::cerr << "\n[TEST] test_passocf_faddeeva_overflow_branch\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF() - Faddeeva::w() branch\n"
              << "  Scenario:      alpha so large that exp(e1) overflows.\n"
              << "  Pass criteria: result is finite, matches the erfcx-based\n"
              << "                 reference, and lies just below the fully\n"
              << "                 absorbing limit erfc(sep)*cof*sigma/r0.\n";

    const double bindRadius = 1.0;
    const double r0 = 1.2;
    const double tCurr = 1.0;
    const double Dtot = 1.0;
    const double alpha = 1.0e6; // e1 ~ 1e12 -> exp() overflows to +inf
    const double cof = 1.0;

    // Confirm we really are on the overflow code path.
    const bool onFaddeevaBranch
        = passocf_hits_faddeeva_branch(r0, tCurr, Dtot, bindRadius, alpha);
    std::cerr << "    Faddeeva branch exercised: "
              << (onFaddeevaBranch ? "yes" : "no") << '\n';
    EXPECT_TRUE(onFaddeevaBranch)
        << "Test setup should force exp(e1) to overflow so the Faddeeva branch "
           "is taken";

    const double p = passocF(r0, tCurr, Dtot, bindRadius, alpha, cof);
    const double ref = passocf_reference(r0, tCurr, Dtot, bindRadius, alpha, cof);

    // Fully absorbing (alpha -> infinity) limit for this geometry.
    const double sep = (r0 - bindRadius) / std::sqrt(4.0 * Dtot * tCurr);
    const double absorbingLimit = std::erfc(sep) * cof * bindRadius / r0;

    std::cerr << "    p = " << p << ", reference = " << ref
              << ", absorbing limit = " << absorbingLimit << '\n';

    EXPECT_TRUE(std::isfinite(p))
        << "The Faddeeva branch must return a finite number, not inf/nan";
    EXPECT_NEAR(p, ref, 1e-12)
        << "Faddeeva branch must agree with the erfcx-based reference";
    EXPECT_LT(p, absorbingLimit)
        << "Result must stay strictly below the fully absorbing limit";
    EXPECT_NEAR(p, absorbingLimit, 1e-5)
        << "For alpha = 1e6 the result should be very close to the fully "
           "absorbing limit";
}

// -----------------------------------------------------------------------------
// Test 9: Continuity across the exp() overflow threshold.
//         Two alpha values straddling the overflow point (exp(709) is finite,
//         exp(710) is inf) must produce nearly the same answer - i.e. the two
//         code paths implement the same mathematics.
// -----------------------------------------------------------------------------
void test_passocf_branch_continuity_across_overflow()
{
    std::cerr << "\n[TEST] test_passocf_branch_continuity_across_overflow\n"
              << "  Source file:   src/reactions/passocF.cpp\n"
              << "  Function:      passocF() - direct vs Faddeeva branch\n"
              << "  Scenario:      sep = 0, tCurr = 1, alpha = 26 (exp finite)\n"
              << "                 and alpha = 27 (exp overflows).\n"
              << "  Pass criteria: both branches match the erfcx reference and\n"
              << "                 differ from each other by < 1e-2.\n";

    const double bindRadius = 1.0;
    const double r0 = bindRadius; // sep = 0 -> e1 = alpha^2 * tCurr
    const double tCurr = 1.0;
    const double Dtot = 1.0;
    const double cof = 1.0;

    const double alphaLow = 26.0; // 26^2 = 676  -> exp finite
    const double alphaHigh = 27.0; // 27^2 = 729 -> exp overflows

    const bool lowIsFaddeeva
        = passocf_hits_faddeeva_branch(r0, tCurr, Dtot, bindRadius, alphaLow);
    const bool highIsFaddeeva
        = passocf_hits_faddeeva_branch(r0, tCurr, Dtot, bindRadius, alphaHigh);

    std::cerr << "    alpha = " << alphaLow << " uses Faddeeva branch: "
              << (lowIsFaddeeva ? "yes" : "no") << '\n'
              << "    alpha = " << alphaHigh << " uses Faddeeva branch: "
              << (highIsFaddeeva ? "yes" : "no") << '\n';

    EXPECT_FALSE(lowIsFaddeeva)
        << "alpha = 26 should stay on the direct exp()/erfc() path";
    EXPECT_TRUE(highIsFaddeeva)
        << "alpha = 27 should overflow and switch to the Faddeeva path";

    const double pLow = passocF(r0, tCurr, Dtot, bindRadius, alphaLow, cof);
    const double pHigh = passocF(r0, tCurr, Dtot, bindRadius, alphaHigh, cof);
    const double refLow
        = passocf_reference(r0, tCurr, Dtot, bindRadius, alphaLow, cof);
    const double refHigh
        = passocf_reference(r0, tCurr, Dtot, bindRadius, alphaHigh, cof);

    std::cerr << "    p(alpha=26) = " << pLow << " (ref " << refLow << ")\n"
              << "    p(alpha=27) = " << pHigh << " (ref " << refHigh << ")\n";

    EXPECT_NEAR(pLow, refLow, 1e-12)
        << "Direct branch must match the reference at alpha = 26";
    EXPECT_NEAR(pHigh, refHigh, 1e-12)
        << "Faddeeva branch must match the reference at alpha = 27";

    // The exact values differ slightly because alpha differs, but the jump must
    // be small - a broken branch would produce an O(1) discontinuity.
    EXPECT_NEAR(pHigh, pLow, 1e-2)
        << "The two branches must be continuous across the overflow threshold";
    EXPECT_GT(pHigh, pLow)
        << "A larger intrinsic rate (alpha) must give a larger probability";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario runs in its own TEST so that a failure in
// one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(PassocFTest, ZeroAlphaGivesZeroProbability) { test_passocf_zero_alpha_gives_zero_probability(); }
TEST(PassocFTest, AtContactMatchesClosedForm) { test_passocf_at_contact_matches_closed_form(); }
TEST(PassocFTest, MatchesReferenceFormula) { test_passocf_matches_reference_formula(); }
TEST(PassocFTest, BoundsAreRespected) { test_passocf_bounds_are_respected(); }
TEST(PassocFTest, IncreasesWithTime) { test_passocf_increases_with_time(); }
TEST(PassocFTest, DecreasesWithSeparation) { test_passocf_decreases_with_separation(); }
TEST(PassocFTest, IsLinearInCof) { test_passocf_is_linear_in_cof(); }
TEST(PassocFTest, FaddeevaOverflowBranch) { test_passocf_faddeeva_overflow_branch(); }
TEST(PassocFTest, BranchContinuityAcrossOverflow) { test_passocf_branch_continuity_across_overflow(); }