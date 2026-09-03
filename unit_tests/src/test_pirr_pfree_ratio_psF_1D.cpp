/*! \file test_pirr_pfree_ratio_psF_1D.cpp
 *
 * ### Unit tests for src/reactions/pirr_pfree_ratio_psF_1D.cpp
 *
 * The function under test is
 *
 *     double pirr_pfree_ratio_psF_1D(double rCurr, double r0, double tCurr,
 *                                    double Dtot, double bindrad, double ka,
 *                                    double ps_prev);
 *
 * It returns the reweighting ratio  p_irr / (ps_prev * p_free_norm)  used by the
 * 1D (fiber / rod) bimolecular reaction machinery.  The implementation contains
 * two completely separate analytic branches:
 *
 *   * bindrad == 0  -> "semi-permeable"/delta-sink solution
 *          p_free_norm = exp(-(r-r0)^2/4Dt) / sqrt(4 pi D t)
 *          p_irr       = p_free_norm
 *                        - (ka/D) exp(-(r+r0)^2/4Dt) erfcx((r+r0)/sqrt(4Dt) + ka sqrt(t/D))
 *
 *   * bindrad  > 0  -> radiation boundary condition at x = bindrad
 *          I_free      = sqrt(pi D t) [ erf((r0-a)/sqrt(4Dt)) + erf((r0+a)/sqrt(4Dt)) ]
 *          p_free_norm = [exp(-(r-r0)^2/4Dt) - exp(-(r+r0)^2/4Dt)] / I_free
 *          p_irr       = [exp(-(r-r0)^2/4Dt) + exp(-(r+r0-2a)^2/4Dt)] / sqrt(4 pi D t)
 *                        - (ka/D) exp(-(r+r0-2a)^2/4Dt) erfcx((r+r0-2a)/sqrt(4Dt) + ka sqrt(t/D))
 *
 * The tests below verify the implementation against
 *   (a) independently coded closed forms that use only <cmath> (erf / erfc / exp)
 *       instead of the Faddeeva erfcx routine used inside the function, and
 *   (b) analytic limits (no reaction, strong absorption, large separation)
 *       that can be derived on paper without re-implementing the formula.
 *
 * The function is a pure numeric routine: it has no side effects, does not use
 * the global GSL RNG, and never calls exit()/abort(), so every case below is
 * safe to run inside the shared gtest binary.
 */

#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <limits>

// -----------------------------------------------------------------------------
// Local helpers (internal linkage; names additionally prefixed to stay unique).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Scaled complementary error function, computed independently of the
 *         Faddeeva library used inside the function under test.
 *
 * erfcx(x) = exp(x^2) * erfc(x).  This naive form is only valid while
 * exp(x^2) does not overflow (x < ~26), which is true for every argument used
 * in the closed-form comparison tests below.
 */
double pirr1d_erfcx_naive(double x) { return std::exp(x * x) * std::erfc(x); }

/*! \brief Independent reference implementation for the bindrad == 0 branch. */
double pirr1d_ref_zero_bindrad(
    double rCurr, double r0, double tCurr, double Dtot, double ka, double ps_prev)
{
    const double fDt = 4.0 * Dtot * tCurr;          // 4 D t
    const double sqFDt = std::sqrt(fDt);            // sqrt(4 D t)
    const double sqPiFDt = std::sqrt(4.0 * M_PI * Dtot * tCurr);
    const double dxMinus = rCurr - r0;
    const double dxPlus = rCurr + r0;

    const double e1 = std::exp(-dxMinus * dxMinus / fDt);
    const double e2 = std::exp(-dxPlus * dxPlus / fDt);
    const double arg = dxPlus / sqFDt + ka * std::sqrt(tCurr / Dtot);

    const double pFree = e1 / sqPiFDt;
    const double pIrr = pFree - (ka / Dtot) * e2 * pirr1d_erfcx_naive(arg);

    return pIrr / ps_prev / pFree;
}

/*! \brief Independent reference implementation for the bindrad > 0 branch. */
double pirr1d_ref_finite_bindrad(double rCurr, double r0, double tCurr, double Dtot,
    double bindrad, double ka, double ps_prev)
{
    const double fDt = 4.0 * Dtot * tCurr;
    const double sqFDt = std::sqrt(fDt);
    const double sqPiFDt = std::sqrt(4.0 * M_PI * Dtot * tCurr);

    const double dxMinus = rCurr - r0;
    const double dxPlus = rCurr + r0;
    const double dxSigma = rCurr + r0 - 2.0 * bindrad;
    const double xSigmaMinus = r0 - bindrad;
    const double xSigmaPlus = r0 + bindrad;

    const double e1 = std::exp(-dxMinus * dxMinus / fDt);
    const double e2 = std::exp(-dxPlus * dxPlus / fDt);
    const double e3 = std::exp(-dxSigma * dxSigma / fDt);

    const double iFree = std::sqrt(M_PI * Dtot * tCurr)
        * (std::erf(xSigmaMinus / sqFDt) + std::erf(xSigmaPlus / sqFDt));

    const double pFree = (e1 - e2) / iFree;
    const double sep = dxSigma / sqFDt + ka * std::sqrt(tCurr / Dtot);
    const double pIrr = (e1 + e3) / sqPiFDt - (ka / Dtot) * e3 * pirr1d_erfcx_naive(sep);

    return pIrr / ps_prev / pFree;
}

/*! \brief Relative-plus-absolute tolerance helper for floating point compares. */
double pirr1d_tol(double expected, double relTol = 1e-10)
{
    return std::abs(expected) * relTol + 1e-14;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: bindrad == 0 with no reactivity (ka == 0).
//
// With ka == 0 the correction term is multiplied by exactly zero, so p_irr and
// p_free_norm are the identical expression and the ratio must collapse to
// 1 / ps_prev.
// -----------------------------------------------------------------------------
void test_pirr1d_zero_bindrad_no_reaction()
{
    std::cerr << "\n[TEST] test_pirr1d_zero_bindrad_no_reaction\n"
              << "  Source file: src/reactions/pirr_pfree_ratio_psF_1D.cpp\n"
              << "  Function:    pirr_pfree_ratio_psF_1D (bindrad == 0 branch)\n"
              << "  Scenario:    ka = 0, so p_irr == p_free_norm exactly.\n"
              << "  Pass:        ratio == 1 / ps_prev.\n";

    // ps_prev = 1 -> ratio must be exactly unity.
    double ratio = pirr_pfree_ratio_psF_1D(/*rCurr*/ 2.0, /*r0*/ 1.0, /*tCurr*/ 1.0,
        /*Dtot*/ 1.0, /*bindrad*/ 0.0, /*ka*/ 0.0, /*ps_prev*/ 1.0);
    std::cerr << "    ka=0, ps_prev=1.0 -> ratio = " << ratio << " (expect 1)\n";
    EXPECT_DOUBLE_EQ(ratio, 1.0)
        << "With no reactivity the survival-normalised ratio must be unity";

    // ps_prev = 0.5 -> ratio must be 2.
    ratio = pirr_pfree_ratio_psF_1D(2.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.5);
    std::cerr << "    ka=0, ps_prev=0.5 -> ratio = " << ratio << " (expect 2)\n";
    EXPECT_DOUBLE_EQ(ratio, 2.0) << "Ratio must scale as 1/ps_prev";

    // ps_prev = 2.0 -> ratio must be 0.5.
    ratio = pirr_pfree_ratio_psF_1D(2.0, 1.0, 1.0, 1.0, 0.0, 0.0, 2.0);
    std::cerr << "    ka=0, ps_prev=2.0 -> ratio = " << ratio << " (expect 0.5)\n";
    EXPECT_DOUBLE_EQ(ratio, 0.5) << "Ratio must scale as 1/ps_prev";
}

// -----------------------------------------------------------------------------
// Test 2: bindrad == 0 compared against an independently coded closed form.
//
// The reference uses exp(x^2)*erfc(x) from <cmath> rather than the Faddeeva
// erfcx routine used inside the function, so agreement validates both the
// algebra and the erfcx call.
// -----------------------------------------------------------------------------
void test_pirr1d_zero_bindrad_matches_closed_form()
{
    std::cerr << "\n[TEST] test_pirr1d_zero_bindrad_matches_closed_form\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (bindrad == 0 branch)\n"
              << "  Scenario:  several (r, r0, t, D, ka, ps) sets with ka > 0.\n"
              << "  Pass:      matches exp(x^2)erfc(x)-based reference to 1e-10 rel.\n";

    struct Case {
        double rCurr, r0, tCurr, Dtot, ka, ps_prev;
    };
    // Arguments of erfcx stay below ~10 so the naive reference cannot overflow.
    const Case cases[] = {
        { 2.0, 1.0, 1.0, 1.0, 1.0, 1.0 },
        { 1.0, 1.0, 1.0, 1.0, 0.5, 1.0 },
        { 3.0, 2.0, 0.5, 2.0, 2.0, 0.9 },
        { 0.5, 0.4, 2.0, 0.25, 0.1, 0.75 },
        { 5.0, 1.0, 1.0, 4.0, 3.0, 1.0 },
    };

    for (const auto& c : cases) {
        const double got = pirr_pfree_ratio_psF_1D(
            c.rCurr, c.r0, c.tCurr, c.Dtot, /*bindrad*/ 0.0, c.ka, c.ps_prev);
        const double want
            = pirr1d_ref_zero_bindrad(c.rCurr, c.r0, c.tCurr, c.Dtot, c.ka, c.ps_prev);

        std::cerr << "    r=" << c.rCurr << " r0=" << c.r0 << " t=" << c.tCurr
                  << " D=" << c.Dtot << " ka=" << c.ka << " ps=" << c.ps_prev
                  << " -> got " << got << ", reference " << want << '\n';

        EXPECT_NEAR(got, want, pirr1d_tol(want))
            << "bindrad==0 branch disagrees with the independently coded closed form";
    }
}

// -----------------------------------------------------------------------------
// Test 3: bindrad == 0, monotone dependence on the intrinsic rate ka.
//
// The correction term  (ka/D) * exp(-(r+r0)^2/4Dt) * erfcx(...)  is strictly
// positive for ka > 0, so the ratio must be strictly below 1/ps_prev and must
// decrease as ka grows (more reactivity -> more depletion near the sink).
// -----------------------------------------------------------------------------
void test_pirr1d_zero_bindrad_ka_monotonic()
{
    std::cerr << "\n[TEST] test_pirr1d_zero_bindrad_ka_monotonic\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (bindrad == 0 branch)\n"
              << "  Scenario:  sweep ka upward at fixed geometry.\n"
              << "  Pass:      ratio < 1/ps_prev and strictly decreasing in ka.\n";

    const double rCurr = 2.0, r0 = 1.0, tCurr = 1.0, Dtot = 1.0, ps_prev = 1.0;
    const double kaList[] = { 0.0, 0.1, 0.5, 1.0, 5.0, 20.0 };

    double previous = std::numeric_limits<double>::infinity();
    for (double ka : kaList) {
        const double ratio
            = pirr_pfree_ratio_psF_1D(rCurr, r0, tCurr, Dtot, 0.0, ka, ps_prev);
        std::cerr << "    ka = " << ka << " -> ratio = " << ratio << '\n';

        // Every reactive case must sit strictly below the non-reactive value.
        if (ka > 0.0) {
            EXPECT_LT(ratio, 1.0 / ps_prev)
                << "Reactive sink must deplete the irreversible density below p_free";
        }
        EXPECT_LT(ratio, previous)
            << "Ratio must decrease monotonically with increasing intrinsic rate ka";
        previous = ratio;
    }
}

// -----------------------------------------------------------------------------
// Test 4: bindrad == 0, strong-absorption (ka -> infinity) analytic limit.
//
// For ka -> infinity, erfcx(x) ~ 1/(x sqrt(pi)) with x ~ ka sqrt(t/D), hence the
// correction tends to exp(-(r+r0)^2/4Dt)/sqrt(pi D t) and
//
//        ratio -> 1 - 2 exp(-(r+r0)^2/4Dt) / exp(-(r-r0)^2/4Dt)
//
// which is derived on paper, independent of the implementation.
// -----------------------------------------------------------------------------
void test_pirr1d_zero_bindrad_strong_absorption_limit()
{
    std::cerr << "\n[TEST] test_pirr1d_zero_bindrad_strong_absorption_limit\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (bindrad == 0 branch)\n"
              << "  Scenario:  ka = 1e8 (effectively perfect absorption).\n"
              << "  Pass:      ratio -> 1 - 2 e2/e1 (analytic limit).\n";

    const double rCurr = 2.0, r0 = 1.0, tCurr = 1.0, Dtot = 1.0, ps_prev = 1.0;
    const double fDt = 4.0 * Dtot * tCurr;
    const double e1 = std::exp(-(rCurr - r0) * (rCurr - r0) / fDt);
    const double e2 = std::exp(-(rCurr + r0) * (rCurr + r0) / fDt);
    const double expected = 1.0 - 2.0 * e2 / e1;

    const double ratio
        = pirr_pfree_ratio_psF_1D(rCurr, r0, tCurr, Dtot, 0.0, 1.0e8, ps_prev);

    std::cerr << "    ratio = " << ratio << ", analytic limit = " << expected << '\n';

    // Leading correction to the asymptotic expansion is O((r+r0)/(sqrt(4Dt) ka)),
    // i.e. ~1e-8 here, so a 1e-6 window is a comfortable pass criterion.
    EXPECT_NEAR(ratio, expected, 1.0e-6)
        << "Strong-absorption limit of the delta-sink branch is incorrect";
}

// -----------------------------------------------------------------------------
// Test 5: bindrad == 0, large-separation limit.
//
// When the pair is many diffusion lengths apart the sink is invisible, so the
// irreversible and free propagators coincide and ratio -> 1/ps_prev.
// -----------------------------------------------------------------------------
void test_pirr1d_zero_bindrad_large_separation()
{
    std::cerr << "\n[TEST] test_pirr1d_zero_bindrad_large_separation\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (bindrad == 0 branch)\n"
              << "  Scenario:  rCurr = 40 with sqrt(4Dt) = 2 (20 diffusion lengths).\n"
              << "  Pass:      ratio -> 1/ps_prev to within 1e-12.\n";

    const double ratio = pirr_pfree_ratio_psF_1D(
        /*rCurr*/ 40.0, /*r0*/ 1.0, /*tCurr*/ 1.0, /*Dtot*/ 1.0,
        /*bindrad*/ 0.0, /*ka*/ 1.0, /*ps_prev*/ 1.0);

    std::cerr << "    ratio = " << ratio << " (expect ~1)\n";
    EXPECT_NEAR(ratio, 1.0, 1.0e-12)
        << "Far from the sink the reweighting factor must reduce to 1/ps_prev";

    // Approach to the limit must be monotone from below in rCurr.
    double previous = -std::numeric_limits<double>::infinity();
    for (double r : { 1.5, 2.0, 4.0, 8.0, 16.0 }) {
        const double val = pirr_pfree_ratio_psF_1D(r, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0);
        std::cerr << "    rCurr = " << r << " -> ratio = " << val << '\n';
        EXPECT_GT(val, previous) << "Ratio should rise toward 1 as the pair separates";
        EXPECT_LT(val, 1.0) << "Ratio must stay below the non-reactive value of 1";
        previous = val;
    }
}

// -----------------------------------------------------------------------------
// Test 6: bindrad > 0 compared against an independently coded closed form.
// -----------------------------------------------------------------------------
void test_pirr1d_finite_bindrad_matches_closed_form()
{
    std::cerr << "\n[TEST] test_pirr1d_finite_bindrad_matches_closed_form\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (bindrad > 0 branch)\n"
              << "  Scenario:  radiation boundary condition at x = bindrad.\n"
              << "  Pass:      matches erf/erfc-based reference to 1e-10 rel.\n";

    struct Case {
        double rCurr, r0, tCurr, Dtot, bindrad, ka, ps_prev;
    };
    const Case cases[] = {
        { 3.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0 },   // generic reactive case
        { 3.0, 2.0, 1.0, 1.0, 1.0, 0.0, 1.0 },   // ka = 0 -> no erfcx contribution
        { 1.5, 1.2, 0.5, 2.0, 1.0, 3.0, 0.8 },   // faster diffusion, shorter time
        { 6.0, 2.5, 2.0, 0.5, 2.0, 0.25, 0.95 }, // larger binding radius
        { 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0 },   // rCurr == r0
    };

    for (const auto& c : cases) {
        const double got = pirr_pfree_ratio_psF_1D(
            c.rCurr, c.r0, c.tCurr, c.Dtot, c.bindrad, c.ka, c.ps_prev);
        const double want = pirr1d_ref_finite_bindrad(
            c.rCurr, c.r0, c.tCurr, c.Dtot, c.bindrad, c.ka, c.ps_prev);

        std::cerr << "    r=" << c.rCurr << " r0=" << c.r0 << " t=" << c.tCurr
                  << " D=" << c.Dtot << " a=" << c.bindrad << " ka=" << c.ka
                  << " ps=" << c.ps_prev << " -> got " << got << ", reference " << want
                  << '\n';

        EXPECT_NEAR(got, want, pirr1d_tol(want))
            << "bindrad>0 branch disagrees with the independently coded closed form";
    }
}

// -----------------------------------------------------------------------------
// Test 7: bindrad > 0, large-separation limit.
//
// For r >> sqrt(4Dt) both image terms vanish relative to exp(-(r-r0)^2/4Dt) and
//
//        ratio -> I_free / sqrt(4 pi D t) / ps_prev
//               = [erf((r0-a)/sqrt(4Dt)) + erf((r0+a)/sqrt(4Dt))] / 2 / ps_prev
//
// Note that, unlike the delta-sink branch, this is NOT 1/ps_prev: the radiation
// branch normalises p_free over the half-line x > bindrad.
// -----------------------------------------------------------------------------
void test_pirr1d_finite_bindrad_large_separation()
{
    std::cerr << "\n[TEST] test_pirr1d_finite_bindrad_large_separation\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (bindrad > 0 branch)\n"
              << "  Scenario:  rCurr = 40, r0 = 2, bindrad = 1, sqrt(4Dt) = 2.\n"
              << "  Pass:      ratio -> [erf((r0-a)/sqrt(4Dt))+erf((r0+a)/sqrt(4Dt))]/2.\n";

    const double rCurr = 40.0, r0 = 2.0, tCurr = 1.0, Dtot = 1.0;
    const double bindrad = 1.0, ka = 1.0, ps_prev = 1.0;

    const double sqFDt = std::sqrt(4.0 * Dtot * tCurr);
    const double expected
        = 0.5 * (std::erf((r0 - bindrad) / sqFDt) + std::erf((r0 + bindrad) / sqFDt))
        / ps_prev;

    const double ratio = pirr_pfree_ratio_psF_1D(
        rCurr, r0, tCurr, Dtot, bindrad, ka, ps_prev);

    std::cerr << "    ratio = " << ratio << ", analytic limit = " << expected << '\n';
    EXPECT_NEAR(ratio, expected, 1.0e-10)
        << "Large-separation limit of the radiation-BC branch is incorrect";
}

// -----------------------------------------------------------------------------
// Test 8: both branches scale exactly as 1 / ps_prev.
//
// ps_prev only enters as a division, so ratio(ps) * ps must be independent of
// ps for fixed geometry and rates.
// -----------------------------------------------------------------------------
void test_pirr1d_ps_prev_scaling()
{
    std::cerr << "\n[TEST] test_pirr1d_ps_prev_scaling\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (both branches)\n"
              << "  Scenario:  vary only ps_prev.\n"
              << "  Pass:      ratio * ps_prev is constant.\n";

    // --- delta-sink branch -------------------------------------------------
    const double baseZero
        = pirr_pfree_ratio_psF_1D(2.0, 1.0, 1.0, 1.0, /*bindrad*/ 0.0, 1.0, 1.0);
    for (double ps : { 0.25, 0.5, 0.75, 0.99 }) {
        const double val
            = pirr_pfree_ratio_psF_1D(2.0, 1.0, 1.0, 1.0, 0.0, 1.0, ps);
        std::cerr << "    [bindrad=0] ps_prev = " << ps << " -> ratio = " << val
                  << ", ratio*ps = " << val * ps << " (base " << baseZero << ")\n";
        EXPECT_NEAR(val * ps, baseZero, pirr1d_tol(baseZero))
            << "Delta-sink branch must scale exactly as 1/ps_prev";
    }

    // --- radiation-BC branch ----------------------------------------------
    const double baseFinite
        = pirr_pfree_ratio_psF_1D(3.0, 2.0, 1.0, 1.0, /*bindrad*/ 1.0, 1.0, 1.0);
    for (double ps : { 0.25, 0.5, 0.75, 0.99 }) {
        const double val
            = pirr_pfree_ratio_psF_1D(3.0, 2.0, 1.0, 1.0, 1.0, 1.0, ps);
        std::cerr << "    [bindrad=1] ps_prev = " << ps << " -> ratio = " << val
                  << ", ratio*ps = " << val * ps << " (base " << baseFinite << ")\n";
        EXPECT_NEAR(val * ps, baseFinite, pirr1d_tol(baseFinite))
            << "Radiation-BC branch must scale exactly as 1/ps_prev";
    }
}

// -----------------------------------------------------------------------------
// Test 9: the two branches are genuinely different models.
//
// This documents (rather than "fixes") the behaviour: bindrad == 0 selects a
// delta-sink solution whose free propagator is normalised over the whole line,
// while bindrad = 1e-9 selects the radiation-BC solution whose free propagator
// is normalised over the half-line.  The two therefore do NOT agree in the
// bindrad -> 0 limit, and callers must not rely on continuity there.
// -----------------------------------------------------------------------------
void test_pirr1d_branch_discontinuity_at_zero_bindrad()
{
    std::cerr << "\n[TEST] test_pirr1d_branch_discontinuity_at_zero_bindrad\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (branch selection)\n"
              << "  Scenario:  bindrad = 0 versus bindrad = 1e-9, all else equal.\n"
              << "  Pass:      the two branches differ noticeably (different PDEs).\n";

    const double atZero
        = pirr_pfree_ratio_psF_1D(2.0, 1.0, 1.0, 1.0, /*bindrad*/ 0.0, 1.0, 1.0);
    const double nearZero
        = pirr_pfree_ratio_psF_1D(2.0, 1.0, 1.0, 1.0, /*bindrad*/ 1.0e-9, 1.0, 1.0);

    std::cerr << "    bindrad = 0    -> ratio = " << atZero << '\n';
    std::cerr << "    bindrad = 1e-9 -> ratio = " << nearZero << '\n';
    std::cerr << "    |difference|   = " << std::abs(atZero - nearZero) << '\n';

    // Both must be finite numbers ...
    EXPECT_TRUE(std::isfinite(atZero)) << "bindrad==0 branch produced a non-finite value";
    EXPECT_TRUE(std::isfinite(nearZero)) << "bindrad>0 branch produced a non-finite value";

    // ... but they are not the same, because the free-propagator normalisation
    // differs by the image term exp(-(r+r0)^2/4Dt).
    EXPECT_GT(std::abs(atZero - nearZero), 0.1)
        << "The two analytic branches are expected to differ; if they now agree "
           "the normalisation convention has changed";
}

// -----------------------------------------------------------------------------
// Test 10: broad parameter sweep - the routine must never produce NaN or Inf
//          for physically sensible inputs.
// -----------------------------------------------------------------------------
void test_pirr1d_finiteness_sweep()
{
    std::cerr << "\n[TEST] test_pirr1d_finiteness_sweep\n"
              << "  Function:  pirr_pfree_ratio_psF_1D (both branches)\n"
              << "  Scenario:  sweep separation, time, diffusion and rate.\n"
              << "  Pass:      every returned value is finite (no NaN/Inf).\n";

    int evaluated = 0;
    for (double bindrad : { 0.0, 1.0 }) {
        for (double r0 : { 1.1, 2.0, 5.0 }) {
            for (double rCurr : { 1.05, 2.0, 6.0, 20.0 }) {
                for (double tCurr : { 1.0e-3, 1.0, 10.0 }) {
                    for (double Dtot : { 0.05, 1.0, 25.0 }) {
                        for (double ka : { 0.0, 1.0, 100.0 }) {
                            // Keep rCurr and r0 outside the binding radius, which
                            // is the regime the caller always supplies.
                            if (bindrad > 0.0 && (rCurr < bindrad || r0 < bindrad))
                                continue;

                            const double val = pirr_pfree_ratio_psF_1D(
                                rCurr, r0, tCurr, Dtot, bindrad, ka, /*ps_prev*/ 1.0);
                            ++evaluated;

                            EXPECT_TRUE(std::isfinite(val))
                                << "Non-finite ratio for r=" << rCurr << " r0=" << r0
                                << " t=" << tCurr << " D=" << Dtot
                                << " a=" << bindrad << " ka=" << ka;
                        }
                    }
                }
            }
        }
    }

    std::cerr << "    evaluated " << evaluated << " parameter combinations, "
              << "all finite.\n";
    EXPECT_GT(evaluated, 0) << "Sweep should have evaluated at least one combination";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is executed inside its own TEST so a
// failure in one scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(PirrPfreeRatioPsF1D, ZeroBindradNoReaction) { test_pirr1d_zero_bindrad_no_reaction(); }
TEST(PirrPfreeRatioPsF1D, ZeroBindradMatchesClosedForm) { test_pirr1d_zero_bindrad_matches_closed_form(); }
TEST(PirrPfreeRatioPsF1D, ZeroBindradKaMonotonic) { test_pirr1d_zero_bindrad_ka_monotonic(); }
TEST(PirrPfreeRatioPsF1D, ZeroBindradStrongAbsorptionLimit) { test_pirr1d_zero_bindrad_strong_absorption_limit(); }
TEST(PirrPfreeRatioPsF1D, ZeroBindradLargeSeparation) { test_pirr1d_zero_bindrad_large_separation(); }
TEST(PirrPfreeRatioPsF1D, FiniteBindradMatchesClosedForm) { test_pirr1d_finite_bindrad_matches_closed_form(); }
TEST(PirrPfreeRatioPsF1D, FiniteBindradLargeSeparation) { test_pirr1d_finite_bindrad_large_separation(); }
TEST(PirrPfreeRatioPsF1D, PsPrevScaling) { test_pirr1d_ps_prev_scaling(); }
TEST(PirrPfreeRatioPsF1D, BranchDiscontinuityAtZeroBindrad) { test_pirr1d_branch_discontinuity_at_zero_bindrad(); }
TEST(PirrPfreeRatioPsF1D, FinitenessSweep) { test_pirr1d_finiteness_sweep(); }
