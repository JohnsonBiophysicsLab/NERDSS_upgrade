/*! \file test_DDpirr_pfree_ratio_ps.cpp
 *
 * ### Unit tests for src/reactions/DDpirr_pfree_ratio_ps.cpp
 *
 * The single function under test is
 *
 * \code
 * double DDpirr_pfree_ratio_ps(gsl_matrix* pirMatrix, gsl_matrix* survMatrix,
 *                              gsl_matrix* normMatrix, double r, double Dtot,
 *                              double deltaT, double r0, double ps_prev,
 *                              double rTol, double bindRadius);
 * \endcode
 *
 * What the implementation does, line by line:
 *   1. It evaluates the analytic 2D free-diffusion Green's function
 *
 *          p_free(r, dt | r0) = 1/(4 pi D dt) * exp(-(r^2 + r0^2)/(4 D dt))
 *                                            * I0( r*r0 / (2 D dt) )
 *
 *      but written in the numerically stable "scaled Bessel" form
 *      temp2 * temp3 * I0_scaled(temp), where
 *          temp  = r*r0/(2 D dt)
 *          temp2 = 1/(4 pi D dt)
 *          temp3 = exp(temp - (r0^2 + r^2)/(4 D dt))
 *      Note temp - (r^2+r0^2)/(4Ddt) == -(r-r0)^2/(4Ddt), and
 *      I0_scaled(x) == exp(-x) * I0(x), so the two forms are algebraically
 *      identical.  The tests below exploit this to check the formula with an
 *      *independent* implementation that uses the unscaled gsl_sf_bessel_I0().
 *   2. It looks up the normalisation and the irreversible pair distribution
 *      from the precomputed 2D tables via get_prevNorm() and calc_pirr(),
 *      always using RstepSize = sqrt(Dtot*deltaT)/50.
 *   3. It returns pirr / ((pFree/pNorm) * ps_prev), except that it short
 *      circuits to exactly 1.0 whenever |pirr - (pFree/pNorm)*ps_prev| < rTol.
 *
 * Test strategy
 * -------------
 * The lookup helpers get_prevNorm() and calc_pirr() live in another
 * translation unit; rather than guessing at their internals the tests simply
 * *call them* to build the expected value, exactly as the function under test
 * does.  That still gives real coverage of:
 *   - the analytic p_free expression (checked against an independent
 *     closed-form evaluation using the unscaled Bessel function),
 *   - the RstepSize definition sqrt(Dtot*deltaT)/50 that is handed to the
 *     lookup helpers,
 *   - the ordering/plumbing of the arguments,
 *   - the rTol short-circuit branch and its threshold,
 *   - the 1/ps_prev scaling of the returned ratio.
 *
 * The lookup tables are synthetic (smooth, strictly positive) matrices that
 * are far larger than any index the helpers can generate for the chosen
 * physical parameters, so the reads stay in range.  As a belt-and-braces
 * measure the GSL error handler is switched off for the duration of a test
 * (and restored afterwards) so that an unexpected out-of-range read would
 * produce a failing number instead of aborting the whole test binary.
 */

#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Fixed physical parameters shared by every test in this file.
// With Dtot*deltaT = 1 the internal RstepSize is sqrt(1)/50 = 0.02 nm, so the
// largest table index reachable from r,r0 <= 1.5 nm is (1.5-1.0)/0.02 = 25.
// The 512x512 tables below therefore leave an enormous safety margin.
// -----------------------------------------------------------------------------
constexpr double kDdpBindRadius = 1.0; //!< sigma, nm
constexpr double kDdpDtot = 1.0;       //!< total diffusion constant, nm^2/us
constexpr double kDdpDeltaT = 1.0;     //!< time step, us
constexpr size_t kDdpDim = 512;        //!< rows == cols of each synthetic table

/*! \brief Small RAII helper that silences the GSL error handler for a scope.
 *
 * Without this, a hypothetical out-of-range gsl_matrix_get() inside one of the
 * lookup helpers would abort() the entire gtest binary.  With the handler off
 * the read simply returns 0 and the assertions below fail loudly instead.
 */
struct DdpErrorHandlerGuard {
    gsl_error_handler_t* previous { nullptr };
    DdpErrorHandlerGuard() { previous = gsl_set_error_handler_off(); }
    ~DdpErrorHandlerGuard() { gsl_set_error_handler(previous); }
};

/*! \brief Container for the three lookup tables the function needs. */
struct DdpMatrixSet {
    gsl_matrix* surv { nullptr };
    gsl_matrix* norm { nullptr };
    gsl_matrix* pir { nullptr };
};

/*! \brief Allocate and fill the three synthetic lookup tables.
 *
 * The numbers are arbitrary but deliberately smooth, finite and strictly
 * positive so that no division by zero can occur downstream.  Both the code
 * under test and the reference computation see exactly the same tables, so the
 * actual values never need to be physically meaningful.
 */
DdpMatrixSet ddp_alloc_tables()
{
    DdpMatrixSet tables;
    tables.surv = gsl_matrix_alloc(kDdpDim, kDdpDim);
    tables.norm = gsl_matrix_alloc(kDdpDim, kDdpDim);
    tables.pir = gsl_matrix_alloc(kDdpDim, kDdpDim);

    for (size_t i = 0; i < kDdpDim; ++i) {
        for (size_t j = 0; j < kDdpDim; ++j) {
            const double fi = static_cast<double>(i);
            const double fj = static_cast<double>(j);
            // normalisation table: O(1), never zero
            gsl_matrix_set(tables.norm, i, j, 1.0 + 0.010 * fi + 0.001 * fj);
            // survival table: a decaying probability strictly inside (0, 1)
            gsl_matrix_set(tables.surv, i, j, 0.9 / (1.0 + 0.001 * (fi + fj)));
            // irreversible pair distribution: small and positive
            gsl_matrix_set(tables.pir, i, j, 0.05 + 0.0007 * fi + 0.0003 * fj);
        }
    }
    return tables;
}

/*! \brief Release the three lookup tables. */
void ddp_free_tables(DdpMatrixSet& tables)
{
    if (tables.surv != nullptr)
        gsl_matrix_free(tables.surv);
    if (tables.norm != nullptr)
        gsl_matrix_free(tables.norm);
    if (tables.pir != nullptr)
        gsl_matrix_free(tables.pir);
    tables.surv = nullptr;
    tables.norm = nullptr;
    tables.pir = nullptr;
}

/*! \brief Reproduce the RstepSize definition used inside the function.
 *
 * If somebody ever edits the "/50" in the source, the reference values
 * computed with this helper will disagree with the function's result because
 * the lookup helpers will be handed a different step size.
 */
double ddp_step_size(double Dtot, double deltaT)
{
    return std::sqrt(Dtot * deltaT) / 50.0;
}

/*! \brief Independent, closed-form evaluation of the 2D free Green's function.
 *
 * Uses the *unscaled* modified Bessel function I0 and the plain exponential,
 *      p_free = 1/(4 pi D dt) * exp(-(r^2+r0^2)/(4 D dt)) * I0(r r0/(2 D dt)),
 * i.e. a different code path from the source's temp/temp2/temp3 +
 * gsl_sf_bessel_I0_scaled() arrangement.  Agreement between the two is a real
 * check of the algebra in the function under test.
 */
double ddp_reference_pfree(double r, double r0, double Dtot, double deltaT)
{
    const double fourDt = 4.0 * Dtot * deltaT;
    const double prefactor = 1.0 / (M_PI * fourDt);
    const double gaussian = std::exp(-(r * r + r0 * r0) / fourDt);
    const double bessel = gsl_sf_bessel_I0(r * r0 / (2.0 * Dtot * deltaT));
    return prefactor * gaussian * bessel;
}

/*! \brief Rebuild the whole expected return value of DDpirr_pfree_ratio_ps.
 *
 * Uses the independent p_free above plus the very same lookup helpers the
 * source calls, so that only the arithmetic/branching of the function under
 * test is being verified.
 */
double ddp_reference_ratio(const DdpMatrixSet& tables, double r, double Dtot, double deltaT,
    double r0, double ps_prev, double rTol, double bindRadius)
{
    const double step = ddp_step_size(Dtot, deltaT);
    const double pFree = ddp_reference_pfree(r, r0, Dtot, deltaT);
    const double pNorm = get_prevNorm(tables.norm, step, r0, bindRadius);
    const double pirr = calc_pirr(tables.pir, tables.surv, step, r, r0, bindRadius);

    const double pfreeN = pFree / pNorm;
    if (std::fabs(pirr - pfreeN * ps_prev) < rTol)
        return 1.0;
    return pirr / (pfreeN * ps_prev);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 0: sanity check of the inputs the function depends on.
//
// This is not strictly a test of DDpirr_pfree_ratio_ps() itself, but if the
// lookup helpers return 0 or a non-finite value with our synthetic tables then
// every later assertion would be meaningless, so report that up front.
// -----------------------------------------------------------------------------
void test_ddpirr_lookup_inputs_are_usable()
{
    std::cerr << "\n[TEST] test_ddpirr_lookup_inputs_are_usable\n"
              << "  Source file : src/reactions/DDpirr_pfree_ratio_ps.cpp\n"
              << "  Checking    : the get_prevNorm()/calc_pirr() lookups that\n"
              << "                DDpirr_pfree_ratio_ps() consumes return finite,\n"
              << "                non-zero numbers for our synthetic tables.\n"
              << "  Pass if     : pNorm != 0, and pNorm and pirr are both finite.\n";

    DdpErrorHandlerGuard guard;
    DdpMatrixSet tables = ddp_alloc_tables();

    const double step = ddp_step_size(kDdpDtot, kDdpDeltaT);
    const double r = 1.2;
    const double r0 = 1.1;

    const double pNorm = get_prevNorm(tables.norm, step, r0, kDdpBindRadius);
    const double pirr = calc_pirr(tables.pir, tables.surv, step, r, r0, kDdpBindRadius);

    std::cerr << "  RstepSize = sqrt(Dtot*deltaT)/50 = " << step << '\n'
              << "  get_prevNorm(r0=" << r0 << ") = " << pNorm << '\n'
              << "  calc_pirr(r=" << r << ", r0=" << r0 << ") = " << pirr << '\n';

    EXPECT_TRUE(std::isfinite(pNorm)) << "get_prevNorm returned a non-finite value";
    EXPECT_TRUE(std::isfinite(pirr)) << "calc_pirr returned a non-finite value";
    EXPECT_NE(pNorm, 0.0) << "get_prevNorm returned 0 -> pFree/pNorm would be infinite";

    ddp_free_tables(tables);
}

// -----------------------------------------------------------------------------
// Test 1: the returned ratio equals pirr / ((pFree/pNorm) * ps_prev), with
//         pFree evaluated independently from the closed-form Green's function.
// -----------------------------------------------------------------------------
void test_ddpirr_matches_analytic_greens_function()
{
    std::cerr << "\n[TEST] test_ddpirr_matches_analytic_greens_function\n"
              << "  Source file : src/reactions/DDpirr_pfree_ratio_ps.cpp\n"
              << "  Function    : DDpirr_pfree_ratio_ps (else branch)\n"
              << "  Scenario    : several (r, r0, ps_prev) triples with a tiny rTol\n"
              << "                so the short-circuit never fires.\n"
              << "  Pass if     : returned value == pirr/((pFree/pNorm)*ps_prev)\n"
              << "                where pFree is recomputed independently as\n"
              << "                1/(4 pi D dt) exp(-(r^2+r0^2)/4Ddt) I0(r r0/2Ddt).\n";

    DdpErrorHandlerGuard guard;
    DdpMatrixSet tables = ddp_alloc_tables();

    // (r, r0) pairs kept inside [bindRadius, bindRadius + 0.5] so the table
    // indices stay tiny, plus one case with r == r0 (temp3 == 1 exactly).
    const std::vector<std::pair<double, double>> rPairs {
        { 1.00, 1.10 }, // r sitting exactly on the binding radius
        { 1.05, 1.05 }, // r == r0
        { 1.20, 1.10 }, // generic case
        { 1.50, 1.30 }, // largest separation used here
    };
    const std::vector<double> psPrevList { 1.0, 0.7 };

    const double rTol = 1.0e-15; // effectively disables the |diff| < rTol branch

    for (const auto& pair : rPairs) {
        for (double psPrev : psPrevList) {
            const double r = pair.first;
            const double r0 = pair.second;

            const double actual = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm,
                r, kDdpDtot, kDdpDeltaT, r0, psPrev, rTol, kDdpBindRadius);
            const double expected = ddp_reference_ratio(tables, r, kDdpDtot, kDdpDeltaT, r0,
                psPrev, rTol, kDdpBindRadius);

            std::cerr << "    r = " << r << ", r0 = " << r0 << ", ps_prev = " << psPrev
                      << " -> returned " << actual << ", expected " << expected << '\n';

            // Relative tolerance: the two paths differ only by
            // I0_scaled(x) vs exp(-x)*I0(x), accurate to ~1e-15 relative.
            const double tol = std::fabs(expected) * 1.0e-10 + 1.0e-14;
            EXPECT_NEAR(actual, expected, tol)
                << "ratio mismatch for r=" << r << ", r0=" << r0 << ", ps_prev=" << psPrev;
            EXPECT_TRUE(std::isfinite(actual))
                << "ratio should be finite for r=" << r << ", r0=" << r0;
        }
    }

    ddp_free_tables(tables);
}

// -----------------------------------------------------------------------------
// Test 2: the rTol short-circuit.
//
// The source returns exactly 1.0 when |pirr - pfreeN*ps_prev| < rTol.  We
// compute that difference here and probe just above and just below it.
// -----------------------------------------------------------------------------
void test_ddpirr_returns_unity_within_tolerance()
{
    std::cerr << "\n[TEST] test_ddpirr_returns_unity_within_tolerance\n"
              << "  Source file : src/reactions/DDpirr_pfree_ratio_ps.cpp\n"
              << "  Function    : DDpirr_pfree_ratio_ps (rTol branch)\n"
              << "  Scenario    : same inputs, three different rTol values.\n"
              << "  Pass if     : rTol > |pirr - pfreeN*ps_prev| gives exactly 1.0,\n"
              << "                and rTol just below that gives the true ratio.\n";

    DdpErrorHandlerGuard guard;
    DdpMatrixSet tables = ddp_alloc_tables();

    const double r = 1.2;
    const double r0 = 1.1;
    const double psPrev = 0.8;

    // Rebuild the two quantities the branch compares.
    const double step = ddp_step_size(kDdpDtot, kDdpDeltaT);
    const double pFree = ddp_reference_pfree(r, r0, kDdpDtot, kDdpDeltaT);
    const double pNorm = get_prevNorm(tables.norm, step, r0, kDdpBindRadius);
    const double pirr = calc_pirr(tables.pir, tables.surv, step, r, r0, kDdpBindRadius);
    const double pfreeN = pFree / pNorm;
    const double diff = std::fabs(pirr - pfreeN * psPrev);

    std::cerr << "  pFree = " << pFree << ", pNorm = " << pNorm
              << ", pirr = " << pirr << '\n'
              << "  |pirr - (pFree/pNorm)*ps_prev| = " << diff << '\n';

    // Sanity: the branch test is only meaningful if the difference is non-zero.
    EXPECT_GT(diff, 0.0) << "the two quantities coincide, cannot probe the rTol threshold";

    // (a) A very generous rTol must short-circuit to exactly 1.0.
    const double hugeTol = 1.0e6;
    const double ratioHuge = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm, r,
        kDdpDtot, kDdpDeltaT, r0, psPrev, hugeTol, kDdpBindRadius);
    std::cerr << "  rTol = " << hugeTol << " -> " << ratioHuge << " (expect exactly 1.0)\n";
    EXPECT_DOUBLE_EQ(ratioHuge, 1.0) << "a huge rTol must force the unity short-circuit";

    // (b) rTol just above the difference: still 1.0.
    if (diff > 0.0) {
        const double justAbove = diff * 1.001;
        const double ratioAbove = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm, r,
            kDdpDtot, kDdpDeltaT, r0, psPrev, justAbove, kDdpBindRadius);
        std::cerr << "  rTol = " << justAbove << " (1.001x diff) -> " << ratioAbove
                  << " (expect 1.0)\n";
        EXPECT_DOUBLE_EQ(ratioAbove, 1.0)
            << "rTol slightly larger than the difference should short-circuit";

        // (c) rTol just below the difference: the real ratio comes back.
        const double justBelow = diff * 0.999;
        const double ratioBelow = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm, r,
            kDdpDtot, kDdpDeltaT, r0, psPrev, justBelow, kDdpBindRadius);
        const double expectedBelow = pirr / (pfreeN * psPrev);
        std::cerr << "  rTol = " << justBelow << " (0.999x diff) -> " << ratioBelow
                  << ", expected " << expectedBelow << '\n';
        EXPECT_NEAR(ratioBelow, expectedBelow,
            std::fabs(expectedBelow) * 1.0e-10 + 1.0e-14)
            << "rTol slightly smaller than the difference should give the true ratio";
        EXPECT_NE(ratioBelow, 1.0) << "the computed ratio should not coincidentally be 1.0";
    }

    // (d) rTol == 0: |x| < 0 is never true, so the else branch is always taken.
    const double ratioZeroTol = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm, r,
        kDdpDtot, kDdpDeltaT, r0, psPrev, 0.0, kDdpBindRadius);
    const double expectedZeroTol = pirr / (pfreeN * psPrev);
    std::cerr << "  rTol = 0 -> " << ratioZeroTol << ", expected " << expectedZeroTol << '\n';
    EXPECT_NEAR(ratioZeroTol, expectedZeroTol,
        std::fabs(expectedZeroTol) * 1.0e-10 + 1.0e-14)
        << "rTol == 0 must always take the division branch";

    ddp_free_tables(tables);
}

// -----------------------------------------------------------------------------
// Test 3: the returned ratio scales as 1/ps_prev.
//
// With rTol == 0 the division branch is guaranteed, so
//     ratio(ps_prev) * ps_prev == pirr / pfreeN
// must be independent of ps_prev.  This isolates the role of ps_prev without
// needing to know anything about the lookup tables.
// -----------------------------------------------------------------------------
void test_ddpirr_scales_inversely_with_ps_prev()
{
    std::cerr << "\n[TEST] test_ddpirr_scales_inversely_with_ps_prev\n"
              << "  Source file : src/reactions/DDpirr_pfree_ratio_ps.cpp\n"
              << "  Function    : DDpirr_pfree_ratio_ps (else branch, rTol = 0)\n"
              << "  Scenario    : identical geometry, ps_prev varied.\n"
              << "  Pass if     : ratio * ps_prev is the same constant for every\n"
              << "                ps_prev (i.e. ratio is proportional to 1/ps_prev).\n";

    DdpErrorHandlerGuard guard;
    DdpMatrixSet tables = ddp_alloc_tables();

    const double r = 1.3;
    const double r0 = 1.2;
    const double rTol = 0.0; // forces the division branch for every call

    const std::vector<double> psPrevList { 1.0, 0.5, 0.25, 0.125 };

    double invariant = 0.0;
    bool haveInvariant = false;

    for (double psPrev : psPrevList) {
        const double ratio = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm, r,
            kDdpDtot, kDdpDeltaT, r0, psPrev, rTol, kDdpBindRadius);
        const double product = ratio * psPrev;

        std::cerr << "    ps_prev = " << psPrev << " -> ratio = " << ratio
                  << ", ratio*ps_prev = " << product << '\n';

        EXPECT_TRUE(std::isfinite(ratio)) << "ratio must be finite for ps_prev=" << psPrev;

        if (!haveInvariant) {
            invariant = product;
            haveInvariant = true;
        } else {
            EXPECT_NEAR(product, invariant, std::fabs(invariant) * 1.0e-12 + 1.0e-15)
                << "ratio*ps_prev should be independent of ps_prev (ps_prev=" << psPrev << ')';
        }
    }

    ddp_free_tables(tables);
}

// -----------------------------------------------------------------------------
// Test 4: purity / determinism.
//
// The function must not mutate the lookup tables and must return exactly the
// same double for repeated identical calls.
// -----------------------------------------------------------------------------
void test_ddpirr_is_deterministic_and_leaves_tables_intact()
{
    std::cerr << "\n[TEST] test_ddpirr_is_deterministic_and_leaves_tables_intact\n"
              << "  Source file : src/reactions/DDpirr_pfree_ratio_ps.cpp\n"
              << "  Function    : DDpirr_pfree_ratio_ps\n"
              << "  Scenario    : the same call is issued twice in a row.\n"
              << "  Pass if     : both calls return bit-identical values and the\n"
              << "                three lookup tables are unmodified.\n";

    DdpErrorHandlerGuard guard;
    DdpMatrixSet tables = ddp_alloc_tables();

    // Snapshot a handful of table entries so we can prove they are read-only.
    const double normBefore = gsl_matrix_get(tables.norm, 3, 4);
    const double survBefore = gsl_matrix_get(tables.surv, 5, 6);
    const double pirBefore = gsl_matrix_get(tables.pir, 7, 8);

    const double r = 1.15;
    const double r0 = 1.25;
    const double psPrev = 0.9;
    const double rTol = 1.0e-12;

    const double first = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm, r,
        kDdpDtot, kDdpDeltaT, r0, psPrev, rTol, kDdpBindRadius);
    const double second = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm, r,
        kDdpDtot, kDdpDeltaT, r0, psPrev, rTol, kDdpBindRadius);

    std::cerr << "  first call  = " << first << '\n'
              << "  second call = " << second << '\n';

    EXPECT_DOUBLE_EQ(first, second) << "repeated identical calls must be deterministic";
    EXPECT_TRUE(std::isfinite(first)) << "the returned ratio must be finite";

    EXPECT_DOUBLE_EQ(gsl_matrix_get(tables.norm, 3, 4), normBefore)
        << "normMatrix must not be modified";
    EXPECT_DOUBLE_EQ(gsl_matrix_get(tables.surv, 5, 6), survBefore)
        << "survMatrix must not be modified";
    EXPECT_DOUBLE_EQ(gsl_matrix_get(tables.pir, 7, 8), pirBefore)
        << "pirMatrix must not be modified";

    ddp_free_tables(tables);
}

// -----------------------------------------------------------------------------
// Test 5: sensitivity to the diffusion/time-step inputs.
//
// Changing Dtot*deltaT changes both the Green's function and the RstepSize
// handed to the lookup helpers, so the result must change too (a function that
// silently ignored Dtot or deltaT would fail here).  The exact value is still
// checked against the independent reference.
// -----------------------------------------------------------------------------
void test_ddpirr_responds_to_dtot_and_timestep()
{
    std::cerr << "\n[TEST] test_ddpirr_responds_to_dtot_and_timestep\n"
              << "  Source file : src/reactions/DDpirr_pfree_ratio_ps.cpp\n"
              << "  Function    : DDpirr_pfree_ratio_ps\n"
              << "  Scenario    : Dtot and deltaT are each doubled in turn.\n"
              << "  Pass if     : each result still matches the independent\n"
              << "                reference, and the results actually differ from\n"
              << "                the baseline (arguments are really used).\n";

    DdpErrorHandlerGuard guard;
    DdpMatrixSet tables = ddp_alloc_tables();

    const double r = 1.2;
    const double r0 = 1.1;
    const double psPrev = 1.0;
    const double rTol = 1.0e-15;

    struct Case {
        const char* label;
        double dTot;
        double deltaT;
    };
    const std::vector<Case> cases {
        { "baseline (D=1, dt=1)", 1.0, 1.0 },
        { "D doubled  (D=2, dt=1)", 2.0, 1.0 },
        { "dt doubled (D=1, dt=2)", 1.0, 2.0 },
    };

    double baseline = 0.0;
    bool haveBaseline = false;

    for (const auto& oneCase : cases) {
        const double actual = DDpirr_pfree_ratio_ps(tables.pir, tables.surv, tables.norm, r,
            oneCase.dTot, oneCase.deltaT, r0, psPrev, rTol, kDdpBindRadius);
        const double expected = ddp_reference_ratio(tables, r, oneCase.dTot, oneCase.deltaT, r0,
            psPrev, rTol, kDdpBindRadius);

        std::cerr << "    " << oneCase.label << " -> returned " << actual
                  << ", expected " << expected
                  << " (RstepSize = " << ddp_step_size(oneCase.dTot, oneCase.deltaT) << ")\n";

        EXPECT_NEAR(actual, expected, std::fabs(expected) * 1.0e-10 + 1.0e-14)
            << "ratio mismatch for case: " << oneCase.label;

        if (!haveBaseline) {
            baseline = actual;
            haveBaseline = true;
        } else {
            EXPECT_NE(actual, baseline)
                << "changing " << oneCase.label << " should change the result";
        }
    }

    ddp_free_tables(tables);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* helper runs inside its own TEST so a
// failure in one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(DDpirrPfreeRatioPs, LookupInputsAreUsable) { test_ddpirr_lookup_inputs_are_usable(); }
TEST(DDpirrPfreeRatioPs, MatchesAnalyticGreensFunction) { test_ddpirr_matches_analytic_greens_function(); }
TEST(DDpirrPfreeRatioPs, ReturnsUnityWithinTolerance) { test_ddpirr_returns_unity_within_tolerance(); }
TEST(DDpirrPfreeRatioPs, ScalesInverselyWithPsPrev) { test_ddpirr_scales_inversely_with_ps_prev(); }
TEST(DDpirrPfreeRatioPs, DeterministicAndTablesIntact) { test_ddpirr_is_deterministic_and_leaves_tables_intact(); }
TEST(DDpirrPfreeRatioPs, RespondsToDtotAndTimestep) { test_ddpirr_responds_to_dtot_and_timestep(); }