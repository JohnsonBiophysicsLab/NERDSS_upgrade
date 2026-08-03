/*! \file test_functions_implicitlipid.cpp
 *
 * ### Unit tests for src/reactions/functions_implicitlipid.cpp
 *
 * The file under test provides the closed-form / quadrature based
 * probabilities used by the implicit-lipid (and implicit compartment)
 * reaction model:
 *
 *   - double dissociate2D(paramsIL&)                       2D unbinding prob.
 *   - double function2D(double u, void* parameter)          Bessel integrand
 *   - double integral_for_blockdistance2D(paramsIL&)        integral of the above
 *   - void   block_distance(paramsIL&)                      bisection for R2D
 *   - double pimplicitlipid_2D(paramsIL&)                   2D binding prob.
 *   - double dissociate3D(dt, D, sigma, ka, kb)             3D unbinding prob.
 *   - double pimplicitlipid_3D(double z, paramsIL&)         3D binding prob.
 *   - double prob_entering_compartment(double dr, paramsIL&)
 *   - double prob_exiting_compartment(double dr, paramsIL&)
 *
 * Since these are physical probability expressions, the tests check
 *   (a) the documented "short-circuit" behaviour (rate == 0 -> prob == 0),
 *   (b) exact agreement with an independently written closed form where one
 *       exists (dissociate2D / dissociate3D),
 *   (c) mathematical invariants: finiteness, sign, monotonicity, the
 *       clamping of separations below sigma, and linear scaling with the
 *       binding-site density,
 *   (d) that the numerical bisection in block_distance() lands inside the
 *       physically meaningful window [sigma, Rmax].
 *
 * Every test writes to stderr what is being exercised and why it passes.
 */

#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <limits>

// -----------------------------------------------------------------------------
// Helpers: build well conditioned parameter sets.  paramsIL has *no* default
// member initializers, so every field that the functions read must be set.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Parameter block for the 2D (membrane) routines.
 *
 * Units follow NERDSS conventions: lengths in nm, time in us, D in nm^2/us,
 * 2D ka in nm^2/us, kb in s^-1 (the source divides it by 1e6 internally).
 */
paramsIL ilf_make_2d_params()
{
    paramsIL p;
    p.R2D = 1.0; // will be overwritten by block_distance()
    p.sigma = 1.0; // binding radius, nm
    p.Dtot = 0.1; // 2D diffusion constant, nm^2/us
    p.ka = 10.0; // intrinsic association rate
    p.kb = 1.0; // intrinsic dissociation rate, s^-1
    p.area = 1.0e6; // membrane area, nm^2
    p.dt = 1.0; // time step, us
    p.compartmentR = 0.0; // unused in 2D
    p.compartSiteRho = 0.0; // unused in 2D
    p.Na = 100; // proteins in solution
    p.Nlipid = 1000; // lipids on the surface
    return p;
}

/*! \brief Parameter block for the 3D (solution -> surface) routine. */
paramsIL ilf_make_3d_params()
{
    paramsIL p;
    p.R2D = 0.0;
    p.sigma = 1.0;
    p.Dtot = 1.0; // 3D diffusion constant, nm^2/us
    p.ka = 10.0; // 3D intrinsic ka, nm^3/us
    p.kb = 1.0;
    p.area = 1.0e6;
    p.dt = 1.0;
    p.compartmentR = 0.0;
    p.compartSiteRho = 0.0;
    p.Na = 100;
    p.Nlipid = 1000;
    return p;
}

/*! \brief Parameter block for the spherical compartment (droplet) routines. */
paramsIL ilf_make_compartment_params()
{
    paramsIL p = ilf_make_3d_params();
    p.compartmentR = 100.0; // compartment radius, nm
    p.compartSiteRho = 0.01; // surface density of binding sites, nm^-2
    return p;
}

/*! \brief Independent re-implementation of the dissociate2D closed form.
 *
 * Written directly from the analytic expression so that a typo in the
 * production code would show up as a mismatch.
 */
double ilf_reference_dissociate2D(const paramsIL& p)
{
    const double kb = p.kb / 1.0e6; // s^-1 -> us^-1
    if (kb < 1e-15)
        return 0.0;

    const double KD = kb / p.ka;
    const double maxN = (p.Na > p.Nlipid) ? static_cast<double>(p.Na)
                                          : static_cast<double>(p.Nlipid);
    const double s = p.sigma;
    const double b = 2.0 * std::sqrt(p.area / M_PI / maxN + s * s);
    const double ratio2 = (s / b) * (s / b);
    const double kon = 1.0
        / (1.0 / p.ka
            + 1.0 / (8.0 * M_PI * p.Dtot)
                * (4.0 * std::log(b / s) / ((1.0 - ratio2) * (1.0 - ratio2))
                    - 2.0 / (1.0 - ratio2) - 1.0));
    const double koff = kon * KD;
    return 1.0 - std::exp(-koff * p.dt);
}

/*! \brief Independent re-implementation of the dissociate3D closed form. */
double ilf_reference_dissociate3D(double h, double D, double sigma, double ka, double kbSec)
{
    const double kb = kbSec / 1.0e6;
    if (kb < 1e-15)
        return 0.0;

    const double KD = 2.0 * kb / ka;
    const double kon = 0.5 / (1.0 / ka + 1.0 / (4.0 * M_PI * D * sigma));
    const double koff = kon * KD;
    return 1.0 - std::exp(-koff * h);
}

} // namespace

// -----------------------------------------------------------------------------
// dissociate2D()
// -----------------------------------------------------------------------------
void test_il_dissociate2D()
{
    std::cerr << "\n[TEST] test_il_dissociate2D\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    dissociate2D(paramsIL&)\n"
              << "  Checks:      (1) zero/negligible kb short-circuits to 0,\n"
              << "               (2) exact match with an independently coded\n"
              << "                   closed form,\n"
              << "               (3) probability stays in [0,1),\n"
              << "               (4) probability grows with kb.\n";

    // (1) kb == 0 must short-circuit before any math is done.
    paramsIL zeroKb = ilf_make_2d_params();
    zeroKb.kb = 0.0;
    const double pZero = dissociate2D(zeroKb);
    std::cerr << "  -> kb = 0 gives p = " << pZero << " (expected exactly 0)\n";
    EXPECT_DOUBLE_EQ(pZero, 0.0) << "dissociate2D must return 0 when kb is negligible";

    // (2) Nominal parameters compared against the reference implementation.
    paramsIL p = ilf_make_2d_params();
    const double pCode = dissociate2D(p);
    const double pRef = ilf_reference_dissociate2D(p);
    std::cerr << "  -> kb = " << p.kb << " s^-1 gives p = " << pCode
              << ", reference = " << pRef << "\n";
    EXPECT_TRUE(std::isfinite(pCode)) << "2D unbinding probability must be finite";
    EXPECT_NEAR(pCode, pRef, 1e-14)
        << "dissociate2D should reproduce the analytic koff expression";

    // (3) It is a probability.
    EXPECT_GE(pCode, 0.0) << "probability cannot be negative";
    EXPECT_LT(pCode, 1.0) << "1 - exp(-koff*dt) is strictly below 1 for finite koff";

    // (4) Larger intrinsic off rate -> larger unbinding probability.
    paramsIL fast = ilf_make_2d_params();
    fast.kb = 100.0;
    const double pFast = dissociate2D(fast);
    std::cerr << "  -> kb = 100 s^-1 gives p = " << pFast
              << " (must exceed the kb = 1 s^-1 value)\n";
    EXPECT_GT(pFast, pCode) << "unbinding probability must increase with kb";

    // The input parameter block must not be mutated by dissociate2D.
    EXPECT_DOUBLE_EQ(p.kb, 1.0) << "dissociate2D should not modify kb";
    EXPECT_DOUBLE_EQ(p.R2D, 1.0) << "dissociate2D should not modify R2D";
}

// -----------------------------------------------------------------------------
// dissociate3D()
// -----------------------------------------------------------------------------
void test_il_dissociate3D()
{
    std::cerr << "\n[TEST] test_il_dissociate3D\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    dissociate3D(dt, Dtot, sigma, ka, kb)\n"
              << "  Checks:      zero kb -> 0, agreement with the analytic\n"
              << "               1 - exp(-koff*dt) form, range, monotonicity.\n";

    const double dt = 1.0;
    const double D = 1.0;
    const double sigma = 1.0;
    const double ka = 10.0;

    // kb = 0 must short-circuit.
    const double pZero = dissociate3D(dt, D, sigma, ka, 0.0);
    std::cerr << "  -> kb = 0 gives p = " << pZero << " (expected exactly 0)\n";
    EXPECT_DOUBLE_EQ(pZero, 0.0) << "dissociate3D must return 0 when kb is negligible";

    // Nominal case against the reference formula.
    const double pCode = dissociate3D(dt, D, sigma, ka, 1.0);
    const double pRef = ilf_reference_dissociate3D(dt, D, sigma, ka, 1.0);
    std::cerr << "  -> kb = 1 s^-1 gives p = " << pCode << ", reference = " << pRef << "\n";
    EXPECT_TRUE(std::isfinite(pCode)) << "3D unbinding probability must be finite";
    EXPECT_NEAR(pCode, pRef, 1e-14) << "dissociate3D should match the analytic expression";
    EXPECT_GE(pCode, 0.0) << "probability cannot be negative";
    EXPECT_LT(pCode, 1.0) << "probability must remain below unity";

    // Monotonic in kb and in dt.
    const double pBigKb = dissociate3D(dt, D, sigma, ka, 50.0);
    const double pBigDt = dissociate3D(10.0 * dt, D, sigma, ka, 1.0);
    std::cerr << "  -> larger kb gives " << pBigKb << ", larger dt gives " << pBigDt << "\n";
    EXPECT_GT(pBigKb, pCode) << "unbinding probability must grow with kb";
    EXPECT_GT(pBigDt, pCode) << "unbinding probability must grow with the time step";
}

// -----------------------------------------------------------------------------
// function2D() -- the Bessel-function integrand used by the 2D quadratures.
// -----------------------------------------------------------------------------
void test_il_function2D()
{
    std::cerr << "\n[TEST] test_il_function2D\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    function2D(double u, void* parameter)\n"
              << "  Checks:      the integrand is finite for a range of u and\n"
              << "               decays to (numerically) zero at large u, which\n"
              << "               is the assumption made by the fallback\n"
              << "               quadrature inside the file under test.\n";

    paramsIL p = ilf_make_2d_params();
    p.R2D = p.sigma; // integrand needs a valid lower separation

    // Sample a spread of u values; all must be finite numbers.
    const double uSamples[] = { 1.0e-3, 1.0e-2, 0.1, 1.0, 10.0, 100.0 };
    for (double u : uSamples) {
        const double val = function2D(u, &p);
        std::cerr << "  -> function2D(u = " << u << ") = " << val << "\n";
        EXPECT_TRUE(std::isfinite(val))
            << "integrand must be finite at u = " << u << " (no NaN/Inf allowed)";
    }

    // The production code assumes |function2D(u)| <= 1e-5 by u ~ 1e4 so that
    // the finite-range fallback integration is valid.  Verify that here.
    const double uTail = 1.0e4;
    const double tail = function2D(uTail, &p);
    std::cerr << "  -> tail value function2D(u = " << uTail << ") = " << tail
              << " (must be tiny, |.| < 1e-5)\n";
    EXPECT_TRUE(std::isfinite(tail)) << "tail value must be finite";
    EXPECT_LT(std::fabs(tail), 1.0e-5)
        << "integrand must have decayed by u = 1e4, as assumed by the fallback loop";

    // function2D must not modify the parameters it is handed.
    EXPECT_DOUBLE_EQ(p.R2D, p.sigma) << "function2D must treat its parameters as read-only";
}

// -----------------------------------------------------------------------------
// integral_for_blockdistance2D()
// -----------------------------------------------------------------------------
void test_il_integral_for_blockdistance2D()
{
    std::cerr << "\n[TEST] test_il_integral_for_blockdistance2D\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    integral_for_blockdistance2D(paramsIL&)\n"
              << "  Checks:      the quadrature returns a finite positive value\n"
              << "               and decreases as the lower separation R2D grows\n"
              << "               (this monotonicity is exactly what the bisection\n"
              << "               inside block_distance() relies on).\n";

    paramsIL p = ilf_make_2d_params();
    const double Rmax = p.sigma + 3.0 * std::sqrt(4.0 * p.Dtot * p.dt);
    std::cerr << "  -> sigma = " << p.sigma << ", Rmax = " << Rmax << "\n";

    // Evaluate at the smallest admissible separation.
    p.R2D = p.sigma;
    const double atSigma = integral_for_blockdistance2D(p);
    std::cerr << "  -> integral at R2D = sigma  : " << atSigma << "\n";
    EXPECT_TRUE(std::isfinite(atSigma)) << "integral must be a finite number at R2D = sigma";
    EXPECT_GT(atSigma, 0.0)
        << "the integral is proportional to a binding probability and must be positive";

    // Evaluate midway and at the outer cutoff.
    p.R2D = 0.5 * (p.sigma + Rmax);
    const double atMid = integral_for_blockdistance2D(p);
    std::cerr << "  -> integral at R2D = middle : " << atMid << "\n";
    EXPECT_TRUE(std::isfinite(atMid)) << "integral must be finite at the mid separation";

    p.R2D = Rmax;
    const double atRmax = integral_for_blockdistance2D(p);
    std::cerr << "  -> integral at R2D = Rmax   : " << atRmax << "\n";
    EXPECT_TRUE(std::isfinite(atRmax)) << "integral must be finite at R2D = Rmax";

    // Monotone decreasing in R2D.
    EXPECT_GT(atSigma, atMid) << "integral must decrease as the lower separation increases";
    EXPECT_GT(atMid, atRmax) << "integral must decrease as the lower separation increases";
}

// -----------------------------------------------------------------------------
// block_distance()
// -----------------------------------------------------------------------------
void test_il_block_distance()
{
    std::cerr << "\n[TEST] test_il_block_distance\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    block_distance(paramsIL&)\n"
              << "  Checks:      the bisection writes an R2D inside the physical\n"
              << "               window [sigma, Rmax]; with kb = 0 (no unbinding)\n"
              << "               the block distance collapses onto sigma.\n";

    // Nominal case: R2D must be bracketed by sigma and Rmax.
    paramsIL p = ilf_make_2d_params();
    const double Rmax = p.sigma + 3.0 * std::sqrt(4.0 * p.Dtot * p.dt);
    p.R2D = -1.0; // poison the field so we can see it being written
    block_distance(p);
    std::cerr << "  -> block_distance produced R2D = " << p.R2D
              << " (window [" << p.sigma << ", " << Rmax << "])\n";
    EXPECT_TRUE(std::isfinite(p.R2D)) << "block distance must be a finite number";
    EXPECT_GE(p.R2D, p.sigma - 1e-6) << "block distance cannot be smaller than sigma";
    EXPECT_LE(p.R2D, Rmax + 1e-6) << "block distance cannot exceed Rmax";

    // kb = 0 -> left-hand side of the bisection is 0, so the search always
    // moves the upper bound down and R2D converges to sigma.
    paramsIL noOff = ilf_make_2d_params();
    noOff.kb = 0.0;
    noOff.R2D = -1.0;
    block_distance(noOff);
    std::cerr << "  -> with kb = 0, R2D = " << noOff.R2D
              << " (expected to converge to sigma = " << noOff.sigma << ")\n";
    EXPECT_NEAR(noOff.R2D, noOff.sigma, 1.0e-4)
        << "with no unbinding the block distance must collapse to sigma";
}

// -----------------------------------------------------------------------------
// pimplicitlipid_2D()
// -----------------------------------------------------------------------------
void test_il_pimplicitlipid_2D()
{
    std::cerr << "\n[TEST] test_il_pimplicitlipid_2D\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    pimplicitlipid_2D(paramsIL&)\n"
              << "  Checks:      ka = 0 short-circuits to 0 without touching the\n"
              << "               quadrature; a physical ka gives a finite,\n"
              << "               positive value and leaves R2D inside\n"
              << "               [sigma, Rmax] (block_distance is called first).\n";

    // ka == 0: cheap early return, and R2D must be untouched.
    paramsIL noOn = ilf_make_2d_params();
    noOn.ka = 0.0;
    noOn.R2D = 42.0; // sentinel
    const double pZero = pimplicitlipid_2D(noOn);
    std::cerr << "  -> ka = 0 gives p = " << pZero << " (expected exactly 0)\n";
    EXPECT_DOUBLE_EQ(pZero, 0.0) << "pimplicitlipid_2D must return 0 for a zero ka";
    EXPECT_DOUBLE_EQ(noOn.R2D, 42.0)
        << "the ka == 0 early return must happen before block_distance()";

    // Nominal case.
    paramsIL p = ilf_make_2d_params();
    const double Rmax = p.sigma + 3.0 * std::sqrt(4.0 * p.Dtot * p.dt);
    const double prob = pimplicitlipid_2D(p);
    std::cerr << "  -> ka = " << p.ka << " gives p (per unit lipid density) = " << prob
              << ", R2D set to " << p.R2D << "\n";
    EXPECT_TRUE(std::isfinite(prob)) << "2D binding probability must be finite";
    EXPECT_GT(prob, 0.0) << "2D binding probability must be positive for a positive ka";
    EXPECT_GE(p.R2D, p.sigma - 1e-6)
        << "pimplicitlipid_2D must leave a valid block distance (>= sigma)";
    EXPECT_LE(p.R2D, Rmax + 1e-6)
        << "pimplicitlipid_2D must leave a valid block distance (<= Rmax)";
}

// -----------------------------------------------------------------------------
// pimplicitlipid_3D()
// -----------------------------------------------------------------------------
void test_il_pimplicitlipid_3D()
{
    std::cerr << "\n[TEST] test_il_pimplicitlipid_3D\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    pimplicitlipid_3D(double z, paramsIL&)\n"
              << "  Checks:      ka = 0 -> 0; positive and finite for z >= sigma;\n"
              << "               monotone decreasing with the height z; the two\n"
              << "               analytic branches agree at z = sigma; heights\n"
              << "               below sigma are clamped to sigma.\n";

    // ka == 0 short-circuit.
    paramsIL noOn = ilf_make_3d_params();
    noOn.ka = 0.0;
    const double pZero = pimplicitlipid_3D(2.0, noOn);
    std::cerr << "  -> ka = 0 gives p = " << pZero << " (expected exactly 0)\n";
    EXPECT_DOUBLE_EQ(pZero, 0.0) << "pimplicitlipid_3D must return 0 for a zero ka";

    paramsIL p = ilf_make_3d_params();

    // Positive and finite just above contact, and monotone decreasing with z.
    const double pNear = pimplicitlipid_3D(1.5 * p.sigma, p);
    const double pMid = pimplicitlipid_3D(3.0 * p.sigma, p);
    const double pFar = pimplicitlipid_3D(20.0 * p.sigma, p);
    std::cerr << "  -> p(z = 1.5 sigma) = " << pNear << "\n"
              << "  -> p(z = 3.0 sigma) = " << pMid << "\n"
              << "  -> p(z = 20  sigma) = " << pFar << "\n";
    EXPECT_TRUE(std::isfinite(pNear)) << "probability must be finite near contact";
    EXPECT_TRUE(std::isfinite(pMid)) << "probability must be finite at intermediate z";
    EXPECT_TRUE(std::isfinite(pFar)) << "probability must be finite far from the surface";
    EXPECT_GT(pNear, 0.0) << "binding probability must be positive near the surface";
    EXPECT_GT(pNear, pMid) << "binding probability must fall off with height";
    EXPECT_GT(pMid, pFar) << "binding probability must fall off with height";
    EXPECT_LT(pFar, pNear * 1.0e-2)
        << "far from the surface the binding probability must be strongly suppressed";

    // The z > sigma branch and the z <= sigma branch must agree at contact:
    // both reduce to the same prefactor times the same bracket.
    const double pAtSigma = pimplicitlipid_3D(p.sigma, p);
    const double pJustAbove = pimplicitlipid_3D(p.sigma * (1.0 + 1.0e-9), p);
    std::cerr << "  -> p(z = sigma)        = " << pAtSigma << "\n"
              << "  -> p(z = sigma + eps)  = " << pJustAbove
              << " (branches must agree at contact)\n";
    EXPECT_NEAR(pAtSigma, pJustAbove, 1.0e-8 * std::fabs(pAtSigma) + 1.0e-12)
        << "the two analytic branches of pimplicitlipid_3D must join at z = sigma";

    // Heights below sigma are clamped to sigma internally.
    const double pBelow = pimplicitlipid_3D(0.25 * p.sigma, p);
    std::cerr << "  -> p(z = 0.25 sigma)   = " << pBelow
              << " (must equal the z = sigma value because z is clamped)\n";
    EXPECT_NEAR(pBelow, pAtSigma, 1.0e-12)
        << "separations below sigma must be clamped to sigma";
}

// -----------------------------------------------------------------------------
// prob_entering_compartment()
// -----------------------------------------------------------------------------
void test_il_prob_entering_compartment()
{
    std::cerr << "\n[TEST] test_il_prob_entering_compartment\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    prob_entering_compartment(double dr, paramsIL&)\n"
              << "  Checks:      ka = 0 -> 0; positive, finite and below unity;\n"
              << "               decreasing with the distance dr outside the\n"
              << "               compartment; clamped for dr < sigma; exactly\n"
              << "               linear in the surface site density.\n";

    // ka == 0 short-circuit.
    paramsIL noOn = ilf_make_compartment_params();
    noOn.ka = 0.0;
    const double pZero = prob_entering_compartment(2.0, noOn);
    std::cerr << "  -> ka = 0 gives p = " << pZero << " (expected exactly 0)\n";
    EXPECT_DOUBLE_EQ(pZero, 0.0) << "prob_entering_compartment must return 0 for a zero ka";

    paramsIL p = ilf_make_compartment_params();

    // Finite, positive, sub-unity probabilities that decay with distance.
    const double pNear = prob_entering_compartment(p.sigma, p);
    const double pMid = prob_entering_compartment(2.0 * p.sigma, p);
    const double pFar = prob_entering_compartment(5.0 * p.sigma, p);
    std::cerr << "  -> p(dr = sigma)   = " << pNear << "\n"
              << "  -> p(dr = 2 sigma) = " << pMid << "\n"
              << "  -> p(dr = 5 sigma) = " << pFar << "\n";
    EXPECT_TRUE(std::isfinite(pNear)) << "probability must be finite at contact";
    EXPECT_TRUE(std::isfinite(pMid)) << "probability must be finite at 2 sigma";
    EXPECT_TRUE(std::isfinite(pFar)) << "probability must be finite at 5 sigma";
    EXPECT_GT(pNear, 0.0) << "entering probability must be positive at contact";
    EXPECT_LE(pNear, 1.0) << "entering probability must not exceed unity";
    EXPECT_GT(pNear, pMid) << "entering probability must decay with distance";
    EXPECT_GT(pMid, pFar) << "entering probability must decay with distance";

    // Distances below sigma are clamped to sigma inside the function.
    const double pBelow = prob_entering_compartment(0.25 * p.sigma, p);
    std::cerr << "  -> p(dr = 0.25 sigma) = " << pBelow
              << " (must equal the dr = sigma value because dr is clamped)\n";
    EXPECT_NEAR(pBelow, pNear, 1.0e-12)
        << "distances below sigma must be clamped to sigma";

    // The prefactor is strictly linear in compartSiteRho.
    paramsIL doubled = ilf_make_compartment_params();
    doubled.compartSiteRho = 2.0 * p.compartSiteRho;
    const double pDoubled = prob_entering_compartment(p.sigma, doubled);
    std::cerr << "  -> doubling rho gives " << pDoubled << " vs 2 * " << pNear
              << " = " << 2.0 * pNear << "\n";
    EXPECT_NEAR(pDoubled, 2.0 * pNear, 1.0e-10 * std::fabs(pNear) + 1.0e-15)
        << "entering probability must scale linearly with the site density";
}

// -----------------------------------------------------------------------------
// prob_exiting_compartment()
// -----------------------------------------------------------------------------
void test_il_prob_exiting_compartment()
{
    std::cerr << "\n[TEST] test_il_prob_exiting_compartment\n"
              << "  Source file: functions_implicitlipid.cpp\n"
              << "  Function:    prob_exiting_compartment(double dr, paramsIL&)\n"
              << "  Checks:      ka = 0 -> 0; positive, finite and below unity;\n"
              << "               decreasing with the depth dr inside the\n"
              << "               compartment; clamped for dr < sigma; exactly\n"
              << "               linear in the surface site density.\n";

    // ka == 0 short-circuit.
    paramsIL noOn = ilf_make_compartment_params();
    noOn.ka = 0.0;
    const double pZero = prob_exiting_compartment(2.0, noOn);
    std::cerr << "  -> ka = 0 gives p = " << pZero << " (expected exactly 0)\n";
    EXPECT_DOUBLE_EQ(pZero, 0.0) << "prob_exiting_compartment must return 0 for a zero ka";

    paramsIL p = ilf_make_compartment_params();

    // Finite, positive, sub-unity probabilities that decay with depth.
    const double pNear = prob_exiting_compartment(p.sigma, p);
    const double pMid = prob_exiting_compartment(2.0 * p.sigma, p);
    const double pFar = prob_exiting_compartment(5.0 * p.sigma, p);
    std::cerr << "  -> p(dr = sigma)   = " << pNear << "\n"
              << "  -> p(dr = 2 sigma) = " << pMid << "\n"
              << "  -> p(dr = 5 sigma) = " << pFar << "\n";
    EXPECT_TRUE(std::isfinite(pNear)) << "probability must be finite at contact";
    EXPECT_TRUE(std::isfinite(pMid)) << "probability must be finite at 2 sigma";
    EXPECT_TRUE(std::isfinite(pFar)) << "probability must be finite at 5 sigma";
    EXPECT_GT(pNear, 0.0) << "exiting probability must be positive at contact";
    EXPECT_LE(pNear, 1.0) << "exiting probability must not exceed unity";
    EXPECT_GT(pNear, pMid) << "exiting probability must decay with depth";
    EXPECT_GT(pMid, pFar) << "exiting probability must decay with depth";

    // Depths below sigma are clamped to sigma inside the function.
    const double pBelow = prob_exiting_compartment(0.25 * p.sigma, p);
    std::cerr << "  -> p(dr = 0.25 sigma) = " << pBelow
              << " (must equal the dr = sigma value because dr is clamped)\n";
    EXPECT_NEAR(pBelow, pNear, 1.0e-12)
        << "depths below sigma must be clamped to sigma";

    // The prefactor is strictly linear in compartSiteRho.
    paramsIL doubled = ilf_make_compartment_params();
    doubled.compartSiteRho = 2.0 * p.compartSiteRho;
    const double pDoubled = prob_exiting_compartment(p.sigma, doubled);
    std::cerr << "  -> doubling rho gives " << pDoubled << " vs 2 * " << pNear
              << " = " << 2.0 * pNear << "\n";
    EXPECT_NEAR(pDoubled, 2.0 * pNear, 1.0e-10 * std::fabs(pNear) + 1.0e-15)
        << "exiting probability must scale linearly with the site density";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each helper is invoked from its own TEST so that the
// framework reports the results individually while still running all of them
// (only non-fatal EXPECT_* assertions are used above).
// -----------------------------------------------------------------------------
TEST(FunctionsImplicitLipid, Dissociate2D) { test_il_dissociate2D(); }
TEST(FunctionsImplicitLipid, Dissociate3D) { test_il_dissociate3D(); }
TEST(FunctionsImplicitLipid, Function2D) { test_il_function2D(); }
TEST(FunctionsImplicitLipid, IntegralForBlockDistance2D) { test_il_integral_for_blockdistance2D(); }
TEST(FunctionsImplicitLipid, BlockDistance) { test_il_block_distance(); }
TEST(FunctionsImplicitLipid, PImplicitLipid2D) { test_il_pimplicitlipid_2D(); }
TEST(FunctionsImplicitLipid, PImplicitLipid3D) { test_il_pimplicitlipid_3D(); }
TEST(FunctionsImplicitLipid, ProbEnteringCompartment) { test_il_prob_entering_compartment(); }
TEST(FunctionsImplicitLipid, ProbExitingCompartment) { test_il_prob_exiting_compartment(); }