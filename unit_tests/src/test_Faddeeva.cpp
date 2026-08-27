/*! \file test_Faddeeva.cpp
 *
 * ### Unit tests for ../src/math/Faddeeva.cpp
 *
 * The file under test provides the `Faddeeva::` namespace, an implementation of
 * the scaled complex error function ("Faddeeva function")
 *
 *      w(z) = exp(-z^2) * erfc(-i z)
 *
 * together with a family of related special functions built on top of it:
 *
 *   - std::complex<double> Faddeeva::w(std::complex<double> z, double relerr)
 *   - double               Faddeeva::w_im(double x)          (= 2*Dawson(x)/sqrt(pi))
 *   - std::complex<double> Faddeeva::erfcx(std::complex<double> z, double relerr)
 *   - double               Faddeeva::erfcx(double x)
 *   - std::complex<double> Faddeeva::erf(std::complex<double> z, double relerr)
 *   - double               Faddeeva::erf(double x)
 *   - std::complex<double> Faddeeva::erfi(std::complex<double> z, double relerr)
 *   - double               Faddeeva::erfi(double x)
 *   - std::complex<double> Faddeeva::erfc(std::complex<double> z, double relerr)
 *   - double               Faddeeva::erfc(double x)
 *   - std::complex<double> Faddeeva::Dawson(std::complex<double> z, double relerr)
 *   - double               Faddeeva::Dawson(double x)
 *
 * The reference values used below are the high-precision (Maple / WolframAlpha)
 * values that ship inside the `#ifdef TEST_FADDEEVA` block of the source file
 * itself, plus independent cross-checks against the C++11 `std::erf`/`std::erfc`
 * from <cmath> and against the analytic identities that relate the functions.
 *
 * NOTE: NaN inputs are deliberately NOT exercised.  The lookup-table helpers do
 * `switch ((int) y100)` and converting a NaN to `int` is undefined behaviour, so
 * a NaN test would be testing the compiler rather than this code.  Infinities on
 * the real axis are well defined and are tested.
 */

#include "math/Faddeeva.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

// Convenience alias for the complex type used throughout the Faddeeva API.
using FadComplex = std::complex<double>;

// sqrt(pi)/2 and 2/sqrt(pi), used by several analytic identities below.
const double kFadSpi2 = 0.8862269254527580136490837416705725913990; // sqrt(pi)/2
const double kFad2ospi = 1.1283791670955125738961589031215451716769; // 2/sqrt(pi)
const double kFadIspi = 0.56418958354775628694807945156; // 1/sqrt(pi)

/*! \brief Compute the relative error of `actual` with respect to `expected`.
 *
 * If `expected` is exactly zero the absolute error is returned instead, so the
 * printed diagnostics remain meaningful.
 */
double fad_relerr(double expected, double actual)
{
    if (expected == 0.0)
        return std::fabs(actual);
    return std::fabs((actual - expected) / expected);
}

/*! \brief Non-fatal check of one real value against a reference value.
 *
 * The pass criterion is  |actual - expected| <= max(absFloor, relTol*|expected|).
 * The value, the reference and the achieved relative error are echoed to stderr
 * so a reader can see exactly what was compared.
 */
void fad_expect_scalar(const std::string& label, double expected, double actual,
    double relTol, double absFloor = 0.0)
{
    const double tol = std::max(absFloor, relTol * std::fabs(expected));
    EXPECT_NEAR(actual, expected, tol) << label;
    std::cerr << "      " << label << ": got " << std::setprecision(17) << actual
              << ", expected " << expected << "  (rel.err = " << std::setprecision(3)
              << fad_relerr(expected, actual) << ", tol = " << relTol << ")\n";
}

/*! \brief Non-fatal check of one complex value against a reference value.
 *
 * Each component is compared with a tolerance that is scaled both by that
 * component and by the magnitude of the whole complex number.  Scaling by the
 * magnitude is what allows components whose reference value is exactly zero
 * (e.g. Im w(iy) == 0) to be checked sensibly rather than requiring bit
 * equality with 0.
 */
void fad_expect_complex(const std::string& label, FadComplex expected,
    FadComplex actual, double relTol)
{
    const double scale = std::abs(expected);
    const double tolRe = relTol * (std::fabs(expected.real()) + scale);
    const double tolIm = relTol * (std::fabs(expected.imag()) + scale);

    EXPECT_NEAR(actual.real(), expected.real(), tolRe) << label << " [real part]";
    EXPECT_NEAR(actual.imag(), expected.imag(), tolIm) << label << " [imag part]";

    std::cerr << "      " << label << ": got (" << std::setprecision(17) << actual.real()
              << ", " << actual.imag() << "i)\n"
              << "        expected  (" << expected.real() << ", " << expected.imag()
              << "i)  rel.err = " << std::setprecision(3)
              << fad_relerr(expected.real(), actual.real()) << " / "
              << fad_relerr(expected.imag(), actual.imag()) << "\n";
}

} // namespace

// -----------------------------------------------------------------------------
// 1. Faddeeva::w(complex) against high-precision reference values.
//
//    The sample points are chosen to walk through every major branch of the
//    implementation: the Algorithm-916 sum (|z| small/moderate), the
//    continued-fraction expansion (|z| large), the pure real axis, the pure
//    imaginary axis, and negative imaginary parts (which use the reflection
//    w(z) = 2 exp(-z^2) - w(-z)).
// -----------------------------------------------------------------------------
void faddeeva_test_w_known_values()
{
    std::cerr << "\n[TEST] faddeeva_test_w_known_values\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::w(std::complex<double>, double)\n"
              << "  Checking w(z) at reference points spanning the Algorithm-916,\n"
              << "  continued-fraction, real-axis and imaginary-axis branches.\n"
              << "  Pass criterion: componentwise relative error < 1e-13.\n";

    struct Sample {
        const char* name;
        FadComplex z;
        FadComplex w;
    };

    // Reference values reproduced from the (Maple/WolframAlpha generated) table
    // embedded in Faddeeva.cpp itself.
    const std::vector<Sample> samples = {
        { "w(-1+1i)   [Alg.916 branch]",
            FadComplex(-1.0, 1.0),
            FadComplex(0.30474420525691259245713884106959496013413834051768,
                -0.20821893820283162728743734725471561394145872072738) },
        { "w(0.6+2i)  [Alg.916 branch]",
            FadComplex(0.6, 2.0),
            FadComplex(0.2410250715772692146133539023007113781272362309451,
                0.06087579663428089745895459735240964093522265589350) },
        { "w(-0.4+3i) [Alg.916 branch]",
            FadComplex(-0.4, 3.0),
            FadComplex(0.1764906227004816847297495349730234591778719532788,
                -0.02146550539468457616788719893991501311573031095617) },
        { "w(-0.7-0.7i) [negative Im, reflection]",
            FadComplex(-0.7, -0.7),
            FadComplex(0.69504753678406939989115375989939096800793577783885,
                -1.8916411171103639136680830887017670616339912024317) },
        { "w(-1-9i)   [large |exp(-z^2)| reflection]",
            FadComplex(-1.0, -9.0),
            FadComplex(7.317131068972378096865595229600561710140617977e34,
                8.321873499714402777186848353320412813066170427e34) },
        { "w(11+1i)   [continued fraction]",
            FadComplex(11.0, 1.0),
            FadComplex(0.00468190164965444174367477874864366058339647648741,
                0.0510735563901306197993676329845149741675029197050) },
        { "w(-22-2i)  [continued fraction, negative Im]",
            FadComplex(-22.0, -2.0),
            FadComplex(-0.0023193175200187620902125853834909543869428763219,
                -0.025460054739731556004902057663500272721780776336) },
        { "w(1e5+1e5i) [asymptotic i/(sqrt(pi) z)]",
            FadComplex(1e5, 1e5),
            FadComplex(2.820947917809305132678577516325951485807107151e-6,
                2.820947917668257736791638444590253942253354058e-6) },
        { "w(1+0i)    [real axis]",
            FadComplex(1.0, 0.0),
            FadComplex(0.36787944117144232159552377016146086744581113103176,
                0.60715770584139372911503823580074492116122092866515) },
        { "w(0+0.12345i) [imaginary axis]",
            FadComplex(0.0, 0.12345),
            FadComplex(0.8746342859608052666092782112565360755791467973338452, 0.0) },
    };

    for (const auto& s : samples) {
        const FadComplex got = Faddeeva::w(s.z, 0.0);
        fad_expect_complex(s.name, s.w, got, 1e-13);
    }
}

// -----------------------------------------------------------------------------
// 2. Structural properties of w(z): the special-cased axes and the mirror
//    symmetry w(-conj(z)) = conj(w(z)).
// -----------------------------------------------------------------------------
void faddeeva_test_w_axes_and_symmetry()
{
    std::cerr << "\n[TEST] faddeeva_test_w_axes_and_symmetry\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::w (axis special cases + mirror symmetry)\n"
              << "  1) w(x + 0i) must be exactly (exp(-x^2), w_im(x)).\n"
              << "  2) w(0 + iy) must be exactly (erfcx(y), 0).\n"
              << "  3) w(-conj z) must equal conj(w(z)) to ~1e-13 relative.\n";

    // 1) On the real axis the implementation short-circuits to
    //    (exp(-x^2), w_im(x)) -- these must agree bit-for-bit (4 ULP).
    for (double x : { 0.25, 1.0, 3.5, 12.0, -2.0 }) {
        const FadComplex got = Faddeeva::w(FadComplex(x, 0.0));
        std::cerr << "      real axis x = " << x << " -> w = (" << std::setprecision(12)
                  << got.real() << ", " << got.imag() << "i)\n";
        EXPECT_DOUBLE_EQ(got.real(), std::exp(-x * x))
            << "Re w(x) must equal exp(-x*x) for real x = " << x;
        EXPECT_DOUBLE_EQ(got.imag(), Faddeeva::w_im(x))
            << "Im w(x) must equal w_im(x) for real x = " << x;
    }

    // 2) On the imaginary axis w(iy) = erfcx(y) with a zero imaginary part.
    for (double y : { 0.12345, 1.0, 8.0, -0.5 }) {
        const FadComplex got = Faddeeva::w(FadComplex(0.0, y));
        std::cerr << "      imag axis y = " << y << " -> w = (" << std::setprecision(12)
                  << got.real() << ", " << got.imag() << "i)\n";
        EXPECT_DOUBLE_EQ(got.real(), Faddeeva::erfcx(y))
            << "Re w(iy) must equal erfcx(y) for y = " << y;
        EXPECT_DOUBLE_EQ(got.imag(), 0.0) << "Im w(iy) must vanish for y = " << y;
    }

    // 3) Mirror symmetry of the Faddeeva function.
    const std::vector<FadComplex> pts = { FadComplex(1.0, 2.0), FadComplex(0.3, 0.7),
        FadComplex(6.5, 0.05), FadComplex(2.0, -1.5) };
    for (const auto& z : pts) {
        const FadComplex lhs = Faddeeva::w(FadComplex(-z.real(), z.imag())); // w(-conj z)
        const FadComplex rhs = std::conj(Faddeeva::w(z));
        std::cerr << "      symmetry check at z = (" << z.real() << ", " << z.imag()
                  << "i)\n";
        fad_expect_complex("        w(-conj z) vs conj(w(z))", rhs, lhs, 1e-13);
    }
}

// -----------------------------------------------------------------------------
// 3. The `relerr` argument of w(): requesting a coarser tolerance must still
//    return a result that is within (roughly) that tolerance of the fully
//    converged answer.  This exercises the "relerr != DBL_EPSILON" code path
//    which recomputes exp(-a2*n*n) on the fly instead of using the lookup table.
// -----------------------------------------------------------------------------
void faddeeva_test_w_relerr_parameter()
{
    std::cerr << "\n[TEST] faddeeva_test_w_relerr_parameter\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::w(z, relerr) with a coarse relerr\n"
              << "  Pass criterion: |w(z,1e-6) - w(z,0)| / |w(z,0)| < 1e-4\n"
              << "  (the argument is a requested accuracy, not a guarantee, so the\n"
              << "   comparison is deliberately looser than the request).\n";

    const std::vector<FadComplex> pts
        = { FadComplex(-0.7, -0.7), FadComplex(1.0, 2.0), FadComplex(0.05, 0.05) };

    for (const auto& z : pts) {
        const FadComplex exact = Faddeeva::w(z, 0.0);
        const FadComplex coarse = Faddeeva::w(z, 1e-6);
        const double err = std::abs(coarse - exact) / std::abs(exact);
        std::cerr << "      z = (" << z.real() << ", " << z.imag()
                  << "i): |dw|/|w| = " << std::setprecision(4) << err << '\n';
        EXPECT_LT(err, 1e-4)
            << "Coarse-tolerance w(z, 1e-6) deviates too much from w(z, 0)";
    }

    // Requesting an absurdly loose tolerance is clamped internally to 0.1; the
    // result must still be finite and of the right order of magnitude.
    const FadComplex clamped = Faddeeva::w(FadComplex(0.5, 0.5), 10.0);
    const FadComplex exact = Faddeeva::w(FadComplex(0.5, 0.5), 0.0);
    std::cerr << "      relerr=10 (clamped to 0.1) gives (" << clamped.real() << ", "
              << clamped.imag() << "i)\n";
    EXPECT_TRUE(std::isfinite(clamped.real()) && std::isfinite(clamped.imag()))
        << "w() with an over-large relerr must still return a finite value";
    EXPECT_LT(std::abs(clamped - exact) / std::abs(exact), 0.5)
        << "Even a clamped relerr should stay in the right ballpark";
}

// -----------------------------------------------------------------------------
// 4. Faddeeva::w_im(double) -- the scaled Dawson integral for real arguments.
// -----------------------------------------------------------------------------
void faddeeva_test_w_im()
{
    std::cerr << "\n[TEST] faddeeva_test_w_im\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::w_im(double)\n"
              << "  Checks: w_im(0)=0, reference value at x=1, oddness, the small-x\n"
              << "  Taylor branch, and the large-x continued-fraction asymptotics.\n";

    // w_im(0) is exactly zero (the Taylor branch evaluates x*(...)).
    std::cerr << "    -> w_im(0) should be exactly 0\n";
    EXPECT_DOUBLE_EQ(Faddeeva::w_im(0.0), 0.0) << "w_im(0) must be 0";

    // Reference value: Im w(1) from the table in Faddeeva.cpp.
    fad_expect_scalar("w_im(1)", 0.60715770584139372911503823580074492116122092866515,
        Faddeeva::w_im(1.0), 1e-13);

    // Oddness: w_im(-x) == -w_im(x) (implemented explicitly in the source).
    std::cerr << "    -> checking oddness w_im(-x) == -w_im(x)\n";
    for (double x : { 0.001, 0.5, 1.0, 12.0, 60.0 }) {
        EXPECT_DOUBLE_EQ(Faddeeva::w_im(-x), -Faddeeva::w_im(x))
            << "w_im must be an odd function; failed at x = " << x;
    }

    // Small-x Taylor branch: w_im(x) -> 2/sqrt(pi) * x as x -> 0.
    std::cerr << "    -> small-x Taylor branch (x = 1e-8) should give 2x/sqrt(pi)\n";
    fad_expect_scalar("w_im(1e-8)", kFad2ospi * 1e-8, Faddeeva::w_im(1e-8), 1e-14);

    // Large-x continued fraction vs the asymptotic series
    //     w_im(x) ~ 1/(x sqrt(pi)) [1 + 1/(2x^2) + 3/(4x^4) + 15/(8x^6)]
    std::cerr << "    -> large-x branch (x = 50) vs asymptotic series\n";
    {
        const double x = 50.0;
        const double x2 = x * x;
        const double asym = kFadIspi / x
            * (1.0 + 0.5 / x2 + 0.75 / (x2 * x2) + 1.875 / (x2 * x2 * x2));
        fad_expect_scalar("w_im(50)", asym, Faddeeva::w_im(x), 1e-10);
    }

    // Extremely large x uses the one-term expansion 1/(x sqrt(pi)) exactly.
    std::cerr << "    -> huge-x branch (x = 1e8) uses the exact 1-term expansion\n";
    EXPECT_DOUBLE_EQ(Faddeeva::w_im(1e8), kFadIspi / 1e8)
        << "w_im(1e8) must reduce to 1/(x*sqrt(pi))";
}

// -----------------------------------------------------------------------------
// 5. Faddeeva::erfcx(double) -- the scaled complementary error function.
// -----------------------------------------------------------------------------
void faddeeva_test_erfcx_real()
{
    std::cerr << "\n[TEST] faddeeva_test_erfcx_real\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::erfcx(double)\n"
              << "  Checks: erfcx(0)=1, erfcx(x)=exp(x^2)erfc(x) against std::erfc,\n"
              << "  the negative-x mirror formula, and the large-|x| expansions.\n";

    // erfcx(0) falls out of the Chebyshev table as exactly 1.
    std::cerr << "    -> erfcx(0) must be exactly 1\n";
    EXPECT_DOUBLE_EQ(Faddeeva::erfcx(0.0), 1.0) << "erfcx(0) must be 1";

    // Reference values.
    fad_expect_scalar("erfcx(1)", 0.42758357615580700442, Faddeeva::erfcx(1.0), 1e-13);
    fad_expect_scalar("erfcx(-1)", 2.0 * std::exp(1.0) - 0.42758357615580700442,
        Faddeeva::erfcx(-1.0), 1e-13);

    // Independent cross-check against <cmath>: erfcx(x) == exp(x^2)*erfc(x) for
    // arguments where exp(x^2) does not overflow.
    std::cerr << "    -> cross-checking erfcx(x) == exp(x*x)*std::erfc(x)\n";
    for (double x : { -3.0, -1.5, -0.25, 0.0, 0.25, 1.5, 3.0, 8.0, 20.0 }) {
        const double ref = std::exp(x * x) * std::erfc(x);
        fad_expect_scalar("erfcx(" + std::to_string(x) + ")", ref, Faddeeva::erfcx(x),
            1e-12);
    }

    // 5-term continued fraction region (x > 50), checked against the asymptotic
    // series erfcx(x) ~ 1/(x sqrt(pi)) [1 - 1/(2x^2) + 3/(4x^4) - 15/(8x^6)].
    std::cerr << "    -> continued-fraction region (x = 60) vs asymptotic series\n";
    {
        const double x = 60.0;
        const double x2 = x * x;
        const double asym = kFadIspi / x
            * (1.0 - 0.5 / x2 + 0.75 / (x2 * x2) - 1.875 / (x2 * x2 * x2));
        fad_expect_scalar("erfcx(60)", asym, Faddeeva::erfcx(x), 1e-11);
    }

    // Very large x uses the single-term expansion exactly.
    std::cerr << "    -> huge-x branch (x = 1e8) uses the exact 1-term expansion\n";
    EXPECT_DOUBLE_EQ(Faddeeva::erfcx(1e8), kFadIspi / 1e8)
        << "erfcx(1e8) must reduce to 1/(x*sqrt(pi))";

    // Deeply negative x: the implementation returns the leading 2*exp(x^2) term,
    // and saturates to +inf below x = -26.7.
    std::cerr << "    -> deeply negative x uses 2*exp(x*x), and overflows below -26.7\n";
    EXPECT_DOUBLE_EQ(Faddeeva::erfcx(-10.0), 2.0 * std::exp(100.0))
        << "erfcx(-10) must be the leading 2*exp(x*x) term";
    EXPECT_TRUE(std::isinf(Faddeeva::erfcx(-27.0)))
        << "erfcx(-27) must overflow to infinity";
}

// -----------------------------------------------------------------------------
// 6. Faddeeva::erfcx(complex) -- must be exactly w(i z), and must (by design)
//    ignore the relerr argument entirely.
// -----------------------------------------------------------------------------
void faddeeva_test_erfcx_complex()
{
    std::cerr << "\n[TEST] faddeeva_test_erfcx_complex\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::erfcx(std::complex<double>, double)\n"
              << "  Checks: reference value, the identity erfcx(z) == w(i z), and\n"
              << "  the documented fact that the relerr argument is ignored.\n";

    // Reference value from the table inside Faddeeva.cpp.
    fad_expect_complex("erfcx(1.234+0.5678i)",
        FadComplex(0.3382187479799972294747793561190487832579,
            -0.1116077470811648467464927471872945833154),
        Faddeeva::erfcx(FadComplex(1.234, 0.5678), 0.0), 1e-13);

    // erfcx(z) is implemented as w(i z) -- results must be bit identical.
    std::cerr << "    -> erfcx(z) must equal w(i*z) exactly\n";
    const std::vector<FadComplex> pts
        = { FadComplex(1.0, 2.0), FadComplex(-0.5, 0.25), FadComplex(0.0, 3.0) };
    for (const auto& z : pts) {
        const FadComplex viaW = Faddeeva::w(FadComplex(-z.imag(), z.real()));
        const FadComplex direct = Faddeeva::erfcx(z, 0.0);
        EXPECT_DOUBLE_EQ(direct.real(), viaW.real())
            << "Re erfcx(z) != Re w(i z) at z = (" << z.real() << "," << z.imag() << ")";
        EXPECT_DOUBLE_EQ(direct.imag(), viaW.imag())
            << "Im erfcx(z) != Im w(i z) at z = (" << z.real() << "," << z.imag() << ")";
    }

    // The relerr argument is unused by the complex erfcx implementation.
    std::cerr << "    -> passing relerr = 0.1 must give the identical result\n";
    const FadComplex a = Faddeeva::erfcx(FadComplex(1.0, 1.0), 0.0);
    const FadComplex b = Faddeeva::erfcx(FadComplex(1.0, 1.0), 0.1);
    EXPECT_DOUBLE_EQ(a.real(), b.real()) << "erfcx(complex) should ignore relerr";
    EXPECT_DOUBLE_EQ(a.imag(), b.imag()) << "erfcx(complex) should ignore relerr";
}

// -----------------------------------------------------------------------------
// 7. Faddeeva::erf(double).
// -----------------------------------------------------------------------------
void faddeeva_test_erf_real()
{
    std::cerr << "\n[TEST] faddeeva_test_erf_real\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::erf(double)\n"
              << "  Checks: exact zero at the origin, reference value at x=1,\n"
              << "  oddness, saturation at large |x| and agreement with std::erf\n"
              << "  (independent <cmath> implementation) to 1e-12 relative.\n";

    std::cerr << "    -> erf(0) must be exactly 0\n";
    EXPECT_DOUBLE_EQ(Faddeeva::erf(0.0), 0.0) << "erf(0) must be 0";

    fad_expect_scalar("erf(1)", 0.84270079294971486934122063508260925929606699796630,
        Faddeeva::erf(1.0), 1e-13);
    fad_expect_scalar("erf(-1)", -0.84270079294971486934122063508260925929606699796630,
        Faddeeva::erf(-1.0), 1e-13);

    // Independent cross-check on a grid that spans the |x| < 5e-3 Taylor branch
    // and the general branch on both sides of the origin.
    std::cerr << "    -> grid comparison against std::erf\n";
    const std::vector<double> grid = { -6.0, -3.0, -1.0, -0.5, -0.004, -1e-8, 1e-8,
        0.004, 0.5, 1.0, 3.0, 6.0 };
    for (double x : grid) {
        fad_expect_scalar("erf(" + std::to_string(x) + ")", std::erf(x),
            Faddeeva::erf(x), 1e-12);
    }

    // Saturation / underflow handling.
    std::cerr << "    -> saturation at very large |x| and at infinity\n";
    EXPECT_DOUBLE_EQ(Faddeeva::erf(30.0), 1.0) << "erf(30) must saturate to 1";
    EXPECT_DOUBLE_EQ(Faddeeva::erf(-30.0), -1.0) << "erf(-30) must saturate to -1";
    EXPECT_DOUBLE_EQ(Faddeeva::erf(std::numeric_limits<double>::infinity()), 1.0)
        << "erf(+inf) must be 1";
    EXPECT_DOUBLE_EQ(Faddeeva::erf(-std::numeric_limits<double>::infinity()), -1.0)
        << "erf(-inf) must be -1";
}

// -----------------------------------------------------------------------------
// 8. Faddeeva::erf(complex).
// -----------------------------------------------------------------------------
void faddeeva_test_erf_complex()
{
    std::cerr << "\n[TEST] faddeeva_test_erf_complex\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::erf(std::complex<double>, double)\n"
              << "  Reference points cover the general branch, the small-|z| Taylor\n"
              << "  branch, the 'taylor_erfi' branch (small x, small x*y) and the\n"
              << "  pure real / pure imaginary axes.\n"
              << "  Pass criterion: componentwise relative error < 1e-13.\n";

    struct Sample {
        const char* name;
        FadComplex z;
        FadComplex e;
    };

    const std::vector<Sample> samples = {
        { "erf(1+2i)      [general branch, x>0]", FadComplex(1.0, 2.0),
            FadComplex(-0.5366435657785650339917955593141927494421,
                -5.049143703447034669543036958614140565553) },
        { "erf(-1+2i)     [general branch, x<0]", FadComplex(-1.0, 2.0),
            FadComplex(0.5366435657785650339917955593141927494421,
                -5.049143703447034669543036958614140565553) },
        { "erf(5.1e-3+1e-8i) [just above the Taylor cutoff]",
            FadComplex(5.1e-3, 1e-8),
            FadComplex(0.005754683859034800134412990541076554934877,
                0.1128349818335058741511924929801267822634e-7) },
        { "erf(-4.9e-3+4.95e-3i) [small-|z| Taylor branch]",
            FadComplex(-4.9e-3, 4.95e-3),
            FadComplex(-0.005529149142341821193633460286828381876955,
                0.005585388387864706679609092447916333443570) },
        { "erf(4.9e-3+0.5i)  [taylor_erfi branch]", FadComplex(4.9e-3, 0.5),
            FadComplex(0.007099365669981359632319829148438283865814,
                0.6149347012854211635026981277569074001219) },
        { "erf(1e-6+2e-6i)   [tiny-|z| Taylor branch]", FadComplex(1e-6, 2e-6),
            FadComplex(0.1128379167099649964175513742247082845155e-5,
                0.2256758334191777400570377193451519478895e-5) },
        { "erf(0+2i)         [imaginary axis]", FadComplex(0.0, 2.0),
            FadComplex(0.0, 18.56480241457555259870429191324101719886) },
    };

    for (const auto& s : samples) {
        fad_expect_complex(s.name, s.e, Faddeeva::erf(s.z, 0.0), 1e-13);
    }

    // On the real axis the complex routine short-circuits to the real routine.
    std::cerr << "    -> erf(x + 0i) must reproduce the real erf(x) exactly\n";
    for (double x : { -2.0, -0.001, 0.0, 0.001, 2.0 }) {
        const FadComplex got = Faddeeva::erf(FadComplex(x, 0.0), 0.0);
        EXPECT_DOUBLE_EQ(got.real(), Faddeeva::erf(x))
            << "Re erf(x+0i) must equal erf(x) at x = " << x;
        EXPECT_DOUBLE_EQ(got.imag(), 0.0)
            << "Im erf(x+0i) must vanish at x = " << x;
    }

    // Oddness of the complex error function: erf(-z) == -erf(z).
    std::cerr << "    -> erf(-z) must equal -erf(z)\n";
    for (const auto& z : { FadComplex(1.0, 2.0), FadComplex(0.3, -0.4) }) {
        fad_expect_complex("        erf(-z) vs -erf(z)", -Faddeeva::erf(z, 0.0),
            Faddeeva::erf(-z, 0.0), 1e-13);
    }
}

// -----------------------------------------------------------------------------
// 9. Faddeeva::erfc(double) and the complementary identity erf + erfc == 1.
// -----------------------------------------------------------------------------
void faddeeva_test_erfc_real()
{
    std::cerr << "\n[TEST] faddeeva_test_erfc_real\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::erfc(double)\n"
              << "  Checks: erfc(0)=1, reference value at x=1, agreement with\n"
              << "  std::erfc, the identity erf(x)+erfc(x)=1, and the underflow /\n"
              << "  saturation limits at |x| > 27.\n";

    std::cerr << "    -> erfc(0) must be exactly 1\n";
    EXPECT_DOUBLE_EQ(Faddeeva::erfc(0.0), 1.0) << "erfc(0) must be 1";

    fad_expect_scalar("erfc(1)", 0.15729920705028513065877936491739074070393300203369,
        Faddeeva::erfc(1.0), 1e-13);

    std::cerr << "    -> grid comparison against std::erfc\n";
    for (double x : { -3.0, -1.0, -0.25, 0.0, 0.25, 1.0, 3.0, 10.0, 20.0 }) {
        fad_expect_scalar("erfc(" + std::to_string(x) + ")", std::erfc(x),
            Faddeeva::erfc(x), 1e-11);
    }

    std::cerr << "    -> identity erf(x) + erfc(x) == 1\n";
    for (double x : { -2.0, -0.5, 0.0, 0.5, 2.0, 4.0 }) {
        const double sum = Faddeeva::erf(x) + Faddeeva::erfc(x);
        EXPECT_NEAR(sum, 1.0, 1e-14)
            << "erf(x)+erfc(x) must be 1; failed at x = " << x;
    }

    std::cerr << "    -> underflow (large +x) and saturation (large -x)\n";
    EXPECT_DOUBLE_EQ(Faddeeva::erfc(30.0), 0.0)
        << "erfc(30) must underflow to 0";
    EXPECT_DOUBLE_EQ(Faddeeva::erfc(-30.0), 2.0)
        << "erfc(-30) must saturate to 2";
}

// -----------------------------------------------------------------------------
// 10. Faddeeva::erfc(complex).
// -----------------------------------------------------------------------------
void faddeeva_test_erfc_complex()
{
    std::cerr << "\n[TEST] faddeeva_test_erfc_complex\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::erfc(std::complex<double>, double)\n"
              << "  Checks: reference points (general branch, both axes, underflow)\n"
              << "  and the identity erf(z) + erfc(z) == 1.\n";

    struct Sample {
        const char* name;
        FadComplex z;
        FadComplex e;
    };

    const std::vector<Sample> samples = {
        { "erfc(1+2i)   [general branch]", FadComplex(1.0, 2.0),
            FadComplex(1.536643565778565033991795559314192749442,
                5.049143703447034669543036958614140565553) },
        { "erfc(0+2i)   [imaginary axis]", FadComplex(0.0, 2.0),
            FadComplex(1.0, -18.56480241457555259870429191324101719886) },
        { "erfc(2+0i)   [real axis]", FadComplex(2.0, 0.0),
            FadComplex(0.004677734981047265837930743632747071389108, 0.0) },
        { "erfc(88+0i)  [underflow to zero]", FadComplex(88.0, 0.0),
            FadComplex(0.0, 0.0) },
    };

    for (const auto& s : samples) {
        const FadComplex got = Faddeeva::erfc(s.z, 0.0);
        // The underflow entry is exactly (0,0); a relative test is meaningless
        // there, so it is checked with an absolute tolerance.
        if (s.e == FadComplex(0.0, 0.0)) {
            std::cerr << "      " << s.name << ": got (" << got.real() << ", "
                      << got.imag() << "i), expected (0, 0i)\n";
            EXPECT_NEAR(got.real(), 0.0, 1e-300) << s.name << " [real part]";
            EXPECT_NEAR(got.imag(), 0.0, 1e-300) << s.name << " [imag part]";
        } else {
            fad_expect_complex(s.name, s.e, got, 1e-13);
        }
    }

    std::cerr << "    -> identity erf(z) + erfc(z) == 1 for complex z\n";
    for (const auto& z : { FadComplex(1.0, 2.0), FadComplex(-0.5, 0.5),
             FadComplex(0.0, 1.0), FadComplex(3.0, 0.0) }) {
        const FadComplex sum = Faddeeva::erf(z, 0.0) + Faddeeva::erfc(z, 0.0);
        const double scale = std::max(1.0, std::abs(Faddeeva::erf(z, 0.0)));
        std::cerr << "      z = (" << z.real() << ", " << z.imag() << "i): sum = ("
                  << std::setprecision(12) << sum.real() << ", " << sum.imag() << "i)\n";
        EXPECT_NEAR(sum.real(), 1.0, 1e-13 * scale)
            << "Re[erf(z)+erfc(z)] must be 1";
        EXPECT_NEAR(sum.imag(), 0.0, 1e-13 * scale)
            << "Im[erf(z)+erfc(z)] must be 0";
    }
}

// -----------------------------------------------------------------------------
// 11. Faddeeva::erfi (real and complex): erfi(z) = -i erf(i z).
// -----------------------------------------------------------------------------
void faddeeva_test_erfi()
{
    std::cerr << "\n[TEST] faddeeva_test_erfi\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Functions:   Faddeeva::erfi(double) and\n"
              << "               Faddeeva::erfi(std::complex<double>, double)\n"
              << "  Checks: reference values, oddness, the definition\n"
              << "  erfi(x) = exp(x^2) w_im(x), and erfi(z) = -i erf(i z).\n";

    // Reference value.
    fad_expect_scalar("erfi(1)", 1.65042575879754287602533772956136462951,
        Faddeeva::erfi(1.0), 1e-13);

    std::cerr << "    -> erfi(0) must be exactly 0 and erfi must be odd\n";
    EXPECT_DOUBLE_EQ(Faddeeva::erfi(0.0), 0.0) << "erfi(0) must be 0";
    for (double x : { 0.25, 1.0, 2.0 }) {
        EXPECT_DOUBLE_EQ(Faddeeva::erfi(-x), -Faddeeva::erfi(x))
            << "erfi must be odd; failed at x = " << x;
    }

    // The real implementation is exactly exp(x^2)*w_im(x) below the overflow cut.
    std::cerr << "    -> erfi(x) == exp(x*x)*w_im(x) for moderate x\n";
    for (double x : { 0.5, 1.0, 2.5 }) {
        EXPECT_DOUBLE_EQ(Faddeeva::erfi(x), std::exp(x * x) * Faddeeva::w_im(x))
            << "erfi(x) must equal exp(x*x)*w_im(x) at x = " << x;
    }

    // Overflow guard: for x*x > 720 the routine returns +/- infinity.
    std::cerr << "    -> erfi(+/-40) must overflow to +/- infinity\n";
    EXPECT_TRUE(std::isinf(Faddeeva::erfi(40.0)) && Faddeeva::erfi(40.0) > 0)
        << "erfi(40) must be +inf";
    EXPECT_TRUE(std::isinf(Faddeeva::erfi(-40.0)) && Faddeeva::erfi(-40.0) < 0)
        << "erfi(-40) must be -inf";

    // Complex reference value from the table in Faddeeva.cpp.
    fad_expect_complex("erfi(1.234+0.5678i)",
        FadComplex(1.081032284405373149432716643834106923212,
            1.926775520840916645838949402886591180834),
        Faddeeva::erfi(FadComplex(1.234, 0.5678), 0.0), 1e-13);

    // Definition check: erfi(z) == -i * erf(i z).
    std::cerr << "    -> erfi(z) == -i*erf(i*z)\n";
    for (const auto& z : { FadComplex(1.0, 2.0), FadComplex(-0.3, 0.7) }) {
        const FadComplex iz(-z.imag(), z.real());
        const FadComplex e = Faddeeva::erf(iz, 0.0);
        const FadComplex ref(e.imag(), -e.real()); // -i * e
        fad_expect_complex("        erfi(z) vs -i erf(iz)", ref,
            Faddeeva::erfi(z, 0.0), 1e-13);
    }

    // The real and complex overloads must agree on the real axis.
    std::cerr << "    -> erfi(x + 0i) must reproduce the real erfi(x)\n";
    for (double x : { 0.5, 1.234 }) {
        const FadComplex got = Faddeeva::erfi(FadComplex(x, 0.0), 0.0);
        fad_expect_scalar("        Re erfi(" + std::to_string(x) + "+0i)",
            Faddeeva::erfi(x), got.real(), 1e-13);
        EXPECT_NEAR(got.imag(), 0.0, 1e-13 * std::fabs(got.real()))
            << "Im erfi(x+0i) must vanish at x = " << x;
    }
}

// -----------------------------------------------------------------------------
// 12. Faddeeva::Dawson(double).
// -----------------------------------------------------------------------------
void faddeeva_test_dawson_real()
{
    std::cerr << "\n[TEST] faddeeva_test_dawson_real\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::Dawson(double)\n"
              << "  Checks: the defining relation D(x) = sqrt(pi)/2 * w_im(x),\n"
              << "  reference values at x = 2e-6, 2, 20, 200, oddness, and the\n"
              << "  location/height of the global maximum near x = 0.9241.\n";

    // The implementation is literally spi2*w_im(x): verify that contract.
    std::cerr << "    -> D(x) must equal sqrt(pi)/2 * w_im(x) exactly\n";
    for (double x : { 0.1, 1.0, 5.0, -2.0 }) {
        EXPECT_DOUBLE_EQ(Faddeeva::Dawson(x), kFadSpi2 * Faddeeva::w_im(x))
            << "Dawson(x) must equal sqrt(pi)/2*w_im(x) at x = " << x;
    }

    // Reference values (real parts of the y == 0 entries of the Dawson table).
    fad_expect_scalar("Dawson(2e-6)", 0.1999999999994666666666675199999999990248e-5,
        Faddeeva::Dawson(2e-6), 1e-13);
    fad_expect_scalar("Dawson(2)", 0.3013403889237919660346644392864226952119,
        Faddeeva::Dawson(2.0), 1e-13);
    fad_expect_scalar("Dawson(20)", 0.02503136792640367194699495234782353186858,
        Faddeeva::Dawson(20.0), 1e-13);
    fad_expect_scalar("Dawson(200)", 0.002500031251171948248596912483183760683918,
        Faddeeva::Dawson(200.0), 1e-13);

    // Oddness and the value at the origin.
    std::cerr << "    -> D(0) = 0 and D(-x) = -D(x)\n";
    EXPECT_DOUBLE_EQ(Faddeeva::Dawson(0.0), 0.0) << "Dawson(0) must be 0";
    for (double x : { 0.5, 2.0, 30.0 }) {
        EXPECT_DOUBLE_EQ(Faddeeva::Dawson(-x), -Faddeeva::Dawson(x))
            << "Dawson must be odd; failed at x = " << x;
    }

    // The well-known global maximum: D(0.9241388730...) = 0.5410442246...
    std::cerr << "    -> global maximum D(0.9241388730) ~ 0.5410442246\n";
    fad_expect_scalar("Dawson(0.9241388730)", 0.5410442246, Faddeeva::Dawson(0.9241388730),
        1e-9);
    EXPECT_GT(Faddeeva::Dawson(0.9241388730), Faddeeva::Dawson(0.85))
        << "0.9241 should be past the rising flank of the Dawson maximum";
    EXPECT_GT(Faddeeva::Dawson(0.9241388730), Faddeeva::Dawson(1.0))
        << "0.9241 should be before the falling flank of the Dawson maximum";
}

// -----------------------------------------------------------------------------
// 13. Faddeeva::Dawson(complex).
// -----------------------------------------------------------------------------
void faddeeva_test_dawson_complex()
{
    std::cerr << "\n[TEST] faddeeva_test_dawson_complex\n"
              << "  Source file: src/math/Faddeeva.cpp\n"
              << "  Function:    Faddeeva::Dawson(std::complex<double>, double)\n"
              << "  Reference points exercise the general branch, the small-|z|\n"
              << "  Taylor branch, both near-real-axis continued-fraction branches\n"
              << "  (|x| > 40 and |x| > 5e7), and the two coordinate axes.\n"
              << "  Pass criterion: componentwise relative error < 1e-12.\n";

    struct Sample {
        const char* name;
        FadComplex z;
        FadComplex d;
    };

    const std::vector<Sample> samples = {
        { "Dawson(2+1i)      [general branch]", FadComplex(2.0, 1.0),
            FadComplex(0.1635394094345355614904345232875688576839,
                -0.1531245755371229803585918112683241066853) },
        { "Dawson(0.5+4.9e-3i) [near-real-axis Taylor]", FadComplex(0.5, 4.9e-3),
            FadComplex(0.4244534840871830045021143490355372016428,
                0.002820278933186814021399602648373095266538) },
        { "Dawson(1e-6+2e-6i)  [small-|z| Taylor]", FadComplex(1e-6, 2e-6),
            FadComplex(0.1000000000007333333333344266666666664457e-5,
                0.2000000000001333333333323199999999978819e-5) },
        { "Dawson(39+6.4e-5i)  [|x|>40 continued fraction]", FadComplex(39.0, 6.4e-5),
            FadComplex(0.01282473148489433743567240624939698290584,
                -0.2105957276516618621447832572909153498104e-7) },
        { "Dawson(4.9e7+5e-11i)[|x|>5e7 continued fraction]", FadComplex(4.9e7, 5e-11),
            FadComplex(0.1020408163265306334945473399689037886997e-7,
                -0.1041232819658476285651490827866174985330e-25) },
        { "Dawson(0-2i)        [imaginary axis]", FadComplex(0.0, -2.0),
            FadComplex(0.0, -48.16001211429122974789822893525016528191) },
        { "Dawson(2e-6+0i)     [real axis]", FadComplex(2e-6, 0.0),
            FadComplex(0.1999999999994666666666675199999999990248e-5, 0.0) },
        { "Dawson(2+0i)        [real axis]", FadComplex(2.0, 0.0),
            FadComplex(0.3013403889237919660346644392864226952119, 0.0) },
    };

    for (const auto& s : samples) {
        fad_expect_complex(s.name, s.d, Faddeeva::Dawson(s.z, 0.0), 1e-12);
    }

    // On the real axis the complex routine short-circuits to spi2*w_im(x).
    std::cerr << "    -> Dawson(x + 0i) must reproduce the real Dawson(x) exactly\n";
    for (double x : { -3.0, -0.25, 0.0, 0.25, 3.0 }) {
        const FadComplex got = Faddeeva::Dawson(FadComplex(x, 0.0), 0.0);
        EXPECT_DOUBLE_EQ(got.real(), Faddeeva::Dawson(x))
            << "Re Dawson(x+0i) must equal Dawson(x) at x = " << x;
        EXPECT_DOUBLE_EQ(std::fabs(got.imag()), 0.0)
            << "Im Dawson(x+0i) must vanish at x = " << x;
    }

    // Oddness of the complex Dawson function: D(-z) == -D(z).
    std::cerr << "    -> Dawson(-z) must equal -Dawson(z)\n";
    for (const auto& z : { FadComplex(2.0, 1.0), FadComplex(0.4, -0.9) }) {
        fad_expect_complex("        Dawson(-z) vs -Dawson(z)",
            -Faddeeva::Dawson(z, 0.0), Faddeeva::Dawson(-z, 0.0), 1e-12);
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* helper is invoked from its own TEST so
// that a failure in one group does not prevent the remaining groups from running
// (all assertions used above are non-fatal EXPECT_* macros).
// -----------------------------------------------------------------------------
TEST(FaddeevaTest, WKnownValues) { faddeeva_test_w_known_values(); }
TEST(FaddeevaTest, WAxesAndSymmetry) { faddeeva_test_w_axes_and_symmetry(); }
TEST(FaddeevaTest, WRelerrParameter) { faddeeva_test_w_relerr_parameter(); }
TEST(FaddeevaTest, WImaginaryPartRealArg) { faddeeva_test_w_im(); }
TEST(FaddeevaTest, ErfcxReal) { faddeeva_test_erfcx_real(); }
TEST(FaddeevaTest, ErfcxComplex) { faddeeva_test_erfcx_complex(); }
TEST(FaddeevaTest, ErfReal) { faddeeva_test_erf_real(); }
TEST(FaddeevaTest, ErfComplex) { faddeeva_test_erf_complex(); }
TEST(FaddeevaTest, ErfcReal) { faddeeva_test_erfc_real(); }
TEST(FaddeevaTest, ErfcComplex) { faddeeva_test_erfc_complex(); }
TEST(FaddeevaTest, Erfi) { faddeeva_test_erfi(); }
TEST(FaddeevaTest, DawsonReal) { faddeeva_test_dawson_real(); }
TEST(FaddeevaTest, DawsonComplex) { faddeeva_test_dawson_complex(); }