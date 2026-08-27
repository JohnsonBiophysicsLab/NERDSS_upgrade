/*! \file test_math_functions.cpp
 *
 * ### Unit tests for ../src/math/math_functions.cpp
 *
 * The translation unit under test provides three functions inside the
 * `MathFuncs` namespace:
 *
 *   - `long double MathFuncs::factorial(unsigned n)`
 *        A simple recursive factorial returning a `long double`.
 *
 *   - `double MathFuncs::gammln(double n)`
 *        The Numerical Recipes (ch. 6) Lanczos approximation of
 *        \f$\ln \Gamma(n)\f$.  Accurate to roughly 1e-10 relative for
 *        positive arguments.
 *
 *   - `double MathFuncs::gammFactorial(int n)`
 *        Returns \f$n!\f$.  For `n <= 32` the value is taken from (and
 *        progressively cached into) a `static float` lookup table; for
 *        `n > 32` it falls back on `exp(gammln(n+1))`.
 *
 * IMPORTANT NOTES ABOUT WHAT IS *NOT* TESTED HERE:
 *   - `gammFactorial()` calls `exit(1)` when handed a negative argument.
 *     Exercising that path would tear down the whole gtest binary, so the
 *     negative-input branch is deliberately never invoked.
 *   - `gammln()` is only mathematically defined for arguments > 0; passing
 *     0 divides by zero inside the routine, so only positive arguments are
 *     used below.
 *   - Because the lookup table inside `gammFactorial()` is `static float`,
 *     values with more than 24 significant bits (roughly n >= 14) carry a
 *     small rounding error.  Assertions therefore use *relative*
 *     tolerances rather than exact equality for those inputs.
 */

#include "math/math_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace so they cannot collide with other TUs).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Compute the relative difference |a-b| / max(1,|b|).
 *
 * Used so that assertions on very large factorials (up to ~1e35) can be
 * expressed as a fractional error rather than an absolute one.
 */
double mathfuncs_rel_err(double a, double b)
{
    const double denom = (std::fabs(b) > 1.0) ? std::fabs(b) : 1.0;
    return std::fabs(a - b) / denom;
}

/*! \brief Reference factorial computed iteratively in long double. */
long double mathfuncs_ref_factorial(unsigned n)
{
    long double result = 1.0L;
    for (unsigned i = 2; i <= n; ++i)
        result *= static_cast<long double>(i);
    return result;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: MathFuncs::factorial() base cases.
// -----------------------------------------------------------------------------
void test_mathfuncs_factorial_base_cases()
{
    std::cerr << "\n[TEST] test_mathfuncs_factorial_base_cases\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::factorial(unsigned)\n"
              << "  Scenario:    the recursion terminating cases n = 0 and n = 1.\n"
              << "  Criteria:    0! == 1 and 1! == 1 exactly.\n";

    const long double f0 = MathFuncs::factorial(0);
    const long double f1 = MathFuncs::factorial(1);

    std::cerr << "  factorial(0) = " << static_cast<double>(f0) << '\n';
    std::cerr << "  factorial(1) = " << static_cast<double>(f1) << '\n';

    // Both of these are small integers, so exact double comparison is valid.
    EXPECT_DOUBLE_EQ(static_cast<double>(f0), 1.0)
        << "factorial(0) must be exactly 1 (recursion terminator)";
    EXPECT_DOUBLE_EQ(static_cast<double>(f1), 1.0)
        << "factorial(1) must be exactly 1";
}

// -----------------------------------------------------------------------------
// Test 2: MathFuncs::factorial() against hand-written known values.
// -----------------------------------------------------------------------------
void test_mathfuncs_factorial_known_values()
{
    std::cerr << "\n[TEST] test_mathfuncs_factorial_known_values\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::factorial(unsigned)\n"
              << "  Scenario:    n = 2..12, all of which are exactly\n"
              << "               representable in a double (12! = 479001600).\n"
              << "  Criteria:    exact equality with the tabulated value.\n";

    // Hand-tabulated factorials 2! through 12!.
    const double expected[] = {
        2.0,        // 2!
        6.0,        // 3!
        24.0,       // 4!
        120.0,      // 5!
        720.0,      // 6!
        5040.0,     // 7!
        40320.0,    // 8!
        362880.0,   // 9!
        3628800.0,  // 10!
        39916800.0, // 11!
        479001600.0 // 12!
    };

    for (unsigned n = 2; n <= 12; ++n) {
        const double got = static_cast<double>(MathFuncs::factorial(n));
        const double want = expected[n - 2];
        std::cerr << "  factorial(" << n << ") = " << got
                  << " (expected " << want << ")\n";
        EXPECT_DOUBLE_EQ(got, want)
            << "factorial(" << n << ") should be exactly " << want;
    }
}

// -----------------------------------------------------------------------------
// Test 3: MathFuncs::factorial() for larger n, and the defining recurrence.
// -----------------------------------------------------------------------------
void test_mathfuncs_factorial_recurrence_and_large()
{
    std::cerr << "\n[TEST] test_mathfuncs_factorial_recurrence_and_large\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::factorial(unsigned)\n"
              << "  Scenario:    (a) verify n! == n * (n-1)! for n = 2..20,\n"
              << "               (b) compare n = 13..20 against an iterative\n"
              << "                   long-double reference.\n"
              << "  Criteria:    relative error below 1e-15 (values remain\n"
              << "               integral and exactly representable).\n";

    // (a) The recurrence relation the implementation is built on.
    for (unsigned n = 2; n <= 20; ++n) {
        const long double fn = MathFuncs::factorial(n);
        const long double fnm1 = MathFuncs::factorial(n - 1);
        const double lhs = static_cast<double>(fn);
        const double rhs = static_cast<double>(static_cast<long double>(n) * fnm1);
        EXPECT_LT(mathfuncs_rel_err(lhs, rhs), 1e-15)
            << "factorial(" << n << ") should equal " << n << " * factorial("
            << n - 1 << ")";
    }
    std::cerr << "  Recurrence n! == n*(n-1)! verified for n = 2..20\n";

    // (b) Larger values checked against an independent iterative computation.
    for (unsigned n = 13; n <= 20; ++n) {
        const double got = static_cast<double>(MathFuncs::factorial(n));
        const double want = static_cast<double>(mathfuncs_ref_factorial(n));
        std::cerr << "  factorial(" << n << ") = " << got
                  << " (reference " << want << ")\n";
        EXPECT_LT(mathfuncs_rel_err(got, want), 1e-15)
            << "factorial(" << n << ") disagrees with the iterative reference";
    }
}

// -----------------------------------------------------------------------------
// Test 4: MathFuncs::gammln() against closed-form values of ln(Gamma(x)).
// -----------------------------------------------------------------------------
void test_mathfuncs_gammln_known_values()
{
    std::cerr << "\n[TEST] test_mathfuncs_gammln_known_values\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::gammln(double)\n"
              << "  Scenario:    integer and half-integer arguments where\n"
              << "               ln(Gamma(x)) has a closed form.\n"
              << "  Criteria:    absolute error below 1e-9 (the Lanczos\n"
              << "               series used here is good to ~2e-10).\n";

    const double tol = 1e-9;

    // Gamma(1) = 0! = 1  ->  ln = 0
    std::cerr << "  gammln(1.0) = " << MathFuncs::gammln(1.0) << " (expect 0)\n";
    EXPECT_NEAR(MathFuncs::gammln(1.0), 0.0, tol)
        << "ln(Gamma(1)) = ln(1) = 0";

    // Gamma(2) = 1! = 1  ->  ln = 0
    std::cerr << "  gammln(2.0) = " << MathFuncs::gammln(2.0) << " (expect 0)\n";
    EXPECT_NEAR(MathFuncs::gammln(2.0), 0.0, tol)
        << "ln(Gamma(2)) = ln(1) = 0";

    // Gamma(3) = 2! = 2  ->  ln = ln(2)
    std::cerr << "  gammln(3.0) = " << MathFuncs::gammln(3.0)
              << " (expect " << std::log(2.0) << ")\n";
    EXPECT_NEAR(MathFuncs::gammln(3.0), std::log(2.0), tol)
        << "ln(Gamma(3)) = ln(2)";

    // Gamma(4) = 3! = 6  ->  ln = ln(6)
    std::cerr << "  gammln(4.0) = " << MathFuncs::gammln(4.0)
              << " (expect " << std::log(6.0) << ")\n";
    EXPECT_NEAR(MathFuncs::gammln(4.0), std::log(6.0), tol)
        << "ln(Gamma(4)) = ln(6)";

    // Gamma(1/2) = sqrt(pi)  ->  ln = 0.5*ln(pi)
    const double halfLogPi = 0.5 * std::log(M_PI);
    std::cerr << "  gammln(0.5) = " << MathFuncs::gammln(0.5)
              << " (expect " << halfLogPi << ")\n";
    EXPECT_NEAR(MathFuncs::gammln(0.5), halfLogPi, tol)
        << "ln(Gamma(0.5)) = ln(sqrt(pi))";

    // Gamma(3/2) = sqrt(pi)/2  ->  ln = 0.5*ln(pi) - ln(2)
    const double expect15 = halfLogPi - std::log(2.0);
    std::cerr << "  gammln(1.5) = " << MathFuncs::gammln(1.5)
              << " (expect " << expect15 << ")\n";
    EXPECT_NEAR(MathFuncs::gammln(1.5), expect15, tol)
        << "ln(Gamma(1.5)) = ln(sqrt(pi)/2)";
}

// -----------------------------------------------------------------------------
// Test 5: MathFuncs::gammln() recurrence and agreement with std::lgamma.
// -----------------------------------------------------------------------------
void test_mathfuncs_gammln_recurrence_and_stdlib()
{
    std::cerr << "\n[TEST] test_mathfuncs_gammln_recurrence_and_stdlib\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::gammln(double)\n"
              << "  Scenario:    (a) the functional equation\n"
              << "                   ln Gamma(x+1) - ln Gamma(x) == ln(x),\n"
              << "               (b) cross-check against std::lgamma over a\n"
              << "                   range of arguments including large ones.\n"
              << "  Criteria:    relative error below 1e-9.\n";

    // (a) Functional equation Gamma(x+1) = x * Gamma(x).
    const double xs[] = { 0.25, 0.5, 1.0, 2.7, 5.0, 12.5, 40.0 };
    for (double x : xs) {
        const double diff = MathFuncs::gammln(x + 1.0) - MathFuncs::gammln(x);
        std::cerr << "  gammln(" << x + 1.0 << ") - gammln(" << x << ") = "
                  << diff << " (expect ln(" << x << ") = " << std::log(x) << ")\n";
        EXPECT_NEAR(diff, std::log(x), 1e-9)
            << "Functional equation failed at x = " << x;
    }

    // (b) Compare directly to the standard library implementation. lgamma()
    //     returns ln|Gamma(x)|, which equals ln Gamma(x) for x > 0.
    const double xs2[] = { 0.5, 1.0, 2.0, 3.3, 10.0, 50.0, 100.0, 500.0 };
    for (double x : xs2) {
        const double got = MathFuncs::gammln(x);
        const double want = std::lgamma(x);
        std::cerr << "  gammln(" << x << ") = " << got
                  << " vs std::lgamma = " << want << '\n';
        EXPECT_LT(mathfuncs_rel_err(got, want), 1e-9)
            << "gammln(" << x << ") disagrees with std::lgamma";
    }
}

// -----------------------------------------------------------------------------
// Test 6: MathFuncs::gammFactorial() for the pre-initialised table entries.
// -----------------------------------------------------------------------------
void test_mathfuncs_gammfactorial_preloaded_entries()
{
    std::cerr << "\n[TEST] test_mathfuncs_gammfactorial_preloaded_entries\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::gammFactorial(int)\n"
              << "  Scenario:    n = 0..4, which are hard-coded in the static\n"
              << "               lookup table {1, 1, 2, 6, 24}.\n"
              << "  Criteria:    exact equality (all are small integers).\n";

    const double expected[] = { 1.0, 1.0, 2.0, 6.0, 24.0 };
    for (int n = 0; n <= 4; ++n) {
        const double got = MathFuncs::gammFactorial(n);
        std::cerr << "  gammFactorial(" << n << ") = " << got
                  << " (expected " << expected[n] << ")\n";
        EXPECT_DOUBLE_EQ(got, expected[n])
            << "gammFactorial(" << n << ") should return the pre-loaded value";
    }
}

// -----------------------------------------------------------------------------
// Test 7: MathFuncs::gammFactorial() table extension (5 <= n <= 32).
// -----------------------------------------------------------------------------
void test_mathfuncs_gammfactorial_table_extension()
{
    std::cerr << "\n[TEST] test_mathfuncs_gammfactorial_table_extension\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::gammFactorial(int)\n"
              << "  Scenario:    n = 5..32, values produced by the internal\n"
              << "               while-loop that extends the static table.\n"
              << "  Criteria:    n <= 13 is exact in a float, so use exact\n"
              << "               equality there; beyond that the float table\n"
              << "               accumulates rounding, so require only a 1e-5\n"
              << "               relative agreement with the long-double\n"
              << "               reference factorial.\n";

    for (int n = 5; n <= 32; ++n) {
        const double got = MathFuncs::gammFactorial(n);
        const double want = static_cast<double>(
            mathfuncs_ref_factorial(static_cast<unsigned>(n)));

        if (n <= 13) {
            // 13! = 6227020800 needs 23 mantissa bits -> exact as a float.
            std::cerr << "  gammFactorial(" << n << ") = " << got
                      << " (exact expected " << want << ")\n";
            EXPECT_DOUBLE_EQ(got, want)
                << "gammFactorial(" << n << ") should be exact for n <= 13";
        } else {
            const double rel = mathfuncs_rel_err(got, want);
            std::cerr << "  gammFactorial(" << n << ") = " << got
                      << " (reference " << want << ", rel err " << rel << ")\n";
            EXPECT_LT(rel, 1e-5)
                << "gammFactorial(" << n
                << ") drifts more than the float table should allow";
        }
    }
}

// -----------------------------------------------------------------------------
// Test 8: MathFuncs::gammFactorial() for n > 32 (the gammln fallback branch).
// -----------------------------------------------------------------------------
void test_mathfuncs_gammfactorial_large_branch()
{
    std::cerr << "\n[TEST] test_mathfuncs_gammfactorial_large_branch\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::gammFactorial(int)\n"
              << "  Scenario:    n > 32, where the implementation returns\n"
              << "               exp(gammln(n+1)) instead of a table lookup.\n"
              << "  Criteria:    relative agreement with exp(std::lgamma(n+1))\n"
              << "               better than 1e-8.\n";

    const int ns[] = { 33, 40, 50, 60, 80, 100 };
    for (int n : ns) {
        const double got = MathFuncs::gammFactorial(n);
        const double want = std::exp(std::lgamma(static_cast<double>(n) + 1.0));
        const double rel = mathfuncs_rel_err(got, want);
        std::cerr << "  gammFactorial(" << n << ") = " << got
                  << " (reference " << want << ", rel err " << rel << ")\n";
        EXPECT_GT(got, 0.0)
            << "gammFactorial(" << n << ") must be strictly positive";
        EXPECT_LT(rel, 1e-8)
            << "gammFactorial(" << n
            << ") from the exp(gammln()) branch disagrees with the reference";
    }

    // The two branches should meet smoothly at the n = 32/33 boundary:
    // 33! must equal 32! * 33 to within the float table's precision.
    const double f32 = MathFuncs::gammFactorial(32); // table branch
    const double f33 = MathFuncs::gammFactorial(33); // gammln branch
    const double rel = mathfuncs_rel_err(f33, f32 * 33.0);
    std::cerr << "  Boundary check: 33! (=" << f33 << ") vs 32!*33 (="
              << f32 * 33.0 << "), rel err " << rel << '\n';
    EXPECT_LT(rel, 1e-5)
        << "The table branch (n=32) and gammln branch (n=33) should be "
           "consistent to within float-table precision";
}

// -----------------------------------------------------------------------------
// Test 9: MathFuncs::gammFactorial() static-state safety.
// -----------------------------------------------------------------------------
void test_mathfuncs_gammfactorial_static_state()
{
    std::cerr << "\n[TEST] test_mathfuncs_gammfactorial_static_state\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Function:    MathFuncs::gammFactorial(int)\n"
              << "  Scenario:    the routine caches results in a function-local\n"
              << "               static table, so repeated and out-of-order\n"
              << "               calls must not corrupt earlier entries.\n"
              << "  Criteria:    (a) repeated calls are bit-identical,\n"
              << "               (b) a small n queried after a large n still\n"
              << "                   returns the correct value.\n";

    // (a) Idempotence: asking twice must give exactly the same answer.
    const double first = MathFuncs::gammFactorial(20);
    const double second = MathFuncs::gammFactorial(20);
    std::cerr << "  gammFactorial(20) called twice: " << first << " and "
              << second << '\n';
    EXPECT_DOUBLE_EQ(first, second)
        << "Repeated calls with the same n must return identical values";

    // (b) Warm the cache to a high index, then ask for a low index again.
    MathFuncs::gammFactorial(30);
    const double f6 = MathFuncs::gammFactorial(6);
    std::cerr << "  gammFactorial(6) after warming the table to n=30: " << f6
              << " (expected 720)\n";
    EXPECT_DOUBLE_EQ(f6, 720.0)
        << "Extending the static table must not corrupt earlier entries";

    const double f10 = MathFuncs::gammFactorial(10);
    std::cerr << "  gammFactorial(10) after warming the table: " << f10
              << " (expected 3628800)\n";
    EXPECT_DOUBLE_EQ(f10, 3628800.0)
        << "Extending the static table must not corrupt earlier entries";
}

// -----------------------------------------------------------------------------
// Test 10: cross-consistency between factorial() and gammFactorial().
// -----------------------------------------------------------------------------
void test_mathfuncs_factorial_vs_gammfactorial()
{
    std::cerr << "\n[TEST] test_mathfuncs_factorial_vs_gammfactorial\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Functions:   MathFuncs::factorial(unsigned) and\n"
              << "               MathFuncs::gammFactorial(int)\n"
              << "  Scenario:    both compute n!, so they must agree for\n"
              << "               n = 0..20.\n"
              << "  Criteria:    relative difference below 1e-5 (limited by\n"
              << "               the float lookup table used by\n"
              << "               gammFactorial()).\n";

    for (int n = 0; n <= 20; ++n) {
        const double viaRecursion =
            static_cast<double>(MathFuncs::factorial(static_cast<unsigned>(n)));
        const double viaTable = MathFuncs::gammFactorial(n);
        const double rel = mathfuncs_rel_err(viaTable, viaRecursion);
        std::cerr << "  n = " << n << ": factorial = " << viaRecursion
                  << ", gammFactorial = " << viaTable << ", rel err " << rel
                  << '\n';
        EXPECT_LT(rel, 1e-5)
            << "factorial(" << n << ") and gammFactorial(" << n
            << ") should describe the same quantity";
    }
}

// -----------------------------------------------------------------------------
// Test 11: gammln() as the log of the factorial, tying the two together.
// -----------------------------------------------------------------------------
void test_mathfuncs_gammln_matches_log_factorial()
{
    std::cerr << "\n[TEST] test_mathfuncs_gammln_matches_log_factorial\n"
              << "  Source file: src/math/math_functions.cpp\n"
              << "  Functions:   MathFuncs::gammln(double) and\n"
              << "               MathFuncs::factorial(unsigned)\n"
              << "  Scenario:    for integer n, ln Gamma(n+1) == ln(n!).\n"
              << "  Criteria:    absolute error below 1e-9.\n";

    for (unsigned n = 1; n <= 15; ++n) {
        const double lnFact =
            std::log(static_cast<double>(MathFuncs::factorial(n)));
        const double got = MathFuncs::gammln(static_cast<double>(n) + 1.0);
        std::cerr << "  gammln(" << n + 1 << ") = " << got
                  << " vs ln(" << n << "!) = " << lnFact << '\n';
        EXPECT_NEAR(got, lnFact, 1e-9)
            << "gammln(n+1) should equal ln(n!) for n = " << n;
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is executed inside its own TEST so a
// failure in one does not prevent the remaining checks from running.
// -----------------------------------------------------------------------------
TEST(MathFunctionsTest, FactorialBaseCases) { test_mathfuncs_factorial_base_cases(); }
TEST(MathFunctionsTest, FactorialKnownValues) { test_mathfuncs_factorial_known_values(); }
TEST(MathFunctionsTest, FactorialRecurrenceAndLarge) { test_mathfuncs_factorial_recurrence_and_large(); }
TEST(MathFunctionsTest, GammlnKnownValues) { test_mathfuncs_gammln_known_values(); }
TEST(MathFunctionsTest, GammlnRecurrenceAndStdlib) { test_mathfuncs_gammln_recurrence_and_stdlib(); }
TEST(MathFunctionsTest, GammFactorialPreloadedEntries) { test_mathfuncs_gammfactorial_preloaded_entries(); }
TEST(MathFunctionsTest, GammFactorialTableExtension) { test_mathfuncs_gammfactorial_table_extension(); }
TEST(MathFunctionsTest, GammFactorialLargeBranch) { test_mathfuncs_gammfactorial_large_branch(); }
TEST(MathFunctionsTest, GammFactorialStaticState) { test_mathfuncs_gammfactorial_static_state(); }
TEST(MathFunctionsTest, FactorialVsGammFactorial) { test_mathfuncs_factorial_vs_gammfactorial(); }
TEST(MathFunctionsTest, GammlnMatchesLogFactorial) { test_mathfuncs_gammln_matches_log_factorial(); }