#include "math/math_functions.hpp"

#include <cmath>

#include <gtest/gtest.h>

namespace {

constexpr double kTolerance = 1.0e-10;

TEST(MathFunctionsTest, ComputesFactorialBaseAndRecursiveCases) {
  EXPECT_EQ(MathFuncs::factorial(0), 1);
  EXPECT_EQ(MathFuncs::factorial(1), 1);
  EXPECT_EQ(MathFuncs::factorial(5), 120);
  EXPECT_EQ(MathFuncs::factorial(10), 3628800);
}

TEST(MathFunctionsTest, GammaLogMatchesStandardLibraryForPositiveInputs) {
  EXPECT_NEAR(MathFuncs::gammln(0.5), std::lgamma(0.5), kTolerance);
  EXPECT_NEAR(MathFuncs::gammln(1.0), std::lgamma(1.0), kTolerance);
  EXPECT_NEAR(MathFuncs::gammln(5.5), std::lgamma(5.5), kTolerance);
  EXPECT_NEAR(MathFuncs::gammln(33.0), std::lgamma(33.0), kTolerance);
}

TEST(MathFunctionsTest, GammaFactorialUsesCacheAndGammaFallback) {
  EXPECT_NEAR(MathFuncs::gammFactorial(0), 1.0, kTolerance);
  EXPECT_NEAR(MathFuncs::gammFactorial(1), 1.0, kTolerance);
  EXPECT_NEAR(MathFuncs::gammFactorial(5), 120.0, kTolerance);
  EXPECT_NEAR(MathFuncs::gammFactorial(10), 3628800.0, kTolerance);

  EXPECT_NEAR(MathFuncs::gammFactorial(33), std::tgamma(34.0),
              std::tgamma(34.0) * 1.0e-10);
}

}  // namespace
