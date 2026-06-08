#include "reactions/association/association.hpp"

#include <cmath>

#include <gtest/gtest.h>

namespace {

constexpr double kTolerance = 1.0e-12;

double Pi() { return std::acos(-1.0); }

Vector VectorWithMagnitude(double x, double y, double z) {
  Vector vector{x, y, z};
  vector.calc_magnitude();
  return vector;
}

TEST(AssociationAngleHelpersTest, ComparesAnglesWithLegacyTolerance) {
  EXPECT_TRUE(areSameAngle(1.0, 1.0 + 5.0e-5));
  EXPECT_FALSE(areSameAngle(1.0, 1.0 + 2.0e-4));
  EXPECT_FALSE(areSameAngle(Pi(), -Pi()));
}

TEST(AssociationAngleHelpersTest, DetectsOnlyExactParallelSentinels) {
  EXPECT_TRUE(areParallel(0.0));
  EXPECT_TRUE(areParallel(Pi()));
  EXPECT_FALSE(areParallel(Pi() - 1.0e-12));
}

TEST(AssociationAngleHelpersTest, ChecksAngleSignInXYAndXZPlanes) {
  EXPECT_TRUE(angleSignIsCorrect(Vector{1.0, 0.0, 0.0},
                                 Vector{0.0, -1.0, 0.0}));
  EXPECT_FALSE(angleSignIsCorrect(Vector{1.0, 0.0, 0.0},
                                  Vector{0.0, 1.0, 0.0}));

  EXPECT_TRUE(angleSignIsCorrect(Vector{1.0, 0.0, 1.0},
                                 Vector{0.0, 0.0, 1.0}));
  EXPECT_FALSE(angleSignIsCorrect(Vector{0.0, 0.0, 1.0},
                                  Vector{1.0, 0.0, 1.0}));
}

TEST(AssociationAngleHelpersTest, CreatesOrthogonalFallbackVector) {
  Vector x_axis = VectorWithMagnitude(1.0, 0.0, 0.0);
  Vector fallback = create_arbitrary_vector(x_axis);
  EXPECT_NEAR(fallback.dot(x_axis), 0.0, kTolerance);
  EXPECT_NEAR(fallback.magnitude, 1.0, kTolerance);

  Vector diagonal = VectorWithMagnitude(1.0, 1.0, 0.0);
  Vector orthogonal = create_arbitrary_vector(diagonal);
  EXPECT_NEAR(orthogonal.dot(diagonal), 0.0, kTolerance);
  EXPECT_NEAR(orthogonal.magnitude, 1.0, kTolerance);
}

TEST(AssociationAngleHelpersTest, PreservesRequiresSignFlipProjectionCases) {
  EXPECT_FALSE(requiresSignFlip(VectorWithMagnitude(0.0, 0.0, 1.0),
                                VectorWithMagnitude(1.0, 0.0, 0.0),
                                VectorWithMagnitude(0.0, 1.0, 0.0)));
  EXPECT_TRUE(requiresSignFlip(VectorWithMagnitude(0.0, 0.0, 1.0),
                               VectorWithMagnitude(1.0, 0.0, 0.0),
                               VectorWithMagnitude(0.0, -1.0, 0.0)));

  EXPECT_TRUE(requiresSignFlip(VectorWithMagnitude(1.0, 0.0, 0.0),
                               VectorWithMagnitude(0.0, 1.0, 0.0),
                               VectorWithMagnitude(0.0, 0.0, 1.0)));
  EXPECT_FALSE(requiresSignFlip(VectorWithMagnitude(1.0, 0.0, 0.0),
                                VectorWithMagnitude(0.0, 1.0, 0.0),
                                VectorWithMagnitude(0.0, 0.0, -1.0)));
}

}  // namespace
