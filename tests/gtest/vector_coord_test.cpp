#include "classes/class_Coord.hpp"
#include "classes/class_Vector.hpp"

#include <cmath>

#include <gtest/gtest.h>

namespace {

constexpr double kTolerance = 1.0e-12;

TEST(CoordTest, RoundsCoordinatesAndDetectsColinearity) {
  Coord rounded = round({1.23456, -1.23456, 0.00004});
  EXPECT_NEAR(rounded.x, 1.2346, kTolerance);
  EXPECT_NEAR(rounded.y, -1.2346, kTolerance);
  EXPECT_NEAR(rounded.z, 0.0, kTolerance);

  Coord first{0.0, 0.0, 0.0};
  Coord second{1.0, 1.0, 1.0};
  Coord third{2.0, 2.0, 2.0};
  EXPECT_TRUE(is_co_linear(first, second, third));
}

TEST(VectorTest, ComputesMagnitudeDotAndNormalization) {
  Vector vector{3.0, 4.0, 12.0};
  vector.calc_magnitude();
  EXPECT_NEAR(vector.magnitude, 13.0, kTolerance);

  Vector x_axis{1.0, 0.0, 0.0};
  Vector y_axis{0.0, 1.0, 0.0};
  EXPECT_NEAR(x_axis.dot(y_axis), 0.0, kTolerance);

  vector.normalize();
  EXPECT_NEAR(vector.magnitude, 1.0, kTolerance);
  EXPECT_NEAR(vector.x, 3.0 / 13.0, kTolerance);
  EXPECT_NEAR(vector.y, 4.0 / 13.0, kTolerance);
  EXPECT_NEAR(vector.z, 12.0 / 13.0, kTolerance);
}

TEST(VectorTest, ComputesCrossProjectionAndAngle) {
  Vector x_axis{1.0, 0.0, 0.0};
  Vector y_axis{0.0, 1.0, 0.0};
  Vector cross = x_axis.cross(y_axis);
  EXPECT_NEAR(cross.x, 0.0, kTolerance);
  EXPECT_NEAR(cross.y, 0.0, kTolerance);
  EXPECT_NEAR(cross.z, 1.0, kTolerance);
  EXPECT_NEAR(cross.magnitude, 1.0, kTolerance);

  Vector original{2.0, 3.0, 4.0};
  Vector normal{0.0, 0.0, 1.0};
  Vector projected = original.vector_projection(normal);
  EXPECT_NEAR(projected.x, 2.0, kTolerance);
  EXPECT_NEAR(projected.y, 3.0, kTolerance);
  EXPECT_NEAR(projected.z, 0.0, kTolerance);

  x_axis.calc_magnitude();
  y_axis.calc_magnitude();
  EXPECT_NEAR(x_axis.dot_theta(y_axis), std::acos(0.0), kTolerance);
}

}  // namespace
