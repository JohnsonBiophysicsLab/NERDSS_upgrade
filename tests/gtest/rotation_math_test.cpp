#include "classes/class_Quat.hpp"
#include "classes/class_Vector.hpp"
#include "math/matrix.hpp"

#include <array>
#include <cmath>

#include <gtest/gtest.h>

namespace {

constexpr double kTolerance = 1.0e-12;

double HalfPi() { return std::acos(-1.0) / 2.0; }

void ExpectVectorNear(const Vector& vector, double x, double y, double z) {
  EXPECT_NEAR(vector.x, x, kTolerance);
  EXPECT_NEAR(vector.y, y, kTolerance);
  EXPECT_NEAR(vector.z, z, kTolerance);
}

TEST(QuatTest, ComputesAlgebraAndUnitOperations) {
  Quat first{1.0, 2.0, 3.0, 4.0};
  Quat second{5.0, 6.0, 7.0, 8.0};

  Quat product = first * second;
  EXPECT_NEAR(product.w, -60.0, kTolerance);
  EXPECT_NEAR(product.x, 12.0, kTolerance);
  EXPECT_NEAR(product.y, 30.0, kTolerance);
  EXPECT_NEAR(product.z, 24.0, kTolerance);

  EXPECT_NEAR(first.norm(), 30.0, kTolerance);
  EXPECT_NEAR(first.mag(), std::sqrt(30.0), kTolerance);

  Quat scaled = first.scale(0.5);
  EXPECT_NEAR(scaled.w, 0.5, kTolerance);
  EXPECT_NEAR(scaled.x, 1.0, kTolerance);
  EXPECT_NEAR(scaled.y, 1.5, kTolerance);
  EXPECT_NEAR(scaled.z, 2.0, kTolerance);

  Quat conjugate = first.conjugate();
  EXPECT_NEAR(conjugate.w, 1.0, kTolerance);
  EXPECT_NEAR(conjugate.x, -2.0, kTolerance);
  EXPECT_NEAR(conjugate.y, -3.0, kTolerance);
  EXPECT_NEAR(conjugate.z, -4.0, kTolerance);

  Quat inverse = first.inverse();
  Quat identity = first * inverse;
  EXPECT_NEAR(identity.w, 1.0, kTolerance);
  EXPECT_NEAR(identity.x, 0.0, kTolerance);
  EXPECT_NEAR(identity.y, 0.0, kTolerance);
  EXPECT_NEAR(identity.z, 0.0, kTolerance);

  EXPECT_NEAR(first.unit().mag(), 1.0, kTolerance);
}

TEST(QuatTest, RotatesVectorAroundZAxis) {
  Quat rotation{std::cos(HalfPi() / 2.0), 0.0, 0.0,
                std::sin(HalfPi() / 2.0)};
  Vector vector{1.0, 0.0, 0.0};

  rotation.rotate(vector);

  ExpectVectorNear(vector, 0.0, 1.0, 0.0);
}

TEST(MatrixMathTest, RotatesVectorWithEulerMatrix) {
  std::array<double, 9> matrix =
      create_euler_rotation_matrix(0.0, 0.0, HalfPi());
  Vector vector{1.0, 0.0, 0.0};

  Vector rotated = matrix_rotate(vector, matrix);

  ExpectVectorNear(rotated, 0.0, 1.0, 0.0);
}

TEST(MatrixMathTest, CoordAndScalarEulerMatrixOverloadsMatch) {
  Coord angles{0.25, -0.5, 0.75};

  std::array<double, 9> from_coord = create_euler_rotation_matrix(angles);
  std::array<double, 9> from_scalars =
      create_euler_rotation_matrix(angles.x, angles.y, angles.z);

  for (std::size_t index = 0; index < from_coord.size(); ++index) {
    EXPECT_NEAR(from_coord[index], from_scalars[index], kTolerance);
  }
}

}  // namespace
