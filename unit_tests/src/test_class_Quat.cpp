// test_class_Quat.cpp
// Unit tests for class_Quat.cpp
// Tests the Quat quaternion class functions: operator*, operator<<, norm,
// mag, scale, conjugate, inverse, unit, and rotate.

#include "classes/class_Quat.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <sstream>

// Small tolerance for floating point comparisons.
static const double TOL = 1e-9;

// Test the quaternion multiplication operator*.
// Uses known Hamilton product identities: i*j = k, etc.
void test_Quat_multiply()
{
    std::cerr << "[class_Quat.cpp] Testing Quat::operator* (quaternion multiplication)\n";

    // i * j should equal k.
    // i = (0,1,0,0), j = (0,0,1,0), expected k = (0,0,0,1)
    Quat i { 0, 1, 0, 0 };
    Quat j { 0, 0, 1, 0 };
    Quat result = i * j;
    std::cerr << "  i * j = " << result << " (expected [0, 0i, 0j, 1k])\n";
    EXPECT_NEAR(result.w, 0.0, TOL);
    EXPECT_NEAR(result.x, 0.0, TOL);
    EXPECT_NEAR(result.y, 0.0, TOL);
    EXPECT_NEAR(result.z, 1.0, TOL);

    // Multiplying by the identity quaternion (1,0,0,0) should return the same quaternion.
    Quat identity { 1, 0, 0, 0 };
    Quat q { 2, 3, 4, 5 };
    Quat qId = q * identity;
    std::cerr << "  q * identity = " << qId << " (expected [2, 3i, 4j, 5k])\n";
    EXPECT_NEAR(qId.w, 2.0, TOL);
    EXPECT_NEAR(qId.x, 3.0, TOL);
    EXPECT_NEAR(qId.y, 4.0, TOL);
    EXPECT_NEAR(qId.z, 5.0, TOL);
}

// Test the stream insertion operator<<.
// Verifies that a Quat is formatted as "[w, xi, yj, zk]".
void test_Quat_ostream()
{
    std::cerr << "[class_Quat.cpp] Testing operator<< (stream output)\n";

    Quat q { 1, 2, 3, 4 };
    std::ostringstream oss;
    oss << q;
    std::string expected = "[1, 2i, 3j, 4k]";
    std::cerr << "  Output: \"" << oss.str() << "\" (expected \"" << expected << "\")\n";
    EXPECT_EQ(oss.str(), expected);
}

// Test the norm() function which returns the squared magnitude.
void test_Quat_norm()
{
    std::cerr << "[class_Quat.cpp] Testing Quat::norm (squared magnitude)\n";

    // norm = w^2 + x^2 + y^2 + z^2 = 1+4+9+16 = 30
    Quat q { 1, 2, 3, 4 };
    double n = q.norm();
    std::cerr << "  norm of [1,2,3,4] = " << n << " (expected 30)\n";
    EXPECT_NEAR(n, 30.0, TOL);
}

// Test the mag() function which returns the magnitude (sqrt of norm).
void test_Quat_mag()
{
    std::cerr << "[class_Quat.cpp] Testing Quat::mag (magnitude)\n";

    // mag = sqrt(30)
    Quat q { 1, 2, 3, 4 };
    double m = q.mag();
    std::cerr << "  mag of [1,2,3,4] = " << m << " (expected " << std::sqrt(30.0) << ")\n";
    EXPECT_NEAR(m, std::sqrt(30.0), TOL);
}

// Test the scale() function which multiplies all components by a scalar.
void test_Quat_scale()
{
    std::cerr << "[class_Quat.cpp] Testing Quat::scale (scalar multiplication)\n";

    Quat q { 1, 2, 3, 4 };
    Quat scaled = q.scale(2.0);
    std::cerr << "  [1,2,3,4] scaled by 2 = " << scaled << " (expected [2, 4i, 6j, 8k])\n";
    EXPECT_NEAR(scaled.w, 2.0, TOL);
    EXPECT_NEAR(scaled.x, 4.0, TOL);
    EXPECT_NEAR(scaled.y, 6.0, TOL);
    EXPECT_NEAR(scaled.z, 8.0, TOL);
}

// Test the conjugate() function which negates the vector components.
void test_Quat_conjugate()
{
    std::cerr << "[class_Quat.cpp] Testing Quat::conjugate\n";

    Quat q { 1, 2, 3, 4 };
    Quat conj = q.conjugate();
    std::cerr << "  conjugate of [1,2,3,4] = " << conj << " (expected [1, -2i, -3j, -4k])\n";
    EXPECT_NEAR(conj.w, 1.0, TOL);
    EXPECT_NEAR(conj.x, -2.0, TOL);
    EXPECT_NEAR(conj.y, -3.0, TOL);
    EXPECT_NEAR(conj.z, -4.0, TOL);
}

// Test the inverse() function which is conjugate scaled by 1/norm.
void test_Quat_inverse()
{
    std::cerr << "[class_Quat.cpp] Testing Quat::inverse\n";

    // For q=[1,2,3,4], norm=30, conjugate=[1,-2,-3,-4]
    // inverse = [1/30, -2/30, -3/30, -4/30]
    Quat q { 1, 2, 3, 4 };
    Quat inv = q.inverse();
    std::cerr << "  inverse of [1,2,3,4] = " << inv << "\n";
    EXPECT_NEAR(inv.w, 1.0 / 30.0, TOL);
    EXPECT_NEAR(inv.x, -2.0 / 30.0, TOL);
    EXPECT_NEAR(inv.y, -3.0 / 30.0, TOL);
    EXPECT_NEAR(inv.z, -4.0 / 30.0, TOL);

    // Additionally, q * q.inverse() should give the identity quaternion [1,0,0,0].
    Quat prod = q * q.inverse();
    std::cerr << "  q * q.inverse() = " << prod << " (expected [1, 0i, 0j, 0k])\n";
    EXPECT_NEAR(prod.w, 1.0, TOL);
    EXPECT_NEAR(prod.x, 0.0, TOL);
    EXPECT_NEAR(prod.y, 0.0, TOL);
    EXPECT_NEAR(prod.z, 0.0, TOL);
}

// Test the unit() function which normalizes the quaternion to magnitude 1.
void test_Quat_unit()
{
    std::cerr << "[class_Quat.cpp] Testing Quat::unit (normalization)\n";

    Quat q { 1, 2, 3, 4 };
    Quat u = q.unit();
    // The magnitude of the resulting unit quaternion should be 1.
    double mag = u.mag();
    std::cerr << "  unit of [1,2,3,4] = " << u << ", magnitude = " << mag << " (expected 1)\n";
    EXPECT_NEAR(mag, 1.0, TOL);

    // Verify individual components equal original / sqrt(30).
    double d = std::sqrt(30.0);
    EXPECT_NEAR(u.w, 1.0 / d, TOL);
    EXPECT_NEAR(u.x, 2.0 / d, TOL);
    EXPECT_NEAR(u.y, 3.0 / d, TOL);
    EXPECT_NEAR(u.z, 4.0 / d, TOL);
}

// Test the rotate() function which rotates a Vector using this quaternion.
void test_Quat_rotate()
{
    std::cerr << "[class_Quat.cpp] Testing Quat::rotate (vector rotation)\n";

    // A 90-degree rotation about the z-axis.
    // Unit quaternion for angle theta about z: [cos(theta/2), 0, 0, sin(theta/2)]
    double theta = M_PI / 2.0;
    Quat rot { std::cos(theta / 2.0), 0.0, 0.0, std::sin(theta / 2.0) };

    // Rotating the x unit vector (1,0,0) by 90 deg about z should give (0,1,0).
    Vector vec;
    vec.x = 1.0;
    vec.y = 0.0;
    vec.z = 0.0;

    rot.rotate(vec);
    std::cerr << "  Rotated (1,0,0) by 90deg about z -> ("
              << vec.x << ", " << vec.y << ", " << vec.z << ") (expected (0,1,0))\n";
    EXPECT_NEAR(vec.x, 0.0, TOL);
    EXPECT_NEAR(vec.y, 1.0, TOL);
    EXPECT_NEAR(vec.z, 0.0, TOL);

    // Rotation about z should not change the z component.
    Vector vec2;
    vec2.x = 0.0;
    vec2.y = 0.0;
    vec2.z = 5.0;

    rot.rotate(vec2);
    std::cerr << "  Rotated (0,0,5) by 90deg about z -> ("
              << vec2.x << ", " << vec2.y << ", " << vec2.z << ") (expected (0,0,5))\n";
    EXPECT_NEAR(vec2.x, 0.0, TOL);
    EXPECT_NEAR(vec2.y, 0.0, TOL);
    EXPECT_NEAR(vec2.z, 5.0, TOL);
}

// GoogleTest wrappers that invoke each test function.
// Grouping under a single test suite name "QuatTest".

TEST(QuatTest, Multiply)   { test_Quat_multiply(); }
TEST(QuatTest, Ostream)    { test_Quat_ostream(); }
TEST(QuatTest, Norm)       { test_Quat_norm(); }
TEST(QuatTest, Mag)        { test_Quat_mag(); }
TEST(QuatTest, Scale)      { test_Quat_scale(); }
TEST(QuatTest, Conjugate)  { test_Quat_conjugate(); }
TEST(QuatTest, Inverse)    { test_Quat_inverse(); }
TEST(QuatTest, Unit)       { test_Quat_unit(); }
TEST(QuatTest, Rotate)     { test_Quat_rotate(); }