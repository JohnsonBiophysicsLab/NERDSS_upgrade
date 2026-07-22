/*! \file test_class_Coord.cpp
 *
 * Unit tests for class_Coord.cpp
 *
 * Tests all constructors, operators, free functions, and member functions
 * of the Coord class using the GoogleTest framework.
 */

#include "classes/class_Coord.hpp"
#include <classes/class_Parameters.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

// -----------------------------------------------------------------------------
// Test the three Coord constructors.
// -----------------------------------------------------------------------------
void coord_test_constructors()
{
    std::cerr << "[class_Coord.cpp] Testing Coord constructors...\n";

    // Test the (double, double, double) constructor.
    std::cerr << "  -> Testing Coord(double, double, double)\n";
    Coord c1(1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(c1.x, 1.0) << "x should be 1.0";
    EXPECT_DOUBLE_EQ(c1.y, 2.0) << "y should be 2.0";
    EXPECT_DOUBLE_EQ(c1.z, 3.0) << "z should be 3.0";

    // Test the std::array<double,3>& constructor.
    std::cerr << "  -> Testing Coord(std::array<double,3>&)\n";
    std::array<double, 3> arr { 4.0, 5.0, 6.0 };
    Coord c2(arr);
    EXPECT_DOUBLE_EQ(c2.x, 4.0) << "x should be 4.0";
    EXPECT_DOUBLE_EQ(c2.y, 5.0) << "y should be 5.0";
    EXPECT_DOUBLE_EQ(c2.z, 6.0) << "z should be 6.0";

    // Test the std::vector<double> constructor (valid size).
    std::cerr << "  -> Testing Coord(std::vector<double>)\n";
    std::vector<double> vals { 7.0, 8.0, 9.0 };
    Coord c3(vals);
    EXPECT_DOUBLE_EQ(c3.x, 7.0) << "x should be 7.0";
    EXPECT_DOUBLE_EQ(c3.y, 8.0) << "y should be 8.0";
    EXPECT_DOUBLE_EQ(c3.z, 9.0) << "z should be 9.0";
}

// -----------------------------------------------------------------------------
// Test the free function roundv() which rounds to 4 decimal places.
// -----------------------------------------------------------------------------
void test_roundv()
{
    std::cerr << "[class_Coord.cpp] Testing roundv()...\n";

    // Positive value rounding.
    std::cerr << "  -> Testing roundv() with positive value 1.234567\n";
    EXPECT_DOUBLE_EQ(roundv(1.234567), 1.2346) << "1.234567 should round to 1.2346";

    // Negative value rounding.
    std::cerr << "  -> Testing roundv() with negative value -1.234567\n";
    EXPECT_DOUBLE_EQ(roundv(-1.234567), -1.2346) << "-1.234567 should round to -1.2346";

    // Zero rounding.
    std::cerr << "  -> Testing roundv() with 0.0\n";
    EXPECT_DOUBLE_EQ(roundv(0.0), 0.0) << "0.0 should round to 0.0";
}

// -----------------------------------------------------------------------------
// Test the free function round(Coord) which rounds each component.
// -----------------------------------------------------------------------------
void test_round()
{
    std::cerr << "[class_Coord.cpp] Testing round(const Coord&)...\n";

    Coord c(1.234567, -2.345678, 3.456789);
    Coord r = round(c);
    std::cerr << "  -> Testing round on Coord(1.234567, -2.345678, 3.456789)\n";
    EXPECT_DOUBLE_EQ(r.x, 1.2346) << "x should round to 1.2346";
    EXPECT_DOUBLE_EQ(r.y, -2.3457) << "y should round to -2.3457";
    EXPECT_DOUBLE_EQ(r.z, 3.4568) << "z should round to 3.4568";
}

// -----------------------------------------------------------------------------
// Test the is_co_linear() function using Heron's formula.
// -----------------------------------------------------------------------------
void test_is_co_linear()
{
    std::cerr << "[class_Coord.cpp] Testing is_co_linear()...\n";

    // Three points on a straight line along the x-axis are co-linear.
    std::cerr << "  -> Testing three co-linear points along x-axis\n";
    Coord a(0.0, 0.0, 0.0);
    Coord b(1.0, 0.0, 0.0);
    Coord c(2.0, 0.0, 0.0);
    EXPECT_TRUE(is_co_linear(a, b, c)) << "Points (0,0,0), (1,0,0), (2,0,0) should be co-linear";

    // Three points forming a triangle are not co-linear.
    std::cerr << "  -> Testing three non-co-linear points (triangle)\n";
    Coord d(0.0, 0.0, 0.0);
    Coord e(1.0, 0.0, 0.0);
    Coord f(0.0, 1.0, 0.0);
    EXPECT_FALSE(is_co_linear(d, e, f)) << "Points forming a triangle should not be co-linear";
}

// -----------------------------------------------------------------------------
// Test the comparison operators == and !=.
// -----------------------------------------------------------------------------
void test_comparison_operators()
{
    std::cerr << "[class_Coord.cpp] Testing operator== and operator!=...\n";

    // Two identical coordinates should compare equal.
    std::cerr << "  -> Testing operator== with equal coordinates\n";
    Coord c1(1.0, 2.0, 3.0);
    Coord c2(1.0, 2.0, 3.0);
    EXPECT_TRUE(c1 == c2) << "Identical coordinates should be equal";
    EXPECT_FALSE(c1 != c2) << "Identical coordinates should not be unequal";

    // Two different coordinates should not compare equal.
    std::cerr << "  -> Testing operator!= with different coordinates\n";
    Coord c3(1.0, 2.0, 4.0);
    EXPECT_FALSE(c1 == c3) << "Different coordinates should not be equal";
    EXPECT_TRUE(c1 != c3) << "Different coordinates should be unequal";
}

// -----------------------------------------------------------------------------
// Test the stream insertion operator <<.
// -----------------------------------------------------------------------------
void test_ostream_operator()
{
    std::cerr << "[class_Coord.cpp] Testing operator<<...\n";

    Coord c(1.5, 2.5, 3.5);
    std::ostringstream oss;
    oss << c;
    std::cerr << "  -> Streaming Coord(1.5, 2.5, 3.5) produced: \"" << oss.str() << "\"\n";
    // The output should be non-empty and contain the coordinate values.
    EXPECT_FALSE(oss.str().empty()) << "Streamed output should not be empty";
    EXPECT_NE(oss.str().find("1.5"), std::string::npos) << "Output should contain 1.5";
    EXPECT_NE(oss.str().find("2.5"), std::string::npos) << "Output should contain 2.5";
    EXPECT_NE(oss.str().find("3.5"), std::string::npos) << "Output should contain 3.5";
}

// -----------------------------------------------------------------------------
// Test the addition operators.
// -----------------------------------------------------------------------------
void test_addition_operators()
{
    std::cerr << "[class_Coord.cpp] Testing addition operators...\n";

    // Test Coord + Coord.
    std::cerr << "  -> Testing operator+(Coord, Coord)\n";
    Coord c1(1.0, 2.0, 3.0);
    Coord c2(4.0, 5.0, 6.0);
    Coord sum = c1 + c2;
    EXPECT_DOUBLE_EQ(sum.x, 5.0) << "x sum should be 5.0";
    EXPECT_DOUBLE_EQ(sum.y, 7.0) << "y sum should be 7.0";
    EXPECT_DOUBLE_EQ(sum.z, 9.0) << "z sum should be 9.0";

    // Test std::array + Coord.
    std::cerr << "  -> Testing operator+(std::array<double,3>, Coord)\n";
    std::array<double, 3> arr { 10.0, 20.0, 30.0 };
    Coord sum2 = arr + c1;
    EXPECT_DOUBLE_EQ(sum2.x, 11.0) << "x sum should be 11.0";
    EXPECT_DOUBLE_EQ(sum2.y, 22.0) << "y sum should be 22.0";
    EXPECT_DOUBLE_EQ(sum2.z, 33.0) << "z sum should be 33.0";
}

// -----------------------------------------------------------------------------
// Test the subtraction operator (Coord - double).
// -----------------------------------------------------------------------------
void test_subtraction_operator()
{
    std::cerr << "[class_Coord.cpp] Testing operator-(Coord, double)...\n";

    Coord c(10.0, 20.0, 30.0);
    Coord diff = c - 5.0;
    std::cerr << "  -> Subtracting 5.0 from Coord(10,20,30)\n";
    EXPECT_DOUBLE_EQ(diff.x, 5.0) << "x should be 5.0";
    EXPECT_DOUBLE_EQ(diff.y, 15.0) << "y should be 15.0";
    EXPECT_DOUBLE_EQ(diff.z, 25.0) << "z should be 25.0";
}

// -----------------------------------------------------------------------------
// Test the compound assignment operators += and /=.
// -----------------------------------------------------------------------------
void test_compound_assignment_operators()
{
    std::cerr << "[class_Coord.cpp] Testing compound assignment operators...\n";

    // Test Coord::operator+=(Coord).
    std::cerr << "  -> Testing Coord::operator+=(Coord)\n";
    Coord c1(1.0, 2.0, 3.0);
    Coord c2(4.0, 5.0, 6.0);
    c1 += c2;
    EXPECT_DOUBLE_EQ(c1.x, 5.0) << "x should be 5.0 after +=";
    EXPECT_DOUBLE_EQ(c1.y, 7.0) << "y should be 7.0 after +=";
    EXPECT_DOUBLE_EQ(c1.z, 9.0) << "z should be 9.0 after +=";

    // Test operator/=(Coord, double).
    std::cerr << "  -> Testing operator/=(Coord, double)\n";
    Coord c3(10.0, 20.0, 30.0);
    double scal = 2.0;
    c3 /= scal;
    EXPECT_DOUBLE_EQ(c3.x, 5.0) << "x should be 5.0 after /=";
    EXPECT_DOUBLE_EQ(c3.y, 10.0) << "y should be 10.0 after /=";
    EXPECT_DOUBLE_EQ(c3.z, 15.0) << "z should be 15.0 after /=";

    // Test operator+=(Coord, std::array) which returns a new Coord.
    std::cerr << "  -> Testing operator+=(Coord, std::array) (returns new Coord)\n";
    Coord c4(1.0, 1.0, 1.0);
    std::array<double, 3> arr { 2.0, 3.0, 4.0 };
    Coord result = (c4 += arr);
    EXPECT_DOUBLE_EQ(result.x, 3.0) << "x should be 3.0";
    EXPECT_DOUBLE_EQ(result.y, 4.0) << "y should be 4.0";
    EXPECT_DOUBLE_EQ(result.z, 5.0) << "z should be 5.0";
}

// -----------------------------------------------------------------------------
// Test the member function zero_crds().
// -----------------------------------------------------------------------------
void test_zero_crds()
{
    std::cerr << "[class_Coord.cpp] Testing Coord::zero_crds()...\n";

    Coord c(1.0, 2.0, 3.0);
    c.zero_crds();
    std::cerr << "  -> After zero_crds() all components should be 0\n";
    EXPECT_DOUBLE_EQ(c.x, 0.0) << "x should be 0.0";
    EXPECT_DOUBLE_EQ(c.y, 0.0) << "y should be 0.0";
    EXPECT_DOUBLE_EQ(c.z, 0.0) << "z should be 0.0";
}

// -----------------------------------------------------------------------------
// Test the member function get_magnitude().
// -----------------------------------------------------------------------------
void test_get_magnitude()
{
    std::cerr << "[class_Coord.cpp] Testing Coord::get_magnitude()...\n";

    // Magnitude of (3,4,0) should be 5.
    std::cerr << "  -> Testing get_magnitude on Coord(3,4,0)\n";
    Coord c(3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(c.get_magnitude(), 5.0) << "Magnitude of (3,4,0) should be 5.0";

    // Magnitude of (0,0,0) should be 0.
    std::cerr << "  -> Testing get_magnitude on Coord(0,0,0)\n";
    Coord zero(0.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(zero.get_magnitude(), 0.0) << "Magnitude of (0,0,0) should be 0.0";
}

// -----------------------------------------------------------------------------
// Test the member function isOutOfBox().
// -----------------------------------------------------------------------------
void test_isOutOfBox()
{
    std::cerr << "[class_Coord.cpp] Testing Coord::isOutOfBox()...\n";

    // Set up a Membrane object with a water box.
    Membrane membraneObject;
    membraneObject.waterBox.x = 10.0;
    membraneObject.waterBox.y = 10.0;
    membraneObject.waterBox.z = 10.0;

    // A coordinate inside the box should return false.
    std::cerr << "  -> Testing coordinate inside the water box\n";
    Coord inside(1.0, 1.0, 1.0);
    EXPECT_FALSE(inside.isOutOfBox(membraneObject)) << "Coordinate inside box should not be out of box";

    // A coordinate outside the box in x should return true.
    std::cerr << "  -> Testing coordinate outside the water box (x too large)\n";
    Coord outsideX(6.0, 1.0, 1.0);
    EXPECT_TRUE(outsideX.isOutOfBox(membraneObject)) << "Coordinate with x=6 should be out of box";

    // A coordinate outside the box in y (negative) should return true.
    std::cerr << "  -> Testing coordinate outside the water box (y too small)\n";
    Coord outsideY(1.0, -6.0, 1.0);
    EXPECT_TRUE(outsideY.isOutOfBox(membraneObject)) << "Coordinate with y=-6 should be out of box";

    // A coordinate outside the box in z should return true.
    std::cerr << "  -> Testing coordinate outside the water box (z too large)\n";
    Coord outsideZ(1.0, 1.0, 6.0);
    EXPECT_TRUE(outsideZ.isOutOfBox(membraneObject)) << "Coordinate with z=6 should be out of box";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers that invoke each test function.
// -----------------------------------------------------------------------------
TEST(ClassCoordTest, Constructors) { coord_test_constructors(); }
TEST(ClassCoordTest, Roundv) { test_roundv(); }
TEST(ClassCoordTest, Round) { test_round(); }
TEST(ClassCoordTest, IsCoLinear) { test_is_co_linear(); }
TEST(ClassCoordTest, ComparisonOperators) { test_comparison_operators(); }
TEST(ClassCoordTest, OstreamOperator) { test_ostream_operator(); }
TEST(ClassCoordTest, AdditionOperators) { test_addition_operators(); }
TEST(ClassCoordTest, SubtractionOperator) { test_subtraction_operator(); }
TEST(ClassCoordTest, CompoundAssignmentOperators) { test_compound_assignment_operators(); }
TEST(ClassCoordTest, ZeroCrds) { test_zero_crds(); }
TEST(ClassCoordTest, GetMagnitude) { test_get_magnitude(); }
TEST(ClassCoordTest, IsOutOfBox) { test_isOutOfBox(); }