#include "classes/class_Coord.hpp"
#include "classes/class_Vector.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

constexpr double kTolerance = 1.0e-12;

void require_close(double actual, double expected, const std::string& label)
{
    if (std::abs(actual - expected) > kTolerance) {
        std::cerr << label << ": expected " << expected << ", got " << actual << '\n';
        std::exit(1);
    }
}

void require_true(bool condition, const std::string& label)
{
    if (!condition) {
        std::cerr << label << '\n';
        std::exit(1);
    }
}

void test_coord_rounding_and_colinearity()
{
    Coord rounded = round({ 1.23456, -1.23456, 0.00004 });
    require_close(rounded.x, 1.2346, "round positive coordinate");
    require_close(rounded.y, -1.2346, "round negative coordinate");
    require_close(rounded.z, 0.0, "round near-zero coordinate");

    Coord first { 0.0, 0.0, 0.0 };
    Coord second { 1.0, 1.0, 1.0 };
    Coord third { 2.0, 2.0, 2.0 };
    require_true(is_co_linear(first, second, third), "expected diagonal points to be co-linear");
}

void test_vector_magnitude_dot_and_normalize()
{
    Vector vector { 3.0, 4.0, 12.0 };
    vector.calc_magnitude();
    require_close(vector.magnitude, 13.0, "vector magnitude");

    Vector x_axis { 1.0, 0.0, 0.0 };
    Vector y_axis { 0.0, 1.0, 0.0 };
    require_close(x_axis.dot(y_axis), 0.0, "orthogonal dot product");

    vector.normalize();
    require_close(vector.magnitude, 1.0, "normalized magnitude");
    require_close(vector.x, 3.0 / 13.0, "normalized x");
    require_close(vector.y, 4.0 / 13.0, "normalized y");
    require_close(vector.z, 12.0 / 13.0, "normalized z");
}

void test_vector_cross_projection_and_angle()
{
    Vector x_axis { 1.0, 0.0, 0.0 };
    Vector y_axis { 0.0, 1.0, 0.0 };
    Vector cross = x_axis.cross(y_axis);
    require_close(cross.x, 0.0, "cross product x");
    require_close(cross.y, 0.0, "cross product y");
    require_close(cross.z, 1.0, "cross product z");
    require_close(cross.magnitude, 1.0, "cross product magnitude");

    Vector original { 2.0, 3.0, 4.0 };
    Vector normal { 0.0, 0.0, 1.0 };
    Vector projected = original.vector_projection(normal);
    require_close(projected.x, 2.0, "projection x");
    require_close(projected.y, 3.0, "projection y");
    require_close(projected.z, 0.0, "projection removes normal component");

    x_axis.calc_magnitude();
    y_axis.calc_magnitude();
    require_close(x_axis.dot_theta(y_axis), std::acos(0.0), "right angle between axes");
}

} // namespace

int main()
{
    test_coord_rounding_and_colinearity();
    test_vector_magnitude_dot_and_normalize();
    test_vector_cross_projection_and_angle();
    return 0;
}
