
// test_class_Vector.cpp
// Unit tests for the Vector class defined in class_Vector.cpp
// These tests exercise the constructors, operators, and member functions
// of the Vector class, printing verbose information about each test.

#include "classes/class_Vector.hpp" // The header for the Vector (and Coord) classes

#include <gtest/gtest.h> // Google Test framework
#include <cmath>         // For std::sqrt, std::acos, M_PI, std::abs
#include <array>         // For std::array
#include <vector>        // For std::vector
#include <sstream>       // For std::ostringstream (to test operator<<)

// A small tolerance used for floating-point comparisons throughout the tests.
static const double kTol = 1e-9;

/*
 * test_constructors()
 * -------------------
 * Verifies that the various Vector constructors correctly initialize the
 * x, y, and z members.
 */
void test_constructors()
{
    std::cerr << "\n[TEST] test_constructors in class_Vector.cpp\n";

    // --- 2-argument constructor Vector(double x, double y) ---
    std::cerr << "  Testing Vector(double x, double y) constructor...\n";
    Vector v2(3.0, 4.0);
    // Only x and y are set explicitly by this constructor.
    EXPECT_NEAR(v2.x, 3.0, kTol) << "Vector(x,y): x should be 3.0";
    EXPECT_NEAR(v2.y, 4.0, kTol) << "Vector(x,y): y should be 4.0";

    // --- 3-argument constructor Vector(const double&, const double&, const double&) ---
    std::cerr << "  Testing Vector(x, y, z) constructor...\n";
    Vector v3(1.0, 2.0, 3.0);
    EXPECT_NEAR(v3.x, 1.0, kTol) << "Vector(x,y,z): x should be 1.0";
    EXPECT_NEAR(v3.y, 2.0, kTol) << "Vector(x,y,z): y should be 2.0";
    EXPECT_NEAR(v3.z, 3.0, kTol) << "Vector(x,y,z): z should be 3.0";

    // --- Constructor from two Coords: Vector(end, start) ---
    std::cerr << "  Testing Vector(Coord end, Coord start) constructor...\n";
    Coord start; // default-constructed origin (assumed 0,0,0)
    start.x = 1.0; start.y = 1.0; start.z = 1.0;
    Coord end;
    end.x = 4.0; end.y = 5.0; end.z = 6.0;
    Vector vEndStart(end, start);
    // Resulting vector components should be the difference end - start.
    EXPECT_NEAR(vEndStart.x, 3.0, kTol) << "Vector(end,start): x = end.x - start.x";
    EXPECT_NEAR(vEndStart.y, 4.0, kTol) << "Vector(end,start): y = end.y - start.y";
    EXPECT_NEAR(vEndStart.z, 5.0, kTol) << "Vector(end,start): z = end.z - start.z";

    // --- Constructor from a single Coord (implied origin start) ---
    std::cerr << "  Testing Vector(Coord) constructor...\n";
    Coord single;
    single.x = 7.0; single.y = 8.0; single.z = 9.0;
    Vector vSingle(single);
    EXPECT_NEAR(vSingle.x, 7.0, kTol) << "Vector(Coord): x should copy coord.x";
    EXPECT_NEAR(vSingle.y, 8.0, kTol) << "Vector(Coord): y should copy coord.y";
    EXPECT_NEAR(vSingle.z, 9.0, kTol) << "Vector(Coord): z should copy coord.z";

    // --- Constructor from std::array<double,3> ---
    std::cerr << "  Testing Vector(std::array<double,3>&) constructor...\n";
    std::array<double, 3> arr = { 10.0, 11.0, 12.0 };
    Vector vArr(arr);
    EXPECT_NEAR(vArr.x, 10.0, kTol) << "Vector(array): x should be arr[0]";
    EXPECT_NEAR(vArr.y, 11.0, kTol) << "Vector(array): y should be arr[1]";
    EXPECT_NEAR(vArr.z, 12.0, kTol) << "Vector(array): z should be arr[2]";

    // --- Constructor from std::vector<double> (size == 3) ---
    std::cerr << "  Testing Vector(std::vector<double>) constructor...\n";
    std::vector<double> vec = { 13.0, 14.0, 15.0 };
    Vector vVec(vec);
    EXPECT_NEAR(vVec.x, 13.0, kTol) << "Vector(vector): x should be vec[0]";
    EXPECT_NEAR(vVec.y, 14.0, kTol) << "Vector(vector): y should be vec[1]";
    EXPECT_NEAR(vVec.z, 15.0, kTol) << "Vector(vector): z should be vec[2]";
}

/*
 * test_output_operator()
 * ----------------------
 * Verifies the overloaded operator<< produces the expected formatted string
 * "[xi + yj + zk]".
 */
void test_output_operator()
{
    std::cerr << "\n[TEST] test_output_operator (operator<<) in class_Vector.cpp\n";

    Vector v(1.0, 2.0, 3.0);
    std::ostringstream oss;
    oss << v; // Uses the overloaded stream insertion operator.

    std::string expected = "[1i + 2j + 3k]";
    std::cerr << "  Produced string: " << oss.str() << "\n";
    std::cerr << "  Expected string: " << expected << "\n";
    EXPECT_EQ(oss.str(), expected) << "operator<< should format as [xi + yj + zk]";
}

/*
 * test_plus_operator()
 * --------------------
 * Verifies the operator+ between a Vector and a Coord returns a Coord that is
 * the component-wise sum.
 */
void test_plus_operator()
{
    std::cerr << "\n[TEST] test_plus_operator (Vector + Coord) in class_Vector.cpp\n";

    Vector v(1.0, 2.0, 3.0);
    Coord c;
    c.x = 4.0; c.y = 5.0; c.z = 6.0;

    Coord result = v + c; // Component-wise addition.
    std::cerr << "  (1,2,3) + (4,5,6) = (" << result.x << "," << result.y << "," << result.z << ")\n";
    EXPECT_NEAR(result.x, 5.0, kTol) << "Vector+Coord: x should be 5.0";
    EXPECT_NEAR(result.y, 7.0, kTol) << "Vector+Coord: y should be 7.0";
    EXPECT_NEAR(result.z, 9.0, kTol) << "Vector+Coord: z should be 9.0";
}

/*
 * test_cross()
 * ------------
 * Verifies the cross product. Note the implementation normalizes the result,
 * so we compare against the normalized cross product direction.
 */
void test_cross()
{
    std::cerr << "\n[TEST] test_cross in class_Vector.cpp\n";

    // x-axis cross y-axis = z-axis (already unit length, so normalization is a no-op).
    Vector a(1.0, 0.0, 0.0);
    Vector b(0.0, 1.0, 0.0);
    Vector result = a.cross(b);
    std::cerr << "  (1,0,0) x (0,1,0) normalized = ("
              << result.x << "," << result.y << "," << result.z << ")\n";
    EXPECT_NEAR(result.x, 0.0, kTol) << "cross: x component should be 0";
    EXPECT_NEAR(result.y, 0.0, kTol) << "cross: y component should be 0";
    EXPECT_NEAR(result.z, 1.0, kTol) << "cross: z component should be 1 (unit z)";

    // A non-unit example: (1,2,3) x (4,5,6) = (-3,6,-3), which normalizes.
    Vector c(1.0, 2.0, 3.0);
    Vector d(4.0, 5.0, 6.0);
    Vector result2 = c.cross(d);
    // Unnormalized cross is (-3, 6, -3), magnitude = sqrt(9+36+9) = sqrt(54).
    double mag = std::sqrt(54.0);
    std::cerr << "  (1,2,3) x (4,5,6) normalized = ("
              << result2.x << "," << result2.y << "," << result2.z << ")\n";
    EXPECT_NEAR(result2.x, -3.0 / mag, kTol) << "cross: normalized x component";
    EXPECT_NEAR(result2.y, 6.0 / mag, kTol) << "cross: normalized y component";
    EXPECT_NEAR(result2.z, -3.0 / mag, kTol) << "cross: normalized z component";
}

/*
 * test_dot()
 * ----------
 * Verifies the dot product returns the correct scalar value.
 */
void test_dot()
{
    std::cerr << "\n[TEST] test_dot in class_Vector.cpp\n";

    Vector a(1.0, 2.0, 3.0);
    Vector b(4.0, 5.0, 6.0);
    // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    double dp = a.dot(b);
    std::cerr << "  (1,2,3) . (4,5,6) = " << dp << " (expected 32)\n";
    EXPECT_NEAR(dp, 32.0, kTol) << "dot: should be 32";

    // Orthogonal vectors should have a dot product of zero.
    Vector x(1.0, 0.0, 0.0);
    Vector y(0.0, 1.0, 0.0);
    double dpOrtho = x.dot(y);
    std::cerr << "  orthogonal (1,0,0).(0,1,0) = " << dpOrtho << " (expected 0)\n";
    EXPECT_NEAR(dpOrtho, 0.0, kTol) << "dot: orthogonal vectors give 0";
}

/*
 * test_dot_theta()
 * ----------------
 * Verifies the angle between vectors is computed correctly. This requires
 * magnitudes to be set, so calc_magnitude() is called first.
 */
void test_dot_theta()
{
    std::cerr << "\n[TEST] test_dot_theta in class_Vector.cpp\n";

    // 90 degree angle between x-axis and y-axis => pi/2.
    Vector a(1.0, 0.0, 0.0);
    Vector b(0.0, 1.0, 0.0);
    a.calc_magnitude(); // dot_theta relies on magnitude members being populated.
    b.calc_magnitude();
    double angle = a.dot_theta(b);
    std::cerr << "  angle between (1,0,0) and (0,1,0) = " << angle
              << " (expected " << M_PI / 2 << ")\n";
    EXPECT_NEAR(angle, M_PI / 2, 1e-6) << "dot_theta: perpendicular vectors -> pi/2";

    // Parallel vectors => angle 0.
    Vector c(1.0, 0.0, 0.0);
    Vector d(2.0, 0.0, 0.0);
    c.calc_magnitude();
    d.calc_magnitude();
    double angleParallel = c.dot_theta(d);
    std::cerr << "  angle between (1,0,0) and (2,0,0) = " << angleParallel
              << " (expected 0)\n";
    EXPECT_NEAR(angleParallel, 0.0, 1e-6) << "dot_theta: parallel vectors -> 0";

    // Anti-parallel vectors => angle pi.
    Vector e(1.0, 0.0, 0.0);
    Vector f(-1.0, 0.0, 0.0);
    e.calc_magnitude();
    f.calc_magnitude();
    double angleAnti = e.dot_theta(f);
    std::cerr << "  angle between (1,0,0) and (-1,0,0) = " << angleAnti
              << " (expected " << M_PI << ")\n";
    EXPECT_NEAR(angleAnti, M_PI, 1e-6) << "dot_theta: anti-parallel vectors -> pi";

    // Zero-magnitude vector should return 0.0 (with a printed warning).
    Vector zero(0.0, 0.0, 0.0);
    Vector g(1.0, 0.0, 0.0);
    zero.calc_magnitude();
    g.calc_magnitude();
    double angleZero = zero.dot_theta(g);
    std::cerr << "  angle with a zero-magnitude vector = " << angleZero
              << " (expected 0)\n";
    EXPECT_NEAR(angleZero, 0.0, kTol) << "dot_theta: zero magnitude -> returns 0";
}

/*
 * test_calc_magnitude()
 * ---------------------
 * Verifies that calc_magnitude() correctly computes the Euclidean norm.
 */
void test_calc_magnitude()
{
    std::cerr << "\n[TEST] test_calc_magnitude in class_Vector.cpp\n";

    // 3-4-5 style in 3D: sqrt(3^2 + 4^2 + 0^2) = 5.
    Vector v(3.0, 4.0, 0.0);
    v.calc_magnitude();
    std::cerr << "  magnitude of (3,4,0) = " << v.magnitude << " (expected 5)\n";
    EXPECT_NEAR(v.magnitude, 5.0, kTol) << "calc_magnitude: should be 5";

    // Full 3D: sqrt(1+4+4) = 3.
    Vector w(1.0, 2.0, 2.0);
    w.calc_magnitude();
    std::cerr << "  magnitude of (1,2,2) = " << w.magnitude << " (expected 3)\n";
    EXPECT_NEAR(w.magnitude, 3.0, kTol) << "calc_magnitude: should be 3";
}

/*
 * test_normalize()
 * ----------------
 * Verifies that normalize() scales the vector to unit length and leaves the
 * direction unchanged. Also verifies the zero-vector special case.
 */
void test_normalize()
{
    std::cerr << "\n[TEST] test_normalize in class_Vector.cpp\n";

    // Normalizing (3,4,0) should give (0.6, 0.8, 0) with magnitude 1.
    Vector v(3.0, 4.0, 0.0);
    v.normalize();
    std::cerr << "  normalized (3,4,0) = (" << v.x << "," << v.y << "," << v.z
              << "), magnitude = " << v.magnitude << "\n";
    EXPECT_NEAR(v.x, 0.6, kTol) << "normalize: x should be 0.6";
    EXPECT_NEAR(v.y, 0.8, kTol) << "normalize: y should be 0.8";
    EXPECT_NEAR(v.z, 0.0, kTol) << "normalize: z should be 0.0";
    EXPECT_NEAR(v.magnitude, 1.0, kTol) << "normalize: resulting magnitude should be 1";

    // Normalizing a zero vector: the implementation sets magnitude to 1 to avoid
    // NaNs, so the components remain 0 and the recomputed magnitude is 0.
    Vector zero(0.0, 0.0, 0.0);
    zero.normalize();
    std::cerr << "  normalized zero vector = (" << zero.x << "," << zero.y << ","
              << zero.z << "), magnitude = " << zero.magnitude << "\n";
    EXPECT_NEAR(zero.x, 0.0, kTol) << "normalize: zero vector x stays 0";
    EXPECT_NEAR(zero.y, 0.0, kTol) << "normalize: zero vector y stays 0";
    EXPECT_NEAR(zero.z, 0.0, kTol) << "normalize: zero vector z stays 0";
    EXPECT_NEAR(zero.magnitude, 0.0, kTol) << "normalize: zero vector magnitude recomputed to 0";
}

/*
 * test_vector_projection()
 * ------------------------
 * Verifies the vector_projection() method, which returns the component of the
 * current vector that is orthogonal to the supplied normal vector.
 * result = this - proj_onto_normal(this)
 */
void test_vector_projection()
{
    std::cerr << "\n[TEST] test_vector_projection in class_Vector.cpp\n";

    // Project (2,3,0) with normal along x-axis (1,0,0).
    // proj onto x-axis is (2,0,0), so result = (2,3,0) - (2,0,0) = (0,3,0).
    Vector v(2.0, 3.0, 0.0);
    Vector normal(1.0, 0.0, 0.0);
    Vector result = v.vector_projection(normal);
    std::cerr << "  projection-removed of (2,3,0) w.r.t normal (1,0,0) = ("
              << result.x << "," << result.y << "," << result.z << ")\n";
    EXPECT_NEAR(result.x, 0.0, kTol) << "vector_projection: x component removed";
    EXPECT_NEAR(result.y, 3.0, kTol) << "vector_projection: y component preserved";
    EXPECT_NEAR(result.z, 0.0, kTol) << "vector_projection: z component preserved";

    // A vector already orthogonal to the normal should be unchanged.
    Vector ortho(0.0, 5.0, 0.0);
    Vector n2(1.0, 0.0, 0.0);
    Vector result2 = ortho.vector_projection(n2);
    std::cerr << "  projection-removed of (0,5,0) w.r.t normal (1,0,0) = ("
              << result2.x << "," << result2.y << "," << result2.z << ")\n";
    EXPECT_NEAR(result2.x, 0.0, kTol) << "vector_projection: orthogonal x unchanged";
    EXPECT_NEAR(result2.y, 5.0, kTol) << "vector_projection: orthogonal y unchanged";
    EXPECT_NEAR(result2.z, 0.0, kTol) << "vector_projection: orthogonal z unchanged";
}

/* ----------------------------------------------------------------------------
 * Google Test wrappers. Each TEST simply calls the corresponding test_* helper
 * so that assertions are captured by the framework while keeping the verbose,
 * grouped structure above.
 * -------------------------------------------------------------------------- */

TEST(VectorTest, Constructors)       { test_constructors(); }
TEST(VectorTest, OutputOperator)     { test_output_operator(); }
TEST(VectorTest, PlusOperator)       { test_plus_operator(); }
TEST(VectorTest, Cross)              { test_cross(); }
TEST(VectorTest, Dot)                { test_dot(); }
TEST(VectorTest, DotTheta)           { test_dot_theta(); }
TEST(VectorTest, CalcMagnitude)      { test_calc_magnitude(); }
TEST(VectorTest, Normalize)          { test_normalize(); }
TEST(VectorTest, VectorProjection)   { test_vector_projection(); }