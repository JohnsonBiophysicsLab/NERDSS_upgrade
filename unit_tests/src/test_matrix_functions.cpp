/*! \file test_matrix_functions.cpp
 *
 * ### Unit test for ../src/math/matrix_functions.cpp
 *
 * The file under test provides three free functions declared in
 * include/math/matrix.hpp:
 *
 *   1. Vector matrix_rotate(Vector& vec, std::array<double,9>& M)
 *        - Applies a 3x3 row-major matrix M to a Vector.
 *
 *   2. std::array<double,9> create_euler_rotation_matrix(double x, double y, double z)
 *        - Builds a Tait-Bryan (Euler) rotation matrix from three angles.
 *
 *   3. std::array<double,9> create_euler_rotation_matrix(const Coord& angles)
 *        - Overload of (2) that reads the angles from a Coord's x/y/z members.
 *
 * The tests below check:
 *   - the exact row-major arithmetic performed by matrix_rotate,
 *   - identity / basis-vector behaviour,
 *   - the algebraic form of each of the 9 matrix entries for the Euler builder,
 *   - the special cases of zero angles (identity) and single-axis rotations,
 *   - that the Coord overload produces bit-identical results to the
 *     three-double overload,
 *   - that the generated matrix is orthonormal with determinant +1 (a proper
 *     rotation), and that rotating a vector preserves its magnitude.
 *
 * Verbose progress and criteria are printed to stderr so a reader can follow
 * exactly which function is under test and what is being asserted.
 */

#include "classes/class_Coord.hpp"
#include "classes/class_Vector.hpp"
#include "math/matrix.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with mtxfn_ to avoid collisions in the combined suite)
// -----------------------------------------------------------------------------
namespace {

//! Tolerance used for floating point comparisons of trigonometric results.
constexpr double kMtxfnTol = 1e-12;

/*! \brief Pretty-print a 3x3 row-major matrix to stderr for diagnostics. */
void mtxfn_print_matrix(const char* label, const std::array<double, 9>& M)
{
    std::cerr << "    " << label << " =\n";
    for (int row = 0; row < 3; ++row) {
        std::cerr << "      [ " << M[row * 3 + 0] << ", " << M[row * 3 + 1] << ", "
                  << M[row * 3 + 2] << " ]\n";
    }
}

/*! \brief Print a Vector to stderr for diagnostics. */
void mtxfn_print_vector(const char* label, const Vector& v)
{
    std::cerr << "    " << label << " = (" << v.x << ", " << v.y << ", " << v.z << ")\n";
}

/*! \brief Compute the determinant of a 3x3 row-major matrix. */
double mtxfn_determinant(const std::array<double, 9>& M)
{
    return M[0] * (M[4] * M[8] - M[5] * M[7])
        - M[1] * (M[3] * M[8] - M[5] * M[6])
        + M[2] * (M[3] * M[7] - M[4] * M[6]);
}

/*! \brief Dot product of two rows (0,1,2) of a 3x3 row-major matrix. */
double mtxfn_row_dot(const std::array<double, 9>& M, int r1, int r2)
{
    return M[r1 * 3 + 0] * M[r2 * 3 + 0]
        + M[r1 * 3 + 1] * M[r2 * 3 + 1]
        + M[r1 * 3 + 2] * M[r2 * 3 + 2];
}

} // namespace

// -----------------------------------------------------------------------------
// Test: matrix_rotate performs the exact row-major matrix-vector product.
//
// matrix_rotate returns:
//    x' = M[0]*x + M[1]*y + M[2]*z
//    y' = M[3]*x + M[4]*y + M[5]*z
//    z' = M[6]*x + M[7]*y + M[8]*z
//
// We use an arbitrary (non-rotation) matrix with distinct integer entries so
// that any mis-ordering of indices would immediately produce a wrong answer.
// -----------------------------------------------------------------------------
void mtxfn_test_matrix_rotate_arithmetic()
{
    std::cerr << "\n[TEST] mtxfn_test_matrix_rotate_arithmetic\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Function:    matrix_rotate(Vector&, std::array<double,9>&)\n"
              << "  Criteria:    result must equal the row-major product M * vec,\n"
              << "               computed by hand from distinct integer entries.\n";

    // A deliberately non-symmetric matrix so index ordering matters.
    std::array<double, 9> M { 1.0, 2.0, 3.0,
                              4.0, 5.0, 6.0,
                              7.0, 8.0, 9.0 };
    Vector vec(1.0, 2.0, 3.0);

    mtxfn_print_matrix("M", M);
    mtxfn_print_vector("input vec", vec);

    Vector out = matrix_rotate(vec, M);
    mtxfn_print_vector("matrix_rotate(vec, M)", out);

    // Hand-computed expectations:
    //   x' = 1*1 + 2*2 + 3*3 = 14
    //   y' = 4*1 + 5*2 + 6*3 = 32
    //   z' = 7*1 + 8*2 + 9*3 = 50
    EXPECT_DOUBLE_EQ(out.x, 14.0) << "x' should be M[0]*x + M[1]*y + M[2]*z = 14";
    EXPECT_DOUBLE_EQ(out.y, 32.0) << "y' should be M[3]*x + M[4]*y + M[5]*z = 32";
    EXPECT_DOUBLE_EQ(out.z, 50.0) << "z' should be M[6]*x + M[7]*y + M[8]*z = 50";

    // matrix_rotate takes its arguments by non-const reference but must not
    // modify them; confirm the input vector is untouched.
    std::cerr << "  Also checking the input Vector is not mutated by the call.\n";
    EXPECT_DOUBLE_EQ(vec.x, 1.0) << "input vec.x should be unmodified";
    EXPECT_DOUBLE_EQ(vec.y, 2.0) << "input vec.y should be unmodified";
    EXPECT_DOUBLE_EQ(vec.z, 3.0) << "input vec.z should be unmodified";
}

// -----------------------------------------------------------------------------
// Test: matrix_rotate with the identity matrix returns the vector unchanged,
// and applied to the cartesian basis vectors it extracts matrix columns.
// -----------------------------------------------------------------------------
void mtxfn_test_matrix_rotate_identity_and_basis()
{
    std::cerr << "\n[TEST] mtxfn_test_matrix_rotate_identity_and_basis\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Function:    matrix_rotate(Vector&, std::array<double,9>&)\n"
              << "  Criteria:    (a) identity matrix leaves the vector unchanged;\n"
              << "               (b) rotating the unit basis vectors returns the\n"
              << "                   corresponding column of M.\n";

    // (a) Identity matrix must act as a no-op.
    std::cerr << "  -> Sub-check (a): identity matrix.\n";
    std::array<double, 9> I { 1.0, 0.0, 0.0,
                              0.0, 1.0, 0.0,
                              0.0, 0.0, 1.0 };
    Vector v(-3.5, 7.25, 0.125);
    Vector rotI = matrix_rotate(v, I);
    mtxfn_print_vector("identity * v", rotI);

    EXPECT_DOUBLE_EQ(rotI.x, -3.5) << "identity should preserve x";
    EXPECT_DOUBLE_EQ(rotI.y, 7.25) << "identity should preserve y";
    EXPECT_DOUBLE_EQ(rotI.z, 0.125) << "identity should preserve z";

    // (b) Rotating e_x, e_y, e_z pulls out columns 0, 1, 2 of M respectively.
    std::cerr << "  -> Sub-check (b): basis vectors extract columns of M.\n";
    std::array<double, 9> M { 10.0, 11.0, 12.0,
                              20.0, 21.0, 22.0,
                              30.0, 31.0, 32.0 };
    mtxfn_print_matrix("M", M);

    Vector ex(1.0, 0.0, 0.0);
    Vector ey(0.0, 1.0, 0.0);
    Vector ez(0.0, 0.0, 1.0);

    Vector colX = matrix_rotate(ex, M);
    Vector colY = matrix_rotate(ey, M);
    Vector colZ = matrix_rotate(ez, M);

    mtxfn_print_vector("M * e_x (column 0)", colX);
    mtxfn_print_vector("M * e_y (column 1)", colY);
    mtxfn_print_vector("M * e_z (column 2)", colZ);

    // Column 0 = (M[0], M[3], M[6])
    EXPECT_DOUBLE_EQ(colX.x, M[0]) << "M*e_x should give column 0 entry M[0]";
    EXPECT_DOUBLE_EQ(colX.y, M[3]) << "M*e_x should give column 0 entry M[3]";
    EXPECT_DOUBLE_EQ(colX.z, M[6]) << "M*e_x should give column 0 entry M[6]";

    // Column 1 = (M[1], M[4], M[7])
    EXPECT_DOUBLE_EQ(colY.x, M[1]) << "M*e_y should give column 1 entry M[1]";
    EXPECT_DOUBLE_EQ(colY.y, M[4]) << "M*e_y should give column 1 entry M[4]";
    EXPECT_DOUBLE_EQ(colY.z, M[7]) << "M*e_y should give column 1 entry M[7]";

    // Column 2 = (M[2], M[5], M[8])
    EXPECT_DOUBLE_EQ(colZ.x, M[2]) << "M*e_z should give column 2 entry M[2]";
    EXPECT_DOUBLE_EQ(colZ.y, M[5]) << "M*e_z should give column 2 entry M[5]";
    EXPECT_DOUBLE_EQ(colZ.z, M[8]) << "M*e_z should give column 2 entry M[8]";

    // A zero vector must always map to the zero vector.
    std::cerr << "  -> Sub-check (c): zero vector maps to zero vector.\n";
    Vector zero(0.0, 0.0, 0.0);
    Vector rotZero = matrix_rotate(zero, M);
    EXPECT_DOUBLE_EQ(rotZero.x, 0.0) << "M * 0 should have x == 0";
    EXPECT_DOUBLE_EQ(rotZero.y, 0.0) << "M * 0 should have y == 0";
    EXPECT_DOUBLE_EQ(rotZero.z, 0.0) << "M * 0 should have z == 0";
}

// -----------------------------------------------------------------------------
// Test: create_euler_rotation_matrix(double,double,double) with all-zero angles
// must produce the 3x3 identity matrix (sin=0, cos=1 everywhere).
// -----------------------------------------------------------------------------
void mtxfn_test_euler_zero_angles_is_identity()
{
    std::cerr << "\n[TEST] mtxfn_test_euler_zero_angles_is_identity\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Function:    create_euler_rotation_matrix(double, double, double)\n"
              << "  Criteria:    with x = y = z = 0 all sines vanish and all cosines\n"
              << "               are 1, so the result must be exactly the identity.\n";

    std::array<double, 9> M = create_euler_rotation_matrix(0.0, 0.0, 0.0);
    mtxfn_print_matrix("create_euler_rotation_matrix(0,0,0)", M);

    const std::array<double, 9> expected { 1.0, 0.0, 0.0,
                                           0.0, 1.0, 0.0,
                                           0.0, 0.0, 1.0 };

    for (int i = 0; i < 9; ++i) {
        EXPECT_DOUBLE_EQ(M[i], expected[i])
            << "Element M[" << i << "] should equal identity element " << expected[i];
    }
}

// -----------------------------------------------------------------------------
// Test: verify every one of the nine matrix entries against the closed-form
// expressions used in the implementation, for a set of generic (non-special)
// angles. This catches any transposition or sign error in a single element.
// -----------------------------------------------------------------------------
void mtxfn_test_euler_entries_match_formula()
{
    std::cerr << "\n[TEST] mtxfn_test_euler_entries_match_formula\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Function:    create_euler_rotation_matrix(double, double, double)\n"
              << "  Criteria:    each of the 9 entries must match the documented\n"
              << "               Tait-Bryan formula, evaluated independently here.\n";

    // Generic angles chosen so no sine or cosine is 0, 1 or -1.
    const double ax = 0.3;
    const double ay = -0.7;
    const double az = 1.1;

    std::cerr << "  Angles used: x = " << ax << ", y = " << ay << ", z = " << az << "\n";

    std::array<double, 9> M = create_euler_rotation_matrix(ax, ay, az);
    mtxfn_print_matrix("M", M);

    // Recompute the expected entries independently of the implementation.
    const double sx = std::sin(ax);
    const double cx = std::cos(ax);
    const double sy = std::sin(ay);
    const double cy = std::cos(ay);
    const double sz = std::sin(az);
    const double cz = std::cos(az);

    EXPECT_NEAR(M[0], cz * cy, kMtxfnTol) << "M[0] should be cz*cy";
    EXPECT_NEAR(M[1], cz * sx * sy - sz * cx, kMtxfnTol) << "M[1] should be cz*sx*sy - sz*cx";
    EXPECT_NEAR(M[2], cz * sy * cx + sz * sx, kMtxfnTol) << "M[2] should be cz*sy*cx + sz*sx";
    EXPECT_NEAR(M[3], sz * cy, kMtxfnTol) << "M[3] should be sz*cy";
    EXPECT_NEAR(M[4], sz * sx * sy + cz * cx, kMtxfnTol) << "M[4] should be sz*sx*sy + cz*cx";
    EXPECT_NEAR(M[5], sz * sy * cx - cz * sx, kMtxfnTol) << "M[5] should be sz*sy*cx - cz*sx";
    EXPECT_NEAR(M[6], -sy, kMtxfnTol) << "M[6] should be -sy";
    EXPECT_NEAR(M[7], cy * sx, kMtxfnTol) << "M[7] should be cy*sx";
    EXPECT_NEAR(M[8], cy * cx, kMtxfnTol) << "M[8] should be cy*cx";
}

// -----------------------------------------------------------------------------
// Test: single-axis rotations should reduce to the familiar 2D rotation blocks.
// This validates the sign convention of each axis independently.
// -----------------------------------------------------------------------------
void mtxfn_test_euler_single_axis_rotations()
{
    std::cerr << "\n[TEST] mtxfn_test_euler_single_axis_rotations\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Function:    create_euler_rotation_matrix(double, double, double)\n"
              << "  Criteria:    rotating about a single axis must leave that axis\n"
              << "               fixed and rotate the other two by the expected angle.\n";

    const double ang = M_PI / 2.0; // 90 degrees: sin = 1, cos ~ 0.
    const double s = std::sin(ang);
    const double c = std::cos(ang);

    // ---- Rotation about x only (y = z = 0) ----------------------------------
    // Expected form: [ 1  0   0 ; 0  cx -sx ; 0  sx  cx ]
    std::cerr << "  -> Sub-check: pure X rotation of +90 degrees.\n";
    std::array<double, 9> Mx = create_euler_rotation_matrix(ang, 0.0, 0.0);
    mtxfn_print_matrix("Mx", Mx);

    EXPECT_NEAR(Mx[0], 1.0, kMtxfnTol) << "pure-X: M[0] should be 1";
    EXPECT_NEAR(Mx[1], 0.0, kMtxfnTol) << "pure-X: M[1] should be 0";
    EXPECT_NEAR(Mx[2], 0.0, kMtxfnTol) << "pure-X: M[2] should be 0";
    EXPECT_NEAR(Mx[3], 0.0, kMtxfnTol) << "pure-X: M[3] should be 0";
    EXPECT_NEAR(Mx[4], c, kMtxfnTol) << "pure-X: M[4] should be cos(x)";
    EXPECT_NEAR(Mx[5], -s, kMtxfnTol) << "pure-X: M[5] should be -sin(x)";
    EXPECT_NEAR(Mx[6], 0.0, kMtxfnTol) << "pure-X: M[6] should be 0";
    EXPECT_NEAR(Mx[7], s, kMtxfnTol) << "pure-X: M[7] should be sin(x)";
    EXPECT_NEAR(Mx[8], c, kMtxfnTol) << "pure-X: M[8] should be cos(x)";

    // The x axis itself must be invariant under a pure x rotation.
    Vector ex(1.0, 0.0, 0.0);
    Vector rotEx = matrix_rotate(ex, Mx);
    mtxfn_print_vector("Mx * e_x", rotEx);
    EXPECT_NEAR(rotEx.x, 1.0, kMtxfnTol) << "pure-X rotation must fix e_x (x component)";
    EXPECT_NEAR(rotEx.y, 0.0, kMtxfnTol) << "pure-X rotation must fix e_x (y component)";
    EXPECT_NEAR(rotEx.z, 0.0, kMtxfnTol) << "pure-X rotation must fix e_x (z component)";

    // e_y must map onto +e_z for a +90 degree x rotation.
    Vector ey(0.0, 1.0, 0.0);
    Vector rotEy = matrix_rotate(ey, Mx);
    mtxfn_print_vector("Mx * e_y", rotEy);
    EXPECT_NEAR(rotEy.x, 0.0, kMtxfnTol) << "Mx*e_y should have x == 0";
    EXPECT_NEAR(rotEy.y, 0.0, kMtxfnTol) << "Mx*e_y should have y == 0 (cos 90 = 0)";
    EXPECT_NEAR(rotEy.z, 1.0, kMtxfnTol) << "Mx*e_y should have z == 1 (maps to +z)";

    // ---- Rotation about y only (x = z = 0) ----------------------------------
    // Expected form: [ cy 0 sy ; 0 1 0 ; -sy 0 cy ]
    std::cerr << "  -> Sub-check: pure Y rotation of +90 degrees.\n";
    std::array<double, 9> My = create_euler_rotation_matrix(0.0, ang, 0.0);
    mtxfn_print_matrix("My", My);

    EXPECT_NEAR(My[0], c, kMtxfnTol) << "pure-Y: M[0] should be cos(y)";
    EXPECT_NEAR(My[1], 0.0, kMtxfnTol) << "pure-Y: M[1] should be 0";
    EXPECT_NEAR(My[2], s, kMtxfnTol) << "pure-Y: M[2] should be sin(y)";
    EXPECT_NEAR(My[3], 0.0, kMtxfnTol) << "pure-Y: M[3] should be 0";
    EXPECT_NEAR(My[4], 1.0, kMtxfnTol) << "pure-Y: M[4] should be 1";
    EXPECT_NEAR(My[5], 0.0, kMtxfnTol) << "pure-Y: M[5] should be 0";
    EXPECT_NEAR(My[6], -s, kMtxfnTol) << "pure-Y: M[6] should be -sin(y)";
    EXPECT_NEAR(My[7], 0.0, kMtxfnTol) << "pure-Y: M[7] should be 0";
    EXPECT_NEAR(My[8], c, kMtxfnTol) << "pure-Y: M[8] should be cos(y)";

    Vector rotEyAxis = matrix_rotate(ey, My);
    mtxfn_print_vector("My * e_y", rotEyAxis);
    EXPECT_NEAR(rotEyAxis.x, 0.0, kMtxfnTol) << "pure-Y rotation must fix e_y (x component)";
    EXPECT_NEAR(rotEyAxis.y, 1.0, kMtxfnTol) << "pure-Y rotation must fix e_y (y component)";
    EXPECT_NEAR(rotEyAxis.z, 0.0, kMtxfnTol) << "pure-Y rotation must fix e_y (z component)";

    // ---- Rotation about z only (x = y = 0) ----------------------------------
    // Expected form: [ cz -sz 0 ; sz cz 0 ; 0 0 1 ]
    std::cerr << "  -> Sub-check: pure Z rotation of +90 degrees.\n";
    std::array<double, 9> Mz = create_euler_rotation_matrix(0.0, 0.0, ang);
    mtxfn_print_matrix("Mz", Mz);

    EXPECT_NEAR(Mz[0], c, kMtxfnTol) << "pure-Z: M[0] should be cos(z)";
    EXPECT_NEAR(Mz[1], -s, kMtxfnTol) << "pure-Z: M[1] should be -sin(z)";
    EXPECT_NEAR(Mz[2], 0.0, kMtxfnTol) << "pure-Z: M[2] should be 0";
    EXPECT_NEAR(Mz[3], s, kMtxfnTol) << "pure-Z: M[3] should be sin(z)";
    EXPECT_NEAR(Mz[4], c, kMtxfnTol) << "pure-Z: M[4] should be cos(z)";
    EXPECT_NEAR(Mz[5], 0.0, kMtxfnTol) << "pure-Z: M[5] should be 0";
    EXPECT_NEAR(Mz[6], 0.0, kMtxfnTol) << "pure-Z: M[6] should be 0";
    EXPECT_NEAR(Mz[7], 0.0, kMtxfnTol) << "pure-Z: M[7] should be 0";
    EXPECT_NEAR(Mz[8], 1.0, kMtxfnTol) << "pure-Z: M[8] should be 1";

    // e_x must map onto +e_y for a +90 degree z rotation.
    Vector rotExZ = matrix_rotate(ex, Mz);
    mtxfn_print_vector("Mz * e_x", rotExZ);
    EXPECT_NEAR(rotExZ.x, 0.0, kMtxfnTol) << "Mz*e_x should have x == 0 (cos 90 = 0)";
    EXPECT_NEAR(rotExZ.y, 1.0, kMtxfnTol) << "Mz*e_x should have y == 1 (maps to +y)";
    EXPECT_NEAR(rotExZ.z, 0.0, kMtxfnTol) << "Mz*e_x should have z == 0";
}

// -----------------------------------------------------------------------------
// Test: the Coord overload must produce results identical to the three-double
// overload, since both implementations share the same formula.
// -----------------------------------------------------------------------------
void mtxfn_test_euler_coord_overload_matches()
{
    std::cerr << "\n[TEST] mtxfn_test_euler_coord_overload_matches\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Function:    create_euler_rotation_matrix(const Coord&)\n"
              << "  Criteria:    the Coord overload must reproduce, bit-for-bit,\n"
              << "               the three-double overload given the same angles.\n";

    // A handful of angle triplets, including negatives and values > pi.
    const double angleSets[4][3] = {
        { 0.0, 0.0, 0.0 },
        { 0.25, -1.3, 2.7 },
        { -M_PI / 3.0, M_PI / 4.0, -M_PI / 6.0 },
        { 5.0, -4.2, 3.9 }
    };

    for (int setItr = 0; setItr < 4; ++setItr) {
        const double ax = angleSets[setItr][0];
        const double ay = angleSets[setItr][1];
        const double az = angleSets[setItr][2];

        std::cerr << "  -> Angle set " << setItr << ": (" << ax << ", " << ay << ", "
                  << az << ")\n";

        std::array<double, 9> fromDoubles = create_euler_rotation_matrix(ax, ay, az);

        // Build the same angles as a Coord and call the overload.
        Coord angles(ax, ay, az);
        std::array<double, 9> fromCoord = create_euler_rotation_matrix(angles);

        for (int i = 0; i < 9; ++i) {
            EXPECT_DOUBLE_EQ(fromCoord[i], fromDoubles[i])
                << "Angle set " << setItr << ": element " << i
                << " from the Coord overload must equal the double overload";
        }
    }
}

// -----------------------------------------------------------------------------
// Test: the generated Euler matrix must be a proper rotation, i.e. orthonormal
// with determinant +1. This is a strong structural property that would break
// under nearly any algebraic mistake in the construction.
// -----------------------------------------------------------------------------
void mtxfn_test_euler_matrix_is_proper_rotation()
{
    std::cerr << "\n[TEST] mtxfn_test_euler_matrix_is_proper_rotation\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Function:    create_euler_rotation_matrix(double, double, double)\n"
              << "  Criteria:    rows must be unit length and mutually orthogonal,\n"
              << "               and the determinant must be +1 (proper rotation).\n";

    const double angleSets[3][3] = {
        { 0.4, 0.9, -1.2 },
        { -2.0, 0.15, 3.0 },
        { 1.0, -1.0, 0.5 }
    };

    for (int setItr = 0; setItr < 3; ++setItr) {
        const double ax = angleSets[setItr][0];
        const double ay = angleSets[setItr][1];
        const double az = angleSets[setItr][2];

        std::cerr << "  -> Angle set " << setItr << ": (" << ax << ", " << ay << ", "
                  << az << ")\n";

        std::array<double, 9> M = create_euler_rotation_matrix(ax, ay, az);

        // Each row must have unit length.
        for (int row = 0; row < 3; ++row) {
            const double normSq = mtxfn_row_dot(M, row, row);
            std::cerr << "     |row " << row << "|^2 = " << normSq << "\n";
            EXPECT_NEAR(normSq, 1.0, 1e-10)
                << "Row " << row << " of the rotation matrix should be a unit vector";
        }

        // Distinct rows must be orthogonal (dot product zero).
        const double dot01 = mtxfn_row_dot(M, 0, 1);
        const double dot02 = mtxfn_row_dot(M, 0, 2);
        const double dot12 = mtxfn_row_dot(M, 1, 2);
        std::cerr << "     row0.row1 = " << dot01 << ", row0.row2 = " << dot02
                  << ", row1.row2 = " << dot12 << "\n";

        EXPECT_NEAR(dot01, 0.0, 1e-10) << "Rows 0 and 1 should be orthogonal";
        EXPECT_NEAR(dot02, 0.0, 1e-10) << "Rows 0 and 2 should be orthogonal";
        EXPECT_NEAR(dot12, 0.0, 1e-10) << "Rows 1 and 2 should be orthogonal";

        // Determinant of a proper rotation is exactly +1.
        const double det = mtxfn_determinant(M);
        std::cerr << "     det(M) = " << det << "\n";
        EXPECT_NEAR(det, 1.0, 1e-10)
            << "Determinant should be +1 for a proper (non-reflecting) rotation";
    }
}

// -----------------------------------------------------------------------------
// Test: because the Euler matrix is a rotation, applying it via matrix_rotate
// must preserve the length of any vector.
// -----------------------------------------------------------------------------
void mtxfn_test_rotation_preserves_magnitude()
{
    std::cerr << "\n[TEST] mtxfn_test_rotation_preserves_magnitude\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Functions:   create_euler_rotation_matrix + matrix_rotate\n"
              << "  Criteria:    rotating a vector must leave its magnitude\n"
              << "               unchanged (rotations are isometries).\n";

    // A rotation with all three angles non-trivial.
    std::array<double, 9> M = create_euler_rotation_matrix(0.6, -1.1, 2.2);
    mtxfn_print_matrix("M", M);

    // Several test vectors, including one along a single axis.
    Vector testVecs[3] = {
        Vector(3.0, 4.0, 12.0), // magnitude 13
        Vector(-1.0, 0.0, 0.0), // magnitude 1
        Vector(0.5, -2.5, 7.25)
    };

    for (int i = 0; i < 3; ++i) {
        Vector v = testVecs[i];
        v.calc_magnitude();
        const double beforeMag = v.magnitude;

        Vector rotated = matrix_rotate(v, M);
        rotated.calc_magnitude();
        const double afterMag = rotated.magnitude;

        std::cerr << "  -> vector " << i << ": |v| before = " << beforeMag
                  << ", |Mv| after = " << afterMag << "\n";

        EXPECT_NEAR(afterMag, beforeMag, 1e-10)
            << "Rotation should preserve the magnitude of test vector " << i;
    }

    // A rotation must also preserve the dot product between two vectors.
    std::cerr << "  -> Also checking that the dot product of two vectors is preserved.\n";
    Vector a(1.0, 2.0, 3.0);
    Vector b(-2.0, 0.5, 4.0);
    const double dotBefore = a.dot(b);

    Vector aRot = matrix_rotate(a, M);
    Vector bRot = matrix_rotate(b, M);
    const double dotAfter = aRot.dot(bRot);

    std::cerr << "     a.b before = " << dotBefore << ", (Ma).(Mb) after = "
              << dotAfter << "\n";
    EXPECT_NEAR(dotAfter, dotBefore, 1e-10)
        << "Rotation should preserve the dot product between two vectors";
}

// -----------------------------------------------------------------------------
// Test: an Euler matrix built with the negated angles applied in the reverse
// composition order acts as the inverse of the original. We verify this the
// simple way: rotating by M and then by its transpose returns the original
// vector, since R^T == R^-1 for an orthonormal rotation matrix.
// -----------------------------------------------------------------------------
void mtxfn_test_transpose_is_inverse()
{
    std::cerr << "\n[TEST] mtxfn_test_transpose_is_inverse\n"
              << "  Source file: src/math/matrix_functions.cpp\n"
              << "  Functions:   create_euler_rotation_matrix + matrix_rotate\n"
              << "  Criteria:    since the generated matrix is orthonormal,\n"
              << "               applying M then M^T must recover the original\n"
              << "               vector to within floating point tolerance.\n";

    std::array<double, 9> M = create_euler_rotation_matrix(-0.85, 1.4, 0.33);
    mtxfn_print_matrix("M", M);

    // Build the transpose explicitly (row-major swap of off-diagonal entries).
    std::array<double, 9> Mt { M[0], M[3], M[6],
                               M[1], M[4], M[7],
                               M[2], M[5], M[8] };
    mtxfn_print_matrix("M^T", Mt);

    Vector original(2.0, -5.0, 0.75);
    mtxfn_print_vector("original", original);

    Vector rotated = matrix_rotate(original, M);
    mtxfn_print_vector("M * original", rotated);

    Vector restored = matrix_rotate(rotated, Mt);
    mtxfn_print_vector("M^T * (M * original)", restored);

    EXPECT_NEAR(restored.x, original.x, 1e-10)
        << "x component should be recovered after M then M^T";
    EXPECT_NEAR(restored.y, original.y, 1e-10)
        << "y component should be recovered after M then M^T";
    EXPECT_NEAR(restored.z, original.z, 1e-10)
        << "z component should be recovered after M then M^T";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is executed inside its own TEST so a
// failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(MatrixFunctions, MatrixRotateArithmetic) { mtxfn_test_matrix_rotate_arithmetic(); }
TEST(MatrixFunctions, MatrixRotateIdentityAndBasis) { mtxfn_test_matrix_rotate_identity_and_basis(); }
TEST(MatrixFunctions, EulerZeroAnglesIsIdentity) { mtxfn_test_euler_zero_angles_is_identity(); }
TEST(MatrixFunctions, EulerEntriesMatchFormula) { mtxfn_test_euler_entries_match_formula(); }
TEST(MatrixFunctions, EulerSingleAxisRotations) { mtxfn_test_euler_single_axis_rotations(); }
TEST(MatrixFunctions, EulerCoordOverloadMatches) { mtxfn_test_euler_coord_overload_matches(); }
TEST(MatrixFunctions, EulerMatrixIsProperRotation) { mtxfn_test_euler_matrix_is_proper_rotation(); }
TEST(MatrixFunctions, RotationPreservesMagnitude) { mtxfn_test_rotation_preserves_magnitude(); }
TEST(MatrixFunctions, TransposeIsInverse) { mtxfn_test_transpose_is_inverse(); }