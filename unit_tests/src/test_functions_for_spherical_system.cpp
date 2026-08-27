/*! \file test_functions_for_spherical_system.cpp
 *
 * ### Unit tests for src/reactions/functions_for_spherical_system.cpp
 *
 * The file under test contains all of the geometry helpers that NERDSS uses
 * when the simulation boundary is a sphere.  The functions exercised here are:
 *
 *   - double radius(Coord)
 *   - Coord  find_spherical_coords(Coord)
 *   - Coord  find_cardesian_coords(Coord)
 *   - double theta_plus(double, double)
 *   - double phi_plus(double, double)
 *   - Coord  angle_plus(Coord, Coord)
 *   - Coord  find_position_after_association(double, Coord, Coord, double, double)
 *   - std::array<double,9> inner_coord_set(Coord, Coord)
 *   - std::array<double,9> inner_coord_set_new(Coord, Coord)
 *   - std::array<double,3> calculate_inner_coord_coefficients(Coord, Coord, std::array<double,9>)
 *   - Coord  translate_on_sphere(Coord, Coord, Coord, std::array<double,9>, std::array<double,9>)
 *   - Coord  rotate_on_sphere(Coord, Coord, std::array<double,9>, double)
 *   - double calc_bindRadius2D(double, Coord)
 *   - void   set_memProtein_sphere(Complex, Molecule&, std::vector<Molecule>, const Membrane)
 *   - void   find_Lipid_sphere(Complex, Molecule&, std::vector<Molecule>, const Membrane)
 *
 * Notes on paths that are deliberately NOT exercised:
 *   - find_position_after_association() calls exit(1) if it produces a NaN
 *     position (which happens for degenerate input where the plane through the
 *     origin and the two interfaces is undefined).  All inputs below are
 *     chosen so the plane is well defined.
 *   - set_memProtein_sphere() and find_Lipid_sphere() call exit(1) when no
 *     lipid / implicit lipid can be located.  Every fixture below always
 *     supplies one.
 *   - calculate_inner_coord_coefficients() divides by quantities that vanish
 *     for axis-aligned bases in its "general" branch, so the general branch is
 *     only exercised with a generic (no zero component) orthonormal basis, and
 *     the special branches are exercised with the identity basis.
 *   - rotate_on_sphere() has an early-out when the projection of the target on
 *     the rotation axis has magnitude ~1.0; all targets below avoid that value.
 */

#include "reactions/association/functions_for_spherical_system.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Small local helpers (anonymous namespace + fss_ prefix so they cannot clash
// with anything else linked into the unit-test binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Euclidean length of a Coord. */
double fss_mag(const Coord& c) { return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z); }

/*! \brief Dot product of two Coords. */
double fss_dot(const Coord& a, const Coord& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/*! \brief Build a Coord of length R pointing along (x,y,z). */
Coord fss_on_sphere(double x, double y, double z, double R)
{
    const double m = std::sqrt(x * x + y * y + z * z);
    return Coord { R * x / m, R * y / m, R * z / m };
}

/*! \brief Great-circle arc length between two points assumed to lie on a sphere of radius R. */
double fss_arc(const Coord& a, const Coord& b, double R)
{
    double c = fss_dot(a, b) / (R * R);
    if (c > 1.0)
        c = 1.0;
    if (c < -1.0)
        c = -1.0;
    return R * std::acos(c);
}

/*! \brief Pull one of the three basis vectors (0=i, 1=j, 2=k) out of a 9 element crdset. */
Coord fss_basis(const std::array<double, 9>& s, int which)
{
    return Coord { s[3 * which], s[3 * which + 1], s[3 * which + 2] };
}

/*! \brief Assert that a 9-element crdset really is a right-handed orthonormal triad. */
void fss_expect_orthonormal(const std::array<double, 9>& s, const char* label)
{
    const Coord i = fss_basis(s, 0);
    const Coord j = fss_basis(s, 1);
    const Coord k = fss_basis(s, 2);

    std::cerr << "    " << label << " i = (" << i.x << ", " << i.y << ", " << i.z << ")\n";
    std::cerr << "    " << label << " j = (" << j.x << ", " << j.y << ", " << j.z << ")\n";
    std::cerr << "    " << label << " k = (" << k.x << ", " << k.y << ", " << k.z << ")\n";

    // Every basis vector is normalised at the end of the routine.
    EXPECT_NEAR(fss_mag(i), 1.0, 1e-12) << label << ": |i| must be 1";
    EXPECT_NEAR(fss_mag(j), 1.0, 1e-12) << label << ": |j| must be 1";
    EXPECT_NEAR(fss_mag(k), 1.0, 1e-12) << label << ": |k| must be 1";

    // Mutually perpendicular.
    EXPECT_NEAR(fss_dot(i, j), 0.0, 1e-12) << label << ": i.j must be 0";
    EXPECT_NEAR(fss_dot(i, k), 0.0, 1e-12) << label << ": i.k must be 0";
    EXPECT_NEAR(fss_dot(j, k), 0.0, 1e-12) << label << ": j.k must be 0";

    // Right handed: i x j == k.
    const Coord cross { i.y * j.z - i.z * j.y, i.z * j.x - i.x * j.z, i.x * j.y - i.y * j.x };
    EXPECT_NEAR(cross.x, k.x, 1e-12) << label << ": (i x j).x must equal k.x";
    EXPECT_NEAR(cross.y, k.y, 1e-12) << label << ": (i x j).y must equal k.y";
    EXPECT_NEAR(cross.z, k.z, 1e-12) << label << ": (i x j).z must equal k.z";
}

/*! \brief Build a plain, non-lipid molecule with a single interface. */
Molecule fss_make_protein(const Coord& com, const Coord& iface)
{
    Molecule mol;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.comCoord = com;
    mol.interfaceList.resize(1);
    mol.interfaceList[0].coord = iface;
    mol.tmpComCoord = com;
    mol.tmpICoords.clear();
    mol.tmpICoords.push_back(iface);
    return mol;
}

/*! \brief Build an explicit lipid molecule with a single interface. */
Molecule fss_make_lipid(const Coord& com, const Coord& iface)
{
    Molecule mol = fss_make_protein(com, iface);
    mol.isLipid = true;
    return mol;
}

} // namespace

// -----------------------------------------------------------------------------
// radius()
// -----------------------------------------------------------------------------
void fss_test_radius()
{
    std::cerr << "\n[TEST] fss_test_radius\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    radius(Coord)\n"
              << "  Criteria:    returns sqrt(x^2+y^2+z^2).\n";

    // Classic 3-4-0 right triangle -> length 5.
    std::cerr << "  -> radius((3,4,0)) should be 5\n";
    EXPECT_DOUBLE_EQ(radius(Coord { 3.0, 4.0, 0.0 }), 5.0) << "radius of (3,4,0) must be 5";

    // 1-2-2 -> length 3.
    std::cerr << "  -> radius((1,2,2)) should be 3\n";
    EXPECT_DOUBLE_EQ(radius(Coord { 1.0, 2.0, 2.0 }), 3.0) << "radius of (1,2,2) must be 3";

    // The origin has zero length.
    std::cerr << "  -> radius((0,0,0)) should be 0\n";
    EXPECT_DOUBLE_EQ(radius(Coord { 0.0, 0.0, 0.0 }), 0.0) << "radius of the origin must be 0";

    // Sign should not matter.
    std::cerr << "  -> radius((-3,0,-4)) should be 5 (signs ignored)\n";
    EXPECT_DOUBLE_EQ(radius(Coord { -3.0, 0.0, -4.0 }), 5.0) << "radius must be sign independent";
}

// -----------------------------------------------------------------------------
// find_spherical_coords()
// -----------------------------------------------------------------------------
void fss_test_find_spherical_coords()
{
    std::cerr << "\n[TEST] fss_test_find_spherical_coords\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    find_spherical_coords(Coord)\n"
              << "  Criteria:    returns (theta, phi, R) with the special north/south\n"
              << "               pole handling coded in the source.\n";

    // North pole: the source short circuits when z == R exactly.
    std::cerr << "  -> north pole (0,0,5): expect theta=0, phi=0, R=5\n";
    Coord north = find_spherical_coords(Coord { 0.0, 0.0, 5.0 });
    EXPECT_DOUBLE_EQ(north.x, 0.0) << "theta at the north pole must be 0";
    EXPECT_DOUBLE_EQ(north.y, 0.0) << "phi at the north pole must be 0";
    EXPECT_DOUBLE_EQ(north.z, 5.0) << "R at the north pole must be 5";

    // South pole: the source stores theta = -PI (not +PI) for z == -R.
    std::cerr << "  -> south pole (0,0,-5): expect theta=-PI, phi=0, R=5 (as coded)\n";
    Coord south = find_spherical_coords(Coord { 0.0, 0.0, -5.0 });
    EXPECT_DOUBLE_EQ(south.x, -M_PI) << "theta at the south pole is coded as -PI";
    EXPECT_DOUBLE_EQ(south.y, 0.0) << "phi at the south pole must be 0";
    EXPECT_DOUBLE_EQ(south.z, 5.0) << "R at the south pole must be 5";

    // +x axis: theta = PI/2, phi = 0.
    std::cerr << "  -> (1,0,0): expect theta=PI/2, phi=0, R=1\n";
    Coord px = find_spherical_coords(Coord { 1.0, 0.0, 0.0 });
    EXPECT_NEAR(px.x, M_PI / 2.0, 1e-12) << "theta on the +x axis must be PI/2";
    EXPECT_NEAR(px.y, 0.0, 1e-12) << "phi on the +x axis must be 0";
    EXPECT_NEAR(px.z, 1.0, 1e-12) << "R must be 1";

    // +y axis: phi = PI/2 (y >= 0 so no 2*PI - phi correction).
    std::cerr << "  -> (0,1,0): expect theta=PI/2, phi=PI/2\n";
    Coord py = find_spherical_coords(Coord { 0.0, 1.0, 0.0 });
    EXPECT_NEAR(py.x, M_PI / 2.0, 1e-12) << "theta on the +y axis must be PI/2";
    EXPECT_NEAR(py.y, M_PI / 2.0, 1e-12) << "phi on the +y axis must be PI/2";

    // -y axis: y < 0 triggers phi -> 2*PI - phi.
    std::cerr << "  -> (0,-1,0): expect phi=3*PI/2 because y<0 flips phi\n";
    Coord ny = find_spherical_coords(Coord { 0.0, -1.0, 0.0 });
    EXPECT_NEAR(ny.x, M_PI / 2.0, 1e-12) << "theta on the -y axis must be PI/2";
    EXPECT_NEAR(ny.y, 3.0 * M_PI / 2.0, 1e-12) << "phi on the -y axis must be 3*PI/2";

    // -x axis: phi = PI.
    std::cerr << "  -> (-1,0,0): expect phi=PI\n";
    Coord nx = find_spherical_coords(Coord { -1.0, 0.0, 0.0 });
    EXPECT_NEAR(nx.y, M_PI, 1e-12) << "phi on the -x axis must be PI";
}

// -----------------------------------------------------------------------------
// find_cardesian_coords() and a spherical <-> cartesian round trip.
// -----------------------------------------------------------------------------
void fss_test_find_cardesian_coords()
{
    std::cerr << "\n[TEST] fss_test_find_cardesian_coords\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    find_cardesian_coords(Coord)\n"
              << "  Criteria:    (theta,phi,R) -> (x,y,z) and a round trip through\n"
              << "               find_spherical_coords() returns the original point.\n";

    // theta = PI/2, phi = 0 puts the point on the +x axis at distance R.
    std::cerr << "  -> (theta=PI/2, phi=0, R=3) should map to (3,0,0)\n";
    Coord xyz = find_cardesian_coords(Coord { M_PI / 2.0, 0.0, 3.0 });
    EXPECT_NEAR(xyz.x, 3.0, 1e-12) << "x must be 3";
    EXPECT_NEAR(xyz.y, 0.0, 1e-12) << "y must be 0";
    EXPECT_NEAR(xyz.z, 0.0, 1e-12) << "z must be 0";

    // theta = 0 puts the point on the +z axis.
    std::cerr << "  -> (theta=0, phi=1.234, R=4) should map to (0,0,4)\n";
    Coord pole = find_cardesian_coords(Coord { 0.0, 1.234, 4.0 });
    EXPECT_NEAR(pole.x, 0.0, 1e-12) << "x must be 0 on the polar axis";
    EXPECT_NEAR(pole.y, 0.0, 1e-12) << "y must be 0 on the polar axis";
    EXPECT_NEAR(pole.z, 4.0, 1e-12) << "z must be R on the polar axis";

    // Round trip: cartesian -> spherical -> cartesian.
    std::cerr << "  -> round trip of a generic point (2, -3, 6)\n";
    Coord original { 2.0, -3.0, 6.0 };
    Coord back = find_cardesian_coords(find_spherical_coords(original));
    EXPECT_NEAR(back.x, original.x, 1e-10) << "round trip must restore x";
    EXPECT_NEAR(back.y, original.y, 1e-10) << "round trip must restore y";
    EXPECT_NEAR(back.z, original.z, 1e-10) << "round trip must restore z";
}

// -----------------------------------------------------------------------------
// theta_plus()
// -----------------------------------------------------------------------------
void fss_test_theta_plus()
{
    std::cerr << "\n[TEST] fss_test_theta_plus\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    theta_plus(double, double)\n"
              << "  Criteria:    plain sum, reflected back into [0,PI] when it leaves.\n";

    // In-range sum is untouched.
    std::cerr << "  -> 0.5 + 0.5 stays 1.0 (inside [0,PI])\n";
    EXPECT_NEAR(theta_plus(0.5, 0.5), 1.0, 1e-14) << "in-range sum must be unchanged";

    // Overshooting PI reflects: sum -> 2*PI - sum.
    std::cerr << "  -> 3.0 + 1.0 = 4.0 > PI, so expect 2*PI - 4.0\n";
    EXPECT_NEAR(theta_plus(3.0, 1.0), 2.0 * M_PI - 4.0, 1e-14) << "sum > PI must reflect to 2*PI - sum";

    // Undershooting 0 reflects: sum -> -sum.
    std::cerr << "  -> 0.2 + (-0.5) = -0.3 < 0, so expect 0.3\n";
    EXPECT_NEAR(theta_plus(0.2, -0.5), 0.3, 1e-14) << "sum < 0 must reflect to -sum";

    // Exactly PI is left alone (not > PI).
    std::cerr << "  -> 1.0 + (PI - 1.0) == PI stays PI\n";
    EXPECT_NEAR(theta_plus(1.0, M_PI - 1.0), M_PI, 1e-14) << "sum exactly PI must be unchanged";
}

// -----------------------------------------------------------------------------
// phi_plus()
// -----------------------------------------------------------------------------
void fss_test_phi_plus()
{
    std::cerr << "\n[TEST] fss_test_phi_plus\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    phi_plus(double, double)\n"
              << "  Criteria:    plain sum wrapped into [0,2*PI).\n";

    // In-range sum is untouched.
    std::cerr << "  -> 1.0 + 2.0 stays 3.0\n";
    EXPECT_NEAR(phi_plus(1.0, 2.0), 3.0, 1e-14) << "in-range sum must be unchanged";

    // Wrap down when past 2*PI.
    std::cerr << "  -> 6.0 + 1.0 = 7.0 > 2*PI, so expect 7.0 - 2*PI\n";
    EXPECT_NEAR(phi_plus(6.0, 1.0), 7.0 - 2.0 * M_PI, 1e-14) << "sum > 2*PI must wrap down";

    // Wrap up when negative.
    std::cerr << "  -> 0.5 + (-1.0) = -0.5 < 0, so expect 2*PI - 0.5\n";
    EXPECT_NEAR(phi_plus(0.5, -1.0), 2.0 * M_PI - 0.5, 1e-14) << "negative sum must wrap up";
}

// -----------------------------------------------------------------------------
// angle_plus()
// -----------------------------------------------------------------------------
void fss_test_angle_plus()
{
    std::cerr << "\n[TEST] fss_test_angle_plus\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    angle_plus(Coord, Coord)\n"
              << "  Criteria:    theta/phi combined with pole crossing corrections;\n"
              << "               the radius component is taken from the first argument.\n";

    // Simple in-range addition, radius comes from angle1.
    std::cerr << "  -> generic case (0.5,1.0,10) + (0.3,0.5,99) -> (0.8,1.5,10)\n";
    Coord simple = angle_plus(Coord { 0.5, 1.0, 10.0 }, Coord { 0.3, 0.5, 99.0 });
    EXPECT_NEAR(simple.x, 0.8, 1e-14) << "theta must simply add";
    EXPECT_NEAR(simple.y, 1.5, 1e-14) << "phi must simply add";
    EXPECT_DOUBLE_EQ(simple.z, 10.0) << "radius must be taken from the first argument";

    // theta overshoots PI: reflect theta and shift phi by PI.
    std::cerr << "  -> theta overshoot (3.0,1.0,10) + (1.0,0.5,0): theta=4.0>PI\n";
    Coord over = angle_plus(Coord { 3.0, 1.0, 10.0 }, Coord { 1.0, 0.5, 0.0 });
    EXPECT_NEAR(over.x, 2.0 * M_PI - 4.0, 1e-14) << "theta > PI must reflect to 2*PI - theta";
    EXPECT_NEAR(over.y, 1.5 + M_PI, 1e-13) << "phi must be shifted by +PI when crossing the pole";
    EXPECT_DOUBLE_EQ(over.z, 10.0) << "radius must be preserved";

    // theta undershoots 0: reflect theta and shift phi by -PI.
    std::cerr << "  -> theta undershoot (0.2,1.0,10) + (-0.5,0.5,0): theta=-0.3<0\n";
    Coord under = angle_plus(Coord { 0.2, 1.0, 10.0 }, Coord { -0.5, 0.5, 0.0 });
    EXPECT_NEAR(under.x, 0.3, 1e-14) << "theta < 0 must reflect to -theta";
    EXPECT_NEAR(under.y, 2.0 * M_PI + (1.5 - M_PI), 1e-13) << "phi must be shifted by -PI and wrapped";
    EXPECT_DOUBLE_EQ(under.z, 10.0) << "radius must be preserved";

    // theta lands exactly on 0: phi is zeroed by the source.
    std::cerr << "  -> theta lands exactly on 0: phi is forced to 0\n";
    Coord onPole = angle_plus(Coord { 0.5, 1.0, 7.0 }, Coord { -0.5, 0.5, 0.0 });
    EXPECT_NEAR(onPole.x, 0.0, 1e-14) << "theta must be exactly 0";
    EXPECT_NEAR(onPole.y, 0.0, 1e-14) << "phi must be forced to 0 on the pole";
    EXPECT_DOUBLE_EQ(onPole.z, 7.0) << "radius must be preserved";
}

// -----------------------------------------------------------------------------
// find_position_after_association()
// -----------------------------------------------------------------------------
void fss_test_find_position_after_association()
{
    std::cerr << "\n[TEST] fss_test_find_position_after_association\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    find_position_after_association(arc1, Iface1, Iface2, arc_total, bindRadius)\n"
              << "  Criteria:    the returned point stays on the sphere, sits at arc\n"
              << "               distance |arc1| from Iface1 along the great circle through\n"
              << "               Iface1 and Iface2, and the root chosen depends on whether\n"
              << "               bindRadius is smaller or larger than arc_total.\n";

    const double R = 10.0;

    // Two interfaces on the same sphere, separated only in azimuth so the plane
    // through the origin and both interfaces is well defined (no NaN exit path).
    const double dphi = 0.2;
    Coord iface1 { R * std::sin(M_PI / 4.0), 0.0, R * std::cos(M_PI / 4.0) };
    Coord iface2 { R * std::sin(M_PI / 4.0) * std::cos(dphi),
        R * std::sin(M_PI / 4.0) * std::sin(dphi),
        R * std::cos(M_PI / 4.0) };

    const double arcTotal = fss_arc(iface1, iface2, R);
    std::cerr << "  Setup: R = " << R << ", arc between the two interfaces = " << arcTotal << '\n';

    const double arc1 = 0.5;

    // --- Case 1: bindRadius < arc_total -> the root nearest Iface2 is taken,
    //     i.e. the moving interface travels *toward* its partner.
    std::cerr << "  -> case bindRadius(1.0) < arc_total: expect motion toward Iface2\n";
    Coord toward = find_position_after_association(arc1, iface1, iface2, arcTotal, 1.0);
    std::cerr << "     new position = (" << toward.x << ", " << toward.y << ", " << toward.z << ")\n";

    EXPECT_NEAR(fss_mag(toward), R, 1e-6) << "result must remain on the sphere of radius R";
    EXPECT_NEAR(fss_arc(toward, iface1, R), arc1, 1e-6)
        << "result must be arc1 away from Iface1 along the great circle";
    EXPECT_NEAR(fss_arc(toward, iface2, R), arcTotal - arc1, 1e-6)
        << "moving toward the partner shortens the remaining arc by arc1";

    // --- Case 2: bindRadius > arc_total -> the far root is taken, i.e. the
    //     moving interface travels *away* from its partner.
    std::cerr << "  -> case bindRadius(2.0) > arc_total: expect motion away from Iface2\n";
    Coord away = find_position_after_association(arc1, iface1, iface2, arcTotal, 2.0);
    std::cerr << "     new position = (" << away.x << ", " << away.y << ", " << away.z << ")\n";

    EXPECT_NEAR(fss_mag(away), R, 1e-6) << "result must remain on the sphere of radius R";
    EXPECT_NEAR(fss_arc(away, iface1, R), arc1, 1e-6)
        << "result must still be arc1 away from Iface1";
    EXPECT_NEAR(fss_arc(away, iface2, R), arcTotal + arc1, 1e-6)
        << "moving away from the partner lengthens the arc by arc1";

    // The two branches must produce different points.
    EXPECT_GT(std::sqrt(std::pow(toward.x - away.x, 2.0) + std::pow(toward.y - away.y, 2.0)
                  + std::pow(toward.z - away.z, 2.0)),
        1e-3)
        << "the two roots must be distinct points";

    // --- Case 3: the source takes std::abs(arc1), so a negative arc must give
    //     exactly the same answer as the positive one.
    std::cerr << "  -> negative arc1 (-0.5) must give the same point as +0.5\n";
    Coord negArc = find_position_after_association(-arc1, iface1, iface2, arcTotal, 1.0);
    EXPECT_NEAR(negArc.x, toward.x, 1e-10) << "arc1 sign must not matter (x)";
    EXPECT_NEAR(negArc.y, toward.y, 1e-10) << "arc1 sign must not matter (y)";
    EXPECT_NEAR(negArc.z, toward.z, 1e-10) << "arc1 sign must not matter (z)";
}

// -----------------------------------------------------------------------------
// inner_coord_set()
// -----------------------------------------------------------------------------
void fss_test_inner_coord_set()
{
    std::cerr << "\n[TEST] fss_test_inner_coord_set\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    inner_coord_set(com, comnew)\n"
              << "  Criteria:    a right-handed orthonormal triad whose i axis points\n"
              << "               along com; when com != comnew the j axis lies in the\n"
              << "               com/comnew plane pointing toward comnew.\n";

    const double R = 10.0;

    // --- Degenerate branch: com == comnew, non-polar. Fallback v = (0,0,1).
    std::cerr << "  -> degenerate branch, com == comnew == (10,0,0)\n";
    std::array<double, 9> same = inner_coord_set(Coord { R, 0.0, 0.0 }, Coord { R, 0.0, 0.0 });
    fss_expect_orthonormal(same, "same");
    EXPECT_NEAR(same[0], 1.0, 1e-12) << "i must point along +x for com = (10,0,0)";
    EXPECT_NEAR(same[1], 0.0, 1e-12) << "i.y must be 0";
    EXPECT_NEAR(same[2], 0.0, 1e-12) << "i.z must be 0";

    // --- Degenerate branch on the pole: the fallback vector becomes (-1,0,0).
    std::cerr << "  -> degenerate branch on the polar axis, com == comnew == (0,0,10)\n";
    std::array<double, 9> polar = inner_coord_set(Coord { 0.0, 0.0, R }, Coord { 0.0, 0.0, R });
    fss_expect_orthonormal(polar, "polar");
    EXPECT_NEAR(polar[2], 1.0, 1e-12) << "i must point along +z for a polar com";

    // --- General branch: com != comnew.
    std::cerr << "  -> general branch, com and comnew are distinct points on the sphere\n";
    Coord com = fss_on_sphere(1.0, 1.0, 1.0, R);
    Coord comnew = fss_on_sphere(1.05, 1.1, 0.95, R);
    std::array<double, 9> crdset = inner_coord_set(com, comnew);
    fss_expect_orthonormal(crdset, "crdset");

    // i must be the unit vector along com.
    const Coord i = fss_basis(crdset, 0);
    EXPECT_NEAR(i.x, com.x / R, 1e-12) << "i.x must be com.x/|com|";
    EXPECT_NEAR(i.y, com.y / R, 1e-12) << "i.y must be com.y/|com|";
    EXPECT_NEAR(i.z, com.z / R, 1e-12) << "i.z must be com.z/|com|";

    // k is normal to the plane, so it is perpendicular to both com and comnew.
    const Coord k = fss_basis(crdset, 2);
    EXPECT_NEAR(fss_dot(k, com), 0.0, 1e-10) << "k must be normal to com";
    EXPECT_NEAR(fss_dot(k, comnew), 0.0, 1e-10) << "k must be normal to comnew";

    // j lies in the plane and points along the direction of travel.
    const Coord j = fss_basis(crdset, 1);
    EXPECT_GT(fss_dot(j, comnew), 0.0) << "j must point toward comnew (direction of travel)";
}

// -----------------------------------------------------------------------------
// inner_coord_set_new()
// -----------------------------------------------------------------------------
void fss_test_inner_coord_set_new()
{
    std::cerr << "\n[TEST] fss_test_inner_coord_set_new\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    inner_coord_set_new(com, comnew)\n"
              << "  Criteria:    an orthonormal triad anchored on comnew; k stays normal\n"
              << "               to the com/comnew plane so the pair of triads describes\n"
              << "               a pure rotation about k.\n";

    const double R = 10.0;

    // --- Degenerate branch: com == comnew, i must follow comnew.
    std::cerr << "  -> degenerate branch, com == comnew == (0,10,0)\n";
    std::array<double, 9> same = inner_coord_set_new(Coord { 0.0, R, 0.0 }, Coord { 0.0, R, 0.0 });
    fss_expect_orthonormal(same, "same_new");
    EXPECT_NEAR(same[1], 1.0, 1e-12) << "i must point along +y for comnew = (0,10,0)";

    // --- General branch.
    std::cerr << "  -> general branch, distinct com and comnew\n";
    Coord com = fss_on_sphere(1.0, 1.0, 1.0, R);
    Coord comnew = fss_on_sphere(1.05, 1.1, 0.95, R);

    std::array<double, 9> crdset = inner_coord_set(com, comnew);
    std::array<double, 9> crdsetnew = inner_coord_set_new(com, comnew);
    fss_expect_orthonormal(crdsetnew, "crdsetnew");

    // The new i axis must be the unit vector along comnew.
    const Coord iNew = fss_basis(crdsetnew, 0);
    EXPECT_NEAR(iNew.x, comnew.x / R, 1e-10) << "new i.x must be comnew.x/|comnew|";
    EXPECT_NEAR(iNew.y, comnew.y / R, 1e-10) << "new i.y must be comnew.y/|comnew|";
    EXPECT_NEAR(iNew.z, comnew.z / R, 1e-10) << "new i.z must be comnew.z/|comnew|";

    // Both triads must share the same k axis: they only differ by a rotation
    // about the normal of the com/comnew plane.
    const Coord kOld = fss_basis(crdset, 2);
    const Coord kNew = fss_basis(crdsetnew, 2);
    std::cerr << "     k(old) . k(new) = " << fss_dot(kOld, kNew) << " (should be ~1)\n";
    EXPECT_NEAR(fss_dot(kOld, kNew), 1.0, 1e-8)
        << "both triads must share the same plane normal k";
}

// -----------------------------------------------------------------------------
// calculate_inner_coord_coefficients()
// -----------------------------------------------------------------------------
void fss_test_calculate_inner_coord_coefficients()
{
    std::cerr << "\n[TEST] fss_test_calculate_inner_coord_coefficients\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    calculate_inner_coord_coefficients(TARG, COM, crdset)\n"
              << "  Criteria:    decomposes (TARG-COM) into alpha*i + beta*j + gama*k.\n";

    // The identity triad exercises the axis-aligned special branches exactly.
    std::array<double, 9> ident { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    Coord com { 1.0, 2.0, 3.0 };

    // TARG == COM: the source returns all-zero coefficients immediately.
    std::cerr << "  -> TARG == COM returns (0,0,0)\n";
    std::array<double, 3> zero = calculate_inner_coord_coefficients(com, com, ident);
    EXPECT_DOUBLE_EQ(zero[0], 0.0) << "alpha must be 0 when TARG == COM";
    EXPECT_DOUBLE_EQ(zero[1], 0.0) << "beta must be 0 when TARG == COM";
    EXPECT_DOUBLE_EQ(zero[2], 0.0) << "gama must be 0 when TARG == COM";

    // Pure +i offset -> alpha only.
    std::cerr << "  -> offset (5,0,0) along i gives alpha=5, beta=0, gama=0\n";
    std::array<double, 3> alongI
        = calculate_inner_coord_coefficients(Coord { com.x + 5.0, com.y, com.z }, com, ident);
    EXPECT_NEAR(alongI[0], 5.0, 1e-12) << "alpha must be 5";
    EXPECT_NEAR(alongI[1], 0.0, 1e-12) << "beta must be 0";
    EXPECT_NEAR(alongI[2], 0.0, 1e-12) << "gama must be 0";

    // Pure -j offset -> negative beta only (the source restores the sign).
    std::cerr << "  -> offset (0,-7,0) along -j gives alpha=0, beta=-7, gama=0\n";
    std::array<double, 3> alongJ
        = calculate_inner_coord_coefficients(Coord { com.x, com.y - 7.0, com.z }, com, ident);
    EXPECT_NEAR(alongJ[0], 0.0, 1e-12) << "alpha must be 0";
    EXPECT_NEAR(alongJ[1], -7.0, 1e-12) << "beta must be -7";
    EXPECT_NEAR(alongJ[2], 0.0, 1e-12) << "gama must be 0";

    // Pure +k offset -> gama only.
    std::cerr << "  -> offset (0,0,3) along k gives alpha=0, beta=0, gama=3\n";
    std::array<double, 3> alongK
        = calculate_inner_coord_coefficients(Coord { com.x, com.y, com.z + 3.0 }, com, ident);
    EXPECT_NEAR(alongK[0], 0.0, 1e-12) << "alpha must be 0";
    EXPECT_NEAR(alongK[1], 0.0, 1e-12) << "beta must be 0";
    EXPECT_NEAR(alongK[2], 3.0, 1e-12) << "gama must be 3";

    // --- General (non-degenerate) branch with a generic triad, i.e. one with
    //     no zero components, produced by inner_coord_set().
    std::cerr << "  -> general branch with a generic (no zero component) triad\n";
    const double R = 10.0;
    Coord comG = fss_on_sphere(1.0, 1.0, 1.0, R);
    Coord comnewG = fss_on_sphere(1.05, 1.1, 0.95, R);
    std::array<double, 9> crdset = inner_coord_set(comG, comnewG);

    const Coord bi = fss_basis(crdset, 0);
    const Coord bj = fss_basis(crdset, 1);
    const Coord bk = fss_basis(crdset, 2);

    // Choose the target so that no projection lands on the degenerate branches.
    const double alphaIn = 1.5, betaIn = -2.25, gamaIn = 3.75;
    Coord targRel { alphaIn * bi.x + betaIn * bj.x + gamaIn * bk.x,
        alphaIn * bi.y + betaIn * bj.y + gamaIn * bk.y,
        alphaIn * bi.z + betaIn * bj.z + gamaIn * bk.z };
    Coord targ { comG.x + targRel.x, comG.y + targRel.y, comG.z + targRel.z };

    std::array<double, 3> M = calculate_inner_coord_coefficients(targ, comG, crdset);
    std::cerr << "     recovered (alpha,beta,gama) = (" << M[0] << ", " << M[1] << ", " << M[2] << ")\n";
    EXPECT_NEAR(M[0], alphaIn, 1e-8) << "alpha must be recovered";
    EXPECT_NEAR(M[1], betaIn, 1e-8) << "beta must be recovered";
    EXPECT_NEAR(M[2], gamaIn, 1e-8) << "gama must be recovered";

    // Reconstructing from the coefficients must return the original offset.
    Coord rebuilt { M[0] * bi.x + M[1] * bj.x + M[2] * bk.x,
        M[0] * bi.y + M[1] * bj.y + M[2] * bk.y,
        M[0] * bi.z + M[1] * bj.z + M[2] * bk.z };
    EXPECT_NEAR(rebuilt.x, targRel.x, 1e-8) << "reconstruction must restore the offset (x)";
    EXPECT_NEAR(rebuilt.y, targRel.y, 1e-8) << "reconstruction must restore the offset (y)";
    EXPECT_NEAR(rebuilt.z, targRel.z, 1e-8) << "reconstruction must restore the offset (z)";
}

// -----------------------------------------------------------------------------
// translate_on_sphere()
// -----------------------------------------------------------------------------
void fss_test_translate_on_sphere()
{
    std::cerr << "\n[TEST] fss_test_translate_on_sphere\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    translate_on_sphere(targ, COM, COMnew, crdset, crdsetnew)\n"
              << "  Criteria:    a rigid rotation about the plane normal - distances to\n"
              << "               the complex center and to the sphere origin are conserved,\n"
              << "               and a target sitting on COM lands exactly on COMnew.\n";

    const double R = 10.0;
    Coord com = fss_on_sphere(1.0, 1.0, 1.0, R);
    Coord comnew = fss_on_sphere(1.05, 1.1, 0.95, R);

    std::array<double, 9> crdset = inner_coord_set(com, comnew);
    std::array<double, 9> crdsetnew = inner_coord_set_new(com, comnew);

    // --- No motion at all: the source returns targ untouched.
    std::cerr << "  -> COM == COMnew short circuits and returns the target unchanged\n";
    Coord staticTarg { 3.0, -4.0, 5.0 };
    Coord unchanged = translate_on_sphere(staticTarg, com, com, crdset, crdsetnew);
    EXPECT_DOUBLE_EQ(unchanged.x, staticTarg.x) << "x must be untouched when COM does not move";
    EXPECT_DOUBLE_EQ(unchanged.y, staticTarg.y) << "y must be untouched when COM does not move";
    EXPECT_DOUBLE_EQ(unchanged.z, staticTarg.z) << "z must be untouched when COM does not move";

    // --- A target sitting exactly on the complex center follows it exactly.
    std::cerr << "  -> a target sitting on COM must be mapped onto COMnew\n";
    Coord center = translate_on_sphere(com, com, comnew, crdset, crdsetnew);
    EXPECT_NEAR(center.x, comnew.x, 1e-9) << "the center must move to COMnew (x)";
    EXPECT_NEAR(center.y, comnew.y, 1e-9) << "the center must move to COMnew (y)";
    EXPECT_NEAR(center.z, comnew.z, 1e-9) << "the center must move to COMnew (z)";

    // --- A generic member of the complex, also on the sphere surface.
    std::cerr << "  -> a generic surface target keeps its radius and its distance to COM\n";
    Coord targ = fss_on_sphere(0.2, 1.0, 0.5, R);
    Coord moved = translate_on_sphere(targ, com, comnew, crdset, crdsetnew);
    std::cerr << "     targ  = (" << targ.x << ", " << targ.y << ", " << targ.z << ")\n";
    std::cerr << "     moved = (" << moved.x << ", " << moved.y << ", " << moved.z << ")\n";

    const double dOld = fss_mag(Coord { targ.x - com.x, targ.y - com.y, targ.z - com.z });
    const double dNew = fss_mag(Coord { moved.x - comnew.x, moved.y - comnew.y, moved.z - comnew.z });
    EXPECT_NEAR(dNew, dOld, 1e-8) << "the distance to the complex center must be conserved";
    EXPECT_NEAR(fss_mag(moved), R, 1e-8) << "a surface target must stay on the sphere";
    EXPECT_NEAR(fss_dot(moved, comnew), fss_dot(targ, com), 1e-6)
        << "the angle between the target and the complex center must be conserved";

    // The point really did move (COM moved, so the target must too).
    EXPECT_GT(fss_mag(Coord { moved.x - targ.x, moved.y - targ.y, moved.z - targ.z }), 1e-6)
        << "a non-zero COM displacement must displace the target";
}

// -----------------------------------------------------------------------------
// rotate_on_sphere()
// -----------------------------------------------------------------------------
void fss_test_rotate_on_sphere()
{
    std::cerr << "\n[TEST] fss_test_rotate_on_sphere\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    rotate_on_sphere(Targ, COM, crdset, dangle)\n"
              << "  Criteria:    rotation about the O-COM axis (the i axis of crdset):\n"
              << "               the i component and the distance to COM are conserved,\n"
              << "               dangle == 0 and dangle == 2*PI are no-ops, and the j/k\n"
              << "               components follow the analytic rotation formula.\n";

    const double R = 10.0;
    Coord com = fss_on_sphere(1.0, 1.0, 1.0, R);

    // Use the "no motion" flavour of the triad: i is along COM, j/k span the
    // tangent plane. This is exactly what the caller uses for pure rotation.
    std::array<double, 9> crdset = inner_coord_set(com, com);
    const Coord bi = fss_basis(crdset, 0);
    const Coord bj = fss_basis(crdset, 1);
    const Coord bk = fss_basis(crdset, 2);

    // Build a target whose axial component is 0.5 (deliberately NOT 1.0, which
    // would hit the early-out branch coded in the source).
    const double a = 0.5, b = 2.0, c = 3.0;
    Coord targ { com.x + a * bi.x + b * bj.x + c * bk.x,
        com.y + a * bi.y + b * bj.y + c * bk.y,
        com.z + a * bi.z + b * bj.z + c * bk.z };

    // --- dangle = 0 must be a no-op.
    std::cerr << "  -> dangle = 0 must leave the target where it is\n";
    Coord noRot = rotate_on_sphere(targ, com, crdset, 0.0);
    EXPECT_NEAR(noRot.x, targ.x, 1e-9) << "zero rotation must not move x";
    EXPECT_NEAR(noRot.y, targ.y, 1e-9) << "zero rotation must not move y";
    EXPECT_NEAR(noRot.z, targ.z, 1e-9) << "zero rotation must not move z";

    // --- A finite rotation must match the analytic result.
    const double dangle = 0.7;
    std::cerr << "  -> dangle = " << dangle << " must match the analytic j/k rotation\n";
    Coord rot = rotate_on_sphere(targ, com, crdset, dangle);

    const double m = std::sqrt(b * b + c * c);
    const double phi0 = std::atan2(c, b);
    const double bExp = m * std::cos(phi0 + dangle);
    const double cExp = m * std::sin(phi0 + dangle);
    Coord expected { com.x + a * bi.x + bExp * bj.x + cExp * bk.x,
        com.y + a * bi.y + bExp * bj.y + cExp * bk.y,
        com.z + a * bi.z + bExp * bj.z + cExp * bk.z };

    std::cerr << "     rotated  = (" << rot.x << ", " << rot.y << ", " << rot.z << ")\n";
    std::cerr << "     expected = (" << expected.x << ", " << expected.y << ", " << expected.z << ")\n";
    EXPECT_NEAR(rot.x, expected.x, 1e-8) << "rotated x must match the analytic value";
    EXPECT_NEAR(rot.y, expected.y, 1e-8) << "rotated y must match the analytic value";
    EXPECT_NEAR(rot.z, expected.z, 1e-8) << "rotated z must match the analytic value";

    // --- Invariants of a rotation about the i axis.
    const Coord relOld { targ.x - com.x, targ.y - com.y, targ.z - com.z };
    const Coord relNew { rot.x - com.x, rot.y - com.y, rot.z - com.z };
    EXPECT_NEAR(fss_mag(relNew), fss_mag(relOld), 1e-9)
        << "the distance to COM must be conserved by a rotation";
    EXPECT_NEAR(fss_dot(relNew, bi), fss_dot(relOld, bi), 1e-9)
        << "the component along the rotation axis must be conserved";

    // --- Rotating forward then backward returns to the start.
    std::cerr << "  -> rotating by +dangle then -dangle must return the original point\n";
    Coord roundTrip = rotate_on_sphere(rot, com, crdset, -dangle);
    EXPECT_NEAR(roundTrip.x, targ.x, 1e-8) << "round trip must restore x";
    EXPECT_NEAR(roundTrip.y, targ.y, 1e-8) << "round trip must restore y";
    EXPECT_NEAR(roundTrip.z, targ.z, 1e-8) << "round trip must restore z";

    // --- A target sitting exactly on COM has no tangential component and is
    //     returned unchanged by the early-out branch.
    std::cerr << "  -> a target coincident with COM is returned unchanged\n";
    Coord onCom = rotate_on_sphere(com, com, crdset, 1.234);
    EXPECT_DOUBLE_EQ(onCom.x, com.x) << "a target on COM must not move (x)";
    EXPECT_DOUBLE_EQ(onCom.y, com.y) << "a target on COM must not move (y)";
    EXPECT_DOUBLE_EQ(onCom.z, com.z) << "a target on COM must not move (z)";
}

// -----------------------------------------------------------------------------
// calc_bindRadius2D()
// -----------------------------------------------------------------------------
void fss_test_calc_bindRadius2D()
{
    std::cerr << "\n[TEST] fss_test_calc_bindRadius2D\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    calc_bindRadius2D(bindRadius, iFace)\n"
              << "  Criteria:    converts a chord of length bindRadius into the\n"
              << "               corresponding great-circle arc, 2*R*asin(sigma/(2R)).\n";

    // Explicit formula check on a sphere of radius 10.
    std::cerr << "  -> R = 10, sigma = 1: expect 20*asin(0.05)\n";
    Coord iface { 0.0, 0.0, 10.0 };
    const double expected = 2.0 * 10.0 * std::asin(0.5 * 1.0 / 10.0);
    EXPECT_NEAR(calc_bindRadius2D(1.0, iface), expected, 1e-12)
        << "must equal 2*R*asin(sigma/(2R))";

    // The arc is always at least as long as the chord it subtends.
    std::cerr << "  -> the arc must never be shorter than the chord\n";
    EXPECT_GE(calc_bindRadius2D(1.0, iface), 1.0) << "arc length >= chord length";

    // On a much larger sphere the arc converges to the chord.
    std::cerr << "  -> on a very large sphere (R = 1e6) the arc ~ the chord\n";
    Coord bigIface { 1.0e6, 0.0, 0.0 };
    EXPECT_NEAR(calc_bindRadius2D(1.0, bigIface), 1.0, 1e-6)
        << "for R >> sigma the arc must converge to the chord";

    // A smaller sphere curves more, so the arc grows.
    std::cerr << "  -> a smaller sphere must give a longer arc for the same chord\n";
    Coord smallIface { 0.0, 2.0, 0.0 };
    EXPECT_GT(calc_bindRadius2D(1.0, smallIface), calc_bindRadius2D(1.0, iface))
        << "a tighter curvature must stretch the arc";
}

// -----------------------------------------------------------------------------
// set_memProtein_sphere() -- explicit lipid case
// -----------------------------------------------------------------------------
void fss_test_set_memProtein_sphere_explicit()
{
    std::cerr << "\n[TEST] fss_test_set_memProtein_sphere_explicit\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    set_memProtein_sphere(reactCom, memProtein, moleculeList, membrane)\n"
              << "  Scenario:    an EXPLICIT lipid is a member of the complex.\n"
              << "  Criteria:    memProtein is the lipid, its comCoord is taken from the\n"
              << "               temporary association coordinates, and its interface is\n"
              << "               pushed radially inward by the COM-interface bond length.\n";

    const double R = 10.0;
    const double bond = 0.5;

    Membrane membraneObject;
    membraneObject.implicitLipid = false;

    // Molecule 0: an ordinary protein (must be ignored by the search).
    // Molecule 1: the explicit lipid sitting on the sphere surface.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(fss_make_protein(Coord { 0.0, 0.0, R + 1.0 }, Coord { 0.0, 0.0, R + 0.5 }));
    moleculeList.push_back(fss_make_lipid(Coord { 0.0, 0.0, R }, Coord { 0.0, 0.0, R - bond }));
    moleculeList[0].index = 0;
    moleculeList[1].index = 1;

    Complex reactCom;
    reactCom.memberList = { 0, 1 };

    Molecule memProtein;
    set_memProtein_sphere(reactCom, memProtein, moleculeList, membraneObject);

    std::cerr << "  memProtein.comCoord = (" << memProtein.comCoord.x << ", " << memProtein.comCoord.y
              << ", " << memProtein.comCoord.z << ")\n";
    std::cerr << "  memProtein.iface[0] = (" << memProtein.interfaceList[0].coord.x << ", "
              << memProtein.interfaceList[0].coord.y << ", " << memProtein.interfaceList[0].coord.z << ")\n";

    // The lipid, not the protein, must have been selected.
    EXPECT_TRUE(memProtein.isLipid) << "the selected molecule must be the explicit lipid";

    // comCoord is copied from tmpComCoord.
    EXPECT_NEAR(memProtein.comCoord.x, 0.0, 1e-12) << "comCoord.x must come from tmpComCoord";
    EXPECT_NEAR(memProtein.comCoord.y, 0.0, 1e-12) << "comCoord.y must come from tmpComCoord";
    EXPECT_NEAR(memProtein.comCoord.z, R, 1e-12) << "comCoord.z must come from tmpComCoord";

    // The interface is re-seated at (R - bond)/R * comCoord.
    EXPECT_NEAR(memProtein.interfaceList[0].coord.x, 0.0, 1e-12) << "interface x must be 0";
    EXPECT_NEAR(memProtein.interfaceList[0].coord.y, 0.0, 1e-12) << "interface y must be 0";
    EXPECT_NEAR(memProtein.interfaceList[0].coord.z, R - bond, 1e-12)
        << "interface must be pushed inward by the bond length";
}

// -----------------------------------------------------------------------------
// set_memProtein_sphere() -- implicit lipid case
// -----------------------------------------------------------------------------
void fss_test_set_memProtein_sphere_implicit()
{
    std::cerr << "\n[TEST] fss_test_set_memProtein_sphere_implicit\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    set_memProtein_sphere(...)\n"
              << "  Scenario:    IMPLICIT lipid model - the complex holds only a protein\n"
              << "               that is bound to an implicit lipid outside the member list.\n"
              << "  Criteria:    memProtein becomes the implicit lipid, its comCoord is the\n"
              << "               protein's bound interface, and its own interface is placed\n"
              << "               one bond length radially inward.\n";

    const double R = 10.0;
    const double bond = 0.5;

    Membrane membraneObject;
    membraneObject.implicitLipid = true;

    // Molecule 0: the protein, bound through interface 0 to molecule 1.
    Molecule protein = fss_make_protein(Coord { 0.0, 0.0, R + bond }, Coord { 0.0, 0.0, R });
    protein.index = 0;
    protein.interfaceList[0].isBound = true;
    protein.interfaceList[0].interaction.partnerIndex = 1;

    // Molecule 1: the implicit lipid; it needs one interface slot to be filled in.
    Molecule implicitLipid;
    implicitLipid.index = 1;
    implicitLipid.isImplicitLipid = true;
    implicitLipid.interfaceList.resize(1);

    std::vector<Molecule> moleculeList { protein, implicitLipid };

    Complex reactCom;
    reactCom.memberList = { 0 }; // the implicit lipid is NOT a member

    Molecule memProtein;
    set_memProtein_sphere(reactCom, memProtein, moleculeList, membraneObject);

    std::cerr << "  memProtein.comCoord = (" << memProtein.comCoord.x << ", " << memProtein.comCoord.y
              << ", " << memProtein.comCoord.z << ")\n";
    std::cerr << "  memProtein.iface[0] = (" << memProtein.interfaceList[0].coord.x << ", "
              << memProtein.interfaceList[0].coord.y << ", " << memProtein.interfaceList[0].coord.z << ")\n";

    EXPECT_TRUE(memProtein.isImplicitLipid) << "the implicit lipid partner must be selected";

    // The implicit lipid's COM is placed on the protein's bound interface.
    EXPECT_NEAR(memProtein.comCoord.x, 0.0, 1e-12) << "comCoord.x must be the bound interface x";
    EXPECT_NEAR(memProtein.comCoord.y, 0.0, 1e-12) << "comCoord.y must be the bound interface y";
    EXPECT_NEAR(memProtein.comCoord.z, R, 1e-12) << "comCoord.z must be the bound interface z";

    // Its interface is (R - bond)/R * comCoord with bond = |iface - com|.
    EXPECT_NEAR(memProtein.interfaceList[0].coord.z, R - bond, 1e-12)
        << "the implicit lipid interface must sit one bond length inward";
}

// -----------------------------------------------------------------------------
// find_Lipid_sphere() -- explicit lipid case
// -----------------------------------------------------------------------------
void fss_test_find_Lipid_sphere_explicit()
{
    std::cerr << "\n[TEST] fss_test_find_Lipid_sphere_explicit\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    find_Lipid_sphere(reactCom, Lipid, moleculeList, membrane)\n"
              << "  Scenario:    two EXPLICIT lipids in the complex at different radii.\n"
              << "  Criteria:    the lipid with the largest |tmpComCoord| is returned,\n"
              << "               copied verbatim (no coordinate rewriting).\n";

    Membrane membraneObject;
    membraneObject.implicitLipid = false;

    // Molecule 0: a protein (ignored).
    // Molecule 1: a lipid at radius 9   (should lose the "largest radius" race).
    // Molecule 2: a lipid at radius 10  (should be selected).
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(fss_make_protein(Coord { 1.0, 1.0, 1.0 }, Coord { 1.0, 1.0, 1.5 }));
    moleculeList.push_back(fss_make_lipid(Coord { 0.0, 0.0, 9.0 }, Coord { 0.0, 0.0, 8.5 }));
    moleculeList.push_back(fss_make_lipid(Coord { 0.0, 0.0, 10.0 }, Coord { 0.0, 0.0, 9.5 }));
    for (int i = 0; i < 3; ++i)
        moleculeList[i].index = i;

    Complex reactCom;
    reactCom.memberList = { 0, 1, 2 };

    Molecule lipid;
    find_Lipid_sphere(reactCom, lipid, moleculeList, membraneObject);

    std::cerr << "  Lipid.tmpComCoord = (" << lipid.tmpComCoord.x << ", " << lipid.tmpComCoord.y << ", "
              << lipid.tmpComCoord.z << ")\n";

    EXPECT_TRUE(lipid.isLipid) << "an explicit lipid must be returned";
    EXPECT_EQ(lipid.index, 2) << "the lipid furthest from the origin must win the search";

    // The explicit branch copies the molecule verbatim - coordinates untouched.
    EXPECT_NEAR(lipid.tmpComCoord.z, 10.0, 1e-12) << "tmpComCoord must be copied unchanged";
    EXPECT_NEAR(lipid.comCoord.z, 10.0, 1e-12) << "comCoord must be copied unchanged";
    ASSERT_FALSE(lipid.interfaceList.empty()) << "the copied lipid must keep its interface list";
    EXPECT_NEAR(lipid.interfaceList[0].coord.z, 9.5, 1e-12) << "interface coord must be copied unchanged";
}

// -----------------------------------------------------------------------------
// find_Lipid_sphere() -- implicit lipid case
// -----------------------------------------------------------------------------
void fss_test_find_Lipid_sphere_implicit()
{
    std::cerr << "\n[TEST] fss_test_find_Lipid_sphere_implicit\n"
              << "  Source file: functions_for_spherical_system.cpp\n"
              << "  Function:    find_Lipid_sphere(...)\n"
              << "  Scenario:    IMPLICIT lipid model - the complex only holds a protein\n"
              << "               bound to an implicit lipid.\n"
              << "  Criteria:    the implicit lipid is returned with com/tmpCom set to the\n"
              << "               protein's bound interface and its own interface set to the\n"
              << "               protein's COM (both real and temporary coordinates).\n";

    const double R = 10.0;
    const double bond = 0.5;

    Membrane membraneObject;
    membraneObject.implicitLipid = true;

    // Protein bound through its single interface to the implicit lipid (index 1).
    Molecule protein = fss_make_protein(Coord { 0.0, 0.0, R + bond }, Coord { 0.0, 0.0, R });
    protein.index = 0;
    protein.interfaceList[0].isBound = true;
    protein.interfaceList[0].interaction.partnerIndex = 1;

    // The implicit lipid: one interface, and an EMPTY tmpICoords so that
    // set_tmp_association_coords() (called inside the function) produces
    // exactly one temporary interface coordinate.
    Molecule implicitLipid;
    implicitLipid.index = 1;
    implicitLipid.isImplicitLipid = true;
    implicitLipid.interfaceList.resize(1);
    implicitLipid.tmpICoords.clear();

    std::vector<Molecule> moleculeList { protein, implicitLipid };

    Complex reactCom;
    reactCom.memberList = { 0 };

    Molecule lipid;
    find_Lipid_sphere(reactCom, lipid, moleculeList, membraneObject);

    std::cerr << "  Lipid.comCoord    = (" << lipid.comCoord.x << ", " << lipid.comCoord.y << ", "
              << lipid.comCoord.z << ")\n";
    std::cerr << "  Lipid.tmpComCoord = (" << lipid.tmpComCoord.x << ", " << lipid.tmpComCoord.y << ", "
              << lipid.tmpComCoord.z << ")\n";

    EXPECT_TRUE(lipid.isImplicitLipid) << "the implicit lipid partner must be returned";

    // The implicit lipid is anchored on the protein's bound interface.
    EXPECT_NEAR(lipid.comCoord.z, R, 1e-12) << "comCoord must be the protein's bound interface";
    EXPECT_NEAR(lipid.tmpComCoord.z, R, 1e-12) << "tmpComCoord must be the protein's bound interface";

    // Its interface is set to the protein's temporary COM.
    ASSERT_FALSE(lipid.interfaceList.empty()) << "the implicit lipid must keep an interface slot";
    EXPECT_NEAR(lipid.interfaceList[0].coord.z, R + bond, 1e-12)
        << "the implicit lipid interface must be the protein COM";

    // set_tmp_association_coords() must have created exactly one temporary
    // interface coordinate, which is then overwritten with the protein COM.
    ASSERT_EQ(lipid.tmpICoords.size(), 1u) << "exactly one temporary interface coordinate is expected";
    EXPECT_NEAR(lipid.tmpICoords[0].z, R + bond, 1e-12)
        << "the temporary interface coordinate must also be the protein COM";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each named helper is run inside its own TEST so that a
// failure in one does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(FunctionsForSphericalSystem, Radius) { fss_test_radius(); }
TEST(FunctionsForSphericalSystem, FindSphericalCoords) { fss_test_find_spherical_coords(); }
TEST(FunctionsForSphericalSystem, FindCardesianCoords) { fss_test_find_cardesian_coords(); }
TEST(FunctionsForSphericalSystem, ThetaPlus) { fss_test_theta_plus(); }
TEST(FunctionsForSphericalSystem, PhiPlus) { fss_test_phi_plus(); }
TEST(FunctionsForSphericalSystem, AnglePlus) { fss_test_angle_plus(); }
TEST(FunctionsForSphericalSystem, FindPositionAfterAssociation) { fss_test_find_position_after_association(); }
TEST(FunctionsForSphericalSystem, InnerCoordSet) { fss_test_inner_coord_set(); }
TEST(FunctionsForSphericalSystem, InnerCoordSetNew) { fss_test_inner_coord_set_new(); }
TEST(FunctionsForSphericalSystem, CalculateInnerCoordCoefficients) { fss_test_calculate_inner_coord_coefficients(); }
TEST(FunctionsForSphericalSystem, TranslateOnSphere) { fss_test_translate_on_sphere(); }
TEST(FunctionsForSphericalSystem, RotateOnSphere) { fss_test_rotate_on_sphere(); }
TEST(FunctionsForSphericalSystem, CalcBindRadius2D) { fss_test_calc_bindRadius2D(); }
TEST(FunctionsForSphericalSystem, SetMemProteinSphereExplicit) { fss_test_set_memProtein_sphere_explicit(); }
TEST(FunctionsForSphericalSystem, SetMemProteinSphereImplicit) { fss_test_set_memProtein_sphere_implicit(); }
TEST(FunctionsForSphericalSystem, FindLipidSphereExplicit) { fss_test_find_Lipid_sphere_explicit(); }
TEST(FunctionsForSphericalSystem, FindLipidSphereImplicit) { fss_test_find_Lipid_sphere_implicit(); }