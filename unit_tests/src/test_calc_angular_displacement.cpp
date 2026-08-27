/*! \file test_calc_angular_displacement.cpp
 *
 * ### Unit test for src/reactions/calc_angular_displacement.cpp
 *
 * Function under test:
 *
 *     void calc_angular_displacement(int ifaceIndex1, int ifaceIndex2,
 *                                    Molecule& reactMol1, Molecule& reactMol2,
 *                                    Complex& reactCom1, Complex& reactCom2,
 *                                    std::vector<Molecule>& moleculeList);
 *
 * IMPORTANT NOTE ABOUT WHAT CAN BE ASSERTED
 * -----------------------------------------
 * In the current source every `std::cout` statement inside
 * calc_angular_displacement() is commented out, the computed angle is stored in
 * a local variable (`currTheta`) that is never returned, and none of the
 * arguments are written to.  The function therefore has **no observable output**
 * other than:
 *
 *   1. it must not crash / must not read out of bounds, and
 *   2. it must leave every argument it is handed completely unmodified
 *      (it is a pure "diagnostic" routine),
 *   3. it must take its early-return path when the interface-to-COM vector has
 *      (near) zero magnitude, i.e. for POINT particles.
 *
 * So the tests below:
 *   * build fully-initialised Molecules/Complexes (an under-initialised
 *     Molecule would segfault the whole suite),
 *   * snapshot every coordinate the function touches, call the function, and
 *     verify the snapshot is bit-identical afterwards,
 *   * independently re-compute, with Vector::calc_magnitude()/Vector::dot_theta(),
 *     exactly the four angles the function computes internally, so the geometry
 *     the routine inspects is documented and checked, and
 *   * exercise each of the four early-return (POINT particle) branches so all
 *     code paths in the file are executed.
 *
 * Verbose progress information is written to stderr so the reader can follow
 * which source file / function is under test and what each assertion checks.
 */

#include "reactions/association/association.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vector.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers.  All names are prefixed with "cad_" (calc_angular_displacement)
// so they cannot collide with helpers from other files in the test suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully-initialised Molecule with a given COM and interfaces.
 *
 * The temporary association coordinates (tmpComCoord / tmpICoords) are seeded
 * with the real coordinates via Molecule::set_tmp_association_coords(), which is
 * exactly what association() does before it starts moving things around.
 *
 * \param[in] com          center of mass of the molecule
 * \param[in] ifaceCoords  absolute coordinates of each interface
 */
Molecule cad_make_molecule(const Coord& com, const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.index = 0;
    mol.myComIndex = 0;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.comCoord = com;

    mol.interfaceList.clear();
    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface;
        iface.coord = ifaceCoords[i];
        iface.relIndex = static_cast<int>(i);
        iface.index = static_cast<int>(i);
        iface.molTypeIndex = 0;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
    }

    // Seed tmpComCoord and tmpICoords so both are the same size as
    // interfaceList; calc_angular_displacement indexes tmpICoords directly.
    mol.set_tmp_association_coords();
    return mol;
}

/*! \brief Build a minimal but fully-initialised Complex holding one molecule.
 *
 * calc_angular_displacement only reads comCoord and tmpComCoord, but we fill in
 * the other members so the object is safe to pass anywhere.
 */
Complex cad_make_complex(const Coord& com, int memberMolIndex)
{
    Complex com1;
    com1.index = 0;
    com1.comCoord = com;
    com1.tmpComCoord = com; // association starts with tmp == real
    com1.mass = 1.0;
    com1.radius = 1.0;
    com1.D = Coord{ 1.0, 1.0, 1.0 };
    com1.Dr = Coord{ 0.01, 0.01, 0.01 };
    com1.memberList.clear();
    com1.memberList.push_back(memberMolIndex);
    com1.numEachMol.clear();
    com1.numEachMol.push_back(1);
    com1.isEmpty = false;
    com1.OnSurface = false;
    return com1;
}

/*! \brief Rotate a molecule's *temporary* coordinates about the z axis.
 *
 * Rotating the tmp coordinates around the molecule's own tmpComCoord leaves the
 * COM in place and rotates the interface-to-COM vector by exactly `angle`,
 * which is the situation calc_angular_displacement is meant to measure.
 */
void cad_rotate_tmp_about_z(Molecule& mol, const Coord& pivot, double angle)
{
    const double ca = std::cos(angle);
    const double sa = std::sin(angle);

    auto rotate_one = [&](Coord& c) {
        const double dx = c.x - pivot.x;
        const double dy = c.y - pivot.y;
        c.x = pivot.x + dx * ca - dy * sa;
        c.y = pivot.y + dx * sa + dy * ca;
    };

    rotate_one(mol.tmpComCoord);
    for (auto& ic : mol.tmpICoords)
        rotate_one(ic);
}

/*! \brief Angle (radians) between (a - b) and (c - d), using the same Vector
 *         calls the function under test uses internally.
 *
 * This mirrors, line for line, what calc_angular_displacement computes so the
 * test documents the geometry the routine is inspecting.
 */
double cad_angle_between(const Coord& tipTmp, const Coord& baseTmp,
    const Coord& tipOrig, const Coord& baseOrig)
{
    Vector vTmp{ tipTmp - baseTmp };
    Vector vOrig{ tipOrig - baseOrig };
    vTmp.calc_magnitude();
    vOrig.calc_magnitude();
    return vOrig.dot_theta(vTmp);
}

/*! \brief Small container used to snapshot everything the function reads. */
struct CadSnapshot {
    Coord mol1Com, mol1Tmp, mol1Iface, mol1TmpIface;
    Coord mol2Com, mol2Tmp, mol2Iface, mol2TmpIface;
    Coord com1Com, com1Tmp, com2Com, com2Tmp;
};

/*! \brief Take a snapshot of all coordinates calc_angular_displacement reads. */
CadSnapshot cad_snapshot(const Molecule& m1, const Molecule& m2, const Complex& c1,
    const Complex& c2, int i1, int i2)
{
    CadSnapshot s;
    s.mol1Com = m1.comCoord;
    s.mol1Tmp = m1.tmpComCoord;
    s.mol1Iface = m1.interfaceList[i1].coord;
    s.mol1TmpIface = m1.tmpICoords[i1];
    s.mol2Com = m2.comCoord;
    s.mol2Tmp = m2.tmpComCoord;
    s.mol2Iface = m2.interfaceList[i2].coord;
    s.mol2TmpIface = m2.tmpICoords[i2];
    s.com1Com = c1.comCoord;
    s.com1Tmp = c1.tmpComCoord;
    s.com2Com = c2.comCoord;
    s.com2Tmp = c2.tmpComCoord;
    return s;
}

/*! \brief Assert that nothing the function reads has been written to.
 *
 * Coord::operator== compares components after roundv() (4 decimal-place
 * truncation with a sign-dependent rule), which is more than enough resolution
 * for "did anything change at all"; we additionally compare the raw doubles
 * exactly, because the function is expected to be a strict read-only routine.
 */
void cad_expect_unchanged(const CadSnapshot& before, const Molecule& m1,
    const Molecule& m2, const Complex& c1, const Complex& c2, int i1, int i2)
{
    EXPECT_DOUBLE_EQ(m1.comCoord.x, before.mol1Com.x) << "reactMol1.comCoord.x must not change";
    EXPECT_DOUBLE_EQ(m1.comCoord.y, before.mol1Com.y) << "reactMol1.comCoord.y must not change";
    EXPECT_DOUBLE_EQ(m1.comCoord.z, before.mol1Com.z) << "reactMol1.comCoord.z must not change";

    EXPECT_DOUBLE_EQ(m1.tmpComCoord.x, before.mol1Tmp.x) << "reactMol1.tmpComCoord.x must not change";
    EXPECT_DOUBLE_EQ(m1.tmpComCoord.y, before.mol1Tmp.y) << "reactMol1.tmpComCoord.y must not change";
    EXPECT_DOUBLE_EQ(m1.tmpComCoord.z, before.mol1Tmp.z) << "reactMol1.tmpComCoord.z must not change";

    EXPECT_DOUBLE_EQ(m1.interfaceList[i1].coord.x, before.mol1Iface.x) << "mol1 iface coord must not change";
    EXPECT_DOUBLE_EQ(m1.interfaceList[i1].coord.y, before.mol1Iface.y) << "mol1 iface coord must not change";
    EXPECT_DOUBLE_EQ(m1.interfaceList[i1].coord.z, before.mol1Iface.z) << "mol1 iface coord must not change";

    EXPECT_DOUBLE_EQ(m1.tmpICoords[i1].x, before.mol1TmpIface.x) << "mol1 tmp iface coord must not change";
    EXPECT_DOUBLE_EQ(m1.tmpICoords[i1].y, before.mol1TmpIface.y) << "mol1 tmp iface coord must not change";
    EXPECT_DOUBLE_EQ(m1.tmpICoords[i1].z, before.mol1TmpIface.z) << "mol1 tmp iface coord must not change";

    EXPECT_DOUBLE_EQ(m2.comCoord.x, before.mol2Com.x) << "reactMol2.comCoord.x must not change";
    EXPECT_DOUBLE_EQ(m2.comCoord.y, before.mol2Com.y) << "reactMol2.comCoord.y must not change";
    EXPECT_DOUBLE_EQ(m2.comCoord.z, before.mol2Com.z) << "reactMol2.comCoord.z must not change";

    EXPECT_DOUBLE_EQ(m2.tmpComCoord.x, before.mol2Tmp.x) << "reactMol2.tmpComCoord.x must not change";
    EXPECT_DOUBLE_EQ(m2.tmpComCoord.y, before.mol2Tmp.y) << "reactMol2.tmpComCoord.y must not change";
    EXPECT_DOUBLE_EQ(m2.tmpComCoord.z, before.mol2Tmp.z) << "reactMol2.tmpComCoord.z must not change";

    EXPECT_DOUBLE_EQ(m2.interfaceList[i2].coord.x, before.mol2Iface.x) << "mol2 iface coord must not change";
    EXPECT_DOUBLE_EQ(m2.interfaceList[i2].coord.y, before.mol2Iface.y) << "mol2 iface coord must not change";
    EXPECT_DOUBLE_EQ(m2.interfaceList[i2].coord.z, before.mol2Iface.z) << "mol2 iface coord must not change";

    EXPECT_DOUBLE_EQ(m2.tmpICoords[i2].x, before.mol2TmpIface.x) << "mol2 tmp iface coord must not change";
    EXPECT_DOUBLE_EQ(m2.tmpICoords[i2].y, before.mol2TmpIface.y) << "mol2 tmp iface coord must not change";
    EXPECT_DOUBLE_EQ(m2.tmpICoords[i2].z, before.mol2TmpIface.z) << "mol2 tmp iface coord must not change";

    EXPECT_DOUBLE_EQ(c1.comCoord.x, before.com1Com.x) << "reactCom1.comCoord must not change";
    EXPECT_DOUBLE_EQ(c1.comCoord.y, before.com1Com.y) << "reactCom1.comCoord must not change";
    EXPECT_DOUBLE_EQ(c1.comCoord.z, before.com1Com.z) << "reactCom1.comCoord must not change";

    EXPECT_DOUBLE_EQ(c1.tmpComCoord.x, before.com1Tmp.x) << "reactCom1.tmpComCoord must not change";
    EXPECT_DOUBLE_EQ(c1.tmpComCoord.y, before.com1Tmp.y) << "reactCom1.tmpComCoord must not change";
    EXPECT_DOUBLE_EQ(c1.tmpComCoord.z, before.com1Tmp.z) << "reactCom1.tmpComCoord must not change";

    EXPECT_DOUBLE_EQ(c2.comCoord.x, before.com2Com.x) << "reactCom2.comCoord must not change";
    EXPECT_DOUBLE_EQ(c2.comCoord.y, before.com2Com.y) << "reactCom2.comCoord must not change";
    EXPECT_DOUBLE_EQ(c2.comCoord.z, before.com2Com.z) << "reactCom2.comCoord must not change";

    EXPECT_DOUBLE_EQ(c2.tmpComCoord.x, before.com2Tmp.x) << "reactCom2.tmpComCoord must not change";
    EXPECT_DOUBLE_EQ(c2.tmpComCoord.y, before.com2Tmp.y) << "reactCom2.tmpComCoord must not change";
    EXPECT_DOUBLE_EQ(c2.tmpComCoord.z, before.com2Tmp.z) << "reactCom2.tmpComCoord must not change";
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: no association move has been made yet (tmp coords == real coords).
//         All four internal angles must therefore be exactly 0, and the
//         function must not touch any of its arguments.
// -----------------------------------------------------------------------------
void cad_test_zero_displacement()
{
    std::cerr << "\n[TEST] cad_test_zero_displacement\n"
              << "  Source file: src/reactions/calc_angular_displacement.cpp\n"
              << "  Function:    calc_angular_displacement()\n"
              << "  Scenario:    tmp coordinates identical to the real coordinates\n"
              << "               (i.e. no association rotation has happened yet).\n"
              << "  Pass:        the four interface->COM angles the routine measures\n"
              << "               are all 0 rad, and no argument is modified.\n";

    // Molecule 1: COM at origin, single interface 1 nm along +x.
    Molecule mol1 = cad_make_molecule(Coord{ 0.0, 0.0, 0.0 }, { Coord{ 1.0, 0.0, 0.0 } });
    // Molecule 2: COM at (5,0,0), single interface at (4,0,0) (pointing back at mol1).
    Molecule mol2 = cad_make_molecule(Coord{ 5.0, 0.0, 0.0 }, { Coord{ 4.0, 0.0, 0.0 } });
    mol2.index = 1;
    mol2.myComIndex = 1;

    Complex com1 = cad_make_complex(mol1.comCoord, 0);
    Complex com2 = cad_make_complex(mol2.comCoord, 1);
    com2.index = 1;

    std::vector<Molecule> moleculeList{ mol1, mol2 };

    const int i1 = 0;
    const int i2 = 0;

    // Independently reproduce the four angles the function computes internally.
    const double a1 = cad_angle_between(mol1.tmpICoords[i1], mol1.tmpComCoord,
        mol1.interfaceList[i1].coord, mol1.comCoord);
    const double a2 = cad_angle_between(mol2.tmpICoords[i2], mol2.tmpComCoord,
        mol2.interfaceList[i2].coord, mol2.comCoord);
    const double a3 = cad_angle_between(mol1.tmpICoords[i1], com1.tmpComCoord,
        mol1.interfaceList[i1].coord, com1.comCoord);
    const double a4 = cad_angle_between(mol2.tmpICoords[i2], com2.tmpComCoord,
        mol2.interfaceList[i2].coord, com2.comCoord);

    std::cerr << "  Independently computed angles (rad): mol1->molCOM=" << a1
              << ", mol2->molCOM=" << a2 << ", mol1->complexCOM=" << a3
              << ", mol2->complexCOM=" << a4 << '\n';

    EXPECT_NEAR(a1, 0.0, 1e-12) << "mol1 interface->molecule-COM angle should be 0";
    EXPECT_NEAR(a2, 0.0, 1e-12) << "mol2 interface->molecule-COM angle should be 0";
    EXPECT_NEAR(a3, 0.0, 1e-12) << "mol1 interface->complex-COM angle should be 0";
    EXPECT_NEAR(a4, 0.0, 1e-12) << "mol2 interface->complex-COM angle should be 0";

    const CadSnapshot before = cad_snapshot(mol1, mol2, com1, com2, i1, i2);

    std::cerr << "  Calling calc_angular_displacement()...\n";
    calc_angular_displacement(i1, i2, mol1, mol2, com1, com2, moleculeList);
    std::cerr << "  Returned normally.\n";

    cad_expect_unchanged(before, mol1, mol2, com1, com2, i1, i2);
}

// -----------------------------------------------------------------------------
// Test 2: both molecules' temporary coordinates rotated by a known angle about
//         their own COM.  Documents that the routine is measuring exactly that
//         rotation, and again checks nothing is written back.
// -----------------------------------------------------------------------------
void cad_test_known_rotation()
{
    std::cerr << "\n[TEST] cad_test_known_rotation\n"
              << "  Source file: src/reactions/calc_angular_displacement.cpp\n"
              << "  Function:    calc_angular_displacement()\n"
              << "  Scenario:    mol1 tmp coords rotated +90 deg about z around its\n"
              << "               own COM, mol2 tmp coords rotated +180 deg.\n"
              << "  Pass:        the interface->COM angles equal pi/2 and pi, and\n"
              << "               the function leaves all arguments untouched.\n";

    Molecule mol1 = cad_make_molecule(Coord{ 0.0, 0.0, 0.0 }, { Coord{ 2.0, 0.0, 0.0 } });
    Molecule mol2 = cad_make_molecule(Coord{ 10.0, 0.0, 0.0 }, { Coord{ 12.0, 0.0, 0.0 } });
    mol2.index = 1;
    mol2.myComIndex = 1;

    // Rotate ONLY the temporary (association) coordinates about each molecule's
    // own tmp COM, so the molecule COM stays put and the iface vector turns.
    cad_rotate_tmp_about_z(mol1, mol1.tmpComCoord, M_PI / 2.0);
    cad_rotate_tmp_about_z(mol2, mol2.tmpComCoord, M_PI);

    // Complexes: put the complex COM at the molecule COM so the complex-level
    // angles are the same rotations as the molecule-level ones.
    Complex com1 = cad_make_complex(mol1.comCoord, 0);
    Complex com2 = cad_make_complex(mol2.comCoord, 1);
    com2.index = 1;

    std::vector<Molecule> moleculeList{ mol1, mol2 };

    const int i1 = 0;
    const int i2 = 0;

    const double a1 = cad_angle_between(mol1.tmpICoords[i1], mol1.tmpComCoord,
        mol1.interfaceList[i1].coord, mol1.comCoord);
    const double a2 = cad_angle_between(mol2.tmpICoords[i2], mol2.tmpComCoord,
        mol2.interfaceList[i2].coord, mol2.comCoord);
    const double a3 = cad_angle_between(mol1.tmpICoords[i1], com1.tmpComCoord,
        mol1.interfaceList[i1].coord, com1.comCoord);
    const double a4 = cad_angle_between(mol2.tmpICoords[i2], com2.tmpComCoord,
        mol2.interfaceList[i2].coord, com2.comCoord);

    std::cerr << "  Independently computed angles (rad): mol1->molCOM=" << a1
              << " (expect " << M_PI / 2.0 << "), mol2->molCOM=" << a2
              << " (expect " << M_PI << ")\n";

    EXPECT_NEAR(a1, M_PI / 2.0, 1e-10) << "mol1 was rotated 90 degrees";
    EXPECT_NEAR(a2, M_PI, 1e-10) << "mol2 was rotated 180 degrees (dot_theta clamps cos to -1)";
    EXPECT_NEAR(a3, M_PI / 2.0, 1e-10) << "complex COM == molecule COM here, so same 90 deg";
    EXPECT_NEAR(a4, M_PI, 1e-10) << "complex COM == molecule COM here, so same 180 deg";

    const CadSnapshot before = cad_snapshot(mol1, mol2, com1, com2, i1, i2);

    std::cerr << "  Calling calc_angular_displacement()...\n";
    calc_angular_displacement(i1, i2, mol1, mol2, com1, com2, moleculeList);
    std::cerr << "  Returned normally.\n";

    cad_expect_unchanged(before, mol1, mol2, com1, com2, i1, i2);
}

// -----------------------------------------------------------------------------
// Test 3: mol1 is a POINT particle (interface coincides with its COM).
//         This drives the first `if (v1.magnitude < 1E-12) return;` branch.
// -----------------------------------------------------------------------------
void cad_test_point_particle_mol1()
{
    std::cerr << "\n[TEST] cad_test_point_particle_mol1\n"
              << "  Source file: src/reactions/calc_angular_displacement.cpp\n"
              << "  Function:    calc_angular_displacement()\n"
              << "  Scenario:    molecule 1's interface sits exactly on its COM, so\n"
              << "               |v1| == 0 and the routine takes its FIRST early return.\n"
              << "  Pass:        no crash, arguments unchanged, and the measured\n"
              << "               v1 magnitude is below the 1E-12 cut-off.\n";

    // POINT particle: interface coordinate == COM coordinate.
    Molecule mol1 = cad_make_molecule(Coord{ 3.0, -2.0, 1.0 }, { Coord{ 3.0, -2.0, 1.0 } });
    Molecule mol2 = cad_make_molecule(Coord{ 0.0, 0.0, 0.0 }, { Coord{ 0.0, 1.0, 0.0 } });
    mol2.index = 1;
    mol2.myComIndex = 1;

    Complex com1 = cad_make_complex(mol1.comCoord, 0);
    Complex com2 = cad_make_complex(mol2.comCoord, 1);
    com2.index = 1;

    std::vector<Molecule> moleculeList{ mol1, mol2 };

    const int i1 = 0;
    const int i2 = 0;

    // Confirm the branch condition the function tests really is satisfied.
    Vector v1{ mol1.tmpICoords[i1] - mol1.tmpComCoord };
    v1.calc_magnitude();
    std::cerr << "  |v1| (mol1 tmp iface - tmp COM) = " << v1.magnitude << '\n';
    EXPECT_LT(v1.magnitude, 1e-12) << "mol1 must be a POINT particle for this branch";

    const CadSnapshot before = cad_snapshot(mol1, mol2, com1, com2, i1, i2);

    std::cerr << "  Calling calc_angular_displacement() (expects early return)...\n";
    calc_angular_displacement(i1, i2, mol1, mol2, com1, com2, moleculeList);
    std::cerr << "  Returned normally from the POINT-particle branch.\n";

    cad_expect_unchanged(before, mol1, mol2, com1, com2, i1, i2);
}

// -----------------------------------------------------------------------------
// Test 4: mol1 is normal but mol2 is a POINT particle.  This drives the SECOND
//         early return (`if (v3.magnitude < 1E-12) return;`), i.e. the routine
//         must first successfully evaluate molecule 1 and then bail out.
// -----------------------------------------------------------------------------
void cad_test_point_particle_mol2()
{
    std::cerr << "\n[TEST] cad_test_point_particle_mol2\n"
              << "  Source file: src/reactions/calc_angular_displacement.cpp\n"
              << "  Function:    calc_angular_displacement()\n"
              << "  Scenario:    molecule 1 is a normal rigid body (angle is computed),\n"
              << "               molecule 2 is a POINT particle -> SECOND early return.\n"
              << "  Pass:        no crash, arguments unchanged, |v1| > 0 but |v3| == 0.\n";

    Molecule mol1 = cad_make_molecule(Coord{ 0.0, 0.0, 0.0 }, { Coord{ 0.0, 0.0, 1.5 } });
    cad_rotate_tmp_about_z(mol1, mol1.tmpComCoord, M_PI / 4.0); // a rotation about z of a +z vector is a no-op angle

    // POINT particle for molecule 2.
    Molecule mol2 = cad_make_molecule(Coord{ -4.0, 4.0, 0.0 }, { Coord{ -4.0, 4.0, 0.0 } });
    mol2.index = 1;
    mol2.myComIndex = 1;

    Complex com1 = cad_make_complex(mol1.comCoord, 0);
    Complex com2 = cad_make_complex(mol2.comCoord, 1);
    com2.index = 1;

    std::vector<Molecule> moleculeList{ mol1, mol2 };

    const int i1 = 0;
    const int i2 = 0;

    Vector v1{ mol1.tmpICoords[i1] - mol1.tmpComCoord };
    v1.calc_magnitude();
    Vector v3{ mol2.tmpICoords[i2] - mol2.tmpComCoord };
    v3.calc_magnitude();
    std::cerr << "  |v1| = " << v1.magnitude << " (should be non-zero), |v3| = "
              << v3.magnitude << " (should be zero)\n";

    EXPECT_GT(v1.magnitude, 1e-12) << "molecule 1 must NOT be a point particle here";
    EXPECT_LT(v3.magnitude, 1e-12) << "molecule 2 must be a POINT particle here";

    const CadSnapshot before = cad_snapshot(mol1, mol2, com1, com2, i1, i2);

    std::cerr << "  Calling calc_angular_displacement() (expects 2nd early return)...\n";
    calc_angular_displacement(i1, i2, mol1, mol2, com1, com2, moleculeList);
    std::cerr << "  Returned normally.\n";

    cad_expect_unchanged(before, mol1, mol2, com1, com2, i1, i2);
}

// -----------------------------------------------------------------------------
// Test 5: both molecules are normal, but molecule 1's interface lies exactly on
//         the *complex* tmp COM, so the third vector (v5) has zero magnitude and
//         the THIRD early return is taken.
// -----------------------------------------------------------------------------
void cad_test_zero_complex_vector()
{
    std::cerr << "\n[TEST] cad_test_zero_complex_vector\n"
              << "  Source file: src/reactions/calc_angular_displacement.cpp\n"
              << "  Function:    calc_angular_displacement()\n"
              << "  Scenario:    mol1's tmp interface coincides with reactCom1.tmpComCoord,\n"
              << "               so |v5| == 0 and the THIRD early return fires.\n"
              << "  Pass:        no crash, arguments unchanged, |v5| below 1E-12.\n";

    Molecule mol1 = cad_make_molecule(Coord{ 0.0, 0.0, 0.0 }, { Coord{ 1.0, 0.0, 0.0 } });
    Molecule mol2 = cad_make_molecule(Coord{ 6.0, 0.0, 0.0 }, { Coord{ 5.0, 0.0, 0.0 } });
    mol2.index = 1;
    mol2.myComIndex = 1;

    // Complex 1's temporary COM is deliberately put on top of mol1's tmp interface.
    Complex com1 = cad_make_complex(Coord{ 0.5, 0.0, 0.0 }, 0);
    com1.tmpComCoord = mol1.tmpICoords[0];
    Complex com2 = cad_make_complex(mol2.comCoord, 1);
    com2.index = 1;

    std::vector<Molecule> moleculeList{ mol1, mol2 };

    const int i1 = 0;
    const int i2 = 0;

    Vector v5{ mol1.tmpICoords[i1] - com1.tmpComCoord };
    v5.calc_magnitude();
    std::cerr << "  |v5| (mol1 tmp iface - complex1 tmp COM) = " << v5.magnitude << '\n';
    EXPECT_LT(v5.magnitude, 1e-12) << "v5 must be degenerate to hit the third early return";

    const CadSnapshot before = cad_snapshot(mol1, mol2, com1, com2, i1, i2);

    std::cerr << "  Calling calc_angular_displacement() (expects 3rd early return)...\n";
    calc_angular_displacement(i1, i2, mol1, mol2, com1, com2, moleculeList);
    std::cerr << "  Returned normally.\n";

    cad_expect_unchanged(before, mol1, mol2, com1, com2, i1, i2);
}

// -----------------------------------------------------------------------------
// Test 6: realistic multi-molecule complexes whose COM differs from the
//         molecule COM, and a non-zero interface index, so the routine runs all
//         the way to the end (no early return) and indexes tmpICoords/
//         interfaceList at index 1 rather than 0.
// -----------------------------------------------------------------------------
void cad_test_offset_complex_and_nonzero_iface_index()
{
    std::cerr << "\n[TEST] cad_test_offset_complex_and_nonzero_iface_index\n"
              << "  Source file: src/reactions/calc_angular_displacement.cpp\n"
              << "  Function:    calc_angular_displacement()\n"
              << "  Scenario:    two-interface molecules, ifaceIndex1 = ifaceIndex2 = 1,\n"
              << "               complex COM offset from the molecule COM, tmp coords\n"
              << "               rotated -> the routine runs to completion.\n"
              << "  Pass:        the four angles are finite and non-negative, and the\n"
              << "               function leaves every argument untouched.\n";

    // Two interfaces each; we deliberately use index 1 in the call.
    Molecule mol1 = cad_make_molecule(Coord{ 0.0, 0.0, 0.0 },
        { Coord{ 1.0, 0.0, 0.0 }, Coord{ 0.0, 1.0, 0.0 } });
    Molecule mol2 = cad_make_molecule(Coord{ 8.0, 0.0, 0.0 },
        { Coord{ 9.0, 0.0, 0.0 }, Coord{ 8.0, -1.0, 0.0 } });
    mol2.index = 1;
    mol2.myComIndex = 1;

    // Rotate the temporary coordinates of both molecules about the origin,
    // which moves both the molecule tmp COM and its tmp interfaces.
    cad_rotate_tmp_about_z(mol1, Coord{ 0.0, 0.0, 0.0 }, M_PI / 6.0);
    cad_rotate_tmp_about_z(mol2, Coord{ 0.0, 0.0, 0.0 }, -M_PI / 6.0);

    // Complex COMs offset from the molecule COMs (as in a multi-molecule complex).
    Complex com1 = cad_make_complex(Coord{ -1.0, -1.0, 0.0 }, 0);
    com1.tmpComCoord = Coord{ -0.8, -1.1, 0.0 };
    Complex com2 = cad_make_complex(Coord{ 9.0, 1.0, 0.0 }, 1);
    com2.index = 1;
    com2.tmpComCoord = Coord{ 9.2, 0.9, 0.0 };

    std::vector<Molecule> moleculeList{ mol1, mol2 };

    const int i1 = 1; // exercise a non-zero interface index
    const int i2 = 1;

    // Independently reproduce all four angles and sanity-check them.
    const double a1 = cad_angle_between(mol1.tmpICoords[i1], mol1.tmpComCoord,
        mol1.interfaceList[i1].coord, mol1.comCoord);
    const double a2 = cad_angle_between(mol2.tmpICoords[i2], mol2.tmpComCoord,
        mol2.interfaceList[i2].coord, mol2.comCoord);
    const double a3 = cad_angle_between(mol1.tmpICoords[i1], com1.tmpComCoord,
        mol1.interfaceList[i1].coord, com1.comCoord);
    const double a4 = cad_angle_between(mol2.tmpICoords[i2], com2.tmpComCoord,
        mol2.interfaceList[i2].coord, com2.comCoord);

    std::cerr << "  Angles (rad): mol1->molCOM=" << a1 << ", mol2->molCOM=" << a2
              << ", mol1->complexCOM=" << a3 << ", mol2->complexCOM=" << a4 << '\n';

    // Every angle returned by Vector::dot_theta() is an acos() result, so it must
    // lie in [0, pi] and must be finite.
    EXPECT_TRUE(std::isfinite(a1)) << "mol1 molecule-level angle must be finite";
    EXPECT_TRUE(std::isfinite(a2)) << "mol2 molecule-level angle must be finite";
    EXPECT_TRUE(std::isfinite(a3)) << "mol1 complex-level angle must be finite";
    EXPECT_TRUE(std::isfinite(a4)) << "mol2 complex-level angle must be finite";
    EXPECT_GE(a1, 0.0);
    EXPECT_LE(a1, M_PI);
    EXPECT_GE(a2, 0.0);
    EXPECT_LE(a2, M_PI);
    EXPECT_GE(a3, 0.0);
    EXPECT_LE(a3, M_PI);
    EXPECT_GE(a4, 0.0);
    EXPECT_LE(a4, M_PI);

    // The molecule-level rotations were exactly +/- 30 degrees about z applied to
    // vectors lying in the xy-plane, so both should come out as pi/6.
    EXPECT_NEAR(a1, M_PI / 6.0, 1e-10) << "mol1 tmp frame was rotated by 30 degrees";
    EXPECT_NEAR(a2, M_PI / 6.0, 1e-10) << "mol2 tmp frame was rotated by 30 degrees";

    const CadSnapshot before = cad_snapshot(mol1, mol2, com1, com2, i1, i2);

    std::cerr << "  Calling calc_angular_displacement() (full path, no early return)...\n";
    calc_angular_displacement(i1, i2, mol1, mol2, com1, com2, moleculeList);
    std::cerr << "  Returned normally.\n";

    cad_expect_unchanged(before, mol1, mol2, com1, com2, i1, i2);

    // The moleculeList argument is likewise never used for writing; verify the
    // copies stored there are untouched too.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 0.0) << "moleculeList must not be modified";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 8.0) << "moleculeList must not be modified";
    EXPECT_EQ(moleculeList.size(), 2u) << "moleculeList size must not change";
}

// -----------------------------------------------------------------------------
// Test 7: aliasing case - the very common situation in association where both
//         reacting molecules belong to the *same* complex, so reactCom1 and
//         reactCom2 refer to identical data.  Just make sure that is handled.
// -----------------------------------------------------------------------------
void cad_test_same_complex_for_both_molecules()
{
    std::cerr << "\n[TEST] cad_test_same_complex_for_both_molecules\n"
              << "  Source file: src/reactions/calc_angular_displacement.cpp\n"
              << "  Function:    calc_angular_displacement()\n"
              << "  Scenario:    both molecules are members of one complex, so the\n"
              << "               same Complex object is passed for reactCom1/reactCom2.\n"
              << "  Pass:        no crash and the shared Complex is left unmodified.\n";

    Molecule mol1 = cad_make_molecule(Coord{ -2.0, 0.0, 0.0 }, { Coord{ -1.0, 0.0, 0.0 } });
    Molecule mol2 = cad_make_molecule(Coord{ 2.0, 0.0, 0.0 }, { Coord{ 1.0, 0.0, 0.0 } });
    mol2.index = 1;
    mol2.myComIndex = 0; // same complex as mol1

    // Small rotations of the temporary frames so the angles are non-trivial.
    cad_rotate_tmp_about_z(mol1, mol1.tmpComCoord, 0.2);
    cad_rotate_tmp_about_z(mol2, mol2.tmpComCoord, -0.35);

    // One shared complex containing both molecules.
    Complex sharedCom = cad_make_complex(Coord{ 0.0, 0.0, 0.0 }, 0);
    sharedCom.memberList.push_back(1);
    sharedCom.tmpComCoord = Coord{ 0.05, 0.02, 0.0 };

    std::vector<Molecule> moleculeList{ mol1, mol2 };

    const int i1 = 0;
    const int i2 = 0;

    const CadSnapshot before = cad_snapshot(mol1, mol2, sharedCom, sharedCom, i1, i2);

    std::cerr << "  Calling calc_angular_displacement() with the same Complex twice...\n";
    calc_angular_displacement(i1, i2, mol1, mol2, sharedCom, sharedCom, moleculeList);
    std::cerr << "  Returned normally.\n";

    cad_expect_unchanged(before, mol1, mol2, sharedCom, sharedCom, i1, i2);

    // The shared complex's member list must be intact.
    ASSERT_EQ(sharedCom.memberList.size(), 2u) << "shared complex memberList must not change";
    EXPECT_EQ(sharedCom.memberList[0], 0) << "member 0 preserved";
    EXPECT_EQ(sharedCom.memberList[1], 1) << "member 1 preserved";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so the
// framework reports individual results and every test still runs even if an
// earlier one records a (non-fatal) failure.
// -----------------------------------------------------------------------------
TEST(CalcAngularDisplacement, ZeroDisplacement) { cad_test_zero_displacement(); }
TEST(CalcAngularDisplacement, KnownRotation) { cad_test_known_rotation(); }
TEST(CalcAngularDisplacement, PointParticleMol1) { cad_test_point_particle_mol1(); }
TEST(CalcAngularDisplacement, PointParticleMol2) { cad_test_point_particle_mol2(); }
TEST(CalcAngularDisplacement, ZeroComplexVector) { cad_test_zero_complex_vector(); }
TEST(CalcAngularDisplacement, OffsetComplexAndNonZeroIfaceIndex)
{
    cad_test_offset_complex_and_nonzero_iface_index();
}
TEST(CalcAngularDisplacement, SameComplexForBothMolecules)
{
    cad_test_same_complex_for_both_molecules();
}