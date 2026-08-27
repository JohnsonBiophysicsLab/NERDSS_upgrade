/*! \file test_get_distance_to_surface.cpp
 *
 * ### Unit test for ../src/reactions/get_distance_to_surface.cpp
 *
 * Function under test:
 *
 * \code
 * bool get_distance_to_surface(int pro1, int pro2, int iface1, int iface2,
 *                              int rxnIndex, int rateIndex,
 *                              bool isStateChangeBackRxn,
 *                              double& sep, double& R1, double Rmax,
 *                              std::vector<Complex>& complexList,
 *                              const ForwardRxn& currRxn,
 *                              std::vector<Molecule>& moleculeList,
 *                              const Membrane& membraneObject);
 * \endcode
 *
 * The routine measures the distance R1 between the reacting interface of
 * molecule `pro1` and the implicit-lipid surface.  Three distinct branches
 * exist in the implementation:
 *
 *   1. The parent complex is already ON the surface (`Complex::OnSurface`
 *      is true)                       -> R1 is forced to exactly 0.
 *   2. Spherical membrane geometry     -> R1 = |sphereR - |r||  where |r| is the
 *                                         full 3D distance of the interface
 *                                         from the origin.
 *   3. Rectangular (box) geometry      -> R1 = |z - (-boxZ/2) - lipidLength|,
 *                                         i.e. only the z-component matters.
 *
 * If (and only if) R1 < Rmax the function registers a possible reaction by
 * appending to `crossbase`, `mycrossint` and `crossrxn` of molecule `pro1`,
 * incrementing `ncross` on its parent complex, and returning true.
 * Otherwise nothing is touched and false is returned.
 *
 * Every assertion below is written directly against those documented
 * behaviours; verbose console output describes what is being checked.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers.  All names carry the "gdts_" prefix (get_distance_to_surface)
// so they cannot collide with helpers from other files in the test suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal, fully-initialised Molecule with one interface.
 *
 * The function under test only reads:
 *   - moleculeList[pro1].interfaceList[iface1].coord
 *   - moleculeList[pro1].myComIndex
 * and writes to crossbase / mycrossint / crossrxn.  Everything else is still
 * initialised here so the object is in a sane state.
 *
 * \param[in] index         index of this molecule inside moleculeList
 * \param[in] comIndex      index of its parent complex inside complexList
 * \param[in] ifaceCoord    absolute coordinate given to the single interface
 */
Molecule gdts_make_molecule(int index, int comIndex, const Coord& ifaceCoord)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.comCoord = ifaceCoord;

    // Single interface placed exactly on the molecule centre of mass.
    Molecule::Iface iface;
    iface.coord = ifaceCoord;
    iface.relIndex = 0;
    iface.index = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // Encounter bookkeeping vectors start empty.
    mol.crossbase.clear();
    mol.mycrossint.clear();
    mol.crossrxn.clear();

    return mol;
}

/*! \brief Build a minimal Complex owning exactly one member molecule. */
Complex gdts_make_complex(int index, int memberMolIndex, bool onSurface)
{
    Complex com;
    com.index = index;
    com.isEmpty = false;
    com.OnSurface = onSurface;
    com.ncross = 0;
    com.mass = 1.0;
    com.radius = 1.0;
    com.D = Coord { 1.0, 1.0, 1.0 };
    com.Dr = Coord { 0.01, 0.01, 0.01 };
    com.memberList.clear();
    com.memberList.push_back(memberMolIndex);
    return com;
}

/*! \brief Build a box (non-spherical) Membrane with a given z-height. */
Membrane gdts_make_box_membrane(double boxZ, double lipidLength)
{
    Membrane mem;
    mem.isSphere = false;
    mem.isBox = true;
    mem.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, boxZ });
    mem.lipidLength = lipidLength;
    mem.implicitLipid = true;
    return mem;
}

/*! \brief Build a spherical Membrane of the requested radius. */
Membrane gdts_make_sphere_membrane(double sphereR)
{
    Membrane mem;
    mem.isSphere = true;
    mem.isBox = false;
    mem.sphereR = sphereR;
    mem.lipidLength = 0.0;
    mem.implicitLipid = true;
    // A water box is still filled in (unused by the sphere branch).
    mem.waterBox = Membrane::WaterBox(std::vector<double> { 2.0 * sphereR, 2.0 * sphereR, 2.0 * sphereR });
    return mem;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: complex already on the surface -> R1 is forced to exactly zero.
// -----------------------------------------------------------------------------
void test_gdts_on_surface_returns_zero_distance()
{
    std::cerr << "\n[TEST] test_gdts_on_surface_returns_zero_distance\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (OnSurface branch)\n"
              << "  Scenario:      parent Complex::OnSurface == true, interface far\n"
              << "                 away from the origin in x/y/z.\n"
              << "  Pass criteria: R1 == 0 and sep == 0 regardless of coordinates,\n"
              << "                 the reaction is registered and true is returned.\n";

    // Interface deliberately placed far from the membrane; the OnSurface flag
    // must short-circuit any geometric calculation.
    std::vector<Molecule> moleculeList { gdts_make_molecule(0, 0, Coord { 10.0, 20.0, 30.0 }) };
    std::vector<Complex> complexList { gdts_make_complex(0, 0, /*onSurface=*/true) };

    Membrane membraneObject = gdts_make_box_membrane(/*boxZ=*/100.0, /*lipidLength=*/2.0);
    ForwardRxn currRxn; // unused by the function, default constructed

    double sep = -1.0;
    double R1 = -1.0;
    const double Rmax = 1.0;

    std::cerr << "  Calling get_distance_to_surface with Rmax = " << Rmax << "\n";
    bool result = get_distance_to_surface(/*pro1=*/0, /*pro2=*/7, /*iface1=*/0, /*iface2=*/3,
        /*rxnIndex=*/5, /*rateIndex=*/2, /*isStateChangeBackRxn=*/false,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    std::cerr << "  Returned " << std::boolalpha << result << ", R1 = " << R1
              << ", sep = " << sep << "\n";

    // The OnSurface branch hard-codes R1 to 0, and sep is a copy of R1.
    EXPECT_DOUBLE_EQ(R1, 0.0) << "R1 must be exactly 0 when the complex is on the surface";
    EXPECT_DOUBLE_EQ(sep, 0.0) << "sep must mirror R1 (0) when the complex is on the surface";

    // 0 < Rmax, so the encounter must have been registered.
    EXPECT_TRUE(result) << "R1 (0) < Rmax (1) so the function must return true";
    ASSERT_EQ(moleculeList[0].crossbase.size(), 1u) << "crossbase should hold exactly one partner";
    EXPECT_EQ(moleculeList[0].crossbase[0], 7) << "crossbase must record pro2 == 7";
    ASSERT_EQ(moleculeList[0].mycrossint.size(), 1u) << "mycrossint should hold exactly one entry";
    EXPECT_EQ(moleculeList[0].mycrossint[0], 0) << "mycrossint must record iface1 == 0";
    ASSERT_EQ(moleculeList[0].crossrxn.size(), 1u) << "crossrxn should hold exactly one entry";
    EXPECT_EQ(moleculeList[0].crossrxn[0][0], 5) << "crossrxn[0] must record rxnIndex == 5";
    EXPECT_EQ(moleculeList[0].crossrxn[0][1], 2) << "crossrxn[1] must record rateIndex == 2";
    EXPECT_EQ(moleculeList[0].crossrxn[0][2], 0) << "crossrxn[2] must record false as 0";
    EXPECT_EQ(complexList[0].ncross, 1) << "ncross must be incremented once";
}

// -----------------------------------------------------------------------------
// Test 2: box geometry, interface within Rmax of the lipid plane.
// -----------------------------------------------------------------------------
void test_gdts_box_within_rmax()
{
    std::cerr << "\n[TEST] test_gdts_box_within_rmax\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (box branch)\n"
              << "  Scenario:      box of height 100 (bottom at z = -50), lipidLength 2,\n"
              << "                 interface at z = -45.\n"
              << "  Pass criteria: R1 == |(-45) - (-50) - 2| == 3, reaction registered.\n";

    std::vector<Molecule> moleculeList { gdts_make_molecule(0, 0, Coord { 12.0, -7.0, -45.0 }) };
    std::vector<Complex> complexList { gdts_make_complex(0, 0, /*onSurface=*/false) };

    Membrane membraneObject = gdts_make_box_membrane(/*boxZ=*/100.0, /*lipidLength=*/2.0);
    ForwardRxn currRxn;

    double sep = -1.0;
    double R1 = -1.0;
    const double Rmax = 5.0;

    std::cerr << "  Calling get_distance_to_surface with Rmax = " << Rmax << "\n";
    bool result = get_distance_to_surface(0, /*pro2=*/11, 0, 0,
        /*rxnIndex=*/1, /*rateIndex=*/0, /*isStateChangeBackRxn=*/true,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    std::cerr << "  Returned " << std::boolalpha << result << ", R1 = " << R1
              << " (expected 3), sep = " << sep << "\n";

    // dz = -45 - (-50) - 2 = 3; R1 = sqrt(dz*dz) = 3.
    EXPECT_DOUBLE_EQ(R1, 3.0) << "R1 must equal the z-offset above the lipid plane";
    EXPECT_DOUBLE_EQ(sep, R1) << "sep is assigned directly from R1";

    EXPECT_TRUE(result) << "R1 (3) < Rmax (5) so the function must return true";
    ASSERT_EQ(moleculeList[0].crossrxn.size(), 1u) << "one crossrxn entry expected";
    EXPECT_EQ(moleculeList[0].crossrxn[0][0], 1) << "crossrxn[0] must record rxnIndex == 1";
    EXPECT_EQ(moleculeList[0].crossrxn[0][1], 0) << "crossrxn[1] must record rateIndex == 0";
    EXPECT_EQ(moleculeList[0].crossrxn[0][2], 1) << "crossrxn[2] must record true as 1";
    EXPECT_EQ(complexList[0].ncross, 1) << "ncross must be incremented once";
}

// -----------------------------------------------------------------------------
// Test 3: box geometry, interface too far away -> nothing is recorded.
// -----------------------------------------------------------------------------
void test_gdts_box_beyond_rmax_no_side_effects()
{
    std::cerr << "\n[TEST] test_gdts_box_beyond_rmax_no_side_effects\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (box branch, R1 >= Rmax)\n"
              << "  Scenario:      same geometry as before (R1 == 3) but Rmax == 2.\n"
              << "  Pass criteria: false is returned, R1/sep are still filled in, and\n"
              << "                 crossbase/mycrossint/crossrxn/ncross are untouched.\n";

    std::vector<Molecule> moleculeList { gdts_make_molecule(0, 0, Coord { 0.0, 0.0, -45.0 }) };
    std::vector<Complex> complexList { gdts_make_complex(0, 0, /*onSurface=*/false) };

    Membrane membraneObject = gdts_make_box_membrane(/*boxZ=*/100.0, /*lipidLength=*/2.0);
    ForwardRxn currRxn;

    double sep = -1.0;
    double R1 = -1.0;
    const double Rmax = 2.0;

    std::cerr << "  Calling get_distance_to_surface with Rmax = " << Rmax << "\n";
    bool result = get_distance_to_surface(0, 4, 0, 0, 9, 9, false,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    std::cerr << "  Returned " << std::boolalpha << result << ", R1 = " << R1 << "\n";

    // The distance is still computed even when the reaction is not registered.
    EXPECT_DOUBLE_EQ(R1, 3.0) << "R1 is computed before the Rmax comparison";
    EXPECT_DOUBLE_EQ(sep, 3.0) << "sep is computed before the Rmax comparison";

    EXPECT_FALSE(result) << "R1 (3) >= Rmax (2) so the function must return false";
    EXPECT_TRUE(moleculeList[0].crossbase.empty()) << "crossbase must remain empty";
    EXPECT_TRUE(moleculeList[0].mycrossint.empty()) << "mycrossint must remain empty";
    EXPECT_TRUE(moleculeList[0].crossrxn.empty()) << "crossrxn must remain empty";
    EXPECT_EQ(complexList[0].ncross, 0) << "ncross must not be incremented";
}

// -----------------------------------------------------------------------------
// Test 4: box geometry ignores the x and y coordinates entirely.
// -----------------------------------------------------------------------------
void test_gdts_box_ignores_xy()
{
    std::cerr << "\n[TEST] test_gdts_box_ignores_xy\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (box branch)\n"
              << "  Scenario:      two interfaces with identical z but wildly different\n"
              << "                 x and y coordinates.\n"
              << "  Pass criteria: both produce exactly the same R1 (box branch uses\n"
              << "                 only the z component).\n";

    Membrane membraneObject = gdts_make_box_membrane(/*boxZ=*/40.0, /*lipidLength=*/1.0);
    ForwardRxn currRxn;
    const double Rmax = 100.0; // large, so both calls register (irrelevant here)

    // First molecule: interface at (0, 0, -14).
    std::vector<Molecule> molA { gdts_make_molecule(0, 0, Coord { 0.0, 0.0, -14.0 }) };
    std::vector<Complex> comA { gdts_make_complex(0, 0, false) };
    double sepA = -1.0, R1A = -1.0;
    get_distance_to_surface(0, 1, 0, 0, 0, 0, false, sepA, R1A, Rmax, comA, currRxn, molA, membraneObject);

    // Second molecule: same z, but very different x/y.
    std::vector<Molecule> molB { gdts_make_molecule(0, 0, Coord { 999.0, -321.0, -14.0 }) };
    std::vector<Complex> comB { gdts_make_complex(0, 0, false) };
    double sepB = -1.0, R1B = -1.0;
    get_distance_to_surface(0, 1, 0, 0, 0, 0, false, sepB, R1B, Rmax, comB, currRxn, molB, membraneObject);

    std::cerr << "  R1 with (x,y) = (0,0)        : " << R1A << "\n";
    std::cerr << "  R1 with (x,y) = (999,-321)   : " << R1B << "\n";

    // dz = -14 - (-20) - 1 = 5 for both.
    EXPECT_DOUBLE_EQ(R1A, 5.0) << "R1 should be |z + boxZ/2 - lipidLength| == 5";
    EXPECT_DOUBLE_EQ(R1B, R1A) << "The box branch must be independent of x and y";
}

// -----------------------------------------------------------------------------
// Test 5: spherical geometry, interface inside the sphere.
// -----------------------------------------------------------------------------
void test_gdts_sphere_inside()
{
    std::cerr << "\n[TEST] test_gdts_sphere_inside\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (sphere branch)\n"
              << "  Scenario:      sphereR = 100, interface at (0, 0, 90) so |r| = 90.\n"
              << "  Pass criteria: R1 == |100 - 90| == 10 and the reaction is registered\n"
              << "                 because 10 < Rmax (15).\n";

    std::vector<Molecule> moleculeList { gdts_make_molecule(0, 0, Coord { 0.0, 0.0, 90.0 }) };
    std::vector<Complex> complexList { gdts_make_complex(0, 0, /*onSurface=*/false) };

    Membrane membraneObject = gdts_make_sphere_membrane(/*sphereR=*/100.0);
    ForwardRxn currRxn;

    double sep = -1.0;
    double R1 = -1.0;
    const double Rmax = 15.0;

    std::cerr << "  Calling get_distance_to_surface with Rmax = " << Rmax << "\n";
    bool result = get_distance_to_surface(0, /*pro2=*/2, 0, 0,
        /*rxnIndex=*/3, /*rateIndex=*/1, /*isStateChangeBackRxn=*/false,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    std::cerr << "  Returned " << std::boolalpha << result << ", R1 = " << R1
              << " (expected 10)\n";

    EXPECT_NEAR(R1, 10.0, 1e-12) << "R1 must be |sphereR - |r||";
    EXPECT_DOUBLE_EQ(sep, R1) << "sep is assigned directly from R1";
    EXPECT_TRUE(result) << "R1 (10) < Rmax (15) so the function must return true";
    EXPECT_EQ(complexList[0].ncross, 1) << "ncross must be incremented once";
    ASSERT_EQ(moleculeList[0].crossbase.size(), 1u) << "one crossbase entry expected";
    EXPECT_EQ(moleculeList[0].crossbase[0], 2) << "crossbase must record pro2 == 2";
}

// -----------------------------------------------------------------------------
// Test 6: spherical geometry, interface OUTSIDE the sphere -> absolute value.
// -----------------------------------------------------------------------------
void test_gdts_sphere_outside_uses_absolute_value()
{
    std::cerr << "\n[TEST] test_gdts_sphere_outside_uses_absolute_value\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (sphere branch)\n"
              << "  Scenario:      sphereR = 100, interface at (110, 0, 0) so |r| = 110,\n"
              << "                 which lies outside the sphere.\n"
              << "  Pass criteria: R1 == |100 - 110| == 10 (a positive distance), i.e.\n"
              << "                 std::abs is applied to the signed difference.\n";

    std::vector<Molecule> moleculeList { gdts_make_molecule(0, 0, Coord { 110.0, 0.0, 0.0 }) };
    std::vector<Complex> complexList { gdts_make_complex(0, 0, /*onSurface=*/false) };

    Membrane membraneObject = gdts_make_sphere_membrane(/*sphereR=*/100.0);
    ForwardRxn currRxn;

    double sep = -1.0;
    double R1 = -1.0;
    const double Rmax = 20.0;

    std::cerr << "  Calling get_distance_to_surface with Rmax = " << Rmax << "\n";
    bool result = get_distance_to_surface(0, 1, 0, 0, 0, 0, false,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    std::cerr << "  Returned " << std::boolalpha << result << ", R1 = " << R1
              << " (expected 10, and must be non-negative)\n";

    EXPECT_NEAR(R1, 10.0, 1e-12) << "R1 must be the absolute distance to the shell";
    EXPECT_GE(R1, 0.0) << "R1 must never be negative";
    EXPECT_TRUE(result) << "R1 (10) < Rmax (20) so the function must return true";
}

// -----------------------------------------------------------------------------
// Test 7: spherical geometry, interface too far from the shell.
// -----------------------------------------------------------------------------
void test_gdts_sphere_beyond_rmax()
{
    std::cerr << "\n[TEST] test_gdts_sphere_beyond_rmax\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (sphere branch, R1 >= Rmax)\n"
              << "  Scenario:      sphereR = 100, interface at the origin so |r| = 0 and\n"
              << "                 R1 == 100, with Rmax == 5.\n"
              << "  Pass criteria: false returned and no bookkeeping performed.\n";

    std::vector<Molecule> moleculeList { gdts_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }) };
    std::vector<Complex> complexList { gdts_make_complex(0, 0, /*onSurface=*/false) };

    Membrane membraneObject = gdts_make_sphere_membrane(/*sphereR=*/100.0);
    ForwardRxn currRxn;

    double sep = -1.0;
    double R1 = -1.0;
    const double Rmax = 5.0;

    std::cerr << "  Calling get_distance_to_surface with Rmax = " << Rmax << "\n";
    bool result = get_distance_to_surface(0, 1, 0, 0, 0, 0, false,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    std::cerr << "  Returned " << std::boolalpha << result << ", R1 = " << R1
              << " (expected 100)\n";

    EXPECT_NEAR(R1, 100.0, 1e-12) << "A particle at the centre is sphereR away from the shell";
    EXPECT_FALSE(result) << "R1 (100) >= Rmax (5) so the function must return false";
    EXPECT_TRUE(moleculeList[0].crossbase.empty()) << "crossbase must remain empty";
    EXPECT_TRUE(moleculeList[0].mycrossint.empty()) << "mycrossint must remain empty";
    EXPECT_TRUE(moleculeList[0].crossrxn.empty()) << "crossrxn must remain empty";
    EXPECT_EQ(complexList[0].ncross, 0) << "ncross must not be incremented";
}

// -----------------------------------------------------------------------------
// Test 8: repeated successful calls accumulate encounters on the same molecule.
// -----------------------------------------------------------------------------
void test_gdts_repeated_calls_accumulate()
{
    std::cerr << "\n[TEST] test_gdts_repeated_calls_accumulate\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (bookkeeping)\n"
              << "  Scenario:      the same molecule is evaluated twice against two\n"
              << "                 different partners, both times within Rmax.\n"
              << "  Pass criteria: crossbase/mycrossint/crossrxn each grow to size 2 and\n"
              << "                 ncross ends at 2, preserving call order.\n";

    // Complex flagged as OnSurface so R1 is always 0 and both calls succeed.
    std::vector<Molecule> moleculeList { gdts_make_molecule(0, 0, Coord { 1.0, 2.0, 3.0 }) };
    std::vector<Complex> complexList { gdts_make_complex(0, 0, /*onSurface=*/true) };

    Membrane membraneObject = gdts_make_box_membrane(/*boxZ=*/50.0, /*lipidLength=*/0.0);
    ForwardRxn currRxn;

    double sep = 0.0;
    double R1 = 0.0;
    const double Rmax = 1.0;

    std::cerr << "  First call  (pro2 = 21, iface1 = 0, rxn 4, rate 1, backRxn true)\n";
    bool first = get_distance_to_surface(0, /*pro2=*/21, /*iface1=*/0, 0,
        /*rxnIndex=*/4, /*rateIndex=*/1, /*isStateChangeBackRxn=*/true,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    std::cerr << "  Second call (pro2 = 33, iface1 = 0, rxn 6, rate 0, backRxn false)\n";
    bool second = get_distance_to_surface(0, /*pro2=*/33, /*iface1=*/0, 0,
        /*rxnIndex=*/6, /*rateIndex=*/0, /*isStateChangeBackRxn=*/false,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    EXPECT_TRUE(first) << "First call should register an encounter";
    EXPECT_TRUE(second) << "Second call should register an encounter";

    ASSERT_EQ(moleculeList[0].crossbase.size(), 2u) << "crossbase should hold two partners";
    EXPECT_EQ(moleculeList[0].crossbase[0], 21) << "first partner recorded first";
    EXPECT_EQ(moleculeList[0].crossbase[1], 33) << "second partner recorded second";

    ASSERT_EQ(moleculeList[0].mycrossint.size(), 2u) << "mycrossint should hold two entries";
    EXPECT_EQ(moleculeList[0].mycrossint[0], 0) << "iface1 of the first call";
    EXPECT_EQ(moleculeList[0].mycrossint[1], 0) << "iface1 of the second call";

    ASSERT_EQ(moleculeList[0].crossrxn.size(), 2u) << "crossrxn should hold two entries";
    EXPECT_EQ(moleculeList[0].crossrxn[0][0], 4) << "first call rxnIndex";
    EXPECT_EQ(moleculeList[0].crossrxn[0][1], 1) << "first call rateIndex";
    EXPECT_EQ(moleculeList[0].crossrxn[0][2], 1) << "first call isStateChangeBackRxn == true -> 1";
    EXPECT_EQ(moleculeList[0].crossrxn[1][0], 6) << "second call rxnIndex";
    EXPECT_EQ(moleculeList[0].crossrxn[1][1], 0) << "second call rateIndex";
    EXPECT_EQ(moleculeList[0].crossrxn[1][2], 0) << "second call isStateChangeBackRxn == false -> 0";

    EXPECT_EQ(complexList[0].ncross, 2) << "ncross must count both successful evaluations";

    std::cerr << "  Final ncross = " << complexList[0].ncross
              << ", crossbase size = " << moleculeList[0].crossbase.size() << "\n";
}

// -----------------------------------------------------------------------------
// Test 9: the correct molecule/complex pair is used when the lists hold several
//         entries (index plumbing through pro1 and myComIndex).
// -----------------------------------------------------------------------------
void test_gdts_uses_correct_indices()
{
    std::cerr << "\n[TEST] test_gdts_uses_correct_indices\n"
              << "  Source file:   src/reactions/get_distance_to_surface.cpp\n"
              << "  Function:      get_distance_to_surface (index plumbing)\n"
              << "  Scenario:      moleculeList and complexList hold three entries; the\n"
              << "                 call targets pro1 == 2 whose myComIndex == 2.\n"
              << "  Pass criteria: R1 is computed from molecule 2's interface and only\n"
              << "                 molecule 2 / complex 2 are modified.\n";

    // Three molecules at very different heights, each with its own complex.
    std::vector<Molecule> moleculeList {
        gdts_make_molecule(0, 0, Coord { 0.0, 0.0, 40.0 }),   // far from the plane
        gdts_make_molecule(1, 1, Coord { 0.0, 0.0, 20.0 }),   // also far
        gdts_make_molecule(2, 2, Coord { 0.0, 0.0, -47.0 })   // close to the plane
    };
    std::vector<Complex> complexList {
        gdts_make_complex(0, 0, false),
        gdts_make_complex(1, 1, false),
        gdts_make_complex(2, 2, false)
    };

    // Box of height 100 -> bottom at z = -50, lipidLength 1 -> plane at z = -49.
    Membrane membraneObject = gdts_make_box_membrane(/*boxZ=*/100.0, /*lipidLength=*/1.0);
    ForwardRxn currRxn;

    double sep = -1.0;
    double R1 = -1.0;
    const double Rmax = 5.0;

    std::cerr << "  Calling get_distance_to_surface for pro1 = 2\n";
    bool result = get_distance_to_surface(/*pro1=*/2, /*pro2=*/0, /*iface1=*/0, 0,
        /*rxnIndex=*/0, /*rateIndex=*/0, /*isStateChangeBackRxn=*/false,
        sep, R1, Rmax, complexList, currRxn, moleculeList, membraneObject);

    std::cerr << "  Returned " << std::boolalpha << result << ", R1 = " << R1
              << " (expected 2)\n";

    // dz = -47 - (-50) - 1 = 2 for molecule 2.
    EXPECT_DOUBLE_EQ(R1, 2.0) << "R1 must come from molecule 2's interface, not 0 or 1";
    EXPECT_TRUE(result) << "R1 (2) < Rmax (5) so the function must return true";

    // Only molecule 2 / complex 2 should have been touched.
    EXPECT_TRUE(moleculeList[0].crossbase.empty()) << "molecule 0 must be untouched";
    EXPECT_TRUE(moleculeList[1].crossbase.empty()) << "molecule 1 must be untouched";
    ASSERT_EQ(moleculeList[2].crossbase.size(), 1u) << "molecule 2 should have one encounter";
    EXPECT_EQ(moleculeList[2].crossbase[0], 0) << "molecule 2 records pro2 == 0";

    EXPECT_EQ(complexList[0].ncross, 0) << "complex 0 ncross must stay 0";
    EXPECT_EQ(complexList[1].ncross, 0) << "complex 1 ncross must stay 0";
    EXPECT_EQ(complexList[2].ncross, 1) << "complex 2 ncross must be incremented";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is executed inside its own TEST so
// that a failure in one scenario does not stop the remaining scenarios.
// -----------------------------------------------------------------------------
TEST(GetDistanceToSurface, OnSurfaceReturnsZeroDistance) { test_gdts_on_surface_returns_zero_distance(); }
TEST(GetDistanceToSurface, BoxWithinRmax) { test_gdts_box_within_rmax(); }
TEST(GetDistanceToSurface, BoxBeyondRmaxNoSideEffects) { test_gdts_box_beyond_rmax_no_side_effects(); }
TEST(GetDistanceToSurface, BoxIgnoresXY) { test_gdts_box_ignores_xy(); }
TEST(GetDistanceToSurface, SphereInside) { test_gdts_sphere_inside(); }
TEST(GetDistanceToSurface, SphereOutsideUsesAbsoluteValue) { test_gdts_sphere_outside_uses_absolute_value(); }
TEST(GetDistanceToSurface, SphereBeyondRmax) { test_gdts_sphere_beyond_rmax(); }
TEST(GetDistanceToSurface, RepeatedCallsAccumulate) { test_gdts_repeated_calls_accumulate(); }
TEST(GetDistanceToSurface, UsesCorrectIndices) { test_gdts_uses_correct_indices(); }