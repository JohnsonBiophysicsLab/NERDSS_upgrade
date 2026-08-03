/*! \file test_get_distance.cpp
 *
 * ### Unit test for src/reactions/get_distance.cpp
 *
 * The single (active) function in that translation unit is:
 *
 *     bool get_distance(int pro1, int pro2, int iface1, int iface2,
 *                       int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
 *                       double& sep, double& R1, double Rmax,
 *                       std::vector<Complex>& complexList,
 *                       const ForwardRxn& currRxn,
 *                       std::vector<Molecule>& moleculeList, bool isSphere)
 *
 * It measures the separation between two interfaces on two molecules and, if
 * they are closer than Rmax, records the encounter on both molecules
 * (crossbase / mycrossint / crossrxn) and increments the ncross counter of both
 * parent complexes.  There are four distinct geometry branches:
 *
 *   1. 1D  ("onFiber" for both complexes)  -> distance uses only the x-axis
 *   2. spherical 2D (isSphere && both complexes OnSurface) -> geodesic arc length
 *   3. planar   2D (both complexes OnSurface, not a sphere) -> dz forced to zero
 *   4. 3D (default) -> full Cartesian distance
 *
 * Each branch is exercised below, along with the "too far away" early-out and
 * the bookkeeping side effects.  Verbose output describes each scenario and the
 * criteria used for passing.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gtest/gtest.h>

namespace {

/*! \brief A tiny two-molecule / two-complex system used by every test below. */
struct GdSystem {
    std::vector<Molecule> moleculeList; //!< exactly two molecules, indices 0 and 1
    std::vector<Complex> complexList; //!< exactly two complexes, indices 0 and 1
};

/*! \brief Build a molecule with one interface, both COM and interface at \p crd.
 *
 * \param[in] crd    coordinate for the COM and the single interface
 * \param[in] comIdx index of the parent complex in complexList
 */
Molecule gd_make_molecule(const Coord& crd, int comIdx)
{
    Molecule mol;
    mol.comCoord = crd;
    mol.myComIndex = comIdx;

    Molecule::Iface iface;
    iface.coord = crd; // interface sits on top of the COM: keeps geometry trivial
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a complex owning a single member molecule. */
Complex gd_make_complex(const Coord& crd, int memberMolIndex)
{
    Complex com;
    com.comCoord = crd;
    com.memberList.clear();
    com.memberList.push_back(memberMolIndex);
    com.ncross = 0; // encounters counted by get_distance start at zero
    com.OnSurface = false; // default: free 3D diffusion
    com.onFiber = false; // default: not a 1D fiber-bound complex
    return com;
}

/*! \brief Assemble a complete two-body system with molecule 0 at \p c1 and 1 at \p c2. */
GdSystem gd_make_system(const Coord& c1, const Coord& c2)
{
    GdSystem sys;
    sys.moleculeList.push_back(gd_make_molecule(c1, 0));
    sys.moleculeList.push_back(gd_make_molecule(c2, 1));
    sys.complexList.push_back(gd_make_complex(c1, 0));
    sys.complexList.push_back(gd_make_complex(c2, 1));
    return sys;
}

/*! \brief A minimal reaction carrying only the binding radius used by get_distance. */
ForwardRxn gd_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.bindRadius = bindRadius;
    return rxn;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: plain 3D geometry, molecules inside Rmax.
// -----------------------------------------------------------------------------
void test_gd_3d_within_rmax()
{
    std::cerr << "\n[TEST] test_gd_3d_within_rmax\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (3D branch)\n"
              << "  Scenario:      interfaces at (0,0,0) and (3,4,0), bindRadius = 1,\n"
              << "                 Rmax = 20 so the pair is a valid encounter.\n"
              << "  Pass criteria: returns true, R1 == 5, sep == 5 - bindRadius,\n"
              << "                 and both molecules/complexes record the encounter.\n";

    GdSystem sys = gd_make_system(Coord{ 0.0, 0.0, 0.0 }, Coord{ 3.0, 4.0, 0.0 });
    ForwardRxn rxn = gd_make_rxn(1.0);

    double sep = -1.0;
    double R1 = -1.0;

    const bool didEncounter = get_distance(/*pro1*/ 0, /*pro2*/ 1, /*iface1*/ 0, /*iface2*/ 0,
        /*rxnIndex*/ 7, /*rateIndex*/ 2, /*isStateChangeBackRxn*/ false, sep, R1,
        /*Rmax*/ 20.0, sys.complexList, rxn, sys.moleculeList, /*isSphere*/ false);

    std::cerr << "  Returned " << std::boolalpha << didEncounter
              << ", R1 = " << R1 << ", sep = " << sep << '\n';

    // Geometry: sqrt(3^2 + 4^2) == 5.
    EXPECT_TRUE(didEncounter) << "Pair inside Rmax must be flagged as an encounter";
    EXPECT_NEAR(R1, 5.0, 1e-12) << "3D distance between (0,0,0) and (3,4,0) must be 5";
    EXPECT_NEAR(sep, 5.0 - rxn.bindRadius, 1e-12) << "sep must be R1 minus the binding radius";

    // Bookkeeping on molecule 0.
    EXPECT_EQ(sys.moleculeList[0].crossbase.size(), 1u) << "molecule 0 should record one partner";
    if (!sys.moleculeList[0].crossbase.empty())
        EXPECT_EQ(sys.moleculeList[0].crossbase[0], 1) << "molecule 0's partner should be molecule 1";
    EXPECT_EQ(sys.moleculeList[0].mycrossint.size(), 1u) << "molecule 0 should record one interface";

    // Bookkeeping on molecule 1.
    EXPECT_EQ(sys.moleculeList[1].crossbase.size(), 1u) << "molecule 1 should record one partner";
    if (!sys.moleculeList[1].crossbase.empty())
        EXPECT_EQ(sys.moleculeList[1].crossbase[0], 0) << "molecule 1's partner should be molecule 0";

    // Both parent complexes must have their encounter counter incremented.
    EXPECT_EQ(sys.complexList[0].ncross, 1) << "complex 0 ncross should be incremented";
    EXPECT_EQ(sys.complexList[1].ncross, 1) << "complex 1 ncross should be incremented";
}

// -----------------------------------------------------------------------------
// Test 2: plain 3D geometry, molecules farther apart than Rmax.
// -----------------------------------------------------------------------------
void test_gd_3d_outside_rmax()
{
    std::cerr << "\n[TEST] test_gd_3d_outside_rmax\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (3D branch, early-out)\n"
              << "  Scenario:      interfaces separated by 100 nm but Rmax is only 10.\n"
              << "  Pass criteria: returns false, R1/sep still filled in, and no\n"
              << "                 encounter data is appended to either molecule.\n";

    GdSystem sys = gd_make_system(Coord{ 0.0, 0.0, 0.0 }, Coord{ 100.0, 0.0, 0.0 });
    ForwardRxn rxn = gd_make_rxn(1.0);

    double sep = -1.0;
    double R1 = -1.0;

    const bool didEncounter = get_distance(0, 1, 0, 0, 3, 0, false, sep, R1,
        /*Rmax*/ 10.0, sys.complexList, rxn, sys.moleculeList, /*isSphere*/ false);

    std::cerr << "  Returned " << std::boolalpha << didEncounter
              << ", R1 = " << R1 << ", sep = " << sep << '\n';

    EXPECT_FALSE(didEncounter) << "Pair outside Rmax must not be flagged as an encounter";
    EXPECT_NEAR(R1, 100.0, 1e-12) << "R1 should still hold the measured distance";
    EXPECT_NEAR(sep, 99.0, 1e-12) << "sep should still hold R1 - bindRadius";

    // No side effects when the pair is too far apart.
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty()) << "molecule 0 must not record a partner";
    EXPECT_TRUE(sys.moleculeList[1].crossbase.empty()) << "molecule 1 must not record a partner";
    EXPECT_TRUE(sys.moleculeList[0].crossrxn.empty()) << "molecule 0 must not record a reaction";
    EXPECT_EQ(sys.complexList[0].ncross, 0) << "complex 0 ncross must stay at zero";
    EXPECT_EQ(sys.complexList[1].ncross, 0) << "complex 1 ncross must stay at zero";
}

// -----------------------------------------------------------------------------
// Test 3: planar 2D geometry (both complexes OnSurface) ignores the z offset.
// -----------------------------------------------------------------------------
void test_gd_2d_ignores_z()
{
    std::cerr << "\n[TEST] test_gd_2d_ignores_z\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (planar 2D branch)\n"
              << "  Scenario:      both complexes flagged OnSurface; interfaces are\n"
              << "                 (0,0,0) and (3,4,50) - a large z offset.\n"
              << "  Pass criteria: dz is forced to zero, so R1 == 5 (not sqrt(25+2500)).\n";

    GdSystem sys = gd_make_system(Coord{ 0.0, 0.0, 0.0 }, Coord{ 3.0, 4.0, 50.0 });
    // Mark both complexes as membrane bound -> triggers the 2D branch.
    sys.complexList[0].OnSurface = true;
    sys.complexList[1].OnSurface = true;

    ForwardRxn rxn = gd_make_rxn(2.0);

    double sep = -1.0;
    double R1 = -1.0;

    const bool didEncounter = get_distance(0, 1, 0, 0, 1, 0, false, sep, R1,
        /*Rmax*/ 20.0, sys.complexList, rxn, sys.moleculeList, /*isSphere*/ false);

    std::cerr << "  Returned " << std::boolalpha << didEncounter
              << ", R1 = " << R1 << " (3D distance would have been "
              << std::sqrt(3.0 * 3.0 + 4.0 * 4.0 + 50.0 * 50.0) << ")\n";

    EXPECT_TRUE(didEncounter) << "In-plane distance of 5 is inside Rmax = 20";
    EXPECT_NEAR(R1, 5.0, 1e-12) << "2D branch must drop the z component entirely";
    EXPECT_NEAR(sep, 5.0 - 2.0, 1e-12) << "sep must be the in-plane R1 minus bindRadius";
}

// -----------------------------------------------------------------------------
// Test 4: spherical 2D geometry uses the geodesic arc length.
// -----------------------------------------------------------------------------
void test_gd_sphere_2d_geodesic()
{
    std::cerr << "\n[TEST] test_gd_sphere_2d_geodesic\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (spherical 2D branch)\n"
              << "  Scenario:      isSphere = true and both complexes OnSurface;\n"
              << "                 interfaces at (50,0,0) and (0,50,0) on a sphere of\n"
              << "                 radius 50 -> subtended angle pi/2.\n"
              << "  Pass criteria: R1 == 50 * pi/2 (the arc length), not the chord.\n";

    const double sphereRadius = 50.0;
    GdSystem sys = gd_make_system(Coord{ sphereRadius, 0.0, 0.0 }, Coord{ 0.0, sphereRadius, 0.0 });
    sys.complexList[0].OnSurface = true;
    sys.complexList[1].OnSurface = true;

    ForwardRxn rxn = gd_make_rxn(1.0);

    double sep = -1.0;
    double R1 = -1.0;

    const bool didEncounter = get_distance(0, 1, 0, 0, 0, 0, false, sep, R1,
        /*Rmax*/ 500.0, sys.complexList, rxn, sys.moleculeList, /*isSphere*/ true);

    const double expectedArc = sphereRadius * (M_PI / 2.0);
    const double chord = sphereRadius * std::sqrt(2.0);

    std::cerr << "  Returned " << std::boolalpha << didEncounter
              << ", R1 = " << R1 << " (expected arc " << expectedArc
              << ", straight-line chord would be " << chord << ")\n";

    EXPECT_TRUE(didEncounter) << "Arc length is well within the generous Rmax";
    EXPECT_NEAR(R1, expectedArc, 1e-9) << "Spherical branch must return the geodesic arc length";
    EXPECT_NEAR(sep, expectedArc - rxn.bindRadius, 1e-9) << "sep must be arc length minus bindRadius";
    // Sanity check that the geodesic really differs from the chord for this configuration.
    EXPECT_GT(R1, chord) << "For a 90-degree separation the arc must exceed the chord";
}

// -----------------------------------------------------------------------------
// Test 5: isSphere = true but the complexes are NOT on the surface -> 3D branch.
// -----------------------------------------------------------------------------
void test_gd_sphere_but_3d_uses_cartesian()
{
    std::cerr << "\n[TEST] test_gd_sphere_but_3d_uses_cartesian\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (isSphere with 3D complexes)\n"
              << "  Scenario:      isSphere = true, but neither complex is OnSurface,\n"
              << "                 so the spherical geodesic branch must NOT be taken.\n"
              << "  Pass criteria: R1 equals the plain Cartesian distance (5).\n";

    GdSystem sys = gd_make_system(Coord{ 0.0, 0.0, 0.0 }, Coord{ 0.0, 3.0, 4.0 });
    ForwardRxn rxn = gd_make_rxn(1.0);

    double sep = -1.0;
    double R1 = -1.0;

    const bool didEncounter = get_distance(0, 1, 0, 0, 0, 0, false, sep, R1,
        /*Rmax*/ 50.0, sys.complexList, rxn, sys.moleculeList, /*isSphere*/ true);

    std::cerr << "  Returned " << std::boolalpha << didEncounter << ", R1 = " << R1 << '\n';

    EXPECT_TRUE(didEncounter) << "Distance of 5 is inside Rmax = 50";
    EXPECT_NEAR(R1, 5.0, 1e-12) << "Non-surface complexes must use the Cartesian distance";
    EXPECT_NEAR(sep, 4.0, 1e-12) << "sep must be the Cartesian R1 minus bindRadius";
}

// -----------------------------------------------------------------------------
// Test 6: 1D fiber geometry uses only the x coordinate.
// -----------------------------------------------------------------------------
void test_gd_1d_fiber_uses_x_only()
{
    std::cerr << "\n[TEST] test_gd_1d_fiber_uses_x_only\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (1D fiber branch)\n"
              << "  Scenario:      both complexes flagged onFiber; interfaces are\n"
              << "                 (10,7,9) and (14,-3,100).\n"
              << "  Pass criteria: R1 == |14 - 10| == 4 (y and z ignored) and\n"
              << "                 sep == R1 - bindRadius since neither is a promoter.\n";

    GdSystem sys = gd_make_system(Coord{ 10.0, 7.0, 9.0 }, Coord{ 14.0, -3.0, 100.0 });
    // Both complexes live on the same 1D fiber.
    sys.complexList[0].onFiber = true;
    sys.complexList[1].onFiber = true;

    ForwardRxn rxn = gd_make_rxn(1.0);

    double sep = -1.0;
    double R1 = -1.0;

    const bool didEncounter = get_distance(0, 1, 0, 0, 4, 1, true, sep, R1,
        /*Rmax*/ 10.0, sys.complexList, rxn, sys.moleculeList, /*isSphere*/ false);

    std::cerr << "  Returned " << std::boolalpha << didEncounter
              << ", R1 = " << R1 << ", sep = " << sep << '\n';

    EXPECT_TRUE(didEncounter) << "An axial separation of 4 is inside Rmax = 10";
    EXPECT_NEAR(R1, 4.0, 1e-12) << "1D branch must use only the x-axis separation";
    EXPECT_NEAR(sep, 4.0 - rxn.bindRadius, 1e-12)
        << "Non-promoter pairs on a fiber get sep = R1 - bindRadius";

    // The reaction indices passed in must be forwarded verbatim into crossrxn.
    ASSERT_EQ(sys.moleculeList[0].crossrxn.size(), 1u) << "molecule 0 should record one reaction";
    EXPECT_EQ(sys.moleculeList[0].crossrxn[0][0], 4) << "crossrxn[0] should store rxnIndex = 4";
    EXPECT_EQ(sys.moleculeList[0].crossrxn[0][1], 1) << "crossrxn[1] should store rateIndex = 1";
    EXPECT_EQ(sys.moleculeList[0].crossrxn[0][2], 1)
        << "crossrxn[2] should store isStateChangeBackRxn (true -> 1)";
}

// -----------------------------------------------------------------------------
// Test 7: 1D geometry where exactly one partner is a promoter.
// -----------------------------------------------------------------------------
void test_gd_1d_promoter_skips_bindradius()
{
    std::cerr << "\n[TEST] test_gd_1d_promoter_skips_bindradius\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (1D fiber branch, promoter special case)\n"
              << "  Scenario:      both complexes onFiber, molecule 0 is a promoter and\n"
              << "                 molecule 1 is not; axial separation is 4.\n"
              << "  Pass criteria: sep == R1 (the binding radius is NOT subtracted), and\n"
              << "                 the same holds when the roles are swapped.\n";

    ForwardRxn rxn = gd_make_rxn(1.0);

    // --- Case A: molecule 0 is the promoter -------------------------------------
    {
        GdSystem sys = gd_make_system(Coord{ 10.0, 0.0, 0.0 }, Coord{ 14.0, 0.0, 0.0 });
        sys.complexList[0].onFiber = true;
        sys.complexList[1].onFiber = true;
        sys.moleculeList[0].isPromoter = true; // co-localized promoter
        sys.moleculeList[1].isPromoter = false;

        double sep = -1.0;
        double R1 = -1.0;
        const bool didEncounter = get_distance(0, 1, 0, 0, 0, 0, false, sep, R1,
            10.0, sys.complexList, rxn, sys.moleculeList, false);

        std::cerr << "  Case A (mol 0 promoter): R1 = " << R1 << ", sep = " << sep << '\n';
        EXPECT_TRUE(didEncounter) << "Separation of 4 is inside Rmax = 10";
        EXPECT_NEAR(R1, 4.0, 1e-12) << "R1 is still the axial separation";
        EXPECT_NEAR(sep, R1, 1e-12) << "With one promoter, sep must equal R1 exactly";
    }

    // --- Case B: molecule 1 is the promoter (symmetric branch) ------------------
    {
        GdSystem sys = gd_make_system(Coord{ 10.0, 0.0, 0.0 }, Coord{ 14.0, 0.0, 0.0 });
        sys.complexList[0].onFiber = true;
        sys.complexList[1].onFiber = true;
        sys.moleculeList[0].isPromoter = false;
        sys.moleculeList[1].isPromoter = true;

        double sep = -1.0;
        double R1 = -1.0;
        const bool didEncounter = get_distance(0, 1, 0, 0, 0, 0, false, sep, R1,
            10.0, sys.complexList, rxn, sys.moleculeList, false);

        std::cerr << "  Case B (mol 1 promoter): R1 = " << R1 << ", sep = " << sep << '\n';
        EXPECT_TRUE(didEncounter) << "Separation of 4 is inside Rmax = 10";
        EXPECT_NEAR(sep, R1, 1e-12) << "Symmetric promoter branch must also give sep == R1";
    }

    // --- Case C: both are promoters -> normal 1D rule applies -------------------
    {
        GdSystem sys = gd_make_system(Coord{ 10.0, 0.0, 0.0 }, Coord{ 14.0, 0.0, 0.0 });
        sys.complexList[0].onFiber = true;
        sys.complexList[1].onFiber = true;
        sys.moleculeList[0].isPromoter = true;
        sys.moleculeList[1].isPromoter = true;

        double sep = -1.0;
        double R1 = -1.0;
        get_distance(0, 1, 0, 0, 0, 0, false, sep, R1, 10.0, sys.complexList, rxn,
            sys.moleculeList, false);

        std::cerr << "  Case C (both promoters): R1 = " << R1 << ", sep = " << sep << '\n';
        EXPECT_NEAR(sep, R1 - rxn.bindRadius, 1e-12)
            << "When both are promoters the binding radius is subtracted again";
    }
}

// -----------------------------------------------------------------------------
// Test 8: repeated calls accumulate encounters instead of overwriting them.
// -----------------------------------------------------------------------------
void test_gd_accumulates_multiple_encounters()
{
    std::cerr << "\n[TEST] test_gd_accumulates_multiple_encounters\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (bookkeeping side effects)\n"
              << "  Scenario:      call get_distance twice for the same close pair.\n"
              << "  Pass criteria: crossbase/mycrossint/crossrxn each grow to length 2\n"
              << "                 and both complexes report ncross == 2.\n";

    GdSystem sys = gd_make_system(Coord{ 0.0, 0.0, 0.0 }, Coord{ 1.0, 0.0, 0.0 });
    ForwardRxn rxn = gd_make_rxn(0.5);

    double sep = 0.0;
    double R1 = 0.0;

    // Two successive evaluations, e.g. two different reactions for the same pair.
    const bool first = get_distance(0, 1, 0, 0, 11, 0, false, sep, R1, 10.0,
        sys.complexList, rxn, sys.moleculeList, false);
    const bool second = get_distance(0, 1, 0, 0, 12, 1, false, sep, R1, 10.0,
        sys.complexList, rxn, sys.moleculeList, false);

    std::cerr << "  First call returned " << std::boolalpha << first
              << ", second call returned " << second << '\n';
    std::cerr << "  molecule 0 crossbase size = " << sys.moleculeList[0].crossbase.size()
              << ", complex 0 ncross = " << sys.complexList[0].ncross << '\n';

    EXPECT_TRUE(first) << "First close-pair evaluation should register an encounter";
    EXPECT_TRUE(second) << "Second close-pair evaluation should register an encounter";

    EXPECT_EQ(sys.moleculeList[0].crossbase.size(), 2u) << "crossbase should accumulate two entries";
    EXPECT_EQ(sys.moleculeList[0].mycrossint.size(), 2u) << "mycrossint should accumulate two entries";
    ASSERT_EQ(sys.moleculeList[0].crossrxn.size(), 2u) << "crossrxn should accumulate two entries";
    EXPECT_EQ(sys.moleculeList[0].crossrxn[0][0], 11) << "first recorded rxnIndex should be 11";
    EXPECT_EQ(sys.moleculeList[0].crossrxn[1][0], 12) << "second recorded rxnIndex should be 12";

    EXPECT_EQ(sys.complexList[0].ncross, 2) << "complex 0 should count two encounters";
    EXPECT_EQ(sys.complexList[1].ncross, 2) << "complex 1 should count two encounters";
}

// -----------------------------------------------------------------------------
// Test 9: a pair exactly at the binding radius yields sep == 0.
// -----------------------------------------------------------------------------
void test_gd_sep_zero_at_bindradius()
{
    std::cerr << "\n[TEST] test_gd_sep_zero_at_bindradius\n"
              << "  Source file:   src/reactions/get_distance.cpp\n"
              << "  Function:      get_distance (3D branch, boundary value)\n"
              << "  Scenario:      interfaces exactly bindRadius apart (2.0 nm).\n"
              << "  Pass criteria: returns true and sep is exactly 0.\n";

    GdSystem sys = gd_make_system(Coord{ 0.0, 0.0, 0.0 }, Coord{ 2.0, 0.0, 0.0 });
    ForwardRxn rxn = gd_make_rxn(2.0);

    double sep = -1.0;
    double R1 = -1.0;

    const bool didEncounter = get_distance(0, 1, 0, 0, 0, 0, false, sep, R1, 10.0,
        sys.complexList, rxn, sys.moleculeList, false);

    std::cerr << "  Returned " << std::boolalpha << didEncounter
              << ", R1 = " << R1 << ", sep = " << sep << '\n';

    EXPECT_TRUE(didEncounter) << "Contact distance is inside Rmax";
    EXPECT_NEAR(R1, 2.0, 1e-12) << "R1 should equal the binding radius here";
    EXPECT_NEAR(sep, 0.0, 1e-12) << "sep should vanish exactly at contact";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario runs in its own TEST so a failure in one
// does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(GetDistance, ThreeDWithinRmax) { test_gd_3d_within_rmax(); }
TEST(GetDistance, ThreeDOutsideRmax) { test_gd_3d_outside_rmax(); }
TEST(GetDistance, TwoDIgnoresZ) { test_gd_2d_ignores_z(); }
TEST(GetDistance, SphereTwoDGeodesic) { test_gd_sphere_2d_geodesic(); }
TEST(GetDistance, SphereButThreeDUsesCartesian) { test_gd_sphere_but_3d_uses_cartesian(); }
TEST(GetDistance, OneDFiberUsesXOnly) { test_gd_1d_fiber_uses_x_only(); }
TEST(GetDistance, OneDPromoterSkipsBindRadius) { test_gd_1d_promoter_skips_bindradius(); }
TEST(GetDistance, AccumulatesMultipleEncounters) { test_gd_accumulates_multiple_encounters(); }
TEST(GetDistance, SepZeroAtBindRadius) { test_gd_sep_zero_at_bindradius(); }