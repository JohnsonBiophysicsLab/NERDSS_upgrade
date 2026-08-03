/*! \file test_com_of_two_tmp_complexes.cpp
 *
 * ### Unit test for src/reactions/com_of_two_tmp_complexes.cpp
 *
 * Function under test:
 *
 *     void com_of_two_tmp_complexes(Complex& reactCom1, Complex& reactCom2,
 *                                  Coord& vectorCOM,
 *                                  std::vector<Molecule>& moleculeList)
 *
 * The routine computes the mass weighted center of mass (COM) of the union of
 * the two supplied complexes, using the *temporary* association coordinates
 * (Molecule::tmpComCoord) rather than the real coordinates.  The result is
 * written into `vectorCOM`:
 *
 *     vectorCOM = sum_i ( m_i * tmpComCoord_i ) / sum_i ( m_i )
 *
 * where i runs over every molecule index listed in reactCom1.memberList
 * followed by every index listed in reactCom2.memberList.
 *
 * The tests below verify:
 *   1. the simple symmetric case (equal masses -> arithmetic midpoint),
 *   2. proper mass weighting (heavy molecule pulls the COM towards itself),
 *   3. that the result depends on tmpComCoord and NOT on comCoord,
 *   4. that any pre-existing value in vectorCOM is fully overwritten,
 *   5. multi-member complexes with several molecules per complex,
 *   6. that molecules not referenced by either memberList are ignored,
 *   7. the degenerate all-zero-mass case (documents the 0/0 behaviour).
 *
 * Verbose progress information is written to stderr so a reader of the test log
 * can see exactly which scenario is being exercised and what the pass criteria
 * are.  Only non-fatal EXPECT_* assertions are used so every test runs to
 * completion even if one check fails.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

namespace {

// -----------------------------------------------------------------------------
// Helpers (file-local, anonymous namespace so no symbol collisions occur).
// -----------------------------------------------------------------------------

/*! \brief Build a minimal Molecule carrying a mass and a temporary COM coord.
 *
 * The function under test only reads Molecule::mass and Molecule::tmpComCoord,
 * so those are the only fields that need to be meaningful.  `comCoord` is
 * deliberately set to a *different* value than `tmpComCoord` so that a test can
 * prove the routine reads the temporary coordinates.
 *
 * \param[in] tmpCom  value assigned to Molecule::tmpComCoord
 * \param[in] mass    value assigned to Molecule::mass
 * \param[in] realCom value assigned to Molecule::comCoord (decoy)
 */
Molecule cotc_make_molecule(const Coord& tmpCom, double mass, const Coord& realCom)
{
    Molecule mol;
    mol.mass = mass;
    mol.tmpComCoord = tmpCom;
    mol.comCoord = realCom; // must NOT influence the computed COM
    return mol;
}

/*! \brief Build a Complex whose memberList is exactly the supplied indices. */
Complex cotc_make_complex(const std::vector<int>& memberIndices)
{
    Complex com;
    com.memberList = memberIndices;
    return com;
}

/*! \brief Convenience printer used by the verbose logging below. */
void cotc_print_coord(const char* label, const Coord& c)
{
    std::cerr << "    " << label << " = (" << c.x << ", " << c.y << ", " << c.z << ")\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: two single-molecule complexes of equal mass.
//         The COM must be the exact midpoint of the two tmp coordinates.
// -----------------------------------------------------------------------------
void test_cotc_equal_masses_midpoint()
{
    std::cerr << "\n[TEST] test_cotc_equal_masses_midpoint\n"
              << "  Source file:   src/reactions/com_of_two_tmp_complexes.cpp\n"
              << "  Function:      com_of_two_tmp_complexes()\n"
              << "  Scenario:      two complexes, one molecule each, identical mass\n"
              << "                 at (0,0,0) and (10,20,30).\n"
              << "  Pass criteria: vectorCOM equals the arithmetic midpoint\n"
              << "                 (5,10,15).\n";

    // Two molecules of identical mass 1.0.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cotc_make_molecule(Coord{ 0.0, 0.0, 0.0 }, 1.0, Coord{ 999.0, 999.0, 999.0 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ 10.0, 20.0, 30.0 }, 1.0, Coord{ -999.0, -999.0, -999.0 }));

    // Each complex owns exactly one of them.
    Complex reactCom1 = cotc_make_complex({ 0 });
    Complex reactCom2 = cotc_make_complex({ 1 });

    Coord vectorCOM;
    std::cerr << "  Calling com_of_two_tmp_complexes...\n";
    com_of_two_tmp_complexes(reactCom1, reactCom2, vectorCOM, moleculeList);
    cotc_print_coord("vectorCOM", vectorCOM);

    // Midpoint of (0,0,0) and (10,20,30).
    EXPECT_DOUBLE_EQ(vectorCOM.x, 5.0) << "x of the COM should be the midpoint 5.0";
    EXPECT_DOUBLE_EQ(vectorCOM.y, 10.0) << "y of the COM should be the midpoint 10.0";
    EXPECT_DOUBLE_EQ(vectorCOM.z, 15.0) << "z of the COM should be the midpoint 15.0";
}

// -----------------------------------------------------------------------------
// Test 2: unequal masses -> the COM must be weighted towards the heavy body.
// -----------------------------------------------------------------------------
void test_cotc_mass_weighting()
{
    std::cerr << "\n[TEST] test_cotc_mass_weighting\n"
              << "  Source file:   src/reactions/com_of_two_tmp_complexes.cpp\n"
              << "  Function:      com_of_two_tmp_complexes()\n"
              << "  Scenario:      mass 3 at x=0 and mass 1 at x=8.\n"
              << "  Pass criteria: vectorCOM.x == (3*0 + 1*8)/4 == 2.0, i.e. the\n"
              << "                 COM is pulled towards the heavier molecule.\n";

    std::vector<Molecule> moleculeList;
    // Heavy molecule (mass 3) sitting at the origin.
    moleculeList.push_back(cotc_make_molecule(Coord{ 0.0, 0.0, 0.0 }, 3.0, Coord{ 0.0, 0.0, 0.0 }));
    // Light molecule (mass 1) sitting at x = 8.
    moleculeList.push_back(cotc_make_molecule(Coord{ 8.0, 4.0, -4.0 }, 1.0, Coord{ 0.0, 0.0, 0.0 }));

    Complex reactCom1 = cotc_make_complex({ 0 });
    Complex reactCom2 = cotc_make_complex({ 1 });

    Coord vectorCOM;
    std::cerr << "  Calling com_of_two_tmp_complexes...\n";
    com_of_two_tmp_complexes(reactCom1, reactCom2, vectorCOM, moleculeList);
    cotc_print_coord("vectorCOM", vectorCOM);

    // Expected: (3*0 + 1*8)/4 = 2, (3*0 + 1*4)/4 = 1, (3*0 + 1*-4)/4 = -1.
    EXPECT_DOUBLE_EQ(vectorCOM.x, 2.0) << "mass-weighted x should be 2.0";
    EXPECT_DOUBLE_EQ(vectorCOM.y, 1.0) << "mass-weighted y should be 1.0";
    EXPECT_DOUBLE_EQ(vectorCOM.z, -1.0) << "mass-weighted z should be -1.0";

    // Sanity check: heavier side must dominate, so the COM is closer to x=0.
    EXPECT_LT(vectorCOM.x, 4.0) << "COM must lie closer to the heavy molecule than the midpoint";
}

// -----------------------------------------------------------------------------
// Test 3: the routine must use tmpComCoord, never comCoord.
// -----------------------------------------------------------------------------
void test_cotc_uses_tmp_coords_only()
{
    std::cerr << "\n[TEST] test_cotc_uses_tmp_coords_only\n"
              << "  Source file:   src/reactions/com_of_two_tmp_complexes.cpp\n"
              << "  Function:      com_of_two_tmp_complexes()\n"
              << "  Scenario:      tmpComCoord values are (1,1,1) and (3,3,3) while\n"
              << "                 the real comCoord values are absurd (1e6).\n"
              << "  Pass criteria: vectorCOM == (2,2,2), showing the real\n"
              << "                 coordinates were ignored.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cotc_make_molecule(Coord{ 1.0, 1.0, 1.0 }, 1.0, Coord{ 1.0e6, 1.0e6, 1.0e6 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ 3.0, 3.0, 3.0 }, 1.0, Coord{ -1.0e6, -1.0e6, -1.0e6 }));

    Complex reactCom1 = cotc_make_complex({ 0 });
    Complex reactCom2 = cotc_make_complex({ 1 });

    Coord vectorCOM;
    std::cerr << "  Calling com_of_two_tmp_complexes...\n";
    com_of_two_tmp_complexes(reactCom1, reactCom2, vectorCOM, moleculeList);
    cotc_print_coord("vectorCOM", vectorCOM);

    EXPECT_DOUBLE_EQ(vectorCOM.x, 2.0) << "x must come from tmpComCoord, not comCoord";
    EXPECT_DOUBLE_EQ(vectorCOM.y, 2.0) << "y must come from tmpComCoord, not comCoord";
    EXPECT_DOUBLE_EQ(vectorCOM.z, 2.0) << "z must come from tmpComCoord, not comCoord";

    // Guard against the (wrong) use of the huge real coordinates.
    EXPECT_LT(std::fabs(vectorCOM.x), 100.0) << "COM magnitude indicates comCoord leaked into the sum";
}

// -----------------------------------------------------------------------------
// Test 4: any pre-existing content of vectorCOM must be discarded.
// -----------------------------------------------------------------------------
void test_cotc_overwrites_output_argument()
{
    std::cerr << "\n[TEST] test_cotc_overwrites_output_argument\n"
              << "  Source file:   src/reactions/com_of_two_tmp_complexes.cpp\n"
              << "  Function:      com_of_two_tmp_complexes()\n"
              << "  Scenario:      vectorCOM is pre-loaded with garbage (-77,88,-99)\n"
              << "                 before the call; both molecules sit at (2,4,6).\n"
              << "  Pass criteria: vectorCOM == (2,4,6), i.e. the routine zeroes\n"
              << "                 the accumulator before summing.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cotc_make_molecule(Coord{ 2.0, 4.0, 6.0 }, 2.0, Coord{ 0.0, 0.0, 0.0 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ 2.0, 4.0, 6.0 }, 5.0, Coord{ 0.0, 0.0, 0.0 }));

    Complex reactCom1 = cotc_make_complex({ 0 });
    Complex reactCom2 = cotc_make_complex({ 1 });

    // Deliberately poison the output argument.
    Coord vectorCOM{ -77.0, 88.0, -99.0 };
    cotc_print_coord("vectorCOM before call", vectorCOM);

    std::cerr << "  Calling com_of_two_tmp_complexes...\n";
    com_of_two_tmp_complexes(reactCom1, reactCom2, vectorCOM, moleculeList);
    cotc_print_coord("vectorCOM after call", vectorCOM);

    // Both molecules share the same position, so the COM is that position
    // regardless of the masses used.
    EXPECT_DOUBLE_EQ(vectorCOM.x, 2.0) << "stale x value was not cleared";
    EXPECT_DOUBLE_EQ(vectorCOM.y, 4.0) << "stale y value was not cleared";
    EXPECT_DOUBLE_EQ(vectorCOM.z, 6.0) << "stale z value was not cleared";
}

// -----------------------------------------------------------------------------
// Test 5: multi-member complexes (2 molecules + 3 molecules).
// -----------------------------------------------------------------------------
void test_cotc_multi_member_complexes()
{
    std::cerr << "\n[TEST] test_cotc_multi_member_complexes\n"
              << "  Source file:   src/reactions/com_of_two_tmp_complexes.cpp\n"
              << "  Function:      com_of_two_tmp_complexes()\n"
              << "  Scenario:      complex 1 has 2 members, complex 2 has 3 members,\n"
              << "                 all with different masses and positions.\n"
              << "  Pass criteria: vectorCOM matches the hand computed mass\n"
              << "                 weighted average over all five molecules.\n";

    // Positions and masses chosen so the expected answer is easy to verify.
    // idx : tmpComCoord            mass
    //  0  : ( 0,  0,  0)            1
    //  1  : ( 2,  0,  0)            1
    //  2  : ( 0,  4,  0)            2
    //  3  : ( 0,  0,  6)            2
    //  4  : (-4, -4, -4)            4
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cotc_make_molecule(Coord{ 0.0, 0.0, 0.0 }, 1.0, Coord{ 0.0, 0.0, 0.0 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ 2.0, 0.0, 0.0 }, 1.0, Coord{ 0.0, 0.0, 0.0 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ 0.0, 4.0, 0.0 }, 2.0, Coord{ 0.0, 0.0, 0.0 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ 0.0, 0.0, 6.0 }, 2.0, Coord{ 0.0, 0.0, 0.0 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ -4.0, -4.0, -4.0 }, 4.0, Coord{ 0.0, 0.0, 0.0 }));

    Complex reactCom1 = cotc_make_complex({ 0, 1 });
    Complex reactCom2 = cotc_make_complex({ 2, 3, 4 });

    // Reference computation performed independently of the routine.
    double totalMass = 0.0;
    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    for (int idx : { 0, 1, 2, 3, 4 }) {
        totalMass += moleculeList[idx].mass;
        sx += moleculeList[idx].tmpComCoord.x * moleculeList[idx].mass;
        sy += moleculeList[idx].tmpComCoord.y * moleculeList[idx].mass;
        sz += moleculeList[idx].tmpComCoord.z * moleculeList[idx].mass;
    }
    const double expectedX = sx / totalMass;
    const double expectedY = sy / totalMass;
    const double expectedZ = sz / totalMass;

    Coord vectorCOM;
    std::cerr << "  Calling com_of_two_tmp_complexes...\n";
    com_of_two_tmp_complexes(reactCom1, reactCom2, vectorCOM, moleculeList);
    cotc_print_coord("vectorCOM", vectorCOM);
    std::cerr << "    expected  = (" << expectedX << ", " << expectedY << ", " << expectedZ << ")\n";

    EXPECT_NEAR(vectorCOM.x, expectedX, 1e-12) << "x should match the reference mass-weighted mean";
    EXPECT_NEAR(vectorCOM.y, expectedY, 1e-12) << "y should match the reference mass-weighted mean";
    EXPECT_NEAR(vectorCOM.z, expectedZ, 1e-12) << "z should match the reference mass-weighted mean";
}

// -----------------------------------------------------------------------------
// Test 6: molecules that are not members of either complex must be ignored.
// -----------------------------------------------------------------------------
void test_cotc_ignores_non_member_molecules()
{
    std::cerr << "\n[TEST] test_cotc_ignores_non_member_molecules\n"
              << "  Source file:   src/reactions/com_of_two_tmp_complexes.cpp\n"
              << "  Function:      com_of_two_tmp_complexes()\n"
              << "  Scenario:      moleculeList contains a very heavy 'spectator'\n"
              << "                 molecule that appears in neither memberList.\n"
              << "  Pass criteria: vectorCOM is unaffected by the spectator, i.e.\n"
              << "                 it equals the midpoint of the two members.\n";

    std::vector<Molecule> moleculeList;
    // index 0: spectator - huge mass, far away, listed in NO memberList.
    moleculeList.push_back(cotc_make_molecule(Coord{ 1000.0, 1000.0, 1000.0 }, 1000.0, Coord{ 0.0, 0.0, 0.0 }));
    // index 1 and 2: the actual reactants (equal mass).
    moleculeList.push_back(cotc_make_molecule(Coord{ -1.0, -2.0, -3.0 }, 1.0, Coord{ 0.0, 0.0, 0.0 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ 1.0, 2.0, 3.0 }, 1.0, Coord{ 0.0, 0.0, 0.0 }));

    Complex reactCom1 = cotc_make_complex({ 1 });
    Complex reactCom2 = cotc_make_complex({ 2 });

    Coord vectorCOM;
    std::cerr << "  Calling com_of_two_tmp_complexes...\n";
    com_of_two_tmp_complexes(reactCom1, reactCom2, vectorCOM, moleculeList);
    cotc_print_coord("vectorCOM", vectorCOM);

    // Symmetric reactants about the origin -> COM at the origin.
    EXPECT_DOUBLE_EQ(vectorCOM.x, 0.0) << "spectator molecule must not contribute to x";
    EXPECT_DOUBLE_EQ(vectorCOM.y, 0.0) << "spectator molecule must not contribute to y";
    EXPECT_DOUBLE_EQ(vectorCOM.z, 0.0) << "spectator molecule must not contribute to z";
}

// -----------------------------------------------------------------------------
// Test 7: degenerate case, every member has zero mass.
//         totalmass becomes 0 so the final division is 0/0.  This test simply
//         documents/locks in the observed (non-finite) behaviour rather than
//         claiming a numeric answer.
// -----------------------------------------------------------------------------
void test_cotc_zero_total_mass_is_not_finite()
{
    std::cerr << "\n[TEST] test_cotc_zero_total_mass_is_not_finite\n"
              << "  Source file:   src/reactions/com_of_two_tmp_complexes.cpp\n"
              << "  Function:      com_of_two_tmp_complexes()\n"
              << "  Scenario:      all member molecules have mass 0, so the routine\n"
              << "                 divides a zero numerator by a zero total mass.\n"
              << "  Pass criteria: the resulting components are NOT finite (NaN),\n"
              << "                 documenting that callers must guarantee a\n"
              << "                 non-zero total mass.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cotc_make_molecule(Coord{ 1.0, 2.0, 3.0 }, 0.0, Coord{ 0.0, 0.0, 0.0 }));
    moleculeList.push_back(cotc_make_molecule(Coord{ 4.0, 5.0, 6.0 }, 0.0, Coord{ 0.0, 0.0, 0.0 }));

    Complex reactCom1 = cotc_make_complex({ 0 });
    Complex reactCom2 = cotc_make_complex({ 1 });

    Coord vectorCOM;
    std::cerr << "  Calling com_of_two_tmp_complexes...\n";
    com_of_two_tmp_complexes(reactCom1, reactCom2, vectorCOM, moleculeList);
    cotc_print_coord("vectorCOM", vectorCOM);

    // 0.0 / 0.0 -> NaN on any IEEE-754 platform.
    EXPECT_FALSE(std::isfinite(vectorCOM.x)) << "x should be non-finite when total mass is zero";
    EXPECT_FALSE(std::isfinite(vectorCOM.y)) << "y should be non-finite when total mass is zero";
    EXPECT_FALSE(std::isfinite(vectorCOM.z)) << "z should be non-finite when total mass is zero";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* helper is invoked from its own TEST so
// the framework reports individual pass/fail results while still running every
// scenario (all assertions above are non-fatal EXPECT_*).
// -----------------------------------------------------------------------------
TEST(ComOfTwoTmpComplexes, EqualMassesMidpoint) { test_cotc_equal_masses_midpoint(); }
TEST(ComOfTwoTmpComplexes, MassWeighting) { test_cotc_mass_weighting(); }
TEST(ComOfTwoTmpComplexes, UsesTmpCoordsOnly) { test_cotc_uses_tmp_coords_only(); }
TEST(ComOfTwoTmpComplexes, OverwritesOutputArgument) { test_cotc_overwrites_output_argument(); }
TEST(ComOfTwoTmpComplexes, MultiMemberComplexes) { test_cotc_multi_member_complexes(); }
TEST(ComOfTwoTmpComplexes, IgnoresNonMemberMolecules) { test_cotc_ignores_non_member_molecules(); }
TEST(ComOfTwoTmpComplexes, ZeroTotalMassIsNotFinite) { test_cotc_zero_total_mass_is_not_finite(); }