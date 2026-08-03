/*! \file test_conservedMags.cpp
 *
 * ### Unit test for src/reactions/conservedMags.cpp
 *
 * Function under test:
 *
 *     bool conservedMags(const Complex& targCom,
 *                        const std::vector<Molecule>& moleculeList)
 *
 * What the function does
 * ----------------------
 * During association NERDSS keeps two sets of coordinates for every Molecule:
 *   - the "real" coordinates (Molecule::comCoord and Molecule::interfaceList[i].coord)
 *   - the temporary association coordinates (Molecule::tmpComCoord and
 *     Molecule::tmpICoords[i])
 *
 * A rigid-body association may translate and rotate a Complex, but it must never
 * change the length of any interface-to-center-of-mass (IFACE-COM) vector.
 * conservedMags() walks every member Molecule of the Complex and compares the
 * magnitude of each real IFACE-COM vector with the magnitude of the
 * corresponding temporary IFACE-COM vector.  The comparison is done through
 * roundv(), i.e. the magnitudes only have to agree to 4 decimal places.
 *
 * It returns:
 *   - true  if every IFACE-COM magnitude was conserved
 *   - false as soon as one magnitude differs (early return)
 *
 * The tests below build tiny Complex/Molecule structures by hand so we can
 * control the geometry exactly, and verify:
 *   1. pure translation                            -> conserved (true)
 *   2. pure rotation                               -> conserved (true)
 *   3. a stretched interface                       -> not conserved (false)
 *   4. differences smaller than the roundv() cutoff -> still conserved (true)
 *   5. differences larger than the roundv() cutoff  -> not conserved (false)
 *   6. only the second molecule is bad             -> false (loop covers all members)
 *   7. no tmpICoords / no members                  -> vacuously true
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (uniquely prefixed with "cm_" so they cannot collide with other
// translation units in the shared gtest binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a Molecule with matching real and temporary coordinate sets.
 *
 * \param[in] com         real center of mass
 * \param[in] ifaceCoords real interface coordinates (absolute)
 * \param[in] tmpCom      temporary (association) center of mass
 * \param[in] tmpICoords  temporary (association) interface coordinates (absolute)
 *
 * The two coordinate lists are expected to be the same length; conservedMags()
 * iterates over tmpICoords and indexes interfaceList with the same counter.
 */
Molecule cm_make_molecule(const Coord& com, const std::vector<Coord>& ifaceCoords,
    const Coord& tmpCom, const std::vector<Coord>& tmpICoords)
{
    Molecule mol;
    mol.comCoord = com;
    mol.tmpComCoord = tmpCom;

    // Real interfaces.
    mol.interfaceList.clear();
    for (const auto& c : ifaceCoords) {
        Molecule::Iface iface;
        iface.coord = c;
        mol.interfaceList.push_back(iface);
    }

    // Temporary (association) interfaces.
    mol.tmpICoords = tmpICoords;

    return mol;
}

/*! \brief Build a Complex whose memberList points at the given molecule indices. */
Complex cm_make_complex(const std::vector<int>& members, int index = 0)
{
    Complex targCom;
    targCom.index = index;
    targCom.memberList = members;
    return targCom;
}

/*! \brief Convenience: magnitude of (a - b), used only for verbose reporting. */
double cm_mag(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: pure translation of a molecule -> all magnitudes conserved.
// -----------------------------------------------------------------------------
void test_cm_pure_translation_conserved()
{
    std::cerr << "\n[TEST] test_cm_pure_translation_conserved\n"
              << "  Source file:   src/reactions/conservedMags.cpp\n"
              << "  Function:      conservedMags()\n"
              << "  Scenario:      one molecule with two interfaces is shifted by\n"
              << "                 (+10, -5, +2); COM and interfaces move together.\n"
              << "  Pass criteria: conservedMags() returns true.\n";

    // Real geometry: COM at origin, interfaces at (3,0,0) and (0,0,-4).
    const Coord com { 0.0, 0.0, 0.0 };
    const std::vector<Coord> ifaces { Coord { 3.0, 0.0, 0.0 }, Coord { 0.0, 0.0, -4.0 } };

    // Temporary geometry: rigid translation by (10, -5, 2).
    const Coord tmpCom { 10.0, -5.0, 2.0 };
    const std::vector<Coord> tmpIfaces { Coord { 13.0, -5.0, 2.0 }, Coord { 10.0, -5.0, -2.0 } };

    std::vector<Molecule> moleculeList { cm_make_molecule(com, ifaces, tmpCom, tmpIfaces) };
    Complex targCom = cm_make_complex({ 0 });

    std::cerr << "  Real magnitudes:  " << cm_mag(ifaces[0], com) << ", "
              << cm_mag(ifaces[1], com) << "\n"
              << "  Tmp  magnitudes:  " << cm_mag(tmpIfaces[0], tmpCom) << ", "
              << cm_mag(tmpIfaces[1], tmpCom) << "\n";

    const bool result = conservedMags(targCom, moleculeList);
    std::cerr << "  conservedMags() returned " << std::boolalpha << result << "\n";

    EXPECT_TRUE(result) << "A pure translation must conserve every IFACE-COM magnitude";
}

// -----------------------------------------------------------------------------
// Test 2: pure rotation of a molecule -> all magnitudes conserved.
// -----------------------------------------------------------------------------
void test_cm_pure_rotation_conserved()
{
    std::cerr << "\n[TEST] test_cm_pure_rotation_conserved\n"
              << "  Source file:   src/reactions/conservedMags.cpp\n"
              << "  Function:      conservedMags()\n"
              << "  Scenario:      the interface vector (3,0,0) is rotated 90 degrees\n"
              << "                 about z to (0,3,0); the COM does not move.\n"
              << "  Pass criteria: conservedMags() returns true (length unchanged).\n";

    const Coord com { 1.0, 1.0, 1.0 };
    const std::vector<Coord> ifaces { Coord { 4.0, 1.0, 1.0 } }; // offset (3,0,0), mag 3

    const Coord tmpCom { 1.0, 1.0, 1.0 };
    const std::vector<Coord> tmpIfaces { Coord { 1.0, 4.0, 1.0 } }; // offset (0,3,0), mag 3

    std::vector<Molecule> moleculeList { cm_make_molecule(com, ifaces, tmpCom, tmpIfaces) };
    Complex targCom = cm_make_complex({ 0 });

    std::cerr << "  Real magnitude: " << cm_mag(ifaces[0], com)
              << "  Tmp magnitude: " << cm_mag(tmpIfaces[0], tmpCom) << "\n";

    const bool result = conservedMags(targCom, moleculeList);
    std::cerr << "  conservedMags() returned " << std::boolalpha << result << "\n";

    EXPECT_TRUE(result) << "A rigid rotation must conserve the IFACE-COM magnitude";
}

// -----------------------------------------------------------------------------
// Test 3: an interface that is stretched -> magnitude not conserved.
// -----------------------------------------------------------------------------
void test_cm_stretched_iface_not_conserved()
{
    std::cerr << "\n[TEST] test_cm_stretched_iface_not_conserved\n"
              << "  Source file:   src/reactions/conservedMags.cpp\n"
              << "  Function:      conservedMags()\n"
              << "  Scenario:      the temporary interface sits 5 nm from the COM\n"
              << "                 while the real one sits 3 nm away.\n"
              << "  Pass criteria: conservedMags() returns false.\n";

    const Coord com { 0.0, 0.0, 0.0 };
    const std::vector<Coord> ifaces { Coord { 3.0, 0.0, 0.0 } }; // mag 3

    const Coord tmpCom { 0.0, 0.0, 0.0 };
    const std::vector<Coord> tmpIfaces { Coord { 5.0, 0.0, 0.0 } }; // mag 5 -> broken

    std::vector<Molecule> moleculeList { cm_make_molecule(com, ifaces, tmpCom, tmpIfaces) };
    Complex targCom = cm_make_complex({ 0 }, 7);

    std::cerr << "  Real magnitude: " << cm_mag(ifaces[0], com)
              << "  Tmp magnitude: " << cm_mag(tmpIfaces[0], tmpCom) << "\n";

    const bool result = conservedMags(targCom, moleculeList);
    std::cerr << "  conservedMags() returned " << std::boolalpha << result << "\n";

    EXPECT_FALSE(result) << "A stretched IFACE-COM vector must be reported as not conserved";
}

// -----------------------------------------------------------------------------
// Test 4: tiny difference below the roundv() resolution -> still "conserved".
// -----------------------------------------------------------------------------
void test_cm_tiny_difference_within_rounding()
{
    std::cerr << "\n[TEST] test_cm_tiny_difference_within_rounding\n"
              << "  Source file:   src/reactions/conservedMags.cpp\n"
              << "  Function:      conservedMags()\n"
              << "  Scenario:      magnitudes differ by 1e-6 nm, far below the\n"
              << "                 4-decimal-place resolution of roundv().\n"
              << "  Pass criteria: conservedMags() returns true (difference rounds away).\n";

    const Coord com { 0.0, 0.0, 0.0 };
    const std::vector<Coord> ifaces { Coord { 3.0, 0.0, 0.0 } }; // mag 3.0

    const Coord tmpCom { 0.0, 0.0, 0.0 };
    const std::vector<Coord> tmpIfaces { Coord { 3.000001, 0.0, 0.0 } }; // mag 3.000001

    std::vector<Molecule> moleculeList { cm_make_molecule(com, ifaces, tmpCom, tmpIfaces) };
    Complex targCom = cm_make_complex({ 0 });

    std::cerr << "  roundv(real) = " << roundv(cm_mag(ifaces[0], com))
              << "  roundv(tmp)  = " << roundv(cm_mag(tmpIfaces[0], tmpCom)) << "\n";

    const bool result = conservedMags(targCom, moleculeList);
    std::cerr << "  conservedMags() returned " << std::boolalpha << result << "\n";

    EXPECT_TRUE(result) << "Sub-1e-4 deviations must be tolerated because roundv() is used";
}

// -----------------------------------------------------------------------------
// Test 5: difference just above the roundv() resolution -> not conserved.
// -----------------------------------------------------------------------------
void test_cm_small_difference_beyond_rounding()
{
    std::cerr << "\n[TEST] test_cm_small_difference_beyond_rounding\n"
              << "  Source file:   src/reactions/conservedMags.cpp\n"
              << "  Function:      conservedMags()\n"
              << "  Scenario:      magnitudes are 3.0 vs 3.001, which roundv() keeps\n"
              << "                 distinct at 4 decimal places.\n"
              << "  Pass criteria: conservedMags() returns false.\n";

    const Coord com { 0.0, 0.0, 0.0 };
    const std::vector<Coord> ifaces { Coord { 3.0, 0.0, 0.0 } }; // mag 3.0000

    const Coord tmpCom { 0.0, 0.0, 0.0 };
    const std::vector<Coord> tmpIfaces { Coord { 3.001, 0.0, 0.0 } }; // mag 3.0010

    std::vector<Molecule> moleculeList { cm_make_molecule(com, ifaces, tmpCom, tmpIfaces) };
    Complex targCom = cm_make_complex({ 0 });

    std::cerr << "  roundv(real) = " << roundv(cm_mag(ifaces[0], com))
              << "  roundv(tmp)  = " << roundv(cm_mag(tmpIfaces[0], tmpCom)) << "\n";

    const bool result = conservedMags(targCom, moleculeList);
    std::cerr << "  conservedMags() returned " << std::boolalpha << result << "\n";

    EXPECT_FALSE(result) << "Deviations resolvable by roundv() must be flagged";
}

// -----------------------------------------------------------------------------
// Test 6: multi-molecule complex where only the *second* molecule is broken.
//         Confirms the loop really visits every member of the complex.
// -----------------------------------------------------------------------------
void test_cm_multi_molecule_second_bad()
{
    std::cerr << "\n[TEST] test_cm_multi_molecule_second_bad\n"
              << "  Source file:   src/reactions/conservedMags.cpp\n"
              << "  Function:      conservedMags()\n"
              << "  Scenario:      a complex with two member molecules; the first is\n"
              << "                 rigid, the second has a shrunken interface vector.\n"
              << "  Pass criteria: conservedMags() returns false (all members checked).\n";

    // Molecule 0: rigid translation, magnitudes fine.
    Molecule mol0 = cm_make_molecule(
        Coord { 0.0, 0.0, 0.0 }, { Coord { 2.0, 0.0, 0.0 } },
        Coord { 5.0, 5.0, 5.0 }, { Coord { 7.0, 5.0, 5.0 } });

    // Molecule 1: interface collapsed from magnitude 4 to magnitude 1.
    Molecule mol1 = cm_make_molecule(
        Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 4.0, 0.0 } },
        Coord { 0.0, 0.0, 0.0 }, { Coord { 0.0, 1.0, 0.0 } });

    std::vector<Molecule> moleculeList { mol0, mol1 };
    Complex targCom = cm_make_complex({ 0, 1 });

    std::cerr << "  Molecule 0 magnitudes: real 2  tmp 2 (ok)\n"
              << "  Molecule 1 magnitudes: real 4  tmp 1 (broken)\n";

    const bool result = conservedMags(targCom, moleculeList);
    std::cerr << "  conservedMags() returned " << std::boolalpha << result << "\n";

    EXPECT_FALSE(result) << "A violation on any member molecule must return false";

    // Sanity companion check: with both molecules rigid it must return true.
    std::cerr << "  Re-testing with molecule 1 repaired (tmp interface back at mag 4)...\n";
    moleculeList[1].tmpICoords[0] = Coord { 0.0, 4.0, 0.0 };
    const bool repaired = conservedMags(targCom, moleculeList);
    std::cerr << "  conservedMags() returned " << std::boolalpha << repaired << "\n";
    EXPECT_TRUE(repaired) << "With both molecules rigid the function must return true";
}

// -----------------------------------------------------------------------------
// Test 7: degenerate inputs -- no temporary coordinates and no members at all.
//         Both loops simply never execute, so the answer is vacuously true.
// -----------------------------------------------------------------------------
void test_cm_degenerate_inputs_return_true()
{
    std::cerr << "\n[TEST] test_cm_degenerate_inputs_return_true\n"
              << "  Source file:   src/reactions/conservedMags.cpp\n"
              << "  Function:      conservedMags()\n"
              << "  Scenario A:    member molecule has an empty tmpICoords list.\n"
              << "  Scenario B:    the complex has an empty memberList.\n"
              << "  Pass criteria: conservedMags() returns true in both cases.\n";

    // Scenario A: molecule exists but carries no temporary interface coordinates,
    // so the inner loop body never runs.
    Molecule noTmp = cm_make_molecule(
        Coord { 0.0, 0.0, 0.0 }, { Coord { 1.0, 0.0, 0.0 } },
        Coord { 0.0, 0.0, 0.0 }, {} /* tmpICoords intentionally empty */);
    std::vector<Molecule> listA { noTmp };
    Complex comA = cm_make_complex({ 0 });

    const bool resultA = conservedMags(comA, listA);
    std::cerr << "  Scenario A (empty tmpICoords) -> " << std::boolalpha << resultA << "\n";
    EXPECT_TRUE(resultA) << "With no temporary coordinates there is nothing to violate";

    // Scenario B: complex with no member molecules at all.
    std::vector<Molecule> listB {};
    Complex comB = cm_make_complex({});

    const bool resultB = conservedMags(comB, listB);
    std::cerr << "  Scenario B (empty memberList) -> " << std::boolalpha << resultB << "\n";
    EXPECT_TRUE(resultB) << "An empty complex must trivially satisfy magnitude conservation";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each helper is invoked from its own TEST so a failure
// in one scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(ConservedMags, PureTranslationConserved) { test_cm_pure_translation_conserved(); }
TEST(ConservedMags, PureRotationConserved) { test_cm_pure_rotation_conserved(); }
TEST(ConservedMags, StretchedIfaceNotConserved) { test_cm_stretched_iface_not_conserved(); }
TEST(ConservedMags, TinyDifferenceWithinRounding) { test_cm_tiny_difference_within_rounding(); }
TEST(ConservedMags, SmallDifferenceBeyondRounding) { test_cm_small_difference_beyond_rounding(); }
TEST(ConservedMags, MultiMoleculeSecondBad) { test_cm_multi_molecule_second_bad(); }
TEST(ConservedMags, DegenerateInputsReturnTrue) { test_cm_degenerate_inputs_return_true(); }