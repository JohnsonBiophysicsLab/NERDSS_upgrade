/*! \file test_are_molecules_in_vicinity.cpp
 *
 * ### Unit test for ../src/system_setup/are_molecules_in_vicinity.cpp
 *
 * The single function under test is:
 *
 *     bool are_molecules_in_vicinity(const Molecule& mol1,
 *                                    const Molecule& mol2,
 *                                    const std::vector<MolTemplate>& molTemplates)
 *
 * Behaviour implemented by that source file:
 *   1. Build the vector (mol2.comCoord - mol1.comCoord) and take its magnitude,
 *      i.e. the centre-of-mass to centre-of-mass distance (a double).
 *   2. Compute `radiusSum` as a **float**:
 *          radius(template of mol1) + radius(template of mol2) + 2.0
 *   3. Return `distance < radiusSum`   (STRICTLY less-than).
 *
 * So the pass/fail criteria used throughout this file are:
 *   - distance strictly below the sum of the two template radii plus the
 *     hard-coded safety margin of 2.0  ->  expect true
 *   - distance exactly equal to, or greater than, that sum -> expect false
 *
 * Note that the function indexes `molTemplates` with `mol.molTypeIndex` without
 * any bounds checking, so every Molecule built here is given a valid
 * molTypeIndex and the template list is always fully populated.
 */

#include "system_setup/system_setup.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with amiv_ so they cannot collide with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal MolTemplate carrying only the fields the function reads.
 *
 * `are_molecules_in_vicinity` only ever touches MolTemplate::radius, but we also
 * set molTypeIndex/molName so the object is self-consistent and easy to debug.
 *
 * \param[in] typeIndex index this template will occupy in the template list
 * \param[in] radius    bounding-sphere radius used by the vicinity test
 * \param[in] name      human readable name (for console output only)
 */
MolTemplate amiv_make_template(int typeIndex, double radius, const std::string& name)
{
    MolTemplate temp;
    temp.molTypeIndex = typeIndex;
    temp.radius = radius;
    temp.molName = name;
    return temp;
}

/*! \brief Build a minimal Molecule positioned at a given centre of mass.
 *
 * Only comCoord and molTypeIndex are consulted by the function under test, but
 * index/myComIndex are filled in as well so the object is not half-initialised.
 *
 * \param[in] index      position in a hypothetical moleculeList
 * \param[in] typeIndex  index into the MolTemplate list
 * \param[in] com        centre-of-mass coordinate
 */
Molecule amiv_make_molecule(int index, int typeIndex, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = index;
    mol.molTypeIndex = typeIndex;
    mol.comCoord = com;
    // A single interface coincident with the COM keeps the object valid without
    // affecting the result (interfaces are not read by the function).
    Molecule::Iface iface;
    iface.coord = com;
    iface.relIndex = 0;
    iface.molTypeIndex = typeIndex;
    mol.interfaceList.push_back(iface);
    return mol;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: two molecules sitting nearly on top of each other are "in vicinity".
// -----------------------------------------------------------------------------
void test_amiv_overlapping_molecules_are_in_vicinity()
{
    std::cerr << "\n[TEST] test_amiv_overlapping_molecules_are_in_vicinity\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    both molecules share the same centre of mass.\n"
              << "  Criteria:    distance (0) < radius1+radius2+2 -> expect TRUE\n";

    // Two templates, radius 1.0 each => threshold = 1 + 1 + 2 = 4.0
    std::vector<MolTemplate> molTemplates;
    molTemplates.push_back(amiv_make_template(0, 1.0, "A"));
    molTemplates.push_back(amiv_make_template(1, 1.0, "B"));

    Molecule mol1 = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
    Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 0.0, 0.0, 0.0 });

    bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
    std::cerr << "  distance = 0, threshold = 4 -> returned "
              << std::boolalpha << result << '\n';

    EXPECT_TRUE(result)
        << "Coincident molecules must always be reported as in vicinity";
}

// -----------------------------------------------------------------------------
// Test 2: a clearly separated pair is NOT in vicinity.
// -----------------------------------------------------------------------------
void test_amiv_distant_molecules_are_not_in_vicinity()
{
    std::cerr << "\n[TEST] test_amiv_distant_molecules_are_not_in_vicinity\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    molecules separated by 100 nm along x.\n"
              << "  Criteria:    distance (100) >= threshold (4) -> expect FALSE\n";

    std::vector<MolTemplate> molTemplates;
    molTemplates.push_back(amiv_make_template(0, 1.0, "A"));
    molTemplates.push_back(amiv_make_template(1, 1.0, "B"));

    Molecule mol1 = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
    Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 100.0, 0.0, 0.0 });

    bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
    std::cerr << "  distance = 100, threshold = 4 -> returned "
              << std::boolalpha << result << '\n';

    EXPECT_FALSE(result)
        << "Molecules 100 nm apart with a threshold of 4 nm must not be in vicinity";
}

// -----------------------------------------------------------------------------
// Test 3: exercise the strict "<" comparison right at the threshold.
//         Radii are chosen so the float `radiusSum` is exactly representable
//         (1.0 + 2.0 + 2.0 == 5.0) and no float/double rounding can blur the
//         boundary.
// -----------------------------------------------------------------------------
void test_amiv_threshold_is_strictly_exclusive()
{
    std::cerr << "\n[TEST] test_amiv_threshold_is_strictly_exclusive\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    radii 1.0 and 2.0 -> threshold exactly 5.0.\n"
              << "  Criteria:    distance == 5.0        -> FALSE (comparison is '<')\n"
              << "               distance slightly < 5  -> TRUE\n"
              << "               distance slightly > 5  -> FALSE\n";

    std::vector<MolTemplate> molTemplates;
    molTemplates.push_back(amiv_make_template(0, 1.0, "small"));
    molTemplates.push_back(amiv_make_template(1, 2.0, "large"));

    Molecule mol1 = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });

    // (a) Exactly on the threshold: 1 + 2 + 2 = 5, distance = 5.
    {
        Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 5.0, 0.0, 0.0 });
        bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
        std::cerr << "  distance = 5.0 (exactly threshold) -> returned "
                  << std::boolalpha << result << '\n';
        EXPECT_FALSE(result)
            << "Distance exactly equal to the threshold must return false "
               "because the source uses a strict less-than comparison";
    }

    // (b) Just inside the threshold.
    {
        Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 4.99, 0.0, 0.0 });
        bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
        std::cerr << "  distance = 4.99 (just inside)      -> returned "
                  << std::boolalpha << result << '\n';
        EXPECT_TRUE(result) << "Distance just below the threshold must return true";
    }

    // (c) Just outside the threshold.
    {
        Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 5.01, 0.0, 0.0 });
        bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
        std::cerr << "  distance = 5.01 (just outside)     -> returned "
                  << std::boolalpha << result << '\n';
        EXPECT_FALSE(result) << "Distance just above the threshold must return false";
    }
}

// -----------------------------------------------------------------------------
// Test 4: verify the hard-coded safety margin of 2.0 really is applied, by
//         placing the pair further apart than the bare sum of the radii but
//         still inside radius1 + radius2 + 2.
// -----------------------------------------------------------------------------
void test_amiv_safety_margin_of_two_is_applied()
{
    std::cerr << "\n[TEST] test_amiv_safety_margin_of_two_is_applied\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    radii 1.0 and 1.0 (bare sum 2.0), distance 3.0.\n"
              << "  Criteria:    3.0 > 2.0 but 3.0 < 2.0+2.0 = 4.0 -> expect TRUE,\n"
              << "               proving the +2 margin is included.\n";

    std::vector<MolTemplate> molTemplates;
    molTemplates.push_back(amiv_make_template(0, 1.0, "A"));
    molTemplates.push_back(amiv_make_template(1, 1.0, "B"));

    Molecule mol1 = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
    Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 3.0, 0.0, 0.0 });

    bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
    std::cerr << "  distance = 3.0, bare radius sum = 2.0, threshold = 4.0 -> returned "
              << std::boolalpha << result << '\n';

    EXPECT_TRUE(result)
        << "The function must add a safety margin of 2.0 to the sum of radii";
}

// -----------------------------------------------------------------------------
// Test 5: the distance really is the full 3-D Euclidean magnitude, not just a
//         single component. A (3,4,0) offset has magnitude 5.
// -----------------------------------------------------------------------------
void test_amiv_uses_full_3d_distance()
{
    std::cerr << "\n[TEST] test_amiv_uses_full_3d_distance\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    offset (3,4,0) -> Euclidean distance 5.\n"
              << "  Criteria:    with threshold 6 -> TRUE; with threshold 4 -> FALSE.\n"
              << "               (If only one component were used the second case\n"
              << "                would wrongly report TRUE.)\n";

    // (a) Radii 1.0 + 3.0 -> threshold 6.0, so a distance of 5 is inside.
    {
        std::vector<MolTemplate> molTemplates;
        molTemplates.push_back(amiv_make_template(0, 1.0, "A"));
        molTemplates.push_back(amiv_make_template(1, 3.0, "B"));

        Molecule mol1 = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
        Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 3.0, 4.0, 0.0 });

        bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
        std::cerr << "  |(3,4,0)| = 5, threshold = 6 -> returned "
                  << std::boolalpha << result << '\n';
        EXPECT_TRUE(result) << "3-D distance of 5 should be inside a threshold of 6";
    }

    // (b) Radii 1.0 + 1.0 -> threshold 4.0, so a distance of 5 is outside even
    //     though the individual x (3) and y (4) components are each < 4.
    {
        std::vector<MolTemplate> molTemplates;
        molTemplates.push_back(amiv_make_template(0, 1.0, "A"));
        molTemplates.push_back(amiv_make_template(1, 1.0, "B"));

        Molecule mol1 = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
        Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 3.0, 4.0, 0.0 });

        bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
        std::cerr << "  |(3,4,0)| = 5, threshold = 4 -> returned "
                  << std::boolalpha << result << '\n';
        EXPECT_FALSE(result)
            << "3-D distance of 5 must be outside a threshold of 4, proving the "
               "full vector magnitude (not a single component) is used";
    }

    // (c) A fully 3-D offset (1,2,2) has magnitude 3; confirm with threshold 4.
    {
        std::vector<MolTemplate> molTemplates;
        molTemplates.push_back(amiv_make_template(0, 1.0, "A"));
        molTemplates.push_back(amiv_make_template(1, 1.0, "B"));

        Molecule mol1 = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
        Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 1.0, 2.0, 2.0 });

        const double expectedDist = std::sqrt(1.0 + 4.0 + 4.0); // == 3.0
        bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
        std::cerr << "  |(1,2,2)| = " << expectedDist
                  << ", threshold = 4 -> returned " << std::boolalpha << result << '\n';
        EXPECT_DOUBLE_EQ(expectedDist, 3.0) << "Sanity check on the test geometry";
        EXPECT_TRUE(result) << "A 3-D distance of 3 should be inside a threshold of 4";
    }
}

// -----------------------------------------------------------------------------
// Test 6: the result is symmetric in the two molecule arguments, and negative
//         coordinate offsets are handled correctly (magnitude, not signed sum).
// -----------------------------------------------------------------------------
void test_amiv_is_symmetric_and_sign_independent()
{
    std::cerr << "\n[TEST] test_amiv_is_symmetric_and_sign_independent\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    swap the argument order, and mirror the offset.\n"
              << "  Criteria:    (mol1,mol2) and (mol2,mol1) return the same value;\n"
              << "               +d and -d offsets return the same value.\n";

    // Different radii on purpose so a mis-indexed template would show up.
    std::vector<MolTemplate> molTemplates;
    molTemplates.push_back(amiv_make_template(0, 0.5, "small"));
    molTemplates.push_back(amiv_make_template(1, 4.0, "large"));

    // threshold = 0.5 + 4.0 + 2.0 = 6.5

    // Inside case: distance 3.
    {
        Molecule molA = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
        Molecule molB = amiv_make_molecule(1, 1, Coord{ 3.0, 0.0, 0.0 });

        bool forward = are_molecules_in_vicinity(molA, molB, molTemplates);
        bool reverse = are_molecules_in_vicinity(molB, molA, molTemplates);
        std::cerr << "  inside case (d=3, thr=6.5): forward=" << std::boolalpha
                  << forward << ", reverse=" << reverse << '\n';
        EXPECT_TRUE(forward) << "Distance 3 is inside threshold 6.5";
        EXPECT_EQ(forward, reverse)
            << "Swapping the two molecule arguments must not change the result";
    }

    // Outside case: distance 10.
    {
        Molecule molA = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
        Molecule molB = amiv_make_molecule(1, 1, Coord{ 10.0, 0.0, 0.0 });

        bool forward = are_molecules_in_vicinity(molA, molB, molTemplates);
        bool reverse = are_molecules_in_vicinity(molB, molA, molTemplates);
        std::cerr << "  outside case (d=10, thr=6.5): forward=" << std::boolalpha
                  << forward << ", reverse=" << reverse << '\n';
        EXPECT_FALSE(forward) << "Distance 10 is outside threshold 6.5";
        EXPECT_EQ(forward, reverse)
            << "Swapping the two molecule arguments must not change the result";
    }

    // Mirrored (negative) offset must behave identically to the positive one.
    {
        Molecule molA = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });
        Molecule molPos = amiv_make_molecule(1, 1, Coord{ 3.0, 0.0, 0.0 });
        Molecule molNeg = amiv_make_molecule(1, 1, Coord{ -3.0, 0.0, 0.0 });

        bool pos = are_molecules_in_vicinity(molA, molPos, molTemplates);
        bool neg = are_molecules_in_vicinity(molA, molNeg, molTemplates);
        std::cerr << "  mirrored offsets (+3 vs -3): pos=" << std::boolalpha << pos
                  << ", neg=" << neg << '\n';
        EXPECT_EQ(pos, neg)
            << "The vector magnitude is used, so the sign of the offset is irrelevant";
    }
}

// -----------------------------------------------------------------------------
// Test 7: two molecules of the SAME type index share one template entry, and a
//         list containing several templates is indexed correctly (the function
//         must pick entry [molTypeIndex], not entry [0]).
// -----------------------------------------------------------------------------
void test_amiv_indexes_correct_templates()
{
    std::cerr << "\n[TEST] test_amiv_indexes_correct_templates\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    a 3-entry template list; the pair uses entry 2 twice.\n"
              << "  Criteria:    threshold must come from template[2] (radius 10),\n"
              << "               i.e. 10+10+2 = 22, so a distance of 20 is TRUE while\n"
              << "               the tiny template[0] (radius 0.001) would give FALSE.\n";

    std::vector<MolTemplate> molTemplates;
    molTemplates.push_back(amiv_make_template(0, 0.001, "tiny"));   // decoy
    molTemplates.push_back(amiv_make_template(1, 0.5, "medium"));   // decoy
    molTemplates.push_back(amiv_make_template(2, 10.0, "huge"));    // actually used

    // Both molecules are of type 2 -> threshold = 10 + 10 + 2 = 22.
    Molecule mol1 = amiv_make_molecule(0, 2, Coord{ 0.0, 0.0, 0.0 });
    Molecule mol2 = amiv_make_molecule(1, 2, Coord{ 20.0, 0.0, 0.0 });

    bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
    std::cerr << "  distance = 20, threshold from template[2] = 22 -> returned "
              << std::boolalpha << result << '\n';
    EXPECT_TRUE(result)
        << "Both molecules are type 2 (radius 10), so distance 20 < 22 is in vicinity";

    // Now push them beyond the type-2 threshold to confirm the boundary again.
    Molecule mol3 = amiv_make_molecule(1, 2, Coord{ 25.0, 0.0, 0.0 });
    bool farResult = are_molecules_in_vicinity(mol1, mol3, molTemplates);
    std::cerr << "  distance = 25, threshold = 22 -> returned "
              << std::boolalpha << farResult << '\n';
    EXPECT_FALSE(farResult) << "Distance 25 exceeds the type-2 threshold of 22";
}

// -----------------------------------------------------------------------------
// Test 8: molecules using templates with the default (very small) radius.
//         MolTemplate::radius defaults to 0.0001, so the threshold collapses to
//         essentially the bare safety margin of 2.
// -----------------------------------------------------------------------------
void test_amiv_default_radius_templates()
{
    std::cerr << "\n[TEST] test_amiv_default_radius_templates\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    both templates keep MolTemplate's default radius\n"
              << "               (0.0001), so the threshold is ~2.0002.\n"
              << "  Criteria:    distance 1.0 -> TRUE, distance 2.5 -> FALSE.\n";

    // Deliberately do NOT touch .radius so the struct default (0.0001) is used.
    std::vector<MolTemplate> molTemplates;
    {
        MolTemplate t0;
        t0.molTypeIndex = 0;
        t0.molName = "defaultA";
        molTemplates.push_back(t0);

        MolTemplate t1;
        t1.molTypeIndex = 1;
        t1.molName = "defaultB";
        molTemplates.push_back(t1);
    }

    Molecule mol1 = amiv_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 });

    {
        Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 0.0, 1.0, 0.0 });
        bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
        std::cerr << "  distance = 1.0, threshold ~= 2.0002 -> returned "
                  << std::boolalpha << result << '\n';
        EXPECT_TRUE(result) << "Distance 1.0 is inside the ~2.0002 threshold";
    }

    {
        Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 0.0, 2.5, 0.0 });
        bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
        std::cerr << "  distance = 2.5, threshold ~= 2.0002 -> returned "
                  << std::boolalpha << result << '\n';
        EXPECT_FALSE(result) << "Distance 2.5 is outside the ~2.0002 threshold";
    }
}

// -----------------------------------------------------------------------------
// Test 9: calling the function must not modify either molecule or the template
//         list (all parameters are const references in the signature).
// -----------------------------------------------------------------------------
void test_amiv_does_not_modify_inputs()
{
    std::cerr << "\n[TEST] test_amiv_does_not_modify_inputs\n"
              << "  Source file: are_molecules_in_vicinity.cpp\n"
              << "  Function:    are_molecules_in_vicinity()\n"
              << "  Scenario:    call the function and re-inspect the inputs.\n"
              << "  Criteria:    coordinates, type indices and template radii are\n"
              << "               all unchanged after the call.\n";

    std::vector<MolTemplate> molTemplates;
    molTemplates.push_back(amiv_make_template(0, 1.25, "A"));
    molTemplates.push_back(amiv_make_template(1, 2.75, "B"));

    Molecule mol1 = amiv_make_molecule(0, 0, Coord{ -1.0, 2.0, -3.0 });
    Molecule mol2 = amiv_make_molecule(1, 1, Coord{ 4.0, -5.0, 6.0 });

    // Snapshot the interesting state before the call.
    const Coord com1Before = mol1.comCoord;
    const Coord com2Before = mol2.comCoord;
    const int type1Before = mol1.molTypeIndex;
    const int type2Before = mol2.molTypeIndex;
    const double rad0Before = molTemplates[0].radius;
    const double rad1Before = molTemplates[1].radius;

    bool result = are_molecules_in_vicinity(mol1, mol2, molTemplates);
    std::cerr << "  call returned " << std::boolalpha << result
              << "; now verifying inputs are untouched\n";

    // Coord::operator== compares rounded (4 dp) values, which is exactly what we
    // want here: the values should be bit-identical, so rounded equality holds.
    EXPECT_TRUE(mol1.comCoord == com1Before) << "mol1.comCoord must be unchanged";
    EXPECT_TRUE(mol2.comCoord == com2Before) << "mol2.comCoord must be unchanged";
    EXPECT_EQ(mol1.molTypeIndex, type1Before) << "mol1.molTypeIndex must be unchanged";
    EXPECT_EQ(mol2.molTypeIndex, type2Before) << "mol2.molTypeIndex must be unchanged";
    EXPECT_DOUBLE_EQ(molTemplates[0].radius, rad0Before)
        << "template[0].radius must be unchanged";
    EXPECT_DOUBLE_EQ(molTemplates[1].radius, rad1Before)
        << "template[1].radius must be unchanged";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is invoked from its own TEST so that a
// failure in one scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(AreMoleculesInVicinity, OverlappingMoleculesAreInVicinity)
{
    test_amiv_overlapping_molecules_are_in_vicinity();
}

TEST(AreMoleculesInVicinity, DistantMoleculesAreNotInVicinity)
{
    test_amiv_distant_molecules_are_not_in_vicinity();
}

TEST(AreMoleculesInVicinity, ThresholdIsStrictlyExclusive)
{
    test_amiv_threshold_is_strictly_exclusive();
}

TEST(AreMoleculesInVicinity, SafetyMarginOfTwoIsApplied)
{
    test_amiv_safety_margin_of_two_is_applied();
}

TEST(AreMoleculesInVicinity, UsesFull3DDistance)
{
    test_amiv_uses_full_3d_distance();
}

TEST(AreMoleculesInVicinity, IsSymmetricAndSignIndependent)
{
    test_amiv_is_symmetric_and_sign_independent();
}

TEST(AreMoleculesInVicinity, IndexesCorrectTemplates)
{
    test_amiv_indexes_correct_templates();
}

TEST(AreMoleculesInVicinity, DefaultRadiusTemplates)
{
    test_amiv_default_radius_templates();
}

TEST(AreMoleculesInVicinity, DoesNotModifyInputs)
{
    test_amiv_does_not_modify_inputs();
}