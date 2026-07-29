/*! \file test_init_counterCopyNums.cpp
 *
 * ### Unit test for src/io/init_counterCopyNums.cpp
 *
 * Function under test:
 *
 *     void init_counterCopyNums(copyCounters& counterArrays,
 *                               std::vector<Molecule>& moleculeList,
 *                               std::vector<Complex>& complexList,
 *                               std::vector<MolTemplate>& molTemplateList,
 *                               const Membrane& membraneObject,
 *                               int totalSpeciesNum,
 *                               Parameters& params);
 *
 * What the function does (and therefore what we verify):
 *   1. It resizes counterArrays.copyNumSpecies to totalSpeciesNum, all zeroed.
 *   2. It walks over every Molecule (skipping "ghosted" MPI copies) and, for
 *      every interface, increments copyNumSpecies[absolute interface index].
 *   3. For implicit-lipid molecules it *assigns* the copy numbers straight from
 *      Membrane::numberOfFreeLipidsEachState instead of incrementing.
 *   4. When the simulation is NOT restarted and the species can dissociate and
 *      the parent complex has no links to the surface, it records exactly one
 *      of the two partners of a bound pair in counterArrays.bindPairList.
 *   5. Finally, species flagged as "double" (two proteins per species) and not
 *      implicit are halved, since the loop above counted both partners.
 *
 * Every test prints what it is doing and what the pass criteria are, so the
 * console output can be read as a description of the behaviour being checked.
 */

#include <iostream>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Forward declaration of the function under test.
//
// We declare it here rather than including io/io.hpp so this translation unit
// does not have to pull in the whole simulation-volume / parser header chain.
// -----------------------------------------------------------------------------
void init_counterCopyNums(copyCounters& counterArrays,
                          std::vector<Molecule>& moleculeList,
                          std::vector<Complex>& complexList,
                          std::vector<MolTemplate>& molTemplateList,
                          const Membrane& membraneObject, int totalSpeciesNum,
                          Parameters& params);

namespace {

// -----------------------------------------------------------------------------
// Helper: build a copyCounters object whose bookkeeping arrays are all sized
// for `totalSpecies` species and set to benign defaults.
//
// NOTE: copyNumSpecies MUST start empty, because init_counterCopyNums appends
// (push_back) the zeroed entries rather than overwriting them.
// -----------------------------------------------------------------------------
copyCounters iccn_make_counters(int totalSpecies)
{
    copyCounters counters;
    counters.copyNumSpecies.clear();                                  // filled by the function
    counters.singleDouble.assign(totalSpecies, 1);                    // 1 == single molecule species
    counters.implicitDouble.assign(totalSpecies, false);              // not an implicit-lipid product
    counters.canDissociate.assign(totalSpecies, false);               // no dissociation bookkeeping
    counters.bindPairList.assign(totalSpecies, std::vector<int>{});   // explicit bound-pair pool
    counters.bindPairListIL2D.assign(totalSpecies, std::vector<int>{});
    counters.bindPairListIL3D.assign(totalSpecies, std::vector<int>{});
    return counters;
}

// -----------------------------------------------------------------------------
// Helper: build a minimal Molecule belonging to complex `comIndex` and owning
// one interface per entry of `ifaceIndices` (the absolute species indices).
// -----------------------------------------------------------------------------
Molecule iccn_make_molecule(int comIndex, const std::vector<int>& ifaceIndices)
{
    Molecule mol;
    mol.myComIndex = comIndex;
    mol.isGhosted = false;        // owned by this rank (not an MPI ghost copy)
    mol.isImplicitLipid = false;  // a normal explicit molecule
    mol.interfaceList.clear();

    for (std::size_t i = 0; i < ifaceIndices.size(); ++i) {
        Molecule::Iface iface;
        iface.index = ifaceIndices[i];      // absolute species index -> copyNumSpecies slot
        iface.relIndex = static_cast<int>(i);
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

// -----------------------------------------------------------------------------
// Helper: build a Complex owning the given member molecule indices.
// -----------------------------------------------------------------------------
Complex iccn_make_complex(const std::vector<int>& members, int linksToSurface)
{
    Complex com;
    com.memberList = members;
    com.linksToSurface = linksToSurface;
    return com;
}

// -----------------------------------------------------------------------------
// Helper: a Membrane with no implicit-lipid states (the common case).
// -----------------------------------------------------------------------------
Membrane iccn_make_plain_membrane()
{
    Membrane membraneObject;
    membraneObject.nStates = 0;
    membraneObject.numberOfFreeLipidsEachState.clear();
    return membraneObject;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: plain counting of free (single) species.
// -----------------------------------------------------------------------------
void test_iccn_counts_single_species()
{
    std::cerr << "\n[TEST] test_iccn_counts_single_species\n"
              << "  Source file:   src/io/init_counterCopyNums.cpp\n"
              << "  Function:      init_counterCopyNums\n"
              << "  Scenario:      3 explicit molecules, each with interfaces at\n"
              << "                 absolute species indices 0 and 1.\n"
              << "  Pass criteria: copyNumSpecies has size 4, entries 0 and 1 are\n"
              << "                 3 each, entries 2 and 3 remain 0.\n";

    const int totalSpecies = 4;
    copyCounters counters = iccn_make_counters(totalSpecies);

    // Three independent molecules, each in its own complex.
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    for (int i = 0; i < 3; ++i) {
        moleculeList.push_back(iccn_make_molecule(i, { 0, 1 }));
        complexList.push_back(iccn_make_complex({ i }, 0));
    }

    std::vector<MolTemplate> molTemplateList;   // unused by the function
    Membrane membraneObject = iccn_make_plain_membrane();
    Parameters params;
    params.fromRestart = false;

    std::cerr << "  Calling init_counterCopyNums with totalSpeciesNum = "
              << totalSpecies << "...\n";
    init_counterCopyNums(counters, moleculeList, complexList, molTemplateList,
                         membraneObject, totalSpecies, params);

    // The array must be exactly totalSpecies long.
    EXPECT_EQ(counters.copyNumSpecies.size(), static_cast<size_t>(totalSpecies))
        << "copyNumSpecies should be resized to totalSpeciesNum";

    if (counters.copyNumSpecies.size() == static_cast<size_t>(totalSpecies)) {
        EXPECT_EQ(counters.copyNumSpecies[0], 3)
            << "species 0 appears once on each of the 3 molecules";
        EXPECT_EQ(counters.copyNumSpecies[1], 3)
            << "species 1 appears once on each of the 3 molecules";
        EXPECT_EQ(counters.copyNumSpecies[2], 0)
            << "species 2 is not present on any molecule";
        EXPECT_EQ(counters.copyNumSpecies[3], 0)
            << "species 3 is not present on any molecule";

        std::cerr << "  Result copyNumSpecies = ["
                  << counters.copyNumSpecies[0] << ", "
                  << counters.copyNumSpecies[1] << ", "
                  << counters.copyNumSpecies[2] << ", "
                  << counters.copyNumSpecies[3] << "]\n";
    }
}

// -----------------------------------------------------------------------------
// Test 2: ghosted molecules (MPI copies from a neighbouring rank) are skipped.
// -----------------------------------------------------------------------------
void test_iccn_skips_ghosted_molecules()
{
    std::cerr << "\n[TEST] test_iccn_skips_ghosted_molecules\n"
              << "  Source file:   src/io/init_counterCopyNums.cpp\n"
              << "  Function:      init_counterCopyNums\n"
              << "  Scenario:      2 molecules with species index 0, one of which\n"
              << "                 is flagged isGhosted = true.\n"
              << "  Pass criteria: copyNumSpecies[0] == 1 (the ghost is ignored).\n";

    const int totalSpecies = 2;
    copyCounters counters = iccn_make_counters(totalSpecies);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(iccn_make_molecule(0, { 0 }));   // owned molecule
    moleculeList.push_back(iccn_make_molecule(1, { 0 }));   // ghost copy
    moleculeList[1].isGhosted = true;

    std::vector<Complex> complexList;
    complexList.push_back(iccn_make_complex({ 0 }, 0));
    complexList.push_back(iccn_make_complex({ 1 }, 0));

    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject = iccn_make_plain_membrane();
    Parameters params;
    params.fromRestart = false;

    std::cerr << "  Calling init_counterCopyNums...\n";
    init_counterCopyNums(counters, moleculeList, complexList, molTemplateList,
                         membraneObject, totalSpecies, params);

    ASSERT_EQ(counters.copyNumSpecies.size(), static_cast<size_t>(totalSpecies));
    EXPECT_EQ(counters.copyNumSpecies[0], 1)
        << "only the non-ghosted molecule should be counted";
    std::cerr << "  Result copyNumSpecies[0] = " << counters.copyNumSpecies[0]
              << " (expected 1)\n";
}

// -----------------------------------------------------------------------------
// Test 3: bound-pair (double) species are halved, unless flagged implicit.
// -----------------------------------------------------------------------------
void test_iccn_halves_double_species()
{
    std::cerr << "\n[TEST] test_iccn_halves_double_species\n"
              << "  Source file:   src/io/init_counterCopyNums.cpp\n"
              << "  Function:      init_counterCopyNums\n"
              << "  Scenario:      species 0 is a two-protein product counted on\n"
              << "                 both partners; species 1 is also a two-protein\n"
              << "                 product but flagged implicitDouble = true.\n"
              << "  Pass criteria: species 0 count is halved (2 -> 1) while the\n"
              << "                 implicit species 1 keeps its raw count (2).\n";

    const int totalSpecies = 2;
    copyCounters counters = iccn_make_counters(totalSpecies);
    counters.singleDouble[0] = 2;       // bound complex specie -> double counted
    counters.implicitDouble[0] = false; // explicit partner, so halve it
    counters.singleDouble[1] = 2;       // also a "double" specie ...
    counters.implicitDouble[1] = true;  // ... but implicit, so do NOT halve

    // Two molecules bound to each other: each carries species 0 and species 1.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(iccn_make_molecule(0, { 0, 1 }));
    moleculeList.push_back(iccn_make_molecule(0, { 0, 1 }));

    std::vector<Complex> complexList;
    complexList.push_back(iccn_make_complex({ 0, 1 }, 0));

    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject = iccn_make_plain_membrane();
    Parameters params;
    params.fromRestart = false;

    std::cerr << "  Calling init_counterCopyNums...\n";
    init_counterCopyNums(counters, moleculeList, complexList, molTemplateList,
                         membraneObject, totalSpecies, params);

    ASSERT_EQ(counters.copyNumSpecies.size(), static_cast<size_t>(totalSpecies));
    EXPECT_EQ(counters.copyNumSpecies[0], 1)
        << "explicit two-protein species should be halved from 2 to 1";
    EXPECT_EQ(counters.copyNumSpecies[1], 2)
        << "implicitDouble species must not be halved";

    std::cerr << "  Result copyNumSpecies = ["
              << counters.copyNumSpecies[0] << ", "
              << counters.copyNumSpecies[1] << "] (expected [1, 2])\n";
}

// -----------------------------------------------------------------------------
// Test 4: the bound-pair dissociation pool records exactly one partner.
// -----------------------------------------------------------------------------
void test_iccn_populates_bindPairList()
{
    std::cerr << "\n[TEST] test_iccn_populates_bindPairList\n"
              << "  Source file:   src/io/init_counterCopyNums.cpp\n"
              << "  Function:      init_counterCopyNums\n"
              << "  Scenario:      molecules 0 and 1 are bound to each other via a\n"
              << "                 dissociable species (index 0), parent complex has\n"
              << "                 linksToSurface == 0, params.fromRestart == false.\n"
              << "  Pass criteria: bindPairList[0] holds exactly one entry (the first\n"
              << "                 partner, molecule 0); the second partner is skipped\n"
              << "                 because its partner is already in the list.\n";

    const int totalSpecies = 2;
    copyCounters counters = iccn_make_counters(totalSpecies);
    counters.singleDouble[0] = 2;        // bound-pair specie
    counters.canDissociate[0] = true;    // eligible for the dissociation pool

    // Molecule 0 <-> Molecule 1, both in complex 0, both carrying species 0.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(iccn_make_molecule(0, { 0 }));
    moleculeList.push_back(iccn_make_molecule(0, { 0 }));
    moleculeList[0].interfaceList[0].isBound = true;
    moleculeList[0].interfaceList[0].interaction.partnerIndex = 1;
    moleculeList[1].interfaceList[0].isBound = true;
    moleculeList[1].interfaceList[0].interaction.partnerIndex = 0;

    std::vector<Complex> complexList;
    complexList.push_back(iccn_make_complex({ 0, 1 }, 0));  // no surface links

    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject = iccn_make_plain_membrane();
    Parameters params;
    params.fromRestart = false;

    std::cerr << "  Calling init_counterCopyNums...\n";
    init_counterCopyNums(counters, moleculeList, complexList, molTemplateList,
                         membraneObject, totalSpecies, params);

    EXPECT_EQ(counters.bindPairList[0].size(), static_cast<size_t>(1))
        << "only one of the two bound partners should be added to the pool";
    if (counters.bindPairList[0].size() == 1) {
        EXPECT_EQ(counters.bindPairList[0][0], 0)
            << "the first molecule of the pair (index 0) should be the recorded one";
    }

    // The species itself must still be halved (2 raw counts -> 1 bound pair).
    ASSERT_EQ(counters.copyNumSpecies.size(), static_cast<size_t>(totalSpecies));
    EXPECT_EQ(counters.copyNumSpecies[0], 1)
        << "one bound pair should be reported for species 0";

    std::cerr << "  bindPairList[0].size() = " << counters.bindPairList[0].size()
              << ", copyNumSpecies[0] = " << counters.copyNumSpecies[0] << "\n";
}

// -----------------------------------------------------------------------------
// Test 5: on a restart the bindPairList is left alone (it is read from file).
// -----------------------------------------------------------------------------
void test_iccn_restart_skips_bindPairList()
{
    std::cerr << "\n[TEST] test_iccn_restart_skips_bindPairList\n"
              << "  Source file:   src/io/init_counterCopyNums.cpp\n"
              << "  Function:      init_counterCopyNums\n"
              << "  Scenario:      same bound pair as before, but with\n"
              << "                 params.fromRestart = true.\n"
              << "  Pass criteria: bindPairList[0] stays empty while copy numbers\n"
              << "                 are still counted (and halved).\n";

    const int totalSpecies = 2;
    copyCounters counters = iccn_make_counters(totalSpecies);
    counters.singleDouble[0] = 2;
    counters.canDissociate[0] = true;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(iccn_make_molecule(0, { 0 }));
    moleculeList.push_back(iccn_make_molecule(0, { 0 }));
    moleculeList[0].interfaceList[0].interaction.partnerIndex = 1;
    moleculeList[1].interfaceList[0].interaction.partnerIndex = 0;

    std::vector<Complex> complexList;
    complexList.push_back(iccn_make_complex({ 0, 1 }, 0));

    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject = iccn_make_plain_membrane();
    Parameters params;
    params.fromRestart = true;   // <- key difference for this test

    std::cerr << "  Calling init_counterCopyNums (fromRestart = true)...\n";
    init_counterCopyNums(counters, moleculeList, complexList, molTemplateList,
                         membraneObject, totalSpecies, params);

    EXPECT_TRUE(counters.bindPairList[0].empty())
        << "on a restart the dissociation pool must not be rebuilt here";

    ASSERT_EQ(counters.copyNumSpecies.size(), static_cast<size_t>(totalSpecies));
    EXPECT_EQ(counters.copyNumSpecies[0], 1)
        << "copy numbers are still counted and halved on a restart";

    std::cerr << "  bindPairList[0].size() = " << counters.bindPairList[0].size()
              << " (expected 0), copyNumSpecies[0] = "
              << counters.copyNumSpecies[0] << " (expected 1)\n";
}

// -----------------------------------------------------------------------------
// Test 6: implicit-lipid molecules take their copy numbers from the Membrane.
// -----------------------------------------------------------------------------
void test_iccn_implicit_lipid_states()
{
    std::cerr << "\n[TEST] test_iccn_implicit_lipid_states\n"
              << "  Source file:   src/io/init_counterCopyNums.cpp\n"
              << "  Function:      init_counterCopyNums\n"
              << "  Scenario:      molecule 0 is an implicit lipid whose interface\n"
              << "                 starts at absolute index 0, and the Membrane\n"
              << "                 declares 2 states with 100 and 50 free lipids.\n"
              << "  Pass criteria: copyNumSpecies[0] == 100 and copyNumSpecies[1] == 50\n"
              << "                 (assigned, not incremented), and an explicit\n"
              << "                 molecule at species 2 is still counted normally.\n";

    const int totalSpecies = 3;
    copyCounters counters = iccn_make_counters(totalSpecies);

    // Membrane carrying two implicit-lipid states.
    Membrane membraneObject;
    membraneObject.nStates = 2;
    membraneObject.numberOfFreeLipidsEachState = { 100, 50 };

    // Molecule 0: the implicit lipid (must be first in the list, per the code).
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(iccn_make_molecule(0, { 0 }));
    moleculeList[0].isImplicitLipid = true;

    // Molecule 1: a plain explicit molecule carrying species index 2.
    moleculeList.push_back(iccn_make_molecule(1, { 2 }));

    std::vector<Complex> complexList;
    complexList.push_back(iccn_make_complex({ 0 }, 0));
    complexList.push_back(iccn_make_complex({ 1 }, 0));

    std::vector<MolTemplate> molTemplateList;
    Parameters params;
    params.fromRestart = false;

    std::cerr << "  Calling init_counterCopyNums...\n";
    init_counterCopyNums(counters, moleculeList, complexList, molTemplateList,
                         membraneObject, totalSpecies, params);

    ASSERT_EQ(counters.copyNumSpecies.size(), static_cast<size_t>(totalSpecies));
    EXPECT_EQ(counters.copyNumSpecies[0], 100)
        << "implicit lipid state 0 should be taken from the Membrane";
    EXPECT_EQ(counters.copyNumSpecies[1], 50)
        << "implicit lipid state 1 should be taken from the Membrane";
    EXPECT_EQ(counters.copyNumSpecies[2], 1)
        << "the explicit molecule's interface should still be counted";

    std::cerr << "  Result copyNumSpecies = ["
              << counters.copyNumSpecies[0] << ", "
              << counters.copyNumSpecies[1] << ", "
              << counters.copyNumSpecies[2] << "] (expected [100, 50, 1])\n";
}

// -----------------------------------------------------------------------------
// Test 7: an empty molecule list still produces a fully zeroed array.
// -----------------------------------------------------------------------------
void test_iccn_empty_system()
{
    std::cerr << "\n[TEST] test_iccn_empty_system\n"
              << "  Source file:   src/io/init_counterCopyNums.cpp\n"
              << "  Function:      init_counterCopyNums\n"
              << "  Scenario:      no molecules and no complexes at all.\n"
              << "  Pass criteria: copyNumSpecies has size totalSpeciesNum and every\n"
              << "                 entry is zero.\n";

    const int totalSpecies = 5;
    copyCounters counters = iccn_make_counters(totalSpecies);

    std::vector<Molecule> moleculeList;   // empty system
    std::vector<Complex> complexList;
    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject = iccn_make_plain_membrane();
    Parameters params;
    params.fromRestart = false;

    std::cerr << "  Calling init_counterCopyNums on an empty system...\n";
    init_counterCopyNums(counters, moleculeList, complexList, molTemplateList,
                         membraneObject, totalSpecies, params);

    EXPECT_EQ(counters.copyNumSpecies.size(), static_cast<size_t>(totalSpecies))
        << "the array must be allocated even with no molecules present";

    // Every slot must be zero.
    for (int i = 0; i < static_cast<int>(counters.copyNumSpecies.size()); ++i) {
        EXPECT_EQ(counters.copyNumSpecies[i], 0)
            << "species " << i << " should have a copy number of zero";
    }
    std::cerr << "  All " << counters.copyNumSpecies.size()
              << " entries verified to be zero.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper runs inside its own TEST so a
// failure in one does not stop the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(InitCounterCopyNums, CountsSingleSpecies) { test_iccn_counts_single_species(); }
TEST(InitCounterCopyNums, SkipsGhostedMolecules) { test_iccn_skips_ghosted_molecules(); }
TEST(InitCounterCopyNums, HalvesDoubleSpecies) { test_iccn_halves_double_species(); }
TEST(InitCounterCopyNums, PopulatesBindPairList) { test_iccn_populates_bindPairList(); }
TEST(InitCounterCopyNums, RestartSkipsBindPairList) { test_iccn_restart_skips_bindPairList(); }
TEST(InitCounterCopyNums, ImplicitLipidStates) { test_iccn_implicit_lipid_states(); }
TEST(InitCounterCopyNums, EmptySystem) { test_iccn_empty_system(); }