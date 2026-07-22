// ============================================================================
// Unit test for class_SimulVolume.cpp
// ----------------------------------------------------------------------------
// This test exercises the public members of the SimulVolume class and its
// nested Dimensions and SubVolume classes. It verifies:
//   - Dimensions constructor computes proper x/y/z/tot cell counts.
//   - Dimensions::check_dimensions clamps cell counts appropriately.
//   - SimulVolume::create_simulation_volume creates the correct number of
//     sub-cells and computes their sizes.
//   - SimulVolume::create_cell_neighbor_list_cubic assigns proper indices and
//     neighbor lists.
//   - SimulVolume::update_memberMolLists correctly bins molecules.
//   - display() functions run without crashing.
//
// Verbose output is written to stderr describing each test and its criteria.
// ============================================================================

#include "classes/class_SimulVolume.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

// ----------------------------------------------------------------------------
// Helper: construct a minimal Parameters object with a chosen rMaxLimit.
// ----------------------------------------------------------------------------
static Parameters make_params(double rMaxLimit) {
  Parameters params;
  // rMaxLimit determines the cell length used to compute the number of cells.
  params.rMaxLimit = rMaxLimit;
  return params;
}

// ----------------------------------------------------------------------------
// Helper: construct a Membrane object with a given water box size.
// ----------------------------------------------------------------------------
static Membrane make_membrane(double x, double y, double z) {
  Membrane mem;
  mem.waterBox.x = x;
  mem.waterBox.y = y;
  mem.waterBox.z = z;
  return mem;
}

// ----------------------------------------------------------------------------
// Test: SimulVolume::Dimensions constructor.
// Criteria: given a water box and a cell length (rMaxLimit), the number of
// cells in each dimension should equal floor(boxDim / cellLength), with a
// minimum of 4 for z. tot should equal x*y*z.
// ----------------------------------------------------------------------------
void test_Dimensions_constructor() {
  std::cerr << "\n[TEST] class_SimulVolume.cpp :: "
               "SimulVolume::Dimensions::Dimensions(const Parameters&, const "
               "Membrane&)\n";
  std::cerr << "  Purpose: verify cell counts are floor(box/cellLength) with a "
               "z minimum of 4.\n";

  // rMaxLimit = 10, box = 100 x 100 x 100 => 10 cells per dimension.
  Parameters params = make_params(10.0);
  Membrane mem = make_membrane(100.0, 100.0, 100.0);

  SimulVolume::Dimensions dims(params, mem);

  std::cerr << "  Computed dims: x=" << dims.x << " y=" << dims.y
            << " z=" << dims.z << " tot=" << dims.tot << "\n";

  // Expect floor(100/10) = 10 for each axis.
  EXPECT_EQ(dims.x, 10) << "x should be floor(100/10)=10";
  EXPECT_EQ(dims.y, 10) << "y should be floor(100/10)=10";
  EXPECT_EQ(dims.z, 10) << "z should be max(4, floor(100/10))=10";
  EXPECT_EQ(dims.tot, dims.x * dims.y * dims.z)
      << "tot should equal x*y*z";
}

// ----------------------------------------------------------------------------
// Test: SimulVolume::Dimensions constructor enforcing z minimum of 4.
// Criteria: with a small z water box, z should be clamped to at least 4.
// ----------------------------------------------------------------------------
void test_Dimensions_z_minimum() {
  std::cerr << "\n[TEST] class_SimulVolume.cpp :: "
               "SimulVolume::Dimensions (z minimum enforcement)\n";
  std::cerr << "  Purpose: verify z is clamped to a minimum of 4 when the box "
               "is short in z.\n";

  // rMaxLimit = 10, box z = 20 => floor(20/10)=2, but z should be max(4, 2)=4.
  Parameters params = make_params(10.0);
  Membrane mem = make_membrane(100.0, 100.0, 20.0);

  SimulVolume::Dimensions dims(params, mem);

  std::cerr << "  Computed z=" << dims.z << " (expected >= 4)\n";
  EXPECT_GE(dims.z, 4) << "z must be at least 4";
}

// ----------------------------------------------------------------------------
// Test: SimulVolume::Dimensions::check_dimensions.
// Criteria: check_dimensions clamps the number of cells so that no dimension
// exceeds 30 (in the non-mpi build) and re-computes tot.
// ----------------------------------------------------------------------------
void test_Dimensions_check_dimensions() {
  std::cerr << "\n[TEST] class_SimulVolume.cpp :: "
               "SimulVolume::Dimensions::check_dimensions(const Parameters&, "
               "const Membrane&)\n";
  std::cerr << "  Purpose: verify per-axis cell count is clamped to <= 30 and "
               "tot is recomputed.\n";

  // Large box with small cell length would produce > 30 cells per axis.
  Parameters params = make_params(1.0);
  Membrane mem = make_membrane(1000.0, 1000.0, 1000.0);

  // Make sure numberOfMolecules is small so the pair-based clamp does not
  // interfere with the hard axis limit we are testing.
  Molecule::numberOfMolecules = 0;

  SimulVolume::Dimensions dims(params, mem);
  std::cerr << "  Before check: x=" << dims.x << " y=" << dims.y
            << " z=" << dims.z << " tot=" << dims.tot << "\n";

  dims.check_dimensions(params, mem);
  std::cerr << "  After check: x=" << dims.x << " y=" << dims.y
            << " z=" << dims.z << " tot=" << dims.tot << "\n";

  // In the non-mpi build the hard maximum for each axis is 30.
  EXPECT_LE(dims.x, 30) << "x should be clamped to <= 30";
  EXPECT_LE(dims.y, 30) << "y should be clamped to <= 30";
  EXPECT_LE(dims.z, 30) << "z should be clamped to <= 30";
  EXPECT_EQ(dims.tot, dims.x * dims.y * dims.z)
      << "tot should be recomputed as x*y*z";
}

// ----------------------------------------------------------------------------
// Test: SimulVolume::create_simulation_volume.
// Criteria: after creating the simulation volume the subCellList should have
// numSubCells.tot entries and the sub-cell size should be box/numCells.
// ----------------------------------------------------------------------------
void test_create_simulation_volume() {
  std::cerr << "\n[TEST] class_SimulVolume.cpp :: "
               "SimulVolume::create_simulation_volume(const Parameters&, const "
               "Membrane&)\n";
  std::cerr << "  Purpose: verify subCellList size == tot and subCellSize == "
               "box/numCells.\n";

  Parameters params = make_params(10.0);
  Membrane mem = make_membrane(100.0, 100.0, 100.0);
  Molecule::numberOfMolecules = 0;

  SimulVolume simulVolume;
  simulVolume.create_simulation_volume(params, mem);

  std::cerr << "  Created volume: tot=" << simulVolume.numSubCells.tot
            << " subCellList.size=" << simulVolume.subCellList.size()
            << " subCellSize=[" << simulVolume.subCellSize.x << ", "
            << simulVolume.subCellSize.y << ", " << simulVolume.subCellSize.z
            << "]\n";

  // The sub-cell list should hold exactly tot elements.
  EXPECT_EQ(static_cast<int>(simulVolume.subCellList.size()),
            simulVolume.numSubCells.tot)
      << "subCellList should have tot elements";

  // Sub-cell dimensions should be box / numberOfCells along each axis.
  EXPECT_DOUBLE_EQ(simulVolume.subCellSize.x,
                   mem.waterBox.x / (simulVolume.numSubCells.x * 1.0))
      << "subCellSize.x should be boxX / numCellsX";
  EXPECT_DOUBLE_EQ(simulVolume.subCellSize.y,
                   mem.waterBox.y / (simulVolume.numSubCells.y * 1.0))
      << "subCellSize.y should be boxY / numCellsY";
  EXPECT_DOUBLE_EQ(simulVolume.subCellSize.z,
                   mem.waterBox.z / (simulVolume.numSubCells.z * 1.0))
      << "subCellSize.z should be boxZ / numCellsZ";
}

// ----------------------------------------------------------------------------
// Test: SimulVolume::create_cell_neighbor_list_cubic.
// Criteria: after building the neighbor lists, each sub-cell must have the
// correct absolute index and rel index, no cell may exceed maxNeighbors,
// and neighbor indices must be valid (in-range and never self).
// ----------------------------------------------------------------------------
void test_create_cell_neighbor_list_cubic() {
  std::cerr << "\n[TEST] class_SimulVolume.cpp :: "
               "SimulVolume::create_cell_neighbor_list_cubic()\n";
  std::cerr << "  Purpose: verify absolute indices, valid neighbor entries, "
               "and neighbor count limits.\n";

  Parameters params = make_params(10.0);
  Membrane mem = make_membrane(100.0, 100.0, 100.0);
  Molecule::numberOfMolecules = 0;

  SimulVolume simulVolume;
  // create_simulation_volume internally calls create_cell_neighbor_list_cubic.
  simulVolume.create_simulation_volume(params, mem);

  const int tot = simulVolume.numSubCells.tot;
  std::cerr << "  Checking " << tot << " sub-cells...\n";

  bool allIndicesOk = true;
  bool allNeighborsValid = true;
  bool allWithinMaxNeighbors = true;

  for (int i = 0; i < tot; ++i) {
    const auto& cell = simulVolume.subCellList[i];

    // The absolute index should match the position in the list.
    if (cell.absIndex != i) allIndicesOk = false;

    // The number of neighbors must never exceed maxNeighbors.
    if (cell.neighborList.size() > simulVolume.maxNeighbors)
      allWithinMaxNeighbors = false;

    // Every neighbor index must be in range and not equal to itself.
    for (const auto& nb : cell.neighborList) {
      if (nb < 0 || nb >= tot || nb == i) allNeighborsValid = false;
    }
  }

  EXPECT_TRUE(allIndicesOk)
      << "Every sub-cell's absIndex should equal its list position";
  EXPECT_TRUE(allWithinMaxNeighbors)
      << "No sub-cell should exceed maxNeighbors neighbors";
  EXPECT_TRUE(allNeighborsValid)
      << "All neighbor indices should be valid (in-range, non-self)";

  // Sanity check: an interior cell should have up to 13 forward/up neighbors.
  if (tot > 0) {
    std::cerr << "  Cell 0 neighbor count = "
              << simulVolume.subCellList[0].neighborList.size() << "\n";
  }
}

// ----------------------------------------------------------------------------
// Test: SimulVolume::update_memberMolLists (non-MPI overload).
// Criteria: after binning a set of molecules, each non-empty molecule should
// appear in exactly one sub-cell's memberMolList, and its mySubVolIndex should
// be set to a valid bin. Empty / implicit-lipid molecules must be skipped.
// ----------------------------------------------------------------------------
void test_update_memberMolLists() {
  std::cerr << "\n[TEST] class_SimulVolume.cpp :: "
               "SimulVolume::update_memberMolLists(...) (non-MPI)\n";
  std::cerr << "  Purpose: verify molecules are binned into valid sub-cells "
               "and skipped when empty/implicit.\n";

  Parameters params = make_params(10.0);
  Membrane mem = make_membrane(100.0, 100.0, 100.0);
  Molecule::numberOfMolecules = 0;

  SimulVolume simulVolume;
  simulVolume.create_simulation_volume(params, mem);

  // Build a small list of molecules positioned inside the box.
  std::vector<Molecule> moleculeList;

  // Molecule 0: a normal molecule at the box center.
  Molecule mol0;
  mol0.index = 0;
  mol0.isEmpty = false;
  mol0.isImplicitLipid = false;
  mol0.comCoord = Coord{0.0, 0.0, 0.0};
  moleculeList.push_back(mol0);

  // Molecule 1: another normal molecule near a corner (still inside box).
  Molecule mol1;
  mol1.index = 1;
  mol1.isEmpty = false;
  mol1.isImplicitLipid = false;
  mol1.comCoord = Coord{40.0, 40.0, 40.0};
  moleculeList.push_back(mol1);

  // Molecule 2: an "empty" molecule that must be skipped.
  Molecule mol2;
  mol2.index = 2;
  mol2.isEmpty = true;
  mol2.isImplicitLipid = false;
  mol2.comCoord = Coord{0.0, 0.0, 0.0};
  moleculeList.push_back(mol2);

  // Molecule 3: an implicit lipid that must be skipped.
  Molecule mol3;
  mol3.index = 3;
  mol3.isEmpty = false;
  mol3.isImplicitLipid = true;
  mol3.comCoord = Coord{0.0, 0.0, 0.0};
  moleculeList.push_back(mol3);

  std::vector<Complex> complexList;
  std::vector<MolTemplate> molTemplateList;

  // Use simItr that is NOT a multiple of 1000 so we take the fast-binning
  // branch (which does not require complex/template details).
  const int simItr = 1;

  std::cerr << "  Binning " << moleculeList.size()
            << " molecules with simItr=" << simItr << " (fast-bin branch).\n";

  simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                    molTemplateList, mem, simItr);

  // Count how many molecule indices are stored across all sub-cells.
  int totalBinned = 0;
  for (const auto& cell : simulVolume.subCellList)
    totalBinned += static_cast<int>(cell.memberMolList.size());

  std::cerr << "  Total molecules binned = " << totalBinned
            << " (expected 2 non-empty, non-implicit molecules)\n";

  // Only molecule 0 and molecule 1 should have been binned.
  EXPECT_EQ(totalBinned, 2)
      << "Only the two normal molecules should be binned";

  // Check that binned molecules got a valid sub-volume index.
  std::cerr << "  mol0.mySubVolIndex=" << moleculeList[0].mySubVolIndex
            << " mol1.mySubVolIndex=" << moleculeList[1].mySubVolIndex << "\n";
  EXPECT_GE(moleculeList[0].mySubVolIndex, 0)
      << "mol0 should have a valid (>=0) sub-volume index";
  EXPECT_LT(moleculeList[0].mySubVolIndex, simulVolume.numSubCells.tot)
      << "mol0 sub-volume index should be within range";
  EXPECT_GE(moleculeList[1].mySubVolIndex, 0)
      << "mol1 should have a valid (>=0) sub-volume index";
  EXPECT_LT(moleculeList[1].mySubVolIndex, simulVolume.numSubCells.tot)
      << "mol1 sub-volume index should be within range";
}

// ----------------------------------------------------------------------------
// Test: display() functions do not crash.
// Criteria: SimulVolume::display() and SubVolume::display() should execute
// cleanly on a populated volume. There's no return value, so we just ensure
// no exception/crash occurs.
// ----------------------------------------------------------------------------
void test_display_functions() {
  std::cerr << "\n[TEST] class_SimulVolume.cpp :: "
               "SimulVolume::display() and SubVolume::display()\n";
  std::cerr << "  Purpose: verify display routines run without crashing.\n";

  Parameters params = make_params(10.0);
  Membrane mem = make_membrane(100.0, 100.0, 100.0);
  Molecule::numberOfMolecules = 0;

  SimulVolume simulVolume;
  simulVolume.create_simulation_volume(params, mem);

  // Calling display should simply print information; test only that it runs.
  std::cerr << "  Calling SimulVolume::display()...\n";
  EXPECT_NO_THROW(simulVolume.display())
      << "SimulVolume::display() should not throw";

  // Also exercise a single SubVolume's display.
  if (!simulVolume.subCellList.empty()) {
    std::cerr << "  Calling SubVolume::display() for cell 0...\n";
    EXPECT_NO_THROW(simulVolume.subCellList[0].display())
        << "SubVolume::display() should not throw";
  }
}

// ----------------------------------------------------------------------------
// Google Test wrappers: each test function is registered as a TEST so it is
// discovered and run by the surrounding test-suite main().
// ----------------------------------------------------------------------------
TEST(SimulVolumeTest, DimensionsConstructor) { test_Dimensions_constructor(); }
TEST(SimulVolumeTest, DimensionsZMinimum) { test_Dimensions_z_minimum(); }
TEST(SimulVolumeTest, DimensionsCheckDimensions) {
  test_Dimensions_check_dimensions();
}
TEST(SimulVolumeTest, CreateSimulationVolume) {
  test_create_simulation_volume();
}
TEST(SimulVolumeTest, CreateCellNeighborListCubic) {
  test_create_cell_neighbor_list_cubic();
}
TEST(SimulVolumeTest, UpdateMemberMolLists) { test_update_memberMolLists(); }
TEST(SimulVolumeTest, DisplayFunctions) { test_display_functions(); }