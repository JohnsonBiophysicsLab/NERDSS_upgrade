/*! \file test_check_perform_zeroth_first_order_reactions.cpp
 *
 * ### Unit test for src/reactions/check_perform_zeroth_first_order_reactions.cpp
 *
 * The single function defined in that translation unit is
 *
 *     void check_perform_zeroth_first_order_reactions(
 *         unsigned simItr, Parameters&, std::vector<Molecule>&,
 *         std::vector<Complex>&, SimulVolume&, std::vector<ForwardRxn>&,
 *         std::vector<BackRxn>&, std::vector<CreateDestructRxn>&,
 *         std::vector<MolTemplate>&, std::map<std::string,int>&,
 *         copyCounters&, Membrane&, std::vector<double>&,
 *         std::vector<double>&, std::vector<double>&, MpiContext&);
 *
 * It is a *driver* / dispatcher: it opens the association-dissociation log file
 * `DATA/assoc_dissoc_time.dat` and then calls, in order,
 *
 *   1. check_for_unimolecular_reactions()              (always)
 *   2. check_for_zeroth_order_creation()               (always)
 *   3. check_for_unimolstatechange_reactions()         (if params.hasUniMolStateChange)
 *   4. check_dissociation_implicitlipid() per molecule (if params.implicitLipid,
 *      and only for molecules whose complex is OnSurface and which are not
 *      empty / implicit lipids / ghosted)
 *
 * Because every one of those children is a no-op when the reaction lists are
 * empty, the observable, deterministic behaviour that can be unit tested is:
 *
 *   - the routine runs to completion without touching the system when no
 *     zeroth/first order reactions have been defined (molecule and complex
 *     bookkeeping is preserved),
 *   - the log file `DATA/assoc_dissoc_time.dat` is (re)created on every call,
 *   - the optional branches (state change, implicit lipid) are entered without
 *     crashing and still leave an inert system untouched,
 *   - molecules that are flagged empty or ghosted are never revived/altered.
 *
 * Verbose progress information is written to stderr so the reader can follow
 * which source file / function is exercised and what each assertion checks.
 */

#include <gtest/gtest.h>
#include <sys/stat.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Declaration of the function under test lives in this header.  It also drags
// in Parameters, Membrane, SimulVolume, Molecule/Complex, the reaction classes,
// copyCounters and the (forward declared) MpiContext type.
#include "reactions/shared_reaction_functions.hpp"

// For the GSL random number generator handle `r` and srand_gsl().
#include "math/rand_gsl.hpp"

namespace {

// -----------------------------------------------------------------------------
// Small utilities used by every test below.
// -----------------------------------------------------------------------------

/*! \brief Make sure the GSL RNG is usable.
 *
 * `r` is defined (as nullptr) in unit_tests/src/gtest_main.cpp.  Any child
 * routine that draws a random number would dereference it, so we lazily
 * initialise it once with a fixed seed to keep the suite deterministic.
 */
void czfr_init_rng() {
  if (r == nullptr) {
    std::cerr << "  [setup] GSL RNG handle was null -> srand_gsl(1)\n";
    srand_gsl(1);
  }
}

/*! \brief The routine writes DATA/assoc_dissoc_time.dat; make sure DATA exists. */
void czfr_ensure_data_dir() {
  // mkdir failing because the directory already exists is fine.
  mkdir("DATA", 0775);
}

/*! \brief Absolute-ish path of the log file opened by the function under test. */
const char* const kCzfrAssocDissocFile = "DATA/assoc_dissoc_time.dat";

/*! \brief Does a file exist and is it openable? */
bool czfr_file_exists(const char* name) {
  std::ifstream f(name);
  return f.good();
}

// -----------------------------------------------------------------------------
// Minimal, self-consistent simulation system.
// -----------------------------------------------------------------------------

/*! \brief Bundle of every argument the function under test needs. */
struct CzfrSystem {
  Parameters params {};
  Membrane membrane {};
  SimulVolume simulVolume {};
  std::vector<Molecule> moleculeList {};
  std::vector<Complex> complexList {};
  std::vector<MolTemplate> molTemplateList {};
  std::vector<ForwardRxn> forwardRxns {};
  std::vector<BackRxn> backRxns {};
  std::vector<CreateDestructRxn> createDestructRxns {};
  std::map<std::string, int> observablesList {};
  copyCounters counterArrays {};
  std::vector<double> IL2DbindingVec {};
  std::vector<double> IL2DUnbindingVec {};
  std::vector<double> ILTableIDs {};
};

/*! \brief Build the one and only MolTemplate used by these tests.
 *
 * A single point-like molecule type with one stateless, unbound interface.  It
 * takes part in no reaction whatsoever, so every child routine of the function
 * under test must treat it as inert.
 */
MolTemplate czfr_make_template() {
  MolTemplate temp;
  temp.molName = "A";
  temp.molTypeIndex = 0;
  temp.copies = 0;
  temp.mass = 1.0;
  temp.radius = 1.0;
  temp.comCoord = Coord(0.0, 0.0, 0.0);
  temp.D = Coord(1.0, 1.0, 1.0);
  temp.Dr = Coord(0.01, 0.01, 0.01);
  temp.isPoint = true;   // interfaces coincide with the COM
  temp.isLipid = false;
  temp.isRod = false;
  temp.isImplicitLipid = false;
  temp.canDestroy = false;      // no destruction reaction
  temp.checkOverlap = false;
  temp.countTransition = false;

  // One interface named "a" with a single (default) state.  All of the state's
  // reaction lists stay empty, which is what makes the molecule inert.
  Interface iface;
  iface.index = 0;
  iface.name = "a";
  iface.iCoord = Coord(0.0, 0.0, 0.0);
  iface.stateList.push_back(Interface::State('\0', 0));
  temp.interfaceList.push_back(iface);

  // No interface carries an explicit state -> no state change reactions.
  temp.ifacesWithStates.clear();
  return temp;
}

/*! \brief Build one molecule of the template above sitting at \p com. */
Molecule czfr_make_molecule(int index, const Coord& com) {
  Molecule mol;
  mol.index = index;
  mol.id = index;
  mol.myComIndex = index;        // one molecule per complex
  mol.molTypeIndex = 0;
  mol.mySubVolIndex = 0;
  mol.mass = 1.0;
  mol.comCoord = com;
  mol.isEmpty = false;
  mol.isLipid = false;
  mol.isImplicitLipid = false;
  mol.isGhosted = false;
  mol.trajStatus = TrajStatus::none;

  Molecule::Iface iface;
  iface.coord = com;
  iface.index = 0;               // absolute interface/state index
  iface.relIndex = 0;
  iface.molTypeIndex = 0;
  iface.stateIndex = 0;
  iface.stateIden = '\0';
  iface.isBound = false;
  mol.interfaceList.push_back(iface);

  // Interface 0 is free, nothing is bound.
  mol.freelist.push_back(0);
  mol.bndlist.clear();
  mol.bndpartner.clear();
  return mol;
}

/*! \brief Build the parent Complex of molecule \p memberIndex. */
Complex czfr_make_complex(int index, int memberIndex, const Coord& com) {
  Complex com1;
  com1.index = index;
  com1.id = index;
  com1.ownerRank = 0;
  com1.comCoord = com;
  com1.radius = 1.0;
  com1.mass = 1.0;
  com1.D = Coord(1.0, 1.0, 1.0);
  com1.Dr = Coord(0.01, 0.01, 0.01);
  com1.memberList.push_back(memberIndex);
  com1.numEachMol = std::vector<int>(1, 1);
  com1.lastNumberUpdateItrEachMol = std::vector<long long int>(1, 0);
  com1.isEmpty = false;
  com1.OnSurface = false;        // not on the implicit-lipid membrane
  com1.linksToSurface = 0;
  com1.trajStatus = TrajStatus::none;
  return com1;
}

/*! \brief Create a single-cell SimulVolume that owns every molecule.
 *
 * The function under test never rebuilds the member lists (that call is
 * commented out in the source), so a hand-made, trivially valid volume is
 * enough and avoids depending on create_simulation_volume().
 */
void czfr_make_volume(SimulVolume& simulVolume, const Membrane& membrane,
                      std::size_t numMolecules) {
  simulVolume.maxNeighbors = 13;
  simulVolume.numSubCells.x = 1;
  simulVolume.numSubCells.y = 1;
  simulVolume.numSubCells.z = 1;
  simulVolume.numSubCells.tot = 1;
  simulVolume.subCellSize =
      Coord(membrane.waterBox.x, membrane.waterBox.y, membrane.waterBox.z);

  SimulVolume::SubVolume cell;
  cell.absIndex = 0;
  cell.xIndex = 0;
  cell.yIndex = 0;
  cell.zIndex = 0;
  for (std::size_t i = 0; i < numMolecules; ++i)
    cell.memberMolList.push_back(static_cast<int>(i));
  cell.neighborList.clear();

  simulVolume.subCellList.clear();
  simulVolume.subCellList.push_back(cell);
}

/*! \brief Assemble a complete, inert system with \p numMolecules molecules. */
void czfr_build_system(CzfrSystem& sys, int numMolecules) {
  // ---- Parameters -----------------------------------------------------------
  sys.params.rank = 0;
  sys.params.numMolTypes = 1;
  sys.params.numTotalSpecies = 4;   // generous slack for the counter arrays
  sys.params.numTotalComplex = numMolecules;
  sys.params.nItr = 10;
  sys.params.timeStep = 0.1;        // microseconds
  Parameters::dt = sys.params.timeStep;
  Parameters::lastUpdateTransition = std::vector<long long int>(1, 0);
  sys.params.implicitLipid = false;
  sys.params.hasUniMolStateChange = false;
  sys.params.hasCreationDestruction = false;
  sys.params.checkUnimoleculeReactionPopulation = false;
  sys.params.rMaxLimit = 10.0;
  sys.params.rMaxRadius = 5.0;
  sys.params.name = "czfr_unit_test";
  sys.params.debugParams.verbosity = 0;
  sys.params.assocDissocWrite = false;

  // ---- Boundary -------------------------------------------------------------
  sys.membrane.isBox = true;
  sys.membrane.isSphere = false;
  sys.membrane.implicitLipid = false;
  sys.membrane.waterBox = Membrane::WaterBox(std::vector<double>{ 100.0, 100.0, 100.0 });
  sys.membrane.implicitlipidIndex = -1;
  sys.membrane.nStates = 0;
  sys.membrane.No_free_lipids = 0;
  sys.membrane.No_protein = numMolecules;
  sys.membrane.numberOfFreeLipidsEachState.clear();
  sys.membrane.numberOfProteinEachState.clear();
  sys.membrane.xBCtype = "reflect";
  sys.membrane.yBCtype = "reflect";
  sys.membrane.zBCtype = "reflect";

  // ---- Templates / statics --------------------------------------------------
  sys.molTemplateList.clear();
  sys.molTemplateList.push_back(czfr_make_template());
  MolTemplate::numMolTypes = 1;
  MolTemplate::numEachMolType = std::vector<int>(1, numMolecules);
  MolTemplate::absToRelIface = std::vector<int>(1, 0);

  // ---- Molecules and complexes ---------------------------------------------
  sys.moleculeList.clear();
  sys.complexList.clear();
  for (int i = 0; i < numMolecules; ++i) {
    // Spread the molecules out so nothing overlaps.
    Coord com(-20.0 + 10.0 * i, 0.0, 0.0);
    sys.moleculeList.push_back(czfr_make_molecule(i, com));
    sys.complexList.push_back(czfr_make_complex(i, i, com));
  }
  Molecule::numberOfMolecules = numMolecules;
  Molecule::emptyMolList.clear();
  Complex::numberOfComplexes = numMolecules;
  Complex::currNumberMolTypes = 1;
  Complex::currNumberComTypes = 1;
  Complex::emptyComList.clear();
  Complex::obs = std::vector<int>(1, 0);

  // ---- Simulation volume ---------------------------------------------------
  czfr_make_volume(sys.simulVolume, sys.membrane, sys.moleculeList.size());

  // ---- Species counters ----------------------------------------------------
  sys.counterArrays.copyNumSpecies = std::vector<int>(sys.params.numTotalSpecies, 0);
  if (!sys.counterArrays.copyNumSpecies.empty())
    sys.counterArrays.copyNumSpecies[0] = numMolecules;  // free interfaces
  sys.counterArrays.nBoundPairs = std::vector<int>(1, 0);
  sys.counterArrays.proPairlist = std::vector<int>(1, 0);
  sys.counterArrays.singleDouble =
      std::vector<int>(sys.params.numTotalSpecies, 1);
  sys.counterArrays.implicitDouble =
      std::vector<bool>(sys.params.numTotalSpecies, false);
  sys.counterArrays.canDissociate =
      std::vector<bool>(sys.params.numTotalSpecies, false);
  sys.counterArrays.bindPairList.assign(sys.params.numTotalSpecies, {});
  sys.counterArrays.bindPairListIL2D.assign(sys.params.numTotalSpecies, {});
  sys.counterArrays.bindPairListIL3D.assign(sys.params.numTotalSpecies, {});
  sys.counterArrays.events3D.assign(sys.counterArrays.eventArraySize, 0);
  sys.counterArrays.events2D.assign(sys.counterArrays.eventArraySize, 0);
  sys.counterArrays.events3Dto2D.assign(sys.counterArrays.eventArraySize, 0);

  // No reactions at all -> every child routine must be a no-op.
  sys.forwardRxns.clear();
  sys.backRxns.clear();
  sys.createDestructRxns.clear();
  sys.observablesList.clear();
  sys.IL2DbindingVec.clear();
  sys.IL2DUnbindingVec.clear();
  sys.ILTableIDs.clear();
}

/*! \brief Call the function under test on \p sys.
 *
 * NOTE about the MpiContext argument: in the serial build the parameter is
 * completely unused by the implementation -- the only statements that mention
 * it are the commented-out simulVolume.update_memberMolLists() call and the
 * DEBUG_FIND_MOL() macro that sits inside `if (DEBUG)` with `#define DEBUG
 * false`.  MpiContext is merely forward declared in the public headers, so
 * instead of trying to construct an object of a type whose definition is an
 * implementation detail we bind the reference through a null pointer.  It is
 * never dereferenced at run time.
 */
void czfr_invoke(CzfrSystem& sys, unsigned simItr) {
  MpiContext* nullContext = nullptr;
  check_perform_zeroth_first_order_reactions(
      simItr, sys.params, sys.moleculeList, sys.complexList, sys.simulVolume,
      sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
      sys.molTemplateList, sys.observablesList, sys.counterArrays,
      sys.membrane, sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs,
      *nullContext);
}

}  // namespace

// -----------------------------------------------------------------------------
// Test 1: completely empty system.
// -----------------------------------------------------------------------------
void test_czfr_empty_system_is_noop() {
  std::cerr << "\n[TEST] test_czfr_empty_system_is_noop\n"
            << "  Source file:   check_perform_zeroth_first_order_reactions.cpp\n"
            << "  Function:      check_perform_zeroth_first_order_reactions\n"
            << "  Scenario:      no molecules, no complexes, no reactions.\n"
            << "  Pass criteria: the call returns and all containers stay empty.\n";

  czfr_init_rng();
  czfr_ensure_data_dir();

  CzfrSystem sys;
  czfr_build_system(sys, /*numMolecules=*/0);

  std::cerr << "  Calling check_perform_zeroth_first_order_reactions(simItr=0)...\n";
  czfr_invoke(sys, 0);

  // Nothing can be created or destroyed without reactions.
  EXPECT_EQ(sys.moleculeList.size(), 0u)
      << "moleculeList must remain empty when no reactions are defined";
  EXPECT_EQ(sys.complexList.size(), 0u)
      << "complexList must remain empty when no reactions are defined";
  EXPECT_EQ(sys.observablesList.size(), 0u)
      << "no observable may be touched when no reaction fires";

  std::cerr << "  Result: moleculeList.size()=" << sys.moleculeList.size()
            << ", complexList.size()=" << sys.complexList.size() << "\n";
}

// -----------------------------------------------------------------------------
// Test 2: inert molecules must survive the sweep untouched.
// -----------------------------------------------------------------------------
void test_czfr_inert_molecules_unchanged() {
  std::cerr << "\n[TEST] test_czfr_inert_molecules_unchanged\n"
            << "  Source file:   check_perform_zeroth_first_order_reactions.cpp\n"
            << "  Function:      check_perform_zeroth_first_order_reactions\n"
            << "  Scenario:      3 unbound, non-reactive molecules.\n"
            << "  Pass criteria: molecule/complex counts, coordinates, empty\n"
            << "                 flags and species counters are all preserved.\n";

  czfr_init_rng();
  czfr_ensure_data_dir();

  CzfrSystem sys;
  czfr_build_system(sys, /*numMolecules=*/3);

  // Snapshot everything we expect to be invariant.
  const std::size_t numMolBefore = sys.moleculeList.size();
  const std::size_t numComBefore = sys.complexList.size();
  const std::vector<int> countsBefore = sys.counterArrays.copyNumSpecies;
  std::vector<Coord> comBefore;
  for (const auto& mol : sys.moleculeList) comBefore.push_back(mol.comCoord);

  std::cerr << "  Calling check_perform_zeroth_first_order_reactions(simItr=5)...\n";
  czfr_invoke(sys, 5);

  EXPECT_EQ(sys.moleculeList.size(), numMolBefore)
      << "no molecule may be created or removed";
  EXPECT_EQ(sys.complexList.size(), numComBefore)
      << "no complex may be created or removed";

  // Each molecule must still be alive, unbound and at its original position.
  for (std::size_t i = 0; i < sys.moleculeList.size(); ++i) {
    const Molecule& mol = sys.moleculeList[i];
    EXPECT_FALSE(mol.isEmpty) << "molecule " << i << " must not be destroyed";
    EXPECT_EQ(mol.myComIndex, static_cast<int>(i))
        << "molecule " << i << " must keep its parent complex";
    EXPECT_DOUBLE_EQ(mol.comCoord.x, comBefore[i].x)
        << "molecule " << i << " x coordinate must be unchanged";
    EXPECT_DOUBLE_EQ(mol.comCoord.y, comBefore[i].y)
        << "molecule " << i << " y coordinate must be unchanged";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, comBefore[i].z)
        << "molecule " << i << " z coordinate must be unchanged";
    ASSERT_EQ(mol.interfaceList.size(), 1u)
        << "molecule " << i << " should still expose exactly one interface";
    EXPECT_FALSE(mol.interfaceList[0].isBound)
        << "molecule " << i << " interface must remain unbound";
  }

  // Species bookkeeping is only modified by reactions that actually fire.
  ASSERT_EQ(sys.counterArrays.copyNumSpecies.size(), countsBefore.size());
  for (std::size_t i = 0; i < countsBefore.size(); ++i) {
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[i], countsBefore[i])
        << "copyNumSpecies[" << i << "] must be unchanged";
  }

  std::cerr << "  Result: " << sys.moleculeList.size() << " molecules and "
            << sys.complexList.size() << " complexes survived unchanged\n";
}

// -----------------------------------------------------------------------------
// Test 3: the routine (re)creates the association/dissociation log file.
// -----------------------------------------------------------------------------
void test_czfr_creates_assoc_dissoc_log() {
  std::cerr << "\n[TEST] test_czfr_creates_assoc_dissoc_log\n"
            << "  Source file:   check_perform_zeroth_first_order_reactions.cpp\n"
            << "  Function:      check_perform_zeroth_first_order_reactions\n"
            << "  Scenario:      delete DATA/assoc_dissoc_time.dat, then call.\n"
            << "  Pass criteria: the file exists again after the call, because\n"
            << "                 the routine opens it with an std::ofstream.\n";

  czfr_init_rng();
  czfr_ensure_data_dir();

  // Remove any stale file so the check below is meaningful.
  std::remove(kCzfrAssocDissocFile);
  const bool existedBefore = czfr_file_exists(kCzfrAssocDissocFile);
  std::cerr << "  File present before the call? "
            << (existedBefore ? "yes" : "no") << "\n";
  EXPECT_FALSE(existedBefore)
      << "the log file should have been removed before the call";

  CzfrSystem sys;
  czfr_build_system(sys, /*numMolecules=*/1);

  std::cerr << "  Calling check_perform_zeroth_first_order_reactions(simItr=1)...\n";
  czfr_invoke(sys, 1);

  const bool existsAfter = czfr_file_exists(kCzfrAssocDissocFile);
  std::cerr << "  File present after the call?  "
            << (existsAfter ? "yes" : "no") << "\n";
  EXPECT_TRUE(existsAfter)
      << "check_perform_zeroth_first_order_reactions must open (and therefore "
         "create) DATA/assoc_dissoc_time.dat";
}

// -----------------------------------------------------------------------------
// Test 4: the unimolecular state-change branch.
// -----------------------------------------------------------------------------
void test_czfr_unimol_state_change_branch() {
  std::cerr << "\n[TEST] test_czfr_unimol_state_change_branch\n"
            << "  Source file:   check_perform_zeroth_first_order_reactions.cpp\n"
            << "  Function:      check_perform_zeroth_first_order_reactions\n"
            << "  Scenario:      params.hasUniMolStateChange = true but no\n"
            << "                 interface declares a state, so\n"
            << "                 check_for_unimolstatechange_reactions() has\n"
            << "                 nothing to do.\n"
            << "  Pass criteria: branch is entered and the system is unchanged\n"
            << "                 (interface state identity/index preserved).\n";

  czfr_init_rng();
  czfr_ensure_data_dir();

  CzfrSystem sys;
  czfr_build_system(sys, /*numMolecules=*/2);
  sys.params.hasUniMolStateChange = true;   // <- enables the optional branch

  const char stateIdenBefore = sys.moleculeList[0].interfaceList[0].stateIden;
  const int stateIndexBefore = sys.moleculeList[0].interfaceList[0].stateIndex;

  std::cerr << "  Calling check_perform_zeroth_first_order_reactions(simItr=2) "
            << "with hasUniMolStateChange = true...\n";
  czfr_invoke(sys, 2);

  EXPECT_EQ(sys.moleculeList.size(), 2u)
      << "state change reactions never change the molecule count";
  EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, stateIdenBefore)
      << "interface state identity must be preserved when no state change "
         "reaction exists";
  EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, stateIndexBefore)
      << "interface state index must be preserved when no state change "
         "reaction exists";

  std::cerr << "  Result: state identity/index still ('"
            << sys.moleculeList[0].interfaceList[0].stateIden << "', "
            << sys.moleculeList[0].interfaceList[0].stateIndex << ")\n";
}

// -----------------------------------------------------------------------------
// Test 5: the implicit-lipid dissociation branch is guarded by OnSurface.
// -----------------------------------------------------------------------------
void test_czfr_implicit_lipid_branch_skips_off_surface() {
  std::cerr << "\n[TEST] test_czfr_implicit_lipid_branch_skips_off_surface\n"
            << "  Source file:   check_perform_zeroth_first_order_reactions.cpp\n"
            << "  Function:      check_perform_zeroth_first_order_reactions\n"
            << "  Scenario:      params.implicitLipid = true, but every complex\n"
            << "                 has OnSurface == false, so the per-molecule\n"
            << "                 check_dissociation_implicitlipid() call must\n"
            << "                 be skipped by the continue-guard.\n"
            << "  Pass criteria: no crash, no molecule/complex modification and\n"
            << "                 linksToSurface stays zero.\n";

  czfr_init_rng();
  czfr_ensure_data_dir();

  CzfrSystem sys;
  czfr_build_system(sys, /*numMolecules=*/2);
  sys.params.implicitLipid = true;          // <- enables the optional branch
  sys.membrane.implicitLipid = true;
  sys.membrane.implicitlipidIndex = -1;     // no implicit lipid molecule exists
  for (auto& com : sys.complexList) {
    com.OnSurface = false;                  // <- the guard we want to exercise
    com.linksToSurface = 0;
  }

  std::cerr << "  Calling check_perform_zeroth_first_order_reactions(simItr=3) "
            << "with implicitLipid = true...\n";
  czfr_invoke(sys, 3);

  EXPECT_EQ(sys.moleculeList.size(), 2u)
      << "no molecule may be created/removed by the skipped IL branch";
  EXPECT_EQ(sys.complexList.size(), 2u)
      << "no complex may be created/removed by the skipped IL branch";
  for (std::size_t i = 0; i < sys.complexList.size(); ++i) {
    EXPECT_FALSE(sys.complexList[i].OnSurface)
        << "complex " << i << " must still be off the membrane";
    EXPECT_EQ(sys.complexList[i].linksToSurface, 0)
        << "complex " << i << " must not acquire surface links";
    EXPECT_EQ(sys.moleculeList[i].linksToSurface, 0)
        << "molecule " << i << " must not acquire surface links";
  }

  std::cerr << "  Result: implicit-lipid loop entered and correctly skipped "
            << sys.moleculeList.size() << " off-surface molecules\n";
}

// -----------------------------------------------------------------------------
// Test 6: empty (destroyed) and ghosted molecules are ignored.
// -----------------------------------------------------------------------------
void test_czfr_skips_empty_and_ghosted_molecules() {
  std::cerr << "\n[TEST] test_czfr_skips_empty_and_ghosted_molecules\n"
            << "  Source file:   check_perform_zeroth_first_order_reactions.cpp\n"
            << "  Function:      check_perform_zeroth_first_order_reactions\n"
            << "  Scenario:      molecule 1 is flagged isEmpty, molecule 2 is\n"
            << "                 flagged isGhosted; implicit lipids are on so\n"
            << "                 the guarded per-molecule loop runs.\n"
            << "  Pass criteria: flags are preserved, the live molecule is\n"
            << "                 untouched, and nothing is added/removed.\n";

  czfr_init_rng();
  czfr_ensure_data_dir();

  CzfrSystem sys;
  czfr_build_system(sys, /*numMolecules=*/3);
  sys.params.implicitLipid = true;
  sys.membrane.implicitLipid = true;
  sys.membrane.implicitlipidIndex = -1;

  // Molecule 1 has been destroyed, molecule 2 is a copy owned by another rank.
  sys.moleculeList[1].isEmpty = true;
  sys.complexList[1].isEmpty = true;
  sys.moleculeList[2].isGhosted = true;

  const Coord liveComBefore = sys.moleculeList[0].comCoord;

  std::cerr << "  Calling check_perform_zeroth_first_order_reactions(simItr=4)...\n";
  czfr_invoke(sys, 4);

  EXPECT_EQ(sys.moleculeList.size(), 3u)
      << "the molecule list must not be resized";
  EXPECT_TRUE(sys.moleculeList[1].isEmpty)
      << "a destroyed molecule must stay flagged empty";
  EXPECT_TRUE(sys.complexList[1].isEmpty)
      << "the complex of a destroyed molecule must stay flagged empty";
  EXPECT_TRUE(sys.moleculeList[2].isGhosted)
      << "a ghosted molecule must stay flagged ghosted";
  EXPECT_FALSE(sys.moleculeList[0].isEmpty)
      << "the live molecule must not be destroyed";
  EXPECT_DOUBLE_EQ(sys.moleculeList[0].comCoord.x, liveComBefore.x)
      << "the live molecule must not be moved";
  EXPECT_DOUBLE_EQ(sys.moleculeList[0].comCoord.y, liveComBefore.y)
      << "the live molecule must not be moved";
  EXPECT_DOUBLE_EQ(sys.moleculeList[0].comCoord.z, liveComBefore.z)
      << "the live molecule must not be moved";

  std::cerr << "  Result: empty/ghosted flags preserved, live molecule intact\n";
}

// -----------------------------------------------------------------------------
// Test 7: repeated invocations remain stable (idempotent for an inert system).
// -----------------------------------------------------------------------------
void test_czfr_repeated_calls_are_stable() {
  std::cerr << "\n[TEST] test_czfr_repeated_calls_are_stable\n"
            << "  Source file:   check_perform_zeroth_first_order_reactions.cpp\n"
            << "  Function:      check_perform_zeroth_first_order_reactions\n"
            << "  Scenario:      call the routine for several consecutive\n"
            << "                 iterations on the same inert system.\n"
            << "  Pass criteria: the system state is identical after every\n"
            << "                 call (no drift, no leak of molecules).\n";

  czfr_init_rng();
  czfr_ensure_data_dir();

  CzfrSystem sys;
  czfr_build_system(sys, /*numMolecules=*/2);

  const std::size_t numMolBefore = sys.moleculeList.size();
  const std::size_t numComBefore = sys.complexList.size();

  for (unsigned simItr = 0; simItr < 5; ++simItr) {
    std::cerr << "  Iteration " << simItr << ": invoking the routine...\n";
    czfr_invoke(sys, simItr);

    EXPECT_EQ(sys.moleculeList.size(), numMolBefore)
        << "molecule count changed at iteration " << simItr;
    EXPECT_EQ(sys.complexList.size(), numComBefore)
        << "complex count changed at iteration " << simItr;
    for (std::size_t i = 0; i < sys.moleculeList.size(); ++i) {
      EXPECT_FALSE(sys.moleculeList[i].isEmpty)
          << "molecule " << i << " was destroyed at iteration " << simItr;
    }
  }

  std::cerr << "  Result: 5 consecutive sweeps left the system unchanged\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* function is executed inside its own
// TEST so that a failure in one does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(CheckPerformZerothFirstOrderReactions, EmptySystemIsNoop) {
  test_czfr_empty_system_is_noop();
}
TEST(CheckPerformZerothFirstOrderReactions, InertMoleculesUnchanged) {
  test_czfr_inert_molecules_unchanged();
}
TEST(CheckPerformZerothFirstOrderReactions, CreatesAssocDissocLog) {
  test_czfr_creates_assoc_dissoc_log();
}
TEST(CheckPerformZerothFirstOrderReactions, UniMolStateChangeBranch) {
  test_czfr_unimol_state_change_branch();
}
TEST(CheckPerformZerothFirstOrderReactions, ImplicitLipidBranchSkipsOffSurface) {
  test_czfr_implicit_lipid_branch_skips_off_surface();
}
TEST(CheckPerformZerothFirstOrderReactions, SkipsEmptyAndGhostedMolecules) {
  test_czfr_skips_empty_and_ghosted_molecules();
}
TEST(CheckPerformZerothFirstOrderReactions, RepeatedCallsAreStable) {
  test_czfr_repeated_calls_are_stable();
}