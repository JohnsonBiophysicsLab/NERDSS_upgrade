/*! \file test_perform_bimolecular_reactions.cpp
 *
 * ### Unit test for src/reactions/perform_bimolecular_reactions.cpp
 *
 * The function under test is
 *
 *     void perform_bimolecular_reactions(unsigned simItr, Parameters&,
 *              std::vector<Molecule>&, std::vector<Complex>&, SimulVolume&,
 *              std::vector<ForwardRxn>&, std::vector<BackRxn>&,
 *              std::vector<CreateDestructRxn>&, std::vector<MolTemplate>&,
 *              std::map<std::string,int>&, copyCounters&, Membrane&,
 *              std::vector<int>& region)
 *
 * It loops over the molecule indices contained in `region` and, for every
 * molecule, decides between four mutually exclusive outcomes:
 *
 *   1. The molecule is skipped entirely (empty / implicit lipid / ghosted /
 *      already propagated).
 *   2. The molecule has no encounter partners (`crossbase` empty) and has not
 *      moved yet -> its whole complex is given a random propagation vector and
 *      every member molecule is flagged `TrajStatus::canBeResampled`.
 *   3. The molecule has encounter partners but no reaction is drawn -> the
 *      complex diffuses (as in 2) AND the reaction probabilities that its
 *      partners hold *for this molecule* are zeroed out so the partners cannot
 *      react with it later in the same time step.
 *   4. A reaction is drawn -> associate()/state-change is performed. This last
 *      branch involves the full association machinery and is deliberately NOT
 *      exercised here; it is covered by the association unit tests. All the
 *      scenarios below are constructed with zero reaction probability so the
 *      deterministic bookkeeping branches (1-3) can be verified in isolation.
 *
 * NOTE: the function unconditionally opens "DATA/assoc_dissoc_time.dat". If the
 * DATA directory does not exist the stream simply fails to open, which is
 * harmless for the code paths tested here.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"

// The GSL random number generator is a global defined by gtest_main.cpp.
extern gsl_rng* r;

// -----------------------------------------------------------------------------
// Small helpers used to build a minimal but self-consistent system.
// All helper names carry the "pbr_" prefix so they cannot collide with helpers
// from other translation units of the test suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Make sure the GSL RNG exists (diffusion uses GaussV()/rand_gsl()). */
void pbr_ensure_rng() {
  if (r == nullptr) {
    std::cerr << "  [setup] GSL RNG was null -> initializing with srand_gsl(1)\n";
    srand_gsl(1);
  }
}

/*! \brief Minimal simulation parameters: a single molecule type, small dt. */
Parameters pbr_make_params() {
  Parameters params;
  params.timeStep = 0.1;      // microseconds
  params.numMolTypes = 1;
  params.numTotalSpecies = 1;
  params.nItr = 1;
  params.rMaxLimit = 10.0;
  params.rMaxRadius = 10.0;
  params.overlapSepLimit = 0.1;
  params.scaleMaxDisplace = 1e6;  // never veto a displacement
  return params;
}

/*! \brief A large reflecting cubic box so that diffusion never hits a wall. */
Membrane pbr_make_membrane() {
  Membrane membraneObject;
  membraneObject.isBox = true;
  membraneObject.isSphere = false;
  membraneObject.implicitLipid = false;
  membraneObject.TwoD = false;
  membraneObject.waterBox.x = 500.0;
  membraneObject.waterBox.y = 500.0;
  membraneObject.waterBox.z = 500.0;
  membraneObject.waterBox.volume = 500.0 * 500.0 * 500.0;
  membraneObject.waterBox.xLeft = -250.0;
  membraneObject.waterBox.xRight = 250.0;
  membraneObject.xBCtype = "reflect";
  membraneObject.yBCtype = "reflect";
  membraneObject.zBCtype = "reflect";
  // Reflecting-surface lookup table: give it a few benign entries so any
  // indexed access inside the boundary routines is well defined.
  membraneObject.RS3Dvect.assign(16, 0.0);
  membraneObject.numberOfFreeLipidsEachState.assign(1, 0);
  membraneObject.numberOfProteinEachState.assign(1, 0);
  return membraneObject;
}

/*! \brief One molecule template ("A") with a single interface.
 *
 * The static members of MolTemplate (numMolTypes, numEachMolType) are
 * intentionally left untouched: they are shared by the whole test binary.
 */
std::vector<MolTemplate> pbr_make_templates() {
  MolTemplate temp;
  temp.molName = "A";
  temp.molTypeIndex = 0;
  temp.mass = 1.0;
  temp.radius = 1.0;
  temp.copies = 2;
  temp.D = Coord(10.0, 10.0, 10.0);
  temp.Dr = Coord(0.01, 0.01, 0.01);
  temp.isLipid = false;
  temp.isImplicitLipid = false;
  temp.isRod = false;
  temp.isPoint = false;
  temp.checkOverlap = false;

  Interface iface;
  iface.index = 0;
  iface.name = "a";
  iface.iCoord = Coord(1.0, 0.0, 0.0);
  iface.stateList.push_back(Interface::State('\0', 0));
  temp.interfaceList.push_back(iface);

  return std::vector<MolTemplate>{temp};
}

/*! \brief Build a single-interface molecule belonging to complex `comIndex`. */
Molecule pbr_make_molecule(int index, int comIndex, const Coord& com) {
  Molecule mol;
  mol.index = index;
  mol.id = index;
  mol.myComIndex = comIndex;
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
  iface.coord = Coord(com.x + 1.0, com.y, com.z);
  iface.index = 0;
  iface.relIndex = 0;
  iface.molTypeIndex = 0;
  iface.isBound = false;
  mol.interfaceList.push_back(iface);

  mol.freelist.push_back(0);
  return mol;
}

/*! \brief Build a complex owning the molecules listed in `members`. */
Complex pbr_make_complex(int index, const Coord& com,
                         const std::vector<int>& members) {
  Complex targCom;
  targCom.index = index;
  targCom.id = index;
  targCom.comCoord = com;
  targCom.radius = 1.0;
  targCom.mass = static_cast<double>(members.size());
  targCom.D = Coord(10.0, 10.0, 10.0);
  targCom.Dr = Coord(0.01, 0.01, 0.01);
  targCom.memberList = members;
  targCom.numEachMol.assign(1, static_cast<int>(members.size()));
  targCom.lastNumberUpdateItrEachMol.assign(1, 0);
  targCom.isEmpty = false;
  targCom.OnSurface = false;
  targCom.onFiber = false;
  targCom.trajStatus = TrajStatus::none;
  targCom.ncross = 0;
  // trajTrans / trajRot default to (0,0,0), which is what the assertions use
  // as the "did not move" reference value.
  return targCom;
}

/*! \brief One well-formed bimolecular forward reaction A(a) + A(a) -> A(a!1).A(a!1). */
std::vector<ForwardRxn> pbr_make_forward_rxns() {
  ForwardRxn rxn;
  rxn.rxnType = ReactionType::bimolecular;
  rxn.absRxnIndex = 0;
  rxn.relRxnIndex = 0;
  rxn.bindRadius = 1.0;
  rxn.isSymmetric = false;
  rxn.isOnMem = false;
  rxn.isReversible = false;
  rxn.rxnLabel = "pbr_test_rxn";

  RxnIface reactant1("a", 0, 0, 0, '\0', false);
  RxnIface reactant2("a", 0, 0, 0, '\0', false);
  rxn.reactantListNew.push_back(reactant1);
  rxn.reactantListNew.push_back(reactant2);
  rxn.productListNew.push_back(reactant1);
  rxn.productListNew.push_back(reactant2);

  RxnBase::RateState state;
  state.rate = 1.0;
  state.prob = 0.0;
  rxn.rateList.push_back(state);

  return std::vector<ForwardRxn>{rxn};
}

/*! \brief Counter arrays sized generously so no indexed access can go wild. */
copyCounters pbr_make_counters() {
  copyCounters counterArrays;
  counterArrays.copyNumSpecies.assign(8, 0);
  counterArrays.nBoundPairs.assign(4, 0);
  counterArrays.proPairlist.assign(4, 0);
  counterArrays.singleDouble.assign(8, 0);
  counterArrays.implicitDouble.assign(8, false);
  counterArrays.canDissociate.assign(8, false);
  counterArrays.events3D.assign(counterArrays.eventArraySize, 0);
  counterArrays.events2D.assign(counterArrays.eventArraySize, 0);
  counterArrays.events3Dto2D.assign(counterArrays.eventArraySize, 0);
  return counterArrays;
}

/*! \brief Magnitude of a complex's translational trajectory vector. */
double pbr_traj_magnitude(const Complex& targCom) {
  return std::sqrt(targCom.trajTrans.x * targCom.trajTrans.x +
                   targCom.trajTrans.y * targCom.trajTrans.y +
                   targCom.trajTrans.z * targCom.trajTrans.z);
}

/*! \brief Human readable TrajStatus for the console log. */
const char* pbr_traj_status_name(TrajStatus status) {
  switch (status) {
    case TrajStatus::none:
      return "none";
    case TrajStatus::propagated:
      return "propagated";
    case TrajStatus::canBeResampled:
      return "canBeResampled";
  }
  return "unknown";
}

}  // namespace

// -----------------------------------------------------------------------------
// Test 1: an empty region must be a complete no-op.
// -----------------------------------------------------------------------------
void test_pbr_empty_region_is_noop() {
  std::cerr << "\n[TEST] test_pbr_empty_region_is_noop\n"
            << "  Source file:   src/reactions/perform_bimolecular_reactions.cpp\n"
            << "  Function:      perform_bimolecular_reactions\n"
            << "  Scenario:      region contains no molecule indices.\n"
            << "  Pass criteria: molecule TrajStatus stays 'none' and the\n"
            << "                 complex trajectory vector stays exactly zero.\n";

  pbr_ensure_rng();

  Parameters params = pbr_make_params();
  Membrane membraneObject = pbr_make_membrane();
  std::vector<MolTemplate> molTemplateList = pbr_make_templates();
  std::vector<ForwardRxn> forwardRxns = pbr_make_forward_rxns();
  std::vector<BackRxn> backRxns;
  std::vector<CreateDestructRxn> createDestructRxns;
  std::map<std::string, int> observablesList;
  copyCounters counterArrays = pbr_make_counters();
  SimulVolume simulVolume;  // unused by the function, but part of the signature

  std::vector<Molecule> moleculeList{pbr_make_molecule(0, 0, Coord(0, 0, 0))};
  std::vector<Complex> complexList{pbr_make_complex(0, Coord(0, 0, 0), {0})};

  std::vector<int> region;  // deliberately empty

  std::cerr << "  Calling perform_bimolecular_reactions with an empty region...\n";
  perform_bimolecular_reactions(0, params, moleculeList, complexList, simulVolume,
                                forwardRxns, backRxns, createDestructRxns,
                                molTemplateList, observablesList, counterArrays,
                                membraneObject, region);

  // Nothing at all should have been touched.
  EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::none)
      << "A molecule outside the region must not be flagged for resampling";
  EXPECT_DOUBLE_EQ(pbr_traj_magnitude(complexList[0]), 0.0)
      << "A complex outside the region must not receive a propagation vector";

  std::cerr << "  molecule 0 trajStatus = "
            << pbr_traj_status_name(moleculeList[0].trajStatus)
            << ", |trajTrans| = " << pbr_traj_magnitude(complexList[0]) << "\n";
}

// -----------------------------------------------------------------------------
// Test 2: a free molecule (no encounter partners) diffuses, and every member of
//         its complex is flagged canBeResampled - even members not in region.
// -----------------------------------------------------------------------------
void test_pbr_free_molecule_diffuses_whole_complex() {
  std::cerr << "\n[TEST] test_pbr_free_molecule_diffuses_whole_complex\n"
            << "  Source file:   src/reactions/perform_bimolecular_reactions.cpp\n"
            << "  Function:      perform_bimolecular_reactions (crossbase empty branch)\n"
            << "  Scenario:      a two-molecule complex whose first member is in\n"
            << "                 the region; nobody has encounter partners.\n"
            << "  Pass criteria: a non-zero propagation vector is generated and\n"
            << "                 BOTH member molecules become canBeResampled.\n";

  pbr_ensure_rng();

  Parameters params = pbr_make_params();
  Membrane membraneObject = pbr_make_membrane();
  std::vector<MolTemplate> molTemplateList = pbr_make_templates();
  std::vector<ForwardRxn> forwardRxns = pbr_make_forward_rxns();
  std::vector<BackRxn> backRxns;
  std::vector<CreateDestructRxn> createDestructRxns;
  std::map<std::string, int> observablesList;
  copyCounters counterArrays = pbr_make_counters();
  SimulVolume simulVolume;

  // One complex, two member molecules, both sitting near the box centre.
  std::vector<Molecule> moleculeList{pbr_make_molecule(0, 0, Coord(0, 0, 0)),
                                     pbr_make_molecule(1, 0, Coord(3, 0, 0))};
  std::vector<Complex> complexList{pbr_make_complex(0, Coord(1.5, 0, 0), {0, 1})};

  // Only the first member is handed to the routine.
  std::vector<int> region{0};

  std::cerr << "  Calling perform_bimolecular_reactions for molecule 0 only...\n";
  perform_bimolecular_reactions(0, params, moleculeList, complexList, simulVolume,
                                forwardRxns, backRxns, createDestructRxns,
                                molTemplateList, observablesList, counterArrays,
                                membraneObject, region);

  const double mag = pbr_traj_magnitude(complexList[0]);
  std::cerr << "  |trajTrans| after the call = " << mag << "\n"
            << "  molecule 0 trajStatus = "
            << pbr_traj_status_name(moleculeList[0].trajStatus) << "\n"
            << "  molecule 1 trajStatus = "
            << pbr_traj_status_name(moleculeList[1].trajStatus) << "\n";

  // A Gaussian displacement is essentially never identically zero.
  EXPECT_GT(mag, 0.0)
      << "A free complex should have been given a random translation";
  EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::canBeResampled)
      << "The processed molecule must be flagged canBeResampled";
  EXPECT_EQ(moleculeList[1].trajStatus, TrajStatus::canBeResampled)
      << "Every member of the propagated complex must be flagged, not just the "
         "molecule listed in region";
}

// -----------------------------------------------------------------------------
// Test 3: molecules that must be skipped (empty / implicit lipid / ghosted /
//         already propagated) are never propagated.
// -----------------------------------------------------------------------------
void test_pbr_skips_ineligible_molecules() {
  std::cerr << "\n[TEST] test_pbr_skips_ineligible_molecules\n"
            << "  Source file:   src/reactions/perform_bimolecular_reactions.cpp\n"
            << "  Function:      perform_bimolecular_reactions (skip guard)\n"
            << "  Scenario:      four molecules, each in its own complex, that\n"
            << "                 are respectively isEmpty, isImplicitLipid,\n"
            << "                 isGhosted and already TrajStatus::propagated.\n"
            << "  Pass criteria: no complex receives a propagation vector and no\n"
            << "                 TrajStatus is modified.\n";

  pbr_ensure_rng();

  Parameters params = pbr_make_params();
  Membrane membraneObject = pbr_make_membrane();
  std::vector<MolTemplate> molTemplateList = pbr_make_templates();
  std::vector<ForwardRxn> forwardRxns = pbr_make_forward_rxns();
  std::vector<BackRxn> backRxns;
  std::vector<CreateDestructRxn> createDestructRxns;
  std::map<std::string, int> observablesList;
  copyCounters counterArrays = pbr_make_counters();
  SimulVolume simulVolume;

  std::vector<Molecule> moleculeList;
  std::vector<Complex> complexList;
  for (int i = 0; i < 4; ++i) {
    moleculeList.push_back(pbr_make_molecule(i, i, Coord(5.0 * i, 0, 0)));
    complexList.push_back(pbr_make_complex(i, Coord(5.0 * i, 0, 0), {i}));
  }

  // Each molecule trips a different clause of the skip condition.
  moleculeList[0].isEmpty = true;
  moleculeList[1].isImplicitLipid = true;
  moleculeList[2].isGhosted = true;
  moleculeList[3].trajStatus = TrajStatus::propagated;

  std::vector<int> region{0, 1, 2, 3};

  std::cerr << "  Calling perform_bimolecular_reactions on all four molecules...\n";
  perform_bimolecular_reactions(0, params, moleculeList, complexList, simulVolume,
                                forwardRxns, backRxns, createDestructRxns,
                                molTemplateList, observablesList, counterArrays,
                                membraneObject, region);

  for (int i = 0; i < 4; ++i) {
    const double mag = pbr_traj_magnitude(complexList[i]);
    std::cerr << "  complex " << i << " |trajTrans| = " << mag
              << ", molecule trajStatus = "
              << pbr_traj_status_name(moleculeList[i].trajStatus) << "\n";
    EXPECT_DOUBLE_EQ(mag, 0.0)
        << "Skipped molecule " << i << " must not have its complex propagated";
  }

  // The first three keep their original status, the fourth keeps "propagated".
  EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::none)
      << "isEmpty molecule must keep TrajStatus::none";
  EXPECT_EQ(moleculeList[1].trajStatus, TrajStatus::none)
      << "implicit lipid must keep TrajStatus::none";
  EXPECT_EQ(moleculeList[2].trajStatus, TrajStatus::none)
      << "ghosted molecule must keep TrajStatus::none";
  EXPECT_EQ(moleculeList[3].trajStatus, TrajStatus::propagated)
      << "already-propagated molecule must keep TrajStatus::propagated";
}

// -----------------------------------------------------------------------------
// Test 4: a molecule that has already been moved this step (canBeResampled) and
//         has no encounter partners must not be propagated a second time.
// -----------------------------------------------------------------------------
void test_pbr_already_moved_molecule_not_repropagated() {
  std::cerr << "\n[TEST] test_pbr_already_moved_molecule_not_repropagated\n"
            << "  Source file:   src/reactions/perform_bimolecular_reactions.cpp\n"
            << "  Function:      perform_bimolecular_reactions (crossbase empty branch)\n"
            << "  Scenario:      molecule already carries TrajStatus::canBeResampled\n"
            << "                 (its complex moved earlier in the time step).\n"
            << "  Pass criteria: the trajectory vector is left exactly at zero,\n"
            << "                 i.e. create_complex_propagation_vectors is not\n"
            << "                 called a second time.\n";

  pbr_ensure_rng();

  Parameters params = pbr_make_params();
  Membrane membraneObject = pbr_make_membrane();
  std::vector<MolTemplate> molTemplateList = pbr_make_templates();
  std::vector<ForwardRxn> forwardRxns = pbr_make_forward_rxns();
  std::vector<BackRxn> backRxns;
  std::vector<CreateDestructRxn> createDestructRxns;
  std::map<std::string, int> observablesList;
  copyCounters counterArrays = pbr_make_counters();
  SimulVolume simulVolume;

  std::vector<Molecule> moleculeList{pbr_make_molecule(0, 0, Coord(0, 0, 0))};
  std::vector<Complex> complexList{pbr_make_complex(0, Coord(0, 0, 0), {0})};

  // Pretend the complex was already displaced: status is canBeResampled but we
  // reset trajTrans to zero so that any *new* propagation becomes visible.
  moleculeList[0].trajStatus = TrajStatus::canBeResampled;
  complexList[0].trajTrans.x = 0.0;
  complexList[0].trajTrans.y = 0.0;
  complexList[0].trajTrans.z = 0.0;

  std::vector<int> region{0};

  std::cerr << "  Calling perform_bimolecular_reactions on an already moved molecule...\n";
  perform_bimolecular_reactions(0, params, moleculeList, complexList, simulVolume,
                                forwardRxns, backRxns, createDestructRxns,
                                molTemplateList, observablesList, counterArrays,
                                membraneObject, region);

  const double mag = pbr_traj_magnitude(complexList[0]);
  std::cerr << "  |trajTrans| after the call = " << mag
            << " (expected exactly 0)\n";

  EXPECT_DOUBLE_EQ(mag, 0.0)
      << "A molecule whose complex already moved must not be propagated again";
  EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::canBeResampled)
      << "TrajStatus should remain canBeResampled";
}

// -----------------------------------------------------------------------------
// Test 5: encounter partners but zero reaction probability -> the molecule
//         diffuses and the partner's probability entry for it is zeroed.
// -----------------------------------------------------------------------------
void test_pbr_no_reaction_zeros_partner_probability() {
  std::cerr << "\n[TEST] test_pbr_no_reaction_zeros_partner_probability\n"
            << "  Source file:   src/reactions/perform_bimolecular_reactions.cpp\n"
            << "  Function:      perform_bimolecular_reactions (no-reaction branch)\n"
            << "  Scenario:      molecules 0 and 1 list each other as encounter\n"
            << "                 partners. Molecule 0 has probability 0 (so no\n"
            << "                 reaction can be drawn) while molecule 1 still\n"
            << "                 holds a probability of 0.5 for molecule 0.\n"
            << "  Pass criteria: molecule 0's complex diffuses, and molecule 1's\n"
            << "                 probvec entry that points back at molecule 0 is\n"
            << "                 reset to exactly 0.\n";

  pbr_ensure_rng();

  Parameters params = pbr_make_params();
  Membrane membraneObject = pbr_make_membrane();
  std::vector<MolTemplate> molTemplateList = pbr_make_templates();
  std::vector<ForwardRxn> forwardRxns = pbr_make_forward_rxns();
  std::vector<BackRxn> backRxns;
  std::vector<CreateDestructRxn> createDestructRxns;
  std::map<std::string, int> observablesList;
  copyCounters counterArrays = pbr_make_counters();
  SimulVolume simulVolume;

  // Two separate single-molecule complexes.
  std::vector<Molecule> moleculeList{pbr_make_molecule(0, 0, Coord(-2, 0, 0)),
                                     pbr_make_molecule(1, 1, Coord(2, 0, 0))};
  std::vector<Complex> complexList{pbr_make_complex(0, Coord(-2, 0, 0), {0}),
                                   pbr_make_complex(1, Coord(2, 0, 0), {1})};

  // Mutual encounter records: interface 0 of each, forward reaction 0.
  const std::array<int, 3> rxnRecord{0, 0, 0};

  moleculeList[0].crossbase.push_back(1);
  moleculeList[0].mycrossint.push_back(0);
  moleculeList[0].crossrxn.push_back(rxnRecord);
  moleculeList[0].probvec.push_back(0.0);  // cannot react

  moleculeList[1].crossbase.push_back(0);
  moleculeList[1].mycrossint.push_back(0);
  moleculeList[1].crossrxn.push_back(rxnRecord);
  moleculeList[1].probvec.push_back(0.5);  // must be cleared by the routine

  // Only molecule 0 is processed so the effect on molecule 1 is unambiguous.
  std::vector<int> region{0};

  std::cerr << "  Before: molecule 1 probvec[0] = " << moleculeList[1].probvec[0]
            << "\n  Calling perform_bimolecular_reactions for molecule 0...\n";
  perform_bimolecular_reactions(0, params, moleculeList, complexList, simulVolume,
                                forwardRxns, backRxns, createDestructRxns,
                                molTemplateList, observablesList, counterArrays,
                                membraneObject, region);

  const double mag0 = pbr_traj_magnitude(complexList[0]);
  std::cerr << "  After:  molecule 1 probvec[0] = " << moleculeList[1].probvec[0]
            << "\n          complex 0 |trajTrans| = " << mag0
            << "\n          complex 1 |trajTrans| = "
            << pbr_traj_magnitude(complexList[1])
            << "\n          molecule 0 trajStatus = "
            << pbr_traj_status_name(moleculeList[0].trajStatus) << "\n";

  // The non-reacting molecule diffuses ...
  EXPECT_GT(mag0, 0.0)
      << "A molecule that did not react should still diffuse";
  EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::canBeResampled)
      << "The non-reacting molecule must be flagged canBeResampled";

  // ... and its partner can no longer react with it.
  EXPECT_DOUBLE_EQ(moleculeList[1].probvec[0], 0.0)
      << "The partner's reaction probability for this molecule must be zeroed";

  // Molecule 1 was not in the region, so it must not have been propagated.
  EXPECT_DOUBLE_EQ(pbr_traj_magnitude(complexList[1]), 0.0)
      << "Molecule 1 was not in region and must not have been propagated";
  EXPECT_EQ(moleculeList[1].trajStatus, TrajStatus::none)
      << "Molecule 1 was not in region so its TrajStatus must be unchanged";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each scenario runs independently so that a failure in
// one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(PerformBimolecularReactions, EmptyRegionIsNoop) {
  test_pbr_empty_region_is_noop();
}
TEST(PerformBimolecularReactions, FreeMoleculeDiffusesWholeComplex) {
  test_pbr_free_molecule_diffuses_whole_complex();
}
TEST(PerformBimolecularReactions, SkipsIneligibleMolecules) {
  test_pbr_skips_ineligible_molecules();
}
TEST(PerformBimolecularReactions, AlreadyMovedMoleculeNotRepropagated) {
  test_pbr_already_moved_molecule_not_repropagated();
}
TEST(PerformBimolecularReactions, NoReactionZerosPartnerProbability) {
  test_pbr_no_reaction_zeros_partner_probability();
}