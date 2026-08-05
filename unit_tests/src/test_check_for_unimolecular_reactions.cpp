/*! \file test_check_for_unimolecular_reactions.cpp
 *
 * ### Unit test for src/reactions/check_for_unimolecular_reactions.cpp
 *
 * The function under test is
 *
 * \code
 * void check_for_unimolecular_reactions(unsigned simItr, Parameters&,
 *         std::vector<Molecule>&, std::vector<Complex>&, SimulVolume&,
 *         const std::vector<ForwardRxn>&, const std::vector<BackRxn>&,
 *         const std::vector<CreateDestructRxn>&, std::vector<MolTemplate>&,
 *         std::map<std::string,int>&, copyCounters&, Membrane&,
 *         std::vector<double>&, std::vector<double>&, std::vector<double>&,
 *         std::ofstream&);
 * \endcode
 *
 * It walks the whole moleculeList once per time step and, for every *live*
 * molecule, it
 *   1. skips molecules that are empty (destroyed) or ghosted (MPI),
 *   2. optionally attempts a unimolecular state change (only when
 *      params.hasUniMolStateChange is true and the MolTemplate declares
 *      interfaces with states),
 *   3. attempts unimolecular destruction / creation reactions from the
 *      createDestructRxns list, and
 *   4. finally delegates to check_dissociation() for everything that was not
 *      destroyed, is not an implicit lipid and did not just dissociate.
 *
 * All reactions fire with probability  p = 1 - exp(-rate * dt),  therefore the
 * tests below drive the code with either
 *   - rate = 0            -> p = 0  -> the reaction must NEVER fire, or
 *   - a very large rate   -> p ~ 1  -> the reaction must (essentially) always
 *                                      fire.
 *
 * That makes the outcome of a single call deterministic even though the
 * function consumes random numbers from the GSL generator.
 *
 * Verbose progress information is printed to stderr so that a reader of the
 * test log can follow exactly which behaviour is being probed and what the
 * pass criteria are.
 */

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

// The random number generator used by rand_gsl() is a global owned by
// gtest_main.cpp (declared extern in math/rand_gsl.hpp). It starts out as
// nullptr, so this translation unit lazily allocates it before use.
namespace {

/*! \brief Make sure the global GSL random number generator is usable. */
void cfur_init_rng()
{
    if (r == nullptr) {
          // random generator
        const gsl_rng_type *T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
        std::cerr << "  [setup] allocated global GSL RNG (seed 42)\n";
    }
}

/*! \brief Everything check_for_unimolecular_reactions() needs, in one bag. */
struct CfurSystem {
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

/*! \brief Reset the static bookkeeping members that Molecule/Complex/MolTemplate share.
 *
 * Complex::destroy() decrements these counters, so they have to be consistent
 * with the little system we hand to the function under test.
 */
void cfur_reset_statics(int numMolecules)
{
    Molecule::numberOfMolecules = numMolecules;
    Molecule::emptyMolList.clear();
    // Molecule::mapIdToIndex.clear();

    Complex::numberOfComplexes = numMolecules;
    Complex::emptyComList.clear();
    // Complex::mapIdToIndex.clear();
    Complex::currNumberComTypes = 1;
    Complex::currNumberMolTypes = 1;
    Complex::obs.assign(8, 0);

    MolTemplate::numMolTypes = 1;
    MolTemplate::numEachMolType = std::vector<int>(1, numMolecules);
    MolTemplate::absToRelIface = std::vector<int>(4, 0);

    Parameters::dt = 1.0;
    Parameters::lastUpdateTransition = std::vector<long long int>(1, 0);
}

/*! \brief Build a minimal system: `numMolecules` identical monomers of type "A".
 *
 * Each molecule owns exactly one unbound interface, lives alone in its own
 * complex and is registered in sub-volume 0 of a single-cell SimulVolume.
 */
CfurSystem cfur_build_system(int numMolecules)
{
    cfur_init_rng();
    cfur_reset_statics(numMolecules);

    CfurSystem sys;

    // ---- Parameters -------------------------------------------------------
    sys.params.rank = 0;
    sys.params.numMolTypes = 1;
    sys.params.numTotalSpecies = 4;
    sys.params.numTotalComplex = numMolecules;
    sys.params.nItr = 10;
    sys.params.timeStep = 1.0; // microseconds
    sys.params.hasUniMolStateChange = false;
    sys.params.hasCreationDestruction = false;

    // ---- Membrane / boundary ---------------------------------------------
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    sys.membrane.isBox = true;
    sys.membrane.isSphere = false;
    sys.membrane.implicitLipid = false;
    sys.membrane.nStates = 0;
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";
    sys.membrane.RS3Dvect.assign(100, 0.0);

    // ---- MolTemplate for species "A" -------------------------------------
    MolTemplate temp;
    temp.molName = "A";
    temp.molTypeIndex = 0;
    temp.copies = numMolecules;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord { 10.0, 10.0, 10.0 };
    temp.Dr = Coord { 0.1, 0.1, 0.1 };
    temp.checkOverlap = false;
    temp.countTransition = false;
    temp.canDestroy = true;

    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord { 1.0, 0.0, 0.0 };
    iface.stateList.push_back(Interface::State('\0', 0)); // one stateless state
    temp.interfaceList.push_back(iface);
    sys.molTemplateList.push_back(temp);

    // ---- Molecules and their (monomeric) complexes ------------------------
    for (int i = 0; i < numMolecules; ++i) {
        Molecule mol;
        mol.index = i;
        mol.id = i;
        mol.molTypeIndex = 0;
        mol.myComIndex = i;
        mol.complexId = i;
        mol.mySubVolIndex = 0;
        mol.mass = 1.0;
        mol.comCoord = Coord { 2.0 * i, 0.0, 0.0 };
        mol.trajStatus = TrajStatus::none;
        mol.isEmpty = false;
        mol.isGhosted = false;
        mol.isImplicitLipid = false;
        mol.isDissociated = false;

        Molecule::Iface mIface;
        mIface.coord = Coord { 2.0 * i + 1.0, 0.0, 0.0 };
        mIface.index = 0; // absolute state index
        mIface.relIndex = 0;
        mIface.stateIndex = 0;
        mIface.stateIden = '\0';
        mIface.molTypeIndex = 0;
        mIface.isBound = false;
        mol.interfaceList.push_back(mIface);
        mol.freelist.push_back(0);

        sys.moleculeList.push_back(mol);

        Complex com;
        com.index = i;
        com.id = i;
        com.comCoord = mol.comCoord;
        com.memberList = std::vector<int> { i };
        com.numEachMol = std::vector<int> { 1 };
        com.lastNumberUpdateItrEachMol = std::vector<long long int> { 0 };
        com.D = temp.D;
        com.Dr = temp.Dr;
        com.radius = 1.0;
        com.mass = 1.0;
        com.isEmpty = false;
        com.OnSurface = false;
        com.trajStatus = TrajStatus::none;
        sys.complexList.push_back(com);
    }

    // ---- Single-cell SimulVolume ----------------------------------------
    sys.simulVolume.numSubCells.x = 1;
    sys.simulVolume.numSubCells.y = 1;
    sys.simulVolume.numSubCells.z = 1;
    sys.simulVolume.numSubCells.tot = 1;
    sys.simulVolume.subCellSize = Coord { 100.0, 100.0, 100.0 };
    SimulVolume::SubVolume cell;
    cell.absIndex = 0;
    cell.xIndex = 0;
    cell.yIndex = 0;
    cell.zIndex = 0;
    for (int i = 0; i < numMolecules; ++i)
        cell.memberMolList.push_back(i);
    sys.simulVolume.subCellList.push_back(cell);

    // ---- copyCounters ----------------------------------------------------
    sys.counterArrays.copyNumSpecies.assign(4, numMolecules);
    sys.counterArrays.nBoundPairs.assign(1, 0);
    sys.counterArrays.singleDouble.assign(4, 1);
    sys.counterArrays.implicitDouble.assign(4, false);
    sys.counterArrays.canDissociate.assign(4, false);

    return sys;
}

/*! \brief Create a "A -> 0" destruction reaction with the requested rate.
 *
 * The observable label is always "Adestroyed"; the caller can pre-seed
 * observablesList to check that the counter is decremented on success.
 */
CreateDestructRxn cfur_make_destruction_rxn(double rate)
{
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::destruction;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;

    CreateDestructRxn::CreateDestructMol reactant;
    reactant.molTypeIndex = 0;
    reactant.molName = "A";
    reactant.interfaceList.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    rxn.reactantMolList.push_back(reactant);

    rxn.intReactantList = std::vector<int> { 0 };
    rxn.rateList.push_back(RxnBase::RateState());
    rxn.rateList[0].rate = rate;
    rxn.isObserved = true;
    rxn.observeLabel = "Adestroyed";
    return rxn;
}

/*! \brief Create an "A -> A + A" unimolecular creation reaction. */
CreateDestructRxn cfur_make_creation_rxn(double rate)
{
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::uniMolCreation;
    rxn.absRxnIndex = 1;
    rxn.relRxnIndex = 1;
    rxn.creationRadius = 2.0;

    CreateDestructRxn::CreateDestructMol reactant;
    reactant.molTypeIndex = 0;
    reactant.molName = "A";
    reactant.interfaceList.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    rxn.reactantMolList.push_back(reactant);
    rxn.productMolList.push_back(reactant); // the product is another "A"

    rxn.intReactantList = std::vector<int> { 0 };
    rxn.intProductList = std::vector<int> { 0 };
    rxn.rateList.push_back(RxnBase::RateState());
    rxn.rateList[0].rate = rate;
    rxn.isObserved = false;
    return rxn;
}

/*! \brief Turn the basic system into one that supports A(a~U) <-> A(a~P).
 *
 * Adds two states to the MolTemplate interface, registers the state-change
 * reaction pair on both states, and creates the matching ForwardRxn/BackRxn.
 *
 * \param[in,out] sys       system to modify
 * \param[in]     fwdRate   rate of U -> P
 * \param[in]     backRate  rate of P -> U
 */
void cfur_add_state_change(CfurSystem& sys, double fwdRate, double backRate)
{
    sys.params.hasUniMolStateChange = true;

    // Two states, 'U' (absolute index 0) and 'P' (absolute index 1).
    Interface::State stateU('U', 0);
    stateU.ifaceAndStateName = "a~U";
    stateU.stateChangeRxns.push_back(std::make_pair(0, 0));
    stateU.myForwardRxns.push_back(0u);
    Interface::State stateP('P', 1);
    stateP.ifaceAndStateName = "a~P";
    stateP.stateChangeRxns.push_back(std::make_pair(0, 0));
    stateP.myForwardRxns.push_back(0u);

    sys.molTemplateList[0].interfaceList[0].stateList.clear();
    sys.molTemplateList[0].interfaceList[0].stateList.push_back(stateU);
    sys.molTemplateList[0].interfaceList[0].stateList.push_back(stateP);
    sys.molTemplateList[0].ifacesWithStates = std::vector<int> { 0 };

    // Every molecule starts in state 'U'.
    for (auto& mol : sys.moleculeList) {
        mol.interfaceList[0].stateIden = 'U';
        mol.interfaceList[0].stateIndex = 0;
        mol.interfaceList[0].index = 0;
    }

    // Forward reaction: A(a~U) -> A(a~P)
    ForwardRxn fwd;
    fwd.rxnType = ReactionType::uniMolStateChange;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.hasStateChange = true;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.reactantListNew.push_back(RxnIface("a~U", 0, 0, 0, 'U', false));
    fwd.productListNew.push_back(RxnIface("a~P", 0, 1, 0, 'P', false));
    fwd.stateChangeIface = std::make_pair(fwd.reactantListNew[0], fwd.productListNew[0]);
    fwd.intReactantList = std::vector<int> { 0 };
    fwd.intProductList = std::vector<int> { 1 };
    fwd.rateList.push_back(RxnBase::RateState());
    fwd.rateList[0].rate = fwdRate;
    fwd.rateList[0].otherIfaceLists.resize(1);
    fwd.isObserved = false;
    sys.forwardRxns.push_back(fwd);

    // Back reaction: A(a~P) -> A(a~U)
    BackRxn back;
    back.rxnType = ReactionType::uniMolStateChange;
    back.absRxnIndex = 0;
    back.relRxnIndex = 0;
    back.hasStateChange = true;
    back.conjForwardRxnIndex = 0;
    back.reactantListNew.push_back(RxnIface("a~P", 0, 1, 0, 'P', false));
    back.productListNew.push_back(RxnIface("a~U", 0, 0, 0, 'U', false));
    back.stateChangeIface = std::make_pair(back.reactantListNew[0], back.productListNew[0]);
    back.intReactantList = std::vector<int> { 1 };
    back.intProductList = std::vector<int> { 0 };
    back.rateList.push_back(RxnBase::RateState());
    back.rateList[0].rate = backRate;
    back.rateList[0].otherIfaceLists.resize(1);
    back.isObserved = false;
    sys.backRxns.push_back(back);
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: with no reactions at all the function must be a pure no-op.
// -----------------------------------------------------------------------------
void test_cfur_no_reactions_is_noop()
{
    std::cerr << "\n[TEST] test_cfur_no_reactions_is_noop\n"
              << "  Source file:   src/reactions/check_for_unimolecular_reactions.cpp\n"
              << "  Function:      check_for_unimolecular_reactions()\n"
              << "  Scenario:      one free monomer, empty reaction lists.\n"
              << "  Pass criteria: molecule/complex are untouched, no crash,\n"
              << "                 moleculeList size unchanged.\n";

    CfurSystem sys = cfur_build_system(1);
    std::ofstream assocDissocFile("/dev/null");

    const std::size_t nMolBefore = sys.moleculeList.size();

    std::cerr << "  Calling check_for_unimolecular_reactions (simItr = 1)...\n";
    check_for_unimolecular_reactions(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
        sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);

    EXPECT_EQ(sys.moleculeList.size(), nMolBefore)
        << "No creation reaction exists, so moleculeList must not grow";
    EXPECT_FALSE(sys.moleculeList[0].isEmpty)
        << "No destruction reaction exists, so the molecule must survive";
    EXPECT_FALSE(sys.complexList[0].isEmpty)
        << "The parent complex must survive as well";
    EXPECT_EQ(sys.moleculeList[0].mySubVolIndex, 0)
        << "The molecule must stay registered in its sub-volume";
    EXPECT_EQ(sys.moleculeList[0].trajStatus, TrajStatus::none)
        << "Nothing happened, so trajStatus must remain 'none'";
    EXPECT_EQ(sys.simulVolume.subCellList[0].memberMolList.size(), 1u)
        << "Sub-volume membership list must be unchanged";

    std::cerr << "  Result: molecule isEmpty=" << std::boolalpha
              << sys.moleculeList[0].isEmpty
              << ", mySubVolIndex=" << sys.moleculeList[0].mySubVolIndex << "\n";
}

// -----------------------------------------------------------------------------
// Test 2: destruction with rate 0 => probability 0 => must never fire.
// -----------------------------------------------------------------------------
void test_cfur_destruction_zero_rate_never_fires()
{
    std::cerr << "\n[TEST] test_cfur_destruction_zero_rate_never_fires\n"
              << "  Function:      check_for_unimolecular_reactions() (destruction branch)\n"
              << "  Scenario:      destruction reaction A -> 0 with rate = 0.\n"
              << "  Pass criteria: p = 1-exp(0) = 0, so the molecule must NOT be\n"
              << "                 destroyed and the observable must not change.\n";

    CfurSystem sys = cfur_build_system(1);
    sys.createDestructRxns.push_back(cfur_make_destruction_rxn(0.0));
    sys.observablesList["Adestroyed"] = 7;
    std::ofstream assocDissocFile("/dev/null");

    // Repeat a few times: with p = 0 the outcome must be deterministic.
    for (int rep = 0; rep < 5; ++rep) {
        check_for_unimolecular_reactions(rep, sys.params, sys.moleculeList, sys.complexList,
            sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
            sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
            sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);
    }

    EXPECT_FALSE(sys.moleculeList[0].isEmpty)
        << "A zero-rate destruction reaction must never destroy the molecule";
    EXPECT_FALSE(sys.complexList[0].isEmpty)
        << "A zero-rate destruction reaction must never destroy the complex";
    EXPECT_EQ(sys.moleculeList[0].mySubVolIndex, 0)
        << "The molecule must remain assigned to its sub-volume";
    EXPECT_EQ(sys.observablesList["Adestroyed"], 7)
        << "Observable counter must be untouched when the reaction does not fire";

    std::cerr << "  Result after 5 calls: isEmpty=" << std::boolalpha
              << sys.moleculeList[0].isEmpty
              << ", observable=" << sys.observablesList["Adestroyed"] << "\n";
}

// -----------------------------------------------------------------------------
// Test 3: destruction with an enormous rate => probability ~1 => must fire.
// -----------------------------------------------------------------------------
void test_cfur_destruction_high_rate_destroys_molecule()
{
    std::cerr << "\n[TEST] test_cfur_destruction_high_rate_destroys_molecule\n"
              << "  Function:      check_for_unimolecular_reactions() (destruction branch)\n"
              << "  Scenario:      destruction reaction A -> 0 with rate = 1e9 /s and\n"
              << "                 dt = 1 us, i.e. p = 1 - exp(-1000) ~ 1.\n"
              << "  Pass criteria: the molecule is destroyed (isEmpty), it is removed\n"
              << "                 from the sub-volume member list, mySubVolIndex == -1,\n"
              << "                 and the observable counter is decremented.\n";

    CfurSystem sys = cfur_build_system(1);
    sys.createDestructRxns.push_back(cfur_make_destruction_rxn(1.0e9));
    sys.observablesList["Adestroyed"] = 7;
    std::ofstream assocDissocFile("/dev/null");

    const double prob = 1.0 - std::exp(-1.0e9 * sys.params.timeStep * 1e-6);
    std::cerr << "  Reaction probability for this time step = " << prob << "\n";
    EXPECT_GT(prob, 0.999) << "Sanity check on the test setup: p must be ~1";

    check_for_unimolecular_reactions(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
        sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);

    const bool destroyed = sys.moleculeList[0].isEmpty;
    std::cerr << "  Result: molecule destroyed = " << std::boolalpha << destroyed
              << ", mySubVolIndex = " << sys.moleculeList[0].mySubVolIndex
              << ", subCell members = " << sys.simulVolume.subCellList[0].memberMolList.size()
              << ", observable = " << sys.observablesList["Adestroyed"] << "\n";

    EXPECT_TRUE(destroyed)
        << "With p ~ 1 the lone monomer should be destroyed by the destruction reaction";

    if (destroyed) {
        // Bookkeeping that the function itself is responsible for.
        EXPECT_EQ(sys.moleculeList[0].mySubVolIndex, -1)
            << "A destroyed molecule must have its sub-volume index cleared to -1";
        EXPECT_TRUE(sys.simulVolume.subCellList[0].memberMolList.empty())
            << "A destroyed molecule must be erased from the sub-volume member list";
        EXPECT_EQ(sys.observablesList["Adestroyed"], 6)
            << "The observable of an observed destruction reaction must be decremented";
        EXPECT_TRUE(sys.complexList[0].isEmpty)
            << "Destruction acts on the whole complex, so the complex must be empty";
    }
}

// -----------------------------------------------------------------------------
// Test 4: empty (already destroyed) molecules are skipped.
// -----------------------------------------------------------------------------
void test_cfur_skips_empty_molecule()
{
    std::cerr << "\n[TEST] test_cfur_skips_empty_molecule\n"
              << "  Function:      check_for_unimolecular_reactions() (skip logic)\n"
              << "  Scenario:      two monomers; molecule 0 is flagged isEmpty.\n"
              << "                 A destruction reaction with p ~ 1 is active.\n"
              << "  Pass criteria: molecule 0 is never touched (sub-volume index stays 0,\n"
              << "                 its complex is not destroyed) while molecule 1 IS\n"
              << "                 processed, proving the loop continues past the skip.\n";

    CfurSystem sys = cfur_build_system(2);
    sys.moleculeList[0].isEmpty = true; // pretend it was destroyed earlier
    sys.createDestructRxns.push_back(cfur_make_destruction_rxn(1.0e9));
    sys.observablesList["Adestroyed"] = 2;
    std::ofstream assocDissocFile("/dev/null");

    check_for_unimolecular_reactions(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
        sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);

    std::cerr << "  Result: mol0.mySubVolIndex=" << sys.moleculeList[0].mySubVolIndex
              << ", com0.isEmpty=" << std::boolalpha << sys.complexList[0].isEmpty
              << " | mol1.isEmpty=" << sys.moleculeList[1].isEmpty
              << ", mol1.mySubVolIndex=" << sys.moleculeList[1].mySubVolIndex << "\n";

    // The skipped molecule must not have been processed at all.
    EXPECT_EQ(sys.moleculeList[0].mySubVolIndex, 0)
        << "An 'isEmpty' molecule must be skipped, so its sub-volume index stays 0";
    EXPECT_FALSE(sys.complexList[0].isEmpty)
        << "The complex of a skipped molecule must not be destroyed";

    // The live molecule must have been processed.
    EXPECT_TRUE(sys.moleculeList[1].isEmpty)
        << "The live molecule should still be destroyed (loop must not stop at the skip)";
}

// -----------------------------------------------------------------------------
// Test 5: ghosted molecules (MPI copies) are skipped.
// -----------------------------------------------------------------------------
void test_cfur_skips_ghosted_molecule()
{
    std::cerr << "\n[TEST] test_cfur_skips_ghosted_molecule\n"
              << "  Function:      check_for_unimolecular_reactions() (skip logic)\n"
              << "  Scenario:      single monomer flagged isGhosted, destruction p ~ 1.\n"
              << "  Pass criteria: the ghosted molecule is left completely alone\n"
              << "                 (not empty, still in its sub-volume).\n";

    CfurSystem sys = cfur_build_system(1);
    sys.moleculeList[0].isGhosted = true;
    sys.createDestructRxns.push_back(cfur_make_destruction_rxn(1.0e9));
    std::ofstream assocDissocFile("/dev/null");

    check_for_unimolecular_reactions(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
        sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);

    EXPECT_FALSE(sys.moleculeList[0].isEmpty)
        << "A ghosted molecule must not be destroyed by the owning rank's loop";
    EXPECT_EQ(sys.moleculeList[0].mySubVolIndex, 0)
        << "A ghosted molecule must keep its sub-volume index";
    EXPECT_EQ(sys.simulVolume.subCellList[0].memberMolList.size(), 1u)
        << "A ghosted molecule must remain in the sub-volume member list";

    std::cerr << "  Result: isEmpty=" << std::boolalpha << sys.moleculeList[0].isEmpty
              << ", mySubVolIndex=" << sys.moleculeList[0].mySubVolIndex << "\n";
}

// -----------------------------------------------------------------------------
// Test 6: unimolecular creation with rate 0 must not create anything.
// -----------------------------------------------------------------------------
void test_cfur_creation_zero_rate_creates_nothing()
{
    std::cerr << "\n[TEST] test_cfur_creation_zero_rate_creates_nothing\n"
              << "  Function:      check_for_unimolecular_reactions() (uniMolCreation branch)\n"
              << "  Scenario:      creation reaction A -> A + A with rate = 0.\n"
              << "  Pass criteria: moleculeList/complexList sizes are unchanged and the\n"
              << "                 parent molecule keeps trajStatus 'none'.\n";

    CfurSystem sys = cfur_build_system(1);
    sys.createDestructRxns.push_back(cfur_make_creation_rxn(0.0));
    std::ofstream assocDissocFile("/dev/null");

    const std::size_t nMol = sys.moleculeList.size();
    const std::size_t nCom = sys.complexList.size();

    check_for_unimolecular_reactions(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
        sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);

    EXPECT_EQ(sys.moleculeList.size(), nMol)
        << "A zero-rate creation reaction must not add molecules";
    EXPECT_EQ(sys.complexList.size(), nCom)
        << "A zero-rate creation reaction must not add complexes";
    EXPECT_EQ(sys.moleculeList[0].trajStatus, TrajStatus::none)
        << "The parent molecule must not be marked as propagated";

    std::cerr << "  Result: moleculeList.size()=" << sys.moleculeList.size()
              << ", complexList.size()=" << sys.complexList.size() << "\n";
}

// -----------------------------------------------------------------------------
// Test 7: unimolecular state change with rate 0 must leave the state alone.
// -----------------------------------------------------------------------------
void test_cfur_state_change_zero_rate_keeps_state()
{
    std::cerr << "\n[TEST] test_cfur_state_change_zero_rate_keeps_state\n"
              << "  Function:      check_for_unimolecular_reactions() (state change branch)\n"
              << "  Scenario:      A(a~U) <-> A(a~P) with both rates = 0 and\n"
              << "                 params.hasUniMolStateChange = true.\n"
              << "  Pass criteria: interface stays in state 'U' with absolute index 0 and\n"
              << "                 trajStatus stays 'none'.\n";

    CfurSystem sys = cfur_build_system(1);
    cfur_add_state_change(sys, 0.0, 0.0);
    std::ofstream assocDissocFile("/dev/null");

    for (int rep = 0; rep < 5; ++rep) {
        check_for_unimolecular_reactions(rep, sys.params, sys.moleculeList, sys.complexList,
            sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
            sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
            sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);
    }

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIden, 'U')
        << "With rate 0 the interface must stay in its initial state 'U'";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 0)
        << "With rate 0 the absolute state index must stay 0";
    EXPECT_EQ(sys.moleculeList[0].trajStatus, TrajStatus::none)
        << "No state change occurred, so trajStatus must remain 'none'";

    std::cerr << "  Result after 5 calls: stateIden='"
              << sys.moleculeList[0].interfaceList[0].stateIden
              << "', absIndex=" << sys.moleculeList[0].interfaceList[0].index << "\n";
}

// -----------------------------------------------------------------------------
// Test 8: unimolecular state change with a very large rate.
//
// Whether the reaction is actually located depends on the reaction bookkeeping
// (find_which_state_change_reaction), so this test verifies *self consistency*:
// either nothing changed, or the interface moved to the product state and the
// molecule/complex were both flagged as propagated.
// -----------------------------------------------------------------------------
void test_cfur_state_change_high_rate_is_consistent()
{
    std::cerr << "\n[TEST] test_cfur_state_change_high_rate_is_consistent\n"
              << "  Function:      check_for_unimolecular_reactions() (state change branch)\n"
              << "  Scenario:      A(a~U) -> A(a~P) with forward rate 1e9 /s (p ~ 1)\n"
              << "                 and back rate 0.\n"
              << "  Pass criteria: the resulting state is one of the two declared states,\n"
              << "                 and IF the state changed to 'P' then the absolute state\n"
              << "                 index is 1 and molecule + complex are 'propagated'.\n";

    CfurSystem sys = cfur_build_system(1);
    cfur_add_state_change(sys, 1.0e9, 0.0);
    std::ofstream assocDissocFile("/dev/null");

    check_for_unimolecular_reactions(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
        sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);

    const char finalState = sys.moleculeList[0].interfaceList[0].stateIden;
    std::cerr << "  Result: stateIden='" << finalState
              << "', absIndex=" << sys.moleculeList[0].interfaceList[0].index
              << ", molecule trajStatus="
              << static_cast<int>(sys.moleculeList[0].trajStatus) << "\n";

    // The state must remain one of the declared states in any case.
    EXPECT_TRUE(finalState == 'U' || finalState == 'P')
        << "The interface must end in one of the two declared states ('U' or 'P')";

    if (finalState == 'P') {
        std::cerr << "    -> the state change fired; verifying the bookkeeping\n";
        EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 1)
            << "After U -> P the absolute state index must be the product index (1)";
        EXPECT_EQ(sys.moleculeList[0].interfaceList[0].stateIndex, 1)
            << "After U -> P the relative state index must be 1";
        EXPECT_EQ(sys.moleculeList[0].trajStatus, TrajStatus::propagated)
            << "A molecule that changed state must be marked as propagated";
        EXPECT_EQ(sys.complexList[0].trajStatus, TrajStatus::propagated)
            << "The parent complex of a changed molecule must be marked as propagated";
    } else {
        std::cerr << "    -> no state change fired; verifying nothing moved\n";
        EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 0)
            << "If no state change occurred the absolute state index must stay 0";
    }

    // In neither case may the molecule be destroyed - no destruction reaction exists.
    EXPECT_FALSE(sys.moleculeList[0].isEmpty)
        << "No destruction reaction exists, so the molecule must survive";
}

// -----------------------------------------------------------------------------
// Test 9: implicit lipids skip the dissociation call and their state counters
//         are left alone when no reactions are defined.
// -----------------------------------------------------------------------------
void test_cfur_implicit_lipid_is_skipped_for_dissociation()
{
    std::cerr << "\n[TEST] test_cfur_implicit_lipid_is_skipped_for_dissociation\n"
              << "  Function:      check_for_unimolecular_reactions() (implicit lipid path)\n"
              << "  Scenario:      the single molecule is an implicit lipid with one state\n"
              << "                 and 10 free lipids; no reactions are defined.\n"
              << "  Pass criteria: the free-lipid counter is unchanged, the molecule is not\n"
              << "                 destroyed and the call returns without crashing.\n";

    CfurSystem sys = cfur_build_system(1);
    sys.moleculeList[0].isImplicitLipid = true;
    sys.molTemplateList[0].isImplicitLipid = true;
    sys.molTemplateList[0].isLipid = true;
    sys.membrane.implicitLipid = true;
    sys.membrane.implicitlipidIndex = 0;
    sys.membrane.nStates = 1;
    sys.membrane.numberOfFreeLipidsEachState = std::vector<int> { 10 };
    sys.membrane.numberOfProteinEachState = std::vector<int> { 0 };
    std::ofstream assocDissocFile("/dev/null");

    check_for_unimolecular_reactions(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
        sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);

    ASSERT_EQ(sys.membrane.numberOfFreeLipidsEachState.size(), 1u)
        << "The free-lipid state vector must keep its size";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], 10)
        << "Without reactions the number of free implicit lipids must stay 10";
    EXPECT_FALSE(sys.moleculeList[0].isEmpty)
        << "The implicit lipid molecule must not be destroyed";
    EXPECT_EQ(sys.moleculeList[0].mySubVolIndex, 0)
        << "The implicit lipid must keep its sub-volume index";

    std::cerr << "  Result: free lipids = " << sys.membrane.numberOfFreeLipidsEachState[0]
              << ", isEmpty = " << std::boolalpha << sys.moleculeList[0].isEmpty << "\n";
}

// -----------------------------------------------------------------------------
// Test 10: a molecule that already dissociated this step is skipped before the
//          dissociation check (and must otherwise be left intact).
// -----------------------------------------------------------------------------
void test_cfur_just_dissociated_molecule_is_skipped()
{
    std::cerr << "\n[TEST] test_cfur_just_dissociated_molecule_is_skipped\n"
              << "  Function:      check_for_unimolecular_reactions() (isDissociated skip)\n"
              << "  Scenario:      the molecule has isDissociated = true; no reactions.\n"
              << "  Pass criteria: no crash, flags/indices untouched, molecule alive.\n";

    CfurSystem sys = cfur_build_system(1);
    sys.moleculeList[0].isDissociated = true;
    sys.moleculeList[0].trajStatus = TrajStatus::propagated;
    std::ofstream assocDissocFile("/dev/null");

    check_for_unimolecular_reactions(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.membrane,
        sys.IL2DbindingVec, sys.IL2DUnbindingVec, sys.ILTableIDs, assocDissocFile);

    EXPECT_TRUE(sys.moleculeList[0].isDissociated)
        << "The isDissociated flag must not be cleared by this function";
    EXPECT_FALSE(sys.moleculeList[0].isEmpty)
        << "The molecule must remain alive (no destruction reaction defined)";
    EXPECT_EQ(sys.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "trajStatus must be preserved for a molecule that already moved";
    EXPECT_EQ(sys.moleculeList[0].mySubVolIndex, 0)
        << "The sub-volume index must be preserved";

    std::cerr << "  Result: isDissociated=" << std::boolalpha
              << sys.moleculeList[0].isDissociated
              << ", isEmpty=" << sys.moleculeList[0].isEmpty << "\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is executed inside its own TEST so
// that every scenario is reported separately and a failure in one does not
// prevent the others from running.
// -----------------------------------------------------------------------------
TEST(CheckForUnimolecularReactions, NoReactionsIsNoop) { test_cfur_no_reactions_is_noop(); }
TEST(CheckForUnimolecularReactions, DestructionZeroRateNeverFires) { test_cfur_destruction_zero_rate_never_fires(); }
TEST(CheckForUnimolecularReactions, DestructionHighRateDestroysMolecule) { test_cfur_destruction_high_rate_destroys_molecule(); }
TEST(CheckForUnimolecularReactions, SkipsEmptyMolecule) { test_cfur_skips_empty_molecule(); }
TEST(CheckForUnimolecularReactions, SkipsGhostedMolecule) { test_cfur_skips_ghosted_molecule(); }
TEST(CheckForUnimolecularReactions, CreationZeroRateCreatesNothing) { test_cfur_creation_zero_rate_creates_nothing(); }
TEST(CheckForUnimolecularReactions, StateChangeZeroRateKeepsState) { test_cfur_state_change_zero_rate_keeps_state(); }
TEST(CheckForUnimolecularReactions, StateChangeHighRateIsConsistent) { test_cfur_state_change_high_rate_is_consistent(); }
TEST(CheckForUnimolecularReactions, ImplicitLipidIsSkippedForDissociation) { test_cfur_implicit_lipid_is_skipped_for_dissociation(); }
TEST(CheckForUnimolecularReactions, JustDissociatedMoleculeIsSkipped) { test_cfur_just_dissociated_molecule_is_skipped(); }