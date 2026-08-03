/*! \file test_check_for_unimolecular_reactions_population.cpp
 *
 * ### Unit test for src/reactions/check_for_unimolecular_reactions_population.cpp
 *
 * Function under test:
 * \code
 * void check_for_unimolecular_reactions_population(
 *          long long int simItr, Parameters&, std::vector<Molecule>&,
 *          std::vector<Complex>&, SimulVolume&,
 *          const std::vector<ForwardRxn>&, const std::vector<BackRxn>&,
 *          const std::vector<CreateDestructRxn>&, std::vector<MolTemplate>&,
 *          std::map<std::string,int>&, copyCounters&, Membrane&,
 *          std::ofstream&);
 * \endcode
 *
 * The routine performs *population based* (binomially sampled) unimolecular
 * events once per timestep.  It contains four independent blocks:
 *
 *   1. Explicit destruction         (A -> null,  A is a normal molecule)
 *   2. Implicit-lipid destruction   (A -> null,  A is the implicit lipid)
 *   3. Unimolecular creation        (A -> A + B)
 *   4. Explicit dissociation        (A-B -> A + B, driven by backRxns)
 *
 * The number of events per block is drawn from
 *      numEvents = Binomial(1 - exp(-rate*dt), N_available)
 * so the *deterministic* limits are easy to test:
 *
 *   - rate == 0                  -> prob == 0 -> numEvents == 0 (nothing happens)
 *   - rate * dt  >> 1            -> prob == 1 -> numEvents == N (everything happens)
 *   - rate == -1 (destruction)   -> special "strong destruction" branch where
 *                                   numEvents = monomerList.size() - copies
 *
 * Those three limits are what the tests below exercise.  The verbose stderr
 * output states which source file / function / branch is being probed and what
 * the pass criteria are.
 *
 * NOTE: the actual dissociation machinery (break_interaction()) requires a fully
 * bonded, self-consistent molecular system; it is therefore only exercised here
 * in the "no events" limit (rate == 0), which still covers the NAB counting, the
 * off-rate-consistency scan and the early-exit logic of block 4.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

// `r` (the global GSL random number generator) is defined by gtest_main.cpp.
// It is used by rand_gsl() and by gsl_ran_binomial() inside the function under
// test, so it must point at a valid generator before we call anything.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Small helper: make sure the global GSL RNG exists (gtest_main leaves it null
// unless a test initialises it).  Without this, gsl_ran_binomial() would
// dereference a null generator.
// -----------------------------------------------------------------------------
void cfurp_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] global GSL rng 'r' was null -> allocating mt19937 (seed 12345)\n";
        r = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(r, 12345);
    }
}

/*! \brief Bundles every argument the function under test needs.
 *
 * Keeping them together makes each test read as "build a system, fire one
 * reaction block at it, inspect the system".
 */
struct CfurpSystem {
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
    // Never written to in these tests (no dissociation is actually performed),
    // but the signature demands a stream.
    std::ofstream assocDissocFile {};
};

// -----------------------------------------------------------------------------
// Reset the class-level statics that Molecule::destroy() / Complex::destroy()
// touch.  Tests share one process, so stale statics from a previous test would
// otherwise leak in (or index out of range).
// -----------------------------------------------------------------------------
void cfurp_reset_statics(int numMolTypes, int numMolecules, int numComplexes)
{
    MolTemplate::numMolTypes = static_cast<unsigned>(numMolTypes);
    MolTemplate::numEachMolType.assign(numMolTypes, numMolecules);
    MolTemplate::absToRelIface.assign(1, 0);

    Molecule::numberOfMolecules = numMolecules;
    Molecule::emptyMolList.clear();

    Complex::numberOfComplexes = numComplexes;
    Complex::emptyComList.clear();
    Complex::currNumberMolTypes = numMolTypes;
    Complex::currNumberComTypes = numMolTypes;
    Complex::obs.clear();
}

/*! \brief Build a minimal but self-consistent system of two free monomers.
 *
 * Layout:
 *   - one MolTemplate ("A", molTypeIndex 0) with one interface (absolute
 *     species index 0) and monomerList = {0, 1}
 *   - two molecules (indices 0 and 1), each alone in its own complex
 *   - one SubVolume containing both molecules
 *   - copyNumSpecies = {2, 0, 0}  (species 0 = free interface 'a')
 *
 * \param[out] sys           system to populate
 * \param[in]  implicitLipid mark the template as the implicit lipid (selects the
 *                           second destruction branch of the function)
 */
void cfurp_build_basic_system(CfurpSystem& sys, bool implicitLipid)
{
    cfurp_reset_statics(/*numMolTypes*/ 1, /*numMolecules*/ 2, /*numComplexes*/ 2);

    // ---- Parameters ---------------------------------------------------------
    sys.params.timeStep = 1.0; // microseconds; rate*dt*1e-6 == lambda
    sys.params.numMolTypes = 1;
    sys.params.numTotalSpecies = 3;
    sys.params.nItr = 1;

    // ---- Boundary ----------------------------------------------------------
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    sys.membrane.isBox = true;
    sys.membrane.isSphere = false;
    sys.membrane.implicitlipidIndex = -1;
    sys.membrane.numberOfFreeLipidsEachState = std::vector<int> { 0 };
    sys.membrane.numberOfProteinEachState = std::vector<int> { 0 };
    sys.membrane.nStates = 1;

    // ---- MolTemplate -------------------------------------------------------
    MolTemplate temp;
    temp.molName = "A";
    temp.molTypeIndex = 0;
    temp.copies = 2;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord(10.0, 10.0, 10.0);
    temp.Dr = Coord(0.1, 0.1, 0.1);
    temp.isImplicitLipid = implicitLipid;
    temp.canDestroy = true;

    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord(0.0, 0.0, 0.0);
    iface.stateList.push_back(Interface::State('\0', 0)); // header-inline ctor
    temp.interfaceList.push_back(iface);

    // Both molecules are free monomers and therefore candidates for
    // destruction / creation events.
    temp.monomerList = std::vector<int> { 0, 1 };
    sys.molTemplateList.push_back(temp);

    // ---- Molecules ---------------------------------------------------------
    for (int i = 0; i < 2; ++i) {
        Molecule mol;
        mol.index = i;
        mol.id = i;
        mol.molTypeIndex = 0;
        mol.myComIndex = i;
        mol.complexId = i;
        mol.mySubVolIndex = 0; // lives in the single sub-volume
        mol.mass = 1.0;
        mol.comCoord = Coord(5.0 * i, 0.0, 0.0);
        mol.trajStatus = TrajStatus::none;
        mol.isImplicitLipid = false;

        Molecule::Iface mIface;
        mIface.coord = mol.comCoord;
        mIface.index = 0; // absolute species index of the free interface
        mIface.relIndex = 0;
        mIface.molTypeIndex = 0;
        mIface.isBound = false;
        mol.interfaceList.push_back(mIface);

        sys.moleculeList.push_back(mol);
    }

    // ---- Complexes ---------------------------------------------------------
    for (int i = 0; i < 2; ++i) {
        Complex com;
        com.index = i;
        com.id = i;
        com.comCoord = sys.moleculeList[i].comCoord;
        com.memberList = std::vector<int> { i };
        com.numEachMol = std::vector<int> { 1 };
        com.lastNumberUpdateItrEachMol = std::vector<long long int> { 0 };
        com.mass = 1.0;
        com.radius = 1.0;
        com.D = Coord(10.0, 10.0, 10.0);
        com.Dr = Coord(0.1, 0.1, 0.1);
        com.trajStatus = TrajStatus::none;
        sys.complexList.push_back(com);
    }

    // ---- SimulVolume (one cell holding both molecules) ---------------------
    sys.simulVolume.numSubCells.x = 1;
    sys.simulVolume.numSubCells.y = 1;
    sys.simulVolume.numSubCells.z = 1;
    sys.simulVolume.numSubCells.tot = 1;
    sys.simulVolume.subCellSize = Coord(100.0, 100.0, 100.0);
    SimulVolume::SubVolume sub;
    sub.absIndex = 0;
    sub.xIndex = 0;
    sub.yIndex = 0;
    sub.zIndex = 0;
    sub.memberMolList = std::vector<int> { 0, 1 };
    sys.simulVolume.subCellList.push_back(sub);

    // ---- copyCounters ------------------------------------------------------
    // species 0 -> free interface 'a' (2 copies), species 1 -> bound pair.
    sys.counterArrays.copyNumSpecies = std::vector<int> { 2, 0, 0 };
    sys.counterArrays.bindPairList.resize(3);
    sys.counterArrays.nBoundPairs = std::vector<int> { 0 };
    sys.counterArrays.proPairlist = std::vector<int> { 0 };
    sys.counterArrays.nLoops = 0;
}

/*! \brief Build a destruction reaction (A -> null) with the requested rate. */
CreateDestructRxn cfurp_make_destruction_rxn(double rate, bool isObserved, const std::string& label)
{
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::destruction;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isObserved = isObserved;
    rxn.observeLabel = label;

    CreateDestructRxn::CreateDestructMol reactant;
    reactant.molTypeIndex = 0;
    reactant.molName = "A";
    reactant.interfaceList.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    rxn.reactantMolList.push_back(reactant);

    rxn.rateList.push_back(RxnBase::RateState(rate, {}));
    return rxn;
}

/*! \brief Build a unimolecular creation reaction (A -> A + A). */
CreateDestructRxn cfurp_make_uni_creation_rxn(double rate)
{
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::uniMolCreation;
    rxn.absRxnIndex = 1;
    rxn.relRxnIndex = 1;
    rxn.isObserved = false;
    rxn.creationRadius = 2.0;

    CreateDestructRxn::CreateDestructMol reactant;
    reactant.molTypeIndex = 0;
    reactant.molName = "A";
    reactant.interfaceList.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    rxn.reactantMolList.push_back(reactant);

    // Product list: parent + newly created molecule (same type here).
    rxn.productMolList.push_back(reactant);
    rxn.productMolList.push_back(reactant);

    rxn.rateList.push_back(RxnBase::RateState(rate, {}));
    return rxn;
}

/*! \brief Build a matched ForwardRxn/BackRxn pair describing A(a)+A(a) <-> A-A.
 *
 * The bound species has absolute index 1, the free interfaces index 0.
 */
void cfurp_make_dissociation_pair(CfurpSystem& sys, double offRate)
{
    ForwardRxn fwd;
    fwd.rxnType = ReactionType::bimolecular;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.bindRadius = 1.0;
    fwd.reactantListNew.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    fwd.reactantListNew.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    fwd.productListNew.push_back(RxnIface("a", 0, 1, 0, '\0', true));
    fwd.productListNew.push_back(RxnIface("a", 0, 1, 0, '\0', true));
    fwd.rateList.push_back(RxnBase::RateState(1.0, {}));
    sys.forwardRxns.push_back(fwd);

    BackRxn back;
    back.rxnType = ReactionType::bimolecular;
    back.absRxnIndex = 1;
    back.relRxnIndex = 0;
    back.conjForwardRxnIndex = 0;
    back.isCoupled = false;
    back.isObserved = false;
    // Reactant of the back reaction == the bound species (absIfaceIndex 1)
    back.reactantListNew.push_back(RxnIface("a", 0, 1, 0, '\0', true));
    back.reactantListNew.push_back(RxnIface("a", 0, 1, 0, '\0', true));
    // Products are the two free interfaces (absIfaceIndex 0)
    back.productListNew.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    back.productListNew.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    back.rateList.push_back(RxnBase::RateState(offRate, {}));
    sys.backRxns.push_back(back);
}

} // namespace

// =============================================================================
// Test 1: destruction block, rate == 0  -> no molecule may be destroyed.
// =============================================================================
void test_cfurp_destruction_zero_rate()
{
    std::cerr << "\n[TEST] test_cfurp_destruction_zero_rate\n"
              << "  Source file: check_for_unimolecular_reactions_population.cpp\n"
              << "  Branch:      explicit destruction (A -> null), rate = 0\n"
              << "  Criteria:    Binomial(p=0, N=2) == 0 events, so the two monomers,\n"
              << "               their complexes, the monomerList and copyNumSpecies\n"
              << "               must all be untouched.\n";

    cfurp_ensure_rng();

    CfurpSystem sys;
    cfurp_build_basic_system(sys, /*implicitLipid*/ false);
    sys.createDestructRxns.push_back(cfurp_make_destruction_rxn(0.0, false, ""));

    std::cerr << "  Calling check_for_unimolecular_reactions_population(simItr=1)...\n";
    check_for_unimolecular_reactions_population(1, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.membrane, sys.assocDissocFile);

    EXPECT_EQ(sys.molTemplateList[0].monomerList.size(), 2u)
        << "monomerList must still contain both monomers at zero rate";
    EXPECT_FALSE(sys.moleculeList[0].isEmpty) << "molecule 0 must not be destroyed";
    EXPECT_FALSE(sys.moleculeList[1].isEmpty) << "molecule 1 must not be destroyed";
    EXPECT_FALSE(sys.complexList[0].isEmpty) << "complex 0 must not be destroyed";
    EXPECT_FALSE(sys.complexList[1].isEmpty) << "complex 1 must not be destroyed";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 2)
        << "copy number of the free interface must stay 2";
    EXPECT_EQ(sys.simulVolume.subCellList[0].memberMolList.size(), 2u)
        << "no molecule may be removed from the sub-volume member list";

    std::cerr << "  monomerList size = " << sys.molTemplateList[0].monomerList.size()
              << ", copyNumSpecies[0] = " << sys.counterArrays.copyNumSpecies[0] << "\n";
}

// =============================================================================
// Test 2: destruction block, huge rate -> all monomers destroyed.
// =============================================================================
void test_cfurp_destruction_all_events()
{
    std::cerr << "\n[TEST] test_cfurp_destruction_all_events\n"
              << "  Source file: check_for_unimolecular_reactions_population.cpp\n"
              << "  Branch:      explicit destruction (A -> null), rate*dt >> 1\n"
              << "  Criteria:    Binomial(p~1, N=2) == 2 events, so both molecules and\n"
              << "               complexes are flagged empty, copyNumSpecies[0] -> 0,\n"
              << "               monomerList empties and mySubVolIndex is reset to -1.\n";

    cfurp_ensure_rng();

    CfurpSystem sys;
    cfurp_build_basic_system(sys, /*implicitLipid*/ false);
    // lambda = 1e12 (1/s) * 1 us * 1e-6 = 1e6 -> p = 1 - exp(-1e6) = 1
    sys.createDestructRxns.push_back(cfurp_make_destruction_rxn(1.0e12, true, "Acount"));
    sys.observablesList["Acount"] = 2;

    std::cerr << "  Calling check_for_unimolecular_reactions_population(simItr=2)...\n";
    check_for_unimolecular_reactions_population(2, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.membrane, sys.assocDissocFile);

    EXPECT_TRUE(sys.molTemplateList[0].monomerList.empty())
        << "monomerList must be emptied once every monomer has been destroyed";
    EXPECT_TRUE(sys.moleculeList[0].isEmpty) << "molecule 0 should be flagged empty";
    EXPECT_TRUE(sys.moleculeList[1].isEmpty) << "molecule 1 should be flagged empty";
    EXPECT_TRUE(sys.complexList[0].isEmpty) << "complex 0 should be flagged empty";
    EXPECT_TRUE(sys.complexList[1].isEmpty) << "complex 1 should be flagged empty";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 0)
        << "both interface copies must be removed from copyNumSpecies";
    EXPECT_EQ(sys.moleculeList[0].mySubVolIndex, -1)
        << "destroyed molecule 0 must have its sub-volume index invalidated";
    EXPECT_EQ(sys.moleculeList[1].mySubVolIndex, -1)
        << "destroyed molecule 1 must have its sub-volume index invalidated";
    EXPECT_TRUE(sys.simulVolume.subCellList[0].memberMolList.empty())
        << "destroyed molecules must be erased from the sub-volume member list";
    EXPECT_EQ(sys.observablesList["Acount"], 0)
        << "the observable coupled to the destruction reaction must be decremented twice";

    std::cerr << "  monomerList size = " << sys.molTemplateList[0].monomerList.size()
              << ", copyNumSpecies[0] = " << sys.counterArrays.copyNumSpecies[0]
              << ", subCell members = " << sys.simulVolume.subCellList[0].memberMolList.size()
              << ", observable = " << sys.observablesList["Acount"] << "\n";
}

// =============================================================================
// Test 3: destruction block, rate == -1 -> "strong destruction" titration mode.
// =============================================================================
void test_cfurp_destruction_strong_titration()
{
    std::cerr << "\n[TEST] test_cfurp_destruction_strong_titration\n"
              << "  Source file: check_for_unimolecular_reactions_population.cpp\n"
              << "  Branch:      explicit destruction with the special rate == -1 flag\n"
              << "  Criteria:    numEvents = monomerList.size() - template.copies\n"
              << "               = 2 - 1 = 1, i.e. exactly ONE monomer is destroyed and\n"
              << "               the monomer pool shrinks from 2 to 1.\n";

    cfurp_ensure_rng();

    CfurpSystem sys;
    cfurp_build_basic_system(sys, /*implicitLipid*/ false);
    sys.molTemplateList[0].copies = 1; // target population is a single monomer
    sys.createDestructRxns.push_back(cfurp_make_destruction_rxn(-1.0, false, ""));

    std::cerr << "  Calling check_for_unimolecular_reactions_population(simItr=3)...\n";
    check_for_unimolecular_reactions_population(3, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.membrane, sys.assocDissocFile);

    EXPECT_EQ(sys.molTemplateList[0].monomerList.size(), 1u)
        << "strong destruction should trim the monomer pool down to 'copies' (1)";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 1)
        << "exactly one interface copy should have been removed";

    // Exactly one of the two molecules must have been destroyed.
    const int numDestroyed = (sys.moleculeList[0].isEmpty ? 1 : 0) + (sys.moleculeList[1].isEmpty ? 1 : 0);
    EXPECT_EQ(numDestroyed, 1) << "exactly one molecule must be destroyed in titration mode";
    EXPECT_EQ(sys.simulVolume.subCellList[0].memberMolList.size(), 1u)
        << "the destroyed molecule should have left the sub-volume member list";

    std::cerr << "  monomerList size = " << sys.molTemplateList[0].monomerList.size()
              << ", molecules destroyed = " << numDestroyed
              << ", copyNumSpecies[0] = " << sys.counterArrays.copyNumSpecies[0] << "\n";
}

// =============================================================================
// Test 4: implicit-lipid destruction branch (no explicit molecules involved).
// =============================================================================
void test_cfurp_implicit_lipid_destruction()
{
    std::cerr << "\n[TEST] test_cfurp_implicit_lipid_destruction\n"
              << "  Source file: check_for_unimolecular_reactions_population.cpp\n"
              << "  Branch:      destruction where the reactant IS the implicit lipid\n"
              << "  Criteria:    with rate*dt >> 1 all 5 free implicit lipids are removed,\n"
              << "               decrementing both membrane.numberOfFreeLipidsEachState[0]\n"
              << "               and counterArrays.copyNumSpecies[0] to zero, and the\n"
              << "               attached observable by 5.\n";

    cfurp_ensure_rng();

    CfurpSystem sys;
    cfurp_build_basic_system(sys, /*implicitLipid*/ true);
    sys.membrane.implicitLipid = true;
    sys.membrane.numberOfFreeLipidsEachState[0] = 5;
    sys.counterArrays.copyNumSpecies[0] = 5;

    sys.createDestructRxns.push_back(cfurp_make_destruction_rxn(1.0e12, true, "ILcount"));
    sys.observablesList["ILcount"] = 5;

    std::cerr << "  Before: freeLipids = " << sys.membrane.numberOfFreeLipidsEachState[0]
              << ", copyNumSpecies[0] = " << sys.counterArrays.copyNumSpecies[0] << "\n";
    std::cerr << "  Calling check_for_unimolecular_reactions_population(simItr=4)...\n";
    check_for_unimolecular_reactions_population(4, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.membrane, sys.assocDissocFile);

    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], 0)
        << "all free implicit lipids must be consumed at effectively unit probability";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 0)
        << "implicit lipid species copy number must be driven to zero";
    EXPECT_EQ(sys.observablesList["ILcount"], 0)
        << "the observable must be decremented once per destruction event";

    // Nothing explicit should have been touched by the implicit-lipid branch.
    EXPECT_FALSE(sys.moleculeList[0].isEmpty)
        << "the implicit-lipid branch must not destroy explicit molecules";
    EXPECT_EQ(sys.simulVolume.subCellList[0].memberMolList.size(), 2u)
        << "sub-volume membership must be untouched by the implicit-lipid branch";

    std::cerr << "  After:  freeLipids = " << sys.membrane.numberOfFreeLipidsEachState[0]
              << ", copyNumSpecies[0] = " << sys.counterArrays.copyNumSpecies[0]
              << ", observable = " << sys.observablesList["ILcount"] << "\n";
}

// =============================================================================
// Test 5: unimolecular creation block at zero rate -> no new molecules.
// =============================================================================
void test_cfurp_uni_creation_zero_rate()
{
    std::cerr << "\n[TEST] test_cfurp_uni_creation_zero_rate\n"
              << "  Source file: check_for_unimolecular_reactions_population.cpp\n"
              << "  Branch:      unimolecular creation (A -> A + B), rate = 0\n"
              << "  Criteria:    the candidate pool is still scanned via isReactant(), but\n"
              << "               Binomial(p=0,N) == 0 so moleculeList/complexList sizes,\n"
              << "               copy numbers and TrajStatus flags stay unchanged.\n";

    cfurp_ensure_rng();

    CfurpSystem sys;
    cfurp_build_basic_system(sys, /*implicitLipid*/ false);
    sys.createDestructRxns.push_back(cfurp_make_uni_creation_rxn(0.0));

    const size_t nMolBefore = sys.moleculeList.size();
    const size_t nComBefore = sys.complexList.size();

    std::cerr << "  Calling check_for_unimolecular_reactions_population(simItr=5)...\n";
    check_for_unimolecular_reactions_population(5, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.membrane, sys.assocDissocFile);

    EXPECT_EQ(sys.moleculeList.size(), nMolBefore)
        << "no molecule may be created when the creation rate is zero";
    EXPECT_EQ(sys.complexList.size(), nComBefore)
        << "no complex may be created when the creation rate is zero";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 2)
        << "copy numbers must be unchanged when no creation event fires";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "the parent molecule must not be marked propagated if it did not react";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[1].trajStatus), static_cast<int>(TrajStatus::none))
        << "the parent molecule must not be marked propagated if it did not react";

    std::cerr << "  moleculeList size = " << sys.moleculeList.size()
              << ", complexList size = " << sys.complexList.size() << "\n";
}

// =============================================================================
// Test 6: dissociation block with a bound pair present but zero off-rate.
// =============================================================================
void test_cfurp_dissociation_zero_rate()
{
    std::cerr << "\n[TEST] test_cfurp_dissociation_zero_rate\n"
              << "  Source file: check_for_unimolecular_reactions_population.cpp\n"
              << "  Branch:      explicit dissociation driven by backRxns, offRate = 0\n"
              << "  Criteria:    NAB (= bindPairList size) is 1 so the block is entered and\n"
              << "               the off-rate consistency scan runs, but Binomial(p=0,1) == 0\n"
              << "               events, so the bound pair, copy numbers, nLoops and the\n"
              << "               isDissociated flags remain untouched.\n";

    cfurp_ensure_rng();

    CfurpSystem sys;
    cfurp_build_basic_system(sys, /*implicitLipid*/ false);
    cfurp_make_dissociation_pair(sys, /*offRate*/ 0.0);

    // Bind molecule 0's interface 0 to molecule 1's interface 0 (species index 1).
    sys.moleculeList[0].interfaceList[0].isBound = true;
    sys.moleculeList[0].interfaceList[0].index = 1;
    sys.moleculeList[0].interfaceList[0].interaction = Molecule::Interaction(1, 0, 0);
    sys.moleculeList[1].interfaceList[0].isBound = true;
    sys.moleculeList[1].interfaceList[0].index = 1;
    sys.moleculeList[1].interfaceList[0].interaction = Molecule::Interaction(0, 0, 0);

    // Both molecules now live in one complex (index 0).
    sys.moleculeList[1].myComIndex = 0;
    sys.complexList[0].memberList = std::vector<int> { 0, 1 };
    sys.complexList[1].isEmpty = true;
    sys.complexList[1].memberList.clear();

    // Population bookkeeping: one A-A pair, no free 'a' interfaces.
    sys.counterArrays.copyNumSpecies[0] = 0;
    sys.counterArrays.copyNumSpecies[1] = 1;
    sys.counterArrays.bindPairList[1] = std::vector<int> { 0 };
    sys.counterArrays.nBoundPairs[0] = 1;

    std::cerr << "  Calling check_for_unimolecular_reactions_population(simItr=6)...\n";
    check_for_unimolecular_reactions_population(6, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.membrane, sys.assocDissocFile);

    EXPECT_EQ(sys.counterArrays.bindPairList[1].size(), 1u)
        << "the bound pair must remain in bindPairList when the off-rate is zero";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], 1)
        << "bound species copy number must be unchanged";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 0)
        << "no free interfaces may appear without a dissociation event";
    EXPECT_EQ(sys.counterArrays.nBoundPairs[0], 1)
        << "nBoundPairs must be unchanged when no dissociation occurs";
    EXPECT_EQ(sys.counterArrays.nLoops, 0) << "loop counter must be unchanged";
    EXPECT_FALSE(sys.moleculeList[0].isDissociated) << "molecule 0 must not be flagged dissociated";
    EXPECT_FALSE(sys.moleculeList[1].isDissociated) << "molecule 1 must not be flagged dissociated";
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound)
        << "interface of molecule 0 must still be bound";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, 1)
        << "the interaction partner bookkeeping must be preserved";

    std::cerr << "  bindPairList[1] size = " << sys.counterArrays.bindPairList[1].size()
              << ", copyNumSpecies = {" << sys.counterArrays.copyNumSpecies[0] << ", "
              << sys.counterArrays.copyNumSpecies[1] << "}\n";
}

// =============================================================================
// Test 7: nothing configured at all -> the routine must be a harmless no-op.
// =============================================================================
void test_cfurp_no_reactions_is_noop()
{
    std::cerr << "\n[TEST] test_cfurp_no_reactions_is_noop\n"
              << "  Source file: check_for_unimolecular_reactions_population.cpp\n"
              << "  Branch:      all four reaction lists empty\n"
              << "  Criteria:    the function must return without altering any molecule,\n"
              << "               complex, counter or membrane state.\n";

    cfurp_ensure_rng();

    CfurpSystem sys;
    cfurp_build_basic_system(sys, /*implicitLipid*/ false);
    // Deliberately leave createDestructRxns / backRxns / forwardRxns empty.

    std::cerr << "  Calling check_for_unimolecular_reactions_population(simItr=7)...\n";
    check_for_unimolecular_reactions_population(7, sys.params, sys.moleculeList, sys.complexList,
        sys.simulVolume, sys.forwardRxns, sys.backRxns, sys.createDestructRxns, sys.molTemplateList,
        sys.observablesList, sys.counterArrays, sys.membrane, sys.assocDissocFile);

    EXPECT_EQ(sys.moleculeList.size(), 2u) << "molecule count must not change";
    EXPECT_EQ(sys.complexList.size(), 2u) << "complex count must not change";
    EXPECT_EQ(sys.molTemplateList[0].monomerList.size(), 2u) << "monomer pool must not change";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 2) << "copy numbers must not change";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], 0)
        << "membrane lipid bookkeeping must not change";
    EXPECT_FALSE(sys.moleculeList[0].isEmpty) << "molecule 0 must survive";
    EXPECT_FALSE(sys.moleculeList[1].isEmpty) << "molecule 1 must survive";

    std::cerr << "  System unchanged: molecules = " << sys.moleculeList.size()
              << ", complexes = " << sys.complexList.size() << "\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so a failure
// in one scenario never prevents the others from executing.
// -----------------------------------------------------------------------------
TEST(CheckForUnimolecularReactionsPopulation, DestructionZeroRate) { test_cfurp_destruction_zero_rate(); }
TEST(CheckForUnimolecularReactionsPopulation, DestructionAllEvents) { test_cfurp_destruction_all_events(); }
TEST(CheckForUnimolecularReactionsPopulation, DestructionStrongTitration) { test_cfurp_destruction_strong_titration(); }
TEST(CheckForUnimolecularReactionsPopulation, ImplicitLipidDestruction) { test_cfurp_implicit_lipid_destruction(); }
TEST(CheckForUnimolecularReactionsPopulation, UniCreationZeroRate) { test_cfurp_uni_creation_zero_rate(); }
TEST(CheckForUnimolecularReactionsPopulation, DissociationZeroRate) { test_cfurp_dissociation_zero_rate(); }
TEST(CheckForUnimolecularReactionsPopulation, NoReactionsIsNoOp) { test_cfurp_no_reactions_is_noop(); }