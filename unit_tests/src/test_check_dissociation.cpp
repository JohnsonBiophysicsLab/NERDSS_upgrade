/*! \file test_check_dissociation.cpp
 *
 * ### Unit test for src/reactions/check_dissociation.cpp
 *
 * Function under test:
 *
 *     void check_dissociation(unsigned int simItr, const Parameters&, SimulVolume&,
 *                             std::vector<MolTemplate>&, std::map<std::string,int>&,
 *                             unsigned int molItr, std::vector<Molecule>&,
 *                             std::vector<Complex>&, const std::vector<BackRxn>&,
 *                             const std::vector<ForwardRxn>&,
 *                             const std::vector<CreateDestructRxn>&,
 *                             copyCounters&, const Membrane&, std::ofstream&)
 *
 * check_dissociation() walks the bound-interface list of one molecule and, for
 * every bound partner, decides whether the bond should be broken this time step.
 * There are a large number of early-exit ("guard") branches before any bond is
 * actually broken:
 *
 *   1. the molecule already dissociated this step        -> break out of the loop
 *   2. the interface is not actually bound               -> skip
 *   3. molItr > partnerIndex (pair handled by partner)   -> skip
 *   4. partner is an implicit lipid                      -> skip
 *   5. partner already dissociated this step             -> skip
 *   6. the partner does not point back at us             -> skip
 *   7. conjBackRxn == -1 (irreversible forward rxn)      -> skip
 *   8. find_reaction_rate_state() returns -1             -> skip
 *
 * Only when all guards pass is the dissociation probability evaluated and, if it
 * fires, break_interaction() is called and the bookkeeping arrays
 * (copyNumSpecies, nBoundPairs, observables, trajStatus, isDissociated) updated.
 *
 * The tests below build a small, fully self-consistent two-molecule / one-complex
 * system (mol 0 of type 0 bound through interface 0 to mol 1 of type 1) and then:
 *   - drive every guard branch and verify that *nothing* changes, and
 *   - force a dissociation (params.debugParams.forceDissoc = true so that
 *     prob == 1.0) and verify all of the bookkeeping side effects.
 *
 * Verbose progress information is written to stderr so a reader of the test log
 * can see which source file/function is exercised and what the pass criteria are.
 */

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

#include <gtest/gtest.h>

namespace {

// -----------------------------------------------------------------------------
// Absolute interface-state indices used by the toy reaction
//   A(a) + B(b) <-> A(a!1).B(b!1)
// free A.a  -> 0     free B.b  -> 1
// bound A.a -> 2     bound B.b -> 3
// -----------------------------------------------------------------------------
constexpr int kCdFreeA = 0;
constexpr int kCdFreeB = 1;
constexpr int kCdBoundA = 2;
constexpr int kCdBoundB = 3;

/*! \brief RAII guard that saves/restores the process-wide static counters that
 *         the Molecule/Complex/MolTemplate classes keep.
 *
 * All generated tests are linked into one binary, so we must not leak modified
 * static state into the other test translation units.
 */
struct CdStaticGuard {
    unsigned savedNumMolTypes;
    std::vector<int> savedNumEachMolType;
    int savedNumberOfMolecules;
    std::vector<int> savedEmptyMolList;
    int savedNumberOfComplexes;
    std::vector<int> savedEmptyComList;
    int savedCurrNumberMolTypes;

    CdStaticGuard()
        : savedNumMolTypes(MolTemplate::numMolTypes)
        , savedNumEachMolType(MolTemplate::numEachMolType)
        , savedNumberOfMolecules(Molecule::numberOfMolecules)
        , savedEmptyMolList(Molecule::emptyMolList)
        , savedNumberOfComplexes(Complex::numberOfComplexes)
        , savedEmptyComList(Complex::emptyComList)
        , savedCurrNumberMolTypes(Complex::currNumberMolTypes)
    {
    }

    ~CdStaticGuard()
    {
        MolTemplate::numMolTypes = savedNumMolTypes;
        MolTemplate::numEachMolType = savedNumEachMolType;
        Molecule::numberOfMolecules = savedNumberOfMolecules;
        Molecule::emptyMolList = savedEmptyMolList;
        Complex::numberOfComplexes = savedNumberOfComplexes;
        Complex::emptyComList = savedEmptyComList;
        Complex::currNumberMolTypes = savedCurrNumberMolTypes;
    }
};

/*! \brief Everything check_dissociation() needs, kept alive for the call. */
struct CdSystem {
    Parameters params;
    Membrane membraneObject;
    SimulVolume simulVolume;
    std::vector<MolTemplate> molTemplateList;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<BackRxn> backRxns;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<CreateDestructRxn> createDestructRxns;
    copyCounters counterArrays;
    std::map<std::string, int> observablesList;
    std::ofstream assocDissocFile;
};

/*! \brief Make sure the global GSL RNG exists.
 *
 * check_dissociation() always calls rand_gsl() once it reaches the probability
 * test, and the suite-wide `r` pointer starts out as nullptr, so it must be
 * initialised before the first call or the process would crash.
 */
void cd_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  (initialising global GSL RNG for rand_gsl())\n";
        srand_gsl(1);
    }
}

/*! \brief Build one MolTemplate with a single interface that has a single state. */
MolTemplate cd_make_template(int typeIndex, const std::string& molName, const std::string& ifaceName, int freeAbsIndex)
{
    MolTemplate mt;
    mt.molName = molName;
    mt.molTypeIndex = typeIndex;
    mt.copies = 1;
    mt.mass = 1.0;
    mt.radius = 1.0;
    mt.D = Coord { 1.0, 1.0, 1.0 };
    mt.Dr = Coord { 0.01, 0.01, 0.01 };
    mt.isLipid = false;
    mt.isImplicitLipid = false;
    mt.isRod = false;
    mt.isPoint = false;
    mt.checkOverlap = false;
    mt.excludeVolumeBound = false;

    Interface iface;
    iface.index = 0;
    iface.name = ifaceName;
    iface.iCoord = Coord { 1.0, 0.0, 0.0 };

    Interface::State state;
    state.index = freeAbsIndex; // absolute index of the *free* state
    state.iden = '\0';
    state.ifaceAndStateName = ifaceName;
    iface.stateList.push_back(state);

    mt.interfaceList.push_back(iface);
    return mt;
}

/*! \brief Build a molecule whose only interface is bound to `partnerIndex`. */
Molecule cd_make_bound_molecule(
    int index, int molTypeIndex, const Coord& com, int partnerIndex, int boundAbsIndex, int backRxnIndex)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = 0; // both molecules start in complex 0
    mol.mySubVolIndex = 0;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isLipid = false;
    mol.isEmpty = false;
    mol.isImplicitLipid = false;
    mol.isDissociated = false;
    mol.trajStatus = TrajStatus::none;
    mol.linksToSurface = 0;

    Molecule::Iface iface;
    iface.coord = Coord { com.x + 1.0, com.y, com.z };
    iface.index = boundAbsIndex; // currently in its *bound* state
    iface.relIndex = 0;
    iface.molTypeIndex = molTypeIndex;
    iface.stateIndex = 0;
    iface.stateIden = '\0';
    iface.isBound = true;
    iface.excludeVolume = false;
    iface.interaction = Molecule::Interaction(partnerIndex, 0, backRxnIndex);
    mol.interfaceList.push_back(iface);

    // legacy bookkeeping lists used by break_interaction()
    mol.bndlist.push_back(0);
    mol.bndpartner.push_back(partnerIndex);
    mol.bndRxnList.push_back(backRxnIndex);
    mol.freelist.clear();

    return mol;
}

/*! \brief Populate a CdSystem with the A(a!1).B(b!1) dimer described above. */
void cd_build_system(CdSystem& sys)
{
    // ---- statics -----------------------------------------------------------
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    Molecule::numberOfMolecules = 2;
    Molecule::emptyMolList.clear();
    Complex::numberOfComplexes = 1;
    Complex::currNumberMolTypes = 2;

    // ---- parameters --------------------------------------------------------
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 4;
    sys.params.timeStep = 1.0; // microseconds
    sys.params.debugParams.forceDissoc = false;
    sys.params.assocDissocWrite = false;

    // ---- boundary ----------------------------------------------------------
    sys.membraneObject.isBox = true;
    sys.membraneObject.isSphere = false;
    sys.membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    sys.membraneObject.implicitlipidIndex = -1;
    sys.membraneObject.implicitLipid = false;

    // ---- one (unused but valid) sub volume ---------------------------------
    SimulVolume::SubVolume subVol;
    subVol.absIndex = 0;
    subVol.xIndex = 0;
    subVol.yIndex = 0;
    subVol.zIndex = 0;
    subVol.memberMolList = std::vector<int> { 0, 1 };
    sys.simulVolume.subCellList.clear();
    sys.simulVolume.subCellList.push_back(subVol);

    // ---- templates ---------------------------------------------------------
    sys.molTemplateList.clear();
    sys.molTemplateList.push_back(cd_make_template(0, "A", "a", kCdFreeA));
    sys.molTemplateList.push_back(cd_make_template(1, "B", "b", kCdFreeB));

    // ---- molecules ---------------------------------------------------------
    sys.moleculeList.clear();
    sys.moleculeList.push_back(cd_make_bound_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, 1, kCdBoundA, 0));
    sys.moleculeList.push_back(cd_make_bound_molecule(1, 1, Coord { 2.0, 0.0, 0.0 }, 0, kCdBoundB, 0));

    // ---- complexes ---------------------------------------------------------
    sys.complexList.clear();

    Complex dimer;
    dimer.index = 0;
    dimer.id = 0;
    dimer.comCoord = Coord { 1.0, 0.0, 0.0 };
    dimer.memberList = std::vector<int> { 0, 1 };
    dimer.numEachMol = std::vector<int> { 1, 1 };
    dimer.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    dimer.D = Coord { 0.5, 0.5, 0.5 };
    dimer.Dr = Coord { 0.005, 0.005, 0.005 };
    dimer.radius = 2.0;
    dimer.mass = 2.0;
    dimer.isEmpty = false;
    dimer.OnSurface = false;
    dimer.trajStatus = TrajStatus::none;
    dimer.linksToSurface = 0;
    sys.complexList.push_back(dimer);

    // A spare, empty complex slot so that break_interaction() can re-use it when
    // the dimer splits (this mirrors what the real simulation looks like).
    Complex spare;
    spare.index = 1;
    spare.id = 1;
    spare.isEmpty = true;
    spare.memberList.clear();
    spare.numEachMol = std::vector<int> { 0, 0 };
    spare.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    sys.complexList.push_back(spare);
    Complex::emptyComList = std::vector<int> { 1 };

    // ---- reactions ---------------------------------------------------------
    std::vector<std::vector<RxnIface>> noOtherIfaces {};

    ForwardRxn fwd;
    fwd.rxnType = ReactionType::bimolecular;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.bindRadius = 1.0;
    fwd.isObserved = false;
    fwd.reactantListNew.push_back(RxnIface("a", 0, kCdFreeA, 0, '\0', false));
    fwd.reactantListNew.push_back(RxnIface("b", 1, kCdFreeB, 0, '\0', false));
    fwd.productListNew.push_back(RxnIface("a", 0, kCdBoundA, 0, '\0', true));
    fwd.productListNew.push_back(RxnIface("b", 1, kCdBoundB, 0, '\0', true));
    fwd.rateList.emplace_back(1.0, noOtherIfaces);
    fwd.intReactantList = std::vector<int> { kCdFreeA, kCdFreeB };
    fwd.intProductList = std::vector<int> { kCdBoundA, kCdBoundB };
    sys.forwardRxns.clear();
    sys.forwardRxns.push_back(fwd);

    BackRxn back;
    back.rxnType = ReactionType::bimolecular;
    back.absRxnIndex = 1;
    back.relRxnIndex = 0;
    back.conjForwardRxnIndex = 0;
    back.isCoupled = false;
    back.isObserved = false;
    back.isOnMem = false;
    back.isSymmetric = false;
    back.hasStateChange = false;
    // For a BackRxn the reactants are the *bound* species and the products the
    // two free species.
    back.reactantListNew.push_back(RxnIface("a", 0, kCdBoundA, 0, '\0', true));
    back.reactantListNew.push_back(RxnIface("b", 1, kCdBoundB, 0, '\0', true));
    back.productListNew.push_back(RxnIface("a", 0, kCdFreeA, 0, '\0', false));
    back.productListNew.push_back(RxnIface("b", 1, kCdFreeB, 0, '\0', false));
    back.rateList.emplace_back(1.0, noOtherIfaces); // 1 s^-1
    back.intReactantList = std::vector<int> { kCdBoundA, kCdBoundB };
    back.intProductList = std::vector<int> { kCdFreeA, kCdFreeB };
    sys.backRxns.clear();
    sys.backRxns.push_back(back);

    sys.createDestructRxns.clear();

    // ---- counters ----------------------------------------------------------
    sys.counterArrays.copyNumSpecies.assign(8, 0);
    sys.counterArrays.copyNumSpecies[kCdBoundA] = 1;
    sys.counterArrays.copyNumSpecies[kCdBoundB] = 1;
    sys.counterArrays.nBoundPairs.assign(4, 1);
    sys.counterArrays.proPairlist.assign(4, 0);
    sys.counterArrays.singleDouble.assign(8, 0);
    sys.counterArrays.implicitDouble.assign(8, false);
    sys.counterArrays.canDissociate.assign(8, false);
    sys.counterArrays.bindPairList.resize(8);
    sys.counterArrays.nLoops = 0;

    sys.observablesList.clear();

    // ---- output stream (contents are irrelevant for the test) -------------
    if (!sys.assocDissocFile.is_open())
        sys.assocDissocFile.open("/dev/null");
}

/*! \brief Convenience wrapper so every call site reads identically. */
void cd_call(CdSystem& sys, unsigned int molItr)
{
    check_dissociation(/*simItr=*/1, sys.params, sys.simulVolume, sys.molTemplateList, sys.observablesList, molItr,
        sys.moleculeList, sys.complexList, sys.backRxns, sys.forwardRxns, sys.createDestructRxns, sys.counterArrays,
        sys.membraneObject, sys.assocDissocFile);
}

/*! \brief Assert that no dissociation happened.
 *
 * Pass criteria: neither molecule is flagged dissociated, neither molecule was
 * marked as propagated, and the species / bound-pair counters are byte-for-byte
 * the same as the pre-call snapshot.
 */
void cd_expect_no_dissociation(const CdSystem& sys, const std::vector<int>& speciesBefore,
    const std::vector<int>& pairsBefore, const std::string& context)
{
    EXPECT_FALSE(sys.moleculeList[0].isDissociated) << context << ": molecule 0 must not be flagged dissociated";
    EXPECT_FALSE(sys.moleculeList[1].isDissociated) << context << ": molecule 1 must not be flagged dissociated";
    EXPECT_TRUE(sys.moleculeList[0].trajStatus == TrajStatus::none)
        << context << ": molecule 0 trajStatus must stay 'none' (got " << static_cast<int>(sys.moleculeList[0].trajStatus)
        << ')';
    EXPECT_TRUE(sys.moleculeList[1].trajStatus == TrajStatus::none)
        << context << ": molecule 1 trajStatus must stay 'none' (got " << static_cast<int>(sys.moleculeList[1].trajStatus)
        << ')';
    EXPECT_EQ(sys.counterArrays.copyNumSpecies, speciesBefore) << context << ": copyNumSpecies must be untouched";
    EXPECT_EQ(sys.counterArrays.nBoundPairs, pairsBefore) << context << ": nBoundPairs must be untouched";
}

} // namespace

// -----------------------------------------------------------------------------
// Guard 1: an empty bndlist means the for-loop body never runs.
// -----------------------------------------------------------------------------
void cd_test_empty_bndlist_is_a_noop()
{
    std::cerr << "\n[TEST] cd_test_empty_bndlist_is_a_noop\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : molecule 0 has no entries in bndlist.\n"
              << "  Criteria    : loop body never executes -> system unchanged.\n";

    CdStaticGuard guard;
    CdSystem sys;
    cd_build_system(sys);
    cd_ensure_rng();

    // Empty the bound-interface list; the interface itself stays bound so we can
    // be sure it is the list (and not the isBound flag) that stops the loop.
    sys.moleculeList[0].bndlist.clear();
    sys.params.debugParams.forceDissoc = true; // would fire if we got that far

    const std::vector<int> speciesBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> pairsBefore = sys.counterArrays.nBoundPairs;

    cd_call(sys, 0);

    cd_expect_no_dissociation(sys, speciesBefore, pairsBefore, "empty bndlist");
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "interface must remain bound";
}

// -----------------------------------------------------------------------------
// Guard 2: the interface listed in bndlist is not flagged isBound.
// -----------------------------------------------------------------------------
void cd_test_unbound_interface_is_skipped()
{
    std::cerr << "\n[TEST] cd_test_unbound_interface_is_skipped\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : bndlist references interface 0, but isBound==false.\n"
              << "  Criteria    : the 'if (isBound)' guard rejects it -> no changes.\n";

    CdStaticGuard guard;
    CdSystem sys;
    cd_build_system(sys);
    cd_ensure_rng();

    sys.moleculeList[0].interfaceList[0].isBound = false;
    sys.params.debugParams.forceDissoc = true;

    const std::vector<int> speciesBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> pairsBefore = sys.counterArrays.nBoundPairs;

    cd_call(sys, 0);

    cd_expect_no_dissociation(sys, speciesBefore, pairsBefore, "interface not bound");
}

// -----------------------------------------------------------------------------
// Guard 3: molItr > partnerIndex -- the pair is handled from the other side.
// -----------------------------------------------------------------------------
void cd_test_higher_index_molecule_defers_to_partner()
{
    std::cerr << "\n[TEST] cd_test_higher_index_molecule_defers_to_partner\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : called with molItr = 1 whose partner index is 0.\n"
              << "  Criteria    : 'molItr > pro2Index' -> continue, so the pair is\n"
              << "                only ever considered once (from molecule 0).\n";

    CdStaticGuard guard;
    CdSystem sys;
    cd_build_system(sys);
    cd_ensure_rng();

    sys.params.debugParams.forceDissoc = true;

    const std::vector<int> speciesBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> pairsBefore = sys.counterArrays.nBoundPairs;

    cd_call(sys, 1); // higher index -> must defer

    cd_expect_no_dissociation(sys, speciesBefore, pairsBefore, "molItr > partner");
    EXPECT_TRUE(sys.moleculeList[1].interfaceList[0].isBound) << "bond must survive";
}

// -----------------------------------------------------------------------------
// Guard 4: the partner is an implicit lipid (handled elsewhere).
// -----------------------------------------------------------------------------
void cd_test_implicit_lipid_partner_is_skipped()
{
    std::cerr << "\n[TEST] cd_test_implicit_lipid_partner_is_skipped\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : partner molecule has isImplicitLipid == true.\n"
              << "  Criteria    : implicit-lipid bonds are handled by\n"
              << "                check_dissociation_implicitlipid() -> no changes here.\n";

    CdStaticGuard guard;
    CdSystem sys;
    cd_build_system(sys);
    cd_ensure_rng();

    sys.moleculeList[1].isImplicitLipid = true;
    sys.params.debugParams.forceDissoc = true;

    const std::vector<int> speciesBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> pairsBefore = sys.counterArrays.nBoundPairs;

    cd_call(sys, 0);

    cd_expect_no_dissociation(sys, speciesBefore, pairsBefore, "implicit lipid partner");
}

// -----------------------------------------------------------------------------
// Guard 1b/5: molecules already flagged as dissociated this time step.
// -----------------------------------------------------------------------------
void cd_test_already_dissociated_molecules_are_skipped()
{
    std::cerr << "\n[TEST] cd_test_already_dissociated_molecules_are_skipped\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : (a) the queried molecule already dissociated,\n"
              << "                (b) the partner already dissociated.\n"
              << "  Criteria    : both cases leave all counters untouched.\n";

    CdStaticGuard guard;

    // (a) the molecule itself already dissociated -> 'break' at loop top
    {
        CdSystem sys;
        cd_build_system(sys);
        cd_ensure_rng();
        sys.moleculeList[0].isDissociated = true;
        sys.params.debugParams.forceDissoc = true;

        const std::vector<int> speciesBefore = sys.counterArrays.copyNumSpecies;
        const std::vector<int> pairsBefore = sys.counterArrays.nBoundPairs;

        std::cerr << "  -> case (a): moleculeList[0].isDissociated = true\n";
        cd_call(sys, 0);

        EXPECT_EQ(sys.counterArrays.copyNumSpecies, speciesBefore) << "case (a): copyNumSpecies untouched";
        EXPECT_EQ(sys.counterArrays.nBoundPairs, pairsBefore) << "case (a): nBoundPairs untouched";
        EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "case (a): bond must survive";
        EXPECT_FALSE(sys.moleculeList[1].isDissociated) << "case (a): partner untouched";
    }

    // (b) the partner already dissociated -> 'continue'
    {
        CdSystem sys;
        cd_build_system(sys);
        cd_ensure_rng();
        sys.moleculeList[1].isDissociated = true;
        sys.params.debugParams.forceDissoc = true;

        const std::vector<int> speciesBefore = sys.counterArrays.copyNumSpecies;
        const std::vector<int> pairsBefore = sys.counterArrays.nBoundPairs;

        std::cerr << "  -> case (b): moleculeList[1].isDissociated = true\n";
        cd_call(sys, 0);

        EXPECT_EQ(sys.counterArrays.copyNumSpecies, speciesBefore) << "case (b): copyNumSpecies untouched";
        EXPECT_EQ(sys.counterArrays.nBoundPairs, pairsBefore) << "case (b): nBoundPairs untouched";
        EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "case (b): bond must survive";
        EXPECT_FALSE(sys.moleculeList[0].isDissociated) << "case (b): queried molecule untouched";
    }
}

// -----------------------------------------------------------------------------
// Guard 6: the partner's interaction record does not point back at us.
// -----------------------------------------------------------------------------
void cd_test_inconsistent_partner_record_prevents_dissociation()
{
    std::cerr << "\n[TEST] cd_test_inconsistent_partner_record_prevents_dissociation\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : molecule 1 does not list molecule 0 as its partner.\n"
              << "  Criteria    : candissociate stays false -> no dissociation.\n";

    CdStaticGuard guard;
    CdSystem sys;
    cd_build_system(sys);
    cd_ensure_rng();

    // Break the reciprocity of the interaction records.
    sys.moleculeList[1].interfaceList[0].interaction.partnerIndex = -1;
    sys.params.debugParams.forceDissoc = true;

    const std::vector<int> speciesBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> pairsBefore = sys.counterArrays.nBoundPairs;

    cd_call(sys, 0);

    cd_expect_no_dissociation(sys, speciesBefore, pairsBefore, "non-reciprocal interaction");
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "bond must survive";
}

// -----------------------------------------------------------------------------
// Guard 7: irreversible forward reaction (conjBackRxn == -1).
// -----------------------------------------------------------------------------
void cd_test_irreversible_bond_never_dissociates()
{
    std::cerr << "\n[TEST] cd_test_irreversible_bond_never_dissociates\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : interaction.conjBackRxn == -1 (irreversible bond).\n"
              << "  Criteria    : 'if (mu == -1) continue' -> nothing happens even\n"
              << "                with params.debugParams.forceDissoc == true.\n";

    CdStaticGuard guard;
    CdSystem sys;
    cd_build_system(sys);
    cd_ensure_rng();

    sys.moleculeList[0].interfaceList[0].interaction.conjBackRxn = -1;
    sys.params.debugParams.forceDissoc = true;

    const std::vector<int> speciesBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> pairsBefore = sys.counterArrays.nBoundPairs;

    cd_call(sys, 0);

    cd_expect_no_dissociation(sys, speciesBefore, pairsBefore, "irreversible bond");
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "irreversible bond must survive";
}

// -----------------------------------------------------------------------------
// The full dissociation path, forced with debugParams.forceDissoc.
// -----------------------------------------------------------------------------
void cd_test_forced_dissociation_updates_all_bookkeeping()
{
    std::cerr << "\n[TEST] cd_test_forced_dissociation_updates_all_bookkeeping\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : a valid, reversible A(a!1).B(b!1) bond with\n"
              << "                params.debugParams.forceDissoc == true (prob = 1).\n"
              << "  Criteria    : - both interfaces become unbound\n"
              << "                - both molecules flagged isDissociated\n"
              << "                - both molecules/complexes flagged 'propagated'\n"
              << "                - copyNumSpecies: bound-- , free++ , free++\n"
              << "                - nBoundPairs for the (A,B) pair decremented\n";

    CdStaticGuard guard;
    CdSystem sys;
    cd_build_system(sys);
    cd_ensure_rng();

    sys.params.debugParams.forceDissoc = true; // prob forced to 1.0

    std::cerr << "  Before: copyNumSpecies[boundA]=" << sys.counterArrays.copyNumSpecies[kCdBoundA]
              << " copyNumSpecies[freeA]=" << sys.counterArrays.copyNumSpecies[kCdFreeA]
              << " copyNumSpecies[freeB]=" << sys.counterArrays.copyNumSpecies[kCdFreeB]
              << " nBoundPairs[0*2+1]=" << sys.counterArrays.nBoundPairs[1] << '\n';

    cd_call(sys, 0);

    std::cerr << "  After : copyNumSpecies[boundA]=" << sys.counterArrays.copyNumSpecies[kCdBoundA]
              << " copyNumSpecies[freeA]=" << sys.counterArrays.copyNumSpecies[kCdFreeA]
              << " copyNumSpecies[freeB]=" << sys.counterArrays.copyNumSpecies[kCdFreeB]
              << " nBoundPairs[0*2+1]=" << sys.counterArrays.nBoundPairs[1] << '\n';
    std::cerr << "  Complex indices now: mol0 -> " << sys.moleculeList[0].myComIndex << ", mol1 -> "
              << sys.moleculeList[1].myComIndex << " (complexList.size() = " << sys.complexList.size() << ")\n";
    std::cerr << "  counterArrays.nLoops = " << sys.counterArrays.nLoops << " (informational only)\n";

    // --- interface state ---------------------------------------------------
    EXPECT_FALSE(sys.moleculeList[0].interfaceList[0].isBound) << "molecule 0 interface must be released";
    EXPECT_FALSE(sys.moleculeList[1].interfaceList[0].isBound) << "molecule 1 interface must be released";

    // --- per-molecule flags set directly by check_dissociation -------------
    EXPECT_TRUE(sys.moleculeList[0].isDissociated) << "molecule 0 must be flagged isDissociated";
    EXPECT_TRUE(sys.moleculeList[1].isDissociated) << "molecule 1 must be flagged isDissociated";
    EXPECT_TRUE(sys.moleculeList[0].trajStatus == TrajStatus::propagated)
        << "molecule 0 must be flagged propagated (got " << static_cast<int>(sys.moleculeList[0].trajStatus) << ')';
    EXPECT_TRUE(sys.moleculeList[1].trajStatus == TrajStatus::propagated)
        << "molecule 1 must be flagged propagated (got " << static_cast<int>(sys.moleculeList[1].trajStatus) << ')';

    // --- complex flags -----------------------------------------------------
    const int com0 = sys.moleculeList[0].myComIndex;
    const int com1 = sys.moleculeList[1].myComIndex;
    ASSERT_GE(com0, 0);
    ASSERT_LT(com0, static_cast<int>(sys.complexList.size()));
    ASSERT_GE(com1, 0);
    ASSERT_LT(com1, static_cast<int>(sys.complexList.size()));
    EXPECT_TRUE(sys.complexList[com0].trajStatus == TrajStatus::propagated)
        << "complex of molecule 0 must be flagged propagated";
    EXPECT_TRUE(sys.complexList[com1].trajStatus == TrajStatus::propagated)
        << "complex of molecule 1 must be flagged propagated";
    EXPECT_NE(com0, com1) << "the dimer was only held together by this one bond, so the two "
                             "molecules must now live in different complexes";

    // --- copy-number bookkeeping ------------------------------------------
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCdBoundA], 0) << "bound species (reactantListNew[0]) must be decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCdFreeA], 1) << "free A.a (productListNew[0]) must be incremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kCdFreeB], 1) << "free B.b (productListNew[1]) must be incremented";

    // --- bound-pair bookkeeping -------------------------------------------
    EXPECT_EQ(sys.counterArrays.nBoundPairs[1], 0)
        << "update_Nboundpairs() must remove the single (type 0, type 1) bound pair";
}

// -----------------------------------------------------------------------------
// The reaction observable must be decremented when the bond breaks.
// -----------------------------------------------------------------------------
void cd_test_forced_dissociation_decrements_observable()
{
    std::cerr << "\n[TEST] cd_test_forced_dissociation_decrements_observable\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : backRxn.isObserved == true with a label present in\n"
              << "                observablesList, dissociation forced.\n"
              << "  Criteria    : the observable count drops by exactly one; an\n"
              << "                unrelated observable is left alone.\n";

    CdStaticGuard guard;
    CdSystem sys;
    cd_build_system(sys);
    cd_ensure_rng();

    sys.backRxns[0].isObserved = true;
    sys.backRxns[0].observeLabel = "AB";
    sys.observablesList["AB"] = 5;
    sys.observablesList["other"] = 7;
    sys.params.debugParams.forceDissoc = true;

    cd_call(sys, 0);

    std::cerr << "  observable \"AB\" = " << sys.observablesList["AB"] << " (expected 4), \"other\" = "
              << sys.observablesList["other"] << " (expected 7)\n";

    EXPECT_EQ(sys.observablesList["AB"], 4) << "observed product count must be decremented by the dissociation";
    EXPECT_EQ(sys.observablesList["other"], 7) << "unrelated observables must not be touched";
    EXPECT_TRUE(sys.moleculeList[0].isDissociated) << "sanity: the dissociation really did happen";
}

// -----------------------------------------------------------------------------
// A missing observable label must only produce a warning, never a crash and
// never a spurious map entry change.
// -----------------------------------------------------------------------------
void cd_test_forced_dissociation_with_unknown_observable_label()
{
    std::cerr << "\n[TEST] cd_test_forced_dissociation_with_unknown_observable_label\n"
              << "  Source file : src/reactions/check_dissociation.cpp\n"
              << "  Function    : check_dissociation()\n"
              << "  Scenario    : backRxn.isObserved == true but the label is not in\n"
              << "                observablesList.\n"
              << "  Criteria    : dissociation still completes, observablesList keeps\n"
              << "