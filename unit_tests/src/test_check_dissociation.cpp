/*! \file test_check_dissociation.cpp
 *
 * ### Unit test for ../src/reactions/check_dissociation.cpp
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
 * The routine walks the bound-interface list (`bndlist`) of one molecule and, for
 * each bond, decides whether the bond should be broken this time step.  The
 * decision logic that is fully visible in the source (and therefore what we
 * assert on) is:
 *
 *   * a bond is only considered if the interface really is flagged `isBound`;
 *   * the *lower indexed* molecule of the pair drives the reaction
 *     (`molItr > pro2Index` -> skip), so each pair is only attempted once;
 *   * implicit-lipid partners and partners that already dissociated this step
 *     are skipped;
 *   * an irreversible bond (`interaction.conjBackRxn == -1`) is skipped;
 *   * probability is `1 - exp(-kb * dt * 1e-6)`, forced to 1.0 when
 *     `params.debugParams.forceDissoc` is set;
 *   * on success: `break_interaction()` is called, `nBoundPairs` is decremented,
 *     the bound species copy number is decremented, both free species copy
 *     numbers are incremented, every member molecule/complex of both products is
 *     marked `TrajStatus::propagated`, both molecules get `isDissociated = true`,
 *     an optionally coupled reaction is performed, and an observable of the back
 *     reaction is decremented.
 *
 * NOTE: the "coupled destruction" branch is intentionally *not* exercised here.
 * It destroys an entire complex and mutates SimulVolume sub-cell membership, and
 * a mis-specified fixture would reach an `exit(1)` inside the library which would
 * abort the whole gtest binary rather than fail a single case.  The coupled
 * *uniMolStateChange* branch is exercised, with the fixture deliberately built so
 * that the `stateChangeProIndex == -1 -> exit(1)` path can never be taken.
 */

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

// The GSL random number generator lives in gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Absolute interface-state ("species") indices used by the fixture below.
//
//   0 : A(a)                 -- free binding site      (back-reaction product)
//   1 : A(a!1).A(a!1)        -- bound binding site      (back-reaction reactant)
//   2 : (unused padding)
//   3 : A(b~U)               -- free reporter iface, state U
//   4 : A(b~P)               -- free reporter iface, state P (coupled product)
// -----------------------------------------------------------------------------
constexpr int kCdFreeSite = 0;
constexpr int kCdBoundSite = 1;
constexpr int kCdStateU = 3;
constexpr int kCdStateP = 4;
constexpr int kCdNumSpecies = 5;

/*! \brief Everything check_dissociation() needs, kept together for convenience. */
struct CdFixture {
    Parameters params {};
    SimulVolume simulVolume {};
    std::vector<MolTemplate> molTemplateList {};
    std::map<std::string, int> observablesList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<BackRxn> backRxns {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {};
    copyCounters counterArrays {};
    Membrane membraneObject {};

    // break_interaction() logs association/dissociation events to this stream,
    // so hand it a real (temporary) file rather than an unopened stream.
    std::ofstream assocDissocFile {};
    std::string assocDissocFileName { "test_check_dissociation_events.tmp" };

    ~CdFixture()
    {
        if (assocDissocFile.is_open())
            assocDissocFile.close();
        std::remove(assocDissocFileName.c_str());
    }
};

/*! \brief Initialise the global GSL generator with a fixed seed (reproducible). */
void cd_init_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief Put all shared static bookkeeping into a state consistent with the
 *         two-molecule / one-complex fixture built below.
 *
 * These statics are shared by the whole test binary, so every test re-initialises
 * them instead of assuming a particular starting value.
 */
void cd_reset_statics()
{
    MolTemplate::numMolTypes = 1;
    MolTemplate::numEachMolType = std::vector<int> { 2 };

    Molecule::numberOfMolecules = 2;
    Molecule::emptyMolList.clear();
    Molecule::maxID = 2;

    Complex::numberOfComplexes = 1;
    Complex::obs.clear();
    Complex::maxID = 2;
    // Complex index 1 is a spare, already-empty complex which the library may
    // recycle when the dissociation splits the pair into two complexes.
    Complex::emptyComList.clear();
    Complex::emptyComList.push_back(1);
}

/*! \brief Build one molecule of type 0 that is bound to \p partnerIndex.
 *
 * Interface 0 is the binding site (bound, abs. index kCdBoundSite).
 * Interface 1 is a free "reporter" interface carrying states U/P; it exists so
 * that the coupled uniMolStateChange test always finds a matching interface.
 */
Molecule cd_make_molecule(int index, int partnerIndex, const Coord& com)
{
    Molecule mol {};
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = 0;
    mol.myComIndex = 0;
    mol.complexId = 0;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isDissociated = false;
    mol.trajStatus = TrajStatus::none;
    mol.mySubVolIndex = 0;
    mol.linksToSurface = 0;

    // --- interface 0: the bond that check_dissociation may break ---
    Molecule::Iface bindSite {};
    bindSite.coord = Coord(0.0, 0.0, 0.0); // both partners meet at the origin
    bindSite.index = kCdBoundSite;
    bindSite.relIndex = 0;
    bindSite.molTypeIndex = 0;
    bindSite.stateIndex = 0;
    bindSite.stateIden = '\0';
    bindSite.isBound = true;
    bindSite.interaction.partnerIndex = partnerIndex;
    bindSite.interaction.partnerIfaceIndex = 0;
    bindSite.interaction.partnerId = partnerIndex;
    bindSite.interaction.conjBackRxn = 0; // -> backRxns[0]

    // --- interface 1: free, has states, never touched by the dissociation ---
    Molecule::Iface reporter {};
    reporter.coord = com;
    reporter.index = kCdStateU;
    reporter.relIndex = 1;
    reporter.molTypeIndex = 0;
    reporter.stateIndex = 0;
    reporter.stateIden = 'U';
    reporter.isBound = false;

    mol.interfaceList = { bindSite, reporter };

    // legacy parallel bookkeeping lists (must stay mutually consistent)
    mol.bndlist = { 0 };
    mol.freelist = { 1 };
    mol.bndpartner = { partnerIndex };
    mol.bndRxnList = { 0 };

    return mol;
}

/*! \brief Build the complete fixture: 2 molecules of one type bound in 1 complex,
 *         one reversible bimolecular forward/back reaction pair with off rate
 *         \p offRate (units of s^-1).
 */
void cd_build_fixture(CdFixture& f, double offRate)
{
    cd_init_rng();
    cd_reset_statics();

    // ---------------- Parameters ----------------
    f.params.numMolTypes = 1;
    f.params.numTotalSpecies = kCdNumSpecies;
    f.params.timeStep = 1.0; // microseconds
    f.params.nItr = 10;
    f.params.numTotalComplex = 1;
    f.params.debugParams.forceDissoc = false;
    Parameters::dt = f.params.timeStep;

    // ---------------- Membrane (plain reflecting box, no implicit lipid) -----
    f.membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    f.membraneObject.isBox = true;
    f.membraneObject.isSphere = false;
    f.membraneObject.implicitLipid = false;
    f.membraneObject.implicitlipidIndex = -1;

    // ---------------- MolTemplate ----------------
    MolTemplate temp {};
    temp.molName = "A";
    temp.molTypeIndex = 0;
    temp.copies = 2;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord(10.0, 10.0, 10.0);
    temp.Dr = Coord(0.1, 0.1, 0.1);
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPoint = false;
    temp.isRod = false;
    temp.checkOverlap = false;
    temp.countTransition = false; // keeps the transition-matrix code out of play
    temp.canDestroy = false;

    Interface bindIface {};
    bindIface.name = "a";
    bindIface.index = 0;
    bindIface.iCoord = Coord(0.5, 0.0, 0.0);
    bindIface.stateList.emplace_back('\0', kCdFreeSite);

    Interface reportIface {};
    reportIface.name = "b";
    reportIface.index = 1;
    reportIface.iCoord = Coord(-0.5, 0.0, 0.0);
    reportIface.stateList.emplace_back('U', kCdStateU);
    reportIface.stateList.emplace_back('P', kCdStateP);

    temp.interfaceList = { bindIface, reportIface };
    f.molTemplateList.push_back(temp);

    // ---------------- Molecules ----------------
    f.moleculeList.push_back(cd_make_molecule(0, 1, Coord(-0.5, 0.0, 0.0)));
    f.moleculeList.push_back(cd_make_molecule(1, 0, Coord(0.5, 0.0, 0.0)));

    // ---------------- Complexes ----------------
    Complex bound {};
    bound.index = 0;
    bound.id = 0;
    bound.comCoord = Coord(0.0, 0.0, 0.0);
    bound.memberList = { 0, 1 };
    bound.numEachMol = std::vector<int>(MolTemplate::numMolTypes, 0);
    bound.numEachMol[0] = 2;
    bound.lastNumberUpdateItrEachMol = std::vector<long long int>(MolTemplate::numMolTypes, 0);
    bound.D = Coord(5.0, 5.0, 5.0);
    bound.Dr = Coord(0.05, 0.05, 0.05);
    bound.mass = 2.0;
    bound.radius = 2.0;
    bound.isEmpty = false;
    bound.OnSurface = false;
    bound.onFiber = false;
    bound.linksToSurface = 0;
    bound.iLipidIndex = 0;
    bound.trajStatus = TrajStatus::none;
    f.complexList.push_back(bound);

    // Spare, already-empty complex that the library can recycle for the product.
    Complex spare {};
    spare.index = 1;
    spare.id = 1;
    spare.isEmpty = true;
    spare.numEachMol = std::vector<int>(MolTemplate::numMolTypes, 0);
    spare.lastNumberUpdateItrEachMol = std::vector<long long int>(MolTemplate::numMolTypes, 0);
    spare.trajStatus = TrajStatus::none;
    f.complexList.push_back(spare);

    // ---------------- SimulVolume (single sub-volume holding both molecules) --
    f.simulVolume.numSubCells.x = 1;
    f.simulVolume.numSubCells.y = 1;
    f.simulVolume.numSubCells.z = 1;
    f.simulVolume.numSubCells.tot = 1;
    f.simulVolume.subCellSize = Coord(100.0, 100.0, 100.0);
    SimulVolume::SubVolume sub {};
    sub.absIndex = 0;
    sub.xIndex = 0;
    sub.yIndex = 0;
    sub.zIndex = 0;
    sub.memberMolList = { 0, 1 };
    f.simulVolume.subCellList.push_back(sub);

    // ---------------- copyCounters ----------------
    f.counterArrays.copyNumSpecies = std::vector<int>(kCdNumSpecies, 0);
    f.counterArrays.copyNumSpecies[kCdBoundSite] = 2; // two bound sites
    f.counterArrays.copyNumSpecies[kCdStateU] = 2;    // two reporter ifaces in U
    f.counterArrays.nBoundPairs = std::vector<int>(4, 0);
    f.counterArrays.nBoundPairs[0] = 1; // the single A-A bond
    f.counterArrays.proPairlist = std::vector<int>(4, 0);
    f.counterArrays.singleDouble = std::vector<int>(kCdNumSpecies, 0);
    f.counterArrays.implicitDouble = std::vector<bool>(kCdNumSpecies, false);
    f.counterArrays.canDissociate = std::vector<bool>(kCdNumSpecies, false);
    f.counterArrays.bindPairList.resize(kCdNumSpecies);
    f.counterArrays.bindPairListIL2D.resize(kCdNumSpecies);
    f.counterArrays.bindPairListIL3D.resize(kCdNumSpecies);
    f.counterArrays.nLoops = 0;

    // ---------------- Reactions ----------------
    // Two empty "other interface" lists: satisfies every implementation of
    // hasIntangibles() (no ancillary interfaces are required by the rate).
    std::vector<std::vector<RxnIface>> noOtherIfaces(2);

    RxnIface freeSite("a", 0, kCdFreeSite, 0, '\0', false);
    RxnIface boundSite("a", 0, kCdBoundSite, 0, '\0', true);

    ForwardRxn fwd {};
    fwd.rxnType = ReactionType::bimolecular;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.isSymmetric = true;
    fwd.isOnMem = false;
    fwd.bindRadius = 1.0;
    fwd.bindRadSameCom = 1.1;
    fwd.length3Dto2D = 2.0;
    fwd.irrevRingClosure = false;
    fwd.productName = "A(a!1).A(a!1)";
    fwd.reactantListNew = { freeSite, freeSite };
    fwd.productListNew = { boundSite, boundSite };
    fwd.intReactantList = { kCdFreeSite, kCdFreeSite };
    fwd.intProductList = { kCdBoundSite, kCdBoundSite };
    fwd.rateList.emplace_back(1.0, noOtherIfaces);
    f.forwardRxns.push_back(fwd);

    BackRxn back {};
    back.rxnType = ReactionType::bimolecular;
    back.absRxnIndex = 1;
    back.relRxnIndex = 0;
    back.conjForwardRxnIndex = 0;
    back.isSymmetric = true;
    back.isOnMem = false;
    back.length3Dto2D = 2.0;
    back.reactantListNew = { boundSite, boundSite }; // bound species
    back.productListNew = { freeSite, freeSite };    // free species
    back.intReactantList = { kCdBoundSite, kCdBoundSite };
    back.intProductList = { kCdFreeSite, kCdFreeSite };
    back.rateList.emplace_back(offRate, noOtherIfaces);
    back.isCoupled = false;
    back.isObserved = true;
    back.observeLabel = "AAdimer";
    f.backRxns.push_back(back);

    f.observablesList["AAdimer"] = 1;

    // ---------------- output stream used by break_interaction() -------------
    f.assocDissocFile.open(f.assocDissocFileName.c_str());
}

/*! \brief Attach a coupled uniMolStateChange reaction (A(b~U) -> A(b~P)) to the
 *         back reaction, as forwardRxns[1].
 *
 * Both molecules own an interface whose absolute index is kCdStateU, so
 * check_dissociation() is guaranteed to find a state-change target and can never
 * hit its `exit(1)` "products do not match" path.
 */
void cd_attach_coupled_state_change(CdFixture& f)
{
    std::vector<std::vector<RxnIface>> noOtherIfaces(1);

    RxnIface stateU("b", 0, kCdStateU, 1, 'U', false);
    RxnIface stateP("b", 0, kCdStateP, 1, 'P', false);

    ForwardRxn stateChange {};
    stateChange.rxnType = ReactionType::uniMolStateChange;
    stateChange.absRxnIndex = 2;
    stateChange.relRxnIndex = 1;
    stateChange.isReversible = false;
    stateChange.conjBackRxnIndex = -1;
    stateChange.hasStateChange = true;
    stateChange.reactantListNew = { stateU };
    stateChange.productListNew = { stateP };
    stateChange.intReactantList = { kCdStateU };
    stateChange.intProductList = { kCdStateP };
    stateChange.stateChangeIface = std::pair<RxnIface, RxnIface> { stateU, stateP };
    stateChange.rateList.emplace_back(1.0, noOtherIfaces);
    stateChange.isObserved = true;
    stateChange.observeLabel = "Bphospho";
    f.forwardRxns.push_back(stateChange); // index 1

    f.backRxns[0].isCoupled = true;
    f.backRxns[0].coupledRxn.rxnType = ReactionType::uniMolStateChange;
    f.backRxns[0].coupledRxn.relRxnIndex = 1; // -> forwardRxns[1]
    f.backRxns[0].coupledRxn.absRxnIndex = 2;
    f.backRxns[0].coupledRxn.probCoupled = 1.0; // always perform it

    f.observablesList["Bphospho"] = 0;
}

/*! \brief Convenience wrapper around the (long) call signature. */
void cd_call(CdFixture& f, unsigned int molItr, unsigned int simItr = 1)
{
    check_dissociation(simItr, f.params, f.simulVolume, f.molTemplateList, f.observablesList, molItr,
        f.moleculeList, f.complexList, f.backRxns, f.forwardRxns, f.createDestructRxns, f.counterArrays,
        f.membraneObject, f.assocDissocFile);
}

/*! \brief Assert that nothing about the bond changed. Used by all "skip" cases. */
void cd_expect_no_dissociation(const CdFixture& f, const char* why)
{
    std::cerr << "    checking that no dissociation happened (" << why << ")\n";
    EXPECT_FALSE(f.moleculeList[0].isDissociated) << "molecule 0 should not be flagged dissociated";
    EXPECT_FALSE(f.moleculeList[1].isDissociated) << "molecule 1 should not be flagged dissociated";
    EXPECT_EQ(f.counterArrays.copyNumSpecies[kCdBoundSite], 2)
        << "bound-species copy number must be untouched";
    EXPECT_EQ(f.counterArrays.copyNumSpecies[kCdFreeSite], 0)
        << "free-species copy number must be untouched";
    EXPECT_EQ(f.counterArrays.nBoundPairs[0], 1) << "nBoundPairs must be untouched";
    EXPECT_EQ(f.observablesList.at("AAdimer"), 1) << "observable must be untouched";
    EXPECT_EQ(f.moleculeList[0].myComIndex, f.moleculeList[1].myComIndex)
        << "the two molecules must still share one complex";
    EXPECT_EQ(static_cast<int>(f.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must not be promoted to propagated";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: an interface listed in bndlist but flagged unbound must be ignored.
// -----------------------------------------------------------------------------
void test_check_dissoc_unbound_interface_ignored()
{
    std::cerr << "\n[TEST] test_check_dissoc_unbound_interface_ignored\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation()\n"
              << "  Scenario:      relIface 0 is in bndlist but Iface::isBound == false,\n"
              << "                 and forceDissoc is on so any accepted bond would break.\n"
              << "  Pass criteria: the `isBound` guard rejects the bond; no state changes.\n";

    CdFixture f;
    cd_build_fixture(f, 1.0);
    f.params.debugParams.forceDissoc = true;      // would dissociate if considered
    f.moleculeList[0].interfaceList[0].isBound = false; // ...but it is not bound

    cd_call(f, 0);
    cd_expect_no_dissociation(f, "interface not flagged bound");
}

// -----------------------------------------------------------------------------
// Test 2: an empty bndlist means the loop body never executes.
// -----------------------------------------------------------------------------
void test_check_dissoc_empty_bndlist()
{
    std::cerr << "\n[TEST] test_check_dissoc_empty_bndlist\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation()\n"
              << "  Scenario:      the molecule reports no bound interfaces at all.\n"
              << "  Pass criteria: the for-loop never runs, nothing is modified.\n";

    CdFixture f;
    cd_build_fixture(f, 1.0);
    f.params.debugParams.forceDissoc = true;
    f.moleculeList[0].bndlist.clear();

    cd_call(f, 0);
    cd_expect_no_dissociation(f, "bndlist is empty");
}

// -----------------------------------------------------------------------------
// Test 3: irreversible bond (conjBackRxn == -1) is skipped.
// -----------------------------------------------------------------------------
void test_check_dissoc_irreversible_bond_skipped()
{
    std::cerr << "\n[TEST] test_check_dissoc_irreversible_bond_skipped\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation()\n"
              << "  Scenario:      interaction.conjBackRxn == -1 (no back reaction exists),\n"
              << "                 forceDissoc is on.\n"
              << "  Pass criteria: `mu == -1` triggers `continue`; nothing is modified.\n";

    CdFixture f;
    cd_build_fixture(f, 1.0);
    f.params.debugParams.forceDissoc = true;
    f.moleculeList[0].interfaceList[0].interaction.conjBackRxn = -1;

    cd_call(f, 0);
    cd_expect_no_dissociation(f, "reaction is irreversible");
}

// -----------------------------------------------------------------------------
// Test 4: only the lower-indexed molecule of the pair drives dissociation.
// -----------------------------------------------------------------------------
void test_check_dissoc_higher_index_molecule_skipped()
{
    std::cerr << "\n[TEST] test_check_dissoc_higher_index_molecule_skipped\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation()\n"
              << "  Scenario:      called for molecule 1 whose partner is molecule 0,\n"
              << "                 so `molItr > pro2Index`; forceDissoc is on.\n"
              << "  Pass criteria: the pair is skipped so that a bond is only ever\n"
              << "                 attempted once per time step.\n";

    CdFixture f;
    cd_build_fixture(f, 1.0);
    f.params.debugParams.forceDissoc = true;

    cd_call(f, 1); // molItr = 1 > pro2Index = 0
    cd_expect_no_dissociation(f, "molItr > pro2Index");
}

// -----------------------------------------------------------------------------
// Test 5: implicit-lipid partners are handled elsewhere and skipped here.
// -----------------------------------------------------------------------------
void test_check_dissoc_implicit_lipid_partner_skipped()
{
    std::cerr << "\n[TEST] test_check_dissoc_implicit_lipid_partner_skipped\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation()\n"
              << "  Scenario:      the partner molecule is an implicit lipid;\n"
              << "                 forceDissoc is on.\n"
              << "  Pass criteria: skipped (implicit lipids use\n"
              << "                 check_dissociation_implicitlipid instead).\n";

    CdFixture f;
    cd_build_fixture(f, 1.0);
    f.params.debugParams.forceDissoc = true;
    f.moleculeList[1].isImplicitLipid = true;

    cd_call(f, 0);
    cd_expect_no_dissociation(f, "partner is an implicit lipid");
}

// -----------------------------------------------------------------------------
// Test 6: a partner that already dissociated this step is skipped.
// -----------------------------------------------------------------------------
void test_check_dissoc_partner_already_dissociated_skipped()
{
    std::cerr << "\n[TEST] test_check_dissoc_partner_already_dissociated_skipped\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation()\n"
              << "  Scenario:      the partner already carries isDissociated == true;\n"
              << "                 forceDissoc is on.\n"
              << "  Pass criteria: skipped so a molecule dissociates at most once per step.\n";

    CdFixture f;
    cd_build_fixture(f, 1.0);
    f.params.debugParams.forceDissoc = true;
    f.moleculeList[1].isDissociated = true; // pretend it dissociated earlier

    cd_call(f, 0);

    std::cerr << "    checking that molecule 0 was left alone\n";
    EXPECT_FALSE(f.moleculeList[0].isDissociated) << "molecule 0 should not be flagged dissociated";
    EXPECT_EQ(f.counterArrays.copyNumSpecies[kCdBoundSite], 2) << "bound copy number untouched";
    EXPECT_EQ(f.counterArrays.copyNumSpecies[kCdFreeSite], 0) << "free copy number untouched";
    EXPECT_EQ(f.counterArrays.nBoundPairs[0], 1) << "nBoundPairs untouched";
    EXPECT_EQ(f.observablesList.at("AAdimer"), 1) << "observable untouched";
    EXPECT_EQ(f.moleculeList[0].myComIndex, f.moleculeList[1].myComIndex)
        << "molecules must still share a complex";
}

// -----------------------------------------------------------------------------
// Test 7: a zero off-rate gives probability exactly 0 -> never dissociates.
// -----------------------------------------------------------------------------
void test_check_dissoc_zero_rate_never_dissociates()
{
    std::cerr << "\n[TEST] test_check_dissoc_zero_rate_never_dissociates\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation()\n"
              << "  Scenario:      off rate kb = 0 so prob = 1 - exp(0) = 0 exactly,\n"
              << "                 forceDissoc is off.\n"
              << "  Pass criteria: `prob > rand_gsl()` is always false; nothing changes.\n";

    CdFixture f;
    cd_build_fixture(f, 0.0); // kb = 0 s^-1
    f.params.debugParams.forceDissoc = false;

    cd_call(f, 0);
    cd_expect_no_dissociation(f, "dissociation probability is exactly zero");
}

// -----------------------------------------------------------------------------
// Test 8: the full accepted-dissociation bookkeeping (forceDissoc == true).
// -----------------------------------------------------------------------------
void test_check_dissoc_forced_dissociation_bookkeeping()
{
    std::cerr << "\n[TEST] test_check_dissoc_forced_dissociation_bookkeeping\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation()\n"
              << "  Scenario:      debugParams.forceDissoc = true forces prob = 1.0 for the\n"
              << "                 single A(a!1).A(a!1) bond between molecules 0 and 1.\n"
              << "  Pass criteria: bound species -1, both free species +1, nBoundPairs -1,\n"
              << "                 both molecules flagged dissociated & propagated, both\n"
              << "                 complexes flagged propagated, interfaces unbound, the\n"
              << "                 pair no longer shares a complex, observable decremented.\n";

    CdFixture f;
    cd_build_fixture(f, 1.0);
    f.params.debugParams.forceDissoc = true;

    const int boundBefore = f.counterArrays.copyNumSpecies[kCdBoundSite];
    const int freeBefore = f.counterArrays.copyNumSpecies[kCdFreeSite];
    const int pairsBefore = f.counterArrays.nBoundPairs[0];

    std::cerr << "    calling check_dissociation(simItr=1, molItr=0)...\n";
    cd_call(f, 0);

    // --- species copy numbers -------------------------------------------------
    std::cerr << "    copyNumSpecies[bound] " << boundBefore << " -> "
              << f.counterArrays.copyNumSpecies[kCdBoundSite] << ", copyNumSpecies[free] "
              << freeBefore << " -> " << f.counterArrays.copyNumSpecies[kCdFreeSite] << '\n';
    EXPECT_EQ(f.counterArrays.copyNumSpecies[kCdBoundSite], boundBefore - 1)
        << "the bound species copy number must be decremented once";
    EXPECT_EQ(f.counterArrays.copyNumSpecies[kCdFreeSite], freeBefore + 2)
        << "both free product interfaces must be incremented";

    // --- bound pair counter ---------------------------------------------------
    std::cerr << "    nBoundPairs[0] " << pairsBefore << " -> " << f.counterArrays.nBoundPairs[0] << '\n';
    EXPECT_EQ(f.counterArrays.nBoundPairs[0], pairsBefore - 1)
        << "update_Nboundpairs(...,-1,...) must remove the A-A pair";

    // --- dissociation flags ---------------------------------------------------
    EXPECT_TRUE(f.moleculeList[0].isDissociated) << "molecule 0 must be flagged dissociated";
    EXPECT_TRUE(f.moleculeList[1].isDissociated) << "molecule 1 must be flagged dissociated";

    // --- trajectory status ----------------------------------------------------
    EXPECT_EQ(static_cast<int>(f.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "molecule 0 must be marked propagated so it is not moved again this step";
    EXPECT_EQ(static_cast<int>(f.moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "molecule 1 must be marked propagated so it is not moved again this step";

    const int com0 = f.moleculeList[0].myComIndex;
    const int com1 = f.moleculeList[1].myComIndex;
    std::cerr << "    product complex indices: mol0 -> " << com0 << ", mol1 -> " << com1 << '\n';
    ASSERT_GE(com0, 0);
    ASSERT_LT(com0, static_cast<int>(f.complexList.size()));
    ASSERT_GE(com1, 0);
    ASSERT_LT(com1, static_cast<int>(f.complexList.size()));
    EXPECT_NE(com0, com1) << "breaking the only bond must split the complex in two";
    EXPECT_EQ(static_cast<int>(f.complexList[com0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "product complex of molecule 0 must be marked propagated";
    EXPECT_EQ(static_cast<int>(f.complexList[com1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "product complex of molecule 1 must be marked propagated";

    // --- the bond itself really is gone --------------------------------------
    EXPECT_FALSE(f.moleculeList[0].interfaceList[0].isBound)
        << "molecule 0's binding site must be free after break_interaction()";
    EXPECT_FALSE(f.moleculeList[1].interfaceList[0].isBound)
        << "molecule 1's binding site must be free after break_interaction()";

    // --- observable of the back reaction ------------------------------------
    EXPECT_EQ(f.observablesList.at("AAdimer"), 0)
        << "the observed back-reaction product count must be decremented";
}

// -----------------------------------------------------------------------------
// Test 9: dissociation coupled to a unimolecular state change.
// -----------------------------------------------------------------------------
void test_check_dissoc_coupled_state_change()
{
    std::cerr << "\n[TEST] test_check_dissoc_coupled_state_change\n"
              << "  Source file:   src/reactions/check_dissociation.cpp\n"
              << "  Function:      check_dissociation() -- isCoupled / uniMolStateChange branch\n"
              << "  Scenario:      the back reaction is coupled (probCoupled = 1.0) to\n"
              << "                 A(b~U) -> A(b~P); forceDissoc forces the dissociation.\n"
              << "                 Both products own an interface with absolute index "
              << kCdStateU << ",\n"
              << "                 and the source keeps the *last* match (the partner),\n"
              << "                 so molecule 1 is the state-change target.\n"
              << "  Pass criteria: molecule 1's reporter interface becomes state 'P'\n"
              << "                 (abs. index " << kCdStateP << ", rel. state index 1),\n"
              << "                 copyNumSpecies U -1 / P +1, observable incremented,\n"
              << "                 molecule 0's reporter interface untouched.\n";

    CdFixture f;
    cd_build_fixture(f, 1.0);
    cd_attach_coupled_state_change(f);
    f.params.debugParams.forceDissoc = true;

    const int uBefore = f.counterArrays.copyNumSpecies[kCdStateU];
    const int pBefore = f.counterArrays.copyNumSpecies[kCdStateP];

    std::cerr << "    calling check_dissociation(simItr=1, molItr=0)...\n";
    cd_call(f, 0);

    // The dissociation itself must still have happened.
    EXPECT_TRUE(f.moleculeList[0].isDissociated) << "molecule 0 must be flagged dissociated";
    EXPECT_TRUE(f.moleculeList[1].isDissociated) << "molecule 1 must be flagged dissociated";

    // Coupled state change applied to the *partner* (last match wins in the source).
    std::cerr << "    mol1 reporter iface: absIndex=" << f.moleculeList[1].interfaceList[1].index
              << ", stateIden='" << f.moleculeList[1].interfaceList[1].stateIden
              << "', stateIndex=" << f.moleculeList[1].interfaceList[1].stateIndex << '\n';
    EXPECT_EQ(f.moleculeList[1].interfaceList[1].index, kCdStateP)
        << "partner's reporter interface must take the coupled product's absolute index";
    EXPECT_EQ(f.moleculeList[1].interfaceList[1].stateIden, 'P')
        << "partner's reporter interface must carry the new state identity";
    EXPECT_EQ(f.moleculeList[1].interfaceList[1].stateIndex, 1)
        << "relative state index must point at stateList[1] (the 'P' state)";

    // Molecule 0 must be untouched by the state change.
    std::cerr << "    mol0 reporter iface: absIndex=" << f.moleculeList[0].interfaceList[1].index
              << ", stateIden='" << f.moleculeList[0].interfaceList[1].stateIden << "'\n";
    EXPECT_EQ(f.moleculeList[0].interfaceList[1].index, kCdStateU)
        << "only one molecule may undergo the coupled state change";
    EXPECT_EQ(f.moleculeList[0].interfaceList[1].stateIden, 'U')
        << "molecule 0's reporter state identity must be unchanged";

    // Copy numbers for the coupled reaction.
    std::cerr << "    copyNumSpecies[U] " << uBefore << " -> "
              << f.counterArrays.copyNumSpecies[kCdStateU] << ", copyNumSpecies[P] " << pBefore
              << " -> " << f.counterArrays.copyNumSpecies[kCdStateP] << '\n';
    EXPECT_EQ(f.counterArrays.copyNumSpecies[kCdStateU], uBefore - 1)
        << "coupled reactant species must be decremented exactly once";
    EXPECT_EQ(f.counterArrays.copyNumSpecies[kCdStateP], pBefore + 1)
        << "coupled product species must be incremented exactly once";

    // Observables of both the back reaction and the coupled reaction.
    EXPECT_EQ(f.observablesList.at("Bphospho"), 1)
        << "coupled reaction observable must be incremented";
    EXPECT_EQ(f.observablesList.at("AAdimer"), 0)
        << "back-reaction observable must be decremented";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- one per scenario so every case is reported separately
// and a failure in one does not stop the others from running.
// -----------------------------------------------------------------------------
TEST(CheckDissociation, UnboundInterfaceIgnored) { test_check_dissoc_unbound_interface_ignored(); }
TEST(CheckDissociation, EmptyBndList) { test_check_dissoc_empty_bndlist(); }
TEST(CheckDissociation, IrreversibleBondSkipped) { test_check_dissoc_irreversible_bond_skipped(); }
TEST(CheckDissociation, HigherIndexMoleculeSkipped) { test_check_dissoc_higher_index_molecule_skipped(); }
TEST(CheckDissociation, ImplicitLipidPartnerSkipped) { test_check_dissoc_implicit_lipid_partner_skipped(); }
TEST(CheckDissociation, PartnerAlreadyDissociatedSkipped) { test_check_dissoc_partner_already_dissociated_skipped(); }
TEST(CheckDissociation, ZeroRateNeverDissociates) { test_check_dissoc_zero_rate_never_dissociates(); }
TEST(CheckDissociation, ForcedDissociationBookkeeping) { test_check_dissoc_forced_dissociation_bookkeeping(); }
TEST(CheckDissociation, CoupledStateChange) { test_check_dissoc_coupled_state_change(); }