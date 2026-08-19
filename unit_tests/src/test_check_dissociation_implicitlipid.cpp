/*! \file test_check_dissociation_implicitlipid.cpp
 *
 * ### Unit test for src/reactions/check_dissociation_implicitlipid.cpp
 *
 * Function under test:
 * \code
 * void check_dissociation_implicitlipid(long long int simItr, const Parameters&,
 *          SimulVolume&, std::vector<MolTemplate>&, std::map<std::string,int>&,
 *          unsigned int molItr, std::vector<Molecule>&, std::vector<Complex>&,
 *          const std::vector<BackRxn>&, const std::vector<ForwardRxn>&,
 *          const std::vector<CreateDestructRxn>&, copyCounters&, Membrane&,
 *          std::vector<double>&, std::vector<double>&, std::vector<double>&,
 *          std::ofstream&)
 * \endcode
 *
 * The routine walks the bound-interface list (`bndlist`) of one molecule and, for
 * every bond that connects it to the *implicit lipid*, evaluates a dissociation
 * probability and (if the reaction fires) breaks the bond, updates all copy
 * number bookkeeping, releases the complex from the membrane and repositions it
 * at the 3D reflecting surface.
 *
 * The tests below build a minimal but *fully initialised* system consisting of
 *   - molTemplate 0 : "prot" (a soluble protein with one interface "a")
 *   - molTemplate 1 : "il"   (the implicit lipid, one interface "il")
 *   - one bound pair prot(a!1).il(il!1) sitting on the bottom membrane
 * and then exercise:
 *   1. the four early-exit / skip conditions (ghosted molecule, complex not on
 *      the surface, partner is not an implicit lipid, irreversible bond,
 *      interface flagged unbound),
 *   2. the "probability is far too small to fire" path (state must be untouched),
 *   3. the forced-dissociation path (`debugParams.forceDissoc == true`), where
 *      every documented side effect is verified.
 *
 * NOTE: the 2D->2D branch (`linksToSurface != 1` or `membraneObject.TwoD`) calls
 * dissociate2D(), which performs GSL numerical integration.  A GSL integration
 * error aborts the whole process, which would take down the entire test binary,
 * so that branch is deliberately not exercised here.
 */

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gtest/gtest.h>
#include <gsl/gsl_rng.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Constants describing the miniature reaction network used by every test.
// -----------------------------------------------------------------------------
constexpr int kCdilProtFreeState = 0; //!< absolute index: prot(a) free
constexpr int kCdilIlFreeState = 1; //!< absolute index: il(il) free
constexpr int kCdilBoundState = 2; //!< absolute index: prot(a!1).il(il!1)

constexpr double kCdilBindRadius = 1.0; //!< sigma [nm]
constexpr double kCdilOnRate = 10.0; //!< ka [nm^3 us^-1]
constexpr double kCdilOffRate = 1.0; //!< kb [us^-1]
constexpr double kCdilRS3D = 0.5; //!< 3D reflecting surface offset [nm]
constexpr double kCdilBoxLen = 100.0; //!< cubic water box side length [nm]

// -----------------------------------------------------------------------------
// Initialise (or re-seed) the global GSL random number generator.  The generator
// is only consulted at the very end of the routine (prob > rand_gsl()), so a
// fixed seed keeps the tests reproducible.
// -----------------------------------------------------------------------------
void cdil_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    } else {
        gsl_rng_set(r, 42);
    }
}

// -----------------------------------------------------------------------------
// Container holding every argument required by check_dissociation_implicitlipid.
// -----------------------------------------------------------------------------
struct CdilSystem {
    Parameters params {};
    SimulVolume simulVolume {};
    Membrane membraneObject {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
    std::vector<double> IL2DbindingVec {};
    std::vector<double> IL2DUnbindingVec {};
    std::vector<double> ILTableIDs {};
    std::ofstream assocDissocFile {};
};

// -----------------------------------------------------------------------------
// Build the complete, self-consistent test system:
//   moleculeList[0] : protein, complexList[0], one interface bound to the IL
//   moleculeList[1] : implicit lipid, complexList[1]
// -----------------------------------------------------------------------------
void cdil_build(CdilSystem& s)
{
    /* ---------------- simulation parameters ---------------- */
    s.params.numMolTypes = 2;
    s.params.numTotalSpecies = 3;
    s.params.timeStep = 0.1; // us
    s.params.debugParams.forceDissoc = false;
    s.params.rank = 0;

    /* ---------------- molecule templates ---------------- */
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };

    MolTemplate protTemp;
    protTemp.molName = "prot";
    protTemp.molTypeIndex = 0;
    protTemp.copies = 1;
    protTemp.mass = 1.0;
    protTemp.radius = 1.0;
    protTemp.D = Coord(1.0, 1.0, 1.0); // freely diffusing in 3D
    protTemp.Dr = Coord(0.1, 0.1, 0.1);
    protTemp.isLipid = false;
    protTemp.isImplicitLipid = false;
    protTemp.canDestroy = false;
    {
        Interface iface;
        iface.index = 0;
        iface.name = "a";
        iface.iCoord = Coord(0.0, 0.0, -1.0);
        iface.stateList.emplace_back(Interface::State("a", '\0', kCdilProtFreeState));
        protTemp.interfaceList.push_back(iface);
    }
    s.molTemplateList.push_back(protTemp);

    MolTemplate ilTemp;
    ilTemp.molName = "il";
    ilTemp.molTypeIndex = 1;
    ilTemp.copies = 1;
    ilTemp.mass = 1.0;
    ilTemp.radius = 1.0;
    ilTemp.D = Coord(1.0, 1.0, 0.0); // confined to the membrane
    ilTemp.Dr = Coord(0.0, 0.0, 0.0);
    ilTemp.isLipid = true;
    ilTemp.isImplicitLipid = true;
    ilTemp.canDestroy = false;
    {
        Interface iface;
        iface.index = 0;
        iface.name = "il";
        iface.iCoord = Coord(0.0, 0.0, 0.0);
        iface.stateList.emplace_back(Interface::State("il", '\0', kCdilIlFreeState));
        ilTemp.interfaceList.push_back(iface);
    }
    s.molTemplateList.push_back(ilTemp);

    /* ---------------- forward (binding) reaction ---------------- */
    // prot(a) + il(il) <-> prot(a!1).il(il!1)
    ForwardRxn fwd;
    fwd.rxnType = ReactionType::bimolecular;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.bindRadius = kCdilBindRadius;
    fwd.length3Dto2D = 2.0 * kCdilBindRadius;
    fwd.isOnMem = false;
    fwd.isCoupled = false;
    fwd.isObserved = false;
    fwd.reactantListNew.emplace_back("a", 0, kCdilProtFreeState, 0, '\0', false);
    fwd.reactantListNew.emplace_back("il", 1, kCdilIlFreeState, 0, '\0', false);
    fwd.productListNew.emplace_back("a", 0, kCdilBoundState, 0, '\0', true);
    fwd.productListNew.emplace_back("il", 1, kCdilBoundState, 0, '\0', true);
    fwd.intReactantList = std::vector<int> { kCdilProtFreeState, kCdilIlFreeState };
    fwd.intProductList = std::vector<int> { kCdilBoundState, kCdilBoundState };
    // Two *empty* ancillary interface lists: this makes hasIntangibles() return
    // true (no ancillary interfaces required) while still being safely indexable.
    fwd.rateList.emplace_back(kCdilOnRate, std::vector<std::vector<RxnIface>> { {}, {} });
    s.forwardRxns.push_back(fwd);

    /* ---------------- back (unbinding) reaction ---------------- */
    BackRxn back;
    back.rxnType = ReactionType::bimolecular;
    back.absRxnIndex = 1;
    back.relRxnIndex = 0;
    back.conjForwardRxnIndex = 0;
    back.isCoupled = false; // keep the coupled-reaction branches out of play
    back.isObserved = false;
    back.reactantListNew.emplace_back("a", 0, kCdilBoundState, 0, '\0', true);
    back.reactantListNew.emplace_back("il", 1, kCdilBoundState, 0, '\0', true);
    back.productListNew.emplace_back("a", 0, kCdilProtFreeState, 0, '\0', false);
    back.productListNew.emplace_back("il", 1, kCdilIlFreeState, 0, '\0', false);
    back.intReactantList = std::vector<int> { kCdilBoundState, kCdilBoundState };
    back.intProductList = std::vector<int> { kCdilProtFreeState, kCdilIlFreeState };
    back.rateList.emplace_back(kCdilOffRate, std::vector<std::vector<RxnIface>> { {}, {} });
    s.backRxns.push_back(back);

    /* ---------------- membrane / boundary ---------------- */
    s.membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { kCdilBoxLen, kCdilBoxLen, kCdilBoxLen });
    s.membraneObject.isBox = true;
    s.membraneObject.isSphere = false;
    s.membraneObject.TwoD = false; // -> the 2D->3D (safe) branch
    s.membraneObject.hasCompartment = false;
    s.membraneObject.implicitLipid = true;
    s.membraneObject.implicitlipidIndex = 1; // moleculeList index of the IL
    s.membraneObject.nStates = 1;
    s.membraneObject.nSites = 1;
    s.membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 100 };
    s.membraneObject.numberOfProteinEachState = std::vector<int> { 1 };
    s.membraneObject.No_free_lipids = 100;
    s.membraneObject.No_protein = 1;
    s.membraneObject.totalSA = kCdilBoxLen * kCdilBoxLen;

    // The RS3D lookup table: entries [i], [i+100], [i+200] must match the
    // binding radius, the forward rate and the reactants' summed diffusion
    // constant *exactly* (tolerance 1e-15) for [i+300] to be picked up as RS3D.
    const Coord& d0 = s.molTemplateList[fwd.reactantListNew[0].molTypeIndex].D;
    const Coord& d1 = s.molTemplateList[fwd.reactantListNew[1].molTypeIndex].D;
    const double dSum = 1.0 / 3.0 * (d0.x + d1.x) + 1.0 / 3.0 * (d0.y + d1.y) + 1.0 / 3.0 * (d0.z + d1.z);
    s.membraneObject.RS3Dvect.assign(500, 0.0);
    s.membraneObject.RS3Dvect[0] = kCdilBindRadius;
    s.membraneObject.RS3Dvect[100] = kCdilOnRate;
    s.membraneObject.RS3Dvect[200] = dSum;
    s.membraneObject.RS3Dvect[300] = kCdilRS3D;
    s.membraneObject.RS3Dvect[400] = 0.0; // molTypeIndex of the bound protein

    /* ---------------- molecules ---------------- */
    Molecule::numberOfMolecules = 2;
    Complex::numberOfComplexes = 2;

    Molecule prot;
    prot.index = 0;
    prot.id = 0;
    prot.molTypeIndex = 0;
    prot.myComIndex = 0;
    prot.mass = 1.0;
    prot.isLipid = false;
    prot.isImplicitLipid = false;
    prot.isGhosted = false;
    prot.isEmpty = false;
    prot.linksToSurface = 1; // exactly one link -> 2D->3D branch
    prot.trajStatus = TrajStatus::none;
    prot.comCoord = Coord(0.0, 0.0, -49.0);
    {
        Molecule::Iface pIface;
        pIface.coord = Coord(0.0, 0.0, -kCdilBoxLen / 2.0); // sitting on the membrane
        pIface.index = kCdilBoundState; // currently the bound species
        pIface.relIndex = 0;
        pIface.molTypeIndex = 0;
        pIface.stateIndex = 0;
        pIface.stateIden = '\0';
        pIface.isBound = true;
        pIface.interaction.partnerIndex = 1; // the implicit lipid
        pIface.interaction.partnerIfaceIndex = 0;
        pIface.interaction.conjBackRxn = 0; // reversible -> backRxns[0]
        prot.interfaceList.push_back(pIface);
    }
    prot.bndlist = std::vector<int> { 0 };
    prot.freelist = std::vector<int> {};
    prot.bndpartner = std::vector<int> { 1 };
    prot.bndRxnList = std::vector<int> { 0 };
    s.moleculeList.push_back(prot);

    Molecule il;
    il.index = 1;
    il.id = 1;
    il.molTypeIndex = 1;
    il.myComIndex = 1;
    il.mass = 1.0;
    il.isLipid = true;
    il.isImplicitLipid = true;
    il.isGhosted = false;
    il.isEmpty = false;
    il.linksToSurface = 0;
    il.trajStatus = TrajStatus::none;
    il.comCoord = Coord(0.0, 0.0, -kCdilBoxLen / 2.0);
    {
        Molecule::Iface iIface;
        iIface.coord = Coord(0.0, 0.0, -kCdilBoxLen / 2.0);
        iIface.index = kCdilBoundState;
        iIface.relIndex = 0;
        iIface.molTypeIndex = 1;
        iIface.stateIndex = 0;
        iIface.stateIden = '\0';
        iIface.isBound = true;
        iIface.interaction.partnerIndex = 0;
        iIface.interaction.partnerIfaceIndex = 0;
        iIface.interaction.conjBackRxn = 0;
        il.interfaceList.push_back(iIface);
    }
    // Give the IL both list entries so that a symmetric "break" implementation
    // always finds what it looks for.
    il.bndlist = std::vector<int> { 0 };
    il.freelist = std::vector<int> {};
    il.bndpartner = std::vector<int> { 0 };
    il.bndRxnList = std::vector<int> { 0 };
    s.moleculeList.push_back(il);

    /* ---------------- complexes ---------------- */
    Complex protCom;
    protCom.index = 0;
    protCom.id = 0;
    protCom.comCoord = s.moleculeList[0].comCoord;
    protCom.memberList = std::vector<int> { 0 };
    protCom.numEachMol = std::vector<int> { 1, 0 };
    protCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    protCom.mass = 1.0;
    protCom.radius = 1.0;
    protCom.D = Coord(0.5, 0.5, 0.0); // membrane-bound complex
    protCom.Dr = Coord(0.05, 0.05, 0.05);
    protCom.OnSurface = true; // required to even consider dissociation
    protCom.linksToSurface = 1;
    protCom.iLipidIndex = 1;
    protCom.isEmpty = false;
    protCom.trajStatus = TrajStatus::none;
    s.complexList.push_back(protCom);

    Complex ilCom;
    ilCom.index = 1;
    ilCom.id = 1;
    ilCom.comCoord = s.moleculeList[1].comCoord;
    ilCom.memberList = std::vector<int> { 1 };
    ilCom.numEachMol = std::vector<int> { 0, 1 };
    ilCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    ilCom.mass = 1.0;
    ilCom.radius = 1.0;
    ilCom.D = Coord(1.0, 1.0, 0.0);
    ilCom.Dr = Coord(0.0, 0.0, 0.0);
    ilCom.OnSurface = true;
    ilCom.linksToSurface = 0;
    ilCom.iLipidIndex = 1;
    ilCom.isEmpty = false;
    ilCom.trajStatus = TrajStatus::none;
    s.complexList.push_back(ilCom);

    /* ---------------- copy-number bookkeeping ---------------- */
    // Oversize the arrays so that any plausible index formula stays in bounds.
    s.counterArrays.copyNumSpecies.assign(16, 0);
    s.counterArrays.copyNumSpecies[kCdilProtFreeState] = 0;
    s.counterArrays.copyNumSpecies[kCdilIlFreeState] = 100;
    s.counterArrays.copyNumSpecies[kCdilBoundState] = 1;
    s.counterArrays.nBoundPairs.assign(16, 0);
    s.counterArrays.nBoundPairs[0 * s.params.numMolTypes + 1] = 1;
    s.counterArrays.nBoundPairs[1 * s.params.numMolTypes + 0] = 1;
    s.counterArrays.proPairlist.assign(16, 0);
    s.counterArrays.singleDouble.assign(16, 0);
    s.counterArrays.implicitDouble.assign(16, false);
    s.counterArrays.canDissociate.assign(16, false);

    /* ---------------- association/dissociation log ---------------- */
    s.assocDissocFile.open("test_cdil_assoc_dissoc.dat");
}

// -----------------------------------------------------------------------------
// Convenience wrapper: run the routine on molecule 0 at iteration 1.
// -----------------------------------------------------------------------------
void cdil_call(CdilSystem& s)
{
    check_dissociation_implicitlipid(1, s.params, s.simulVolume, s.molTemplateList, s.observablesList,
        /*molItr=*/0, s.moleculeList, s.complexList, s.backRxns, s.forwardRxns, s.createDestructRxns,
        s.counterArrays, s.membraneObject, s.IL2DbindingVec, s.IL2DUnbindingVec, s.ILTableIDs,
        s.assocDissocFile);
}

// -----------------------------------------------------------------------------
// A compact snapshot of everything the routine is allowed to touch.  Used by the
// "nothing must change" tests.
// -----------------------------------------------------------------------------
struct CdilSnapshot {
    bool isBound { false };
    int partnerIndex { -99 };
    int ifaceAbsIndex { -99 };
    int molLinksToSurface { -99 };
    int comLinksToSurface { -99 };
    bool comOnSurface { false };
    int freeLipids { -99 };
    int numBoundSpecies { -99 };
    int numFreeProtSpecies { -99 };
    bool isDissociated { false };
    std::size_t bndlistSize { 0 };
    double ifaceZ { 0.0 };
};

CdilSnapshot cdil_snapshot(const CdilSystem& s)
{
    CdilSnapshot snap;
    snap.isBound = s.moleculeList[0].interfaceList[0].isBound;
    snap.partnerIndex = s.moleculeList[0].interfaceList[0].interaction.partnerIndex;
    snap.ifaceAbsIndex = s.moleculeList[0].interfaceList[0].index;
    snap.molLinksToSurface = s.moleculeList[0].linksToSurface;
    snap.comLinksToSurface = s.complexList[0].linksToSurface;
    snap.comOnSurface = s.complexList[0].OnSurface;
    snap.freeLipids = s.membraneObject.numberOfFreeLipidsEachState[0];
    snap.numBoundSpecies = s.counterArrays.copyNumSpecies[kCdilBoundState];
    snap.numFreeProtSpecies = s.counterArrays.copyNumSpecies[kCdilProtFreeState];
    snap.isDissociated = s.moleculeList[0].isDissociated;
    snap.bndlistSize = s.moleculeList[0].bndlist.size();
    snap.ifaceZ = s.moleculeList[0].interfaceList[0].coord.z;
    return snap;
}

//! Assert that the system is bit-identical (in the interesting fields) to \p before.
void cdil_expect_unchanged(const CdilSnapshot& before, const CdilSnapshot& after, const std::string& why)
{
    EXPECT_EQ(after.isBound, before.isBound) << why << ": interface must stay bound";
    EXPECT_EQ(after.partnerIndex, before.partnerIndex) << why << ": partner index must be untouched";
    EXPECT_EQ(after.ifaceAbsIndex, before.ifaceAbsIndex) << why << ": interface species index must be untouched";
    EXPECT_EQ(after.molLinksToSurface, before.molLinksToSurface) << why << ": molecule linksToSurface must be untouched";
    EXPECT_EQ(after.comLinksToSurface, before.comLinksToSurface) << why << ": complex linksToSurface must be untouched";
    EXPECT_EQ(after.comOnSurface, before.comOnSurface) << why << ": Complex::OnSurface must be untouched";
    EXPECT_EQ(after.freeLipids, before.freeLipids) << why << ": free implicit-lipid count must be untouched";
    EXPECT_EQ(after.numBoundSpecies, before.numBoundSpecies) << why << ": bound species copy number must be untouched";
    EXPECT_EQ(after.numFreeProtSpecies, before.numFreeProtSpecies) << why << ": free protein copy number must be untouched";
    EXPECT_FALSE(after.isDissociated) << why << ": Molecule::isDissociated must stay false";
    EXPECT_EQ(after.bndlistSize, before.bndlistSize) << why << ": bndlist must be untouched";
    EXPECT_DOUBLE_EQ(after.ifaceZ, before.ifaceZ) << why << ": interface coordinate must be untouched";
}

} // namespace

// =============================================================================
// Test 1: a ghosted molecule (MPI copy owned by another rank) causes an
//         immediate return before anything is inspected.
// =============================================================================
void test_cdil_ghosted_molecule_returns_early()
{
    std::cerr << "\n[TEST] test_cdil_ghosted_molecule_returns_early\n"
              << "  Source file:   check_dissociation_implicitlipid.cpp\n"
              << "  Function:      check_dissociation_implicitlipid\n"
              << "  Scenario:      moleculeList[0].isGhosted == true\n"
              << "  Pass criteria: the function returns immediately; the bond, the\n"
              << "                 copy numbers and the free-lipid count are all\n"
              << "                 left exactly as they were.\n";

    CdilSystem s;
    cdil_build(s);
    cdil_init_rng();

    // Even forcing dissociation must not matter: the ghost check comes first.
    s.params.debugParams.forceDissoc = true;
    s.moleculeList[0].isGhosted = true;

    const CdilSnapshot before = cdil_snapshot(s);
    std::cerr << "  Calling check_dissociation_implicitlipid (molItr = 0)...\n";
    cdil_call(s);
    const CdilSnapshot after = cdil_snapshot(s);

    cdil_expect_unchanged(before, after, "ghosted molecule");
    std::cerr << "  Interface still bound? " << std::boolalpha << after.isBound << '\n';
}

// =============================================================================
// Test 2: the parent complex is not on the membrane surface -> every bond is
//         skipped (only membrane-bound complexes can unbind from the IL).
// =============================================================================
void test_cdil_skips_complex_not_on_surface()
{
    std::cerr << "\n[TEST] test_cdil_skips_complex_not_on_surface\n"
              << "  Source file:   check_dissociation_implicitlipid.cpp\n"
              << "  Function:      check_dissociation_implicitlipid\n"
              << "  Scenario:      complexList[0].OnSurface == false\n"
              << "  Pass criteria: the bond is skipped even with forceDissoc set.\n";

    CdilSystem s;
    cdil_build(s);
    cdil_init_rng();

    s.params.debugParams.forceDissoc = true;
    s.complexList[0].OnSurface = false; // pretend the complex is in solution

    const CdilSnapshot before = cdil_snapshot(s);
    std::cerr << "  Calling check_dissociation_implicitlipid...\n";
    cdil_call(s);
    const CdilSnapshot after = cdil_snapshot(s);

    cdil_expect_unchanged(before, after, "complex not on surface");
}

// =============================================================================
// Test 3: the bound partner is a normal (explicit) molecule, not the implicit
//         lipid -> this routine must ignore the bond entirely.
// =============================================================================
void test_cdil_skips_non_implicitlipid_partner()
{
    std::cerr << "\n[TEST] test_cdil_skips_non_implicitlipid_partner\n"
              << "  Source file:   check_dissociation_implicitlipid.cpp\n"
              << "  Function:      check_dissociation_implicitlipid\n"
              << "  Scenario:      partner molecule has isImplicitLipid == false\n"
              << "  Pass criteria: explicit-partner bonds are handled elsewhere, so\n"
              << "                 nothing here may change.\n";

    CdilSystem s;
    cdil_build(s);
    cdil_init_rng();

    s.params.debugParams.forceDissoc = true;
    s.moleculeList[1].isImplicitLipid = false; // partner is now explicit

    const CdilSnapshot before = cdil_snapshot(s);
    std::cerr << "  Calling check_dissociation_implicitlipid...\n";
    cdil_call(s);
    const CdilSnapshot after = cdil_snapshot(s);

    cdil_expect_unchanged(before, after, "explicit partner");
}

// =============================================================================
// Test 4: an irreversible bond (conjBackRxn == -1) can never dissociate.
// =============================================================================
void test_cdil_skips_irreversible_bond()
{
    std::cerr << "\n[TEST] test_cdil_skips_irreversible_bond\n"
              << "  Source file:   check_dissociation_implicitlipid.cpp\n"
              << "  Function:      check_dissociation_implicitlipid\n"
              << "  Scenario:      interaction.conjBackRxn == -1 (irreversible)\n"
              << "  Pass criteria: no back reaction exists, so the bond survives.\n";

    CdilSystem s;
    cdil_build(s);
    cdil_init_rng();

    s.params.debugParams.forceDissoc = true;
    s.moleculeList[0].interfaceList[0].interaction.conjBackRxn = -1;

    const CdilSnapshot before = cdil_snapshot(s);
    std::cerr << "  Calling check_dissociation_implicitlipid...\n";
    cdil_call(s);
    const CdilSnapshot after = cdil_snapshot(s);

    cdil_expect_unchanged(before, after, "irreversible bond");
}

// =============================================================================
// Test 5: an interface listed in bndlist but flagged isBound == false is
//         defensively skipped.
// =============================================================================
void test_cdil_skips_interface_flagged_unbound()
{
    std::cerr << "\n[TEST] test_cdil_skips_interface_flagged_unbound\n"
              << "  Source file:   check_dissociation_implicitlipid.cpp\n"
              << "  Function:      check_dissociation_implicitlipid\n"
              << "  Scenario:      bndlist contains iface 0 but iface.isBound == false\n"
              << "  Pass criteria: the stale bndlist entry is ignored (no crash, no\n"
              << "                 counter updates).\n";

    CdilSystem s;
    cdil_build(s);
    cdil_init_rng();

    s.params.debugParams.forceDissoc = true;
    s.moleculeList[0].interfaceList[0].isBound = false; // stale bndlist entry

    const CdilSnapshot before = cdil_snapshot(s);
    std::cerr << "  Calling check_dissociation_implicitlipid...\n";
    cdil_call(s);
    const CdilSnapshot after = cdil_snapshot(s);

    // isBound is already false here, so compare the remaining fields.
    EXPECT_EQ(after.partnerIndex, before.partnerIndex) << "partner index must be untouched";
    EXPECT_EQ(after.ifaceAbsIndex, before.ifaceAbsIndex) << "interface species index must be untouched";
    EXPECT_EQ(after.molLinksToSurface, before.molLinksToSurface) << "linksToSurface must be untouched";
    EXPECT_EQ(after.freeLipids, before.freeLipids) << "free implicit-lipid count must be untouched";
    EXPECT_EQ(after.numBoundSpecies, before.numBoundSpecies) << "bound copy number must be untouched";
    EXPECT_FALSE(after.isDissociated) << "Molecule::isDissociated must stay false";
    EXPECT_EQ(after.bndlistSize, before.bndlistSize) << "bndlist must be untouched";
}

// =============================================================================
// Test 6: with a vanishingly small off rate *and* time step the computed
//         dissociation probability cannot beat the uniform random number, so the
//         bond must survive.  This exercises the full 2D->3D probability
//         evaluation path (dissociate3D) without firing the reaction.
// =============================================================================
void test_cdil_negligible_probability_keeps_bond()
{
    std::cerr << "\n[TEST] test_cdil_negligible_probability_keeps_bond\n"
              << "  Source file:   check_dissociation_implicitlipid.cpp\n"
              << "  Function:      check_dissociation_implicitlipid (2D->3D branch)\n"
              << "  Scenario:      kb = 1e-6 us^-1 and timeStep = 1e-9 us, so the\n"
              << "                 dissociation probability is ~1e-15.\n"
              << "  Pass criteria: prob <= rand_gsl() (seed 42) and therefore the\n"
              << "                 bond, the counters and the coordinates survive.\n";

    CdilSystem s;
    cdil_build(s);
    cdil_init_rng();

    s.params.debugParams.forceDissoc = false;
    s.params.timeStep = 1e-9; // absurdly small step
    s.backRxns[0].rateList[0].rate = 1e-6; // absurdly slow unbinding

    const CdilSnapshot before = cdil_snapshot(s);
    std::cerr << "  Calling check_dissociation_implicitlipid...\n";
    cdil_call(s);
    const CdilSnapshot after = cdil_snapshot(s);

    cdil_expect_unchanged(before, after, "negligible dissociation probability");
    std::cerr << "  Bound species copy number is still " << after.numBoundSpecies << '\n';
}

// =============================================================================
// Test 7: forced dissociation (debugParams.forceDissoc == true) sets prob = 1 so
//         the reaction always fires.  Verify every documented side effect.
// =============================================================================
void test_cdil_forced_dissociation_updates_everything()
{
    std::cerr << "\n[TEST] test_cdil_forced_dissociation_updates_everything\n"
              << "  Source file:   check_dissociation_implicitlipid.cpp\n"
              << "  Function:      check_dissociation_implicitlipid (2D->3D branch)\n"
              << "  Scenario:      params.debugParams.forceDissoc == true, one link\n"
              << "                 to the surface, box geometry.\n"
              << "  Pass criteria: the bond is broken, the interface returns to its\n"
              << "                 free state, copy numbers / bound-pair counts /\n"
              << "                 free-lipid counts are updated, the complex leaves\n"
              << "                 the membrane (linksToSurface 1 -> 0, OnSurface\n"
              << "                 false, D.z > 0), the complex is lifted to the 3D\n"
              << "                 reflecting surface and the traj status is set.\n";

    CdilSystem s;
    cdil_build(s);
    cdil_init_rng();

    s.params.debugParams.forceDissoc = true;

    const CdilSnapshot before = cdil_snapshot(s);
    const int boundPairSumBefore
        = std::accumulate(s.counterArrays.nBoundPairs.begin(), s.counterArrays.nBoundPairs.end(), 0);
    const double ifaceZBefore = s.moleculeList[0].interfaceList[0].coord.z;
    const double comZBefore = s.moleculeList[0].comCoord.z;

    std::cerr << "  Before: iface z = " << ifaceZBefore << ", COM z = " << comZBefore
              << ", complex linksToSurface = " << before.comLinksToSurface
              << ", free lipids = " << before.freeLipids << '\n';

    std::cerr << "  Calling check_dissociation_implicitlipid...\n";
    cdil_call(s);

    const CdilSnapshot after = cdil_snapshot(s);
    const int boundPairSumAfter
        = std::accumulate(s.counterArrays.nBoundPairs.begin(), s.counterArrays.nBoundPairs.end(), 0);

    std::cerr << "  After:  iface z = " << s.moleculeList[0].interfaceList[0].coord.z
              << ", COM z = " << s.moleculeList[0].comCoord.z
              << ", complex linksToSurface = " << after.comLinksToSurface
              << ", free lipids = " << after.freeLipids << '\n';

    /* --- the bond itself --- */
    EXPECT_FALSE(after.isBound) << "the protein interface must no longer be bound";
    EXPECT_EQ(after.partnerIndex, -1) << "Interaction::clear() must reset partnerIndex to -1";
    EXPECT_EQ(s.moleculeList[0].interfaceList[0].interaction.partnerIfaceIndex, -1)
        << "Interaction::clear() must reset partnerIfaceIndex to -1";
    EXPECT_EQ(after.ifaceAbsIndex, kCdilProtFreeState)
        << "the interface must revert to the free product species index";
    EXPECT_TRUE(s.moleculeList[0].bndlist.empty())
        << "the dissociated interface must be removed from bndlist";

    /* --- species copy numbers --- */
    EXPECT_EQ(after.numBoundSpecies, before.numBoundSpecies - 1)
        << "the bound species copy number must be decremented";
    EXPECT_EQ(after.numFreeProtSpecies, before.numFreeProtSpecies + 1)
        << "the free protein species copy number must be incremented";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kCdilIlFreeState], 101)
        << "the free implicit-lipid species copy number must be incremented";

    /* --- bound-pair bookkeeping (index scheme agnostic) --- */
    EXPECT_LT(boundPairSumAfter, boundPairSumBefore)
        << "update_Nboundpairs(-1) must reduce the total number of bound pairs";
    std::cerr << "  nBoundPairs sum: " << boundPairSumBefore << " -> " << boundPairSumAfter << '\n';

    /* --- membrane bookkeeping --- */
    EXPECT_EQ(after.freeLipids, before.freeLipids + 1)
        << "one implicit lipid must be returned to the free pool";
    EXPECT_EQ(after.molLinksToSurface, 0) << "the molecule must have no links to the surface left";
    EXPECT_EQ(after.comLinksToSurface, 0) << "the complex must have no links to the surface left";
    EXPECT_FALSE(after.comOnSurface) << "the complex must no longer be flagged as membrane bound";
    EXPECT_GT(s.complexList[0].D.z, 0.0)
        << "after leaving the membrane the complex must regain z-diffusion";

    /* --- geometry: lifted to the 3D reflecting surface, still inside the box --- */
    const double ifaceZAfter = s.moleculeList[0].interfaceList[0].coord.z;
    EXPECT_GT(ifaceZAfter, ifaceZBefore)
        << "the interface must be displaced away from the membrane (toward RS3D)";
    EXPECT_GE(ifaceZAfter, -kCdilBoxLen / 2.0)
        << "the interface must remain inside the water box";
    EXPECT_NEAR(ifaceZAfter, -kCdilBoxLen / 2.0 + kCdilRS3D, 1e-3)
        << "the interface should sit at the 3D reflecting surface (-z/2 + RS3D)";
    EXPECT_NEAR(s.moleculeList[0].comCoord.z - ifaceZAfter, comZBefore - ifaceZBefore, 1e-6)
        << "the whole complex must be translated rigidly (COM-iface vector preserved)";
    EXPECT_TRUE(s.moleculeList[0].tmpICoords.empty())
        << "temporary association coordinates must be cleared afterwards";

    /* --- status flags --- */
    EXPECT_TRUE(after.isDissociated) << "Molecule::isDissociated must be set";
    EXPECT_EQ(s.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "the molecule trajectory status must be 'propagated'";
    EXPECT_EQ(s.complexList[0].trajStatus, TrajStatus::propagated)
        << "the complex trajectory status must be 'propagated'";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario builds its own independent system, so a
// failure in one test cannot corrupt the others.
// -----------------------------------------------------------------------------
TEST(CheckDissociationImplicitLipid, GhostedMoleculeReturnsEarly) { test_cdil_ghosted_molecule_returns_early(); }
TEST(CheckDissociationImplicitLipid, SkipsComplexNotOnSurface) { test_cdil_skips_complex_not_on_surface(); }
TEST(CheckDissociationImplicitLipid, SkipsNonImplicitLipidPartner) { test_cdil_skips_non_implicitlipid_partner(); }
TEST(CheckDissociationImplicitLipid, SkipsIrreversibleBond) { test_cdil_skips_irreversible_bond(); }
TEST(CheckDissociationImplicitLipid, SkipsInterfaceFlaggedUnbound) { test_cdil_skips_interface_flagged_unbound(); }
TEST(CheckDissociationImplicitLipid, NegligibleProbabilityKeepsBond) { test_cdil_negligible_probability_keeps_bond(); }
TEST(CheckDissociationImplicitLipid, ForcedDissociationUpdatesEverything)
{
    test_cdil_forced_dissociation_updates_everything();
}