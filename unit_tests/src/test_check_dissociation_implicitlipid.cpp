/*! \file test_check_dissociation_implicitlipid.cpp
 *
 * ### Unit test for src/reactions/check_dissociation_implicitlipid.cpp
 *
 * The single function under test is
 *
 * \code
 * void check_dissociation_implicitlipid(long long int simItr, const Parameters&,
 *          SimulVolume&, std::vector<MolTemplate>&, std::map<std::string,int>&,
 *          unsigned int molItr, std::vector<Molecule>&, std::vector<Complex>&,
 *          const std::vector<BackRxn>&, const std::vector<ForwardRxn>&,
 *          const std::vector<CreateDestructRxn>&, copyCounters&, Membrane&,
 *          std::vector<double>&, std::vector<double>&, std::vector<double>&,
 *          std::ofstream&);
 * \endcode
 *
 * It walks the bound-interface list (`bndlist`) of one molecule and, for every
 * interface bound to the *implicit lipid*, evaluates an unbinding probability
 * and (stochastically) performs the dissociation.
 *
 * The function contains a number of early-exit / "skip this interface" guards
 * that can be probed completely deterministically:
 *
 *   1. the molecule is ghosted (MPI)                    -> immediate return
 *   2. the molecule has no bound interfaces             -> nothing happens
 *   3. the parent complex is not on the surface         -> interface skipped
 *   4. the bound partner is not an implicit lipid       -> interface skipped
 *   5. the bond came from an irreversible reaction      -> interface skipped
 *
 * plus two "productive" paths:
 *
 *   6. forced dissociation (`params.debugParams.forceDissoc == true`) which
 *      makes the reaction probability unity, so the 2D->3D unbinding is
 *      guaranteed to be executed and every book-keeping side effect can be
 *      checked, and
 *   7. a vanishingly small off rate, for which the reaction probability is
 *      essentially zero so that no dissociation is expected.
 *
 * A minimal but self-consistent system is built for each test: one protein
 * (molTypeIndex 0) whose single interface is bound to an implicit lipid
 * (molTypeIndex 1), a reversible forward/back reaction pair, and a flat
 * (box) membrane geometry.
 */

#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

// The random number generator lives in gtest_main.cpp (as does totMatches which
// the reaction-matching helpers use); we only ever *use* them here.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Container holding every argument required by check_dissociation_implicitlipid.
// -----------------------------------------------------------------------------
struct CdilSystem {
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
    std::vector<double> IL2DbindingVec {};
    std::vector<double> IL2DUnbindingVec {};
    std::vector<double> ILTableIDs {};
};

//! Species (absolute interface state) indices used throughout the tests.
constexpr int kFreeProteinIface { 0 }; //!< free protein interface species
constexpr int kFreeLipidIface { 1 }; //!< free implicit-lipid interface species
constexpr int kBoundPairIface { 2 }; //!< protein--lipid bound species

//! Molecule indices.
constexpr int kProteinMol { 0 };
constexpr int kLipidMol { 1 };

/*! \brief Make sure the global GSL random number generator is usable.
 *
 * check_dissociation_implicitlipid() calls rand_gsl(); the generator pointer is
 * defined (but not initialized) in gtest_main.cpp, so initialize it once.
 */
void cdil_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] GSL RNG was null -> calling srand_gsl(42)\n";
        srand_gsl(42);
    }
}

/*! \brief Build the protein / implicit-lipid MolTemplates.
 *
 * Template 0: soluble protein "A" with a single interface (free state index 0).
 * Template 1: implicit lipid  "IL" with a single interface (free state index 1).
 */
void cdil_build_templates(CdilSystem& s)
{
    // ---------------- protein template ----------------
    MolTemplate proTemp {};
    proTemp.molName = "A";
    proTemp.molTypeIndex = 0;
    proTemp.mass = 1.0;
    proTemp.radius = 1.0;
    proTemp.copies = 1;
    proTemp.D = Coord(1.0, 1.0, 0.0); // membrane bound -> Dz == 0
    proTemp.Dr = Coord(0.01, 0.01, 0.01);
    proTemp.isImplicitLipid = false;
    proTemp.isLipid = false;

    Interface proIface {};
    proIface.name = "a";
    proIface.index = 0;
    proIface.iCoord = Coord(0.0, 0.0, -1.0);
    Interface::State proState {};
    proState.index = kFreeProteinIface;
    proState.iden = '\0';
    proState.ifaceAndStateName = "A(a)";
    proIface.stateList.push_back(proState);
    proTemp.interfaceList.push_back(proIface);

    // ---------------- implicit lipid template ----------------
    MolTemplate ilTemp {};
    ilTemp.molName = "IL";
    ilTemp.molTypeIndex = 1;
    ilTemp.mass = 1.0;
    ilTemp.radius = 1.0;
    ilTemp.copies = 1;
    ilTemp.D = Coord(0.5, 0.5, 0.0);
    ilTemp.Dr = Coord(0.0, 0.0, 0.0);
    ilTemp.isImplicitLipid = true;
    ilTemp.isLipid = true;

    Interface ilIface {};
    ilIface.name = "il";
    ilIface.index = 0;
    ilIface.iCoord = Coord(0.0, 0.0, 0.0);
    Interface::State ilState {};
    ilState.index = kFreeLipidIface;
    ilState.iden = '\0';
    ilState.ifaceAndStateName = "IL(il)";
    ilIface.stateList.push_back(ilState);
    ilTemp.interfaceList.push_back(ilIface);

    s.molTemplateList.clear();
    s.molTemplateList.push_back(proTemp);
    s.molTemplateList.push_back(ilTemp);

    // Statics the simulation classes rely on.
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };
    MolTemplate::absToRelIface = std::vector<int>(8, 0);
    Interface::State::totalNumOfStates = 2;
}

/*! \brief Build the reversible forward/back reaction pair A(a) + IL(il) <-> A(a!1).IL(il!1) */
void cdil_build_reactions(CdilSystem& s, double offRate)
{
    // ------------- forward (binding) reaction -------------
    ForwardRxn fwd {};
    fwd.rxnType = ReactionType::bimolecular;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.isOnMem = true;
    fwd.bindRadius = 1.0;
    fwd.length3Dto2D = 2.0;
    fwd.rateList.push_back(RxnBase::RateState(10.0, std::vector<std::vector<RxnIface>>(1)));
    fwd.reactantListNew.emplace_back("a", 0, kFreeProteinIface, 0, '\0', false);
    fwd.reactantListNew.emplace_back("il", 1, kFreeLipidIface, 0, '\0', false);
    fwd.productListNew.emplace_back("a", 0, kBoundPairIface, 0, '\0', true);
    fwd.productListNew.emplace_back("il", 1, kBoundPairIface, 0, '\0', true);

    // ------------- back (unbinding) reaction -------------
    BackRxn back {};
    back.rxnType = ReactionType::bimolecular;
    back.absRxnIndex = 0;
    back.relRxnIndex = 0;
    back.isOnMem = true;
    back.isCoupled = false;
    back.isObserved = false;
    back.conjForwardRxnIndex = 0;
    back.rateList.push_back(RxnBase::RateState(offRate, std::vector<std::vector<RxnIface>>(1)));
    back.reactantListNew.emplace_back("a", 0, kBoundPairIface, 0, '\0', true);
    back.reactantListNew.emplace_back("il", 1, kBoundPairIface, 0, '\0', true);
    back.productListNew.emplace_back("a", 0, kFreeProteinIface, 0, '\0', false);
    back.productListNew.emplace_back("il", 1, kFreeLipidIface, 0, '\0', false);

    s.forwardRxns.clear();
    s.forwardRxns.push_back(fwd);
    s.backRxns.clear();
    s.backRxns.push_back(back);
}

/*! \brief Build the two molecules (protein + implicit lipid) and their complexes. */
void cdil_build_molecules(CdilSystem& s)
{
    // ------------- protein, bound to the implicit lipid -------------
    Molecule pro {};
    pro.index = kProteinMol;
    pro.id = kProteinMol;
    pro.myComIndex = 0;
    pro.molTypeIndex = 0;
    pro.mass = 1.0;
    pro.isImplicitLipid = false;
    pro.isLipid = false;
    pro.isGhosted = false;
    pro.comCoord = Coord(0.0, 0.0, -49.0);

    Molecule::Iface proIface {};
    proIface.coord = Coord(0.0, 0.0, -50.0); // sitting on the membrane
    proIface.index = kBoundPairIface; // currently in the bound state
    proIface.relIndex = 0;
    proIface.stateIndex = 0;
    proIface.molTypeIndex = 0;
    proIface.isBound = true;
    proIface.interaction = Molecule::Interaction(kLipidMol, 0, 0); // partner, its iface, conjBackRxn
    pro.interfaceList.push_back(proIface);
    pro.bndlist.push_back(0); // interface 0 is bound
    pro.linksToSurface = 1;
    pro.trajStatus = TrajStatus::none;

    // ------------- the implicit lipid -------------
    Molecule lip {};
    lip.index = kLipidMol;
    lip.id = kLipidMol;
    lip.myComIndex = 1;
    lip.molTypeIndex = 1;
    lip.mass = 1.0;
    lip.isImplicitLipid = true;
    lip.isLipid = true;
    lip.isGhosted = false;
    lip.comCoord = Coord(0.0, 0.0, -50.0);

    Molecule::Iface lipIface {};
    lipIface.coord = Coord(0.0, 0.0, -50.0);
    lipIface.index = kBoundPairIface;
    lipIface.relIndex = 0;
    lipIface.stateIndex = 0;
    lipIface.molTypeIndex = 1;
    lipIface.isBound = true;
    lipIface.interaction = Molecule::Interaction(kProteinMol, 0, 0);
    lip.interfaceList.push_back(lipIface);
    lip.bndlist.push_back(0);

    s.moleculeList.clear();
    s.moleculeList.push_back(pro);
    s.moleculeList.push_back(lip);
    Molecule::numberOfMolecules = 2;
    Molecule::emptyMolList.clear();

    // ------------- complexes -------------
    Complex proCom {};
    proCom.index = 0;
    proCom.id = 0;
    proCom.ownerRank = 0;
    proCom.comCoord = pro.comCoord;
    proCom.memberList = std::vector<int> { kProteinMol };
    proCom.numEachMol = std::vector<int> { 1, 0 };
    proCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    proCom.mass = 1.0;
    proCom.radius = 1.0;
    // Deliberately equal to the implicit lipid's D so the source's
    // "std::abs(Dcomplex - DIL) < 1E-10" shortcut is taken (no 1/0 blow-up).
    proCom.D = Coord(0.5, 0.5, 0.0);
    proCom.Dr = Coord(0.01, 0.01, 0.01);
    proCom.OnSurface = true;
    proCom.linksToSurface = 1;
    proCom.iLipidIndex = kLipidMol;
    proCom.trajStatus = TrajStatus::none;

    Complex lipCom {};
    lipCom.index = 1;
    lipCom.id = 1;
    lipCom.ownerRank = 0;
    lipCom.comCoord = lip.comCoord;
    lipCom.memberList = std::vector<int> { kLipidMol };
    lipCom.numEachMol = std::vector<int> { 0, 1 };
    lipCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    lipCom.mass = 1.0;
    lipCom.radius = 1.0;
    lipCom.D = Coord(0.5, 0.5, 0.0);
    lipCom.Dr = Coord(0.0, 0.0, 0.0);
    lipCom.OnSurface = true;
    lipCom.linksToSurface = 1;

    s.complexList.clear();
    s.complexList.push_back(proCom);
    s.complexList.push_back(lipCom);
    Complex::numberOfComplexes = 2;
    Complex::currNumberMolTypes = 2;
    Complex::emptyComList.clear();
}

/*! \brief Fill the copy-number/bound-pair counters. */
void cdil_build_counters(CdilSystem& s)
{
    s.counterArrays.copyNumSpecies.assign(8, 0);
    s.counterArrays.copyNumSpecies[kBoundPairIface] = 1; // one bound pair exists
    s.counterArrays.nBoundPairs.assign(16, 5); // generous, indexing scheme agnostic
    s.counterArrays.proPairlist.assign(16, 0);
    s.counterArrays.singleDouble.assign(8, 0);
    s.counterArrays.implicitDouble.assign(8, false);
    s.counterArrays.canDissociate.assign(8, false);
    s.counterArrays.bindPairList.assign(8, std::vector<int> {});
    s.counterArrays.bindPairListIL2D.assign(8, std::vector<int> {});
    s.counterArrays.bindPairListIL3D.assign(8, std::vector<int> {});
}

/*! \brief Build the whole system.
 *
 * \param[in] offRate  the back (unbinding) rate used for the reaction.
 */
void cdil_build_system(CdilSystem& s, double offRate)
{
    // ---------------- parameters ----------------
    s.params.timeStep = 0.1;
    Parameters::dt = 0.1;
    s.params.numMolTypes = 2;
    s.params.numTotalSpecies = 8;
    s.params.numTotalComplex = 2;
    s.params.implicitLipid = true;
    s.params.debugParams.forceDissoc = false;
    s.params.assocDissocWrite = false;

    // ---------------- membrane ----------------
    s.membraneObject.implicitLipid = true;
    s.membraneObject.TwoD = false;
    s.membraneObject.isSphere = false;
    s.membraneObject.isBox = true;
    s.membraneObject.hasCompartment = false;
    s.membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    s.membraneObject.totalSA = 100.0 * 100.0;
    s.membraneObject.xBCtype = "reflect";
    s.membraneObject.yBCtype = "reflect";
    s.membraneObject.zBCtype = "reflect";
    s.membraneObject.implicitlipidIndex = kLipidMol;
    s.membraneObject.nStates = 1;
    s.membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 10 };
    s.membraneObject.No_free_lipids = 10;
    s.membraneObject.No_protein = 1;
    // The RS3D lookup table is scanned as [i], [i+100], [i+200], [i+300] for
    // i in [0,100) so it must hold at least 400 doubles.
    s.membraneObject.RS3Dvect.assign(400, 0.0);

    cdil_build_templates(s);
    cdil_build_reactions(s, offRate);
    cdil_build_molecules(s);
    cdil_build_counters(s);

    // Make the RS3D lookup succeed: bindRadius / rate / Dtot must match exactly.
    const Coord& d0 = s.molTemplateList[s.forwardRxns[0].reactantListNew[0].molTypeIndex].D;
    const Coord& d1 = s.molTemplateList[s.forwardRxns[0].reactantListNew[1].molTypeIndex].D;
    const double dTot = 1.0 / 3.0 * (d0.x + d1.x) + 1.0 / 3.0 * (d0.y + d1.y) + 1.0 / 3.0 * (d0.z + d1.z);
    s.membraneObject.RS3Dvect[0] = s.forwardRxns[0].bindRadius;
    s.membraneObject.RS3Dvect[100] = s.forwardRxns[0].rateList[0].rate;
    s.membraneObject.RS3Dvect[200] = dTot;
    s.membraneObject.RS3Dvect[300] = 1.0; // the reflecting surface offset

    // A single sub-volume so any book-keeping that touches it is in range.
    s.simulVolume.subCellList.resize(1);
    s.simulVolume.subCellList[0].absIndex = 0;
    s.simulVolume.subCellList[0].memberMolList = std::vector<int> { kProteinMol, kLipidMol };
    s.moleculeList[kProteinMol].mySubVolIndex = 0;
    s.moleculeList[kLipidMol].mySubVolIndex = 0;

    // The dissociation code translates the *temporary* association coordinates
    // and then copies them back, so make sure they exist before the call.
    for (auto& mol : s.moleculeList)
        mol.set_tmp_association_coords();
}

//! Convenience: total number of bound pairs currently recorded.
int cdil_bound_pair_sum(const CdilSystem& s)
{
    return std::accumulate(s.counterArrays.nBoundPairs.begin(), s.counterArrays.nBoundPairs.end(), 0);
}

/*! \brief Assert that absolutely nothing about the bond changed.
 *
 * Pass criteria for all of the "guard" tests: the interface is still bound, the
 * link counters, species copy numbers and free-lipid counts are untouched, and
 * the molecule was not flagged as dissociated.
 */
void cdil_expect_no_dissociation(const CdilSystem& s, const char* context)
{
    std::cerr << "    checking that no dissociation happened (" << context << ")\n";
    EXPECT_TRUE(s.moleculeList[kProteinMol].interfaceList[0].isBound)
        << context << ": interface must still be bound";
    EXPECT_EQ(s.moleculeList[kProteinMol].linksToSurface, 1)
        << context << ": molecule linksToSurface must be unchanged";
    EXPECT_EQ(s.complexList[0].linksToSurface, 1)
        << context << ": complex linksToSurface must be unchanged";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kBoundPairIface], 1)
        << context << ": bound species copy number must be unchanged";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kFreeProteinIface], 0)
        << context << ": free protein copy number must be unchanged";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kFreeLipidIface], 0)
        << context << ": free lipid copy number must be unchanged";
    EXPECT_EQ(s.membraneObject.numberOfFreeLipidsEachState[0], 10)
        << context << ": free implicit lipid count must be unchanged";
    EXPECT_FALSE(s.moleculeList[kProteinMol].isDissociated)
        << context << ": molecule must not be flagged as dissociated";
    EXPECT_EQ(cdil_bound_pair_sum(s), 16 * 5)
        << context << ": nBoundPairs must be unchanged";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a ghosted molecule (MPI copy from a neighbouring rank) must cause an
//         immediate return, even though the setup would otherwise dissociate.
// -----------------------------------------------------------------------------
void test_cdil_ghosted_molecule_returns_early()
{
    std::cerr << "\n[TEST] test_cdil_ghosted_molecule_returns_early\n"
              << "  Source file:   src/reactions/check_dissociation_implicitlipid.cpp\n"
              << "  Function:      check_dissociation_implicitlipid()\n"
              << "  Scenario:      protein bound to the implicit lipid, forceDissoc = true,\n"
              << "                 but the molecule is flagged isGhosted = true.\n"
              << "  Pass criteria: the function returns before doing anything.\n";

    cdil_ensure_rng();
    CdilSystem s;
    cdil_build_system(s, 1.0);
    s.params.debugParams.forceDissoc = true; // would guarantee dissociation
    s.moleculeList[kProteinMol].isGhosted = true; // ... but we are ghosted

    std::ofstream assocDissocFile("/dev/null");

    std::cerr << "  Calling check_dissociation_implicitlipid(molItr = 0)...\n";
    check_dissociation_implicitlipid(1, s.params, s.simulVolume, s.molTemplateList, s.observablesList, 0,
        s.moleculeList, s.complexList, s.backRxns, s.forwardRxns, s.createDestructRxns, s.counterArrays,
        s.membraneObject, s.IL2DbindingVec, s.IL2DUnbindingVec, s.ILTableIDs, assocDissocFile);

    cdil_expect_no_dissociation(s, "ghosted molecule");
}

// -----------------------------------------------------------------------------
// Test 2: a molecule with no bound interfaces must be a complete no-op.
// -----------------------------------------------------------------------------
void test_cdil_no_bound_interfaces_is_noop()
{
    std::cerr << "\n[TEST] test_cdil_no_bound_interfaces_is_noop\n"
              << "  Function:      check_dissociation_implicitlipid()\n"
              << "  Scenario:      the protein's bndlist is empty (nothing is bound).\n"
              << "  Pass criteria: the bndlist loop never executes; state unchanged.\n";

    cdil_ensure_rng();
    CdilSystem s;
    cdil_build_system(s, 1.0);
    s.params.debugParams.forceDissoc = true;
    s.moleculeList[kProteinMol].bndlist.clear(); // nothing to iterate over

    std::ofstream assocDissocFile("/dev/null");

    std::cerr << "  Calling check_dissociation_implicitlipid with an empty bndlist...\n";
    check_dissociation_implicitlipid(1, s.params, s.simulVolume, s.molTemplateList, s.observablesList, 0,
        s.moleculeList, s.complexList, s.backRxns, s.forwardRxns, s.createDestructRxns, s.counterArrays,
        s.membraneObject, s.IL2DbindingVec, s.IL2DUnbindingVec, s.ILTableIDs, assocDissocFile);

    cdil_expect_no_dissociation(s, "empty bndlist");
}

// -----------------------------------------------------------------------------
// Test 3: an interface which is listed in bndlist but flagged as not bound must
//         be skipped by the "isBound" guard.
// -----------------------------------------------------------------------------
void test_cdil_unbound_interface_is_skipped()
{
    std::cerr << "\n[TEST] test_cdil_unbound_interface_is_skipped\n"
              << "  Function:      check_dissociation_implicitlipid()\n"
              << "  Scenario:      interface 0 is in bndlist but isBound == false.\n"
              << "  Pass criteria: the interface is skipped, nothing is modified.\n";

    cdil_ensure_rng();
    CdilSystem s;
    cdil_build_system(s, 1.0);
    s.params.debugParams.forceDissoc = true;
    s.moleculeList[kProteinMol].interfaceList[0].isBound = false; // stale bndlist entry

    std::ofstream assocDissocFile("/dev/null");

    std::cerr << "  Calling check_dissociation_implicitlipid with isBound = false...\n";
    check_dissociation_implicitlipid(1, s.params, s.simulVolume, s.molTemplateList, s.observablesList, 0,
        s.moleculeList, s.complexList, s.backRxns, s.forwardRxns, s.createDestructRxns, s.counterArrays,
        s.membraneObject, s.IL2DbindingVec, s.IL2DUnbindingVec, s.ILTableIDs, assocDissocFile);

    // The interface is (still) not bound - check the remaining counters only.
    std::cerr << "    checking counters were not touched\n";
    EXPECT_EQ(s.counterArrays.copyNumSpecies[kBoundPairIface], 1)
        << "bound species copy number must be unchanged";
    EXPECT_EQ(s.membraneObject.numberOfFreeLipidsEachState[0], 10)
        << "free implicit lipid count must be unchanged";
    EXPECT_EQ(s.moleculeList[kProteinMol].linksToSurface, 1)
        << "linksToSurface must be unchanged";
    EXPECT_FALSE(s.moleculeList[kProteinMol].isDissociated)
        << "molecule must not be flagged as dissociated";
}

// -----------------------------------------------------------------------------
// Test 4: the parent complex is not on the membrane surface -> skip.
// -----------------------------------------------------------------------------
void test_cdil_complex_not_on_surface_is_skipped()
{
    std::cerr << "\n[TEST] test_cdil_complex_not_on_surface_is_skipped\n"
              << "  Function:      check_dissociation_implicitlipid()\n"
              << "  Scenario:      complex.OnSurface == false (only surface-bound\n"
              << "                 complexes can unbind from the implicit lipid).\n"
              << "  Pass criteria: the interface is skipped by the OnSurface guard.\n";

    cdil_ensure_rng();
    CdilSystem s;
    cdil_build_system(s, 