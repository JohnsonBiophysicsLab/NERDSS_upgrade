/*! \file test_associate_ImplicitLipid_box.cpp
 *
 * ### Unit test for src/reactions/associate_ImplicitLipid_box.cpp
 *
 * The single function under test is
 *
 *     void associate_implicitlipid_box(long long int iter,
 *                                      int ifaceIndex1, int ifaceIndex2,
 *                                      Molecule& reactMol1, Molecule& reactMol2,
 *                                      Complex& reactCom1, Complex& reactCom2,
 *                                      const Parameters& params, ForwardRxn& currRxn,
 *                                      std::vector<Molecule>& moleculeList,
 *                                      std::vector<MolTemplate>& molTemplateList,
 *                                      std::map<std::string,int>& observablesList,
 *                                      copyCounters& counterArrays,
 *                                      std::vector<Complex>& complexList,
 *                                      Membrane& membraneObject,
 *                                      const std::vector<ForwardRxn>& forwardRxns,
 *                                      const std::vector<BackRxn>& backRxns,
 *                                      std::ofstream& assocDissocFile)
 *
 * The routine binds a real molecule to the *implicit lipid* (IL) in a
 * rectangular ("box") geometry.  It has three interesting code paths:
 *
 *   1. Both complexes are already on the membrane (`OnSurface == true` for
 *      both) -> `isOnMembrane == true`, the whole re-orientation /
 *      reflection / rejection block is skipped and only the bookkeeping is
 *      performed.
 *   2. The binder is in solution -> `transitionToSurface == true`; the binder
 *      is translated down to the membrane so that the interface separation
 *      equals `bindRadius + RS3D`.  If the binder is a *point* molecule all
 *      angular rotations are skipped, which is the variant exercised here.
 *   3. Which of the two reactants is the implicit lipid selects one of two
 *      nearly mirror-image bookkeeping blocks (only the second one writes to
 *      the association/dissociation log file).
 *
 * Every test below builds a *complete* two-molecule system (protein + implicit
 * lipid) so that no member read by the routine (or by the helpers it calls)
 * is left uninitialized.
 */

#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Fixed indices / sizes shared by all tests in this file.
// -----------------------------------------------------------------------------
constexpr int kAilProtTypeIdx = 0;   //!< molTemplateList index of the protein
constexpr int kAilIlTypeIdx = 1;     //!< molTemplateList index of the implicit lipid
constexpr int kAilProtMolIdx = 0;    //!< moleculeList index of the protein
constexpr int kAilIlMolIdx = 1;      //!< moleculeList index of the implicit lipid
constexpr int kAilProtAbsIface = 0;  //!< absolute species index, free protein iface
constexpr int kAilIlAbsIface = 1;    //!< absolute species index, free IL iface
constexpr int kAilProdAbsIface = 2;  //!< absolute species index of the bound product
constexpr double kAilBoxLen = 100.0; //!< cubic water box edge length (nm)
constexpr double kAilBindRadius = 1.0;
constexpr double kAilRate = 10.0;
constexpr int kAilConjBackRxn = 7;   //!< arbitrary "conjugate back reaction" index
constexpr int kAilStartFreeLipids = 100;

/*! \brief Everything the routine under test needs, kept together so that a
 *         test can be written in a handful of lines.
 */
struct AilSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
    ForwardRxn rxn {}; //!< the reaction object actually handed to the routine
};

/*! \brief Initialise the global GSL generator once (some helper routines called
 *         from the association code may resample random numbers).
 */
void ailbox_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}

/*! \brief Sum of all association-event histograms.
 *
 * `track_association_events()` bins a successful association into one of
 * events2D / events3D / events3Dto2D.  Summing all three makes the assertion
 * independent of the exact binning rule.
 */
long ailbox_total_events(const copyCounters& c)
{
    long tot { 0 };
    for (auto v : c.events2D)
        tot += v;
    for (auto v : c.events3D)
        tot += v;
    for (auto v : c.events3Dto2D)
        tot += v;
    return tot;
}

/*! \brief Build a fully initialised protein + implicit-lipid system.
 *
 * \param[out] sys                     system to fill in
 * \param[in]  proteinComplexOnSurface if true the protein's complex is flagged
 *                                     as being on the membrane, which selects
 *                                     the 2D (no re-orientation) code path
 * \param[in]  proteinIsPoint          marks the protein MolTemplate as a point
 *                                     particle (skips all angular rotations)
 * \param[in]  ilIsFirstReactant       ordering of currRxn.reactantListNew
 * \param[in]  proteinCom              starting COM of the protein
 */
void ailbox_build_system(AilSystem& sys, bool proteinComplexOnSurface, bool proteinIsPoint,
    bool ilIsFirstReactant, const Coord& proteinCom)
{
    const double nanv { std::numeric_limits<double>::quiet_NaN() };

    // ---------------------------------------------------------------------
    // Simulation parameters.  scaleMaxDisplace is set huge so that
    // measure_complex_displacement() can never reject the move.
    // ---------------------------------------------------------------------
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 6;
    sys.params.numTotalComplex = 2;
    sys.params.timeStep = 0.1;
    sys.params.overlapSepLimit = 0.0;
    sys.params.scaleMaxDisplace = 1.0e10;
    sys.params.implicitLipid = true;
    sys.params.debugParams.verbosity = 0;

    // ---------------------------------------------------------------------
    // Molecule templates: [0] = protein "A", [1] = implicit lipid "IL".
    // ---------------------------------------------------------------------
    MolTemplate prot;
    prot.molName = "A";
    prot.molTypeIndex = kAilProtTypeIdx;
    prot.copies = 1;
    prot.mass = 1.0;
    prot.radius = 1.0;
    prot.D = Coord { 1.0, 1.0, 1.0 };
    prot.Dr = Coord { 0.01, 0.01, 0.01 };
    prot.checkOverlap = false; // keeps the system-wide overlap check trivial
    prot.canDestroy = false;   // keeps the monomerList branch inactive
    prot.isPoint = proteinIsPoint;
    prot.isLipid = false;
    prot.isImplicitLipid = false;
    {
        Interface iface;
        iface.index = 0;
        iface.name = "a";
        iface.iCoord = Coord { 0.0, 0.0, 0.0 };
        Interface::State st;
        st.index = kAilProtAbsIface;
        st.iden = '\0';
        st.ifaceAndStateName = "a";
        iface.stateList.push_back(st);
        prot.interfaceList.push_back(iface);
    }
    sys.molTemplateList.push_back(prot);

    MolTemplate il;
    il.molName = "IL";
    il.molTypeIndex = kAilIlTypeIdx;
    il.copies = 1;
    il.mass = 1.0;
    il.radius = 1.0;
    il.D = Coord { 0.1, 0.1, 0.0 }; // lipid: no motion normal to the membrane
    il.Dr = Coord { 0.01, 0.01, 0.01 };
    il.checkOverlap = false;
    il.canDestroy = false;
    il.isPoint = false;
    il.isLipid = true;
    il.isImplicitLipid = true;
    {
        Interface iface;
        iface.index = 0;
        iface.name = "il";
        iface.iCoord = Coord { 0.0, 0.0, 0.0 };
        Interface::State st;
        // NOTE: the routine looks the implicit-lipid state up by absolute index
        // at the very end; this value MUST match the IL reactant of the
        // reaction or the state search would return -1.
        st.index = kAilIlAbsIface;
        st.iden = '\0';
        st.ifaceAndStateName = "il";
        iface.stateList.push_back(st);
        il.interfaceList.push_back(iface);
    }
    sys.molTemplateList.push_back(il);

    // ---------------------------------------------------------------------
    // Membrane: flat reflecting box with an implicit lipid.
    // ---------------------------------------------------------------------
    sys.membrane.isSphere = false;
    sys.membrane.isBox = true;
    sys.membrane.hasCompartment = false;
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { kAilBoxLen, kAilBoxLen, kAilBoxLen });
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";
    sys.membrane.implicitLipid = true;
    sys.membrane.implicitlipidIndex = kAilIlMolIdx;
    sys.membrane.nStates = 1;
    sys.membrane.nSites = 0;
    sys.membrane.No_free_lipids = kAilStartFreeLipids;
    sys.membrane.No_protein = 1;
    sys.membrane.numberOfFreeLipidsEachState = std::vector<int> { kAilStartFreeLipids };
    sys.membrane.numberOfProteinEachState = std::vector<int> { 1 };
    sys.membrane.totalSA = kAilBoxLen * kAilBoxLen;
    sys.membrane.Dx = 0.0;
    sys.membrane.Dy = 0.0;
    sys.membrane.Dz = 0.0;
    sys.membrane.Drx = 0.0;
    sys.membrane.Dry = 0.0;
    sys.membrane.Drz = 0.0;
    sys.membrane.offset = 0.0;
    sys.membrane.lipidLength = 0.0;

    // The RS3D lookup table.  The routine scans entries [0..99] for a matching
    // (bindRadius, rate, Dtot) triple and then reads the reflecting distance
    // from entry+300.  We compute Dtot with *exactly* the same expression so
    // that the 1E-15 comparisons succeed and RS3D becomes 0.0 (instead of the
    // "not found" sentinel of -1).
    sys.membrane.RS3Dvect.assign(500, 0.0);
    const double dTot = 1.0 / 3.0 * (prot.D.x + il.D.x) + 1.0 / 3.0 * (prot.D.y + il.D.y)
        + 1.0 / 3.0 * (prot.D.z + il.D.z);
    sys.membrane.RS3Dvect[0] = kAilBindRadius;
    sys.membrane.RS3Dvect[100] = kAilRate;
    sys.membrane.RS3Dvect[200] = dTot;
    sys.membrane.RS3Dvect[300] = 0.0; // <- RS3D that will be picked up
    sys.membrane.RS3Dvect[400] = static_cast<double>(kAilProtTypeIdx);

    // ---------------------------------------------------------------------
    // Molecules: [0] = protein, [1] = implicit lipid.
    // ---------------------------------------------------------------------
    Molecule protMol;
    protMol.index = kAilProtMolIdx;
    protMol.id = kAilProtMolIdx;
    protMol.myComIndex = 0;
    protMol.complexId = 0;
    protMol.molTypeIndex = kAilProtTypeIdx;
    protMol.mass = 1.0;
    protMol.isLipid = false;
    protMol.isImplicitLipid = false;
    protMol.isEmpty = false;
    protMol.linksToSurface = 0;
    protMol.trajStatus = TrajStatus::none;
    protMol.comCoord = proteinCom;
    {
        Molecule::Iface pif;
        pif.coord = proteinCom; // point molecule: interface sits on the COM
        pif.index = kAilProtAbsIface;
        pif.relIndex = 0;
        pif.stateIndex = 0;
        pif.stateIden = '\0';
        pif.molTypeIndex = kAilProtTypeIdx;
        pif.isBound = false;
        protMol.interfaceList.push_back(pif);
    }
    protMol.freelist.push_back(0); // interface 0 is free -> required, the
                                   // routine pops it out of freelist
    sys.moleculeList.push_back(protMol);

    Molecule ilMol;
    ilMol.index = kAilIlMolIdx;
    ilMol.id = kAilIlMolIdx;
    ilMol.myComIndex = 1;
    ilMol.complexId = 1;
    ilMol.molTypeIndex = kAilIlTypeIdx;
    ilMol.mass = 1.0;
    ilMol.isLipid = true;
    ilMol.isImplicitLipid = true;
    ilMol.isEmpty = false;
    ilMol.linksToSurface = 0;
    ilMol.trajStatus = TrajStatus::none;
    ilMol.comCoord = Coord { 0.0, 0.0, -kAilBoxLen / 2.0 };
    {
        Molecule::Iface iif;
        iif.coord = ilMol.comCoord;
        iif.index = kAilIlAbsIface;
        iif.relIndex = 0;
        iif.stateIndex = 0;
        iif.stateIden = '\0';
        iif.molTypeIndex = kAilIlTypeIdx;
        iif.isBound = false;
        ilMol.interfaceList.push_back(iif);
    }
    ilMol.freelist.push_back(0);
    sys.moleculeList.push_back(ilMol);

    // ---------------------------------------------------------------------
    // Complexes: [0] holds the protein, [1] holds the implicit lipid.
    // ---------------------------------------------------------------------
    Complex protCom;
    protCom.index = 0;
    protCom.id = 0;
    protCom.ownerRank = 0;
    protCom.comCoord = proteinCom;
    protCom.memberList = std::vector<int> { kAilProtMolIdx };
    protCom.numEachMol = std::vector<int> { 1, 0 };
    protCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    protCom.mass = 1.0;
    protCom.radius = prot.radius;
    protCom.D = prot.D;
    protCom.Dr = prot.Dr;
    protCom.isEmpty = false;
    protCom.OnSurface = proteinComplexOnSurface;
    protCom.onFiber = false;
    protCom.linksToSurface = 0;
    protCom.iLipidIndex = 0;
    protCom.ncross = 0;
    protCom.trajStatus = TrajStatus::none;
    sys.complexList.push_back(protCom);

    Complex ilCom;
    ilCom.index = 1;
    ilCom.id = 1;
    ilCom.ownerRank = 0;
    ilCom.comCoord = sys.moleculeList[kAilIlMolIdx].comCoord;
    ilCom.memberList = std::vector<int> { kAilIlMolIdx };
    ilCom.numEachMol = std::vector<int> { 0, 1 };
    ilCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    ilCom.mass = 1.0;
    ilCom.radius = il.radius;
    ilCom.D = il.D;
    ilCom.Dr = il.Dr;
    ilCom.isEmpty = false;
    ilCom.OnSurface = true; // the implicit lipid is always on the membrane
    ilCom.onFiber = false;
    ilCom.linksToSurface = 0;
    ilCom.iLipidIndex = kAilIlMolIdx;
    ilCom.ncross = 0;
    ilCom.trajStatus = TrajStatus::none;
    sys.complexList.push_back(ilCom);

    // ---------------------------------------------------------------------
    // The forward reaction A(a) + IL(il) -> A(a!1).IL(il!1)
    // ---------------------------------------------------------------------
    sys.rxn.rxnType = ReactionType::bimolecular;
    sys.rxn.absRxnIndex = 0;
    sys.rxn.relRxnIndex = 0;
    sys.rxn.bindRadius = kAilBindRadius;
    sys.rxn.bindRadius2D = kAilBindRadius;
    sys.rxn.isOnMem = true;
    sys.rxn.isReversible = true;
    sys.rxn.conjBackRxnIndex = kAilConjBackRxn;
    sys.rxn.isObserved = true;
    sys.rxn.observeLabel = "ILbond";
    sys.rxn.norm1 = Vector(0.0, 0.0, 1.0);
    sys.rxn.norm2 = Vector(0.0, 0.0, 1.0);
    // The association angles are never consumed by the code paths tested here
    // (2D binding and point-particle 3D->2D binding both skip all rotations),
    // but they are given sane values anyway.
    sys.rxn.assocAngles = ForwardRxn::Angles(M_PI, M_PI, nanv, nanv, nanv);
    sys.rxn.rateList.emplace_back();
    sys.rxn.rateList.back().rate = kAilRate;

    if (ilIsFirstReactant) {
        sys.rxn.reactantListNew.emplace_back("il", kAilIlTypeIdx, kAilIlAbsIface, 0, '\0', false);
        sys.rxn.reactantListNew.emplace_back("a", kAilProtTypeIdx, kAilProtAbsIface, 0, '\0', false);
    } else {
        sys.rxn.reactantListNew.emplace_back("a", kAilProtTypeIdx, kAilProtAbsIface, 0, '\0', false);
        sys.rxn.reactantListNew.emplace_back("il", kAilIlTypeIdx, kAilIlAbsIface, 0, '\0', false);
    }
    sys.rxn.productListNew.emplace_back("a", kAilProtTypeIdx, kAilProdAbsIface, 0, '\0', true);
    sys.rxn.productListNew.emplace_back("il", kAilIlTypeIdx, kAilProdAbsIface, 0, '\0', true);
    sys.rxn.intReactantList = std::vector<int> { kAilProtAbsIface, kAilIlAbsIface };
    sys.rxn.intProductList = std::vector<int> { kAilProdAbsIface };

    sys.forwardRxns.push_back(sys.rxn); // used only by the overlap checks
    sys.observablesList["ILbond"] = 0;

    // ---------------------------------------------------------------------
    // Counters.  Everything the routine (or update_Nboundpairs /
    // track_association_events) may index into is generously sized.
    // ---------------------------------------------------------------------
    sys.counterArrays.copyNumSpecies = std::vector<int> { 10, kAilStartFreeLipids, 0, 0, 0, 0 };
    sys.counterArrays.nBoundPairs.assign(16, 0);
    sys.counterArrays.proPairlist.assign(16, 0);
    sys.counterArrays.singleDouble.assign(6, 0);
    sys.counterArrays.implicitDouble.assign(6, false);
    sys.counterArrays.canDissociate.assign(6, false);
    sys.counterArrays.bindPairList.resize(6);
    sys.counterArrays.bindPairListIL2D.resize(6);
    sys.counterArrays.bindPairListIL3D.resize(6);
    sys.counterArrays.eventArraySize = 20;
    init_association_events(sys.counterArrays);
    const size_t needed { static_cast<size_t>(sys.counterArrays.eventArraySize) + 1 };
    if (sys.counterArrays.events2D.size() < needed)
        sys.counterArrays.events2D.resize(needed, 0);
    if (sys.counterArrays.events3D.size() < needed)
        sys.counterArrays.events3D.resize(needed, 0);
    if (sys.counterArrays.events3Dto2D.size() < needed)
        sys.counterArrays.events3Dto2D.resize(needed, 0);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: 2D binding, the SECOND reactant is the implicit lipid.
//
//   * both complexes are flagged OnSurface -> isOnMembrane == true, so the
//     re-orientation / reflection / rejection block is bypassed entirely.
//   * the first bookkeeping branch (reactMol2.isImplicitLipid == true) runs.
// -----------------------------------------------------------------------------
void test_ailbox_2d_mol2_is_implicit_lipid()
{
    std::cerr << "\n[TEST] test_ailbox_2d_mol2_is_implicit_lipid\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_box.cpp\n"
              << "  Function:    associate_implicitlipid_box()\n"
              << "  Scenario:    membrane-bound protein binds the implicit lipid\n"
              << "               (reactMol2 == IL, isOnMembrane == true)\n"
              << "  Checks:      bond bookkeeping, copy numbers, free-lipid count,\n"
              << "               complex properties and observable increment.\n";

    ailbox_init_rng();

    AilSystem sys;
    const Coord protStart { 0.0, 0.0, -45.0 };
    // proteinComplexOnSurface = true  -> 2D path
    // proteinIsPoint          = false -> irrelevant here (no rotations happen)
    // ilIsFirstReactant       = false -> reactantListNew = {protein, IL}
    ailbox_build_system(sys, /*onSurface*/ true, /*isPoint*/ false, /*ilFirst*/ false, protStart);

    const long eventsBefore = ailbox_total_events(sys.counterArrays);
    std::ofstream noLogFile; // deliberately closed: this branch never writes

    std::cerr << "  Calling associate_implicitlipid_box()...\n";
    associate_implicitlipid_box(/*iter*/ 1, /*ifaceIndex1*/ 0, /*ifaceIndex2*/ 0,
        sys.moleculeList[kAilProtMolIdx], sys.moleculeList[kAilIlMolIdx],
        sys.complexList[0], sys.complexList[1], sys.params, sys.rxn, sys.moleculeList,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.complexList,
        sys.membrane, sys.forwardRxns, sys.backRxns, noLogFile);

    const Molecule& protein = sys.moleculeList[kAilProtMolIdx];
    const Molecule& lipid = sys.moleculeList[kAilIlMolIdx];
    const Complex& protComplex = sys.complexList[0];

    // --- the association must have been accepted --------------------------
    std::cerr << "  nAssocSuccess = " << sys.counterArrays.nAssocSuccess << '\n';
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "one successful association expected";
    EXPECT_EQ(sys.counterArrays.nCancelSpanBox, 0) << "2D path must not run the span check";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapSystem, 0);
    EXPECT_EQ(sys.counterArrays.nCancelDisplace2D, 0);
    EXPECT_EQ(ailbox_total_events(sys.counterArrays) - eventsBefore, 1)
        << "track_association_events() should have binned exactly one event";

    // --- interface / interaction bookkeeping on the protein ---------------
    EXPECT_TRUE(protein.interfaceList[0].isBound) << "protein interface must be flagged bound";
    EXPECT_EQ(protein.interfaceList[0].interaction.partnerIndex, kAilIlMolIdx)
        << "partner must be the implicit lipid molecule";
    EXPECT_EQ(protein.interfaceList[0].interaction.partnerIfaceIndex, 0);
    EXPECT_EQ(protein.interfaceList[0].interaction.conjBackRxn, kAilConjBackRxn)
        << "reversible reaction stores its conjugate back-reaction index";
    EXPECT_EQ(protein.interfaceList[0].index, kAilProdAbsIface)
        << "interface must adopt the product species index";
    ASSERT_EQ(protein.bndlist.size(), 1u);
    EXPECT_EQ(protein.bndlist[0], 0);
    ASSERT_EQ(protein.bndpartner.size(), 1u);
    EXPECT_EQ(protein.bndpartner[0], kAilIlMolIdx);
    EXPECT_TRUE(protein.freelist.empty()) << "the bound interface is removed from freelist";
    EXPECT_EQ(protein.linksToSurface, 1) << "protein gains one link to the surface";
    EXPECT_EQ(protein.trajStatus, TrajStatus::propagated);

    // --- complex level bookkeeping ---------------------------------------
    EXPECT_EQ(protComplex.iLipidIndex, kAilIlMolIdx);
    EXPECT_EQ(protComplex.linksToSurface, 1);
    EXPECT_TRUE(protComplex.OnSurface) << "complex must be marked as membrane bound";
    EXPECT_DOUBLE_EQ(protComplex.D.z, 0.0) << "membrane bound complexes cannot diffuse in z";
    EXPECT_EQ(protComplex.ncross, -1) << "complex is removed from the overlap sweep";
    EXPECT_TRUE(protein.crossbase.empty());

    // --- temporary association coordinates must always be released -------
    EXPECT_TRUE(protein.tmpICoords.empty()) << "protein tmp coords must be cleared";
    EXPECT_TRUE(lipid.tmpICoords.empty()) << "implicit lipid tmp coords must be cleared";

    // --- species / lipid counters ----------------------------------------
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilProtAbsIface], 9)
        << "free protein interface count decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilIlAbsIface], kAilStartFreeLipids - 1)
        << "free IL interface count decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilProdAbsIface], 1)
        << "bound product count incremented";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], kAilStartFreeLipids - 1)
        << "one implicit lipid of the reacting state is consumed";

    long boundPairs { 0 };
    for (auto v : sys.counterArrays.nBoundPairs)
        boundPairs += v;
    EXPECT_GT(boundPairs, 0) << "update_Nboundpairs() must record the new pair";

    // --- observable tracking ---------------------------------------------
    EXPECT_EQ(sys.observablesList["ILbond"], 1) << "observed reaction increments its label";

    // --- geometry: the implicit lipid was repositioned underneath the
    //     binder by create_position_implicit_lipid() (shift of +/-0.1 nm in
    //     x/y and pinned to the bottom of the box in z).
    std::cerr << "  IL COM after call = (" << lipid.comCoord.x << ", " << lipid.comCoord.y
              << ", " << lipid.comCoord.z << ")\n";
    EXPECT_NEAR(lipid.comCoord.x, protStart.x + 0.1, 1e-12);
    EXPECT_NEAR(lipid.comCoord.y, protStart.y - 0.1, 1e-12);
    EXPECT_NEAR(lipid.comCoord.z, -kAilBoxLen / 2.0, 1e-12)
        << "the implicit lipid sits on the bottom membrane of the box";

    // --- the binder itself is not re-oriented on the 2D path; it must at
    //     the very least still be inside the simulation box.
    std::cerr << "  protein COM after call = (" << protein.comCoord.x << ", "
              << protein.comCoord.y << ", " << protein.comCoord.z << ")\n";
    EXPECT_NEAR(protein.comCoord.x, protStart.x, 1e-9) << "no lateral motion expected in 2D";
    EXPECT_NEAR(protein.comCoord.y, protStart.y, 1e-9) << "no lateral motion expected in 2D";
    EXPECT_GE(protein.comCoord.z, -kAilBoxLen / 2.0 - 1e-9) << "protein must stay in the box";
    EXPECT_LE(protein.comCoord.z, kAilBoxLen / 2.0 + 1e-9) << "protein must stay in the box";
}

// -----------------------------------------------------------------------------
// Test 2: 2D binding, the FIRST reactant is the implicit lipid.
//
// This exercises the mirror bookkeeping branch, which is the only one that
// writes a line into the association/dissociation log file.
// -----------------------------------------------------------------------------
void test_ailbox_2d_mol1_is_implicit_lipid()
{
    std::cerr << "\n[TEST] test_ailbox_2d_mol1_is_implicit_lipid\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_box.cpp\n"
              << "  Function:    associate_implicitlipid_box()\n"
              << "  Scenario:    reactMol1 == implicit lipid (mirror branch),\n"
              << "               assocDissocFile is open so a BOND record is written.\n"
              << "  Checks:      the partner molecule is bound and the log line has\n"
              << "               the documented 'ITR:<iter>,BOND,...' format.\n";

    ailbox_init_rng();

    AilSystem sys;
    const Coord protStart { 5.0, -5.0, -40.0 };
    // ilIsFirstReactant = true -> reactantListNew = {IL, protein}
    ailbox_build_system(sys, /*onSurface*/ true, /*isPoint*/ false, /*ilFirst*/ true, protStart);

    const char* logName = "test_associate_implicitlipid_box_assocdissoc.dat";
    std::ofstream logFile(logName);
    ASSERT_TRUE(logFile.is_open()) << "could not open temporary log file";

    std::cerr << "  Calling associate_implicitlipid_box() with reactMol1 = IL...\n";
    // NOTE the swapped argument order: reactMol1/reactCom1 are the implicit
    // lipid and its complex, reactMol2/reactCom2 are the protein.
    associate_implicitlipid_box(/*iter*/ 5, /*ifaceIndex1*/ 0, /*ifaceIndex2*/ 0,
        sys.moleculeList[kAilIlMolIdx], sys.moleculeList[kAilProtMolIdx],
        sys.complexList[1], sys.complexList[0], sys.params, sys.rxn, sys.moleculeList,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.complexList,
        sys.membrane, sys.forwardRxns, sys.backRxns, logFile);

    logFile.close();

    const Molecule& protein = sys.moleculeList[kAilProtMolIdx];
    const Complex& protComplex = sys.complexList[0];

    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "association should succeed";

    // The protein (reactMol2 here) receives all of the bond bookkeeping.
    EXPECT_TRUE(protein.interfaceList[0].isBound);
    EXPECT_EQ(protein.interfaceList[0].interaction.partnerIndex, kAilIlMolIdx);
    EXPECT_EQ(protein.interfaceList[0].interaction.partnerIfaceIndex, 0);
    EXPECT_EQ(protein.interfaceList[0].interaction.conjBackRxn, kAilConjBackRxn);
    EXPECT_EQ(protein.interfaceList[0].index, kAilProdAbsIface);
    EXPECT_EQ(protein.linksToSurface, 1);
    EXPECT_TRUE(protein.freelist.empty());
    ASSERT_EQ(protein.bndpartner.size(), 1u);
    EXPECT_EQ(protein.bndpartner[0], kAilIlMolIdx);

    EXPECT_EQ(protComplex.iLipidIndex, kAilIlMolIdx);
    EXPECT_EQ(protComplex.linksToSurface, 1);
    EXPECT_TRUE(protComplex.OnSurface);
    EXPECT_DOUBLE_EQ(protComplex.D.z, 0.0);
    EXPECT_EQ(protComplex.ncross, -1);

    // Even with the reactants swapped the implicit-lipid pool is reduced.
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], kAilStartFreeLipids - 1);
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilProdAbsIface], 1);

    // --- verify the log record -------------------------------------------
    std::ifstream in(logName);
    std::stringstream buffer;
    buffer << in.rdbuf();
    in.close();
    const std::string contents { buffer.str() };
    std::remove(logName);

    std::cerr << "  log file contents: " << contents;
    EXPECT_NE(contents.find("ITR:5"), std::string::npos) << "iteration must be logged";
    EXPECT_NE(contents.find("BOND"), std::string::npos) << "event type must be logged";
    // Format: ITR:<iter>,BOND,<name1>,<idx1>,<iface1>,<name2>,<idx2>,<iface2>
    EXPECT_NE(contents.find("ITR:5,BOND,IL,1,0,A,0,0"), std::string::npos)
        << "the whole record should follow the documented layout";
}

// -----------------------------------------------------------------------------
// Test 3: 3D -> 2D transition using a POINT binder.
//
// The protein complex starts in solution (OnSurface == false) so
// transitionToSurface becomes true and the binder is translated down to the
// membrane.  Because the protein template is flagged isPoint, every angular
// rotation is skipped, which is exactly the "no orientation" branch of the
// routine.
// -----------------------------------------------------------------------------
void test_ailbox_3d_to_2d_point_binder()
{
    std::cerr << "\n[TEST] test_ailbox_3d_to_2d_point_binder\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_box.cpp\n"
              << "  Function:    associate_implicitlipid_box()\n"
              << "  Scenario:    point protein in solution binds the implicit lipid\n"
              << "               (isOnMembrane == false, transitionToSurface == true)\n"
              << "  Checks:      the binder is translated down to sigma above the\n"
              << "               membrane, association is not cancelled and the\n"
              << "               complex becomes membrane bound.\n";

    ailbox_init_rng();

    AilSystem sys;
    const Coord protStart { 0.0, 0.0, 0.0 }; // middle of the box, in solution
    // proteinComplexOnSurface = false -> 3D->2D path
    // proteinIsPoint          = true  -> skip theta/phi/omega rotations
    ailbox_build_system(sys, /*onSurface*/ false, /*isPoint*/ true, /*ilFirst*/ false, protStart);

    const long eventsBefore = ailbox_total_events(sys.counterArrays);
    std::ofstream noLogFile;

    std::cerr << "  protein COM before call = (" << protStart.x << ", " << protStart.y
              << ", " << protStart.z << ")\n";
    std::cerr << "  Calling associate_implicitlipid_box()...\n";
    associate_implicitlipid_box(/*iter*/ 2, /*ifaceIndex1*/ 0, /*ifaceIndex2*/ 0,
        sys.moleculeList[kAilProtMolIdx], sys.moleculeList[kAilIlMolIdx],
        sys.complexList[0], sys.complexList[1], sys.params, sys.rxn, sys.moleculeList,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.complexList,
        sys.membrane, sys.forwardRxns, sys.backRxns, noLogFile);

    const Molecule& protein = sys.moleculeList[kAilProtMolIdx];
    const Molecule& lipid = sys.moleculeList[kAilIlMolIdx];
    const Complex& protComplex = sys.complexList[0];

    // --- the move must be accepted (none of the rejection counters fire) --
    std::cerr << "  cancel counters: span=" << sys.counterArrays.nCancelSpanBox
              << " overlapSystem=" << sys.counterArrays.nCancelOverlapSystem
              << " displace3Dto2D=" << sys.counterArrays.nCancelDisplace3Dto2D << '\n';
    EXPECT_EQ(sys.counterArrays.nCancelSpanBox, 0) << "small complexes cannot span the box";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapSystem, 0)
        << "overlap checking is disabled for both templates";
    EXPECT_EQ(sys.counterArrays.nCancelDisplace3Dto2D, 0)
        << "scaleMaxDisplace is huge, so displacement cannot reject the move";
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "association should be accepted";
    EXPECT_EQ(ailbox_total_events(sys.counterArrays) - eventsBefore, 1);

    // --- the binder ends up at the membrane ------------------------------
    std::cerr << "  protein COM after call  = (" << protein.comCoord.x << ", "
              << protein.comCoord.y << ", " << protein.comCoord.z << ")\n";
    EXPECT_LT(protein.comCoord.z, protStart.z - 40.0)
        << "the binder must have been dragged down towards the membrane";
    EXPECT_GE(protein.comCoord.z, -kAilBoxLen / 2.0 - 1e-6)
        << "the binder may not end up below the membrane";
    EXPECT_LE(protein.comCoord.z, -kAilBoxLen / 2.0 + kAilBindRadius + 1e-6)
        << "the binder should sit no higher than sigma above the membrane";

    // Interface separation should be of the order of the binding radius
    // (the implicit lipid was placed directly beneath the binder).
    const double dx = protein.interfaceList[0].coord.x - lipid.interfaceList[0].coord.x;
    const double dy = protein.interfaceList[0].coord.y - lipid.interfaceList[0].coord.y;
    const double dz = protein.interfaceList[0].coord.z - lipid.interfaceList[0].coord.z;
    const double sep = std::sqrt(dx * dx + dy * dy + dz * dz);
    std::cerr << "  interface separation after association = " << sep
              << " (bindRadius = " << kAilBindRadius << ")\n";
    EXPECT_LT(sep, kAilBindRadius + 0.2) << "interfaces should be at (or inside) sigma";

    // --- bookkeeping identical to the 2D case ----------------------------
    EXPECT_TRUE(protein.interfaceList[0].isBound);
    EXPECT_EQ(protein.interfaceList[0].index, kAilProdAbsIface);
    EXPECT_EQ(protein.linksToSurface, 1);
    EXPECT_TRUE(protComplex.OnSurface) << "complex transitioned onto the membrane";
    EXPECT_DOUBLE_EQ(protComplex.D.z, 0.0);
    EXPECT_EQ(protComplex.linksToSurface, 1);
    EXPECT_EQ(protComplex.iLipidIndex, kAilIlMolIdx);
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], kAilStartFreeLipids - 1);
    EXPECT_TRUE(protein.tmpICoords.empty()) << "tmp coords must be released";
    EXPECT_TRUE(lipid.tmpICoords.empty()) << "tmp coords must be released";
}

// -----------------------------------------------------------------------------
// Test 4: irreversible and unobserved reaction.
//
// Confirms the two optional branches at the end of the routine:
//   * currRxn.isReversible == false -> conjBackRxn is left at its default (-1)
//   * currRxn.isObserved  == false -> the observable map is untouched
// while the rest of the bookkeeping still happens.
// -----------------------------------------------------------------------------
void test_ailbox_irreversible_and_unobserved()
{
    std::cerr << "\n[TEST] test_ailbox_irreversible_and_unobserved\n"
              << "  Source file: src/reactions/associate_ImplicitLipid_box.cpp\n"
              << "  Function:    associate_implicitlipid_box()\n"
              << "  Scenario:    irreversible, unobserved 2D binding to the IL.\n"
              << "  Checks:      conjBackRxn stays -1, observables untouched, but\n"
              << "               copy numbers / bond lists are still updated.\n";

    ailbox_init_rng();

    AilSystem sys;
    const Coord protStart { -10.0, 10.0, -30.0 };
    ailbox_build_system(sys, /*onSurface*/ true, /*isPoint*/ false, /*ilFirst*/ false, protStart);

    // Flip the two optional flags off.
    sys.rxn.isReversible = false;
    sys.rxn.conjBackRxnIndex = -1;
    sys.rxn.isObserved = false;

    std::ofstream noLogFile;

    std::cerr << "  Calling associate_implicitlipid_box()...\n";
    associate_implicitlipid_box(/*iter*/ 3, /*ifaceIndex1*/ 0, /*ifaceIndex2*/ 0,
        sys.moleculeList[kAilProtMolIdx], sys.moleculeList[kAilIlMolIdx],
        sys.complexList[0], sys.complexList[1], sys.params, sys.rxn, sys.moleculeList,
        sys.molTemplateList, sys.observablesList, sys.counterArrays, sys.complexList,
        sys.membrane, sys.forwardRxns, sys.backRxns, noLogFile);

    const Molecule& protein = sys.moleculeList[kAilProtMolIdx];

    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1) << "association should still succeed";
    EXPECT_TRUE(protein.interfaceList[0].isBound);
    EXPECT_EQ(protein.interfaceList[0].interaction.conjBackRxn, -1)
        << "irreversible reactions must not store a back-reaction index";
    EXPECT_EQ(sys.observablesList["ILbond"], 0)
        << "an unobserved reaction must not touch the observable map";

    // The rest of the bookkeeping is unaffected by those two flags.
    EXPECT_EQ(protein.interfaceList[0].index, kAilProdAbsIface);
    EXPECT_EQ(protein.linksToSurface, 1);
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilProtAbsIface], 9);
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kAilProdAbsIface], 1);
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], kAilStartFreeLipids - 1);
    EXPECT_TRUE(protein.freelist.empty());
    EXPECT_TRUE(protein.tmpICoords.empty());
    EXPECT_TRUE(sys.moleculeList[kAilIlMolIdx].tmpICoords.empty());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario is reported separately, and because the
// helpers use EXPECT_* (not ASSERT_*, except where continuing would be
// meaningless) every test runs to completion even when one check fails.
// -----------------------------------------------------------------------------
TEST(AssociateImplicitLipidBox, TwoDBindingMol2IsImplicitLipid) { test_ailbox_2d_mol2_is_implicit_lipid(); }
TEST(AssociateImplicitLipidBox, TwoDBindingMol1IsImplicitLipid) { test_ailbox_2d_mol1_is_implicit_lipid(); }
TEST(AssociateImplicitLipidBox, ThreeDToTwoDPointBinder) { test_ailbox_3d_to_2d_point_binder(); }
TEST(AssociateImplicitLipidBox, IrreversibleAndUnobserved) { test_ailbox_irreversible_and_unobserved(); }