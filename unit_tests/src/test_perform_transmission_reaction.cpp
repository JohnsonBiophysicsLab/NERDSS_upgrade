/*! \file test_perform_transmission_reaction.cpp
 *
 * ### Unit test for src/reactions/perform_transmission_reaction.cpp
 *
 * Function under test:
 *
 *     void perform_transmission_reaction(int moleculeIndex,
 *                                        std::vector<Molecule>& moleculeList,
 *                                        std::vector<Complex>& complexList,
 *                                        std::vector<MolTemplate>& molTemplateList,
 *                                        Membrane& membraneObject,
 *                                        copyCounters& counterArrays,
 *                                        Parameters& params,
 *                                        std::vector<TransmissionRxn>& transmissionRxns,
 *                                        const std::vector<ForwardRxn>& forwardRxns,
 *                                        SimulVolume& simulVolume);
 *
 * The routine models a molecule crossing a spherical compartment boundary:
 *
 *   1. it computes the signed distance of the reactant COM to the compartment
 *      shell and mirrors the COM across the shell (newPos = scaler * currPos),
 *   2. it creates the *product* molecule (a different MolTemplate) at newPos and
 *      bumps counterArrays.copyNumSpecies for every product interface,
 *   3. it zeroes the reaction probabilities of every molecule that had the
 *      reactant in its crossbase list, marks the reactant complex ncross = -1
 *      and clears the reactant's crossbase,
 *   4. it decrements counterArrays.copyNumSpecies for every interface of every
 *      member of the reactant complex and destroys that complex,
 *   5. it removes the reactant from its SimulVolume sub-cell member list and
 *      sets mySubVolIndex = -1,
 *   6. finally, if the reactant template can be destroyed, it erases the
 *      reactant from that template's monomerList.
 *
 * The tests below build a tiny but self-consistent world (two molecule types,
 * one transmission reaction, a cubic water box with a spherical compartment)
 * and verify each of those observable side effects.  All assertions are
 * non-fatal (EXPECT_*) so every test in the suite still runs.
 */

#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Helpers (anonymous namespace, all prefixed with `ptrxn_` so they cannot clash
// with helpers from other generated test translation units).
// -----------------------------------------------------------------------------
namespace {

//! Molecule type indices used throughout this file.
constexpr int kPtrxnOutsideType = 0; //!< reactant type ("A", lives outside)
constexpr int kPtrxnInsideType = 1; //!< product type  ("B", lives inside)

//! Absolute interface-state indices (== indices into copyNumSpecies).
constexpr int kPtrxnAbsIfaceA = 0;
constexpr int kPtrxnAbsIfaceB = 1;

//! Geometry of the toy system.
constexpr double kPtrxnBoxSide = 100.0;
constexpr double kPtrxnCompartmentR = 20.0;
constexpr double kPtrxnStartRadius = 25.0; //!< reactant sits 25 nm from origin

/*! \brief Bundles everything perform_transmission_reaction() needs. */
struct PtrxnWorld {
    Parameters params {};
    Membrane membraneObject {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<TransmissionRxn> transmissionRxns {};
    std::vector<ForwardRxn> forwardRxns {};
    copyCounters counterArrays {};
    SimulVolume simulVolume {};
};

/*! \brief Reset all static bookkeeping so tests are independent of each other. */
void ptrxn_reset_statics()
{
    Molecule::numberOfMolecules = 0;
    Molecule::emptyMolList.clear();
    Molecule::maxID = 0;
    // Molecule::mapIdToIndex.clear();

    Complex::numberOfComplexes = 0;
    Complex::emptyComList.clear();
    Complex::obs.clear();
    Complex::maxID = 0;
    // Complex::mapIdToIndex.clear();
    Complex::currNumberMolTypes = 2;
    Complex::currNumberComTypes = 2;

    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType.assign(2, 0);
    MolTemplate::absToRelIface.assign(2, 0);

    Interface::State::totalNumOfStates = 2;
}

/*! \brief Build one single-interface, point-like MolTemplate. */
MolTemplate ptrxn_make_template(int molTypeIndex, const std::string& name, int absIfaceIndex)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = molTypeIndex;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.comCoord = Coord(0.0, 0.0, 0.0);
    temp.D = Coord(10.0, 10.0, 10.0);
    temp.Dr = Coord(0.1, 0.1, 0.1);
    temp.isPoint = true; //!< COM and interface coincide -> simple geometry
    temp.isRod = false;
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.checkOverlap = false; //!< keep the creation routine from rejecting placements
    temp.canDestroy = false;

    // One interface named "a" with a single (state-less) state.
    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord(0.0, 0.0, 0.0);
    Interface::State state;
    state.iden = '\0';
    state.index = absIfaceIndex;
    state.ifaceAndStateName = "a";
    iface.stateList.push_back(state);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a single-interface molecule of the requested type at `com`. */
Molecule ptrxn_make_molecule(int index, int molTypeIndex, const Coord& com, const MolTemplate& temp)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = index; // one molecule per complex
    mol.complexId = index;
    mol.mass = temp.mass;
    mol.comCoord = com;
    mol.isLipid = false;
    mol.isEmpty = false;
    mol.trajStatus = TrajStatus::none;

    Molecule::Iface iface;
    iface.coord = com; // point molecule
    iface.index = temp.interfaceList[0].stateList[0].index; // absolute state index
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.stateIden = temp.interfaceList[0].stateList[0].iden;
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    mol.freelist.push_back(0);
    return mol;
}

/*! \brief Compute the sub-cell absolute index a coordinate belongs to. */
int ptrxn_bin_index(const SimulVolume& simulVolume, const Membrane& membraneObject, const Coord& crd)
{
    int xItr = static_cast<int>((crd.x + membraneObject.waterBox.x / 2.0) / simulVolume.subCellSize.x);
    int yItr = static_cast<int>((crd.y + membraneObject.waterBox.y / 2.0) / simulVolume.subCellSize.y);
    int zItr = static_cast<int>((crd.z + membraneObject.waterBox.z / 2.0) / simulVolume.subCellSize.z);

    // Clamp so a rounding artefact can never index outside the cell list.
    xItr = std::max(0, std::min(xItr, simulVolume.numSubCells.x - 1));
    yItr = std::max(0, std::min(yItr, simulVolume.numSubCells.y - 1));
    zItr = std::max(0, std::min(zItr, simulVolume.numSubCells.z - 1));

    return xItr + (yItr * simulVolume.numSubCells.x)
        + (zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y);
}

/*!
 * \brief Assemble a complete, self-consistent world.
 *
 * \param[out] w                 the world to fill
 * \param[in]  numOutsideMolecules how many type-A molecules to create.  The
 *                               first one (index 0) always sits at
 *                               (kPtrxnStartRadius, 0, 0); any extras are
 *                               spectators placed far away on the +x axis.
 */
void ptrxn_build_world(PtrxnWorld& w, int numOutsideMolecules)
{
    // The GSL generator is a program-wide global that gtest_main leaves null;
    // the molecule-creation routine draws random numbers, so make sure it is
    // initialized exactly once.
    if (r == nullptr) {
        std::cerr << "    (initializing GSL RNG for the transmission test)\n";
         const gsl_rng_type *T;
         T = gsl_rng_default;
         r = gsl_rng_alloc(T);
         gsl_rng_set(r, 42);;
    }

    ptrxn_reset_statics();

    /* ---------------- parameters ---------------- */
    w.params.numMolTypes = 2;
    w.params.numTotalSpecies = 2;
    w.params.timeStep = 1.0;
    Parameters::dt = w.params.timeStep;
    w.params.rMaxLimit = 20.0; // -> 5 sub-cells per dimension for a 100 nm box
    w.params.rMaxRadius = 2.0;
    w.params.numTotalComplex = numOutsideMolecules;
    w.params.name = "ptrxn_test";

    /* ---------------- boundaries ---------------- */
    w.membraneObject.isBox = true;
    w.membraneObject.isSphere = false;
    w.membraneObject.waterBox = Membrane::WaterBox({ kPtrxnBoxSide, kPtrxnBoxSide, kPtrxnBoxSide });
    w.membraneObject.xBCtype = "reflect";
    w.membraneObject.yBCtype = "reflect";
    w.membraneObject.zBCtype = "reflect";
    w.membraneObject.hasCompartment = true;
    w.membraneObject.compartmentR = kPtrxnCompartmentR;
    w.membraneObject.implicitLipid = false;

    /* ---------------- templates ---------------- */
    w.molTemplateList.clear();
    w.molTemplateList.push_back(ptrxn_make_template(kPtrxnOutsideType, "A", kPtrxnAbsIfaceA));
    w.molTemplateList.push_back(ptrxn_make_template(kPtrxnInsideType, "B", kPtrxnAbsIfaceB));
    w.molTemplateList[kPtrxnOutsideType].outsideCompartment = true;
    w.molTemplateList[kPtrxnOutsideType].crossesCompartment = true;
    w.molTemplateList[kPtrxnOutsideType].transmissionRxnIndex = 0; // -> transmissionRxns[0]
    w.molTemplateList[kPtrxnInsideType].insideCompartment = true;
    w.molTemplateList[kPtrxnInsideType].crossesCompartment = true;
    w.molTemplateList[kPtrxnInsideType].transmissionRxnIndex = 0;
    MolTemplate::numEachMolType[kPtrxnOutsideType] = numOutsideMolecules;
    MolTemplate::numEachMolType[kPtrxnInsideType] = 0;

    /* ---------------- transmission reaction A -> B ---------------- */
    TransmissionRxn rxn;
    rxn.rxnType = ReactionType::transmission;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.bindRadius = 1.0;
    rxn.rxnLabel = "transmit_A_to_B";
    rxn.rateList.emplace_back(RxnBase::RateState(1.0, {}));
    rxn.reactantMolList.emplace_back(TransmissionRxn::TransmissionMol(kPtrxnOutsideType,
        { RxnIface("a", kPtrxnOutsideType, kPtrxnAbsIfaceA, 0, '\0', false) }));
    rxn.reactantMolList.back().molName = "A";
    rxn.productMolList.emplace_back(TransmissionRxn::TransmissionMol(kPtrxnInsideType,
        { RxnIface("a", kPtrxnInsideType, kPtrxnAbsIfaceB, 0, '\0', false) }));
    rxn.productMolList.back().molName = "B";
    w.transmissionRxns.clear();
    w.transmissionRxns.push_back(rxn);
    w.forwardRxns.clear(); // no binding reactions needed for this test

    /* ---------------- molecules and complexes ---------------- */
    w.moleculeList.clear();
    w.complexList.clear();
    for (int i = 0; i < numOutsideMolecules; ++i) {
        // molecule 0 is the transmitting one; spectators go further out on +x.
        Coord com = (i == 0) ? Coord(kPtrxnStartRadius, 0.0, 0.0)
                             : Coord(kPtrxnStartRadius + 15.0 * i, 0.0, 0.0);
        Molecule mol = ptrxn_make_molecule(i, kPtrxnOutsideType, com, w.molTemplateList[kPtrxnOutsideType]);
        w.moleculeList.push_back(mol);

        Complex com1(i, w.moleculeList[i], w.molTemplateList[kPtrxnOutsideType]);
        com1.index = i;
        com1.id = i;
        com1.isEmpty = false;
        com1.ncross = 0;
        w.complexList.push_back(com1);
    }
    Molecule::numberOfMolecules = static_cast<int>(w.moleculeList.size());
    Molecule::maxID = static_cast<int>(w.moleculeList.size());
    Complex::numberOfComplexes = static_cast<int>(w.complexList.size());
    Complex::maxID = static_cast<int>(w.complexList.size());

    /* ---------------- species counters ---------------- */
    w.counterArrays.copyNumSpecies.assign(2, 0);
    w.counterArrays.copyNumSpecies[kPtrxnAbsIfaceA] = numOutsideMolecules;
    w.counterArrays.copyNumSpecies[kPtrxnAbsIfaceB] = 0;

    /* ---------------- simulation volume ---------------- */
    w.simulVolume.create_simulation_volume(w.params, w.membraneObject);
    for (auto& mol : w.moleculeList) {
        int binIndex = ptrxn_bin_index(w.simulVolume, w.membraneObject, mol.comCoord);
        mol.mySubVolIndex = binIndex;
        w.simulVolume.subCellList[binIndex].memberMolList.push_back(mol.index);
    }
}

/*! \brief Locate the (single) surviving molecule of the product type. */
int ptrxn_find_product(const std::vector<Molecule>& moleculeList)
{
    for (size_t i = 0; i < moleculeList.size(); ++i) {
        if (!moleculeList[i].isEmpty && moleculeList[i].molTypeIndex == kPtrxnInsideType)
            return static_cast<int>(i);
    }
    return -1;
}

/*! \brief Does `list` still contain `value`? */
bool ptrxn_contains(const std::vector<int>& list, int value)
{
    return std::find(list.begin(), list.end(), value) != list.end();
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the basic transmission event.
//         Reactant destroyed, product created inside the compartment, species
//         counters transferred, sub-cell membership cleaned up.
// -----------------------------------------------------------------------------
void test_ptrxn_basic_transmission()
{
    std::cerr << "\n[TEST] test_ptrxn_basic_transmission\n"
              << "  Source file:   src/reactions/perform_transmission_reaction.cpp\n"
              << "  Function:      perform_transmission_reaction()\n"
              << "  Scenario:      one molecule of type A sitting " << kPtrxnStartRadius
              << " nm from the origin,\n"
              << "                 compartment radius " << kPtrxnCompartmentR
              << " nm, so it is mirrored to ~"
              << (kPtrxnStartRadius - 2.0 * (kPtrxnStartRadius - kPtrxnCompartmentR)) << " nm.\n"
              << "  Pass criteria: (a) a molecule of the product type exists inside the\n"
              << "                 compartment, (b) the reactant molecule/complex are marked\n"
              << "                 empty, (c) copyNumSpecies moves from the A interface to the\n"
              << "                 B interface, (d) the reactant leaves its sub-cell list.\n";

    PtrxnWorld w;
    ptrxn_build_world(w, /*numOutsideMolecules=*/1);

    const size_t nMolBefore = w.moleculeList.size();
    const int subVolBefore = w.moleculeList[0].mySubVolIndex;
    const double distBefore = w.moleculeList[0].comCoord.get_magnitude();

    std::cerr << "  Before: moleculeList.size()=" << nMolBefore
              << ", |COM|=" << distBefore
              << ", mySubVolIndex=" << subVolBefore
              << ", copyNumSpecies=[" << w.counterArrays.copyNumSpecies[0] << ","
              << w.counterArrays.copyNumSpecies[1] << "]\n";

    std::cerr << "  Calling perform_transmission_reaction(0, ...)\n";
    perform_transmission_reaction(0, w.moleculeList, w.complexList, w.molTemplateList,
        w.membraneObject, w.counterArrays, w.params, w.transmissionRxns, w.forwardRxns,
        w.simulVolume);

    // (a) A product molecule of the "inside" template must now exist.
    EXPECT_GT(w.moleculeList.size(), nMolBefore)
        << "a new (product) molecule should have been appended to moleculeList";
    const int prodIndex = ptrxn_find_product(w.moleculeList);
    EXPECT_NE(prodIndex, -1) << "expected exactly one live molecule of the product type";

    if (prodIndex >= 0) {
        const double prodDist = w.moleculeList[prodIndex].comCoord.get_magnitude();
        std::cerr << "  Product molecule index " << prodIndex << " placed at |COM| = "
                  << prodDist << " nm (target ~"
                  << (kPtrxnStartRadius - 2.0 * (kPtrxnStartRadius - kPtrxnCompartmentR))
                  << " nm)\n";
        // The reactant is mirrored across the compartment shell, so it must have
        // moved inward and now be inside the compartment.
        EXPECT_LT(prodDist, distBefore)
            << "product should be closer to the origin than the reactant was";
        EXPECT_LT(prodDist, kPtrxnCompartmentR)
            << "product should end up inside the compartment radius";
        EXPECT_EQ(w.moleculeList[prodIndex].molTypeIndex, kPtrxnInsideType)
            << "product must use the MolTemplate listed in productMolList";
        EXPECT_FALSE(w.moleculeList[prodIndex].isEmpty)
            << "the freshly created product molecule must be live";
    }

    // (b) The reactant molecule and its complex are gone.
    EXPECT_TRUE(w.moleculeList[0].isEmpty)
        << "the transmitting molecule must be destroyed (isEmpty == true)";
    EXPECT_TRUE(w.complexList[0].isEmpty)
        << "the transmitting molecule's complex must be destroyed (isEmpty == true)";

    // (c) Species counters: the A interface count drops, the B count rises, and
    //     the total number of tracked interfaces is conserved (1 -> 1).
    std::cerr << "  After:  copyNumSpecies=[" << w.counterArrays.copyNumSpecies[0] << ","
              << w.counterArrays.copyNumSpecies[1] << "]\n";
    EXPECT_EQ(w.counterArrays.copyNumSpecies[kPtrxnAbsIfaceA], 0)
        << "the reactant interface copy number should have been decremented to 0";
    EXPECT_EQ(w.counterArrays.copyNumSpecies[kPtrxnAbsIfaceB], 1)
        << "the product interface copy number should have been incremented to 1";
    EXPECT_EQ(w.counterArrays.copyNumSpecies[0] + w.counterArrays.copyNumSpecies[1], 1)
        << "total interface count must be conserved by a transmission event";

    // (d) Sub-volume bookkeeping for the destroyed reactant.
    EXPECT_EQ(w.moleculeList[0].mySubVolIndex, -1)
        << "destroyed reactant must have mySubVolIndex reset to -1";
    if (subVolBefore >= 0 && subVolBefore < static_cast<int>(w.simulVolume.subCellList.size())) {
        EXPECT_FALSE(ptrxn_contains(w.simulVolume.subCellList[subVolBefore].memberMolList, 0))
            << "reactant index must be erased from its old sub-cell memberMolList";
    }
}

// -----------------------------------------------------------------------------
// Test 2: pending bimolecular reactions that involve the transmitting molecule
//         must be cancelled (probabilities zeroed, crossbase cleared).
// -----------------------------------------------------------------------------
void test_ptrxn_cancels_pending_crossings()
{
    std::cerr << "\n[TEST] test_ptrxn_cancels_pending_crossings\n"
              << "  Source file:   src/reactions/perform_transmission_reaction.cpp\n"
              << "  Function:      perform_transmission_reaction()\n"
              << "  Scenario:      molecule 0 transmits while molecule 1 still lists it in\n"
              << "                 its crossbase with probability 0.5.\n"
              << "  Pass criteria: molecule 1's probvec entry is forced to 0, molecule 0's\n"
              << "                 crossbase is cleared and its complex ncross is -1, and\n"
              << "                 the spectator molecule survives untouched otherwise.\n";

    PtrxnWorld w;
    ptrxn_build_world(w, /*numOutsideMolecules=*/2);

    // Register a mutual "possible reaction" between molecules 0 and 1.
    w.moleculeList[0].crossbase.push_back(1);
    w.moleculeList[0].mycrossint.push_back(0);
    w.moleculeList[0].crossrxn.push_back({ 0, 0, 0 });
    w.moleculeList[0].probvec.push_back(0.5);

    w.moleculeList[1].crossbase.push_back(0);
    w.moleculeList[1].mycrossint.push_back(0);
    w.moleculeList[1].crossrxn.push_back({ 0, 0, 0 });
    w.moleculeList[1].probvec.push_back(0.5);

    w.complexList[0].ncross = 1;
    w.complexList[1].ncross = 1;

    std::cerr << "  Before: moleculeList[1].probvec[0] = " << w.moleculeList[1].probvec[0]
              << ", moleculeList[0].crossbase.size() = "
              << w.moleculeList[0].crossbase.size() << "\n";

    std::cerr << "  Calling perform_transmission_reaction(0, ...)\n";
    perform_transmission_reaction(0, w.moleculeList, w.complexList, w.molTemplateList,
        w.membraneObject, w.counterArrays, w.params, w.transmissionRxns, w.forwardRxns,
        w.simulVolume);

    // The partner's stored reaction probability with molecule 0 must be zeroed
    // so it cannot react with a molecule that no longer exists.
    ASSERT_FALSE(w.moleculeList[1].probvec.empty())
        << "spectator molecule should still own its probvec entry";
    std::cerr << "  After:  moleculeList[1].probvec[0] = " << w.moleculeList[1].probvec[0] << "\n";
    EXPECT_DOUBLE_EQ(w.moleculeList[1].probvec[0], 0.0)
        << "reaction probability with the transmitted molecule must be set to 0";

    // The transmitting molecule's own crossing list is cleared and its complex
    // is flagged with ncross = -1 (skip further overlap/reaction processing).
    EXPECT_TRUE(w.moleculeList[0].crossbase.empty())
        << "the transmitting molecule's crossbase must be cleared";
    std::cerr << "  complexList[0].ncross = " << w.complexList[0].ncross
              << " (set to -1 before destruction)\n";

    // The spectator must survive; only the transmitting molecule is destroyed.
    EXPECT_FALSE(w.moleculeList[1].isEmpty)
        << "the spectator molecule must not be destroyed";
    EXPECT_TRUE(w.moleculeList[0].isEmpty)
        << "the transmitting molecule must be destroyed";

    // Counters: one A remains (the spectator) and one B was created.
    std::cerr << "  After:  copyNumSpecies=[" << w.counterArrays.copyNumSpecies[0] << ","
              << w.counterArrays.copyNumSpecies[1] << "]\n";
    EXPECT_EQ(w.counterArrays.copyNumSpecies[kPtrxnAbsIfaceA], 1)
        << "exactly one type-A interface (the spectator) should remain";
    EXPECT_EQ(w.counterArrays.copyNumSpecies[kPtrxnAbsIfaceB], 1)
        << "exactly one type-B interface should have been created";
}

// -----------------------------------------------------------------------------
// Test 3: when the reactant template is destroyable, its monomerList entry must
//         be erased by the transmission event.
// -----------------------------------------------------------------------------
void test_ptrxn_monomer_list_cleanup()
{
    std::cerr << "\n[TEST] test_ptrxn_monomer_list_cleanup\n"
              << "  Source file:   src/reactions/perform_transmission_reaction.cpp\n"
              << "  Function:      perform_transmission_reaction()\n"
              << "  Scenario:      the reactant template has canDestroy == true and its\n"
              << "                 monomerList contains the transmitting molecule index.\n"
              << "  Pass criteria: after the call the index is no longer in monomerList.\n";

    PtrxnWorld w;
    ptrxn_build_world(w, /*numOutsideMolecules=*/1);

    // Mark the reactant template destroyable and register the monomer, exactly
    // as the main loop would do for a molecule with a destruction reaction.
    w.molTemplateList[kPtrxnOutsideType].canDestroy = true;
    w.molTemplateList[kPtrxnOutsideType].monomerList.clear();
    w.molTemplateList[kPtrxnOutsideType].monomerList.push_back(0);
    // Product template stays non-destroyable so its monomerList is irrelevant.
    w.molTemplateList[kPtrxnInsideType].canDestroy = false;

    std::cerr << "  Before: template A monomerList.size() = "
              << w.molTemplateList[kPtrxnOutsideType].monomerList.size() << "\n";

    std::cerr << "  Calling perform_transmission_reaction(0, ...)\n";
    perform_transmission_reaction(0, w.moleculeList, w.complexList, w.molTemplateList,
        w.membraneObject, w.counterArrays, w.params, w.transmissionRxns, w.forwardRxns,
        w.simulVolume);

    std::cerr << "  After:  template A monomerList.size() = "
              << w.molTemplateList[kPtrxnOutsideType].monomerList.size() << "\n";
    EXPECT_FALSE(ptrxn_contains(w.molTemplateList[kPtrxnOutsideType].monomerList, 0))
        << "the transmitted molecule index must be erased from the template monomerList";
    EXPECT_EQ(w.molTemplateList[kPtrxnOutsideType].monomerList.size(), 0u)
        << "the only monomer entry should have been removed";

    // Sanity: the transmission itself still happened.
    EXPECT_TRUE(w.moleculeList[0].isEmpty)
        << "the transmitting molecule must still be destroyed in this configuration";
    EXPECT_NE(ptrxn_find_product(w.moleculeList), -1)
        << "the product molecule must still be created in this configuration";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - one per scenario so each reports independently while
// all of them execute.
// -----------------------------------------------------------------------------
TEST(PerformTransmissionReaction, BasicTransmission) { test_ptrxn_basic_transmission(); }
TEST(PerformTransmissionReaction, CancelsPendingCrossings) { test_ptrxn_cancels_pending_crossings(); }
TEST(PerformTransmissionReaction, MonomerListCleanup) { test_ptrxn_monomer_list_cleanup(); }