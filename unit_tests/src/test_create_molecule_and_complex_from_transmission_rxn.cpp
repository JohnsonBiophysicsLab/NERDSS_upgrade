/*! \file test_create_molecule_and_complex_from_transmission_rxn.cpp
 *
 * ### Unit test for
 *     src/reactions/create_molecule_and_complex_from_transmission_rxn.cpp
 *
 * Function under test:
 *
 *     void create_molecule_and_complex_from_transmission_rxn(
 *              int parentMolIndex, int& newMolIndex, int& newComIndex,
 *              MolTemplate& createdMolTemp, Parameters& params,
 *              const TransmissionRxn& currRxn, SimulVolume& simulVolume,
 *              std::vector<Molecule>& moleculeList,
 *              std::vector<Complex>& complexList,
 *              std::vector<MolTemplate>& molTemplateList,
 *              const std::vector<ForwardRxn>& forwardRxns,
 *              const Membrane& membraneObject, const Coord& newPos)
 *
 * The routine is a "book-keeping" function.  Its observable behaviour is:
 *
 *   1. It picks a slot for the new Molecule.  If Molecule::emptyMolList holds
 *      indices of recycled (isEmpty == true) molecules it re-uses the last one,
 *      otherwise it appends a brand new element to moleculeList.
 *   2. The same logic is applied to Complex::emptyComList / complexList.
 *   3. Stale entries (indices of molecules/complexes that are *not* actually
 *      empty any more) sitting on top of the recycle lists are popped off
 *      before a slot is accepted.
 *   4. The new Molecule is initialised (via
 *      initialize_molecule_after_transmission_reaction) and re-sampled while
 *      moleculeOverlaps() reports an overlap.
 *   5. Book-keeping is finished: myComIndex, TrajStatus::propagated and
 *      isDissociated are set on the Molecule, a fresh Complex is constructed in
 *      place, its TrajStatus is set, Complex::numberOfComplexes is incremented
 *      and - when the MolTemplate is destroyable - the new molecule index is
 *      pushed onto MolTemplate::monomerList.
 *
 * The tests below build a very small but fully initialised world (one molecule
 * type with one interface, one parent molecule/complex, a real SimulVolume and
 * an empty ForwardRxn list so that the overlap check can never succeed) and then
 * verify each of the points above.  Every test prints what it is doing to
 * stderr so that the console log explains itself.
 */

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

// `r` (the global GSL generator) and `totMatches` are defined in
// unit_tests/src/gtest_main.cpp; they are only declared here through the
// project headers above.

namespace {

// -----------------------------------------------------------------------------
// RNG helper.  The creation routine can call rand_gsl()/GaussV() when it has to
// re-sample an overlapping position, therefore the global generator must exist.
// A fixed seed keeps the test deterministic.
// -----------------------------------------------------------------------------
void cmacftr_seed_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, 42);
}

// -----------------------------------------------------------------------------
// A minimal, but completely initialised, simulation "world".
//
// Everything the function under test touches is filled in explicitly: an
// under-initialised Molecule / Complex / MolTemplate / Parameters would make the
// whole gtest binary crash rather than fail a single assertion.
// -----------------------------------------------------------------------------
struct CmacftrWorld {
    Parameters params {};
    Membrane membrane {};
    SimulVolume simulVolume {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {}; //!< deliberately empty: no pair of
                                            //!< molecules can ever "overlap"
    TransmissionRxn rxn {};

    //! Where the transmitted molecule is supposed to appear (on the compartment
    //! surface, which sits at radius 15 nm around the origin).
    Coord newPos { 15.0, 0.0, 0.0 };

    /*! \brief Fill the world.
     *  \param[in] canDestroy value given to MolTemplate::canDestroy - controls
     *                        whether the new molecule is added to monomerList.
     */
    void build(bool canDestroy);
};

void CmacftrWorld::build(bool canDestroy)
{
    cmacftr_seed_rng();

    // ---- reset the process-wide static book-keeping ------------------------
    // These statics are shared with every other test in the suite, so they are
    // reset here to make this test self-contained and deterministic.
    Molecule::emptyMolList.clear();
    Complex::emptyComList.clear();
    MolTemplate::numMolTypes = 1;
    MolTemplate::numEachMolType = std::vector<int>(1, 1);
    Molecule::numberOfMolecules = 1;
    Complex::numberOfComplexes = 1;

    // ---- simulation parameters --------------------------------------------
    params.numMolTypes = 1;
    params.numTotalSpecies = 1;
    params.numTotalComplex = 1;
    params.nItr = 10;
    params.timeStep = 0.1;
    params.overlapSepLimit = 0.1;
    params.scaleMaxDisplace = 100.0;
    params.rMaxLimit = 25.0; // gives a 4 x 4 x 4 sub-volume grid in a 100 nm box
    params.rMaxRadius = 2.0;
    params.name = "cmacftr_test";

    // ---- boundary ----------------------------------------------------------
    membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    membrane.isBox = true;
    membrane.isSphere = false;
    membrane.implicitLipid = false;
    membrane.hasCompartment = true;
    membrane.compartmentR = 15.0;
    membrane.xBCtype = "reflect";
    membrane.yBCtype = "reflect";
    membrane.zBCtype = "reflect";

    // ---- one molecule type, one interface ---------------------------------
    MolTemplate temp {};
    temp.molName = "A";
    temp.molTypeIndex = 0;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord { 10.0, 10.0, 10.0 };
    temp.Dr = Coord { 0.1, 0.1, 0.1 };
    temp.comCoord = Coord { 0.0, 0.0, 0.0 };
    temp.isLipid = false;
    temp.isPoint = false;
    temp.isRod = false;
    temp.isImplicitLipid = false;
    temp.isPromoter = false;
    temp.checkOverlap = false;
    temp.canDestroy = canDestroy;
    temp.crossesCompartment = true;
    temp.insideCompartment = true;
    temp.transmissionRxnIndex = 0;

    Interface iface {};
    iface.name = "a";
    iface.index = 0;
    iface.iCoord = Coord { 1.0, 0.0, 0.0 };
    iface.stateList.emplace_back(std::string("a"), 0); // single, unnamed state
    temp.interfaceList.push_back(iface);

    molTemplateList.push_back(temp);

    // ---- the parent molecule (the one that is "transmitted") ---------------
    // Kept far away from newPos so that no overlap can be detected no matter
    // which distance criterion moleculeOverlaps() applies.
    Molecule parent {};
    parent.index = 0;
    parent.id = 0;
    parent.molTypeIndex = 0;
    parent.myComIndex = 0;
    parent.mass = 1.0;
    parent.isEmpty = false;
    parent.isLipid = false;
    parent.comCoord = Coord { -30.0, -30.0, -30.0 };

    Molecule::Iface parentIface {};
    parentIface.coord = parent.comCoord + Coord { 1.0, 0.0, 0.0 };
    parentIface.index = 0;
    parentIface.relIndex = 0;
    parentIface.stateIndex = 0;
    parentIface.stateIden = '\0';
    parentIface.molTypeIndex = 0;
    parentIface.isBound = false;
    parent.interfaceList.push_back(parentIface);
    parent.freelist.push_back(0);

    moleculeList.push_back(parent);

    // ---- the parent complex -------------------------------------------------
    complexList.emplace_back(0, moleculeList[0], molTemplateList[0]);

    // ---- the transmission reaction -----------------------------------------
    rxn.rxnType = ReactionType::transmission;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.bindRadius = 1.0;
    rxn.rxnLabel = "transmit";
    rxn.isOnMem = false;

    TransmissionRxn::TransmissionMol tMol {};
    tMol.molTypeIndex = 0;
    tMol.molName = "A";
    tMol.interfaceList.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.reactantMolList.push_back(tMol);
    rxn.productMolList.push_back(tMol);

    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.productListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.intReactantList.push_back(0);
    rxn.intProductList.push_back(0);
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = 1.0;
    rxn.rateList.back().prob = 0.5;

    // ---- the sub-volume grid -----------------------------------------------
    simulVolume.create_simulation_volume(params, membrane);
    // simItr == 1 selects the cheap "just assign bins" branch (no boundary
    // sanity checks / no exit() calls).
    simulVolume.update_memberMolLists(params, moleculeList, complexList,
        molTemplateList, membrane, 1);
}

//! Thin wrapper so every test calls the function under test the same way.
void cmacftr_run(CmacftrWorld& w, int& newMolIndex, int& newComIndex)
{
    create_molecule_and_complex_from_transmission_rxn(
        /*parentMolIndex=*/0, newMolIndex, newComIndex, w.molTemplateList[0],
        w.params, w.rxn, w.simulVolume, w.moleculeList, w.complexList,
        w.molTemplateList, w.forwardRxns, w.membrane, w.newPos);
}

//! Non-fatal bounds check so that a bad index reports a failure instead of
//! taking the whole test binary down with an out-of-range access.
bool cmacftr_indices_valid(const CmacftrWorld& w, int newMolIndex, int newComIndex)
{
    if (newMolIndex < 0 || newMolIndex >= static_cast<int>(w.moleculeList.size())) {
        ADD_FAILURE() << "newMolIndex " << newMolIndex
                      << " is outside moleculeList (size " << w.moleculeList.size() << ")";
        return false;
    }
    if (newComIndex < 0 || newComIndex >= static_cast<int>(w.complexList.size())) {
        ADD_FAILURE() << "newComIndex " << newComIndex
                      << " is outside complexList (size " << w.complexList.size() << ")";
        return false;
    }
    return true;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: with empty recycle lists a brand new Molecule and Complex must be
//         appended to the back of the two containers and fully initialised.
// -----------------------------------------------------------------------------
void test_cmacftr_creates_new_molecule_and_complex()
{
    std::cerr << "\n[TEST] test_cmacftr_creates_new_molecule_and_complex\n"
              << "  Source file: create_molecule_and_complex_from_transmission_rxn.cpp\n"
              << "  Function:    create_molecule_and_complex_from_transmission_rxn\n"
              << "  Scenario:    Molecule::emptyMolList and Complex::emptyComList are\n"
              << "               empty, so new slots have to be appended.\n"
              << "  Pass:        indices == old container sizes, containers grew by one,\n"
              << "               molecule/complex book-keeping fields are set.\n";

    CmacftrWorld w;
    w.build(/*canDestroy=*/false);

    const std::size_t molsBefore = w.moleculeList.size();
    const std::size_t comsBefore = w.complexList.size();
    const int numComplexesBefore = Complex::numberOfComplexes;

    std::cerr << "  moleculeList.size() before = " << molsBefore
              << ", complexList.size() before = " << comsBefore
              << ", Complex::numberOfComplexes before = " << numComplexesBefore << '\n';

    int newMolIndex = -1;
    int newComIndex = -1;
    cmacftr_run(w, newMolIndex, newComIndex);

    std::cerr << "  -> newMolIndex = " << newMolIndex
              << ", newComIndex = " << newComIndex << '\n';

    // The new objects must have been appended at the back of the containers.
    EXPECT_EQ(newMolIndex, static_cast<int>(molsBefore))
        << "A fresh molecule slot must be the old back of moleculeList";
    EXPECT_EQ(newComIndex, static_cast<int>(comsBefore))
        << "A fresh complex slot must be the old back of complexList";
    EXPECT_EQ(w.moleculeList.size(), molsBefore + 1) << "moleculeList must grow by one";
    EXPECT_EQ(w.complexList.size(), comsBefore + 1) << "complexList must grow by one";

    if (!cmacftr_indices_valid(w, newMolIndex, newComIndex))
        return;

    const Molecule& newMol = w.moleculeList[newMolIndex];
    Complex& newCom = w.complexList[newComIndex];

    std::cerr << "  -> new molecule COM = (" << newMol.comCoord.x << ", "
              << newMol.comCoord.y << ", " << newMol.comCoord.z << ")\n";

    // Book-keeping performed by the function itself.
    EXPECT_EQ(newMol.myComIndex, newComIndex)
        << "The new molecule must point at the newly created complex";
    EXPECT_EQ(static_cast<int>(newMol.trajStatus), static_cast<int>(TrajStatus::propagated))
        << "The new molecule must be flagged as already propagated this step";
    EXPECT_TRUE(newMol.isDissociated)
        << "A transmitted molecule is flagged as dissociated for this step";
    EXPECT_FALSE(newMol.isEmpty) << "The new molecule slot must not be empty";
    EXPECT_EQ(newMol.molTypeIndex, w.molTemplateList[0].molTypeIndex)
        << "The new molecule must have the type of the created template";

    // Sanity of the produced coordinate: finite and inside the simulation box.
    EXPECT_FALSE(std::isnan(newMol.comCoord.x)) << "COM x must be a number";
    EXPECT_FALSE(std::isnan(newMol.comCoord.y)) << "COM y must be a number";
    EXPECT_FALSE(std::isnan(newMol.comCoord.z)) << "COM z must be a number";
    Coord comCopy { newMol.comCoord };
    EXPECT_FALSE(comCopy.isOutOfBox(w.membrane))
        << "The transmitted molecule must be placed inside the water box";

    // The complex constructed in place.
    EXPECT_EQ(newCom.index, newComIndex) << "Complex::index must equal its slot index";
    EXPECT_FALSE(newCom.isEmpty) << "The new complex slot must not be empty";
    EXPECT_EQ(static_cast<int>(newCom.trajStatus), static_cast<int>(TrajStatus::propagated))
        << "The new complex must be flagged as already propagated this step";
    EXPECT_EQ(newCom.memberList.size(), 1u) << "A newly transmitted complex holds one molecule";
    if (!newCom.memberList.empty()) {
        EXPECT_EQ(newCom.memberList[0], newMolIndex)
            << "The single member must be the newly created molecule";
    }
    EXPECT_EQ(Complex::numberOfComplexes, numComplexesBefore + 1)
        << "Complex::numberOfComplexes must be incremented exactly once";

    std::cerr << "  -> Complex::numberOfComplexes after = " << Complex::numberOfComplexes << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: recycled (empty) slots must be re-used instead of growing the arrays.
// -----------------------------------------------------------------------------
void test_cmacftr_reuses_empty_slots()
{
    std::cerr << "\n[TEST] test_cmacftr_reuses_empty_slots\n"
              << "  Source file: create_molecule_and_complex_from_transmission_rxn.cpp\n"
              << "  Function:    create_molecule_and_complex_from_transmission_rxn\n"
              << "  Scenario:    moleculeList[1] and complexList[1] are empty shells whose\n"
              << "               indices sit on Molecule::emptyMolList / Complex::emptyComList.\n"
              << "  Pass:        indices 1 are handed out, containers do NOT grow and the\n"
              << "               recycle lists are consumed.\n";

    CmacftrWorld w;
    w.build(/*canDestroy=*/false);

    // Append one empty molecule / complex shell and register them for recycling.
    Molecule emptyMol {};
    emptyMol.index = 1;
    emptyMol.isEmpty = true;
    w.moleculeList.push_back(emptyMol);

    Complex emptyCom {};
    emptyCom.index = 1;
    emptyCom.isEmpty = true;
    w.complexList.push_back(emptyCom);

    Molecule::emptyMolList.push_back(1);
    Complex::emptyComList.push_back(1);

    const std::size_t molsBefore = w.moleculeList.size(); // 2
    const std::size_t comsBefore = w.complexList.size(); // 2

    std::cerr << "  moleculeList.size() before = " << molsBefore
              << " (index 1 is an empty shell)\n";

    int newMolIndex = -1;
    int newComIndex = -1;
    cmacftr_run(w, newMolIndex, newComIndex);

    std::cerr << "  -> newMolIndex = " << newMolIndex
              << ", newComIndex = " << newComIndex << '\n';

    EXPECT_EQ(newMolIndex, 1) << "The recycled molecule slot 1 must be re-used";
    EXPECT_EQ(newComIndex, 1) << "The recycled complex slot 1 must be re-used";
    EXPECT_EQ(w.moleculeList.size(), molsBefore) << "moleculeList must not grow when recycling";
    EXPECT_EQ(w.complexList.size(), comsBefore) << "complexList must not grow when recycling";

    EXPECT_TRUE(Molecule::emptyMolList.empty())
        << "The consumed index must be popped off Molecule::emptyMolList";
    EXPECT_TRUE(Complex::emptyComList.empty())
        << "The consumed index must be popped off Complex::emptyComList";

    if (!cmacftr_indices_valid(w, newMolIndex, newComIndex))
        return;

    EXPECT_FALSE(w.moleculeList[newMolIndex].isEmpty)
        << "The recycled molecule slot must now be occupied";
    EXPECT_EQ(w.moleculeList[newMolIndex].myComIndex, newComIndex)
        << "The recycled molecule must point at the recycled complex";
    EXPECT_FALSE(w.complexList[newComIndex].isEmpty)
        << "The recycled complex slot must now be occupied";
    EXPECT_EQ(w.complexList[newComIndex].index, newComIndex)
        << "The recycled complex must carry its own slot index";
}

// -----------------------------------------------------------------------------
// Test 3: stale entries (indices of molecules/complexes which are no longer
//         empty) sitting on top of the recycle lists must be discarded.
// -----------------------------------------------------------------------------
void test_cmacftr_skips_stale_empty_entries()
{
    std::cerr << "\n[TEST] test_cmacftr_skips_stale_empty_entries\n"
              << "  Source file: create_molecule_and_complex_from_transmission_rxn.cpp\n"
              << "  Function:    create_molecule_and_complex_from_transmission_rxn\n"
              << "  Scenario:    recycle lists are {1, 2} where slot 2 is *occupied* (stale)\n"
              << "               and slot 1 is genuinely empty.\n"
              << "  Pass:        the stale index 2 is popped, index 1 is handed out and both\n"
              << "               recycle lists end up empty.\n";

    CmacftrWorld w;
    w.build(/*canDestroy=*/false);

    // index 1 : genuinely empty shells
    Molecule emptyMol {};
    emptyMol.index = 1;
    emptyMol.isEmpty = true;
    w.moleculeList.push_back(emptyMol);

    Complex emptyCom {};
    emptyCom.index = 1;
    emptyCom.isEmpty = true;
    w.complexList.push_back(emptyCom);

    // index 2 : occupied objects whose indices were left behind on the lists.
    Molecule occupiedMol { w.moleculeList[0] }; // copy of the fully valid parent
    occupiedMol.index = 2;
    occupiedMol.myComIndex = 2;
    occupiedMol.isEmpty = false;
    occupiedMol.comCoord = Coord { 40.0, 40.0, 40.0 };
    w.moleculeList.push_back(occupiedMol);

    Complex occupiedCom { w.complexList[0] }; // copy of the fully valid parent complex
    occupiedCom.index = 2;
    occupiedCom.isEmpty = false;
    occupiedCom.memberList = std::vector<int> { 2 };
    w.complexList.push_back(occupiedCom);

    // The stale index has to be on top (back) of the list.
    Molecule::emptyMolList = std::vector<int> { 1, 2 };
    Complex::emptyComList = std::vector<int> { 1, 2 };

    const std::size_t molsBefore = w.moleculeList.size(); // 3
    const std::size_t comsBefore = w.complexList.size(); // 3

    std::cerr << "  emptyMolList = {1, 2}; slot 2 is occupied -> must be discarded\n";

    int newMolIndex = -1;
    int newComIndex = -1;
    cmacftr_run(w, newMolIndex, newComIndex);

    std::cerr << "  -> newMolIndex = " << newMolIndex
              << ", newComIndex = " << newComIndex
              << ", emptyMolList.size() = " << Molecule::emptyMolList.size()
              << ", emptyComList.size() = " << Complex::emptyComList.size() << '\n';

    EXPECT_EQ(newMolIndex, 1) << "The stale index 2 must be skipped and slot 1 used";
    EXPECT_EQ(newComIndex, 1) << "The stale index 2 must be skipped and slot 1 used";
    EXPECT_EQ(w.moleculeList.size(), molsBefore) << "moleculeList must not grow";
    EXPECT_EQ(w.complexList.size(), comsBefore) << "complexList must not grow";
    EXPECT_TRUE(Molecule::emptyMolList.empty())
        << "Both the stale and the used index must be removed from emptyMolList";
    EXPECT_TRUE(Complex::emptyComList.empty())
        << "Both the stale and the used index must be removed from emptyComList";

    if (!cmacftr_indices_valid(w, newMolIndex, newComIndex))
        return;

    // The previously occupied slot 2 must be untouched.
    EXPECT_FALSE(w.moleculeList[2].isEmpty) << "The occupied molecule 2 must be left alone";
    EXPECT_EQ(w.moleculeList[2].index, 2) << "The occupied molecule 2 must keep its index";
    EXPECT_FALSE(w.complexList[2].isEmpty) << "The occupied complex 2 must be left alone";
}

// -----------------------------------------------------------------------------
// Test 4: monomerList tracking, driven by MolTemplate::canDestroy.
// -----------------------------------------------------------------------------
void test_cmacftr_monomer_list_tracking()
{
    std::cerr << "\n[TEST] test_cmacftr_monomer_list_tracking\n"
              << "  Source file: create_molecule_and_complex_from_transmission_rxn.cpp\n"
              << "  Function:    create_molecule_and_complex_from_transmission_rxn\n"
              << "  Scenario:    run once with MolTemplate::canDestroy == true and once with\n"
              << "               canDestroy == false.\n"
              << "  Pass:        the new molecule index is appended to monomerList only when\n"
              << "               the template is destroyable.\n";

    // --- destroyable template: the index must be recorded -------------------
    {
        CmacftrWorld w;
        w.build(/*canDestroy=*/true);
        w.molTemplateList[0].monomerList.clear();

        int newMolIndex = -1;
        int newComIndex = -1;
        cmacftr_run(w, newMolIndex, newComIndex);

        std::cerr << "  canDestroy == true  -> monomerList.size() = "
                  << w.molTemplateList[0].monomerList.size() << '\n';

        EXPECT_EQ(w.molTemplateList[0].monomerList.size(), 1u)
            << "A destroyable template must record the new monomer";
        if (!w.molTemplateList[0].monomerList.empty()) {
            EXPECT_EQ(w.molTemplateList[0].monomerList[0], newMolIndex)
                << "The recorded monomer must be the newly created molecule";
        }
    }

    // --- non-destroyable template: nothing must be recorded -----------------
    {
        CmacftrWorld w;
        w.build(/*canDestroy=*/false);
        w.molTemplateList[0].monomerList.clear();

        int newMolIndex = -1;
        int newComIndex = -1;
        cmacftr_run(w, newMolIndex, newComIndex);
        (void)newComIndex;

        std::cerr << "  canDestroy == false -> monomerList.size() = "
                  << w.molTemplateList[0].monomerList.size() << '\n';

        EXPECT_TRUE(w.molTemplateList[0].monomerList.empty())
            << "A non-destroyable template must not record monomers";
        EXPECT_GE(newMolIndex, 0) << "A molecule must still have been created";
    }
}

// -----------------------------------------------------------------------------
// Test 5: two consecutive transmissions must yield two distinct molecules and
//         complexes and bump the complex counter twice.
// -----------------------------------------------------------------------------
void test_cmacftr_two_successive_creations()
{
    std::cerr << "\n[TEST] test_cmacftr_two_successive_creations\n"
              << "  Source file: create_molecule_and_complex_from_transmission_rxn.cpp\n"
              << "  Function:    create_molecule_and_complex_from_transmission_rxn\n"
              << "  Scenario:    call the routine twice in a row on the same world.\n"
              << "  Pass:        two different molecule/complex slots are produced, the\n"
              << "               containers grow by two and numberOfComplexes grows by two.\n";

    CmacftrWorld w;
    w.build(/*canDestroy=*/false);

    const std::size_t molsBefore = w.moleculeList.size();
    const std::size_t comsBefore = w.complexList.size();
    const int numComplexesBefore = Complex::numberOfComplexes;

    int molA = -1;
    int comA = -1;
    cmacftr_run(w, molA, comA);

    int molB = -1;
    int comB = -1;
    cmacftr_run(w, molB, comB);

    std::cerr << "  -> first  (mol, com) = (" << molA << ", " << comA << ")\n"
              << "  -> second (mol, com) = (" << molB << ", " << comB << ")\n"
              << "  -> moleculeList.size() = " << w.moleculeList.size()
              << ", complexList.size() = " << w.complexList.size() << '\n';

    EXPECT_NE(molA, molB) << "Two creations must not share a molecule slot";
    EXPECT_NE(comA, comB) << "Two creations must not share a complex slot";
    EXPECT_EQ(w.moleculeList.size(), molsBefore + 2) << "moleculeList must grow by two";
    EXPECT_EQ(w.complexList.size(), comsBefore + 2) << "complexList must grow by two";
    EXPECT_EQ(Complex::numberOfComplexes, numComplexesBefore + 2)
        << "Complex::numberOfComplexes must be incremented once per creation";

    if (!cmacftr_indices_valid(w, molA, comA) || !cmacftr_indices_valid(w, molB, comB))
        return;

    // Each molecule must belong to its own complex.
    EXPECT_EQ(w.moleculeList[molA].myComIndex, comA)
        << "First molecule must point at the first new complex";
    EXPECT_EQ(w.moleculeList[molB].myComIndex, comB)
        << "Second molecule must point at the second new complex";
    EXPECT_EQ(w.complexList[comA].memberList.size(), 1u)
        << "First new complex must contain exactly one molecule";
    EXPECT_EQ(w.complexList[comB].memberList.size(), 1u)
        << "Second new complex must contain exactly one molecule";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - every named test_* routine is executed inside its own
// TEST() so results are reported individually while all of them still run.
// -----------------------------------------------------------------------------
TEST(CreateMolAndComplexFromTransmissionRxn, CreatesNewMoleculeAndComplex)
{
    test_cmacftr_creates_new_molecule_and_complex();
}

TEST(CreateMolAndComplexFromTransmissionRxn, ReusesEmptySlots)
{
    test_cmacftr_reuses_empty_slots();
}

TEST(CreateMolAndComplexFromTransmissionRxn, SkipsStaleEmptyEntries)
{
    test_cmacftr_skips_stale_empty_entries();
}

TEST(CreateMolAndComplexFromTransmissionRxn, MonomerListTracking)
{
    test_cmacftr_monomer_list_tracking();
}

TEST(CreateMolAndComplexFromTransmissionRxn, TwoSuccessiveCreations)
{
    test_cmacftr_two_successive_creations();
}