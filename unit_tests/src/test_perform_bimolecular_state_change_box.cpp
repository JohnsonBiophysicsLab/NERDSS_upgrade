/*! \file test_perform_bimolecular_state_change_box.cpp
 *
 * ### Unit test for src/reactions/perform_bimolecular_state_change_box.cpp
 *
 * Function under test:
 *
 *     void perform_bimolecular_state_change_box(int stateChangeIface, int facilitatorIface,
 *         std::array<int,3>& rxnItr, Molecule& stateChangeMol, Molecule& facilitatorMol,
 *         Complex& stateChangeCom, Complex& facilitatorCom, copyCounters& counterArrays,
 *         const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
 *         std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList,
 *         std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList,
 *         std::map<std::string,int>& observablesList, Membrane& membraneObject)
 *
 * This routine performs a *bimolecular state change* reaction inside a rectangular
 * (box) simulation volume.  Its job is to:
 *
 *   1. Look up the binding radius / product state / association angles of the
 *      reaction (either from the ForwardRxn or, for a back reaction, from the
 *      conjugate ForwardRxn plus the BackRxn product list).
 *   2. Copy both reacting complexes into temporary "association" coordinates and
 *      translate them so the two reacting interfaces are exactly `bindRadius`
 *      apart, splitting the displacement according to the diffusion constants.
 *   3. Orient the two complexes (skipped entirely when both molecules are points).
 *   4. Restore the pre-reaction center of mass of the complex pair, reflect off
 *      the box walls, and run the overlap / box-spanning sanity checks.
 *   5. If any check fails -> cancel: throw away the temporary coordinates and
 *      leave the system *completely unchanged*.
 *      If all checks pass -> commit: change the interface state, write the
 *      temporary coordinates into the real coordinates, mark the molecules as
 *      propagated, update the species copy numbers and observables, and clear
 *      the "crossed molecule" bookkeeping.
 *
 * ### Test strategy
 * A minimal two-molecule / two-complex system is constructed.  Both molecule
 * templates are flagged `isPoint = true` so that the orientation (theta / phi /
 * omega) machinery is bypassed - this makes the resulting geometry exactly
 * predictable and keeps the test independent of the association angle code.
 *
 *   * molecule 0 ("A") -> the molecule whose interface changes state, at (5,0,0)
 *   * molecule 1 ("B") -> the facilitator, at (0,0,0)
 *   * bindRadius = 1 nm, equal diffusion constants -> each complex moves 2 nm.
 *
 * Four scenarios are exercised:
 *   1. Forward reaction commits: state, copy numbers, observables, trajStatus,
 *      probability/crossbase bookkeeping.
 *   2. Forward reaction geometry: the pair ends up exactly bindRadius apart with
 *      the pair center of mass conserved.
 *   3. Cancelled reaction (structure overlap): nothing in the system changes.
 *   4. Back reaction: the state is reverted and the observable is decremented.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

namespace {

// -----------------------------------------------------------------------------
// Constants describing the two states of the reacting interface.
//
// The state list of interface "a" holds two states:
//    relative index 0 -> 'U' (absolute species index 2)
//    relative index 1 -> 'P' (absolute species index 3)
// -----------------------------------------------------------------------------
constexpr int kPbscbAbsStateU = 2; //!< absolute species index of the 'U' state
constexpr int kPbscbAbsStateP = 3; //!< absolute species index of the 'P' state
constexpr double kPbscbBindRadius = 1.0; //!< binding radius used by the reaction
constexpr double kPbscbBoxSide = 100.0; //!< cubic water box side length (nm)

/*! \brief Container for the whole miniature simulation system.
 *
 * Everything the function under test needs is grouped here so a test can build
 * it with one call and then hand references into the vectors to the function
 * (exactly like the production caller does).
 */
struct PbscbSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
    std::array<int, 3> rxnItr { 0, 0, 0 }; //!< {rxnIndex, rateIndex, isBackRxn}
};

/*! \brief Make sure the global GSL RNG exists.
 *
 * Some boundary-condition helpers may draw random numbers when they have to
 * resample a move.  The suite's main() defines `r` but we defensively seed it
 * here so a null generator can never crash this test.
 */
void pbscb_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  (global GSL RNG was null - seeding it so boundary helpers are safe)\n";
        srand_gsl(1);
    }
}

/*! \brief Build one MolTemplate describing a single-interface *point* molecule.
 *
 * `isPoint = true` is essential: it makes perform_bimolecular_state_change_box
 * skip all theta/phi/omega rotations, so the final coordinates are analytic.
 */
MolTemplate pbscb_make_template(int molTypeIndex, const std::string& name, bool checkOverlap)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0; // non-zero so complex radius updates are well defined
    temp.isPoint = true; // <- bypasses the orientation machinery
    temp.isRod = false;
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.checkOverlap = checkOverlap;

    // Non-zero 3D diffusion in every direction -> the reaction is treated as a
    // solution (3D) event, not a membrane (2D) event.
    temp.D.x = 10.0;
    temp.D.y = 10.0;
    temp.D.z = 10.0;
    temp.Dr.x = 0.1;
    temp.Dr.y = 0.1;
    temp.Dr.z = 0.1;

    // A single interface "a" with two possible states, 'U' and 'P'.
    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord(0.0, 0.0, 0.0); // point molecule: interface sits on the COM
    iface.stateList.emplace_back('U', kPbscbAbsStateU);
    iface.stateList.emplace_back('P', kPbscbAbsStateP);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a single-interface molecule whose interface lies on its COM. */
Molecule pbscb_make_molecule(int index, int molTypeIndex, int comIndex, const Coord& com, char stateIden,
    int absStateIndex, int relStateIndex)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = comIndex;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isEmpty = false;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = com; // point molecule
    iface.stateIden = stateIden;
    iface.stateIndex = relStateIndex;
    iface.index = absStateIndex;
    iface.relIndex = 0;
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a complex owning exactly one member molecule. */
Complex pbscb_make_complex(int index, const Coord& com, int memberMolIndex, int molTypeIndex)
{
    Complex targCom;
    targCom.index = index;
    targCom.comCoord = com;
    targCom.tmpComCoord = com;
    targCom.mass = 1.0;
    targCom.radius = 1.0;

    targCom.D.x = 10.0;
    targCom.D.y = 10.0;
    targCom.D.z = 10.0;
    targCom.Dr.x = 0.1;
    targCom.Dr.y = 0.1;
    targCom.Dr.z = 0.1;

    targCom.memberList.clear();
    targCom.memberList.push_back(memberMolIndex);
    targCom.numEachMol.assign(2, 0);
    targCom.numEachMol[molTypeIndex] = 1;
    targCom.lastNumberUpdateItrEachMol.assign(2, 0);

    targCom.isEmpty = false;
    targCom.OnSurface = false;
    targCom.tmpOnSurface = false;
    targCom.ncross = 0;
    targCom.trajStatus = TrajStatus::none;
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot.x = 0.0;
    targCom.trajRot.y = 0.0;
    targCom.trajRot.z = 0.0;

    return targCom;
}

/*! \brief Assemble the whole miniature system.
 *
 * \param[out] sys              system to populate
 * \param[in]  checkOverlap     value of MolTemplate::checkOverlap for both types.
 *                              `true` activates the structure-overlap check and
 *                              is used to force a cancellation.
 * \param[in]  startInStateP    if true molecule 0 starts in state 'P' (used for
 *                              the back-reaction test), otherwise in state 'U'.
 */
void pbscb_build_system(PbscbSystem& sys, bool checkOverlap, bool startInStateP)
{
    // ---- statics that Complex/MolTemplate bookkeeping relies on --------------
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };

    // ---- simulation parameters ---------------------------------------------
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 8;
    sys.params.numTotalComplex = 2;
    sys.params.timeStep = 1.0;
    sys.params.overlapSepLimit = 0.1; // well below the binding radius -> no cancel
    sys.params.scaleMaxDisplace = 1.0e12; // never reject on displacement grounds
    sys.params.nItr = 10;

    // ---- rectangular (box) boundary ---------------------------------------
    sys.membrane.isSphere = false;
    sys.membrane.isBox = true;
    sys.membrane.sphereR = 0.0;
    sys.membrane.hasCompartment = false;
    sys.membrane.implicitLipid = false;
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { kPbscbBoxSide, kPbscbBoxSide, kPbscbBoxSide });
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";

    // ---- molecule templates ------------------------------------------------
    sys.molTemplateList.push_back(pbscb_make_template(0, "A", checkOverlap));
    sys.molTemplateList.push_back(pbscb_make_template(1, "B", checkOverlap));

    // ---- molecules ---------------------------------------------------------
    // molecule 0: the state-changing molecule, 5 nm away on the +x axis
    sys.moleculeList.push_back(pbscb_make_molecule(0, 0, 0, Coord(5.0, 0.0, 0.0),
        startInStateP ? 'P' : 'U', startInStateP ? kPbscbAbsStateP : kPbscbAbsStateU, startInStateP ? 1 : 0));
    // molecule 1: the facilitator, at the origin
    sys.moleculeList.push_back(
        pbscb_make_molecule(1, 1, 1, Coord(0.0, 0.0, 0.0), 'U', kPbscbAbsStateU, 0));

    // ---- complexes (one molecule each) -------------------------------------
    sys.complexList.push_back(pbscb_make_complex(0, Coord(5.0, 0.0, 0.0), 0, 0));
    sys.complexList.push_back(pbscb_make_complex(1, Coord(0.0, 0.0, 0.0), 1, 1));

    // ---- forward reaction: B(a) + A(a~U) -> B(a) + A(a~P) ------------------
    ForwardRxn fwd;
    fwd.rxnType = ReactionType::biMolStateChange;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.bindRadius = kPbscbBindRadius;
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.hasStateChange = true;
    fwd.isObserved = true;
    fwd.observeLabel = "A_P";
    fwd.reactantListNew.emplace_back("B(a)", 1, kPbscbAbsStateU, 0, '\0', false);
    fwd.reactantListNew.emplace_back("A(a~U)", 0, kPbscbAbsStateU, 0, 'U', false);
    // productListNew[1] is the interface that actually changes state
    fwd.productListNew.emplace_back("B(a)", 1, kPbscbAbsStateU, 0, '\0', false);
    fwd.productListNew.emplace_back("A(a~P)", 0, kPbscbAbsStateP, 0, 'P', false);
    // assocAngles stay NaN - legal because both molecules are points.
    sys.forwardRxns.push_back(fwd);

    // ---- conjugate back reaction: A(a~P) -> A(a~U) -------------------------
    BackRxn back;
    back.rxnType = ReactionType::biMolStateChange;
    back.absRxnIndex = 0;
    back.relRxnIndex = 0;
    back.conjForwardRxnIndex = 0;
    back.hasStateChange = true;
    back.isObserved = true;
    back.observeLabel = "A_P";
    back.productListNew.emplace_back("B(a)", 1, kPbscbAbsStateU, 0, '\0', false);
    back.productListNew.emplace_back("A(a~U)", 0, kPbscbAbsStateU, 0, 'U', false);
    sys.backRxns.push_back(back);

    // ---- observables and species counters ---------------------------------
    sys.observablesList["A_P"] = 3;
    sys.counterArrays.copyNumSpecies.assign(8, 5);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: forward reaction commits - state, counters, observables, bookkeeping.
// -----------------------------------------------------------------------------
void test_pbscb_forward_changes_state_and_counters()
{
    std::cerr << "\n[TEST] test_pbscb_forward_changes_state_and_counters\n"
              << "  Source file:   perform_bimolecular_state_change_box.cpp\n"
              << "  Function:      perform_bimolecular_state_change_box (forward branch)\n"
              << "  Scenario:      two point molecules 5 nm apart in a 100 nm box undergo a\n"
              << "                 bimolecular state change (A(a~U) -> A(a~P)) with all\n"
              << "                 overlap checks disabled, so the reaction must commit.\n"
              << "  Pass criteria: interface state flips to 'P', copy numbers shift by -1/+1,\n"
              << "                 the observable is incremented, molecules are marked\n"
              << "                 propagated, temporary coords are cleared, and the\n"
              << "                 crossed-molecule bookkeeping is zeroed/cleared.\n";

    pbscb_ensure_rng();

    PbscbSystem sys;
    pbscb_build_system(sys, /*checkOverlap=*/false, /*startInStateP=*/false);
    sys.rxnItr = { 0, 0, 0 }; // forward reaction 0, rate state 0, NOT a back reaction

    // Bookkeeping we expect the function to clean up: each molecule "crossed"
    // the other with a non-zero reaction probability.
    sys.moleculeList[0].crossbase = std::vector<int> { 1 };
    sys.moleculeList[0].probvec = std::vector<double> { 0.75 };
    sys.moleculeList[1].crossbase = std::vector<int> { 0 };
    sys.moleculeList[1].probvec = std::vector<double> { 0.75 };

    // Snapsh