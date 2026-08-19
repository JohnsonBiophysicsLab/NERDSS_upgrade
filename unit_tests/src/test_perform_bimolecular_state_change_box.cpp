/*! \file test_perform_bimolecular_state_change_box.cpp
 *
 * ### Unit tests for src/reactions/perform_bimolecular_state_change_box.cpp
 *
 * Function under test:
 *
 *     void perform_bimolecular_state_change_box(int stateChangeIface,
 *              int facilitatorIface, std::array<int,3>& rxnItr,
 *              Molecule& stateChangeMol, Molecule& facilitatorMol,
 *              Complex& stateChangeCom, Complex& facilitatorCom,
 *              copyCounters& counterArrays, const Parameters& params,
 *              std::vector<ForwardRxn>& forwardRxns,
 *              std::vector<BackRxn>& backRxns,
 *              std::vector<Molecule>& moleculeList,
 *              std::vector<Complex>& complexList,
 *              std::vector<MolTemplate>& molTemplateList,
 *              std::map<std::string,int>& observablesList,
 *              Membrane& membraneObject)
 *
 * The routine performs a *bimolecular state change* inside a rectangular
 * (box) simulation volume:
 *   1. Both complexes are copied into temporary association coordinates.
 *   2. The two reacting interfaces are brought to contact (separation ==
 *      ForwardRxn::bindRadius) using the diffusion-constant weighted split.
 *   3. If neither molecule is a "point" molecule the association angles
 *      (theta/phi/omega) are enforced by rotation.
 *   4. The center of mass of the *pair* is forced back to its pre-reaction
 *      value and the pair is reflected off of the box walls.
 *   5. Overlap / box-spanning checks may CANCEL the whole event, in which
 *      case the temporary coordinates are discarded and nothing changes.
 *   6. On success the state of the state-changing interface is switched to
 *      the product state, the temporary coordinates are committed, species
 *      copy numbers and observables are updated, and the crossing
 *      bookkeeping (probvec / crossbase / ncross) is cleaned up.
 *
 * The tests below build a deliberately tiny, fully initialized two-molecule
 * system so that every branch of interest can be exercised and the expected
 * geometry can be computed by hand.
 */

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers / fixture data.  Everything is prefixed with "pbscb_" so no
// symbol collides with the rest of the generated test suite.
// -----------------------------------------------------------------------------
namespace {

//! Absolute interface-state indices of the tiny system used by these tests.
constexpr int kPbscbIfaceA_Abs { 0 }; //!< A(a)    -- facilitator, stateless
constexpr int kPbscbIfaceB_U_Abs { 1 }; //!< B(b~U) -- state changer, reactant state
constexpr int kPbscbIfaceB_P_Abs { 2 }; //!< B(b~P) -- state changer, product state

//! Relative state indices, i.e. positions inside Interface::stateList.
//! NOTE: the implementation indexes copyCounters::copyNumSpecies with the
//! *relative* stateIndex (not the absolute one), so the expectations below
//! use these values on purpose -- they mirror what the code actually does.
constexpr int kPbscbRelStateU { 0 };
constexpr int kPbscbRelStateP { 1 };

//! Binding radius of the state-change reaction (nm).
constexpr double kPbscbBindRadius { 5.0 };

//! Starting |x| of each reactant, so their separation is 16 nm.
constexpr double kPbscbStartX { 8.0 };

/*! \brief Make sure the shared GSL random number generator is usable.
 *
 * Some of the boundary-condition helpers called downstream can resample a
 * trajectory, which requires the global generator. The suite's main() only
 * declares `r`, so allocate/seed it here if nobody has done so yet.
 */
void pbscb_ensure_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}

/*! \brief Container holding a complete, self-consistent mini simulation. */
struct PbscbSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    copyCounters counterArrays {};
    std::map<std::string, int> observablesList {};
};

/*! \brief Build the two-molecule test system.
 *
 * \param[in] bothArePoints  when true both MolTemplates are flagged isPoint,
 *                           which makes the routine skip all rotations.
 * \param[in] checkOverlap   when true both templates request overlap checking
 *                           and overlapSepLimit is set larger than the binding
 *                           radius, so the association must be cancelled.
 * \param[in] initialBState  'U' or 'P': the starting state of the interface
 *                           that is going to change.
 * \param[in] observeLabel   observable label stored on the reactions.
 */
PbscbSystem pbscb_build_system(bool bothArePoints, bool checkOverlap, char initialBState,
    const std::string& observeLabel)
{
    PbscbSystem sys;

    /* ------------------------- parameters ------------------------- */
    sys.params.rank = 0;
    sys.params.timeStep = 0.1; // us
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 3;
    sys.params.numTotalComplex = 2;
    // COM-COM distance below which association is cancelled; only applied to
    // molecules whose template sets checkOverlap.
    sys.params.overlapSepLimit = checkOverlap ? 10.0 : 0.1;
    sys.params.scaleMaxDisplace = 1.0e6; // never reject on displacement here

    /* ------------------------- boundary --------------------------- */
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    sys.membrane.isBox = true;
    sys.membrane.isSphere = false;
    sys.membrane.implicitLipid = false;
    sys.membrane.hasCompartment = false;
    sys.membrane.TwoD = false;
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";
    sys.membrane.nSites = 0;
    sys.membrane.nStates = 0;
    sys.membrane.No_free_lipids = 0;
    sys.membrane.No_protein = 0;
    sys.membrane.totalSA = 0.0;
    sys.membrane.Dx = 0.0;
    sys.membrane.Dy = 0.0;
    sys.membrane.Dz = 0.0;
    sys.membrane.Drx = 0.0;
    sys.membrane.Dry = 0.0;
    sys.membrane.Drz = 0.0;
    sys.membrane.offset = 0.0;
    sys.membrane.lipidLength = 0.0;

    /* --------------------- molecule templates -------------------- */
    sys.molTemplateList.resize(2);

    // Template 0: "A", the facilitator. One stateless interface "a".
    MolTemplate& tempA = sys.molTemplateList[0];
    tempA.molName = "A";
    tempA.molTypeIndex = 0;
    tempA.copies = 1;
    tempA.mass = 1.0;
    tempA.radius = 2.0;
    tempA.D = Coord { 10.0, 10.0, 10.0 };
    tempA.Dr = Coord { 0.1, 0.1, 0.1 };
    tempA.isPoint = bothArePoints;
    tempA.isRod = false;
    tempA.isLipid = false;
    tempA.isImplicitLipid = false;
    tempA.isPromoter = false;
    tempA.checkOverlap = checkOverlap;
    {
        Interface ifaceA;
        ifaceA.name = "a";
        ifaceA.index = 0;
        ifaceA.iCoord = bothArePoints ? Coord { 0.0, 0.0, 0.0 } : Coord { 1.0, 0.5, 0.0 };
        ifaceA.stateList.emplace_back("a", '\0', kPbscbIfaceA_Abs);
        tempA.interfaceList.push_back(ifaceA);
    }

    // Template 1: "B", the state changer. One interface "b" with states U and P.
    MolTemplate& tempB = sys.molTemplateList[1];
    tempB.molName = "B";
    tempB.molTypeIndex = 1;
    tempB.copies = 1;
    tempB.mass = 1.0;
    tempB.radius = 2.0;
    tempB.D = Coord { 10.0, 10.0, 10.0 };
    tempB.Dr = Coord { 0.1, 0.1, 0.1 };
    tempB.isPoint = bothArePoints;
    tempB.isRod = false;
    tempB.isLipid = false;
    tempB.isImplicitLipid = false;
    tempB.isPromoter = false;
    tempB.checkOverlap = checkOverlap;
    {
        Interface ifaceB;
        ifaceB.name = "b";
        ifaceB.index = 0;
        ifaceB.iCoord = bothArePoints ? Coord { 0.0, 0.0, 0.0 } : Coord { -1.0, -0.5, 0.0 };
        ifaceB.stateList.emplace_back("b~U", 'U', kPbscbIfaceB_U_Abs);
        ifaceB.stateList.emplace_back("b~P", 'P', kPbscbIfaceB_P_Abs);
        tempB.interfaceList.push_back(ifaceB);
        tempB.ifacesWithStates.push_back(0);
    }

    /* --------------------------- molecules ------------------------ */
    const int bRelState { (initialBState == 'P') ? kPbscbRelStateP : kPbscbRelStateU };
    const int bAbsState { (initialBState == 'P') ? kPbscbIfaceB_P_Abs : kPbscbIfaceB_U_Abs };

    const Coord comA { -kPbscbStartX, 0.0, 0.0 };
    const Coord comB { kPbscbStartX, 0.0, 0.0 };
    const Coord ifaceCrdA = bothArePoints ? comA : Coord { comA.x + 1.0, comA.y + 0.5, comA.z };
    const Coord ifaceCrdB = bothArePoints ? comB : Coord { comB.x - 1.0, comB.y - 0.5, comB.z };

    // Molecule 0: the facilitator A, sole member of complex 0.
    {
        Molecule mol;
        mol.index = 0;
        mol.id = 0;
        mol.myComIndex = 0;
        mol.complexId = 0;
        mol.molTypeIndex = 0;
        mol.mass = 1.0;
        mol.isEmpty = false;
        mol.isLipid = false;
        mol.isImplicitLipid = false;
        mol.isPromoter = false;
        mol.linksToSurface = 0;
        mol.trajStatus = TrajStatus::none;
        mol.comCoord = comA;

        Molecule::Iface iface;
        iface.coord = ifaceCrdA;
        iface.index = kPbscbIfaceA_Abs;
        iface.relIndex = 0;
        iface.stateIden = '\0';
        iface.stateIndex = 0;
        iface.molTypeIndex = 0;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
        mol.freelist.push_back(0);

        sys.moleculeList.push_back(mol);
    }

    // Molecule 1: the state changer B, sole member of complex 1.
    {
        Molecule mol;
        mol.index = 1;
        mol.id = 1;
        mol.myComIndex = 1;
        mol.complexId = 1;
        mol.molTypeIndex = 1;
        mol.mass = 1.0;
        mol.isEmpty = false;
        mol.isLipid = false;
        mol.isImplicitLipid = false;
        mol.isPromoter = false;
        mol.linksToSurface = 0;
        mol.trajStatus = TrajStatus::none;
        mol.comCoord = comB;

        Molecule::Iface iface;
        iface.coord = ifaceCrdB;
        iface.index = bAbsState;
        iface.relIndex = 0;
        iface.stateIden = initialBState;
        iface.stateIndex = bRelState;
        iface.molTypeIndex = 1;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
        mol.freelist.push_back(0);

        sys.moleculeList.push_back(mol);
    }

    /* --------------------------- complexes ------------------------ */
    for (int comItr = 0; comItr < 2; ++comItr) {
        Complex com;
        com.index = comItr;
        com.id = comItr;
        com.ownerRank = 0;
        com.comCoord = sys.moleculeList[comItr].comCoord;
        com.mass = 1.0;
        com.radius = 2.0;
        com.D = Coord { 10.0, 10.0, 10.0 };
        com.Dr = Coord { 0.1, 0.1, 0.1 };
        com.memberList.push_back(comItr);
        com.numEachMol = std::vector<int>(2, 0);
        com.numEachMol[sys.moleculeList[comItr].molTypeIndex] = 1;
        com.lastNumberUpdateItrEachMol = std::vector<long long int>(2, 0);
        com.isEmpty = false;
        com.OnSurface = false;
        com.onFiber = false;
        com.linksToSurface = 0;
        com.ncross = 0;
        com.trajStatus = TrajStatus::none;
        sys.complexList.push_back(com);
    }

    /* --------------------------- reactions ------------------------ */
    {
        ForwardRxn fwd;
        fwd.rxnType = ReactionType::biMolStateChange;
        fwd.absRxnIndex = 0;
        fwd.relRxnIndex = 0;
        fwd.isReversible = true;
        fwd.conjBackRxnIndex = 0;
        fwd.hasStateChange = true;
        fwd.isOnMem = false;
        fwd.bindRadius = kPbscbBindRadius;
        fwd.norm1 = Vector(0.0, 0.0, 1.0);
        fwd.norm2 = Vector(0.0, 0.0, 1.0);
        fwd.isObserved = !observeLabel.empty();
        fwd.observeLabel = observeLabel;
        fwd.rxnLabel = "stateChange";

        // Angles: for point molecules everything stays NaN so all rotations are
        // skipped. Otherwise use a generic (non-degenerate) theta and leave
        // phi/omega undefined so only theta_rotation() runs.
        if (!bothArePoints) {
            fwd.assocAngles = ForwardRxn::Angles(1.5, 1.5,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN());
        }

        fwd.reactantListNew.emplace_back("a", 0, kPbscbIfaceA_Abs, 0, '\0', false);
        fwd.reactantListNew.emplace_back("b", 1, kPbscbIfaceB_U_Abs, 0, 'U', false);
        fwd.productListNew.emplace_back("a", 0, kPbscbIfaceA_Abs, 0, '\0', false);
        // productListNew[1] is what the routine reads as the new state.
        fwd.productListNew.emplace_back("b", 1, kPbscbIfaceB_P_Abs, 0, 'P', false);
        fwd.stateChangeIface = std::make_pair(fwd.reactantListNew[1], fwd.productListNew[1]);
        fwd.intReactantList = { kPbscbIfaceA_Abs, kPbscbIfaceB_U_Abs };
        fwd.intProductList = { kPbscbIfaceA_Abs, kPbscbIfaceB_P_Abs };
        fwd.rateList.emplace_back();
        fwd.rateList.back().rate = 10.0;

        sys.forwardRxns.push_back(fwd);
    }
    {
        // The conjugate back reaction simply swaps reactants and products, so
        // its productListNew[1] restores state 'U'.
        BackRxn bck;
        bck.rxnType = ReactionType::biMolStateChange;
        bck.absRxnIndex = 1;
        bck.relRxnIndex = 0;
        bck.conjForwardRxnIndex = 0;
        bck.hasStateChange = true;
        bck.isOnMem = false;
        bck.isObserved = !observeLabel.empty();
        bck.observeLabel = observeLabel;
        bck.reactantListNew = sys.forwardRxns[0].productListNew;
        bck.productListNew = sys.forwardRxns[0].reactantListNew;
        bck.stateChangeIface = std::make_pair(bck.reactantListNew[1], bck.productListNew[1]);
        bck.rateList.emplace_back();
        bck.rateList.back().rate = 1.0;

        sys.backRxns.push_back(bck);
    }

    /* ----------------------- counters / observables --------------- */
    // Indexed by the *relative* state index (see note at the top of the file).
    sys.counterArrays.copyNumSpecies = std::vector<int> { 7, 3, 11, 0, 0 };
    if (!observeLabel.empty())
        sys.observablesList[observeLabel] = 20;

    return sys;
}

/*! \brief Append an unrelated "spectator" molecule + complex to the system.
 *
 * The spectator sits far away from the reaction and is used to verify that the
 * routine zeroes the reaction probabilities that other molecules hold against
 * the two reacting molecules.
 */
void pbscb_add_spectator(PbscbSystem& sys)
{
    Molecule mol;
    mol.index = 2;
    mol.id = 2;
    mol.myComIndex = 2;
    mol.complexId = 2;
    mol.molTypeIndex = 0; // same template as the facilitator
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = Coord { 30.0, 30.0, 30.0 };

    Molecule::Iface iface;
    iface.coord = mol.comCoord;
    iface.index = kPbscbIfaceA_Abs;
    iface.relIndex = 0;
    iface.stateIden = '\0';
    iface.stateIndex = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    // The spectator "encountered" both reactants this step.
    mol.crossbase = { 0, 1 };
    mol.mycrossint = { 0, 0 };
    mol.crossrxn = { std::array<int, 3> { 0, 0, 0 }, std::array<int, 3> { 0, 0, 0 } };
    mol.probvec = { 0.70, 0.30 };
    sys.moleculeList.push_back(mol);

    Complex com;
    com.index = 2;
    com.id = 2;
    com.ownerRank = 0;
    com.comCoord = mol.comCoord;
    com.mass = 1.0;
    com.radius = 2.0;
    com.D = Coord { 10.0, 10.0, 10.0 };
    com.Dr = Coord { 0.1, 0.1, 0.1 };
    com.memberList.push_back(2);
    com.numEachMol = std::vector<int> { 1, 0 };
    com.lastNumberUpdateItrEachMol = std::vector<long long int>(2, 0);
    com.isEmpty = false;
    com.OnSurface = false;
    com.ncross = 2;
    com.trajStatus = TrajStatus::none;
    sys.complexList.push_back(com);

    // ...and each reactant "encountered" the spectator.
    sys.moleculeList[0].crossbase = { 2 };
    sys.moleculeList[0].mycrossint = { 0 };
    sys.moleculeList[0].crossrxn = { std::array<int, 3> { 0, 0, 0 } };
    sys.moleculeList[0].probvec = { 0.50 };

    sys.moleculeList[1].crossbase = { 2 };
    sys.moleculeList[1].mycrossint = { 0 };
    sys.moleculeList[1].crossrxn = { std::array<int, 3> { 0, 0, 0 } };
    sys.moleculeList[1].probvec = { 0.40 };
}

/*! \brief Distance between the two reacting interfaces (molecules 0 and 1). */
double pbscb_iface_separation(const PbscbSystem& sys)
{
    const Coord& c1 = sys.moleculeList[0].interfaceList[0].coord;
    const Coord& c2 = sys.moleculeList[1].interfaceList[0].coord;
    const double dx { c1.x - c2.x };
    const double dy { c1.y - c2.y };
    const double dz { c1.z - c2.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Convenience wrapper calling the function under test on a system.
 *
 * \param[in] rxnItr {reaction index, rate index, isBackRxn}
 */
void pbscb_invoke(PbscbSystem& sys, std::array<int, 3> rxnItr)
{
    const int stateChangeIface { 0 }; // relative index of B(b)
    const int facilitatorIface { 0 }; // relative index of A(a)

    perform_bimolecular_state_change_box(stateChangeIface, facilitatorIface, rxnItr,
        sys.moleculeList[1], // stateChangeMol  (B)
        sys.moleculeList[0], // facilitatorMol  (A)
        sys.complexList[1], // stateChangeCom
        sys.complexList[0], // facilitatorCom
        sys.counterArrays, sys.params, sys.forwardRxns, sys.backRxns, sys.moleculeList,
        sys.complexList, sys.molTemplateList, sys.observablesList, sys.membrane);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: forward state change with two point molecules.
//         Because both templates are points every rotation is skipped, so the
//         final geometry can be computed analytically.
// -----------------------------------------------------------------------------
void test_pbscb_forward_point_molecules()
{
    std::cerr << "\n[TEST] test_pbscb_forward_point_molecules\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_box.cpp\n"
              << "  Function:    perform_bimolecular_state_change_box\n"
              << "  Scenario:    two point molecules 16 nm apart, forward bimolecular\n"
              << "               state change with bindRadius = " << kPbscbBindRadius << " nm.\n"
              << "  Criteria:    (a) interfaces end up exactly bindRadius apart,\n"
              << "               (b) the state-change interface switches U -> P,\n"
              << "               (c) copyNumSpecies / observable / bookkeeping updated.\n";

    pbscb_ensure_rng();
    PbscbSystem sys = pbscb_build_system(/*bothArePoints=*/true, /*checkOverlap=*/false,
        /*initialBState=*/'U', /*observeLabel=*/"Bphos");

    const std::vector<int> copyBefore { sys.counterArrays.copyNumSpecies };
    const int obsBefore { sys.observablesList["Bphos"] };

    std::cerr << "  Pre-reaction: A.com.x = " << sys.moleculeList[0].comCoord.x
              << ", B.com.x = " << sys.moleculeList[1].comCoord.x
              << ", separation = " << pbscb_iface_separation(sys) << " nm\n";

    // rxnItr = {forward reaction 0, rate 0, isBackRxn = 0}
    pbscb_invoke(sys, std::array<int, 3> { 0, 0, 0 });

    std::cerr << "  Post-reaction: A.com.x = " << sys.moleculeList[0].comCoord.x
              << ", B.com.x = " << sys.moleculeList[1].comCoord.x
              << ", separation = " << pbscb_iface_separation(sys) << " nm\n";

    // (a) Geometry. Both complexes have identical diffusion constants, so each
    //     is displaced by half of (16 - 5) = 5.5 nm along x, and the pair COM
    //     (the origin) is unchanged, so no extra shift is applied.
    EXPECT_NEAR(pbscb_iface_separation(sys), kPbscbBindRadius, 1e-9)
        << "Reacting interfaces must be exactly bindRadius apart after association";
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, -kPbscbBindRadius / 2.0, 1e-9)
        << "Facilitator should be moved to -bindRadius/2 along x";
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, kPbscbBindRadius / 2.0, 1e-9)
        << "State changer should be moved to +bindRadius/2 along x";
    EXPECT_NEAR(sys.moleculeList[0].comCoord.y, 0.0, 1e-9) << "y should be untouched";
    EXPECT_NEAR(sys.moleculeList[0].comCoord.z, 0.0, 1e-9) << "z should be untouched";

    // (b) The state change itself.
    const Molecule::Iface& bIface = sys.moleculeList[1].interfaceList[0];
    EXPECT_EQ(bIface.stateIden, 'P') << "State identity should now be the product state 'P'";
    EXPECT_EQ(bIface.index, kPbscbIfaceB_P_Abs) << "Absolute interface index should be the product index";
    EXPECT_EQ(bIface.stateIndex, kPbscbRelStateP) << "Relative state index should point at stateList[1]";

    // The facilitator's interface must be untouched.
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, kPbscbIfaceA_Abs)
        << "Facilitator interface index must not change";

    // (c) Species counters: the reactant state loses one copy, the product gains one.
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscbRelStateU], copyBefore[kPbscbRelStateU] - 1)
        << "Reactant state copy number should be decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscbRelStateP], copyBefore[kPbscbRelStateP] + 1)
        << "Product state copy number should be incremented";

    // Observable for a forward reaction is incremented.
    EXPECT_EQ(sys.observablesList["Bphos"], obsBefore + 1)
        << "Observable of an observed forward reaction should be incremented";

    // Temporary association coordinates must have been cleared and committed.
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty())
        << "Facilitator temporary interface coordinates should be cleared";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty())
        << "State changer temporary interface coordinates should be cleared";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "Facilitator should be marked as propagated";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "State changer should be marked as propagated";

    // Crossing bookkeeping is reset so these molecules react only once.
    EXPECT_EQ(sys.complexList[0].ncross, -1) << "Facilitator complex ncross should be set to -1";
    EXPECT_EQ(sys.complexList[1].ncross, -1) << "State changer complex ncross should be set to -1";
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty()) << "Facilitator crossbase should be cleared";
    EXPECT_TRUE(sys.moleculeList[1].crossbase.empty()) << "State changer crossbase should be cleared";

    // Complex properties are recomputed, so the single-member complex COM must
    // track its member molecule.
    EXPECT_NEAR(sys.complexList[0].comCoord.x, sys.moleculeList[0].comCoord.x, 1e-9)
        << "Facilitator complex COM should be updated from its member";
    EXPECT_NEAR(sys.complexList[1].comCoord.x, sys.moleculeList[1].comCoord.x, 1e-9)
        << "State changer complex COM should be updated from its member";
}

// -----------------------------------------------------------------------------
// Test 2: the conjugate BACK reaction (P -> U), selected with rxnItr[2] == 1.
//         The routine must look up bindRadius/angles through
//         BackRxn::conjForwardRxnIndex and take the new state from
//         BackRxn::productListNew[1].
// -----------------------------------------------------------------------------
void test_pbscb_back_reaction_state_change()
{
    std::cerr << "\n[TEST] test_pbscb_back_reaction_state_change\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_box.cpp\n"
              << "  Function:    perform_bimolecular_state_change_box (back-reaction branch)\n"
              << "  Scenario:    B starts in state 'P'; rxnItr = {0,0,1} selects backRxns[0].\n"
              << "  Criteria:    state reverts P -> U, counters swap, observed observable\n"
              << "               is DECREMENTED, and geometry still lands at bindRadius.\n";

    pbscb_ensure_rng();
    PbscbSystem sys = pbscb_build_system(/*bothArePoints=*/true, /*checkOverlap=*/false,
        /*initialBState=*/'P', /*observeLabel=*/"Bphos");

    const std::vector<int> copyBefore { sys.counterArrays.copyNumSpecies };
    const int obsBefore { sys.observablesList["Bphos"] };

    // rxnItr = {back reaction 0, rate 0, isBackRxn = 1}
    pbscb_invoke(sys, std::array<int, 3> { 0, 0, 1 });

    const Molecule::Iface& bIface = sys.moleculeList[1].interfaceList[0];
    std::cerr << "  Post-reaction state: '" << bIface.stateIden << "' (abs index "
              << bIface.index << ", rel index " << bIface.stateIndex << ")\n";

    EXPECT_EQ(bIface.stateIden, 'U') << "Back reaction should restore state 'U'";
    EXPECT_EQ(bIface.index, kPbscbIfaceB_U_Abs) << "Absolute index should be the 'U' state index";
    EXPECT_EQ(bIface.stateIndex, kPbscbRelStateU) << "Relative state index should point at stateList[0]";

    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscbRelStateP], copyBefore[kPbscbRelStateP] - 1)
        << "'P' copy number should be decremented for the back reaction";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[kPbscbRelStateU], copyBefore[kPbscbRelStateU] + 1)
        << "'U' copy number should be incremented for the back reaction";

    EXPECT_EQ(sys.observablesList["Bphos"], obsBefore - 1)
        << "Observable of an observed BACK reaction should be decremented";

    // The binding radius came from the conjugate forward reaction.
    EXPECT_NEAR(pbscb_iface_separation(sys), kPbscbBindRadius, 1e-9)
        << "Back reaction must still place interfaces at the forward reaction's bindRadius";
}

// -----------------------------------------------------------------------------
// Test 3: the association is CANCELLED by the structure-overlap check.
//         Both templates set checkOverlap and overlapSepLimit (10 nm) exceeds
//         the binding radius (5 nm), so bringing the molecules to contact is
//         necessarily an overlap.  Nothing at all may change.
// -----------------------------------------------------------------------------
void test_pbscb_cancelled_by_overlap()
{
    std::cerr << "\n[TEST] test_pbscb_cancelled_by_overlap\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_box.cpp\n"
              << "  Function:    perform_bimolecular_state_change_box (cancel branch)\n"
              << "  Scenario:    checkOverlap = true and overlapSepLimit (10 nm) >\n"
              << "               bindRadius (5 nm), so the event must be rejected.\n"
              << "  Criteria:    coordinates, interface state, copy numbers, observables,\n"
              << "               trajStatus, ncross and crossbase are all unchanged, and\n"
              << "               the temporary coordinates are discarded.\n";

    pbscb_ensure_rng();
    PbscbSystem sys = pbscb_build_system(/*bothArePoints=*/true, /*checkOverlap=*/true,
        /*initialBState=*/'U', /*observeLabel=*/"Bphos");

    const std::vector<int> copyBefore { sys.counterArrays.copyNumSpecies };
    const int obsBefore { sys.observablesList["Bphos"] };
    const double sepBefore { pbscb_iface_separation(sys) };

    pbscb_invoke(sys, std::array<int, 3> { 0, 0, 0 });

    std::cerr << "  Separation before = " << sepBefore
              << " nm, after = " << pbscb_iface_separation(sys) << " nm\n";

    // Coordinates untouched.
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, -kPbscbStartX, 1e-12)
        << "Facilitator must not move when association is cancelled";
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, kPbscbStartX, 1e-12)
        << "State changer must not move when association is cancelled";
    EXPECT_NEAR(pbscb_iface_separation(sys), sepBefore, 1e-12)
        << "Interface separation must be unchanged when association is cancelled";

    // State untouched.
    const Molecule::Iface& bIface = sys.moleculeList[1].interfaceList[0];
    EXPECT_EQ(bIface.stateIden, 'U') << "State must remain 'U' after cancellation";
    EXPECT_EQ(bIface.index, kPbscbIfaceB_U_Abs) << "Absolute index must remain the 'U' index";
    EXPECT_EQ(bIface.stateIndex, kPbscbRelStateU) << "Relative state index must remain 0";

    // Counters and observables untouched.
    EXPECT_EQ(sys.counterArrays.copyNumSpecies, copyBefore)
        << "copyNumSpecies must be untouched when association is cancelled";
    EXPECT_EQ(sys.observablesList["Bphos"], obsBefore)
        << "Observables must be untouched when association is cancelled";

    // Temporary coordinates are cleared even on the cancel path.
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty())
        << "Temporary coordinates must be cleared on cancellation (facilitator)";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty())
        << "Temporary coordinates must be cleared on cancellation (state changer)";

    // Bookkeeping is NOT reset when the event is rejected.
    EXPECT_EQ(static_cast<int>(sys.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must remain 'none' on cancellation (facilitator)";
    EXPECT_EQ(static_cast<int>(sys.moleculeList[1].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must remain 'none' on cancellation (state changer)";
    EXPECT_EQ(sys.complexList[0].ncross, 0) << "ncross must be untouched on cancellation";
    EXPECT_EQ(sys.complexList[1].ncross, 0) << "ncross must be untouched on cancellation";
}

// -----------------------------------------------------------------------------
// Test 4: non-point molecules, so the theta rotations are actually executed.
//         phi/omega are left as NaN so only theta_rotation() runs.
// -----------------------------------------------------------------------------
void test_pbscb_rotation_path_non_point()
{
    std::cerr << "\n[TEST] test_pbscb_rotation_path_non_point\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_box.cpp\n"
              << "  Function:    perform_bimolecular_state_change_box (rotation path)\n"
              << "  Scenario:    both molecules have off-center interfaces (isPoint = false)\n"
              << "               and theta1 = theta2 = 1.5 rad, phi/omega undefined.\n"
              << "  Criteria:    interfaces still end at bindRadius, the theta angle between\n"
              << "               (COM - iface) and sigma matches the requested value, and the\n"
              << "               state change plus bookkeeping still happen.\n";

    pbscb_ensure_rng();
    PbscbSystem sys = pbscb_build_system(/*bothArePoints=*/false, /*checkOverlap=*/false,
        /*initialBState=*/'U', /*observeLabel=*/"Bphos");

    pbscb_invoke(sys, std::array<int, 3> { 0, 0, 0 });

    const double sep { pbscb_iface_separation(sys) };
    std::cerr << "  Final interface separation = " << sep << " nm (target "
              << kPbscbBindRadius << " nm)\n";

    // Rotations pivot about the reacting interfaces, so sigma is conserved.
    EXPECT_NEAR(sep, kPbscbBindRadius, 1e-6)
        << "Rotations must preserve the binding-radius separation";

    // Measure theta1: angle between (COM1 - iface1) and sigma = iface1 - iface2.
    const Coord& com1 = sys.moleculeList[0].comCoord;
    const Coord& if1 = sys.moleculeList[0].interfaceList[0].coord;
    const Coord& if2 = sys.moleculeList[1].interfaceList[0].coord;
    Vector v1 { com1 - if1 };
    Vector sigma { if1 - if2 };
    v1.calc_magnitude();
    sigma.calc_magnitude();
    const double theta1 { v1.dot_theta(sigma) };
    std::cerr << "  Measured theta1 = " << theta1 << " rad (requested 1.5 rad; the "
              << "sigma sign convention may report pi - theta)\n";

    // Depending on the sign convention used for sigma, the measured angle is
    // either the requested theta or its supplement -- both are acceptable.
    EXPECT_TRUE(std::abs(theta1 - 1.5) < 1e-6 || std::abs((M_PI - theta1) - 1.5) < 1e-6)
        << "theta1 should equal the requested association angle (or its supplement); got "
        << theta1;

    // The state change and bookkeeping still take place.
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].stateIden, 'P')
        << "State should change to 'P' on the rotation path as well";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].stateIndex, kPbscbRelStateP)
        << "Relative state index should be updated on the rotation path";
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty())
        << "Temporary coordinates should be cleared after a successful rotation path";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty())
        << "Temporary coordinates should be cleared after a successful rotation path";
    EXPECT_EQ(sys.complexList[0].ncross, -1) << "ncross should be reset on success";
    EXPECT_EQ(sys.complexList[1].ncross, -1) << "ncross should be reset on success";
}

// -----------------------------------------------------------------------------
// Test 5: observed reaction whose label is missing from observablesList.
//         The routine warns and continues; the state change must still occur
//         and no new observable may be created.
// -----------------------------------------------------------------------------
void test_pbscb_missing_observable_label()
{
    std::cerr << "\n[TEST] test_pbscb_missing_observable_label\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_box.cpp\n"
              << "  Function:    perform_bimolecular_state_change_box (observable lookup)\n"
              << "  Scenario:    the reaction is flagged isObserved but its label is not\n"
              << "               present in observablesList.\n"
              << "  Criteria:    no new map entry is created and the state change still\n"
              << "               completes successfully.\n";

    pbscb_ensure_rng();
    PbscbSystem sys = pbscb_build_system(/*bothArePoints=*/true, /*checkOverlap=*/false,
        /*initialBState=*/'U', /*observeLabel=*/"Bphos");

    // Rename the label on the reaction so the lookup fails.
    sys.forwardRxns[0].observeLabel = "labelThatDoesNotExist";
    const std::size_t mapSizeBefore { sys.observablesList.size() };
    const int knownObsBefore { sys.observablesList["Bphos"] };

    pbscb_invoke(sys, std::array<int, 3> { 0, 0, 0 });

    std::cerr << "  observablesList size before = " << mapSizeBefore
              << ", after = " << sys.observablesList.size() << '\n';

    EXPECT_EQ(sys.observablesList.size(), mapSizeBefore)
        << "A missing observable label must not insert a new map entry";
    EXPECT_EQ(sys.observablesList["Bphos"], knownObsBefore)
        << "Unrelated observables must not be modified";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].stateIden, 'P')
        << "The state change must still happen when the observable is unknown";
}

// -----------------------------------------------------------------------------
// Test 6: crossing bookkeeping. A third "spectator" molecule holds non-zero
//         reaction probabilities against both reactants; those probabilities
//         must be zeroed so the reacted molecules cannot react again this step
//         (while the spectator still avoids overlapping them).
// -----------------------------------------------------------------------------
void test_pbscb_zeroes_partner_probabilities()
{
    std::cerr << "\n[TEST] test_pbscb_zeroes_partner_probabilities\n"
              << "  Source file: src/reactions/perform_bimolecular_state_change_box.cpp\n"
              << "  Function:    perform_bimolecular_state_change_box (crossing cleanup)\n"
              << "  Scenario:    a spectator molecule lists both reactants in its crossbase\n"
              << "               with probvec = {0.70, 0.30}.\n"
              << "  Criteria:    both spectator probabilities are set to 0, the reactants'\n"
              << "               crossbase lists are cleared and ncross becomes -1, while\n"
              << "               the spectator's own ncross is untouched.\n";

    pbscb_ensure_rng();
    PbscbSystem sys = pbscb_build_system(/*bothArePoints=*/true, /*checkOverlap=*/false,
        /*initialBState=*/'U', /*observeLabel=*/"Bphos");
    pbscb_add_spectator(sys);

    std::cerr << "  Spectator probvec before = {" << sys.moleculeList[2].probvec[0] << ", "
              << sys.moleculeList[2].probvec[1] << "}\n";

    pbscb_invoke(sys, std::array<int, 3> { 0, 0, 0 });

    std::cerr << "  Spectator probvec after  = {" << sys.moleculeList[2].probvec[0] << ", "
              << sys.moleculeList[2].probvec[1] << "}\n";

    ASSERT_EQ(sys.moleculeList[2].probvec.size(), 2u)
        << "Spectator probvec should keep both entries (only the values change)";
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[0], 0.0)
        << "Spectator probability against the facilitator should be zeroed";
    EXPECT_DOUBLE_EQ(sys.moleculeList[2].probvec[1], 0.0)
        << "Spectator probability against the state changer should be zeroed";

    // The spectator keeps its crossbase (it must still avoid overlap).
    EXPECT_EQ(sys.moleculeList[2].crossbase.size(), 2u)
        << "Spectator crossbase should not be cleared";
    EXPECT_EQ(sys.complexList[2].ncross, 2)
        << "Spectator complex ncross should be untouched";

    // The reactants are removed from further consideration this step.
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty())
        << "Facilitator crossbase should be cleared after the reaction";
    EXPECT_TRUE(sys.moleculeList[1].crossbase.empty())
        << "State changer crossbase should be cleared after the reaction";
    EXPECT_EQ(sys.complexList[0].ncross, -1) << "Facilitator complex ncross should be -1";
    EXPECT_EQ(sys.complexList[1].ncross, -1) << "State changer complex ncross should be -1";

    // Sanity: the reaction itself succeeded.
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].stateIden, 'P')
        << "The state change should have been applied";
    EXPECT_NEAR(pbscb_iface_separation(sys), kPbscbBindRadius, 1e-9)
        << "Interfaces should be at contact after the successful reaction";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each scenario is reported separately, and no fatal
// assertions are used at file scope so every scenario runs.
// -----------------------------------------------------------------------------
TEST(PerformBimolecularStateChangeBox, ForwardPointMolecules) { test_pbscb_forward_point_molecules(); }
TEST(PerformBimolecularStateChangeBox, BackReactionStateChange) { test_pbscb_back_reaction_state_change(); }
TEST(PerformBimolecularStateChangeBox, CancelledByOverlap) { test_pbscb_cancelled_by_overlap(); }
TEST(PerformBimolecularStateChangeBox, RotationPathNonPoint) { test_pbscb_rotation_path_non_point(); }
TEST(PerformBimolecularStateChangeBox, MissingObservableLabel) { test_pbscb_missing_observable_label(); }
TEST(PerformBimolecularStateChangeBox, ZeroesPartnerProbabilities) { test_pbscb_zeroes_partner_probabilities(); }