/*! \file test_perform_implicitlipid_state_change_sphere.cpp
 *
 * ### Unit test for src/reactions/perform_implicitlipid_state_change_sphere.cpp
 *
 * The file under test contains exactly one function:
 *
 * \code
 * void perform_implicitlipid_state_change_sphere(int stateChangeIface, int facilitatorIface,
 *      std::array<int,3>& rxnItr, Molecule& stateChangeMol, Molecule& facilitatorMol,
 *      Complex& stateChangeCom, Complex& facilitatorCom, copyCounters& counterArrays,
 *      const Parameters& params, std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
 *      std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
 *      std::vector<MolTemplate>& molTemplateList, std::map<std::string,int>& observablesList,
 *      Membrane& membraneObject);
 * \endcode
 *
 * It performs a *bimolecular state change* between a soluble "facilitator" molecule and the
 * implicit lipid, for a system whose boundary is a sphere.  The routine
 *
 *   1. looks up the reflecting distance RS3D in `membraneObject.RS3Dvect`,
 *   2. creates a position for the implicit lipid next to the facilitator
 *      (`Molecule::create_position_implicit_lipid`),
 *   3. builds temporary association coordinates for both complexes, translates the facilitator
 *      to sigma (only when it is *not* already on the surface), optionally rotates the
 *      molecules into the association angles (skipped when the facilitator is a point),
 *   4. re-orients the lipid so that it is perpendicular to the sphere and reflects the
 *      facilitator back inside the sphere,
 *   5. cancels the whole event if a structural overlap is detected, otherwise commits the
 *      temporary coordinates, updates `counterArrays.copyNumSpecies`,
 *      `membraneObject.numberOfFreeLipidsEachState`, the observables map and the book-keeping
 *      fields (`ncross`, `crossbase`, `trajStatus`).
 *
 * Because the routine writes into references that must alias entries of `moleculeList`
 * (the temporary association coordinates are filled through `moleculeList`), the tests below
 * always pass `moleculeList[i]` / `complexList[i]` as the molecule / complex arguments.
 *
 * The tests exercise:
 *   - the forward-reaction branch (rxnItr[2] == 0),
 *   - the back-reaction branch (rxnItr[2] == 1, which resolves the conjugate forward reaction),
 *   - the "already on the membrane" branch of the move-to-sigma block,
 *   - the cancel-on-overlap early return.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace -> no symbol collisions with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! Geometry / reaction constants shared by every test in this file. */
constexpr double kIlscsSphereR { 100.0 }; //!< radius of the spherical boundary, nm
constexpr double kIlscsBindRadius { 1.0 }; //!< sigma of the reaction, nm
constexpr double kIlscsRate { 100.0 }; //!< micro on-rate stored in rateList[0]
constexpr double kIlscsRS3D { 0.5 }; //!< value stored in the RS3D look-up table
constexpr double kIlscsLipidLen { 1.0 }; //!< |COM->interface| of the implicit lipid

/*! \brief Everything the function under test needs, kept together so that references
 *         into the vectors stay valid for the whole test (never resize after taking them!).
 */
struct IlscsSystem {
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

/*! \brief Initialise the global GSL random number generator once.
 *
 * `r` is defined (as nullptr) in gtest_main.cpp; several routines reachable from the
 * function under test may draw random numbers, so it must point at a live generator.
 */
void ilscs_init_rng()
{
    if (r != nullptr)
        return;
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief Distance of a coordinate from the origin (centre of the spherical boundary). */
double ilscs_radius(const Coord& c) { return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z); }

/*! \brief Reproduce, verbatim, the Dtot expression the source uses for the RS3D look-up.
 *
 * The source compares `RS3Dvect[i+200]` against
 * 1/3*(Dx1+Dx2) + 1/3*(Dy1+Dy2) + 1/3*(Dz1+Dz2) of the two reactant templates.
 * Storing the same number lets the look-up succeed; if it ever failed the routine would
 * simply keep RS3D = -1.0, so no assertion in this file depends on the exact value.
 */
double ilscs_expected_dtot(const std::vector<MolTemplate>& t)
{
    return 1.0 / 3.0 * (t[0].D.x + t[1].D.x) + 1.0 / 3.0 * (t[0].D.y + t[1].D.y)
        + 1.0 / 3.0 * (t[0].D.z + t[1].D.z);
}

/*! \brief Build a complete, self-consistent two-species system on a spherical membrane.
 *
 * Species layout
 *   - molTemplate 0 "A"  : soluble point protein, one interface with states {abs 0 ('U'), abs 1 ('P')}
 *   - molTemplate 1 "IL" : implicit lipid, one interface with states {abs 2 ('U'), abs 3 ('P')}
 *
 * Reaction 0 (forward, bimolecular state change): A(a~U) + IL(il~U) -> A(a~P) + IL(il~P)
 *   reactantListNew = [absIface 0, absIface 2], productListNew = [absIface 1, absIface 3].
 *
 * \param[in] forceOverlapCancel when true the overlap machinery is switched on
 *            (checkOverlap == true on every template, an absurdly large overlapSepLimit and an
 *            extra bystander molecule sitting on top of the reactants) so that the routine is
 *            guaranteed to take the "cancel association" early return.
 */
IlscsSystem ilscs_build_system(bool forceOverlapCancel)
{
    IlscsSystem s;

    // ---------------- simulation parameters ----------------
    s.params.timeStep = 0.1;
    s.params.numMolTypes = 2;
    s.params.numTotalSpecies = 4;
    // Any COM-COM distance below overlapSepLimit cancels the event (only relevant when
    // the templates have checkOverlap == true).
    s.params.overlapSepLimit = forceOverlapCancel ? 1.0e6 : 0.1;
    s.params.scaleMaxDisplace = 1.0e12;

    // ---------------- molecule templates ----------------
    // Template 0: the facilitator.  isPoint == true makes the routine skip the
    // theta/phi/omega rotations, so its interface is genuinely placed at the COM.
    MolTemplate prot {};
    prot.molTypeIndex = 0;
    prot.molName = "A";
    prot.copies = 1;
    prot.mass = 1.0;
    prot.radius = 1.0;
    prot.D = Coord { 3.0, 3.0, 3.0 };
    prot.Dr = Coord { 0.1, 0.1, 0.1 };
    prot.isPoint = true;
    prot.isRod = true;
    prot.isLipid = false;
    prot.isImplicitLipid = false;
    prot.checkOverlap = forceOverlapCancel;
    {
        Interface iface {};
        iface.index = 0;
        iface.name = "a";
        iface.iCoord = Coord { 0.0, 0.0, 0.0 }; // point molecule: interface sits on the COM
        iface.stateList.emplace_back(std::string("a~U"), 'U', 0);
        iface.stateList.emplace_back(std::string("a~P"), 'P', 1);
        prot.interfaceList.push_back(iface);
    }
    s.molTemplateList.push_back(prot);

    // Template 1: the implicit lipid.  One interface -> a "rod", which keeps
    // save_mem_orientation() on its single-vector path.
    MolTemplate lipid {};
    lipid.molTypeIndex = 1;
    lipid.molName = "IL";
    lipid.copies = 1;
    lipid.mass = 1.0;
    lipid.radius = 1.0;
    lipid.D = Coord { 0.0, 0.0, 0.0 };
    lipid.Dr = Coord { 0.0, 0.0, 0.0 };
    lipid.isPoint = false;
    lipid.isRod = true;
    lipid.isLipid = true;
    lipid.isImplicitLipid = true;
    lipid.checkOverlap = forceOverlapCancel;
    {
        Interface iface {};
        iface.index = 0;
        iface.name = "il";
        iface.iCoord = Coord { 0.0, 0.0, kIlscsLipidLen };
        iface.stateList.emplace_back(std::string("il~U"), 'U', 2);
        iface.stateList.emplace_back(std::string("il~P"), 'P', 3);
        lipid.interfaceList.push_back(iface);
    }
    s.molTemplateList.push_back(lipid);

    // ---------------- membrane / boundary ----------------
    s.membrane.isSphere = true;
    s.membrane.isBox = false;
    s.membrane.sphereR = kIlscsSphereR;
    s.membrane.implicitLipid = true;
    s.membrane.implicitlipidIndex = 1;
    s.membrane.nStates = 2;
    s.membrane.lipidLength = kIlscsLipidLen;
    s.membrane.waterBox = Membrane::WaterBox(
        std::vector<double> { 2.0 * kIlscsSphereR, 2.0 * kIlscsSphereR, 2.0 * kIlscsSphereR });
    // Two lipid states: [0] == 'U' (abs 2), [1] == 'P' (abs 3)
    s.membrane.numberOfFreeLipidsEachState = std::vector<int> { 1000, 7 };
    s.membrane.numberOfProteinEachState = std::vector<int> { 1, 0 };

    // RS3D look-up table: 5 blocks of 100 entries (bindRadius | rate | Dtot | RS3D | molType).
    s.membrane.RS3Dvect.assign(500, 0.0);
    s.membrane.RS3Dvect[0] = kIlscsBindRadius;
    s.membrane.RS3Dvect[100] = kIlscsRate;
    s.membrane.RS3Dvect[200] = ilscs_expected_dtot(s.molTemplateList);
    s.membrane.RS3Dvect[300] = kIlscsRS3D;
    s.membrane.RS3Dvect[400] = 1.0; // molTypeIndex of the implicit lipid

    // ---------------- molecules ----------------
    // molecule 0 : facilitator, in solution, 10 nm inside the spherical membrane
    Molecule facil {};
    facil.index = 0;
    facil.id = 0;
    facil.molTypeIndex = 0;
    facil.myComIndex = 0;
    facil.mass = 1.0;
    facil.isEmpty = false;
    facil.isLipid = false;
    facil.isImplicitLipid = false;
    facil.trajStatus = TrajStatus::none;
    facil.comCoord = Coord { kIlscsSphereR - 10.0, 0.0, 0.0 };
    {
        Molecule::Iface i {};
        i.coord = facil.comCoord; // point molecule
        i.index = 0; // absolute index of state 'U'
        i.relIndex = 0;
        i.molTypeIndex = 0;
        i.stateIndex = 0;
        i.stateIden = 'U';
        i.isBound = false;
        facil.interfaceList.push_back(i);
    }
    facil.freelist.push_back(0);
    s.moleculeList.push_back(facil);

    // molecule 1 : the implicit lipid.  Its coordinates are overwritten by
    // create_position_implicit_lipid(); only |COM->iface| (== 1 nm) matters here.
    Molecule lip {};
    lip.index = 1;
    lip.id = 1;
    lip.molTypeIndex = 1;
    lip.myComIndex = 1;
    lip.mass = 1.0;
    lip.isEmpty = false;
    lip.isLipid = true;
    lip.isImplicitLipid = true;
    lip.trajStatus = TrajStatus::none;
    lip.comCoord = Coord { 0.0, 0.0, -(kIlscsSphereR + kIlscsLipidLen) };
    {
        Molecule::Iface i {};
        i.coord = Coord { 0.0, 0.0, -kIlscsSphereR };
        i.index = 2; // absolute index of lipid state 'U'
        i.relIndex = 0;
        i.molTypeIndex = 1;
        i.stateIndex = 0;
        i.stateIden = 'U';
        i.isBound = false;
        lip.interfaceList.push_back(i);
    }
    lip.freelist.push_back(0);
    s.moleculeList.push_back(lip);

    // molecule 2 (only for the cancel test): a bystander protein sitting essentially on top of
    // the facilitator so that check_for_structure_overlap_system() also flags an overlap.
    if (forceOverlapCancel) {
        Molecule bystander {};
        bystander.index = 2;
        bystander.id = 2;
        bystander.molTypeIndex = 0;
        bystander.myComIndex = 2;
        bystander.mass = 1.0;
        bystander.isEmpty = false;
        bystander.trajStatus = TrajStatus::none;
        bystander.comCoord = Coord { kIlscsSphereR - 12.0, 0.0, 0.0 };
        Molecule::Iface i {};
        i.coord = bystander.comCoord;
        i.index = 0;
        i.relIndex = 0;
        i.molTypeIndex = 0;
        i.stateIndex = 0;
        i.stateIden = 'U';
        bystander.interfaceList.push_back(i);
        bystander.freelist.push_back(0);
        s.moleculeList.push_back(bystander);
    }

    // ---------------- complexes ----------------
    Complex facilCom {};
    facilCom.index = 0;
    facilCom.id = 0;
    facilCom.comCoord = s.moleculeList[0].comCoord;
    facilCom.memberList = std::vector<int> { 0 };
    facilCom.numEachMol = std::vector<int> { 1, 0 };
    facilCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    facilCom.mass = 1.0;
    facilCom.radius = 1.0;
    facilCom.D = Coord { 3.0, 3.0, 3.0 };
    facilCom.Dr = Coord { 0.1, 0.1, 0.1 };
    facilCom.OnSurface = false; // in solution -> the 3D->2D "transition" branch
    facilCom.ncross = 0;
    facilCom.trajStatus = TrajStatus::none;
    s.complexList.push_back(facilCom);

    Complex lipCom {};
    lipCom.index = 1;
    lipCom.id = 1;
    lipCom.comCoord = s.moleculeList[1].comCoord;
    lipCom.memberList = std::vector<int> { 1 };
    lipCom.numEachMol = std::vector<int> { 0, 1 };
    lipCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
    lipCom.mass = 1.0;
    lipCom.radius = 1.0;
    lipCom.D = Coord { 0.0, 0.0, 0.0 };
    lipCom.Dr = Coord { 0.0, 0.0, 0.0 };
    lipCom.OnSurface = true;
    lipCom.iLipidIndex = 1;
    lipCom.ncross = 0;
    lipCom.trajStatus = TrajStatus::none;
    s.complexList.push_back(lipCom);

    if (forceOverlapCancel) {
        Complex byCom {};
        byCom.index = 2;
        byCom.id = 2;
        byCom.comCoord = s.moleculeList[2].comCoord;
        byCom.memberList = std::vector<int> { 2 };
        byCom.numEachMol = std::vector<int> { 1, 0 };
        byCom.lastNumberUpdateItrEachMol = std::vector<long long int> { 0, 0 };
        byCom.mass = 1.0;
        byCom.radius = 1.0;
        byCom.D = Coord { 3.0, 3.0, 3.0 };
        byCom.Dr = Coord { 0.1, 0.1, 0.1 };
        byCom.ncross = 0;
        s.complexList.push_back(byCom);
    }

    // ---------------- reactions ----------------
    ForwardRxn fr {};
    fr.rxnType = ReactionType::biMolStateChange;
    fr.absRxnIndex = 0;
    fr.relRxnIndex = 0;
    fr.isReversible = true;
    fr.conjBackRxnIndex = 0;
    fr.hasStateChange = true;
    fr.isOnMem = false;
    fr.bindRadius = kIlscsBindRadius;
    fr.norm1 = Vector { 0.0, 0.0, 1.0 };
    fr.norm2 = Vector { 0.0, 0.0, 1.0 };
    // The facilitator is a point, so the angles are never used; give them finite values anyway.
    fr.assocAngles = ForwardRxn::Angles(M_PI, M_PI, M_PI, M_PI, M_PI);
    fr.isObserved = true;
    fr.observeLabel = "ILstateChange";
    fr.rateList.emplace_back();
    fr.rateList.back().rate = kIlscsRate;
    fr.reactantListNew.emplace_back("a", 0, 0, 0, 'U', false);
    fr.reactantListNew.emplace_back("il", 1, 2, 0, 'U', false);
    fr.productListNew.emplace_back("a", 0, 1, 0, 'P', false);
    fr.productListNew.emplace_back("il", 1, 3, 0, 'P', false);
    s.forwardRxns.push_back(fr);

    BackRxn br {};
    br.rxnType = ReactionType::biMolStateChange;
    br.absRxnIndex = 1;
    br.relRxnIndex = 0;
    br.conjForwardRxnIndex = 0; // resolves back to forwardRxns[0]
    br.hasStateChange = true;
    br.isObserved = true;
    br.observeLabel = "ILstateChangeBack";
    br.rateList.emplace_back();
    br.rateList.back().rate = 1.0;
    br.reactantListNew = s.forwardRxns[0].productListNew;
    br.productListNew = s.forwardRxns[0].reactantListNew;
    s.backRxns.push_back(br);

    // ---------------- counters / observables ----------------
    // absolute interface-state indices 0..3 -> four species
    s.counterArrays.copyNumSpecies = std::vector<int> { 10, 0, 1000, 0 };
    s.observablesList["ILstateChange"] = 0;
    s.observablesList["ILstateChangeBack"] = 5;

    return s;
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: forward reaction, no overlap -> the full state change is committed.
// -----------------------------------------------------------------------------
void test_ilscs_forward_state_change_commits()
{
    std::cerr << "\n[TEST] test_ilscs_forward_state_change_commits\n"
              << "  Source file:   src/reactions/perform_implicitlipid_state_change_sphere.cpp\n"
              << "  Function:      perform_implicitlipid_state_change_sphere()\n"
              << "  Scenario:      a soluble point protein at r=90 nm changes the state of the\n"
              << "                 implicit lipid on a sphere of radius 100 nm; overlap checking\n"
              << "                 is disabled so the event is NOT cancelled.\n"
              << "  Pass criteria: copyNumSpecies[0]-1, [1]+1, [2]-1, [3]+1;\n"
              << "                 free-lipid state counts move one lipid from state 0 to state 1;\n"
              << "                 the observable is incremented; temporary association coords are\n"
              << "                 cleared; the facilitator stays inside the sphere and its complex\n"
              << "                 book-keeping (ncross/crossbase/trajStatus) is reset.\n";

    ilscs_init_rng();
    IlscsSystem sys = ilscs_build_system(/*forceOverlapCancel=*/false);

    // rxnItr = {rxnIndex, rateIndex, isStateChangeBackRxn}
    std::array<int, 3> rxnItr { 0, 0, 0 };

    // NOTE: the routine fills the temporary association coordinates through moleculeList, so the
    // Molecule/Complex references handed to it must alias the container entries.
    Molecule& facilMol = sys.moleculeList[0];
    Molecule& lipidMol = sys.moleculeList[1];
    Complex& facilCom = sys.complexList[0];
    Complex& lipidCom = sys.complexList[1];

    const std::vector<int> copyBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> freeBefore = sys.membrane.numberOfFreeLipidsEachState;
    const Coord facilPosBefore = facilMol.comCoord;

    std::cerr << "  Before: copyNumSpecies = [" << copyBefore[0] << ", " << copyBefore[1] << ", "
              << copyBefore[2] << ", " << copyBefore[3] << "], freeLipids = [" << freeBefore[0]
              << ", " << freeBefore[1] << "]\n";
    std::cerr << "  Calling perform_implicitlipid_state_change_sphere...\n";

    perform_implicitlipid_state_change_sphere(/*stateChangeIface=*/0, /*facilitatorIface=*/0, rxnItr,
        lipidMol, facilMol, lipidCom, facilCom, sys.counterArrays, sys.params, sys.forwardRxns,
        sys.backRxns, sys.moleculeList, sys.complexList, sys.molTemplateList, sys.observablesList,
        sys.membrane);

    std::cerr << "  After : copyNumSpecies = [" << sys.counterArrays.copyNumSpecies[0] << ", "
              << sys.counterArrays.copyNumSpecies[1] << ", " << sys.counterArrays.copyNumSpecies[2]
              << ", " << sys.counterArrays.copyNumSpecies[3] << "], freeLipids = ["
              << sys.membrane.numberOfFreeLipidsEachState[0] << ", "
              << sys.membrane.numberOfFreeLipidsEachState[1] << "]\n";

    // --- species copy numbers: reactant states decrement, product states increment.
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], copyBefore[0] - 1)
        << "reactant state A(a~U) (abs index 0) should lose one copy";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], copyBefore[1] + 1)
        << "product state A(a~P) (abs index 1) should gain one copy";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[2], copyBefore[2] - 1)
        << "reactant state IL(il~U) (abs index 2) should lose one copy";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[3], copyBefore[3] + 1)
        << "product state IL(il~P) (abs index 3) should gain one copy";

    // --- implicit lipid book-keeping: one free lipid moves from state 0 to state 1.
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], freeBefore[0] - 1)
        << "one free lipid should leave state 0 ('U')";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[1], freeBefore[1] + 1)
        << "one free lipid should enter state 1 ('P')";

    // --- observable of the forward reaction is incremented.
    EXPECT_EQ(sys.observablesList["ILstateChange"], 1)
        << "forward reaction observable should have been incremented once";
    EXPECT_EQ(sys.observablesList["ILstateChangeBack"], 5)
        << "the back-reaction observable must not be touched by a forward event";

    // --- temporary association coordinates must have been released again.
    EXPECT_TRUE(facilMol.tmpICoords.empty()) << "facilitator tmpICoords should be cleared";
    EXPECT_TRUE(lipidMol.tmpICoords.empty()) << "implicit lipid tmpICoords should be cleared";

    // --- the facilitator was committed from its temporary coordinates and flagged as moved.
    EXPECT_EQ(static_cast<int>(facilMol.trajStatus), static_cast<int>(TrajStatus::propagated))
        << "the facilitator should be marked as propagated";
    const double facilR = ilscs_radius(facilMol.comCoord);
    std::cerr << "  Facilitator moved from r=" << ilscs_radius(facilPosBefore) << " to r=" << facilR
              << " (sphere radius " << kIlscsSphereR << ")\n";
    EXPECT_LE(facilR, kIlscsSphereR + 1e-6)
        << "boundary conditions must keep the facilitator inside the sphere";
    EXPECT_GT(facilR, ilscs_radius(facilPosBefore))
        << "the facilitator should have been pulled out towards the lipid on the membrane";

    // --- the implicit lipid was repositioned by create_position_implicit_lipid():
    //     |COM| = sphereR + |COM->iface| + bindRadius (it is placed just outside the sphere).
    const double expectedLipidR = kIlscsSphereR + kIlscsLipidLen + kIlscsBindRadius;
    std::cerr << "  Implicit lipid COM radius = " << ilscs_radius(lipidMol.comCoord)
              << " (expected " << expectedLipidR << ")\n";
    EXPECT_NEAR(ilscs_radius(lipidMol.comCoord), expectedLipidR, 1e-6)
        << "the implicit lipid should be placed just outside the sphere along the facilitator";

    // --- the implicit lipid complex is reduced to the single lipid molecule.
    ASSERT_EQ(lipidCom.memberList.size(), 1u)
        << "the implicit lipid complex must contain exactly one member";
    EXPECT_EQ(lipidCom.memberList[0], lipidMol.index)
        << "the only member of the lipid complex should be the implicit lipid itself";

    // --- crossing book-keeping is reset so neither complex reacts again this step.
    EXPECT_EQ(facilCom.ncross, -1) << "facilitator complex ncross should be set to -1";
    EXPECT_EQ(lipidCom.ncross, -1) << "lipid complex ncross should be set to -1";
    EXPECT_TRUE(facilMol.crossbase.empty()) << "facilitator crossbase should be cleared";
    EXPECT_TRUE(lipidMol.crossbase.empty()) << "implicit lipid crossbase should be cleared";
}

// -----------------------------------------------------------------------------
// Test 2: back reaction (rxnItr[2] == 1) -> the conjugate forward reaction is used
//         for geometry and counters, the *back* observable is decremented.
// -----------------------------------------------------------------------------
void test_ilscs_back_reaction_uses_conjugate_forward()
{
    std::cerr << "\n[TEST] test_ilscs_back_reaction_uses_conjugate_forward\n"
              << "  Source file:   src/reactions/perform_implicitlipid_state_change_sphere.cpp\n"
              << "  Function:      perform_implicitlipid_state_change_sphere()\n"
              << "  Scenario:      same system, but rxnItr[2] == 1 so the reaction index refers to\n"
              << "                 backRxns[0], whose conjForwardRxnIndex points at forwardRxns[0].\n"
              << "  Pass criteria: the geometry/counter updates come from the conjugate FORWARD\n"
              << "                 reaction (identical copy-number arithmetic as test 1) while the\n"
              << "                 BACK reaction observable is decremented instead of incremented.\n";

    ilscs_init_rng();
    IlscsSystem sys = ilscs_build_system(/*forceOverlapCancel=*/false);

    std::array<int, 3> rxnItr { 0, 0, 1 }; // index 0 of backRxns, flagged as a back reaction

    Molecule& facilMol = sys.moleculeList[0];
    Molecule& lipidMol = sys.moleculeList[1];
    Complex& facilCom = sys.complexList[0];
    Complex& lipidCom = sys.complexList[1];

    const std::vector<int> copyBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> freeBefore = sys.membrane.numberOfFreeLipidsEachState;

    std::cerr << "  Calling perform_implicitlipid_state_change_sphere (back-reaction branch)...\n";
    perform_implicitlipid_state_change_sphere(0, 0, rxnItr, lipidMol, facilMol, lipidCom, facilCom,
        sys.counterArrays, sys.params, sys.forwardRxns, sys.backRxns, sys.moleculeList,
        sys.complexList, sys.molTemplateList, sys.observablesList, sys.membrane);

    std::cerr << "  After : copyNumSpecies = [" << sys.counterArrays.copyNumSpecies[0] << ", "
              << sys.counterArrays.copyNumSpecies[1] << ", " << sys.counterArrays.copyNumSpecies[2]
              << ", " << sys.counterArrays.copyNumSpecies[3] << "], observable(back) = "
              << sys.observablesList["ILstateChangeBack"] << '\n';

    // The source always updates copyNumSpecies with the *forward* reaction's lists.
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], copyBefore[0] - 1)
        << "forward reactant list is used for the counters even on the back branch";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], copyBefore[1] + 1)
        << "forward product list is used for the counters even on the back branch";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[2], copyBefore[2] - 1)
        << "implicit lipid reactant state should lose one copy";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[3], copyBefore[3] + 1)
        << "implicit lipid product state should gain one copy";

    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[0], freeBefore[0] - 1)
        << "free lipids of state 0 should be decremented";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState[1], freeBefore[1] + 1)
        << "free lipids of state 1 should be incremented";

    // Observables: back reactions decrement their own label and leave the forward label alone.
    EXPECT_EQ(sys.observablesList["ILstateChangeBack"], 4)
        << "the back-reaction observable should have been decremented from 5 to 4";
    EXPECT_EQ(sys.observablesList["ILstateChange"], 0)
        << "the forward observable must be untouched on a back reaction";

    // The geometric part still ran: the lipid was repositioned and coordinates were released.
    EXPECT_NEAR(ilscs_radius(lipidMol.comCoord), kIlscsSphereR + kIlscsLipidLen + kIlscsBindRadius,
        1e-6)
        << "the conjugate forward bindRadius should drive create_position_implicit_lipid()";
    EXPECT_TRUE(facilMol.tmpICoords.empty()) << "facilitator tmpICoords should be cleared";
    EXPECT_TRUE(lipidMol.tmpICoords.empty()) << "implicit lipid tmpICoords should be cleared";
}

// -----------------------------------------------------------------------------
// Test 3: the facilitator complex is already on the membrane (OnSurface == true),
//         so the "move to sigma" translation is skipped.
// -----------------------------------------------------------------------------
void test_ilscs_facilitator_already_on_surface()
{
    std::cerr << "\n[TEST] test_ilscs_facilitator_already_on_surface\n"
              << "  Source file:   src/reactions/perform_implicitlipid_state_change_sphere.cpp\n"
              << "  Function:      perform_implicitlipid_state_change_sphere()\n"
              << "  Scenario:      facilitatorCom.OnSurface == true and the facilitator already sits\n"
              << "                 on the sphere, so the 3D->2D translation to sigma is skipped\n"
              << "                 (isOnMembrane branch).\n"
              << "  Pass criteria: the state change is still committed (counters/observables) and\n"
              << "                 the facilitator remains on/inside the spherical boundary.\n";

    ilscs_init_rng();
    IlscsSystem sys = ilscs_build_system(/*forceOverlapCancel=*/false);

    // Put the facilitator directly on the membrane and flag its complex as surface bound.
    sys.moleculeList[0].comCoord = Coord { kIlscsSphereR, 0.0, 0.0 };
    sys.moleculeList[0].interfaceList[0].coord = sys.moleculeList[0].comCoord;
    sys.complexList[0].comCoord = sys.moleculeList[0].comCoord;
    sys.complexList[0].OnSurface = true;
    sys.complexList[0].D = Coord { 3.0, 3.0, 0.0 };

    std::array<int, 3> rxnItr { 0, 0, 0 };
    Molecule& facilMol = sys.moleculeList[0];
    Molecule& lipidMol = sys.moleculeList[1];
    Complex& facilCom = sys.complexList[0];
    Complex& lipidCom = sys.complexList[1];

    const std::vector<int> copyBefore = sys.counterArrays.copyNumSpecies;

    std::cerr << "  Calling perform_implicitlipid_state_change_sphere (isOnMembrane branch)...\n";
    perform_implicitlipid_state_change_sphere(0, 0, rxnItr, lipidMol, facilMol, lipidCom, facilCom,
        sys.counterArrays, sys.params, sys.forwardRxns, sys.backRxns, sys.moleculeList,
        sys.complexList, sys.molTemplateList, sys.observablesList, sys.membrane);

    const double facilR = ilscs_radius(facilMol.comCoord);
    std::cerr << "  Facilitator radius after the call = " << facilR << " (sphere radius "
              << kIlscsSphereR << ")\n";

    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], copyBefore[0] - 1)
        << "the state change should still be committed when starting on the membrane";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[3], copyBefore[3] + 1)
        << "the lipid product state should still be incremented";
    EXPECT_EQ(sys.observablesList["ILstateChange"], 1)
        << "the forward observable should be incremented once";
    EXPECT_LE(facilR, kIlscsSphereR + 1e-6)
        << "the facilitator must not be left outside the spherical boundary";
    EXPECT_TRUE(facilMol.tmpICoords.empty())
        << "temporary association coordinates must be released";
    EXPECT_EQ(facilCom.ncross, -1) << "facilitator complex ncross should be set to -1";
}

// -----------------------------------------------------------------------------
// Test 4: overlap detected -> association is cancelled and nothing is committed.
// -----------------------------------------------------------------------------
void test_ilscs_cancel_on_structure_overlap()
{
    std::cerr << "\n[TEST] test_ilscs_cancel_on_structure_overlap\n"
              << "  Source file:   src/reactions/perform_implicitlipid_state_change_sphere.cpp\n"
              << "  Function:      perform_implicitlipid_state_change_sphere()\n"
              << "  Scenario:      every template has checkOverlap == true, overlapSepLimit is set to\n"
              << "                 1e6 nm and a bystander molecule sits next to the reactants, so the\n"
              << "                 structure-overlap checks must flag the event.\n"
              << "  Pass criteria: the routine returns early -> copy numbers, free-lipid counts and\n"
              << "                 observables are untouched, temporary coordinates are released, the\n"
              << "                 facilitator has NOT moved and ncross is left alone. (The implicit\n"
              << "                 lipid position is still updated: create_position_implicit_lipid()\n"
              << "                 runs before the overlap test.)\n";

    ilscs_init_rng();
    IlscsSystem sys = ilscs_build_system(/*forceOverlapCancel=*/true);

    std::array<int, 3> rxnItr { 0, 0, 0 };
    Molecule& facilMol = sys.moleculeList[0];
    Molecule& lipidMol = sys.moleculeList[1];
    Complex& facilCom = sys.complexList[0];
    Complex& lipidCom = sys.complexList[1];

    const std::vector<int> copyBefore = sys.counterArrays.copyNumSpecies;
    const std::vector<int> freeBefore = sys.membrane.numberOfFreeLipidsEachState;
    const Coord facilPosBefore = facilMol.comCoord;
    const int ncrossBefore = facilCom.ncross;

    std::cerr << "  Calling perform_implicitlipid_state_change_sphere (expected cancel)...\n";
    perform_implicitlipid_state_change_sphere(0, 0, rxnItr, lipidMol, facilMol, lipidCom, facilCom,
        sys.counterArrays, sys.params, sys.forwardRxns, sys.backRxns, sys.moleculeList,
        sys.complexList, sys.molTemplateList, sys.observablesList, sys.membrane);

    std::cerr << "  After : copyNumSpecies = [" << sys.counterArrays.copyNumSpecies[0] << ", "
              << sys.counterArrays.copyNumSpecies[1] << ", " << sys.counterArrays.copyNumSpecies[2]
              << ", " << sys.counterArrays.copyNumSpecies[3] << "], facilitator COM = ("
              << facilMol.comCoord.x << ", " << facilMol.comCoord.y << ", " << facilMol.comCoord.z
              << ")\n";

    // Nothing may be committed on the cancel path.
    EXPECT_EQ(sys.counterArrays.copyNumSpecies, copyBefore)
        << "copy numbers must be unchanged when the association is cancelled";
    EXPECT_EQ(sys.membrane.numberOfFreeLipidsEachState, freeBefore)
        << "free lipid counts must be unchanged when the association is cancelled";
    EXPECT_EQ(sys.observablesList["ILstateChange"], 0)
        << "no observable should be recorded for a cancelled event";

    // Temporary association coordinates are still released before the early return.
    EXPECT_TRUE(facilMol.tmpICoords.empty())
        << "facilitator tmpICoords should be cleared on cancel";
    EXPECT_TRUE(lipidMol.tmpICoords.empty())
        << "implicit lipid tmpICoords should be cleared on cancel";

    // Real coordinates of the facilitator are untouched (only temporaries were manipulated).
    EXPECT_DOUBLE_EQ(facilMol.comCoord.x, facilPosBefore.x)
        << "facilitator x must not move when the event is cancelled";
    EXPECT_DOUBLE_EQ(facilMol.comCoord.y, facilPosBefore.y)
        << "facilitator y must not move when the event is cancelled";
    EXPECT_DOUBLE_EQ(facilMol.comCoord.z, facilPosBefore.z)
        << "facilitator z must not move when the event is cancelled";

    // ncross/crossbase resetting happens only after a successful event.
    EXPECT_EQ(facilCom.ncross, ncrossBefore)
        << "ncross should not be reset when the association is cancelled";

    // The implicit lipid was nevertheless repositioned next to the facilitator, because
    // create_position_implicit_lipid() runs before the overlap checks.
    EXPECT_NEAR(ilscs_radius(lipidMol.comCoord), kIlscsSphereR + kIlscsLipidLen + kIlscsBindRadius,
        1e-6)
        << "the implicit lipid placement happens before the overlap test and is kept";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each scenario is reported individually but all of them run,
// because only non-fatal EXPECT_* assertions are used inside.
// -----------------------------------------------------------------------------
TEST(PerformImplicitLipidStateChangeSphere, ForwardStateChangeCommits)
{
    test_ilscs_forward_state_change_commits();
}
TEST(PerformImplicitLipidStateChangeSphere, BackReactionUsesConjugateForward)
{
    test_ilscs_back_reaction_uses_conjugate_forward();
}
TEST(PerformImplicitLipidStateChangeSphere, FacilitatorAlreadyOnSurface)
{
    test_ilscs_facilitator_already_on_surface();
}
TEST(PerformImplicitLipidStateChangeSphere, CancelOnStructureOverlap)
{
    test_ilscs_cancel_on_structure_overlap();
}