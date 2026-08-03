/*! \file test_determine_3D_implicitlipid_reaction_probability.cpp
 *
 * ### Unit test for src/reactions/determine_3D_implicitlipid_reaction_probability.cpp
 *
 * Function under test:
 *
 *     void determine_3D_implicitlipid_reaction_probability(
 *         int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
 *         BiMolData& biMolData, const Parameters& params,
 *         std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
 *         const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
 *         Membrane& membraneObject, const int& relStateIndex)
 *
 * What the routine does (and therefore what we verify):
 *
 *  1. It augments `biMolData.Dtot` with a rotational-diffusion contribution for
 *     each of the two participating complexes.  The pre-factor depends on
 *     whether the complex has a non-zero z translational diffusion constant
 *     (3D, factor 6*dt) or not (membrane bound, factor 4*dt).
 *  2. It computes Rmax = 3*sqrt(6*Dtot*dt) + bindRadius and asks
 *     get_distance_to_surface() whether the protein interface is close enough
 *     to the implicit-lipid membrane to react.
 *  3. If it is inside Rmax a *placeholder* zero probability is appended to
 *     moleculeList[pro1].probvec (note: nothing is ever appended to the
 *     implicit lipid's probvec - that line is commented out in the source).
 *  4. Only if the protein did not just dissociate AND the reaction rate is
 *     positive does it evaluate the binding probability
 *         prob = rho * pimplicitlipid_3D(z, params3D)
 *     write it into probvec.back(), and record the reweighting book-keeping
 *     lists (currprevsep, currlist, currmyface, currpface, currprevnorm,
 *     currps_prev).
 *  5. If the separation is negative (interface already inside the membrane)
 *     the probability is forced to unity.
 *
 * All tests are written with non-fatal EXPECT_* assertions and every access to
 * a possibly-empty container is guarded, so a single failure never aborts the
 * remainder of the suite.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (all prefixed with d3il_ so they cannot collide with any other
// translation unit in the combined gtest binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Everything the function under test needs, kept together so a test can
 *         build a complete, self consistent little "simulation" in one call.
 *
 * Layout of the system:
 *   - molecule 0 / complex 0 : the protein whose binding probability we want
 *   - molecule 1 / complex 1 : the implicit lipid (a stand-in for the membrane)
 *   - the membrane surface lives at the bottom of the cubic water box,
 *     i.e. at z = -waterBox.z / 2.
 */
struct D3ILScenario {
    Parameters params {};
    Membrane membrane {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    BiMolData biMolData {};

    // Constants used both when building and when checking the scenario.
    double boxLength { 100.0 }; //!< cubic water box side length, nm
    double bindRadius { 1.0 }; //!< sigma of the forward reaction, nm
    double dtotTranslational { 10.0 }; //!< translational Dtot handed in, nm^2/us
    double magMol1 { 1.0 }; //!< |COM -> iface| of the protein, nm
    double magMol2 { 0.0 }; //!< |COM -> iface| of the implicit lipid, nm
};

/*! \brief Reproduces (independently of the implementation) the rotational
 *         diffusion contribution that the function adds to Dtot.
 *
 * \param[in] Dz  the complex' z translational diffusion constant
 * \param[in] Drz the complex' z rotational diffusion constant
 * \param[in] mag the COM-to-interface vector magnitude for that complex
 * \param[in] dt  the simulation time step
 */
double d3il_rot_contribution(double Dz, double Drz, double mag, double dt)
{
    if (std::abs(Dz - 0.0) < 1E-15) {
        // Membrane bound complex: only 2 rotational degrees of freedom matter.
        const double cf = std::cos(std::sqrt(2.0 * Drz * dt));
        return (2.0 * mag * (1.0 - cf)) / (4.0 * dt);
    }
    // Freely diffusing (3D) complex.
    const double cf = std::cos(std::sqrt(4.0 * Drz * dt));
    return (2.0 * mag * (1.0 - cf)) / (6.0 * dt);
}

/*! \brief Builds a minimal but complete scenario.
 *
 * \param[in] ifaceZ            z coordinate of the protein's reacting interface
 * \param[in] rate              forward (intrinsic) association rate
 * \param[in] nFreeLipids       number of free implicit lipids of the state
 * \param[in] proteinOnMembrane if true the protein complex has D.z == 0, which
 *                              selects the "membrane" branch of the Dtot update
 */
D3ILScenario d3il_make_scenario(double ifaceZ, double rate, int nFreeLipids, bool proteinOnMembrane = false)
{
    D3ILScenario s;

    /* ---- simulation parameters --------------------------------------- */
    s.params.timeStep = 0.1; // us; small so Rmax stays small
    s.params.numMolTypes = 2;
    s.params.numTotalSpecies = 2;

    /* ---- membrane / boundary ----------------------------------------- */
    s.membrane.isBox = true;
    s.membrane.isSphere = false;
    s.membrane.implicitLipid = true;
    s.membrane.waterBox = Membrane::WaterBox(std::vector<double> { s.boxLength, s.boxLength, s.boxLength });
    s.membrane.totalSA = s.boxLength * s.boxLength; // area of the flat membrane
    s.membrane.nStates = 1;
    s.membrane.numberOfFreeLipidsEachState = std::vector<int> { nFreeLipids };
    s.membrane.numberOfProteinEachState = std::vector<int> { 1 };
    s.membrane.No_free_lipids = nFreeLipids;
    s.membrane.No_protein = 1;
    s.membrane.implicitlipidIndex = 1;
    s.membrane.xBCtype = "reflect";
    s.membrane.yBCtype = "reflect";
    s.membrane.zBCtype = "reflect";

    /* ---- molecule 0: the protein ------------------------------------- */
    Molecule protein;
    protein.index = 0;
    protein.myComIndex = 0;
    protein.molTypeIndex = 0;
    protein.mass = 1.0;
    protein.isImplicitLipid = false;
    protein.isLipid = false;
    protein.trajStatus = TrajStatus::none;
    protein.isDissociated = false;
    // COM sits |magMol1| above its own interface so that magMol1 is consistent.
    protein.comCoord = Coord(0.0, 0.0, ifaceZ + s.magMol1);
    Molecule::Iface pIface;
    pIface.coord = Coord(0.0, 0.0, ifaceZ);
    pIface.index = 0;
    pIface.relIndex = 0;
    pIface.molTypeIndex = 0;
    pIface.isBound = false;
    protein.interfaceList.push_back(pIface);
    protein.freelist.push_back(0);
    s.moleculeList.push_back(protein);

    /* ---- molecule 1: the implicit lipid ------------------------------ */
    Molecule lipid;
    lipid.index = 1;
    lipid.myComIndex = 1;
    lipid.molTypeIndex = 1;
    lipid.mass = 1.0;
    lipid.isImplicitLipid = true;
    lipid.isLipid = true;
    lipid.trajStatus = TrajStatus::none;
    lipid.comCoord = Coord(0.0, 0.0, -s.boxLength / 2.0);
    Molecule::Iface lIface;
    lIface.coord = Coord(0.0, 0.0, -s.boxLength / 2.0);
    lIface.index = 1;
    lIface.relIndex = 0;
    lIface.molTypeIndex = 1;
    lIface.isBound = false;
    lipid.interfaceList.push_back(lIface);
    lipid.freelist.push_back(0);
    s.moleculeList.push_back(lipid);

    /* ---- complex 0: protein complex ---------------------------------- */
    Complex proteinCom;
    proteinCom.index = 0;
    proteinCom.comCoord = s.moleculeList[0].comCoord;
    proteinCom.memberList.push_back(0);
    proteinCom.radius = 1.0;
    proteinCom.mass = 1.0;
    // D.z decides which branch of the Dtot update is exercised.
    proteinCom.D = Coord(10.0, 10.0, proteinOnMembrane ? 0.0 : 10.0);
    proteinCom.Dr = Coord(0.05, 0.05, 0.05);
    proteinCom.OnSurface = proteinOnMembrane;
    proteinCom.linksToSurface = 0;
    proteinCom.trajStatus = TrajStatus::none;
    s.complexList.push_back(proteinCom);

    /* ---- complex 1: implicit lipid complex --------------------------- */
    Complex lipidCom;
    lipidCom.index = 1;
    lipidCom.comCoord = s.moleculeList[1].comCoord;
    lipidCom.memberList.push_back(1);
    lipidCom.radius = 1.0;
    lipidCom.mass = 1.0;
    lipidCom.D = Coord(0.0, 0.0, 0.0); // implicit lipid does not move
    lipidCom.Dr = Coord(0.0, 0.0, 0.0);
    lipidCom.OnSurface = true;
    lipidCom.trajStatus = TrajStatus::none;
    s.complexList.push_back(lipidCom);

    /* ---- the forward reaction --------------------------------------- */
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.bindRadius = s.bindRadius;
    rxn.bindRadius2D = s.bindRadius;
    rxn.isOnMem = false;
    rxn.isReversible = true;
    rxn.conjBackRxnIndex = 0;
    rxn.intReactantList = std::vector<int> { 0, 1 };
    rxn.intProductList = std::vector<int> { 0 };
    rxn.reactantListNew.push_back(RxnIface("prot", 0, 0, 0, '\0', false));
    rxn.reactantListNew.push_back(RxnIface("lipid", 1, 1, 0, '\0', false));
    rxn.productListNew.push_back(RxnIface("prot", 0, 0, 0, '\0', true));
    RxnBase::RateState state;
    state.rate = rate;
    rxn.rateList.push_back(state);
    s.forwardRxns.push_back(rxn);

    /* ---- the (unused, but referenced) back reaction ------------------ */
    BackRxn back;
    back.rxnType = ReactionType::bimolecular;
    back.absRxnIndex = 0;
    back.relRxnIndex = 0;
    back.conjForwardRxnIndex = 0;
    RxnBase::RateState backState;
    backState.rate = 1.0;
    back.rateList.push_back(backState);
    s.backRxns.push_back(back);

    /* ---- the pair data --------------------------------------------- */
    s.biMolData = BiMolData(/*pro1*/ 0, /*pro2*/ 1, /*com1*/ 0, /*com2*/ 1,
        /*relIface1*/ 0, /*relIface2*/ 0, /*absIface1*/ 0, /*absIface2*/ 1,
        /*Dtot*/ s.dtotTranslational, /*magMol1*/ s.magMol1, /*magMol2*/ s.magMol2);

    return s;
}

/*! \brief Convenience wrapper that calls the function under test on a scenario. */
void d3il_run(D3ILScenario& s, int simItr = 1, int relStateIndex = 0)
{
    determine_3D_implicitlipid_reaction_probability(simItr, /*rxnIndex*/ 0, /*rateIndex*/ 0,
        /*isStateChangeBackRxn*/ false, s.biMolData, s.params, s.moleculeList, s.complexList,
        s.forwardRxns, s.backRxns, s.membrane, relStateIndex);
}

/*! \brief Dumps the reweighting book-keeping of the protein for readability. */
void d3il_report(const D3ILScenario& s)
{
    const Molecule& p = s.moleculeList[0];
    std::cerr << "    probvec.size()=" << p.probvec.size();
    if (!p.probvec.empty())
        std::cerr << " probvec.back()=" << p.probvec.back();
    std::cerr << " currlist=" << p.currlist.size() << " currmyface=" << p.currmyface.size()
              << " currpface=" << p.currpface.size() << " currprevsep=" << p.currprevsep.size()
              << " currprevnorm=" << p.currprevnorm.size() << " currps_prev=" << p.currps_prev.size()
              << " Dtot=" << s.biMolData.Dtot << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the rotational diffusion contributions added to Dtot for a freely
//         diffusing (3D) protein complex plus a stationary implicit lipid.
// -----------------------------------------------------------------------------
void test_d3il_dtot_update_3d_complex()
{
    std::cerr << "\n[TEST] test_d3il_dtot_update_3d_complex\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    protein complex has D.z != 0 (3D branch), lipid is static.\n"
              << "  Criteria:    biMolData.Dtot == Dtot_in + rot(protein) + rot(lipid),\n"
              << "               using the 1/(6*dt) pre-factor for the protein.\n";

    // Interface placed far from the membrane so the reaction branch is skipped;
    // the Dtot update however always happens.
    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ 0.0, /*rate*/ 10.0, /*nFreeLipids*/ 1000);

    const double dtotIn = s.biMolData.Dtot;
    const double dt = s.params.timeStep;
    const double expected = dtotIn
        + d3il_rot_contribution(s.complexList[0].D.z, s.complexList[0].Dr.z, s.magMol1, dt)
        + d3il_rot_contribution(s.complexList[1].D.z, s.complexList[1].Dr.z, s.magMol2, dt);

    std::cerr << "  Dtot before call = " << dtotIn << ", expected after = " << expected << '\n';
    d3il_run(s);
    std::cerr << "  Dtot after  call = " << s.biMolData.Dtot << '\n';

    EXPECT_NEAR(s.biMolData.Dtot, expected, 1e-12)
        << "Dtot must be incremented by the rotational contribution of both complexes";
    EXPECT_GT(s.biMolData.Dtot, dtotIn)
        << "A rotating complex must strictly increase the total diffusion constant";
    d3il_report(s);
}

// -----------------------------------------------------------------------------
// Test 2: same, but for a membrane bound protein complex (D.z == 0), which
//         selects the 1/(4*dt) pre-factor.
// -----------------------------------------------------------------------------
void test_d3il_dtot_update_membrane_complex()
{
    std::cerr << "\n[TEST] test_d3il_dtot_update_membrane_complex\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    protein complex has D.z == 0 (membrane branch).\n"
              << "  Criteria:    the 2D pre-factor 1/(4*dt) with cos(sqrt(2*Dr*dt)) is used.\n";

    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ 0.0, /*rate*/ 10.0, /*nFreeLipids*/ 1000,
        /*proteinOnMembrane*/ true);

    const double dtotIn = s.biMolData.Dtot;
    const double dt = s.params.timeStep;
    const double expected = dtotIn
        + d3il_rot_contribution(0.0, s.complexList[0].Dr.z, s.magMol1, dt)
        + d3il_rot_contribution(0.0, s.complexList[1].Dr.z, s.magMol2, dt);

    // Sanity check that the two branches really differ for these inputs.
    const double threeDValue = dtotIn
        + d3il_rot_contribution(10.0, s.complexList[0].Dr.z, s.magMol1, dt)
        + d3il_rot_contribution(0.0, s.complexList[1].Dr.z, s.magMol2, dt);
    std::cerr << "  Membrane-branch expectation = " << expected
              << ", 3D-branch value would be " << threeDValue << '\n';

    d3il_run(s);
    std::cerr << "  Dtot after call = " << s.biMolData.Dtot << '\n';

    EXPECT_NEAR(s.biMolData.Dtot, expected, 1e-12)
        << "For D.z == 0 the 1/(4*dt) rotational pre-factor must be used";
    EXPECT_NE(expected, threeDValue)
        << "The two branches should give different results for the chosen inputs";
    d3il_report(s);
}

// -----------------------------------------------------------------------------
// Test 3: protein interface far away from the membrane -> nothing is recorded.
// -----------------------------------------------------------------------------
void test_d3il_no_probability_when_outside_rmax()
{
    std::cerr << "\n[TEST] test_d3il_no_probability_when_outside_rmax\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    interface sits 50 nm above the membrane while\n"
              << "               Rmax = 3*sqrt(6*Dtot*dt)+sigma is only a few nm.\n"
              << "  Criteria:    no entry appended to probvec and no reweighting data.\n";

    // Membrane is at z = -50, interface at z = 0 => 50 nm separation.
    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ 0.0, /*rate*/ 10.0, /*nFreeLipids*/ 1000);
    d3il_run(s);

    const double rmax = 3.0 * std::sqrt(6.0 * s.biMolData.Dtot * s.params.timeStep) + s.bindRadius;
    std::cerr << "  Rmax computed from the updated Dtot = " << rmax
              << " nm (actual separation is 50 nm)\n";

    EXPECT_TRUE(s.moleculeList[0].probvec.empty())
        << "Outside Rmax no reaction probability slot may be created";
    EXPECT_TRUE(s.moleculeList[0].currlist.empty())
        << "Outside Rmax no partner should be recorded in currlist";
    EXPECT_TRUE(s.moleculeList[0].currprevsep.empty())
        << "Outside Rmax no separation should be recorded";
    d3il_report(s);
}

// -----------------------------------------------------------------------------
// Test 4: protein interface within Rmax -> probability plus complete
//         reweighting book-keeping is stored on the protein only.
// -----------------------------------------------------------------------------
void test_d3il_probability_recorded_when_inside_rmax()
{
    std::cerr << "\n[TEST] test_d3il_probability_recorded_when_inside_rmax\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    interface 2 nm above the membrane, sigma = 1 nm,\n"
              << "               1000 free implicit lipids, positive rate.\n"
              << "  Criteria:    probvec gains exactly one entry in [0,1]; currlist,\n"
              << "               currmyface, currpface, currprevsep, currprevnorm and\n"
              << "               currps_prev each gain exactly one, consistent entry.\n";

    // Membrane at z = -50, interface at z = -48 => R1 = 2 nm, sep = 1 nm.
    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ -48.0, /*rate*/ 10.0, /*nFreeLipids*/ 1000);
    d3il_run(s);
    d3il_report(s);

    const Molecule& prot = s.moleculeList[0];

    // (a) exactly one probability entry was created for the protein
    EXPECT_EQ(prot.probvec.size(), 1u)
        << "Being inside Rmax must append exactly one probability to the protein";

    // (b) the implicit lipid must NOT receive a probability (that push is
    //     commented out in the source under test)
    EXPECT_EQ(s.moleculeList[1].probvec.size(), 0u)
        << "The implicit lipid partner must never get a probvec entry";

    if (!prot.probvec.empty()) {
        const double prob = prot.probvec.back();
        std::cerr << "  Reaction probability = " << prob << '\n';
        EXPECT_GE(prob, 0.0) << "Probability may not be negative";
        EXPECT_LE(prob, 1.000001) << "Probability may not exceed unity";
        EXPECT_TRUE(std::isfinite(prob)) << "Probability must be a finite number";

        // (c) currps_prev must be the complement of the probability
        if (!prot.currps_prev.empty()) {
            EXPECT_NEAR(prot.currps_prev.back(), 1.0 - prob, 1e-12)
                << "currps_prev must equal 1 - prob (currnorm is 1 for implicit lipids)";
        }
    }

    // (d) all reweighting lists get exactly one element with the expected value
    EXPECT_EQ(prot.currlist.size(), 1u) << "currlist should hold the single partner index";
    if (!prot.currlist.empty())
        EXPECT_EQ(prot.currlist.back(), s.biMolData.pro2Index) << "partner index must be pro2Index";

    EXPECT_EQ(prot.currmyface.size(), 1u) << "currmyface should hold one interface index";
    if (!prot.currmyface.empty())
        EXPECT_EQ(prot.currmyface.back(), s.biMolData.relIface1) << "currmyface must be relIface1";

    EXPECT_EQ(prot.currpface.size(), 1u) << "currpface should hold one interface index";
    if (!prot.currpface.empty())
        EXPECT_EQ(prot.currpface.back(), s.biMolData.relIface2) << "currpface must be relIface2";

    EXPECT_EQ(prot.currprevnorm.size(), 1u) << "currprevnorm should hold one entry";
    if (!prot.currprevnorm.empty())
        EXPECT_DOUBLE_EQ(prot.currprevnorm.back(), 1.0)
            << "currnorm is hard-coded to 1.0 for the implicit lipid model";

    EXPECT_EQ(prot.currprevsep.size(), 1u) << "currprevsep should hold the stored distance";
    if (!prot.currprevsep.empty()) {
        std::cerr << "  Stored distance to surface R1 = " << prot.currprevsep.back() << " nm\n";
        EXPECT_GT(prot.currprevsep.back(), 0.0)
            << "The recorded distance to the membrane should be positive here";
    }
}

// -----------------------------------------------------------------------------
// Test 5: a zero reaction rate leaves the placeholder probability at zero and
//         records no reweighting data.
// -----------------------------------------------------------------------------
void test_d3il_zero_rate_gives_zero_probability()
{
    std::cerr << "\n[TEST] test_d3il_zero_rate_gives_zero_probability\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    interface inside Rmax but forward rate == 0.\n"
              << "  Criteria:    the placeholder 0 probability stays 0 and no\n"
              << "               reweighting entries are created.\n";

    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ -48.0, /*rate*/ 0.0, /*nFreeLipids*/ 1000);
    d3il_run(s);
    d3il_report(s);

    const Molecule& prot = s.moleculeList[0];
    EXPECT_EQ(prot.probvec.size(), 1u)
        << "A placeholder probability is still appended because we are inside Rmax";
    if (!prot.probvec.empty())
        EXPECT_DOUBLE_EQ(prot.probvec.back(), 0.0) << "Zero rate must leave the probability at 0";
    EXPECT_TRUE(prot.currlist.empty()) << "No reweighting data for a zero-rate reaction";
    EXPECT_TRUE(prot.currprevsep.empty()) << "No stored separation for a zero-rate reaction";
    EXPECT_TRUE(prot.currps_prev.empty()) << "No stored survival prob for a zero-rate reaction";
}

// -----------------------------------------------------------------------------
// Test 6: a molecule that just dissociated is skipped entirely.
// -----------------------------------------------------------------------------
void test_d3il_dissociated_molecule_is_skipped()
{
    std::cerr << "\n[TEST] test_d3il_dissociated_molecule_is_skipped\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    identical to the 'inside Rmax' case but the protein\n"
              << "               carries isDissociated == true.\n"
              << "  Criteria:    the placeholder 0 probability is appended but never\n"
              << "               overwritten, and no reweighting data is stored.\n";

    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ -48.0, /*rate*/ 10.0, /*nFreeLipids*/ 1000);
    s.moleculeList[0].isDissociated = true; // just came off the membrane this step
    d3il_run(s);
    d3il_report(s);

    const Molecule& prot = s.moleculeList[0];
    EXPECT_EQ(prot.probvec.size(), 1u)
        << "The placeholder entry is pushed before the isDissociated check";
    if (!prot.probvec.empty())
        EXPECT_DOUBLE_EQ(prot.probvec.back(), 0.0)
            << "A dissociated molecule must not be given a binding probability";
    EXPECT_TRUE(prot.currlist.empty()) << "No reweighting data for a dissociated molecule";
    EXPECT_TRUE(prot.currmyface.empty()) << "No reweighting data for a dissociated molecule";
    EXPECT_TRUE(prot.currps_prev.empty()) << "No reweighting data for a dissociated molecule";
}

// -----------------------------------------------------------------------------
// Test 7: without free lipids the surface density is zero, so the probability
//         must vanish.
// -----------------------------------------------------------------------------
void test_d3il_zero_free_lipids_gives_zero_probability()
{
    std::cerr << "\n[TEST] test_d3il_zero_free_lipids_gives_zero_probability\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    interface inside Rmax, positive rate, but the state\n"
              << "               has 0 free implicit lipids => rho = 0.\n"
              << "  Criteria:    prob = rho * pimplicitlipid_3D(...) == 0.\n";

    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ -48.0, /*rate*/ 10.0, /*nFreeLipids*/ 0);
    d3il_run(s);
    d3il_report(s);

    const Molecule& prot = s.moleculeList[0];
    EXPECT_EQ(prot.probvec.size(), 1u) << "Inside Rmax a probability slot is created";
    if (!prot.probvec.empty()) {
        std::cerr << "  Probability with zero free lipids = " << prot.probvec.back() << '\n';
        EXPECT_NEAR(prot.probvec.back(), 0.0, 1e-9)
            << "With zero lipid density the binding probability must be zero";
    }
    // The reweighting lists are still filled (the rate is positive).
    EXPECT_EQ(prot.currlist.size(), 1u) << "The encounter is still recorded for reweighting";
    if (!prot.currps_prev.empty())
        EXPECT_NEAR(prot.currps_prev.back(), 1.0, 1e-9)
            << "Survival probability must be 1 when the reaction probability is 0";
}

// -----------------------------------------------------------------------------
// Test 8: monotonic behaviour - a larger intrinsic rate (and a larger lipid
//         density) can never lower the binding probability.
// -----------------------------------------------------------------------------
void test_d3il_probability_monotonic_in_rate_and_density()
{
    std::cerr << "\n[TEST] test_d3il_probability_monotonic_in_rate_and_density\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    same geometry evaluated with (rate, lipids) =\n"
              << "               (5,500), (50,500) and (5,5000).\n"
              << "  Criteria:    probability grows (or at least does not shrink) with\n"
              << "               both the intrinsic rate and the lipid density.\n";

    D3ILScenario low = d3il_make_scenario(-48.0, 5.0, 500);
    D3ILScenario fastRate = d3il_make_scenario(-48.0, 50.0, 500);
    D3ILScenario denser = d3il_make_scenario(-48.0, 5.0, 5000);

    d3il_run(low);
    d3il_run(fastRate);
    d3il_run(denser);

    const bool haveAll = !low.moleculeList[0].probvec.empty()
        && !fastRate.moleculeList[0].probvec.empty() && !denser.moleculeList[0].probvec.empty();

    EXPECT_TRUE(haveAll) << "All three scenarios should be inside Rmax and produce a probability";

    if (haveAll) {
        const double pLow = low.moleculeList[0].probvec.back();
        const double pFast = fastRate.moleculeList[0].probvec.back();
        const double pDense = denser.moleculeList[0].probvec.back();
        std::cerr << "  p(rate=5,rho_low)=" << pLow << "  p(rate=50,rho_low)=" << pFast
                  << "  p(rate=5,rho_high)=" << pDense << '\n';

        EXPECT_GE(pFast, pLow - 1e-12) << "A larger intrinsic rate must not lower the probability";
        EXPECT_GE(pDense, pLow - 1e-12) << "A larger lipid density must not lower the probability";
        EXPECT_LE(pFast, 1.000001) << "Probability stays bounded by unity";
        EXPECT_LE(pDense, 1.000001) << "Probability stays bounded by unity";
    }
}

// -----------------------------------------------------------------------------
// Test 9: an interface that already overlaps the membrane (negative separation)
//         is forced to react with probability one.
// -----------------------------------------------------------------------------
void test_d3il_negative_separation_forces_unit_probability()
{
    std::cerr << "\n[TEST] test_d3il_negative_separation_forces_unit_probability\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    the interface is pushed below the membrane plane so\n"
              << "               that the separation is negative.\n"
              << "  Criteria:    the routine clamps sep to 0 and sets prob = 1.0, so\n"
              << "               currps_prev becomes 0.\n";

    // Membrane at z = -50; put the interface at z = -50.5 (inside the membrane).
    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ -50.5, /*rate*/ 10.0, /*nFreeLipids*/ 1000);
    d3il_run(s);
    d3il_report(s);

    const Molecule& prot = s.moleculeList[0];
    EXPECT_EQ(prot.probvec.size(), 1u) << "An overlapping interface is certainly within Rmax";
    if (!prot.probvec.empty()) {
        std::cerr << "  Probability for the overlapping interface = " << prot.probvec.back() << '\n';
        EXPECT_DOUBLE_EQ(prot.probvec.back(), 1.0)
            << "A negative separation must force the reaction probability to unity";
    }
    if (!prot.currps_prev.empty()) {
        EXPECT_NEAR(prot.currps_prev.back(), 0.0, 1e-12)
            << "With prob == 1 the stored survival probability must be 0";
    }
    if (!prot.currprevsep.empty()) {
        EXPECT_DOUBLE_EQ(prot.currprevsep.back(), s.bindRadius)
            << "R1 is reset to the binding radius when the separation was negative";
    }
}

// -----------------------------------------------------------------------------
// Test 10: repeated calls accumulate independent entries (the routine appends,
//          it never clears the lists).
// -----------------------------------------------------------------------------
void test_d3il_repeated_calls_accumulate_entries()
{
    std::cerr << "\n[TEST] test_d3il_repeated_calls_accumulate_entries\n"
              << "  Source file: determine_3D_implicitlipid_reaction_probability.cpp\n"
              << "  Function:    determine_3D_implicitlipid_reaction_probability\n"
              << "  Scenario:    the same pair is evaluated twice in one time step.\n"
              << "  Criteria:    probvec and every reweighting list grow by one per\n"
              << "               call, and Dtot keeps accumulating rotational terms.\n";

    D3ILScenario s = d3il_make_scenario(/*ifaceZ*/ -48.0, /*rate*/ 10.0, /*nFreeLipids*/ 1000);

    d3il_run(s);
    const size_t sizeAfterFirst = s.moleculeList[0].probvec.size();
    const double dtotAfterFirst = s.biMolData.Dtot;
    std::cerr << "  after 1st call: probvec=" << sizeAfterFirst << " Dtot=" << dtotAfterFirst << '\n';

    d3il_run(s);
    const size_t sizeAfterSecond = s.moleculeList[0].probvec.size();
    std::cerr << "  after 2nd call: probvec=" << sizeAfterSecond << " Dtot=" << s.biMolData.Dtot << '\n';
    d3il_report(s);

    EXPECT_EQ(sizeAfterFirst, 1u) << "First call should append exactly one probability";
    EXPECT_EQ(sizeAfterSecond, 2u) << "Second call should append a second probability";
    EXPECT_EQ(s.moleculeList[0].currlist.size(), 2u) << "currlist should also have grown to two";
    EXPECT_EQ(s.moleculeList[0].currps_prev.size(), 2u) << "currps_prev should also have grown to two";
    EXPECT_GT(s.biMolData.Dtot, dtotAfterFirst)
        << "Dtot is modified in place, so a second call adds the rotational terms again";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - one per scenario so failures are reported individually
// while the whole set always runs.
// -----------------------------------------------------------------------------
TEST(Determine3DImplicitLipidProbability, DtotUpdate3DComplex) { test_d3il_dtot_update_3d_complex(); }
TEST(Determine3DImplicitLipidProbability, DtotUpdateMembraneComplex) { test_d3il_dtot_update_membrane_complex(); }
TEST(Determine3DImplicitLipidProbability, NoProbabilityOutsideRmax) { test_d3il_no_probability_when_outside_rmax(); }
TEST(Determine3DImplicitLipidProbability, ProbabilityRecordedInsideRmax) { test_d3il_probability_recorded_when_inside_rmax(); }
TEST(Determine3DImplicitLipidProbability, ZeroRateGivesZeroProbability) { test_d3il_zero_rate_gives_zero_probability(); }
TEST(Determine3DImplicitLipidProbability, DissociatedMoleculeSkipped) { test_d3il_dissociated_molecule_is_skipped(); }
TEST(Determine3DImplicitLipidProbability, ZeroFreeLipidsGivesZeroProbability) { test_d3il_zero_free_lipids_gives_zero_probability(); }
TEST(Determine3DImplicitLipidProbability, ProbabilityMonotonicInRateAndDensity) { test_d3il_probability_monotonic_in_rate_and_density(); }
TEST(Determine3DImplicitLipidProbability, NegativeSeparationForcesUnitProbability) { test_d3il_negative_separation_forces_unit_probability(); }
TEST(Determine3DImplicitLipidProbability, RepeatedCallsAccumulateEntries) { test_d3il_repeated_calls_accumulate_entries(); }