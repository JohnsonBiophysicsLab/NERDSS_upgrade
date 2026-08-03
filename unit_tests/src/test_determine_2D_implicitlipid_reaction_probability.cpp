/*! \file test_determine_2D_implicitlipid_reaction_probability.cpp
 *
 * ### Unit test for src/reactions/determine_2D_implicitlipid_reaction_probability.cpp
 *
 * Function under test:
 *
 *     void determine_2D_implicitlipid_reaction_probability(
 *         int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
 *         std::vector<double>& ILTableIDs, BiMolData& biMolData,
 *         const Parameters& params, std::vector<Molecule>& moleculeList,
 *         std::vector<Complex>& complexList,
 *         const std::vector<ForwardRxn>& forwardRxns,
 *         const std::vector<BackRxn>& backRxns,
 *         std::vector<double>& IL2DbindingVec,
 *         std::vector<double>& IL2DUnbindingVec,
 *         Membrane& membraneObject, const int& relStateIndex)
 *
 * The routine is responsible for four distinct things, and each one is
 * covered by its own named test function below:
 *
 *   1. Augmenting `biMolData.Dtot` with the rotational contribution of both
 *      complexes and then *discretizing* the result (1-3 significant figures
 *      depending on magnitude) so that only a limited number of 2D lookup
 *      tables ever need to be built.
 *   2. Looking up an existing 2D binding-probability table entry keyed on
 *      (ka, Dtot, kb) inside ILTableIDs / IL2DbindingVec, or creating a new
 *      one (three doubles pushed onto ILTableIDs, one probability pushed onto
 *      IL2DbindingVec).
 *   3. Writing the resulting reaction probability (rho * tabulated value) into
 *      moleculeList[pro1].probvec and registering the encounter in the
 *      crossbase / mycrossint / crossrxn / curr* bookkeeping vectors.
 *   4. Short-circuiting entirely when the molecule just dissociated, or when
 *      the reaction rate is zero (only a placeholder 0 is pushed to probvec).
 *
 * Verbose progress information is written to stderr so that a reader of the
 * test log can follow exactly which behaviour is being exercised.
 */

#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// A small container that owns every object the function under test needs.
// -----------------------------------------------------------------------------
struct D2ilrpSystem {
    Parameters params;
    Membrane membrane;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<double> ilTableIDs;
    std::vector<double> il2DbindingVec;
    std::vector<double> il2DUnbindingVec;
    BiMolData biMolData;
};

// Fixed properties shared by all the scenarios below.
const double kD2ilrpTimeStep = 0.1; //!< us
const double kD2ilrpMagMol = 1.0; //!< COM-to-interface distance of both partners (nm)
const double kD2ilrpBindRadius = 1.0; //!< sigma (nm)
const double kD2ilrpLength3Dto2D = 10.0; //!< nm, converts the 3D rate into a 2D rate
const double kD2ilrpTotalSA = 10000.0; //!< nm^2 of membrane surface
const int kD2ilrpFreeLipids = 1000; //!< free implicit lipids in state 0
const int kD2ilrpProteins = 100; //!< binding-competent protein interfaces

/*! \brief Build a minimal but self-consistent system for the routine.
 *
 * \param[in] dtotInput   Translational Dtot handed to the function.
 * \param[in] drz         Rotational diffusion constant (z) of *both* complexes.
 * \param[in] rate        Forward (3D) reaction rate.
 * \param[in] kb          Back-reaction rate stored in backRxns.
 * \param[in] reversible  Whether the forward reaction advertises a back reaction.
 */
D2ilrpSystem d2ilrp_build_system(double dtotInput, double drz, double rate, double kb, bool reversible)
{
    D2ilrpSystem sys;

    // ---- simulation parameters -------------------------------------------
    sys.params.timeStep = kD2ilrpTimeStep;

    // ---- membrane / implicit lipid bookkeeping ---------------------------
    sys.membrane.implicitLipid = true;
    sys.membrane.totalSA = kD2ilrpTotalSA;
    sys.membrane.numberOfFreeLipidsEachState = std::vector<int> { kD2ilrpFreeLipids };
    sys.membrane.numberOfProteinEachState = std::vector<int> { kD2ilrpProteins };

    // ---- molecule 0: the real protein ------------------------------------
    Molecule protein;
    protein.index = 0;
    protein.myComIndex = 0;
    protein.molTypeIndex = 0;
    protein.comCoord = Coord { 0.0, 0.0, 0.0 };
    protein.interfaceList.push_back(Molecule::Iface(Coord { 0.0, 0.0, -kD2ilrpMagMol }));
    protein.isDissociated = false;

    // ---- molecule 1: the implicit lipid ----------------------------------
    Molecule lipid;
    lipid.index = 1;
    lipid.myComIndex = 1;
    lipid.molTypeIndex = 1;
    lipid.isImplicitLipid = true;
    lipid.isLipid = true;
    lipid.comCoord = Coord { 0.0, 0.0, 0.0 };
    lipid.interfaceList.push_back(Molecule::Iface(Coord { 0.0, 0.0, 0.0 }));

    sys.moleculeList.push_back(protein);
    sys.moleculeList.push_back(lipid);

    // ---- the two parent complexes ----------------------------------------
    Complex proteinCom;
    proteinCom.index = 0;
    proteinCom.memberList.push_back(0);
    proteinCom.Dr = Coord { drz, drz, drz };
    proteinCom.ncross = 0;

    Complex lipidCom;
    lipidCom.index = 1;
    lipidCom.memberList.push_back(1);
    lipidCom.Dr = Coord { drz, drz, drz };
    lipidCom.ncross = 0;

    sys.complexList.push_back(proteinCom);
    sys.complexList.push_back(lipidCom);

    // ---- the forward reaction --------------------------------------------
    ForwardRxn fwd;
    fwd.bindRadius = kD2ilrpBindRadius;
    fwd.length3Dto2D = kD2ilrpLength3Dto2D;
    fwd.isReversible = reversible;
    fwd.conjBackRxnIndex = 0;
    fwd.rateList.push_back(RxnBase::RateState());
    fwd.rateList[0].rate = rate;
    sys.forwardRxns.push_back(fwd);

    // ---- the (possibly unused) back reaction -----------------------------
    BackRxn back;
    back.rateList.push_back(RxnBase::RateState());
    back.rateList[0].rate = kb;
    sys.backRxns.push_back(back);

    // ---- the pairwise data block -----------------------------------------
    sys.biMolData = BiMolData(/*pro1*/ 0, /*pro2*/ 1, /*com1*/ 0, /*com2*/ 1,
        /*relIface1*/ 0, /*relIface2*/ 0, /*absIface1*/ 0, /*absIface2*/ 1,
        /*Dtot*/ dtotInput, /*magMol1*/ kD2ilrpMagMol, /*magMol2*/ kD2ilrpMagMol);

    return sys;
}

/*! \brief Reference implementation of the documented Dtot transformation.
 *
 * Mirrors the rule stated in the source: add the rotational contribution of
 * both complexes, then keep only one significant figure below 0.1 and two
 * significant figures above (with the corresponding scaling for the very
 * small decades), collapsing anything below 1e-50 to exactly zero.
 */
double d2ilrp_expected_dtot(double dtotInput, double drz, double dt, double magMol1, double magMol2)
{
    const double dr1 = 2.0 * magMol1 * (1.0 - std::cos(std::sqrt(2.0 * drz * dt)));
    const double dr2 = 2.0 * magMol2 * (1.0 - std::cos(std::sqrt(2.0 * drz * dt)));
    double dtot = dtotInput + (dr1 + dr2) / (4.0 * dt);

    double dtmp = 0.0;
    if (dtot < 0.0001)
        dtmp = dtot * 100000;
    else if (dtot < 0.001)
        dtmp = dtot * 10000;
    else if (dtot < 0.01)
        dtmp = dtot * 1000;
    else if (dtot < 0.1)
        dtmp = dtot * 100;
    else
        dtmp = dtot * 100;

    const int dOnes = int(std::round(dtmp));

    double out = 0.0;
    if (dtot < 0.0001)
        out = dOnes * 0.00001;
    else if (dtot < 0.001)
        out = dOnes * 0.0001;
    else if (dtot < 0.01)
        out = dOnes * 0.001;
    else if (dtot < 0.1)
        out = dOnes * 0.01;
    else
        out = dOnes * 0.01;

    if (out < 1E-50)
        out = 0.0;

    return out;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Dtot receives the rotational contribution and is then discretized.
// -----------------------------------------------------------------------------
void test_d2ilrp_dtot_rotation_and_discretization()
{
    std::cerr << "\n[TEST] test_d2ilrp_dtot_rotation_and_discretization\n"
              << "  Source file: determine_2D_implicitlipid_reaction_probability.cpp\n"
              << "  Checking:    biMolData.Dtot += rotational terms, then is snapped\n"
              << "               onto the coarse grid used to key the 2D tables.\n";

    // ---- case (a): no rotation, value that must be truncated to 2 sig figs.
    {
        const double dtotIn = 0.12345;
        D2ilrpSystem sys = d2ilrp_build_system(dtotIn, /*drz*/ 0.0, /*rate*/ 10.0, /*kb*/ 1.0, true);

        const double expected = d2ilrp_expected_dtot(dtotIn, 0.0, kD2ilrpTimeStep, kD2ilrpMagMol, kD2ilrpMagMol);
        std::cerr << "  (a) Dtot in = " << dtotIn << ", expected snapped Dtot = " << expected << "\n";

        // Pre-seed the lookup table with the *expected* key so that the heavy
        // 2D probability computation is skipped; the test only cares about Dtot.
        const double ktemp = sys.forwardRxns[0].rateList[0].rate / sys.forwardRxns[0].length3Dto2D;
        sys.ilTableIDs = { ktemp, expected, sys.backRxns[0].rateList[0].rate };
        sys.il2DbindingVec = { 0.25 };

        determine_2D_implicitlipid_reaction_probability(0, 0, 0, false, sys.ilTableIDs, sys.biMolData, sys.params,
            sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec,
            sys.il2DUnbindingVec, sys.membrane, 0);

        EXPECT_NEAR(sys.biMolData.Dtot, expected, 1e-12)
            << "Dtot should be snapped to 0.12 (2 significant figures) with no rotation";
        EXPECT_DOUBLE_EQ(sys.biMolData.Dtot, 0.12)
            << "0.12345 must discretize to exactly 0.12";
    }

    // ---- case (b): non-zero rotational diffusion increases Dtot -----------
    {
        const double dtotIn = 0.05;
        const double drz = 0.01;
        D2ilrpSystem sys = d2ilrp_build_system(dtotIn, drz, /*rate*/ 10.0, /*kb*/ 1.0, true);

        const double expected = d2ilrp_expected_dtot(dtotIn, drz, kD2ilrpTimeStep, kD2ilrpMagMol, kD2ilrpMagMol);
        std::cerr << "  (b) Dtot in = " << dtotIn << ", Dr.z = " << drz
                  << ", expected snapped Dtot = " << expected << "\n";

        const double ktemp = sys.forwardRxns[0].rateList[0].rate / sys.forwardRxns[0].length3Dto2D;
        sys.ilTableIDs = { ktemp, expected, sys.backRxns[0].rateList[0].rate };
        sys.il2DbindingVec = { 0.25 };

        determine_2D_implicitlipid_reaction_probability(0, 0, 0, false, sys.ilTableIDs, sys.biMolData, sys.params,
            sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec,
            sys.il2DUnbindingVec, sys.membrane, 0);

        EXPECT_NEAR(sys.biMolData.Dtot, expected, 1e-12)
            << "Dtot should include the rotational contribution of both complexes";
        EXPECT_GT(sys.biMolData.Dtot, dtotIn)
            << "Rotational diffusion must strictly increase Dtot";
    }

    // ---- case (c): an absurdly small Dtot must collapse to exactly zero ---
    {
        const double dtotIn = 1e-60;
        D2ilrpSystem sys = d2ilrp_build_system(dtotIn, /*drz*/ 0.0, /*rate*/ 10.0, /*kb*/ 1.0, true);

        std::cerr << "  (c) Dtot in = " << dtotIn << ", expected snapped Dtot = 0\n";

        const double ktemp = sys.forwardRxns[0].rateList[0].rate / sys.forwardRxns[0].length3Dto2D;
        sys.ilTableIDs = { ktemp, 0.0, sys.backRxns[0].rateList[0].rate };
        sys.il2DbindingVec = { 0.25 };

        determine_2D_implicitlipid_reaction_probability(0, 0, 0, false, sys.ilTableIDs, sys.biMolData, sys.params,
            sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec,
            sys.il2DUnbindingVec, sys.membrane, 0);

        EXPECT_DOUBLE_EQ(sys.biMolData.Dtot, 0.0)
            << "Anything below 1e-50 must be forced to exactly zero";
    }
}

// -----------------------------------------------------------------------------
// Test 2: an existing (ka, Dtot, kb) table entry is reused and the resulting
//         probability plus all bookkeeping vectors are filled in correctly.
// -----------------------------------------------------------------------------
void test_d2ilrp_reuses_cached_table_entry()
{
    std::cerr << "\n[TEST] test_d2ilrp_reuses_cached_table_entry\n"
              << "  Source file: determine_2D_implicitlipid_reaction_probability.cpp\n"
              << "  Checking:    a matching cache key is reused (tables do not grow),\n"
              << "               probvec.back() == rho * cachedValue, and the crossing\n"
              << "               bookkeeping vectors are all updated once.\n";

    const double dtotIn = 0.05; // 0.05 -> discretizes back onto itself
    const double rate = 10.0;
    const double kb = 1.0;
    D2ilrpSystem sys = d2ilrp_build_system(dtotIn, /*drz*/ 0.0, rate, kb, /*reversible*/ true);

    const double ktemp = rate / kD2ilrpLength3Dto2D; // intrinsic 2D rate
    const double cachedProb = 0.25; // pretend tabulated 2D probability
    sys.ilTableIDs = { ktemp, dtotIn, kb };
    sys.il2DbindingVec = { cachedProb };
    sys.il2DUnbindingVec = { 0.75 }; // must be left untouched by this routine

    // rho = (free lipids) / (total surface area)
    const double rho = static_cast<double>(kD2ilrpFreeLipids) / kD2ilrpTotalSA;
    const double expectedProb = rho * cachedProb;
    std::cerr << "  rho = " << rho << ", expected reaction probability = " << expectedProb << "\n";

    determine_2D_implicitlipid_reaction_probability(/*simItr*/ 7, /*rxnIndex*/ 0, /*rateIndex*/ 0,
        /*isStateChangeBackRxn*/ false, sys.ilTableIDs, sys.biMolData, sys.params, sys.moleculeList,
        sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec, sys.il2DUnbindingVec,
        sys.membrane, /*relStateIndex*/ 0);

    // --- the cache must NOT have grown ------------------------------------
    EXPECT_EQ(sys.ilTableIDs.size(), 3u) << "A matching key must be reused, not appended";
    EXPECT_EQ(sys.il2DbindingVec.size(), 1u) << "No new 2D table should be generated";
    EXPECT_EQ(sys.il2DUnbindingVec.size(), 1u) << "IL2DUnbindingVec is not used by this routine";

    // --- probability written to the *first* molecule only ------------------
    const Molecule& protein = sys.moleculeList[0];
    const Molecule& lipid = sys.moleculeList[1];
    ASSERT_EQ(protein.probvec.size(), 1u) << "Exactly one probability entry should be pushed";
    EXPECT_NEAR(protein.probvec.back(), expectedProb, 1e-12)
        << "probvec.back() should equal rho * tabulated probability (currnorm == 1)";
    EXPECT_TRUE(lipid.probvec.empty())
        << "The implicit lipid partner gets no probvec entry in this routine";

    // --- crossing bookkeeping --------------------------------------------
    ASSERT_EQ(protein.crossbase.size(), 1u) << "crossbase should record the lipid partner";
    EXPECT_EQ(protein.crossbase[0], 1) << "crossbase must hold pro2Index (the implicit lipid)";
    ASSERT_EQ(protein.mycrossint.size(), 1u);
    EXPECT_EQ(protein.mycrossint[0], 0) << "mycrossint must hold this molecule's relative iface";
    ASSERT_EQ(protein.crossrxn.size(), 1u);
    EXPECT_EQ(protein.crossrxn[0][0], 0) << "crossrxn[0] should be the reaction index";
    EXPECT_EQ(protein.crossrxn[0][1], 0) << "crossrxn[1] should be the rate index";
    EXPECT_EQ(protein.crossrxn[0][2], 0) << "crossrxn[2] should be isStateChangeBackRxn (false)";
    EXPECT_EQ(sys.complexList[0].ncross, 1) << "ncross of the parent complex must be incremented";

    // --- reweighting bookkeeping -----------------------------------------
    ASSERT_EQ(protein.currprevsep.size(), 1u);
    EXPECT_DOUBLE_EQ(protein.currprevsep[0], 0.0) << "Separation is stored as 0 for implicit lipids";
    ASSERT_EQ(protein.currlist.size(), 1u);
    EXPECT_EQ(protein.currlist[0], 1) << "currlist must hold the partner index";
    ASSERT_EQ(protein.currmyface.size(), 1u);
    EXPECT_EQ(protein.currmyface[0], 0);
    ASSERT_EQ(protein.currpface.size(), 1u);
    EXPECT_EQ(protein.currpface[0], 0) << "currpface must hold the partner's relative iface";
    ASSERT_EQ(protein.currprevnorm.size(), 1u);
    EXPECT_DOUBLE_EQ(protein.currprevnorm[0], 1.0) << "currnorm is hard-coded to 1.0";
    ASSERT_EQ(protein.currps_prev.size(), 1u);
    EXPECT_NEAR(protein.currps_prev[0], 1.0 - expectedProb, 1e-12)
        << "currps_prev should be 1 - probability";
}

// -----------------------------------------------------------------------------
// Test 3: with an empty cache the routine must build (and store) a brand new
//         2D table entry keyed on (ka, Dtot, kb).
// -----------------------------------------------------------------------------
void test_d2ilrp_creates_new_table_entry()
{
    std::cerr << "\n[TEST] test_d2ilrp_creates_new_table_entry\n"
              << "  Source file: determine_2D_implicitlipid_reaction_probability.cpp\n"
              << "  Checking:    an empty cache leads to three keys appended to\n"
              << "               ILTableIDs (ka, Dtot, kb) and one probability appended\n"
              << "               to IL2DbindingVec; the resulting probability is finite\n"
              << "               and non-negative.\n";

    const double dtotIn = 0.05;
    const double rate = 10.0;
    const double kb = 1.0;
    D2ilrpSystem sys = d2ilrp_build_system(dtotIn, /*drz*/ 0.0, rate, kb, /*reversible*/ true);

    // The cache starts empty, so the expensive 2D probability is computed.
    ASSERT_TRUE(sys.ilTableIDs.empty());
    ASSERT_TRUE(sys.il2DbindingVec.empty());

    determine_2D_implicitlipid_reaction_probability(0, 0, 0, false, sys.ilTableIDs, sys.biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec,
        sys.il2DUnbindingVec, sys.membrane, 0);

    const double expectedKa = rate / kD2ilrpLength3Dto2D;

    ASSERT_EQ(sys.ilTableIDs.size(), 3u) << "Exactly three key values must be appended";
    EXPECT_NEAR(sys.ilTableIDs[0], expectedKa, 1e-12) << "Key 0 is the intrinsic 2D rate ka";
    EXPECT_NEAR(sys.ilTableIDs[1], sys.biMolData.Dtot, 1e-12) << "Key 1 is the discretized Dtot";
    EXPECT_NEAR(sys.ilTableIDs[2], kb, 1e-12) << "Key 2 is the back-reaction rate kb";

    ASSERT_EQ(sys.il2DbindingVec.size(), 1u) << "One tabulated probability must be appended";
    std::cerr << "  Tabulated 2D binding value = " << sys.il2DbindingVec[0] << "\n";
    EXPECT_TRUE(std::isfinite(sys.il2DbindingVec[0]))
        << "The tabulated 2D binding probability must be a finite number";

    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u);
    const double prob = sys.moleculeList[0].probvec.back();
    std::cerr << "  Resulting reaction probability = " << prob << "\n";
    EXPECT_TRUE(std::isfinite(prob)) << "The reaction probability must be finite";
    EXPECT_GE(prob, 0.0) << "The reaction probability must be non-negative";

    // The encounter must still be registered.
    EXPECT_EQ(sys.moleculeList[0].crossbase.size(), 1u) << "The encounter must be registered";
    EXPECT_EQ(sys.complexList[0].ncross, 1) << "ncross must be incremented";
    EXPECT_TRUE(sys.il2DUnbindingVec.empty()) << "IL2DUnbindingVec must remain untouched";
}

// -----------------------------------------------------------------------------
// Test 4: an irreversible forward reaction must use kb = 0 for the table key,
//         even when a back reaction with a non-zero rate exists in backRxns.
// -----------------------------------------------------------------------------
void test_d2ilrp_irreversible_uses_zero_kb()
{
    std::cerr << "\n[TEST] test_d2ilrp_irreversible_uses_zero_kb\n"
              << "  Source file: determine_2D_implicitlipid_reaction_probability.cpp\n"
              << "  Checking:    when ForwardRxn::isReversible is false the kb key must\n"
              << "               be 0 regardless of backRxns content. Verified by seeding\n"
              << "               the cache with kb = 0 and confirming it is reused.\n";

    const double dtotIn = 0.05;
    const double rate = 10.0;
    const double storedKb = 5.0; // deliberately non-zero and *must be ignored*
    D2ilrpSystem sys = d2ilrp_build_system(dtotIn, /*drz*/ 0.0, rate, storedKb, /*reversible*/ false);

    const double ktemp = rate / kD2ilrpLength3Dto2D;
    sys.ilTableIDs = { ktemp, dtotIn, 0.0 }; // kb key == 0
    sys.il2DbindingVec = { 0.4 };

    determine_2D_implicitlipid_reaction_probability(0, 0, 0, false, sys.ilTableIDs, sys.biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec,
        sys.il2DUnbindingVec, sys.membrane, 0);

    // If kb had been taken from backRxns (5.0) the key would not match and the
    // routine would have appended a new entry.
    EXPECT_EQ(sys.ilTableIDs.size(), 3u)
        << "kb must be 0 for an irreversible reaction, so the cached entry matches";
    EXPECT_EQ(sys.il2DbindingVec.size(), 1u)
        << "No new table should be created for an irreversible reaction with a cached key";

    const double rho = static_cast<double>(kD2ilrpFreeLipids) / kD2ilrpTotalSA;
    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u);
    EXPECT_NEAR(sys.moleculeList[0].probvec.back(), rho * 0.4, 1e-12)
        << "The cached probability (0.4) scaled by rho should be used";
}

// -----------------------------------------------------------------------------
// Test 5: a molecule that just dissociated must be skipped entirely (only the
//         placeholder 0 probability is pushed).
// -----------------------------------------------------------------------------
void test_d2ilrp_dissociated_molecule_is_skipped()
{
    std::cerr << "\n[TEST] test_d2ilrp_dissociated_molecule_is_skipped\n"
              << "  Source file: determine_2D_implicitlipid_reaction_probability.cpp\n"
              << "  Checking:    Molecule::isDissociated == true leaves probvec.back() == 0,\n"
              << "               leaves the tables untouched and registers no crossing.\n";

    D2ilrpSystem sys = d2ilrp_build_system(/*dtot*/ 0.05, /*drz*/ 0.0, /*rate*/ 10.0, /*kb*/ 1.0, true);
    sys.moleculeList[0].isDissociated = true; // molecule just came apart this step

    determine_2D_implicitlipid_reaction_probability(0, 0, 0, false, sys.ilTableIDs, sys.biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec,
        sys.il2DUnbindingVec, sys.membrane, 0);

    // Dtot is still updated (it happens before the isDissociated guard).
    EXPECT_DOUBLE_EQ(sys.biMolData.Dtot, 0.05) << "Dtot is discretized before the guard";

    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u) << "Placeholder 0 is still pushed";
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].probvec.back(), 0.0)
        << "A dissociated molecule must keep a zero reaction probability";
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty()) << "No crossing should be registered";
    EXPECT_TRUE(sys.moleculeList[0].currlist.empty()) << "No reweighting entry should be created";
    EXPECT_EQ(sys.complexList[0].ncross, 0) << "ncross must not be incremented";
    EXPECT_TRUE(sys.ilTableIDs.empty()) << "No 2D table should be generated";
    EXPECT_TRUE(sys.il2DbindingVec.empty()) << "No 2D table should be generated";
}

// -----------------------------------------------------------------------------
// Test 6: a zero forward rate must also short-circuit the routine.
// -----------------------------------------------------------------------------
void test_d2ilrp_zero_rate_is_skipped()
{
    std::cerr << "\n[TEST] test_d2ilrp_zero_rate_is_skipped\n"
              << "  Source file: determine_2D_implicitlipid_reaction_probability.cpp\n"
              << "  Checking:    rate == 0 leaves probvec.back() == 0, builds no table\n"
              << "               and registers no crossing.\n";

    D2ilrpSystem sys = d2ilrp_build_system(/*dtot*/ 0.05, /*drz*/ 0.0, /*rate*/ 0.0, /*kb*/ 1.0, true);

    determine_2D_implicitlipid_reaction_probability(0, 0, 0, false, sys.ilTableIDs, sys.biMolData, sys.params,
        sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec,
        sys.il2DUnbindingVec, sys.membrane, 0);

    ASSERT_EQ(sys.moleculeList[0].probvec.size(), 1u) << "Placeholder 0 is still pushed";
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].probvec.back(), 0.0)
        << "A zero-rate reaction must produce a zero probability";
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty()) << "No crossing should be registered";
    EXPECT_EQ(sys.complexList[0].ncross, 0) << "ncross must not be incremented";
    EXPECT_TRUE(sys.ilTableIDs.empty()) << "No 2D table keys should be stored";
    EXPECT_TRUE(sys.il2DbindingVec.empty()) << "No 2D probability should be stored";
}

// -----------------------------------------------------------------------------
// Test 7: repeated calls append one probvec entry each time and keep reusing
//         the same single cache entry.
// -----------------------------------------------------------------------------
void test_d2ilrp_repeated_calls_accumulate()
{
    std::cerr << "\n[TEST] test_d2ilrp_repeated_calls_accumulate\n"
              << "  Source file: determine_2D_implicitlipid_reaction_probability.cpp\n"
              << "  Checking:    calling the routine twice on the same pair appends two\n"
              << "               probvec / crossbase entries, increments ncross twice, and\n"
              << "               still uses exactly one cached 2D table entry.\n";

    const double dtotIn = 0.05;
    const double rate = 10.0;
    const double kb = 1.0;
    D2ilrpSystem sys = d2ilrp_build_system(dtotIn, /*drz*/ 0.0, rate, kb, true);

    const double ktemp = rate / kD2ilrpLength3Dto2D;
    sys.ilTableIDs = { ktemp, dtotIn, kb };
    sys.il2DbindingVec = { 0.25 };

    for (int call = 0; call < 2; ++call) {
        std::cerr << "  -> call #" << (call + 1) << "\n";
        determine_2D_implicitlipid_reaction_probability(call, 0, 0, false, sys.ilTableIDs, sys.biMolData, sys.params,
            sys.moleculeList, sys.complexList, sys.forwardRxns, sys.backRxns, sys.il2DbindingVec,
            sys.il2DUnbindingVec, sys.membrane, 0);
    }

    const double rho = static_cast<double>(kD2ilrpFreeLipids) / kD2ilrpTotalSA;
    const double expectedProb = rho * 0.25;

    EXPECT_EQ(sys.moleculeList[0].probvec.size(), 2u) << "One probability per call";
    EXPECT_EQ(sys.moleculeList[0].crossbase.size(), 2u) << "One crossing record per call";
    EXPECT_EQ(sys.moleculeList[0].currlist.size(), 2u) << "One reweighting record per call";
    EXPECT_EQ(sys.complexList[0].ncross, 2) << "ncross must be incremented once per call";
    EXPECT_EQ(sys.ilTableIDs.size(), 3u) << "The single cached key must be reused";
    EXPECT_EQ(sys.il2DbindingVec.size(), 1u) << "The single cached probability must be reused";

    for (std::size_t i = 0; i < sys.moleculeList[0].probvec.size(); ++i) {
        EXPECT_NEAR(sys.moleculeList[0].probvec[i], expectedProb, 1e-12)
            << "Every stored probability should equal rho * cached value (entry " << i << ")";
    }

    // The second call re-discretizes the already discretized Dtot: the value
    // must be idempotent under the transformation.
    EXPECT_DOUBLE_EQ(sys.biMolData.Dtot, dtotIn)
        << "Discretizing an already discretized Dtot must be a no-op";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each scenario is its own TEST so that a failure in one
// scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(Determine2DImplicitLipidReactionProbability, DtotRotationAndDiscretization)
{
    test_d2ilrp_dtot_rotation_and_discretization();
}
TEST(Determine2DImplicitLipidReactionProbability, ReusesCachedTableEntry)
{
    test_d2ilrp_reuses_cached_table_entry();
}
TEST(Determine2DImplicitLipidReactionProbability, CreatesNewTableEntry)
{
    test_d2ilrp_creates_new_table_entry();
}
TEST(Determine2DImplicitLipidReactionProbability, IrreversibleUsesZeroKb)
{
    test_d2ilrp_irreversible_uses_zero_kb();
}
TEST(Determine2DImplicitLipidReactionProbability, DissociatedMoleculeIsSkipped)
{
    test_d2ilrp_dissociated_molecule_is_skipped();
}
TEST(Determine2DImplicitLipidReactionProbability, ZeroRateIsSkipped)
{
    test_d2ilrp_zero_rate_is_skipped();
}
TEST(Determine2DImplicitLipidReactionProbability, RepeatedCallsAccumulate)
{
    test_d2ilrp_repeated_calls_accumulate();
}