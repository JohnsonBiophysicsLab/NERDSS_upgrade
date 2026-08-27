/*! \file test_determine_entering_compartment_probability.cpp
 *
 * ### Unit test for src/reactions/determine_entering_compartment_probability.cpp
 *
 * Function under test:
 *
 *     void determine_entering_compartment_probability(
 *         double distToCompartment,
 *         const std::vector<TransmissionRxn>& transmissionRxns,
 *         int rxnIndex, int pro1Index,
 *         std::vector<Molecule>& moleculeList,
 *         double Dtot, const Parameters& parameters,
 *         Membrane& membraneObject);
 *
 * The routine is a thin wrapper whose *only* observable side effect is writing
 * a probability into `moleculeList[pro1Index].transmissionProb`.  Reading the
 * implementation, the contract is:
 *
 *   1. If the target molecule has `isDissociated == true`, the function returns
 *      immediately and `transmissionProb` is left untouched (default -1).
 *   2. If `transmissionRxns[rxnIndex].rateList[0].rate <= 0`, nothing happens
 *      either (the `if (rate > 0)` guard fails).
 *   3. Otherwise a `paramsIL` struct is filled in from the reaction, the
 *      Parameters and the Membrane, using an intrinsic rate of
 *      `ka = 2.0 * rate` (note the factor of two!), `R2D = 0.0`,
 *      `sigma = rxn.bindRadius`, `Dtot`, `dt = params.timeStep`,
 *      `compartmentR = membraneObject.compartmentR` and
 *      `compartSiteRho = membraneObject.droplet.rho`, and the result of
 *      `prob_entering_compartment(distToCompartment, params3D)` is stored on
 *      the molecule.
 *   4. Only the molecule at index `pro1Index` is modified.
 *   5. A probability larger than 1.000001 only prints a warning to std::cerr —
 *      the code does **not** exit — so the function is always safe to call.
 *
 * Every assertion below therefore checks one of these documented behaviours.
 * The "expected" probability is obtained by calling `prob_entering_compartment`
 * directly with an identically-constructed `paramsIL`, which lets us verify the
 * parameter marshalling (in particular the ka = 2*rate doubling) exactly,
 * without hard-coding a magic number from the closed-form expression.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers.  All names are prefixed with "decp_" (Determine Entering
// Compartment Probability) so they cannot collide with helpers from any other
// translation unit in the combined gtest binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal, fully-initialised TransmissionRxn.
 *
 * The function under test only reads `rateList[0].rate` and `bindRadius`, but
 * we still set the reaction type so the object is self-consistent.
 *
 * \param[in] rate       Intrinsic transmission rate stored in rateList[0].
 * \param[in] bindRadius Binding radius (sigma) used to fill paramsIL.
 */
TransmissionRxn decp_make_transmission_rxn(double rate, double bindRadius)
{
    TransmissionRxn rxn {};
    rxn.rxnType = ReactionType::transmission;
    rxn.bindRadius = bindRadius;
    // rateList must have at least one entry: the code indexes rateList[0].
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = rate;
    return rxn;
}

/*! \brief Build a minimal Molecule.
 *
 * Only `isDissociated` (read) and `transmissionProb` (written) matter here, but
 * we fill in the index/type so the object is not left half-initialised.
 */
Molecule decp_make_molecule(int index, bool isDissociated)
{
    Molecule mol {};
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = 0;
    mol.myComIndex = index;
    mol.isDissociated = isDissociated;
    mol.trajStatus = TrajStatus::none;
    // transmissionProb defaults to -1 (the "not computed" sentinel).
    return mol;
}

/*! \brief Build a Membrane carrying a compartment/droplet description. */
Membrane decp_make_membrane(double compartmentR, double dropletRho)
{
    Membrane membraneObject {};
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = compartmentR;
    membraneObject.droplet.rho = dropletRho;
    membraneObject.droplet.D = 1.0;
    // A water box big enough to contain the compartment (unused by the call,
    // but keeps the object physically sensible).
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 500.0, 500.0, 500.0 });
    return membraneObject;
}

/*! \brief Build Parameters with only the fields the routine consumes. */
Parameters decp_make_parameters(double timeStep)
{
    Parameters params {};
    params.timeStep = timeStep;
    params.nItr = 1;
    params.numMolTypes = 1;
    return params;
}

/*! \brief Reproduce, by hand, the paramsIL the function builds internally.
 *
 * This mirrors the source exactly, including the `ka = 2.0 * rate` doubling,
 * so we can call prob_entering_compartment() ourselves and compare.
 */
paramsIL decp_expected_params(const TransmissionRxn& rxn, double Dtot,
    const Parameters& params, const Membrane& membraneObject)
{
    paramsIL p {};
    p.R2D = 0.0;
    p.sigma = rxn.bindRadius;
    p.Dtot = Dtot;
    p.ka = 2.0 * rxn.rateList[0].rate; // NOTE: factor of two, as in the source
    p.dt = params.timeStep;
    p.compartmentR = membraneObject.compartmentR;
    p.compartSiteRho = membraneObject.droplet.rho;
    return p;
}

// Common, physically-sensible test inputs shared by several tests.
constexpr double kDecpBindRadius = 1.0;   // nm
constexpr double kDecpDtot = 1.0;         // nm^2 us^-1
constexpr double kDecpTimeStep = 0.1;     // us
constexpr double kDecpCompartmentR = 50.0; // nm
constexpr double kDecpRho = 0.1;          // sites nm^-2
constexpr double kDecpRate = 1.0;         // reaction rate (> 0 => active)
constexpr double kDecpDist = 2.0;         // nm from the compartment surface

} // namespace

// -----------------------------------------------------------------------------
// Test 1: an active reaction (rate > 0, molecule not dissociated) must write a
//         probability onto the target molecule.
// -----------------------------------------------------------------------------
void test_decp_sets_probability_for_active_reaction()
{
    std::cerr << "\n[TEST] test_decp_sets_probability_for_active_reaction\n"
              << "  Source file:   determine_entering_compartment_probability.cpp\n"
              << "  Function:      determine_entering_compartment_probability\n"
              << "  Scenario:      rate > 0 and the molecule has NOT just dissociated.\n"
              << "  Pass criteria: molecule.transmissionProb is overwritten (no longer\n"
              << "                 the -1 sentinel), is finite, and is a sane probability.\n";

    // --- arrange -------------------------------------------------------------
    std::vector<TransmissionRxn> transmissionRxns;
    transmissionRxns.push_back(decp_make_transmission_rxn(kDecpRate, kDecpBindRadius));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(decp_make_molecule(0, /*isDissociated=*/false));

    Parameters params = decp_make_parameters(kDecpTimeStep);
    Membrane membraneObject = decp_make_membrane(kDecpCompartmentR, kDecpRho);

    std::cerr << "  Initial transmissionProb (sentinel) = "
              << moleculeList[0].transmissionProb << '\n';
    EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, -1.0)
        << "Pre-condition: a fresh Molecule must start with transmissionProb == -1";

    // --- act -----------------------------------------------------------------
    std::cerr << "  Calling determine_entering_compartment_probability(dist=" << kDecpDist
              << ", rxnIndex=0, pro1Index=0, Dtot=" << kDecpDtot << ")\n";
    determine_entering_compartment_probability(kDecpDist, transmissionRxns, /*rxnIndex=*/0,
        /*pro1Index=*/0, moleculeList, kDecpDtot, params, membraneObject);

    // --- assert --------------------------------------------------------------
    const double prob = moleculeList[0].transmissionProb;
    std::cerr << "  Resulting transmissionProb = " << prob << '\n';

    EXPECT_NE(prob, -1.0)
        << "transmissionProb must have been overwritten for an active reaction";
    EXPECT_TRUE(std::isfinite(prob))
        << "transmissionProb must be a finite number, got " << prob;
    EXPECT_GE(prob, 0.0)
        << "A probability can never be negative";
    // The source itself treats anything above 1.000001 as an error condition,
    // so that is exactly the bound we check against here.
    EXPECT_LE(prob, 1.000001)
        << "Probability should not exceed unity for these (small dt) parameters";
}

// -----------------------------------------------------------------------------
// Test 2: the stored value must equal prob_entering_compartment() evaluated with
//         the paramsIL the function is documented to build - in particular the
//         intrinsic rate must be DOUBLED (ka = 2 * rate).
// -----------------------------------------------------------------------------
void test_decp_matches_direct_prob_entering_compartment()
{
    std::cerr << "\n[TEST] test_decp_matches_direct_prob_entering_compartment\n"
              << "  Source file:   determine_entering_compartment_probability.cpp\n"
              << "  Function:      determine_entering_compartment_probability\n"
              << "  Scenario:      compare the stored value against a hand-built\n"
              << "                 paramsIL + direct prob_entering_compartment() call.\n"
              << "  Pass criteria: exact bit-for-bit agreement with ka = 2*rate, and\n"
              << "                 disagreement with the (wrong) ka = rate variant.\n";

    std::vector<TransmissionRxn> transmissionRxns;
    transmissionRxns.push_back(decp_make_transmission_rxn(kDecpRate, kDecpBindRadius));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(decp_make_molecule(0, false));

    Parameters params = decp_make_parameters(kDecpTimeStep);
    Membrane membraneObject = decp_make_membrane(kDecpCompartmentR, kDecpRho);

    // Build the reference parameter set exactly as the source does.
    paramsIL expectedParams
        = decp_expected_params(transmissionRxns[0], kDecpDtot, params, membraneObject);
    const double expectedProb = prob_entering_compartment(kDecpDist, expectedParams);

    // And a deliberately-wrong variant with the un-doubled rate, used to prove
    // the ka doubling really is exercised (skipped if the two coincide).
    paramsIL singleRateParams = expectedParams;
    singleRateParams.ka = transmissionRxns[0].rateList[0].rate;
    const double singleRateProb = prob_entering_compartment(kDecpDist, singleRateParams);

    std::cerr << "  Reference prob (ka = 2*rate) = " << expectedProb << '\n';
    std::cerr << "  Variant   prob (ka =   rate) = " << singleRateProb << '\n';

    determine_entering_compartment_probability(kDecpDist, transmissionRxns, 0, 0,
        moleculeList, kDecpDtot, params, membraneObject);

    const double prob = moleculeList[0].transmissionProb;
    std::cerr << "  Stored    prob               = " << prob << '\n';

    EXPECT_DOUBLE_EQ(prob, expectedProb)
        << "Stored probability must equal prob_entering_compartment() evaluated with "
           "R2D=0, sigma=bindRadius, Dtot, ka=2*rate, dt=timeStep, "
           "compartmentR and compartSiteRho taken from the Membrane";

    if (expectedProb != singleRateProb) {
        EXPECT_NE(prob, singleRateProb)
            << "The intrinsic rate handed to prob_entering_compartment must be the "
               "DOUBLED rate (2*rate), not the raw reaction rate";
    } else {
        std::cerr << "  NOTE: doubling the rate did not change the probability for these\n"
                     "        inputs, so the ka-doubling discrimination check is skipped.\n";
    }
}

// -----------------------------------------------------------------------------
// Test 3: a molecule flagged as just-dissociated must be skipped entirely.
// -----------------------------------------------------------------------------
void test_decp_skips_dissociated_molecule()
{
    std::cerr << "\n[TEST] test_decp_skips_dissociated_molecule\n"
              << "  Source file:   determine_entering_compartment_probability.cpp\n"
              << "  Function:      determine_entering_compartment_probability\n"
              << "  Scenario:      molecule.isDissociated == true (it just came apart).\n"
              << "  Pass criteria: transmissionProb keeps whatever value it already had.\n";

    std::vector<TransmissionRxn> transmissionRxns;
    transmissionRxns.push_back(decp_make_transmission_rxn(kDecpRate, kDecpBindRadius));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(decp_make_molecule(0, /*isDissociated=*/true));

    // Stamp a recognisable marker so we can tell "untouched" from "recomputed".
    const double marker = 0.4242;
    moleculeList[0].transmissionProb = marker;

    Parameters params = decp_make_parameters(kDecpTimeStep);
    Membrane membraneObject = decp_make_membrane(kDecpCompartmentR, kDecpRho);

    std::cerr << "  transmissionProb before call = " << moleculeList[0].transmissionProb << '\n';
    determine_entering_compartment_probability(kDecpDist, transmissionRxns, 0, 0,
        moleculeList, kDecpDtot, params, membraneObject);
    std::cerr << "  transmissionProb after  call = " << moleculeList[0].transmissionProb << '\n';

    EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, marker)
        << "A dissociated molecule must be skipped, leaving transmissionProb untouched";
}

// -----------------------------------------------------------------------------
// Test 4: a reaction with a non-positive rate must be skipped (the `rate > 0`
//         guard).  Both a zero rate and a negative rate are exercised.
// -----------------------------------------------------------------------------
void test_decp_skips_non_positive_rate()
{
    std::cerr << "\n[TEST] test_decp_skips_non_positive_rate\n"
              << "  Source file:   determine_entering_compartment_probability.cpp\n"
              << "  Function:      determine_entering_compartment_probability\n"
              << "  Scenario:      rateList[0].rate is 0.0 and then -1.0.\n"
              << "  Pass criteria: the `rate > 0` guard rejects both, so\n"
              << "                 transmissionProb stays at its -1 sentinel.\n";

    Parameters params = decp_make_parameters(kDecpTimeStep);
    Membrane membraneObject = decp_make_membrane(kDecpCompartmentR, kDecpRho);

    // --- zero rate -----------------------------------------------------------
    {
        std::vector<TransmissionRxn> transmissionRxns;
        transmissionRxns.push_back(decp_make_transmission_rxn(0.0, kDecpBindRadius));

        std::vector<Molecule> moleculeList;
        moleculeList.push_back(decp_make_molecule(0, false));

        std::cerr << "  -> rate = 0.0\n";
        determine_entering_compartment_probability(kDecpDist, transmissionRxns, 0, 0,
            moleculeList, kDecpDtot, params, membraneObject);
        std::cerr << "     transmissionProb = " << moleculeList[0].transmissionProb << '\n';

        EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, -1.0)
            << "A zero-rate transmission reaction must not produce a probability";
    }

    // --- negative rate -------------------------------------------------------
    {
        std::vector<TransmissionRxn> transmissionRxns;
        transmissionRxns.push_back(decp_make_transmission_rxn(-1.0, kDecpBindRadius));

        std::vector<Molecule> moleculeList;
        moleculeList.push_back(decp_make_molecule(0, false));

        std::cerr << "  -> rate = -1.0\n";
        determine_entering_compartment_probability(kDecpDist, transmissionRxns, 0, 0,
            moleculeList, kDecpDtot, params, membraneObject);
        std::cerr << "     transmissionProb = " << moleculeList[0].transmissionProb << '\n';

        EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, -1.0)
            << "A negative-rate transmission reaction must not produce a probability";
    }
}

// -----------------------------------------------------------------------------
// Test 5: only the molecule at pro1Index may be modified; its neighbours in the
//         moleculeList must be left completely alone.
// -----------------------------------------------------------------------------
void test_decp_only_modifies_target_molecule()
{
    std::cerr << "\n[TEST] test_decp_only_modifies_target_molecule\n"
              << "  Source file:   determine_entering_compartment_probability.cpp\n"
              << "  Function:      determine_entering_compartment_probability\n"
              << "  Scenario:      three molecules in the list, pro1Index = 1.\n"
              << "  Pass criteria: molecule 1 gets a probability; molecules 0 and 2\n"
              << "                 keep their untouched -1 sentinel.\n";

    std::vector<TransmissionRxn> transmissionRxns;
    transmissionRxns.push_back(decp_make_transmission_rxn(kDecpRate, kDecpBindRadius));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(decp_make_molecule(0, false));
    moleculeList.push_back(decp_make_molecule(1, false)); // <-- the target
    moleculeList.push_back(decp_make_molecule(2, false));

    Parameters params = decp_make_parameters(kDecpTimeStep);
    Membrane membraneObject = decp_make_membrane(kDecpCompartmentR, kDecpRho);

    std::cerr << "  Calling with pro1Index = 1\n";
    determine_entering_compartment_probability(kDecpDist, transmissionRxns, 0, /*pro1Index=*/1,
        moleculeList, kDecpDtot, params, membraneObject);

    for (std::size_t i = 0; i < moleculeList.size(); ++i) {
        std::cerr << "  molecule " << i
                  << " transmissionProb = " << moleculeList[i].transmissionProb << '\n';
    }

    EXPECT_NE(moleculeList[1].transmissionProb, -1.0)
        << "The targeted molecule (index 1) should have received a probability";
    EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, -1.0)
        << "Molecule 0 is not the target and must be untouched";
    EXPECT_DOUBLE_EQ(moleculeList[2].transmissionProb, -1.0)
        << "Molecule 2 is not the target and must be untouched";
}

// -----------------------------------------------------------------------------
// Test 6: rxnIndex must select the correct reaction out of the list.
// -----------------------------------------------------------------------------
void test_decp_uses_requested_reaction_index()
{
    std::cerr << "\n[TEST] test_decp_uses_requested_reaction_index\n"
              << "  Source file:   determine_entering_compartment_probability.cpp\n"
              << "  Function:      determine_entering_compartment_probability\n"
              << "  Scenario:      two transmission reactions with very different rates\n"
              << "                 and binding radii; rxnIndex = 1 is requested.\n"
              << "  Pass criteria: the stored probability matches reaction 1's parameters\n"
              << "                 and not reaction 0's.\n";

    std::vector<TransmissionRxn> transmissionRxns;
    transmissionRxns.push_back(decp_make_transmission_rxn(0.01, 1.0));  // index 0
    transmissionRxns.push_back(decp_make_transmission_rxn(5.00, 3.0));  // index 1 (target)

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(decp_make_molecule(0, false));

    Parameters params = decp_make_parameters(kDecpTimeStep);
    Membrane membraneObject = decp_make_membrane(kDecpCompartmentR, kDecpRho);

    // Reference probabilities for each candidate reaction.
    paramsIL p0 = decp_expected_params(transmissionRxns[0], kDecpDtot, params, membraneObject);
    paramsIL p1 = decp_expected_params(transmissionRxns[1], kDecpDtot, params, membraneObject);
    const double prob0 = prob_entering_compartment(kDecpDist, p0);
    const double prob1 = prob_entering_compartment(kDecpDist, p1);

    std::cerr << "  Reference prob for reaction 0 = " << prob0 << '\n';
    std::cerr << "  Reference prob for reaction 1 = " << prob1 << '\n';

    determine_entering_compartment_probability(kDecpDist, transmissionRxns, /*rxnIndex=*/1, 0,
        moleculeList, kDecpDtot, params, membraneObject);

    const double prob = moleculeList[0].transmissionProb;
    std::cerr << "  Stored prob (rxnIndex = 1)    = " << prob << '\n';

    EXPECT_DOUBLE_EQ(prob, prob1)
        << "The reaction pointed at by rxnIndex must supply rate and bindRadius";

    if (prob0 != prob1) {
        EXPECT_NE(prob, prob0)
            << "The wrong reaction (index 0) must not have been used";
    } else {
        std::cerr << "  NOTE: both reactions gave identical probabilities, so the\n"
                     "        discrimination check is skipped.\n";
    }
}

// -----------------------------------------------------------------------------
// Test 7: the routine is a pure function of its inputs - calling it twice with
//         the same arguments must give the same answer, while changing the
//         distance must be reflected in the stored value.
// -----------------------------------------------------------------------------
void test_decp_is_deterministic_and_distance_dependent()
{
    std::cerr << "\n[TEST] test_decp_is_deterministic_and_distance_dependent\n"
              << "  Source file:   determine_entering_compartment_probability.cpp\n"
              << "  Function:      determine_entering_compartment_probability\n"
              << "  Scenario:      repeated calls with identical arguments, then a call\n"
              << "                 with a different distToCompartment.\n"
              << "  Pass criteria: identical inputs -> identical output (no RNG use), and\n"
              << "                 each distance reproduces its own direct reference value.\n";

    std::vector<TransmissionRxn> transmissionRxns;
    transmissionRxns.push_back(decp_make_transmission_rxn(kDecpRate, kDecpBindRadius));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(decp_make_molecule(0, false));

    Parameters params = decp_make_parameters(kDecpTimeStep);
    Membrane membraneObject = decp_make_membrane(kDecpCompartmentR, kDecpRho);

    // --- determinism ---------------------------------------------------------
    determine_entering_compartment_probability(kDecpDist, transmissionRxns, 0, 0,
        moleculeList, kDecpDtot, params, membraneObject);
    const double firstCall = moleculeList[0].transmissionProb;

    moleculeList[0].transmissionProb = -1.0; // reset the sentinel
    determine_entering_compartment_probability(kDecpDist, transmissionRxns, 0, 0,
        moleculeList, kDecpDtot, params, membraneObject);
    const double secondCall = moleculeList[0].transmissionProb;

    std::cerr << "  Call 1 (dist = " << kDecpDist << ") -> " << firstCall << '\n';
    std::cerr << "  Call 2 (dist = " << kDecpDist << ") -> " << secondCall << '\n';

    EXPECT_DOUBLE_EQ(firstCall, secondCall)
        << "The routine draws no random numbers, so repeated calls must agree exactly";

    // --- distance dependence -------------------------------------------------
    const double farDist = 20.0;
    paramsIL nearParams = decp_expected_params(transmissionRxns[0], kDecpDtot, params, membraneObject);
    const double nearRef = prob_entering_compartment(kDecpDist, nearParams);
    const double farRef = prob_entering_compartment(farDist, nearParams);

    moleculeList[0].transmissionProb = -1.0;
    determine_entering_compartment_probability(farDist, transmissionRxns, 0, 0,
        moleculeList, kDecpDtot, params, membraneObject);
    const double farProb = moleculeList[0].transmissionProb;

    std::cerr << "  Call 3 (dist = " << farDist << ") -> " << farProb
              << "  (direct reference " << farRef << ")\n";

    EXPECT_DOUBLE_EQ(firstCall, nearRef)
        << "distToCompartment = " << kDecpDist
        << " must be forwarded verbatim to prob_entering_compartment";
    EXPECT_DOUBLE_EQ(farProb, farRef)
        << "distToCompartment = " << farDist
        << " must be forwarded verbatim to prob_entering_compartment";
    EXPECT_TRUE(std::isfinite(farProb))
        << "The probability at a large separation must still be finite, got " << farProb;
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* helper is invoked from its own TEST
// so that a failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(DetermineEnteringCompartmentProbability, SetsProbabilityForActiveReaction)
{
    test_decp_sets_probability_for_active_reaction();
}

TEST(DetermineEnteringCompartmentProbability, MatchesDirectProbEnteringCompartment)
{
    test_decp_matches_direct_prob_entering_compartment();
}

TEST(DetermineEnteringCompartmentProbability, SkipsDissociatedMolecule)
{
    test_decp_skips_dissociated_molecule();
}

TEST(DetermineEnteringCompartmentProbability, SkipsNonPositiveRate)
{
    test_decp_skips_non_positive_rate();
}

TEST(DetermineEnteringCompartmentProbability, OnlyModifiesTargetMolecule)
{
    test_decp_only_modifies_target_molecule();
}

TEST(DetermineEnteringCompartmentProbability, UsesRequestedReactionIndex)
{
    test_decp_uses_requested_reaction_index();
}

TEST(DetermineEnteringCompartmentProbability, IsDeterministicAndDistanceDependent)
{
    test_decp_is_deterministic_and_distance_dependent();
}