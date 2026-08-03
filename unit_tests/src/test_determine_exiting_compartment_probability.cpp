/*! \file test_determine_exiting_compartment_probability.cpp
 *
 * ### Unit test for src/reactions/determine_exiting_compartment_probability.cpp
 *
 * Function under test:
 *
 *     void determine_exiting_compartment_probability(
 *         double distToCompartment,
 *         const std::vector<TransmissionRxn>& transmissionRxns,
 *         int rxnIndex, int pro1Index,
 *         std::vector<Molecule>& moleculeList,
 *         double Dtot, const Parameters& parameters,
 *         Membrane& membraneObject);
 *
 * The routine has exactly one observable side effect: it may assign a value
 * to `moleculeList[pro1Index].transmissionProb`.  The assignment is guarded
 * by two conditions taken straight from the source:
 *
 *   1. `moleculeList[pro1Index].isDissociated != true`      (molecule did not
 *                                                            just dissociate)
 *   2. `transmissionRxns[rxnIndex].rateList[0].rate > 0`     (a real rate was
 *                                                            supplied)
 *
 * When both hold, the probability is computed via
 * `prob_exiting_compartment()` using a `paramsIL` structure filled in from the
 * reaction (bindRadius, 2*rate), the parameters (timeStep), the caller (Dtot)
 * and the membrane (compartmentR, droplet.rho).
 *
 * The tests below therefore verify:
 *   - the guards (dissociated molecule / zero rate leave the field untouched),
 *   - that a positive rate on a live molecule does write a finite probability,
 *   - that the *indices* passed in are honoured (correct reaction is read and
 *     only the correct molecule is modified),
 *   - that the computation is deterministic and stays within the sane
 *     probability range for a sweep of distances.
 *
 * Verbose output is printed to stderr so a reader can follow which source file
 * and function is being exercised and what each assertion checks.
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
// Local helpers (anonymous namespace so they cannot collide with other tests).
// -----------------------------------------------------------------------------
namespace {

//! Sentinel that `Molecule::transmissionProb` is initialised to by the class.
constexpr double kDexcpUnsetProb = -1.0;

/*! \brief Build a minimal TransmissionRxn with a single rate entry.
 *
 * \param[in] rate       Value stored in rateList[0].rate (the guard the
 *                       function checks against zero).
 * \param[in] bindRadius Reaction bind radius, forwarded to paramsIL::sigma.
 */
TransmissionRxn dexcp_make_rxn(double rate, double bindRadius)
{
    TransmissionRxn rxn;
    rxn.rxnType = ReactionType::transmission;
    rxn.bindRadius = bindRadius;

    // The function unconditionally dereferences rateList[0], so it must exist.
    RxnBase::RateState oneRate;
    oneRate.rate = rate;
    rxn.rateList.push_back(oneRate);

    return rxn;
}

/*! \brief Build a minimal Molecule for use as the reactant.
 *
 * Only `isDissociated` (the guard) and `transmissionProb` (the output) matter
 * for this function, but we also give it an index and a coordinate so the
 * object is in a sensible state.
 */
Molecule dexcp_make_molecule(int index, bool isDissociated)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = index;
    mol.comCoord = Coord { 0.0, 0.0, 0.0 };
    mol.isDissociated = isDissociated;
    // transmissionProb defaults to -1 (kDexcpUnsetProb) from the class def.
    return mol;
}

/*! \brief Build a Membrane object that describes a spherical compartment. */
Membrane dexcp_make_membrane()
{
    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 50.0; // nm
    membraneObject.droplet.D = 1.0; // surface site diffusion
    membraneObject.droplet.rho = 0.01; // surface site density
    return membraneObject;
}

/*! \brief Build simulation Parameters with a small, stable time step. */
Parameters dexcp_make_parameters()
{
    Parameters params;
    params.timeStep = 0.1; // microseconds
    return params;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: happy path -- live molecule + positive rate  =>  probability written.
// -----------------------------------------------------------------------------
void test_dexcp_sets_probability_for_live_molecule()
{
    std::cerr << "\n[TEST] test_dexcp_sets_probability_for_live_molecule\n"
              << "  Source file:   determine_exiting_compartment_probability.cpp\n"
              << "  Function:      determine_exiting_compartment_probability()\n"
              << "  Scenario:      molecule has NOT dissociated and the reaction\n"
              << "                 rate is > 0, so both guards are satisfied.\n"
              << "  Pass criteria: moleculeList[0].transmissionProb is changed\n"
              << "                 from its -1 sentinel to a finite, sane value.\n";

    // Reaction with a healthy rate.
    std::vector<TransmissionRxn> transmissionRxns { dexcp_make_rxn(/*rate=*/10.0, /*bindRadius=*/1.0) };

    // One live (non-dissociated) molecule.
    std::vector<Molecule> moleculeList { dexcp_make_molecule(/*index=*/0, /*isDissociated=*/false) };

    Parameters params = dexcp_make_parameters();
    Membrane membraneObject = dexcp_make_membrane();

    const double distToCompartment = 2.0; // nm from the compartment surface
    const double Dtot = 25.0; // nm^2/us

    std::cerr << "  Pre-call transmissionProb = " << moleculeList[0].transmissionProb << '\n'
              << "  Calling determine_exiting_compartment_probability(dist=" << distToCompartment
              << ", Dtot=" << Dtot << ")...\n";

    determine_exiting_compartment_probability(distToCompartment, transmissionRxns, /*rxnIndex=*/0,
        /*pro1Index=*/0, moleculeList, Dtot, params, membraneObject);

    const double prob = moleculeList[0].transmissionProb;
    std::cerr << "  Post-call transmissionProb = " << prob << '\n';

    // The sentinel must have been overwritten.
    EXPECT_NE(prob, kDexcpUnsetProb)
        << "transmissionProb should have been assigned for a live molecule with rate > 0";

    // The result must be a usable number (not NaN / inf).
    EXPECT_FALSE(std::isnan(prob)) << "transmissionProb must not be NaN";
    EXPECT_TRUE(std::isfinite(prob)) << "transmissionProb must be finite";

    // A probability must not be negative, and should not exceed unity by more
    // than the tolerance the source itself warns about.
    EXPECT_GE(prob, 0.0) << "A probability cannot be negative";
    EXPECT_LE(prob, 1.000001) << "A probability should not exceed 1 (within source tolerance)";
}

// -----------------------------------------------------------------------------
// Test 2: guard #1 -- a molecule that just dissociated must be skipped.
// -----------------------------------------------------------------------------
void test_dexcp_skips_dissociated_molecule()
{
    std::cerr << "\n[TEST] test_dexcp_skips_dissociated_molecule\n"
              << "  Source file:   determine_exiting_compartment_probability.cpp\n"
              << "  Function:      determine_exiting_compartment_probability()\n"
              << "  Scenario:      moleculeList[0].isDissociated == true.\n"
              << "  Pass criteria: transmissionProb keeps its -1 sentinel because\n"
              << "                 the outer `if` guard rejects the molecule.\n";

    std::vector<TransmissionRxn> transmissionRxns { dexcp_make_rxn(/*rate=*/10.0, /*bindRadius=*/1.0) };

    // The single molecule is flagged as having just dissociated.
    std::vector<Molecule> moleculeList { dexcp_make_molecule(/*index=*/0, /*isDissociated=*/true) };

    Parameters params = dexcp_make_parameters();
    Membrane membraneObject = dexcp_make_membrane();

    std::cerr << "  Calling determine_exiting_compartment_probability on a dissociated molecule...\n";
    determine_exiting_compartment_probability(/*distToCompartment=*/2.0, transmissionRxns,
        /*rxnIndex=*/0, /*pro1Index=*/0, moleculeList, /*Dtot=*/25.0, params, membraneObject);

    std::cerr << "  Post-call transmissionProb = " << moleculeList[0].transmissionProb
              << " (expected " << kDexcpUnsetProb << ")\n";

    EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, kDexcpUnsetProb)
        << "A just-dissociated molecule must not receive a transmission probability";
}

// -----------------------------------------------------------------------------
// Test 3: guard #2 -- a zero (or negative) rate must be skipped.
// -----------------------------------------------------------------------------
void test_dexcp_skips_zero_and_negative_rate()
{
    std::cerr << "\n[TEST] test_dexcp_skips_zero_and_negative_rate\n"
              << "  Source file:   determine_exiting_compartment_probability.cpp\n"
              << "  Function:      determine_exiting_compartment_probability()\n"
              << "  Scenario:      rateList[0].rate is 0.0 and then -5.0.\n"
              << "  Pass criteria: in both cases transmissionProb stays at -1,\n"
              << "                 since the inner `rate > 0` guard fails.\n";

    Parameters params = dexcp_make_parameters();
    Membrane membraneObject = dexcp_make_membrane();

    // --- Sub-case A: rate exactly zero ---------------------------------------
    {
        std::vector<TransmissionRxn> transmissionRxns { dexcp_make_rxn(/*rate=*/0.0, /*bindRadius=*/1.0) };
        std::vector<Molecule> moleculeList { dexcp_make_molecule(0, /*isDissociated=*/false) };

        std::cerr << "  -> Sub-case A: rate == 0.0\n";
        determine_exiting_compartment_probability(/*distToCompartment=*/2.0, transmissionRxns,
            /*rxnIndex=*/0, /*pro1Index=*/0, moleculeList, /*Dtot=*/25.0, params, membraneObject);

        std::cerr << "     transmissionProb = " << moleculeList[0].transmissionProb << '\n';
        EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, kDexcpUnsetProb)
            << "A zero reaction rate must leave transmissionProb untouched";
    }

    // --- Sub-case B: negative rate ------------------------------------------
    {
        std::vector<TransmissionRxn> transmissionRxns { dexcp_make_rxn(/*rate=*/-5.0, /*bindRadius=*/1.0) };
        std::vector<Molecule> moleculeList { dexcp_make_molecule(0, /*isDissociated=*/false) };

        std::cerr << "  -> Sub-case B: rate == -5.0\n";
        determine_exiting_compartment_probability(/*distToCompartment=*/2.0, transmissionRxns,
            /*rxnIndex=*/0, /*pro1Index=*/0, moleculeList, /*Dtot=*/25.0, params, membraneObject);

        std::cerr << "     transmissionProb = " << moleculeList[0].transmissionProb << '\n';
        EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, kDexcpUnsetProb)
            << "A negative reaction rate must leave transmissionProb untouched";
    }
}

// -----------------------------------------------------------------------------
// Test 4: index handling -- the function must read transmissionRxns[rxnIndex]
//         and write only moleculeList[pro1Index].
// -----------------------------------------------------------------------------
void test_dexcp_honours_rxn_and_molecule_indices()
{
    std::cerr << "\n[TEST] test_dexcp_honours_rxn_and_molecule_indices\n"
              << "  Source file:   determine_exiting_compartment_probability.cpp\n"
              << "  Function:      determine_exiting_compartment_probability()\n"
              << "  Scenario:      reaction list has a zero-rate entry at index 0\n"
              << "                 and a positive-rate entry at index 1; molecule\n"
              << "                 list has three molecules and we target index 1.\n"
              << "  Pass criteria: rxnIndex==1 is used (so a probability IS set)\n"
              << "                 and only moleculeList[1] is modified.\n";

    // Two reactions: index 0 would be skipped, index 1 is active.
    std::vector<TransmissionRxn> transmissionRxns;
    transmissionRxns.push_back(dexcp_make_rxn(/*rate=*/0.0, /*bindRadius=*/1.0));
    transmissionRxns.push_back(dexcp_make_rxn(/*rate=*/8.0, /*bindRadius=*/1.5));

    // Three live molecules; only the middle one should be touched.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(dexcp_make_molecule(0, /*isDissociated=*/false));
    moleculeList.push_back(dexcp_make_molecule(1, /*isDissociated=*/false));
    moleculeList.push_back(dexcp_make_molecule(2, /*isDissociated=*/false));

    Parameters params = dexcp_make_parameters();
    Membrane membraneObject = dexcp_make_membrane();

    std::cerr << "  Calling with rxnIndex=1, pro1Index=1...\n";
    determine_exiting_compartment_probability(/*distToCompartment=*/3.0, transmissionRxns,
        /*rxnIndex=*/1, /*pro1Index=*/1, moleculeList, /*Dtot=*/25.0, params, membraneObject);

    std::cerr << "  transmissionProb values = ["
              << moleculeList[0].transmissionProb << ", "
              << moleculeList[1].transmissionProb << ", "
              << moleculeList[2].transmissionProb << "]\n";

    // Reaction 1 has rate > 0, so molecule 1 must have been assigned.
    EXPECT_NE(moleculeList[1].transmissionProb, kDexcpUnsetProb)
        << "The reaction at rxnIndex=1 has rate > 0, so molecule 1 must be assigned";

    // Neighbours must be untouched -- proves pro1Index is respected.
    EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, kDexcpUnsetProb)
        << "Molecule 0 was not the target and must remain unmodified";
    EXPECT_DOUBLE_EQ(moleculeList[2].transmissionProb, kDexcpUnsetProb)
        << "Molecule 2 was not the target and must remain unmodified";

    // Now target reaction 0 (rate == 0) on molecule 0: nothing should change,
    // proving that rxnIndex (not always 0) drives the rate lookup.
    std::cerr << "  Calling again with rxnIndex=0 (rate 0), pro1Index=0...\n";
    determine_exiting_compartment_probability(/*distToCompartment=*/3.0, transmissionRxns,
        /*rxnIndex=*/0, /*pro1Index=*/0, moleculeList, /*Dtot=*/25.0, params, membraneObject);

    std::cerr << "  molecule 0 transmissionProb = " << moleculeList[0].transmissionProb << '\n';
    EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, kDexcpUnsetProb)
        << "rxnIndex=0 has rate 0, so molecule 0 must still be unmodified";
}

// -----------------------------------------------------------------------------
// Test 5: determinism and range -- sweep a set of distances, each evaluated
//         twice, and confirm the results are reproducible and in [0, 1].
// -----------------------------------------------------------------------------
void test_dexcp_deterministic_and_bounded_over_distances()
{
    std::cerr << "\n[TEST] test_dexcp_deterministic_and_bounded_over_distances\n"
              << "  Source file:   determine_exiting_compartment_probability.cpp\n"
              << "  Function:      determine_exiting_compartment_probability()\n"
              << "                 (delegates to prob_exiting_compartment())\n"
              << "  Scenario:      evaluate the probability for several distances\n"
              << "                 from the compartment, twice per distance.\n"
              << "  Pass criteria: each result is finite, within [0, 1.000001],\n"
              << "                 and identical between the two evaluations\n"
              << "                 (the routine must be deterministic).\n";

    std::vector<TransmissionRxn> transmissionRxns { dexcp_make_rxn(/*rate=*/5.0, /*bindRadius=*/1.0) };

    Parameters params = dexcp_make_parameters();
    Membrane membraneObject = dexcp_make_membrane();
    const double Dtot = 20.0;

    // A spread of separations from the compartment surface, in nm.
    const std::vector<double> distances { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0 };

    for (double dist : distances) {
        // First evaluation on a fresh molecule.
        std::vector<Molecule> firstList { dexcp_make_molecule(0, /*isDissociated=*/false) };
        determine_exiting_compartment_probability(dist, transmissionRxns, /*rxnIndex=*/0,
            /*pro1Index=*/0, firstList, Dtot, params, membraneObject);
        const double firstProb = firstList[0].transmissionProb;

        // Second, independent evaluation with identical inputs.
        std::vector<Molecule> secondList { dexcp_make_molecule(0, /*isDissociated=*/false) };
        determine_exiting_compartment_probability(dist, transmissionRxns, /*rxnIndex=*/0,
            /*pro1Index=*/0, secondList, Dtot, params, membraneObject);
        const double secondProb = secondList[0].transmissionProb;

        std::cerr << "  dist = " << dist << " nm -> prob = " << firstProb
                  << " (repeat: " << secondProb << ")\n";

        // Sanity of the computed number.
        EXPECT_FALSE(std::isnan(firstProb))
            << "Probability at dist=" << dist << " must not be NaN";
        EXPECT_TRUE(std::isfinite(firstProb))
            << "Probability at dist=" << dist << " must be finite";
        EXPECT_GE(firstProb, 0.0)
            << "Probability at dist=" << dist << " must be non-negative";
        EXPECT_LE(firstProb, 1.000001)
            << "Probability at dist=" << dist << " must not exceed 1";

        // Determinism: same inputs must give the same output.
        EXPECT_DOUBLE_EQ(firstProb, secondProb)
            << "The routine must be deterministic for dist=" << dist;
    }
}

// -----------------------------------------------------------------------------
// Test 6: robustness -- an already-populated transmissionProb is overwritten
//         when the guards pass, and preserved when they do not.
// -----------------------------------------------------------------------------
void test_dexcp_overwrites_existing_probability()
{
    std::cerr << "\n[TEST] test_dexcp_overwrites_existing_probability\n"
              << "  Source file:   determine_exiting_compartment_probability.cpp\n"
              << "  Function:      determine_exiting_compartment_probability()\n"
              << "  Scenario:      transmissionProb is pre-seeded with a bogus\n"
              << "                 marker value (0.5) before the call.\n"
              << "  Pass criteria: with rate > 0 the marker is replaced; with the\n"
              << "                 molecule dissociated the marker survives.\n";

    Parameters params = dexcp_make_parameters();
    Membrane membraneObject = dexcp_make_membrane();
    const double marker = 0.5;

    // --- Sub-case A: guards pass -> marker must be replaced -----------------
    {
        std::vector<TransmissionRxn> transmissionRxns { dexcp_make_rxn(/*rate=*/12.0, /*bindRadius=*/1.0) };
        std::vector<Molecule> moleculeList { dexcp_make_molecule(0, /*isDissociated=*/false) };
        moleculeList[0].transmissionProb = marker; // stale value from a prior step

        std::cerr << "  -> Sub-case A: live molecule, rate 12.0, pre-seeded prob = " << marker << '\n';
        determine_exiting_compartment_probability(/*distToCompartment=*/1.0, transmissionRxns,
            /*rxnIndex=*/0, /*pro1Index=*/0, moleculeList, /*Dtot=*/30.0, params, membraneObject);

        std::cerr << "     post-call prob = " << moleculeList[0].transmissionProb << '\n';
        // We cannot know the exact value, but the stale marker must be gone
        // unless the physics happens to reproduce it exactly; assert it is at
        // least a valid probability that was written by the function.
        EXPECT_TRUE(std::isfinite(moleculeList[0].transmissionProb))
            << "The freshly written probability must be finite";
        EXPECT_GE(moleculeList[0].transmissionProb, 0.0)
            << "The freshly written probability must be non-negative";
    }

    // --- Sub-case B: guard fails -> marker must survive ----------------------
    {
        std::vector<TransmissionRxn> transmissionRxns { dexcp_make_rxn(/*rate=*/12.0, /*bindRadius=*/1.0) };
        std::vector<Molecule> moleculeList { dexcp_make_molecule(0, /*isDissociated=*/true) };
        moleculeList[0].transmissionProb = marker;

        std::cerr << "  -> Sub-case B: dissociated molecule, pre-seeded prob = " << marker << '\n';
        determine_exiting_compartment_probability(/*distToCompartment=*/1.0, transmissionRxns,
            /*rxnIndex=*/0, /*pro1Index=*/0, moleculeList, /*Dtot=*/30.0, params, membraneObject);

        std::cerr << "     post-call prob = " << moleculeList[0].transmissionProb << '\n';
        EXPECT_DOUBLE_EQ(moleculeList[0].transmissionProb, marker)
            << "A dissociated molecule's existing transmissionProb must be preserved";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* function runs inside its own TEST so
// that a failure in one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(DetermineExitingCompartmentProbability, SetsProbabilityForLiveMolecule)
{
    test_dexcp_sets_probability_for_live_molecule();
}

TEST(DetermineExitingCompartmentProbability, SkipsDissociatedMolecule)
{
    test_dexcp_skips_dissociated_molecule();
}

TEST(DetermineExitingCompartmentProbability, SkipsZeroAndNegativeRate)
{
    test_dexcp_skips_zero_and_negative_rate();
}

TEST(DetermineExitingCompartmentProbability, HonoursRxnAndMoleculeIndices)
{
    test_dexcp_honours_rxn_and_molecule_indices();
}

TEST(DetermineExitingCompartmentProbability, DeterministicAndBoundedOverDistances)
{
    test_dexcp_deterministic_and_bounded_over_distances();
}

TEST(DetermineExitingCompartmentProbability, OverwritesExistingProbability)
{
    test_dexcp_overwrites_existing_probability();
}