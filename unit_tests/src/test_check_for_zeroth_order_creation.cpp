/*! \file test_check_for_zeroth_order_creation.cpp
 *
 * ### Unit test for src/reactions/check_for_zeroth_order_creation.cpp
 *
 * Function under test:
 * \code
 * void check_for_zeroth_order_creation(unsigned simItr, Parameters&, SimulVolume&,
 *                                      const std::vector<ForwardRxn>&,
 *                                      const std::vector<CreateDestructRxn>&,
 *                                      std::vector<Molecule>&, std::vector<Complex>&,
 *                                      std::vector<MolTemplate>&,
 *                                      std::map<std::string,int>&, copyCounters&,
 *                                      Membrane&);
 * \endcode
 *
 * The routine walks the list of creation/destruction reactions and, for every
 * reaction of type ReactionType::zerothOrderCreation, draws the number of
 * creation events from a Poisson distribution with
 *
 *     lambda = N_A * rate * Volume * nm3ToLiters * timeStep * usToSeconds
 *
 * where the Volume used depends on the boundary description held by the
 * Membrane object (sphere, compartment interior, compartment exterior, or the
 * plain water box).  Two very different code paths then follow:
 *
 *   * implicit-lipid products -> only counters are bumped
 *     (copyCounters::copyNumSpecies and Membrane::numberOfFreeLipidsEachState)
 *   * explicit products       -> real Molecules/Complexes are instantiated by
 *     create_molecule_and_complex_from_rxn()
 *
 * Testing strategy
 * ----------------
 * Because the number of events is random, the assertions below rely either on
 *   (a) deterministic corner cases (rate == 0, rate == -1 titration already
 *       satisfied, non-creation reaction types) where the event count is
 *       provably zero, or
 *   (b) invariants that must hold for *any* drawn event count (the lipid
 *       counter, the species counter and the observable counter must all move
 *       by exactly the same amount), or
 *   (c) lambda values chosen so extreme that the outcome is statistically
 *       certain (lambda ~ 60 -> P(0 events) ~ 1e-26;
 *                lambda ~ 2.5e-7 -> P(>0 events) ~ 2.5e-7).
 *
 * Only the implicit-lipid path is exercised with a non-zero event count so
 * that the test never has to build a fully initialised SimulVolume for the
 * real molecule-creation machinery.
 */

#include "math/constants.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (all prefixed with czoc_ so they cannot collide with other
// generated test translation units).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Make sure the global GSL random number generator exists.
 *
 * gtest_main.cpp defines `gsl_rng* r = nullptr;`.  check_for_zeroth_order_creation
 * calls rand_gsl(), which dereferences that generator, so we lazily initialise
 * it with a fixed seed if nobody else has done so yet.
 */
void czoc_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] global gsl_rng was null -> calling srand_gsl(1)\n";
         const gsl_rng_type *T;
         T = gsl_rng_default;
         r = gsl_rng_alloc(T);
         gsl_rng_set(r, 42);
    }
}

/*! \brief Parameters with a 1 microsecond time step (units used by the code). */
Parameters czoc_make_params()
{
    Parameters params;
    params.timeStep = 1.0; // microseconds
    return params;
}

/*! \brief A cubic 100x100x100 nm water box (volume = 1e6 nm^3), no compartment. */
Membrane czoc_make_box_membrane()
{
    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    membraneObject.isSphere = false;
    membraneObject.hasCompartment = false;
    membraneObject.compartmentR = 0.0;
    membraneObject.sphereVol = 0.0;
    // Two implicit lipid states, both starting empty.
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };
    return membraneObject;
}

/*! \brief Build a minimal MolTemplate with a single interface (abs index 0). */
MolTemplate czoc_make_moltemplate(const std::string& name, bool isImplicitLipid)
{
    MolTemplate oneTemp;
    oneTemp.molName = name;
    oneTemp.molTypeIndex = 0;
    oneTemp.isImplicitLipid = isImplicitLipid;
    oneTemp.isLipid = isImplicitLipid;
    oneTemp.copies = 0;
    oneTemp.mass = 1.0;
    oneTemp.radius = 1.0;
    oneTemp.D.x = 1.0;
    oneTemp.D.y = 1.0;
    oneTemp.D.z = 1.0;
    oneTemp.Dr.x = 0.01;
    oneTemp.Dr.y = 0.01;
    oneTemp.Dr.z = 0.01;

    // A single interface carrying a single (default) state with absolute index 0.
    Interface iface;
    iface.name = "a";
    iface.index = 0;
    iface.stateList.push_back(Interface::State('\0', 0));
    oneTemp.interfaceList.push_back(iface);

    return oneTemp;
}

/*! \brief Build a zeroth-order creation reaction producing molTypeIndex.
 *
 * \param[in] molTypeIndex   index into molTemplateList of the product
 * \param[in] absIfaceIndex  absolute interface/state index of the product iface
 * \param[in] rate           macroscopic creation rate in M/s
 */
CreateDestructRxn czoc_make_creation_rxn(int molTypeIndex, int absIfaceIndex, double rate)
{
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::zerothOrderCreation;
    rxn.creationRadius = 1.0;
    rxn.isObserved = false;
    rxn.observeLabel = "";

    RxnIface prodIface("a", molTypeIndex, absIfaceIndex, 0, '\0', false);
    CreateDestructRxn::CreateDestructMol prodMol(molTypeIndex, std::vector<RxnIface> { prodIface });
    prodMol.molName = "A";
    rxn.productMolList.push_back(prodMol);

    rxn.rateList.emplace_back(rate, std::vector<std::vector<RxnIface>> {});
    return rxn;
}

/*! \brief copyCounters with a two-element species counter array, all zero. */
copyCounters czoc_make_counters()
{
    copyCounters counterArrays;
    counterArrays.copyNumSpecies = std::vector<int> { 0, 0 };
    return counterArrays;
}

/*! \brief Reproduce the lambda the routine computes, for diagnostic printing. */
long double czoc_expected_lambda(double rate, double volume, double timeStep)
{
    return Constants::avogadro * rate * volume * Constants::nm3ToLiters * timeStep
        * Constants::usToSeconds;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: an empty reaction list must be a complete no-op.
// -----------------------------------------------------------------------------
void test_czoc_empty_reaction_list_is_noop()
{
    std::cerr << "\n[TEST] test_czoc_empty_reaction_list_is_noop\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation\n"
              << "  Scenario:      createDestructRxns is empty.\n"
              << "  Pass criteria: no molecules/complexes appear and no counter\n"
              << "                 or observable is modified.\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    SimulVolume simulVolume; // never used on this path
    copyCounters counterArrays = czoc_make_counters();

    std::vector<MolTemplate> molTemplateList { czoc_make_moltemplate("A", false) };
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {}; // <- empty
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};

    check_for_zeroth_order_creation(0, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    EXPECT_EQ(moleculeList.size(), 0u) << "No Molecule should have been created";
    EXPECT_EQ(complexList.size(), 0u) << "No Complex should have been created";
    EXPECT_EQ(counterArrays.copyNumSpecies[0], 0) << "Species counter must stay at zero";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 0)
        << "Free-lipid counter must stay at zero";
    EXPECT_TRUE(observablesList.empty()) << "Observables map must stay empty";
}

// -----------------------------------------------------------------------------
// Test 2: reactions that are not zerothOrderCreation must be skipped, even when
//         their rate is enormous.
// -----------------------------------------------------------------------------
void test_czoc_ignores_non_creation_reactions()
{
    std::cerr << "\n[TEST] test_czoc_ignores_non_creation_reactions\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation\n"
              << "  Scenario:      the only reaction present is a destruction\n"
              << "                 reaction with a very large rate.\n"
              << "  Pass criteria: the rxnType filter rejects it, so every counter\n"
              << "                 is untouched.\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    // Implicit lipid template so that, if the filter failed, the lipid counter
    // would visibly move (and no real molecule creation would be attempted).
    std::vector<MolTemplate> molTemplateList { czoc_make_moltemplate("IL", true) };

    // Same reaction data, but flagged as a destruction reaction.
    CreateDestructRxn destroyRxn = czoc_make_creation_rxn(0, 0, 1000.0);
    destroyRxn.rxnType = ReactionType::destruction;

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns { destroyRxn };
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};

    check_for_zeroth_order_creation(1, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    EXPECT_EQ(counterArrays.copyNumSpecies[0], 0)
        << "A destruction reaction must not increment the species counter";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 0)
        << "A destruction reaction must not increment the free-lipid counter";
    EXPECT_EQ(moleculeList.size(), 0u) << "No Molecule should have been created";
}

// -----------------------------------------------------------------------------
// Test 3: rate == 0 -> lambda == 0 -> prob == 1 -> zero events (deterministic).
//         Uses an explicit (non implicit-lipid) product, so this also proves
//         create_molecule_and_complex_from_rxn() is never reached.
// -----------------------------------------------------------------------------
void test_czoc_zero_rate_creates_nothing()
{
    std::cerr << "\n[TEST] test_czoc_zero_rate_creates_nothing\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation (explicit branch)\n"
              << "  Scenario:      zerothOrderCreation reaction with rate = 0 M/s.\n"
              << "  Pass criteria: lambda = 0 => P(0 events) = 1, so moleculeList,\n"
              << "                 complexList and the counters stay empty/zero.\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    std::vector<MolTemplate> molTemplateList { czoc_make_moltemplate("A", false) };
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns { czoc_make_creation_rxn(0, 0, 0.0) };
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};

    std::cerr << "  expected lambda = "
              << static_cast<double>(czoc_expected_lambda(0.0, membraneObject.waterBox.volume,
                     params.timeStep))
              << '\n';

    check_for_zeroth_order_creation(2, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    EXPECT_EQ(moleculeList.size(), 0u) << "Zero rate must create no Molecule";
    EXPECT_EQ(complexList.size(), 0u) << "Zero rate must create no Complex";
    EXPECT_EQ(counterArrays.copyNumSpecies[0], 0) << "Zero rate must not bump copyNumSpecies";
}

// -----------------------------------------------------------------------------
// Test 4: implicit lipid product with rate == 0 -> lipid counters untouched.
// -----------------------------------------------------------------------------
void test_czoc_implicit_lipid_zero_rate()
{
    std::cerr << "\n[TEST] test_czoc_implicit_lipid_zero_rate\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation (implicit-lipid branch)\n"
              << "  Scenario:      implicit lipid creation reaction with rate = 0 M/s.\n"
              << "  Pass criteria: both copyNumSpecies[state] and\n"
              << "                 numberOfFreeLipidsEachState[state] remain zero.\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    std::vector<MolTemplate> molTemplateList { czoc_make_moltemplate("IL", true) };
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns { czoc_make_creation_rxn(0, 0, 0.0) };
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};

    check_for_zeroth_order_creation(3, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    EXPECT_EQ(counterArrays.copyNumSpecies[0], 0)
        << "copyNumSpecies must be unchanged for a zero-rate lipid creation";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 0)
        << "numberOfFreeLipidsEachState must be unchanged for a zero-rate lipid creation";
    EXPECT_EQ(moleculeList.size(), 0u)
        << "Implicit lipid creation never instantiates explicit Molecules";
}

// -----------------------------------------------------------------------------
// Test 5: implicit lipid product with a large lambda.  Verifies both that
//         lipids really get created and that the two counters (plus the
//         observable) all advance by exactly the same amount.
// -----------------------------------------------------------------------------
void test_czoc_implicit_lipid_creates_lipids_and_keeps_counters_consistent()
{
    std::cerr << "\n[TEST] test_czoc_implicit_lipid_creates_lipids_and_keeps_counters_consistent\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation (implicit-lipid branch)\n"
              << "  Scenario:      rate = 100 M/s in a 1e6 nm^3 box, dt = 1 us,\n"
              << "                 reaction is observed with a registered label.\n"
              << "  Pass criteria: (a) at least one lipid is created (lambda ~ 60,\n"
              << "                     P(0 events) ~ 1e-26),\n"
              << "                 (b) copyNumSpecies delta == free-lipid delta,\n"
              << "                 (c) observable delta == copyNumSpecies delta,\n"
              << "                 (d) the unrelated state index 1 is untouched.\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    std::vector<MolTemplate> molTemplateList { czoc_make_moltemplate("IL", true) };
    std::vector<ForwardRxn> forwardRxns {};

    CreateDestructRxn rxn = czoc_make_creation_rxn(0, /*absIfaceIndex=*/0, /*rate=*/100.0);
    rxn.isObserved = true;
    rxn.observeLabel = "lipidObs";
    std::vector<CreateDestructRxn> createDestructRxns { rxn };

    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList { { "lipidObs", 0 } };

    const long double lambda
        = czoc_expected_lambda(100.0, membraneObject.waterBox.volume, params.timeStep);
    std::cerr << "  expected lambda = " << static_cast<double>(lambda) << '\n';

    check_for_zeroth_order_creation(4, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    const int speciesDelta = counterArrays.copyNumSpecies[0];
    const int lipidDelta = membraneObject.numberOfFreeLipidsEachState[0];
    const int observableDelta = observablesList["lipidObs"];

    std::cerr << "  drawn events: copyNumSpecies=" << speciesDelta
              << ", freeLipids=" << lipidDelta << ", observable=" << observableDelta << '\n';

    EXPECT_GT(speciesDelta, 0) << "With lambda ~ 60 at least one lipid must be created";
    EXPECT_EQ(speciesDelta, lipidDelta)
        << "copyNumSpecies and numberOfFreeLipidsEachState must advance together";
    EXPECT_EQ(observableDelta, speciesDelta)
        << "The observable must be incremented by exactly the number of events";
    EXPECT_EQ(counterArrays.copyNumSpecies[1], 0)
        << "A different species index must not be modified";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[1], 0)
        << "A different lipid state must not be modified";
    EXPECT_EQ(moleculeList.size(), 0u)
        << "Implicit lipid creation must not push entries into moleculeList";
}

// -----------------------------------------------------------------------------
// Test 6: an observed reaction whose label is missing from the observables map
//         must not crash and must not insert a new map entry.
// -----------------------------------------------------------------------------
void test_czoc_unregistered_observable_label_is_tolerated()
{
    std::cerr << "\n[TEST] test_czoc_unregistered_observable_label_is_tolerated\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation (implicit-lipid branch)\n"
              << "  Scenario:      isObserved == true but observeLabel is not in\n"
              << "                 the observables map.\n"
              << "  Pass criteria: lipids are still created and the observables map\n"
              << "                 remains empty (no accidental insertion).\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    std::vector<MolTemplate> molTemplateList { czoc_make_moltemplate("IL", true) };
    std::vector<ForwardRxn> forwardRxns {};

    CreateDestructRxn rxn = czoc_make_creation_rxn(0, 0, 100.0);
    rxn.isObserved = true;
    rxn.observeLabel = "neverDefined";
    std::vector<CreateDestructRxn> createDestructRxns { rxn };

    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {}; // deliberately empty

    check_for_zeroth_order_creation(5, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    EXPECT_TRUE(observablesList.empty())
        << "An undefined observable label must not be inserted into the map";
    EXPECT_GT(counterArrays.copyNumSpecies[0], 0)
        << "Creation must still happen even if the observable is undefined";
}

// -----------------------------------------------------------------------------
// Test 7: sphere boundary -> the sphere volume, not the water box, drives lambda.
// -----------------------------------------------------------------------------
void test_czoc_uses_sphere_volume_when_membrane_is_sphere()
{
    std::cerr << "\n[TEST] test_czoc_uses_sphere_volume_when_membrane_is_sphere\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation (volume selection)\n"
              << "  Scenario:      Membrane::isSphere = true, sphereVol = 1e6 nm^3,\n"
              << "                 waterBox volume deliberately set to 0.\n"
              << "  Pass criteria: lipids ARE created, proving sphereVol (and not the\n"
              << "                 zero-sized water box) was used for lambda.\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 62.0;
    membraneObject.sphereVol = 1.0e6; // nm^3
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 0.0, 0.0, 0.0 });
    membraneObject.hasCompartment = false;
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    std::vector<MolTemplate> molTemplateList { czoc_make_moltemplate("IL", true) };
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns { czoc_make_creation_rxn(0, 0, 100.0) };
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};

    std::cerr << "  expected lambda (sphere) = "
              << static_cast<double>(
                     czoc_expected_lambda(100.0, membraneObject.sphereVol, params.timeStep))
              << ", lambda if waterBox were used = "
              << static_cast<double>(
                     czoc_expected_lambda(100.0, membraneObject.waterBox.volume, params.timeStep))
              << '\n';

    check_for_zeroth_order_creation(6, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    std::cerr << "  drawn events = " << counterArrays.copyNumSpecies[0] << '\n';
    EXPECT_GT(counterArrays.copyNumSpecies[0], 0)
        << "Sphere volume should give lambda ~ 60, so events are essentially certain";
    EXPECT_EQ(counterArrays.copyNumSpecies[0], membraneObject.numberOfFreeLipidsEachState[0])
        << "Species and lipid counters must stay in lockstep";
}

// -----------------------------------------------------------------------------
// Test 8: compartment interior -> tiny compartment volume suppresses creation.
// -----------------------------------------------------------------------------
void test_czoc_uses_compartment_interior_volume()
{
    std::cerr << "\n[TEST] test_czoc_uses_compartment_interior_volume\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation (volume selection)\n"
              << "  Scenario:      hasCompartment = true, template insideCompartment,\n"
              << "                 compartmentR = 0.1 nm (V ~ 4.2e-3 nm^3) inside a\n"
              << "                 1e6 nm^3 box, rate = 100 M/s.\n"
              << "  Pass criteria: lambda ~ 2.5e-7 so NO event is drawn; had the box\n"
              << "                 volume been used lambda would have been ~60.\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 0.1; // nm -> essentially no volume

    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    MolTemplate lipidTemp = czoc_make_moltemplate("IL", true);
    lipidTemp.insideCompartment = true;
    std::vector<MolTemplate> molTemplateList { lipidTemp };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns { czoc_make_creation_rxn(0, 0, 100.0) };
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};

    const double vCompartment
        = 4.0 / 3.0 * M_PI * std::pow(membraneObject.compartmentR, 3.0);
    std::cerr << "  compartment volume = " << vCompartment << " nm^3, expected lambda = "
              << static_cast<double>(czoc_expected_lambda(100.0, vCompartment, params.timeStep))
              << '\n';

    check_for_zeroth_order_creation(7, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    std::cerr << "  drawn events = " << counterArrays.copyNumSpecies[0] << '\n';
    EXPECT_EQ(counterArrays.copyNumSpecies[0], 0)
        << "A ~4e-3 nm^3 compartment must give lambda ~ 2.5e-7 -> no creation events";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 0)
        << "No lipid should be added when the compartment interior volume is used";
}

// -----------------------------------------------------------------------------
// Test 9: compartment exterior -> (box - compartment) volume drives lambda.
// -----------------------------------------------------------------------------
void test_czoc_uses_compartment_exterior_volume()
{
    std::cerr << "\n[TEST] test_czoc_uses_compartment_exterior_volume\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation (volume selection)\n"
              << "  Scenario:      hasCompartment = true, template outsideCompartment,\n"
              << "                 compartmentR = 0.1 nm, box volume = 1e6 nm^3.\n"
              << "  Pass criteria: exterior volume ~ 1e6 nm^3 -> lambda ~ 60, so at\n"
              << "                 least one lipid is created.\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 0.1;

    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    MolTemplate lipidTemp = czoc_make_moltemplate("IL", true);
    lipidTemp.outsideCompartment = true;
    std::vector<MolTemplate> molTemplateList { lipidTemp };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns { czoc_make_creation_rxn(0, 0, 100.0) };
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};

    const double vCompartment
        = 4.0 / 3.0 * M_PI * std::pow(membraneObject.compartmentR, 3.0);
    const double vExterior = membraneObject.waterBox.volume - vCompartment;
    std::cerr << "  exterior volume = " << vExterior << " nm^3, expected lambda = "
              << static_cast<double>(czoc_expected_lambda(100.0, vExterior, params.timeStep))
              << '\n';

    check_for_zeroth_order_creation(8, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    std::cerr << "  drawn events = " << counterArrays.copyNumSpecies[0] << '\n';
    EXPECT_GT(counterArrays.copyNumSpecies[0], 0)
        << "The compartment exterior volume should yield lambda ~ 60 -> events certain";
    EXPECT_EQ(counterArrays.copyNumSpecies[0], membraneObject.numberOfFreeLipidsEachState[0])
        << "Species and lipid counters must stay in lockstep";
}

// -----------------------------------------------------------------------------
// Test 10: the special rate == -1 "strong titration" branch.  When the target
//          copy number is already met the computed event count is clamped to 0.
// -----------------------------------------------------------------------------
void test_czoc_titration_rate_with_target_already_met()
{
    std::cerr << "\n[TEST] test_czoc_titration_rate_with_target_already_met\n"
              << "  Source file:   check_for_zeroth_order_creation.cpp\n"
              << "  Function:      check_for_zeroth_order_creation (explicit branch,\n"
              << "                 rate == -1 titration special case)\n"
              << "  Scenario:      rate = -1 with copies == monomerList.size() (== 3).\n"
              << "  Pass criteria: numEvents = copies - monomers = 0, so nothing is\n"
              << "                 created and no counter moves.  (A negative lambda\n"
              << "                 also makes prob > 1 so the Poisson loop is skipped.)\n";

    czoc_ensure_rng();

    Parameters params = czoc_make_params();
    Membrane membraneObject = czoc_make_box_membrane();
    SimulVolume simulVolume;
    copyCounters counterArrays = czoc_make_counters();

    // Explicit (non-lipid) template whose monomer population already equals the
    // requested copy number, so the titration branch asks for zero new molecules.
    MolTemplate oneTemp = czoc_make_moltemplate("A", false);
    oneTemp.copies = 3;
    oneTemp.monomerList = std::vector<int> { 0, 1, 2 };
    std::vector<MolTemplate> molTemplateList { oneTemp };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<CreateDestructRxn> createDestructRxns { czoc_make_creation_rxn(0, 0, -1.0) };
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::map<std::string, int> observablesList {};

    check_for_zeroth_order_creation(9, params, simulVolume, forwardRxns, createDestructRxns,
        moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);

    EXPECT_EQ(moleculeList.size(), 0u)
        << "Titration must create nothing when copies == monomerList.size()";
    EXPECT_EQ(complexList.size(), 0u)
        << "Titration must create no Complex when the target is already met";
    EXPECT_EQ(counterArrays.copyNumSpecies[0], 0)
        << "Species counters must not move when zero molecules are titrated in";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario is its own TEST so that a failure in one
// scenario does not stop the remaining ones from running.
// -----------------------------------------------------------------------------
TEST(CheckForZerothOrderCreation, EmptyReactionListIsNoop)
{
    test_czoc_empty_reaction_list_is_noop();
}
TEST(CheckForZerothOrderCreation, IgnoresNonCreationReactions)
{
    test_czoc_ignores_non_creation_reactions();
}
TEST(CheckForZerothOrderCreation, ZeroRateCreatesNothing)
{
    test_czoc_zero_rate_creates_nothing();
}
TEST(CheckForZerothOrderCreation, ImplicitLipidZeroRate)
{
    test_czoc_implicit_lipid_zero_rate();
}
TEST(CheckForZerothOrderCreation, ImplicitLipidCountersConsistent)
{
    test_czoc_implicit_lipid_creates_lipids_and_keeps_counters_consistent();
}
TEST(CheckForZerothOrderCreation, UnregisteredObservableLabelTolerated)
{
    test_czoc_unregistered_observable_label_is_tolerated();
}
TEST(CheckForZerothOrderCreation, UsesSphereVolume)
{
    test_czoc_uses_sphere_volume_when_membrane_is_sphere();
}
TEST(CheckForZerothOrderCreation, UsesCompartmentInteriorVolume)
{
    test_czoc_uses_compartment_interior_volume();
}
TEST(CheckForZerothOrderCreation, UsesCompartmentExteriorVolume)
{
    test_czoc_uses_compartment_exterior_volume();
}
TEST(CheckForZerothOrderCreation, TitrationTargetAlreadyMet)
{
    test_czoc_titration_rate_with_target_already_met();
}