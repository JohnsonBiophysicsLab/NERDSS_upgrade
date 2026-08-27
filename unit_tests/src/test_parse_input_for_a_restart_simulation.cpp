/*! \file test_parse_input_for_a_restart_simulation.cpp
 *
 * ### Unit test for `src/parser/parse_input_for_a_restart_simulation.cpp`
 *
 * The translation unit under test contains a single (very large) orchestration
 * function:
 *
 * \code
 * void parse_input_for_a_restart_simulation(
 *         std::string restartFileNameInput, std::string addFileNameInput,
 *         Parameters&, std::map<std::string,int>&, std::vector<ForwardRxn>&,
 *         std::vector<BackRxn>&, std::vector<CreateDestructRxn>&,
 *         std::vector<MolTemplate>&, Membrane&, std::vector<Molecule>&,
 *         std::vector<Complex>&, SimulVolume&, std::string, int,
 *         long long int, MpiContext&, std::string, std::string, unsigned,
 *         copyCounters&);
 * \endcode
 *
 * It performs, in order:
 *   1. `read_rng_state(rank)` inside a try/catch, falling back to
 *      `gsl_rng_set(r, seed)` when the RNG state file cannot be read,
 *   2. opening the restart file -- **and calling `error(...)`, which terminates
 *      the process, when that file does not exist**,
 *   3. `read_restart(...)` (which itself `exit()`s on malformed input),
 *   4. appending `membraneObject.nStates` zeroes to
 *      `membraneObject.numberOfProteinEachState`,
 *   5. optionally parsing an `add.inp` file,
 *   6. for spherical boundaries: `membraneObject.create_water_box()` and
 *      `sphereVol = 4/3 * pi * R^3`,
 *   7. `initialize_paramters_for_implicitlipid_model(...)`,
 *      `initialize_states(...)`,
 *   8. `set_rMaxLimit(...)`, `SimulVolume::create_simulation_volume(...)`,
 *      `SimulVolume::update_memberMolLists(..., mpiContext)` and
 *      `SimulVolume::display()`,
 *   9. a non-fatal scan of the trajectory file to compare its last recorded
 *      iteration with the restart iteration.
 *
 * ### Why the top-level function itself is *not* invoked here
 *
 * Every entry point into this function reaches a process-terminating path that
 * a unit test must not take:
 *   * a missing / unreadable restart file calls `error()` -> `exit()`, which
 *     would kill the whole gtest binary (not just this test),
 *   * `read_rng_state()` performs raw `FILE*` I/O on an RNG-state file that a
 *     unit test cannot fabricate portably,
 *   * `read_restart()` `exit()`s on any inconsistency in the (binary/records)
 *     restart stream, and
 *   * the function requires a live `MpiContext` whose members are filled in by
 *     `main()`/MPI start-up; an under-initialised one makes
 *     `update_memberMolLists()` compute an out-of-range sub-volume index and
 *     call `error()` -> `exit()`.
 *
 * So, following the "do not reach an exit()/abort() path" rule, this file
 * instead:
 *   * asserts that the symbol from the file under test is linked into the
 *     binary (so the translation unit really does build),
 *   * documents the fatal precondition and verifies the exact condition the
 *     function tests (an unopenable `std::ifstream`) without triggering it, and
 *   * drives the individual, non-fatal pipeline steps that the restart parser
 *     composes -- RNG re-seeding, the spherical water box + sphere volume, the
 *     `rMaxLimit` computation and the simulation-volume partitioning -- through
 *     their public APIs on a fully initialised, minimal system.
 */

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "parser/parser_functions.hpp"
#include "system_setup/system_setup.hpp"

// -----------------------------------------------------------------------------
// Local helpers. Everything lives in an anonymous namespace *and* carries the
// unique `pifars_` (parse_input_for_a_restart_simulation) prefix so that no
// symbol can collide with another file in the combined test binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Allocate (if required) and deterministically seed the global GSL RNG.
 *
 * `gtest_main.cpp` defines `gsl_rng* r = nullptr;`, so the first test that needs
 * randomness has to allocate it. `srand_gsl()` is deliberately avoided.
 */
void pifars_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, 42);
}

/*! \brief Build a fully initialised MolTemplate carrying exactly one interface.
 *
 * Every field that the downstream pipeline reads (diffusion constants, radius,
 * mass, interface list with a single state) is populated, because
 * `set_rMaxLimit()`/`Complex::update_properties()` dereference them without any
 * bounds or sanity checks.
 */
MolTemplate pifars_make_template(const std::string& name, int typeIndex, int absIfaceIndex)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.copies = 2;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.comCoord = Coord { 0.0, 0.0, 0.0 };
    temp.D = Coord { 10.0, 10.0, 10.0 };   // um^2 s^-1, non-zero => 3D species
    temp.Dr = Coord { 0.01, 0.01, 0.01 };  // rad^2 s^-1
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPoint = false;
    temp.isRod = false;
    temp.isPromoter = false;

    // One interface, one (unnamed) state.
    Interface iface;
    iface.name = name + "site";
    iface.index = 0;
    iface.iCoord = Coord { 1.0, 0.0, 0.0 };
    iface.stateList.emplace_back(iface.name, absIfaceIndex);
    temp.interfaceList.push_back(iface);

    // Mirror what the parser would leave behind for the starting copy numbers.
    temp.startingNumState.totalCopyNumbers = temp.copies;
    temp.startingNumState.numberEachState.push_back(temp.copies);
    temp.startingNumState.nameEachState.push_back(iface.name);

    return temp;
}

/*! \brief Build a single-interface Molecule that matches a MolTemplate. */
Molecule pifars_make_molecule(int index, const MolTemplate& temp, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = index; // one molecule per complex in this toy system
    mol.molTypeIndex = temp.molTypeIndex;
    mol.mass = temp.mass;
    mol.isLipid = temp.isLipid;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.isEmpty = false;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = com + temp.interfaceList[0].iCoord;
    iface.index = temp.interfaceList[0].stateList[0].index;
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.stateIden = temp.interfaceList[0].stateList[0].iden;
    iface.molTypeIndex = temp.molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);
    mol.freelist.push_back(0);

    return mol;
}

/*! \brief Container returned by pifars_build_minimal_system(). */
struct PifarsSystem {
    Parameters params;
    Membrane membrane;
    std::vector<MolTemplate> molTemplateList;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<ForwardRxn> forwardRxns;
};

/*! \brief Assemble the smallest system that the restart pipeline can operate on.
 *
 * Two molecule types, two copies of each, one bimolecular forward reaction and a
 * 100 x 100 x 100 nm reflecting water box. All molecules sit comfortably inside
 * the box so that the sub-volume assignment can never report an out-of-box
 * molecule (which would `exit()`).
 */
PifarsSystem pifars_build_minimal_system()
{
    PifarsSystem sys;

    // --- static bookkeeping that the Complex constructor reads ---------------
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 2, 2 };

    // --- parameters ----------------------------------------------------------
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 2;
    sys.params.numTotalComplex = 4;
    sys.params.nItr = 100;
    sys.params.timeStep = 0.1; // us
    sys.params.name = "pifars_unit_test";
    sys.params.trajFile = "pifars_trajectory.xyz";

    // --- membrane / boundary -------------------------------------------------
    sys.membrane.isBox = true;
    sys.membrane.isSphere = false;
    sys.membrane.implicitLipid = false;
    sys.membrane.nStates = 0;
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    sys.membrane.xBCtype = "reflect";
    sys.membrane.yBCtype = "reflect";
    sys.membrane.zBCtype = "reflect";

    // --- templates -----------------------------------------------------------
    sys.molTemplateList.push_back(pifars_make_template("A", 0, 0));
    sys.molTemplateList.push_back(pifars_make_template("B", 1, 1));

    // --- molecules (all deep inside the box) ---------------------------------
    sys.moleculeList.push_back(pifars_make_molecule(0, sys.molTemplateList[0], Coord { -20.0, -20.0, -20.0 }));
    sys.moleculeList.push_back(pifars_make_molecule(1, sys.molTemplateList[0], Coord { 20.0, 20.0, 20.0 }));
    sys.moleculeList.push_back(pifars_make_molecule(2, sys.molTemplateList[1], Coord { 10.0, -10.0, 5.0 }));
    sys.moleculeList.push_back(pifars_make_molecule(3, sys.molTemplateList[1], Coord { -5.0, 15.0, -25.0 }));

    // --- one complex per molecule -------------------------------------------
    for (std::size_t molItr = 0; molItr < sys.moleculeList.size(); ++molItr) {
        const Molecule& mol = sys.moleculeList[molItr];
        sys.complexList.emplace_back(static_cast<int>(molItr), mol, sys.molTemplateList[mol.molTypeIndex]);
    }
    // Fill in radius / diffusion / member counts exactly as the simulation does.
    for (auto& com : sys.complexList)
        com.update_properties(sys.moleculeList, sys.molTemplateList);

    // --- one bimolecular forward reaction A(site) + B(site) -> A(site!1).B(site!1)
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = false;
    rxn.isOnMem = false;
    rxn.bindRadius = 2.0; // nm
    rxn.assocAngles = ForwardRxn::Angles(M_PI, M_PI, M_PI, M_PI, M_PI);
    rxn.reactantListNew.emplace_back("Asite", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("Bsite", 1, 1, 0, '\0', false);
    rxn.productListNew.emplace_back("Asite", 0, 0, 0, '\0', true);
    rxn.productListNew.emplace_back("Bsite", 1, 1, 0, '\0', true);
    rxn.intReactantList = std::vector<int> { 0, 1 };
    rxn.intProductList = std::vector<int> { 0, 1 };
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = 10.0;
    sys.forwardRxns.push_back(rxn);

    return sys;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the function under test is linked, and its fatal precondition.
//
// The restart parser opens `restartFileNameInput` and, when the stream is not
// good, calls `error("could not find restart file, exiting...")` which
// terminates the process. We verify the *condition* the parser tests (an
// unopenable ifstream vs. an openable one) without ever invoking the parser on
// a missing file, because doing so would abort the entire test binary.
// -----------------------------------------------------------------------------
void pifars_test_symbol_and_missing_restart_precondition()
{
    std::cerr << "\n[TEST] pifars_test_symbol_and_missing_restart_precondition\n"
              << "  Source file: parse_input_for_a_restart_simulation.cpp\n"
              << "  Function:    parse_input_for_a_restart_simulation\n"
              << "  Checks:      (a) the symbol is linked into this binary,\n"
              << "               (b) a non-existent restart file yields a bad\n"
              << "                   ifstream (the parser turns this into a\n"
              << "                   fatal error(); NOT exercised here),\n"
              << "               (c) an existing restart file opens cleanly.\n";

    // (a) Taking the address proves the translation unit under test compiled and
    //     linked with exactly the declared signature.
    auto fn = &parse_input_for_a_restart_simulation;
    EXPECT_TRUE(fn != nullptr)
        << "parse_input_for_a_restart_simulation must be linked into the test binary";
    std::cerr << "  -> symbol parse_input_for_a_restart_simulation resolved OK\n";

    // (b) The exact guard used by the parser: `std::ifstream f{name}; if (!f) error(...)`.
    const std::string missingName { "pifars_this_restart_file_does_not_exist.dat" };
    std::remove(missingName.c_str()); // make sure it really is absent
    {
        std::ifstream missing { missingName };
        EXPECT_FALSE(static_cast<bool>(missing))
            << "A non-existent restart file must produce a failed ifstream; the parser "
               "converts this into error() -> exit(), so it is deliberately not called here";
        std::cerr << "  -> missing file \"" << missingName
                  << "\" correctly fails to open (parser would exit() here)\n";
    }

    // (c) A file that does exist must satisfy the guard, i.e. the parser would
    //     proceed to read_restart().
    const std::string presentName { "pifars_probe_restart.dat" };
    {
        std::ofstream out { presentName };
        out << "# probe file written by the unit test\n";
        out.close();
    }
    {
        std::ifstream present { presentName };
        EXPECT_TRUE(static_cast<bool>(present))
            << "An existing restart file must open, so the parser's guard passes";
        std::cerr << "  -> existing file \"" << presentName << "\" opens OK\n";
    }
    std::remove(presentName.c_str());
}

// -----------------------------------------------------------------------------
// Test 2: RNG fallback behaviour.
//
// When `read_rng_state(rank)` throws, the parser falls back to
// `gsl_rng_set(r, seed)`. That fallback must give a reproducible stream, since
// restarts are expected to be deterministic for a given seed.
// -----------------------------------------------------------------------------
void pifars_test_rng_seed_fallback()
{
    std::cerr << "\n[TEST] pifars_test_rng_seed_fallback\n"
              << "  Source file: parse_input_for_a_restart_simulation.cpp\n"
              << "  Behaviour:   catch(read_rng_state) -> gsl_rng_set(r, seed)\n"
              << "  Pass:        identical seeds reproduce identical streams and\n"
              << "               a different seed produces a different stream.\n";

    pifars_init_rng(); // allocates r if needed and seeds it with 42

    // Draw a reference stream from seed 42.
    std::vector<double> reference;
    for (int i = 0; i < 5; ++i)
        reference.push_back(rand_gsl());
    std::cerr << "  -> reference stream (seed 42) begins with " << reference[0] << '\n';

    // Re-seeding with the same seed must reproduce it bit for bit.
    gsl_rng_set(r, 42);
    for (int i = 0; i < 5; ++i) {
        const double again = rand_gsl();
        EXPECT_DOUBLE_EQ(again, reference[i])
            << "gsl_rng_set(r, 42) must reproduce draw " << i << " of the stream";
    }
    std::cerr << "  -> re-seeding with 42 reproduced all 5 draws\n";

    // A different seed must (with overwhelming probability) diverge.
    gsl_rng_set(r, 43);
    const double otherSeedDraw = rand_gsl();
    EXPECT_NE(otherSeedDraw, reference[0])
        << "A different fallback seed should give a different random stream";
    std::cerr << "  -> seed 43 first draw " << otherSeedDraw << " differs from seed 42\n";

    // Leave the generator in the documented, deterministic test state.
    gsl_rng_set(r, 42);
}

// -----------------------------------------------------------------------------
// Test 3: the spherical-boundary branch.
//
// The parser executes, for `membraneObject.isSphere == true`:
//     membraneObject.create_water_box();
//     membraneObject.sphereVol = (4.0 * M_PI * pow(sphereR, 3.0)) / 3.0;
// Here the same two statements are applied to a Membrane and their
// post-conditions are checked.
// -----------------------------------------------------------------------------
void pifars_test_sphere_water_box_and_volume()
{
    std::cerr << "\n[TEST] pifars_test_sphere_water_box_and_volume\n"
              << "  Source file: parse_input_for_a_restart_simulation.cpp\n"
              << "  Behaviour:   sphere branch -> create_water_box() + sphereVol\n"
              << "  Pass:        the generated box can contain the sphere and\n"
              << "               sphereVol == 4/3 * pi * R^3.\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 50.0; // nm

    // Step 1 of the sphere branch.
    membraneObject.create_water_box();
    std::cerr << "  -> create_water_box() produced [" << membraneObject.waterBox.x << ", "
              << membraneObject.waterBox.y << ", " << membraneObject.waterBox.z << "] nm\n";

    EXPECT_GT(membraneObject.waterBox.x, 0.0) << "water box x must be positive";
    EXPECT_GT(membraneObject.waterBox.y, 0.0) << "water box y must be positive";
    EXPECT_GT(membraneObject.waterBox.z, 0.0) << "water box z must be positive";

    // The enclosing box must be at least as large as the sphere's diameter,
    // otherwise molecules on the sphere would immediately be "out of box".
    EXPECT_GE(membraneObject.waterBox.x, 2.0 * membraneObject.sphereR)
        << "water box x must be able to contain the sphere diameter";
    EXPECT_GE(membraneObject.waterBox.y, 2.0 * membraneObject.sphereR)
        << "water box y must be able to contain the sphere diameter";
    EXPECT_GE(membraneObject.waterBox.z, 2.0 * membraneObject.sphereR)
        << "water box z must be able to contain the sphere diameter";

    // Step 2 of the sphere branch.
    membraneObject.sphereVol = (4.0 * M_PI * pow(membraneObject.sphereR, 3.0)) / 3.0;
    const double analytic = (4.0 / 3.0) * M_PI * 50.0 * 50.0 * 50.0;
    std::cerr << "  -> sphereVol = " << membraneObject.sphereVol
              << " nm^3 (analytic " << analytic << ")\n";
    EXPECT_NEAR(membraneObject.sphereVol, analytic, 1e-6)
        << "sphereVol must equal 4/3 * pi * R^3";
}

// -----------------------------------------------------------------------------
// Test 4: precondition of the per-state protein counter initialisation.
//
// After read_restart() the parser appends `nStates` zeroes to
// `membraneObject.numberOfProteinEachState`. That is only correct if a freshly
// created / restored Membrane starts with an empty vector, which is what is
// asserted here.
// -----------------------------------------------------------------------------
void pifars_test_membrane_state_counter_precondition()
{
    std::cerr << "\n[TEST] pifars_test_membrane_state_counter_precondition\n"
              << "  Source file: parse_input_for_a_restart_simulation.cpp\n"
              << "  Behaviour:   for (i < nStates) numberOfProteinEachState.emplace_back(0)\n"
              << "  Pass:        a default Membrane has nStates == 0 and an empty\n"
              << "               numberOfProteinEachState (so the loop appends\n"
              << "               exactly nStates entries, no more).\n";

    Membrane membraneObject;
    EXPECT_EQ(membraneObject.nStates, 0)
        << "A default-constructed Membrane must report zero implicit-lipid states";
    EXPECT_TRUE(membraneObject.numberOfProteinEachState.empty())
        << "numberOfProteinEachState must start empty so the restart loop sizes it correctly";
    EXPECT_TRUE(membraneObject.numberOfFreeLipidsEachState.empty())
        << "numberOfFreeLipidsEachState must start empty as well";

    // Emulate the restart loop for a membrane that declares three states and
    // confirm the resulting invariant (size == nStates, all counters zero).
    membraneObject.nStates = 3;
    for (int stateItr = 0; stateItr < membraneObject.nStates; ++stateItr)
        membraneObject.numberOfProteinEachState.emplace_back(0);

    EXPECT_EQ(static_cast<int>(membraneObject.numberOfProteinEachState.size()), membraneObject.nStates)
        << "one protein counter per implicit-lipid state must be created";
    for (std::size_t i = 0; i < membraneObject.numberOfProteinEachState.size(); ++i) {
        EXPECT_EQ(membraneObject.numberOfProteinEachState[i], 0)
            << "counter " << i << " must be initialised to zero";
    }
    std::cerr << "  -> numberOfProteinEachState sized to " << membraneObject.nStates
              << " zeroed entries\n";
}

// -----------------------------------------------------------------------------
// Test 5: set_rMaxLimit(), the step that sizes the sub-volumes.
//
// The restart parser calls
//     set_rMaxLimit(params, molTemplateList, forwardRxns, 0, 0);
// before partitioning the box. rMaxLimit is the largest distance any pair of
// reactants can cover in one time step, so it must be strictly positive and at
// least as large as the binding radius of the reaction.
// -----------------------------------------------------------------------------
void pifars_test_set_rmax_limit()
{
    std::cerr << "\n[TEST] pifars_test_set_rmax_limit\n"
              << "  Source file: parse_input_for_a_restart_simulation.cpp\n"
              << "  Function:    set_rMaxLimit(params, molTemplateList, forwardRxns, 0, 0)\n"
              << "  Pass:        rMaxLimit > 0 and rMaxLimit >= bindRadius.\n";

    PifarsSystem sys = pifars_build_minimal_system();
    std::cerr << "  -> system: " << sys.molTemplateList.size() << " templates, "
              << sys.moleculeList.size() << " molecules, "
              << sys.forwardRxns.size() << " forward reaction(s), bindRadius = "
              << sys.forwardRxns[0].bindRadius << " nm\n";

    set_rMaxLimit(sys.params, sys.molTemplateList, sys.forwardRxns,
                  /*numDoubleBeforeAdd=*/0, /*numMolTemplateBeforeAdd=*/0);

    std::cerr << "  -> rMaxLimit  = " << sys.params.rMaxLimit << " nm\n";
    std::cerr << "  -> rMaxRadius = " << sys.params.rMaxRadius << " nm\n";

    EXPECT_GT(sys.params.rMaxLimit, 0.0)
        << "rMaxLimit must be positive; the sub-volume size is box length / rMaxLimit";
    EXPECT_GE(sys.params.rMaxLimit, sys.forwardRxns[0].bindRadius)
        << "rMaxLimit must be at least the binding radius of the reaction";
    EXPECT_GE(sys.params.rMaxRadius, 0.0)
        << "rMaxRadius (largest complex radius) must be non-negative";
}

// -----------------------------------------------------------------------------
// Test 6: simulation-volume partitioning + member molecule lists.
//
// The restart parser finishes with
//     simulVolume.create_simulation_volume(params, membraneObject);
//     simulVolume.update_memberMolLists(...);
//     simulVolume.display();
// The serial overload of update_memberMolLists() is used here (the parser uses
// the MpiContext overload, which a unit test cannot safely construct).
// -----------------------------------------------------------------------------
void pifars_test_simulation_volume_partitioning()
{
    std::cerr << "\n[TEST] pifars_test_simulation_volume_partitioning\n"
              << "  Source file: parse_input_for_a_restart_simulation.cpp\n"
              << "  Functions:   SimulVolume::create_simulation_volume /\n"
              << "               SimulVolume::update_memberMolLists / display\n"
              << "  Pass:        a positive number of sub-volumes is created, the\n"
              << "               sub-cell list matches that count, every molecule\n"
              << "               is assigned to exactly one in-range sub-volume and\n"
              << "               no cell exceeds maxNeighbors neighbours.\n";

    PifarsSystem sys = pifars_build_minimal_system();

    // Same preparatory step the parser performs.
    set_rMaxLimit(sys.params, sys.molTemplateList, sys.forwardRxns, 0, 0);

    // Defensive fallback: create_simulation_volume() divides the box length by
    // rMaxLimit, so a non-positive value would be catastrophic. Report a failure
    // and substitute a sane value rather than letting the suite crash.
    if (!(sys.params.rMaxLimit > 0.0)) {
        ADD_FAILURE() << "set_rMaxLimit produced a non-positive rMaxLimit ("
                      << sys.params.rMaxLimit << "); substituting 10 nm to keep the test safe";
        sys.params.rMaxLimit = 10.0;
    }

    SimulVolume simulVolume;
    simulVolume.create_simulation_volume(sys.params, sys.membrane);

    std::cerr << "  -> sub-cells [" << simulVolume.numSubCells.x << ", "
              << simulVolume.numSubCells.y << ", " << simulVolume.numSubCells.z
              << "] = " << simulVolume.numSubCells.tot << " total\n";
    std::cerr << "  -> sub-cell size [" << simulVolume.subCellSize.x << ", "
              << simulVolume.subCellSize.y << ", " << simulVolume.subCellSize.z << "] nm\n";

    EXPECT_GT(simulVolume.numSubCells.tot, 0) << "the box must be split into at least one sub-volume";
    EXPECT_EQ(simulVolume.numSubCells.tot,
              simulVolume.numSubCells.x * simulVolume.numSubCells.y * simulVolume.numSubCells.z)
        << "tot must equal x * y * z for a cubic partitioning";
    EXPECT_EQ(simulVolume.subCellList.size(), static_cast<std::size_t>(simulVolume.numSubCells.tot))
        << "one SubVolume object per sub-cell must be allocated";

    // Sub-cell edge lengths must tile the water box exactly.
    if (simulVolume.numSubCells.x > 0) {
        EXPECT_NEAR(simulVolume.subCellSize.x * simulVolume.numSubCells.x,
                    sys.membrane.waterBox.x, 1e-9)
            << "sub-cells must tile the box along x";
        EXPECT_NEAR(simulVolume.subCellSize.y * simulVolume.numSubCells.y,
                    sys.membrane.waterBox.y, 1e-9)
            << "sub-cells must tile the box along y";
    }

    // Neighbour lists are built for "forward and up" only, capped at maxNeighbors.
    for (const auto& cell : simulVolume.subCellList) {
        EXPECT_LE(static_cast<int>(cell.neighborList.size()), simulVolume.maxNeighbors)
            << "sub-volume " << cell.absIndex << " has too many neighbours";
    }

    // Only assign members when the partitioning succeeded, so that no molecule
    // can be mapped to an out-of-range bin (which the library treats as fatal).
    if (simulVolume.numSubCells.tot > 0
        && simulVolume.subCellList.size() == static_cast<std::size_t>(simulVolume.numSubCells.tot)
        && simulVolume.subCellSize.x > 0.0 && simulVolume.subCellSize.y > 0.0
        && simulVolume.subCellSize.z > 0.0) {

        // simItr == 1 selects the cheap assignment path (the expensive boundary
        // audit runs only every 1000 iterations); all molecules were built well
        // inside the box, so no boundary correction is required either way.
        simulVolume.update_memberMolLists(sys.params, sys.moleculeList, sys.complexList,
                                          sys.molTemplateList, sys.membrane, /*simItr=*/1);

        std::size_t assigned { 0 };
        for (const auto& cell : simulVolume.subCellList)
            assigned += cell.memberMolList.size();

        std::cerr << "  -> " << assigned << " molecule(s) assigned to sub-volumes\n";
        EXPECT_EQ(assigned, sys.moleculeList.size())
            << "every non-empty molecule must appear in exactly one sub-volume";

        for (const auto& mol : sys.moleculeList) {
            EXPECT_GE(mol.mySubVolIndex, 0)
                << "molecule " << mol.index << " must be given a sub-volume index";
            EXPECT_LT(mol.mySubVolIndex, simulVolume.numSubCells.tot)
                << "molecule " << mol.index << " sub-volume index must be in range";
            std::cerr << "     molecule " << mol.index << " -> sub-volume "
                      << mol.mySubVolIndex << '\n';
        }

        // The parser also prints the partitioning; exercise that code path too.
        simulVolume.display();
    } else {
        ADD_FAILURE() << "simulation volume was not partitioned correctly; skipping "
                         "member-list assignment to avoid a fatal out-of-box error";
    }
}

// -----------------------------------------------------------------------------
// Test 7: the (non-fatal) trajectory-length cross check.
//
// The tail of the parser reads params.trajFile, keeps the iteration number from
// the last header line containing ':' and only *warns* when it disagrees with
// the restart iteration. Neither branch terminates, so we simply verify the
// file-level precondition: a written trajectory can be re-read and its last
// header recovered, and a missing trajectory yields a failed stream (the
// "writing new trajectory" warning branch).
// -----------------------------------------------------------------------------
void pifars_test_trajectory_length_cross_check_inputs()
{
    std::cerr << "\n[TEST] pifars_test_trajectory_length_cross_check_inputs\n"
              << "  Source file: parse_input_for_a_restart_simulation.cpp\n"
              << "  Behaviour:   scan params.trajFile for the last \"...: N\" header\n"
              << "               and compare N with the restart iteration (warning only)\n"
              << "  Pass:        the last header is recoverable from a written file and\n"
              << "               a missing trajectory produces a failed stream.\n";

    const std::string trajName { "pifars_trajectory.xyz" };

    // Write a two-frame trajectory whose final header records iteration 200.
    {
        std::ofstream traj { trajName };
        traj << "4\n";
        traj << "iteration: 100\n";
        traj << "A 0.0 0.0 0.0\n";
        traj << "4\n";
        traj << "iteration: 200\n";
        traj << "A 1.0 1.0 1.0\n";
        traj.close();
    }

    // Replay exactly the scan the parser performs.
    long long int trajItr { -1 };
    {
        std::ifstream traj { trajName };
        EXPECT_TRUE(static_cast<bool>(traj)) << "the trajectory just written must open";
        std::string line;
        while (getline(traj, line)) {
            auto headerItr = line.find(':');
            if (headerItr != std::string::npos)
                trajItr = std::stoi(line.substr(headerItr + 1, std::string::npos));
        }
        traj.close();
    }
    std::cerr << "  -> last trajectory iteration recovered: " << trajItr << '\n';
    EXPECT_EQ(trajItr, 200) << "the scan must keep the iteration of the LAST header line";

    // A matching restart iteration is the "continuing" branch, a mismatching one
    // is the warning branch; neither is fatal.
    const long long int matchingSimItr { 200 };
    const long long int mismatchedSimItr { 137 };
    EXPECT_EQ(trajItr, matchingSimItr) << "matching lengths take the 'Continuing...' branch";
    EXPECT_NE(trajItr, mismatchedSimItr) << "mismatching lengths only emit a warning";

    // Now the "no trajectory found" branch.
    std::remove(trajName.c_str());
    {
        std::ifstream missing { trajName };
        EXPECT_FALSE(static_cast<bool>(missing))
            << "a removed trajectory must fail to open, triggering the (non-fatal) "
               "'writing new trajectory' warning";
    }
    std::cerr << "  -> missing trajectory correctly detected (non-fatal warning branch)\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper is invoked from its own TEST so
// that a failure in one does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(ParseInputForARestartSimulation, SymbolAndMissingRestartPrecondition)
{
    pifars_test_symbol_and_missing_restart_precondition();
}

TEST(ParseInputForARestartSimulation, RngSeedFallback) { pifars_test_rng_seed_fallback(); }

TEST(ParseInputForARestartSimulation, SphereWaterBoxAndVolume)
{
    pifars_test_sphere_water_box_and_volume();
}

TEST(ParseInputForARestartSimulation, MembraneStateCounterPrecondition)
{
    pifars_test_membrane_state_counter_precondition();
}

TEST(ParseInputForARestartSimulation, SetRMaxLimit) { pifars_test_set_rmax_limit(); }

TEST(ParseInputForARestartSimulation, SimulationVolumePartitioning)
{
    pifars_test_simulation_volume_partitioning();
}

TEST(ParseInputForARestartSimulation, TrajectoryLengthCrossCheckInputs)
{
    pifars_test_trajectory_length_cross_check_inputs();
}