/*! \file test_parse_input_for_a_new_simulation.cpp
 *
 * ### Unit tests for src/parser/parse_input_for_a_new_simulation.cpp
 *
 * The file under test contains exactly one function:
 *
 *     void parse_input_for_a_new_simulation(std::string paramFile,
 *              Parameters&, std::map<std::string,int>&,
 *              std::vector<ForwardRxn>&, std::vector<BackRxn>&,
 *              std::vector<CreateDestructRxn>&, std::vector<MolTemplate>&,
 *              Membrane&, std::vector<Molecule>&, std::vector<Complex>&,
 *              SimulVolume&, std::string, int, long long int,
 *              MpiContext&, std::string, std::string);
 *
 * It is a *driver*: it chains together parse_input(), generate_coordinates(),
 * write_psf(), initialize_paramters_for_implicitlipid_model(),
 * initialize_states(), set_rMaxLimit(), SimulVolume::create_simulation_volume(),
 * write_traj() and write_transition(), and in between performs a handful of
 * small, self-contained bookkeeping steps.
 *
 * ### Why the driver itself is not invoked here
 *
 * Calling parse_input_for_a_new_simulation() directly is *not* safe inside a
 * shared gtest binary:
 *
 *   1. parse_input() reads a parameter file plus one `<molname>.mol` file per
 *      molecule from the current working directory.  Every malformed/missing
 *      token path in the parser ends in `error()` / `exit(1)`, which would tear
 *      down the whole test executable (not just this test).
 *   2. It mutates program-wide statics that other translation units in this
 *      suite depend on (MolTemplate::numMolTypes, MolTemplate::numEachMolType,
 *      Molecule::numberOfMolecules, Complex::numberOfComplexes,
 *      RxnBase::numberOfRxns, Interface::State::totalNumOfStates, ...).
 *   3. It requires a fully populated MpiContext (only defined for the parallel
 *      build) for SimulVolume::update_memberMolLists().
 *
 * So instead each *observable step* the driver performs is exercised here in
 * isolation:
 *   - the implicit-lipid state bookkeeping block (Membrane::nStates and the
 *     per-state free-lipid / protein counters),
 *   - the "implicit lipid must be molecule type 0" validation predicate,
 *   - MolTemplate::numMolTypes assignment and the moleculeList/complexList
 *     reservation arithmetic,
 *   - the spherical-boundary water box + sphere volume set-up (calls the real
 *     Membrane::create_water_box()),
 *   - the observables-file header block (a faithful transcription, so that the
 *     resulting on-disk format is pinned down),
 *   - the transition-matrix allocation loop followed by a real
 *     write_transition() call,
 *   - a real write_traj() call on a minimal system,
 *   - a real SimulVolume::create_simulation_volume() +
 *     update_memberMolLists() partitioning of the simulation box.
 *
 * Every static that is touched is saved and restored so that the rest of the
 * suite is unaffected.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "io/io.hpp"

// -----------------------------------------------------------------------------
// Local helpers.  Everything is file-local and prefixed with `pifans_`
// (Parse Input For A New Simulation) so no symbol can collide with another
// translation unit of the combined test binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Read an entire file into a std::string (empty string if absent). */
std::string pifans_slurp(const std::string& fileName)
{
    std::ifstream in(fileName);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

/*! \brief Count the number of '\n' terminated lines in a blob of text. */
std::size_t pifans_count_lines(const std::string& text)
{
    std::size_t count { 0 };
    for (char c : text) {
        if (c == '\n')
            ++count;
    }
    return count;
}

/*! \brief Faithful transcription of the observables-header block of the
 *         function under test.
 *
 * The production code writes:
 *   - nothing at all when the observables map is empty,
 *   - a header line plus a "0,0" data line when there is exactly one observable,
 *   - only a header line when there is more than one observable (the "0,..."
 *     data row is sent to std::cout, *not* to the file -- reproduced faithfully
 *     below so the assertions document the real on-disk result).
 */
void pifans_write_observables_header(const std::map<std::string, int>& observablesList,
    const std::string& observablesFileName)
{
    std::ofstream observablesFile { observablesFileName };
    if (observablesList.size() == 1) {
        observablesFile << "Time (s)," << observablesList.begin()->first << '\n';
        observablesFile << "0,0\n";
    } else if (observablesList.size() > 1) {
        observablesFile << "Time (s)";
        for (auto obsItr = observablesList.begin(); obsItr != observablesList.end(); ++obsItr) {
            observablesFile << ',' << obsItr->first;
        }
        // NOTE: production sends the following row to std::cout, not to the file.
        std::cout << '\n' << "0";
        for (auto obsItr = observablesList.begin(); obsItr != observablesList.end(); ++obsItr) {
            std::cout << ',' << obsItr->second;
        }
        observablesFile << '\n' << std::flush;
    }
    observablesFile.close();
}

/*! \brief Build a fully initialised MolTemplate with one interface.
 *
 * Every field that the collaborators dereference (name, indices, mass, radius,
 * diffusion constants, interface state list) is populated so that no callee can
 * read uninitialised memory.
 *
 * \param[in] molName     name of the molecule type
 * \param[in] molTypeIndex index of this template in molTemplateList
 * \param[in] copies      starting copy number
 * \param[in] numStates   number of states on the single interface
 */
MolTemplate pifans_make_template(const std::string& molName, int molTypeIndex, int copies, int numStates)
{
    MolTemplate temp;
    temp.molName = molName;
    temp.molTypeIndex = molTypeIndex;
    temp.copies = copies;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord(10.0, 10.0, 10.0);
    temp.Dr = Coord(0.1, 0.1, 0.1);
    temp.comCoord = Coord(0.0, 0.0, 0.0);

    // Build the state list of the single interface "a".
    std::vector<Interface::State> stateList;
    for (int stateItr = 0; stateItr < numStates; ++stateItr) {
        // Interface::State(char iden, int index)
        stateList.emplace_back(static_cast<char>('U' + stateItr), stateItr);
    }

    Interface iface { "a", stateList, Coord(0.0, 0.0, 1.0) };
    iface.index = 0;
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a fully initialised Molecule that matches a MolTemplate. */
Molecule pifans_make_molecule(int index, const Coord& com, const MolTemplate& temp)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = temp.molTypeIndex;
    mol.myComIndex = index;
    mol.mySubVolIndex = -1;
    mol.mass = temp.mass;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.trajStatus = TrajStatus::none;

    for (std::size_t ifaceItr = 0; ifaceItr < temp.interfaceList.size(); ++ifaceItr) {
        Molecule::Iface iface;
        iface.coord = com + temp.interfaceList[ifaceItr].iCoord;
        iface.relIndex = static_cast<int>(ifaceItr);
        iface.index = temp.interfaceList[ifaceItr].stateList[0].index;
        iface.stateIden = temp.interfaceList[ifaceItr].stateList[0].iden;
        iface.stateIndex = 0;
        iface.molTypeIndex = temp.molTypeIndex;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Build a Complex owning exactly one member Molecule.
 *
 * Built by hand (rather than with Complex(const Molecule&, const MolTemplate&))
 * so that the Complex::maxID / Complex::numberOfComplexes statics used by the
 * rest of the suite are left untouched.
 */
Complex pifans_make_complex(int index, const Molecule& mol, const MolTemplate& temp, int numMolTypes)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.comCoord = mol.comCoord;
    com.D = temp.D;
    com.Dr = temp.Dr;
    com.radius = temp.radius;
    com.mass = mol.mass;
    com.isEmpty = false;
    com.OnSurface = false;
    com.memberList.push_back(mol.index);
    com.numEachMol.resize(numMolTypes, 0);
    if (mol.molTypeIndex >= 0 && mol.molTypeIndex < numMolTypes)
        ++com.numEachMol[mol.molTypeIndex];
    com.lastNumberUpdateItrEachMol.resize(numMolTypes, 0);
    return com;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the implicit-lipid bookkeeping performed right after parse_input().
//
//   for (auto& t : molTemplateList)
//       if (t.isImplicitLipid) { membraneObject.nStates = t.interfaceList[0].stateList.size(); break; }
//   for (i < nStates) { numberOfFreeLipidsEachState.emplace_back(0);
//                       numberOfProteinEachState.emplace_back(0); }
//   for (auto& t : molTemplateList)
//       if (t.isImplicitLipid && t.molTypeIndex != 0) error(...)
//
// Pass criteria: nStates equals the number of states on interface 0 of the
// implicit lipid, both per-state counter vectors are that long and zeroed, and
// the "implicit lipid must be first" predicate correctly separates the legal
// and illegal layouts.
// -----------------------------------------------------------------------------
void test_pifans_implicit_lipid_state_setup()
{
    std::cerr << "\n[TEST] parse_input_for_a_new_simulation.cpp -> implicit-lipid state bookkeeping\n"
              << "  Checking Membrane::nStates derivation and the per-state counter vectors.\n";

    // --- molTemplateList with the implicit lipid first (the legal layout) ----
    std::vector<MolTemplate> molTemplateList;
    MolTemplate implicitLipid = pifans_make_template("IL", 0, 100, 3); // 3 states
    implicitLipid.isImplicitLipid = true;
    implicitLipid.isLipid = true;
    molTemplateList.push_back(implicitLipid);
    molTemplateList.push_back(pifans_make_template("A", 1, 10, 1)); // ordinary protein

    Membrane membraneObject;
    std::cerr << "  -> Membrane::nStates default before the loop = " << membraneObject.nStates << '\n';
    EXPECT_EQ(membraneObject.nStates, 0) << "Membrane::nStates must default to 0";

    // Block under test (verbatim transcription).
    for (auto& molTemplateTmp : molTemplateList) {
        if (molTemplateTmp.isImplicitLipid == true) {
            membraneObject.nStates = static_cast<int>(molTemplateTmp.interfaceList[0].stateList.size());
            break;
        }
    }
    std::cerr << "  -> nStates after scanning the templates = " << membraneObject.nStates << '\n';
    EXPECT_EQ(membraneObject.nStates, 3) << "nStates must equal the state count of the implicit lipid's iface 0";

    for (int tmpStateIndex = 0; tmpStateIndex < membraneObject.nStates; tmpStateIndex++) {
        membraneObject.numberOfFreeLipidsEachState.emplace_back(0);
        membraneObject.numberOfProteinEachState.emplace_back(0);
    }
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState.size(), 3u)
        << "one free-lipid counter must be created per implicit-lipid state";
    EXPECT_EQ(membraneObject.numberOfProteinEachState.size(), 3u)
        << "one protein counter must be created per implicit-lipid state";
    for (std::size_t i = 0; i < membraneObject.numberOfFreeLipidsEachState.size(); ++i) {
        EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[i], 0) << "free-lipid counters start at zero";
        EXPECT_EQ(membraneObject.numberOfProteinEachState[i], 0) << "protein counters start at zero";
    }

    // --- the "implicit lipid must be molecule type 0" validation -------------
    // The production code calls error() (which exits) on a violation, so only
    // the *predicate* is exercised here -- never the error path itself.
    bool violationDetected { false };
    for (auto& tempMolTemplate : molTemplateList) {
        if (tempMolTemplate.isImplicitLipid == true && tempMolTemplate.molTypeIndex != 0)
            violationDetected = true;
    }
    std::cerr << "  -> legal layout (implicit lipid at index 0): violation = " << std::boolalpha
              << violationDetected << '\n';
    EXPECT_FALSE(violationDetected) << "implicit lipid at molTypeIndex 0 must be accepted";

    // Now move the implicit lipid off index 0 and re-run the predicate.
    molTemplateList[0].molTypeIndex = 2;
    violationDetected = false;
    for (auto& tempMolTemplate : molTemplateList) {
        if (tempMolTemplate.isImplicitLipid == true && tempMolTemplate.molTypeIndex != 0)
            violationDetected = true;
    }
    std::cerr << "  -> illegal layout (implicit lipid at index 2): violation = " << violationDetected << '\n';
    EXPECT_TRUE(violationDetected) << "an implicit lipid that is not molecule type 0 must be flagged";

    // --- a system with no implicit lipid at all leaves nStates untouched -----
    std::vector<MolTemplate> noLipidList { pifans_make_template("A", 0, 5, 2) };
    Membrane noLipidMembrane;
    for (auto& molTemplateTmp : noLipidList) {
        if (molTemplateTmp.isImplicitLipid == true) {
            noLipidMembrane.nStates = static_cast<int>(molTemplateTmp.interfaceList[0].stateList.size());
            break;
        }
    }
    EXPECT_EQ(noLipidMembrane.nStates, 0) << "without an implicit lipid nStates must remain 0";
    EXPECT_TRUE(noLipidMembrane.numberOfFreeLipidsEachState.empty())
        << "without an implicit lipid no per-state counters are created";
}

// -----------------------------------------------------------------------------
// Test 2: MolTemplate::numMolTypes assignment and the container reservation.
//
//   MolTemplate::numMolTypes = molTemplateList.size();
//   unsigned long reservation{};
//   for (auto& molTemp : molTemplateList) reservation += molTemp.copies;
//   moleculeList.reserve(reservation); complexList.reserve(reservation);
//
// Pass criteria: numMolTypes matches the template count, the reservation is the
// sum of the copy numbers, and both containers end up with at least that
// capacity while still being logically empty.
// -----------------------------------------------------------------------------
void test_pifans_moltype_count_and_reservation()
{
    std::cerr << "\n[TEST] parse_input_for_a_new_simulation.cpp -> numMolTypes + list reservation\n";

    // Save the process-wide static so the rest of the suite is unaffected.
    const unsigned savedNumMolTypes = MolTemplate::numMolTypes;

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(pifans_make_template("IL", 0, 100, 1));
    molTemplateList.push_back(pifans_make_template("A", 1, 7, 1));
    molTemplateList.push_back(pifans_make_template("B", 2, 3, 1));

    MolTemplate::numMolTypes = molTemplateList.size();
    std::cerr << "  -> MolTemplate::numMolTypes set to " << MolTemplate::numMolTypes << '\n';
    EXPECT_EQ(MolTemplate::numMolTypes, 3u) << "numMolTypes must equal molTemplateList.size()";

    unsigned long reservation {};
    for (auto& molTemp : molTemplateList) {
        reservation += molTemp.copies;
    }
    std::cerr << "  -> reservation (sum of copies 100 + 7 + 3) = " << reservation << '\n';
    EXPECT_EQ(reservation, 110ul) << "reservation must be the sum of every template's copy number";

    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    moleculeList.reserve(reservation);
    complexList.reserve(reservation);

    EXPECT_GE(moleculeList.capacity(), reservation) << "moleculeList must reserve at least `reservation` slots";
    EXPECT_GE(complexList.capacity(), reservation) << "complexList must reserve at least `reservation` slots";
    EXPECT_TRUE(moleculeList.empty()) << "reserve() must not create any Molecule";
    EXPECT_TRUE(complexList.empty()) << "reserve() must not create any Complex";

    // Restore the static.
    MolTemplate::numMolTypes = savedNumMolTypes;
    std::cerr << "  -> restored MolTemplate::numMolTypes to " << MolTemplate::numMolTypes << '\n';
}

// -----------------------------------------------------------------------------
// Test 3: the spherical-boundary branch.
//
//   if (membraneObject.isSphere) {
//       membraneObject.create_water_box();
//       membraneObject.sphereVol = (4.0 * M_PI * pow(sphereR, 3.0)) / 3.0;
//   }
//
// Pass criteria: sphereVol is the analytic volume of the sphere, and
// create_water_box() (the real implementation from parse_input.cpp) produces a
// water box with strictly positive dimensions.  A non-spherical membrane must
// be left completely alone.
// -----------------------------------------------------------------------------
void test_pifans_sphere_water_box_and_volume()
{
    std::cerr << "\n[TEST] parse_input_for_a_new_simulation.cpp -> spherical boundary set-up\n";

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 50.0;

    if (membraneObject.isSphere) {
        membraneObject.create_water_box();
        membraneObject.sphereVol = (4.0 * M_PI * pow(membraneObject.sphereR, 3.0)) / 3.0;
    }

    const double expectedVol = (4.0 / 3.0) * M_PI * 50.0 * 50.0 * 50.0;
    std::cerr << "  -> sphereR = " << membraneObject.sphereR << ", sphereVol = " << membraneObject.sphereVol
              << " (expected " << expectedVol << ")\n";
    EXPECT_NEAR(membraneObject.sphereVol, expectedVol, 1e-6)
        << "sphereVol must be 4/3 pi R^3";

    std::cerr << "  -> create_water_box() produced [" << membraneObject.waterBox.x << ", "
              << membraneObject.waterBox.y << ", " << membraneObject.waterBox.z
              << "], volume = " << membraneObject.waterBox.volume << '\n';
    EXPECT_GT(membraneObject.waterBox.x, 0.0) << "the enclosing water box must have a positive x extent";
    EXPECT_GT(membraneObject.waterBox.y, 0.0) << "the enclosing water box must have a positive y extent";
    EXPECT_GT(membraneObject.waterBox.z, 0.0) << "the enclosing water box must have a positive z extent";

    // A non-spherical membrane must skip the branch entirely.
    Membrane boxMembrane;
    boxMembrane.isSphere = false;
    boxMembrane.waterBox = Membrane::WaterBox(std::vector<double> { 20.0, 20.0, 20.0 });
    const double preservedVolume = boxMembrane.waterBox.volume;
    if (boxMembrane.isSphere) { // deliberately false
        boxMembrane.create_water_box();
        boxMembrane.sphereVol = (4.0 * M_PI * pow(boxMembrane.sphereR, 3.0)) / 3.0;
    }
    EXPECT_DOUBLE_EQ(boxMembrane.sphereVol, 0.0) << "a box membrane must keep sphereVol == 0";
    EXPECT_DOUBLE_EQ(boxMembrane.waterBox.volume, preservedVolume)
        << "a box membrane's water box must not be recreated";
    EXPECT_DOUBLE_EQ(boxMembrane.waterBox.volume, 8000.0) << "20 x 20 x 20 = 8000 nm^3";
}

// -----------------------------------------------------------------------------
// Test 4: the observables file header block.
//
// Pass criteria (exactly what the production branches do):
//   - 0 observables  -> the file is created but stays empty,
//   - 1 observable   -> "Time (s),<name>\n0,0\n"  (two lines),
//   - >1 observables -> only the header line ends up on disk; the "0,..." data
//                       row is sent to std::cout by the production code.
// -----------------------------------------------------------------------------
void test_pifans_observables_file_writing()
{
    std::cerr << "\n[TEST] parse_input_for_a_new_simulation.cpp -> observables file header\n";

    const std::string fileName { "unittest_pifans_observables.csv" };

    // --- case 1: no observables at all --------------------------------------
    {
        std::map<std::string, int> observablesList;
        pifans_write_observables_header(observablesList, fileName);
        const std::string content = pifans_slurp(fileName);
        std::cerr << "  -> 0 observables: file length = " << content.size() << " bytes\n";
        EXPECT_TRUE(content.empty()) << "with no observables the file must be created but left empty";
    }

    // --- case 2: exactly one observable -------------------------------------
    {
        std::map<std::string, int> observablesList;
        observablesList["dimer"] = 0;
        pifans_write_observables_header(observablesList, fileName);
        const std::string content = pifans_slurp(fileName);
        std::cerr << "  -> 1 observable: file contents =\n\"" << content << "\"\n";
        EXPECT_EQ(content, std::string("Time (s),dimer\n0,0\n"))
            << "a single observable must produce a header line and a zeroed data line";
        EXPECT_EQ(pifans_count_lines(content), 2u) << "single-observable file must contain exactly 2 lines";
    }

    // --- case 3: more than one observable ------------------------------------
    {
        std::map<std::string, int> observablesList;
        observablesList["obsA"] = 3;
        observablesList["obsB"] = 7;
        pifans_write_observables_header(observablesList, fileName);
        const std::string content = pifans_slurp(fileName);
        std::cerr << "  -> 2 observables: file contents =\n\"" << content << "\"\n";
        // std::map iterates in sorted key order, so the column order is obsA, obsB.
        EXPECT_EQ(content, std::string("Time (s),obsA,obsB\n"))
            << "with >1 observables only the header reaches the file (the data row goes to std::cout)";
        EXPECT_EQ(pifans_count_lines(content), 1u)
            << "multi-observable file must contain exactly 1 line -- documents the std::cout quirk";
        EXPECT_NE(content.find("obsA"), std::string::npos) << "header must list obsA";
        EXPECT_NE(content.find("obsB"), std::string::npos) << "header must list obsB";
    }

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 5: the transition-matrix allocation loop plus a real write_transition().
//
//   for (auto& molTemp : molTemplateList)
//       if (molTemp.countTransition) {
//           molTemp.transitionMatrix.resize(N);
//           molTemp.lifeTime.resize(N);
//           for (i < N) molTemp.transitionMatrix[i].resize(N);
//       }
//
// Pass criteria: only templates that request transition counting get storage;
// transitionMatrix becomes a zeroed N x N matrix; lifeTime gets N rows but the
// rows are deliberately left empty by the production code; write_transition()
// then succeeds and emits a non-empty file.
// -----------------------------------------------------------------------------
void test_pifans_transition_matrix_allocation()
{
    std::cerr << "\n[TEST] parse_input_for_a_new_simulation.cpp -> transition matrix allocation\n";

    const unsigned savedNumMolTypes = MolTemplate::numMolTypes;

    std::vector<MolTemplate> molTemplateList;
    MolTemplate counted = pifans_make_template("A", 0, 4, 1);
    counted.countTransition = true;
    counted.transitionMatrixSize = 4;
    MolTemplate uncounted = pifans_make_template("B", 1, 4, 1);
    uncounted.countTransition = false;
    uncounted.transitionMatrixSize = 4;
    molTemplateList.push_back(counted);
    molTemplateList.push_back(uncounted);
    MolTemplate::numMolTypes = molTemplateList.size();

    // Block under test (verbatim transcription).
    for (auto& molTemp : molTemplateList) {
        if (molTemp.countTransition == true) {
            molTemp.transitionMatrix.resize(molTemp.transitionMatrixSize);
            molTemp.lifeTime.resize(molTemp.transitionMatrixSize);
            for (int indexOne = 0; indexOne < molTemp.transitionMatrixSize; ++indexOne) {
                molTemp.transitionMatrix[indexOne].resize(molTemp.transitionMatrixSize);
            }
        }
    }

    std::cerr << "  -> template A (countTransition=true): transitionMatrix rows = "
              << molTemplateList[0].transitionMatrix.size()
              << ", lifeTime rows = " << molTemplateList[0].lifeTime.size() << '\n';
    EXPECT_EQ(molTemplateList[0].transitionMatrix.size(), 4u) << "transitionMatrix must have N rows";
    EXPECT_EQ(molTemplateList[0].lifeTime.size(), 4u) << "lifeTime must have N rows";
    for (const auto& row : molTemplateList[0].transitionMatrix) {
        EXPECT_EQ(row.size(), 4u) << "every transitionMatrix row must be resized to N columns";
        for (long long int value : row)
            EXPECT_EQ(value, 0) << "a freshly allocated transition matrix must be zeroed";
    }
    // Documented behaviour: the inner vectors of lifeTime are *not* resized.
    EXPECT_TRUE(molTemplateList[0].lifeTime[0].empty())
        << "lifeTime rows are intentionally left empty by the allocation loop";

    std::cerr << "  -> template B (countTransition=false): transitionMatrix rows = "
              << molTemplateList[1].transitionMatrix.size() << '\n';
    EXPECT_TRUE(molTemplateList[1].transitionMatrix.empty())
        << "templates that do not request transition counting must get no storage";
    EXPECT_TRUE(molTemplateList[1].lifeTime.empty())
        << "templates that do not request transition counting must get no lifetime storage";

    // --- now drive the real write_transition() exactly as production does ----
    const std::string transitionFileName { "unittest_pifans_transition.dat" };
    {
        std::ofstream transitionFile { transitionFileName };
        ASSERT_TRUE(transitionFile.is_open()) << "could not open the temporary transition file";
        write_transition(0, transitionFile, molTemplateList);
        transitionFile.close();
    }
    const std::string content = pifans_slurp(transitionFileName);
    std::cerr << "  -> write_transition() wrote " << content.size() << " bytes\n";
    EXPECT_FALSE(content.empty())
        << "write_transition() must emit the iteration-0 record for the counted template";

    std::remove(transitionFileName.c_str());
    MolTemplate::numMolTypes = savedNumMolTypes;
}

// -----------------------------------------------------------------------------
// Test 6: the "write beginning of trajectory" step.
//
//   std::ofstream trajFile{trajFileName};
//   write_traj(0, trajFile, params, moleculeList, molTemplateList, membraneObject);
//
// Pass criteria: the trajectory file is created, is non-empty, and mentions the
// molecule type name of the single molecule in the system.
// -----------------------------------------------------------------------------
void test_pifans_write_traj_output()
{
    std::cerr << "\n[TEST] parse_input_for_a_new_simulation.cpp -> initial trajectory frame\n";

    const unsigned savedNumMolTypes = MolTemplate::numMolTypes;

    Parameters params;
    params.timeStep = 0.1;
    params.nItr = 10;
    params.numMolTypes = 1;
    params.numTotalComplex = 1;
    params.name = "unittest_pifans";

    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });

    std::vector<MolTemplate> molTemplateList { pifans_make_template("Amol", 0, 1, 1) };
    MolTemplate::numMolTypes = molTemplateList.size();

    std::vector<Molecule> moleculeList { pifans_make_molecule(0, Coord(1.0, 2.0, 3.0), molTemplateList[0]) };

    const std::string trajFileName { "unittest_pifans_traj.xyz" };
    {
        std::ofstream trajFile { trajFileName };
        ASSERT_TRUE(trajFile.is_open()) << "could not open the temporary trajectory file";
        write_traj(0, trajFile, params, moleculeList, molTemplateList, membraneObject);
        trajFile.close();
    }

    const std::string content = pifans_slurp(trajFileName);
    std::cerr << "  -> write_traj() wrote " << content.size() << " bytes, first line: \""
              << content.substr(0, content.find('\n')) << "\"\n";
    EXPECT_FALSE(content.empty()) << "the initial trajectory frame must not be empty";
    EXPECT_NE(content.find("Amol"), std::string::npos)
        << "the frame must reference the molecule type name written by write_traj()";
    EXPECT_GE(pifans_count_lines(content), 2u)
        << "an xyz frame contains at least a count/comment header plus one record";

    std::remove(trajFileName.c_str());
    MolTemplate::numMolTypes = savedNumMolTypes;
}

// -----------------------------------------------------------------------------
// Test 7: the "partition the simulation box into sub-boxes" step.
//
//   simulVolume.create_simulation_volume(params, membraneObject);
//   simulVolume.update_memberMolLists(...);
//
// params.rMaxLimit is normally produced by set_rMaxLimit(); it is set here
// explicitly (to a physically sensible 10 nm) because a zero rMaxLimit makes
// the Dimensions constructor divide by zero, which would blow up the whole
// suite rather than fail a single assertion.
//
// Pass criteria: tot == x*y*z, subCellList has exactly `tot` entries, the
// sub-cell size tiles the box, and the single molecule lands in exactly one
// sub-volume whose memberMolList contains its index.
// -----------------------------------------------------------------------------
void test_pifans_simulation_volume_partitioning()
{
    std::cerr << "\n[TEST] parse_input_for_a_new_simulation.cpp -> simulation volume partitioning\n";

    const unsigned savedNumMolTypes = MolTemplate::numMolTypes;

    Parameters params;
    params.timeStep = 0.1;
    params.rMaxLimit = 10.0; // stand-in for the value set_rMaxLimit() would compute
    params.name = "unittest_pifans";

    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });

    std::vector<MolTemplate> molTemplateList { pifans_make_template("Amol", 0, 1, 1) };
    MolTemplate::numMolTypes = molTemplateList.size();

    std::vector<Molecule> moleculeList { pifans_make_molecule(0, Coord(0.0, 0.0, 0.0), molTemplateList[0]) };
    std::vector<Complex> complexList { pifans_make_complex(0, moleculeList[0], molTemplateList[0],
        static_cast<int>(MolTemplate::numMolTypes)) };

    SimulVolume simulVolume;
    simulVolume.create_simulation_volume(params, membraneObject);

    std::cerr << "  -> sub-cell grid = [" << simulVolume.numSubCells.x << ", " << simulVolume.numSubCells.y
              << ", " << simulVolume.numSubCells.z << "], tot = " << simulVolume.numSubCells.tot << '\n';
    std::cerr << "  -> sub-cell size = [" << simulVolume.subCellSize.x << ", " << simulVolume.subCellSize.y
              << ", " << simulVolume.subCellSize.z << "]\n";

    EXPECT_GT(simulVolume.numSubCells.tot, 0) << "the box must be split into at least one sub-volume";
    EXPECT_EQ(simulVolume.numSubCells.tot,
        simulVolume.numSubCells.x * simulVolume.numSubCells.y * simulVolume.numSubCells.z)
        << "tot must be the product of the three grid dimensions";
    EXPECT_EQ(static_cast<int>(simulVolume.subCellList.size()), simulVolume.numSubCells.tot)
        << "one SubVolume object must exist per grid cell";
    EXPECT_NEAR(simulVolume.subCellSize.x, membraneObject.waterBox.x / (simulVolume.numSubCells.x * 1.0), 1e-9)
        << "sub-cells must tile the box exactly in x";
    EXPECT_NEAR(simulVolume.subCellSize.y, membraneObject.waterBox.y / (simulVolume.numSubCells.y * 1.0), 1e-9)
        << "sub-cells must tile the box exactly in y";
    EXPECT_NEAR(simulVolume.subCellSize.z, membraneObject.waterBox.z / (simulVolume.numSubCells.z * 1.0), 1e-9)
        << "sub-cells must tile the box exactly in z";

    // Every sub-volume knows its own absolute index and has a neighbour list
    // that never exceeds maxNeighbors (13 for the cubic decomposition).
    for (std::size_t cellItr = 0; cellItr < simulVolume.subCellList.size(); ++cellItr) {
        ASSERT_EQ(simulVolume.subCellList[cellItr].absIndex, static_cast<int>(cellItr))
            << "SubVolume " << cellItr << " has an inconsistent absIndex";
        EXPECT_LE(static_cast<int>(simulVolume.subCellList[cellItr].neighborList.size()),
            simulVolume.maxNeighbors)
            << "SubVolume " << cellItr << " exceeded the maximum neighbour count";
    }

    // simItr == 1 selects the fast "just assign bins" path (the expensive
    // boundary audit only runs when simItr % 1000 == 0).
    simulVolume.update_memberMolLists(params, moleculeList, complexList, molTemplateList, membraneObject, 1);

    std::cerr << "  -> molecule 0 assigned to sub-volume " << moleculeList[0].mySubVolIndex << '\n';
    EXPECT_GE(moleculeList[0].mySubVolIndex, 0) << "the molecule must be assigned to a sub-volume";
    EXPECT_LT(moleculeList[0].mySubVolIndex, simulVolume.numSubCells.tot)
        << "the assigned sub-volume index must be inside the grid";

    if (moleculeList[0].mySubVolIndex >= 0
        && moleculeList[0].mySubVolIndex < static_cast<int>(simulVolume.subCellList.size())) {
        const auto& members = simulVolume.subCellList[moleculeList[0].mySubVolIndex].memberMolList;
        EXPECT_EQ(members.size(), 1u) << "the owning sub-volume must list exactly the one molecule";
        if (!members.empty())
            EXPECT_EQ(members[0], moleculeList[0].index) << "the listed member must be molecule 0";
    }

    // The molecule must appear in exactly one sub-volume across the whole grid.
    std::size_t totalMembers { 0 };
    for (const auto& subVol : simulVolume.subCellList)
        totalMembers += subVol.memberMolList.size();
    EXPECT_EQ(totalMembers, 1u) << "each molecule must belong to exactly one sub-volume";

    MolTemplate::numMolTypes = savedNumMolTypes;
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is called from its own TEST so that a
// failure in one step still lets every other step run (all assertions are the
// non-fatal EXPECT_* form, except the few ASSERT_* guards protecting a
// subsequent dereference).
// -----------------------------------------------------------------------------
TEST(ParseInputForANewSimulation, ImplicitLipidStateSetup) { test_pifans_implicit_lipid_state_setup(); }
TEST(ParseInputForANewSimulation, MolTypeCountAndReservation) { test_pifans_moltype_count_and_reservation(); }
TEST(ParseInputForANewSimulation, SphereWaterBoxAndVolume) { test_pifans_sphere_water_box_and_volume(); }
TEST(ParseInputForANewSimulation, ObservablesFileWriting) { test_pifans_observables_file_writing(); }
TEST(ParseInputForANewSimulation, TransitionMatrixAllocation) { test_pifans_transition_matrix_allocation(); }
TEST(ParseInputForANewSimulation, WriteTrajOutput) { test_pifans_write_traj_output(); }
TEST(ParseInputForANewSimulation, SimulationVolumePartitioning) { test_pifans_simulation_volume_partitioning(); }