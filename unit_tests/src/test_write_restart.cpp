/*! \file test_write_restart.cpp
 *
 * ### Unit test for ../src/io/write_restart.cpp
 *
 * This test exercises the single function defined in that file:
 *
 *     void write_restart(long long int simItr, std::ofstream& restartFile,
 *                        const Parameters&, const SimulVolume&,
 *                        const std::vector<Molecule>&, const std::vector<Complex>&,
 *                        const std::vector<MolTemplate>&, const std::vector<ForwardRxn>&,
 *                        const std::vector<BackRxn>&, const std::vector<CreateDestructRxn>&,
 *                        const std::vector<TransmissionRxn>&,
 *                        const std::map<std::string,int>&, const Membrane&,
 *                        const copyCounters&)
 *
 * `write_restart()` is a pure serializer: it takes the current simulation state
 * and dumps it, as plain text, into an already-open std::ofstream.  Therefore the
 * strategy of this test is:
 *
 *   1. Build small, fully-controlled simulation state objects (parameters,
 *      molecules, complexes, templates, reactions, observables, counters).
 *   2. Call write_restart() with a temporary file stream.
 *   3. Read the produced file back into memory and verify that the expected
 *      section headers and the expected serialized values are present, on the
 *      expected lines relative to their section header.
 *
 * NOTE ON STATICS: write_restart() also serializes several *static* class
 * members (MolTemplate::numMolTypes, Molecule::numberOfMolecules, ...).  Since
 * these are global to the whole test binary, every test here saves them on entry
 * and restores them on exit through the WrStaticsGuard RAII helper, so that this
 * file cannot disturb other tests in the suite.
 *
 * Verbose progress messages are printed to stderr so that the reader can follow
 * which function is being exercised and what each check verifies.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace keeps every symbol private to this file,
// which also guarantees there can be no name collisions with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! \brief RAII guard that snapshots and restores every static member that
 *         write_restart() reads.
 *
 * Constructing the guard stores the current values; destroying it puts them
 * back.  Tests may therefore freely overwrite the statics.
 */
struct WrStaticsGuard {
    unsigned savedNumMolTypes;
    std::vector<int> savedNumEachMolType;
    std::vector<int> savedAbsToRelIface;
    int savedTotalNumOfStates;
    unsigned savedNumberOfRxns;
    int savedTotRxnSpecies;
    int savedNumberOfMolecules;
    std::vector<int> savedEmptyMolList;
    int savedNumberOfComplexes;
    std::vector<int> savedEmptyComList;
    std::vector<long long int> savedLastUpdateTransition;

    WrStaticsGuard()
        : savedNumMolTypes(MolTemplate::numMolTypes)
        , savedNumEachMolType(MolTemplate::numEachMolType)
        , savedAbsToRelIface(MolTemplate::absToRelIface)
        , savedTotalNumOfStates(Interface::State::totalNumOfStates)
        , savedNumberOfRxns(RxnBase::numberOfRxns)
        , savedTotRxnSpecies(RxnBase::totRxnSpecies)
        , savedNumberOfMolecules(Molecule::numberOfMolecules)
        , savedEmptyMolList(Molecule::emptyMolList)
        , savedNumberOfComplexes(Complex::numberOfComplexes)
        , savedEmptyComList(Complex::emptyComList)
        , savedLastUpdateTransition(Parameters::lastUpdateTransition)
    {
    }

    ~WrStaticsGuard()
    {
        MolTemplate::numMolTypes = savedNumMolTypes;
        MolTemplate::numEachMolType = savedNumEachMolType;
        MolTemplate::absToRelIface = savedAbsToRelIface;
        Interface::State::totalNumOfStates = savedTotalNumOfStates;
        RxnBase::numberOfRxns = savedNumberOfRxns;
        RxnBase::totRxnSpecies = savedTotRxnSpecies;
        Molecule::numberOfMolecules = savedNumberOfMolecules;
        Molecule::emptyMolList = savedEmptyMolList;
        Complex::numberOfComplexes = savedNumberOfComplexes;
        Complex::emptyComList = savedEmptyComList;
        Parameters::lastUpdateTransition = savedLastUpdateTransition;
    }
};

/*! \brief Bundle of every argument write_restart() needs, with sane defaults. */
struct WrScenario {
    long long int simItr { 0 };
    Parameters params {};
    SimulVolume simulVolume {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {};
    std::vector<TransmissionRxn> transmissionRxns {};
    std::map<std::string, int> observablesList {};
    Membrane membraneObject {};
    copyCounters counterArrays {};
};

/*! \brief Fill a Membrane with fully-defined (non garbage) values.
 *
 * Membrane has several members without in-class initializers (nSites,
 * No_free_lipids, No_protein, totalSA) which write_restart() prints, so they
 * must be set explicitly to keep the test deterministic.
 */
void wr_init_membrane(Membrane& membraneObject)
{
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    membraneObject.implicitlipidIndex = -1;
    membraneObject.nSites = 0;
    membraneObject.nStates = 0;
    membraneObject.No_free_lipids = 0;
    membraneObject.No_protein = 0;
    membraneObject.totalSA = 0.0;
    membraneObject.implicitLipid = false;
    membraneObject.TwoD = false;
    membraneObject.isBox = true;
    membraneObject.isSphere = false;
    membraneObject.sphereR = 0.0;
    membraneObject.hasCompartment = false;
    membraneObject.compartmentR = 0.0;
    membraneObject.numberOfFreeLipidsEachState.clear();
}

/*! \brief Fill Parameters with easily recognizable integral values. */
void wr_init_params(Parameters& params)
{
    params.nItr = 1000;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;
    params.numMolTypes = 1;
    params.numTotalSpecies = 3;
    params.numTotalComplex = 2;
    params.numTotalUnits = 4;
    params.numLipids = 5;
    params.timeStep = 1.0;
    params.max2DRxns = 1234;
    params.overlapSepLimit = 0.1;
    params.rMaxLimit = 2.0;
    params.timeWrite = 11;
    params.trajWrite = 22;
    params.restartWrite = 33;
    params.pdbWrite = 44;
    params.assocDissocWrite = true;
    params.checkPoint = 55;
    params.scaleMaxDisplace = 100.0;
    params.transitionWrite = 66;
    params.clusterOverlapCheck = true;
    params.rngwrite = false;
    params.bondedComplexWrite = 77;
}

/*! \brief Run write_restart() on a scenario and return the whole file content.
 *
 * A uniquely named temporary file is used and removed afterwards so nothing is
 * left behind in the working directory.
 */
std::string wr_run_and_read(WrScenario& sc)
{
    const char* fileName = "unit_test_write_restart_tmp.dat";

    std::ofstream restartFile(fileName);
    // Sanity: if the file cannot be opened the rest of the test is meaningless,
    // but we use a non-fatal expectation so the suite keeps going.
    EXPECT_TRUE(restartFile.is_open()) << "Could not open temporary restart file for writing";

    write_restart(sc.simItr, restartFile, sc.params, sc.simulVolume, sc.moleculeList, sc.complexList,
        sc.molTemplateList, sc.forwardRxns, sc.backRxns, sc.createDestructRxns, sc.transmissionRxns,
        sc.observablesList, sc.membraneObject, sc.counterArrays);
    restartFile.close();

    std::ifstream in(fileName);
    std::stringstream ss;
    ss << in.rdbuf();
    in.close();
    std::remove(fileName);

    return ss.str();
}

/*! \brief Split a blob of text into individual lines (newline separated). */
std::vector<std::string> wr_split_lines(const std::string& content)
{
    std::vector<std::string> lines;
    std::string line;
    std::istringstream iss(content);
    while (std::getline(iss, line))
        lines.push_back(line);
    return lines;
}

/*! \brief Index of the first line beginning with `prefix`, or -1 if not found. */
int wr_index_of_line_with_prefix(const std::vector<std::string>& lines, const std::string& prefix)
{
    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].compare(0, prefix.size(), prefix) == 0)
            return static_cast<int>(i);
    }
    return -1;
}

/*! \brief Fetch line `offset` lines after the header identified by `prefix`.
 *
 * Returns an empty string if the header (or the requested offset) is missing.
 */
std::string wr_line_after_header(
    const std::vector<std::string>& lines, const std::string& prefix, size_t offset)
{
    int headerIdx = wr_index_of_line_with_prefix(lines, prefix);
    if (headerIdx < 0)
        return std::string();
    size_t target = static_cast<size_t>(headerIdx) + offset;
    if (target >= lines.size())
        return std::string();
    return lines[target];
}

/*! \brief True when the file contains a line exactly equal to `exact`. */
bool wr_has_exact_line(const std::vector<std::string>& lines, const std::string& exact)
{
    for (const auto& line : lines) {
        if (line == exact)
            return true;
    }
    return false;
}

/*! \brief Build a tiny single-interface MolTemplate named "A". */
MolTemplate wr_make_moltemplate()
{
    MolTemplate oneTemp;
    oneTemp.molTypeIndex = 0;
    oneTemp.molName = "A";
    oneTemp.copies = 1;
    oneTemp.mass = 1.0;
    oneTemp.radius = 1.0;
    oneTemp.comCoord = Coord(0.0, 0.0, 0.0);
    oneTemp.D = Coord(10.0, 10.0, 10.0);
    oneTemp.Dr = Coord(0.1, 0.1, 0.1);

    // One interface named "a" carrying a single state 'U'.
    Interface oneIface;
    oneIface.index = 0;
    oneIface.name = "a";
    oneIface.iCoord = Coord(1.0, 0.0, 0.0);
    Interface::State oneState('U', 0);
    oneIface.stateList.push_back(oneState);
    oneTemp.interfaceList.push_back(oneIface);

    // Auxiliary lists that write_restart() serializes.
    oneTemp.rxnPartners.push_back(0);
    oneTemp.bondList.push_back(std::array<int, 2> { 0, 1 });
    oneTemp.ifacesWithStates.push_back(0);
    oneTemp.monomerList.push_back(0);

    return oneTemp;
}

/*! \brief Build a single unbound molecule that belongs to complex 0. */
Molecule wr_make_molecule()
{
    Molecule oneMol;
    oneMol.index = 0;
    oneMol.isEmpty = false;
    oneMol.myComIndex = 0;
    oneMol.molTypeIndex = 0;
    oneMol.mySubVolIndex = 3;
    oneMol.mass = 1.0;
    oneMol.isLipid = false;
    oneMol.isImplicitLipid = false;
    oneMol.linksToSurface = 0;
    oneMol.isPromoter = false;
    oneMol.comCoord = Coord(1.0, 2.0, 3.0);

    Molecule::Iface iface;
    iface.index = 0;
    iface.relIndex = 0;
    iface.molTypeIndex = 0;
    iface.stateIndex = 0;
    iface.stateIden = 'U';
    iface.isBound = false;
    iface.coord = Coord(2.0, 2.0, 3.0);
    oneMol.interfaceList.push_back(iface);

    oneMol.freelist.push_back(0);
    return oneMol;
}

/*! \brief Build a single complex holding molecule 0. */
Complex wr_make_complex()
{
    Complex oneCom;
    oneCom.index = 0;
    oneCom.isEmpty = false;
    oneCom.radius = 1.0;
    oneCom.mass = 1.0;
    oneCom.linksToSurface = 0;
    oneCom.iLipidIndex = 0;
    oneCom.OnSurface = false;
    oneCom.onFiber = false;
    oneCom.comCoord = Coord(1.0, 2.0, 3.0);
    oneCom.D = Coord(10.0, 10.0, 10.0);
    oneCom.Dr = Coord(0.1, 0.1, 0.1);
    oneCom.memberList.push_back(0);
    oneCom.numEachMol.push_back(1);
    oneCom.lastNumberUpdateItrEachMol.push_back(0);
    return oneCom;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: The "#Parameters" section.
// -----------------------------------------------------------------------------
void test_wr_parameters_section()
{
    std::cerr << "\n[TEST] test_wr_parameters_section\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart()\n"
              << "  Scenario:      minimal state, only Parameters/Membrane populated.\n"
              << "  Pass criteria: the '#Parameters' header is the first line and each\n"
              << "                 'key = value' pair matches the value we supplied.\n";

    WrStaticsGuard guard; // protect global statics
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);
    sc.simItr = 42;
    Parameters::lastUpdateTransition = std::vector<long long int> { 0, 10 };

    std::cerr << "  Calling write_restart() with simItr = " << sc.simItr << "...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    // The file must not be empty and must start with the parameters header.
    EXPECT_FALSE(content.empty()) << "write_restart() produced an empty file";
    ASSERT_FALSE(lines.empty()) << "no lines were written";
    EXPECT_EQ(lines[0].compare(0, 11, "#Parameters"), 0)
        << "first line should be the '#Parameters' header, got: " << lines[0];

    // Integral parameters are written verbatim, so we can compare whole lines.
    std::cerr << "  -> Checking integral parameter lines\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "numItr = 1000")) << "numItr line missing/incorrect";
    EXPECT_TRUE(wr_has_exact_line(lines, "currItr = 42")) << "currItr should echo simItr";
    EXPECT_TRUE(wr_has_exact_line(lines, "numMolTypes = 1")) << "numMolTypes line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "numTotalSpecies = 3")) << "numTotalSpecies line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "numComplexs = 2")) << "numComplexs line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "numTotalUnits = 4")) << "numTotalUnits line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "numLipids = 5")) << "numLipids line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "max2DRxns = 1234")) << "max2DRxns line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "timeWrite = 11")) << "timeWrite line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "trajWrite = 22")) << "trajWrite line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "restartWrite = 33")) << "restartWrite line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "pdbWrite = 44")) << "pdbWrite line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "accocDissocWrite = 1"))
        << "assocDissocWrite should be written as boolean 1";
    EXPECT_TRUE(wr_has_exact_line(lines, "checkPoint = 55")) << "checkPoint line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "transitionWrite = 66")) << "transitionWrite line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "clusterOverlapCheck = 1"))
        << "clusterOverlapCheck should be written as boolean 1";
    EXPECT_TRUE(wr_has_exact_line(lines, "RNGwrite = 0")) << "RNGwrite should be written as boolean 0";
    EXPECT_TRUE(wr_has_exact_line(lines, "bondedComplexWrite = 77")) << "bondedComplexWrite line missing";

    // Floating point values are written with 20 digits, so only check prefixes.
    std::cerr << "  -> Checking floating point parameter lines (prefix match only)\n";
    EXPECT_NE(content.find("currSimTime (s) = "), std::string::npos) << "currSimTime label missing";
    EXPECT_NE(content.find("timestep = 1.0"), std::string::npos) << "timestep value should start with 1.0";
    EXPECT_NE(content.find("ifaceOverlapSepLimit = 0.1"), std::string::npos)
        << "overlapSepLimit value should start with 0.1";
    EXPECT_NE(content.find("rMaxLimit = 2.0"), std::string::npos) << "rMaxLimit value should start with 2.0";
    EXPECT_NE(content.find("scaleMaxDisplace = 100.0"), std::string::npos)
        << "scaleMaxDisplace value should start with 100.0";
    EXPECT_NE(content.find("simulDimensions = 100.0"), std::string::npos)
        << "simulDimensions should echo the water box x dimension";

    // Parameters::lastUpdateTransition is written as "size v0 v1 ...".
    std::cerr << "  -> Checking Parameters::lastUpdateTransition serialization\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "2 0 10"))
        << "lastUpdateTransition should be written as its size followed by its elements";
}

// -----------------------------------------------------------------------------
// Test 2: Membrane / implicit lipid state serialization.
// -----------------------------------------------------------------------------
void test_wr_membrane_and_lipid_states()
{
    std::cerr << "\n[TEST] test_wr_membrane_and_lipid_states\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (membrane portion)\n"
              << "  Scenario:      spherical membrane with 2 implicit lipid states.\n"
              << "  Pass criteria: 'membrane = ...' echoes the membrane counters and\n"
              << "                 'implicitLipidStates = ' lists both state populations.\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);

    // Configure a spherical, implicit-lipid membrane with two lipid states.
    sc.membraneObject.implicitlipidIndex = 7;
    sc.membraneObject.nSites = 3;
    sc.membraneObject.nStates = 2;
    sc.membraneObject.No_free_lipids = 90;
    sc.membraneObject.No_protein = 11;
    sc.membraneObject.totalSA = 500.0;
    sc.membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 40, 50 };
    sc.membraneObject.implicitLipid = true;
    sc.membraneObject.isBox = false;
    sc.membraneObject.isSphere = true;
    sc.membraneObject.sphereR = 25.0;
    sc.membraneObject.hasCompartment = true;
    sc.membraneObject.compartmentR = 5.0;

    std::cerr << "  Calling write_restart()...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    // "membrane = implicitlipidIndex nSites nStates No_free_lipids No_protein totalSA"
    std::cerr << "  -> Checking the 'membrane = ' summary line\n";
    EXPECT_NE(content.find("membrane = 7 3 2 90 11 500"), std::string::npos)
        << "membrane summary line does not contain the expected counters";

    // The state populations are space separated on one line.
    std::cerr << "  -> Checking 'implicitLipidStates = 40 50'\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "implicitLipidStates = 40 50"))
        << "per-state free lipid numbers were not serialized correctly";

    // "implicitLipidsParams = implicitLipid TwoD isBox isSphere sphereR hasCompartment compartmentR"
    std::cerr << "  -> Checking 'implicitLipidsParams = ' booleans and radii\n";
    EXPECT_NE(content.find("implicitLipidsParams = 1 0 0 1 25"), std::string::npos)
        << "implicit lipid flags / sphere radius were not serialized correctly";
}

// -----------------------------------------------------------------------------
// Test 3: The "#MolTemplates" section for one simple template.
// -----------------------------------------------------------------------------
void test_wr_moltemplate_section()
{
    std::cerr << "\n[TEST] test_wr_moltemplate_section\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (MolTemplate portion)\n"
              << "  Scenario:      one MolTemplate 'A' with a single interface 'a~U'.\n"
              << "  Pass criteria: the three static-member lines directly follow the\n"
              << "                 '#MolTemplates' header, and the template's index,\n"
              << "                 name, flags, interface and state are all present.\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);
    sc.molTemplateList.push_back(wr_make_moltemplate());

    // Static template bookkeeping that write_restart() serializes.
    MolTemplate::numMolTypes = 1;
    MolTemplate::numEachMolType = std::vector<int> { 1 };
    MolTemplate::absToRelIface = std::vector<int> { 0 };
    Interface::State::totalNumOfStates = 1;

    std::cerr << "  Calling write_restart()...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    // The section header must exist.
    const int headerIdx = wr_index_of_line_with_prefix(lines, "#MolTemplates");
    EXPECT_GE(headerIdx, 0) << "'#MolTemplates' header was not written";

    // Line 1 after header: numMolTypes then numEachMolType contents.
    std::cerr << "  -> Line 1 after header should be 'numMolTypes numEachMolType...'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#MolTemplates", 1), "1 1")
        << "static numMolTypes/numEachMolType line incorrect";

    // Line 2 after header: absToRelIface size then contents.
    std::cerr << "  -> Line 2 after header should be 'absToRelIface.size() contents'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#MolTemplates", 2), "1 0")
        << "static absToRelIface line incorrect";

    // Line 3 after header: total number of interface states.
    std::cerr << "  -> Line 3 after header should be Interface::State::totalNumOfStates\n";
    EXPECT_EQ(wr_line_after_header(lines, "#MolTemplates", 3), "1")
        << "totalNumOfStates line incorrect";

    // Line 4 after header: "molTypeIndex molName".
    std::cerr << "  -> Line 4 after header should be 'molTypeIndex molName'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#MolTemplates", 4), "0 A")
        << "template index/name line incorrect";

    // Flag line: isLipid isImplicitLipid isRod isPoint checkOverlap countTransition
    //            transitionMatrixSize outsideCompartment insideCompartment
    //            crossesCompartment transmissionRxnIndex
    std::cerr << "  -> Checking the MolTemplate boolean-flag line\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 0 0 0 0 0 500 0 0 0 -1"))
        << "MolTemplate flag line (with default transitionMatrixSize 500) not found";

    // Interface and state descriptors.
    std::cerr << "  -> Checking interface ('0 a') and state ('0 U') descriptor lines\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 a")) << "interface index/name line not found";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 U")) << "interface state index/identity line not found";

    // The optional bond (0,1) supplied in the template.
    std::cerr << "  -> Checking the optional bond line '0 1'\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 1")) << "bondList entry was not serialized";
}

// -----------------------------------------------------------------------------
// Test 4: Transition matrix / lifetime output when countTransition == true.
// -----------------------------------------------------------------------------
void test_wr_moltemplate_transition_matrix()
{
    std::cerr << "\n[TEST] test_wr_moltemplate_transition_matrix\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (transition-matrix branch)\n"
              << "  Scenario:      MolTemplate with countTransition = true and a 2x2\n"
              << "                 transition matrix / lifetime table.\n"
              << "  Pass criteria: each transition matrix row and each lifetime row is\n"
              << "                 written out (this branch is skipped otherwise).\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);

    MolTemplate oneTemp = wr_make_moltemplate();
    oneTemp.countTransition = true;
    oneTemp.transitionMatrixSize = 2; // keep the matrix tiny for the test
    oneTemp.transitionMatrix = std::vector<std::vector<long long int>> { { 1, 2 }, { 3, 4 } };
    oneTemp.lifeTime = std::vector<std::vector<double>> { { 0.5, 1.5 }, { 2.5, 3.5 } };
    sc.molTemplateList.push_back(oneTemp);

    MolTemplate::numMolTypes = 1;
    MolTemplate::numEachMolType = std::vector<int> { 1 };
    MolTemplate::absToRelIface = std::vector<int> { 0 };
    Interface::State::totalNumOfStates = 1;

    std::cerr << "  Calling write_restart()...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    // Each transition matrix row is written as " v0 v1" (leading space per value).
    std::cerr << "  -> Checking transition matrix rows ' 1 2' and ' 3 4'\n";
    EXPECT_TRUE(wr_has_exact_line(lines, " 1 2")) << "first transition matrix row missing";
    EXPECT_TRUE(wr_has_exact_line(lines, " 3 4")) << "second transition matrix row missing";

    // Each lifetime row is written as "count v0 v1"; doubles use 20 digits so we
    // only require the row to begin with the count and the first value.
    std::cerr << "  -> Checking lifetime rows begin with '2 0.5' and '2 2.5'\n";
    EXPECT_NE(content.find("2 0.5"), std::string::npos) << "first lifeTime row missing";
    EXPECT_NE(content.find("2 2.5"), std::string::npos) << "second lifeTime row missing";

    // The flag line should now show countTransition = 1 and the small matrix size.
    std::cerr << "  -> Checking the flag line reports countTransition=1, size=2\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 0 0 0 0 1 2 0 0 0 -1"))
        << "MolTemplate flag line should report countTransition=1 and transitionMatrixSize=2";
}

// -----------------------------------------------------------------------------
// Test 5: The "#Reactions" section (forward, back, create/destruct, transmission).
// -----------------------------------------------------------------------------
void test_wr_reactions_section()
{
    std::cerr << "\n[TEST] test_wr_reactions_section\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (reaction portion)\n"
              << "  Scenario:      one of each reaction flavour (forward, back,\n"
              << "                 create/destruct, transmission).\n"
              << "  Pass criteria: the counts line right after '#Reactions' matches the\n"
              << "                 container sizes, and each reaction's unique observe\n"
              << "                 label / reactant descriptors appear in the output.\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);

    RxnBase::numberOfRxns = 4;
    RxnBase::totRxnSpecies = 3;

    // ---- one forward (bimolecular) reaction --------------------------------
    ForwardRxn fwd;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.rxnLabel = "bindAB";
    fwd.rxnType = ReactionType::bimolecular;
    fwd.isSymmetric = false;
    fwd.isOnMem = false;
    fwd.hasStateChange = false;
    fwd.isObserved = true;
    fwd.observeLabel = "fwdObs";
    fwd.productName = "AB";
    fwd.isReversible = true;
    fwd.conjBackRxnIndex = 0;
    fwd.irrevRingClosure = false;
    fwd.bindRadSameCom = 1.1;
    fwd.loopCoopFactor = 1.0;
    fwd.length3Dto2D = 2.0;
    fwd.area3Dto1D = 4.0;
    fwd.bindRadius = 1.0;
    // Use finite angles so the output stays human readable (defaults are NaN).
    fwd.assocAngles = ForwardRxn::Angles(M_PI, M_PI, 0.0, 0.0, 0.0);
    fwd.norm1 = Vector(0.0, 0.0, 1.0);
    fwd.norm2 = Vector(0.0, 0.0, 1.0);
    fwd.excludeVolumeBound = false;
    fwd.isCoupled = false;
    fwd.intReactantList = std::vector<int> { 0, 1 };
    fwd.intProductList = std::vector<int> { 2 };
    fwd.reactantListNew.push_back(RxnIface("a", 0, 0, 0, 'U', false));
    fwd.reactantListNew.push_back(RxnIface("b", 1, 1, 0, 'U', false));
    fwd.productListNew.push_back(RxnIface("a", 0, 2, 0, 'U', true));
    {
        RxnBase::RateState oneRate;
        oneRate.rate = 10.0;
        // One ancillary-interface requirement, to cover that nested loop.
        oneRate.otherIfaceLists = std::vector<std::vector<RxnIface>> {
            std::vector<RxnIface> { RxnIface("c", 0, 3, 1, 'U', true) }
        };
        fwd.rateList.push_back(oneRate);
    }
    sc.forwardRxns.push_back(fwd);

    // ---- the conjugate back reaction --------------------------------------
    BackRxn bck;
    bck.absRxnIndex = 1;
    bck.relRxnIndex = 0;
    bck.rxnType = ReactionType::bimolecular;
    bck.isObserved = false;
    bck.observeLabel = "backObs";
    bck.conjForwardRxnIndex = 0;
    bck.isCoupled = false;
    bck.intReactantList = std::vector<int> { 2 };
    bck.intProductList = std::vector<int> { 0, 1 };
    bck.reactantListNew.push_back(RxnIface("a", 0, 2, 0, 'U', true));
    bck.productListNew.push_back(RxnIface("a", 0, 0, 0, 'U', false));
    {
        RxnBase::RateState oneRate;
        oneRate.rate = 1.0;
        bck.rateList.push_back(oneRate);
    }
    sc.backRxns.push_back(bck);

    // ---- one destruction reaction -----------------------------------------
    CreateDestructRxn cd;
    cd.absRxnIndex = 2;
    cd.relRxnIndex = 0;
    cd.rxnType = ReactionType::destruction;
    cd.isOnMem = false;
    cd.isObserved = false;
    cd.observeLabel = "cdObs";
    cd.creationRadius = 1.0;
    cd.intReactantList = std::vector<int> { 0 };
    {
        CreateDestructRxn::CreateDestructMol oneMol(
            0, std::vector<RxnIface> { RxnIface("a", 0, 0, 0, 'U', false) });
        oneMol.molName = "Adestroy";
        cd.reactantMolList.push_back(oneMol);
    }
    {
        RxnBase::RateState oneRate;
        oneRate.rate = 0.5;
        cd.rateList.push_back(oneRate);
    }
    sc.createDestructRxns.push_back(cd);

    // ---- one transmission reaction ----------------------------------------
    TransmissionRxn tr;
    tr.absRxnIndex = 3;
    tr.relRxnIndex = 0;
    tr.rxnType = ReactionType::transmission;
    tr.isOnMem = false;
    tr.isObserved = false;
    tr.observeLabel = "trObs";
    tr.intReactantList = std::vector<int> { 0 };
    {
        TransmissionRxn::TransmissionMol oneMol(
            0, std::vector<RxnIface> { RxnIface("a", 0, 0, 0, 'U', false) });
        oneMol.molName = "Atransmit";
        tr.reactantMolList.push_back(oneMol);
    }
    {
        RxnBase::RateState oneRate;
        oneRate.rate = 0.25;
        tr.rateList.push_back(oneRate);
    }
    sc.transmissionRxns.push_back(tr);

    std::cerr << "  Calling write_restart() with 1 forward / 1 back / 1 createDestruct / 1 transmission...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    // Header + counts line.
    EXPECT_GE(wr_index_of_line_with_prefix(lines, "#Reactions"), 0) << "'#Reactions' header missing";
    std::cerr << "  -> Line right after '#Reactions' should be "
              << "'numberOfRxns nFwd nBack nCreateDestruct nTransmission totRxnSpecies'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#Reactions", 1), "4 1 1 1 1 3")
        << "reaction count line does not match the container sizes";

    // Forward reaction identifiers.
    std::cerr << "  -> Checking forward reaction identifier and observable lines\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 0 bindAB"))
        << "forward reaction index/label line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "1 fwdObs AB"))
        << "forward reaction isObserved/observeLabel/productName line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "1 0 0 0"))
        << "forward reaction type/flags line (bimolecular=1) missing";

    // Ancillary interface of the forward reaction's rate.
    std::cerr << "  -> Checking the ancillary interface line of the rate state\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 c 3 1 U 1"))
        << "ancillary (otherIfaceLists) interface line missing";

    // Back reaction identifiers.
    std::cerr << "  -> Checking back reaction observable line\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 backObs")) << "back reaction observeLabel line missing";

    // Create/destruct and transmission reaction identifiers.
    std::cerr << "  -> Checking create/destruct and transmission reaction lines\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 cdObs")) << "createDestruct observeLabel line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 Adestroy 1"))
        << "createDestruct reactant molecule descriptor missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 trObs")) << "transmission observeLabel line missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 Atransmit 1"))
        << "transmission reactant molecule descriptor missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "4 0"))
        << "destruction reaction type line (destruction=4, isOnMem=0) missing";
    EXPECT_TRUE(wr_has_exact_line(lines, "7 0"))
        << "transmission reaction type line (transmission=7, isOnMem=0) missing";
}

// -----------------------------------------------------------------------------
// Test 6: The "#All Molecules and coordinates" section.
// -----------------------------------------------------------------------------
void test_wr_molecules_section()
{
    std::cerr << "\n[TEST] test_wr_molecules_section\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (Molecule portion)\n"
              << "  Scenario:      one molecule, Molecule::numberOfMolecules = 1,\n"
              << "                 Molecule::emptyMolList = {5}.\n"
              << "  Pass criteria: counts line, molecule descriptor line, coordinate\n"
              << "                 line and the trailing emptyMolList line all match.\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);
    sc.moleculeList.push_back(wr_make_molecule());
    Molecule::numberOfMolecules = 1;
    Molecule::emptyMolList = std::vector<int> { 5 };

    std::cerr << "  Calling write_restart()...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    EXPECT_GE(wr_index_of_line_with_prefix(lines, "#All Molecules"), 0)
        << "'#All Molecules and coordinates' header missing";

    // "moleculeList.size() Molecule::numberOfMolecules"
    std::cerr << "  -> Line after header should be 'listSize numberOfMolecules'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#All Molecules", 1), "1 1")
        << "molecule count line incorrect";

    // "index isEmpty myComIndex molTypeIndex mySubVolIndex"
    std::cerr << "  -> Checking molecule descriptor line '0 0 0 0 3'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#All Molecules", 2), "0 0 0 0 3")
        << "molecule index/complex/type/subvolume line incorrect";

    // Center-of-mass coordinates are printed with 20 fixed digits: prefix check.
    std::cerr << "  -> Checking the center-of-mass coordinate line starts with '1.0'\n";
    const std::string comLine = wr_line_after_header(lines, "#All Molecules", 4);
    EXPECT_EQ(comLine.compare(0, 3, "1.0"), 0)
        << "molecule COM coordinate line should start with x = 1.0, got: " << comLine;

    // The interface descriptor: "index relIndex molTypeIndex stateIndex stateIden isBound"
    std::cerr << "  -> Checking the interface descriptor line '0 0 0 0 U 0'\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "0 0 0 0 U 0"))
        << "molecule interface descriptor line missing";

    // Molecule::emptyMolList is written as "size contents".
    std::cerr << "  -> Checking the emptyMolList line '1 5'\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "1 5")) << "Molecule::emptyMolList line missing";
}

// -----------------------------------------------------------------------------
// Test 7: The "#All Complexes and their components" section.
// -----------------------------------------------------------------------------
void test_wr_complexes_section()
{
    std::cerr << "\n[TEST] test_wr_complexes_section\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (Complex portion)\n"
              << "  Scenario:      one complex containing molecule 0,\n"
              << "                 Complex::numberOfComplexes = 1, emptyComList = {7}.\n"
              << "  Pass criteria: counts line, surface/link line, membership lines and\n"
              << "                 the trailing emptyComList line all match.\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);
    sc.moleculeList.push_back(wr_make_molecule());
    sc.complexList.push_back(wr_make_complex());
    Molecule::numberOfMolecules = 1;
    Complex::numberOfComplexes = 1;
    Complex::emptyComList = std::vector<int> { 7 };

    std::cerr << "  Calling write_restart()...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    EXPECT_GE(wr_index_of_line_with_prefix(lines, "#All Complexes"), 0)
        << "'#All Complexes and their components' header missing";

    // "complexList.size() Complex::numberOfComplexes"
    std::cerr << "  -> Line after header should be 'listSize numberOfComplexes'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#All Complexes", 1), "1 1")
        << "complex count line incorrect";

    // "index isEmpty radius mass" (radius/mass are fixed with 20 digits).
    std::cerr << "  -> Checking the complex descriptor line starts with '0 0 1.0'\n";
    const std::string descLine = wr_line_after_header(lines, "#All Complexes", 2);
    EXPECT_EQ(descLine.compare(0, 6, "0 0 1."), 0)
        << "complex index/isEmpty/radius line unexpected, got: " << descLine;

    // "linksToSurface iLipidIndex OnSurface"
    std::cerr << "  -> Checking the links/lipid/OnSurface line '0 0 0'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#All Complexes", 3), "0 0 0")
        << "complex linksToSurface/iLipidIndex/OnSurface line incorrect";

    // onFiber on its own line.
    std::cerr << "  -> Checking the onFiber line '0'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#All Complexes", 4), "0")
        << "complex onFiber line incorrect";

    // Complex::emptyComList is written as "size contents".
    std::cerr << "  -> Checking the emptyComList line '1 7'\n";
    EXPECT_TRUE(wr_has_exact_line(lines, "1 7")) << "Complex::emptyComList line missing";
}

// -----------------------------------------------------------------------------
// Test 8: The "#Observables" section.
// -----------------------------------------------------------------------------
void test_wr_observables_section()
{
    std::cerr << "\n[TEST] test_wr_observables_section\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (observables portion)\n"
              << "  Scenario:      two observables, 'dimer' = 3 and 'trimer' = 4.\n"
              << "  Pass criteria: the observable count is followed by one\n"
              << "                 'name count' line per map entry (map order).\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);
    sc.observablesList["dimer"] = 3;
    sc.observablesList["trimer"] = 4;

    std::cerr << "  Calling write_restart()...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    EXPECT_GE(wr_index_of_line_with_prefix(lines, "#Observables"), 0)
        << "'#Observables' header missing";

    std::cerr << "  -> Line after header should be the observable count\n";
    EXPECT_EQ(wr_line_after_header(lines, "#Observables", 1), "2")
        << "observable count line incorrect";

    // std::map iterates in sorted key order, so 'dimer' precedes 'trimer'.
    std::cerr << "  -> Checking the two 'name count' lines in map (sorted) order\n";
    EXPECT_EQ(wr_line_after_header(lines, "#Observables", 2), "dimer 3")
        << "first observable line incorrect";
    EXPECT_EQ(wr_line_after_header(lines, "#Observables", 3), "trimer 4")
        << "second observable line incorrect";
}

// -----------------------------------------------------------------------------
// Test 9: The two counterArrays sections.
// -----------------------------------------------------------------------------
void test_wr_counter_arrays_section()
{
    std::cerr << "\n[TEST] test_wr_counter_arrays_section\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (copyCounters portion)\n"
              << "  Scenario:      distinct values for every cancel counter, small\n"
              << "                 event histograms and a two-entry bindPairList.\n"
              << "  Pass criteria: the counter summary line, the three event-histogram\n"
              << "                 lines and the bindPairList block all match.\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);

    // Distinct values make the single summary line unambiguous.
    sc.counterArrays.nLoops = 1;
    sc.counterArrays.nCancelOverlapPartner = 2;
    sc.counterArrays.nCancelOverlapSystem = 3;
    sc.counterArrays.nCancelDisplace2D = 4;
    sc.counterArrays.nCancelDisplace3D = 5;
    sc.counterArrays.nCancelDisplace3Dto2D = 6;
    sc.counterArrays.nCancelSpanBox = 7;
    sc.counterArrays.nAssocSuccess = 8;
    sc.counterArrays.eventArraySize = 9;
    sc.counterArrays.events3D = std::vector<int> { 10, 11 };
    sc.counterArrays.events3Dto2D = std::vector<int> { 12 };
    sc.counterArrays.events2D = std::vector<int> { 13, 14 };
    sc.counterArrays.bindPairList = std::vector<std::vector<int>> { { 1, 2 }, { 3 } };

    std::cerr << "  Calling write_restart()...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    EXPECT_GE(wr_index_of_line_with_prefix(lines, "#counterArrays.NLoops"), 0)
        << "'#counterArrays.NLoops .nCancels' header missing";

    // Summary line of all cancel/success counters plus the event array size.
    std::cerr << "  -> Checking the counter summary line '1 2 3 4 5 6 7 8 9'\n";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.NLoops", 1), "1 2 3 4 5 6 7 8 9")
        << "counterArrays summary line incorrect";

    // Each event histogram is written as "v0 v1 ... " -- note the trailing space
    // produced by "restartFile << event << ' '".
    std::cerr << "  -> Checking the three event-histogram lines (trailing space expected)\n";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.NLoops", 2), "10 11 ")
        << "events3D histogram line incorrect";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.NLoops", 3), "12 ")
        << "events3Dto2D histogram line incorrect";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.NLoops", 4), "13 14 ")
        << "events2D histogram line incorrect";

    // bindPairList: outer size, then for each entry its size and its elements
    // (each element prefixed with a space).
    EXPECT_GE(wr_index_of_line_with_prefix(lines, "#counterArrays.bindPairList"), 0)
        << "'#counterArrays.bindPairList' header missing";
    std::cerr << "  -> Checking the bindPairList block (outer size, then size/elements)\n";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.bindPairList", 1), "2")
        << "bindPairList outer size line incorrect";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.bindPairList", 2), "2")
        << "first bindPair size line incorrect";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.bindPairList", 3), " 1 2")
        << "first bindPair contents line incorrect";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.bindPairList", 4), "1")
        << "second bindPair size line incorrect";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.bindPairList", 5), " 3")
        << "second bindPair contents line incorrect";
}

// -----------------------------------------------------------------------------
// Test 10: Smoke test -- all section headers are written, in the documented
//          order, even when every container is empty.
// -----------------------------------------------------------------------------
void test_wr_all_headers_in_order_for_empty_system()
{
    std::cerr << "\n[TEST] test_wr_all_headers_in_order_for_empty_system\n"
              << "  Source file:   src/io/write_restart.cpp\n"
              << "  Function:      write_restart() (whole-file structure)\n"
              << "  Scenario:      completely empty system (no molecules, complexes,\n"
              << "                 templates, reactions or observables).\n"
              << "  Pass criteria: every section header is present and the header line\n"
              << "                 indices are strictly increasing (correct order).\n";

    WrStaticsGuard guard;
    WrScenario sc;
    wr_init_params(sc.params);
    wr_init_membrane(sc.membraneObject);

    // Explicitly zero every static so the empty system is self-consistent.
    MolTemplate::numMolTypes = 0;
    MolTemplate::numEachMolType.clear();
    MolTemplate::absToRelIface.clear();
    Interface::State::totalNumOfStates = 0;
    RxnBase::numberOfRxns = 0;
    RxnBase::totRxnSpecies = 0;
    Molecule::numberOfMolecules = 0;
    Molecule::emptyMolList.clear();
    Complex::numberOfComplexes = 0;
    Complex::emptyComList.clear();
    Parameters::lastUpdateTransition.clear();

    std::cerr << "  Calling write_restart() on an empty system...\n";
    const std::string content = wr_run_and_read(sc);
    const std::vector<std::string> lines = wr_split_lines(content);

    // Collect the line index of each expected header.
    const std::vector<std::string> expectedHeaders {
        "#Parameters",
        "#MolTemplates",
        "#Reactions",
        "#All Molecules",
        "#All Complexes",
        "#Observables",
        "#counterArrays.NLoops",
        "#counterArrays.bindPairList",
    };

    int previousIdx = -1;
    for (const auto& header : expectedHeaders) {
        const int idx = wr_index_of_line_with_prefix(lines, header);
        std::cerr << "  -> Header \"" << header << "\" found at line " << idx << '\n';
        EXPECT_GE(idx, 0) << "expected section header not written: " << header;
        // Header order must be strictly increasing.
        EXPECT_GT(idx, previousIdx) << "section header out of order: " << header;
        if (idx > previousIdx)
            previousIdx = idx;
    }

    // Even with nothing in the system, the container-size lines must be "0"-ish.
    std::cerr << "  -> Checking that empty containers serialize zero sizes\n";
    EXPECT_EQ(wr_line_after_header(lines, "#Reactions", 1), "0 0 0 0 0 0")
        << "empty reaction section should report all zero counts";
    EXPECT_EQ(wr_line_after_header(lines, "#All Molecules", 1), "0 0")
        << "empty molecule section should report zero molecules";
    EXPECT_EQ(wr_line_after_header(lines, "#All Complexes", 1), "0 0")
        << "empty complex section should report zero complexes";
    EXPECT_EQ(wr_line_after_header(lines, "#Observables", 1), "0")
        << "empty observables section should report zero observables";
    EXPECT_EQ(wr_line_after_header(lines, "#counterArrays.bindPairList", 1), "0")
        << "empty bindPairList should report size zero";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: one TEST per named test_* helper so that failures are
// reported individually while every test still runs.
// -----------------------------------------------------------------------------
TEST(WriteRestart, ParametersSection) { test_wr_parameters_section(); }
TEST(WriteRestart, MembraneAndLipidStates) { test_wr_membrane_and_lipid_states(); }
TEST(WriteRestart, MolTemplateSection) { test_wr_moltemplate_section(); }
TEST(WriteRestart, MolTemplateTransitionMatrix) { test_wr_moltemplate_transition_matrix(); }
TEST(WriteRestart, ReactionsSection) { test_wr_reactions_section(); }
TEST(WriteRestart, MoleculesSection) { test_wr_molecules_section(); }
TEST(WriteRestart, ComplexesSection) { test_wr_complexes_section(); }
TEST(WriteRestart, ObservablesSection) { test_wr_observables_section(); }
TEST(WriteRestart, CounterArraysSection) { test_wr_counter_arrays_section(); }
TEST(WriteRestart, AllHeadersInOrderForEmptySystem) { test_wr_all_headers_in_order_for_empty_system(); }