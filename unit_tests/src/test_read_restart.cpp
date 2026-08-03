/*! \file test_read_restart.cpp
 *
 * ### Unit test for src/io/read_restart.cpp
 *
 * The single function under test is:
 *
 *     void read_restart(long long int& simItr, std::ifstream& restartFile,
 *                       Parameters& params, SimulVolume& simulVolume,
 *                       std::vector<Molecule>& moleculeList,
 *                       std::vector<Complex>& complexList,
 *                       std::vector<MolTemplate>& molTemplateList,
 *                       std::vector<ForwardRxn>& forwardRxns,
 *                       std::vector<BackRxn>& backRxns,
 *                       std::vector<CreateDestructRxn>& createDestructRxns,
 *                       std::vector<TransmissionRxn>& transmissionRxns,
 *                       std::map<std::string,int>& observablesList,
 *                       Membrane& membraneObject, copyCounters& counterArrays)
 *
 * read_restart() is a pure *parser*: it consumes a plain-text NERDSS restart
 * file and fills in every simulation container.  It has no return value, so the
 * only way to test it is to
 *
 *   1. synthesise a restart file whose byte layout matches exactly what the
 *      parser expects (this is done once by rr_write_test_restart_file()),
 *   2. run the parser on that file, and
 *   3. assert that every container/field received the value that was written.
 *
 * Because the parser reads the whole file in one pass, every test below has to
 * re-read the complete file; each test then verifies one logical *section* of
 * the restart format:
 *
 *   - simulation Parameters (nItr, timeStep, write intervals, ...)
 *   - Membrane / water box description
 *   - MolTemplates (including interfaces and interface states)
 *   - Reactions (one ForwardRxn + its BackRxn)
 *   - Molecules (including a bound interface / Interaction)
 *   - Complexes
 *   - Observables and the copyCounters block
 *
 * Everything is checked with EXPECT_* so a mis-parse in one section does not
 * abort the remaining tests.
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

namespace {

//! Name of the scratch restart file created/destroyed by these tests.
const char* kRRRestartFileName = "test_read_restart_tmp.restart";

//! Size of the association-event histograms stored in the restart file.
const int kRREventArraySize = 20;

// -----------------------------------------------------------------------------
// Helper: write a synthetic - but format-exact - restart file.
//
// Every line below corresponds 1:1 to a read in read_restart().  Comments give
// the field(s) each line supplies.  The parser skips a "header" line before
// each major section (it calls ignore(...,'\n')), so the '#...' lines below are
// deliberately present as throw-away separators.
// -----------------------------------------------------------------------------
void rr_write_test_restart_file(const std::string& path)
{
    std::ofstream out(path);

    /* ------------------------- PARAMETERS SECTION ------------------------- */
    out << "#Parameters\n"; // skipped by the very first ignore()
    out << "nItr = 1000\n"; // params.nItr
    out << "simItr = 250\n"; // simItr (and params.itrRestartFrom)
    out << "timeRestartFrom = 0.0025\n"; // params.timeRestartFrom
    out << "numMolTypes = 1\n"; // params.numMolTypes
    out << "numTotalSpecies = 3\n"; // params.numTotalSpecies
    out << "numTotalComplex = 1\n"; // params.numTotalComplex
    out << "numTotalUnits = 4\n"; // params.numTotalUnits
    out << "numLipids = 0\n"; // params.numLipids
    out << "timeStep = 0.1\n"; // params.timeStep
    out << "max2DRxns = 100\n"; // params.max2DRxns
    out << "waterBox = 20 30 40\n"; // membrane waterBox x,y,z
    // implicitlipidIndex nSites nStates No_free_lipids No_protein totalSA
    out << "implicitLipid = -1 0 1 12 3 1234.5\n";
    // one entry per membrane state (nStates == 1 above)
    out << "freeLipidsEachState = 12\n";
    // implicitLipid TwoD isBox isSphere sphereR hasCompartment compartmentR
    out << "membraneFlags = 0 0 1 0 0 0 0\n";
    out << "overlapSepLimit = 0.2\n"; // params.overlapSepLimit
    out << "rMaxLimit = 5.5\n"; // params.rMaxLimit
    out << "timeWrite = 10\n"; // params.timeWrite
    out << "trajWrite = 20\n"; // params.trajWrite
    out << "restartWrite = 30\n"; // params.restartWrite
    out << "pdbWrite = 40\n"; // params.pdbWrite
    out << "assocDissocWrite = 1\n"; // params.assocDissocWrite
    out << "checkPoint = 50\n"; // params.checkPoint
    out << "scaleMaxDisplace = 60.5\n"; // params.scaleMaxDisplace
    out << "transitionWrite = 70\n"; // params.transitionWrite
    out << "clusterOverlapCheck = 1\n"; // params.clusterOverlapCheck
    // optional '='-containing lines consumed by the getline() loop
    out << "RNGwrite = 0\n"; // params.rngwrite
    out << "bondedComplexWrite = -1\n"; // params.bondedComplexWrite
    out << "0\n"; // lastUpdateTransition size (== 0 -> empty)

    /* ------------------------ MOL TEMPLATE SECTION ------------------------ */
    out << "#MolTemplates\n"; // skipped
    out << "1 1\n"; // numMolTypes, numEachMolType[...]
    out << "1 0\n"; // absToRelIface size, entries
    out << "1\n"; // Interface::State::totalNumOfStates
    out << "0 testmol\n"; // molTypeIndex, molName
    out << "5 1000 2.5\n"; // copies, mass, radius
    // isLipid isImplicitLipid isRod isPoint checkOverlap countTransition
    // transitionMatrixSize outside inside crosses transmissionRxnIndex
    out << "0 0 0 0 0 0 500 0 0 0 -1\n";
    out << "0 0 0\n"; // template COM
    out << "25 25 25\n"; // translational D
    out << "0.5 0.5 0.5\n"; // rotational Dr
    out << "1 0\n"; // rxnPartners size + entries
    out << "0\n"; // bondList size (none)
    out << "1\n"; // number of interfaces on this template
    out << "0 iface0\n"; // interface index, name
    out << "1 0 0\n"; // interface internal coordinate
    out << "1\n"; // number of states on this interface
    out << "0 U\n"; // state index, state identity char
    out << "1 1\n"; // rxnPartners size + entries
    out << "1 0\n"; // myForwardRxns size + entries
    out << "0\n"; // myCreateDestructRxns size
    out << "0\n"; // stateChangeRxns size
    out << "0\n"; // ifacesWithStates size
    out << "2 3 4\n"; // monomerList size + entries

    /* -------------------------- REACTION SECTION -------------------------- */
    out << "#Reactions\n"; // skipped
    // numberOfRxns, #forward, #back, #createDestruct, #transmission, totRxnSpecies
    out << "2 1 1 0 0 4\n";
    /* --- one ForwardRxn --- */
    out << "0 0 label0\n"; // absRxnIndex, relRxnIndex, rxnLabel
    out << "1 0 0 0\n"; // rxnType(=bimolecular), isSymmetric, isOnMem, hasStateChange
    out << "0 A_B\n"; // isObserved, productName
    // isReversible conjBackRxnIndex irrevRingClosure bindRadSameCom
    // loopCoopFactor length3Dto2D area3Dto1D
    out << "1 0 0 1.1 1 2 4\n";
    // bindRadius theta1 theta2 phi1 phi2 omega ("nan" -> quiet_NaN)
    out << "1.5 1.5708 1.5708 nan nan 3.14159\n";
    out << "0 0 1\n"; // norm1
    out << "0 0 1\n"; // norm2
    out << "0\n"; // excludeVolumeBound
    out << "0\n"; // isCoupled
    out << "2 0 1\n"; // intReactantList size + entries
    out << "1 2\n"; // intProductList size + entries
    out << "2\n"; // reactantListNew size
    out << "0\n" << "iface0 0 0\n" << "U 0\n"; // reactant 1
    out << "0\n" << "iface0 1 0\n" << "U 0\n"; // reactant 2
    out << "1\n"; // productListNew size
    out << "0\n" << "iface0 2 0\n" << "U 1\n"; // product 1
    out << "1\n"; // rateList size
    out << "10.5\n"; // rate
    out << "0\n"; // otherIfaceLists size
    /* --- the conjugate BackRxn --- */
    out << "1 0\n"; // absRxnIndex, relRxnIndex
    out << "1 0 0 0\n"; // rxnType, isSymmetric, isOnMem, hasStateChange
    out << "0\n"; // isObserved
    out << "0\n"; // conjForwardRxnIndex
    out << "0\n"; // isCoupled
    out << "1 2\n"; // intReactantList
    out << "2 0 1\n"; // intProductList
    out << "1\n"; // reactantListNew size
    out << "0\n" << "iface0 2 0\n" << "U 1\n";
    out << "2\n"; // productListNew size
    out << "0\n" << "iface0 0 0\n" << "U 0\n";
    out << "0\n" << "iface0 1 0\n" << "U 0\n";
    out << "1\n"; // rateList size
    out << "0.5\n"; // rate
    out << "0\n"; // otherIfaceLists size

    /* -------------------------- MOLECULE SECTION -------------------------- */
    out << "#Molecules\n"; // skipped
    out << "2 2\n"; // moleculeList size, Molecule::numberOfMolecules
    /* --- molecule 0: free interface --- */
    out << "0 0 0 0 5\n"; // index isEmpty myComIndex molTypeIndex mySubVolIndex
    out << "1000 0 0 0 0 0\n"; // mass isLipid isImplicitLipid links isPromoter isEmpty
    out << "1.0 2.0 3.0\n"; // COM
    out << "1 0\n"; // freelist
    out << "0\n"; // bndlist
    out << "0\n"; // bndpartner
    out << "1\n"; // interfaceList size
    out << "0 0 0 0 U 0\n"; // index relIndex molTypeIndex stateIndex iden isBound
    out << "2.0 2.0 3.0\n"; // interface coordinate
    out << "0\n0\n0\n0\n0\n0\n"; // six (empty) reweighting lists
    /* --- molecule 1: bound interface (exercises the Interaction branch) --- */
    out << "1 0 0 0 6\n";
    out << "1000 0 0 1 0 0\n";
    out << "-1.0 -2.0 -3.0\n";
    out << "0\n"; // freelist (empty)
    out << "1 0\n"; // bndlist
    out << "1 0\n"; // bndpartner
    out << "1\n"; // interfaceList size
    out << "0 0 0 0 U 1\n"; // isBound == 1
    out << "-2.0 -2.0 -3.0\n";
    out << "0 0 0\n"; // partnerIndex partnerIfaceIndex conjBackRxn
    out << "0\n0\n0\n0\n0\n0\n"; // reweighting lists
    out << "0\n"; // Molecule::emptyMolList size

    /* --------------------------- COMPLEX SECTION -------------------------- */
    out << "#Complexes\n"; // skipped
    out << "1 1\n"; // complexList size, Complex::numberOfComplexes
    out << "0 0 3.5 2000\n"; // index isEmpty radius mass
    out << "0 -1 0\n"; // linksToSurface iLipidIndex OnSurface
    out << "0\n"; // onFiber
    out << "0.5 0.5 0.5\n"; // COM
    out << "25 25 25\n"; // D
    out << "0.5 0.5 0.5\n"; // Dr
    out << "2 0 1\n"; // memberList
    out << "1 2\n"; // numEachMol
    out << "1 0\n"; // lastNumberUpdateItrEachMol
    out << "0\n"; // Complex::emptyComList size

    /* ------------------------- OBSERVABLE SECTION ------------------------- */
    out << "#Observables\n"; // skipped
    out << "2\n"; // number of observables
    out << "obsA 5\n";
    out << "obsB 7\n";

    /* -------------------------- COUNTERS SECTION -------------------------- */
    out << "#Counters\n"; // skipped
    // nLoops nCancelOverlapPartner nCancelOverlapSystem nCancelDisplace2D
    // nCancelDisplace3D nCancelDisplace3Dto2D nCancelSpanBox nAssocSuccess
    // eventArraySize
    out << "3 4 5 6 7 8 9 10 " << kRREventArraySize << '\n';
    for (int i = 0; i < kRREventArraySize; ++i) // events3D
        out << i << ' ';
    out << '\n';
    for (int i = 0; i < kRREventArraySize; ++i) // events3Dto2D
        out << (100 + i) << ' ';
    out << '\n';
    for (int i = 0; i < kRREventArraySize; ++i) // events2D
        out << (200 + i) << ' ';
    out << '\n';
    out << "#Species\n"; // extra line skipped before numSpecies
    out << "2\n"; // number of species with bindPairLists
    out << "2\n"; // species 0: two bound pairs
    out << "0 1\n";
    out << "0\n"; // species 1: no bound pairs
    out << "\n"; // (empty value line for species 1)

    out.close();
}

// -----------------------------------------------------------------------------
// Helper: reset every *static* member that read_restart() overwrites/appends to,
// so repeated reads in this file (and other tests in the suite) start clean.
// -----------------------------------------------------------------------------
void rr_reset_statics()
{
    MolTemplate::numMolTypes = 0;
    MolTemplate::numEachMolType.clear();
    MolTemplate::absToRelIface.clear();
    Interface::State::totalNumOfStates = 0;
    Molecule::numberOfMolecules = 0;
    Molecule::emptyMolList.clear();
    Complex::numberOfComplexes = 0;
    Complex::emptyComList.clear();
    RxnBase::numberOfRxns = 0;
    RxnBase::totRxnSpecies = 0;
    Parameters::lastUpdateTransition.clear();
}

//! Bundle of every output container read_restart() needs.
struct RRSystem {
    long long int simItr { -1 };
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

// -----------------------------------------------------------------------------
// Helper: create the file, run read_restart() on it and report success.
//
// NOTE: the parser writes directly into counterArrays.events*[i] without
// resizing, exactly like the production code expects init_association_events()
// to have pre-allocated them.  We therefore size them here to avoid UB.
// -----------------------------------------------------------------------------
bool rr_load_system(RRSystem& sys)
{
    rr_write_test_restart_file(kRRRestartFileName);
    rr_reset_statics();

    // Pre-allocate the event histograms (see note above).
    sys.counterArrays.eventArraySize = kRREventArraySize;
    sys.counterArrays.events3D.assign(kRREventArraySize, -1);
    sys.counterArrays.events3Dto2D.assign(kRREventArraySize, -1);
    sys.counterArrays.events2D.assign(kRREventArraySize, -1);

    std::ifstream restartFile(kRRRestartFileName);
    if (!restartFile.is_open()) {
        ADD_FAILURE() << "Could not open the synthetic restart file "
                      << kRRRestartFileName;
        return false;
    }

    bool ok = true;
    try {
        read_restart(sys.simItr, restartFile, sys.params, sys.simulVolume,
            sys.moleculeList, sys.complexList, sys.molTemplateList,
            sys.forwardRxns, sys.backRxns, sys.createDestructRxns,
            sys.transmissionRxns, sys.observablesList, sys.membraneObject,
            sys.counterArrays);
    } catch (const std::exception& e) {
        // read_restart() only catches std::string / std::length_error, so guard
        // against anything else escaping (e.g. std::stod on a bad token).
        ADD_FAILURE() << "read_restart() threw an unexpected exception: " << e.what();
        ok = false;
    }

    restartFile.close();
    std::remove(kRRRestartFileName);
    return ok;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Parameters block.
// -----------------------------------------------------------------------------
void test_rr_parameters_block()
{
    std::cerr << "\n[TEST] test_rr_parameters_block\n"
              << "  Source file:   src/io/read_restart.cpp\n"
              << "  Function:      read_restart()\n"
              << "  Scenario:      parse the leading '#Parameters' block.\n"
              << "  Pass criteria: every scalar parameter equals the value that\n"
              << "                 was written into the synthetic restart file.\n";

    RRSystem sys;
    ASSERT_TRUE(rr_load_system(sys)) << "restart file could not be parsed";

    std::cerr << "  -> checking iteration counters and time step\n";
    EXPECT_EQ(sys.params.nItr, 1000) << "nItr should be 1000";
    EXPECT_EQ(sys.simItr, 250) << "simItr should be 250";
    EXPECT_EQ(sys.params.itrRestartFrom, 250)
        << "itrRestartFrom must be set from simItr";
    EXPECT_DOUBLE_EQ(sys.params.timeRestartFrom, 0.0025)
        << "timeRestartFrom should be 0.0025 s";
    EXPECT_DOUBLE_EQ(sys.params.timeStep, 0.1) << "timeStep should be 0.1 us";

    std::cerr << "  -> checking system-size parameters\n";
    EXPECT_EQ(sys.params.numMolTypes, 1) << "numMolTypes should be 1";
    EXPECT_EQ(sys.params.numTotalSpecies, 3) << "numTotalSpecies should be 3";
    EXPECT_EQ(sys.params.numTotalComplex, 1) << "numTotalComplex should be 1";
    EXPECT_EQ(sys.params.numTotalUnits, 4u) << "numTotalUnits should be 4";
    EXPECT_EQ(sys.params.numLipids, 0) << "numLipids should be 0";
    EXPECT_EQ(sys.params.max2DRxns, 100) << "max2DRxns should be 100";

    std::cerr << "  -> checking geometry/limit parameters\n";
    EXPECT_DOUBLE_EQ(sys.params.overlapSepLimit, 0.2)
        << "overlapSepLimit should be 0.2 nm";
    EXPECT_DOUBLE_EQ(sys.params.rMaxLimit, 5.5) << "rMaxLimit should be 5.5 nm";
    EXPECT_DOUBLE_EQ(sys.params.scaleMaxDisplace, 60.5)
        << "scaleMaxDisplace should be 60.5";

    std::cerr << "  -> checking IO intervals and flags\n";
    EXPECT_EQ(sys.params.timeWrite, 10) << "timeWrite should be 10";
    EXPECT_EQ(sys.params.trajWrite, 20) << "trajWrite should be 20";
    EXPECT_EQ(sys.params.restartWrite, 30) << "restartWrite should be 30";
    EXPECT_EQ(sys.params.pdbWrite, 40) << "pdbWrite should be 40";
    EXPECT_TRUE(sys.params.assocDissocWrite) << "assocDissocWrite should be true";
    EXPECT_EQ(sys.params.checkPoint, 50) << "checkPoint should be 50";
    EXPECT_EQ(sys.params.transitionWrite, 70) << "transitionWrite should be 70";
    EXPECT_TRUE(sys.params.clusterOverlapCheck)
        << "clusterOverlapCheck should be true";

    std::cerr << "  -> checking the optional (keyword) trailing parameters\n";
    EXPECT_FALSE(sys.params.rngwrite) << "RNGwrite was written as 0";
    EXPECT_EQ(sys.params.bondedComplexWrite, -1)
        << "bondedComplexWrite was written as -1";
    EXPECT_TRUE(Parameters::lastUpdateTransition.empty())
        << "lastUpdateTransition list length was 0 in the file";
}

// -----------------------------------------------------------------------------
// Test 2: Membrane / water box block.
// -----------------------------------------------------------------------------
void test_rr_membrane_block()
{
    std::cerr << "\n[TEST] test_rr_membrane_block\n"
              << "  Source file:   src/io/read_restart.cpp\n"
              << "  Function:      read_restart()\n"
              << "  Scenario:      parse the water box + implicit-lipid fields.\n"
              << "  Pass criteria: box dimensions, derived volume, membrane state\n"
              << "                 counts and boundary flags all match the file.\n";

    RRSystem sys;
    ASSERT_TRUE(rr_load_system(sys)) << "restart file could not be parsed";

    std::cerr << "  -> checking water box dimensions and derived volume\n";
    EXPECT_DOUBLE_EQ(sys.membraneObject.waterBox.x, 20.0) << "box x should be 20";
    EXPECT_DOUBLE_EQ(sys.membraneObject.waterBox.y, 30.0) << "box y should be 30";
    EXPECT_DOUBLE_EQ(sys.membraneObject.waterBox.z, 40.0) << "box z should be 40";
    EXPECT_DOUBLE_EQ(sys.membraneObject.waterBox.volume, 20.0 * 30.0 * 40.0)
        << "volume must be recomputed as x*y*z by read_restart";

    std::cerr << "  -> checking implicit lipid bookkeeping\n";
    EXPECT_EQ(sys.membraneObject.implicitlipidIndex, -1)
        << "implicitlipidIndex should be -1 (no implicit lipid)";
    EXPECT_EQ(sys.membraneObject.nSites, 0) << "nSites should be 0";
    EXPECT_EQ(sys.membraneObject.nStates, 1) << "nStates should be 1";
    EXPECT_EQ(sys.membraneObject.No_free_lipids, 12)
        << "No_free_lipids should be 12";
    EXPECT_EQ(sys.membraneObject.No_protein, 3) << "No_protein should be 3";
    EXPECT_DOUBLE_EQ(sys.membraneObject.totalSA, 1234.5)
        << "totalSA should be 1234.5";

    std::cerr << "  -> checking per-state free lipid list (nStates entries)\n";
    ASSERT_EQ(sys.membraneObject.numberOfFreeLipidsEachState.size(), 1u)
        << "one entry expected because nStates == 1";
    EXPECT_EQ(sys.membraneObject.numberOfFreeLipidsEachState[0], 12)
        << "the single state holds 12 free lipids";

    std::cerr << "  -> checking boundary-condition flags\n";
    EXPECT_FALSE(sys.membraneObject.implicitLipid) << "implicitLipid should be false";
    EXPECT_FALSE(sys.membraneObject.TwoD) << "TwoD should be false";
    EXPECT_TRUE(sys.membraneObject.isBox) << "isBox should be true";
    EXPECT_FALSE(sys.membraneObject.isSphere) << "isSphere should be false";
    EXPECT_DOUBLE_EQ(sys.membraneObject.sphereR, 0.0) << "sphereR should be 0";
    EXPECT_FALSE(sys.membraneObject.hasCompartment)
        << "hasCompartment should be false";
    EXPECT_DOUBLE_EQ(sys.membraneObject.compartmentR, 0.0)
        << "compartmentR should be 0";
}

// -----------------------------------------------------------------------------
// Test 3: MolTemplate block (template, interface and interface state).
// -----------------------------------------------------------------------------
void test_rr_moltemplate_block()
{
    std::cerr << "\n[TEST] test_rr_moltemplate_block\n"
              << "  Source file:   src/io/read_restart.cpp\n"
              << "  Function:      read_restart()\n"
              << "  Scenario:      parse one MolTemplate with one interface that\n"
              << "                 carries one state.\n"
              << "  Pass criteria: statics (numMolTypes, absToRelIface, ...) and\n"
              << "                 every template/interface/state field match.\n";

    RRSystem sys;
    ASSERT_TRUE(rr_load_system(sys)) << "restart file could not be parsed";

    std::cerr << "  -> checking MolTemplate statics\n";
    EXPECT_EQ(MolTemplate::numMolTypes, 1u) << "one molecule type in the file";
    ASSERT_EQ(MolTemplate::numEachMolType.size(), 1u)
        << "numEachMolType must have numMolTypes entries";
    EXPECT_EQ(MolTemplate::numEachMolType[0], 1) << "one copy of type 0";
    ASSERT_EQ(MolTemplate::absToRelIface.size(), 1u)
        << "absToRelIface list length was 1 in the file";
    EXPECT_EQ(MolTemplate::absToRelIface[0], 0) << "absolute iface 0 -> relative 0";
    EXPECT_EQ(Interface::State::totalNumOfStates, 1)
        << "totalNumOfStates should be 1";

    std::cerr << "  -> checking the parsed template itself\n";
    ASSERT_EQ(sys.molTemplateList.size(), 1u) << "exactly one template expected";
    const MolTemplate& tmp = sys.molTemplateList[0];
    EXPECT_EQ(tmp.molTypeIndex, 0) << "molTypeIndex should be 0";
    EXPECT_EQ(tmp.molName, std::string("testmol")) << "molName should be 'testmol'";
    EXPECT_EQ(tmp.copies, 5) << "copies should be 5";
    EXPECT_DOUBLE_EQ(tmp.mass, 1000.0) << "mass should be 1000";
    EXPECT_DOUBLE_EQ(tmp.radius, 2.5) << "radius should be 2.5";
    EXPECT_FALSE(tmp.isLipid) << "isLipid should be false";
    EXPECT_FALSE(tmp.isImplicitLipid) << "isImplicitLipid should be false";
    EXPECT_FALSE(tmp.isRod) << "isRod should be false";
    EXPECT_FALSE(tmp.isPoint) << "isPoint should be false";
    EXPECT_FALSE(tmp.checkOverlap) << "checkOverlap should be false";
    EXPECT_FALSE(tmp.countTransition) << "countTransition should be false";
    EXPECT_EQ(tmp.transitionMatrixSize, 500) << "transitionMatrixSize should be 500";
    EXPECT_EQ(tmp.transmissionRxnIndex, -1)
        << "transmissionRxnIndex should be -1 (unused)";

    std::cerr << "  -> checking template coordinates and diffusion constants\n";
    EXPECT_DOUBLE_EQ(tmp.comCoord.x, 0.0) << "template COM x should be 0";
    EXPECT_DOUBLE_EQ(tmp.D.x, 25.0) << "D.x should be 25";
    EXPECT_DOUBLE_EQ(tmp.D.z, 25.0) << "D.z should be 25";
    EXPECT_DOUBLE_EQ(tmp.Dr.y, 0.5) << "Dr.y should be 0.5";

    std::cerr << "  -> checking rxnPartners / bondList / monomerList\n";
    ASSERT_EQ(tmp.rxnPartners.size(), 1u) << "one reaction partner expected";
    EXPECT_EQ(tmp.rxnPartners[0], 0) << "partner molecule type should be 0";
    EXPECT_TRUE(tmp.bondList.empty()) << "no optional bonds were written";
    ASSERT_EQ(tmp.monomerList.size(), 2u) << "monomerList had two entries";
    EXPECT_EQ(tmp.monomerList[0], 3) << "first monomer index should be 3";
    EXPECT_EQ(tmp.monomerList[1], 4) << "second monomer index should be 4";
    EXPECT_TRUE(tmp.ifacesWithStates.empty())
        << "ifacesWithStates list was empty in the file";

    std::cerr << "  -> checking the interface and its single state\n";
    ASSERT_EQ(tmp.interfaceList.size(), 1u) << "one interface expected";
    const Interface& iface = tmp.interfaceList[0];
    EXPECT_EQ(iface.index, 0) << "interface index should be 0";
    EXPECT_EQ(iface.name, std::string("iface0")) << "interface name should be iface0";
    EXPECT_DOUBLE_EQ(iface.iCoord.x, 1.0) << "interface x coordinate should be 1";
    EXPECT_DOUBLE_EQ(iface.iCoord.y, 0.0) << "interface y coordinate should be 0";
    ASSERT_EQ(iface.stateList.size(), 1u) << "one state expected";
    const Interface::State& state = iface.stateList[0];
    EXPECT_EQ(state.index, 0) << "state index should be 0";
    EXPECT_EQ(state.iden, 'U') << "state identity char should be 'U'";
    ASSERT_EQ(state.rxnPartners.size(), 1u) << "one rxn partner for this state";
    EXPECT_EQ(state.rxnPartners[0], 1u) << "state partner index should be 1";
    ASSERT_EQ(state.myForwardRxns.size(), 1u) << "one forward reaction for this state";
    EXPECT_EQ(state.myForwardRxns[0], 0u) << "forward reaction index should be 0";
    EXPECT_TRUE(state.myCreateDestructRxns.empty())
        << "no create/destruct reactions were written";
    EXPECT_TRUE(state.stateChangeRxns.empty())
        << "no state-change reactions were written";
}

// -----------------------------------------------------------------------------
// Test 4: Reaction block (one ForwardRxn plus its BackRxn).
// -----------------------------------------------------------------------------
void test_rr_reaction_block()
{
    std::cerr << "\n[TEST] test_rr_reaction_block\n"
              << "  Source file:   src/io/read_restart.cpp\n"
              << "  Function:      read_restart()\n"
              << "  Scenario:      parse one reversible bimolecular reaction and\n"
              << "                 its conjugate back reaction.\n"
              << "  Pass criteria: reaction indices, angles (including the 'nan'\n"
              << "                 tokens), reactant/product lists and rates match.\n";

    RRSystem sys;
    ASSERT_TRUE(rr_load_system(sys)) << "restart file could not be parsed";

    std::cerr << "  -> checking reaction counters / list sizes\n";
    EXPECT_EQ(RxnBase::numberOfRxns, 2u) << "two reactions in total";
    EXPECT_EQ(RxnBase::totRxnSpecies, 4) << "totRxnSpecies should be 4";
    ASSERT_EQ(sys.forwardRxns.size(), 1u) << "one forward reaction expected";
    ASSERT_EQ(sys.backRxns.size(), 1u) << "one back reaction expected";
    EXPECT_TRUE(sys.createDestructRxns.empty())
        << "no create/destruct reactions were written";
    EXPECT_TRUE(sys.transmissionRxns.empty())
        << "no transmission reactions were written";

    const ForwardRxn& fwd = sys.forwardRxns[0];
    std::cerr << "  -> checking ForwardRxn scalar fields\n";
    EXPECT_EQ(fwd.absRxnIndex, 0) << "absRxnIndex should be 0";
    EXPECT_EQ(fwd.relRxnIndex, 0) << "relRxnIndex should be 0";
    EXPECT_EQ(fwd.rxnLabel, std::string("label0")) << "rxnLabel should be label0";
    EXPECT_EQ(fwd.rxnType, ReactionType::bimolecular)
        << "rxnType 1 must map to ReactionType::bimolecular";
    EXPECT_FALSE(fwd.isSymmetric) << "isSymmetric should be false";
    EXPECT_FALSE(fwd.isOnMem) << "isOnMem should be false";
    EXPECT_FALSE(fwd.hasStateChange) << "hasStateChange should be false";
    EXPECT_FALSE(fwd.isObserved) << "isObserved should be false";
    EXPECT_EQ(fwd.productName, std::string("A_B")) << "productName should be A_B";
    EXPECT_TRUE(fwd.isReversible) << "isReversible should be true";
    EXPECT_EQ(fwd.conjBackRxnIndex, 0) << "conjBackRxnIndex should be 0";
    EXPECT_FALSE(fwd.irrevRingClosure) << "irrevRingClosure should be false";
    EXPECT_DOUBLE_EQ(fwd.bindRadSameCom, 1.1) << "bindRadSameCom should be 1.1";
    EXPECT_DOUBLE_EQ(fwd.loopCoopFactor, 1.0) << "loopCoopFactor should be 1";
    EXPECT_DOUBLE_EQ(fwd.length3Dto2D, 2.0) << "length3Dto2D should be 2";
    EXPECT_DOUBLE_EQ(fwd.area3Dto1D, 4.0) << "area3Dto1D should be 4";
    EXPECT_FALSE(fwd.excludeVolumeBound) << "excludeVolumeBound should be false";
    EXPECT_FALSE(fwd.isCoupled) << "isCoupled should be false";

    std::cerr << "  -> checking binding radius, association angles and normals\n";
    EXPECT_DOUBLE_EQ(fwd.bindRadius, 1.5) << "bindRadius should be 1.5";
    EXPECT_DOUBLE_EQ(fwd.assocAngles.theta1, 1.5708) << "theta1 should be 1.5708";
    EXPECT_DOUBLE_EQ(fwd.assocAngles.theta2, 1.5708) << "theta2 should be 1.5708";
    EXPECT_TRUE(std::isnan(fwd.assocAngles.phi1))
        << "the 'nan' token must become quiet_NaN for phi1";
    EXPECT_TRUE(std::isnan(fwd.assocAngles.phi2))
        << "the 'nan' token must become quiet_NaN for phi2";
    EXPECT_DOUBLE_EQ(fwd.assocAngles.omega, 3.14159) << "omega should be 3.14159";
    EXPECT_DOUBLE_EQ(fwd.norm1.z, 1.0) << "norm1 should be (0,0,1)";
    EXPECT_DOUBLE_EQ(fwd.norm2.z, 1.0) << "norm2 should be (0,0,1)";

    std::cerr << "  -> checking ForwardRxn reactant/product/rate lists\n";
    ASSERT_EQ(fwd.intReactantList.size(), 2u) << "two integer reactants";
    EXPECT_EQ(fwd.intReactantList[0], 0) << "first integer reactant is 0";
    EXPECT_EQ(fwd.intReactantList[1], 1) << "second integer reactant is 1";
    ASSERT_EQ(fwd.intProductList.size(), 1u) << "one integer product";
    EXPECT_EQ(fwd.intProductList[0], 2) << "integer product is 2";
    ASSERT_EQ(fwd.reactantListNew.size(), 2u) << "two RxnIface reactants";
    EXPECT_EQ(fwd.reactantListNew[0].ifaceName, std::string("iface0"))
        << "reactant 0 interface name";
    EXPECT_EQ(fwd.reactantListNew[1].absIfaceIndex, 1)
        << "reactant 1 absolute interface index should be 1";
    EXPECT_EQ(fwd.reactantListNew[0].requiresState, 'U')
        << "reactant 0 requiresState should be 'U'";
    EXPECT_FALSE(fwd.reactantListNew[0].requiresInteraction)
        << "reactant 0 requiresInteraction should be false";
    ASSERT_EQ(fwd.productListNew.size(), 1u) << "one RxnIface product";
    EXPECT_EQ(fwd.productListNew[0].absIfaceIndex, 2)
        << "product absolute interface index should be 2";
    EXPECT_TRUE(fwd.productListNew[0].requiresInteraction)
        << "product requiresInteraction should be true";
    ASSERT_EQ(fwd.rateList.size(), 1u) << "one rate state";
    EXPECT_DOUBLE_EQ(fwd.rateList[0].rate, 10.5) << "forward rate should be 10.5";
    EXPECT_TRUE(fwd.rateList[0].otherIfaceLists.empty())
        << "no ancillary interface lists were written";

    const BackRxn& back = sys.backRxns[0];
    std::cerr << "  -> checking BackRxn fields\n";
    EXPECT_EQ(back.absRxnIndex, 1) << "back reaction absRxnIndex should be 1";
    EXPECT_EQ(back.relRxnIndex, 0) << "back reaction relRxnIndex should be 0";
    EXPECT_EQ(back.conjForwardRxnIndex, 0u)
        << "conjForwardRxnIndex should point at forward reaction 0";
    EXPECT_FALSE(back.isCoupled) << "back reaction isCoupled should be false";
    ASSERT_EQ(back.reactantListNew.size(), 1u) << "one back-reaction reactant";
    ASSERT_EQ(back.productListNew.size(), 2u) << "two back-reaction products";
    ASSERT_EQ(back.rateList.size(), 1u) << "one back-reaction rate";
    EXPECT_DOUBLE_EQ(back.rateList[0].rate, 0.5) << "back rate should be 0.5";
}

// -----------------------------------------------------------------------------
// Test 5: Molecule block, including the bound-interface (Interaction) branch.
// -----------------------------------------------------------------------------
void test_rr_molecule_block()
{
    std::cerr << "\n[TEST] test_rr_mol