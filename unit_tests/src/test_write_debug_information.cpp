/*! \file test_write_debug_information.cpp
 *
 * ### Unit test for ../src/debug/write_debug_information.cpp
 *
 * The single function defined in that translation unit is:
 *
 * \code
 * void write_debug_information(MpiContext& mpiContext, long long int simItr,
 *                              ofstream& debugFile,
 *                              vector<Molecule>& moleculeList,
 *                              vector<Complex>& complexList,
 *                              vector<MolTemplate>& molTemplateList,
 *                              copyCounters& counterArrays, double simTime,
 *                              string s);
 * \endcode
 *
 * In the current version of the code base the body of this function is
 * intentionally *empty* -- it is a hook that developers fill in temporarily
 * when they need to dump per-iteration state while chasing a bug.  Its
 * contract, as shipped, is therefore very specific and completely testable:
 *
 *   1. It must be callable with fully populated simulation containers and
 *      must return normally (no crash, no exception).
 *   2. It must be callable with *empty* containers (it must not index into
 *      moleculeList / complexList / molTemplateList).
 *   3. It must not mutate any of the objects passed by (non-const) reference.
 *   4. It must not write anything to the supplied debug stream (the stream is
 *      passed by reference but the empty body never touches it).
 *   5. It must tolerate arbitrary scalar arguments (negative / huge simItr,
 *      negative simTime, empty label string).
 *
 * Every test below prints, to stderr, which source file and function is under
 * test, exactly what data is being handed to it, and what the pass criteria
 * are.
 */

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// Pulls in Molecule, Complex, MolTemplate, copyCounters, SimulVolume,
// Membrane, MpiContext (via split.cpp) and the declaration of
// write_debug_information() itself.
#include "debug/debug.hpp"

// -----------------------------------------------------------------------------
// Local helpers.  Everything is given the "wdi_" prefix so that no symbol in
// this translation unit can collide with another generated test file that is
// linked into the same gtest binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully-initialised MolTemplate with a single interface.
 *
 * Every field that the simulation code normally reads is filled in so that
 * nothing downstream can dereference an uninitialised member.
 *
 * \param[in] typeIndex The molTypeIndex (also used as the interface state index).
 * \param[in] name      Human readable molecule name.
 */
MolTemplate wdi_make_mol_template(int typeIndex, const std::string& name)
{
    MolTemplate molTemplate;
    molTemplate.molName = name;
    molTemplate.molTypeIndex = typeIndex;
    molTemplate.copies = 1;
    molTemplate.mass = 1.0;
    molTemplate.radius = 1.0;
    molTemplate.comCoord = Coord{ 0.0, 0.0, 0.0 };
    molTemplate.D = Coord{ 10.0, 10.0, 10.0 };
    molTemplate.Dr = Coord{ 0.01, 0.01, 0.01 };
    molTemplate.isLipid = false;
    molTemplate.isImplicitLipid = false;
    molTemplate.isPoint = false;
    molTemplate.isRod = false;

    // One interface, offset 1 nm along +x, with one (default) state.
    Interface iface(name + "_iface", Coord{ 1.0, 0.0, 0.0 });
    iface.index = 0;
    iface.stateList.emplace_back(name + "_iface", typeIndex);
    molTemplate.interfaceList.push_back(iface);

    return molTemplate;
}

/*! \brief Build a fully-initialised Molecule owning a single interface.
 *
 * \param[in] index      Index of this molecule inside moleculeList.
 * \param[in] comIndex   Index of the parent complex inside complexList.
 * \param[in] typeIndex  molTypeIndex, matching a MolTemplate.
 * \param[in] com        Centre-of-mass coordinate.
 */
Molecule wdi_make_molecule(int index, int comIndex, int typeIndex, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.id = index; // unique id, only meaningful for the MPI build
    mol.complexId = comIndex;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = typeIndex;
    mol.mySubVolIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    // A single unbound interface located 1 nm along +x from the COM.
    Molecule::Iface iface;
    iface.coord = Coord{ com.x + 1.0, com.y, com.z };
    iface.relIndex = 0;
    iface.index = typeIndex;
    iface.stateIndex = 0;
    iface.molTypeIndex = typeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    // Legacy encounter/bookkeeping lists, left empty but explicitly cleared.
    mol.freelist.clear();
    mol.bndlist.clear();
    mol.bndpartner.clear();
    mol.bndRxnList.clear();

    return mol;
}

/*! \brief Build a Complex that owns exactly one member molecule.
 *
 * The Complex is assembled field-by-field rather than through the
 * Complex(Molecule, MolTemplate) constructor on purpose: that constructor
 * mutates the static counters Complex::maxID / MolTemplate::numMolTypes, and a
 * unit test should not perturb global simulation state shared with the rest of
 * the test binary.
 *
 * \param[in] index    Index of this complex inside complexList.
 * \param[in] memberMol Index (in moleculeList) of the single member molecule.
 * \param[in] com      Centre-of-mass coordinate.
 */
Complex wdi_make_complex(int index, int memberMol, const Coord& com)
{
    Complex targCom;
    targCom.index = index;
    targCom.id = index;
    targCom.comCoord = com;
    targCom.mass = 1.0;
    targCom.radius = 1.0;
    targCom.D = Coord{ 10.0, 10.0, 10.0 };
    targCom.Dr = Coord{ 0.01, 0.01, 0.01 };
    targCom.isEmpty = false;
    targCom.OnSurface = false;
    targCom.onFiber = false;
    targCom.trajStatus = TrajStatus::none;
    targCom.memberList.clear();
    targCom.memberList.push_back(memberMol);
    targCom.numEachMol.assign(1, 1);
    targCom.lastNumberUpdateItrEachMol.assign(1, 0);
    return targCom;
}

/*! \brief Build a copyCounters object with a few non-default entries.
 *
 * The specific numbers are arbitrary; they exist purely so we can verify that
 * write_debug_information() leaves them untouched.
 */
copyCounters wdi_make_counters()
{
    copyCounters counterArrays;
    counterArrays.nBoundPairs = std::vector<int>{ 3, 5 };
    counterArrays.proPairlist = std::vector<int>{ 0, 1 };
    counterArrays.copyNumSpecies = std::vector<int>{ 7, 11, 13 };
    counterArrays.singleDouble = std::vector<int>{ 1, 0, 1 };
    counterArrays.nLoops = 2;
    counterArrays.nAssocSuccess = 17;
    counterArrays.nCancelSpanBox = 4;
    counterArrays.eventArraySize = 20;
    counterArrays.events3D.assign(counterArrays.eventArraySize, 0);
    counterArrays.events2D.assign(counterArrays.eventArraySize, 0);
    counterArrays.events3Dto2D.assign(counterArrays.eventArraySize, 0);
    return counterArrays;
}

/*! \brief Return the size, in bytes, of a file on disk (0 if it cannot be read). */
std::streamoff wdi_file_size(const std::string& fileName)
{
    std::ifstream in(fileName, std::ios::binary | std::ios::ate);
    if (!in.is_open())
        return 0;
    return in.tellg();
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the function can be called with a fully populated simulation state
//         and returns normally.
// -----------------------------------------------------------------------------
void test_wdi_callable_with_populated_state()
{
    std::cerr << "\n[TEST] test_wdi_callable_with_populated_state\n"
              << "  Source file:   src/debug/write_debug_information.cpp\n"
              << "  Function:      write_debug_information()\n"
              << "  Scenario:      2 MolTemplates, 2 Molecules, 2 Complexes, a\n"
              << "                 populated copyCounters and an open debug file.\n"
              << "  Pass criteria: the call returns normally (no crash / throw)\n"
              << "                 and the containers still have their original\n"
              << "                 sizes afterwards.\n";

    // --- Build the simulation state -----------------------------------------
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wdi_make_mol_template(0, "A"));
    molTemplateList.push_back(wdi_make_mol_template(1, "B"));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wdi_make_molecule(0, 0, 0, Coord{ 0.0, 0.0, 0.0 }));
    moleculeList.push_back(wdi_make_molecule(1, 1, 1, Coord{ 5.0, -5.0, 2.5 }));

    std::vector<Complex> complexList;
    complexList.push_back(wdi_make_complex(0, 0, Coord{ 0.0, 0.0, 0.0 }));
    complexList.push_back(wdi_make_complex(1, 1, Coord{ 5.0, -5.0, 2.5 }));

    copyCounters counterArrays = wdi_make_counters();

    // Supporting objects that the (currently empty) debug helper would use if
    // its body were filled in.  A default MpiContext is created and only the
    // fields that debug.hpp / mpi_function.hpp actually read are populated.
    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double>{ 100.0, 100.0, 100.0 });

    Parameters params;
    params.timeStep = 1.0;
    params.rMaxLimit = 10.0;

    SimulVolume simulVolume;
    simulVolume.create_simulation_volume(params, membraneObject);

    MpiContext mpiContext;
    mpiContext.rank = 0;
    mpiContext.nprocs = 1;
    mpiContext.xOffset = 0;
    mpiContext.membraneObject = &membraneObject;
    mpiContext.simulVolume = &simulVolume;

    const std::string fileName = "wdi_test_populated.debug";
    std::ofstream debugFile(fileName, std::ios::trunc);
    ASSERT_TRUE(debugFile.is_open()) << "Could not open scratch debug file " << fileName;

    std::cerr << "  Calling write_debug_information(simItr=42, simTime=0.5, s=\"populated\")...\n";
    write_debug_information(mpiContext, 42, debugFile, moleculeList, complexList,
                            molTemplateList, counterArrays, 0.5, std::string("populated"));
    std::cerr << "  ...returned normally.\n";

    // The function must not have added or removed elements from any container.
    EXPECT_EQ(moleculeList.size(), static_cast<size_t>(2))
        << "moleculeList must still hold 2 molecules";
    EXPECT_EQ(complexList.size(), static_cast<size_t>(2))
        << "complexList must still hold 2 complexes";
    EXPECT_EQ(molTemplateList.size(), static_cast<size_t>(2))
        << "molTemplateList must still hold 2 templates";

    // The stream must still be usable (not put into a fail state).
    EXPECT_TRUE(debugFile.good()) << "debugFile stream should still be in a good state";

    debugFile.close();
    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 2: the function must not mutate any of its by-reference arguments.
// -----------------------------------------------------------------------------
void test_wdi_does_not_modify_inputs()
{
    std::cerr << "\n[TEST] test_wdi_does_not_modify_inputs\n"
              << "  Source file:   src/debug/write_debug_information.cpp\n"
              << "  Function:      write_debug_information()\n"
              << "  Scenario:      snapshot molecule/complex coordinates and the\n"
              << "                 copyCounters values, call the function, then\n"
              << "                 compare everything again.\n"
              << "  Pass criteria: every observed field is bit-for-bit identical\n"
              << "                 before and after the call (the shipped body is\n"
              << "                 empty, so it is a pure no-op).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wdi_make_mol_template(0, "A"));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wdi_make_molecule(0, 0, 0, Coord{ 1.25, -2.5, 3.75 }));

    std::vector<Complex> complexList;
    complexList.push_back(wdi_make_complex(0, 0, Coord{ 1.25, -2.5, 3.75 }));

    copyCounters counterArrays = wdi_make_counters();

    // --- Snapshot the state we are going to compare against ------------------
    const Coord molComBefore = moleculeList[0].comCoord;
    const Coord ifaceBefore = moleculeList[0].interfaceList[0].coord;
    const int molIndexBefore = moleculeList[0].index;
    const int molComIndexBefore = moleculeList[0].myComIndex;
    const bool molEmptyBefore = moleculeList[0].isEmpty;

    const Coord comComBefore = complexList[0].comCoord;
    const double comRadiusBefore = complexList[0].radius;
    const size_t comMembersBefore = complexList[0].memberList.size();

    const std::vector<int> boundPairsBefore = counterArrays.nBoundPairs;
    const std::vector<int> copyNumBefore = counterArrays.copyNumSpecies;
    const int nLoopsBefore = counterArrays.nLoops;
    const int nAssocBefore = counterArrays.nAssocSuccess;

    const std::string templateNameBefore = molTemplateList[0].molName;
    const double templateRadiusBefore = molTemplateList[0].radius;

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double>{ 50.0, 50.0, 50.0 });

    Parameters params;
    params.timeStep = 1.0;
    params.rMaxLimit = 10.0;

    SimulVolume simulVolume;
    simulVolume.create_simulation_volume(params, membraneObject);

    MpiContext mpiContext;
    mpiContext.rank = 0;
    mpiContext.nprocs = 1;
    mpiContext.xOffset = 0;
    mpiContext.membraneObject = &membraneObject;
    mpiContext.simulVolume = &simulVolume;

    const std::string fileName = "wdi_test_nomutate.debug";
    std::ofstream debugFile(fileName, std::ios::trunc);
    ASSERT_TRUE(debugFile.is_open()) << "Could not open scratch debug file " << fileName;

    std::cerr << "  Calling write_debug_information(simItr=7, simTime=1.5e-3, s=\"snapshot\")...\n";
    write_debug_information(mpiContext, 7, debugFile, moleculeList, complexList,
                            molTemplateList, counterArrays, 1.5e-3, std::string("snapshot"));

    // --- Molecule must be untouched -----------------------------------------
    std::cerr << "  Comparing Molecule state before/after...\n";
    EXPECT_TRUE(moleculeList[0].comCoord == molComBefore)
        << "Molecule centre of mass must be unchanged";
    EXPECT_TRUE(moleculeList[0].interfaceList[0].coord == ifaceBefore)
        << "Molecule interface coordinate must be unchanged";
    EXPECT_EQ(moleculeList[0].index, molIndexBefore) << "Molecule::index must be unchanged";
    EXPECT_EQ(moleculeList[0].myComIndex, molComIndexBefore)
        << "Molecule::myComIndex must be unchanged";
    EXPECT_EQ(moleculeList[0].isEmpty, molEmptyBefore) << "Molecule::isEmpty must be unchanged";

    // --- Complex must be untouched ------------------------------------------
    std::cerr << "  Comparing Complex state before/after...\n";
    EXPECT_TRUE(complexList[0].comCoord == comComBefore)
        << "Complex centre of mass must be unchanged";
    EXPECT_DOUBLE_EQ(complexList[0].radius, comRadiusBefore)
        << "Complex::radius must be unchanged";
    EXPECT_EQ(complexList[0].memberList.size(), comMembersBefore)
        << "Complex::memberList size must be unchanged";

    // --- copyCounters must be untouched -------------------------------------
    std::cerr << "  Comparing copyCounters state before/after...\n";
    EXPECT_EQ(counterArrays.nBoundPairs, boundPairsBefore)
        << "copyCounters::nBoundPairs must be unchanged";
    EXPECT_EQ(counterArrays.copyNumSpecies, copyNumBefore)
        << "copyCounters::copyNumSpecies must be unchanged";
    EXPECT_EQ(counterArrays.nLoops, nLoopsBefore) << "copyCounters::nLoops must be unchanged";
    EXPECT_EQ(counterArrays.nAssocSuccess, nAssocBefore)
        << "copyCounters::nAssocSuccess must be unchanged";

    // --- MolTemplate must be untouched --------------------------------------
    std::cerr << "  Comparing MolTemplate state before/after...\n";
    EXPECT_EQ(molTemplateList[0].molName, templateNameBefore)
        << "MolTemplate::molName must be unchanged";
    EXPECT_DOUBLE_EQ(molTemplateList[0].radius, templateRadiusBefore)
        << "MolTemplate::radius must be unchanged";

    debugFile.close();
    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 3: the empty body must not emit a single byte to the debug stream.
// -----------------------------------------------------------------------------
void test_wdi_writes_nothing_to_debug_file()
{
    std::cerr << "\n[TEST] test_wdi_writes_nothing_to_debug_file\n"
              << "  Source file:   src/debug/write_debug_information.cpp\n"
              << "  Function:      write_debug_information()\n"
              << "  Scenario:      hand the function a freshly truncated file\n"
              << "                 stream and call it three times in a row.\n"
              << "  Pass criteria: the stream write position stays at 0 and the\n"
              << "                 file on disk is 0 bytes long afterwards\n"
              << "                 (the shipped body writes nothing).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wdi_make_mol_template(0, "A"));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wdi_make_molecule(0, 0, 0, Coord{ 0.0, 0.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(wdi_make_complex(0, 0, Coord{ 0.0, 0.0, 0.0 }));

    copyCounters counterArrays = wdi_make_counters();

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double>{ 50.0, 50.0, 50.0 });

    Parameters params;
    params.timeStep = 1.0;
    params.rMaxLimit = 10.0;

    SimulVolume simulVolume;
    simulVolume.create_simulation_volume(params, membraneObject);

    MpiContext mpiContext;
    mpiContext.rank = 0;
    mpiContext.nprocs = 1;
    mpiContext.xOffset = 0;
    mpiContext.membraneObject = &membraneObject;
    mpiContext.simulVolume = &simulVolume;

    const std::string fileName = "wdi_test_output.debug";
    std::ofstream debugFile(fileName, std::ios::trunc);
    ASSERT_TRUE(debugFile.is_open()) << "Could not open scratch debug file " << fileName;

    // Call it repeatedly with different iteration numbers / labels.
    for (long long int simItr = 0; simItr < 3; ++simItr) {
        std::cerr << "  Calling write_debug_information(simItr=" << simItr << ")...\n";
        write_debug_information(mpiContext, simItr, debugFile, moleculeList, complexList,
                                molTemplateList, counterArrays,
                                static_cast<double>(simItr) * 1e-6,
                                std::string("iteration ") + std::to_string(simItr));
    }

    // The put pointer must not have advanced.
    const std::streamoff putPos = debugFile.tellp();
    std::cerr << "  Stream put position after 3 calls = " << putPos << " (expected 0)\n";
    EXPECT_EQ(putPos, static_cast<std::streamoff>(0))
        << "write_debug_information() must not advance the stream put position";

    debugFile.flush();
    debugFile.close();

    const std::streamoff sizeOnDisk = wdi_file_size(fileName);
    std::cerr << "  File size on disk = " << sizeOnDisk << " bytes (expected 0)\n";
    EXPECT_EQ(sizeOnDisk, static_cast<std::streamoff>(0))
        << "The debug file should still be empty after calling the no-op hook";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 4: the function must survive completely empty containers.
// -----------------------------------------------------------------------------
void test_wdi_handles_empty_containers()
{
    std::cerr << "\n[TEST] test_wdi_handles_empty_containers\n"
              << "  Source file:   src/debug/write_debug_information.cpp\n"
              << "  Function:      write_debug_information()\n"
              << "  Scenario:      empty moleculeList, complexList, molTemplateList\n"
              << "                 and a default-constructed copyCounters.\n"
              << "  Pass criteria: the call returns normally and all containers\n"
              << "                 are still empty (no indexing into element 0).\n";

    std::vector<Molecule> moleculeList;      // deliberately empty
    std::vector<Complex> complexList;        // deliberately empty
    std::vector<MolTemplate> molTemplateList; // deliberately empty
    copyCounters counterArrays;              // all vectors empty, counters at 0

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double>{ 20.0, 20.0, 20.0 });

    Parameters params;
    params.timeStep = 1.0;
    params.rMaxLimit = 5.0;

    SimulVolume simulVolume;
    simulVolume.create_simulation_volume(params, membraneObject);

    MpiContext mpiContext;
    mpiContext.rank = 0;
    mpiContext.nprocs = 1;
    mpiContext.xOffset = 0;
    mpiContext.membraneObject = &membraneObject;
    mpiContext.simulVolume = &simulVolume;

    const std::string fileName = "wdi_test_empty.debug";
    std::ofstream debugFile(fileName, std::ios::trunc);
    ASSERT_TRUE(debugFile.is_open()) << "Could not open scratch debug file " << fileName;

    std::cerr << "  Calling write_debug_information with empty containers...\n";
    write_debug_information(mpiContext, 0, debugFile, moleculeList, complexList,
                            molTemplateList, counterArrays, 0.0, std::string("empty"));
    std::cerr << "  ...returned normally.\n";

    EXPECT_TRUE(moleculeList.empty()) << "moleculeList must remain empty";
    EXPECT_TRUE(complexList.empty()) << "complexList must remain empty";
    EXPECT_TRUE(molTemplateList.empty()) << "molTemplateList must remain empty";
    EXPECT_TRUE(counterArrays.copyNumSpecies.empty())
        << "copyCounters::copyNumSpecies must remain empty";
    EXPECT_EQ(counterArrays.nLoops, 0) << "copyCounters::nLoops must remain 0";

    debugFile.close();
    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test 5: extreme / degenerate scalar arguments must be tolerated.
// -----------------------------------------------------------------------------
void test_wdi_handles_varied_scalar_arguments()
{
    std::cerr << "\n[TEST] test_wdi_handles_varied_scalar_arguments\n"
              << "  Source file:   src/debug/write_debug_information.cpp\n"
              << "  Function:      write_debug_information()\n"
              << "  Scenario:      call the hook with negative, zero and very large\n"
              << "                 simItr values, negative/huge simTime values and\n"
              << "                 an empty label string.\n"
              << "  Pass criteria: every call returns normally, the debug stream is\n"
              << "                 still good, and nothing was written.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(wdi_make_mol_template(0, "A"));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wdi_make_molecule(0, 0, 0, Coord{ -1.0, 2.0, -3.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(wdi_make_complex(0, 0, Coord{ -1.0, 2.0, -3.0 }));

    copyCounters counterArrays = wdi_make_counters();

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double>{ 30.0, 30.0, 30.0 });

    Parameters params;
    params.timeStep = 1.0;
    params.rMaxLimit = 5.0;

    SimulVolume simulVolume;
    simulVolume.create_simulation_volume(params, membraneObject);

    MpiContext mpiContext;
    mpiContext.rank = 0;
    mpiContext.nprocs = 1;
    mpiContext.xOffset = 0;
    mpiContext.membraneObject = &membraneObject;
    mpiContext.simulVolume = &simulVolume;

    const std::string fileName = "wdi_test_scalars.debug";
    std::ofstream debugFile(fileName, std::ios::trunc);
    ASSERT_TRUE(debugFile.is_open()) << "Could not open scratch debug file " << fileName;

    // A small table of degenerate (simItr, simTime, label) triples.
    struct Case {
        long long int simItr;
        double simTime;
        const char* label;
    };
    const Case cases[] = {
        { -1, -1.0, "negative iteration" },
        { 0, 0.0, "" },
        { std::numeric_limits<long long int>::max(), 1.0e12, "huge iteration" },
        { 123456789LL, -0.0, "negative zero time" },
    };

    for (const Case& oneCase : cases) {
        std::cerr << "  Calling write_debug_information(simItr=" << oneCase.simItr
                  << ", simTime=" << oneCase.simTime << ", s=\"" << oneCase.label << "\")...\n";
        write_debug_information(mpiContext, oneCase.simItr, debugFile, moleculeList,
                                complexList, molTemplateList, counterArrays,
                                oneCase.simTime, std::string(oneCase.label));
        // The stream must never be pushed into a fail/bad state by the hook.
        EXPECT_TRUE(debugFile.good())
            << "debugFile stream should still be good after simItr=" << oneCase.simItr;
    }

    // Still nothing written for any of the degenerate calls.
    EXPECT_EQ(debugFile.tellp(), static_cast<std::streamoff>(0))
        << "No bytes should have been written for any of the degenerate calls";

    debugFile.close();
    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so a
// failure in one scenario does not prevent the remaining scenarios from running
// (all assertions above are non-fatal EXPECT_* except for the file-open guards).
// -----------------------------------------------------------------------------
TEST(WriteDebugInformationTest, CallableWithPopulatedState)
{
    test_wdi_callable_with_populated_state();
}

TEST(WriteDebugInformationTest, DoesNotModifyInputs)
{
    test_wdi_does_not_modify_inputs();
}

TEST(WriteDebugInformationTest, WritesNothingToDebugFile)
{
    test_wdi_writes_nothing_to_debug_file();
}

TEST(WriteDebugInformationTest, HandlesEmptyContainers)
{
    test_wdi_handles_empty_containers();
}

TEST(WriteDebugInformationTest, HandlesVariedScalarArguments)
{
    test_wdi_handles_varied_scalar_arguments();
}