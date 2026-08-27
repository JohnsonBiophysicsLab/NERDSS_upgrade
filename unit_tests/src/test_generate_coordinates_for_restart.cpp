/*! \file test_generate_coordinates_for_restart.cpp
 *
 * ### Unit tests for src/system_setup/generate_coordinates_for_restart.cpp
 *
 * The file under test contains four free functions that are used when a
 * simulation is restarted and *new* molecules (declared in `add.inp`) have to
 * be inserted into an already existing system:
 *
 *   - generate_coordinates_for_restart()      : loops over the newly added
 *                                               MolTemplates and creates the
 *                                               requested number of copies.
 *   - create_molecule_and_complex_for_restart(): allocates (or recycles) a slot
 *                                               in moleculeList/complexList and
 *                                               fills it with a new species.
 *   - initialize_molecule_for_restart()       : builds a single Molecule with
 *                                               random coordinates.
 *   - moleculeOverlapsForRestart()            : decides whether a freshly placed
 *                                               molecule has to be resampled.
 *
 * All four are exercised below.  Because the code under test manipulates a
 * number of *static* counters that are shared with the rest of the test suite
 * (Molecule::numberOfMolecules, Complex::emptyComList, MolTemplate::numEachMolType,
 * ...), every test snapshots those statics on entry and restores them on exit.
 *
 * Verbose progress information is written to stderr so that a reader of the
 * console log can follow which function is being probed and what the pass
 * criterion for every assertion is.
 */

#include "system_setup/system_setup.hpp"
#include "math/rand_gsl.hpp"

#include <gtest/gtest.h>
#include <gsl/gsl_rng.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// The GSL random number generator is defined in unit_tests/src/gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// RNG helper.  create_random_coords() (called from initialize_molecule_for_restart)
// uses rand_gsl(), so the generator must exist and be seeded deterministically.
// -----------------------------------------------------------------------------
void gcfr_seed_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, 42); // fixed seed -> reproducible coordinates
}

// -----------------------------------------------------------------------------
// RAII guard that snapshots every static counter touched by the code under test
// and puts it back the way it was, so neighbouring tests in the suite are not
// disturbed.
// -----------------------------------------------------------------------------
struct GcfrStaticsGuard {
    int savedMolCount;
    int savedMolMaxID;
    int savedComCount;
    int savedComMaxID;
    unsigned savedNumMolTypes;
    std::vector<int> savedNumEachMolType;
    std::vector<int> savedEmptyMolList;
    std::vector<int> savedEmptyComList;

    GcfrStaticsGuard()
    {
        savedMolCount = Molecule::numberOfMolecules;
        savedMolMaxID = Molecule::maxID;
        savedComCount = Complex::numberOfComplexes;
        savedComMaxID = Complex::maxID;
        savedNumMolTypes = MolTemplate::numMolTypes;
        savedNumEachMolType = MolTemplate::numEachMolType;
        savedEmptyMolList = Molecule::emptyMolList;
        savedEmptyComList = Complex::emptyComList;
    }

    ~GcfrStaticsGuard()
    {
        Molecule::numberOfMolecules = savedMolCount;
        Molecule::maxID = savedMolMaxID;
        Complex::numberOfComplexes = savedComCount;
        Complex::maxID = savedComMaxID;
        MolTemplate::numMolTypes = savedNumMolTypes;
        MolTemplate::numEachMolType = savedNumEachMolType;
        Molecule::emptyMolList = savedEmptyMolList;
        Complex::emptyComList = savedEmptyComList;
    }
};

/*! \brief Puts the shared statics into a clean, well defined starting state.
 *
 * \param[in] numTypes number of MolTemplates the test will use.  Both
 *            MolTemplate::numEachMolType and Complex::numEachMol are indexed by
 *            molTypeIndex, so these have to be large enough or the code under
 *            test writes out of bounds.
 */
void gcfr_prepare_statics(unsigned numTypes)
{
    MolTemplate::numMolTypes = numTypes;
    MolTemplate::numEachMolType.assign(numTypes, 0);
    Molecule::emptyMolList.clear();
    Complex::emptyComList.clear();
    Molecule::numberOfMolecules = 0;
    Complex::numberOfComplexes = 0;
}

/*! \brief Cubic, non-spherical water box centred on the origin. */
Membrane gcfr_make_membrane(double side)
{
    Membrane mem;
    mem.waterBox = Membrane::WaterBox(std::vector<double> { side, side, side });
    mem.isSphere = false;
    mem.implicitLipid = false;
    mem.hasCompartment = false;
    return mem;
}

/*! \brief Fully initialised MolTemplate for a simple soluble (3D) molecule.
 *
 * Every interface is placed `ifaceLen` nm away from the centre of mass and owns
 * exactly one state, which is what create_random_coords() expects
 * (it dereferences stateList[0] unconditionally).
 */
MolTemplate gcfr_make_template(
    const std::string& name, int typeIndex, int copies, int nIfaces, double ifaceLen)
{
    MolTemplate temp {};
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.copies = copies;
    temp.mass = 1.0;
    temp.radius = 2.0;
    temp.D = Coord(10.0, 10.0, 10.0);
    temp.Dr = Coord(0.1, 0.1, 0.1);
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPromoter = false;
    temp.isRod = false;
    temp.isPoint = false;
    temp.checkOverlap = false;

    for (int i = 0; i < nIfaces; ++i) {
        Interface iface {};
        iface.index = i;
        iface.name = std::string("i") + std::to_string(i);
        // spread the interfaces over the three axes, all at distance ifaceLen
        if (i % 3 == 0)
            iface.iCoord = Coord(ifaceLen, 0.0, 0.0);
        else if (i % 3 == 1)
            iface.iCoord = Coord(0.0, ifaceLen, 0.0);
        else
            iface.iCoord = Coord(0.0, 0.0, ifaceLen);
        iface.stateList.emplace_back(i); // Interface::State(int index), iden == '\0'
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Hand-built, fully initialised Molecule with exactly one interface. */
Molecule gcfr_make_mol(int index, int comIndex, int molTypeIndex, const Coord& com,
    const Coord& ifaceCoord, int absIfaceIndex)
{
    Molecule mol {};
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.comCoord = com;

    Molecule::Iface ifc {};
    ifc.coord = ifaceCoord;
    ifc.index = absIfaceIndex; // absolute (state) index
    ifc.relIndex = 0;          // relative index within the molecule
    ifc.stateIndex = 0;
    ifc.stateIden = '\0';
    ifc.molTypeIndex = molTypeIndex;
    ifc.isBound = false;
    mol.interfaceList.push_back(ifc);
    mol.freelist.push_back(0);
    return mol;
}

/*! \brief Hand-built Complex owning a single member molecule. */
Complex gcfr_make_com(int index, const Coord& com, double radius, int memberMolIndex)
{
    Complex cplx {};
    cplx.index = index;
    cplx.comCoord = com;
    cplx.radius = radius;
    cplx.mass = 1.0;
    cplx.isEmpty = false;
    cplx.memberList.push_back(memberMolIndex);
    cplx.numEachMol.assign(MolTemplate::numMolTypes, 0);
    cplx.lastNumberUpdateItrEachMol.assign(MolTemplate::numMolTypes, 0);
    return cplx;
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: initialize_molecule_for_restart()
// -----------------------------------------------------------------------------
void test_gcfr_initialize_molecule_for_restart()
{
    std::cerr << "\n[TEST] test_gcfr_initialize_molecule_for_restart\n"
              << "  Source file: generate_coordinates_for_restart.cpp\n"
              << "  Function:    initialize_molecule_for_restart()\n"
              << "  Scenario:    build one soluble molecule (2 interfaces) inside a\n"
              << "               200 nm cubic water box.\n"
              << "  Criteria:    template data copied over, freelist/interfaceList sized\n"
              << "               to the template, statics incremented, coordinates in box.\n";

    GcfrStaticsGuard guard;      // restores the shared statics on scope exit
    gcfr_prepare_statics(1);     // one molecule type in this test
    gcfr_seed_rng();

    Parameters params {};
    params.numTotalUnits = 0;

    MolTemplate temp = gcfr_make_template("A", /*typeIndex*/ 0, /*copies*/ 1,
        /*nIfaces*/ 2, /*ifaceLen*/ 1.0);
    Membrane mem = gcfr_make_membrane(200.0);

    std::cerr << "  Calling initialize_molecule_for_restart(index = 7)...\n";
    Molecule mol = initialize_molecule_for_restart(7, params, temp, mem);

    // ---- identity / template data --------------------------------------------
    EXPECT_EQ(mol.index, 7) << "the index argument must be stored on the molecule";
    EXPECT_EQ(mol.molTypeIndex, 0) << "molTypeIndex must be copied from the template";
    EXPECT_DOUBLE_EQ(mol.mass, temp.mass) << "mass must be copied from the template";
    EXPECT_FALSE(mol.isLipid) << "isLipid must be copied from the template (false here)";
    EXPECT_FALSE(mol.isPromoter) << "isPromoter must be copied from the template (false here)";
    EXPECT_FALSE(mol.isEmpty) << "a freshly initialised molecule is not an empty slot";

    // ---- interface / freelist bookkeeping -------------------------------------
    ASSERT_EQ(mol.interfaceList.size(), static_cast<size_t>(2))
        << "one Iface per template interface is required for the checks below";
    ASSERT_EQ(mol.freelist.size(), static_cast<size_t>(2))
        << "freelist is sized to the number of interfaces";
    EXPECT_EQ(mol.freelist[0], 0) << "freelist is filled by std::iota starting at 0";
    EXPECT_EQ(mol.freelist[1], 1) << "freelist is filled by std::iota starting at 0";

    for (size_t i = 0; i < mol.interfaceList.size(); ++i) {
        EXPECT_EQ(mol.interfaceList[i].relIndex, static_cast<int>(i))
            << "relative interface index " << i << " should equal its position";
        EXPECT_EQ(mol.interfaceList[i].molTypeIndex, 0)
            << "each interface records its parent molecule type";
        EXPECT_EQ(mol.interfaceList[i].stateIndex, 0)
            << "the default (first) state of the template is used";
        EXPECT_EQ(mol.interfaceList[i].index, temp.interfaceList[i].stateList[0].index)
            << "absolute index comes from the template's first state";
        EXPECT_FALSE(mol.interfaceList[i].isBound)
            << "a newly created molecule has no bound interfaces";

        // The molecule is rigidly rotated, so the |COM - interface| distance
        // must still be the template value (1.0 nm).
        Vector armVec { mol.interfaceList[i].coord - mol.comCoord };
        arm_magnitude_check:
        armVec.calc_magnitude();
        EXPECT_NEAR(armVec.magnitude, 1.0, 1e-9)
            << "rigid rotation must preserve the COM-interface distance";
    }

    // ---- coordinates must be inside the box -----------------------------------
    const double halfBox = 100.0;
    EXPECT_LE(std::abs(mol.comCoord.x), halfBox + 1e-9) << "COM x inside the water box";
    EXPECT_LE(std::abs(mol.comCoord.y), halfBox + 1e-9) << "COM y inside the water box";
    EXPECT_LE(std::abs(mol.comCoord.z), halfBox + 1e-9) << "COM z inside the water box";
    for (const auto& ifc : mol.interfaceList) {
        EXPECT_LE(std::abs(ifc.coord.x), halfBox + 1e-9) << "interface x inside the water box";
        EXPECT_LE(std::abs(ifc.coord.y), halfBox + 1e-9) << "interface y inside the water box";
        EXPECT_LE(std::abs(ifc.coord.z), halfBox + 1e-9) << "interface z inside the water box";
    }
    std::cerr << "  Generated COM = (" << mol.comCoord.x << ", " << mol.comCoord.y
              << ", " << mol.comCoord.z << ")\n";

    // ---- static / parameter side effects --------------------------------------
    EXPECT_EQ(Molecule::numberOfMolecules, 1)
        << "the global molecule counter must be incremented once";
    ASSERT_EQ(MolTemplate::numEachMolType.size(), static_cast<size_t>(1));
    EXPECT_EQ(MolTemplate::numEachMolType[0], 1)
        << "per-type molecule counter must be incremented once";
    EXPECT_EQ(params.numTotalUnits, static_cast<unsigned>(3))
        << "numTotalUnits grows by (#interfaces + 1) = 3";
}

// -----------------------------------------------------------------------------
// Test 2: create_molecule_and_complex_for_restart() appending to empty lists
// -----------------------------------------------------------------------------
void test_gcfr_create_molecule_and_complex_append()
{
    std::cerr << "\n[TEST] test_gcfr_create_molecule_and_complex_append\n"
              << "  Source file: generate_coordinates_for_restart.cpp\n"
              << "  Function:    create_molecule_and_complex_for_restart()\n"
              << "  Scenario:    empty moleculeList/complexList and empty recycle lists,\n"
              << "               so two consecutive calls must append at index 0 and 1.\n"
              << "  Criteria:    lists grow, molecule<->complex cross links are correct,\n"
              << "               counters increment and the template's monomerList grows.\n";

    GcfrStaticsGuard guard;
    gcfr_prepare_statics(1);
    gcfr_seed_rng();

    Parameters params {};
    params.numTotalUnits = 0;

    std::vector<MolTemplate> molTemplateList { gcfr_make_template("A", 0, 1, 1, 1.0) };
    Membrane mem = gcfr_make_membrane(200.0);
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {}; // no reactions -> no overlap can be detected

    std::cerr << "  Calling create_molecule_and_complex_for_restart() twice...\n";
    create_molecule_and_complex_for_restart(molTemplateList[0], params, moleculeList,
        complexList, molTemplateList, forwardRxns, mem);
    create_molecule_and_complex_for_restart(molTemplateList[0], params, moleculeList,
        complexList, molTemplateList, forwardRxns, mem);

    ASSERT_EQ(moleculeList.size(), static_cast<size_t>(2))
        << "each call must append exactly one Molecule";
    ASSERT_EQ(complexList.size(), static_cast<size_t>(2))
        << "each call must append exactly one Complex";

    for (int i = 0; i < 2; ++i) {
        EXPECT_EQ(moleculeList[i].index, i) << "molecule " << i << " stores its own index";
        EXPECT_EQ(moleculeList[i].myComIndex, i)
            << "molecule " << i << " must point at the complex created for it";
        EXPECT_FALSE(moleculeList[i].isEmpty) << "molecule " << i << " is a real molecule";

        EXPECT_EQ(complexList[i].index, i) << "complex " << i << " stores its own index";
        ASSERT_EQ(complexList[i].memberList.size(), static_cast<size_t>(1));
        EXPECT_EQ(complexList[i].memberList[0], i)
            << "complex " << i << " owns exactly the molecule that created it";
        EXPECT_DOUBLE_EQ(complexList[i].radius, molTemplateList[0].radius)
            << "complex radius comes from the MolTemplate";
        EXPECT_DOUBLE_EQ(complexList[i].mass, moleculeList[i].mass)
            << "complex mass comes from the member molecule";
        ASSERT_EQ(complexList[i].numEachMol.size(), static_cast<size_t>(1))
            << "numEachMol is sized by MolTemplate::numMolTypes";
        EXPECT_EQ(complexList[i].numEachMol[0], 1)
            << "the complex contains one molecule of type 0";
    }

    EXPECT_EQ(Molecule::numberOfMolecules, 2) << "two molecules were created";
    EXPECT_EQ(Complex::numberOfComplexes, 2) << "two complexes were created";
    EXPECT_EQ(MolTemplate::numEachMolType[0], 2) << "two molecules of type 0 exist";
    EXPECT_EQ(params.numTotalUnits, static_cast<unsigned>(4))
        << "two molecules x (1 interface + 1 COM) = 4 units";

    // The template is passed by reference here, so the monomer bookkeeping is visible.
    ASSERT_EQ(molTemplateList[0].monomerList.size(), static_cast<size_t>(2))
        << "each created molecule is appended to the template's monomerList";
    EXPECT_EQ(molTemplateList[0].monomerList[0], 0);
    EXPECT_EQ(molTemplateList[0].monomerList[1], 1);
}

// -----------------------------------------------------------------------------
// Test 3: create_molecule_and_complex_for_restart() recycling empty slots
// -----------------------------------------------------------------------------
void test_gcfr_create_molecule_and_complex_recycles_empty_slots()
{
    std::cerr << "\n[TEST] test_gcfr_create_molecule_and_complex_recycles_empty_slots\n"
              << "  Source file: generate_coordinates_for_restart.cpp\n"
              << "  Function:    create_molecule_and_complex_for_restart()\n"
              << "  Scenario:    moleculeList/complexList each contain one live entry and\n"
              << "               one void entry that is registered in the recycle lists.\n"
              << "  Criteria:    the void slot (index 1) is reused, the containers do NOT\n"
              << "               grow, and the recycle lists are emptied.\n";

    GcfrStaticsGuard guard;
    gcfr_prepare_statics(1);
    gcfr_seed_rng();

    Parameters params {};
    params.numTotalUnits = 0;

    std::vector<MolTemplate> molTemplateList { gcfr_make_template("A", 0, 1, 1, 1.0) };
    Membrane mem = gcfr_make_membrane(200.0);
    std::vector<ForwardRxn> forwardRxns {};

    // --- one live molecule/complex pair, parked far away in the corner of the box
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    moleculeList.push_back(gcfr_make_mol(0, 0, 0, Coord(90.0, 90.0, 90.0),
        Coord(91.0, 90.0, 90.0), 0));
    complexList.push_back(gcfr_make_com(0, Coord(90.0, 90.0, 90.0), 1.0, 0));

    // --- one void molecule/complex pair that must be recycled
    Molecule voidMol {};
    voidMol.index = 1;
    voidMol.isEmpty = true;
    moleculeList.push_back(voidMol);

    Complex voidCom {};
    voidCom.index = 1;
    voidCom.isEmpty = true;
    complexList.push_back(voidCom);

    Molecule::emptyMolList.push_back(1);
    Complex::emptyComList.push_back(1);
    Molecule::numberOfMolecules = 1; // one live molecule already in the system
    Complex::numberOfComplexes = 1;

    std::cerr << "  Calling create_molecule_and_complex_for_restart()...\n";
    create_molecule_and_complex_for_restart(molTemplateList[0], params, moleculeList,
        complexList, molTemplateList, forwardRxns, mem);

    EXPECT_EQ(moleculeList.size(), static_cast<size_t>(2))
        << "an existing void slot must be reused instead of appending";
    EXPECT_EQ(complexList.size(), static_cast<size_t>(2))
        << "an existing void slot must be reused instead of appending";

    EXPECT_FALSE(moleculeList[1].isEmpty) << "slot 1 now holds a real molecule";
    EXPECT_EQ(moleculeList[1].index, 1) << "the recycled molecule keeps slot index 1";
    EXPECT_EQ(moleculeList[1].myComIndex, 1) << "and points at the recycled complex slot";
    ASSERT_EQ(complexList[1].memberList.size(), static_cast<size_t>(1));
    EXPECT_EQ(complexList[1].memberList[0], 1) << "the recycled complex owns molecule 1";

    EXPECT_TRUE(Molecule::emptyMolList.empty())
        << "the consumed molecule slot is popped from emptyMolList";
    EXPECT_TRUE(Complex::emptyComList.empty())
        << "the consumed complex slot is popped from emptyComList";

    EXPECT_EQ(Molecule::numberOfMolecules, 2) << "one more live molecule";
    EXPECT_EQ(Complex::numberOfComplexes, 2) << "one more live complex";
}

// -----------------------------------------------------------------------------
// Test 4: generate_coordinates_for_restart()
// -----------------------------------------------------------------------------
void test_gcfr_generate_coordinates_for_restart()
{
    std::cerr << "\n[TEST] test_gcfr_generate_coordinates_for_restart\n"
              << "  Source file: generate_coordinates_for_restart.cpp\n"
              << "  Function:    generate_coordinates_for_restart()\n"
              << "  Scenario:    two templates, numMolTemplateBeforeAdd = 1, so only the\n"
              << "               second template (3 copies) may be instantiated.\n"
              << "  Criteria:    exactly 3 molecules/complexes of type 1 appear, none of\n"
              << "               type 0, and (a known quirk) the template's monomerList is\n"
              << "               NOT updated because the loop works on a local copy.\n";

    GcfrStaticsGuard guard;
    gcfr_prepare_statics(2); // two molecule types
    gcfr_seed_rng();

    Parameters params {};
    params.numTotalUnits = 0;

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gcfr_make_template("A", /*typeIndex*/ 0, /*copies*/ 5,
        /*nIfaces*/ 1, /*ifaceLen*/ 1.0)); // pre-existing type, must be skipped
    molTemplateList.push_back(gcfr_make_template("B", /*typeIndex*/ 1, /*copies*/ 3,
        /*nIfaces*/ 2, /*ifaceLen*/ 1.0)); // newly added type

    Membrane mem = gcfr_make_membrane(200.0);
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};

    std::cerr << "  Calling generate_coordinates_for_restart(numMolTemplateBeforeAdd = 1)...\n";
    generate_coordinates_for_restart(params, moleculeList, complexList, molTemplateList,
        forwardRxns, mem, /*numMolTemplateBeforeAdd*/ 1, /*numForwardRxnBdeforeAdd*/ 0);

    ASSERT_EQ(moleculeList.size(), static_cast<size_t>(3))
        << "only the 3 copies of template index 1 are created";
    ASSERT_EQ(complexList.size(), static_cast<size_t>(3))
        << "each new molecule gets its own complex";

    for (size_t i = 0; i < moleculeList.size(); ++i) {
        EXPECT_EQ(moleculeList[i].molTypeIndex, 1)
            << "molecule " << i << " must be of the newly added type";
        EXPECT_EQ(moleculeList[i].myComIndex, static_cast<int>(i))
            << "molecule " << i << " links to complex " << i;
        EXPECT_EQ(moleculeList[i].interfaceList.size(), static_cast<size_t>(2))
            << "template B declares two interfaces";
        EXPECT_FALSE(moleculeList[i].isEmpty);
    }

    EXPECT_EQ(MolTemplate::numEachMolType[0], 0)
        << "no copies of the pre-existing type may be added";
    EXPECT_EQ(MolTemplate::numEachMolType[1], 3)
        << "three copies of the newly added type were requested";
    EXPECT_EQ(Molecule::numberOfMolecules, 3);
    EXPECT_EQ(Complex::numberOfComplexes, 3);
    EXPECT_EQ(params.numTotalUnits, static_cast<unsigned>(9))
        << "3 molecules x (2 interfaces + 1 COM) = 9 units";

    // Documented behaviour: generate_coordinates_for_restart() copies the template
    // into a local `oneTemp`, so the monomerList bookkeeping performed inside
    // create_molecule_and_complex_for_restart() is discarded.
    EXPECT_TRUE(molTemplateList[1].monomerList.empty())
        << "monomerList of the stored template stays empty (local copy is modified)";
    EXPECT_TRUE(molTemplateList[0].monomerList.empty())
        << "the skipped template is untouched";
}

// -----------------------------------------------------------------------------
// Test 5: moleculeOverlapsForRestart() - molecule compared only with itself
// -----------------------------------------------------------------------------
void test_gcfr_overlap_self_only()
{
    std::cerr << "\n[TEST] test_gcfr_overlap_self_only\n"
              << "  Source file: generate_coordinates_for_restart.cpp\n"
              << "  Function:    moleculeOverlapsForRestart()\n"
              << "  Scenario:    the only molecule in the system is the created one, so the\n"
              << "               scan loop `continue`s immediately.\n"
              << "  Criteria:    returns false (no overlap, no resampling needed).\n";

    GcfrStaticsGuard guard;
    gcfr_prepare_statics(1);

    Parameters params {};
    Membrane mem = gcfr_make_membrane(200.0);
    std::vector<MolTemplate> molTemplateList { gcfr_make_template("A", 0, 1, 1, 1.0) };
    std::vector<ForwardRxn> forwardRxns {};

    std::vector<Molecule> moleculeList { gcfr_make_mol(0, 0, 0, Coord(0.0, 0.0, 0.0),
        Coord(1.0, 0.0, 0.0), 0) };
    std::vector<Complex> complexList { gcfr_make_com(0, Coord(0.0, 0.0, 0.0), 2.0, 0) };

    std::cerr << "  Calling moleculeOverlapsForRestart()...\n";
    bool overlaps = moleculeOverlapsForRestart(params, moleculeList[0], moleculeList,
        complexList, forwardRxns, molTemplateList, mem);

    EXPECT_FALSE(overlaps) << "a molecule is never considered to overlap with itself";
}

// -----------------------------------------------------------------------------
// Test 6: moleculeOverlapsForRestart() - bounding spheres do not touch
// -----------------------------------------------------------------------------
void test_gcfr_overlap_bounding_sphere_miss()
{
    std::cerr << "\n[TEST] test_gcfr_overlap_bounding_sphere_miss\n"
              << "  Source file: generate_coordinates_for_restart.cpp\n"
              << "  Function:    moleculeOverlapsForRestart()\n"
              << "  Scenario:    the partner complex is 50 nm away while the bounding radii\n"
              << "               only sum to 4 nm.\n"
              << "  Criteria:    returns false.  NOTE: the implementation *returns* false\n"
              << "               here rather than continuing the scan, so the very first\n"
              << "               distant partner short-circuits the whole search.\n";

    GcfrStaticsGuard guard;
    gcfr_prepare_statics(2);

    Parameters params {};
    Membrane mem = gcfr_make_membrane(200.0);
    std::vector<MolTemplate> molTemplateList {
        gcfr_make_template("A", 0, 1, 1, 1.0),
        gcfr_make_template("B", 1, 1, 1, 1.0)
    };
    std::vector<ForwardRxn> forwardRxns {};

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(gcfr_make_mol(0, 0, 0, Coord(0.0, 0.0, 0.0), Coord(1.0, 0.0, 0.0), 0));
    moleculeList.push_back(gcfr_make_mol(1, 1, 1, Coord(50.0, 0.0, 0.0), Coord(51.0, 0.0, 0.0), 1));

    std::vector<Complex> complexList;
    complexList.push_back(gcfr_make_com(0, Coord(0.0, 0.0, 0.0), 2.0, 0));
    complexList.push_back(gcfr_make_com(1, Coord(50.0, 0.0, 0.0), 2.0, 1));

    std::cerr << "  Calling moleculeOverlapsForRestart()...\n";
    bool overlaps = moleculeOverlapsForRestart(params, moleculeList[0], moleculeList,
        complexList, forwardRxns, molTemplateList, mem);

    EXPECT_FALSE(overlaps) << "well separated bounding spheres cannot overlap";
}

// -----------------------------------------------------------------------------
// Test 7/8: moleculeOverlapsForRestart() with a matching bimolecular reaction
//
// The two molecules are inside each other's bounding spheres and their reacting
// interfaces match reactantListNew[0] / reactantListNew[1] of a ForwardRxn.  The
// implementation then reports "overlap" when the interface separation is LARGER
// than the binding radius (and reports no overlap when it is smaller) - the tests
// below simply pin down that behaviour as written.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Builds the shared fixture used by the two reaction-driven overlap tests.
 *
 * \param[in] partnerIfaceX x coordinate of the partner molecule's interface; the
 *            created molecule's interface always sits at x = 0.5.
 */
bool gcfr_run_reaction_overlap(double partnerIfaceX, std::vector<MolTemplate>& molTemplateList)
{
    Parameters params {};
    Membrane mem = gcfr_make_membrane(200.0);

    // A single bimolecular reaction with a 1 nm binding radius.
    ForwardRxn rxn {};
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = 1.0;
    // reactant 0: molecule type 0, relative/absolute interface index 0, free, no state
    rxn.reactantListNew.emplace_back("i0", 0, 0, 0, '\0', false);
    // reactant 1: molecule type 1, absolute interface index 1, relative index 0
    rxn.reactantListNew.emplace_back("i0", 1, 1, 0, '\0', false);
    std::vector<ForwardRxn> forwardRxns { rxn };

    // The two molecules are 2 nm apart; bounding radii sum to 4 nm, so the
    // bounding-sphere pre-filter is passed and the reaction loop is reached.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(gcfr_make_mol(0, 0, /*type*/ 0, Coord(0.0, 0.0, 0.0),
        Coord(0.5, 0.0, 0.0), /*absIfaceIndex*/ 0));
    moleculeList.push_back(gcfr_make_mol(1, 1, /*type*/ 1, Coord(2.0, 0.0, 0.0),
        Coord(partnerIfaceX, 0.0, 0.0), /*absIfaceIndex*/ 1));

    std::vector<Complex> complexList;
    complexList.push_back(gcfr_make_com(0, Coord(0.0, 0.0, 0.0), 2.0, 0));
    complexList.push_back(gcfr_make_com(1, Coord(2.0, 0.0, 0.0), 2.0, 1));

    return moleculeOverlapsForRestart(params, moleculeList[0], moleculeList, complexList,
        forwardRxns, molTemplateList, mem);
}

} // anonymous namespace

void test_gcfr_overlap_reactive_pair_far_interfaces()
{
    std::cerr << "\n[TEST] test_gcfr_overlap_reactive_pair_far_interfaces\n"
              << "  Source file: generate_coordinates_for_restart.cpp\n"
              << "  Function:    moleculeOverlapsForRestart()\n"
              << "  Scenario:    bounding spheres intersect, the interfaces are recognised\n"
              << "               as reactants of a ForwardRxn (bindRadius = 1 nm) and are\n"
              << "               separated by 2 nm.\n"
              << "  Criteria:    returns true, because the code flags the pair when the\n"
              << "               interface separation exceeds the binding radius.\n";

    GcfrStaticsGuard guard;
    gcfr_prepare_statics(2);

    std::vector<MolTemplate> molTemplateList {
        gcfr_make_template("A", 0, 1, 1, 0.5),
        gcfr_make_template("B", 1, 1, 1, 0.5)
    };

    // partner interface at x = 2.5 -> separation from x = 0.5 is 2.0 nm > bindRadius
    std::cerr << "  Interface separation = 2.0 nm, bindRadius = 1.0 nm\n";
    bool overlaps = gcfr_run_reaction_overlap(2.5, molTemplateList);

    EXPECT_TRUE(overlaps)
        << "separation (2 nm) > bindRadius (1 nm) makes the function report an overlap";
}

void test_gcfr_overlap_reactive_pair_close_interfaces()
{
    std::cerr << "\n[TEST] test_gcfr_overlap_reactive_pair_close_interfaces\n"
              << "  Source file: generate_coordinates_for_restart.cpp\n"
              << "  Function:    moleculeOverlapsForRestart()\n"
              << "  Scenario:    identical to the previous test but the two reacting\n"
              << "               interfaces are only 0.4 nm apart.\n"
              << "  Criteria:    returns false, i.e. separations below the binding radius\n"
              << "               are accepted and the scan simply runs to completion.\n";

    GcfrStaticsGuard guard;
    gcfr_prepare_statics(2);

    std::vector<MolTemplate> molTemplateList {
        gcfr_make_template("A", 0, 1, 1, 0.5),
        gcfr_make_template("B", 1, 1, 1, 0.5)
    };

    // partner interface at x = 0.9 -> separation from x = 0.5 is 0.4 nm < bindRadius
    std::cerr << "  Interface separation = 0.4 nm, bindRadius = 1.0 nm\n";
    bool overlaps = gcfr_run_reaction_overlap(0.9, molTemplateList);

    EXPECT_FALSE(overlaps)
        << "separation (0.4 nm) < bindRadius (1 nm) is not reported as an overlap";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Every check above uses non-fatal EXPECT_* (or ASSERT_*
// only where continuing would dereference out of range), so all tests run even
// if some of them fail.
// -----------------------------------------------------------------------------
TEST(GenerateCoordinatesForRestart, InitializeMoleculeForRestart)
{
    test_gcfr_initialize_molecule_for_restart();
}

TEST(GenerateCoordinatesForRestart, CreateMoleculeAndComplexAppend)
{
    test_gcfr_create_molecule_and_complex_append();
}

TEST(GenerateCoordinatesForRestart, CreateMoleculeAndComplexRecyclesEmptySlots)
{
    test_gcfr_create_molecule_and_complex_recycles_empty_slots();
}

TEST(GenerateCoordinatesForRestart, GenerateCoordinatesForRestart)
{
    test_gcfr_generate_coordinates_for_restart();
}

TEST(GenerateCoordinatesForRestart, OverlapSelfOnly) { test_gcfr_overlap_self_only(); }

TEST(GenerateCoordinatesForRestart, OverlapBoundingSphereMiss)
{
    test_gcfr_overlap_bounding_sphere_miss();
}

TEST(GenerateCoordinatesForRestart, OverlapReactivePairFarInterfaces)
{
    test_gcfr_overlap_reactive_pair_far_interfaces();
}

TEST(GenerateCoordinatesForRestart, OverlapReactivePairCloseInterfaces)
{
    test_gcfr_overlap_reactive_pair_close_interfaces();
}