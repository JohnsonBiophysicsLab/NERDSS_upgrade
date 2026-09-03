/*! \file test_generate_coordinates.cpp
 *
 * ### Unit tests for src/system_setup/generate_coordinates.cpp
 *
 * The translation unit under test contains four functions:
 *
 *   1. generate_coordinates()        - top level driver: creates every molecule
 *                                      and complex from the MolTemplate copy
 *                                      numbers, removes overlaps, optionally
 *                                      applies coordinates read from a PDB file
 *                                      and finally dumps an XYZ file.
 *   2. updateMoleculeCoordinates()   - reads "ATOM ... COM ..." records out of a
 *                                      PDB file and rigidly translates the first
 *                                      N molecules of each type onto them.
 *   3. checkAndResolveOverlap()      - decides whether two molecules whose
 *                                      interfaces can react are closer than the
 *                                      binding radius and, if so, either warns
 *                                      (PDB supplied coordinates) or relocates
 *                                      the second molecule at random.
 *   4. fixOverlappingMolecules()     - iterates checkAndResolveOverlap() over
 *                                      either every pair or only the pairs that
 *                                      involve a molecule read from the PDB.
 *
 * The helper routines are ordinary (non-static) free functions, so they are
 * re-declared below and exercised directly.  Every test prints what it is doing
 * and what the pass criterion is, so a reader of the console log can follow
 * along.
 *
 * NOTE: the implicit-lipid branch of generate_coordinates() is intentionally not
 * exercised here; it mutates global lipid bookkeeping shared with the rest of
 * the suite, so it is covered by the implicit lipid tests instead.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "system_setup/system_setup.hpp"

// -----------------------------------------------------------------------------
// Declarations of the helper routines that live in generate_coordinates.cpp.
// The default arguments are deliberately omitted; every call below passes all
// arguments explicitly.
// -----------------------------------------------------------------------------
int fixOverlappingMolecules(std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject,
    const std::vector<int>& changedMoleculeIndices, bool checkOnlyChanged);

std::vector<int> updateMoleculeCoordinates(std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList, const std::string& coordinateFileName);

bool checkAndResolveOverlap(const Molecule& mol1, Molecule& mol2,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const Membrane& membraneObject, std::vector<Complex>& complexList, bool mol2Changed);

namespace {

// -----------------------------------------------------------------------------
// Small fixtures / helpers, all prefixed with "gc_" so they cannot collide with
// helpers defined by any other test file in the suite.
// -----------------------------------------------------------------------------

/*! \brief Allocate (once) and deterministically seed the global GSL generator.
 *
 * `r` is defined in gtest_main.cpp; we only ever declare it extern (via
 * math/rand_gsl.hpp).  Seeding with a fixed value keeps the random molecule
 * placement reproducible from run to run.
 */
void gc_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, 42);
}

/*! \brief Reset every static counter the code under test relies on.
 *
 * generate_coordinates() feeds Complex::numberOfComplexes to
 * initialize_molecule(), and Complex's constructor indexes
 * MolTemplate::numMolTypes-sized vectors, so these have to be consistent before
 * anything is built.
 */
void gc_reset_globals(int numMolTypes)
{
    Molecule::numberOfMolecules = 0;
    Molecule::maxID = 0;
    Molecule::emptyMolList.clear();

    Complex::numberOfComplexes = 0;
    Complex::maxID = 0;
    Complex::emptyComList.clear();
    Complex::obs.clear();

    MolTemplate::numMolTypes = static_cast<unsigned>(numMolTypes);
    MolTemplate::numEachMolType.assign(numMolTypes, 0);
    MolTemplate::absToRelIface.clear();
}

/*! \brief Build a cubic, non-spherical membrane centred on the origin.
 *
 * The WaterBox(vector) constructor is used on purpose: it is the only one that
 * fills in xLeft/xRight, which Molecule::create_random_coords() reads when it
 * samples a new x coordinate.
 */
Membrane gc_make_membrane(double side)
{
    Membrane membrane;
    membrane.waterBox = Membrane::WaterBox(std::vector<double> { side, side, side });
    membrane.isSphere = false;
    membrane.isBox = true;
    membrane.implicitLipid = false;
    membrane.hasCompartment = false;
    membrane.sphereR = 0.0;
    return membrane;
}

/*! \brief Create a simple MolTemplate with `numIfaces` single-state interfaces.
 *
 * \param absStateIndex running counter used to hand out unique absolute state
 *        indices across all templates (incremented in place).
 */
MolTemplate gc_make_template(
    const std::string& name, int typeIndex, int numIfaces, double radius, int& absStateIndex)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.mass = 1.0;
    temp.radius = radius;
    temp.copies = 0;
    temp.comCoord = Coord { 0.0, 0.0, 0.0 };
    temp.D = Coord { 10.0, 10.0, 10.0 };
    temp.Dr = Coord { 0.01, 0.01, 0.01 };
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPoint = false;
    temp.isRod = false;
    temp.isPromoter = false;

    for (int i = 0; i < numIfaces; ++i) {
        Interface iface;
        iface.index = i;
        iface.name = name + "_i" + std::to_string(i);
        // Interfaces sit one nanometre away from the centre of mass.
        iface.iCoord = Coord { (i == 0) ? 1.0 : 0.0, (i == 1) ? 1.0 : 0.0, 0.0 };

        Interface::State state;
        state.index = absStateIndex++;
        state.iden = '\0';
        state.ifaceAndStateName = iface.name;
        iface.stateList.push_back(state);

        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Hand-build a fully initialised Molecule that matches a MolTemplate. */
Molecule gc_make_molecule(int index, int comIndex, const MolTemplate& temp, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.complexId = comIndex;
    mol.molTypeIndex = temp.molTypeIndex;
    mol.mass = temp.mass;
    mol.isLipid = temp.isLipid;
    mol.isImplicitLipid = temp.isImplicitLipid;
    mol.isEmpty = false;
    mol.comCoord = com;

    for (size_t i = 0; i < temp.interfaceList.size(); ++i) {
        Molecule::Iface iface;
        iface.coord = com + temp.interfaceList[i].iCoord;
        iface.relIndex = static_cast<int>(i);
        iface.index = temp.interfaceList[i].stateList[0].index;
        iface.stateIndex = 0;
        iface.stateIden = temp.interfaceList[i].stateList[0].iden;
        iface.molTypeIndex = temp.molTypeIndex;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Hand-build a single-member Complex for the given Molecule. */
Complex gc_make_complex(int index, const Molecule& mol)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.comCoord = mol.comCoord;
    com.mass = mol.mass;
    com.radius = 1.0;
    com.D = Coord { 10.0, 10.0, 10.0 };
    com.Dr = Coord { 0.01, 0.01, 0.01 };
    com.memberList.push_back(mol.index);
    com.numEachMol.assign(MolTemplate::numMolTypes, 0);
    if (mol.molTypeIndex >= 0 && mol.molTypeIndex < static_cast<int>(MolTemplate::numMolTypes))
        ++com.numEachMol[mol.molTypeIndex];
    com.lastNumberUpdateItrEachMol.assign(MolTemplate::numMolTypes, 0);
    return com;
}

/*! \brief Register a bimolecular reaction between two molecule types.
 *
 * The reaction index is pushed onto the *first* interface of `templateA` so that
 * checkAndResolveOverlap() can find it while scanning mol1's forward reactions.
 */
void gc_add_reaction(std::vector<ForwardRxn>& forwardRxns, std::vector<MolTemplate>& molTemplateList,
    int typeA, int typeB, double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.reactantListNew.emplace_back(molTemplateList[typeA].interfaceList[0].name, typeA,
        molTemplateList[typeA].interfaceList[0].stateList[0].index, 0, '\0', false);
    rxn.reactantListNew.emplace_back(molTemplateList[typeB].interfaceList[0].name, typeB,
        molTemplateList[typeB].interfaceList[0].stateList[0].index, 0, '\0', false);
    forwardRxns.push_back(rxn);

    const unsigned rxnIndex = static_cast<unsigned>(forwardRxns.size() - 1);
    molTemplateList[typeA].interfaceList[0].stateList[0].myForwardRxns.push_back(rxnIndex);
    molTemplateList[typeB].interfaceList[0].stateList[0].myForwardRxns.push_back(rxnIndex);
}

/*! \brief Write one PDB "ATOM" record with the exact column layout the parser
 *         in updateMoleculeCoordinates() expects.
 *
 * substr(12,4) -> atom name, substr(17,4) -> residue (molecule) name,
 * substr(30,8)/substr(38,8)/substr(46,8) -> x/y/z.
 */
void gc_write_pdb_line(std::ofstream& os, int serial, const std::string& atomName,
    const std::string& resName, int resSeq, double x, double y, double z)
{
    char buf[128];
    std::snprintf(buf, sizeof(buf), "ATOM  %5d %-4s %-4s%4d     %8.3f%8.3f%8.3f  1.00  0.00", serial,
        atomName.c_str(), resName.c_str(), resSeq, x, y, z);
    os << buf << '\n';
}

/*! \brief Squared distance between two coordinates (matches the code's metric). */
double gc_sq_dist(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

} // anonymous namespace

// =============================================================================
// updateMoleculeCoordinates()
// =============================================================================

/*! Reads two COM records for molecule type "A" and checks that exactly the first
 *  two A molecules are rigidly translated onto them. */
void gc_test_update_coordinates_from_pdb()
{
    std::cerr << "\n[TEST] gc_test_update_coordinates_from_pdb\n"
              << "  Source file : generate_coordinates.cpp\n"
              << "  Function    : updateMoleculeCoordinates()\n"
              << "  Scenario    : PDB supplies 2 COM records for type A while the\n"
              << "                system holds 3 A molecules and 1 B molecule.\n"
              << "  Pass criteria: only molecules 0 and 1 are reported as changed,\n"
              << "                 their COM lands exactly on the PDB value and the\n"
              << "                 interface offsets are rigidly carried along.\n";

    gc_init_rng();
    gc_reset_globals(2);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 2.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 2.0, absState));

    // Three A molecules and one B molecule at known starting positions.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(gc_make_molecule(1, 1, molTemplateList[0], Coord { 5.0, 5.0, 5.0 }));
    moleculeList.push_back(gc_make_molecule(2, 2, molTemplateList[0], Coord { -7.0, 1.0, 2.0 }));
    moleculeList.push_back(gc_make_molecule(3, 3, molTemplateList[1], Coord { 9.0, 9.0, 9.0 }));

    // Remember the untouched state so we can verify what must NOT move.
    const Coord thirdA = moleculeList[2].comCoord;
    const Coord onlyB = moleculeList[3].comCoord;

    // Build the PDB file.  A stray non-COM atom record and a comment line are
    // included to prove they are ignored.
    const std::string fileName = "test_gc_update_coords.pdb";
    {
        std::ofstream out(fileName);
        out << "REMARK generated by unit test\n";
        gc_write_pdb_line(out, 1, "COM", "A", 1, 10.0, -20.0, 5.0);
        gc_write_pdb_line(out, 2, "A_i0", "A", 1, 11.0, -20.0, 5.0); // must be ignored
        gc_write_pdb_line(out, 3, "COM", "A", 2, -30.0, 15.0, -5.0);
        out.close();
    }

    std::cerr << "  Calling updateMoleculeCoordinates()...\n";
    std::vector<int> changed = updateMoleculeCoordinates(moleculeList, molTemplateList, fileName);

    // Exactly two molecules should have been reported as changed.
    EXPECT_EQ(changed.size(), 2u) << "Only the first two A molecules should be updated";
    if (changed.size() == 2u) {
        EXPECT_EQ(changed[0], 0) << "First updated molecule should be index 0";
        EXPECT_EQ(changed[1], 1) << "Second updated molecule should be index 1";
    }

    // COM coordinates must match the PDB values exactly.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 10.0);
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, -20.0);
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 5.0);
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, -30.0);
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.y, 15.0);
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.z, -5.0);

    // The translation is rigid: the interface keeps its offset of (1,0,0).
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x - moleculeList[0].comCoord.x, 1.0);
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.y - moleculeList[0].comCoord.y, 0.0);
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.z - moleculeList[0].comCoord.z, 0.0);

    // The third A molecule and the B molecule must be untouched.
    EXPECT_DOUBLE_EQ(moleculeList[2].comCoord.x, thirdA.x) << "Extra A molecule must not move";
    EXPECT_DOUBLE_EQ(moleculeList[3].comCoord.x, onlyB.x) << "B molecule must not move";

    std::cerr << "  Molecule 0 COM = (" << moleculeList[0].comCoord.x << ", "
              << moleculeList[0].comCoord.y << ", " << moleculeList[0].comCoord.z << ")\n";

    std::remove(fileName.c_str());
}

/*! A missing file must be reported and produce an empty "changed" list. */
void gc_test_update_coordinates_missing_file()
{
    std::cerr << "\n[TEST] gc_test_update_coordinates_missing_file\n"
              << "  Function    : updateMoleculeCoordinates()\n"
              << "  Scenario    : the coordinate file does not exist.\n"
              << "  Pass criteria: an empty index vector is returned and no molecule\n"
              << "                 is modified (no crash).\n";

    gc_reset_globals(1);
    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 2.0, absState));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(gc_make_molecule(0, 0, molTemplateList[0], Coord { 3.0, 3.0, 3.0 }));

    std::vector<int> changed = updateMoleculeCoordinates(
        moleculeList, molTemplateList, "this_file_should_not_exist_gc_test.pdb");

    EXPECT_TRUE(changed.empty()) << "No molecules can change when the file cannot be opened";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 3.0) << "Coordinates must be untouched";
}

/*! Coordinates for an unknown molecule type are skipped; surplus coordinates for
 *  a known type beyond the number of molecules present are ignored. */
void gc_test_update_coordinates_unknown_and_surplus()
{
    std::cerr << "\n[TEST] gc_test_update_coordinates_unknown_and_surplus\n"
              << "  Function    : updateMoleculeCoordinates()\n"
              << "  Scenario    : PDB names a molecule type absent from the system and\n"
              << "                supplies more A records than there are A molecules.\n"
              << "  Pass criteria: unknown records are skipped, and only min(#records,\n"
              << "                 #molecules) molecules are updated.\n";

    gc_init_rng();
    gc_reset_globals(1);
    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 2.0, absState));

    // Only a single A molecule exists.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 }));

    const std::string fileName = "test_gc_unknown_coords.pdb";
    {
        std::ofstream out(fileName);
        gc_write_pdb_line(out, 1, "COM", "ZZZ", 1, 1.0, 1.0, 1.0);   // unknown type
        gc_write_pdb_line(out, 2, "COM", "A", 1, 4.0, -4.0, 8.0);    // used
        gc_write_pdb_line(out, 3, "COM", "A", 2, 40.0, 40.0, 40.0);  // surplus, ignored
        out.close();
    }

    std::vector<int> changed = updateMoleculeCoordinates(moleculeList, molTemplateList, fileName);

    EXPECT_EQ(changed.size(), 1u) << "Only one A molecule exists, so only one can be updated";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 4.0) << "First A record should be applied";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, -4.0);
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 8.0);

    std::remove(fileName.c_str());
}

// =============================================================================
// checkAndResolveOverlap()
// =============================================================================

/*! Without a registered reaction between the two interfaces no overlap can be
 *  detected, even when the interfaces are coincident. */
void gc_test_overlap_no_reaction()
{
    std::cerr << "\n[TEST] gc_test_overlap_no_reaction\n"
              << "  Function    : checkAndResolveOverlap()\n"
              << "  Scenario    : two molecules sit on top of each other but no\n"
              << "                forward reaction connects their interfaces.\n"
              << "  Pass criteria: returns false and leaves the second molecule where\n"
              << "                 it was (overlap is only defined for reactive pairs).\n";

    gc_init_rng();
    gc_reset_globals(2);
    Membrane membrane = gc_make_membrane(100.0);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 3.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 3.0, absState));

    std::vector<ForwardRxn> forwardRxns; // deliberately empty

    Molecule mol1 = gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 });
    Molecule mol2 = gc_make_molecule(1, 1, molTemplateList[1], Coord { 0.0, 0.0, 0.0 });

    std::vector<Complex> complexList;
    complexList.push_back(gc_make_complex(0, mol1));
    complexList.push_back(gc_make_complex(1, mol2));

    const Coord before = mol2.comCoord;

    bool overlap = checkAndResolveOverlap(
        mol1, mol2, molTemplateList, forwardRxns, membrane, complexList, false);

    EXPECT_FALSE(overlap) << "No reaction => no overlap can be flagged";
    EXPECT_DOUBLE_EQ(mol2.comCoord.x, before.x) << "Molecule must not be relocated";
    EXPECT_DOUBLE_EQ(mol2.comCoord.y, before.y);
    EXPECT_DOUBLE_EQ(mol2.comCoord.z, before.z);
}

/*! Molecules that are not "in vicinity" are rejected before any interface pair
 *  distance is computed. */
void gc_test_overlap_far_apart()
{
    std::cerr << "\n[TEST] gc_test_overlap_far_apart\n"
              << "  Function    : checkAndResolveOverlap()\n"
              << "  Scenario    : a reactive pair whose centres of mass are far apart\n"
              << "                (well beyond radius1 + radius2 + 2).\n"
              << "  Pass criteria: the vicinity short-circuit returns false and the\n"
              << "                 second molecule is untouched.\n";

    gc_init_rng();
    gc_reset_globals(2);
    Membrane membrane = gc_make_membrane(100.0);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 3.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 3.0, absState));

    std::vector<ForwardRxn> forwardRxns;
    gc_add_reaction(forwardRxns, molTemplateList, 0, 1, 5.0);

    Molecule mol1 = gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 });
    Molecule mol2 = gc_make_molecule(1, 1, molTemplateList[1], Coord { 40.0, 0.0, 0.0 });

    std::vector<Complex> complexList;
    complexList.push_back(gc_make_complex(0, mol1));
    complexList.push_back(gc_make_complex(1, mol2));

    bool overlap = checkAndResolveOverlap(
        mol1, mol2, molTemplateList, forwardRxns, membrane, complexList, false);

    EXPECT_FALSE(overlap) << "Molecules 40 nm apart cannot overlap";
    EXPECT_DOUBLE_EQ(mol2.comCoord.x, 40.0) << "Distant molecule must not be relocated";
}

/*! A genuine overlap for a randomly generated configuration must be reported and
 *  the second molecule must be given brand new random coordinates whose complex
 *  centre of mass is kept in sync. */
void gc_test_overlap_relocates_molecule()
{
    std::cerr << "\n[TEST] gc_test_overlap_relocates_molecule\n"
              << "  Function    : checkAndResolveOverlap()\n"
              << "  Scenario    : reactive interfaces 2 nm apart with a binding radius\n"
              << "                of 5 nm, mol2Changed == false.\n"
              << "  Pass criteria: returns true, mol2 receives new coordinates, its\n"
              << "                 complex COM is synchronised, and the rigid\n"
              << "                 interface-to-COM distance (1 nm) is preserved.\n";

    gc_init_rng();
    gc_reset_globals(2);
    Membrane membrane = gc_make_membrane(100.0);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 3.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 3.0, absState));

    std::vector<ForwardRxn> forwardRxns;
    gc_add_reaction(forwardRxns, molTemplateList, 0, 1, 5.0);

    // Interfaces end up at (1,0,0) and (3,0,0): squared separation 4 < 25.
    Molecule mol1 = gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 });
    Molecule mol2 = gc_make_molecule(1, 1, molTemplateList[1], Coord { 2.0, 0.0, 0.0 });

    std::vector<Complex> complexList;
    complexList.push_back(gc_make_complex(0, mol1));
    complexList.push_back(gc_make_complex(1, mol2));

    const Coord before = mol2.comCoord;

    bool overlap = checkAndResolveOverlap(
        mol1, mol2, molTemplateList, forwardRxns, membrane, complexList, false);

    EXPECT_TRUE(overlap) << "Interfaces closer than the binding radius must overlap";

    // The molecule was re-randomised somewhere in the 100 nm box.
    const bool moved = (mol2.comCoord.x != before.x) || (mol2.comCoord.y != before.y)
        || (mol2.comCoord.z != before.z);
    EXPECT_TRUE(moved) << "Overlapping molecule should receive new random coordinates";

    // The parent complex must have been updated to the new position.
    EXPECT_DOUBLE_EQ(complexList[mol2.myComIndex].comCoord.x, mol2.comCoord.x)
        << "Complex COM must track the relocated molecule";
    EXPECT_DOUBLE_EQ(complexList[mol2.myComIndex].comCoord.y, mol2.comCoord.y);
    EXPECT_DOUBLE_EQ(complexList[mol2.myComIndex].comCoord.z, mol2.comCoord.z);

    // The molecule stays rigid: |iface - COM| is still the template value 1 nm.
    const double offset = std::sqrt(gc_sq_dist(mol2.interfaceList[0].coord, mol2.comCoord));
    EXPECT_NEAR(offset, 1.0, 1e-9) << "Random rotation must preserve the interface offset";

    // The new position stays inside the simulation box.
    EXPECT_LE(std::abs(mol2.comCoord.x), 50.0);
    EXPECT_LE(std::abs(mol2.comCoord.y), 50.0);
    EXPECT_LE(std::abs(mol2.comCoord.z), 50.0);

    std::cerr << "  Relocated COM = (" << mol2.comCoord.x << ", " << mol2.comCoord.y << ", "
              << mol2.comCoord.z << ")\n";
}

/*! When the second molecule's coordinates came from the user's PDB file the
 *  routine only prints a warning: it returns false and moves nothing. */
void gc_test_overlap_warns_for_pdb_molecule()
{
    std::cerr << "\n[TEST] gc_test_overlap_warns_for_pdb_molecule\n"
              << "  Function    : checkAndResolveOverlap()\n"
              << "  Scenario    : identical overlap as before but mol2Changed == true\n"
              << "                (coordinates supplied by the user).\n"
              << "  Pass criteria: returns false (no overlap fix), only a warning is\n"
              << "                 printed and the coordinates stay exactly as given.\n";

    gc_init_rng();
    gc_reset_globals(2);
    Membrane membrane = gc_make_membrane(100.0);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 3.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 3.0, absState));

    std::vector<ForwardRxn> forwardRxns;
    gc_add_reaction(forwardRxns, molTemplateList, 0, 1, 5.0);

    Molecule mol1 = gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 });
    Molecule mol2 = gc_make_molecule(1, 1, molTemplateList[1], Coord { 2.0, 0.0, 0.0 });

    std::vector<Complex> complexList;
    complexList.push_back(gc_make_complex(0, mol1));
    complexList.push_back(gc_make_complex(1, mol2));

    bool overlap = checkAndResolveOverlap(
        mol1, mol2, molTemplateList, forwardRxns, membrane, complexList, true);

    EXPECT_FALSE(overlap) << "User supplied coordinates are never treated as an overlap";
    EXPECT_DOUBLE_EQ(mol2.comCoord.x, 2.0) << "PDB coordinates must be preserved";
    EXPECT_DOUBLE_EQ(mol2.comCoord.y, 0.0);
    EXPECT_DOUBLE_EQ(mol2.comCoord.z, 0.0);
}

// =============================================================================
// fixOverlappingMolecules()
// =============================================================================

/*! No reactions at all => the very first sweep finds nothing and the loop exits
 *  after a single iteration. */
void gc_test_fix_overlaps_no_reactions()
{
    std::cerr << "\n[TEST] gc_test_fix_overlaps_no_reactions\n"
              << "  Function    : fixOverlappingMolecules()  (full sweep)\n"
              << "  Scenario    : four molecules, no reactions defined.\n"
              << "  Pass criteria: returns 1 (single iteration, no overlaps) and no\n"
              << "                 molecule is displaced.\n";

    gc_init_rng();
    gc_reset_globals(1);
    Membrane membrane = gc_make_membrane(100.0);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 2.0, absState));

    std::vector<ForwardRxn> forwardRxns; // none

    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    for (int i = 0; i < 4; ++i) {
        moleculeList.push_back(
            gc_make_molecule(i, i, molTemplateList[0], Coord { i * 1.0, 0.0, 0.0 }));
        complexList.push_back(gc_make_complex(i, moleculeList.back()));
    }

    int iterations = fixOverlappingMolecules(
        moleculeList, molTemplateList, complexList, forwardRxns, membrane, std::vector<int> {}, false);

    EXPECT_EQ(iterations, 1) << "One clean sweep is enough when nothing can overlap";
    for (int i = 0; i < 4; ++i)
        EXPECT_DOUBLE_EQ(moleculeList[i].comCoord.x, i * 1.0) << "Molecule " << i << " must not move";
}

/*! Two genuinely overlapping molecules must be separated within the 50 iteration
 *  budget. */
void gc_test_fix_overlaps_resolves()
{
    std::cerr << "\n[TEST] gc_test_fix_overlaps_resolves\n"
              << "  Function    : fixOverlappingMolecules()  (full sweep)\n"
              << "  Scenario    : two reactive molecules whose interfaces overlap.\n"
              << "  Pass criteria: more than one iteration is required, fewer than the\n"
              << "                 50 iteration cap are used, and the final interface\n"
              << "                 separation is >= the binding radius.\n";

    gc_init_rng();
    gc_reset_globals(2);
    Membrane membrane = gc_make_membrane(200.0);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 3.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 3.0, absState));

    const double bindRadius = 5.0;
    std::vector<ForwardRxn> forwardRxns;
    gc_add_reaction(forwardRxns, molTemplateList, 0, 1, bindRadius);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(gc_make_molecule(1, 1, molTemplateList[1], Coord { 2.0, 0.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(gc_make_complex(0, moleculeList[0]));
    complexList.push_back(gc_make_complex(1, moleculeList[1]));

    int iterations = fixOverlappingMolecules(
        moleculeList, molTemplateList, complexList, forwardRxns, membrane, std::vector<int> {}, false);

    std::cerr << "  Iterations performed: " << iterations << '\n';
    EXPECT_GE(iterations, 2) << "The first sweep must detect the overlap, so >= 2 sweeps run";
    EXPECT_LE(iterations, 50) << "The routine must never exceed its iteration cap";

    if (iterations < 50) {
        // The loop only exits early when a full sweep found no overlapping pair.
        const double sq = gc_sq_dist(
            moleculeList[0].interfaceList[0].coord, moleculeList[1].interfaceList[0].coord);
        std::cerr << "  Final interface separation = " << std::sqrt(sq) << " nm (binding radius "
                  << bindRadius << ")\n";
        EXPECT_GE(sq, bindRadius * bindRadius)
            << "After a clean sweep the interfaces must be at least a binding radius apart";
    }
}

/*! In "check only changed" mode a pair where both partners came from the PDB is
 *  only warned about, so the sweep terminates immediately with no motion. */
void gc_test_fix_overlaps_only_changed_warns()
{
    std::cerr << "\n[TEST] gc_test_fix_overlaps_only_changed_warns\n"
              << "  Function    : fixOverlappingMolecules()  (checkOnlyChanged = true)\n"
              << "  Scenario    : both overlapping molecules are listed as changed.\n"
              << "  Pass criteria: returns 1 (nothing to fix, only warnings) and both\n"
              << "                 molecules keep their user supplied coordinates.\n";

    gc_init_rng();
    gc_reset_globals(2);
    Membrane membrane = gc_make_membrane(100.0);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 3.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 3.0, absState));

    std::vector<ForwardRxn> forwardRxns;
    gc_add_reaction(forwardRxns, molTemplateList, 0, 1, 5.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(gc_make_molecule(1, 1, molTemplateList[1], Coord { 2.0, 0.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(gc_make_complex(0, moleculeList[0]));
    complexList.push_back(gc_make_complex(1, moleculeList[1]));

    const std::vector<int> changed { 0, 1 };
    int iterations = fixOverlappingMolecules(
        moleculeList, molTemplateList, complexList, forwardRxns, membrane, changed, true);

    EXPECT_EQ(iterations, 1) << "Warnings only => the sweep reports no overlap and stops";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 0.0) << "PDB molecule 0 must stay put";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 2.0) << "PDB molecule 1 must stay put";
}

/*! In "check only changed" mode a pair where the partner was NOT read from the
 *  PDB is repaired by moving that partner. */
void gc_test_fix_overlaps_only_changed_moves_partner()
{
    std::cerr << "\n[TEST] gc_test_fix_overlaps_only_changed_moves_partner\n"
              << "  Function    : fixOverlappingMolecules()  (checkOnlyChanged = true)\n"
              << "  Scenario    : molecule 0 comes from the PDB, molecule 1 does not.\n"
              << "  Pass criteria: more than one iteration runs, molecule 0 stays where\n"
              << "                 the user put it and molecule 1 is relocated.\n";

    gc_init_rng();
    gc_reset_globals(2);
    Membrane membrane = gc_make_membrane(200.0);

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 3.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 3.0, absState));

    const double bindRadius = 5.0;
    std::vector<ForwardRxn> forwardRxns;
    gc_add_reaction(forwardRxns, molTemplateList, 0, 1, bindRadius);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(gc_make_molecule(0, 0, molTemplateList[0], Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(gc_make_molecule(1, 1, molTemplateList[1], Coord { 2.0, 0.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(gc_make_complex(0, moleculeList[0]));
    complexList.push_back(gc_make_complex(1, moleculeList[1]));

    const Coord partnerBefore = moleculeList[1].comCoord;

    const std::vector<int> changed { 0 };
    int iterations = fixOverlappingMolecules(
        moleculeList, molTemplateList, complexList, forwardRxns, membrane, changed, true);

    std::cerr << "  Iterations performed: " << iterations << '\n';
    EXPECT_GE(iterations, 2) << "An overlap is found on the first sweep, so a retry happens";
    EXPECT_LE(iterations, 50) << "The iteration cap must be respected";

    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 0.0)
        << "The molecule taken from the PDB must never be relocated";

    const bool partnerMoved = (moleculeList[1].comCoord.x != partnerBefore.x)
        || (moleculeList[1].comCoord.y != partnerBefore.y)
        || (moleculeList[1].comCoord.z != partnerBefore.z);
    EXPECT_TRUE(partnerMoved) << "The non-PDB partner should have been re-randomised";
}

// =============================================================================
// generate_coordinates()
// =============================================================================

/*! The driver must create exactly `copies` molecules (and one complex each) for
 *  every non implicit-lipid template and register them in monomerList. */
void gc_test_generate_coordinates_builds_system()
{
    std::cerr << "\n[TEST] gc_test_generate_coordinates_builds_system\n"
              << "  Function    : generate_coordinates()\n"
              << "  Scenario    : two molecule templates with 3 and 2 copies, no\n"
              << "                reactions and no coordinate file.\n"
              << "  Pass criteria: 5 molecules and 5 complexes are created with the\n"
              << "                 right types, monomerLists are filled, molecule and\n"
              << "                 complex indices agree and everything lies in the box.\n";

    gc_init_rng();
    gc_reset_globals(2);
    // write_xyz() dumps into DATA/, make sure the directory exists.
    (void)std::system("mkdir -p DATA");

    Membrane membrane = gc_make_membrane(100.0);

    Parameters params;
    params.timeStep = 1.0;
    params.numMolTypes = 2;
    params.numTotalUnits = 5;

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 2.0, absState));
    molTemplateList.push_back(gc_make_template("B", 1, 1, 2.0, absState));
    molTemplateList[0].copies = 3;
    molTemplateList[1].copies = 2;

    std::vector<ForwardRxn> forwardRxns; // no reactions -> no overlap handling
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::string coordinateFileName; // empty -> PDB branch skipped

    std::cerr << "  Calling generate_coordinates()...\n";
    generate_coordinates(params, moleculeList, complexList, molTemplateList, forwardRxns, membrane,
        coordinateFileName);

    EXPECT_EQ(moleculeList.size(), 5u) << "3 copies of A + 2 copies of B = 5 molecules";
    EXPECT_EQ(complexList.size(), 5u) << "Every free molecule starts in its own complex";

    EXPECT_EQ(molTemplateList[0].monomerList.size(), 3u) << "Template A monomerList size";
    EXPECT_EQ(molTemplateList[1].monomerList.size(), 2u) << "Template B monomerList size";

    // Types are laid down template by template: indices 0-2 are A, 3-4 are B.
    for (size_t i = 0; i < moleculeList.size(); ++i) {
        const int expectedType = (i < 3) ? 0 : 1;
        EXPECT_EQ(moleculeList[i].molTypeIndex, expectedType)
            << "Molecule " << i << " should belong to template " << expectedType;

        // Each molecule got one interface, matching the template.
        EXPECT_EQ(moleculeList[i].interfaceList.size(), 1u)
            << "Molecule " << i << " should have one interface";

        // Coordinates must be inside the 100 nm box.
        EXPECT_LE(std::abs(moleculeList[i].comCoord.x), 50.0) << "Molecule " << i << " x in box";
        EXPECT_LE(std::abs(moleculeList[i].comCoord.y), 50.0) << "Molecule " << i << " y in box";
        EXPECT_LE(std::abs(moleculeList[i].comCoord.z), 50.0) << "Molecule " << i << " z in box";

        // The complex bookkeeping created alongside the molecule must line up.
        ASSERT_LT(moleculeList[i].myComIndex, static_cast<int>(complexList.size()))
            << "Molecule " << i << " points at a valid complex";
        const Complex& com = complexList[moleculeList[i].myComIndex];
        EXPECT_EQ(com.memberList.size(), 1u) << "Fresh complexes hold exactly one molecule";
        if (!com.memberList.empty())
            EXPECT_EQ(com.memberList[0], moleculeList[i].index)
                << "Complex member index should be the molecule index";
        EXPECT_EQ(moleculeList[i].complexId, com.id) << "complexId must be copied from the complex";
    }

    std::cerr << "  Created " << moleculeList.size() << " molecules and " << complexList.size()
              << " complexes.\n";
}

/*! An empty template list with numTotalUnits == 0 takes the early-return branch. */
void gc_test_generate_coordinates_empty_system()
{
    std::cerr << "\n[TEST] gc_test_generate_coordinates_empty_system\n"
              << "  Function    : generate_coordinates()\n"
              << "  Scenario    : no templates at all and params.numTotalUnits == 0.\n"
              << "  Pass criteria: the routine returns early, leaving both lists empty\n"
              << "                 and writing no coordinate file.\n";

    gc_init_rng();
    gc_reset_globals(0);

    Membrane membrane = gc_make_membrane(100.0);

    Parameters params;
    params.timeStep = 1.0;
    params.numMolTypes = 0;
    params.numTotalUnits = 0;

    std::vector<MolTemplate> molTemplateList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::string coordinateFileName;

    generate_coordinates(params, moleculeList, complexList, molTemplateList, forwardRxns, membrane,
        coordinateFileName);

    EXPECT_TRUE(moleculeList.empty()) << "No templates => no molecules";
    EXPECT_TRUE(complexList.empty()) << "No templates => no complexes";
}

/*! With a coordinate file, generate_coordinates() must place the first molecules
 *  of the named type exactly on the supplied COM positions. */
void gc_test_generate_coordinates_applies_pdb()
{
    std::cerr << "\n[TEST] gc_test_generate_coordinates_applies_pdb\n"
              << "  Function    : generate_coordinates() (PDB branch)\n"
              << "  Scenario    : 3 copies of A, a PDB file with 2 COM records for A.\n"
              << "  Pass criteria: molecules 0 and 1 land exactly on the PDB values,\n"
              << "                 molecule 2 keeps its random position, and every\n"
              << "                 molecule stays rigid (1 nm interface offset).\n";

    gc_init_rng();
    gc_reset_globals(1);
    (void)std::system("mkdir -p DATA");

    Membrane membrane = gc_make_membrane(100.0);

    Parameters params;
    params.timeStep = 1.0;
    params.numMolTypes = 1;
    params.numTotalUnits = 3;

    int absState = 0;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(gc_make_template("A", 0, 1, 2.0, absState));
    molTemplateList[0].copies = 3;

    std::vector<ForwardRxn> forwardRxns; // no reactions -> the PDB positions survive
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;

    const std::string fileName = "test_gc_generate_coords.pdb";
    {
        std::ofstream out(fileName);
        gc_write_pdb_line(out, 1, "COM", "A", 1, 10.0, -20.0, 5.0);
        gc_write_pdb_line(out, 2, "COM", "A", 2, -30.0, 15.0, -5.0);
        out.close();
    }
    std::string coordinateFileName = fileName;

    generate_coordinates(params, moleculeList, complexList, molTemplateList, forwardRxns, membrane,
        coordinateFileName);

    ASSERT_EQ(moleculeList.size(), 3u) << "Three copies of A must exist";

    // The first two molecules must have been snapped onto the PDB coordinates.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 10.0) << "Molecule 0 x from PDB";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, -20.0) << "Molecule 0 y from PDB";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 5.0) << "Molecule 0 z from PDB";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, -30.0) << "Molecule 1 x from PDB";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.y, 15.0) << "Molecule 1 y from PDB";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.z, -5.0) << "Molecule 1 z from PDB";

    // The third molecule keeps whatever random position it was given, which must
    // still be inside the box.
    EXPECT_LE(std::abs(moleculeList[2].comCoord.x), 50.0) << "Molecule 2 x remains in the box";
    EXPECT_LE(std::abs(moleculeList[2].comCoord.y), 50.0) << "Molecule 2 y remains in the box";
    EXPECT_LE(std::abs(moleculeList[2].comCoord.z), 50.0) << "Molecule 2 z remains in the box";

    // The PDB update is a rigid translation, so the interface offset is intact.
    for (size_t i = 0; i < moleculeList.size(); ++i) {
        const double offset
            = std::sqrt(gc_sq_dist(moleculeList[i].interfaceList[0].coord, moleculeList[i].comCoord));
        EXPECT_NEAR(offset, 1.0, 1e-9)
            << "Molecule " << i << " must keep its 1 nm interface offset";
    }

    std::cerr << "  Molecule 2 (untouched by the PDB) COM = (" << moleculeList[2].comCoord.x << ", "
              << moleculeList[2].comCoord.y << ", " << moleculeList[2].comCoord.z << ")\n";

    std::remove(fileName.c_str());
}

// =============================================================================
// GoogleTest wrappers - one per named test function so that a failure in one
// case never prevents the remaining cases from running.
// =============================================================================
TEST(GenerateCoordinatesTest, UpdateMoleculeCoordinatesFromPdb) { gc_test_update_coordinates_from_pdb(); }
TEST(GenerateCoordinatesTest, UpdateMoleculeCoordinatesMissingFile) { gc_test_update_coordinates_missing_file(); }
TEST(GenerateCoordinatesTest, UpdateMoleculeCoordinatesUnknownAndSurplus) { gc_test_update_coordinates_unknown_and_surplus(); }
TEST(GenerateCoordinatesTest, CheckAndResolveOverlapNoReaction) { gc_test_overlap_no_reaction(); }
TEST(GenerateCoordinatesTest, CheckAndResolveOverlapFarApart) { gc_test_overlap_far_apart(); }
TEST(GenerateCoordinatesTest, CheckAndResolveOverlapRelocates) { gc_test_overlap_relocates_molecule(); }
TEST(GenerateCoordinatesTest, CheckAndResolveOverlapWarnsForPdbMolecule) { gc_test_overlap_warns_for_pdb_molecule(); }
TEST(GenerateCoordinatesTest, FixOverlappingNoReactions) { gc_test_fix_overlaps_no_reactions(); }
TEST(GenerateCoordinatesTest, FixOverlappingResolves) { gc_test_fix_overlaps_resolves(); }
TEST(GenerateCoordinatesTest, FixOverlappingOnlyChangedWarns) { gc_test_fix_overlaps_only_changed_warns(); }
TEST(GenerateCoordinatesTest, FixOverlappingOnlyChangedMovesPartner) { gc_test_fix_overlaps_only_changed_moves_partner(); }
TEST(GenerateCoordinatesTest, GenerateCoordinatesBuildsSystem) { gc_test_generate_coordinates_builds_system(); }
TEST(GenerateCoordinatesTest, GenerateCoordinatesEmptySystem) { gc_test_generate_coordinates_empty_system(); }
TEST(GenerateCoordinatesTest, GenerateCoordinatesAppliesPdb) { gc_test_generate_coordinates_applies_pdb(); }
