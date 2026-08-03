/*! \file test_create_molecule_and_complex_from_rxn.cpp
 *
 * ### Unit tests for src/reactions/create_molecule_and_complex_from_rxn.cpp
 *
 * The file under test provides three free functions:
 *
 *   1. bool moleculeOverlaps(...)
 *        - decides whether a freshly placed Molecule is invalid, either because
 *          it left the simulation volume / membrane, or because it landed inside
 *          the binding radius of a molecule it can react with.
 *
 *   2. void create_molecule_and_complex_from_rxn(...)
 *        - allocates (or recycles) a Molecule slot and a Complex slot, places the
 *          new molecule (randomly in the box, or in the vicinity of a parent),
 *          resampling while moleculeOverlaps() reports a bad placement, and then
 *          wires up the molecule <-> complex bookkeeping.
 *
 *   3. void MPI_create_molecule_and_complex_on_rank(...)
 *        - the MPI variant which only allocates the slots and links them, without
 *          generating coordinates.
 *
 * Every test below prints which source file / function is exercised, what the
 * scenario is, and what the pass criteria are, so the console log is readable on
 * its own.  All assertions are non-fatal (EXPECT_*) so a single failure never
 * prevents the remaining checks from running.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "classes/class_SimulVolume.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

// -----------------------------------------------------------------------------
// The MPI helper is defined in the file under test but is not declared in any
// project header, so declare it here to be able to link against it.
// -----------------------------------------------------------------------------
void MPI_create_molecule_and_complex_on_rank(Molecule& mol, int& newMolIndex,
    int& newComIndex, MolTemplate& createdMolTemp, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);

// =============================================================================
// Small fixture-style helpers (kept in an anonymous namespace so they cannot
// collide with helpers from other translation units in the test binary).
// =============================================================================
namespace {

//! Side length of the cubic water box used by every test below.
constexpr double kBoxSide = 100.0;

/*! \brief Make sure the GSL random number generator used by the placement
 *  routines has been seeded. `r` itself is defined in gtest_main.cpp.
 */
void cmacr_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] GSL rng pointer was null -> calling srand_gsl(1)\n";
        srand_gsl(1);
    }
}

/*! \brief Reset the static bookkeeping shared by Molecule/Complex/MolTemplate.
 *
 * The functions under test read and write these statics (empty-slot lists,
 * complex counters, number of molecule types), so each test starts from a known
 * state.
 *
 * \param[in] numMolTypes number of MolTemplates the test will use.
 */
void cmacr_reset_statics(int numMolTypes)
{
    MolTemplate::numMolTypes = static_cast<unsigned>(numMolTypes);
    MolTemplate::numEachMolType.assign(numMolTypes, 0);
    Molecule::numberOfMolecules = 0;
    Molecule::emptyMolList.clear();
    Complex::numberOfComplexes = 0;
    Complex::emptyComList.clear();
    Complex::currNumberMolTypes = numMolTypes;
    Complex::currNumberComTypes = numMolTypes;
}

/*! \brief Build a cubic, non-spherical Membrane of side kBoxSide centered on 0. */
Membrane cmacr_make_membrane()
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.waterBox.x = kBoxSide;
    membraneObject.waterBox.y = kBoxSide;
    membraneObject.waterBox.z = kBoxSide;
    membraneObject.waterBox.volume = kBoxSide * kBoxSide * kBoxSide;
    return membraneObject;
}

/*! \brief Build a SimulVolume consisting of exactly one sub-cell spanning the
 *  whole water box.  Using a single cell makes the bin arithmetic inside
 *  moleculeOverlaps() completely deterministic (every in-box molecule -> bin 0).
 */
SimulVolume cmacr_make_simul_volume()
{
    SimulVolume simulVolume;
    simulVolume.numSubCells.x = 1;
    simulVolume.numSubCells.y = 1;
    simulVolume.numSubCells.z = 1;
    simulVolume.numSubCells.tot = 1;
    simulVolume.subCellSize = Coord { kBoxSide, kBoxSide, kBoxSide };
    simulVolume.subCellList.clear();
    simulVolume.subCellList.resize(1);
    simulVolume.subCellList[0].absIndex = 0;
    simulVolume.subCellList[0].xIndex = 0;
    simulVolume.subCellList[0].yIndex = 0;
    simulVolume.subCellList[0].zIndex = 0;
    return simulVolume;
}

/*! \brief Build a minimal MolTemplate with a single, single-state interface.
 *
 * \param[in] molTypeIndex index of the template in molTemplateList.
 * \param[in] name         template name (only used for diagnostics).
 * \param[in] isLipid      flag the template as a (membrane bound) lipid.
 */
MolTemplate cmacr_make_template(int molTypeIndex, const std::string& name, bool isLipid)
{
    MolTemplate molTemplate;
    molTemplate.molTypeIndex = molTypeIndex;
    molTemplate.molName = name;
    molTemplate.mass = 1.0;
    molTemplate.radius = 2.0;
    molTemplate.comCoord = Coord { 0.0, 0.0, 0.0 };
    molTemplate.D = Coord { 10.0, 10.0, isLipid ? 0.0 : 10.0 };
    molTemplate.Dr = Coord { 0.1, 0.1, 0.1 };
    molTemplate.isLipid = isLipid;
    molTemplate.isImplicitLipid = false;
    molTemplate.copies = 1;

    // One interface, offset 1 nm along +x from the center of mass, with one state.
    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord { 1.0, 0.0, 0.0 };
    iface.stateList.push_back(Interface::State('\0', 0));
    molTemplate.interfaceList.push_back(iface);

    return molTemplate;
}

/*! \brief Build a Molecule with a single interface at an explicit coordinate. */
Molecule cmacr_make_molecule(int index, int comIndex, int molTypeIndex,
    const Coord& com, const Coord& ifaceCoord, int ifaceAbsIndex)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = ifaceCoord;
    iface.index = ifaceAbsIndex;
    iface.relIndex = 0;
    iface.molTypeIndex = molTypeIndex;
    iface.stateIden = '\0';
    iface.stateIndex = 0;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a Complex owning one member molecule. */
Complex cmacr_make_complex(int index, const Coord& com, double radius, int memberIndex)
{
    Complex targCom;
    targCom.index = index;
    targCom.comCoord = com;
    targCom.radius = radius;
    targCom.mass = 1.0;
    targCom.D = Coord { 10.0, 10.0, 10.0 };
    targCom.Dr = Coord { 0.1, 0.1, 0.1 };
    targCom.memberList.push_back(memberIndex);
    targCom.numEachMol.assign(MolTemplate::numMolTypes, 0);
    targCom.lastNumberUpdateItrEachMol.assign(MolTemplate::numMolTypes, 0);
    return targCom;
}

/*! \brief Build a symmetric bimolecular ForwardRxn A(a) + A(a) -> A(a!1).A(a!1).
 *
 * Both entries of reactantListNew must exist because moleculeOverlaps() indexes
 * reactantListNew[0] and reactantListNew[1] unconditionally.
 */
ForwardRxn cmacr_make_forward_rxn(double bindRadius)
{
    ForwardRxn oneRxn;
    oneRxn.rxnType = ReactionType::bimolecular;
    oneRxn.bindRadius = bindRadius;
    oneRxn.isSymmetric = true;

    RxnIface reactant;
    reactant.ifaceName = "a";
    reactant.molTypeIndex = 0;
    reactant.absIfaceIndex = 0;
    reactant.relIfaceIndex = 0;
    reactant.requiresState = '\0';
    reactant.requiresInteraction = false;

    oneRxn.reactantListNew.push_back(reactant);
    oneRxn.reactantListNew.push_back(reactant);
    oneRxn.productListNew.push_back(reactant);
    oneRxn.rateList.emplace_back();
    return oneRxn;
}

/*! \brief Build a zeroth-order creation reaction that produces molType 0. */
CreateDestructRxn cmacr_make_creation_rxn(double creationRadius)
{
    CreateDestructRxn currRxn;
    currRxn.rxnType = ReactionType::zerothOrderCreation;
    currRxn.creationRadius = creationRadius;
    currRxn.absRxnIndex = 0;
    currRxn.relRxnIndex = 0;
    currRxn.rateList.emplace_back();
    currRxn.rateList.back().rate = 1.0;

    // The product molecule: molType 0 with its single interface in state 0.
    RxnIface prodIface;
    prodIface.ifaceName = "a";
    prodIface.molTypeIndex = 0;
    prodIface.absIfaceIndex = 0;
    prodIface.relIfaceIndex = 0;
    prodIface.requiresState = '\0';
    prodIface.requiresInteraction = false;

    CreateDestructRxn::CreateDestructMol prodMol;
    prodMol.molTypeIndex = 0;
    prodMol.molName = "A";
    prodMol.interfaceList.push_back(prodIface);
    currRxn.productMolList.push_back(prodMol);

    return currRxn;
}

/*! \brief Simple Parameters object good enough for the placement routines. */
Parameters cmacr_make_params()
{
    Parameters params;
    params.timeStep = 1.0;
    params.numMolTypes = 1;
    params.numTotalSpecies = 1;
    params.rMaxLimit = 10.0;
    params.rMaxRadius = 10.0;
    return params;
}

} // namespace

// =============================================================================
// moleculeOverlaps() tests
// =============================================================================

/*! Molecule sits at the origin of an empty box: no overlap, and the molecule is
 *  registered into the sub-volume it landed in.
 */
void test_cmacr_overlap_inside_empty_bin()
{
    std::cerr << "\n[TEST] test_cmacr_overlap_inside_empty_bin\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      moleculeOverlaps()\n"
              << "  Scenario:      lone molecule at the center of an empty box.\n"
              << "  Pass criteria: returns false (valid placement), mySubVolIndex\n"
              << "                 is set to bin 0 and the molecule index is added\n"
              << "                 to that sub-volume's memberMolList.\n";

    cmacr_reset_statics(1);
    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };
    std::vector<ForwardRxn> forwardRxns {}; // no reactions -> only geometry matters
    std::vector<Complex> complexList {};
    std::vector<Molecule> moleculeList {};

    Molecule createdMol = cmacr_make_molecule(
        0, -1, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 0);

    bool overlaps = moleculeOverlaps(params, simulVolume, createdMol, moleculeList,
        complexList, forwardRxns, molTemplateList, membraneObject);

    std::cerr << "  moleculeOverlaps() returned " << std::boolalpha << overlaps << '\n';
    EXPECT_FALSE(overlaps) << "A molecule alone at the box center must be a valid placement";
    EXPECT_EQ(createdMol.mySubVolIndex, 0) << "The molecule should be assigned to sub-volume 0";
    ASSERT_NO_FATAL_FAILURE(); // keep going even if the above failed
    EXPECT_EQ(simulVolume.subCellList[0].memberMolList.size(), 1u)
        << "The molecule should have been appended to the sub-volume member list";
    if (!simulVolume.subCellList[0].memberMolList.empty()) {
        EXPECT_EQ(simulVolume.subCellList[0].memberMolList[0], createdMol.index)
            << "The registered index should be the created molecule's index";
    }
}

/*! Molecules pushed outside the box in x, y or z must be rejected. */
void test_cmacr_overlap_outside_box()
{
    std::cerr << "\n[TEST] test_cmacr_overlap_outside_box\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      moleculeOverlaps()\n"
              << "  Scenario:      molecule placed beyond the +x, -y and +z walls\n"
              << "                 of a 100 nm cubic box.\n"
              << "  Pass criteria: returns true for each dimension and never\n"
              << "                 registers itself in a sub-volume.\n";

    cmacr_reset_statics(1);
    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<Complex> complexList {};
    std::vector<Molecule> moleculeList {};

    // --- outside in +x -------------------------------------------------------
    Molecule outX = cmacr_make_molecule(
        0, -1, 0, Coord { 60.0, 0.0, 0.0 }, Coord { 61.0, 0.0, 0.0 }, 0);
    bool resX = moleculeOverlaps(params, simulVolume, outX, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);
    std::cerr << "  x = 60 (wall at 50) -> " << std::boolalpha << resX << '\n';
    EXPECT_TRUE(resX) << "x = 60 nm lies outside the +x wall and must be rejected";

    // --- outside in -y -------------------------------------------------------
    Molecule outY = cmacr_make_molecule(
        0, -1, 0, Coord { 0.0, -60.0, 0.0 }, Coord { 1.0, -60.0, 0.0 }, 0);
    bool resY = moleculeOverlaps(params, simulVolume, outY, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);
    std::cerr << "  y = -60 (wall at -50) -> " << resY << '\n';
    EXPECT_TRUE(resY) << "y = -60 nm lies outside the -y wall and must be rejected";

    // --- outside in +z -------------------------------------------------------
    Molecule outZ = cmacr_make_molecule(
        0, -1, 0, Coord { 0.0, 0.0, 60.0 }, Coord { 1.0, 0.0, 60.0 }, 0);
    bool resZ = moleculeOverlaps(params, simulVolume, outZ, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);
    std::cerr << "  z = 60 (wall at 50) -> " << resZ << '\n';
    EXPECT_TRUE(resZ) << "z = 60 nm lies outside the +z wall and must be rejected";

    // Nothing should have been added to the (single) sub-volume.
    EXPECT_TRUE(simulVolume.subCellList[0].memberMolList.empty())
        << "Rejected molecules must not be registered in a sub-volume";
}

/*! Lipids must remain on the membrane plane (|z| == box.z/2). */
void test_cmacr_overlap_lipid_membrane_check()
{
    std::cerr << "\n[TEST] test_cmacr_overlap_lipid_membrane_check\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      moleculeOverlaps()\n"
              << "  Scenario:      a lipid-flagged molecule exactly on the lower\n"
              << "                 membrane (z = -50) and one below it (z = -60).\n"
              << "  Pass criteria: on-membrane lipid accepted (false), lipid that\n"
              << "                 has left the membrane rejected (true).\n";

    cmacr_reset_statics(2);
    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    // template 0 -> normal protein, template 1 -> lipid
    std::vector<MolTemplate> molTemplateList {
        cmacr_make_template(0, "A", false),
        cmacr_make_template(1, "LIPID", true),
    };
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<Complex> complexList {};
    std::vector<Molecule> moleculeList {};

    // Lipid exactly on the membrane -> accepted.
    Molecule onMem = cmacr_make_molecule(
        0, -1, 1, Coord { 0.0, 0.0, -50.0 }, Coord { 1.0, 0.0, -50.0 }, 0);
    onMem.isLipid = true;
    bool resOn = moleculeOverlaps(params, simulVolume, onMem, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);
    std::cerr << "  lipid at z = -50 -> " << std::boolalpha << resOn << '\n';
    EXPECT_FALSE(resOn) << "A lipid sitting on the membrane is a legal placement";

    // Lipid that fell off the membrane -> rejected.
    Molecule offMem = cmacr_make_molecule(
        1, -1, 1, Coord { 0.0, 0.0, -60.0 }, Coord { 1.0, 0.0, -60.0 }, 0);
    offMem.isLipid = true;
    bool resOff = moleculeOverlaps(params, simulVolume, offMem, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);
    std::cerr << "  lipid at z = -60 -> " << resOff << '\n';
    EXPECT_TRUE(resOff) << "A lipid below the membrane must be rejected";
}

/*! A reactive partner inside the binding radius counts as an overlap. */
void test_cmacr_overlap_reactive_partner_inside_bindradius()
{
    std::cerr << "\n[TEST] test_cmacr_overlap_reactive_partner_inside_bindradius\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      moleculeOverlaps()\n"
              << "  Scenario:      new molecule lands 0.5 nm from a molecule it can\n"
              << "                 bind to, with a binding radius of 5 nm.\n"
              << "  Pass criteria: returns true (placement rejected).\n";

    cmacr_reset_statics(1);
    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };
    std::vector<ForwardRxn> forwardRxns { cmacr_make_forward_rxn(5.0) };

    // Existing molecule / complex, close to the origin and inside the same bin.
    std::vector<Molecule> moleculeList { cmacr_make_molecule(
        0, 0, 0, Coord { 1.0, 0.0, 0.0 }, Coord { 1.5, 0.0, 0.0 }, 0) };
    std::vector<Complex> complexList { cmacr_make_complex(
        0, Coord { 1.0, 0.0, 0.0 }, 10.0, 0) };
    simulVolume.subCellList[0].memberMolList.push_back(0);

    // New molecule: interface at (1,0,0) -> 0.5 nm from the partner interface.
    Molecule createdMol = cmacr_make_molecule(
        1, -1, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 0);

    bool overlaps = moleculeOverlaps(params, simulVolume, createdMol, moleculeList,
        complexList, forwardRxns, molTemplateList, membraneObject);

    std::cerr << "  interface separation = 0.5 nm, bindRadius = 5 nm -> "
              << std::boolalpha << overlaps << '\n';
    EXPECT_TRUE(overlaps) << "Interfaces closer than the binding radius must be an overlap";
}

/*! The same geometry with a tiny binding radius is a legal placement. */
void test_cmacr_overlap_reactive_partner_outside_bindradius()
{
    std::cerr << "\n[TEST] test_cmacr_overlap_reactive_partner_outside_bindradius\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      moleculeOverlaps()\n"
              << "  Scenario:      same neighbour as before, but the binding radius\n"
              << "                 is only 0.1 nm (separation is 0.5 nm).\n"
              << "  Pass criteria: returns false and the molecule is registered in\n"
              << "                 the sub-volume member list.\n";

    cmacr_reset_statics(1);
    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };
    std::vector<ForwardRxn> forwardRxns { cmacr_make_forward_rxn(0.1) };

    std::vector<Molecule> moleculeList { cmacr_make_molecule(
        0, 0, 0, Coord { 1.0, 0.0, 0.0 }, Coord { 1.5, 0.0, 0.0 }, 0) };
    std::vector<Complex> complexList { cmacr_make_complex(
        0, Coord { 1.0, 0.0, 0.0 }, 10.0, 0) };
    simulVolume.subCellList[0].memberMolList.push_back(0);

    Molecule createdMol = cmacr_make_molecule(
        1, -1, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 0);

    bool overlaps = moleculeOverlaps(params, simulVolume, createdMol, moleculeList,
        complexList, forwardRxns, molTemplateList, membraneObject);

    std::cerr << "  interface separation = 0.5 nm, bindRadius = 0.1 nm -> "
              << std::boolalpha << overlaps << '\n';
    EXPECT_FALSE(overlaps) << "Interfaces farther apart than the binding radius do not overlap";
    EXPECT_EQ(createdMol.mySubVolIndex, 0) << "Accepted molecule should be binned into cell 0";
    EXPECT_EQ(simulVolume.subCellList[0].memberMolList.size(), 2u)
        << "The accepted molecule should be appended after the pre-existing member";
}

/*! Neighbours outside the combined bounding spheres short-circuit to "no overlap"
 *  without registering the molecule - this documents the current behaviour.
 */
void test_cmacr_overlap_bounding_sphere_shortcut()
{
    std::cerr << "\n[TEST] test_cmacr_overlap_bounding_sphere_shortcut\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      moleculeOverlaps()\n"
              << "  Scenario:      the only neighbour in the bin is 80 nm away while\n"
              << "                 the bounding spheres sum to 12 nm.\n"
              << "  Pass criteria: returns false immediately (bounding-sphere\n"
              << "                 shortcut) and leaves mySubVolIndex untouched.\n";

    cmacr_reset_statics(1);
    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };
    // A reaction with a huge binding radius: it must never even be evaluated.
    std::vector<ForwardRxn> forwardRxns { cmacr_make_forward_rxn(1000.0) };

    std::vector<Molecule> moleculeList { cmacr_make_molecule(
        0, 0, 0, Coord { 40.0, 0.0, 0.0 }, Coord { 41.0, 0.0, 0.0 }, 0) };
    std::vector<Complex> complexList { cmacr_make_complex(
        0, Coord { 40.0, 0.0, 0.0 }, 10.0, 0) };
    simulVolume.subCellList[0].memberMolList.push_back(0);

    Molecule createdMol = cmacr_make_molecule(
        1, -1, 0, Coord { -40.0, 0.0, 0.0 }, Coord { -39.0, 0.0, 0.0 }, 0);

    bool overlaps = moleculeOverlaps(params, simulVolume, createdMol, moleculeList,
        complexList, forwardRxns, molTemplateList, membraneObject);

    std::cerr << "  COM separation = 80 nm, radius sum = 12 nm -> "
              << std::boolalpha << overlaps << '\n';
    EXPECT_FALSE(overlaps) << "A distant neighbour cannot cause an overlap";
    EXPECT_EQ(createdMol.mySubVolIndex, -1)
        << "The bounding-sphere shortcut returns before the molecule is binned";
}

// =============================================================================
// create_molecule_and_complex_from_rxn() tests
// =============================================================================

/*! Zeroth-order creation: molecule placed randomly in the box, new complex made. */
void test_cmacr_create_zeroth_order()
{
    std::cerr << "\n[TEST] test_cmacr_create_zeroth_order\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      create_molecule_and_complex_from_rxn()\n"
              << "  Scenario:      createInVicinity == false, empty system, no\n"
              << "                 reactions, so the first random placement is kept.\n"
              << "  Pass criteria: new Molecule and Complex slots appended, indices\n"
              << "                 cross-linked, trajStatus propagated, isDissociated\n"
              << "                 set, coordinates inside the water box, and the\n"
              << "                 molecule pushed onto the template monomerList.\n";

    cmacr_ensure_rng();
    cmacr_reset_statics(1);

    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };
    molTemplateList[0].canDestroy = true; // so the monomerList bookkeeping runs
    std::vector<ForwardRxn> forwardRxns {};
    CreateDestructRxn currRxn = cmacr_make_creation_rxn(3.0);

    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};

    int newMolIndex = -1;
    int newComIndex = -1;
    const int complexesBefore = Complex::numberOfComplexes;

    std::cerr << "  Calling create_molecule_and_complex_from_rxn (zeroth order)...\n";
    create_molecule_and_complex_from_rxn(/*parentMolIndex=*/0, newMolIndex, newComIndex,
        /*createInVicinity=*/false, molTemplateList[0], params, currRxn, simulVolume,
        moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);

    std::cerr << "  newMolIndex = " << newMolIndex << ", newComIndex = " << newComIndex << '\n';

    EXPECT_EQ(newMolIndex, 0) << "First molecule in an empty system should take slot 0";
    EXPECT_EQ(newComIndex, 0) << "First complex in an empty system should take slot 0";
    EXPECT_EQ(moleculeList.size(), 1u) << "One molecule slot should have been appended";
    EXPECT_EQ(complexList.size(), 1u) << "One complex slot should have been appended";

    if (newMolIndex >= 0 && newMolIndex < static_cast<int>(moleculeList.size())
        && newComIndex >= 0 && newComIndex < static_cast<int>(complexList.size())) {
        const Molecule& newMol = moleculeList[newMolIndex];
        const Complex& newCom = complexList[newComIndex];

        std::cerr << "  created COM = (" << newMol.comCoord.x << ", " << newMol.comCoord.y
                  << ", " << newMol.comCoord.z << ")\n";

        EXPECT_EQ(newMol.myComIndex, newComIndex) << "Molecule must point at its new complex";
        EXPECT_EQ(newMol.molTypeIndex, 0) << "Molecule type should come from the template";
        EXPECT_EQ(static_cast<int>(newMol.trajStatus), static_cast<int>(TrajStatus::propagated))
            << "New molecules are flagged as already propagated this step";
        EXPECT_TRUE(newMol.isDissociated) << "New molecules are flagged isDissociated";
        EXPECT_FALSE(newMol.isGhosted) << "New molecules are owned by this rank";
        EXPECT_FALSE(newMol.interfaceList.empty()) << "Molecule should have its interface built";

        EXPECT_EQ(static_cast<int>(newCom.trajStatus), static_cast<int>(TrajStatus::propagated))
            << "New complexes are flagged as already propagated this step";
        EXPECT_EQ(newCom.memberList.size(), 1u) << "New complex holds exactly one molecule";
        if (!newCom.memberList.empty()) {
            EXPECT_EQ(newCom.memberList[0], newMolIndex)
                << "The complex's member should be the created molecule";
        }
        EXPECT_EQ(newMol.complexId, newCom.id) << "Molecule complexId must match the complex id";

        // Placement must respect the water box.
        EXPECT_LE(std::abs(newMol.comCoord.x), kBoxSide / 2.0 + 1e-6)
            << "x coordinate must be inside the water box";
        EXPECT_LE(std::abs(newMol.comCoord.y), kBoxSide / 2.0 + 1e-6)
            << "y coordinate must be inside the water box";
        EXPECT_LE(std::abs(newMol.comCoord.z), kBoxSide / 2.0 + 1e-6)
            << "z coordinate must be inside the water box";
    }

    EXPECT_EQ(Complex::numberOfComplexes, complexesBefore + 1)
        << "Complex::numberOfComplexes should be incremented once";
    EXPECT_EQ(molTemplateList[0].monomerList.size(), 1u)
        << "canDestroy templates track new monomers in monomerList";
    if (!molTemplateList[0].monomerList.empty()) {
        EXPECT_EQ(molTemplateList[0].monomerList[0], newMolIndex)
            << "monomerList should record the new molecule index";
    }
}

/*! Unimolecular creation: the product is placed within creationRadius of a parent. */
void test_cmacr_create_in_vicinity()
{
    std::cerr << "\n[TEST] test_cmacr_create_in_vicinity\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      create_molecule_and_complex_from_rxn()\n"
              << "  Scenario:      createInVicinity == true with a parent molecule at\n"
              << "                 the origin and creationRadius = 3 nm.\n"
              << "  Pass criteria: the product is appended after the parent and its\n"
              << "                 COM lies within creationRadius of the parent.\n";

    cmacr_ensure_rng();
    cmacr_reset_statics(1);

    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };
    std::vector<ForwardRxn> forwardRxns {}; // keep the resampling loop deterministic
    const double creationRadius = 3.0;
    CreateDestructRxn currRxn = cmacr_make_creation_rxn(creationRadius);

    // Parent molecule / complex at the origin.
    std::vector<Molecule> moleculeList { cmacr_make_molecule(
        0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 0) };
    std::vector<Complex> complexList { cmacr_make_complex(
        0, Coord { 0.0, 0.0, 0.0 }, 2.0, 0) };

    int newMolIndex = -1;
    int newComIndex = -1;

    std::cerr << "  Calling create_molecule_and_complex_from_rxn (in vicinity)...\n";
    create_molecule_and_complex_from_rxn(/*parentMolIndex=*/0, newMolIndex, newComIndex,
        /*createInVicinity=*/true, molTemplateList[0], params, currRxn, simulVolume,
        moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);

    EXPECT_EQ(newMolIndex, 1) << "Product should occupy the next free molecule slot";
    EXPECT_EQ(newComIndex, 1) << "Product should occupy the next free complex slot";

    if (newMolIndex >= 0 && newMolIndex < static_cast<int>(moleculeList.size())) {
        const Molecule& newMol = moleculeList[newMolIndex];
        const double dist = std::sqrt(newMol.comCoord.x * newMol.comCoord.x
            + newMol.comCoord.y * newMol.comCoord.y + newMol.comCoord.z * newMol.comCoord.z);
        std::cerr << "  distance from parent = " << dist
                  << " (creationRadius = " << creationRadius << ")\n";
        EXPECT_LE(dist, creationRadius + 1e-6)
            << "Product must be created within creationRadius of the parent";
        EXPECT_EQ(newMol.myComIndex, newComIndex) << "Molecule must point at its new complex";
        EXPECT_TRUE(newMol.isDissociated) << "New molecules are flagged isDissociated";
    }

    // The parent must be left alone.
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "Parent complex membership must not change";
}

/*! An index in Molecule::emptyMolList is recycled instead of growing the list. */
void test_cmacr_create_recycles_empty_molecule_slot()
{
    std::cerr << "\n[TEST] test_cmacr_create_recycles_empty_molecule_slot\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      create_molecule_and_complex_from_rxn()\n"
              << "  Scenario:      moleculeList slot 1 is marked empty and pushed on\n"
              << "                 Molecule::emptyMolList before the call.\n"
              << "  Pass criteria: newMolIndex == 1, moleculeList does not grow, and\n"
              << "                 emptyMolList is consumed.\n";

    cmacr_ensure_rng();
    cmacr_reset_statics(1);

    Parameters params = cmacr_make_params();
    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();

    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };
    std::vector<ForwardRxn> forwardRxns {};
    CreateDestructRxn currRxn = cmacr_make_creation_rxn(3.0);

    // Slot 0: a live parent. Slot 1: a destroyed (empty) molecule to be reused.
    std::vector<Molecule> moleculeList {
        cmacr_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 0),
        cmacr_make_molecule(1, -1, 0, Coord { 5.0, 0.0, 0.0 }, Coord { 6.0, 0.0, 0.0 }, 0),
    };
    moleculeList[1].isEmpty = true;
    Molecule::emptyMolList.push_back(1);

    std::vector<Complex> complexList { cmacr_make_complex(
        0, Coord { 0.0, 0.0, 0.0 }, 2.0, 0) };

    const size_t molCountBefore = moleculeList.size();
    int newMolIndex = -1;
    int newComIndex = -1;

    std::cerr << "  Calling create_molecule_and_complex_from_rxn (recycling slot 1)...\n";
    create_molecule_and_complex_from_rxn(/*parentMolIndex=*/0, newMolIndex, newComIndex,
        /*createInVicinity=*/false, molTemplateList[0], params, currRxn, simulVolume,
        moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);

    std::cerr << "  newMolIndex = " << newMolIndex
              << ", moleculeList.size() = " << moleculeList.size() << '\n';

    EXPECT_EQ(newMolIndex, 1) << "The empty slot at index 1 should have been reused";
    EXPECT_EQ(moleculeList.size(), molCountBefore)
        << "Reusing an empty slot must not grow moleculeList";
    EXPECT_TRUE(Molecule::emptyMolList.empty())
        << "The consumed empty-slot index should be popped from emptyMolList";
    if (newMolIndex >= 0 && newMolIndex < static_cast<int>(moleculeList.size())) {
        EXPECT_FALSE(moleculeList[newMolIndex].isEmpty)
            << "The recycled slot now holds a live molecule";
        EXPECT_EQ(moleculeList[newMolIndex].myComIndex, newComIndex)
            << "Recycled molecule must point at the new complex";
    }
}

// =============================================================================
// MPI_create_molecule_and_complex_on_rank() tests
// =============================================================================

/*! The MPI variant only allocates and links slots: no coordinates are generated. */
void test_cmacr_mpi_create_on_rank()
{
    std::cerr << "\n[TEST] test_cmacr_mpi_create_on_rank\n"
              << "  Source file:   create_molecule_and_complex_from_rxn.cpp\n"
              << "  Function:      MPI_create_molecule_and_complex_on_rank()\n"
              << "  Scenario:      one empty molecule slot is available; no complex\n"
              << "                 slots are free.\n"
              << "  Pass criteria: the empty molecule slot is recycled, a complex is\n"
              << "                 appended, the two are cross-linked, trajStatus is\n"
              << "                 propagated, isGhosted is cleared and the global\n"
              << "                 complex counter is incremented.\n";

    cmacr_reset_statics(1);

    Membrane membraneObject = cmacr_make_membrane();
    SimulVolume simulVolume = cmacr_make_simul_volume();
    std::vector<MolTemplate> molTemplateList { cmacr_make_template(0, "A", false) };

    // Slot 0: a live molecule (also used as the unused `mol` argument).
    // Slot 1: an empty-but-typed slot which the routine should recycle.
    std::vector<Molecule> moleculeList {
        cmacr_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 0),
        cmacr_make_molecule(1, -1, 0, Coord { 4.0, 0.0, 0.0 }, Coord { 5.0, 0.0, 0.0 }, 0),
    };
    moleculeList[1].isEmpty = true;
    moleculeList[1].isGhosted = true; // must be cleared by the routine
    Molecule::emptyMolList.push_back(1);

    std::vector<Complex> complexList { cmacr_make_complex(
        0, Coord { 0.0, 0.0, 0.0 }, 2.0, 0) };

    const size_t comCountBefore = complexList.size();
    const int complexesBefore = Complex::numberOfComplexes;

    int newMolIndex = -1;
    int newComIndex = -1;

    std::cerr << "  Calling MPI_create_molecule_and_complex_on_rank...\n";
    MPI_create_molecule_and_complex_on_rank(moleculeList[0], newMolIndex, newComIndex,
        molTemplateList[0], simulVolume, moleculeList, complexList, molTemplateList,
        membraneObject);

    std::cerr << "  newMolIndex = " << newMolIndex << ", newComIndex = " << newComIndex
              << ", complexList.size() = " << complexList.size() << '\n';

    EXPECT_EQ(newMolIndex, 1) << "The empty molecule slot should have been recycled";
    EXPECT_EQ(newComIndex, static_cast<int>(comCountBefore))
        << "With no free complex slots a new complex must be appended";
    EXPECT_EQ(complexList.size(), comCountBefore + 1) << "complexList should grow by one";

    if (newMolIndex >= 0 && newMolIndex < static_cast<int>(moleculeList.size())
        && newComIndex >= 0 && newComIndex < static_cast<int>(complexList.size())) {
        const Molecule& newMol = moleculeList[newMolIndex];
        const Complex& newCom = complexList[newComIndex];

        EXPECT_EQ(newMol.myComIndex, newComIndex) << "Molecule must point at its new complex";
        EXPECT_EQ(static_cast<int>(newMol.trajStatus), static_cast<int>(TrajStatus::propagated))
            << "The molecule is flagged as propagated";
        EXPECT_FALSE(newMol.isGhosted) << "isGhosted must be cleared for a locally owned molecule";
        EXPECT_EQ(newMol.complexId, newCom.id) << "Molecule complexId must match the complex id";
        EXPECT_EQ(static_cast<int>(newCom.trajStatus), static_cast<int>(TrajStatus::propagated))
            << "The complex is flagged as propagated";
        EXPECT_EQ(newCom.memberList.size(), 1u) << "The new complex holds one member molecule";
        if (!newCom.memberList.empty()) {
            EXPECT_EQ(newCom.memberList[0], newMolIndex)
                << "The member should be the recycled molecule";
        }
    }

    EXPECT_EQ(Complex::numberOfComplexes, complexesBefore + 1)
        << "Complex::numberOfComplexes should be incremented once";
    EXPECT_TRUE(Molecule::emptyMolList.empty())
        << "The consumed empty-slot index should be popped from emptyMolList";
}

// =============================================================================
// GoogleTest wrappers - one per named test function so results are reported
// individually while every test still runs.
// =============================================================================
TEST(CreateMoleculeAndComplexFromRxn, OverlapInsideEmptyBin) { test_cmacr_overlap_inside_empty_bin(); }
TEST(CreateMoleculeAndComplexFromRxn, OverlapOutsideBox) { test_cmacr_overlap_outside_box(); }
TEST(CreateMoleculeAndComplexFromRxn, OverlapLipidMembraneCheck) { test_cmacr_overlap_lipid_membrane_check(); }
TEST(CreateMoleculeAndComplexFromRxn, OverlapReactivePartnerInsideBindRadius) { test_cmacr_overlap_reactive_partner_inside_bindradius(); }
TEST(CreateMoleculeAndComplexFromRxn, OverlapReactivePartnerOutsideBindRadius) { test_cmacr_overlap_reactive_partner_outside_bindradius(); }
TEST(CreateMoleculeAndComplexFromRxn, OverlapBoundingSphereShortcut) { test_cmacr_overlap_bounding_sphere_shortcut(); }
TEST(CreateMoleculeAndComplexFromRxn, CreateZerothOrder) { test_cmacr_create_zeroth_order(); }
TEST(CreateMoleculeAndComplexFromRxn, CreateInVicinity) { test_cmacr_create_in_vicinity(); }
TEST(CreateMoleculeAndComplexFromRxn, CreateRecyclesEmptyMoleculeSlot) { test_cmacr_create_recycles_empty_molecule_slot(); }
TEST(CreateMoleculeAndComplexFromRxn, MpiCreateOnRank) { test_cmacr_mpi_create_on_rank(); }