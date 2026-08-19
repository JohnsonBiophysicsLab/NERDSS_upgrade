/*! \file test_initialize_molecule_after_reaction.cpp
 *
 * ### Unit tests for src/reactions/initialize_molecule_after_reaction.cpp
 *
 * The file under test provides two factory functions which build a brand new
 * Molecule as the product of a creation reaction (CreateDestructRxn):
 *
 *   1. Molecule initialize_molecule_after_zeroth_reaction(
 *          int index, Parameters&, MolTemplate&, const CreateDestructRxn&,
 *          const Membrane&)
 *      -> "0 -> X" creation from concentration. The new molecule is given
 *         *random* coordinates inside the simulation volume.
 *
 *   2. Molecule initialize_molecule_after_uni_reaction(
 *          int index, const Molecule& parentMol, Parameters&, MolTemplate&,
 *          const CreateDestructRxn&)
 *      -> "X -> X + Y" creation from a parent molecule. The new molecule is
 *         placed on a sphere of radius currRxn.creationRadius centered on the
 *         parent molecule's center of mass.
 *
 * Both functions additionally:
 *   - copy identifying fields out of the MolTemplate (molTypeIndex, mass,
 *     isLipid, isPromoter),
 *   - build a freelist of every interface index (std::iota),
 *   - size the interfaceList to the template's interface count,
 *   - optionally override interface states using the reaction's product,
 *   - bump the global counters (Molecule::numberOfMolecules, Molecule::maxID,
 *     MolTemplate::numEachMolType) and params.numTotalUnits.
 *
 * The tests below verify each of those behaviours and print what is being
 * checked to stderr so a failure is easy to localise.
 */

#include "classes/class_Molecule_Complex.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with imar_ = Initialize Molecule After Reaction)
// -----------------------------------------------------------------------------
namespace {

/*! \brief Makes sure the GSL random number generator used by rand_gsl() exists.
 *
 * gtest_main.cpp defines `gsl_rng* r = nullptr;`, so the very first test that
 * needs randomness has to seed it. Seeding is done only once so repeated calls
 * from several tests are harmless.
 */
void imar_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] GSL rng was null -> calling srand_gsl(1)\n";
         const gsl_rng_type *T;
         T = gsl_rng_default;
         r = gsl_rng_alloc(T);
         gsl_rng_set(r, 42);
    }
}

/*! \brief Builds a small two-interface MolTemplate used by every test.
 *
 * Interface 0 ("a") has two states, 'U' (index 10) and 'P' (index 11).
 * Interface 1 ("b") has a single state 'X' (index 12).
 * The internal coordinates are deliberately simple so that the expected
 * absolute interface coordinates are trivial to compute by hand.
 */
MolTemplate imar_make_mol_template(int molTypeIndex)
{
    MolTemplate mt;
    mt.molTypeIndex = molTypeIndex;
    mt.molName = "TestMol";
    mt.mass = 42.0;
    mt.radius = 1.0;
    mt.isLipid = false; // keep it a solution molecule (simplest placement)
    mt.isPromoter = true; // non-default so we can verify it is copied over
    mt.isImplicitLipid = false;
    mt.isPoint = false;
    mt.isRod = false;
    mt.comCoord = Coord(0.0, 0.0, 0.0);
    mt.D = Coord(1.0, 1.0, 1.0);
    mt.Dr = Coord(0.01, 0.01, 0.01);

    // ---- interface 0: two possible states -----------------------------------
    Interface iface0;
    iface0.index = 0;
    iface0.name = "a";
    iface0.iCoord = Coord(1.0, 0.0, 0.0);
    iface0.stateList.emplace_back('U', 10);
    iface0.stateList.emplace_back('P', 11);

    // ---- interface 1: single state ------------------------------------------
    Interface iface1;
    iface1.index = 1;
    iface1.name = "b";
    iface1.iCoord = Coord(0.0, 1.0, 0.0);
    iface1.stateList.emplace_back('X', 12);

    mt.interfaceList.push_back(iface0);
    mt.interfaceList.push_back(iface1);

    // The functions under test increment MolTemplate::numEachMolType[molTypeIndex],
    // so that static vector must be large enough. Only ever grow it, so other
    // tests in the suite are not disturbed.
    if (MolTemplate::numEachMolType.size() < static_cast<size_t>(molTypeIndex) + 1)
        MolTemplate::numEachMolType.resize(molTypeIndex + 1, 0);

    return mt;
}

/*! \brief Builds a creation reaction whose product declares one interface state.
 *
 * \param[in] productMolTypeIndex molTypeIndex stored on the product molecule.
 * \param[in] relIfaceIndex       which interface of the product is declared.
 * \param[in] requiresState       the state character requested for it.
 * \param[in] creationRadius      distance used by the unimolecular variant.
 */
CreateDestructRxn imar_make_creation_rxn(
    int productMolTypeIndex, int relIfaceIndex, char requiresState, double creationRadius)
{
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::zerothOrderCreation;
    rxn.creationRadius = creationRadius;

    CreateDestructRxn::CreateDestructMol product;
    product.molTypeIndex = productMolTypeIndex;
    product.molName = "TestMol";
    product.interfaceList.push_back(
        RxnIface("a", productMolTypeIndex, /*absIfaceIndex*/ 11, relIfaceIndex, requiresState, false));

    rxn.productMolList.push_back(product);
    return rxn;
}

/*! \brief A cubic, reflecting, non-spherical water box centred on the origin. */
Membrane imar_make_box_membrane(double side)
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.implicitLipid = false;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { side, side, side });
    membraneObject.xBCtype = "reflect";
    membraneObject.yBCtype = "reflect";
    membraneObject.zBCtype = "reflect";
    return membraneObject;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: zeroth-order creation copies template fields and bumps all counters.
// -----------------------------------------------------------------------------
void test_imar_zeroth_basic_fields_and_counters()
{
    std::cerr << "\n[TEST] test_imar_zeroth_basic_fields_and_counters\n"
              << "  Source file: initialize_molecule_after_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_zeroth_reaction()\n"
              << "  Checks:      molTypeIndex/mass/isLipid/isPromoter copied,\n"
              << "               freelist = {0,1}, interfaceList sized to 2,\n"
              << "               isEmpty==false, index set, and the global\n"
              << "               counters/params.numTotalUnits are incremented.\n";

    imar_ensure_rng();

    MolTemplate mt = imar_make_mol_template(0);
    CreateDestructRxn rxn = imar_make_creation_rxn(0, 0, 'U', 1.0);
    Membrane membraneObject = imar_make_box_membrane(100.0);

    Parameters params;
    params.numTotalUnits = 5; // arbitrary starting value

    // Snapshot the globals so we can assert on the *delta* rather than absolute
    // values (other tests in the suite share these statics).
    const int molsBefore = Molecule::numberOfMolecules;
    const int idBefore = Molecule::maxID;
    const int typeCountBefore = MolTemplate::numEachMolType[0];

    std::cerr << "  Calling initialize_molecule_after_zeroth_reaction(index=7, ...)\n";
    Molecule mol = initialize_molecule_after_zeroth_reaction(7, params, mt, rxn, membraneObject);

    // --- identifying fields copied from the template -------------------------
    EXPECT_EQ(mol.index, 7) << "the index argument should be stored on the molecule";
    EXPECT_EQ(mol.molTypeIndex, 0) << "molTypeIndex should come from the MolTemplate";
    EXPECT_DOUBLE_EQ(mol.mass, 42.0) << "mass should come from the MolTemplate";
    EXPECT_FALSE(mol.isLipid) << "isLipid should mirror the MolTemplate (false here)";
    EXPECT_TRUE(mol.isPromoter) << "isPromoter should mirror the MolTemplate (true here)";
    EXPECT_FALSE(mol.isEmpty) << "a freshly created molecule must not be flagged empty";

    // --- freelist should be std::iota over every interface -------------------
    ASSERT_EQ(mol.freelist.size(), 2u) << "freelist must have one entry per interface";
    EXPECT_EQ(mol.freelist[0], 0) << "freelist should be filled by std::iota (0)";
    EXPECT_EQ(mol.freelist[1], 1) << "freelist should be filled by std::iota (1)";

    // --- interfaceList sized to match the template ---------------------------
    EXPECT_EQ(mol.interfaceList.size(), 2u) << "interfaceList must have one entry per template interface";

    // --- global bookkeeping --------------------------------------------------
    EXPECT_EQ(Molecule::numberOfMolecules, molsBefore + 1)
        << "Molecule::numberOfMolecules should increase by exactly one";
    EXPECT_EQ(mol.id, idBefore) << "the new molecule should take the previous maxID";
    EXPECT_EQ(Molecule::maxID, idBefore + 1) << "Molecule::maxID should be post-incremented";
    EXPECT_EQ(MolTemplate::numEachMolType[0], typeCountBefore + 1)
        << "the per-type molecule counter should increase by one";

    // numTotalUnits grows by (#interfaces + 1 for the COM) = 3 here.
    EXPECT_EQ(params.numTotalUnits, 8u)
        << "numTotalUnits should grow by interfaceList.size() + 1 (5 + 2 + 1 = 8)";

    std::cerr << "  Created COM = (" << mol.comCoord.x << ", " << mol.comCoord.y << ", "
              << mol.comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: zeroth-order creation places the molecule inside the water box.
// -----------------------------------------------------------------------------
void test_imar_zeroth_random_coords_inside_box()
{
    std::cerr << "\n[TEST] test_imar_zeroth_random_coords_inside_box\n"
              << "  Source file: initialize_molecule_after_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_zeroth_reaction()\n"
              << "  Checks:      create_random_coords() puts the COM inside the\n"
              << "               100x100x100 water box for several draws.\n"
              << "  Pass:        |x|,|y|,|z| <= 50 (box half-width) for each draw.\n";

    imar_ensure_rng();

    MolTemplate mt = imar_make_mol_template(0);
    CreateDestructRxn rxn = imar_make_creation_rxn(0, 0, 'U', 1.0);
    Membrane membraneObject = imar_make_box_membrane(100.0);
    const double half = membraneObject.waterBox.x / 2.0;

    Parameters params;
    params.numTotalUnits = 0;

    // Draw a handful of molecules; every one must land inside the box.
    for (int trial = 0; trial < 5; ++trial) {
        Molecule mol = initialize_molecule_after_zeroth_reaction(trial, params, mt, rxn, membraneObject);
        std::cerr << "  draw " << trial << ": COM = (" << mol.comCoord.x << ", " << mol.comCoord.y
                  << ", " << mol.comCoord.z << ")\n";
        EXPECT_LE(std::abs(mol.comCoord.x), half + 1e-6) << "x must be inside the water box";
        EXPECT_LE(std::abs(mol.comCoord.y), half + 1e-6) << "y must be inside the water box";
        EXPECT_LE(std::abs(mol.comCoord.z), half + 1e-6) << "z must be inside the water box";
    }
}

// -----------------------------------------------------------------------------
// Test 3: zeroth-order creation applies the product's requested state.
// -----------------------------------------------------------------------------
void test_imar_zeroth_applies_product_state()
{
    std::cerr << "\n[TEST] test_imar_zeroth_applies_product_state\n"
              << "  Source file: initialize_molecule_after_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_zeroth_reaction()\n"
              << "  Scenario:    the reaction product asks interface 0 to be in\n"
              << "               state 'P', which is stateList[1] on the template.\n"
              << "  Pass:        stateIden == 'P' and stateIndex == 1.\n";

    imar_ensure_rng();

    MolTemplate mt = imar_make_mol_template(0);
    // Product molTypeIndex matches the template, interface 0 requires state 'P'.
    CreateDestructRxn rxn = imar_make_creation_rxn(0, 0, 'P', 1.0);
    Membrane membraneObject = imar_make_box_membrane(100.0);

    Parameters params;
    params.numTotalUnits = 0;

    Molecule mol = initialize_molecule_after_zeroth_reaction(0, params, mt, rxn, membraneObject);

    ASSERT_EQ(mol.interfaceList.size(), 2u) << "need two interfaces to test state assignment";
    EXPECT_EQ(mol.interfaceList[0].stateIden, 'P')
        << "interface 0 should adopt the product's requested state character";
    EXPECT_EQ(mol.interfaceList[0].stateIndex, 1)
        << "stateIndex should be the position of 'P' in the template stateList (1)";

    std::cerr << "  interface0 stateIden='" << mol.interfaceList[0].stateIden
              << "', stateIndex=" << mol.interfaceList[0].stateIndex << "\n";
}

// -----------------------------------------------------------------------------
// Test 4: zeroth-order creation ignores a product whose molTypeIndex differs.
// -----------------------------------------------------------------------------
void test_imar_zeroth_ignores_mismatched_product()
{
    std::cerr << "\n[TEST] test_imar_zeroth_ignores_mismatched_product\n"
              << "  Source file: initialize_molecule_after_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_zeroth_reaction()\n"
              << "  Scenario:    the reaction product's molTypeIndex (1) does NOT\n"
              << "               match the template being created (0), so the\n"
              << "               state-assignment branch must be skipped.\n"
              << "  Pass:        interface 0 is not set to the requested 'Z'.\n";

    imar_ensure_rng();

    MolTemplate mt = imar_make_mol_template(0);
    // Product belongs to a *different* molecule type and asks for a state 'Z'
    // that does not even exist on the template -> must be ignored entirely.
    CreateDestructRxn rxn = imar_make_creation_rxn(1, 0, 'Z', 1.0);
    Membrane membraneObject = imar_make_box_membrane(100.0);

    Parameters params;
    params.numTotalUnits = 0;

    Molecule mol = initialize_molecule_after_zeroth_reaction(0, params, mt, rxn, membraneObject);

    ASSERT_EQ(mol.interfaceList.size(), 2u) << "interfaceList should still be sized from the template";
    EXPECT_NE(mol.interfaceList[0].stateIden, 'Z')
        << "a product of a different molTypeIndex must not set the state";

    std::cerr << "  interface0 stateIden='" << mol.interfaceList[0].stateIden
              << "' (expected: anything but 'Z')\n";
}

// -----------------------------------------------------------------------------
// Test 5: unimolecular creation places the product on a sphere of
//         creationRadius around the parent, and fills the interface fields.
// -----------------------------------------------------------------------------
void test_imar_uni_position_and_interfaces()
{
    std::cerr << "\n[TEST] test_imar_uni_position_and_interfaces\n"
              << "  Source file: initialize_molecule_after_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_uni_reaction()\n"
              << "  Checks:      |newCOM - parentCOM| == creationRadius, and every\n"
              << "               interface gets coord = iCoord + newCOM plus the\n"
              << "               correct index/relIndex/stateIden/stateIndex/molTypeIndex.\n";

    imar_ensure_rng();

    MolTemplate mt = imar_make_mol_template(0);
    const double creationRadius = 5.0;
    // Product molTypeIndex intentionally different so the default template
    // states (set in the interface loop) survive and can be checked here.
    CreateDestructRxn rxn = imar_make_creation_rxn(1, 0, 'P', creationRadius);

    Parameters params;
    params.numTotalUnits = 0;

    // Parent molecule sitting at an off-origin position.
    Molecule parentMol;
    parentMol.comCoord = Coord(10.0, 20.0, 30.0);

    std::cerr << "  Parent COM = (10, 20, 30), creationRadius = " << creationRadius << "\n";
    Molecule mol = initialize_molecule_after_uni_reaction(3, parentMol, params, mt, rxn);

    // --- placement: exactly creationRadius away from the parent --------------
    const double dx = mol.comCoord.x - parentMol.comCoord.x;
    const double dy = mol.comCoord.y - parentMol.comCoord.y;
    const double dz = mol.comCoord.z - parentMol.comCoord.z;
    const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    std::cerr << "  New COM    = (" << mol.comCoord.x << ", " << mol.comCoord.y << ", "
              << mol.comCoord.z << "), distance from parent = " << dist << "\n";
    EXPECT_NEAR(dist, creationRadius, 1e-8)
        << "the product must be placed on a sphere of radius creationRadius";

    // --- the parent must be untouched ---------------------------------------
    EXPECT_DOUBLE_EQ(parentMol.comCoord.x, 10.0) << "parent COM x must not change";
    EXPECT_DOUBLE_EQ(parentMol.comCoord.y, 20.0) << "parent COM y must not change";
    EXPECT_DOUBLE_EQ(parentMol.comCoord.z, 30.0) << "parent COM z must not change";

    // --- interface bookkeeping ----------------------------------------------
    ASSERT_EQ(mol.interfaceList.size(), 2u) << "one interface entry per template interface";

    // interface 0: iCoord (1,0,0), first state 'U' with absolute index 10
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.x, mol.comCoord.x + 1.0)
        << "iface0 x should be COM.x + iCoord.x";
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.y, mol.comCoord.y)
        << "iface0 y should be COM.y + 0";
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.z, mol.comCoord.z)
        << "iface0 z should be COM.z + 0";
    EXPECT_EQ(mol.interfaceList[0].index, 10) << "iface0 absolute index = stateList[0].index";
    EXPECT_EQ(mol.interfaceList[0].relIndex, 0) << "iface0 relative index should be 0";
    EXPECT_EQ(mol.interfaceList[0].stateIden, 'U') << "iface0 default state is stateList[0] ('U')";
    EXPECT_EQ(mol.interfaceList[0].stateIndex, 0) << "iface0 default stateIndex should be 0";
    EXPECT_EQ(mol.interfaceList[0].molTypeIndex, 0) << "iface0 molTypeIndex mirrors the template";

    // interface 1: iCoord (0,1,0), single state 'X' with absolute index 12
    EXPECT_DOUBLE_EQ(mol.interfaceList[1].coord.x, mol.comCoord.x)
        << "iface1 x should be COM.x + 0";
    EXPECT_DOUBLE_EQ(mol.interfaceList[1].coord.y, mol.comCoord.y + 1.0)
        << "iface1 y should be COM.y + iCoord.y";
    EXPECT_DOUBLE_EQ(mol.interfaceList[1].coord.z, mol.comCoord.z)
        << "iface1 z should be COM.z + 0";
    EXPECT_EQ(mol.interfaceList[1].index, 12) << "iface1 absolute index = stateList[0].index";
    EXPECT_EQ(mol.interfaceList[1].relIndex, 1) << "iface1 relative index should be 1";
    EXPECT_EQ(mol.interfaceList[1].stateIden, 'X') << "iface1 state should be 'X'";
    EXPECT_EQ(mol.interfaceList[1].stateIndex, 0) << "iface1 stateIndex should be 0";
    EXPECT_EQ(mol.interfaceList[1].molTypeIndex, 0) << "iface1 molTypeIndex mirrors the template";
}

// -----------------------------------------------------------------------------
// Test 6: unimolecular creation overrides interface states from the product.
// -----------------------------------------------------------------------------
void test_imar_uni_applies_product_state()
{
    std::cerr << "\n[TEST] test_imar_uni_applies_product_state\n"
              << "  Source file: initialize_molecule_after_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_uni_reaction()\n"
              << "  Scenario:    the product (same molTypeIndex) asks interface 0\n"
              << "               to be in state 'P' (stateList[1]).\n"
              << "  Pass:        iface0 becomes ('P', 1) while iface1 keeps ('X', 0).\n";

    imar_ensure_rng();

    MolTemplate mt = imar_make_mol_template(0);
    CreateDestructRxn rxn = imar_make_creation_rxn(0, 0, 'P', 2.0);

    Parameters params;
    params.numTotalUnits = 0;

    Molecule parentMol;
    parentMol.comCoord = Coord(0.0, 0.0, 0.0);

    Molecule mol = initialize_molecule_after_uni_reaction(1, parentMol, params, mt, rxn);

    ASSERT_EQ(mol.interfaceList.size(), 2u) << "need two interfaces for this check";
    EXPECT_EQ(mol.interfaceList[0].stateIden, 'P')
        << "iface0 should be overridden with the product's state 'P'";
    EXPECT_EQ(mol.interfaceList[0].stateIndex, 1)
        << "iface0 stateIndex should point at 'P' (position 1)";
    EXPECT_EQ(mol.interfaceList[1].stateIden, 'X')
        << "iface1 was not declared by the product and keeps its default 'X'";
    EXPECT_EQ(mol.interfaceList[1].stateIndex, 0)
        << "iface1 stateIndex should remain the default 0";

    std::cerr << "  iface0=('" << mol.interfaceList[0].stateIden << "',"
              << mol.interfaceList[0].stateIndex << ")  iface1=('"
              << mol.interfaceList[1].stateIden << "," << mol.interfaceList[1].stateIndex << ")\n";
}

// -----------------------------------------------------------------------------
// Test 7: unimolecular creation performs the same bookkeeping as the zeroth
//         order variant (counters, freelist, flags, index/id).
// -----------------------------------------------------------------------------
void test_imar_uni_counters_and_flags()
{
    std::cerr << "\n[TEST] test_imar_uni_counters_and_flags\n"
              << "  Source file: initialize_molecule_after_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_uni_reaction()\n"
              << "  Checks:      numberOfMolecules/maxID/numEachMolType and\n"
              << "               params.numTotalUnits are all updated, freelist\n"
              << "               is iota, isEmpty is false and index/id are set.\n";

    imar_ensure_rng();

    MolTemplate mt = imar_make_mol_template(0);
    CreateDestructRxn rxn = imar_make_creation_rxn(0, 0, 'U', 3.0);

    Parameters params;
    params.numTotalUnits = 100;

    Molecule parentMol;
    parentMol.comCoord = Coord(-4.0, 6.5, 0.25);

    const int molsBefore = Molecule::numberOfMolecules;
    const int idBefore = Molecule::maxID;
    const int typeCountBefore = MolTemplate::numEachMolType[0];

    Molecule mol = initialize_molecule_after_uni_reaction(11, parentMol, params, mt, rxn);

    EXPECT_EQ(mol.index, 11) << "the index argument must be stored on the molecule";
    EXPECT_FALSE(mol.isEmpty) << "a newly created molecule is not empty";
    EXPECT_TRUE(mol.isPromoter) << "isPromoter should be copied from the MolTemplate";
    EXPECT_FALSE(mol.isLipid) << "isLipid should be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(mol.mass, 42.0) << "mass should be copied from the MolTemplate";

    ASSERT_EQ(mol.freelist.size(), 2u) << "freelist should have one entry per interface";
    EXPECT_EQ(mol.freelist[0], 0) << "freelist[0] should be 0 (std::iota)";
    EXPECT_EQ(mol.freelist[1], 1) << "freelist[1] should be 1 (std::iota)";

    EXPECT_EQ(Molecule::numberOfMolecules, molsBefore + 1)
        << "Molecule::numberOfMolecules must increase by one";
    EXPECT_EQ(mol.id, idBefore) << "the new molecule takes the previous maxID";
    EXPECT_EQ(Molecule::maxID, idBefore + 1) << "Molecule::maxID must be post-incremented";
    EXPECT_EQ(MolTemplate::numEachMolType[0], typeCountBefore + 1)
        << "the per-type counter must increase by one";
    EXPECT_EQ(params.numTotalUnits, 103u)
        << "numTotalUnits should grow by interfaces + 1 (100 + 2 + 1 = 103)";

    std::cerr << "  numberOfMolecules " << molsBefore << " -> " << Molecule::numberOfMolecules
              << ", maxID " << idBefore << " -> " << Molecule::maxID << "\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each named helper runs inside its own TEST so that a
// single failure does not stop the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(InitializeMoleculeAfterReaction, ZerothBasicFieldsAndCounters)
{
    test_imar_zeroth_basic_fields_and_counters();
}
TEST(InitializeMoleculeAfterReaction, ZerothRandomCoordsInsideBox)
{
    test_imar_zeroth_random_coords_inside_box();
}
TEST(InitializeMoleculeAfterReaction, ZerothAppliesProductState)
{
    test_imar_zeroth_applies_product_state();
}
TEST(InitializeMoleculeAfterReaction, ZerothIgnoresMismatchedProduct)
{
    test_imar_zeroth_ignores_mismatched_product();
}
TEST(InitializeMoleculeAfterReaction, UniPositionAndInterfaces)
{
    test_imar_uni_position_and_interfaces();
}
TEST(InitializeMoleculeAfterReaction, UniAppliesProductState)
{
    test_imar_uni_applies_product_state();
}
TEST(InitializeMoleculeAfterReaction, UniCountersAndFlags)
{
    test_imar_uni_counters_and_flags();
}