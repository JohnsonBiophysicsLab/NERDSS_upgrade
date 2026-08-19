/*! \file test_initialize_molecule_after_transmission_reaction.cpp
 *
 * ### Unit test for src/reactions/initialize_molecule_after_transmission_reaction.cpp
 *
 * Function under test:
 *
 *     Molecule initialize_molecule_after_transmission_reaction(
 *         int index, const Molecule& parentMol, Parameters& params,
 *         MolTemplate& molTemplate, const TransmissionRxn& currRxn,
 *         const Coord& newPos, const bool plusRand,
 *         const Membrane& membraneObject);
 *
 * The routine builds a brand new Molecule that has just been "transmitted"
 * across a compartment boundary. Its responsibilities are:
 *
 *   1. Copy identity data out of the MolTemplate (molTypeIndex, mass, isLipid).
 *   2. Create a freelist [0, 1, 2, ...] and a matching interfaceList.
 *   3. Place the center of mass at `newPos` (or, when `plusRand == true`, at a
 *      randomly oriented point exactly `molTemplate.radius` away from `newPos`
 *      that lies on the correct side of the compartment sphere).
 *   4. Position every interface at (template iCoord + new COM) and default all
 *      interface states to the first state in the template's stateList.
 *   5. Override interface states with the states demanded by the product of the
 *      TransmissionRxn, but only when the product's molTypeIndex matches.
 *   6. Bump the bookkeeping counters: Molecule::numberOfMolecules,
 *      params.numTotalUnits (+= nInterfaces + 1) and
 *      MolTemplate::numEachMolType[molTypeIndex].
 *
 * Each test below prints what is being exercised and what the pass criteria
 * are, then uses non-fatal EXPECT_* assertions so every test runs to
 * completion even when one fails.
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (all prefixed with imatr_ == Initialize Molecule After
// Transmission Reaction) so they cannot collide with other suite files.
// -----------------------------------------------------------------------------
namespace {

//! Molecule type index used by all the tests in this file.
constexpr int kImatrMolType = 0;

/*! \brief Make sure the static per-type counter vector is large enough.
 *
 * initialize_molecule_after_transmission_reaction() does
 * `++MolTemplate::numEachMolType[molTypeIndex]`, so the static vector must
 * already have an entry for our molecule type or we would index out of range.
 */
void imatr_ensure_static_storage(int molTypeIndex)
{
    if (static_cast<int>(MolTemplate::numEachMolType.size()) <= molTypeIndex)
        MolTemplate::numEachMolType.resize(molTypeIndex + 1, 0);
}

/*! \brief Build a small, fully-populated MolTemplate.
 *
 * \param[in] molTypeIndex index this template identifies as
 * \param[in] numIfaces    how many interfaces the template owns
 * \param[in] twoStates    if true each interface gets states 'U' (index 0) and
 *                         'P' (index 1); otherwise only 'U'
 */
MolTemplate imatr_make_moltemplate(int molTypeIndex, int numIfaces, bool twoStates)
{
    MolTemplate molTemplate;
    molTemplate.molName = "TRANSMOL";
    molTemplate.molTypeIndex = molTypeIndex;
    molTemplate.mass = 3.5;
    molTemplate.isLipid = false;
    molTemplate.radius = 2.0; // used as the random displacement magnitude
    molTemplate.comCoord = Coord(0.0, 0.0, 0.0);

    for (int ifaceItr = 0; ifaceItr < numIfaces; ++ifaceItr) {
        Interface iface;
        iface.index = ifaceItr;
        iface.name = "i" + std::to_string(ifaceItr);
        // Distinct, easy-to-verify internal coordinates for each interface.
        iface.iCoord = Coord(1.0 * (ifaceItr + 1), 0.5 * (ifaceItr + 1), -0.25 * (ifaceItr + 1));

        // State 0 is always 'U'; absolute index chosen to be recognizable.
        iface.stateList.push_back(Interface::State('U', 10 + 2 * ifaceItr));
        if (twoStates)
            iface.stateList.push_back(Interface::State('P', 11 + 2 * ifaceItr));

        molTemplate.interfaceList.push_back(iface);
    }
    return molTemplate;
}

/*! \brief Build a TransmissionRxn whose *product* demands one specific state.
 *
 * \param[in] productMolTypeIndex molTypeIndex stored on the product molecule.
 *            When this differs from the template's index the function under
 *            test must NOT override any interface state.
 * \param[in] relIfaceIndex which interface of the product is constrained
 * \param[in] requiresState the state character demanded for that interface
 */
TransmissionRxn imatr_make_rxn(int productMolTypeIndex, int relIfaceIndex, char requiresState)
{
    TransmissionRxn rxn;

    TransmissionRxn::TransmissionMol product;
    product.molTypeIndex = productMolTypeIndex;
    product.molName = "TRANSMOL";

    RxnIface rxnIface;
    rxnIface.ifaceName = "i" + std::to_string(relIfaceIndex);
    rxnIface.molTypeIndex = productMolTypeIndex;
    rxnIface.relIfaceIndex = relIfaceIndex;
    rxnIface.absIfaceIndex = 0;
    rxnIface.requiresState = requiresState;
    product.interfaceList.push_back(rxnIface);

    // The function only inspects productMolList.back(), but the list must be
    // non-empty for that to be legal.
    rxn.productMolList.push_back(product);
    return rxn;
}

/*! \brief A trivial "parent" molecule; the function ignores it entirely. */
Molecule imatr_make_parent(const Coord& com)
{
    Molecule parent;
    parent.comCoord = com;
    parent.index = 7;
    parent.molTypeIndex = kImatrMolType;
    return parent;
}

/*! \brief Lazily initialize the global GSL RNG (needed for plusRand == true). */
void imatr_init_rng_if_needed()
{
    if (r == nullptr) {
        std::cerr << "  (global gsl_rng* r was null -> calling srand_gsl(1))\n";
         const gsl_rng_type *T;
         T = gsl_rng_default;
         r = gsl_rng_alloc(T);
         gsl_rng_set(r, 42);
    }
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: deterministic placement (plusRand == false) and identity copy.
// -----------------------------------------------------------------------------
void test_imatr_basic_placement_and_identity()
{
    std::cerr << "\n[TEST] test_imatr_basic_placement_and_identity\n"
              << "  Source file: initialize_molecule_after_transmission_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_transmission_reaction\n"
              << "  Scenario:    plusRand == false, product molTypeIndex does not\n"
              << "               match the template (so no state override).\n"
              << "  Pass criteria: identity fields copied from the template, COM\n"
              << "               placed exactly at newPos, freelist == {0,1},\n"
              << "               interfaces placed at (iCoord + COM) with default\n"
              << "               state 0 information.\n";

    imatr_ensure_static_storage(kImatrMolType);

    Parameters params;
    params.numTotalUnits = 100;

    Membrane membraneObject;
    membraneObject.compartmentR = 10.0;

    MolTemplate molTemplate = imatr_make_moltemplate(kImatrMolType, 2, /*twoStates=*/true);
    // Product belongs to a *different* molecule type -> no state override.
    TransmissionRxn rxn = imatr_make_rxn(/*productMolTypeIndex=*/99, 0, 'P');

    const Coord parentPos { -20.0, -20.0, -20.0 };
    const Coord newPos { 5.0, -3.0, 2.0 };
    Molecule parent = imatr_make_parent(parentPos);

    std::cerr << "  Calling with index=42, newPos=(" << newPos.x << ", " << newPos.y
              << ", " << newPos.z << ")\n";
    Molecule mol = initialize_molecule_after_transmission_reaction(
        42, parent, params, molTemplate, rxn, newPos, /*plusRand=*/false, membraneObject);

    // ---- identity data lifted straight out of the template -------------------
    EXPECT_EQ(mol.molTypeIndex, molTemplate.molTypeIndex)
        << "molTypeIndex must be copied from the MolTemplate";
    EXPECT_DOUBLE_EQ(mol.mass, molTemplate.mass) << "mass must be copied from the MolTemplate";
    EXPECT_EQ(mol.isLipid, molTemplate.isLipid) << "isLipid must be copied from the MolTemplate";
    EXPECT_EQ(mol.index, 42) << "index must be the value passed in by the caller";
    EXPECT_FALSE(mol.isEmpty) << "a freshly created molecule must not be flagged empty";

    // ---- placement: exactly at newPos, NOT at the parent position ------------
    EXPECT_DOUBLE_EQ(mol.comCoord.x, newPos.x) << "COM.x should equal newPos.x when plusRand is false";
    EXPECT_DOUBLE_EQ(mol.comCoord.y, newPos.y) << "COM.y should equal newPos.y when plusRand is false";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, newPos.z) << "COM.z should equal newPos.z when plusRand is false";
    EXPECT_NE(mol.comCoord.x, parentPos.x)
        << "the parent molecule's coordinate must not be used for placement";

    // ---- freelist should be a simple iota over the interfaces ----------------
    ASSERT_EQ(mol.freelist.size(), molTemplate.interfaceList.size())
        << "freelist must have one entry per template interface";
    for (size_t i = 0; i < mol.freelist.size(); ++i) {
        EXPECT_EQ(mol.freelist[i], static_cast<int>(i))
            << "freelist should be filled with std::iota starting at 0";
    }

    // ---- interface geometry / default state information ---------------------
    ASSERT_EQ(mol.interfaceList.size(), molTemplate.interfaceList.size())
        << "interfaceList must have one entry per template interface";
    for (size_t i = 0; i < mol.interfaceList.size(); ++i) {
        const Coord expected = molTemplate.interfaceList[i].iCoord + mol.comCoord;
        std::cerr << "    iface " << i << " coord = (" << mol.interfaceList[i].coord.x << ", "
                  << mol.interfaceList[i].coord.y << ", " << mol.interfaceList[i].coord.z
                  << "), expected (" << expected.x << ", " << expected.y << ", " << expected.z << ")\n";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.x, expected.x) << "interface x = iCoord.x + COM.x";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.y, expected.y) << "interface y = iCoord.y + COM.y";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.z, expected.z) << "interface z = iCoord.z + COM.z";

        EXPECT_EQ(mol.interfaceList[i].relIndex, static_cast<int>(i))
            << "relIndex should be the interface's position in the list";
        EXPECT_EQ(mol.interfaceList[i].index, molTemplate.interfaceList[i].stateList[0].index)
            << "absolute index should come from stateList[0] of the template interface";
        EXPECT_EQ(mol.interfaceList[i].stateIden, 'U')
            << "default state identity should be stateList[0].iden ('U')";
        EXPECT_EQ(mol.interfaceList[i].stateIndex, 0)
            << "default stateIndex should be 0 when the reaction imposes nothing";
        EXPECT_EQ(mol.interfaceList[i].molTypeIndex, molTemplate.molTypeIndex)
            << "interface molTypeIndex should match the parent template";
    }
}

// -----------------------------------------------------------------------------
// Test 2: the bookkeeping counters must all advance by the right amount.
// -----------------------------------------------------------------------------
void test_imatr_counter_updates()
{
    std::cerr << "\n[TEST] test_imatr_counter_updates\n"
              << "  Source file: initialize_molecule_after_transmission_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_transmission_reaction\n"
              << "  Scenario:    create two molecules from a 3-interface template.\n"
              << "  Pass criteria: Molecule::numberOfMolecules += 1 per call,\n"
              << "               MolTemplate::numEachMolType[type] += 1 per call,\n"
              << "               params.numTotalUnits += (nIfaces + 1) per call.\n";

    imatr_ensure_static_storage(kImatrMolType);

    Parameters params;
    params.numTotalUnits = 0;

    Membrane membraneObject;
    membraneObject.compartmentR = 10.0;

    MolTemplate molTemplate = imatr_make_moltemplate(kImatrMolType, 3, /*twoStates=*/false);
    TransmissionRxn rxn = imatr_make_rxn(/*productMolTypeIndex=*/99, 0, 'U');
    Molecule parent = imatr_make_parent(Coord(0.0, 0.0, 0.0));
    const Coord newPos { 1.0, 1.0, 1.0 };

    const int molsBefore = Molecule::numberOfMolecules;
    const int typeBefore = MolTemplate::numEachMolType[kImatrMolType];
    const unsigned unitsBefore = params.numTotalUnits;
    const unsigned expectedUnitsPerCall
        = static_cast<unsigned>(molTemplate.interfaceList.size()) + 1u; // ifaces + COM

    std::cerr << "  Before: numberOfMolecules=" << molsBefore
              << " numEachMolType[" << kImatrMolType << "]=" << typeBefore
              << " numTotalUnits=" << unitsBefore << "\n";

    initialize_molecule_after_transmission_reaction(
        0, parent, params, molTemplate, rxn, newPos, false, membraneObject);
    initialize_molecule_after_transmission_reaction(
        1, parent, params, molTemplate, rxn, newPos, false, membraneObject);

    std::cerr << "  After : numberOfMolecules=" << Molecule::numberOfMolecules
              << " numEachMolType[" << kImatrMolType << "]="
              << MolTemplate::numEachMolType[kImatrMolType]
              << " numTotalUnits=" << params.numTotalUnits << "\n";

    EXPECT_EQ(Molecule::numberOfMolecules, molsBefore + 2)
        << "each call must increment the global molecule counter exactly once";
    EXPECT_EQ(MolTemplate::numEachMolType[kImatrMolType], typeBefore + 2)
        << "each call must increment the per-type molecule counter exactly once";
    EXPECT_EQ(params.numTotalUnits, unitsBefore + 2u * expectedUnitsPerCall)
        << "each call must add (numInterfaces + 1) units to params.numTotalUnits";
}

// -----------------------------------------------------------------------------
// Test 3: the product of the reaction overrides the interface state.
// -----------------------------------------------------------------------------
void test_imatr_product_state_override()
{
    std::cerr << "\n[TEST] test_imatr_product_state_override\n"
              << "  Source file: initialize_molecule_after_transmission_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_transmission_reaction\n"
              << "  Scenario:    product molTypeIndex MATCHES the template and the\n"
              << "               product requires state 'P' on interface 1.\n"
              << "  Pass criteria: interface 1 -> stateIden 'P', stateIndex 1;\n"
              << "               interface 0 keeps the default 'U'/0.\n";

    imatr_ensure_static_storage(kImatrMolType);

    Parameters params;
    Membrane membraneObject;
    membraneObject.compartmentR = 10.0;

    MolTemplate molTemplate = imatr_make_moltemplate(kImatrMolType, 2, /*twoStates=*/true);
    // Matching molTypeIndex -> the override branch must be taken.
    TransmissionRxn rxn = imatr_make_rxn(kImatrMolType, /*relIfaceIndex=*/1, /*requiresState=*/'P');
    Molecule parent = imatr_make_parent(Coord(0.0, 0.0, 0.0));

    Molecule mol = initialize_molecule_after_transmission_reaction(
        3, parent, params, molTemplate, rxn, Coord(0.0, 0.0, 0.0), false, membraneObject);

    ASSERT_EQ(mol.interfaceList.size(), 2u) << "template has two interfaces";

    std::cerr << "    iface0 state = '" << mol.interfaceList[0].stateIden << "' index "
              << mol.interfaceList[0].stateIndex << "\n"
              << "    iface1 state = '" << mol.interfaceList[1].stateIden << "' index "
              << mol.interfaceList[1].stateIndex << "\n";

    EXPECT_EQ(mol.interfaceList[1].stateIden, 'P')
        << "interface 1 identity should be overridden by the reaction product";
    EXPECT_EQ(mol.interfaceList[1].stateIndex, 1)
        << "interface 1 stateIndex should be the position of 'P' in the template stateList";

    EXPECT_EQ(mol.interfaceList[0].stateIden, 'U')
        << "interface 0 was not mentioned by the product and must keep state 'U'";
    EXPECT_EQ(mol.interfaceList[0].stateIndex, 0)
        << "interface 0 was not mentioned by the product and must keep stateIndex 0";
}

// -----------------------------------------------------------------------------
// Test 4: no override happens when the product is a different molecule type.
// -----------------------------------------------------------------------------
void test_imatr_no_override_for_other_moltype()
{
    std::cerr << "\n[TEST] test_imatr_no_override_for_other_moltype\n"
              << "  Source file: initialize_molecule_after_transmission_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_transmission_reaction\n"
              << "  Scenario:    product requires state 'P' but belongs to a\n"
              << "               different molTypeIndex than the template.\n"
              << "  Pass criteria: every interface keeps the default 'U'/0 state.\n";

    imatr_ensure_static_storage(kImatrMolType);

    Parameters params;
    Membrane membraneObject;
    membraneObject.compartmentR = 10.0;

    MolTemplate molTemplate = imatr_make_moltemplate(kImatrMolType, 2, /*twoStates=*/true);
    TransmissionRxn rxn = imatr_make_rxn(/*productMolTypeIndex=*/kImatrMolType + 5, 1, 'P');
    Molecule parent = imatr_make_parent(Coord(0.0, 0.0, 0.0));

    Molecule mol = initialize_molecule_after_transmission_reaction(
        4, parent, params, molTemplate, rxn, Coord(0.0, 0.0, 0.0), false, membraneObject);

    ASSERT_EQ(mol.interfaceList.size(), 2u) << "template has two interfaces";
    for (size_t i = 0; i < mol.interfaceList.size(); ++i) {
        EXPECT_EQ(mol.interfaceList[i].stateIden, 'U')
            << "interface " << i << " must keep the default state when molTypeIndex differs";
        EXPECT_EQ(mol.interfaceList[i].stateIndex, 0)
            << "interface " << i << " must keep stateIndex 0 when molTypeIndex differs";
    }
}

// -----------------------------------------------------------------------------
// Test 5: a required state that is absent from the template stateList only
//         changes the identity character, leaving stateIndex at its default.
// -----------------------------------------------------------------------------
void test_imatr_unknown_required_state()
{
    std::cerr << "\n[TEST] test_imatr_unknown_required_state\n"
              << "  Source file: initialize_molecule_after_transmission_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_transmission_reaction\n"
              << "  Scenario:    product requires state 'X' which does not exist in\n"
              << "               the template's stateList (only 'U' exists).\n"
              << "  Pass criteria: stateIden becomes 'X' but stateIndex stays 0\n"
              << "               because the search loop never finds a match.\n";

    imatr_ensure_static_storage(kImatrMolType);

    Parameters params;
    Membrane membraneObject;
    membraneObject.compartmentR = 10.0;

    // Only one state ('U') per interface here.
    MolTemplate molTemplate = imatr_make_moltemplate(kImatrMolType, 1, /*twoStates=*/false);
    TransmissionRxn rxn = imatr_make_rxn(kImatrMolType, /*relIfaceIndex=*/0, /*requiresState=*/'X');
    Molecule parent = imatr_make_parent(Coord(0.0, 0.0, 0.0));

    Molecule mol = initialize_molecule_after_transmission_reaction(
        5, parent, params, molTemplate, rxn, Coord(0.0, 0.0, 0.0), false, membraneObject);

    ASSERT_EQ(mol.interfaceList.size(), 1u) << "template has one interface";
    std::cerr << "    iface0 state = '" << mol.interfaceList[0].stateIden << "' index "
              << mol.interfaceList[0].stateIndex << "\n";

    EXPECT_EQ(mol.interfaceList[0].stateIden, 'X')
        << "the requested state character is assigned unconditionally";
    EXPECT_EQ(mol.interfaceList[0].stateIndex, 0)
        << "stateIndex should remain the default 0 when 'X' is not in the stateList";
}

// -----------------------------------------------------------------------------
// Test 6: plusRand == true, molecule must stay OUTSIDE the compartment sphere.
// -----------------------------------------------------------------------------
void test_imatr_plus_rand_outside_compartment()
{
    std::cerr << "\n[TEST] test_imatr_plus_rand_outside_compartment\n"
              << "  Source file: initialize_molecule_after_transmission_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_transmission_reaction\n"
              << "  Scenario:    plusRand == true, template flagged\n"
              << "               outsideCompartment, newPos far from the origin so\n"
              << "               the first random draw is always accepted.\n"
              << "  Pass criteria: |COM - newPos| == molTemplate.radius, the COM is\n"
              << "               outside the compartment radius, and interfaces are\n"
              << "               rebuilt relative to the displaced COM.\n";

    imatr_ensure_static_storage(kImatrMolType);
    imatr_init_rng_if_needed();

    Parameters params;
    Membrane membraneObject;
    membraneObject.compartmentR = 1.0; // tiny compartment, newPos is far outside

    MolTemplate molTemplate = imatr_make_moltemplate(kImatrMolType, 2, /*twoStates=*/false);
    molTemplate.outsideCompartment = true; // accepted when |COM| > compartmentR
    molTemplate.insideCompartment = false;

    TransmissionRxn rxn = imatr_make_rxn(/*productMolTypeIndex=*/99, 0, 'U');
    Molecule parent = imatr_make_parent(Coord(0.0, 0.0, 0.0));
    const Coord newPos { 50.0, 0.0, 0.0 };

    Molecule mol = initialize_molecule_after_transmission_reaction(
        6, parent, params, molTemplate, rxn, newPos, /*plusRand=*/true, membraneObject);

    const double dx = mol.comCoord.x - newPos.x;
    const double dy = mol.comCoord.y - newPos.y;
    const double dz = mol.comCoord.z - newPos.z;
    const double displacement = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double comRadius = std::sqrt(mol.comCoord.x * mol.comCoord.x
        + mol.comCoord.y * mol.comCoord.y + mol.comCoord.z * mol.comCoord.z);

    std::cerr << "    COM = (" << mol.comCoord.x << ", " << mol.comCoord.y << ", "
              << mol.comCoord.z << ")\n"
              << "    |COM - newPos| = " << displacement << " (expected "
              << molTemplate.radius << ")\n"
              << "    |COM| = " << comRadius << " (compartmentR = "
              << membraneObject.compartmentR << ")\n";

    EXPECT_NEAR(displacement, molTemplate.radius, 1e-9)
        << "the random displacement magnitude must equal the template radius";
    EXPECT_GT(comRadius, membraneObject.compartmentR)
        << "an outsideCompartment molecule must land outside the compartment sphere";

    // The interfaces must track the *displaced* COM, not the raw newPos.
    ASSERT_EQ(mol.interfaceList.size(), molTemplate.interfaceList.size())
        << "interfaceList size must follow the template";
    for (size_t i = 0; i < mol.interfaceList.size(); ++i) {
        const Coord expected = molTemplate.interfaceList[i].iCoord + mol.comCoord;
        EXPECT_NEAR(mol.interfaceList[i].coord.x, expected.x, 1e-12)
            << "interface " << i << " x must be relative to the displaced COM";
        EXPECT_NEAR(mol.interfaceList[i].coord.y, expected.y, 1e-12)
            << "interface " << i << " y must be relative to the displaced COM";
        EXPECT_NEAR(mol.interfaceList[i].coord.z, expected.z, 1e-12)
            << "interface " << i << " z must be relative to the displaced COM";
    }
}

// -----------------------------------------------------------------------------
// Test 7: plusRand == true, molecule must stay INSIDE the compartment sphere.
// -----------------------------------------------------------------------------
void test_imatr_plus_rand_inside_compartment()
{
    std::cerr << "\n[TEST] test_imatr_plus_rand_inside_compartment\n"
              << "  Source file: initialize_molecule_after_transmission_reaction.cpp\n"
              << "  Function:    initialize_molecule_after_transmission_reaction\n"
              << "  Scenario:    plusRand == true, template flagged\n"
              << "               insideCompartment, newPos at the origin of a large\n"
              << "               compartment so the first draw is always accepted.\n"
              << "  Pass criteria: |COM - newPos| == molTemplate.radius and the COM\n"
              << "               lies inside the compartment radius.\n";

    imatr_ensure_static_storage(kImatrMolType);
    imatr_init_rng_if_needed();

    Parameters params;
    Membrane membraneObject;
    membraneObject.compartmentR = 100.0; // huge compartment, newPos at its center

    MolTemplate molTemplate = imatr_make_moltemplate(kImatrMolType, 1, /*twoStates=*/false);
    molTemplate.outsideCompartment = false;
    molTemplate.insideCompartment = true; // accepted when |COM| < compartmentR

    TransmissionRxn rxn = imatr_make_rxn(/*productMolTypeIndex=*/99, 0, 'U');
    Molecule parent = imatr_make_parent(Coord(0.0, 0.0, 0.0));
    const Coord newPos { 0.0, 0.0, 0.0 };

    Molecule mol = initialize_molecule_after_transmission_reaction(
        7, parent, params, molTemplate, rxn, newPos, /*plusRand=*/true, membraneObject);

    const double displacement = std::sqrt(mol.comCoord.x * mol.comCoord.x
        + mol.comCoord.y * mol.comCoord.y + mol.comCoord.z * mol.comCoord.z);

    std::cerr << "    COM = (" << mol.comCoord.x << ", " << mol.comCoord.y << ", "
              << mol.comCoord.z << "), |COM| = " << displacement << "\n";

    EXPECT_NEAR(displacement, molTemplate.radius, 1e-9)
        << "the random displacement magnitude must equal the template radius";
    EXPECT_LT(displacement, membraneObject.compartmentR)
        << "an insideCompartment molecule must land inside the compartment sphere";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named helper runs inside its own TEST so that a
// failure in one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(InitializeMoleculeAfterTransmissionReaction, BasicPlacementAndIdentity)
{
    test_imatr_basic_placement_and_identity();
}
TEST(InitializeMoleculeAfterTransmissionReaction, CounterUpdates)
{
    test_imatr_counter_updates();
}
TEST(InitializeMoleculeAfterTransmissionReaction, ProductStateOverride)
{
    test_imatr_product_state_override();
}
TEST(InitializeMoleculeAfterTransmissionReaction, NoOverrideForOtherMolType)
{
    test_imatr_no_override_for_other_moltype();
}
TEST(InitializeMoleculeAfterTransmissionReaction, UnknownRequiredState)
{
    test_imatr_unknown_required_state();
}
TEST(InitializeMoleculeAfterTransmissionReaction, PlusRandOutsideCompartment)
{
    test_imatr_plus_rand_outside_compartment();
}
TEST(InitializeMoleculeAfterTransmissionReaction, PlusRandInsideCompartment)
{
    test_imatr_plus_rand_inside_compartment();
}