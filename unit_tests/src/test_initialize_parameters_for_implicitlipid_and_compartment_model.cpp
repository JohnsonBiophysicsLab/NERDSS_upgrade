/*! \file test_initialize_parameters_for_implicitlipid_and_compartment_model.cpp
 *
 * ### Unit tests for
 *     src/system_setup/initialize_parameters_for_implicitlipid_and_compartment_model.cpp
 *
 * The file under test contains two setup routines:
 *
 *   1. initialize_paramters_for_implicitlipid_and_compartment_model(...)
 *   2. initialize_paramters_for_implicitlipid_model(...)
 *
 * Both walk the molecule / template / reaction lists and fill in the
 * Membrane bookkeeping needed by the implicit-lipid (and, for the first
 * routine, the compartment/transmission) machinery:
 *
 *   - membraneObject.implicitLipid / lipidLength / TwoD
 *   - implicitlipidIndex (and membraneObject.implicitlipidIndex)
 *   - Molecule::isImplicitLipid flags
 *   - MolTemplate::bindToSurface flags
 *   - membraneObject.totalSA (box vs sphere)
 *   - membraneObject.nSites, numberOfFreeLipidsEachState,
 *     numberOfProteinEachState
 *   - the 500-entry look-up table membraneObject.RS3Dvect, laid out as
 *         [i]      sigma (bindRadius)
 *         [i+100]  ka    (rate)
 *         [i+200]  Dtot
 *         [i+300]  RS3D = sigma*ka*2 / (ka*2 + 4*pi*sigma*Dtot)
 *         [i+400]  molTypeIndex of the non-implicit-lipid partner
 *
 * The one behavioural difference exercised below: the *combined* routine only
 * allocates RS3Dvect when there is an implicit lipid **or** a compartment,
 * whereas the implicit-lipid-only routine always allocates it.
 *
 * NOTE: these routines index into moleculeList/molTemplateList without bounds
 * checks, so every object built here is fully initialised (interfaces, states,
 * rate lists, per-state counter vectors sized to Membrane::nStates).  No test
 * below drives the code down an exit()/abort() path.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "system_setup/system_setup.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Small builders shared by the tests below.  All names carry the "ilcm_" prefix
// (Implicit Lipid / Compartment Model) so they cannot collide with helpers from
// other translation units in the combined gtest binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a MolTemplate carrying exactly one interface with one state.
 *
 * \param[in] name           molecule name (also used to name the interface)
 * \param[in] molTypeIndex   index of this template inside molTemplateList
 * \param[in] absIfaceIndex  absolute index of the single interface state
 * \param[in] isImplicitLipid  flags this template as the implicit lipid
 * \param[in] diffusion      translational diffusion constants
 * \param[in] radius         template radius (becomes Membrane::lipidLength for IL)
 * \param[in] copies         starting copy number (becomes Membrane::nSites for IL)
 */
MolTemplate ilcm_make_template(const std::string& name, int molTypeIndex, int absIfaceIndex,
    bool isImplicitLipid, const Coord& diffusion, double radius, int copies)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = molTypeIndex;
    temp.copies = copies;
    temp.radius = radius;
    temp.mass = 1.0;
    temp.isImplicitLipid = isImplicitLipid;
    temp.isLipid = isImplicitLipid; // an implicit lipid is always a lipid
    temp.D = diffusion;
    temp.Dr = Coord { 0.0, 0.0, 0.0 };
    temp.bindToSurface = false;

    // One interface, one state.  The state's absolute index is what the
    // reaction reactant lists are compared against.
    Interface iface;
    iface.name = name + "_iface";
    iface.index = 0;
    Interface::State state;
    state.index = absIfaceIndex;
    state.iden = '\0';
    iface.stateList.push_back(state);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a single-interface Molecule that matches ilcm_make_template(). */
Molecule ilcm_make_molecule(int index, int molTypeIndex, int absIfaceIndex)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = index;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.comCoord = Coord { 0.0, 0.0, 0.0 };

    Molecule::Iface iface;
    iface.coord = mol.comCoord;
    iface.index = absIfaceIndex; // absolute index of the current state
    iface.relIndex = 0;
    iface.stateIndex = 0; // index into the template's stateList
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a two-reactant bimolecular ForwardRxn with a single rate. */
ForwardRxn ilcm_make_bimolecular_rxn(
    int molType1, int absIface1, int molType2, int absIface2, double bindRadius, double rate)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.reactantListNew.emplace_back("r1", molType1, absIface1, 0, '\0', false);
    rxn.reactantListNew.emplace_back("r2", molType2, absIface2, 0, '\0', false);
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = rate;
    return rxn;
}

/*! \brief Build a one-reactant TransmissionRxn with a single rate. */
TransmissionRxn ilcm_make_transmission_rxn(
    int molTypeIndex, int absIfaceIndex, double bindRadius, double rate)
{
    TransmissionRxn rxn;
    rxn.rxnType = ReactionType::transmission;
    rxn.bindRadius = bindRadius;
    rxn.reactantListNew.emplace_back("t1", molTypeIndex, absIfaceIndex, 0, '\0', false);
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = rate;
    return rxn;
}

/*! \brief Reference implementation of the RS3D closed form used by the code. */
double ilcm_expected_RS3D(double sigma, double ka, double Dtot)
{
    return sigma * ka * 2.0 / (ka * 2.0 + 4 * M_PI * sigma * Dtot);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: combined routine, box geometry, one implicit lipid + one protein.
//         Checks every field the routine is responsible for.
// -----------------------------------------------------------------------------
void test_ilcm_combined_with_implicit_lipid_box()
{
    std::cerr << "\n[TEST] test_ilcm_combined_with_implicit_lipid_box\n"
              << "  Source file: initialize_parameters_for_implicitlipid_and_compartment_model.cpp\n"
              << "  Function:    initialize_paramters_for_implicitlipid_and_compartment_model()\n"
              << "  Scenario:    rectangular water box, molecule 0 is the implicit lipid,\n"
              << "               molecule 1 is a soluble protein that binds it.\n"
              << "  Pass:        IL flags/index set, bindToSurface set on the protein,\n"
              << "               totalSA = x*y, nSites = IL copies, per-state counters\n"
              << "               filled, and RS3Dvect[0/100/200/300/400] populated.\n";

    // --- templates: 0 = implicit lipid (immobile), 1 = protein (D = 10 in x,y,z)
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(ilcm_make_template(
        "IL", 0, /*absIface*/ 0, /*isIL*/ true, Coord { 0.0, 0.0, 0.0 }, /*radius*/ 1.0, /*copies*/ 100));
    molTemplateList.push_back(ilcm_make_template(
        "A", 1, /*absIface*/ 1, /*isIL*/ false, Coord { 10.0, 10.0, 10.0 }, /*radius*/ 2.0, /*copies*/ 5));

    // --- one bimolecular reaction IL(iface,abs 0) + A(iface,abs 1)
    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(ilcm_make_bimolecular_rxn(
        /*molType1*/ 0, /*abs1*/ 0, /*molType2*/ 1, /*abs2*/ 1, /*sigma*/ 1.0, /*ka*/ 10.0));
    // The protein's interface state must point at that reaction so the
    // "number of proteins per implicit-lipid state" counter can be built.
    molTemplateList[1].interfaceList[0].stateList[0].myForwardRxns.push_back(0);

    std::vector<BackRxn> backRxns;
    std::vector<TransmissionRxn> transmissionRxns;

    // --- molecules: index 0 is the implicit lipid, index 1 the protein
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(ilcm_make_molecule(0, 0, 0));
    moleculeList.push_back(ilcm_make_molecule(1, 1, 1));

    std::vector<Complex> complexList; // unused by the routine, kept empty

    // --- membrane: 100 x 100 x 100 box, one implicit-lipid state
    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    membraneObject.isSphere = false;
    membraneObject.hasCompartment = false;
    membraneObject.nStates = 1;
    membraneObject.nSites = 0;
    membraneObject.numberOfFreeLipidsEachState.assign(membraneObject.nStates, 0);
    membraneObject.numberOfProteinEachState.assign(membraneObject.nStates, 0);

    Parameters params;
    params.fromRestart = false;

    int implicitlipidIndex { -1 };

    std::cerr << "  Calling initialize_paramters_for_implicitlipid_and_compartment_model...\n";
    initialize_paramters_for_implicitlipid_and_compartment_model(implicitlipidIndex, params, forwardRxns,
        backRxns, transmissionRxns, moleculeList, molTemplateList, complexList, membraneObject);

    // --- implicit-lipid discovery -------------------------------------------
    EXPECT_TRUE(membraneObject.implicitLipid) << "membrane should be flagged as having an implicit lipid";
    EXPECT_DOUBLE_EQ(membraneObject.lipidLength, 1.0) << "lipidLength should equal the IL template radius";
    EXPECT_EQ(implicitlipidIndex, 0) << "implicitlipidIndex should be the moleculeList index of the IL";
    EXPECT_EQ(membraneObject.implicitlipidIndex, 0) << "membraneObject.implicitlipidIndex should match";
    EXPECT_TRUE(moleculeList[0].isImplicitLipid) << "molecule 0 is the implicit lipid";
    EXPECT_FALSE(moleculeList[1].isImplicitLipid) << "molecule 1 is a normal protein";

    // The protein diffuses in z, so the system is NOT flagged as 2D.
    EXPECT_FALSE(membraneObject.TwoD) << "a protein with D.z > 0 in a z>0 box means the system is 3D";

    // --- bindToSurface ------------------------------------------------------
    EXPECT_TRUE(molTemplateList[1].bindToSurface)
        << "protein template must be flagged bindToSurface (it reacts with the IL)";
    EXPECT_FALSE(molTemplateList[0].bindToSurface)
        << "the implicit lipid template itself should not be flagged bindToSurface";

    // --- surface area (box branch) ------------------------------------------
    EXPECT_DOUBLE_EQ(membraneObject.totalSA, 100.0 * 100.0) << "box totalSA should be waterBox.x * waterBox.y";

    // --- copy-number bookkeeping --------------------------------------------
    EXPECT_EQ(membraneObject.nSites, 100) << "nSites should be the implicit lipid's copy number";
    ASSERT_FALSE(membraneObject.numberOfFreeLipidsEachState.empty());
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 100)
        << "all lipids start free in state 0 for a non-restart run";
    ASSERT_FALSE(membraneObject.numberOfProteinEachState.empty());
    EXPECT_EQ(membraneObject.numberOfProteinEachState[0], 1)
        << "exactly one protein interface can bind implicit-lipid state 0";

    // --- RS3D lookup table --------------------------------------------------
    ASSERT_EQ(membraneObject.RS3Dvect.size(), 500u) << "RS3Dvect must be allocated with 500 slots";

    const double sigma = 1.0;
    const double ka = 10.0;
    // Dtot = (1/3)(Dx1+Dx2) + (1/3)(Dy1+Dy2) + (1/3)(Dz1+Dz2)
    const double expectedDtot = 1.0 / 3.0 * (0.0 + 10.0) + 1.0 / 3.0 * (0.0 + 10.0) + 1.0 / 3.0 * (0.0 + 10.0);
    const double expectedRS3D = ilcm_expected_RS3D(sigma, ka, expectedDtot);

    std::cerr << "  RS3Dvect[0]=" << membraneObject.RS3Dvect[0] << " [100]=" << membraneObject.RS3Dvect[100]
              << " [200]=" << membraneObject.RS3Dvect[200] << " [300]=" << membraneObject.RS3Dvect[300]
              << " [400]=" << membraneObject.RS3Dvect[400] << '\n';

    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[0], sigma) << "slot 0 stores the binding radius";
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[100], ka) << "slot 100 stores the association rate";
    EXPECT_NEAR(membraneObject.RS3Dvect[200], expectedDtot, 1e-12) << "slot 200 stores the total diffusion";
    EXPECT_NEAR(membraneObject.RS3Dvect[300], expectedRS3D, 1e-12) << "slot 300 stores the RS3D closed form";
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[400], 1.0)
        << "slot 400 stores the molTypeIndex of the non-implicit-lipid partner";

    // Only one entry was written; the next slot of each block stays at -1.
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[1], -1.0) << "unused RS3D slots stay at the -1 sentinel";
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[401], -1.0) << "unused RS3D slots stay at the -1 sentinel";
}

// -----------------------------------------------------------------------------
// Test 2: combined routine with neither implicit lipid nor compartment.
//         RS3Dvect must stay completely unallocated.
// -----------------------------------------------------------------------------
void test_ilcm_combined_without_il_or_compartment()
{
    std::cerr << "\n[TEST] test_ilcm_combined_without_il_or_compartment\n"
              << "  Function:    initialize_paramters_for_implicitlipid_and_compartment_model()\n"
              << "  Scenario:    two soluble proteins, no implicit lipid, no compartment.\n"
              << "  Pass:        RS3Dvect is left empty, no IL flags set, totalSA still computed.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(
        ilcm_make_template("A", 0, 0, /*isIL*/ false, Coord { 5.0, 5.0, 5.0 }, 1.0, 10));
    molTemplateList.push_back(
        ilcm_make_template("B", 1, 1, /*isIL*/ false, Coord { 3.0, 3.0, 3.0 }, 1.0, 10));

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(ilcm_make_bimolecular_rxn(0, 0, 1, 1, 1.5, 20.0));

    std::vector<BackRxn> backRxns;
    std::vector<TransmissionRxn> transmissionRxns;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(ilcm_make_molecule(0, 0, 0));
    moleculeList.push_back(ilcm_make_molecule(1, 1, 1));

    std::vector<Complex> complexList;

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 50.0, 40.0, 30.0 });
    membraneObject.isSphere = false;
    membraneObject.hasCompartment = false;
    // No implicit-lipid states exist, so the RS3D scan over states is a no-op.
    membraneObject.nStates = 0;
    membraneObject.nSites = 0;

    Parameters params;
    params.fromRestart = false;

    int implicitlipidIndex { -1 };

    std::cerr << "  Calling initialize_paramters_for_implicitlipid_and_compartment_model...\n";
    initialize_paramters_for_implicitlipid_and_compartment_model(implicitlipidIndex, params, forwardRxns,
        backRxns, transmissionRxns, moleculeList, molTemplateList, complexList, membraneObject);

    EXPECT_FALSE(membraneObject.implicitLipid) << "no implicit lipid present";
    EXPECT_EQ(implicitlipidIndex, -1) << "implicitlipidIndex must be left untouched";
    EXPECT_EQ(membraneObject.implicitlipidIndex, -1) << "membrane IL index must be left at its default";
    EXPECT_FALSE(moleculeList[0].isImplicitLipid) << "molecule 0 is not an implicit lipid";
    EXPECT_FALSE(moleculeList[1].isImplicitLipid) << "molecule 1 is not an implicit lipid";

    EXPECT_FALSE(molTemplateList[0].bindToSurface) << "no reactant is an implicit lipid, so no surface binding";
    EXPECT_FALSE(molTemplateList[1].bindToSurface) << "no reactant is an implicit lipid, so no surface binding";

    EXPECT_DOUBLE_EQ(membraneObject.totalSA, 50.0 * 40.0) << "box totalSA should still be computed";

    // This is the distinguishing behaviour of the combined routine.
    EXPECT_TRUE(membraneObject.RS3Dvect.empty())
        << "RS3Dvect is only allocated when an implicit lipid or a compartment exists";

    std::cerr << "  RS3Dvect size = " << membraneObject.RS3Dvect.size() << " (expected 0)\n";
}

// -----------------------------------------------------------------------------
// Test 3: 2D detection (every mobile species has D.z == 0) and the spherical
//         surface-area branch.
// -----------------------------------------------------------------------------
void test_ilcm_two_d_detection_and_sphere_area()
{
    std::cerr << "\n[TEST] test_ilcm_two_d_detection_and_sphere_area\n"
              << "  Function:    initialize_paramters_for_implicitlipid_and_compartment_model()\n"
              << "  Scenario:    every molecule has D.z == 0 and the boundary is a sphere.\n"
              << "  Pass:        membraneObject.TwoD becomes true and\n"
              << "               totalSA = 4*pi*R^2 instead of x*y.\n";

    std::vector<MolTemplate> molTemplateList;
    // D.z == 0 for both templates -> the routine's "is2D" flag stays true.
    molTemplateList.push_back(
        ilcm_make_template("A", 0, 0, /*isIL*/ false, Coord { 1.0, 1.0, 0.0 }, 1.0, 4));
    molTemplateList.push_back(
        ilcm_make_template("B", 1, 1, /*isIL*/ false, Coord { 2.0, 2.0, 0.0 }, 1.0, 4));

    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<TransmissionRxn> transmissionRxns;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(ilcm_make_molecule(0, 0, 0));
    moleculeList.push_back(ilcm_make_molecule(1, 1, 1));

    std::vector<Complex> complexList;

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 20.0, 20.0, 20.0 });
    membraneObject.isSphere = true;
    membraneObject.sphereR = 50.0;
    membraneObject.hasCompartment = false;
    membraneObject.nStates = 0;
    membraneObject.nSites = 0;

    Parameters params;
    params.fromRestart = false;

    int implicitlipidIndex { -1 };

    std::cerr << "  Calling initialize_paramters_for_implicitlipid_and_compartment_model...\n";
    initialize_paramters_for_implicitlipid_and_compartment_model(implicitlipidIndex, params, forwardRxns,
        backRxns, transmissionRxns, moleculeList, molTemplateList, complexList, membraneObject);

    EXPECT_TRUE(membraneObject.TwoD) << "all species confined to z => system flagged as 2D";

    const double expectedSA = 4.0 * M_PI * std::pow(50.0, 2.0);
    std::cerr << "  totalSA = " << membraneObject.totalSA << " (expected " << expectedSA << ")\n";
    EXPECT_NEAR(membraneObject.totalSA, expectedSA, 1e-9) << "sphere totalSA should be 4*pi*R^2";
}

// -----------------------------------------------------------------------------
// Test 4: compartment-only path.  No implicit lipid, but hasCompartment == true
//         and a transmission reaction, so RS3Dvect is built from that reaction.
// -----------------------------------------------------------------------------
void test_ilcm_compartment_transmission_table()
{
    std::cerr << "\n[TEST] test_ilcm_compartment_transmission_table\n"
              << "  Function:    initialize_paramters_for_implicitlipid_and_compartment_model()\n"
              << "  Scenario:    no implicit lipid, hasCompartment = true, one transmission rxn.\n"
              << "  Pass:        RS3Dvect allocated (500 slots) and filled from the\n"
              << "               transmission reaction using the droplet diffusion constant.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(
        ilcm_make_template("A", 0, 0, /*isIL*/ false, Coord { 2.0, 2.0, 2.0 }, 1.0, 7));

    std::vector<ForwardRxn> forwardRxns; // none needed for this path
    std::vector<BackRxn> backRxns;

    const double sigma = 3.0;
    const double ka = 4.0;
    std::vector<TransmissionRxn> transmissionRxns;
    transmissionRxns.push_back(ilcm_make_transmission_rxn(/*molTypeIndex*/ 0, /*absIface*/ 0, sigma, ka));

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(ilcm_make_molecule(0, 0, 0));

    std::vector<Complex> complexList;

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 80.0, 60.0, 40.0 });
    membraneObject.isSphere = false;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = 10.0;
    membraneObject.droplet.D = 5.0; // diffusion of the compartment surface sites
    // No implicit-lipid states -> the forward-reaction RS3D scan is skipped
    // entirely, so no implicit-lipid molecule is dereferenced.
    membraneObject.nStates = 0;
    membraneObject.nSites = 0;

    Parameters params;
    params.fromRestart = false;

    int implicitlipidIndex { -1 };

    std::cerr << "  Calling initialize_paramters_for_implicitlipid_and_compartment_model...\n";
    initialize_paramters_for_implicitlipid_and_compartment_model(implicitlipidIndex, params, forwardRxns,
        backRxns, transmissionRxns, moleculeList, molTemplateList, complexList, membraneObject);

    ASSERT_EQ(membraneObject.RS3Dvect.size(), 500u)
        << "a compartment alone is enough to allocate the RS3D table";

    // Dtot = (1/3)(Dx + dropletD) + (1/3)(Dy + dropletD) + (1/3)(Dz + 0)
    const double expectedDtot
        = 1.0 / 3.0 * (2.0 + 5.0) + 1.0 / 3.0 * (2.0 + 5.0) + 1.0 / 3.0 * (2.0 + 0.0);
    const double expectedRS3D = ilcm_expected_RS3D(sigma, ka, expectedDtot);

    std::cerr << "  RS3Dvect[0]=" << membraneObject.RS3Dvect[0] << " [100]=" << membraneObject.RS3Dvect[100]
              << " [200]=" << membraneObject.RS3Dvect[200] << " [300]=" << membraneObject.RS3Dvect[300]
              << " [400]=" << membraneObject.RS3Dvect[400] << '\n';

    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[0], sigma) << "transmission binding radius stored at slot 0";
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[100], ka) << "transmission rate stored at slot 100";
    EXPECT_NEAR(membraneObject.RS3Dvect[200], expectedDtot, 1e-12)
        << "Dtot mixes the molecule D with the droplet D (z term uses 0)";
    EXPECT_NEAR(membraneObject.RS3Dvect[300], expectedRS3D, 1e-12) << "RS3D closed form for the compartment";
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[400], 0.0)
        << "slot 400 stores the reactant molTypeIndex of the transmission reaction";

    // No implicit lipid was present, so nothing IL-related may be flagged.
    EXPECT_FALSE(membraneObject.implicitLipid) << "no implicit lipid in a compartment-only setup";
    EXPECT_EQ(implicitlipidIndex, -1) << "implicitlipidIndex untouched without an implicit lipid";
}

// -----------------------------------------------------------------------------
// Test 5: implicit-lipid-only routine with NO implicit lipid.
//         Unlike the combined routine, RS3Dvect is still allocated (all -1).
// -----------------------------------------------------------------------------
void test_ilcm_il_only_allocates_table_without_il()
{
    std::cerr << "\n[TEST] test_ilcm_il_only_allocates_table_without_il\n"
              << "  Function:    initialize_paramters_for_implicitlipid_model()\n"
              << "  Scenario:    two soluble proteins, no implicit lipid.\n"
              << "  Pass:        RS3Dvect still has 500 slots, every one at the -1 sentinel.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(
        ilcm_make_template("A", 0, 0, /*isIL*/ false, Coord { 6.0, 6.0, 6.0 }, 1.0, 3));
    molTemplateList.push_back(
        ilcm_make_template("B", 1, 1, /*isIL*/ false, Coord { 6.0, 6.0, 6.0 }, 1.0, 3));

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(ilcm_make_bimolecular_rxn(0, 0, 1, 1, 2.0, 8.0));
    std::vector<BackRxn> backRxns;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(ilcm_make_molecule(0, 0, 0));
    moleculeList.push_back(ilcm_make_molecule(1, 1, 1));

    std::vector<Complex> complexList;

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 10.0, 10.0, 10.0 });
    membraneObject.isSphere = false;
    membraneObject.nStates = 0; // no implicit-lipid states exist
    membraneObject.nSites = 0;

    Parameters params;
    params.fromRestart = false;

    int implicitlipidIndex { -1 };

    std::cerr << "  Calling initialize_paramters_for_implicitlipid_model...\n";
    initialize_paramters_for_implicitlipid_model(implicitlipidIndex, params, forwardRxns, backRxns,
        moleculeList, molTemplateList, complexList, membraneObject);

    ASSERT_EQ(membraneObject.RS3Dvect.size(), 500u)
        << "the implicit-lipid-only routine always allocates the RS3D table";

    // Every slot must remain at the -1 sentinel because nothing was written.
    bool allSentinel = true;
    for (std::size_t i = 0; i < membraneObject.RS3Dvect.size(); ++i) {
        if (membraneObject.RS3Dvect[i] != -1.0) {
            allSentinel = false;
            std::cerr << "  Unexpected non-sentinel value at index " << i << ": "
                      << membraneObject.RS3Dvect[i] << '\n';
            break;
        }
    }
    EXPECT_TRUE(allSentinel) << "with no implicit lipid every RS3D slot should still be -1";

    EXPECT_FALSE(membraneObject.implicitLipid) << "no implicit lipid should be detected";
    EXPECT_DOUBLE_EQ(membraneObject.totalSA, 10.0 * 10.0) << "box totalSA should be x*y";
}

// -----------------------------------------------------------------------------
// Test 6: implicit-lipid-only routine WITH an implicit lipid.  Should produce
//         the same bookkeeping as the combined routine for that scenario.
// -----------------------------------------------------------------------------
void test_ilcm_il_only_with_implicit_lipid()
{
    std::cerr << "\n[TEST] test_ilcm_il_only_with_implicit_lipid\n"
              << "  Function:    initialize_paramters_for_implicitlipid_model()\n"
              << "  Scenario:    identical IL + protein system as the combined-routine test.\n"
              << "  Pass:        same IL flags, nSites, per-state counters and RS3D entries.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(ilcm_make_template(
        "IL", 0, /*absIface*/ 0, /*isIL*/ true, Coord { 0.0, 0.0, 0.0 }, /*radius*/ 1.5, /*copies*/ 250));
    molTemplateList.push_back(ilcm_make_template(
        "A", 1, /*absIface*/ 1, /*isIL*/ false, Coord { 4.0, 4.0, 4.0 }, /*radius*/ 2.0, /*copies*/ 2));

    const double sigma = 2.0;
    const double ka = 12.0;
    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(ilcm_make_bimolecular_rxn(0, 0, 1, 1, sigma, ka));
    molTemplateList[1].interfaceList[0].stateList[0].myForwardRxns.push_back(0);

    std::vector<BackRxn> backRxns;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(ilcm_make_molecule(0, 0, 0)); // implicit lipid
    moleculeList.push_back(ilcm_make_molecule(1, 1, 1)); // protein
    moleculeList.push_back(ilcm_make_molecule(2, 1, 1)); // a second protein

    std::vector<Complex> complexList;

    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 30.0, 20.0, 10.0 });
    membraneObject.isSphere = false;
    membraneObject.nStates = 1;
    membraneObject.nSites = 0;
    membraneObject.numberOfFreeLipidsEachState.assign(membraneObject.nStates, 0);
    membraneObject.numberOfProteinEachState.assign(membraneObject.nStates, 0);

    Parameters params;
    params.fromRestart = false;

    int implicitlipidIndex { -1 };

    std::cerr << "  Calling initialize_paramters_for_implicitlipid_model...\n";
    initialize_paramters_for_implicitlipid_model(implicitlipidIndex, params, forwardRxns, backRxns,
        moleculeList, molTemplateList, complexList, membraneObject);

    EXPECT_TRUE(membraneObject.implicitLipid) << "implicit lipid must be detected";
    EXPECT_DOUBLE_EQ(membraneObject.lipidLength, 1.5) << "lipidLength comes from the IL template radius";
    EXPECT_EQ(implicitlipidIndex, 0) << "the IL sits at moleculeList index 0";
    EXPECT_TRUE(moleculeList[0].isImplicitLipid) << "molecule 0 flagged as the implicit lipid";
    EXPECT_FALSE(moleculeList[1].isImplicitLipid) << "molecule 1 is a protein";
    EXPECT_FALSE(moleculeList[2].isImplicitLipid) << "molecule 2 is a protein";

    EXPECT_TRUE(molTemplateList[1].bindToSurface) << "the protein template binds the surface";

    EXPECT_EQ(membraneObject.nSites, 250) << "nSites is the implicit lipid copy number";
    ASSERT_FALSE(membraneObject.numberOfFreeLipidsEachState.empty());
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 250) << "all lipids start free in state 0";
    ASSERT_FALSE(membraneObject.numberOfProteinEachState.empty());
    EXPECT_EQ(membraneObject.numberOfProteinEachState[0], 2)
        << "both protein molecules contribute one bindable interface each";

    ASSERT_EQ(membraneObject.RS3Dvect.size(), 500u) << "RS3D table allocated";

    const double expectedDtot
        = 1.0 / 3.0 * (0.0 + 4.0) + 1.0 / 3.0 * (0.0 + 4.0) + 1.0 / 3.0 * (0.0 + 4.0);
    const double expectedRS3D = ilcm_expected_RS3D(sigma, ka, expectedDtot);

    std::cerr << "  RS3Dvect[0]=" << membraneObject.RS3Dvect[0] << " [100]=" << membraneObject.RS3Dvect[100]
              << " [200]=" << membraneObject.RS3Dvect[200] << " [300]=" << membraneObject.RS3Dvect[300]
              << " [400]=" << membraneObject.RS3Dvect[400] << '\n';

    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[0], sigma) << "binding radius stored at slot 0";
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[100], ka) << "association rate stored at slot 100";
    EXPECT_NEAR(membraneObject.RS3Dvect[200], expectedDtot, 1e-12) << "total diffusion stored at slot 200";
    EXPECT_NEAR(membraneObject.RS3Dvect[300], expectedRS3D, 1e-12) << "RS3D stored at slot 300";
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[400], 1.0) << "partner molTypeIndex stored at slot 400";
    EXPECT_DOUBLE_EQ(membraneObject.RS3Dvect[1], -1.0) << "only one RS3D entry should have been written";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so a failure
// in one scenario does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(InitializeImplicitLipidCompartmentModel, CombinedWithImplicitLipidBox)
{
    test_ilcm_combined_with_implicit_lipid_box();
}

TEST(InitializeImplicitLipidCompartmentModel, CombinedWithoutIlOrCompartment)
{
    test_ilcm_combined_without_il_or_compartment();
}

TEST(InitializeImplicitLipidCompartmentModel, TwoDDetectionAndSphereArea)
{
    test_ilcm_two_d_detection_and_sphere_area();
}

TEST(InitializeImplicitLipidCompartmentModel, CompartmentTransmissionTable)
{
    test_ilcm_compartment_transmission_table();
}

TEST(InitializeImplicitLipidModel, AllocatesTableWithoutImplicitLipid)
{
    test_ilcm_il_only_allocates_table_without_il();
}

TEST(InitializeImplicitLipidModel, WithImplicitLipid) { test_ilcm_il_only_with_implicit_lipid(); }