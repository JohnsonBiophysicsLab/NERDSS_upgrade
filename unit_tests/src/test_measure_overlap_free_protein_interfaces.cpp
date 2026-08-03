/*! \file test_measure_overlap_free_protein_interfaces.cpp
 *
 * ### Unit test for src/reactions/measure_overlap_free_protein_interfaces.cpp
 *
 * Function under test:
 *
 *     void measure_overlap_free_protein_interfaces(Molecule base1, Molecule baseTmp,
 *                                                  bool& flagCancel,
 *                                                  const std::vector<MolTemplate>& molTemplateList,
 *                                                  const std::vector<ForwardRxn>& forwardRxns,
 *                                                  const std::vector<BackRxn>& backRxns)
 *
 * Semantics of the routine (as documented in the source file):
 *   - `base1`   is a molecule that already lives in the system (real coordinates are used).
 *   - `baseTmp` is a molecule that is currently performing association, so its *temporary*
 *               coordinates (`tmpICoords`) are the ones that must be inspected.
 *   - The routine only does work if the two molecule types are declared reaction partners
 *     (via MolTemplate::rxnPartners) and `baseTmp` is not an implicit lipid.
 *   - For every pair of *free* interfaces that are declared reaction partners and for which
 *     find_which_reaction() locates a valid reaction, the separation between
 *     `base1.interfaceList[i].coord` and `baseTmp.tmpICoords[j]` is compared against the
 *     reaction's binding radius.  If the separation is strictly smaller than the binding
 *     radius, `flagCancel` is set to true (association will be cancelled).
 *   - `flagCancel` is never reset to false by the routine.
 *
 * The tests below construct a minimal but self-consistent world:
 *   molType 0 ("A") with one interface, absolute interface index 0
 *   molType 1 ("B") with one interface, absolute interface index 1
 *   one bimolecular ForwardRxn  A(a) + B(b) -> A(a!1).B(b!1)  with a known binding radius.
 *
 * Verbose progress information is printed to stderr so a reader of the test log can follow
 * exactly which source file / function is being exercised and what each check verifies.
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "classes/class_Rxns.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

namespace {

//! Absolute interface index of molType 0's single interface.
constexpr int kMofpiAbsIfaceA = 0;
//! Absolute interface index of molType 1's single interface.
constexpr int kMofpiAbsIfaceB = 1;
//! Binding radius used by the single test reaction (nm).
constexpr double kMofpiBindRadius = 2.0;

/*! \brief Build a single MolTemplate with exactly one interface and one state.
 *
 * \param[in] name             Molecule name (only used for readability of failure output).
 * \param[in] molTypeIndex     Index of this template in molTemplateList.
 * \param[in] absIfaceIndex    Absolute index of the single interface's single state.
 * \param[in] rxnPartners      Molecule-type indices this template can react with.
 * \param[in] stateRxnPartners Absolute interface indices the state can react with.
 * \param[in] myForwardRxns    Indices into forwardRxns in which the state participates.
 */
MolTemplate mofpi_make_template(const std::string& name, int molTypeIndex, int absIfaceIndex,
    const std::vector<int>& rxnPartners, const std::vector<unsigned>& stateRxnPartners,
    const std::vector<unsigned>& myForwardRxns)
{
    MolTemplate molTemplate;
    molTemplate.molName = name;
    molTemplate.molTypeIndex = molTypeIndex;
    molTemplate.copies = 1;
    molTemplate.radius = 1.0;
    molTemplate.mass = 1.0;
    molTemplate.D = Coord(1.0, 1.0, 1.0);
    molTemplate.Dr = Coord(0.01, 0.01, 0.01);
    molTemplate.isLipid = false;
    molTemplate.isImplicitLipid = false;
    molTemplate.isPoint = false;
    molTemplate.isRod = false;
    molTemplate.rxnPartners = rxnPartners;

    // Build the single interface with its single (default) state.
    Interface iface;
    iface.index = 0; // relative index within the template
    iface.name = "iface";
    iface.iCoord = Coord(0.0, 0.0, 0.0);

    Interface::State state;
    state.ifaceAndStateName = "iface";
    state.iden = '\0'; // no explicit state identifier
    state.index = absIfaceIndex; // absolute index of this state
    state.myForwardRxns = myForwardRxns;
    state.rxnPartners = stateRxnPartners;
    iface.stateList.push_back(state);

    molTemplate.interfaceList.push_back(iface);
    return molTemplate;
}

/*! \brief Build the two-template list used by all tests.
 *
 * \param[in] aPartners       rxnPartners of template 0 (allows switching interaction on/off).
 * \param[in] aStatePartners  rxnPartners of template 0's state.
 */
std::vector<MolTemplate> mofpi_make_template_list(
    const std::vector<int>& aPartners, const std::vector<unsigned>& aStatePartners)
{
    std::vector<MolTemplate> molTemplateList;
    // Template 0: "A", interface absolute index 0, participates in forward reaction 0.
    molTemplateList.push_back(
        mofpi_make_template("A", 0, kMofpiAbsIfaceA, aPartners, aStatePartners, std::vector<unsigned> { 0 }));
    // Template 1: "B", interface absolute index 1, participates in forward reaction 0.
    molTemplateList.push_back(mofpi_make_template("B", 1, kMofpiAbsIfaceB, std::vector<int> { 0 },
        std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceA) }, std::vector<unsigned> { 0 }));
    return molTemplateList;
}

/*! \brief Build the single bimolecular ForwardRxn A(a) + B(b) -> A(a!1).B(b!1). */
std::vector<ForwardRxn> mofpi_make_forward_rxns(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = false;
    rxn.isSymmetric = false;
    rxn.isOnMem = false;
    rxn.hasStateChange = false;
    rxn.isCoupled = false;
    rxn.bindRadius = bindRadius;
    rxn.bindRadius2D = bindRadius;
    rxn.rxnLabel = "testRxn";
    rxn.productName = "A(a!1).B(b!1)";

    // Reactants: (name, molTypeIndex, absIfaceIndex, relIfaceIndex, requiresState, requiresInteraction)
    rxn.reactantListNew.emplace_back("A(a)", 0, kMofpiAbsIfaceA, 0, '\0', false);
    rxn.reactantListNew.emplace_back("B(b)", 1, kMofpiAbsIfaceB, 0, '\0', false);
    // Products: the same interfaces, now bound (absolute indices 2 and 3).
    rxn.productListNew.emplace_back("A(a!1)", 0, 2, 0, '\0', true);
    rxn.productListNew.emplace_back("B(b!1)", 1, 3, 0, '\0', true);

    rxn.intReactantList = { kMofpiAbsIfaceA, kMofpiAbsIfaceB };
    rxn.intProductList = { 2, 3 };

    // A single rate with no ancillary ("other") interface requirements.
    rxn.rateList.emplace_back(1.0, std::vector<std::vector<RxnIface>> {});

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(rxn);
    return forwardRxns;
}

/*! \brief Build a minimal Molecule with a single free interface.
 *
 * \param[in] molTypeIndex  Index of the molecule's MolTemplate.
 * \param[in] absIfaceIndex Absolute index of the molecule's single interface state.
 * \param[in] realCoord     "Real" (system) coordinate of the COM and the interface.
 * \param[in] tmpCoord      Temporary (association) coordinate of the COM and the interface.
 */
Molecule mofpi_make_molecule(int molTypeIndex, int absIfaceIndex, const Coord& realCoord, const Coord& tmpCoord)
{
    Molecule mol;
    mol.index = molTypeIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isEmpty = false;
    mol.comCoord = realCoord;
    mol.tmpComCoord = tmpCoord;

    Molecule::Iface iface;
    iface.coord = realCoord;
    iface.stateIden = '\0';
    iface.stateIndex = 0; // position inside MolTemplate's stateList
    iface.index = absIfaceIndex; // absolute state index
    iface.relIndex = 0;
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    // The single interface is free (available for binding).
    mol.freelist.push_back(0);

    // Temporary association coordinate for the single interface.
    mol.tmpICoords.push_back(tmpCoord);

    return mol;
}

/*! \brief Convenience: Euclidean distance between two Coords (for log messages). */
double mofpi_distance(const Coord& c1, const Coord& c2)
{
    const double dx = c1.x - c2.x;
    const double dy = c1.y - c2.y;
    const double dz = c1.z - c2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: molecule types are not reaction partners -> routine must do nothing.
// -----------------------------------------------------------------------------
void test_mofpi_non_interacting_types()
{
    std::cerr << "\n[TEST] test_mofpi_non_interacting_types\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario:    template 0 lists NO reaction partners, so the two molecules\n"
              << "               cannot interact even though their interfaces overlap.\n"
              << "  Criteria:    flagCancel must remain false (early exit before any distance\n"
              << "               calculation).\n";

    // Template 0 declares no reaction partners -> canInteract == false.
    std::vector<MolTemplate> molTemplateList
        = mofpi_make_template_list(std::vector<int> {}, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    // Deliberately place the interfaces right on top of one another.
    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    Molecule baseTmp = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(0.1, 0.0, 0.0), Coord(0.1, 0.0, 0.0));

    std::cerr << "  Interface separation (tmp) = " << mofpi_distance(base1.interfaceList[0].coord, baseTmp.tmpICoords[0])
              << " nm, bindRadius = " << kMofpiBindRadius << " nm\n";

    bool flagCancel = false;
    measure_overlap_free_protein_interfaces(base1, baseTmp, flagCancel, molTemplateList, forwardRxns, backRxns);

    EXPECT_FALSE(flagCancel) << "Non-partner molecule types must never cancel association";
    std::cerr << "  flagCancel after call = " << std::boolalpha << flagCancel << "\n";
}

// -----------------------------------------------------------------------------
// Test 2: the associating molecule is an implicit lipid -> routine must do nothing.
// -----------------------------------------------------------------------------
void test_mofpi_implicit_lipid_partner()
{
    std::cerr << "\n[TEST] test_mofpi_implicit_lipid_partner\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario:    baseTmp.isImplicitLipid == true (partner is the implicit\n"
              << "               membrane, which has no real interface coordinates).\n"
              << "  Criteria:    flagCancel must remain false.\n";

    std::vector<MolTemplate> molTemplateList = mofpi_make_template_list(
        std::vector<int> { 1 }, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    Molecule baseTmp = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(0.1, 0.0, 0.0), Coord(0.1, 0.0, 0.0));
    baseTmp.isImplicitLipid = true; // this is the property being tested

    bool flagCancel = false;
    measure_overlap_free_protein_interfaces(base1, baseTmp, flagCancel, molTemplateList, forwardRxns, backRxns);

    EXPECT_FALSE(flagCancel) << "Implicit-lipid partners must be skipped, so no cancellation";
    std::cerr << "  flagCancel after call = " << std::boolalpha << flagCancel << "\n";
}

// -----------------------------------------------------------------------------
// Test 3: interacting types but interfaces far apart -> no cancellation.
// -----------------------------------------------------------------------------
void test_mofpi_far_apart_no_cancel()
{
    std::cerr << "\n[TEST] test_mofpi_far_apart_no_cancel\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario:    the two reactive interfaces are 20 nm apart, far outside the\n"
              << "               2 nm binding radius.\n"
              << "  Criteria:    flagCancel must remain false.\n";

    std::vector<MolTemplate> molTemplateList = mofpi_make_template_list(
        std::vector<int> { 1 }, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    Molecule baseTmp = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(20.0, 0.0, 0.0), Coord(20.0, 0.0, 0.0));

    std::cerr << "  Interface separation (tmp) = " << mofpi_distance(base1.interfaceList[0].coord, baseTmp.tmpICoords[0])
              << " nm, bindRadius = " << kMofpiBindRadius << " nm\n";

    bool flagCancel = false;
    measure_overlap_free_protein_interfaces(base1, baseTmp, flagCancel, molTemplateList, forwardRxns, backRxns);

    EXPECT_FALSE(flagCancel) << "Well-separated interfaces must not cancel association";
    std::cerr << "  flagCancel after call = " << std::boolalpha << flagCancel << "\n";
}

// -----------------------------------------------------------------------------
// Test 4: interacting types and overlapping interfaces -> cancellation expected.
// -----------------------------------------------------------------------------
void test_mofpi_overlap_cancels()
{
    std::cerr << "\n[TEST] test_mofpi_overlap_cancels\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario:    partner interfaces separated by only 0.5 nm, well inside the\n"
              << "               2 nm binding radius of the single test reaction.\n"
              << "  Criteria:    flagCancel must be set to true (association cancelled).\n";

    std::vector<MolTemplate> molTemplateList = mofpi_make_template_list(
        std::vector<int> { 1 }, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    Molecule baseTmp = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(0.5, 0.0, 0.0), Coord(0.5, 0.0, 0.0));

    const double sep = mofpi_distance(base1.interfaceList[0].coord, baseTmp.tmpICoords[0]);
    std::cerr << "  Interface separation (tmp) = " << sep << " nm, bindRadius = " << kMofpiBindRadius << " nm\n";

    bool flagCancel = false;
    measure_overlap_free_protein_interfaces(base1, baseTmp, flagCancel, molTemplateList, forwardRxns, backRxns);

    EXPECT_TRUE(flagCancel) << "Overlapping free interfaces (separation " << sep << " nm < bindRadius "
                            << kMofpiBindRadius << " nm) must cancel association";
    std::cerr << "  flagCancel after call = " << std::boolalpha << flagCancel << "\n";
}

// -----------------------------------------------------------------------------
// Test 5: exactly at the binding radius -> strict "<" comparison means no cancel.
// -----------------------------------------------------------------------------
void test_mofpi_exactly_at_bindradius()
{
    std::cerr << "\n[TEST] test_mofpi_exactly_at_bindradius\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario:    separation is exactly the binding radius (2.0 nm) so that\n"
              << "               d2 == bindrad2 exactly in floating point.\n"
              << "  Criteria:    the source uses a strict 'd2 < bindrad2' test, therefore\n"
              << "               flagCancel must remain false.\n";

    std::vector<MolTemplate> molTemplateList = mofpi_make_template_list(
        std::vector<int> { 1 }, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    // 2.0 nm along x: both 2.0 and 4.0 are exactly representable in binary floating point.
    Molecule baseTmp = mofpi_make_molecule(
        1, kMofpiAbsIfaceB, Coord(kMofpiBindRadius, 0.0, 0.0), Coord(kMofpiBindRadius, 0.0, 0.0));

    std::cerr << "  Interface separation (tmp) = " << mofpi_distance(base1.interfaceList[0].coord, baseTmp.tmpICoords[0])
              << " nm, bindRadius = " << kMofpiBindRadius << " nm\n";

    bool flagCancel = false;
    measure_overlap_free_protein_interfaces(base1, baseTmp, flagCancel, molTemplateList, forwardRxns, backRxns);

    EXPECT_FALSE(flagCancel) << "Separation exactly equal to the binding radius must not cancel";
    std::cerr << "  flagCancel after call = " << std::boolalpha << flagCancel << "\n";
}

// -----------------------------------------------------------------------------
// Test 6: the routine must inspect baseTmp's *temporary* coordinates, not its real ones.
// -----------------------------------------------------------------------------
void test_mofpi_uses_tmp_coords_of_second_molecule()
{
    std::cerr << "\n[TEST] test_mofpi_uses_tmp_coords_of_second_molecule\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario A:  baseTmp real coord far away (30 nm) but tmp coord overlapping\n"
              << "               (0.25 nm) -> the tmp coord must be the one used.\n"
              << "  Scenario B:  the mirror image: real coord overlapping, tmp coord far away.\n"
              << "  Criteria:    A cancels, B does not.\n";

    std::vector<MolTemplate> molTemplateList = mofpi_make_template_list(
        std::vector<int> { 1 }, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));

    // --- Scenario A: real far / tmp close -------------------------------------
    Molecule tmpClose = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(30.0, 0.0, 0.0), Coord(0.25, 0.0, 0.0));
    std::cerr << "  Scenario A: real separation = "
              << mofpi_distance(base1.interfaceList[0].coord, tmpClose.interfaceList[0].coord) << " nm, tmp separation = "
              << mofpi_distance(base1.interfaceList[0].coord, tmpClose.tmpICoords[0]) << " nm\n";

    bool flagCancelA = false;
    measure_overlap_free_protein_interfaces(base1, tmpClose, flagCancelA, molTemplateList, forwardRxns, backRxns);
    EXPECT_TRUE(flagCancelA) << "Temporary coordinates of the associating molecule must be used (overlap missed)";
    std::cerr << "    flagCancel (A) = " << std::boolalpha << flagCancelA << "\n";

    // --- Scenario B: real close / tmp far -------------------------------------
    Molecule tmpFar = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(0.25, 0.0, 0.0), Coord(30.0, 0.0, 0.0));
    std::cerr << "  Scenario B: real separation = "
              << mofpi_distance(base1.interfaceList[0].coord, tmpFar.interfaceList[0].coord) << " nm, tmp separation = "
              << mofpi_distance(base1.interfaceList[0].coord, tmpFar.tmpICoords[0]) << " nm\n";

    bool flagCancelB = false;
    measure_overlap_free_protein_interfaces(base1, tmpFar, flagCancelB, molTemplateList, forwardRxns, backRxns);
    EXPECT_FALSE(flagCancelB) << "Real (non-temporary) coordinates of the associating molecule must be ignored";
    std::cerr << "    flagCancel (B) = " << std::boolalpha << flagCancelB << "\n";
}

// -----------------------------------------------------------------------------
// Test 7: no free interfaces -> the interface loops never execute.
// -----------------------------------------------------------------------------
void test_mofpi_empty_freelists()
{
    std::cerr << "\n[TEST] test_mofpi_empty_freelists\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario 1:  base1 has an empty freelist (all interfaces already bound).\n"
              << "  Scenario 2:  baseTmp has an empty freelist.\n"
              << "  Criteria:    flagCancel must remain false in both cases even though the\n"
              << "               interfaces physically overlap.\n";

    std::vector<MolTemplate> molTemplateList = mofpi_make_template_list(
        std::vector<int> { 1 }, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    // Scenario 1: base1 has no free interfaces.
    Molecule base1NoFree = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    base1NoFree.freelist.clear();
    Molecule baseTmp = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(0.1, 0.0, 0.0), Coord(0.1, 0.0, 0.0));

    bool flagCancel1 = false;
    measure_overlap_free_protein_interfaces(base1NoFree, baseTmp, flagCancel1, molTemplateList, forwardRxns, backRxns);
    EXPECT_FALSE(flagCancel1) << "base1 with no free interfaces cannot trigger a cancellation";
    std::cerr << "  flagCancel (base1 freelist empty) = " << std::boolalpha << flagCancel1 << "\n";

    // Scenario 2: baseTmp has no free interfaces.
    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    Molecule baseTmpNoFree = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(0.1, 0.0, 0.0), Coord(0.1, 0.0, 0.0));
    baseTmpNoFree.freelist.clear();

    bool flagCancel2 = false;
    measure_overlap_free_protein_interfaces(base1, baseTmpNoFree, flagCancel2, molTemplateList, forwardRxns, backRxns);
    EXPECT_FALSE(flagCancel2) << "baseTmp with no free interfaces cannot trigger a cancellation";
    std::cerr << "  flagCancel (baseTmp freelist empty) = " << std::boolalpha << flagCancel2 << "\n";
}

// -----------------------------------------------------------------------------
// Test 8: an already-true flagCancel is never cleared by the routine.
// -----------------------------------------------------------------------------
void test_mofpi_preserves_existing_flag()
{
    std::cerr << "\n[TEST] test_mofpi_preserves_existing_flag\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario:    flagCancel is passed in as true while the molecules are far\n"
              << "               apart (and, in a second call, non-interacting).\n"
              << "  Criteria:    flagCancel must still be true afterwards - the routine only\n"
              << "               ever raises the flag, never lowers it.\n";

    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    // Case 1: interacting types, but far apart.
    std::vector<MolTemplate> interacting = mofpi_make_template_list(
        std::vector<int> { 1 }, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    Molecule farTmp = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(50.0, 0.0, 0.0), Coord(50.0, 0.0, 0.0));

    bool flagCancel1 = true;
    measure_overlap_free_protein_interfaces(base1, farTmp, flagCancel1, interacting, forwardRxns, backRxns);
    EXPECT_TRUE(flagCancel1) << "A previously raised cancel flag must survive the call";
    std::cerr << "  flagCancel (interacting, far apart) = " << std::boolalpha << flagCancel1 << "\n";

    // Case 2: non-interacting types (early exit path).
    std::vector<MolTemplate> nonInteracting = mofpi_make_template_list(
        std::vector<int> {}, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    bool flagCancel2 = true;
    measure_overlap_free_protein_interfaces(base1, farTmp, flagCancel2, nonInteracting, forwardRxns, backRxns);
    EXPECT_TRUE(flagCancel2) << "The early-exit path must not clear a previously raised cancel flag";
    std::cerr << "  flagCancel (non-interacting early exit) = " << std::boolalpha << flagCancel2 << "\n";
}

// -----------------------------------------------------------------------------
// Test 9: interfaces are reaction partners by type, but the state's rxnPartners list
//         does not contain the partner interface -> the inner comparison never happens.
// -----------------------------------------------------------------------------
void test_mofpi_state_partner_mismatch()
{
    std::cerr << "\n[TEST] test_mofpi_state_partner_mismatch\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario:    molecule types are declared partners, but the state's\n"
              << "               rxnPartners list points at an unrelated absolute interface\n"
              << "               index (99) so 'absIface2 == statePartner' is never true.\n"
              << "  Criteria:    flagCancel must remain false despite the physical overlap.\n";

    // rxnPartners of the *type* allow interaction, but the interface state lists a
    // partner interface index (99) that does not exist on baseTmp.
    std::vector<MolTemplate> molTemplateList
        = mofpi_make_template_list(std::vector<int> { 1 }, std::vector<unsigned> { 99u });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    Molecule baseTmp = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(0.1, 0.0, 0.0), Coord(0.1, 0.0, 0.0));

    std::cerr << "  Interface separation (tmp) = " << mofpi_distance(base1.interfaceList[0].coord, baseTmp.tmpICoords[0])
              << " nm (overlapping), but the state partner index does not match\n";

    bool flagCancel = false;
    measure_overlap_free_protein_interfaces(base1, baseTmp, flagCancel, molTemplateList, forwardRxns, backRxns);

    EXPECT_FALSE(flagCancel) << "Interfaces that are not declared reaction partners must be ignored";
    std::cerr << "  flagCancel after call = " << std::boolalpha << flagCancel << "\n";
}

// -----------------------------------------------------------------------------
// Test 10: the routine takes its molecules by value -> the caller's copies must be
//          left untouched (only flagCancel is an out-parameter).
// -----------------------------------------------------------------------------
void test_mofpi_inputs_unmodified()
{
    std::cerr << "\n[TEST] test_mofpi_inputs_unmodified\n"
              << "  Source file: measure_overlap_free_protein_interfaces.cpp\n"
              << "  Function:    measure_overlap_free_protein_interfaces\n"
              << "  Scenario:    call the routine with overlapping interfaces and then verify\n"
              << "               that neither molecule's coordinates nor freelists changed.\n"
              << "  Criteria:    all inspected fields are bitwise-identical afterwards.\n";

    std::vector<MolTemplate> molTemplateList = mofpi_make_template_list(
        std::vector<int> { 1 }, std::vector<unsigned> { static_cast<unsigned>(kMofpiAbsIfaceB) });
    std::vector<ForwardRxn> forwardRxns = mofpi_make_forward_rxns(kMofpiBindRadius);
    std::vector<BackRxn> backRxns;

    Molecule base1 = mofpi_make_molecule(0, kMofpiAbsIfaceA, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0));
    Molecule baseTmp = mofpi_make_molecule(1, kMofpiAbsIfaceB, Coord(0.4, 0.0, 0.0), Coord(0.4, 0.0, 0.0));

    // Snapshot the values we care about.
    const Coord base1IfaceBefore = base1.interfaceList[0].coord;
    const Coord tmpIfaceBefore = baseTmp.tmpICoords[0];
    const size_t base1FreeBefore = base1.freelist.size();
    const size_t tmpFreeBefore = baseTmp.freelist.size();

    bool flagCancel = false;
    measure_overlap_free_protein_interfaces(base1, baseTmp, flagCancel, molTemplateList, forwardRxns, backRxns);
    std::cerr << "  flagCancel after call = " << std::boolalpha << flagCancel << " (overlap expected)\n";

    // Inputs are passed by value, so the caller's objects must be unchanged.
    EXPECT_DOUBLE_EQ(base1.interfaceList[0].coord.x, base1IfaceBefore.x) << "base1 interface x must be unchanged";
    EXPECT_DOUBLE_EQ(base1.interfaceList[0].coord.y, base1IfaceBefore.y) << "base1 interface y must be unchanged";
    EXPECT_DOUBLE_EQ(base1.interfaceList[0].coord.z, base1IfaceBefore.z) << "base1 interface z must be unchanged";
    EXPECT_DOUBLE_EQ(baseTmp.tmpICoords[0].x, tmpIfaceBefore.x) << "baseTmp tmp interface x must be unchanged";
    EXPECT_DOUBLE_EQ(baseTmp.tmpICoords[0].y, tmpIfaceBefore.y) << "baseTmp tmp interface y must be unchanged";
    EXPECT_DOUBLE_EQ(baseTmp.tmpICoords[0].z, tmpIfaceBefore.z) << "baseTmp tmp interface z must be unchanged";
    EXPECT_EQ(base1.freelist.size(), base1FreeBefore) << "base1 freelist size must be unchanged";
    EXPECT_EQ(baseTmp.freelist.size(), tmpFreeBefore) << "baseTmp freelist size must be unchanged";

    std::cerr << "  Input molecules verified unmodified.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each named helper runs inside its own TEST so that a
// failure in one scenario does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(MeasureOverlapFreeProteinInterfaces, NonInteractingTypes) { test_mofpi_non_interacting_types(); }
TEST(MeasureOverlapFreeProteinInterfaces, ImplicitLipidPartner) { test_mofpi_implicit_lipid_partner(); }
TEST(MeasureOverlapFreeProteinInterfaces, FarApartNoCancel) { test_mofpi_far_apart_no_cancel(); }
TEST(MeasureOverlapFreeProteinInterfaces, OverlapCancels) { test_mofpi_overlap_cancels(); }
TEST(MeasureOverlapFreeProteinInterfaces, ExactlyAtBindRadius) { test_mofpi_exactly_at_bindradius(); }
TEST(MeasureOverlapFreeProteinInterfaces, UsesTmpCoordsOfSecondMolecule)
{
    test_mofpi_uses_tmp_coords_of_second_molecule();
}
TEST(MeasureOverlapFreeProteinInterfaces, EmptyFreeLists) { test_mofpi_empty_freelists(); }
TEST(MeasureOverlapFreeProteinInterfaces, PreservesExistingFlag) { test_mofpi_preserves_existing_flag(); }
TEST(MeasureOverlapFreeProteinInterfaces, StatePartnerMismatch) { test_mofpi_state_partner_mismatch(); }
TEST(MeasureOverlapFreeProteinInterfaces, InputsUnmodified) { test_mofpi_inputs_unmodified(); }