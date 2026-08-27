/*! \file test_populate_reaction_lists.cpp
 *
 * ### Unit tests for src/parser/populate_reaction_lists.cpp
 *
 * Functions under test:
 *   - populate_reaction_lists(forwardRxns, backRxns, createDestructRxns, molTemplateList)
 *   - populate_reaction_lists_for_add(forwardRxns, backRxns, createDestructRxns,
 *                                     molTemplateList, addForwardRxnNum,
 *                                     addBackRxnNum, addCreateDestructRxnNum)
 *
 * Both functions walk a list of ForwardRxn objects and back-annotate the
 * MolTemplate interface State objects with:
 *   - State::myForwardRxns    (relative index of the reaction)
 *   - State::rxnPartners      (absolute interface index of the partner)
 *   - State::stateChangeRxns  (pair of forward/back relative indices)
 *
 * The tests below construct tiny, fully-initialised MolTemplate and ForwardRxn
 * objects by hand and then assert exactly which entries land in which list.
 *
 * NOTE ON HEADERS: parser/parser_functions.hpp transitively pulls in
 * classes/class_SimulVolume.hpp which does `#include "split.cpp"`. Including a
 * .cpp from more than one translation unit in the same test binary risks
 * duplicate-symbol link errors, so instead we include only the class headers we
 * need and declare the two functions under test ourselves. The declarations
 * below are byte-for-byte identical to the ones in parser_functions.hpp.
 */

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Rxns.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------
// Declarations of the functions under test (see note above about headers).
// -----------------------------------------------------------------------------
void populate_reaction_lists(const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList);

void populate_reaction_lists_for_add(const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList, int addForwardRxnNum,
    int addBackRxnNum, int addCreateDestructRxnNum);

// -----------------------------------------------------------------------------
// Small builders used by every test. Everything the functions under test read
// must be filled in explicitly -- the code indexes into interfaceList and
// stateList without any bounds checks.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build one Interface carrying an explicit list of states.
 *
 * \param[in] name      Interface name (only used for readable output).
 * \param[in] relIndex  Relative index of this interface inside its MolTemplate.
 * \param[in] states    List of (state identity, absolute state index) pairs.
 */
Interface prl_make_interface(const std::string& name, int relIndex,
    const std::vector<std::pair<char, int>>& states)
{
    Interface iface;
    iface.name = name;
    iface.index = relIndex;
    iface.iCoord = Coord { 0.0, 0.0, 0.0 };
    for (const auto& oneState : states) {
        // Interface::State(char iden, int index)
        iface.stateList.emplace_back(oneState.first, oneState.second);
    }
    return iface;
}

/*! \brief Build a MolTemplate holding a single interface. */
MolTemplate prl_make_template(const std::string& molName, int molTypeIndex, const Interface& iface)
{
    MolTemplate molTemplate;
    molTemplate.molName = molName;
    molTemplate.molTypeIndex = molTypeIndex;
    molTemplate.interfaceList.push_back(iface);
    return molTemplate;
}

/*! \brief Build the standard three-template system used by most of the tests.
 *
 * Index 0: "A" with interface "a"  -> one   state, absolute index 0
 * Index 1: "B" with interface "b"  -> one   state, absolute index 1
 * Index 2: "C" with interface "c"  -> two   states, absolute indices 2 ('U')
 *                                             and 3 ('P')
 */
std::vector<MolTemplate> prl_make_templates()
{
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(
        prl_make_template("A", 0, prl_make_interface("a", 0, { { '\0', 0 } })));
    molTemplateList.push_back(
        prl_make_template("B", 1, prl_make_interface("b", 0, { { '\0', 1 } })));
    molTemplateList.push_back(
        prl_make_template("C", 2, prl_make_interface("c", 0, { { 'U', 2 }, { 'P', 3 } })));
    return molTemplateList;
}

/*! \brief Convenience wrapper for building a reactant/product descriptor. */
RxnIface prl_make_rxn_iface(const std::string& name, int molTypeIndex, int absIfaceIndex, int relIfaceIndex)
{
    return RxnIface(name, molTypeIndex, absIfaceIndex, relIfaceIndex, '\0', false);
}

/*! \brief Build a bimolecular ForwardRxn between two given interfaces.
 *
 * A non-empty rateList is mandatory: populate_reaction_lists_for_add prints
 * forwardRxns[...].rateList[0].rate for every reaction already registered on
 * the first state of the first reactant interface.
 */
ForwardRxn prl_make_bimolecular_rxn(int relRxnIndex, const RxnIface& react1, const RxnIface& react2,
    bool isSymmetric, double rate)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.relRxnIndex = relRxnIndex;
    rxn.absRxnIndex = relRxnIndex;
    rxn.isSymmetric = isSymmetric;
    rxn.isReversible = false;
    rxn.conjBackRxnIndex = -1;
    rxn.reactantListNew.push_back(react1);
    rxn.reactantListNew.push_back(react2);
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = rate;
    return rxn;
}

/*! \brief Build a unimolecular state change reaction (react -> product). */
ForwardRxn prl_make_unimol_state_change_rxn(int relRxnIndex, int conjBackRxnIndex, bool isReversible,
    const RxnIface& reactant, const RxnIface& product)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::uniMolStateChange;
    rxn.relRxnIndex = relRxnIndex;
    rxn.absRxnIndex = relRxnIndex;
    rxn.isReversible = isReversible;
    rxn.conjBackRxnIndex = conjBackRxnIndex;
    rxn.hasStateChange = true;
    rxn.reactantListNew.push_back(reactant);
    rxn.productListNew.push_back(product);
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = 1.0;
    return rxn;
}

/*! \brief Build a bimolecular state change reaction (facilitator + target). */
ForwardRxn prl_make_bimol_state_change_rxn(int relRxnIndex, int conjBackRxnIndex, bool isReversible,
    const RxnIface& facilitatorReact, const RxnIface& targetReact,
    const RxnIface& facilitatorProd, const RxnIface& targetProd)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::biMolStateChange;
    rxn.relRxnIndex = relRxnIndex;
    rxn.absRxnIndex = relRxnIndex;
    rxn.isReversible = isReversible;
    rxn.conjBackRxnIndex = conjBackRxnIndex;
    rxn.hasStateChange = true;
    rxn.reactantListNew.push_back(facilitatorReact);
    rxn.reactantListNew.push_back(targetReact);
    rxn.productListNew.push_back(facilitatorProd);
    rxn.productListNew.push_back(targetProd);
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = 1.0;
    return rxn;
}

/*! \brief Dump a state's bookkeeping lists to stderr for readability. */
void prl_report_state(const std::string& label, const Interface::State& state)
{
    std::cerr << "    " << label << " (abs index " << state.index << "): myForwardRxns={";
    for (const auto& rxn : state.myForwardRxns)
        std::cerr << rxn << ' ';
    std::cerr << "} rxnPartners={";
    for (const auto& partner : state.rxnPartners)
        std::cerr << partner << ' ';
    std::cerr << "} stateChangeRxns={";
    for (const auto& pair : state.stateChangeRxns)
        std::cerr << '(' << pair.first << ',' << pair.second << ") ";
    std::cerr << "}\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: bimolecular reaction between two single-state interfaces.
// -----------------------------------------------------------------------------
void test_prl_bimolecular_single_states()
{
    std::cerr << "\n[TEST] test_prl_bimolecular_single_states\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists (bimolecular branch)\n"
              << "  Scenario:    A(a) + B(b), both interfaces have a single state.\n"
              << "  Pass criteria: each state records the reaction index once and\n"
              << "                 records the partner's absolute interface index.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    // A(a) has abs index 0, B(b) has abs index 1.
    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_bimolecular_rxn(0, prl_make_rxn_iface("a", 0, 0, 0),
        prl_make_rxn_iface("b", 1, 1, 0), /*isSymmetric=*/false, /*rate=*/10.0));

    const std::vector<BackRxn> backRxns {};
    const std::vector<CreateDestructRxn> createDestructRxns {};

    std::cerr << "  Calling populate_reaction_lists...\n";
    populate_reaction_lists(forwardRxns, backRxns, createDestructRxns, molTemplateList);

    const Interface::State& stateA = molTemplateList[0].interfaceList[0].stateList[0];
    const Interface::State& stateB = molTemplateList[1].interfaceList[0].stateList[0];
    prl_report_state("A(a)", stateA);
    prl_report_state("B(b)", stateB);

    // Both interfaces should have exactly one forward reaction registered.
    ASSERT_EQ(stateA.myForwardRxns.size(), 1u) << "A(a) should list exactly one forward reaction";
    EXPECT_EQ(stateA.myForwardRxns[0], 0u) << "A(a) should list reaction relRxnIndex 0";
    ASSERT_EQ(stateB.myForwardRxns.size(), 1u) << "B(b) should list exactly one forward reaction";
    EXPECT_EQ(stateB.myForwardRxns[0], 0u) << "B(b) should list reaction relRxnIndex 0";

    // The partner lists are cross-linked by absolute interface index.
    ASSERT_EQ(stateA.rxnPartners.size(), 1u) << "A(a) should have one reaction partner";
    EXPECT_EQ(stateA.rxnPartners[0], 1u) << "A(a)'s partner is B(b) with absolute index 1";
    ASSERT_EQ(stateB.rxnPartners.size(), 1u) << "B(b) should have one reaction partner";
    EXPECT_EQ(stateB.rxnPartners[0], 0u) << "B(b)'s partner is A(a) with absolute index 0";

    // A plain bimolecular reaction is not a state change.
    EXPECT_TRUE(stateA.stateChangeRxns.empty()) << "Bimolecular binding is not a state change for A";
    EXPECT_TRUE(stateB.stateChangeRxns.empty()) << "Bimolecular binding is not a state change for B";

    // The untouched template must remain pristine.
    EXPECT_TRUE(molTemplateList[2].interfaceList[0].stateList[0].myForwardRxns.empty())
        << "C(c~U) takes no part in this reaction";
    EXPECT_TRUE(molTemplateList[2].interfaceList[0].stateList[1].myForwardRxns.empty())
        << "C(c~P) takes no part in this reaction";
}

// -----------------------------------------------------------------------------
// Test 2: symmetric (homodimerisation) bimolecular reaction.
// -----------------------------------------------------------------------------
void test_prl_bimolecular_symmetric()
{
    std::cerr << "\n[TEST] test_prl_bimolecular_symmetric\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists (symmetric bimolecular)\n"
              << "  Scenario:    A(a) + A(a) with isSymmetric == true.\n"
              << "  Pass criteria: only the first reactant interface is touched, and\n"
              << "                 its partner is its own absolute interface index.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_bimolecular_rxn(4, prl_make_rxn_iface("a", 0, 0, 0),
        prl_make_rxn_iface("a", 0, 0, 0), /*isSymmetric=*/true, /*rate=*/2.5));

    populate_reaction_lists(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    const Interface::State& stateA = molTemplateList[0].interfaceList[0].stateList[0];
    prl_report_state("A(a)", stateA);

    // The symmetric branch pushes the reaction exactly once (not twice).
    ASSERT_EQ(stateA.myForwardRxns.size(), 1u)
        << "A symmetric reaction should only be registered once on the interface";
    EXPECT_EQ(stateA.myForwardRxns[0], 4u) << "The stored index is the reaction's relRxnIndex (4)";

    // Its partner is itself.
    ASSERT_EQ(stateA.rxnPartners.size(), 1u) << "A symmetric reaction adds a single partner entry";
    EXPECT_EQ(stateA.rxnPartners[0], 0u) << "The partner of a symmetric reaction is its own absolute index (0)";

    // The second molecule template must not be touched at all.
    EXPECT_TRUE(molTemplateList[1].interfaceList[0].stateList[0].myForwardRxns.empty())
        << "B is not involved in the symmetric A+A reaction";
}

// -----------------------------------------------------------------------------
// Test 3: bimolecular reaction where the FIRST reactant has multiple states.
// -----------------------------------------------------------------------------
void test_prl_bimolecular_multistate_first_reactant()
{
    std::cerr << "\n[TEST] test_prl_bimolecular_multistate_first_reactant\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists (multi-state reactant 1)\n"
              << "  Scenario:    C(c~P) + A(a), where C(c) has states abs 2 and abs 3.\n"
              << "  Pass criteria: only the matching state (abs 3) is annotated; the\n"
              << "                 other state (abs 2) is left untouched.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    // Reactant 1 = C(c) in its SECOND state (absolute index 3).
    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_bimolecular_rxn(1, prl_make_rxn_iface("c", 2, 3, 0),
        prl_make_rxn_iface("a", 0, 0, 0), /*isSymmetric=*/false, /*rate=*/7.0));

    populate_reaction_lists(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    const Interface::State& stateCU = molTemplateList[2].interfaceList[0].stateList[0]; // abs 2
    const Interface::State& stateCP = molTemplateList[2].interfaceList[0].stateList[1]; // abs 3
    const Interface::State& stateA = molTemplateList[0].interfaceList[0].stateList[0];  // abs 0
    prl_report_state("C(c~U)", stateCU);
    prl_report_state("C(c~P)", stateCP);
    prl_report_state("A(a)", stateA);

    // Only the state whose absolute index matches gets the reaction.
    EXPECT_TRUE(stateCU.myForwardRxns.empty()) << "State abs index 2 does not match the reactant (abs 3)";
    EXPECT_TRUE(stateCU.rxnPartners.empty()) << "State abs index 2 gains no partners";
    ASSERT_EQ(stateCP.myForwardRxns.size(), 1u) << "State abs index 3 should record the reaction";
    EXPECT_EQ(stateCP.myForwardRxns[0], 1u) << "The recorded index is relRxnIndex 1";

    // The single-state partner is linked both ways.
    ASSERT_EQ(stateA.myForwardRxns.size(), 1u) << "A(a) should record the reaction once";
    ASSERT_EQ(stateA.rxnPartners.size(), 1u) << "A(a) should gain a single partner";
    EXPECT_EQ(stateA.rxnPartners[0], 3u) << "A(a)'s partner is the C(c~P) state, absolute index 3";
    ASSERT_EQ(stateCP.rxnPartners.size(), 1u) << "C(c~P) should gain a single partner";
    EXPECT_EQ(stateCP.rxnPartners[0], 0u) << "C(c~P)'s partner is A(a), absolute index 0";
}

// -----------------------------------------------------------------------------
// Test 4: bimolecular reaction where the SECOND reactant has multiple states.
// -----------------------------------------------------------------------------
void test_prl_bimolecular_multistate_second_reactant()
{
    std::cerr << "\n[TEST] test_prl_bimolecular_multistate_second_reactant\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists (multi-state reactant 2)\n"
              << "  Scenario:    A(a) + C(c~U), matching the FIRST state of C(c).\n"
              << "  Pass criteria: the abs-2 state is annotated and cross-linked, the\n"
              << "                 abs-3 state stays empty.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    // Reactant 2 = C(c) in its FIRST state (absolute index 2).
    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_bimolecular_rxn(3, prl_make_rxn_iface("a", 0, 0, 0),
        prl_make_rxn_iface("c", 2, 2, 0), /*isSymmetric=*/false, /*rate=*/5.0));

    populate_reaction_lists(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    const Interface::State& stateA = molTemplateList[0].interfaceList[0].stateList[0];
    const Interface::State& stateCU = molTemplateList[2].interfaceList[0].stateList[0]; // abs 2
    const Interface::State& stateCP = molTemplateList[2].interfaceList[0].stateList[1]; // abs 3
    prl_report_state("A(a)", stateA);
    prl_report_state("C(c~U)", stateCU);
    prl_report_state("C(c~P)", stateCP);

    ASSERT_EQ(stateA.myForwardRxns.size(), 1u) << "A(a) records the reaction once";
    EXPECT_EQ(stateA.myForwardRxns[0], 3u) << "A(a) records relRxnIndex 3";
    ASSERT_EQ(stateCU.myForwardRxns.size(), 1u) << "C(c~U) records the reaction once";
    EXPECT_EQ(stateCU.myForwardRxns[0], 3u) << "C(c~U) records relRxnIndex 3";
    EXPECT_TRUE(stateCP.myForwardRxns.empty()) << "C(c~P) is not a reactant here";

    ASSERT_EQ(stateA.rxnPartners.size(), 1u) << "A(a) gains one partner";
    EXPECT_EQ(stateA.rxnPartners[0], 2u) << "A(a)'s partner is C(c~U), absolute index 2";
    ASSERT_EQ(stateCU.rxnPartners.size(), 1u) << "C(c~U) gains one partner";
    EXPECT_EQ(stateCU.rxnPartners[0], 0u) << "C(c~U)'s partner is A(a), absolute index 0";
    EXPECT_TRUE(stateCP.rxnPartners.empty()) << "C(c~P) gains no partners";
}

// -----------------------------------------------------------------------------
// Test 5: reversible unimolecular state change.
// -----------------------------------------------------------------------------
void test_prl_unimolecular_state_change()
{
    std::cerr << "\n[TEST] test_prl_unimolecular_state_change\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists (uniMolStateChange branch)\n"
              << "  Scenario:    C(c~U) <-> C(c~P), relRxnIndex 5, conjBackRxnIndex 2.\n"
              << "  Pass criteria: the reactant state stores (5,2), the product state\n"
              << "                 stores the reversed pair (2,5), and myForwardRxns\n"
              << "                 stays empty for both.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_unimol_state_change_rxn(/*relRxnIndex=*/5, /*conjBackRxnIndex=*/2,
        /*isReversible=*/true, prl_make_rxn_iface("c", 2, 2, 0), prl_make_rxn_iface("c", 2, 3, 0)));

    populate_reaction_lists(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    const Interface::State& stateCU = molTemplateList[2].interfaceList[0].stateList[0]; // abs 2 (reactant)
    const Interface::State& stateCP = molTemplateList[2].interfaceList[0].stateList[1]; // abs 3 (product)
    prl_report_state("C(c~U)", stateCU);
    prl_report_state("C(c~P)", stateCP);

    // Reactant side stores (forward, back).
    ASSERT_EQ(stateCU.stateChangeRxns.size(), 1u) << "The reactant state stores one state-change pair";
    EXPECT_EQ(stateCU.stateChangeRxns[0].first, 5) << "Forward index for the reactant state is 5";
    EXPECT_EQ(stateCU.stateChangeRxns[0].second, 2) << "Back index for the reactant state is 2";

    // Product side stores the reversed pair (back, forward) since it is reversible.
    ASSERT_EQ(stateCP.stateChangeRxns.size(), 1u) << "The product state stores one state-change pair";
    EXPECT_EQ(stateCP.stateChangeRxns[0].first, 2) << "Forward index for the product state is the back index 2";
    EXPECT_EQ(stateCP.stateChangeRxns[0].second, 5) << "Back index for the product state is the forward index 5";

    // Unimolecular state changes are deliberately NOT put into myForwardRxns.
    EXPECT_TRUE(stateCU.myForwardRxns.empty()) << "uniMolStateChange must not touch myForwardRxns";
    EXPECT_TRUE(stateCP.myForwardRxns.empty()) << "uniMolStateChange must not touch myForwardRxns";
    EXPECT_TRUE(stateCU.rxnPartners.empty()) << "uniMolStateChange must not add reaction partners";
    EXPECT_TRUE(stateCP.rxnPartners.empty()) << "uniMolStateChange must not add reaction partners";
}

// -----------------------------------------------------------------------------
// Test 6: irreversible unimolecular state change only annotates the reactant.
// -----------------------------------------------------------------------------
void test_prl_unimolecular_state_change_irreversible()
{
    std::cerr << "\n[TEST] test_prl_unimolecular_state_change_irreversible\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists (uniMolStateChange branch)\n"
              << "  Scenario:    C(c~U) -> C(c~P) with isReversible == false.\n"
              << "  Pass criteria: only the reactant state gains a stateChangeRxns\n"
              << "                 entry; the product state remains empty.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_unimol_state_change_rxn(/*relRxnIndex=*/9, /*conjBackRxnIndex=*/-1,
        /*isReversible=*/false, prl_make_rxn_iface("c", 2, 2, 0), prl_make_rxn_iface("c", 2, 3, 0)));

    populate_reaction_lists(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    const Interface::State& stateCU = molTemplateList[2].interfaceList[0].stateList[0];
    const Interface::State& stateCP = molTemplateList[2].interfaceList[0].stateList[1];
    prl_report_state("C(c~U)", stateCU);
    prl_report_state("C(c~P)", stateCP);

    ASSERT_EQ(stateCU.stateChangeRxns.size(), 1u) << "The reactant state records the irreversible reaction";
    EXPECT_EQ(stateCU.stateChangeRxns[0].first, 9) << "Forward index is 9";
    EXPECT_EQ(stateCU.stateChangeRxns[0].second, -1) << "Back index is -1 for an irreversible reaction";
    EXPECT_TRUE(stateCP.stateChangeRxns.empty())
        << "The product state must stay empty when the reaction is irreversible";
}

// -----------------------------------------------------------------------------
// Test 7: reversible bimolecular state change (facilitator + target).
// -----------------------------------------------------------------------------
void test_prl_bimolecular_state_change()
{
    std::cerr << "\n[TEST] test_prl_bimolecular_state_change\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists (biMolStateChange branch)\n"
              << "  Scenario:    A(a) + C(c~U) <-> A(a) + C(c~P), relRxnIndex 7,\n"
              << "               conjBackRxnIndex 3.\n"
              << "  Pass criteria: both C states and the facilitator A state gain\n"
              << "                 stateChangeRxns, myForwardRxns and rxnPartners.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_bimol_state_change_rxn(/*relRxnIndex=*/7, /*conjBackRxnIndex=*/3,
        /*isReversible=*/true,
        /*facilitatorReact=*/prl_make_rxn_iface("a", 0, 0, 0),
        /*targetReact=*/prl_make_rxn_iface("c", 2, 2, 0),
        /*facilitatorProd=*/prl_make_rxn_iface("a", 0, 0, 0),
        /*targetProd=*/prl_make_rxn_iface("c", 2, 3, 0)));

    populate_reaction_lists(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    const Interface::State& stateA = molTemplateList[0].interfaceList[0].stateList[0];  // facilitator, abs 0
    const Interface::State& stateCU = molTemplateList[2].interfaceList[0].stateList[0]; // abs 2 (reactant)
    const Interface::State& stateCP = molTemplateList[2].interfaceList[0].stateList[1]; // abs 3 (product)
    prl_report_state("A(a) facilitator", stateA);
    prl_report_state("C(c~U)", stateCU);
    prl_report_state("C(c~P)", stateCP);

    // Target reactant state: forward pair, forward reaction, facilitator partner.
    ASSERT_EQ(stateCU.stateChangeRxns.size(), 1u) << "C(c~U) records one state-change pair";
    EXPECT_EQ(stateCU.stateChangeRxns[0].first, 7) << "C(c~U) forward index is 7";
    EXPECT_EQ(stateCU.stateChangeRxns[0].second, 3) << "C(c~U) back index is 3";
    ASSERT_EQ(stateCU.myForwardRxns.size(), 1u) << "C(c~U) records the forward reaction";
    EXPECT_EQ(stateCU.myForwardRxns[0], 7u) << "C(c~U) records relRxnIndex 7";
    ASSERT_EQ(stateCU.rxnPartners.size(), 1u) << "C(c~U) records the facilitator as its partner";
    EXPECT_EQ(stateCU.rxnPartners[0], 0u) << "The facilitator has absolute index 0";

    // Target product state: reversed pair because the reaction is reversible.
    ASSERT_EQ(stateCP.stateChangeRxns.size(), 1u) << "C(c~P) records one state-change pair";
    EXPECT_EQ(stateCP.stateChangeRxns[0].first, 3) << "C(c~P) forward index is the back index 3";
    EXPECT_EQ(stateCP.stateChangeRxns[0].second, 7) << "C(c~P) back index is the forward index 7";
    ASSERT_EQ(stateCP.myForwardRxns.size(), 1u) << "C(c~P) also records the forward reaction index";
    EXPECT_EQ(stateCP.myForwardRxns[0], 7u) << "C(c~P) records relRxnIndex 7";
    ASSERT_EQ(stateCP.rxnPartners.size(), 1u) << "C(c~P) records the product-side facilitator partner";
    EXPECT_EQ(stateCP.rxnPartners[0], 0u) << "The product-side facilitator has absolute index 0";

    // Facilitator: the reactant branch matches first, so exactly one entry each.
    ASSERT_EQ(stateA.stateChangeRxns.size(), 1u) << "The facilitator gets exactly one pair (first branch matches)";
    EXPECT_EQ(stateA.stateChangeRxns[0].first, 7) << "Facilitator forward index is 7";
    EXPECT_EQ(stateA.stateChangeRxns[0].second, 3) << "Facilitator back index is 3";
    ASSERT_EQ(stateA.myForwardRxns.size(), 1u) << "The facilitator records the forward reaction once";
    ASSERT_EQ(stateA.rxnPartners.size(), 1u) << "The facilitator records one partner";
    EXPECT_EQ(stateA.rxnPartners[0], 2u) << "The facilitator's partner is the target reactant, absolute index 2";
}

// -----------------------------------------------------------------------------
// Test 8: unhandled reaction types leave every template untouched, and the
//         backRxns / createDestructRxns arguments are ignored entirely.
// -----------------------------------------------------------------------------
void test_prl_ignored_reaction_types()
{
    std::cerr << "\n[TEST] test_prl_ignored_reaction_types\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists (default / no-op path)\n"
              << "  Scenario:    a zerothOrderCreation and a destruction reaction, plus\n"
              << "               non-empty backRxns and createDestructRxns lists.\n"
              << "  Pass criteria: no MolTemplate state list is modified at all.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    // Neither of these reaction types is handled by the function.
    std::vector<ForwardRxn> forwardRxns;
    {
        ForwardRxn creationRxn;
        creationRxn.rxnType = ReactionType::zerothOrderCreation;
        creationRxn.relRxnIndex = 0;
        creationRxn.reactantListNew.push_back(prl_make_rxn_iface("a", 0, 0, 0));
        creationRxn.productListNew.push_back(prl_make_rxn_iface("a", 0, 0, 0));
        creationRxn.rateList.emplace_back();
        forwardRxns.push_back(creationRxn);

        ForwardRxn destructionRxn;
        destructionRxn.rxnType = ReactionType::destruction;
        destructionRxn.relRxnIndex = 1;
        destructionRxn.reactantListNew.push_back(prl_make_rxn_iface("b", 1, 1, 0));
        destructionRxn.rateList.emplace_back();
        forwardRxns.push_back(destructionRxn);
    }

    // These two lists are accepted but never read by the implementation.
    std::vector<BackRxn> backRxns(1);
    std::vector<CreateDestructRxn> createDestructRxns(1);

    populate_reaction_lists(forwardRxns, backRxns, createDestructRxns, molTemplateList);

    // Walk every state of every template and confirm nothing was recorded.
    for (const auto& molTemplate : molTemplateList) {
        for (const auto& iface : molTemplate.interfaceList) {
            for (const auto& state : iface.stateList) {
                prl_report_state(molTemplate.molName + "(" + iface.name + ")", state);
                EXPECT_TRUE(state.myForwardRxns.empty())
                    << "Unhandled reaction types must not add forward reactions";
                EXPECT_TRUE(state.rxnPartners.empty())
                    << "Unhandled reaction types must not add reaction partners";
                EXPECT_TRUE(state.stateChangeRxns.empty())
                    << "Unhandled reaction types must not add state-change reactions";
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 9: populate_reaction_lists_for_add only processes reactions at or after
//         the supplied starting index.
// -----------------------------------------------------------------------------
void test_prl_for_add_skips_existing_reactions()
{
    std::cerr << "\n[TEST] test_prl_for_add_skips_existing_reactions\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists_for_add\n"
              << "  Scenario:    reaction 0 (A+B) was already registered by a previous\n"
              << "               call; reaction 1 (A+C~U) is newly added, so the\n"
              << "               function is called with addForwardRxnNum == 1.\n"
              << "  Pass criteria: reaction 0's bookkeeping is not duplicated and\n"
              << "                 reaction 1 is appended to the existing lists.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    // The full reaction list as it looks AFTER the add-file has been parsed.
    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_bimolecular_rxn(0, prl_make_rxn_iface("a", 0, 0, 0),
        prl_make_rxn_iface("b", 1, 1, 0), /*isSymmetric=*/false, /*rate=*/11.0));
    forwardRxns.push_back(prl_make_bimolecular_rxn(1, prl_make_rxn_iface("a", 0, 0, 0),
        prl_make_rxn_iface("c", 2, 2, 0), /*isSymmetric=*/false, /*rate=*/22.0));

    // Step 1: emulate the original run, which only knew about reaction 0.
    std::vector<ForwardRxn> originalRxns { forwardRxns[0] };
    std::cerr << "  Pre-populating with reaction 0 only...\n";
    populate_reaction_lists(originalRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    // Step 2: the restart/add path processes only reactions from index 1 on.
    std::cerr << "  Calling populate_reaction_lists_for_add with addForwardRxnNum = 1...\n";
    populate_reaction_lists_for_add(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {},
        molTemplateList, /*addForwardRxnNum=*/1, /*addBackRxnNum=*/0, /*addCreateDestructRxnNum=*/0);

    const Interface::State& stateA = molTemplateList[0].interfaceList[0].stateList[0];
    const Interface::State& stateB = molTemplateList[1].interfaceList[0].stateList[0];
    const Interface::State& stateCU = molTemplateList[2].interfaceList[0].stateList[0];
    const Interface::State& stateCP = molTemplateList[2].interfaceList[0].stateList[1];
    prl_report_state("A(a)", stateA);
    prl_report_state("B(b)", stateB);
    prl_report_state("C(c~U)", stateCU);
    prl_report_state("C(c~P)", stateCP);

    // A participates in both reactions -> two entries, in order.
    ASSERT_EQ(stateA.myForwardRxns.size(), 2u) << "A(a) participates in reactions 0 and 1";
    EXPECT_EQ(stateA.myForwardRxns[0], 0u) << "The pre-existing entry (reaction 0) is preserved";
    EXPECT_EQ(stateA.myForwardRxns[1], 1u) << "The newly added reaction 1 is appended";
    ASSERT_EQ(stateA.rxnPartners.size(), 2u) << "A(a) now has two reaction partners";
    EXPECT_EQ(stateA.rxnPartners[0], 1u) << "First partner is B(b), absolute index 1";
    EXPECT_EQ(stateA.rxnPartners[1], 2u) << "Second partner is C(c~U), absolute index 2";

    // B is only in reaction 0 and must not have been re-processed.
    ASSERT_EQ(stateB.myForwardRxns.size(), 1u) << "B(b) must not be re-registered by the add pass";
    EXPECT_EQ(stateB.myForwardRxns[0], 0u) << "B(b) still only knows about reaction 0";
    ASSERT_EQ(stateB.rxnPartners.size(), 1u) << "B(b) must not gain a duplicate partner";
    EXPECT_EQ(stateB.rxnPartners[0], 0u) << "B(b)'s only partner is A(a), absolute index 0";

    // The newly added C(c~U) state picks up reaction 1; C(c~P) stays empty.
    ASSERT_EQ(stateCU.myForwardRxns.size(), 1u) << "C(c~U) is registered by the add pass";
    EXPECT_EQ(stateCU.myForwardRxns[0], 1u) << "C(c~U) records relRxnIndex 1";
    ASSERT_EQ(stateCU.rxnPartners.size(), 1u) << "C(c~U) gains one partner";
    EXPECT_EQ(stateCU.rxnPartners[0], 0u) << "C(c~U)'s partner is A(a), absolute index 0";
    EXPECT_TRUE(stateCP.myForwardRxns.empty()) << "C(c~P) is not a reactant of the added reaction";
    EXPECT_TRUE(stateCP.rxnPartners.empty()) << "C(c~P) gains no partners";
}

// -----------------------------------------------------------------------------
// Test 10: populate_reaction_lists_for_add handles state-change reactions in the
//          added range exactly like populate_reaction_lists does.
// -----------------------------------------------------------------------------
void test_prl_for_add_state_change_reaction()
{
    std::cerr << "\n[TEST] test_prl_for_add_state_change_reaction\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists_for_add (state change path)\n"
              << "  Scenario:    reaction 0 (A+B, already handled) plus a newly added\n"
              << "               reversible unimolecular state change at index 1.\n"
              << "  Pass criteria: the state change pairs are recorded on C(c) and the\n"
              << "                 already-processed reaction 0 is not duplicated.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_bimolecular_rxn(0, prl_make_rxn_iface("a", 0, 0, 0),
        prl_make_rxn_iface("b", 1, 1, 0), /*isSymmetric=*/false, /*rate=*/1.0));
    forwardRxns.push_back(prl_make_unimol_state_change_rxn(/*relRxnIndex=*/1, /*conjBackRxnIndex=*/0,
        /*isReversible=*/true, prl_make_rxn_iface("c", 2, 2, 0), prl_make_rxn_iface("c", 2, 3, 0)));

    // Pre-populate with reaction 0 only, mimicking the original simulation run.
    std::vector<ForwardRxn> originalRxns { forwardRxns[0] };
    populate_reaction_lists(originalRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    std::cerr << "  Calling populate_reaction_lists_for_add with addForwardRxnNum = 1...\n";
    populate_reaction_lists_for_add(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {},
        molTemplateList, /*addForwardRxnNum=*/1, /*addBackRxnNum=*/0, /*addCreateDestructRxnNum=*/0);

    const Interface::State& stateA = molTemplateList[0].interfaceList[0].stateList[0];
    const Interface::State& stateCU = molTemplateList[2].interfaceList[0].stateList[0];
    const Interface::State& stateCP = molTemplateList[2].interfaceList[0].stateList[1];
    prl_report_state("A(a)", stateA);
    prl_report_state("C(c~U)", stateCU);
    prl_report_state("C(c~P)", stateCP);

    // Reaction 0 must not be replayed.
    ASSERT_EQ(stateA.myForwardRxns.size(), 1u) << "A(a) must still only list the pre-existing reaction 0";
    EXPECT_EQ(stateA.myForwardRxns[0], 0u) << "A(a)'s entry is relRxnIndex 0";

    // The added state change is recorded on both states of C(c).
    ASSERT_EQ(stateCU.stateChangeRxns.size(), 1u) << "C(c~U) gains the added state-change pair";
    EXPECT_EQ(stateCU.stateChangeRxns[0].first, 1) << "C(c~U) forward index is 1";
    EXPECT_EQ(stateCU.stateChangeRxns[0].second, 0) << "C(c~U) back index is 0";
    ASSERT_EQ(stateCP.stateChangeRxns.size(), 1u) << "C(c~P) gains the reversed state-change pair";
    EXPECT_EQ(stateCP.stateChangeRxns[0].first, 0) << "C(c~P) forward index is the back index 0";
    EXPECT_EQ(stateCP.stateChangeRxns[0].second, 1) << "C(c~P) back index is the forward index 1";

    // Unimolecular state changes never populate myForwardRxns.
    EXPECT_TRUE(stateCU.myForwardRxns.empty()) << "uniMolStateChange must not populate myForwardRxns";
    EXPECT_TRUE(stateCP.myForwardRxns.empty()) << "uniMolStateChange must not populate myForwardRxns";
}

// -----------------------------------------------------------------------------
// Test 11: populate_reaction_lists_for_add with a start index equal to the list
//          size is a complete no-op.
// -----------------------------------------------------------------------------
void test_prl_for_add_noop_when_nothing_added()
{
    std::cerr << "\n[TEST] test_prl_for_add_noop_when_nothing_added\n"
              << "  Source file: src/parser/populate_reaction_lists.cpp\n"
              << "  Function:    populate_reaction_lists_for_add\n"
              << "  Scenario:    addForwardRxnNum equals forwardRxns.size(), i.e. the\n"
              << "               add-file introduced no new reactions.\n"
              << "  Pass criteria: the MolTemplate state lists are byte-identical to\n"
              << "                 their state before the call.\n";

    std::vector<MolTemplate> molTemplateList = prl_make_templates();

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(prl_make_bimolecular_rxn(0, prl_make_rxn_iface("a", 0, 0, 0),
        prl_make_rxn_iface("b", 1, 1, 0), /*isSymmetric=*/false, /*rate=*/3.0));

    // Populate normally first.
    populate_reaction_lists(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {}, molTemplateList);

    // Snapshot the lists we care about.
    const std::vector<unsigned> beforeA = molTemplateList[0].interfaceList[0].stateList[0].myForwardRxns;
    const std::vector<unsigned> beforeAPartners = molTemplateList[0].interfaceList[0].stateList[0].rxnPartners;
    const std::vector<unsigned> beforeB = molTemplateList[1].interfaceList[0].stateList[0].myForwardRxns;

    std::cerr << "  Calling populate_reaction_lists_for_add with addForwardRxnNum = "
              << forwardRxns.size() << " (nothing new)...\n";
    populate_reaction_lists_for_add(forwardRxns, std::vector<BackRxn> {}, std::vector<CreateDestructRxn> {},
        molTemplateList, static_cast<int>(forwardRxns.size()), 0, 0);

    const Interface::State& stateA = molTemplateList[0].interfaceList[0].stateList[0];
    const Interface::State& stateB = molTemplateList[1].interfaceList[0].stateList[0];
    prl_report_state("A(a)", stateA);
    prl_report_state("B(b)", stateB);

    EXPECT_EQ(stateA.myForwardRxns, beforeA) << "A(a)'s forward reaction list must be unchanged";
    EXPECT_EQ(stateA.rxnPartners, beforeAPartners) << "A(a)'s partner list must be unchanged";
    EXPECT_EQ(stateB.myForwardRxns, beforeB) << "B(b)'s forward reaction list must be unchanged";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper runs inside its own TEST so that a
// failure in one scenario does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(PopulateReactionLists, BimolecularSingleStates) { test_prl_bimolecular_single_states(); }
TEST(PopulateReactionLists, BimolecularSymmetric) { test_prl_bimolecular_symmetric(); }
TEST(PopulateReactionLists, BimolecularMultistateFirstReactant) { test_prl_bimolecular_multistate_first_reactant(); }
TEST(PopulateReactionLists, BimolecularMultistateSecondReactant) { test_prl_bimolecular_multistate_second_reactant(); }
TEST(PopulateReactionLists, UnimolecularStateChange) { test_prl_unimolecular_state_change(); }
TEST(PopulateReactionLists, UnimolecularStateChangeIrreversible) { test_prl_unimolecular_state_change_irreversible(); }
TEST(PopulateReactionLists, BimolecularStateChange) { test_prl_bimolecular_state_change(); }
TEST(PopulateReactionLists, IgnoredReactionTypes) { test_prl_ignored_reaction_types(); }
TEST(PopulateReactionListsForAdd, SkipsExistingReactions) { test_prl_for_add_skips_existing_reactions(); }
TEST(PopulateReactionListsForAdd, StateChangeReaction) { test_prl_for_add_state_change_reaction(); }
TEST(PopulateReactionListsForAdd, NoopWhenNothingAdded) { test_prl_for_add_noop_when_nothing_added(); }