/*! \file test_check_for_state_change.cpp
 *
 * ### Unit test for ../src/parser/check_for_state_change.cpp
 *
 * Function under test:
 *
 *     void check_for_state_change(ParsedMol::IfaceInfo& targetIface,
 *                                 ParsedMol&            targetMol,
 *                                 ParsedRxn&            parsedRxn)
 *
 * Behaviour of the function (read directly from the implementation):
 *
 *   1. It walks over every ParsedMol in `parsedRxn.productList`.
 *   2. For every product whose `molName` equals `targetMol.molName` it copies
 *      `targetMol.molTypeIndex` onto the product (pure bookkeeping) -- this
 *      happens whether or not a state change is eventually detected.
 *   3. For each interface of such a product it calls `areSameExceptState()`
 *      with the reactant interface.  When that predicate is true the reactant
 *      interface and the product interface are both flagged with
 *      `Involvement::stateChange` and `parsedRxn.hasStateChange` is set true.
 *   4. Products whose `molName` differs are skipped entirely -- their
 *      `molTypeIndex` is left untouched.
 *
 * `areSameExceptState()` (declared in parser/parser_functions.hpp) compares two
 * ParsedMol::IfaceInfo objects and reports true only when they agree in every
 * identifying field but carry *different* state characters.  All test fixtures
 * below therefore build interface pairs that are byte-for-byte identical except
 * for the `state` member, so the predicate is exercised unambiguously.
 *
 * No path in this function calls exit()/abort(), so all cases below are safe.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (file-scope, uniquely named to avoid collisions in the suite)
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully-initialised ParsedMol::IfaceInfo.
 *
 * Every field that `areSameExceptState()` can possibly inspect is set
 * explicitly, so two interfaces built with the same arguments except `state`
 * differ *only* in their state character.
 *
 * \param[in] name         interface name (e.g. "a")
 * \param[in] state        the state identifier character (e.g. 'U')
 * \param[in] isBound      is the interface bound in this reaction side
 * \param[in] bondIndex    bond label (-1 when unbound)
 * \param[in] speciesIndex which reaction species the interface sits on
 * \param[in] molTypeIndex parent MolTemplate index
 * \param[in] relIndex     relative interface index on the molecule
 * \param[in] absIndex     absolute interface/state index
 */
ParsedMol::IfaceInfo cfsc_make_iface(const std::string& name, char state, bool isBound,
    int bondIndex, int speciesIndex, int molTypeIndex, int relIndex, int absIndex)
{
    // Use the fullest available constructor, then fill in the remaining members.
    ParsedMol::IfaceInfo info(name, state, isBound, bondIndex, Involvement::possible, speciesIndex);
    info.molTypeIndex = molTypeIndex;
    info.relIndex = relIndex;
    info.absIndex = absIndex;
    return info;
}

/*! \brief Build a ParsedMol carrying a single interface. */
ParsedMol cfsc_make_mol(const std::string& molName, int molTypeIndex, int specieIndex,
    const ParsedMol::IfaceInfo& iface)
{
    ParsedMol mol;
    mol.molName = molName;
    mol.molTypeIndex = molTypeIndex;
    mol.specieIndex = specieIndex;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);
    return mol;
}

/*! \brief Pretty-print the Involvement enum as an int for verbose output. */
int cfsc_involvement_as_int(Involvement inv)
{
    return static_cast<std::underlying_type<Involvement>::type>(inv);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a genuine state change is detected and flagged on BOTH sides.
// -----------------------------------------------------------------------------
void test_cfsc_detects_state_change()
{
    std::cerr << "\n[TEST] test_cfsc_detects_state_change\n"
              << "  Source file: src/parser/check_for_state_change.cpp\n"
              << "  Function:    check_for_state_change()\n"
              << "  Scenario:    reactant A(a~U) and product A(a~P) -- identical\n"
              << "               interfaces apart from the state character.\n"
              << "  Pass criteria:\n"
              << "     * parsedRxn.hasStateChange becomes true\n"
              << "     * reactant iface status  == Involvement::stateChange\n"
              << "     * product  iface status  == Involvement::stateChange\n"
              << "     * product.molTypeIndex is copied from the reactant molecule\n";

    // --- Build the reactant side -------------------------------------------
    ParsedMol::IfaceInfo reactIface
        = cfsc_make_iface("a", 'U', /*isBound*/ false, /*bondIndex*/ -1,
            /*speciesIndex*/ 0, /*molTypeIndex*/ 3, /*relIndex*/ 0, /*absIndex*/ 7);
    ParsedMol reactMol = cfsc_make_mol("A", /*molTypeIndex*/ 3, /*specieIndex*/ 0, reactIface);

    // --- Build the product side (same everything except the state) ---------
    ParsedMol::IfaceInfo prodIface
        = cfsc_make_iface("a", 'P', /*isBound*/ false, /*bondIndex*/ -1,
            /*speciesIndex*/ 0, /*molTypeIndex*/ 3, /*relIndex*/ 0, /*absIndex*/ 7);
    ParsedMol prodMol = cfsc_make_mol("A", /*molTypeIndex*/ -1, /*specieIndex*/ 0, prodIface);

    ParsedRxn parsedRxn;
    parsedRxn.hasStateChange = false;
    parsedRxn.productList.push_back(prodMol);

    // Sanity check on the helper predicate before invoking the function.
    std::cerr << "  Pre-check: areSameExceptState(reactIface, prodIface) = "
              << std::boolalpha
              << areSameExceptState(reactMol.interfaceList[0], parsedRxn.productList[0].interfaceList[0])
              << "\n";
    EXPECT_TRUE(areSameExceptState(reactMol.interfaceList[0], parsedRxn.productList[0].interfaceList[0]))
        << "Fixture is wrong: the two interfaces must differ only by state";

    std::cerr << "  Calling check_for_state_change...\n";
    check_for_state_change(reactMol.interfaceList[0], reactMol, parsedRxn);

    // --- Assertions ---------------------------------------------------------
    EXPECT_TRUE(parsedRxn.hasStateChange)
        << "hasStateChange must be set when a state change is detected";

    EXPECT_EQ(cfsc_involvement_as_int(reactMol.interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::stateChange))
        << "The reactant interface must be flagged as Involvement::stateChange";

    EXPECT_EQ(cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::stateChange))
        << "The product interface must be flagged as Involvement::stateChange";

    EXPECT_EQ(parsedRxn.productList[0].molTypeIndex, 3)
        << "The product molTypeIndex must be copied from the reactant molecule";

    std::cerr << "  Result: hasStateChange=" << std::boolalpha << parsedRxn.hasStateChange
              << ", reactant status=" << cfsc_involvement_as_int(reactMol.interfaceList[0].ifaceRxnStatus)
              << ", product status="
              << cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[0].ifaceRxnStatus)
              << ", product molTypeIndex=" << parsedRxn.productList[0].molTypeIndex << "\n";
}

// -----------------------------------------------------------------------------
// Test 2: matching molecule but IDENTICAL state -> no state change is reported,
//         yet the molTypeIndex bookkeeping still happens.
// -----------------------------------------------------------------------------
void test_cfsc_same_state_is_not_a_change()
{
    std::cerr << "\n[TEST] test_cfsc_same_state_is_not_a_change\n"
              << "  Source file: src/parser/check_for_state_change.cpp\n"
              << "  Function:    check_for_state_change()\n"
              << "  Scenario:    reactant A(a~U) and product A(a~U) -- the state is\n"
              << "               unchanged, so areSameExceptState() is false.\n"
              << "  Pass criteria:\n"
              << "     * hasStateChange stays false\n"
              << "     * both interface statuses remain Involvement::possible\n"
              << "     * product.molTypeIndex is STILL updated (bookkeeping runs\n"
              << "       for every name-matching product regardless of outcome)\n";

    ParsedMol::IfaceInfo reactIface
        = cfsc_make_iface("a", 'U', false, -1, 0, 5, 0, 11);
    ParsedMol reactMol = cfsc_make_mol("A", /*molTypeIndex*/ 5, /*specieIndex*/ 0, reactIface);

    // Product interface is a perfect copy - identical state included.
    ParsedMol::IfaceInfo prodIface
        = cfsc_make_iface("a", 'U', false, -1, 0, 5, 0, 11);
    ParsedMol prodMol = cfsc_make_mol("A", /*molTypeIndex*/ -1, /*specieIndex*/ 0, prodIface);

    ParsedRxn parsedRxn;
    parsedRxn.hasStateChange = false;
    parsedRxn.productList.push_back(prodMol);

    std::cerr << "  Pre-check: areSameExceptState(...) = " << std::boolalpha
              << areSameExceptState(reactMol.interfaceList[0], parsedRxn.productList[0].interfaceList[0])
              << " (expected false)\n";

    std::cerr << "  Calling check_for_state_change...\n";
    check_for_state_change(reactMol.interfaceList[0], reactMol, parsedRxn);

    EXPECT_FALSE(parsedRxn.hasStateChange)
        << "hasStateChange must stay false when the state is unchanged";

    EXPECT_EQ(cfsc_involvement_as_int(reactMol.interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::possible))
        << "The reactant interface status must be untouched";

    EXPECT_EQ(cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::possible))
        << "The product interface status must be untouched";

    // The molTypeIndex assignment sits outside the state-change branch.
    EXPECT_EQ(parsedRxn.productList[0].molTypeIndex, 5)
        << "molTypeIndex bookkeeping runs for every name-matching product";

    std::cerr << "  Result: hasStateChange=" << std::boolalpha << parsedRxn.hasStateChange
              << ", product molTypeIndex=" << parsedRxn.productList[0].molTypeIndex << "\n";
}

// -----------------------------------------------------------------------------
// Test 3: products belonging to a DIFFERENT molecule are skipped entirely.
// -----------------------------------------------------------------------------
void test_cfsc_ignores_other_molecules()
{
    std::cerr << "\n[TEST] test_cfsc_ignores_other_molecules\n"
              << "  Source file: src/parser/check_for_state_change.cpp\n"
              << "  Function:    check_for_state_change()\n"
              << "  Scenario:    the only product is molecule B while the reactant is\n"
              << "               molecule A, so the name comparison never matches.\n"
              << "  Pass criteria:\n"
              << "     * hasStateChange stays false\n"
              << "     * the product's molTypeIndex is NOT overwritten (-1 preserved)\n"
              << "     * neither interface status changes\n";

    ParsedMol::IfaceInfo reactIface
        = cfsc_make_iface("a", 'U', false, -1, 0, 2, 0, 4);
    ParsedMol reactMol = cfsc_make_mol("A", /*molTypeIndex*/ 2, /*specieIndex*/ 0, reactIface);

    // Same interface name and a different state, but on a molecule called "B".
    ParsedMol::IfaceInfo prodIface
        = cfsc_make_iface("a", 'P', false, -1, 0, 2, 0, 4);
    ParsedMol prodMol = cfsc_make_mol("B", /*molTypeIndex*/ -1, /*specieIndex*/ 0, prodIface);

    ParsedRxn parsedRxn;
    parsedRxn.hasStateChange = false;
    parsedRxn.productList.push_back(prodMol);

    std::cerr << "  Calling check_for_state_change...\n";
    check_for_state_change(reactMol.interfaceList[0], reactMol, parsedRxn);

    EXPECT_FALSE(parsedRxn.hasStateChange)
        << "A product on a different molecule must not trigger a state change";

    EXPECT_EQ(parsedRxn.productList[0].molTypeIndex, -1)
        << "molTypeIndex must be left alone for non-matching molecule names";

    EXPECT_EQ(cfsc_involvement_as_int(reactMol.interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::possible))
        << "The reactant interface status must be untouched";

    EXPECT_EQ(cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::possible))
        << "The product interface status must be untouched";

    std::cerr << "  Result: hasStateChange=" << std::boolalpha << parsedRxn.hasStateChange
              << ", product molTypeIndex=" << parsedRxn.productList[0].molTypeIndex << "\n";
}

// -----------------------------------------------------------------------------
// Test 4: an empty product list is handled gracefully (loop body never runs).
// -----------------------------------------------------------------------------
void test_cfsc_empty_product_list()
{
    std::cerr << "\n[TEST] test_cfsc_empty_product_list\n"
              << "  Source file: src/parser/check_for_state_change.cpp\n"
              << "  Function:    check_for_state_change()\n"
              << "  Scenario:    parsedRxn.productList is empty.\n"
              << "  Pass criteria: the call is a no-op -- nothing is flagged and\n"
              << "                 hasStateChange remains false.\n";

    ParsedMol::IfaceInfo reactIface
        = cfsc_make_iface("a", 'U', false, -1, 0, 0, 0, 0);
    ParsedMol reactMol = cfsc_make_mol("A", /*molTypeIndex*/ 0, /*specieIndex*/ 0, reactIface);

    ParsedRxn parsedRxn;
    parsedRxn.hasStateChange = false;
    parsedRxn.productList.clear();

    std::cerr << "  Calling check_for_state_change on an empty product list...\n";
    check_for_state_change(reactMol.interfaceList[0], reactMol, parsedRxn);

    EXPECT_FALSE(parsedRxn.hasStateChange)
        << "An empty product list cannot produce a state change";
    EXPECT_EQ(cfsc_involvement_as_int(reactMol.interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::possible))
        << "The reactant interface status must be untouched";
    EXPECT_TRUE(parsedRxn.productList.empty())
        << "The product list must remain empty";

    std::cerr << "  Result: no changes made, as expected.\n";
}

// -----------------------------------------------------------------------------
// Test 5: with several interfaces on the matching product, only the one that
//         truly corresponds to the target reactant interface is flagged.
// -----------------------------------------------------------------------------
void test_cfsc_flags_only_the_matching_interface()
{
    std::cerr << "\n[TEST] test_cfsc_flags_only_the_matching_interface\n"
              << "  Source file: src/parser/check_for_state_change.cpp\n"
              << "  Function:    check_for_state_change()\n"
              << "  Scenario:    product A carries two interfaces: 'a' (state 'U'->'P')\n"
              << "               and 'b' (state 'U' on both sides). The reactant\n"
              << "               interface under test is 'a'.\n"
              << "  Pass criteria:\n"
              << "     * only product interface 'a' gets Involvement::stateChange\n"
              << "     * product interface 'b' keeps Involvement::possible\n"
              << "     * hasStateChange becomes true\n";

    // Reactant interface under test: A(a~U)
    ParsedMol::IfaceInfo reactIfaceA
        = cfsc_make_iface("a", 'U', false, -1, 0, 4, 0, 21);
    ParsedMol reactMol = cfsc_make_mol("A", /*molTypeIndex*/ 4, /*specieIndex*/ 0, reactIfaceA);

    // Product molecule A with two interfaces.
    ParsedMol prodMol;
    prodMol.molName = "A";
    prodMol.molTypeIndex = -1;
    prodMol.specieIndex = 0;
    prodMol.interfaceList.clear();
    // 'a' differs only by state -> should be flagged.
    prodMol.interfaceList.push_back(cfsc_make_iface("a", 'P', false, -1, 0, 4, 0, 21));
    // 'b' has a different interface name -> should never match.
    prodMol.interfaceList.push_back(cfsc_make_iface("b", 'U', false, -1, 0, 4, 1, 22));

    ParsedRxn parsedRxn;
    parsedRxn.hasStateChange = false;
    parsedRxn.productList.push_back(prodMol);

    std::cerr << "  Calling check_for_state_change...\n";
    check_for_state_change(reactMol.interfaceList[0], reactMol, parsedRxn);

    EXPECT_TRUE(parsedRxn.hasStateChange)
        << "The 'a' interface state change must be detected";

    EXPECT_EQ(cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::stateChange))
        << "Product interface 'a' must be flagged as a state change";

    EXPECT_EQ(cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[1].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::possible))
        << "Product interface 'b' must be left as Involvement::possible";

    EXPECT_EQ(cfsc_involvement_as_int(reactMol.interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::stateChange))
        << "The reactant interface 'a' must be flagged as a state change";

    std::cerr << "  Result: iface 'a' status="
              << cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[0].ifaceRxnStatus)
              << ", iface 'b' status="
              << cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[1].ifaceRxnStatus)
              << "\n";
}

// -----------------------------------------------------------------------------
// Test 6: several products share the reactant's molecule name -- every one of
//         them receives the molTypeIndex, and each qualifying interface is
//         flagged (the function does not stop at the first hit).
// -----------------------------------------------------------------------------
void test_cfsc_multiple_matching_products()
{
    std::cerr << "\n[TEST] test_cfsc_multiple_matching_products\n"
              << "  Source file: src/parser/check_for_state_change.cpp\n"
              << "  Function:    check_for_state_change()\n"
              << "  Scenario:    productList = [A(a~P), B(a~P), A(a~P)].\n"
              << "  Pass criteria:\n"
              << "     * both 'A' products get molTypeIndex copied over\n"
              << "     * the 'B' product keeps its original molTypeIndex\n"
              << "     * both 'A' interfaces get Involvement::stateChange\n"
              << "     * hasStateChange becomes true\n";

    ParsedMol::IfaceInfo reactIface
        = cfsc_make_iface("a", 'U', false, -1, 0, 9, 0, 33);
    ParsedMol reactMol = cfsc_make_mol("A", /*molTypeIndex*/ 9, /*specieIndex*/ 0, reactIface);

    ParsedRxn parsedRxn;
    parsedRxn.hasStateChange = false;
    parsedRxn.productList.push_back(
        cfsc_make_mol("A", -1, 0, cfsc_make_iface("a", 'P', false, -1, 0, 9, 0, 33)));
    parsedRxn.productList.push_back(
        cfsc_make_mol("B", -7, 0, cfsc_make_iface("a", 'P', false, -1, 0, 9, 0, 33)));
    parsedRxn.productList.push_back(
        cfsc_make_mol("A", -1, 0, cfsc_make_iface("a", 'P', false, -1, 0, 9, 0, 33)));

    std::cerr << "  Calling check_for_state_change...\n";
    check_for_state_change(reactMol.interfaceList[0], reactMol, parsedRxn);

    EXPECT_TRUE(parsedRxn.hasStateChange)
        << "A state change exists on the matching 'A' products";

    // Bookkeeping on the two 'A' products.
    EXPECT_EQ(parsedRxn.productList[0].molTypeIndex, 9)
        << "First 'A' product must receive the reactant molTypeIndex";
    EXPECT_EQ(parsedRxn.productList[2].molTypeIndex, 9)
        << "Second 'A' product must receive the reactant molTypeIndex";
    EXPECT_EQ(parsedRxn.productList[1].molTypeIndex, -7)
        << "The 'B' product must keep its original molTypeIndex";

    // Flags on the interfaces.
    EXPECT_EQ(cfsc_involvement_as_int(parsedRxn.productList[0].interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::stateChange))
        << "First 'A' interface must be flagged";
    EXPECT_EQ(cfsc_involvement_as_int(parsedRxn.productList[2].interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::stateChange))
        << "Second 'A' interface must be flagged";
    EXPECT_EQ(cfsc_involvement_as_int(parsedRxn.productList[1].interfaceList[0].ifaceRxnStatus),
        cfsc_involvement_as_int(Involvement::possible))
        << "The 'B' interface must not be flagged";

    std::cerr << "  Result: molTypeIndex = [" << parsedRxn.productList[0].molTypeIndex << ", "
              << parsedRxn.productList[1].molTypeIndex << ", "
              << parsedRxn.productList[2].molTypeIndex << "]\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* helper runs inside its own TEST so the
// framework reports them individually while still executing every case.
// -----------------------------------------------------------------------------
TEST(CheckForStateChange, DetectsStateChange) { test_cfsc_detects_state_change(); }
TEST(CheckForStateChange, SameStateIsNotAChange) { test_cfsc_same_state_is_not_a_change(); }
TEST(CheckForStateChange, IgnoresOtherMolecules) { test_cfsc_ignores_other_molecules(); }
TEST(CheckForStateChange, EmptyProductList) { test_cfsc_empty_product_list(); }
TEST(CheckForStateChange, FlagsOnlyMatchingInterface) { test_cfsc_flags_only_the_matching_interface(); }
TEST(CheckForStateChange, MultipleMatchingProducts) { test_cfsc_multiple_matching_products(); }