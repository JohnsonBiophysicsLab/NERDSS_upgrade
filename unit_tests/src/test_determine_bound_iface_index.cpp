/*! \file test_determine_bound_iface_index.cpp
 *
 * ### Unit tests for src/parser/determine_bound_iface_index.cpp
 *
 * The function under test is
 *
 * \code
 * void determine_bound_iface_index(int& totSpecies,
 *                                  ParsedMol::IfaceInfo& targIface,
 *                                  ParsedRxn& parsedRxn,
 *                                  const std::vector<ForwardRxn>& forwardRxns,
 *                                  const std::vector<MolTemplate>& molTemplateList);
 * \endcode
 *
 * Behaviour, straight from the implementation:
 *   - It walks every ParsedMol in `parsedRxn.productList` and every
 *     ParsedMol::IfaceInfo in that molecule's `interfaceList`.
 *   - The *first* product interface that is (a) NOT equal to `targIface` and
 *     (b) has the SAME `bondIndex` as `targIface` is treated as the partner
 *     of the newly formed bond.
 *   - When such a partner is found the function
 *       * increments `totSpecies` (a new bound species is created),
 *       * calls `change_ifaceRxnStatus(totSpecies, Involvement::interactionChange)`
 *         on BOTH interfaces (which also sets `indexFound = true` and
 *         `absIndex = totSpecies`),
 *       * back-fills `prodIface.molTypeIndex` from the parent ParsedMol when it
 *         is still the default -1,
 *       * appends `totSpecies` to `parsedRxn.intProductList`,
 *       * appends the product interface FIRST and the target interface SECOND
 *         to `parsedRxn.rxnProducts`,
 *       * and returns immediately (only one partner is ever processed).
 *   - When no partner is found nothing at all is modified.
 *
 * Note on equality: `ParsedMol::IfaceInfo::operator!=` only compares
 * `absIndex`, `state`, `isBound` and `bondIndex` -- names and molTypeIndex are
 * ignored. The tests below are written against those exact semantics.
 *
 * The `forwardRxns` and `molTemplateList` arguments are never read by the
 * function, so empty vectors are passed for them.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (unique "dbii_" prefix so they cannot collide with other tests)
// -----------------------------------------------------------------------------
namespace {

/*! \brief Convert an Involvement enumerator to its underlying int.
 *
 * We deliberately do NOT stream Involvement values directly: the overloaded
 * operator<< for Involvement calls exit() for enumerators it does not know
 * about (e.g. Involvement::none), which would tear down the whole test binary.
 */
int dbii_involvement_value(Involvement inv) { return static_cast<int>(inv); }

/*! \brief Build a fully initialised ParsedMol::IfaceInfo.
 *
 * Every field that `operator!=` inspects (absIndex, state, isBound, bondIndex)
 * is set explicitly so the tests are not at the mercy of default values.
 */
ParsedMol::IfaceInfo dbii_make_iface(const std::string& name, char state, bool isBound, int bondIndex, int absIndex,
    int relIndex, int molTypeIndex, Involvement status)
{
    // The 5-argument constructor sets name/state/isBound/bondIndex/status and
    // composes ifaceAndStateName for us.
    ParsedMol::IfaceInfo info(name, state, isBound, bondIndex, status);
    info.absIndex = absIndex;
    info.relIndex = relIndex;
    info.molTypeIndex = molTypeIndex;
    return info;
}

/*! \brief Build a ParsedMol shell (no interfaces yet). */
ParsedMol dbii_make_mol(const std::string& molName, int molTypeIndex, int specieIndex)
{
    ParsedMol mol;
    mol.molName = molName;
    mol.molTypeIndex = molTypeIndex;
    mol.specieIndex = specieIndex;
    return mol;
}

/*! \brief Dump an interface to stderr so failures are easy to diagnose. */
void dbii_print_iface(const std::string& label, const ParsedMol::IfaceInfo& iface)
{
    std::cerr << "    " << label << ": name=" << iface.ifaceName << ", absIndex=" << iface.absIndex
              << ", relIndex=" << iface.relIndex << ", molTypeIndex=" << iface.molTypeIndex
              << ", bondIndex=" << iface.bondIndex << ", isBound=" << std::boolalpha << iface.isBound
              << ", indexFound=" << iface.indexFound
              << ", involvement(int)=" << dbii_involvement_value(iface.ifaceRxnStatus) << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a product interface sharing the bond index is found -> a brand new
//         bound species is created and both interfaces are updated.
// -----------------------------------------------------------------------------
void test_dbii_matching_partner_creates_new_species()
{
    std::cerr << "\n[TEST] test_dbii_matching_partner_creates_new_species\n"
              << "  Source file: determine_bound_iface_index.cpp\n"
              << "  Function:    determine_bound_iface_index()\n"
              << "  Scenario:    target interface (bondIndex 1) has exactly one\n"
              << "               distinct partner in the product list with bondIndex 1.\n"
              << "  Pass criteria: totSpecies incremented once; both interfaces get\n"
              << "               absIndex == new totSpecies, indexFound == true and\n"
              << "               Involvement::interactionChange; intProductList gains the\n"
              << "               new index; rxnProducts gains {product, target} in that order.\n";

    int totSpecies = 7; // pretend 7 species already exist in the system

    // --- build the parsed reaction's product side: B(b1!1) ---
    ParsedRxn parsedRxn;
    ParsedMol prodMol = dbii_make_mol("B", 1, 1);
    prodMol.interfaceList.push_back(
        dbii_make_iface("b1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 20, /*relIndex*/ 0,
            /*molTypeIndex*/ 1, Involvement::possible));
    parsedRxn.productList.push_back(prodMol);

    // --- the interface whose bound index we are trying to determine: A(a1!1) ---
    ParsedMol::IfaceInfo targIface = dbii_make_iface("a1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 10,
        /*relIndex*/ 0, /*molTypeIndex*/ 0, Involvement::possible);

    // unused-by-the-function arguments
    std::vector<ForwardRxn> forwardRxns;
    std::vector<MolTemplate> molTemplateList;

    std::cerr << "  Before the call (totSpecies = " << totSpecies << "):\n";
    dbii_print_iface("targIface", targIface);
    dbii_print_iface("prodIface", parsedRxn.productList[0].interfaceList[0]);

    determine_bound_iface_index(totSpecies, targIface, parsedRxn, forwardRxns, molTemplateList);

    std::cerr << "  After the call (totSpecies = " << totSpecies << "):\n";
    dbii_print_iface("targIface", targIface);
    dbii_print_iface("prodIface", parsedRxn.productList[0].interfaceList[0]);

    // The species counter must have advanced by exactly one.
    EXPECT_EQ(totSpecies, 8) << "One new bound species should have been created";

    // The target interface must now point at the new species and be flagged as
    // an interaction-changing interface.
    EXPECT_EQ(targIface.absIndex, 8) << "targIface.absIndex should be the new species index";
    EXPECT_TRUE(targIface.indexFound) << "change_ifaceRxnStatus() sets indexFound";
    EXPECT_EQ(dbii_involvement_value(targIface.ifaceRxnStatus),
        dbii_involvement_value(Involvement::interactionChange))
        << "targIface should be marked as interactionChange";

    // ...and so must the matched product interface (modified in place).
    const ParsedMol::IfaceInfo& updatedProd = parsedRxn.productList[0].interfaceList[0];
    EXPECT_EQ(updatedProd.absIndex, 8) << "product interface should share the new species index";
    EXPECT_TRUE(updatedProd.indexFound) << "product interface should have indexFound set";
    EXPECT_EQ(dbii_involvement_value(updatedProd.ifaceRxnStatus),
        dbii_involvement_value(Involvement::interactionChange))
        << "product interface should be marked as interactionChange";

    // The new species index is recorded in the integer product list.
    EXPECT_EQ(parsedRxn.intProductList.size(), size_t { 1 }) << "exactly one product index should be recorded";
    if (parsedRxn.intProductList.size() == 1)
        EXPECT_EQ(parsedRxn.intProductList[0], 8) << "recorded product index should be the new species index";

    // rxnProducts gets the product interface first, then the target interface.
    EXPECT_EQ(parsedRxn.rxnProducts.size(), size_t { 2 }) << "two product interfaces should be recorded";
    if (parsedRxn.rxnProducts.size() == 2) {
        EXPECT_EQ(parsedRxn.rxnProducts[0].ifaceName, std::string { "b1" })
            << "the matched product interface must be pushed first";
        EXPECT_EQ(parsedRxn.rxnProducts[1].ifaceName, std::string { "a1" })
            << "the target interface must be pushed second";
        EXPECT_EQ(parsedRxn.rxnProducts[0].absIndex, 8) << "stored copy must carry the new index";
        EXPECT_EQ(parsedRxn.rxnProducts[1].absIndex, 8) << "stored copy must carry the new index";
    }
}

// -----------------------------------------------------------------------------
// Test 2: no product interface shares the target's bond index -> nothing at all
//         is touched.
// -----------------------------------------------------------------------------
void test_dbii_no_matching_bond_index_is_a_no_op()
{
    std::cerr << "\n[TEST] test_dbii_no_matching_bond_index_is_a_no_op\n"
              << "  Source file: determine_bound_iface_index.cpp\n"
              << "  Function:    determine_bound_iface_index()\n"
              << "  Scenario:    target has bondIndex 1 while the only product interface\n"
              << "               carries bondIndex 2.\n"
              << "  Pass criteria: totSpecies, both interfaces and both product lists are\n"
              << "               left exactly as they were.\n";

    int totSpecies = 42;

    ParsedRxn parsedRxn;
    ParsedMol prodMol = dbii_make_mol("B", 1, 1);
    prodMol.interfaceList.push_back(
        dbii_make_iface("b1", '\0', true, /*bondIndex*/ 2, /*absIndex*/ 20, 0, 1, Involvement::possible));
    parsedRxn.productList.push_back(prodMol);

    ParsedMol::IfaceInfo targIface
        = dbii_make_iface("a1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 10, 0, 0, Involvement::possible);

    std::vector<ForwardRxn> forwardRxns;
    std::vector<MolTemplate> molTemplateList;

    determine_bound_iface_index(totSpecies, targIface, parsedRxn, forwardRxns, molTemplateList);

    std::cerr << "  After the call (totSpecies = " << totSpecies << "):\n";
    dbii_print_iface("targIface", targIface);
    dbii_print_iface("prodIface", parsedRxn.productList[0].interfaceList[0]);

    EXPECT_EQ(totSpecies, 42) << "totSpecies must not change when no partner is found";
    EXPECT_EQ(targIface.absIndex, 10) << "targIface.absIndex must be untouched";
    EXPECT_FALSE(targIface.indexFound) << "targIface.indexFound must remain false";
    EXPECT_EQ(dbii_involvement_value(targIface.ifaceRxnStatus), dbii_involvement_value(Involvement::possible))
        << "targIface involvement must remain 'possible'";

    const ParsedMol::IfaceInfo& prodIface = parsedRxn.productList[0].interfaceList[0];
    EXPECT_EQ(prodIface.absIndex, 20) << "product interface absIndex must be untouched";
    EXPECT_FALSE(prodIface.indexFound) << "product interface indexFound must remain false";

    EXPECT_TRUE(parsedRxn.intProductList.empty()) << "intProductList must stay empty";
    EXPECT_TRUE(parsedRxn.rxnProducts.empty()) << "rxnProducts must stay empty";
}

// -----------------------------------------------------------------------------
// Test 3: an interface that compares EQUAL to the target is skipped, and the
//         next distinct interface with the same bond index is used instead.
// -----------------------------------------------------------------------------
void test_dbii_identical_interface_is_skipped()
{
    std::cerr << "\n[TEST] test_dbii_identical_interface_is_skipped\n"
              << "  Source file: determine_bound_iface_index.cpp\n"
              << "  Function:    determine_bound_iface_index()\n"
              << "  Scenario:    the first product interface is identical to the target\n"
              << "               (same absIndex/state/isBound/bondIndex); a second, distinct\n"
              << "               interface also carries bondIndex 1.\n"
              << "  Pass criteria: the identical interface is left alone and the second one\n"
              << "               becomes the reaction partner.\n";

    int totSpecies = 3;

    ParsedRxn parsedRxn;
    ParsedMol prodMol = dbii_make_mol("A", 0, 0);
    // interface 0: an exact clone of the target -> operator!= is false -> skipped
    prodMol.interfaceList.push_back(
        dbii_make_iface("a1clone", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 10, 0, 0, Involvement::possible));
    // interface 1: distinct absIndex but the same bondIndex -> this is the partner
    prodMol.interfaceList.push_back(
        dbii_make_iface("b1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 20, 1, 0, Involvement::possible));
    parsedRxn.productList.push_back(prodMol);

    ParsedMol::IfaceInfo targIface
        = dbii_make_iface("a1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 10, 0, 0, Involvement::possible);

    std::vector<ForwardRxn> forwardRxns;
    std::vector<MolTemplate> molTemplateList;

    determine_bound_iface_index(totSpecies, targIface, parsedRxn, forwardRxns, molTemplateList);

    const ParsedMol::IfaceInfo& clone = parsedRxn.productList[0].interfaceList[0];
    const ParsedMol::IfaceInfo& partner = parsedRxn.productList[0].interfaceList[1];

    std::cerr << "  After the call (totSpecies = " << totSpecies << "):\n";
    dbii_print_iface("clone   ", clone);
    dbii_print_iface("partner ", partner);

    EXPECT_EQ(totSpecies, 4) << "one new species should have been created";

    // The clone must be untouched: it compares equal to targIface so it is skipped.
    EXPECT_EQ(clone.absIndex, 10) << "the identical interface must not be re-indexed";
    EXPECT_FALSE(clone.indexFound) << "the identical interface must not be flagged";
    EXPECT_EQ(dbii_involvement_value(clone.ifaceRxnStatus), dbii_involvement_value(Involvement::possible))
        << "the identical interface must keep its original involvement";

    // The second interface is the one that becomes the bound partner.
    EXPECT_EQ(partner.absIndex, 4) << "the distinct interface should get the new species index";
    EXPECT_TRUE(partner.indexFound) << "the distinct interface should be flagged as found";
    EXPECT_EQ(dbii_involvement_value(partner.ifaceRxnStatus),
        dbii_involvement_value(Involvement::interactionChange))
        << "the distinct interface should be marked interactionChange";

    EXPECT_EQ(parsedRxn.rxnProducts.size(), size_t { 2 }) << "exactly one pair should be recorded";
    if (parsedRxn.rxnProducts.size() == 2)
        EXPECT_EQ(parsedRxn.rxnProducts[0].ifaceName, std::string { "b1" })
            << "the recorded partner must be the second (distinct) interface";
}

// -----------------------------------------------------------------------------
// Test 4: a product interface with molTypeIndex == -1 inherits the parent
//         ParsedMol's molTypeIndex before being stored.
// -----------------------------------------------------------------------------
void test_dbii_backfills_missing_moltype_index()
{
    std::cerr << "\n[TEST] test_dbii_backfills_missing_moltype_index\n"
              << "  Source file: determine_bound_iface_index.cpp\n"
              << "  Function:    determine_bound_iface_index()\n"
              << "  Scenario:    the matched product interface still has the default\n"
              << "               molTypeIndex of -1 while its parent molecule has index 3.\n"
              << "  Pass criteria: the interface (and the copy stored in rxnProducts) end up\n"
              << "               with molTypeIndex == 3.\n";

    int totSpecies = 0;

    ParsedRxn parsedRxn;
    ParsedMol prodMol = dbii_make_mol("C", /*molTypeIndex*/ 3, /*specieIndex*/ 1);
    prodMol.interfaceList.push_back(
        dbii_make_iface("c1", '\0', true, /*bondIndex*/ 5, /*absIndex*/ 20, 0, /*molTypeIndex*/ -1,
            Involvement::possible));
    parsedRxn.productList.push_back(prodMol);

    ParsedMol::IfaceInfo targIface
        = dbii_make_iface("a1", '\0', true, /*bondIndex*/ 5, /*absIndex*/ 10, 0, 0, Involvement::possible);

    std::vector<ForwardRxn> forwardRxns;
    std::vector<MolTemplate> molTemplateList;

    determine_bound_iface_index(totSpecies, targIface, parsedRxn, forwardRxns, molTemplateList);

    const ParsedMol::IfaceInfo& partner = parsedRxn.productList[0].interfaceList[0];
    std::cerr << "  After the call (totSpecies = " << totSpecies << "):\n";
    dbii_print_iface("partner", partner);

    EXPECT_EQ(totSpecies, 1) << "the species counter starts at 0 and becomes 1";
    EXPECT_EQ(partner.molTypeIndex, 3) << "missing molTypeIndex must be copied from the parent ParsedMol";

    EXPECT_EQ(parsedRxn.rxnProducts.size(), size_t { 2 }) << "one interface pair should be recorded";
    if (parsedRxn.rxnProducts.size() == 2)
        EXPECT_EQ(parsedRxn.rxnProducts[0].molTypeIndex, 3)
            << "the stored copy is taken after the molTypeIndex fix-up";
}

// -----------------------------------------------------------------------------
// Test 5: the function returns after the first match -- later candidates,
//         even in other product molecules, are ignored.
// -----------------------------------------------------------------------------
void test_dbii_returns_after_first_match()
{
    std::cerr << "\n[TEST] test_dbii_returns_after_first_match\n"
              << "  Source file: determine_bound_iface_index.cpp\n"
              << "  Function:    determine_bound_iface_index()\n"
              << "  Scenario:    two different product molecules both expose an interface\n"
              << "               with the target's bondIndex.\n"
              << "  Pass criteria: only the first candidate is consumed, totSpecies is bumped\n"
              << "               once, and the second candidate is untouched.\n";

    int totSpecies = 100;

    ParsedRxn parsedRxn;

    ParsedMol firstMol = dbii_make_mol("B", 1, 1);
    firstMol.interfaceList.push_back(
        dbii_make_iface("b1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 20, 0, 1, Involvement::possible));
    parsedRxn.productList.push_back(firstMol);

    ParsedMol secondMol = dbii_make_mol("C", 2, 2);
    secondMol.interfaceList.push_back(
        dbii_make_iface("c1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 30, 0, 2, Involvement::possible));
    parsedRxn.productList.push_back(secondMol);

    ParsedMol::IfaceInfo targIface
        = dbii_make_iface("a1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 10, 0, 0, Involvement::possible);

    std::vector<ForwardRxn> forwardRxns;
    std::vector<MolTemplate> molTemplateList;

    determine_bound_iface_index(totSpecies, targIface, parsedRxn, forwardRxns, molTemplateList);

    const ParsedMol::IfaceInfo& first = parsedRxn.productList[0].interfaceList[0];
    const ParsedMol::IfaceInfo& second = parsedRxn.productList[1].interfaceList[0];

    std::cerr << "  After the call (totSpecies = " << totSpecies << "):\n";
    dbii_print_iface("first candidate ", first);
    dbii_print_iface("second candidate", second);

    EXPECT_EQ(totSpecies, 101) << "the counter must advance by exactly one";

    EXPECT_EQ(first.absIndex, 101) << "the first candidate is the one that reacts";
    EXPECT_TRUE(first.indexFound) << "the first candidate should be flagged";

    EXPECT_EQ(second.absIndex, 30) << "the second candidate must be left alone";
    EXPECT_FALSE(second.indexFound) << "the second candidate must not be flagged";
    EXPECT_EQ(dbii_involvement_value(second.ifaceRxnStatus), dbii_involvement_value(Involvement::possible))
        << "the second candidate keeps its original involvement";

    EXPECT_EQ(parsedRxn.intProductList.size(), size_t { 1 }) << "only one product index is added";
    EXPECT_EQ(parsedRxn.rxnProducts.size(), size_t { 2 }) << "only one interface pair is added";
}

// -----------------------------------------------------------------------------
// Test 6: the matching test uses bondIndex ONLY -- it never inspects isBound.
//         Two "unbound" interfaces (bondIndex == -1) that differ in some other
//         field are therefore treated as partners. This documents the current
//         behaviour of the code rather than an intended physical rule.
// -----------------------------------------------------------------------------
void test_dbii_matches_on_bond_index_only()
{
    std::cerr << "\n[TEST] test_dbii_matches_on_bond_index_only\n"
              << "  Source file: determine_bound_iface_index.cpp\n"
              << "  Function:    determine_bound_iface_index()\n"
              << "  Scenario:    target and product interface both have bondIndex -1 (the\n"
              << "               'unbound' default) but differ in their state character.\n"
              << "  Pass criteria: the function still pairs them, proving the match relies\n"
              << "               purely on bondIndex equality plus interface inequality.\n";

    int totSpecies = 0;

    ParsedRxn parsedRxn;
    ParsedMol prodMol = dbii_make_mol("B", 1, 1);
    prodMol.interfaceList.push_back(
        dbii_make_iface("b1", 'P', /*isBound*/ false, /*bondIndex*/ -1, /*absIndex*/ 20, 0, 1,
            Involvement::possible));
    parsedRxn.productList.push_back(prodMol);

    // Same bondIndex (-1) but a different state character -> operator!= is true.
    ParsedMol::IfaceInfo targIface = dbii_make_iface("a1", 'U', /*isBound*/ false, /*bondIndex*/ -1,
        /*absIndex*/ 10, 0, 0, Involvement::possible);

    std::vector<ForwardRxn> forwardRxns;
    std::vector<MolTemplate> molTemplateList;

    determine_bound_iface_index(totSpecies, targIface, parsedRxn, forwardRxns, molTemplateList);

    std::cerr << "  After the call (totSpecies = " << totSpecies << "):\n";
    dbii_print_iface("targIface", targIface);
    dbii_print_iface("prodIface", parsedRxn.productList[0].interfaceList[0]);

    EXPECT_EQ(totSpecies, 1) << "matching is done on bondIndex alone, so a species is created";
    EXPECT_EQ(targIface.absIndex, 1) << "target picks up the new species index";
    EXPECT_EQ(parsedRxn.productList[0].interfaceList[0].absIndex, 1)
        << "product picks up the same new species index";
    EXPECT_EQ(parsedRxn.rxnProducts.size(), size_t { 2 }) << "the pair is recorded in rxnProducts";
}

// -----------------------------------------------------------------------------
// Test 7: an empty product list must be handled gracefully (no crash, no work).
// -----------------------------------------------------------------------------
void test_dbii_empty_product_list_is_safe()
{
    std::cerr << "\n[TEST] test_dbii_empty_product_list_is_safe\n"
              << "  Source file: determine_bound_iface_index.cpp\n"
              << "  Function:    determine_bound_iface_index()\n"
              << "  Scenario:    parsedRxn.productList is empty.\n"
              << "  Pass criteria: the call returns without touching anything.\n";

    int totSpecies = 5;

    ParsedRxn parsedRxn; // no products at all
    ParsedMol::IfaceInfo targIface
        = dbii_make_iface("a1", '\0', true, /*bondIndex*/ 1, /*absIndex*/ 10, 0, 0, Involvement::possible);

    std::vector<ForwardRxn> forwardRxns;
    std::vector<MolTemplate> molTemplateList;

    determine_bound_iface_index(totSpecies, targIface, parsedRxn, forwardRxns, molTemplateList);

    std::cerr << "  After the call (totSpecies = " << totSpecies << "):\n";
    dbii_print_iface("targIface", targIface);

    EXPECT_EQ(totSpecies, 5) << "no species may be created with an empty product list";
    EXPECT_EQ(targIface.absIndex, 10) << "target interface must be untouched";
    EXPECT_FALSE(targIface.indexFound) << "target interface must not be flagged";
    EXPECT_TRUE(parsedRxn.intProductList.empty()) << "intProductList must stay empty";
    EXPECT_TRUE(parsedRxn.rxnProducts.empty()) << "rxnProducts must stay empty";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each scenario runs independently so a failure in one
// does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(DetermineBoundIfaceIndex, MatchingPartnerCreatesNewSpecies) { test_dbii_matching_partner_creates_new_species(); }
TEST(DetermineBoundIfaceIndex, NoMatchingBondIndexIsANoOp) { test_dbii_no_matching_bond_index_is_a_no_op(); }
TEST(DetermineBoundIfaceIndex, IdenticalInterfaceIsSkipped) { test_dbii_identical_interface_is_skipped(); }
TEST(DetermineBoundIfaceIndex, BackfillsMissingMolTypeIndex) { test_dbii_backfills_missing_moltype_index(); }
TEST(DetermineBoundIfaceIndex, ReturnsAfterFirstMatch) { test_dbii_returns_after_first_match(); }
TEST(DetermineBoundIfaceIndex, MatchesOnBondIndexOnly) { test_dbii_matches_on_bond_index_only(); }
TEST(DetermineBoundIfaceIndex, EmptyProductListIsSafe) { test_dbii_empty_product_list_is_safe(); }