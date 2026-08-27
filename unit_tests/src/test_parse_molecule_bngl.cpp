/*! \file test_parse_molecule_bngl.cpp
 *
 * ### Unit tests for src/parser/parse_molecule_bngl.cpp
 *
 * The file under test provides two free functions used while parsing the BNGL
 * style reaction / concentration blocks of a NERDSS input file:
 *
 *   1. `ParsedMol parse_molecule_bngl(int& totSpecies, bool isProductSide,
 *                                     std::pair<std::string,int> oneMol)`
 *      Turns a single molecule string such as `"A(a1,a2~u,a3!1)"` into a
 *      `ParsedMol` holding the molecule name and a list of `IfaceInfo` entries.
 *
 *   2. `ParsedMolNumState parse_number_bngl(std::string oneLine)`
 *      Turns a copy-number line such as `"100(a~U),200(a~P)"` into a
 *      `ParsedMolNumState` holding the total copy number, the copy number of
 *      each state and the name of each state.
 *
 * Important behaviours that the assertions below encode (read straight from the
 * implementation):
 *   - Interface *states* are up-cased (`~u` becomes `'U'`).
 *   - The 5-argument `IfaceInfo` constructor is used for *unbound* interfaces,
 *     so their `bondIndex` keeps its default value of -1.  The 6-argument
 *     constructor is used for *bound* interfaces, so wildcards get `bondIndex`
 *     0 and explicit bonds get the parsed integer.
 *   - On the reactant side only wildcard bonds (`!*`) are legal; an explicit
 *     bond index there calls `exit(1)`, therefore that path is NOT exercised
 *     (it would kill the whole gtest binary).
 *   - Any character that is not alphanumeric and not one of
 *     `_ ~ ( ! , )` also reaches an `exit(1)` default branch, so no spaces,
 *     dots or plus signs appear in the strings handed to the parser here.
 *   - `totSpecies` is taken by reference but is never modified.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

/*! \brief Pretty-print a ParsedMol to stderr so the test log shows what the
 *         parser actually produced.
 */
void pmb_dump_parsed_mol(const ParsedMol& mol)
{
    std::cerr << "    -> molName=\"" << mol.molName << "\""
              << ", specieIndex=" << mol.specieIndex
              << ", molTypeIndex=" << mol.molTypeIndex
              << ", #interfaces=" << mol.interfaceList.size() << '\n';
    for (const auto& iface : mol.interfaceList) {
        std::cerr << "       iface name=\"" << iface.ifaceName << "\""
                  << " fullName=\"" << iface.ifaceAndStateName << "\""
                  << " state='" << (iface.state == '\0' ? '0' : iface.state) << "'"
                  << " isBound=" << std::boolalpha << iface.isBound
                  << " bondIndex=" << iface.bondIndex
                  << " involvement(int)=" << static_cast<int>(iface.ifaceRxnStatus)
                  << " speciesIndex=" << iface.speciesIndex << '\n';
    }
}

/*! \brief Pretty-print a ParsedMolNumState to stderr. */
void pnb_dump_parsed_num(const ParsedMolNumState& num)
{
    std::cerr << "    -> totalCopyNumbers=" << num.totalCopyNumbers
              << ", #numberEachState=" << num.numberEachState.size()
              << ", #nameEachState=" << num.nameEachState.size() << '\n';
    for (size_t i = 0; i < num.numberEachState.size(); ++i)
        std::cerr << "       numberEachState[" << i << "] = " << num.numberEachState[i] << '\n';
    for (size_t i = 0; i < num.nameEachState.size(); ++i)
        std::cerr << "       nameEachState[" << i << "]   = \"" << num.nameEachState[i] << "\"\n";
}

} // namespace

// -----------------------------------------------------------------------------
// parse_molecule_bngl: a molecule with an empty interface list.
// -----------------------------------------------------------------------------
void parse_bngl_test_molecule_no_interfaces()
{
    std::cerr << "\n[TEST] parse_bngl_test_molecule_no_interfaces\n"
              << "  Source file: src/parser/parse_molecule_bngl.cpp\n"
              << "  Function:    parse_molecule_bngl()\n"
              << "  Input:       \"A()\" (reactant side, specieIndex 0)\n"
              << "  Criteria:    molName == \"A\", no interfaces parsed, and the\n"
              << "               by-reference totSpecies argument is untouched.\n";

    int totSpecies { 7 }; // the routine documents that it never changes this
    ParsedMol mol = parse_molecule_bngl(totSpecies, false, std::make_pair(std::string { "A()" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.molName, "A") << "The text before '(' becomes the molecule name";
    EXPECT_EQ(mol.interfaceList.size(), 0u) << "An empty parenthesis pair yields no interfaces";
    EXPECT_EQ(mol.specieIndex, 0) << "specieIndex is copied from the pair's second element";
    EXPECT_EQ(mol.molTypeIndex, -1) << "molTypeIndex is only resolved later, must stay at -1";
    EXPECT_EQ(totSpecies, 7) << "parse_molecule_bngl must not modify totSpecies";
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: several plain (unbound, stateless) interfaces.
// -----------------------------------------------------------------------------
void parse_bngl_test_unbound_interfaces()
{
    std::cerr << "\n[TEST] parse_bngl_test_unbound_interfaces\n"
              << "  Function:  parse_molecule_bngl()\n"
              << "  Input:     \"pip2(head,tail)\" (reactant side)\n"
              << "  Criteria:  two interfaces named head/tail, no state, not bound,\n"
              << "             bondIndex left at the default -1 (5-arg ctor),\n"
              << "             Involvement::possible for both.\n";

    int totSpecies { 0 };
    ParsedMol mol
        = parse_molecule_bngl(totSpecies, false, std::make_pair(std::string { "pip2(head,tail)" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.molName, "pip2") << "Alphanumeric name before '(' is kept verbatim";
    EXPECT_EQ(mol.interfaceList.size(), 2u) << "Comma and closing paren each flush one interface";

    if (mol.interfaceList.size() >= 2) {
        const auto& head = mol.interfaceList[0];
        const auto& tail = mol.interfaceList[1];

        EXPECT_EQ(head.ifaceName, "head") << "First interface name";
        EXPECT_EQ(head.state, '\0') << "No '~' means no state character";
        EXPECT_FALSE(head.isBound) << "No '!' means the interface is free";
        EXPECT_EQ(head.bondIndex, -1) << "Unbound ifaces use the 5-arg ctor so bondIndex stays -1";
        EXPECT_EQ(static_cast<int>(head.ifaceRxnStatus), static_cast<int>(Involvement::possible))
            << "Free interfaces start out as 'possible'";
        EXPECT_EQ(head.ifaceAndStateName, "head") << "With no state the full name equals the name";

        EXPECT_EQ(tail.ifaceName, "tail") << "Second interface name (flushed by ')')";
        EXPECT_FALSE(tail.isBound) << "Second interface is also free";
        EXPECT_EQ(static_cast<int>(tail.ifaceRxnStatus), static_cast<int>(Involvement::possible))
            << "Second interface is also 'possible'";
    }
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: '~' introduces a state and the state is up-cased.
// -----------------------------------------------------------------------------
void parse_bngl_test_state_is_uppercased()
{
    std::cerr << "\n[TEST] parse_bngl_test_state_is_uppercased\n"
              << "  Function:  parse_molecule_bngl()\n"
              << "  Input:     \"A(a~u)\" (reactant side)\n"
              << "  Criteria:  one interface named \"a\", state stored as 'U'\n"
              << "             (std::toupper is applied), unbound, 'possible',\n"
              << "             and ifaceAndStateName == \"a~U\".\n";

    int totSpecies { 0 };
    ParsedMol mol = parse_molecule_bngl(totSpecies, false, std::make_pair(std::string { "A(a~u)" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.molName, "A") << "Molecule name";
    EXPECT_EQ(mol.interfaceList.size(), 1u) << "One state-bearing interface expected";

    if (!mol.interfaceList.empty()) {
        const auto& iface = mol.interfaceList[0];
        EXPECT_EQ(iface.ifaceName, "a") << "Only the text before '~' is the interface name";
        EXPECT_EQ(iface.state, 'U') << "The parser up-cases the state character";
        EXPECT_FALSE(iface.isBound) << "No bond indicator present";
        EXPECT_EQ(iface.bondIndex, -1) << "Unbound iface keeps default bondIndex";
        EXPECT_EQ(static_cast<int>(iface.ifaceRxnStatus), static_cast<int>(Involvement::possible))
            << "A stated but unbound iface is 'possible'";
        EXPECT_EQ(iface.ifaceAndStateName, "a~U")
            << "IfaceInfo ctor builds name+'~'+state when a state exists";
    }
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: reactant side wildcard bond, no state -> ancillary.
// -----------------------------------------------------------------------------
void parse_bngl_test_reactant_wildcard_bond()
{
    std::cerr << "\n[TEST] parse_bngl_test_reactant_wildcard_bond\n"
              << "  Function:  parse_molecule_bngl() with isProductSide == false\n"
              << "  Input:     \"A(a!*)\"\n"
              << "  Criteria:  interface \"a\" is bound, bondIndex forced to 0 and\n"
              << "             its Involvement is 'ancillary' (it cannot react).\n"
              << "  Note:      an explicit bond index on the reactant side calls\n"
              << "             exit(1), so that branch is deliberately not tested.\n";

    int totSpecies { 0 };
    ParsedMol mol = parse_molecule_bngl(totSpecies, false, std::make_pair(std::string { "A(a!*)" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.interfaceList.size(), 1u) << "One bound interface expected";
    if (!mol.interfaceList.empty()) {
        const auto& iface = mol.interfaceList[0];
        EXPECT_EQ(iface.ifaceName, "a") << "Interface name is the text before '!'";
        EXPECT_TRUE(iface.isBound) << "'!' marks the interface as bound";
        EXPECT_EQ(iface.bondIndex, 0) << "Wildcard bonds are stored with bondIndex 0";
        EXPECT_EQ(iface.state, '\0') << "No state was given";
        EXPECT_EQ(static_cast<int>(iface.ifaceRxnStatus), static_cast<int>(Involvement::ancillary))
            << "Reactant wildcard bond without a state is ancillary";
    }
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: reactant side wildcard bond WITH a state -> possible.
// -----------------------------------------------------------------------------
void parse_bngl_test_reactant_wildcard_with_state()
{
    std::cerr << "\n[TEST] parse_bngl_test_reactant_wildcard_with_state\n"
              << "  Function:  parse_molecule_bngl() with isProductSide == false\n"
              << "  Input:     \"A(a~p!*)\"\n"
              << "  Criteria:  the '~' branch does NOT emit an interface because the\n"
              << "             next character is '!', so the '!' branch emits a single\n"
              << "             interface named \"a\" with state 'P', bound, bondIndex 0\n"
              << "             and Involvement::possible (its state may change).\n";

    int totSpecies { 0 };
    ParsedMol mol = parse_molecule_bngl(totSpecies, false, std::make_pair(std::string { "A(a~p!*)" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.interfaceList.size(), 1u) << "Exactly one interface must be produced, not two";
    if (!mol.interfaceList.empty()) {
        const auto& iface = mol.interfaceList[0];
        EXPECT_EQ(iface.ifaceName, "a") << "The '~state' suffix is stripped from the name";
        EXPECT_EQ(iface.state, 'P') << "State is recovered from the buffer and up-cased";
        EXPECT_TRUE(iface.isBound) << "Wildcard bond means bound";
        EXPECT_EQ(iface.bondIndex, 0) << "Wildcard bond index is 0";
        EXPECT_EQ(static_cast<int>(iface.ifaceRxnStatus), static_cast<int>(Involvement::possible))
            << "A stated wildcard-bound reactant iface is 'possible'";
        EXPECT_EQ(iface.ifaceAndStateName, "a~P") << "Full name is rebuilt as name~STATE";
    }
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: product side explicit bond index -> interactionChange.
// -----------------------------------------------------------------------------
void parse_bngl_test_product_indexed_bond()
{
    std::cerr << "\n[TEST] parse_bngl_test_product_indexed_bond\n"
              << "  Function:  parse_molecule_bngl() with isProductSide == true\n"
              << "  Input:     \"A(a!1)\"\n"
              << "  Criteria:  bound interface with bondIndex parsed from the digits\n"
              << "             (1) and Involvement::interactionChange.\n";

    int totSpecies { 0 };
    ParsedMol mol = parse_molecule_bngl(totSpecies, true, std::make_pair(std::string { "A(a!1)" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.interfaceList.size(), 1u) << "One bound product interface expected";
    if (!mol.interfaceList.empty()) {
        const auto& iface = mol.interfaceList[0];
        EXPECT_EQ(iface.ifaceName, "a") << "Interface name";
        EXPECT_TRUE(iface.isBound) << "Product interface is bound";
        EXPECT_EQ(iface.bondIndex, 1) << "std::stoi of the characters after '!' gives 1";
        EXPECT_EQ(iface.state, '\0') << "No state was declared";
        EXPECT_EQ(static_cast<int>(iface.ifaceRxnStatus), static_cast<int>(Involvement::interactionChange))
            << "An indexed product bond marks an interaction change";
    }
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: product side explicit bond index plus a state.
// -----------------------------------------------------------------------------
void parse_bngl_test_product_indexed_bond_with_state()
{
    std::cerr << "\n[TEST] parse_bngl_test_product_indexed_bond_with_state\n"
              << "  Function:  parse_molecule_bngl() with isProductSide == true\n"
              << "  Input:     \"A(a~p!2)\"\n"
              << "  Criteria:  single interface \"a\", state 'P', bound, bondIndex 2,\n"
              << "             Involvement::interactionChange.\n";

    int totSpecies { 0 };
    ParsedMol mol = parse_molecule_bngl(totSpecies, true, std::make_pair(std::string { "A(a~p!2)" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.interfaceList.size(), 1u) << "The '~' branch must not emit a duplicate interface";
    if (!mol.interfaceList.empty()) {
        const auto& iface = mol.interfaceList[0];
        EXPECT_EQ(iface.ifaceName, "a") << "State suffix stripped from the name";
        EXPECT_EQ(iface.state, 'P') << "State up-cased and preserved through the bond branch";
        EXPECT_TRUE(iface.isBound) << "Bound because of '!'";
        EXPECT_EQ(iface.bondIndex, 2) << "Bond index parsed as 2";
        EXPECT_EQ(static_cast<int>(iface.ifaceRxnStatus), static_cast<int>(Involvement::interactionChange))
            << "Stated + indexed bond on the product side is an interaction change";
        EXPECT_EQ(iface.ifaceAndStateName, "a~P") << "Full name rebuilt as a~P";
    }
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: product side wildcard bonds (with and without state).
// -----------------------------------------------------------------------------
void parse_bngl_test_product_wildcard_bond()
{
    std::cerr << "\n[TEST] parse_bngl_test_product_wildcard_bond\n"
              << "  Function:  parse_molecule_bngl() with isProductSide == true\n"
              << "  Inputs:    \"A(a!*)\"   -> ancillary, bondIndex 0\n"
              << "             \"A(a~p!*)\" -> possible,  bondIndex 0, state 'P'\n"
              << "  Criteria:  a wildcard bond on the product side cannot have been\n"
              << "             created by this reaction, so it is never an\n"
              << "             interaction change.\n";

    int totSpecies { 0 };

    ParsedMol plain = parse_molecule_bngl(totSpecies, true, std::make_pair(std::string { "A(a!*)" }, 0));
    std::cerr << "  Case 1: \"A(a!*)\"\n";
    pmb_dump_parsed_mol(plain);
    EXPECT_EQ(plain.interfaceList.size(), 1u) << "One interface expected";
    if (!plain.interfaceList.empty()) {
        EXPECT_TRUE(plain.interfaceList[0].isBound) << "Wildcard means bound";
        EXPECT_EQ(plain.interfaceList[0].bondIndex, 0) << "Wildcard bondIndex is 0";
        EXPECT_EQ(static_cast<int>(plain.interfaceList[0].ifaceRxnStatus),
            static_cast<int>(Involvement::ancillary))
            << "Stateless product wildcard is ancillary";
    }

    ParsedMol stated = parse_molecule_bngl(totSpecies, true, std::make_pair(std::string { "A(a~p!*)" }, 0));
    std::cerr << "  Case 2: \"A(a~p!*)\"\n";
    pmb_dump_parsed_mol(stated);
    EXPECT_EQ(stated.interfaceList.size(), 1u) << "One interface expected";
    if (!stated.interfaceList.empty()) {
        EXPECT_EQ(stated.interfaceList[0].ifaceName, "a") << "Name without the state suffix";
        EXPECT_EQ(stated.interfaceList[0].state, 'P') << "State kept and up-cased";
        EXPECT_EQ(stated.interfaceList[0].bondIndex, 0) << "Wildcard bondIndex is 0";
        EXPECT_EQ(static_cast<int>(stated.interfaceList[0].ifaceRxnStatus),
            static_cast<int>(Involvement::possible))
            << "Stated product wildcard is 'possible' (its state may have changed)";
    }
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: a mixture of free / stated / bonded interfaces.
// -----------------------------------------------------------------------------
void parse_bngl_test_mixed_interface_list()
{
    std::cerr << "\n[TEST] parse_bngl_test_mixed_interface_list\n"
              << "  Function:  parse_molecule_bngl() with isProductSide == true\n"
              << "  Input:     \"clat(a,b~u,c!1)\"\n"
              << "  Criteria:  three interfaces are produced, in source order, with\n"
              << "             the correct per-interface flags.\n";

    int totSpecies { 0 };
    ParsedMol mol
        = parse_molecule_bngl(totSpecies, true, std::make_pair(std::string { "clat(a,b~u,c!1)" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.molName, "clat") << "Molecule name comes from before the '('";
    EXPECT_EQ(mol.interfaceList.size(), 3u) << "Three interfaces separated by commas";

    if (mol.interfaceList.size() >= 3) {
        // 1: plain free interface flushed by the comma
        EXPECT_EQ(mol.interfaceList[0].ifaceName, "a") << "First interface name";
        EXPECT_FALSE(mol.interfaceList[0].isBound) << "First interface is free";
        EXPECT_EQ(mol.interfaceList[0].state, '\0') << "First interface has no state";

        // 2: stated free interface flushed by the '~' branch
        EXPECT_EQ(mol.interfaceList[1].ifaceName, "b") << "Second interface name";
        EXPECT_EQ(mol.interfaceList[1].state, 'U') << "Second interface state up-cased";
        EXPECT_FALSE(mol.interfaceList[1].isBound) << "Second interface is free";
        EXPECT_EQ(static_cast<int>(mol.interfaceList[1].ifaceRxnStatus),
            static_cast<int>(Involvement::possible))
            << "Stated free interface is 'possible'";

        // 3: bonded interface flushed by the '!' branch
        EXPECT_EQ(mol.interfaceList[2].ifaceName, "c") << "Third interface name";
        EXPECT_TRUE(mol.interfaceList[2].isBound) << "Third interface is bound";
        EXPECT_EQ(mol.interfaceList[2].bondIndex, 1) << "Third interface bond index";
        EXPECT_EQ(static_cast<int>(mol.interfaceList[2].ifaceRxnStatus),
            static_cast<int>(Involvement::interactionChange))
            << "Indexed product bond is an interaction change";
    }
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: the species index from the input pair is propagated.
// -----------------------------------------------------------------------------
void parse_bngl_test_species_index_propagation()
{
    std::cerr << "\n[TEST] parse_bngl_test_species_index_propagation\n"
              << "  Function:  parse_molecule_bngl()\n"
              << "  Input:     \"B(b1,b2~u)\" with specieIndex 1\n"
              << "  Criteria:  ParsedMol::specieIndex and every IfaceInfo::speciesIndex\n"
              << "             equal the second element of the input pair.\n";

    int totSpecies { 3 };
    ParsedMol mol = parse_molecule_bngl(totSpecies, false, std::make_pair(std::string { "B(b1,b2~u)" }, 1));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.specieIndex, 1) << "ParsedMol carries the species index it was given";
    EXPECT_EQ(mol.interfaceList.size(), 2u) << "Two interfaces expected";
    for (const auto& iface : mol.interfaceList) {
        EXPECT_EQ(iface.speciesIndex, 1)
            << "Every interface must inherit the same species index (" << iface.ifaceName << ")";
    }
    EXPECT_EQ(totSpecies, 3) << "totSpecies is still untouched";
}

// -----------------------------------------------------------------------------
// parse_molecule_bngl: documents the '_' case fall-through into the '~' case.
// -----------------------------------------------------------------------------
void parse_bngl_test_underscore_falls_through_to_state()
{
    std::cerr << "\n[TEST] parse_bngl_test_underscore_falls_through_to_state\n"
              << "  Function:  parse_molecule_bngl()\n"
              << "  Input:     \"A(a_b)\"\n"
              << "  Criteria:  this test DOCUMENTS current behaviour: the `case '_'`\n"
              << "             branch has no `break`, so it falls through into the\n"
              << "             `case '~'` branch. The result is one interface named\n"
              << "             \"a_\" whose state is the up-cased following character\n"
              << "             ('B'), not an interface literally named \"a_b\".\n";

    int totSpecies { 0 };
    ParsedMol mol = parse_molecule_bngl(totSpecies, false, std::make_pair(std::string { "A(a_b)" }, 0));
    pmb_dump_parsed_mol(mol);

    EXPECT_EQ(mol.molName, "A") << "Molecule name is unaffected by the fall-through";
    EXPECT_EQ(mol.interfaceList.size(), 1u) << "The fall-through emits exactly one interface";
    if (!mol.interfaceList.empty()) {
        EXPECT_EQ(mol.interfaceList[0].ifaceName, "a_")
            << "Name is the buffer up to and including the underscore";
        EXPECT_EQ(mol.interfaceList[0].state, 'B')
            << "The character after '_' is consumed as an up-cased state";
        EXPECT_FALSE(mol.interfaceList[0].isBound) << "Still an unbound interface";
    }
}

// -----------------------------------------------------------------------------
// parse_number_bngl: a bare copy number with no per-state breakdown.
// -----------------------------------------------------------------------------
void parse_bngl_test_number_plain_copy_number()
{
    std::cerr << "\n[TEST] parse_bngl_test_number_plain_copy_number\n"
              << "  Source file: src/parser/parse_molecule_bngl.cpp\n"
              << "  Function:    parse_number_bngl()\n"
              << "  Input:       \"500\"\n"
              << "  Criteria:    the trailing-buffer branch adds 500 to the total and\n"
              << "               both per-state vectors stay empty.\n";

    ParsedMolNumState num = parse_number_bngl(std::string { "500" });
    pnb_dump_parsed_num(num);

    EXPECT_EQ(num.totalCopyNumbers, 500) << "The whole string is one copy number";
    EXPECT_EQ(num.numberEachState.size(), 0u) << "No '(' was seen so no per-state counts";
    EXPECT_EQ(num.nameEachState.size(), 0u) << "No ')' was seen so no state names";
}

// -----------------------------------------------------------------------------
// parse_number_bngl: two explicit states, comma separated.
// -----------------------------------------------------------------------------
void parse_bngl_test_number_per_state_copy_numbers()
{
    std::cerr << "\n[TEST] parse_bngl_test_number_per_state_copy_numbers\n"
              << "  Function:  parse_number_bngl()\n"
              << "  Input:     \"100(a~U),200(a~P)\"\n"
              << "  Criteria:  totalCopyNumbers == 300, numberEachState == {100,200},\n"
              << "             nameEachState == {\"a~U\",\"a~P\"}.\n";

    ParsedMolNumState num = parse_number_bngl(std::string { "100(a~U),200(a~P)" });
    pnb_dump_parsed_num(num);

    EXPECT_EQ(num.totalCopyNumbers, 300) << "Copy numbers accumulate at every '('";
    EXPECT_EQ(num.numberEachState.size(), 2u) << "One entry pushed per '('";
    EXPECT_EQ(num.nameEachState.size(), 2u) << "One entry pushed per ')'";

    if (num.numberEachState.size() >= 2) {
        EXPECT_EQ(num.numberEachState[0], 100) << "First state copy number";
        EXPECT_EQ(num.numberEachState[1], 200) << "Second state copy number";
    }
    if (num.nameEachState.size() >= 2) {
        EXPECT_EQ(num.nameEachState[0], "a~U") << "First state name";
        EXPECT_EQ(num.nameEachState[1], "a~P") << "Second state name";
    }
}

// -----------------------------------------------------------------------------
// parse_number_bngl: state characters are up-cased here as well.
// -----------------------------------------------------------------------------
void parse_bngl_test_number_state_name_uppercased()
{
    std::cerr << "\n[TEST] parse_bngl_test_number_state_name_uppercased\n"
              << "  Function:  parse_number_bngl()\n"
              << "  Input:     \"50(x~y)\"\n"
              << "  Criteria:  the stored state name is \"x~Y\" (state up-cased) and\n"
              << "             totalCopyNumbers == 50.\n";

    ParsedMolNumState num = parse_number_bngl(std::string { "50(x~y)" });
    pnb_dump_parsed_num(num);

    EXPECT_EQ(num.totalCopyNumbers, 50) << "Single copy number";
    EXPECT_EQ(num.numberEachState.size(), 1u) << "Exactly one state group";
    EXPECT_EQ(num.nameEachState.size(), 1u) << "Exactly one state name";
    if (!num.numberEachState.empty())
        EXPECT_EQ(num.numberEachState[0], 50) << "Copy number for the single state";
    if (!num.nameEachState.empty())
        EXPECT_EQ(num.nameEachState[0], "x~Y") << "std::toupper is applied to the state character";
}

// -----------------------------------------------------------------------------
// parse_number_bngl: several interfaces inside one parenthesis group.
// -----------------------------------------------------------------------------
void parse_bngl_test_number_multiple_states_in_one_group()
{
    std::cerr << "\n[TEST] parse_bngl_test_number_multiple_states_in_one_group\n"
              << "  Function:  parse_number_bngl()\n"
              << "  Input:     \"100(a~u,b~p)\"\n"
              << "  Criteria:  the comma is appended to the (non-empty) buffer, so the\n"
              << "             whole group is stored as the single name \"a~U,b~P\"\n"
              << "             with one count of 100.\n";

    ParsedMolNumState num = parse_number_bngl(std::string { "100(a~u,b~p)" });
    pnb_dump_parsed_num(num);

    EXPECT_EQ(num.totalCopyNumbers, 100) << "Only one '(' so only one count added";
    EXPECT_EQ(num.numberEachState.size(), 1u) << "One count pushed at the single '('";
    EXPECT_EQ(num.nameEachState.size(), 1u) << "One name pushed at the single ')'";
    if (!num.numberEachState.empty())
        EXPECT_EQ(num.numberEachState[0], 100) << "The parsed count for the group";
    if (!num.nameEachState.empty())
        EXPECT_EQ(num.nameEachState[0], "a~U,b~P")
            << "Interior commas are kept and both states are up-cased";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each named helper is invoked from its own TEST so that
// a failure in one does not prevent the remaining checks from running.
// -----------------------------------------------------------------------------
TEST(ParseMoleculeBngl, MoleculeWithNoInterfaces) { parse_bngl_test_molecule_no_interfaces(); }
TEST(ParseMoleculeBngl, UnboundInterfaces) { parse_bngl_test_unbound_interfaces(); }
TEST(ParseMoleculeBngl, StateIsUppercased) { parse_bngl_test_state_is_uppercased(); }
TEST(ParseMoleculeBngl, ReactantWildcardBond) { parse_bngl_test_reactant_wildcard_bond(); }
TEST(ParseMoleculeBngl, ReactantWildcardBondWithState) { parse_bngl_test_reactant_wildcard_with_state(); }
TEST(ParseMoleculeBngl, ProductIndexedBond) { parse_bngl_test_product_indexed_bond(); }
TEST(ParseMoleculeBngl, ProductIndexedBondWithState) { parse_bngl_test_product_indexed_bond_with_state(); }
TEST(ParseMoleculeBngl, ProductWildcardBond) { parse_bngl_test_product_wildcard_bond(); }
TEST(ParseMoleculeBngl, MixedInterfaceList) { parse_bngl_test_mixed_interface_list(); }
TEST(ParseMoleculeBngl, SpeciesIndexPropagation) { parse_bngl_test_species_index_propagation(); }
TEST(ParseMoleculeBngl, UnderscoreFallsThroughToState) { parse_bngl_test_underscore_falls_through_to_state(); }

TEST(ParseNumberBngl, PlainCopyNumber) { parse_bngl_test_number_plain_copy_number(); }
TEST(ParseNumberBngl, PerStateCopyNumbers) { parse_bngl_test_number_per_state_copy_numbers(); }
TEST(ParseNumberBngl, StateNameUppercased) { parse_bngl_test_number_state_name_uppercased(); }
TEST(ParseNumberBngl, MultipleStatesInOneGroup) { parse_bngl_test_number_multiple_states_in_one_group(); }