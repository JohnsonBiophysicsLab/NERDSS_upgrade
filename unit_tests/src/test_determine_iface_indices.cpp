/*! \file test_determine_iface_indices.cpp
 *
 * ### Unit test for src/parser/determine_iface_indices.cpp
 *
 * Function under test:
 * \code
 *   void determine_iface_indices(int specieIndex, int& totSpecies, ParsedMol& targMol,
 *                                ParsedRxn& parsedRxn,
 *                                const std::vector<ForwardRxn>& forwardRxns,
 *                                const std::vector<MolTemplate>& molTemplateList);
 * \endcode
 *
 * The routine walks the list of MolTemplates looking for the template whose
 * `molName` matches the ParsedMol being parsed.  For every interface of the
 * ParsedMol that has a name matching a template interface it:
 *
 *   1. copies `targMol.molTypeIndex` onto the interface (note: it copies the
 *      value stored on the *ParsedMol*, NOT the one on the MolTemplate),
 *   2. stores the `specieIndex` argument as the interface's `speciesIndex`,
 *   3. stores the template interface's relative index as `relIndex`,
 *   4. and -- only when the absolute index is still unknown (`absIndex == -1`) --
 *      resolves the absolute (state) index:
 *        * template interface with exactly one state  -> that state's index,
 *        * template interface with several states and an explicit state char
 *          -> the matching state's index,
 *        * several states but no explicit state char  -> for association type
 *          reactions the reaction is flagged `willBeMultipleRxns = true`; for
 *          creation/destruction reactions a warning is printed instead,
 *        * a bound, non-wildcard interface additionally delegates to
 *          determine_bound_iface_index().
 *
 * NOTE ON COVERAGE: the `isBound && bondIndex != 0` branch calls
 * determine_bound_iface_index(), which terminates the process (exit(1)) when the
 * bond partner cannot be resolved from the list of previously parsed reactions.
 * Because a fatal exit would kill the whole gtest binary, that branch is *not*
 * exercised here; instead we verify the complementary case (wildcard bond,
 * `bondIndex == 0`) where the delegation is deliberately skipped.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Small builders shared by all the tests below.  Everything is fully
// initialized, since the parser reads members the caller is expected to have
// filled in (names, relative indices, state lists, ...).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a template Interface with an explicit relative index and states. */
Interface dii_make_interface(const std::string& name, int relIndex, std::vector<Interface::State> stateList)
{
    Interface iface; // default ctor: index = -1, empty name/state list
    iface.name = name;
    iface.index = relIndex; // relative index inside MolTemplate::interfaceList
    iface.stateList = std::move(stateList);
    iface.set_ifaceAndStateNames(); // only relevant when there is > 1 state
    return iface;
}

/*! \brief Build a MolTemplate carrying the given interfaces. */
MolTemplate dii_make_template(const std::string& molName, int molTypeIndex, std::vector<Interface> ifaces)
{
    MolTemplate temp;
    temp.molName = molName;
    temp.molTypeIndex = molTypeIndex;
    temp.interfaceList = std::move(ifaces);
    return temp;
}

/*!
 * \brief Build the MolTemplate list used by every test.
 *
 * Template "a" (molTypeIndex 0):
 *   - interface "a1", relative index 0, ONE unnamed state with absolute index 0
 *   - interface "a2", relative index 1, TWO states: 'U' (abs 1) and 'P' (abs 2)
 * Template "b" (molTypeIndex 1):
 *   - interface "b1", relative index 0, ONE unnamed state with absolute index 3
 */
std::vector<MolTemplate> dii_make_template_list()
{
    std::vector<Interface> aIfaces;
    aIfaces.push_back(dii_make_interface("a1", 0, { Interface::State(0) }));
    aIfaces.push_back(dii_make_interface("a2", 1, { Interface::State('U', 1), Interface::State('P', 2) }));

    std::vector<Interface> bIfaces;
    bIfaces.push_back(dii_make_interface("b1", 0, { Interface::State(3) }));

    std::vector<MolTemplate> list;
    list.push_back(dii_make_template("a", 0, std::move(aIfaces)));
    list.push_back(dii_make_template("b", 1, std::move(bIfaces)));
    return list;
}

/*!
 * \brief Build a parsed interface exactly as the BNGL parser would hand it over:
 *        name/state/bond information known, all indices still at their -1
 *        "unknown" defaults.
 */
ParsedMol::IfaceInfo dii_make_iface_info(const std::string& name, char state, bool isBound, int bondIndex)
{
    ParsedMol::IfaceInfo info; // defaults: absIndex = relIndex = molTypeIndex = speciesIndex = -1
    info.ifaceName = name;
    info.state = state;
    info.isBound = isBound;
    info.bondIndex = bondIndex;
    return info;
}

/*! \brief Build a ParsedMol with the supplied name/type index and interfaces. */
ParsedMol dii_make_parsed_mol(const std::string& molName, int molTypeIndex, std::vector<ParsedMol::IfaceInfo> ifaces)
{
    ParsedMol mol;
    mol.molName = molName;
    mol.molTypeIndex = molTypeIndex;
    mol.interfaceList = std::move(ifaces);
    return mol;
}

/*! \brief Build a ParsedRxn of a given reaction type, everything else default. */
ParsedRxn dii_make_parsed_rxn(ReactionType type)
{
    ParsedRxn rxn; // willBeMultipleRxns defaults to false
    rxn.rxnType = type;
    return rxn;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: an interface with a single state gets its absolute index from that
//         single state, and all bookkeeping indices are filled in.
// -----------------------------------------------------------------------------
void test_dii_single_state_interface()
{
    std::cerr << "\n[TEST] test_dii_single_state_interface\n"
              << "  Source file: src/parser/determine_iface_indices.cpp\n"
              << "  Function:    determine_iface_indices()\n"
              << "  Scenario:    mol \"a\", interface \"a1\" which has exactly one state.\n"
              << "  Pass:        absIndex == 0 (the state index), relIndex == 0,\n"
              << "               molTypeIndex/speciesIndex copied in, totSpecies untouched.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{}; // no previously parsed reactions needed
    ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
    int totSpecies = 5;

    // molTypeIndex deliberately set to 7 (not the template's 0) to show which
    // value the parser actually propagates onto the interface.
    ParsedMol mol{ dii_make_parsed_mol("a", 7, { dii_make_iface_info("a1", '\0', false, -1) }) };

    std::cerr << "  Calling determine_iface_indices(specieIndex=3, ...)\n";
    determine_iface_indices(3, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

    const ParsedMol::IfaceInfo& iface = mol.interfaceList.front();
    std::cerr << "  Result: absIndex=" << iface.absIndex << " relIndex=" << iface.relIndex
              << " molTypeIndex=" << iface.molTypeIndex << " speciesIndex=" << iface.speciesIndex << '\n';

    EXPECT_EQ(iface.absIndex, 0) << "single-state interface must take stateList[0].index";
    EXPECT_EQ(iface.relIndex, 0) << "relative index must come from the template interface";
    EXPECT_EQ(iface.molTypeIndex, 7)
        << "molTypeIndex is copied from the ParsedMol (7), not from the MolTemplate (0)";
    EXPECT_EQ(iface.speciesIndex, 3) << "speciesIndex must be the specieIndex argument";
    EXPECT_FALSE(parsedRxn.willBeMultipleRxns) << "a fully determined state must not split the reaction";
    EXPECT_EQ(totSpecies, 5) << "totSpecies is only modified through determine_bound_iface_index()";
}

// -----------------------------------------------------------------------------
// Test 2: a multi-state interface with an explicitly named state resolves to
//         that state's absolute index.
// -----------------------------------------------------------------------------
void test_dii_multi_state_explicit_state()
{
    std::cerr << "\n[TEST] test_dii_multi_state_explicit_state\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: interface \"a2\" has states 'U'(abs 1) and 'P'(abs 2);\n"
              << "            the parsed molecule names each of them in turn.\n"
              << "  Pass:     absIndex matches the named state, relIndex == 1.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};

    // --- state 'U' -> absolute index 1 ---
    {
        ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
        int totSpecies = 0;
        ParsedMol mol{ dii_make_parsed_mol("a", 0, { dii_make_iface_info("a2", 'U', false, -1) }) };

        determine_iface_indices(0, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

        std::cerr << "  state 'U' -> absIndex=" << mol.interfaceList.front().absIndex << '\n';
        EXPECT_EQ(mol.interfaceList.front().absIndex, 1) << "state 'U' has absolute index 1";
        EXPECT_EQ(mol.interfaceList.front().relIndex, 1) << "\"a2\" is relative index 1 on the template";
        EXPECT_FALSE(parsedRxn.willBeMultipleRxns) << "an explicit state must not split the reaction";
    }

    // --- state 'P' -> absolute index 2 ---
    {
        ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
        int totSpecies = 0;
        ParsedMol mol{ dii_make_parsed_mol("a", 0, { dii_make_iface_info("a2", 'P', false, -1) }) };

        determine_iface_indices(0, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

        std::cerr << "  state 'P' -> absIndex=" << mol.interfaceList.front().absIndex << '\n';
        EXPECT_EQ(mol.interfaceList.front().absIndex, 2) << "state 'P' has absolute index 2";
        EXPECT_EQ(mol.interfaceList.front().relIndex, 1) << "\"a2\" is relative index 1 on the template";
        EXPECT_FALSE(parsedRxn.willBeMultipleRxns);
    }
}

// -----------------------------------------------------------------------------
// Test 3: a multi-state interface with an *unknown* state character leaves the
//         absolute index unresolved (the state search simply finds nothing).
// -----------------------------------------------------------------------------
void test_dii_multi_state_unknown_state_char()
{
    std::cerr << "\n[TEST] test_dii_multi_state_unknown_state_char\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: interface \"a2\" asked for state 'X', which the template\n"
              << "            does not define.\n"
              << "  Pass:     absIndex stays -1, willBeMultipleRxns stays false,\n"
              << "            but relIndex/speciesIndex are still filled in.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};
    ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
    int totSpecies = 0;

    ParsedMol mol{ dii_make_parsed_mol("a", 0, { dii_make_iface_info("a2", 'X', false, -1) }) };

    determine_iface_indices(2, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

    const ParsedMol::IfaceInfo& iface = mol.interfaceList.front();
    std::cerr << "  Result: absIndex=" << iface.absIndex << " relIndex=" << iface.relIndex
              << " speciesIndex=" << iface.speciesIndex << '\n';

    EXPECT_EQ(iface.absIndex, -1) << "no template state matches 'X', so absIndex stays unknown";
    EXPECT_EQ(iface.relIndex, 1) << "relIndex is assigned before the absolute-index search";
    EXPECT_EQ(iface.speciesIndex, 2) << "speciesIndex is assigned before the absolute-index search";
    EXPECT_FALSE(parsedRxn.willBeMultipleRxns)
        << "the split flag is only raised when NO state character was given";
}

// -----------------------------------------------------------------------------
// Test 4: a multi-state interface with no state given in an association type
//         reaction flags the reaction for splitting into multiple reactions.
// -----------------------------------------------------------------------------
void test_dii_multi_state_no_state_flags_split()
{
    std::cerr << "\n[TEST] test_dii_multi_state_no_state_flags_split\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: bimolecular reaction naming \"a2\" with no explicit state.\n"
              << "  Pass:     parsedRxn.willBeMultipleRxns becomes true and absIndex\n"
              << "            is left at -1 for later resolution.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};
    ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
    int totSpecies = 11;

    ParsedMol mol{ dii_make_parsed_mol("a", 0, { dii_make_iface_info("a2", '\0', false, -1) }) };

    ASSERT_FALSE(parsedRxn.willBeMultipleRxns) << "sanity: the flag starts out false";
    determine_iface_indices(1, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

    std::cerr << "  Result: willBeMultipleRxns=" << std::boolalpha << parsedRxn.willBeMultipleRxns
              << " absIndex=" << mol.interfaceList.front().absIndex << '\n';

    EXPECT_TRUE(parsedRxn.willBeMultipleRxns)
        << "an unspecified state on a multi-state interface must request a reaction split";
    EXPECT_EQ(mol.interfaceList.front().absIndex, -1) << "absIndex cannot be resolved yet";
    EXPECT_EQ(mol.interfaceList.front().relIndex, 1) << "relIndex is still assigned";
    EXPECT_EQ(totSpecies, 11) << "totSpecies must not change on this path";
}

// -----------------------------------------------------------------------------
// Test 5: the same "no state given" situation in a creation/destruction reaction
//         only prints a warning -- the split flag must stay down.
// -----------------------------------------------------------------------------
void test_dii_creation_rxn_no_state_only_warns()
{
    std::cerr << "\n[TEST] test_dii_creation_rxn_no_state_only_warns\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: zerothOrderCreation / uniMolCreation / destruction reactions\n"
              << "            naming \"a2\" without a state.  (A WARNING line printed by\n"
              << "            the parser on stderr below is expected and is not a failure.)\n"
              << "  Pass:     willBeMultipleRxns stays false and absIndex stays -1.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};

    // All three reaction types take the "warn instead of split" branch.
    const ReactionType typesToCheck[]
        = { ReactionType::zerothOrderCreation, ReactionType::uniMolCreation, ReactionType::destruction };

    for (const ReactionType type : typesToCheck) {
        ParsedRxn parsedRxn{ dii_make_parsed_rxn(type) };
        int totSpecies = 4;
        ParsedMol mol{ dii_make_parsed_mol("a", 0, { dii_make_iface_info("a2", '\0', false, -1) }) };

        determine_iface_indices(0, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

        EXPECT_FALSE(parsedRxn.willBeMultipleRxns)
            << "creation/destruction reactions warn instead of splitting the reaction";
        EXPECT_EQ(mol.interfaceList.front().absIndex, -1) << "the default state is used later, not here";
        EXPECT_EQ(totSpecies, 4) << "totSpecies must not change on this path";
    }
}

// -----------------------------------------------------------------------------
// Test 6: an interface whose absolute index is already known is not re-resolved,
//         but its bookkeeping indices are still refreshed.
// -----------------------------------------------------------------------------
void test_dii_existing_absIndex_is_preserved()
{
    std::cerr << "\n[TEST] test_dii_existing_absIndex_is_preserved\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: interface \"a2\" already carries absIndex = 99 and no state\n"
              << "            character; the whole absolute-index block must be skipped.\n"
              << "  Pass:     absIndex stays 99, willBeMultipleRxns stays false, and\n"
              << "            relIndex/molTypeIndex/speciesIndex are still updated.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};
    ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
    int totSpecies = 0;

    ParsedMol::IfaceInfo preset{ dii_make_iface_info("a2", '\0', false, -1) };
    preset.absIndex = 99; // pretend a previous parsing pass already resolved this
    ParsedMol mol{ dii_make_parsed_mol("a", 5, { preset }) };

    determine_iface_indices(8, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

    const ParsedMol::IfaceInfo& iface = mol.interfaceList.front();
    std::cerr << "  Result: absIndex=" << iface.absIndex << " relIndex=" << iface.relIndex
              << " molTypeIndex=" << iface.molTypeIndex << " speciesIndex=" << iface.speciesIndex << '\n';

    EXPECT_EQ(iface.absIndex, 99) << "a known absolute index must never be overwritten";
    EXPECT_FALSE(parsedRxn.willBeMultipleRxns)
        << "the split flag lives inside the skipped block, so it must stay false";
    EXPECT_EQ(iface.relIndex, 1) << "relIndex is assigned outside the skipped block";
    EXPECT_EQ(iface.molTypeIndex, 5) << "molTypeIndex is assigned outside the skipped block";
    EXPECT_EQ(iface.speciesIndex, 8) << "speciesIndex is assigned outside the skipped block";
}

// -----------------------------------------------------------------------------
// Test 7: a bound interface with a wildcard bond (bondIndex == 0) resolves its
//         index normally and does NOT delegate to determine_bound_iface_index().
// -----------------------------------------------------------------------------
void test_dii_wildcard_bond_skips_bound_lookup()
{
    std::cerr << "\n[TEST] test_dii_wildcard_bond_skips_bound_lookup\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: interface \"a1\" is bound but with a wildcard bond index 0.\n"
              << "  Pass:     absIndex resolves from the single state and totSpecies is\n"
              << "            unchanged, proving determine_bound_iface_index() was skipped.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};
    ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
    int totSpecies = 42;

    // isBound == true, bondIndex == 0 -> wildcard, the bound lookup is skipped.
    ParsedMol mol{ dii_make_parsed_mol("a", 0, { dii_make_iface_info("a1", '\0', true, 0) }) };

    determine_iface_indices(0, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

    std::cerr << "  Result: absIndex=" << mol.interfaceList.front().absIndex << " totSpecies=" << totSpecies << '\n';

    EXPECT_EQ(mol.interfaceList.front().absIndex, 0) << "the single state still supplies the absolute index";
    EXPECT_EQ(mol.interfaceList.front().relIndex, 0) << "relative index of \"a1\" is 0";
    EXPECT_EQ(totSpecies, 42)
        << "a wildcard bond must not reach determine_bound_iface_index(), which would grow totSpecies";
}

// -----------------------------------------------------------------------------
// Test 8: a ParsedMol whose name matches no template is left completely alone.
// -----------------------------------------------------------------------------
void test_dii_unknown_molecule_name_is_ignored()
{
    std::cerr << "\n[TEST] test_dii_unknown_molecule_name_is_ignored\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: molecule \"zzz\" does not exist in the MolTemplate list.\n"
              << "  Pass:     every index on the parsed interface stays at its -1 default.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};
    ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
    int totSpecies = 0;

    ParsedMol mol{ dii_make_parsed_mol("zzz", 3, { dii_make_iface_info("a1", '\0', false, -1) }) };

    determine_iface_indices(6, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

    const ParsedMol::IfaceInfo& iface = mol.interfaceList.front();
    std::cerr << "  Result: absIndex=" << iface.absIndex << " relIndex=" << iface.relIndex
              << " molTypeIndex=" << iface.molTypeIndex << " speciesIndex=" << iface.speciesIndex << '\n';

    EXPECT_EQ(iface.absIndex, -1) << "no template matched, so nothing may be resolved";
    EXPECT_EQ(iface.relIndex, -1) << "no template matched, so relIndex stays unknown";
    EXPECT_EQ(iface.molTypeIndex, -1) << "no template matched, so molTypeIndex stays unknown";
    EXPECT_EQ(iface.speciesIndex, -1) << "no template matched, so speciesIndex stays unknown";
    EXPECT_FALSE(parsedRxn.willBeMultipleRxns);
}

// -----------------------------------------------------------------------------
// Test 9: an interface name that the matching template does not define is left
//         untouched, while its sibling interfaces are still processed.
// -----------------------------------------------------------------------------
void test_dii_unknown_interface_name_is_ignored()
{
    std::cerr << "\n[TEST] test_dii_unknown_interface_name_is_ignored\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: molecule \"a\" lists interfaces \"a1\" (valid) and\n"
              << "            \"nosuchiface\" (not on the template).\n"
              << "  Pass:     \"a1\" resolves, \"nosuchiface\" keeps all -1 defaults.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};
    ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
    int totSpecies = 0;

    ParsedMol mol{ dii_make_parsed_mol("a", 0,
        { dii_make_iface_info("a1", '\0', false, -1), dii_make_iface_info("nosuchiface", '\0', false, -1) }) };

    determine_iface_indices(0, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

    const ParsedMol::IfaceInfo& good = mol.interfaceList[0];
    const ParsedMol::IfaceInfo& bad = mol.interfaceList[1];
    std::cerr << "  Known iface \"a1\":        absIndex=" << good.absIndex << " relIndex=" << good.relIndex << '\n';
    std::cerr << "  Unknown iface \"nosuch\": absIndex=" << bad.absIndex << " relIndex=" << bad.relIndex << '\n';

    EXPECT_EQ(good.absIndex, 0) << "the recognised interface must still be resolved";
    EXPECT_EQ(good.relIndex, 0) << "the recognised interface must still get its relative index";
    EXPECT_EQ(bad.absIndex, -1) << "an unrecognised interface name must be left alone";
    EXPECT_EQ(bad.relIndex, -1) << "an unrecognised interface name must be left alone";
    EXPECT_EQ(bad.molTypeIndex, -1) << "an unrecognised interface name must be left alone";
    EXPECT_EQ(bad.speciesIndex, -1) << "an unrecognised interface name must be left alone";
}

// -----------------------------------------------------------------------------
// Test 10: several interfaces on the same molecule are all resolved in one call,
//          and a molecule matching the *second* template is found as well.
// -----------------------------------------------------------------------------
void test_dii_multiple_ifaces_and_second_template()
{
    std::cerr << "\n[TEST] test_dii_multiple_ifaces_and_second_template\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: (a) molecule \"a\" carrying both \"a1\" and \"a2~P\",\n"
              << "            (b) molecule \"b\" which matches the second MolTemplate.\n"
              << "  Pass:     all interfaces get the indices defined by their template.\n";

    const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
    const std::vector<ForwardRxn> forwardRxns{};

    // --- (a) two interfaces resolved in a single call ---
    {
        ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
        int totSpecies = 0;
        ParsedMol mol{ dii_make_parsed_mol(
            "a", 0, { dii_make_iface_info("a1", '\0', false, -1), dii_make_iface_info("a2", 'P', false, -1) }) };

        determine_iface_indices(0, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

        std::cerr << "  mol \"a\": a1 -> abs " << mol.interfaceList[0].absIndex << ", rel "
                  << mol.interfaceList[0].relIndex << " | a2~P -> abs " << mol.interfaceList[1].absIndex << ", rel "
                  << mol.interfaceList[1].relIndex << '\n';

        EXPECT_EQ(mol.interfaceList[0].absIndex, 0) << "\"a1\" single state -> absolute index 0";
        EXPECT_EQ(mol.interfaceList[0].relIndex, 0) << "\"a1\" is relative index 0";
        EXPECT_EQ(mol.interfaceList[1].absIndex, 2) << "\"a2\" state 'P' -> absolute index 2";
        EXPECT_EQ(mol.interfaceList[1].relIndex, 1) << "\"a2\" is relative index 1";
        EXPECT_FALSE(parsedRxn.willBeMultipleRxns);
    }

    // --- (b) the search continues past the first template ---
    {
        ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
        int totSpecies = 0;
        ParsedMol mol{ dii_make_parsed_mol("b", 1, { dii_make_iface_info("b1", '\0', false, -1) }) };

        determine_iface_indices(1, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

        std::cerr << "  mol \"b\": b1 -> abs " << mol.interfaceList[0].absIndex << ", rel "
                  << mol.interfaceList[0].relIndex << ", speciesIndex " << mol.interfaceList[0].speciesIndex << '\n';

        EXPECT_EQ(mol.interfaceList[0].absIndex, 3) << "\"b1\" single state -> absolute index 3";
        EXPECT_EQ(mol.interfaceList[0].relIndex, 0) << "\"b1\" is relative index 0 on template \"b\"";
        EXPECT_EQ(mol.interfaceList[0].molTypeIndex, 1) << "molTypeIndex comes from the ParsedMol (1)";
        EXPECT_EQ(mol.interfaceList[0].speciesIndex, 1) << "speciesIndex is the specieIndex argument (1)";
    }
}

// -----------------------------------------------------------------------------
// Test 11: an empty interface list (and an empty template list) must be handled
//          gracefully -- no crash, no state changes.
// -----------------------------------------------------------------------------
void test_dii_empty_inputs_are_safe()
{
    std::cerr << "\n[TEST] test_dii_empty_inputs_are_safe\n"
              << "  Function: determine_iface_indices()\n"
              << "  Scenario: (a) parsed molecule with zero interfaces,\n"
              << "            (b) completely empty MolTemplate list.\n"
              << "  Pass:     the call returns without touching totSpecies or the split flag.\n";

    const std::vector<ForwardRxn> forwardRxns{};

    // (a) molecule without any interfaces, but with a matching template
    {
        const std::vector<MolTemplate> molTemplateList{ dii_make_template_list() };
        ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
        int totSpecies = 13;
        ParsedMol mol{ dii_make_parsed_mol("a", 0, {}) };

        determine_iface_indices(0, totSpecies, mol, parsedRxn, forwardRxns, molTemplateList);

        EXPECT_TRUE(mol.interfaceList.empty()) << "the interface list must remain empty";
        EXPECT_EQ(totSpecies, 13) << "nothing to resolve, so totSpecies must not change";
        EXPECT_FALSE(parsedRxn.willBeMultipleRxns);
    }

    // (b) no templates at all
    {
        const std::vector<MolTemplate> emptyTemplates{};
        ParsedRxn parsedRxn{ dii_make_parsed_rxn(ReactionType::bimolecular) };
        int totSpecies = 13;
        ParsedMol mol{ dii_make_parsed_mol("a", 0, { dii_make_iface_info("a1", '\0', false, -1) }) };

        determine_iface_indices(0, totSpecies, mol, parsedRxn, forwardRxns, emptyTemplates);

        EXPECT_EQ(mol.interfaceList.front().absIndex, -1) << "with no templates nothing can be resolved";
        EXPECT_EQ(mol.interfaceList.front().relIndex, -1) << "with no templates nothing can be resolved";
        EXPECT_EQ(totSpecies, 13) << "totSpecies must not change";
        EXPECT_FALSE(parsedRxn.willBeMultipleRxns);
    }

    std::cerr << "  Both empty-input cases completed without crashing.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each named helper runs inside its own TEST so that a
// failure in one case does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(DetermineIfaceIndicesTest, SingleStateInterface) { test_dii_single_state_interface(); }
TEST(DetermineIfaceIndicesTest, MultiStateExplicitState) { test_dii_multi_state_explicit_state(); }
TEST(DetermineIfaceIndicesTest, MultiStateUnknownStateChar) { test_dii_multi_state_unknown_state_char(); }
TEST(DetermineIfaceIndicesTest, MultiStateNoStateFlagsSplit) { test_dii_multi_state_no_state_flags_split(); }
TEST(DetermineIfaceIndicesTest, CreationRxnNoStateOnlyWarns) { test_dii_creation_rxn_no_state_only_warns(); }
TEST(DetermineIfaceIndicesTest, ExistingAbsIndexPreserved) { test_dii_existing_absIndex_is_preserved(); }
TEST(DetermineIfaceIndicesTest, WildcardBondSkipsBoundLookup) { test_dii_wildcard_bond_skips_bound_lookup(); }
TEST(DetermineIfaceIndicesTest, UnknownMoleculeNameIgnored) { test_dii_unknown_molecule_name_is_ignored(); }
TEST(DetermineIfaceIndicesTest, UnknownInterfaceNameIgnored) { test_dii_unknown_interface_name_is_ignored(); }
TEST(DetermineIfaceIndicesTest, MultipleIfacesAndSecondTemplate) { test_dii_multiple_ifaces_and_second_template(); }
TEST(DetermineIfaceIndicesTest, EmptyInputsAreSafe) { test_dii_empty_inputs_are_safe(); }