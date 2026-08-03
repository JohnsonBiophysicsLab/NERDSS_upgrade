/*! \file test_init_speciesFile.cpp
 *
 * ### Unit test for src/io/init_speciesFile.cpp
 *
 * The single function under test is
 *
 *     int init_speciesFile(std::ofstream& speciesFile,
 *                          copyCounters& counterArrays,
 *                          std::vector<MolTemplate>& molTemplateList,
 *                          std::vector<ForwardRxn>& forwardRxns,
 *                          Parameters& params);
 *
 * It has two observable effects:
 *
 *   1. It writes a single CSV header line into `speciesFile`, beginning with
 *      "Time (s)" and followed by one column per tracked species:
 *        - "Mol(iface)"            for interfaces with exactly one state
 *        - "Mol(iface~S)"          for every state S of a multi-state interface
 *        - the ForwardRxn productName for every *bimolecular* forward reaction
 *   2. It grows the bookkeeping vectors inside `copyCounters`:
 *        - singleDouble    : 1 for monomer species, 2 for bimolecular products
 *        - implicitDouble  : true only if a bimolecular reactant is an implicit lipid
 *        - canDissociate   : true only for reversible bimolecular reactions
 *        - bindPairList    : one empty slot per species, *only* when
 *                            params.fromRestart == false
 *
 * The return value is the total number of species written.
 *
 * Every test below builds a minimal molTemplateList / forwardRxns pair,
 * runs the function against a scratch file, then checks both the header text
 * and the counter vectors. Verbose progress is printed to stderr so the
 * reader can follow which behaviour is being exercised.
 */

#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Declaration of the function under test (declared in include/io/io.hpp, but
// re-declared here so that this translation unit does not need to pull in the
// whole IO header chain).
// -----------------------------------------------------------------------------
int init_speciesFile(std::ofstream& speciesFile, copyCounters& counterArrays,
    std::vector<MolTemplate>& molTemplateList,
    std::vector<ForwardRxn>& forwardRxns, Parameters& params);

namespace {

/*! \brief Small container holding everything a single run produced. */
struct IsfRunResult {
    int nSpecies { 0 }; //!< value returned by init_speciesFile
    std::string header {}; //!< first (and only) line written to the file
};

/*! \brief Split a comma separated string into its individual tokens.
 *
 * Used to inspect the CSV header produced by init_speciesFile.
 */
std::vector<std::string> isf_split_csv(const std::string& line)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream iss(line);
    while (std::getline(iss, token, ','))
        tokens.push_back(token);
    return tokens;
}

/*! \brief Build a MolTemplate with the requested interfaces / states.
 *
 * \param[in] molName          name used in the header column, e.g. "A"
 * \param[in] molTypeIndex     index of this template inside molTemplateList
 * \param[in] ifaceNames       one entry per interface
 * \param[in] ifaceStates      state identifiers per interface. An empty inner
 *                             vector means "interface with NO states at all"
 *                             (an edge case handled by the function), a single
 *                             entry means a single-state interface.
 * \param[in] isImplicitLipid  marks the template as an implicit lipid
 */
MolTemplate isf_make_mol_template(const std::string& molName, int molTypeIndex,
    const std::vector<std::string>& ifaceNames,
    const std::vector<std::vector<char>>& ifaceStates,
    bool isImplicitLipid = false)
{
    MolTemplate temp;
    temp.molName = molName;
    temp.molTypeIndex = molTypeIndex;
    temp.isImplicitLipid = isImplicitLipid;

    for (size_t ifaceItr = 0; ifaceItr < ifaceNames.size(); ++ifaceItr) {
        Interface iface;
        iface.name = ifaceNames[ifaceItr];
        iface.index = static_cast<int>(ifaceItr);

        // Attach one Interface::State per requested state identifier.
        for (size_t stateItr = 0; stateItr < ifaceStates[ifaceItr].size(); ++stateItr) {
            iface.stateList.emplace_back(
                ifaceStates[ifaceItr][stateItr], static_cast<int>(stateItr));
        }
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Build a ForwardRxn suitable for the species file.
 *
 * \param[in] productName    string written into the header
 * \param[in] molTypeIndex1  molTemplateList index of the first reactant
 * \param[in] molTypeIndex2  molTemplateList index of the second reactant
 * \param[in] rxnType        reaction type (only bimolecular is counted)
 * \param[in] isReversible   controls the expected canDissociate flag
 */
ForwardRxn isf_make_forward_rxn(const std::string& productName, int molTypeIndex1,
    int molTypeIndex2, ReactionType rxnType, bool isReversible)
{
    ForwardRxn rxn;
    rxn.productName = productName;
    rxn.rxnType = rxnType;
    rxn.isReversible = isReversible;

    // The function indexes reactantListNew[0] and [1], so both must exist.
    RxnIface react1;
    react1.molTypeIndex = molTypeIndex1;
    RxnIface react2;
    react2.molTypeIndex = molTypeIndex2;
    rxn.reactantListNew.push_back(react1);
    rxn.reactantListNew.push_back(react2);

    return rxn;
}

/*! \brief Run init_speciesFile against a scratch file and read the header back.
 *
 * The scratch file is deleted before returning so repeated runs cannot
 * interfere with one another.
 */
IsfRunResult isf_run(copyCounters& counterArrays,
    std::vector<MolTemplate>& molTemplateList,
    std::vector<ForwardRxn>& forwardRxns, Parameters& params,
    const std::string& fileName)
{
    IsfRunResult result;

    std::ofstream speciesFile(fileName);
    // Sanity check: without an open stream the test cannot be meaningful.
    EXPECT_TRUE(speciesFile.is_open()) << "Could not open scratch file " << fileName;

    result.nSpecies = init_speciesFile(
        speciesFile, counterArrays, molTemplateList, forwardRxns, params);
    speciesFile.close();

    std::ifstream in(fileName);
    std::getline(in, result.header);
    in.close();
    std::remove(fileName.c_str());

    return result;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: no molecules and no reactions -> header only, nothing counted.
// -----------------------------------------------------------------------------
void test_isf_empty_inputs()
{
    std::cerr << "\n[TEST] test_isf_empty_inputs\n"
              << "  Source file:   src/io/init_speciesFile.cpp\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      empty molTemplateList and empty forwardRxns.\n"
              << "  Pass criteria: returns 0, header is exactly \"Time (s)\",\n"
              << "                 and all copyCounters vectors stay empty.\n";

    copyCounters counterArrays;
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    Parameters params;
    params.fromRestart = false;

    IsfRunResult res = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_empty.csv");

    std::cerr << "  header  = \"" << res.header << "\"\n"
              << "  species = " << res.nSpecies << '\n';

    EXPECT_EQ(res.nSpecies, 0) << "No molecules/reactions means zero species";
    EXPECT_EQ(res.header, "Time (s)") << "Header must contain only the time column";
    EXPECT_TRUE(counterArrays.singleDouble.empty()) << "singleDouble should stay empty";
    EXPECT_TRUE(counterArrays.implicitDouble.empty()) << "implicitDouble should stay empty";
    EXPECT_TRUE(counterArrays.canDissociate.empty()) << "canDissociate should stay empty";
    EXPECT_TRUE(counterArrays.bindPairList.empty()) << "bindPairList should stay empty";
}

// -----------------------------------------------------------------------------
// Test 2: single-state interfaces across two molecule templates.
// -----------------------------------------------------------------------------
void test_isf_single_state_interfaces()
{
    std::cerr << "\n[TEST] test_isf_single_state_interfaces\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      template A with ifaces a1,a2 and template B with b1,\n"
              << "                 each interface carrying exactly one state.\n"
              << "  Pass criteria: 3 species named Mol(iface); singleDouble==1,\n"
              << "                 implicitDouble/canDissociate==false, and one\n"
              << "                 bindPairList slot per species.\n";

    copyCounters counterArrays;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(isf_make_mol_template(
        "A", 0, { "a1", "a2" }, { { 'U' }, { 'U' } }));
    molTemplateList.push_back(isf_make_mol_template(
        "B", 1, { "b1" }, { { 'U' } }));
    std::vector<ForwardRxn> forwardRxns {};
    Parameters params;
    params.fromRestart = false;

    IsfRunResult res = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_single.csv");

    std::cerr << "  header  = \"" << res.header << "\"\n"
              << "  species = " << res.nSpecies << '\n';

    EXPECT_EQ(res.nSpecies, 3) << "Three single-state interfaces => three species";

    // Check the exact header tokens (order follows template/interface order).
    std::vector<std::string> tokens = isf_split_csv(res.header);
    ASSERT_EQ(tokens.size(), 4u) << "Header must be time column plus three species";
    EXPECT_EQ(tokens[0], "Time (s)") << "First column is the time column";
    EXPECT_EQ(tokens[1], "A(a1)") << "Single-state interfaces omit the ~state suffix";
    EXPECT_EQ(tokens[2], "A(a2)") << "Interfaces are written in declaration order";
    EXPECT_EQ(tokens[3], "B(b1)") << "Templates are written in list order";

    // Every monomer species is flagged as a 'single' species.
    ASSERT_EQ(counterArrays.singleDouble.size(), 3u);
    for (size_t i = 0; i < counterArrays.singleDouble.size(); ++i) {
        EXPECT_EQ(counterArrays.singleDouble[i], 1)
            << "Monomer species " << i << " should be flagged singleDouble==1";
        EXPECT_FALSE(counterArrays.implicitDouble[i])
            << "Monomer species " << i << " cannot be an implicit pair";
        EXPECT_FALSE(counterArrays.canDissociate[i])
            << "Monomer species " << i << " cannot dissociate";
    }
    EXPECT_EQ(counterArrays.bindPairList.size(), 3u)
        << "fromRestart==false means one bindPairList slot per species";
}

// -----------------------------------------------------------------------------
// Test 3: a multi-state interface produces one column per state.
// -----------------------------------------------------------------------------
void test_isf_multi_state_interface()
{
    std::cerr << "\n[TEST] test_isf_multi_state_interface\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      one template with a single interface holding two\n"
              << "                 states (U and P).\n"
              << "  Pass criteria: 2 species, columns \"M(s~U)\" and \"M(s~P)\".\n";

    copyCounters counterArrays;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(isf_make_mol_template(
        "M", 0, { "s" }, { { 'U', 'P' } }));
    std::vector<ForwardRxn> forwardRxns {};
    Parameters params;
    params.fromRestart = false;

    IsfRunResult res = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_multistate.csv");

    std::cerr << "  header  = \"" << res.header << "\"\n"
              << "  species = " << res.nSpecies << '\n';

    EXPECT_EQ(res.nSpecies, 2) << "Two states => two tracked species";

    std::vector<std::string> tokens = isf_split_csv(res.header);
    ASSERT_EQ(tokens.size(), 3u) << "Time column plus one column per state";
    EXPECT_EQ(tokens[1], "M(s~U)") << "State identifier must follow a tilde";
    EXPECT_EQ(tokens[2], "M(s~P)") << "States appear in stateList order";

    ASSERT_EQ(counterArrays.singleDouble.size(), 2u);
    EXPECT_EQ(counterArrays.singleDouble[0], 1) << "State species are still monomers";
    EXPECT_EQ(counterArrays.singleDouble[1], 1) << "State species are still monomers";
    EXPECT_EQ(counterArrays.bindPairList.size(), 2u)
        << "One bindPairList slot per state species";
}

// -----------------------------------------------------------------------------
// Test 4: interface with NO states contributes nothing at all.
// -----------------------------------------------------------------------------
void test_isf_interface_without_states()
{
    std::cerr << "\n[TEST] test_isf_interface_without_states\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      an interface whose stateList is empty (neither the\n"
              << "                 size==1 branch nor any loop iteration applies).\n"
              << "  Pass criteria: that interface adds no column and no counter entry.\n";

    copyCounters counterArrays;
    std::vector<MolTemplate> molTemplateList;
    // First interface has one state (counted); second has none (skipped).
    molTemplateList.push_back(isf_make_mol_template(
        "E", 0, { "good", "empty" }, { { 'U' }, {} }));
    std::vector<ForwardRxn> forwardRxns {};
    Parameters params;
    params.fromRestart = false;

    IsfRunResult res = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_nostate.csv");

    std::cerr << "  header  = \"" << res.header << "\"\n"
              << "  species = " << res.nSpecies << '\n';

    EXPECT_EQ(res.nSpecies, 1) << "Only the interface with a state is counted";
    std::vector<std::string> tokens = isf_split_csv(res.header);
    ASSERT_EQ(tokens.size(), 2u) << "Time column plus exactly one species column";
    EXPECT_EQ(tokens[1], "E(good)") << "The stateless interface must be absent";
    EXPECT_EQ(counterArrays.singleDouble.size(), 1u)
        << "Stateless interfaces must not push counter entries";
}

// -----------------------------------------------------------------------------
// Test 5: a reversible bimolecular reaction adds a 'double' species.
// -----------------------------------------------------------------------------
void test_isf_bimolecular_reversible_reaction()
{
    std::cerr << "\n[TEST] test_isf_bimolecular_reversible_reaction\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      templates A(a) and B(b) plus one reversible\n"
              << "                 bimolecular reaction with productName A(a!1).B(b!1).\n"
              << "  Pass criteria: 3 species; last entry has singleDouble==2,\n"
              << "                 canDissociate==true, implicitDouble==false.\n";

    copyCounters counterArrays;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(isf_make_mol_template("A", 0, { "a" }, { { 'U' } }));
    molTemplateList.push_back(isf_make_mol_template("B", 1, { "b" }, { { 'U' } }));

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(isf_make_forward_rxn(
        "A(a!1).B(b!1)", 0, 1, ReactionType::bimolecular, /*isReversible=*/true));

    Parameters params;
    params.fromRestart = false;

    IsfRunResult res = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_bimol.csv");

    std::cerr << "  header  = \"" << res.header << "\"\n"
              << "  species = " << res.nSpecies << '\n';

    EXPECT_EQ(res.nSpecies, 3) << "Two monomers plus one bimolecular product";

    std::vector<std::string> tokens = isf_split_csv(res.header);
    ASSERT_EQ(tokens.size(), 4u) << "Time column plus three species";
    EXPECT_EQ(tokens[3], "A(a!1).B(b!1)")
        << "The bimolecular product column uses ForwardRxn::productName";

    ASSERT_EQ(counterArrays.singleDouble.size(), 3u);
    EXPECT_EQ(counterArrays.singleDouble.back(), 2)
        << "Bimolecular product contains two species => singleDouble==2";
    ASSERT_EQ(counterArrays.canDissociate.size(), 3u);
    EXPECT_TRUE(counterArrays.canDissociate.back())
        << "Reversible reactions mark the product as dissociable";
    ASSERT_EQ(counterArrays.implicitDouble.size(), 3u);
    EXPECT_FALSE(counterArrays.implicitDouble.back())
        << "Neither reactant is an implicit lipid";
    EXPECT_EQ(counterArrays.bindPairList.size(), 3u)
        << "One bindPairList slot per species, product included";
}

// -----------------------------------------------------------------------------
// Test 6: irreversible reaction involving an implicit lipid.
// -----------------------------------------------------------------------------
void test_isf_implicit_lipid_irreversible_reaction()
{
    std::cerr << "\n[TEST] test_isf_implicit_lipid_irreversible_reaction\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      second reactant template is an implicit lipid and\n"
              << "                 the bimolecular reaction is irreversible.\n"
              << "  Pass criteria: product entry has implicitDouble==true and\n"
              << "                 canDissociate==false.\n";

    copyCounters counterArrays;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(isf_make_mol_template("P", 0, { "p" }, { { 'U' } }));
    molTemplateList.push_back(isf_make_mol_template(
        "IL", 1, { "il" }, { { 'U' } }, /*isImplicitLipid=*/true));

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(isf_make_forward_rxn(
        "P(p!1).IL(il!1)", 0, 1, ReactionType::bimolecular, /*isReversible=*/false));

    Parameters params;
    params.fromRestart = false;

    IsfRunResult res = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_implicit.csv");

    std::cerr << "  header  = \"" << res.header << "\"\n"
              << "  species = " << res.nSpecies << '\n';

    EXPECT_EQ(res.nSpecies, 3) << "Two monomers plus the implicit-lipid product";
    ASSERT_EQ(counterArrays.implicitDouble.size(), 3u);
    EXPECT_TRUE(counterArrays.implicitDouble.back())
        << "An implicit lipid reactant must set implicitDouble true";
    ASSERT_EQ(counterArrays.canDissociate.size(), 3u);
    EXPECT_FALSE(counterArrays.canDissociate.back())
        << "Irreversible reactions cannot dissociate";
    EXPECT_EQ(counterArrays.singleDouble.back(), 2)
        << "The product still counts as a two-species complex";
}

// -----------------------------------------------------------------------------
// Test 7: non-bimolecular reactions are ignored entirely.
// -----------------------------------------------------------------------------
void test_isf_non_bimolecular_reactions_ignored()
{
    std::cerr << "\n[TEST] test_isf_non_bimolecular_reactions_ignored\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      forwardRxns contains a uniMolStateChange and a\n"
              << "                 zerothOrderCreation reaction plus one bimolecular.\n"
              << "  Pass criteria: only the bimolecular reaction produces a column.\n";

    copyCounters counterArrays;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(isf_make_mol_template("A", 0, { "a" }, { { 'U' } }));
    molTemplateList.push_back(isf_make_mol_template("B", 1, { "b" }, { { 'U' } }));

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(isf_make_forward_rxn(
        "IGNORED_UNI", 0, 0, ReactionType::uniMolStateChange, true));
    forwardRxns.push_back(isf_make_forward_rxn(
        "A(a!1).B(b!1)", 0, 1, ReactionType::bimolecular, true));
    forwardRxns.push_back(isf_make_forward_rxn(
        "IGNORED_CREATE", 1, 1, ReactionType::zerothOrderCreation, false));

    Parameters params;
    params.fromRestart = false;

    IsfRunResult res = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_mixedrxn.csv");

    std::cerr << "  header  = \"" << res.header << "\"\n"
              << "  species = " << res.nSpecies << '\n';

    EXPECT_EQ(res.nSpecies, 3) << "Only the bimolecular reaction adds a species";
    EXPECT_EQ(res.header.find("IGNORED_UNI"), std::string::npos)
        << "uniMolStateChange products must not appear in the header";
    EXPECT_EQ(res.header.find("IGNORED_CREATE"), std::string::npos)
        << "zerothOrderCreation products must not appear in the header";
    EXPECT_NE(res.header.find("A(a!1).B(b!1)"), std::string::npos)
        << "The bimolecular product must appear in the header";
    EXPECT_EQ(counterArrays.singleDouble.size(), 3u)
        << "Counter vectors only grow for tracked species";
}

// -----------------------------------------------------------------------------
// Test 8: params.fromRestart == true suppresses bindPairList growth.
// -----------------------------------------------------------------------------
void test_isf_from_restart_skips_bind_pair_list()
{
    std::cerr << "\n[TEST] test_isf_from_restart_skips_bind_pair_list\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      identical inputs run with params.fromRestart==true.\n"
              << "  Pass criteria: species count/header/other counters unchanged, but\n"
              << "                 bindPairList is left untouched (restart supplies it).\n";

    // Build identical inputs twice so the two runs can be compared directly.
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(isf_make_mol_template("A", 0, { "a" }, { { 'U' } }));
    molTemplateList.push_back(isf_make_mol_template("B", 1, { "b" }, { { 'U' } }));
    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(isf_make_forward_rxn(
        "A(a!1).B(b!1)", 0, 1, ReactionType::bimolecular, true));

    // --- run with fromRestart == false (reference behaviour) ---
    copyCounters freshCounters;
    Parameters freshParams;
    freshParams.fromRestart = false;
    IsfRunResult freshRes = isf_run(freshCounters, molTemplateList, forwardRxns,
        freshParams, "isf_tmp_fresh.csv");

    // --- run with fromRestart == true ---
    copyCounters restartCounters;
    Parameters restartParams;
    restartParams.fromRestart = true;
    IsfRunResult restartRes = isf_run(restartCounters, molTemplateList, forwardRxns,
        restartParams, "isf_tmp_restart.csv");

    std::cerr << "  fresh   species = " << freshRes.nSpecies
              << ", bindPairList = " << freshCounters.bindPairList.size() << '\n'
              << "  restart species = " << restartRes.nSpecies
              << ", bindPairList = " << restartCounters.bindPairList.size() << '\n';

    // The header and species count must not depend on the restart flag.
    EXPECT_EQ(restartRes.nSpecies, freshRes.nSpecies)
        << "The restart flag must not change the number of species";
    EXPECT_EQ(restartRes.header, freshRes.header)
        << "The restart flag must not change the header text";

    // The classification vectors still grow ...
    EXPECT_EQ(restartCounters.singleDouble.size(), freshCounters.singleDouble.size())
        << "singleDouble grows regardless of the restart flag";
    EXPECT_EQ(restartCounters.implicitDouble.size(), freshCounters.implicitDouble.size())
        << "implicitDouble grows regardless of the restart flag";
    EXPECT_EQ(restartCounters.canDissociate.size(), freshCounters.canDissociate.size())
        << "canDissociate grows regardless of the restart flag";

    // ... but bindPairList must stay empty on a restart.
    EXPECT_EQ(freshCounters.bindPairList.size(), 3u)
        << "Non-restart runs allocate one bindPairList slot per species";
    EXPECT_TRUE(restartCounters.bindPairList.empty())
        << "Restart runs must not allocate bindPairList slots";
}

// -----------------------------------------------------------------------------
// Test 9: counters accumulate across successive calls (no internal reset).
// -----------------------------------------------------------------------------
void test_isf_counters_accumulate_between_calls()
{
    std::cerr << "\n[TEST] test_isf_counters_accumulate_between_calls\n"
              << "  Function:      init_speciesFile()\n"
              << "  Scenario:      the same copyCounters object is passed to two\n"
              << "                 successive calls.\n"
              << "  Pass criteria: the function appends (does not clear), so the\n"
              << "                 vectors end up twice as long while the returned\n"
              << "                 species count stays the same each call.\n";

    copyCounters counterArrays;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(isf_make_mol_template("A", 0, { "a" }, { { 'U' } }));
    std::vector<ForwardRxn> forwardRxns {};
    Parameters params;
    params.fromRestart = false;

    IsfRunResult first = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_accum1.csv");
    const size_t sizeAfterFirst = counterArrays.singleDouble.size();

    IsfRunResult second = isf_run(counterArrays, molTemplateList, forwardRxns, params,
        "isf_tmp_accum2.csv");
    const size_t sizeAfterSecond = counterArrays.singleDouble.size();

    std::cerr << "  first call  -> species=" << first.nSpecies
              << ", singleDouble=" << sizeAfterFirst << '\n'
              << "  second call -> species=" << second.nSpecies
              << ", singleDouble=" << sizeAfterSecond << '\n';

    EXPECT_EQ(first.nSpecies, 1) << "One single-state interface => one species";
    EXPECT_EQ(second.nSpecies, first.nSpecies)
        << "The return value counts only species written by this call";
    EXPECT_EQ(sizeAfterFirst, 1u) << "First call pushes one classification entry";
    EXPECT_EQ(sizeAfterSecond, 2u)
        << "The function appends to counters without clearing them";
    EXPECT_EQ(counterArrays.bindPairList.size(), 2u)
        << "bindPairList likewise accumulates across calls";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each named helper is executed in its own TEST so all of
// them run even if one reports failures (only non-fatal EXPECT_* are used).
// -----------------------------------------------------------------------------
TEST(InitSpeciesFileTest, EmptyInputs) { test_isf_empty_inputs(); }
TEST(InitSpeciesFileTest, SingleStateInterfaces) { test_isf_single_state_interfaces(); }
TEST(InitSpeciesFileTest, MultiStateInterface) { test_isf_multi_state_interface(); }
TEST(InitSpeciesFileTest, InterfaceWithoutStates) { test_isf_interface_without_states(); }
TEST(InitSpeciesFileTest, BimolecularReversibleReaction) { test_isf_bimolecular_reversible_reaction(); }
TEST(InitSpeciesFileTest, ImplicitLipidIrreversibleReaction) { test_isf_implicit_lipid_irreversible_reaction(); }
TEST(InitSpeciesFileTest, NonBimolecularReactionsIgnored) { test_isf_non_bimolecular_reactions_ignored(); }
TEST(InitSpeciesFileTest, FromRestartSkipsBindPairList) { test_isf_from_restart_skips_bind_pair_list(); }
TEST(InitSpeciesFileTest, CountersAccumulateBetweenCalls) { test_isf_counters_accumulate_between_calls(); }