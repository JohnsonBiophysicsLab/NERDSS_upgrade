/*! \file test_parse_observable.cpp
 *
 * ### Unit test for src/parser/parse_observable.cpp
 *
 * Function under test:
 * \code
 *   SpeciesTracker::Observable parse_observable(const std::string& line,
 *       const std::vector<MolTemplate>& molTemplateList,
 *       const std::vector<ForwardRxn>& forwardRxns,
 *       const std::vector<CreateDestructRxn>& createDestructRxns);
 * \endcode
 *
 * The routine takes one line from an `start observables` block of the input
 * file. The expected format is three whitespace separated tokens:
 *
 *      <observableType> <observableName> <BNGL specie string>
 *
 * e.g.   `molecule freeA a(b~U)`
 *        `complex  dimer d(e!1).f(g!1)`
 *
 * What the function does (and therefore what we verify here):
 *   1. It reads the type keyword, lower-cases it and maps it to
 *      SpeciesTracker::ObservableType (only "molecule" and "complex" are legal).
 *   2. It stores the observable's name.
 *   3. It splits the specie string on '.' and hands each piece to
 *      parse_molecule_bngl() to obtain a ParsedMol.
 *   4. For every parsed molecule it looks up the matching MolTemplate by name
 *      and back-fills:
 *          - ParsedMol::molTypeIndex   from MolTemplate::molTypeIndex
 *          - IfaceInfo::relIndex       from Interface::index
 *          - IfaceInfo::absIndex       from the matching Interface::State::index
 *      A molecule (or interface, or state) that cannot be found is simply left
 *      at whatever the BNGL parser produced (-1 by default).
 *   5. The `complex` branch loops over the forward reactions but its inner
 *      `if` body is empty, so it must have no observable effect. It does,
 *      however, unconditionally index parsedMols[0]/[1] and
 *      productListNew[0]/[1] for every *bimolecular* reaction, so those tests
 *      always supply two molecules and two products.
 *   6. Finally it copies every parsed molecule/interface into
 *      Observable::constituentList.
 *
 * NOTE: an unrecognised observable type keyword makes the function call
 * exit(1). That path is deliberately *not* exercised, because terminating the
 * process would kill the whole gtest binary.
 *
 * NOTE: the `createDestructRxns` argument is never read by the implementation;
 * we assert that supplying a non-empty list changes nothing.
 */

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Observable.hpp"
#include "classes/class_Rxns.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// parse_observable() is defined in src/parser/parse_observable.cpp but is not
// declared in any project header, so declare it here with the exact signature
// (global scope, C++ linkage) so the test links against the real definition.
// -----------------------------------------------------------------------------
SpeciesTracker::Observable parse_observable(const std::string& line,
    const std::vector<MolTemplate>& molTemplateList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns);

namespace {

/*! \brief Build a small, fully initialised MolTemplate list used by all tests.
 *
 * Molecule 0: "a" -- iface 0 "b" with two states 'U'(abs 11) and 'P'(abs 12)
 *                 -- iface 1 "c" with a single state-less state (abs 13)
 * Molecule 1: "d" -- iface 0 "e" with a single state-less state (abs 20)
 * Molecule 2: "f" -- iface 0 "g" with a single state-less state (abs 30)
 *
 * All names are lower case and all state identifiers upper case, which matches
 * the convention the BNGL parser enforces on the strings it reads.
 */
std::vector<MolTemplate> po_make_mol_template_list()
{
    std::vector<MolTemplate> molTemplateList;

    // ---- molecule "a" ----------------------------------------------------
    MolTemplate molA;
    molA.molName = "a";
    molA.molTypeIndex = 0;
    {
        // interface "b": two explicit states
        Interface ifaceB;
        ifaceB.name = "b";
        ifaceB.index = 0; // relative index within the molecule
        ifaceB.stateList.emplace_back('U', 11); // State(char iden, int index)
        ifaceB.stateList.emplace_back('P', 12);
        molA.interfaceList.push_back(ifaceB);

        // interface "c": one state with no identity ('\0')
        Interface ifaceC;
        ifaceC.name = "c";
        ifaceC.index = 1;
        ifaceC.stateList.emplace_back(13); // explicit State(int index) -> iden '\0'
        molA.interfaceList.push_back(ifaceC);
    }
    molTemplateList.push_back(molA);

    // ---- molecule "d" ----------------------------------------------------
    MolTemplate molD;
    molD.molName = "d";
    molD.molTypeIndex = 1;
    {
        Interface ifaceE;
        ifaceE.name = "e";
        ifaceE.index = 0;
        ifaceE.stateList.emplace_back(20);
        molD.interfaceList.push_back(ifaceE);
    }
    molTemplateList.push_back(molD);

    // ---- molecule "f" ----------------------------------------------------
    MolTemplate molF;
    molF.molName = "f";
    molF.molTypeIndex = 2;
    {
        Interface ifaceG;
        ifaceG.name = "g";
        ifaceG.index = 0;
        ifaceG.stateList.emplace_back(30);
        molF.interfaceList.push_back(ifaceG);
    }
    molTemplateList.push_back(molF);

    return molTemplateList;
}

/*! \brief Dump an Observable to stderr so the test log shows what was produced. */
void po_dump_observable(const SpeciesTracker::Observable& obs)
{
    std::cerr << "    -> parsed observable: name=\"" << obs.name << "\", type="
              << static_cast<int>(obs.observableType)
              << " (1=molecule, 2=complex), currNum=" << obs.currNum
              << ", constituents=" << obs.constituentList.size() << '\n';
    for (std::size_t c = 0; c < obs.constituentList.size(); ++c) {
        const auto& cons = obs.constituentList[c];
        std::cerr << "       constituent " << c << ": molTypeIndex=" << cons.molTypeIndex
                  << ", ifaces=" << cons.interfaceList.size() << '\n';
        for (std::size_t i = 0; i < cons.interfaceList.size(); ++i) {
            const auto& ifc = cons.interfaceList[i];
            std::cerr << "         iface " << i << ": relIndex=" << ifc.relIndex
                      << ", absIndex=" << ifc.absIndex << ", state='"
                      << (ifc.state == '\0' ? '0' : ifc.state) << "', isBound=" << std::boolalpha
                      << ifc.isBound << ", bondIndex=" << ifc.bondIndex << '\n';
        }
    }
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a "molecule" observable naming one interface that carries a state.
//         Pass criteria: type/name are read correctly and the interface is
//         resolved to the relative index of the template interface and the
//         absolute index of the *matching state* ('U' -> 11, not 'P' -> 12).
// -----------------------------------------------------------------------------
void po_test_molecule_with_state()
{
    std::cerr << "\n[TEST] po_test_molecule_with_state\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    line = \"molecule freeA a(b~U)\"\n"
              << "  Criteria:    observableType==molecule, name==\"freeA\", one\n"
              << "               constituent of molTypeIndex 0 whose single iface\n"
              << "               resolves to relIndex 0 / absIndex 11 (state 'U').\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<ForwardRxn> forwardRxns {}; // no reactions needed
    const std::vector<CreateDestructRxn> createDestructRxns {};

    SpeciesTracker::Observable obs
        = parse_observable("molecule freeA a(b~U)", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(obs);

    // Header information
    EXPECT_EQ(obs.name, std::string("freeA")) << "the second token is the observable name";
    EXPECT_TRUE(obs.observableType == SpeciesTracker::ObservableType::molecule)
        << "the \"molecule\" keyword must map to ObservableType::molecule";
    EXPECT_EQ(obs.currNum, 0u) << "a freshly parsed observable has a zero population";

    // One molecule in the specie string -> exactly one constituent
    ASSERT_EQ(obs.constituentList.size(), 1u) << "\"a(b~U)\" contains a single molecule";
    EXPECT_EQ(obs.constituentList[0].molTypeIndex, 0u)
        << "molecule \"a\" is molTypeIndex 0 in the template list";

    ASSERT_EQ(obs.constituentList[0].interfaceList.size(), 1u) << "one interface was declared";
    const auto& ifc = obs.constituentList[0].interfaceList[0];
    EXPECT_EQ(ifc.relIndex, 0u) << "interface \"b\" has relative index 0";
    EXPECT_EQ(ifc.absIndex, 11u) << "state 'U' of interface \"b\" has absolute index 11";
    EXPECT_EQ(ifc.state, 'U') << "the parsed state identifier should be preserved";
    EXPECT_FALSE(ifc.isBound) << "no '!' bond was given, so the interface is free";
}

// -----------------------------------------------------------------------------
// Test 2: a "molecule" observable naming a state-less interface.
//         Pass criteria: the interface whose template state has iden '\0'
//         is matched and receives that state's absolute index.
// -----------------------------------------------------------------------------
void po_test_molecule_stateless_iface()
{
    std::cerr << "\n[TEST] po_test_molecule_stateless_iface\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    line = \"molecule freeC a(c)\" (interface has no state)\n"
              << "  Criteria:    relIndex 1 (second interface of \"a\") and absIndex 13.\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<ForwardRxn> forwardRxns {};
    const std::vector<CreateDestructRxn> createDestructRxns {};

    SpeciesTracker::Observable obs
        = parse_observable("molecule freeC a(c)", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(obs);

    EXPECT_EQ(obs.name, std::string("freeC"));
    ASSERT_EQ(obs.constituentList.size(), 1u);
    EXPECT_EQ(obs.constituentList[0].molTypeIndex, 0u) << "still molecule \"a\"";

    ASSERT_EQ(obs.constituentList[0].interfaceList.size(), 1u);
    const auto& ifc = obs.constituentList[0].interfaceList[0];
    EXPECT_EQ(ifc.relIndex, 1u) << "interface \"c\" is the second interface of \"a\"";
    EXPECT_EQ(ifc.absIndex, 13u) << "the '\\0' state of interface \"c\" has absolute index 13";
    EXPECT_EQ(ifc.state, '\0') << "no state was declared, so the state identifier stays null";
    EXPECT_FALSE(ifc.isBound) << "interface is free";
}

// -----------------------------------------------------------------------------
// Test 3: two interfaces on one molecule are both resolved independently.
// -----------------------------------------------------------------------------
void po_test_molecule_two_interfaces()
{
    std::cerr << "\n[TEST] po_test_molecule_two_interfaces\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    line = \"molecule bothIfaces a(b~P,c)\"\n"
              << "  Criteria:    both interfaces appear in the constituent, resolved\n"
              << "               to (relIndex 0, absIndex 12) and (relIndex 1, absIndex 13).\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<ForwardRxn> forwardRxns {};
    const std::vector<CreateDestructRxn> createDestructRxns {};

    SpeciesTracker::Observable obs
        = parse_observable("molecule bothIfaces a(b~P,c)", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(obs);

    ASSERT_EQ(obs.constituentList.size(), 1u) << "still only one molecule";
    ASSERT_EQ(obs.constituentList[0].interfaceList.size(), 2u)
        << "two comma-separated interfaces should both be parsed";

    // first interface: b~P
    EXPECT_EQ(obs.constituentList[0].interfaceList[0].relIndex, 0u) << "\"b\" is relative index 0";
    EXPECT_EQ(obs.constituentList[0].interfaceList[0].absIndex, 12u)
        << "state 'P' of \"b\" has absolute index 12";
    EXPECT_EQ(obs.constituentList[0].interfaceList[0].state, 'P');

    // second interface: c
    EXPECT_EQ(obs.constituentList[0].interfaceList[1].relIndex, 1u) << "\"c\" is relative index 1";
    EXPECT_EQ(obs.constituentList[0].interfaceList[1].absIndex, 13u)
        << "the state-less state of \"c\" has absolute index 13";
    EXPECT_EQ(obs.constituentList[0].interfaceList[1].state, '\0');
}

// -----------------------------------------------------------------------------
// Test 4: the observable type keyword is matched case-insensitively (it is
//         lower-cased by std::transform before the map lookup).
// -----------------------------------------------------------------------------
void po_test_type_keyword_is_case_insensitive()
{
    std::cerr << "\n[TEST] po_test_type_keyword_is_case_insensitive\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    the keyword is written \"MoLeCuLe\" / \"COMPLEX\"\n"
              << "  Criteria:    both are still mapped to the correct ObservableType.\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<ForwardRxn> forwardRxns {};
    const std::vector<CreateDestructRxn> createDestructRxns {};

    SpeciesTracker::Observable obsMol
        = parse_observable("MoLeCuLe mixedCase a(b~U)", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(obsMol);
    EXPECT_TRUE(obsMol.observableType == SpeciesTracker::ObservableType::molecule)
        << "\"MoLeCuLe\" must be lower-cased and recognised";
    EXPECT_EQ(obsMol.name, std::string("mixedCase")) << "the name itself is NOT case folded";

    // the complex branch needs two molecules (it indexes parsedMols[0] and [1])
    SpeciesTracker::Observable obsCom
        = parse_observable("COMPLEX upperCase d(e!1).f(g!1)", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(obsCom);
    EXPECT_TRUE(obsCom.observableType == SpeciesTracker::ObservableType::complex)
        << "\"COMPLEX\" must be lower-cased and recognised";
}

// -----------------------------------------------------------------------------
// Test 5: a "complex" observable made of two bound molecules, with an empty
//         reaction list so the (empty-bodied) reaction loop is not entered.
// -----------------------------------------------------------------------------
void po_test_complex_two_molecules()
{
    std::cerr << "\n[TEST] po_test_complex_two_molecules\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    line = \"complex dimer d(e!1).f(g!1)\", no reactions\n"
              << "  Criteria:    the '.' splits the specie into two constituents whose\n"
              << "               molTypeIndex/relIndex/absIndex are resolved and whose\n"
              << "               interfaces are flagged as bound.\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<ForwardRxn> forwardRxns {};
    const std::vector<CreateDestructRxn> createDestructRxns {};

    SpeciesTracker::Observable obs
        = parse_observable("complex dimer d(e!1).f(g!1)", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(obs);

    EXPECT_EQ(obs.name, std::string("dimer"));
    EXPECT_TRUE(obs.observableType == SpeciesTracker::ObservableType::complex);

    ASSERT_EQ(obs.constituentList.size(), 2u) << "the specie string contains one '.' -> two molecules";

    // First constituent: d(e!1)
    EXPECT_EQ(obs.constituentList[0].molTypeIndex, 1u) << "molecule \"d\" is molTypeIndex 1";
    ASSERT_EQ(obs.constituentList[0].interfaceList.size(), 1u);
    EXPECT_EQ(obs.constituentList[0].interfaceList[0].relIndex, 0u) << "\"e\" is relative index 0";
    EXPECT_EQ(obs.constituentList[0].interfaceList[0].absIndex, 20u) << "\"e\" has absolute index 20";
    EXPECT_TRUE(obs.constituentList[0].interfaceList[0].isBound)
        << "the '!' bond marker must set isBound";

    // Second constituent: f(g!1)
    EXPECT_EQ(obs.constituentList[1].molTypeIndex, 2u) << "molecule \"f\" is molTypeIndex 2";
    ASSERT_EQ(obs.constituentList[1].interfaceList.size(), 1u);
    EXPECT_EQ(obs.constituentList[1].interfaceList[0].relIndex, 0u) << "\"g\" is relative index 0";
    EXPECT_EQ(obs.constituentList[1].interfaceList[0].absIndex, 30u) << "\"g\" has absolute index 30";
    EXPECT_TRUE(obs.constituentList[1].interfaceList[0].isBound)
        << "the '!' bond marker must set isBound";
}

// -----------------------------------------------------------------------------
// Test 6: the same "complex" observable, but now with a populated reaction
//         list. The bimolecular branch in parse_observable computes four
//         booleans and then does nothing with them, so the returned Observable
//         must be byte-for-byte the same as in test 5.
//
//         The reaction list deliberately contains
//           - one uniMolStateChange reaction (must be skipped entirely), and
//           - one bimolecular reaction with two products (the branch that is
//             actually taken, and which indexes productListNew[0] and [1]).
// -----------------------------------------------------------------------------
void po_test_complex_with_forward_reactions()
{
    std::cerr << "\n[TEST] po_test_complex_with_forward_reactions\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    same complex observable, but forwardRxns holds one\n"
              << "               uniMolStateChange (skipped) and one bimolecular rxn.\n"
              << "  Criteria:    the reaction loop has no side effect -- the result is\n"
              << "               identical to the no-reaction case.\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<CreateDestructRxn> createDestructRxns {};

    std::vector<ForwardRxn> forwardRxns;

    // A unimolecular reaction: the `if (rxnType == bimolecular)` guard must skip
    // it, so it is safe for it to carry an empty product list.
    {
        ForwardRxn uniRxn;
        uniRxn.rxnType = ReactionType::uniMolStateChange;
        forwardRxns.push_back(uniRxn);
    }

    // A bimolecular reaction with exactly the two products that make up our
    // observed complex. productListNew[0] and [1] are both dereferenced by the
    // function, so both must exist.
    {
        ForwardRxn biRxn;
        biRxn.rxnType = ReactionType::bimolecular;
        biRxn.isSymmetric = false;
        // RxnIface(name, molTypeIndex, absIfaceIndex, relIfaceIndex, requiresState, requiresInteraction)
        biRxn.productListNew.emplace_back("e", 1, 20, 0, '\0', true);
        biRxn.productListNew.emplace_back("g", 2, 30, 0, '\0', true);
        forwardRxns.push_back(biRxn);
    }

    SpeciesTracker::Observable withRxns
        = parse_observable("complex dimer d(e!1).f(g!1)", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(withRxns);

    // Reference result computed with an empty reaction list.
    const std::vector<ForwardRxn> noRxns {};
    SpeciesTracker::Observable withoutRxns
        = parse_observable("complex dimer d(e!1).f(g!1)", molTemplateList, noRxns, createDestructRxns);

    EXPECT_EQ(withRxns.name, withoutRxns.name) << "reactions must not alter the name";
    EXPECT_TRUE(withRxns.observableType == withoutRxns.observableType)
        << "reactions must not alter the observable type";
    ASSERT_EQ(withRxns.constituentList.size(), withoutRxns.constituentList.size())
        << "reactions must not add or remove constituents";

    for (std::size_t c = 0; c < withRxns.constituentList.size(); ++c) {
        EXPECT_EQ(withRxns.constituentList[c].molTypeIndex, withoutRxns.constituentList[c].molTypeIndex)
            << "constituent " << c << " molTypeIndex must be unchanged";
        ASSERT_EQ(withRxns.constituentList[c].interfaceList.size(),
            withoutRxns.constituentList[c].interfaceList.size())
            << "constituent " << c << " interface count must be unchanged";
        for (std::size_t i = 0; i < withRxns.constituentList[c].interfaceList.size(); ++i) {
            EXPECT_EQ(withRxns.constituentList[c].interfaceList[i].relIndex,
                withoutRxns.constituentList[c].interfaceList[i].relIndex)
                << "relIndex must be unchanged (constituent " << c << ", iface " << i << ')';
            EXPECT_EQ(withRxns.constituentList[c].interfaceList[i].absIndex,
                withoutRxns.constituentList[c].interfaceList[i].absIndex)
                << "absIndex must be unchanged (constituent " << c << ", iface " << i << ')';
            EXPECT_EQ(withRxns.constituentList[c].interfaceList[i].isBound,
                withoutRxns.constituentList[c].interfaceList[i].isBound)
                << "isBound must be unchanged (constituent " << c << ", iface " << i << ')';
        }
    }
}

// -----------------------------------------------------------------------------
// Test 7: the createDestructRxns argument is never read by the implementation,
//         so passing a non-empty list must not change the parsed result.
// -----------------------------------------------------------------------------
void po_test_createdestruct_list_is_ignored()
{
    std::cerr << "\n[TEST] po_test_createdestruct_list_is_ignored\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    identical line parsed with an empty and a non-empty\n"
              << "               createDestructRxns list.\n"
              << "  Criteria:    both calls yield the same Observable (the argument is\n"
              << "               unused by the implementation).\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<ForwardRxn> forwardRxns {};

    std::vector<CreateDestructRxn> cdRxns;
    cdRxns.emplace_back(); // default constructed, contents irrelevant

    const std::vector<CreateDestructRxn> emptyCdRxns {};

    SpeciesTracker::Observable a
        = parse_observable("molecule freeA a(b~U)", molTemplateList, forwardRxns, emptyCdRxns);
    SpeciesTracker::Observable b
        = parse_observable("molecule freeA a(b~U)", molTemplateList, forwardRxns, cdRxns);
    po_dump_observable(b);

    EXPECT_EQ(a.name, b.name) << "name must not depend on createDestructRxns";
    EXPECT_TRUE(a.observableType == b.observableType)
        << "observable type must not depend on createDestructRxns";
    ASSERT_EQ(a.constituentList.size(), b.constituentList.size());
    ASSERT_EQ(a.constituentList[0].interfaceList.size(), b.constituentList[0].interfaceList.size());
    EXPECT_EQ(a.constituentList[0].molTypeIndex, b.constituentList[0].molTypeIndex);
    EXPECT_EQ(a.constituentList[0].interfaceList[0].absIndex, b.constituentList[0].interfaceList[0].absIndex);
}

// -----------------------------------------------------------------------------
// Test 8: a molecule name that does not exist in the template list is left
//         unresolved. ParsedMol::molTypeIndex stays at its default of -1, which
//         is copied into Constituent::molTypeIndex -- an *unsigned* field --
//         and therefore wraps to std::numeric_limits<unsigned>::max().
//         This documents the current (unchecked) behaviour of the function.
// -----------------------------------------------------------------------------
void po_test_unknown_molecule_is_left_unresolved()
{
    std::cerr << "\n[TEST] po_test_unknown_molecule_is_left_unresolved\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    line names molecule \"zzz\" which is not in the\n"
              << "               MolTemplate list.\n"
              << "  Criteria:    the function does not abort; the constituent keeps the\n"
              << "               unresolved index (-1 stored in an unsigned field).\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<ForwardRxn> forwardRxns {};
    const std::vector<CreateDestructRxn> createDestructRxns {};

    SpeciesTracker::Observable obs
        = parse_observable("molecule ghost zzz(q)", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(obs);

    EXPECT_EQ(obs.name, std::string("ghost")) << "the header is parsed regardless of the specie";
    ASSERT_EQ(obs.constituentList.size(), 1u) << "one molecule was named, so one constituent is stored";
    EXPECT_EQ(obs.constituentList[0].molTypeIndex, static_cast<unsigned>(-1))
        << "an unmatched molecule keeps ParsedMol's default molTypeIndex of -1";
    ASSERT_EQ(obs.constituentList[0].interfaceList.size(), 1u)
        << "the interface is still copied across even though it was not resolved";
    EXPECT_EQ(obs.constituentList[0].interfaceList[0].relIndex, static_cast<unsigned>(-1))
        << "an unmatched interface keeps the default relative index of -1";
}

// -----------------------------------------------------------------------------
// Test 9: only the first three whitespace separated tokens are consumed, so any
//         trailing text on the line is ignored.
// -----------------------------------------------------------------------------
void po_test_trailing_tokens_are_ignored()
{
    std::cerr << "\n[TEST] po_test_trailing_tokens_are_ignored\n"
              << "  Source file: src/parser/parse_observable.cpp\n"
              << "  Function:    parse_observable()\n"
              << "  Scenario:    line = \"molecule freeA a(b~U) extra tokens here\"\n"
              << "  Criteria:    the extra tokens are never read; the result matches the\n"
              << "               plain three-token line.\n";

    const std::vector<MolTemplate> molTemplateList = po_make_mol_template_list();
    const std::vector<ForwardRxn> forwardRxns {};
    const std::vector<CreateDestructRxn> createDestructRxns {};

    SpeciesTracker::Observable padded = parse_observable(
        "molecule freeA a(b~U) extra tokens here", molTemplateList, forwardRxns, createDestructRxns);
    po_dump_observable(padded);

    SpeciesTracker::Observable plain
        = parse_observable("molecule freeA a(b~U)", molTemplateList, forwardRxns, createDestructRxns);

    EXPECT_EQ(padded.name, plain.name) << "trailing tokens must not affect the name";
    EXPECT_TRUE(padded.observableType == plain.observableType)
        << "trailing tokens must not affect the type";
    ASSERT_EQ(padded.constituentList.size(), plain.constituentList.size())
        << "trailing tokens must not add constituents";
    ASSERT_EQ(padded.constituentList[0].interfaceList.size(), plain.constituentList[0].interfaceList.size());
    EXPECT_EQ(padded.constituentList[0].interfaceList[0].absIndex,
        plain.constituentList[0].interfaceList[0].absIndex)
        << "the resolved absolute index must be identical";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each scenario is run in its own TEST so a failure in
// one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(ParseObservableTest, MoleculeWithState) { po_test_molecule_with_state(); }
TEST(ParseObservableTest, MoleculeStatelessInterface) { po_test_molecule_stateless_iface(); }
TEST(ParseObservableTest, MoleculeTwoInterfaces) { po_test_molecule_two_interfaces(); }
TEST(ParseObservableTest, TypeKeywordCaseInsensitive) { po_test_type_keyword_is_case_insensitive(); }
TEST(ParseObservableTest, ComplexTwoMolecules) { po_test_complex_two_molecules(); }
TEST(ParseObservableTest, ComplexWithForwardReactions) { po_test_complex_with_forward_reactions(); }
TEST(ParseObservableTest, CreateDestructListIgnored) { po_test_createdestruct_list_is_ignored(); }
TEST(ParseObservableTest, UnknownMoleculeUnresolved) { po_test_unknown_molecule_is_left_unresolved(); }
TEST(ParseObservableTest, TrailingTokensIgnored) { po_test_trailing_tokens_are_ignored(); }