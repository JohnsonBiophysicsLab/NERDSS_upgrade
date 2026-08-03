/*! \file test_init_NboundPairs.cpp
 *
 * ### Unit test for src/io/init_NboundPairs.cpp
 *
 * This test exercises the single function defined in that file:
 *
 *     void init_NboundPairs(copyCounters& counterArrays, std::ofstream& outfile,
 *                           Parameters params,
 *                           std::vector<MolTemplate>& molTemplateList,
 *                           std::vector<Molecule>& moleculeList)
 *
 * What the function is supposed to do:
 *   1. Write a single header line to `outfile` consisting of "TIME(s)", one
 *      column per unique reacting protein pair ('A,B'), and then the fixed
 *      counter columns (Nloops, nOverlapPartner, ... nAssocSuccess).
 *   2. Fill `counterArrays.proPairlist` with the linear index
 *      (p1 * params.numMolTypes + p2) of every unique pair, only storing a pair
 *      once, i.e. only when p2 >= p1.
 *   3. Size `counterArrays.nBoundPairs` to numMolTypes * numMolTypes and
 *      initialize every entry to 0.
 *   4. If params.fromRestart is true, walk `moleculeList` and increment
 *      nBoundPairs[p1 * numMolTypes + p2] for every bond, counting each bond
 *      exactly once (by requiring molecule index bounder1 >= bounder2).
 *
 * Each test below constructs a small synthetic system, calls the function, and
 * checks those four behaviours. Verbose progress is printed to stderr.
 */

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_copyCounters.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Declaration of the function under test.
//
// We declare it by hand (rather than including io/io.hpp) so that this
// translation unit does not pull in headers that #include "split.cpp", which
// would risk duplicate symbol definitions when every generated test is linked
// into one binary.
// -----------------------------------------------------------------------------
void init_NboundPairs(copyCounters& counterArrays, std::ofstream& outfile,
    Parameters params, std::vector<MolTemplate>& molTemplateList,
    std::vector<Molecule>& moleculeList);

namespace {

/*! \brief Build a minimal MolTemplate with a name, type index and partner list.
 *
 * init_NboundPairs() only reads `molName` and `rxnPartners`, so nothing else
 * has to be set up.
 */
MolTemplate inbp_make_template(const std::string& name, int typeIndex,
    const std::vector<int>& partners)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.rxnPartners = partners;
    return temp;
}

/*! \brief Build a minimal Molecule with an index, a type and a bound-partner list.
 *
 * Only `index`, `molTypeIndex` and `bndpartner` are consulted by the restart
 * branch of init_NboundPairs().
 */
Molecule inbp_make_molecule(int index, int molTypeIndex,
    const std::vector<int>& bndpartner)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = molTypeIndex;
    mol.bndpartner = bndpartner;
    return mol;
}

/*! \brief Read the entire contents of a file into a std::string. */
std::string inbp_read_file(const std::string& fileName)
{
    std::ifstream in(fileName);
    std::ostringstream oss;
    oss << in.rdbuf();
    return oss.str();
}

/*! \brief Convenience: run init_NboundPairs writing into a named temp file and
 *         return the resulting file contents.
 */
std::string inbp_run(copyCounters& counterArrays, const Parameters& params,
    std::vector<MolTemplate>& molTemplateList, std::vector<Molecule>& moleculeList,
    const std::string& fileName)
{
    std::ofstream outfile(fileName);
    init_NboundPairs(counterArrays, outfile, params, molTemplateList, moleculeList);
    outfile.close();
    const std::string contents = inbp_read_file(fileName);
    std::remove(fileName.c_str());
    return contents;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Two molecule types that react with each other.
//
// A (index 0) has partner B (index 1); B (index 1) has partner A (index 0).
// Only the pair with p2 >= p1 (i.e. 0,1) may be stored, so proPairlist must
// hold exactly one entry with value 0 * 2 + 1 == 1.
// The header line must name the pair as 'A,B' and end with the fixed columns.
// -----------------------------------------------------------------------------
void test_inbp_header_and_pairs_two_types()
{
    std::cerr << "\n[TEST] test_inbp_header_and_pairs_two_types\n"
              << "  Source file: src/io/init_NboundPairs.cpp\n"
              << "  Function:    init_NboundPairs()\n"
              << "  Scenario:    two mol types A and B that bind each other.\n"
              << "  Pass criteria: proPairlist == {1}; nBoundPairs has 4 zeros;\n"
              << "                 header contains TIME(s), 'A,B' and nAssocSuccess.\n";

    // Two molecule types, each listing the other as a reaction partner.
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(inbp_make_template("A", 0, { 1 }));
    molTemplateList.push_back(inbp_make_template("B", 1, { 0 }));

    Parameters params;
    params.numMolTypes = 2;
    params.fromRestart = false; // no restart bookkeeping in this test

    std::vector<Molecule> moleculeList; // unused when fromRestart == false
    copyCounters counterArrays;

    std::cerr << "  Calling init_NboundPairs...\n";
    const std::string header
        = inbp_run(counterArrays, params, molTemplateList, moleculeList,
            "test_inbp_two_types.tmp");

    // --- proPairlist: only the unique (0,1) pair should be recorded ----------
    EXPECT_EQ(counterArrays.proPairlist.size(), 1u)
        << "Exactly one unique protein pair should be stored";
    if (counterArrays.proPairlist.size() == 1u) {
        EXPECT_EQ(counterArrays.proPairlist[0], 1)
            << "Pair index for (p1=0, p2=1) must be p1*numMolTypes+p2 == 1";
    }

    // --- nBoundPairs: numMolTypes^2 entries, all initialized to zero ---------
    EXPECT_EQ(counterArrays.nBoundPairs.size(), 4u)
        << "nBoundPairs should be sized numMolTypes*numMolTypes == 4";
    for (size_t i = 0; i < counterArrays.nBoundPairs.size(); ++i) {
        EXPECT_EQ(counterArrays.nBoundPairs[i], 0)
            << "nBoundPairs[" << i << "] should start at zero";
    }

    // --- header line contents ------------------------------------------------
    std::cerr << "  Header written: " << header;
    EXPECT_NE(header.find("TIME(s)"), std::string::npos)
        << "Header must begin with the TIME(s) column";
    EXPECT_NE(header.find("'A,B'"), std::string::npos)
        << "Header must contain the pair column 'A,B'";
    EXPECT_NE(header.find("Nloops"), std::string::npos)
        << "Header must contain the Nloops column";
    EXPECT_NE(header.find("nAssocSuccess"), std::string::npos)
        << "Header must end with the nAssocSuccess column";
}

// -----------------------------------------------------------------------------
// Test 2: A single molecule type that binds to itself (homodimer).
//
// p1 == p2 == 0 satisfies p2 >= p1, so the self pair must be recorded with
// index 0 * 1 + 0 == 0, and the header should show 'A,A'.
// -----------------------------------------------------------------------------
void test_inbp_self_binding_single_type()
{
    std::cerr << "\n[TEST] test_inbp_self_binding_single_type\n"
              << "  Source file: src/io/init_NboundPairs.cpp\n"
              << "  Function:    init_NboundPairs()\n"
              << "  Scenario:    one mol type A whose partner is itself.\n"
              << "  Pass criteria: proPairlist == {0}; nBoundPairs has 1 zero;\n"
              << "                 header contains 'A,A'.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(inbp_make_template("A", 0, { 0 }));

    Parameters params;
    params.numMolTypes = 1;
    params.fromRestart = false;

    std::vector<Molecule> moleculeList;
    copyCounters counterArrays;

    std::cerr << "  Calling init_NboundPairs...\n";
    const std::string header
        = inbp_run(counterArrays, params, molTemplateList, moleculeList,
            "test_inbp_self_binding.tmp");

    // A self pair is a legal, unique pair and must be stored.
    EXPECT_EQ(counterArrays.proPairlist.size(), 1u)
        << "The self pair (A,A) should be stored once";
    if (counterArrays.proPairlist.size() == 1u) {
        EXPECT_EQ(counterArrays.proPairlist[0], 0)
            << "Self-pair index must be 0*1+0 == 0";
    }

    // With one molecule type the bound-pair array has a single element.
    EXPECT_EQ(counterArrays.nBoundPairs.size(), 1u)
        << "nBoundPairs should be sized 1*1 == 1";
    if (counterArrays.nBoundPairs.size() == 1u) {
        EXPECT_EQ(counterArrays.nBoundPairs[0], 0)
            << "The only bound-pair counter should start at zero";
    }

    std::cerr << "  Header written: " << header;
    EXPECT_NE(header.find("'A,A'"), std::string::npos)
        << "Header must contain the self-pair column 'A,A'";
}

// -----------------------------------------------------------------------------
// Test 3: Molecule types with no reaction partners.
//
// No pair columns are written and proPairlist stays empty, but nBoundPairs must
// still be allocated with numMolTypes^2 zeros.
// -----------------------------------------------------------------------------
void test_inbp_no_reaction_partners()
{
    std::cerr << "\n[TEST] test_inbp_no_reaction_partners\n"
              << "  Source file: src/io/init_NboundPairs.cpp\n"
              << "  Function:    init_NboundPairs()\n"
              << "  Scenario:    three mol types, none with reaction partners.\n"
              << "  Pass criteria: proPairlist empty; nBoundPairs has 9 zeros;\n"
              << "                 header has no pair column but keeps fixed columns.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(inbp_make_template("A", 0, {}));
    molTemplateList.push_back(inbp_make_template("B", 1, {}));
    molTemplateList.push_back(inbp_make_template("C", 2, {}));

    Parameters params;
    params.numMolTypes = 3;
    params.fromRestart = false;

    std::vector<Molecule> moleculeList;
    copyCounters counterArrays;

    std::cerr << "  Calling init_NboundPairs...\n";
    const std::string header
        = inbp_run(counterArrays, params, molTemplateList, moleculeList,
            "test_inbp_no_partners.tmp");

    EXPECT_TRUE(counterArrays.proPairlist.empty())
        << "No reaction partners means no pairs should be recorded";

    EXPECT_EQ(counterArrays.nBoundPairs.size(), 9u)
        << "nBoundPairs should still be sized numMolTypes*numMolTypes == 9";
    for (size_t i = 0; i < counterArrays.nBoundPairs.size(); ++i) {
        EXPECT_EQ(counterArrays.nBoundPairs[i], 0)
            << "nBoundPairs[" << i << "] should start at zero";
    }

    std::cerr << "  Header written: " << header;
    // No pair columns, but the fixed trailing columns must still be present.
    EXPECT_EQ(header.find("','"), std::string::npos)
        << "No protein-pair column should be written";
    EXPECT_NE(header.find("nOverlapPartner"), std::string::npos)
        << "Fixed counter columns must still be written";
}

// -----------------------------------------------------------------------------
// Test 4: Three molecule types with multiple partners.
//
// A(0) partners with B(1) and C(2)  -> stores indices 0*3+1 = 1 and 0*3+2 = 2
// B(1) partners with A(0) and B(1)  -> A skipped (p2 < p1), self pair 1*3+1 = 4
// C(2) partners with A(0)           -> skipped (p2 < p1)
// Expected proPairlist (in order): {1, 2, 4}
// -----------------------------------------------------------------------------
void test_inbp_multiple_partners_ordering()
{
    std::cerr << "\n[TEST] test_inbp_multiple_partners_ordering\n"
              << "  Source file: src/io/init_NboundPairs.cpp\n"
              << "  Function:    init_NboundPairs()\n"
              << "  Scenario:    A binds B and C, B binds A and itself, C binds A.\n"
              << "  Pass criteria: proPairlist == {1, 2, 4} (only p2 >= p1 kept).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(inbp_make_template("A", 0, { 1, 2 }));
    molTemplateList.push_back(inbp_make_template("B", 1, { 0, 1 }));
    molTemplateList.push_back(inbp_make_template("C", 2, { 0 }));

    Parameters params;
    params.numMolTypes = 3;
    params.fromRestart = false;

    std::vector<Molecule> moleculeList;
    copyCounters counterArrays;

    std::cerr << "  Calling init_NboundPairs...\n";
    const std::string header
        = inbp_run(counterArrays, params, molTemplateList, moleculeList,
            "test_inbp_multi_partners.tmp");

    // Exactly the three unique pairs (A,B), (A,C) and (B,B) should be stored.
    ASSERT_EQ(counterArrays.proPairlist.size(), 3u)
        << "Three unique pairs expected: (A,B), (A,C), (B,B)";
    EXPECT_EQ(counterArrays.proPairlist[0], 1) << "First pair should be (0,1) -> 1";
    EXPECT_EQ(counterArrays.proPairlist[1], 2) << "Second pair should be (0,2) -> 2";
    EXPECT_EQ(counterArrays.proPairlist[2], 4) << "Third pair should be (1,1) -> 4";

    // The header must name each of those pairs.
    std::cerr << "  Header written: " << header;
    EXPECT_NE(header.find("'A,B'"), std::string::npos) << "Header should list 'A,B'";
    EXPECT_NE(header.find("'A,C'"), std::string::npos) << "Header should list 'A,C'";
    EXPECT_NE(header.find("'B,B'"), std::string::npos) << "Header should list 'B,B'";
}

// -----------------------------------------------------------------------------
// Test 5: Restart bookkeeping for a heterodimer A-B.
//
// moleculeList holds molecule 0 (type A) bound to molecule 1 (type B). The
// function only counts a bond when bounder1 >= bounder2, so the bond is counted
// exactly once, from molecule 1's perspective. After swapping so p1 <= p2, the
// incremented slot is index 0*2 + 1 == 1.
// -----------------------------------------------------------------------------
void test_inbp_restart_counts_heterodimer()
{
    std::cerr << "\n[TEST] test_inbp_restart_counts_heterodimer\n"
              << "  Source file: src/io/init_NboundPairs.cpp\n"
              << "  Function:    init_NboundPairs() (fromRestart branch)\n"
              << "  Scenario:    molecule 0 (type A) bound to molecule 1 (type B).\n"
              << "  Pass criteria: nBoundPairs[1] == 1, all other entries zero.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(inbp_make_template("A", 0, { 1 }));
    molTemplateList.push_back(inbp_make_template("B", 1, { 0 }));

    Parameters params;
    params.numMolTypes = 2;
    params.fromRestart = true; // exercise the restart counting branch

    // A single A-B bond, recorded reciprocally on both molecules.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(inbp_make_molecule(0, 0, { 1 }));
    moleculeList.push_back(inbp_make_molecule(1, 1, { 0 }));

    copyCounters counterArrays;

    std::cerr << "  Calling init_NboundPairs...\n";
    inbp_run(counterArrays, params, molTemplateList, moleculeList,
        "test_inbp_restart_hetero.tmp");

    ASSERT_EQ(counterArrays.nBoundPairs.size(), 4u)
        << "nBoundPairs should be sized numMolTypes*numMolTypes == 4";

    std::cerr << "  nBoundPairs = (" << counterArrays.nBoundPairs[0] << ", "
              << counterArrays.nBoundPairs[1] << ", " << counterArrays.nBoundPairs[2]
              << ", " << counterArrays.nBoundPairs[3] << ")\n";

    EXPECT_EQ(counterArrays.nBoundPairs[1], 1)
        << "The single A-B bond should be counted once in slot (0,1)";
    EXPECT_EQ(counterArrays.nBoundPairs[0], 0) << "Slot (0,0) should remain zero";
    EXPECT_EQ(counterArrays.nBoundPairs[2], 0)
        << "Slot (1,0) is never used because pairs are normalized to p1 <= p2";
    EXPECT_EQ(counterArrays.nBoundPairs[3], 0) << "Slot (1,1) should remain zero";
}

// -----------------------------------------------------------------------------
// Test 6: Restart bookkeeping for a homodimer A-A plus an unbound molecule.
//
// Molecules 0 and 1 are both type A and bound to each other; molecule 2 is a
// free type A. Only molecule 1 satisfies bounder1 >= bounder2, so the (A,A)
// slot (index 0) is incremented exactly once.
// -----------------------------------------------------------------------------
void test_inbp_restart_counts_homodimer()
{
    std::cerr << "\n[TEST] test_inbp_restart_counts_homodimer\n"
              << "  Source file: src/io/init_NboundPairs.cpp\n"
              << "  Function:    init_NboundPairs() (fromRestart branch)\n"
              << "  Scenario:    molecules 0 and 1 (both type A) bound; molecule 2 free.\n"
              << "  Pass criteria: nBoundPairs[0] == 1 (counted once, not twice).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(inbp_make_template("A", 0, { 0 }));

    Parameters params;
    params.numMolTypes = 1;
    params.fromRestart = true;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(inbp_make_molecule(0, 0, { 1 })); // bound to mol 1
    moleculeList.push_back(inbp_make_molecule(1, 0, { 0 })); // bound to mol 0
    moleculeList.push_back(inbp_make_molecule(2, 0, {})); // free monomer

    copyCounters counterArrays;

    std::cerr << "  Calling init_NboundPairs...\n";
    inbp_run(counterArrays, params, molTemplateList, moleculeList,
        "test_inbp_restart_homo.tmp");

    ASSERT_EQ(counterArrays.nBoundPairs.size(), 1u)
        << "nBoundPairs should be sized 1*1 == 1 for a single mol type";

    std::cerr << "  nBoundPairs[0] = " << counterArrays.nBoundPairs[0] << "\n";
    EXPECT_EQ(counterArrays.nBoundPairs[0], 1)
        << "A reciprocal A-A bond must be counted exactly once, not twice";
}

// -----------------------------------------------------------------------------
// Test 7: Restart bookkeeping with several bonds spread over two types.
//
// System: mol0(A)-mol1(A) homodimer, mol2(B)-mol3(A) heterodimer.
//   - mol1 sees partner 0 (1 >= 0)  -> (A,A) slot 0*2+0 = 0
//   - mol3 sees partner 2 (3 >= 2)  -> types A and B normalized to (0,1) -> 1
// Every other visit is skipped by the bounder1 >= bounder2 filter.
// -----------------------------------------------------------------------------
void test_inbp_restart_multiple_bonds()
{
    std::cerr << "\n[TEST] test_inbp_restart_multiple_bonds\n"
              << "  Source file: src/io/init_NboundPairs.cpp\n"
              << "  Function:    init_NboundPairs() (fromRestart branch)\n"
              << "  Scenario:    one A-A bond and one A-B bond in the same system.\n"
              << "  Pass criteria: nBoundPairs == {1, 1, 0, 0}.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(inbp_make_template("A", 0, { 0, 1 }));
    molTemplateList.push_back(inbp_make_template("B", 1, { 0 }));

    Parameters params;
    params.numMolTypes = 2;
    params.fromRestart = true;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(inbp_make_molecule(0, 0, { 1 })); // A bound to mol 1
    moleculeList.push_back(inbp_make_molecule(1, 0, { 0 })); // A bound to mol 0
    moleculeList.push_back(inbp_make_molecule(2, 1, { 3 })); // B bound to mol 3
    moleculeList.push_back(inbp_make_molecule(3, 0, { 2 })); // A bound to mol 2

    copyCounters counterArrays;

    std::cerr << "  Calling init_NboundPairs...\n";
    inbp_run(counterArrays, params, molTemplateList, moleculeList,
        "test_inbp_restart_multi.tmp");

    ASSERT_EQ(counterArrays.nBoundPairs.size(), 4u)
        << "nBoundPairs should be sized numMolTypes*numMolTypes == 4";

    std::cerr << "  nBoundPairs = (" << counterArrays.nBoundPairs[0] << ", "
              << counterArrays.nBoundPairs[1] << ", " << counterArrays.nBoundPairs[2]
              << ", " << counterArrays.nBoundPairs[3] << ")\n";

    EXPECT_EQ(counterArrays.nBoundPairs[0], 1) << "One A-A bond expected in slot (0,0)";
    EXPECT_EQ(counterArrays.nBoundPairs[1], 1) << "One A-B bond expected in slot (0,1)";
    EXPECT_EQ(counterArrays.nBoundPairs[2], 0) << "Slot (1,0) must stay unused";
    EXPECT_EQ(counterArrays.nBoundPairs[3], 0) << "No B-B bonds exist";
}

// -----------------------------------------------------------------------------
// Test 8: The function appends to (rather than replaces) an existing
//         proPairlist, which is the documented push_back behaviour.
// -----------------------------------------------------------------------------
void test_inbp_appends_to_existing_pairlist()
{
    std::cerr << "\n[TEST] test_inbp_appends_to_existing_pairlist\n"
              << "  Source file: src/io/init_NboundPairs.cpp\n"
              << "  Function:    init_NboundPairs()\n"
              << "  Scenario:    proPairlist already holds a sentinel value.\n"
              << "  Pass criteria: the sentinel survives and the new pair index\n"
              << "                 is appended after it (push_back semantics).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(inbp_make_template("A", 0, { 1 }));
    molTemplateList.push_back(inbp_make_template("B", 1, { 0 }));

    Parameters params;
    params.numMolTypes = 2;
    params.fromRestart = false;

    std::vector<Molecule> moleculeList;

    // Pre-load a sentinel so we can detect append vs. overwrite behaviour.
    copyCounters counterArrays;
    counterArrays.proPairlist.push_back(-99);

    std::cerr << "  Calling init_NboundPairs with a pre-seeded proPairlist...\n";
    inbp_run(counterArrays, params, molTemplateList, moleculeList,
        "test_inbp_append.tmp");

    ASSERT_EQ(counterArrays.proPairlist.size(), 2u)
        << "The pair index should be appended, leaving the sentinel in place";
    EXPECT_EQ(counterArrays.proPairlist[0], -99)
        << "Pre-existing sentinel entry must not be modified";
    EXPECT_EQ(counterArrays.proPairlist[1], 1)
        << "Newly appended pair index for (0,1) must be 1";

    std::cerr << "  proPairlist = (" << counterArrays.proPairlist[0] << ", "
              << counterArrays.proPairlist[1] << ")\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is invoked from its own TEST so that a
// failure in one scenario does not prevent the remaining scenarios from running
// (all assertions used are non-fatal except the size preconditions, which are
// local to a single test).
// -----------------------------------------------------------------------------
TEST(InitNboundPairs, HeaderAndPairsTwoTypes) { test_inbp_header_and_pairs_two_types(); }
TEST(InitNboundPairs, SelfBindingSingleType) { test_inbp_self_binding_single_type(); }
TEST(InitNboundPairs, NoReactionPartners) { test_inbp_no_reaction_partners(); }
TEST(InitNboundPairs, MultiplePartnersOrdering) { test_inbp_multiple_partners_ordering(); }
TEST(InitNboundPairs, RestartCountsHeterodimer) { test_inbp_restart_counts_heterodimer(); }
TEST(InitNboundPairs, RestartCountsHomodimer) { test_inbp_restart_counts_homodimer(); }
TEST(InitNboundPairs, RestartMultipleBonds) { test_inbp_restart_multiple_bonds(); }
TEST(InitNboundPairs, AppendsToExistingPairlist) { test_inbp_appends_to_existing_pairlist(); }