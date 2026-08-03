/*! \file test_determine_if_reaction_occurs.cpp
 *
 * ### Unit test for src/reactions/determine_if_reaction_occurs.cpp
 *
 * Function under test:
 *
 *     bool determine_if_reaction_occurs(int& crossIndex1, int& crossIndex2,
 *                                       const double maxRandInt, Molecule& mol,
 *                                       std::vector<Molecule>& moleculeList,
 *                                       const std::vector<ForwardRxn>& forwardRxns)
 *
 * The routine walks the list of potential reaction partners stored on the
 * Molecule (`crossbase`, `probvec`, `crossrxn`, `mycrossint`).  For each
 * candidate it draws a uniform random number in [0,1) and, if that number is
 * smaller than the stored reaction probability, it tries to locate the
 * matching entry on the *partner* molecule's cross lists.  It returns true
 * (and fills in crossIndex1 / crossIndex2) only when a fully consistent pair
 * of cross entries is found, or immediately when the partner is an implicit
 * lipid.
 *
 * Because the acceptance test is `rand < probvec`, the stochastic behaviour
 * can be made deterministic for testing purposes:
 *   - probvec entry of 2.0  => always accepted (rand_gsl64() < 1.0 always)
 *   - probvec entry of -1.0 => never accepted
 *
 * Every test below therefore has a completely deterministic pass/fail
 * criterion even though the function uses the GSL random number generator.
 */

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (unique dirocc_ prefix so they cannot clash with other tests)
// -----------------------------------------------------------------------------
namespace {

//! Probability value guaranteed to be accepted by `rand < prob` (rand < 1).
constexpr double kDiroccAlwaysReact = 2.0;
//! Probability value guaranteed to be rejected by `rand < prob` (rand >= 0).
constexpr double kDiroccNeverReact = -1.0;

/*! \brief Make sure the global GSL random number generator is usable.
 *
 * gtest_main.cpp defines `gsl_rng* r = nullptr;`.  The function under test
 * calls rand_gsl64(), which dereferences `r`, so it must be seeded once
 * before any test runs.
 */
void dirocc_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] global gsl_rng* r was null -> calling srand_gsl(1)\n";
        srand_gsl(1);
    }
    // If seeding somehow failed the tests below would crash, so report it.
    EXPECT_NE(r, nullptr) << "GSL random number generator must be initialized";
}

/*! \brief Build a minimal Molecule with a single interface.
 *
 * \param[in] index         Index of this molecule inside moleculeList.
 * \param[in] absIfaceIndex Absolute interface index stored on interface 0.
 *                          This is what the function compares against
 *                          ForwardRxn::reactantListNew[*].absIfaceIndex.
 */
Molecule dirocc_make_molecule(int index, int absIfaceIndex)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = index;
    mol.isImplicitLipid = false;

    Molecule::Iface iface;
    iface.index = absIfaceIndex; // absolute interface (state) index
    iface.relIndex = 0;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // Empty encounter lists by default; each test fills them in.
    mol.crossbase.clear();
    mol.mycrossint.clear();
    mol.crossrxn.clear();
    mol.probvec.clear();
    return mol;
}

/*! \brief Append one "encounter" entry to a molecule's cross lists.
 *
 * \param[in,out] mol           Molecule receiving the encounter.
 * \param[in] partnerMolIndex   moleculeList index of the encountered partner.
 * \param[in] relIfaceIndex     Which of this molecule's interfaces encountered.
 * \param[in] rxnIndex          Index into forwardRxns for this encounter.
 * \param[in] prob              Reaction probability for this encounter.
 */
void dirocc_add_encounter(Molecule& mol, int partnerMolIndex, int relIfaceIndex, int rxnIndex, double prob)
{
    mol.crossbase.push_back(partnerMolIndex);
    mol.mycrossint.push_back(relIfaceIndex);
    mol.crossrxn.push_back(std::array<int, 3> { rxnIndex, 0, 0 }); // [rxnIndex, rateIndex, isBackRxn]
    mol.probvec.push_back(prob);
}

/*! \brief Build a single bimolecular ForwardRxn with the two given reactants.
 *
 * Only `reactantListNew[0/1].absIfaceIndex` is read by the function under
 * test, so nothing else needs to be filled in.
 */
std::vector<ForwardRxn> dirocc_make_forward_rxns(int absIface1, int absIface2)
{
    ForwardRxn rxn;
    RxnIface react1;
    react1.absIfaceIndex = absIface1;
    RxnIface react2;
    react2.absIfaceIndex = absIface2;
    rxn.reactantListNew.clear();
    rxn.reactantListNew.push_back(react1);
    rxn.reactantListNew.push_back(react2);

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(rxn);
    return forwardRxns;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: molecule has no encounters at all -> must return false and must not
//         touch the caller's crossIndex outputs.
// -----------------------------------------------------------------------------
void test_dirocc_no_encounters_returns_false()
{
    std::cerr << "\n[TEST] test_dirocc_no_encounters_returns_false\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      mol.crossbase is empty (no potential partners).\n"
              << "  Pass criteria: returns false, crossIndex1/2 left untouched.\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10) };

    // Sentinel values so we can prove they were not overwritten.
    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_FALSE(willReact) << "With no encounters no reaction can be chosen";
    EXPECT_EQ(crossIndex1, -99) << "crossIndex1 must not be written when no reaction occurs";
    EXPECT_EQ(crossIndex2, -99) << "crossIndex2 must not be written when no reaction occurs";
}

// -----------------------------------------------------------------------------
// Test 2: encounter exists but its probability can never be accepted.
// -----------------------------------------------------------------------------
void test_dirocc_zero_probability_returns_false()
{
    std::cerr << "\n[TEST] test_dirocc_zero_probability_returns_false\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      single encounter with probability -1 (rand is\n"
              << "                 always >= 0, so acceptance is impossible).\n"
              << "  Pass criteria: returns false, crossIndex1/2 left untouched.\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 20) };

    // molecule 0 sees molecule 1, but with an impossible probability.
    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccNeverReact);
    dirocc_add_encounter(moleculeList[1], 0, 0, 0, kDiroccNeverReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_FALSE(willReact) << "A negative probability can never be accepted";
    EXPECT_EQ(crossIndex1, -99) << "crossIndex1 must remain unmodified";
    EXPECT_EQ(crossIndex2, -99) << "crossIndex2 must remain unmodified";
}

// -----------------------------------------------------------------------------
// Test 3: partner is an implicit lipid -> shortcut return of true with
//         crossIndex2 forced to 0.
// -----------------------------------------------------------------------------
void test_dirocc_implicit_lipid_partner()
{
    std::cerr << "\n[TEST] test_dirocc_implicit_lipid_partner\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      accepted encounter whose partner molecule has\n"
              << "                 isImplicitLipid == true.\n"
              << "  Pass criteria: returns true immediately, crossIndex1 == 0 (the\n"
              << "                 accepted encounter) and crossIndex2 == 0.\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 20) };

    // Mark molecule 1 as the implicit lipid; it needs no cross lists at all.
    moleculeList[1].isImplicitLipid = true;

    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccAlwaysReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_TRUE(willReact) << "Implicit lipid partners always accept the drawn reaction";
    EXPECT_EQ(crossIndex1, 0) << "crossIndex1 should point at the accepted encounter";
    EXPECT_EQ(crossIndex2, 0) << "crossIndex2 is hard-coded to 0 for implicit lipids";
}

// -----------------------------------------------------------------------------
// Test 4: fully consistent explicit pair, with `mol` acting as reactant 1.
// -----------------------------------------------------------------------------
void test_dirocc_matched_pair_mol_is_reactant1()
{
    std::cerr << "\n[TEST] test_dirocc_matched_pair_mol_is_reactant1\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      mol(iface 10) and partner(iface 20) point at each\n"
              << "                 other, share probability and crossrxn, and the\n"
              << "                 interfaces match reactantListNew[0]/[1].\n"
              << "  Pass criteria: returns true with crossIndex1 == 0, crossIndex2 == 0.\n";

    dirocc_ensure_rng();

    // Reaction: absolute iface 10 (reactant 1) + absolute iface 20 (reactant 2).
    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 20) };

    // Symmetric, fully consistent encounter entries.
    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccAlwaysReact);
    dirocc_add_encounter(moleculeList[1], 0, 0, 0, kDiroccAlwaysReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_TRUE(willReact) << "A fully consistent pair must be accepted";
    EXPECT_EQ(crossIndex1, 0) << "crossIndex1 must index mol's accepted encounter";
    EXPECT_EQ(crossIndex2, 0) << "crossIndex2 must index the partner's matching encounter";
}

// -----------------------------------------------------------------------------
// Test 5: the same consistent pair, but with `mol` playing the role of
//         reactant 2 (its interface matches reactantListNew[1]).
// -----------------------------------------------------------------------------
void test_dirocc_matched_pair_mol_is_reactant2()
{
    std::cerr << "\n[TEST] test_dirocc_matched_pair_mol_is_reactant2\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      mol carries iface 20 (reactant 2) and the partner\n"
              << "                 carries iface 10 (reactant 1) - the else-branch of\n"
              << "                 the reactant identification.\n"
              << "  Pass criteria: returns true with crossIndex1 == 0, crossIndex2 == 0.\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    // Note the swapped interface assignment relative to test 4.
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 20), dirocc_make_molecule(1, 10) };

    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccAlwaysReact);
    dirocc_add_encounter(moleculeList[1], 0, 0, 0, kDiroccAlwaysReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_TRUE(willReact) << "Reactant order must not matter for a consistent pair";
    EXPECT_EQ(crossIndex1, 0) << "crossIndex1 must index mol's accepted encounter";
    EXPECT_EQ(crossIndex2, 0) << "crossIndex2 must index the partner's matching encounter";
}

// -----------------------------------------------------------------------------
// Test 6: the partner's interface does not belong to the reaction -> no match.
// -----------------------------------------------------------------------------
void test_dirocc_interface_mismatch_returns_false()
{
    std::cerr << "\n[TEST] test_dirocc_interface_mismatch_returns_false\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      probabilities and crossrxn agree and the partner\n"
              << "                 points back at mol, but the partner's interface\n"
              << "                 (99) is not a reactant of the reaction.\n"
              << "  Pass criteria: returns false (reaction cannot be identified).\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    // Partner carries interface 99, which is not in reactantListNew.
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 99) };

    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccAlwaysReact);
    dirocc_add_encounter(moleculeList[1], 0, 0, 0, kDiroccAlwaysReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_FALSE(willReact) << "A non-reactant partner interface must not be accepted";
}

// -----------------------------------------------------------------------------
// Test 7: the partner does not list `mol` as its encountered molecule.
// -----------------------------------------------------------------------------
void test_dirocc_partner_does_not_point_back()
{
    std::cerr << "\n[TEST] test_dirocc_partner_does_not_point_back\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      partner's crossbase entry references a different\n"
              << "                 molecule (7) instead of mol (0).\n"
              << "  Pass criteria: returns false.\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 20) };

    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccAlwaysReact);
    // Partner's single encounter is with molecule 7, not with molecule 0.
    dirocc_add_encounter(moleculeList[1], 7, 0, 0, kDiroccAlwaysReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_FALSE(willReact) << "Both molecules must reference each other for a match";
}

// -----------------------------------------------------------------------------
// Test 8: probabilities on the two molecules disagree by more than 1e-10.
// -----------------------------------------------------------------------------
void test_dirocc_probability_mismatch_returns_false()
{
    std::cerr << "\n[TEST] test_dirocc_probability_mismatch_returns_false\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      mol stores probability 2.0 for the encounter while\n"
              << "                 the partner stores 1.0 (difference > 1e-10).\n"
              << "  Pass criteria: returns false (probabilities must agree).\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 20) };

    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccAlwaysReact); // 2.0
    dirocc_add_encounter(moleculeList[1], 0, 0, 0, 1.0); // deliberately different

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_FALSE(willReact) << "Mismatched stored probabilities must not be paired";
}

// -----------------------------------------------------------------------------
// Test 9: the crossrxn triples differ between the two molecules.
// -----------------------------------------------------------------------------
void test_dirocc_crossrxn_mismatch_returns_false()
{
    std::cerr << "\n[TEST] test_dirocc_crossrxn_mismatch_returns_false\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      identical probabilities, but the partner's crossrxn\n"
              << "                 triple names a different reaction index.\n"
              << "  Pass criteria: returns false (crossrxn triples must be equal).\n";

    dirocc_ensure_rng();

    // Two reactions so that a second index is meaningful.
    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    forwardRxns.push_back(forwardRxns[0]); // reaction index 1, same reactants

    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 20) };

    dirocc_add_encounter(moleculeList[0], 1, 0, /*rxnIndex=*/0, kDiroccAlwaysReact);
    dirocc_add_encounter(moleculeList[1], 0, 0, /*rxnIndex=*/1, kDiroccAlwaysReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_FALSE(willReact) << "Different crossrxn triples cannot describe the same event";
}

// -----------------------------------------------------------------------------
// Test 10: several encounters on `mol`; only the second one can be accepted, so
//          crossIndex1 must be 1 rather than 0.
// -----------------------------------------------------------------------------
void test_dirocc_selects_second_encounter()
{
    std::cerr << "\n[TEST] test_dirocc_selects_second_encounter\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      mol has two encounters; the first has probability\n"
              << "                 -1 (never accepted) and the second 2.0 (always).\n"
              << "  Pass criteria: returns true with crossIndex1 == 1, crossIndex2 == 0.\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    std::vector<Molecule> moleculeList {
        dirocc_make_molecule(0, 10), // mol under test, reactant 1
        dirocc_make_molecule(1, 20), // rejected candidate
        dirocc_make_molecule(2, 20) // accepted candidate
    };

    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccNeverReact); // encounter 0
    dirocc_add_encounter(moleculeList[0], 2, 0, 0, kDiroccAlwaysReact); // encounter 1

    // Both candidates mirror the encounter, but only #2 will be reached.
    dirocc_add_encounter(moleculeList[1], 0, 0, 0, kDiroccNeverReact);
    dirocc_add_encounter(moleculeList[2], 0, 0, 0, kDiroccAlwaysReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_TRUE(willReact) << "The second, always-accepted encounter must fire";
    EXPECT_EQ(crossIndex1, 1) << "crossIndex1 must be the index of the accepted encounter";
    EXPECT_EQ(crossIndex2, 0) << "crossIndex2 must index the partner's only encounter";
}

// -----------------------------------------------------------------------------
// Test 11: the matching entry is not the first one on the partner molecule, so
//          the inner search loop must advance crossIndex2.
// -----------------------------------------------------------------------------
void test_dirocc_partner_second_entry_matches()
{
    std::cerr << "\n[TEST] test_dirocc_partner_second_entry_matches\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      the partner's first encounter is with an unrelated\n"
              << "                 molecule; its second encounter is the mirror of\n"
              << "                 mol's accepted encounter.\n"
              << "  Pass criteria: returns true with crossIndex1 == 0, crossIndex2 == 1.\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };
    std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 20) };

    dirocc_add_encounter(moleculeList[0], 1, 0, 0, kDiroccAlwaysReact);

    // Partner: entry 0 references molecule 7 (irrelevant), entry 1 references mol.
    dirocc_add_encounter(moleculeList[1], 7, 0, 0, kDiroccAlwaysReact);
    dirocc_add_encounter(moleculeList[1], 0, 0, 0, kDiroccAlwaysReact);

    int crossIndex1 { -99 };
    int crossIndex2 { -99 };

    const bool willReact = determine_if_reaction_occurs(
        crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

    std::cerr << "  result: willReact=" << std::boolalpha << willReact
              << " crossIndex1=" << crossIndex1 << " crossIndex2=" << crossIndex2 << '\n';

    EXPECT_TRUE(willReact) << "The partner's later matching entry must still be found";
    EXPECT_EQ(crossIndex1, 0) << "crossIndex1 must index mol's accepted encounter";
    EXPECT_EQ(crossIndex2, 1) << "crossIndex2 must advance to the matching partner entry";
}

// -----------------------------------------------------------------------------
// Test 12: repeated calls with a genuinely random probability must always
//          produce a self-consistent answer (true only for a valid pair).
// -----------------------------------------------------------------------------
void test_dirocc_stochastic_consistency()
{
    std::cerr << "\n[TEST] test_dirocc_stochastic_consistency\n"
              << "  Source file:   determine_if_reaction_occurs.cpp\n"
              << "  Function:      determine_if_reaction_occurs\n"
              << "  Scenario:      a valid pair with probability 0.5 is queried 200\n"
              << "                 times using the real GSL random stream.\n"
              << "  Pass criteria: whenever true is returned the indices point at the\n"
              << "                 valid pair; at least one acceptance and one\n"
              << "                 rejection are observed over 200 draws.\n";

    dirocc_ensure_rng();

    std::vector<ForwardRxn> forwardRxns { dirocc_make_forward_rxns(10, 20) };

    int numTrue { 0 };
    int numFalse { 0 };

    for (int trial = 0; trial < 200; ++trial) {
        // Rebuild fresh molecules each trial: the function may modify the lists'
        // consumers elsewhere, so keep every trial independent.
        std::vector<Molecule> moleculeList { dirocc_make_molecule(0, 10), dirocc_make_molecule(1, 20) };
        dirocc_add_encounter(moleculeList[0], 1, 0, 0, 0.5);
        dirocc_add_encounter(moleculeList[1], 0, 0, 0, 0.5);

        int crossIndex1 { -99 };
        int crossIndex2 { -99 };
        const bool willReact = determine_if_reaction_occurs(
            crossIndex1, crossIndex2, 0.0, moleculeList[0], moleculeList, forwardRxns);

        if (willReact) {
            ++numTrue;
            // When accepted, the indices must describe the one valid pairing.
            EXPECT_EQ(crossIndex1, 0) << "accepted reaction must use mol encounter 0";
            EXPECT_EQ(crossIndex2, 0) << "accepted reaction must use partner encounter 0";
        } else {
            ++numFalse;
        }
    }

    std::cerr << "  observed " << numTrue << " acceptances and " << numFalse
              << " rejections out of 200 draws (expected ~50/50)\n";

    // With p = 0.5 and 200 draws, seeing zero of either outcome is
    // astronomically unlikely, so this is a safe determinism-free check.
    EXPECT_GT(numTrue, 0) << "Some draws should accept the p=0.5 reaction";
    EXPECT_GT(numFalse, 0) << "Some draws should reject the p=0.5 reaction";
    EXPECT_EQ(numTrue + numFalse, 200) << "Every trial must return a definite answer";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: one TEST per named helper so all of them execute even if
// an earlier expectation fails (only non-fatal EXPECT_* assertions are used).
// -----------------------------------------------------------------------------
TEST(DetermineIfReactionOccurs, NoEncountersReturnsFalse) { test_dirocc_no_encounters_returns_false(); }
TEST(DetermineIfReactionOccurs, ZeroProbabilityReturnsFalse) { test_dirocc_zero_probability_returns_false(); }
TEST(DetermineIfReactionOccurs, ImplicitLipidPartner) { test_dirocc_implicit_lipid_partner(); }
TEST(DetermineIfReactionOccurs, MatchedPairMolIsReactant1) { test_dirocc_matched_pair_mol_is_reactant1(); }
TEST(DetermineIfReactionOccurs, MatchedPairMolIsReactant2) { test_dirocc_matched_pair_mol_is_reactant2(); }
TEST(DetermineIfReactionOccurs, InterfaceMismatchReturnsFalse) { test_dirocc_interface_mismatch_returns_false(); }
TEST(DetermineIfReactionOccurs, PartnerDoesNotPointBack) { test_dirocc_partner_does_not_point_back(); }
TEST(DetermineIfReactionOccurs, ProbabilityMismatchReturnsFalse) { test_dirocc_probability_mismatch_returns_false(); }
TEST(DetermineIfReactionOccurs, CrossRxnMismatchReturnsFalse) { test_dirocc_crossrxn_mismatch_returns_false(); }
TEST(DetermineIfReactionOccurs, SelectsSecondEncounter) { test_dirocc_selects_second_encounter(); }
TEST(DetermineIfReactionOccurs, PartnerSecondEntryMatches) { test_dirocc_partner_second_entry_matches(); }
TEST(DetermineIfReactionOccurs, StochasticConsistency) { test_dirocc_stochastic_consistency(); }