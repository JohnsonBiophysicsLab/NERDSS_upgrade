/*! \file test_hasIntangibles.cpp
 *
 * ### Unit test for src/reactions/hasIntangibles.cpp
 *
 * The file under test provides two overloads:
 *
 *   bool hasIntangibles(int reactantIndex, const Molecule& reactMol,
 *                       const RxnBase::RateState& currRxnState);
 *
 *   bool hasIntangibles(int reactIndex1, int reactIndex2,
 *                       const Molecule& reactMol1, const Molecule& reactMol2,
 *                       const RxnBase::RateState& currRxnState);
 *
 * Both functions answer the question: "does the reacting Molecule (or pair of
 * Molecules) carry every ancillary ('intangible') interface that this reaction
 * rate state requires?"
 *
 * An ancillary interface requirement (an RxnIface) is considered satisfied by a
 * Molecule::Iface only when ALL FOUR of the following match:
 *   - molTypeIndex        == molTypeIndex
 *   - relIfaceIndex       == relIndex
 *   - requiresInteraction == isBound
 *   - requiresState       == stateIden
 *
 * The tests below therefore build small synthetic Molecules and RateStates and
 * check each of those four criteria individually, plus the special-case early
 * return in the two-reactant overload when both ancillary lists are empty.
 *
 * Verbose progress messages are written to stderr so a reader can follow which
 * function is being exercised and what criteria decide pass/fail.
 */

#include "reactions/shared_reaction_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (prefixed hi_ so they cannot collide with other test files).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a single Molecule::Iface with the four fields used for matching.
 *
 * \param[in] molTypeIndex Index of the parent Molecule's MolTemplate.
 * \param[in] relIndex     Relative index of the interface on the Molecule.
 * \param[in] stateIden    State character of the interface ('\0' == no state).
 * \param[in] isBound      Whether the interface is currently bound.
 */
Molecule::Iface hi_make_iface(int molTypeIndex, int relIndex, char stateIden, bool isBound)
{
    Molecule::Iface iface;
    iface.molTypeIndex = molTypeIndex;
    iface.relIndex = relIndex;
    iface.stateIden = stateIden;
    iface.isBound = isBound;
    return iface;
}

/*! \brief Build an ancillary reaction-interface requirement (RxnIface).
 *
 * The name string is irrelevant to the matching logic but is required by the
 * constructor, so a descriptive placeholder is used.
 */
RxnIface hi_make_rxn_iface(int molTypeIndex, int relIfaceIndex, char requiresState, bool requiresInteraction)
{
    // absIfaceIndex is not used by hasIntangibles(); pass a dummy value.
    return RxnIface(std::string("ancillary"), molTypeIndex, /*absIfaceIndex=*/0, relIfaceIndex, requiresState,
        requiresInteraction);
}

/*! \brief Build a Molecule holding the provided list of interfaces. */
Molecule hi_make_molecule(int molTypeIndex, const std::vector<Molecule::Iface>& ifaces)
{
    Molecule mol;
    mol.molTypeIndex = molTypeIndex;
    mol.interfaceList = ifaces;
    return mol;
}

/*! \brief Build a RateState whose otherIfaceLists holds exactly two lists.
 *
 * Both overloads of hasIntangibles index into otherIfaceLists, and the
 * two-reactant overload always touches indices 0 and 1, so every RateState used
 * in these tests carries two (possibly empty) lists.
 */
RxnBase::RateState hi_make_rate_state(
    const std::vector<RxnIface>& list0, const std::vector<RxnIface>& list1)
{
    RxnBase::RateState rateState;
    rateState.rate = 1.0;
    rateState.otherIfaceLists.clear();
    rateState.otherIfaceLists.push_back(list0);
    rateState.otherIfaceLists.push_back(list1);
    return rateState;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: single-reactant overload, empty ancillary list.
// -----------------------------------------------------------------------------
void test_hi_unimol_empty_requirement_returns_true()
{
    std::cerr << "\n[TEST] test_hi_unimol_empty_requirement_returns_true\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, const Molecule&, const RateState&)\n"
              << "  Scenario:    the rate state requires no ancillary interfaces.\n"
              << "  Pass:        function returns true (nothing to check).\n";

    // A molecule with one arbitrary interface.
    Molecule mol = hi_make_molecule(0, { hi_make_iface(0, 0, 'U', false) });

    // Reactant 0 has an empty ancillary requirement list.
    RxnBase::RateState rateState = hi_make_rate_state({}, {});

    const bool result = hasIntangibles(0, mol, rateState);
    std::cerr << "  hasIntangibles returned " << std::boolalpha << result << '\n';
    EXPECT_TRUE(result) << "With no ancillary requirements the molecule trivially qualifies";
}

// -----------------------------------------------------------------------------
// Test 2: single-reactant overload, all requirements satisfied.
// -----------------------------------------------------------------------------
void test_hi_unimol_all_requirements_present_returns_true()
{
    std::cerr << "\n[TEST] test_hi_unimol_all_requirements_present_returns_true\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, const Molecule&, const RateState&)\n"
              << "  Scenario:    molecule carries every required ancillary interface\n"
              << "               (two requirements, both matched exactly).\n"
              << "  Pass:        function returns true.\n";

    // Molecule of type 3 with three interfaces.
    Molecule mol = hi_make_molecule(3,
        { hi_make_iface(3, 0, 'U', false),   // unbound, state U
            hi_make_iface(3, 1, 'P', true),  // bound, state P
            hi_make_iface(3, 2, '\0', false) });

    // Require iface 1 bound in state P, and iface 2 unbound with no state.
    RxnBase::RateState rateState = hi_make_rate_state(
        { hi_make_rxn_iface(3, 1, 'P', true), hi_make_rxn_iface(3, 2, '\0', false) }, {});

    const bool result = hasIntangibles(0, mol, rateState);
    std::cerr << "  hasIntangibles returned " << std::boolalpha << result << '\n';
    EXPECT_TRUE(result) << "Both ancillary requirements are present on the molecule";
}

// -----------------------------------------------------------------------------
// Test 3: single-reactant overload, each of the four match criteria broken.
// -----------------------------------------------------------------------------
void test_hi_unimol_mismatch_each_field_returns_false()
{
    std::cerr << "\n[TEST] test_hi_unimol_mismatch_each_field_returns_false\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, const Molecule&, const RateState&)\n"
              << "  Scenario:    break one matching field at a time (molTypeIndex,\n"
              << "               relIndex, state identity, bound status).\n"
              << "  Pass:        each broken requirement makes the function return false.\n";

    // Reference molecule: type 3, iface 1, state 'P', bound.
    Molecule mol = hi_make_molecule(3, { hi_make_iface(3, 1, 'P', true) });

    // (a) wrong molTypeIndex.
    {
        RxnBase::RateState rs = hi_make_rate_state({ hi_make_rxn_iface(9, 1, 'P', true) }, {});
        const bool result = hasIntangibles(0, mol, rs);
        std::cerr << "  wrong molTypeIndex  -> returned " << std::boolalpha << result << '\n';
        EXPECT_FALSE(result) << "A different molTypeIndex must not satisfy the requirement";
    }

    // (b) wrong relative interface index.
    {
        RxnBase::RateState rs = hi_make_rate_state({ hi_make_rxn_iface(3, 7, 'P', true) }, {});
        const bool result = hasIntangibles(0, mol, rs);
        std::cerr << "  wrong relIfaceIndex -> returned " << std::boolalpha << result << '\n';
        EXPECT_FALSE(result) << "A different relative interface index must not match";
    }

    // (c) wrong required state.
    {
        RxnBase::RateState rs = hi_make_rate_state({ hi_make_rxn_iface(3, 1, 'U', true) }, {});
        const bool result = hasIntangibles(0, mol, rs);
        std::cerr << "  wrong state         -> returned " << std::boolalpha << result << '\n';
        EXPECT_FALSE(result) << "A different state identity must not match";
    }

    // (d) wrong bound requirement.
    {
        RxnBase::RateState rs = hi_make_rate_state({ hi_make_rxn_iface(3, 1, 'P', false) }, {});
        const bool result = hasIntangibles(0, mol, rs);
        std::cerr << "  wrong bound status  -> returned " << std::boolalpha << result << '\n';
        EXPECT_FALSE(result) << "Requiring an unbound iface must not match a bound one";
    }
}

// -----------------------------------------------------------------------------
// Test 4: single-reactant overload reads the list selected by reactantIndex.
// -----------------------------------------------------------------------------
void test_hi_unimol_uses_requested_list_index()
{
    std::cerr << "\n[TEST] test_hi_unimol_uses_requested_list_index\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, const Molecule&, const RateState&)\n"
              << "  Scenario:    list 0 is satisfiable, list 1 is not; call with each\n"
              << "               index in turn.\n"
              << "  Pass:        index 0 -> true, index 1 -> false, proving the\n"
              << "               reactantIndex argument selects the list.\n";

    Molecule mol = hi_make_molecule(0, { hi_make_iface(0, 0, 'U', false) });

    // List 0 matches the molecule; list 1 demands an interface it does not have.
    RxnBase::RateState rateState = hi_make_rate_state(
        { hi_make_rxn_iface(0, 0, 'U', false) }, { hi_make_rxn_iface(0, 5, 'X', true) });

    const bool resultList0 = hasIntangibles(0, mol, rateState);
    const bool resultList1 = hasIntangibles(1, mol, rateState);
    std::cerr << "  reactantIndex=0 -> " << std::boolalpha << resultList0
              << ", reactantIndex=1 -> " << resultList1 << '\n';

    EXPECT_TRUE(resultList0) << "List 0 requirement is present on the molecule";
    EXPECT_FALSE(resultList1) << "List 1 requirement is absent from the molecule";
}

// -----------------------------------------------------------------------------
// Test 5: two-reactant overload, both ancillary lists empty (early return).
// -----------------------------------------------------------------------------
void test_hi_bimol_both_lists_empty_returns_true()
{
    std::cerr << "\n[TEST] test_hi_bimol_both_lists_empty_returns_true\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, int, const Molecule&, const Molecule&,\n"
              << "                              const RateState&)\n"
              << "  Scenario:    otherIfaceLists[0] and [1] are both empty, which hits the\n"
              << "               dedicated early-return branch.\n"
              << "  Pass:        function returns true even for molecules with no interfaces.\n";

    // Deliberately give the molecules no interfaces at all.
    Molecule mol1 = hi_make_molecule(0, {});
    Molecule mol2 = hi_make_molecule(1, {});

    RxnBase::RateState rateState = hi_make_rate_state({}, {});

    const bool result = hasIntangibles(0, 1, mol1, mol2, rateState);
    std::cerr << "  hasIntangibles returned " << std::boolalpha << result << '\n';
    EXPECT_TRUE(result) << "Empty ancillary lists should short-circuit to true";
}

// -----------------------------------------------------------------------------
// Test 6: two-reactant overload, both reactants satisfy their requirements.
// -----------------------------------------------------------------------------
void test_hi_bimol_both_reactants_satisfied_returns_true()
{
    std::cerr << "\n[TEST] test_hi_bimol_both_reactants_satisfied_returns_true\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, int, ...)\n"
              << "  Scenario:    reactant 1 needs (type 0, iface 1, state 'P', bound) and\n"
              << "               reactant 2 needs (type 1, iface 2, no state, unbound);\n"
              << "               both are present.\n"
              << "  Pass:        function returns true.\n";

    Molecule mol1 = hi_make_molecule(0,
        { hi_make_iface(0, 0, 'U', false), hi_make_iface(0, 1, 'P', true) });
    Molecule mol2 = hi_make_molecule(1,
        { hi_make_iface(1, 2, '\0', false), hi_make_iface(1, 3, 'A', true) });

    RxnBase::RateState rateState = hi_make_rate_state(
        { hi_make_rxn_iface(0, 1, 'P', true) }, { hi_make_rxn_iface(1, 2, '\0', false) });

    const bool result = hasIntangibles(0, 1, mol1, mol2, rateState);
    std::cerr << "  hasIntangibles returned " << std::boolalpha << result << '\n';
    EXPECT_TRUE(result) << "Both reactants carry their required ancillary interfaces";
}

// -----------------------------------------------------------------------------
// Test 7: two-reactant overload, failure caused by the first reactant.
// -----------------------------------------------------------------------------
void test_hi_bimol_first_reactant_missing_returns_false()
{
    std::cerr << "\n[TEST] test_hi_bimol_first_reactant_missing_returns_false\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, int, ...)\n"
              << "  Scenario:    reactant 2 satisfies its list, but reactant 1 lacks a\n"
              << "               required interface state.\n"
              << "  Pass:        function returns false.\n";

    // mol1's iface 1 is in state 'U', but 'P' is demanded.
    Molecule mol1 = hi_make_molecule(0, { hi_make_iface(0, 1, 'U', true) });
    Molecule mol2 = hi_make_molecule(1, { hi_make_iface(1, 2, '\0', false) });

    RxnBase::RateState rateState = hi_make_rate_state(
        { hi_make_rxn_iface(0, 1, 'P', true) }, { hi_make_rxn_iface(1, 2, '\0', false) });

    const bool result = hasIntangibles(0, 1, mol1, mol2, rateState);
    std::cerr << "  hasIntangibles returned " << std::boolalpha << result << '\n';
    EXPECT_FALSE(result) << "A missing requirement on reactant 1 must veto the reaction";
}

// -----------------------------------------------------------------------------
// Test 8: two-reactant overload, failure caused by the second reactant.
// -----------------------------------------------------------------------------
void test_hi_bimol_second_reactant_missing_returns_false()
{
    std::cerr << "\n[TEST] test_hi_bimol_second_reactant_missing_returns_false\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, int, ...)\n"
              << "  Scenario:    reactant 1 satisfies its list, but reactant 2 lacks the\n"
              << "               required bound interface.\n"
              << "  Pass:        function returns false (second loop rejects).\n";

    Molecule mol1 = hi_make_molecule(0, { hi_make_iface(0, 1, 'P', true) });
    // mol2's iface 2 is unbound while the requirement demands it be bound.
    Molecule mol2 = hi_make_molecule(1, { hi_make_iface(1, 2, '\0', false) });

    RxnBase::RateState rateState = hi_make_rate_state(
        { hi_make_rxn_iface(0, 1, 'P', true) }, { hi_make_rxn_iface(1, 2, '\0', true) });

    const bool result = hasIntangibles(0, 1, mol1, mol2, rateState);
    std::cerr << "  hasIntangibles returned " << std::boolalpha << result << '\n';
    EXPECT_FALSE(result) << "A missing requirement on reactant 2 must veto the reaction";
}

// -----------------------------------------------------------------------------
// Test 9: two-reactant overload with swapped list indices.
// -----------------------------------------------------------------------------
void test_hi_bimol_swapped_indices_map_lists_correctly()
{
    std::cerr << "\n[TEST] test_hi_bimol_swapped_indices_map_lists_correctly\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, int, ...)\n"
              << "  Scenario:    the same molecules/lists are checked with (0,1) and with\n"
              << "               (1,0); only one mapping is self-consistent.\n"
              << "  Pass:        (1,0) -> true (correct mapping), (0,1) -> false.\n";

    // mol1 owns the interface described by list 1; mol2 owns the one in list 0.
    Molecule mol1 = hi_make_molecule(0, { hi_make_iface(0, 1, 'P', true) });
    Molecule mol2 = hi_make_molecule(1, { hi_make_iface(1, 2, 'Q', false) });

    RxnBase::RateState rateState = hi_make_rate_state(
        { hi_make_rxn_iface(1, 2, 'Q', false) },   // list 0 describes mol2
        { hi_make_rxn_iface(0, 1, 'P', true) });   // list 1 describes mol1

    // Correct mapping: mol1 <-> list 1, mol2 <-> list 0.
    const bool resultSwapped = hasIntangibles(1, 0, mol1, mol2, rateState);
    // Incorrect mapping: mol1 <-> list 0, mol2 <-> list 1.
    const bool resultStraight = hasIntangibles(0, 1, mol1, mol2, rateState);

    std::cerr << "  indices (1,0) -> " << std::boolalpha << resultSwapped
              << ", indices (0,1) -> " << resultStraight << '\n';

    EXPECT_TRUE(resultSwapped) << "With indices swapped each molecule matches its own list";
    EXPECT_FALSE(resultStraight) << "With the straight mapping neither molecule matches its list";
}

// -----------------------------------------------------------------------------
// Test 10: two-reactant overload, one empty list plus one populated list.
// -----------------------------------------------------------------------------
void test_hi_bimol_one_empty_list_still_checks_other()
{
    std::cerr << "\n[TEST] test_hi_bimol_one_empty_list_still_checks_other\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, int, ...)\n"
              << "  Scenario:    list 0 is empty so the early return does NOT fire; list 1\n"
              << "               holds a requirement that is first satisfied and then not.\n"
              << "  Pass:        true when reactant 2 has the interface, false when it does not.\n";

    Molecule mol1 = hi_make_molecule(0, { hi_make_iface(0, 0, 'U', false) });

    // Case A: reactant 2 carries the required interface.
    {
        Molecule mol2 = hi_make_molecule(1, { hi_make_iface(1, 4, 'Z', true) });
        RxnBase::RateState rateState = hi_make_rate_state({}, { hi_make_rxn_iface(1, 4, 'Z', true) });
        const bool result = hasIntangibles(0, 1, mol1, mol2, rateState);
        std::cerr << "  requirement present -> returned " << std::boolalpha << result << '\n';
        EXPECT_TRUE(result) << "The only non-empty list is satisfied, so result should be true";
    }

    // Case B: reactant 2 does not carry the required interface.
    {
        Molecule mol2 = hi_make_molecule(1, { hi_make_iface(1, 4, 'Z', false) });
        RxnBase::RateState rateState = hi_make_rate_state({}, { hi_make_rxn_iface(1, 4, 'Z', true) });
        const bool result = hasIntangibles(0, 1, mol1, mol2, rateState);
        std::cerr << "  requirement absent  -> returned " << std::boolalpha << result << '\n';
        EXPECT_FALSE(result) << "An unsatisfied non-empty list must return false";
    }
}

// -----------------------------------------------------------------------------
// Test 11: multiple requirements on one reactant, only one of them missing.
// -----------------------------------------------------------------------------
void test_hi_bimol_partial_match_returns_false()
{
    std::cerr << "\n[TEST] test_hi_bimol_partial_match_returns_false\n"
              << "  Source file: src/reactions/hasIntangibles.cpp\n"
              << "  Function:    hasIntangibles(int, int, ...)\n"
              << "  Scenario:    reactant 1 must satisfy TWO ancillary requirements but\n"
              << "               only the first is present.\n"
              << "  Pass:        function returns false (ALL requirements are mandatory).\n";

    // mol1 has iface 0 (state U, unbound) but no iface 5.
    Molecule mol1 = hi_make_molecule(0, { hi_make_iface(0, 0, 'U', false) });
    Molecule mol2 = hi_make_molecule(1, { hi_make_iface(1, 0, 'U', false) });

    RxnBase::RateState rateState = hi_make_rate_state(
        { hi_make_rxn_iface(0, 0, 'U', false), hi_make_rxn_iface(0, 5, 'U', false) }, {});

    const bool result = hasIntangibles(0, 1, mol1, mol2, rateState);
    std::cerr << "  hasIntangibles returned " << std::boolalpha << result << '\n';
    EXPECT_FALSE(result) << "Every ancillary interface in the list must be found, not just one";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper runs in its own TEST so that all of
// them execute even if an earlier one records failures (only EXPECT_* is used).
// -----------------------------------------------------------------------------
TEST(HasIntangibles, UnimolEmptyRequirement) { test_hi_unimol_empty_requirement_returns_true(); }
TEST(HasIntangibles, UnimolAllRequirementsPresent) { test_hi_unimol_all_requirements_present_returns_true(); }
TEST(HasIntangibles, UnimolMismatchEachField) { test_hi_unimol_mismatch_each_field_returns_false(); }
TEST(HasIntangibles, UnimolUsesRequestedListIndex) { test_hi_unimol_uses_requested_list_index(); }
TEST(HasIntangibles, BimolBothListsEmpty) { test_hi_bimol_both_lists_empty_returns_true(); }
TEST(HasIntangibles, BimolBothReactantsSatisfied) { test_hi_bimol_both_reactants_satisfied_returns_true(); }
TEST(HasIntangibles, BimolFirstReactantMissing) { test_hi_bimol_first_reactant_missing_returns_false(); }
TEST(HasIntangibles, BimolSecondReactantMissing) { test_hi_bimol_second_reactant_missing_returns_false(); }
TEST(HasIntangibles, BimolSwappedIndices) { test_hi_bimol_swapped_indices_map_lists_correctly(); }
TEST(HasIntangibles, BimolOneEmptyList) { test_hi_bimol_one_empty_list_still_checks_other(); }
TEST(HasIntangibles, BimolPartialMatch) { test_hi_bimol_partial_match_returns_false(); }