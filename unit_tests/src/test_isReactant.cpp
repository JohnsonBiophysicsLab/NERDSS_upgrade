/*! \file test_isReactant.cpp
 *
 * ### Unit tests for src/reactions/isReactant.cpp
 *
 * That translation unit provides two overloaded predicates:
 *
 *   1. bool isReactant(const Molecule::Iface&, const Molecule&, const RxnIface&)
 *        - Compares a single concrete molecule interface against a reaction's
 *          templated reactant interface. It requires that:
 *            * bound-status matches RxnIface::requiresInteraction
 *            * the parent molecule type matches RxnIface::molTypeIndex
 *            * the absolute interface index matches RxnIface::absIfaceIndex
 *            * the interface state identity matches RxnIface::requiresState
 *
 *   2. bool isReactant(const Molecule&, const Complex&,
 *                      const CreateDestructRxn&, const std::vector<Molecule>&)
 *        - Determines whether a Molecule (possibly bound in a Complex) is the
 *          reactant of a creation/destruction reaction. Handles both the
 *          single-reactant case and the bound two-reactant case.
 *
 * Each test below builds up the minimal Molecule / Complex / CreateDestructRxn
 * state required to hit one specific code path and prints, verbosely, what is
 * being exercised and what the pass criterion is.
 */

#include <iostream>
#include <string>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "reactions/shared_reaction_functions.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with isr_ so they cannot collide with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a Molecule::Iface with fully-specified comparison fields.
 *
 * \param[in] absIndex      absolute interface (state) index -> Iface::index
 * \param[in] molTypeIndex  parent molecule template index
 * \param[in] state         state identifier character ('\0' means "no state")
 * \param[in] isBound       whether the interface is currently bound
 * \param[in] partnerIndex  moleculeList index of the bound partner (-1 if free)
 */
Molecule::Iface isr_make_iface(int absIndex, int molTypeIndex, char state, bool isBound,
    int partnerIndex = -1)
{
    Molecule::Iface iface;
    iface.index = absIndex; // absolute index of the current state
    iface.molTypeIndex = molTypeIndex;
    iface.stateIden = state;
    iface.isBound = isBound;
    iface.interaction.partnerIndex = partnerIndex;
    iface.interaction.partnerIfaceIndex = isBound ? 0 : -1;
    return iface;
}

/*! \brief Build a templated reactant interface (RxnIface) for a reaction. */
RxnIface isr_make_rxn_iface(int molTypeIndex, int absIfaceIndex, char requiresState,
    bool requiresInteraction)
{
    RxnIface rxnIface;
    rxnIface.ifaceName = "testIface";
    rxnIface.molTypeIndex = molTypeIndex;
    rxnIface.absIfaceIndex = absIfaceIndex;
    rxnIface.relIfaceIndex = 0;
    rxnIface.requiresState = requiresState;
    rxnIface.requiresInteraction = requiresInteraction;
    return rxnIface;
}

/*! \brief Build a Molecule with a given type, index, and interface list. */
Molecule isr_make_molecule(int molTypeIndex, int index, const std::vector<Molecule::Iface>& ifaces)
{
    Molecule mol;
    mol.molTypeIndex = molTypeIndex;
    mol.index = index;
    mol.interfaceList = ifaces;
    mol.myComIndex = 0;
    return mol;
}

/*! \brief Build a CreateDestructRxn::CreateDestructMol reactant descriptor. */
CreateDestructRxn::CreateDestructMol isr_make_cd_mol(int molTypeIndex,
    const std::vector<RxnIface>& ifaces)
{
    CreateDestructRxn::CreateDestructMol cdMol;
    cdMol.molTypeIndex = molTypeIndex;
    cdMol.molName = "testMol";
    cdMol.interfaceList = ifaces;
    return cdMol;
}

/*! \brief Build a Complex whose memberList contains the provided indices. */
Complex isr_make_complex(const std::vector<int>& members)
{
    Complex com;
    com.memberList = members;
    com.index = 0;
    return com;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: single-interface overload, matching free (unbound) interface.
// -----------------------------------------------------------------------------
void test_isr_iface_overload_matches_free_iface()
{
    std::cerr << "\n[TEST] test_isr_iface_overload_matches_free_iface\n"
              << "  Source file:   src/reactions/isReactant.cpp\n"
              << "  Function:      isReactant(const Molecule::Iface&, const Molecule&,\n"
              << "                            const RxnIface&)\n"
              << "  Scenario:      free interface, reaction does NOT require an interaction,\n"
              << "                 molTypeIndex / absIfaceIndex / state all agree.\n"
              << "  Pass criteria: returns true.\n";

    // A free interface of molecule type 2, absolute state index 7, state 'U'.
    Molecule::Iface iface = isr_make_iface(/*absIndex=*/7, /*molTypeIndex=*/2, /*state=*/'U',
        /*isBound=*/false);
    Molecule mol = isr_make_molecule(/*molTypeIndex=*/2, /*index=*/0, { iface });

    // Reaction template asking for exactly that free interface.
    RxnIface tempReactant = isr_make_rxn_iface(/*molTypeIndex=*/2, /*absIfaceIndex=*/7,
        /*requiresState=*/'U', /*requiresInteraction=*/false);

    const bool result = isReactant(iface, mol, tempReactant);
    std::cerr << "  isReactant returned: " << std::boolalpha << result << '\n';
    EXPECT_TRUE(result) << "A fully-matching free interface must be recognized as a reactant";
}

// -----------------------------------------------------------------------------
// Test 2: single-interface overload, matching bound interface.
// -----------------------------------------------------------------------------
void test_isr_iface_overload_matches_bound_iface()
{
    std::cerr << "\n[TEST] test_isr_iface_overload_matches_bound_iface\n"
              << "  Source file:   src/reactions/isReactant.cpp\n"
              << "  Function:      isReactant(Iface, Molecule, RxnIface)\n"
              << "  Scenario:      bound interface, reaction requires an interaction.\n"
              << "  Pass criteria: returns true.\n";

    Molecule::Iface iface = isr_make_iface(/*absIndex=*/3, /*molTypeIndex=*/1, /*state=*/'\0',
        /*isBound=*/true, /*partnerIndex=*/5);
    Molecule mol = isr_make_molecule(/*molTypeIndex=*/1, /*index=*/4, { iface });

    RxnIface tempReactant = isr_make_rxn_iface(/*molTypeIndex=*/1, /*absIfaceIndex=*/3,
        /*requiresState=*/'\0', /*requiresInteraction=*/true);

    const bool result = isReactant(iface, mol, tempReactant);
    std::cerr << "  isReactant returned: " << std::boolalpha << result << '\n';
    EXPECT_TRUE(result) << "A matching bound interface must be recognized as a reactant";
}

// -----------------------------------------------------------------------------
// Test 3: single-interface overload, all of the individual mismatch paths.
// -----------------------------------------------------------------------------
void test_isr_iface_overload_rejects_mismatches()
{
    std::cerr << "\n[TEST] test_isr_iface_overload_rejects_mismatches\n"
              << "  Source file:   src/reactions/isReactant.cpp\n"
              << "  Function:      isReactant(Iface, Molecule, RxnIface)\n"
              << "  Scenario:      exercise every rejection branch one at a time.\n"
              << "  Pass criteria: each mismatched comparison returns false.\n";

    // Baseline: a free interface that would otherwise match.
    const Molecule::Iface freeIface = isr_make_iface(7, 2, 'U', false);
    const Molecule mol = isr_make_molecule(2, 0, { freeIface });

    // (a) Reaction requires a bond, but the interface is free.
    std::cerr << "  -> (a) requiresInteraction=true but interface is unbound\n";
    RxnIface needsBond = isr_make_rxn_iface(2, 7, 'U', true);
    EXPECT_FALSE(isReactant(freeIface, mol, needsBond))
        << "Unbound interface cannot satisfy a reaction requiring an interaction";

    // (b) Reaction requires a free interface, but the interface is bound.
    std::cerr << "  -> (b) requiresInteraction=false but interface is bound\n";
    const Molecule::Iface boundIface = isr_make_iface(7, 2, 'U', true, 9);
    RxnIface needsFree = isr_make_rxn_iface(2, 7, 'U', false);
    EXPECT_FALSE(isReactant(boundIface, mol, needsFree))
        << "Bound interface cannot satisfy a reaction requiring a free interface";

    // (c) Wrong parent molecule type.
    std::cerr << "  -> (c) molTypeIndex mismatch (mol type 2 vs rxn type 3)\n";
    RxnIface wrongMolType = isr_make_rxn_iface(3, 7, 'U', false);
    EXPECT_FALSE(isReactant(freeIface, mol, wrongMolType))
        << "Different molecule type must not be a reactant";

    // (d) Wrong absolute interface index.
    std::cerr << "  -> (d) absIfaceIndex mismatch (iface index 7 vs rxn index 8)\n";
    RxnIface wrongAbsIndex = isr_make_rxn_iface(2, 8, 'U', false);
    EXPECT_FALSE(isReactant(freeIface, mol, wrongAbsIndex))
        << "Different absolute interface index must not be a reactant";

    // (e) Wrong state identifier.
    std::cerr << "  -> (e) requiresState mismatch (iface 'U' vs rxn 'P')\n";
    RxnIface wrongState = isr_make_rxn_iface(2, 7, 'P', false);
    EXPECT_FALSE(isReactant(freeIface, mol, wrongState))
        << "Different interface state must not be a reactant";
}

// -----------------------------------------------------------------------------
// Test 4: CreateDestructRxn overload, single reactant, everything matches.
// -----------------------------------------------------------------------------
void test_isr_createdestruct_single_reactant_match()
{
    std::cerr << "\n[TEST] test_isr_createdestruct_single_reactant_match\n"
              << "  Source file:   src/reactions/isReactant.cpp\n"
              << "  Function:      isReactant(Molecule, Complex, CreateDestructRxn,\n"
              << "                            std::vector<Molecule>)\n"
              << "  Scenario:      reaction has exactly one reactant molecule, and the\n"
              << "                 candidate molecule matches type plus both interfaces.\n"
              << "  Pass criteria: returns true.\n";

    // Candidate molecule: type 0, two free interfaces (abs indices 0 and 1).
    Molecule mol = isr_make_molecule(0, 0,
        { isr_make_iface(0, 0, '\0', false), isr_make_iface(1, 0, '\0', false) });
    std::vector<Molecule> moleculeList { mol };
    Complex com = isr_make_complex({ 0 });

    // Reaction: destroys a single free molecule of type 0 with two free interfaces.
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::destruction;
    rxn.reactantMolList.push_back(isr_make_cd_mol(0,
        { isr_make_rxn_iface(0, 0, '\0', false), isr_make_rxn_iface(0, 1, '\0', false) }));

    const bool result = isReactant(mol, com, rxn, moleculeList);
    std::cerr << "  isReactant returned: " << std::boolalpha << result << '\n';
    EXPECT_TRUE(result) << "A free monomer whose interfaces all match must be a reactant";
}

// -----------------------------------------------------------------------------
// Test 5: CreateDestructRxn overload, single reactant, rejection branches.
// -----------------------------------------------------------------------------
void test_isr_createdestruct_single_reactant_mismatch()
{
    std::cerr << "\n[TEST] test_isr_createdestruct_single_reactant_mismatch\n"
              << "  Source file:   src/reactions/isReactant.cpp\n"
              << "  Function:      isReactant(Molecule, Complex, CreateDestructRxn, list)\n"
              << "  Scenario:      wrong molecule type, wrong interface count, and a\n"
              << "                 mismatched interface state.\n"
              << "  Pass criteria: all three cases return false.\n";

    Complex com = isr_make_complex({ 0 });

    // Reference reaction: one reactant of type 0 with two free, stateless ifaces.
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::destruction;
    rxn.reactantMolList.push_back(isr_make_cd_mol(0,
        { isr_make_rxn_iface(0, 0, '\0', false), isr_make_rxn_iface(0, 1, '\0', false) }));

    // (a) Wrong molecule type index.
    std::cerr << "  -> (a) candidate molecule is of type 1, reaction expects type 0\n";
    Molecule wrongType = isr_make_molecule(1, 0,
        { isr_make_iface(0, 1, '\0', false), isr_make_iface(1, 1, '\0', false) });
    std::vector<Molecule> listA { wrongType };
    EXPECT_FALSE(isReactant(wrongType, com, rxn, listA))
        << "Molecule of a different type must not be the reactant";

    // (b) Wrong number of interfaces.
    std::cerr << "  -> (b) candidate molecule has 1 interface, reaction expects 2\n";
    Molecule wrongCount = isr_make_molecule(0, 0, { isr_make_iface(0, 0, '\0', false) });
    std::vector<Molecule> listB { wrongCount };
    EXPECT_FALSE(isReactant(wrongCount, com, rxn, listB))
        << "Interface-count mismatch must reject the molecule";

    // (c) Correct type/count, but one interface is in the wrong state.
    std::cerr << "  -> (c) second interface has state 'P', reaction expects '\\0'\n";
    Molecule wrongState = isr_make_molecule(0, 0,
        { isr_make_iface(0, 0, '\0', false), isr_make_iface(1, 0, 'P', false) });
    std::vector<Molecule> listC { wrongState };
    EXPECT_FALSE(isReactant(wrongState, com, rxn, listC))
        << "State mismatch on any interface must reject the molecule";
}

// -----------------------------------------------------------------------------
// Test 6: CreateDestructRxn overload, bound dimer reactant (two reactant mols).
// -----------------------------------------------------------------------------
void test_isr_createdestruct_bound_dimer_match()
{
    std::cerr << "\n[TEST] test_isr_createdestruct_bound_dimer_match\n"
              << "  Source file:   src/reactions/isReactant.cpp\n"
              << "  Function:      isReactant(Molecule, Complex, CreateDestructRxn, list)\n"
              << "  Scenario:      reaction has two reactant molecules; the candidate is\n"
              << "                 molecule 0 of a bound A(0)-B(1) dimer.\n"
              << "  Pass criteria: returns true when queried with the lower-index member\n"
              << "                 (the implementation only walks partners whose index is\n"
              << "                 larger, avoiding double counting).\n";

    // mol0: type 0, index 0, one interface bound to molecule 1.
    Molecule mol0 = isr_make_molecule(0, 0, { isr_make_iface(0, 0, '\0', true, /*partner=*/1) });
    // mol1: type 1, index 1, one interface bound back to molecule 0.
    Molecule mol1 = isr_make_molecule(1, 1, { isr_make_iface(1, 1, '\0', true, /*partner=*/0) });
    std::vector<Molecule> moleculeList { mol0, mol1 };

    // Both molecules live in the same two-member complex.
    Complex com = isr_make_complex({ 0, 1 });

    // Reaction destroys the bound dimer: both reactant interfaces require a bond.
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::destruction;
    rxn.reactantMolList.push_back(isr_make_cd_mol(0, { isr_make_rxn_iface(0, 0, '\0', true) }));
    rxn.reactantMolList.push_back(isr_make_cd_mol(1, { isr_make_rxn_iface(1, 1, '\0', true) }));

    // Querying with the lower-index member should succeed.
    const bool resultLow = isReactant(moleculeList[0], com, rxn, moleculeList);
    std::cerr << "  isReactant(mol0, ...) returned: " << std::boolalpha << resultLow << '\n';
    EXPECT_TRUE(resultLow) << "Bound dimer should be detected from its lower-index member";

    // Querying with the higher-index member is expected to return false, since
    // the implementation only inspects partners with a strictly larger index.
    const bool resultHigh = isReactant(moleculeList[1], com, rxn, moleculeList);
    std::cerr << "  isReactant(mol1, ...) returned: " << std::boolalpha << resultHigh
              << " (documented one-way behaviour)\n";
    EXPECT_FALSE(resultHigh)
        << "Higher-index member yields false because only larger partner indices are walked";
}

// -----------------------------------------------------------------------------
// Test 7: CreateDestructRxn overload, bound-dimer rejection branches.
// -----------------------------------------------------------------------------
void test_isr_createdestruct_bound_dimer_mismatch()
{
    std::cerr << "\n[TEST] test_isr_createdestruct_bound_dimer_mismatch\n"
              << "  Source file:   src/reactions/isReactant.cpp\n"
              << "  Function:      isReactant(Molecule, Complex, CreateDestructRxn, list)\n"
              << "  Scenario:      (a) candidate type matches neither reactant,\n"
              << "                 (b) candidate interface is free while the reaction needs\n"
              << "                     a bond,\n"
              << "                 (c) a two-reactant reaction is queried against a\n"
              << "                     single-member complex.\n"
              << "  Pass criteria: all three cases return false.\n";

    // Two-reactant destruction reaction of a type-0 / type-1 dimer.
    CreateDestructRxn rxn;
    rxn.rxnType = ReactionType::destruction;
    rxn.reactantMolList.push_back(isr_make_cd_mol(0, { isr_make_rxn_iface(0, 0, '\0', true) }));
    rxn.reactantMolList.push_back(isr_make_cd_mol(1, { isr_make_rxn_iface(1, 1, '\0', true) }));

    // (a) Candidate molecule is type 5, which is neither reactant type.
    std::cerr << "  -> (a) candidate molecule type 5 matches neither reactant (0 or 1)\n";
    Molecule strangerA = isr_make_molecule(5, 0, { isr_make_iface(9, 5, '\0', true, 1) });
    Molecule partnerA = isr_make_molecule(1, 1, { isr_make_iface(1, 1, '\0', true, 0) });
    std::vector<Molecule> listA { strangerA, partnerA };
    Complex comA = isr_make_complex({ 0, 1 });
    EXPECT_FALSE(isReactant(listA[0], comA, rxn, listA))
        << "A molecule matching neither reactant type must be rejected";

    // (b) Candidate is the right type, but its interface is unbound.
    std::cerr << "  -> (b) candidate type 0 has an unbound interface, reaction needs a bond\n";
    Molecule freeMol = isr_make_molecule(0, 0, { isr_make_iface(0, 0, '\0', false) });
    Molecule otherMol = isr_make_molecule(1, 1, { isr_make_iface(1, 1, '\0', true, 0) });
    std::vector<Molecule> listB { freeMol, otherMol };
    Complex comB = isr_make_complex({ 0, 1 });
    EXPECT_FALSE(isReactant(listB[0], comB, rxn, listB))
        << "Unbound interface cannot satisfy a bound-dimer reactant";

    // (c) The complex contains only one molecule, so the bound branch is skipped.
    std::cerr << "  -> (c) two-reactant reaction queried against a monomeric complex\n";
    Molecule loneMol = isr_make_molecule(0, 0, { isr_make_iface(0, 0, '\0', false) });
    std::vector<Molecule> listC { loneMol };
    Complex comC = isr_make_complex({ 0 });
    EXPECT_FALSE(isReactant(listC[0], comC, rxn, listC))
        << "A monomeric complex cannot satisfy a two-molecule reactant";
}

// -----------------------------------------------------------------------------
// Test 8: CreateDestructRxn overload, degenerate reaction with no reactants.
// -----------------------------------------------------------------------------
void test_isr_createdestruct_no_reactants()
{
    std::cerr << "\n[TEST] test_isr_createdestruct_no_reactants\n"
              << "  Source file:   src/reactions/isReactant.cpp\n"
              << "  Function:      isReactant(Molecule, Complex, CreateDestructRxn, list)\n"
              << "  Scenario:      zeroth-order creation reaction (empty reactantMolList).\n"
              << "  Pass criteria: returns false because neither size branch is taken.\n";

    Molecule mol = isr_make_molecule(0, 0, { isr_make_iface(0, 0, '\0', false) });
    std::vector<Molecule> moleculeList { mol };
    Complex com = isr_make_complex({ 0 });

    CreateDestructRxn rxn; // reactantMolList intentionally left empty
    rxn.rxnType = ReactionType::zerothOrderCreation;

    const bool result = isReactant(moleculeList[0], com, rxn, moleculeList);
    std::cerr << "  isReactant returned: " << std::boolalpha << result << '\n';
    EXPECT_FALSE(result) << "A reaction with no reactant molecules can have no reactant match";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper is invoked from its own TEST so
// that a failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(IsReactantIfaceOverload, MatchesFreeIface) { test_isr_iface_overload_matches_free_iface(); }
TEST(IsReactantIfaceOverload, MatchesBoundIface) { test_isr_iface_overload_matches_bound_iface(); }
TEST(IsReactantIfaceOverload, RejectsMismatches) { test_isr_iface_overload_rejects_mismatches(); }
TEST(IsReactantCreateDestruct, SingleReactantMatch) { test_isr_createdestruct_single_reactant_match(); }
TEST(IsReactantCreateDestruct, SingleReactantMismatch) { test_isr_createdestruct_single_reactant_mismatch(); }
TEST(IsReactantCreateDestruct, BoundDimerMatch) { test_isr_createdestruct_bound_dimer_match(); }
TEST(IsReactantCreateDestruct, BoundDimerMismatch) { test_isr_createdestruct_bound_dimer_mismatch(); }
TEST(IsReactantCreateDestruct, NoReactants) { test_isr_createdestruct_no_reactants(); }