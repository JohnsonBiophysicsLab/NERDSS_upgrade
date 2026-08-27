/*! \file test_check_implicit_reactions.cpp
 *
 * ### Unit test for src/reactions/check_implicit_reactions.cpp
 *
 * The file under test contains exactly one function:
 *
 * \code
 * void check_implicit_reactions(int pro1Index, int pro2Index, int simItr,
 *                               const Parameters&, std::vector<Molecule>&,
 *                               std::vector<Complex>&,
 *                               const std::vector<MolTemplate>&,
 *                               const std::vector<ForwardRxn>&,
 *                               const std::vector<BackRxn>&, copyCounters&,
 *                               Membrane&, std::vector<double>&,
 *                               std::vector<double>&, std::vector<double>&);
 * \endcode
 *
 * The routine is a *guard/dispatch* function: it walks the free interfaces of
 * protein `pro1Index`, checks whether any of them can react with the implicit
 * lipid `pro2Index`, and - only when an actual reaction is found - delegates to
 * determine_2D_implicitlipid_reaction_probability() or
 * determine_3D_implicitlipid_reaction_probability().
 *
 * Those delegate routines are the only place where observable state is written
 * (they append to Molecule::probvec / crossbase / mycrossint / crossrxn and
 * increment Complex::ncross).  Therefore every guard branch of
 * check_implicit_reactions() can be verified by asserting that *no* reaction
 * bookkeeping was recorded.
 *
 * The branches exercised here are:
 *   1. pro2 is not an implicit lipid                       -> immediate return
 *   2. no free implicit lipids left in the membrane        -> canInteract=false
 *   3. pro1 has no free interfaces                         -> outer loop empty
 *   4. pro1's interface state lists no reaction partners   -> middle loop empty
 *   5. the listed partner index does not match the lipid   -> inner `if` false
 *   6. the partner index matches, but no reaction is
 *      registered in forwardRxns/backRxns                  -> rxnIndex stays -1
 *      (this also verifies the temporary +=/-= of the
 *       implicit-lipid interface index is undone again)
 *
 * NOTE: the "reaction really happens" path is deliberately *not* driven here.
 * It would require a fully populated ForwardRxn/2D lookup-table machinery, and
 * an incomplete setup there would abort the whole gtest binary rather than fail
 * a single case.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// The global GSL random number generator is defined in gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Make sure the global RNG exists.  check_implicit_reactions() itself does not
// draw random numbers on the paths tested here, but any helper it might reach
// would, so we guarantee a deterministic, seeded generator.
// -----------------------------------------------------------------------------
void cir_ensure_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}

/*! \brief Build a minimal soluble protein Molecule with one free interface.
 *
 * \param[in] index          index of this Molecule inside moleculeList
 * \param[in] comIndex       index of its parent Complex inside complexList
 * \param[in] molTypeIndex   index of its MolTemplate
 * \param[in] com            centre-of-mass coordinate
 * \param[in] absIfaceIndex  absolute (state) index of its single interface
 */
Molecule cir_make_protein(int index, int comIndex, int molTypeIndex, const Coord& com, int absIfaceIndex)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;

    Molecule::Iface iface;
    iface.coord = Coord { com.x, com.y, com.z - 1.0 }; // interface offset from COM
    iface.index = absIfaceIndex;
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.stateIden = '\0';
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    // interface 0 is free and therefore available for binding
    mol.freelist.push_back(0);
    return mol;
}

/*! \brief Build the implicit-lipid Molecule (single interface, always free). */
Molecule cir_make_implicit_lipid(int index, int comIndex, int molTypeIndex, int absIfaceIndex)
{
    Molecule mol = cir_make_protein(index, comIndex, molTypeIndex, Coord { 0.0, 0.0, -50.0 }, absIfaceIndex);
    mol.isImplicitLipid = true;
    mol.isLipid = true;
    // The implicit lipid interface sits exactly on the membrane plane.
    mol.interfaceList[0].coord = mol.comCoord;
    return mol;
}

/*! \brief Build a Complex owning exactly one member Molecule. */
Complex cir_make_complex(int index, int memberMolIndex, const Coord& com, const Coord& D)
{
    Complex targCom;
    targCom.index = index;
    targCom.id = index;
    targCom.comCoord = com;
    targCom.D = D;
    targCom.Dr = Coord { 0.01, 0.01, 0.01 };
    targCom.mass = 1.0;
    targCom.radius = 1.0;
    targCom.isEmpty = false;
    targCom.ncross = 0;
    targCom.memberList.push_back(memberMolIndex);
    targCom.numEachMol = std::vector<int>(2, 0);
    targCom.lastNumberUpdateItrEachMol = std::vector<long long int>(2, 0);
    targCom.OnSurface = (D.z < 1E-10);
    return targCom;
}

/*! \brief Build a MolTemplate with a single interface carrying a single state.
 *
 * \param[in] typeIndex      molTypeIndex of the template
 * \param[in] name           molecule name
 * \param[in] absIfaceIndex  absolute index given to the single interface state
 * \param[in] isIL           is this template the implicit lipid?
 * \param[in] partners       list of absolute interface indices this state can
 *                           react with (Interface::State::rxnPartners)
 */
MolTemplate cir_make_template(int typeIndex, const std::string& name, int absIfaceIndex, bool isIL,
    const std::vector<unsigned>& partners)
{
    MolTemplate temp;
    temp.molTypeIndex = typeIndex;
    temp.molName = name;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.copies = 1;
    temp.isImplicitLipid = isIL;
    temp.isLipid = isIL;
    temp.D = isIL ? Coord { 1.0, 1.0, 0.0 } : Coord { 10.0, 10.0, 10.0 };
    temp.Dr = Coord { 0.01, 0.01, 0.01 };

    Interface iface;
    iface.index = 0;
    iface.name = name + "site";
    iface.iCoord = Coord { 0.0, 0.0, -1.0 };

    Interface::State state;
    state.index = absIfaceIndex;
    state.iden = '\0';
    state.ifaceAndStateName = iface.name;
    state.rxnPartners = partners; // <- drives the middle loop of the function
    iface.stateList.push_back(state);

    temp.interfaceList.push_back(iface);
    return temp;
}

/*! \brief Assert that no reaction bookkeeping was recorded anywhere.
 *
 * determine_2D/3D_implicitlipid_reaction_probability() are the only writers of
 * these containers, so "all empty / ncross == 0" is exactly the criterion for
 * "check_implicit_reactions() bailed out without evaluating a reaction".
 */
void cir_expect_no_reaction_recorded(const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList, const std::string& context)
{
    for (const auto& mol : moleculeList) {
        EXPECT_TRUE(mol.probvec.empty())
            << context << ": probvec of molecule " << mol.index << " should remain empty";
        EXPECT_TRUE(mol.crossbase.empty())
            << context << ": crossbase of molecule " << mol.index << " should remain empty";
        EXPECT_TRUE(mol.mycrossint.empty())
            << context << ": mycrossint of molecule " << mol.index << " should remain empty";
        EXPECT_TRUE(mol.crossrxn.empty())
            << context << ": crossrxn of molecule " << mol.index << " should remain empty";
    }
    for (const auto& com : complexList) {
        EXPECT_EQ(com.ncross, 0)
            << context << ": ncross of complex " << com.index << " should remain 0";
    }
}

/*! \brief Container holding a complete, self-consistent mini-system.
 *
 * moleculeList[0] : soluble protein   (pro1)
 * moleculeList[1] : implicit lipid    (pro2)
 */
struct CirSystem {
    Parameters params {};
    Membrane membraneObject {};
    copyCounters counterArrays {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<double> IL2DbindingVec {};
    std::vector<double> IL2DUnbindingVec {};
    std::vector<double> ILTableIDs {};
};

/*! \brief Build the default mini-system.
 *
 * \param[in] partners    rxnPartners of the protein's interface state
 * \param[in] nStates     number of implicit-lipid states
 * \param[in] freeLipids  number of free lipids for each state
 * \param[in] pro2IsIL    should molecule 1 be flagged as an implicit lipid?
 */
CirSystem cir_build_system(const std::vector<unsigned>& partners, int nStates, int freeLipids, bool pro2IsIL)
{
    CirSystem sys;

    // --- simulation parameters -------------------------------------------
    sys.params.timeStep = 1.0;
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 2;

    // --- boundary: a 100 nm cubic water box, flat membrane ---------------
    sys.membraneObject.isSphere = false;
    sys.membraneObject.isBox = true;
    sys.membraneObject.implicitLipid = true;
    sys.membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    sys.membraneObject.nStates = nStates;
    sys.membraneObject.numberOfFreeLipidsEachState = std::vector<int>(nStates, freeLipids);
    sys.membraneObject.numberOfProteinEachState = std::vector<int>(nStates, 0);
    sys.membraneObject.implicitlipidIndex = 1;
    sys.membraneObject.RS3Dvect = std::vector<double>(500, 0.0);

    // --- templates: protein "A" (abs iface 0), implicit lipid "IL" (abs 1)
    sys.molTemplateList.push_back(cir_make_template(0, "A", 0, false, partners));
    sys.molTemplateList.push_back(cir_make_template(1, "IL", 1, true, std::vector<unsigned> { 0u }));
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType = std::vector<int> { 1, 1 };

    // --- molecules --------------------------------------------------------
    sys.moleculeList.push_back(cir_make_protein(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, 0));
    Molecule lipid = cir_make_implicit_lipid(1, 1, 1, 1);
    lipid.isImplicitLipid = pro2IsIL; // toggled for the "not a lipid" test
    sys.moleculeList.push_back(lipid);

    // --- complexes --------------------------------------------------------
    // Complex 0 is in solution (D.z != 0)  -> would take the 3D branch.
    sys.complexList.push_back(cir_make_complex(0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 10.0, 10.0, 10.0 }));
    // Complex 1 is the implicit-lipid surface (D.z == 0).
    sys.complexList.push_back(cir_make_complex(1, 1, Coord { 0.0, 0.0, -50.0 }, Coord { 1.0, 1.0, 0.0 }));

    return sys;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: pro2 is NOT an implicit lipid -> the function must return immediately.
// -----------------------------------------------------------------------------
void test_cir_returns_when_partner_is_not_implicit_lipid()
{
    std::cerr << "\n[TEST] test_cir_returns_when_partner_is_not_implicit_lipid\n"
              << "  Source file:   src/reactions/check_implicit_reactions.cpp\n"
              << "  Function:      check_implicit_reactions()\n"
              << "  Scenario:      molecule 1 is an ordinary protein, not the implicit lipid.\n"
              << "  Pass criteria: the guard `if (!isImplicitLipid) return;` fires, so no\n"
              << "                 reaction bookkeeping (probvec/crossbase/ncross) is written.\n";

    cir_ensure_rng();

    // rxnPartners deliberately DOES match the lipid interface (absIndex 1) so
    // that the only thing preventing evaluation is the isImplicitLipid flag.
    CirSystem sys = cir_build_system(std::vector<unsigned> { 1u }, 1, 100, /*pro2IsIL=*/false);

    std::cerr << "  Calling check_implicit_reactions(pro1=0, pro2=1, simItr=1)...\n";
    check_implicit_reactions(0, 1, 1, sys.params, sys.moleculeList, sys.complexList, sys.molTemplateList,
        sys.forwardRxns, sys.backRxns, sys.counterArrays, sys.membraneObject, sys.IL2DbindingVec,
        sys.IL2DUnbindingVec, sys.ILTableIDs);

    cir_expect_no_reaction_recorded(sys.moleculeList, sys.complexList, "non-implicit-lipid partner");

    // The interface indices must be untouched by the early return.
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 0)
        << "protein interface absolute index must be unchanged";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].index, 1)
        << "partner interface absolute index must be unchanged";

    std::cerr << "  Observed: probvec size = " << sys.moleculeList[0].probvec.size()
              << ", complex0.ncross = " << sys.complexList[0].ncross << "\n";
}

// -----------------------------------------------------------------------------
// Test 2: implicit lipid present, but no free lipids left -> canInteract=false.
// -----------------------------------------------------------------------------
void test_cir_no_free_lipids_blocks_interaction()
{
    std::cerr << "\n[TEST] test_cir_no_free_lipids_blocks_interaction\n"
              << "  Source file:   src/reactions/check_implicit_reactions.cpp\n"
              << "  Function:      check_implicit_reactions()\n"
              << "  Scenario:      membraneObject.numberOfFreeLipidsEachState sums to 0.\n"
              << "  Pass criteria: canInteract becomes false, so the association-probability\n"
              << "                 block is skipped entirely (no bookkeeping written).\n";

    cir_ensure_rng();

    // Matching partner, valid implicit lipid, but ZERO free lipids.
    CirSystem sys = cir_build_system(std::vector<unsigned> { 1u }, 1, /*freeLipids=*/0, /*pro2IsIL=*/true);

    std::cerr << "  Free lipid count = " << sys.membraneObject.numberOfFreeLipidsEachState[0] << "\n";
    std::cerr << "  Calling check_implicit_reactions(pro1=0, pro2=1, simItr=1)...\n";
    check_implicit_reactions(0, 1, 1, sys.params, sys.moleculeList, sys.complexList, sys.molTemplateList,
        sys.forwardRxns, sys.backRxns, sys.counterArrays, sys.membraneObject, sys.IL2DbindingVec,
        sys.IL2DUnbindingVec, sys.ILTableIDs);

    cir_expect_no_reaction_recorded(sys.moleculeList, sys.complexList, "no free lipids");

    // The routine must not consume/modify the free-lipid bookkeeping either.
    EXPECT_EQ(sys.membraneObject.numberOfFreeLipidsEachState[0], 0)
        << "free lipid counter must not be modified by check_implicit_reactions";

    std::cerr << "  Observed: probvec size = " << sys.moleculeList[0].probvec.size() << " (expected 0)\n";
}

// -----------------------------------------------------------------------------
// Test 3: protein has no free interfaces -> the outer `freelist` loop is empty.
// -----------------------------------------------------------------------------
void test_cir_protein_without_free_interfaces()
{
    std::cerr << "\n[TEST] test_cir_protein_without_free_interfaces\n"
              << "  Source file:   src/reactions/check_implicit_reactions.cpp\n"
              << "  Function:      check_implicit_reactions()\n"
              << "  Scenario:      moleculeList[0].freelist is empty (all sites occupied).\n"
              << "  Pass criteria: the loop over pro1's free interfaces never executes, so\n"
              << "                 nothing is recorded.\n";

    cir_ensure_rng();

    CirSystem sys = cir_build_system(std::vector<unsigned> { 1u }, 1, 100, true);
    sys.moleculeList[0].freelist.clear(); // no free binding sites on the protein
    sys.moleculeList[0].interfaceList[0].isBound = true;

    std::cerr << "  pro1 freelist size = " << sys.moleculeList[0].freelist.size() << "\n";
    std::cerr << "  Calling check_implicit_reactions(pro1=0, pro2=1, simItr=5)...\n";
    check_implicit_reactions(0, 1, 5, sys.params, sys.moleculeList, sys.complexList, sys.molTemplateList,
        sys.forwardRxns, sys.backRxns, sys.counterArrays, sys.membraneObject, sys.IL2DbindingVec,
        sys.IL2DUnbindingVec, sys.ILTableIDs);

    cir_expect_no_reaction_recorded(sys.moleculeList, sys.complexList, "protein has no free interfaces");

    std::cerr << "  Observed: complex0.ncross = " << sys.complexList[0].ncross << " (expected 0)\n";
}

// -----------------------------------------------------------------------------
// Test 4: the protein's interface state declares no reaction partners.
// -----------------------------------------------------------------------------
void test_cir_state_without_reaction_partners()
{
    std::cerr << "\n[TEST] test_cir_state_without_reaction_partners\n"
              << "  Source file:   src/reactions/check_implicit_reactions.cpp\n"
              << "  Function:      check_implicit_reactions()\n"
              << "  Scenario:      Interface::State::rxnPartners of pro1 is empty.\n"
              << "  Pass criteria: the `for (auto statePartner : state.rxnPartners)` loop is\n"
              << "                 never entered, so no probability is evaluated.\n";

    cir_ensure_rng();

    CirSystem sys = cir_build_system(std::vector<unsigned> {}, 1, 100, true);

    std::cerr << "  rxnPartners size = "
              << sys.molTemplateList[0].interfaceList[0].stateList[0].rxnPartners.size() << "\n";
    std::cerr << "  Calling check_implicit_reactions(pro1=0, pro2=1, simItr=7)...\n";
    check_implicit_reactions(0, 1, 7, sys.params, sys.moleculeList, sys.complexList, sys.molTemplateList,
        sys.forwardRxns, sys.backRxns, sys.counterArrays, sys.membraneObject, sys.IL2DbindingVec,
        sys.IL2DUnbindingVec, sys.ILTableIDs);

    cir_expect_no_reaction_recorded(sys.moleculeList, sys.complexList, "state has no rxnPartners");

    std::cerr << "  Observed: probvec size = " << sys.moleculeList[0].probvec.size() << " (expected 0)\n";
}

// -----------------------------------------------------------------------------
// Test 5: the declared partner index does not match the implicit-lipid state.
// -----------------------------------------------------------------------------
void test_cir_partner_index_mismatch()
{
    std::cerr << "\n[TEST] test_cir_partner_index_mismatch\n"
              << "  Source file:   src/reactions/check_implicit_reactions.cpp\n"
              << "  Function:      check_implicit_reactions()\n"
              << "  Scenario:      rxnPartners = {99} while the implicit lipid interface has\n"
              << "                 absolute index 1 and nStates = 1.\n"
              << "  Pass criteria: `absIface2 == statePartner - stateIdx` is false for every\n"
              << "                 state, so find_which_reaction() is never called and the\n"
              << "                 lipid interface index is left exactly as it was.\n";

    cir_ensure_rng();

    CirSystem sys = cir_build_system(std::vector<unsigned> { 99u }, 1, 100, true);

    const int lipidIfaceIndexBefore = sys.moleculeList[1].interfaceList[0].index;
    std::cerr << "  lipid interface absolute index before call = " << lipidIfaceIndexBefore << "\n";
    std::cerr << "  Calling check_implicit_reactions(pro1=0, pro2=1, simItr=11)...\n";
    check_implicit_reactions(0, 1, 11, sys.params, sys.moleculeList, sys.complexList, sys.molTemplateList,
        sys.forwardRxns, sys.backRxns, sys.counterArrays, sys.membraneObject, sys.IL2DbindingVec,
        sys.IL2DUnbindingVec, sys.ILTableIDs);

    cir_expect_no_reaction_recorded(sys.moleculeList, sys.complexList, "partner index mismatch");

    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].index, lipidIfaceIndexBefore)
        << "implicit-lipid interface index must be unchanged when no partner matches";

    std::cerr << "  lipid interface absolute index after call  = "
              << sys.moleculeList[1].interfaceList[0].index << "\n";
}

// -----------------------------------------------------------------------------
// Test 6: partner index matches for the second lipid state, but there is no
//         registered reaction.  This also checks that the temporary
//         `interfaceList[relIface2].index += stateIdx` is undone afterwards.
// -----------------------------------------------------------------------------
void test_cir_matching_partner_without_registered_reaction()
{
    std::cerr << "\n[TEST] test_cir_matching_partner_without_registered_reaction\n"
              << "  Source file:   src/reactions/check_implicit_reactions.cpp\n"
              << "  Function:      check_implicit_reactions()\n"
              << "  Scenario:      nStates = 2, rxnPartners = {2}; the match therefore occurs\n"
              << "                 for implicit-lipid state 1 (2 - 1 == absIface2 == 1) and the\n"
              << "                 lipid interface index is temporarily shifted by +1.\n"
              << "                 forwardRxns/backRxns are empty and the protein state lists\n"
              << "                 no reactions, so find_which_reaction() must report -1.\n"
              << "  Pass criteria: rxnIndex stays -1 -> no probability evaluated, AND the\n"
              << "                 lipid interface index is restored to its original value.\n";

    cir_ensure_rng();

    CirSystem sys = cir_build_system(std::vector<unsigned> { 2u }, /*nStates=*/2, /*freeLipids=*/50, true);

    // The protein state deliberately declares no forward/state-change reactions,
    // so find_which_reaction() has nothing to find.
    EXPECT_TRUE(sys.molTemplateList[0].interfaceList[0].stateList[0].myForwardRxns.empty())
        << "test setup: the protein state must not own any forward reactions";
    EXPECT_TRUE(sys.forwardRxns.empty()) << "test setup: forwardRxns must be empty";
    EXPECT_TRUE(sys.backRxns.empty()) << "test setup: backRxns must be empty";

    const int lipidIfaceIndexBefore = sys.moleculeList[1].interfaceList[0].index;
    std::cerr << "  lipid interface absolute index before call = " << lipidIfaceIndexBefore << "\n";
    std::cerr << "  Calling check_implicit_reactions(pro1=0, pro2=1, simItr=13)...\n";
    check_implicit_reactions(0, 1, 13, sys.params, sys.moleculeList, sys.complexList, sys.molTemplateList,
        sys.forwardRxns, sys.backRxns, sys.counterArrays, sys.membraneObject, sys.IL2DbindingVec,
        sys.IL2DUnbindingVec, sys.ILTableIDs);

    cir_expect_no_reaction_recorded(sys.moleculeList, sys.complexList, "matching partner, no reaction defined");

    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].index, lipidIfaceIndexBefore)
        << "the temporary state-index offset applied to the implicit-lipid interface "
           "must be subtracted again before returning";

    // The free lipid populations are read-only for this routine.
    EXPECT_EQ(sys.membraneObject.numberOfFreeLipidsEachState[0], 50)
        << "free lipid count of state 0 must be untouched";
    EXPECT_EQ(sys.membraneObject.numberOfFreeLipidsEachState[1], 50)
        << "free lipid count of state 1 must be untouched";

    std::cerr << "  lipid interface absolute index after call  = "
              << sys.moleculeList[1].interfaceList[0].index << "\n";
    std::cerr << "  Observed: probvec size = " << sys.moleculeList[0].probvec.size()
              << ", complex0.ncross = " << sys.complexList[0].ncross << " (both expected 0)\n";
}

// -----------------------------------------------------------------------------
// Test 7: the routine must not disturb the surrounding system state when it
//         bails out (defensive regression check on molecule/complex geometry).
// -----------------------------------------------------------------------------
void test_cir_leaves_system_state_intact()
{
    std::cerr << "\n[TEST] test_cir_leaves_system_state_intact\n"
              << "  Source file:   src/reactions/check_implicit_reactions.cpp\n"
              << "  Function:      check_implicit_reactions()\n"
              << "  Scenario:      a bail-out call (no free lipids) on a fully populated system.\n"
              << "  Pass criteria: coordinates, complex membership and free-lists are all\n"
              << "                 identical before and after the call.\n";

    cir_ensure_rng();

    CirSystem sys = cir_build_system(std::vector<unsigned> { 1u }, 1, /*freeLipids=*/0, true);

    // Snapshot of everything we expect to remain invariant.
    const Coord proComBefore = sys.moleculeList[0].comCoord;
    const Coord proIfaceBefore = sys.moleculeList[0].interfaceList[0].coord;
    const Coord lipidComBefore = sys.moleculeList[1].comCoord;
    const size_t proFreeListBefore = sys.moleculeList[0].freelist.size();
    const size_t lipidFreeListBefore = sys.moleculeList[1].freelist.size();
    const int proComIndexBefore = sys.moleculeList[0].myComIndex;

    std::cerr << "  Calling check_implicit_reactions(pro1=0, pro2=1, simItr=17)...\n";
    check_implicit_reactions(0, 1, 17, sys.params, sys.moleculeList, sys.complexList, sys.molTemplateList,
        sys.forwardRxns, sys.backRxns, sys.counterArrays, sys.membraneObject, sys.IL2DbindingVec,
        sys.IL2DUnbindingVec, sys.ILTableIDs);

    // Coord::operator== compares components after roundv() (4 decimal places).
    EXPECT_TRUE(sys.moleculeList[0].comCoord == proComBefore)
        << "protein centre of mass must not move";
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].coord == proIfaceBefore)
        << "protein interface coordinate must not move";
    EXPECT_TRUE(sys.moleculeList[1].comCoord == lipidComBefore)
        << "implicit-lipid centre of mass must not move";
    EXPECT_EQ(sys.moleculeList[0].freelist.size(), proFreeListBefore)
        << "protein free-list must be unchanged";
    EXPECT_EQ(sys.moleculeList[1].freelist.size(), lipidFreeListBefore)
        << "implicit-lipid free-list must be unchanged";
    EXPECT_EQ(sys.moleculeList[0].myComIndex, proComIndexBefore)
        << "protein complex membership must be unchanged";
    EXPECT_FALSE(sys.moleculeList[0].interfaceList[0].isBound)
        << "no binding may be performed by check_implicit_reactions itself";

    std::cerr << "  Geometry and connectivity verified unchanged.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario is registered separately so that a failure
// in one branch does not hide the results of the others.
// -----------------------------------------------------------------------------
TEST(CheckImplicitReactions, ReturnsWhenPartnerIsNotImplicitLipid)
{
    test_cir_returns_when_partner_is_not_implicit_lipid();
}
TEST(CheckImplicitReactions, NoFreeLipidsBlocksInteraction)
{
    test_cir_no_free_lipids_blocks_interaction();
}
TEST(CheckImplicitReactions, ProteinWithoutFreeInterfaces)
{
    test_cir_protein_without_free_interfaces();
}
TEST(CheckImplicitReactions, StateWithoutReactionPartners)
{
    test_cir_state_without_reaction_partners();
}
TEST(CheckImplicitReactions, PartnerIndexMismatch)
{
    test_cir_partner_index_mismatch();
}
TEST(CheckImplicitReactions, MatchingPartnerWithoutRegisteredReaction)
{
    test_cir_matching_partner_without_registered_reaction();
}
TEST(CheckImplicitReactions, LeavesSystemStateIntact)
{
    test_cir_leaves_system_state_intact();
}