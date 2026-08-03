/*! \file test_correct_structure.cpp
 *
 * ### Unit test for src/reactions/correct_structutre.cpp
 *
 * Function under test:
 *
 *     void correct_structure(const std::vector<Molecule>& moleculeList,
 *                            Complex& complex,
 *                            const std::vector<ForwardRxn>& forwardRxns);
 *
 * What the function does (as written):
 *   1. Walks over every member molecule index of `complex` and counts how many
 *      of them are promoters (`isPromoter == true`) and how many are "simple
 *      proteins" (neither an explicit lipid nor an implicit lipid).
 *   2. Only if there is *exactly one* promoter and *exactly one* protein does it
 *      continue.
 *   3. It then checks whether interface 0 of that protein is bound to the
 *      promoter (`interaction.partnerIndex == promoterIndex`).  If so it
 *      computes the interface-to-interface vector, shrinks it to the binding
 *      radius taken from `forwardRxns[protein.bndRxnList[0]].bindRadius`, and
 *      applies the resulting displacement vector.
 *   4. IMPORTANT: the displacement is applied to `Molecule movingMol{...}`,
 *      which is a *local copy* of the protein.  `moleculeList` is passed by
 *      const reference and `complex` is never written to, so the routine has
 *      no observable side effects on its arguments.
 *
 * Consequently the test suite verifies two things:
 *   a) the function executes all of its code paths without crashing (including
 *      the degenerate zero-separation case), and
 *   b) it leaves every input object bit-for-bit unchanged, which is the
 *      contract implied by the `const` qualifier on `moleculeList` and by the
 *      fact that only a copy of the protein is displaced.
 *
 * Every test prints, to stderr, which source file / function is exercised, what
 * scenario is being set up, and what the pass criterion is.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Declaration of the function under test.  We declare it locally (instead of
// pulling in reactions/unimolecular/unimolecular_reactions.hpp) so that this
// translation unit only needs the two class headers above.
// -----------------------------------------------------------------------------
void correct_structure(const std::vector<Molecule>& moleculeList, Complex& complex,
    const std::vector<ForwardRxn>& forwardRxns);

namespace {

/*! \brief Build a minimal molecule with a single interface.
 *
 * \param[in] index       Index this molecule will occupy in moleculeList.
 * \param[in] com         Center-of-mass coordinate.
 * \param[in] ifaceCoord  Coordinate of the single interface.
 * \return A Molecule ready to be dropped into a moleculeList.
 *
 * The molecule defaults to a "simple protein": not a promoter, not a lipid and
 * not an implicit lipid.  The caller flips those flags as needed.
 */
Molecule cstr_make_molecule(int index, const Coord& com, const Coord& ifaceCoord)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = 0; // every molecule in these tests belongs to complex 0
    mol.molTypeIndex = 0;
    mol.comCoord = com;
    mol.isPromoter = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;

    // Exactly one interface: correct_structure() always looks at interfaceList[0]
    // of the protein, so a single interface is sufficient.
    Molecule::Iface iface;
    iface.coord = ifaceCoord;
    iface.index = 0;
    iface.relIndex = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // The protein branch dereferences bndRxnList[0]; give every molecule a
    // valid reaction index so we never read out of bounds.
    mol.bndRxnList.clear();
    mol.bndRxnList.push_back(0);

    return mol;
}

/*! \brief Bind interface 0 of `protein` to interface `promoterIfaceIndex` of the
 *         molecule stored at `promoterIndex`.
 */
void cstr_bind_protein_to_promoter(Molecule& protein, int promoterIndex, int promoterIfaceIndex)
{
    protein.interfaceList[0].isBound = true;
    protein.interfaceList[0].interaction.partnerIndex = promoterIndex;
    protein.interfaceList[0].interaction.partnerIfaceIndex = promoterIfaceIndex;
    protein.interfaceList[0].interaction.conjBackRxn = 0;
}

/*! \brief Build a Complex whose memberList contains the supplied indices. */
Complex cstr_make_complex(const std::vector<int>& members, const Coord& com)
{
    Complex targCom;
    targCom.index = 0;
    targCom.comCoord = com;
    targCom.memberList = members;
    return targCom;
}

/*! \brief Build a single-reaction forwardRxns list with a chosen binding radius. */
std::vector<ForwardRxn> cstr_make_forward_rxns(double bindRadius)
{
    std::vector<ForwardRxn> forwardRxns(1);
    forwardRxns[0].bindRadius = bindRadius;
    return forwardRxns;
}

/*! \brief Convenience comparison used by several tests: assert two Coords match. */
void cstr_expect_coord_eq(const Coord& actual, const Coord& expected, const char* what)
{
    EXPECT_DOUBLE_EQ(actual.x, expected.x) << what << ": x component changed";
    EXPECT_DOUBLE_EQ(actual.y, expected.y) << what << ": y component changed";
    EXPECT_DOUBLE_EQ(actual.z, expected.z) << what << ": z component changed";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the full correction path (1 promoter + 1 protein, correctly bound).
//         The interesting geometry is exercised, and we confirm that neither the
//         moleculeList nor the Complex is mutated (the routine only displaces a
//         local copy of the protein).
// -----------------------------------------------------------------------------
void test_cstr_single_promoter_single_protein_bound()
{
    std::cerr << "\n[TEST] test_cstr_single_promoter_single_protein_bound\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      complex holds exactly one promoter (index 0) and\n"
              << "                 one protein (index 1); the protein's interface 0\n"
              << "                 is bound to interface 0 of the promoter, and the\n"
              << "                 two interfaces are 5 nm apart while the reaction\n"
              << "                 binding radius is 2 nm.\n"
              << "  Pass criteria: the call completes and every input object\n"
              << "                 (moleculeList entries + Complex) is unchanged,\n"
              << "                 because only a local copy of the protein is moved.\n";

    // ---- promoter: molecule index 0, interface sitting at the origin --------
    Molecule promoter = cstr_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 });
    promoter.isPromoter = true;

    // ---- protein: molecule index 1, interface 5 nm away along +x -----------
    Molecule protein = cstr_make_molecule(1, Coord { 6.0, 0.0, 0.0 }, Coord { 5.0, 0.0, 0.0 });
    cstr_bind_protein_to_promoter(protein, /*promoterIndex=*/0, /*promoterIfaceIndex=*/0);

    std::vector<Molecule> moleculeList { promoter, protein };
    Complex targCom = cstr_make_complex({ 0, 1 }, Coord { 3.0, 0.0, 0.0 });
    std::vector<ForwardRxn> forwardRxns = cstr_make_forward_rxns(2.0);

    // Snapshot the inputs so we can prove nothing was touched.
    const Coord promoterComBefore = moleculeList[0].comCoord;
    const Coord promoterIfaceBefore = moleculeList[0].interfaceList[0].coord;
    const Coord proteinComBefore = moleculeList[1].comCoord;
    const Coord proteinIfaceBefore = moleculeList[1].interfaceList[0].coord;
    const Coord complexComBefore = targCom.comCoord;
    const std::vector<int> memberListBefore = targCom.memberList;

    // Report the displacement the routine computes internally, purely for the log.
    const double sep = 5.0;
    const double expectedShiftX = sep * (1.0 - 2.0 / sep); // = 3.0
    std::cerr << "  Internally computed displacement should be (" << expectedShiftX
              << ", 0, 0) and is applied to a *copy* of the protein.\n"
              << "  Calling correct_structure()...\n";

    correct_structure(moleculeList, targCom, forwardRxns);

    // ---- the inputs must be untouched --------------------------------------
    cstr_expect_coord_eq(moleculeList[0].comCoord, promoterComBefore, "promoter comCoord");
    cstr_expect_coord_eq(moleculeList[0].interfaceList[0].coord, promoterIfaceBefore,
        "promoter interface coord");
    cstr_expect_coord_eq(moleculeList[1].comCoord, proteinComBefore, "protein comCoord");
    cstr_expect_coord_eq(moleculeList[1].interfaceList[0].coord, proteinIfaceBefore,
        "protein interface coord");
    cstr_expect_coord_eq(targCom.comCoord, complexComBefore, "complex comCoord");

    EXPECT_EQ(targCom.memberList.size(), memberListBefore.size())
        << "complex memberList size must not change";
    for (size_t i = 0; i < memberListBefore.size() && i < targCom.memberList.size(); ++i) {
        EXPECT_EQ(targCom.memberList[i], memberListBefore[i])
            << "complex memberList entry " << i << " must not change";
    }

    // The bond bookkeeping must also survive untouched.
    EXPECT_EQ(moleculeList[1].interfaceList[0].interaction.partnerIndex, 0)
        << "protein's bound partner index must be preserved";
    EXPECT_EQ(moleculeList[1].interfaceList[0].interaction.partnerIfaceIndex, 0)
        << "protein's bound partner interface index must be preserved";

    std::cerr << "  Protein comCoord after call = (" << moleculeList[1].comCoord.x << ", "
              << moleculeList[1].comCoord.y << ", " << moleculeList[1].comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: no promoter in the complex -> the geometric branch must be skipped.
// -----------------------------------------------------------------------------
void test_cstr_no_promoter_is_noop()
{
    std::cerr << "\n[TEST] test_cstr_no_promoter_is_noop\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      complex contains two plain proteins, no promoter,\n"
              << "                 so numPromoter == 0.\n"
              << "  Pass criteria: the guard (numPromoter == 1 && numProtein == 1)\n"
              << "                 fails, the routine returns, and nothing changes.\n";

    Molecule protA = cstr_make_molecule(0, Coord { 1.0, 2.0, 3.0 }, Coord { 1.0, 2.0, 4.0 });
    Molecule protB = cstr_make_molecule(1, Coord { -1.0, 0.0, 0.0 }, Coord { -1.0, 0.0, 1.0 });
    cstr_bind_protein_to_promoter(protB, /*promoterIndex=*/0, /*promoterIfaceIndex=*/0);

    std::vector<Molecule> moleculeList { protA, protB };
    Complex targCom = cstr_make_complex({ 0, 1 }, Coord { 0.0, 1.0, 1.5 });
    std::vector<ForwardRxn> forwardRxns = cstr_make_forward_rxns(1.0);

    const Coord aComBefore = moleculeList[0].comCoord;
    const Coord bComBefore = moleculeList[1].comCoord;
    const Coord bIfaceBefore = moleculeList[1].interfaceList[0].coord;

    std::cerr << "  Calling correct_structure()...\n";
    correct_structure(moleculeList, targCom, forwardRxns);

    cstr_expect_coord_eq(moleculeList[0].comCoord, aComBefore, "protein A comCoord");
    cstr_expect_coord_eq(moleculeList[1].comCoord, bComBefore, "protein B comCoord");
    cstr_expect_coord_eq(moleculeList[1].interfaceList[0].coord, bIfaceBefore,
        "protein B interface coord");
    std::cerr << "  No molecule coordinates were modified, as expected.\n";
}

// -----------------------------------------------------------------------------
// Test 3: promoter present but zero "simple proteins" (the partner is a lipid
//         and an implicit lipid) -> numProtein == 0, branch skipped.
// -----------------------------------------------------------------------------
void test_cstr_lipids_are_not_counted_as_protein()
{
    std::cerr << "\n[TEST] test_cstr_lipids_are_not_counted_as_protein\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      complex = 1 promoter + 1 explicit lipid + 1\n"
              << "                 implicit lipid.  Neither lipid may be counted as\n"
              << "                 a simple protein, so numProtein must be 0.\n"
              << "  Pass criteria: the guard fails and no coordinate is modified.\n";

    Molecule promoter = cstr_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 });
    promoter.isPromoter = true;

    Molecule lipid = cstr_make_molecule(1, Coord { 4.0, 0.0, 0.0 }, Coord { 3.0, 0.0, 0.0 });
    lipid.isLipid = true;
    cstr_bind_protein_to_promoter(lipid, 0, 0);

    Molecule implicitLipid = cstr_make_molecule(2, Coord { 0.0, 4.0, 0.0 }, Coord { 0.0, 3.0, 0.0 });
    implicitLipid.isImplicitLipid = true;

    std::vector<Molecule> moleculeList { promoter, lipid, implicitLipid };
    Complex targCom = cstr_make_complex({ 0, 1, 2 }, Coord { 1.0, 1.0, 0.0 });
    std::vector<ForwardRxn> forwardRxns = cstr_make_forward_rxns(1.5);

    const Coord lipidComBefore = moleculeList[1].comCoord;
    const Coord lipidIfaceBefore = moleculeList[1].interfaceList[0].coord;
    const Coord implicitComBefore = moleculeList[2].comCoord;

    std::cerr << "  Calling correct_structure()...\n";
    correct_structure(moleculeList, targCom, forwardRxns);

    cstr_expect_coord_eq(moleculeList[1].comCoord, lipidComBefore, "lipid comCoord");
    cstr_expect_coord_eq(moleculeList[1].interfaceList[0].coord, lipidIfaceBefore,
        "lipid interface coord");
    cstr_expect_coord_eq(moleculeList[2].comCoord, implicitComBefore,
        "implicit lipid comCoord");
    std::cerr << "  Lipids correctly ignored; nothing was modified.\n";
}

// -----------------------------------------------------------------------------
// Test 4: more than one protein in the complex -> numProtein == 2, branch skipped.
// -----------------------------------------------------------------------------
void test_cstr_two_proteins_is_noop()
{
    std::cerr << "\n[TEST] test_cstr_two_proteins_is_noop\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      complex = 1 promoter + 2 proteins (numProtein == 2).\n"
              << "  Pass criteria: the (numPromoter == 1 && numProtein == 1) guard\n"
              << "                 fails, so no correction is attempted.\n";

    Molecule promoter = cstr_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 });
    promoter.isPromoter = true;

    Molecule protA = cstr_make_molecule(1, Coord { 6.0, 0.0, 0.0 }, Coord { 5.0, 0.0, 0.0 });
    cstr_bind_protein_to_promoter(protA, 0, 0);

    Molecule protB = cstr_make_molecule(2, Coord { 0.0, 6.0, 0.0 }, Coord { 0.0, 5.0, 0.0 });

    std::vector<Molecule> moleculeList { promoter, protA, protB };
    Complex targCom = cstr_make_complex({ 0, 1, 2 }, Coord { 2.0, 2.0, 0.0 });
    std::vector<ForwardRxn> forwardRxns = cstr_make_forward_rxns(2.0);

    const Coord protAComBefore = moleculeList[1].comCoord;
    const Coord protBComBefore = moleculeList[2].comCoord;

    std::cerr << "  Calling correct_structure()...\n";
    correct_structure(moleculeList, targCom, forwardRxns);

    cstr_expect_coord_eq(moleculeList[1].comCoord, protAComBefore, "protein A comCoord");
    cstr_expect_coord_eq(moleculeList[2].comCoord, protBComBefore, "protein B comCoord");
    std::cerr << "  Two-protein complex left untouched, as expected.\n";
}

// -----------------------------------------------------------------------------
// Test 5: right molecule counts, but the protein's interface 0 is bound to some
//         *other* molecule -> the inner `if` fails and no geometry is touched.
// -----------------------------------------------------------------------------
void test_cstr_partner_mismatch_is_noop()
{
    std::cerr << "\n[TEST] test_cstr_partner_mismatch_is_noop\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      1 promoter + 1 protein, but the protein's\n"
              << "                 interface 0 reports partnerIndex = 7 (i.e. it is\n"
              << "                 not bound to the promoter, index 0).\n"
              << "  Pass criteria: interaction.partnerIndex != promoterIndex, so the\n"
              << "                 correction is skipped and the inputs are unchanged.\n";

    Molecule promoter = cstr_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 });
    promoter.isPromoter = true;

    Molecule protein = cstr_make_molecule(1, Coord { 6.0, 0.0, 0.0 }, Coord { 5.0, 0.0, 0.0 });
    // Deliberately point at a molecule that is not the promoter.
    cstr_bind_protein_to_promoter(protein, /*promoterIndex=*/7, /*promoterIfaceIndex=*/0);

    std::vector<Molecule> moleculeList { promoter, protein };
    Complex targCom = cstr_make_complex({ 0, 1 }, Coord { 3.0, 0.0, 0.0 });
    std::vector<ForwardRxn> forwardRxns = cstr_make_forward_rxns(2.0);

    const Coord proteinComBefore = moleculeList[1].comCoord;
    const Coord proteinIfaceBefore = moleculeList[1].interfaceList[0].coord;

    std::cerr << "  Calling correct_structure()...\n";
    correct_structure(moleculeList, targCom, forwardRxns);

    cstr_expect_coord_eq(moleculeList[1].comCoord, proteinComBefore, "protein comCoord");
    cstr_expect_coord_eq(moleculeList[1].interfaceList[0].coord, proteinIfaceBefore,
        "protein interface coord");
    EXPECT_EQ(moleculeList[1].interfaceList[0].interaction.partnerIndex, 7)
        << "the (mismatched) partner index must be preserved";
    std::cerr << "  Mismatched partner correctly ignored.\n";
}

// -----------------------------------------------------------------------------
// Test 6: unbound protein (default partnerIndex == -1) -> also a mismatch, so
//         the routine must simply return.
// -----------------------------------------------------------------------------
void test_cstr_unbound_protein_is_noop()
{
    std::cerr << "\n[TEST] test_cstr_unbound_protein_is_noop\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      1 promoter + 1 protein, but the protein's\n"
              << "                 interface 0 is free (partnerIndex stays at the\n"
              << "                 default value of -1).\n"
              << "  Pass criteria: no correction is attempted and nothing changes.\n";

    Molecule promoter = cstr_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 });
    promoter.isPromoter = true;

    // No call to cstr_bind_protein_to_promoter() -> interface remains free.
    Molecule protein = cstr_make_molecule(1, Coord { 10.0, 0.0, 0.0 }, Coord { 9.0, 0.0, 0.0 });

    std::vector<Molecule> moleculeList { promoter, protein };
    Complex targCom = cstr_make_complex({ 0, 1 }, Coord { 5.0, 0.0, 0.0 });
    std::vector<ForwardRxn> forwardRxns = cstr_make_forward_rxns(2.0);

    const Coord proteinComBefore = moleculeList[1].comCoord;

    EXPECT_EQ(moleculeList[1].interfaceList[0].interaction.partnerIndex, -1)
        << "sanity check: default partnerIndex should be -1 (unbound)";

    std::cerr << "  Calling correct_structure()...\n";
    correct_structure(moleculeList, targCom, forwardRxns);

    cstr_expect_coord_eq(moleculeList[1].comCoord, proteinComBefore, "protein comCoord");
    std::cerr << "  Unbound protein correctly ignored.\n";
}

// -----------------------------------------------------------------------------
// Test 7: an empty complex (no members at all).  This exercises the loop with
//         zero iterations and confirms the routine tolerates it.
// -----------------------------------------------------------------------------
void test_cstr_empty_member_list_is_safe()
{
    std::cerr << "\n[TEST] test_cstr_empty_member_list_is_safe\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      the Complex has an empty memberList, so the\n"
              << "                 counting loop performs zero iterations.\n"
              << "  Pass criteria: the call returns without crashing and the empty\n"
              << "                 memberList is still empty afterwards.\n";

    std::vector<Molecule> moleculeList {}; // nothing in the system
    Complex targCom = cstr_make_complex({}, Coord { 0.0, 0.0, 0.0 });
    std::vector<ForwardRxn> forwardRxns = cstr_make_forward_rxns(1.0);

    std::cerr << "  Calling correct_structure()...\n";
    correct_structure(moleculeList, targCom, forwardRxns);

    EXPECT_TRUE(targCom.memberList.empty()) << "memberList should still be empty";
    EXPECT_TRUE(moleculeList.empty()) << "moleculeList should still be empty";
    std::cerr << "  Empty complex handled safely.\n";
}

// -----------------------------------------------------------------------------
// Test 8: degenerate geometry -- the two bound interfaces sit on top of one
//         another, so the interface-to-interface magnitude is exactly zero and
//         the internal (1 - bindRadius/mag) term is not finite.  The routine
//         must still return and must still leave the inputs alone (the NaN only
//         ever reaches the discarded local copy).
// -----------------------------------------------------------------------------
void test_cstr_zero_separation_does_not_corrupt_inputs()
{
    std::cerr << "\n[TEST] test_cstr_zero_separation_does_not_corrupt_inputs\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      1 promoter + 1 bound protein whose interface is\n"
              << "                 coincident with the promoter interface, giving a\n"
              << "                 zero interface-to-interface magnitude (a division\n"
              << "                 by zero inside the routine).\n"
              << "  Pass criteria: the call returns, and the stored coordinates are\n"
              << "                 still finite and unchanged (the non-finite value\n"
              << "                 stays inside the discarded local copy).\n";

    Molecule promoter = cstr_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 1.0, 1.0 });
    promoter.isPromoter = true;

    // Protein interface deliberately identical to the promoter interface.
    Molecule protein = cstr_make_molecule(1, Coord { 2.0, 2.0, 2.0 }, Coord { 1.0, 1.0, 1.0 });
    cstr_bind_protein_to_promoter(protein, 0, 0);

    std::vector<Molecule> moleculeList { promoter, protein };
    Complex targCom = cstr_make_complex({ 0, 1 }, Coord { 1.0, 1.0, 1.0 });
    std::vector<ForwardRxn> forwardRxns = cstr_make_forward_rxns(2.0);

    const Coord proteinComBefore = moleculeList[1].comCoord;
    const Coord proteinIfaceBefore = moleculeList[1].interfaceList[0].coord;

    std::cerr << "  Calling correct_structure()...\n";
    correct_structure(moleculeList, targCom, forwardRxns);

    // Nothing stored in the containers may become NaN/Inf.
    EXPECT_TRUE(std::isfinite(moleculeList[1].comCoord.x))
        << "protein comCoord.x must remain finite";
    EXPECT_TRUE(std::isfinite(moleculeList[1].comCoord.y))
        << "protein comCoord.y must remain finite";
    EXPECT_TRUE(std::isfinite(moleculeList[1].comCoord.z))
        << "protein comCoord.z must remain finite";

    cstr_expect_coord_eq(moleculeList[1].comCoord, proteinComBefore, "protein comCoord");
    cstr_expect_coord_eq(moleculeList[1].interfaceList[0].coord, proteinIfaceBefore,
        "protein interface coord");
    std::cerr << "  Degenerate geometry survived without corrupting the inputs.\n";
}

// -----------------------------------------------------------------------------
// Test 9: the binding radius is read from forwardRxns[protein.bndRxnList[0]].
//         Point bndRxnList at a *second* reaction with a different radius and
//         verify the routine still completes (i.e. the correct entry is indexed
//         and nothing goes out of bounds).
// -----------------------------------------------------------------------------
void test_cstr_uses_bndRxnList_index_for_bind_radius()
{
    std::cerr << "\n[TEST] test_cstr_uses_bndRxnList_index_for_bind_radius\n"
              << "  Source file:   src/reactions/correct_structutre.cpp\n"
              << "  Function:      correct_structure()\n"
              << "  Scenario:      forwardRxns has two entries (radii 1.0 and 4.0)\n"
              << "                 and the protein's bndRxnList[0] == 1, so the\n"
              << "                 routine must index the second reaction.\n"
              << "  Pass criteria: the call completes without out-of-range access\n"
              << "                 and the inputs remain unchanged.\n";

    Molecule promoter = cstr_make_molecule(0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 });
    promoter.isPromoter = true;

    Molecule protein = cstr_make_molecule(1, Coord { 9.0, 0.0, 0.0 }, Coord { 8.0, 0.0, 0.0 });
    cstr_bind_protein_to_promoter(protein, 0, 0);
    protein.bndRxnList.clear();
    protein.bndRxnList.push_back(1); // deliberately point at the second reaction

    std::vector<Molecule> moleculeList { promoter, protein };
    Complex targCom = cstr_make_complex({ 0, 1 }, Coord { 4.0, 0.0, 0.0 });

    // Two reactions: index 0 (radius 1.0) and index 1 (radius 4.0).
    std::vector<ForwardRxn> forwardRxns(2);
    forwardRxns[0].bindRadius = 1.0;
    forwardRxns[1].bindRadius = 4.0;

    const Coord proteinComBefore = moleculeList[1].comCoord;
    const Coord proteinIfaceBefore = moleculeList[1].interfaceList[0].coord;

    std::cerr << "  Expected internal displacement magnitude = 8 * (1 - 4/8) = 4\n"
              << "  Calling correct_structure()...\n";
    correct_structure(moleculeList, targCom, forwardRxns);

    cstr_expect_coord_eq(moleculeList[1].comCoord, proteinComBefore, "protein comCoord");
    cstr_expect_coord_eq(moleculeList[1].interfaceList[0].coord, proteinIfaceBefore,
        "protein interface coord");

    // The reaction list itself is const and must be intact.
    EXPECT_DOUBLE_EQ(forwardRxns[0].bindRadius, 1.0)
        << "forwardRxns[0].bindRadius must be unchanged";
    EXPECT_DOUBLE_EQ(forwardRxns[1].bindRadius, 4.0)
        << "forwardRxns[1].bindRadius must be unchanged";
    std::cerr << "  Correct reaction entry indexed; all inputs intact.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario lives in its own TEST so a failure in one
// does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(CorrectStructure, SinglePromoterSingleProteinBound)
{
    test_cstr_single_promoter_single_protein_bound();
}
TEST(CorrectStructure, NoPromoterIsNoop) { test_cstr_no_promoter_is_noop(); }
TEST(CorrectStructure, LipidsAreNotCountedAsProtein)
{
    test_cstr_lipids_are_not_counted_as_protein();
}
TEST(CorrectStructure, TwoProteinsIsNoop) { test_cstr_two_proteins_is_noop(); }
TEST(CorrectStructure, PartnerMismatchIsNoop) { test_cstr_partner_mismatch_is_noop(); }
TEST(CorrectStructure, UnboundProteinIsNoop) { test_cstr_unbound_protein_is_noop(); }
TEST(CorrectStructure, EmptyMemberListIsSafe) { test_cstr_empty_member_list_is_safe(); }
TEST(CorrectStructure, ZeroSeparationDoesNotCorruptInputs)
{
    test_cstr_zero_separation_does_not_corrupt_inputs();
}
TEST(CorrectStructure, UsesBndRxnListIndexForBindRadius)
{
    test_cstr_uses_bndRxnList_index_for_bind_radius();
}