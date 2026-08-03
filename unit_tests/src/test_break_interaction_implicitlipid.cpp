/*! \file test_break_interaction_implicitlipid.cpp
 *
 * ### Unit test for src/reactions/break_interaction_implicitlipid.cpp
 *
 * Function under test:
 *
 *     void break_interaction_implicitlipid(long long int iter,
 *                                          size_t relIface1, size_t relIface2,
 *                                          Molecule& reactMol1, Molecule& reactMol2,
 *                                          const BackRxn& currRxn,
 *                                          std::vector<Molecule>& moleculeList,
 *                                          std::vector<Complex>& complexList,
 *                                          std::vector<MolTemplate>& molTemplateList,
 *                                          std::ofstream& assocDissocFile)
 *
 * The routine breaks the bond bookkeeping for reactMol1 only (reactMol2 is the
 * implicit lipid and is left alone).  Concretely it:
 *
 *   1. Determines the new absolute interface index for reactMol1 from the
 *      BackRxn product list.  If the reaction is symmetric, product[0] is used
 *      for mol1 and product[1] for mol2.  Otherwise it matches on
 *      (molTypeIndex, relIfaceIndex, requiresState) against product[0], then
 *      product[1] (in which case the two absolute indices are swapped).
 *   2. Clears the Interaction on reactMol1's reacting interface, marks it
 *      unbound, and assigns the new absolute index.
 *   3. Optionally logs a "BREAK" record to the assoc/dissoc file.
 *   4. Moves the interface from bndlist to freelist and erases reactMol2 from
 *      bndpartner.
 *   5. If reactMol1 has become a monomer and its MolTemplate canDestroy, adds
 *      it to that template's monomerList.
 *
 * Notes on what is *not* tested: the "no product matched" branch calls exit(1),
 * which cannot be exercised without killing the whole test binary, so every
 * scenario below is constructed to match a product.
 */

#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Small local helpers used to build the minimal objects the function needs.
// All helpers are file-local (anonymous namespace) so they cannot collide with
// other translation units in the suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a molecule with \p numIfaces interfaces, all in state \p state.
 *
 * \param[in] index        Index of the molecule in moleculeList.
 * \param[in] molTypeIndex Index of the molecule's MolTemplate.
 * \param[in] numIfaces    How many interfaces to create.
 * \param[in] state        State identifier stamped onto every interface.
 */
Molecule bii_make_molecule(int index, int molTypeIndex, int numIfaces, char state)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = index; // one molecule per complex; not used by the function
    mol.interfaceList.resize(static_cast<size_t>(numIfaces));
    for (int i = 0; i < numIfaces; ++i) {
        mol.interfaceList[i].relIndex = i;
        mol.interfaceList[i].molTypeIndex = molTypeIndex;
        mol.interfaceList[i].stateIden = state;
        mol.interfaceList[i].index = -1; // absolute index assigned by the function
        mol.interfaceList[i].isBound = false;
    }
    return mol;
}

/*! \brief Mark interface \p relIface of \p mol as bound to \p partnerIndex.
 *
 * Also adds the interface to bndlist and the partner to bndpartner, both of
 * which the function under test erases from (erasing an unfound element would
 * be undefined behaviour, so the setup must be consistent).
 */
void bii_bind_interface(Molecule& mol, size_t relIface, int partnerIndex, int partnerIface, int conjBackRxn)
{
    mol.interfaceList[relIface].isBound = true;
    mol.interfaceList[relIface].interaction.partnerIndex = partnerIndex;
    mol.interfaceList[relIface].interaction.partnerIfaceIndex = partnerIface;
    mol.interfaceList[relIface].interaction.conjBackRxn = conjBackRxn;
    mol.bndlist.push_back(static_cast<int>(relIface));
    mol.bndpartner.push_back(partnerIndex);
}

/*! \brief Convenience wrapper for constructing a RxnIface product entry. */
RxnIface bii_make_rxn_iface(
    const std::string& ifaceName, int molTypeIndex, int absIfaceIndex, int relIfaceIndex, char requiresState)
{
    return RxnIface(ifaceName, molTypeIndex, absIfaceIndex, relIfaceIndex, requiresState, false);
}

/*! \brief Build a MolTemplate with just the fields the function reads. */
MolTemplate bii_make_template(const std::string& name, int molTypeIndex, bool canDestroy)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = molTypeIndex;
    temp.canDestroy = canDestroy;
    temp.monomerList.clear();
    return temp;
}

/*! \brief Returns true if \p vec contains \p val. */
bool bii_contains(const std::vector<int>& vec, int val)
{
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: symmetric BackRxn.
//   absIface1 <- productListNew[0].absIfaceIndex
//   absIface2 <- productListNew[1].absIfaceIndex (unused for mol1, but the
//                branch must still be taken without any state matching).
// -----------------------------------------------------------------------------
void test_bii_symmetric_reaction_updates_mol1()
{
    std::cerr << "\n[TEST] test_bii_symmetric_reaction_updates_mol1\n"
              << "  Source file: break_interaction_implicitlipid.cpp\n"
              << "  Function:    break_interaction_implicitlipid (isSymmetric == true)\n"
              << "  Scenario:    A(a) bound to implicit lipid; symmetric back reaction.\n"
              << "  Pass criteria: mol1's reacting interface becomes unbound, gets the\n"
              << "                 absolute index of productListNew[0], is moved from\n"
              << "                 bndlist to freelist, and the partner is erased from\n"
              << "                 bndpartner. mol2 must be untouched.\n";

    // --- templates -----------------------------------------------------------
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(bii_make_template("A", 0, /*canDestroy=*/false));
    molTemplateList.push_back(bii_make_template("IL", 1, /*canDestroy=*/false));

    // --- molecules -----------------------------------------------------------
    // mol1 is a protein with two interfaces; interface 0 is the bound one.
    Molecule mol1 = bii_make_molecule(/*index=*/0, /*molTypeIndex=*/0, /*numIfaces=*/2, '\0');
    mol1.freelist.push_back(1); // interface 1 is already free
    bii_bind_interface(mol1, /*relIface=*/0, /*partnerIndex=*/1, /*partnerIface=*/0, /*conjBackRxn=*/3);

    // mol2 stands in for the implicit lipid.
    Molecule mol2 = bii_make_molecule(/*index=*/1, /*molTypeIndex=*/1, /*numIfaces=*/1, '\0');
    bii_bind_interface(mol2, /*relIface=*/0, /*partnerIndex=*/0, /*partnerIface=*/0, /*conjBackRxn=*/3);

    std::vector<Molecule> moleculeList { mol1, mol2 };
    std::vector<Complex> complexList {}; // unused by the function

    // --- reaction ------------------------------------------------------------
    BackRxn currRxn;
    currRxn.isSymmetric = true;
    currRxn.productListNew.push_back(bii_make_rxn_iface("a", 0, /*abs=*/5, /*rel=*/0, '\0'));
    currRxn.productListNew.push_back(bii_make_rxn_iface("il", 1, /*abs=*/6, /*rel=*/0, '\0'));

    // --- call ----------------------------------------------------------------
    std::ofstream closedFile; // deliberately not opened -> no logging path
    std::cerr << "  Calling break_interaction_implicitlipid with a closed log stream...\n";
    break_interaction_implicitlipid(/*iter=*/42, /*relIface1=*/0, /*relIface2=*/0, moleculeList[0], moleculeList[1],
        currRxn, moleculeList, complexList, molTemplateList, closedFile);

    const Molecule& out1 = moleculeList[0];
    const Molecule& out2 = moleculeList[1];

    // The interface must now report the product's absolute index.
    std::cerr << "  mol1.interfaceList[0].index = " << out1.interfaceList[0].index << " (expected 5)\n";
    EXPECT_EQ(out1.interfaceList[0].index, 5) << "absIface1 should come from productListNew[0]";

    // It must be flagged unbound.
    EXPECT_FALSE(out1.interfaceList[0].isBound) << "Broken interface should no longer be bound";

    // The Interaction must have been cleared (implementation resets the indices;
    // accept either the 0 or -1 convention, but it must no longer point at mol2).
    std::cerr << "  cleared interaction.partnerIndex = " << out1.interfaceList[0].interaction.partnerIndex << '\n';
    EXPECT_LE(out1.interfaceList[0].interaction.partnerIndex, 0)
        << "Interaction should have been cleared (no longer pointing at partner index 1)";

    // Bookkeeping lists.
    EXPECT_TRUE(bii_contains(out1.freelist, 0)) << "Broken interface 0 should be pushed onto freelist";
    EXPECT_EQ(out1.freelist.size(), 2u) << "freelist should now hold both interfaces (1 pre-existing + 1 new)";
    EXPECT_TRUE(out1.bndlist.empty()) << "bndlist should be empty after erasing interface 0";
    EXPECT_TRUE(out1.bndpartner.empty()) << "bndpartner should be empty after erasing partner 1";

    // mol2 (the implicit lipid) is not modified by this function.
    EXPECT_TRUE(out2.interfaceList[0].isBound) << "reactMol2 must be left untouched by this function";
    EXPECT_EQ(out2.bndpartner.size(), 1u) << "reactMol2's bndpartner must be left untouched";
}

// -----------------------------------------------------------------------------
// Test 2: non-symmetric BackRxn where mol1 matches productListNew[0].
// -----------------------------------------------------------------------------
void test_bii_nonsymmetric_matches_first_product()
{
    std::cerr << "\n[TEST] test_bii_nonsymmetric_matches_first_product\n"
              << "  Source file: break_interaction_implicitlipid.cpp\n"
              << "  Function:    break_interaction_implicitlipid (isSymmetric == false)\n"
              << "  Scenario:    mol1's (molTypeIndex, relIface, state) matches product[0].\n"
              << "  Pass criteria: mol1's interface receives product[0]'s absIfaceIndex (11).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(bii_make_template("A", 0, /*canDestroy=*/false));
    molTemplateList.push_back(bii_make_template("IL", 1, /*canDestroy=*/false));

    // mol1's interface is in state 'U'; the product requires state 'U'.
    Molecule mol1 = bii_make_molecule(0, 0, 1, 'U');
    bii_bind_interface(mol1, 0, /*partnerIndex=*/1, 0, 3);
    Molecule mol2 = bii_make_molecule(1, 1, 1, '\0');
    bii_bind_interface(mol2, 0, /*partnerIndex=*/0, 0, 3);

    std::vector<Molecule> moleculeList { mol1, mol2 };
    std::vector<Complex> complexList {};

    BackRxn currRxn;
    currRxn.isSymmetric = false;
    // product[0] describes mol1 (type 0, rel iface 0, state 'U').
    currRxn.productListNew.push_back(bii_make_rxn_iface("a", 0, /*abs=*/11, /*rel=*/0, 'U'));
    // product[1] describes the implicit lipid.
    currRxn.productListNew.push_back(bii_make_rxn_iface("il", 1, /*abs=*/12, /*rel=*/0, '\0'));

    std::ofstream closedFile;
    std::cerr << "  Calling break_interaction_implicitlipid...\n";
    break_interaction_implicitlipid(7, 0, 0, moleculeList[0], moleculeList[1], currRxn, moleculeList, complexList,
        molTemplateList, closedFile);

    std::cerr << "  mol1.interfaceList[0].index = " << moleculeList[0].interfaceList[0].index << " (expected 11)\n";
    EXPECT_EQ(moleculeList[0].interfaceList[0].index, 11)
        << "When mol1 matches product[0], absIface1 must be product[0].absIfaceIndex";
    EXPECT_FALSE(moleculeList[0].interfaceList[0].isBound) << "Interface should be unbound after the break";
    EXPECT_TRUE(bii_contains(moleculeList[0].freelist, 0)) << "Interface 0 should now be free";
}

// -----------------------------------------------------------------------------
// Test 3: non-symmetric BackRxn where mol1 matches productListNew[1].
//         In that branch the implementation swaps the two absolute indices, so
//         mol1 must receive product[1]'s absIfaceIndex.
// -----------------------------------------------------------------------------
void test_bii_nonsymmetric_matches_second_product()
{
    std::cerr << "\n[TEST] test_bii_nonsymmetric_matches_second_product\n"
              << "  Source file: break_interaction_implicitlipid.cpp\n"
              << "  Function:    break_interaction_implicitlipid (isSymmetric == false)\n"
              << "  Scenario:    the product list is listed in the reverse order, so mol1\n"
              << "               matches product[1] instead of product[0].\n"
              << "  Pass criteria: mol1's interface receives product[1]'s absIfaceIndex (22),\n"
              << "                 demonstrating the index swap in that branch.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(bii_make_template("A", 0, /*canDestroy=*/false));
    molTemplateList.push_back(bii_make_template("IL", 1, /*canDestroy=*/false));

    Molecule mol1 = bii_make_molecule(0, /*molTypeIndex=*/0, 1, 'P');
    bii_bind_interface(mol1, 0, /*partnerIndex=*/1, 0, 3);
    Molecule mol2 = bii_make_molecule(1, /*molTypeIndex=*/1, 1, '\0');
    bii_bind_interface(mol2, 0, /*partnerIndex=*/0, 0, 3);

    std::vector<Molecule> moleculeList { mol1, mol2 };
    std::vector<Complex> complexList {};

    BackRxn currRxn;
    currRxn.isSymmetric = false;
    // product[0] now describes the implicit lipid (different molTypeIndex).
    currRxn.productListNew.push_back(bii_make_rxn_iface("il", 1, /*abs=*/21, /*rel=*/0, '\0'));
    // product[1] describes mol1 (type 0, rel iface 0, state 'P').
    currRxn.productListNew.push_back(bii_make_rxn_iface("a", 0, /*abs=*/22, /*rel=*/0, 'P'));

    std::ofstream closedFile;
    std::cerr << "  Calling break_interaction_implicitlipid...\n";
    break_interaction_implicitlipid(9, 0, 0, moleculeList[0], moleculeList[1], currRxn, moleculeList, complexList,
        molTemplateList, closedFile);

    std::cerr << "  mol1.interfaceList[0].index = " << moleculeList[0].interfaceList[0].index << " (expected 22)\n";
    EXPECT_EQ(moleculeList[0].interfaceList[0].index, 22)
        << "When mol1 matches product[1], absIface1 must be product[1].absIfaceIndex";
    EXPECT_TRUE(moleculeList[0].bndlist.empty()) << "bndlist should be empty after the only bond is broken";
    EXPECT_TRUE(moleculeList[0].bndpartner.empty()) << "bndpartner should be empty after the only bond is broken";
}

// -----------------------------------------------------------------------------
// Test 4: monomerList maintenance.
//   (a) canDestroy == true and mol1 becomes a monomer  -> index appended.
//   (b) canDestroy == true but mol1 keeps a partner    -> nothing appended.
//   (c) canDestroy == false and mol1 becomes a monomer -> nothing appended.
// -----------------------------------------------------------------------------
void test_bii_monomer_list_update()
{
    std::cerr << "\n[TEST] test_bii_monomer_list_update\n"
              << "  Source file: break_interaction_implicitlipid.cpp\n"
              << "  Function:    break_interaction_implicitlipid (monomerList block)\n"
              << "  Pass criteria: the molecule index is appended to its MolTemplate's\n"
              << "                 monomerList only when the molecule has become a\n"
              << "                 monomer AND its template canDestroy.\n";

    // ---- case (a): becomes a monomer, canDestroy == true --------------------
    {
        std::cerr << "  (a) canDestroy=true and mol1 becomes a monomer -> expect append\n";
        std::vector<MolTemplate> molTemplateList;
        molTemplateList.push_back(bii_make_template("A", 0, /*canDestroy=*/true));
        molTemplateList.push_back(bii_make_template("IL", 1, /*canDestroy=*/false));

        Molecule mol1 = bii_make_molecule(/*index=*/4, 0, 1, '\0');
        bii_bind_interface(mol1, 0, /*partnerIndex=*/1, 0, 3);
        Molecule mol2 = bii_make_molecule(1, 1, 1, '\0');
        bii_bind_interface(mol2, 0, /*partnerIndex=*/4, 0, 3);

        std::vector<Molecule> moleculeList { mol1, mol2 };
        std::vector<Complex> complexList {};

        BackRxn currRxn;
        currRxn.isSymmetric = true;
        currRxn.productListNew.push_back(bii_make_rxn_iface("a", 0, 30, 0, '\0'));
        currRxn.productListNew.push_back(bii_make_rxn_iface("il", 1, 31, 0, '\0'));

        std::ofstream closedFile;
        break_interaction_implicitlipid(1, 0, 0, moleculeList[0], moleculeList[1], currRxn, moleculeList, complexList,
            molTemplateList, closedFile);

        std::cerr << "      monomerList size = " << molTemplateList[0].monomerList.size() << " (expected 1)\n";
        EXPECT_EQ(molTemplateList[0].monomerList.size(), 1u)
            << "A newly produced destroyable monomer should be recorded once";
        if (!molTemplateList[0].monomerList.empty()) {
            EXPECT_EQ(molTemplateList[0].monomerList[0], 4)
                << "The recorded monomer should be the molecule's index (4)";
        }
    }

    // ---- case (b): still has a partner, canDestroy == true ------------------
    {
        std::cerr << "  (b) canDestroy=true but mol1 keeps another partner -> expect no append\n";
        std::vector<MolTemplate> molTemplateList;
        molTemplateList.push_back(bii_make_template("A", 0, /*canDestroy=*/true));
        molTemplateList.push_back(bii_make_template("IL", 1, /*canDestroy=*/false));

        // mol1 has two bonds: iface 0 -> mol index 1 (the one we break),
        //                     iface 1 -> mol index 2 (kept).
        Molecule mol1 = bii_make_molecule(/*index=*/0, 0, 2, '\0');
        bii_bind_interface(mol1, 0, /*partnerIndex=*/1, 0, 3);
        bii_bind_interface(mol1, 1, /*partnerIndex=*/2, 0, 3);
        Molecule mol2 = bii_make_molecule(1, 1, 1, '\0');
        bii_bind_interface(mol2, 0, 0, 0, 3);

        std::vector<Molecule> moleculeList { mol1, mol2 };
        std::vector<Complex> complexList {};

        BackRxn currRxn;
        currRxn.isSymmetric = true;
        currRxn.productListNew.push_back(bii_make_rxn_iface("a", 0, 40, 0, '\0'));
        currRxn.productListNew.push_back(bii_make_rxn_iface("il", 1, 41, 0, '\0'));

        std::ofstream closedFile;
        break_interaction_implicitlipid(2, 0, 0, moleculeList[0], moleculeList[1], currRxn, moleculeList, complexList,
            molTemplateList, closedFile);

        std::cerr << "      remaining bndpartner size = " << moleculeList[0].bndpartner.size()
                  << ", monomerList size = " << molTemplateList[0].monomerList.size() << " (expected 1 and 0)\n";
        EXPECT_EQ(moleculeList[0].bndpartner.size(), 1u) << "Only the broken partner should have been erased";
        EXPECT_TRUE(molTemplateList[0].monomerList.empty())
            << "A molecule that still has partners is not a monomer and must not be recorded";
    }

    // ---- case (c): becomes a monomer, canDestroy == false -------------------
    {
        std::cerr << "  (c) canDestroy=false and mol1 becomes a monomer -> expect no append\n";
        std::vector<MolTemplate> molTemplateList;
        molTemplateList.push_back(bii_make_template("A", 0, /*canDestroy=*/false));
        molTemplateList.push_back(bii_make_template("IL", 1, /*canDestroy=*/false));

        Molecule mol1 = bii_make_molecule(0, 0, 1, '\0');
        bii_bind_interface(mol1, 0, 1, 0, 3);
        Molecule mol2 = bii_make_molecule(1, 1, 1, '\0');
        bii_bind_interface(mol2, 0, 0, 0, 3);

        std::vector<Molecule> moleculeList { mol1, mol2 };
        std::vector<Complex> complexList {};

        BackRxn currRxn;
        currRxn.isSymmetric = true;
        currRxn.productListNew.push_back(bii_make_rxn_iface("a", 0, 50, 0, '\0'));
        currRxn.productListNew.push_back(bii_make_rxn_iface("il", 1, 51, 0, '\0'));

        std::ofstream closedFile;
        break_interaction_implicitlipid(3, 0, 0, moleculeList[0], moleculeList[1], currRxn, moleculeList, complexList,
            molTemplateList, closedFile);

        std::cerr << "      monomerList size = " << molTemplateList[0].monomerList.size() << " (expected 0)\n";
        EXPECT_TRUE(moleculeList[0].bndpartner.empty()) << "Molecule should be a monomer after the break";
        EXPECT_TRUE(molTemplateList[0].monomerList.empty())
            << "Templates that cannot be destroyed must not accumulate monomers";
    }
}

// -----------------------------------------------------------------------------
// Test 5: assoc/dissoc logging.  When the supplied ofstream is open, a single
//         "BREAK" record containing the iteration, both molecule names/indices
//         and both relative interface indices must be written.
// -----------------------------------------------------------------------------
void test_bii_writes_assoc_dissoc_record()
{
    std::cerr << "\n[TEST] test_bii_writes_assoc_dissoc_record\n"
              << "  Source file: break_interaction_implicitlipid.cpp\n"
              << "  Function:    break_interaction_implicitlipid (logging block)\n"
              << "  Scenario:    an open ofstream is passed in.\n"
              << "  Pass criteria: exactly one line is written containing ITR, BREAK,\n"
              << "                 both molecule names and both interface indices.\n";

    const std::string fileName = "test_bii_assoc_dissoc.tmp";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(bii_make_template("ProtA", 0, /*canDestroy=*/false));
    molTemplateList.push_back(bii_make_template("Lipid", 1, /*canDestroy=*/false));

    Molecule mol1 = bii_make_molecule(/*index=*/8, /*molTypeIndex=*/0, 1, '\0');
    bii_bind_interface(mol1, 0, /*partnerIndex=*/9, 0, 3);
    Molecule mol2 = bii_make_molecule(/*index=*/9, /*molTypeIndex=*/1, 1, '\0');
    bii_bind_interface(mol2, 0, /*partnerIndex=*/8, 0, 3);

    std::vector<Molecule> moleculeList { mol1, mol2 };
    std::vector<Complex> complexList {};

    BackRxn currRxn;
    currRxn.isSymmetric = true;
    currRxn.productListNew.push_back(bii_make_rxn_iface("a", 0, 60, 0, '\0'));
    currRxn.productListNew.push_back(bii_make_rxn_iface("il", 1, 61, 0, '\0'));

    // Open the log file, run the break, then close so the data is flushed.
    std::ofstream logFile(fileName);
    ASSERT_TRUE(logFile.is_open()) << "Could not open temporary log file " << fileName;
    std::cerr << "  Calling break_interaction_implicitlipid with an open log stream (iter=123)...\n";
    break_interaction_implicitlipid(/*iter=*/123, /*relIface1=*/0, /*relIface2=*/0, moleculeList[0], moleculeList[1],
        currRxn, moleculeList, complexList, molTemplateList, logFile);
    logFile.close();

    // Read the log back in.
    std::ifstream in(fileName);
    EXPECT_TRUE(in.is_open()) << "Could not reopen the temporary log file for reading";
    std::string line;
    std::getline(in, line);
    std::string extra;
    const bool hasSecondLine = static_cast<bool>(std::getline(in, extra));
    in.close();
    std::remove(fileName.c_str()); // tidy up regardless of assertion outcome

    std::cerr << "  Logged line: \"" << line << "\"\n";

    // Content checks: the record format is
    //   ITR:<iter>,BREAK,<name1>,<idx1>,<relIface1>,<name2>,<idx2>,<relIface2>
    EXPECT_NE(line.find("ITR:123"), std::string::npos) << "Log line should record the iteration number";
    EXPECT_NE(line.find("BREAK"), std::string::npos) << "Log line should be tagged BREAK";
    EXPECT_NE(line.find("ProtA"), std::string::npos) << "Log line should contain reactMol1's template name";
    EXPECT_NE(line.find("Lipid"), std::string::npos) << "Log line should contain reactMol2's template name";
    EXPECT_NE(line.find(",8,"), std::string::npos) << "Log line should contain reactMol1's index (8)";
    EXPECT_NE(line.find(",9,"), std::string::npos) << "Log line should contain reactMol2's index (9)";
    EXPECT_FALSE(hasSecondLine) << "Exactly one record should be written per break";

    // The bookkeeping should still have happened alongside the logging.
    EXPECT_EQ(moleculeList[0].interfaceList[0].index, 60) << "Interface index should still be updated when logging";
    EXPECT_FALSE(moleculeList[0].interfaceList[0].isBound) << "Interface should still be unbound when logging";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario lives in its own TEST so that a failure
// in one does not prevent the others from running (all assertions above are
// non-fatal EXPECT_* apart from the unavoidable file-open ASSERT).
// -----------------------------------------------------------------------------
TEST(BreakInteractionImplicitLipid, SymmetricReactionUpdatesMol1) { test_bii_symmetric_reaction_updates_mol1(); }
TEST(BreakInteractionImplicitLipid, NonSymmetricMatchesFirstProduct) { test_bii_nonsymmetric_matches_first_product(); }
TEST(BreakInteractionImplicitLipid, NonSymmetricMatchesSecondProduct) { test_bii_nonsymmetric_matches_second_product(); }
TEST(BreakInteractionImplicitLipid, MonomerListUpdate) { test_bii_monomer_list_update(); }
TEST(BreakInteractionImplicitLipid, WritesAssocDissocRecord) { test_bii_writes_assoc_dissoc_record(); }