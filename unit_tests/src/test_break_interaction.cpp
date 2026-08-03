/*! \file test_break_interaction.cpp
 *
 * ### Unit test for src/reactions/break_interaction.cpp
 *
 * Function under test:
 *
 *     bool break_interaction(long long int iter, size_t relIface1, size_t relIface2,
 *                            Molecule& reactMol1, Molecule& reactMol2,
 *                            const BackRxn& currRxn,
 *                            std::vector<Molecule>& moleculeList,
 *                            std::vector<Complex>& complexList,
 *                            std::vector<MolTemplate>& molTemplateList,
 *                            int ILindexMol, const ForwardRxn& conjForwardRxn,
 *                            bool& breakLinkComplex, double timeStep,
 *                            std::ofstream& assocDissocFile)
 *
 * break_interaction() performs the book-keeping half of a dissociation event:
 *   1. the bond is removed from both partners' bndpartner lists,
 *   2. a slot for a possible new Complex is reserved (either an empty slot or a
 *      brand new entry at the end of complexList),
 *   3. determine_parent_complex_IL() decides whether the two molecules are still
 *      connected through some other path (a "loop"),
 *   4. if they are still connected, a rebinding (loop closure) correction factor
 *      is used to stochastically CANCEL the dissociation,
 *   5. otherwise the interfaces are freed, absolute interface indices are
 *      swapped to the (unbound) product indices, the complexes' properties are
 *      recomputed and (for a true split) a new Complex is created.
 *
 * The tests below build tiny, fully self-contained systems and check each of
 * these behaviours:
 *
 *   - A simple bound dimer  A(a!1).B(b!1)  ->  A(a) + B(b)            (true split)
 *   - The same dimer, with the roles of reactMol1/reactMol2 swapped
 *     (exercises the "which product matches which molecule" logic)
 *   - A closed triangle 0-1-2-0 where breaking 0-1 leaves them connected:
 *       * with a huge rebinding rate  -> dissociation is CANCELLED
 *       * with a tiny rebinding rate  -> dissociation proceeds but only the
 *         link is broken (breakLinkComplex == true, no new Complex)
 *
 * Verbose commentary is written to stderr so that the reader can follow which
 * source file / function is under test and what each assertion checks.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

// -----------------------------------------------------------------------------
// Local helpers.  Everything lives in an anonymous namespace and additionally
// carries a "bi_" (break_interaction) prefix so that nothing can collide with
// other translation units in the combined test binary.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Small container holding one complete miniature system. */
struct BiSystem {
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<MolTemplate> molTemplateList {};
    BackRxn backRxn {}; //!< the dissociation (back) reaction being performed
    ForwardRxn forwardRxn {}; //!< its conjugate forward (binding) reaction
};

/*! \brief Make sure the global GSL random number generator exists.
 *
 * break_interaction() calls rand_gsl() whenever the two molecules remain in the
 * same complex (loop breaking).  The suite-wide global `r` starts out as
 * nullptr, so initialise it here if nobody else did.
 */
void bi_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  (initialising GSL RNG for rand_gsl())\n";
        srand_gsl(12345);
    }
}

/*! \brief Reset all statics that break_interaction() reads or modifies.
 *
 * These are process-wide, so each test starts from a known state.
 */
void bi_reset_statics(size_t numMolTypes)
{
    Complex::emptyComList.clear();
    Complex::numberOfComplexes = 0;
    Complex::maxID = 0;
    Molecule::emptyMolList.clear();
    Molecule::numberOfMolecules = 0;
    Molecule::maxID = 0;
    MolTemplate::numMolTypes = static_cast<unsigned>(numMolTypes);
    MolTemplate::numEachMolType.assign(numMolTypes, 0);
    Parameters::dt = 1.0;
    Parameters::lastUpdateTransition.assign(numMolTypes, 0);
}

/*! \brief Build a MolTemplate with `numIfaces` interfaces.
 *
 * Non-zero mass/radius/diffusion constants are required because
 * Complex::update_properties() recomputes the complex diffusion constants from
 * these values.
 */
MolTemplate bi_make_template(const std::string& name, int typeIndex, int numIfaces, bool canDestroy)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord { 10.0, 10.0, 10.0 };
    temp.Dr = Coord { 0.1, 0.1, 0.1 };
    temp.canDestroy = canDestroy;
    temp.countTransition = false; // keep the transition-matrix code path out of the way
    temp.isImplicitLipid = false;
    temp.copies = 1;

    for (int i = 0; i < numIfaces; ++i) {
        Interface iface;
        iface.index = i;
        iface.name = "i" + std::to_string(i);
        iface.iCoord = Coord { (i == 0) ? 1.0 : -1.0, 0.0, 0.0 };
        iface.stateList.emplace_back(Interface::State('\0', i));
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Build a Molecule with `numIfaces` (initially free) interfaces. */
Molecule bi_make_molecule(int index, int molTypeIndex, int comIndex, const Coord& com, int numIfaces)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = comIndex;
    mol.complexId = comIndex;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isImplicitLipid = false;
    mol.isLipid = false;

    for (int i = 0; i < numIfaces; ++i) {
        Molecule::Iface iface;
        // interface coordinates are only used for the tiny "nudge apart" vector
        iface.coord = Coord { com.x + ((i == 0) ? 1.0 : -1.0), com.y, com.z };
        iface.index = 100 * molTypeIndex + i; // arbitrary "bound" absolute index
        iface.relIndex = i;
        iface.molTypeIndex = molTypeIndex;
        iface.stateIndex = 0;
        iface.stateIden = '\0';
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
        mol.freelist.push_back(i); // everything starts free
    }
    return mol;
}

/*! \brief Bind interface i1 of molecule m1 to interface i2 of molecule m2.
 *
 * Fills in exactly the fields that break_interaction() (and
 * determine_parent_complex_IL()) inspect: the Iface::interaction, the
 * bndlist/bndpartner/bndRxnList triplet and the freelist.
 */
void bi_bind(std::vector<Molecule>& molList, int m1, int i1, int m2, int i2, int backRxnIndex)
{
    Molecule& a = molList[m1];
    Molecule& b = molList[m2];

    a.interfaceList[i1].isBound = true;
    a.interfaceList[i1].interaction = Molecule::Interaction(m2, i2, backRxnIndex);
    a.bndlist.push_back(i1);
    a.bndpartner.push_back(m2);
    a.bndRxnList.push_back(backRxnIndex);
    a.freelist.erase(std::remove(a.freelist.begin(), a.freelist.end(), i1), a.freelist.end());

    b.interfaceList[i2].isBound = true;
    b.interfaceList[i2].interaction = Molecule::Interaction(m1, i1, backRxnIndex);
    b.bndlist.push_back(i2);
    b.bndpartner.push_back(m1);
    b.bndRxnList.push_back(backRxnIndex);
    b.freelist.erase(std::remove(b.freelist.begin(), b.freelist.end(), i2), b.freelist.end());
}

/*! \brief Build one Complex holding the given member molecules.
 *
 * numEachMol and lastNumberUpdateItrEachMol must already have one entry per
 * molecule type, because break_interaction() snapshots them before any call to
 * Complex::update_properties().
 */
Complex bi_make_complex(int index, const std::vector<int>& members,
    const std::vector<Molecule>& molList, size_t numMolTypes)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.isEmpty = false;
    com.memberList = members;
    com.numEachMol.assign(numMolTypes, 0);
    com.lastNumberUpdateItrEachMol.assign(numMolTypes, 0);

    Coord sum {};
    for (auto memMol : members) {
        sum.x += molList[memMol].comCoord.x;
        sum.y += molList[memMol].comCoord.y;
        sum.z += molList[memMol].comCoord.z;
        ++com.numEachMol[molList[memMol].molTypeIndex];
    }
    const double n = static_cast<double>(members.size());
    com.comCoord = Coord { sum.x / n, sum.y / n, sum.z / n };
    com.mass = n;
    com.radius = 2.0;
    com.D = Coord { 10.0, 10.0, 10.0 };
    com.Dr = Coord { 0.1, 0.1, 0.1 };
    return com;
}

/*! \brief Two molecules of two different types, bound through interface 0.
 *
 * Layout:  A(index 0, type 0, 2 ifaces) -- B(index 1, type 1, 2 ifaces)
 * The bond joins A.iface0 to B.iface0; iface1 of each is left free.
 *
 * Product (unbound) absolute interface indices are deliberately distinct (7 for
 * A, 9 for B) so we can verify which one is written back to which molecule.
 */
BiSystem bi_build_dimer_system(bool canDestroy)
{
    BiSystem sys;
    sys.molTemplateList.push_back(bi_make_template("A", 0, 2, canDestroy));
    sys.molTemplateList.push_back(bi_make_template("B", 1, 2, canDestroy));
    bi_reset_statics(sys.molTemplateList.size());

    sys.moleculeList.push_back(bi_make_molecule(0, 0, 0, Coord { -1.0, 0.0, 0.0 }, 2));
    sys.moleculeList.push_back(bi_make_molecule(1, 1, 0, Coord { 1.0, 0.0, 0.0 }, 2));
    bi_bind(sys.moleculeList, 0, 0, 1, 0, 0);

    sys.complexList.push_back(
        bi_make_complex(0, { 0, 1 }, sys.moleculeList, sys.molTemplateList.size()));
    Complex::numberOfComplexes = 1;
    Complex::maxID = 1;
    Molecule::numberOfMolecules = 2;

    // Back (dissociation) reaction:  A(a!1).B(b!1) -> A(a) + B(b)
    sys.backRxn.isSymmetric = false;
    sys.backRxn.rxnType = ReactionType::bimolecular;
    sys.backRxn.productListNew.push_back(RxnIface("a", /*molTypeIndex*/ 0, /*abs*/ 7, /*rel*/ 0, '\0', false));
    sys.backRxn.productListNew.push_back(RxnIface("b", /*molTypeIndex*/ 1, /*abs*/ 9, /*rel*/ 0, '\0', false));

    // Conjugate forward reaction: only consulted when a loop is broken.
    sys.forwardRxn.loopCoopFactor = 1.0;
    sys.forwardRxn.rateList.emplace_back(1.0, std::vector<std::vector<RxnIface>> {});
    return sys;
}

/*! \brief Three identical molecules bound in a closed triangle 0-1, 1-2, 2-0.
 *
 * Breaking the 0-1 bond leaves 0 and 1 connected through 2, i.e. this is the
 * "loop breaking" case which triggers the rebinding correction.
 *
 * \param[in] forwardRate rate of the conjugate binding reaction (drives the
 *                        loop-closure correction ratio).
 */
BiSystem bi_build_triangle_system(double forwardRate)
{
    BiSystem sys;
    sys.molTemplateList.push_back(bi_make_template("A", 0, 2, /*canDestroy*/ false));
    bi_reset_statics(sys.molTemplateList.size());

    sys.moleculeList.push_back(bi_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, 2));
    sys.moleculeList.push_back(bi_make_molecule(1, 0, 0, Coord { 2.0, 0.0, 0.0 }, 2));
    sys.moleculeList.push_back(bi_make_molecule(2, 0, 0, Coord { 1.0, 2.0, 0.0 }, 2));

    bi_bind(sys.moleculeList, 0, 0, 1, 0, 0); // the bond we will break
    bi_bind(sys.moleculeList, 1, 1, 2, 0, 0);
    bi_bind(sys.moleculeList, 2, 1, 0, 1, 0);

    sys.complexList.push_back(
        bi_make_complex(0, { 0, 1, 2 }, sys.moleculeList, sys.molTemplateList.size()));
    Complex::numberOfComplexes = 1;
    Complex::maxID = 1;
    Molecule::numberOfMolecules = 3;

    // Homodimer style reaction -> symmetric, so no product matching is needed.
    sys.backRxn.isSymmetric = true;
    sys.backRxn.rxnType = ReactionType::bimolecular;
    sys.backRxn.productListNew.push_back(RxnIface("a", 0, /*abs*/ 5, /*rel*/ 0, '\0', false));
    sys.backRxn.productListNew.push_back(RxnIface("a", 0, /*abs*/ 5, /*rel*/ 0, '\0', false));

    sys.forwardRxn.loopCoopFactor = 1.0;
    sys.forwardRxn.rateList.emplace_back(forwardRate, std::vector<std::vector<RxnIface>> {});
    return sys;
}

/*! \brief Convenience predicate: is `value` inside `vec`? */
bool bi_contains(const std::vector<int>& vec, int value)
{
    return std::find(vec.begin(), vec.end(), value) != vec.end();
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a plain dimer dissociates into two independent complexes.
// -----------------------------------------------------------------------------
void test_bi_dimer_full_dissociation()
{
    std::cerr << "\n[TEST] test_bi_dimer_full_dissociation\n"
              << "  Source file:   src/reactions/break_interaction.cpp\n"
              << "  Function:      break_interaction()\n"
              << "  Scenario:      bound dimer A(a!1).B(b!1) dissociates; the two\n"
              << "                 molecules are not connected by any other bond.\n"
              << "  Pass criteria: returns false (not cancelled), breakLinkComplex\n"
              << "                 false, both interfaces freed and given their\n"
              << "                 unbound absolute indices, molecules end up in\n"
              << "                 two different complexes, and the dissociation is\n"
              << "                 logged to the assoc/dissoc file.\n";

    bi_ensure_rng();
    BiSystem sys = bi_build_dimer_system(/*canDestroy*/ true);

    const int complexesBefore = Complex::numberOfComplexes;
    const size_t complexListSizeBefore = sys.complexList.size();

    // Open a temporary log file so the "assocDissocFile.is_open()" branch runs.
    const std::string logName = "test_break_interaction_log.tmp";
    std::ofstream logFile(logName);

    bool breakLinkComplex = true; // deliberately wrong; must be reset to false
    std::cerr << "  Calling break_interaction(iter=42, relIface1=0, relIface2=0)...\n";
    const bool cancelled = break_interaction(/*iter*/ 42, /*relIface1*/ 0, /*relIface2*/ 0,
        sys.moleculeList[0], sys.moleculeList[1], sys.backRxn, sys.moleculeList, sys.complexList,
        sys.molTemplateList, /*ILindexMol*/ -1, sys.forwardRxn, breakLinkComplex,
        /*timeStep*/ 1.0, logFile);
    logFile.close();

    // --- return values -------------------------------------------------------
    EXPECT_FALSE(cancelled) << "A non-loop dissociation must never be cancelled";
    EXPECT_FALSE(breakLinkComplex) << "A true split is not a 'link only' break";

    // --- interface state ----------------------------------------------------
    const Molecule& molA = sys.moleculeList[0];
    const Molecule& molB = sys.moleculeList[1];
    EXPECT_FALSE(molA.interfaceList[0].isBound) << "A.iface0 must be unbound after the break";
    EXPECT_FALSE(molB.interfaceList[0].isBound) << "B.iface0 must be unbound after the break";
    EXPECT_EQ(molA.interfaceList[0].interaction.partnerIndex, -1)
        << "A.iface0 interaction must be cleared";
    EXPECT_EQ(molB.interfaceList[0].interaction.partnerIndex, -1)
        << "B.iface0 interaction must be cleared";
    EXPECT_EQ(molA.interfaceList[0].index, 7)
        << "A.iface0 should get productListNew[0].absIfaceIndex (7)";
    EXPECT_EQ(molB.interfaceList[0].index, 9)
        << "B.iface0 should get productListNew[1].absIfaceIndex (9)";

    // --- bookkeeping lists --------------------------------------------------
    EXPECT_TRUE(molA.bndpartner.empty()) << "A should have no bound partners left";
    EXPECT_TRUE(molB.bndpartner.empty()) << "B should have no bound partners left";
    EXPECT_TRUE(molA.bndlist.empty()) << "A.bndlist should be empty";
    EXPECT_TRUE(molB.bndlist.empty()) << "B.bndlist should be empty";
    EXPECT_TRUE(bi_contains(molA.freelist, 0)) << "A.iface0 should be back on the freelist";
    EXPECT_TRUE(bi_contains(molB.freelist, 0)) << "B.iface0 should be back on the freelist";

    // --- complex bookkeeping ------------------------------------------------
    EXPECT_GE(sys.complexList.size(), complexListSizeBefore + 1)
        << "A brand new Complex slot should have been used for the split";
    EXPECT_NE(molA.myComIndex, molB.myComIndex)
        << "The two molecules must now live in different complexes";
    EXPECT_EQ(Complex::numberOfComplexes, complexesBefore + 1)
        << "The total number of complexes should grow by exactly one";
    std::cerr << "  A is now in complex " << molA.myComIndex << ", B in complex "
              << molB.myComIndex << " (complexList size " << sys.complexList.size() << ")\n";

    // --- monomerList update (canDestroy == true for both templates) ---------
    EXPECT_EQ(sys.molTemplateList[0].monomerList.size(), 1u)
        << "Template A should have recorded molecule 0 as a new monomer";
    EXPECT_EQ(sys.molTemplateList[1].monomerList.size(), 1u)
        << "Template B should have recorded molecule 1 as a new monomer";
    if (!sys.molTemplateList[0].monomerList.empty())
        EXPECT_EQ(sys.molTemplateList[0].monomerList[0], 0) << "monomerList should hold index 0";
    if (!sys.molTemplateList[1].monomerList.empty())
        EXPECT_EQ(sys.molTemplateList[1].monomerList[0], 1) << "monomerList should hold index 1";

    // --- assoc/dissoc logging ----------------------------------------------
    std::ifstream in(logName);
    std::string line;
    std::getline(in, line);
    in.close();
    std::remove(logName.c_str());
    std::cerr << "  assocDissocFile line: \"" << line << "\"\n";
    EXPECT_NE(line.find("BREAK"), std::string::npos) << "Log line should contain the BREAK tag";
    EXPECT_NE(line.find("ITR:42"), std::string::npos) << "Log line should contain the iteration";
    EXPECT_NE(line.find('A'), std::string::npos) << "Log line should name molecule type A";
    EXPECT_NE(line.find('B'), std::string::npos) << "Log line should name molecule type B";
}

// -----------------------------------------------------------------------------
// Test 2: the same dimer, but with reactMol1/reactMol2 swapped.  This exercises
//         the "which entry of productListNew matches reactMol1" branch.
// -----------------------------------------------------------------------------
void test_bi_dimer_reversed_reactant_order()
{
    std::cerr << "\n[TEST] test_bi_dimer_reversed_reactant_order\n"
              << "  Source file:   src/reactions/break_interaction.cpp\n"
              << "  Function:      break_interaction()\n"
              << "  Scenario:      identical dimer, but reactMol1 = B and reactMol2 = A,\n"
              << "                 so reactMol1 matches productListNew[1].\n"
              << "  Pass criteria: each molecule still receives ITS OWN product\n"
              << "                 absolute interface index (A->7, B->9).\n";

    bi_ensure_rng();
    BiSystem sys = bi_build_dimer_system(/*canDestroy*/ false);

    std::ofstream closedFile; // never opened -> the logging branch is skipped
    bool breakLinkComplex = false;

    std::cerr << "  Calling break_interaction with reactMol1 = molecule 1 (type B)...\n";
    const bool cancelled = break_interaction(/*iter*/ 7, /*relIface1*/ 0, /*relIface2*/ 0,
        sys.moleculeList[1], sys.moleculeList[0], sys.backRxn, sys.moleculeList, sys.complexList,
        sys.molTemplateList, -1, sys.forwardRxn, breakLinkComplex, 1.0, closedFile);

    EXPECT_FALSE(cancelled) << "A non-loop dissociation must never be cancelled";
    EXPECT_FALSE(breakLinkComplex) << "A true split is not a 'link only' break";

    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].index, 9)
        << "B (reactMol1) must get productListNew[1].absIfaceIndex (9)";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 7)
        << "A (reactMol2) must get productListNew[0].absIfaceIndex (7)";
    EXPECT_FALSE(sys.moleculeList[0].interfaceList[0].isBound) << "A.iface0 must be unbound";
    EXPECT_FALSE(sys.moleculeList[1].interfaceList[0].isBound) << "B.iface0 must be unbound";
    EXPECT_NE(sys.moleculeList[0].myComIndex, sys.moleculeList[1].myComIndex)
        << "The molecules must end up in different complexes";

    // canDestroy is false here, so nothing may be pushed onto monomerList.
    EXPECT_TRUE(sys.molTemplateList[0].monomerList.empty())
        << "monomerList must stay empty when canDestroy is false";
    EXPECT_TRUE(sys.molTemplateList[1].monomerList.empty())
        << "monomerList must stay empty when canDestroy is false";
}

// -----------------------------------------------------------------------------
// Test 3: loop breaking with an enormous rebinding rate -> dissociation is
//         cancelled and the system must be restored exactly.
// -----------------------------------------------------------------------------
void test_bi_loop_dissociation_cancelled()
{
    std::cerr << "\n[TEST] test_bi_loop_dissociation_cancelled\n"
              << "  Source file:   src/reactions/break_interaction.cpp\n"
              << "  Function:      break_interaction()\n"
              << "  Scenario:      closed triangle 0-1-2-0; the 0-1 bond is broken but\n"
              << "                 the molecules stay connected through molecule 2.\n"
              << "                 A huge conjugate binding rate makes the loop-closure\n"
              << "                 correction ratio ~1e-12, so the event is rejected.\n"
              << "  Pass criteria: returns true, bndpartner entries are restored, the\n"
              << "                 interfaces remain bound and no Complex is added.\n";

    bi_ensure_rng();
    // rateClose = rate * 0.602 * coop, poisson = timeStep * rateClose.
    // rate = 1e12 -> poisson ~ 6e11 -> correctionRatio ~ 1.7e-12 -> cancel.
    BiSystem sys = bi_build_triangle_system(/*forwardRate*/ 1.0e12);

    const size_t complexListSizeBefore = sys.complexList.size();
    const int complexesBefore = Complex::numberOfComplexes;

    std::ofstream closedFile;
    bool breakLinkComplex = true; // must be reset to false on entry

    std::cerr << "  Calling break_interaction (expecting cancellation)...\n";
    const bool cancelled = break_interaction(/*iter*/ 100, /*relIface1*/ 0, /*relIface2*/ 0,
        sys.moleculeList[0], sys.moleculeList[1], sys.backRxn, sys.moleculeList, sys.complexList,
        sys.molTemplateList, -1, sys.forwardRxn, breakLinkComplex, /*timeStep*/ 1.0, closedFile);

    EXPECT_TRUE(cancelled) << "With a huge rebinding rate the dissociation must be cancelled";
    EXPECT_FALSE(breakLinkComplex) << "breakLinkComplex must remain false when cancelled";

    // The only state break_interaction() touched before the decision was the
    // bndpartner list, which must have been restored.
    EXPECT_TRUE(bi_contains(sys.moleculeList[0].bndpartner, 1))
        << "molecule 0 must still list molecule 1 as a partner";
    EXPECT_TRUE(bi_contains(sys.moleculeList[1].bndpartner, 0))
        << "molecule 1 must still list molecule 0 as a partner";
    EXPECT_EQ(sys.moleculeList[0].bndpartner.size(), 2u)
        << "molecule 0 should still have two bound partners";
    EXPECT_EQ(sys.moleculeList[1].bndpartner.size(), 2u)
        << "molecule 1 should still have two bound partners";

    // Interfaces untouched.
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound)
        << "molecule 0 iface0 must still be bound";
    EXPECT_TRUE(sys.moleculeList[1].interfaceList[0].isBound)
        << "molecule 1 iface0 must still be bound";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, 1)
        << "molecule 0 iface0 interaction must be intact";
    EXPECT_TRUE(sys.moleculeList[0].bndlist.size() == 2u)
        << "molecule 0 bndlist must be untouched (2 bonds)";
    EXPECT_FALSE(bi_contains(sys.moleculeList[0].freelist, 0))
        << "molecule 0 iface0 must NOT be added to the freelist";

    // Complex bookkeeping restored.
    EXPECT_EQ(sys.complexList.size(), complexListSizeBefore)
        << "The speculatively added Complex slot must be removed again";
    EXPECT_EQ(Complex::numberOfComplexes, complexesBefore)
        << "The complex counter must not change on a cancelled dissociation";
    EXPECT_TRUE(Complex::emptyComList.empty())
        << "The speculative slot was at the end, so it is popped, not recycled";
    EXPECT_EQ(sys.moleculeList[0].myComIndex, sys.moleculeList[1].myComIndex)
        << "All molecules must remain in the original complex";
    std::cerr << "  complexList size = " << sys.complexList.size()
              << ", numberOfComplexes = " << Complex::numberOfComplexes << '\n';
}

// -----------------------------------------------------------------------------
// Test 4: loop breaking with a negligible rebinding rate -> the link is broken
//         but no new complex is created.
// -----------------------------------------------------------------------------
void test_bi_loop_link_break_proceeds()
{
    std::cerr << "\n[TEST] test_bi_loop_link_break_proceeds\n"
              << "  Source file:   src/reactions/break_interaction.cpp\n"
              << "  Function:      break_interaction()\n"
              << "  Scenario:      same closed triangle, but the conjugate binding rate\n"
              << "                 is ~0 so the correction ratio is ~1 and the event is\n"
              << "                 accepted.\n"
              << "  Pass criteria: returns false with breakLinkComplex == true, the 0-1\n"
              << "                 interfaces are freed, the remaining bonds survive and\n"
              << "                 NO new Complex is created.\n";

    bi_ensure_rng();
    // rate = 1e-12 -> poisson ~ 6e-13 -> correctionRatio ~ 1 -> accept.
    BiSystem sys = bi_build_triangle_system(/*forwardRate*/ 1.0e-12);

    const size_t complexListSizeBefore = sys.complexList.size();
    const int complexesBefore = Complex::numberOfComplexes;

    std::ofstream closedFile;
    bool breakLinkComplex = false;

    std::cerr << "  Calling break_interaction (expecting a link-only break)...\n";
    const bool cancelled = break_interaction(/*iter*/ 200, /*relIface1*/ 0, /*relIface2*/ 0,
        sys.moleculeList[0], sys.moleculeList[1], sys.backRxn, sys.moleculeList, sys.complexList,
        sys.molTemplateList, -1, sys.forwardRxn, breakLinkComplex, /*timeStep*/ 1.0, closedFile);

    EXPECT_FALSE(cancelled) << "With a negligible rebinding rate the break must proceed";
    EXPECT_TRUE(breakLinkComplex) << "Only the link is broken, so breakLinkComplex must be true";

    // The 0-1 bond is gone...
    EXPECT_FALSE(sys.moleculeList[0].interfaceList[0].isBound)
        << "molecule 0 iface0 must be unbound";
    EXPECT_FALSE(sys.moleculeList[1].interfaceList[0].isBound)
        << "molecule 1 iface0 must be unbound";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, -1)
        << "molecule 0 iface0 interaction must be cleared";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 5)
        << "molecule 0 iface0 should get the symmetric product abs index (5)";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].index, 5)
        << "molecule 1 iface0 should get the symmetric product abs index (5)";
    EXPECT_TRUE(bi_contains(sys.moleculeList[0].freelist, 0))
        << "molecule 0 iface0 must be on the freelist";
    EXPECT_TRUE(bi_contains(sys.moleculeList[1].freelist, 0))
        << "molecule 1 iface0 must be on the freelist";
    EXPECT_FALSE(bi_contains(sys.moleculeList[0].bndlist, 0))
        << "molecule 0 iface0 must be removed from bndlist";

    // ...but the other two bonds of the (former) ring survive.
    EXPECT_TRUE(bi_contains(sys.moleculeList[0].bndpartner, 2))
        << "molecule 0 must still be bound to molecule 2";
    EXPECT_TRUE(bi_contains(sys.moleculeList[1].bndpartner, 2))
        << "molecule 1 must still be bound to molecule 2";
    EXPECT_EQ(sys.moleculeList[2].bndpartner.size(), 2u)
        << "molecule 2 keeps both of its bonds";

    // No new complex: everything is still one connected object.
    EXPECT_EQ(sys.complexList.size(), complexListSizeBefore)
        << "The speculative Complex slot must be released again";
    EXPECT_EQ(Complex::numberOfComplexes, complexesBefore)
        << "No new complex is created when only a link is broken";
    EXPECT_EQ(sys.moleculeList[0].myComIndex, sys.moleculeList[1].myComIndex)
        << "Molecules 0 and 1 remain in the same complex";
    EXPECT_EQ(sys.moleculeList[0].myComIndex, sys.moleculeList[2].myComIndex)
        << "Molecule 2 remains in the same complex as well";

    // canDestroy is false and neither molecule became a monomer anyway.
    EXPECT_TRUE(sys.molTemplateList[0].monomerList.empty())
        << "No monomers are produced by a link-only break";
    std::cerr << "  complexList size = " << sys.complexList.size()
              << ", molecule complex indices = (" << sys.moleculeList[0].myComIndex << ", "
              << sys.moleculeList[1].myComIndex << ", " << sys.moleculeList[2].myComIndex << ")\n";
}

// -----------------------------------------------------------------------------
// Test 5: an empty Complex slot present in Complex::emptyComList must be reused
//         instead of growing complexList.
// -----------------------------------------------------------------------------
void test_bi_reuses_empty_complex_slot()
{
    std::cerr << "\n[TEST] test_bi_reuses_empty_complex_slot\n"
              << "  Source file:   src/reactions/break_interaction.cpp\n"
              << "  Function:      break_interaction()\n"
              << "  Scenario:      an empty (isEmpty == true) Complex already exists and\n"
              << "                 is registered in Complex::emptyComList.\n"
              << "  Pass criteria: the recycled slot is used (complexList does not grow,\n"
              << "                 emptyComList is emptied) and the split still happens.\n";

    bi_ensure_rng();
    BiSystem sys = bi_build_dimer_system(/*canDestroy*/ false);

    // Append an explicitly empty complex and advertise it as recyclable.
    Complex emptyCom;
    emptyCom.index = 1;
    emptyCom.isEmpty = true;
    emptyCom.numEachMol.assign(sys.molTemplateList.size(), 0);
    emptyCom.lastNumberUpdateItrEachMol.assign(sys.molTemplateList.size(), 0);
    sys.complexList.push_back(emptyCom);
    Complex::emptyComList.push_back(1);

    const size_t complexListSizeBefore = sys.complexList.size();

    std::ofstream closedFile;
    bool breakLinkComplex = false;

    std::cerr << "  Calling break_interaction with one recyclable complex slot...\n";
    const bool cancelled = break_interaction(/*iter*/ 5, 0, 0, sys.moleculeList[0],
        sys.moleculeList[1], sys.backRxn, sys.moleculeList, sys.complexList, sys.molTemplateList,
        -1, sys.forwardRxn, breakLinkComplex, 1.0, closedFile);

    EXPECT_FALSE(cancelled) << "A non-loop dissociation must never be cancelled";
    EXPECT_EQ(sys.complexList.size(), complexListSizeBefore)
        << "complexList must NOT grow when an empty slot is recycled";
    EXPECT_TRUE(Complex::emptyComList.empty())
        << "The recycled index must be removed from Complex::emptyComList";
    EXPECT_NE(sys.moleculeList[0].myComIndex, sys.moleculeList[1].myComIndex)
        << "The two molecules still have to end up in different complexes";
    std::cerr << "  complexList size = " << sys.complexList.size() << ", emptyComList size = "
              << Complex::emptyComList.size() << '\n';
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_bi_* helper runs inside its own TEST so
// that a failure in one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(BreakInteractionTest, DimerFullDissociation) { test_bi_dimer_full_dissociation(); }
TEST(BreakInteractionTest, DimerReversedReactantOrder) { test_bi_dimer_reversed_reactant_order(); }
TEST(BreakInteractionTest, LoopDissociationCancelled) { test_bi_loop_dissociation_cancelled(); }
TEST(BreakInteractionTest, LoopLinkBreakProceeds) { test_bi_loop_link_break_proceeds(); }
TEST(BreakInteractionTest, ReusesEmptyComplexSlot) { test_bi_reuses_empty_complex_slot(); }