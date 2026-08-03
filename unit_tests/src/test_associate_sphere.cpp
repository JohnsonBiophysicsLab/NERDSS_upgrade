/*! \file test_associate_sphere.cpp
 *
 * ### Unit test for src/reactions/associate_sphere.cpp
 *
 * Function under test:
 * \code
 * void associate_sphere(long long int iter, int ifaceIndex1, int ifaceIndex2,
 *                       Molecule& reactMol1, Molecule& reactMol2,
 *                       Complex& reactCom1, Complex& reactCom2,
 *                       const Parameters& params, ForwardRxn& currRxn,
 *                       std::vector<Molecule>& moleculeList,
 *                       std::vector<MolTemplate>& molTemplateList,
 *                       std::map<std::string,int>& observablesList,
 *                       copyCounters& counterArrays,
 *                       std::vector<Complex>& complexList,
 *                       Membrane& membraneObject,
 *                       const std::vector<ForwardRxn>& forwardRxns,
 *                       const std::vector<BackRxn>& backRxns,
 *                       std::ofstream& assocDissocFile);
 * \endcode
 *
 * associate_sphere() has four clearly distinguishable behaviours which are each
 * exercised below:
 *
 *   1. **Loop closure** - the two reacting molecules already live in the same
 *      Complex.  No geometry is changed; only the loop counter, the trajectory
 *      status and the bookkeeping tail (bond records, copy numbers, ...) run.
 *   2. **Successful association in solution** - two separate complexes are
 *      pushed to sigma (bindRadius), merged into one complex, and all
 *      bookkeeping is performed.
 *   3. **Cancelled association** - a steric-overlap check rejects the move, the
 *      temporary coordinates are discarded, the real coordinates are untouched
 *      and the function returns early (no bond is formed).
 *   4. **Misuse guard** - if either reactant is an implicit lipid the routine
 *      prints an error and terminates the process with exit code 1.
 *
 * All molecules are declared "point" molecules (MolTemplate::isPoint == true) so
 * that the theta/phi/omega rotations - which need meaningful association angles
 * - are skipped and the test focuses on the deterministic parts of the routine.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Small geometry helper: Euclidean distance between two coordinates.
// -----------------------------------------------------------------------------
double as_sphere_dist(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// -----------------------------------------------------------------------------
// A bundle holding every piece of state that associate_sphere() touches, so the
// individual tests stay readable.
// -----------------------------------------------------------------------------
struct AsSphereSystem {
    Parameters params {};
    Membrane membrane {};
    ForwardRxn rxn {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
};

/*! \brief Build a minimal single-interface molecule.
 *
 * The interface sits 1 nm away from the centre of mass so that the
 * interface-to-COM vectors used by the routine are well defined (non zero).
 */
Molecule as_sphere_make_molecule(int index, int molTypeIndex, int comIndex,
    const Coord& com, const Coord& ifaceCrd)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = comIndex;
    mol.complexId = comIndex;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isEmpty = false;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    // A single, free interface.
    Molecule::Iface oneIface {};
    oneIface.coord = ifaceCrd;
    oneIface.index = molTypeIndex; // absolute (state) index of the free iface
    oneIface.relIndex = 0;
    oneIface.stateIndex = 0;
    oneIface.molTypeIndex = molTypeIndex;
    oneIface.isBound = false;
    mol.interfaceList.push_back(oneIface);

    // The interface is currently free -> it must appear on the free list, since
    // associate_sphere() removes it from there when the bond forms.
    mol.freelist.push_back(0);

    return mol;
}

/*! \brief Build a Complex containing the supplied member molecule indices. */
Complex as_sphere_make_complex(int index, const std::vector<int>& members,
    const Coord& com, int numMolTypes)
{
    Complex targCom;
    targCom.index = index;
    targCom.id = index;
    targCom.comCoord = com;
    targCom.tmpComCoord = com;
    targCom.mass = static_cast<double>(members.size());
    targCom.radius = 1.0;
    targCom.D = Coord(10.0, 10.0, 10.0);
    targCom.Dr = Coord(0.01, 0.01, 0.01);
    targCom.memberList = members;
    targCom.numEachMol.assign(numMolTypes, 0);
    targCom.lastNumberUpdateItrEachMol.assign(numMolTypes, 0);
    targCom.OnSurface = false; // both complexes live in solution
    targCom.tmpOnSurface = false;
    targCom.isEmpty = false;
    targCom.ncross = 0;
    targCom.trajStatus = TrajStatus::none;
    return targCom;
}

/*! \brief Build a "point" MolTemplate so no orientational rotations happen. */
MolTemplate as_sphere_make_template(int molTypeIndex, const std::string& name)
{
    MolTemplate oneTemp;
    oneTemp.molName = name;
    oneTemp.molTypeIndex = molTypeIndex;
    oneTemp.isPoint = true; // -> associate_sphere skips theta/phi/omega
    oneTemp.isRod = false;
    oneTemp.isLipid = false;
    oneTemp.isImplicitLipid = false;
    oneTemp.checkOverlap = false; // no steric checks unless a test asks for them
    oneTemp.countTransition = false; // keeps the transition-matrix code inert
    oneTemp.canDestroy = false; // keeps the monomerList code inert
    oneTemp.mass = 1.0;
    oneTemp.radius = 1.0;
    oneTemp.D = Coord(10.0, 10.0, 10.0);
    oneTemp.Dr = Coord(0.01, 0.01, 0.01);

    // one interface with one state
    Interface oneIface("a", Coord(1.0, 0.0, 0.0));
    oneIface.index = 0;
    oneIface.stateList.push_back(Interface::State(molTypeIndex));
    oneTemp.interfaceList.push_back(oneIface);
    return oneTemp;
}

/*! \brief Reset the static members the routine (indirectly) relies on. */
void as_sphere_reset_statics()
{
    MolTemplate::numMolTypes = 2;
    MolTemplate::numEachMolType.assign(2, 1);
    MolTemplate::absToRelIface.assign(16, 0);

    Molecule::numberOfMolecules = 2;
    Molecule::emptyMolList.clear();

    Complex::numberOfComplexes = 2;
    Complex::currNumberComTypes = 2;
    Complex::currNumberMolTypes = 2;
    Complex::emptyComList.clear();

    Parameters::dt = 0.1;
    Parameters::lastUpdateTransition.assign(2, 0);
}

/*! \brief Populate a full, self-consistent two molecule system.
 *
 * \param[in]  sameComplex if true both molecules belong to one Complex (the
 *                         loop-closure branch), otherwise each molecule gets
 *                         its own Complex.
 */
void as_sphere_setup(AsSphereSystem& sys, bool sameComplex)
{
    as_sphere_reset_statics();

    // ---- simulation parameters -------------------------------------------
    sys.params.numMolTypes = 2;
    sys.params.numTotalSpecies = 3;
    sys.params.timeStep = 0.1;
    sys.params.overlapSepLimit = 0.1; // permissive by default
    sys.params.scaleMaxDisplace = 100.0; // permissive by default

    // ---- spherical boundary, large enough to never interfere -------------
    sys.membrane.isSphere = true;
    sys.membrane.isBox = false;
    sys.membrane.sphereR = 100.0;
    sys.membrane.implicitLipid = false;
    sys.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 400.0, 400.0, 400.0 });

    // ---- templates -------------------------------------------------------
    sys.molTemplateList.push_back(as_sphere_make_template(0, "A"));
    sys.molTemplateList.push_back(as_sphere_make_template(1, "B"));

    // ---- molecules -------------------------------------------------------
    // mol0: COM (0,0,0) iface (1,0,0)   mol1: COM (5,0,0) iface (4,0,0)
    // => interface separation is 3 nm, bindRadius will be 1 nm.
    const int com1Index = 0;
    const int com2Index = sameComplex ? 0 : 1;
    sys.moleculeList.push_back(
        as_sphere_make_molecule(0, 0, com1Index, Coord(0.0, 0.0, 0.0), Coord(1.0, 0.0, 0.0)));
    sys.moleculeList.push_back(
        as_sphere_make_molecule(1, 1, com2Index, Coord(5.0, 0.0, 0.0), Coord(4.0, 0.0, 0.0)));

    // ---- complexes -------------------------------------------------------
    if (sameComplex) {
        Complex oneCom = as_sphere_make_complex(0, { 0, 1 }, Coord(2.5, 0.0, 0.0), 2);
        oneCom.numEachMol[0] = 1;
        oneCom.numEachMol[1] = 1;
        sys.complexList.push_back(oneCom);
    } else {
        Complex com0 = as_sphere_make_complex(0, { 0 }, Coord(0.0, 0.0, 0.0), 2);
        com0.numEachMol[0] = 1;
        Complex com1 = as_sphere_make_complex(1, { 1 }, Coord(5.0, 0.0, 0.0), 2);
        com1.numEachMol[1] = 1;
        sys.complexList.push_back(com0);
        sys.complexList.push_back(com1);
    }

    // ---- the reaction ----------------------------------------------------
    sys.rxn.rxnType = ReactionType::bimolecular;
    sys.rxn.bindRadius = 1.0;
    sys.rxn.absRxnIndex = 0;
    sys.rxn.relRxnIndex = 0;
    sys.rxn.isReversible = true;
    sys.rxn.conjBackRxnIndex = 7; // arbitrary, we only check it is copied
    sys.rxn.isObserved = true;
    sys.rxn.observeLabel = "AB";
    // reactants: A(a) absIface 0, B(b) absIface 1 ; product absIface 2
    sys.rxn.reactantListNew.push_back(RxnIface("a", 0, 0, 0, '\0', false));
    sys.rxn.reactantListNew.push_back(RxnIface("b", 1, 1, 0, '\0', false));
    sys.rxn.productListNew.push_back(RxnIface("a-b", 0, 2, 0, '\0', false));
    sys.rxn.rateList.push_back(RxnBase::RateState(1.0, {}));
    // association angles stay NaN (default) - legal because both mols are points

    sys.observablesList["AB"] = 0;

    // ---- copy counters ---------------------------------------------------
    sys.counterArrays.copyNumSpecies.assign(16, 10);
    sys.counterArrays.nBoundPairs.assign(4, 0);
    sys.counterArrays.proPairlist.assign(4, 0);
    sys.counterArrays.singleDouble.assign(16, 0);
    sys.counterArrays.implicitDouble.assign(16, false);
    sys.counterArrays.canDissociate.assign(16, false);
    sys.counterArrays.bindPairList.assign(16, std::vector<int>());
    sys.counterArrays.bindPairListIL2D.assign(16, std::vector<int>());
    sys.counterArrays.bindPairListIL3D.assign(16, std::vector<int>());
    sys.counterArrays.eventArraySize = 20;
    sys.counterArrays.events3D.assign(400, 0);
    sys.counterArrays.events2D.assign(400, 0);
    sys.counterArrays.events3Dto2D.assign(400, 0);
}

} // namespace

// -----------------------------------------------------------------------------
// TEST 1: loop closure - both reactants already in the same Complex.
// -----------------------------------------------------------------------------
void test_as_sphere_loop_closure()
{
    std::cerr << "\n[TEST] test_as_sphere_loop_closure\n"
              << "  Source file : src/reactions/associate_sphere.cpp\n"
              << "  Function    : associate_sphere() (reactCom1.index == reactCom2.index)\n"
              << "  Scenario    : two molecules of ONE complex bind (ring closure).\n"
              << "  Criteria    : nLoops incremented, coordinates untouched, bond\n"
              << "                bookkeeping (partner indices, bndlist, freelist,\n"
              << "                copy numbers, observables, assoc/dissoc record) done.\n";

    AsSphereSystem sys;
    as_sphere_setup(sys, /*sameComplex=*/true);

    // The product species should be flagged dissociable so the bindPairList
    // branch is exercised as well.
    sys.counterArrays.canDissociate[2] = true;

    // Remember the coordinates: the loop-closure branch must not move anything.
    const Coord com0Before = sys.moleculeList[0].comCoord;
    const Coord com1Before = sys.moleculeList[1].comCoord;

    // Open a real assoc/dissoc log so the "BOND" record can be inspected.
    const std::string fileName = "test_associate_sphere_bond_record.txt";
    std::ofstream assocDissocFile(fileName);

    std::cerr << "  Calling associate_sphere (iter = 5, ifaces 0 and 0)...\n";
    associate_sphere(5, 0, 0, sys.moleculeList[0], sys.moleculeList[1],
        sys.complexList[0], sys.complexList[0], sys.params, sys.rxn,
        sys.moleculeList, sys.molTemplateList, sys.observablesList,
        sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, assocDissocFile);
    assocDissocFile.close();

    // --- loop counter ----------------------------------------------------
    EXPECT_EQ(sys.counterArrays.nLoops, 1)
        << "Binding inside one complex must increment the closed-loop counter";
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 0)
        << "Loop closure takes the short-circuit branch, so no 'assoc success' is counted";

    // --- coordinates must be untouched -----------------------------------
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].comCoord.x, com0Before.x)
        << "Loop closure performs no rotations/translations";
    EXPECT_DOUBLE_EQ(sys.moleculeList[1].comCoord.x, com1Before.x)
        << "Loop closure performs no rotations/translations";

    // --- trajectory status ------------------------------------------------
    EXPECT_TRUE(sys.moleculeList[0].trajStatus == TrajStatus::propagated)
        << "Members of the closing complex are marked as propagated";
    EXPECT_TRUE(sys.moleculeList[1].trajStatus == TrajStatus::propagated)
        << "Members of the closing complex are marked as propagated";

    // --- interaction bookkeeping -----------------------------------------
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound)
        << "Interface 0 of molecule 0 should now be bound";
    EXPECT_TRUE(sys.moleculeList[1].interfaceList[0].isBound)
        << "Interface 0 of molecule 1 should now be bound";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, 1)
        << "Molecule 0's partner index must be molecule 1";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].interaction.partnerIndex, 0)
        << "Molecule 1's partner index must be molecule 0";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIfaceIndex, 0)
        << "Partner interface index must be the reacting interface of the partner";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.conjBackRxn, 7)
        << "Reversible reaction must store its conjugate back-reaction index";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].interaction.conjBackRxn, 7)
        << "Reversible reaction must store its conjugate back-reaction index";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].index, 2)
        << "Bound interface must adopt the product's absolute interface index";

    // --- bound / free lists ----------------------------------------------
    EXPECT_EQ(sys.moleculeList[0].bndlist.size(), 1u) << "Interface added to bndlist";
    EXPECT_EQ(sys.moleculeList[1].bndlist.size(), 1u) << "Interface added to bndlist";
    EXPECT_EQ(sys.moleculeList[0].bndpartner.size(), 1u) << "Partner added to bndpartner";
    EXPECT_TRUE(sys.moleculeList[0].freelist.empty())
        << "The only free interface was consumed, freelist must be empty";
    EXPECT_TRUE(sys.moleculeList[1].freelist.empty())
        << "The only free interface was consumed, freelist must be empty";

    // --- copy numbers ----------------------------------------------------
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 9) << "Reactant 0 copy number decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], 9) << "Reactant 1 copy number decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[2], 11) << "Product copy number incremented";

    // --- dissociation pool ------------------------------------------------
    ASSERT_LT(2u, sys.counterArrays.bindPairList.size());
    EXPECT_EQ(sys.counterArrays.bindPairList[2].size(), 1u)
        << "A dissociable product must register one molecule in bindPairList";
    if (!sys.counterArrays.bindPairList[2].empty()) {
        EXPECT_EQ(sys.counterArrays.bindPairList[2][0], 0)
            << "The molecule holding reactant 0 (molecule index 0) is registered";
    }

    // --- observables -----------------------------------------------------
    EXPECT_EQ(sys.observablesList["AB"], 1)
        << "Observed reaction must increment its observable counter";

    // --- crossing lists ---------------------------------------------------
    EXPECT_EQ(sys.complexList[0].ncross, -1)
        << "The reacting complex must be flagged with ncross = -1";
    EXPECT_TRUE(sys.moleculeList[0].crossbase.empty()) << "crossbase must be cleared";
    EXPECT_TRUE(sys.moleculeList[1].crossbase.empty()) << "crossbase must be cleared";

    // --- assoc/dissoc log -------------------------------------------------
    std::ifstream in(fileName);
    std::stringstream buffer;
    buffer << in.rdbuf();
    in.close();
    const std::string logContents = buffer.str();
    std::cerr << "  assoc/dissoc log contents: " << logContents;
    EXPECT_NE(logContents.find("BOND"), std::string::npos)
        << "An open assoc/dissoc file must receive a BOND record";
    EXPECT_NE(logContents.find("ITR:5"), std::string::npos)
        << "The BOND record must contain the iteration number";
    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// TEST 2: successful association of two complexes in solution.
// -----------------------------------------------------------------------------
void test_as_sphere_solution_association_success()
{
    std::cerr << "\n[TEST] test_as_sphere_solution_association_success\n"
              << "  Source file : src/reactions/associate_sphere.cpp\n"
              << "  Function    : associate_sphere() (two distinct complexes)\n"
              << "  Scenario    : two point molecules 3 nm apart bind with a 1 nm\n"
              << "                binding radius inside a spherical boundary.\n"
              << "  Criteria    : interfaces end up exactly bindRadius apart, the\n"
              << "                second complex is destroyed and merged into the\n"
              << "                first, temporary coords are cleared and the\n"
              << "                success counter is incremented.\n";

    AsSphereSystem sys;
    as_sphere_setup(sys, /*sameComplex=*/false);
    sys.counterArrays.canDissociate[2] = true;

    std::ofstream noFile; // closed stream -> no assoc/dissoc record written

    std::cerr << "  Interface separation before: "
              << as_sphere_dist(sys.moleculeList[0].interfaceList[0].coord,
                     sys.moleculeList[1].interfaceList[0].coord)
              << " nm (bindRadius = " << sys.rxn.bindRadius << " nm)\n";

    associate_sphere(11, 0, 0, sys.moleculeList[0], sys.moleculeList[1],
        sys.complexList[0], sys.complexList[1], sys.params, sys.rxn,
        sys.moleculeList, sys.molTemplateList, sys.observablesList,
        sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, noFile);

    const double finalSep = as_sphere_dist(sys.moleculeList[0].interfaceList[0].coord,
        sys.moleculeList[1].interfaceList[0].coord);
    std::cerr << "  Interface separation after : " << finalSep << " nm\n";

    // --- the move was accepted -------------------------------------------
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 1)
        << "A geometrically legal association must be counted as a success";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapPartner, 0) << "No overlap cancellation expected";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapSystem, 0) << "No system overlap expected";
    EXPECT_EQ(sys.counterArrays.nCancelDisplace3D, 0) << "Displacement is small, no cancellation";
    EXPECT_EQ(sys.counterArrays.nLoops, 0) << "Distinct complexes must not count as a loop";

    // --- geometry: pushed to sigma ---------------------------------------
    EXPECT_NEAR(finalSep, sys.rxn.bindRadius, 1e-8)
        << "The two reacting interfaces must sit exactly at the binding radius";
    // Equal diffusion constants => each complex moves half of the excess (1 nm).
    EXPECT_NEAR(sys.moleculeList[0].comCoord.x, 1.0, 1e-8)
        << "Molecule 0 should have moved +1 nm along x (half of the 2 nm excess)";
    EXPECT_NEAR(sys.moleculeList[1].comCoord.x, 4.0, 1e-8)
        << "Molecule 1 should have moved -1 nm along x (half of the 2 nm excess)";

    // --- temporary association coordinates cleaned up --------------------
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty())
        << "Temporary association coordinates must be cleared after committing";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty())
        << "Temporary association coordinates must be cleared after committing";

    // --- complexes merged -------------------------------------------------
    EXPECT_EQ(sys.moleculeList[1].myComIndex, sys.complexList[0].index)
        << "Molecule 1 must be adopted by complex 0";
    EXPECT_EQ(sys.complexList[0].memberList.size(), 2u)
        << "Complex 0 must now contain both molecules";
    EXPECT_TRUE(sys.complexList[1].isEmpty)
        << "Complex 1 must be destroyed (marked empty) after the merge";
    EXPECT_TRUE(sys.complexList[1].memberList.empty())
        << "The destroyed complex must not keep any members";
    if (sys.complexList[0].numEachMol.size() >= 2) {
        EXPECT_EQ(sys.complexList[0].numEachMol[0], 1)
            << "update_properties() should count one molecule of type 0";
        EXPECT_EQ(sys.complexList[0].numEachMol[1], 1)
            << "update_properties() should count one molecule of type 1";
    }

    // --- trajectory status ------------------------------------------------
    EXPECT_TRUE(sys.moleculeList[0].trajStatus == TrajStatus::propagated)
        << "Associated molecules must be flagged as propagated";
    EXPECT_TRUE(sys.moleculeList[1].trajStatus == TrajStatus::propagated)
        << "Associated molecules must be flagged as propagated";

    // --- bond bookkeeping (same tail as the loop-closure branch) ----------
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound) << "Interface must be bound";
    EXPECT_TRUE(sys.moleculeList[1].interfaceList[0].isBound) << "Interface must be bound";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, 1)
        << "Partner index of molecule 0 must be 1";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 9) << "Reactant 0 decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[1], 9) << "Reactant 1 decremented";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[2], 11) << "Product incremented";
    EXPECT_EQ(sys.observablesList["AB"], 1) << "Observable must be incremented";
    EXPECT_EQ(sys.complexList[0].ncross, -1) << "Surviving complex flagged with ncross = -1";

    // --- bound pair counter (index layout is implementation defined, only the
    //     total number of registered bound pairs is checked) ---------------
    int totalBoundPairs = 0;
    for (auto n : sys.counterArrays.nBoundPairs)
        totalBoundPairs += n;
    EXPECT_EQ(totalBoundPairs, 1)
        << "Exactly one new bound pair (A-B) must be registered";
}

// -----------------------------------------------------------------------------
// TEST 3: association cancelled by a steric overlap check.
// -----------------------------------------------------------------------------
void test_as_sphere_cancel_on_overlap()
{
    std::cerr << "\n[TEST] test_as_sphere_cancel_on_overlap\n"
              << "  Source file : src/reactions/associate_sphere.cpp\n"
              << "  Function    : associate_sphere() (cancellation path)\n"
              << "  Scenario    : both molecule types request overlap checking and\n"
              << "                params.overlapSepLimit (5 nm) is larger than the\n"
              << "                post-association COM separation (3 nm).\n"
              << "  Criteria    : a cancel counter is incremented, no bond is made,\n"
              << "                real coordinates are unchanged and the temporary\n"
              << "                coordinates are discarded.\n";

    AsSphereSystem sys;
    as_sphere_setup(sys, /*sameComplex=*/false);

    // Force the structure-overlap check to reject the move.
    sys.molTemplateList[0].checkOverlap = true;
    sys.molTemplateList[1].checkOverlap = true;
    sys.params.overlapSepLimit = 5.0; // > final COM separation of 3 nm

    const Coord com0Before = sys.moleculeList[0].comCoord;
    const Coord com1Before = sys.moleculeList[1].comCoord;
    const Coord iface0Before = sys.moleculeList[0].interfaceList[0].coord;

    std::ofstream noFile;

    associate_sphere(3, 0, 0, sys.moleculeList[0], sys.moleculeList[1],
        sys.complexList[0], sys.complexList[1], sys.params, sys.rxn,
        sys.moleculeList, sys.molTemplateList, sys.observablesList,
        sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, noFile);

    const int totalCancels = sys.counterArrays.nCancelOverlapPartner
        + sys.counterArrays.nCancelOverlapSystem + sys.counterArrays.nCancelSpanBox
        + sys.counterArrays.nCancelDisplace2D + sys.counterArrays.nCancelDisplace3D
        + sys.counterArrays.nCancelDisplace3Dto2D;
    std::cerr << "  Cancel counters total = " << totalCancels
              << " (overlapPartner = " << sys.counterArrays.nCancelOverlapPartner
              << ", spanBox = " << sys.counterArrays.nCancelSpanBox << ")\n";

    EXPECT_GT(totalCancels, 0)
        << "An overlapping association must increment at least one cancel counter";
    EXPECT_EQ(sys.counterArrays.nCancelOverlapPartner, 1)
        << "Overlap between the two associating structures is counted as an "
           "'overlap partner' cancellation";
    EXPECT_EQ(sys.counterArrays.nAssocSuccess, 0)
        << "A cancelled association must not be counted as a success";

    // --- nothing may have changed ----------------------------------------
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].comCoord.x, com0Before.x)
        << "Real coordinates must be restored/untouched on cancellation";
    EXPECT_DOUBLE_EQ(sys.moleculeList[1].comCoord.x, com1Before.x)
        << "Real coordinates must be restored/untouched on cancellation";
    EXPECT_DOUBLE_EQ(sys.moleculeList[0].interfaceList[0].coord.x, iface0Before.x)
        << "Interface coordinates must be untouched on cancellation";
    EXPECT_TRUE(sys.moleculeList[0].tmpICoords.empty())
        << "Temporary association coordinates must be cleared on cancellation";
    EXPECT_TRUE(sys.moleculeList[1].tmpICoords.empty())
        << "Temporary association coordinates must be cleared on cancellation";

    // --- no bond was formed ----------------------------------------------
    EXPECT_FALSE(sys.moleculeList[0].interfaceList[0].isBound)
        << "No bond may be created when association is cancelled";
    EXPECT_FALSE(sys.moleculeList[1].interfaceList[0].isBound)
        << "No bond may be created when association is cancelled";
    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.partnerIndex, -1)
        << "Partner index must stay unset when association is cancelled";
    EXPECT_TRUE(sys.moleculeList[0].bndlist.empty()) << "bndlist must stay empty";
    EXPECT_EQ(sys.moleculeList[0].freelist.size(), 1u) << "freelist must be preserved";
    EXPECT_FALSE(sys.complexList[1].isEmpty)
        << "The second complex must survive a cancelled association";
    EXPECT_EQ(sys.moleculeList[1].myComIndex, 1)
        << "Molecule 1 must remain in its original complex";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[0], 10)
        << "Copy numbers must not change on cancellation";
    EXPECT_EQ(sys.counterArrays.copyNumSpecies[2], 10)
        << "Copy numbers must not change on cancellation";
    EXPECT_EQ(sys.observablesList["AB"], 0)
        << "Observables must not change on cancellation";
}

// -----------------------------------------------------------------------------
// TEST 4: irreversible + unobserved reaction (loop-closure branch reused).
// -----------------------------------------------------------------------------
void test_as_sphere_irreversible_unobserved()
{
    std::cerr << "\n[TEST] test_as_sphere_irreversible_unobserved\n"
              << "  Source file : src/reactions/associate_sphere.cpp\n"
              << "  Function    : associate_sphere() bookkeeping tail\n"
              << "  Scenario    : irreversible, unobserved, non dissociable product.\n"
              << "  Criteria    : conjBackRxn stays at its default (-1), the\n"
              << "                observable map is untouched and no entry is added\n"
              << "                to bindPairList.\n";

    AsSphereSystem sys;
    as_sphere_setup(sys, /*sameComplex=*/true);
    sys.rxn.isReversible = false; // -> conjBackRxn must not be written
    sys.rxn.isObserved = false; // -> observable must not be incremented
    sys.counterArrays.canDissociate[2] = false; // -> bindPairList untouched

    std::ofstream noFile;

    associate_sphere(1, 0, 0, sys.moleculeList[0], sys.moleculeList[1],
        sys.complexList[0], sys.complexList[0], sys.params, sys.rxn,
        sys.moleculeList, sys.molTemplateList, sys.observablesList,
        sys.counterArrays, sys.complexList, sys.membrane, sys.forwardRxns,
        sys.backRxns, noFile);

    EXPECT_EQ(sys.moleculeList[0].interfaceList[0].interaction.conjBackRxn, -1)
        << "Irreversible reactions must leave conjBackRxn at its default value";
    EXPECT_EQ(sys.moleculeList[1].interfaceList[0].interaction.conjBackRxn, -1)
        << "Irreversible reactions must leave conjBackRxn at its default value";
    EXPECT_EQ(sys.observablesList["AB"], 0)
        << "An unobserved reaction must not touch the observables map";
    ASSERT_LT(2u, sys.counterArrays.bindPairList.size());
    EXPECT_TRUE(sys.counterArrays.bindPairList[2].empty())
        << "A non dissociable product must not be added to the dissociation pool";

    // The bond itself must still be formed.
    EXPECT_TRUE(sys.moleculeList[0].interfaceList[0].isBound)
        << "The bond is still created for an irreversible reaction";
    EXPECT_EQ(sys.counterArrays.nLoops, 1) << "Loop counter still incremented";
}

// -----------------------------------------------------------------------------
// TEST 5: misuse guard - implicit lipids must not reach associate_sphere().
//
// This is a death test: associate_sphere() prints a diagnostic and calls
// exit(1) when either reactant is an implicit lipid.
// -----------------------------------------------------------------------------
#if GTEST_HAS_DEATH_TEST
void test_as_sphere_implicit_lipid_aborts()
{
    std::cerr << "\n[TEST] test_as_sphere_implicit_lipid_aborts\n"
              << "  Source file : src/reactions/associate_sphere.cpp\n"
              << "  Function    : associate_sphere() implicit-lipid guard\n"
              << "  Scenario    : reactant 2 is flagged as an implicit lipid while\n"
              << "                the two complexes differ.\n"
              << "  Criteria    : the process terminates with exit code 1.\n";

    AsSphereSystem sys;
    as_sphere_setup(sys, /*sameComplex=*/false);
    sys.moleculeList[1].isImplicitLipid = true; // wrong routine for this reactant

    std::ofstream noFile;

    // The empty regex matches any output; the exit status is what matters here
    // (the diagnostic itself is written to stdout, not stderr).
    EXPECT_EXIT(
        {
            associate_sphere(0, 0, 0, sys.moleculeList[0], sys.moleculeList[1],
                sys.complexList[0], sys.complexList[1], sys.params, sys.rxn,
                sys.moleculeList, sys.molTemplateList, sys.observablesList,
                sys.counterArrays, sys.complexList, sys.membrane,
                sys.forwardRxns, sys.backRxns, noFile);
        },
        ::testing::ExitedWithCode(1), "")
        << "associate_sphere must refuse implicit-lipid reactants with exit(1)";
}
#endif // GTEST_HAS_DEATH_TEST

// -----------------------------------------------------------------------------
// GoogleTest wrappers - one TEST per scenario so every scenario is reported
// individually and all of them run even if an earlier one fails.
// -----------------------------------------------------------------------------
TEST(AssociateSphere, LoopClosure) { test_as_sphere_loop_closure(); }
TEST(AssociateSphere, SolutionAssociationSuccess) { test_as_sphere_solution_association_success(); }
TEST(AssociateSphere, CancelOnOverlap) { test_as_sphere_cancel_on_overlap(); }
TEST(AssociateSphere, IrreversibleUnobserved) { test_as_sphere_irreversible_unobserved(); }
#if GTEST_HAS_DEATH_TEST
TEST(AssociateSphereDeathTest, ImplicitLipidAborts) { test_as_sphere_implicit_lipid_aborts(); }
#endif