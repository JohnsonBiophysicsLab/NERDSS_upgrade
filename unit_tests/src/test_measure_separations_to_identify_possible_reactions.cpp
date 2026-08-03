/*! \file test_measure_separations_to_identify_possible_reactions.cpp
 *
 * ### Unit test for src/reactions/measure_separations_to_identify_possible_reactions.cpp
 *
 * Function under test:
 *
 *     void measure_separations_to_identify_possible_reactions(
 *         unsigned simItr, Parameters&, std::vector<Molecule>&,
 *         std::vector<Complex>&, SimulVolume&, std::vector<ForwardRxn>&,
 *         std::vector<BackRxn>&, std::vector<CreateDestructRxn>&,
 *         std::vector<MolTemplate>&, std::map<std::string,int>&,
 *         copyCounters&, Membrane&, std::vector<double>&, std::vector<double>&,
 *         std::vector<double>&, std::vector<gsl_matrix*>&,
 *         std::vector<gsl_matrix*>&, std::vector<gsl_matrix*>&, int, double*,
 *         unsigned&)
 *
 * The routine walks every SimulVolume sub-cell, and for every molecule in that
 * sub-cell it
 *
 *   1. skips implicit-lipid molecules (`continue`),
 *   2. skips ghosted (MPI) molecules (`continue`),
 *   3. only continues if the molecule has free interfaces OR its MolTemplate
 *      requests bound-interface volume exclusion,
 *   4. optionally evaluates implicit-lipid binding (only when the MolTemplate
 *      has `bindToSurface == true`),
 *   5. evaluates every *later* molecule in the same sub-cell, and
 *   6. evaluates every molecule in each neighbouring sub-cell.
 *
 * The observable side effect of a "possible reaction" being identified is that
 * `check_bimolecular_reactions()` appends to the participating molecules'
 * `probvec` / `crossbase` / `mycrossint` / `crossrxn` lists and increments the
 * parent complexes' `ncross` counters.  The tests below therefore build tiny,
 * fully self-consistent systems and inspect those bookkeeping vectors.
 *
 * All console output goes to stderr so a reader can follow exactly which branch
 * is being probed and what the pass criteria are.  Only EXPECT_* macros are
 * used so that every test runs to completion even when one fails.
 */

#include <gsl/gsl_matrix.h>
#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"

// -----------------------------------------------------------------------------
// Helpers, kept in an anonymous namespace so nothing collides with the rest of
// the generated unit-test suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Make sure the global GSL RNG exists.
 *
 * `gtest_main.cpp` defines `gsl_rng* r = nullptr;`.  Some of the reaction
 * probability helpers may draw random numbers, so initialise the generator once
 * if nobody else has done it yet.
 */
void msipr_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  (global GSL RNG was null -> initialising with seed 1)\n";
        srand_gsl(1);
    }
}

/*! \brief Build a MolTemplate carrying exactly one interface with one state.
 *
 * \param[in] molTypeIndex  Index of this template in molTemplateList.
 * \param[in] name          Human readable molecule name.
 * \param[in] absIfaceIndex Absolute (state) index of the single interface.
 * \param[in] partnerAbsIfaceIndex Absolute index of the interface this state
 *                                 can react with (-1 => no partners at all).
 * \param[in] hasForwardRxn Whether the state points at forward reaction 0.
 */
MolTemplate msipr_make_template(int molTypeIndex, const std::string& name, int absIfaceIndex,
    int partnerAbsIfaceIndex, bool hasForwardRxn)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;

    // Diffusion constants in nm^2/us.  All three components are non-zero so
    // that the reaction is treated as a 3D reaction (no 2D lookup tables are
    // needed, which keeps the test fast and free of gsl matrix allocation).
    temp.D = Coord { 10.0, 10.0, 10.0 };
    temp.Dr = Coord { 0.1, 0.1, 0.1 };

    // A single interface sitting on the center of mass.
    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord { 0.0, 0.0, 0.0 };

    // One (default, unnamed) state for that interface.
    Interface::State state('\0', absIfaceIndex);
    if (partnerAbsIfaceIndex >= 0)
        state.rxnPartners.push_back(static_cast<unsigned>(partnerAbsIfaceIndex));
    if (hasForwardRxn)
        state.myForwardRxns.push_back(0u);
    iface.stateList.push_back(state);

    temp.interfaceList.push_back(iface);

    // Nothing exotic: no implicit lipid, no surface adsorption, no exclusion.
    temp.isImplicitLipid = false;
    temp.isLipid = false;
    temp.bindToSurface = false;
    temp.excludeVolumeBound = false;
    temp.checkOverlap = false;

    return temp;
}

/*! \brief Build a single-interface Molecule at a given coordinate.
 *
 * \param[in] index         Index of the molecule in moleculeList.
 * \param[in] comIndex      Index of its parent complex in complexList.
 * \param[in] molTypeIndex  Index of its MolTemplate.
 * \param[in] absIfaceIndex Absolute (state) index of its single interface.
 * \param[in] com           Center-of-mass coordinate (interface is coincident).
 * \param[in] isFree        If true the interface is put on the free list.
 */
Molecule msipr_make_molecule(int index, int comIndex, int molTypeIndex, int absIfaceIndex,
    const Coord& com, bool isFree)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isImplicitLipid = false;
    mol.isGhosted = false;
    mol.trajStatus = TrajStatus::none;

    // One interface, unbound, coincident with the COM.
    Molecule::Iface iface;
    iface.coord = com;
    iface.index = absIfaceIndex;
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.stateIden = '\0';
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    iface.excludeVolume = false;
    mol.interfaceList.push_back(iface);

    // Only free interfaces make the molecule a candidate for binding.
    if (isFree)
        mol.freelist.push_back(0);

    return mol;
}

/*! \brief Build a Complex owning exactly one member molecule. */
Complex msipr_make_complex(int index, int memberMolIndex, const Coord& com)
{
    Complex com1;
    com1.index = index;
    com1.comCoord = com;
    com1.radius = 1.0;
    com1.mass = 1.0;
    com1.memberList.push_back(memberMolIndex);
    com1.numEachMol.assign(2, 0);
    com1.lastNumberUpdateItrEachMol.assign(2, 0);
    com1.D = Coord { 10.0, 10.0, 10.0 };
    com1.Dr = Coord { 0.1, 0.1, 0.1 };
    com1.isEmpty = false;
    com1.OnSurface = false;
    com1.ncross = 0;
    com1.trajStatus = TrajStatus::none;
    return com1;
}

/*! \brief Build one bimolecular ForwardRxn between two single-interface types. */
ForwardRxn msipr_make_forward_rxn(int molType1, int absIface1, int molType2, int absIface2,
    double bindRadius, double rate)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = false;
    rxn.conjBackRxnIndex = -1;
    rxn.isOnMem = false;
    rxn.isSymmetric = false;
    rxn.hasStateChange = false;
    rxn.bindRadius = bindRadius;
    rxn.bindRadius2D = bindRadius;
    rxn.length3Dto2D = 2.0 * bindRadius;
    rxn.rxnLabel = "testBimolecular";

    // Reactants: (molType1, absIface1) + (molType2, absIface2), both unbound.
    rxn.reactantListNew.emplace_back("a", molType1, absIface1, 0, '\0', false);
    rxn.reactantListNew.emplace_back("a", molType2, absIface2, 0, '\0', false);

    // Products: the same interfaces, now bound.
    rxn.productListNew.emplace_back("a", molType1, absIface1 + 2, 0, '\0', true);
    rxn.productListNew.emplace_back("a", molType2, absIface2 + 2, 0, '\0', true);

    rxn.intReactantList = { absIface1, absIface2 };
    rxn.intProductList = { absIface1 + 2, absIface2 + 2 };

    // One rate, and one (empty) ancillary-interface list per reactant so that
    // hasIntangibles() can safely index both entries.
    rxn.rateList.emplace_back(rate, std::vector<std::vector<RxnIface>>(2));

    return rxn;
}

/*!
 * \brief Container bundling everything the function under test needs.
 *
 * Default construction produces a 100 nm reflecting cubic box, a 0.1 us
 * timestep, two molecule templates (A and B) able to react through forward
 * reaction 0, an empty SimulVolume, and empty molecule/complex lists.  Each
 * test then adds exactly the molecules and sub-cells it needs.
 */
struct MsiprSystem {
    Parameters params {};
    Membrane membrane {};
    SimulVolume simulVolume {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
    std::vector<double> IL2DbindingVec {};
    std::vector<double> IL2DUnbindingVec {};
    std::vector<double> ILTableIDs {};
    std::vector<gsl_matrix*> normMatrices {};
    std::vector<gsl_matrix*> survMatrices {};
    std::vector<gsl_matrix*> pirMatrices {};
    double tableIDs[100] {};
    unsigned DDTableIndex { 0 };
    int implicitlipidIndex { -1 };

    MsiprSystem()
    {
        // ---- Simulation parameters -------------------------------------
        params.numMolTypes = 2;
        params.numTotalSpecies = 4;
        params.timeStep = 0.1; // microseconds
        params.nItr = 10;
        params.numTotalComplex = 2;
        params.debugParams.verbosity = 0;

        // ---- Boundary: reflecting cubic box, no implicit lipid ----------
        membrane.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
        membrane.isBox = true;
        membrane.isSphere = false;
        membrane.implicitLipid = false;
        membrane.TwoD = false;
        membrane.xBCtype = "reflect";
        membrane.yBCtype = "reflect";
        membrane.zBCtype = "reflect";
        membrane.No_free_lipids = 0;
        membrane.nStates = 0;
        membrane.RS3Dvect.assign(10, 0.0);

        // ---- Molecule templates: A(a) + B(a) <-> A(a!1).B(a!1) ---------
        // A has absolute interface index 0, B has absolute index 1.
        molTemplateList.push_back(msipr_make_template(0, "A", 0, /*partner*/ 1, true));
        molTemplateList.push_back(msipr_make_template(1, "B", 1, /*partner*/ 0, true));

        // ---- The bimolecular reaction itself ---------------------------
        forwardRxns.push_back(msipr_make_forward_rxn(0, 0, 1, 1, /*bindRadius*/ 1.0, /*rate*/ 10.0));

        // ---- Species counters (sized generously; unused in 3D path) ----
        counterArrays.copyNumSpecies.assign(8, 0);
        counterArrays.nBoundPairs.assign(4, 0);
        counterArrays.proPairlist.assign(4, 0);
        counterArrays.singleDouble.assign(8, 0);
        counterArrays.implicitDouble.assign(8, false);
        counterArrays.canDissociate.assign(8, false);
        counterArrays.bindPairList.assign(8, std::vector<int> {});
        counterArrays.bindPairListIL2D.assign(8, std::vector<int> {});
        counterArrays.bindPairListIL3D.assign(8, std::vector<int> {});
        counterArrays.events3D.assign(counterArrays.eventArraySize, 0);
        counterArrays.events2D.assign(counterArrays.eventArraySize, 0);
        counterArrays.events3Dto2D.assign(counterArrays.eventArraySize, 0);

        // ---- SimulVolume skeleton (cells added per test) ----------------
        simulVolume.subCellSize = Coord { 50.0, 50.0, 50.0 };
        simulVolume.numSubCells.x = 0;
        simulVolume.numSubCells.y = 0;
        simulVolume.numSubCells.z = 0;
        simulVolume.numSubCells.tot = 0;
    }

    /*! \brief Append a sub-cell with the given members and neighbour cells. */
    void add_cell(int absIndex, const std::vector<int>& members, const std::vector<int>& neighbors)
    {
        SimulVolume::SubVolume cell;
        cell.absIndex = absIndex;
        cell.xIndex = absIndex;
        cell.yIndex = 0;
        cell.zIndex = 0;
        cell.memberMolList = members;
        cell.neighborList = neighbors;
        simulVolume.subCellList.push_back(cell);
        simulVolume.numSubCells.x = static_cast<int>(simulVolume.subCellList.size());
        simulVolume.numSubCells.y = 1;
        simulVolume.numSubCells.z = 1;
        simulVolume.numSubCells.tot = static_cast<int>(simulVolume.subCellList.size());

        // Keep each molecule's cached sub-volume index consistent.
        for (int molIndex : members) {
            if (molIndex >= 0 && molIndex < static_cast<int>(moleculeList.size()))
                moleculeList[molIndex].mySubVolIndex = absIndex;
        }
    }

    /*! \brief Add an A/B pair of molecules (and their complexes) at two points. */
    void add_AB_pair(const Coord& posA, const Coord& posB)
    {
        moleculeList.push_back(msipr_make_molecule(0, 0, /*molType*/ 0, /*absIface*/ 0, posA, true));
        moleculeList.push_back(msipr_make_molecule(1, 1, /*molType*/ 1, /*absIface*/ 1, posB, true));
        complexList.push_back(msipr_make_complex(0, 0, posA));
        complexList.push_back(msipr_make_complex(1, 1, posB));
    }

    /*! \brief Invoke the function under test. */
    void run(unsigned simItr)
    {
        measure_separations_to_identify_possible_reactions(simItr, params, moleculeList, complexList,
            simulVolume, forwardRxns, backRxns, createDestructRxns, molTemplateList, observablesList,
            counterArrays, membrane, IL2DbindingVec, IL2DUnbindingVec, ILTableIDs, normMatrices,
            survMatrices, pirMatrices, implicitlipidIndex, tableIDs, DDTableIndex);
    }

    /*! \brief Total number of candidate reactions recorded on all molecules. */
    size_t total_recorded_candidates() const
    {
        size_t total { 0 };
        for (const auto& mol : moleculeList)
            total += mol.probvec.size();
        return total;
    }

    /*! \brief Sum of ncross over all complexes. */
    int total_ncross() const
    {
        int total { 0 };
        for (const auto& com : complexList)
            total += com.ncross;
        return total;
    }
};

/*! \brief Dump the per-molecule bookkeeping vectors for a system to stderr. */
void msipr_report(const MsiprSystem& sys)
{
    for (size_t i { 0 }; i < sys.moleculeList.size(); ++i) {
        const auto& mol = sys.moleculeList[i];
        std::cerr << "    mol[" << i << "]: probvec=" << mol.probvec.size()
                  << " crossbase=" << mol.crossbase.size()
                  << " mycrossint=" << mol.mycrossint.size()
                  << " crossrxn=" << mol.crossrxn.size();
        if (!mol.probvec.empty())
            std::cerr << " prob[0]=" << mol.probvec[0];
        std::cerr << '\n';
    }
    for (size_t i { 0 }; i < sys.complexList.size(); ++i)
        std::cerr << "    complex[" << i << "]: ncross=" << sys.complexList[i].ncross << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: an empty SimulVolume must be a complete no-op.
// -----------------------------------------------------------------------------
void msipr_test_empty_volume_is_noop()
{
    std::cerr << "\n[TEST] msipr_test_empty_volume_is_noop\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      outer loop over simulVolume.subCellList (zero cells)\n"
              << "  Scenario:    two reactive molecules exist 2 nm apart but no sub-cell\n"
              << "               references them, so the loops never execute.\n"
              << "  Criteria:    no probvec entries, ncross == 0, DDTableIndex untouched.\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    sys.add_AB_pair(Coord { 0.0, 0.0, 0.0 }, Coord { 2.0, 0.0, 0.0 });
    // Deliberately add NO sub-cells.

    std::cerr << "  Calling measure_separations_to_identify_possible_reactions...\n";
    sys.run(0);
    msipr_report(sys);

    EXPECT_EQ(sys.total_recorded_candidates(), 0u)
        << "No molecule may be examined when subCellList is empty";
    EXPECT_EQ(sys.total_ncross(), 0)
        << "No complex may register a crossing when subCellList is empty";
    EXPECT_EQ(sys.DDTableIndex, 0u)
        << "DDTableIndex must be untouched when nothing is evaluated";
}

// -----------------------------------------------------------------------------
// Test 2: implicit-lipid and ghosted molecules are skipped by the two
//         `continue` statements at the top of the member loop.
// -----------------------------------------------------------------------------
void msipr_test_skips_implicit_lipid_and_ghosted()
{
    std::cerr << "\n[TEST] msipr_test_skips_implicit_lipid_and_ghosted\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      'if (isImplicitLipid) continue;' and\n"
              << "               'if (isGhosted) continue;'\n"
              << "  Scenario:    the only sub-cell holds one implicit lipid and one\n"
              << "               ghosted molecule, both with free interfaces and only\n"
              << "               2 nm apart.\n"
              << "  Criteria:    both are skipped, so nothing is recorded at all.\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    sys.add_AB_pair(Coord { 0.0, 0.0, 0.0 }, Coord { 2.0, 0.0, 0.0 });

    // Molecule 0 pretends to be the implicit lipid, molecule 1 lives on a
    // neighbouring MPI rank (ghosted).
    sys.moleculeList[0].isImplicitLipid = true;
    sys.moleculeList[1].isGhosted = true;

    sys.add_cell(0, { 0, 1 }, {});

    std::cerr << "  Calling measure_separations_to_identify_possible_reactions...\n";
    sys.run(0);
    msipr_report(sys);

    EXPECT_EQ(sys.total_recorded_candidates(), 0u)
        << "Implicit-lipid and ghosted molecules must be skipped by the outer loop";
    EXPECT_EQ(sys.total_ncross(), 0)
        << "No crossings may be registered for skipped molecules";
}

// -----------------------------------------------------------------------------
// Test 3: molecules with no free interfaces and no volume-exclusion request
//         never reach check_bimolecular_reactions().
// -----------------------------------------------------------------------------
void msipr_test_no_free_interfaces_is_noop()
{
    std::cerr << "\n[TEST] msipr_test_no_free_interfaces_is_noop\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      'if (freelist.size() > 0 || excludeVolumeBound)'\n"
              << "  Scenario:    two molecules 2 nm apart, both with an EMPTY freelist and\n"
              << "               templates that do not request bound volume exclusion.\n"
              << "  Criteria:    the guard is false, so nothing is recorded.\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    // Build the pair, then strip the free lists.
    sys.add_AB_pair(Coord { 0.0, 0.0, 0.0 }, Coord { 2.0, 0.0, 0.0 });
    sys.moleculeList[0].freelist.clear();
    sys.moleculeList[1].freelist.clear();
    sys.molTemplateList[0].excludeVolumeBound = false;
    sys.molTemplateList[1].excludeVolumeBound = false;

    sys.add_cell(0, { 0, 1 }, {});

    std::cerr << "  Calling measure_separations_to_identify_possible_reactions...\n";
    sys.run(0);
    msipr_report(sys);

    EXPECT_EQ(sys.total_recorded_candidates(), 0u)
        << "Molecules without free interfaces must not be evaluated for binding";
    EXPECT_EQ(sys.total_ncross(), 0)
        << "No crossings may be registered for fully bound molecules";
}

// -----------------------------------------------------------------------------
// Test 4: a close pair inside the SAME sub-cell is identified as a possible
//         reaction (the `memItr2 = memItr + 1` inner loop).
// -----------------------------------------------------------------------------
void msipr_test_same_cell_close_pair_is_found()
{
    std::cerr << "\n[TEST] msipr_test_same_cell_close_pair_is_found\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      inner loop over later members of the SAME sub-cell\n"
              << "  Scenario:    A(a) and B(a) 3 nm apart (bindRadius 1 nm, Dtot 20 nm^2/us,\n"
              << "               dt 0.1 us => Rmax ~ 11.4 nm) in one sub-cell, with a\n"
              << "               matching bimolecular forward reaction.\n"
              << "  Criteria:    both molecules record a candidate reaction, they point at\n"
              << "               each other via crossbase, ncross is incremented and the\n"
              << "               probability lies in [0,1].\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    sys.add_AB_pair(Coord { 0.0, 0.0, 0.0 }, Coord { 3.0, 0.0, 0.0 });
    sys.add_cell(0, { 0, 1 }, {});

    std::cerr << "  Calling measure_separations_to_identify_possible_reactions...\n";
    sys.run(0);
    msipr_report(sys);

    // Both partners should have registered exactly one candidate reaction.
    EXPECT_FALSE(sys.moleculeList[0].probvec.empty())
        << "Molecule 0 should record a candidate bimolecular reaction";
    EXPECT_FALSE(sys.moleculeList[1].probvec.empty())
        << "Molecule 1 should record a candidate bimolecular reaction";

    // crossbase must reference the partner molecule index.
    if (!sys.moleculeList[0].crossbase.empty()) {
        EXPECT_EQ(sys.moleculeList[0].crossbase[0], 1)
            << "Molecule 0's crossbase entry should point at molecule 1";
    }
    if (!sys.moleculeList[1].crossbase.empty()) {
        EXPECT_EQ(sys.moleculeList[1].crossbase[0], 0)
            << "Molecule 1's crossbase entry should point at molecule 0";
    }

    // The interface recorded must be the (only) relative interface 0.
    if (!sys.moleculeList[0].mycrossint.empty()) {
        EXPECT_EQ(sys.moleculeList[0].mycrossint[0], 0)
            << "Only relative interface 0 exists on molecule 0";
    }

    // Both parent complexes must know they have a pending crossing.
    EXPECT_GT(sys.total_ncross(), 0)
        << "At least one complex must register a crossing for a close reactive pair";

    // The stored association probability must be a valid probability.
    if (!sys.moleculeList[0].probvec.empty()) {
        const double prob = sys.moleculeList[0].probvec[0];
        EXPECT_GE(prob, 0.0) << "Association probability must be >= 0";
        EXPECT_LE(prob, 1.0) << "Association probability must be <= 1";
    }
}

// -----------------------------------------------------------------------------
// Test 5: a pair that is far beyond Rmax is examined but not recorded.
// -----------------------------------------------------------------------------
void msipr_test_far_pair_is_not_recorded()
{
    std::cerr << "\n[TEST] msipr_test_far_pair_is_not_recorded\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      inner same-cell loop, separation > Rmax\n"
              << "  Scenario:    identical setup to the previous test but the molecules are\n"
              << "               45 nm apart, far outside Rmax (~11.4 nm).\n"
              << "  Criteria:    the pair is inspected but no candidate reaction is stored.\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    sys.add_AB_pair(Coord { -22.5, 0.0, 0.0 }, Coord { 22.5, 0.0, 0.0 });
    sys.add_cell(0, { 0, 1 }, {});

    std::cerr << "  Calling measure_separations_to_identify_possible_reactions...\n";
    sys.run(0);
    msipr_report(sys);

    EXPECT_EQ(sys.total_recorded_candidates(), 0u)
        << "A pair separated by 45 nm is well outside Rmax and must not be recorded";
    EXPECT_EQ(sys.total_ncross(), 0)
        << "No crossing may be registered for a pair outside Rmax";
}

// -----------------------------------------------------------------------------
// Test 6: a close pair split across two NEIGHBOURING sub-cells is found via the
//         neighbour-cell loop.
// -----------------------------------------------------------------------------
void msipr_test_neighbor_cell_pair_is_found()
{
    std::cerr << "\n[TEST] msipr_test_neighbor_cell_pair_is_found\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      loop over simulVolume.subCellList[cellItr].neighborList\n"
              << "  Scenario:    molecule 0 lives in cell 0, molecule 1 in cell 1; cell 0\n"
              << "               lists cell 1 as a neighbour and the two molecules are\n"
              << "               only 4 nm apart.\n"
              << "  Criteria:    exactly one candidate reaction is recorded per molecule\n"
              << "               (the pair is not double counted because cell 1 has no\n"
              << "               neighbours of its own).\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    sys.add_AB_pair(Coord { -2.0, 0.0, 0.0 }, Coord { 2.0, 0.0, 0.0 });
    sys.add_cell(0, { 0 }, { 1 }); // cell 0 sees cell 1
    sys.add_cell(1, { 1 }, {}); // cell 1 has no forward neighbours

    std::cerr << "  Calling measure_separations_to_identify_possible_reactions...\n";
    sys.run(0);
    msipr_report(sys);

    EXPECT_EQ(sys.moleculeList[0].probvec.size(), 1u)
        << "Molecule 0 should record exactly one candidate via the neighbour cell";
    EXPECT_EQ(sys.moleculeList[1].probvec.size(), 1u)
        << "Molecule 1 should record exactly one candidate via the neighbour cell";
    EXPECT_GT(sys.total_ncross(), 0)
        << "Complexes must register the crossing found across neighbouring cells";
}

// -----------------------------------------------------------------------------
// Test 7: with no reaction partners declared on the interface state, close
//         molecules are inspected but no reaction is possible.
// -----------------------------------------------------------------------------
void msipr_test_no_reaction_partners_records_nothing()
{
    std::cerr << "\n[TEST] msipr_test_no_reaction_partners_records_nothing\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      same-cell inner loop, but the interface states declare no\n"
              << "               reaction partners\n"
              << "  Scenario:    two free molecules 3 nm apart whose templates were built\n"
              << "               without rxnPartners / myForwardRxns.\n"
              << "  Criteria:    the molecules are evaluated but nothing is recorded.\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    // Replace the reactive templates with inert ones (no partners, no rxns).
    sys.molTemplateList.clear();
    sys.molTemplateList.push_back(msipr_make_template(0, "A", 0, /*partner*/ -1, false));
    sys.molTemplateList.push_back(msipr_make_template(1, "B", 1, /*partner*/ -1, false));

    sys.add_AB_pair(Coord { 0.0, 0.0, 0.0 }, Coord { 3.0, 0.0, 0.0 });
    sys.add_cell(0, { 0, 1 }, {});

    std::cerr << "  Calling measure_separations_to_identify_possible_reactions...\n";
    sys.run(0);
    msipr_report(sys);

    EXPECT_EQ(sys.total_recorded_candidates(), 0u)
        << "Without declared reaction partners no candidate reaction may be stored";
    EXPECT_EQ(sys.total_ncross(), 0)
        << "Without declared reaction partners no crossing may be stored";
}

// -----------------------------------------------------------------------------
// Test 8: the excludeVolumeBound half of the guard is honoured even when the
//         molecule has no free interfaces.
// -----------------------------------------------------------------------------
void msipr_test_exclude_volume_bound_guard()
{
    std::cerr << "\n[TEST] msipr_test_exclude_volume_bound_guard\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      'if (freelist.size() > 0 || excludeVolumeBound == true)'\n"
              << "  Scenario:    both molecules have EMPTY free lists but their templates\n"
              << "               set excludeVolumeBound = true, so the body is entered.\n"
              << "               Because no interface is actually bound, there is nothing\n"
              << "               to exclude.\n"
              << "  Criteria:    the call completes safely and records no reactions.\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    sys.add_AB_pair(Coord { 0.0, 0.0, 0.0 }, Coord { 3.0, 0.0, 0.0 });
    sys.moleculeList[0].freelist.clear();
    sys.moleculeList[1].freelist.clear();
    sys.molTemplateList[0].excludeVolumeBound = true;
    sys.molTemplateList[1].excludeVolumeBound = true;
    sys.add_cell(0, { 0, 1 }, {});

    std::cerr << "  Calling measure_separations_to_identify_possible_reactions...\n";
    sys.run(0);
    msipr_report(sys);

    // The guard lets execution in, but with no bound interfaces there is no
    // exclusion pair and no binding candidate.
    EXPECT_EQ(sys.total_recorded_candidates(), 0u)
        << "No bound interfaces exist, so no exclusion/binding candidate may be stored";
    EXPECT_EQ(sys.total_ncross(), 0)
        << "No crossings expected when nothing is bound and nothing is free";
}

// -----------------------------------------------------------------------------
// Test 9: repeated calls accumulate (the routine never clears the lists itself,
//         that is the caller's job at the top of each timestep).
// -----------------------------------------------------------------------------
void msipr_test_repeated_calls_accumulate()
{
    std::cerr << "\n[TEST] msipr_test_repeated_calls_accumulate\n"
              << "  Source file: measure_separations_to_identify_possible_reactions.cpp\n"
              << "  Branch:      whole routine, invoked twice in a row\n"
              << "  Scenario:    a close reactive pair is evaluated for simItr 0 and again\n"
              << "               for simItr 1 without the caller clearing probvec.\n"
              << "  Criteria:    the routine does not clear state itself, so the number of\n"
              << "               recorded candidates after the 2nd call is >= after the 1st.\n";

    msipr_ensure_rng();

    MsiprSystem sys;
    sys.add_AB_pair(Coord { 0.0, 0.0, 0.0 }, Coord { 3.0, 0.0, 0.0 });
    sys.add_cell(0, { 0, 1 }, {});

    std::cerr << "  First call (simItr = 0)...\n";
    sys.run(0);
    const size_t afterFirst = sys.total_recorded_candidates();
    const int ncrossFirst = sys.total_ncross();
    std::cerr << "    candidates after first call = " << afterFirst
              << ", total ncross = " << ncrossFirst << '\n';

    std::cerr << "  Second call (simItr = 1), lists intentionally NOT cleared...\n";
    sys.run(1);
    const size_t afterSecond = sys.total_recorded_candidates();
    const int ncrossSecond = sys.total_ncross();
    std::cerr << "    candidates after second call = " << afterSecond
              << ", total ncross = " << ncrossSecond << '\n';
    msipr_report(sys);

    EXPECT_GE(afterSecond, afterFirst)
        << "The routine must not remove previously recorded candidates";
    EXPECT_GE(ncrossSecond, ncrossFirst)
        << "The routine must not decrement previously recorded crossings";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so that a
// failure in one scenario still lets the remaining scenarios execute.
// -----------------------------------------------------------------------------
TEST(MeasureSeparationsToIdentifyPossibleReactions, EmptyVolumeIsNoop)
{
    msipr_test_empty_volume_is_noop();
}
TEST(MeasureSeparationsToIdentifyPossibleReactions, SkipsImplicitLipidAndGhosted)
{
    msipr_test_skips_implicit_lipid_and_ghosted();
}
TEST(MeasureSeparationsToIdentifyPossibleReactions, NoFreeInterfacesIsNoop)
{
    msipr_test_no_free_interfaces_is_noop();
}
TEST(MeasureSeparationsToIdentifyPossibleReactions, SameCellClosePairIsFound)
{
    msipr_test_same_cell_close_pair_is_found();
}
TEST(MeasureSeparationsToIdentifyPossibleReactions, FarPairIsNotRecorded)
{
    msipr_test_far_pair_is_not_recorded();
}
TEST(MeasureSeparationsToIdentifyPossibleReactions, NeighborCellPairIsFound)
{
    msipr_test_neighbor_cell_pair_is_found();
}
TEST(MeasureSeparationsToIdentifyPossibleReactions, NoReactionPartnersRecordsNothing)
{
    msipr_test_no_reaction_partners_records_nothing();
}
TEST(MeasureSeparationsToIdentifyPossibleReactions, ExcludeVolumeBoundGuard)
{
    msipr_test_exclude_volume_bound_guard();
}
TEST(MeasureSeparationsToIdentifyPossibleReactions, RepeatedCallsAccumulate)
{
    msipr_test_repeated_calls_accumulate();
}