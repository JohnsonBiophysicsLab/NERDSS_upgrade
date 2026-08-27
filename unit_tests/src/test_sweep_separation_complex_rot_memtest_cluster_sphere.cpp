/*! \file test_sweep_separation_complex_rot_memtest_cluster_sphere.cpp
 *
 * ### Unit tests for
 * ###   src/trajectory_functions/sweep_separation_complex_rot_memtest_cluster_sphere.cpp
 *
 * Function under test:
 *
 *   void sweep_separation_complex_rot_memtest_cluster_sphere(
 *            int simItr, int pro1Index, Parameters& params,
 *            std::vector<Molecule>& moleculeList,
 *            std::vector<Complex>& complexList,
 *            const std::vector<ForwardRxn>& forwardRxns,
 *            const std::vector<MolTemplate>& molTemplateList,
 *            const Membrane& membraneObject);
 *
 * The routine performs the "cluster sweep" overlap resolution used for a
 * spherical simulation boundary:
 *
 *   1. It builds the list of interacting molecule pairs of the cluster that
 *      contains `pro1Index` (via define_cluster_pairs()).
 *   2. If the cluster is large (> 30 pairs) it subdivides the time step into
 *      5 sub-steps, resampling every involved trajectory accordingly.
 *   3. For every pair it computes the *post-move* interface separation.  If
 *      the pair would end up closer than the reaction binding radius the
 *      trajectories are resampled (up to 100 attempts).
 *   4. All involved complexes are then physically propagated and every
 *      involved molecule is flagged TrajStatus::propagated.
 *   5. params.timeStep is always restored to its incoming value.
 *
 * The tests below build tiny, fully initialised systems (molecules, complexes,
 * templates, a spherical membrane and one bimolecular reaction) and check the
 * observable consequences of the above.
 *
 * NOTE: the "cluster sweep failed after 100 attempts" branch is deliberately
 *       *not* exercised - it delegates to sweep_separation_complex_rot_memtest_sphere()
 *       and can only be reached by constructing a system whose overlap can never
 *       be resolved (e.g. zero diffusion constants), which would make the test
 *       both slow and non-deterministic.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (file-local, prefixed sscrmcs_ so they cannot collide with any
// other translation unit in the combined gtest binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief (Re)initialise the shared GSL generator so the random trajectory
 *         resampling inside the function under test is reproducible.
 *
 *  `r` is the global generator defined in gtest_main.cpp; we simply point it at
 *  a freshly seeded generator.
 */
void sscrmcs_seed_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief Euclidean distance between two coordinates. */
double sscrmcs_dist(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Build a minimal, fully initialised MolTemplate with one interface.
 *
 *  Every field that Complex::update_properties() reads (radius, D, Dr) is set,
 *  because propagate() calls update_properties() on each moved complex.
 */
MolTemplate sscrmcs_make_template(int molTypeIndex, const std::string& name, double D, double Dr, bool isLipid)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.radius = 1.0;
    temp.mass = 1.0;
    temp.copies = 1;
    temp.isLipid = isLipid;
    temp.isImplicitLipid = isLipid;
    temp.isPromoter = false;
    temp.isRod = false;
    temp.isPoint = false;
    temp.checkOverlap = false;
    temp.D = Coord { D, D, D };
    temp.Dr = Coord { Dr, Dr, Dr };

    // one interface, sitting on top of the center of mass
    Interface iface { "a", Coord { 0.0, 0.0, 0.0 } };
    iface.index = 0;
    iface.stateList.emplace_back(std::string { "a" }, 0);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a molecule with a single interface offset from its COM. */
Molecule sscrmcs_make_molecule(
    int index, int comIndex, int molTypeIndex, const Coord& com, const Coord& ifaceOffset, bool isImplicitLipid)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = isImplicitLipid;
    mol.isImplicitLipid = isImplicitLipid;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.comCoord = com;
    mol.trajStatus = TrajStatus::none;

    Molecule::Iface iface;
    iface.coord = com + ifaceOffset;
    iface.relIndex = 0;
    iface.index = 0;
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a single-member complex around a molecule. */
Complex sscrmcs_make_complex(int index, const Molecule& mol, double D, double Dr, bool onSurface, size_t numMolTypes)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.comCoord = mol.comCoord;
    com.mass = mol.mass;
    com.radius = 1.0;
    com.memberList.push_back(mol.index);
    com.numEachMol = std::vector<int>(numMolTypes, 0);
    ++com.numEachMol[mol.molTypeIndex];
    com.lastNumberUpdateItrEachMol = std::vector<long long int>(numMolTypes, 0);
    com.D = Coord { D, D, D };
    com.Dr = Coord { Dr, Dr, Dr };
    com.isEmpty = false;
    com.OnSurface = onSurface;
    com.onFiber = false;
    com.tmpOnSurface = false;
    com.linksToSurface = 0;
    com.iLipidIndex = 0;
    com.ncross = 0;
    com.trajStatus = TrajStatus::none;
    com.trajTrans.x = 0.0;
    com.trajTrans.y = 0.0;
    com.trajTrans.z = 0.0;
    com.trajRot.x = 0.0;
    com.trajRot.y = 0.0;
    com.trajRot.z = 0.0;
    com.tmpComCoord = mol.comCoord;

    return com;
}

/*! \brief One symmetric bimolecular reaction between interface 0 of two molecules.
 *
 *  cluster_one_complex() reads only bindRadius and reactantListNew[*].relIfaceIndex.
 */
ForwardRxn sscrmcs_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.bindRadius = bindRadius;
    rxn.bindRadius2D = bindRadius;
    rxn.isReversible = false;
    rxn.conjBackRxnIndex = -1;
    rxn.reactantListNew.emplace_back(std::string { "a" }, 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back(std::string { "a" }, 0, 0, 0, '\0', false);
    rxn.productListNew = rxn.reactantListNew;
    rxn.rateList.emplace_back();
    return rxn;
}

/*! \brief A spherical membrane large enough that nothing in these tests reflects. */
Membrane sscrmcs_make_membrane(double sphereR, bool implicitLipid)
{
    Membrane mem;
    mem.isSphere = true;
    mem.isBox = false;
    mem.TwoD = false;
    mem.sphereR = sphereR;
    mem.sphereVol = (4.0 / 3.0) * M_PI * sphereR * sphereR * sphereR;
    mem.waterBox = Membrane::WaterBox(std::vector<double> { 4 * sphereR, 4 * sphereR, 4 * sphereR });
    mem.implicitLipid = implicitLipid;
    mem.hasCompartment = false;
    mem.compartmentR = 0.0;
    mem.nSites = 0;
    mem.nStates = 0;
    mem.No_free_lipids = 0;
    mem.No_protein = 0;
    mem.implicitlipidIndex = -1;
    mem.totalSA = 0.0;
    mem.Dx = 0.0;
    mem.Dy = 0.0;
    mem.Dz = 0.0;
    mem.Drx = 0.0;
    mem.Dry = 0.0;
    mem.Drz = 0.0;
    mem.offset = 0.0;
    mem.lipidLength = 0.0;
    mem.xBCtype = "reflect";
    mem.yBCtype = "reflect";
    mem.zBCtype = "reflect";
    // The RS3D look-up table is indexed at [300..399] and [400..499] whenever
    // implicitLipid == true, so it must be at least 500 elements long.
    mem.RS3Dvect.assign(500, 0.0);
    return mem;
}

/*! \brief Minimal Parameters object. */
Parameters sscrmcs_make_params(double timeStep)
{
    Parameters params;
    params.rank = 0;
    params.timeStep = timeStep;
    params.nItr = 10;
    params.itrRestartFrom = 0;
    params.numMolTypes = 1;
    params.numTotalComplex = 2;
    return params;
}

/*! \brief Register molecule `a` as a potential binding partner of molecule `b`. */
void sscrmcs_add_crossing(Molecule& mol, int partnerIndex, int myIfaceRelIndex, int rxnIndex)
{
    mol.crossbase.push_back(partnerIndex);
    mol.mycrossint.push_back(myIfaceRelIndex);
    mol.crossrxn.push_back(std::array<int, 3> { rxnIndex, 0, 0 });
    mol.probvec.push_back(0.0);
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: a molecule with no encountered partners produces an empty pair list,
//         so absolutely nothing is touched.
// -----------------------------------------------------------------------------
void test_sscrmcs_no_pairs_leaves_system_unchanged()
{
    std::cerr << "\n[TEST] test_sscrmcs_no_pairs_leaves_system_unchanged\n"
              << "  Source file: sweep_separation_complex_rot_memtest_cluster_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_cluster_sphere\n"
              << "  Scenario:    the target molecule has an empty crossbase list, so\n"
              << "               define_cluster_pairs() returns no pairs at all.\n"
              << "  Criteria:    coordinates, trajTrans and trajStatus are untouched and\n"
              << "               params.timeStep is left at its incoming value.\n";

    sscrmcs_seed_rng();

    std::vector<MolTemplate> molTemplateList { sscrmcs_make_template(0, "A", 1.0, 0.0, false) };
    std::vector<ForwardRxn> forwardRxns { sscrmcs_make_rxn(1.0) };
    Membrane membraneObject = sscrmcs_make_membrane(100.0, false);
    Parameters params = sscrmcs_make_params(1.0);

    // A single lonely molecule/complex - no crossings registered.
    std::vector<Molecule> moleculeList { sscrmcs_make_molecule(
        0, 0, 0, Coord { 3.0, -4.0, 5.0 }, Coord { 0.0, 0.0, 0.0 }, false) };
    std::vector<Complex> complexList { sscrmcs_make_complex(0, moleculeList[0], 1.0, 0.0, false, molTemplateList.size()) };

    // Give the complex a non-zero trajectory: it must NOT be applied because the
    // molecule is never entered into the pair list.
    complexList[0].trajTrans.x = 7.0;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_sphere(...)\n";
    sweep_separation_complex_rot_memtest_cluster_sphere(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 3.0) << "COM x must not move when no pairs exist";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, -4.0) << "COM y must not move when no pairs exist";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 5.0) << "COM z must not move when no pairs exist";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 7.0) << "Unused trajectory must not be consumed";
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must stay 'none' - the molecule was never swept";
    EXPECT_DOUBLE_EQ(params.timeStep, 1.0) << "params.timeStep must be restored/unchanged";

    std::cerr << "  Final COM = (" << moleculeList[0].comCoord.x << ", " << moleculeList[0].comCoord.y << ", "
              << moleculeList[0].comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: a well separated pair is propagated exactly by its own trajectory
//         (no resampling is triggered) and is flagged as propagated afterwards.
// -----------------------------------------------------------------------------
void test_sscrmcs_separated_pair_propagates_exactly()
{
    std::cerr << "\n[TEST] test_sscrmcs_separated_pair_propagates_exactly\n"
              << "  Source file: sweep_separation_complex_rot_memtest_cluster_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_cluster_sphere\n"
              << "  Scenario:    two single-molecule complexes 40 nm apart with a 1 nm\n"
              << "               binding radius; rotational diffusion is zero so\n"
              << "               Complex::propagate() takes its exact translation path.\n"
              << "  Criteria:    each COM/interface is displaced by exactly the preset\n"
              << "               trajTrans, trajTrans is zeroed, and both molecules end\n"
              << "               up flagged TrajStatus::propagated.\n";

    sscrmcs_seed_rng();

    std::vector<MolTemplate> molTemplateList { sscrmcs_make_template(0, "A", 1.0, 0.0, false) };
    std::vector<ForwardRxn> forwardRxns { sscrmcs_make_rxn(1.0) }; // bindRadius = 1 nm
    Membrane membraneObject = sscrmcs_make_membrane(100.0, false);
    Parameters params = sscrmcs_make_params(1.0);

    std::vector<Molecule> moleculeList {
        sscrmcs_make_molecule(0, 0, 0, Coord { -20.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }, false),
        sscrmcs_make_molecule(1, 1, 0, Coord { 20.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }, false)
    };
    // The two molecules "encountered" each other through interface 0 / reaction 0.
    sscrmcs_add_crossing(moleculeList[0], 1, 0, 0);
    sscrmcs_add_crossing(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList {
        sscrmcs_make_complex(0, moleculeList[0], 1.0, 0.0, false, molTemplateList.size()),
        sscrmcs_make_complex(1, moleculeList[1], 1.0, 0.0, false, molTemplateList.size())
    };

    // Deterministic, non-overlapping displacement: they approach by 2 nm only.
    complexList[0].trajTrans.x = 1.0;
    complexList[1].trajTrans.x = -1.0;

    std::cerr << "  Initial separation = " << sscrmcs_dist(moleculeList[0].comCoord, moleculeList[1].comCoord) << " nm\n";
    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_sphere(...)\n";
    sweep_separation_complex_rot_memtest_cluster_sphere(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Positions must be exactly the pre-set translation - no resampling happened.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, -19.0) << "molecule 0 should move exactly +1 nm in x";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 19.0) << "molecule 1 should move exactly -1 nm in x";
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, -19.0) << "interface 0 follows its molecule";
    EXPECT_DOUBLE_EQ(moleculeList[1].interfaceList[0].coord.x, 19.0) << "interface 1 follows its molecule";

    // update_properties() inside propagate() re-derives the complex COM.
    EXPECT_DOUBLE_EQ(complexList[0].comCoord.x, -19.0) << "complex 0 COM tracks its single member";
    EXPECT_DOUBLE_EQ(complexList[1].comCoord.x, 19.0) << "complex 1 COM tracks its single member";

    // propagate() zeroes the trajectory once it has been consumed.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.x, 0.0) << "trajTrans must be zeroed after propagation";

    // Everything that was swept is marked as propagated.
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "molecule 0 must be flagged propagated";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "molecule 1 must be flagged propagated";

    EXPECT_DOUBLE_EQ(params.timeStep, 1.0) << "params.timeStep must be restored";

    std::cerr << "  Final separation = " << sscrmcs_dist(moleculeList[0].comCoord, moleculeList[1].comCoord) << " nm\n";
}

// -----------------------------------------------------------------------------
// Test 3: an overlapping pair has its trajectories resampled until the
//         post-move interface separation exceeds the binding radius.
// -----------------------------------------------------------------------------
void test_sscrmcs_overlapping_pair_is_swept_apart()
{
    std::cerr << "\n[TEST] test_sscrmcs_overlapping_pair_is_swept_apart\n"
              << "  Source file: sweep_separation_complex_rot_memtest_cluster_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_cluster_sphere\n"
              << "  Scenario:    two complexes start 0.2 nm apart while the reaction\n"
              << "               binding radius is 2 nm, i.e. deeply overlapping.\n"
              << "               Diffusion is large enough that resampling resolves the\n"
              << "               overlap well within the 100 attempt limit.\n"
              << "  Criteria:    after the call the two interfaces are at least\n"
              << "               bindRadius apart and both molecules have moved.\n";

    sscrmcs_seed_rng();

    const double bindRadius = 2.0;

    std::vector<MolTemplate> molTemplateList { sscrmcs_make_template(0, "A", 10.0, 0.0, false) };
    std::vector<ForwardRxn> forwardRxns { sscrmcs_make_rxn(bindRadius) };
    Membrane membraneObject = sscrmcs_make_membrane(200.0, false);
    Parameters params = sscrmcs_make_params(1.0);

    const Coord start0 { -0.1, 0.0, 0.0 };
    const Coord start1 { 0.1, 0.0, 0.0 };

    std::vector<Molecule> moleculeList {
        sscrmcs_make_molecule(0, 0, 0, start0, Coord { 0.0, 0.0, 0.0 }, false),
        sscrmcs_make_molecule(1, 1, 0, start1, Coord { 0.0, 0.0, 0.0 }, false)
    };
    sscrmcs_add_crossing(moleculeList[0], 1, 0, 0);
    sscrmcs_add_crossing(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList {
        sscrmcs_make_complex(0, moleculeList[0], 10.0, 0.0, false, molTemplateList.size()),
        sscrmcs_make_complex(1, moleculeList[1], 10.0, 0.0, false, molTemplateList.size())
    };

    std::cerr << "  Initial interface separation = "
              << sscrmcs_dist(moleculeList[0].interfaceList[0].coord, moleculeList[1].interfaceList[0].coord)
              << " nm (bindRadius = " << bindRadius << " nm)\n";
    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_sphere(...)\n";
    sweep_separation_complex_rot_memtest_cluster_sphere(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    const double finalSep
        = sscrmcs_dist(moleculeList[0].interfaceList[0].coord, moleculeList[1].interfaceList[0].coord);
    std::cerr << "  Final interface separation   = " << finalSep << " nm\n";

    // The swept configuration is the one that was actually accepted, therefore
    // the final geometry must respect the binding radius.
    EXPECT_GE(finalSep, bindRadius - 1e-8) << "Overlap must be resolved: separation >= bindRadius";

    // Both molecules had to be displaced to escape the overlap.
    EXPECT_GT(sscrmcs_dist(moleculeList[0].comCoord, start0), 0.0) << "molecule 0 should have been displaced";
    EXPECT_GT(sscrmcs_dist(moleculeList[1].comCoord, start1), 0.0) << "molecule 1 should have been displaced";

    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "molecule 0 must be flagged propagated";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "molecule 1 must be flagged propagated";
    EXPECT_DOUBLE_EQ(params.timeStep, 1.0) << "params.timeStep must be restored";
}

// -----------------------------------------------------------------------------
// Test 4: partners that are implicit lipids are skipped by the overlap check,
//         and the implicit-lipid RS3D look-up branch is exercised.
// -----------------------------------------------------------------------------
void test_sscrmcs_implicit_lipid_partner_is_skipped()
{
    std::cerr << "\n[TEST] test_sscrmcs_implicit_lipid_partner_is_skipped\n"
              << "  Source file: sweep_separation_complex_rot_memtest_cluster_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_cluster_sphere\n"
              << "  Scenario:    membraneObject.implicitLipid == true (so the RS3D\n"
              << "               look-up table branch runs) and the only partner of the\n"
              << "               target molecule is an implicit lipid placed *inside*\n"
              << "               the binding radius.\n"
              << "  Criteria:    the pair is skipped (no resampling), so the target is\n"
              << "               translated by exactly its preset trajTrans; the implicit\n"
              << "               lipid (already propagated) never moves.\n";

    sscrmcs_seed_rng();

    std::vector<MolTemplate> molTemplateList {
        sscrmcs_make_template(0, "A", 1.0, 0.0, false), // ordinary protein
        sscrmcs_make_template(1, "IL", 0.0, 0.0, true) // implicit lipid
    };
    std::vector<ForwardRxn> forwardRxns { sscrmcs_make_rxn(5.0) }; // large bindRadius
    Membrane membraneObject = sscrmcs_make_membrane(100.0, true); // implicit lipid ON
    // molTypeIndex 0 maps to RS3Dvect[400]; the corresponding RS3D value lives at
    // RS3Dvect[300].  Both are zero here, so the reflection behaviour is unchanged.
    membraneObject.RS3Dvect[400] = 0.0;
    membraneObject.RS3Dvect[300] = 0.0;

    Parameters params = sscrmcs_make_params(1.0);

    const Coord ilStart { 0.5, 0.0, 0.0 };
    std::vector<Molecule> moleculeList {
        sscrmcs_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }, false),
        sscrmcs_make_molecule(1, 1, 1, ilStart, Coord { 0.0, 0.0, 0.0 }, true)
    };
    // Only the protein records the encounter; the implicit lipid keeps an empty
    // crossbase so no reverse pair is generated.
    sscrmcs_add_crossing(moleculeList[0], 1, 0, 0);
    // The implicit lipid is already done moving for this step.
    moleculeList[1].trajStatus = TrajStatus::propagated;

    std::vector<Complex> complexList {
        sscrmcs_make_complex(0, moleculeList[0], 1.0, 0.0, false, molTemplateList.size()),
        sscrmcs_make_complex(1, moleculeList[1], 0.0, 0.0, true, molTemplateList.size())
    };

    // A deliberate displacement that keeps the protein well inside the (huge)
    // binding radius of the implicit lipid.
    complexList[0].trajTrans.y = 0.25;

    std::cerr << "  Initial protein/lipid separation = "
              << sscrmcs_dist(moleculeList[0].comCoord, moleculeList[1].comCoord)
              << " nm (bindRadius = 5 nm -> would normally count as overlap)\n";
    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_sphere(...)\n";
    sweep_separation_complex_rot_memtest_cluster_sphere(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Because the partner is an implicit lipid the separation test is bypassed,
    // so no resampling occurred and the preset translation was applied verbatim.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 0.0) << "no x displacement was requested";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, 0.25) << "protein moves by exactly its preset trajTrans.y";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 0.0) << "no z displacement was requested";

    // The implicit lipid was already flagged propagated, so its complex is not moved.
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, ilStart.x) << "implicit lipid must not move (x)";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.y, ilStart.y) << "implicit lipid must not move (y)";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.z, ilStart.z) << "implicit lipid must not move (z)";

    // Both members of the pair are flagged propagated at the end of the routine.
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "protein must be flagged propagated";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "implicit lipid must be flagged propagated";

    EXPECT_DOUBLE_EQ(params.timeStep, 1.0) << "params.timeStep must be restored";

    std::cerr << "  Final protein COM = (" << moleculeList[0].comCoord.x << ", " << moleculeList[0].comCoord.y << ", "
              << moleculeList[0].comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 5: a large cluster (> 30 pairs) triggers the 5 sub-step path; the time
//         step is temporarily divided but must be restored on exit.
// -----------------------------------------------------------------------------
void test_sscrmcs_large_cluster_substeps_and_restores_timestep()
{
    std::cerr << "\n[TEST] test_sscrmcs_large_cluster_substeps_and_restores_timestep\n"
              << "  Source file: sweep_separation_complex_rot_memtest_cluster_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_cluster_sphere\n"
              << "  Scenario:    one hub molecule encounters 31 partners, so the pair\n"
              << "               list exceeds the 30 pair threshold and the routine\n"
              << "               splits the step into 5 smaller sub-steps.\n"
              << "  Criteria:    params.timeStep is back to its original value, every\n"
              << "               involved molecule is flagged propagated, all molecules\n"
              << "               moved by a small (sub-step sized) amount and no hub/\n"
              << "               partner pair ends up closer than the binding radius.\n";

    sscrmcs_seed_rng();

    const int numPartners = 31; // -> 31 pairs, one more than the maxPairs = 30 threshold
    const double bindRadius = 1.0;
    const double D = 0.001; // tiny, so the geometry stays well separated
    const double origTimeStep = 1.0;

    std::vector<MolTemplate> molTemplateList { sscrmcs_make_template(0, "A", D, 0.0, false) };
    std::vector<ForwardRxn> forwardRxns { sscrmcs_make_rxn(bindRadius) };
    Membrane membraneObject = sscrmcs_make_membrane(300.0, false);
    Parameters params = sscrmcs_make_params(origTimeStep);

    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;

    // Hub molecule (index 0), sitting far above the partner plane.
    moleculeList.push_back(
        sscrmcs_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 60.0 }, Coord { 0.0, 0.0, 0.0 }, false));

    // Partners laid out on a coarse grid in the z = 0 plane, 12 nm apart.
    for (int i = 1; i <= numPartners; ++i) {
        const double px = ((i - 1) % 8) * 12.0 - 42.0;
        const double py = ((i - 1) / 8) * 12.0 - 18.0;
        moleculeList.push_back(
            sscrmcs_make_molecule(i, i, 0, Coord { px, py, 0.0 }, Coord { 0.0, 0.0, 0.0 }, false));
        // Only the hub records the encounters, so each partner complex contributes
        // no additional pairs of its own.
        sscrmcs_add_crossing(moleculeList[0], i, 0, 0);
    }

    for (int i = 0; i <= numPartners; ++i)
        complexList.push_back(sscrmcs_make_complex(i, moleculeList[i], D, 0.0, false, molTemplateList.size()));

    // Remember where everybody started.
    std::vector<Coord> startCoords;
    for (const auto& mol : moleculeList)
        startCoords.push_back(mol.comCoord);

    std::cerr << "  Built " << moleculeList.size() << " molecules -> " << numPartners
              << " cluster pairs (threshold is 30)\n";
    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_sphere(...)\n";
    sweep_separation_complex_rot_memtest_cluster_sphere(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // The routine divides params.timeStep by nStep internally; it must put it back.
    EXPECT_DOUBLE_EQ(params.timeStep, origTimeStep)
        << "params.timeStep must be restored after the sub-stepped sweep";

    // Everything that participated must be flagged as propagated.
    for (int i = 0; i <= numPartners; ++i) {
        EXPECT_EQ(static_cast<int>(moleculeList[i].trajStatus), static_cast<int>(TrajStatus::propagated))
            << "molecule " << i << " must be flagged propagated";
    }

    // Every molecule should have been displaced, but only by a small, sub-step
    // sized amount (5 sub-steps of sqrt(2*dt/5*D) ~ 0.02 nm each).
    for (int i = 0; i <= numPartners; ++i) {
        const double disp = sscrmcs_dist(moleculeList[i].comCoord, startCoords[i]);
        EXPECT_GT(disp, 0.0) << "molecule " << i << " should have been propagated somewhere";
        EXPECT_LT(disp, 5.0) << "molecule " << i << " displacement should stay small for D = " << D;
    }

    // And no hub/partner pair may violate the binding radius.
    double minSep = 1e30;
    for (int i = 1; i <= numPartners; ++i)
        minSep = std::min(minSep, sscrmcs_dist(moleculeList[0].comCoord, moleculeList[i].comCoord));
    std::cerr << "  Minimum hub-partner separation after sweep = " << minSep << " nm\n";
    EXPECT_GE(minSep, bindRadius) << "no hub/partner pair may end up inside the binding radius";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - one per named test function so a failure in one does
// not stop the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(SweepSeparationClusterSphere, NoPairsLeavesSystemUnchanged)
{
    test_sscrmcs_no_pairs_leaves_system_unchanged();
}

TEST(SweepSeparationClusterSphere, SeparatedPairPropagatesExactly)
{
    test_sscrmcs_separated_pair_propagates_exactly();
}

TEST(SweepSeparationClusterSphere, OverlappingPairIsSweptApart)
{
    test_sscrmcs_overlapping_pair_is_swept_apart();
}

TEST(SweepSeparationClusterSphere, ImplicitLipidPartnerIsSkipped)
{
    test_sscrmcs_implicit_lipid_partner_is_skipped();
}

TEST(SweepSeparationClusterSphere, LargeClusterSubstepsAndRestoresTimestep)
{
    test_sscrmcs_large_cluster_substeps_and_restores_timestep();
}