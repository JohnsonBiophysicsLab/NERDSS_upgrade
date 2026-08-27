/*! \file test_sweep_separation_complex_rot_memtest_cluster_box.cpp
 *
 * ### Unit test for
 * ###   src/trajectory_functions/sweep_separation_complex_rot_memtest_cluster_box.cpp
 *
 * The function under test is:
 *
 *     void sweep_separation_complex_rot_memtest_cluster_box(
 *              int simItr, int pro1Index, Parameters& params,
 *              std::vector<Molecule>& moleculeList,
 *              std::vector<Complex>& complexList,
 *              const std::vector<ForwardRxn>& forwardRxns,
 *              const std::vector<MolTemplate>& molTemplateList,
 *              const Membrane& membraneObject)
 *
 * What the routine does (and therefore what we test):
 *
 *   1. It builds a *cluster* of interacting pairs starting from the complex of
 *      molecule `pro1Index` (via define_cluster_pairs()).
 *   2. For every pair it predicts the interface-interface separation that the
 *      already-sampled trajectories (trajTrans / trajRot) would produce.  If
 *      the predicted separation is smaller than the reaction binding radius the
 *      trajectories are *resampled* (resample_traj()) until the overlap is
 *      gone, or until 100 attempts have been made.
 *   3. Pairs whose partner is an implicit lipid are skipped entirely by the
 *      overlap test.
 *   4. All involved complexes are then physically propagated, their trajectory
 *      vectors are zeroed, and their molecules are flagged
 *      TrajStatus::propagated.
 *   5. If the cluster contains more than 30 pairs the time step is temporarily
 *      divided into 10 sub-steps; the original params.timeStep must always be
 *      restored on exit.
 *
 * Everything is built from minimal, fully-initialised objects so that no
 * uninitialised member is ever dereferenced by the production code.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

// The RNG pointer is defined in unit_tests/src/gtest_main.cpp; we only use it.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Small helpers used to build a completely initialised, minimal system.
// All names carry the "sscrmcb_" prefix (sweep_separation_complex_rot_
// memtest_cluster_box) so they cannot collide with anything else in the suite.
// -----------------------------------------------------------------------------

/*! \brief Make sure the global GSL generator exists and is seeded reproducibly.
 *
 * resample_traj() (called indirectly by the function under test) uses GaussV(),
 * which requires the global generator `r` to be allocated.
 */
void sscrmcb_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}

/*! \brief Build a one-interface molecule whose interface sits on its COM.
 *
 * Keeping the interface coincident with the centre of mass means the
 * interface-interface separation used by the routine is exactly the
 * COM-COM separation, which makes the expected behaviour trivial to reason
 * about.
 */
Molecule sscrmcb_make_molecule(int index, int comIndex, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0; // single MolTemplate in every test
    mol.mass = 1.0;       // must be non-zero: update_properties() divides by it
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isEmpty = false;
    mol.trajStatus = TrajStatus::none; // "not yet propagated" => may be moved
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = com;
    iface.relIndex = 0;
    iface.index = 0;
    iface.molTypeIndex = 0;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a single-member complex holding molecule \p memberMolIndex.
 *
 * Rotational diffusion is deliberately set to zero so that trajRot stays
 * exactly (0,0,0); Complex::propagate() then takes its pure-translation fast
 * path, which makes the expected coordinates bit-for-bit predictable.
 */
Complex sscrmcb_make_complex(int index, int memberMolIndex, const Coord& com, double diffusion)
{
    Complex com_;
    com_.index = index;
    com_.id = index;
    com_.comCoord = com;
    com_.mass = 1.0;
    com_.radius = 1.0;
    com_.memberList.push_back(memberMolIndex);
    com_.numEachMol = std::vector<int>(1, 1);
    com_.lastNumberUpdateItrEachMol.resize(1);
    com_.D = Coord(diffusion, diffusion, diffusion);
    com_.Dr = Coord(0.0, 0.0, 0.0);
    com_.trajTrans.x = 0.0;
    com_.trajTrans.y = 0.0;
    com_.trajTrans.z = 0.0;
    com_.trajRot = Coord(0.0, 0.0, 0.0);
    com_.isEmpty = false;
    com_.OnSurface = false; // 3D complex => memtest == 0 => full 3D distance
    com_.onFiber = false;
    com_.linksToSurface = 0;
    com_.trajStatus = TrajStatus::none;
    return com_;
}

/*! \brief The one MolTemplate all molecules in these tests refer to. */
MolTemplate sscrmcb_make_template(double diffusion)
{
    MolTemplate temp;
    temp.molTypeIndex = 0;
    temp.molName = "A";
    temp.mass = 1.0;
    temp.radius = 1.0; // must be non-zero for the diffusion update
    temp.D = Coord(diffusion, diffusion, diffusion);
    temp.Dr = Coord(0.0, 0.0, 0.0);
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPoint = true;

    Interface iface("a", Coord(0.0, 0.0, 0.0));
    iface.index = 0;
    iface.stateList.emplace_back(std::string("a"), 0);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief One bimolecular reaction supplying the binding radius of every pair. */
ForwardRxn sscrmcb_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    // Both reactants use relative interface index 0, so the partner interface
    // index resolved inside cluster_one_complex() is 0 as well.
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    return rxn;
}

/*! \brief A reflecting cubic water box centred on the origin (side = 100 nm). */
Membrane sscrmcb_make_membrane()
{
    Membrane mem;
    mem.waterBox = Membrane::WaterBox(std::vector<double> { 100.0, 100.0, 100.0 });
    mem.isSphere = false;
    mem.isBox = true;
    mem.implicitLipid = false; // => RS3Dinput stays 0 and RS3Dvect is never read
    mem.hasCompartment = false;
    mem.xBCtype = "reflect";
    mem.yBCtype = "reflect";
    mem.zBCtype = "reflect";
    return mem;
}

/*! \brief Register molecule \p partner as a crossing partner of molecule \p mol.
 *
 * define_cluster_pairs() walks crossbase / mycrossint / crossrxn to discover
 * which pairs must be swept, so these three vectors must stay in lock-step.
 */
void sscrmcb_add_cross_partner(Molecule& mol, int partnerIndex)
{
    mol.crossbase.push_back(partnerIndex);
    mol.mycrossint.push_back(0);                     // our interface index
    mol.crossrxn.push_back(std::array<int, 3> { 0, 0, 0 }); // forward rxn 0
    mol.probvec.push_back(0.0);
}

/*! \brief Euclidean distance between the (single) interfaces of two molecules. */
double sscrmcb_iface_distance(const Molecule& m1, const Molecule& m2)
{
    const double dx = m1.interfaceList[0].coord.x - m2.interfaceList[0].coord.x;
    const double dy = m1.interfaceList[0].coord.y - m2.interfaceList[0].coord.y;
    const double dz = m1.interfaceList[0].coord.z - m2.interfaceList[0].coord.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a molecule with no crossing partners produces an empty cluster and the
//         routine must be a complete no-op (apart from restoring the timestep).
// -----------------------------------------------------------------------------
void test_sscrmcb_no_pairs_is_noop()
{
    std::cerr << "\n[TEST] test_sscrmcb_no_pairs_is_noop\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster_box.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster_box\n"
              << "  Scenario    : molecule 0 has an empty crossbase => no cluster pairs.\n"
              << "  Pass criteria: coordinates, trajStatus and trajTrans are all left\n"
              << "                 untouched and params.timeStep is unchanged.\n";

    sscrmcb_init_rng();

    Parameters params;
    params.timeStep = 0.1;
    const double origTimeStep = params.timeStep;

    Membrane membraneObject = sscrmcb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmcb_make_template(1.0e-4) };
    std::vector<ForwardRxn> forwardRxns { sscrmcb_make_rxn(1.0) };

    std::vector<Molecule> moleculeList { sscrmcb_make_molecule(0, 0, Coord(0.0, 0.0, 0.0)) };
    std::vector<Complex> complexList { sscrmcb_make_complex(0, 0, Coord(0.0, 0.0, 0.0), 1.0e-4) };

    // A non-zero trajectory that must NOT be consumed, because there is nothing
    // to sweep and therefore nothing to propagate.
    complexList[0].trajTrans.x = 0.25;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_box(pro1Index = 0)...\n";
    sweep_separation_complex_rot_memtest_cluster_box(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 0.0)
        << "With no cluster pairs the molecule must not be propagated";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.25)
        << "The unused trajectory vector must be preserved";
    EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::none)
        << "trajStatus must stay 'none' when no pair was processed";
    EXPECT_DOUBLE_EQ(params.timeStep, origTimeStep)
        << "params.timeStep must be restored to its original value";

    std::cerr << "  Result: comCoord.x = " << moleculeList[0].comCoord.x
              << ", trajTrans.x = " << complexList[0].trajTrans.x
              << ", timeStep = " << params.timeStep << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: a well separated pair keeps its sampled trajectory and is propagated
//         exactly by trajTrans (no resampling is needed).
// -----------------------------------------------------------------------------
void test_sscrmcb_separated_pair_propagates()
{
    std::cerr << "\n[TEST] test_sscrmcb_separated_pair_propagates\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster_box.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster_box\n"
              << "  Scenario    : two complexes 4 nm apart move 0.5 nm towards each other;\n"
              << "                the predicted separation (3 nm) exceeds the binding\n"
              << "                radius (1 nm), so no overlap correction is required.\n"
              << "  Pass criteria: each molecule is displaced by exactly its trajTrans,\n"
              << "                 trajStatus becomes 'propagated' and trajTrans is zeroed.\n";

    sscrmcb_init_rng();

    Parameters params;
    params.timeStep = 0.1;
    const double origTimeStep = params.timeStep;

    Membrane membraneObject = sscrmcb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmcb_make_template(1.0e-4) };
    std::vector<ForwardRxn> forwardRxns { sscrmcb_make_rxn(1.0) }; // bindRadius = 1 nm

    std::vector<Molecule> moleculeList {
        sscrmcb_make_molecule(0, 0, Coord(-2.0, 0.0, 0.0)),
        sscrmcb_make_molecule(1, 1, Coord(2.0, 0.0, 0.0))
    };
    std::vector<Complex> complexList {
        sscrmcb_make_complex(0, 0, Coord(-2.0, 0.0, 0.0), 1.0e-4),
        sscrmcb_make_complex(1, 1, Coord(2.0, 0.0, 0.0), 1.0e-4)
    };

    // Only molecule 0 lists a crossing partner, which yields exactly one pair.
    sscrmcb_add_cross_partner(moleculeList[0], 1);

    // Deterministic trajectories: pure translation towards each other.
    complexList[0].trajTrans.x = 0.5;
    complexList[1].trajTrans.x = -0.5;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_box(pro1Index = 0)...\n";
    sweep_separation_complex_rot_memtest_cluster_box(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // The trajectory should have been applied verbatim (no rotation, no resample).
    EXPECT_NEAR(moleculeList[0].comCoord.x, -1.5, 1e-12)
        << "Molecule 0 should have translated by exactly +0.5 nm";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 1.5, 1e-12)
        << "Molecule 1 should have translated by exactly -0.5 nm";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.x, -1.5, 1e-12)
        << "Interface of molecule 0 must move together with its COM";
    EXPECT_NEAR(moleculeList[1].interfaceList[0].coord.x, 1.5, 1e-12)
        << "Interface of molecule 1 must move together with its COM";
    EXPECT_NEAR(complexList[0].comCoord.x, -1.5, 1e-12)
        << "Complex 0 COM is recomputed from its members after propagation";

    // Book-keeping performed at the very end of the routine.
    EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::propagated)
        << "Molecule 0 must be flagged as propagated";
    EXPECT_EQ(moleculeList[1].trajStatus, TrajStatus::propagated)
        << "Molecule 1 must be flagged as propagated";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0)
        << "Complex 0 trajTrans must be zeroed after the sweep";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.x, 0.0)
        << "Complex 1 trajTrans must be zeroed after the sweep";
    EXPECT_DOUBLE_EQ(params.timeStep, origTimeStep)
        << "params.timeStep must be restored";

    std::cerr << "  Result: x0 = " << moleculeList[0].comCoord.x
              << ", x1 = " << moleculeList[1].comCoord.x
              << ", final interface separation = "
              << sscrmcb_iface_distance(moleculeList[0], moleculeList[1]) << '\n';
}

// -----------------------------------------------------------------------------
// Test 3: an overlapping pair triggers resampling; the huge, overlap-producing
//         displacement is thrown away and replaced by a fresh (tiny) sample.
// -----------------------------------------------------------------------------
void test_sscrmcb_overlapping_pair_is_resampled()
{
    std::cerr << "\n[TEST] test_sscrmcb_overlapping_pair_is_resampled\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster_box.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster_box\n"
              << "  Scenario    : two complexes 2 nm apart each try to move 0.95 nm\n"
              << "                towards the other, i.e. down to 0.1 nm which is well\n"
              << "                inside the 1 nm binding radius.  The diffusion constant\n"
              << "                is tiny so any resampled step is <0.02 nm.\n"
              << "  Pass criteria: the offending displacement is discarded (molecules stay\n"
              << "                 close to their starting positions) and the final\n"
              << "                 separation is no longer inside the binding radius.\n";

    sscrmcb_init_rng();

    Parameters params;
    params.timeStep = 0.1;
    const double origTimeStep = params.timeStep;

    const double bindRadius = 1.0;
    const double smallD = 1.0e-4; // sqrt(2*dt*D) ~ 0.0045 nm per axis

    Membrane membraneObject = sscrmcb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmcb_make_template(smallD) };
    std::vector<ForwardRxn> forwardRxns { sscrmcb_make_rxn(bindRadius) };

    std::vector<Molecule> moleculeList {
        sscrmcb_make_molecule(0, 0, Coord(-1.0, 0.0, 0.0)),
        sscrmcb_make_molecule(1, 1, Coord(1.0, 0.0, 0.0))
    };
    std::vector<Complex> complexList {
        sscrmcb_make_complex(0, 0, Coord(-1.0, 0.0, 0.0), smallD),
        sscrmcb_make_complex(1, 1, Coord(1.0, 0.0, 0.0), smallD)
    };

    sscrmcb_add_cross_partner(moleculeList[0], 1);

    // These two trajectories would put the interfaces 0.1 nm apart => overlap.
    complexList[0].trajTrans.x = 0.95;
    complexList[1].trajTrans.x = -0.95;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_box(pro1Index = 0)...\n";
    sweep_separation_complex_rot_memtest_cluster_box(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    const double finalSep = sscrmcb_iface_distance(moleculeList[0], moleculeList[1]);
    std::cerr << "  Result: x0 = " << moleculeList[0].comCoord.x
              << ", x1 = " << moleculeList[1].comCoord.x
              << ", final separation = " << finalSep
              << " (binding radius " << bindRadius << ")\n";

    EXPECT_GT(finalSep, bindRadius)
        << "After the sweep the pair must no longer be inside the binding radius";

    // The original 0.95 nm steps must have been discarded, so both molecules
    // remain within a few hundredths of a nm of where they started.
    EXPECT_NEAR(moleculeList[0].comCoord.x, -1.0, 0.1)
        << "Molecule 0 should have been resampled to a very small displacement";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 1.0, 0.1)
        << "Molecule 1 should have been resampled to a very small displacement";

    EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::propagated)
        << "Molecule 0 must be flagged as propagated";
    EXPECT_EQ(moleculeList[1].trajStatus, TrajStatus::propagated)
        << "Molecule 1 must be flagged as propagated";
    EXPECT_DOUBLE_EQ(params.timeStep, origTimeStep)
        << "params.timeStep must be restored";
}

// -----------------------------------------------------------------------------
// Test 4: pairs whose partner is an implicit lipid are skipped by the overlap
//         test, so no resampling happens even at "overlapping" separations.
// -----------------------------------------------------------------------------
void test_sscrmcb_implicit_lipid_partner_skipped()
{
    std::cerr << "\n[TEST] test_sscrmcb_implicit_lipid_partner_skipped\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster_box.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster_box\n"
              << "  Scenario    : the crossing partner is flagged isImplicitLipid, and the\n"
              << "                sampled trajectories would put the two interfaces only\n"
              << "                0.1 nm apart (inside the 1 nm binding radius).\n"
              << "  Pass criteria: the pair is skipped by the separation test, so the\n"
              << "                 original trajectory survives and is applied verbatim.\n";

    sscrmcb_init_rng();

    Parameters params;
    params.timeStep = 0.1;
    const double origTimeStep = params.timeStep;

    Membrane membraneObject = sscrmcb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmcb_make_template(1.0e-4) };
    std::vector<ForwardRxn> forwardRxns { sscrmcb_make_rxn(1.0) };

    std::vector<Molecule> moleculeList {
        sscrmcb_make_molecule(0, 0, Coord(-1.0, 0.0, 0.0)),
        sscrmcb_make_molecule(1, 1, Coord(1.0, 0.0, 0.0))
    };
    // Only the isImplicitLipid flag matters for the branch under test.
    moleculeList[1].isImplicitLipid = true;

    std::vector<Complex> complexList {
        sscrmcb_make_complex(0, 0, Coord(-1.0, 0.0, 0.0), 1.0e-4),
        sscrmcb_make_complex(1, 1, Coord(1.0, 0.0, 0.0), 1.0e-4)
    };

    sscrmcb_add_cross_partner(moleculeList[0], 1);

    // Would create a 0.1 nm separation - normally an overlap.
    complexList[0].trajTrans.x = 0.95;
    complexList[1].trajTrans.x = -0.95;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster_box(pro1Index = 0)...\n";
    sweep_separation_complex_rot_memtest_cluster_box(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Because the pair was skipped, the trajectory was applied unchanged.
    EXPECT_NEAR(moleculeList[0].comCoord.x, -0.05, 1e-12)
        << "Molecule 0 must keep its original +0.95 nm step (pair was skipped)";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 0.05, 1e-12)
        << "The implicit lipid partner is still propagated by its own trajectory";
    EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::propagated)
        << "Molecule 0 must be flagged as propagated";
    EXPECT_DOUBLE_EQ(params.timeStep, origTimeStep)
        << "params.timeStep must be restored";

    std::cerr << "  Result: x0 = " << moleculeList[0].comCoord.x
              << ", x1 = " << moleculeList[1].comCoord.x
              << ", final separation = "
              << sscrmcb_iface_distance(moleculeList[0], moleculeList[1]) << '\n';
}

// -----------------------------------------------------------------------------
// Test 5: more than 30 pairs switch the routine into its 10 sub-step mode.  The
//         key contract is that params.timeStep is divided *internally only* and
//         is restored exactly on exit, and that every molecule ends propagated.
// -----------------------------------------------------------------------------
void test_sscrmcb_large_cluster_substeps_restore_timestep()
{
    std::cerr << "\n[TEST] test_sscrmcb_large_cluster_substeps_restore_timestep\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster_box.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster_box\n"
              << "  Scenario    : one central molecule crosses 31 partners (>30 pairs),\n"
              << "                which activates the nStep = 10 sub-stepping branch.\n"
              << "  Pass criteria: params.timeStep is restored to its original value, all\n"
              << "                 involved molecules end as 'propagated', every complex\n"
              << "                 trajectory vector is zeroed, and displacements stay\n"
              << "                 small (tiny diffusion constant).\n";

    sscrmcb_init_rng();

    Parameters params;
    params.timeStep = 0.1;
    const double origTimeStep = params.timeStep;

    const int numPartners = 31; // > maxPairs (30) => sub-stepping is triggered
    const double smallD = 1.0e-4;

    Membrane membraneObject = sscrmcb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmcb_make_template(smallD) };
    std::vector<ForwardRxn> forwardRxns { sscrmcb_make_rxn(1.0) };

    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;

    // Central molecule at the origin.
    moleculeList.push_back(sscrmcb_make_molecule(0, 0, Coord(0.0, 0.0, 0.0)));
    complexList.push_back(sscrmcb_make_complex(0, 0, Coord(0.0, 0.0, 0.0), smallD));

    // Partners on a 20 nm circle: far from the centre (no overlap) and safely
    // inside the 100 nm reflecting box.
    for (int i = 0; i < numPartners; ++i) {
        const double angle = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(numPartners);
        const Coord pos(20.0 * std::cos(angle), 20.0 * std::sin(angle), 0.0);
        const int molIndex = i + 1;
        moleculeList.push_back(sscrmcb_make_molecule(molIndex, molIndex, pos));
        complexList.push_back(sscrmcb_make_complex(molIndex, molIndex, pos, smallD));
        sscrmcb_add_cross_partner(moleculeList[0], molIndex);
    }

    std::cerr << "  Built " << numPartners << " crossing partners; calling the sweep...\n";
    sweep_separation_complex_rot_memtest_cluster_box(
        0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    EXPECT_DOUBLE_EQ(params.timeStep, origTimeStep)
        << "The temporarily sub-divided timeStep must be restored exactly";

    // Every molecule of the cluster must be marked as propagated and every
    // complex trajectory must have been consumed (zeroed).
    for (std::size_t i = 0; i < moleculeList.size(); ++i) {
        EXPECT_EQ(moleculeList[i].trajStatus, TrajStatus::propagated)
            << "Molecule " << i << " should be flagged as propagated";
        EXPECT_DOUBLE_EQ(complexList[i].trajTrans.x, 0.0)
            << "Complex " << i << " trajTrans.x should be zeroed";
        EXPECT_DOUBLE_EQ(complexList[i].trajTrans.y, 0.0)
            << "Complex " << i << " trajTrans.y should be zeroed";
        EXPECT_DOUBLE_EQ(complexList[i].trajTrans.z, 0.0)
            << "Complex " << i << " trajTrans.z should be zeroed";
    }

    // With D = 1e-4 nm^2/us and 10 sub-steps of 0.01 us the accumulated random
    // walk is far below 1 nm, so the central molecule must stay near the origin.
    const double displacement = std::sqrt(moleculeList[0].comCoord.x * moleculeList[0].comCoord.x
        + moleculeList[0].comCoord.y * moleculeList[0].comCoord.y
        + moleculeList[0].comCoord.z * moleculeList[0].comCoord.z);
    std::cerr << "  Central molecule displacement after sub-stepping = " << displacement << " nm\n";
    EXPECT_LT(displacement, 1.0)
        << "Sub-stepped diffusion of a slow molecule must remain sub-nanometre";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario is reported individually, and none of them
// uses a fatal assertion so the whole set always runs to completion.
// -----------------------------------------------------------------------------
TEST(SweepSeparationComplexRotMemtestClusterBox, NoPairsIsNoop)
{
    test_sscrmcb_no_pairs_is_noop();
}

TEST(SweepSeparationComplexRotMemtestClusterBox, SeparatedPairPropagates)
{
    test_sscrmcb_separated_pair_propagates();
}

TEST(SweepSeparationComplexRotMemtestClusterBox, OverlappingPairIsResampled)
{
    test_sscrmcb_overlapping_pair_is_resampled();
}

TEST(SweepSeparationComplexRotMemtestClusterBox, ImplicitLipidPartnerSkipped)
{
    test_sscrmcb_implicit_lipid_partner_skipped();
}

TEST(SweepSeparationComplexRotMemtestClusterBox, LargeClusterSubstepsRestoreTimestep)
{
    test_sscrmcb_large_cluster_substeps_restore_timestep();
}