/*! \file test_sweep_separation_complex_rot_fiber_box.cpp
 *
 * ### Unit test for src/trajectory_functions/sweep_separation_complex_rot_fiber_box.cpp
 *
 * Function under test:
 *
 * \code
 * void sweep_separation_complex_rot_fiber_box(int simItr, int selfIndex,
 *          Parameters& params, std::vector<Molecule>& moleculeList,
 *          std::vector<Complex>& complexList,
 *          const std::vector<ForwardRxn>& forwardRxns,
 *          const std::vector<MolTemplate>& molTemplateList,
 *          const Membrane& membraneObject);
 * \endcode
 *
 * The routine sweeps the sampled displacement of the complex that owns
 * molecule `selfIndex` against every "crossed" partner molecule found during
 * the neighbour search of this timestep. If the proposed move would place a
 * reacting interface pair closer than the reaction binding radius (or, for
 * partners living on a fiber, would make the two objects hop past each other
 * along the fiber axis), the trajectory of the moving complex - and of the
 * offending partner complexes - is resampled. Finally the complex is always
 * propagated and its trajectory vectors are zeroed.
 *
 * ### Test strategy
 * Random resampling makes the routine non-deterministic in general, so every
 * scenario below is built with **zero diffusion constants**. With
 * `D == Dr == 0` the resampled trajectory `sqrt(2*dt*D)*GaussV()` is exactly
 * 0.0, which turns the resampling step into an observable, deterministic
 * event:
 *
 *   - if an overlap *is* detected, the initially non-zero trajectory is
 *     replaced by zero and the molecule ends the call where it started;
 *   - if no overlap is detected (or the pair is skipped), the molecule is
 *     propagated by exactly the trajectory we handed in.
 *
 * So "did the molecule move?" is a clean, exact yes/no probe of the overlap
 * logic. All coordinates are chosen so that no complex ever approaches a wall
 * of the (very large) reflecting simulation box.
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#include <gtest/gtest.h>

namespace {

// -----------------------------------------------------------------------------
// Container holding a complete, minimal simulation state for one scenario.
// -----------------------------------------------------------------------------
struct SscrfbWorld {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
};

/*! \brief (Re)seed the global GSL generator used by GaussV().
 *
 * `r` is defined in unit_tests/src/gtest_main.cpp; srand_gsl() must not be
 * used, so the generator is allocated and seeded explicitly here.
 */
void sscrfb_seed_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief Build a minimal molecule with a single interface sitting on its COM.
 *
 * Every field that Complex::update_properties() (called from
 * Complex::propagate()) reads is initialised: mass, molTypeIndex, isLipid,
 * isPromoter and linksToSurface.
 */
Molecule sscrfb_make_molecule(int index, int comIndex, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0; // the single MolTemplate of the world
    mol.mass = 1.0; // update_properties divides by the total mass
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    // One interface, coincident with the centre of mass: this makes the
    // interface-interface separation identical to the COM-COM separation.
    Molecule::Iface iface;
    iface.coord = com;
    iface.relIndex = 0;
    iface.index = 0;
    iface.molTypeIndex = 0;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a complex holding the listed member molecules.
 *
 * Diffusion constants are deliberately zero (see the file-level comment).
 */
Complex sscrfb_make_complex(int index, const Coord& com, const std::vector<int>& members)
{
    Complex cplx;
    cplx.index = index;
    cplx.id = index;
    cplx.comCoord = com;
    cplx.memberList = members;
    cplx.mass = static_cast<double>(members.size());
    cplx.radius = 1.0;
    cplx.isEmpty = false;
    cplx.OnSurface = false;
    cplx.onFiber = false;
    cplx.linksToSurface = 0;
    cplx.iLipidIndex = 0;
    cplx.ncross = 0;
    cplx.trajStatus = TrajStatus::none;

    // Frozen complex: every resample yields sqrt(2*dt*0)*GaussV() == 0.
    cplx.D = Coord { 0.0, 0.0, 0.0 };
    cplx.Dr = Coord { 0.0, 0.0, 0.0 };

    cplx.trajTrans = Vector { 0.0, 0.0, 0.0 };
    cplx.trajRot = Coord { 0.0, 0.0, 0.0 };

    cplx.numEachMol = std::vector<int>(1, static_cast<int>(members.size()));
    cplx.lastNumberUpdateItrEachMol = std::vector<long long int>(1, 0);

    return cplx;
}

/*! \brief Record that `mol` crossed paths with molecule `partnerIndex`.
 *
 * The three parallel vectors crossbase / mycrossint / crossrxn are what the
 * routine walks to rebuild the reacting interface pair.
 */
void sscrfb_link_partner(Molecule& mol, int partnerIndex)
{
    mol.crossbase.push_back(partnerIndex);
    mol.mycrossint.push_back(0); // our molecules only own relative interface 0
    mol.crossrxn.push_back(std::array<int, 3> { 0, 0, 0 }); // forwardRxns[0]
}

/*! \brief Assemble a world with two single-molecule complexes.
 *
 * Molecule 0 (complex 0) is the "self" molecule that will be swept; molecule 1
 * (complex 1) is its crossed partner. `bindRadius` sets the overlap threshold
 * used by the routine.
 */
SscrfbWorld sscrfb_build_world(const Coord& pos0, const Coord& pos1, double bindRadius)
{
    SscrfbWorld w;

    // --- parameters -------------------------------------------------------
    w.params.timeStep = 1.0;
    w.params.overlapSepLimit = 0.1;

    // --- boundary: a large reflecting cube, no compartment, no implicit lipid
    // (implicitLipid must stay false, otherwise the routine reads RS3Dvect).
    w.membrane.isSphere = false;
    w.membrane.implicitLipid = false;
    w.membrane.hasCompartment = false;
    w.membrane.waterBox = Membrane::WaterBox(std::vector<double> { 500.0, 500.0, 500.0 });
    w.membrane.xBCtype = "reflect";
    w.membrane.yBCtype = "reflect";
    w.membrane.zBCtype = "reflect";

    // --- one molecule template, with zero diffusion so that the diffusion
    //     constants recomputed inside propagate() stay zero as well ---------
    w.molTemplateList.emplace_back();
    MolTemplate& temp = w.molTemplateList.back();
    temp.molTypeIndex = 0;
    temp.molName = "A";
    temp.radius = 1.0;
    temp.mass = 1.0;
    temp.D = Coord { 0.0, 0.0, 0.0 };
    temp.Dr = Coord { 0.0, 0.0, 0.0 };
    temp.interfaceList.emplace_back(std::string("a"), Coord { 0.0, 0.0, 0.0 });
    temp.interfaceList.back().index = 0;

    // --- one bimolecular reaction supplying the binding radius -------------
    w.forwardRxns.emplace_back();
    ForwardRxn& rxn = w.forwardRxns.back();
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.reactantListNew.resize(2);
    rxn.reactantListNew[0].molTypeIndex = 0;
    rxn.reactantListNew[0].relIfaceIndex = 0;
    rxn.reactantListNew[1].molTypeIndex = 0;
    rxn.reactantListNew[1].relIfaceIndex = 0;

    // --- molecules and complexes ------------------------------------------
    w.moleculeList.push_back(sscrfb_make_molecule(0, 0, pos0));
    w.moleculeList.push_back(sscrfb_make_molecule(1, 1, pos1));

    // molecule 0 crossed molecule 1 during this timestep
    sscrfb_link_partner(w.moleculeList[0], 1);

    w.complexList.push_back(sscrfb_make_complex(0, pos0, std::vector<int> { 0 }));
    w.complexList.push_back(sscrfb_make_complex(1, pos1, std::vector<int> { 1 }));

    return w;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: partners far apart -> no overlap, the complex is propagated by
//         exactly the trajectory it was given and the trajectory is zeroed.
// -----------------------------------------------------------------------------
void test_sscrfb_no_overlap_propagates()
{
    std::cerr << "\n[TEST] test_sscrfb_no_overlap_propagates\n"
              << "  Source file: sweep_separation_complex_rot_fiber_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_fiber_box\n"
              << "  Scenario:    self at (0,0,0), partner at (10,0,0), sigma = 2.\n"
              << "               The proposed step keeps them ~9 nm apart.\n"
              << "  Pass:        molecule 0 lands exactly at COM + trajTrans and\n"
              << "               the complex trajectory vectors are zeroed.\n";

    sscrfb_seed_rng();
    SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 0.0 }, 2.0);

    // A small, purely translational displacement (no rotation keeps the
    // propagate() maths exact).
    w.complexList[0].trajTrans = Vector { 1.0, 0.5, -0.25 };
    w.complexList[0].trajRot = Coord { 0.0, 0.0, 0.0 };

    std::cerr << "  Calling sweep_separation_complex_rot_fiber_box(simItr=1, selfIndex=0)...\n";
    sweep_separation_complex_rot_fiber_box(1, 0, w.params, w.moleculeList, w.complexList,
        w.forwardRxns, w.molTemplateList, w.membrane);

    std::cerr << "  molecule 0 COM after sweep = (" << w.moleculeList[0].comCoord.x << ", "
              << w.moleculeList[0].comCoord.y << ", " << w.moleculeList[0].comCoord.z << ")\n";

    // The molecule (and its interface) must have moved by the whole step.
    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 1.0) << "x should advance by trajTrans.x";
    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.y, 0.5) << "y should advance by trajTrans.y";
    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.z, -0.25) << "z should advance by trajTrans.z";
    EXPECT_DOUBLE_EQ(w.moleculeList[0].interfaceList[0].coord.x, 1.0)
        << "the interface must be translated together with the COM";

    // The complex COM is recomputed from its (single) member.
    EXPECT_DOUBLE_EQ(w.complexList[0].comCoord.x, 1.0) << "complex COM should follow its member";

    // The routine always clears the trajectory before returning.
    EXPECT_DOUBLE_EQ(w.complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed on exit";
    EXPECT_DOUBLE_EQ(w.complexList[0].trajTrans.y, 0.0) << "trajTrans.y must be zeroed on exit";
    EXPECT_DOUBLE_EQ(w.complexList[0].trajTrans.z, 0.0) << "trajTrans.z must be zeroed on exit";
    EXPECT_DOUBLE_EQ(w.complexList[0].trajRot.x, 0.0) << "trajRot.x must be zeroed on exit";
    EXPECT_DOUBLE_EQ(w.complexList[0].trajRot.y, 0.0) << "trajRot.y must be zeroed on exit";
    EXPECT_DOUBLE_EQ(w.complexList[0].trajRot.z, 0.0) << "trajRot.z must be zeroed on exit";

    // The partner molecule is never propagated by this routine.
    EXPECT_DOUBLE_EQ(w.moleculeList[1].comCoord.x, 10.0) << "the partner must not be propagated";
}

// -----------------------------------------------------------------------------
// Test 2: genuine 3D overlap -> the moving complex AND the offending partner
//         complex both get their trajectories resampled (to zero here).
// -----------------------------------------------------------------------------
void test_sscrfb_overlap_resamples_self_and_partner()
{
    std::cerr << "\n[TEST] test_sscrfb_overlap_resamples_self_and_partner\n"
              << "  Source file: sweep_separation_complex_rot_fiber_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_fiber_box\n"
              << "  Scenario:    partner sits 1 nm away while sigma = 2 -> overlap.\n"
              << "               Both complexes have D = Dr = 0, so a resample sets\n"
              << "               every trajectory component to exactly 0.\n"
              << "  Pass:        molecule 0 does NOT move (its 5 nm step was thrown\n"
              << "               away) and the partner's pre-set trajectory is wiped.\n";

    sscrfb_seed_rng();
    SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 2.0);

    w.complexList[0].trajTrans = Vector { 5.0, 0.0, 0.0 }; // would be kept if no overlap
    w.complexList[1].trajTrans = Vector { 3.0, 3.0, 3.0 }; // must be resampled away
    w.moleculeList[1].trajStatus = TrajStatus::none; // partner is resamplable

    std::cerr << "  Calling sweep_separation_complex_rot_fiber_box(simItr=2, selfIndex=0)...\n";
    sweep_separation_complex_rot_fiber_box(2, 0, w.params, w.moleculeList, w.complexList,
        w.forwardRxns, w.molTemplateList, w.membrane);

    std::cerr << "  molecule 0 COM after sweep = (" << w.moleculeList[0].comCoord.x << ", "
              << w.moleculeList[0].comCoord.y << ", " << w.moleculeList[0].comCoord.z << ")\n"
              << "  partner complex trajTrans  = (" << w.complexList[1].trajTrans.x << ", "
              << w.complexList[1].trajTrans.y << ", " << w.complexList[1].trajTrans.z << ")\n";

    // Overlap detected -> the 5 nm step was replaced by the frozen resample.
    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 0.0)
        << "the overlapping step must have been discarded (x unchanged)";
    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.y, 0.0) << "y unchanged";
    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.z, 0.0) << "z unchanged";

    // The partner complex is resampled in place (never propagated here).
    EXPECT_DOUBLE_EQ(w.complexList[1].trajTrans.x, 0.0)
        << "partner trajTrans.x must have been resampled";
    EXPECT_DOUBLE_EQ(w.complexList[1].trajTrans.y, 0.0)
        << "partner trajTrans.y must have been resampled";
    EXPECT_DOUBLE_EQ(w.complexList[1].trajTrans.z, 0.0)
        << "partner trajTrans.z must have been resampled";
    EXPECT_DOUBLE_EQ(w.moleculeList[1].comCoord.x, 1.0)
        << "the partner molecule itself is not propagated by this routine";
}

// -----------------------------------------------------------------------------
// Test 3: a partner whose trajStatus is `propagated` has already moved this
//         step and must therefore not be resampled.
// -----------------------------------------------------------------------------
void test_sscrfb_already_propagated_partner_not_resampled()
{
    std::cerr << "\n[TEST] test_sscrfb_already_propagated_partner_not_resampled\n"
              << "  Source file: sweep_separation_complex_rot_fiber_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_fiber_box\n"
              << "  Scenario:    overlapping partner carries TrajStatus::propagated.\n"
              << "  Pass:        the partner keeps its (3,3,3) trajectory, while the\n"
              << "               swept complex is still frozen by the overlap.\n";

    sscrfb_seed_rng();
    SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 2.0);

    w.complexList[0].trajTrans = Vector { 5.0, 0.0, 0.0 };
    w.complexList[1].trajTrans = Vector { 3.0, 3.0, 3.0 };
    w.moleculeList[1].trajStatus = TrajStatus::propagated; // already moved this step

    std::cerr << "  Calling sweep_separation_complex_rot_fiber_box(simItr=3, selfIndex=0)...\n";
    sweep_separation_complex_rot_fiber_box(3, 0, w.params, w.moleculeList, w.complexList,
        w.forwardRxns, w.molTemplateList, w.membrane);

    std::cerr << "  partner complex trajTrans = (" << w.complexList[1].trajTrans.x << ", "
              << w.complexList[1].trajTrans.y << ", " << w.complexList[1].trajTrans.z << ")\n";

    EXPECT_DOUBLE_EQ(w.complexList[1].trajTrans.x, 3.0)
        << "an already-propagated partner must keep its trajectory (x)";
    EXPECT_DOUBLE_EQ(w.complexList[1].trajTrans.y, 3.0)
        << "an already-propagated partner must keep its trajectory (y)";
    EXPECT_DOUBLE_EQ(w.complexList[1].trajTrans.z, 3.0)
        << "an already-propagated partner must keep its trajectory (z)";

    // The swept complex still had its own overlapping step rejected.
    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 0.0)
        << "the swept complex is still frozen by the unresolved overlap";
}

// -----------------------------------------------------------------------------
// Test 4: implicit-lipid partners are skipped outright.
// -----------------------------------------------------------------------------
void test_sscrfb_implicit_lipid_partner_skipped()
{
    std::cerr << "\n[TEST] test_sscrfb_implicit_lipid_partner_skipped\n"
              << "  Source file: sweep_separation_complex_rot_fiber_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_fiber_box\n"
              << "  Scenario:    the (geometrically overlapping) partner is flagged\n"
              << "               isImplicitLipid.\n"
              << "  Pass:        the pair is skipped, so the full step is kept.\n";

    sscrfb_seed_rng();
    SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 2.0);

    w.moleculeList[1].isImplicitLipid = true; // must be ignored by the sweep
    w.complexList[0].trajTrans = Vector { 5.0, 0.0, 0.0 };

    std::cerr << "  Calling sweep_separation_complex_rot_fiber_box(simItr=4, selfIndex=0)...\n";
    sweep_separation_complex_rot_fiber_box(4, 0, w.params, w.moleculeList, w.complexList,
        w.forwardRxns, w.molTemplateList, w.membrane);

    std::cerr << "  molecule 0 x after sweep = " << w.moleculeList[0].comCoord.x << " (expected 5)\n";
    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 5.0)
        << "an implicit-lipid partner must never trigger a resample";
}

// -----------------------------------------------------------------------------
// Test 5: mixed promoter / non-promoter pairs are skipped in both orders.
// -----------------------------------------------------------------------------
void test_sscrfb_promoter_mismatch_skipped()
{
    std::cerr << "\n[TEST] test_sscrfb_promoter_mismatch_skipped\n"
              << "  Source file: sweep_separation_complex_rot_fiber_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_fiber_box\n"
              << "  Scenario:    exactly one of the two overlapping molecules is a\n"
              << "               promoter (co-localized 1D object).\n"
              << "  Pass:        the pair is skipped in BOTH orderings, so the swept\n"
              << "               complex keeps its whole 5 nm step.\n";

    // --- case A: self is the promoter, partner is not ----------------------
    {
        sscrfb_seed_rng();
        SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 2.0);
        w.moleculeList[0].isPromoter = true;
        w.moleculeList[1].isPromoter = false;
        w.complexList[0].trajTrans = Vector { 5.0, 0.0, 0.0 };

        std::cerr << "  case A (self promoter, partner not): calling sweep...\n";
        sweep_separation_complex_rot_fiber_box(5, 0, w.params, w.moleculeList, w.complexList,
            w.forwardRxns, w.molTemplateList, w.membrane);

        std::cerr << "    molecule 0 x = " << w.moleculeList[0].comCoord.x << " (expected 5)\n";
        EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 5.0)
            << "promoter vs non-promoter must be skipped (self is the promoter)";
    }

    // --- case B: partner is the promoter, self is not ----------------------
    {
        sscrfb_seed_rng();
        SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 2.0);
        w.moleculeList[0].isPromoter = false;
        w.moleculeList[1].isPromoter = true;
        w.complexList[0].trajTrans = Vector { 5.0, 0.0, 0.0 };

        std::cerr << "  case B (partner promoter, self not): calling sweep...\n";
        sweep_separation_complex_rot_fiber_box(5, 0, w.params, w.moleculeList, w.complexList,
            w.forwardRxns, w.molTemplateList, w.membrane);

        std::cerr << "    molecule 0 x = " << w.moleculeList[0].comCoord.x << " (expected 5)\n";
        EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 5.0)
            << "promoter vs non-promoter must be skipped (partner is the promoter)";
    }
}

// -----------------------------------------------------------------------------
// Test 6: partners inside the *same* complex cannot move relative to each
//         other and must therefore be ignored.
// -----------------------------------------------------------------------------
void test_sscrfb_same_complex_partner_ignored()
{
    std::cerr << "\n[TEST] test_sscrfb_same_complex_partner_ignored\n"
              << "  Source file: sweep_separation_complex_rot_fiber_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_fiber_box\n"
              << "  Scenario:    both crossed molecules belong to complex 0 and sit\n"
              << "               1 nm apart (inside sigma = 2).\n"
              << "  Pass:        no resample happens and both members translate by\n"
              << "               the full trajectory.\n";

    sscrfb_seed_rng();
    SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }, 2.0);

    // Move molecule 1 into complex 0 and make complex 1 an empty shell.
    w.moleculeList[1].myComIndex = 0;
    w.complexList[0].memberList = std::vector<int> { 0, 1 };
    w.complexList[0].mass = 2.0;
    w.complexList[1].memberList.clear();
    w.complexList[1].isEmpty = true;

    // Cross-link the two members both ways, as the neighbour search would.
    sscrfb_link_partner(w.moleculeList[1], 0);

    w.complexList[0].trajTrans = Vector { 5.0, 0.0, 0.0 };

    std::cerr << "  Calling sweep_separation_complex_rot_fiber_box(simItr=6, selfIndex=0)...\n";
    sweep_separation_complex_rot_fiber_box(6, 0, w.params, w.moleculeList, w.complexList,
        w.forwardRxns, w.molTemplateList, w.membrane);

    std::cerr << "  member x positions after sweep = " << w.moleculeList[0].comCoord.x << ", "
              << w.moleculeList[1].comCoord.x << " (expected 5 and 6)\n";

    EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 5.0)
        << "intra-complex neighbours must not veto the move (member 0)";
    EXPECT_DOUBLE_EQ(w.moleculeList[1].comCoord.x, 6.0)
        << "intra-complex neighbours must not veto the move (member 1)";
    // The recomputed, mass-weighted complex COM is the midpoint of both members.
    EXPECT_DOUBLE_EQ(w.complexList[0].comCoord.x, 5.5)
        << "the complex COM should be the mass-weighted mean of its members";
}

// -----------------------------------------------------------------------------
// Test 7: for a partner living on a fiber only the x separation counts.
// -----------------------------------------------------------------------------
void test_sscrfb_fiber_uses_only_x_separation()
{
    std::cerr << "\n[TEST] test_sscrfb_fiber_uses_only_x_separation\n"
              << "  Source file: sweep_separation_complex_rot_fiber_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_fiber_box\n"
              << "  Scenario:    partner at (1, 50, 50), sigma = 2. In 3D the pair is\n"
              << "               ~70 nm apart, along the fiber axis only 1 nm apart.\n"
              << "  Pass:        onFiber partner  -> overlap detected -> step rejected;\n"
              << "               off-fiber partner -> no overlap      -> step kept.\n";

    // --- case A: partner complex flagged onFiber (1D distance is used) -----
    {
        sscrfb_seed_rng();
        SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 50.0, 50.0 }, 2.0);
        w.complexList[0].onFiber = true;
        w.complexList[1].onFiber = true; // this is what memCheckList reads
        w.complexList[0].trajTrans = Vector { 0.0, 0.0, 0.0 };
        // A tiny step that neither separates them nor makes them hop.
        w.complexList[0].trajTrans.y = 4.0;

        std::cerr << "  case A (partner onFiber): calling sweep...\n";
        sweep_separation_complex_rot_fiber_box(7, 0, w.params, w.moleculeList, w.complexList,
            w.forwardRxns, w.molTemplateList, w.membrane);

        std::cerr << "    molecule 0 y = " << w.moleculeList[0].comCoord.y << " (expected 0)\n";
        EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.y, 0.0)
            << "on a fiber the 1 nm x-separation is inside sigma -> step rejected";
    }

    // --- case B: identical geometry but the partner is a 3D object ---------
    {
        sscrfb_seed_rng();
        SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 50.0, 50.0 }, 2.0);
        w.complexList[0].onFiber = true;
        w.complexList[1].onFiber = false; // full 3D separation is used
        w.complexList[0].trajTrans = Vector { 0.0, 4.0, 0.0 };

        std::cerr << "  case B (partner off fiber): calling sweep...\n";
        sweep_separation_complex_rot_fiber_box(7, 0, w.params, w.moleculeList, w.complexList,
            w.forwardRxns, w.molTemplateList, w.membrane);

        std::cerr << "    molecule 0 y = " << w.moleculeList[0].comCoord.y << " (expected 4)\n";
        EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.y, 4.0)
            << "in 3D the pair is far apart -> the whole step must be kept";
    }
}

// -----------------------------------------------------------------------------
// Test 8: hop detection along the fiber (only for non-promoter pairs).
// -----------------------------------------------------------------------------
void test_sscrfb_fiber_hop_detection()
{
    std::cerr << "\n[TEST] test_sscrfb_fiber_hop_detection\n"
              << "  Source file: sweep_separation_complex_rot_fiber_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_fiber_box\n"
              << "  Scenario:    self at x = 0, fiber partner at x = 10, step = +20\n"
              << "               so the two would swap sides (a 'hop'), even though\n"
              << "               the final separation (10 nm) exceeds sigma = 2.\n"
              << "  Pass:        ordinary fiber pair -> hop rejected (no motion);\n"
              << "               promoter/promoter pair -> hop allowed (moves 20 nm).\n";

    // --- case A: two ordinary (non-promoter) objects on the fiber ----------
    {
        sscrfb_seed_rng();
        SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 0.0 }, 2.0);
        w.complexList[0].onFiber = true;
        w.complexList[1].onFiber = true;
        w.complexList[0].trajTrans = Vector { 20.0, 0.0, 0.0 };

        std::cerr << "  case A (non-promoter fiber pair): calling sweep...\n";
        sweep_separation_complex_rot_fiber_box(8, 0, w.params, w.moleculeList, w.complexList,
            w.forwardRxns, w.molTemplateList, w.membrane);

        std::cerr << "    molecule 0 x = " << w.moleculeList[0].comCoord.x << " (expected 0)\n";
        EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 0.0)
            << "a sign change of the x separation is a hop and must be rejected";
    }

    // --- case B: both molecules are promoters (DNA sites) ------------------
    // For co-localized DNA sites hopping is legal; only true overlap counts.
    {
        sscrfb_seed_rng();
        SscrfbWorld w = sscrfb_build_world(Coord { 0.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 0.0 }, 2.0);
        w.complexList[0].onFiber = true;
        w.complexList[1].onFiber = true;
        w.moleculeList[0].isPromoter = true;
        w.moleculeList[1].isPromoter = true;
        w.complexList[0].trajTrans = Vector { 20.0, 0.0, 0.0 };

        std::cerr << "  case B (promoter/promoter fiber pair): calling sweep...\n";
        sweep_separation_complex_rot_fiber_box(8, 0, w.params, w.moleculeList, w.complexList,
            w.forwardRxns, w.molTemplateList, w.membrane);

        std::cerr << "    molecule 0 x = " << w.moleculeList[0].comCoord.x << " (expected 20)\n";
        EXPECT_DOUBLE_EQ(w.moleculeList[0].comCoord.x, 20.0)
            << "DNA-site pairs are only checked for overlap, hops are allowed";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario is reported individually but all of them
// run even if an earlier one fails (only EXPECT_* assertions are used).
// -----------------------------------------------------------------------------
TEST(SweepSeparationComplexRotFiberBox, NoOverlapPropagates) { test_sscrfb_no_overlap_propagates(); }
TEST(SweepSeparationComplexRotFiberBox, OverlapResamplesSelfAndPartner) { test_sscrfb_overlap_resamples_self_and_partner(); }
TEST(SweepSeparationComplexRotFiberBox, AlreadyPropagatedPartnerNotResampled) { test_sscrfb_already_propagated_partner_not_resampled(); }
TEST(SweepSeparationComplexRotFiberBox, ImplicitLipidPartnerSkipped) { test_sscrfb_implicit_lipid_partner_skipped(); }
TEST(SweepSeparationComplexRotFiberBox, PromoterMismatchSkipped) { test_sscrfb_promoter_mismatch_skipped(); }
TEST(SweepSeparationComplexRotFiberBox, SameComplexPartnerIgnored) { test_sscrfb_same_complex_partner_ignored(); }
TEST(SweepSeparationComplexRotFiberBox, FiberUsesOnlyXSeparation) { test_sscrfb_fiber_uses_only_x_separation(); }
TEST(SweepSeparationComplexRotFiberBox, FiberHopDetection) { test_sscrfb_fiber_hop_detection(); }