/*! \file test_sweep_separation_complex_rot_memtest_sphere.cpp
 *
 * ### Unit test for src/trajectory_functions/sweep_separation_complex_rot_memtest_sphere.cpp
 *
 * The single function under test is
 *
 *     void sweep_separation_complex_rot_memtest_sphere(int simItr, int pro1Index,
 *              Parameters& params, std::vector<Molecule>& moleculeList,
 *              std::vector<Complex>& complexList,
 *              const std::vector<ForwardRxn>& forwardRxns,
 *              const std::vector<MolTemplate>& molTemplateList,
 *              const Membrane& membraneObject);
 *
 * The routine sweeps the "cross" (encounter) list of every molecule inside the
 * complex that owns `pro1Index`.  For each encountered partner it predicts the
 * interface-interface separation *after* the currently sampled trajectory would
 * be applied.  If that separation is smaller than the reaction binding radius,
 * the trajectories of the target complex (and of eligible partner complexes)
 * are re-sampled, up to ten times.  Finally the target complex is propagated and
 * its stored translation/rotation vectors are zeroed.
 *
 * The branches exercised here are:
 *
 *   1. no cross partners at all            -> trajectory applied verbatim
 *   2. partner inside the *same* complex   -> pair is skipped (no re-sampling)
 *   3. partner is an implicit lipid        -> pair is skipped (no re-sampling)
 *   4. genuine overlap in another complex  -> both trajectories re-sampled
 *   5. partner complex flagged OnSurface   -> spherical-surface code path
 *
 * All console chatter is written to stderr so the reader can follow exactly what
 * is being exercised and what the pass criteria are.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

// `r` (the GSL random number generator) and `totMatches` are defined in
// unit_tests/src/gtest_main.cpp; rand_gsl.hpp already declares `r` as extern.

namespace {

// -----------------------------------------------------------------------------
// Small helpers used to build a completely initialised, minimal simulation
// state.  Every field that the code under test (or anything it calls) reads is
// filled in explicitly -- an under-initialised Molecule/Complex would segfault
// the whole test binary.
// -----------------------------------------------------------------------------

/*! \brief Initialise (once) and re-seed the global GSL generator.
 *
 * GaussV()/rand_gsl() are used whenever a trajectory has to be re-sampled, so
 * the generator must exist before the function under test is called.
 */
void sscrms_seed_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, 42);
}

/*! \brief Build a point-like molecule with exactly one interface.
 *
 * \param[in] index    index of this molecule inside moleculeList
 * \param[in] comIndex index of the owning complex inside complexList
 * \param[in] com      center of mass coordinate
 * \param[in] ifaceCrd absolute coordinate of the single interface
 */
Molecule sscrms_make_molecule(int index, int comIndex, const Coord& com, const Coord& ifaceCrd)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0; // must be non-zero: Complex::update_properties divides by total mass
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = ifaceCrd;
    iface.relIndex = 0;
    iface.index = 0;
    iface.stateIndex = 0;
    iface.stateIden = 'A';
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // no encounters by default
    mol.crossbase.clear();
    mol.mycrossint.clear();
    mol.crossrxn.clear();
    return mol;
}

/*! \brief Register molecule `partnerIndex` as an encounter of `mol`.
 *
 * The relative interface index of `mol` is 0 (its only interface) and the
 * reaction used is reaction 0 of the forwardRxns list.
 */
void sscrms_add_cross_partner(Molecule& mol, int partnerIndex)
{
    mol.crossbase.push_back(partnerIndex);
    mol.mycrossint.push_back(0);
    mol.crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
}

/*! \brief Build a complex containing the listed member molecules. */
Complex sscrms_make_complex(int index, const Coord& com, const std::vector<int>& members,
    bool onSurface, double dTrans)
{
    Complex targCom;
    targCom.index = index;
    targCom.id = index;
    targCom.comCoord = com;
    targCom.memberList = members;
    targCom.numEachMol = std::vector<int>(1, static_cast<int>(members.size()));
    targCom.lastNumberUpdateItrEachMol = std::vector<long long int>(1, 0);
    targCom.mass = static_cast<double>(members.size());
    targCom.radius = 1.0;

    // Translational diffusion is finite; rotational diffusion is intentionally
    // zero so that any re-sampled rotation is exactly zero.  That keeps the
    // propagation step a pure translation and therefore easy to predict.
    targCom.D.x = dTrans;
    targCom.D.y = dTrans;
    targCom.D.z = onSurface ? 0.0 : dTrans;
    targCom.Dr.x = 0.0;
    targCom.Dr.y = 0.0;
    targCom.Dr.z = 0.0;

    targCom.isEmpty = false;
    targCom.OnSurface = onSurface;
    targCom.onFiber = false;
    targCom.tmpOnSurface = false;
    targCom.linksToSurface = 0;
    targCom.iLipidIndex = 0;
    targCom.ncross = 0;
    targCom.trajStatus = TrajStatus::none;

    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot.x = 0.0;
    targCom.trajRot.y = 0.0;
    targCom.trajRot.z = 0.0;
    return targCom;
}

/*! \brief Minimal MolTemplate; only radius / D / Dr are read by update_properties. */
MolTemplate sscrms_make_template()
{
    MolTemplate temp;
    temp.molTypeIndex = 0;
    temp.molName = "A";
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.copies = 0;
    temp.D = Coord { 5.0, 5.0, 5.0 };
    temp.Dr = Coord { 0.0, 0.0, 0.0 };
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isRod = false;
    temp.isPoint = true;
    temp.isPromoter = false;

    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord { 0.0, 0.0, 0.0 };
    iface.stateList.emplace_back('A', 0);
    temp.interfaceList.push_back(iface);
    return temp;
}

/*! \brief One bimolecular reaction whose two reactants are both interface 0. */
ForwardRxn sscrms_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = false;
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    return rxn;
}

/*! \brief Spherical membrane.
 *
 * NOTE: the routine unconditionally reads membraneObject.RS3Dvect at indices
 * 300..499, so the vector must be at least 500 elements long even when implicit
 * lipids are not in use.
 */
Membrane sscrms_make_sphere_membrane(double radius)
{
    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.isBox = false;
    membraneObject.implicitLipid = false;
    membraneObject.sphereR = radius;
    membraneObject.sphereVol = (4.0 / 3.0) * M_PI * radius * radius * radius;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 2 * radius, 2 * radius, 2 * radius });
    membraneObject.RS3Dvect.assign(500, 0.0);
    membraneObject.xBCtype = "reflect";
    membraneObject.yBCtype = "reflect";
    membraneObject.zBCtype = "reflect";
    return membraneObject;
}

/*! \brief Convenience: distance between two coordinates. */
double sscrms_distance(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: a complex whose only molecule has an empty cross list.
//         Nothing can overlap, therefore the sampled trajectory must be applied
//         exactly as it was, and afterwards the stored vectors must be zeroed.
// -----------------------------------------------------------------------------
void test_sscrms_no_cross_partners()
{
    std::cerr << "\n[TEST] test_sscrms_no_cross_partners\n"
              << "  Source file:   sweep_separation_complex_rot_memtest_sphere.cpp\n"
              << "  Function:      sweep_separation_complex_rot_memtest_sphere\n"
              << "  Scenario:      single molecule / single complex, no encounters.\n"
              << "  Pass criteria: molecule is displaced by exactly trajTrans,\n"
              << "                 complex trajTrans/trajRot are zeroed and both\n"
              << "                 molecule and complex are flagged 'propagated'.\n";

    sscrms_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrms_make_sphere_membrane(100.0);
    std::vector<MolTemplate> molTemplateList { sscrms_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrms_make_rxn(1.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrms_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(sscrms_make_complex(0, Coord { 0.0, 0.0, 0.0 }, std::vector<int> { 0 }, false, 5.0));

    // A deterministic displacement; no rotation so propagate() takes its exact
    // translation-only fast path.
    complexList[0].trajTrans.x = 1.0;
    complexList[0].trajTrans.y = 2.0;
    complexList[0].trajTrans.z = -3.0;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_sphere(simItr=0, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_sphere(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);

    std::cerr << "  Molecule 0 now at (" << moleculeList[0].comCoord.x << ", "
              << moleculeList[0].comCoord.y << ", " << moleculeList[0].comCoord.z << ")\n";

    // The trajectory must have been applied verbatim (no overlap = no re-sampling).
    EXPECT_NEAR(moleculeList[0].comCoord.x, 1.0, 1e-12) << "x should be old x + trajTrans.x";
    EXPECT_NEAR(moleculeList[0].comCoord.y, 2.0, 1e-12) << "y should be old y + trajTrans.y";
    EXPECT_NEAR(moleculeList[0].comCoord.z, -3.0, 1e-12) << "z should be old z + trajTrans.z";

    // The complex center of mass is recomputed from its (single) member.
    EXPECT_NEAR(complexList[0].comCoord.x, moleculeList[0].comCoord.x, 1e-12)
        << "complex COM should track its only member";
    EXPECT_NEAR(complexList[0].comCoord.y, moleculeList[0].comCoord.y, 1e-12)
        << "complex COM should track its only member";
    EXPECT_NEAR(complexList[0].comCoord.z, moleculeList[0].comCoord.z, 1e-12)
        << "complex COM should track its only member";

    // Trajectory storage is cleared at the very end of the routine.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "trajTrans.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "trajTrans.z must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.x, 0.0) << "trajRot.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.y, 0.0) << "trajRot.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.z, 0.0) << "trajRot.z must be zeroed";

    // propagate() marks everything it moved.
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "molecule should be flagged as propagated";
    EXPECT_EQ(static_cast<int>(complexList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "complex should be flagged as propagated";
}

// -----------------------------------------------------------------------------
// Test 2: the encountered partner lives inside the *same* complex.  Members of a
//         complex cannot move relative to one another, so the pair must be
//         skipped even though the interfaces are much closer than the binding
//         radius.  Detection is by checking that the trajectory was NOT
//         re-sampled (both molecules move by exactly the original trajTrans).
// -----------------------------------------------------------------------------
void test_sscrms_same_complex_partner_is_skipped()
{
    std::cerr << "\n[TEST] test_sscrms_same_complex_partner_is_skipped\n"
              << "  Source file:   sweep_separation_complex_rot_memtest_sphere.cpp\n"
              << "  Function:      sweep_separation_complex_rot_memtest_sphere\n"
              << "  Scenario:      two molecules of the SAME complex sit 0.1 nm apart\n"
              << "                 while the binding radius is 1.0 nm.\n"
              << "  Pass criteria: no re-sampling occurs, i.e. both members are\n"
              << "                 translated by exactly the original trajTrans.\n";

    sscrms_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrms_make_sphere_membrane(100.0);
    std::vector<MolTemplate> molTemplateList { sscrms_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrms_make_rxn(1.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrms_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrms_make_molecule(1, 0, Coord { 0.1, 0.0, 0.0 }, Coord { 0.1, 0.0, 0.0 }));

    // They "see" each other -- but they are in the same complex.
    sscrms_add_cross_partner(moleculeList[0], 1);
    sscrms_add_cross_partner(moleculeList[1], 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrms_make_complex(0, Coord { 0.05, 0.0, 0.0 },
        std::vector<int> { 0, 1 }, false, 5.0));

    complexList[0].trajTrans.x = 2.0;
    complexList[0].trajTrans.y = -1.0;
    complexList[0].trajTrans.z = 0.5;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_sphere(simItr=1, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_sphere(1, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);

    std::cerr << "  Molecule 0 -> (" << moleculeList[0].comCoord.x << ", " << moleculeList[0].comCoord.y
              << ", " << moleculeList[0].comCoord.z << ")\n";
    std::cerr << "  Molecule 1 -> (" << moleculeList[1].comCoord.x << ", " << moleculeList[1].comCoord.y
              << ", " << moleculeList[1].comCoord.z << ")\n";

    // Original trajectory (2, -1, 0.5) must have been applied unchanged.
    EXPECT_NEAR(moleculeList[0].comCoord.x, 2.0, 1e-12) << "member 0 x moved by the original trajTrans";
    EXPECT_NEAR(moleculeList[0].comCoord.y, -1.0, 1e-12) << "member 0 y moved by the original trajTrans";
    EXPECT_NEAR(moleculeList[0].comCoord.z, 0.5, 1e-12) << "member 0 z moved by the original trajTrans";

    EXPECT_NEAR(moleculeList[1].comCoord.x, 2.1, 1e-12) << "member 1 x moved by the original trajTrans";
    EXPECT_NEAR(moleculeList[1].comCoord.y, -1.0, 1e-12) << "member 1 y moved by the original trajTrans";
    EXPECT_NEAR(moleculeList[1].comCoord.z, 0.5, 1e-12) << "member 1 z moved by the original trajTrans";

    // Their relative geometry inside the rigid complex is untouched.
    EXPECT_NEAR(sscrms_distance(moleculeList[0].comCoord, moleculeList[1].comCoord), 0.1, 1e-12)
        << "intra-complex separation must be conserved";
}

// -----------------------------------------------------------------------------
// Test 3: the encountered partner is an implicit lipid.  The routine `continue`s
//         before it ever computes a separation, so the trajectory must again be
//         applied verbatim and the partner complex must be left untouched.
// -----------------------------------------------------------------------------
void test_sscrms_implicit_lipid_partner_is_skipped()
{
    std::cerr << "\n[TEST] test_sscrms_implicit_lipid_partner_is_skipped\n"
              << "  Source file:   sweep_separation_complex_rot_memtest_sphere.cpp\n"
              << "  Function:      sweep_separation_complex_rot_memtest_sphere\n"
              << "  Scenario:      the cross partner (in a different complex) is\n"
              << "                 flagged isImplicitLipid and sits 0.1 nm away.\n"
              << "  Pass criteria: no overlap is registered -- target moves by the\n"
              << "                 original trajTrans, partner keeps zero trajTrans\n"
              << "                 and does not move.\n";

    sscrms_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrms_make_sphere_membrane(100.0);
    std::vector<MolTemplate> molTemplateList { sscrms_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrms_make_rxn(1.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrms_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrms_make_molecule(1, 1, Coord { 0.1, 0.0, 0.0 }, Coord { 0.1, 0.0, 0.0 }));
    moleculeList[1].isImplicitLipid = true; // <- the branch under test

    sscrms_add_cross_partner(moleculeList[0], 1);

    std::vector<Complex> complexList;
    complexList.push_back(sscrms_make_complex(0, Coord { 0.0, 0.0, 0.0 }, std::vector<int> { 0 }, false, 5.0));
    complexList.push_back(sscrms_make_complex(1, Coord { 0.1, 0.0, 0.0 }, std::vector<int> { 1 }, false, 5.0));

    complexList[0].trajTrans.x = 1.5;
    complexList[0].trajTrans.y = 0.25;
    complexList[0].trajTrans.z = -0.75;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_sphere(simItr=2, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_sphere(2, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);

    std::cerr << "  Molecule 0 -> (" << moleculeList[0].comCoord.x << ", " << moleculeList[0].comCoord.y
              << ", " << moleculeList[0].comCoord.z << ")\n";

    EXPECT_NEAR(moleculeList[0].comCoord.x, 1.5, 1e-12) << "target x moved by the original trajTrans";
    EXPECT_NEAR(moleculeList[0].comCoord.y, 0.25, 1e-12) << "target y moved by the original trajTrans";
    EXPECT_NEAR(moleculeList[0].comCoord.z, -0.75, 1e-12) << "target z moved by the original trajTrans";

    // The implicit lipid partner must not have been re-sampled or moved.
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.x, 0.0) << "partner trajTrans.x must stay zero";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.y, 0.0) << "partner trajTrans.y must stay zero";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.z, 0.0) << "partner trajTrans.z must stay zero";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 0.1, 1e-12) << "partner molecule must not move";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::none))
        << "partner molecule must not be propagated by this call";
}

// -----------------------------------------------------------------------------
// Test 4: a genuine overlap with a partner in a *different* complex.  Both
//         trajectories must be re-sampled; only the target complex is
//         propagated at the end.
// -----------------------------------------------------------------------------
void test_sscrms_overlap_triggers_resampling()
{
    std::cerr << "\n[TEST] test_sscrms_overlap_triggers_resampling\n"
              << "  Source file:   sweep_separation_complex_rot_memtest_sphere.cpp\n"
              << "  Function:      sweep_separation_complex_rot_memtest_sphere\n"
              << "  Scenario:      two single-molecule complexes whose interfaces\n"
              << "                 would end 0.1 nm apart, binding radius 1.0 nm.\n"
              << "  Pass criteria: (a) the partner complex trajectory is re-sampled\n"
              << "                     (becomes non-zero), (b) the target complex is\n"
              << "                     propagated (its molecule moves) and its stored\n"
              << "                     trajectory vectors end up zeroed, (c) the\n"
              << "                     partner molecule itself is NOT moved here.\n";

    sscrms_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrms_make_sphere_membrane(100.0);
    std::vector<MolTemplate> molTemplateList { sscrms_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrms_make_rxn(1.0) }; // binding radius 1.0 nm

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrms_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrms_make_molecule(1, 1, Coord { 0.1, 0.0, 0.0 }, Coord { 0.1, 0.0, 0.0 }));

    // Partner index 1 > pro1Index 0 and its trajStatus is 'none' -> eligible for
    // re-sampling inside the routine.
    sscrms_add_cross_partner(moleculeList[0], 1);

    std::vector<Complex> complexList;
    complexList.push_back(sscrms_make_complex(0, Coord { 0.0, 0.0, 0.0 }, std::vector<int> { 0 }, false, 5.0));
    complexList.push_back(sscrms_make_complex(1, Coord { 0.1, 0.0, 0.0 }, std::vector<int> { 1 }, false, 5.0));

    // Both trajectories start at zero -> the predicted separation is the current
    // 0.1 nm, which is smaller than the 1.0 nm binding radius -> overlap.
    const Coord mol0Start = moleculeList[0].comCoord;
    const Coord mol1Start = moleculeList[1].comCoord;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_sphere(simItr=3, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_sphere(3, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);

    const double partnerTrajMag = std::sqrt(complexList[1].trajTrans.x * complexList[1].trajTrans.x
        + complexList[1].trajTrans.y * complexList[1].trajTrans.y
        + complexList[1].trajTrans.z * complexList[1].trajTrans.z);
    const double targetDisplacement = sscrms_distance(moleculeList[0].comCoord, mol0Start);

    std::cerr << "  Partner (complex 1) re-sampled |trajTrans| = " << partnerTrajMag << '\n';
    std::cerr << "  Target molecule displacement          = " << targetDisplacement << '\n';

    // (a) the partner complex had a fresh trajectory drawn for it
    EXPECT_GT(partnerTrajMag, 0.0)
        << "overlapping partner complex should have had its trajectory re-sampled";

    // (b) the target complex was propagated with the (re-sampled) trajectory
    EXPECT_GT(targetDisplacement, 0.0)
        << "target molecule should have been displaced by the re-sampled trajectory";
    EXPECT_TRUE(std::isfinite(moleculeList[0].comCoord.x)) << "target x must remain finite";
    EXPECT_TRUE(std::isfinite(moleculeList[0].comCoord.y)) << "target y must remain finite";
    EXPECT_TRUE(std::isfinite(moleculeList[0].comCoord.z)) << "target z must remain finite";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "target trajTrans.x must be zeroed at the end";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "target trajTrans.y must be zeroed at the end";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "target trajTrans.z must be zeroed at the end";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.x, 0.0) << "target trajRot.x must be zeroed at the end";

    // (c) the partner complex is not propagated by this call: it only stores the
    //     new trajectory for when its own sweep is performed.
    EXPECT_NEAR(moleculeList[1].comCoord.x, mol1Start.x, 1e-12) << "partner x must not move here";
    EXPECT_NEAR(moleculeList[1].comCoord.y, mol1Start.y, 1e-12) << "partner y must not move here";
    EXPECT_NEAR(moleculeList[1].comCoord.z, mol1Start.z, 1e-12) << "partner z must not move here";

    // The target must still be inside the spherical boundary.
    const double radial = moleculeList[0].comCoord.get_magnitude();
    std::cerr << "  Target radial distance after sweep    = " << radial
              << " (sphere radius " << membraneObject.sphereR << ")\n";
    EXPECT_LE(radial, membraneObject.sphereR + 1e-6)
        << "propagated molecule should remain inside the spherical boundary";
}

// -----------------------------------------------------------------------------
// Test 5: the partner complex is flagged OnSurface, which makes the routine take
//         the spherical-surface code path (calculate_update_position_interface)
//         when predicting the partner interface position.  This test mainly
//         checks that the surface branch is reachable and leaves the system in a
//         sane state; the exact spherical geometry is validated elsewhere.
// -----------------------------------------------------------------------------
void test_sscrms_partner_on_sphere_surface()
{
    std::cerr << "\n[TEST] test_sscrms_partner_on_sphere_surface\n"
              << "  Source file:   sweep_separation_complex_rot_memtest_sphere.cpp\n"
              << "  Function:      sweep_separation_complex_rot_memtest_sphere\n"
              << "  Scenario:      target complex is in solution, the encountered\n"
              << "                 partner complex is flagged OnSurface so the\n"
              << "                 spherical-surface interface prediction is used.\n"
              << "  Pass criteria: the routine completes, the target complex is\n"
              << "                 propagated with finite coordinates and its stored\n"
              << "                 trajectory vectors end up zeroed.\n";

    sscrms_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    const double sphereR = 100.0;
    Membrane membraneObject = sscrms_make_sphere_membrane(sphereR);
    std::vector<MolTemplate> molTemplateList { sscrms_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrms_make_rxn(1.0) };

    // Target sits just inside the shell; the surface partner sits exactly on it.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrms_make_molecule(0, 0, Coord { sphereR - 0.5, 0.0, 0.0 },
        Coord { sphereR - 0.5, 0.0, 0.0 }));
    moleculeList.push_back(sscrms_make_molecule(1, 1, Coord { sphereR, 0.0, 0.0 },
        Coord { sphereR, 0.0, 0.0 }));
    moleculeList[1].isLipid = true; // membrane-bound partner

    sscrms_add_cross_partner(moleculeList[0], 1);

    std::vector<Complex> complexList;
    complexList.push_back(sscrms_make_complex(0, Coord { sphereR - 0.5, 0.0, 0.0 },
        std::vector<int> { 0 }, false, 1.0));
    // OnSurface == true -> the "complex on sphere surface" branch of the sweep.
    complexList.push_back(sscrms_make_complex(1, Coord { sphereR, 0.0, 0.0 },
        std::vector<int> { 1 }, true, 1.0));

    // Give the surface partner a small tangential displacement so the spherical
    // translation helper has a well defined direction to work with.
    complexList[1].trajTrans.x = 0.0;
    complexList[1].trajTrans.y = 0.05;
    complexList[1].trajTrans.z = 0.0;

    // Small displacement for the target so it cannot leave the sphere on its own.
    complexList[0].trajTrans.x = -0.1;
    complexList[0].trajTrans.y = 0.0;
    complexList[0].trajTrans.z = 0.0;

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_sphere(simItr=4, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_sphere(4, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membraneObject);

    std::cerr << "  Molecule 0 -> (" << moleculeList[0].comCoord.x << ", " << moleculeList[0].comCoord.y
              << ", " << moleculeList[0].comCoord.z << ")\n";

    // Regardless of whether the surface prediction registered an overlap, the
    // routine must leave the target complex propagated with a cleared trajectory
    // and finite coordinates.
    EXPECT_TRUE(std::isfinite(moleculeList[0].comCoord.x)) << "target x must remain finite";
    EXPECT_TRUE(std::isfinite(moleculeList[0].comCoord.y)) << "target y must remain finite";
    EXPECT_TRUE(std::isfinite(moleculeList[0].comCoord.z)) << "target z must remain finite";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "target trajTrans.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "target trajTrans.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "target trajTrans.z must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.z, 0.0) << "target trajRot.z must be zeroed";
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "target molecule should be flagged as propagated";

    // The surface partner is never propagated by this routine.
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::none))
        << "surface partner must not be propagated by this call";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so a failure
// in one scenario does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(SweepSepComplexRotMemtestSphere, NoCrossPartners) { test_sscrms_no_cross_partners(); }
TEST(SweepSepComplexRotMemtestSphere, SameComplexPartnerSkipped) { test_sscrms_same_complex_partner_is_skipped(); }
TEST(SweepSepComplexRotMemtestSphere, ImplicitLipidPartnerSkipped) { test_sscrms_implicit_lipid_partner_is_skipped(); }
TEST(SweepSepComplexRotMemtestSphere, OverlapTriggersResampling) { test_sscrms_overlap_triggers_resampling(); }
TEST(SweepSepComplexRotMemtestSphere, PartnerOnSphereSurface) { test_sscrms_partner_on_sphere_surface(); }