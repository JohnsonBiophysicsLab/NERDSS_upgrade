/*! \file test_sweep_separation_complex_rot_memtest_box.cpp
 *
 * ### Unit test for ../src/trajectory_functions/sweep_separation_complex_rot_memtest_box.cpp
 *
 * Function under test:
 *
 *     void sweep_separation_complex_rot_memtest_box(int simItr, int pro1Index,
 *              Parameters& params, std::vector<Molecule>& moleculeList,
 *              std::vector<Complex>& complexList,
 *              const std::vector<ForwardRxn>& forwardRxns,
 *              const std::vector<MolTemplate>& molTemplateList,
 *              const Membrane& membraneObject)
 *
 * The routine performs an overlap "sweep" for the complex that owns molecule
 * `pro1Index`:
 *   * For every molecule that was flagged as a possible reaction partner
 *     (Molecule::crossbase / mycrossint / crossrxn) it computes the separation
 *     that the two reacting interfaces *would* have after the currently stored
 *     trajectories (trajTrans + trajRot) are applied.
 *   * If the partner complex is flagged `OnSurface` (the "memtest" part) the
 *     z-component of the separation is ignored, i.e. only the xy distance is
 *     compared to the reaction binding radius.
 *   * Partners inside the *same* complex are never checked (they cannot move
 *     relative to each other), and implicit lipids are skipped entirely.
 *   * When an overlap is detected, new random trajectories are drawn (via
 *     GaussV()) for this complex and for every overlapping partner complex
 *     whose molecules may still be resampled, then reflected back into the box.
 *   * Finally the complex is propagated and its trajTrans / trajRot are zeroed.
 *
 * The tests below build tiny, fully-initialised systems and check exactly those
 * behaviours.  Because a fresh random trajectory drawn from a Gaussian is
 * essentially never identically zero, "was resampled" is verified by checking
 * that a partner complex whose trajectory started at (0,0,0) is no longer zero.
 */

#include <gtest/gtest.h>

#include <gsl/gsl_rng.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

// The GSL random number generator lives in gtest_main.cpp; we only reference it.
extern gsl_rng* r;

// -----------------------------------------------------------------------------
// Local helpers.  All names carry the `sscrmb_` prefix (Sweep Separation Complex
// Rot Memtest Box) so that they cannot collide with other tests in the suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief (Re)initialise the global GSL random number generator deterministically.
 *
 * The sweep routine calls GaussV() whenever it has to resample a trajectory, so
 * the generator must be alive before the function is invoked.
 */
void sscrmb_seed_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief Build a minimal MolTemplate with a single interface at the COM.
 *
 * Rotational diffusion is deliberately left at zero so that every resampled
 * rotation is exactly zero.  This keeps Complex::propagate() on its
 * translation-only fast path, which makes the position assertions exact.
 */
MolTemplate sscrmb_make_template()
{
    MolTemplate temp;
    temp.molTypeIndex = 0;
    temp.molName = "sweeper";
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.copies = 0;
    temp.D = Coord { 1.0, 1.0, 1.0 };   // non-zero translational diffusion
    temp.Dr = Coord { 0.0, 0.0, 0.0 };  // zero rotational diffusion -> no rotation
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPromoter = false;
    temp.isRod = false;
    temp.isPoint = true;
    temp.insideCompartment = false;
    temp.outsideCompartment = false;

    // One interface, coincident with the center of mass.
    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord { 0.0, 0.0, 0.0 };
    iface.stateList.emplace_back(std::string("a"), 0);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a single-interface molecule sitting at `com`. */
Molecule sscrmb_make_molecule(int index, int comIndex, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = 0;
    mol.myComIndex = comIndex;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;

    // Interface 0 lies exactly on the COM, so the interface-to-COM vector is
    // the zero vector and the computed separations are just COM separations.
    Molecule::Iface iface;
    iface.coord = com;
    iface.relIndex = 0;
    iface.index = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a complex holding the listed member molecules.
 *
 * The center of mass is the (equal-mass) average of its members, matching what
 * Complex::update_properties() would compute.
 */
Complex sscrmb_make_complex(int index, const std::vector<int>& members,
    const std::vector<Molecule>& moleculeList)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.memberList = members;

    Coord center { 0.0, 0.0, 0.0 };
    for (int m : members)
        center += moleculeList[m].comCoord;
    double n = static_cast<double>(members.size());
    center.x /= n;
    center.y /= n;
    center.z /= n;
    com.comCoord = center;

    com.mass = n;
    com.radius = 1.0;
    com.D = Coord { 1.0, 1.0, 1.0 };
    com.Dr = Coord { 0.0, 0.0, 0.0 };
    com.isEmpty = false;
    com.OnSurface = false;
    com.onFiber = false;
    com.linksToSurface = 0;
    com.trajStatus = TrajStatus::none;
    com.trajTrans = Vector(0.0, 0.0, 0.0);
    com.trajRot = Coord { 0.0, 0.0, 0.0 };
    com.numEachMol = std::vector<int>(1, static_cast<int>(members.size()));
    com.lastNumberUpdateItrEachMol = std::vector<long long int>(1, 0);

    return com;
}

/*! \brief Register `partnerIndex` as a possible reaction partner of `mol`. */
void sscrmb_add_cross(Molecule& mol, int partnerIndex, int myIfaceRelIndex, int rxnIndex)
{
    mol.crossbase.push_back(partnerIndex);
    mol.mycrossint.push_back(myIfaceRelIndex);
    mol.crossrxn.push_back(std::array<int, 3> { rxnIndex, 0, 0 });
    mol.probvec.push_back(1.0);
}

/*! \brief Build a bimolecular ForwardRxn whose two reactants both use iface 0. */
ForwardRxn sscrmb_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    // (name, molTypeIndex, absIfaceIndex, relIfaceIndex, requiresState, requiresInteraction)
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.productListNew.emplace_back("a", 0, 0, 0, '\0', true);
    rxn.productListNew.emplace_back("a", 0, 0, 0, '\0', true);
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = 1.0;
    return rxn;
}

/*! \brief Build a simple reflecting cubic water box, no sphere, no compartment. */
Membrane sscrmb_make_membrane()
{
    Membrane membrane;
    membrane.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
    membrane.isSphere = false;
    membrane.isBox = true;
    membrane.implicitLipid = false;   // avoids the RS3Dvect lookup table entirely
    membrane.hasCompartment = false;
    membrane.xBCtype = "reflect";
    membrane.yBCtype = "reflect";
    membrane.zBCtype = "reflect";
    membrane.TwoD = false;
    return membrane;
}

/*! \brief Convenience predicate: is this translation vector exactly zero? */
bool sscrmb_traj_is_zero(const Complex& com)
{
    return com.trajTrans.x == 0.0 && com.trajTrans.y == 0.0 && com.trajTrans.z == 0.0;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: no partner at all -> the complex is simply propagated by the
//         trajectory that was already stored, and its trajectory is zeroed.
// -----------------------------------------------------------------------------
void test_sscrmb_no_partners_just_propagates()
{
    std::cerr << "\n[TEST] test_sscrmb_no_partners_just_propagates\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    a lone complex with an empty crossbase list.\n"
              << "  Pass criteria: molecule moves by exactly trajTrans, complex\n"
              << "                 trajTrans/trajRot are zeroed, trajStatus is\n"
              << "                 set to propagated.\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 10.0, 10.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0 }, moleculeList));

    // Pre-set displacement; there is nothing to collide with so it must be kept.
    complexList[0].trajTrans = Vector(2.0, 3.0, -1.0);

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_box(pro1Index = 0)...\n";
    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Molecule 0 COM after sweep = (" << moleculeList[0].comCoord.x << ", "
              << moleculeList[0].comCoord.y << ", " << moleculeList[0].comCoord.z << ")\n";

    // Translation-only fast path (Dr == 0) -> exact arithmetic.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 12.0) << "x should be 10 + 2";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, 13.0) << "y should be 10 + 3";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, -1.0) << "z should be 0 - 1";

    // The complex COM is recomputed inside propagate()/update_properties().
    EXPECT_DOUBLE_EQ(complexList[0].comCoord.x, 12.0) << "complex COM should track its lone member";
    EXPECT_DOUBLE_EQ(complexList[0].comCoord.y, 13.0) << "complex COM should track its lone member";

    // Trajectories must be cleared at the end of the routine.
    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[0])) << "trajTrans must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.x, 0.0) << "trajRot.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.y, 0.0) << "trajRot.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.z, 0.0) << "trajRot.z must be zeroed";

    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "member molecule should be flagged as propagated";
}

// -----------------------------------------------------------------------------
// Test 2: partner exists but is far away -> no overlap, so neither complex has
//         its trajectory resampled and the stored move is applied verbatim.
// -----------------------------------------------------------------------------
void test_sscrmb_far_partner_no_resample()
{
    std::cerr << "\n[TEST] test_sscrmb_far_partner_no_resample\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    partner complex is 50 nm away, binding radius 5 nm.\n"
              << "  Pass criteria: partner trajTrans stays exactly (0,0,0) and the\n"
              << "                 swept complex moves by its original trajTrans.\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrmb_make_molecule(1, 1, Coord { 50.0, 0.0, 0.0 }));

    // Molecule 0 lists molecule 1 as a potential reaction partner (and back).
    sscrmb_add_cross(moleculeList[0], 1, 0, 0);
    sscrmb_add_cross(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0 }, moleculeList));
    complexList.push_back(sscrmb_make_complex(1, std::vector<int> { 1 }, moleculeList));

    complexList[0].trajTrans = Vector(1.0, 0.0, 0.0);

    std::cerr << "  Separation before sweep = 50 nm, bindRadius = 5 nm -> no overlap expected.\n";
    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Molecule 0 x after sweep = " << moleculeList[0].comCoord.x
              << " (expected 1)\n";

    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 1.0)
        << "no overlap -> the original translation must be applied unchanged";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, 0.0) << "y should be untouched";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 0.0) << "z should be untouched";

    // Partner complex must not have been touched at all.
    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[1]))
        << "partner complex trajectory must not be resampled when there is no overlap";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 50.0)
        << "partner molecule is never propagated by this routine";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::none))
        << "partner molecule trajStatus must stay untouched";

    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[0])) << "swept complex trajTrans is zeroed at the end";
}

// -----------------------------------------------------------------------------
// Test 3: overlapping partner in a *different* complex -> both trajectories are
//         resampled (partner trajectory becomes non-zero).
// -----------------------------------------------------------------------------
void test_sscrmb_overlap_resamples_both()
{
    std::cerr << "\n[TEST] test_sscrmb_overlap_resamples_both\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    two separate complexes 1 nm apart, binding radius 5 nm.\n"
              << "  Pass criteria: the overlapping partner complex gets a freshly\n"
              << "                 sampled (non-zero) trajTrans, the swept complex is\n"
              << "                 moved away from its original position, and the\n"
              << "                 partner molecule itself is not propagated.\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrmb_make_molecule(1, 1, Coord { 1.0, 0.0, 0.0 }));

    sscrmb_add_cross(moleculeList[0], 1, 0, 0);
    sscrmb_add_cross(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0 }, moleculeList));
    complexList.push_back(sscrmb_make_complex(1, std::vector<int> { 1 }, moleculeList));

    const Coord partnerStart = moleculeList[1].comCoord;

    std::cerr << "  Separation before sweep = 1 nm < bindRadius 5 nm -> overlap expected.\n";
    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Partner complex trajTrans after sweep = ("
              << complexList[1].trajTrans.x << ", " << complexList[1].trajTrans.y << ", "
              << complexList[1].trajTrans.z << ")\n";

    // The partner complex started with a zero trajectory; the routine must have
    // drawn a new Gaussian displacement for it.
    EXPECT_FALSE(sscrmb_traj_is_zero(complexList[1]))
        << "overlapping partner complex should have been given a resampled trajectory";

    // The partner is only *scheduled* to move: this routine propagates the swept
    // complex only.
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, partnerStart.x)
        << "partner molecule coordinates are not propagated here";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.y, partnerStart.y)
        << "partner molecule coordinates are not propagated here";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.z, partnerStart.z)
        << "partner molecule coordinates are not propagated here";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::none))
        << "partner molecule keeps trajStatus none so it can still be resampled later";

    // The swept complex was moved by a randomly resampled (non-zero) displacement.
    const double moved = std::sqrt(moleculeList[0].comCoord.x * moleculeList[0].comCoord.x
        + moleculeList[0].comCoord.y * moleculeList[0].comCoord.y
        + moleculeList[0].comCoord.z * moleculeList[0].comCoord.z);
    std::cerr << "  Swept molecule displacement magnitude = " << moved << '\n';
    EXPECT_GT(moved, 0.0) << "the swept complex should have been displaced by the resampled move";

    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[0])) << "swept complex trajTrans is zeroed at the end";
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "swept molecule should be flagged as propagated";
}

// -----------------------------------------------------------------------------
// Test 4: overlapping partner inside the *same* complex -> the check is skipped
//         entirely because members of a complex cannot diffuse relative to each
//         other, so the stored trajectory is applied verbatim.
// -----------------------------------------------------------------------------
void test_sscrmb_same_complex_partner_is_skipped()
{
    std::cerr << "\n[TEST] test_sscrmb_same_complex_partner_is_skipped\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    both 'partners' belong to the same complex and sit\n"
              << "               1 nm apart with a 5 nm binding radius.\n"
              << "  Pass criteria: no resampling occurs; both members translate by\n"
              << "                 exactly the pre-set trajTrans.\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrmb_make_molecule(1, 0, Coord { 1.0, 0.0, 0.0 })); // same complex 0

    sscrmb_add_cross(moleculeList[0], 1, 0, 0);
    sscrmb_add_cross(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0, 1 }, moleculeList));

    complexList[0].trajTrans = Vector(4.0, 0.0, 0.0);

    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Member positions after sweep: x0 = " << moleculeList[0].comCoord.x
              << ", x1 = " << moleculeList[1].comCoord.x << " (expected 4 and 5)\n";

    // Because the overlap test is bypassed the original 4 nm translation survives.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 4.0)
        << "same-complex partners must not trigger a resample";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 5.0)
        << "the second member is translated by the same vector";
    EXPECT_DOUBLE_EQ(complexList[0].comCoord.x, 4.5)
        << "complex COM is the mass-weighted average of its two members";

    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[0])) << "trajTrans is zeroed after propagation";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "every member of the swept complex is flagged as propagated";
}

// -----------------------------------------------------------------------------
// Test 5: an implicit-lipid partner is skipped by the `continue` statement, so
//         even a perfectly coincident partner produces no resampling.
// -----------------------------------------------------------------------------
void test_sscrmb_implicit_lipid_partner_skipped()
{
    std::cerr << "\n[TEST] test_sscrmb_implicit_lipid_partner_skipped\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    the only cross-partner is flagged isImplicitLipid.\n"
              << "  Pass criteria: the partner is skipped, so the stored trajectory\n"
              << "                 is applied unchanged and the partner complex keeps\n"
              << "                 its zero trajectory.\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }));
    // Coincident partner that would definitely overlap - but it is implicit.
    moleculeList.push_back(sscrmb_make_molecule(1, 1, Coord { 0.0, 0.0, 0.0 }));
    moleculeList[1].isImplicitLipid = true;

    sscrmb_add_cross(moleculeList[0], 1, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0 }, moleculeList));
    complexList.push_back(sscrmb_make_complex(1, std::vector<int> { 1 }, moleculeList));

    complexList[0].trajTrans = Vector(0.0, 2.0, 0.0);

    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Molecule 0 y after sweep = " << moleculeList[0].comCoord.y
              << " (expected 2 - i.e. the move was NOT resampled)\n";

    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, 2.0)
        << "implicit lipid partners are skipped, so the original move survives";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 0.0) << "x should be untouched";
    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[1]))
        << "the implicit lipid complex must never be resampled";
}

// -----------------------------------------------------------------------------
// Test 6a: "memtest" behaviour - when the partner complex is OnSurface the
//          z separation is ignored, so a large z gap still counts as an overlap.
// -----------------------------------------------------------------------------
void test_sscrmb_membrane_partner_ignores_z()
{
    std::cerr << "\n[TEST] test_sscrmb_membrane_partner_ignores_z\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    partner complex has OnSurface = true and is 20 nm\n"
              << "               away in z but only 1 nm away in xy (bindRadius 5 nm).\n"
              << "  Pass criteria: the z component is ignored -> overlap detected ->\n"
              << "                 the partner complex trajectory is resampled.\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrmb_make_molecule(1, 1, Coord { 0.0, 0.0, 20.0 }));

    sscrmb_add_cross(moleculeList[0], 1, 0, 0);
    sscrmb_add_cross(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0 }, moleculeList));
    complexList.push_back(sscrmb_make_complex(1, std::vector<int> { 1 }, moleculeList));

    // Both complexes are treated as membrane-bound: the partner's OnSurface flag
    // is what switches on the xy-only distance test.
    complexList[0].OnSurface = true;
    complexList[1].OnSurface = true;

    complexList[0].trajTrans = Vector(1.0, 0.0, 0.0); // xy separation would be 1 nm

    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Partner trajTrans after sweep = (" << complexList[1].trajTrans.x << ", "
              << complexList[1].trajTrans.y << ", " << complexList[1].trajTrans.z << ")\n";

    EXPECT_FALSE(sscrmb_traj_is_zero(complexList[1]))
        << "with OnSurface partners only xy is compared, so this must be an overlap";
    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[0]))
        << "the swept complex trajectory is still zeroed at the end";
}

// -----------------------------------------------------------------------------
// Test 6b: the mirror image of 6a - if the partner is NOT OnSurface the z gap is
//          included, the pair is far apart, and nothing is resampled.
// -----------------------------------------------------------------------------
void test_sscrmb_solution_partner_counts_z()
{
    std::cerr << "\n[TEST] test_sscrmb_solution_partner_counts_z\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    identical geometry to the previous test but the\n"
              << "               partner complex has OnSurface = false.\n"
              << "  Pass criteria: z is included in the separation (20 nm) -> no\n"
              << "                 overlap -> no resampling, exact translation.\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrmb_make_molecule(1, 1, Coord { 0.0, 0.0, 20.0 }));

    sscrmb_add_cross(moleculeList[0], 1, 0, 0);
    sscrmb_add_cross(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0 }, moleculeList));
    complexList.push_back(sscrmb_make_complex(1, std::vector<int> { 1 }, moleculeList));

    complexList[0].OnSurface = true;
    complexList[1].OnSurface = false; // <-- solution partner: z matters

    complexList[0].trajTrans = Vector(1.0, 0.0, 0.0);

    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Molecule 0 x after sweep = " << moleculeList[0].comCoord.x
              << " (expected 1 - move kept), partner trajTrans zero? "
              << std::boolalpha << sscrmb_traj_is_zero(complexList[1]) << '\n';

    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[1]))
        << "including z gives a 20 nm separation, far beyond the 5 nm binding radius";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 1.0)
        << "no overlap -> the pre-set translation is applied unchanged";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 0.0) << "z should be untouched";
}

// -----------------------------------------------------------------------------
// Test 7: several partners, one overlapping and one far away -> only the
//         overlapping partner complex is resampled.
// -----------------------------------------------------------------------------
void test_sscrmb_only_overlapping_partner_resampled()
{
    std::cerr << "\n[TEST] test_sscrmb_only_overlapping_partner_resampled\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    two cross-partners: one 1 nm away (overlapping) and\n"
              << "               one 60 nm away (not overlapping), bindRadius 5 nm.\n"
              << "  Pass criteria: only the close partner's complex trajectory is\n"
              << "                 resampled; the distant one keeps trajTrans (0,0,0).\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrmb_make_molecule(1, 1, Coord { 1.0, 0.0, 0.0 }));   // close
    moleculeList.push_back(sscrmb_make_molecule(2, 2, Coord { 60.0, 0.0, 0.0 }));  // far

    sscrmb_add_cross(moleculeList[0], 1, 0, 0);
    sscrmb_add_cross(moleculeList[0], 2, 0, 0);
    sscrmb_add_cross(moleculeList[1], 0, 0, 0);
    sscrmb_add_cross(moleculeList[2], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0 }, moleculeList));
    complexList.push_back(sscrmb_make_complex(1, std::vector<int> { 1 }, moleculeList));
    complexList.push_back(sscrmb_make_complex(2, std::vector<int> { 2 }, moleculeList));

    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Close partner trajTrans zero? " << std::boolalpha
              << sscrmb_traj_is_zero(complexList[1])
              << ", far partner trajTrans zero? " << sscrmb_traj_is_zero(complexList[2]) << '\n';

    EXPECT_FALSE(sscrmb_traj_is_zero(complexList[1]))
        << "the overlapping partner complex must be resampled";
    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[2]))
        << "the distant partner complex must be left alone";

    // The swept complex was resampled and then propagated; its trajectory is cleared.
    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[0]))
        << "swept complex trajTrans is zeroed after propagation";
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "swept molecule should be flagged as propagated";
}

// -----------------------------------------------------------------------------
// Test 8: a partner whose molecules have already been propagated must NOT be
//         resampled (only the swept complex is moved).
// -----------------------------------------------------------------------------
void test_sscrmb_propagated_partner_not_resampled()
{
    std::cerr << "\n[TEST] test_sscrmb_propagated_partner_not_resampled\n"
              << "  Source file: sweep_separation_complex_rot_memtest_box.cpp\n"
              << "  Function:    sweep_separation_complex_rot_memtest_box\n"
              << "  Scenario:    overlapping partner whose molecule already has\n"
              << "               trajStatus == propagated.\n"
              << "  Pass criteria: the partner's trajectory is left at (0,0,0) while\n"
              << "                 the swept complex is still resampled and moved.\n";

    sscrmb_seed_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membrane = sscrmb_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrmb_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrmb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmb_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrmb_make_molecule(1, 1, Coord { 1.0, 0.0, 0.0 }));
    // This partner has already moved this timestep and may not be resampled.
    moleculeList[1].trajStatus = TrajStatus::propagated;

    sscrmb_add_cross(moleculeList[0], 1, 0, 0);
    sscrmb_add_cross(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmb_make_complex(0, std::vector<int> { 0 }, moleculeList));
    complexList.push_back(sscrmb_make_complex(1, std::vector<int> { 1 }, moleculeList));

    sweep_separation_complex_rot_memtest_box(0, 0, params, moleculeList, complexList,
        forwardRxns, molTemplateList, membrane);

    std::cerr << "  Partner (already propagated) trajTrans = (" << complexList[1].trajTrans.x
              << ", " << complexList[1].trajTrans.y << ", " << complexList[1].trajTrans.z << ")\n";

    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[1]))
        << "an already-propagated partner must keep its trajectory";

    // The swept complex still had its own trajectory resampled, so it moved.
    const double moved = std::sqrt(moleculeList[0].comCoord.x * moleculeList[0].comCoord.x
        + moleculeList[0].comCoord.y * moleculeList[0].comCoord.y
        + moleculeList[0].comCoord.z * moleculeList[0].comCoord.z);
    std::cerr << "  Swept molecule displacement magnitude = " << moved << '\n';
    EXPECT_GT(moved, 0.0) << "the swept complex should still have been resampled and moved";

    EXPECT_TRUE(sscrmb_traj_is_zero(complexList[0]))
        << "swept complex trajTrans is zeroed after propagation";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is executed inside its own TEST so a
// failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(SweepSepComplexRotMemtestBox, NoPartnersJustPropagates) { test_sscrmb_no_partners_just_propagates(); }
TEST(SweepSepComplexRotMemtestBox, FarPartnerNoResample) { test_sscrmb_far_partner_no_resample(); }
TEST(SweepSepComplexRotMemtestBox, OverlapResamplesBoth) { test_sscrmb_overlap_resamples_both(); }
TEST(SweepSepComplexRotMemtestBox, SameComplexPartnerIsSkipped) { test_sscrmb_same_complex_partner_is_skipped(); }
TEST(SweepSepComplexRotMemtestBox, ImplicitLipidPartnerSkipped) { test_sscrmb_implicit_lipid_partner_skipped(); }
TEST(SweepSepComplexRotMemtestBox, MembranePartnerIgnoresZ) { test_sscrmb_membrane_partner_ignores_z(); }
TEST(SweepSepComplexRotMemtestBox, SolutionPartnerCountsZ) { test_sscrmb_solution_partner_counts_z(); }
TEST(SweepSepComplexRotMemtestBox, OnlyOverlappingPartnerResampled) { test_sscrmb_only_overlapping_partner_resampled(); }
TEST(SweepSepComplexRotMemtestBox, PropagatedPartnerNotResampled) { test_sscrmb_propagated_partner_not_resampled(); }