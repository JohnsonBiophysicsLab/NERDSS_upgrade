/*! \file test_sweep_separation_complex_rot_box.cpp
 *
 * ### Unit test for ../src/trajectory_functions/sweep_separation_complex_rot_box.cpp
 *
 * Function under test:
 *
 *     void sweep_separation_complex_rot_box(int simItr, int pro1Index,
 *                                           Parameters& params,
 *                                           std::vector<Molecule>& moleculeList,
 *                                           std::vector<Complex>& complexList,
 *                                           const std::vector<ForwardRxn>& forwardRxns,
 *                                           const std::vector<MolTemplate>& molTemplateList,
 *                                           const Membrane& membraneObject)
 *
 * The routine walks over every molecule of the complex that owns `pro1Index`,
 * looks at each of that molecule's "crossing" partners (crossbase / mycrossint /
 * crossrxn), and computes the projected interface-interface separation *after*
 * the currently-stored trajectory (trajTrans + trajRot) would be applied.
 *
 *   * If any partner pair ends up closer than the reaction binding radius the
 *     move is rejected: a brand new Gaussian trajectory is drawn for the target
 *     complex and (if they have not moved yet this step) for each overlapping
 *     partner complex. Up to 10 attempts are made.
 *   * Partners that live in the *same* complex are skipped, because they cannot
 *     diffuse relative to one another.
 *   * Implicit-lipid partners are skipped entirely.
 *   * Finally the target complex is propagated and its trajTrans/trajRot are
 *     zeroed.
 *
 * The tests below build tiny, fully-initialised systems (one interface per
 * molecule) so that every one of those code paths can be exercised, and so the
 * "no overlap" paths are completely deterministic (no random numbers are drawn
 * at all when nothing overlaps).
 *
 * NOTE: `membraneObject.implicitLipid` is deliberately left false in every test.
 *       When it is true the routine indexes `membraneObject.RS3Dvect` at offsets
 *       300..499, which would read out of bounds for a default-constructed
 *       Membrane.
 */

#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// The GSL random number generator is defined once in gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Small helpers used by all of the tests in this file.  Names are prefixed with
// "sscrb_" (sweep_separation_complex_rot_box) so they cannot collide with the
// helpers of any other generated test translation unit.
// -----------------------------------------------------------------------------

/*! \brief Deterministically (re)seed the GSL RNG used by GaussV().
 *
 * The instructions for this test suite require this exact initialisation
 * sequence instead of calling srand_gsl().
 */
void sscrb_init_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief Build a MolTemplate with exactly one interface.
 *
 * `Complex::update_properties()` (called from `Complex::propagate()`) reads the
 * template radius and the translational/rotational diffusion constants, so all
 * of those must be filled in.
 */
MolTemplate sscrb_make_template(int typeIndex, const std::string& name, double transD)
{
    MolTemplate temp;
    temp.molTypeIndex = typeIndex;
    temp.molName = name;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord { transD, transD, transD };
    temp.Dr = Coord { 0.0, 0.0, 0.0 }; // no rotational diffusion -> no rotation resampling
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPoint = false;
    temp.isRod = false;
    temp.isPromoter = false;
    temp.insideCompartment = false;
    temp.outsideCompartment = false;
    temp.crossesCompartment = false;

    // One interface, one (default) state.
    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord { 1.0, 0.0, 0.0 };
    iface.stateList.emplace_back(std::string { "a" }, 0);
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a fully-initialised Molecule with a single interface. */
Molecule sscrb_make_molecule(int index, int comIndex, int molTypeIndex, const Coord& com, const Coord& ifaceCoord)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;

    Molecule::Iface iface;
    iface.coord = ifaceCoord;
    iface.index = 0;
    iface.relIndex = 0;
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a Complex whose center of mass is the mass-weighted average of
 *         its member molecules (which is exactly what update_properties() will
 *         recompute later).
 */
Complex sscrb_make_complex(int index, const std::vector<int>& members, const std::vector<Molecule>& moleculeList,
    const Coord& D, const Coord& Dr, size_t numMolTypes)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.memberList = members;
    com.D = D;
    com.Dr = Dr;
    com.radius = 1.0;
    com.isEmpty = false;
    com.OnSurface = false;
    com.onFiber = false;
    com.linksToSurface = 0;
    com.iLipidIndex = 0;
    com.ncross = 0;
    com.trajStatus = TrajStatus::none;

    double totMass { 0.0 };
    Coord c { 0.0, 0.0, 0.0 };
    for (int m : members) {
        totMass += moleculeList[m].mass;
        c.x += moleculeList[m].comCoord.x * moleculeList[m].mass;
        c.y += moleculeList[m].comCoord.y * moleculeList[m].mass;
        c.z += moleculeList[m].comCoord.z * moleculeList[m].mass;
    }
    c.x /= totMass;
    c.y /= totMass;
    c.z /= totMass;
    com.mass = totMass;
    com.comCoord = c;

    com.numEachMol = std::vector<int>(numMolTypes, 0);
    com.numEachMol[0] = static_cast<int>(members.size());
    com.lastNumberUpdateItrEachMol = std::vector<long long int>(numMolTypes, 0);

    // Start with a zero trajectory; individual tests overwrite this.
    com.trajTrans.x = 0.0;
    com.trajTrans.y = 0.0;
    com.trajTrans.z = 0.0;
    com.trajRot.x = 0.0;
    com.trajRot.y = 0.0;
    com.trajRot.z = 0.0;

    return com;
}

/*! \brief Build a symmetric bimolecular ForwardRxn between interface 0 and
 *         interface 0 with the requested binding radius.
 *
 * `sweep_separation_complex_rot_box()` only reads `bindRadius` and the two
 * `reactantListNew[*].relIfaceIndex` entries from the reaction.
 */
ForwardRxn sscrb_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.reactantListNew.emplace_back(std::string { "a" }, 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back(std::string { "a" }, 0, 0, 0, '\0', false);
    rxn.productListNew.emplace_back(std::string { "a" }, 0, 0, 0, '\0', true);
    rxn.productListNew.emplace_back(std::string { "a" }, 0, 0, 0, '\0', true);
    return rxn;
}

/*! \brief Build a rectangular (non-spherical, no compartment) membrane. */
Membrane sscrb_make_box(double side)
{
    Membrane mem;
    mem.isBox = true;
    mem.isSphere = false;
    mem.hasCompartment = false;
    mem.implicitLipid = false; // avoids the RS3Dvect lookup, see file header
    mem.waterBox = Membrane::WaterBox(std::vector<double> { side, side, side });
    return mem;
}

/*! \brief Register molecule `p2` as a crossing partner of molecule `p1`. */
void sscrb_add_crossing(Molecule& p1, int p2Index, int myIfaceRelIndex, int rxnIndex)
{
    p1.crossbase.push_back(p2Index);
    p1.mycrossint.push_back(myIfaceRelIndex);
    p1.crossrxn.push_back(std::array<int, 3> { rxnIndex, 0, 0 });
}

/*! \brief Euclidean distance helper for the assertions below. */
double sscrb_dist(const Coord& a, const Coord& b)
{
    const double dx { a.x - b.x };
    const double dy { a.y - b.y };
    const double dz { a.z - b.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: two well-separated complexes -> no overlap is detected, so the stored
//         trajectory must be applied verbatim (this path draws no random
//         numbers at all, so the result is bit-for-bit predictable).
// -----------------------------------------------------------------------------
void test_sscrb_no_overlap_applies_trajectory()
{
    std::cerr << "\n[TEST] test_sscrb_no_overlap_applies_trajectory\n"
              << "  Source file:   sweep_separation_complex_rot_box.cpp\n"
              << "  Function:      sweep_separation_complex_rot_box()\n"
              << "  Scenario:      complex 0 and complex 1 are 50 nm apart while the\n"
              << "                 binding radius is only 1 nm, so nothing overlaps.\n"
              << "  Pass criteria: complex 0 (and its molecule) is translated by exactly\n"
              << "                 the pre-set trajTrans, trajTrans/trajRot are zeroed,\n"
              << "                 the molecule is flagged propagated, and complex 1 is\n"
              << "                 left completely untouched.\n";

    Parameters params;
    params.timeStep = 0.1;

    Membrane membraneObject = sscrb_make_box(500.0);

    std::vector<MolTemplate> molTemplateList { sscrb_make_template(0, "A", 1.0) };
    std::vector<ForwardRxn> forwardRxns { sscrb_make_rxn(1.0) };

    // Molecule 0 sits at the origin, molecule 1 sits 50 nm away in x.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrb_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrb_make_molecule(1, 1, 0, Coord { 50.0, 0.0, 0.0 }, Coord { 49.0, 0.0, 0.0 }));

    // Each molecule knows about the other as a possible reaction partner.
    sscrb_add_crossing(moleculeList[0], 1, 0, 0);
    sscrb_add_crossing(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrb_make_complex(0, { 0 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));
    complexList.push_back(
        sscrb_make_complex(1, { 1 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));

    // The move we want to see applied unchanged.
    complexList[0].trajTrans.x = 2.0;
    complexList[0].trajTrans.y = -1.0;
    complexList[0].trajTrans.z = 0.5;

    const Coord com1Before = complexList[1].comCoord;

    std::cerr << "  Calling sweep_separation_complex_rot_box(simItr=0, pro1Index=0)...\n";
    sweep_separation_complex_rot_box(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // The stored trajectory must have been applied exactly (no resampling).
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 2.0) << "molecule 0 should translate by trajTrans.x = 2.0";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, -1.0) << "molecule 0 should translate by trajTrans.y = -1.0";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.z, 0.5) << "molecule 0 should translate by trajTrans.z = 0.5";

    // Its interface rides along rigidly with the molecule.
    EXPECT_DOUBLE_EQ(moleculeList[0].interfaceList[0].coord.x, 3.0) << "interface should translate with the molecule";

    // The complex COM is recomputed from its (single) member.
    EXPECT_DOUBLE_EQ(complexList[0].comCoord.x, 2.0) << "complex 0 COM should follow its only member";

    // The trajectory buffers are always cleared on exit.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "trajTrans.y must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "trajTrans.z must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.x, 0.0) << "trajRot.x must be zeroed after propagation";

    // Bookkeeping performed by Complex::propagate().
    EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::propagated) << "molecule 0 should be flagged as propagated";
    EXPECT_EQ(complexList[0].trajStatus, TrajStatus::propagated) << "complex 0 should be flagged as propagated";

    // The partner complex was never resampled and never propagated.
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.x, 0.0) << "complex 1 must not be resampled when nothing overlaps";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.y, 0.0) << "complex 1 must not be resampled when nothing overlaps";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.z, 0.0) << "complex 1 must not be resampled when nothing overlaps";
    EXPECT_DOUBLE_EQ(complexList[1].comCoord.x, com1Before.x) << "only the target complex may be propagated";
    EXPECT_EQ(moleculeList[1].trajStatus, TrajStatus::none) << "partner molecule must keep trajStatus none";

    std::cerr << "  molecule 0 final COM = (" << moleculeList[0].comCoord.x << ", " << moleculeList[0].comCoord.y << ", "
              << moleculeList[0].comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: overlapping interfaces in two *different* complexes -> both complexes
//         must have their trajectories resampled.
// -----------------------------------------------------------------------------
void test_sscrb_overlap_resamples_both_complexes()
{
    std::cerr << "\n[TEST] test_sscrb_overlap_resamples_both_complexes\n"
              << "  Source file:   sweep_separation_complex_rot_box.cpp\n"
              << "  Function:      sweep_separation_complex_rot_box()\n"
              << "  Scenario:      the two reacting interfaces coincide and the binding\n"
              << "                 radius is 5 nm, so overlap is guaranteed.\n"
              << "  Pass criteria: the partner complex (which has not moved yet this\n"
              << "                 step) receives a freshly sampled non-zero trajectory,\n"
              << "                 but is NOT itself propagated; the target complex is\n"
              << "                 moved and its trajectory buffers are cleared.\n";

    sscrb_init_rng();

    Parameters params;
    params.timeStep = 0.1;

    Membrane membraneObject = sscrb_make_box(500.0);

    std::vector<MolTemplate> molTemplateList { sscrb_make_template(0, "A", 1.0) };
    std::vector<ForwardRxn> forwardRxns { sscrb_make_rxn(5.0) }; // big binding radius -> overlap

    // Both interfaces are placed at exactly the same point (separation == 0).
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrb_make_molecule(0, 0, 0, Coord { -1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrb_make_molecule(1, 1, 0, Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));

    sscrb_add_crossing(moleculeList[0], 1, 0, 0);
    sscrb_add_crossing(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrb_make_complex(0, { 0 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));
    complexList.push_back(
        sscrb_make_complex(1, { 1 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));

    const Coord mol0Before = moleculeList[0].comCoord;
    const Coord com1Before = complexList[1].comCoord;
    const Coord mol1Before = moleculeList[1].comCoord;

    std::cerr << "  Calling sweep_separation_complex_rot_box(simItr=0, pro1Index=0)...\n";
    sweep_separation_complex_rot_box(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // The partner complex must have been given a brand new (non-zero) trajectory.
    const double partnerStep = std::sqrt(complexList[1].trajTrans.x * complexList[1].trajTrans.x
        + complexList[1].trajTrans.y * complexList[1].trajTrans.y
        + complexList[1].trajTrans.z * complexList[1].trajTrans.z);
    std::cerr << "  Resampled |trajTrans| of the partner complex = " << partnerStep << '\n';
    EXPECT_GT(partnerStep, 0.0) << "the overlapping partner complex should have been resampled";

    // ...but it must not have been moved: only the target complex is propagated.
    EXPECT_DOUBLE_EQ(complexList[1].comCoord.x, com1Before.x) << "partner complex must not be propagated here";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, mol1Before.x) << "partner molecule must not be moved here";
    EXPECT_EQ(moleculeList[1].trajStatus, TrajStatus::none) << "partner molecule must still be un-propagated";

    // The target complex was resampled (D > 0 so the draw is non-zero) and moved.
    EXPECT_NE(moleculeList[0].comCoord.x, mol0Before.x) << "target molecule should have been displaced";
    EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::propagated) << "target molecule should be flagged propagated";

    // Trajectory buffers of the target are always cleared before returning.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "target trajTrans.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "target trajTrans.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "target trajTrans.z must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.z, 0.0) << "target trajRot.z must be zeroed";
}

// -----------------------------------------------------------------------------
// Test 3: overlapping partners that already moved this step must NOT be
//         resampled (trajStatus gate).
// -----------------------------------------------------------------------------
void test_sscrb_propagated_partner_is_not_resampled()
{
    std::cerr << "\n[TEST] test_sscrb_propagated_partner_is_not_resampled\n"
              << "  Source file:   sweep_separation_complex_rot_box.cpp\n"
              << "  Function:      sweep_separation_complex_rot_box()\n"
              << "  Scenario:      identical overlap to the previous test, but the\n"
              << "                 partner molecule is marked TrajStatus::propagated,\n"
              << "                 i.e. it has already been moved this timestep.\n"
              << "  Pass criteria: the partner's trajTrans stays exactly zero while the\n"
              << "                 target complex is still resampled and propagated.\n";

    sscrb_init_rng();

    Parameters params;
    params.timeStep = 0.1;

    Membrane membraneObject = sscrb_make_box(500.0);

    std::vector<MolTemplate> molTemplateList { sscrb_make_template(0, "A", 1.0) };
    std::vector<ForwardRxn> forwardRxns { sscrb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrb_make_molecule(0, 0, 0, Coord { -1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrb_make_molecule(1, 1, 0, Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));

    // The partner has already been propagated this step and may not be resampled.
    moleculeList[1].trajStatus = TrajStatus::propagated;

    sscrb_add_crossing(moleculeList[0], 1, 0, 0);
    sscrb_add_crossing(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrb_make_complex(0, { 0 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));
    complexList.push_back(
        sscrb_make_complex(1, { 1 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));

    const Coord mol0Before = moleculeList[0].comCoord;

    std::cerr << "  Calling sweep_separation_complex_rot_box(simItr=0, pro1Index=0)...\n";
    sweep_separation_complex_rot_box(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // The trajStatus gate must have blocked the resampling of complex 1.
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.x, 0.0) << "already-propagated partner must keep trajTrans.x == 0";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.y, 0.0) << "already-propagated partner must keep trajTrans.y == 0";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.z, 0.0) << "already-propagated partner must keep trajTrans.z == 0";

    // The target complex still gets a new move and is propagated.
    EXPECT_NE(moleculeList[0].comCoord.x, mol0Before.x) << "target molecule should still have been displaced";
    EXPECT_EQ(complexList[0].trajStatus, TrajStatus::propagated) << "target complex should be flagged propagated";

    std::cerr << "  partner trajTrans = (" << complexList[1].trajTrans.x << ", " << complexList[1].trajTrans.y << ", "
              << complexList[1].trajTrans.z << ") (expected all zero)\n";
}

// -----------------------------------------------------------------------------
// Test 4: partners inside the *same* complex are skipped -> no resampling even
//         though their interfaces are on top of each other.
// -----------------------------------------------------------------------------
void test_sscrb_same_complex_partners_skipped()
{
    std::cerr << "\n[TEST] test_sscrb_same_complex_partners_skipped\n"
              << "  Source file:   sweep_separation_complex_rot_box.cpp\n"
              << "  Function:      sweep_separation_complex_rot_box()\n"
              << "  Scenario:      both crossing partners belong to complex 0 and their\n"
              << "                 interfaces coincide (binding radius 5 nm).\n"
              << "  Pass criteria: because molecules inside one complex cannot diffuse\n"
              << "                 relative to each other, the overlap test is skipped\n"
              << "                 and the pre-set trajectory is applied verbatim to\n"
              << "                 BOTH member molecules.\n";

    Parameters params;
    params.timeStep = 0.1;

    Membrane membraneObject = sscrb_make_box(500.0);

    std::vector<MolTemplate> molTemplateList { sscrb_make_template(0, "A", 1.0) };
    std::vector<ForwardRxn> forwardRxns { sscrb_make_rxn(5.0) };

    // Two molecules, both members of complex 0, interfaces at the same point.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrb_make_molecule(0, 0, 0, Coord { -1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrb_make_molecule(1, 0, 0, Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));

    sscrb_add_crossing(moleculeList[0], 1, 0, 0);
    sscrb_add_crossing(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrb_make_complex(0, { 0, 1 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));

    // A deterministic move that must survive untouched.
    complexList[0].trajTrans.x = 3.0;
    complexList[0].trajTrans.y = 0.0;
    complexList[0].trajTrans.z = 0.0;

    std::cerr << "  Calling sweep_separation_complex_rot_box(simItr=0, pro1Index=0)...\n";
    sweep_separation_complex_rot_box(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Both members translate by exactly +3 in x.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, 2.0) << "member 0 should move from -1 to +2";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 4.0) << "member 1 should move from +1 to +4";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, 0.0) << "no y displacement was requested";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.z, 0.0) << "no z displacement was requested";

    // Mass-weighted COM of the two members.
    EXPECT_DOUBLE_EQ(complexList[0].comCoord.x, 3.0) << "complex COM should be the average of (2, 4)";

    // Both members are flagged, and the buffers are cleared.
    EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::propagated) << "member 0 should be flagged propagated";
    EXPECT_EQ(moleculeList[1].trajStatus, TrajStatus::propagated) << "member 1 should be flagged propagated";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans must be zeroed on exit";

    std::cerr << "  member COMs after the call: " << moleculeList[0].comCoord.x << " and " << moleculeList[1].comCoord.x
              << " (expected 2 and 4)\n";
}

// -----------------------------------------------------------------------------
// Test 5: implicit-lipid partners are skipped outright.
// -----------------------------------------------------------------------------
void test_sscrb_implicit_lipid_partner_skipped()
{
    std::cerr << "\n[TEST] test_sscrb_implicit_lipid_partner_skipped\n"
              << "  Source file:   sweep_separation_complex_rot_box.cpp\n"
              << "  Function:      sweep_separation_complex_rot_box()\n"
              << "  Scenario:      the crossing partner is flagged isImplicitLipid and\n"
              << "                 its interface coincides with the target's interface.\n"
              << "  Pass criteria: the `continue` guard skips the pair entirely, so the\n"
              << "                 pre-set trajectory is applied verbatim and the lipid\n"
              << "                 complex is left alone.\n";

    Parameters params;
    params.timeStep = 0.1;

    Membrane membraneObject = sscrb_make_box(500.0);

    std::vector<MolTemplate> molTemplateList { sscrb_make_template(0, "A", 1.0) };
    std::vector<ForwardRxn> forwardRxns { sscrb_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrb_make_molecule(0, 0, 0, Coord { -1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrb_make_molecule(1, 1, 0, Coord { 1.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 0.0 }));
    moleculeList[1].isImplicitLipid = true; // <- the branch under test

    sscrb_add_crossing(moleculeList[0], 1, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrb_make_complex(0, { 0 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));
    complexList.push_back(
        sscrb_make_complex(1, { 1 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));

    complexList[0].trajTrans.x = 0.0;
    complexList[0].trajTrans.y = 4.0;
    complexList[0].trajTrans.z = 0.0;

    std::cerr << "  Calling sweep_separation_complex_rot_box(simItr=0, pro1Index=0)...\n";
    sweep_separation_complex_rot_box(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Deterministic +4 nm displacement in y proves no resampling happened.
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, 4.0) << "target should translate by exactly trajTrans.y = 4.0";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, -1.0) << "target x should be unchanged";

    // The implicit lipid complex is untouched.
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.y, 0.0) << "implicit lipid partner must never be resampled";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 1.0) << "implicit lipid molecule must not move";

    std::cerr << "  target COM after the call = (" << moleculeList[0].comCoord.x << ", " << moleculeList[0].comCoord.y
              << ", " << moleculeList[0].comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 6: a complex with no crossing partners at all is simply propagated, and
//         a non-zero rotation is applied rigidly and then cleared.
// -----------------------------------------------------------------------------
void test_sscrb_no_partners_rigid_propagation()
{
    std::cerr << "\n[TEST] test_sscrb_no_partners_rigid_propagation\n"
              << "  Source file:   sweep_separation_complex_rot_box.cpp\n"
              << "  Function:      sweep_separation_complex_rot_box()\n"
              << "  Scenario:      the target complex has an empty crossbase list and a\n"
              << "                 non-zero trajRot; its molecule sits exactly on the\n"
              << "                 complex COM.\n"
              << "  Pass criteria: the COM translates by exactly trajTrans (rotation about\n"
              << "                 the COM cannot move it), the interface stays rigidly at\n"
              << "                 1 nm from the COM, and both trajectory buffers end up\n"
              << "                 zeroed.\n";

    Parameters params;
    params.timeStep = 0.1;

    Membrane membraneObject = sscrb_make_box(500.0);

    std::vector<MolTemplate> molTemplateList { sscrb_make_template(0, "A", 1.0) };
    std::vector<ForwardRxn> forwardRxns { sscrb_make_rxn(1.0) };

    // Molecule COM == complex COM, interface exactly 1 nm away along +x.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrb_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrb_make_complex(0, { 0 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.01, 0.01, 0.01 }, molTemplateList.size()));

    complexList[0].trajTrans.x = 2.0;
    complexList[0].trajTrans.y = 3.0;
    complexList[0].trajTrans.z = -1.0;
    complexList[0].trajRot.x = 0.0;
    complexList[0].trajRot.y = 0.0;
    complexList[0].trajRot.z = 0.5; // a real rotation about z

    std::cerr << "  Calling sweep_separation_complex_rot_box(simItr=0, pro1Index=0)...\n";
    sweep_separation_complex_rot_box(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Rotating about the COM leaves the COM itself where it is, so the COM only
    // picks up the translation.
    EXPECT_NEAR(moleculeList[0].comCoord.x, 2.0, 1e-12) << "COM x should be the pure translation";
    EXPECT_NEAR(moleculeList[0].comCoord.y, 3.0, 1e-12) << "COM y should be the pure translation";
    EXPECT_NEAR(moleculeList[0].comCoord.z, -1.0, 1e-12) << "COM z should be the pure translation";

    // The rotation is rigid: the COM->interface distance is conserved.
    const double ifaceDist = sscrb_dist(moleculeList[0].interfaceList[0].coord, moleculeList[0].comCoord);
    std::cerr << "  COM->interface distance after rigid rotation = " << ifaceDist << " (expected 1.0)\n";
    EXPECT_NEAR(ifaceDist, 1.0, 1e-9) << "rotation must preserve the internal geometry";

    // A rotation about z really did happen (the interface left the +x axis).
    EXPECT_GT(std::abs(moleculeList[0].interfaceList[0].coord.y - moleculeList[0].comCoord.y), 1e-6)
        << "a non-zero trajRot.z should tilt the interface out of the x direction";

    // Both buffers zeroed on exit.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.z, 0.0) << "trajRot.z must be zeroed";
    EXPECT_EQ(moleculeList[0].trajStatus, TrajStatus::propagated) << "the lone molecule should be flagged propagated";
}

// -----------------------------------------------------------------------------
// Test 7: the target complex may be addressed through any of its members, and
//         the correct partner interface index is looked up from the reaction.
// -----------------------------------------------------------------------------
void test_sscrb_entry_via_second_member()
{
    std::cerr << "\n[TEST] test_sscrb_entry_via_second_member\n"
              << "  Source file:   sweep_separation_complex_rot_box.cpp\n"
              << "  Function:      sweep_separation_complex_rot_box()\n"
              << "  Scenario:      a two-molecule complex is entered through its SECOND\n"
              << "                 member; only the first member has a (far away)\n"
              << "                 crossing partner.\n"
              << "  Pass criteria: the routine still resolves the owning complex from\n"
              << "                 pro1Index, finds no overlap, and translates every\n"
              << "                 member by exactly the stored trajTrans.\n";

    Parameters params;
    params.timeStep = 0.1;

    Membrane membraneObject = sscrb_make_box(500.0);

    std::vector<MolTemplate> molTemplateList { sscrb_make_template(0, "A", 1.0) };
    std::vector<ForwardRxn> forwardRxns { sscrb_make_rxn(1.0) };

    // Molecules 0 and 1 form complex 0; molecule 2 is a distant partner (complex 1).
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrb_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 1.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrb_make_molecule(1, 0, 0, Coord { 4.0, 0.0, 0.0 }, Coord { 5.0, 0.0, 0.0 }));
    moleculeList.push_back(sscrb_make_molecule(2, 1, 0, Coord { 80.0, 0.0, 0.0 }, Coord { 79.0, 0.0, 0.0 }));

    sscrb_add_crossing(moleculeList[0], 2, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrb_make_complex(0, { 0, 1 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));
    complexList.push_back(
        sscrb_make_complex(1, { 2 }, moleculeList, Coord { 1.0, 1.0, 1.0 }, Coord { 0.0, 0.0, 0.0 }, molTemplateList.size()));

    complexList[0].trajTrans.x = -2.0;
    complexList[0].trajTrans.y = 1.0;
    complexList[0].trajTrans.z = 0.0;

    std::cerr << "  Calling sweep_separation_complex_rot_box(simItr=0, pro1Index=1)...\n";
    // Entering through member index 1 rather than 0.
    sweep_separation_complex_rot_box(0, 1, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.x, -2.0) << "member 0 should move from 0 to -2";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.x, 2.0) << "member 1 should move from 4 to 2";
    EXPECT_DOUBLE_EQ(moleculeList[0].comCoord.y, 1.0) << "member 0 should move +1 in y";
    EXPECT_DOUBLE_EQ(moleculeList[1].comCoord.y, 1.0) << "member 1 should move +1 in y";

    // The distant partner is untouched.
    EXPECT_DOUBLE_EQ(moleculeList[2].comCoord.x, 80.0) << "the distant partner must not move";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.x, 0.0) << "the distant partner must not be resampled";

    // Complex COM is the average of the two members.
    EXPECT_DOUBLE_EQ(complexList[0].comCoord.x, 0.0) << "complex COM should be the average of (-2, 2)";

    std::cerr << "  complex 0 COM after the call = (" << complexList[0].comCoord.x << ", " << complexList[0].comCoord.y
              << ", " << complexList[0].comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so that
// a failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(SweepSeparationComplexRotBox, NoOverlapAppliesTrajectory) { test_sscrb_no_overlap_applies_trajectory(); }
TEST(SweepSeparationComplexRotBox, OverlapResamplesBothComplexes) { test_sscrb_overlap_resamples_both_complexes(); }
TEST(SweepSeparationComplexRotBox, PropagatedPartnerIsNotResampled) { test_sscrb_propagated_partner_is_not_resampled(); }
TEST(SweepSeparationComplexRotBox, SameComplexPartnersSkipped) { test_sscrb_same_complex_partners_skipped(); }
TEST(SweepSeparationComplexRotBox, ImplicitLipidPartnerSkipped) { test_sscrb_implicit_lipid_partner_skipped(); }
TEST(SweepSeparationComplexRotBox, NoPartnersRigidPropagation) { test_sscrb_no_partners_rigid_propagation(); }
TEST(SweepSeparationComplexRotBox, EntryViaSecondMember) { test_sscrb_entry_via_second_member(); }