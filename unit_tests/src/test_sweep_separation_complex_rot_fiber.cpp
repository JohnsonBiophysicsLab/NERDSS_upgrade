/*! \file test_sweep_separation_complex_rot_fiber.cpp
 *
 * ### Unit test for src/trajectory_functions/sweep_separation_complex_rot_fiber.cpp
 *
 * The file under test contains a single, very thin dispatcher:
 *
 *     void sweep_separation_complex_rot_fiber(int simItr, int pro1Index,
 *          Parameters&, std::vector<Molecule>&, std::vector<Complex>&,
 *          const std::vector<ForwardRxn>&, const std::vector<MolTemplate>&,
 *          const Membrane&)
 *
 * The sphere branch is currently commented out in the source, so every call is
 * forwarded to sweep_separation_complex_rot_fiber_box().  The contract of that
 * routine is: "protein pro1Index is about to be displaced by its complex'
 * trajTrans/trajRot.  If, at the new position, it would overlap (get closer
 * than the binding radius of the relevant reaction) any of the partners listed
 * in its crossbase list, the displacement is rejected and new displacements are
 * drawn until no overlap remains."
 *
 * Because the routine is stochastic, the tests below are written around the
 * pieces of behaviour that are deterministic and observable from outside:
 *
 *   1. A molecule with *no* crossing partners and a zero trajectory must be a
 *      no-op (nothing to sweep against, nothing to move).
 *   2. Two well separated partners whose trajectories have already been fixed
 *      (TrajStatus::propagated, zero displacement) must be left alone.
 *   3. Two partners that currently overlap and whose trajectories may still be
 *      resampled (TrajStatus::none) must end up with a *proposed* final
 *      separation that is at least the binding radius.
 *   4. A multi-molecule complex must survive the sweep with its bookkeeping
 *      (memberList / myComIndex / isEmpty) intact and with finite coordinates.
 *
 * Note on how "final position" is measured: depending on whether the routine
 * only fixes up the trajectory vectors or actually propagates the complexes,
 * the answer lives either in comCoord+trajTrans or in comCoord alone (with
 * trajTrans zeroed).  Every assertion below therefore uses (comCoord+trajTrans)
 * which is correct in both cases.
 */

#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// The RNG used by GaussV()/rand_gsl() is a global owned by gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// RNG bootstrap.  The sweep resamples displacements with GaussV(), which needs a
// live gsl_rng.  Allocate one if the suite has not done so yet and always seed
// deterministically so the test is reproducible.
// -----------------------------------------------------------------------------
void sscrf_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    } else {
        gsl_rng_set(r, 42);
    }
}

// -----------------------------------------------------------------------------
// Build the single MolTemplate shared by all tests: a point-like protein that
// sits on a fiber, i.e. it only diffuses along x and does not rotate.
// -----------------------------------------------------------------------------
MolTemplate sscrf_make_template()
{
    MolTemplate temp;
    temp.molTypeIndex = 0;
    temp.molName = "fiberProtein";
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord { 10.0, 0.0, 0.0 }; // 1-D diffusion along the fiber (nm^2/us)
    temp.Dr = Coord { 0.0, 0.0, 0.0 }; // a point on a fiber does not rotate
    temp.isPromoter = true;            // "lives on a fiber"
    temp.isPoint = true;
    temp.isLipid = false;
    temp.isRod = false;
    temp.checkOverlap = false;
    temp.copies = 2;

    // One interface, coincident with the centre of mass (keeps the geometry
    // trivial: COM separation == interface separation).
    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord { 0.0, 0.0, 0.0 };
    iface.stateList.push_back(Interface::State(std::string("a"), 0));
    temp.interfaceList.push_back(iface);

    return temp;
}

// -----------------------------------------------------------------------------
// Build a fully initialised Molecule sitting at (x, 0, 0) on the fiber.
// Every field the sweep (or a downstream propagate()) might read is filled in.
// -----------------------------------------------------------------------------
Molecule sscrf_make_molecule(int index, int comIndex, double x)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mySubVolIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = true;
    mol.linksToSurface = 0;
    mol.comCoord = Coord { x, 0.0, 0.0 };

    Molecule::Iface iface;
    iface.coord = mol.comCoord; // interface sits on the COM
    iface.index = 0;
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    mol.freelist.push_back(0); // the single interface is free
    mol.trajStatus = TrajStatus::none;
    return mol;
}

// -----------------------------------------------------------------------------
// Build a Complex holding the supplied member molecules.  The complex starts
// with a zero trajectory so that, even if the routine ends by propagating, the
// coordinates are unchanged.
// -----------------------------------------------------------------------------
Complex sscrf_make_complex(int index, const std::vector<Molecule>& moleculeList,
    const std::vector<int>& members)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.isEmpty = false;
    com.OnSurface = false;
    com.onFiber = true;
    com.linksToSurface = 0;
    com.iLipidIndex = 0;
    com.ncross = 0;
    com.radius = 1.0;
    com.D = Coord { 10.0, 0.0, 0.0 };
    com.Dr = Coord { 0.0, 0.0, 0.0 };
    com.trajStatus = TrajStatus::none;

    com.trajTrans.x = 0.0;
    com.trajTrans.y = 0.0;
    com.trajTrans.z = 0.0;
    com.trajRot.x = 0.0;
    com.trajRot.y = 0.0;
    com.trajRot.z = 0.0;

    // mass weighted centre of mass over the member molecules
    double totMass { 0.0 };
    Coord sum { 0.0, 0.0, 0.0 };
    for (int m : members) {
        com.memberList.push_back(m);
        totMass += moleculeList[m].mass;
        sum.x += moleculeList[m].comCoord.x * moleculeList[m].mass;
        sum.y += moleculeList[m].comCoord.y * moleculeList[m].mass;
        sum.z += moleculeList[m].comCoord.z * moleculeList[m].mass;
    }
    com.mass = totMass;
    com.comCoord = Coord { sum.x / totMass, sum.y / totMass, sum.z / totMass };

    com.numEachMol.assign(1, static_cast<int>(members.size()));
    com.lastNumberUpdateItrEachMol.assign(1, 0);

    return com;
}

// -----------------------------------------------------------------------------
// Register two molecules as mutual "crossing partners", i.e. they are close
// enough that the sweep has to consider them.  Reaction 0 (index stored in
// crossrxn) supplies the binding radius used as the overlap criterion.
// -----------------------------------------------------------------------------
void sscrf_link_cross_partners(std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, int p1, int p2, double prob)
{
    moleculeList[p1].crossbase.push_back(p2);
    moleculeList[p1].mycrossint.push_back(0); // our only interface
    moleculeList[p1].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    moleculeList[p1].probvec.push_back(prob);

    moleculeList[p2].crossbase.push_back(p1);
    moleculeList[p2].mycrossint.push_back(0);
    moleculeList[p2].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    moleculeList[p2].probvec.push_back(prob);

    ++complexList[moleculeList[p1].myComIndex].ncross;
    ++complexList[moleculeList[p2].myComIndex].ncross;
}

// -----------------------------------------------------------------------------
// A single bimolecular reaction whose binding radius is the overlap threshold.
// -----------------------------------------------------------------------------
ForwardRxn sscrf_make_forward_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = false;
    rxn.conjBackRxnIndex = -1;
    rxn.bindRadius = bindRadius;
    rxn.bindRadius2D = bindRadius;
    rxn.isOnMem = false;
    rxn.isSymmetric = true;

    // reactant/product interface descriptors (name, molType, abs, rel, state, bound)
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.productListNew.emplace_back("a", 0, 0, 0, '\0', true);
    rxn.productListNew.emplace_back("a", 0, 0, 0, '\0', true);

    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = 1.0;

    return rxn;
}

// -----------------------------------------------------------------------------
// A reflecting cubic water box of side `side` centred on the origin.
// -----------------------------------------------------------------------------
Membrane sscrf_make_membrane(double side)
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.implicitLipid = false;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { side, side, side });
    membraneObject.xBCtype = "reflect";
    membraneObject.yBCtype = "reflect";
    membraneObject.zBCtype = "reflect";
    // The 3D reflecting-surface lookup table is only consulted for implicit
    // lipids, but give it a valid size so no stray index can run off the end.
    membraneObject.RS3Dvect.assign(500, 0.0);
    membraneObject.nSites = 0;
    membraneObject.No_free_lipids = 0;
    membraneObject.No_protein = 0;
    return membraneObject;
}

// -----------------------------------------------------------------------------
// Simulation parameters: 1 us timestep gives sqrt(2*D*dt) ~ 4.5 nm steps, which
// is comfortably larger than the 1 nm binding radius used in the overlap test.
// -----------------------------------------------------------------------------
Parameters sscrf_make_params()
{
    Parameters params;
    params.rank = 0;
    params.numMolTypes = 1;
    params.numTotalSpecies = 1;
    params.numTotalComplex = 2;
    params.timeStep = 1.0;
    params.nItr = 10;
    params.overlapSepLimit = 0.1;
    params.scaleMaxDisplace = 100.0;
    params.clusterOverlapCheck = false;
    params.rMaxLimit = 20.0;
    params.rMaxRadius = 2.0;
    return params;
}

//! \brief Proposed final position of a molecule = current position + pending translation.
Coord sscrf_final_position(const Molecule& mol, const Complex& com)
{
    return Coord { mol.comCoord.x + com.trajTrans.x, mol.comCoord.y + com.trajTrans.y,
        mol.comCoord.z + com.trajTrans.z };
}

//! \brief Euclidean distance helper.
double sscrf_distance(const Coord& a, const Coord& b)
{
    const double dx { a.x - b.x };
    const double dy { a.y - b.y };
    const double dz { a.z - b.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

//! \brief True if every component of the coordinate is a finite number.
bool sscrf_is_finite(const Coord& c)
{
    return std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z);
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: an isolated molecule (no crossing partners at all).
//
// There is nothing to sweep against, so the dispatcher must simply return and
// leave the system exactly as it found it.
// -----------------------------------------------------------------------------
void test_sscrf_isolated_molecule_is_noop()
{
    std::cerr << "\n[TEST] test_sscrf_isolated_molecule_is_noop\n"
              << "  Source file:   sweep_separation_complex_rot_fiber.cpp\n"
              << "  Function:      sweep_separation_complex_rot_fiber\n"
              << "  Scenario:      one molecule, empty crossbase, zero trajectory.\n"
              << "  Pass criteria: coordinates unchanged and finite, complex/molecule\n"
              << "                 bookkeeping untouched.\n";

    sscrf_init_rng();

    Parameters params { sscrf_make_params() };
    Membrane membraneObject { sscrf_make_membrane(100.0) };
    std::vector<MolTemplate> molTemplateList { sscrf_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrf_make_forward_rxn(1.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrf_make_molecule(0, 0, 0.0));

    std::vector<Complex> complexList;
    complexList.push_back(sscrf_make_complex(0, moleculeList, std::vector<int> { 0 }));

    // No crossing partners at all.
    EXPECT_TRUE(moleculeList[0].crossbase.empty())
        << "precondition: the isolated molecule must have no crossing partners";

    const Coord before { moleculeList[0].comCoord };
    std::cerr << "  Before: COM = (" << before.x << ", " << before.y << ", " << before.z << ")\n";

    sweep_separation_complex_rot_fiber(0, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    const Coord after { sscrf_final_position(moleculeList[0], complexList[0]) };
    std::cerr << "  After : proposed COM = (" << after.x << ", " << after.y << ", " << after.z << ")\n";

    // With no partners and a zero trajectory nothing may move.
    EXPECT_TRUE(sscrf_is_finite(after)) << "coordinates must remain finite";
    EXPECT_NEAR(after.x, before.x, 1e-9) << "x must not change for an isolated molecule";
    EXPECT_NEAR(after.y, before.y, 1e-9) << "y must not change for an isolated molecule";
    EXPECT_NEAR(after.z, before.z, 1e-9) << "z must not change for an isolated molecule";

    // Bookkeeping must be preserved.
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "molecule must still belong to complex 0";
    EXPECT_EQ(complexList[0].memberList.size(), 1u) << "complex must still own exactly one molecule";
    EXPECT_FALSE(moleculeList[0].isEmpty) << "the molecule must not be marked empty";
    EXPECT_FALSE(complexList[0].isEmpty) << "the complex must not be marked empty";
}

// -----------------------------------------------------------------------------
// Test 2: two crossing partners that are nowhere near each other.
//
// Both trajectories have already been decided (TrajStatus::propagated) and are
// zero, so the sweep has no reason to resample anything: the pair is 40 nm
// apart while the binding radius is only 1 nm.
// -----------------------------------------------------------------------------
void test_sscrf_separated_pair_left_alone()
{
    std::cerr << "\n[TEST] test_sscrf_separated_pair_left_alone\n"
              << "  Source file:   sweep_separation_complex_rot_fiber.cpp\n"
              << "  Function:      sweep_separation_complex_rot_fiber\n"
              << "  Scenario:      two cross partners 40 nm apart, bindRadius = 1 nm,\n"
              << "                 trajectories already fixed and equal to zero.\n"
              << "  Pass criteria: no displacement is introduced; separation unchanged.\n";

    sscrf_init_rng();

    const double bindRadius { 1.0 };

    Parameters params { sscrf_make_params() };
    Membrane membraneObject { sscrf_make_membrane(100.0) };
    std::vector<MolTemplate> molTemplateList { sscrf_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrf_make_forward_rxn(bindRadius) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrf_make_molecule(0, 0, -20.0));
    moleculeList.push_back(sscrf_make_molecule(1, 1, +20.0));

    std::vector<Complex> complexList;
    complexList.push_back(sscrf_make_complex(0, moleculeList, std::vector<int> { 0 }));
    complexList.push_back(sscrf_make_complex(1, moleculeList, std::vector<int> { 1 }));

    sscrf_link_cross_partners(moleculeList, complexList, 0, 1, 0.0);

    // Freeze the trajectories: "already propagated" means the sweep must not
    // draw new random displacements for these complexes.
    for (auto& mol : moleculeList)
        mol.trajStatus = TrajStatus::propagated;
    for (auto& com : complexList)
        com.trajStatus = TrajStatus::propagated;

    const double sepBefore { sscrf_distance(moleculeList[0].comCoord, moleculeList[1].comCoord) };
    std::cerr << "  Before: separation = " << sepBefore << " nm (bindRadius = " << bindRadius << " nm)\n";

    sweep_separation_complex_rot_fiber(0, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    const Coord p0 { sscrf_final_position(moleculeList[0], complexList[0]) };
    const Coord p1 { sscrf_final_position(moleculeList[1], complexList[1]) };
    const double sepAfter { sscrf_distance(p0, p1) };
    std::cerr << "  After : separation = " << sepAfter << " nm\n";
    std::cerr << "          trajTrans[0] = (" << complexList[0].trajTrans.x << ", "
              << complexList[0].trajTrans.y << ", " << complexList[0].trajTrans.z << ")\n";

    EXPECT_TRUE(sscrf_is_finite(p0)) << "molecule 0 coordinates must be finite";
    EXPECT_TRUE(sscrf_is_finite(p1)) << "molecule 1 coordinates must be finite";

    // Nothing overlapped, so the well separated pair must be untouched.
    EXPECT_NEAR(sepAfter, sepBefore, 1e-9)
        << "a non-overlapping pair with fixed trajectories must not be moved";
    EXPECT_NEAR(p0.x, -20.0, 1e-9) << "molecule 0 must stay at x = -20";
    EXPECT_NEAR(p1.x, +20.0, 1e-9) << "molecule 1 must stay at x = +20";

    // And, trivially, the pair still satisfies the non-overlap condition.
    EXPECT_GE(sepAfter, bindRadius) << "separation must remain at least the binding radius";
}

// -----------------------------------------------------------------------------
// Test 3: two crossing partners that currently overlap.
//
// Their trajectories may still be resampled (TrajStatus::none), which is
// exactly the situation the sweep exists for.  After the call the *proposed*
// positions must no longer be inside the binding radius of each other.
// -----------------------------------------------------------------------------
void test_sscrf_overlapping_pair_is_separated()
{
    std::cerr << "\n[TEST] test_sscrf_overlapping_pair_is_separated\n"
              << "  Source file:   sweep_separation_complex_rot_fiber.cpp\n"
              << "  Function:      sweep_separation_complex_rot_fiber\n"
              << "  Scenario:      two cross partners only 0.4 nm apart while the\n"
              << "                 binding radius is 1.0 nm; trajectories are still\n"
              << "                 resampleable (TrajStatus::none).\n"
              << "  Pass criteria: the routine terminates, coordinates stay finite and\n"
              << "                 inside the box, and the proposed final separation is\n"
              << "                 at least the binding radius.\n";

    sscrf_init_rng();

    const double bindRadius { 1.0 };
    const double boxSide { 100.0 };

    Parameters params { sscrf_make_params() };
    Membrane membraneObject { sscrf_make_membrane(boxSide) };
    std::vector<MolTemplate> molTemplateList { sscrf_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrf_make_forward_rxn(bindRadius) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrf_make_molecule(0, 0, 0.0));
    moleculeList.push_back(sscrf_make_molecule(1, 1, 0.4)); // overlapping!

    std::vector<Complex> complexList;
    complexList.push_back(sscrf_make_complex(0, moleculeList, std::vector<int> { 0 }));
    complexList.push_back(sscrf_make_complex(1, moleculeList, std::vector<int> { 1 }));

    sscrf_link_cross_partners(moleculeList, complexList, 0, 1, 0.0);

    // Leave both molecules/complexes resampleable so the sweep is allowed to
    // draw new displacements for them.
    for (auto& mol : moleculeList)
        mol.trajStatus = TrajStatus::none;
    for (auto& com : complexList)
        com.trajStatus = TrajStatus::none;

    const double sepBefore { sscrf_distance(moleculeList[0].comCoord, moleculeList[1].comCoord) };
    std::cerr << "  Before: separation = " << sepBefore << " nm (< bindRadius = " << bindRadius << " nm)\n";
    EXPECT_LT(sepBefore, bindRadius) << "precondition: the pair must start out overlapping";

    sweep_separation_complex_rot_fiber(0, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    const Coord p0 { sscrf_final_position(moleculeList[0], complexList[0]) };
    const Coord p1 { sscrf_final_position(moleculeList[1], complexList[1]) };
    const double sepAfter { sscrf_distance(p0, p1) };

    std::cerr << "  After : molecule 0 proposed COM = (" << p0.x << ", " << p0.y << ", " << p0.z << ")\n";
    std::cerr << "          molecule 1 proposed COM = (" << p1.x << ", " << p1.y << ", " << p1.z << ")\n";
    std::cerr << "          separation = " << sepAfter << " nm\n";

    // Basic sanity: nothing blew up numerically.
    EXPECT_TRUE(sscrf_is_finite(p0)) << "molecule 0 coordinates must be finite after the sweep";
    EXPECT_TRUE(sscrf_is_finite(p1)) << "molecule 1 coordinates must be finite after the sweep";

    // Reflecting boundaries: the proposed positions must remain in the box.
    const double half { boxSide / 2.0 };
    EXPECT_LE(std::abs(p0.x), half + 1e-6) << "molecule 0 must stay inside the box in x";
    EXPECT_LE(std::abs(p1.x), half + 1e-6) << "molecule 1 must stay inside the box in x";

    // The whole point of the routine: the overlap must be resolved.
    EXPECT_GE(sepAfter, bindRadius - 1e-6)
        << "after sweeping, the proposed separation must be at least the binding radius";
}

// -----------------------------------------------------------------------------
// Test 4: a two-molecule complex crossing a molecule in a second complex.
//
// This exercises the loop over the members of the swept complex.  The partner
// is far away and every trajectory is frozen, so the expected outcome is again
// "nothing changes" - what we are really checking is that the routine walks a
// multi-member memberList without corrupting the molecule/complex mapping.
// -----------------------------------------------------------------------------
void test_sscrf_multimember_complex_state_preserved()
{
    std::cerr << "\n[TEST] test_sscrf_multimember_complex_state_preserved\n"
              << "  Source file:   sweep_separation_complex_rot_fiber.cpp\n"
              << "  Function:      sweep_separation_complex_rot_fiber\n"
              << "  Scenario:      complex 0 holds two molecules, both listed as cross\n"
              << "                 partners of a distant molecule in complex 1.\n"
              << "  Pass criteria: memberList / myComIndex / isEmpty are preserved and\n"
              << "                 all coordinates remain finite and unchanged.\n";

    sscrf_init_rng();

    const double bindRadius { 1.0 };

    Parameters params { sscrf_make_params() };
    Membrane membraneObject { sscrf_make_membrane(100.0) };
    std::vector<MolTemplate> molTemplateList { sscrf_make_template() };
    std::vector<ForwardRxn> forwardRxns { sscrf_make_forward_rxn(bindRadius) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrf_make_molecule(0, 0, -10.0)); // complex 0, member 1
    moleculeList.push_back(sscrf_make_molecule(1, 0, -8.0));  // complex 0, member 2
    moleculeList.push_back(sscrf_make_molecule(2, 1, 25.0));  // complex 1, far away

    std::vector<Complex> complexList;
    complexList.push_back(sscrf_make_complex(0, moleculeList, std::vector<int> { 0, 1 }));
    complexList.push_back(sscrf_make_complex(1, moleculeList, std::vector<int> { 2 }));

    // Both members of complex 0 "see" the distant molecule 2.
    sscrf_link_cross_partners(moleculeList, complexList, 0, 2, 0.0);
    sscrf_link_cross_partners(moleculeList, complexList, 1, 2, 0.0);

    // Freeze all trajectories so the expected result is deterministic.
    for (auto& mol : moleculeList)
        mol.trajStatus = TrajStatus::propagated;
    for (auto& com : complexList)
        com.trajStatus = TrajStatus::propagated;

    const Coord before0 { moleculeList[0].comCoord };
    const Coord before1 { moleculeList[1].comCoord };
    const Coord before2 { moleculeList[2].comCoord };
    std::cerr << "  Before: complex 0 members at x = " << before0.x << " and " << before1.x
              << ", complex 1 member at x = " << before2.x << '\n';

    sweep_separation_complex_rot_fiber(0, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    const Coord after0 { sscrf_final_position(moleculeList[0], complexList[0]) };
    const Coord after1 { sscrf_final_position(moleculeList[1], complexList[0]) };
    const Coord after2 { sscrf_final_position(moleculeList[2], complexList[1]) };
    std::cerr << "  After : complex 0 members at x = " << after0.x << " and " << after1.x
              << ", complex 1 member at x = " << after2.x << '\n';

    // Coordinates: finite and, because everything was frozen and separated,
    // exactly where they started.
    EXPECT_TRUE(sscrf_is_finite(after0)) << "molecule 0 coordinates must be finite";
    EXPECT_TRUE(sscrf_is_finite(after1)) << "molecule 1 coordinates must be finite";
    EXPECT_TRUE(sscrf_is_finite(after2)) << "molecule 2 coordinates must be finite";
    EXPECT_NEAR(after0.x, before0.x, 1e-9) << "member 0 must not move";
    EXPECT_NEAR(after1.x, before1.x, 1e-9) << "member 1 must not move";
    EXPECT_NEAR(after2.x, before2.x, 1e-9) << "the distant partner must not move";

    // Bookkeeping: the sweep must never re-assign molecules between complexes.
    ASSERT_EQ(complexList[0].memberList.size(), 2u)
        << "complex 0 must still own exactly two molecules";
    EXPECT_EQ(complexList[0].memberList[0], 0) << "complex 0 must still own molecule 0";
    EXPECT_EQ(complexList[0].memberList[1], 1) << "complex 0 must still own molecule 1";
    EXPECT_EQ(complexList[1].memberList.size(), 1u)
        << "complex 1 must still own exactly one molecule";
    EXPECT_EQ(moleculeList[0].myComIndex, 0) << "molecule 0 must still point at complex 0";
    EXPECT_EQ(moleculeList[1].myComIndex, 0) << "molecule 1 must still point at complex 0";
    EXPECT_EQ(moleculeList[2].myComIndex, 1) << "molecule 2 must still point at complex 1";

    for (const auto& mol : moleculeList)
        EXPECT_FALSE(mol.isEmpty) << "no molecule may be destroyed by the sweep";
    for (const auto& com : complexList)
        EXPECT_FALSE(com.isEmpty) << "no complex may be destroyed by the sweep";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario lives in its own TEST so that a failure in
// one does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(SweepSeparationComplexRotFiberTest, IsolatedMoleculeIsNoop)
{
    test_sscrf_isolated_molecule_is_noop();
}

TEST(SweepSeparationComplexRotFiberTest, SeparatedPairLeftAlone)
{
    test_sscrf_separated_pair_left_alone();
}

TEST(SweepSeparationComplexRotFiberTest, OverlappingPairIsSeparated)
{
    test_sscrf_overlapping_pair_is_separated();
}

TEST(SweepSeparationComplexRotFiberTest, MultiMemberComplexStatePreserved)
{
    test_sscrf_multimember_complex_state_preserved();
}