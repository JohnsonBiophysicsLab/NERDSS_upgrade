/*! \file test_sweep_separation_complex_rot.cpp
 *
 * ### Unit test for src/trajectory_functions/sweep_separation_complex_rot.cpp
 *
 * The single function defined in that file is
 *
 *     void sweep_separation_complex_rot(int simItr, int pro1Index, Parameters& params,
 *                                       std::vector<Molecule>& moleculeList,
 *                                       std::vector<Complex>& complexList,
 *                                       const std::vector<ForwardRxn>& forwardRxns,
 *                                       const std::vector<MolTemplate>& molTemplateList,
 *                                       const Membrane& membraneObject)
 *
 * It is a dispatcher: when the Membrane describes a sphere it forwards to
 * sweep_separation_complex_rot_sphere(), otherwise to
 * sweep_separation_complex_rot_box().  Both of those routines
 *
 *   1. look at every "cross" partner recorded on the molecules of the target
 *      complex,
 *   2. reject and re-sample the trajectory of the target complex (and of any
 *      partner complex that has not moved yet) whenever the proposed step
 *      would place two reacting interfaces closer than the reaction binding
 *      radius,
 *   3. finally propagate the target complex and zero out its stored
 *      translation/rotation vectors.
 *
 * The tests below therefore verify observable behaviour of the dispatcher:
 *   * a complex with no cross partners is simply propagated by its stored
 *     trajTrans (box branch and sphere branch),
 *   * partners that live inside the *same* complex are ignored (they cannot
 *     diffuse relative to one another) so no re-sampling happens,
 *   * genuinely overlapping partners in *different* complexes cause the
 *     trajectories to be re-sampled until the interfaces are separated by at
 *     least the binding radius,
 *   * only the target complex is propagated; the partner complex keeps its
 *     (re-sampled but not yet applied) trajectory.
 *
 * Verbose progress information is written to stderr so a reader can follow
 * exactly which source file / function is under test and what each assertion
 * is checking.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// The RNG pointer is defined in unit_tests/src/gtest_main.cpp -- only declare it here.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Small helpers used to build the fully initialised objects that the routines
// under test dereference.  Every field that the sweep / propagate / reflect
// code path reads is set explicitly; under-initialised objects would crash the
// whole gtest binary.
// -----------------------------------------------------------------------------

/*! \brief Initialise (or re-seed) the GSL random number generator.
 *
 * The overlap-resolution branch of the sweep calls GaussV(), which uses the
 * global generator `r`.  A fixed seed keeps the test reproducible.
 */
void sscr_init_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief One molecule type: a point-like molecule with a single interface
 *         that sits on top of its centre of mass.
 *
 * Rotational diffusion is deliberately zero so that Complex::propagate() takes
 * its pure-translation fast path; this makes the expected coordinates exact.
 */
MolTemplate sscr_make_moltemplate()
{
    MolTemplate temp;
    temp.molName = "A";
    temp.molTypeIndex = 0;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord { 10.0, 10.0, 10.0 }; // large enough for overlap to be resolved fast
    temp.Dr = Coord { 0.0, 0.0, 0.0 }; // no rotation -> deterministic translation
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPromoter = false;
    temp.isPoint = true;

    Interface iface { "a", Coord { 0.0, 0.0, 0.0 } };
    iface.index = 0;
    iface.stateList.emplace_back('\0', 0); // a single, default (unbound) state
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a molecule of type 0 with one interface coincident with the COM. */
Molecule sscr_make_molecule(int index, int comIndex, const Coord& com)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.comCoord = com;
    mol.trajStatus = TrajStatus::none;

    Molecule::Iface iface;
    iface.coord = com;
    iface.index = 0;
    iface.relIndex = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    mol.freelist.push_back(0);

    return mol;
}

/*! \brief Build a complex owning the listed member molecules. */
Complex sscr_make_complex(int index, const std::vector<int>& members, const Coord& com, double transD)
{
    Complex targCom;
    targCom.index = index;
    targCom.id = index;
    targCom.comCoord = com;
    targCom.radius = 1.0;
    targCom.mass = static_cast<double>(members.size());
    targCom.memberList = members;
    targCom.numEachMol = std::vector<int> { static_cast<int>(members.size()) };
    targCom.lastNumberUpdateItrEachMol = std::vector<long long int>(1, 0);
    targCom.D = Coord { transD, transD, transD };
    targCom.Dr = Coord { 0.0, 0.0, 0.0 }; // no rotation
    targCom.isEmpty = false;
    targCom.OnSurface = false;
    targCom.onFiber = false;
    targCom.linksToSurface = 0;
    targCom.ncross = 0;
    targCom.trajStatus = TrajStatus::none;
    targCom.trajTrans = Vector { 0.0, 0.0, 0.0 };
    targCom.trajRot = Coord { 0.0, 0.0, 0.0 };
    return targCom;
}

/*! \brief A single symmetric bimolecular reaction between interface 0 of both
 *         partners, with the requested binding radius. */
ForwardRxn sscr_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.isReversible = false;
    // RxnIface(name, molTypeIndex, absIfaceIndex, relIfaceIndex, requiresState, requiresInteraction)
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.productListNew = rxn.reactantListNew;
    return rxn;
}

/*! \brief A large cubic water box centred on the origin (reflecting box case). */
Membrane sscr_make_box_membrane()
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.implicitLipid = false;
    membraneObject.hasCompartment = false;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
    membraneObject.xBCtype = "reflect";
    membraneObject.yBCtype = "reflect";
    membraneObject.zBCtype = "reflect";
    // The sweep looks up a reflecting-surface value in this table; make sure it
    // is large enough to be indexed safely (entries 300-499 are read).
    membraneObject.RS3Dvect.assign(500, 0.0);
    return membraneObject;
}

/*! \brief A large spherical boundary centred on the origin (sphere case). */
Membrane sscr_make_sphere_membrane()
{
    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.isBox = false;
    membraneObject.implicitLipid = false;
    membraneObject.hasCompartment = false;
    membraneObject.sphereR = 100.0;
    membraneObject.sphereVol = (4.0 / 3.0) * M_PI * std::pow(membraneObject.sphereR, 3.0);
    // Several sphere routines still consult the enclosing water box.
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
    membraneObject.RS3Dvect.assign(500, 0.0);
    return membraneObject;
}

/*! \brief Euclidean distance between two coordinates. */
double sscr_distance(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: box boundary, target complex has no cross partners.
//         The complex must simply be propagated by its stored trajTrans and
//         the trajectory vectors must be zeroed afterwards.
// -----------------------------------------------------------------------------
void test_sscr_box_no_partners_translates_complex()
{
    std::cerr << "\n[TEST] test_sscr_box_no_partners_translates_complex\n"
              << "  Source file: sweep_separation_complex_rot.cpp\n"
              << "  Function:    sweep_separation_complex_rot (box branch)\n"
              << "  Scenario:    one complex, one molecule, no cross partners.\n"
              << "  Criteria:    molecule + complex move by exactly trajTrans,\n"
              << "               trajTrans/trajRot are zeroed, molecule is flagged\n"
              << "               as propagated.\n";

    sscr_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscr_make_box_membrane();
    std::vector<MolTemplate> molTemplateList { sscr_make_moltemplate() };
    std::vector<ForwardRxn> forwardRxns { sscr_make_rxn(1.0) };

    // One molecule at the origin inside its own complex.
    std::vector<Molecule> moleculeList { sscr_make_molecule(0, 0, Coord { 0.0, 0.0, 0.0 }) };
    std::vector<Complex> complexList { sscr_make_complex(0, std::vector<int> { 0 }, Coord { 0.0, 0.0, 0.0 }, 10.0) };

    // A modest step that stays well inside the 200 nm box.
    const Vector step { 2.0, -3.0, 1.5 };
    complexList[0].trajTrans = step;

    std::cerr << "  Proposed step = (" << step.x << ", " << step.y << ", " << step.z << ")\n";
    std::cerr << "  Calling sweep_separation_complex_rot...\n";
    sweep_separation_complex_rot(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // The molecule (and hence the complex COM) should have shifted by the step.
    EXPECT_NEAR(moleculeList[0].comCoord.x, step.x, 1e-10)
        << "molecule COM x should be displaced by trajTrans.x";
    EXPECT_NEAR(moleculeList[0].comCoord.y, step.y, 1e-10)
        << "molecule COM y should be displaced by trajTrans.y";
    EXPECT_NEAR(moleculeList[0].comCoord.z, step.z, 1e-10)
        << "molecule COM z should be displaced by trajTrans.z";

    // The interface rides along with the molecule.
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.x, step.x, 1e-10)
        << "interface x should follow the molecule";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.y, step.y, 1e-10)
        << "interface y should follow the molecule";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.z, step.z, 1e-10)
        << "interface z should follow the molecule";

    // Complex COM is recomputed from its members during propagation.
    EXPECT_NEAR(complexList[0].comCoord.x, step.x, 1e-10) << "complex COM x should track its member";
    EXPECT_NEAR(complexList[0].comCoord.y, step.y, 1e-10) << "complex COM y should track its member";
    EXPECT_NEAR(complexList[0].comCoord.z, step.z, 1e-10) << "complex COM z should track its member";

    // Trajectory vectors must be reset so the position does not change again.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "trajTrans.y must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "trajTrans.z must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.x, 0.0) << "trajRot.x must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.y, 0.0) << "trajRot.y must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.z, 0.0) << "trajRot.z must be zeroed after propagation";

    // Propagation marks the members so they are not moved twice in one step.
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "member molecule should be flagged TrajStatus::propagated";

    std::cerr << "  Final molecule COM = (" << moleculeList[0].comCoord.x << ", "
              << moleculeList[0].comCoord.y << ", " << moleculeList[0].comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: spherical boundary, target complex has no cross partners.
//         The sphere branch of the dispatcher must also propagate the complex
//         and leave it inside the spherical volume.
// -----------------------------------------------------------------------------
void test_sscr_sphere_no_partners_translates_complex()
{
    std::cerr << "\n[TEST] test_sscr_sphere_no_partners_translates_complex\n"
              << "  Source file: sweep_separation_complex_rot.cpp\n"
              << "  Function:    sweep_separation_complex_rot (sphere branch)\n"
              << "  Scenario:    isSphere == true, one complex, no cross partners.\n"
              << "  Criteria:    the complex is displaced, remains inside the\n"
              << "               sphere, and its trajectory vectors are zeroed.\n";

    sscr_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscr_make_sphere_membrane();
    std::vector<MolTemplate> molTemplateList { sscr_make_moltemplate() };
    std::vector<ForwardRxn> forwardRxns { sscr_make_rxn(1.0) };

    // Start near the centre of the sphere so no reflection is required.
    std::vector<Molecule> moleculeList { sscr_make_molecule(0, 0, Coord { 1.0, 2.0, -1.0 }) };
    std::vector<Complex> complexList { sscr_make_complex(0, std::vector<int> { 0 }, Coord { 1.0, 2.0, -1.0 }, 10.0) };

    const Vector step { 2.0, 0.0, 0.0 };
    complexList[0].trajTrans = step;

    std::cerr << "  Sphere radius = " << membraneObject.sphereR
              << ", proposed step = (" << step.x << ", " << step.y << ", " << step.z << ")\n";
    std::cerr << "  Calling sweep_separation_complex_rot...\n";
    sweep_separation_complex_rot(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Expected position: the molecule started at (1,2,-1) and should be shifted.
    EXPECT_NEAR(moleculeList[0].comCoord.x, 1.0 + step.x, 1e-10)
        << "molecule COM x should be displaced by trajTrans.x on the sphere branch";
    EXPECT_NEAR(moleculeList[0].comCoord.y, 2.0 + step.y, 1e-10)
        << "molecule COM y should be displaced by trajTrans.y on the sphere branch";
    EXPECT_NEAR(moleculeList[0].comCoord.z, -1.0 + step.z, 1e-10)
        << "molecule COM z should be displaced by trajTrans.z on the sphere branch";

    // Sanity: still inside the spherical volume.
    const double radial = moleculeList[0].comCoord.get_magnitude();
    std::cerr << "  Final radial distance = " << radial << " (radius " << membraneObject.sphereR << ")\n";
    EXPECT_LE(radial, membraneObject.sphereR + 1e-6)
        << "the propagated complex must stay inside the spherical boundary";

    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "trajTrans.y must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "trajTrans.z must be zeroed after propagation";
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "member molecule should be flagged TrajStatus::propagated";
}

// -----------------------------------------------------------------------------
// Test 3: two molecules that are cross partners but belong to the SAME complex.
//         They cannot diffuse relative to one another, so the sweep must skip
//         them entirely: the originally proposed step survives untouched even
//         though the two interfaces are much closer than the binding radius.
// -----------------------------------------------------------------------------
void test_sscr_same_complex_partners_are_ignored()
{
    std::cerr << "\n[TEST] test_sscr_same_complex_partners_are_ignored\n"
              << "  Source file: sweep_separation_complex_rot.cpp\n"
              << "  Function:    sweep_separation_complex_rot (box branch)\n"
              << "  Scenario:    two cross partners inside one complex, separated\n"
              << "               by 0.4 nm while the binding radius is 1.0 nm.\n"
              << "  Criteria:    no re-sampling occurs (partners in the same complex\n"
              << "               cannot move relative to one another), so both\n"
              << "               molecules move by exactly the proposed trajTrans.\n";

    sscr_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscr_make_box_membrane();
    std::vector<MolTemplate> molTemplateList { sscr_make_moltemplate() };
    const double bindRadius = 1.0;
    std::vector<ForwardRxn> forwardRxns { sscr_make_rxn(bindRadius) };

    // Two molecules, 0.4 nm apart, both members of complex 0.
    const Coord pos0 { 0.0, 0.0, 0.0 };
    const Coord pos1 { 0.4, 0.0, 0.0 };
    std::vector<Molecule> moleculeList {
        sscr_make_molecule(0, 0, pos0),
        sscr_make_molecule(1, 0, pos1),
    };
    std::vector<Complex> complexList {
        sscr_make_complex(0, std::vector<int> { 0, 1 }, Coord { 0.2, 0.0, 0.0 }, 10.0)
    };
    complexList[0].numEachMol = std::vector<int> { 2 };
    complexList[0].mass = 2.0;

    // Register each molecule as the other's "cross" partner via reaction 0.
    moleculeList[0].crossbase.push_back(1);
    moleculeList[0].mycrossint.push_back(0);
    moleculeList[0].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    moleculeList[0].probvec.push_back(1.0);

    moleculeList[1].crossbase.push_back(0);
    moleculeList[1].mycrossint.push_back(0);
    moleculeList[1].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    moleculeList[1].probvec.push_back(1.0);

    const Vector step { 1.0, 1.0, 0.0 };
    complexList[0].trajTrans = step;

    std::cerr << "  Interface separation before the call = "
              << sscr_distance(moleculeList[0].interfaceList[0].coord, moleculeList[1].interfaceList[0].coord)
              << " nm (bindRadius = " << bindRadius << " nm)\n";
    std::cerr << "  Calling sweep_separation_complex_rot...\n";
    sweep_separation_complex_rot(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Both members must have moved by exactly the originally proposed step:
    // any re-sampling would have produced a random (different) displacement.
    EXPECT_NEAR(moleculeList[0].comCoord.x, pos0.x + step.x, 1e-10)
        << "molecule 0 x should be moved by the *original* trajTrans (no re-sampling)";
    EXPECT_NEAR(moleculeList[0].comCoord.y, pos0.y + step.y, 1e-10)
        << "molecule 0 y should be moved by the *original* trajTrans (no re-sampling)";
    EXPECT_NEAR(moleculeList[1].comCoord.x, pos1.x + step.x, 1e-10)
        << "molecule 1 x should be moved by the *original* trajTrans (no re-sampling)";
    EXPECT_NEAR(moleculeList[1].comCoord.y, pos1.y + step.y, 1e-10)
        << "molecule 1 y should be moved by the *original* trajTrans (no re-sampling)";

    // Rigid-body motion: the internal separation is unchanged.
    const double sepAfter
        = sscr_distance(moleculeList[0].interfaceList[0].coord, moleculeList[1].interfaceList[0].coord);
    std::cerr << "  Interface separation after the call = " << sepAfter << " nm\n";
    EXPECT_NEAR(sepAfter, 0.4, 1e-10)
        << "members of the same complex must keep their relative geometry";

    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed after propagation";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "trajTrans.y must be zeroed after propagation";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "every member of the propagated complex should be flagged propagated";
}

// -----------------------------------------------------------------------------
// Test 4: two overlapping molecules in DIFFERENT complexes.
//         The sweep must re-sample trajectories until the two reacting
//         interfaces are no longer inside the binding radius.
// -----------------------------------------------------------------------------
void test_sscr_resolves_overlap_between_complexes()
{
    std::cerr << "\n[TEST] test_sscr_resolves_overlap_between_complexes\n"
              << "  Source file: sweep_separation_complex_rot.cpp\n"
              << "  Function:    sweep_separation_complex_rot (box branch)\n"
              << "  Scenario:    two separate complexes whose reacting interfaces\n"
              << "               are 0.5 nm apart while the binding radius is 1.0 nm\n"
              << "               and the proposed steps are zero (i.e. overlapping).\n"
              << "  Criteria:    trajectories are re-sampled so the projected\n"
              << "               interface separation is >= the binding radius.\n";

    sscr_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscr_make_box_membrane();
    std::vector<MolTemplate> molTemplateList { sscr_make_moltemplate() };
    const double bindRadius = 1.0;
    std::vector<ForwardRxn> forwardRxns { sscr_make_rxn(bindRadius) };

    // Molecule 0 in complex 0, molecule 1 in complex 1, 0.5 nm apart.
    const Coord pos0 { 0.0, 0.0, 0.0 };
    const Coord pos1 { 0.5, 0.0, 0.0 };
    std::vector<Molecule> moleculeList {
        sscr_make_molecule(0, 0, pos0),
        sscr_make_molecule(1, 1, pos1),
    };
    std::vector<Complex> complexList {
        sscr_make_complex(0, std::vector<int> { 0 }, pos0, 10.0),
        sscr_make_complex(1, std::vector<int> { 1 }, pos1, 10.0),
    };

    // Cross lists: molecule 0 sees molecule 1 through reaction 0 (and back).
    moleculeList[0].crossbase.push_back(1);
    moleculeList[0].mycrossint.push_back(0);
    moleculeList[0].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    moleculeList[0].probvec.push_back(1.0);

    moleculeList[1].crossbase.push_back(0);
    moleculeList[1].mycrossint.push_back(0);
    moleculeList[1].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    moleculeList[1].probvec.push_back(1.0);

    // Both complexes start with a zero step, so as-is they overlap.
    complexList[0].trajTrans = Vector { 0.0, 0.0, 0.0 };
    complexList[1].trajTrans = Vector { 0.0, 0.0, 0.0 };

    const Coord mol1PosBefore = moleculeList[1].comCoord;

    std::cerr << "  Initial interface separation = "
              << sscr_distance(moleculeList[0].interfaceList[0].coord, moleculeList[1].interfaceList[0].coord)
              << " nm (bindRadius = " << bindRadius << " nm) -> overlap expected\n";
    std::cerr << "  Calling sweep_separation_complex_rot...\n";
    sweep_separation_complex_rot(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Complex 0 has been propagated; complex 1 has only had its trajectory
    // re-sampled (it will be applied when molecule 1 is swept in turn).
    const Coord projectedIface1 { moleculeList[1].interfaceList[0].coord.x + complexList[1].trajTrans.x,
        moleculeList[1].interfaceList[0].coord.y + complexList[1].trajTrans.y,
        moleculeList[1].interfaceList[0].coord.z + complexList[1].trajTrans.z };
    const double finalSep = sscr_distance(moleculeList[0].interfaceList[0].coord, projectedIface1);

    std::cerr << "  Final (projected) interface separation = " << finalSep << " nm\n";
    EXPECT_GE(finalSep, bindRadius - 1e-9)
        << "after sweeping, the reacting interfaces must not be closer than the binding radius";

    // The target complex must have actually moved away from its start.
    EXPECT_GT(sscr_distance(moleculeList[0].comCoord, pos0), 0.0)
        << "the target complex should have been given a new (re-sampled) displacement";

    // Only the target complex is propagated during this call.
    EXPECT_NEAR(moleculeList[1].comCoord.x, mol1PosBefore.x, 1e-10)
        << "the partner complex must not be propagated by this call";
    EXPECT_NEAR(moleculeList[1].comCoord.y, mol1PosBefore.y, 1e-10)
        << "the partner complex must not be propagated by this call";
    EXPECT_NEAR(moleculeList[1].comCoord.z, mol1PosBefore.z, 1e-10)
        << "the partner complex must not be propagated by this call";

    // The target complex's stored trajectory is always cleared at the end.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "target trajTrans.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "target trajTrans.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "target trajTrans.z must be zeroed";
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "the swept molecule should be flagged TrajStatus::propagated";
}

// -----------------------------------------------------------------------------
// Test 5: two cross partners in different complexes that are already far
//         enough apart.  No overlap is detected, therefore the originally
//         proposed step must survive.
// -----------------------------------------------------------------------------
void test_sscr_no_overlap_keeps_proposed_step()
{
    std::cerr << "\n[TEST] test_sscr_no_overlap_keeps_proposed_step\n"
              << "  Source file: sweep_separation_complex_rot.cpp\n"
              << "  Function:    sweep_separation_complex_rot (box branch)\n"
              << "  Scenario:    two complexes 20 nm apart with a 1 nm binding\n"
              << "               radius and small proposed steps.\n"
              << "  Criteria:    no re-sampling -> the target moves by exactly the\n"
              << "               proposed trajTrans and the partner's trajTrans is\n"
              << "               left untouched.\n";

    sscr_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscr_make_box_membrane();
    std::vector<MolTemplate> molTemplateList { sscr_make_moltemplate() };
    const double bindRadius = 1.0;
    std::vector<ForwardRxn> forwardRxns { sscr_make_rxn(bindRadius) };

    const Coord pos0 { 0.0, 0.0, 0.0 };
    const Coord pos1 { 20.0, 0.0, 0.0 };
    std::vector<Molecule> moleculeList {
        sscr_make_molecule(0, 0, pos0),
        sscr_make_molecule(1, 1, pos1),
    };
    std::vector<Complex> complexList {
        sscr_make_complex(0, std::vector<int> { 0 }, pos0, 10.0),
        sscr_make_complex(1, std::vector<int> { 1 }, pos1, 10.0),
    };

    moleculeList[0].crossbase.push_back(1);
    moleculeList[0].mycrossint.push_back(0);
    moleculeList[0].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    moleculeList[0].probvec.push_back(1.0);

    moleculeList[1].crossbase.push_back(0);
    moleculeList[1].mycrossint.push_back(0);
    moleculeList[1].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    moleculeList[1].probvec.push_back(1.0);

    // Small steps that do not bring the two interfaces within the binding radius.
    const Vector step0 { -1.0, 0.5, 0.0 };
    const Vector step1 { 1.0, 0.0, 0.0 };
    complexList[0].trajTrans = step0;
    complexList[1].trajTrans = step1;

    std::cerr << "  Initial separation = " << sscr_distance(pos0, pos1) << " nm\n";
    std::cerr << "  Calling sweep_separation_complex_rot...\n";
    sweep_separation_complex_rot(0, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // Target complex moved by exactly the proposed step (no re-sampling).
    EXPECT_NEAR(moleculeList[0].comCoord.x, pos0.x + step0.x, 1e-10)
        << "no overlap -> target x must move by the originally proposed step";
    EXPECT_NEAR(moleculeList[0].comCoord.y, pos0.y + step0.y, 1e-10)
        << "no overlap -> target y must move by the originally proposed step";
    EXPECT_NEAR(moleculeList[0].comCoord.z, pos0.z + step0.z, 1e-10)
        << "no overlap -> target z must move by the originally proposed step";

    // Partner trajectory untouched.
    EXPECT_NEAR(complexList[1].trajTrans.x, step1.x, 1e-10)
        << "the partner complex's trajTrans.x should not be re-sampled";
    EXPECT_NEAR(complexList[1].trajTrans.y, step1.y, 1e-10)
        << "the partner complex's trajTrans.y should not be re-sampled";
    EXPECT_NEAR(complexList[1].trajTrans.z, step1.z, 1e-10)
        << "the partner complex's trajTrans.z should not be re-sampled";

    // Partner has not been propagated.
    EXPECT_NEAR(moleculeList[1].comCoord.x, pos1.x, 1e-10)
        << "the partner complex must not be propagated by this call";

    std::cerr << "  Final target COM = (" << moleculeList[0].comCoord.x << ", "
              << moleculeList[0].comCoord.y << ", " << moleculeList[0].comCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each scenario lives in its own TEST so that a failure
// in one does not prevent the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(SweepSeparationComplexRot, BoxNoPartnersTranslatesComplex) { test_sscr_box_no_partners_translates_complex(); }
TEST(SweepSeparationComplexRot, SphereNoPartnersTranslatesComplex) { test_sscr_sphere_no_partners_translates_complex(); }
TEST(SweepSeparationComplexRot, SameComplexPartnersAreIgnored) { test_sscr_same_complex_partners_are_ignored(); }
TEST(SweepSeparationComplexRot, ResolvesOverlapBetweenComplexes) { test_sscr_resolves_overlap_between_complexes(); }
TEST(SweepSeparationComplexRot, NoOverlapKeepsProposedStep) { test_sscr_no_overlap_keeps_proposed_step(); }