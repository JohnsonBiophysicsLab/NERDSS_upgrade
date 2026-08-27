/*! \file test_sweep_separation_complex_rot_sphere.cpp
 *
 * ### Unit test for src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp
 *
 * The single function under test is
 *
 *     void sweep_separation_complex_rot_sphere(int simItr, int pro1Index,
 *                                              Parameters& params,
 *                                              std::vector<Molecule>& moleculeList,
 *                                              std::vector<Complex>& complexList,
 *                                              const std::vector<ForwardRxn>& forwardRxns,
 *                                              const std::vector<MolTemplate>& molTemplateList,
 *                                              const Membrane& membraneObject);
 *
 * Behaviour that the tests below pin down (taken directly from the implementation):
 *
 *   1. For every member molecule of the complex that owns `pro1Index`, each entry of
 *      `crossbase` is examined.  Partners that are implicit lipids, or that belong to
 *      the *same* complex, are skipped entirely.
 *   2. If the trial step would bring a reacting interface pair closer than the
 *      reaction binding radius, both complexes get a freshly sampled trajectory
 *      (Gaussian, scaled by sqrt(2*dt*D)), and the partner is only resampled when its
 *      molecule's TrajStatus is `none` or `canBeResampled`.
 *   3. Whatever happens in the sweep loop, the *first* complex (and only that one) is
 *      finally propagated, and its trajTrans/trajRot are zeroed afterwards.
 *
 * Notes on set-up requirements discovered by reading the implementation:
 *   - `membraneObject.RS3Dvect` is indexed unconditionally at [300..499], so the
 *     vector must contain at least 500 entries or the call is undefined behaviour.
 *   - `Complex::propagate()` calls `update_properties()`, which divides by the total
 *     mass of the member molecules; every molecule therefore needs a non-zero mass and
 *     a valid `molTypeIndex` into molTemplateList.
 *   - With `trajRot == (0,0,0)` propagate takes the pure-translation fast path, so the
 *     resulting coordinates are exact and can be compared with tight tolerances.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

#include "classes/class_Molecule_Complex.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

// -----------------------------------------------------------------------------
// Small helpers used to build a fully-initialized, minimal simulation system.
// Everything here is prefixed with "sscrs_" so it cannot collide with the rest
// of the generated test suite.
// -----------------------------------------------------------------------------
namespace {

//! Tolerance used for the deterministic (no-resampling) comparisons.
constexpr double kSscrsTol = 1e-9;

/*! \brief Make sure the global GSL generator used by GaussV()/rand_gsl() exists.
 *
 * `r` is defined (as nullptr) in unit_tests/src/gtest_main.cpp, so we allocate it
 * once here if nobody else did, and re-seed it every call so that the random
 * branches of the sweep are reproducible from run to run.
 */
void sscrs_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, 42);
}

/*! \brief Build a one-interface molecule owned by complex `comIndex`. */
Molecule sscrs_make_molecule(int index, int comIndex, const Coord& com, const Coord& ifaceCoord)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.molTypeIndex = 0;          // single molecule type in these tests
    mol.myComIndex = comIndex;
    mol.mass = 1.0;                // must be non-zero: update_properties divides by it
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = ifaceCoord;
    iface.index = 0;
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // No encounters recorded yet; individual tests fill these in.
    mol.crossbase.clear();
    mol.mycrossint.clear();
    mol.crossrxn.clear();

    return mol;
}

/*! \brief Register molecule `partnerIndex` as a possible reaction partner of `mol`. */
void sscrs_add_cross(Molecule& mol, int partnerIndex, int myIfaceRelIndex, int rxnIndex)
{
    mol.crossbase.push_back(partnerIndex);
    mol.mycrossint.push_back(myIfaceRelIndex);
    mol.crossrxn.push_back(std::array<int, 3> { rxnIndex, 0, 0 });
}

/*! \brief Build a complex holding the given member molecules, with isotropic D/Dr. */
Complex sscrs_make_complex(int index, const Coord& com, const std::vector<int>& members,
    double transD, double rotD)
{
    Complex targCom;
    targCom.index = index;
    targCom.id = index;
    targCom.comCoord = com;
    targCom.memberList = members;
    targCom.numEachMol = std::vector<int> { static_cast<int>(members.size()) };
    targCom.lastNumberUpdateItrEachMol.resize(1);
    targCom.mass = static_cast<double>(members.size());
    targCom.radius = 1.0;
    targCom.D.x = transD;
    targCom.D.y = transD;
    targCom.D.z = transD;
    targCom.Dr.x = rotD;
    targCom.Dr.y = rotD;
    targCom.Dr.z = rotD;
    targCom.isEmpty = false;
    targCom.OnSurface = false;   // these complexes live in the interior of the sphere
    targCom.onFiber = false;
    targCom.linksToSurface = 0;
    targCom.trajStatus = TrajStatus::none;
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot.x = 0.0;
    targCom.trajRot.y = 0.0;
    targCom.trajRot.z = 0.0;
    return targCom;
}

/*! \brief One molecule template; its D/Dr are what update_properties() re-derives. */
MolTemplate sscrs_make_template(double transD, double rotD)
{
    MolTemplate molTemplate;
    molTemplate.molTypeIndex = 0;
    molTemplate.molName = "A";
    molTemplate.mass = 1.0;
    molTemplate.radius = 1.0;
    molTemplate.copies = 1;
    molTemplate.isLipid = false;
    molTemplate.isImplicitLipid = false;
    molTemplate.isRod = false;
    molTemplate.isPoint = false;
    molTemplate.isPromoter = false;
    molTemplate.D.x = transD;
    molTemplate.D.y = transD;
    molTemplate.D.z = transD;
    molTemplate.Dr.x = rotD;
    molTemplate.Dr.y = rotD;
    molTemplate.Dr.z = rotD;

    Interface iface { "a", Coord { 0.0, 0.0, 0.0 } };
    iface.index = 0;
    iface.stateList.emplace_back('\0', 0);
    molTemplate.interfaceList.push_back(iface);

    return molTemplate;
}

/*! \brief A bimolecular reaction between interface 0 of the (single) molecule type. */
ForwardRxn sscrs_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.bindRadius = bindRadius;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    // Both reactants use relative interface index 0, so the partner interface
    // resolved inside the sweep is also 0.
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("a", 0, 0, 0, '\0', false);
    return rxn;
}

/*! \brief A spherical membrane of radius 100 with a legal-sized RS3D lookup table. */
Membrane sscrs_make_membrane()
{
    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.isBox = false;
    membraneObject.sphereR = 100.0;
    membraneObject.implicitLipid = false;
    membraneObject.hasCompartment = false;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
    // The sweep reads RS3Dvect[300..499] unconditionally -> must be >= 500 long.
    membraneObject.RS3Dvect.assign(500, 0.0);
    return membraneObject;
}

/*! \brief Euclidean distance between two coordinates (test-side convenience). */
double sscrs_distance(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a complex with no recorded encounters is simply propagated.
// -----------------------------------------------------------------------------
void test_sscrs_no_partners_pure_propagation()
{
    std::cerr << "\n[TEST] test_sscrs_no_partners_pure_propagation\n"
              << "  Source file: src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_sphere\n"
              << "  Scenario:    lone complex, crossbase empty, trajRot == 0.\n"
              << "  Pass:        molecule + interface are translated by trajTrans exactly,\n"
              << "               the complex COM follows, and trajTrans/trajRot are zeroed.\n";

    sscrs_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrs_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrs_make_template(1.0, 0.0) };
    std::vector<ForwardRxn> forwardRxns { sscrs_make_rxn(1.0) };

    // A single molecule sitting at (10,0,0) with its interface one nm further out.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(
        sscrs_make_molecule(0, 0, Coord { 10.0, 0.0, 0.0 }, Coord { 11.0, 0.0, 0.0 }));

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrs_make_complex(0, Coord { 10.0, 0.0, 0.0 }, std::vector<int> { 0 }, 1.0, 0.0));

    // Trial displacement that must be applied verbatim (no partners -> no overlap).
    complexList[0].trajTrans.x = 2.0;
    complexList[0].trajTrans.y = -1.0;
    complexList[0].trajTrans.z = 0.5;

    std::cerr << "  Calling sweep_separation_complex_rot_sphere(simItr=1, pro1Index=0)...\n";
    sweep_separation_complex_rot_sphere(1, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    std::cerr << "  Molecule COM after sweep = (" << moleculeList[0].comCoord.x << ", "
              << moleculeList[0].comCoord.y << ", " << moleculeList[0].comCoord.z << ")\n";

    // Molecule and interface must have moved by exactly the requested displacement.
    EXPECT_NEAR(moleculeList[0].comCoord.x, 12.0, kSscrsTol) << "COM x should be 10 + 2";
    EXPECT_NEAR(moleculeList[0].comCoord.y, -1.0, kSscrsTol) << "COM y should be 0 - 1";
    EXPECT_NEAR(moleculeList[0].comCoord.z, 0.5, kSscrsTol) << "COM z should be 0 + 0.5";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.x, 13.0, kSscrsTol)
        << "Interface x should be 11 + 2";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.y, -1.0, kSscrsTol)
        << "Interface y should be 0 - 1";
    EXPECT_NEAR(moleculeList[0].interfaceList[0].coord.z, 0.5, kSscrsTol)
        << "Interface z should be 0 + 0.5";

    // update_properties() recomputes the complex COM from its (single) member.
    EXPECT_NEAR(complexList[0].comCoord.x, 12.0, kSscrsTol) << "Complex COM x follows the member";
    EXPECT_NEAR(complexList[0].comCoord.y, -1.0, kSscrsTol) << "Complex COM y follows the member";
    EXPECT_NEAR(complexList[0].comCoord.z, 0.5, kSscrsTol) << "Complex COM z follows the member";

    // Trajectories are cleared at the very end of the routine.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "trajTrans.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "trajTrans.z must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.x, 0.0) << "trajRot.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.y, 0.0) << "trajRot.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajRot.z, 0.0) << "trajRot.z must be zeroed";

    // Both the molecule and the complex are flagged as propagated.
    EXPECT_EQ(static_cast<int>(moleculeList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "Molecule trajStatus should be 'propagated' after the sweep";
    EXPECT_EQ(static_cast<int>(complexList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "Complex trajStatus should be 'propagated' after the sweep";
}

// -----------------------------------------------------------------------------
// Test 2: partners inside the *same* complex are never checked for overlap.
// -----------------------------------------------------------------------------
void test_sscrs_same_complex_partner_is_skipped()
{
    std::cerr << "\n[TEST] test_sscrs_same_complex_partner_is_skipped\n"
              << "  Source file: src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_sphere\n"
              << "  Scenario:    two molecules of the SAME complex list each other in\n"
              << "               crossbase and sit only 0.5 nm apart while bindRadius = 5.\n"
              << "  Pass:        no resampling happens (members cannot diffuse relative to\n"
              << "               each other), so the displacement is applied exactly.\n";

    sscrs_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrs_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrs_make_template(1.0, 0.0) };
    std::vector<ForwardRxn> forwardRxns { sscrs_make_rxn(5.0) }; // deliberately large

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(
        sscrs_make_molecule(0, 0, Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 0.0 }));
    moleculeList.push_back(
        sscrs_make_molecule(1, 0, Coord { 10.5, 0.0, 0.0 }, Coord { 10.5, 0.0, 0.0 }));
    // They "see" each other, but they belong to the same complex.
    sscrs_add_cross(moleculeList[0], 1, 0, 0);
    sscrs_add_cross(moleculeList[1], 0, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrs_make_complex(0, Coord { 10.25, 0.0, 0.0 }, std::vector<int> { 0, 1 }, 1.0, 0.0));
    complexList[0].trajTrans.x = 1.0;

    std::cerr << "  Calling sweep_separation_complex_rot_sphere(simItr=2, pro1Index=0)...\n";
    sweep_separation_complex_rot_sphere(2, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    std::cerr << "  Member COM x after sweep = " << moleculeList[0].comCoord.x << " and "
              << moleculeList[1].comCoord.x << " (complex COM x = " << complexList[0].comCoord.x
              << ")\n";

    // A deterministic +1 nm shift proves that the trajectory was never resampled.
    EXPECT_NEAR(moleculeList[0].comCoord.x, 11.0, kSscrsTol)
        << "Member 0 should translate by exactly +1 nm";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 11.5, kSscrsTol)
        << "Member 1 should translate by exactly +1 nm";
    EXPECT_NEAR(complexList[0].comCoord.x, 11.25, kSscrsTol)
        << "Complex COM is the mass-weighted mean of both members";

    // Rigid-body translation: the internal separation must be preserved.
    EXPECT_NEAR(sscrs_distance(moleculeList[0].comCoord, moleculeList[1].comCoord), 0.5, kSscrsTol)
        << "Intra-complex separation must be unchanged by a pure translation";
}

// -----------------------------------------------------------------------------
// Test 3: implicit-lipid partners are skipped before any distance evaluation.
// -----------------------------------------------------------------------------
void test_sscrs_implicit_lipid_partner_is_skipped()
{
    std::cerr << "\n[TEST] test_sscrs_implicit_lipid_partner_is_skipped\n"
              << "  Source file: src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_sphere\n"
              << "  Scenario:    the only crossbase entry is an implicit lipid that sits\n"
              << "               right on top of the moving interface.\n"
              << "  Pass:        the `continue` for implicit lipids means no overlap is\n"
              << "               registered, so the trial move survives untouched.\n";

    sscrs_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrs_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrs_make_template(1.0, 0.0) };
    std::vector<ForwardRxn> forwardRxns { sscrs_make_rxn(5.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(
        sscrs_make_molecule(0, 0, Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 0.0 }));
    // Partner is an implicit lipid at (almost) the same place -> would overlap badly.
    moleculeList.push_back(
        sscrs_make_molecule(1, 1, Coord { 10.1, 0.0, 0.0 }, Coord { 10.1, 0.0, 0.0 }));
    moleculeList[1].isImplicitLipid = true;
    sscrs_add_cross(moleculeList[0], 1, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrs_make_complex(0, Coord { 10.0, 0.0, 0.0 }, std::vector<int> { 0 }, 1.0, 0.0));
    complexList.push_back(
        sscrs_make_complex(1, Coord { 10.1, 0.0, 0.0 }, std::vector<int> { 1 }, 0.0, 0.0));

    complexList[0].trajTrans.y = 3.0;

    std::cerr << "  Calling sweep_separation_complex_rot_sphere(simItr=3, pro1Index=0)...\n";
    sweep_separation_complex_rot_sphere(3, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    std::cerr << "  Molecule 0 COM after sweep = (" << moleculeList[0].comCoord.x << ", "
              << moleculeList[0].comCoord.y << ", " << moleculeList[0].comCoord.z << ")\n";

    EXPECT_NEAR(moleculeList[0].comCoord.x, 10.0, kSscrsTol) << "x should be untouched";
    EXPECT_NEAR(moleculeList[0].comCoord.y, 3.0, kSscrsTol)
        << "y should be shifted by exactly the requested +3 nm (no resampling)";
    EXPECT_NEAR(moleculeList[0].comCoord.z, 0.0, kSscrsTol) << "z should be untouched";

    // The implicit lipid is never moved by this routine.
    EXPECT_NEAR(moleculeList[1].comCoord.x, 10.1, kSscrsTol) << "Implicit lipid must not move";
    EXPECT_NEAR(moleculeList[1].comCoord.y, 0.0, kSscrsTol) << "Implicit lipid must not move";
}

// -----------------------------------------------------------------------------
// Test 4: a distant partner in another complex leaves the trajectory alone.
// -----------------------------------------------------------------------------
void test_sscrs_distant_partner_no_resample()
{
    std::cerr << "\n[TEST] test_sscrs_distant_partner_no_resample\n"
              << "  Source file: src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_sphere\n"
              << "  Scenario:    partner belongs to a different complex but is ~40 nm away\n"
              << "               while bindRadius = 1 nm.\n"
              << "  Pass:        separation^2 >= bindRadius^2 so nothing is resampled; the\n"
              << "               partner complex keeps its own stored trajectory.\n";

    sscrs_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrs_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrs_make_template(1.0, 0.0) };
    std::vector<ForwardRxn> forwardRxns { sscrs_make_rxn(1.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(
        sscrs_make_molecule(0, 0, Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 0.0 }));
    moleculeList.push_back(
        sscrs_make_molecule(1, 1, Coord { 50.0, 0.0, 0.0 }, Coord { 50.0, 0.0, 0.0 }));
    sscrs_add_cross(moleculeList[0], 1, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrs_make_complex(0, Coord { 10.0, 0.0, 0.0 }, std::vector<int> { 0 }, 1.0, 0.0));
    complexList.push_back(
        sscrs_make_complex(1, Coord { 50.0, 0.0, 0.0 }, std::vector<int> { 1 }, 1.0, 0.0));

    complexList[0].trajTrans.x = 1.5;
    // Sentinel trajectory on the partner: it must survive verbatim.
    complexList[1].trajTrans.z = -0.25;

    std::cerr << "  Calling sweep_separation_complex_rot_sphere(simItr=4, pro1Index=0)...\n";
    sweep_separation_complex_rot_sphere(4, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    std::cerr << "  Molecule 0 COM x = " << moleculeList[0].comCoord.x
              << ", partner trajTrans.z = " << complexList[1].trajTrans.z << '\n';

    EXPECT_NEAR(moleculeList[0].comCoord.x, 11.5, kSscrsTol)
        << "No overlap -> the exact trial displacement is applied";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.z, -0.25)
        << "Partner complex trajectory must not be resampled when there is no overlap";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 50.0, kSscrsTol)
        << "Only the first complex is propagated by this routine";
    EXPECT_EQ(static_cast<int>(moleculeList[1].trajStatus), static_cast<int>(TrajStatus::none))
        << "The partner molecule is not propagated, so its trajStatus stays 'none'";
}

// -----------------------------------------------------------------------------
// Test 5: an overlapping partner in another complex triggers resampling of both.
// -----------------------------------------------------------------------------
void test_sscrs_overlap_resamples_both_complexes()
{
    std::cerr << "\n[TEST] test_sscrs_overlap_resamples_both_complexes\n"
              << "  Source file: src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_sphere\n"
              << "  Scenario:    two single-molecule complexes 0.5 nm apart with\n"
              << "               bindRadius = 2 nm -> guaranteed overlap on the first pass.\n"
              << "  Pass:        the moving complex ends up somewhere other than its\n"
              << "               original trial destination (its step was resampled), the\n"
              << "               partner complex receives a fresh non-zero trajectory, and\n"
              << "               the partner molecule itself is left where it was.\n";

    sscrs_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrs_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrs_make_template(1.0, 0.0) };
    std::vector<ForwardRxn> forwardRxns { sscrs_make_rxn(2.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(
        sscrs_make_molecule(0, 0, Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 0.0 }));
    moleculeList.push_back(
        sscrs_make_molecule(1, 1, Coord { 10.5, 0.0, 0.0 }, Coord { 10.5, 0.0, 0.0 }));
    sscrs_add_cross(moleculeList[0], 1, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrs_make_complex(0, Coord { 10.0, 0.0, 0.0 }, std::vector<int> { 0 }, 1.0, 0.0));
    complexList.push_back(
        sscrs_make_complex(1, Coord { 10.5, 0.0, 0.0 }, std::vector<int> { 1 }, 1.0, 0.0));

    // Zero trial displacements: the initial 0.5 nm separation is inside bindRadius.
    const Coord originalPartnerPos = moleculeList[1].comCoord;

    std::cerr << "  Calling sweep_separation_complex_rot_sphere(simItr=5, pro1Index=0)...\n";
    sweep_separation_complex_rot_sphere(5, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    const double movedDist = sscrs_distance(moleculeList[0].comCoord, Coord { 10.0, 0.0, 0.0 });
    const double partnerTrajMag = std::sqrt(complexList[1].trajTrans.x * complexList[1].trajTrans.x
        + complexList[1].trajTrans.y * complexList[1].trajTrans.y
        + complexList[1].trajTrans.z * complexList[1].trajTrans.z);

    std::cerr << "  Moving molecule displaced by " << movedDist << " nm; partner |trajTrans| = "
              << partnerTrajMag << '\n';

    // The trial step was (0,0,0); any non-zero displacement proves it was resampled.
    EXPECT_GT(movedDist, 0.0)
        << "Overlap must force a fresh Gaussian trajectory for the moving complex";

    // The partner molecule had TrajStatus::none, so its complex is resampled too,
    // but the partner is NOT propagated inside this routine.
    EXPECT_GT(partnerTrajMag, 0.0)
        << "Partner complex should receive a resampled (non-zero) trajTrans";
    EXPECT_NEAR(moleculeList[1].comCoord.x, originalPartnerPos.x, kSscrsTol)
        << "Partner molecule coordinates must not change";
    EXPECT_NEAR(moleculeList[1].comCoord.y, originalPartnerPos.y, kSscrsTol)
        << "Partner molecule coordinates must not change";
    EXPECT_NEAR(moleculeList[1].comCoord.z, originalPartnerPos.z, kSscrsTol)
        << "Partner molecule coordinates must not change";

    // The moving complex is always propagated and cleaned up at the end.
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.x, 0.0) << "trajTrans.x must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.y, 0.0) << "trajTrans.y must be zeroed";
    EXPECT_DOUBLE_EQ(complexList[0].trajTrans.z, 0.0) << "trajTrans.z must be zeroed";
    EXPECT_EQ(static_cast<int>(complexList[0].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "The swept complex is propagated regardless of the overlap outcome";

    // The complex COM must still track its single member molecule.
    EXPECT_NEAR(complexList[0].comCoord.x, moleculeList[0].comCoord.x, 1e-8)
        << "update_properties() keeps the COM consistent with the member";
    EXPECT_NEAR(complexList[0].comCoord.y, moleculeList[0].comCoord.y, 1e-8)
        << "update_properties() keeps the COM consistent with the member";
    EXPECT_NEAR(complexList[0].comCoord.z, moleculeList[0].comCoord.z, 1e-8)
        << "update_properties() keeps the COM consistent with the member";

    // Sanity: the resampled step is a diffusive move, so it should stay small
    // compared with the sphere radius, and the complex must remain inside.
    const double radial = moleculeList[0].comCoord.get_magnitude();
    EXPECT_LE(radial, membraneObject.sphereR)
        << "The resampled/reflected position must remain inside the spherical boundary";
}

// -----------------------------------------------------------------------------
// Test 6: an already-propagated partner is never resampled.
// -----------------------------------------------------------------------------
void test_sscrs_propagated_partner_not_resampled()
{
    std::cerr << "\n[TEST] test_sscrs_propagated_partner_not_resampled\n"
              << "  Source file: src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp\n"
              << "  Function:    sweep_separation_complex_rot_sphere\n"
              << "  Scenario:    overlapping partner whose molecule already carries\n"
              << "               TrajStatus::propagated (it moved earlier this step).\n"
              << "  Pass:        only the swept complex is resampled; the partner keeps the\n"
              << "               exact sentinel trajectory it was given.\n";

    sscrs_init_rng();

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = sscrs_make_membrane();
    std::vector<MolTemplate> molTemplateList { sscrs_make_template(1.0, 0.0) };
    std::vector<ForwardRxn> forwardRxns { sscrs_make_rxn(2.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(
        sscrs_make_molecule(0, 0, Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 0.0 }));
    moleculeList.push_back(
        sscrs_make_molecule(1, 1, Coord { 10.5, 0.0, 0.0 }, Coord { 10.5, 0.0, 0.0 }));
    // The partner has already been moved this timestep and may not be resampled.
    moleculeList[1].trajStatus = TrajStatus::propagated;
    sscrs_add_cross(moleculeList[0], 1, 0, 0);

    std::vector<Complex> complexList;
    complexList.push_back(
        sscrs_make_complex(0, Coord { 10.0, 0.0, 0.0 }, std::vector<int> { 0 }, 1.0, 0.0));
    complexList.push_back(
        sscrs_make_complex(1, Coord { 10.5, 0.0, 0.0 }, std::vector<int> { 1 }, 1.0, 0.0));

    // Tiny sentinel step so the pair still overlaps on the first evaluation.
    const double sentinel = 0.01;
    complexList[1].trajTrans.x = sentinel;

    std::cerr << "  Calling sweep_separation_complex_rot_sphere(simItr=6, pro1Index=0)...\n";
    sweep_separation_complex_rot_sphere(6, 0, params, moleculeList, complexList, forwardRxns,
        molTemplateList, membraneObject);

    const double movedDist = sscrs_distance(moleculeList[0].comCoord, Coord { 10.0, 0.0, 0.0 });
    std::cerr << "  Moving molecule displaced by " << movedDist
              << " nm; partner trajTrans = (" << complexList[1].trajTrans.x << ", "
              << complexList[1].trajTrans.y << ", " << complexList[1].trajTrans.z << ")\n";

    // The swept complex started with a zero trial step, so it must have been resampled.
    EXPECT_GT(movedDist, 0.0) << "The swept complex must still be resampled on overlap";

    // The already-propagated partner keeps its exact trajectory.
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.x, sentinel)
        << "A propagated partner must keep its trajTrans.x untouched";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.y, 0.0)
        << "A propagated partner must keep its trajTrans.y untouched";
    EXPECT_DOUBLE_EQ(complexList[1].trajTrans.z, 0.0)
        << "A propagated partner must keep its trajTrans.z untouched";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 10.5, kSscrsTol)
        << "The partner molecule is not propagated here";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario runs independently so that a failure in
// one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(SweepSeparationComplexRotSphere, NoPartnersPurePropagation)
{
    test_sscrs_no_partners_pure_propagation();
}
TEST(SweepSeparationComplexRotSphere, SameComplexPartnerIsSkipped)
{
    test_sscrs_same_complex_partner_is_skipped();
}
TEST(SweepSeparationComplexRotSphere, ImplicitLipidPartnerIsSkipped)
{
    test_sscrs_implicit_lipid_partner_is_skipped();
}
TEST(SweepSeparationComplexRotSphere, DistantPartnerNoResample)
{
    test_sscrs_distant_partner_no_resample();
}
TEST(SweepSeparationComplexRotSphere, OverlapResamplesBothComplexes)
{
    test_sscrs_overlap_resamples_both_complexes();
}
TEST(SweepSeparationComplexRotSphere, PropagatedPartnerNotResampled)
{
    test_sscrs_propagated_partner_not_resampled();
}