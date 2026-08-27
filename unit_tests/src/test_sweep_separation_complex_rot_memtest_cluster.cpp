/*! \file test_sweep_separation_complex_rot_memtest_cluster.cpp
 *
 * ### Unit tests for src/trajectory_functions/sweep_separation_complex_rot_memtest_cluster.cpp
 *
 * The function under test is a thin dispatcher:
 *
 *     void sweep_separation_complex_rot_memtest_cluster(int simItr, int pro1Index,
 *              Parameters& params, std::vector<Molecule>& moleculeList,
 *              std::vector<Complex>& complexList,
 *              const std::vector<ForwardRxn>& forwardRxns,
 *              const std::vector<MolTemplate>& molTemplateList,
 *              const Membrane& membraneObject)
 *
 * It forwards to
 *   - sweep_separation_complex_rot_memtest_cluster_sphere() when membraneObject.isSphere
 *   - sweep_separation_complex_rot_memtest_cluster_box()    otherwise
 *
 * The sweep routines resample the propagation vectors of a *cluster* of complexes
 * until no pair of reacting interfaces is closer than its binding radius, and then
 * physically propagate the complexes.  Proteins that live inside the *same* complex
 * are never resampled against each other (they cannot diffuse relative to one another).
 *
 * The tests below therefore check invariants that the routine must preserve no matter
 * which internal branch is taken:
 *
 *   1. The "position after this time step" of two paired interfaces, i.e.
 *      (interface coordinate + complex trajTrans), is never closer than the binding
 *      radius once the sweep has returned.  This quantity is correct whether or not
 *      the complex has already been physically propagated, because propagation both
 *      updates the coordinate and zeroes trajTrans.
 *   2. The motion is a rigid-body motion: the |interface - centre-of-mass| distance of
 *      each molecule and the distance between two molecules of the same complex are
 *      preserved.
 *   3. Nothing is pushed outside of the simulation volume (box or sphere), and no
 *      coordinate becomes NaN/inf.
 *
 * To make invariant (2) exact, all diffusion constants for rotation (Dr) are set to
 * zero so every resampled rotation angle is identically zero and propagation reduces
 * to a pure translation.
 */

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

// The random number generator is a global owned by unit_tests/src/gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Small, self-contained builders for the objects the sweep routine touches.
// Every field the sweep (and the Complex::propagate it calls) reads is filled in,
// because these structs are indexed without bounds checks deep inside NERDSS.
// -----------------------------------------------------------------------------

/*! \brief Deterministically (re)initialise the shared GSL RNG so resampling is
 *         reproducible from run to run. */
void sscrmc_init_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

/*! \brief One molecule type: a single interface sitting 1 nm along +x from the COM.
 *
 * \param[in] onMembrane when true the type is a lipid with D.z == 0 (2D diffusion).
 *            Rotational diffusion is always zero so propagation is a pure translation.
 */
MolTemplate sscrmc_make_moltemplate(bool onMembrane)
{
    MolTemplate molTemplate;
    molTemplate.molTypeIndex = 0;
    molTemplate.molName = "sweepTestMol";
    molTemplate.mass = 1.0;
    molTemplate.radius = 1.0;
    molTemplate.copies = 2;
    molTemplate.isLipid = onMembrane;
    molTemplate.isImplicitLipid = false;
    molTemplate.isRod = false;
    molTemplate.isPoint = false;
    molTemplate.isPromoter = false;
    molTemplate.checkOverlap = false;

    // Translational diffusion: 1 nm^2/us laterally, no z motion when on the membrane.
    molTemplate.D = Coord(1.0, 1.0, onMembrane ? 0.0 : 1.0);
    // No rotational diffusion -> resampled rotation angles are exactly zero.
    molTemplate.Dr = Coord(0.0, 0.0, 0.0);

    Interface iface("site", Coord(1.0, 0.0, 0.0));
    iface.index = 0;
    iface.stateList.emplace_back("site", 0);
    molTemplate.interfaceList.push_back(iface);

    return molTemplate;
}

/*! \brief Build a molecule of the single template type at a given centre of mass. */
Molecule sscrmc_make_molecule(int index, int comIndex, const Coord& com, bool isLipid)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = isLipid;
    mol.isImplicitLipid = false;
    mol.isPromoter = false;
    mol.linksToSurface = 0;
    mol.mySubVolIndex = 0;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    // Single free interface, 1 nm along +x from the COM.
    Molecule::Iface iface;
    iface.coord = Coord(com.x + 1.0, com.y, com.z);
    iface.relIndex = 0;
    iface.index = 0;
    iface.stateIndex = 0;
    iface.molTypeIndex = 0;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);
    mol.freelist.push_back(0);

    return mol;
}

/*! \brief Build a complex that owns the supplied member molecules. */
Complex sscrmc_make_complex(int index, const std::vector<int>& members,
    const std::vector<Molecule>& moleculeList, const MolTemplate& molTemplate, bool onSurface)
{
    Complex com;
    com.index = index;
    com.id = index;
    com.memberList = members;
    com.numEachMol = std::vector<int>(1, static_cast<int>(members.size()));
    com.lastNumberUpdateItrEachMol = std::vector<long long int>(1, 0);

    // Mass weighted centre of mass of the members (the invariant Complex keeps).
    double totMass{ 0.0 };
    Coord accum{ 0.0, 0.0, 0.0 };
    for (auto memMol : members) {
        totMass += moleculeList[memMol].mass;
        accum.x += moleculeList[memMol].comCoord.x * moleculeList[memMol].mass;
        accum.y += moleculeList[memMol].comCoord.y * moleculeList[memMol].mass;
        accum.z += moleculeList[memMol].comCoord.z * moleculeList[memMol].mass;
    }
    com.mass = totMass;
    com.comCoord = Coord(accum.x / totMass, accum.y / totMass, accum.z / totMass);
    com.tmpComCoord = com.comCoord;

    com.radius = molTemplate.radius + 1.0;
    com.D = molTemplate.D;
    com.Dr = molTemplate.Dr; // zero -> pure translation
    com.OnSurface = onSurface;
    com.onFiber = false;
    com.isEmpty = false;
    com.linksToSurface = 0;
    com.iLipidIndex = 0;
    com.ncross = 0;
    com.trajStatus = TrajStatus::none;
    com.trajTrans.x = 0.0;
    com.trajTrans.y = 0.0;
    com.trajTrans.z = 0.0;
    com.trajRot.zero_crds();

    return com;
}

/*! \brief Register molecule \p b as a possible reaction partner of molecule \p a. */
void sscrmc_link_cross(Molecule& a, const Molecule& b, int rxnIndex)
{
    a.crossbase.push_back(b.index);   // partner molecule index
    a.mycrossint.push_back(0);        // my relative interface index
    a.crossrxn.push_back(std::array<int, 3>{ rxnIndex, 0, 0 });
    a.probvec.push_back(0.0);
}

/*! \brief A minimal bimolecular ForwardRxn: interface 0 binds interface 0. */
ForwardRxn sscrmc_make_rxn(double bindRadius)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.absRxnIndex = 0;
    rxn.relRxnIndex = 0;
    rxn.isReversible = false;
    rxn.isOnMem = false;
    rxn.bindRadius = bindRadius;
    rxn.bindRadius2D = bindRadius;
    rxn.reactantListNew.emplace_back("site", 0, 0, 0, '\0', false);
    rxn.reactantListNew.emplace_back("site", 0, 0, 0, '\0', false);
    rxn.productListNew = rxn.reactantListNew;
    rxn.rateList.emplace_back();
    rxn.rateList.back().rate = 1.0;
    return rxn;
}

/*! \brief A 100 nm cubic reflecting water box.
 *
 * RS3Dvect is sized to 500 entries because the sweep routines scan
 * RS3Dvect[400..499] looking for the reflecting-surface value of a molecule type;
 * filling it with -1 guarantees no match and hence RS3Dinput == 0. */
Membrane sscrmc_make_box_membrane()
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.implicitLipid = false;
    membraneObject.TwoD = false;
    membraneObject.hasCompartment = false;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double>{ 100.0, 100.0, 100.0 });
    membraneObject.xBCtype = "reflect";
    membraneObject.yBCtype = "reflect";
    membraneObject.zBCtype = "reflect";
    membraneObject.RS3Dvect.assign(500, -1.0);
    return membraneObject;
}

/*! \brief A spherical boundary of radius 100 nm (with the enclosing box set as well). */
Membrane sscrmc_make_sphere_membrane()
{
    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.isBox = false;
    membraneObject.implicitLipid = false;
    membraneObject.hasCompartment = false;
    membraneObject.sphereR = 100.0;
    membraneObject.sphereVol = (4.0 / 3.0) * M_PI * std::pow(membraneObject.sphereR, 3.0);
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double>{ 200.0, 200.0, 200.0 });
    membraneObject.xBCtype = "reflect";
    membraneObject.yBCtype = "reflect";
    membraneObject.zBCtype = "reflect";
    membraneObject.RS3Dvect.assign(500, -1.0);
    return membraneObject;
}

/*! \brief Simulation parameters, 1 us time step. */
Parameters sscrmc_make_params()
{
    Parameters params;
    params.rank = 0;
    params.timeStep = 1.0;
    Parameters::dt = 1.0;
    params.nItr = 10;
    params.numMolTypes = 1;
    params.numTotalSpecies = 1;
    params.numTotalComplex = 2;
    params.overlapSepLimit = 0.1;
    params.scaleMaxDisplace = 100.0;
    params.rMaxLimit = 10.0;
    params.rMaxRadius = 1.0;
    params.clusterOverlapCheck = true;
    return params;
}

/*! \brief true if every component of the coordinate is a finite number. */
bool sscrmc_is_finite(const Coord& c)
{
    return std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z);
}

/*! \brief Distance between molecule \p molIndex's interface 0 and its own COM. */
double sscrmc_iface_offset(const std::vector<Molecule>& moleculeList, int molIndex)
{
    const Coord& com = moleculeList[molIndex].comCoord;
    const Coord& ifc = moleculeList[molIndex].interfaceList[0].coord;
    const double dx{ ifc.x - com.x };
    const double dy{ ifc.y - com.y };
    const double dz{ ifc.z - com.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief "Position after this time step" of a molecule's interface 0.
 *
 * If the sweep already propagated the parent complex, trajTrans is zero and the
 * coordinate is the new one.  If it has not, the coordinate is the old one and
 * trajTrans still holds the (resampled) displacement.  Either way the sum is the
 * position the algorithm guaranteed to be overlap free. */
Coord sscrmc_effective_iface(const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList, int molIndex)
{
    const Complex& com = complexList[moleculeList[molIndex].myComIndex];
    const Coord& ifc = moleculeList[molIndex].interfaceList[0].coord;
    return Coord(ifc.x + com.trajTrans.x, ifc.y + com.trajTrans.y, ifc.z + com.trajTrans.z);
}

/*! \brief Separation of the two paired interfaces at the end of the time step. */
double sscrmc_effective_separation(const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList, int mol1, int mol2)
{
    const Coord a{ sscrmc_effective_iface(moleculeList, complexList, mol1) };
    const Coord b{ sscrmc_effective_iface(moleculeList, complexList, mol2) };
    const double dx{ a.x - b.x };
    const double dy{ a.y - b.y };
    const double dz{ a.z - b.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief Distance between the centres of mass of two molecules. */
double sscrmc_mol_mol_distance(const std::vector<Molecule>& moleculeList, int mol1, int mol2)
{
    const double dx{ moleculeList[mol1].comCoord.x - moleculeList[mol2].comCoord.x };
    const double dy{ moleculeList[mol1].comCoord.y - moleculeList[mol2].comCoord.y };
    const double dz{ moleculeList[mol1].comCoord.z - moleculeList[mol2].comCoord.z };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: box geometry, two membrane complexes that are nowhere near each other.
//         The sweep must find no overlap and must not corrupt anything.
// -----------------------------------------------------------------------------
void test_sscrmc_box_no_overlap()
{
    std::cerr << "\n[TEST] test_sscrmc_box_no_overlap\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster (box branch)\n"
              << "  Scenario    : two single-molecule complexes on the membrane, 18 nm apart,\n"
              << "                each carrying a small translational trajectory.\n"
              << "  Pass criteria: coordinates stay finite and inside the 100 nm box, the\n"
              << "                molecules stay rigid, and the end-of-step interface\n"
              << "                separation is >= the binding radius (1 nm).\n";

    sscrmc_init_rng();

    Parameters params{ sscrmc_make_params() };
    Membrane membraneObject{ sscrmc_make_box_membrane() };

    std::vector<MolTemplate> molTemplateList{ sscrmc_make_moltemplate(true) };
    MolTemplate::numMolTypes = 1;

    const double bindRadius{ 1.0 };
    std::vector<ForwardRxn> forwardRxns{ sscrmc_make_rxn(bindRadius) };

    // Two lipid molecules on the membrane plane (z = -49 nm), well separated in x.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmc_make_molecule(0, 0, Coord(-10.0, 0.0, -49.0), true));
    moleculeList.push_back(sscrmc_make_molecule(1, 1, Coord(10.0, 0.0, -49.0), true));

    // Each molecule sees the other as a possible binding partner.
    sscrmc_link_cross(moleculeList[0], moleculeList[1], 0);
    sscrmc_link_cross(moleculeList[1], moleculeList[0], 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmc_make_complex(0, { 0 }, moleculeList, molTemplateList[0], true));
    complexList.push_back(sscrmc_make_complex(1, { 1 }, moleculeList, molTemplateList[0], true));
    complexList[0].ncross = 1;
    complexList[1].ncross = 1;

    // Small in-plane displacements that cannot bring the interfaces together.
    complexList[0].trajTrans.x = 0.50;
    complexList[0].trajTrans.y = 0.25;
    complexList[1].trajTrans.x = -0.50;
    complexList[1].trajTrans.y = 0.10;

    const double offset0Before{ sscrmc_iface_offset(moleculeList, 0) };
    const double offset1Before{ sscrmc_iface_offset(moleculeList, 1) };
    std::cerr << "  Initial |iface - COM| = " << offset0Before << " and " << offset1Before << " nm\n";

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster(simItr=1, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_cluster(
        1, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // --- finiteness -----------------------------------------------------------
    for (int m = 0; m < 2; ++m) {
        EXPECT_TRUE(sscrmc_is_finite(moleculeList[m].comCoord))
            << "molecule " << m << " COM must remain a finite coordinate";
        EXPECT_TRUE(sscrmc_is_finite(moleculeList[m].interfaceList[0].coord))
            << "molecule " << m << " interface must remain a finite coordinate";
        EXPECT_TRUE(sscrmc_is_finite(complexList[m].comCoord))
            << "complex " << m << " COM must remain a finite coordinate";
    }

    // --- rigid body -----------------------------------------------------------
    EXPECT_NEAR(sscrmc_iface_offset(moleculeList, 0), offset0Before, 1e-9)
        << "molecule 0 must move rigidly (Dr = 0 -> pure translation)";
    EXPECT_NEAR(sscrmc_iface_offset(moleculeList, 1), offset1Before, 1e-9)
        << "molecule 1 must move rigidly (Dr = 0 -> pure translation)";

    // --- complex COM tracks its (single) member -------------------------------
    for (int c = 0; c < 2; ++c) {
        EXPECT_NEAR(complexList[c].comCoord.x, moleculeList[c].comCoord.x, 1e-9)
            << "complex " << c << " COM.x must equal its single member's COM.x";
        EXPECT_NEAR(complexList[c].comCoord.y, moleculeList[c].comCoord.y, 1e-9)
            << "complex " << c << " COM.y must equal its single member's COM.y";
        EXPECT_NEAR(complexList[c].comCoord.z, moleculeList[c].comCoord.z, 1e-9)
            << "complex " << c << " COM.z must equal its single member's COM.z";
    }

    // --- overlap free ---------------------------------------------------------
    const double sep{ sscrmc_effective_separation(moleculeList, complexList, 0, 1) };
    std::cerr << "  End-of-step interface separation = " << sep
              << " nm (binding radius " << bindRadius << " nm)\n";
    EXPECT_GE(sep, bindRadius - 1e-6)
        << "well separated complexes must remain further apart than the binding radius";

    // --- still inside the simulation volume -----------------------------------
    for (int m = 0; m < 2; ++m) {
        const Coord eff{ sscrmc_effective_iface(moleculeList, complexList, m) };
        EXPECT_LE(std::abs(eff.x), membraneObject.waterBox.x / 2.0 + 1e-6)
            << "molecule " << m << " must stay inside the box in x";
        EXPECT_LE(std::abs(eff.y), membraneObject.waterBox.y / 2.0 + 1e-6)
            << "molecule " << m << " must stay inside the box in y";
        EXPECT_LE(std::abs(eff.z), membraneObject.waterBox.z / 2.0 + 1e-6)
            << "molecule " << m << " must stay inside the box in z";
    }
}

// -----------------------------------------------------------------------------
// Test 2: box geometry, two complexes whose reacting interfaces start *inside*
//         the binding radius.  The sweep must resample until they are separated.
// -----------------------------------------------------------------------------
void test_sscrmc_box_resolves_overlap()
{
    std::cerr << "\n[TEST] test_sscrmc_box_resolves_overlap\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster (box branch)\n"
              << "  Scenario    : the two paired interfaces start only 0.5 nm apart while the\n"
              << "                binding radius is 1.0 nm, so the initial (zero) trajectories\n"
              << "                overlap and must be resampled.\n"
              << "  Pass criteria: after the sweep the end-of-step interface separation is\n"
              << "                >= the binding radius, and at least one complex moved.\n";

    sscrmc_init_rng();

    Parameters params{ sscrmc_make_params() };
    Membrane membraneObject{ sscrmc_make_box_membrane() };

    std::vector<MolTemplate> molTemplateList{ sscrmc_make_moltemplate(true) };
    MolTemplate::numMolTypes = 1;

    const double bindRadius{ 1.0 };
    std::vector<ForwardRxn> forwardRxns{ sscrmc_make_rxn(bindRadius) };

    // Interfaces sit at x = -9.0 and x = -8.5 -> only 0.5 nm apart (overlapping).
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmc_make_molecule(0, 0, Coord(-10.0, 0.0, -49.0), true));
    moleculeList.push_back(sscrmc_make_molecule(1, 1, Coord(-9.5, 0.0, -49.0), true));

    sscrmc_link_cross(moleculeList[0], moleculeList[1], 0);
    sscrmc_link_cross(moleculeList[1], moleculeList[0], 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmc_make_complex(0, { 0 }, moleculeList, molTemplateList[0], true));
    complexList.push_back(sscrmc_make_complex(1, { 1 }, moleculeList, molTemplateList[0], true));
    complexList[0].ncross = 1;
    complexList[1].ncross = 1;
    // Both complexes are free to be resampled (trajStatus == none on their molecules).

    const double sepBefore{ sscrmc_effective_separation(moleculeList, complexList, 0, 1) };
    std::cerr << "  Initial end-of-step interface separation = " << sepBefore
              << " nm (binding radius " << bindRadius << " nm) -> overlapping\n";
    EXPECT_LT(sepBefore, bindRadius) << "test setup must actually start in an overlapping state";

    const Coord eff0Before{ sscrmc_effective_iface(moleculeList, complexList, 0) };
    const Coord eff1Before{ sscrmc_effective_iface(moleculeList, complexList, 1) };
    const double offset0Before{ sscrmc_iface_offset(moleculeList, 0) };

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster(simItr=1, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_cluster(
        1, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    const double sepAfter{ sscrmc_effective_separation(moleculeList, complexList, 0, 1) };
    std::cerr << "  Resolved end-of-step interface separation = " << sepAfter << " nm\n";

    // The whole purpose of the routine: no residual overlap of the reacting pair.
    EXPECT_GE(sepAfter, bindRadius - 1e-6)
        << "the sweep must resample trajectories until the pair is no longer overlapping";

    // Something has to have changed - either a resampled trajectory or a propagation.
    const Coord eff0After{ sscrmc_effective_iface(moleculeList, complexList, 0) };
    const Coord eff1After{ sscrmc_effective_iface(moleculeList, complexList, 1) };
    const double moved0{ std::sqrt(std::pow(eff0After.x - eff0Before.x, 2)
        + std::pow(eff0After.y - eff0Before.y, 2) + std::pow(eff0After.z - eff0Before.z, 2)) };
    const double moved1{ std::sqrt(std::pow(eff1After.x - eff1Before.x, 2)
        + std::pow(eff1After.y - eff1Before.y, 2) + std::pow(eff1After.z - eff1Before.z, 2)) };
    std::cerr << "  End-of-step displacement caused by resampling: mol0 = " << moved0
              << " nm, mol1 = " << moved1 << " nm\n";
    EXPECT_GT(moved0 + moved1, 1e-9)
        << "at least one of the overlapping complexes must have been resampled";

    // Rigid body motion is still preserved while resampling.
    EXPECT_NEAR(sscrmc_iface_offset(moleculeList, 0), offset0Before, 1e-9)
        << "resampling must not deform the molecule (Dr = 0 -> translation only)";

    // Still inside the box.
    for (int m = 0; m < 2; ++m) {
        const Coord eff{ sscrmc_effective_iface(moleculeList, complexList, m) };
        EXPECT_TRUE(sscrmc_is_finite(eff)) << "molecule " << m << " coordinate must stay finite";
        EXPECT_LE(std::abs(eff.x), membraneObject.waterBox.x / 2.0 + 1e-6)
            << "molecule " << m << " must stay inside the box in x";
        EXPECT_LE(std::abs(eff.y), membraneObject.waterBox.y / 2.0 + 1e-6)
            << "molecule " << m << " must stay inside the box in y";
        EXPECT_LE(std::abs(eff.z), membraneObject.waterBox.z / 2.0 + 1e-6)
            << "molecule " << m << " must stay inside the box in z";
    }
}

// -----------------------------------------------------------------------------
// Test 3: two cross-listed molecules that belong to the SAME complex.
//         Such pairs are deliberately ignored (they cannot diffuse relative to
//         one another) so the internal geometry must be untouched.
// -----------------------------------------------------------------------------
void test_sscrmc_same_complex_pairs_ignored()
{
    std::cerr << "\n[TEST] test_sscrmc_same_complex_pairs_ignored\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster (box branch)\n"
              << "  Scenario    : one complex holding two molecules that list each other as\n"
              << "                cross partners, closer than the binding radius.\n"
              << "  Pass criteria: the routine returns (no attempt to fix an unfixable overlap)\n"
              << "                and the intra-complex geometry is exactly preserved.\n";

    sscrmc_init_rng();

    Parameters params{ sscrmc_make_params() };
    Membrane membraneObject{ sscrmc_make_box_membrane() };

    std::vector<MolTemplate> molTemplateList{ sscrmc_make_moltemplate(true) };
    MolTemplate::numMolTypes = 1;

    // Binding radius larger than the intra-complex separation, so a naive
    // implementation would try (and fail) to resolve the "overlap".
    const double bindRadius{ 3.0 };
    std::vector<ForwardRxn> forwardRxns{ sscrmc_make_rxn(bindRadius) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmc_make_molecule(0, 0, Coord(-1.0, 0.0, -49.0), true));
    moleculeList.push_back(sscrmc_make_molecule(1, 0, Coord(1.0, 0.0, -49.0), true));

    sscrmc_link_cross(moleculeList[0], moleculeList[1], 0);
    sscrmc_link_cross(moleculeList[1], moleculeList[0], 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmc_make_complex(0, { 0, 1 }, moleculeList, molTemplateList[0], true));
    complexList[0].numEachMol = std::vector<int>(1, 2);
    complexList[0].ncross = 2;
    complexList[0].trajTrans.x = 0.5;
    complexList[0].trajTrans.y = 0.5;

    const double molMolBefore{ sscrmc_mol_mol_distance(moleculeList, 0, 1) };
    const double offset0Before{ sscrmc_iface_offset(moleculeList, 0) };
    const double offset1Before{ sscrmc_iface_offset(moleculeList, 1) };
    std::cerr << "  Intra-complex molecule-molecule distance before = " << molMolBefore << " nm\n";

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster(simItr=1, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_cluster(
        1, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    const double molMolAfter{ sscrmc_mol_mol_distance(moleculeList, 0, 1) };
    std::cerr << "  Intra-complex molecule-molecule distance after  = " << molMolAfter << " nm\n";

    EXPECT_NEAR(molMolAfter, molMolBefore, 1e-9)
        << "members of the same complex must not move relative to one another";
    EXPECT_NEAR(sscrmc_iface_offset(moleculeList, 0), offset0Before, 1e-9)
        << "molecule 0 must stay rigid";
    EXPECT_NEAR(sscrmc_iface_offset(moleculeList, 1), offset1Before, 1e-9)
        << "molecule 1 must stay rigid";

    EXPECT_TRUE(sscrmc_is_finite(moleculeList[0].comCoord)) << "molecule 0 COM must stay finite";
    EXPECT_TRUE(sscrmc_is_finite(moleculeList[1].comCoord)) << "molecule 1 COM must stay finite";
    EXPECT_TRUE(sscrmc_is_finite(complexList[0].comCoord)) << "complex COM must stay finite";

    // The complex COM must remain the (mass weighted) mean of its two members.
    const double meanX{ 0.5 * (moleculeList[0].comCoord.x + moleculeList[1].comCoord.x) };
    const double meanY{ 0.5 * (moleculeList[0].comCoord.y + moleculeList[1].comCoord.y) };
    EXPECT_NEAR(complexList[0].comCoord.x, meanX, 1e-9)
        << "complex COM.x must be the mass weighted mean of its members";
    EXPECT_NEAR(complexList[0].comCoord.y, meanY, 1e-9)
        << "complex COM.y must be the mass weighted mean of its members";
}

// -----------------------------------------------------------------------------
// Test 4: a complex with no cross partners at all (degenerate cluster).
// -----------------------------------------------------------------------------
void test_sscrmc_no_cross_partners()
{
    std::cerr << "\n[TEST] test_sscrmc_no_cross_partners\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster (box branch)\n"
              << "  Scenario    : the target molecule has an empty crossbase, so the cluster\n"
              << "                pair list is empty and there is nothing to resolve.\n"
              << "  Pass criteria: the call returns cleanly, coordinates stay finite, rigid\n"
              << "                and inside the box.\n";

    sscrmc_init_rng();

    Parameters params{ sscrmc_make_params() };
    Membrane membraneObject{ sscrmc_make_box_membrane() };

    std::vector<MolTemplate> molTemplateList{ sscrmc_make_moltemplate(true) };
    MolTemplate::numMolTypes = 1;

    std::vector<ForwardRxn> forwardRxns{ sscrmc_make_rxn(1.0) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmc_make_molecule(0, 0, Coord(0.0, 0.0, -49.0), true));
    // Intentionally no crossbase entries.

    std::vector<Complex> complexList;
    complexList.push_back(sscrmc_make_complex(0, { 0 }, moleculeList, molTemplateList[0], true));
    complexList[0].ncross = 0;
    complexList[0].trajTrans.x = 1.0;
    complexList[0].trajTrans.y = -1.0;

    const double offsetBefore{ sscrmc_iface_offset(moleculeList, 0) };

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster(simItr=1, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_cluster(
        1, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    EXPECT_TRUE(sscrmc_is_finite(moleculeList[0].comCoord)) << "COM must stay finite";
    EXPECT_TRUE(sscrmc_is_finite(moleculeList[0].interfaceList[0].coord))
        << "interface coordinate must stay finite";
    EXPECT_NEAR(sscrmc_iface_offset(moleculeList, 0), offsetBefore, 1e-9)
        << "the lone molecule must remain rigid";

    const Coord eff{ sscrmc_effective_iface(moleculeList, complexList, 0) };
    std::cerr << "  End-of-step interface position = " << eff << '\n';
    EXPECT_LE(std::abs(eff.x), membraneObject.waterBox.x / 2.0 + 1e-6)
        << "the molecule must stay inside the box in x";
    EXPECT_LE(std::abs(eff.y), membraneObject.waterBox.y / 2.0 + 1e-6)
        << "the molecule must stay inside the box in y";
    EXPECT_LE(std::abs(eff.z), membraneObject.waterBox.z / 2.0 + 1e-6)
        << "the molecule must stay inside the box in z";
}

// -----------------------------------------------------------------------------
// Test 5: sphere geometry - exercises the other branch of the dispatcher.
// -----------------------------------------------------------------------------
void test_sscrmc_sphere_branch()
{
    std::cerr << "\n[TEST] test_sscrmc_sphere_branch\n"
              << "  Source file : sweep_separation_complex_rot_memtest_cluster.cpp\n"
              << "  Function    : sweep_separation_complex_rot_memtest_cluster (sphere branch)\n"
              << "  Scenario    : membraneObject.isSphere == true (R = 100 nm) with two well\n"
              << "                separated complexes diffusing in solution inside the sphere.\n"
              << "  Pass criteria: the sphere routine runs, keeps coordinates finite and inside\n"
              << "                the sphere, keeps the molecules rigid, and leaves the pair\n"
              << "                further apart than the binding radius.\n";

    sscrmc_init_rng();

    Parameters params{ sscrmc_make_params() };
    Membrane membraneObject{ sscrmc_make_sphere_membrane() };

    // 3D diffusion (not membrane bound) so the sphere-surface propagation math is
    // not exercised with a degenerate (zero length) displacement.
    std::vector<MolTemplate> molTemplateList{ sscrmc_make_moltemplate(false) };
    MolTemplate::numMolTypes = 1;

    const double bindRadius{ 1.0 };
    std::vector<ForwardRxn> forwardRxns{ sscrmc_make_rxn(bindRadius) };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(sscrmc_make_molecule(0, 0, Coord(-10.0, 0.0, 0.0), false));
    moleculeList.push_back(sscrmc_make_molecule(1, 1, Coord(10.0, 0.0, 0.0), false));

    sscrmc_link_cross(moleculeList[0], moleculeList[1], 0);
    sscrmc_link_cross(moleculeList[1], moleculeList[0], 0);

    std::vector<Complex> complexList;
    complexList.push_back(sscrmc_make_complex(0, { 0 }, moleculeList, molTemplateList[0], false));
    complexList.push_back(sscrmc_make_complex(1, { 1 }, moleculeList, molTemplateList[0], false));
    complexList[0].ncross = 1;
    complexList[1].ncross = 1;

    // Non-degenerate displacements (needed by the spherical geometry helpers).
    complexList[0].trajTrans.x = 0.20;
    complexList[0].trajTrans.y = 0.30;
    complexList[0].trajTrans.z = 0.10;
    complexList[1].trajTrans.x = -0.20;
    complexList[1].trajTrans.y = 0.10;
    complexList[1].trajTrans.z = -0.30;

    const double offset0Before{ sscrmc_iface_offset(moleculeList, 0) };
    const double offset1Before{ sscrmc_iface_offset(moleculeList, 1) };

    std::cerr << "  Calling sweep_separation_complex_rot_memtest_cluster(simItr=1, pro1Index=0)...\n";
    sweep_separation_complex_rot_memtest_cluster(
        1, 0, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    for (int m = 0; m < 2; ++m) {
        EXPECT_TRUE(sscrmc_is_finite(moleculeList[m].comCoord))
            << "molecule " << m << " COM must remain finite in the spherical branch";
        EXPECT_TRUE(sscrmc_is_finite(moleculeList[m].interfaceList[0].coord))
            << "molecule " << m << " interface must remain finite in the spherical branch";
    }

    EXPECT_NEAR(sscrmc_iface_offset(moleculeList, 0), offset0Before, 1e-9)
        << "molecule 0 must remain rigid (Dr = 0)";
    EXPECT_NEAR(sscrmc_iface_offset(moleculeList, 1), offset1Before, 1e-9)
        << "molecule 1 must remain rigid (Dr = 0)";

    // Both molecules must remain inside the spherical boundary.
    for (int m = 0; m < 2; ++m) {
        const Coord eff{ sscrmc_effective_iface(moleculeList, complexList, m) };
        const double radial{ std::sqrt(eff.x * eff.x + eff.y * eff.y + eff.z * eff.z) };
        std::cerr << "  molecule " << m << " end-of-step radial distance = " << radial
                  << " nm (sphere radius " << membraneObject.sphereR << " nm)\n";
        EXPECT_LE(radial, membraneObject.sphereR + 1e-6)
            << "molecule " << m << " must remain inside the spherical boundary";
    }

    const double sep{ sscrmc_effective_separation(moleculeList, complexList, 0, 1) };
    std::cerr << "  End-of-step interface separation = " << sep << " nm\n";
    EXPECT_GE(sep, bindRadius - 1e-6)
        << "the sphere branch must also leave the pair further apart than the binding radius";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario is run inside its own TEST so that a
// failure in one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(SweepSeparationComplexRotMemtestCluster, BoxNoOverlap) { test_sscrmc_box_no_overlap(); }
TEST(SweepSeparationComplexRotMemtestCluster, BoxResolvesOverlap) { test_sscrmc_box_resolves_overlap(); }
TEST(SweepSeparationComplexRotMemtestCluster, SameComplexPairsIgnored) { test_sscrmc_same_complex_pairs_ignored(); }
TEST(SweepSeparationComplexRotMemtestCluster, NoCrossPartners) { test_sscrmc_no_cross_partners(); }
TEST(SweepSeparationComplexRotMemtestCluster, SphereBranch) { test_sscrmc_sphere_branch(); }