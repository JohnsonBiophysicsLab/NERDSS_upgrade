/*! \file test_check_overlap.cpp
 *
 * ### Unit tests for src/reactions/check_overlap.cpp
 *
 * The file under test contains exactly one function:
 *
 * \code
 * void check_overlap(std::vector<int>& region, unsigned simItr,
 *                    Parameters& params, std::vector<Molecule>& moleculeList,
 *                    std::vector<Complex>& complexList, SimulVolume& simulVolume,
 *                    std::vector<ForwardRxn>& forwardRxns,
 *                    std::vector<BackRxn>& backRxns,
 *                    std::vector<CreateDestructRxn>& createDestructRxns,
 *                    std::vector<MolTemplate>& molTemplateList,
 *                    std::map<std::string,int>& observablesList,
 *                    copyCounters& counterArrays, Membrane& membraneObject,
 *                    MpiContext& mpiContext);
 * \endcode
 *
 * `check_overlap()` walks a list of molecule indices (`region`) and, for every
 * molecule that is still eligible to move, finalises its position:
 *
 *   1. molecules that are empty / implicit lipids / already propagated /
 *      ghosted are skipped outright,
 *   2. molecules whose parent complex registered crossings (`ncross > 0`) are
 *      handed to the `sweep_separation_complex_rot*` overlap sweepers,
 *   3. molecules whose parent complex has `ncross == 0` are simply propagated
 *      (after first generating a propagation vector if the trajectory has not
 *      been sampled yet), with a spherical-boundary reflection when the system
 *      uses a sphere.
 *
 * These tests exercise branch (1) and branch (3) - the branches whose observable
 * effect is a deterministic change (or non-change) of the molecule/complex
 * coordinates and TrajStatus flags.
 *
 * Branch (2) (`ncross > 0`) delegates to `sweep_separation_complex_rot*`, which
 * requires a *fully* populated crossing bookkeeping (crossbase / mycrossint /
 * crossrxn / probvec on every member molecule, matching reaction entries, ...).
 * Feeding those routines a partially initialised system would read past the end
 * of those vectors and take the whole gtest binary down with it, so instead we
 * only verify here that the guard ordering in `check_overlap()` prevents the
 * sweep from being reached for molecules that were already propagated.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "classes/class_SimulVolume.hpp"            // SimulVolume (+ MpiContext via split.cpp)
#include "classes/class_copyCounters.hpp"           // copyCounters
#include "math/rand_gsl.hpp"                        // extern gsl_rng* r
#include "reactions/shared_reaction_functions.hpp"  // declaration of check_overlap

// The RNG instance and the match counter are defined in unit_tests/src/gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// RNG bootstrap.
//
// `check_overlap()` can call create_complex_propagation_vectors(), which draws
// Gaussian random numbers through the global GSL generator, so the generator has
// to exist before we call in. We re-seed deterministically on every use so the
// tests are reproducible no matter what ran before them.
// -----------------------------------------------------------------------------
void co_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, 42);
}

// -----------------------------------------------------------------------------
// A stand-in MpiContext.
//
// Inside check_overlap() `mpiContext` is only ever touched from `if (DEBUG)`
// blocks, and macro.hpp defines `DEBUG` to `false`, so none of those statements
// execute. Handing over a zero-initialised buffer of the right size/alignment
// therefore gives the function a perfectly valid reference while avoiding any
// dependency on the (MPI-build specific) MpiContext constructor signature.
// -----------------------------------------------------------------------------
MpiContext& co_dummy_mpi_context()
{
    static std::aligned_storage<sizeof(MpiContext), alignof(MpiContext)>::type storage {};
    return *reinterpret_cast<MpiContext*>(&storage);
}

// -----------------------------------------------------------------------------
// Builders for the minimal - but fully initialised - system objects.
// -----------------------------------------------------------------------------

/*! \brief Build a one-interface MolTemplate.
 *
 * `Complex::update_properties()` (called from `Complex::propagate()`) reads the
 * template's radius and diffusion constants, so both must be meaningful.
 * A diffusion constant of exactly zero makes the propagation vector generation
 * deterministic (sqrt(2*dt*0)*GaussV() == 0), which we exploit below.
 */
MolTemplate co_make_mol_template(int typeIndex, const std::string& name, double diffusion)
{
    MolTemplate tmpl {};
    tmpl.molTypeIndex = typeIndex;
    tmpl.molName = name;
    tmpl.mass = 1.0;
    tmpl.radius = 1.0;
    tmpl.copies = 1;
    tmpl.D = Coord { diffusion, diffusion, diffusion };
    tmpl.Dr = Coord { diffusion * 0.01, diffusion * 0.01, diffusion * 0.01 };
    tmpl.isLipid = false;
    tmpl.isImplicitLipid = false;
    tmpl.isPromoter = false;
    tmpl.isRod = false;
    tmpl.isPoint = true;

    // A single interface, coincident with the centre of mass.
    Interface iface {};
    iface.index = 0;
    iface.name = "i0";
    iface.iCoord = Coord { 0.0, 0.0, 0.0 };
    iface.stateList.emplace_back(std::string("i0"), 0);
    tmpl.interfaceList.push_back(iface);

    return tmpl;
}

/*! \brief Build a molecule with one interface sitting on its centre of mass. */
Molecule co_make_molecule(int index, int comIndex, int molTypeIndex, const Coord& com)
{
    Molecule mol {};
    mol.index = index;
    mol.id = index;
    mol.complexId = comIndex;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mySubVolIndex = 0;
    mol.mass = 1.0;  // must be non-zero: update_properties() divides by total mass
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isImplicitLipid = false;
    mol.isLipid = false;
    mol.isPromoter = false;
    mol.isGhosted = false;
    mol.isShared = false;
    mol.linksToSurface = 0;
    mol.trajStatus = TrajStatus::none;

    Molecule::Iface iface {};
    iface.coord = com;
    iface.index = 0;
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a complex owning the given member molecules. */
Complex co_make_complex(int index, const std::vector<int>& members,
    const std::vector<Molecule>& moleculeList, const MolTemplate& tmpl, std::size_t numMolTypes)
{
    Complex com {};
    com.index = index;
    com.id = index;
    com.ownerRank = 0;
    com.memberList = members;
    com.numEachMol.assign(numMolTypes, 0);
    com.lastNumberUpdateItrEachMol.assign(numMolTypes, 0);

    // Mass-weighted centre of mass of the member molecules.
    double totMass { 0.0 };
    Coord weighted { 0.0, 0.0, 0.0 };
    for (int memMol : members) {
        totMass += moleculeList[memMol].mass;
        weighted.x += moleculeList[memMol].comCoord.x * moleculeList[memMol].mass;
        weighted.y += moleculeList[memMol].comCoord.y * moleculeList[memMol].mass;
        weighted.z += moleculeList[memMol].comCoord.z * moleculeList[memMol].mass;
        ++com.numEachMol[moleculeList[memMol].molTypeIndex];
    }
    com.mass = totMass;
    com.comCoord = Coord { weighted.x / totMass, weighted.y / totMass, weighted.z / totMass };

    com.D = tmpl.D;
    com.Dr = tmpl.Dr;
    com.radius = tmpl.radius;
    com.ncross = 0;
    com.trajStatus = TrajStatus::none;
    com.trajTrans.x = 0.0;
    com.trajTrans.y = 0.0;
    com.trajTrans.z = 0.0;
    com.trajRot = Coord { 0.0, 0.0, 0.0 };
    com.OnSurface = false;
    com.onFiber = false;
    com.tmpOnSurface = false;
    com.isEmpty = false;
    com.linksToSurface = 0;
    com.iLipidIndex = 0;

    return com;
}

/*! \brief Cubic reflecting water box centred on the origin.
 *
 * NOTE: `check_overlap()` scans `RS3Dvect[400 .. 499]` looking for the molecule
 * type index, so the vector *must* be at least 500 entries long. All-zeros means
 * "type 0 matches, RS3D offset 0".
 */
Membrane co_make_box_membrane(double side)
{
    Membrane membrane {};
    membrane.waterBox = Membrane::WaterBox(std::vector<double> { side, side, side });
    membrane.isBox = true;
    membrane.isSphere = false;
    membrane.implicitLipid = false;
    membrane.hasCompartment = false;
    membrane.sphereR = 0.0;
    membrane.xBCtype = "reflect";
    membrane.yBCtype = "reflect";
    membrane.zBCtype = "reflect";
    membrane.RS3Dvect.assign(500, 0.0);
    membrane.totalSA = 6.0 * side * side;
    return membrane;
}

/*! \brief Spherical boundary of the given radius, centred on the origin. */
Membrane co_make_sphere_membrane(double radius)
{
    Membrane membrane {};
    membrane.waterBox
        = Membrane::WaterBox(std::vector<double> { 4.0 * radius, 4.0 * radius, 4.0 * radius });
    membrane.isBox = false;
    membrane.isSphere = true;
    membrane.sphereR = radius;
    membrane.sphereVol = (4.0 / 3.0) * M_PI * radius * radius * radius;
    membrane.implicitLipid = false;
    membrane.hasCompartment = false;
    membrane.xBCtype = "reflect";
    membrane.yBCtype = "reflect";
    membrane.zBCtype = "reflect";
    membrane.RS3Dvect.assign(500, 0.0);
    membrane.totalSA = 4.0 * M_PI * radius * radius;
    return membrane;
}

/*! \brief Radial distance of a coordinate from the origin. */
double co_radial(const Coord& c) { return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z); }

// -----------------------------------------------------------------------------
// Container bundling everything check_overlap() needs, so each test only has to
// tweak the pieces it actually cares about.
// -----------------------------------------------------------------------------
struct CheckOverlapWorld {
    Parameters params {};
    Membrane membrane {};
    SimulVolume simulVolume {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};
    std::vector<int> region {};

    /*! \brief Invoke the function under test. */
    void run(unsigned simItr)
    {
        check_overlap(region, simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
            backRxns, createDestructRxns, molTemplateList, observablesList, counterArrays, membrane,
            co_dummy_mpi_context());
    }
};

} // namespace

// -----------------------------------------------------------------------------
// Test 1: molecules flagged as empty / implicit lipid / already propagated /
//         ghosted must be skipped, leaving positions AND pending trajectories
//         completely untouched.
// -----------------------------------------------------------------------------
void test_check_overlap_skips_flagged_molecules()
{
    std::cerr << "\n[TEST] test_check_overlap_skips_flagged_molecules\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap()\n"
              << "  Scenario:      region holds four molecules, each flagged with one of the\n"
              << "                 four skip conditions (isEmpty, isImplicitLipid,\n"
              << "                 TrajStatus::propagated, isGhosted). Every parent complex\n"
              << "                 carries a large pending translation of (5,5,5).\n"
              << "  Pass criteria: no coordinate moves and no pending trajTrans is consumed,\n"
              << "                 proving the `continue` guard fired for each molecule.\n";

    co_init_rng();

    CheckOverlapWorld world;
    world.params.timeStep = 0.1;
    world.params.clusterOverlapCheck = false;
    world.membrane = co_make_box_membrane(100.0);
    world.molTemplateList.push_back(co_make_mol_template(0, "A", 1.0));

    // Four independent single-molecule complexes at distinct, easy-to-check spots.
    const std::vector<Coord> starts {
        Coord { 0.0, 0.0, 0.0 }, Coord { 5.0, 0.0, 0.0 },
        Coord { 0.0, 5.0, 0.0 }, Coord { 0.0, 0.0, 5.0 }
    };
    for (int i = 0; i < 4; ++i)
        world.moleculeList.push_back(co_make_molecule(i, i, 0, starts[i]));
    for (int i = 0; i < 4; ++i) {
        world.complexList.push_back(co_make_complex(
            i, std::vector<int> { i }, world.moleculeList, world.molTemplateList[0], 1));
        world.complexList[i].ncross = 0;
        world.complexList[i].trajTrans.x = 5.0;
        world.complexList[i].trajTrans.y = 5.0;
        world.complexList[i].trajTrans.z = 5.0;
    }

    // Apply one skip condition per molecule.
    world.moleculeList[0].isEmpty = true;
    world.moleculeList[1].isImplicitLipid = true;
    world.moleculeList[2].trajStatus = TrajStatus::propagated;
    world.moleculeList[3].isGhosted = true;

    world.region = { 0, 1, 2, 3 };

    std::cerr << "  Calling check_overlap()...\n";
    world.run(1u);

    for (int i = 0; i < 4; ++i) {
        std::cerr << "  molecule " << i << " comCoord = (" << world.moleculeList[i].comCoord.x
                  << ", " << world.moleculeList[i].comCoord.y << ", "
                  << world.moleculeList[i].comCoord.z << ")\n";
        EXPECT_DOUBLE_EQ(world.moleculeList[i].comCoord.x, starts[i].x)
            << "skipped molecule " << i << " must not move in x";
        EXPECT_DOUBLE_EQ(world.moleculeList[i].comCoord.y, starts[i].y)
            << "skipped molecule " << i << " must not move in y";
        EXPECT_DOUBLE_EQ(world.moleculeList[i].comCoord.z, starts[i].z)
            << "skipped molecule " << i << " must not move in z";

        // A skipped molecule never reaches Complex::propagate(), so the pending
        // translation must still be sitting on the complex.
        EXPECT_DOUBLE_EQ(world.complexList[i].trajTrans.x, 5.0)
            << "pending trajTrans of complex " << i << " must be preserved";
        EXPECT_DOUBLE_EQ(world.complexList[i].trajTrans.y, 5.0)
            << "pending trajTrans of complex " << i << " must be preserved";
        EXPECT_DOUBLE_EQ(world.complexList[i].trajTrans.z, 5.0)
            << "pending trajTrans of complex " << i << " must be preserved";
    }

    // The skip must not silently rewrite the status flags either.
    EXPECT_EQ(world.moleculeList[0].trajStatus, TrajStatus::none)
        << "empty molecule keeps TrajStatus::none";
    EXPECT_EQ(world.moleculeList[2].trajStatus, TrajStatus::propagated)
        << "already-propagated molecule keeps TrajStatus::propagated";
}

// -----------------------------------------------------------------------------
// Test 2: ncross == 0 and TrajStatus::canBeResampled in a box -> the already
//         sampled trajectory is applied verbatim by Complex::propagate().
// -----------------------------------------------------------------------------
void test_check_overlap_applies_sampled_trajectory_box()
{
    std::cerr << "\n[TEST] test_check_overlap_applies_sampled_trajectory_box\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() -> Complex::propagate() (box branch)\n"
              << "  Scenario:      single molecule, complex.ncross == 0, molecule already has\n"
              << "                 TrajStatus::canBeResampled and a pre-sampled translation\n"
              << "                 of (1, -2, 0.5) with zero rotation.\n"
              << "  Pass criteria: molecule COM and interface are shifted by exactly that\n"
              << "                 translation, the trajectory buffers are zeroed and both\n"
              << "                 molecule and complex end up TrajStatus::propagated.\n";

    co_init_rng();

    CheckOverlapWorld world;
    world.params.timeStep = 0.1;
    world.params.clusterOverlapCheck = false;
    world.membrane = co_make_box_membrane(100.0);
    world.molTemplateList.push_back(co_make_mol_template(0, "A", 1.0));

    world.moleculeList.push_back(co_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }));
    world.complexList.push_back(
        co_make_complex(0, std::vector<int> { 0 }, world.moleculeList, world.molTemplateList[0], 1));

    // The trajectory was "already sampled" elsewhere in the timestep.
    world.moleculeList[0].trajStatus = TrajStatus::canBeResampled;
    world.complexList[0].ncross = 0;
    world.complexList[0].trajTrans.x = 1.0;
    world.complexList[0].trajTrans.y = -2.0;
    world.complexList[0].trajTrans.z = 0.5;
    world.complexList[0].trajRot = Coord { 0.0, 0.0, 0.0 }; // keeps the move purely translational

    world.region = { 0 };

    std::cerr << "  Calling check_overlap()...\n";
    world.run(1u);

    std::cerr << "  molecule COM after propagation = (" << world.moleculeList[0].comCoord.x << ", "
              << world.moleculeList[0].comCoord.y << ", " << world.moleculeList[0].comCoord.z
              << "), expected (1, -2, 0.5)\n";

    // Exact equality is legitimate here: with trajRot == 0 the propagate() fast
    // path performs a plain coordinate + trajTrans addition.
    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.x, 1.0) << "COM x must advance by trajTrans.x";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.y, -2.0) << "COM y must advance by trajTrans.y";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.z, 0.5) << "COM z must advance by trajTrans.z";

    EXPECT_DOUBLE_EQ(world.moleculeList[0].interfaceList[0].coord.x, 1.0)
        << "interface must be carried along with the COM (x)";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].interfaceList[0].coord.y, -2.0)
        << "interface must be carried along with the COM (y)";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].interfaceList[0].coord.z, 0.5)
        << "interface must be carried along with the COM (z)";

    // update_properties() recomputes the complex COM from its single member.
    EXPECT_DOUBLE_EQ(world.complexList[0].comCoord.x, 1.0) << "complex COM tracks its member";
    EXPECT_DOUBLE_EQ(world.complexList[0].comCoord.y, -2.0) << "complex COM tracks its member";
    EXPECT_DOUBLE_EQ(world.complexList[0].comCoord.z, 0.5) << "complex COM tracks its member";

    // The trajectory buffers are consumed (zeroed) by propagate().
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.x, 0.0) << "trajTrans.x zeroed after propagate";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.y, 0.0) << "trajTrans.y zeroed after propagate";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.z, 0.0) << "trajTrans.z zeroed after propagate";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajRot.x, 0.0) << "trajRot.x zeroed after propagate";

    EXPECT_EQ(world.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "molecule must be marked propagated";
    EXPECT_EQ(world.complexList[0].trajStatus, TrajStatus::propagated)
        << "complex must be marked propagated";
}

// -----------------------------------------------------------------------------
// Test 3: two molecules of the *same* complex both listed in `region`. The
//         complex must be propagated exactly once - the second visit has to be
//         short-circuited by the TrajStatus::propagated guard.
// -----------------------------------------------------------------------------
void test_check_overlap_moves_shared_complex_only_once()
{
    std::cerr << "\n[TEST] test_check_overlap_moves_shared_complex_only_once\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() -> Complex::propagate()\n"
              << "  Scenario:      one complex with two member molecules, both indices present\n"
              << "                 in region, pending translation (1,0,0).\n"
              << "  Pass criteria: both members shift by exactly (1,0,0) - i.e. the shift is\n"
              << "                 applied once, not twice - and both are marked propagated.\n";

    co_init_rng();

    CheckOverlapWorld world;
    world.params.timeStep = 0.1;
    world.params.clusterOverlapCheck = false;
    world.membrane = co_make_box_membrane(100.0);
    world.molTemplateList.push_back(co_make_mol_template(0, "A", 1.0));

    world.moleculeList.push_back(co_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }));
    world.moleculeList.push_back(co_make_molecule(1, 0, 0, Coord { 2.0, 0.0, 0.0 }));
    world.complexList.push_back(co_make_complex(
        0, std::vector<int> { 0, 1 }, world.moleculeList, world.molTemplateList[0], 1));

    world.moleculeList[0].trajStatus = TrajStatus::canBeResampled;
    world.moleculeList[1].trajStatus = TrajStatus::canBeResampled;
    world.complexList[0].ncross = 0;
    world.complexList[0].trajTrans.x = 1.0;
    world.complexList[0].trajRot = Coord { 0.0, 0.0, 0.0 };

    world.region = { 0, 1 };

    std::cerr << "  Calling check_overlap()...\n";
    world.run(1u);

    std::cerr << "  member 0 x = " << world.moleculeList[0].comCoord.x << " (expected 1)\n";
    std::cerr << "  member 1 x = " << world.moleculeList[1].comCoord.x << " (expected 3)\n";

    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.x, 1.0)
        << "first member must be shifted exactly once";
    EXPECT_DOUBLE_EQ(world.moleculeList[1].comCoord.x, 3.0)
        << "second member must be shifted exactly once";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.y, 0.0) << "no motion expected in y";
    EXPECT_DOUBLE_EQ(world.moleculeList[1].comCoord.z, 0.0) << "no motion expected in z";

    // Mass-weighted midpoint of x = 1 and x = 3.
    EXPECT_DOUBLE_EQ(world.complexList[0].comCoord.x, 2.0)
        << "complex COM must be the mass-weighted centroid of its members";

    EXPECT_EQ(world.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "member 0 marked propagated";
    EXPECT_EQ(world.moleculeList[1].trajStatus, TrajStatus::propagated)
        << "member 1 marked propagated";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.x, 0.0)
        << "trajTrans consumed by the single propagate() call";
}

// -----------------------------------------------------------------------------
// Test 4: ncross == 0 and TrajStatus::none -> check_overlap() must first sample
//         a propagation vector (create_complex_propagation_vectors) and then
//         propagate. With all diffusion constants set to zero the sampled vector
//         is identically zero, which makes the outcome deterministic.
// -----------------------------------------------------------------------------
void test_check_overlap_samples_trajectory_for_unmoved_complex()
{
    std::cerr << "\n[TEST] test_check_overlap_samples_trajectory_for_unmoved_complex\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() -> create_complex_propagation_vectors()\n"
              << "                 followed by Complex::propagate()\n"
              << "  Scenario:      complex.ncross == 0 and the molecule is still\n"
              << "                 TrajStatus::none, so a fresh trajectory is drawn. All\n"
              << "                 diffusion constants are zero => the draw is exactly zero.\n"
              << "  Pass criteria: the molecule does not move, yet it is promoted to\n"
              << "                 TrajStatus::propagated and the trajectory buffers are clear,\n"
              << "                 proving the sample+propagate path really ran.\n";

    co_init_rng();

    CheckOverlapWorld world;
    world.params.timeStep = 0.1;
    world.params.clusterOverlapCheck = false;
    world.membrane = co_make_box_membrane(100.0);
    // Zero diffusion: sqrt(2*dt*D)*GaussV() == 0 for translation and rotation.
    world.molTemplateList.push_back(co_make_mol_template(0, "A", 0.0));

    world.moleculeList.push_back(co_make_molecule(0, 0, 0, Coord { 3.0, -4.0, 2.0 }));
    world.complexList.push_back(
        co_make_complex(0, std::vector<int> { 0 }, world.moleculeList, world.molTemplateList[0], 1));

    world.moleculeList[0].trajStatus = TrajStatus::none;  // nothing sampled yet
    world.complexList[0].ncross = 0;
    world.complexList[0].D = Coord { 0.0, 0.0, 0.0 };
    world.complexList[0].Dr = Coord { 0.0, 0.0, 0.0 };

    world.region = { 0 };

    std::cerr << "  Calling check_overlap()...\n";
    world.run(1u);

    std::cerr << "  molecule COM after zero-diffusion propagation = ("
              << world.moleculeList[0].comCoord.x << ", " << world.moleculeList[0].comCoord.y << ", "
              << world.moleculeList[0].comCoord.z << "), expected (3, -4, 2)\n";

    // A tiny tolerance guards against harmless round-trip arithmetic inside
    // update_properties(); the physical displacement must still be zero.
    EXPECT_NEAR(world.moleculeList[0].comCoord.x, 3.0, 1e-12)
        << "zero diffusion must produce zero displacement in x";
    EXPECT_NEAR(world.moleculeList[0].comCoord.y, -4.0, 1e-12)
        << "zero diffusion must produce zero displacement in y";
    EXPECT_NEAR(world.moleculeList[0].comCoord.z, 2.0, 1e-12)
        << "zero diffusion must produce zero displacement in z";

    EXPECT_EQ(world.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "molecule must be promoted from none to propagated";
    EXPECT_EQ(world.complexList[0].trajStatus, TrajStatus::propagated)
        << "complex must be marked propagated";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.x, 0.0) << "trajTrans cleared after propagate";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.y, 0.0) << "trajTrans cleared after propagate";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.z, 0.0) << "trajTrans cleared after propagate";
}

// -----------------------------------------------------------------------------
// Test 5: spherical system, complex deep inside the sphere -> the sphere branch
//         propagates and then calls reflect_complex_rad_rot(), which must be a
//         no-op for an interior complex.
// -----------------------------------------------------------------------------
void test_check_overlap_sphere_interior_translation()
{
    std::cerr << "\n[TEST] test_check_overlap_sphere_interior_translation\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() (isSphere == true branch)\n"
              << "  Scenario:      complex sits at the centre of a sphere of radius 100 with a\n"
              << "                 pre-sampled translation of (1,1,1).\n"
              << "  Pass criteria: the translation is applied and the subsequent spherical\n"
              << "                 reflection leaves the interior complex alone.\n";

    co_init_rng();

    CheckOverlapWorld world;
    world.params.timeStep = 0.1;
    world.params.clusterOverlapCheck = false;
    world.membrane = co_make_sphere_membrane(100.0);
    world.molTemplateList.push_back(co_make_mol_template(0, "A", 1.0));

    world.moleculeList.push_back(co_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }));
    world.complexList.push_back(
        co_make_complex(0, std::vector<int> { 0 }, world.moleculeList, world.molTemplateList[0], 1));

    world.moleculeList[0].trajStatus = TrajStatus::canBeResampled;
    world.complexList[0].ncross = 0;
    world.complexList[0].OnSurface = false;  // free 3D complex, not membrane bound
    world.complexList[0].trajTrans.x = 1.0;
    world.complexList[0].trajTrans.y = 1.0;
    world.complexList[0].trajTrans.z = 1.0;
    world.complexList[0].trajRot = Coord { 0.0, 0.0, 0.0 };

    world.region = { 0 };

    std::cerr << "  Calling check_overlap()...\n";
    world.run(1u);

    const Coord& com = world.moleculeList[0].comCoord;
    std::cerr << "  molecule COM = (" << com.x << ", " << com.y << ", " << com.z
              << "), radial distance = " << co_radial(com) << " (sphere R = "
              << world.membrane.sphereR << ")\n";

    EXPECT_NEAR(com.x, 1.0, 1e-9) << "interior complex must simply translate in x";
    EXPECT_NEAR(com.y, 1.0, 1e-9) << "interior complex must simply translate in y";
    EXPECT_NEAR(com.z, 1.0, 1e-9) << "interior complex must simply translate in z";
    EXPECT_LE(co_radial(com), world.membrane.sphereR + 1e-6)
        << "complex must remain inside the spherical boundary";
    EXPECT_EQ(world.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "molecule must be marked propagated";
}

// -----------------------------------------------------------------------------
// Test 6: spherical system, the pre-sampled step would push the complex outside
//         the sphere -> reflect_complex_rad_rot() must pull it back in.
// -----------------------------------------------------------------------------
void test_check_overlap_sphere_reflects_escaping_complex()
{
    std::cerr << "\n[TEST] test_check_overlap_sphere_reflects_escaping_complex\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() -> reflect_complex_rad_rot() (sphere)\n"
              << "  Scenario:      complex starts at x = 98 inside a sphere of radius 100 and\n"
              << "                 carries a translation of +4 in x, which would place it at\n"
              << "                 x = 102, outside the boundary.\n"
              << "  Pass criteria: after the call the complex is back inside the sphere\n"
              << "                 (radial distance <= sphere radius) and it is no longer at\n"
              << "                 the naive, unreflected position.\n";

    co_init_rng();

    CheckOverlapWorld world;
    world.params.timeStep = 0.1;
    world.params.clusterOverlapCheck = false;
    world.membrane = co_make_sphere_membrane(100.0);
    world.molTemplateList.push_back(co_make_mol_template(0, "A", 1.0));

    world.moleculeList.push_back(co_make_molecule(0, 0, 0, Coord { 98.0, 0.0, 0.0 }));
    world.complexList.push_back(
        co_make_complex(0, std::vector<int> { 0 }, world.moleculeList, world.molTemplateList[0], 1));

    world.moleculeList[0].trajStatus = TrajStatus::canBeResampled;
    world.complexList[0].ncross = 0;
    world.complexList[0].OnSurface = false;
    world.complexList[0].trajTrans.x = 4.0;  // 98 + 4 = 102 > R = 100
    world.complexList[0].trajRot = Coord { 0.0, 0.0, 0.0 };

    world.region = { 0 };

    std::cerr << "  Calling check_overlap()...\n";
    world.run(1u);

    const Coord& com = world.moleculeList[0].comCoord;
    const double radial = co_radial(com);
    std::cerr << "  molecule COM = (" << com.x << ", " << com.y << ", " << com.z
              << "), radial distance = " << radial << " (sphere R = " << world.membrane.sphereR
              << ")\n";

    EXPECT_LE(radial, world.membrane.sphereR + 1e-6)
        << "the escaping complex must be reflected back inside the sphere";
    EXPECT_LT(com.x, 102.0)
        << "the naive unreflected position (x = 102) must have been corrected";
    EXPECT_EQ(world.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "molecule must still be marked propagated";
}

// -----------------------------------------------------------------------------
// Test 7: guard ordering. A molecule whose complex reports crossings
//         (ncross > 0) but which was already propagated must be skipped by the
//         very first `continue`, so the overlap sweepers are never reached.
//
//         This is deliberately the only ncross > 0 case exercised here: the
//         sweep routines require complete crossbase/mycrossint/crossrxn/probvec
//         bookkeeping on every member molecule, which cannot be faked safely.
// -----------------------------------------------------------------------------
void test_check_overlap_ncross_guard_ordering()
{
    std::cerr << "\n[TEST] test_check_overlap_ncross_guard_ordering\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() (skip guard vs. ncross > 0 branch)\n"
              << "  Scenario:      complex.ncross == 3 (crossings were registered) but the\n"
              << "                 molecule is already TrajStatus::propagated.\n"
              << "  Pass criteria: the molecule is skipped before the sweep branch, so its\n"
              << "                 coordinates, status and the complex's ncross are unchanged.\n";

    co_init_rng();

    CheckOverlapWorld world;
    world.params.timeStep = 0.1;
    world.params.clusterOverlapCheck = false;
    world.membrane = co_make_box_membrane(100.0);
    world.molTemplateList.push_back(co_make_mol_template(0, "A", 1.0));

    world.moleculeList.push_back(co_make_molecule(0, 0, 0, Coord { 7.0, 8.0, 9.0 }));
    world.complexList.push_back(
        co_make_complex(0, std::vector<int> { 0 }, world.moleculeList, world.molTemplateList[0], 1));

    world.moleculeList[0].trajStatus = TrajStatus::propagated;  // fires the first guard
    world.complexList[0].ncross = 3;                            // would otherwise enter the sweep
    world.complexList[0].trajTrans.x = 12.0;

    world.region = { 0 };

    std::cerr << "  Calling check_overlap()...\n";
    world.run(1u);

    std::cerr << "  molecule COM = (" << world.moleculeList[0].comCoord.x << ", "
              << world.moleculeList[0].comCoord.y << ", " << world.moleculeList[0].comCoord.z
              << "), expected (7, 8, 9)\n";

    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.x, 7.0) << "propagated molecule must not move";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.y, 8.0) << "propagated molecule must not move";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.z, 9.0) << "propagated molecule must not move";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.x, 12.0)
        << "pending trajTrans must not be consumed";
    EXPECT_EQ(world.complexList[0].ncross, 3) << "check_overlap must not alter ncross";
    EXPECT_EQ(world.moleculeList[0].trajStatus, TrajStatus::propagated)
        << "TrajStatus must be left as-is";
}

// -----------------------------------------------------------------------------
// Test 8: an empty region is a well-defined no-op (loop body never executes).
// -----------------------------------------------------------------------------
void test_check_overlap_empty_region_is_noop()
{
    std::cerr << "\n[TEST] test_check_overlap_empty_region_is_noop\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap()\n"
              << "  Scenario:      region is empty although a movable molecule with a pending\n"
              << "                 trajectory exists in the system.\n"
              << "  Pass criteria: nothing at all happens - the function iterates over region,\n"
              << "                 not over moleculeList.\n";

    co_init_rng();

    CheckOverlapWorld world;
    world.params.timeStep = 0.1;
    world.params.clusterOverlapCheck = false;
    world.membrane = co_make_box_membrane(100.0);
    world.molTemplateList.push_back(co_make_mol_template(0, "A", 1.0));

    world.moleculeList.push_back(co_make_molecule(0, 0, 0, Coord { -1.0, -2.0, -3.0 }));
    world.complexList.push_back(
        co_make_complex(0, std::vector<int> { 0 }, world.moleculeList, world.molTemplateList[0], 1));

    world.moleculeList[0].trajStatus = TrajStatus::canBeResampled;
    world.complexList[0].ncross = 0;
    world.complexList[0].trajTrans.x = 9.0;

    world.region.clear();  // nothing to process

    std::cerr << "  Calling check_overlap() with an empty region...\n";
    world.run(1u);

    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.x, -1.0)
        << "molecule outside region must not move";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.y, -2.0)
        << "molecule outside region must not move";
    EXPECT_DOUBLE_EQ(world.moleculeList[0].comCoord.z, -3.0)
        << "molecule outside region must not move";
    EXPECT_DOUBLE_EQ(world.complexList[0].trajTrans.x, 9.0)
        << "pending trajTrans must survive an empty-region call";
    EXPECT_EQ(world.moleculeList[0].trajStatus, TrajStatus::canBeResampled)
        << "TrajStatus must be untouched for a molecule outside the region";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - one per scenario so every case is reported separately
// and a failure in one does not stop the others from running.
// -----------------------------------------------------------------------------
TEST(CheckOverlapTest, SkipsFlaggedMolecules) { test_check_overlap_skips_flagged_molecules(); }
TEST(CheckOverlapTest, AppliesSampledTrajectoryBox)
{
    test_check_overlap_applies_sampled_trajectory_box();
}
TEST(CheckOverlapTest, MovesSharedComplexOnlyOnce)
{
    test_check_overlap_moves_shared_complex_only_once();
}
TEST(CheckOverlapTest, SamplesTrajectoryForUnmovedComplex)
{
    test_check_overlap_samples_trajectory_for_unmoved_complex();
}
TEST(CheckOverlapTest, SphereInteriorTranslation) { test_check_overlap_sphere_interior_translation(); }
TEST(CheckOverlapTest, SphereReflectsEscapingComplex)
{
    test_check_overlap_sphere_reflects_escaping_complex();
}
TEST(CheckOverlapTest, NcrossGuardOrdering) { test_check_overlap_ncross_guard_ordering(); }
TEST(CheckOverlapTest, EmptyRegionIsNoop) { test_check_overlap_empty_region_is_noop(); }