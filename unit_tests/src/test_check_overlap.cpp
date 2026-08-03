/*! \file test_check_overlap.cpp
 *
 * ### Unit test for src/reactions/check_overlap.cpp
 *
 * The file under test contains a single function:
 *
 * \code
 * void check_overlap(std::vector<int>& region, unsigned simItr,
 *                    Parameters& params, std::vector<Molecule>& moleculeList,
 *                    std::vector<Complex>& complexList, SimulVolume& simulVolume,
 *                    std::vector<ForwardRxn>& forwardRxns,
 *                    std::vector<BackRxn>& backRxns,
 *                    std::vector<CreateDestructRxn>& createDestructRxns,
 *                    std::vector<MolTemplate>& molTemplateList,
 *                    std::map<std::string, int>& observablesList,
 *                    copyCounters& counterArrays, Membrane& membraneObject,
 *                    MpiContext& mpiContext);
 * \endcode
 *
 * `check_overlap()` is the final position-update stage of one simulation
 * iteration.  For every molecule index listed in `region` it:
 *
 *   1. skips molecules that are empty / implicit lipid / already propagated /
 *      ghosted,
 *   2. looks up the reflecting-surface value `RS3Dinput` from
 *      `membraneObject.RS3Dvect` (entries 300..499 of that table),
 *   3. if the parent complex had potential reaction partners (`ncross > 0`)
 *      delegates to one of the `sweep_separation_complex_rot*` routines
 *      (membrane / cluster / 3D flavours), and
 *   4. otherwise creates propagation vectors (only when the trajectory has not
 *      been sampled yet) and propagates the complex, reflecting it back into a
 *      spherical volume when the boundary is a sphere.
 *
 * The tests below therefore verify observable *behaviour* of each of those
 * branches:
 *   - inactive molecules are left untouched,
 *   - an empty region is a no-op,
 *   - a molecule with a pre-sampled trajectory is translated by exactly that
 *     trajectory (fully deterministic, no RNG involved),
 *   - a molecule whose trajectory has not been sampled gets one (RNG driven)
 *     and ends up inside the simulation volume,
 *   - members of the same complex are only propagated once,
 *   - the sphere branch never leaves the spherical boundary,
 *   - the "ncross > 0" sweep branches (2D membrane, 2D cluster, 3D) run to
 *     completion and leave the system in a sane state.
 *
 * Every test prints, to stderr, which source file / function is under test,
 * what the scenario is, and what the pass criteria are.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace => no symbol collisions with the rest of
// the generated test suite; names additionally carry a "co_" prefix).
// -----------------------------------------------------------------------------
namespace {

//! check_overlap() indexes membraneObject.RS3Dvect up to element 499, so the
//! lookup table must be at least 500 entries long.
constexpr std::size_t kCoRS3DSize = 501;

//! Storage used to fabricate the MpiContext reference required by the
//! signature.  check_overlap() only ever touches `mpiContext` inside blocks
//! guarded by the compile-time DEBUG flag (which is `false` in this build), so
//! the reference is never dereferenced.  A raw byte buffer is used because
//! MpiContext is only forward declared in the public headers.
alignas(16) unsigned char co_mpi_context_storage[8192] = { 0 };

/*! \brief Returns a never-dereferenced placeholder MpiContext reference. */
MpiContext& co_dummy_mpi_context()
{
    return *reinterpret_cast<MpiContext*>(co_mpi_context_storage);
}

/*! \brief Make sure the global GSL RNG exists before any propagation happens.
 *
 * gtest_main.cpp defines `gsl_rng* r = nullptr;`, so the first test that needs
 * random numbers has to seed the generator.  A fixed seed keeps runs
 * reproducible.
 */
void co_ensure_rng()
{
    if (r == nullptr) {
        std::cerr << "  [setup] global GSL RNG was null -> seeding with srand_gsl(17)\n";
        srand_gsl(17);
    }
}

/*! \brief Everything check_overlap() needs, bundled so tests stay readable. */
struct CoScenario {
    Parameters params {};
    Membrane membraneObject {};
    SimulVolume simulVolume {};
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList {};
    std::vector<MolTemplate> molTemplateList {};
    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};
    std::vector<CreateDestructRxn> createDestructRxns {};
    std::map<std::string, int> observablesList {};
    copyCounters counterArrays {};

    /*! \brief Thin wrapper so each test calls the function under test once. */
    void run(std::vector<int>& region, unsigned simItr = 1)
    {
        check_overlap(region, simItr, params, moleculeList, complexList, simulVolume,
            forwardRxns, backRxns, createDestructRxns, molTemplateList,
            observablesList, counterArrays, membraneObject, co_dummy_mpi_context());
    }
};

/*! \brief Build a minimal but self-consistent simulation state.
 *
 * \param[in] isSphere true -> spherical boundary of radius 100 nm,
 *                     false -> reflecting cubic box of side 200 nm.
 */
CoScenario co_make_scenario(bool isSphere)
{
    co_ensure_rng();

    CoScenario s;

    // --- simulation parameters -----------------------------------------------
    s.params.timeStep = 0.1; // microseconds
    s.params.numMolTypes = 1;
    s.params.clusterOverlapCheck = false;
    s.params.overlapSepLimit = 0.1;
    s.params.rMaxLimit = 10.0;
    s.params.scaleMaxDisplace = 1e6; // never veto a displacement

    // --- boundary ------------------------------------------------------------
    s.membraneObject.isSphere = isSphere;
    s.membraneObject.isBox = !isSphere;
    s.membraneObject.implicitLipid = false;
    s.membraneObject.xBCtype = "reflect";
    s.membraneObject.yBCtype = "reflect";
    s.membraneObject.zBCtype = "reflect";
    s.membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 200.0, 200.0, 200.0 });
    if (isSphere) {
        s.membraneObject.sphereR = 100.0;
        s.membraneObject.sphereVol = 4.0 / 3.0 * M_PI * 100.0 * 100.0 * 100.0;
    }

    // RS3D lookup table: entries 400..499 hold molTypeIndex tags, 300..399 hold
    // the corresponding reflecting-surface distances.  Tag with -1 so that no
    // molecule type matches and RS3Dinput stays 0 (no implicit lipid here).
    s.membraneObject.RS3Dvect.assign(kCoRS3DSize, 0.0);
    for (std::size_t i = 400; i < 500; ++i)
        s.membraneObject.RS3Dvect[i] = -1.0;

    // --- one molecule template ----------------------------------------------
    MolTemplate temp;
    temp.molName = "A";
    temp.molTypeIndex = 0;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord(10.0, 10.0, 10.0);
    temp.Dr = Coord(0.01, 0.01, 0.01);
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPoint = false;
    temp.isRod = false;
    temp.checkOverlap = false;

    Interface iface;
    iface.index = 0;
    iface.name = "a";
    iface.iCoord = Coord(1.0, 0.0, 0.0);
    iface.stateList.push_back(Interface::State(0));
    temp.interfaceList.push_back(iface);
    s.molTemplateList.push_back(temp);

    // --- one bimolecular reaction so the sweep routines can look up a
    //     binding radius for any recorded crossing pair --------------------
    ForwardRxn fwd;
    fwd.bindRadius = 1.0;
    fwd.bindRadius2D = 1.0;
    fwd.absRxnIndex = 0;
    fwd.relRxnIndex = 0;
    fwd.rxnType = ReactionType::bimolecular;
    fwd.isReversible = false;
    fwd.rateList.emplace_back();
    fwd.rateList[0].rate = 1.0;
    s.forwardRxns.push_back(fwd);

    return s;
}

/*! \brief Append a single-interface molecule at \p com to the scenario. */
int co_add_molecule(CoScenario& s, int comIndex, const Coord& com)
{
    Molecule mol;
    mol.index = static_cast<int>(s.moleculeList.size());
    mol.id = mol.index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.comCoord = com;
    mol.isEmpty = false;
    mol.isImplicitLipid = false;
    mol.isGhosted = false;
    mol.isLipid = false;
    mol.trajStatus = TrajStatus::none;
    mol.mySubVolIndex = 0;

    Molecule::Iface mIface;
    mIface.coord = Coord(com.x + 1.0, com.y, com.z);
    mIface.index = 0;
    mIface.relIndex = 0;
    mIface.stateIndex = 0;
    mIface.molTypeIndex = 0;
    mIface.isBound = false;
    mol.interfaceList.push_back(mIface);
    mol.freelist.push_back(0);

    s.moleculeList.push_back(mol);
    return mol.index;
}

/*! \brief Append a complex owning \p members, with z-diffusion \p dz.
 *
 * `dz == 0` marks a membrane-bound (2D) complex, which is what selects the
 * `sweep_separation_complex_rot_memtest*` branch inside check_overlap().
 */
int co_add_complex(CoScenario& s, const Coord& com, const std::vector<int>& members, double dz)
{
    Complex cplx;
    cplx.index = static_cast<int>(s.complexList.size());
    cplx.id = cplx.index;
    cplx.ownerRank = 0;
    cplx.comCoord = com;
    cplx.radius = 2.0;
    cplx.mass = static_cast<double>(members.size());
    cplx.memberList = members;
    cplx.numEachMol = std::vector<int> { static_cast<int>(members.size()) };
    cplx.lastNumberUpdateItrEachMol = std::vector<long long int> { 0 };
    cplx.D = Coord(10.0, 10.0, dz);
    cplx.Dr = Coord(0.01, 0.01, 0.01);
    cplx.isEmpty = false;
    cplx.OnSurface = false;
    cplx.tmpOnSurface = false;
    cplx.linksToSurface = 0;
    cplx.iLipidIndex = 0;
    cplx.ncross = 0;
    cplx.trajStatus = TrajStatus::none;
    cplx.trajTrans.x = 0.0;
    cplx.trajTrans.y = 0.0;
    cplx.trajTrans.z = 0.0;
    cplx.trajRot.x = 0.0;
    cplx.trajRot.y = 0.0;
    cplx.trajRot.z = 0.0;

    s.complexList.push_back(cplx);
    return cplx.index;
}

/*! \brief Record a (non-reacted) crossing between two molecules.
 *
 * This is what makes `complexList[...].ncross > 0`, i.e. it selects the
 * "sweep separation" branch of check_overlap().  The two molecules are placed
 * far apart by the caller so no real overlap has to be resolved.
 */
void co_add_crossing_pair(CoScenario& s, int molA, int molB)
{
    s.moleculeList[molA].crossbase.push_back(molB);
    s.moleculeList[molA].mycrossint.push_back(0);
    s.moleculeList[molA].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    s.moleculeList[molA].probvec.push_back(0.0);

    s.moleculeList[molB].crossbase.push_back(molA);
    s.moleculeList[molB].mycrossint.push_back(0);
    s.moleculeList[molB].crossrxn.push_back(std::array<int, 3> { 0, 0, 0 });
    s.moleculeList[molB].probvec.push_back(0.0);

    s.complexList[s.moleculeList[molA].myComIndex].ncross = 1;
    s.complexList[s.moleculeList[molB].myComIndex].ncross = 1;
}

/*! \brief Is \p c inside the (reflecting) water box, within a tolerance? */
bool co_inside_box(const Coord& c, const Membrane& m)
{
    const double tol = 1e-6;
    return std::abs(c.x) <= m.waterBox.x / 2.0 + tol
        && std::abs(c.y) <= m.waterBox.y / 2.0 + tol
        && std::abs(c.z) <= m.waterBox.z / 2.0 + tol;
}

/*! \brief Euclidean distance between two coordinates. */
double co_distance(const Coord& a, const Coord& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/*! \brief True when no component of \p c is NaN or infinite. */
bool co_is_finite(const Coord& c)
{
    return std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z);
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the guard clause at the top of the loop.
// -----------------------------------------------------------------------------
void test_co_skips_inactive_molecules()
{
    std::cerr << "\n[TEST] test_co_skips_inactive_molecules\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap()\n"
              << "  Scenario:      four molecules in the region, each with a\n"
              << "                 different reason to be skipped (isEmpty,\n"
              << "                 isImplicitLipid, already propagated, ghosted).\n"
              << "  Pass criteria: no coordinate and no trajStatus is modified.\n";

    CoScenario s = co_make_scenario(/*isSphere=*/false);

    // Build four independent single-molecule complexes spaced 20 nm apart.
    for (int i = 0; i < 4; ++i) {
        Coord pos(-30.0 + 20.0 * i, 0.0, 0.0);
        int molIdx = co_add_molecule(s, i, pos);
        co_add_complex(s, pos, std::vector<int> { molIdx }, /*dz=*/10.0);
    }

    // Flag each molecule so that the guard clause fires.
    s.moleculeList[0].isEmpty = true;
    s.moleculeList[1].isImplicitLipid = true;
    s.moleculeList[2].trajStatus = TrajStatus::propagated;
    s.moleculeList[3].isGhosted = true;

    // Snapshot the starting coordinates.
    std::vector<Coord> before;
    for (const auto& mol : s.moleculeList)
        before.push_back(mol.comCoord);

    std::vector<int> region { 0, 1, 2, 3 };
    std::cerr << "  Calling check_overlap() with region = {0,1,2,3}...\n";
    s.run(region);

    for (std::size_t i = 0; i < s.moleculeList.size(); ++i) {
        std::cerr << "    mol " << i << " COM before (" << before[i].x << ") after ("
                  << s.moleculeList[i].comCoord.x << ")\n";
        EXPECT_DOUBLE_EQ(s.moleculeList[i].comCoord.x, before[i].x)
            << "skipped molecule " << i << " must not move in x";
        EXPECT_DOUBLE_EQ(s.moleculeList[i].comCoord.y, before[i].y)
            << "skipped molecule " << i << " must not move in y";
        EXPECT_DOUBLE_EQ(s.moleculeList[i].comCoord.z, before[i].z)
            << "skipped molecule " << i << " must not move in z";
    }

    // trajStatus of the skipped molecules must also be untouched.
    EXPECT_EQ(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "empty molecule trajStatus should remain 'none'";
    EXPECT_EQ(static_cast<int>(s.moleculeList[1].trajStatus), static_cast<int>(TrajStatus::none))
        << "implicit lipid trajStatus should remain 'none'";
    EXPECT_EQ(static_cast<int>(s.moleculeList[2].trajStatus), static_cast<int>(TrajStatus::propagated))
        << "already propagated molecule should stay 'propagated'";
    EXPECT_EQ(static_cast<int>(s.moleculeList[3].trajStatus), static_cast<int>(TrajStatus::none))
        << "ghosted molecule trajStatus should remain 'none'";
}

// -----------------------------------------------------------------------------
// Test 2: an empty region list must be a complete no-op.
// -----------------------------------------------------------------------------
void test_co_empty_region_is_noop()
{
    std::cerr << "\n[TEST] test_co_empty_region_is_noop\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap()\n"
              << "  Scenario:      an active, un-propagated molecule exists but the\n"
              << "                 region list is empty.\n"
              << "  Pass criteria: the molecule is not propagated (check_overlap only\n"
              << "                 walks the region, never the whole moleculeList).\n";

    CoScenario s = co_make_scenario(/*isSphere=*/false);
    int molIdx = co_add_molecule(s, 0, Coord(0.0, 0.0, 0.0));
    co_add_complex(s, Coord(0.0, 0.0, 0.0), std::vector<int> { molIdx }, /*dz=*/10.0);

    const Coord before = s.moleculeList[0].comCoord;

    std::vector<int> region {}; // deliberately empty
    std::cerr << "  Calling check_overlap() with an empty region...\n";
    s.run(region);

    EXPECT_DOUBLE_EQ(s.moleculeList[0].comCoord.x, before.x) << "x must be unchanged";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].comCoord.y, before.y) << "y must be unchanged";
    EXPECT_DOUBLE_EQ(s.moleculeList[0].comCoord.z, before.z) << "z must be unchanged";
    EXPECT_EQ(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must stay 'none' when the molecule is not in the region";
}

// -----------------------------------------------------------------------------
// Test 3: ncross == 0 and trajStatus == canBeResampled -> deterministic
//         translation by the already sampled trajTrans.
// -----------------------------------------------------------------------------
void test_co_box_presampled_trajectory_is_applied()
{
    std::cerr << "\n[TEST] test_co_box_presampled_trajectory_is_applied\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() (ncross==0, box branch)\n"
              << "  Scenario:      trajStatus == canBeResampled with trajTrans =\n"
              << "                 (2,-3,1) and zero rotation, so check_overlap must\n"
              << "                 NOT resample and simply propagate.\n"
              << "  Pass criteria: COM and interface both shift by exactly (2,-3,1).\n";

    CoScenario s = co_make_scenario(/*isSphere=*/false);
    const Coord start(5.0, -5.0, 2.0);
    int molIdx = co_add_molecule(s, 0, start);
    co_add_complex(s, start, std::vector<int> { molIdx }, /*dz=*/10.0);

    // Pretend the trajectory was already sampled earlier in the iteration.
    s.moleculeList[0].trajStatus = TrajStatus::canBeResampled;
    s.complexList[0].trajStatus = TrajStatus::canBeResampled;
    s.complexList[0].trajTrans.x = 2.0;
    s.complexList[0].trajTrans.y = -3.0;
    s.complexList[0].trajTrans.z = 1.0;
    // zero rotation keeps the expected result analytic

    const Coord ifaceStart = s.moleculeList[0].interfaceList[0].coord;

    std::vector<int> region { molIdx };
    std::cerr << "  Calling check_overlap()...\n";
    s.run(region);

    std::cerr << "    COM   : (" << s.moleculeList[0].comCoord.x << ", "
              << s.moleculeList[0].comCoord.y << ", " << s.moleculeList[0].comCoord.z << ")\n"
              << "    iface : (" << s.moleculeList[0].interfaceList[0].coord.x << ", "
              << s.moleculeList[0].interfaceList[0].coord.y << ", "
              << s.moleculeList[0].interfaceList[0].coord.z << ")\n";

    EXPECT_NEAR(s.moleculeList[0].comCoord.x, start.x + 2.0, 1e-9) << "COM x shifted by +2";
    EXPECT_NEAR(s.moleculeList[0].comCoord.y, start.y - 3.0, 1e-9) << "COM y shifted by -3";
    EXPECT_NEAR(s.moleculeList[0].comCoord.z, start.z + 1.0, 1e-9) << "COM z shifted by +1";

    EXPECT_NEAR(s.moleculeList[0].interfaceList[0].coord.x, ifaceStart.x + 2.0, 1e-9)
        << "interface x must follow the rigid-body translation";
    EXPECT_NEAR(s.moleculeList[0].interfaceList[0].coord.y, ifaceStart.y - 3.0, 1e-9)
        << "interface y must follow the rigid-body translation";
    EXPECT_NEAR(s.moleculeList[0].interfaceList[0].coord.z, ifaceStart.z + 1.0, 1e-9)
        << "interface z must follow the rigid-body translation";

    // Complex COM must track its single member.
    EXPECT_NEAR(s.complexList[0].comCoord.x, s.moleculeList[0].comCoord.x, 1e-9)
        << "complex COM should coincide with its single member";
    EXPECT_NE(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "molecule should no longer be flagged 'none' after propagation";
}

// -----------------------------------------------------------------------------
// Test 4: ncross == 0 and trajStatus == none -> a trajectory is sampled and
//         applied (RNG driven), and the molecule stays inside the box.
// -----------------------------------------------------------------------------
void test_co_box_samples_trajectory_when_none()
{
    std::cerr << "\n[TEST] test_co_box_samples_trajectory_when_none\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() (ncross==0, box branch)\n"
              << "  Scenario:      trajStatus == none, so check_overlap must call\n"
              << "                 create_complex_propagation_vectors() and propagate.\n"
              << "  Pass criteria: molecule moved, coordinates finite, still in box,\n"
              << "                 trajStatus no longer 'none'.\n";

    CoScenario s = co_make_scenario(/*isSphere=*/false);
    const Coord start(0.0, 0.0, 0.0);
    int molIdx = co_add_molecule(s, 0, start);
    co_add_complex(s, start, std::vector<int> { molIdx }, /*dz=*/10.0);

    std::vector<int> region { molIdx };
    std::cerr << "  Calling check_overlap()...\n";
    s.run(region);

    const Coord after = s.moleculeList[0].comCoord;
    const double displacement = co_distance(after, start);
    std::cerr << "    displacement = " << displacement << " nm, new COM = ("
              << after.x << ", " << after.y << ", " << after.z << ")\n";

    EXPECT_TRUE(co_is_finite(after)) << "propagated coordinates must be finite";
    EXPECT_GT(displacement, 0.0)
        << "a Gaussian displacement with D=10, dt=0.1 must be non-zero";
    EXPECT_TRUE(co_inside_box(after, s.membraneObject))
        << "propagated molecule must remain inside the reflecting box";
    EXPECT_NE(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must be updated once the trajectory has been created";
    EXPECT_NEAR(co_distance(after, s.complexList[0].comCoord), 0.0, 1e-9)
        << "single-member complex COM must track the molecule COM";
}

// -----------------------------------------------------------------------------
// Test 5: two molecules of the SAME complex both listed in the region; the
//         complex must only be propagated once and stay rigid.
// -----------------------------------------------------------------------------
void test_co_shared_complex_propagates_once()
{
    std::cerr << "\n[TEST] test_co_shared_complex_propagates_once\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() (ncross==0, box branch)\n"
              << "  Scenario:      one complex with two member molecules, both listed\n"
              << "                 in the region.  The first pass propagates the\n"
              << "                 complex and marks members 'propagated', so the\n"
              << "                 second member hits the guard clause.\n"
              << "  Pass criteria: internal distance conserved (rigid body) and both\n"
              << "                 molecules end up inside the box.\n";

    CoScenario s = co_make_scenario(/*isSphere=*/false);
    const Coord comA(-1.0, 0.0, 0.0);
    const Coord comB(1.0, 0.0, 0.0);
    int molA = co_add_molecule(s, 0, comA);
    int molB = co_add_molecule(s, 0, comB);
    co_add_complex(s, Coord(0.0, 0.0, 0.0), std::vector<int> { molA, molB }, /*dz=*/10.0);

    const double distBefore = co_distance(comA, comB);

    std::vector<int> region { molA, molB };
    std::cerr << "  Calling check_overlap() with both members in the region...\n";
    s.run(region);

    const double distAfter = co_distance(s.moleculeList[molA].comCoord,
        s.moleculeList[molB].comCoord);
    std::cerr << "    internal distance before = " << distBefore
              << ", after = " << distAfter << "\n";

    EXPECT_NEAR(distAfter, distBefore, 1e-6)
        << "rigid-body propagation must conserve the intra-complex distance";
    EXPECT_TRUE(co_is_finite(s.moleculeList[molA].comCoord)) << "member A coords finite";
    EXPECT_TRUE(co_is_finite(s.moleculeList[molB].comCoord)) << "member B coords finite";
    EXPECT_TRUE(co_inside_box(s.moleculeList[molA].comCoord, s.membraneObject))
        << "member A must stay inside the box";
    EXPECT_TRUE(co_inside_box(s.moleculeList[molB].comCoord, s.membraneObject))
        << "member B must stay inside the box";
    EXPECT_NE(static_cast<int>(s.moleculeList[molB].trajStatus), static_cast<int>(TrajStatus::none))
        << "the second member must have been moved along with its complex";
}

// -----------------------------------------------------------------------------
// Test 6: spherical boundary branch -> reflect_complex_rad_rot() must keep the
//         complex inside the sphere.
// -----------------------------------------------------------------------------
void test_co_sphere_stays_inside_boundary()
{
    std::cerr << "\n[TEST] test_co_sphere_stays_inside_boundary\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() (ncross==0, sphere branch)\n"
              << "  Scenario:      spherical boundary R=100, complex parked near the\n"
              << "                 shell with a pre-sampled outward trajectory.\n"
              << "  Pass criteria: final radial distance <= sphere radius and all\n"
              << "                 coordinates finite.\n";

    CoScenario s = co_make_scenario(/*isSphere=*/true);
    const Coord start(90.0, 0.0, 0.0);
    int molIdx = co_add_molecule(s, 0, start);
    co_add_complex(s, start, std::vector<int> { molIdx }, /*dz=*/10.0);

    // Already-sampled trajectory pointing straight at the shell.
    s.moleculeList[0].trajStatus = TrajStatus::canBeResampled;
    s.complexList[0].trajStatus = TrajStatus::canBeResampled;
    s.complexList[0].trajTrans.x = 5.0;

    std::vector<int> region { molIdx };
    std::cerr << "  Calling check_overlap()...\n";
    s.run(region);

    const Coord after = s.moleculeList[0].comCoord;
    const double radial = std::sqrt(after.x * after.x + after.y * after.y + after.z * after.z);
    std::cerr << "    final radial distance = " << radial
              << " (sphere radius " << s.membraneObject.sphereR << ")\n";

    EXPECT_TRUE(co_is_finite(after)) << "coordinates must stay finite on the sphere branch";
    EXPECT_LE(radial, s.membraneObject.sphereR + 1e-6)
        << "the complex must be reflected back inside the spherical boundary";
    EXPECT_NE(static_cast<int>(s.moleculeList[0].trajStatus), static_cast<int>(TrajStatus::none))
        << "trajStatus must be updated on the sphere branch as well";
}

// -----------------------------------------------------------------------------
// Test 7: ncross > 0 with a 3D complex (|D.z| > 0) -> the
//         sweep_separation_complex_rot() branch.
// -----------------------------------------------------------------------------
void test_co_cross_3d_sweep_branch()
{
    std::cerr << "\n[TEST] test_co_cross_3d_sweep_branch\n"
              << "  Source file:   src/reactions/check_overlap.cpp\n"
              << "  Function:      check_overlap() (ncross>0, |D.z|>0 -> 3D sweep)\n"
              << "  Scenario:      two complexes recorded a crossing but did not\n"
              << "                 react; they are 40 nm apart so no real overlap has\n"
              << "                 to be resolved.\n"
              << "  Pass criteria: the sweep completes, coordinates stay finite and\n"
              << "                 inside the box, and both molecules are moved.\n";

    CoScenario s = co_make_scenario(/*isSphere=*/false);
    const Coord posA(-20.0, 0.0, 0.0);
    const Coord posB(20.0, 0.0, 0.0);
    int molA = co_add_molecule(s, 0, posA);
    co_add_complex(s, posA, std::vector<int> { molA }, /*dz=*/10.0);
    