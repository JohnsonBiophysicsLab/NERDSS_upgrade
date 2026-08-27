/*! \file test_create_complex_propagation_vectors.cpp
 *
 * ### Unit test for src/trajectory_functions/create_complex_propagation_vectors.cpp
 *
 * Function under test:
 *
 *     void create_complex_propagation_vectors(const Parameters& params,
 *                                             Complex& targCom,
 *                                             std::vector<Molecule>& moleculeList,
 *                                             std::vector<Complex>& complexList,
 *                                             const std::vector<MolTemplate>& molTemplateList,
 *                                             const Membrane& membraneObject)
 *
 * Behaviour that is being verified:
 *
 *   1. In a box (or inside a sphere, i.e. Complex::OnSurface == false) the
 *      translational displacement is  sqrt(2*dt*D_i)*GaussV()  and the
 *      rotational displacement is  sqrt(2*dt*Dr_i)*GaussV(), drawn in the
 *      order x, y, z, rx, ry, rz.
 *   2. When a diffusion constant is exactly zero the corresponding component
 *      is set to exactly 0.0 *without* consuming a random number (this is the
 *      optimisation documented in the source file).  We verify that by
 *      re-seeding the RNG and comparing against an explicitly drawn sequence.
 *   3. The generated vectors are statistically correct (zero mean, standard
 *      deviation sqrt(2*dt*D)).
 *   4. After the vectors are generated, reflecting boundary conditions are
 *      applied, so a complex sitting next to a wall never ends up outside the
 *      simulation box.
 *   5. The implicit-lipid look-up table (RS3Dvect) is indexed safely and does
 *      not change the result for a complex far away from the membrane.
 *   6. The "on the sphere surface" branch is taken when
 *      Complex::OnSurface == true and Membrane::isSphere == true, and produces
 *      finite displacements (with an exactly zero rotation when Dr == 0).
 *
 * All console output goes to stderr so it interleaves with gtest output.
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// The random number generator handle is owned by gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Helper: (re)initialise the GSL random number generator with a fixed seed so
// every test is completely reproducible.
// -----------------------------------------------------------------------------
void ccpv_seed_rng()
{
    const gsl_rng_type* T;
    T = gsl_rng_default;
    if (r == nullptr)
        r = gsl_rng_alloc(T);
    gsl_rng_set(r, 42);
}

// -----------------------------------------------------------------------------
// Helper: build a minimal, fully initialised Molecule.
//
// The molecule owns exactly one interface which is coincident with its centre
// of mass, which keeps all the geometry used by the reflecting boundary
// conditions trivially simple.
// -----------------------------------------------------------------------------
Molecule ccpv_make_molecule(const Coord& com, int molTypeIndex)
{
    Molecule mol;
    mol.index = 0;
    mol.id = 0;
    mol.myComIndex = 0;
    mol.complexId = 0;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.isEmpty = false;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = com;
    iface.index = 0;
    iface.relIndex = 0;
    iface.stateIndex = 0;
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    return mol;
}

// -----------------------------------------------------------------------------
// Helper: build a Complex made of the single molecule stored at index 0 of the
// moleculeList.  Diffusion constants are supplied isotropically.
// -----------------------------------------------------------------------------
Complex ccpv_make_complex(const Coord& com, double D, double Dr)
{
    Complex targCom;
    targCom.index = 0;
    targCom.id = 0;
    targCom.isEmpty = false;
    targCom.comCoord = com;
    targCom.mass = 1.0;
    targCom.radius = 1.0;

    targCom.D = Coord { D, D, D };
    targCom.Dr = Coord { Dr, Dr, Dr };

    targCom.memberList.clear();
    targCom.memberList.push_back(0);
    targCom.numEachMol.assign(1, 1);
    targCom.lastNumberUpdateItrEachMol.assign(1, 0);

    targCom.OnSurface = false;
    targCom.onFiber = false;
    targCom.linksToSurface = 0;

    // start with zeroed propagation vectors
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot.zero_crds();

    return targCom;
}

// -----------------------------------------------------------------------------
// Helper: a cubic, reflecting water box centred on the origin.
// -----------------------------------------------------------------------------
Membrane ccpv_make_box_membrane(double side)
{
    Membrane membraneObject;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { side, side, side });
    membraneObject.isBox = true;
    membraneObject.isSphere = false;
    membraneObject.implicitLipid = false;
    membraneObject.hasCompartment = false;
    membraneObject.xBCtype = "reflect";
    membraneObject.yBCtype = "reflect";
    membraneObject.zBCtype = "reflect";
    return membraneObject;
}

// -----------------------------------------------------------------------------
// Helper: one MolTemplate matching the molecules created above.  Only
// MolTemplate::insideCompartment is read by the function under test, but every
// field that could be dereferenced is initialised anyway.
// -----------------------------------------------------------------------------
std::vector<MolTemplate> ccpv_make_templates()
{
    std::vector<MolTemplate> molTemplateList;
    MolTemplate tmp;
    tmp.molTypeIndex = 0;
    tmp.molName = "A";
    tmp.copies = 1;
    tmp.mass = 1.0;
    tmp.radius = 1.0;
    tmp.D = Coord { 10.0, 10.0, 10.0 };
    tmp.Dr = Coord { 0.01, 0.01, 0.01 };
    tmp.isLipid = false;
    tmp.isImplicitLipid = false;
    tmp.isPoint = true;
    tmp.insideCompartment = false;
    tmp.outsideCompartment = false;
    tmp.crossesCompartment = false;
    tmp.interfaceList.emplace_back("a", Coord { 0.0, 0.0, 0.0 });
    tmp.interfaceList.back().index = 0;
    molTemplateList.push_back(tmp);
    return molTemplateList;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the exact formula and the order in which random numbers are consumed.
//
// Criteria: with the identical seed, the six components must be exactly
//           sqrt(2*dt*D_i) * g_i where g_0..g_5 are the first six GaussV()
//           draws (order: x, y, z, rot x, rot y, rot z).
// -----------------------------------------------------------------------------
void test_ccpv_exact_gaussian_formula_and_draw_order()
{
    std::cerr << "\n[TEST] test_ccpv_exact_gaussian_formula_and_draw_order\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Function:    create_complex_propagation_vectors (box branch)\n"
              << "  Checking:    trajTrans/trajRot == sqrt(2*dt*D)*GaussV() and\n"
              << "               that draws are consumed in the order x,y,z,rx,ry,rz.\n";

    Parameters params;
    params.timeStep = 0.1; // us
    const double D = 10.0; // nm^2/us
    const double Dr = 0.05; // rad^2/us

    // Reference: six Gaussian deviates from the freshly seeded generator.
    ccpv_seed_rng();
    double g[6];
    for (int i = 0; i < 6; ++i)
        g[i] = GaussV();
    std::cerr << "  Reference GaussV() draws: " << g[0] << ", " << g[1] << ", " << g[2] << ", "
              << g[3] << ", " << g[4] << ", " << g[5] << '\n';

    // Now re-seed and let the function under test do the drawing.
    ccpv_seed_rng();

    Membrane membraneObject = ccpv_make_box_membrane(1000.0); // huge box => no reflection
    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();
    std::vector<Molecule> moleculeList { ccpv_make_molecule(Coord { 0.0, 0.0, 0.0 }, 0) };
    Complex targCom = ccpv_make_complex(Coord { 0.0, 0.0, 0.0 }, D, Dr);
    std::vector<Complex> complexList { targCom };

    create_complex_propagation_vectors(
        params, targCom, moleculeList, complexList, molTemplateList, membraneObject);

    const double transScale = std::sqrt(2.0 * params.timeStep * D);
    const double rotScale = std::sqrt(2.0 * params.timeStep * Dr);

    std::cerr << "  Produced trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";
    std::cerr << "  Produced trajRot   = (" << targCom.trajRot.x << ", " << targCom.trajRot.y
              << ", " << targCom.trajRot.z << ")\n";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, transScale * g[0])
        << "trajTrans.x must be sqrt(2*dt*D.x) times the 1st Gaussian draw";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, transScale * g[1])
        << "trajTrans.y must be sqrt(2*dt*D.y) times the 2nd Gaussian draw";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, transScale * g[2])
        << "trajTrans.z must be sqrt(2*dt*D.z) times the 3rd Gaussian draw";
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, rotScale * g[3])
        << "trajRot.x must be sqrt(2*dt*Dr.x) times the 4th Gaussian draw";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, rotScale * g[4])
        << "trajRot.y must be sqrt(2*dt*Dr.y) times the 5th Gaussian draw";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, rotScale * g[5])
        << "trajRot.z must be sqrt(2*dt*Dr.z) times the 6th Gaussian draw";
}

// -----------------------------------------------------------------------------
// Test 2: all diffusion constants zero.
//
// Criteria: every component is exactly 0.0, and no random number is consumed
//           (the next GaussV() must equal the very first draw of the seeded
//           sequence).
// -----------------------------------------------------------------------------
void test_ccpv_zero_diffusion_is_exactly_zero_and_draws_nothing()
{
    std::cerr << "\n[TEST] test_ccpv_zero_diffusion_is_exactly_zero_and_draws_nothing\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Checking:    D == 0 and Dr == 0 produce exactly zero motion AND\n"
              << "               skip the GaussV() calls entirely.\n";

    Parameters params;
    params.timeStep = 0.1;

    // What the first draw of the seeded sequence looks like.
    ccpv_seed_rng();
    const double firstDraw = GaussV();
    std::cerr << "  First GaussV() of the seeded stream = " << firstDraw << '\n';

    ccpv_seed_rng();

    Membrane membraneObject = ccpv_make_box_membrane(1000.0);
    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();
    std::vector<Molecule> moleculeList { ccpv_make_molecule(Coord { 0.0, 0.0, 0.0 }, 0) };
    Complex targCom = ccpv_make_complex(Coord { 0.0, 0.0, 0.0 }, 0.0, 0.0);
    std::vector<Complex> complexList { targCom };

    create_complex_propagation_vectors(
        params, targCom, moleculeList, complexList, molTemplateList, membraneObject);

    std::cerr << "  Produced trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";
    std::cerr << "  Produced trajRot   = (" << targCom.trajRot.x << ", " << targCom.trajRot.y
              << ", " << targCom.trajRot.z << ")\n";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, 0.0) << "frozen x translation must be exactly 0";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "frozen y translation must be exactly 0";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, 0.0) << "frozen z translation must be exactly 0";
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, 0.0) << "frozen x rotation must be exactly 0";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, 0.0) << "frozen y rotation must be exactly 0";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, 0.0) << "frozen z rotation must be exactly 0";

    // If no random numbers were consumed, the generator is still at the start of
    // the stream and the next draw must reproduce firstDraw.
    const double afterCall = GaussV();
    std::cerr << "  GaussV() after the call             = " << afterCall
              << " (expected " << firstDraw << ")\n";
    EXPECT_DOUBLE_EQ(afterCall, firstDraw)
        << "no random deviate should be consumed when D and Dr are all zero";
}

// -----------------------------------------------------------------------------
// Test 3: mixed zero / non-zero diffusion constants.
//
// Criteria: only the non-zero components consume a draw, and they do so in
//           order.  Here D = (Dx, 0, Dz), Dr = (0, Dry, 0), so the expected
//           draw order is: Dx -> g0, Dz -> g1, Dry -> g2.
// -----------------------------------------------------------------------------
void test_ccpv_mixed_zero_diffusion_skips_draws()
{
    std::cerr << "\n[TEST] test_ccpv_mixed_zero_diffusion_skips_draws\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Checking:    components with D == 0 are skipped so the remaining\n"
              << "               components use consecutive Gaussian deviates.\n";

    Parameters params;
    params.timeStep = 0.25;
    const double Dx = 4.0;
    const double Dz = 9.0;
    const double Dry = 0.16;

    // Reference draws for the three "live" degrees of freedom.
    ccpv_seed_rng();
    const double g0 = GaussV();
    const double g1 = GaussV();
    const double g2 = GaussV();
    std::cerr << "  Reference draws: " << g0 << ", " << g1 << ", " << g2 << '\n';

    ccpv_seed_rng();

    Membrane membraneObject = ccpv_make_box_membrane(1000.0);
    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();
    std::vector<Molecule> moleculeList { ccpv_make_molecule(Coord { 0.0, 0.0, 0.0 }, 0) };
    Complex targCom = ccpv_make_complex(Coord { 0.0, 0.0, 0.0 }, 0.0, 0.0);
    targCom.D = Coord { Dx, 0.0, Dz };
    targCom.Dr = Coord { 0.0, Dry, 0.0 };
    std::vector<Complex> complexList { targCom };

    create_complex_propagation_vectors(
        params, targCom, moleculeList, complexList, molTemplateList, membraneObject);

    std::cerr << "  Produced trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";
    std::cerr << "  Produced trajRot   = (" << targCom.trajRot.x << ", " << targCom.trajRot.y
              << ", " << targCom.trajRot.z << ")\n";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, std::sqrt(2.0 * params.timeStep * Dx) * g0)
        << "x translation should use the 1st deviate";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, 0.0) << "y translation is frozen (D.y == 0)";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, std::sqrt(2.0 * params.timeStep * Dz) * g1)
        << "z translation should use the 2nd deviate (y consumed none)";
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, 0.0) << "x rotation is frozen (Dr.x == 0)";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, std::sqrt(2.0 * params.timeStep * Dry) * g2)
        << "y rotation should use the 3rd deviate";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, 0.0) << "z rotation is frozen (Dr.z == 0)";
}

// -----------------------------------------------------------------------------
// Test 4: statistics of the generated displacements.
//
// Criteria: over many samples in a very large box (so reflection essentially
//           never fires) the sample mean is ~0 and the sample standard
//           deviation is ~sqrt(2*dt*D) to within a few percent.
// -----------------------------------------------------------------------------
void test_ccpv_displacement_statistics()
{
    std::cerr << "\n[TEST] test_ccpv_displacement_statistics\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Checking:    <dx> ~ 0 and std(dx) ~ sqrt(2*dt*D) over many samples.\n";

    ccpv_seed_rng();

    Parameters params;
    params.timeStep = 0.1;
    const double D = 10.0;
    const double expectedSigma = std::sqrt(2.0 * params.timeStep * D);

    Membrane membraneObject = ccpv_make_box_membrane(10000.0); // effectively unbounded
    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();
    std::vector<Molecule> moleculeList { ccpv_make_molecule(Coord { 0.0, 0.0, 0.0 }, 0) };
    Complex targCom = ccpv_make_complex(Coord { 0.0, 0.0, 0.0 }, D, 0.0);
    std::vector<Complex> complexList { targCom };

    const int nSamples = 5000;
    double sum[3] = { 0.0, 0.0, 0.0 };
    double sumSq[3] = { 0.0, 0.0, 0.0 };

    for (int i = 0; i < nSamples; ++i) {
        create_complex_propagation_vectors(
            params, targCom, moleculeList, complexList, molTemplateList, membraneObject);
        const double v[3] = { targCom.trajTrans.x, targCom.trajTrans.y, targCom.trajTrans.z };
        for (int k = 0; k < 3; ++k) {
            sum[k] += v[k];
            sumSq[k] += v[k] * v[k];
        }
    }

    const char* label[3] = { "x", "y", "z" };
    for (int k = 0; k < 3; ++k) {
        const double mean = sum[k] / nSamples;
        const double var = (sumSq[k] / nSamples) - (mean * mean);
        const double sigma = std::sqrt(var);
        std::cerr << "  component " << label[k] << ": mean = " << mean << " (expect ~0), sigma = "
                  << sigma << " (expect " << expectedSigma << ")\n";

        // Mean: standard error is sigma/sqrt(N) ~ 0.02, allow 10% of sigma.
        EXPECT_NEAR(mean, 0.0, 0.10 * expectedSigma)
            << "sample mean of the " << label[k] << " displacement should be ~0";
        // Std dev: relative error is ~1/sqrt(2N) ~ 1%, allow 8%.
        EXPECT_NEAR(sigma, expectedSigma, 0.08 * expectedSigma)
            << "sample sigma of the " << label[k] << " displacement should be sqrt(2*dt*D)";
    }
}

// -----------------------------------------------------------------------------
// Test 5: reflecting boundary conditions are applied after the draw.
//
// Criteria: a complex parked next to the +x wall, with a step size comparable
//           to its distance from the wall, must never end up outside the box.
// -----------------------------------------------------------------------------
void test_ccpv_reflects_at_box_wall()
{
    std::cerr << "\n[TEST] test_ccpv_reflects_at_box_wall\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Checking:    reflect_traj_complex_rad_rot is applied, so\n"
              << "               comCoord + trajTrans always stays inside the box.\n";

    ccpv_seed_rng();

    Parameters params;
    params.timeStep = 1.0;
    const double D = 50.0; // sigma = 10 nm, comparable to the 5 nm gap to the wall

    const double side = 100.0;
    const double halfBox = side / 2.0;
    Membrane membraneObject = ccpv_make_box_membrane(side);
    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();

    const Coord startCrd { 45.0, 0.0, 0.0 }; // 5 nm from the +x wall
    std::vector<Molecule> moleculeList { ccpv_make_molecule(startCrd, 0) };
    Complex targCom = ccpv_make_complex(startCrd, D, 0.0); // no rotation: keeps the geometry simple
    std::vector<Complex> complexList { targCom };

    const int nTrials = 200;
    double maxFinalX = -1e300;
    double minFinalX = 1e300;
    int nBeyondWallBeforeCheck = 0;

    for (int i = 0; i < nTrials; ++i) {
        // The function never moves the complex, only sets its trajectory, so the
        // starting coordinate stays put between calls.
        create_complex_propagation_vectors(
            params, targCom, moleculeList, complexList, molTemplateList, membraneObject);

        const double finalX = targCom.comCoord.x + targCom.trajTrans.x;
        const double finalY = targCom.comCoord.y + targCom.trajTrans.y;
        const double finalZ = targCom.comCoord.z + targCom.trajTrans.z;

        if (finalX > maxFinalX)
            maxFinalX = finalX;
        if (finalX < minFinalX)
            minFinalX = finalX;
        if (std::abs(targCom.trajTrans.x) > 5.0)
            ++nBeyondWallBeforeCheck;

        EXPECT_LE(finalX, halfBox + 1e-6) << "trial " << i << ": complex escaped through +x wall";
        EXPECT_GE(finalX, -halfBox - 1e-6) << "trial " << i << ": complex escaped through -x wall";
        EXPECT_LE(finalY, halfBox + 1e-6) << "trial " << i << ": complex escaped through +y wall";
        EXPECT_GE(finalY, -halfBox - 1e-6) << "trial " << i << ": complex escaped through -y wall";
        EXPECT_LE(finalZ, halfBox + 1e-6) << "trial " << i << ": complex escaped through +z wall";
        EXPECT_GE(finalZ, -halfBox - 1e-6) << "trial " << i << ": complex escaped through -z wall";
    }

    std::cerr << "  " << nTrials << " trials, final x in [" << minFinalX << ", " << maxFinalX
              << "], wall at +/-" << halfBox << '\n';
    std::cerr << "  " << nBeyondWallBeforeCheck
              << " trials drew |dx| > 5 nm (i.e. would have crossed the wall without reflection)\n";

    // Sanity check on the test itself: the step size must actually challenge the wall.
    EXPECT_GT(nBeyondWallBeforeCheck, 0)
        << "test is not exercising the reflection code; increase D or move closer to the wall";
}

// -----------------------------------------------------------------------------
// Test 6: implicit-lipid RS3D look-up table.
//
// Criteria: with membraneObject.implicitLipid == true the RS3Dvect table is
//           scanned for the molecule type; the look-up must not crash and, for
//           a complex far away from the membrane, must not alter the result
//           relative to the non-implicit-lipid case with the same seed.
// -----------------------------------------------------------------------------
void test_ccpv_implicit_lipid_lookup()
{
    std::cerr << "\n[TEST] test_ccpv_implicit_lipid_lookup\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Checking:    the RS3Dvect look-up (indices 300..499) is safe and\n"
              << "               leaves an interior complex's trajectory untouched.\n";

    Parameters params;
    params.timeStep = 0.1;
    const double D = 10.0;
    const double Dr = 0.05;

    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();
    const Coord startCrd { 0.0, 0.0, 0.0 }; // dead centre of a large box

    // --- reference run without implicit lipid -------------------------------
    ccpv_seed_rng();
    Membrane plainMembrane = ccpv_make_box_membrane(1000.0);
    std::vector<Molecule> molsA { ccpv_make_molecule(startCrd, 0) };
    Complex comA = ccpv_make_complex(startCrd, D, Dr);
    std::vector<Complex> listA { comA };
    create_complex_propagation_vectors(params, comA, molsA, listA, molTemplateList, plainMembrane);

    // --- same seed, implicit lipid table populated --------------------------
    ccpv_seed_rng();
    Membrane ilMembrane = ccpv_make_box_membrane(1000.0);
    ilMembrane.implicitLipid = true;
    // The table is read at [i + 400] (molecule type) and [i + 300] (RS3D value)
    // for i in [0,100), so at least 500 entries must exist.
    ilMembrane.RS3Dvect.assign(500, -1.0);
    ilMembrane.RS3Dvect[400] = 0.0; // molTypeIndex 0 lives in the first slot
    ilMembrane.RS3Dvect[300] = 0.75; // its reflecting-surface offset

    std::vector<Molecule> molsB { ccpv_make_molecule(startCrd, 0) };
    Complex comB = ccpv_make_complex(startCrd, D, Dr);
    std::vector<Complex> listB { comB };
    create_complex_propagation_vectors(params, comB, molsB, listB, molTemplateList, ilMembrane);

    std::cerr << "  no-IL  trajTrans = (" << comA.trajTrans.x << ", " << comA.trajTrans.y << ", "
              << comA.trajTrans.z << ")\n";
    std::cerr << "  with-IL trajTrans = (" << comB.trajTrans.x << ", " << comB.trajTrans.y << ", "
              << comB.trajTrans.z << ")\n";

    EXPECT_DOUBLE_EQ(comB.trajTrans.x, comA.trajTrans.x)
        << "RS3D look-up must not change x motion of an interior complex";
    EXPECT_DOUBLE_EQ(comB.trajTrans.y, comA.trajTrans.y)
        << "RS3D look-up must not change y motion of an interior complex";
    EXPECT_DOUBLE_EQ(comB.trajTrans.z, comA.trajTrans.z)
        << "RS3D look-up must not change z motion of an interior complex";
    EXPECT_DOUBLE_EQ(comB.trajRot.x, comA.trajRot.x) << "rotation x must be identical";
    EXPECT_DOUBLE_EQ(comB.trajRot.y, comA.trajRot.y) << "rotation y must be identical";
    EXPECT_DOUBLE_EQ(comB.trajRot.z, comA.trajRot.z) << "rotation z must be identical";
}

// -----------------------------------------------------------------------------
// Test 7: a complex inside a spherical volume but NOT on its surface still uses
//         the plain Gaussian branch.
//
// Criteria: OnSurface == false with isSphere == true must reproduce exactly the
//           same sqrt(2*dt*D)*GaussV() values as the box case.
// -----------------------------------------------------------------------------
void test_ccpv_inside_sphere_uses_cartesian_branch()
{
    std::cerr << "\n[TEST] test_ccpv_inside_sphere_uses_cartesian_branch\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Checking:    isSphere == true but OnSurface == false takes the\n"
              << "               ordinary Cartesian Gaussian branch.\n";

    Parameters params;
    params.timeStep = 0.1;
    const double D = 5.0;
    const double Dr = 0.02;

    // reference draws
    ccpv_seed_rng();
    double g[6];
    for (int i = 0; i < 6; ++i)
        g[i] = GaussV();

    ccpv_seed_rng();

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 500.0; // large sphere so no reflection occurs
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 1000.0, 1000.0, 1000.0 });
    membraneObject.implicitLipid = false;
    membraneObject.hasCompartment = false;

    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();
    std::vector<Molecule> moleculeList { ccpv_make_molecule(Coord { 0.0, 0.0, 0.0 }, 0) };
    Complex targCom = ccpv_make_complex(Coord { 0.0, 0.0, 0.0 }, D, Dr);
    targCom.OnSurface = false; // freely diffusing inside the sphere
    std::vector<Complex> complexList { targCom };

    create_complex_propagation_vectors(
        params, targCom, moleculeList, complexList, molTemplateList, membraneObject);

    const double transScale = std::sqrt(2.0 * params.timeStep * D);
    const double rotScale = std::sqrt(2.0 * params.timeStep * Dr);

    std::cerr << "  Produced trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, transScale * g[0]) << "x must use the Cartesian formula";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, transScale * g[1]) << "y must use the Cartesian formula";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, transScale * g[2]) << "z must use the Cartesian formula";
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, rotScale * g[3]) << "rot x must use the Cartesian formula";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, rotScale * g[4]) << "rot y must use the Cartesian formula";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, rotScale * g[5]) << "rot z must use the Cartesian formula";
}

// -----------------------------------------------------------------------------
// Test 8: the "on the sphere surface" branch.
//
// Criteria: with OnSurface == true and isSphere == true the surface propagator
//           is used.  We cannot predict the exact deviates consumed by that
//           helper, so we assert that (a) the translation is finite and
//           non-degenerate, and (b) with Dr == 0 the rotation is exactly zero.
// -----------------------------------------------------------------------------
void test_ccpv_on_sphere_surface_branch()
{
    std::cerr << "\n[TEST] test_ccpv_on_sphere_surface_branch\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Checking:    OnSurface == true && isSphere == true delegates to\n"
              << "               create_complex_propagation_vectors_on_sphere and\n"
              << "               yields finite motion (zero rotation when Dr == 0).\n";

    ccpv_seed_rng();

    Parameters params;
    params.timeStep = 0.1;
    const double D = 1.0;

    const double R = 100.0;
    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = R;
    membraneObject.sphereVol = (4.0 / 3.0) * M_PI * R * R * R;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 2 * R, 2 * R, 2 * R });
    membraneObject.implicitLipid = false;
    membraneObject.hasCompartment = false;

    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();

    // The complex sits exactly on the spherical shell (required: the surface
    // propagator works in spherical coordinates, so the radius must be > 0).
    const Coord onShell { 0.0, 0.0, -R };
    std::vector<Molecule> moleculeList { ccpv_make_molecule(onShell, 0) };
    Complex targCom = ccpv_make_complex(onShell, D, 0.0); // Dr == 0 on purpose
    targCom.OnSurface = true;
    targCom.D.z = 0.0; // membrane-bound complexes have no motion normal to the surface
    std::vector<Complex> complexList { targCom };

    create_complex_propagation_vectors(
        params, targCom, moleculeList, complexList, molTemplateList, membraneObject);

    std::cerr << "  Produced trajTrans = (" << targCom.trajTrans.x << ", " << targCom.trajTrans.y
              << ", " << targCom.trajTrans.z << ")\n";
    std::cerr << "  Produced trajRot   = (" << targCom.trajRot.x << ", " << targCom.trajRot.y
              << ", " << targCom.trajRot.z << ")\n";

    EXPECT_TRUE(std::isfinite(targCom.trajTrans.x)) << "surface x displacement must be finite";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.y)) << "surface y displacement must be finite";
    EXPECT_TRUE(std::isfinite(targCom.trajTrans.z)) << "surface z displacement must be finite";

    // Dr == 0 means sqrt(2*dt*0)*GaussV() == 0 for every rotational component.
    EXPECT_DOUBLE_EQ(targCom.trajRot.x, 0.0) << "rotation must vanish when Dr.x == 0";
    EXPECT_DOUBLE_EQ(targCom.trajRot.y, 0.0) << "rotation must vanish when Dr.y == 0";
    EXPECT_DOUBLE_EQ(targCom.trajRot.z, 0.0) << "rotation must vanish when Dr.z == 0";

    // The step must be a genuine displacement, not identically zero.
    const double stepMag = std::sqrt(targCom.trajTrans.x * targCom.trajTrans.x
        + targCom.trajTrans.y * targCom.trajTrans.y + targCom.trajTrans.z * targCom.trajTrans.z);
    std::cerr << "  |trajTrans| on the sphere = " << stepMag << '\n';
    EXPECT_GT(stepMag, 0.0) << "a diffusing surface complex should actually move";

    // The complex must remain inside/on the spherical boundary after reflection.
    const double fx = targCom.comCoord.x + targCom.trajTrans.x;
    const double fy = targCom.comCoord.y + targCom.trajTrans.y;
    const double fz = targCom.comCoord.z + targCom.trajTrans.z;
    const double radial = std::sqrt(fx * fx + fy * fy + fz * fz);
    std::cerr << "  final radial distance = " << radial << " (sphere R = " << R << ")\n";
    EXPECT_LE(radial, R + 1e-6) << "surface complex must not leave the sphere";
}

// -----------------------------------------------------------------------------
// Test 9: repeated calls keep producing fresh, different vectors, and the
//         function is fully reproducible for a fixed seed.
// -----------------------------------------------------------------------------
void test_ccpv_repeatability_and_variation()
{
    std::cerr << "\n[TEST] test_ccpv_repeatability_and_variation\n"
              << "  Source file: create_complex_propagation_vectors.cpp\n"
              << "  Checking:    consecutive calls give different vectors, while the\n"
              << "               same seed reproduces the same sequence exactly.\n";

    Parameters params;
    params.timeStep = 0.1;
    const double D = 10.0;
    const double Dr = 0.05;

    Membrane membraneObject = ccpv_make_box_membrane(1000.0);
    std::vector<MolTemplate> molTemplateList = ccpv_make_templates();

    // --- first pass ---------------------------------------------------------
    ccpv_seed_rng();
    std::vector<Molecule> mols1 { ccpv_make_molecule(Coord { 0.0, 0.0, 0.0 }, 0) };
    Complex com1 = ccpv_make_complex(Coord { 0.0, 0.0, 0.0 }, D, Dr);
    std::vector<Complex> list1 { com1 };

    create_complex_propagation_vectors(params, com1, mols1, list1, molTemplateList, membraneObject);
    const double firstX = com1.trajTrans.x;
    const double firstRotX = com1.trajRot.x;

    create_complex_propagation_vectors(params, com1, mols1, list1, molTemplateList, membraneObject);
    const double secondX = com1.trajTrans.x;

    std::cerr << "  call 1 dx = " << firstX << ", call 2 dx = " << secondX << '\n';
    EXPECT_NE(firstX, secondX)
        << "two consecutive calls should draw different random displacements";

    // --- second pass with the same seed ------------------------------------
    ccpv_seed_rng();
    std::vector<Molecule> mols2 { ccpv_make_molecule(Coord { 0.0, 0.0, 0.0 }, 0) };
    Complex com2 = ccpv_make_complex(Coord { 0.0, 0.0, 0.0 }, D, Dr);
    std::vector<Complex> list2 { com2 };

    create_complex_propagation_vectors(params, com2, mols2, list2, molTemplateList, membraneObject);

    std::cerr << "  reseeded call 1 dx = " << com2.trajTrans.x << " (expected " << firstX << ")\n";
    EXPECT_DOUBLE_EQ(com2.trajTrans.x, firstX)
        << "the same seed must reproduce the same translation";
    EXPECT_DOUBLE_EQ(com2.trajRot.x, firstRotX)
        << "the same seed must reproduce the same rotation";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* routine is executed inside its own
// TEST so a failure in one does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(CreateComplexPropagationVectors, ExactGaussianFormulaAndDrawOrder)
{
    test_ccpv_exact_gaussian_formula_and_draw_order();
}
TEST(CreateComplexPropagationVectors, ZeroDiffusionIsExactlyZeroAndDrawsNothing)
{
    test_ccpv_zero_diffusion_is_exactly_zero_and_draws_nothing();
}
TEST(CreateComplexPropagationVectors, MixedZeroDiffusionSkipsDraws)
{
    test_ccpv_mixed_zero_diffusion_skips_draws();
}
TEST(CreateComplexPropagationVectors, DisplacementStatistics)
{
    test_ccpv_displacement_statistics();
}
TEST(CreateComplexPropagationVectors, ReflectsAtBoxWall) { test_ccpv_reflects_at_box_wall(); }
TEST(CreateComplexPropagationVectors, ImplicitLipidLookup) { test_ccpv_implicit_lipid_lookup(); }
TEST(CreateComplexPropagationVectors, InsideSphereUsesCartesianBranch)
{
    test_ccpv_inside_sphere_uses_cartesian_branch();
}
TEST(CreateComplexPropagationVectors, OnSphereSurfaceBranch)
{
    test_ccpv_on_sphere_surface_branch();
}
TEST(CreateComplexPropagationVectors, RepeatabilityAndVariation)
{
    test_ccpv_repeatability_and_variation();
}