/*! \file test_create_complex_propagation_vectors_on_sphere.cpp
 *
 * ### Unit test for
 *     src/trajectory_functions/create_complex_propagation_vectors_on_sphere.cpp
 *
 * Function under test:
 *
 *     Coord create_complex_propagation_vectors_on_sphere(const Parameters& params,
 *                                                        Complex& targCom);
 *
 * What the function does (read from the implementation):
 *   1. Projects the complex centre of mass onto a sphere of radius
 *      `params.sphereR` (so only the *direction* of `targCom.comCoord` matters).
 *   2. Draws two independent Gaussian displacements
 *          dx = sqrt(2*dt*D.x) * GaussV()
 *          dy = sqrt(2*dt*D.y) * GaussV()
 *      and builds a geodesic step of length dl = sqrt(dx^2 + dy^2) in a
 *      uniformly random direction within the local tangent plane.
 *   3. Moves the projected COM along that geodesic (angle dtheta = dl / R) and
 *      returns the *chord* vector `COMnew - COMprojected`.
 *
 * Consequences that make good, deterministic assertions:
 *   - The end point `COMprojected + trajTrans` always lies exactly on the
 *     sphere of radius R (rotations/geodesic moves conserve the radius).
 *   - The chord always points slightly *inward*, i.e.
 *     dot(unit(COM), trajTrans) = R*(cos(dtheta) - 1) <= 0.
 *   - The result depends only on the *direction* of comCoord, never on its
 *     magnitude (because of the projection in step 1).
 *   - With the same RNG seed the function is bit-for-bit reproducible.
 *   - Statistically E[dl^2] = 2*dt*(D.x + D.y), and the step direction is
 *     isotropic in the tangent plane so E[trajTrans] = 0.
 *
 * DEGENERATE INPUTS THAT ARE DELIBERATELY *NOT* EXERCISED (they produce NaN /
 * divide-by-zero and would poison the whole suite):
 *   - `targCom.comCoord == (0,0,0)`  -> COM_norm == 0 -> division by zero.
 *   - `params.sphereR == 0`          -> dtheta = dl / 0.
 *   - `D.x == D.y == 0`              -> dl == 0 -> acos(dx/dl) == acos(NaN).
 */

#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// `r` (the GSL random number generator) is defined in gtest_main.cpp; we only
// declare/refer to it here. Never define it in this translation unit.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Helper: (re)seed the global GSL generator in a reproducible way.
//
// NOTE: srand_gsl() must not be used. We allocate the generator once (if the
// suite has not already done so) and then simply re-seed it for every test so
// that each test starts from a known, identical random stream.
// -----------------------------------------------------------------------------
void ccpvs_seed_rng(unsigned long seed)
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, seed);
}

//! \brief Euclidean length of a Coord.
double ccpvs_norm(const Coord& c) { return std::sqrt(c.x * c.x + c.y * c.y + c.z * c.z); }

/*! \brief Build a minimal Complex sitting on (or along) a given direction.
 *
 * Only `comCoord` and the translational diffusion constants `D` are read by the
 * function under test, but we fill in everything that is cheap to initialise so
 * the object is never left half-constructed.
 *
 * \param[in] dir  Direction of the centre of mass (need not be normalised).
 * \param[in] len  Length to scale `dir` to (the function projects onto sphereR
 *                 anyway, so this only matters for the projection test).
 * \param[in] Dx   Translational diffusion constant in x.
 * \param[in] Dy   Translational diffusion constant in y.
 */
Complex ccpvs_make_complex(const Coord& dir, double len, double Dx, double Dy)
{
    const double mag = ccpvs_norm(dir);
    Complex targCom;
    targCom.comCoord = Coord { len * dir.x / mag, len * dir.y / mag, len * dir.z / mag };

    // Translational diffusion: only x and y are used by the routine, but z is
    // set explicitly so nothing is read uninitialised.
    targCom.D = Coord { Dx, Dy, 0.0 };
    targCom.Dr = Coord { 0.0, 0.0, 0.0 };

    // Housekeeping members (unused by the function, but keep the object sane).
    targCom.index = 0;
    targCom.radius = 1.0;
    targCom.mass = 1.0;
    targCom.memberList.clear();
    targCom.numEachMol.clear();
    targCom.OnSurface = true; // this routine is only used for surface complexes
    targCom.trajTrans = Vector { 0.0, 0.0, 0.0 };
    targCom.trajRot = Coord { 0.0, 0.0, 0.0 };

    return targCom;
}

/*! \brief Build a Parameters object with the only two fields the routine reads. */
Parameters ccpvs_make_params(double sphereR, double timeStep)
{
    Parameters params;
    params.sphereR = sphereR;
    params.timeStep = timeStep;
    return params;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: the propagated point never leaves the sphere.
//
// Criterion: |COMprojected + trajTrans| == params.sphereR for every sample.
// -----------------------------------------------------------------------------
void ccpvs_test_stays_on_sphere()
{
    std::cerr << "\n[TEST] ccpvs_test_stays_on_sphere\n"
              << "  Source file: create_complex_propagation_vectors_on_sphere.cpp\n"
              << "  Function:    create_complex_propagation_vectors_on_sphere()\n"
              << "  Scenario:    500 random propagation steps of a surface complex.\n"
              << "  Criterion:   the end point stays on the sphere (|COM+traj| == R)\n"
              << "               and the chord points inward (dot(COMhat,traj) <= 0).\n";

    const double R = 100.0;
    Parameters params = ccpvs_make_params(R, 1.0);

    // Direction (1,1,1) - deliberately away from the poles so that the
    // spherical-coordinate conversion is well behaved.
    const Coord dir { 1.0, 1.0, 1.0 };

    ccpvs_seed_rng(42);

    int worstIdx = -1;
    double worstRadialErr = 0.0;
    double worstInwardDot = 0.0;

    for (int i = 0; i < 500; ++i) {
        Complex targCom = ccpvs_make_complex(dir, R, 1.0, 1.0);

        // COM as the routine sees it: projected onto the sphere of radius R.
        const Coord comProj = targCom.comCoord;

        const Coord traj = create_complex_propagation_vectors_on_sphere(params, targCom);

        // Result must be a finite vector.
        ASSERT_TRUE(std::isfinite(traj.x) && std::isfinite(traj.y) && std::isfinite(traj.z))
            << "Sample " << i << " produced a non-finite trajectory vector";

        const Coord endPt { comProj.x + traj.x, comProj.y + traj.y, comProj.z + traj.z };
        const double endRad = ccpvs_norm(endPt);

        if (std::fabs(endRad - R) > worstRadialErr) {
            worstRadialErr = std::fabs(endRad - R);
            worstIdx = i;
        }

        // The chord of a geodesic step always cuts inside the sphere, so the
        // radial component of trajTrans is <= 0 (it is R*(cos(dtheta)-1)).
        const double radialDot = (comProj.x * traj.x + comProj.y * traj.y + comProj.z * traj.z) / R;
        if (radialDot > worstInwardDot)
            worstInwardDot = radialDot;

        EXPECT_NEAR(endRad, R, 1e-6)
            << "Sample " << i << ": propagated point left the sphere (radius " << endRad << ")";
        EXPECT_LE(radialDot, 1e-9)
            << "Sample " << i << ": chord should not point outward (radial dot = " << radialDot << ")";
    }

    std::cerr << "  Worst radial deviation from R: " << worstRadialErr
              << " (sample " << worstIdx << ")\n"
              << "  Largest (least negative) radial dot: " << worstInwardDot << "\n";
}

// -----------------------------------------------------------------------------
// Test 2: only the *direction* of comCoord matters.
//
// The routine normalises comCoord onto the sphere before doing anything else,
// so two complexes lying on the same ray must produce the identical step for
// the same RNG stream.
// -----------------------------------------------------------------------------
void ccpvs_test_projection_independent_of_com_magnitude()
{
    std::cerr << "\n[TEST] ccpvs_test_projection_independent_of_com_magnitude\n"
              << "  Function:  create_complex_propagation_vectors_on_sphere()\n"
              << "  Scenario:  same direction, |comCoord| = R vs. |comCoord| = R/3\n"
              << "             vs. |comCoord| = 5R, identical RNG seed.\n"
              << "  Criterion: the three returned vectors agree (the COM is\n"
              << "             projected onto the sphere before use).\n";

    const double R = 50.0;
    Parameters params = ccpvs_make_params(R, 2.0);
    const Coord dir { 0.3, -0.7, 1.2 };

    // Case A: exactly on the sphere.
    ccpvs_seed_rng(42);
    Complex comA = ccpvs_make_complex(dir, R, 1.5, 1.5);
    const Coord trajA = create_complex_propagation_vectors_on_sphere(params, comA);

    // Case B: well inside the sphere, same ray.
    ccpvs_seed_rng(42);
    Complex comB = ccpvs_make_complex(dir, R / 3.0, 1.5, 1.5);
    const Coord trajB = create_complex_propagation_vectors_on_sphere(params, comB);

    // Case C: well outside the sphere, same ray.
    ccpvs_seed_rng(42);
    Complex comC = ccpvs_make_complex(dir, 5.0 * R, 1.5, 1.5);
    const Coord trajC = create_complex_propagation_vectors_on_sphere(params, comC);

    std::cerr << "  trajA = (" << trajA.x << ", " << trajA.y << ", " << trajA.z << ")\n"
              << "  trajB = (" << trajB.x << ", " << trajB.y << ", " << trajB.z << ")\n"
              << "  trajC = (" << trajC.x << ", " << trajC.y << ", " << trajC.z << ")\n";

    // Only tiny round-off differences are allowed (the projection divides by
    // slightly different magnitudes).
    EXPECT_NEAR(trajA.x, trajB.x, 1e-9) << "x differs between |COM|=R and |COM|=R/3";
    EXPECT_NEAR(trajA.y, trajB.y, 1e-9) << "y differs between |COM|=R and |COM|=R/3";
    EXPECT_NEAR(trajA.z, trajB.z, 1e-9) << "z differs between |COM|=R and |COM|=R/3";

    EXPECT_NEAR(trajA.x, trajC.x, 1e-9) << "x differs between |COM|=R and |COM|=5R";
    EXPECT_NEAR(trajA.y, trajC.y, 1e-9) << "y differs between |COM|=R and |COM|=5R";
    EXPECT_NEAR(trajA.z, trajC.z, 1e-9) << "z differs between |COM|=R and |COM|=5R";
}

// -----------------------------------------------------------------------------
// Test 3: reproducibility.
//
// Same seed -> same random draws -> bit-for-bit identical result.
// Different seed -> (with overwhelming probability) a different result.
// -----------------------------------------------------------------------------
void ccpvs_test_deterministic_with_same_seed()
{
    std::cerr << "\n[TEST] ccpvs_test_deterministic_with_same_seed\n"
              << "  Function:  create_complex_propagation_vectors_on_sphere()\n"
              << "  Criterion: identical seed -> identical vector;\n"
              << "             different seed -> different vector.\n";

    const double R = 80.0;
    Parameters params = ccpvs_make_params(R, 0.5);
    const Coord dir { 1.0, 0.0, 0.5 };

    ccpvs_seed_rng(42);
    Complex com1 = ccpvs_make_complex(dir, R, 1.0, 1.0);
    const Coord traj1 = create_complex_propagation_vectors_on_sphere(params, com1);

    ccpvs_seed_rng(42);
    Complex com2 = ccpvs_make_complex(dir, R, 1.0, 1.0);
    const Coord traj2 = create_complex_propagation_vectors_on_sphere(params, com2);

    std::cerr << "  seed 42 (1st) = (" << traj1.x << ", " << traj1.y << ", " << traj1.z << ")\n"
              << "  seed 42 (2nd) = (" << traj2.x << ", " << traj2.y << ", " << traj2.z << ")\n";

    EXPECT_DOUBLE_EQ(traj1.x, traj2.x) << "Same seed must reproduce x exactly";
    EXPECT_DOUBLE_EQ(traj1.y, traj2.y) << "Same seed must reproduce y exactly";
    EXPECT_DOUBLE_EQ(traj1.z, traj2.z) << "Same seed must reproduce z exactly";

    // A different seed should give a different displacement.
    ccpvs_seed_rng(12345);
    Complex com3 = ccpvs_make_complex(dir, R, 1.0, 1.0);
    const Coord traj3 = create_complex_propagation_vectors_on_sphere(params, com3);

    std::cerr << "  seed 12345    = (" << traj3.x << ", " << traj3.y << ", " << traj3.z << ")\n";

    const double diff = std::sqrt((traj3.x - traj1.x) * (traj3.x - traj1.x)
        + (traj3.y - traj1.y) * (traj3.y - traj1.y) + (traj3.z - traj1.z) * (traj3.z - traj1.z));
    EXPECT_GT(diff, 1e-12) << "A different RNG seed should give a different step";
}

// -----------------------------------------------------------------------------
// Test 4: the geodesic step length obeys the diffusion statistics.
//
// dl = sqrt(dx^2 + dy^2) with dx ~ N(0, 2*dt*D.x), dy ~ N(0, 2*dt*D.y), hence
//     E[dl^2] = 2*dt*(D.x + D.y).
// We recover dl from the returned chord:
//     chord   = 2 R sin(dtheta/2)  ->  dtheta = 2 asin(chord/(2R)),  dl = R*dtheta.
//
// The step direction is isotropic in the tangent plane, so the mean of each
// Cartesian component of trajTrans must be ~0.
// -----------------------------------------------------------------------------
void ccpvs_test_step_length_statistics()
{
    std::cerr << "\n[TEST] ccpvs_test_step_length_statistics\n"
              << "  Function:  create_complex_propagation_vectors_on_sphere()\n"
              << "  Scenario:  20000 samples, R=100, dt=1, D.x=D.y=1.\n"
              << "  Criterion: mean geodesic dl^2 ~ 2*dt*(D.x+D.y) = 4, and the\n"
              << "             mean displacement vector is ~0 (isotropic direction).\n";

    const double R = 100.0;
    const double dt = 1.0;
    const double Dx = 1.0;
    const double Dy = 1.0;
    const int nSamples = 20000;

    Parameters params = ccpvs_make_params(R, dt);
    const Coord dir { 1.0, 1.0, 1.0 };

    ccpvs_seed_rng(42);

    double sumSq = 0.0;
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    double maxChord = 0.0;

    for (int i = 0; i < nSamples; ++i) {
        Complex targCom = ccpvs_make_complex(dir, R, Dx, Dy);
        const Coord traj = create_complex_propagation_vectors_on_sphere(params, targCom);

        const double chord = ccpvs_norm(traj);
        if (chord > maxChord)
            maxChord = chord;

        // Recover the geodesic arc length from the chord.
        double s = chord / (2.0 * R);
        if (s > 1.0)
            s = 1.0; // guard against round-off just above 1
        const double dl = R * 2.0 * std::asin(s);

        sumSq += dl * dl;
        sumX += traj.x;
        sumY += traj.y;
        sumZ += traj.z;
    }

    const double meanSq = sumSq / nSamples;
    const double expectedMeanSq = 2.0 * dt * (Dx + Dy); // = 4.0
    const double meanX = sumX / nSamples;
    const double meanY = sumY / nSamples;
    const double meanZ = sumZ / nSamples;

    std::cerr << "  mean dl^2 = " << meanSq << " (expected " << expectedMeanSq << ")\n"
              << "  mean displacement = (" << meanX << ", " << meanY << ", " << meanZ << ")\n"
              << "  largest chord observed = " << maxChord << " (must be <= 2R = " << 2 * R << ")\n";

    // Standard error of the mean of dl^2 is ~0.03 here, so 0.3 is a very safe
    // (~10 sigma) tolerance that still catches a wrong scale factor.
    EXPECT_NEAR(meanSq, expectedMeanSq, 0.3)
        << "Mean squared geodesic step should be 2*dt*(D.x+D.y)";

    // Isotropic direction => each component averages to zero. The per-sample
    // component sd is ~1.4, so the SE of the mean is ~0.01; 0.15 is generous.
    EXPECT_NEAR(meanX, 0.0, 0.15) << "Mean x displacement should vanish (isotropic direction)";
    EXPECT_NEAR(meanY, 0.0, 0.15) << "Mean y displacement should vanish (isotropic direction)";
    EXPECT_NEAR(meanZ, 0.0, 0.15) << "Mean z displacement should vanish (isotropic direction)";

    // A chord on a sphere can never exceed the diameter.
    EXPECT_LE(maxChord, 2.0 * R + 1e-9) << "Chord length cannot exceed the sphere diameter";
}

// -----------------------------------------------------------------------------
// Test 5: anisotropic diffusion constants.
//
// The statistics must follow 2*dt*(D.x + D.y) even when D.x != D.y, because the
// two Gaussian draws are combined into a single radial step length.
// -----------------------------------------------------------------------------
void ccpvs_test_anisotropic_diffusion_statistics()
{
    std::cerr << "\n[TEST] ccpvs_test_anisotropic_diffusion_statistics\n"
              << "  Function:  create_complex_propagation_vectors_on_sphere()\n"
              << "  Scenario:  20000 samples with D.x = 2.0, D.y = 0.5, dt = 1.\n"
              << "  Criterion: mean dl^2 ~ 2*dt*(D.x+D.y) = 5, and every end\n"
              << "             point still lies on the sphere.\n";

    const double R = 120.0;
    const double dt = 1.0;
    const double Dx = 2.0;
    const double Dy = 0.5;
    const int nSamples = 20000;

    Parameters params = ccpvs_make_params(R, dt);
    const Coord dir { -0.4, 0.9, 0.2 };

    ccpvs_seed_rng(42);

    double sumSq = 0.0;
    int offSphere = 0;

    for (int i = 0; i < nSamples; ++i) {
        Complex targCom = ccpvs_make_complex(dir, R, Dx, Dy);
        const Coord comProj = targCom.comCoord;
        const Coord traj = create_complex_propagation_vectors_on_sphere(params, targCom);

        const Coord endPt { comProj.x + traj.x, comProj.y + traj.y, comProj.z + traj.z };
        if (std::fabs(ccpvs_norm(endPt) - R) > 1e-6)
            ++offSphere;

        double s = ccpvs_norm(traj) / (2.0 * R);
        if (s > 1.0)
            s = 1.0;
        const double dl = R * 2.0 * std::asin(s);
        sumSq += dl * dl;
    }

    const double meanSq = sumSq / nSamples;
    const double expectedMeanSq = 2.0 * dt * (Dx + Dy); // = 5.0

    std::cerr << "  mean dl^2 = " << meanSq << " (expected " << expectedMeanSq << ")\n"
              << "  samples that left the sphere = " << offSphere << " (expected 0)\n";

    EXPECT_NEAR(meanSq, expectedMeanSq, 0.4)
        << "Anisotropic D must still give mean dl^2 = 2*dt*(D.x+D.y)";
    EXPECT_EQ(offSphere, 0) << "No sample may leave the spherical surface";
}

// -----------------------------------------------------------------------------
// Test 6: scaling with the time step.
//
// dl^2 is linear in params.timeStep, so quadrupling dt must quadruple the mean
// squared step length.
// -----------------------------------------------------------------------------
void ccpvs_test_timestep_scaling()
{
    std::cerr << "\n[TEST] ccpvs_test_timestep_scaling\n"
              << "  Function:  create_complex_propagation_vectors_on_sphere()\n"
              << "  Scenario:  mean dl^2 measured for dt = 1 and dt = 4.\n"
              << "  Criterion: the ratio of the two means is ~4 (dl^2 ~ dt).\n";

    const double R = 200.0;
    const double Dx = 1.0;
    const double Dy = 1.0;
    const int nSamples = 20000;
    const Coord dir { 1.0, 0.2, -0.3 };

    // Small helper lambda measuring E[dl^2] for a given time step.
    auto meanSqStep = [&](double dt) -> double {
        Parameters params = ccpvs_make_params(R, dt);
        ccpvs_seed_rng(42); // identical random stream for both measurements
        double sumSq = 0.0;
        for (int i = 0; i < nSamples; ++i) {
            Complex targCom = ccpvs_make_complex(dir, R, Dx, Dy);
            const Coord traj = create_complex_propagation_vectors_on_sphere(params, targCom);
            double s = ccpvs_norm(traj) / (2.0 * R);
            if (s > 1.0)
                s = 1.0;
            const double dl = R * 2.0 * std::asin(s);
            sumSq += dl * dl;
        }
        return sumSq / nSamples;
    };

    const double m1 = meanSqStep(1.0);
    const double m4 = meanSqStep(4.0);
    const double ratio = m4 / m1;

    std::cerr << "  mean dl^2 (dt=1) = " << m1 << "\n"
              << "  mean dl^2 (dt=4) = " << m4 << "\n"
              << "  ratio            = " << ratio << " (expected 4)\n";

    // Because both runs use the same seed the ratio is essentially exact; a
    // loose tolerance still guards against an outright wrong dt dependence.
    EXPECT_NEAR(ratio, 4.0, 0.2) << "Mean squared step length must be proportional to timeStep";
    EXPECT_NEAR(m1, 4.0, 0.3) << "dt=1, D.x=D.y=1 should give E[dl^2] = 4";
}

// -----------------------------------------------------------------------------
// Test 7: the routine must not mutate the Complex it is handed.
//
// It takes `Complex&` but is documented/implemented as a pure computation: the
// caller is expected to assign the returned Coord into targCom.trajTrans.
// -----------------------------------------------------------------------------
void ccpvs_test_complex_is_not_modified()
{
    std::cerr << "\n[TEST] ccpvs_test_complex_is_not_modified\n"
              << "  Function:  create_complex_propagation_vectors_on_sphere()\n"
              << "  Criterion: comCoord, D and trajTrans of the input Complex are\n"
              << "             untouched (the step is only *returned*).\n";

    const double R = 60.0;
    Parameters params = ccpvs_make_params(R, 1.0);
    const Coord dir { 0.0, 1.0, 1.0 };

    ccpvs_seed_rng(42);
    Complex targCom = ccpvs_make_complex(dir, R, 1.0, 1.0);

    // Snapshot of everything the caller cares about.
    const Coord comBefore = targCom.comCoord;
    const Coord dBefore = targCom.D;
    const double transXBefore = targCom.trajTrans.x;
    const double transYBefore = targCom.trajTrans.y;
    const double transZBefore = targCom.trajTrans.z;

    const Coord traj = create_complex_propagation_vectors_on_sphere(params, targCom);

    std::cerr << "  returned traj = (" << traj.x << ", " << traj.y << ", " << traj.z << ")\n"
              << "  comCoord after = (" << targCom.comCoord.x << ", " << targCom.comCoord.y
              << ", " << targCom.comCoord.z << ")\n";

    EXPECT_DOUBLE_EQ(targCom.comCoord.x, comBefore.x) << "comCoord.x must not be modified";
    EXPECT_DOUBLE_EQ(targCom.comCoord.y, comBefore.y) << "comCoord.y must not be modified";
    EXPECT_DOUBLE_EQ(targCom.comCoord.z, comBefore.z) << "comCoord.z must not be modified";

    EXPECT_DOUBLE_EQ(targCom.D.x, dBefore.x) << "D.x must not be modified";
    EXPECT_DOUBLE_EQ(targCom.D.y, dBefore.y) << "D.y must not be modified";
    EXPECT_DOUBLE_EQ(targCom.D.z, dBefore.z) << "D.z must not be modified";

    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, transXBefore)
        << "trajTrans.x must be left for the caller to assign";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, transYBefore)
        << "trajTrans.y must be left for the caller to assign";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, transZBefore)
        << "trajTrans.z must be left for the caller to assign";
}

// -----------------------------------------------------------------------------
// Test 8: a very small sphere still yields a well-formed (finite, on-sphere)
// step, even when the geodesic angle dtheta = dl/R wraps past pi.
// -----------------------------------------------------------------------------
void ccpvs_test_small_sphere_large_angle()
{
    std::cerr << "\n[TEST] ccpvs_test_small_sphere_large_angle\n"
              << "  Function:  create_complex_propagation_vectors_on_sphere()\n"
              << "  Scenario:  R = 1 with dt*D large enough that dl > pi*R, i.e.\n"
              << "             the geodesic angle wraps around the sphere.\n"
              << "  Criterion: the result is finite and the end point is still on\n"
              << "             the sphere (chord <= 2R).\n";

    const double R = 1.0;
    Parameters params = ccpvs_make_params(R, 5.0); // dt*D large compared with R
    const Coord dir { 0.0, 0.0, 1.0 + 0.3 };       // slightly off the pole in x/y? keep simple

    ccpvs_seed_rng(42);

    for (int i = 0; i < 200; ++i) {
        Complex targCom = ccpvs_make_complex(Coord { 0.6, 0.5, 0.7 }, R, 2.0, 2.0);
        const Coord comProj = targCom.comCoord;
        const Coord traj = create_complex_propagation_vectors_on_sphere(params, targCom);

        ASSERT_TRUE(std::isfinite(traj.x) && std::isfinite(traj.y) && std::isfinite(traj.z))
            << "Sample " << i << ": non-finite result on a small sphere";

        const Coord endPt { comProj.x + traj.x, comProj.y + traj.y, comProj.z + traj.z };
        EXPECT_NEAR(ccpvs_norm(endPt), R, 1e-8)
            << "Sample " << i << ": wrapped geodesic step left the sphere";
        EXPECT_LE(ccpvs_norm(traj), 2.0 * R + 1e-9)
            << "Sample " << i << ": chord exceeded the sphere diameter";
    }

    // Silence "unused variable" warnings for the illustrative direction above.
    (void)dir;

    std::cerr << "  200 wrapped-angle samples all stayed on the unit sphere.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is executed inside its own TEST so a
// failure in one does not stop the others from running.
// -----------------------------------------------------------------------------
TEST(CreateComplexPropagationVectorsOnSphere, StaysOnSphere) { ccpvs_test_stays_on_sphere(); }
TEST(CreateComplexPropagationVectorsOnSphere, ProjectionIndependentOfComMagnitude)
{
    ccpvs_test_projection_independent_of_com_magnitude();
}
TEST(CreateComplexPropagationVectorsOnSphere, DeterministicWithSameSeed)
{
    ccpvs_test_deterministic_with_same_seed();
}
TEST(CreateComplexPropagationVectorsOnSphere, StepLengthStatistics)
{
    ccpvs_test_step_length_statistics();
}
TEST(CreateComplexPropagationVectorsOnSphere, AnisotropicDiffusionStatistics)
{
    ccpvs_test_anisotropic_diffusion_statistics();
}
TEST(CreateComplexPropagationVectorsOnSphere, TimestepScaling) { ccpvs_test_timestep_scaling(); }
TEST(CreateComplexPropagationVectorsOnSphere, ComplexIsNotModified)
{
    ccpvs_test_complex_is_not_modified();
}
TEST(CreateComplexPropagationVectorsOnSphere, SmallSphereLargeAngle)
{
    ccpvs_test_small_sphere_large_angle();
}