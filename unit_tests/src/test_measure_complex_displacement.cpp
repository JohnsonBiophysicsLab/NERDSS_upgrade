/*! \file test_measure_complex_displacement.cpp
 *
 * ### Unit test for src/reactions/measure_complex_displacement.cpp
 *
 * Function under test:
 *
 *     void measure_complex_displacement(bool& flag, Complex& reactCom1, Complex& reactCom2,
 *                                       std::vector<Molecule>& moleculeList,
 *                                       const Parameters& params,
 *                                       const std::vector<MolTemplate>& molTemplateList,
 *                                       const std::vector<Complex>& complexList)
 *
 * The routine is called during association.  `Molecule::comCoord` holds the
 * pre-association coordinates while `Molecule::tmpComCoord` holds the proposed
 * (post-association) coordinates.  The routine:
 *
 *   1. Determines a dimensionality for each complex (3 by default, 2 if the
 *      complex is `OnSurface`, 1 if it is `onFiber`).
 *   2. Adds a rotational contribution to the translational diffusion constant.
 *   3. Computes the average Einstein displacement
 *          avgDisp = sqrt(2 * dim * Dtot * timeStep)
 *      and a maximum allowed displacement
 *          maxDisp = params.scaleMaxDisplace * avgDisp
 *   4. Sets `flag = true` (i.e. "cancel association") if either complex COM, or
 *      any member molecule COM, moved further than that maximum.
 *
 * Every test below prints what is being exercised and what the pass criterion
 * is, and all assertions are non-fatal EXPECT_* so the whole file always runs.
 */

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (file-internal, so they cannot clash with other test files).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal Molecule.
 *
 * \param[in] index    index of this Molecule inside moleculeList
 * \param[in] comIndex index of the parent Complex
 * \param[in] com      pre-association centre of mass
 * \param[in] disp     displacement applied to obtain the temporary (proposed)
 *                     centre of mass, i.e. tmpComCoord = com + disp
 */
Molecule mcd_make_molecule(int index, int comIndex, const Coord& com, const Coord& disp)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;

    // "Real" (pre-association) coordinates.
    mol.comCoord = com;
    Molecule::Iface iface;
    iface.coord = com;
    mol.interfaceList.push_back(iface);

    // Temporary (proposed) coordinates used by the routine under test.
    mol.tmpComCoord = Coord(com.x + disp.x, com.y + disp.y, com.z + disp.z);
    mol.tmpICoords.push_back(mol.tmpComCoord);

    return mol;
}

/*! \brief Build a Complex that owns the given member molecules.
 *
 * The Complex centre of mass is set to the (unweighted) average of the member
 * molecules' pre-association coordinates so that the pre/post comparison
 * performed inside the routine is meaningful.
 */
Complex mcd_make_complex(int index, const std::vector<int>& members,
    const std::vector<Molecule>& moleculeList, double diffConst, double rotDiffConst,
    double radius)
{
    Complex com;
    com.index = index;
    com.memberList = members;
    com.numEachMol = std::vector<int>(1, static_cast<int>(members.size()));

    com.D = Coord(diffConst, diffConst, diffConst);
    com.Dr = Coord(rotDiffConst, rotDiffConst, rotDiffConst);
    com.radius = radius;
    com.mass = static_cast<double>(members.size());

    // Average the member coordinates to get the complex COM.
    Coord avg(0.0, 0.0, 0.0);
    for (auto memIdx : members)
        avg += moleculeList[memIdx].comCoord;
    double n = static_cast<double>(members.size());
    avg /= n;
    com.comCoord = avg;
    com.tmpComCoord = avg;

    com.OnSurface = false;
    com.onFiber = false;
    return com;
}

/*! \brief Re-implementation of the maximum-allowed-displacement formula.
 *
 * This mirrors the arithmetic inside measure_complex_displacement() so the
 * tests can probe just below / just above the acceptance threshold.
 *
 * \param[in] com       complex whose threshold is wanted
 * \param[in] dim       dimensionality the routine will use for this complex
 * \param[in] timeStep  Parameters::timeStep
 * \param[in] scale     Parameters::scaleMaxDisplace
 */
double mcd_expected_max_disp(const Complex& com, double dim, double timeStep, double scale)
{
    const double cf = std::cos(std::sqrt(2.0 * (dim - 1.0) * com.Dr.z * timeStep));
    const double Dr = 2.0 * com.radius * com.radius * (1.0 - cf);
    const double Dtot = com.D.x + Dr / (2.0 * dim * timeStep);
    const double avgDisp = std::sqrt(2.0 * dim * Dtot * timeStep);
    return scale * avgDisp;
}

/*! \brief Convenience: default Parameters used by most tests. */
Parameters mcd_make_params(double timeStep, double scaleMaxDisplace)
{
    Parameters params;
    params.timeStep = timeStep;
    params.scaleMaxDisplace = scaleMaxDisplace;
    return params;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: no displacement at all -> association must NOT be cancelled.
// -----------------------------------------------------------------------------
void test_mcd_zero_displacement_keeps_flag_false()
{
    std::cerr << "\n[TEST] test_mcd_zero_displacement_keeps_flag_false\n"
              << "  Source file:   src/reactions/measure_complex_displacement.cpp\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      tmpComCoord == comCoord for every molecule and\n"
              << "                 both complexes (a completely static proposal).\n"
              << "  Pass criteria: flag remains false (association not cancelled).\n";

    Parameters params = mcd_make_params(1.0, 100.0);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList; // unused by the routine

    // Two single-molecule complexes, zero displacement.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));
    moleculeList.push_back(mcd_make_molecule(1, 1, Coord(10.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

    Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.01, 1.0);
    Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.01, 1.0);

    bool flag = false;
    std::cerr << "  Calling measure_complex_displacement()...\n";
    measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);

    EXPECT_FALSE(flag) << "A zero displacement must never cancel association";
    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: complex 1 translated by an enormous distance -> flag set.
// -----------------------------------------------------------------------------
void test_mcd_huge_displacement_complex1_sets_flag()
{
    std::cerr << "\n[TEST] test_mcd_huge_displacement_complex1_sets_flag\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      every member of complex 1 is moved 1e6 nm while\n"
              << "                 complex 2 stays put.\n"
              << "  Pass criteria: flag becomes true (association cancelled).\n";

    Parameters params = mcd_make_params(1.0, 100.0);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    std::vector<Molecule> moleculeList;
    // Complex 1 member: displaced by a physically absurd amount.
    moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(1.0e6, 0.0, 0.0)));
    // Complex 2 member: not displaced.
    moleculeList.push_back(mcd_make_molecule(1, 1, Coord(10.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

    Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.01, 1.0);
    Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.01, 1.0);

    bool flag = false;
    std::cerr << "  Calling measure_complex_displacement()...\n";
    measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);

    EXPECT_TRUE(flag) << "A 1e6 nm displacement of complex 1 must cancel association";
    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
}

// -----------------------------------------------------------------------------
// Test 3: complex 2 translated by an enormous distance -> flag set.
// -----------------------------------------------------------------------------
void test_mcd_huge_displacement_complex2_sets_flag()
{
    std::cerr << "\n[TEST] test_mcd_huge_displacement_complex2_sets_flag\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      complex 1 stays put, complex 2 is moved 1e6 nm.\n"
              << "  Pass criteria: flag becomes true; this exercises the second\n"
              << "                 (reactCom2) branch of the COM check.\n";

    Parameters params = mcd_make_params(1.0, 100.0);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));
    moleculeList.push_back(mcd_make_molecule(1, 1, Coord(10.0, 0.0, 0.0), Coord(0.0, 1.0e6, 0.0)));

    Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.01, 1.0);
    Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.01, 1.0);

    bool flag = false;
    std::cerr << "  Calling measure_complex_displacement()...\n";
    measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);

    EXPECT_TRUE(flag) << "A 1e6 nm displacement of complex 2 must cancel association";
    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
}

// -----------------------------------------------------------------------------
// Test 4: complex COM is conserved but individual molecules move a lot.
//         This isolates the per-member-molecule loops.
// -----------------------------------------------------------------------------
void test_mcd_member_molecule_loop_sets_flag()
{
    std::cerr << "\n[TEST] test_mcd_member_molecule_loop_sets_flag\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      complex 1 has two members displaced by equal and\n"
              << "                 opposite vectors, so the complex COM does NOT move\n"
              << "                 but each molecule moves 500 nm.\n"
              << "  Pass criteria: flag becomes true, proving the per-molecule loop\n"
              << "                 (not just the COM check) rejects the move.\n";

    Parameters params = mcd_make_params(1.0, 10.0);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    std::vector<Molecule> moleculeList;
    // Complex 1: two molecules symmetric about the origin, swapped outward.
    moleculeList.push_back(mcd_make_molecule(0, 0, Coord(-5.0, 0.0, 0.0), Coord(500.0, 0.0, 0.0)));
    moleculeList.push_back(mcd_make_molecule(1, 0, Coord(5.0, 0.0, 0.0), Coord(-500.0, 0.0, 0.0)));
    // Complex 2: a single static molecule.
    moleculeList.push_back(mcd_make_molecule(2, 1, Coord(50.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

    Complex com1 = mcd_make_complex(0, { 0, 1 }, moleculeList, 1.0, 0.01, 5.0);
    Complex com2 = mcd_make_complex(1, { 2 }, moleculeList, 1.0, 0.01, 1.0);

    bool flag = false;
    std::cerr << "  Calling measure_complex_displacement()...\n";
    measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);

    // Informational: the routine recomputes the temporary complex COM internally.
    std::cerr << "  complex 1 pre-assoc COM  = (" << com1.comCoord.x << ", "
              << com1.comCoord.y << ", " << com1.comCoord.z << ")\n"
              << "  complex 1 temporary COM  = (" << com1.tmpComCoord.x << ", "
              << com1.tmpComCoord.y << ", " << com1.tmpComCoord.z << ")\n";

    EXPECT_TRUE(flag) << "Large per-molecule displacements must cancel association even "
                         "when the complex COM is unchanged";
    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
}

// -----------------------------------------------------------------------------
// Test 5: probe the acceptance threshold from below and from above.
// -----------------------------------------------------------------------------
void test_mcd_threshold_boundary_3d()
{
    std::cerr << "\n[TEST] test_mcd_threshold_boundary_3d\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      3D complexes; the maximum allowed displacement is\n"
              << "                 recomputed in the test from the documented formula\n"
              << "                 maxDisp = scaleMaxDisplace*sqrt(2*dim*Dtot*dt).\n"
              << "  Pass criteria: 0.9*maxDisp is accepted (flag false) and\n"
              << "                 1.05*maxDisp is rejected (flag true).\n";

    const double timeStep = 1.0;
    const double scale = 1.0;
    Parameters params = mcd_make_params(timeStep, scale);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    // A probe complex used only to evaluate the threshold formula.
    std::vector<Molecule> probeList;
    probeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));
    Complex probe = mcd_make_complex(0, { 0 }, probeList, 1.0, 0.01, 1.0);
    const double maxDisp = mcd_expected_max_disp(probe, 3.0, timeStep, scale);
    std::cerr << "  Computed 3D maxDisp = " << maxDisp << " nm\n";

    // ---- Sub-case A: displacement safely below the threshold. --------------
    {
        std::vector<Molecule> moleculeList;
        moleculeList.push_back(
            mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.9 * maxDisp, 0.0, 0.0)));
        moleculeList.push_back(
            mcd_make_molecule(1, 1, Coord(20.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

        Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.01, 1.0);
        Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.01, 1.0);

        bool flag = false;
        std::cerr << "  Sub-case A: displacement = " << 0.9 * maxDisp << " nm (below limit)\n";
        measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);
        EXPECT_FALSE(flag) << "Displacement of 0.9*maxDisp must be accepted";
        std::cerr << "    flag = " << std::boolalpha << flag << '\n';
    }

    // ---- Sub-case B: displacement just above the threshold. ----------------
    {
        std::vector<Molecule> moleculeList;
        moleculeList.push_back(
            mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(1.05 * maxDisp, 0.0, 0.0)));
        moleculeList.push_back(
            mcd_make_molecule(1, 1, Coord(20.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

        Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.01, 1.0);
        Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.01, 1.0);

        bool flag = false;
        std::cerr << "  Sub-case B: displacement = " << 1.05 * maxDisp << " nm (above limit)\n";
        measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);
        EXPECT_TRUE(flag) << "Displacement of 1.05*maxDisp must be rejected";
        std::cerr << "    flag = " << std::boolalpha << flag << '\n';
    }
}

// -----------------------------------------------------------------------------
// Test 6: Complex::OnSurface lowers the dimensionality from 3 to 2, which
//         tightens the allowed displacement.
// -----------------------------------------------------------------------------
void test_mcd_on_surface_tightens_threshold()
{
    std::cerr << "\n[TEST] test_mcd_on_surface_tightens_threshold\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      the rotational diffusion constant is zero, so the\n"
              << "                 threshold is scale*sqrt(2*dim*D*dt).  A displacement\n"
              << "                 chosen strictly between the 2D and 3D limits is used.\n"
              << "  Pass criteria: accepted when complex 1 is a 3D complex, rejected\n"
              << "                 when complex 1 is flagged OnSurface (dim = 2).\n";

    const double timeStep = 1.0;
    const double scale = 1.0;
    Parameters params = mcd_make_params(timeStep, scale);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    std::vector<Molecule> probeList;
    probeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));
    Complex probe = mcd_make_complex(0, { 0 }, probeList, 1.0, 0.0 /*no rotation*/, 1.0);

    const double maxDisp3D = mcd_expected_max_disp(probe, 3.0, timeStep, scale);
    const double maxDisp2D = mcd_expected_max_disp(probe, 2.0, timeStep, scale);
    const double testDisp = 0.5 * (maxDisp2D + maxDisp3D); // between the two limits
    std::cerr << "  maxDisp(2D) = " << maxDisp2D << ", maxDisp(3D) = " << maxDisp3D
              << ", test displacement = " << testDisp << " nm\n";

    // Sanity check on the construction of this test.
    EXPECT_LT(maxDisp2D, maxDisp3D) << "With zero rotation the 2D limit must be smaller";

    // ---- Sub-case A: complex 1 is a normal 3D complex -> accepted. ---------
    {
        std::vector<Molecule> moleculeList;
        moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(testDisp, 0.0, 0.0)));
        moleculeList.push_back(mcd_make_molecule(1, 1, Coord(20.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

        Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.0, 1.0);
        Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.0, 1.0);

        bool flag = false;
        std::cerr << "  Sub-case A: complex 1 OnSurface = false (dim = 3)\n";
        measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);
        EXPECT_FALSE(flag) << "Displacement below the 3D limit must be accepted";
        std::cerr << "    flag = " << std::boolalpha << flag << '\n';
    }

    // ---- Sub-case B: same displacement but complex 1 is on the membrane. ---
    {
        std::vector<Molecule> moleculeList;
        moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(testDisp, 0.0, 0.0)));
        moleculeList.push_back(mcd_make_molecule(1, 1, Coord(20.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

        Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.0, 1.0);
        Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.0, 1.0);
        com1.OnSurface = true; // forces dim1 = 2 inside the routine

        bool flag = false;
        std::cerr << "  Sub-case B: complex 1 OnSurface = true (dim = 2)\n";
        measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);
        EXPECT_TRUE(flag) << "Displacement above the 2D limit must be rejected for a "
                             "membrane-bound complex";
        std::cerr << "    flag = " << std::boolalpha << flag << '\n';
    }
}

// -----------------------------------------------------------------------------
// Test 7: Complex::onFiber lowers the dimensionality of complex 1 to 1.
// -----------------------------------------------------------------------------
void test_mcd_on_fiber_tightens_threshold()
{
    std::cerr << "\n[TEST] test_mcd_on_fiber_tightens_threshold\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      complex 1 has onFiber = true, so the routine uses\n"
              << "                 dim1 = 1 and the allowed displacement shrinks to\n"
              << "                 scale*sqrt(2*D*dt).\n"
              << "  Pass criteria: a displacement between the 1D and 3D limits is\n"
              << "                 rejected only when onFiber is set.\n";

    const double timeStep = 1.0;
    const double scale = 1.0;
    Parameters params = mcd_make_params(timeStep, scale);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    std::vector<Molecule> probeList;
    probeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));
    Complex probe = mcd_make_complex(0, { 0 }, probeList, 1.0, 0.0, 1.0);

    const double maxDisp1D = mcd_expected_max_disp(probe, 1.0, timeStep, scale);
    const double maxDisp3D = mcd_expected_max_disp(probe, 3.0, timeStep, scale);
    const double testDisp = 0.5 * (maxDisp1D + maxDisp3D);
    std::cerr << "  maxDisp(1D) = " << maxDisp1D << ", maxDisp(3D) = " << maxDisp3D
              << ", test displacement = " << testDisp << " nm\n";

    EXPECT_LT(maxDisp1D, maxDisp3D) << "The 1D limit must be smaller than the 3D limit";

    // ---- Sub-case A: no fiber -> accepted. ---------------------------------
    {
        std::vector<Molecule> moleculeList;
        moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(testDisp, 0.0, 0.0)));
        moleculeList.push_back(mcd_make_molecule(1, 1, Coord(20.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

        Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.0, 1.0);
        Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.0, 1.0);

        bool flag = false;
        std::cerr << "  Sub-case A: complex 1 onFiber = false (dim = 3)\n";
        measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);
        EXPECT_FALSE(flag) << "Displacement below the 3D limit must be accepted";
        std::cerr << "    flag = " << std::boolalpha << flag << '\n';
    }

    // ---- Sub-case B: complex 1 on a fiber -> rejected. ---------------------
    {
        std::vector<Molecule> moleculeList;
        moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(testDisp, 0.0, 0.0)));
        moleculeList.push_back(mcd_make_molecule(1, 1, Coord(20.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

        Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.0, 1.0);
        Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.0, 1.0);
        com1.onFiber = true; // forces dim1 = 1 inside the routine

        bool flag = false;
        std::cerr << "  Sub-case B: complex 1 onFiber = true (dim = 1)\n";
        measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);
        EXPECT_TRUE(flag) << "Displacement above the 1D limit must be rejected for a "
                             "fiber-bound complex";
        std::cerr << "    flag = " << std::boolalpha << flag << '\n';
    }
}

// -----------------------------------------------------------------------------
// Test 8: documents the current behaviour that reactCom2.onFiber modifies
//         *dim1* (the dimensionality of complex 1) in the implementation.
// -----------------------------------------------------------------------------
void test_mcd_com2_on_fiber_affects_com1_dimensionality()
{
    std::cerr << "\n[TEST] test_mcd_com2_on_fiber_affects_com1_dimensionality\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      only complex 2 is flagged onFiber.  In the current\n"
              << "                 implementation that branch assigns dim1 = 1 (rather\n"
              << "                 than dim2), so COMPLEX 1's limit is the one that\n"
              << "                 shrinks.  This test documents that behaviour.\n"
              << "  Pass criteria: a displacement of complex 1 that lies between the\n"
              << "                 1D and 3D limits is rejected.\n";

    const double timeStep = 1.0;
    const double scale = 1.0;
    Parameters params = mcd_make_params(timeStep, scale);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    std::vector<Molecule> probeList;
    probeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));
    Complex probe = mcd_make_complex(0, { 0 }, probeList, 1.0, 0.0, 1.0);

    const double maxDisp1D = mcd_expected_max_disp(probe, 1.0, timeStep, scale);
    const double maxDisp3D = mcd_expected_max_disp(probe, 3.0, timeStep, scale);
    const double testDisp = 0.5 * (maxDisp1D + maxDisp3D);
    std::cerr << "  displacing complex 1 by " << testDisp << " nm (between "
              << maxDisp1D << " and " << maxDisp3D << ")\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(testDisp, 0.0, 0.0)));
    moleculeList.push_back(mcd_make_molecule(1, 1, Coord(20.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

    Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.0, 1.0);
    Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.0, 1.0);
    com2.onFiber = true; // implementation sets dim1 = 1 in this branch

    bool flag = false;
    std::cerr << "  Calling measure_complex_displacement()...\n";
    measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);

    EXPECT_TRUE(flag) << "Current implementation reduces dim1 when reactCom2.onFiber is "
                         "set, so complex 1's displacement is rejected";
    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
}

// -----------------------------------------------------------------------------
// Test 9: the routine never clears an already-raised flag.
// -----------------------------------------------------------------------------
void test_mcd_preexisting_flag_is_preserved()
{
    std::cerr << "\n[TEST] test_mcd_preexisting_flag_is_preserved\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      flag is already true (an earlier check cancelled\n"
              << "                 association) and the proposed move is perfectly\n"
              << "                 acceptable (zero displacement).\n"
              << "  Pass criteria: flag is still true; the routine only ever raises\n"
              << "                 the flag, it never resets it.\n";

    Parameters params = mcd_make_params(1.0, 100.0);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(mcd_make_molecule(0, 0, Coord(0.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));
    moleculeList.push_back(mcd_make_molecule(1, 1, Coord(10.0, 0.0, 0.0), Coord(0.0, 0.0, 0.0)));

    Complex com1 = mcd_make_complex(0, { 0 }, moleculeList, 1.0, 0.01, 1.0);
    Complex com2 = mcd_make_complex(1, { 1 }, moleculeList, 1.0, 0.01, 1.0);

    bool flag = true;
    std::cerr << "  Calling measure_complex_displacement() with flag pre-set to true...\n";
    measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);

    EXPECT_TRUE(flag) << "An already-raised cancellation flag must not be cleared";
    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
}

// -----------------------------------------------------------------------------
// Test 10: multi-molecule complexes with a small rigid-body translation are
//          accepted (regression style sanity check of the whole code path).
// -----------------------------------------------------------------------------
void test_mcd_small_rigid_translation_accepted()
{
    std::cerr << "\n[TEST] test_mcd_small_rigid_translation_accepted\n"
              << "  Function:      measure_complex_displacement()\n"
              << "  Scenario:      both complexes contain several molecules and are\n"
              << "                 rigidly translated by a small vector (0.1 nm), which\n"
              << "                 is far below the allowed limit.\n"
              << "  Pass criteria: flag stays false and the recomputed temporary complex\n"
              << "                 COM is offset from the original COM by that same\n"
              << "                 small vector.\n";

    Parameters params = mcd_make_params(1.0, 100.0);
    std::vector<MolTemplate> molTemplateList(1);
    std::vector<Complex> complexList;

    const Coord shift(0.1, -0.1, 0.05);

    std::vector<Molecule> moleculeList;
    // Complex 1: three molecules, all shifted by the same small vector.
    moleculeList.push_back(mcd_make_molecule(0, 0, Coord(-2.0, 0.0, 0.0), shift));
    moleculeList.push_back(mcd_make_molecule(1, 0, Coord(0.0, 0.0, 0.0), shift));
    moleculeList.push_back(mcd_make_molecule(2, 0, Coord(2.0, 0.0, 0.0), shift));
    // Complex 2: two molecules, also rigidly shifted.
    moleculeList.push_back(mcd_make_molecule(3, 1, Coord(20.0, 0.0, 0.0), shift));
    moleculeList.push_back(mcd_make_molecule(4, 1, Coord(22.0, 0.0, 0.0), shift));

    Complex com1 = mcd_make_complex(0, { 0, 1, 2 }, moleculeList, 1.0, 0.01, 2.0);
    Complex com2 = mcd_make_complex(1, { 3, 4 }, moleculeList, 1.0, 0.01, 1.0);

    const Coord origCom1 = com1.comCoord;

    bool flag = false;
    std::cerr << "  Calling measure_complex_displacement()...\n";
    measure_complex_displacement(flag, com1, com2, moleculeList, params, molTemplateList, complexList);

    EXPECT_FALSE(flag) << "A 0.1 nm rigid translation must be accepted";

    // The routine calls update_complex_tmp_com_crds(), so the temporary complex
    // COM should now reflect the shifted member coordinates.
    std::cerr << "  complex 1 tmp COM = (" << com1.tmpComCoord.x << ", "
              << com1.tmpComCoord.y << ", " << com1.tmpComCoord.z << ")\n";
    EXPECT_NEAR(com1.tmpComCoord.x, origCom1.x + shift.x, 1e-9)
        << "temporary COM x should track the rigid shift";
    EXPECT_NEAR(com1.tmpComCoord.y, origCom1.y + shift.y, 1e-9)
        << "temporary COM y should track the rigid shift";
    EXPECT_NEAR(com1.tmpComCoord.z, origCom1.z + shift.z, 1e-9)
        << "temporary COM z should track the rigid shift";

    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: one TEST per named test_* function so that failures are
// reported individually while every scenario still executes.
// -----------------------------------------------------------------------------
TEST(MeasureComplexDisplacement, ZeroDisplacementKeepsFlagFalse) { test_mcd_zero_displacement_keeps_flag_false(); }
TEST(MeasureComplexDisplacement, HugeDisplacementComplex1SetsFlag) { test_mcd_huge_displacement_complex1_sets_flag(); }
TEST(MeasureComplexDisplacement, HugeDisplacementComplex2SetsFlag) { test_mcd_huge_displacement_complex2_sets_flag(); }
TEST(MeasureComplexDisplacement, MemberMoleculeLoopSetsFlag) { test_mcd_member_molecule_loop_sets_flag(); }
TEST(MeasureComplexDisplacement, ThresholdBoundary3D) { test_mcd_threshold_boundary_3d(); }
TEST(MeasureComplexDisplacement, OnSurfaceTightensThreshold) { test_mcd_on_surface_tightens_threshold(); }
TEST(MeasureComplexDisplacement, OnFiberTightensThreshold) { test_mcd_on_fiber_tightens_threshold(); }
TEST(MeasureComplexDisplacement, Com2OnFiberAffectsCom1Dim) { test_mcd_com2_on_fiber_affects_com1_dimensionality(); }
TEST(MeasureComplexDisplacement, PreexistingFlagPreserved) { test_mcd_preexisting_flag_is_preserved(); }
TEST(MeasureComplexDisplacement, SmallRigidTranslationAccepted) { test_mcd_small_rigid_translation_accepted(); }