/*! \file test_reflect_traj_tmp_crds_compartment.cpp
 *
 * ### Unit test for src/boundary_conditions/reflect_traj_tmp_crds_compartment.cpp
 *
 * Function under test:
 *
 *     void reflect_traj_tmp_crds_compartment(const Parameters& params,
 *                                            std::vector<Molecule>& moleculeList,
 *                                            Complex& targCom,
 *                                            std::array<double, 3>& traj,
 *                                            const Membrane& membraneObject,
 *                                            double RS3Dinput);
 *
 * ### What the function does (as implemented)
 *
 *  1. Computes an effective compartment radius
 *         sphereR = membraneObject.compartmentR + RS3D,
 *     where RS3D == 0 if the complex is on (or temporarily on) the surface and
 *     RS3D == RS3Dinput otherwise.
 *  2. Performs a cheap pre-check: only if the *whole* complex could fit inside
 *     the compartment sphere
 *         |tmpComCoord + traj| + targCom.radius  <  sphereR
 *     does it look any closer.  Otherwise `traj` is left completely alone.
 *  3. If the pre-check passes, every member molecule COM and every member
 *     interface (using the *temporary* association coordinates) is measured
 *     against the compartment sphere.  The point with the most negative
 *     (magnitude - sphereR), i.e. the point furthest *inside* the compartment,
 *     is kept.
 *  4. `traj` is then overwritten with the reflection that pushes that deepest
 *     point back out through the compartment surface:
 *         lamda = -2 * (|p| - sphereR) / |p|
 *         traj  = lamda * p
 *
 *     Note: only `traj` is modified; neither the Complex nor the Molecules are
 *     touched.
 *
 * The tests below verify each of those behaviours with exact, hand-computable
 * geometry so the pass/fail criteria are unambiguous.  Verbose output is sent
 * to stderr so the reader can follow which scenario is being exercised.
 */

#include "boundary_conditions/reflect_functions.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Small helpers used by all of the tests in this file.  Everything is prefixed
// with `rttcc_` (Reflect Traj Tmp Crds Compartment) so it cannot collide with
// helpers defined in other translation units of the suite.
// -----------------------------------------------------------------------------

/*! \brief Build a Molecule whose temporary (association) coordinates are set.
 *
 * The function under test reads only `tmpComCoord`, `tmpICoords` and
 * `interfaceList.size()`, so the interface list is sized to match the temporary
 * interface coordinates.
 *
 * \param[in] tmpCom     temporary center-of-mass coordinate
 * \param[in] tmpIfaces  temporary interface coordinates (may be empty)
 * \return a Molecule ready to be dropped into a moleculeList
 */
Molecule rttcc_make_molecule(const Coord& tmpCom, const std::vector<Coord>& tmpIfaces)
{
    Molecule mol;
    mol.comCoord = tmpCom; // real coords, unused by the function
    mol.tmpComCoord = tmpCom; // temporary coords, this is what is read
    mol.interfaceList.clear();
    mol.tmpICoords.clear();
    for (const auto& oneIface : tmpIfaces) {
        Molecule::Iface iface;
        iface.coord = oneIface; // real coord, unused
        mol.interfaceList.push_back(iface);
        mol.tmpICoords.push_back(oneIface); // temporary coord, this is read
    }
    mol.myComIndex = 0;
    return mol;
}

/*! \brief Build a Complex with the given temporary COM, bounding radius and members. */
Complex rttcc_make_complex(const Coord& tmpCom, double radius, const std::vector<int>& members)
{
    Complex targCom;
    targCom.comCoord = tmpCom;
    targCom.tmpComCoord = tmpCom; // the function only uses tmpComCoord
    targCom.radius = radius;

    // Diffusion constants: not used by this routine (no resampling here), but
    // set to sane non-zero values anyway.
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Coord { 0.01, 0.01, 0.01 };

    targCom.memberList = members;
    targCom.OnSurface = false;
    targCom.tmpOnSurface = false;
    return targCom;
}

/*! \brief Build a Membrane that owns a spherical compartment of the given radius. */
Membrane rttcc_make_membrane(double compartmentR)
{
    Membrane membraneObject;
    membraneObject.hasCompartment = true;
    membraneObject.compartmentR = compartmentR;
    // A water box is not consulted by this routine but is filled in for realism.
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { 1000.0, 1000.0, 1000.0 });
    return membraneObject;
}

/*! \brief Reference implementation of the reflection the routine should apply.
 *
 * \param[in] deepestPoint the point furthest inside the compartment
 * \param[in] sphereR      effective compartment radius
 * \return the expected trajectory vector
 */
std::array<double, 3> rttcc_expected_traj(const Coord& deepestPoint, double sphereR)
{
    Coord p { deepestPoint };
    const double mag = p.get_magnitude();
    const double lamda = -2.0 * (mag - sphereR) / mag;
    return { lamda * p.x, lamda * p.y, lamda * p.z };
}

/*! \brief Convenience printer for a trajectory triple. */
void rttcc_print_traj(const char* label, const std::array<double, 3>& traj)
{
    std::cerr << "    " << label << " = (" << traj[0] << ", " << traj[1] << ", " << traj[2] << ")\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: complex nowhere near the compartment -> trajectory untouched.
// -----------------------------------------------------------------------------
void test_rttcc_far_outside_leaves_traj_unchanged()
{
    std::cerr << "\n[TEST] test_rttcc_far_outside_leaves_traj_unchanged\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      complex sits at |r| = 200 while the compartment\n"
              << "                 radius is only 100, so the cheap pre-check\n"
              << "                 (|com+traj| + radius < sphereR) fails.\n"
              << "  Pass criteria: traj is bit-for-bit unchanged.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);

    // One molecule far outside the compartment.
    std::vector<Molecule> moleculeList { rttcc_make_molecule(Coord { 200.0, 0.0, 0.0 },
        { Coord { 201.0, 0.0, 0.0 } }) };
    Complex targCom = rttcc_make_complex(Coord { 200.0, 0.0, 0.0 }, 1.0, { 0 });

    std::array<double, 3> traj { 1.0, 2.0, 3.0 };
    const std::array<double, 3> orig = traj;

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);

    EXPECT_DOUBLE_EQ(traj[0], orig[0]) << "traj[0] must not change for a complex outside the compartment";
    EXPECT_DOUBLE_EQ(traj[1], orig[1]) << "traj[1] must not change for a complex outside the compartment";
    EXPECT_DOUBLE_EQ(traj[2], orig[2]) << "traj[2] must not change for a complex outside the compartment";
}

// -----------------------------------------------------------------------------
// Test 2: complex fully inside the compartment -> reflected back out.
// -----------------------------------------------------------------------------
void test_rttcc_inside_reflects_outward()
{
    std::cerr << "\n[TEST] test_rttcc_inside_reflects_outward\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      single molecule at (10,0,0), compartment R = 100,\n"
              << "                 initial traj = (0,0,0). The deepest point is\n"
              << "                 (10,0,0) so lamda = -2*(10-100)/10 = 18.\n"
              << "  Pass criteria: traj == (180,0,0) and the final position is at\n"
              << "                 |2*sphereR - 10| = 190 from the origin.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);

    // Molecule COM and its single interface both at (10,0,0) so the deepest
    // point is unambiguous.
    std::vector<Molecule> moleculeList { rttcc_make_molecule(Coord { 10.0, 0.0, 0.0 },
        { Coord { 10.0, 0.0, 0.0 } }) };
    Complex targCom = rttcc_make_complex(Coord { 10.0, 0.0, 0.0 }, 1.0, { 0 });

    std::array<double, 3> traj { 0.0, 0.0, 0.0 };

    const std::array<double, 3> expected = rttcc_expected_traj(Coord { 10.0, 0.0, 0.0 }, 100.0);

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);
    rttcc_print_traj("traj expected  ", expected);

    EXPECT_NEAR(traj[0], expected[0], 1e-9) << "traj[0] should equal lamda * deepestPoint.x";
    EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] should equal lamda * deepestPoint.y";
    EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] should equal lamda * deepestPoint.z";

    // Sanity check that the reflection actually pushed us outside the compartment.
    const double finalMag = std::sqrt(std::pow(targCom.tmpComCoord.x + traj[0], 2)
        + std::pow(targCom.tmpComCoord.y + traj[1], 2)
        + std::pow(targCom.tmpComCoord.z + traj[2], 2));
    std::cerr << "    final |com + traj| = " << finalMag << " (compartment R = 100)\n";
    EXPECT_GT(finalMag, 100.0) << "After reflection the complex must be outside the compartment";
    EXPECT_NEAR(finalMag, 190.0, 1e-9) << "Reflection should mirror the point across the compartment surface";
}

// -----------------------------------------------------------------------------
// Test 3: an interface coordinate deeper than the COM must be the one selected.
// -----------------------------------------------------------------------------
void test_rttcc_deepest_interface_is_selected()
{
    std::cerr << "\n[TEST] test_rttcc_deepest_interface_is_selected\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      molecule COM at (10,0,0) but one interface sits at\n"
              << "                 (2,0,0), i.e. deeper inside the compartment.\n"
              << "  Pass criteria: reflection is computed from (2,0,0), giving\n"
              << "                 lamda = -2*(2-100)/2 = 98 -> traj = (196,0,0).\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);

    // Two interfaces: one at the COM, one much deeper inside.
    std::vector<Molecule> moleculeList { rttcc_make_molecule(Coord { 10.0, 0.0, 0.0 },
        { Coord { 10.0, 0.0, 0.0 }, Coord { 2.0, 0.0, 0.0 } }) };
    Complex targCom = rttcc_make_complex(Coord { 10.0, 0.0, 0.0 }, 8.0, { 0 });

    std::array<double, 3> traj { 0.0, 0.0, 0.0 };

    const std::array<double, 3> expected = rttcc_expected_traj(Coord { 2.0, 0.0, 0.0 }, 100.0);

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);
    rttcc_print_traj("traj expected  ", expected);

    EXPECT_NEAR(traj[0], expected[0], 1e-9) << "Deepest interface (2,0,0) should drive the reflection";
    EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] should stay zero for an on-axis geometry";
    EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] should stay zero for an on-axis geometry";
}

// -----------------------------------------------------------------------------
// Test 4: with several member molecules the deepest COM must win.
// -----------------------------------------------------------------------------
void test_rttcc_deepest_member_molecule_is_selected()
{
    std::cerr << "\n[TEST] test_rttcc_deepest_member_molecule_is_selected\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      two-member complex, molecule 0 at (10,0,0) and\n"
              << "                 molecule 1 at (5,0,0) (deeper inside).\n"
              << "  Pass criteria: reflection is computed from (5,0,0), giving\n"
              << "                 lamda = -2*(5-100)/5 = 38 -> traj = (190,0,0).\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);

    // Both molecules keep their single interface on top of their own COM so the
    // COM comparison is what decides the deepest point.
    std::vector<Molecule> moleculeList {
        rttcc_make_molecule(Coord { 10.0, 0.0, 0.0 }, { Coord { 10.0, 0.0, 0.0 } }),
        rttcc_make_molecule(Coord { 5.0, 0.0, 0.0 }, { Coord { 5.0, 0.0, 0.0 } })
    };
    // Complex COM taken as molecule 0's position; radius small enough to pass
    // the cheap pre-check (10 + 6 < 100).
    Complex targCom = rttcc_make_complex(Coord { 10.0, 0.0, 0.0 }, 6.0, { 0, 1 });
    moleculeList[1].myComIndex = 0;

    std::array<double, 3> traj { 0.0, 0.0, 0.0 };

    const std::array<double, 3> expected = rttcc_expected_traj(Coord { 5.0, 0.0, 0.0 }, 100.0);

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);
    rttcc_print_traj("traj expected  ", expected);

    EXPECT_NEAR(traj[0], expected[0], 1e-9) << "Deepest member COM (5,0,0) should drive the reflection";
    EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] should remain zero";
    EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] should remain zero";
}

// -----------------------------------------------------------------------------
// Test 5: a non-zero incoming trajectory is folded into the test point.
// -----------------------------------------------------------------------------
void test_rttcc_incoming_traj_is_included()
{
    std::cerr << "\n[TEST] test_rttcc_incoming_traj_is_included\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      molecule at (10,0,0) with an incoming trajectory\n"
              << "                 of (5,0,0); the tested point is therefore (15,0,0).\n"
              << "  Pass criteria: traj is overwritten with lamda*(15,0,0) where\n"
              << "                 lamda = -2*(15-100)/15.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);

    std::vector<Molecule> moleculeList { rttcc_make_molecule(Coord { 10.0, 0.0, 0.0 },
        { Coord { 10.0, 0.0, 0.0 } }) };
    Complex targCom = rttcc_make_complex(Coord { 10.0, 0.0, 0.0 }, 1.0, { 0 });

    std::array<double, 3> traj { 5.0, 0.0, 0.0 };

    // The candidate point is tmpComCoord + traj = (15,0,0).
    const std::array<double, 3> expected = rttcc_expected_traj(Coord { 15.0, 0.0, 0.0 }, 100.0);

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);
    rttcc_print_traj("traj expected  ", expected);

    EXPECT_NEAR(traj[0], expected[0], 1e-9) << "Incoming traj must shift the tested point before reflection";
    EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] should remain zero";
    EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] should remain zero";
}

// -----------------------------------------------------------------------------
// Test 6: RS3Dinput enlarges the effective radius, unless the complex is on the
//         surface, in which case RS3D is forced to zero.
// -----------------------------------------------------------------------------
void test_rttcc_rs3d_and_on_surface_handling()
{
    std::cerr << "\n[TEST] test_rttcc_rs3d_and_on_surface_handling\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario A:    OnSurface == false, RS3Dinput = 5 -> sphereR = 105.\n"
              << "  Scenario B:    OnSurface == true,  RS3Dinput = 5 -> sphereR = 100.\n"
              << "  Scenario C:    tmpOnSurface == true, RS3Dinput = 5 -> sphereR = 100.\n"
              << "  Pass criteria: each reflection matches the analytic value for the\n"
              << "                 corresponding effective radius.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);
    const double rs3dInput = 5.0;
    const Coord deepest { 10.0, 0.0, 0.0 };

    // --- Scenario A: solution complex, RS3D added to the compartment radius ---
    {
        std::vector<Molecule> moleculeList { rttcc_make_molecule(deepest, { deepest }) };
        Complex targCom = rttcc_make_complex(deepest, 1.0, { 0 });
        targCom.OnSurface = false;
        targCom.tmpOnSurface = false;

        std::array<double, 3> traj { 0.0, 0.0, 0.0 };
        const std::array<double, 3> expected = rttcc_expected_traj(deepest, 100.0 + rs3dInput);

        std::cerr << "  Scenario A: calling with OnSurface=false, RS3Dinput=5...\n";
        reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, rs3dInput);
        rttcc_print_traj("traj after call", traj);
        rttcc_print_traj("traj expected  ", expected);

        EXPECT_NEAR(traj[0], expected[0], 1e-9) << "RS3Dinput should enlarge sphereR for a 3D complex";
        EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] should remain zero";
        EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] should remain zero";
    }

    // --- Scenario B: complex already on the surface -> RS3D forced to zero ---
    {
        std::vector<Molecule> moleculeList { rttcc_make_molecule(deepest, { deepest }) };
        Complex targCom = rttcc_make_complex(deepest, 1.0, { 0 });
        targCom.OnSurface = true;

        std::array<double, 3> traj { 0.0, 0.0, 0.0 };
        const std::array<double, 3> expected = rttcc_expected_traj(deepest, 100.0);

        std::cerr << "  Scenario B: calling with OnSurface=true, RS3Dinput=5...\n";
        reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, rs3dInput);
        rttcc_print_traj("traj after call", traj);
        rttcc_print_traj("traj expected  ", expected);

        EXPECT_NEAR(traj[0], expected[0], 1e-9) << "OnSurface complexes must ignore RS3Dinput";
        EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] should remain zero";
        EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] should remain zero";
    }

    // --- Scenario C: complex temporarily on the surface -> RS3D forced to zero ---
    {
        std::vector<Molecule> moleculeList { rttcc_make_molecule(deepest, { deepest }) };
        Complex targCom = rttcc_make_complex(deepest, 1.0, { 0 });
        targCom.OnSurface = false;
        targCom.tmpOnSurface = true;

        std::array<double, 3> traj { 0.0, 0.0, 0.0 };
        const std::array<double, 3> expected = rttcc_expected_traj(deepest, 100.0);

        std::cerr << "  Scenario C: calling with tmpOnSurface=true, RS3Dinput=5...\n";
        reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, rs3dInput);
        rttcc_print_traj("traj after call", traj);
        rttcc_print_traj("traj expected  ", expected);

        EXPECT_NEAR(traj[0], expected[0], 1e-9) << "tmpOnSurface complexes must ignore RS3Dinput";
        EXPECT_NEAR(traj[1], expected[1], 1e-9) << "traj[1] should remain zero";
        EXPECT_NEAR(traj[2], expected[2], 1e-9) << "traj[2] should remain zero";
    }
}

// -----------------------------------------------------------------------------
// Test 7: a complex whose bounding radius spans the compartment fails the cheap
//         pre-check, so no reflection is applied even though its COM is inside.
// -----------------------------------------------------------------------------
void test_rttcc_large_radius_skips_check()
{
    std::cerr << "\n[TEST] test_rttcc_large_radius_skips_check\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      COM is at |r| = 10 (inside a compartment of R = 100)\n"
              << "                 but targCom.radius = 95, so |r| + radius >= R and the\n"
              << "                 cheap pre-check `canBeInsphere` is false.\n"
              << "  Pass criteria: traj is left unchanged.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);

    std::vector<Molecule> moleculeList { rttcc_make_molecule(Coord { 10.0, 0.0, 0.0 },
        { Coord { 10.0, 0.0, 0.0 } }) };
    Complex targCom = rttcc_make_complex(Coord { 10.0, 0.0, 0.0 }, 95.0, { 0 });

    std::array<double, 3> traj { 0.25, -0.5, 0.75 };
    const std::array<double, 3> orig = traj;

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);

    EXPECT_DOUBLE_EQ(traj[0], orig[0]) << "A complex spanning the compartment must not be reflected";
    EXPECT_DOUBLE_EQ(traj[1], orig[1]) << "A complex spanning the compartment must not be reflected";
    EXPECT_DOUBLE_EQ(traj[2], orig[2]) << "A complex spanning the compartment must not be reflected";
}

// -----------------------------------------------------------------------------
// Test 8: an empty member list means no point is ever inspected -> no change.
// -----------------------------------------------------------------------------
void test_rttcc_empty_member_list_no_change()
{
    std::cerr << "\n[TEST] test_rttcc_empty_member_list_no_change\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      the complex passes the cheap pre-check but owns no\n"
              << "                 member molecules, so `inside` never becomes true.\n"
              << "  Pass criteria: traj is left unchanged (and the call does not crash).\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);

    std::vector<Molecule> moleculeList {}; // deliberately empty
    Complex targCom = rttcc_make_complex(Coord { 10.0, 0.0, 0.0 }, 1.0, {});

    std::array<double, 3> traj { 1.0, 1.0, 1.0 };
    const std::array<double, 3> orig = traj;

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);

    EXPECT_DOUBLE_EQ(traj[0], orig[0]) << "No members -> nothing measured -> traj unchanged";
    EXPECT_DOUBLE_EQ(traj[1], orig[1]) << "No members -> nothing measured -> traj unchanged";
    EXPECT_DOUBLE_EQ(traj[2], orig[2]) << "No members -> nothing measured -> traj unchanged";
}

// -----------------------------------------------------------------------------
// Test 9: the routine must be side-effect free apart from `traj`.
// -----------------------------------------------------------------------------
void test_rttcc_does_not_modify_coordinates()
{
    std::cerr << "\n[TEST] test_rttcc_does_not_modify_coordinates\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      a reflection *is* performed (complex inside the\n"
              << "                 compartment) and we then inspect the Complex and\n"
              << "                 Molecule coordinates.\n"
              << "  Pass criteria: Complex::tmpComCoord, Complex::trajTrans, and all\n"
              << "                 molecule temporary coordinates are untouched.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(100.0);

    const Coord molCom { 12.0, -3.0, 4.0 };
    const Coord molIface { 13.0, -3.0, 4.0 };

    std::vector<Molecule> moleculeList { rttcc_make_molecule(molCom, { molIface }) };
    Complex targCom = rttcc_make_complex(molCom, 2.0, { 0 });

    // Remember what should not change.
    const Coord origComTmp = targCom.tmpComCoord;
    const Vector origTrajTrans = targCom.trajTrans;
    const Coord origMolTmpCom = moleculeList[0].tmpComCoord;
    const Coord origMolTmpIface = moleculeList[0].tmpICoords[0];

    std::array<double, 3> traj { 0.0, 0.0, 0.0 };

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);

    // First confirm the routine actually did something, otherwise the
    // "no side effects" assertions below would be vacuous.
    const double trajMag = std::sqrt(traj[0] * traj[0] + traj[1] * traj[1] + traj[2] * traj[2]);
    std::cerr << "    |traj| = " << trajMag << " (should be > 0, i.e. a reflection happened)\n";
    EXPECT_GT(trajMag, 0.0) << "A reflection should have been applied for this geometry";

    // Complex coordinates must be untouched.
    EXPECT_DOUBLE_EQ(targCom.tmpComCoord.x, origComTmp.x) << "Complex tmpComCoord.x must not be modified";
    EXPECT_DOUBLE_EQ(targCom.tmpComCoord.y, origComTmp.y) << "Complex tmpComCoord.y must not be modified";
    EXPECT_DOUBLE_EQ(targCom.tmpComCoord.z, origComTmp.z) << "Complex tmpComCoord.z must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.x, origTrajTrans.x) << "Complex trajTrans.x must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.y, origTrajTrans.y) << "Complex trajTrans.y must not be modified";
    EXPECT_DOUBLE_EQ(targCom.trajTrans.z, origTrajTrans.z) << "Complex trajTrans.z must not be modified";

    // Molecule temporary coordinates must be untouched.
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpComCoord.x, origMolTmpCom.x) << "Molecule tmpComCoord.x must not change";
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpComCoord.y, origMolTmpCom.y) << "Molecule tmpComCoord.y must not change";
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpComCoord.z, origMolTmpCom.z) << "Molecule tmpComCoord.z must not change";
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpICoords[0].x, origMolTmpIface.x) << "Molecule tmpICoords[0].x must not change";
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpICoords[0].y, origMolTmpIface.y) << "Molecule tmpICoords[0].y must not change";
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpICoords[0].z, origMolTmpIface.z) << "Molecule tmpICoords[0].z must not change";
}

// -----------------------------------------------------------------------------
// Test 10: an off-axis geometry, to confirm all three components are reflected.
// -----------------------------------------------------------------------------
void test_rttcc_off_axis_reflection_components()
{
    std::cerr << "\n[TEST] test_rttcc_off_axis_reflection_components\n"
              << "  Source file:   reflect_traj_tmp_crds_compartment.cpp\n"
              << "  Function:      reflect_traj_tmp_crds_compartment\n"
              << "  Scenario:      deepest point at (6,8,0) (|r| = 10) inside a\n"
              << "                 compartment of R = 50.\n"
              << "  Pass criteria: traj == lamda*(6,8,0) with lamda = -2*(10-50)/10 = 8\n"
              << "                 and the reflected point lies radially outward.\n";

    Parameters params;
    params.timeStep = 1.0;

    Membrane membraneObject = rttcc_make_membrane(50.0);

    const Coord deepest { 6.0, 8.0, 0.0 }; // magnitude exactly 10
    std::vector<Molecule> moleculeList { rttcc_make_molecule(deepest, { deepest }) };
    Complex targCom = rttcc_make_complex(deepest, 1.0, { 0 });

    std::array<double, 3> traj { 0.0, 0.0, 0.0 };
    const std::array<double, 3> expected = rttcc_expected_traj(deepest, 50.0);

    std::cerr << "  Calling reflect_traj_tmp_crds_compartment...\n";
    reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, 0.0);
    rttcc_print_traj("traj after call", traj);
    rttcc_print_traj("traj expected  ", expected);

    EXPECT_NEAR(traj[0], expected[0], 1e-9) << "x component of the reflection is wrong";
    EXPECT_NEAR(traj[1], expected[1], 1e-9) << "y component of the reflection is wrong";
    EXPECT_NEAR(traj[2], expected[2], 1e-9) << "z component of the reflection is wrong";

    // The reflected point should be collinear with, and outside of, the original.
    const double fx = deepest.x + traj[0];
    const double fy = deepest.y + traj[1];
    const double fz = deepest.z + traj[2];
    const double finalMag = std::sqrt(fx * fx + fy * fy + fz * fz);
    std::cerr << "    final |point| = " << finalMag << " (expected 2*50 - 10 = 90)\n";
    EXPECT_NEAR(finalMag, 90.0, 1e-9) << "Reflected point should be mirrored across the compartment surface";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each scenario is run in its own TEST so a failure in one
// does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(ReflectTrajTmpCrdsCompartment, FarOutsideLeavesTrajUnchanged) { test_rttcc_far_outside_leaves_traj_unchanged(); }
TEST(ReflectTrajTmpCrdsCompartment, InsideReflectsOutward) { test_rttcc_inside_reflects_outward(); }
TEST(ReflectTrajTmpCrdsCompartment, DeepestInterfaceIsSelected) { test_rttcc_deepest_interface_is_selected(); }
TEST(ReflectTrajTmpCrdsCompartment, DeepestMemberMoleculeIsSelected) { test_rttcc_deepest_member_molecule_is_selected(); }
TEST(ReflectTrajTmpCrdsCompartment, IncomingTrajIsIncluded) { test_rttcc_incoming_traj_is_included(); }
TEST(ReflectTrajTmpCrdsCompartment, RS3DAndOnSurfaceHandling) { test_rttcc_rs3d_and_on_surface_handling(); }
TEST(ReflectTrajTmpCrdsCompartment, LargeRadiusSkipsCheck) { test_rttcc_large_radius_skips_check(); }
TEST(ReflectTrajTmpCrdsCompartment, EmptyMemberListNoChange) { test_rttcc_empty_member_list_no_change(); }
TEST(ReflectTrajTmpCrdsCompartment, DoesNotModifyCoordinates) { test_rttcc_does_not_modify_coordinates(); }
TEST(ReflectTrajTmpCrdsCompartment, OffAxisReflectionComponents) { test_rttcc_off_axis_reflection_components(); }