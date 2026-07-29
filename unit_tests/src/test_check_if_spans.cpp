/*! \file test_check_if_spans.cpp
 *
 * ### Unit test for ../src/boundary_conditions/check_if_spans.cpp
 *
 * The file under test contains a single function:
 *
 *     void check_if_spans(bool& cancelAssoc, const Parameters& params,
 *                         Complex& reactCom1, Complex& reactCom2,
 *                         std::vector<Molecule>& moleculeList,
 *                         const Membrane& membraneObject)
 *
 * It is a thin dispatcher used during association:
 *   - When the boundary is a *box* (membraneObject.isSphere == false) it
 *     forwards the call to check_if_spans_box(), which decides whether the
 *     newly formed (still temporary) complex is so large that it spans the
 *     entire simulation volume.  If it does, the association is cancelled
 *     (cancelAssoc is set to true).  If the complex only pokes out of one
 *     wall, it is pushed back inside and the association proceeds.
 *   - When the boundary is a *sphere* (membraneObject.isSphere == true) the
 *     function does nothing at all: no flags and no coordinates change.
 *
 * All geometry checks in the association machinery are performed on the
 * *temporary* association coordinates (Molecule::tmpComCoord and
 * Molecule::tmpICoords), so the fixtures below always populate those via
 * Molecule::set_tmp_association_coords().
 *
 * Verbose output is written to stderr so a reader can see which source file
 * and function are exercised, what geometry is used, and what the pass
 * criteria are for each assertion.
 */

#include "boundary_conditions/reflect_functions.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local fixture helpers (prefixed with cis_ == check_if_spans) so they cannot
// collide with helpers in other translation units of the test suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Create a minimal Molecule for the span check.
 *
 * The molecule owns one interface that sits exactly on its center of mass,
 * which keeps the geometry (and therefore the expected bounding box) trivial.
 * set_tmp_association_coords() is called so that tmpComCoord / tmpICoords are
 * populated, since check_if_spans_box() works on the temporary coordinates.
 *
 * \param[in] com       Center-of-mass (and interface) coordinate.
 * \param[in] index     Index this molecule will occupy in moleculeList.
 * \param[in] comIndex  Index of the parent Complex.
 */
Molecule cis_make_molecule(const Coord& com, int index, int comIndex)
{
    Molecule mol;
    mol.index = index;
    mol.myComIndex = comIndex;
    mol.molTypeIndex = 0;
    mol.mass = 1.0;
    mol.comCoord = com;

    // One interface, coincident with the COM.
    Molecule::Iface iface;
    iface.coord = com;
    iface.index = 0;
    iface.relIndex = 0;
    iface.molTypeIndex = 0;
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    // Copy comCoord/interfaceList into tmpComCoord/tmpICoords.
    mol.set_tmp_association_coords();
    return mol;
}

/*! \brief Create a Complex that owns the given member molecules.
 *
 * \param[in] index    Index of this Complex in complexList.
 * \param[in] com      Center-of-mass of the complex (also copied to tmpComCoord).
 * \param[in] members  Indices into moleculeList belonging to this complex.
 */
Complex cis_make_complex(int index, const Coord& com, const std::vector<int>& members)
{
    Complex targCom;
    targCom.index = index;
    targCom.comCoord = com;
    targCom.tmpComCoord = com;
    targCom.memberList = members;
    targCom.radius = 1.0;
    targCom.mass = static_cast<double>(members.size());

    // Non-zero diffusion constants in case the callee needs to resample a move.
    targCom.D = Coord { 1.0, 1.0, 1.0 };
    targCom.Dr = Coord { 0.01, 0.01, 0.01 };

    // No pending propagation.
    targCom.trajTrans.x = 0.0;
    targCom.trajTrans.y = 0.0;
    targCom.trajTrans.z = 0.0;
    targCom.trajRot.x = 0.0;
    targCom.trajRot.y = 0.0;
    targCom.trajRot.z = 0.0;

    targCom.OnSurface = false;
    return targCom;
}

/*! \brief Build a cubic, reflecting box Membrane of the requested side length. */
Membrane cis_make_box_membrane(double side)
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.waterBox = Membrane::WaterBox(std::vector<double> { side, side, side });
    membraneObject.xBCtype = "reflect";
    membraneObject.yBCtype = "reflect";
    membraneObject.zBCtype = "reflect";
    return membraneObject;
}

/*! \brief Return the largest tmpComCoord.x over all molecules in the list. */
double cis_max_tmp_x(const std::vector<Molecule>& moleculeList)
{
    double maxX = -std::numeric_limits<double>::infinity();
    for (const auto& mol : moleculeList)
        maxX = std::max(maxX, mol.tmpComCoord.x);
    return maxX;
}

/*! \brief Return the smallest tmpComCoord.x over all molecules in the list. */
double cis_min_tmp_x(const std::vector<Molecule>& moleculeList)
{
    double minX = std::numeric_limits<double>::infinity();
    for (const auto& mol : moleculeList)
        minX = std::min(minX, mol.tmpComCoord.x);
    return minX;
}

/*! \brief Simple default Parameters object (only timeStep matters here). */
Parameters cis_make_params()
{
    Parameters params;
    params.timeStep = 1.0;
    params.numMolTypes = 1;
    return params;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Spherical boundary -> the dispatcher must do absolutely nothing.
//
// Pass criteria: cancelAssoc stays false AND every temporary coordinate is
// bit-for-bit unchanged, even though the two "complexes" are placed far
// outside any reasonable volume.
// -----------------------------------------------------------------------------
void test_cis_sphere_is_a_noop()
{
    std::cerr << "\n[TEST] test_cis_sphere_is_a_noop\n"
              << "  Source file:   check_if_spans.cpp\n"
              << "  Function:      check_if_spans (sphere branch)\n"
              << "  Scenario:      membraneObject.isSphere == true, molecules placed\n"
              << "                 at x = -500 and x = +500 (grossly out of bounds).\n"
              << "  Pass criteria: cancelAssoc remains false and no temporary\n"
              << "                 coordinate is modified (the function returns early).\n";

    Parameters params = cis_make_params();

    // Spherical boundary: check_if_spans() should not call any child function.
    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 50.0;

    // Deliberately extreme positions -- they must be ignored.
    std::vector<Molecule> moleculeList {
        cis_make_molecule(Coord { -500.0, 0.0, 0.0 }, 0, 0),
        cis_make_molecule(Coord { 500.0, 0.0, 0.0 }, 1, 1)
    };
    Complex reactCom1 = cis_make_complex(0, Coord { -500.0, 0.0, 0.0 }, std::vector<int> { 0 });
    Complex reactCom2 = cis_make_complex(1, Coord { 500.0, 0.0, 0.0 }, std::vector<int> { 1 });

    bool cancelAssoc { false };

    std::cerr << "  Calling check_if_spans...\n";
    check_if_spans(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    // The sphere branch is a no-op, so nothing may change.
    EXPECT_FALSE(cancelAssoc)
        << "cancelAssoc must stay false: the sphere branch performs no checks";
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpComCoord.x, -500.0)
        << "molecule 0 tmpComCoord.x must be untouched on the sphere branch";
    EXPECT_DOUBLE_EQ(moleculeList[1].tmpComCoord.x, 500.0)
        << "molecule 1 tmpComCoord.x must be untouched on the sphere branch";
    EXPECT_DOUBLE_EQ(reactCom1.tmpComCoord.x, -500.0)
        << "complex 1 tmpComCoord.x must be untouched on the sphere branch";
    EXPECT_DOUBLE_EQ(reactCom2.tmpComCoord.x, 500.0)
        << "complex 2 tmpComCoord.x must be untouched on the sphere branch";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc
              << ", molecule x positions = (" << moleculeList[0].tmpComCoord.x << ", "
              << moleculeList[1].tmpComCoord.x << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: Spherical boundary with an already-raised cancel flag.
//
// Pass criteria: the flag is passed through untouched (the function neither
// sets nor clears it when the boundary is a sphere).
// -----------------------------------------------------------------------------
void test_cis_sphere_preserves_cancel_flag()
{
    std::cerr << "\n[TEST] test_cis_sphere_preserves_cancel_flag\n"
              << "  Source file:   check_if_spans.cpp\n"
              << "  Function:      check_if_spans (sphere branch)\n"
              << "  Scenario:      cancelAssoc is already true when the function is\n"
              << "                 called with a spherical boundary.\n"
              << "  Pass criteria: cancelAssoc is still true afterwards (no reset).\n";

    Parameters params = cis_make_params();

    Membrane membraneObject;
    membraneObject.isSphere = true;
    membraneObject.sphereR = 100.0;

    std::vector<Molecule> moleculeList {
        cis_make_molecule(Coord { 0.0, 0.0, 0.0 }, 0, 0),
        cis_make_molecule(Coord { 1.0, 0.0, 0.0 }, 1, 1)
    };
    Complex reactCom1 = cis_make_complex(0, Coord { 0.0, 0.0, 0.0 }, std::vector<int> { 0 });
    Complex reactCom2 = cis_make_complex(1, Coord { 1.0, 0.0, 0.0 }, std::vector<int> { 1 });

    // A previous check (e.g. overlap detection) already cancelled association.
    bool cancelAssoc { true };

    std::cerr << "  Calling check_if_spans with cancelAssoc pre-set to true...\n";
    check_if_spans(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_TRUE(cancelAssoc)
        << "an already-raised cancelAssoc flag must survive the sphere branch";

    std::cerr << "  cancelAssoc after call = " << std::boolalpha << cancelAssoc << '\n';
}

// -----------------------------------------------------------------------------
// Test 3: Box boundary, both complexes comfortably inside the volume.
//
// Pass criteria: association is NOT cancelled, and because nothing sticks out
// of a wall the temporary molecule coordinates are left exactly as supplied.
// -----------------------------------------------------------------------------
void test_cis_box_inside_no_cancel()
{
    std::cerr << "\n[TEST] test_cis_box_inside_no_cancel\n"
              << "  Source file:   check_if_spans.cpp\n"
              << "  Function:      check_if_spans -> check_if_spans_box\n"
              << "  Scenario:      100 nm cubic box, two one-molecule complexes at\n"
              << "                 x = -2 and x = +2 (well inside every wall).\n"
              << "  Pass criteria: cancelAssoc stays false and the temporary molecule\n"
              << "                 coordinates are unchanged (nothing to push back).\n";

    Parameters params = cis_make_params();
    Membrane membraneObject = cis_make_box_membrane(100.0);

    std::vector<Molecule> moleculeList {
        cis_make_molecule(Coord { -2.0, 0.0, 0.0 }, 0, 0),
        cis_make_molecule(Coord { 2.0, 0.0, 0.0 }, 1, 1)
    };
    Complex reactCom1 = cis_make_complex(0, Coord { -2.0, 0.0, 0.0 }, std::vector<int> { 0 });
    Complex reactCom2 = cis_make_complex(1, Coord { 2.0, 0.0, 0.0 }, std::vector<int> { 1 });

    bool cancelAssoc { false };

    std::cerr << "  Calling check_if_spans...\n";
    check_if_spans(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    // A small complex in the middle of a large box can never span it.
    EXPECT_FALSE(cancelAssoc)
        << "a tiny complex at the box center must not be flagged as spanning";

    // Nothing is out of bounds, so no repositioning should occur.
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpComCoord.x, -2.0)
        << "interior molecule 0 must not be translated";
    EXPECT_DOUBLE_EQ(moleculeList[1].tmpComCoord.x, 2.0)
        << "interior molecule 1 must not be translated";
    EXPECT_DOUBLE_EQ(moleculeList[0].tmpComCoord.y, 0.0)
        << "interior molecule 0 y-coordinate must not change";
    EXPECT_DOUBLE_EQ(moleculeList[1].tmpComCoord.z, 0.0)
        << "interior molecule 1 z-coordinate must not change";

    // Sanity check: everything is still inside the half-box in x.
    const double halfBox = membraneObject.waterBox.x / 2.0;
    EXPECT_LE(cis_max_tmp_x(moleculeList), halfBox)
        << "all molecules should remain inside the +x wall";
    EXPECT_GE(cis_min_tmp_x(moleculeList), -halfBox)
        << "all molecules should remain inside the -x wall";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc
              << ", tmp x range = [" << cis_min_tmp_x(moleculeList) << ", "
              << cis_max_tmp_x(moleculeList) << "] with walls at +/-" << halfBox << '\n';
}

// -----------------------------------------------------------------------------
// Test 4: Box boundary, the merged complex sticks out of BOTH x walls.
//
// Pass criteria: the complex cannot be fitted into the volume, so the
// association must be cancelled (cancelAssoc == true).
// -----------------------------------------------------------------------------
void test_cis_box_spanning_cancels_assoc()
{
    std::cerr << "\n[TEST] test_cis_box_spanning_cancels_assoc\n"
              << "  Source file:   check_if_spans.cpp\n"
              << "  Function:      check_if_spans -> check_if_spans_box\n"
              << "  Scenario:      10 nm cubic box (walls at +/-5) with molecules at\n"
              << "                 x = -30 and x = +30, i.e. the merged complex\n"
              << "                 protrudes through the -x AND the +x wall.\n"
              << "  Pass criteria: cancelAssoc is set to true, since the complex\n"
              << "                 spans the whole simulation volume.\n";

    Parameters params = cis_make_params();
    Membrane membraneObject = cis_make_box_membrane(10.0);

    std::vector<Molecule> moleculeList {
        cis_make_molecule(Coord { -30.0, 0.0, 0.0 }, 0, 0),
        cis_make_molecule(Coord { 30.0, 0.0, 0.0 }, 1, 1)
    };
    Complex reactCom1 = cis_make_complex(0, Coord { -30.0, 0.0, 0.0 }, std::vector<int> { 0 });
    Complex reactCom2 = cis_make_complex(1, Coord { 30.0, 0.0, 0.0 }, std::vector<int> { 1 });

    bool cancelAssoc { false };

    std::cerr << "  Calling check_if_spans...\n";
    check_if_spans(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_TRUE(cancelAssoc)
        << "a complex protruding through both x walls must cancel association";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc
              << " (expected true), tmp x range = [" << cis_min_tmp_x(moleculeList)
              << ", " << cis_max_tmp_x(moleculeList) << "]\n";
}

// -----------------------------------------------------------------------------
// Test 5: Box boundary, the merged complex only pokes out of the +x wall.
//
// A complex that fits in the box but is currently out of bounds on one side is
// supposed to be shifted back inside rather than rejected.
//
// Pass criteria: cancelAssoc stays false, and the molecules are certainly not
// pushed FURTHER outside than they started (monotonic, implementation-agnostic
// check).  The final extent is also reported for inspection.
// -----------------------------------------------------------------------------
void test_cis_box_out_one_side_not_cancelled()
{
    std::cerr << "\n[TEST] test_cis_box_out_one_side_not_cancelled\n"
              << "  Source file:   check_if_spans.cpp\n"
              << "  Function:      check_if_spans -> check_if_spans_box\n"
              << "  Scenario:      100 nm cubic box (walls at +/-50) with molecules at\n"
              << "                 x = +58 and x = +60, i.e. out of bounds on the +x\n"
              << "                 side only; the complex easily fits in the box.\n"
              << "  Pass criteria: cancelAssoc stays false (the complex does not span\n"
              << "                 the volume) and no molecule is moved further out\n"
              << "                 than where it started.\n";

    Parameters params = cis_make_params();
    Membrane membraneObject = cis_make_box_membrane(100.0);

    std::vector<Molecule> moleculeList {
        cis_make_molecule(Coord { 58.0, 0.0, 0.0 }, 0, 0),
        cis_make_molecule(Coord { 60.0, 0.0, 0.0 }, 1, 1)
    };
    Complex reactCom1 = cis_make_complex(0, Coord { 58.0, 0.0, 0.0 }, std::vector<int> { 0 });
    Complex reactCom2 = cis_make_complex(1, Coord { 60.0, 0.0, 0.0 }, std::vector<int> { 1 });

    const double initialMaxX = cis_max_tmp_x(moleculeList);
    const double halfBox = membraneObject.waterBox.x / 2.0;

    bool cancelAssoc { false };

    std::cerr << "  Calling check_if_spans (initial max x = " << initialMaxX
              << ", +x wall = " << halfBox << ")...\n";
    check_if_spans(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    const double finalMaxX = cis_max_tmp_x(moleculeList);

    // The complex is only 2 nm long, so it can never span a 100 nm box.
    EXPECT_FALSE(cancelAssoc)
        << "a small complex out of bounds on one side only must not be cancelled";

    // Whatever repositioning happens, it may not make the situation worse.
    EXPECT_LE(finalMaxX, initialMaxX + 1e-6)
        << "molecules must not be pushed further outside the +x wall";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc
              << ", max tmp x before = " << initialMaxX
              << ", after = " << finalMaxX << '\n';
}

// -----------------------------------------------------------------------------
// Test 6: Box boundary, multi-molecule complexes that span the y dimension.
//
// This confirms the dispatcher also works when each Complex owns more than one
// member molecule and the offending dimension is not x.
//
// Pass criteria: cancelAssoc is set to true.
// -----------------------------------------------------------------------------
void test_cis_box_multimolecule_spans_y()
{
    std::cerr << "\n[TEST] test_cis_box_multimolecule_spans_y\n"
              << "  Source file:   check_if_spans.cpp\n"
              << "  Function:      check_if_spans -> check_if_spans_box\n"
              << "  Scenario:      20 nm cubic box (walls at +/-10); complex 1 owns two\n"
              << "                 molecules near y = -40 and complex 2 owns two\n"
              << "                 molecules near y = +40, so the merged structure\n"
              << "                 protrudes through both y walls.\n"
              << "  Pass criteria: cancelAssoc is set to true.\n";

    Parameters params = cis_make_params();
    Membrane membraneObject = cis_make_box_membrane(20.0);

    // Two molecules per complex, all four sharing the same y extremes.
    std::vector<Molecule> moleculeList {
        cis_make_molecule(Coord { 0.0, -40.0, 0.0 }, 0, 0),
        cis_make_molecule(Coord { 1.0, -41.0, 0.0 }, 1, 0),
        cis_make_molecule(Coord { 0.0, 40.0, 0.0 }, 2, 1),
        cis_make_molecule(Coord { 1.0, 41.0, 0.0 }, 3, 1)
    };
    Complex reactCom1 = cis_make_complex(0, Coord { 0.5, -40.5, 0.0 }, std::vector<int> { 0, 1 });
    Complex reactCom2 = cis_make_complex(1, Coord { 0.5, 40.5, 0.0 }, std::vector<int> { 2, 3 });

    bool cancelAssoc { false };

    std::cerr << "  Calling check_if_spans...\n";
    check_if_spans(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_TRUE(cancelAssoc)
        << "a multi-molecule complex protruding through both y walls must cancel association";

    // The member lists must not have been altered by a pure geometry check.
    EXPECT_EQ(reactCom1.memberList.size(), 2u)
        << "complex 1 member list size must be preserved by the span check";
    EXPECT_EQ(reactCom2.memberList.size(), 2u)
        << "complex 2 member list size must be preserved by the span check";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc
              << " (expected true); member counts = " << reactCom1.memberList.size()
              << " and " << reactCom2.memberList.size() << '\n';
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each scenario runs in its own TEST so that a failure in
// one does not prevent the remaining scenarios from executing.
// -----------------------------------------------------------------------------
TEST(CheckIfSpans, SphereIsANoop) { test_cis_sphere_is_a_noop(); }
TEST(CheckIfSpans, SpherePreservesCancelFlag) { test_cis_sphere_preserves_cancel_flag(); }
TEST(CheckIfSpans, BoxInsideNoCancel) { test_cis_box_inside_no_cancel(); }
TEST(CheckIfSpans, BoxSpanningCancelsAssoc) { test_cis_box_spanning_cancels_assoc(); }
TEST(CheckIfSpans, BoxOutOneSideNotCancelled) { test_cis_box_out_one_side_not_cancelled(); }
TEST(CheckIfSpans, BoxMultiMoleculeSpansY) { test_cis_box_multimolecule_spans_y(); }