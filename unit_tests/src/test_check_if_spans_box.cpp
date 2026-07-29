/*! \file test_check_if_spans_box.cpp
 *
 * ### Unit test for ../src/boundary_conditions/check_if_spans_box.cpp
 *
 * Function under test:
 *
 *     void check_if_spans_box(bool& cancelAssoc, const Parameters& params,
 *                             Complex& reactCom1, Complex& reactCom2,
 *                             std::vector<Molecule>& moleculeList,
 *                             const Membrane& membraneObject)
 *
 * The routine is called right after two associating complexes have been placed
 * in contact (using the *temporary* association coordinates `tmpComCoord` and
 * `tmpICoords`). It then, independently for z, y and x:
 *
 *   1. Skips the dimension entirely unless
 *      (reactCom1.radius + reactCom2.radius) > waterBox.dim / 2
 *      (i.e. unless the merged complex is even *capable* of leaving the box).
 *   2. Finds the largest protrusion past the positive wall (`posWall`) and past
 *      the negative wall (`negWall`), scanning both complexes' member molecule
 *      centers of mass **and** all of their temporary interface coordinates.
 *   3. Cancels the association (cancelAssoc = true) if the structure sticks out
 *      of both walls at once, or if pushing it back in from one wall would make
 *      it poke out of the opposite wall (posWall + negWall > 0).
 *   4. Otherwise translates *both* complexes (and every member molecule and
 *      interface) so that the outermost point sits exactly on the wall.
 *
 * All the tests below verify these four behaviours and print what is being
 * checked so the console log is self-describing.
 */

#include "boundary_conditions/reflect_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Small helpers used to build a minimal but valid two-complex system.
// -----------------------------------------------------------------------------

/*! \brief Build a Molecule that only carries association (temporary) coords.
 *
 * `check_if_spans_box` reads/writes exclusively `tmpComCoord` and `tmpICoords`,
 * so those are the fields we populate. `comCoord` is set to the same value so
 * we can later assert that the *real* coordinates were left untouched.
 *
 * \param[in] tmpCom   temporary center-of-mass coordinate
 * \param[in] ifaces   temporary interface coordinates
 * \param[in] comIndex index of the parent complex
 */
Molecule cisb_make_molecule(const Coord& tmpCom, const std::vector<Coord>& ifaces, int comIndex)
{
    Molecule mol;
    mol.myComIndex = comIndex;
    mol.comCoord = tmpCom; // real coordinate, must not be modified
    mol.tmpComCoord = tmpCom; // temporary coordinate, may be translated
    mol.tmpICoords = ifaces; // temporary interface coordinates
    return mol;
}

/*! \brief Build a Complex with a bounding radius and a single member molecule.
 *
 * \param[in] tmpCom    temporary center-of-mass coordinate of the complex
 * \param[in] radius    bounding-sphere radius (drives the "canBeOutside" test)
 * \param[in] memberIdx index of its member molecule in moleculeList
 */
Complex cisb_make_complex(const Coord& tmpCom, double radius, int memberIdx)
{
    Complex com;
    com.comCoord = tmpCom;
    com.tmpComCoord = tmpCom;
    com.radius = radius;
    com.memberList.clear();
    com.memberList.push_back(memberIdx);
    return com;
}

/*! \brief Build a rectangular (non-spherical) water box of the given size. */
Membrane cisb_make_box(double x, double y, double z)
{
    Membrane membraneObject;
    membraneObject.isSphere = false;
    membraneObject.isBox = true;
    membraneObject.waterBox.x = x;
    membraneObject.waterBox.y = y;
    membraneObject.waterBox.z = z;
    membraneObject.waterBox.volume = x * y * z;
    return membraneObject;
}

// Tolerance used for all floating point comparisons of coordinates.
const double kCisbTol = 1e-9;

} // namespace

// -----------------------------------------------------------------------------
// Test 1: When the two complexes are small (radius1 + radius2 <= boxDim/2), the
//         dimension is never examined -- even coordinates outside the box are
//         left alone and association is not cancelled.
// -----------------------------------------------------------------------------
void test_cisb_no_check_when_radii_small()
{
    std::cerr << "\n[TEST] test_cisb_no_check_when_radii_small\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box\n"
              << "  Scenario:      radius1 + radius2 (2.0) <= boxDim/2 (50) in every\n"
              << "                 dimension, so no dimension is inspected at all.\n"
              << "  Pass criteria: cancelAssoc stays false and no coordinate moves,\n"
              << "                 even though the molecules lie outside the box.\n";

    Parameters params; // unused by the function, but part of the signature
    Membrane membraneObject = cisb_make_box(100.0, 100.0, 100.0);

    // Molecules deliberately placed well outside the +/-50 walls.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ 90.0, 90.0, 90.0 }, { Coord{ 91.0, 91.0, 91.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ -90.0, -90.0, -90.0 }, { Coord{ -91.0, -91.0, -91.0 } }, 1));

    // Tiny radii: 1.0 + 1.0 = 2.0, far below 100/2 = 50 in every dimension.
    Complex reactCom1 = cisb_make_complex(Coord{ 90.0, 90.0, 90.0 }, 1.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ -90.0, -90.0, -90.0 }, 1.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    // Nothing at all should have happened.
    EXPECT_FALSE(cancelAssoc) << "Association must not be cancelled when the size test is not triggered";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, 90.0, kCisbTol) << "mol0 x must be untouched";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.y, 90.0, kCisbTol) << "mol0 y must be untouched";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.z, 90.0, kCisbTol) << "mol0 z must be untouched";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.x, -90.0, kCisbTol) << "mol1 x must be untouched";
    EXPECT_NEAR(moleculeList[0].tmpICoords[0].x, 91.0, kCisbTol) << "mol0 interface x must be untouched";
    EXPECT_NEAR(reactCom1.tmpComCoord.x, 90.0, kCisbTol) << "complex1 tmpComCoord must be untouched";
    EXPECT_NEAR(reactCom2.tmpComCoord.x, -90.0, kCisbTol) << "complex2 tmpComCoord must be untouched";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc
              << ", mol0 tmp x = " << moleculeList[0].tmpComCoord.x << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: Protrusion past the +X wall only -> everything is shifted back by
//         exactly the protrusion so the outermost point lands on the wall.
// -----------------------------------------------------------------------------
void test_cisb_translate_positive_x()
{
    std::cerr << "\n[TEST] test_cisb_translate_positive_x\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (x branch, +X protrusion)\n"
              << "  Scenario:      box x=100 (+X wall at 50), radii 30+30 = 60 > 50, so\n"
              << "                 the x dimension is inspected. Farthest point x = 55.\n"
              << "  Pass criteria: every temporary coordinate shifts by -5 (55 -> 50)\n"
              << "                 and cancelAssoc stays false.\n";

    Parameters params;
    Membrane membraneObject = cisb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ 55.0, 0.0, 0.0 }, { Coord{ 54.0, 0.0, 0.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ 45.0, 0.0, 0.0 }, { Coord{ 46.0, 0.0, 0.0 } }, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ 55.0, 0.0, 0.0 }, 30.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ 45.0, 0.0, 0.0 }, 30.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    // Protrusion is 55 - 50 = 5, so the whole structure shifts by -5 in x.
    EXPECT_FALSE(cancelAssoc) << "One-sided protrusion should be fixed by translation, not cancelled";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, 50.0, kCisbTol) << "mol0 COM should sit exactly on the +X wall";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.x, 40.0, kCisbTol) << "mol1 COM should shift by the same -5";
    EXPECT_NEAR(moleculeList[0].tmpICoords[0].x, 49.0, kCisbTol) << "mol0 interface should shift by -5";
    EXPECT_NEAR(moleculeList[1].tmpICoords[0].x, 41.0, kCisbTol) << "mol1 interface should shift by -5";
    EXPECT_NEAR(reactCom1.tmpComCoord.x, 50.0, kCisbTol) << "complex1 tmpComCoord should shift by -5";
    EXPECT_NEAR(reactCom2.tmpComCoord.x, 40.0, kCisbTol) << "complex2 tmpComCoord should shift by -5";

    // The permanent coordinates must never be modified by this routine.
    EXPECT_NEAR(moleculeList[0].comCoord.x, 55.0, kCisbTol) << "real comCoord must be left alone";
    EXPECT_NEAR(moleculeList[1].comCoord.x, 45.0, kCisbTol) << "real comCoord must be left alone";

    std::cerr << "  Shifted x positions: mol0 = " << moleculeList[0].tmpComCoord.x
              << ", mol1 = " << moleculeList[1].tmpComCoord.x << '\n';
}

// -----------------------------------------------------------------------------
// Test 3: Protrusion past the -X wall only -> shift in the +x direction.
// -----------------------------------------------------------------------------
void test_cisb_translate_negative_x()
{
    std::cerr << "\n[TEST] test_cisb_translate_negative_x\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (x branch, -X protrusion)\n"
              << "  Scenario:      box x=100 (-X wall at -50), radii 30+30 > 50.\n"
              << "                 Farthest point x = -55, i.e. 5 nm past the wall.\n"
              << "  Pass criteria: every temporary coordinate shifts by +5 (-55 -> -50).\n";

    Parameters params;
    Membrane membraneObject = cisb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ -55.0, 0.0, 0.0 }, { Coord{ -54.0, 0.0, 0.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ -45.0, 0.0, 0.0 }, { Coord{ -46.0, 0.0, 0.0 } }, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ -55.0, 0.0, 0.0 }, 30.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ -45.0, 0.0, 0.0 }, 30.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_FALSE(cancelAssoc) << "One-sided protrusion should be fixed by translation, not cancelled";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, -50.0, kCisbTol) << "mol0 COM should sit exactly on the -X wall";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.x, -40.0, kCisbTol) << "mol1 COM should shift by the same +5";
    EXPECT_NEAR(moleculeList[0].tmpICoords[0].x, -49.0, kCisbTol) << "mol0 interface should shift by +5";
    EXPECT_NEAR(reactCom1.tmpComCoord.x, -50.0, kCisbTol) << "complex1 tmpComCoord should shift by +5";
    EXPECT_NEAR(reactCom2.tmpComCoord.x, -40.0, kCisbTol) << "complex2 tmpComCoord should shift by +5";

    std::cerr << "  Shifted x positions: mol0 = " << moleculeList[0].tmpComCoord.x
              << ", mol1 = " << moleculeList[1].tmpComCoord.x << '\n';
}

// -----------------------------------------------------------------------------
// Test 4: Structure pokes out of BOTH x walls -> association is cancelled and
//         no translation is applied.
// -----------------------------------------------------------------------------
void test_cisb_cancel_when_spanning_x()
{
    std::cerr << "\n[TEST] test_cisb_cancel_when_spanning_x\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (x branch, spans the box)\n"
              << "  Scenario:      narrow box x=20 (walls at +/-10), radii 20+20 > 10.\n"
              << "                 Points at x = +15 and x = -15 poke out of both walls.\n"
              << "  Pass criteria: cancelAssoc becomes true and coordinates are NOT moved.\n";

    Parameters params;
    // y and z are big enough (radii sum 40 <= 50) that only x is inspected.
    Membrane membraneObject = cisb_make_box(20.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ 15.0, 0.0, 0.0 }, { Coord{ 14.0, 0.0, 0.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ -15.0, 0.0, 0.0 }, { Coord{ -14.0, 0.0, 0.0 } }, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ 15.0, 0.0, 0.0 }, 20.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ -15.0, 0.0, 0.0 }, 20.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_TRUE(cancelAssoc) << "Association must be cancelled when the structure spans the box in x";
    // With both walls violated neither translation branch runs.
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, 15.0, kCisbTol) << "no translation should occur when cancelling";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.x, -15.0, kCisbTol) << "no translation should occur when cancelling";
    EXPECT_NEAR(reactCom1.tmpComCoord.x, 15.0, kCisbTol) << "complex1 tmpComCoord should be unchanged";
    EXPECT_NEAR(reactCom2.tmpComCoord.x, -15.0, kCisbTol) << "complex2 tmpComCoord should be unchanged";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc << '\n';
}

// -----------------------------------------------------------------------------
// Test 5: Only the +X wall is violated, but the structure is so wide that
//         pushing it back in would make it exit the -X wall
//         (posWall + negWall > 0) -> association is cancelled.
// -----------------------------------------------------------------------------
void test_cisb_cancel_when_pushback_would_span()
{
    std::cerr << "\n[TEST] test_cisb_cancel_when_pushback_would_span\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (x branch, posWall + negWall > 0)\n"
              << "  Scenario:      box x=4 (walls at +/-2), radii 1.5+1.5 = 3 > 2.\n"
              << "                 Points at x = 3 (posWall = 1) and x = -1.5\n"
              << "                 (negWall = -0.5); sum = +0.5 > 0, so putting it back\n"
              << "                 inside would push the other end out of -X.\n"
              << "  Pass criteria: cancelAssoc becomes true.\n";

    Parameters params;
    Membrane membraneObject = cisb_make_box(4.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ 3.0, 0.0, 0.0 }, {}, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ -1.5, 0.0, 0.0 }, {}, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ 3.0, 0.0, 0.0 }, 1.5, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ -1.5, 0.0, 0.0 }, 1.5, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_TRUE(cancelAssoc) << "Cancel expected because pushing back in would exit the opposite wall";

    // Documented side effect: the one-sided translation branch still executes
    // (the code does not return early), so coordinates are also shifted by
    // posWall = 1.0. Verify the documented behaviour so a change is noticed.
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, 2.0, kCisbTol)
        << "the +X translation branch still runs after cancelAssoc is set";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.x, -2.5, kCisbTol)
        << "the +X translation branch still runs after cancelAssoc is set";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc
              << ", mol0 tmp x = " << moleculeList[0].tmpComCoord.x << '\n';
}

// -----------------------------------------------------------------------------
// Test 6: The y branch works the same way as the x branch.
// -----------------------------------------------------------------------------
void test_cisb_translate_y()
{
    std::cerr << "\n[TEST] test_cisb_translate_y\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (y branch)\n"
              << "  Scenario:      box y=100 (+Y wall at 50), radii 30+30 > 50, farthest\n"
              << "                 point y = 57.\n"
              << "  Pass criteria: all y coordinates shift by -7; x and z untouched.\n";

    Parameters params;
    Membrane membraneObject = cisb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ 0.0, 57.0, 0.0 }, { Coord{ 0.0, 56.0, 0.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ 0.0, 47.0, 0.0 }, { Coord{ 0.0, 48.0, 0.0 } }, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ 0.0, 57.0, 0.0 }, 30.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ 0.0, 47.0, 0.0 }, 30.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_FALSE(cancelAssoc) << "One-sided +Y protrusion should be translated, not cancelled";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.y, 50.0, kCisbTol) << "mol0 y should land on the +Y wall";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.y, 40.0, kCisbTol) << "mol1 y should shift by the same -7";
    EXPECT_NEAR(moleculeList[0].tmpICoords[0].y, 49.0, kCisbTol) << "mol0 interface y should shift by -7";
    EXPECT_NEAR(reactCom1.tmpComCoord.y, 50.0, kCisbTol) << "complex1 tmpComCoord.y should shift by -7";
    // x and z were inside the box, so they must be untouched.
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, 0.0, kCisbTol) << "x must not change";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.z, 0.0, kCisbTol) << "z must not change";

    std::cerr << "  Shifted y positions: mol0 = " << moleculeList[0].tmpComCoord.y
              << ", mol1 = " << moleculeList[1].tmpComCoord.y << '\n';
}

// -----------------------------------------------------------------------------
// Test 7: The z branch (evaluated first inside the function) works too, here
//         for a protrusion past the -Z wall.
// -----------------------------------------------------------------------------
void test_cisb_translate_z()
{
    std::cerr << "\n[TEST] test_cisb_translate_z\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (z branch, -Z protrusion)\n"
              << "  Scenario:      box z=100 (-Z wall at -50), radii 30+30 > 50, farthest\n"
              << "                 point z = -56.\n"
              << "  Pass criteria: all z coordinates shift by +6 (-56 -> -50).\n";

    Parameters params;
    Membrane membraneObject = cisb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ 0.0, 0.0, -56.0 }, { Coord{ 0.0, 0.0, -55.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ 0.0, 0.0, -46.0 }, { Coord{ 0.0, 0.0, -47.0 } }, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ 0.0, 0.0, -56.0 }, 30.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ 0.0, 0.0, -46.0 }, 30.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_FALSE(cancelAssoc) << "One-sided -Z protrusion should be translated, not cancelled";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.z, -50.0, kCisbTol) << "mol0 z should land on the -Z wall";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.z, -40.0, kCisbTol) << "mol1 z should shift by the same +6";
    EXPECT_NEAR(moleculeList[0].tmpICoords[0].z, -49.0, kCisbTol) << "mol0 interface z should shift by +6";
    EXPECT_NEAR(reactCom1.tmpComCoord.z, -50.0, kCisbTol) << "complex1 tmpComCoord.z should shift by +6";
    EXPECT_NEAR(reactCom2.tmpComCoord.z, -40.0, kCisbTol) << "complex2 tmpComCoord.z should shift by +6";

    std::cerr << "  Shifted z positions: mol0 = " << moleculeList[0].tmpComCoord.z
              << ", mol1 = " << moleculeList[1].tmpComCoord.z << '\n';
}

// -----------------------------------------------------------------------------
// Test 8: The scan includes interface coordinates, not just centers of mass --
//         a molecule whose COM is inside but whose interface is outside must
//         still drive the translation.
// -----------------------------------------------------------------------------
void test_cisb_interface_determines_shift()
{
    std::cerr << "\n[TEST] test_cisb_interface_determines_shift\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (interface scan)\n"
              << "  Scenario:      all centers of mass are inside the box, but one\n"
              << "                 interface sits at x = 58 (8 nm past the +X wall).\n"
              << "  Pass criteria: the shift is driven by the interface: -8 in x, so the\n"
              << "                 offending interface ends up exactly at x = 50.\n";

    Parameters params;
    Membrane membraneObject = cisb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    // COM at 45 (inside) but interface at 58 (outside by 8).
    moleculeList.push_back(cisb_make_molecule(Coord{ 45.0, 0.0, 0.0 },
        { Coord{ 58.0, 0.0, 0.0 }, Coord{ 40.0, 0.0, 0.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ 30.0, 0.0, 0.0 }, { Coord{ 31.0, 0.0, 0.0 } }, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ 45.0, 0.0, 0.0 }, 30.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ 30.0, 0.0, 0.0 }, 30.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_FALSE(cancelAssoc) << "A single-sided interface protrusion should be translated";
    EXPECT_NEAR(moleculeList[0].tmpICoords[0].x, 50.0, kCisbTol)
        << "the offending interface should end up exactly on the +X wall";
    EXPECT_NEAR(moleculeList[0].tmpICoords[1].x, 32.0, kCisbTol) << "the other interface shifts by -8 too";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, 37.0, kCisbTol) << "mol0 COM shifts by -8 (45 -> 37)";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.x, 22.0, kCisbTol) << "mol1 COM shifts by -8 (30 -> 22)";
    EXPECT_NEAR(moleculeList[1].tmpICoords[0].x, 23.0, kCisbTol) << "mol1 interface shifts by -8 (31 -> 23)";

    std::cerr << "  Interface after shift = " << moleculeList[0].tmpICoords[0].x
              << ", mol0 COM = " << moleculeList[0].tmpComCoord.x << '\n';
}

// -----------------------------------------------------------------------------
// Test 9: Two dimensions can be corrected in a single call; the z branch runs
//         before the x branch and each is independent.
// -----------------------------------------------------------------------------
void test_cisb_multiple_dimensions_corrected()
{
    std::cerr << "\n[TEST] test_cisb_multiple_dimensions_corrected\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (z branch then x branch)\n"
              << "  Scenario:      box 100^3, radii 30+30 > 50 in every dimension.\n"
              << "                 Structure pokes 7 nm past -Z and 5 nm past +X.\n"
              << "  Pass criteria: z shifts by +7 and x shifts by -5 in the same call;\n"
              << "                 y (all zeros) is untouched and cancelAssoc is false.\n";

    Parameters params;
    Membrane membraneObject = cisb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ 55.0, 0.0, -57.0 }, { Coord{ 55.0, 0.0, -57.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ 45.0, 0.0, -40.0 }, { Coord{ 45.0, 0.0, -40.0 } }, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ 55.0, 0.0, -57.0 }, 30.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ 45.0, 0.0, -40.0 }, 30.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_FALSE(cancelAssoc) << "Two independent one-sided protrusions should both be translated";

    // z: farthest point -57, so +7 shift.
    EXPECT_NEAR(moleculeList[0].tmpComCoord.z, -50.0, kCisbTol) << "mol0 z should land on the -Z wall";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.z, -33.0, kCisbTol) << "mol1 z should shift by +7";
    // x: farthest point 55, so -5 shift.
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, 50.0, kCisbTol) << "mol0 x should land on the +X wall";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.x, 40.0, kCisbTol) << "mol1 x should shift by -5";
    // interfaces follow their molecules.
    EXPECT_NEAR(moleculeList[0].tmpICoords[0].x, 50.0, kCisbTol) << "mol0 interface x should shift by -5";
    EXPECT_NEAR(moleculeList[0].tmpICoords[0].z, -50.0, kCisbTol) << "mol0 interface z should shift by +7";
    // y never left the box.
    EXPECT_NEAR(moleculeList[0].tmpComCoord.y, 0.0, kCisbTol) << "y must remain unchanged";
    EXPECT_NEAR(moleculeList[1].tmpComCoord.y, 0.0, kCisbTol) << "y must remain unchanged";
    // complexes tracked both shifts.
    EXPECT_NEAR(reactCom1.tmpComCoord.x, 50.0, kCisbTol) << "complex1 x should shift by -5";
    EXPECT_NEAR(reactCom1.tmpComCoord.z, -50.0, kCisbTol) << "complex1 z should shift by +7";
    EXPECT_NEAR(reactCom2.tmpComCoord.x, 40.0, kCisbTol) << "complex2 x should shift by -5";
    EXPECT_NEAR(reactCom2.tmpComCoord.z, -33.0, kCisbTol) << "complex2 z should shift by +7";

    std::cerr << "  mol0 final tmp coords = (" << moleculeList[0].tmpComCoord.x << ", "
              << moleculeList[0].tmpComCoord.y << ", " << moleculeList[0].tmpComCoord.z << ")\n"
              << "  mol1 final tmp coords = (" << moleculeList[1].tmpComCoord.x << ", "
              << moleculeList[1].tmpComCoord.y << ", " << moleculeList[1].tmpComCoord.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 10: A complex that is entirely inside the box is untouched even when the
//          size test *is* triggered (canBeOutside == true).
// -----------------------------------------------------------------------------
void test_cisb_inside_box_untouched()
{
    std::cerr << "\n[TEST] test_cisb_inside_box_untouched\n"
              << "  Source file:   check_if_spans_box.cpp\n"
              << "  Function:      check_if_spans_box (all branches, nothing outside)\n"
              << "  Scenario:      radii 30+30 > 50 so every dimension IS inspected, but\n"
              << "                 all coordinates are comfortably inside +/-50.\n"
              << "  Pass criteria: cancelAssoc false and no coordinate changes.\n";

    Parameters params;
    Membrane membraneObject = cisb_make_box(100.0, 100.0, 100.0);

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(cisb_make_molecule(Coord{ 5.0, -3.0, 2.0 }, { Coord{ 6.0, -4.0, 3.0 } }, 0));
    moleculeList.push_back(cisb_make_molecule(Coord{ -5.0, 3.0, -2.0 }, { Coord{ -6.0, 4.0, -3.0 } }, 1));

    Complex reactCom1 = cisb_make_complex(Coord{ 5.0, -3.0, 2.0 }, 30.0, 0);
    Complex reactCom2 = cisb_make_complex(Coord{ -5.0, 3.0, -2.0 }, 30.0, 1);

    bool cancelAssoc { false };
    std::cerr << "  Calling check_if_spans_box...\n";
    check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);

    EXPECT_FALSE(cancelAssoc) << "A fully interior structure must never cancel association";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.x, 5.0, kCisbTol) << "interior x must not move";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.y, -3.0, kCisbTol) << "interior y must not move";
    EXPECT_NEAR(moleculeList[0].tmpComCoord.z, 2.0, kCisbTol) << "interior z must not move";
    EXPECT_NEAR(moleculeList[1].tmpICoords[0].x, -6.0, kCisbTol) << "interior interface must not move";
    EXPECT_NEAR(reactCom1.tmpComCoord.z, 2.0, kCisbTol) << "complex1 tmpComCoord must not move";
    EXPECT_NEAR(reactCom2.tmpComCoord.z, -2.0, kCisbTol) << "complex2 tmpComCoord must not move";

    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc << " (no translation applied)\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each scenario is its own TEST so a failure in one does
// not stop the others from running (only non-fatal EXPECT_* macros are used).
// -----------------------------------------------------------------------------
TEST(CheckIfSpansBox, NoCheckWhenRadiiSmall) { test_cisb_no_check_when_radii_small(); }
TEST(CheckIfSpansBox, TranslatePositiveX) { test_cisb_translate_positive_x(); }
TEST(CheckIfSpansBox, TranslateNegativeX) { test_cisb_translate_negative_x(); }
TEST(CheckIfSpansBox, CancelWhenSpanningX) { test_cisb_cancel_when_spanning_x(); }
TEST(CheckIfSpansBox, CancelWhenPushbackWouldSpan) { test_cisb_cancel_when_pushback_would_span(); }
TEST(CheckIfSpansBox, TranslateY) { test_cisb_translate_y(); }
TEST(CheckIfSpansBox, TranslateZ) { test_cisb_translate_z(); }
TEST(CheckIfSpansBox, InterfaceDeterminesShift) { test_cisb_interface_determines_shift(); }
TEST(CheckIfSpansBox, MultipleDimensionsCorrected) { test_cisb_multiple_dimensions_corrected(); }
TEST(CheckIfSpansBox, InsideBoxUntouched) { test_cisb_inside_box_untouched(); }