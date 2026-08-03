/*! \file test_subtract_off_com_position.cpp
 *
 * ### Unit test for src/reactions/subtract_off_com_position.cpp
 *
 * The function under test is:
 *
 *     void subtract_off_com_position(Molecule& base1);
 *
 * Its job is very simple but important for association bookkeeping:
 *   1. Every interface coordinate of the molecule is shifted by the molecule's
 *      center-of-mass (COM) coordinate, so that interface coordinates become
 *      *relative* to the COM.
 *   2. The COM coordinate itself is then set to exactly (0, 0, 0).
 *
 * The function does not appear in any public header, so we forward declare it
 * here with C++ linkage exactly as it is defined in the source file.
 *
 * Every test below prints what is being exercised and what the pass criteria
 * are, and uses only non-fatal EXPECT_* assertions so the whole suite runs
 * even if an individual check fails.
 */

#include <iostream>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Forward declaration of the function under test (defined in
// src/reactions/subtract_off_com_position.cpp, not exposed via a header).
// -----------------------------------------------------------------------------
void subtract_off_com_position(Molecule& base1);

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace so they cannot collide with other TUs).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a molecule with a given COM and an arbitrary list of
 *         absolute interface coordinates.
 *
 * \param[in] com          Absolute center-of-mass coordinate.
 * \param[in] ifaceCoords  Absolute coordinates for each interface.
 * \return A Molecule populated with the requested geometry.
 */
Molecule soc_make_molecule(const Coord& com, const std::vector<Coord>& ifaceCoords)
{
    Molecule mol;
    mol.comCoord = com;
    mol.interfaceList.clear();

    // Add one Iface per requested coordinate. Only the coordinate matters for
    // this function, but we set relIndex so output is easier to interpret.
    for (std::size_t i = 0; i < ifaceCoords.size(); ++i) {
        Molecule::Iface iface;
        iface.coord = ifaceCoords[i];
        iface.relIndex = static_cast<int>(i);
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Print the geometry of a molecule to stderr for traceability. */
void soc_dump_molecule(const Molecule& mol, const std::string& label)
{
    std::cerr << "    " << label << ": COM = (" << mol.comCoord.x << ", "
              << mol.comCoord.y << ", " << mol.comCoord.z << ")\n";
    for (std::size_t i = 0; i < mol.interfaceList.size(); ++i) {
        std::cerr << "      iface[" << i << "] = ("
                  << mol.interfaceList[i].coord.x << ", "
                  << mol.interfaceList[i].coord.y << ", "
                  << mol.interfaceList[i].coord.z << ")\n";
    }
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Typical case -- non-zero COM with several interfaces.
//         Each interface must be shifted by exactly -COM, and the COM zeroed.
// -----------------------------------------------------------------------------
void test_soc_shifts_interfaces_and_zeroes_com()
{
    std::cerr << "\n[TEST] test_soc_shifts_interfaces_and_zeroes_com\n"
              << "  Source file:   src/reactions/subtract_off_com_position.cpp\n"
              << "  Function:      subtract_off_com_position(Molecule&)\n"
              << "  Scenario:      molecule with COM (10, -5, 2) and three\n"
              << "                 absolute interface coordinates.\n"
              << "  Pass criteria: each interface becomes (absolute - COM) and\n"
              << "                 the COM becomes exactly (0, 0, 0).\n";

    const Coord com { 10.0, -5.0, 2.0 };
    const std::vector<Coord> absIfaces {
        Coord { 11.0, -5.0, 2.0 }, // +1 in x relative to COM
        Coord { 10.0, -3.0, 2.0 }, // +2 in y relative to COM
        Coord { 7.5, -5.0, -0.5 } // -2.5 in x, -2.5 in z relative to COM
    };

    Molecule mol = soc_make_molecule(com, absIfaces);
    soc_dump_molecule(mol, "before");

    std::cerr << "  Calling subtract_off_com_position...\n";
    subtract_off_com_position(mol);
    soc_dump_molecule(mol, "after ");

    // The COM must be exactly zero after the call.
    EXPECT_DOUBLE_EQ(mol.comCoord.x, 0.0) << "COM x should be zeroed";
    EXPECT_DOUBLE_EQ(mol.comCoord.y, 0.0) << "COM y should be zeroed";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, 0.0) << "COM z should be zeroed";

    // The interface count must be unchanged.
    EXPECT_EQ(mol.interfaceList.size(), absIfaces.size())
        << "Number of interfaces must not change";

    // Every interface must equal (originalAbsolute - originalCOM).
    for (std::size_t i = 0; i < absIfaces.size() && i < mol.interfaceList.size(); ++i) {
        const double expX = absIfaces[i].x - com.x;
        const double expY = absIfaces[i].y - com.y;
        const double expZ = absIfaces[i].z - com.z;
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.x, expX)
            << "iface[" << i << "].x should be absolute-x minus COM-x";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.y, expY)
            << "iface[" << i << "].y should be absolute-y minus COM-y";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.z, expZ)
            << "iface[" << i << "].z should be absolute-z minus COM-z";
    }
}

// -----------------------------------------------------------------------------
// Test 2: COM already at the origin -- interfaces must be untouched.
// -----------------------------------------------------------------------------
void test_soc_com_already_at_origin_is_noop()
{
    std::cerr << "\n[TEST] test_soc_com_already_at_origin_is_noop\n"
              << "  Source file:   src/reactions/subtract_off_com_position.cpp\n"
              << "  Function:      subtract_off_com_position(Molecule&)\n"
              << "  Scenario:      COM is already (0, 0, 0).\n"
              << "  Pass criteria: interface coordinates are bit-for-bit\n"
              << "                 unchanged and the COM stays at the origin.\n";

    const std::vector<Coord> absIfaces {
        Coord { 1.0, 2.0, 3.0 },
        Coord { -4.0, -5.0, -6.0 }
    };

    Molecule mol = soc_make_molecule(Coord { 0.0, 0.0, 0.0 }, absIfaces);
    soc_dump_molecule(mol, "before");

    std::cerr << "  Calling subtract_off_com_position...\n";
    subtract_off_com_position(mol);
    soc_dump_molecule(mol, "after ");

    // Subtracting zero must leave the interfaces exactly as they were.
    for (std::size_t i = 0; i < absIfaces.size() && i < mol.interfaceList.size(); ++i) {
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.x, absIfaces[i].x)
            << "iface[" << i << "].x should be unchanged when COM is the origin";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.y, absIfaces[i].y)
            << "iface[" << i << "].y should be unchanged when COM is the origin";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.z, absIfaces[i].z)
            << "iface[" << i << "].z should be unchanged when COM is the origin";
    }

    EXPECT_DOUBLE_EQ(mol.comCoord.x, 0.0) << "COM x should remain 0";
    EXPECT_DOUBLE_EQ(mol.comCoord.y, 0.0) << "COM y should remain 0";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, 0.0) << "COM z should remain 0";
}

// -----------------------------------------------------------------------------
// Test 3: Molecule with no interfaces -- the loop body never executes but the
//         COM must still be reset to the origin without crashing.
// -----------------------------------------------------------------------------
void test_soc_empty_interface_list()
{
    std::cerr << "\n[TEST] test_soc_empty_interface_list\n"
              << "  Source file:   src/reactions/subtract_off_com_position.cpp\n"
              << "  Function:      subtract_off_com_position(Molecule&)\n"
              << "  Scenario:      molecule has a non-zero COM but zero interfaces.\n"
              << "  Pass criteria: no crash, interface list stays empty, and the\n"
              << "                 COM is set to (0, 0, 0).\n";

    Molecule mol = soc_make_molecule(Coord { 3.0, -7.0, 11.0 }, {});
    soc_dump_molecule(mol, "before");

    std::cerr << "  Calling subtract_off_com_position...\n";
    subtract_off_com_position(mol);
    soc_dump_molecule(mol, "after ");

    EXPECT_TRUE(mol.interfaceList.empty()) << "Interface list must remain empty";
    EXPECT_DOUBLE_EQ(mol.comCoord.x, 0.0) << "COM x should be zeroed";
    EXPECT_DOUBLE_EQ(mol.comCoord.y, 0.0) << "COM y should be zeroed";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, 0.0) << "COM z should be zeroed";
}

// -----------------------------------------------------------------------------
// Test 4: Single interface coincident with the COM -- it must land exactly on
//         the origin after the shift (a common case for point molecules).
// -----------------------------------------------------------------------------
void test_soc_interface_coincident_with_com()
{
    std::cerr << "\n[TEST] test_soc_interface_coincident_with_com\n"
              << "  Source file:   src/reactions/subtract_off_com_position.cpp\n"
              << "  Function:      subtract_off_com_position(Molecule&)\n"
              << "  Scenario:      one interface located exactly at the COM\n"
              << "                 (the 'point molecule' case).\n"
              << "  Pass criteria: the interface ends up at (0, 0, 0), just\n"
              << "                 like the COM.\n";

    const Coord com { -12.25, 4.75, 0.5 };
    Molecule mol = soc_make_molecule(com, { com });
    soc_dump_molecule(mol, "before");

    std::cerr << "  Calling subtract_off_com_position...\n";
    subtract_off_com_position(mol);
    soc_dump_molecule(mol, "after ");

    ASSERT_FALSE(mol.interfaceList.empty()) << "Test setup should provide one interface";
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.x, 0.0)
        << "Coincident interface x should collapse onto the origin";
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.y, 0.0)
        << "Coincident interface y should collapse onto the origin";
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.z, 0.0)
        << "Coincident interface z should collapse onto the origin";
    EXPECT_DOUBLE_EQ(mol.comCoord.x, 0.0) << "COM x should be zeroed";
    EXPECT_DOUBLE_EQ(mol.comCoord.y, 0.0) << "COM y should be zeroed";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, 0.0) << "COM z should be zeroed";
}

// -----------------------------------------------------------------------------
// Test 5: Internal geometry preservation -- the pairwise distances between the
//         interfaces must not change, since the operation is a pure
//         translation.
// -----------------------------------------------------------------------------
void test_soc_preserves_relative_geometry()
{
    std::cerr << "\n[TEST] test_soc_preserves_relative_geometry\n"
              << "  Source file:   src/reactions/subtract_off_com_position.cpp\n"
              << "  Function:      subtract_off_com_position(Molecule&)\n"
              << "  Scenario:      three interfaces around an off-origin COM.\n"
              << "  Pass criteria: every interface-interface separation vector\n"
              << "                 is identical before and after the call\n"
              << "                 (a translation is rigid).\n";

    const Coord com { 100.0, 200.0, -300.0 };
    const std::vector<Coord> absIfaces {
        Coord { 101.0, 200.0, -300.0 },
        Coord { 100.0, 203.0, -300.0 },
        Coord { 100.0, 200.0, -295.0 }
    };

    Molecule mol = soc_make_molecule(com, absIfaces);
    soc_dump_molecule(mol, "before");

    std::cerr << "  Calling subtract_off_com_position...\n";
    subtract_off_com_position(mol);
    soc_dump_molecule(mol, "after ");

    // Compare all pairwise separation vectors before/after the translation.
    for (std::size_t i = 0; i < absIfaces.size(); ++i) {
        for (std::size_t j = i + 1; j < absIfaces.size(); ++j) {
            if (i >= mol.interfaceList.size() || j >= mol.interfaceList.size())
                continue;
            const double beforeDx = absIfaces[j].x - absIfaces[i].x;
            const double beforeDy = absIfaces[j].y - absIfaces[i].y;
            const double beforeDz = absIfaces[j].z - absIfaces[i].z;
            const double afterDx = mol.interfaceList[j].coord.x - mol.interfaceList[i].coord.x;
            const double afterDy = mol.interfaceList[j].coord.y - mol.interfaceList[i].coord.y;
            const double afterDz = mol.interfaceList[j].coord.z - mol.interfaceList[i].coord.z;

            EXPECT_DOUBLE_EQ(afterDx, beforeDx)
                << "x separation of ifaces " << i << "," << j << " must be preserved";
            EXPECT_DOUBLE_EQ(afterDy, beforeDy)
                << "y separation of ifaces " << i << "," << j << " must be preserved";
            EXPECT_DOUBLE_EQ(afterDz, beforeDz)
                << "z separation of ifaces " << i << "," << j << " must be preserved";
        }
    }
}

// -----------------------------------------------------------------------------
// Test 6: Idempotency -- calling the function a second time must not change
//         anything, because the COM is already the origin after the first call.
// -----------------------------------------------------------------------------
void test_soc_is_idempotent()
{
    std::cerr << "\n[TEST] test_soc_is_idempotent\n"
              << "  Source file:   src/reactions/subtract_off_com_position.cpp\n"
              << "  Function:      subtract_off_com_position(Molecule&)\n"
              << "  Scenario:      call the function twice in a row.\n"
              << "  Pass criteria: the second call leaves the molecule exactly\n"
              << "                 as the first call left it.\n";

    Molecule mol = soc_make_molecule(Coord { 8.0, -2.0, 6.0 },
        { Coord { 9.0, -2.0, 6.0 }, Coord { 8.0, 1.0, 6.0 } });

    std::cerr << "  First call...\n";
    subtract_off_com_position(mol);
    soc_dump_molecule(mol, "after first call ");

    // Snapshot the state produced by the first call.
    const std::vector<Coord> snapshot { mol.interfaceList[0].coord, mol.interfaceList[1].coord };

    std::cerr << "  Second call...\n";
    subtract_off_com_position(mol);
    soc_dump_molecule(mol, "after second call");

    // Nothing should have moved on the second application.
    for (std::size_t i = 0; i < snapshot.size() && i < mol.interfaceList.size(); ++i) {
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.x, snapshot[i].x)
            << "iface[" << i << "].x should be stable across repeated calls";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.y, snapshot[i].y)
            << "iface[" << i << "].y should be stable across repeated calls";
        EXPECT_DOUBLE_EQ(mol.interfaceList[i].coord.z, snapshot[i].z)
            << "iface[" << i << "].z should be stable across repeated calls";
    }
    EXPECT_DOUBLE_EQ(mol.comCoord.x, 0.0) << "COM x should still be 0";
    EXPECT_DOUBLE_EQ(mol.comCoord.y, 0.0) << "COM y should still be 0";
    EXPECT_DOUBLE_EQ(mol.comCoord.z, 0.0) << "COM z should still be 0";
}

// -----------------------------------------------------------------------------
// Test 7: Non-geometric molecule fields must not be disturbed. The function
//         should only touch comCoord and interfaceList coordinates.
// -----------------------------------------------------------------------------
void test_soc_leaves_other_fields_untouched()
{
    std::cerr << "\n[TEST] test_soc_leaves_other_fields_untouched\n"
              << "  Source file:   src/reactions/subtract_off_com_position.cpp\n"
              << "  Function:      subtract_off_com_position(Molecule&)\n"
              << "  Scenario:      molecule carries bookkeeping data (indices,\n"
              << "                 bound flags, temporary coords).\n"
              << "  Pass criteria: only comCoord and iface coords change; all\n"
              << "                 other fields keep their original values.\n";

    Molecule mol = soc_make_molecule(Coord { 5.0, 5.0, 5.0 },
        { Coord { 6.0, 5.0, 5.0 } });

    // Populate bookkeeping fields we expect the function to ignore.
    mol.index = 42;
    mol.myComIndex = 7;
    mol.molTypeIndex = 3;
    mol.mass = 1.5;
    mol.isLipid = true;
    mol.interfaceList[0].isBound = true;
    mol.interfaceList[0].relIndex = 0;
    mol.interfaceList[0].molTypeIndex = 3;
    mol.tmpComCoord = Coord { 99.0, 99.0, 99.0 };

    std::cerr << "  Calling subtract_off_com_position...\n";
    subtract_off_com_position(mol);
    soc_dump_molecule(mol, "after ");

    // Bookkeeping data must be preserved verbatim.
    EXPECT_EQ(mol.index, 42) << "Molecule index must be preserved";
    EXPECT_EQ(mol.myComIndex, 7) << "Parent complex index must be preserved";
    EXPECT_EQ(mol.molTypeIndex, 3) << "MolTemplate index must be preserved";
    EXPECT_DOUBLE_EQ(mol.mass, 1.5) << "Mass must be preserved";
    EXPECT_TRUE(mol.isLipid) << "isLipid flag must be preserved";
    ASSERT_FALSE(mol.interfaceList.empty()) << "Interface should still exist";
    EXPECT_TRUE(mol.interfaceList[0].isBound) << "Interface isBound flag must be preserved";
    EXPECT_EQ(mol.interfaceList[0].relIndex, 0) << "Interface relIndex must be preserved";
    EXPECT_EQ(mol.interfaceList[0].molTypeIndex, 3)
        << "Interface molTypeIndex must be preserved";

    // Temporary association coordinates are not part of this operation.
    EXPECT_DOUBLE_EQ(mol.tmpComCoord.x, 99.0) << "tmpComCoord.x must be untouched";
    EXPECT_DOUBLE_EQ(mol.tmpComCoord.y, 99.0) << "tmpComCoord.y must be untouched";
    EXPECT_DOUBLE_EQ(mol.tmpComCoord.z, 99.0) << "tmpComCoord.z must be untouched";

    // And the geometry itself must have been shifted as expected.
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.x, 1.0)
        << "Interface x should be 6 - 5 = 1";
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.y, 0.0)
        << "Interface y should be 5 - 5 = 0";
    EXPECT_DOUBLE_EQ(mol.interfaceList[0].coord.z, 0.0)
        << "Interface z should be 5 - 5 = 0";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* helper runs inside its own TEST so the
// framework reports individual pass/fail results while still running them all.
// -----------------------------------------------------------------------------
TEST(SubtractOffComPosition, ShiftsInterfacesAndZeroesCom) { test_soc_shifts_interfaces_and_zeroes_com(); }
TEST(SubtractOffComPosition, ComAlreadyAtOriginIsNoop) { test_soc_com_already_at_origin_is_noop(); }
TEST(SubtractOffComPosition, EmptyInterfaceList) { test_soc_empty_interface_list(); }
TEST(SubtractOffComPosition, InterfaceCoincidentWithCom) { test_soc_interface_coincident_with_com(); }
TEST(SubtractOffComPosition, PreservesRelativeGeometry) { test_soc_preserves_relative_geometry(); }
TEST(SubtractOffComPosition, IsIdempotent) { test_soc_is_idempotent(); }
TEST(SubtractOffComPosition, LeavesOtherFieldsUntouched) { test_soc_leaves_other_fields_untouched(); }