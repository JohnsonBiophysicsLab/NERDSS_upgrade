/*! \file test_write_xyz_assoc_cout.cpp
 *
 * ### Unit test for ../src/io/write_xyz_assoc_cout.cpp
 *
 * The function under test is:
 *
 *     void write_xyz_assoc_cout(const Complex& reactCom1,
 *                               const Complex& reactCom2,
 *                               const std::vector<Molecule>& moleculeList)
 *
 * It is a debugging helper that dumps the *temporary* association coordinates
 * (Molecule::tmpComCoord and Molecule::tmpICoords) of the two Complexes taking
 * part in an association event to standard output.
 *
 * Since the routine has no return value and mutates nothing, the only way to
 * test it is to capture std::cout and inspect the produced text.  The tests
 * below therefore:
 *
 *   1. Redirect std::cout into an in-memory buffer (RAII helper).
 *   2. Call write_xyz_assoc_cout() with hand-built Complexes / Molecules.
 *   3. Assert on the structure of the emitted text: the "Complex : <index>"
 *      headers, one "Mol type : <molTypeIndex>" header per member molecule,
 *      the expected number of coordinate lines, the ordering of the two
 *      complexes, and the presence of the actual coordinate values.
 *
 * All assertions are non-fatal (EXPECT_*) so that every test runs to
 * completion even when one of them fails.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <ios>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// RAII helper that redirects std::cout into a string buffer for the lifetime of
// the object, and restores both the original stream buffer *and* the original
// formatting flags on destruction.
//
// Restoring the flags matters here: write_xyz_assoc_cout() applies std::fixed
// to std::cout and never undoes it, so without this cleanup the global stream
// state would leak into other tests in the suite.
// -----------------------------------------------------------------------------
class WxacCoutCapture {
public:
    WxacCoutCapture()
        : buffer_()
        , oldBuf_(std::cout.rdbuf(buffer_.rdbuf()))
        , oldFlags_(std::cout.flags())
        , oldPrecision_(std::cout.precision())
    {
    }

    ~WxacCoutCapture()
    {
        // Put std::cout back exactly the way we found it.
        std::cout.rdbuf(oldBuf_);
        std::cout.flags(oldFlags_);
        std::cout.precision(oldPrecision_);
    }

    //! \brief Everything written to std::cout since construction.
    std::string str() const { return buffer_.str(); }

private:
    std::ostringstream buffer_;
    std::streambuf* oldBuf_;
    std::ios_base::fmtflags oldFlags_;
    std::streamsize oldPrecision_;
};

/*! \brief Build a minimal Molecule carrying temporary association coordinates.
 *
 * write_xyz_assoc_cout() only ever reads molTypeIndex, tmpComCoord and
 * tmpICoords, so nothing else has to be populated.
 *
 * \param[in] molTypeIndex value printed in the "Mol type :" line
 * \param[in] tmpCom       temporary center-of-mass coordinate
 * \param[in] tmpIfaces    temporary interface coordinates
 */
Molecule wxac_make_molecule(int molTypeIndex, const Coord& tmpCom, const std::vector<Coord>& tmpIfaces)
{
    Molecule mol;
    mol.molTypeIndex = molTypeIndex;
    mol.tmpComCoord = tmpCom;
    mol.tmpICoords = tmpIfaces;
    return mol;
}

/*! \brief Build a Complex with a given index and member molecule list. */
Complex wxac_make_complex(int index, const std::vector<int>& memberList)
{
    Complex com;
    com.index = index;
    com.memberList = memberList;
    return com;
}

/*! \brief Count the number of newline-terminated lines in a captured dump. */
size_t wxac_count_lines(const std::string& text)
{
    size_t lines = 0;
    for (char c : text) {
        if (c == '\n')
            ++lines;
    }
    return lines;
}

/*! \brief Count how many times \p needle occurs inside \p haystack. */
size_t wxac_count_occurrences(const std::string& haystack, const std::string& needle)
{
    if (needle.empty())
        return 0;

    size_t count = 0;
    size_t pos = haystack.find(needle);
    while (pos != std::string::npos) {
        ++count;
        pos = haystack.find(needle, pos + needle.size());
    }
    return count;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Typical case - two complexes, each with one molecule that owns a
//         center of mass plus two interfaces.
//
// Pass criteria:
//   * both "Complex : <index>" headers are present, in argument order;
//   * exactly one "Mol type :" header per member molecule;
//   * the total number of printed lines matches
//         2 headers + per molecule (1 type line + 1 COM line + N iface lines).
// -----------------------------------------------------------------------------
void test_wxac_two_complexes_basic_structure()
{
    std::cerr << "\n[TEST] test_wxac_two_complexes_basic_structure\n"
              << "  Source file: src/io/write_xyz_assoc_cout.cpp\n"
              << "  Function:    write_xyz_assoc_cout()\n"
              << "  Scenario:    two complexes, one molecule each, two interfaces each.\n"
              << "  Criteria:    correct headers, header counts and total line count.\n";

    // moleculeList index 0 -> member of complex 7, index 1 -> member of complex 9.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wxac_make_molecule(
        /*molTypeIndex=*/3, Coord(12.5, -3.25, 0.75),
        { Coord(13.5, -3.25, 0.75), Coord(11.5, -3.25, 0.75) }));
    moleculeList.push_back(wxac_make_molecule(
        /*molTypeIndex=*/5, Coord(-20.125, 4.5, -6.25),
        { Coord(-21.125, 4.5, -6.25), Coord(-19.125, 4.5, -6.25) }));

    Complex com1 = wxac_make_complex(7, { 0 });
    Complex com2 = wxac_make_complex(9, { 1 });

    std::string dump;
    {
        // Capture std::cout only for the duration of the call.
        WxacCoutCapture capture;
        write_xyz_assoc_cout(com1, com2, moleculeList);
        dump = capture.str();
    }

    std::cerr << "  Captured output:\n" << dump;

    // --- headers -------------------------------------------------------------
    const size_t posCom1 = dump.find("Complex : 7");
    const size_t posCom2 = dump.find("Complex : 9");
    EXPECT_NE(posCom1, std::string::npos) << "Header for the first complex (index 7) is missing";
    EXPECT_NE(posCom2, std::string::npos) << "Header for the second complex (index 9) is missing";
    if (posCom1 != std::string::npos && posCom2 != std::string::npos) {
        EXPECT_LT(posCom1, posCom2) << "reactCom1 must be printed before reactCom2";
    }
    EXPECT_EQ(wxac_count_occurrences(dump, "Complex :"), 2u)
        << "Exactly two complex headers should be printed";

    // --- molecule type lines -------------------------------------------------
    EXPECT_EQ(wxac_count_occurrences(dump, "Mol type :"), 2u)
        << "One 'Mol type :' line per member molecule (two total) expected";
    EXPECT_NE(dump.find("Mol type : 3"), std::string::npos)
        << "molTypeIndex 3 of the first molecule should be printed";
    EXPECT_NE(dump.find("Mol type : 5"), std::string::npos)
        << "molTypeIndex 5 of the second molecule should be printed";

    // --- total line count ----------------------------------------------------
    // 2 complex headers + 2 * (1 type line + 1 COM line + 2 interface lines) = 10
    const size_t expectedLines = 2 + 2 * (1 + 1 + 2);
    EXPECT_EQ(wxac_count_lines(dump), expectedLines)
        << "Unexpected number of output lines for 2 complexes x 1 molecule x 2 interfaces";
}

// -----------------------------------------------------------------------------
// Test 2: The *temporary* coordinates must be the ones printed, not the
//         permanent comCoord / interfaceList coordinates.
//
// Pass criteria: the distinctive tmp values appear in the dump while the
//                distinctive permanent values do not.
// -----------------------------------------------------------------------------
void test_wxac_prints_temporary_coordinates()
{
    std::cerr << "\n[TEST] test_wxac_prints_temporary_coordinates\n"
              << "  Source file: src/io/write_xyz_assoc_cout.cpp\n"
              << "  Function:    write_xyz_assoc_cout()\n"
              << "  Scenario:    molecule has different permanent and temporary coords.\n"
              << "  Criteria:    only tmpComCoord / tmpICoords values are emitted.\n";

    std::vector<Molecule> moleculeList;

    // Temporary coordinates use easily recognisable values (12.5 / 34.5 / 56.5),
    // whereas the permanent coordinates use completely different ones
    // (777.5 / 888.5 / 999.5) so we can tell which set was printed.
    Molecule mol = wxac_make_molecule(1, Coord(12.5, 34.5, 56.5), { Coord(13.5, 34.5, 56.5) });
    mol.comCoord = Coord(777.5, 888.5, 999.5);
    Molecule::Iface permanentIface;
    permanentIface.coord = Coord(777.5, 888.5, 999.5);
    mol.interfaceList.push_back(permanentIface);
    moleculeList.push_back(mol);

    Complex com1 = wxac_make_complex(0, { 0 });
    Complex com2 = wxac_make_complex(1, {}); // second complex intentionally empty

    std::string dump;
    {
        WxacCoutCapture capture;
        write_xyz_assoc_cout(com1, com2, moleculeList);
        dump = capture.str();
    }

    std::cerr << "  Captured output:\n" << dump;

    // The temporary COM and interface values must show up.
    EXPECT_NE(dump.find("12.5"), std::string::npos)
        << "tmpComCoord.x (12.5) should appear in the output";
    EXPECT_NE(dump.find("34.5"), std::string::npos)
        << "tmpComCoord.y (34.5) should appear in the output";
    EXPECT_NE(dump.find("56.5"), std::string::npos)
        << "tmpComCoord.z (56.5) should appear in the output";
    EXPECT_NE(dump.find("13.5"), std::string::npos)
        << "tmpICoords[0].x (13.5) should appear in the output";

    // The permanent coordinates must NOT show up.
    EXPECT_EQ(dump.find("777.5"), std::string::npos)
        << "Permanent comCoord.x (777.5) must not be printed by this function";
    EXPECT_EQ(dump.find("888.5"), std::string::npos)
        << "Permanent comCoord.y (888.5) must not be printed by this function";
    EXPECT_EQ(dump.find("999.5"), std::string::npos)
        << "Permanent comCoord.z (999.5) must not be printed by this function";
}

// -----------------------------------------------------------------------------
// Test 3: Degenerate inputs - both complexes have empty member lists.
//
// Pass criteria: only the two complex headers are printed (two lines), and no
//                "Mol type :" line is emitted.
// -----------------------------------------------------------------------------
void test_wxac_empty_member_lists()
{
    std::cerr << "\n[TEST] test_wxac_empty_member_lists\n"
              << "  Source file: src/io/write_xyz_assoc_cout.cpp\n"
              << "  Function:    write_xyz_assoc_cout()\n"
              << "  Scenario:    both complexes contain no member molecules.\n"
              << "  Criteria:    only the two 'Complex :' headers are emitted.\n";

    std::vector<Molecule> moleculeList; // deliberately empty

    Complex com1 = wxac_make_complex(11, {});
    Complex com2 = wxac_make_complex(22, {});

    std::string dump;
    {
        WxacCoutCapture capture;
        write_xyz_assoc_cout(com1, com2, moleculeList);
        dump = capture.str();
    }

    std::cerr << "  Captured output:\n" << dump;

    EXPECT_EQ(wxac_count_lines(dump), 2u) << "Only the two complex header lines should be printed";
    EXPECT_NE(dump.find("Complex : 11"), std::string::npos) << "First complex header missing";
    EXPECT_NE(dump.find("Complex : 22"), std::string::npos) << "Second complex header missing";
    EXPECT_EQ(wxac_count_occurrences(dump, "Mol type :"), 0u)
        << "No molecule lines should be printed when the member lists are empty";
}

// -----------------------------------------------------------------------------
// Test 4: A molecule with no interfaces (e.g. a point molecule whose temporary
//         interface list has not been populated) must still print its type and
//         its temporary center-of-mass line, and nothing more.
// -----------------------------------------------------------------------------
void test_wxac_molecule_without_interfaces()
{
    std::cerr << "\n[TEST] test_wxac_molecule_without_interfaces\n"
              << "  Source file: src/io/write_xyz_assoc_cout.cpp\n"
              << "  Function:    write_xyz_assoc_cout()\n"
              << "  Scenario:    single molecule with an empty tmpICoords vector.\n"
              << "  Criteria:    3 lines total (2 headers + type + COM), no iface lines.\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wxac_make_molecule(/*molTypeIndex=*/42, Coord(1.5, 2.5, 3.5), {}));

    Complex com1 = wxac_make_complex(0, { 0 });
    Complex com2 = wxac_make_complex(1, {});

    std::string dump;
    {
        WxacCoutCapture capture;
        write_xyz_assoc_cout(com1, com2, moleculeList);
        dump = capture.str();
    }

    std::cerr << "  Captured output:\n" << dump;

    // 2 complex headers + 1 "Mol type" line + 1 COM line = 4 lines.
    EXPECT_EQ(wxac_count_lines(dump), 4u)
        << "Expected 2 headers + 1 type line + 1 COM line for an interface-less molecule";
    EXPECT_NE(dump.find("Mol type : 42"), std::string::npos)
        << "The molecule's molTypeIndex (42) should be printed";
    EXPECT_NE(dump.find("1.5"), std::string::npos)
        << "The temporary center-of-mass x value (1.5) should be printed";
}

// -----------------------------------------------------------------------------
// Test 5: Multiple molecules per complex, including the case where the two
//         complexes share the same molecule index (which happens in the code
//         base when the same Molecule object is inspected twice).
//
// Pass criteria: one "Mol type :" block per entry in each memberList (i.e. the
//                shared molecule is printed twice), and the aggregate line
//                count matches the sum over both member lists.
// -----------------------------------------------------------------------------
void test_wxac_multiple_and_shared_members()
{
    std::cerr << "\n[TEST] test_wxac_multiple_and_shared_members\n"
              << "  Source file: src/io/write_xyz_assoc_cout.cpp\n"
              << "  Function:    write_xyz_assoc_cout()\n"
              << "  Scenario:    complex 1 has two members, complex 2 re-uses member 0.\n"
              << "  Criteria:    3 molecule blocks printed, line count matches.\n";

    std::vector<Molecule> moleculeList;
    // index 0: one interface
    moleculeList.push_back(wxac_make_molecule(0, Coord(0.5, 0.5, 0.5), { Coord(1.5, 0.5, 0.5) }));
    // index 1: three interfaces
    moleculeList.push_back(wxac_make_molecule(
        1, Coord(-0.5, -0.5, -0.5),
        { Coord(-1.5, -0.5, -0.5), Coord(-0.5, -1.5, -0.5), Coord(-0.5, -0.5, -1.5) }));

    Complex com1 = wxac_make_complex(100, { 0, 1 });
    Complex com2 = wxac_make_complex(200, { 0 }); // deliberately re-uses molecule 0

    std::string dump;
    {
        WxacCoutCapture capture;
        write_xyz_assoc_cout(com1, com2, moleculeList);
        dump = capture.str();
    }

    std::cerr << "  Captured output:\n" << dump;

    // Three molecule blocks in total (two from com1, one from com2).
    EXPECT_EQ(wxac_count_occurrences(dump, "Mol type :"), 3u)
        << "Three molecule blocks expected (2 in reactCom1, 1 in reactCom2)";
    // "Mol type : 0" appears twice because molecule 0 is a member of both.
    EXPECT_EQ(wxac_count_occurrences(dump, "Mol type : 0"), 2u)
        << "The shared molecule (molTypeIndex 0) should be printed once per member list entry";
    EXPECT_EQ(wxac_count_occurrences(dump, "Mol type : 1"), 1u)
        << "Molecule with molTypeIndex 1 should be printed exactly once";

    // Line budget: 2 headers
    //            + molecule 0 twice   -> 2 * (1 type + 1 COM + 1 iface) = 6
    //            + molecule 1 once    -> 1 * (1 type + 1 COM + 3 ifaces) = 5
    const size_t expectedLines = 2 + 6 + 5;
    EXPECT_EQ(wxac_count_lines(dump), expectedLines)
        << "Unexpected total line count for the multi-member / shared-member case";

    // Both complex headers must be present and in order.
    const size_t posCom1 = dump.find("Complex : 100");
    const size_t posCom2 = dump.find("Complex : 200");
    EXPECT_NE(posCom1, std::string::npos) << "Header 'Complex : 100' missing";
    EXPECT_NE(posCom2, std::string::npos) << "Header 'Complex : 200' missing";
    if (posCom1 != std::string::npos && posCom2 != std::string::npos) {
        EXPECT_LT(posCom1, posCom2) << "reactCom1 block must precede the reactCom2 block";
    }
}

// -----------------------------------------------------------------------------
// Test 6: Documented side effect - the function applies std::fixed to the
//         output stream and never restores it.  We verify the behaviour on a
//         local stream-state snapshot so the rest of the suite is unaffected
//         (the capture helper restores the flags for us).
// -----------------------------------------------------------------------------
void test_wxac_applies_fixed_formatting()
{
    std::cerr << "\n[TEST] test_wxac_applies_fixed_formatting\n"
              << "  Source file: src/io/write_xyz_assoc_cout.cpp\n"
              << "  Function:    write_xyz_assoc_cout()\n"
              << "  Scenario:    inspect std::cout formatting flags after the call.\n"
              << "  Criteria:    std::ios::fixed is set (documented side effect).\n";

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(wxac_make_molecule(0, Coord(2.5, 2.5, 2.5), { Coord(3.5, 2.5, 2.5) }));

    Complex com1 = wxac_make_complex(0, { 0 });
    Complex com2 = wxac_make_complex(1, {});

    std::ios_base::fmtflags flagsAfterCall {};
    {
        WxacCoutCapture capture; // restores flags on scope exit
        std::cout.unsetf(std::ios::fixed); // start from a known state
        write_xyz_assoc_cout(com1, com2, moleculeList);
        flagsAfterCall = std::cout.flags();
    }

    // The implementation streams with std::fixed, which is sticky.
    EXPECT_TRUE((flagsAfterCall & std::ios::fixed) != 0)
        << "write_xyz_assoc_cout() is expected to leave std::cout in fixed mode";

    std::cerr << "  std::ios::fixed set after call: "
              << (((flagsAfterCall & std::ios::fixed) != 0) ? "yes" : "no") << "\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each helper is invoked from its own TEST so a failure in
// one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(WriteXyzAssocCout, TwoComplexesBasicStructure) { test_wxac_two_complexes_basic_structure(); }
TEST(WriteXyzAssocCout, PrintsTemporaryCoordinates) { test_wxac_prints_temporary_coordinates(); }
TEST(WriteXyzAssocCout, EmptyMemberLists) { test_wxac_empty_member_lists(); }
TEST(WriteXyzAssocCout, MoleculeWithoutInterfaces) { test_wxac_molecule_without_interfaces(); }
TEST(WriteXyzAssocCout, MultipleAndSharedMembers) { test_wxac_multiple_and_shared_members(); }
TEST(WriteXyzAssocCout, AppliesFixedFormatting) { test_wxac_applies_fixed_formatting(); }