/*! \file test_print_system_information.cpp
 *
 * ### Unit test for src/io/print_system_information.cpp
 *
 * The single function under test is:
 *
 *     void print_system_information(long long int simItr,
 *                                  std::ofstream& systemFile,
 *                                  const std::vector<Molecule>& moleculeList,
 *                                  const std::vector<Complex>& complexList,
 *                                  const std::vector<MolTemplate>& molTemplateList)
 *
 * The function is a pure "dump the system state to a stream" routine; it has no
 * return value and mutates nothing.  The only observable behaviour therefore is
 * the text it writes.  Every test below
 *
 *   1. builds a small, completely synthetic system (molecules / complexes /
 *      molecule templates),
 *   2. calls print_system_information() with an std::ofstream pointed at a
 *      temporary file,
 *   3. reads the file back into a std::string, and
 *   4. asserts that the expected fields appear (or, for "empty" objects, do NOT
 *      appear) in the produced text.
 *
 * All assertions use EXPECT_* so that a single failing field does not abort the
 * remaining checks.  Verbose progress messages are written to stderr so a reader
 * of the test log can follow exactly which source file / function is exercised
 * and what pass criteria are applied.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Small helpers used to build the synthetic system.  All names are prefixed with
// "psi_" (print_system_information) so they cannot collide with helpers defined
// in other translation units of the test suite.
// -----------------------------------------------------------------------------

/*! \brief Build a MolTemplate with a given name and a single named interface.
 *
 * print_system_information() looks up
 *     molTemplateList[mol.molTypeIndex].molName
 * and
 *     molTemplateList[mol.molTypeIndex].interfaceList[iface.relIndex].name
 * so both fields must be populated for a non-empty molecule to be printable.
 *
 * \param[in] molName   Name reported for the molecule type.
 * \param[in] ifaceName Name of the (single) interface of this template.
 * \return A minimally populated MolTemplate.
 */
MolTemplate psi_make_template(const std::string& molName, const std::string& ifaceName)
{
    MolTemplate temp;
    temp.molName = molName;

    // One interface, relative index 0, sitting on the center of mass.
    Interface iface(ifaceName, Coord{ 0.0, 0.0, 0.0 });
    iface.index = 0;
    temp.interfaceList.clear();
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a non-empty Molecule with exactly one interface.
 *
 * \param[in] index        Value reported as "Index:".
 * \param[in] molTypeIndex Index into the molTemplateList (selects the name).
 * \param[in] comCoord     Center-of-mass coordinate.
 * \return A Molecule whose details will all be printed.
 */
Molecule psi_make_molecule(int index, int molTypeIndex, const Coord& comCoord)
{
    Molecule mol;
    mol.index = index;
    mol.isEmpty = false;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = 0;
    mol.mySubVolIndex = 5;
    mol.isLipid = false;
    mol.isPromoter = false;
    mol.comCoord = comCoord;

    // A single interface whose relative index (0) is valid for the template.
    Molecule::Iface iface;
    iface.relIndex = 0;
    iface.index = 0;
    iface.coord = comCoord;
    iface.stateIden = '\0'; // no state by default
    iface.isBound = false;  // not bound by default
    mol.interfaceList.clear();
    mol.interfaceList.push_back(iface);

    return mol;
}

/*! \brief Build a non-empty Complex owning the supplied member molecules.
 *
 * \param[in] index      Value reported as "Index:".
 * \param[in] memberList Indices of the member molecules.
 * \param[in] numEachMol Per-molecule-type counts (must be <= molTemplateList size).
 * \return A Complex whose details will all be printed.
 */
Complex psi_make_complex(int index, const std::vector<int>& memberList, const std::vector<int>& numEachMol)
{
    Complex com;
    com.index = index;
    com.isEmpty = false;
    com.mass = 12.5;
    com.radius = 3.75;
    com.comCoord = Coord{ 4.25, -5.25, 6.25 };
    com.D = Coord{ 1.25, 2.25, 3.25 };
    com.Dr = Coord{ 0.125, 0.25, 0.375 };
    com.memberList = memberList;
    com.numEachMol = numEachMol;
    return com;
}

/*! \brief Run print_system_information() into a temporary file and return the text.
 *
 * The ofstream is closed (and therefore flushed) before the file is read back,
 * and the temporary file is removed afterwards so repeated runs stay clean.
 *
 * \param[in] simItr          Iteration number handed to the function.
 * \param[in] moleculeList    Molecules to print.
 * \param[in] complexList     Complexes to print.
 * \param[in] molTemplateList Templates used for name lookups.
 * \param[in] tag             Unique tag used to build the temporary file name.
 * \return Complete contents of the file written by the function.
 */
std::string psi_run_and_read(long long int simItr, const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList,
    const std::string& tag)
{
    const std::string fileName = "test_print_system_information_" + tag + ".tmp";

    {
        std::ofstream systemFile(fileName);
        // Make sure the scratch file could actually be created; otherwise every
        // downstream assertion would fail for an unrelated reason.
        EXPECT_TRUE(systemFile.is_open()) << "Could not open temporary file " << fileName;
        if (systemFile.is_open()) {
            print_system_information(simItr, systemFile, moleculeList, complexList, molTemplateList);
        }
    } // ofstream destructor flushes and closes the file

    std::ifstream in(fileName);
    std::stringstream ss;
    ss << in.rdbuf();
    in.close();

    // Clean up the scratch file so the test leaves no artifacts behind.
    std::remove(fileName.c_str());

    return ss.str();
}

/*! \brief Convenience predicate: does \p haystack contain \p needle? */
bool psi_contains(const std::string& haystack, const std::string& needle)
{
    return haystack.find(needle) != std::string::npos;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: The iteration number and the two section headers must always be
//         written, even when the system contains nothing at all.
// -----------------------------------------------------------------------------
void test_psi_header_and_iteration()
{
    std::cerr << "\n[TEST] test_psi_header_and_iteration\n"
              << "  Source file:   src/io/print_system_information.cpp\n"
              << "  Function:      print_system_information()\n"
              << "  Scenario:      empty moleculeList and empty complexList.\n"
              << "  Pass criteria: output contains 'Iteration: 1234', a MOLECULES\n"
              << "                 header and a COMPLEXES header.\n";

    const std::vector<Molecule> moleculeList {}; // nothing to print
    const std::vector<Complex> complexList {};   // nothing to print
    const std::vector<MolTemplate> molTemplateList {};

    const std::string out = psi_run_and_read(1234, moleculeList, complexList, molTemplateList, "header");

    std::cerr << "  Wrote " << out.size() << " characters of output.\n";

    // Basic sanity: the function must have produced *something*.
    EXPECT_FALSE(out.empty()) << "print_system_information() should always emit the headers";

    // The iteration number is the very first thing written.
    EXPECT_TRUE(psi_contains(out, "Iteration: 1234"))
        << "Output should report the simulation iteration that was passed in";

    // Both section headers are unconditional.
    EXPECT_TRUE(psi_contains(out, "MOLECULES"))
        << "Output should contain the MOLECULES section header";
    EXPECT_TRUE(psi_contains(out, "COMPLEXES"))
        << "Output should contain the COMPLEXES section header";
}

// -----------------------------------------------------------------------------
// Test 2: A fully populated, non-empty molecule must have every scalar field
//         printed, and its interface block must be printed as well.
// -----------------------------------------------------------------------------
void test_psi_molecule_fields()
{
    std::cerr << "\n[TEST] test_psi_molecule_fields\n"
              << "  Source file:   src/io/print_system_information.cpp\n"
              << "  Function:      print_system_information() - molecule block\n"
              << "  Scenario:      one non-empty molecule of type 'clathrin' with\n"
              << "                 one unbound, stateless interface named 'iface1'.\n"
              << "  Pass criteria: index, empty flag, type name, parent complex,\n"
              << "                 sub volume, lipid/promoter flags, COM coordinate\n"
              << "                 components, interface indices and interface name\n"
              << "                 all appear in the output.\n";

    std::vector<MolTemplate> molTemplateList { psi_make_template("clathrin", "iface1") };

    // Distinctive coordinate values so we can grep for them individually.
    std::vector<Molecule> moleculeList { psi_make_molecule(7, 0, Coord{ 1.5, -2.5, 3.5 }) };
    const std::vector<Complex> complexList {};

    const std::string out = psi_run_and_read(0, moleculeList, complexList, molTemplateList, "molfields");

    // Scalar molecule fields.
    EXPECT_TRUE(psi_contains(out, "Index: 7")) << "Molecule index should be printed";
    EXPECT_TRUE(psi_contains(out, "Is empty: false"))
        << "A non-empty molecule should report 'Is empty: false' (std::boolalpha is set)";
    EXPECT_TRUE(psi_contains(out, "Type: clathrin"))
        << "Molecule type name should be looked up from the MolTemplate";
    EXPECT_TRUE(psi_contains(out, "Parent complex index: 0"))
        << "Parent complex index should be printed";
    EXPECT_TRUE(psi_contains(out, "Sub volume index: 5"))
        << "Sub volume index should be printed";
    EXPECT_TRUE(psi_contains(out, "Is a lipid: false"))
        << "Lipid flag should be printed with boolalpha";
    EXPECT_TRUE(psi_contains(out, "Is a promoter: false"))
        << "Promoter flag should be printed with boolalpha";

    // The center-of-mass coordinate is streamed via Coord's operator<<; we only
    // require that each component value shows up somewhere in the text.
    EXPECT_TRUE(psi_contains(out, "Center of mass coordinate:"))
        << "Center of mass label should be printed";
    EXPECT_TRUE(psi_contains(out, "1.5")) << "COM x component (1.5) should appear";
    EXPECT_TRUE(psi_contains(out, "-2.5")) << "COM y component (-2.5) should appear";
    EXPECT_TRUE(psi_contains(out, "3.5")) << "COM z component (3.5) should appear";

    // Interface block.
    EXPECT_TRUE(psi_contains(out, "Interfaces:")) << "Interface section should be printed";
    EXPECT_TRUE(psi_contains(out, "Relative index: 0"))
        << "Interface relative index should be printed";
    EXPECT_TRUE(psi_contains(out, "Absolute index: 0"))
        << "Interface absolute index should be printed";
    EXPECT_TRUE(psi_contains(out, "Interface name: iface1"))
        << "Interface name should come from the MolTemplate interfaceList";

    // With stateIden == '\0' and isBound == false these optional blocks must be
    // suppressed.
    EXPECT_FALSE(psi_contains(out, "Current state:"))
        << "No state should be printed when stateIden is the null character";
    EXPECT_FALSE(psi_contains(out, "Interaction:"))
        << "No interaction block should be printed for an unbound interface";
}

// -----------------------------------------------------------------------------
// Test 3: Optional interface sub-blocks (state identifier and interaction) must
//         be printed when the corresponding fields are set.
// -----------------------------------------------------------------------------
void test_psi_interface_state_and_interaction()
{
    std::cerr << "\n[TEST] test_psi_interface_state_and_interaction\n"
              << "  Source file:   src/io/print_system_information.cpp\n"
              << "  Function:      print_system_information() - interface block\n"
              << "  Scenario:      the interface has stateIden 'P' and is bound to\n"
              << "                 molecule 9, interface 1.\n"
              << "  Pass criteria: 'Current state: P', the interaction header and\n"
              << "                 both partner indices are printed.\n";

    std::vector<MolTemplate> molTemplateList { psi_make_template("kinase", "site") };

    std::vector<Molecule> moleculeList { psi_make_molecule(3, 0, Coord{ 0.0, 0.0, 0.0 }) };
    // Turn the interface into a phosphorylated, bound interface.
    moleculeList[0].interfaceList[0].stateIden = 'P';
    moleculeList[0].interfaceList[0].isBound = true;
    moleculeList[0].interfaceList[0].interaction.partnerIndex = 9;
    moleculeList[0].interfaceList[0].interaction.partnerIfaceIndex = 1;

    const std::vector<Complex> complexList {};

    const std::string out = psi_run_and_read(0, moleculeList, complexList, molTemplateList, "ifacestate");

    // State identifier block.
    EXPECT_TRUE(psi_contains(out, "Current state: P"))
        << "A non-null stateIden should be printed";

    // Interaction block (note: the source prints the partner interface index
    // without a trailing colon).
    EXPECT_TRUE(psi_contains(out, "Interaction:"))
        << "A bound interface should print the interaction header";
    EXPECT_TRUE(psi_contains(out, "Partner index: 9"))
        << "The bound partner's molecule index should be printed";
    EXPECT_TRUE(psi_contains(out, "Partner interface index 1"))
        << "The bound partner's interface index should be printed";
}

// -----------------------------------------------------------------------------
// Test 4: An empty (destroyed) molecule must report only its index and the empty
//         flag; none of the detail fields may be emitted.
// -----------------------------------------------------------------------------
void test_psi_empty_molecule_skipped()
{
    std::cerr << "\n[TEST] test_psi_empty_molecule_skipped\n"
              << "  Source file:   src/io/print_system_information.cpp\n"
              << "  Function:      print_system_information() - molecule block\n"
              << "  Scenario:      the only molecule has isEmpty == true.\n"
              << "  Pass criteria: 'Index: 11' and 'Is empty: true' are printed but\n"
              << "                 'Type:', 'Interfaces:' and the COM label are not.\n";

    std::vector<MolTemplate> molTemplateList { psi_make_template("ghost", "iface1") };

    // Build a molecule and then mark it as destroyed.
    Molecule emptyMol = psi_make_molecule(11, 0, Coord{ 8.5, 8.5, 8.5 });
    emptyMol.isEmpty = true;
    std::vector<Molecule> moleculeList { emptyMol };
    const std::vector<Complex> complexList {};

    const std::string out = psi_run_and_read(0, moleculeList, complexList, molTemplateList, "emptymol");

    // Always printed, regardless of isEmpty.
    EXPECT_TRUE(psi_contains(out, "Index: 11")) << "Empty molecules still report their index";
    EXPECT_TRUE(psi_contains(out, "Is empty: true"))
        << "An empty molecule should report 'Is empty: true'";

    // Everything guarded by "if (!mol.isEmpty)" must be absent.
    EXPECT_FALSE(psi_contains(out, "Type: ghost"))
        << "The type name must not be printed for an empty molecule";
    EXPECT_FALSE(psi_contains(out, "Interfaces:"))
        << "The interface list must not be printed for an empty molecule";
    EXPECT_FALSE(psi_contains(out, "Center of mass coordinate:"))
        << "The COM coordinate must not be printed for an empty molecule";
}

// -----------------------------------------------------------------------------
// Test 5: A fully populated, non-empty complex must print mass, radius, both
//         diffusion-constant triples, its member list and its per-type counts.
// -----------------------------------------------------------------------------
void test_psi_complex_fields()
{
    std::cerr << "\n[TEST] test_psi_complex_fields\n"
              << "  Source file:   src/io/print_system_information.cpp\n"
              << "  Function:      print_system_information() - complex block\n"
              << "  Scenario:      two molecule templates ('A','B'), one molecule\n"
              << "                 (so boolalpha is already set) and one non-empty\n"
              << "                 complex with members {0,1} and counts {2,1}.\n"
              << "  Pass criteria: mass, radius, Dx/Dy/Dz, Drx/Dry/Drz, the member\n"
              << "                 index list and the per-type counts are printed.\n";

    std::vector<MolTemplate> molTemplateList { psi_make_template("A", "a1"), psi_make_template("B", "b1") };

    // A single molecule is printed first; this both exercises the molecule loop
    // and leaves std::boolalpha enabled on the stream, which is what the complex
    // loop relies on when printing its own isEmpty flag.
    std::vector<Molecule> moleculeList { psi_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }) };

    std::vector<Complex> complexList { psi_make_complex(2, { 0, 1 }, { 2, 1 }) };

    const std::string out = psi_run_and_read(0, moleculeList, complexList, molTemplateList, "complexfields");

    // Complex identity / bulk properties.
    EXPECT_TRUE(psi_contains(out, "Index: 2")) << "Complex index should be printed";
    EXPECT_TRUE(psi_contains(out, "Mass: 12.5")) << "Complex mass should be printed";
    EXPECT_TRUE(psi_contains(out, "Radius: 3.75")) << "Complex radius should be printed";

    // Translational diffusion constants.
    EXPECT_TRUE(psi_contains(out, "Translational diffusion constants:"))
        << "Translational diffusion header should be printed";
    EXPECT_TRUE(psi_contains(out, "Dx = 1.25")) << "Dx should be printed";
    EXPECT_TRUE(psi_contains(out, "Dy = 2.25")) << "Dy should be printed";
    EXPECT_TRUE(psi_contains(out, "Dz = 3.25")) << "Dz should be printed";

    // Rotational diffusion constants.
    EXPECT_TRUE(psi_contains(out, "Rotational diffusion constants:"))
        << "Rotational diffusion header should be printed";
    EXPECT_TRUE(psi_contains(out, "Drx = 0.125")) << "Drx should be printed";
    EXPECT_TRUE(psi_contains(out, "Dry = 0.25")) << "Dry should be printed";
    EXPECT_TRUE(psi_contains(out, "Drz = 0.375")) << "Drz should be printed";

    // Member molecule indices are printed space-separated on one line.
    EXPECT_TRUE(psi_contains(out, "Member molecules: 0 1"))
        << "Member molecule indices should be printed space separated";

    // Per-molecule-type counts, labelled with the template names.
    EXPECT_TRUE(psi_contains(out, "Number of each molecule type:"))
        << "Per-type count header should be printed";
    EXPECT_TRUE(psi_contains(out, "A: 2")) << "Count for molecule type A should be 2";
    EXPECT_TRUE(psi_contains(out, "B: 1")) << "Count for molecule type B should be 1";
}

// -----------------------------------------------------------------------------
// Test 6: An empty (destroyed) complex must report only its index and the empty
//         flag; none of the detail fields may be emitted.
// -----------------------------------------------------------------------------
void test_psi_empty_complex_skipped()
{
    std::cerr << "\n[TEST] test_psi_empty_complex_skipped\n"
              << "  Source file:   src/io/print_system_information.cpp\n"
              << "  Function:      print_system_information() - complex block\n"
              << "  Scenario:      the only complex has isEmpty == true.\n"
              << "  Pass criteria: 'Index: 4' is printed but mass, radius, the\n"
              << "                 diffusion constants and the member list are not.\n";

    std::vector<MolTemplate> molTemplateList { psi_make_template("A", "a1") };
    const std::vector<Molecule> moleculeList {};

    Complex emptyCom = psi_make_complex(4, { 0 }, { 1 });
    emptyCom.isEmpty = true;
    std::vector<Complex> complexList { emptyCom };

    const std::string out = psi_run_and_read(0, moleculeList, complexList, molTemplateList, "emptycomplex");

    // Index is printed unconditionally.
    EXPECT_TRUE(psi_contains(out, "Index: 4")) << "Empty complexes still report their index";
    EXPECT_TRUE(psi_contains(out, "Is empty:"))
        << "The empty flag should be printed for every complex";

    // Everything guarded by "if (!complex.isEmpty)" must be absent.
    EXPECT_FALSE(psi_contains(out, "Mass:"))
        << "Mass must not be printed for an empty complex";
    EXPECT_FALSE(psi_contains(out, "Radius:"))
        << "Radius must not be printed for an empty complex";
    EXPECT_FALSE(psi_contains(out, "Translational diffusion constants:"))
        << "Diffusion constants must not be printed for an empty complex";
    EXPECT_FALSE(psi_contains(out, "Member molecules:"))
        << "Member list must not be printed for an empty complex";
}

// -----------------------------------------------------------------------------
// Test 7: Multiple molecules and multiple complexes are all visited; the loops
//         must not stop after the first element.
// -----------------------------------------------------------------------------
void test_psi_multiple_entries()
{
    std::cerr << "\n[TEST] test_psi_multiple_entries\n"
              << "  Source file:   src/io/print_system_information.cpp\n"
              << "  Function:      print_system_information() - both loops\n"
              << "  Scenario:      two molecules (types 'A' and 'B') and two\n"
              << "                 complexes (indices 0 and 1).\n"
              << "  Pass criteria: both molecule type names appear, and 'Index: 0'\n"
              << "                 and 'Index: 1' both appear (once per loop pass).\n";

    std::vector<MolTemplate> molTemplateList { psi_make_template("A", "a1"), psi_make_template("B", "b1") };

    // Two molecules of different types.
    std::vector<Molecule> moleculeList {
        psi_make_molecule(0, 0, Coord{ 1.0, 1.0, 1.0 }),
        psi_make_molecule(1, 1, Coord{ 2.0, 2.0, 2.0 }),
    };

    // Two complexes, each with a single member.
    std::vector<Complex> complexList {
        psi_make_complex(0, { 0 }, { 1, 0 }),
        psi_make_complex(1, { 1 }, { 0, 1 }),
    };

    const std::string out = psi_run_and_read(99, moleculeList, complexList, molTemplateList, "multi");

    // Both molecule types must have been resolved via molTemplateList.
    EXPECT_TRUE(psi_contains(out, "Type: A")) << "First molecule type should be printed";
    EXPECT_TRUE(psi_contains(out, "Type: B")) << "Second molecule type should be printed";

    // Both index values appear (twice each - once in the molecule loop and once
    // in the complex loop); we simply verify each string is present.
    EXPECT_TRUE(psi_contains(out, "Index: 0")) << "Index 0 should be printed";
    EXPECT_TRUE(psi_contains(out, "Index: 1")) << "Index 1 should be printed";

    // Sanity: both complex member lists should show up.
    EXPECT_TRUE(psi_contains(out, "Member molecules: 0"))
        << "First complex member list should be printed";
    EXPECT_TRUE(psi_contains(out, "Member molecules: 1"))
        << "Second complex member list should be printed";

    // And the iteration number is still correct.
    EXPECT_TRUE(psi_contains(out, "Iteration: 99")) << "Iteration number should be 99";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* helper is executed inside its own TEST
// so that failures are reported per scenario while every scenario still runs.
// -----------------------------------------------------------------------------
TEST(PrintSystemInformation, HeaderAndIteration) { test_psi_header_and_iteration(); }
TEST(PrintSystemInformation, MoleculeFields) { test_psi_molecule_fields(); }
TEST(PrintSystemInformation, InterfaceStateAndInteraction) { test_psi_interface_state_and_interaction(); }
TEST(PrintSystemInformation, EmptyMoleculeSkipped) { test_psi_empty_molecule_skipped(); }
TEST(PrintSystemInformation, ComplexFields) { test_psi_complex_fields(); }
TEST(PrintSystemInformation, EmptyComplexSkipped) { test_psi_empty_complex_skipped(); }
TEST(PrintSystemInformation, MultipleEntries) { test_psi_multiple_entries(); }