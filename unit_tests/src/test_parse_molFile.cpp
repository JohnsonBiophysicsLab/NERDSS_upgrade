/*! \file test_parse_molFile.cpp
 *
 * ### Unit test for src/parser/parse_molFile.cpp
 *
 * Function under test:
 *
 *     MolTemplate parse_molFile(std::string& mol);
 *
 * `parse_molFile()` opens the file `<mol>.mol` in the current working
 * directory and converts the (mostly free-form) molecule information file
 * into a fully populated `MolTemplate`.  It:
 *
 *   - skips whole-line comments (`#...`) and strips trailing comments,
 *   - strips *all* whitespace from every line before parsing,
 *   - dispatches `keyword = value` lines to `MolTemplate::set_value()`,
 *   - detects the internal-coordinate block (a line whose alphabetic prefix
 *     is "com" immediately followed by a digit) and hands it to
 *     `read_internal_coordinates()`,
 *   - delegates `state = ...` to `parse_states()` and `bonds = N` to
 *     `read_bonds()`,
 *   - assigns relative interface indices and *global* (static) state indices
 *     (`Interface::State::totalNumOfStates` / `MolTemplate::absToRelIface`),
 *   - computes the template `radius` as the largest COM-to-interface
 *     distance, and finally
 *   - copies that radius into `mass` when no mass was supplied (mass < 0).
 *
 * Every test below writes its own temporary `*.mol` file, parses it, checks
 * the resulting `MolTemplate`, and then deletes the file.  Note that
 * `parse_molFile()` calls `exit(1)` when the file is missing or a keyword is
 * unknown, so all of the fixtures written here are deliberately valid input
 * (the failure paths are intentionally NOT exercised: they would tear down
 * the whole gtest binary rather than fail a single case).
 */

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "parser/parser_functions.hpp"

namespace {

// -----------------------------------------------------------------------------
// Small helpers shared by every test in this file.  The `pmf_` prefix keeps
// these symbols from colliding with anything else in the unit-test binary.
// -----------------------------------------------------------------------------

/*! \brief Write `contents` to `<basename>.mol` in the current directory.
 *
 * \return true when the file could be created (the tests skip their
 *         assertions when this fails rather than crashing).
 */
bool pmf_write_mol_file(const std::string& basename, const std::string& contents)
{
    std::ofstream out(basename + ".mol");
    if (!out) {
        std::cerr << "  !! could not create fixture file " << basename << ".mol\n";
        return false;
    }
    out << contents;
    out.close();
    std::cerr << "  -> wrote fixture file " << basename << ".mol ("
              << contents.size() << " bytes)\n";
    return true;
}

/*! \brief Delete a fixture file created by pmf_write_mol_file(). */
void pmf_remove_mol_file(const std::string& basename)
{
    std::remove((basename + ".mol").c_str());
    std::cerr << "  -> removed fixture file " << basename << ".mol\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a "typical" molecule file.
//
// Checks name/copies/booleans/diffusion constants, the internal coordinate
// block, the automatically assigned interface + state indices, the derived
// radius (largest |COM - interface| vector) and the derived mass.
// -----------------------------------------------------------------------------
void pmf_test_basic_parsing()
{
    std::cerr << "\n[TEST] pmf_test_basic_parsing\n"
              << "  Source file: src/parser/parse_molFile.cpp\n"
              << "  Function:    parse_molFile()\n"
              << "  Scenario:    name/copies/checkOverlap/D/Dr keywords plus a\n"
              << "               COM + two interface coordinate block.\n"
              << "  Pass:        every field is transferred to the MolTemplate,\n"
              << "               radius == max(|COM-iface|) == 4.0 and, because\n"
              << "               no mass was given, mass == radius.\n";

    const std::string base { "unit_test_pmf_basic" };

    // Note: interface coordinates are (3,0,0) and (0,4,0) so the expected
    // radius is exactly 4.0 (sqrt(16) is exact in IEEE double arithmetic).
    const std::string contents =
        "# unit-test molecule file: basic keyword coverage\n"
        "Name = testmola\n"
        "copies = 25\n"
        "checkOverlap = true\n"
        "D = [10.0,20.0,30.0]\n"
        "Dr = [0.1,0.2,0.3]\n"
        "\n"
        "COM\t0.0000\t0.0000\t0.0000\n"
        "a\t3.0000\t0.0000\t0.0000\n"
        "b\t0.0000\t4.0000\t0.0000\n"
        "\n";

    if (!pmf_write_mol_file(base, contents)) {
        ADD_FAILURE() << "Could not create the .mol fixture; skipping assertions.";
        return;
    }

    // Snapshot the *global* counters that parse_molFile() mutates so the test
    // can assert on the delta instead of on an absolute value (other tests in
    // the suite share these statics).
    const int statesBefore { Interface::State::totalNumOfStates };
    const size_t absToRelBefore { MolTemplate::absToRelIface.size() };
    std::cerr << "  -> Interface::State::totalNumOfStates before parse = "
              << statesBefore << '\n';

    std::string molName { base }; // parse_molFile takes a non-const reference
    MolTemplate tmpl = parse_molFile(molName);

    // --- scalar keywords ------------------------------------------------
    std::cerr << "  Checking scalar keywords (name, copies, checkOverlap)...\n";
    EXPECT_EQ(tmpl.molName, "testmola") << "'Name = testmola' should be stored verbatim";
    EXPECT_EQ(tmpl.copies, 25) << "'copies = 25' should be parsed with stoi";
    EXPECT_TRUE(tmpl.checkOverlap) << "'checkOverlap = true' should set checkOverlap";

    // --- diffusion constants --------------------------------------------
    std::cerr << "  Checking diffusion constants D and Dr...\n";
    EXPECT_DOUBLE_EQ(tmpl.D.x, 10.0) << "D.x should be 10.0";
    EXPECT_DOUBLE_EQ(tmpl.D.y, 20.0) << "D.y should be 20.0";
    EXPECT_DOUBLE_EQ(tmpl.D.z, 30.0) << "D.z should be 30.0";
    EXPECT_DOUBLE_EQ(tmpl.Dr.x, 0.1) << "Dr.x should be 0.1";
    EXPECT_DOUBLE_EQ(tmpl.Dr.y, 0.2) << "Dr.y should be 0.2";
    EXPECT_DOUBLE_EQ(tmpl.Dr.z, 0.3) << "Dr.z should be 0.3";

    // --- internal coordinates -------------------------------------------
    std::cerr << "  Checking the internal coordinate block...\n";
    EXPECT_DOUBLE_EQ(tmpl.comCoord.x, 0.0) << "COM.x should be 0.0";
    EXPECT_DOUBLE_EQ(tmpl.comCoord.y, 0.0) << "COM.y should be 0.0";
    EXPECT_DOUBLE_EQ(tmpl.comCoord.z, 0.0) << "COM.z should be 0.0";

    ASSERT_EQ(tmpl.interfaceList.size(), 2u)
        << "Two interface lines were provided, so two interfaces must exist";

    EXPECT_EQ(tmpl.interfaceList[0].name, "a") << "First interface should be named 'a'";
    EXPECT_EQ(tmpl.interfaceList[1].name, "b") << "Second interface should be named 'b'";
    EXPECT_DOUBLE_EQ(tmpl.interfaceList[0].iCoord.x, 3.0) << "iface a x-coordinate";
    EXPECT_DOUBLE_EQ(tmpl.interfaceList[0].iCoord.y, 0.0) << "iface a y-coordinate";
    EXPECT_DOUBLE_EQ(tmpl.interfaceList[0].iCoord.z, 0.0) << "iface a z-coordinate";
    EXPECT_DOUBLE_EQ(tmpl.interfaceList[1].iCoord.x, 0.0) << "iface b x-coordinate";
    EXPECT_DOUBLE_EQ(tmpl.interfaceList[1].iCoord.y, 4.0) << "iface b y-coordinate";
    EXPECT_DOUBLE_EQ(tmpl.interfaceList[1].iCoord.z, 0.0) << "iface b z-coordinate";

    // --- relative interface indices --------------------------------------
    std::cerr << "  Checking relative interface indices (should be 0,1)...\n";
    EXPECT_EQ(tmpl.interfaceList[0].index, 0) << "First interface gets relative index 0";
    EXPECT_EQ(tmpl.interfaceList[1].index, 1) << "Second interface gets relative index 1";

    // --- automatically created (state-less) states -----------------------
    std::cerr << "  Checking the implicit single state created per interface...\n";
    ASSERT_EQ(tmpl.interfaceList[0].stateList.size(), 1u)
        << "An interface without explicit states gets exactly one default state";
    ASSERT_EQ(tmpl.interfaceList[1].stateList.size(), 1u)
        << "An interface without explicit states gets exactly one default state";
    EXPECT_EQ(tmpl.interfaceList[0].stateList[0].ifaceAndStateName, "a")
        << "The default state name is just the interface name";
    EXPECT_EQ(tmpl.interfaceList[1].stateList[0].ifaceAndStateName, "b")
        << "The default state name is just the interface name";
    EXPECT_EQ(tmpl.interfaceList[0].stateList[0].index, statesBefore)
        << "First default state takes the next free global state index";
    EXPECT_EQ(tmpl.interfaceList[1].stateList[0].index, statesBefore + 1)
        << "Second default state takes the following global state index";
    EXPECT_TRUE(tmpl.ifacesWithStates.empty())
        << "No interface declared explicit states, so ifacesWithStates stays empty";

    // --- global bookkeeping ----------------------------------------------
    std::cerr << "  Checking global static bookkeeping (totalNumOfStates, absToRelIface)...\n";
    EXPECT_EQ(Interface::State::totalNumOfStates, statesBefore + 2)
        << "Parsing two state-less interfaces must add two global states";
    ASSERT_EQ(MolTemplate::absToRelIface.size(), absToRelBefore + 2)
        << "absToRelIface must gain one entry per created state";
    EXPECT_EQ(MolTemplate::absToRelIface[absToRelBefore], 0)
        << "First new absolute index maps to relative interface 0";
    EXPECT_EQ(MolTemplate::absToRelIface[absToRelBefore + 1], 1)
        << "Second new absolute index maps to relative interface 1";

    // --- derived radius / mass -------------------------------------------
    std::cerr << "  Checking derived radius and auto-assigned mass...\n";
    std::cerr << "     radius = " << tmpl.radius << ", mass = " << tmpl.mass << '\n';
    EXPECT_DOUBLE_EQ(tmpl.radius, 4.0)
        << "radius must be the largest COM-to-interface distance (4.0)";
    EXPECT_DOUBLE_EQ(tmpl.mass, 4.0)
        << "mass defaults to -1 and is therefore replaced with the radius";

    // --- untouched defaults ----------------------------------------------
    std::cerr << "  Checking that unspecified flags keep their defaults...\n";
    EXPECT_FALSE(tmpl.isRod) << "isRod defaults to false";
    EXPECT_FALSE(tmpl.isLipid) << "isLipid defaults to false";
    EXPECT_FALSE(tmpl.isPoint) << "isPoint defaults to false";
    EXPECT_FALSE(tmpl.isImplicitLipid) << "isImplicitLipid defaults to false";
    EXPECT_FALSE(tmpl.countTransition) << "countTransition defaults to false";
    EXPECT_EQ(tmpl.transitionMatrixSize, 500) << "transitionMatrixSize defaults to 500";
    EXPECT_TRUE(tmpl.bondList.empty()) << "No bonds block was supplied";

    pmf_remove_mol_file(base);
}

// -----------------------------------------------------------------------------
// Test 2: an explicit mass must NOT be overwritten by the radius, and
//         `isImplicitLipid = true` must also switch `isLipid` on.
// -----------------------------------------------------------------------------
void pmf_test_explicit_mass_and_implicit_lipid()
{
    std::cerr << "\n[TEST] pmf_test_explicit_mass_and_implicit_lipid\n"
              << "  Source file: src/parser/parse_molFile.cpp\n"
              << "  Function:    parse_molFile()\n"
              << "  Scenario:    mass is given explicitly and isImplicitLipid=true.\n"
              << "  Pass:        mass keeps the user value (42.5), radius is still\n"
              << "               derived from the coordinates (1.5), and setting\n"
              << "               isImplicitLipid implies isLipid.\n";

    const std::string base { "unit_test_pmf_mass" };

    const std::string contents =
        "Name = testmolb\n"
        "copies = 3\n"
        "mass = 42.5\n"
        "isImplicitLipid = true\n"
        "D = [0.0,0.0,0.0]\n"
        "Dr = [0.0,0.0,0.0]\n"
        "\n"
        "COM\t0.0000\t0.0000\t0.0000\n"
        "head\t0.0000\t0.0000\t1.5000\n"
        "\n";

    if (!pmf_write_mol_file(base, contents)) {
        ADD_FAILURE() << "Could not create the .mol fixture; skipping assertions.";
        return;
    }

    std::string molName { base };
    MolTemplate tmpl = parse_molFile(molName);

    std::cerr << "  Checking mass/radius interaction...\n";
    std::cerr << "     radius = " << tmpl.radius << ", mass = " << tmpl.mass << '\n';
    EXPECT_DOUBLE_EQ(tmpl.mass, 42.5)
        << "An explicit mass (>= 0) must survive the auto-assignment step";
    EXPECT_DOUBLE_EQ(tmpl.radius, 1.5)
        << "radius is still computed from the single interface at z = 1.5";

    std::cerr << "  Checking implicit-lipid flag propagation...\n";
    EXPECT_TRUE(tmpl.isImplicitLipid) << "'isImplicitLipid = true' should be stored";
    EXPECT_TRUE(tmpl.isLipid)
        << "MolTemplate::set_value() forces isLipid=true for implicit lipids";

    std::cerr << "  Checking remaining fields...\n";
    EXPECT_EQ(tmpl.molName, "testmolb") << "Name should be testmolb";
    EXPECT_EQ(tmpl.copies, 3) << "copies should be 3";
    ASSERT_EQ(tmpl.interfaceList.size(), 1u) << "Exactly one interface was declared";
    EXPECT_EQ(tmpl.interfaceList[0].name, "head") << "Interface should be named 'head'";
    EXPECT_DOUBLE_EQ(tmpl.D.x, 0.0) << "D.x should be 0.0";
    EXPECT_DOUBLE_EQ(tmpl.Dr.z, 0.0) << "Dr.z should be 0.0";

    pmf_remove_mol_file(base);
}

// -----------------------------------------------------------------------------
// Test 3: the remaining boolean/integer keywords, including the compartment
//         keyword that also flips `crossesCompartment`.
// -----------------------------------------------------------------------------
void pmf_test_flag_keywords()
{
    std::cerr << "\n[TEST] pmf_test_flag_keywords\n"
              << "  Source file: src/parser/parse_molFile.cpp\n"
              << "  Function:    parse_molFile()\n"
              << "  Scenario:    isRod/isPoint/isPromoter/countTransition/\n"
              << "               transitionMatrixSize/outsideCompartment keywords.\n"
              << "  Pass:        each keyword lands in the matching MolTemplate\n"
              << "               member; outsideCompartment also sets\n"
              << "               crossesCompartment, insideCompartment stays false.\n";

    const std::string base { "unit_test_pmf_flags" };

    const std::string contents =
        "Name = testmolc\n"
        "isRod = true\n"
        "isPoint = true\n"
        "isPromoter = true\n"
        "countTransition = true\n"
        "transitionMatrixSize = 250\n"
        "outsideCompartment = true\n"
        "\n"
        "COM\t0.0000\t0.0000\t0.0000\n"
        "c1\t1.0000\t0.0000\t0.0000\n"
        "\n";

    if (!pmf_write_mol_file(base, contents)) {
        ADD_FAILURE() << "Could not create the .mol fixture; skipping assertions.";
        return;
    }

    std::string molName { base };
    MolTemplate tmpl = parse_molFile(molName);

    std::cerr << "  Checking boolean keywords...\n";
    EXPECT_TRUE(tmpl.isRod) << "'isRod = true' should set isRod";
    EXPECT_TRUE(tmpl.isPoint) << "'isPoint = true' should set isPoint";
    EXPECT_TRUE(tmpl.isPromoter) << "'isPromoter = true' should set isPromoter";
    EXPECT_TRUE(tmpl.countTransition) << "'countTransition = true' should set countTransition";

    std::cerr << "  Checking integer keyword transitionMatrixSize...\n";
    EXPECT_EQ(tmpl.transitionMatrixSize, 250)
        << "'transitionMatrixSize = 250' should override the default of 500";

    std::cerr << "  Checking compartment keywords...\n";
    EXPECT_TRUE(tmpl.outsideCompartment) << "'outsideCompartment = true' should be stored";
    EXPECT_TRUE(tmpl.crossesCompartment)
        << "outsideCompartment implies crossesCompartment in set_value()";
    EXPECT_FALSE(tmpl.insideCompartment)
        << "insideCompartment was never declared and must stay false";

    std::cerr << "  Checking that geometry was still parsed...\n";
    ASSERT_EQ(tmpl.interfaceList.size(), 1u) << "One interface line was supplied";
    EXPECT_DOUBLE_EQ(tmpl.radius, 1.0) << "radius = |(0,0,0) - (1,0,0)| = 1.0";
    EXPECT_DOUBLE_EQ(tmpl.mass, 1.0) << "mass falls back to the radius";

    pmf_remove_mol_file(base);
}

// -----------------------------------------------------------------------------
// Test 4: comment handling (whole-line and trailing) plus the default values
//         that survive when almost nothing is declared.
// -----------------------------------------------------------------------------
void pmf_test_comments_and_defaults()
{
    std::cerr << "\n[TEST] pmf_test_comments_and_defaults\n"
              << "  Source file: src/parser/parse_molFile.cpp\n"
              << "  Function:    parse_molFile()\n"
              << "  Scenario:    '#' comment lines and trailing '# ...' comments\n"
              << "               on keyword lines; only Name/copies declared.\n"
              << "  Pass:        comments are stripped (Name == 'testmold',\n"
              << "               copies == 7) and every other member keeps the\n"
              << "               MolTemplate default.\n";

    const std::string base { "unit_test_pmf_comments" };

    const std::string contents =
        "##\n"
        "# This whole line is a comment and must be ignored entirely.\n"
        "##\n"
        "Name = testmold # trailing comment after the molecule name\n"
        "copies = 7 # seven copies of this molecule\n"
        "\n"
        "COM\t0.0000\t0.0000\t0.0000\n"
        "d1\t0.0000\t0.0000\t1.0000\n"
        "\n";

    if (!pmf_write_mol_file(base, contents)) {
        ADD_FAILURE() << "Could not create the .mol fixture; skipping assertions.";
        return;
    }

    std::string molName { base };
    MolTemplate tmpl = parse_molFile(molName);

    std::cerr << "  Checking that trailing comments were removed...\n";
    EXPECT_EQ(tmpl.molName, "testmold")
        << "remove_comment() must strip everything from '#' onwards";
    EXPECT_EQ(tmpl.copies, 7)
        << "'copies = 7 # ...' must still parse as the integer 7";

    std::cerr << "  Checking untouched defaults...\n";
    EXPECT_FALSE(tmpl.checkOverlap) << "checkOverlap defaults to false";
    EXPECT_FALSE(tmpl.isLipid) << "isLipid defaults to false";
    EXPECT_FALSE(tmpl.isRod) << "isRod defaults to false";
    EXPECT_FALSE(tmpl.isImplicitLipid) << "isImplicitLipid defaults to false";
    EXPECT_FALSE(tmpl.countTransition) << "countTransition defaults to false";
    EXPECT_EQ(tmpl.transitionMatrixSize, 500) << "transitionMatrixSize defaults to 500";
    EXPECT_EQ(tmpl.molTypeIndex, 0) << "molTypeIndex is not assigned by parse_molFile";
    EXPECT_DOUBLE_EQ(tmpl.D.x, 0.0) << "D defaults to the zero Coord";
    EXPECT_DOUBLE_EQ(tmpl.D.y, 0.0) << "D defaults to the zero Coord";
    EXPECT_DOUBLE_EQ(tmpl.D.z, 0.0) << "D defaults to the zero Coord";
    EXPECT_DOUBLE_EQ(tmpl.Dr.x, 0.0) << "Dr defaults to the zero Coord";

    std::cerr << "  Checking geometry-derived values...\n";
    ASSERT_EQ(tmpl.interfaceList.size(), 1u) << "One interface was declared";
    EXPECT_EQ(tmpl.interfaceList[0].name, "d1") << "Interface should be named 'd1'";
    EXPECT_DOUBLE_EQ(tmpl.radius, 1.0) << "radius = 1.0 from the single interface";
    EXPECT_DOUBLE_EQ(tmpl.mass, 1.0) << "mass falls back to the radius";

    pmf_remove_mol_file(base);
}

// -----------------------------------------------------------------------------
// Test 5: the `bonds = N` branch, which forwards the next N lines to
//         read_bonds().  The fixture uses the repository's conventional bond
//         syntax ("com <ifaceName>") so read_bonds() can resolve every name.
//         Only the number of parsed bonds is asserted; the exact integer
//         encoding of "com" is an implementation detail of read_bonds().
// -----------------------------------------------------------------------------
void pmf_test_bonds_block()
{
    std::cerr << "\n[TEST] pmf_test_bonds_block\n"
              << "  Source file: src/parser/parse_molFile.cpp\n"
              << "  Function:    parse_molFile() -> read_bonds()\n"
              << "  Scenario:    a coordinate block followed by 'bonds = 2' and\n"
              << "               two 'com <iface>' bond definitions.\n"
              << "  Pass:        the coordinate block still yields two interfaces\n"
              << "               and bondList holds exactly two bonds.\n";

    const std::string base { "unit_test_pmf_bonds" };

    const std::string contents =
        "Name = testmole\n"
        "copies = 2\n"
        "\n"
        "COM\t0.0000\t0.0000\t0.0000\n"
        "a\t1.0000\t0.0000\t0.0000\n"
        "b\t0.0000\t2.0000\t0.0000\n"
        "\n"
        "bonds = 2\n"
        "com a\n"
        "com b\n";

    if (!pmf_write_mol_file(base, contents)) {
        ADD_FAILURE() << "Could not create the .mol fixture; skipping assertions.";
        return;
    }

    std::string molName { base };
    MolTemplate tmpl = parse_molFile(molName);

    std::cerr << "  Checking that the coordinate block ended at the blank line...\n";
    ASSERT_EQ(tmpl.interfaceList.size(), 2u)
        << "read_internal_coordinates() must stop before the 'bonds' keyword";
    EXPECT_EQ(tmpl.interfaceList[0].name, "a") << "First interface should be 'a'";
    EXPECT_EQ(tmpl.interfaceList[1].name, "b") << "Second interface should be 'b'";

    std::cerr << "  Checking the parsed bond list...\n";
    EXPECT_EQ(tmpl.bondList.size(), 2u)
        << "'bonds = 2' followed by two bond lines should yield two bonds";
    for (size_t bondItr = 0; bondItr < tmpl.bondList.size(); ++bondItr) {
        std::cerr << "     bond " << bondItr << " = [" << tmpl.bondList[bondItr][0]
                  << ", " << tmpl.bondList[bondItr][1] << "]\n";
    }

    std::cerr << "  Checking keywords that precede the coordinate block...\n";
    EXPECT_EQ(tmpl.molName, "testmole") << "Name should be testmole";
    EXPECT_EQ(tmpl.copies, 2) << "copies should be 2";
    EXPECT_DOUBLE_EQ(tmpl.radius, 2.0)
        << "radius is the larger of the two COM-iface distances (1.0 and 2.0)";

    pmf_remove_mol_file(base);
}

// -----------------------------------------------------------------------------
// Test 6: the global state counter keeps increasing across calls, i.e. parsing
//         the same file twice does not restart the absolute state indices.
// -----------------------------------------------------------------------------
void pmf_test_state_indices_accumulate()
{
    std::cerr << "\n[TEST] pmf_test_state_indices_accumulate\n"
              << "  Source file: src/parser/parse_molFile.cpp\n"
              << "  Function:    parse_molFile()\n"
              << "  Scenario:    parse the same molecule file twice.\n"
              << "  Pass:        Interface::State::totalNumOfStates advances by\n"
              << "               one per interface on *each* call, so the second\n"
              << "               template's state indices are strictly larger.\n";

    const std::string base { "unit_test_pmf_counter" };

    const std::string contents =
        "Name = testmolf\n"
        "copies = 1\n"
        "\n"
        "COM\t0.0000\t0.0000\t0.0000\n"
        "s1\t1.0000\t0.0000\t0.0000\n"
        "s2\t0.0000\t1.0000\t0.0000\n"
        "\n";

    if (!pmf_write_mol_file(base, contents)) {
        ADD_FAILURE() << "Could not create the .mol fixture; skipping assertions.";
        return;
    }

    const int statesBefore { Interface::State::totalNumOfStates };
    std::cerr << "  -> totalNumOfStates before any parse = " << statesBefore << '\n';

    std::string molName { base };
    MolTemplate first = parse_molFile(molName);
    const int statesAfterFirst { Interface::State::totalNumOfStates };
    std::cerr << "  -> totalNumOfStates after first parse = " << statesAfterFirst << '\n';

    MolTemplate second = parse_molFile(molName);
    const int statesAfterSecond { Interface::State::totalNumOfStates };
    std::cerr << "  -> totalNumOfStates after second parse = " << statesAfterSecond << '\n';

    std::cerr << "  Checking counter growth (two interfaces per parse)...\n";
    EXPECT_EQ(statesAfterFirst, statesBefore + 2)
        << "First parse creates one default state per interface";
    EXPECT_EQ(statesAfterSecond, statesAfterFirst + 2)
        << "Second parse continues from the previous counter value";

    ASSERT_EQ(first.interfaceList.size(), 2u) << "First template must have two interfaces";
    ASSERT_EQ(second.interfaceList.size(), 2u) << "Second template must have two interfaces";
    ASSERT_EQ(first.interfaceList[0].stateList.size(), 1u) << "One default state expected";
    ASSERT_EQ(second.interfaceList[0].stateList.size(), 1u) << "One default state expected";

    std::cerr << "  Checking that the second parse used fresh absolute indices...\n";
    EXPECT_EQ(first.interfaceList[0].stateList[0].index, statesBefore)
        << "First template's first state index";
    EXPECT_EQ(first.interfaceList[1].stateList[0].index, statesBefore + 1)
        << "First template's second state index";
    EXPECT_EQ(second.interfaceList[0].stateList[0].index, statesAfterFirst)
        << "Second template's first state index continues the global counter";
    EXPECT_EQ(second.interfaceList[1].stateList[0].index, statesAfterFirst + 1)
        << "Second template's second state index continues the global counter";

    std::cerr << "  Checking that both parses produced identical geometry...\n";
    EXPECT_DOUBLE_EQ(first.radius, second.radius) << "Both parses derive the same radius";
    EXPECT_DOUBLE_EQ(first.mass, second.mass) << "Both parses derive the same mass";
    EXPECT_EQ(first.molName, second.molName) << "Both parses read the same name";

    pmf_remove_mol_file(base);
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Every case is non-fatal (EXPECT_*) wherever possible so
// that a single bad expectation does not hide the remaining checks.
// -----------------------------------------------------------------------------
TEST(ParseMolFile, BasicParsing) { pmf_test_basic_parsing(); }
TEST(ParseMolFile, ExplicitMassAndImplicitLipid) { pmf_test_explicit_mass_and_implicit_lipid(); }
TEST(ParseMolFile, FlagKeywords) { pmf_test_flag_keywords(); }
TEST(ParseMolFile, CommentsAndDefaults) { pmf_test_comments_and_defaults(); }
TEST(ParseMolFile, BondsBlock) { pmf_test_bonds_block(); }
TEST(ParseMolFile, StateIndicesAccumulate) { pmf_test_state_indices_accumulate(); }