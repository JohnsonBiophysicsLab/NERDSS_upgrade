/*! \file test_parse_input.cpp
 *
 * ### Unit tests for src/parser/parse_input.cpp
 *
 * The translation unit under test provides:
 *
 *   - `Membrane::set_value_BC(std::string, BoundaryKeyword)`  -- assigns one
 *      boundary-condition keyword/value pair onto a Membrane object.
 *   - `Membrane::display()`                                   -- prints the
 *      boundary information to std::cout.
 *   - `Membrane::create_water_box()`                          -- builds a cubic
 *      water box that circumscribes the sphere of radius `sphereR`.
 *   - `create_tmp_line(const std::string&)`                   -- normalizes an
 *      input-file line (lower-case, whitespace stripped, comment removed).
 *   - `parse_input(...)` / `parse_input_for_add(...)`         -- top level
 *      input-file parsers.
 *
 * Notes on what is (and is not) exercised here:
 *   * The `start molecules` / `start reactions` blocks of `parse_input()` read
 *     `<molName>.mol` files from disk and call `exit()` on any malformed input,
 *     which would take the whole gtest binary down.  For that reason the
 *     end-to-end parser tests below only use the `start parameters` and
 *     `start boundaries` blocks, which are fully self contained.
 *   * `set_value_BC()` deliberately falls through to a `throw`/`exit(1)` for an
 *     unknown keyword, so that path is documented but never triggered.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "classes/class_Membrane.hpp"
#include "parser/parser_functions.hpp"

// -----------------------------------------------------------------------------
// Small local helpers (file-static so they cannot collide with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Write `contents` verbatim into `fileName` (overwriting it). */
void pi_write_text_file(const std::string& fileName, const std::string& contents)
{
    std::ofstream os(fileName.c_str());
    os << contents;
    os.close();
}

/*! \brief Run `Membrane::display()` while capturing everything it puts on cout. */
std::string pi_capture_display(Membrane& mem)
{
    std::ostringstream captured;
    std::streambuf* oldBuf = std::cout.rdbuf(captured.rdbuf());
    mem.display();
    std::cout.rdbuf(oldBuf); // always restore, even if display() printed nothing
    return captured.str();
}

} // namespace

// -----------------------------------------------------------------------------
// Test: Membrane::set_value_BC() -- water box + implicit lipid keywords.
//
// Criteria: WaterBox = [50,60,70] must populate x/y/z, derive volume = x*y*z,
//           set xLeft/xRight to -x/2 and +x/2, and flip isBox to true.
// -----------------------------------------------------------------------------
void test_pi_set_value_bc_waterbox()
{
    std::cerr << "\n[TEST] test_pi_set_value_bc_waterbox\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    Membrane::set_value_BC (BoundaryKeyword::waterBox,\n"
              << "               BoundaryKeyword::implicitLipid)\n";

    Membrane mem;

    // Sanity check on the defaults before we touch anything.
    EXPECT_FALSE(mem.isBox) << "A default Membrane must not be flagged as a box";
    EXPECT_FALSE(mem.implicitLipid) << "A default Membrane must not use implicit lipid";

    std::cerr << "  -> Feeding 'WaterBox = [50,60,70]' (already normalized) to set_value_BC\n";
    mem.set_value_BC("[50,60,70]", BoundaryKeyword::waterBox);

    EXPECT_DOUBLE_EQ(mem.waterBox.x, 50.0) << "waterBox.x should be 50 nm";
    EXPECT_DOUBLE_EQ(mem.waterBox.y, 60.0) << "waterBox.y should be 60 nm";
    EXPECT_DOUBLE_EQ(mem.waterBox.z, 70.0) << "waterBox.z should be 70 nm";
    // The WaterBox constructor derives these three quantities itself.
    EXPECT_DOUBLE_EQ(mem.waterBox.volume, 50.0 * 60.0 * 70.0)
        << "volume must be x*y*z as computed by the WaterBox constructor";
    EXPECT_DOUBLE_EQ(mem.waterBox.xLeft, -25.0) << "xLeft must be -x/2";
    EXPECT_DOUBLE_EQ(mem.waterBox.xRight, 25.0) << "xRight must be +x/2";
    EXPECT_TRUE(mem.isBox) << "Reading a waterBox must set isBox = true";

    // implicitLipid uses read_boolean(), so both "true" and "false" are valid.
    std::cerr << "  -> Feeding 'implicitLipid = true' then 'false'\n";
    mem.set_value_BC("true", BoundaryKeyword::implicitLipid);
    EXPECT_TRUE(mem.implicitLipid) << "implicitLipid should be true after reading 'true'";
    mem.set_value_BC("false", BoundaryKeyword::implicitLipid);
    EXPECT_FALSE(mem.implicitLipid) << "implicitLipid should be false after reading 'false'";
}

// -----------------------------------------------------------------------------
// Test: Membrane::set_value_BC() -- sphere and compartment keywords.
//
// Criteria: sphereR implicitly turns isSphere on; compartmentR implicitly turns
//           hasCompartment on; the explicit boolean keywords work on their own.
// -----------------------------------------------------------------------------
void test_pi_set_value_bc_sphere_and_compartment()
{
    std::cerr << "\n[TEST] test_pi_set_value_bc_sphere_and_compartment\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    Membrane::set_value_BC (isSphere, sphereR,\n"
              << "               hasCompartment, compartmentR)\n";

    // --- explicit booleans -----------------------------------------------
    Membrane memFlags;
    std::cerr << "  -> isSphere = true / hasCompartment = true\n";
    memFlags.set_value_BC("true", BoundaryKeyword::isSphere);
    memFlags.set_value_BC("true", BoundaryKeyword::hasCompartment);
    EXPECT_TRUE(memFlags.isSphere) << "isSphere should follow the parsed boolean";
    EXPECT_TRUE(memFlags.hasCompartment) << "hasCompartment should follow the parsed boolean";

    // --- numeric radii, which additionally set the corresponding flags ----
    Membrane memRadii;
    EXPECT_FALSE(memRadii.isSphere) << "Fresh Membrane should not be a sphere";
    EXPECT_FALSE(memRadii.hasCompartment) << "Fresh Membrane should have no compartment";

    std::cerr << "  -> sphereR = 100 should also set isSphere = true\n";
    memRadii.set_value_BC("100", BoundaryKeyword::sphereR);
    EXPECT_DOUBLE_EQ(memRadii.sphereR, 100.0) << "sphereR should be 100 nm";
    EXPECT_TRUE(memRadii.isSphere) << "Providing sphereR must imply isSphere = true";

    std::cerr << "  -> compartmentR = 30.5 should also set hasCompartment = true\n";
    memRadii.set_value_BC("30.5", BoundaryKeyword::compartmentR);
    EXPECT_DOUBLE_EQ(memRadii.compartmentR, 30.5) << "compartmentR should be 30.5 nm";
    EXPECT_TRUE(memRadii.hasCompartment) << "Providing compartmentR must imply hasCompartment = true";
}

// -----------------------------------------------------------------------------
// Test: Membrane::set_value_BC() -- boundary condition type strings and the
//       droplet (compartment site) properties.
//
// Criteria: the three BC-type keywords store the raw string; compartmentSiteD
//           and compartmentSiteRho are stored on the nested Droplet struct.
// -----------------------------------------------------------------------------
void test_pi_set_value_bc_types_and_droplet()
{
    std::cerr << "\n[TEST] test_pi_set_value_bc_types_and_droplet\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    Membrane::set_value_BC (xBCtype, yBCtype, zBCtype,\n"
              << "               compartmentSiteD, compartmentSiteRho)\n";

    Membrane mem;

    std::cerr << "  -> Storing the three boundary condition type strings\n";
    mem.set_value_BC("reflect", BoundaryKeyword::xBCtype);
    mem.set_value_BC("pbc", BoundaryKeyword::yBCtype);
    mem.set_value_BC("reflect", BoundaryKeyword::zBCtype);
    EXPECT_EQ(mem.xBCtype, std::string("reflect")) << "xBCtype should be stored verbatim";
    EXPECT_EQ(mem.yBCtype, std::string("pbc")) << "yBCtype should be stored verbatim";
    EXPECT_EQ(mem.zBCtype, std::string("reflect")) << "zBCtype should be stored verbatim";

    std::cerr << "  -> Storing droplet diffusion constant and site density\n";
    EXPECT_DOUBLE_EQ(mem.droplet.D, 0.0) << "droplet.D defaults to 0";
    EXPECT_DOUBLE_EQ(mem.droplet.rho, 0.0) << "droplet.rho defaults to 0";
    mem.set_value_BC("0.25", BoundaryKeyword::compartmentSiteD);
    mem.set_value_BC("0.001", BoundaryKeyword::compartmentSiteRho);
    EXPECT_DOUBLE_EQ(mem.droplet.D, 0.25) << "compartmentSiteD should land on droplet.D";
    EXPECT_DOUBLE_EQ(mem.droplet.rho, 0.001) << "compartmentSiteRho should land on droplet.rho";
}

// -----------------------------------------------------------------------------
// Test: Membrane::create_water_box()
//
// Criteria: the derived box is cubic with side 2*sphereR.  (The function does
//           NOT recompute volume/xLeft/xRight, so those must stay untouched.)
// -----------------------------------------------------------------------------
void test_pi_create_water_box()
{
    std::cerr << "\n[TEST] test_pi_create_water_box\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    Membrane::create_water_box\n";

    Membrane mem;
    mem.sphereR = 42.0;

    // Record the untouched fields so we can prove the function leaves them alone.
    const double volBefore = mem.waterBox.volume;
    const double leftBefore = mem.waterBox.xLeft;
    const double rightBefore = mem.waterBox.xRight;

    std::cerr << "  -> Creating a water box around a sphere of radius 42 nm\n";
    mem.create_water_box();

    EXPECT_DOUBLE_EQ(mem.waterBox.x, 84.0) << "waterBox.x should be 2*sphereR";
    EXPECT_DOUBLE_EQ(mem.waterBox.y, 84.0) << "waterBox.y should be 2*sphereR";
    EXPECT_DOUBLE_EQ(mem.waterBox.z, 84.0) << "waterBox.z should be 2*sphereR";

    EXPECT_DOUBLE_EQ(mem.waterBox.volume, volBefore)
        << "create_water_box() does not recompute volume";
    EXPECT_DOUBLE_EQ(mem.waterBox.xLeft, leftBefore)
        << "create_water_box() does not recompute xLeft";
    EXPECT_DOUBLE_EQ(mem.waterBox.xRight, rightBefore)
        << "create_water_box() does not recompute xRight";
}

// -----------------------------------------------------------------------------
// Test: Membrane::display()
//
// Criteria: the routine must always print the sphere/compartment lines, and the
//           box dimensions only when isBox == true.  We capture std::cout and
//           look for the identifying substrings.
// -----------------------------------------------------------------------------
void test_pi_display()
{
    std::cerr << "\n[TEST] test_pi_display\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    Membrane::display\n";

    // --- case 1: no box has been read in ---------------------------------
    Membrane noBox;
    noBox.sphereR = 12.0;
    noBox.isSphere = true;
    std::cerr << "  -> Displaying a sphere-only Membrane (isBox == false)\n";
    const std::string out1 = pi_capture_display(noBox);
    EXPECT_NE(out1.find("isSphere?"), std::string::npos) << "Output should report isSphere";
    EXPECT_NE(out1.find("sphere Radius"), std::string::npos) << "Output should report the sphere radius";
    EXPECT_NE(out1.find("hasCompartment?"), std::string::npos) << "Output should report hasCompartment";
    EXPECT_NE(out1.find("hasImplicitLipid?"), std::string::npos) << "Output should report implicit lipid usage";
    EXPECT_EQ(out1.find("BOX geometry"), std::string::npos)
        << "Box dimensions must NOT be printed when isBox is false";

    // --- case 2: a water box has been read in ----------------------------
    Membrane withBox;
    withBox.set_value_BC("[10,20,30]", BoundaryKeyword::waterBox); // also sets isBox
    std::cerr << "  -> Displaying a Membrane that owns a water box (isBox == true)\n";
    const std::string out2 = pi_capture_display(withBox);
    EXPECT_NE(out2.find("BOX geometry"), std::string::npos)
        << "Box dimensions must be printed when isBox is true";
    EXPECT_NE(out2.find("10"), std::string::npos) << "Output should contain the x dimension";
    EXPECT_NE(out2.find("20"), std::string::npos) << "Output should contain the y dimension";
    EXPECT_NE(out2.find("30"), std::string::npos) << "Output should contain the z dimension";
}

// -----------------------------------------------------------------------------
// Test: create_tmp_line()
//
// Criteria: the returned line is lower-cased, has every whitespace character
//           removed, and has any trailing '#' comment stripped.  Comment-only
//           lines must be rejected by skipLine(), which is how the parser uses
//           the pair of functions together.
// -----------------------------------------------------------------------------
void test_pi_create_tmp_line()
{
    std::cerr << "\n[TEST] test_pi_create_tmp_line\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    create_tmp_line\n";

    // Block delimiters are recognised only after normalization.
    std::cerr << "  -> '   Start   Parameters  ' should normalize to 'startparameters'\n";
    EXPECT_EQ(create_tmp_line("   Start   Parameters  "), std::string("startparameters"))
        << "Spaces removed and text lower-cased";

    std::cerr << "  -> '\\t end Boundaries \\t' should normalize to 'endboundaries'\n";
    EXPECT_EQ(create_tmp_line("\t end Boundaries \t"), std::string("endboundaries"))
        << "Tabs count as whitespace and must be removed";

    // Trailing comments must be stripped, keeping only the keyword/value pair.
    std::cerr << "  -> 'WaterBox = [10, 20, 30] # comment' should keep only the assignment\n";
    EXPECT_EQ(create_tmp_line("WaterBox = [10, 20, 30] # a trailing comment"),
        std::string("waterbox=[10,20,30]"))
        << "Comment removed, spaces removed, characters lower-cased";

    // An empty line stays empty and is skipped by the parser.
    std::cerr << "  -> The empty string stays empty and is skipped\n";
    EXPECT_TRUE(create_tmp_line("").empty()) << "An empty input line yields an empty result";
    EXPECT_TRUE(skipLine(create_tmp_line(""))) << "skipLine() must reject an empty line";

    // A pure comment line is always discarded by the parser (either because the
    // comment stripper emptied it, or because it still begins with '#').
    std::cerr << "  -> '   # only a comment' must be skipped by the parser\n";
    EXPECT_TRUE(skipLine(create_tmp_line("   # only a comment")))
        << "A comment-only line must be skipped by the caller";
}

// -----------------------------------------------------------------------------
// Test: parse_input() -- parameters + boundaries blocks.
//
// Criteria: a temporary .inp file containing only the (self contained)
//           parameters and boundaries blocks must end up populating both the
//           Parameters object and the Membrane object, and must leave every
//           reaction container untouched (empty).
// -----------------------------------------------------------------------------
void test_pi_parse_input_parameters_and_boundaries()
{
    std::cerr << "\n[TEST] test_pi_parse_input_parameters_and_boundaries\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    parse_input\n"
              << "  Scenario:    a temporary .inp file with only 'start parameters'\n"
              << "               and 'start boundaries' blocks (no .mol files needed).\n";

    std::string fileName { "unit_test_parse_input_basic.inp" };

    // Note: keyword matching is case-insensitive and whitespace-insensitive,
    // which is exactly what create_tmp_line()/parse_paramFile() take care of.
    const std::string contents =
        "# temporary input file generated by the unit test\n"
        "start parameters\n"
        "    nItr = 1234\n"
        "    timeStep = 0.5\n"
        "    timeWrite = 10\n"
        "    trajWrite = 20\n"
        "    restartWrite = 50\n"
        "    name = unittestsim\n"
        "end parameters\n"
        "\n"
        "start boundaries\n"
        "    WaterBox = [50, 60, 70]   # cubic-ish box, in nm\n"
        "    xBCtype = reflect\n"
        "    yBCtype = reflect\n"
        "    zBCtype = reflect\n"
        "    implicitLipid = false\n"
        "end boundaries\n";

    pi_write_text_file(fileName, contents);

    // Everything the parser needs, all empty / default constructed.
    Parameters params;
    std::map<std::string, int> observableList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<CreateDestructRxn> createDestructRxns;
    std::vector<TransmissionRxn> transmissionRxns;
    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject;

    std::cerr << "  Calling parse_input on \"" << fileName << "\"...\n";
    parse_input(fileName, params, observableList, forwardRxns, backRxns, createDestructRxns,
        transmissionRxns, molTemplateList, membraneObject);

    // --- parameters block -------------------------------------------------
    std::cerr << "  Checking values read from the parameters block\n";
    EXPECT_EQ(params.nItr, 1234LL) << "nItr should be 1234 timesteps";
    EXPECT_DOUBLE_EQ(params.timeStep, 0.5) << "timeStep should be 0.5 us";
    EXPECT_EQ(params.timeWrite, 10LL) << "timeWrite should be 10 timesteps";
    EXPECT_EQ(params.trajWrite, 20LL) << "trajWrite should be 20 timesteps";
    EXPECT_EQ(params.restartWrite, 50LL) << "restartWrite should be 50 timesteps";
    EXPECT_EQ(params.name, std::string("unittestsim")) << "name should be 'unittestsim'";

    // --- boundaries block -------------------------------------------------
    std::cerr << "  Checking values read from the boundaries block\n";
    EXPECT_TRUE(membraneObject.isBox) << "A waterBox keyword must set isBox";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.x, 50.0) << "waterBox.x should be 50 nm";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.y, 60.0) << "waterBox.y should be 60 nm";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.z, 70.0) << "waterBox.z should be 70 nm";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.volume, 50.0 * 60.0 * 70.0)
        << "waterBox volume should be x*y*z";
    EXPECT_EQ(membraneObject.xBCtype, std::string("reflect")) << "xBCtype should be 'reflect'";
    EXPECT_EQ(membraneObject.yBCtype, std::string("reflect")) << "yBCtype should be 'reflect'";
    EXPECT_EQ(membraneObject.zBCtype, std::string("reflect")) << "zBCtype should be 'reflect'";
    EXPECT_FALSE(membraneObject.implicitLipid) << "implicitLipid was declared false";
    EXPECT_FALSE(membraneObject.isSphere) << "No sphere keyword was provided";
    EXPECT_FALSE(membraneObject.hasCompartment) << "No compartment keyword was provided";

    // --- containers that must remain untouched ----------------------------
    std::cerr << "  Checking that no reactions / molecules were created\n";
    EXPECT_TRUE(forwardRxns.empty()) << "No reactions block was supplied";
    EXPECT_TRUE(backRxns.empty()) << "No reactions block was supplied";
    EXPECT_TRUE(createDestructRxns.empty()) << "No reactions block was supplied";
    EXPECT_TRUE(transmissionRxns.empty()) << "No reactions block was supplied";
    EXPECT_TRUE(molTemplateList.empty()) << "No molecules block was supplied";
    EXPECT_TRUE(observableList.empty()) << "No observables block was supplied";

    // isNonEQ is derived from the presence of creation/destruction reactions.
    EXPECT_FALSE(params.isNonEQ) << "Without create/destruct reactions isNonEQ must be false";
    // numTotalSpecies is copied straight out of the RxnBase static counter.
    EXPECT_EQ(params.numTotalSpecies, RxnBase::totRxnSpecies)
        << "numTotalSpecies is assigned from RxnBase::totRxnSpecies";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test: parse_input() with a spherical boundary.
//
// Criteria: sphereR read from the file implies isSphere, and create_water_box()
//           can then be used to derive the enclosing cube.
// -----------------------------------------------------------------------------
void test_pi_parse_input_sphere_boundary()
{
    std::cerr << "\n[TEST] test_pi_parse_input_sphere_boundary\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    parse_input (boundaries block, sphere geometry)\n";

    std::string fileName { "unit_test_parse_input_sphere.inp" };

    const std::string contents =
        "start parameters\n"
        "    nItr = 10\n"
        "    timeStep = 1.0\n"
        "end parameters\n"
        "start boundaries\n"
        "    sphereR = 250.0\n"
        "    hasCompartment = true\n"
        "    compartmentR = 100.0\n"
        "    compartmentSiteD = 0.5\n"
        "    compartmentSiteRho = 0.02\n"
        "end boundaries\n";

    pi_write_text_file(fileName, contents);

    Parameters params;
    std::map<std::string, int> observableList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<CreateDestructRxn> createDestructRxns;
    std::vector<TransmissionRxn> transmissionRxns;
    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject;

    std::cerr << "  Calling parse_input on \"" << fileName << "\"...\n";
    parse_input(fileName, params, observableList, forwardRxns, backRxns, createDestructRxns,
        transmissionRxns, molTemplateList, membraneObject);

    std::cerr << "  Checking the sphere / compartment values\n";
    EXPECT_DOUBLE_EQ(membraneObject.sphereR, 250.0) << "sphereR should be 250 nm";
    EXPECT_TRUE(membraneObject.isSphere) << "Reading sphereR implies isSphere = true";
    EXPECT_FALSE(membraneObject.isBox) << "No waterBox keyword was provided, so isBox stays false";
    EXPECT_TRUE(membraneObject.hasCompartment) << "hasCompartment was declared true";
    EXPECT_DOUBLE_EQ(membraneObject.compartmentR, 100.0) << "compartmentR should be 100 nm";
    EXPECT_DOUBLE_EQ(membraneObject.droplet.D, 0.5) << "compartmentSiteD should be 0.5";
    EXPECT_DOUBLE_EQ(membraneObject.droplet.rho, 0.02) << "compartmentSiteRho should be 0.02";

    // The water box for a spherical system is built afterwards from sphereR.
    std::cerr << "  Deriving the enclosing water box via create_water_box()\n";
    membraneObject.create_water_box();
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.x, 500.0) << "Box side should be 2*sphereR";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.y, 500.0) << "Box side should be 2*sphereR";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.z, 500.0) << "Box side should be 2*sphereR";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test: parse_input() ignores unknown boundary keywords instead of aborting.
//
// Criteria: a bogus keyword prints a warning and is skipped, while the valid
//           keywords around it are still applied.
// -----------------------------------------------------------------------------
void test_pi_parse_input_unknown_keyword_is_ignored()
{
    std::cerr << "\n[TEST] test_pi_parse_input_unknown_keyword_is_ignored\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    parse_input (boundaries block keyword lookup)\n"
              << "  Pass criteria: unknown keyword warns but does not abort, and\n"
              << "                 surrounding valid keywords are still honoured.\n";

    std::string fileName { "unit_test_parse_input_unknown.inp" };

    const std::string contents =
        "start boundaries\n"
        "    notARealKeyword = 17\n"
        "    WaterBox = [12, 12, 12]\n"
        "    xBCtype = pbc\n"
        "end boundaries\n";

    pi_write_text_file(fileName, contents);

    Parameters params;
    std::map<std::string, int> observableList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<CreateDestructRxn> createDestructRxns;
    std::vector<TransmissionRxn> transmissionRxns;
    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject;

    std::cerr << "  Calling parse_input on \"" << fileName << "\"...\n";
    parse_input(fileName, params, observableList, forwardRxns, backRxns, createDestructRxns,
        transmissionRxns, molTemplateList, membraneObject);

    EXPECT_TRUE(membraneObject.isBox) << "The valid waterBox line after the bad line must be read";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.x, 12.0) << "waterBox.x should be 12 nm";
    EXPECT_EQ(membraneObject.xBCtype, std::string("pbc")) << "xBCtype should be 'pbc'";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// Test: parse_input_for_add() -- parameters + boundaries blocks.
//
// Criteria: the "add" variant behaves like parse_input() for these two blocks,
//           overwriting the Parameters/Membrane values that a restart had set,
//           and leaving pre-existing reaction lists alone.
// -----------------------------------------------------------------------------
void test_pi_parse_input_for_add_parameters_and_boundaries()
{
    std::cerr << "\n[TEST] test_pi_parse_input_for_add_parameters_and_boundaries\n"
              << "  Source file: src/parser/parse_input.cpp\n"
              << "  Function:    parse_input_for_add\n"
              << "  Scenario:    an 'add' file containing only parameters and\n"
              << "               boundaries blocks (no .mol files required).\n";

    std::string fileName { "unit_test_parse_input_add.inp" };

    const std::string contents =
        "start parameters\n"
        "    nItr = 777\n"
        "    timeStep = 0.1\n"
        "    timeWrite = 5\n"
        "    restartWrite = 25\n"
        "end parameters\n"
        "start boundaries\n"
        "    WaterBox = [80, 80, 80]\n"
        "    implicitLipid = true\n"
        "end boundaries\n";

    pi_write_text_file(fileName, contents);

    // Pretend a restart already gave us some values; the add-file must override.
    Parameters params;
    params.nItr = 1;
    params.timeStep = 999.0;

    std::map<std::string, int> observableList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<CreateDestructRxn> createDestructRxns;
    std::vector<TransmissionRxn> transmissionRxns;
    std::vector<MolTemplate> molTemplateList;
    Membrane membraneObject;

    const std::size_t forwardBefore = forwardRxns.size();
    const std::size_t molTempBefore = molTemplateList.size();

    std::cerr << "  Calling parse_input_for_add on \"" << fileName << "\"...\n";
    parse_input_for_add(fileName, params, observableList, forwardRxns, backRxns, createDestructRxns,
        transmissionRxns, molTemplateList, membraneObject, /*numDoubleBeforeAdd=*/0);

    std::cerr << "  Checking that the add-file overwrote the parameters\n";
    EXPECT_EQ(params.nItr, 777LL) << "nItr should be replaced with the add-file value";
    EXPECT_DOUBLE_EQ(params.timeStep, 0.1) << "timeStep should be replaced with the add-file value";
    EXPECT_EQ(params.timeWrite, 5LL) << "timeWrite should be 5 timesteps";
    EXPECT_EQ(params.restartWrite, 25LL) << "restartWrite should be 25 timesteps";

    std::cerr << "  Checking the boundaries read by the add-file\n";
    EXPECT_TRUE(membraneObject.isBox) << "A waterBox keyword must set isBox";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.x, 80.0) << "waterBox.x should be 80 nm";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.y, 80.0) << "waterBox.y should be 80 nm";
    EXPECT_DOUBLE_EQ(membraneObject.waterBox.z, 80.0) << "waterBox.z should be 80 nm";
    EXPECT_TRUE(membraneObject.implicitLipid) << "implicitLipid was declared true";

    std::cerr << "  Checking that no new reactions or molecule templates appeared\n";
    EXPECT_EQ(forwardRxns.size(), forwardBefore) << "No reactions block was supplied";
    EXPECT_EQ(molTemplateList.size(), molTempBefore) << "No molecules block was supplied";
    EXPECT_FALSE(params.isNonEQ) << "Without create/destruct reactions isNonEQ must be false";

    std::remove(fileName.c_str());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- one per named test_* helper so a failure in one does
// not prevent the remaining checks from running.
// -----------------------------------------------------------------------------
TEST(ParseInput, MembraneSetValueBCWaterBox) { test_pi_set_value_bc_waterbox(); }
TEST(ParseInput, MembraneSetValueBCSphereAndCompartment) { test_pi_set_value_bc_sphere_and_compartment(); }
TEST(ParseInput, MembraneSetValueBCTypesAndDroplet) { test_pi_set_value_bc_types_and_droplet(); }
TEST(ParseInput, MembraneCreateWaterBox) { test_pi_create_water_box(); }
TEST(ParseInput, MembraneDisplay) { test_pi_display(); }
TEST(ParseInput, CreateTmpLine) { test_pi_create_tmp_line(); }
TEST(ParseInput, ParseInputParametersAndBoundaries) { test_pi_parse_input_parameters_and_boundaries(); }
TEST(ParseInput, ParseInputSphereBoundary) { test_pi_parse_input_sphere_boundary(); }
TEST(ParseInput, ParseInputUnknownKeywordIgnored) { test_pi_parse_input_unknown_keyword_is_ignored(); }
TEST(ParseInput, ParseInputForAddParametersAndBoundaries) { test_pi_parse_input_for_add_parameters_and_boundaries(); }