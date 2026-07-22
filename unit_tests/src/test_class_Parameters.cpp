/*! \file test_class_parameter.cpp
 * \brief Unit tests for the Parameters class member functions defined in
 *        class_parameter.cpp
 *
 * Tests cover:
 *   - Parameters::set_value()      : Sets member parameters from string values
 *                                    keyed on the ParamKeyword enumeration.
 *   - Parameters::parse_paramFile(): Parses a parameter block from an ifstream.
 *   - Parameters::display()        : Prints the current parameter values.
 *
 * NOTE: We use gtest and print verbose information to stderr describing which
 *       function is being tested and what each test is doing.
 */

#include "classes/class_Parameters.hpp"
#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <cstdio>

/*!
 * \brief Test the Parameters::set_value function for numeric integer keywords.
 *
 * We call set_value with several ParamKeyword enumerations and verify that the
 * corresponding member variable of the Parameters object is set correctly.
 */
void test_set_value_integers()
{
    std::cerr << "\n[TEST] test_set_value_integers()\n"
              << "  File under test : class_parameter.cpp\n"
              << "  Function        : Parameters::set_value()\n"
              << "  Purpose         : Verify integer/long parameters are parsed correctly.\n";

    Parameters params;

    // numMolTypes (key 0)
    std::cerr << "  -> Setting numMolTypes = \"3\"\n";
    params.set_value("3", ParamKeyword::numMolTypes);
    EXPECT_EQ(params.numMolTypes, 3) << "numMolTypes should be parsed to 3";

    // numTotalSpecies (key 1)
    std::cerr << "  -> Setting numTotalSpecies = \"7\"\n";
    params.set_value("7", ParamKeyword::numTotalSpecies);
    EXPECT_EQ(params.numTotalSpecies, 7) << "numTotalSpecies should be parsed to 7";

    // nItr (key 2) - long long int
    std::cerr << "  -> Setting nItr = \"1000000\"\n";
    params.set_value("1000000", ParamKeyword::nItr);
    EXPECT_EQ(params.nItr, 1000000LL) << "nItr should be parsed to 1000000";

    // numTotalComplex (key 7)
    std::cerr << "  -> Setting numTotalComplex = \"12\"\n";
    params.set_value("12", ParamKeyword::numTotalComplex);
    EXPECT_EQ(params.numTotalComplex, 12) << "numTotalComplex should be parsed to 12";
}

/*!
 * \brief Test the Parameters::set_value function for keywords with long long
 *        timestep-interval values.
 */
void test_set_value_intervals()
{
    std::cerr << "\n[TEST] test_set_value_intervals()\n"
              << "  File under test : class_parameter.cpp\n"
              << "  Function        : Parameters::set_value()\n"
              << "  Purpose         : Verify write-interval (long long) parameters parse correctly.\n";

    Parameters params;

    // timeWrite (key 4)
    std::cerr << "  -> Setting timeWrite = \"500\"\n";
    params.set_value("500", ParamKeyword::timeWrite);
    EXPECT_EQ(params.timeWrite, 500LL) << "timeWrite should be 500";

    // trajWrite (key 5)
    std::cerr << "  -> Setting trajWrite = \"250\"\n";
    params.set_value("250", ParamKeyword::trajWrite);
    EXPECT_EQ(params.trajWrite, 250LL) << "trajWrite should be 250";

    // restartWrite (key 10)
    std::cerr << "  -> Setting restartWrite = \"1000\"\n";
    params.set_value("1000", ParamKeyword::restartWrite);
    EXPECT_EQ(params.restartWrite, 1000LL) << "restartWrite should be 1000";

    // pdbWrite (key 11)
    std::cerr << "  -> Setting pdbWrite = \"100\"\n";
    params.set_value("100", ParamKeyword::pdbWrite);
    EXPECT_EQ(params.pdbWrite, 100LL) << "pdbWrite should be 100";

    // checkPoint (key 14)
    std::cerr << "  -> Setting checkPoint = \"2000\"\n";
    params.set_value("2000", ParamKeyword::checkPoint);
    EXPECT_EQ(params.checkPoint, 2000LL) << "checkPoint should be 2000";

    // transitionWrite (key 16)
    std::cerr << "  -> Setting transitionWrite = \"777\"\n";
    params.set_value("777", ParamKeyword::transitionWrite);
    EXPECT_EQ(params.transitionWrite, 777LL) << "transitionWrite should be 777";

    // bondedComplexWrite (key 20)
    std::cerr << "  -> Setting bondedComplexWrite = \"333\"\n";
    params.set_value("333", ParamKeyword::bondedComplexWrite);
    EXPECT_EQ(params.bondedComplexWrite, 333LL) << "bondedComplexWrite should be 333";
}

/*!
 * \brief Test the Parameters::set_value function for double-valued keywords.
 */
void test_set_value_doubles()
{
    std::cerr << "\n[TEST] test_set_value_doubles()\n"
              << "  File under test : class_parameter.cpp\n"
              << "  Function        : Parameters::set_value()\n"
              << "  Purpose         : Verify floating-point parameters parse correctly.\n";

    Parameters params;

    // timeStep (key 6)
    std::cerr << "  -> Setting timeStep = \"0.5\"\n";
    params.set_value("0.5", ParamKeyword::timeStep);
    EXPECT_DOUBLE_EQ(params.timeStep, 0.5) << "timeStep should be 0.5";

    // mass (key 8)
    std::cerr << "  -> Setting mass = \"2.5\"\n";
    params.set_value("2.5", ParamKeyword::mass);
    EXPECT_DOUBLE_EQ(params.mass, 2.5) << "mass should be 2.5";

    // overlapSepLimit (key 12)
    std::cerr << "  -> Setting overlapSepLimit = \"0.1\"\n";
    params.set_value("0.1", ParamKeyword::overlapSepLimit);
    EXPECT_DOUBLE_EQ(params.overlapSepLimit, 0.1) << "overlapSepLimit should be 0.1";

    // scaleMaxDisplace (key 15)
    std::cerr << "  -> Setting scaleMaxDisplace = \"100.0\"\n";
    params.set_value("100.0", ParamKeyword::scaleMaxDisplace);
    EXPECT_DOUBLE_EQ(params.scaleMaxDisplace, 100.0) << "scaleMaxDisplace should be 100.0";
}

/*!
 * \brief Test the Parameters::set_value function for boolean-valued keywords.
 *
 * These rely on read_boolean() from the parser functions to convert the string.
 */
void test_set_value_booleans()
{
    std::cerr << "\n[TEST] test_set_value_booleans()\n"
              << "  File under test : class_parameter.cpp\n"
              << "  Function        : Parameters::set_value()\n"
              << "  Purpose         : Verify boolean parameters parse correctly.\n";

    Parameters params;

    // fromRestart (key 3)
    std::cerr << "  -> Setting fromRestart = \"true\"\n";
    params.set_value("true", ParamKeyword::fromRestart);
    EXPECT_TRUE(params.fromRestart) << "fromRestart should be true";

    // clusterOverlapCheck (key 17)
    std::cerr << "  -> Setting clusterOverlapCheck = \"false\"\n";
    params.set_value("false", ParamKeyword::clusterOverlapCheck);
    EXPECT_FALSE(params.clusterOverlapCheck) << "clusterOverlapCheck should be false";

    // assocDissocWrite (key 18)
    std::cerr << "  -> Setting assocDissocWrite = \"true\"\n";
    params.set_value("true", ParamKeyword::assocDissocWrite);
    EXPECT_TRUE(params.assocDissocWrite) << "assocDissocWrite should be true";

    // rngwrite (key 19)
    std::cerr << "  -> Setting rngwrite = \"true\"\n";
    params.set_value("true", ParamKeyword::rngwrite);
    EXPECT_TRUE(params.rngwrite) << "rngwrite should be true";
}

/*!
 * \brief Test the Parameters::set_value function for string-valued keywords.
 */
void test_set_value_string()
{
    std::cerr << "\n[TEST] test_set_value_string()\n"
              << "  File under test : class_parameter.cpp\n"
              << "  Function        : Parameters::set_value()\n"
              << "  Purpose         : Verify string parameter (name) is stored correctly.\n";

    Parameters params;

    // name (key 13)
    std::cerr << "  -> Setting name = \"my_simulation\"\n";
    params.set_value("my_simulation", ParamKeyword::name);
    EXPECT_EQ(params.name, "my_simulation") << "name should be 'my_simulation'";
}

/*!
 * \brief Test Parameters::parse_paramFile by writing a temporary parameter
 *        block to a file and parsing it.
 *
 * Verifies that multiple keywords in a block are correctly read and that
 * parsing terminates when the "end parameters" line is encountered.
 */
void test_parse_paramFile()
{
    std::cerr << "\n[TEST] test_parse_paramFile()\n"
              << "  File under test : class_parameter.cpp\n"
              << "  Function        : Parameters::parse_paramFile()\n"
              << "  Purpose         : Verify a parameter block is parsed from an ifstream.\n";

    // Create a temporary parameter file with a small parameter block.
    const std::string fname = "test_param_block.tmp.inp";
    {
        std::ofstream out(fname);
        out << "nItr = 5000\n";               // long long
        out << "timeStep = 0.25\n";            // double
        out << "numMolTypes = 2\n";            // int
        out << "name = testrun\n";             // string
        out << "# this is a comment line\n";   // should be skipped
        out << "unknownKeyword = 42\n";        // should be ignored w/ warning
        out << "end parameters\n";             // terminates the block
        out.close();
    }

    // Open the file and parse it.
    std::ifstream in(fname);
    ASSERT_TRUE(in.is_open()) << "Failed to open temporary parameter file for reading";

    Parameters params;
    std::cerr << "  -> Parsing parameter block from " << fname << "\n";
    params.parse_paramFile(in);
    in.close();

    // Verify that the recognized keywords were parsed.
    EXPECT_EQ(params.nItr, 5000LL) << "nItr should be 5000 after parsing";
    EXPECT_DOUBLE_EQ(params.timeStep, 0.25) << "timeStep should be 0.25 after parsing";
    EXPECT_EQ(params.numMolTypes, 2) << "numMolTypes should be 2 after parsing";
    EXPECT_EQ(params.name, "testrun") << "name should be 'testrun' after parsing";

    // Clean up the temporary file.
    std::remove(fname.c_str());
}

/*!
 * \brief Test Parameters::display simply ensures it runs without crashing.
 *
 * display() only prints to std::cout and has no return value, so we merely
 * exercise it to confirm it does not throw or crash.
 */
void test_display()
{
    std::cerr << "\n[TEST] test_display()\n"
              << "  File under test : class_parameter.cpp\n"
              << "  Function        : Parameters::display()\n"
              << "  Purpose         : Verify display() runs without throwing.\n";

    Parameters params;
    // Populate a few fields so the output is meaningful.
    params.set_value("100", ParamKeyword::nItr);
    params.set_value("0.1", ParamKeyword::timeStep);

    std::cerr << "  -> Calling display() (output goes to std::cout)\n";
    // We wrap the call in EXPECT_NO_THROW to keep the test non-fatal.
    EXPECT_NO_THROW(params.display()) << "display() should not throw";
}

/* ---------------------------------------------------------------------------
 * GoogleTest wrappers.
 * Each TEST invokes one of the named test_* helper functions above so that all
 * assertions run even if some fail (EXPECT_* are non-fatal).
 * ------------------------------------------------------------------------- */

TEST(ParametersTest, SetValueIntegers)   { test_set_value_integers(); }
TEST(ParametersTest, SetValueIntervals)  { test_set_value_intervals(); }
TEST(ParametersTest, SetValueDoubles)    { test_set_value_doubles(); }
TEST(ParametersTest, SetValueBooleans)   { test_set_value_booleans(); }
TEST(ParametersTest, SetValueString)     { test_set_value_string(); }
TEST(ParametersTest, ParseParamFile)     { test_parse_paramFile(); }
TEST(ParametersTest, Display)            { test_display(); }