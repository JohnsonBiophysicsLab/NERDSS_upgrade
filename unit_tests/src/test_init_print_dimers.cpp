/*! \file test_init_print_dimers.cpp
 *
 * ### Unit test for ../src/io/init_print_dimers.cpp
 *
 * This test exercises the single function defined in that translation unit:
 *
 *     void init_print_dimers(std::ofstream& outfile,
 *                            Parameters params,
 *                            std::vector<MolTemplate>& molTemplateList)
 *
 * The function writes the *header line* of the MONO/DIMER observable output
 * file.  Its contract is:
 *
 *   1. It always begins the line with the literal text "TIME (s) " followed by
 *      a tab character.
 *   2. For every MolTemplate in molTemplateList that is **not** an implicit
 *      lipid it appends two tab separated column labels:
 *          "MONO:<molName>"     and     "DIMERS W:<molName>"
 *   3. Implicit lipid templates are skipped entirely (no columns emitted).
 *   4. The line is terminated with a newline (std::endl).
 *
 * Because the function only produces text, every test below writes to a
 * temporary file, reads the file back, and inspects the resulting string.
 * Verbose messages are printed to stderr so the reader can follow which
 * source file / function is being tested and what each assertion checks.
 */

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Parameters.hpp"
#include "io/io.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Helper: build a minimal MolTemplate with just the fields init_print_dimers
//         actually reads (molName and isImplicitLipid).
// -----------------------------------------------------------------------------
MolTemplate ipd_make_template(const std::string& name, bool isImplicitLipid)
{
    MolTemplate temp;
    temp.molName = name;
    temp.isImplicitLipid = isImplicitLipid;
    return temp;
}

// -----------------------------------------------------------------------------
// Helper: run init_print_dimers into a temporary file and return the full text
//         that was written.  The temporary file is removed afterwards so the
//         test leaves no artifacts behind.
// -----------------------------------------------------------------------------
std::string ipd_run_and_capture(std::vector<MolTemplate>& molTemplateList,
    const std::string& tmpFileName)
{
    Parameters params; // The function ignores params, defaults are fine.

    // --- write phase -------------------------------------------------------
    std::ofstream outFile(tmpFileName);
    // Report (non-fatally) if the temporary file could not be opened.
    EXPECT_TRUE(outFile.is_open())
        << "Could not open temporary output file: " << tmpFileName;
    if (!outFile.is_open())
        return std::string();

    std::cerr << "    calling init_print_dimers() with "
              << molTemplateList.size() << " MolTemplate(s)\n";
    init_print_dimers(outFile, params, molTemplateList);
    outFile.close();

    // --- read-back phase ---------------------------------------------------
    std::ifstream inFile(tmpFileName);
    EXPECT_TRUE(inFile.is_open())
        << "Could not re-open temporary output file for reading: " << tmpFileName;
    std::string contents;
    if (inFile.is_open()) {
        std::string line;
        while (std::getline(inFile, line)) {
            contents += line;
            contents += '\n';
        }
        inFile.close();
    }

    // Clean up the temporary file.
    std::remove(tmpFileName.c_str());

    std::cerr << "    captured text (tabs shown as <T>): \"";
    for (char c : contents) {
        if (c == '\t')
            std::cerr << "<T>";
        else if (c == '\n')
            std::cerr << "<NL>";
        else
            std::cerr << c;
    }
    std::cerr << "\"\n";

    return contents;
}

// -----------------------------------------------------------------------------
// Helper: count how many tab characters appear in a string.
// -----------------------------------------------------------------------------
size_t ipd_count_tabs(const std::string& s)
{
    return static_cast<size_t>(std::count(s.begin(), s.end(), '\t'));
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: An empty MolTemplate list should produce nothing but the time column
//         header plus a terminating newline.
// -----------------------------------------------------------------------------
void test_ipd_empty_template_list()
{
    std::cerr << "\n[TEST] test_ipd_empty_template_list\n"
              << "  Source file:   src/io/init_print_dimers.cpp\n"
              << "  Function:      init_print_dimers()\n"
              << "  Scenario:      molTemplateList is empty.\n"
              << "  Pass criteria: output is exactly \"TIME (s) \\t\" plus a newline,\n"
              << "                 i.e. one tab and no MONO/DIMERS columns.\n";

    std::vector<MolTemplate> molTemplateList; // deliberately empty

    const std::string text
        = ipd_run_and_capture(molTemplateList, "test_ipd_empty.tmp.dat");

    // The literal time header must be the very start of the output.
    EXPECT_EQ(text.rfind("TIME (s) ", 0), 0u)
        << "Output must begin with the literal \"TIME (s) \"";

    // Exactly one tab (the one right after "TIME (s) ").
    EXPECT_EQ(ipd_count_tabs(text), 1u)
        << "With no templates only the single time-column tab should be written";

    // No per-molecule column labels should exist.
    EXPECT_EQ(text.find("MONO:"), std::string::npos)
        << "No MONO: column should be written for an empty template list";
    EXPECT_EQ(text.find("DIMERS W:"), std::string::npos)
        << "No DIMERS W: column should be written for an empty template list";

    // The line must be newline terminated (std::endl).
    EXPECT_FALSE(text.empty()) << "Output should not be empty";
    if (!text.empty()) {
        EXPECT_EQ(text.back(), '\n')
            << "init_print_dimers() should terminate the header with a newline";
    }
}

// -----------------------------------------------------------------------------
// Test 2: A single, ordinary (non-implicit-lipid) template should yield one
//         MONO column and one DIMERS W column, both naming that molecule.
// -----------------------------------------------------------------------------
void test_ipd_single_template()
{
    std::cerr << "\n[TEST] test_ipd_single_template\n"
              << "  Source file:   src/io/init_print_dimers.cpp\n"
              << "  Function:      init_print_dimers()\n"
              << "  Scenario:      one explicit molecule template named \"clat\".\n"
              << "  Pass criteria: header contains MONO:clat and DIMERS W:clat and\n"
              << "                 has 3 tabs total (1 time + 2 molecule columns).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(ipd_make_template("clat", false));

    const std::string text
        = ipd_run_and_capture(molTemplateList, "test_ipd_single.tmp.dat");

    // Both expected column headers must be present.
    EXPECT_NE(text.find("MONO:clat"), std::string::npos)
        << "Expected a \"MONO:clat\" column label";
    EXPECT_NE(text.find("DIMERS W:clat"), std::string::npos)
        << "Expected a \"DIMERS W:clat\" column label";

    // 1 tab for the time column + 2 tabs for the two molecule columns.
    EXPECT_EQ(ipd_count_tabs(text), 3u)
        << "Expected 1 + 2*1 = 3 tab separators for one explicit template";

    // MONO must come before DIMERS W for the same molecule (column ordering).
    const size_t monoPos = text.find("MONO:clat");
    const size_t dimerPos = text.find("DIMERS W:clat");
    if (monoPos != std::string::npos && dimerPos != std::string::npos) {
        EXPECT_LT(monoPos, dimerPos)
            << "The MONO column must be written before the DIMERS W column";
    }
}

// -----------------------------------------------------------------------------
// Test 3: Several explicit templates -> columns for every molecule, in the same
//         order as the input vector.
// -----------------------------------------------------------------------------
void test_ipd_multiple_templates_order()
{
    std::cerr << "\n[TEST] test_ipd_multiple_templates_order\n"
              << "  Source file:   src/io/init_print_dimers.cpp\n"
              << "  Function:      init_print_dimers()\n"
              << "  Scenario:      three explicit templates: A, B, C.\n"
              << "  Pass criteria: all six column labels exist, they appear in the\n"
              << "                 input order, and there are 1 + 2*3 = 7 tabs.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(ipd_make_template("A", false));
    molTemplateList.push_back(ipd_make_template("B", false));
    molTemplateList.push_back(ipd_make_template("C", false));

    const std::string text
        = ipd_run_and_capture(molTemplateList, "test_ipd_multi.tmp.dat");

    // Every molecule must contribute both of its columns.
    const size_t monoA = text.find("MONO:A");
    const size_t dimA = text.find("DIMERS W:A");
    const size_t monoB = text.find("MONO:B");
    const size_t dimB = text.find("DIMERS W:B");
    const size_t monoC = text.find("MONO:C");
    const size_t dimC = text.find("DIMERS W:C");

    EXPECT_NE(monoA, std::string::npos) << "Missing MONO:A column";
    EXPECT_NE(dimA, std::string::npos) << "Missing DIMERS W:A column";
    EXPECT_NE(monoB, std::string::npos) << "Missing MONO:B column";
    EXPECT_NE(dimB, std::string::npos) << "Missing DIMERS W:B column";
    EXPECT_NE(monoC, std::string::npos) << "Missing MONO:C column";
    EXPECT_NE(dimC, std::string::npos) << "Missing DIMERS W:C column";

    // Column ordering must follow the order of molTemplateList.
    if (monoA != std::string::npos && monoB != std::string::npos)
        EXPECT_LT(monoA, monoB) << "Molecule A columns must precede molecule B columns";
    if (monoB != std::string::npos && monoC != std::string::npos)
        EXPECT_LT(monoB, monoC) << "Molecule B columns must precede molecule C columns";
    if (dimA != std::string::npos && monoB != std::string::npos)
        EXPECT_LT(dimA, monoB)
            << "Each molecule's DIMERS column must precede the next molecule's MONO column";

    // 1 time-column tab + 2 tabs per explicit molecule.
    EXPECT_EQ(ipd_count_tabs(text), 7u)
        << "Expected 1 + 2*3 = 7 tab separators for three explicit templates";
}

// -----------------------------------------------------------------------------
// Test 4: Implicit lipid templates must be skipped (the `continue` branch).
// -----------------------------------------------------------------------------
void test_ipd_skips_implicit_lipid()
{
    std::cerr << "\n[TEST] test_ipd_skips_implicit_lipid\n"
              << "  Source file:   src/io/init_print_dimers.cpp\n"
              << "  Function:      init_print_dimers()\n"
              << "  Scenario:      list is [pip2 (implicit lipid), ap2 (explicit)].\n"
              << "  Pass criteria: no columns are emitted for the implicit lipid,\n"
              << "                 the explicit molecule still gets both columns,\n"
              << "                 and the tab count is 1 + 2*1 = 3.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(ipd_make_template("pip2", true)); // implicit lipid
    molTemplateList.push_back(ipd_make_template("ap2", false)); // explicit

    const std::string text
        = ipd_run_and_capture(molTemplateList, "test_ipd_implicit.tmp.dat");

    // The implicit lipid must be completely absent from the header.
    EXPECT_EQ(text.find("MONO:pip2"), std::string::npos)
        << "Implicit lipid should not produce a MONO column";
    EXPECT_EQ(text.find("DIMERS W:pip2"), std::string::npos)
        << "Implicit lipid should not produce a DIMERS W column";

    // The explicit molecule is still reported.
    EXPECT_NE(text.find("MONO:ap2"), std::string::npos)
        << "Explicit molecule should still produce a MONO column";
    EXPECT_NE(text.find("DIMERS W:ap2"), std::string::npos)
        << "Explicit molecule should still produce a DIMERS W column";

    // Only the explicit molecule contributes column tabs.
    EXPECT_EQ(ipd_count_tabs(text), 3u)
        << "Expected 1 + 2*1 = 3 tabs when one of two templates is an implicit lipid";
}

// -----------------------------------------------------------------------------
// Test 5: A list made up solely of implicit lipids degenerates to the same
//         output as an empty list.
// -----------------------------------------------------------------------------
void test_ipd_all_implicit_lipids()
{
    std::cerr << "\n[TEST] test_ipd_all_implicit_lipids\n"
              << "  Source file:   src/io/init_print_dimers.cpp\n"
              << "  Function:      init_print_dimers()\n"
              << "  Scenario:      every template in the list is an implicit lipid.\n"
              << "  Pass criteria: output degenerates to just the time column\n"
              << "                 (one tab, no MONO/DIMERS labels).\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(ipd_make_template("pip2", true));
    molTemplateList.push_back(ipd_make_template("ps", true));

    const std::string text
        = ipd_run_and_capture(molTemplateList, "test_ipd_allimplicit.tmp.dat");

    // Header prefix still written.
    EXPECT_EQ(text.rfind("TIME (s) ", 0), 0u)
        << "Output must still begin with the literal \"TIME (s) \"";

    // No molecule columns at all.
    EXPECT_EQ(text.find("MONO:"), std::string::npos)
        << "No MONO: column should be written when all templates are implicit lipids";
    EXPECT_EQ(text.find("DIMERS W:"), std::string::npos)
        << "No DIMERS W: column should be written when all templates are implicit lipids";
    EXPECT_EQ(ipd_count_tabs(text), 1u)
        << "Only the single time-column tab should remain";
}

// -----------------------------------------------------------------------------
// Test 6: Calling the function twice on the same stream appends a second header
//         line; this checks that the routine writes exactly one line per call
//         and does not, for example, seek or truncate the stream.
// -----------------------------------------------------------------------------
void test_ipd_two_calls_produce_two_lines()
{
    std::cerr << "\n[TEST] test_ipd_two_calls_produce_two_lines\n"
              << "  Source file:   src/io/init_print_dimers.cpp\n"
              << "  Function:      init_print_dimers()\n"
              << "  Scenario:      the function is invoked twice on one open stream.\n"
              << "  Pass criteria: the file holds exactly two identical header lines,\n"
              << "                 proving one newline-terminated line per call.\n";

    Parameters params;
    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(ipd_make_template("gag", false));

    const std::string tmpFileName = "test_ipd_twice.tmp.dat";

    std::ofstream outFile(tmpFileName);
    EXPECT_TRUE(outFile.is_open()) << "Could not open " << tmpFileName;
    if (outFile.is_open()) {
        std::cerr << "    calling init_print_dimers() twice on the same stream\n";
        init_print_dimers(outFile, params, molTemplateList);
        init_print_dimers(outFile, params, molTemplateList);
        outFile.close();
    }

    // Read the file line by line so we can count the header lines.
    std::vector<std::string> lines;
    std::ifstream inFile(tmpFileName);
    EXPECT_TRUE(inFile.is_open()) << "Could not re-open " << tmpFileName;
    if (inFile.is_open()) {
        std::string line;
        while (std::getline(inFile, line))
            lines.push_back(line);
        inFile.close();
    }
    std::remove(tmpFileName.c_str());

    std::cerr << "    number of lines written: " << lines.size() << "\n";

    // Two calls => two lines.
    EXPECT_EQ(lines.size(), 2u)
        << "Each call to init_print_dimers() should emit exactly one line";

    // Both lines must be identical and both must carry the molecule columns.
    if (lines.size() == 2u) {
        EXPECT_EQ(lines[0], lines[1])
            << "Repeated calls with identical input should produce identical headers";
        EXPECT_NE(lines[0].find("MONO:gag"), std::string::npos)
            << "First header line should contain MONO:gag";
        EXPECT_NE(lines[1].find("DIMERS W:gag"), std::string::npos)
            << "Second header line should contain DIMERS W:gag";
        EXPECT_EQ(ipd_count_tabs(lines[0]), 3u)
            << "Each header line should carry 1 + 2*1 = 3 tab separators";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* helper is executed inside its own TEST
// so failures are reported individually while every case still runs.
// -----------------------------------------------------------------------------
TEST(InitPrintDimers, EmptyTemplateList) { test_ipd_empty_template_list(); }
TEST(InitPrintDimers, SingleTemplate) { test_ipd_single_template(); }
TEST(InitPrintDimers, MultipleTemplatesOrder) { test_ipd_multiple_templates_order(); }
TEST(InitPrintDimers, SkipsImplicitLipid) { test_ipd_skips_implicit_lipid(); }
TEST(InitPrintDimers, AllImplicitLipids) { test_ipd_all_implicit_lipids(); }
TEST(InitPrintDimers, TwoCallsProduceTwoLines) { test_ipd_two_calls_produce_two_lines(); }