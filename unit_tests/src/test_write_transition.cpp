/*! \file test_write_transition.cpp
 *
 * ### Unit test for src/io/write_transition.cpp
 *
 * Function under test:
 *
 *     void write_transition(double time, std::ofstream& transitionFile,
 *                           const std::vector<MolTemplate>& molTemplateList)
 *
 * The routine dumps, for every MolTemplate that has `countTransition == true`:
 *
 *   1. a "time: <t>" header line,
 *   2. a "transition matrix for each mol type: " section containing, per
 *      molecule name, `transitionMatrixSize` rows of `transitionMatrixSize`
 *      space-prefixed integers taken from MolTemplate::transitionMatrix,
 *   3. a "lifetime for each mol type: " section containing, per molecule name
 *      and per cluster size, all doubles held in MolTemplate::lifeTime.
 *
 * MolTemplates whose `countTransition` flag is false must be skipped entirely.
 *
 * Strategy: write to a temporary file, read it back line by line, and compare
 * the produced text against the exact expected layout.  All checks use
 * non-fatal EXPECT_* assertions so every scenario runs to completion.
 */

#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Helper: build a minimal MolTemplate suitable for write_transition().
//
// The transition matrix is allocated as a size x size matrix of zeros and the
// lifeTime container is allocated with `size` (initially empty) rows, exactly
// as write_transition() expects to find them.
// -----------------------------------------------------------------------------
MolTemplate wt_make_mol_template(const std::string& name, bool countTransition, int size)
{
    MolTemplate molTemp;
    molTemp.molName = name;
    molTemp.countTransition = countTransition;
    molTemp.transitionMatrixSize = size;
    molTemp.transitionMatrix.assign(size, std::vector<long long int>(size, 0));
    molTemp.lifeTime.assign(size, std::vector<double>());
    return molTemp;
}

// -----------------------------------------------------------------------------
// Helper: call write_transition() into a temporary file and read the file back
// as a vector of lines.  The temporary file is removed afterwards.
// -----------------------------------------------------------------------------
std::vector<std::string> wt_run_and_read(const std::string& fileName, double time,
    const std::vector<MolTemplate>& molTemplateList)
{
    // Write phase: open the stream, invoke the function under test, close it.
    std::ofstream transitionFile(fileName);
    if (!transitionFile.is_open()) {
        std::cerr << "  !! Could not open temporary file '" << fileName << "' for writing\n";
        return {};
    }
    write_transition(time, transitionFile, molTemplateList);
    transitionFile.close();

    // Read phase: slurp the file back one line at a time.
    std::vector<std::string> lines;
    std::ifstream inFile(fileName);
    std::string line;
    while (std::getline(inFile, line))
        lines.push_back(line);
    inFile.close();

    // Clean up the scratch file so repeated test runs stay hermetic.
    std::remove(fileName.c_str());

    // Echo the captured output so a human reader can see what was produced.
    std::cerr << "  Captured " << lines.size() << " line(s) of output:\n";
    for (std::size_t i = 0; i < lines.size(); ++i)
        std::cerr << "    [" << i << "] \"" << lines[i] << "\"\n";

    return lines;
}

// -----------------------------------------------------------------------------
// Helper: safe line accessor.  Returns a sentinel instead of throwing so that a
// short/incorrect output cannot abort the whole test binary.
// -----------------------------------------------------------------------------
std::string wt_line(const std::vector<std::string>& lines, std::size_t idx)
{
    return (idx < lines.size()) ? lines[idx] : std::string("<MISSING LINE>");
}

// -----------------------------------------------------------------------------
// Helper: does any line contain the given substring?
// -----------------------------------------------------------------------------
bool wt_contains(const std::vector<std::string>& lines, const std::string& needle)
{
    for (const auto& line : lines) {
        if (line.find(needle) != std::string::npos)
            return true;
    }
    return false;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: an empty molTemplateList must still emit the three fixed header lines
//         (time, transition matrix section, lifetime section) and nothing else.
// -----------------------------------------------------------------------------
void test_wt_empty_mol_template_list()
{
    std::cerr << "\n[TEST] test_wt_empty_mol_template_list\n"
              << "  Source file: src/io/write_transition.cpp\n"
              << "  Function:    write_transition()\n"
              << "  Scenario:    molTemplateList is empty.\n"
              << "  Pass:        exactly 3 lines -> time header + 2 section headers.\n";

    std::vector<MolTemplate> molTemplateList; // deliberately empty

    std::vector<std::string> lines
        = wt_run_and_read("test_write_transition_empty.tmp", 0.0, molTemplateList);

    // Only the three unconditional header lines should be present.
    EXPECT_EQ(lines.size(), 3u) << "Empty template list should produce only 3 header lines";

    // The time header must carry the value that was passed in (0 prints as "0").
    EXPECT_EQ(wt_line(lines, 0), "time: 0") << "First line should be the time header";

    // Section headers (trailing space is part of the literal, so use find()).
    EXPECT_NE(wt_line(lines, 1).find("transition matrix for each mol type:"), std::string::npos)
        << "Second line should be the transition matrix section header";
    EXPECT_NE(wt_line(lines, 2).find("lifetime for each mol type:"), std::string::npos)
        << "Third line should be the lifetime section header";
}

// -----------------------------------------------------------------------------
// Test 2: templates with countTransition == false must be skipped entirely, so
//         neither their name nor any matrix row may appear in the output.
// -----------------------------------------------------------------------------
void test_wt_skips_templates_without_count_transition()
{
    std::cerr << "\n[TEST] test_wt_skips_templates_without_count_transition\n"
              << "  Source file: src/io/write_transition.cpp\n"
              << "  Function:    write_transition()\n"
              << "  Scenario:    one MolTemplate with countTransition == false.\n"
              << "  Pass:        output is header-only; the molecule name is absent.\n";

    // Build a template that carries real data but has counting switched OFF.
    MolTemplate skipped = wt_make_mol_template("SKIPPED_MOL", false, 2);
    skipped.transitionMatrix[0][0] = 111;
    skipped.transitionMatrix[1][1] = 222;
    skipped.lifeTime[0].push_back(9.25);

    std::vector<MolTemplate> molTemplateList { skipped };

    std::vector<std::string> lines
        = wt_run_and_read("test_write_transition_skip.tmp", 2.0, molTemplateList);

    // Same three header lines as the empty case: nothing was emitted for the mol.
    EXPECT_EQ(lines.size(), 3u) << "A non-counted MolTemplate must contribute no lines";

    // Neither the name nor any of its data may leak into the file.
    EXPECT_FALSE(wt_contains(lines, "SKIPPED_MOL"))
        << "Molecule name should not be written when countTransition is false";
    EXPECT_FALSE(wt_contains(lines, "111"))
        << "Transition matrix data should not be written when countTransition is false";
    EXPECT_FALSE(wt_contains(lines, "9.25"))
        << "Lifetime data should not be written when countTransition is false";
    EXPECT_FALSE(wt_contains(lines, "size of the cluster:"))
        << "No cluster-size lines should be written when countTransition is false";
}

// -----------------------------------------------------------------------------
// Test 3: full round-trip for a single counted MolTemplate.  Verifies the exact
//         line-by-line layout: matrix rows, cluster-size markers, and lifetimes
//         (including an empty lifetime row which must still emit a blank line).
// -----------------------------------------------------------------------------
void test_wt_single_template_matrix_and_lifetime()
{
    std::cerr << "\n[TEST] test_wt_single_template_matrix_and_lifetime\n"
              << "  Source file: src/io/write_transition.cpp\n"
              << "  Function:    write_transition()\n"
              << "  Scenario:    one counted MolTemplate, 2x2 matrix {{1,2},{3,4}},\n"
              << "               lifeTime[0] = {0.5, 1.5}, lifeTime[1] = {} (empty).\n"
              << "  Pass:        11 lines matching the documented layout exactly.\n";

    MolTemplate molTemp = wt_make_mol_template("A", true, 2);
    // Fill the transition matrix with easily recognisable integers.
    molTemp.transitionMatrix[0][0] = 1;
    molTemp.transitionMatrix[0][1] = 2;
    molTemp.transitionMatrix[1][0] = 3;
    molTemp.transitionMatrix[1][1] = 4;
    // Cluster size 1 has two recorded lifetimes; cluster size 2 has none.
    molTemp.lifeTime[0].push_back(0.5);
    molTemp.lifeTime[0].push_back(1.5);

    std::vector<MolTemplate> molTemplateList { molTemp };

    std::vector<std::string> lines
        = wt_run_and_read("test_write_transition_single.tmp", 1.5, molTemplateList);

    // Expected layout:
    //   0: time: 1.5
    //   1: transition matrix for each mol type:
    //   2: A
    //   3:  1 2
    //   4:  3 4
    //   5: lifetime for each mol type:
    //   6: A
    //   7: size of the cluster:1
    //   8:  0.5 1.5
    //   9: size of the cluster:2
    //  10: (empty, from the lifetime row with no entries)
    EXPECT_EQ(lines.size(), 11u) << "Expected 11 output lines for this configuration";

    EXPECT_EQ(wt_line(lines, 0), "time: 1.5") << "Time header should print the supplied time";
    EXPECT_NE(wt_line(lines, 1).find("transition matrix for each mol type:"), std::string::npos)
        << "Line 1 should open the transition matrix section";
    EXPECT_EQ(wt_line(lines, 2), "A") << "Line 2 should be the molecule name";

    // Each matrix element is written as ' ' followed by the value.
    EXPECT_EQ(wt_line(lines, 3), " 1 2") << "First matrix row should be ' 1 2'";
    EXPECT_EQ(wt_line(lines, 4), " 3 4") << "Second matrix row should be ' 3 4'";

    EXPECT_NE(wt_line(lines, 5).find("lifetime for each mol type:"), std::string::npos)
        << "Line 5 should open the lifetime section";
    EXPECT_EQ(wt_line(lines, 6), "A") << "Line 6 should repeat the molecule name in lifetime section";

    EXPECT_EQ(wt_line(lines, 7), "size of the cluster:1")
        << "Cluster indices are written 1-based (indexOne + 1)";
    EXPECT_EQ(wt_line(lines, 8), " 0.5 1.5") << "Lifetimes for cluster size 1 should be ' 0.5 1.5'";

    EXPECT_EQ(wt_line(lines, 9), "size of the cluster:2")
        << "Second cluster-size marker should be present";
    EXPECT_EQ(wt_line(lines, 10), "")
        << "A cluster size with no recorded lifetimes should still emit a blank line";
}

// -----------------------------------------------------------------------------
// Test 4: a 1x1 matrix is the smallest meaningful case; confirm the single row
//         and single cluster-size block are produced correctly.
// -----------------------------------------------------------------------------
void test_wt_single_element_matrix()
{
    std::cerr << "\n[TEST] test_wt_single_element_matrix\n"
              << "  Source file: src/io/write_transition.cpp\n"
              << "  Function:    write_transition()\n"
              << "  Scenario:    transitionMatrixSize == 1, one lifetime value.\n"
              << "  Pass:        one matrix row and one cluster-size block.\n";

    MolTemplate molTemp = wt_make_mol_template("MONO", true, 1);
    molTemp.transitionMatrix[0][0] = 42;
    molTemp.lifeTime[0].push_back(7.25);

    std::vector<MolTemplate> molTemplateList { molTemp };

    std::vector<std::string> lines
        = wt_run_and_read("test_write_transition_mono.tmp", 3.0, molTemplateList);

    // Expected: time, matrix header, name, " 42", lifetime header, name,
    //           "size of the cluster:1", " 7.25"  -> 8 lines.
    EXPECT_EQ(lines.size(), 8u) << "A 1x1 matrix with one lifetime should produce 8 lines";

    EXPECT_EQ(wt_line(lines, 0), "time: 3") << "Time 3.0 should print as '3'";
    EXPECT_EQ(wt_line(lines, 2), "MONO") << "Molecule name should head the matrix block";
    EXPECT_EQ(wt_line(lines, 3), " 42") << "Single matrix element should be written as ' 42'";
    EXPECT_EQ(wt_line(lines, 5), "MONO") << "Molecule name should head the lifetime block";
    EXPECT_EQ(wt_line(lines, 6), "size of the cluster:1")
        << "Only one cluster-size marker is expected";
    EXPECT_EQ(wt_line(lines, 7), " 7.25") << "Lifetime value should be written as ' 7.25'";
}

// -----------------------------------------------------------------------------
// Test 5: several templates, mixing counted and non-counted ones.  Confirms
//         section ordering (all matrices first, then all lifetimes) and that
//         only the counted templates appear, in their original list order.
// -----------------------------------------------------------------------------
void test_wt_multiple_templates_ordering()
{
    std::cerr << "\n[TEST] test_wt_multiple_templates_ordering\n"
              << "  Source file: src/io/write_transition.cpp\n"
              << "  Function:    write_transition()\n"
              << "  Scenario:    templates [FIRST(count), MIDDLE(no count), LAST(count)].\n"
              << "  Pass:        matrices of FIRST/LAST precede all lifetime blocks and\n"
              << "               MIDDLE never appears.\n";

    // FIRST: counted, 1x1 matrix.
    MolTemplate first = wt_make_mol_template("FIRST", true, 1);
    first.transitionMatrix[0][0] = 5;
    first.lifeTime[0].push_back(1.0);

    // MIDDLE: not counted, must be invisible in the output.
    MolTemplate middle = wt_make_mol_template("MIDDLE", false, 1);
    middle.transitionMatrix[0][0] = 999;

    // LAST: counted, 2x2 matrix with no lifetime entries at all.
    MolTemplate last = wt_make_mol_template("LAST", true, 2);
    last.transitionMatrix[0][0] = 6;
    last.transitionMatrix[0][1] = 7;
    last.transitionMatrix[1][0] = 8;
    last.transitionMatrix[1][1] = 9;

    std::vector<MolTemplate> molTemplateList { first, middle, last };

    std::vector<std::string> lines
        = wt_run_and_read("test_write_transition_multi.tmp", 10.0, molTemplateList);

    // Expected layout:
    //   0: time: 10
    //   1: transition matrix for each mol type:
    //   2: FIRST
    //   3:  5
    //   4: LAST
    //   5:  6 7
    //   6:  8 9
    //   7: lifetime for each mol type:
    //   8: FIRST
    //   9: size of the cluster:1
    //  10:  1
    //  11: LAST
    //  12: size of the cluster:1
    //  13: (blank)
    //  14: size of the cluster:2
    //  15: (blank)
    EXPECT_EQ(lines.size(), 16u) << "Expected 16 lines for the mixed three-template case";

    EXPECT_EQ(wt_line(lines, 0), "time: 10") << "Time 10.0 should print as '10'";

    // Matrix section: FIRST then LAST, preserving the input order.
    EXPECT_EQ(wt_line(lines, 2), "FIRST") << "FIRST should be the first counted molecule";
    EXPECT_EQ(wt_line(lines, 3), " 5") << "FIRST matrix row should be ' 5'";
    EXPECT_EQ(wt_line(lines, 4), "LAST") << "LAST should follow FIRST in the matrix section";
    EXPECT_EQ(wt_line(lines, 5), " 6 7") << "LAST first matrix row should be ' 6 7'";
    EXPECT_EQ(wt_line(lines, 6), " 8 9") << "LAST second matrix row should be ' 8 9'";

    // Lifetime section starts only after all matrices were written.
    EXPECT_NE(wt_line(lines, 7).find("lifetime for each mol type:"), std::string::npos)
        << "Lifetime section header must come after every matrix block";
    EXPECT_EQ(wt_line(lines, 8), "FIRST") << "Lifetime section should also start with FIRST";
    EXPECT_EQ(wt_line(lines, 9), "size of the cluster:1") << "FIRST has a single cluster size";
    EXPECT_EQ(wt_line(lines, 10), " 1") << "FIRST lifetime value should be ' 1'";
    EXPECT_EQ(wt_line(lines, 11), "LAST") << "LAST lifetime block should follow FIRST";
    EXPECT_EQ(wt_line(lines, 12), "size of the cluster:1") << "LAST cluster size 1 marker";
    EXPECT_EQ(wt_line(lines, 13), "") << "LAST cluster size 1 has no lifetimes -> blank line";
    EXPECT_EQ(wt_line(lines, 14), "size of the cluster:2") << "LAST cluster size 2 marker";
    EXPECT_EQ(wt_line(lines, 15), "") << "LAST cluster size 2 has no lifetimes -> blank line";

    // The non-counted template must be completely absent.
    EXPECT_FALSE(wt_contains(lines, "MIDDLE"))
        << "MIDDLE has countTransition == false and must not be written";
    EXPECT_FALSE(wt_contains(lines, "999"))
        << "MIDDLE matrix data must not be written";
}

// -----------------------------------------------------------------------------
// Test 6: the `time` argument is written verbatim via operator<<; check both a
//         fractional and a negative value to ensure nothing else is injected.
// -----------------------------------------------------------------------------
void test_wt_time_value_is_written()
{
    std::cerr << "\n[TEST] test_wt_time_value_is_written\n"
              << "  Source file: src/io/write_transition.cpp\n"
              << "  Function:    write_transition()\n"
              << "  Scenario:    call with time = 0.125 and then time = -2.5.\n"
              << "  Pass:        the first line reproduces the value via default\n"
              << "               ostream formatting.\n";

    std::vector<MolTemplate> molTemplateList; // no molecules needed for this check

    // Fractional positive time.
    std::vector<std::string> linesA
        = wt_run_and_read("test_write_transition_timeA.tmp", 0.125, molTemplateList);
    EXPECT_EQ(wt_line(linesA, 0), "time: 0.125")
        << "Fractional time should be printed with default precision";

    // Negative time (not physical, but the function must not filter it).
    std::vector<std::string> linesB
        = wt_run_and_read("test_write_transition_timeB.tmp", -2.5, molTemplateList);
    EXPECT_EQ(wt_line(linesB, 0), "time: -2.5")
        << "Negative time should be printed unchanged";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers - each scenario is its own TEST so failures are reported
// individually while every scenario still executes.
// -----------------------------------------------------------------------------
TEST(WriteTransition, EmptyMolTemplateList) { test_wt_empty_mol_template_list(); }
TEST(WriteTransition, SkipsTemplatesWithoutCountTransition) { test_wt_skips_templates_without_count_transition(); }
TEST(WriteTransition, SingleTemplateMatrixAndLifetime) { test_wt_single_template_matrix_and_lifetime(); }
TEST(WriteTransition, SingleElementMatrix) { test_wt_single_element_matrix(); }
TEST(WriteTransition, MultipleTemplatesOrdering) { test_wt_multiple_templates_ordering(); }
TEST(WriteTransition, TimeValueIsWritten) { test_wt_time_value_is_written(); }