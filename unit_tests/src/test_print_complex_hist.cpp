/*! \file test_print_complex_hist.cpp
 *
 * ### Unit test for src/io/print_complex_hist.cpp
 *
 * Function under test:
 *
 *     double print_complex_hist(std::vector<Complex>& complexList,
 *                               std::ofstream& outfile, int it,
 *                               Parameters params,
 *                               std::vector<MolTemplate>& molTemplateList,
 *                               std::vector<Molecule>& moleculeList,
 *                               int nImplicitLipids);
 *
 * The routine builds a histogram of the distinct complex "compositions"
 * present in the system (how many of each molecule type each complex
 * contains), and writes the result to the supplied output stream in the
 * following format:
 *
 *     Time (s): <t>
 *     <count>\t<MolName>: <n>. <MolName2>: <m>. \n
 *     ...
 *
 * Because the only observable effects of the function are (a) the text it
 * writes to the stream and (b) its (always 0.0) return value, every test
 * below writes to a temporary file, reads the file back into a string, and
 * asserts on the presence / absence of specific substrings.
 *
 * Notes on behaviour that is deliberately covered here:
 *   - complexes flagged isEmpty are skipped entirely,
 *   - identical compositions are grouped and counted,
 *   - molecule types with a zero count are omitted from the label,
 *   - a complex whose label would be empty is not counted at all,
 *   - when molTemplateList[0].isImplicitLipid is true, the count of the
 *     first molecule type is replaced by Complex::linksToSurface and an
 *     extra line reporting the number of free implicit lipids is emitted,
 *   - the time stamp uses (it - itrRestartFrom)*timeStep*1e-6 + timeRestartFrom.
 */

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "io/io.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Local helpers (prefixed with pch_ so they cannot collide with other tests).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal MolTemplate that only carries a name/index.
 *
 * print_complex_hist() only ever reads molName and isImplicitLipid from the
 * templates, so nothing else needs to be initialized.
 */
MolTemplate pch_make_mol_template(const std::string& name, int typeIndex)
{
    MolTemplate temp;
    temp.molName = name;
    temp.molTypeIndex = typeIndex;
    temp.isImplicitLipid = false;
    return temp;
}

/*! \brief Build a minimal Molecule belonging to a given complex/type. */
Molecule pch_make_molecule(int index, int molTypeIndex, int comIndex)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = comIndex;
    mol.isEmpty = false;
    return mol;
}

/*! \brief Build a Complex owning the supplied moleculeList indices. */
Complex pch_make_complex(int index, const std::vector<int>& members)
{
    Complex com;
    com.index = index;
    com.memberList = members;
    com.isEmpty = false;
    com.linksToSurface = 0;
    return com;
}

/*! \brief Slurp an entire text file into a std::string. */
std::string pch_read_file(const std::string& fileName)
{
    std::ifstream in(fileName);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

/*! \brief Count non-overlapping occurrences of needle within haystack. */
int pch_count_occurrences(const std::string& haystack, const std::string& needle)
{
    if (needle.empty())
        return 0;

    int count = 0;
    std::string::size_type pos = haystack.find(needle);
    while (pos != std::string::npos) {
        ++count;
        pos = haystack.find(needle, pos + needle.size());
    }
    return count;
}

/*! \brief Convenience wrapper: run the function, close and read back the file. */
std::string pch_run_and_capture(const std::string& fileName, std::vector<Complex>& complexList,
    int it, Parameters params, std::vector<MolTemplate>& molTemplateList,
    std::vector<Molecule>& moleculeList, int nImplicitLipids, double& returnValue)
{
    std::ofstream outFile(fileName, std::ios::trunc);
    returnValue = print_complex_hist(
        complexList, outFile, it, params, molTemplateList, moleculeList, nImplicitLipids);
    outFile.close();

    const std::string contents = pch_read_file(fileName);
    std::remove(fileName.c_str());
    return contents;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: A single complex made of 2 "A" molecules and 1 "B" molecule.
//         Verifies the composition label and the count column.
// -----------------------------------------------------------------------------
void test_pch_single_complex_composition()
{
    std::cerr << "\n[TEST] test_pch_single_complex_composition\n"
              << "  Source file:   src/io/print_complex_hist.cpp\n"
              << "  Function:      print_complex_hist()\n"
              << "  Scenario:      one complex containing 2 x A and 1 x B.\n"
              << "  Pass criteria: file contains a time header and the line\n"
              << "                 \"1<TAB>A: 2. B: 1. \"; return value is 0.0.\n";

    // Two molecule types: A (index 0) and B (index 1).
    std::vector<MolTemplate> molTemplateList {
        pch_make_mol_template("A", 0),
        pch_make_mol_template("B", 1)
    };

    // Three molecules, all in complex 0: A, A, B.
    std::vector<Molecule> moleculeList {
        pch_make_molecule(0, 0, 0),
        pch_make_molecule(1, 0, 0),
        pch_make_molecule(2, 1, 0)
    };

    // One complex owning all three molecules.
    std::vector<Complex> complexList { pch_make_complex(0, { 0, 1, 2 }) };

    // Simple timing parameters -> time stamp is 100 * 1.0 * 1e-6 = 0.0001 s.
    Parameters params;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    double retVal = -1.0;
    const std::string out = pch_run_and_capture(
        "test_pch_single.dat", complexList, 100, params, molTemplateList, moleculeList, 0, retVal);

    std::cerr << "  Captured output:\n" << out;

    // The function always returns 0.0 (it is a legacy signature).
    EXPECT_DOUBLE_EQ(retVal, 0.0) << "print_complex_hist should always return 0.0";

    // The header must be present with the computed simulation time.
    EXPECT_NE(out.find("Time (s): "), std::string::npos)
        << "Output must start with a \"Time (s): \" header";
    EXPECT_NE(out.find("0.0001"), std::string::npos)
        << "Time stamp should be (100 - 0) * 1.0 us = 0.0001 s";

    // The composition histogram line: one complex with A: 2. B: 1.
    EXPECT_NE(out.find("1\tA: 2. B: 1. "), std::string::npos)
        << "Expected the single complex to be reported as \"1\\tA: 2. B: 1. \"";
}

// -----------------------------------------------------------------------------
// Test 2: Identical complexes must be grouped and counted together, and
//         molecule types absent from a complex must be omitted from its label.
// -----------------------------------------------------------------------------
void test_pch_groups_identical_complexes()
{
    std::cerr << "\n[TEST] test_pch_groups_identical_complexes\n"
              << "  Source file:   src/io/print_complex_hist.cpp\n"
              << "  Function:      print_complex_hist()\n"
              << "  Scenario:      two identical A-B dimers plus one lone A monomer.\n"
              << "  Pass criteria: \"2<TAB>A: 1. B: 1. \" and \"1<TAB>A: 1. \" appear,\n"
              << "                 and the monomer label does not mention B.\n";

    std::vector<MolTemplate> molTemplateList {
        pch_make_mol_template("A", 0),
        pch_make_mol_template("B", 1)
    };

    // Molecules 0..3 form two A-B dimers, molecule 4 is a free A monomer.
    std::vector<Molecule> moleculeList {
        pch_make_molecule(0, 0, 0),
        pch_make_molecule(1, 1, 0),
        pch_make_molecule(2, 0, 1),
        pch_make_molecule(3, 1, 1),
        pch_make_molecule(4, 0, 2)
    };

    std::vector<Complex> complexList {
        pch_make_complex(0, { 0, 1 }),
        pch_make_complex(1, { 2, 3 }),
        pch_make_complex(2, { 4 })
    };

    Parameters params;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    double retVal = -1.0;
    const std::string out = pch_run_and_capture(
        "test_pch_group.dat", complexList, 10, params, molTemplateList, moleculeList, 0, retVal);

    std::cerr << "  Captured output:\n" << out;

    EXPECT_DOUBLE_EQ(retVal, 0.0) << "print_complex_hist should always return 0.0";

    // Two identical dimers must be collapsed into a single line with count 2.
    EXPECT_NE(out.find("2\tA: 1. B: 1. "), std::string::npos)
        << "The two identical A-B dimers should be reported with a count of 2";

    // The lone monomer gets its own line, and must not list molecule type B.
    EXPECT_NE(out.find("1\tA: 1. "), std::string::npos)
        << "The free A monomer should be reported with a count of 1";
    EXPECT_EQ(out.find("A: 1. B: 0."), std::string::npos)
        << "Molecule types with zero copies must be omitted from the label";

    // Exactly three data lines' worth of newlines: header + two histogram rows.
    EXPECT_EQ(pch_count_occurrences(out, "\n"), 3)
        << "Expected exactly one header line plus two histogram lines";
}

// -----------------------------------------------------------------------------
// Test 3: Complexes flagged isEmpty are skipped entirely.
// -----------------------------------------------------------------------------
void test_pch_skips_empty_complexes()
{
    std::cerr << "\n[TEST] test_pch_skips_empty_complexes\n"
              << "  Source file:   src/io/print_complex_hist.cpp\n"
              << "  Function:      print_complex_hist()\n"
              << "  Scenario:      two complexes, the second marked isEmpty.\n"
              << "  Pass criteria: only one histogram line is written, with count 1.\n";

    std::vector<MolTemplate> molTemplateList { pch_make_mol_template("A", 0) };

    std::vector<Molecule> moleculeList {
        pch_make_molecule(0, 0, 0),
        pch_make_molecule(1, 0, 1)
    };

    std::vector<Complex> complexList {
        pch_make_complex(0, { 0 }),
        pch_make_complex(1, { 1 })
    };
    // Mark the second complex as a destroyed / void slot.
    complexList[1].isEmpty = true;

    Parameters params;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    double retVal = -1.0;
    const std::string out = pch_run_and_capture(
        "test_pch_empty.dat", complexList, 1, params, molTemplateList, moleculeList, 0, retVal);

    std::cerr << "  Captured output:\n" << out;

    EXPECT_DOUBLE_EQ(retVal, 0.0) << "print_complex_hist should always return 0.0";

    // Only the surviving complex should be counted.
    EXPECT_NE(out.find("1\tA: 1. "), std::string::npos)
        << "The single non-empty complex should be reported with count 1";
    EXPECT_EQ(out.find("2\tA: 1. "), std::string::npos)
        << "The isEmpty complex must not contribute to the histogram count";

    // header + one histogram line == two newlines.
    EXPECT_EQ(pch_count_occurrences(out, "\n"), 2)
        << "Expected one header line plus exactly one histogram line";
}

// -----------------------------------------------------------------------------
// Test 4: Complexes with no members produce an empty label and are dropped,
//         so only the time header should be written.
// -----------------------------------------------------------------------------
void test_pch_no_countable_complexes()
{
    std::cerr << "\n[TEST] test_pch_no_countable_complexes\n"
              << "  Source file:   src/io/print_complex_hist.cpp\n"
              << "  Function:      print_complex_hist()\n"
              << "  Scenario:      (a) an empty complexList and (b) a complex\n"
              << "                 with an empty memberList.\n"
              << "  Pass criteria: only the \"Time (s): \" header is written.\n";

    std::vector<MolTemplate> molTemplateList { pch_make_mol_template("A", 0) };
    std::vector<Molecule> moleculeList {};

    Parameters params;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    // (a) Completely empty complex list.
    std::vector<Complex> noComplexes {};
    double retVal = -1.0;
    std::string out = pch_run_and_capture(
        "test_pch_none_a.dat", noComplexes, 0, params, molTemplateList, moleculeList, 0, retVal);

    std::cerr << "  Captured output for an empty complexList:\n" << out;
    EXPECT_DOUBLE_EQ(retVal, 0.0) << "print_complex_hist should always return 0.0";
    EXPECT_NE(out.find("Time (s): "), std::string::npos)
        << "The header must be written even when there are no complexes";
    EXPECT_EQ(pch_count_occurrences(out, "\n"), 1)
        << "Only the header line should be written for an empty complexList";

    // (b) A complex that exists but owns no molecules -> empty label -> skipped.
    std::vector<Complex> memberlessComplex { pch_make_complex(0, {}) };
    retVal = -1.0;
    out = pch_run_and_capture("test_pch_none_b.dat", memberlessComplex, 0, params,
        molTemplateList, moleculeList, 0, retVal);

    std::cerr << "  Captured output for a member-less complex:\n" << out;
    EXPECT_DOUBLE_EQ(retVal, 0.0) << "print_complex_hist should always return 0.0";
    EXPECT_EQ(pch_count_occurrences(out, "\n"), 1)
        << "A complex with an empty composition label must not be counted";
    EXPECT_EQ(out.find("A:"), std::string::npos)
        << "No composition label should be emitted for a member-less complex";
}

// -----------------------------------------------------------------------------
// Test 5: Restart bookkeeping - the time stamp is computed relative to
//         itrRestartFrom and offset by timeRestartFrom.
// -----------------------------------------------------------------------------
void test_pch_time_header_with_restart()
{
    std::cerr << "\n[TEST] test_pch_time_header_with_restart\n"
              << "  Source file:   src/io/print_complex_hist.cpp\n"
              << "  Function:      print_complex_hist()\n"
              << "  Scenario:      it=200, itrRestartFrom=100, timeStep=2 us,\n"
              << "                 timeRestartFrom=0.5 s.\n"
              << "  Pass criteria: header reports (200-100)*2e-6 + 0.5 = 0.5002 s.\n";

    std::vector<MolTemplate> molTemplateList { pch_make_mol_template("A", 0) };
    std::vector<Molecule> moleculeList { pch_make_molecule(0, 0, 0) };
    std::vector<Complex> complexList { pch_make_complex(0, { 0 }) };

    Parameters params;
    params.timeStep = 2.0;          // microseconds
    params.itrRestartFrom = 100;    // iteration the run restarted from
    params.timeRestartFrom = 0.5;   // seconds already elapsed before restart

    double retVal = -1.0;
    const std::string out = pch_run_and_capture(
        "test_pch_restart.dat", complexList, 200, params, molTemplateList, moleculeList, 0, retVal);

    std::cerr << "  Captured output:\n" << out;

    EXPECT_DOUBLE_EQ(retVal, 0.0) << "print_complex_hist should always return 0.0";
    EXPECT_NE(out.find("Time (s): 0.5002"), std::string::npos)
        << "Time stamp should be (it - itrRestartFrom)*timeStep*1e-6 + timeRestartFrom";
}

// -----------------------------------------------------------------------------
// Test 6: Implicit-lipid mode. When molTemplateList[0] is an implicit lipid the
//         count for type 0 is taken from Complex::linksToSurface, and an extra
//         line reporting the free implicit lipid pool is prepended.
// -----------------------------------------------------------------------------
void test_pch_implicit_lipid_mode()
{
    std::cerr << "\n[TEST] test_pch_implicit_lipid_mode\n"
              << "  Source file:   src/io/print_complex_hist.cpp\n"
              << "  Function:      print_complex_hist()\n"
              << "  Scenario:      molTemplateList[0] is an implicit lipid; one\n"
              << "                 complex has 3 links to the surface plus 1 protein,\n"
              << "                 another complex is a free protein (0 links).\n"
              << "  Pass criteria: a \"<nImplicitLipids><TAB>IL: 1. \" line is written,\n"
              << "                 the bound complex reports \"IL: 3. P: 1. \", and the\n"
              << "                 free protein reports \"P: 1. \" only.\n";

    // Type 0 == implicit lipid ("IL"), type 1 == a normal protein ("P").
    std::vector<MolTemplate> molTemplateList {
        pch_make_mol_template("IL", 0),
        pch_make_mol_template("P", 1)
    };
    molTemplateList[0].isImplicitLipid = true;

    // Molecule 0 is the protein bound to the surface, molecule 1 is free.
    std::vector<Molecule> moleculeList {
        pch_make_molecule(0, 1, 0),
        pch_make_molecule(1, 1, 1)
    };

    std::vector<Complex> complexList {
        pch_make_complex(0, { 0 }),
        pch_make_complex(1, { 1 })
    };
    // The first complex is anchored to the membrane by three implicit lipids.
    complexList[0].linksToSurface = 3;
    // The second complex has no surface links at all.
    complexList[1].linksToSurface = 0;

    Parameters params;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    const int nImplicitLipids = 42;

    double retVal = -1.0;
    const std::string out = pch_run_and_capture("test_pch_implicit.dat", complexList, 0, params,
        molTemplateList, moleculeList, nImplicitLipids, retVal);

    std::cerr << "  Captured output:\n" << out;

    EXPECT_DOUBLE_EQ(retVal, 0.0) << "print_complex_hist should always return 0.0";

    // The pool of free implicit lipids is written on its own line.
    EXPECT_NE(out.find("42\tIL: 1."), std::string::npos)
        << "Expected the free implicit lipid pool line \"42\\tIL: 1.\"";

    // The bound complex reports linksToSurface as the implicit lipid count.
    EXPECT_NE(out.find("1\tIL: 3. P: 1. "), std::string::npos)
        << "The bound complex should report IL: 3 (from linksToSurface) and P: 1";

    // The free protein has zero links, so no IL entry appears for it.
    EXPECT_NE(out.find("1\tP: 1. "), std::string::npos)
        << "The free protein complex should report only \"P: 1. \"";
    EXPECT_EQ(out.find("IL: 0."), std::string::npos)
        << "A complex with zero surface links must not list IL at all";

    // header + implicit lipid line + two histogram lines == four newlines.
    EXPECT_EQ(pch_count_occurrences(out, "\n"), 4)
        << "Expected header + implicit-lipid line + two histogram lines";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_* helper is invoked from its own TEST so
// that a failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(PrintComplexHist, SingleComplexComposition) { test_pch_single_complex_composition(); }
TEST(PrintComplexHist, GroupsIdenticalComplexes) { test_pch_groups_identical_complexes(); }
TEST(PrintComplexHist, SkipsEmptyComplexes) { test_pch_skips_empty_complexes(); }
TEST(PrintComplexHist, NoCountableComplexes) { test_pch_no_countable_complexes(); }
TEST(PrintComplexHist, TimeHeaderWithRestart) { test_pch_time_header_with_restart(); }
TEST(PrintComplexHist, ImplicitLipidMode) { test_pch_implicit_lipid_mode(); }