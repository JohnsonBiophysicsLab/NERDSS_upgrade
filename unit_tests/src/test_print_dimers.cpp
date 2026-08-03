/*! \file test_print_dimers.cpp
 *
 * ### Unit test for src/io/print_dimers.cpp
 *
 * The file under test provides one function:
 *
 *     void print_dimers(std::vector<Complex>& complexList, std::ofstream& outfile,
 *                       int it, Parameters params,
 *                       std::vector<MolTemplate>& molTemplateList);
 *
 * The routine walks over every non-empty Complex in `complexList`, builds a
 * unique "composition index" for each distinct complex composition
 * (using MolTemplate::numEachMolType to build the multipliers), histograms how
 * many complexes share each composition, and then writes one line to `outfile`:
 *
 *     <time>\t<monomers[0]>\t<dimers[0]>\t...<monomers[N-1]>\t<dimers[N-1]>\n
 *
 * where
 *   - time      = (it - params.itrRestartFrom) * params.timeStep * 1e-6
 *                 + params.timeRestartFrom
 *   - monomers[j] = number of complexes whose representative has exactly one
 *                   member molecule and one copy of molecule type j
 *   - dimers[j]   = summed number of complexes whose representative has exactly
 *                   two member molecules and 1 or 2 copies of molecule type j
 *
 * Special case: if `complexList` is empty, the function writes "0\t0\t" once
 * per entry in molTemplateList and returns immediately (no time stamp, no
 * newline).
 *
 * The tests below drive the function with hand-built Complex objects, capture
 * the produced text in a temporary file, and then verify the parsed tokens.
 * Everything is reported verbosely to stderr so the reader can follow along.
 */

#include "classes/class_Molecule_Complex.hpp"
#include "io/io.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (all prefixed with pd_ to stay unique inside the test suite).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal Complex suitable for print_dimers().
 *
 * print_dimers() only ever looks at three fields of a Complex:
 *   - isEmpty      (empty slots are skipped entirely)
 *   - numEachMol   (composition, one entry per molecule type)
 *   - memberList   (its *size* decides monomer vs. dimer vs. larger)
 *
 * \param[in] numEachMol Copy numbers of each molecule type in the complex.
 * \param[in] numMembers How many member molecules the complex has.
 * \param[in] isEmpty    Mark this complex as a destroyed/void slot.
 */
Complex pd_make_complex(const std::vector<int>& numEachMol, int numMembers,
    bool isEmpty = false)
{
    Complex com;
    com.numEachMol = numEachMol;
    com.memberList.clear();
    for (int i = 0; i < numMembers; ++i)
        com.memberList.push_back(i); // actual indices are irrelevant here
    com.isEmpty = isEmpty;
    return com;
}

/*! \brief Build a MolTemplate list of the requested length.
 *
 * Only the number of entries matters for the "empty complexList" branch, so the
 * templates are otherwise left at their defaults (apart from a readable name).
 */
std::vector<MolTemplate> pd_make_mol_templates(int nTypes)
{
    std::vector<MolTemplate> molTemplateList;
    for (int i = 0; i < nTypes; ++i) {
        MolTemplate temp;
        temp.molTypeIndex = i;
        temp.molName = "mol" + std::to_string(i);
        molTemplateList.push_back(temp);
    }
    return molTemplateList;
}

/*! \brief Call print_dimers(), capture its output and return the raw text.
 *
 * The output is written to a small temporary file (which is deleted again) so
 * that we exercise the real std::ofstream interface the function demands.
 */
std::string pd_run_and_capture(std::vector<Complex>& complexList, int it,
    const Parameters& params, std::vector<MolTemplate>& molTemplateList,
    const std::string& tag)
{
    const std::string path = "test_print_dimers_" + tag + ".dat";

    std::ofstream outFile(path);
    // Guard: if we cannot open the file the test cannot be meaningful.
    if (!outFile.is_open()) {
        std::cerr << "    !! could not open temporary file " << path << '\n';
        return std::string {};
    }
    print_dimers(complexList, outFile, it, params, molTemplateList);
    outFile.close();

    std::ifstream inFile(path);
    std::ostringstream contents;
    contents << inFile.rdbuf();
    inFile.close();
    std::remove(path.c_str());

    return contents.str();
}

/*! \brief Split whitespace separated text into a vector of tokens. */
std::vector<std::string> pd_tokenize(const std::string& text)
{
    std::vector<std::string> tokens;
    std::istringstream iss(text);
    std::string tok;
    while (iss >> tok)
        tokens.push_back(tok);
    return tokens;
}

/*! \brief Convenience printer for the captured tokens. */
void pd_dump_tokens(const std::vector<std::string>& tokens)
{
    std::cerr << "    tokens (" << tokens.size() << "):";
    for (const auto& t : tokens)
        std::cerr << " [" << t << ']';
    std::cerr << '\n';
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: empty complexList -> "0\t0\t" per MolTemplate, nothing else.
// -----------------------------------------------------------------------------
void test_pd_empty_complex_list()
{
    std::cerr << "\n[TEST] test_pd_empty_complex_list\n"
              << "  Source file:   src/io/print_dimers.cpp\n"
              << "  Function:      print_dimers()\n"
              << "  Scenario:      complexList is empty, molTemplateList has 3 types.\n"
              << "  Pass criteria: exactly 6 zero tokens are written (\"0\\t0\\t\" per\n"
              << "                 molecule template) and no time stamp/newline.\n";

    // Two molecule types in the parameters, three templates in the list: the
    // early-return branch keys off molTemplateList.size(), not numMolTypes.
    Parameters params;
    params.numMolTypes = 2;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    MolTemplate::numEachMolType = std::vector<int> { 4, 4 };

    std::vector<Complex> complexList {}; // intentionally empty
    std::vector<MolTemplate> molTemplateList = pd_make_mol_templates(3);

    const std::string raw
        = pd_run_and_capture(complexList, 10, params, molTemplateList, "empty");
    const std::vector<std::string> tokens = pd_tokenize(raw);
    pd_dump_tokens(tokens);

    // 3 templates * 2 zeros each = 6 tokens.
    EXPECT_EQ(tokens.size(), 6u)
        << "Empty complexList should emit two zeros per molecule template";
    for (size_t i = 0; i < tokens.size(); ++i)
        EXPECT_EQ(tokens[i], "0") << "Token " << i << " should be the literal 0";

    // The early-return path must not terminate the line.
    EXPECT_EQ(raw.find('\n'), std::string::npos)
        << "Empty-list branch should not write a newline";
}

// -----------------------------------------------------------------------------
// Test 2: single monomer of a single molecule type.
// -----------------------------------------------------------------------------
void test_pd_single_monomer()
{
    std::cerr << "\n[TEST] test_pd_single_monomer\n"
              << "  Source file:   src/io/print_dimers.cpp\n"
              << "  Function:      print_dimers()\n"
              << "  Scenario:      one complex, one member molecule of type 0.\n"
              << "  Pass criteria: line is <time> 1 0 with time = it*timeStep*1e-6.\n";

    Parameters params;
    params.numMolTypes = 1;
    params.timeStep = 1.0; // microseconds
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    // Multipliers are derived from this static list; one type -> mult[0] == 1.
    MolTemplate::numEachMolType = std::vector<int> { 1 };

    std::vector<Complex> complexList;
    complexList.push_back(pd_make_complex(std::vector<int> { 1 }, 1)); // monomer

    std::vector<MolTemplate> molTemplateList = pd_make_mol_templates(1);

    const int it = 100;
    const std::string raw
        = pd_run_and_capture(complexList, it, params, molTemplateList, "monomer");
    const std::vector<std::string> tokens = pd_tokenize(raw);
    pd_dump_tokens(tokens);

    // Expect: time, monomers[0], dimers[0]
    ASSERT_EQ(tokens.size(), 3u)
        << "Expected 3 whitespace separated fields (time, monomers, dimers)";

    const double expectedTime
        = (it - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom;
    const double reportedTime = std::stod(tokens[0]);
    std::cerr << "    expected time = " << expectedTime
              << ", reported time = " << reportedTime << '\n';
    EXPECT_NEAR(reportedTime, expectedTime, 1e-12)
        << "First column must be the simulation time in seconds";

    EXPECT_EQ(tokens[1], "1") << "monomers[0] should be 1 (single monomeric complex)";
    EXPECT_EQ(tokens[2], "0") << "dimers[0] should be 0 (no dimers present)";

    // A populated line must end with a newline (function uses std::endl).
    EXPECT_FALSE(raw.empty());
    if (!raw.empty())
        EXPECT_EQ(raw.back(), '\n') << "Populated output line should be newline terminated";
}

// -----------------------------------------------------------------------------
// Test 3: mixture of monomers, a homodimer and a heterodimer.
// -----------------------------------------------------------------------------
void test_pd_monomers_and_dimers()
{
    std::cerr << "\n[TEST] test_pd_monomers_and_dimers\n"
              << "  Source file:   src/io/print_dimers.cpp\n"
              << "  Function:      print_dimers()\n"
              << "  Scenario:      2 molecule types; 2 identical A monomers, one A-A\n"
              << "                 homodimer and one A-B heterodimer.\n"
              << "  Pass criteria: monomers = {2, 0}; dimers = {2, 1}\n"
              << "                 (dimers[0] sums the homodimer and the heterodimer,\n"
              << "                  dimers[1] counts only the heterodimer).\n";

    Parameters params;
    params.numMolTypes = 2;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    // mult[0] = 1, mult[1] = numEachMolType[0] + 1 = 5, which keeps every
    // composition index used below unique (1, 2 and 6).
    MolTemplate::numEachMolType = std::vector<int> { 4, 2 };

    std::vector<Complex> complexList;
    complexList.push_back(pd_make_complex(std::vector<int> { 1, 0 }, 1)); // A monomer
    complexList.push_back(pd_make_complex(std::vector<int> { 1, 0 }, 1)); // A monomer (same type)
    complexList.push_back(pd_make_complex(std::vector<int> { 2, 0 }, 2)); // A-A homodimer
    complexList.push_back(pd_make_complex(std::vector<int> { 1, 1 }, 2)); // A-B heterodimer

    std::vector<MolTemplate> molTemplateList = pd_make_mol_templates(2);

    const std::string raw
        = pd_run_and_capture(complexList, 0, params, molTemplateList, "mixture");
    const std::vector<std::string> tokens = pd_tokenize(raw);
    pd_dump_tokens(tokens);

    // Expect: time, monomers[0], dimers[0], monomers[1], dimers[1]
    ASSERT_EQ(tokens.size(), 5u)
        << "Expected time plus (monomer,dimer) pair for each of the 2 types";

    EXPECT_EQ(tokens[1], "2") << "monomers[0]: two separate A monomers exist";
    EXPECT_EQ(tokens[2], "2")
        << "dimers[0]: homodimer (2 copies of A) + heterodimer (1 copy of A)";
    EXPECT_EQ(tokens[3], "0") << "monomers[1]: no B monomer complex present";
    EXPECT_EQ(tokens[4], "1") << "dimers[1]: only the A-B heterodimer contains B";
}

// -----------------------------------------------------------------------------
// Test 4: empty complex slots are skipped, including a leading empty slot.
// -----------------------------------------------------------------------------
void test_pd_skips_empty_slots()
{
    std::cerr << "\n[TEST] test_pd_skips_empty_slots\n"
              << "  Source file:   src/io/print_dimers.cpp\n"
              << "  Function:      print_dimers()\n"
              << "  Scenario:      complexList[0] and complexList[2] are marked\n"
              << "                 isEmpty; only two real complexes remain.\n"
              << "  Pass criteria: the empty slots contribute nothing, so the counts\n"
              << "                 match a list containing just the live complexes.\n";

    Parameters params;
    params.numMolTypes = 1;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    MolTemplate::numEachMolType = std::vector<int> { 3 };

    std::vector<Complex> complexList;
    // Leading empty slot: exercises the "while (isEmpty) i++" search.
    complexList.push_back(pd_make_complex(std::vector<int> { 1 }, 1, /*isEmpty=*/true));
    complexList.push_back(pd_make_complex(std::vector<int> { 1 }, 1)); // live monomer
    complexList.push_back(pd_make_complex(std::vector<int> { 2 }, 2, /*isEmpty=*/true));
    complexList.push_back(pd_make_complex(std::vector<int> { 2 }, 2)); // live homodimer

    std::vector<MolTemplate> molTemplateList = pd_make_mol_templates(1);

    const std::string raw
        = pd_run_and_capture(complexList, 0, params, molTemplateList, "emptyslots");
    const std::vector<std::string> tokens = pd_tokenize(raw);
    pd_dump_tokens(tokens);

    ASSERT_EQ(tokens.size(), 3u) << "Expected time, monomers[0] and dimers[0]";
    EXPECT_EQ(tokens[1], "1") << "monomers[0]: only the single live monomer counts";
    EXPECT_EQ(tokens[2], "1") << "dimers[0]: only the single live homodimer counts";
}

// -----------------------------------------------------------------------------
// Test 5: larger assemblies are neither monomers nor dimers.
// -----------------------------------------------------------------------------
void test_pd_larger_assemblies_ignored()
{
    std::cerr << "\n[TEST] test_pd_larger_assemblies_ignored\n"
              << "  Source file:   src/io/print_dimers.cpp\n"
              << "  Function:      print_dimers()\n"
              << "  Scenario:      one trimer (3 members) and one tetramer (4 members)\n"
              << "                 alongside a single monomer.\n"
              << "  Pass criteria: monomers[0] == 1 and dimers[0] == 0; assemblies\n"
              << "                 larger than two members are not reported.\n";

    Parameters params;
    params.numMolTypes = 1;
    params.timeStep = 1.0;
    params.itrRestartFrom = 0;
    params.timeRestartFrom = 0.0;

    // Allow composition indices up to 4 copies of the single molecule type.
    MolTemplate::numEachMolType = std::vector<int> { 8 };

    std::vector<Complex> complexList;
    complexList.push_back(pd_make_complex(std::vector<int> { 1 }, 1)); // monomer
    complexList.push_back(pd_make_complex(std::vector<int> { 3 }, 3)); // trimer
    complexList.push_back(pd_make_complex(std::vector<int> { 4 }, 4)); // tetramer

    std::vector<MolTemplate> molTemplateList = pd_make_mol_templates(1);

    const std::string raw
        = pd_run_and_capture(complexList, 0, params, molTemplateList, "large");
    const std::vector<std::string> tokens = pd_tokenize(raw);
    pd_dump_tokens(tokens);

    ASSERT_EQ(tokens.size(), 3u) << "Expected time, monomers[0] and dimers[0]";
    EXPECT_EQ(tokens[1], "1") << "monomers[0]: only the monomer is counted";
    EXPECT_EQ(tokens[2], "0")
        << "dimers[0]: trimers/tetramers must not be counted as dimers";
}

// -----------------------------------------------------------------------------
// Test 6: the time stamp honours the restart offsets.
// -----------------------------------------------------------------------------
void test_pd_time_column_with_restart()
{
    std::cerr << "\n[TEST] test_pd_time_column_with_restart\n"
              << "  Source file:   src/io/print_dimers.cpp\n"
              << "  Function:      print_dimers()\n"
              << "  Scenario:      itrRestartFrom = 100, timeRestartFrom = 0.5 s,\n"
              << "                 timeStep = 2 us, it = 200.\n"
              << "  Pass criteria: first column == (200-100)*2e-6 + 0.5 = 0.5002 s.\n";

    Parameters params;
    params.numMolTypes = 1;
    params.timeStep = 2.0; // microseconds
    params.itrRestartFrom = 100;
    params.timeRestartFrom = 0.5; // seconds

    MolTemplate::numEachMolType = std::vector<int> { 1 };

    std::vector<Complex> complexList;
    complexList.push_back(pd_make_complex(std::vector<int> { 1 }, 1));

    std::vector<MolTemplate> molTemplateList = pd_make_mol_templates(1);

    const int it = 200;
    const std::string raw
        = pd_run_and_capture(complexList, it, params, molTemplateList, "restart");
    const std::vector<std::string> tokens = pd_tokenize(raw);
    pd_dump_tokens(tokens);

    ASSERT_GE(tokens.size(), 1u) << "At least the time column must be present";

    const double expectedTime
        = (it - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom;
    const double reportedTime = std::stod(tokens[0]);
    std::cerr << "    expected time = " << expectedTime
              << ", reported time = " << reportedTime << '\n';

    // ostream default precision is 6 significant digits, so allow a small slack.
    EXPECT_NEAR(reportedTime, expectedTime, 1e-9)
        << "Time column must include both the iteration and restart offsets";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each scenario runs in its own TEST so a failure in one
// does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(PrintDimersTest, EmptyComplexList) { test_pd_empty_complex_list(); }
TEST(PrintDimersTest, SingleMonomer) { test_pd_single_monomer(); }
TEST(PrintDimersTest, MonomersAndDimers) { test_pd_monomers_and_dimers(); }
TEST(PrintDimersTest, SkipsEmptySlots) { test_pd_skips_empty_slots(); }
TEST(PrintDimersTest, LargerAssembliesIgnored) { test_pd_larger_assemblies_ignored(); }
TEST(PrintDimersTest, TimeColumnWithRestart) { test_pd_time_column_with_restart(); }