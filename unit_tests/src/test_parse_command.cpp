/*! \file test_parse_command.cpp
 *
 * ### Unit test for src/parser/parse_command.cpp
 *
 * Tests the single free function defined in that file:
 *
 * \code
 * void parse_command(int argc, char* argv[], Parameters& params,
 *                    std::string& paramFileName, std::string& restartFileName,
 *                    std::string& addFileName, std::string& coordinateFileName,
 *                    unsigned int& seed);
 * \endcode
 *
 * The function walks the command line and sets:
 *   -f                     -> paramFileName
 *   -s / --seed            -> seed (exits the process on a non-numeric value,
 *                             so that path is deliberately NOT exercised here)
 *   --debug-force-dissoc   -> params.debugParams.forceDissoc
 *   --debug-force-assoc    -> params.debugParams.forceAssoc
 *   --print-system-info    -> params.debugParams.printSystemInfo
 *   -r / --restart         -> restartFileName (+ params.fromRestart = true);
 *                             for params.rank >= 0 the name is forced to
 *                             "restart.dat<rank>" instead of the given argument
 *   -a / --add             -> addFileName
 *   -c / --coordinate      -> coordinateFileName
 *   -v                     -> verbosity = 1
 *   -vv                    -> verbosity = 2
 *   anything else          -> ignored (no state change)
 *
 * Every assertion below states which flag is being fed in and which member the
 * function is expected to have written to.
 */

#include "classes/class_Parameters.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// The function under test lives in src/parser/parse_command.cpp. It is declared
// in include/parser/parser_functions.hpp, but that header drags in a large
// amount of unrelated machinery (and split.cpp), so the prototype is repeated
// here verbatim to keep this translation unit small and self-contained.
// -----------------------------------------------------------------------------
void parse_command(int argc, char* argv[], Parameters& params,
    std::string& paramFileName, std::string& restartFileName,
    std::string& addFileName, std::string& coordinateFileName,
    unsigned int& seed);

namespace {

/*! \brief Small helper that owns the storage backing a C-style argv array.
 *
 * parse_command() takes `char* argv[]`, so we need mutable, NUL-terminated
 * character buffers that stay alive for the duration of the call.
 */
class PcArgv {
public:
    /*! Build the argv storage from a list of command line tokens. */
    explicit PcArgv(const std::vector<std::string>& tokens)
    {
        storage_.reserve(tokens.size());
        pointers_.reserve(tokens.size() + 1);
        for (const auto& tok : tokens) {
            storage_.emplace_back(tok.begin(), tok.end());
            storage_.back().push_back('\0');
            pointers_.push_back(storage_.back().data());
        }
        // Real argv arrays are NULL terminated; mimic that so that a stray
        // read of argv[argc] would not be undefined in our harness.
        pointers_.push_back(nullptr);
    }

    int argc() const { return static_cast<int>(storage_.size()); }
    char** argv() { return pointers_.data(); }

private:
    std::vector<std::vector<char>> storage_; //!< owns each argument's characters
    std::vector<char*> pointers_;            //!< the argv array handed to the function
};

/*! \brief Bundle of all the out-parameters parse_command() writes into.
 *
 * Initialised to obviously-distinct sentinel values so we can tell whether the
 * function touched a field or left it alone.
 */
struct PcOutputs {
    Parameters params {};
    std::string paramFileName { "UNSET_PARAM" };
    std::string restartFileName { "UNSET_RESTART" };
    std::string addFileName { "UNSET_ADD" };
    std::string coordinateFileName { "UNSET_COORD" };
    unsigned int seed { 999999 };

    /*! Parameters::rank has no default member initialiser, so it must always be
     *  set before the -r / --restart branch reads it. */
    explicit PcOutputs(int rank = -1) { params.rank = rank; }
};

/*! \brief Convenience wrapper: build argv then invoke the function under test. */
void pc_run(const std::vector<std::string>& tokens, PcOutputs& out)
{
    PcArgv argvHolder(tokens);
    std::cerr << "  -> invoking parse_command with argc=" << argvHolder.argc() << " :";
    for (const auto& tok : tokens)
        std::cerr << ' ' << tok;
    std::cerr << '\n';
    parse_command(argvHolder.argc(), argvHolder.argv(), out.params,
        out.paramFileName, out.restartFileName, out.addFileName,
        out.coordinateFileName, out.seed);
}

} // namespace

// -----------------------------------------------------------------------------
// Test: program name only. Nothing should be modified.
// -----------------------------------------------------------------------------
void test_pc_no_flags()
{
    std::cerr << "\n[parse_command.cpp] test_pc_no_flags\n"
              << "  Scenario:      only argv[0] is supplied (no flags).\n"
              << "  Pass criteria: every out-parameter keeps its sentinel value\n"
              << "                 and no debug flag is switched on.\n";

    PcOutputs out; // rank = -1 (serial)
    pc_run({ "./nerdss" }, out);

    EXPECT_EQ(out.paramFileName, "UNSET_PARAM") << "paramFileName must not change without -f";
    EXPECT_EQ(out.restartFileName, "UNSET_RESTART") << "restartFileName must not change without -r";
    EXPECT_EQ(out.addFileName, "UNSET_ADD") << "addFileName must not change without -a";
    EXPECT_EQ(out.coordinateFileName, "UNSET_COORD") << "coordinateFileName must not change without -c";
    EXPECT_EQ(out.seed, 999999u) << "seed must not change without -s/--seed";

    EXPECT_FALSE(out.params.fromRestart) << "fromRestart defaults to false";
    EXPECT_FALSE(out.params.debugParams.forceAssoc) << "forceAssoc defaults to false";
    EXPECT_FALSE(out.params.debugParams.forceDissoc) << "forceDissoc defaults to false";
    EXPECT_FALSE(out.params.debugParams.printSystemInfo) << "printSystemInfo defaults to false";
    EXPECT_EQ(out.params.debugParams.verbosity, 0) << "verbosity defaults to 0";
}

// -----------------------------------------------------------------------------
// Test: -f <file> sets the parameter file name and consumes its argument.
// -----------------------------------------------------------------------------
void test_pc_param_file_flag()
{
    std::cerr << "\n[parse_command.cpp] test_pc_param_file_flag\n"
              << "  Flag tested:   -f <paramFile>\n"
              << "  Pass criteria: paramFileName equals the following token and\n"
              << "                 that token is consumed (not re-parsed as a flag).\n";

    PcOutputs out;
    pc_run({ "./nerdss", "-f", "parms.inp" }, out);

    EXPECT_EQ(out.paramFileName, "parms.inp") << "-f should copy the next token into paramFileName";
    // "parms.inp" was consumed as an argument, so it must not have been seen as
    // an unrecognised flag nor changed anything else.
    EXPECT_EQ(out.restartFileName, "UNSET_RESTART") << "-f must not touch restartFileName";
    EXPECT_FALSE(out.params.fromRestart) << "-f must not set fromRestart";
    std::cerr << "  paramFileName = \"" << out.paramFileName << "\"\n";
}

// -----------------------------------------------------------------------------
// Test: -s and --seed both parse an unsigned integer into seed.
// -----------------------------------------------------------------------------
void test_pc_seed_flags()
{
    std::cerr << "\n[parse_command.cpp] test_pc_seed_flags\n"
              << "  Flags tested:  -s <n> and --seed <n>\n"
              << "  Pass criteria: seed holds the parsed unsigned value.\n"
              << "  NOTE:          a non-numeric seed makes the function call\n"
              << "                 exit(1), which would kill the whole test\n"
              << "                 binary, so that path is intentionally untested.\n";

    // Short form.
    PcOutputs shortForm;
    pc_run({ "./nerdss", "-s", "42" }, shortForm);
    EXPECT_EQ(shortForm.seed, 42u) << "-s 42 should store 42 in seed";

    // Long form.
    PcOutputs longForm;
    pc_run({ "./nerdss", "--seed", "1234567" }, longForm);
    EXPECT_EQ(longForm.seed, 1234567u) << "--seed 1234567 should store 1234567 in seed";

    // A later seed flag overwrites an earlier one (flags are processed in order).
    PcOutputs twice;
    pc_run({ "./nerdss", "-s", "7", "--seed", "9" }, twice);
    EXPECT_EQ(twice.seed, 9u) << "the last seed flag on the line should win";

    std::cerr << "  seeds parsed: " << shortForm.seed << ", " << longForm.seed
              << ", " << twice.seed << '\n';
}

// -----------------------------------------------------------------------------
// Test: the three boolean debug flags.
// -----------------------------------------------------------------------------
void test_pc_debug_boolean_flags()
{
    std::cerr << "\n[parse_command.cpp] test_pc_debug_boolean_flags\n"
              << "  Flags tested:  --debug-force-dissoc, --debug-force-assoc,\n"
              << "                 --print-system-info\n"
              << "  Pass criteria: each sets exactly its own Parameters::Debug member.\n";

    // Only forceDissoc.
    PcOutputs dissoc;
    pc_run({ "./nerdss", "--debug-force-dissoc" }, dissoc);
    EXPECT_TRUE(dissoc.params.debugParams.forceDissoc) << "--debug-force-dissoc sets forceDissoc";
    EXPECT_FALSE(dissoc.params.debugParams.forceAssoc) << "--debug-force-dissoc must not set forceAssoc";
    EXPECT_FALSE(dissoc.params.debugParams.printSystemInfo)
        << "--debug-force-dissoc must not set printSystemInfo";

    // Only forceAssoc.
    PcOutputs assoc;
    pc_run({ "./nerdss", "--debug-force-assoc" }, assoc);
    EXPECT_TRUE(assoc.params.debugParams.forceAssoc) << "--debug-force-assoc sets forceAssoc";
    EXPECT_FALSE(assoc.params.debugParams.forceDissoc) << "--debug-force-assoc must not set forceDissoc";

    // Only printSystemInfo.
    PcOutputs info;
    pc_run({ "./nerdss", "--print-system-info" }, info);
    EXPECT_TRUE(info.params.debugParams.printSystemInfo) << "--print-system-info sets printSystemInfo";
    EXPECT_FALSE(info.params.debugParams.forceAssoc) << "--print-system-info must not set forceAssoc";

    // All three together.
    PcOutputs all;
    pc_run({ "./nerdss", "--debug-force-dissoc", "--debug-force-assoc", "--print-system-info" }, all);
    EXPECT_TRUE(all.params.debugParams.forceDissoc) << "combined: forceDissoc set";
    EXPECT_TRUE(all.params.debugParams.forceAssoc) << "combined: forceAssoc set";
    EXPECT_TRUE(all.params.debugParams.printSystemInfo) << "combined: printSystemInfo set";
}

// -----------------------------------------------------------------------------
// Test: -v and -vv set verbosity levels 1 and 2.
// -----------------------------------------------------------------------------
void test_pc_verbosity_flags()
{
    std::cerr << "\n[parse_command.cpp] test_pc_verbosity_flags\n"
              << "  Flags tested:  -v and -vv\n"
              << "  Pass criteria: verbosity becomes 1 and 2 respectively; when\n"
              << "                 both appear the later flag wins.\n";

    PcOutputs v1;
    pc_run({ "./nerdss", "-v" }, v1);
    EXPECT_EQ(v1.params.debugParams.verbosity, 1) << "-v should set verbosity to 1";

    PcOutputs v2;
    pc_run({ "./nerdss", "-vv" }, v2);
    EXPECT_EQ(v2.params.debugParams.verbosity, 2) << "-vv should set verbosity to 2";

    // Flags are applied strictly in order, so the trailing -v resets it to 1.
    PcOutputs both;
    pc_run({ "./nerdss", "-vv", "-v" }, both);
    EXPECT_EQ(both.params.debugParams.verbosity, 1)
        << "the last verbosity flag encountered should determine the value";

    std::cerr << "  verbosity results: " << v1.params.debugParams.verbosity << ", "
              << v2.params.debugParams.verbosity << ", "
              << both.params.debugParams.verbosity << '\n';
}

// -----------------------------------------------------------------------------
// Test: restart flag with a negative rank (the serial case) copies the argument.
// -----------------------------------------------------------------------------
void test_pc_restart_serial()
{
    std::cerr << "\n[parse_command.cpp] test_pc_restart_serial\n"
              << "  Flags tested:  -r <file> and --restart <file> with params.rank < 0\n"
              << "  Pass criteria: restartFileName equals the supplied token and\n"
              << "                 params.fromRestart becomes true.\n";

    // Short form, serial job (rank = -1).
    PcOutputs shortForm(-1);
    pc_run({ "./nerdss", "-r", "my_restart.dat" }, shortForm);
    EXPECT_EQ(shortForm.restartFileName, "my_restart.dat")
        << "serial -r should copy the given restart file name";
    EXPECT_TRUE(shortForm.params.fromRestart) << "-r must set fromRestart to true";

    // Long form, serial job.
    PcOutputs longForm(-5);
    pc_run({ "./nerdss", "--restart", "other_restart.dat" }, longForm);
    EXPECT_EQ(longForm.restartFileName, "other_restart.dat")
        << "serial --restart should copy the given restart file name";
    EXPECT_TRUE(longForm.params.fromRestart) << "--restart must set fromRestart to true";

    std::cerr << "  restartFileName (serial) = \"" << shortForm.restartFileName << "\"\n";
}

// -----------------------------------------------------------------------------
// Test: restart flag with a non-negative rank (the parallel case) ignores the
//       supplied name and synthesises "restart.dat<rank>".
// -----------------------------------------------------------------------------
void test_pc_restart_parallel()
{
    std::cerr << "\n[parse_command.cpp] test_pc_restart_parallel\n"
              << "  Flags tested:  -r <file> with params.rank >= 0\n"
              << "  Pass criteria: the supplied name is discarded and the file name\n"
              << "                 becomes \"restart.dat\" + rank; fromRestart is true.\n";

    // rank 0 -> "restart.dat0"
    PcOutputs rank0(0);
    pc_run({ "./nerdss", "-r", "ignored_name.dat" }, rank0);
    EXPECT_EQ(rank0.restartFileName, "restart.dat0")
        << "rank 0 should produce restart.dat0 regardless of the argument";
    EXPECT_TRUE(rank0.params.fromRestart) << "parallel -r must set fromRestart";

    // rank 3 -> "restart.dat3"
    PcOutputs rank3(3);
    pc_run({ "./nerdss", "--restart", "also_ignored.dat" }, rank3);
    EXPECT_EQ(rank3.restartFileName, "restart.dat3")
        << "rank 3 should produce restart.dat3 regardless of the argument";
    EXPECT_TRUE(rank3.params.fromRestart) << "parallel --restart must set fromRestart";

    std::cerr << "  restartFileName (rank 0) = \"" << rank0.restartFileName << "\"\n"
              << "  restartFileName (rank 3) = \"" << rank3.restartFileName << "\"\n";
}

// -----------------------------------------------------------------------------
// Test: -a/--add and -c/--coordinate capture their following argument.
// -----------------------------------------------------------------------------
void test_pc_add_and_coordinate_flags()
{
    std::cerr << "\n[parse_command.cpp] test_pc_add_and_coordinate_flags\n"
              << "  Flags tested:  -a/--add <file>, -c/--coordinate <file>\n"
              << "  Pass criteria: addFileName / coordinateFileName pick up the\n"
              << "                 following token; neither affects the other.\n";

    // Short forms.
    PcOutputs shortForm;
    pc_run({ "./nerdss", "-a", "add.inp", "-c", "coords.dat" }, shortForm);
    EXPECT_EQ(shortForm.addFileName, "add.inp") << "-a should set addFileName";
    EXPECT_EQ(shortForm.coordinateFileName, "coords.dat") << "-c should set coordinateFileName";
    EXPECT_EQ(shortForm.paramFileName, "UNSET_PARAM") << "-a/-c must not touch paramFileName";

    // Long forms.
    PcOutputs longForm;
    pc_run({ "./nerdss", "--add", "add2.inp", "--coordinate", "coords2.dat" }, longForm);
    EXPECT_EQ(longForm.addFileName, "add2.inp") << "--add should set addFileName";
    EXPECT_EQ(longForm.coordinateFileName, "coords2.dat") << "--coordinate should set coordinateFileName";

    std::cerr << "  addFileName = \"" << shortForm.addFileName
              << "\", coordinateFileName = \"" << shortForm.coordinateFileName << "\"\n";
}

// -----------------------------------------------------------------------------
// Test: unrecognised tokens are ignored and do not disturb known flags.
// -----------------------------------------------------------------------------
void test_pc_unknown_flags_ignored()
{
    std::cerr << "\n[parse_command.cpp] test_pc_unknown_flags_ignored\n"
              << "  Scenario:      unknown flags are interleaved with known ones.\n"
              << "  Pass criteria: unknown tokens change nothing; the known -f flag\n"
              << "                 still resolves correctly.\n";

    PcOutputs out;
    pc_run({ "./nerdss", "--not-a-flag", "-f", "parms.inp", "--another-bogus-flag" }, out);

    EXPECT_EQ(out.paramFileName, "parms.inp") << "-f should still be honoured amid unknown flags";
    EXPECT_EQ(out.addFileName, "UNSET_ADD") << "unknown flags must not set addFileName";
    EXPECT_EQ(out.coordinateFileName, "UNSET_COORD") << "unknown flags must not set coordinateFileName";
    EXPECT_EQ(out.seed, 999999u) << "unknown flags must not set seed";
    EXPECT_FALSE(out.params.fromRestart) << "unknown flags must not set fromRestart";
    EXPECT_EQ(out.params.debugParams.verbosity, 0) << "unknown flags must not change verbosity";
}

// -----------------------------------------------------------------------------
// Test: a realistic full command line, exercising every branch at once.
// -----------------------------------------------------------------------------
void test_pc_full_command_line()
{
    std::cerr << "\n[parse_command.cpp] test_pc_full_command_line\n"
              << "  Scenario:      a complete command line combining every flag.\n"
              << "  Pass criteria: all out-parameters and debug members reflect the\n"
              << "                 flags, with argument-taking flags consuming their\n"
              << "                 value (so file names are never mistaken for flags).\n";

    PcOutputs out(-1); // serial, so -r keeps the given name
    pc_run({ "./nerdss",
                "-f", "parms.inp",
                "-s", "2024",
                "-r", "restart_serial.dat",
                "-a", "add.inp",
                "-c", "coords.dat",
                "--debug-force-assoc",
                "--debug-force-dissoc",
                "--print-system-info",
                "-vv" },
        out);

    EXPECT_EQ(out.paramFileName, "parms.inp") << "combined: -f";
    EXPECT_EQ(out.seed, 2024u) << "combined: -s";
    EXPECT_EQ(out.restartFileName, "restart_serial.dat") << "combined: -r (serial)";
    EXPECT_TRUE(out.params.fromRestart) << "combined: -r sets fromRestart";
    EXPECT_EQ(out.addFileName, "add.inp") << "combined: -a";
    EXPECT_EQ(out.coordinateFileName, "coords.dat") << "combined: -c";
    EXPECT_TRUE(out.params.debugParams.forceAssoc) << "combined: --debug-force-assoc";
    EXPECT_TRUE(out.params.debugParams.forceDissoc) << "combined: --debug-force-dissoc";
    EXPECT_TRUE(out.params.debugParams.printSystemInfo) << "combined: --print-system-info";
    EXPECT_EQ(out.params.debugParams.verbosity, 2) << "combined: -vv";

    std::cerr << "  parsed: param=\"" << out.paramFileName << "\", seed=" << out.seed
              << ", restart=\"" << out.restartFileName << "\", add=\"" << out.addFileName
              << "\", coord=\"" << out.coordinateFileName << "\", verbosity="
              << out.params.debugParams.verbosity << '\n';
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is invoked from its own TEST so that a
// failure in one scenario does not prevent the others from running.
// -----------------------------------------------------------------------------
TEST(ParseCommandTest, NoFlags) { test_pc_no_flags(); }
TEST(ParseCommandTest, ParamFileFlag) { test_pc_param_file_flag(); }
TEST(ParseCommandTest, SeedFlags) { test_pc_seed_flags(); }
TEST(ParseCommandTest, DebugBooleanFlags) { test_pc_debug_boolean_flags(); }
TEST(ParseCommandTest, VerbosityFlags) { test_pc_verbosity_flags(); }
TEST(ParseCommandTest, RestartSerial) { test_pc_restart_serial(); }
TEST(ParseCommandTest, RestartParallel) { test_pc_restart_parallel(); }
TEST(ParseCommandTest, AddAndCoordinateFlags) { test_pc_add_and_coordinate_flags(); }
TEST(ParseCommandTest, UnknownFlagsIgnored) { test_pc_unknown_flags_ignored(); }
TEST(ParseCommandTest, FullCommandLine) { test_pc_full_command_line(); }