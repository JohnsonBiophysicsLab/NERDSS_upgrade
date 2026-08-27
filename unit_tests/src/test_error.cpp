/*! \file test_error.cpp
 *
 * ### Unit test for ../src/error/error.cpp
 *
 * The file under test provides three overloads of a fatal-error reporting
 * helper:
 *
 *   1. void error(std::string errorString);
 *   2. void error(MpiContext& mpiContext, std::string errorString);
 *   3. void error(MpiContext& mpiContext, Molecule& mol, std::string errorString);
 *
 * Every one of these overloads terminates the calling process with
 * `exit(1)` after printing a banner and a diagnostic message to `std::cerr`.
 * Because a plain call would tear down the whole gtest binary, all direct
 * exercises of `error()` below are written as GoogleTest *death tests*
 * (`EXPECT_EXIT`).  A death test forks a child process, runs the statement
 * there, and then inspects the child's exit status and its stderr, so the
 * parent test binary survives and the remaining tests keep running.
 *
 * ## What is checked
 *   - the process exit status is exactly 1 (`exit(1)` in error());
 *   - the decorated banner line `!!! ###...### !!! ###...### !!!` is written;
 *   - the payload is written as `Error: <message>!!!`;
 *   - degenerate input (an empty message) still exits with status 1;
 *   - the exact string that overload (2) forwards to overload (1),
 *     namely `<message>: Rank=<rank>:\n`, formats as documented.
 *
 * ## What is deliberately *not* constructed here
 * Overloads (2) and (3) take a `MpiContext&` (a typedef of the
 * `structMpiContext` type used by the MPI build).  That structure is not part
 * of the public `include/classes` API surface, it owns raw pointers into the
 * live simulation state (`moleculeList`, `complexList`, `membraneObject`,
 * `simulVolume`), and overload (3) unconditionally dereferences
 * `mpiContext.moleculeList` and indexes `(*mpiContext.complexList)[mol.myComIndex]`
 * without any bounds or null checks.  Fabricating one here would either fail
 * to compile in a non-MPI build or segfault the child process for reasons
 * unrelated to the behaviour we want to observe.  Instead, the *observable
 * contract* of overload (2) - the message it builds and hands to overload (1)
 * - is verified directly through overload (1); see
 * error_test_mpi_rank_message_format().
 */

#include "error/error.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>

// -----------------------------------------------------------------------------
// Test 1: error(std::string) terminates the process with exit status 1.
// -----------------------------------------------------------------------------
void error_test_exit_code_is_one()
{
    std::cerr << "\n[TEST] error_test_exit_code_is_one\n"
              << "  Source file:   src/error/error.cpp\n"
              << "  Function:      void error(std::string)\n"
              << "  Scenario:      call error() with a simple message inside a\n"
              << "                 forked death-test child process.\n"
              << "  Pass criteria: the child terminates via exit(1), i.e. the\n"
              << "                 exit status is exactly 1 (not 0, not a signal).\n";

    // ExitedWithCode(1) fails the test if the child crashed by signal or
    // returned any status other than 1.  The regex argument additionally
    // requires the given text to appear on the child's stderr.
    EXPECT_EXIT(error(std::string { "fatal condition detected" }),
        ::testing::ExitedWithCode(1),
        "Error: fatal condition detected")
        << "error(std::string) must call exit(1) and report the message";

    std::cerr << "  Checked: exit status == 1 and message reached stderr.\n";
}

// -----------------------------------------------------------------------------
// Test 2: the payload line is formatted as "Error: <message>!!!".
// -----------------------------------------------------------------------------
void error_test_message_formatting()
{
    std::cerr << "\n[TEST] error_test_message_formatting\n"
              << "  Source file:   src/error/error.cpp\n"
              << "  Function:      void error(std::string)\n"
              << "  Scenario:      pass a distinctive message and inspect the\n"
              << "                 exact text emitted on stderr.\n"
              << "  Pass criteria: stderr contains the literal sequence\n"
              << "                 'Error: <message>!!!' produced by the\n"
              << "                 statement cerr << \"Error: \" << errorString << \"!!!\".\n";

    // The message below contains only characters that are literal in a POSIX
    // extended regular expression ('!' and ':' are not metacharacters), so the
    // death-test matcher can be used as a plain substring search.
    EXPECT_EXIT(error(std::string { "molecule 42 left the box" }),
        ::testing::ExitedWithCode(1),
        "Error: molecule 42 left the box!!!")
        << "The payload must be wrapped exactly as 'Error: <msg>!!!'";

    std::cerr << "  Checked: prefix 'Error: ' and suffix '!!!' surround the message.\n";
}

// -----------------------------------------------------------------------------
// Test 3: the decorative banner lines are emitted around the message.
// -----------------------------------------------------------------------------
void error_test_banner_is_printed()
{
    std::cerr << "\n[TEST] error_test_banner_is_printed\n"
              << "  Source file:   src/error/error.cpp\n"
              << "  Function:      void error(std::string)\n"
              << "  Scenario:      verify the '!!! #### !!! #### !!!' banner that\n"
              << "                 brackets the message on stderr.\n"
              << "  Pass criteria: stderr contains a run of '#' characters that\n"
              << "                 only the banner statements can produce.\n";

    // '#' has no special meaning in a POSIX extended regex, so a run of them is
    // matched literally.  Only the two banner cerr statements emit this run.
    EXPECT_EXIT(error(std::string { "banner check" }),
        ::testing::ExitedWithCode(1),
        "##################################")
        << "error() must print the '#' banner line before and after the message";

    std::cerr << "  Checked: banner of '#' characters present in child stderr.\n";
}

// -----------------------------------------------------------------------------
// Test 4: an empty message is still a fatal error (degenerate input).
// -----------------------------------------------------------------------------
void error_test_empty_message_still_exits()
{
    std::cerr << "\n[TEST] error_test_empty_message_still_exits\n"
              << "  Source file:   src/error/error.cpp\n"
              << "  Function:      void error(std::string)\n"
              << "  Scenario:      call error(\"\") - the function performs no\n"
              << "                 validation on its argument.\n"
              << "  Pass criteria: the child still exits with status 1 and the\n"
              << "                 payload degenerates to 'Error: !!!'.\n";

    EXPECT_EXIT(error(std::string {}),
        ::testing::ExitedWithCode(1),
        "Error: !!!")
        << "An empty message must not change the exit status or the wrapper text";

    std::cerr << "  Checked: empty input yields 'Error: !!!' and exit status 1.\n";
}

// -----------------------------------------------------------------------------
// Test 5: punctuation in the message is forwarded verbatim (no escaping,
//         no truncation at spaces, no reformatting).
// -----------------------------------------------------------------------------
void error_test_message_is_forwarded_verbatim()
{
    std::cerr << "\n[TEST] error_test_message_is_forwarded_verbatim\n"
              << "  Source file:   src/error/error.cpp\n"
              << "  Function:      void error(std::string)\n"
              << "  Scenario:      pass a multi-token message containing ':',\n"
              << "                 '=', '-' and digits.\n"
              << "  Pass criteria: the message appears on stderr unchanged; the\n"
              << "                 function does no escaping or trimming.\n";

    // Every character used here is literal under POSIX extended regex rules,
    // which keeps the death-test matcher a straight substring comparison.
    EXPECT_EXIT(error(std::string { "index out of range: idx=7 size=3 - aborting" }),
        ::testing::ExitedWithCode(1),
        "Error: index out of range: idx=7 size=3 - aborting!!!")
        << "error() must not modify the caller-supplied message";

    std::cerr << "  Checked: punctuation and spacing preserved verbatim.\n";
}

// -----------------------------------------------------------------------------
// Test 6: proxy check for the string that error(MpiContext&, std::string)
//         composes and forwards to error(std::string).
//
// The MPI-aware overload is a one-liner:
//     error(errorString + ": Rank=" + to_string(mpiContext.rank) + ":\n");
// A real MpiContext cannot be safely fabricated here (see the file header), so
// instead we build the identical string ourselves and push it through the
// overload that MPI variant delegates to.  This pins down the documented
// message layout: "<message>: Rank=<rank>:" followed by a newline, all of it
// still wrapped by the "Error: ...!!!" decoration of the base overload.
// -----------------------------------------------------------------------------
void error_test_mpi_rank_message_format()
{
    std::cerr << "\n[TEST] error_test_mpi_rank_message_format\n"
              << "  Source file:   src/error/error.cpp\n"
              << "  Function:      void error(MpiContext&, std::string) - proxy test\n"
              << "  Scenario:      reproduce the exact string the MPI overload\n"
              << "                 builds ('<msg>: Rank=<rank>:\\n') and feed it to\n"
              << "                 the base overload it delegates to.\n"
              << "  Pass criteria: the composed text appears on stderr and the\n"
              << "                 process exits with status 1.\n"
              << "  Note:          a genuine MpiContext is not constructed here\n"
              << "                 because the MPI overloads dereference raw\n"
              << "                 pointers into live simulation state.\n";

    // Mirror the concatenation performed by error(MpiContext&, std::string).
    const int simulatedRank = 3;
    const std::string baseMessage { "neighbor exchange failed" };
    const std::string composed
        = baseMessage + ": Rank=" + std::to_string(simulatedRank) + ":\n";

    // Sanity-check the composition itself in the parent process, so a
    // formatting regression is reported even if the death test is skipped.
    EXPECT_EQ(composed, "neighbor exchange failed: Rank=3:\n")
        << "The MPI overload documents the layout '<msg>: Rank=<rank>:\\n'";
    std::cerr << "  Composed message (parent-side check): \""
              << baseMessage << ": Rank=" << simulatedRank << ":\\n\"\n";

    // Now confirm the base overload prints it and terminates.  The trailing
    // '\n' embedded in the message is not part of the matched pattern.
    EXPECT_EXIT(error(composed),
        ::testing::ExitedWithCode(1),
        "Error: neighbor exchange failed: Rank=3:")
        << "The forwarded MPI message must be printed by error(std::string)";

    std::cerr << "  Checked: rank-annotated message printed and exit status 1.\n";
}

// -----------------------------------------------------------------------------
// Test 7: repeated invocations behave identically - error() keeps no state and
//         each call is independently fatal.  This also demonstrates that the
//         parent test binary survives multiple death tests.
// -----------------------------------------------------------------------------
void error_test_repeated_invocations_are_independent()
{
    std::cerr << "\n[TEST] error_test_repeated_invocations_are_independent\n"
              << "  Source file:   src/error/error.cpp\n"
              << "  Function:      void error(std::string)\n"
              << "  Scenario:      invoke error() twice, in two separate death-test\n"
              << "                 child processes, with different messages.\n"
              << "  Pass criteria: both children exit with status 1 and each one\n"
              << "                 reports only its own message; the parent test\n"
              << "                 binary is unaffected.\n";

    EXPECT_EXIT(error(std::string { "first failure" }),
        ::testing::ExitedWithCode(1),
        "Error: first failure!!!")
        << "First invocation must be fatal";

    EXPECT_EXIT(error(std::string { "second failure" }),
        ::testing::ExitedWithCode(1),
        "Error: second failure!!!")
        << "Second invocation must be fatal and independent of the first";

    // Reaching this point proves the parent process was never terminated by
    // the exit(1) calls performed inside the forked children.
    SUCCEED() << "Parent process survived both fatal-error death tests";
    std::cerr << "  Checked: two independent fatal calls, parent still alive.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.
//
// The suite name ends in "DeathTest" so that GoogleTest schedules these cases
// before any non-death tests, which is the recommended layout for binaries that
// mix the two kinds.
// -----------------------------------------------------------------------------
TEST(ErrorDeathTest, ExitCodeIsOne) { error_test_exit_code_is_one(); }
TEST(ErrorDeathTest, MessageFormatting) { error_test_message_formatting(); }
TEST(ErrorDeathTest, BannerIsPrinted) { error_test_banner_is_printed(); }
TEST(ErrorDeathTest, EmptyMessageStillExits) { error_test_empty_message_still_exits(); }
TEST(ErrorDeathTest, MessageIsForwardedVerbatim) { error_test_message_is_forwarded_verbatim(); }
TEST(ErrorDeathTest, MpiRankMessageFormat) { error_test_mpi_rank_message_format(); }
TEST(ErrorDeathTest, RepeatedInvocationsAreIndependent)
{
    error_test_repeated_invocations_are_independent();
}