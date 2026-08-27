/*! \file test_remove_comment.cpp
 *
 * ### Unit test for ../src/parser/remove_comment.cpp
 *
 * The only active function in that translation unit is:
 *
 *     void remove_comment(std::string& line);
 *
 * Behaviour, as implemented:
 *   - Find the FIRST '#' character in the string.
 *   - If one exists, erase from that position through the end of the string.
 *   - If none exists, leave the string completely untouched.
 *   - The modification is done in place; nothing is returned.
 *
 * Note that the implementation has no notion of quoting or escaping -- a '#'
 * inside quotes, or a '#' that is the first character of the line, is treated
 * exactly the same way as any other '#'. It also does NOT trim whitespace that
 * precedes the comment marker. The tests below pin down each of these details.
 *
 * (The companion function create_tmp_line() is commented out in the source
 * file, so there is nothing to link against and it is deliberately not tested
 * here.)
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>

// -----------------------------------------------------------------------------
// Test 1: A line with no '#' at all must be returned completely unmodified.
// -----------------------------------------------------------------------------
void rc_test_no_comment_marker_leaves_line_untouched()
{
    std::cerr << "\n[TEST] rc_test_no_comment_marker_leaves_line_untouched\n"
              << "  Source file:   src/parser/remove_comment.cpp\n"
              << "  Function:      remove_comment(std::string&)\n"
              << "  Scenario:      input string contains no '#' character.\n"
              << "  Pass criteria: the string is byte-for-byte identical afterwards.\n";

    // A typical parameter-file line with no trailing comment.
    std::string line { "nItr = 1000000" };
    const std::string expected { line }; // keep a pristine copy for comparison

    std::cerr << "  Input : \"" << line << "\"\n";
    remove_comment(line);
    std::cerr << "  Output: \"" << line << "\"\n";

    EXPECT_EQ(line, expected)
        << "remove_comment() must not alter a line that contains no '#'";

    // An empty string is the degenerate no-'#' case; it must survive as well.
    std::string emptyLine {};
    remove_comment(emptyLine);
    std::cerr << "  Empty-string case produced a string of length "
              << emptyLine.size() << '\n';
    EXPECT_TRUE(emptyLine.empty())
        << "An empty string must remain empty after remove_comment()";
}

// -----------------------------------------------------------------------------
// Test 2: A trailing comment is stripped, but everything before the '#' -- and
//         crucially the whitespace immediately before it -- is preserved.
// -----------------------------------------------------------------------------
void rc_test_trailing_comment_is_stripped()
{
    std::cerr << "\n[TEST] rc_test_trailing_comment_is_stripped\n"
              << "  Source file:   src/parser/remove_comment.cpp\n"
              << "  Function:      remove_comment(std::string&)\n"
              << "  Scenario:      'key = value  # explanation'\n"
              << "  Pass criteria: only the '#' and everything after it is erased;\n"
              << "                 preceding whitespace is intentionally kept because\n"
              << "                 the implementation does no trimming.\n";

    std::string line { "timeStep = 0.1  # microseconds" };

    std::cerr << "  Input : \"" << line << "\"\n";
    remove_comment(line);
    std::cerr << "  Output: \"" << line << "\"\n";

    // The two spaces before '#' must still be there -- remove_comment does not trim.
    EXPECT_EQ(line, std::string { "timeStep = 0.1  " })
        << "Text before '#' (including trailing spaces) must be preserved verbatim";

    // Sanity check: no '#' can remain anywhere in the result.
    EXPECT_EQ(line.find('#'), std::string::npos)
        << "No comment marker may survive in the output";
}

// -----------------------------------------------------------------------------
// Test 3: A line that begins with '#' becomes the empty string, and a line that
//         is nothing but "#" likewise becomes empty.
// -----------------------------------------------------------------------------
void rc_test_whole_line_comment_becomes_empty()
{
    std::cerr << "\n[TEST] rc_test_whole_line_comment_becomes_empty\n"
              << "  Source file:   src/parser/remove_comment.cpp\n"
              << "  Function:      remove_comment(std::string&)\n"
              << "  Scenario:      the '#' is the very first character.\n"
              << "  Pass criteria: the resulting string has length 0.\n";

    // Full-line comment such as one would find at the top of an input block.
    std::string fullComment { "# this whole line is a comment" };
    std::cerr << "  Input : \"" << fullComment << "\"\n";
    remove_comment(fullComment);
    std::cerr << "  Output: \"" << fullComment << "\" (length "
              << fullComment.size() << ")\n";

    EXPECT_TRUE(fullComment.empty())
        << "A line starting with '#' must be reduced to the empty string";
    EXPECT_EQ(fullComment.size(), 0u)
        << "Length of a stripped full-line comment must be exactly 0";

    // Degenerate case: the line is a lone '#'.
    std::string loneHash { "#" };
    std::cerr << "  Input : \"" << loneHash << "\"\n";
    remove_comment(loneHash);
    std::cerr << "  Output: \"" << loneHash << "\" (length "
              << loneHash.size() << ")\n";

    EXPECT_TRUE(loneHash.empty())
        << "A string consisting of a single '#' must become empty";
}

// -----------------------------------------------------------------------------
// Test 4: When several '#' characters are present, the erase must begin at the
//         FIRST one -- everything from there onward disappears in one go.
// -----------------------------------------------------------------------------
void rc_test_only_first_hash_matters()
{
    std::cerr << "\n[TEST] rc_test_only_first_hash_matters\n"
              << "  Source file:   src/parser/remove_comment.cpp\n"
              << "  Function:      remove_comment(std::string&)\n"
              << "  Scenario:      multiple '#' characters on one line.\n"
              << "  Pass criteria: text is truncated at the first '#'; the later\n"
              << "                 '#' characters vanish along with it.\n";

    std::string line { "value = 5 # first comment # second comment" };

    std::cerr << "  Input : \"" << line << "\"\n";
    remove_comment(line);
    std::cerr << "  Output: \"" << line << "\"\n";

    EXPECT_EQ(line, std::string { "value = 5 " })
        << "Truncation must occur at the first '#', not the last";
    EXPECT_EQ(line.find('#'), std::string::npos)
        << "Subsequent '#' characters must also be gone (they were after the first)";
}

// -----------------------------------------------------------------------------
// Test 5: The function is idempotent -- running it a second time on an already
//         stripped line changes nothing further.
// -----------------------------------------------------------------------------
void rc_test_repeated_calls_are_idempotent()
{
    std::cerr << "\n[TEST] rc_test_repeated_calls_are_idempotent\n"
              << "  Source file:   src/parser/remove_comment.cpp\n"
              << "  Function:      remove_comment(std::string&)\n"
              << "  Scenario:      call remove_comment() twice on the same string.\n"
              << "  Pass criteria: the second call leaves the string unchanged,\n"
              << "                 because no '#' remains after the first call.\n";

    std::string line { "sigma = 1.0 # binding radius in nm" };

    remove_comment(line);
    const std::string afterFirstCall { line };
    std::cerr << "  After first  call: \"" << afterFirstCall << "\"\n";

    remove_comment(line);
    std::cerr << "  After second call: \"" << line << "\"\n";

    EXPECT_EQ(line, afterFirstCall)
        << "remove_comment() must be idempotent once the comment is gone";
    EXPECT_EQ(line, std::string { "sigma = 1.0 " })
        << "The stripped content must still be the pre-'#' portion of the line";
}

// -----------------------------------------------------------------------------
// Test 6: The implementation is purely textual -- it has no concept of quoting,
//         escaping, or of '#' being embedded inside a word. Document that here
//         so the behaviour is locked in rather than accidentally "fixed".
// -----------------------------------------------------------------------------
void rc_test_hash_is_never_escaped_or_quoted()
{
    std::cerr << "\n[TEST] rc_test_hash_is_never_escaped_or_quoted\n"
              << "  Source file:   src/parser/remove_comment.cpp\n"
              << "  Function:      remove_comment(std::string&)\n"
              << "  Scenario:      '#' appears inside quotes, mid-word, and after a\n"
              << "                 backslash.\n"
              << "  Pass criteria: every one of these is treated as a plain comment\n"
              << "                 marker -- there is no escaping logic in the code.\n";

    // A '#' inside double quotes is still a comment marker to this function.
    std::string quoted { "name = \"run#1\"" };
    std::cerr << "  Input : \"" << quoted << "\"\n";
    remove_comment(quoted);
    std::cerr << "  Output: \"" << quoted << "\"\n";
    EXPECT_EQ(quoted, std::string { "name = \"run" })
        << "A quoted '#' is NOT protected; truncation happens at it";

    // A '#' glued to surrounding characters with no whitespace.
    std::string midWord { "abc#def" };
    std::cerr << "  Input : \"" << midWord << "\"\n";
    remove_comment(midWord);
    std::cerr << "  Output: \"" << midWord << "\"\n";
    EXPECT_EQ(midWord, std::string { "abc" })
        << "A '#' in the middle of a word truncates the word";

    // A backslash before '#' provides no escape either; the backslash itself
    // sits before the '#', so it survives.
    std::string escaped { "path\\#tag" };
    std::cerr << "  Input : \"" << escaped << "\"\n";
    remove_comment(escaped);
    std::cerr << "  Output: \"" << escaped << "\"\n";
    EXPECT_EQ(escaped, std::string { "path\\" })
        << "A preceding backslash does not escape '#'; it is kept, the rest is cut";
}

// -----------------------------------------------------------------------------
// Test 7: Realistic parser-style inputs, mirroring how remove_comment() is used
//         while reading NERDSS parameter / molecule files.
// -----------------------------------------------------------------------------
void rc_test_realistic_parameter_lines()
{
    std::cerr << "\n[TEST] rc_test_realistic_parameter_lines\n"
              << "  Source file:   src/parser/remove_comment.cpp\n"
              << "  Function:      remove_comment(std::string&)\n"
              << "  Scenario:      several representative input-file lines.\n"
              << "  Pass criteria: each line is truncated exactly at its first '#'\n"
              << "                 and untouched lines stay identical.\n";

    // Each entry is {input, expected output after remove_comment}.
    struct Case {
        const char* input;
        const char* expected;
    };

    const Case cases[] = {
        { "start parameters", "start parameters" },              // no comment
        { "nItr = 10 # ten steps", "nItr = 10 " },                // trailing comment
        { "#trajWrite = 100", "" },                              // fully commented out
        { "D = [1.0,1.0,1.0]#diffusion", "D = [1.0,1.0,1.0]" },  // no space before '#'
        { "   ", "   " },                                        // whitespace only, no '#'
        { "  # indented comment", "  " },                        // leading spaces kept
    };

    for (const auto& oneCase : cases) {
        std::string line { oneCase.input };
        std::cerr << "  Input : \"" << line << "\"";
        remove_comment(line);
        std::cerr << "  ->  Output: \"" << line << "\"\n";

        EXPECT_EQ(line, std::string { oneCase.expected })
            << "remove_comment(\"" << oneCase.input << "\") produced an unexpected result";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named helper is invoked from its own TEST so that a
// failure in one scenario does not stop the remaining scenarios from running.
// -----------------------------------------------------------------------------
TEST(RemoveComment, NoCommentMarkerLeavesLineUntouched)
{
    rc_test_no_comment_marker_leaves_line_untouched();
}

TEST(RemoveComment, TrailingCommentIsStripped)
{
    rc_test_trailing_comment_is_stripped();
}

TEST(RemoveComment, WholeLineCommentBecomesEmpty)
{
    rc_test_whole_line_comment_becomes_empty();
}

TEST(RemoveComment, OnlyFirstHashMatters)
{
    rc_test_only_first_hash_matters();
}

TEST(RemoveComment, RepeatedCallsAreIdempotent)
{
    rc_test_repeated_calls_are_idempotent();
}

TEST(RemoveComment, HashIsNeverEscapedOrQuoted)
{
    rc_test_hash_is_never_escaped_or_quoted();
}

TEST(RemoveComment, RealisticParameterLines)
{
    rc_test_realistic_parameter_lines();
}