/*! \file test_read_boolean.cpp
 *
 * ### Unit test for ../src/parser/read_boolean.cpp
 *
 * Function under test:
 *
 *     bool read_boolean(std::string fileLine)
 *
 * Behaviour of the implementation, step by step:
 *   1. The whole line is lower-cased with std::transform(..., ::tolower).
 *   2. remove_comment() strips any trailing comment ('#' and everything after).
 *   3. Every whitespace character (as reported by std::isspace) is erased,
 *      no matter where it appears in the string - not just leading/trailing.
 *   4. "0" or "false"  -> returns false
 *      "1" or "true"   -> returns true
 *      anything else   -> prints "FATAL ERROR: Cannot read boolean." and
 *                         calls exit(1).
 *
 * IMPORTANT: the failure branch calls exit(1), which would terminate the whole
 * gtest binary (not just one test). Therefore this file deliberately never
 * feeds read_boolean() an unrecognised token; that behaviour is documented in
 * test_rb_invalid_input_is_documented_only() instead of being exercised.
 *
 * Verbose progress information is written to stderr so a reader of the test
 * log can see exactly which source file / function is under test and what each
 * assertion is checking.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>

// -----------------------------------------------------------------------------
// Test 1: the two accepted numeric literals, "0" and "1".
// -----------------------------------------------------------------------------
void test_rb_numeric_literals()
{
    std::cerr << "\n[TEST] test_rb_numeric_literals\n"
              << "  Source file:   src/parser/read_boolean.cpp\n"
              << "  Function:      read_boolean(std::string)\n"
              << "  Scenario:      the plain numeric tokens \"0\" and \"1\".\n"
              << "  Pass criteria: \"0\" maps to false and \"1\" maps to true.\n";

    // "0" is one of the two literals explicitly mapped to false.
    std::cerr << "  -> read_boolean(\"0\") should be false\n";
    EXPECT_FALSE(read_boolean("0")) << "The token \"0\" must be parsed as false";

    // "1" is one of the two literals explicitly mapped to true.
    std::cerr << "  -> read_boolean(\"1\") should be true\n";
    EXPECT_TRUE(read_boolean("1")) << "The token \"1\" must be parsed as true";
}

// -----------------------------------------------------------------------------
// Test 2: the two accepted word literals, "false" and "true" (already lower
//         case, so the tolower pass is a no-op here).
// -----------------------------------------------------------------------------
void test_rb_word_literals()
{
    std::cerr << "\n[TEST] test_rb_word_literals\n"
              << "  Source file:   src/parser/read_boolean.cpp\n"
              << "  Function:      read_boolean(std::string)\n"
              << "  Scenario:      the lower-case words \"false\" and \"true\".\n"
              << "  Pass criteria: they map to false and true respectively.\n";

    std::cerr << "  -> read_boolean(\"false\") should be false\n";
    EXPECT_FALSE(read_boolean("false")) << "The word \"false\" must be parsed as false";

    std::cerr << "  -> read_boolean(\"true\") should be true\n";
    EXPECT_TRUE(read_boolean("true")) << "The word \"true\" must be parsed as true";
}

// -----------------------------------------------------------------------------
// Test 3: case insensitivity. The function lower-cases the input before the
//         comparison, so any mix of upper/lower case must still be accepted.
// -----------------------------------------------------------------------------
void test_rb_case_insensitivity()
{
    std::cerr << "\n[TEST] test_rb_case_insensitivity\n"
              << "  Source file:   src/parser/read_boolean.cpp\n"
              << "  Function:      read_boolean(std::string)\n"
              << "  Scenario:      mixed / upper case spellings of the words.\n"
              << "  Pass criteria: std::transform(::tolower) makes every casing\n"
              << "                 variant compare equal to the lower-case form.\n";

    // All-caps variants.
    std::cerr << "  -> read_boolean(\"TRUE\") should be true\n";
    EXPECT_TRUE(read_boolean("TRUE")) << "\"TRUE\" must be lower-cased to \"true\"";

    std::cerr << "  -> read_boolean(\"FALSE\") should be false\n";
    EXPECT_FALSE(read_boolean("FALSE")) << "\"FALSE\" must be lower-cased to \"false\"";

    // Capitalised variants.
    std::cerr << "  -> read_boolean(\"True\") should be true\n";
    EXPECT_TRUE(read_boolean("True")) << "\"True\" must be lower-cased to \"true\"";

    std::cerr << "  -> read_boolean(\"False\") should be false\n";
    EXPECT_FALSE(read_boolean("False")) << "\"False\" must be lower-cased to \"false\"";

    // Deliberately ugly mixed casing.
    std::cerr << "  -> read_boolean(\"TrUe\") should be true\n";
    EXPECT_TRUE(read_boolean("TrUe")) << "\"TrUe\" must be lower-cased to \"true\"";

    std::cerr << "  -> read_boolean(\"fAlSe\") should be false\n";
    EXPECT_FALSE(read_boolean("fAlSe")) << "\"fAlSe\" must be lower-cased to \"false\"";
}

// -----------------------------------------------------------------------------
// Test 4: whitespace handling. The erase(remove_if(isspace)) pass deletes every
//         whitespace character anywhere in the string, not merely the leading
//         and trailing ones.
// -----------------------------------------------------------------------------
void test_rb_whitespace_is_stripped()
{
    std::cerr << "\n[TEST] test_rb_whitespace_is_stripped\n"
              << "  Source file:   src/parser/read_boolean.cpp\n"
              << "  Function:      read_boolean(std::string)\n"
              << "  Scenario:      tokens padded with spaces / tabs / newlines,\n"
              << "                 including whitespace *inside* the word.\n"
              << "  Pass criteria: every std::isspace character is erased before\n"
              << "                 the comparison, so the token still matches.\n";

    // Leading and trailing spaces.
    std::cerr << "  -> read_boolean(\"   true   \") should be true\n";
    EXPECT_TRUE(read_boolean("   true   ")) << "Surrounding spaces must be removed";

    std::cerr << "  -> read_boolean(\"  false \") should be false\n";
    EXPECT_FALSE(read_boolean("  false ")) << "Surrounding spaces must be removed";

    // Tabs, carriage returns and newlines are also std::isspace characters.
    std::cerr << "  -> read_boolean(\"\\t1\\r\\n\") should be true\n";
    EXPECT_TRUE(read_boolean("\t1\r\n")) << "Tabs/CR/LF must be removed around \"1\"";

    std::cerr << "  -> read_boolean(\"\\n0\\t\") should be false\n";
    EXPECT_FALSE(read_boolean("\n0\t")) << "Tabs/LF must be removed around \"0\"";

    // Whitespace in the MIDDLE of the token is removed as well, because
    // remove_if scans the entire string. This is a real (if surprising)
    // property of the implementation, so we assert it explicitly.
    std::cerr << "  -> read_boolean(\"t r u e\") should be true (interior spaces removed)\n";
    EXPECT_TRUE(read_boolean("t r u e"))
        << "Interior whitespace is erased, so \"t r u e\" collapses to \"true\"";

    std::cerr << "  -> read_boolean(\"f a l s e\") should be false (interior spaces removed)\n";
    EXPECT_FALSE(read_boolean("f a l s e"))
        << "Interior whitespace is erased, so \"f a l s e\" collapses to \"false\"";
}

// -----------------------------------------------------------------------------
// Test 5: comment handling. remove_comment() removes a trailing '#' comment
//         before the token comparison happens.
// -----------------------------------------------------------------------------
void test_rb_trailing_comment_is_removed()
{
    std::cerr << "\n[TEST] test_rb_trailing_comment_is_removed\n"
              << "  Source file:   src/parser/read_boolean.cpp\n"
              << "  Function:      read_boolean(std::string) via remove_comment()\n"
              << "  Scenario:      an input-file style line with a trailing\n"
              << "                 '#' comment after the boolean value.\n"
              << "  Pass criteria: the comment text is discarded, leaving only\n"
              << "                 the boolean token to be interpreted.\n";

    // A value followed by an explanatory comment, as it might appear in a
    // parameter/mol input file.
    std::cerr << "  -> read_boolean(\"true #this molecule is a lipid\") should be true\n";
    EXPECT_TRUE(read_boolean("true #this molecule is a lipid"))
        << "Everything from '#' onward must be stripped before comparison";

    std::cerr << "  -> read_boolean(\"0 # not a rod\") should be false\n";
    EXPECT_FALSE(read_boolean("0 # not a rod"))
        << "Everything from '#' onward must be stripped before comparison";

    // The comment is lower-cased first, but that cannot matter because it is
    // removed; verify with an upper-case comment body.
    std::cerr << "  -> read_boolean(\"FALSE  # THIS IS A COMMENT\") should be false\n";
    EXPECT_FALSE(read_boolean("FALSE  # THIS IS A COMMENT"))
        << "Comment casing is irrelevant; the comment is discarded entirely";
}

// -----------------------------------------------------------------------------
// Test 6: the failure branch is documented but intentionally NOT exercised.
// -----------------------------------------------------------------------------
void test_rb_invalid_input_is_documented_only()
{
    std::cerr << "\n[TEST] test_rb_invalid_input_is_documented_only\n"
              << "  Source file:   src/parser/read_boolean.cpp\n"
              << "  Function:      read_boolean(std::string)\n"
              << "  Scenario:      an unrecognised token such as \"yes\", \"2\"\n"
              << "                 or an empty string.\n"
              << "  Behaviour:     the implementation prints\n"
              << "                 \"FATAL ERROR: Cannot read boolean.\" and calls\n"
              << "                 exit(1), which would tear down this entire\n"
              << "                 gtest binary - so it is NOT invoked here.\n"
              << "  Pass criteria: trivially satisfied; this test only records\n"
              << "                 the documented (untestable) failure path.\n";

    // Nothing is called. The single assertion below simply keeps the test body
    // non-empty and makes the intent explicit in the gtest report.
    SUCCEED() << "read_boolean()'s exit(1) failure path is documented, not exercised";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each named test_rb_* routine is run from its own TEST so
// that a failure in one does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(ReadBoolean, NumericLiterals) { test_rb_numeric_literals(); }
TEST(ReadBoolean, WordLiterals) { test_rb_word_literals(); }
TEST(ReadBoolean, CaseInsensitivity) { test_rb_case_insensitivity(); }
TEST(ReadBoolean, WhitespaceIsStripped) { test_rb_whitespace_is_stripped(); }
TEST(ReadBoolean, TrailingCommentIsRemoved) { test_rb_trailing_comment_is_removed(); }
TEST(ReadBoolean, InvalidInputIsDocumentedOnly) { test_rb_invalid_input_is_documented_only(); }