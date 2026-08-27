/*! \file test_parse_input_array.cpp
 *
 * ### Unit test for ../src/parser/parse_input_array.cpp
 *
 * This test exercises the single free function defined in that file:
 *
 *     std::vector<double> parse_input_array(std::string& line);
 *
 * Behaviour of the function, as read from the implementation:
 *
 *   1. It removes the *first* occurrence of each of the characters
 *      '[', ']', '(' and ')' from the input string.  (Only one of each,
 *      because `find()` is called exactly once per bracket character.)
 *   2. It then splits what remains on commas and converts every token to a
 *      double with `std::stod()`.
 *   3. If `std::stod()` throws, the token is lower-cased and compared against
 *      "m_pi"/"pi" (-> M_PI) and "nan" (-> quiet NaN).  Anything else causes a
 *      `std::string` to be thrown out of the function (the trailing
 *      `catch (const std::string&)` clause cannot catch a throw issued from
 *      inside a sibling handler, so the exception escapes to the caller).
 *   4. IMPORTANT SIDE EFFECT: `line` is passed by non-const reference and is
 *      destroyed by the parsing loop -- on return it holds only the final
 *      comma-separated token (with the brackets already stripped).
 *
 * Every test below prints what it is doing to stderr so a reader of the test
 * log can follow which function and which criterion is being checked.
 */

#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Test: a plain bracketed list of three doubles, the canonical "D = [x,y,z]"
//       use case from the molecule information files.
// -----------------------------------------------------------------------------
void pia_test_bracketed_triplet()
{
    std::cerr << "\n[TEST] pia_test_bracketed_triplet\n"
              << "  Source file: src/parser/parse_input_array.cpp\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[1.0,2.0,3.0]\"\n"
              << "  Criteria:    returns exactly {1.0, 2.0, 3.0}.\n";

    std::string line { "[1.0,2.0,3.0]" };
    std::vector<double> vals = parse_input_array(line);

    // Exactly three values must come back.
    ASSERT_EQ(vals.size(), 3u) << "Three comma separated values were supplied";

    EXPECT_DOUBLE_EQ(vals[0], 1.0) << "First element should be 1.0";
    EXPECT_DOUBLE_EQ(vals[1], 2.0) << "Second element should be 2.0";
    EXPECT_DOUBLE_EQ(vals[2], 3.0) << "Third element should be 3.0";

    std::cerr << "  Parsed values: " << vals[0] << ", " << vals[1] << ", " << vals[2] << '\n';
}

// -----------------------------------------------------------------------------
// Test: parentheses are stripped exactly like square brackets.
// -----------------------------------------------------------------------------
void pia_test_parenthesis_delimiters()
{
    std::cerr << "\n[TEST] pia_test_parenthesis_delimiters\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"(4,5,6)\"\n"
              << "  Criteria:    '(' and ')' are removed, values parse to {4,5,6}.\n";

    std::string line { "(4,5,6)" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 3u) << "Three values expected from a parenthesised list";
    EXPECT_DOUBLE_EQ(vals[0], 4.0) << "First element should be 4.0";
    EXPECT_DOUBLE_EQ(vals[1], 5.0) << "Second element should be 5.0";
    EXPECT_DOUBLE_EQ(vals[2], 6.0) << "Third element should be 6.0";

    std::cerr << "  Parsed values: " << vals[0] << ", " << vals[1] << ", " << vals[2] << '\n';
}

// -----------------------------------------------------------------------------
// Test: no brackets at all -- the function must still split on commas.
// -----------------------------------------------------------------------------
void pia_test_no_brackets()
{
    std::cerr << "\n[TEST] pia_test_no_brackets\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"7,8\"  (no surrounding brackets)\n"
              << "  Criteria:    still parses into two values {7, 8}.\n";

    std::string line { "7,8" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 2u) << "Two comma separated values were supplied";
    EXPECT_DOUBLE_EQ(vals[0], 7.0) << "First element should be 7.0";
    EXPECT_DOUBLE_EQ(vals[1], 8.0) << "Second element should be 8.0";

    std::cerr << "  Parsed values: " << vals[0] << ", " << vals[1] << '\n';
}

// -----------------------------------------------------------------------------
// Test: a single value (no comma present at all).  The while-loop never fires
//       and the whole string is pushed as one token.
// -----------------------------------------------------------------------------
void pia_test_single_value()
{
    std::cerr << "\n[TEST] pia_test_single_value\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[42.5]\"\n"
              << "  Criteria:    returns a one-element vector holding 42.5.\n";

    std::string line { "[42.5]" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 1u) << "A single value should give a single element";
    EXPECT_DOUBLE_EQ(vals[0], 42.5) << "The lone element should be 42.5";

    std::cerr << "  Parsed value: " << vals[0] << '\n';
}

// -----------------------------------------------------------------------------
// Test: negative numbers, scientific notation, and integers without a decimal
//       point all flow through std::stod().
// -----------------------------------------------------------------------------
void pia_test_negative_and_scientific()
{
    std::cerr << "\n[TEST] pia_test_negative_and_scientific\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[-1.5,2e3,0,-0.001]\"\n"
              << "  Criteria:    negatives, exponents and bare integers parse correctly.\n";

    std::string line { "[-1.5,2e3,0,-0.001]" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 4u) << "Four comma separated values were supplied";
    EXPECT_DOUBLE_EQ(vals[0], -1.5) << "Negative decimal should parse";
    EXPECT_DOUBLE_EQ(vals[1], 2000.0) << "Scientific notation 2e3 should parse to 2000";
    EXPECT_DOUBLE_EQ(vals[2], 0.0) << "Bare integer 0 should parse to 0.0";
    EXPECT_DOUBLE_EQ(vals[3], -0.001) << "Small negative value should parse";

    std::cerr << "  Parsed values: " << vals[0] << ", " << vals[1] << ", " << vals[2] << ", "
              << vals[3] << '\n';
}

// -----------------------------------------------------------------------------
// Test: whitespace around tokens.  std::stod() skips leading whitespace and
//       simply stops at the first character it cannot consume, so " 1.5 "
//       parses fine.
// -----------------------------------------------------------------------------
void pia_test_whitespace_tolerated()
{
    std::cerr << "\n[TEST] pia_test_whitespace_tolerated\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[ 1.5 , -2.5 , 3.5 ]\"\n"
              << "  Criteria:    std::stod skips leading/trailing blanks -> {1.5,-2.5,3.5}.\n";

    std::string line { "[ 1.5 , -2.5 , 3.5 ]" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 3u) << "Three padded values were supplied";
    EXPECT_DOUBLE_EQ(vals[0], 1.5) << "Leading whitespace must not break parsing";
    EXPECT_DOUBLE_EQ(vals[1], -2.5) << "Whitespace around a negative must not break parsing";
    EXPECT_DOUBLE_EQ(vals[2], 3.5) << "Trailing whitespace must not break parsing";

    std::cerr << "  Parsed values: " << vals[0] << ", " << vals[1] << ", " << vals[2] << '\n';
}

// -----------------------------------------------------------------------------
// Test: the symbolic pi keywords.  "pi" and "M_PI" both make std::stod throw,
//       and the handler substitutes M_PI.  The comparison is case-insensitive
//       because the token is lower-cased before comparison.
// -----------------------------------------------------------------------------
void pia_test_pi_keyword()
{
    std::cerr << "\n[TEST] pia_test_pi_keyword\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[M_PI,pi,Pi,PI]\"\n"
              << "  Criteria:    every spelling maps to M_PI (case-insensitive).\n";

    std::string line { "[M_PI,pi,Pi,PI]" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 4u) << "Four pi spellings were supplied";
    EXPECT_DOUBLE_EQ(vals[0], M_PI) << "\"M_PI\" should become M_PI";
    EXPECT_DOUBLE_EQ(vals[1], M_PI) << "\"pi\" should become M_PI";
    EXPECT_DOUBLE_EQ(vals[2], M_PI) << "\"Pi\" should become M_PI (case-insensitive)";
    EXPECT_DOUBLE_EQ(vals[3], M_PI) << "\"PI\" should become M_PI (case-insensitive)";

    std::cerr << "  Parsed values: " << vals[0] << ", " << vals[1] << ", " << vals[2] << ", "
              << vals[3] << '\n';
}

// -----------------------------------------------------------------------------
// Test: the "nan" keyword.  Note that std::stod() itself already accepts the
//       literal "nan" and returns a quiet NaN, so the explicit isNull branch is
//       not necessarily what supplies the value -- but either path yields NaN,
//       which is exactly what the caller (association angles) relies upon.
// -----------------------------------------------------------------------------
void pia_test_nan_keyword()
{
    std::cerr << "\n[TEST] pia_test_nan_keyword\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[nan,1.0,NAN]\"\n"
              << "  Criteria:    the nan tokens come back as NaN, the numeric one as 1.0.\n";

    std::string line { "[nan,1.0,NAN]" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 3u) << "Three tokens were supplied";
    EXPECT_TRUE(std::isnan(vals[0])) << "Lower-case \"nan\" should produce a NaN";
    EXPECT_DOUBLE_EQ(vals[1], 1.0) << "The numeric token between the NaNs should be 1.0";
    EXPECT_TRUE(std::isnan(vals[2])) << "Upper-case \"NAN\" should also produce a NaN";

    std::cerr << "  Parsed values: " << vals[0] << ", " << vals[1] << ", " << vals[2]
              << " (NaN checks via std::isnan)\n";
}

// -----------------------------------------------------------------------------
// Test: a realistic five-element association-angle array mixing numbers and
//       symbolic keywords, i.e. the way `assocAngles` is written in .inp files.
// -----------------------------------------------------------------------------
void pia_test_assoc_angle_style_input()
{
    std::cerr << "\n[TEST] pia_test_assoc_angle_style_input\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[1.5707963,M_PI,nan,0,pi]\" (assocAngles style)\n"
              << "  Criteria:    five elements, symbolic entries resolved, order preserved.\n";

    std::string line { "[1.5707963,M_PI,nan,0,pi]" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 5u) << "assocAngles always supplies five values";
    EXPECT_DOUBLE_EQ(vals[0], 1.5707963) << "theta1 numeric value preserved";
    EXPECT_DOUBLE_EQ(vals[1], M_PI) << "theta2 given as M_PI";
    EXPECT_TRUE(std::isnan(vals[2])) << "phi1 given as nan should be NaN";
    EXPECT_DOUBLE_EQ(vals[3], 0.0) << "phi2 given as 0";
    EXPECT_DOUBLE_EQ(vals[4], M_PI) << "omega given as pi";

    std::cerr << "  Parsed values: " << vals[0] << ", " << vals[1] << ", " << vals[2] << ", "
              << vals[3] << ", " << vals[4] << '\n';
}

// -----------------------------------------------------------------------------
// Test: the documented side effect on the caller's string.  `line` is taken by
//       non-const reference and is consumed by the split loop, so after the call
//       it contains only the final token with brackets stripped.
// -----------------------------------------------------------------------------
void pia_test_input_string_is_consumed()
{
    std::cerr << "\n[TEST] pia_test_input_string_is_consumed\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[10,20,30]\"\n"
              << "  Criteria:    the by-reference argument is left holding only \"30\"\n"
              << "               (brackets stripped, earlier tokens erased).\n";

    std::string line { "[10,20,30]" };
    std::vector<double> vals = parse_input_array(line);

    ASSERT_EQ(vals.size(), 3u) << "Three values expected";
    EXPECT_DOUBLE_EQ(vals[2], 30.0) << "Final element should be 30.0";

    // Document/verify the destructive side effect on the caller's string.
    EXPECT_EQ(line, std::string("30"))
        << "parse_input_array consumes its argument; only the last token remains";

    std::cerr << "  Residual contents of the input string: \"" << line << "\"\n";
}

// -----------------------------------------------------------------------------
// Test: an unparsable token.  std::stod throws, the token is neither pi nor
//       nan, so a std::string is thrown out of the function.  (The trailing
//       `catch (const std::string&)` clause belongs to the same try block and
//       therefore cannot catch an exception raised inside a sibling handler, so
//       the exit(1) path is NOT taken -- it is safe to check this from a test.)
// -----------------------------------------------------------------------------
void pia_test_invalid_token_throws()
{
    std::cerr << "\n[TEST] pia_test_invalid_token_throws\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[1.0,bogus,3.0]\"\n"
              << "  Criteria:    a std::string exception escapes the function for a\n"
              << "               token that is neither a number, \"pi\" nor \"nan\".\n";

    std::string line { "[1.0,bogus,3.0]" };

    // EXPECT_THROW (not ASSERT_THROW) so that a failure here does not abort the
    // remainder of the suite.
    EXPECT_THROW(
        {
            std::vector<double> unused = parse_input_array(line);
            (void)unused;
        },
        std::string)
        << "An unrecognised token must throw a std::string error message";

    // Grab the message once more so it can be reported to the log.
    std::string line2 { "[1.0,bogus,3.0]" };
    try {
        std::vector<double> unused = parse_input_array(line2);
        (void)unused;
        ADD_FAILURE() << "parse_input_array did not throw for the invalid token \"bogus\"";
    } catch (const std::string& msg) {
        std::cerr << "  Caught expected error message: \"" << msg << "\"\n";
        EXPECT_NE(msg.find("bogus"), std::string::npos)
            << "The thrown message should name the offending token";
    } catch (...) {
        ADD_FAILURE() << "parse_input_array threw an unexpected exception type";
    }
}

// -----------------------------------------------------------------------------
// Test: an empty string.  The split loop produces a single empty token, stod
//       throws, and (since "" is neither pi nor nan) a std::string is thrown.
// -----------------------------------------------------------------------------
void pia_test_empty_string_throws()
{
    std::cerr << "\n[TEST] pia_test_empty_string_throws\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"\" (empty string)\n"
              << "  Criteria:    the empty token cannot be converted, so a std::string\n"
              << "               exception is thrown rather than an empty vector returned.\n";

    std::string line {};

    EXPECT_THROW(
        {
            std::vector<double> unused = parse_input_array(line);
            (void)unused;
        },
        std::string)
        << "An empty input string yields one empty token, which must throw";

    std::cerr << "  Empty input correctly reported an error instead of silently succeeding.\n";
}

// -----------------------------------------------------------------------------
// Test: a trailing comma produces an extra empty token, which likewise throws.
//       This documents that the parser is strict about trailing separators.
// -----------------------------------------------------------------------------
void pia_test_trailing_comma_throws()
{
    std::cerr << "\n[TEST] pia_test_trailing_comma_throws\n"
              << "  Function:    parse_input_array()\n"
              << "  Input:       \"[1.0,2.0,]\" (trailing comma)\n"
              << "  Criteria:    the empty token after the final comma cannot be parsed,\n"
              << "               so a std::string exception is thrown.\n";

    std::string line { "[1.0,2.0,]" };

    EXPECT_THROW(
        {
            std::vector<double> unused = parse_input_array(line);
            (void)unused;
        },
        std::string)
        << "A trailing comma leaves an empty final token, which must throw";

    std::cerr << "  Trailing comma correctly rejected.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so the
// framework reports them individually while still running every one of them.
// -----------------------------------------------------------------------------
TEST(ParseInputArrayTest, BracketedTriplet) { pia_test_bracketed_triplet(); }
TEST(ParseInputArrayTest, ParenthesisDelimiters) { pia_test_parenthesis_delimiters(); }
TEST(ParseInputArrayTest, NoBrackets) { pia_test_no_brackets(); }
TEST(ParseInputArrayTest, SingleValue) { pia_test_single_value(); }
TEST(ParseInputArrayTest, NegativeAndScientific) { pia_test_negative_and_scientific(); }
TEST(ParseInputArrayTest, WhitespaceTolerated) { pia_test_whitespace_tolerated(); }
TEST(ParseInputArrayTest, PiKeyword) { pia_test_pi_keyword(); }
TEST(ParseInputArrayTest, NanKeyword) { pia_test_nan_keyword(); }
TEST(ParseInputArrayTest, AssocAngleStyleInput) { pia_test_assoc_angle_style_input(); }
TEST(ParseInputArrayTest, InputStringIsConsumed) { pia_test_input_string_is_consumed(); }
TEST(ParseInputArrayTest, InvalidTokenThrows) { pia_test_invalid_token_throws(); }
TEST(ParseInputArrayTest, EmptyStringThrows) { pia_test_empty_string_throws(); }
TEST(ParseInputArrayTest, TrailingCommaThrows) { pia_test_trailing_comma_throws(); }