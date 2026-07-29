/*! \file test_eye_candy.cpp
 *
 * ### Unit tests for ../src/io/eye_candy.cpp
 *
 * The file under test provides four stream "manipulator" style free functions
 * that are used throughout NERDSS to prettify console output:
 *
 *     std::ostream& bon(std::ostream& os);        // turn ANSI bold ON
 *     std::ostream& boff(std::ostream& os);       // turn ANSI bold OFF
 *     std::ostream& llinebreak(std::ostream& os); // "long"  dashed line break
 *     std::ostream& linebreak(std::ostream& os);  // "short" dashed line break
 *
 * Because each function simply writes characters into the stream it is handed
 * and returns that same stream, they can be tested deterministically by
 * streaming them into a std::ostringstream and inspecting the resulting bytes.
 *
 * Pass criteria used below:
 *   - bon()        emits exactly the ANSI "bold on"  sequence  ESC [ 1 m
 *   - boff()       emits exactly the ANSI "bold off" sequence  ESC [ 0 m
 *   - llinebreak() emits a 50-wide right-justified ' ' padded with '-'
 *                  (i.e. 49 dashes + a space) followed by '\n'  => 51 chars
 *   - linebreak()  emits a 20-wide right-justified ' ' padded with '-'
 *                  (i.e. 19 dashes + a space) followed by '\n'  => 21 chars
 *   - every function returns a reference to the *same* stream it was given
 *   - the stream fill character is restored to ' ' by the line break functions
 *
 * Verbose progress messages are written to stderr so that a reader of the test
 * log can see exactly which source file / function is being exercised.
 */

#include "io/io.hpp"
#include <iomanip>
#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>

namespace {

//! ANSI escape character used by the bold on/off manipulators.
const char kEyeCandyEsc = '\x1b'; // same as '\e' GNU extension in the source

/*! \brief Turn a raw string into a printable form so escape characters show up
 *         in the test log rather than actually colouring the terminal.
 *
 * \param[in] raw the exact bytes produced by the function under test
 * \return a human readable representation ("<ESC>" for 0x1b, "\\n" for newline)
 */
std::string eyecandy_escape_for_log(const std::string& raw)
{
    std::string out;
    for (char c : raw) {
        if (c == kEyeCandyEsc)
            out += "<ESC>";
        else if (c == '\n')
            out += "\\n";
        else
            out += c;
    }
    return out;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: bon() -- ANSI bold ON manipulator.
// -----------------------------------------------------------------------------
void test_eyecandy_bon()
{
    std::cerr << "\n[TEST] test_eyecandy_bon\n"
              << "  Source file: src/io/eye_candy.cpp\n"
              << "  Function:    std::ostream& bon(std::ostream&)\n"
              << "  Action:      stream bon into an empty ostringstream\n"
              << "  Criteria:    output is exactly ESC '[' '1' 'm' (4 bytes)\n";

#if defined(__APPLE__) || defined(__linux__)
    std::ostringstream oss;

    // Stream the manipulator; this calls bon(oss) via the ostream operator<<
    // overload that accepts a function pointer.
    oss << bon;

    const std::string produced = oss.str();
    std::cerr << "  Produced:    \"" << eyecandy_escape_for_log(produced) << "\""
              << " (" << produced.size() << " bytes)\n";

    // Expected byte sequence for "bold on".
    const std::string expected = std::string(1, kEyeCandyEsc) + "[1m";

    EXPECT_EQ(produced, expected)
        << "bon() must emit the ANSI bold-on escape sequence";
    EXPECT_EQ(produced.size(), static_cast<size_t>(4))
        << "bon() should emit exactly 4 characters";

    // Sanity check the individual bytes so a failure pinpoints the problem.
    if (produced.size() == 4) {
        EXPECT_EQ(produced[0], kEyeCandyEsc) << "first byte must be ESC (0x1b)";
        EXPECT_EQ(produced[1], '[') << "second byte must be '['";
        EXPECT_EQ(produced[2], '1') << "third byte must be '1' (bold)";
        EXPECT_EQ(produced[3], 'm') << "fourth byte must be 'm'";
    }
#else
    std::cerr << "  SKIPPED: bon() is only declared on macOS/Linux builds.\n";
#endif
}

// -----------------------------------------------------------------------------
// Test 2: boff() -- ANSI bold OFF (attribute reset) manipulator.
// -----------------------------------------------------------------------------
void test_eyecandy_boff()
{
    std::cerr << "\n[TEST] test_eyecandy_boff\n"
              << "  Source file: src/io/eye_candy.cpp\n"
              << "  Function:    std::ostream& boff(std::ostream&)\n"
              << "  Action:      stream boff into an empty ostringstream\n"
              << "  Criteria:    output is exactly ESC '[' '0' 'm' (4 bytes)\n";

#if defined(__APPLE__) || defined(__linux__)
    std::ostringstream oss;
    oss << boff;

    const std::string produced = oss.str();
    std::cerr << "  Produced:    \"" << eyecandy_escape_for_log(produced) << "\""
              << " (" << produced.size() << " bytes)\n";

    const std::string expected = std::string(1, kEyeCandyEsc) + "[0m";

    EXPECT_EQ(produced, expected)
        << "boff() must emit the ANSI attribute-reset escape sequence";
    EXPECT_EQ(produced.size(), static_cast<size_t>(4))
        << "boff() should emit exactly 4 characters";

    if (produced.size() == 4) {
        EXPECT_EQ(produced[0], kEyeCandyEsc) << "first byte must be ESC (0x1b)";
        EXPECT_EQ(produced[1], '[') << "second byte must be '['";
        EXPECT_EQ(produced[2], '0') << "third byte must be '0' (reset)";
        EXPECT_EQ(produced[3], 'm') << "fourth byte must be 'm'";
    }

    // A bold-on immediately followed by bold-off must round-trip cleanly, which
    // is how the manipulators are actually used in NERDSS output.
    std::ostringstream pair;
    pair << bon << "text" << boff;
    const std::string pairStr = pair.str();
    std::cerr << "  Combined:    \"" << eyecandy_escape_for_log(pairStr) << "\"\n";
    EXPECT_EQ(pairStr, std::string(1, kEyeCandyEsc) + "[1m" + "text"
                  + std::string(1, kEyeCandyEsc) + "[0m")
        << "bon/boff should wrap the payload text exactly";
#else
    std::cerr << "  SKIPPED: boff() is only declared on macOS/Linux builds.\n";
#endif
}

// -----------------------------------------------------------------------------
// Test 3: llinebreak() -- the "long" 50 column dashed line break.
// -----------------------------------------------------------------------------
void test_eyecandy_llinebreak()
{
    std::cerr << "\n[TEST] test_eyecandy_llinebreak\n"
              << "  Source file: src/io/eye_candy.cpp\n"
              << "  Function:    std::ostream& llinebreak(std::ostream&)\n"
              << "  Action:      stream llinebreak into an empty ostringstream\n"
              << "  Criteria:    49 '-' + ' ' + '\\n'  (51 characters total)\n";

    std::ostringstream oss;
    oss << llinebreak;

    const std::string produced = oss.str();
    std::cerr << "  Produced:    \"" << eyecandy_escape_for_log(produced) << "\""
              << " (" << produced.size() << " bytes)\n";

    // setw(50) with default right-justification pads the single ' ' character
    // with 49 fill characters ('-'), then a newline is appended.
    const std::string expected = std::string(49, '-') + " \n";

    EXPECT_EQ(produced, expected)
        << "llinebreak() should produce a 50 wide dashed field plus newline";
    EXPECT_EQ(produced.size(), static_cast<size_t>(51))
        << "llinebreak() total length should be 50 (field) + 1 (newline)";

    // Verify the character composition explicitly.
    const size_t dashCount = produced.find_first_not_of('-');
    EXPECT_EQ(dashCount, static_cast<size_t>(49))
        << "there should be exactly 49 leading dashes";
    if (produced.size() == 51) {
        EXPECT_EQ(produced[49], ' ') << "character 50 should be the padded space";
        EXPECT_EQ(produced[50], '\n') << "the line break must end with a newline";
    }

    // The function restores the fill character to ' '; confirm the very next
    // width-limited write is space padded, not dash padded.
    oss << std::setw(4) << 'X';
    const std::string after = oss.str().substr(produced.size());
    std::cerr << "  Follow-up write with setw(4): \""
              << eyecandy_escape_for_log(after) << "\"\n";
    EXPECT_EQ(after, "   X")
        << "llinebreak() must restore the stream fill character to a space";
}

// -----------------------------------------------------------------------------
// Test 4: linebreak() -- the "short" 20 column dashed line break.
// -----------------------------------------------------------------------------
void test_eyecandy_linebreak()
{
    std::cerr << "\n[TEST] test_eyecandy_linebreak\n"
              << "  Source file: src/io/eye_candy.cpp\n"
              << "  Function:    std::ostream& linebreak(std::ostream&)\n"
              << "  Action:      stream linebreak into an empty ostringstream\n"
              << "  Criteria:    19 '-' + ' ' + '\\n'  (21 characters total)\n";

    std::ostringstream oss;
    oss << linebreak;

    const std::string produced = oss.str();
    std::cerr << "  Produced:    \"" << eyecandy_escape_for_log(produced) << "\""
              << " (" << produced.size() << " bytes)\n";

    const std::string expected = std::string(19, '-') + " \n";

    EXPECT_EQ(produced, expected)
        << "linebreak() should produce a 20 wide dashed field plus newline";
    EXPECT_EQ(produced.size(), static_cast<size_t>(21))
        << "linebreak() total length should be 20 (field) + 1 (newline)";

    const size_t dashCount = produced.find_first_not_of('-');
    EXPECT_EQ(dashCount, static_cast<size_t>(19))
        << "there should be exactly 19 leading dashes";
    if (produced.size() == 21) {
        EXPECT_EQ(produced[19], ' ') << "character 20 should be the padded space";
        EXPECT_EQ(produced[20], '\n') << "the line break must end with a newline";
    }

    // Fill character restoration check, as above.
    oss << std::setw(3) << 'Y';
    const std::string after = oss.str().substr(produced.size());
    std::cerr << "  Follow-up write with setw(3): \""
              << eyecandy_escape_for_log(after) << "\"\n";
    EXPECT_EQ(after, "  Y")
        << "linebreak() must restore the stream fill character to a space";

    // The short break must be strictly shorter than the long break; this guards
    // against the two width constants being accidentally swapped.
    std::ostringstream longOss;
    longOss << llinebreak;
    std::cerr << "  Short break length = " << produced.size()
              << ", long break length = " << longOss.str().size() << '\n';
    EXPECT_LT(produced.size(), longOss.str().size())
        << "linebreak() must be shorter than llinebreak()";
}

// -----------------------------------------------------------------------------
// Test 5: every manipulator returns a reference to the stream it was handed, so
//         calls can be chained.  This is checked both by comparing addresses
//         from a direct call and by chaining several manipulators together.
// -----------------------------------------------------------------------------
void test_eyecandy_returns_same_stream_and_chaining()
{
    std::cerr << "\n[TEST] test_eyecandy_returns_same_stream_and_chaining\n"
              << "  Source file: src/io/eye_candy.cpp\n"
              << "  Functions:   bon, boff, llinebreak, linebreak\n"
              << "  Action:      call each directly and compare returned address\n"
              << "               to the stream address; then chain them.\n"
              << "  Criteria:    &returned == &original, and chained output is\n"
              << "               the exact concatenation of the pieces.\n";

    std::ostringstream oss;
    std::ostream& base = oss;

    // Direct calls: each function must hand back the very same stream object.
    std::cerr << "  -> checking llinebreak() return identity\n";
    EXPECT_EQ(&llinebreak(base), &base)
        << "llinebreak() must return the stream it was given";

    std::cerr << "  -> checking linebreak() return identity\n";
    EXPECT_EQ(&linebreak(base), &base)
        << "linebreak() must return the stream it was given";

#if defined(__APPLE__) || defined(__linux__)
    std::cerr << "  -> checking bon() return identity\n";
    EXPECT_EQ(&bon(base), &base) << "bon() must return the stream it was given";

    std::cerr << "  -> checking boff() return identity\n";
    EXPECT_EQ(&boff(base), &base) << "boff() must return the stream it was given";
#else
    std::cerr << "  -> SKIPPED bon()/boff() identity checks (non macOS/Linux)\n";
#endif

    // Now verify that chaining produces the concatenation of each piece.
    std::ostringstream chained;
    chained << linebreak << llinebreak;

    const std::string expected = std::string(19, '-') + " \n" + std::string(49, '-') + " \n";
    std::cerr << "  Chained output length = " << chained.str().size()
              << ", expected " << expected.size() << '\n';
    EXPECT_EQ(chained.str(), expected)
        << "chaining linebreak then llinebreak must concatenate their output";

    // Make sure the chained stream is still in a good state (no error bits).
    EXPECT_TRUE(chained.good())
        << "the stream must remain in a good state after the manipulators run";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper is invoked from its own TEST so that
// a failure in one does not prevent the others from running (all assertions
// above are non-fatal EXPECT_* checks).
// -----------------------------------------------------------------------------
TEST(EyeCandy, BoldOn) { test_eyecandy_bon(); }
TEST(EyeCandy, BoldOff) { test_eyecandy_boff(); }
TEST(EyeCandy, LongLineBreak) { test_eyecandy_llinebreak(); }
TEST(EyeCandy, ShortLineBreak) { test_eyecandy_linebreak(); }
TEST(EyeCandy, ReturnsSameStreamAndChaining) { test_eyecandy_returns_same_stream_and_chaining(); }