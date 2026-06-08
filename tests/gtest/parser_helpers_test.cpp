#include "parser/parser_functions.hpp"
#include "parser/parser_diagnostics.hpp"

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

namespace {

constexpr double kTolerance = 1.0e-12;

TEST(ParserHelpersTest, RemovesInlineCommentsAndPreservesPrefixWhitespace) {
  std::string with_comment = "  value = 10  # trailing comment";
  remove_comment(with_comment);
  EXPECT_EQ(with_comment, "  value = 10  ");

  std::string without_comment = "  value = 10  ";
  remove_comment(without_comment);
  EXPECT_EQ(without_comment, "  value = 10  ");
}

TEST(ParserHelpersTest, ReadsBooleanFormsWithWhitespaceAndComments) {
  EXPECT_FALSE(read_boolean("0"));
  EXPECT_FALSE(read_boolean(" false "));
  EXPECT_FALSE(read_boolean("FALSE # disabled"));

  EXPECT_TRUE(read_boolean("1"));
  EXPECT_TRUE(read_boolean(" true "));
  EXPECT_TRUE(read_boolean("TRUE # enabled"));
}

TEST(ParserHelpersTest, ParsesNumericInputArraysAndMutatesLine) {
  std::string line = "[1.5, -2.0, 3]";

  std::vector<double> values = parse_input_array(line);

  ASSERT_EQ(values.size(), 3U);
  EXPECT_NEAR(values[0], 1.5, kTolerance);
  EXPECT_NEAR(values[1], -2.0, kTolerance);
  EXPECT_NEAR(values[2], 3.0, kTolerance);
  EXPECT_EQ(line, " 3");
}

TEST(ParserHelpersTest, ParsesPiAndNanTokens) {
  std::string line = "(pi,M_PI,nan)";

  std::vector<double> values = parse_input_array(line);

  ASSERT_EQ(values.size(), 3U);
  EXPECT_NEAR(values[0], std::acos(-1.0), kTolerance);
  EXPECT_NEAR(values[1], std::acos(-1.0), kTolerance);
  EXPECT_TRUE(std::isnan(values[2]));
}

TEST(ParserHelpersTest, FormatsParserDiagnosticsWithContext) {
  std::string message = nerdss::parser::format_parser_error(
      "read_boolean", "expected boolean token", "maybe");

  EXPECT_NE(message.find("PARSER_ERROR[read_boolean]"), std::string::npos);
  EXPECT_NE(message.find("expected boolean token"), std::string::npos);
  EXPECT_NE(message.find("input: maybe"), std::string::npos);
  EXPECT_NE(message.find("category: input"), std::string::npos);
  EXPECT_NE(message.find("exit_code: 2"), std::string::npos);
}

}  // namespace
