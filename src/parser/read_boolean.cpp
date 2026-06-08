#include "parser/parser_diagnostics.hpp"
#include "parser/parser_functions.hpp"

bool read_boolean(std::string fileLine) {
  const std::string originalLine{fileLine};

  // This function is to read in booleans of any form
  std::transform(fileLine.begin(), fileLine.end(), fileLine.begin(), ::tolower);

  // remove whitespaces
  remove_comment(fileLine);
  fileLine.erase(std::remove_if(fileLine.begin(), fileLine.end(),
                                [](unsigned char x) { return std::isspace(x); }),
                 fileLine.end());

  if (fileLine == "0" || fileLine == "false")
    return false;
  else if (fileLine == "1" || fileLine == "true")
    return true;
  else {
    nerdss::parser::fail_parser_error("read_boolean", "expected boolean token 0, 1, false, or true",
                                      originalLine);
  }
}
