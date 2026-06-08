#pragma once

#include <string>

namespace nerdss {
namespace parser {

std::string format_parser_error(const std::string& parser_name, const std::string& detail,
                                const std::string& input);

[[noreturn]] void fail_parser_error(const std::string& parser_name, const std::string& detail,
                                    const std::string& input);

}  // namespace parser
}  // namespace nerdss
