/*! \file error_codes.hpp
 * \brief Structured error categories and process exit codes.
 */

#pragma once

namespace nerdss {
namespace error {

enum class ErrorCategory {
  success,
  input,
  file_io,
  reaction,
  state,
  numeric,
  dependency,
  mpi,
  internal,
  unsupported,
  interrupted,
  unknown,
};

enum class ExitCode : int {
  success = 0,
  general = 1,
  input = 2,
  coordinate_input = 3,
  invalid_container_value = 4,
  unbalanced_reaction = 5,
  invalid_reaction = 6,
  dependency = 7,
  numeric = 8,
  file_io = 9,
  unsupported = 10,
  internal = 11,
  invariant = 12,
  mpi = 13,
  parser_state = 120,
  unknown = 125,
  interrupted = 130,
};

inline int to_exit_status(ExitCode code) { return static_cast<int>(code); }

inline const char *to_string(ErrorCategory category) {
  switch (category) {
  case ErrorCategory::success:
    return "success";
  case ErrorCategory::input:
    return "input";
  case ErrorCategory::file_io:
    return "file_io";
  case ErrorCategory::reaction:
    return "reaction";
  case ErrorCategory::state:
    return "state";
  case ErrorCategory::numeric:
    return "numeric";
  case ErrorCategory::dependency:
    return "dependency";
  case ErrorCategory::mpi:
    return "mpi";
  case ErrorCategory::internal:
    return "internal";
  case ErrorCategory::unsupported:
    return "unsupported";
  case ErrorCategory::interrupted:
    return "interrupted";
  case ErrorCategory::unknown:
    return "unknown";
  }

  return "unknown";
}

inline const char *to_string(ExitCode code) {
  switch (code) {
  case ExitCode::success:
    return "success";
  case ExitCode::general:
    return "general";
  case ExitCode::input:
    return "input";
  case ExitCode::coordinate_input:
    return "coordinate_input";
  case ExitCode::invalid_container_value:
    return "invalid_container_value";
  case ExitCode::unbalanced_reaction:
    return "unbalanced_reaction";
  case ExitCode::invalid_reaction:
    return "invalid_reaction";
  case ExitCode::dependency:
    return "dependency";
  case ExitCode::numeric:
    return "numeric";
  case ExitCode::file_io:
    return "file_io";
  case ExitCode::unsupported:
    return "unsupported";
  case ExitCode::internal:
    return "internal";
  case ExitCode::invariant:
    return "invariant";
  case ExitCode::mpi:
    return "mpi";
  case ExitCode::parser_state:
    return "parser_state";
  case ExitCode::unknown:
    return "unknown";
  case ExitCode::interrupted:
    return "interrupted";
  }

  return "unknown";
}

inline ExitCode default_exit_code(ErrorCategory category) {
  switch (category) {
  case ErrorCategory::success:
    return ExitCode::success;
  case ErrorCategory::input:
    return ExitCode::input;
  case ErrorCategory::file_io:
    return ExitCode::file_io;
  case ErrorCategory::reaction:
    return ExitCode::invalid_reaction;
  case ErrorCategory::state:
    return ExitCode::invalid_container_value;
  case ErrorCategory::numeric:
    return ExitCode::numeric;
  case ErrorCategory::dependency:
    return ExitCode::dependency;
  case ErrorCategory::mpi:
    return ExitCode::mpi;
  case ErrorCategory::internal:
    return ExitCode::internal;
  case ErrorCategory::unsupported:
    return ExitCode::unsupported;
  case ErrorCategory::interrupted:
    return ExitCode::interrupted;
  case ErrorCategory::unknown:
    return ExitCode::unknown;
  }

  return ExitCode::unknown;
}

} // namespace error
} // namespace nerdss
