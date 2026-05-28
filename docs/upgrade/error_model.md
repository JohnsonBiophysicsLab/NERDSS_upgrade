# Structured Error Model

This document defines the first Phase 5 target for structured NERDSS errors.
It is a planning and compatibility layer only: existing simulation behavior,
stderr formatting, MPI handling, and process exits are not changed by this
slice.

## Goals

- Give future parser, setup, simulation, and I/O errors stable categories.
- Reserve process exit codes so command-line callers can distinguish failure
  classes without parsing human-readable messages.
- Preserve current behavior during migration. Existing `error(...)` calls and
  direct `exit(...)` calls remain valid until their owning code paths are
  intentionally migrated.
- Keep rank-local context available for MPI diagnostics without requiring this
  first slice to change MPI finalization or collective behavior.

## Non-goals

- No broad source rewiring in Phase 5 slice 1.
- No exception-policy change for core simulation loops.
- No change to existing output text, logs, restart files, random seeds, or
  trajectory behavior.
- No requirement that all legacy exit codes immediately match the new table.

## Error Categories

`include/error/error_codes.hpp` defines `nerdss::error::ErrorCategory` as the
common vocabulary for future structured diagnostics:

| Category | Intended use |
| --- | --- |
| `success` | Completed successfully. |
| `input` | Invalid command-line argument, input file syntax, or semantic input validation. |
| `file_io` | Missing, unreadable, unwritable, or malformed external file. |
| `reaction` | Invalid reaction definition or reaction consistency failure. |
| `state` | Invalid molecule, complex, interface, observable, or restart state. |
| `numeric` | Numerical integration, probability, geometry, or floating-point failure. |
| `dependency` | Required external library, runtime feature, or environment dependency is unavailable. |
| `mpi` | MPI rank, communication, decomposition, or parallel consistency failure. |
| `internal` | Failed internal invariant or implementation bug. |
| `unsupported` | Recognized feature or option that is not implemented in this build or mode. |
| `interrupted` | User or system interruption. |
| `unknown` | Fallback when no more specific category is known. |

## Exit Codes

The first structured exit-code set keeps values close to existing NERDSS
conventions where they are already present. Future migrations should prefer
these values over new ad hoc integer literals.

| Exit code | Name | Category | Notes |
| --- | --- | --- | --- |
| 0 | `success` | `success` | Normal completion. |
| 1 | `general` | `unknown` | Legacy fallback used by `error(...)` and many direct exits. |
| 2 | `input` | `input` | General input parse or validation failure. |
| 3 | `coordinate_input` | `file_io` | Coordinate-specific read failure, matching a legacy helper. |
| 4 | `invalid_container_value` | `state` | Required value missing from an internal container. |
| 5 | `unbalanced_reaction` | `reaction` | Reaction definition is not balanced. |
| 6 | `invalid_reaction` | `reaction` | Invalid reaction rule or parameterization. |
| 7 | `dependency` | `dependency` | Missing external dependency or unavailable runtime capability. |
| 8 | `numeric` | `numeric` | Numerical failure. |
| 9 | `file_io` | `file_io` | General file I/O failure. |
| 10 | `unsupported` | `unsupported` | Unsupported but recognized feature or mode. |
| 11 | `internal` | `internal` | Internal bug or invariant failure. |
| 12 | `invariant` | `internal` | Legacy involvement/invariant error. |
| 13 | `mpi` | `mpi` | MPI or parallel consistency failure. |
| 120 | `parser_state` | `input` | Preserves the current state-parser special exit. |
| 125 | `unknown` | `unknown` | Structured fallback when a legacy code is not appropriate. |
| 130 | `interrupted` | `interrupted` | Conventional Ctrl-C/SIGINT-style termination. |

## Migration Rules

1. Add structured metadata at module boundaries first, such as parser entry
   points, file readers, setup validation, and top-level executable handling.
2. Do not change stochastic, geometry, reaction, or scheduler behavior while
   converting error reporting.
3. Keep legacy messages stable unless a migration explicitly updates tests and
   caller documentation.
4. When replacing a direct `exit(...)`, choose the closest `ExitCode` value and
   document any intentional compatibility change in the owning pull request.
5. MPI call sites should retain rank information. Any future collective abort
   policy should be introduced separately from category assignment.

## Current Implementation Surface

`include/error/error_codes.hpp` is header-only and has no side effects. It
provides:

- `ErrorCategory`
- `ExitCode`
- `to_exit_status(ExitCode)`
- `to_string(ErrorCategory)`
- `to_string(ExitCode)`
- `default_exit_code(ErrorCategory)`

The header is not wired into existing simulation code yet.
