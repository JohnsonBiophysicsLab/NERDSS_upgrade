# Google Test Migration Plan

Agent Q owns the Phase 7 unit and integration test framework. The framework is
being migrated to Google Test and C++ tests should follow the Google C++ Style
Guide.

## Current Freeze

Do not add new custom CTest/std-library unit or integration tests while this
migration is active. Existing CTest-era assertions may be preserved on their
branches for reference, but reviewable follow-up PRs should re-spin them as
Google Test cases.

## Migration Order

1. Add a Google Test dependency/bootstrap strategy that works in GitHub Actions
   and local CMake builds without changing the scientific executable.
2. Add a minimal Google Test smoke target and document local commands.
3. Split fast unit tests from slower integration/validation examples.
4. Convert useful CTest-era assertions to Google Test style in focused PRs.
5. Update CI and coverage collection to call the Google Test targets.

## Superseded CTest Work

The custom CTest harness, validation-integration stack PRs, and CTest-oriented
unit/integration commands are superseded by this plan. The branches are kept so
implementation work and useful assertions can be cherry-picked into Google Test
PRs later.

Generated build directories such as `build-unit-tests/` should not be
committed.
