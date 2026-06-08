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

## Initial Harness

The initial CMake hook is opt-in so normal serial and MPI builds do not require
Google Test:

```bash
cmake -S . -B build-gtest -DNERDSS_ENABLE_GTEST=ON
cmake --build build-gtest --target nerdss_gtest_smoke
./build-gtest/nerdss_gtest_smoke
```

`NERDSS_ENABLE_GTEST=ON` requires a discoverable Google Test package. On Ubuntu
CI this is provided by `libgtest-dev`. Local macOS builds can use a package
manager installation or a CMake prefix that provides `GTest::gtest_main`.
If the installed package exposes the older `GTest::Main` imported target, the
CMake harness uses that target instead.

Run the Google Test binary directly. Do not add `ctest` as the default runner
for this migration; a CTest wrapper can be reconsidered later if it adds value
without reintroducing the custom test executable.

## Superseded CTest Work

The custom CTest harness, validation-integration stack PRs, `BUILD_TESTING`
unit target, and CTest-oriented unit/integration commands are superseded by this
plan. The branches are kept so implementation work and useful assertions can be
cherry-picked into Google Test PRs later.

Generated build directories such as `build-unit-tests/` should not be
committed.
