# Static Analysis And Sanitizer Backlog

This backlog records current findings from the Phase 2 sanitizer/static-analysis
guardrail setup. Items here are intentionally not fixed in the tooling PR unless
the fix is limited to build configuration or documentation.

## Open Items

### P1: Serial CMake build fails on undeclared bonded-complex JSON symbols

- Reproduction:
  `cmake --build build/validate-asan --target nerdss --parallel 4`
- Observed on: 2026-05-28, `codex/upgrade-sanitizer-static` based on
  `origin/master`.
- Failure:
  - `EXEs/nerdss.cpp` references `params.bondedComplexWrite`, but
    `include/classes/class_Parameters.hpp` on `origin/master` does not declare
    that member.
  - `EXEs/nerdss.cpp` calls `write_bonded_complex_json`, but
    `include/io/io.hpp` on `origin/master` does not declare that function.
- Impact:
  serial sanitizer builds configure correctly but cannot complete until the
  source/header mismatch is resolved.
- Notes:
  `origin/codex/upgrade-smoke-runner` appears to add the
  `Parameters::bondedComplexWrite` declaration. Coordinate with that branch
  rather than duplicating source/header fixes in this tooling-only PR.

### P2: clang-tidy is not available in the local validation environment

- Reproduction: `tools/run_static_analysis.sh src/math`
- Failure: `clang-tidy was not found on PATH.`
- Impact:
  the focused static-analysis runner is present, but local findings cannot be
  generated until clang-tidy is installed.
- Suggested next step:
  install clang-tidy in the developer/CI image, then run:
  `tools/run_static_analysis.sh src/math src/parser EXEs`

### P3: CMake 4.3 warns about compatibility policy floor

- Reproduction:
  `cmake -S . -B build/validate-asan -DCMAKE_BUILD_TYPE=Debug -DNERDSS_ENABLE_ASAN=ON`
- Warning:
  CMake 4.3 warns that compatibility with CMake versions older than 3.10 will be
  removed in a future CMake release.
- Impact:
  configure succeeds after raising the repository floor from 3.0 to 3.5, but a
  future upgrade should decide whether NERDSS can require CMake 3.10 or newer.
- Suggested next step:
  include the CMake minimum in Agent A's environment policy decision.
