# Upgrade Validation Log

Record validation runs for upgrade work in reverse chronological order. Keep
entries concise, reproducible, and explicit about skipped checks.

## Entry Template

Date:

Branch:

Commit:

Workstream:

Environment:
- OS:
- Compiler:
- GSL:
- CMake:
- Make:
- MPI:

Commands:

Results:

Artifacts:

Notes:

## 2026-06-08: Parser Helpers Google Test Conversion

Date: 2026-06-08

Branch: `codex/gtest-parser-helpers`

Commit: Pending at validation log update

Workstream: Phase 7 Google Test migration and parser normalization
preservation tests

Environment:
- OS: macOS, local Codex workspace
- Compiler: Apple clang via default CMake compiler
- GSL: 2.8 from Homebrew
- CMake: Available on PATH
- Make: GNU Make
- MPI: Local `mpicxx` wrapper remains blocked by missing
  `x86_64-apple-darwin13.4.0-clang++`

Commands:
```sh
git diff --check HEAD
tools/format_changed_files.sh --check --base HEAD
cmake -S . -B /tmp/nerdss-gtest-default
cmake --build /tmp/nerdss-gtest-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-gtest-on -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-gtest-parser-helpers-smoke
```

Results:
- Passed: throwaway `/tmp/nerdss_parser_probe` compile/run for the target source
  set after preserving current special-token behavior (`pi`/`nan` are not
  whitespace-trimmed before token matching).
- Passed: `git diff --check HEAD`.
- Passed: `tools/format_changed_files.sh --check --base HEAD`; the helper
  reported no tracked C/C++ files to format before staging this new test file.
- Passed: default CMake configure in `/tmp/nerdss-gtest-default`.
- Passed: default CMake build target `nerdss`.
- Failed as expected locally: `NERDSS_ENABLE_GTEST=ON` configure could not find
  `GTEST_LIBRARY`, `GTEST_INCLUDE_DIR`, or `GTEST_MAIN_LIBRARY`.
- Passed: `make serial -j4`.
- Passed: serial smoke runner with artifacts in
  `/tmp/nerdss-gtest-parser-helpers-smoke`.
- Not run locally: `clang-format --dry-run --Werror
  tests/gtest/parser_helpers_test.cpp`; `clang-format` is not on PATH in this
  shell.

Artifacts:
- `CMakeLists.txt`
- `.github/workflows/c-cpp.yml`
- `tests/gtest/parser_helpers_test.cpp`
- `docs/upgrade/unit_tests.md`
- This log entry.

Notes:
- This slice adds deterministic coverage for `remove_comment`, `read_boolean`,
  and `parse_input_array`.
- Invalid-input fatal paths remain untested until the GTest migration has a
  settled pattern for death tests and error-message assertions.
- The Google Test target is not expected to build locally unless Google Test is
  installed or exposed through `CMAKE_PREFIX_PATH`.

## 2026-06-08: Association Angle Helpers Google Test Conversion

Date: 2026-06-08

Branch: `codex/gtest-association-angle-helpers`

Commit: Pending at validation log update

Workstream: Phase 7 Google Test migration and association geometry preservation
tests

Environment:
- OS: macOS, local Codex workspace
- Compiler: Apple clang via default CMake compiler
- GSL: 2.8 from Homebrew
- CMake: Available on PATH
- Make: GNU Make
- MPI: Local `mpicxx` wrapper remains blocked by missing
  `x86_64-apple-darwin13.4.0-clang++`

Commands:
```sh
git diff --check HEAD
tools/format_changed_files.sh --check --base HEAD
cmake -S . -B /tmp/nerdss-gtest-default
cmake --build /tmp/nerdss-gtest-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-gtest-on -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-gtest-association-angle-smoke
```

Results:
- Passed: throwaway `/tmp/nerdss_assoc_probe` compile/run to confirm current
  `requiresSignFlip` projection expectations with initialized vector
  magnitudes.
- Passed: `git diff --check HEAD`.
- Passed: `tools/format_changed_files.sh --check --base HEAD`; the helper
  reported no tracked C/C++ files to format before staging this new test file.
- Passed: default CMake configure in `/tmp/nerdss-gtest-default`.
- Passed: default CMake build target `nerdss`.
- Failed as expected locally: `NERDSS_ENABLE_GTEST=ON` configure could not find
  `GTEST_LIBRARY`, `GTEST_INCLUDE_DIR`, or `GTEST_MAIN_LIBRARY`.
- Passed: `make serial -j4`.
- Passed: serial smoke runner with artifacts in
  `/tmp/nerdss-gtest-association-angle-smoke`.
- Not run locally: `clang-format --dry-run --Werror
  tests/gtest/association_angle_helpers_test.cpp`; `clang-format` is not on PATH
  in this shell.

Artifacts:
- `CMakeLists.txt`
- `.github/workflows/c-cpp.yml`
- `tests/gtest/association_angle_helpers_test.cpp`
- `docs/upgrade/unit_tests.md`
- This log entry.

Notes:
- This slice adds deterministic coverage for association angle helper
  predicates used by theta/phi/omega rotation paths.
- The cases preserve current tolerance, exact-parallel sentinel, orthogonal
  fallback, and sign-flip projection behavior.
- The Google Test target is not expected to build locally unless Google Test is
  installed or exposed through `CMAKE_PREFIX_PATH`.

## 2026-06-08: Math Functions Google Test Conversion

Date: 2026-06-08

Branch: `codex/gtest-math-functions`

Commit: Pending at validation log update

Workstream: Phase 7 Google Test migration and MathEngine preservation tests

Environment:
- OS: macOS, local Codex workspace
- Compiler: Apple clang via default CMake compiler
- GSL: 2.8 from Homebrew
- CMake: Available on PATH
- Make: GNU Make
- MPI: Local `mpicxx` wrapper remains blocked by missing
  `x86_64-apple-darwin13.4.0-clang++`

Commands:
```sh
git diff --check HEAD
tools/format_changed_files.sh --check --base HEAD
cmake -S . -B /tmp/nerdss-gtest-default
cmake --build /tmp/nerdss-gtest-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-gtest-on -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-gtest-math-functions-smoke
```

Results:
- Passed: `git diff --check HEAD`.
- Passed: `tools/format_changed_files.sh --check --base HEAD`; the helper
  reported no tracked C/C++ files to format before staging this new test file.
- Passed: default CMake configure in `/tmp/nerdss-gtest-default`.
- Passed: default CMake build target `nerdss`.
- Failed as expected locally: `NERDSS_ENABLE_GTEST=ON` configure could not find
  `GTEST_LIBRARY`, `GTEST_INCLUDE_DIR`, or `GTEST_MAIN_LIBRARY`.
- Passed: `make serial -j4`.
- Passed: serial smoke runner with artifacts in
  `/tmp/nerdss-gtest-math-functions-smoke`.
- Not run locally: `clang-format --dry-run --Werror
  tests/gtest/math_functions_test.cpp`; `clang-format` is not on PATH in this
  shell.

Artifacts:
- `CMakeLists.txt`
- `.github/workflows/c-cpp.yml`
- `tests/gtest/math_functions_test.cpp`
- `docs/upgrade/unit_tests.md`
- This log entry.

Notes:
- This slice adds deterministic coverage for `MathFuncs::factorial`,
  `MathFuncs::gammln`, and `MathFuncs::gammFactorial`.
- The negative-input fatal path remains untested until the GTest migration has
  a settled pattern for death tests and error-message assertions.
- The Google Test target is not expected to build locally unless Google Test is
  installed or exposed through `CMAKE_PREFIX_PATH`.

## 2026-06-08: Rotation Math Google Test Conversion

Date: 2026-06-08

Branch: `codex/gtest-rotation-math`

Commit: Pending at validation log update

Workstream: Phase 7 Google Test migration and MathEngine preservation tests

Environment:
- OS: macOS, local Codex workspace
- Compiler: Apple clang via default CMake compiler
- GSL: 2.8 from Homebrew
- CMake: Available on PATH
- Make: GNU Make
- MPI: Local `mpicxx` wrapper remains blocked by missing
  `x86_64-apple-darwin13.4.0-clang++`

Commands:
```sh
git diff --check HEAD
tools/format_changed_files.sh --check --base HEAD
cmake -S . -B /tmp/nerdss-gtest-default
cmake --build /tmp/nerdss-gtest-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-gtest-on -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-gtest-rotation-math-smoke
```

Results:
- Passed: `git diff --check HEAD`.
- Passed: `tools/format_changed_files.sh --check --base HEAD`; the helper
  reported no tracked C/C++ files to format before staging this new test file.
- Passed: default CMake configure in `/tmp/nerdss-gtest-default`.
- Passed: default CMake build target `nerdss`.
- Failed as expected locally: `NERDSS_ENABLE_GTEST=ON` configure could not find
  `GTEST_LIBRARY`, `GTEST_INCLUDE_DIR`, or `GTEST_MAIN_LIBRARY`.
- Passed: `make serial -j4`.
- Passed: serial smoke runner with artifacts in
  `/tmp/nerdss-gtest-rotation-math-smoke`.
- Not run locally: `clang-format --dry-run --Werror
  tests/gtest/rotation_math_test.cpp`; `clang-format` is not on PATH in this
  shell.

Artifacts:
- `CMakeLists.txt`
- `.github/workflows/c-cpp.yml`
- `tests/gtest/rotation_math_test.cpp`
- `docs/upgrade/unit_tests.md`
- This log entry.

Notes:
- This slice adds deterministic coverage for quaternion algebra, quaternion
  vector rotation, and legacy Euler matrix rotation.
- The Google Test target is not expected to build locally unless Google Test is
  installed or exposed through `CMAKE_PREFIX_PATH`.

## 2026-06-08: First Coord/Vector Google Test Conversion

Date: 2026-06-08

Branch: `codex/gtest-vector-coord-tests`

Commit: Pending at validation log update

Workstream: Phase 7 Google Test migration

Environment:
- OS: macOS, local Codex workspace
- Compiler: Apple clang via default CMake compiler
- GSL: 2.8 from Homebrew
- CMake: Available on PATH
- Make: GNU Make
- MPI: Local `mpicxx` wrapper remains blocked by missing
  `x86_64-apple-darwin13.4.0-clang++`

Commands:
```sh
git diff --check HEAD
tools/format_changed_files.sh --check --base HEAD
cmake -S . -B /tmp/nerdss-gtest-default
cmake --build /tmp/nerdss-gtest-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-gtest-on -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-gtest-vector-coord-smoke
```

Results:
- Passed: `git diff --check HEAD`.
- Passed: `tools/format_changed_files.sh --check --base HEAD`; the helper
  reported no tracked C/C++ files to format before staging this new test file.
- Passed: default CMake configure in `/tmp/nerdss-gtest-default`.
- Passed: default CMake build target `nerdss`.
- Failed as expected locally: `NERDSS_ENABLE_GTEST=ON` configure could not find
  `GTEST_LIBRARY`, `GTEST_INCLUDE_DIR`, or `GTEST_MAIN_LIBRARY`.
- Passed: `make serial -j4`.
- Passed: serial smoke runner with artifacts in
  `/tmp/nerdss-gtest-vector-coord-smoke`.

Artifacts:
- `CMakeLists.txt`
- `.github/workflows/c-cpp.yml`
- `tests/gtest/vector_coord_test.cpp`
- `docs/upgrade/unit_tests.md`
- This log entry.

Notes:
- This slice converts only the prior low-risk Coord/Vector helper assertions.
- The Google Test target is not expected to build locally unless Google Test is
  installed or exposed through `CMAKE_PREFIX_PATH`.

## 2026-06-08: Minimal Google Test Harness

Date: 2026-06-08

Branch: `codex/gtest-minimal-harness`

Commit: Pending at validation log update

Workstream: Phase 7 Google Test migration

Environment:
- OS: macOS, local Codex workspace
- Compiler: Apple clang via default CMake compiler
- GSL: 2.8 from Homebrew
- CMake: Available on PATH
- Make: GNU Make
- MPI: Local `mpicxx` wrapper remains blocked by missing
  `x86_64-apple-darwin13.4.0-clang++`

Commands:
```sh
cmake --find-package -DNAME=GTest -DCOMPILER_ID=GNU -DLANGUAGE=CXX -DMODE=EXIST
cmake -S . -B /tmp/nerdss-gtest-on -DNERDSS_ENABLE_GTEST=ON
cmake -S . -B /tmp/nerdss-gtest-default
cmake --build /tmp/nerdss-gtest-default --target nerdss --parallel 4
tools/format_changed_files.sh --check --base HEAD
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-gtest-harness-smoke
git diff --check HEAD
```

Results:
- Local Google Test discovery reports `GTest not found`, so the new target was
  not built locally.
- `cmake -S . -B /tmp/nerdss-gtest-on -DNERDSS_ENABLE_GTEST=ON` fails locally
  at `find_package(GTest REQUIRED)` because Google Test is not installed.
- Default CMake configure/build and `make serial -j4` passed without Google
  Test because `NERDSS_ENABLE_GTEST` defaults to `OFF`.
- `tools/format_changed_files.sh --check --base HEAD` passed; no tracked C/C++
  files required formatting at that point.
- The direct smoke runner passed with `./bin/nerdss` and wrote artifacts under
  `/tmp/nerdss-gtest-harness-smoke`.
- `git diff --check HEAD` passed.
- GitHub Actions now installs `libgtest-dev`, configures
  `NERDSS_ENABLE_GTEST=ON`, builds `nerdss_gtest_smoke`, and runs the smoke
  binary.

Artifacts:
- `CMakeLists.txt`
- `.github/workflows/c-cpp.yml`
- `tests/gtest/smoke_test.cpp`
- `docs/upgrade/unit_tests.md`
- This log entry.

Notes:
- This slice intentionally adds only a harness smoke test. It does not convert
  existing CTest-era domain assertions.

## 2026-06-08: Google Test Migration Plan Update

Date: 2026-06-08

Branch: `codex/gtest-migration-plan`

Commit: Pending at validation log update

Workstream: Phase 7 test framework planning

Environment:
- OS: macOS, local Codex workspace
- Compiler: Not used
- GSL: Not used
- CMake: Not used
- Make: Not used
- MPI: Not used

Commands:
```sh
git fetch personal --prune
gh pr list --repo yingyue2030699/NERDSS_upgrade --state open --limit 100 --json number,title,headRefName,baseRefName,url
gh pr close 5 8 11 14 17 20 23 26 29 32 35 39 43 46
rg -n "CTest|ctest|custom CTest|Google Test|gtest|Google C\\+\\+" docs/refactor_implementation_plan.md docs/upgrade/unit_tests.md docs/upgrade/decision_log.md
git diff --check HEAD
```

Results:
- Closed the CTest-era validation-integration aggregation PRs as superseded:
  #5, #8, #11, #14, #17, #20, #23, #26, #29, #32, #35, #39, #43, and #46.
- Left focused implementation PRs open for later retargeting or Google
  Test-style assertion conversion.
- Updated Phase 7 planning docs to freeze custom CTest/std-library unit and
  integration test work and adopt Google Test.

Artifacts:
- `docs/refactor_implementation_plan.md`
- `docs/upgrade/unit_tests.md`
- `docs/upgrade/decision_log.md`
- This log entry.

Notes:
- This is a documentation and coordination update only. No test harness,
  executable behavior, parser behavior, or probability kernel code changed.
- Future validation commands should avoid relying on the superseded custom CTest
  unit executable until the Google Test harness lands.

## 2026-05-28: Agent A Phase 0 Environment Baseline

Date: 2026-05-28

Branch: `codex/upgrade-baseline-policy`

Commit: Pending at initial log creation

Workstream: Phase 0 repository baseline and branch policy

Environment:
- OS: macOS 14.5, build 23F79, arm64
- Compiler: Apple clang 16.0.0 via `g++` and `clang++`
- GSL: 2.8 from Homebrew, headers in `/opt/homebrew/Cellar/gsl/2.8/include`
- CMake: Not available on `PATH` in this environment
- Make: GNU Make 3.81
- MPI: MPICH/HYDRA 3.3.2 runtime present through Anaconda; `mpicxx --version`
  fails because the configured `x86_64-apple-darwin13.4.0-clang++` wrapper
  compiler is not available

Commands:
```sh
sw_vers
uname -m
g++ --version
clang++ --version
gsl-config --version
gsl-config --cflags
gsl-config --libs
cmake --version
make --version
mpicxx --version
mpirun --version
git status --short --branch
make serial
```

Results:
- Serial compiler and GSL dependencies are available.
- `make serial` compiles object files but fails when compiling
  `EXEs/nerdss.cpp` because the current baseline references
  `Parameters::bondedComplexWrite` and `write_bonded_complex_json` without
  visible declarations.
- CMake validation is blocked until CMake is installed or added to `PATH`.
- MPI build validation is blocked until the MPI compiler wrapper points to an
  available compiler.
- Repository worktree started clean on `codex/upgrade-baseline-policy`.

Artifacts:
- This log entry.

Notes:
- This Phase 0 branch is documentation and metadata only. It does not modify
  build logic or simulation source code.
- The serial build failure is recorded as a baseline blocker for follow-up by
  the build or smoke-runner workstream.
