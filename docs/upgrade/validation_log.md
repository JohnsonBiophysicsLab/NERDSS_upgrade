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

## 2026-06-09: Input Copy-Count Diagnostics

Date: 2026-06-09

Branch: `codex/parser-input-copy-count-diagnostics`

Commit: Pending at validation log update

Workstream: Parser diagnostics migration

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
cmake -S . -B /tmp/nerdss-parser-input-copy-count-default
cmake --build /tmp/nerdss-parser-input-copy-count-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-parser-input-copy-count-gtest -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-parser-input-copy-count-smoke
/Users/yueying/Documents/NERDSS\ upgrade/NERDSS_upgrade/bin/nerdss -f smoke.inp -s 123
```

Results:
- `git diff --check HEAD`: passed.
- `tools/format_changed_files.sh --check --base HEAD`: failed because the
  touched legacy file `src/parser/parse_input.cpp` has broad pre-existing
  clang-format drift; this slice did not reformat the whole file to keep the
  review focused.
- Default CMake configure: passed.
- Default CMake `nerdss` target build: passed.
- `cmake -S . -B /tmp/nerdss-parser-input-copy-count-gtest
  -DNERDSS_ENABLE_GTEST=ON`: failed as expected in this local workspace because
  Google Test is not installed/discoverable (`GTEST_LIBRARY`,
  `GTEST_INCLUDE_DIR`, and `GTEST_MAIN_LIBRARY` missing).
- `make serial -j4`: passed.
- Smoke runner with `--skip-build`: passed.
- Negative malformed-input probe from
  `/tmp/nerdss-parser-input-copy-count-negative`: passed, exiting with code `2`
  and `PARSER_ERROR[parse_input]` for input `A1`.

Artifacts:
- `src/parser/parse_input.cpp`
- This log entry.

Notes:
- This slice changes malformed molecule copy-count input from a printed
  message plus successful `exit(0)` to the shared parser diagnostics error
  path and input-error exit code.
- No CTest/custom test implementation was added; Google Test migration remains
  the selected path for future unit and integration tests.

## 2026-06-09: Molecule File Blank-Line Handling

Date: 2026-06-09

Branch: `codex/parser-molfile-empty-lines`

Commit: Pending at validation log update

Workstream: Parser diagnostics migration

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
cmake -S . -B /tmp/nerdss-parser-molfile-empty-default
cmake --build /tmp/nerdss-parser-molfile-empty-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-parser-molfile-empty-gtest -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-parser-molfile-empty-smoke
```

Results:
- `git diff --check HEAD`: passed.
- `tools/format_changed_files.sh --check --base HEAD`: failed because the
  touched legacy file `src/parser/parse_molFile.cpp` has broad pre-existing
  clang-format drift; this slice did not reformat the whole file to keep the
  review focused.
- Default CMake configure: passed.
- Default CMake `nerdss` target build: passed.
- `cmake -S . -B /tmp/nerdss-parser-molfile-empty-gtest
  -DNERDSS_ENABLE_GTEST=ON`: failed as expected in this local workspace because
  Google Test is not installed/discoverable (`GTEST_LIBRARY`,
  `GTEST_INCLUDE_DIR`, and `GTEST_MAIN_LIBRARY` missing).
- `make serial -j4`: passed.
- Smoke runner with `--skip-build`: passed.

Artifacts:
- `src/parser/parse_molFile.cpp`
- This log entry.

Notes:
- This slice changes whitespace-only molecule-file lines from undefined
  `line[0]` access to the same skip path used for full-line comments.
- Broader molecule parser syntax diagnostics remain queued for a separate
  focused slice.

## 2026-06-09: Parser File-Open Diagnostics

Date: 2026-06-09

Branch: `codex/parser-file-open-diagnostics`

Commit: Pending at validation log update

Workstream: Parser diagnostics migration

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
cmake -S . -B /tmp/nerdss-parser-file-open-default
cmake --build /tmp/nerdss-parser-file-open-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-parser-file-open-gtest -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-parser-file-open-smoke
```

Results:
- Passed: `git diff --check HEAD`.
- Failed: `tools/format_changed_files.sh --check --base HEAD` because the
  helper checks whole touched files and `parse_input.cpp` / `parse_molFile.cpp`
  still have broad pre-existing clang-format drift. The files were not
  auto-formatted in this focused diagnostics slice to avoid a large unrelated
  parser formatting diff.
- Passed: default CMake configure in `/tmp/nerdss-parser-file-open-default`.
- Passed: default CMake build target `nerdss`.
- Failed as expected locally: `NERDSS_ENABLE_GTEST=ON` configure could not find
  `GTEST_LIBRARY`, `GTEST_INCLUDE_DIR`, or `GTEST_MAIN_LIBRARY`.
- Passed: `make serial -j4`.
- Passed: serial smoke runner with artifacts in
  `/tmp/nerdss-parser-file-open-smoke`.
- Passed: `./bin/nerdss -f /tmp/nerdss-missing-input-file-does-not-exist.inp`
  exited nonzero and printed the new input filename diagnostic.

Artifacts:
- `src/parser/parse_input.cpp`
- `src/parser/parse_molFile.cpp`
- `src/parser/parse_input_for_a_restart_simulation.cpp`
- This log entry.

Notes:
- This slice preserves fatal behavior while adding filename/context to parser
  input, add-input, molecule-file, molecule-keyword, and restart-file open
  diagnostics.
- `parse_molecule_bngl.cpp` syntax diagnostics and empty-line `parse_molFile`
  crash behavior remain queued for separate slices.

## 2026-06-08: Parser Helper Diagnostics

Date: 2026-06-08

Branch: `codex/parser-helper-diagnostics`

Commit: Pending at validation log update

Workstream: Parser diagnostics migration

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
cmake -S . -B /tmp/nerdss-parser-helper-diagnostics-default
cmake --build /tmp/nerdss-parser-helper-diagnostics-default --target nerdss --parallel 4
cmake -S . -B /tmp/nerdss-parser-helper-diagnostics-gtest -DNERDSS_ENABLE_GTEST=ON
make serial -j4
python3 tools/run_smoke_tests.py --skip-build --executable ./bin/nerdss --artifact-dir /tmp/nerdss-parser-helper-diagnostics-smoke
```

Results:
- Passed: `git diff --check HEAD`.
- Fixed and passed: `tools/format_changed_files.sh --check --base HEAD`. The
  helper initially failed on an empty `unique_files` array under `set -u`, then
  successfully checked the changed parser diagnostics files after the guard fix
  and clang-format normalization.
- Passed: default CMake configure in
  `/tmp/nerdss-parser-helper-diagnostics-default`.
- Passed: default CMake build target `nerdss`.
- Failed as expected locally: `NERDSS_ENABLE_GTEST=ON` configure could not find
  `GTEST_LIBRARY`, `GTEST_INCLUDE_DIR`, or `GTEST_MAIN_LIBRARY`.
- Passed: direct `/tmp/nerdss_parser_diagnostics_probe` compile/run for valid
  parser-helper behavior and diagnostic formatting.
- Passed: direct `/tmp/nerdss_parser_diagnostics_failure_probe` compile/run;
  invalid `read_boolean` input printed structured parser diagnostics and exited
  with code 2.
- Passed: `make serial -j4`.
- Passed: serial smoke runner with artifacts in
  `/tmp/nerdss-parser-helper-diagnostics-smoke`.

Artifacts:
- `include/parser/parser_diagnostics.hpp`
- `src/parser/parser_diagnostics.cpp`
- `src/parser/read_boolean.cpp`
- `src/parser/parse_input_array.cpp`
- `tests/gtest/parser_helpers_test.cpp`
- `CMakeLists.txt`
- `docs/upgrade/unit_tests.md`
- `tools/format_changed_files.sh`
- This log entry.

Notes:
- This slice adds a lightweight parser diagnostic formatter/fatal helper and
  applies it to invalid `read_boolean` and `parse_input_array` parser-helper
  paths.
- Valid parser-helper behavior remains covered by `nerdss_gtest_parser_helpers`.
- File-open parser diagnostics remain queued as the next focused slice.

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
c++ -std=c++11 -Iinclude -c src/parser/remove_comment.cpp -o /tmp/nerdss-remove-comment.o
c++ -std=c++11 -Iinclude -c src/parser/read_boolean.cpp -o /tmp/nerdss-read-boolean.o
c++ -std=c++11 -Iinclude -c src/parser/parse_input_array.cpp -o /tmp/nerdss-parse-input-array.o
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
- Passed: standalone object compile probes for `remove_comment.cpp`,
  `read_boolean.cpp`, and `parse_input_array.cpp`.
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
