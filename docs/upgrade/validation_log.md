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
- Make:
- MPI:

Commands:

Results:

Artifacts:

Notes:

## 2026-06-11: Google Test Scope Warning

Date: 2026-06-11

Branch: `codex/cleanup-gtest-integration`

Commit: Pending at log update

Workstream: Google Test cleanup

Environment:
- OS: macOS local Codex workspace
- Compiler: Apple clang via `g++`
- GSL: 2.8 from Homebrew
- Make: GNU Make
- MPI: Not exercised for this warning update

Commands:
```sh
make
make serial
make unittest TEST_VERBOSE=1
```

Results:
- The Google Test cleanup is intentionally Makefile-only.
- Agents must not add or extend CMake/CTest Google Test integration in this
  workstream.
- `make unittest` remains the explicit Google Test entry point.
- CMake-based Google Test harnesses, smoke targets, or CI commands are out of
  scope until a future reviewed request reopens that work.

Artifacts:
- `Makefile`
- `docs/refactor_implementation_plan.md`
- `docs/upgrade/unit_tests.md`
- This log entry.

Notes:
- This warning is intentionally duplicated in the unit-test guidance and the
  implementation plan so future branches do not reintroduce CMake work while
  cleaning up Google Test.

## 2026-05-28: Agent A Phase 0 Environment Baseline

Date: 2026-05-28

Branch: `codex/upgrade-baseline-policy`

Commit: Pending at initial log creation

Workstream: Phase 0 repository baseline and branch policy

Environment:
- OS: macOS 14.5, build 23F79, arm64
- Compiler: Apple clang 16.0.0 via `g++` and `clang++`
- GSL: 2.8 from Homebrew, headers in `/opt/homebrew/Cellar/gsl/2.8/include`
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
