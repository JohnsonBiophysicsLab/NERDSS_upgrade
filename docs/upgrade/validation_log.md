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
