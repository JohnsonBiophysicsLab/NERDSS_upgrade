# Upgrade Contributor Guidance

This page documents the Phase 0 baseline for upgrade work. It is intentionally
conservative: preserve current simulation behavior first, then modernize only
behind validation.

## Environment Baseline

Required for serial development:
- macOS or Linux.
- A C++ compiler with C++11 support.
- GSL available through `gsl-config`.
- GNU Make.

Optional for additional validation:
- CMake for the serial CMake build.
- MPI compiler and runtime for the parallel build.

Current minimum language standard:
- C++11 is the supported minimum for upgrade branches.
- Do not use C++14, C++17, or newer features until a later decision updates the
  build configuration and validation plan.

Reference environment recorded on 2026-05-28:

| Component | Version or status |
| --- | --- |
| OS | macOS 14.5, build 23F79, arm64 |
| Compiler | Apple clang 16.0.0 through `g++` and `clang++` |
| GSL | 2.8 through Homebrew |
| CMake | Not available on `PATH` |
| Make | GNU Make 3.81 |
| MPI | MPICH/HYDRA 3.3.2 runtime present; `mpicxx` wrapper misconfigured |

## Serial Build

From a clean checkout or isolated worktree:

```sh
gsl-config --version
make serial
./bin/nerdss -h
```

Expected result:
- `make serial` creates `bin/nerdss`.
- If GSL is not installed or `gsl-config` is not on `PATH`, the Makefile stops
  before compiling.

Known Phase 0 baseline status:
- On the 2026-05-28 Agent A macOS baseline, `make serial` reached
  `EXEs/nerdss.cpp` and failed because existing references to
  `Parameters::bondedComplexWrite` and `write_bonded_complex_json` were not
  declared in the visible headers. This branch records that blocker but does not
  change build or simulation code.

## Optional Builds

CMake serial build:

```sh
cmake -S . -B build
cmake --build build
```

MPI build:

```sh
make mpi
```

Only report these as validated when the relevant tools are installed and the
commands complete successfully.

## Branch Policy

Use one branch per focused upgrade workstream:
- Branch names use `codex/upgrade-*`.
- Keep branches scoped to one phase, agent, or vertical slice.
- Prefer isolated git worktrees when multiple agents share a repository clone.
- Do not mix generated documentation, formatting sweeps, and behavior changes in
  the same PR.

Examples:
- `codex/upgrade-baseline-policy`
- `codex/upgrade-smoke-runner`
- `codex/upgrade-regression-harness`
- `codex/upgrade-style-tooling`

## PR Policy

Each PR should include:
- Scope summary.
- Files or subsystems changed.
- Build and validation commands run.
- Explicit notes for skipped validation.
- Any scientific behavior risk, even when the expected risk is none.

Additional requirements:
- Behavior-preserving refactors must include baseline validation for the touched
  workflow before review.
- Scientific behavior changes require a separate reviewed decision before code
  changes land.
- Keep old input and output formats compatible unless a documented deprecation
  decision has been accepted.
