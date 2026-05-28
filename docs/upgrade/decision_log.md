# Upgrade Decision Log

This log records coordination decisions for the NERDSS upgrade. Add entries in
reverse chronological order and link follow-up PRs or validation artifacts when
they exist.

## 2026-05-28: Phase 0 Repository Baseline

Status: Accepted

Owners: Agent A, upgrade coordination team

Decision:
- Upgrade branches use the `codex/upgrade-*` naming pattern.
- Refactor PRs must stay scoped to one workstream or vertical slice.
- PRs that touch executable behavior, parser behavior, output formats, RNG call
  order, MPI partitioning, or numerical algorithms must include a validation log
  entry before review.
- Documentation-only and metadata-only PRs may validate with link checks, command
  examples, and repository status checks.

Rationale:
- Multiple agents are working concurrently, so predictable branch names and small
  PRs reduce merge and review conflicts.
- Baseline records make it possible to distinguish mechanical refactors from
  scientific behavior changes.

Consequences:
- Contributors should branch from the current integration branch, then create
  focused branches such as `codex/upgrade-smoke-runner` or
  `codex/upgrade-parser-errors`.
- Generated documentation under `docs/html/` should not be updated incidentally
  in refactor PRs.
- Any behavior-changing proposal needs an explicit scientific review record
  before implementation.

## 2026-05-28: Minimum C++ Standard During Upgrade Baseline

Status: Accepted for Phase 0

Owners: Agent A, upgrade coordination team

Decision:
- Keep C++11 as the minimum supported language standard for the Phase 0 baseline.
- Do not introduce C++14, C++17, or newer syntax in upgrade branches until a
  dedicated modernization decision updates the build files and validation plan.
- Treat the README badge that says C++ >=14 as stale until the team resolves it
  in a documentation or build-system PR.

Rationale:
- `CMakeLists.txt` sets `CMAKE_CXX_STANDARD 11`.
- `Makefile` currently compiles with `-std=c++0x`.
- `docs/devel.md` explicitly asks contributors to write C++11 compatible with
  GCC 4.9.0.
- Preserving the current language floor avoids coupling broad modernization to
  baseline and validation work.

Consequences:
- Phase 0 and dependent refactor work should assume C++11 unless a later entry
  supersedes this decision.
- Modernization proposals must include the minimum compiler matrix, build-system
  edits, and regression validation scope.
