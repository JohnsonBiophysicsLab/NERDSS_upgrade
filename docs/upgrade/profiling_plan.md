# NERDSS Profiling Plan (Per-Function Timing)

Companion to `profile_report.md` (which records measured results). This document
is the *method*: how to get evidence-based, per-function timing for NERDSS on
both macOS (dev) and Linux (authoritative), without changing simulation behavior.

## Context

Phase 6 profiling must identify expensive components from measurement, not
intuition. The repo already had a command guide and empty result template; this
plan turns that into a repeatable workflow backed by scripts, and its first pass
is recorded in `profile_report.md` ("Measured Results — macOS Serial").

## Tooling (reuse these; don't rebuild)

- **`tools/run_profile.sh`** — the driver. Builds (optional), copies a case into
  an isolated `profile-runs/` dir, caps `nItr` in that copy only, measures CPU
  time, runs the chosen profiler, and prints per-function **seconds**. Modes:
  - `--profiler instruments` (macOS, default) — self + inclusive seconds.
  - `--profiler pprof` (macOS) — real call graph via `/usr/bin/sample` →
    `google/instrumentsToPprof` → `go tool pprof`.
  - `--profiler perf` (Linux) — `perf record -g --call-graph dwarf`, per-function
    seconds, and (with `perf_to_profile`) a pprof profile for the same graph tools.
  - `--profiler gperftools` (Linux) — libprofiler via preload.
- **`tools/profile_aggregate.py`** — aggregates an `xctrace` time-profile export
  into self/inclusive time; `--total-seconds` scales sample shares to real CPU time.
- **`tools/profile_commands.sh`** — prints per-tool command blocks (no execution).
- Build: `make serial` (`-O3`); `make serial profile` only for `-pg`/gperftools linkage.

## Strategy

Per-function timing comes from three complementary layers, in order:

1. **Sampling on the optimized build** — the source of truth for which functions
   cost time (no instrumentation skew). Instruments on macOS; `perf` on Linux.
2. **Call graph** — who calls the hot function, and how often. On macOS use the
   `sample`→`instrumentsToPprof`→`pprof` route (gperftools' unwinder is broken on
   arm64 macOS). On Linux use `perf` (+ `perf_to_profile` for the pprof workflow).
3. **Targeted in-code timing** — only when a hot function is inlined/too short for
   the sampler: `TRACE()` (`include/tracing.hpp`, debug builds) or a `std::chrono`
   scope like the main-loop `durationList` in `EXEs/nerdss.cpp`. Revert before merge.

## Cases and knobs

- Small: `homoTrimer` / `parmTri6.inp`, `nItr=2000` (~4 s) — fast sanity, known
  validation artifacts.
- Large: `clathrin` / `parms_clath_kon1uM.inp` — raise `nItr` (~100k for
  Instruments, ~400k for `pprof`/`sample`) so the run is long enough to sample.
- Fixed seed **12345** for every run. `nItr` is capped in the isolated copy only.

## Steps

1. **Baseline** — `make clean && make serial`; optionally
   `python3 benchmarks/run_benchmarks.py` for wall time / peak RSS.
2. **Sample (per-function seconds)** — small then large:
   `tools/run_profile.sh --profiler instruments|perf --case <dir> --input <file> --nitr <n>`.
   Take the top ~15 by self time; note high-inclusive callers.
3. **Call graph** — `--profiler pprof` (macOS) or `--profiler perf` (Linux); then
   `go tool pprof -peek=<fn>` to attribute allocation/math time to callers, and
   `-png/-svg/-http` (needs `brew install graphviz`) for images.
4. **Map to categories** — reaction checks, propagation, overlap/boundary, IO,
   parsing/setup, MPI (per `profile_report.md`).
5. **Record** — fill `profile_report.md`: metadata, artifacts, ranked hotspots
   (with seconds + source lines), call-graph findings. Keep raw traces out of git
   (`profile-runs/` is gitignored).

## Guardrails

- Profile the optimized build; `-pg` only for gprof call graphs.
- Fixed seed 12345; identical input copies per run.
- No hotspot claim without a supporting profiler artifact.
- Link each optimization idea to a benchmark case **and** a fixed-seed regression.
- Preserve scientific behavior — profiling is measurement only.

## Platform notes

- **macOS/arm64:** Instruments is the source of truth; gperftools does **not**
  unwind (mis-attributes ~48% to `libc++`). Use the `instrumentsToPprof` route
  for call graphs (`go install github.com/google/instrumentsToPprof@latest`).
- **Linux:** `perf` for hotspots, gperftools for per-function/MPI. Install
  `perf_to_profile` (google/perf_data_converter) to reuse the `pprof` graph tools;
  `perf record` may need `sysctl kernel.perf_event_paranoid=-1` or sudo.

## Status / next

First pass (macOS serial, small + large) is recorded in `profile_report.md`,
including the finding that ~100% of malloc churn traces to
`check_bimolecular_reactions.cpp:111`. Open: Linux `perf`/gperftools corroboration
and MPI per-rank profiling.
