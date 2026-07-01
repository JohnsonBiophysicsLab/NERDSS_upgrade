# NERDSS Profiling And Hotspot Report

Owner: Agent P

Status: profiler command guide + first measured pass. macOS serial
Instruments results for the small and large cases are recorded under
"Measured Results — macOS Serial (2026-07-01)". Linux `perf`/gperftools
corroboration and MPI profiling still pending.

## Purpose

Phase 6 profiling should identify expensive NERDSS components from
measurements, not intuition. Use this report to record the build, case, profiler
output, ranked hotspots, and the validation guard needed before any
optimization PR touches simulation behavior.

Target attribution categories:

- Reaction checks: bimolecular, unimolecular, create/destruct, transmission,
  implicit-lipid reactions.
- Propagation: translational and rotational trajectory updates.
- Overlap and boundary checks: box, sphere, compartment, membrane, implicit
  lipid, and reflection paths.
- IO: trajectory, restart, PDB/CIF, observables, histograms, transition matrix,
  JSON complex output.
- Parsing and setup: input parsing, molecule/reaction construction, initial
  coordinates, probability table setup.
- MPI communication: rank exchange, neighbor buffers, merge/write paths, and
  per-rank imbalance.

## Current Profile Build Intent

The Makefile advertises profiling through:

```sh
make serial profile
make mpi profile
```

Current behavior to account for in reports:

- `profile` adds `-pg`, which supports `gprof` style call-graph output on
  platforms where `gprof` is available.
- `profile` defines `ENABLE_PROFILING`.
- The MPI executable compiles gperftools hooks under `ENABLE_PROFILING` and
  writes `profile_output_<rank>.prof` files when gperftools headers/libraries
  are installed.
- If a combined MPI goal such as `make mpi profile` selects `g++` instead of
  the MPI wrapper, pass `CC=mpicxx` on the command line and record that in the
  report. This documents the current Makefile behavior without changing build
  logic in this planning PR.
- The serial executable currently has gperftools hooks commented out, so its
  first supported path is `gprof`/external sampling unless the binary is linked
  with `libprofiler` and `CPUPROFILE` is used for whole-process profiling.
- Optimized builds use `-O3`; if symbols are hard to interpret, repeat with a
  local profiling build that adds debug symbols or lowers optimization, and
  record the exact flags. Do not commit source changes just to profile.

## Recommended Cases

Replace or extend this table once Agent O publishes the benchmark harness.

| Size | Case | Command directory | Input | Why profile it |
| --- | --- | --- | --- | --- |
| Small correctness-adjacent | `homoTrimer` | `sample_inputs/VALIDATE_SUITE/homoTrimer` | `parmTri6.inp` | Exercises association/dissociation and restart/output paths with known validation artifacts. |
| Small correctness-adjacent | `hetTrimer` | `sample_inputs/VALIDATE_SUITE/hetTrimer` | `parm_autodiff_hetTri.inp` | Adds multiple molecule types and reaction combinations. |
| Boundary/geometry | `sphere` | `sample_inputs/sphere` | `parms_sphere.inp` | Stresses non-box boundary handling. |
| Representative large | TBD by Agent O | `benchmarks/...` | TBD | Use for optimization decisions; keep out of routine CI if runtime is high. |
| MPI representative | TBD by Agent O | `benchmarks/...` or validation MPI case | TBD | Required before MPI-sensitive refactors. |

## Command Helper

Use `tools/profile_commands.sh` to print platform-specific command blocks:

```sh
tools/profile_commands.sh --mode serial --profiler all \
  --case sample_inputs/VALIDATE_SUITE/homoTrimer \
  --input parmTri6.inp --seed 12345

tools/profile_commands.sh --mode mpi --profiler gperftools \
  --case sample_inputs/VALIDATE_SUITE/homoTrimer \
  --input parmTri6.inp --ranks 4 --seed 12345
```

The helper prints commands only; it does not build or run NERDSS.

## Linux Profiling Guidance

### `gprof`

Use as a first pass when the `make ... profile` binary emits `gmon.out`.

```sh
make clean
make serial profile
mkdir -p profile-runs/linux-gprof-homoTrimer
cp sample_inputs/VALIDATE_SUITE/homoTrimer/*.inp sample_inputs/VALIDATE_SUITE/homoTrimer/*.mol profile-runs/linux-gprof-homoTrimer/
cd profile-runs/linux-gprof-homoTrimer
/usr/bin/time -p ../../bin/nerdss -f parmTri6.inp -s 12345
gprof ../../bin/nerdss gmon.out > gprof.txt
```

For MPI, set `GMON_OUT_PREFIX` so ranks do not overwrite one another:

```sh
make clean
make mpi profile CC=mpicxx
mkdir -p profile-runs/linux-gprof-mpi
cp sample_inputs/VALIDATE_SUITE/homoTrimer/*.inp sample_inputs/VALIDATE_SUITE/homoTrimer/*.mol profile-runs/linux-gprof-mpi/
cd profile-runs/linux-gprof-mpi
GMON_OUT_PREFIX=gmon mpirun -np 4 ../../bin/nerdss_mpi -f parmTri6.inp -s 12345
for f in gmon.*; do gprof ../../bin/nerdss_mpi "$f" > "$f.txt"; done
```

### `perf`

Use `perf` on Linux for sampling optimized binaries without relying on
instrumentation. Prefer call stacks with frame pointers or DWARF call graph
capture.

```sh
make clean
make serial
mkdir -p profile-runs/linux-perf-homoTrimer
cp sample_inputs/VALIDATE_SUITE/homoTrimer/*.inp sample_inputs/VALIDATE_SUITE/homoTrimer/*.mol profile-runs/linux-perf-homoTrimer/
cd profile-runs/linux-perf-homoTrimer
perf record -F 99 -g --call-graph dwarf -o perf.data -- ../../bin/nerdss -f parmTri6.inp -s 12345
perf report -i perf.data
perf script -i perf.data > perf.script.txt
```

For MPI, profile one rank first to keep reports tractable. If profiling all
ranks, write separate output files per rank and record rank imbalance.

`tools/run_profile.sh --profiler perf` wraps the above: it records with perf,
prints per-function overhead scaled to real seconds, and — when Google's
`perf_to_profile` (perf_data_converter) is installed — converts `perf.data` to a
pprof profile so the same `go tool pprof -peek / -png / -http` workflow used on
macOS applies. Without the converter, use `perf report` (interactive call graph)
or FlameGraph (`perf script | stackcollapse-perf.pl | flamegraph.pl > flame.svg`).

```sh
tools/run_profile.sh --profiler perf --build \
  --case sample_inputs/VALIDATE_SUITE/clathrin --input parms_clath_kon1uM.inp --nitr 100000
# then, if converted:
go tool pprof -png -nodefraction=0.01 -output cg.png profile-runs/perf-clathrin/cpu.pb.gz
```

### gperftools

Use when `gperftools` is installed and symbols are readable. MPI profile builds
currently start/stop gperftools internally and emit one `.prof` file per rank.

```sh
make clean
make mpi profile CC=mpicxx
mkdir -p profile-runs/linux-gperftools-mpi
cp sample_inputs/VALIDATE_SUITE/homoTrimer/*.inp sample_inputs/VALIDATE_SUITE/homoTrimer/*.mol profile-runs/linux-gperftools-mpi/
cd profile-runs/linux-gperftools-mpi
mpirun -np 4 ../../bin/nerdss_mpi -f parmTri6.inp -s 12345
for f in profile_output_*.prof; do pprof --text ../../bin/nerdss_mpi "$f" > "$f.txt"; done
```

For serial, use `CPUPROFILE` only if the profile build links `libprofiler`:

```sh
make clean
make serial profile
mkdir -p profile-runs/linux-gperftools-serial
cd profile-runs/linux-gperftools-serial
CPUPROFILE=cpu.prof ../../bin/nerdss -f parmTri6.inp -s 12345
pprof --text ../../bin/nerdss cpu.prof > pprof.txt
```

## macOS Profiling Guidance

### Instruments / `xctrace`

Use Instruments Time Profiler for optimized serial and MPI runs on macOS.

```sh
make clean
make serial
mkdir -p profile-runs/macos-instruments-homoTrimer
cp sample_inputs/VALIDATE_SUITE/homoTrimer/*.inp sample_inputs/VALIDATE_SUITE/homoTrimer/*.mol profile-runs/macos-instruments-homoTrimer/
cd profile-runs/macos-instruments-homoTrimer
xcrun xctrace record --template "Time Profiler" \
  --output TimeProfiler.trace \
  --launch -- ../../bin/nerdss -f parmTri6.inp -s 12345
xcrun xctrace export --input TimeProfiler.trace --xpath '/trace-toc/run/data/table' > xctrace-table.xml
```

For GUI inspection:

```sh
open profile-runs/macos-instruments-homoTrimer/TimeProfiler.trace
```

For MPI, prefer a short run and compare per-rank behavior. Instruments may
attach most cleanly to one rank at a time; record the exact `mpirun` command and
rank selection method used.

### `sample`

Use `sample` for a lightweight sanity pass when Instruments setup is too heavy.
Run NERDSS in one terminal, then sample the process from another:

```sh
pgrep -fl nerdss
sample <PID> 10 -file sample.txt
```

### pprof call graph on macOS (`sample` → instrumentsToPprof)

To get a real pprof call graph on macOS, sample the process and convert with
[`google/instrumentsToPprof`](https://github.com/google/instrumentsToPprof). This
uses the same reliable macOS unwinder as `sample`/Instruments, so stacks are
correctly resolved to source lines (unlike gperftools here).

```sh
go install github.com/google/instrumentsToPprof@latest      # provides the converter
# run NERDSS long enough to sample (raise nItr), grab its pid, then:
sample <PID> 12 -file sample.txt
instrumentsToPprof --format=sample -output cpu.pb.gz sample.txt
go tool pprof -top cpu.pb.gz                 # flat
go tool pprof -tree cpu.pb.gz                # call tree
go tool pprof -peek=free_tiny cpu.pb.gz      # callers of a function
go tool pprof -http=: cpu.pb.gz              # web graph + flame graph (needs graphviz)
```

Call-graph **image** (needs `brew install graphviz`):

```sh
go tool pprof -png -nodefraction=0.01 -output callgraph.png cpu.pb.gz  # PNG (pruned)
go tool pprof -svg -output callgraph.svg cpu.pb.gz                     # SVG (zoomable)
go tool pprof -png -focus=check_bimolecular_reactions -output cg.png cpu.pb.gz  # one subtree
```

The same `-png/-svg/-http` commands work on Linux for the pprof profile produced
by the `perf` path (see Linux `perf`).

`tools/run_profile.sh --profiler pprof ...` automates launch + sample + convert.

### gperftools On macOS

Not viable on Apple Silicon: gperftools' CPU profiler fails to unwind on arm64
(observed mis-attributing ~48% of samples to `libc++` and never resolving NERDSS
symbols, with and without frame pointers). Use the `instrumentsToPprof` path
above on macOS; keep gperftools for Linux. Treat Instruments as the default macOS
source of truth.

## Result Template

### Run Metadata

| Field | Value |
| --- | --- |
| Date | TBD |
| Agent | Agent P |
| Git commit | TBD |
| Branch | `codex/upgrade-profile-plan` or successor branch |
| OS / kernel | TBD |
| CPU / memory | TBD |
| Compiler | TBD |
| GSL version | TBD |
| MPI implementation | TBD or N/A |
| Build command | TBD |
| Case | TBD |
| Input command | TBD |
| Seed | TBD |
| Wall time | TBD |
| Peak RSS | TBD |
| Output artifact directory | TBD |

### Profiler Artifacts

| Tool | Artifact | Notes |
| --- | --- | --- |
| `gprof` | TBD | Include flat profile and call graph summary. |
| `perf` | TBD | Include top symbols, call graph mode, and sampling frequency. |
| Instruments | TBD | Include Time Profiler trace and exported table if available. |
| gperftools | TBD | Include one report per rank for MPI. |

### Ranked Hotspots

| Rank | Component | Evidence | Percent / samples | File or symbol examples | Optimization candidate | Validation guard |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | TBD | TBD | TBD | TBD | TBD | Regression case + benchmark case |
| 2 | TBD | TBD | TBD | TBD | TBD | Regression case + benchmark case |
| 3 | TBD | TBD | TBD | TBD | TBD | Regression case + benchmark case |

### MPI Notes

| Rank group | Observation | Evidence | Follow-up |
| --- | --- | --- | --- |
| Rank 0 | TBD | TBD | TBD |
| Worker ranks | TBD | TBD | TBD |
| Imbalance | TBD | TBD | TBD |

### Optimization Experiments

| Experiment | Expected benefit | Behavior risk | Required validation | Benchmark acceptance |
| --- | --- | --- | --- | --- |
| Reserve or reuse frequently allocated containers in measured hotspot | Lower allocation overhead | Low if ownership and order are unchanged | Fixed-seed regression for affected case | No unexplained slowdown; output unchanged |
| Isolate repeated overlap/boundary checks behind helper with same call order | Reduce duplicated branch work | Medium; geometry-sensitive | Boundary-specific regression plus benchmark | Equal outputs within tolerance; speedup on boundary case |
| Improve data locality in reaction candidate traversal | Reduce cache misses | High; RNG and reaction order sensitive | Fixed-seed and stochastic checks | Speedup with no reaction count drift |

## Measured Results — macOS Serial (2026-07-01)

First populated pass. Function-level sampling with Instruments Time Profiler on
an optimized `make serial` build. Both cases use fixed seed `12345`. Percentages
are share of total on-CPU samples (self / exclusive time, leaf frame).

### Run Metadata

| Field | Small case | Large case |
| --- | --- | --- |
| Date | 2026-07-01 | 2026-07-01 |
| Agent | Agent P | Agent P |
| Git commit | `bb30f57` | `bb30f57` |
| OS / kernel | macOS 15.7.2 (24G325) | macOS 15.7.2 |
| CPU / memory | Apple Silicon MacBook Pro (arm64e, P/E cores) | same |
| Compiler | Apple clang 17.0.0 (`g++` → clang), `-O3` | same |
| GSL version | libgsl.28 | libgsl.28 |
| MPI implementation | N/A (serial) | N/A (serial) |
| Build command | `make clean && make serial` | same |
| Case | `homoTrimer` | `clathrin` |
| Input command | `nerdss -f parmTri6.inp -s 12345` | `nerdss -f parms_clath_kon1uM.inp -s 12345` |
| nItr (profiling) | 2000 (benchmark "short") | 100000 (raised from 2000 so the run lasts long enough to sample) |
| Seed | 12345 | 12345 |
| Wall time | ~3.9 s | ~4.9 s |
| Peak RSS | 6.57 MB | 3.95 MB |
| Samples collected | 1234 (@1 ms) | 1147 (@1 ms) |
| Output artifact directory | `profile-runs/macos-instruments-homoTrimer/` | `profile-runs/macos-instruments-clathrin/` |

Self time by binary (small / large): `nerdss` 44% / 51%, `libsystem_m` (transcendental
math) 19% / 19%, `libsystem_malloc` 17% / 18%, `libgsl` 3% / 7%.

### Profiler Artifacts

| Tool | Artifact | Notes |
| --- | --- | --- |
| Instruments (Time Profiler) | `TimeProfiler.trace`, `time-profile.xml` in each run dir | `xctrace record --template "Time Profiler"`; exported table aggregated to per-function self/inclusive **seconds** by `tools/profile_aggregate.py`. |
| pprof (macOS) | `sample.txt`, `cpu.pb.gz` in each pprof run dir | `/usr/bin/sample` → `instrumentsToPprof --format=sample` → `go tool pprof`. |
| gperftools | — | Attempted; unusable on arm64 macOS (see Notes). |

Both are wrapped by **`tools/run_profile.sh`** (see Reproduction Commands).

### Reproduction Commands

The whole workflow is wrapped by `tools/run_profile.sh`, which builds (optional),
copies the case into an isolated `profile-runs/` dir, caps `nItr` in that copy
only, measures CPU time, runs the profiler, and prints per-function **seconds**.

```sh
# Per-function SELF/INCLUSIVE seconds (Instruments) — small and large
tools/run_profile.sh --case sample_inputs/VALIDATE_SUITE/homoTrimer \
  --input parmTri6.inp --nitr 2000 --build
tools/run_profile.sh --case sample_inputs/VALIDATE_SUITE/clathrin \
  --input parms_clath_kon1uM.inp --nitr 100000

# Call graph (macOS): sample -> instrumentsToPprof -> go tool pprof
#   (needs a longer run so `sample` gets enough on-CPU time)
go install github.com/google/instrumentsToPprof@latest   # one-time
tools/run_profile.sh --profiler pprof --sample-seconds 12 \
  --case sample_inputs/VALIDATE_SUITE/clathrin \
  --input parms_clath_kon1uM.inp --nitr 400000
go tool pprof -peek=free_tiny profile-runs/pprof-clathrin/cpu.pb.gz   # who allocates
```

Underlying primitives (what the script runs): `make serial`; per-case
`xcrun xctrace record --template "Time Profiler"` + `xctrace export
--xpath '.../table[@schema="time-profile"]'` aggregated by
`tools/profile_aggregate.py --total-seconds <measured CPU s>`; and for the graph,
`/usr/bin/sample <pid> <secs>` + `instrumentsToPprof --format=sample`. Print-only
per-tool command blocks are also available via `tools/profile_commands.sh`.

### Ranked Hotspots

Seconds are scaled to measured CPU time (small case ≈ 3.61 s user+sys; percentages
are share of on-CPU samples). "self" = exclusive, "incl" = inclusive of callees.

| Rank | Component | Evidence — small case (self % / self s ; incl) | File or symbol examples (source line) | Optimization candidate | Validation guard |
| --- | --- | --- | --- | --- | --- |
| 1 | Reaction checks — bimolecular candidate eval | ~12% self / ~0.43 s; **incl 43% (1.55 s)** | `check_bimolecular_reactions.cpp` (self 0.18 s), `determine_3D_bimolecular_reaction_probability` (self 0.14 s, incl 0.77 s), `find_which_reaction.cpp:109` (0.11 s), `pirr_pfree_ratio_psF.cpp:70` (incl 0.27 s) | Hoist repeated work / reduce per-candidate branching; cache reaction lookups | Fixed-seed regression per reaction family + benchmark case |
| 2 | Transcendental math in probability/sampling | ~12% self / ~0.42 s (`libsystem_m`) | `exp` (0.23 s), `log` (0.09 s), `__sincos_stret` (0.08 s), `erfc`, `pow` — called from probability + `GaussV` | Precompute/table probability terms; avoid redundant `exp/erfc` per step | Numerical-tolerance regression (probabilities unchanged) |
| 3 | Allocation churn | ~17% self / ~0.61 s (`libsystem_malloc`) | `free_tiny` (0.16 s), `tiny_malloc_should_clear`, `_szone_free` — **100% of `free_tiny` traced to `check_bimolecular_reactions.cpp:111`** | Reserve/reuse the container allocated at check_bimolecular_reactions.cpp:111 (matches Optimization Experiment #1) | Fixed-seed regression; output byte-identical |
| 4 | RNG draws | small ~5% / large ~8% | `GaussV` (`rand_gsl.cpp:130`), `mt_get_double`, `gsl_rng_uniform` | Reduce draw count / batch normals; only if reaction+move order preserved | Fixed-seed + stochastic-count checks (no draw-order drift) |
| 5 | Propagation & boundary geometry | small ~9% / large ~13% | `sweep_separation_complex_rot_box.cpp:247` (incl 0.50 s), `Complex::propagate`, `Quat::rotate` (`class_Quat.cpp:42`), `Coord::Coord` (`class_Coord.cpp:28`), `get_distance.cpp`, `reflect_traj_complex_rad_rot_box` | Improve data locality in complex traversal; grows with complex size | Boundary-specific regression + benchmark (geometry-sensitive) |
| 6 | Spatial binning (large/longer runs) | large ~3–4% self | `SimulVolume::update_memberMolLists` (`class_SimulVolume.cpp:255`) | Only rebin changed cells; avoid full membership rebuild | Fixed-seed regression on a case that grows complexes |

### Call-graph attribution (pprof via instrumentsToPprof)

gperftools' own CPU profiler does **not** unwind on arm64 macOS (it mis-attributes
~48% of samples to `libc++` and never resolves NERDSS symbols, with or without
frame pointers). The working macOS route is to `/usr/bin/sample` the process and
convert with [`google/instrumentsToPprof`](https://github.com/google/instrumentsToPprof)
(`--format=sample`), then render with `go tool pprof`. This yields correctly
unwound, source-line-resolved stacks. Key `-peek` findings:

- **`free_tiny` (malloc churn): 100% called from `check_bimolecular_reactions.cpp:111`** —
  a per-candidate container allocated/freed inside the hot reaction loop. Prime
  target for Optimization Experiment #1 (reserve/reuse).
- `tiny_malloc_should_clear` ← `szone_malloc_should_clear` ← `operator new` on the
  same reaction path; smaller `free_tiny` contributions from
  `sweep_separation_complex_rot_box.cpp:247`.

Reproduce: `tools/run_profile.sh --profiler pprof --case <dir> --input <file> --nitr <large>`,
then `go tool pprof -peek=<function> profile-runs/pprof-<case>/cpu.pb.gz`.

### Notes / Limitations

- Seconds are sampling estimates scaled to a plain (un-profiled) run's measured
  CPU time (`/usr/bin/time`), since raw sample weight under-counts true CPU time.
- The large case was run at `nItr=100000` (not the benchmark `2000`) because at
  2000 it finishes in ~0.13 s — too few samples. The small case at `nItr=2000`
  runs ~3.9 s. pprof-mode runs use a larger `nItr` so `/usr/bin/sample` gets
  ≥ ~9 s of on-CPU samples. All runs are single-threaded serial, fixed seed 12345.
- Attribution is leaf/self time from 1 ms sampling; short/inlined helpers may be
  folded into callers. A few samples resolve to unmapped stubs
  (`0x1a37b60f4` ≈ 3%, `<deduplicated_symbol>`).
- The two cases agree on the top components (bimolecular reaction eval, math,
  allocation), which raises confidence. Propagation/geometry and spatial binning
  rank higher in the clathrin case because complexes grow over the longer run.
- **gperftools does not work on arm64 macOS** (broken unwinding); the call graph
  above comes from `sample` + `instrumentsToPprof`. gperftools remains the right
  tool on Linux.
- Not yet done: Linux `perf`/gperftools corroboration and MPI per-rank profiling.

## Reporting Rules

- Do not claim a hotspot unless a profiler artifact supports it.
- Keep raw profiler outputs outside source control unless they are small,
  review-useful summaries.
- Link every optimization idea to both a benchmark case and a regression guard.
- Record profiler limitations, especially missing symbols, short runs, profiler
  overhead, and MPI rank imbalance.
- Preserve scientific behavior unless a separate reviewed science change
  explicitly approves otherwise.
