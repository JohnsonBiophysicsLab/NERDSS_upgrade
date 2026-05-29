# NERDSS Profiling And Hotspot Report

Owner: Agent P

Status: initial template and profiler command guide. Fill measured results after
Agent O benchmark cases are available.

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

### gperftools On macOS

Homebrew gperftools can be useful, but availability varies by architecture and
compiler. Record `brew --prefix gperftools`, compiler, and linker flags if this
path is used. Treat Instruments as the default macOS source of truth.

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

## Reporting Rules

- Do not claim a hotspot unless a profiler artifact supports it.
- Keep raw profiler outputs outside source control unless they are small,
  review-useful summaries.
- Link every optimization idea to both a benchmark case and a regression guard.
- Record profiler limitations, especially missing symbols, short runs, profiler
  overhead, and MPI rank imbalance.
- Preserve scientific behavior unless a separate reviewed science change
  explicitly approves otherwise.
