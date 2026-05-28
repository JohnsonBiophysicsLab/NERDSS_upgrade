# Benchmark Harness Baseline

This document records the Phase 6 benchmark harness for the NERDSS upgrade.
It is a harness baseline, not a performance baseline: current `master` may not
build in every agent worktree until the build-fix branch is merged.

## Scope

- Runner: `benchmarks/run_benchmarks.py`
- Manifest: `benchmarks/benchmark_manifest.json`
- Generated results: `benchmark_results/<timestamp>/benchmark_results.json`
  and `benchmark_results/<timestamp>/benchmark_results.csv`

The benchmark runner is independent from the regression harness. It does not
compare scientific outputs or enforce tolerances. Its job is to collect timing,
resource, output-volume, git, and system metadata for later performance work.

## Benchmark Cases

The manifest currently defines four serial cases:

| Case | Size | Source input | Purpose |
| --- | --- | --- | --- |
| `small_homotrimer` | small | `sample_inputs/VALIDATE_SUITE/homoTrimer/parmTri6.inp` | Small 3D association workload. |
| `medium_michaelis_menten` | medium | `sample_inputs/VALIDATE_SUITE/michaelis_menten/michaelis.inp` | Moderate reaction-network workload. |
| `representative_implicit_lipid` | representative | `sample_inputs/VALIDATE_SUITE/implicit_lipid/parms.inp` | Implicit-lipid binding workload. |
| `large_clathrin_short` | large | `sample_inputs/VALIDATE_SUITE/clathrin/parms_clath_kon1uM.inp` | Larger model run at short benchmark iteration count. |

Each case copies its source directory into an isolated benchmark work directory
and applies benchmark-only `start parameters` overrides there. The original
sample inputs are not modified.

## Usage

List benchmark cases:

```sh
python3 benchmarks/run_benchmarks.py --list
```

Run all cases with the repository serial executable:

```sh
python3 benchmarks/run_benchmarks.py --nerdss ./bin/nerdss
```

Build first, then run all cases:

```sh
python3 benchmarks/run_benchmarks.py --build
```

Run a single case:

```sh
python3 benchmarks/run_benchmarks.py --nerdss ./bin/nerdss --case small_homotrimer
```

If current `master` cannot build locally, use a known-good executable from
Agent B's build branch once that branch is merged or checked out separately:

```sh
python3 benchmarks/run_benchmarks.py --nerdss /path/to/nerdss
```

## Collected Fields

The JSON result includes:

- Git branch, commit, and dirty status.
- Platform, hostname, CPU count, Python version, and total memory when the OS
  exposes it.
- Optional build command, build exit status, and build logs.
- Per-case command, seed, timeout, wall time, CPU time when available, peak RSS
  when available, RSS source, exit status, timeout status, stdout/stderr log
  paths and byte counts, and produced output file count/bytes.

The CSV result flattens the most useful per-case fields for spreadsheets and
quick trend comparisons.

## Baseline Notes

No numeric performance baseline is committed yet. The first numeric baseline
should be generated after the serial build fix is available on the branch used
for benchmarking, then attached to this document or stored as an external CI
artifact. Avoid committing large raw output directories.
