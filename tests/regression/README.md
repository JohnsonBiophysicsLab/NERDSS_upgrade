# NERDSS Regression Harness

This harness runs selected small validation samples from
`sample_inputs/VALIDATE_SUITE` in isolated temporary directories and compares
their outputs against either:

- a generated two-pass self-baseline from the same binary, or
- a stored baseline directory provided with `--baseline-root`.

It is intentionally separate from the smoke runner. The smoke runner answers
"can this checkout build and run at all?" while this harness answers "did fixed
seed validation outputs change?"

## Quick start

Build the serial executable, then run the default cases:

```sh
make serial
python3 tests/regression/run_regression.py --binary bin/nerdss
```

The runner can build first if desired:

```sh
python3 tests/regression/run_regression.py --build
```

List cases or run only one case:

```sh
python3 tests/regression/run_regression.py --list
python3 tests/regression/run_regression.py --case create_destroy_fresh_small
python3 tests/regression/run_regression.py --dry-run --case create_destroy_stochastic_seed_set
```

Use `--json-output <path>` with `--list`, `--dry-run`, or a real run to emit a
machine-readable report containing selected cases, planned commands, run
directories, stochastic summaries, threshold status, and failures.

## Baselines

Without `--baseline-root`, the runner executes each case twice with the same
fixed seed in two separate temp directories and compares the second run against
the first. This is the Phase 1 "current branch against itself" check.

To write stored baselines:

```sh
python3 tests/regression/run_regression.py --update-baseline /tmp/nerdss-baseline
```

To compare against stored baselines:

```sh
python3 tests/regression/run_regression.py --baseline-root /tmp/nerdss-baseline
```

Each stored case directory contains the captured `stdout.txt`, `stderr.txt`,
`exit_code.txt`, and produced simulation outputs.

## Manifest comparison modes

`manifest.json` defines cases, fixed seeds, staged sample files, optional input
patches for shortening validation samples, and file-level comparisons.

- `exact`: byte-for-byte text comparison after optional ignored-line filtering.
- `tolerant`: text comparison that requires nonnumeric text to match exactly
  while numeric tokens compare with `atol` and `rtol`.

Both modes support `ignore_lines_matching` for volatile text such as dates,
command paths, and timing lines.

## Stochastic seed-set checks

Stochastic cases use `type: "stochastic_ensemble"` and run the same shortened
sample once for each seed in `seed_set`. They do not compare byte-for-byte
outputs. Instead, they read numeric output tables, extract named metrics, compute
summary statistics across the fixed seed set, and validate inclusive thresholds.

Metric fields:

- `path`: output file relative to each run directory, for example
  `DATA/copy_numbers_time.dat`.
- `column`: zero-based numeric column index or an exact header name.
- `value`: row/series selector; supported values are `first`, `last`, `min`,
  `max`, and `mean`.
- `statistics`: any of `count`, `mean`, `stdev`, `min`, `max`, `range`, `sem`.
- `thresholds`: per-statistic checks using `min`, `max`, or `equals`.

The initial Phase 1 stochastic case is
`create_destroy_stochastic_seed_set`. It runs a five-seed ensemble and checks
the final created-species count stays finite, nonnegative, and within the
configured envelope. The JSON report includes per-seed metric values and the
computed ensemble summary.

## Restart support

Restart cases use `mode: "restart"` and `restart_file`. The runner invokes
NERDSS with `--restart <file>` and can stage a matching `rng_state_file` into
the run directory when a restart fixture includes one. The initial wired case
uses `homoTrimer/RESTARTS/restart1000.dat` with shortened `numItr`.
