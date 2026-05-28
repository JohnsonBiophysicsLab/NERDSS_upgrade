#!/usr/bin/env python3
"""Run deterministic NERDSS validation samples and compare outputs."""

from __future__ import annotations

import argparse
import difflib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = Path(__file__).with_name("manifest.json")
NUMBER_RE = re.compile(
    r"[-+]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|nan|inf)",
    re.IGNORECASE,
)


@dataclass
class RunResult:
    case_name: str
    run_dir: Path
    seed: int | None
    returncode: int
    stdout_path: Path
    stderr_path: Path


class RegressionError(Exception):
    """Raised for validation failures with a concise user-facing message."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run fixed-seed NERDSS regression cases from VALIDATE_SUITE."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Path to the regression manifest JSON.",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=REPO_ROOT / "bin" / "nerdss",
        help="NERDSS serial executable to run.",
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="Run `make serial` before executing cases.",
    )
    parser.add_argument(
        "--case",
        action="append",
        dest="case_names",
        help="Case name to run. May be provided more than once.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List manifest cases and exit.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate selected cases and print planned NERDSS commands without running.",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        help="Write a machine-readable JSON report for list, dry-run, or run results.",
    )
    parser.add_argument(
        "--baseline-root",
        type=Path,
        help="Directory containing stored baseline case outputs.",
    )
    parser.add_argument(
        "--update-baseline",
        type=Path,
        help="Generate or replace stored baseline case outputs in this directory.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary run directories for debugging.",
    )
    parser.add_argument(
        "--tmp-root",
        type=Path,
        help="Parent directory for temporary run directories.",
    )
    return parser.parse_args()


def load_manifest(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("version") != 1:
        raise RegressionError(f"Unsupported manifest version in {path}")
    if not isinstance(manifest.get("cases"), list):
        raise RegressionError(f"Manifest {path} must contain a cases list")
    return manifest


def selected_cases(
    manifest: dict[str, Any], requested_names: list[str] | None
) -> list[dict[str, Any]]:
    cases = manifest["cases"]
    if not requested_names:
        return cases

    by_name = {case["name"]: case for case in cases}
    unknown = sorted(set(requested_names) - set(by_name))
    if unknown:
        raise RegressionError(f"Unknown regression case(s): {', '.join(unknown)}")
    return [by_name[name] for name in requested_names]


def case_type(case: dict[str, Any]) -> str:
    return case.get("type", case.get("validation_type", "deterministic"))


def is_stochastic_case(case: dict[str, Any]) -> bool:
    return case_type(case) in {"stochastic", "stochastic_ensemble", "seed_set"}


def seed_set(case: dict[str, Any], defaults: dict[str, Any]) -> list[int]:
    seeds = case.get("seed_set", defaults.get("stochastic_seed_set"))
    if not seeds:
        raise RegressionError(f"Stochastic case {case['name']} needs seed_set")
    try:
        return [int(seed) for seed in seeds]
    except (TypeError, ValueError) as exc:
        raise RegressionError(f"Case {case['name']} has an invalid seed_set") from exc


def resolve_binary_path(binary: Path) -> Path:
    return binary if binary.is_absolute() else REPO_ROOT / binary


def build_binary() -> None:
    print("Building serial NERDSS executable...")
    subprocess.run(["make", "serial"], cwd=REPO_ROOT, check=True)


def require_binary(binary: Path) -> Path:
    resolved = resolve_binary_path(binary)
    if not resolved.exists():
        raise RegressionError(
            f"NERDSS binary not found at {resolved}. Run `make serial` or pass --build."
        )
    if not resolved.is_file():
        raise RegressionError(f"NERDSS binary path is not a file: {resolved}")
    return resolved


def copy_stage_files(case: dict[str, Any], run_dir: Path) -> None:
    source_root = REPO_ROOT / case["sample_dir"]
    if not source_root.is_dir():
        raise RegressionError(f"Sample directory does not exist: {source_root}")

    stage_files = case.get("stage_files")
    if not stage_files:
        stage_files = [path.name for path in source_root.iterdir()]

    for rel_name in stage_files:
        source = source_root / rel_name
        target = run_dir / rel_name
        if not source.exists():
            raise RegressionError(f"Case {case['name']} missing staged file: {source}")
        target.parent.mkdir(parents=True, exist_ok=True)
        if source.is_dir():
            if target.exists():
                shutil.rmtree(target)
            shutil.copytree(source, target)
        else:
            shutil.copy2(source, target)


def apply_patches(case: dict[str, Any], run_dir: Path) -> None:
    for patch in case.get("patches", []):
        path = run_dir / patch["path"]
        text = path.read_text(encoding="latin-1")
        updated, count = re.subn(
            patch["pattern"], patch["replace"], text, count=patch.get("count", 1)
        )
        if count == 0:
            raise RegressionError(
                f"Case {case['name']} patch did not match {patch['path']}: "
                f"{patch['pattern']}"
            )
        path.write_text(updated, encoding="latin-1")


def stage_rng_state(case: dict[str, Any], run_dir: Path) -> None:
    rng_state = case.get("rng_state_file") or case.get("restart_support", {}).get(
        "rng_state_file"
    )
    if not rng_state:
        return
    source = REPO_ROOT / case["sample_dir"] / rng_state
    if not source.exists():
        raise RegressionError(f"Case {case['name']} missing rng state file: {source}")
    shutil.copy2(source, run_dir / "rng_state")


def make_command(
    binary: Path,
    case: dict[str, Any],
    defaults: dict[str, Any],
    seed_override: int | None = None,
) -> list[str]:
    seed = str(seed_override if seed_override is not None else case.get("seed", defaults.get("seed")))
    mode = case.get("mode", "fresh")
    if mode == "fresh":
        input_file = case.get("input_file")
        if not input_file:
            raise RegressionError(f"Fresh case {case['name']} needs input_file")
        command = [str(binary), "-f", input_file, "--seed", seed]
    elif mode == "restart":
        restart_file = case.get("restart_file")
        if not restart_file:
            raise RegressionError(f"Restart case {case['name']} needs restart_file")
        command = [str(binary), "--restart", restart_file, "--seed", seed]
    else:
        raise RegressionError(f"Case {case['name']} has unsupported mode: {mode}")

    command.extend(case.get("extra_args", []))
    return command


def prepare_run_dir(case: dict[str, Any], run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    copy_stage_files(case, run_dir)
    apply_patches(case, run_dir)
    stage_rng_state(case, run_dir)


def run_case(
    binary: Path,
    case: dict[str, Any],
    defaults: dict[str, Any],
    parent_dir: Path,
    run_label: str,
    seed_override: int | None = None,
) -> RunResult:
    run_dir = parent_dir / case["name"] / run_label
    if run_dir.exists():
        shutil.rmtree(run_dir)
    prepare_run_dir(case, run_dir)

    command = make_command(binary, case, defaults, seed_override)
    timeout = int(case.get("timeout_seconds", defaults.get("timeout_seconds", 120)))
    completed = subprocess.run(
        command,
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )

    stdout_path = run_dir / "stdout.txt"
    stderr_path = run_dir / "stderr.txt"
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    (run_dir / "exit_code.txt").write_text(f"{completed.returncode}\n", encoding="utf-8")
    return RunResult(
        case["name"], run_dir, seed_override, completed.returncode, stdout_path, stderr_path
    )


def copy_baseline(run_dir: Path, baseline_case_dir: Path) -> None:
    if baseline_case_dir.exists():
        shutil.rmtree(baseline_case_dir)
    shutil.copytree(run_dir, baseline_case_dir)


def filtered_lines(text: str, ignore_patterns: list[str]) -> list[str]:
    patterns = [re.compile(pattern) for pattern in ignore_patterns]
    lines = text.replace("\r\n", "\n").replace("\r", "\n").splitlines()
    if not patterns:
        return lines
    return [line for line in lines if not any(pattern.search(line) for pattern in patterns)]


def read_filtered(path: Path, comparison: dict[str, Any]) -> list[str]:
    if not path.exists():
        raise RegressionError(f"Expected comparison file is missing: {path}")
    text = path.read_text(encoding="latin-1")
    return filtered_lines(text, comparison.get("ignore_lines_matching", []))


def is_number(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def compare_numeric_line(
    baseline: str, candidate: str, rtol: float, atol: float, label: str
) -> str | None:
    base_matches = list(NUMBER_RE.finditer(baseline))
    cand_matches = list(NUMBER_RE.finditer(candidate))
    if len(base_matches) != len(cand_matches):
        return f"{label}: numeric token count differs"

    base_pos = 0
    cand_pos = 0
    for index, (base_match, cand_match) in enumerate(zip(base_matches, cand_matches)):
        if baseline[base_pos : base_match.start()] != candidate[cand_pos : cand_match.start()]:
            return f"{label}: nonnumeric text differs before numeric token {index + 1}"

        base_token = base_match.group(0)
        cand_token = cand_match.group(0)
        if is_number(base_token) and is_number(cand_token):
            base_value = float(base_token)
            cand_value = float(cand_token)
            if math.isnan(base_value) and math.isnan(cand_value):
                pass
            elif not math.isclose(base_value, cand_value, rel_tol=rtol, abs_tol=atol):
                return (
                    f"{label}: numeric token {index + 1} differs: "
                    f"{base_token} != {cand_token} (atol={atol}, rtol={rtol})"
                )
        elif base_token != cand_token:
            return f"{label}: token {index + 1} differs: {base_token} != {cand_token}"

        base_pos = base_match.end()
        cand_pos = cand_match.end()

    if baseline[base_pos:] != candidate[cand_pos:]:
        return f"{label}: trailing nonnumeric text differs"
    return None


def compare_exact(
    baseline_lines: list[str], candidate_lines: list[str], rel_path: str
) -> list[str]:
    if baseline_lines == candidate_lines:
        return []
    diff = difflib.unified_diff(
        baseline_lines,
        candidate_lines,
        fromfile=f"baseline/{rel_path}",
        tofile=f"candidate/{rel_path}",
        lineterm="",
        n=3,
    )
    return [f"{rel_path} exact comparison failed:", *list(diff)[:40]]


def compare_tolerant(
    baseline_lines: list[str],
    candidate_lines: list[str],
    comparison: dict[str, Any],
    rel_path: str,
) -> list[str]:
    if len(baseline_lines) != len(candidate_lines):
        return [
            f"{rel_path} line count differs: "
            f"{len(baseline_lines)} != {len(candidate_lines)}"
        ]

    rtol = float(comparison.get("rtol", 0.0))
    atol = float(comparison.get("atol", 0.0))
    for line_no, (baseline, candidate) in enumerate(
        zip(baseline_lines, candidate_lines), start=1
    ):
        problem = compare_numeric_line(
            baseline, candidate, rtol, atol, f"{rel_path}:{line_no}"
        )
        if problem:
            return [problem]
    return []


def compare_case_outputs(
    case: dict[str, Any], baseline_dir: Path, candidate_dir: Path
) -> list[str]:
    failures: list[str] = []
    expected_exit_code = case.get("expected_exit_code")
    if expected_exit_code is not None:
        actual = (candidate_dir / "exit_code.txt").read_text(encoding="utf-8").strip()
        if actual != str(expected_exit_code):
            failures.append(
                f"exit_code.txt expected {expected_exit_code}, got {actual}. "
                f"See {candidate_dir / 'stderr.txt'}"
            )

    for comparison in case.get("comparisons", []):
        rel_path = comparison["path"]
        baseline_path = baseline_dir / rel_path
        candidate_path = candidate_dir / rel_path
        baseline_lines = read_filtered(baseline_path, comparison)
        candidate_lines = read_filtered(candidate_path, comparison)

        mode = comparison.get("mode", "exact")
        if mode == "exact":
            failures.extend(compare_exact(baseline_lines, candidate_lines, rel_path))
        elif mode == "tolerant":
            failures.extend(
                compare_tolerant(baseline_lines, candidate_lines, comparison, rel_path)
            )
        else:
            failures.append(f"{rel_path} has unsupported comparison mode: {mode}")
    return failures


def parse_numeric_table(path: Path) -> tuple[list[str], list[list[float]]]:
    if not path.exists():
        raise RegressionError(f"Expected stochastic statistics file is missing: {path}")

    lines = [
        line.strip()
        for line in path.read_text(encoding="latin-1").splitlines()
        if line.strip()
    ]
    if not lines:
        raise RegressionError(f"Statistics file is empty: {path}")

    delimiter = "," if "," in lines[0] else None
    header = [
        token.strip()
        for token in (lines[0].split(delimiter) if delimiter else lines[0].split())
    ]
    rows: list[list[float]] = []
    for line in lines[1:]:
        tokens = [
            token.strip()
            for token in (line.split(delimiter) if delimiter else line.split())
            if token.strip()
        ]
        try:
            rows.append([float(token) for token in tokens])
        except ValueError:
            continue

    if not rows:
        raise RegressionError(f"Statistics file has no numeric rows: {path}")
    return header, rows


def column_index(header: list[str], column: int | str, path: Path) -> int:
    if isinstance(column, int):
        index = column
    elif isinstance(column, str) and column.isdigit():
        index = int(column)
    elif isinstance(column, str):
        normalized = {name.strip(): index for index, name in enumerate(header)}
        if column not in normalized:
            raise RegressionError(
                f"Column {column!r} not found in {path}; available columns: {header}"
            )
        index = normalized[column]
    else:
        raise RegressionError(f"Invalid column selector {column!r} for {path}")

    if index < 0 or index >= len(header):
        raise RegressionError(f"Column index {index} is out of range for {path}")
    return index


def select_metric_value(rows: list[list[float]], index: int, selector: str, path: Path) -> float:
    values = [row[index] for row in rows if index < len(row)]
    if not values:
        raise RegressionError(f"No values found for column {index} in {path}")

    if selector == "last":
        return values[-1]
    if selector == "first":
        return values[0]
    if selector == "min":
        return min(values)
    if selector == "max":
        return max(values)
    if selector == "mean":
        return sum(values) / len(values)
    raise RegressionError(f"Unsupported value selector {selector!r} for {path}")


def summarize_values(values: list[float], requested: list[str]) -> dict[str, float]:
    count = len(values)
    mean = sum(values) / count
    variance = (
        sum((value - mean) ** 2 for value in values) / (count - 1) if count > 1 else 0.0
    )
    available = {
        "count": float(count),
        "mean": mean,
        "stdev": math.sqrt(variance),
        "min": min(values),
        "max": max(values),
        "range": max(values) - min(values),
    }
    if "sem" in requested:
        available["sem"] = available["stdev"] / math.sqrt(count)
    return {name: available[name] for name in requested if name in available}


def threshold_failures(
    metric_name: str, summary: dict[str, float], thresholds: dict[str, Any]
) -> list[str]:
    failures: list[str] = []
    for statistic, limits in thresholds.items():
        if statistic not in summary:
            failures.append(f"{metric_name}: threshold references missing statistic {statistic}")
            continue
        value = summary[statistic]
        if "min" in limits and value < float(limits["min"]):
            failures.append(f"{metric_name}.{statistic}={value} below {limits['min']}")
        if "max" in limits and value > float(limits["max"]):
            failures.append(f"{metric_name}.{statistic}={value} above {limits['max']}")
        if "equals" in limits and value != float(limits["equals"]):
            failures.append(f"{metric_name}.{statistic}={value} != {limits['equals']}")
    return failures


def analyze_stochastic_case(
    case: dict[str, Any], runs: list[RunResult]
) -> tuple[list[dict[str, Any]], list[str]]:
    metrics: list[dict[str, Any]] = []
    failures: list[str] = []
    for metric in case.get("metrics", []):
        metric_name = metric["name"]
        selector = metric.get("value", "last")
        values: list[float] = []
        per_seed: list[dict[str, Any]] = []
        for run in runs:
            path = run.run_dir / metric["path"]
            header, rows = parse_numeric_table(path)
            index = column_index(header, metric["column"], path)
            value = select_metric_value(rows, index, selector, path)
            values.append(value)
            per_seed.append({"seed": run.seed, "value": value})

        requested = metric.get(
            "statistics", ["count", "mean", "stdev", "min", "max", "range"]
        )
        summary = summarize_values(values, requested)
        metric_failures = threshold_failures(
            metric_name, summary, metric.get("thresholds", {})
        )
        failures.extend(metric_failures)
        metrics.append(
            {
                "name": metric_name,
                "path": metric["path"],
                "column": metric["column"],
                "value": selector,
                "per_seed": per_seed,
                "summary": summary,
                "thresholds": metric.get("thresholds", {}),
                "status": "fail" if metric_failures else "pass",
            }
        )
    return metrics, failures


def run_stochastic_case(
    binary: Path,
    case: dict[str, Any],
    defaults: dict[str, Any],
    temp_root: Path,
) -> tuple[dict[str, Any], list[str]]:
    expected_exit_code = int(
        case.get("expected_exit_code", defaults.get("expected_exit_code", 0))
    )
    seeds = seed_set(case, defaults)
    runs: list[RunResult] = []
    failures: list[str] = []

    for seed in seeds:
        print(f"Running {case['name']} seed {seed}...")
        run = run_case(binary, case, defaults, temp_root, f"seed_{seed}", seed)
        runs.append(run)
        if run.returncode != expected_exit_code:
            failures.append(
                f"{case['name']} seed {seed}: exited {run.returncode}; "
                f"see {run.stderr_path}"
            )

    metrics: list[dict[str, Any]] = []
    if not failures:
        metrics, failures = analyze_stochastic_case(case, runs)

    report = {
        "name": case["name"],
        "type": case_type(case),
        "description": case.get("description", ""),
        "seeds": seeds,
        "runs": [
            {
                "seed": run.seed,
                "run_dir": str(run.run_dir),
                "returncode": run.returncode,
            }
            for run in runs
        ],
        "metrics": metrics,
        "status": "fail" if failures else "pass",
    }
    return report, failures


def planned_case_commands(
    binary: Path, case: dict[str, Any], defaults: dict[str, Any]
) -> list[dict[str, Any]]:
    if is_stochastic_case(case):
        return [
            {"seed": seed, "command": make_command(binary, case, defaults, seed)}
            for seed in seed_set(case, defaults)
        ]
    return [{"seed": case.get("seed", defaults.get("seed")), "command": make_command(binary, case, defaults)}]


def write_json_report(path: Path | None, report: dict[str, Any]) -> None:
    if not path:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_with_temp_root(args: argparse.Namespace, fn: Any) -> int:
    if args.keep_temp:
        root = Path(
            tempfile.mkdtemp(
                prefix="nerdss-regression-",
                dir=str(args.tmp_root) if args.tmp_root else None,
            )
        )
        print(f"Keeping temporary runs under {root}")
        return fn(root)

    with tempfile.TemporaryDirectory(
        prefix="nerdss-regression-", dir=str(args.tmp_root) if args.tmp_root else None
    ) as temp_name:
        return fn(Path(temp_name))


def main() -> int:
    args = parse_args()
    manifest = load_manifest(args.manifest)
    cases = selected_cases(manifest, args.case_names)
    defaults = manifest.get("defaults", {})
    binary_path = resolve_binary_path(args.binary)

    if args.list:
        listed_cases: list[dict[str, Any]] = []
        for case in manifest["cases"]:
            case_info = {
                "name": case["name"],
                "type": case_type(case),
                "description": case.get("description", ""),
            }
            if is_stochastic_case(case):
                case_info["seeds"] = seed_set(case, defaults)
            listed_cases.append(case_info)
            type_label = f" [{case_info['type']}]" if case_info["type"] else ""
            print(f"{case['name']}{type_label}: {case.get('description', '')}")
        write_json_report(
            args.json_output,
            {
                "mode": "list",
                "manifest": str(args.manifest),
                "cases": listed_cases,
                "status": "pass",
            },
        )
        return 0

    if args.dry_run:
        planned: list[dict[str, Any]] = []
        for case in cases:
            commands = planned_case_commands(binary_path, case, defaults)
            planned.append(
                {
                    "name": case["name"],
                    "type": case_type(case),
                    "description": case.get("description", ""),
                    "commands": commands,
                }
            )
            print(f"{case['name']} [{case_type(case)}]")
            for command in commands:
                seed = command.get("seed")
                seed_label = f" seed {seed}:" if seed is not None else ":"
                print(f"  {seed_label} {' '.join(str(part) for part in command['command'])}")
        write_json_report(
            args.json_output,
            {
                "mode": "dry-run",
                "manifest": str(args.manifest),
                "binary": str(binary_path),
                "cases": planned,
                "status": "pass",
            },
        )
        return 0

    if args.build:
        build_binary()
    binary = require_binary(args.binary)

    if args.baseline_root and args.update_baseline:
        raise RegressionError("Use either --baseline-root or --update-baseline, not both")

    def run_all(temp_root: Path) -> int:
        failures: list[str] = []
        case_reports: list[dict[str, Any]] = []
        for case in cases:
            if is_stochastic_case(case):
                case_report, case_failures = run_stochastic_case(
                    binary, case, defaults, temp_root
                )
                case_reports.append(case_report)
                if case_failures:
                    failures.append(f"{case['name']}:")
                    failures.extend(f"  {failure}" for failure in case_failures)
                else:
                    print(f"  PASS {case['name']}")
                continue

            expected_exit_code = int(
                case.get("expected_exit_code", defaults.get("expected_exit_code", 0))
            )
            case["expected_exit_code"] = expected_exit_code
            print(f"Running {case['name']}...")
            report_runs: list[dict[str, Any]] = []

            if args.baseline_root:
                baseline_dir = args.baseline_root / case["name"]
                if not baseline_dir.is_dir():
                    failures.append(f"{case['name']}: missing baseline {baseline_dir}")
                    continue
            else:
                baseline = run_case(binary, case, defaults, temp_root, "baseline")
                if baseline.returncode != expected_exit_code:
                    failures.append(
                        f"{case['name']}: baseline exited {baseline.returncode}; "
                        f"see {baseline.stderr_path}"
                    )
                    continue
                baseline_dir = baseline.run_dir
                report_runs.append(
                    {
                        "label": "baseline",
                        "run_dir": str(baseline.run_dir),
                        "returncode": baseline.returncode,
                    }
                )
                if args.update_baseline:
                    copy_baseline(baseline.run_dir, args.update_baseline / case["name"])

            candidate = run_case(binary, case, defaults, temp_root, "candidate")
            report_runs.append(
                {
                    "label": "candidate",
                    "run_dir": str(candidate.run_dir),
                    "returncode": candidate.returncode,
                }
            )
            if candidate.returncode != expected_exit_code:
                failures.append(
                    f"{case['name']}: candidate exited {candidate.returncode}; "
                    f"see {candidate.stderr_path}"
                )
                continue

            case_failures = compare_case_outputs(case, baseline_dir, candidate.run_dir)
            if case_failures:
                failures.append(f"{case['name']}:")
                failures.extend(f"  {failure}" for failure in case_failures)
            else:
                print(f"  PASS {case['name']}")
            case_reports.append(
                {
                    "name": case["name"],
                    "type": case_type(case),
                    "description": case.get("description", ""),
                    "runs": report_runs,
                    "status": "fail" if case_failures else "pass",
                }
            )

        if failures:
            print("\nRegression failures:")
            for failure in failures:
                print(failure)
            write_json_report(
                args.json_output,
                {
                    "mode": "run",
                    "manifest": str(args.manifest),
                    "binary": str(binary),
                    "cases": case_reports,
                    "failures": failures,
                    "status": "fail",
                },
            )
            return 1

        if args.update_baseline:
            print(f"\nUpdated baselines in {args.update_baseline}")
        print(f"\nPASS: {len(cases)} regression case(s)")
        write_json_report(
            args.json_output,
            {
                "mode": "run",
                "manifest": str(args.manifest),
                "binary": str(binary),
                "cases": case_reports,
                "failures": [],
                "status": "pass",
            },
        )
        return 0

    return run_with_temp_root(args, run_all)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.TimeoutExpired as exc:
        print(f"Timed out after {exc.timeout} seconds: {exc.cmd}", file=sys.stderr)
        raise SystemExit(124)
    except (OSError, RegressionError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
