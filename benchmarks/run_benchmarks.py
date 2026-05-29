#!/usr/bin/env python3
"""Run NERDSS benchmark cases and emit machine-readable results.

The harness is intentionally independent from correctness/regression runners.
It copies each benchmark input tree to an isolated output directory, applies
manifest-defined parameter overrides to the copied input file, and runs NERDSS
from that copy so repository fixtures are never modified.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional

try:
    import resource
except ImportError:  # pragma: no cover - resource is unavailable on Windows.
    resource = None  # type: ignore[assignment]


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "benchmarks" / "benchmark_manifest.json"
DEFAULT_RESULTS_ROOT = REPO_ROOT / "benchmark_results"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the NERDSS benchmark manifest and write JSON/CSV results."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=f"Benchmark manifest JSON. Default: {DEFAULT_MANIFEST}",
    )
    parser.add_argument(
        "--nerdss",
        type=Path,
        help="Path to an existing NERDSS executable. Defaults to ./bin/nerdss.",
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="Run the selected build command before benchmarking.",
    )
    parser.add_argument(
        "--build-target",
        default="serial",
        help="Make target used with --build. Default: serial.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for this benchmark run. Defaults under ./benchmark_results/.",
    )
    parser.add_argument(
        "--format",
        choices=("json", "csv", "both"),
        default="both",
        help="Result format to write. Default: both.",
    )
    parser.add_argument(
        "--case",
        action="append",
        dest="case_ids",
        help="Case id to run. May be repeated. Default: all manifest cases.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List cases from the manifest without running them.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        help="Override timeout in seconds for every selected case.",
    )
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Additional argument passed to every NERDSS invocation.",
    )
    parser.add_argument(
        "--keep-workdirs",
        action="store_true",
        help="Keep copied case work directories. They are kept by default on failures.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Prepare commands and write result files without executing NERDSS.",
    )
    return parser.parse_args()


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def timestamp_for_path() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def run_text_command(args: List[str], cwd: Path) -> Dict[str, Any]:
    try:
        completed = subprocess.run(
            args,
            cwd=str(cwd),
            text=True,
            capture_output=True,
            check=False,
        )
        return {
            "command": args,
            "exit_status": completed.returncode,
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
        }
    except OSError as exc:
        return {
            "command": args,
            "exit_status": None,
            "stdout": "",
            "stderr": str(exc),
        }


def git_metadata() -> Dict[str, Any]:
    commit = run_text_command(["git", "rev-parse", "HEAD"], REPO_ROOT)
    branch = run_text_command(["git", "branch", "--show-current"], REPO_ROOT)
    status = run_text_command(["git", "status", "--short"], REPO_ROOT)
    return {
        "commit": commit["stdout"] or None,
        "branch": branch["stdout"] or None,
        "dirty": bool(status["stdout"]),
        "status_short": status["stdout"].splitlines(),
    }


def total_memory_bytes() -> Optional[int]:
    if not hasattr(os, "sysconf"):
        return None
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        pages = os.sysconf("SC_PHYS_PAGES")
    except (OSError, ValueError):
        return None
    return int(page_size * pages)


def system_metadata() -> Dict[str, Any]:
    return {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": sys.version.replace("\n", " "),
        "cpu_count": os.cpu_count(),
        "total_memory_bytes": total_memory_bytes(),
    }


def load_manifest(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if "cases" not in manifest or not isinstance(manifest["cases"], list):
        raise ValueError("manifest must contain a 'cases' list")
    return manifest


def merge_parameter_overrides(
    defaults: Mapping[str, Any], case: Mapping[str, Any]
) -> Dict[str, str]:
    merged: Dict[str, str] = {}
    for source in (
        defaults.get("parameter_overrides", {}),
        case.get("parameter_overrides", {}),
    ):
        if not isinstance(source, Mapping):
            continue
        for key, value in source.items():
            merged[str(key)] = str(value)
    return merged


def resolve_repo_path(path: Path) -> Path:
    if path.is_absolute():
        return path
    return (REPO_ROOT / path).resolve()


def copy_case_inputs(source_dir: Path, destination: Path) -> None:
    if not source_dir.is_dir():
        raise FileNotFoundError(f"case source_dir does not exist: {source_dir}")
    shutil.copytree(source_dir, destination)


PARAMETER_LINE_RE = re.compile(
    r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)(\s*(?:#.*)?)$"
)


def apply_parameter_overrides(input_path: Path, overrides: Mapping[str, str]) -> None:
    if not overrides:
        return
    lines = input_path.read_text(encoding="utf-8").splitlines(keepends=True)
    pending = {key.lower(): (key, value) for key, value in overrides.items()}
    in_parameters = False
    output: List[str] = []

    for line in lines:
        stripped = line.strip().lower()
        if stripped == "start parameters":
            in_parameters = True
            output.append(line)
            continue

        if in_parameters and stripped == "end parameters":
            for original_key, value in pending.values():
                output.append(f"    {original_key} = {value}\n")
            pending.clear()
            in_parameters = False
            output.append(line)
            continue

        if in_parameters:
            line_ending = "\n" if line.endswith("\n") else ""
            body = line[:-1] if line_ending else line
            match = PARAMETER_LINE_RE.match(body)
            if match and match.group(2).lower() in pending:
                key_lower = match.group(2).lower()
                original_key, value = pending.pop(key_lower)
                comment = match.group(4) or ""
                output.append(
                    f"{match.group(1)}{original_key} = {value}{comment}{line_ending}"
                )
                continue

        output.append(line)

    if pending:
        missing = ", ".join(sorted(key for key, _ in pending.values()))
        raise ValueError(f"input file has no parameters block for overrides: {missing}")

    input_path.write_text("".join(output), encoding="utf-8")


def snapshot_files(root: Path) -> Dict[str, int]:
    files: Dict[str, int] = {}
    for path in root.rglob("*"):
        if path.is_file():
            files[str(path.relative_to(root))] = path.stat().st_size
    return files


def produced_outputs(root: Path, before: Mapping[str, int]) -> Dict[str, Any]:
    after = snapshot_files(root)
    produced: List[Dict[str, Any]] = []
    total_bytes = 0
    for rel_path, size in sorted(after.items()):
        if before.get(rel_path) == size:
            continue
        produced.append({"path": rel_path, "bytes": size})
        total_bytes += size
    return {
        "count": len(produced),
        "bytes": total_bytes,
        "files": produced,
    }


def rusage_snapshot() -> Optional[Any]:
    if resource is None:
        return None
    return resource.getrusage(resource.RUSAGE_CHILDREN)


def rusage_delta(before: Optional[Any], after: Optional[Any]) -> Dict[str, Optional[float]]:
    if before is None or after is None:
        return {
            "cpu_time_seconds": None,
            "user_cpu_seconds": None,
            "system_cpu_seconds": None,
            "max_rss_bytes": None,
        }
    user_cpu = max(0.0, after.ru_utime - before.ru_utime)
    system_cpu = max(0.0, after.ru_stime - before.ru_stime)
    max_rss = None
    if after.ru_maxrss > before.ru_maxrss:
        max_rss = after.ru_maxrss
        if platform.system() != "Darwin":
            max_rss *= 1024
    return {
        "cpu_time_seconds": user_cpu + system_cpu,
        "user_cpu_seconds": user_cpu,
        "system_cpu_seconds": system_cpu,
        "max_rss_bytes": int(max_rss) if max_rss is not None else None,
    }


def read_process_rss_bytes(pid: int) -> Optional[int]:
    """Return current resident set size for pid using ps, when available."""
    try:
        completed = subprocess.run(
            ["ps", "-o", "rss=", "-p", str(pid)],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        return None
    if completed.returncode != 0:
        return None
    raw = completed.stdout.strip()
    if not raw:
        return None
    try:
        return int(raw.splitlines()[0].strip()) * 1024
    except ValueError:
        return None


def relative_to_repo(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def run_case(
    case: Mapping[str, Any],
    defaults: Mapping[str, Any],
    executable: Path,
    output_dir: Path,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    case_id = str(case["id"])
    work_dir = output_dir / "workdirs" / case_id
    input_dir = work_dir / "input"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    source_dir = resolve_repo_path(Path(str(case["source_dir"])))
    copy_case_inputs(source_dir, input_dir)

    input_file = input_dir / str(case["input_file"])
    if not input_file.is_file():
        raise FileNotFoundError(f"case input_file does not exist after copy: {input_file}")

    parameter_overrides = merge_parameter_overrides(defaults, case)
    apply_parameter_overrides(input_file, parameter_overrides)

    timeout = float(
        args.timeout
        or case.get("timeout_seconds")
        or defaults.get("timeout_seconds")
        or 300
    )
    seed = case.get("seed", defaults.get("seed"))
    command = [str(executable), "-f", str(case["input_file"])]
    if seed is not None:
        command.extend(["-s", str(seed)])
    command.extend(args.extra_arg)

    before_files = snapshot_files(input_dir)
    stdout_path = work_dir / "stdout.txt"
    stderr_path = work_dir / "stderr.txt"
    started = utc_now()
    wall_start = time.perf_counter()
    usage_before = rusage_snapshot()
    exit_status: Optional[int] = None
    timed_out = False
    error: Optional[str] = None
    peak_rss_bytes: Optional[int] = None

    if args.dry_run:
        exit_status = None
    else:
        try:
            with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
                "w", encoding="utf-8"
            ) as stderr_handle:
                process = subprocess.Popen(
                    command,
                    cwd=str(input_dir),
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    text=True,
                )
                deadline = time.perf_counter() + timeout
                while True:
                    rss_bytes = read_process_rss_bytes(process.pid)
                    if rss_bytes is not None:
                        peak_rss_bytes = max(peak_rss_bytes or 0, rss_bytes)
                    exit_status = process.poll()
                    if exit_status is not None:
                        break
                    if time.perf_counter() >= deadline:
                        timed_out = True
                        process.kill()
                        exit_status = process.wait()
                        with stderr_path.open("a", encoding="utf-8") as stderr_append:
                            stderr_append.write(
                                f"\nBenchmark runner timed out after {timeout} seconds.\n"
                            )
                        break
                    time.sleep(0.2)
        except OSError as exc:
            error = str(exc)
            exit_status = None
            stderr_path.write_text(error + "\n", encoding="utf-8")

    wall_time = time.perf_counter() - wall_start
    usage_after = rusage_snapshot()
    ended = utc_now()
    output_summary = produced_outputs(input_dir, before_files)
    stdout_bytes = stdout_path.stat().st_size if stdout_path.exists() else 0
    stderr_bytes = stderr_path.stat().st_size if stderr_path.exists() else 0

    status = "dry_run"
    if not args.dry_run:
        if timed_out:
            status = "timeout"
        elif exit_status == 0:
            status = "passed"
        else:
            status = "failed"

    usage = rusage_delta(usage_before, usage_after)
    if peak_rss_bytes is not None:
        usage["max_rss_bytes"] = peak_rss_bytes
        max_rss_source = "ps_rss_poll"
    else:
        max_rss_source = "resource_rusage" if usage["max_rss_bytes"] is not None else None

    result: Dict[str, Any] = {
        "case_id": case_id,
        "label": case.get("label"),
        "size": case.get("size"),
        "source_dir": relative_to_repo(source_dir),
        "input_file": str(case["input_file"]),
        "work_dir": str(work_dir),
        "command": command,
        "timeout_seconds": timeout,
        "seed": seed,
        "parameter_overrides": parameter_overrides,
        "started_at": started,
        "ended_at": ended,
        "status": status,
        "exit_status": exit_status,
        "timed_out": timed_out,
        "wall_time_seconds": wall_time,
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
        "stdout_bytes": stdout_bytes,
        "stderr_bytes": stderr_bytes,
        "produced_output_count": output_summary["count"],
        "produced_output_bytes": output_summary["bytes"],
        "produced_outputs": output_summary["files"],
        "error": error,
        "max_rss_source": max_rss_source,
    }
    result.update(usage)

    if status == "passed" and not args.keep_workdirs:
        shutil.rmtree(work_dir)

    return result


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, cases: Iterable[Mapping[str, Any]]) -> None:
    fieldnames = [
        "case_id",
        "label",
        "size",
        "status",
        "exit_status",
        "timed_out",
        "wall_time_seconds",
        "cpu_time_seconds",
        "user_cpu_seconds",
        "system_cpu_seconds",
        "max_rss_bytes",
        "max_rss_source",
        "produced_output_count",
        "produced_output_bytes",
        "stdout_bytes",
        "stderr_bytes",
        "timeout_seconds",
        "seed",
        "work_dir",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in cases:
            writer.writerow({field: case.get(field) for field in fieldnames})


def select_cases(
    manifest: Mapping[str, Any], case_ids: Optional[List[str]]
) -> List[Mapping[str, Any]]:
    cases = manifest["cases"]
    if not case_ids:
        return list(cases)
    wanted = set(case_ids)
    selected = [case for case in cases if case.get("id") in wanted]
    found = {str(case.get("id")) for case in selected}
    missing = sorted(wanted - found)
    if missing:
        raise ValueError(f"unknown benchmark case id(s): {', '.join(missing)}")
    return selected


def build_executable(output_dir: Path, target: str) -> Dict[str, Any]:
    stdout_path = output_dir / "build_stdout.txt"
    stderr_path = output_dir / "build_stderr.txt"
    command = ["make", target]
    started = utc_now()
    wall_start = time.perf_counter()
    try:
        with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
            "w", encoding="utf-8"
        ) as stderr_handle:
            completed = subprocess.run(
                command,
                cwd=str(REPO_ROOT),
                stdout=stdout_handle,
                stderr=stderr_handle,
                text=True,
                check=False,
            )
        exit_status: Optional[int] = completed.returncode
        error = None
    except OSError as exc:
        exit_status = None
        error = str(exc)
        stderr_path.write_text(error + "\n", encoding="utf-8")
    return {
        "command": command,
        "started_at": started,
        "ended_at": utc_now(),
        "wall_time_seconds": time.perf_counter() - wall_start,
        "exit_status": exit_status,
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
        "error": error,
    }


def write_results(output_dir: Path, result_format: str, payload: Mapping[str, Any]) -> None:
    if result_format in ("json", "both"):
        write_json(output_dir / "benchmark_results.json", payload)
    if result_format in ("csv", "both"):
        write_csv(output_dir / "benchmark_results.csv", payload.get("cases", []))


def main() -> int:
    args = parse_args()
    manifest_path = resolve_repo_path(args.manifest)
    manifest = load_manifest(manifest_path)
    selected_cases = select_cases(manifest, args.case_ids)

    if args.list:
        for case in selected_cases:
            print(f"{case['id']}\t{case.get('size', '')}\t{case.get('label', '')}")
        return 0

    output_dir = (args.output_dir or (DEFAULT_RESULTS_ROOT / timestamp_for_path())).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    executable = resolve_repo_path(args.nerdss) if args.nerdss else (REPO_ROOT / "bin" / "nerdss")
    run_payload: MutableMapping[str, Any] = {
        "schema_version": 1,
        "created_at": utc_now(),
        "manifest": str(manifest_path),
        "output_dir": str(output_dir),
        "repo": git_metadata(),
        "system": system_metadata(),
        "runner": {
            "script": str(Path(__file__).resolve()),
            "dry_run": args.dry_run,
            "build_requested": args.build,
            "build_target": args.build_target,
            "executable": str(executable),
            "extra_args": args.extra_arg,
        },
        "build": None,
        "cases": [],
    }

    if args.build:
        build_result = build_executable(output_dir, args.build_target)
        run_payload["build"] = build_result
        if build_result["exit_status"] != 0:
            run_payload["error"] = (
                "build failed; pass --nerdss PATH to use an executable from another build "
                "branch or rerun after the serial build fix lands"
            )
            write_results(output_dir, args.format, run_payload)
            print(run_payload["error"], file=sys.stderr)
            return 2

    if not executable.is_file() and not args.dry_run:
        run_payload["error"] = (
            f"NERDSS executable not found: {executable}. "
            "Build with --build or pass --nerdss PATH to an existing executable."
        )
        write_results(output_dir, args.format, run_payload)
        print(run_payload["error"], file=sys.stderr)
        return 2

    defaults = manifest.get("defaults", {})
    exit_code = 0
    for case in selected_cases:
        try:
            result = run_case(case, defaults, executable, output_dir, args)
        except Exception as exc:  # Keep one broken case from hiding prior results.
            result = {
                "case_id": case.get("id"),
                "label": case.get("label"),
                "size": case.get("size"),
                "status": "harness_error",
                "exit_status": None,
                "timed_out": False,
                "error": str(exc),
            }
        run_payload["cases"].append(result)
        if result.get("status") not in ("passed", "dry_run"):
            exit_code = 1
        write_results(output_dir, args.format, run_payload)

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
