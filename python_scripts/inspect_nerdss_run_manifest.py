#!/usr/bin/env python3
"""Inspect legacy NERDSS outputs and emit a run manifest skeleton."""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
GENERATOR_VERSION = "0.1.0"

KNOWN_FILES: dict[str, dict[str, str]] = {
    "copy_numbers_time.dat": {
        "role": "legacy_text_timeseries",
        "format": "csv",
        "description": "Copy numbers for molecules and bound products over time.",
        "time_units": "s",
    },
    "observables_time.dat": {
        "role": "legacy_text_timeseries",
        "format": "csv",
        "description": "Configured observables over time.",
        "time_units": "s",
    },
    "bound_pair_time.dat": {
        "role": "legacy_text_timeseries",
        "format": "tsv",
        "description": "Bound pair counts and association event counters over time.",
        "time_units": "s",
    },
    "mono_dimer_time.dat": {
        "role": "legacy_text_timeseries",
        "format": "tsv",
        "description": "Monomer and dimer counts by molecule type over time.",
        "time_units": "s",
    },
    "histogram_complexes_time.dat": {
        "role": "legacy_block_timeseries",
        "format": "legacy_text",
        "description": "Repeated time blocks of complex composition counts.",
        "time_units": "s",
    },
    "event_counters_time.dat": {
        "role": "legacy_block_timeseries",
        "format": "legacy_text",
        "description": "Repeated time blocks of event counters.",
        "time_units": "s",
    },
    "transition_matrix_time.dat": {
        "role": "legacy_block_timeseries",
        "format": "legacy_text",
        "description": "Transition matrix and lifetime text blocks.",
        "time_units": "s",
    },
    "assoc_dissoc_time.dat": {
        "role": "legacy_block_timeseries",
        "format": "legacy_text",
        "description": "Association and dissociation event stream.",
        "time_units": "s",
    },
    "smt_reactions_time.dat": {
        "role": "legacy_block_timeseries",
        "format": "legacy_text",
        "description": "Single-molecule reaction event stream.",
        "time_units": "s",
    },
    "trajectory.xyz": {
        "role": "trajectory_xyz",
        "format": "xyz",
        "description": "XYZ trajectory.",
    },
    "initial_crds.xyz": {
        "role": "coordinate_xyz",
        "format": "xyz",
        "description": "Initial molecule coordinates.",
    },
    "system.psf": {
        "role": "topology_psf",
        "format": "psf",
        "description": "PSF topology.",
    },
    "restart.dat": {
        "role": "restart",
        "format": "restart",
        "description": "Restart state.",
    },
    "rng_state": {
        "role": "rng_state",
        "format": "unknown",
        "description": "Random-number generator state.",
    },
    "OUTPUT": {
        "role": "stdout",
        "format": "legacy_text",
        "description": "Captured NERDSS standard output.",
    },
}

RANK_SUFFIX_RE = re.compile(r"^(?P<stem>.+)_(?P<rank>[0-9]+)(?P<suffix>\.[^.]+)$")
TIME_BLOCK_RE = re.compile(r"^time\s*\(s\)\s*:", re.IGNORECASE)


def utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


def relative_path(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def base_name_for_ranked_file(name: str) -> tuple[str, int | None]:
    match = RANK_SUFFIX_RE.match(name)
    if not match:
        return name, None

    base_name = f"{match.group('stem')}{match.group('suffix')}"
    if base_name in KNOWN_FILES:
        return base_name, int(match.group("rank"))
    return name, None


def classify_file(path: Path) -> dict[str, Any]:
    base_name, rank = base_name_for_ranked_file(path.name)
    details: dict[str, Any] = dict(KNOWN_FILES.get(base_name, {}))

    suffix = path.suffix.lower()
    if not details:
        if suffix == ".pdb":
            details = {"role": "pdb_snapshot", "format": "pdb"}
        elif suffix == ".cif":
            details = {"role": "pdb_snapshot", "format": "cif"}
        elif suffix == ".xyz":
            details = {"role": "trajectory_xyz", "format": "xyz"}
        elif suffix == ".psf":
            details = {"role": "topology_psf", "format": "psf"}
        elif suffix == ".json":
            details = {"role": "other", "format": "json"}
        else:
            details = {"role": "other", "format": "unknown"}

    if rank is not None:
        details["generated_by_rank"] = rank

    return details


def read_text_prefix(path: Path, limit: int = 65536) -> str:
    with path.open("rb") as handle:
        data = handle.read(limit)
    return data.decode("utf-8", errors="replace")


def infer_columns_and_samples(path: Path, file_format: str) -> dict[str, Any]:
    if file_format in {"binary", "unknown", "restart", "pdb", "cif", "psf"}:
        return {}

    try:
        text = read_text_prefix(path)
    except OSError:
        return {}

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return {"sample_count": 0}

    metadata: dict[str, Any] = {}
    first = lines[0]
    if file_format == "tsv" or "\t" in first:
        metadata["columns"] = [part for part in first.split("\t") if part]
        metadata["sample_count"] = max(len(lines) - 1, 0)
    elif file_format == "csv" or "," in first:
        metadata["columns"] = next(csv.reader([first]))
        metadata["sample_count"] = max(len(lines) - 1, 0)
    else:
        block_count = sum(
            1
            for line in lines
            if TIME_BLOCK_RE.match(line) or line.lower().startswith("time:")
        )
        if block_count:
            metadata["sample_count"] = block_count

    return metadata


def collect_files(run_dir: Path) -> list[Path]:
    scan_roots: list[Path] = []
    if run_dir.name in {"DATA", "PDB", "RESTARTS"}:
        scan_roots.append(run_dir)
    else:
        scan_roots.append(run_dir)
        for child_name in ("DATA", "PDB", "RESTARTS"):
            child = run_dir / child_name
            if child.is_dir():
                scan_roots.append(child)

    files: list[Path] = []
    seen: set[Path] = set()
    for root in scan_roots:
        for path in sorted(root.iterdir()):
            if not path.is_file() or path in seen:
                continue
            if root == run_dir and path.name not in KNOWN_FILES:
                continue
            seen.add(path)
            files.append(path)
    return files


def collect_directories(run_dir: Path, manifest_root: Path) -> list[dict[str, Any]]:
    roles = {
        "DATA": "data",
        "PDB": "pdb_snapshots",
        "RESTARTS": "restarts",
    }
    directories: list[dict[str, Any]] = []
    for name, role in roles.items():
        path = run_dir / name
        if path.is_dir():
            directories.append(
                {
                    "path": relative_path(path, manifest_root),
                    "role": role,
                    "exists": True,
                    "file_count": sum(1 for child in path.iterdir() if child.is_file()),
                }
            )
    return directories


def build_manifest(args: argparse.Namespace) -> dict[str, Any]:
    run_dir = args.run_directory.resolve()
    manifest_root = run_dir.parent if run_dir.name in {"DATA", "PDB", "RESTARTS"} else run_dir

    files = []
    for path in collect_files(run_dir):
        entry: dict[str, Any] = {
            "path": relative_path(path, manifest_root),
            "exists": True,
            "size_bytes": path.stat().st_size,
        }
        entry.update(classify_file(path))
        entry.update(infer_columns_and_samples(path, entry.get("format", "unknown")))
        if args.hash:
            entry["sha256"] = sha256_file(path)
        files.append(entry)

    run: dict[str, Any] = {
        "output_directory": relative_path(run_dir, manifest_root),
        "working_directory": str(run_dir),
        "status": "unknown",
    }
    data_dir = run_dir / "DATA"
    if data_dir.is_dir():
        run["data_directory"] = relative_path(data_dir, manifest_root)
    if args.run_id:
        run["id"] = args.run_id
    if args.input_file:
        run["input_file"] = args.input_file

    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "manifest_type": "nerdss-run-manifest",
        "generated_at": utc_now(),
        "generator": {
            "name": "inspect_nerdss_run_manifest.py",
            "version": GENERATOR_VERSION,
            "command": sys.argv,
        },
        "run": run,
        "files": files,
        "notes": [
            "Generated from existing legacy outputs; simulation metadata may be incomplete."
        ],
    }

    directories = collect_directories(run_dir, manifest_root)
    if directories:
        manifest["directories"] = directories

    return manifest


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect a NERDSS run directory or DATA directory and emit a run manifest skeleton."
    )
    parser.add_argument(
        "run_directory",
        nargs="?",
        default=Path("."),
        type=Path,
        help="Run directory to inspect, or a DATA/PDB/RESTARTS directory.",
    )
    parser.add_argument(
        "--run-id",
        help="Optional run identifier to include in the manifest.",
    )
    parser.add_argument(
        "--input-file",
        help="Optional input file path to record in run.input_file.",
    )
    parser.add_argument(
        "--hash",
        action="store_true",
        help="Compute sha256 for each discovered file.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Write JSON to this file instead of standard output.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.run_directory.exists():
        print(f"error: path does not exist: {args.run_directory}", file=sys.stderr)
        return 2
    if not args.run_directory.is_dir():
        print(f"error: expected a directory: {args.run_directory}", file=sys.stderr)
        return 2

    manifest = build_manifest(args)
    json_text = json.dumps(
        manifest,
        indent=2 if args.pretty else None,
        sort_keys=True,
    )
    if args.pretty:
        json_text += "\n"

    if args.output:
        args.output.write_text(json_text, encoding="utf-8")
    else:
        print(json_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
