#!/usr/bin/env python3
"""Inspect a legacy NERDSS .inp file and emit schema-oriented JSON.

This is not a replacement parser for NERDSS. It is a small, stdlib-only helper
for Phase 4 validation samples: segment blocks, preserve reaction equations,
parse simple key/value lines, and optionally inspect adjacent .mol files.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "0.1.0"


def strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def normalise_key(key: str) -> str:
    aliases = {
        "nitr": "nItr",
        "timestep": "timeStep",
        "timewrite": "timeWrite",
        "trajwrite": "trajWrite",
        "restartwrite": "restartWrite",
        "pdbwrite": "pdbWrite",
        "transitionwrite": "transitionWrite",
        "checkpoint": "checkPoint",
        "fromrestart": "fromRestart",
        "assocdissocwrite": "assocDissocWrite",
        "clusteroverlapcheck": "clusterOverlapCheck",
        "rngwrite": "rngWrite",
        "overlapseplimit": "overlapSepLimit",
        "scalemaxdisplace": "scaleMaxDisplace",
        "waterbox": "waterBox",
        "implicitlipid": "implicitLipid",
        "xbctype": "xBCtype",
        "ybctype": "yBCtype",
        "zbctype": "zBCtype",
        "issphere": "isSphere",
        "spherer": "sphereR",
        "hascompartment": "hasCompartment",
        "compartmentr": "compartmentR",
        "compartmentsited": "compartmentSiteD",
        "compartmentsiterho": "compartmentSiteRho",
    }
    compact = re.sub(r"\s+", "", key).lower()
    return aliases.get(compact, key.strip())


def parse_scalar(value: str) -> Any:
    text = value.strip()
    lowered = text.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    if lowered in {"nan", "null"}:
        return None
    if re.fullmatch(r"[+-]?\d+", text):
        return int(text)
    if re.fullmatch(r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?", text):
        return float(text)
    return text


def parse_value(value: str) -> Any:
    text = strip_comment(value)
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        if not inner:
            return []
        return [parse_scalar(part) for part in inner.split(",")]
    return parse_scalar(text)


def split_key_value(line: str) -> tuple[str, Any] | None:
    clean = strip_comment(line)
    if "=" not in clean:
        return None
    key, value = clean.split("=", 1)
    return normalise_key(key), parse_value(value)


def collect_blocks(path: Path) -> dict[str, list[str]]:
    blocks: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            clean = strip_comment(raw)
            lowered = re.sub(r"\s+", "", clean).lower()
            if lowered.startswith("start") and len(lowered) > len("start"):
                current = lowered[len("start") :]
                blocks.setdefault(current, [])
                continue
            if lowered.startswith("end") and current == lowered[len("end") :]:
                current = None
                continue
            if current is not None:
                blocks[current].append(raw.rstrip("\n"))
    return blocks


def parse_key_value_block(lines: list[str]) -> dict[str, Any]:
    values: dict[str, Any] = {}
    for line in lines:
        parsed = split_key_value(line)
        if parsed is None:
            continue
        key, value = parsed
        values[key] = value
    return values


def parse_molecule_decl(line: str) -> dict[str, Any] | None:
    clean = strip_comment(line)
    if not clean:
        return None
    if ":" in clean:
        name, rhs = clean.split(":", 1)
        molecule: dict[str, Any] = {
            "name": name.strip(),
            "moleculeFile": f"{name.strip()}.mol",
            "legacy": {"declaration": clean},
        }
        copy_match = re.search(r"\d+", rhs)
        if copy_match:
            molecule["copies"] = int(copy_match.group(0))
        else:
            molecule["legacy"]["copyExpression"] = rhs.strip()
        return molecule
    name = clean.strip()
    return {
        "name": name,
        "moleculeFile": f"{name}.mol",
        "legacy": {"declaration": clean},
    }


def parse_mol_file(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None

    properties: dict[str, Any] = {}
    coordinates: list[dict[str, Any]] = []
    bonds: list[list[str]] = []
    in_bonds = False

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            clean = strip_comment(raw)
            if not clean:
                continue

            if in_bonds:
                parts = clean.split()
                if len(parts) >= 2:
                    bonds.append([parts[0], parts[1]])
                    continue
                in_bonds = False

            parts = clean.split()
            if len(parts) == 4:
                maybe_numbers = [parse_scalar(part) for part in parts[1:]]
                if all(isinstance(value, (int, float)) for value in maybe_numbers):
                    coordinates.append(
                        {"site": parts[0], "position": [float(value) for value in maybe_numbers]}
                    )
                    continue

            parsed = split_key_value(clean)
            if parsed is None:
                if clean.lower().startswith("bonds"):
                    in_bonds = True
                continue
            key, value = parsed
            properties[key] = value
            if key.lower() == "bonds":
                in_bonds = True

    template: dict[str, Any] = {}
    if properties:
        template["properties"] = properties
    if coordinates:
        template["coordinates"] = coordinates
    if bonds:
        template["bonds"] = bonds
    return template or None


def parse_reactions(lines: list[str]) -> list[dict[str, Any]]:
    reactions: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    for raw in lines:
        clean = strip_comment(raw)
        if not clean:
            continue
        is_equation = "->" in clean and "=" not in clean.split("->", 1)[0]
        if is_equation:
            current = {
                "equation": clean,
                "reversible": "<->" in clean,
                "parameters": {},
                "legacy": {"rawLines": [raw]},
            }
            reactions.append(current)
            continue
        parsed = split_key_value(clean)
        if parsed is not None and current is not None:
            key, value = parsed
            current["parameters"][key] = value
            current["legacy"]["rawLines"].append(raw)
        elif current is not None:
            current["legacy"]["rawLines"].append(raw)

    return reactions


def build_document(path: Path, include_mol_files: bool) -> dict[str, Any]:
    blocks = collect_blocks(path)
    base_dir = path.parent
    molecules = [
        molecule
        for molecule in (parse_molecule_decl(line) for line in blocks.get("molecules", []))
        if molecule is not None
    ]

    if include_mol_files:
        for molecule in molecules:
            mol_file = base_dir / str(molecule["moleculeFile"])
            template = parse_mol_file(mol_file)
            if template is not None:
                molecule["template"] = template

    simulation: dict[str, Any] = {
        "parameters": parse_key_value_block(blocks.get("parameters", [])),
        "molecules": molecules,
        "reactions": parse_reactions(blocks.get("reactions", [])),
    }
    boundaries = parse_key_value_block(blocks.get("boundaries", []))
    if boundaries:
        simulation["boundaries"] = boundaries
    observables = [strip_comment(line) for line in blocks.get("observables", []) if strip_comment(line)]
    if observables:
        simulation["observables"] = observables

    return {
        "format": "nerdss-input-json",
        "schemaVersion": SCHEMA_VERSION,
        "source": {
            "kind": "legacy-inspection",
            "path": str(path),
            "generatedBy": Path(__file__).name,
        },
        "compatibility": {
            "legacyInputPath": str(path),
            "legacyWorkingDirectory": str(base_dir),
            "legacyPreserved": True,
            "notes": [
                "Reaction equations are preserved as legacy strings.",
                "This helper is an inspector for validation samples, not a complete parser rewrite.",
            ],
        },
        "simulation": simulation,
    }


def list_blocks(path: Path) -> None:
    blocks = collect_blocks(path)
    for name in sorted(blocks):
        non_empty = sum(1 for line in blocks[name] if strip_comment(line))
        print(f"{name}: {non_empty} non-empty line(s)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect a legacy NERDSS .inp file and emit schema-oriented JSON."
    )
    parser.add_argument("input", type=Path, help="Path to a legacy .inp file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write JSON to this path instead of stdout.",
    )
    parser.add_argument(
        "--include-mol-files",
        action="store_true",
        help="Inspect adjacent legacy .mol files named by the molecules block.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List discovered legacy blocks and exit without emitting JSON.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.list:
        list_blocks(args.input)
        return 0

    document = build_document(args.input, args.include_mol_files)
    payload = json.dumps(document, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(payload + "\n", encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
