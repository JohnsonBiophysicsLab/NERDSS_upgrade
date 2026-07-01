#!/usr/bin/env python3
"""Aggregate an Instruments Time Profiler export into per-function timings.

Input is the XML produced by:
  xcrun xctrace export --input <trace> \
    --xpath '/trace-toc/run[@number="1"]/data/table[@schema="time-profile"]'

It sums sampling weight by leaf frame (self/exclusive time) and by every frame
in the backtrace (inclusive time), resolving the frame/binary/weight `ref` reuse
that xctrace uses to keep the XML small.

Usage:
  tools/profile_aggregate.py <time-profile.xml> [--top N]
"""
from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
from collections import defaultdict


def analyze(path: str, top: int, total_seconds: float | None) -> None:
    root = ET.parse(path).getroot()
    frame_name: dict[str, str] = {}
    frame_binid: dict[str, str | None] = {}
    bin_name: dict[str, str] = {}
    weight_val: dict[str, int] = {}
    self_t: dict[str, int] = defaultdict(int)
    incl_t: dict[str, int] = defaultdict(int)
    self_bin: dict[str, int] = defaultdict(int)
    total = 0
    nsamp = 0

    for row in root.iter("row"):
        w_el = row.find("weight")
        bt = row.find("backtrace")
        if w_el is None or bt is None:
            continue
        ref = w_el.get("ref")
        if ref is not None:
            w = weight_val.get(ref, 0)
        else:
            w = int(w_el.text or 0)
            if w_el.get("id"):
                weight_val[w_el.get("id")] = w

        names: list[str] = []
        bins: list[str] = []
        for f in bt.findall("frame"):
            fr = f.get("ref")
            if fr is not None:
                nm = frame_name.get(fr, "?")
                binid = frame_binid.get(fr)
            else:
                fid = f.get("id")
                nm = f.get("name", "?")
                b = f.find("binary")
                binid = None
                if b is not None:
                    bref = b.get("ref")
                    if bref is not None:
                        binid = bref
                    else:
                        binid = b.get("id")
                        bin_name[binid] = b.get("name", "?")
                if fid:
                    frame_name[fid] = nm
                    frame_binid[fid] = binid
            names.append(nm)
            bins.append(bin_name.get(binid, "?"))

        if not names:
            continue
        total += w
        nsamp += 1
        self_t[names[0]] += w
        self_bin[bins[0]] += w
        for nm in set(names):
            incl_t[nm] += w

    if total == 0:
        print("No samples found. Is this a time-profile export?")
        return

    # Real seconds per function. Sampling weight sums to on-CPU sample time,
    # which under-counts true CPU time; if a measured CPU total is supplied,
    # scale each share to it so the seconds column reflects real wall/CPU time.
    if total_seconds is not None:
        secs = lambda v: total_seconds * v / total  # noqa: E731
        unit = f"(scaled to {total_seconds:.2f}s measured CPU)"
    else:
        secs = lambda v: v / 1e9  # noqa: E731  (weight is ns)
        unit = "(sampling-derived; use --total-seconds for real CPU time)"

    print(f"# samples={nsamp}  on-CPU sample time={total/1e9:.2f}s  seconds {unit}")

    def show(title, table):
        print(f"\n## {title} (top {top})")
        print(f"{'sec':>8}  {'pct':>6}  symbol")
        for nm, v in sorted(table.items(), key=lambda x: -x[1])[:top]:
            print(f"{secs(v):8.3f}  {100*v/total:5.1f}%  {nm[:90]}")

    print("\n## SELF TIME by binary")
    print(f"{'sec':>8}  {'pct':>6}  binary")
    for nm, v in sorted(self_bin.items(), key=lambda x: -x[1])[:10]:
        print(f"{secs(v):8.3f}  {100*v/total:5.1f}%  {nm}")

    show("SELF TIME by function", self_t)
    show("INCLUSIVE TIME by function", incl_t)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("xml", help="time-profile XML exported by xctrace")
    ap.add_argument("--top", type=int, default=25, help="rows to print (default 25)")
    ap.add_argument(
        "--total-seconds",
        type=float,
        default=None,
        help="measured CPU seconds (e.g. from /usr/bin/time); scales the seconds column to real CPU time",
    )
    args = ap.parse_args()
    analyze(args.xml, args.top, args.total_seconds)


if __name__ == "__main__":
    main()
