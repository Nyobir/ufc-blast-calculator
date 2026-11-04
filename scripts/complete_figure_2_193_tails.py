#!/usr/bin/env python3
"""Append missing UFC Fig. 2-193 C-alpha curve tails from official DPlot.

Read the separately downloaded WBDG UFC 3-340-02 DPlot file 02_193.grf;
never commit the source GRF or replace already digitized curve segments.
Example:
    python scripts/complete_figure_2_193_tails.py \
        --grf /path/to/02_193.grf \
        --csv data/free_air/figure_2_193_calpha.csv \
        --output /path/to/proposed_calpha.csv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import defaultdict
from pathlib import Path

PSI_TO_KPA = 6.89476
SOURCE_SHA256 = "dac929c9ac9a7f25823bc47b94d7e0b54466b3c92002ab9d808f6bac6e0e0a56"


def dplot_curves(path: Path) -> dict[float, list[tuple[float, float]]]:
    if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("DPlot source checksum differs from the audited official 02_193.grf")
    lines = path.read_text().splitlines()
    if lines[:2] != ["DPLOT v1.6", "data"] or int(lines[2]) != 20:
        raise ValueError("Unexpected Fig. 2-193 DPlot format or number of curves")
    curves: dict[float, list[tuple[float, float]]] = {}
    i = 3
    for _ in range(20):
        while not lines[i].strip():
            i += 1
        n = int(lines[i]); i += 1
        points = []
        for _ in range(n):
            x, y = lines[i].split(",")
            points.append((float(x), float(y)))
            i += 1
        if not re.fullmatch(r"\s*\d+\s+0\s*", lines[i]):
            raise ValueError(f"Unexpected curve trailer at line {i + 1}")
        i += 1
        psi = float(lines[i].strip()); i += 1
        if psi in curves or any(b[0] <= a[0] for a, b in zip(points, points[1:])):
            raise ValueError(f"Duplicate pressure or unsorted angles: {psi} psi")
        if points[-1] != (90.0, 1.0):
            raise ValueError(f"Unexpected 90-degree endpoint: {psi} psi")
        curves[psi] = points
    return curves


def complete(grf: Path, source: Path, output: Path) -> dict[float, int]:
    curves = dplot_curves(grf)
    with source.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != ["ps0_kpa", "angle_deg", "calpha"]:
            raise ValueError("Unexpected CSV columns")
        rows = list(reader)
    families: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        families[float(row["ps0_kpa"])].append(row)
    if len(families) != len(curves):
        raise ValueError("DPlot and calculator have different numbers of pressure families")
    completed: list[dict[str, str]] = []
    counts: dict[float, int] = {}
    for pressure, family_rows in families.items():
        matches = [psi for psi in curves if abs(pressure - psi * PSI_TO_KPA) < 1e-7]
        if len(matches) != 1:
            raise ValueError(f"Cannot match DPlot family to calculator: {pressure} kPa")
        psi = matches[0]
        last_x = max(float(row["angle_deg"]) for row in family_rows)
        if last_x < 90.0 and any(float(b["angle_deg"]) <= float(a["angle_deg"])
                                 for a, b in zip(family_rows, family_rows[1:])):
            raise ValueError(f"Unsorted truncated calculator curve: {psi} psi")
        if last_x > 90.0:
            raise ValueError(f"Calculator curve exceeds DPlot range: {psi} psi")
        tail = [(x, y) for x, y in curves[psi] if x > last_x]
        if last_x == 90.0 and tail:
            raise ValueError("Cannot add tail after 90 degrees")
        if last_x < 90.0 and not tail:
            raise ValueError(f"No missing tail points found for {psi} psi")
        completed.extend(family_rows)
        completed.extend({"ps0_kpa": family_rows[0]["ps0_kpa"],
                          "angle_deg": format(x, ".15g"), "calpha": format(y, ".15g")}
                         for x, y in tail)
        counts[psi] = len(tail)
    if output.resolve() == source.resolve():
        raise ValueError("Output must be a separate review file; never overwrite input in place")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["ps0_kpa", "angle_deg", "calpha"])
        writer.writeheader()
        writer.writerows(completed)
    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grf", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print("Appended per pressure (psi):", complete(args.grf, args.csv, args.output))
