"""
Command-line interface for the UFC 3-340-02 blast wave calculator.

Subcommands
-----------
point   -- Compute blast parameters at a single surface point.
compute -- Compute blast parameters over a rectangular facade grid.
gui     -- Launch the interactive 3-D pressure visualiser (Task 9).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from typing import Any


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _detect_burst_type(Hc: float, W: float) -> str:
    """Detect burst type and print it."""
    from ufc_blast.core.geometry import determine_burst_type
    burst_type, scaled_hob = determine_burst_type(Hc, W)
    label = "Air burst" if burst_type == "air" else "Surface burst"
    print(f"  Burst type: {label} (Hc/W^(1/3) = {scaled_hob:.4f})")
    return burst_type


def _r_alpha_from_angle(R: float, alpha_deg: float) -> float:
    """Slant distance given perpendicular standoff and angle of incidence."""
    if alpha_deg == 0.0:
        return R
    return R / math.cos(math.radians(alpha_deg))


# ---------------------------------------------------------------------------
# Output formatters
# ---------------------------------------------------------------------------

def _print_point_result(result: Any, W: float, R: float, alpha_deg: float, R_alpha: float) -> None:
    """Print a single BlastPointResult as a formatted table."""
    sep = "-" * 52
    print(sep)
    print(f"  UFC 3-340-02 Blast Point Result")
    print(sep)
    print(f"  Charge mass W          : {W:.1f} kg TNT")
    print(f"  Standoff R             : {R:.2f} m")
    print(f"  Angle of incidence α   : {alpha_deg:.1f} °")
    print(f"  Slant distance R_α     : {R_alpha:.3f} m")
    print(sep)
    print(f"  Peak incident Ps0      : {result.Ps0:.3f} kPa")
    print(f"  Reflection coeff. Cα   : {result.C_alpha:.4f}")
    print(f"  Peak reflected Pr_α    : {result.Pr_alpha:.3f} kPa")
    print(f"  Reflected impulse ir_α : {result.ir_alpha * 1000:.3f} kPa·ms  ({result.ir_alpha:.6f} kPa·s)")
    print(f"  Arrival time tA        : {result.tA * 1000:.3f} ms")
    print(f"  Positive duration t0   : {result.t0 * 1000:.3f} ms")
    print(f"  Friedlander b          : {result.b:.4f}")
    print(sep)


def _result_to_dict(gp: Any, result: Any) -> dict[str, Any]:
    """Convert a GridPoint + BlastPointResult pair to a flat dict."""
    return {
        "dx_m": round(gp.dx, 4),
        "dy_m": round(gp.dy, 4),
        "R_alpha_m": round(gp.R_alpha, 4),
        "alpha_deg": round(gp.alpha_deg, 4),
        "Ps0_kPa": round(result.Ps0, 4),
        "C_alpha": round(result.C_alpha, 6),
        "Pr_alpha_kPa": round(result.Pr_alpha, 4),
        "ir_alpha_kPa_ms": round(result.ir_alpha * 1000.0, 4),
        "tA_ms": round(result.tA * 1000.0, 4),
        "t0_ms": round(result.t0 * 1000.0, 4),
        "b": round(result.b, 6),
    }


_TABLE_COLUMNS = [
    ("dx_m",          "dx(m)",      7),
    ("dy_m",          "dy(m)",      7),
    ("R_alpha_m",     "R_α(m)",     8),
    ("alpha_deg",     "α(°)",       7),
    ("Ps0_kPa",       "Ps0(kPa)",  10),
    ("C_alpha",       "Cα",         8),
    ("Pr_alpha_kPa",  "Pr(kPa)",   10),
    ("ir_alpha_kPa_ms", "ir(kPa·ms)", 11),
    ("tA_ms",         "tA(ms)",     9),
    ("t0_ms",         "t0(ms)",     9),
    ("b",             "b",          8),
]


def _print_table(rows: list[dict[str, Any]]) -> None:
    """Print grid results as a fixed-width column table."""
    header = "  ".join(col.rjust(width) for _, col, width in _TABLE_COLUMNS)
    sep = "-" * len(header)
    print(sep)
    print(header)
    print(sep)
    for row in rows:
        line = "  ".join(
            str(row[key]).rjust(width)
            for key, _, width in _TABLE_COLUMNS
        )
        print(line)
    print(sep)
    print(f"  {len(rows)} point(s)")


def _write_csv(rows: list[dict[str, Any]], path: str) -> None:
    """Write grid results to a CSV file."""
    fieldnames = [key for key, _, _ in _TABLE_COLUMNS]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"CSV written: {path}  ({len(rows)} rows)")


def _print_json(rows: list[dict[str, Any]]) -> None:
    """Print grid results as JSON to stdout."""
    print(json.dumps(rows, indent=2))


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def _cmd_point(args: argparse.Namespace) -> int:
    from ufc_blast.core.blast_params import compute_point, load_ufc_tables

    load_ufc_tables()
    burst_type = _detect_burst_type(args.Hc, args.W)

    R_alpha = _r_alpha_from_angle(args.R, args.alpha)
    try:
        result = compute_point(R_alpha, args.alpha, args.W, burst_type=burst_type)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    _print_point_result(result, args.W, args.R, args.alpha, R_alpha)
    return 0


def _cmd_compute(args: argparse.Namespace) -> int:
    from ufc_blast.core.blast_params import compute_point, load_ufc_tables
    from ufc_blast.core.geometry import generate_grid

    load_ufc_tables()
    burst_type = _detect_burst_type(args.Hc, args.W)

    try:
        grid = generate_grid(
            R=args.R,
            W=args.W,
            Hc=args.Hc,
            width=args.width,
            height=args.height,
            step=args.step,
        )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    rows: list[dict[str, Any]] = []
    errors = 0
    for gp in grid:
        try:
            result = compute_point(gp.R_alpha, gp.alpha_deg, args.W, burst_type=burst_type)
            rows.append(_result_to_dict(gp, result))
        except Exception as exc:
            errors += 1
            print(
                f"WARNING: skipping point (dx={gp.dx:.1f}, dy={gp.dy:.1f}): {exc}",
                file=sys.stderr,
            )

    if not rows:
        print("ERROR: no grid points could be computed.", file=sys.stderr)
        return 1

    if errors:
        print(f"WARNING: {errors} point(s) skipped due to errors.", file=sys.stderr)

    fmt = getattr(args, "format", "table")

    if args.output:
        _write_csv(rows, args.output)
    elif fmt == "json":
        _print_json(rows)
    else:
        _print_table(rows)

    return 0


def _cmd_gui(args: argparse.Namespace) -> int:
    try:
        from ufc_blast.gui import launch_gui
    except ImportError as exc:
        print(
            f"ERROR: GUI module not available ({exc}). "
            "Make sure Task 9 (gui.py) is implemented.",
            file=sys.stderr,
        )
        return 1

    launch_gui(
        W=args.W,
        R=args.R,
        Hc=args.Hc,
        width=args.width,
        height=args.height,
        step=args.step,
    )
    return 0


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ufc-blast",
        description="UFC 3-340-02 blast wave calculator for building facade analysis.",
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    # ------------------------------------------------------------------
    # point
    # ------------------------------------------------------------------
    p_point = sub.add_parser(
        "point",
        help="Compute blast parameters at a single surface point.",
        description=(
            "Compute UFC 3-340-02 reflected blast parameters at a single point "
            "defined by standoff distance R and angle of incidence alpha."
        ),
    )
    p_point.add_argument("--W", type=float, required=True, metavar="kg",
                         help="Charge mass (kg TNT equivalent).")
    p_point.add_argument("--R", type=float, required=True, metavar="m",
                         help="Perpendicular standoff distance (m).")
    p_point.add_argument("--Hc", type=float, required=True, metavar="m",
                         help="Height of burst above ground (m).")
    p_point.add_argument("--alpha", type=float, default=0.0, metavar="deg",
                         help="Angle of incidence in degrees (default: 0).")
    p_point.set_defaults(func=_cmd_point)

    # ------------------------------------------------------------------
    # compute
    # ------------------------------------------------------------------
    p_compute = sub.add_parser(
        "compute",
        help="Compute blast parameters over a rectangular facade grid.",
        description=(
            "Generate a rectangular grid of surface points and compute UFC 3-340-02 "
            "reflected blast parameters at each point."
        ),
    )
    p_compute.add_argument("--W", type=float, required=True, metavar="kg",
                           help="Charge mass (kg TNT equivalent).")
    p_compute.add_argument("--R", type=float, required=True, metavar="m",
                           help="Perpendicular standoff distance (m).")
    p_compute.add_argument("--Hc", type=float, required=True, metavar="m",
                           help="Height of burst above ground (m).")
    p_compute.add_argument("--width", type=float, default=None, metavar="m",
                           help="Facade width (m). Omit for open-field mode.")
    p_compute.add_argument("--height", type=float, default=None, metavar="m",
                           help="Facade height (m). Omit for open-field mode.")
    p_compute.add_argument("--step", type=float, default=1.0, metavar="m",
                           help="Grid spacing (m, default: 1.0).")
    p_compute.add_argument("-o", "--output", type=str, default=None, metavar="FILE",
                           help="Write results to CSV file instead of stdout.")
    p_compute.add_argument("--format", choices=["table", "json"], default="table",
                           help="Output format when writing to stdout (default: table).")
    p_compute.set_defaults(func=_cmd_compute)

    # ------------------------------------------------------------------
    # gui
    # ------------------------------------------------------------------
    p_gui = sub.add_parser(
        "gui",
        help="Launch the interactive 3-D pressure visualiser.",
        description=(
            "Open an interactive 3-D visualisation of the blast pressure "
            "distribution over the facade (requires Task 9 gui.py)."
        ),
    )
    p_gui.add_argument("--W", type=float, default=200.0, metavar="kg",
                       help="Charge mass (kg TNT equivalent, default: 200).")
    p_gui.add_argument("--R", type=float, default=30.0, metavar="m",
                       help="Perpendicular standoff distance (m, default: 30).")
    p_gui.add_argument("--Hc", type=float, default=5.0, metavar="m",
                       help="Height of burst above ground (m, default: 5).")
    p_gui.add_argument("--width", type=float, default=10.0, metavar="m",
                       help="Facade width (m, default: 10).")
    p_gui.add_argument("--height", type=float, default=8.0, metavar="m",
                       help="Facade height (m, default: 8).")
    p_gui.add_argument("--step", type=float, default=1.0, metavar="m",
                       help="Grid spacing (m, default: 1.0).")
    p_gui.set_defaults(func=_cmd_gui)

    return parser


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    # No subcommand → default to gui
    if not hasattr(args, "func"):
        args = parser.parse_args(["gui"])
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
