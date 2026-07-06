"""
export_data_package.py

Assembles the public data package for the Scientific Data Data Descriptor:
SI-normalized UFC 3-340-02 chart tables, raw digitizer exports, and frozen
derived products (benchmark facade grid, multi-distance reference curves,
benchmark anchor values, cross-chart consistency table), plus a README /
data dictionary. The package is fully usable without the calculator code.

Usage:
    python scripts/export_data_package.py [output_dir]

Default output_dir: ./data_package
"""

import csv
import math
import pathlib
import shutil
import subprocess
import sys
import warnings

import numpy as np

from ufc_blast.core import blast_params as bp
from ufc_blast.core.blast_params import (
    compute_facade,
    compute_point,
    friedlander,
    load_ufc_tables,
)

REPO = pathlib.Path(__file__).parent.parent
W, R, HC = 200.0, 30.0, 5.0
DISTANCES = [30.0, 45.0, 60.0, 75.0, 90.0]


def git_describe() -> str:
    try:
        return subprocess.run(
            ["git", "describe", "--tags", "--always"],
            cwd=REPO, capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return "unknown"


def copy_tables(out: pathlib.Path) -> tuple[int, int, int]:
    """Copy SI-normalized tables and raw digitizer exports. Returns
    (n_si_tables, n_raw_files, total_si_rows)."""
    si_dir = out / "tables_si"
    raw_dir = out / "raw_digitizer_exports"
    si_dir.mkdir(parents=True)
    raw_dir.mkdir(parents=True)

    n_si = n_raw = total_rows = 0
    for sub in ["free_air", "surface"]:
        for f in sorted((REPO / "data" / sub).glob("*.csv")):
            dest = si_dir / f"{sub}__{f.name}"
            shutil.copy2(f, dest)
            n_si += 1
            total_rows += sum(1 for _ in open(f)) - 1  # minus header
    for f in sorted((REPO / "data" / "surface" / "raw").glob("*.csv")):
        shutil.copy2(f, raw_dir / f"surface_raw__{f.name}")
        n_raw += 1
    for f in sorted((REPO / "data" / "triple_point").glob("*.csv")):
        shutil.copy2(f, raw_dir / f"triple_point__{f.name.replace(' ', '_')}")
        n_raw += 1
    return n_si, n_raw, total_rows


def export_facade(out: pathlib.Path) -> tuple[int, int]:
    """Frozen benchmark facade grid. Returns (n_points, n_flagged)."""
    fr = compute_facade(R=R, W=W, Hc=HC, width=30.0, height=24.0, step=0.25)
    mach_by_dx = dict(fr.mach_curve)
    path = out / "derived" / "facade_benchmark_200kg_30m.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    n_flagged = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "dx_m", "dy_m", "R_alpha_m", "alpha_deg", "Ps0_kPa", "C_alpha",
            "Pr_alpha_kPa", "ir_alpha_kPa_ms", "tA_ms", "t0_ms",
            "friedlander_b", "mach_zone", "extrapolation_flags",
        ])
        for gp in fr.grid_points:
            res = fr.result_map[(gp.dx, gp.dy)]
            mach_dy = mach_by_dx.get(gp.dx)
            in_mach = mach_dy is not None and gp.dy < mach_dy
            flags = ";".join(sorted(res.extrapolated))
            if flags:
                n_flagged += 1
            w.writerow([
                f"{gp.dx:.2f}", f"{gp.dy:.2f}", f"{gp.R_alpha:.4f}",
                f"{gp.alpha_deg:.4f}", f"{res.Ps0:.4f}", f"{res.C_alpha:.4f}",
                f"{res.Pr_alpha:.4f}", f"{res.ir_alpha * 1000:.4f}",
                f"{res.tA * 1000:.4f}", f"{res.t0 * 1000:.4f}",
                f"{res.b:.6f}", int(in_mach), flags,
            ])
    # Triple-point curve as its own small file
    with open(out / "derived" / "triple_point_curve_200kg_30m.csv", "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["dx_m", "triple_point_dy_m"])
        for dx, dy in fr.mach_curve:
            w.writerow([f"{dx:.2f}", f"{dy:.4f}"])
    return len(fr.grid_points), n_flagged


def export_reference_curves(out: pathlib.Path) -> None:
    """Reference parameters + sampled Friedlander histories, 30-90 m."""
    params_path = out / "derived" / "reference_curve_parameters_200kg.csv"
    hist_path = out / "derived" / "friedlander_histories_200kg.csv"
    results = {}
    with open(params_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["R_m", "Z_m_kg13", "Ps0_kPa", "C_alpha", "Pr_alpha_kPa",
                    "ir_alpha_kPa_ms", "tA_ms", "t0_ms", "friedlander_b"])
        for d in DISTANCES:
            res = compute_point(R_alpha=d, alpha_deg=0.0, W=W, burst_type="air")
            results[d] = res
            w.writerow([
                f"{d:.0f}", f"{d / W ** (1 / 3):.4f}", f"{res.Ps0:.4f}",
                f"{res.C_alpha:.4f}", f"{res.Pr_alpha:.4f}",
                f"{res.ir_alpha * 1000:.4f}", f"{res.tA * 1000:.4f}",
                f"{res.t0 * 1000:.4f}", f"{res.b:.6f}",
            ])
    t = np.linspace(0.0, 0.30, 3001)  # 0-300 ms, 0.1 ms step
    with open(hist_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["R_m", "t_ms", "P_kPa"])
        for d in DISTANCES:
            res = results[d]
            P = friedlander(t, res.Pr_alpha, res.tA, res.t0, res.b)
            for ti, Pi in zip(t, P):
                w.writerow([f"{d:.0f}", f"{ti * 1000:.1f}", f"{Pi:.4f}"])


def export_benchmark_anchor(out: pathlib.Path) -> None:
    point = compute_point(R_alpha=R, alpha_deg=0.0, W=W, burst_type="air")
    fr = compute_facade(R=R, W=W, Hc=HC, width=30.0, height=24.0, step=0.25)
    centre = fr.result_map[(0.0, 0.0)]
    mach_dy = dict(fr.mach_curve)[0.0]
    R_T = math.hypot(R, mach_dy)
    W13 = W ** (1 / 3)
    rows = [
        ("nominal_scaled_distance_Zp_m_kg13", f"{R / W13:.3f}", f"{R / W13:.3f}"),
        ("mach_evaluation_scaled_distance_ZT_m_kg13", "", f"{R_T / W13:.3f}"),
        ("triple_point_height_above_burst_m", "", f"{mach_dy:.3f}"),
        ("Ps0_kPa", f"{point.Ps0:.2f}", f"{centre.Ps0:.2f}"),
        ("C_alpha", f"{point.C_alpha:.3f}", f"{centre.C_alpha:.3f}"),
        ("Pr_alpha_kPa", f"{point.Pr_alpha:.2f}", f"{centre.Pr_alpha:.2f}"),
        ("ir_alpha_kPa_ms", f"{point.ir_alpha * 1e3:.1f}", f"{centre.ir_alpha * 1e3:.1f}"),
        ("tA_ms", f"{point.tA * 1e3:.2f}", f"{centre.tA * 1e3:.2f}"),
        ("t0_ms", f"{point.t0 * 1e3:.2f}", f"{centre.t0 * 1e3:.2f}"),
        ("friedlander_b", f"{point.b:.3f}", f"{centre.b:.3f}"),
    ]
    with open(out / "derived" / "benchmark_anchor_200kg_30m_5m.csv", "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "point_alpha0", "facade_centre_mach"])
        w.writerows(rows)


def export_cross_chart(out: pathlib.Path) -> tuple[int, int, float, float, float]:
    """Cross-route consistency table. Returns (n_total, n_used, median, p95, max)."""
    ps0_tab, ca_tab = bp._tables["ps0"], bp._tables["calpha"]
    zs, prd = np.loadtxt(REPO / "data/free_air/figure_2_7_wpd_pr.csv",
                         delimiter=",", skiprows=1, unpack=True)
    order = np.argsort(zs)
    zs, prd = zs[order], prd[order]
    val_dir = out / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)
    used = []
    with open(val_dir / "cross_chart_check.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Z_m_kg13", "Pr_direct_kPa", "Pr_indirect_kPa",
                    "rel_diff_pct", "chart_supported"])
        for z, pr in zip(zs, prd):
            ps0, oor_z = ps0_tab.lookup_flagged(float(z))
            ca, oor_fam = ca_tab.lookup_flagged(x=0.0, family=ps0)
            if oor_z or oor_fam:
                w.writerow([f"{z:.4f}", f"{pr:.4f}", "", "", 0])
                continue
            indirect = ca * ps0
            rel = abs(indirect - pr) / pr * 100.0
            used.append(rel)
            w.writerow([f"{z:.4f}", f"{pr:.4f}", f"{indirect:.4f}",
                        f"{rel:.4f}", 1])
    rel = np.array(used)
    return len(zs), len(used), float(np.median(rel)), float(
        np.percentile(rel, 95)), float(rel.max())


README = """# UFC 3-340-02 blast-loading chart data (SI-normalized) with spatial facade reference loads

Version: {version} (calculator code tag: manuscript-v1.2)
License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
Authors: Vladyslav D. Kochkarov (ORCID 0009-0002-6245-7025), Petro M. Martyniuk
Affiliation: Department of Computer Science and Applied Mathematics,
National University of Water and Environmental Engineering, Rivne, Ukraine
Contact: v.d.kochkarov@nuwm.edu.ua

## What this is

Machine-readable numerical tables of the air-blast design charts of
UFC 3-340-02, "Structures to Resist the Effects of Accidental Explosions"
(U.S. Department of Defense, 2008, with Change 2, 2014; approved for
public release, distribution unlimited), normalized to SI units, together
with derived spatial facade reference loads for a 200 kg TNT benchmark
configuration. Only numerical values extracted from the published charts
are included; no chart graphics are reproduced.

## Contents

### tables_si/  ({n_si} files, {total_rows} data rows)
SI-normalized tables. File names carry the burst regime prefix
(free_air__/surface__) and the UFC figure number they were extracted from:
- figure_2_7_*   : free-air burst parameters vs scaled distance (UFC Fig. 2-7)
- figure_2_15_*  : hemispherical surface-burst parameters (UFC Fig. 2-15)
- figure_2_193_calpha.csv : reflection coefficient C_alpha(Ps0, angle) (UFC Fig. 2-193; 20 incident-pressure families, 2,191 points)
- figure_2_194_iralpha.csv: scaled reflected impulse vs angle (UFC Fig. 2-194; 21 families, 1,096 points)
- figure_2_13_triple_point.csv : triple-point (Mach stem) height families (UFC Fig. 2-13)

Column naming convention: units are embedded in the column name.
- z_m_kg13           : scaled distance Z = R/W^(1/3) in m/kg^(1/3)
- ps0_kpa, pr_kpa    : peak incident / normally reflected overpressure, kPa
- ir_kpa_ms_kg13, is_kpa_ms_kg13 : scaled reflected / incident impulse, kPa*ms/kg^(1/3)
- tA_scaled_ms_kg13, t0_scaled_ms_kg13 : scaled arrival time / positive-phase duration, ms/kg^(1/3)
- angle_deg          : angle of incidence, degrees
- calpha             : reflection coefficient (dimensionless)
- hc_scaled_m_kg13, rg_scaled_m_kg13, ht_scaled_m_kg13 : scaled charge height,
  ground distance, and triple-point height, m/kg^(1/3)

Extraction routes (see the Data Descriptor Methods for details):
- Z-dependent curves (Figs. 2-7, 2-15) and triple-point families (Fig. 2-13):
  WebPlotDigitizer v5 digitization of the published log-log charts
  (raw exports preserved in raw_digitizer_exports/).
- Angle-dependent tables (Figs. 2-193, 2-194): transcription of a
  chart-reading engineering workbook assembled by the authors, converted
  to SI by the same processing scripts.

### raw_digitizer_exports/  ({n_raw} files)
Unmodified WebPlotDigitizer exports (imperial chart units) retained for
provenance: surface-burst curves and per-family triple-point exports.

### derived/
Frozen products computed from the tables (parameters: W = 200 kg TNT,
R = 30 m, Hc = 5 m unless stated; air burst; positive phase only):
- facade_benchmark_200kg_30m.csv : {n_facade} grid points over a
  30 m x 24 m facade resting on the ground plane (0.25 m step; dy = 0 at
  burst height, ground at dy = -5 m). Columns: geometry (dx_m, dy_m,
  R_alpha_m, alpha_deg), blast parameters (Ps0_kPa, C_alpha, Pr_alpha_kPa,
  ir_alpha_kPa_ms, tA_ms, t0_ms, friedlander_b), mach_zone (1 = below the
  per-column triple point; parameters evaluated at the triple-point height,
  uniform over height - a geometric approximation, not a Mach-front model),
  extrapolation_flags (semicolon-separated names of chart tables whose
  range was exceeded and clamped; {n_flagged} points carry the
  'triple_point' flag in columns beyond the digitized Fig. 2-13 range).
- triple_point_curve_200kg_30m.csv : per-column triple-point height.
- reference_curve_parameters_200kg.csv : normal-incidence parameters at
  R = 30, 45, 60, 75, 90 m.
- friedlander_histories_200kg.csv : sampled pressure-time histories
  P(t) = Pr*(1 - tau)*exp(-b*tau), tau = (t - tA)/t0, 0-300 ms at 0.1 ms.
- benchmark_anchor_200kg_30m_5m.csv : fully specified anchor values
  (regular reflection at alpha = 0 vs Mach-region facade centre).

### validation/
- cross_chart_check.csv : redundancy check of the two extraction routes.
  The digitized normally-reflected-pressure curve (direct route) is
  compared with C_alpha(Ps0, 0 deg) * Ps0(Z) built from independently
  extracted tables (indirect route): {n_used} of {n_cc} points are
  chart-supported; median relative difference {cc_med:.2f}%,
  95th percentile {cc_p95:.1f}%, maximum {cc_max:.1f}%.

### docs/
- unit_conversions.md : imperial-to-SI conversion factors used.

## Reuse notes

- Interpolate 1-D tables in log-log coordinates; interpolate the
  angle-dependent tables along the angle axis within bracketing
  incident-pressure families, then between families in log(Ps0).
- Do not use values outside the tabulated ranges; the derived facade file
  marks all range-clamped points in extrapolation_flags.
- The tables inherit the reading precision of the printed charts and the
  validity envelope of the underlying empirical standard; they are not a
  substitute for CFD or experiments (positive phase only; no clearing,
  shielding, or multi-reflection effects).

## Code availability

The Python calculator used to generate the derived files (110 automated
tests) is not publicly archived because blast-load estimation software is
dual-use; it is available from the corresponding author on reasonable
request for research and protective-design purposes. All files in this
package are usable without the code.

## How to cite

Please cite the Zenodo record of this dataset (DOI on the record page)
and the associated Data Descriptor, and attribute the underlying
empirical content to UFC 3-340-02 (U.S. Department of Defense, 2008).
"""


def main() -> None:
    warnings.filterwarnings("ignore")
    out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "data_package"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    load_ufc_tables()

    n_si, n_raw, total_rows = copy_tables(out)
    n_facade, n_flagged = export_facade(out)
    export_reference_curves(out)
    export_benchmark_anchor(out)
    n_cc, n_used, cc_med, cc_p95, cc_max = export_cross_chart(out)

    (out / "docs").mkdir(exist_ok=True)
    shutil.copy2(REPO / "docs" / "unit-conversions.md",
                 out / "docs" / "unit_conversions.md")

    (out / "README.md").write_text(README.format(
        version=git_describe(), n_si=n_si, n_raw=n_raw,
        total_rows=total_rows, n_facade=n_facade, n_flagged=n_flagged,
        n_cc=n_cc, n_used=n_used, cc_med=cc_med, cc_p95=cc_p95,
        cc_max=cc_max,
    ))

    print(f"Package written to {out}")
    print(f"  SI tables: {n_si} ({total_rows} data rows)")
    print(f"  raw exports: {n_raw}")
    print(f"  facade grid: {n_facade} points ({n_flagged} flagged)")
    print(f"  cross-chart: {n_used}/{n_cc} pts, median {cc_med:.2f}% "
          f"p95 {cc_p95:.2f}% max {cc_max:.2f}%")


if __name__ == "__main__":
    main()
