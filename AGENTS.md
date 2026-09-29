# AGENTS.md

Guidance for AI coding agents (and humans) working in this repository.
Read this before changing code. User-facing documentation is in `README.md`.
`CLAUDE.md` points here, so this file is the single source of agent instructions.

## Quick reference

```bash
uv sync --extra dev                    # install with test dependencies
uv run pytest tests/ -v                # 113 tests, a few seconds
uv run pytest tests/test_X.py -v       # single module
uv run python -m ufc_blast             # launch GUI (default entry point)
uv run ufc-blast point --W 200 --R 30 --Hc 5 --alpha 0
uv run ufc-blast compute --W 200 --R 30 --Hc 5 --width 20 --height 10 --step 0.25 -o facade.csv
uv run ufc-blast batch --W 200 --R 30 --Hc 5 --alphas 0,15,30,45 --format csv
uv run ufc-blast gui-check             # headless GUI smoke test (used in CI)
pyinstaller ufc_blast.spec             # Windows bundle (onedir)
```

## What this project is

A Python implementation of the UFC 3-340-02 blast-wave methodology that computes
incident and reflected pressure, reflected impulse, arrival time, positive-phase
duration and the impulse-consistent Friedlander decay parameter `b`, at a single
point or over a facade grid. Two interfaces:

- **CLI**: scripted generation of reference loads, e.g. for validating explosion
  simulations against UFC data;
- **GUI**: interactive inspection (PySide6 + Matplotlib contour and waveform views).

It is a research and teaching tool, not a certified design tool (see `NOTICE`).

## Architecture: the golden rule

**`compute_facade()` in `ufc_blast/core/blast_params.py` is the single computation
path.** CLI and GUI both call it. Never duplicate calculation logic in `cli.py` or
`gui.py`; `tests/test_compute_facade.py` and `tests/test_gui.py` enforce this.

```
Input (W, R, Hc, width, height, step)
  → determine_burst_type()   # Hc/W^(1/3) > 0.397 → air; else → surface
  → generate_grid()          # GridPoint list with R_alpha, alpha_deg
  → compute_points_batch()   # vectorized: Ps0, C_alpha, ir, tA, t0, b for all points at once
  → apply_mach_stem()        # air bursts only: uniform pressure below the triple point
  → FacadeResult             # grid_points, result_map, mach_curve
```

- Facade grids are processed as NumPy arrays (no per-point Python loop). Keep new
  code vectorized; the scalar `compute_point()` path exists for single queries.
- `b` is solved by Brent's method for single points and by a vectorized safeguarded
  Newton iteration for grids; tests check that both agree.
- `Table1D` / `Table2D` (`core/interpolation.py`) are domain-agnostic: `Table2D`
  uses `family_levels` / `x_tables` / `lookup(x=, family=)`. Call sites map angle,
  Ps0 and Hc_scaled onto these generic names.
- Out-of-range queries are **clamped and flagged** (`extrapolated` set,
  `lookup_flagged`), never silently extrapolated. The CLI warns, the GUI hatches the
  affected region. Keep this behaviour.

## Critical pitfalls (learned from past bugs)

### Unit conversions: the #1 source of errors
All UFC charts are imperial; all internal storage is SI. Full table in
[docs/unit-conversions.md](docs/unit-conversions.md). The most dangerous trap:
**W^(1/3) is in the denominator**, so time/impulse conversion uses
`(1/0.453592)^(1/3)` = **1.3015**, not `(0.453592)^(1/3)` = 0.7683. Inverting it
inflates values by ~1.69×. The impulse factor is **8.974** (not 5.300).

### Coordinate convention
`dy = 0` is at the **burst height Hc**, not at the ground. Ground is at `dy = -Hc`.
The triple-point height HT is measured from the ground, so on the facade
`mach_dy = HT - Hc`. The Mach curve uses this same frame.

### Table2D x-range clamping (fixed bug)
`lookup_batch` used to clamp Rg_scaled to the intersection of *all* family ranges,
which made the Mach line flat. It now uses only the two **bracketing** families.
If you touch Table2D clamping, verify that the Mach line varies per column.

### Interpolation method
- Default is **log-log** (physically appropriate for power-law decay).
- Use `method='linear'` only when cross-validating against the source workbook.
- 2D tables: log on the Ps0/family axis, linear on the angle axis.

### C_alpha physics near 40–45°
C_alpha does **not** decrease monotonically from α = 0; it peaks in the Mach
reflection regime around 40–45° before declining. Tests must compare α = 0 with
α ≥ 70° to see a reliably lower value.

### Mach stem is air-burst only
`apply_mach_stem()` runs only when `burst_type == "air"`; surface bursts skip it.

### Other known issues
- `np.trapz` was removed in NumPy 2.0: use `np.trapezoid`.
- PySide6 `Signal(list, dict, float)` crashes across `QThread`: store results as
  worker attributes and emit a payload-free `Signal()`.
- `ax.cla()` does not remove colorbars: keep a reference and remove explicitly.
- PyInstaller: `console=True` hangs on GitHub Actions (use `console=False`);
  `onefile` is slow and AV-blocked (use `onedir` + zip);
  `matplotlib.backends.backend_qtagg` must be in `hiddenimports`.
- When adding data directories, update the `datas` list in `ufc_blast.spec` and the
  frozen-path logic (`sys._MEIPASS`) in `blast_params.py`.
- Column `pso_kpa` vs `ps0_kpa` (letter o vs digit 0): use `ps0_kpa`.

## Data

The 17 SI tables in `data/` were extracted from UFC 3-340-02 charts (approved for
public release, unlimited distribution; see `NOTICE`). Details and processing in
[docs/data-sources.md](docs/data-sources.md).

| Figure | File | Content |
|--------|------|---------|
| 2-7 | `free_air/figure_2_7_*.csv` | Ps0, tA, t0, ir, is vs Z (free-air burst) |
| 2-13 | `free_air/figure_2_13_triple_point.csv` | triple-point (Mach stem) height |
| 2-15 | `surface/figure_2_15_*.csv` | Ps0, tA, t0 vs Z (surface burst) |
| 2-193 | `free_air/figure_2_193_calpha.csv` | reflection coefficient C_alpha(angle, Ps0) |
| 2-194 | `free_air/figure_2_194_iralpha.csv` | reflected impulse ir_alpha(angle, Ps0) |

Processing scripts take the raw inputs as arguments; the raw WebPlotDigitizer
exports and the chart-reading workbook are kept outside the repository:

```bash
python scripts/process_wpd_export.py <dir-with-Figure-2-7-WPD-exports>
python scripts/extract_excel_data.py <chart-reading-workbook.xlsm>
python scripts/process_wpd_surface.py        # reads data/surface/raw/
python scripts/process_wpd_triple_point.py   # reads data/triple_point/
```

CSV column names always carry units (`z_m_kg13`, `ps0_kpa`, ...) and must match the
`from_csv()` call sites in `blast_params.py`.

## Validation

### Reference values (W = 200 kg, R = 30 m, Hc = 5 m, α = 0°)
- Z = 5.13 m/kg^(1/3), Ps0 = 29.67 kPa, C_alpha = 2.186, Pr = 64.86 kPa
- ir = 468.4 kPa·ms, tA = 54.6 ms, t0 = 19.8 ms, b = 1.034
- 20 m × 10 m facade at 0.25 m step: 3,321 points, centre Pr ≈ 58.3 kPa

If a change moves these numbers, it changes the published results
(tag `manuscript-v1.2`); treat that as a deliberate, documented decision.

### Cross-checks
- `scripts/cross_chart_check.py`: Pr from the reflected-pressure chart vs
  C_alpha × Ps0 over 225 points (median 0.98 %, max 3.1 %).
- `scripts/benchmark_table.py`: reproduces the manuscript benchmark table.
- `scripts/plot_verification.py`: redraws UFC Figure 2-7 in imperial units.
- Match the source workbook with `method='linear'` (expect < 1 % difference).

### Acceptance thresholds for simulation validation (downstream use)
When the calculator is used as the reference for DEM/CFD validation runs:
Pearson r ≥ 0.90, impulse error ≤ 15 %, peak pressure error ≤ 20 %.

## Boundaries

**Always**
- Run `uv run pytest tests/ -v` after any change to core modules.
- Use keyword arguments for `Table2D.from_csv()` and `Table2D.lookup()`.
- Keep `compute_facade()` as the single entry point for facade calculations.
- Keep facade code vectorized and out-of-range values flagged.

**Ask first**
- Modifying CSV data files (carefully digitized reference data).
- Changing the burst-type threshold (0.397) or Hopkinson–Cranz scaling.
- Adding new UFC figure tables (needs `load_ufc_tables()` + CSV + tests).
- Anything that changes the reference values above.

**Never**
- Put calculation logic in `cli.py` or `gui.py`.
- Change interpolation defaults without validating against the UFC charts.
- Commit chart artwork, the source workbook, credentials or local machine paths.

## Detailed documentation

- [docs/calculation-pipeline.md](docs/calculation-pipeline.md): full algorithm walkthrough (Ukrainian, 10 steps)
- [docs/unit-conversions.md](docs/unit-conversions.md): imperial-to-SI factors and traps
- [docs/data-sources.md](docs/data-sources.md): provenance, workbook layout, digitization details
- [docs/design-decisions.md](docs/design-decisions.md): why PySide6, why log-log, rejected approaches
- [docs/friedlander-parameters.md](docs/friedlander-parameters.md): waveform formulas, b-solver, test values
