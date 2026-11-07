# AGENTS.md

Guidance for AI coding agents (and humans) working in this repository.
Read this before changing code. User-facing documentation is in `README.md`.
`CLAUDE.md` points here, so this file is the single source of agent instructions.

## Quick reference

```bash
uv sync --locked --extra dev                           # core/CLI development environment
uv run --locked --extra dev pytest tests/ -q           # GUI tests may skip without PySide6
uv sync --locked --extra dev --extra gui               # optional desktop GUI
QT_QPA_PLATFORM=offscreen uv run --locked --extra dev --extra gui pytest tests/ -q
                                                        # 122 passed on reviewed Pi/Mac candidate
uv run --locked --extra dev pytest tests/test_X.py -v # single module
uv run --locked --extra gui python -m ufc_blast        # launch GUI
uv run --locked --extra dev ufc-blast point --W 200 --R 30 --Hc 5 --alpha 0
uv run --locked --extra dev ufc-blast compute --W 200 --R 30 --Hc 5 --width 20 --height 10 --step 0.25 -o facade.csv
uv run --locked --extra dev ufc-blast batch --W 200 --R 30 --Hc 5 --alphas 0,15,30,45 --format csv
QT_QPA_PLATFORM=offscreen uv run --locked --extra gui ufc-blast gui-check
uv build --no-sources                                   # sdist and wheel
uv run --locked --extra dev python scripts/check_release_wheel.py dist/ufc_blast-1.3.0-py3-none-any.whl
pyinstaller ufc_blast.spec                              # Windows onedir bundle
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

**`compute_facade()` in `ufc_blast/core/blast_params.py` is the shared facade
computation path.** CLI and GUI both call it. `compute_point()` is a separate
regular-reflection query at a specified angle and does not apply the facade's
Mach-region correction. Never duplicate calculation logic in `cli.py` or
`gui.py`; `tests/test_compute_facade.py` and `tests/test_gui.py` enforce this.

```
Input (W, R, Hc, width, height, step)
  → determine_burst_type()   # Hc/W^(1/3) > 0.397 → air; else → surface
  → generate_grid()          # GridPoint list with R_alpha, alpha_deg
  → compute_points_batch()   # vectorized: Ps0, C_alpha, ir, tA, t0, b for all points at once
  → apply_mach_stem()        # air bursts only: below the triple point, UFC §2-13.2 Mach-front values
                             # computed at the ground (compute_mach_front, Figs. 2-9/2-10 → 2-7 → 2-193/2-194)
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
- The 0.397 scaled-height switch is a chart-selection convention, not a physical
  discontinuity; scaled heights 0.3–0.5 receive a sensitivity warning.
- `t0_fictitious` is an informational flag for the Mach-front duration adjustment;
  chart-limit flags and `friedlander_b` still need attention.

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
  frozen-path logic (`sys._MEIPASS`) in `blast_params.py`, plus the wheel payload
  in `pyproject.toml`. Source, wheel and frozen executables resolve data differently.
- Column `pso_kpa` vs `ps0_kpa` (letter o vs digit 0): use `ps0_kpa`.
- On macOS the wheel checker can fail its temporary import-path assertion if
  `TMPDIR` is the symlinked `/var/folders/...`. Re-run the unchanged checker
  with `TMPDIR` set to the real path and inspect independent CSV hashes and
  installed-wheel smoke results before diagnosing a package-data failure.

## Data

The 19 SI-normalized CSV tables in `data/free_air/` (13) and `data/surface/`
(6) derive from UFC 3-340-02 charts (approved for public release, unlimited
distribution; see `NOTICE`). Wheels carry these tables under `ufc_blast/data/`;
raw inputs are not bundled. Details in [docs/data-sources.md](docs/data-sources.md).

| Figure | File | Content |
|--------|------|---------|
| 2-7 | `free_air/figure_2_7_*.csv` | Ps0, tA, t0, ir, is vs Z (free-air burst) |
| 2-9, 2-10 | `free_air/figure_2_9_pra.csv`, `figure_2_10_ira.csv` | ground-evaluated Mach-front pressure and impulse |
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

### Reference values (W = 200 kg, R = 30 m, Hc = 5 m)
- `ufc-blast point --alpha 0` is a **regular-reflection point query**:
  Z = 5.13 m/kg^(1/3), Ps0 = 29.67 kPa, C_alpha = 2.186,
  Pr = 64.864 kPa, ir = 0.468400 kPa·s, tA = 54.598 ms,
  t0 = 19.818 ms, b = 1.0341.
- A 20 m × 10 m facade at 0.25 m step has 3,321 points. Its centre
  `(dx=0, dy=0)` lies in the **ground-evaluated Mach region** and gives
  Pr = 110.081271 kPa, ir = 0.622702163 kPa·s, with no chart-limit flags.
  Do not substitute the point-command value for the corrected facade centre.

The old `manuscript-v1.2` tag describes the previous method, not this
revision. The SoftwareX numerical audit is anchored to commit `9d287e0`;
the reviewed packaging/GUI candidate is `1d69fe4`. Future numerical changes
require explicit documentation and revalidation.

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
- Run `uv run --locked --extra dev --extra gui pytest tests/ -q` after any
  change to core modules (use `QT_QPA_PLATFORM=offscreen` on headless hosts).
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

- [docs/calculation-pipeline.md](docs/calculation-pipeline.md): current algorithm walkthrough (Ukrainian)
- [docs/unit-conversions.md](docs/unit-conversions.md): imperial-to-SI factors and traps
- [docs/data-sources.md](docs/data-sources.md): provenance, workbook layout, digitization details
- [docs/design-decisions.md](docs/design-decisions.md): why PySide6, why log-log, rejected approaches
- [docs/friedlander-parameters.md](docs/friedlander-parameters.md): waveform formulas, solver limits and distinct point/facade checks
