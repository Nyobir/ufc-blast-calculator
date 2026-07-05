# UFC 3-340-02 Blast Wave Calculator

## Quick Reference

```bash
uv run pytest tests/ -v          # 77 tests, ~3s
uv run python -m ufc_blast       # Launch GUI (default)
uv run ufc-blast point --W 200 --R 30 --Hc 5 --alpha 0   # Single point
uv run ufc-blast compute --W 200 --R 30 --Hc 5 --width 10 --height 8 --step 1 -o results.csv
uv run ufc-blast gui-check       # Headless smoke test (CI)
pyinstaller ufc_blast.spec       # Windows EXE (onedir bundle)
```

## What This Project Is

A Python calculator that implements the UFC 3-340-02 blast wave methodology for computing reflected pressures, impulses, and Friedlander waveforms on building facades. Used for:
- **CLI**: Automated validation of explosion simulations against UFC reference data
- **GUI**: Interactive visual inspection (PySide6 + matplotlib contour plots)

The project originated in `/Users/vladyslav/Documents/Phd/Cloud/Research` (session `f7051e3b`, 2026-03-26) as part of PhD research on blast simulation validation.

## Architecture — The Golden Rule

**`compute_facade()` is the single computation path.** Both CLI and GUI call it. Never duplicate calculation logic in `cli.py` or `gui.py`. Tests enforce this.

```
Input (W, R, Hc, width, height, step)
  → determine_burst_type()     # Hc/W^(1/3) > 0.397 → air; else → surface
  → generate_grid()            # GridPoint list with R_alpha, alpha_deg
  → compute_points_batch()     # Vectorized: all Ps0, Calpha, ir, tA, t0, b at once
  → apply_mach_stem()          # Air bursts only: uniform pressure below triple point
  → FacadeResult               # grid_points, result_map, mach_curve
```

## Critical Pitfalls (Learn from Past Bugs)

### Unit Conversions — The #1 Source of Errors
All UFC charts are in imperial. All internal storage is SI. See [docs/unit-conversions.md](docs/unit-conversions.md) for the full table. The most dangerous trap:

**W^(1/3) is in the DENOMINATOR**, so time/impulse conversion uses `(1/0.453592)^(1/3)` = **1.3015**, NOT `(0.453592)^(1/3)` = 0.7683. Getting this inverted inflates values by ~1.69x. The impulse factor is **8.974** (not 5.300).

### Coordinate Convention
`dy=0` is at **burst height Hc**, NOT at ground. Ground is at `dy = -Hc`. Triple point height HT is measured from ground upward, so facade coordinate `mach_dy = HT - Hc`.

### Table2D x-range Clamping Bug (Fixed)
`lookup_batch` was clamping Rg_scaled to the intersection of ALL family ranges, causing the Mach line to appear flat. Fixed to use only the two **bracketing** families' ranges. If you touch Table2D clamping logic, verify the Mach line varies per column.

### Interpolation Method Matters
- **Default: log-log** (physically correct for power-law pressure decay)
- **Use `method='linear'`** only when cross-validating against the Excel spreadsheet
- 2D tables: log on Ps0/family axis, linear on angle axis

### Cα Physics at α ≈ 40-45 deg
Cα does NOT decrease monotonically from α=0. It **peaks** in the Mach reflection regime around 40-45 deg before declining. Tests must compare α=0 vs α=70+ to reliably see a lower value.

### Other Known Issues
- `np.trapz` removed in NumPy 2.0 — use `np.trapezoid`
- PySide6 `Signal(list, dict, float)` crashes across QThread — store results as worker attributes, emit `Signal()` with no payload
- `ax.cla()` does not remove colorbars — store reference and remove explicitly
- `console=True` in PyInstaller hangs on GitHub Actions — always use `console=False`
- PyInstaller `onefile` is slow and AV-blocked — use `onedir` + zip
- `matplotlib.backends.backend_qtagg` must be in `hiddenimports`
- Column `pso_kpa` vs `ps0_kpa`: letter 'o' vs digit '0' — was caught and fixed early

## Data Sources

All data digitized from UFC 3-340-02 charts, stored as SI CSVs in `data/`. See [docs/data-sources.md](docs/data-sources.md) for details.

| Figure | File | What |
|--------|------|------|
| 2-7 | `free_air/figure_2_7_wpd_*.csv` | Ps0, tA, t0, ir, is vs Z (free-air) |
| 2-13 | `free_air/figure_2_13_triple_point.csv` | Triple point height (Mach stem) |
| 2-15 | `surface/figure_2_15_*.csv` | Ps0, tA, t0 vs Z (surface burst) |
| 2-193 | `free_air/figure_2_193_calpha.csv` | Reflection coefficient Cα(angle, Ps0) |
| 2-194 | `free_air/figure_2_194_iralpha.csv` | Reflected impulse irα(angle, Ps0) |

**Source Excel**: `/Users/vladyslav/Documents/Phd/Cloud/Корисні штуки/Таблиця UFC_Кінцевий варіант.xlsm`
**WPD raw exports**: `/Users/vladyslav/Documents/Phd/Cloud/Корисні штуки/DigitalizedUFC/2-7/`

## Validation

### Known-Good Reference Values (W=200kg, R=30m, Hc=5m, α=0 deg)
- Z = 5.13 m/kg^(1/3), Ps0 = 29.7 kPa, Cα = 2.19, Pr = 64.86 kPa
- tA = 54.6 ms, t0 = 19.8 ms, b ≈ 1.0-1.1

### Accuracy Thresholds
- Pearson r >= 0.90
- Impulse error <= 15%
- Peak pressure error <= 20%

### Cross-Validation Approach
Use `method='linear'` to match Excel exactly. Feed same inputs, check outputs match within <1% relative error. `scripts/plot_verification.py` reproduces UFC Figure 2-7 in imperial units for visual confirmation.

## Boundaries

**Always:**
- Run tests after any change to core modules
- Use keyword arguments for `Table2D.from_csv()` and `Table2D.lookup()`
- Keep `compute_facade()` as the single entry point

**Ask first:**
- Modifying CSV data files (carefully digitized reference data)
- Changing burst-type threshold (0.397) or Hopkinson-Cranz scaling
- Adding new UFC figure tables (requires `load_ufc_tables()` + CSV + tests)

**Never:**
- Put calculation logic in `cli.py` or `gui.py`
- Change interpolation defaults without validating against UFC charts
- Commit real UFC data to public repos (copyright-restricted)

## Detailed Documentation

- [docs/calculation-pipeline.md](docs/calculation-pipeline.md) — Full algorithm walkthrough in Ukrainian (10 steps)
- [docs/unit-conversions.md](docs/unit-conversions.md) — Imperial-to-SI conversion factors and traps
- [docs/data-sources.md](docs/data-sources.md) — Data provenance, Excel structure, WPD digitization details
- [docs/design-decisions.md](docs/design-decisions.md) — Why PySide6 not PyVista, why log-log, rejected approaches
- [docs/friedlander-parameters.md](docs/friedlander-parameters.md) — Waveform formulas, b-solver, confirmed test values
