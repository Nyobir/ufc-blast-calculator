# AGENTS.md

## Commands

```bash
uv run pytest tests/ -v          # Run all tests (77 tests, ~3s)
uv run pytest tests/test_X.py -v # Run single module
uv run python -m ufc_blast gui   # Launch GUI
pyinstaller ufc_blast.spec       # Build Windows .exe
```

## Gotchas

- **Placeholder data**: `data/surface/*.csv` (8 rows each) and `data/free_air/figure_2_13_triple_point.csv` (25 rows) are synthetic placeholders. Real digitized data from UFC 3-340-02 will replace them. Do not tune algorithms against placeholder values.
- **Single computation path**: Both CLI and GUI call `compute_facade()` in `blast_params.py`. Never duplicate calculation logic in `cli.py` or `gui.py` — always route through `compute_facade()`. Tests in `test_compute_facade.py` enforce this.
- **Table2D is generic**: `Table2D` uses `family_levels`/`x_tables`/`lookup(x=, family=)` — domain-agnostic names. Call sites map domain concepts (angle, Ps0, Hc_scaled) to these generic parameters.
- **Mach stem is air-burst only**: `apply_mach_stem()` runs only when `burst_type == "air"`. Surface bursts skip it entirely.
- **PyInstaller data bundling**: When adding new data directories, update `ufc_blast.spec` `datas` list and the frozen-path logic in `blast_params.py` (`sys._MEIPASS`).
- **Coordinate convention**: `dy=0` is at burst height, ground is at `dy=-Hc`. The Mach curve uses this same frame.
- **Data pipeline**: CSV tables come from WebPlotDigitizer exports processed by `scripts/process_wpd_export.py`. Column names in CSVs must match the `from_csv()` call sites in `blast_params.py`.

## Boundaries

**Always:**
- Run `uv run pytest tests/ -v` after any change to core modules
- Use keyword arguments for `Table2D.from_csv()` and `Table2D.lookup()` calls
- Keep `compute_facade()` as the single entry point for facade calculations

**Ask first:**
- Modifying CSV data files (they are carefully digitized reference data)
- Changing the burst-type threshold (0.397) or Hopkinson-Cranz scaling
- Adding new UFC figure tables (requires matching `load_ufc_tables()` + CSV + tests)

**Never:**
- Put calculation logic directly in `cli.py` or `gui.py`
- Change Table2D interpolation method defaults without validating against UFC charts
- Commit real UFC data to public repos (copyright-restricted standard)
