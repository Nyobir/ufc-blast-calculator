# Mach Stem + Surface Burst Support — Design Spec

**Date**: 2026-03-27
**Status**: Approved
**Project**: `~/Projects/ufc_calculator`

## Problem

The calculator currently handles only free-air bursts without Mach stem effects. Two gaps:

1. **No Mach stem**: For air bursts, the triple point creates a Mach stem below a certain height on the facade. Below that line, pressure is uniform per column (vertically constant). Our contour plot shows smooth rings everywhere, which is physically incorrect for the lower portion of the facade.

2. **No surface burst support**: When `Hc/W^(1/3) <= 0.397`, the code raises an error instead of switching to surface burst tables (UFC Figure 2-15).

## Solution

### Feature 1: Mach Stem for Free-Air Bursts

Add triple point height calculation from UFC Figure 2-13 and apply a post-processing correction to the pressure grid. Visualize with a hybrid plot: contour rings above the Mach line, flat horizontal bands below.

### Feature 2: Surface Burst Support

Digitize UFC Figure 2-15 data. Auto-detect burst type from `Hc/W^(1/3)` threshold. Use surface tables when applicable. No Mach stem for surface bursts.

## Data Layer

### New CSV Tables

**Figure 2-13 — Triple Point Height** (to be digitized via WebPlotDigitizer):
- File: `data/free_air/figure_2_13_triple_point.csv`
- Columns: `hc_scaled_m_kg13`, `rg_scaled_m_kg13`, `ht_scaled_m_kg13`
- Structure: long-format 2D table (same as `figure_2_193_calpha.csv`)
- Family parameter: `Hc/W^(1/3)` (curves at 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7 in ft/lb^(1/3))
- Raw WPD data in imperial (ft/lb^(1/3)), converted to SI (m/kg^(1/3)) via `process_wpd_export.py`
- Loaded as `Table2D` with `hc_scaled` as the family parameter and `rg_scaled` as the lookup variable

**Figure 2-15 — Surface Burst Parameters** (to be digitized via WebPlotDigitizer):
- Directory: `data/surface/`
- Files:
  - `figure_2_15_wpd_ps0.csv` — columns: `z_m_kg13`, `ps0_kpa`
  - `figure_2_15_tA.csv` — columns: `z_m_kg13`, `tA_scaled_ms_kg13`
  - `figure_2_15_t0.csv` — columns: `z_m_kg13`, `t0_scaled_ms_kg13`
- Same column naming and unit conventions as existing Figure 2-7 tables
- Loaded as `Table1D`

### Table Loading Changes (blast_params.py)

New entries in `_tables` dict:
```python
_tables["triple_point"]  # Table2D: (hc_scaled, rg_scaled) → ht_scaled
_tables["surface_ps0"]   # Table1D: Z → Ps0
_tables["surface_tA"]    # Table1D: Z → tA_scaled
_tables["surface_t0"]    # Table1D: Z → t0_scaled
```

C_alpha and ir_alpha tables (Figures 2-193, 2-194) are shared across burst types — no duplication.

## Burst Type Auto-Detection

The burst type is determined automatically, never specified by the user.

**Rule**: `Hc / W^(1/3) > 0.397` → air burst, otherwise → surface burst.

**Changes**:
- `geometry.py`: remove the `ValueError` for surface bursts. Instead, return burst type information.
- `blast_params.py`: `compute_point()` and `compute_points_batch()` gain a `burst_type: str = "air"` parameter that selects which Ps0/tA/t0 tables to use. The caller determines burst type and passes it in.
- CLI: prints determined burst type in output header.
- GUI: displays burst type as a label, e.g., `"Air burst (Hc/W^⅓ = 0.86)"`.

## Mach Stem Computation

### Function: `apply_mach_stem()`

Located in `blast_params.py`. Called as a post-processing step after `compute_points_batch()`, only for air bursts.

**Signature:**
```python
def apply_mach_stem(
    grid_points: list[GridPoint],
    result_map: dict[tuple[float, float], BlastPointResult],
    W: float,
    Hc: float,
    R: float,
) -> tuple[dict[tuple[float, float], BlastPointResult], list[tuple[float, float]]]:
    """Apply Mach stem correction to grid results.

    Returns
    -------
    result_map : dict
        Modified result map with Mach-corrected values below triple point.
    mach_curve : list[tuple[float, float]]
        List of (dx, mach_dy) points defining the Mach line on the facade.
    """
```

**Algorithm:**
1. Group grid points by unique `dx` values (facade columns).
2. For each column at horizontal offset `dx`:
   a. Compute horizontal ground distance: `Rg = sqrt(R² + dx²)`
   b. Compute scaled values: `Rg_scaled = Rg / W^(1/3)`, `Hc_scaled = Hc / W^(1/3)`
   c. Look up `HT_scaled` from triple point table: `HT_scaled = _tables["triple_point"].lookup(rg_scaled=Rg_scaled, hc_scaled=Hc_scaled)`
   d. Convert to real height: `HT = HT_scaled * W^(1/3)`
   e. Convert to facade coordinate: `mach_dy = HT - Hc` (dy=0 is at burst height, ground is at dy=-Hc)
3. For the Mach-height point at `(dx, mach_dy)`: interpolate blast parameters from the two nearest grid points above and below `mach_dy` in that column.
4. For all points in the column with `dy < mach_dy`: overwrite their `BlastPointResult` with the interpolated Mach-height values.
5. Collect `(dx, mach_dy)` into the Mach curve array.

**Edge cases:**
- `Rg_scaled` outside Figure 2-13 table range → no Mach stem at that column (skip, leave values unchanged)
- `HT` below the bottom of the facade → no Mach zone visible
- `HT` above the top of the facade → entire column is Mach zone (uniform pressure)

### Coordinate Convention

```
dy = +height/2  ┌─────────────────┐  top of facade
                 │   contour rings  │
                 │                  │
dy = mach_dy     │───Mach line───── │
                 │  uniform bands   │
dy = -height/2   └─────────────────┘  bottom of facade

dy = 0           perpendicular point (at burst height Hc)
dy = -Hc         ground level
```

## Visualization

### Contour Plot (ContourCanvas.draw_contours)

**Air burst with Mach stem:**
1. Create two masked arrays from the pressure grid:
   - `Pr_above`: pressure values for points above the Mach line; `NaN` for points below
   - `Pr_below`: Mach-corrected pressure values for points below the Mach line; `NaN` for points above
2. Determine shared color range: `vmin = min(all Pr)`, `vmax = max(all Pr)`
3. Plot upper zone: `ax.contourf(X, Y, Pr_above, levels=15, cmap='YlOrRd', vmin=vmin, vmax=vmax)`
4. Plot lower zone: `ax.pcolormesh(X, Y_below, Pr_below, cmap='YlOrRd', vmin=vmin, vmax=vmax)`
5. Draw Mach curve: `ax.plot(mach_dx, mach_dy, 'k-', linewidth=2.5, label='Mach stem')`
6. Shared colorbar across both zones.

**Surface burst:**
- Pure `contourf` over the entire facade (no Mach line, no split zones).

**Air burst where Mach line is below the facade:**
- Pure `contourf` (Mach stem exists but is not visible on this facade).

### Friedlander Panel

When the clicked point is below the Mach line, show an additional label: `"(Mach zone)"` next to the point coordinates. The Friedlander curve uses the overwritten (uniform) parameters — no special handling needed.

### Burst Type Indicator

A `QLabel` in the top bar, right of the Recalculate button:
- Air burst: blue text `"Air burst (Hc/W^⅓ = 0.86)"`
- Surface burst: orange text `"Surface burst (Hc/W^⅓ = 0.31)"`

Updated on every recalculation.

## File Changes

| File | Change |
|------|--------|
| `data/free_air/figure_2_13_triple_point.csv` | **New** — WPD digitization of triple point height |
| `data/surface/figure_2_15_wpd_ps0.csv` | **New** — surface burst Ps0(Z) |
| `data/surface/figure_2_15_tA.csv` | **New** — surface burst tA(Z) |
| `data/surface/figure_2_15_t0.csv` | **New** — surface burst t0(Z) |
| `scripts/process_wpd_export.py` | **Modify** — add conversion for Figure 2-13 and 2-15 imperial→SI |
| `ufc_blast/core/blast_params.py` | **Modify** — load new tables, `burst_type` parameter, `apply_mach_stem()`, surface table selection |
| `ufc_blast/core/geometry.py` | **Modify** — remove ValueError for surface burst, expose burst type determination |
| `ufc_blast/gui.py` | **Modify** — hybrid contourf+pcolormesh, Mach line curve, burst type label, Mach zone indicator |
| `ufc_blast/cli.py` | **Modify** — print burst type in output, remove surface burst error |
| `tests/test_blast_params.py` | **Modify** — tests for surface burst table selection |
| `tests/test_mach_stem.py` | **New** — tests for triple point lookup and Mach zone application |

No changes to `interpolation.py` — `Table2D` already supports the Figure 2-13 structure.

## Dependencies

No new dependencies. All existing: numpy, scipy, matplotlib, pyside6.
