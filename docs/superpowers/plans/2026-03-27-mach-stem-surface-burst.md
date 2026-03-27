# Mach Stem + Surface Burst Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Mach stem triple-point correction for air bursts and surface burst table support to the UFC blast calculator.

**Architecture:** Post-processing Mach stem correction applied to the grid after batch computation. Burst type auto-detected from Hc/W^(1/3) threshold. Surface bursts use separate 1D tables for Ps0/tA/t0 but share the reflection coefficient tables. GUI uses hybrid contourf (above Mach line) + pcolormesh (below Mach line) visualization.

**Tech Stack:** Python 3.13+, NumPy, SciPy, matplotlib, PySide6. Package manager: uv.

---

## File Structure

| File | Responsibility |
|------|---------------|
| `data/free_air/figure_2_13_triple_point.csv` | **New** — Triple point height 2D table (hc_scaled, rg_scaled → ht_scaled) |
| `data/surface/figure_2_15_wpd_ps0.csv` | **New** — Surface burst incident pressure (Z → Ps0) |
| `data/surface/figure_2_15_tA.csv` | **New** — Surface burst arrival time (Z → tA_scaled) |
| `data/surface/figure_2_15_t0.csv` | **New** — Surface burst positive phase duration (Z → t0_scaled) |
| `scripts/process_wpd_export.py` | **Modify** — Add Figure 2-13 and 2-15 processing functions |
| `ufc_blast/core/blast_params.py` | **Modify** — Load new tables, add `burst_type` parameter, add `apply_mach_stem()` |
| `ufc_blast/core/geometry.py` | **Modify** — Replace ValueError with burst type determination helper |
| `ufc_blast/cli.py` | **Modify** — Print burst type, remove surface burst error |
| `ufc_blast/gui.py` | **Modify** — Hybrid contour+pcolormesh, Mach curve, burst type label |
| `tests/test_geometry.py` | **Modify** — Update surface burst tests from "raises" to "returns type" |
| `tests/test_blast_params.py` | **Modify** — Add surface burst table selection tests |
| `tests/test_mach_stem.py` | **New** — Triple point lookup, Mach zone application, edge cases |

---

### Task 1: Burst Type Detection in geometry.py

Replace the ValueError for surface bursts with a helper function that returns the burst type string.

**Files:**
- Modify: `ufc_blast/core/geometry.py:137-144`
- Modify: `tests/test_geometry.py:134-153`

- [ ] **Step 1: Write tests for the new burst type helper**

Add to `tests/test_geometry.py`, replacing the existing `TestGenerateGridApplicability` class:

```python
from ufc_blast.core.geometry import GridPoint, compute_point_geometry, generate_grid, determine_burst_type


class TestDetermineBurstType:
    """Tests for burst type auto-detection."""

    def test_air_burst(self):
        """Hc/W^(1/3) > 0.397 → 'air'."""
        burst_type, scaled_hob = determine_burst_type(Hc=10.0, W=1.0)
        assert burst_type == "air"
        assert scaled_hob > 0.397

    def test_surface_burst(self):
        """Hc/W^(1/3) <= 0.397 → 'surface'."""
        burst_type, scaled_hob = determine_burst_type(Hc=3.0, W=1000.0)
        assert burst_type == "surface"
        assert scaled_hob <= 0.397

    def test_exact_threshold_is_surface(self):
        """Ratio exactly equal to 0.397 → 'surface' (strict inequality)."""
        W = 8.0
        Hc = 0.397 * (W ** (1.0 / 3.0))
        burst_type, _ = determine_burst_type(Hc=Hc, W=W)
        assert burst_type == "surface"

    def test_zero_hc_is_surface(self):
        """Hc=0 → surface burst."""
        burst_type, scaled_hob = determine_burst_type(Hc=0.0, W=100.0)
        assert burst_type == "surface"
        assert scaled_hob == 0.0


class TestGenerateGridApplicability:
    """generate_grid now works for both air and surface bursts."""

    def test_surface_burst_no_longer_raises(self):
        """Surface burst parameters should not raise ValueError."""
        points = generate_grid(R=30.0, W=1000.0, Hc=3.0, width=10.0, height=8.0, step=2.0)
        assert len(points) > 0

    def test_valid_hob_still_works(self):
        """Air burst parameters continue to work."""
        points = generate_grid(R=30.0, W=1.0, Hc=10.0, width=4.0, height=4.0)
        assert len(points) > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/test_geometry.py::TestDetermineBurstType -v`
Expected: FAIL with `ImportError: cannot import name 'determine_burst_type'`

- [ ] **Step 3: Implement determine_burst_type and remove ValueError from generate_grid**

In `ufc_blast/core/geometry.py`, add the new function after `compute_point_geometry` (before `generate_grid`):

```python
def determine_burst_type(Hc: float, W: float) -> tuple[str, float]:
    """Determine burst type from scaled height-of-burst.

    Parameters
    ----------
    Hc : float
        Height of burst above ground (m).
    W : float
        Charge mass (kg TNT equivalent).

    Returns
    -------
    burst_type : str
        ``'air'`` if Hc/W^(1/3) > 0.397, ``'surface'`` otherwise.
    scaled_hob : float
        The computed Hc/W^(1/3) value.
    """
    scaled_hob = Hc / (W ** (1.0 / 3.0)) if W > 0 else 0.0
    burst_type = "air" if scaled_hob > 0.397 else "surface"
    return burst_type, scaled_hob
```

In `generate_grid`, replace lines 137-144 (the applicability check block):

```python
    # --- Applicability check -------------------------------------------------
    scaled_hob = Hc / (W ** (1.0 / 3.0))
    if scaled_hob <= 0.397:
        raise ValueError(
            f"Scaled height-of-burst Hc/W^(1/3) = {scaled_hob:.4f} does not "
            f"exceed 0.397 — charge is at or below surface-burst threshold. "
            f"UFC 3-340-02 air-burst charts are not applicable."
        )
```

With (simply remove the block — no replacement needed, delete these lines entirely). The function should go straight from the docstring to `points: list[GridPoint] = []`.

- [ ] **Step 4: Run all geometry tests**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/test_geometry.py -v`
Expected: All pass (the old `test_surface_burst_raises` and `test_exact_threshold_raises` tests must be removed/replaced by the new class)

- [ ] **Step 5: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add ufc_blast/core/geometry.py tests/test_geometry.py
git commit -m "feat: replace surface burst ValueError with determine_burst_type helper"
```

---

### Task 2: Surface Burst Data Tables (placeholder CSVs for testing)

Create minimal placeholder CSV files so the loading code can be tested before the real WPD digitization is done. Real data will replace these later.

**Files:**
- Create: `data/surface/figure_2_15_wpd_ps0.csv`
- Create: `data/surface/figure_2_15_tA.csv`
- Create: `data/surface/figure_2_15_t0.csv`

- [ ] **Step 1: Create surface data directory and placeholder CSVs**

```bash
mkdir -p /Users/vladyslav/Projects/ufc_calculator/data/surface
```

Create `data/surface/figure_2_15_wpd_ps0.csv`:
```csv
z_m_kg13,ps0_kpa
0.200000,30000.000000
0.500000,5000.000000
1.000000,1200.000000
2.000000,300.000000
5.000000,50.000000
10.000000,12.000000
20.000000,3.500000
40.000000,1.200000
```

Create `data/surface/figure_2_15_tA.csv`:
```csv
z_m_kg13,tA_scaled_ms_kg13
0.200000,0.050000
0.500000,0.150000
1.000000,0.400000
2.000000,1.200000
5.000000,4.500000
10.000000,12.000000
20.000000,30.000000
40.000000,70.000000
```

Create `data/surface/figure_2_15_t0.csv`:
```csv
z_m_kg13,t0_scaled_ms_kg13
0.200000,0.100000
0.500000,0.350000
1.000000,0.800000
2.000000,1.800000
5.000000,4.000000
10.000000,8.000000
20.000000,14.000000
40.000000,22.000000
```

These are physically plausible placeholder values (surface burst pressures are roughly 1.8× free-air at same Z). They will be replaced with real WPD-digitized data from Figure 2-15 later.

- [ ] **Step 2: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add data/surface/
git commit -m "data: add placeholder surface burst CSVs for Figure 2-15"
```

---

### Task 3: Triple Point Placeholder CSV

Create a placeholder 2D CSV for the triple point height table (Figure 2-13) so Mach stem code can be tested.

**Files:**
- Create: `data/free_air/figure_2_13_triple_point.csv`

- [ ] **Step 1: Create the placeholder CSV**

Create `data/free_air/figure_2_13_triple_point.csv` with the same long-format structure as `figure_2_193_calpha.csv`. The family parameter is `hc_scaled_m_kg13` and the lookup variable is `rg_scaled_m_kg13`:

```csv
hc_scaled_m_kg13,rg_scaled_m_kg13,ht_scaled_m_kg13
0.397,1.200,0.050
0.397,2.000,0.200
0.397,3.000,0.500
0.397,4.000,1.000
0.397,6.000,2.000
0.595,1.500,0.030
0.595,2.500,0.150
0.595,3.500,0.400
0.595,5.000,1.000
0.595,7.000,2.200
0.794,2.000,0.020
0.794,3.000,0.100
0.794,4.500,0.350
0.794,6.000,0.900
0.794,8.000,2.000
1.191,2.500,0.010
1.191,4.000,0.080
1.191,5.500,0.250
1.191,7.000,0.700
1.191,9.000,1.600
1.588,3.500,0.005
1.588,5.000,0.050
1.588,6.500,0.180
1.588,8.000,0.500
1.588,10.000,1.200
```

Note: `hc_scaled` values are the imperial values (1, 1.5, 2, 3, 4 ft/lb^(1/3)) converted to SI via `× 0.39685`. These are placeholder values — real WPD data will replace them.

- [ ] **Step 2: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add data/free_air/figure_2_13_triple_point.csv
git commit -m "data: add placeholder triple point height CSV for Figure 2-13"
```

---

### Task 4: Load Surface Burst and Triple Point Tables in blast_params.py

Extend `load_ufc_tables()` to load the new CSV files and add `burst_type` parameter to `compute_point()` and `compute_points_batch()`.

**Files:**
- Modify: `ufc_blast/core/blast_params.py:23-102, 105-219`
- Modify: `tests/test_blast_params.py`

- [ ] **Step 1: Write tests for burst_type parameter and surface table loading**

Add to `tests/test_blast_params.py`:

```python
from ufc_blast.core.blast_params import (
    BlastPointResult,
    compute_point,
    compute_points_batch,
    friedlander,
    load_ufc_tables,
    solve_friedlander_b,
)


class TestSurfaceBurstComputation:
    """Tests for surface burst table selection."""

    def test_surface_burst_returns_result(self):
        """compute_point with burst_type='surface' returns a valid result."""
        result = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0, burst_type="surface")
        assert isinstance(result, BlastPointResult)
        assert result.Ps0 > 0.0
        assert result.Pr_alpha > 0.0
        assert result.tA > 0.0
        assert result.t0 > 0.0
        assert result.b > 0.0

    def test_surface_uses_different_tables(self):
        """Surface burst should give different Ps0 than air burst at same Z."""
        result_air = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0, burst_type="air")
        result_surface = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0, burst_type="surface")
        # Different tables → different Ps0 (placeholder data differs from air data)
        assert result_air.Ps0 != result_surface.Ps0

    def test_batch_surface_burst(self):
        """compute_points_batch with burst_type='surface' works."""
        import numpy as np
        R_alphas = np.array([30.0, 40.0, 50.0])
        alpha_degs = np.array([0.0, 10.0, 20.0])
        results = compute_points_batch(R_alphas, alpha_degs, W=200.0, burst_type="surface")
        assert len(results) == 3
        for r in results:
            assert isinstance(r, BlastPointResult)
            assert r.Ps0 > 0.0

    def test_default_burst_type_is_air(self):
        """Omitting burst_type defaults to air burst (backward compatible)."""
        result_default = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0)
        result_air = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0, burst_type="air")
        assert result_default.Ps0 == result_air.Ps0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/test_blast_params.py::TestSurfaceBurstComputation -v`
Expected: FAIL — `compute_point()` doesn't accept `burst_type` parameter yet

- [ ] **Step 3: Add SURFACE_DATA_DIR and load new tables**

In `ufc_blast/core/blast_params.py`, after line 30 (`DATA_DIR = _BASE / "data" / "free_air"`), add:

```python
SURFACE_DATA_DIR = _BASE / "data" / "surface"
```

In `load_ufc_tables()`, add after the existing table loads (after line 102):

```python
    # Surface burst tables (Figure 2-15)
    _tables["surface_ps0"] = Table1D.from_csv(
        SURFACE_DATA_DIR / "figure_2_15_wpd_ps0.csv",
        "z_m_kg13",
        "ps0_kpa",
    )
    _tables["surface_tA"] = Table1D.from_csv(
        SURFACE_DATA_DIR / "figure_2_15_tA.csv",
        "z_m_kg13",
        "tA_scaled_ms_kg13",
    )
    _tables["surface_t0"] = Table1D.from_csv(
        SURFACE_DATA_DIR / "figure_2_15_t0.csv",
        "z_m_kg13",
        "t0_scaled_ms_kg13",
    )

    # Triple point height table (Figure 2-13)
    _tables["triple_point"] = Table2D.from_csv(
        DATA_DIR / "figure_2_13_triple_point.csv",
        "hc_scaled_m_kg13",
        "rg_scaled_m_kg13",
        "ht_scaled_m_kg13",
    )
```

- [ ] **Step 4: Add burst_type parameter to compute_point()**

Change the signature of `compute_point` from:

```python
def compute_point(R_alpha: float, alpha_deg: float, W: float) -> BlastPointResult:
```

To:

```python
def compute_point(R_alpha: float, alpha_deg: float, W: float, burst_type: str = "air") -> BlastPointResult:
```

Replace the three 1D lookups (lines 129, 140, 143) with table selection:

```python
    # Select tables based on burst type
    ps0_table = _tables["surface_ps0"] if burst_type == "surface" else _tables["ps0"]
    tA_table = _tables["surface_tA"] if burst_type == "surface" else _tables["tA"]
    t0_table = _tables["surface_t0"] if burst_type == "surface" else _tables["t0"]

    # Incident overpressure
    Ps0 = ps0_table.lookup(Z)

    # Reflection coefficient from Figure 2-193 (shared across burst types)
    C_alpha = _tables["calpha"].lookup(angle=alpha_deg, ps0=Ps0)
    Pr_alpha = C_alpha * Ps0

    # Reflected specific impulse from Figure 2-194 (kPa·ms/kg^1/3 → kPa·s)
    ir_alpha_scaled = _tables["iralpha"].lookup(angle=alpha_deg, ps0=Ps0)
    ir_alpha = ir_alpha_scaled * W_cbrt / 1000.0

    # Arrival time (ms/kg^1/3 → s)
    tA = tA_table.lookup(Z) * W_cbrt / 1000.0

    # Positive phase duration (ms/kg^1/3 → s)
    t0 = t0_table.lookup(Z) * W_cbrt / 1000.0
```

- [ ] **Step 5: Add burst_type parameter to compute_points_batch()**

Change the signature from:

```python
def compute_points_batch(
    R_alphas: np.ndarray,
    alpha_degs: np.ndarray,
    W: float,
) -> list[BlastPointResult]:
```

To:

```python
def compute_points_batch(
    R_alphas: np.ndarray,
    alpha_degs: np.ndarray,
    W: float,
    burst_type: str = "air",
) -> list[BlastPointResult]:
```

Replace the three bulk 1D lookups (lines 194-196) with:

```python
    # Select tables based on burst type
    ps0_table = _tables["surface_ps0"] if burst_type == "surface" else _tables["ps0"]
    tA_table = _tables["surface_tA"] if burst_type == "surface" else _tables["tA"]
    t0_table = _tables["surface_t0"] if burst_type == "surface" else _tables["t0"]

    # Bulk 1D lookups
    Ps0 = ps0_table.lookup_batch(Z)
    tA_vals = tA_table.lookup_batch(Z) * W_cbrt / 1000.0
    t0_vals = t0_table.lookup_batch(Z) * W_cbrt / 1000.0
```

- [ ] **Step 6: Run all tests**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/ -v`
Expected: All pass

- [ ] **Step 7: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add ufc_blast/core/blast_params.py tests/test_blast_params.py
git commit -m "feat: load surface burst + triple point tables, add burst_type parameter"
```

---

### Task 5: Mach Stem apply_mach_stem() Function

Implement the post-processing function that corrects pressure values below the triple point height.

**Files:**
- Modify: `ufc_blast/core/blast_params.py`
- Create: `tests/test_mach_stem.py`

- [ ] **Step 1: Write tests for apply_mach_stem()**

Create `tests/test_mach_stem.py`:

```python
"""Tests for Mach stem triple point correction."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ufc_blast.core.blast_params import (
    BlastPointResult,
    apply_mach_stem,
    compute_points_batch,
    load_ufc_tables,
)
from ufc_blast.core.geometry import GridPoint, generate_grid


@pytest.fixture(autouse=True)
def _load_tables():
    load_ufc_tables()


class TestApplyMachStem:
    """Tests for the Mach stem post-processing correction."""

    def _make_grid_and_results(self, W=200.0, R=30.0, Hc=5.0, width=10.0, height=8.0, step=1.0):
        """Helper: generate grid and compute results."""
        grid_points = generate_grid(R=R, W=W, Hc=Hc, width=width, height=height, step=step)
        R_alphas = np.array([gp.R_alpha for gp in grid_points])
        alpha_degs = np.array([gp.alpha_deg for gp in grid_points])
        results = compute_points_batch(R_alphas, alpha_degs, W)
        result_map = {}
        for gp, res in zip(grid_points, results):
            key = (gp.dx, gp.dy)
            if key not in result_map:
                result_map[key] = res
        return grid_points, result_map

    def test_returns_result_map_and_mach_curve(self):
        """apply_mach_stem returns a dict and a list of (dx, dy) tuples."""
        grid_points, result_map = self._make_grid_and_results()
        new_map, mach_curve = apply_mach_stem(grid_points, result_map, W=200.0, Hc=5.0, R=30.0)
        assert isinstance(new_map, dict)
        assert isinstance(mach_curve, list)
        # Mach curve entries are (dx, mach_dy) tuples
        for entry in mach_curve:
            assert len(entry) == 2

    def test_mach_zone_has_uniform_pressure_per_column(self):
        """Below the Mach line, all points in a column should have the same Pr_alpha."""
        grid_points, result_map = self._make_grid_and_results()
        new_map, mach_curve = apply_mach_stem(grid_points, result_map, W=200.0, Hc=5.0, R=30.0)

        if not mach_curve:
            pytest.skip("No Mach stem visible for these parameters")

        # Group by dx
        from collections import defaultdict
        columns = defaultdict(list)
        for gp in grid_points:
            columns[gp.dx].append(gp.dy)

        # For columns with Mach curve, check uniformity below
        mach_dict = {dx: mach_dy for dx, mach_dy in mach_curve}
        for dx, mach_dy in mach_dict.items():
            if dx not in columns:
                continue
            below_prs = []
            for dy in columns[dx]:
                if dy < mach_dy:
                    below_prs.append(new_map[(dx, dy)].Pr_alpha)
            if len(below_prs) > 1:
                # All should be equal (uniform pressure in Mach zone)
                for pr in below_prs:
                    assert abs(pr - below_prs[0]) < 1e-10, (
                        f"Non-uniform Pr in Mach zone at dx={dx}: {below_prs}"
                    )

    def test_above_mach_line_unchanged(self):
        """Points above the Mach line should be unchanged."""
        grid_points, result_map = self._make_grid_and_results()
        original_map = {k: v for k, v in result_map.items()}
        new_map, mach_curve = apply_mach_stem(grid_points, result_map, W=200.0, Hc=5.0, R=30.0)

        if not mach_curve:
            pytest.skip("No Mach stem visible for these parameters")

        mach_dict = {dx: mach_dy for dx, mach_dy in mach_curve}
        for (dx, dy), res in new_map.items():
            if dx in mach_dict and dy >= mach_dict[dx]:
                orig = original_map[(dx, dy)]
                assert res.Pr_alpha == orig.Pr_alpha, (
                    f"Point ({dx}, {dy}) above Mach line was modified"
                )

    def test_no_mach_stem_when_out_of_table_range(self):
        """If Rg_scaled is outside Figure 2-13 range, return empty Mach curve."""
        # Very large R means Rg_scaled may exceed table bounds
        grid_points, result_map = self._make_grid_and_results(R=500.0, width=2.0, height=2.0)
        new_map, mach_curve = apply_mach_stem(grid_points, result_map, W=200.0, Hc=5.0, R=500.0)
        # Either empty or the values are unchanged
        assert isinstance(mach_curve, list)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/test_mach_stem.py -v`
Expected: FAIL — `cannot import name 'apply_mach_stem'`

- [ ] **Step 3: Implement apply_mach_stem()**

Add to `ufc_blast/core/blast_params.py`, after the `friedlander()` function (at end of file):

```python
def apply_mach_stem(
    grid_points: list,
    result_map: dict[tuple[float, float], BlastPointResult],
    W: float,
    Hc: float,
    R: float,
) -> tuple[dict[tuple[float, float], BlastPointResult], list[tuple[float, float]]]:
    """Apply Mach stem correction to grid results.

    For each facade column (unique dx), compute the triple point height from
    UFC Figure 2-13. All grid points below that height get their blast
    parameters overwritten with the values at the triple point height
    (uniform Mach stem pressure).

    Parameters
    ----------
    grid_points : list[GridPoint]
        Grid points on the facade.
    result_map : dict
        Mapping (dx, dy) → BlastPointResult (modified in-place).
    W : float
        Charge mass (kg TNT equivalent).
    Hc : float
        Height of burst above ground (m).
    R : float
        Perpendicular standoff distance (m).

    Returns
    -------
    result_map : dict
        Modified result map with Mach-corrected values below triple point.
    mach_curve : list[tuple[float, float]]
        List of (dx, mach_dy) points defining the Mach line on the facade.
    """
    if not _tables:
        load_ufc_tables()

    triple_table = _tables["triple_point"]
    W_cbrt = W ** (1.0 / 3.0)
    Hc_scaled = Hc / W_cbrt

    # Group grid points by column (dx)
    from collections import defaultdict
    columns: dict[float, list[float]] = defaultdict(list)
    for gp in grid_points:
        columns[gp.dx].append(gp.dy)

    # Sort dy values in each column
    for dx in columns:
        columns[dx].sort()

    mach_curve: list[tuple[float, float]] = []

    for dx in sorted(columns.keys()):
        dy_list = columns[dx]

        # Horizontal ground distance from charge epicenter
        Rg = math.sqrt(R ** 2 + dx ** 2)
        Rg_scaled = Rg / W_cbrt

        # Look up triple point height
        try:
            HT_scaled = triple_table.lookup(angle=Rg_scaled, ps0=Hc_scaled)
        except (ValueError, KeyError):
            # Rg_scaled or Hc_scaled outside table range — no Mach stem here
            continue

        HT = HT_scaled * W_cbrt
        # Convert to facade coordinate: dy=0 is at burst height, ground at dy=-Hc
        mach_dy = HT - Hc

        mach_curve.append((dx, mach_dy))

        # Find the blast result at the Mach height by interpolating between
        # the two nearest grid points above and below mach_dy
        dy_arr = np.array(dy_list)
        above_mask = dy_arr >= mach_dy
        below_mask = dy_arr < mach_dy

        if not np.any(below_mask):
            # Mach line is below all grid points — no correction needed
            continue

        if not np.any(above_mask):
            # Entire column is in Mach zone — use the topmost point
            top_dy = dy_list[-1]
            mach_result = result_map[(dx, top_dy)]
        else:
            # Interpolate: find nearest point at or above mach_dy
            above_dy = float(dy_arr[above_mask][0])  # smallest dy >= mach_dy
            mach_result = result_map[(dx, above_dy)]

        # Overwrite all points below the Mach line
        for dy in dy_list:
            if dy < mach_dy:
                result_map[(dx, dy)] = mach_result

    return result_map, mach_curve
```

- [ ] **Step 4: Run tests**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/test_mach_stem.py -v`
Expected: All pass

- [ ] **Step 5: Run full test suite**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/ -v`
Expected: All pass

- [ ] **Step 6: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add ufc_blast/core/blast_params.py tests/test_mach_stem.py
git commit -m "feat: implement apply_mach_stem() Mach stem correction"
```

---

### Task 6: Update CLI to Show Burst Type

Update the CLI to auto-detect burst type, pass it to compute functions, and display it in output.

**Files:**
- Modify: `ufc_blast/cli.py:25-38, 141-206`

- [ ] **Step 1: Update _check_hob to return burst type info**

In `ufc_blast/cli.py`, replace the `_check_hob` function (lines 29-38):

```python
def _check_hob(Hc: float, W: float) -> None:
    """Print a warning (but do not abort) when the scaled HOB is borderline."""
    hob = _scaled_hob(Hc, W)
    if hob <= 0.397:
        print(
            f"WARNING: Scaled height-of-burst Hc/W^(1/3) = {hob:.4f} does not "
            f"exceed 0.397 — charge is at or below surface-burst threshold. "
            f"UFC 3-340-02 air-burst charts are not applicable.",
            file=sys.stderr,
        )
```

With:

```python
def _detect_burst_type(Hc: float, W: float) -> str:
    """Detect burst type and print it."""
    from ufc_blast.core.geometry import determine_burst_type
    burst_type, scaled_hob = determine_burst_type(Hc, W)
    label = "Air burst" if burst_type == "air" else "Surface burst"
    print(f"  Burst type: {label} (Hc/W^(1/3) = {scaled_hob:.4f})")
    return burst_type
```

- [ ] **Step 2: Update _cmd_point to use burst type**

In `_cmd_point` (around line 141), replace:

```python
    load_ufc_tables()
    _check_hob(args.Hc, args.W)

    R_alpha = _r_alpha_from_angle(args.R, args.alpha)
    try:
        result = compute_point(R_alpha, args.alpha, args.W)
```

With:

```python
    load_ufc_tables()
    burst_type = _detect_burst_type(args.Hc, args.W)

    R_alpha = _r_alpha_from_angle(args.R, args.alpha)
    try:
        result = compute_point(R_alpha, args.alpha, args.W, burst_type=burst_type)
```

- [ ] **Step 3: Update _cmd_compute to use burst type**

In `_cmd_compute` (around line 158), replace:

```python
    load_ufc_tables()
```

With:

```python
    load_ufc_tables()
    burst_type = _detect_burst_type(args.Hc, args.W)
```

And in the loop where `compute_point` is called (around line 181):

```python
            result = compute_point(gp.R_alpha, gp.alpha_deg, args.W)
```

Replace with:

```python
            result = compute_point(gp.R_alpha, gp.alpha_deg, args.W, burst_type=burst_type)
```

- [ ] **Step 4: Verify CLI works**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run ufc-blast point --W 200 --R 30 --Hc 5 --alpha 0`
Expected: Output includes `Burst type: Air burst (Hc/W^(1/3) = 0.8550)` and normal results.

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run ufc-blast point --W 1000 --R 30 --Hc 3 --alpha 0`
Expected: Output includes `Burst type: Surface burst (Hc/W^(1/3) = 0.3000)` and results from surface tables.

- [ ] **Step 5: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add ufc_blast/cli.py
git commit -m "feat: auto-detect burst type in CLI, display in output"
```

---

### Task 7: Update GUI — Burst Type Label and Mach Stem in Worker

Add burst type indicator in the top bar and integrate Mach stem into the background computation worker.

**Files:**
- Modify: `ufc_blast/gui.py:61-114, 462-618`

- [ ] **Step 1: Update _compute_results to accept burst_type**

In `gui.py`, change the `_compute_results` function (line 61):

```python
def _compute_results(
    grid_points: list[GridPoint],
    W: float,
) -> dict[tuple[float, float], BlastPointResult]:
```

To:

```python
def _compute_results(
    grid_points: list[GridPoint],
    W: float,
    burst_type: str = "air",
) -> dict[tuple[float, float], BlastPointResult]:
```

And change line 71:

```python
    results = compute_points_batch(R_alphas, alpha_degs, W)
```

To:

```python
    results = compute_points_batch(R_alphas, alpha_degs, W, burst_type=burst_type)
```

- [ ] **Step 2: Update _ComputeWorker to handle burst type and Mach stem**

Add imports at top of `gui.py` (add `apply_mach_stem` to the import from blast_params):

```python
from ufc_blast.core.blast_params import (
    BlastPointResult,
    apply_mach_stem,
    compute_point,
    compute_points_batch,
    friedlander,
    load_ufc_tables,
)
from ufc_blast.core.geometry import GridPoint, determine_burst_type, generate_grid
```

Update the `_ComputeWorker` class. Add new attributes and modify `run()`:

```python
class _ComputeWorker(QThread):
    """Background thread for blast computation."""

    finished = Signal()

    def __init__(self, R, W, Hc, width, height, step):
        super().__init__()
        self.R = R
        self.W = W
        self.Hc = Hc
        self.width = width
        self.height = height
        self.step = step
        self.error: str | None = None
        self.grid_points: list | None = None
        self.result_map: dict | None = None
        self.mach_curve: list | None = None
        self.burst_type: str = "air"
        self.scaled_hob: float = 0.0
        self.elapsed: float = 0.0

    def run(self):
        t0 = time.perf_counter()
        try:
            self.burst_type, self.scaled_hob = determine_burst_type(self.Hc, self.W)
            self.grid_points = generate_grid(
                self.R, self.W, self.Hc,
                width=self.width, height=self.height, step=self.step,
            )
            self.result_map = _compute_results(self.grid_points, self.W, burst_type=self.burst_type)

            # Apply Mach stem correction for air bursts only
            if self.burst_type == "air":
                self.result_map, self.mach_curve = apply_mach_stem(
                    self.grid_points, self.result_map,
                    W=self.W, Hc=self.Hc, R=self.R,
                )
            else:
                self.mach_curve = []

            self.elapsed = time.perf_counter() - t0
        except Exception as exc:
            self.error = str(exc)
        self.finished.emit()
```

- [ ] **Step 3: Add burst type label to input bar**

In `BlastWindow._build_input_bar()`, after the Recalculate button (after line 565):

```python
        self._recalc_btn = QPushButton("Recalculate")
        self._recalc_btn.setFixedWidth(110)
        self._recalc_btn.clicked.connect(self._recalculate)
        layout.addWidget(self._recalc_btn)

        # Burst type indicator
        self._burst_type_label = QLabel("")
        self._burst_type_label.setStyleSheet("font-weight: bold; padding-left: 10px;")
        layout.addWidget(self._burst_type_label)

        layout.addStretch()
```

Remove the existing `layout.addStretch()` that was there before.

- [ ] **Step 4: Update _on_compute_done to pass Mach data and set burst label**

Replace the `_on_compute_done` method:

```python
    def _on_compute_done(self) -> None:
        """Called when background computation finishes."""
        self._recalc_btn.setEnabled(True)

        if self._worker.error:
            self._status.showMessage(f"Error: {self._worker.error}")
            return

        self._grid_points = self._worker.grid_points
        self._result_map  = self._worker.result_map
        self._mach_curve  = self._worker.mach_curve
        self._width       = self._worker_width
        self._height      = self._worker_height

        # Update burst type label
        bt = self._worker.burst_type
        hob = self._worker.scaled_hob
        if bt == "air":
            self._burst_type_label.setText(f"Air burst (Hc/W^⅓ = {hob:.2f})")
            self._burst_type_label.setStyleSheet("font-weight: bold; padding-left: 10px; color: #2060c0;")
        else:
            self._burst_type_label.setText(f"Surface burst (Hc/W^⅓ = {hob:.2f})")
            self._burst_type_label.setStyleSheet("font-weight: bold; padding-left: 10px; color: #c06020;")

        self._contour_canvas.draw_contours(
            self._grid_points, self._result_map, self._width, self._height,
            mach_curve=self._mach_curve,
        )
        self._status.showMessage(
            f"Calculated {len(self._grid_points)} points in {self._worker.elapsed:.2f} s"
        )
```

- [ ] **Step 5: Run tests to verify nothing is broken**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/ -v`
Expected: All pass

- [ ] **Step 6: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add ufc_blast/gui.py
git commit -m "feat: burst type label + Mach stem in background worker"
```

---

### Task 8: Hybrid Contour + Pcolormesh Visualization

Update `ContourCanvas.draw_contours()` to render the Mach zone with pcolormesh and draw the Mach curve.

**Files:**
- Modify: `ufc_blast/gui.py:189-246`

- [ ] **Step 1: Update draw_contours() signature**

Change the signature from:

```python
    def draw_contours(
        self,
        grid_points: list[GridPoint],
        result_map: dict[tuple[float, float], BlastPointResult],
        width: float,
        height: float,
    ) -> None:
```

To:

```python
    def draw_contours(
        self,
        grid_points: list[GridPoint],
        result_map: dict[tuple[float, float], BlastPointResult],
        width: float,
        height: float,
        mach_curve: list[tuple[float, float]] | None = None,
    ) -> None:
```

- [ ] **Step 2: Implement hybrid rendering**

Replace the body of `draw_contours()` with:

```python
        """Redraw the contour map with new data, with optional Mach stem zones."""
        self._grid_points = grid_points
        self._result_map = result_map
        self._click_marker = None

        import matplotlib.patches as mpatches

        X, Y, Pr = _build_meshgrid(grid_points, result_map)

        fig = self.figure
        fig.clear()
        self.ax = fig.add_subplot(111)
        ax = self.ax
        self._colorbar = None

        vmin, vmax = float(np.nanmin(Pr)), float(np.nanmax(Pr))

        has_mach = mach_curve is not None and len(mach_curve) > 0

        if has_mach:
            # Build a per-cell Mach height lookup
            mach_dict = {dx: mach_dy for dx, mach_dy in mach_curve}

            # Create masked arrays: above and below the Mach line
            Pr_above = Pr.copy()
            Pr_below = Pr.copy()

            dx_vals = sorted({gp.dx for gp in grid_points})
            dy_vals = sorted({gp.dy for gp in grid_points})

            for j, dy in enumerate(dy_vals):
                for i, dx in enumerate(dx_vals):
                    if dx in mach_dict:
                        if dy < mach_dict[dx]:
                            Pr_above[j, i] = np.nan  # mask from contourf
                        else:
                            Pr_below[j, i] = np.nan  # mask from pcolormesh

            # Upper zone: contourf (rings)
            cf = ax.contourf(X, Y, Pr_above, levels=15, cmap="YlOrRd",
                             vmin=vmin, vmax=vmax, zorder=1)
            cs = ax.contour(X, Y, Pr_above, levels=10, colors="black",
                            linewidths=0.8, zorder=2)
            ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f kPa")

            # Lower zone: pcolormesh (flat bands)
            ax.pcolormesh(X, Y, Pr_below, cmap="YlOrRd",
                          vmin=vmin, vmax=vmax, zorder=1, shading="nearest")

            # Mach curve
            mach_xs = [dx for dx, _ in mach_curve]
            mach_ys = [dy for _, dy in mach_curve]
            ax.plot(mach_xs, mach_ys, "k-", linewidth=2.5, zorder=3, label="Mach stem")

            # Colorbar from contourf
            self._colorbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
            self._colorbar.set_label("Pr_α (kPa)", fontsize=9)
        else:
            # No Mach stem: pure contourf (surface burst or Mach below facade)
            cf = ax.contourf(X, Y, Pr, levels=15, cmap="YlOrRd", zorder=1)
            cs = ax.contour(X, Y, Pr, levels=10, colors="black", linewidths=0.8, zorder=2)
            ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f kPa")
            self._colorbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
            self._colorbar.set_label("Pr_α (kPa)", fontsize=9)

        # Building outline
        half_w = width / 2.0
        half_h = height / 2.0
        ax.add_patch(mpatches.Rectangle(
            (-half_w, -half_h), width, height,
            linewidth=2, edgecolor="black", facecolor="none", zorder=4,
        ))

        # Perpendicular centre marker
        ax.plot(0.0, 0.0, "ko", markersize=6, zorder=5)

        # Clip view to building bounds with small margin
        margin = max(width, height) * 0.1
        ax.set_xlim(-half_w - margin, half_w + margin)
        ax.set_ylim(-half_h - margin, half_h + margin)

        ax.set_xlabel("Horizontal offset (m)", fontsize=9)
        ax.set_ylabel("Vertical offset (m)", fontsize=9)
        ax.set_title("Reflected peak overpressure Pr_α (kPa)", fontsize=10)
        ax.set_aspect("equal", adjustable="box")

        fig.tight_layout()
        self.draw()
```

- [ ] **Step 3: Add Mach zone indicator to Friedlander panel**

In `FriedlanderPanel._update_labels()`, change the coordinate label line to accept an optional `mach_zone` flag. Update the `update_point` method:

```python
    def update_point(self, gp: GridPoint, res: BlastPointResult, mach_zone: bool = False) -> None:
        """Redraw the waveform and update all parameter labels."""
        self._draw_waveform(gp, res)
        self._update_labels(gp, res, mach_zone=mach_zone)
        self._label_placeholder.hide()
        self._param_widget.show()
```

And in `_update_labels`:

```python
    def _update_labels(self, gp: GridPoint, res: BlastPointResult, mach_zone: bool = False) -> None:
        coord_text = f"Point:  dx = {gp.dx:.2f} m,  dy = {gp.dy:.2f} m"
        if mach_zone:
            coord_text += "  (Mach zone)"
        self._param_labels["coord"].setText(coord_text)
```

(rest of `_update_labels` unchanged)

- [ ] **Step 4: Update _on_point_selected to pass Mach zone flag**

In `BlastWindow._on_point_selected`:

```python
    def _on_point_selected(
        self,
        gp: GridPoint,
        res: BlastPointResult,
    ) -> None:
        # Check if point is in Mach zone
        mach_zone = False
        if hasattr(self, '_mach_curve') and self._mach_curve:
            mach_dict = {dx: mach_dy for dx, mach_dy in self._mach_curve}
            if gp.dx in mach_dict and gp.dy < mach_dict[gp.dx]:
                mach_zone = True
        self._friedlander_panel.update_point(gp, res, mach_zone=mach_zone)
        self._status.showMessage(
            f"Selected point  dx={gp.dx:.1f} m  dy={gp.dy:.1f} m  "
            f"Pr_α={res.Pr_alpha:.1f} kPa"
            + ("  [Mach zone]" if mach_zone else "")
        )
```

- [ ] **Step 5: Test the GUI manually**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run ufc-blast gui`
Expected:
- Blue "Air burst" label appears in top bar
- Contour plot renders with Mach line (if visible for default params)
- Clicking below Mach line shows "(Mach zone)" in Friedlander panel

- [ ] **Step 6: Run full test suite**

Run: `cd /Users/vladyslav/Projects/ufc_calculator && uv run pytest tests/ -v`
Expected: All pass

- [ ] **Step 7: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add ufc_blast/gui.py
git commit -m "feat: hybrid contour+pcolormesh Mach zone visualization"
```

---

### Task 9: Update process_wpd_export.py for Figure 2-13 and 2-15

Add processing functions for the new data once the user digitizes it from UFC.

**Files:**
- Modify: `scripts/process_wpd_export.py`

- [ ] **Step 1: Add Figure 2-13 processing function**

Add to `scripts/process_wpd_export.py`, after `process_all()`:

```python
def process_figure_2_13(src_dir: Path, out_dir: Path) -> None:
    """Process Figure 2-13 triple point height data.

    Expects one WPD CSV per Hc/W^(1/3) curve, named like:
      hc_1.0.csv, hc_1.5.csv, hc_2.0.csv, ...

    Each file has two columns: Rg_scaled (ft/lb^1/3), HT_scaled (ft/lb^1/3).
    Both axes use the Z_FACTOR conversion (length/mass^1/3).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "figure_2_13_triple_point.csv"

    all_rows = []
    for csv_file in sorted(src_dir.glob("hc_*.csv")):
        # Parse Hc_scaled from filename: hc_2.5.csv → 2.5
        hc_imp = float(csv_file.stem.split("_", 1)[1])
        hc_si = hc_imp * Z_FACTOR  # ft/lb^(1/3) → m/kg^(1/3)

        raw = read_wpd_csv(csv_file)
        for rg_imp, ht_imp in raw:
            rg_si = rg_imp * Z_FACTOR
            ht_si = ht_imp * Z_FACTOR
            all_rows.append((hc_si, rg_si, ht_si))

    # Write long-format CSV
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["hc_scaled_m_kg13", "rg_scaled_m_kg13", "ht_scaled_m_kg13"])
        for hc, rg, ht in sorted(all_rows):
            writer.writerow([f"{hc:.6f}", f"{rg:.6f}", f"{ht:.6f}"])

    print(f"Figure 2-13: {len(all_rows)} rows → {out_path}")


def process_figure_2_15(src_dir: Path, out_dir: Path) -> None:
    """Process Figure 2-15 surface burst data.

    Same structure as Figure 2-7 processing. Expects:
      Ps0.csv, tA.csv, t0.csv

    Units: same as Figure 2-7 (Z in ft/lb^1/3, pressure in psi, time in ms/lb^1/3).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Ps0
    raw = read_wpd_csv(src_dir / "Ps0.csv")
    converted = [(z * Z_FACTOR, v * KPA_PER_PSI) for z, v in raw]
    write_csv(out_dir / "figure_2_15_wpd_ps0.csv", ["z_m_kg13", "ps0_kpa"], converted)
    print(f"Surface Ps0: {len(converted)} rows → figure_2_15_wpd_ps0.csv")

    # tA
    raw = read_wpd_csv(src_dir / "tA.csv")
    converted = [(z * Z_FACTOR, v * TIME_FACTOR) for z, v in raw]
    write_csv(out_dir / "figure_2_15_tA.csv", ["z_m_kg13", "tA_scaled_ms_kg13"], converted)
    print(f"Surface tA:  {len(converted)} rows → figure_2_15_tA.csv")

    # t0
    raw = read_wpd_csv(src_dir / "t0.csv")
    converted = [(z * Z_FACTOR, v * TIME_FACTOR) for z, v in raw]
    write_csv(out_dir / "figure_2_15_t0.csv", ["z_m_kg13", "t0_scaled_ms_kg13"], converted)
    print(f"Surface t0:  {len(converted)} rows → figure_2_15_t0.csv")
```

- [ ] **Step 2: Update the `if __name__` block to call new processors when directories exist**

Replace the existing `if __name__ == "__main__":` block:

```python
if __name__ == "__main__":
    process_all()

    # Figure 2-13 (triple point) — process if source directory exists
    fig_2_13_src = Path("/Users/vladyslav/Documents/Phd/Cloud/Корисні штуки/DigitalizedUFC/2-13")
    if fig_2_13_src.exists():
        process_figure_2_13(fig_2_13_src, OUT_DIR)

    # Figure 2-15 (surface burst) — process if source directory exists
    fig_2_15_src = Path("/Users/vladyslav/Documents/Phd/Cloud/Корисні штуки/DigitalizedUFC/2-15")
    surface_out = Path("/Users/vladyslav/Projects/ufc_calculator/data/surface")
    if fig_2_15_src.exists():
        process_figure_2_15(fig_2_15_src, surface_out)
```

- [ ] **Step 3: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add scripts/process_wpd_export.py
git commit -m "feat: add WPD processing for Figure 2-13 and 2-15"
```

---

### Task 10: Update PyInstaller Spec for New Data

Ensure the Windows .exe build includes the new `data/surface/` directory.

**Files:**
- Modify: `ufc_blast.spec`

- [ ] **Step 1: Update the spec to include surface burst data**

In `ufc_blast.spec`, after the existing `data_files` line:

```python
data_files = [(str(f), str(Path("data") / "free_air")) for f in DATA_DIR.glob("*.csv")]
```

Add:

```python
SURFACE_DIR = ROOT / "data" / "surface"
if SURFACE_DIR.exists():
    data_files += [(str(f), str(Path("data") / "surface")) for f in SURFACE_DIR.glob("*.csv")]
```

- [ ] **Step 2: Commit**

```bash
cd /Users/vladyslav/Projects/ufc_calculator
git add ufc_blast.spec
git commit -m "build: include surface burst data in PyInstaller bundle"
```

---

## Self-Review Checklist

**Spec coverage:**
- [x] Figure 2-13 triple point table → Task 3 (placeholder), Task 9 (WPD processing)
- [x] Figure 2-15 surface burst tables → Task 2 (placeholder), Task 9 (WPD processing)
- [x] Table loading → Task 4
- [x] Burst type auto-detection → Task 1
- [x] `burst_type` parameter on compute_point/compute_points_batch → Task 4
- [x] `apply_mach_stem()` → Task 5
- [x] CLI burst type display → Task 6
- [x] GUI burst type label → Task 7
- [x] Hybrid contourf+pcolormesh → Task 8
- [x] Mach curve drawing → Task 8
- [x] Friedlander panel Mach zone indicator → Task 8
- [x] PyInstaller update → Task 10
- [x] process_wpd_export.py → Task 9

**Placeholder scan:** No TBD/TODO found. All code blocks are complete.

**Type consistency:**
- `determine_burst_type()` returns `tuple[str, float]` — used consistently in Tasks 1, 6, 7
- `apply_mach_stem()` returns `tuple[dict, list[tuple[float, float]]]` — used consistently in Tasks 5, 7, 8
- `burst_type: str = "air"` parameter — consistent across Tasks 4, 6, 7
- `mach_curve: list[tuple[float, float]]` — consistent across Tasks 5, 7, 8
- `Table2D.lookup(angle=..., ps0=...)` parameter names match existing code — used in Task 5
