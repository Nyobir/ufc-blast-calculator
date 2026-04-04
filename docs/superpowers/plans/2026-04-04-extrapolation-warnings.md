# Extrapolation Warnings Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Flag every blast result point that required extrapolation beyond UFC table data, surface the flag in CLI and GUI so validation pipelines can exclude untrustworthy values.

**Architecture:** Add `lookup_flagged` / `lookup_batch_flagged` methods to `Table2D` that return `(value, out_of_range_bool)` alongside the computed result. Thread a per-point `extrapolated: set[str]` through `BlastPointResult` → `FacadeResult` → CLI/GUI display layers. No changes to existing method signatures.

**Tech Stack:** Python 3.13, numpy, scipy, PySide6, matplotlib, pytest

---

### Task 1: Add `lookup_flagged` and `lookup_batch_flagged` to Table2D

**Files:**
- Modify: `ufc_blast/core/interpolation.py:247-372`
- Test: `tests/test_interpolation.py`

- [ ] **Step 1: Write failing tests for `lookup_flagged`**

Add to the end of `tests/test_interpolation.py`:

```python
class TestTable2DFlagged:
    """Tests for extrapolation-aware lookup methods."""

    @pytest.fixture(scope="class")
    def calpha(self):
        return Table2D.from_csv(CALPHA_CSV, family_col="ps0_kpa", x_col="angle_deg", value_col="calpha")

    # -- lookup_flagged --

    def test_lookup_flagged_in_range_returns_false(self, calpha):
        """In-range family should return (value, False)."""
        val, flag = calpha.lookup_flagged(x=0.0, family=6.89476)
        assert flag is False
        assert val > 0.0

    def test_lookup_flagged_below_range_returns_true(self, calpha):
        """Family below table min should return (value, True), not raise."""
        fam_min = calpha.family_levels[0]
        val, flag = calpha.lookup_flagged(x=0.0, family=fam_min * 0.5)
        assert flag is True
        assert val > 0.0

    def test_lookup_flagged_above_range_returns_true(self, calpha):
        """Family above table max should return (value, True), not raise."""
        fam_max = calpha.family_levels[-1]
        val, flag = calpha.lookup_flagged(x=0.0, family=fam_max * 2.0)
        assert flag is True
        assert val > 0.0

    def test_lookup_flagged_matches_lookup_for_in_range(self, calpha):
        """Flagged and unflagged must return identical values for in-range queries."""
        val_flagged, flag = calpha.lookup_flagged(x=10.0, family=6.89476)
        val_plain = calpha.lookup(x=10.0, family=6.89476)
        assert val_flagged == pytest.approx(val_plain, rel=1e-12)
        assert flag is False

    # -- lookup_batch_flagged --

    def test_lookup_batch_flagged_mixed(self, calpha):
        """Batch with mixed in-range and out-of-range families."""
        fam_min = calpha.family_levels[0]
        fam_max = calpha.family_levels[-1]
        xs = np.array([0.0, 0.0, 0.0])
        families = np.array([fam_min * 0.5, 6.89476, fam_max * 2.0])
        vals, flags = calpha.lookup_batch_flagged(xs, families)
        assert flags[0] is True or flags[0] == True  # below range
        assert flags[1] is False or flags[1] == False  # in range
        assert flags[2] is True or flags[2] == True  # above range
        assert all(v > 0.0 for v in vals)

    def test_lookup_batch_flagged_all_in_range(self, calpha):
        """All in-range families should produce all-False mask."""
        xs = np.array([0.0, 10.0, 20.0])
        families = np.array([6.89476, 68.9476, 689.476])
        vals, flags = calpha.lookup_batch_flagged(xs, families)
        assert not np.any(flags)

    def test_lookup_batch_flagged_values_match_unflagged(self, calpha):
        """Flagged batch values must match unflagged batch for same inputs."""
        xs = np.array([0.0, 10.0, 20.0])
        families = np.array([6.89476, 68.9476, 689.476])
        vals_flagged, _ = calpha.lookup_batch_flagged(xs, families)
        vals_plain = calpha.lookup_batch(xs, families)
        np.testing.assert_allclose(vals_flagged, vals_plain, rtol=1e-12)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_interpolation.py::TestTable2DFlagged -v`
Expected: FAIL with `AttributeError: 'Table2D' object has no attribute 'lookup_flagged'`

- [ ] **Step 3: Implement `lookup_flagged` in Table2D**

Add after the existing `lookup` method (after line 315) in `ufc_blast/core/interpolation.py`:

```python
    def lookup_flagged(self, x: float, family: float, method: str = "log") -> tuple[float, bool]:
        """Bivariate interpolation with out-of-range detection.

        Same as :meth:`lookup` but instead of raising ``ValueError`` when
        *family* is outside the table range, clamps to the nearest boundary
        and returns ``(value, True)``.  In-range queries return
        ``(value, False)``.

        Parameters
        ----------
        x : float
            Lookup variable (e.g. angle in degrees).
        family : float
            Family parameter (e.g. Ps0 in kPa).
        method : str
            ``'log'`` or ``'linear'`` for the family axis.

        Returns
        -------
        value : float
            Interpolated (or extrapolated) value.
        out_of_range : bool
            ``True`` if *family* was outside the table range.
        """
        fam_min, fam_max = self.family_levels[0], self.family_levels[-1]
        out_of_range = family < fam_min or family > fam_max

        if out_of_range:
            family = float(np.clip(family, fam_min, fam_max))

        # Find bracketing family levels
        idx = np.searchsorted(self.family_levels, family)

        if idx == 0:
            tbl = self.x_tables[float(self.family_levels[0])]
            x = self._clamp_x(x, tbl)
            return tbl.lookup(x, method="linear"), out_of_range

        if idx == len(self.family_levels):
            tbl = self.x_tables[float(self.family_levels[-1])]
            x = self._clamp_x(x, tbl)
            return tbl.lookup(x, method="linear"), out_of_range

        fam_lo = float(self.family_levels[idx - 1])
        fam_hi = float(self.family_levels[idx])

        tbl_lo = self.x_tables[fam_lo]
        tbl_hi = self.x_tables[fam_hi]
        x_min = max(tbl_lo.x[0], tbl_hi.x[0])
        x_max = min(tbl_lo.x[-1], tbl_hi.x[-1])
        x = self._clamp_x(x, x_min=x_min, x_max=x_max)

        val_lo = tbl_lo.lookup(x, method="linear")
        val_hi = tbl_hi.lookup(x, method="linear")

        if method == "log":
            t = (np.log(family) - np.log(fam_lo)) / (np.log(fam_hi) - np.log(fam_lo))
            value = float(np.exp((1 - t) * np.log(val_lo) + t * np.log(val_hi)))
        elif method == "linear":
            t = (family - fam_lo) / (fam_hi - fam_lo)
            value = float((1 - t) * val_lo + t * val_hi)
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'log' or 'linear'.")

        return value, out_of_range
```

- [ ] **Step 4: Implement `lookup_batch_flagged` in Table2D**

Add after `lookup_flagged` in `ufc_blast/core/interpolation.py`:

```python
    def lookup_batch_flagged(
        self, xs: np.ndarray, families: np.ndarray, method: str = "log"
    ) -> tuple[np.ndarray, np.ndarray]:
        """Vectorized bivariate interpolation with out-of-range detection.

        Same as :meth:`lookup_batch` but also returns a boolean mask
        indicating which points had *family* outside the table range.

        Returns
        -------
        values : np.ndarray
            Interpolated values, same shape as inputs.
        out_of_range : np.ndarray of bool
            ``True`` where *family* was outside the table range.
        """
        xs = np.asarray(xs, dtype=float)
        families = np.asarray(families, dtype=float)

        fam_min = self.family_levels[0]
        fam_max = self.family_levels[-1]
        out_of_range = (families < fam_min) | (families > fam_max)

        result = np.empty_like(xs)

        # Clamp families for computation (same extrapolation as before)
        families_clamped = np.clip(families, fam_min, fam_max)

        idxs = np.searchsorted(self.family_levels, families_clamped)
        idxs = np.clip(idxs, 1, len(self.family_levels) - 1)

        for lo_idx in range(len(self.family_levels) - 1):
            hi_idx = lo_idx + 1
            mask = idxs == hi_idx
            if not np.any(mask):
                continue

            fam_lo = float(self.family_levels[lo_idx])
            fam_hi = float(self.family_levels[hi_idx])
            tbl_lo = self.x_tables[fam_lo]
            tbl_hi = self.x_tables[fam_hi]

            x_min = max(tbl_lo.x[0], tbl_hi.x[0])
            x_max = min(tbl_lo.x[-1], tbl_hi.x[-1])
            x_sub = np.clip(xs[mask], x_min, x_max)

            val_lo = tbl_lo.lookup_batch(x_sub, method="linear")
            val_hi = tbl_hi.lookup_batch(x_sub, method="linear")

            fam_sub = families_clamped[mask]
            if method == "log":
                t = (np.log(fam_sub) - np.log(fam_lo)) / (np.log(fam_hi) - np.log(fam_lo))
                result[mask] = np.exp((1 - t) * np.log(val_lo) + t * np.log(val_hi))
            else:
                t = (fam_sub - fam_lo) / (fam_hi - fam_lo)
                result[mask] = (1 - t) * val_lo + t * val_hi

        return result, out_of_range
```

- [ ] **Step 5: Refactor `lookup_batch` to delegate to `lookup_batch_flagged`**

Replace the existing `lookup_batch` body (lines 317-372) with:

```python
    def lookup_batch(
        self, xs: np.ndarray, families: np.ndarray, method: str = "log"
    ) -> np.ndarray:
        """Vectorized bivariate interpolation for arrays of (x, family).

        Parameters
        ----------
        xs : np.ndarray
            X-axis query points.
        families : np.ndarray
            Family-axis query points.
        method : str
            ``'log'`` or ``'linear'`` for the family axis.

        Returns
        -------
        np.ndarray
            Interpolated values, same shape as inputs.
        """
        values, _ = self.lookup_batch_flagged(xs, families, method=method)
        return values
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `uv run pytest tests/test_interpolation.py -v`
Expected: ALL PASS (existing tests + new `TestTable2DFlagged`)

- [ ] **Step 7: Commit**

```bash
git add ufc_blast/core/interpolation.py tests/test_interpolation.py
git commit -m "feat: add lookup_flagged and lookup_batch_flagged to Table2D"
```

---

### Task 2: Add `extrapolated` field to `BlastPointResult` and wire it through `compute_point`

**Files:**
- Modify: `ufc_blast/core/blast_params.py:36-63` (dataclass), `ufc_blast/core/blast_params.py:132-188` (compute_point)
- Test: `tests/test_blast_params.py`

- [ ] **Step 1: Write failing tests**

Add to the end of `tests/test_blast_params.py`:

```python
class TestExtrapolationFlags:
    """Tests for the extrapolated field on BlastPointResult."""

    def test_normal_range_has_empty_extrapolated(self):
        """W=200kg, R=30m — well within all table ranges."""
        result = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0)
        assert result.extrapolated == set()

    def test_far_range_flags_iralpha(self):
        """R=120m, W=200kg → Ps0 ≈ 4.24 kPa, below irAlpha min (4.83 kPa)."""
        result = compute_point(R_alpha=120.0, alpha_deg=0.0, W=200.0)
        assert "iralpha" in result.extrapolated

    def test_far_range_does_not_crash(self):
        """compute_point must not raise ValueError for out-of-range Ps0."""
        result = compute_point(R_alpha=120.0, alpha_deg=0.0, W=200.0)
        assert isinstance(result, BlastPointResult)
        assert result.Ps0 > 0.0
        assert result.ir_alpha > 0.0

    def test_batch_matches_single_extrapolation_flag(self):
        """Batch and single compute must agree on extrapolation flags."""
        import numpy as np
        R_alphas = np.array([30.0, 120.0])
        alpha_degs = np.array([0.0, 0.0])
        batch = compute_points_batch(R_alphas, alpha_degs, W=200.0)

        single_near = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0)
        single_far = compute_point(R_alpha=120.0, alpha_deg=0.0, W=200.0)

        assert batch[0].extrapolated == single_near.extrapolated
        assert batch[1].extrapolated == single_far.extrapolated
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_blast_params.py::TestExtrapolationFlags -v`
Expected: FAIL — `BlastPointResult` has no `extrapolated` field, `compute_point` raises `ValueError` for R=120m.

- [ ] **Step 3: Add `extrapolated` field to `BlastPointResult`**

In `ufc_blast/core/blast_params.py`, add import at top (line 6, inside existing imports):

```python
from dataclasses import dataclass, field
```

Then modify the dataclass (after the `b: float` line, line 63):

```python
    extrapolated: set[str] = field(default_factory=set)
```

- [ ] **Step 4: Update `compute_point` to use `lookup_flagged` and collect flags**

Replace the body of `compute_point` (lines 149-188) with:

```python
    if not _tables:
        load_ufc_tables()

    W_cbrt = W ** (1.0 / 3.0)
    Z = R_alpha / W_cbrt

    extrapolated: set[str] = set()

    # Select tables based on burst type
    ps0_table = _tables["surface_ps0"] if burst_type == "surface" else _tables["ps0"]
    tA_table = _tables["surface_tA"] if burst_type == "surface" else _tables["tA"]
    t0_table = _tables["surface_t0"] if burst_type == "surface" else _tables["t0"]

    # Incident overpressure
    Ps0 = ps0_table.lookup(Z)

    # Reflection coefficient from Figure 2-193 (shared across burst types)
    C_alpha, oor_calpha = _tables["calpha"].lookup_flagged(x=alpha_deg, family=Ps0)
    if oor_calpha:
        extrapolated.add("calpha")
    Pr_alpha = C_alpha * Ps0

    # Reflected specific impulse from Figure 2-194 (kPa·ms/kg^1/3 → kPa·s)
    ir_alpha_scaled, oor_iralpha = _tables["iralpha"].lookup_flagged(x=alpha_deg, family=Ps0)
    if oor_iralpha:
        extrapolated.add("iralpha")
    ir_alpha = ir_alpha_scaled * W_cbrt / 1000.0

    # Arrival time (ms/kg^1/3 → s)
    tA = tA_table.lookup(Z) * W_cbrt / 1000.0

    # Positive phase duration (ms/kg^1/3 → s)
    t0 = t0_table.lookup(Z) * W_cbrt / 1000.0

    # Friedlander decay parameter
    b = solve_friedlander_b(Pr_alpha, t0, ir_alpha)

    return BlastPointResult(
        Ps0=Ps0,
        C_alpha=C_alpha,
        Pr_alpha=Pr_alpha,
        ir_alpha=ir_alpha,
        tA=tA,
        t0=t0,
        b=b,
        extrapolated=extrapolated,
    )
```

- [ ] **Step 5: Update `compute_points_batch` to use `lookup_batch_flagged` and collect flags**

Replace the body of `compute_points_batch` (lines 216-257) with:

```python
    if not _tables:
        load_ufc_tables()

    R_alphas = np.asarray(R_alphas, dtype=float)
    alpha_degs = np.asarray(alpha_degs, dtype=float)
    n = len(R_alphas)

    W_cbrt = W ** (1.0 / 3.0)
    Z = R_alphas / W_cbrt

    # Select tables based on burst type
    ps0_table = _tables["surface_ps0"] if burst_type == "surface" else _tables["ps0"]
    tA_table = _tables["surface_tA"] if burst_type == "surface" else _tables["tA"]
    t0_table = _tables["surface_t0"] if burst_type == "surface" else _tables["t0"]

    # Bulk 1D lookups
    Ps0 = ps0_table.lookup_batch(Z)
    tA_vals = tA_table.lookup_batch(Z) * W_cbrt / 1000.0
    t0_vals = t0_table.lookup_batch(Z) * W_cbrt / 1000.0

    # Bulk 2D lookups with flags
    C_alpha, oor_calpha = _tables["calpha"].lookup_batch_flagged(xs=alpha_degs, families=Ps0)
    Pr_alpha = C_alpha * Ps0

    ir_alpha_scaled, oor_iralpha = _tables["iralpha"].lookup_batch_flagged(xs=alpha_degs, families=Ps0)
    ir_alpha = ir_alpha_scaled * W_cbrt / 1000.0

    # Friedlander b — vectorized Newton's method
    b_vals = _solve_friedlander_b_batch(Pr_alpha, t0_vals, ir_alpha)

    return [
        BlastPointResult(
            Ps0=float(Ps0[i]),
            C_alpha=float(C_alpha[i]),
            Pr_alpha=float(Pr_alpha[i]),
            ir_alpha=float(ir_alpha[i]),
            tA=float(tA_vals[i]),
            t0=float(t0_vals[i]),
            b=float(b_vals[i]),
            extrapolated=(
                ({"calpha"} if oor_calpha[i] else set())
                | ({"iralpha"} if oor_iralpha[i] else set())
            ),
        )
        for i in range(n)
    ]
```

- [ ] **Step 6: Fix existing test that checks `__dict__` fields**

In `tests/test_blast_params.py`, method `test_result_fields_are_finite` (line 214) iterates `result.__dict__.items()` and asserts `math.isfinite(value)`. The new `extrapolated` field is a `set`, not a float. Replace:

```python
    def test_result_fields_are_finite(self):
        """All numeric result fields must be finite floats."""
        result = compute_point(R_alpha=45.0, alpha_deg=20.0, W=200.0)
        for field_name, value in result.__dict__.items():
            if field_name == "extrapolated":
                continue
            assert math.isfinite(value), f"Field {field_name} is not finite: {value}"
```

- [ ] **Step 7: Run all tests**

Run: `uv run pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 8: Commit**

```bash
git add ufc_blast/core/blast_params.py tests/test_blast_params.py
git commit -m "feat: add extrapolated flag to BlastPointResult, wire through compute_point and compute_points_batch"
```

---

### Task 3: Add extrapolation warnings to CLI

**Files:**
- Modify: `ufc_blast/cli.py:60-96` (formatters), `ufc_blast/cli.py:150-203` (subcommands)
- Test: manual verification via command line

- [ ] **Step 1: Update `_print_point_result` to show extrapolation warnings**

In `ufc_blast/cli.py`, replace `_print_point_result` (lines 61-79) with:

```python
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
    if result.extrapolated:
        tables = ", ".join(sorted(result.extrapolated))
        print(f"  \u26a0 EXTRAPOLATED: {tables} (Ps0 outside table range)")
        print(sep)
```

- [ ] **Step 2: Update `_result_to_dict` to include `extrapolated` key**

In `ufc_blast/cli.py`, replace `_result_to_dict` (lines 82-96) with:

```python
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
        "extrapolated": ",".join(sorted(result.extrapolated)),
    }
```

- [ ] **Step 3: Add `extrapolated` to `_TABLE_COLUMNS`**

In `ufc_blast/cli.py`, replace the `_TABLE_COLUMNS` list (lines 99-111) with:

```python
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
    ("extrapolated",  "extrap.",    12),
]
```

- [ ] **Step 4: Add extrapolation summary to `_cmd_compute`**

In `ufc_blast/cli.py`, add the summary after the output block. Replace lines 190-203 of `_cmd_compute` with:

```python
    if not rows:
        print("ERROR: no grid points could be computed.", file=sys.stderr)
        return 1

    fmt = getattr(args, "format", "table")

    if args.output:
        _write_csv(rows, args.output)
    elif fmt == "json":
        _print_json(rows)
    else:
        _print_table(rows)

    # Extrapolation summary
    extrap_rows = [r for r in rows if r["extrapolated"]]
    if extrap_rows:
        all_tables: set[str] = set()
        for r in extrap_rows:
            all_tables.update(r["extrapolated"].split(","))
        print(f"  \u26a0 {len(extrap_rows)} of {len(rows)} points extrapolated "
              f"(outside UFC table range)")
        print(f"    Affected tables: {', '.join(sorted(all_tables))}")

    return 0
```

- [ ] **Step 5: Verify CLI point command no longer crashes**

Run: `uv run ufc-blast point --W 200 --R 120 --Hc 5 --alpha 0`
Expected: Prints result with `⚠ EXTRAPOLATED: iralpha (Ps0 outside table range)` instead of crashing.

- [ ] **Step 6: Verify CLI compute command shows summary**

Run: `uv run ufc-blast compute --W 200 --R 110 --Hc 5 --width 10 --height 8 --step 2`
Expected: Table prints, followed by `⚠ N of M points extrapolated (outside UFC table range)`.

- [ ] **Step 7: Run all tests**

Run: `uv run pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 8: Commit**

```bash
git add ufc_blast/cli.py
git commit -m "feat: add extrapolation warnings to CLI point and compute output"
```

---

### Task 4: Add extrapolation indicators to GUI

**Files:**
- Modify: `ufc_blast/gui.py:96-125` (meshgrid builder), `ufc_blast/gui.py:168-224` (contour drawing), `ufc_blast/gui.py:426-439` (label update), `ufc_blast/gui.py:635-651` (click handler)

- [ ] **Step 1: Add `_build_extrap_mask` helper**

In `ufc_blast/gui.py`, add after the `_build_meshgrid` function (after line 125):

```python
def _build_extrap_mask(
    grid_points: list[GridPoint],
    result_map: dict[tuple[float, float], BlastPointResult],
) -> np.ndarray:
    """Build a 2D boolean array: True where a point is extrapolated.

    Returns
    -------
    mask : np.ndarray, shape (ny, nx), dtype bool
    """
    dx_vals = sorted({gp.dx for gp in grid_points})
    dy_vals = sorted({gp.dy for gp in grid_points})
    nx, ny = len(dx_vals), len(dy_vals)

    mask = np.array(
        [bool(result_map[(dx, dy)].extrapolated) for dy in dy_vals for dx in dx_vals],
        dtype=bool,
    ).reshape(ny, nx)

    return mask
```

- [ ] **Step 2: Update `draw_contours` to accept and render extrapolation hatching**

In `ufc_blast/gui.py`, modify the `draw_contours` method signature (line 168) to add `result_map` parameter:

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

Then, after the Mach stem boundary line block (after the `ax.plot(mach_xs, ...)` block, around line 208), add the hatched overlay:

```python
        # Extrapolation hatched overlay
        extrap_mask = _build_extrap_mask(grid_points, result_map)
        if np.any(extrap_mask):
            # Draw hatched contourf over extrapolated regions
            ax.contourf(
                X, Y, extrap_mask.astype(float),
                levels=[0.5, 1.5],
                colors="none",
                hatches=["///"],
                zorder=4,
            )
            ax.contour(
                X, Y, extrap_mask.astype(float),
                levels=[0.5],
                colors="gray",
                linewidths=1.0,
                linestyles="--",
                zorder=4,
            )
```

- [ ] **Step 3: Update the call site in `_on_compute_done`**

The `draw_contours` call in `_on_compute_done` (line 623) already passes `self._result_map` stored on line 608. The existing call is:

```python
        self._contour_canvas.draw_contours(
            self._grid_points, self._result_map, self._width, self._height,
            mach_curve=self._mach_curve,
        )
```

This already passes `result_map` as the second positional argument, which now matches the updated signature. **No change needed** — the parameter was already being passed, it just wasn't declared in the old signature. Wait — looking at the old signature at line 168:

```python
    def draw_contours(self, grid_points, result_map, width, height, mach_curve=None):
```

The old call at line 623 already passes `self._result_map`. So the call site is already correct. The only change is that `draw_contours` now uses `result_map` internally for the extrapolation overlay.

- [ ] **Step 4: Update `_update_labels` in FriedlanderPanel to show extrapolation**

In `ufc_blast/gui.py`, replace `_update_labels` (lines 426-439) with:

```python
    def _update_labels(self, gp: GridPoint, res: BlastPointResult, mach_zone: bool = False) -> None:
        coord_text = f"Point:  dx = {gp.dx:.2f} m,  dy = {gp.dy:.2f} m"
        if mach_zone:
            coord_text += "  (Mach zone)"
        if res.extrapolated:
            tables = ", ".join(sorted(res.extrapolated))
            coord_text += f"  (\u26a0 extrapolated: {tables})"
        self._param_labels["coord"].setText(coord_text)
        self._param_labels["R_alpha"].setText(f"R_α = {gp.R_alpha:.3f} m")
        self._param_labels["alpha"].setText(f"α = {gp.alpha_deg:.2f}°")
        self._param_labels["Ps0"].setText(f"Ps0 = {res.Ps0:.2f} kPa")
        self._param_labels["C_alpha"].setText(f"C_α = {res.C_alpha:.4f}")
        self._param_labels["Pr"].setText(f"Pr_α = {res.Pr_alpha:.2f} kPa")
        self._param_labels["ir"].setText(f"ir_α = {res.ir_alpha*1e3:.2f} kPa·ms")
        self._param_labels["tA"].setText(f"tA = {res.tA*1e3:.3f} ms")
        self._param_labels["t0"].setText(f"t0 = {res.t0*1e3:.3f} ms")
        self._param_labels["b"].setText(f"b = {res.b:.4f}")
```

- [ ] **Step 5: Update status bar message in `_on_point_selected`**

In `ufc_blast/gui.py`, replace `_on_point_selected` (lines 635-651) with:

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
        status = (
            f"Selected point  dx={gp.dx:.1f} m  dy={gp.dy:.1f} m  "
            f"Pr_α={res.Pr_alpha:.1f} kPa"
        )
        if mach_zone:
            status += "  [Mach zone]"
        if res.extrapolated:
            status += f"  [\u26a0 extrapolated: {', '.join(sorted(res.extrapolated))}]"
        self._status.showMessage(status)
```

- [ ] **Step 6: Run all tests**

Run: `uv run pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 7: Manual verification — launch GUI at far range**

Run: `uv run ufc-blast gui --W 200 --R 110 --Hc 5 --width 10 --height 8 --step 1`
Expected: Contour plot shows hatched overlay on points outside UFC range. Clicking a hatched point shows `⚠ extrapolated: iralpha` in the panel.

- [ ] **Step 8: Commit**

```bash
git add ufc_blast/gui.py
git commit -m "feat: add extrapolation hatching and labels to GUI"
```

---

### Task 5: Add facade-level extrapolation tests

**Files:**
- Modify: `tests/test_compute_facade.py`

- [ ] **Step 1: Write extrapolation tests**

Add to the end of `tests/test_compute_facade.py`:

```python
class TestExtrapolationInFacade:
    """Tests for extrapolation flag propagation through compute_facade."""

    def test_far_range_facade_has_extrapolated_points(self):
        """At R=110m, W=200kg, all points should have Ps0 below irAlpha min."""
        facade = compute_facade(R=110.0, W=200.0, Hc=5.0, width=4.0, height=4.0, step=2.0)
        extrap_count = sum(
            1 for res in facade.result_map.values() if res.extrapolated
        )
        assert extrap_count > 0, "Expected some extrapolated points at far range"

    def test_close_range_facade_no_extrapolation(self):
        """At R=30m, W=200kg, no points should be extrapolated."""
        facade = compute_facade(R=30.0, W=200.0, Hc=5.0, width=10.0, height=8.0, step=2.0)
        extrap_count = sum(
            1 for res in facade.result_map.values() if res.extrapolated
        )
        assert extrap_count == 0, "No extrapolation expected at R=30m"

    def test_extrapolated_set_contains_table_names(self):
        """Extrapolated set should contain valid table key strings."""
        facade = compute_facade(R=110.0, W=200.0, Hc=5.0, width=4.0, height=4.0, step=2.0)
        valid_keys = {"calpha", "iralpha", "ps0", "tA", "t0",
                      "surface_ps0", "surface_tA", "surface_t0", "triple_point"}
        for res in facade.result_map.values():
            assert res.extrapolated <= valid_keys, (
                f"Unexpected table keys: {res.extrapolated - valid_keys}"
            )

    def test_mach_stem_propagates_extrapolation(self):
        """Mach-overwritten points should inherit donor's extrapolated set."""
        facade = compute_facade(R=110.0, W=200.0, Hc=5.0, width=4.0, height=4.0, step=2.0)
        if not facade.mach_curve:
            pytest.skip("No Mach stem at this range")

        from collections import defaultdict
        columns = defaultdict(list)
        for gp in facade.grid_points:
            columns[gp.dx].append(gp.dy)

        mach_dict = {dx: mach_dy for dx, mach_dy in facade.mach_curve}
        for dx, mach_dy in mach_dict.items():
            # Find donor (nearest above Mach line)
            above = [dy for dy in columns[dx] if dy >= mach_dy]
            if not above:
                continue
            donor_dy = min(above)
            donor_extrap = facade.result_map[(dx, donor_dy)].extrapolated
            # All below-Mach points must match donor
            for dy in columns[dx]:
                if dy < mach_dy:
                    assert facade.result_map[(dx, dy)].extrapolated == donor_extrap
```

- [ ] **Step 2: Run the new tests**

Run: `uv run pytest tests/test_compute_facade.py::TestExtrapolationInFacade -v`
Expected: ALL PASS

- [ ] **Step 3: Run full test suite**

Run: `uv run pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 4: Commit**

```bash
git add tests/test_compute_facade.py
git commit -m "test: add extrapolation flag propagation tests for compute_facade"
```
