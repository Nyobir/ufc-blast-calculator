# Extrapolation Warnings Design

## Problem

`Table2D.lookup_batch()` silently extrapolates when a family parameter (Ps0) falls outside the UFC table range. `Table2D.lookup()` raises `ValueError` for the same case. This means:

- CLI `compute` and GUI silently return unvalidated data for out-of-range points
- CLI `point` crashes instead of showing the extrapolated value
- No consumer can distinguish UFC-validated results from extrapolated ones

Real-world scenarios hit this:
- **irAlpha underflow**: R > ~110m for W=200kg (Ps0 < 4.83 kPa)
- **Cα overflow**: R < 0.51m for W=200kg (Ps0 > 34,474 kPa)

## Solution

Compute extrapolated values AND flag them per-point, per-table. Every `BlastPointResult` carries a `set[str]` of table names that required extrapolation. Empty set = fully UFC-validated.

## Design

### 1. Interpolation Layer (`interpolation.py`)

Add flagged variants to `Table2D`:

**`lookup_flagged(x, family, method) -> tuple[float, bool]`**
- Same logic as `lookup()` but returns `(value, True)` when `family < fam_min` or `family > fam_max` instead of raising `ValueError`.
- Clamps family to nearest level and computes as if at boundary (same behavior as current `lookup_batch` extrapolation).

**`lookup_batch_flagged(xs, families, method) -> tuple[ndarray, ndarray[bool]]`**
- Returns `(values, out_of_range_mask)`.
- `out_of_range_mask[i]` is `True` when `families[i] < fam_min` or `families[i] > fam_max`.
- Extrapolation uses existing `np.clip` logic — just also records which indices were clipped.

Existing `lookup()` and `lookup_batch()` remain unchanged for backward compatibility. `lookup_batch()` internally delegates to `lookup_batch_flagged()` and discards the flag.

### 2. BlastPointResult (`blast_params.py`)

```python
@dataclass
class BlastPointResult:
    Ps0: float
    C_alpha: float
    Pr_alpha: float
    ir_alpha: float
    tA: float
    t0: float
    b: float
    extrapolated: set[str] = field(default_factory=set)
```

Table keys in the set match `_tables` dict keys: `"calpha"`, `"iralpha"`, `"ps0"`, `"tA"`, `"t0"`, `"surface_ps0"`, `"surface_tA"`, `"surface_t0"`, `"triple_point"`.

**`compute_point()`**: Uses `lookup_flagged()` for each 2D table (calpha, iralpha). Collects flags into the set. No longer crashes for out-of-range Ps0.

**`compute_points_batch()`**: Uses `lookup_batch_flagged()` for each 2D table. Builds per-point extrapolated sets from the boolean masks.

**`apply_mach_stem()`**: When overwriting a point below the Mach line with a donor point's result, the donor's `extrapolated` set propagates naturally (same BlastPointResult object is reused).

### 3. CLI Output (`cli.py`)

**`point` subcommand**: No longer crashes. Prints the result and appends warning lines:

```
  ⚠ EXTRAPOLATED: iralpha (Ps0 outside table range)
```

**`compute` subcommand**:

Table/stdout format — summary after the table:
```
  ⚠ 12 of 99 points extrapolated (outside UFC table range)
    Affected tables: calpha, iralpha
```

CSV format — add `extrapolated` column:
- Empty string for clean points
- Comma-joined table names for flagged points (e.g. `"calpha,iralpha"`)

JSON format — add `"extrapolated": []` or `"extrapolated": ["calpha", "iralpha"]` per point.

### 4. GUI Display (`gui.py`)

**Contour overlay**: After drawing the main `contourf`, overlay a semi-transparent hatched region for points where `extrapolated` is non-empty. Use `contourf` with a binary mask and `hatches='///'` to mark unreliable zones while keeping pressure colors visible underneath.

**Friedlander panel**: When a clicked point has non-empty `extrapolated`, the coord label shows:
```
Point:  dx = 4.00 m,  dy = -3.00 m  (⚠ extrapolated: iralpha)
```

**Status bar**: Includes `[extrapolated]` suffix when applicable.

### 5. Testing

**test_interpolation.py**:
- `lookup_batch_flagged` returns `(values, mask)` with correct mask for mixed in-range/out-of-range families
- `lookup_flagged` returns `(value, True)` for out-of-range family, `(value, False)` for in-range
- Flagged values match unflagged values (same extrapolation math)

**test_blast_params.py**:
- `compute_point` at R=120m, W=200kg populates `extrapolated` with `{"iralpha"}`
- `compute_points_batch` for same scenario matches
- Normal-range scenario (R=30m, W=200kg) has empty `extrapolated`

**test_compute_facade.py**:
- Facade at R=110m has points with non-empty `extrapolated` in `result_map`
- Count of extrapolated points matches expected grid coverage

## Non-Goals

- Negative phase Friedlander (not modeled)
- New interpolation methods for out-of-range data
- GUI modal dialogs or blocking warnings
- Changing existing `lookup()` / `lookup_batch()` return signatures
