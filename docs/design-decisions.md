# Design Decisions

## Architecture: Layered Library + Thin Interfaces (Approach B)

Three approaches were considered during brainstorming:
- **A — Monolithic single-file script**: Rejected. No reuse, hard to test.
- **B — Layered library + thin CLI/GUI**: **Chosen.** Clean separation: `core/` for calculation, `cli.py` and `gui.py` as thin wrappers.
- **C — Config-driven plugin system**: Rejected. Over-engineered for a calculator.

## GUI: PySide6 + matplotlib (not PyVista)

Originally designed with PyVista 3D visualization. Changed because:
- 3D rotation not needed — facade is a flat plane, 2D contour is the natural view
- PyVista/VTK was heavy (~500MB), laggy on large grids
- matplotlib contourf provides publication-quality pressure maps
- PySide6 gives native Qt widgets (spinboxes, layout panels)

**Current layout**: Top input bar + left contour view (65%) + right Friedlander panel (35%). Click on contour to see P(t) waveform for that point.

## TUI: Rejected

Explicitly dropped as redundant — CLI handles automation, GUI handles visualization. A TUI would be neither.

## Interpolation: Log-Log Default

Blast parameters span orders of magnitude (Ps0 from ~1 kPa to ~50,000 kPa). Log-log interpolation is exact for power-law relationships, which is the dominant physics. Linear interpolation is only used for cross-validation against the Excel spreadsheet.

2D tables (Cα, irα): log on Ps0/family axis, linear on angle axis — angle has no power-law behavior.

## Vectorized Newton for b-solver

Single-point b-solver uses `scipy.optimize.brentq` (robust, bracketed). Batch computation uses vectorized Newton's method (converges in ~10 iterations for all points simultaneously). This achieves 2000+ points in ~1ms vs ~2s with a Python loop.

## Mach Stem as Post-Processing

Rather than special-casing the Mach zone in the core loop, `apply_mach_stem()` runs after `compute_points_batch()` and overwrites results below the triple point height with uniform values from above. This keeps the core calculation simple and testable.

## Mach Visualization: Unified Contourf

Several approaches were tried:
- Separate `contourf` (above) + `pcolormesh` (below) — looked disconnected
- Masked NaN for contour lines below Mach — made it worse
- **Final**: Single `contourf` across entire grid (Mach values already overwritten), with a thin white dashed line overlay for the boundary

## No Pandas

All CSV loading uses `np.loadtxt` or stdlib `csv`. Pandas is not a dependency — unnecessary for simple columnar data that loads directly into numpy arrays.

## Package Manager: uv Only

No pip for development. `uv run` for all commands. `uv pip install --system` only in CI (GitHub Actions with `astral-sh/setup-uv@v5`).
