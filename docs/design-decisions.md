# Design Decisions

## Layered Library + Thin Interfaces

The calculation lives in `ufc_blast/core/`; CLI and GUI both call
`compute_facade()` for spatial loads. A monolithic script would make the
numerical path harder to test; a plugin system would add unnecessary
complexity. `compute_point()` remains a separate regular-reflection query for
one specified angle and does not apply the facade's Mach-region correction.

## 2D GUI Instead of a 3D Dependency

A flat facade is naturally represented by a pressure contour and a waveform
panel. PySide6 provides native widgets and Matplotlib produces the contour.
PyVista/VTK was not needed for this view; no benchmark or package-size
comparison with it is claimed. A TUI was not added because CLI covers
scripting and the GUI covers interactive inspection.

## Chart Interpolation and Flags

Positive 1D chart values use log–log interpolation. For 2D reflection tables,
interpolation is linear in angle within each pressure family and logarithmic
between positive pressure families. Requests outside chart ranges are clamped
and flagged, not extrapolated without notice; the chart-selection switch at
`Hc/W^(1/3) = 0.397` is a convention rather than a physical discontinuity.
`compute_facade()` warns for scaled heights from 0.3 through 0.5 and invites
comparison of the air and surface chart branches.

## Impulse-Consistent Friedlander Parameter

The scalar `b` solver uses `scipy.optimize.brentq` on `[0.01, 50]`. Batch
calculation uses a vectorized, bounded Newton iteration rather than a Python
loop over points. If no root in the allowed bracket reproduces the tabulated
impulse, the result is clipped and flagged `friedlander_b`; it must not be
presented as an unflagged impulse-consistent curve. Runtime depends on the
machine, dependencies and workload; no fixed grid speedup is claimed here.
See [friedlander-parameters.md](friedlander-parameters.md).

## Ground-Evaluated Mach Front

For the revised procedure introduced with version 1.3.0, `apply_mach_stem()`
runs after the regular batch calculation on air-burst facades. Figure 2-13
locates the triple-point boundary, but is **not** where the Mach-front blast
parameters are sampled. Below that boundary, `compute_mach_front()` follows
UFC 3-340-02 §2-13.2 and Problems 2A-2/2A-10:

1. For each facade column, calculate ground distance
   `Rg = sqrt(R^2 + dx^2)` and initial-wave incidence angle
   `alpha_g = atan2(Rg, Hc)`.
2. Read ground reflected pressure and scaled reflected impulse from the
   official DPlot-derived Figures 2-9 and 2-10. These become incident
   pressure and impulse of the Mach wave.
3. Invert Figure 2-7 for the distances corresponding to that pressure and
   impulse, then read arrival time and positive-phase duration. An inverse
   chart-range clamp is flagged.
4. Treat the front as travelling horizontally, with wall angle
   `theta = atan2(abs(dx), R)`. Use Figures 2-193 and 2-194 for wall-reflected
   pressure and impulse. The column values are uniform with height within
   this idealized Mach region.
5. If `ir >= Pr*t0/2`, the chart-derived duration cannot support a
   positive-decay Friedlander curve with that peak and impulse. UFC defines
   the fictitious triangular duration `t_rf = 2*ir/Pr` (Eq. 2-11); the
   implementation sets `t0 = 1.005*t_rf` to obtain a finite positive `b`,
   and adds the informational flag `t0_fictitious`. This is **not** an exact
   UFC triangle or the original chart duration.

The earlier triple-point-sampled filling rule was superseded. Chart limits
still matter: four of six cases in the broader sensitivity study had
current-procedure chart-limit flags and must not be called in-range
predictions. The GUI uses one contour plot with a dashed Mach boundary.

## Packaging, Dependencies and Reproduction

PySide6 is an optional `[gui]` extra, not a core dependency. Linux aarch64
selects the 6.7.x line for Debian 12 compatibility; other platforms may
resolve newer wheels. Use `uv sync --locked --extra dev` for reproducible
core/CLI development and add `--extra gui` for desktop tests. Plain `pip`
remains supported for installed users. Windows CI uses `uv pip install
--system` and a PyInstaller onedir bundle. Pandas is not required: the
numeric tables load with NumPy/stdlib CSV handling.

Source checkouts read normalized data from repository `data/`; installed
wheels include 19 normalized CSVs under `ufc_blast/data/`; frozen binaries
use the bundled data directory. Verify sdist and wheel builds, CSV hashes,
and fresh installed-wheel CLI/GUI runs outside the source checkout. On macOS
`check_release_wheel.py` has a known temporary-path assertion issue when
`TMPDIR` passes through `/var` to `/private/var`; an unchanged rerun with a
resolved `TMPDIR` and independent wheel-data and installed-package checks
separate that checker issue from a package-data failure.
