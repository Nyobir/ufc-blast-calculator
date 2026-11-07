# UFC Blast Wave Calculator

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)

An open-source Python tool that turns the chart-based blast-loading method of
**UFC 3-340-02** into a tested, spatially resolved calculator of facade blast
loads.

Given a TNT-equivalent charge mass, a standoff distance and a height of burst,
the calculator returns, at one point or over a whole facade grid:

| Quantity | Symbol | Unit |
|---|---|---|
| Peak incident overpressure | P<sub>s0</sub> | kPa |
| Reflection coefficient | C<sub>α</sub> | – |
| Peak reflected overpressure | P<sub>rα</sub> | kPa |
| Reflected specific impulse | i<sub>rα</sub> | kPa·ms |
| Arrival time | t<sub>A</sub> | ms |
| Positive-phase duration | t<sub>0</sub> | ms |
| Friedlander decay parameter (impulse-consistent) | b | – |

Queries that fall outside the digitized chart range are **flagged, not silently
extrapolated**.

The tool was built as a transparent reference generator for validating
numerical blast simulations (e.g. DEM or CFD) against the UFC design method,
and for teaching and preliminary design studies.

---

## Installation

Requires Python ≥ 3.10.

```bash
git clone https://github.com/Nyobir/ufc-blast-calculator.git
cd ufc-blast-calculator

# recommended: uv
uv sync --extra dev                # core, CLI and non-GUI tests
uv sync --extra dev --extra gui    # additionally install the GUI where PySide6 wheels are available

# or plain pip
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"            # add ",gui" to the extras to run the desktop GUI
```

The v1.2.2 release includes a pre-built Windows bundle. Check the
[GitHub releases](https://github.com/Nyobir/ufc-blast-calculator/releases) for
platform-specific assets of the version you install; a newly tagged version's
Windows executable is available only after its build workflow succeeds and
attaches the asset. On Linux aarch64 the
`[gui]` extra selects PySide6 6.7.x, whose wheel supports Debian 12's glibc;
on other platforms it uses a current compatible PySide6 release. The GUI
requires a working Qt display or `QT_QPA_PLATFORM=offscreen` for smoke tests.

## Quick start

### Command line

```bash
# One point: 200 kg TNT, 30 m standoff, burst 5 m above ground, normal incidence
ufc-blast point --W 200 --R 30 --Hc 5

# Whole facade: 20 m x 10 m, 0.25 m grid, written to CSV
ufc-blast compute --W 200 --R 30 --Hc 5 --width 20 --height 10 --step 0.25 -o facade.csv

# Several angles of incidence at one standoff
ufc-blast batch --W 200 --R 30 --Hc 5 --alphas 0,15,30,45,60

# Interactive GUI (install the [gui] extra first): pressure map + Friedlander curve
ufc-blast gui --W 200 --R 30 --Hc 5 --width 20 --height 10 --step 0.25
```

The `point` command evaluates regular reflection at a single specified angle;
it does **not** apply the facade's Mach-region correction. In the example above
it gives `Pr_α = 64.864 kPa`. The 20 × 10 m facade's centre `(dx=0, dy=0)` is
in the Mach region and gives about `110.081 kPa` (reflected impulse
`0.622702 kPa·s`). These are distinct calculation paths, not conflicting
reference values.

Expected output of the first command (reference case used throughout the tests):

```
Burst type: Air burst (Hc/W^(1/3) = 0.8550)
Peak incident Ps0      : 29.670 kPa
Reflection coeff. Cα   : 2.1862
Peak reflected Pr_α    : 64.864 kPa
Reflected impulse ir_α : 468.400 kPa·ms
Arrival time tA        : 54.598 ms
Positive duration t0   : 19.818 ms
Friedlander b          : 1.0341
```

`compute` writes one row per grid point with the columns
`dx_m, dy_m, R_alpha_m, alpha_deg, Ps0_kPa, C_alpha, Pr_alpha_kPa,
ir_alpha_kPa_ms, tA_ms, t0_ms, b, extrapolated`. Use `--format json` to print
JSON instead of a table.

### Python API

```python
from ufc_blast.core.blast_params import compute_facade

res = compute_facade(R=30, W=200, Hc=5, width=20, height=10, step=0.25)
print(res.burst_type, len(res.grid_points))      # air 3321

p = res.result_map[(0.0, 0.0)]                   # point at burst height, facade centre
print(p.Pr_alpha, p.ir_alpha, p.b, p.extrapolated)
```

`compute_facade()` is the single computation path shared by the CLI and the
GUI, so all three interfaces give identical numbers.

## How it works

```
input (W, R, Hc, facade size, step)
  → burst type        Hc / W^(1/3) > 0.397 → air burst, otherwise surface burst
  → facade grid       slant distance R_α and incidence angle α for every point
  → batch lookup      all points at once (NumPy arrays, no per-point loop)
                      Ps0, tA, t0 from Fig. 2-7 / 2-15; Cα, irα from Fig. 2-193 / 2-194
  → Friedlander b     vectorized Newton solve against the chart impulse
                      (Brent's method for single points)
  → Mach stem         air bursts only: Fig. 2-13 locates the boundary; below it,
                      ground-evaluated Mach-front values follow §2-13.2 using
                      Figs. 2-9/2-10, inverse Fig. 2-7 and wall reflection
                      from Figs. 2-193/2-194
  → results           pressure map, time histories, CSV / JSON
```

Interpolation is log–log by default (power-law decay of the charts); 2-D tables
are logarithmic in P<sub>s0</sub> and linear in angle.

A step-by-step description of the algorithm is in
[`docs/calculation-pipeline.md`](docs/calculation-pipeline.md) (Ukrainian);
the Friedlander solver is described in
[`docs/friedlander-parameters.md`](docs/friedlander-parameters.md) and the unit
conversions in [`docs/unit-conversions.md`](docs/unit-conversions.md).

## Data

`data/` holds 19 SI-normalized CSV tables taken from the UFC 3-340-02 design
charts:

| Chart | Content | Files |
|---|---|---|
| Fig. 2-7 | free-air burst: P<sub>s0</sub>, P<sub>r</sub>, i<sub>s</sub>, i<sub>r</sub>, t<sub>A</sub>, t<sub>0</sub> vs. scaled distance | `data/free_air/figure_2_7_*.csv` |
| Figs. 2-9, 2-10 | Mach-front pressure and impulse at the ground (official DPlot data) | `data/free_air/figure_2_9_pra.csv`, `data/free_air/figure_2_10_ira.csv` |
| Fig. 2-13 | triple-point height (Mach stem) | `data/free_air/figure_2_13_triple_point.csv`, `data/triple_point/` |
| Fig. 2-15 | surface burst: same parameters as Fig. 2-7 | `data/surface/figure_2_15_*.csv` |
| Fig. 2-193 | reflection coefficient C<sub>α</sub>(α, P<sub>s0</sub>) | `data/free_air/figure_2_193_calpha.csv` |
| Fig. 2-194 | reflected scaled impulse i<sub>rα</sub>(α, P<sub>s0</sub>) | `data/free_air/figure_2_194_iralpha.csv` |

Figure 2-7 and 2-15 curves were digitized with WebPlotDigitizer; the
angle-dependent reflection tables were transcribed from a chart-reading
workbook, with 14 Figure 2-193 tail points restored from the official
electronic DPlot curve. Ground-point Figures 2-9 and 2-10 use official UFC
DPlot numerical series. The source `.grf` files are not redistributed here.
Processing scripts in `scripts/` convert source chart units to SI. Provenance
and limits are in [`docs/data-sources.md`](docs/data-sources.md).

UFC 3-340-02 is approved for public release with unlimited distribution. The
repository stores only the extracted numerical values, not the chart artwork.

## Verification

```bash
uv run --locked --extra dev pytest tests/ -q              # GUI tests skip without PySide6
uv run --locked --extra dev --extra gui pytest tests/ -q  # full suite on a supported GUI platform
QT_QPA_PLATFORM=offscreen uv run --locked --extra dev --extra gui ufc-blast gui-check
uv build --no-sources                                      # sdist and wheel
uv run --locked --extra dev python scripts/check_release_wheel.py dist/ufc_blast-1.3.0-py3-none-any.whl
```

The reviewed candidate `1d69fe4` passed **122 tests** and a GUI smoke test on
Pi and Mac; a fresh installed Mac wheel passed CLI/GUI, air/surface and 19-table
hash checks. On macOS, the wheel checker may fail only its temporary import-path
assertion when `TMPDIR` points through `/var/folders` to `/private/var/folders`;
rerun it unchanged with `TMPDIR` set to the resolved path, and independently
verify package data and installed-wheel import. This is a checker portability
issue, not evidence of a wheel payload mismatch. The Windows executable build,
CLI and GUI smoke tests also passed on the candidate branch. None of these
statements claims that `v1.3.0` has already been published.

The test suite covers interpolation, geometry, burst-type selection, the Mach
stem, the Friedlander solvers (Brent and vectorized Newton agree), facade
computation and a headless check that the GUI and CLI share one computation
path.

The scripts that reproduce the figures and tables of the accompanying paper:

| Script | Reproduces |
|---|---|
| `scripts/benchmark_table.py` | benchmark table (W = 200 kg, R = 30 m, H<sub>c</sub> = 5 m) |
| `scripts/cross_chart_check.py` | consistency check between independent charts |
| `scripts/plot_verification.py` | redraw of UFC Fig. 2-7 in imperial units for visual comparison |

The older manuscript snapshot is tagged
[`manuscript-v1.2`](https://github.com/Nyobir/ufc-blast-calculator/tree/manuscript-v1.2).
The revised SoftwareX numerical and Mac verification archive was generated at
commit `9d287e0`. The later `1d69fe4` packaging/dependency candidate had its
own Mac installed-wheel/GUI checks and Windows CI checks; its validation did
**not** repeat the earlier timing benchmark. Do not attribute the `9d287e0`
Mac timing to `1d69fe4`, or imply a like-for-like comparison with `v1.2.2`.

## Scope and limitations

- Covers the **positive phase** of the blast wave from spherical free-air and
  hemispherical surface bursts of TNT-equivalent charges, as tabulated in the
  UFC charts.
- A facade is treated as a flat, rigid, infinite reflecting plane. Clearing,
  confinement, shielding by other buildings and negative-phase effects are not
  modelled.
- Results are only as accurate as the underlying design charts and their
  digitization. Values outside the chart range are clamped and flagged in the
  `extrapolated` field.
- The Mach front is approximated as uniform over each facade column below the
  triple-point boundary. Some broader sensitivity cases carry chart-limit
  flags; they must not be presented as in-range UFC chart predictions.
- The tool is intended for research, teaching and verification of numerical
  models. It does not replace the full UFC 3-340-02 procedure or a qualified
  engineer's assessment.

## Repository layout

```
ufc_blast/
  core/interpolation.py   Table1D / Table2D, log-log interpolation, range flags
  core/geometry.py        burst type, facade grid, slant distance and angle
  core/blast_params.py    compute_point, compute_points_batch, compute_facade,
                          Friedlander solvers, Mach stem
  cli.py                  command-line interface (point, compute, batch, gui)
  gui.py                  PySide6 + Matplotlib interface
data/                     19 normalized SI tables derived from UFC 3-340-02
scripts/                  data processing and paper reproduction scripts
tests/                    pytest suite
docs/                     algorithm, data provenance, design decisions
```

## Contributing and support

Bug reports and questions are welcome in
[GitHub Issues](https://github.com/Nyobir/ufc-blast-calculator/issues). Please
include the exact command or code, the output, and the package version.
Pull requests should keep all tests passing (`uv run pytest`) and route any
calculation change through `compute_facade()`. AI coding agents should read
[`AGENTS.md`](AGENTS.md) first.

Contact: Vladyslav Kochkarov, v.d.kochkarov@nuwm.edu.ua

## License

Code is released under the [MIT License](LICENSE). See [`NOTICE`](NOTICE) for the
note on the UFC-derived data tables and the engineering disclaimer.
