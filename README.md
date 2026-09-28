# UFC Blast Wave Calculator

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)
![Tests](https://img.shields.io/badge/tests-113%20passing-brightgreen)

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
uv sync --extra dev

# or plain pip
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

A pre-built Windows bundle (no Python needed) is attached to each
[GitHub release](https://github.com/Nyobir/ufc-blast-calculator/releases).

## Quick start

### Command line

```bash
# One point: 200 kg TNT, 30 m standoff, burst 5 m above ground, normal incidence
ufc-blast point --W 200 --R 30 --Hc 5

# Whole facade: 20 m x 10 m, 0.25 m grid, written to CSV
ufc-blast compute --W 200 --R 30 --Hc 5 --width 20 --height 10 --step 0.25 -o facade.csv

# Several angles of incidence at one standoff
ufc-blast batch --W 200 --R 30 --Hc 5 --alphas 0,15,30,45,60

# Interactive GUI: pressure map + Friedlander curve for any clicked point
ufc-blast gui --W 200 --R 30 --Hc 5 --width 20 --height 10 --step 0.25
```

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
  → Mach stem         air bursts only, triple-point height from Fig. 2-13
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

`data/` holds 17 SI-normalized CSV tables taken from the UFC 3-340-02 design
charts:

| Chart | Content | Files |
|---|---|---|
| Fig. 2-7 | free-air burst: P<sub>s0</sub>, P<sub>r</sub>, i<sub>s</sub>, i<sub>r</sub>, t<sub>A</sub>, t<sub>0</sub> vs. scaled distance | `data/free_air/figure_2_7_*.csv` |
| Fig. 2-13 | triple-point height (Mach stem) | `data/free_air/figure_2_13_triple_point.csv`, `data/triple_point/` |
| Fig. 2-15 | surface burst: same parameters as Fig. 2-7 | `data/surface/figure_2_15_*.csv` |
| Fig. 2-193 | reflection coefficient C<sub>α</sub>(α, P<sub>s0</sub>) | `data/free_air/figure_2_193_calpha.csv` |
| Fig. 2-194 | reflected scaled impulse i<sub>rα</sub>(α, P<sub>s0</sub>) | `data/free_air/figure_2_194_iralpha.csv` |

Curves were digitized with WebPlotDigitizer; the angle-dependent tables were
transcribed from a chart-reading workbook. The processing scripts in
`scripts/` convert the imperial chart units to SI. Provenance for every table
is in [`docs/data-sources.md`](docs/data-sources.md).

UFC 3-340-02 is approved for public release with unlimited distribution. The
repository stores only the extracted numerical values, not the chart artwork.

## Verification

```bash
uv run pytest tests/ -q        # 113 tests
```

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

The version used in the paper is tagged
[`manuscript-v1.2`](https://github.com/Nyobir/ufc-blast-calculator/tree/manuscript-v1.2).

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
data/                     17 SI tables digitized from UFC 3-340-02
scripts/                  data processing and paper reproduction scripts
tests/                    pytest suite
docs/                     algorithm, data provenance, design decisions
```

## Citation

If you use this software, please cite it. The citation metadata is in
[`CITATION.cff`](CITATION.cff); GitHub shows it under **"Cite this repository"**.
The software paper is currently under review; until it is published please
also cite:

> P.M. Martyniuk, V.D. Kochkarov. Two-dimensional modeling of explosion
> impulses on a structure using the discrete element method. *Mathematical
> Modeling and Computing* 12(4) (2025) 1157–1168.
> https://doi.org/10.23939/mmc2025.04.1157

## Contributing and support

Bug reports and questions are welcome in
[GitHub Issues](https://github.com/Nyobir/ufc-blast-calculator/issues). Please
include the exact command or code, the output, and the package version.
Pull requests should keep all tests passing (`uv run pytest`) and route any
calculation change through `compute_facade()`. AI coding agents should read
[`AGENTS.md`](AGENTS.md) first.

Contact: Vladyslav Kochkarov, v.d.kochkarov@nuwm.edu.ua

## License

Code is released under the [MIT License](LICENSE). See the license file for the
note on the UFC-derived data tables and the engineering disclaimer.
