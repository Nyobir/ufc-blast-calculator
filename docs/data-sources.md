# Data Sources — UFC 3-340-02 Digitized Charts

## Source Excel Workbook

**File**: `Таблиця UFC_Кінцевий варіант.xlsm` (authors' chart-reading workbook, kept outside the repository; its extracted values are in `data/free_air/figure_2_193_calpha.csv` and `figure_2_194_iralpha.csv`)
- 13 sheets, 1160+ rows
- Main sheet: **FreeAirBlast** (434 rows)

### FreeAirBlast Sheet Structure
- **Rows 1-6, A-K**: Calculator inputs (W, R, angle, Hc, Z, Pr, P, ir)
- **Rows 10-51, B-H**: Figure 2-7 Z lookup (53 data points). Columns: Z (ft/lb^1/3), Z (m/kg^1/3), Pr (psi), Pr (kPa), P (psi), P (kPa)
- **Rows 38-60, K-N**: Cα calculator — 20 Ps0 breakpoints in psi: 0.2, 0.5, 1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200, 300, 400, 500, 1000, 2000, 3000, 5000
- **Rows 62-148**: Cα 2D table (Figure 2-193) — 20 Ps0 levels, originally 2191 rows extracted; 14 DPlot tail points were subsequently appended to four truncated families (see below)
- **Rows 226+**: irα 2D table (Figure 2-194) — 22 Ps0 levels x angle points → 1096 rows extracted

### Cα at α=0 deg (Validation Reference)
```
Ps0 [psi]:  0.2   0.5   1     2     5     10    20    30    50    70    100   150   200   300   400   500   1000  2000  3000  5000
Cα(0):      2.0   2.07  2.09  2.12  2.20  2.51  3.00  3.35  4.00  4.45  5.01  5.60  6.00  6.65  6.99  7.80  8.60  10.0  10.80 12.24
```

## WebPlotDigitizer Exports

**Tool**: automeris.io/wpd/
**Raw location**: WebPlotDigitizer project exports kept outside the repository (`DigitalizedUFC/2-7/`); processed values are in `data/free_air/`
**Format**: Semicolon separator, comma decimal, no headers, imperial units

6 files exported from Figure 2-7:
| Curve | Rows | Purpose |
|-------|------|---------|
| t0/W^(1/3) | 167 | Positive phase duration |
| tA/W^(1/3) | 123 | Arrival time |
| Ps0 | 196 | Incident overpressure (cross-validation) |
| Pr | 241 | Reflected pressure (cross-validation) |
| ir | 135 | Reflected impulse |
| is | 176 | Incident impulse |

**Surface burst raw data**: `data/surface/raw/` (Figure 2-15, same format)
**Triple point raw data**: `data/triple_point/` (Figure 2-13, one file per Hc/W^(1/3) level)

## Processing Pipeline

```
WPD export (imperial, semicolons, comma-decimal)
  → scripts/process_wpd_export.py        # Free-air Figure 2-7
  → scripts/process_wpd_surface.py       # Surface Figure 2-15
  → scripts/process_wpd_triple_point.py  # Triple point Figure 2-13
  → data/**/*.csv (SI units, standard CSV)
```

Excel data:
```
Excel .xlsm → scripts/extract_excel_data.py → data/free_air/figure_2_{193,194}_*.csv
Official WBDG DPlot 02_193.grf → scripts/complete_figure_2_193_tails.py →
  14 additional Fig. 2-193 rows for the 5, 10, 30, 70 psi families
```

The additional Fig. 2-193 data come from the UFC 3-340-02 electronic DPlot
package distributed by WBDG (file `02_193.grf`, SHA-256
`dac929c9ac9a7f25823bc47b94d7e0b54466b3c92002ab9d808f6bac6e0e0a56`).
Only points at angles greater than each family's previous maximum were added;
the original 2191 CSV rows were preserved. The separate source GRF remains
outside the repository. Reproduce by passing it explicitly with `--grf` to
`scripts/complete_figure_2_193_tails.py`, together with the original pre-tail
CSV (git history) and a separate `--output` path. This is a DPlot extension,
not a second pass of manual PDF digitization.

The revised ground-evaluated Mach procedure additionally uses official UFC
electronic DPlot numerical series for Figures 2-9 (ground reflected pressure)
and 2-10 (ground reflected scaled impulse). Their SI-normalized outputs are
`data/free_air/figure_2_9_pra.csv` and `figure_2_10_ira.csv`, 2,000 data rows
each. These are distinct from the Figure 2-13 triple-point-height table:
Figure 2-13 locates the Mach-region boundary but does not supply the front
pressure or impulse. Original `.grf` files remain outside this repository;
do not describe the new tables as manual workbook entries or assert that the
original plotting files are bundled.


## Final CSV Files

### Free-Air Burst
| File | Rows | Columns |
|------|------|---------|
| figure_2_7_wpd_ps0.csv | 196 | z_m_kg13, ps0_kpa |
| figure_2_7_wpd_pr.csv | 241 | z_m_kg13, pr_kpa |
| figure_2_7_wpd_ir.csv | 135 | z_m_kg13, ir_kpa_ms_kg13 |
| figure_2_7_wpd_is.csv | 176 | z_m_kg13, is_kpa_ms_kg13 |
| figure_2_7_tA.csv | 123 | z_m_kg13, tA_scaled_ms_kg13 |
| figure_2_7_t0.csv | 167 | z_m_kg13, t0_scaled_ms_kg13 |
| figure_2_7_incident.csv | 53 | z_m_kg13, ps0_kpa (legacy comparison table; not selected by `load_ufc_tables()`) |
| figure_2_7_impulse.csv | 51 | z_m_kg13, is_kpa_ms_kg13 (legacy comparison table; not selected by `load_ufc_tables()`) |
| figure_2_9_pra.csv | 2000 | hc_scaled_m_kg13, angle_deg, pra_kpa (ground Mach pressure) |
| figure_2_10_ira.csv | 2000 | hc_scaled_m_kg13, angle_deg, ira_scaled_kpa_ms_kg13 (ground Mach impulse) |
| figure_2_193_calpha.csv | 2205 (2191 original + 14 DPlot tail) | ps0_kpa, angle_deg, calpha |
| figure_2_194_iralpha.csv | 1096 | ps0_kpa, angle_deg, iralpha_kpa_ms_kg13 |
| figure_2_13_triple_point.csv | 1388 | hc_scaled_m_kg13, rg_scaled_m_kg13, ht_scaled_m_kg13 |

### Surface Burst
| File | Rows | Columns |
|------|------|---------|
| figure_2_15_wpd_ps0.csv | 192 | z_m_kg13, ps0_kpa |
| figure_2_15_tA.csv | 113 | z_m_kg13, tA_scaled_ms_kg13 |
| figure_2_15_t0.csv | 158 | z_m_kg13, t0_scaled_ms_kg13 |
| figure_2_15_wpd_pr.csv | 232 | (extra, not currently loaded) |
| figure_2_15_wpd_ir.csv | 129 | (extra, not currently loaded) |
| figure_2_15_wpd_is.csv | 167 | (extra, not currently loaded) |

The 13 free-air plus 6 surface CSV files above total **19 normalized tables**.
Hatch includes those tables under `ufc_blast/data/` in an installed wheel;
`data/surface/raw/` and `data/triple_point/` are not part of that normalized
wheel payload. Source checkouts, installed wheels and frozen Windows bundles
resolve their data paths separately; the release wheel checker compares all
19 payload hashes against the source tables.

### Column Naming Convention
Units always in column names: `z_m_kg13` = Z in m/kg^(1/3), `ps0_kpa` = Ps0 in kPa.

## Notes
- Figures 2-193 and 2-194 (Cα, irα) are **shared** between air and surface burst — only Ps0/tA/t0 tables differ
- Surface burst Pr/ir/is tables are digitized and stored but **not yet loaded** by the calculator
- Triple point table covers Hc_scaled levels from Figure 2-13 family curves (Hc/W^(1/3) = 1..7 in ft/lb^(1/3), converted to m/kg^(1/3))
- Max Hc_scaled in triple point table: 2.777 m/kg^(1/3) — exceeding this raises ValueError
