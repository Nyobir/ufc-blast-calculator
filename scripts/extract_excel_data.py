#!/usr/bin/env python3
"""Extract UFC 3-340-02 data from the master Excel workbook into CSV files.

Source: Таблиця UFC_Кінцевий варіант.xlsm, sheet "FreeAirBlast"
Output: data/free_air/*.csv

The Excel workbook contains digitised data from UFC 3-340-02 figures:
  - Figure 2-7: Incident overpressure and incident impulse vs scaled distance Z
  - Figure 2-193: Reflection coefficient Calpha vs angle of incidence
  - Figure 2-194: Scaled reflected impulse iralpha vs angle of incidence

Unit conventions:
  - Z: m/kg^(1/3)
  - Pressure: kPa
  - Impulse: kPa-ms/kg^(1/3)
  - Angles: degrees

The Excel provides SI values for 1D tables (Z, pressure, impulse) directly.
For 2D tables (Calpha, iralpha) only imperial values exist; we convert using
the factor 8.953 for impulse (psi-ms/lb^(1/3) -> kPa-ms/kg^(1/3)), which is
the factor hardcoded in the Excel calculator section.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import openpyxl


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

EXCEL_PATH = Path(
    "/Users/vladyslav/Documents/Phd/Cloud/Корисні штуки/"
    "Таблиця UFC_Кінцевий варіант.xlsm"
)
SHEET_NAME = "FreeAirBlast"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "data" / "free_air"

# Unit conversion constants
PSI_TO_KPA = 6.89476
# Impulse scaling: psi-ms/lb^(1/3) -> kPa-ms/kg^(1/3)
# Derived from the Excel calculator section (hardcoded factor in formulas).
IMPULSE_SCALE_FACTOR = 8.953


# ---------------------------------------------------------------------------
# Table layout in the Excel sheet
# ---------------------------------------------------------------------------

# Figure 2-7: Incident overpressure (Pso vs Z)
# Rows 75-127, Col D = Z (m/kg^1/3), Col F = Pso (kPa)
INCIDENT_PRESSURE_ROWS = (75, 127)
INCIDENT_PRESSURE_Z_COL = 4       # D
INCIDENT_PRESSURE_KPA_COL = 6     # F

# Figure 2-7: Incident impulse (is vs Z)
# Rows 179-175... wait, that is reflected. Actually:
# - Reflected pressure: rows 10-71 (Pr vs Z)
# - Incident pressure: rows 75-127 (Pso vs Z)
# - Reflected impulse: rows 131-175 (ir vs Z)
# - Incident impulse: rows 179-? (is vs Z)
#
# The task asks for incident data only.
INCIDENT_IMPULSE_ROWS = (179, 229)  # Data ends at row 229
INCIDENT_IMPULSE_Z_COL = 4         # D  (m/kg^1/3)
INCIDENT_IMPULSE_IMP_COL = 5       # E  (psi-ms/lb^1/3) — need to convert

# Figure 2-193: Calpha (reflection coefficient) vs angle
# Pressure headers at row 62 in specific columns (triplet structure).
# Data rows: 65-148 (varies per pressure).
# Each pressure occupies 3 columns: angle, calpha, interpolation_result.
CALPHA_HEADER_ROW = 62
CALPHA_DATA_START_ROW = 65
CALPHA_DATA_MAX_ROW = 200  # Will stop at first empty pair
CALPHA_PRESSURE_COLS = [
    # (header_col, pressure_psi) — angle at header_col-1, calpha at header_col
    (12, 0.2), (15, 0.5), (18, 1), (21, 2), (24, 5),
    (27, 10), (30, 20), (33, 30), (36, 50), (39, 70),
    (42, 100), (45, 150), (48, 200), (51, 300), (54, 400),
    (57, 500), (60, 1000), (63, 2000), (66, 3000), (69, 5000),
]

# Figure 2-194: iralpha (scaled reflected impulse) vs angle
# Pressure headers at row 226 in specific columns (triplet structure).
# Data rows: 229-309 (varies per pressure).
IRALPHA_HEADER_ROW = 226
IRALPHA_DATA_START_ROW = 229
IRALPHA_DATA_MAX_ROW = 320  # Will stop at first empty pair
IRALPHA_PRESSURE_COLS = [
    # (header_col, pressure_psi) — angle at header_col-2, iralpha at header_col-1
    (13, 0.7), (16, 1), (19, 1.5), (22, 2), (25, 3),
    (28, 5), (31, 10), (34, 20), (37, 50), (40, 100),
    (43, 200), (46, 400), (49, 700), (52, 1000), (55, 1500),
    (58, 2000), (61, 3000), (64, 4000), (67, 5000), (70, 6000),
    (73, 7000),
]


# ---------------------------------------------------------------------------
# Extraction helpers
# ---------------------------------------------------------------------------

def is_valid_number(value) -> bool:
    """Return True if value is a usable numeric cell (not None, not '-')."""
    if value is None:
        return False
    if isinstance(value, str):
        return False  # covers '-' and any other non-numeric strings
    return isinstance(value, (int, float))


def extract_1d_table(ws, row_range, z_col, val_col, convert_factor=1.0):
    """Extract a 1D table (Z vs value) from the worksheet.

    Parameters
    ----------
    ws : Worksheet
    row_range : tuple (start_row, max_row)
    z_col : int — column index for Z (m/kg^1/3)
    val_col : int — column index for the value
    convert_factor : float — multiplied into the value column

    Returns
    -------
    list of (z_m_kg13, value) tuples
    """
    start, max_row = row_range
    rows = []
    for row in range(start, max_row + 1):
        z = ws.cell(row=row, column=z_col).value
        val = ws.cell(row=row, column=val_col).value
        if not is_valid_number(z):
            break  # end of table
        if not is_valid_number(val):
            continue  # skip rows with '-'
        rows.append((float(z), float(val) * convert_factor))
    return rows


def extract_2d_table(ws, pressure_cols, data_start, data_max, psi_to_kpa=True,
                     value_convert=1.0, angle_offset=-2, value_offset=-1):
    """Extract a 2D table (pressure x angle -> value) in long format.

    For each pressure level, reads (angle, value) pairs from the triplet columns.

    Parameters
    ----------
    ws : Worksheet
    pressure_cols : list of (header_col, pressure_psi)
    data_start : int — first data row
    data_max : int — maximum row to scan
    psi_to_kpa : bool — convert pressure from psi to kPa
    value_convert : float — multiplied into the value column
    angle_offset : int — column offset from header_col to angle column
    value_offset : int — column offset from header_col to value column

    Returns
    -------
    list of (ps0_kpa, angle_deg, value) tuples
    """
    rows = []
    for header_col, p_psi in pressure_cols:
        ps0_kpa = p_psi * PSI_TO_KPA if psi_to_kpa else p_psi
        acol = header_col + angle_offset
        vcol = header_col + value_offset

        for row in range(data_start, data_max + 1):
            angle = ws.cell(row=row, column=acol).value
            value = ws.cell(row=row, column=vcol).value
            if not is_valid_number(angle) and not is_valid_number(value):
                break  # end of data for this pressure
            if not is_valid_number(angle) or not is_valid_number(value):
                continue
            rows.append((round(ps0_kpa, 6), float(angle), float(value) * value_convert))
    return rows


def write_csv(path: Path, header: list[str], rows: list[tuple]):
    """Write rows to a CSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
    print(f"  {path.name}: {len(rows)} data rows")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if not EXCEL_PATH.exists():
        print(f"ERROR: Excel file not found: {EXCEL_PATH}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {EXCEL_PATH.name} ...")
    wb = openpyxl.load_workbook(str(EXCEL_PATH), data_only=True)
    ws = wb[SHEET_NAME]
    print(f"  Sheet: {SHEET_NAME}")

    # ------------------------------------------------------------------
    # 1. Figure 2-7: Incident overpressure
    # ------------------------------------------------------------------
    print("\nExtracting incident pressure (Figure 2-7) ...")
    incident_p = extract_1d_table(
        ws, INCIDENT_PRESSURE_ROWS,
        INCIDENT_PRESSURE_Z_COL, INCIDENT_PRESSURE_KPA_COL,
    )
    write_csv(
        OUTPUT_DIR / "figure_2_7_incident.csv",
        ["z_m_kg13", "pso_kpa"],
        incident_p,
    )

    # ------------------------------------------------------------------
    # 2. Figure 2-7: Incident impulse
    # ------------------------------------------------------------------
    print("\nExtracting incident impulse (Figure 2-7) ...")
    # Col E is in psi-ms/lb^(1/3); convert to kPa-ms/kg^(1/3)
    incident_i = extract_1d_table(
        ws, INCIDENT_IMPULSE_ROWS,
        INCIDENT_IMPULSE_Z_COL, INCIDENT_IMPULSE_IMP_COL,
        convert_factor=IMPULSE_SCALE_FACTOR,
    )
    write_csv(
        OUTPUT_DIR / "figure_2_7_impulse.csv",
        ["z_m_kg13", "is_kpa_ms_kg13"],
        incident_i,
    )

    # ------------------------------------------------------------------
    # 3. Figure 2-193: Calpha
    # ------------------------------------------------------------------
    print("\nExtracting Calpha (Figure 2-193) ...")
    calpha = extract_2d_table(
        ws, CALPHA_PRESSURE_COLS,
        CALPHA_DATA_START_ROW, CALPHA_DATA_MAX_ROW,
        psi_to_kpa=True,
        value_convert=1.0,   # Calpha is dimensionless
        angle_offset=-1,     # angle col = header_col - 1
        value_offset=0,      # value col = header_col itself
    )
    write_csv(
        OUTPUT_DIR / "figure_2_193_calpha.csv",
        ["ps0_kpa", "angle_deg", "calpha"],
        calpha,
    )

    # ------------------------------------------------------------------
    # 4. Figure 2-194: iralpha
    # ------------------------------------------------------------------
    print("\nExtracting iralpha (Figure 2-194) ...")
    iralpha = extract_2d_table(
        ws, IRALPHA_PRESSURE_COLS,
        IRALPHA_DATA_START_ROW, IRALPHA_DATA_MAX_ROW,
        psi_to_kpa=True,
        value_convert=IMPULSE_SCALE_FACTOR,  # psi-ms/lb^(1/3) -> kPa-ms/kg^(1/3)
        angle_offset=-2,     # angle col = header_col - 2
        value_offset=-1,     # value col = header_col - 1
    )
    write_csv(
        OUTPUT_DIR / "figure_2_194_iralpha.csv",
        ["ps0_kpa", "angle_deg", "iralpha_kpa_ms_kg13"],
        iralpha,
    )

    wb.close()

    # ------------------------------------------------------------------
    # Verification
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("VERIFICATION")
    print("=" * 60)

    # Incident pressure row count
    print(f"\nIncident pressure: {len(incident_p)} rows")
    print(f"  Z range: {incident_p[0][0]:.4f} - {incident_p[-1][0]:.4f} m/kg^(1/3)")
    print(f"  P range: {incident_p[-1][1]:.4f} - {incident_p[0][1]:.4f} kPa")

    # Incident impulse row count
    print(f"\nIncident impulse: {len(incident_i)} rows")
    print(f"  Z range: {incident_i[0][0]:.4f} - {incident_i[-1][0]:.4f} m/kg^(1/3)")

    # Calpha verification
    print(f"\nCalpha: {len(calpha)} rows (expect ~1600+)")
    # Find Calpha at angle=0, Ps0=0.2 psi (1.379 kPa)
    ca_02 = [r for r in calpha if abs(r[0] - 0.2 * PSI_TO_KPA) < 0.01 and r[1] == 0.0]
    if ca_02:
        print(f"  Calpha(angle=0, Ps0=0.2psi={0.2*PSI_TO_KPA:.3f}kPa) = {ca_02[0][2]:.4f}  (expect ~2.00)")
    # Find Calpha at angle=0, Ps0=5000 psi
    ca_5000 = [r for r in calpha if abs(r[0] - 5000 * PSI_TO_KPA) < 1 and r[1] == 0.0]
    if ca_5000:
        print(f"  Calpha(angle=0, Ps0=5000psi={5000*PSI_TO_KPA:.1f}kPa) = {ca_5000[0][2]:.4f}  (expect ~12.24)")

    # iralpha row count
    print(f"\niralpha: {len(iralpha)} rows")
    # Spot check: iralpha at angle=0, Ps0=0.7 psi
    ira_07 = [r for r in iralpha if abs(r[0] - 0.7 * PSI_TO_KPA) < 0.01 and r[1] == 0.0]
    if ira_07:
        print(f"  iralpha(angle=0, Ps0=0.7psi) = {ira_07[0][2]:.4f} kPa-ms/kg^(1/3)")
        expected = 2.2811156156267 * IMPULSE_SCALE_FACTOR
        print(f"    Expected: {expected:.4f} (2.2811 * {IMPULSE_SCALE_FACTOR})")

    print("\nDone.")


if __name__ == "__main__":
    main()
