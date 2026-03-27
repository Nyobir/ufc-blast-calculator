"""
Process WebPlotDigitizer CSV exports from UFC Figure 2-7 into standardized CSV files.

Source: digitized curves from UFC 3-340-02 Figure 2-7 (free-air burst parameters)
Input format: semicolon-separated, comma as decimal separator (European locale), no headers
Input units:
  - Z: ft/lb^(1/3)
  - Pressure (Ps0, Pr): psi
  - Time (t0, tA): ms/lb^(1/3)
  - Impulse (ir, is): psi-ms/lb^(1/3)

Unit conversions applied:
  - Z [ft/lb^(1/3)] → [m/kg^(1/3)]:
      factor = 0.3048 ft/m / (0.453592 lb/kg)^(1/3) = 0.3048 / 0.76886 ≈ 0.39685
  - Pressure [psi] → [kPa]:
      factor = 6.89476 kPa/psi
  - Time [ms/lb^(1/3)] → [ms/kg^(1/3)]:
      W^(1/3) is in DENOMINATOR, so: factor = 1/(0.453592)^(1/3) = (2.20462)^(1/3) ≈ 1.3015
  - Impulse [psi-ms/lb^(1/3)] → [kPa-ms/kg^(1/3)]:
      factor = 6.89476 / (0.453592)^(1/3) = 6.89476 × (2.20462)^(1/3) ≈ 8.974
      Derivation: (ir/W^1/3)[kPa-ms/kg^1/3] = (ir/W^1/3)[psi-ms/lb^1/3] × (kPa/psi) / (kg/lb)^(1/3)
"""

import csv
import math
from pathlib import Path

# ---------------------------------------------------------------------------
# Conversion constants
# ---------------------------------------------------------------------------
LB_PER_KG = 0.453592          # 1 kg = 1/0.453592 lb => 1 lb = 0.453592 kg
FT_PER_M = 0.3048              # 1 m = 1/0.3048 ft

# Z: ft/lb^(1/3) → m/kg^(1/3)
# Z_SI = Z_imp × ft_per_m / (kg_per_lb)^(1/3)
# But source is ft/lb^(1/3), target is m/kg^(1/3):
#   Z [m/kg^(1/3)] = Z [ft/lb^(1/3)] × (m/ft) / (kg/lb)^(1/3)
#                  = Z × 0.3048 / (0.453592)^(1/3)
CBRT_LB_PER_KG = LB_PER_KG ** (1.0 / 3.0)   # (0.453592)^(1/3) ≈ 0.76886
Z_FACTOR = FT_PER_M / CBRT_LB_PER_KG         # ≈ 0.39685

# Pressure: psi → kPa
KPA_PER_PSI = 6.89476

# Time scaled: ms/lb^(1/3) → ms/kg^(1/3)
# W^(1/3) is in the DENOMINATOR:
# t/W^(1/3) [ms/kg^(1/3)] = t/W^(1/3) [ms/lb^(1/3)] / (kg/lb)^(1/3)
#                          = value / (0.453592)^(1/3) = value × (2.20462)^(1/3)
TIME_FACTOR = 1.0 / CBRT_LB_PER_KG   # ≈ 1.3015

# Impulse scaled: psi-ms/lb^(1/3) → kPa-ms/kg^(1/3)
# W^(1/3) is in the DENOMINATOR:
# i/W^(1/3) [kPa-ms/kg^(1/3)] = i/W^(1/3) [psi-ms/lb^(1/3)] × (kPa/psi) / (kg/lb)^(1/3)
#                               = value × 6.89476 / 0.76886
IMPULSE_FACTOR = KPA_PER_PSI / CBRT_LB_PER_KG  # ≈ 8.974

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SRC_DIR = Path("/Users/vladyslav/Documents/Phd/Cloud/Корисні штуки/DigitalizedUFC/2-7")
OUT_DIR = Path("/Users/vladyslav/Projects/ufc_calculator/data/free_air")


def read_wpd_csv(path: Path) -> list[tuple[float, float]]:
    """Read a WebPlotDigitizer CSV (semicolon sep, comma decimal, no header).

    Returns list of (x, y) float pairs, unsorted.
    """
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=";")
        for line_no, row in enumerate(reader, start=1):
            if not row or all(cell.strip() == "" for cell in row):
                continue
            if len(row) < 2:
                raise ValueError(f"{path.name} line {line_no}: expected 2 columns, got {len(row)}: {row!r}")
            # Replace comma decimal separator with dot
            x_str = row[0].strip().replace(",", ".")
            y_str = row[1].strip().replace(",", ".")
            try:
                x = float(x_str)
                y = float(y_str)
            except ValueError as e:
                raise ValueError(f"{path.name} line {line_no}: cannot parse '{row}': {e}") from e
            rows.append((x, y))
    return rows


def write_csv(path: Path, header: list[str], rows: list[tuple[float, float]]) -> None:
    """Write two-column CSV with header, sorted by first column ascending."""
    sorted_rows = sorted(rows, key=lambda r: r[0])
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in sorted_rows:
            # Write with enough precision to faithfully represent digitized data
            writer.writerow([f"{row[0]:.6f}", f"{row[1]:.6f}"])


def process_all() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Conversion factors:")
    print(f"  Z:       ft/lb^(1/3) → m/kg^(1/3)          × {Z_FACTOR:.6f}")
    print(f"  Pressure: psi → kPa                          × {KPA_PER_PSI:.6f}")
    print(f"  Time:     ms/lb^(1/3) → ms/kg^(1/3)         × {TIME_FACTOR:.6f}")
    print(f"  Impulse:  psi-ms/lb^(1/3) → kPa-ms/kg^(1/3) × {IMPULSE_FACTOR:.6f}")
    print()

    # ------------------------------------------------------------------
    # t0 — arrival time scaled
    # ------------------------------------------------------------------
    raw = read_wpd_csv(SRC_DIR / "t0.csv")
    converted = [(z * Z_FACTOR, v * TIME_FACTOR) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_7_t0.csv", ["z_m_kg13", "t0_scaled_ms_kg13"], converted)
    print(f"t0:  {len(converted)} rows → figure_2_7_t0.csv")
    # Spot-check: smallest Z
    s = sorted(converted, key=lambda r: r[0])
    print(f"     First row: Z={s[0][0]:.4f} m/kg^(1/3), t0={s[0][1]:.4f} ms/kg^(1/3)")
    print(f"     Last row:  Z={s[-1][0]:.4f} m/kg^(1/3), t0={s[-1][1]:.4f} ms/kg^(1/3)")

    # ------------------------------------------------------------------
    # tA — positive phase duration scaled
    # ------------------------------------------------------------------
    raw = read_wpd_csv(SRC_DIR / "tA.csv")
    converted = [(z * Z_FACTOR, v * TIME_FACTOR) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_7_tA.csv", ["z_m_kg13", "tA_scaled_ms_kg13"], converted)
    print(f"tA:  {len(converted)} rows → figure_2_7_tA.csv")
    s = sorted(converted, key=lambda r: r[0])
    print(f"     First row: Z={s[0][0]:.4f} m/kg^(1/3), tA={s[0][1]:.4f} ms/kg^(1/3)")
    print(f"     Last row:  Z={s[-1][0]:.4f} m/kg^(1/3), tA={s[-1][1]:.4f} ms/kg^(1/3)")

    # ------------------------------------------------------------------
    # Ps0 — incident pressure
    # ------------------------------------------------------------------
    raw = read_wpd_csv(SRC_DIR / "Ps0.csv")
    converted = [(z * Z_FACTOR, v * KPA_PER_PSI) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_7_wpd_ps0.csv", ["z_m_kg13", "ps0_kpa"], converted)
    print(f"Ps0: {len(converted)} rows → figure_2_7_wpd_ps0.csv")
    s = sorted(converted, key=lambda r: r[0])
    print(f"     First row: Z={s[0][0]:.4f} m/kg^(1/3), Ps0={s[0][1]:.1f} kPa")
    print(f"     Last row:  Z={s[-1][0]:.4f} m/kg^(1/3), Ps0={s[-1][1]:.1f} kPa")

    # ------------------------------------------------------------------
    # Pr — reflected pressure
    # ------------------------------------------------------------------
    raw = read_wpd_csv(SRC_DIR / "Pr.csv")
    converted = [(z * Z_FACTOR, v * KPA_PER_PSI) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_7_wpd_pr.csv", ["z_m_kg13", "pr_kpa"], converted)
    print(f"Pr:  {len(converted)} rows → figure_2_7_wpd_pr.csv")
    s = sorted(converted, key=lambda r: r[0])
    print(f"     First row: Z={s[0][0]:.4f} m/kg^(1/3), Pr={s[0][1]:.1f} kPa")
    print(f"     Last row:  Z={s[-1][0]:.4f} m/kg^(1/3), Pr={s[-1][1]:.1f} kPa")

    # ------------------------------------------------------------------
    # ir — reflected impulse scaled
    # ------------------------------------------------------------------
    raw = read_wpd_csv(SRC_DIR / "ir.csv")
    converted = [(z * Z_FACTOR, v * IMPULSE_FACTOR) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_7_wpd_ir.csv", ["z_m_kg13", "ir_scaled_kpa_ms_kg13"], converted)
    print(f"ir:  {len(converted)} rows → figure_2_7_wpd_ir.csv")
    s = sorted(converted, key=lambda r: r[0])
    print(f"     First row: Z={s[0][0]:.4f} m/kg^(1/3), ir={s[0][1]:.2f} kPa-ms/kg^(1/3)")
    print(f"     Last row:  Z={s[-1][0]:.4f} m/kg^(1/3), ir={s[-1][1]:.2f} kPa-ms/kg^(1/3)")

    # ------------------------------------------------------------------
    # is — incident impulse scaled
    # ------------------------------------------------------------------
    raw = read_wpd_csv(SRC_DIR / "is.csv")
    converted = [(z * Z_FACTOR, v * IMPULSE_FACTOR) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_7_wpd_is.csv", ["z_m_kg13", "is_scaled_kpa_ms_kg13"], converted)
    print(f"is:  {len(converted)} rows → figure_2_7_wpd_is.csv")
    s = sorted(converted, key=lambda r: r[0])
    print(f"     First row: Z={s[0][0]:.4f} m/kg^(1/3), is={s[0][1]:.2f} kPa-ms/kg^(1/3)")
    print(f"     Last row:  Z={s[-1][0]:.4f} m/kg^(1/3), is={s[-1][1]:.2f} kPa-ms/kg^(1/3)")

    print()
    print("Done. All output files written to:", OUT_DIR)


if __name__ == "__main__":
    process_all()
