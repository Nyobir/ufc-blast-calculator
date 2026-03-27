"""
Process WebPlotDigitizer CSV exports from UFC Figure 2-15 into standardized CSV files.

Source: digitized curves from UFC 3-340-02 Figure 2-15 (surface burst parameters)
Input format: semicolon-separated, comma as decimal separator (European locale), no headers
Input units: same as Figure 2-7 (imperial)
  - Z: ft/lb^(1/3)
  - Pressure (Ps0, Pr): psi
  - Time (t0, tA): ms/lb^(1/3)
  - Impulse (ir, is): psi-ms/lb^(1/3)
"""

from pathlib import Path

# Reuse the reading/writing and conversion constants from the free air script
from process_wpd_export import (
    Z_FACTOR,
    KPA_PER_PSI,
    TIME_FACTOR,
    IMPULSE_FACTOR,
    read_wpd_csv,
    write_csv,
)

SRC_DIR = Path(__file__).parent.parent / "data" / "surface" / "raw"
OUT_DIR = Path(__file__).parent.parent / "data" / "surface"


def process_surface() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Processing UFC Figure 2-15 (surface burst) data")
    print(f"Source: {SRC_DIR}")
    print(f"Output: {OUT_DIR}")
    print()

    # Ps0
    raw = read_wpd_csv(SRC_DIR / "Surface blast - Ps0.csv")
    converted = [(z * Z_FACTOR, v * KPA_PER_PSI) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_15_wpd_ps0.csv", ["z_m_kg13", "ps0_kpa"], converted)
    s = sorted(converted, key=lambda r: r[0])
    print(f"Ps0: {len(converted)} rows → figure_2_15_wpd_ps0.csv")
    print(f"     Z range: [{s[0][0]:.4f}, {s[-1][0]:.4f}] m/kg^(1/3)")
    print(f"     Ps0 range: [{s[-1][1]:.2f}, {s[0][1]:.2f}] kPa")

    # tA
    raw = read_wpd_csv(SRC_DIR / "Surface blast - ta.csv")
    converted = [(z * Z_FACTOR, v * TIME_FACTOR) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_15_tA.csv", ["z_m_kg13", "tA_scaled_ms_kg13"], converted)
    s = sorted(converted, key=lambda r: r[0])
    print(f"tA:  {len(converted)} rows → figure_2_15_tA.csv")
    print(f"     Z range: [{s[0][0]:.4f}, {s[-1][0]:.4f}] m/kg^(1/3)")
    print(f"     tA range: [{s[0][1]:.4f}, {s[-1][1]:.4f}] ms/kg^(1/3)")

    # t0
    raw = read_wpd_csv(SRC_DIR / "Surface blast - t0.csv")
    converted = [(z * Z_FACTOR, v * TIME_FACTOR) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_15_t0.csv", ["z_m_kg13", "t0_scaled_ms_kg13"], converted)
    s = sorted(converted, key=lambda r: r[0])
    print(f"t0:  {len(converted)} rows → figure_2_15_t0.csv")
    print(f"     Z range: [{s[0][0]:.4f}, {s[-1][0]:.4f}] m/kg^(1/3)")
    print(f"     t0 range: [{s[0][1]:.4f}, {s[-1][1]:.4f}] ms/kg^(1/3)")

    # Pr (reflected pressure — not currently used but good to have)
    raw = read_wpd_csv(SRC_DIR / "Surface blast - Pr.csv")
    converted = [(z * Z_FACTOR, v * KPA_PER_PSI) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_15_wpd_pr.csv", ["z_m_kg13", "pr_kpa"], converted)
    s = sorted(converted, key=lambda r: r[0])
    print(f"Pr:  {len(converted)} rows → figure_2_15_wpd_pr.csv")
    print(f"     Z range: [{s[0][0]:.4f}, {s[-1][0]:.4f}] m/kg^(1/3)")

    # ir (reflected impulse)
    raw = read_wpd_csv(SRC_DIR / "Surface Blast - ir.csv")
    converted = [(z * Z_FACTOR, v * IMPULSE_FACTOR) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_15_wpd_ir.csv", ["z_m_kg13", "ir_scaled_kpa_ms_kg13"], converted)
    s = sorted(converted, key=lambda r: r[0])
    print(f"ir:  {len(converted)} rows → figure_2_15_wpd_ir.csv")
    print(f"     Z range: [{s[0][0]:.4f}, {s[-1][0]:.4f}] m/kg^(1/3)")

    # is (incident impulse)
    raw = read_wpd_csv(SRC_DIR / "Surface blast - is.csv")
    converted = [(z * Z_FACTOR, v * IMPULSE_FACTOR) for z, v in raw]
    write_csv(OUT_DIR / "figure_2_15_wpd_is.csv", ["z_m_kg13", "is_scaled_kpa_ms_kg13"], converted)
    s = sorted(converted, key=lambda r: r[0])
    print(f"is:  {len(converted)} rows → figure_2_15_wpd_is.csv")
    print(f"     Z range: [{s[0][0]:.4f}, {s[-1][0]:.4f}] m/kg^(1/3)")

    print()
    print("Done.")


if __name__ == "__main__":
    process_surface()
