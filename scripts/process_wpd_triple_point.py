"""
Process WebPlotDigitizer CSV exports from UFC Figure 2-13 into a single long-format CSV.

Figure 2-13: Scaled Height of Triple Point
  - X-axis: Scaled Horizontal Distance Rg/W^(1/3) [ft/lb^(1/3)]
  - Y-axis: Scaled Height of Triple Point HT/W^(1/3) [ft/lb^(1/3)]
  - Each curve: different Hc/W^(1/3) [ft/lb^(1/3)]

All three dimensions use the same Z_FACTOR conversion (ft/lb^(1/3) → m/kg^(1/3)).
"""

import re
from pathlib import Path

from process_wpd_export import Z_FACTOR, read_wpd_csv, write_csv

SRC_DIR = Path(__file__).parent.parent / "data" / "triple_point"
OUT_FILE = Path(__file__).parent.parent / "data" / "free_air" / "figure_2_13_triple_point.csv"


def process_triple_point() -> None:
    print("Processing UFC Figure 2-13 (triple point height) data")
    print(f"Source: {SRC_DIR}")
    print(f"Output: {OUT_FILE}")
    print(f"Z_FACTOR: {Z_FACTOR:.6f}")
    print()

    # Parse Hc values from filenames like "triple point - 3.5.csv" or "triple-point - 1.csv"
    all_rows: list[tuple[float, float, float]] = []  # (hc_scaled, rg_scaled, ht_scaled)

    for csv_file in sorted(SRC_DIR.glob("*.csv")):
        # Extract the Hc value from filename
        match = re.search(r"(\d+(?:\.\d+)?)\s*\.csv$", csv_file.name)
        if not match:
            print(f"  SKIP: {csv_file.name} (cannot parse Hc value)")
            continue

        hc_imperial = float(match.group(1))  # ft/lb^(1/3)
        hc_si = hc_imperial * Z_FACTOR        # m/kg^(1/3)

        raw = read_wpd_csv(csv_file)
        for rg_imp, ht_imp in raw:
            rg_si = rg_imp * Z_FACTOR
            ht_si = ht_imp * Z_FACTOR
            all_rows.append((hc_si, rg_si, ht_si))

        print(f"  Hc={hc_imperial} ft/lb^(1/3) → {hc_si:.4f} m/kg^(1/3): {len(raw)} points")

    # Sort by (hc, rg) for clean output
    all_rows.sort(key=lambda r: (r[0], r[1]))

    # Write long-format CSV
    with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
        import csv
        writer = csv.writer(f)
        writer.writerow(["hc_scaled_m_kg13", "rg_scaled_m_kg13", "ht_scaled_m_kg13"])
        for hc, rg, ht in all_rows:
            writer.writerow([f"{hc:.6f}", f"{rg:.6f}", f"{ht:.6f}"])

    # Summary
    hc_levels = sorted({r[0] for r in all_rows})
    print(f"\nTotal: {len(all_rows)} rows across {len(hc_levels)} Hc levels")
    print(f"Hc levels (m/kg^(1/3)): {[f'{h:.4f}' for h in hc_levels]}")
    rg_all = [r[1] for r in all_rows]
    ht_all = [r[2] for r in all_rows]
    print(f"Rg range: [{min(rg_all):.4f}, {max(rg_all):.4f}] m/kg^(1/3)")
    print(f"HT range: [{min(ht_all):.4f}, {max(ht_all):.4f}] m/kg^(1/3)")
    print("\nDone.")


if __name__ == "__main__":
    process_triple_point()
