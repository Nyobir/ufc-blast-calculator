"""
benchmark_table.py

Recomputes every number of the manuscript's benchmark anchor case
(Table 1): W = 200 kg TNT, R = 30 m, Hc = 5 m.

Left column:  single-point regular-reflection calculation at alpha = 0.
Right column: the same nominal facade location after the per-column
geometric Mach-region approximation (evaluated at the triple-point
height of the central column).

Usage:
    python scripts/benchmark_table.py
"""

import math
import warnings

from ufc_blast.core.blast_params import (
    compute_facade,
    compute_point,
    load_ufc_tables,
)

W, R, HC = 200.0, 30.0, 5.0


def main() -> None:
    warnings.filterwarnings("ignore")
    load_ufc_tables()

    W_cbrt = W ** (1.0 / 3.0)
    Zp = R / W_cbrt

    point = compute_point(R_alpha=R, alpha_deg=0.0, W=W, burst_type="air")

    facade = compute_facade(R=R, W=W, Hc=HC, width=30.0, height=24.0, step=0.25)
    centre = facade.result_map[(0.0, 0.0)]

    # Triple-point geometry of the central column
    mach_dy = dict(facade.mach_curve)[0.0]
    R_T = math.hypot(R, mach_dy)
    Z_T = R_T / W_cbrt
    alpha_T = math.degrees(math.acos(R / R_T))

    print(f"Benchmark anchor case: W={W:.0f} kg TNT, R={R:.0f} m, Hc={HC:.0f} m")
    print(f"Scaled height of burst Hc/W^(1/3) = {HC / W_cbrt:.3f}  (air burst)")
    print(f"Triple-point height (central column): {mach_dy:.3f} m above burst height")
    print()
    hdr = f"{'Quantity':38s} {'Point (a=0)':>12s} {'Facade centre (Mach)':>22s}"
    print(hdr)
    print("-" * len(hdr))

    def row(label: str, left: str, right: str) -> None:
        print(f"{label:38s} {left:>12s} {right:>22s}")

    row("Nominal scaled distance Z_p (m/kg^1/3)", f"{Zp:.3f}", f"{Zp:.3f}")
    row("Mach evaluation Z_T (m/kg^1/3)", "--", f"{Z_T:.3f}")
    row("  (slant distance / incidence angle)", "--",
        f"{R_T:.2f} m / {alpha_T:.1f} deg")
    row("Incident pressure Ps0 (kPa)", f"{point.Ps0:.2f}", f"{centre.Ps0:.2f}")
    row("Reflection coefficient C_alpha (-)", f"{point.C_alpha:.3f}",
        f"{centre.C_alpha:.3f}")
    row("Reflected pressure Pr_alpha (kPa)", f"{point.Pr_alpha:.2f}",
        f"{centre.Pr_alpha:.2f}")
    row("Reflected impulse ir_alpha (kPa ms)", f"{point.ir_alpha * 1000:.1f}",
        f"{centre.ir_alpha * 1000:.1f}")
    row("Arrival time tA (ms)", f"{point.tA * 1000:.2f}", f"{centre.tA * 1000:.2f}")
    row("Positive-phase duration t0 (ms)", f"{point.t0 * 1000:.2f}",
        f"{centre.t0 * 1000:.2f}")
    row("Friedlander parameter b (-)", f"{point.b:.3f}", f"{centre.b:.3f}")

    # Facade demonstration numbers quoted in the manuscript
    prs = {k: v.Pr_alpha for k, v in facade.result_map.items()}
    k_min = min(prs, key=prs.get)
    n_flagged = sum(1 for v in facade.result_map.values()
                    if "triple_point" in v.extrapolated)
    print()
    print(f"Facade demonstration (30 m x 24 m, step 0.25 m): "
          f"{len(facade.result_map)} points")
    print(f"  max Pr = {max(prs.values()):.2f} kPa, "
          f"min Pr = {prs[k_min]:.2f} kPa at (dx, dy) = {k_min}")
    print(f"  spatial variation = "
          f"{(max(prs.values()) - prs[k_min]) / max(prs.values()) * 100:.1f} %")
    print(f"  triple-point clamp flags: {n_flagged} points "
          f"(columns beyond the digitized Figure 2-13 range)")


if __name__ == "__main__":
    main()
