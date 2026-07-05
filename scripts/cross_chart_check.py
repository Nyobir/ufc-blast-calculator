"""
cross_chart_check.py

Cross-chart consistency check of the digitized data layer (manuscript
Section "Data-layer consistency check").

The normally reflected peak pressure is available through two
independently extracted routes:

  direct:   the digitized reflected-pressure curve of UFC Figure 2-7
            (data/free_air/figure_2_7_wpd_pr.csv, WebPlotDigitizer);
  indirect: C_alpha(Ps0, 0 deg) * Ps0(Z), the product of the digitized
            incident-pressure curve (Figure 2-7) and the reflection-
            coefficient table (Figure 2-193).

The check evaluates the indirect route at every Z of the direct curve
and reports relative-difference statistics. Points where the indirect
route leaves the chart-supported range (Ps0 outside the digitized
C_alpha families, or Z outside the incident curve) are excluded and
counted.

Usage:
    python scripts/cross_chart_check.py
"""

import pathlib
import warnings

import numpy as np

from ufc_blast.core import blast_params as bp

PR_CSV = (
    pathlib.Path(__file__).parent.parent
    / "data" / "free_air" / "figure_2_7_wpd_pr.csv"
)


def main() -> None:
    warnings.filterwarnings("ignore")
    bp.load_ufc_tables()
    ps0_tab = bp._tables["ps0"]
    ca_tab = bp._tables["calpha"]

    zs, pr_direct = np.loadtxt(PR_CSV, delimiter=",", skiprows=1, unpack=True)
    order = np.argsort(zs)
    zs, pr_direct = zs[order], pr_direct[order]

    used, excluded = [], 0
    for z, pr in zip(zs, pr_direct):
        ps0, oor_z = ps0_tab.lookup_flagged(float(z))
        ca, oor_fam = ca_tab.lookup_flagged(x=0.0, family=ps0)
        if oor_z or oor_fam:
            excluded += 1
            continue
        used.append((z, pr, ca * ps0))

    arr = np.array(used)
    rel = np.abs(arr[:, 2] - arr[:, 1]) / arr[:, 1] * 100.0

    print(f"Digitized reflected-pressure curve: {len(zs)} points, "
          f"Z = {zs[0]:.3f} - {zs[-1]:.2f} m/kg^(1/3)")
    print(f"Chart-supported for the indirect route: {len(used)} points "
          f"(Z = {arr[0, 0]:.3f} - {arr[-1, 0]:.2f}); "
          f"{excluded} excluded (C_alpha families do not reach the "
          f"required incident pressures at small Z)")
    print()
    print(f"Relative difference |indirect - direct| / direct:")
    print(f"  median : {np.median(rel):.2f} %")
    print(f"  95th   : {np.percentile(rel, 95):.2f} %")
    print(f"  max    : {rel.max():.2f} %")


if __name__ == "__main__":
    main()
