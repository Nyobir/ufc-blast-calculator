"""
plot_verification.py

Reproduces UFC 3-340-02 Figure 2-7 (Free Air Burst, SI units) from
digitized WebPlotDigitizer data to verify digitization quality.

Usage:
    python scripts/plot_verification.py
"""

import pathlib
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = pathlib.Path(__file__).parent.parent / "data" / "free_air"
OUTPUT_PATH = pathlib.Path(__file__).parent.parent / "verification_plot.png"


def load(filename: str) -> tuple[np.ndarray, np.ndarray]:
    """Load a two-column CSV (with header) and return (x, y) arrays sorted by x."""
    data = np.loadtxt(DATA_DIR / filename, delimiter=",", skiprows=1)
    x, y = data[:, 0], data[:, 1]
    # Sort by x just in case WPD returned points out of order
    idx = np.argsort(x)
    return x[idx], y[idx]


def main() -> None:
    # -----------------------------------------------------------------------
    # Load digitized data
    # -----------------------------------------------------------------------
    z_ps0,  ps0  = load("figure_2_7_wpd_ps0.csv")
    z_pr,   pr   = load("figure_2_7_wpd_pr.csv")
    z_ir,   ir   = load("figure_2_7_wpd_ir.csv")
    z_is,   is_  = load("figure_2_7_wpd_is.csv")
    z_t0,   t0   = load("figure_2_7_t0.csv")
    z_tA,   tA   = load("figure_2_7_tA.csv")

    # -----------------------------------------------------------------------
    # Build figure
    # -----------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 8), dpi=150)

    # Curve style convention matching original UFC Figure 2-7
    curves = [
        # (x,   y,   color,       linestyle, linewidth, label)
        (z_pr,  pr,  "#CC0000",   "-",       1.8,  r"$P_r$ — Peak reflected pressure (kPa)"),
        (z_ps0, ps0, "#CC00CC",   "-",       1.8,  r"$P_{so}$ — Peak incident pressure (kPa)"),
        (z_ir,  ir,  "#CC0000",   ":",       1.8,  r"$i_r/W^{1/3}$ — Reflected impulse (kPa·ms/kg$^{1/3}$)"),
        (z_is,  is_, "#CC00CC",   "--",      1.8,  r"$i_s/W^{1/3}$ — Incident impulse (kPa·ms/kg$^{1/3}$)"),
        (z_tA,  tA,  "#6600CC",   "--",      1.6,  r"$t_A/W^{1/3}$ — Arrival time (ms/kg$^{1/3}$)"),
        (z_t0,  t0,  "#000066",   "-.",      1.6,  r"$t_o/W^{1/3}$ — Duration (ms/kg$^{1/3}$)"),
    ]

    for x, y, color, ls, lw, label in curves:
        ax.plot(x, y, color=color, linestyle=ls, linewidth=lw, label=label)

    # -----------------------------------------------------------------------
    # Axes — log-log
    # -----------------------------------------------------------------------
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel(
        r"Scaled Distance  $Z = R \, / \, W^{1/3}$  (m/kg$^{1/3}$)",
        fontsize=13,
    )
    ax.set_ylabel("Blast Wave Parameter Value (see legend for units)", fontsize=13)
    ax.set_title(
        "Digitized UFC 3-340-02 Figure 2-7 — Free Air Burst (SI units)",
        fontsize=14,
        fontweight="bold",
        pad=12,
    )

    # Match UFC Figure 2-7 axis ranges (converted from imperial):
    #   X: 0.1–100 ft/lb^(1/3) = 0.0397–39.67 m/kg^(1/3)
    #   Y: 0.005–100000 (shared numeric range across all parameters)
    Z_FACTOR = 0.39685
    ax.set_xlim(0.1 * Z_FACTOR, 100 * Z_FACTOR)
    ax.set_ylim(0.005, 100000)

    # -----------------------------------------------------------------------
    # Grid — major and minor
    # -----------------------------------------------------------------------
    ax.grid(which="major", linestyle="-",  linewidth=0.6, color="#AAAAAA", alpha=0.8)
    ax.grid(which="minor", linestyle=":",  linewidth=0.4, color="#CCCCCC", alpha=0.6)
    ax.minorticks_on()

    # -----------------------------------------------------------------------
    # Legend
    # -----------------------------------------------------------------------
    ax.legend(
        loc="upper right",
        fontsize=10,
        framealpha=0.9,
        edgecolor="#888888",
        title="Blast Parameters",
        title_fontsize=10,
    )

    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=150, bbox_inches="tight")
    print(f"Saved verification plot → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
