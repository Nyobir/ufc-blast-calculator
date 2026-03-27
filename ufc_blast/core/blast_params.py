"""
Blast parameter calculator for UFC 3-340-02.

Computes reflected overpressure, reflected impulse, arrival time, positive
phase duration, and Friedlander waveform parameters for a given charge mass,
standoff distance, and angle of incidence.

All scaling follows the cube-root (Hopkinson-Cranz) law.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from ufc_blast.core.interpolation import Table1D, Table2D

DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "free_air"


@dataclass
class BlastPointResult:
    """UFC 3-340-02 blast parameters at a single surface point.

    Attributes
    ----------
    Ps0 : float
        Peak incident (side-on) overpressure (kPa).
    C_alpha : float
        Reflection coefficient (dimensionless).
    Pr_alpha : float
        Peak reflected overpressure (kPa).
    ir_alpha : float
        Reflected specific impulse (kPa·s).
    tA : float
        Blast wave arrival time (s).
    t0 : float
        Positive phase duration (s).
    b : float
        Friedlander exponential decay parameter (dimensionless).
    """

    Ps0: float
    C_alpha: float
    Pr_alpha: float
    ir_alpha: float
    tA: float
    t0: float
    b: float


# Module-level table cache — populated once by load_ufc_tables()
_tables: dict = {}


def load_ufc_tables() -> None:
    """Load all UFC 3-340-02 lookup tables from CSV.

    Idempotent: subsequent calls are no-ops once the tables are loaded.
    """
    if _tables:
        return

    _tables["ps0"] = Table1D.from_csv(
        DATA_DIR / "figure_2_7_wpd_ps0.csv",
        "z_m_kg13",
        "ps0_kpa",
    )
    _tables["calpha"] = Table2D.from_csv(
        DATA_DIR / "figure_2_193_calpha.csv",
        "ps0_kpa",
        "angle_deg",
        "calpha",
    )
    _tables["iralpha"] = Table2D.from_csv(
        DATA_DIR / "figure_2_194_iralpha.csv",
        "ps0_kpa",
        "angle_deg",
        "iralpha_kpa_ms_kg13",
    )
    _tables["t0"] = Table1D.from_csv(
        DATA_DIR / "figure_2_7_t0.csv",
        "z_m_kg13",
        "t0_scaled_ms_kg13",
    )
    _tables["tA"] = Table1D.from_csv(
        DATA_DIR / "figure_2_7_tA.csv",
        "z_m_kg13",
        "tA_scaled_ms_kg13",
    )


def compute_point(R_alpha: float, alpha_deg: float, W: float) -> BlastPointResult:
    """Compute blast parameters at a surface point using UFC 3-340-02 charts.

    Parameters
    ----------
    R_alpha : float
        Slant distance from charge to the surface point (m).
    alpha_deg : float
        Angle of incidence (degrees, 0 = normal incidence).
    W : float
        Charge mass (kg TNT equivalent).

    Returns
    -------
    BlastPointResult
        Full set of blast parameters for this point.
    """
    if not _tables:
        load_ufc_tables()

    W_cbrt = W ** (1.0 / 3.0)
    Z = R_alpha / W_cbrt

    # Incident overpressure from Figure 2-7 (WPD digitization)
    Ps0 = _tables["ps0"].lookup(Z)

    # Reflection coefficient from Figure 2-193
    C_alpha = _tables["calpha"].lookup(angle=alpha_deg, ps0=Ps0)
    Pr_alpha = C_alpha * Ps0

    # Reflected specific impulse from Figure 2-194 (kPa·ms/kg^1/3 → kPa·s)
    ir_alpha_scaled = _tables["iralpha"].lookup(angle=alpha_deg, ps0=Ps0)
    ir_alpha = ir_alpha_scaled * W_cbrt / 1000.0

    # Arrival time from Figure 2-7 (ms/kg^1/3 → s)
    tA = _tables["tA"].lookup(Z) * W_cbrt / 1000.0

    # Positive phase duration from Figure 2-7 (ms/kg^1/3 → s)
    t0 = _tables["t0"].lookup(Z) * W_cbrt / 1000.0

    # Friedlander decay parameter
    b = solve_friedlander_b(Pr_alpha, t0, ir_alpha)

    return BlastPointResult(
        Ps0=Ps0,
        C_alpha=C_alpha,
        Pr_alpha=Pr_alpha,
        ir_alpha=ir_alpha,
        tA=tA,
        t0=t0,
        b=b,
    )


def compute_points_batch(
    R_alphas: np.ndarray,
    alpha_degs: np.ndarray,
    W: float,
) -> list[BlastPointResult]:
    """Vectorized computation of blast parameters for many points.

    Much faster than calling compute_point() in a loop — bulk numpy
    interpolation replaces per-point Python calls.

    Parameters
    ----------
    R_alphas : np.ndarray
        Slant distances (m).
    alpha_degs : np.ndarray
        Angles of incidence (degrees).
    W : float
        Charge mass (kg TNT equivalent).

    Returns
    -------
    list[BlastPointResult]
        One result per input point.
    """
    if not _tables:
        load_ufc_tables()

    R_alphas = np.asarray(R_alphas, dtype=float)
    alpha_degs = np.asarray(alpha_degs, dtype=float)
    n = len(R_alphas)

    W_cbrt = W ** (1.0 / 3.0)
    Z = R_alphas / W_cbrt

    # Bulk 1D lookups
    Ps0 = _tables["ps0"].lookup_batch(Z)
    tA_vals = _tables["tA"].lookup_batch(Z) * W_cbrt / 1000.0
    t0_vals = _tables["t0"].lookup_batch(Z) * W_cbrt / 1000.0

    # Bulk 2D lookups
    C_alpha = _tables["calpha"].lookup_batch(alpha_degs, Ps0)
    Pr_alpha = C_alpha * Ps0

    ir_alpha_scaled = _tables["iralpha"].lookup_batch(alpha_degs, Ps0)
    ir_alpha = ir_alpha_scaled * W_cbrt / 1000.0

    # Friedlander b — vectorized Newton's method
    b_vals = _solve_friedlander_b_batch(Pr_alpha, t0_vals, ir_alpha)

    return [
        BlastPointResult(
            Ps0=float(Ps0[i]),
            C_alpha=float(C_alpha[i]),
            Pr_alpha=float(Pr_alpha[i]),
            ir_alpha=float(ir_alpha[i]),
            tA=float(tA_vals[i]),
            t0=float(t0_vals[i]),
            b=float(b_vals[i]),
        )
        for i in range(n)
    ]


def _solve_friedlander_b_batch(
    Pr: np.ndarray, t0: np.ndarray, ir: np.ndarray,
    tol: float = 1e-10, max_iter: int = 100,
) -> np.ndarray:
    """Vectorized Friedlander b solver using Newton's method.

    Solves: ir = (Pr * t0) / b^2 * (b - 1 + exp(-b))  for each point.
    """
    # Initial guess
    b = np.full_like(Pr, 2.0)

    for _ in range(max_iter):
        eb = np.exp(-b)
        f = (Pr * t0) / (b ** 2) * (b - 1.0 + eb) - ir
        # df/db = Pr*t0 * [(-2/b^3)(b-1+eb) + (1/b^2)(1-eb)]
        #       = (Pr*t0/b^3) * [-(2*(b-1+eb)) + b*(1-eb)]
        #       = (Pr*t0/b^3) * [-2*b + 2 - 2*eb + b - b*eb]
        #       = (Pr*t0/b^3) * [-b + 2 - (2+b)*eb]
        df = (Pr * t0 / (b ** 3)) * (-b + 2.0 - (2.0 + b) * eb)

        # Guard against zero derivative
        safe = np.abs(df) > 1e-30
        delta = np.where(safe, f / df, 0.0)
        b = b - delta
        # Keep b positive
        b = np.clip(b, 0.01, 50.0)

        if np.all(np.abs(delta) < tol):
            break

    return b


def solve_friedlander_b(Pr_alpha: float, t0: float, ir_alpha: float) -> float:
    """Solve the Friedlander impulse equation for the decay parameter *b*.

    The analytical impulse of the Friedlander waveform is:

        ir = (Pr * t0) / b^2 * (b - 1 + exp(-b))

    This equation is solved for *b* using Brent's method on the interval
    (0.01, 50).

    Parameters
    ----------
    Pr_alpha : float
        Peak reflected overpressure (kPa).
    t0 : float
        Positive phase duration (s).
    ir_alpha : float
        Target reflected specific impulse (kPa·s).

    Returns
    -------
    float
        Dimensionless Friedlander decay parameter b > 0.
    """

    def _residual(b: float) -> float:
        return (Pr_alpha * t0) / (b ** 2) * (b - 1.0 + math.exp(-b)) - ir_alpha

    return float(brentq(_residual, 0.01, 50.0))


def friedlander(
    t: np.ndarray,
    Pr_alpha: float,
    tA: float,
    t0: float,
    b: float,
) -> np.ndarray:
    """Evaluate the modified Friedlander pressure waveform.

    The waveform is defined as:

        P(t) = Pr_alpha * (1 - tau) * exp(-b * tau)   for tA <= t <= tA + t0
        P(t) = 0                                        otherwise

    where ``tau = (t - tA) / t0``.

    Parameters
    ----------
    t : np.ndarray
        Time array (s).
    Pr_alpha : float
        Peak reflected overpressure (kPa).
    tA : float
        Arrival time (s).
    t0 : float
        Positive phase duration (s).
    b : float
        Friedlander decay parameter (dimensionless).

    Returns
    -------
    np.ndarray
        Pressure values (kPa) at each time in *t*.
    """
    t = np.asarray(t, dtype=float)
    P = np.zeros_like(t)
    mask = (t >= tA) & (t <= tA + t0)
    tau = (t[mask] - tA) / t0
    P[mask] = Pr_alpha * (1.0 - tau) * np.exp(-b * tau)
    return P
