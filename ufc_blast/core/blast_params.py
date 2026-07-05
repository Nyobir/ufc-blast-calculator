"""
Blast parameter calculator for UFC 3-340-02.

Computes reflected overpressure, reflected impulse, arrival time, positive
phase duration, and Friedlander waveform parameters for a given charge mass,
standoff distance, and angle of incidence.

All scaling follows the cube-root (Hopkinson-Cranz) law.
"""

from __future__ import annotations

import math
import sys
import warnings
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from ufc_blast.core.interpolation import Table1D, Table2D

# When frozen by PyInstaller, data files live next to the executable
# (sys._MEIPASS for --onefile, or the exe directory for --onedir).
if getattr(sys, "frozen", False):
    _BASE = Path(sys._MEIPASS)
else:
    _BASE = Path(__file__).resolve().parent.parent.parent

DATA_DIR = _BASE / "data" / "free_air"
SURFACE_DATA_DIR = _BASE / "data" / "surface"


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
    extrapolated : set[str]
        Names of 2D tables where family (Ps0) was outside the data range
        and the result was extrapolated.  Empty set when all lookups are
        within range.
    """

    Ps0: float
    C_alpha: float
    Pr_alpha: float
    ir_alpha: float
    tA: float
    t0: float
    b: float
    extrapolated: set[str] = field(default_factory=set)


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
        family_col="ps0_kpa",
        x_col="angle_deg",
        value_col="calpha",
    )
    _tables["iralpha"] = Table2D.from_csv(
        DATA_DIR / "figure_2_194_iralpha.csv",
        family_col="ps0_kpa",
        x_col="angle_deg",
        value_col="iralpha_kpa_ms_kg13",
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

    # Surface burst tables (Figure 2-15)
    _tables["surface_ps0"] = Table1D.from_csv(
        SURFACE_DATA_DIR / "figure_2_15_wpd_ps0.csv",
        "z_m_kg13",
        "ps0_kpa",
    )
    _tables["surface_tA"] = Table1D.from_csv(
        SURFACE_DATA_DIR / "figure_2_15_tA.csv",
        "z_m_kg13",
        "tA_scaled_ms_kg13",
    )
    _tables["surface_t0"] = Table1D.from_csv(
        SURFACE_DATA_DIR / "figure_2_15_t0.csv",
        "z_m_kg13",
        "t0_scaled_ms_kg13",
    )

    # Triple point height table (Figure 2-13)
    _tables["triple_point"] = Table2D.from_csv(
        DATA_DIR / "figure_2_13_triple_point.csv",
        family_col="hc_scaled_m_kg13",
        x_col="rg_scaled_m_kg13",
        value_col="ht_scaled_m_kg13",
    )


def compute_point(R_alpha: float, alpha_deg: float, W: float, burst_type: str = "air") -> BlastPointResult:
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

    # Select tables based on burst type
    ps0_table = _tables["surface_ps0"] if burst_type == "surface" else _tables["ps0"]
    tA_table = _tables["surface_tA"] if burst_type == "surface" else _tables["tA"]
    t0_table = _tables["surface_t0"] if burst_type == "surface" else _tables["t0"]

    extrapolated: set[str] = set()

    # Incident overpressure
    Ps0, oor_ps0 = ps0_table.lookup_flagged(Z)
    if oor_ps0:
        extrapolated.add("ps0")

    # Reflection coefficient from Figure 2-193 (shared across burst types)
    C_alpha, oor_calpha = _tables["calpha"].lookup_flagged(x=alpha_deg, family=Ps0)
    if oor_calpha:
        extrapolated.add("calpha")
    Pr_alpha = C_alpha * Ps0

    # Reflected specific impulse from Figure 2-194 (kPa·ms/kg^1/3 → kPa·s)
    ir_alpha_scaled, oor_iralpha = _tables["iralpha"].lookup_flagged(x=alpha_deg, family=Ps0)
    if oor_iralpha:
        extrapolated.add("iralpha")
    ir_alpha = ir_alpha_scaled * W_cbrt / 1000.0

    # Arrival time (ms/kg^1/3 → s)
    tA_scaled, oor_tA = tA_table.lookup_flagged(Z)
    if oor_tA:
        extrapolated.add("tA")
    tA = tA_scaled * W_cbrt / 1000.0

    # Positive phase duration (ms/kg^1/3 → s)
    t0_scaled, oor_t0 = t0_table.lookup_flagged(Z)
    if oor_t0:
        extrapolated.add("t0")
    t0 = t0_scaled * W_cbrt / 1000.0

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
        extrapolated=extrapolated,
    )


def compute_points_batch(
    R_alphas: np.ndarray,
    alpha_degs: np.ndarray,
    W: float,
    burst_type: str = "air",
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

    # Select tables based on burst type
    ps0_table = _tables["surface_ps0"] if burst_type == "surface" else _tables["ps0"]
    tA_table = _tables["surface_tA"] if burst_type == "surface" else _tables["tA"]
    t0_table = _tables["surface_t0"] if burst_type == "surface" else _tables["t0"]

    # Bulk 1D lookups (flagged: out-of-range Z clamps to the chart boundary)
    Ps0, oor_ps0 = ps0_table.lookup_batch_flagged(Z)
    tA_scaled, oor_tA = tA_table.lookup_batch_flagged(Z)
    tA_vals = tA_scaled * W_cbrt / 1000.0
    t0_scaled, oor_t0 = t0_table.lookup_batch_flagged(Z)
    t0_vals = t0_scaled * W_cbrt / 1000.0

    # Bulk 2D lookups (flagged to track out-of-range extrapolation)
    C_alpha, oor_calpha = _tables["calpha"].lookup_batch_flagged(xs=alpha_degs, families=Ps0)
    Pr_alpha = C_alpha * Ps0

    ir_alpha_scaled, oor_iralpha = _tables["iralpha"].lookup_batch_flagged(xs=alpha_degs, families=Ps0)
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
            extrapolated=(
                ({"ps0"} if oor_ps0[i] else set())
                | ({"calpha"} if oor_calpha[i] else set())
                | ({"iralpha"} if oor_iralpha[i] else set())
                | ({"tA"} if oor_tA[i] else set())
                | ({"t0"} if oor_t0[i] else set())
            ),
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

    b_lo, b_hi = 0.01, 50.0
    r_lo, r_hi = _residual(b_lo), _residual(b_hi)
    if r_lo * r_hi > 0.0:
        # No root inside the bracket: the chart-derived (Pr, t0, ir) triple
        # is not impulse-consistent for any admissible b. This only occurs
        # for chart-boundary-clamped (flagged) inputs; return the bracket
        # end with the smaller residual, mirroring the clipping behaviour
        # of the vectorized Newton solver.
        warnings.warn(
            f"No Friedlander b in [{b_lo}, {b_hi}] matches the chart "
            f"impulse (Pr={Pr_alpha:.4g} kPa, t0={t0:.4g} s, "
            f"ir={ir_alpha:.4g} kPa s); clipping to the bracket boundary.",
            UserWarning,
            stacklevel=2,
        )
        return b_lo if abs(r_lo) <= abs(r_hi) else b_hi

    return float(brentq(_residual, b_lo, b_hi))


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


@dataclass
class FacadeResult:
    """Complete facade computation result.

    Single source of truth for the blast parameter pipeline — used by
    both the CLI and the GUI to guarantee identical results.
    """

    burst_type: str
    scaled_hob: float
    grid_points: list  # list[GridPoint] from geometry module
    result_map: dict[tuple[float, float], BlastPointResult]
    mach_curve: list[tuple[float, float]]


def compute_facade(
    R: float,
    W: float,
    Hc: float,
    width: float | None = None,
    height: float | None = None,
    step: float = 1.0,
) -> FacadeResult:
    """Run the full blast-parameter pipeline for a facade.

    This is the **single computation path** shared by the CLI and the GUI.
    Any change to the calculation logic must happen here so that both
    interfaces stay in sync.

    Steps
    -----
    1. Determine burst type from Hc / W^(1/3).
    2. Generate the facade grid.
    3. Batch-compute blast parameters (air or surface tables).
    4. Apply Mach stem correction (air bursts only).

    Parameters
    ----------
    R : float   Perpendicular standoff (m).
    W : float   Charge mass (kg TNT).
    Hc : float  Height of burst (m).
    width, height : float or None   Facade dimensions (building mode).
    step : float   Grid spacing (m).

    Returns
    -------
    FacadeResult
        Contains burst_type, grid_points, result_map, and mach_curve.
    """
    from ufc_blast.core.geometry import determine_burst_type, generate_grid

    if not _tables:
        load_ufc_tables()

    burst_type, scaled_hob = determine_burst_type(Hc, W)

    # Validate that scaled Hc falls within the triple point table range
    # for air bursts (Mach stem correction requires this data)
    if burst_type == "air":
        triple_table = _tables["triple_point"]
        max_hc_scaled = float(triple_table.family_levels[-1])
        if scaled_hob > max_hc_scaled:
            raise ValueError(
                f"No data available: Hc/W^(1/3) = {scaled_hob:.3f} exceeds "
                f"triple point table maximum ({max_hc_scaled:.3f} m/kg^(1/3)). "
                f"Reduce Hc or increase W."
            )

    grid_points = generate_grid(R, W, Hc, width=width, height=height, step=step)

    R_alphas = np.array([gp.R_alpha for gp in grid_points])
    alpha_degs = np.array([gp.alpha_deg for gp in grid_points])
    results = compute_points_batch(R_alphas, alpha_degs, W, burst_type=burst_type)

    result_map: dict[tuple[float, float], BlastPointResult] = {}
    for gp, res in zip(grid_points, results):
        key = (gp.dx, gp.dy)
        if key not in result_map:
            result_map[key] = res

    if burst_type == "air":
        result_map, mach_curve = apply_mach_stem(
            grid_points, result_map, W=W, Hc=Hc, R=R,
        )
    else:
        mach_curve = []

    return FacadeResult(
        burst_type=burst_type,
        scaled_hob=scaled_hob,
        grid_points=grid_points,
        result_map=result_map,
        mach_curve=mach_curve,
    )


def apply_mach_stem(
    grid_points: list,
    result_map: dict[tuple[float, float], BlastPointResult],
    W: float,
    Hc: float,
    R: float,
) -> tuple[dict[tuple[float, float], BlastPointResult], list[tuple[float, float]]]:
    """Apply Mach stem correction to grid results.

    For each facade column (unique dx), compute the triple point height from
    UFC Figure 2-13. All grid points below that height get their blast
    parameters overwritten with the values at the triple point height
    (uniform Mach stem pressure).

    Parameters
    ----------
    grid_points : list[GridPoint]
        Grid points on the facade.
    result_map : dict
        Mapping (dx, dy) -> BlastPointResult (modified in-place).
    W : float
        Charge mass (kg TNT equivalent).
    Hc : float
        Height of burst above ground (m).
    R : float
        Perpendicular standoff distance (m).

    Returns
    -------
    result_map : dict
        Modified result map with Mach-corrected values below triple point.
    mach_curve : list[tuple[float, float]]
        List of (dx, mach_dy) points defining the Mach line on the facade.
    """
    if not _tables:
        load_ufc_tables()

    triple_table = _tables["triple_point"]
    W_cbrt = W ** (1.0 / 3.0)
    Hc_scaled = Hc / W_cbrt

    # Group grid points by column (dx)
    columns: dict[float, list[float]] = defaultdict(list)
    for gp in grid_points:
        columns[gp.dx].append(gp.dy)

    # Sort dy values in each column
    for dx in columns:
        columns[dx].sort()

    mach_curve: list[tuple[float, float]] = []

    for dx in sorted(columns.keys()):
        dy_list = columns[dx]

        # Horizontal ground distance from charge epicenter
        Rg = math.sqrt(R ** 2 + dx ** 2)
        Rg_scaled = Rg / W_cbrt

        # Look up triple point height. An out-of-range scaled ground
        # distance clamps to the digitized Figure 2-13 boundary; the whole
        # column is then flagged, because the clamped triple-point height
        # decides the Mach/regular classification of every point in it.
        # A scaled burst height outside the digitized family range means
        # there is no triple-point data at all — no Mach treatment.
        fam_min = float(triple_table.family_levels[0])
        fam_max = float(triple_table.family_levels[-1])
        if not (fam_min <= Hc_scaled <= fam_max):
            continue
        HT_scaled, tp_clamped = triple_table.lookup_flagged(
            x=Rg_scaled, family=Hc_scaled
        )

        if tp_clamped:
            for dy in dy_list:
                result_map[(dx, dy)].extrapolated.add("triple_point")

        HT = HT_scaled * W_cbrt
        # Convert to facade coordinate: dy=0 is at burst height, ground at dy=-Hc
        mach_dy = HT - Hc

        mach_curve.append((dx, mach_dy))

        dy_arr = np.array(dy_list)
        below_mask = dy_arr < mach_dy

        if not np.any(below_mask):
            # Mach line is below all grid points — no correction needed
            continue

        # Evaluate the merged-front parameters exactly at the triple-point
        # height (dx, mach_dy) through the pointwise pipeline, so the Mach
        # correction is independent of the facade extent and grid step.
        R_alpha_mach = math.sqrt(R ** 2 + dx ** 2 + mach_dy ** 2)
        alpha_mach = math.degrees(math.acos(min(1.0, R / R_alpha_mach)))
        try:
            mach_result = compute_point(
                R_alpha_mach, alpha_mach, W, burst_type="air",
            )
        except ValueError:
            # Triple-point height lies outside the chart-supported range —
            # leave the regular-reflection values for this column.
            continue

        # Overwrite all points below the Mach line, carrying over the
        # column's triple-point clamp flag if it was set.
        if tp_clamped:
            mach_result.extrapolated.add("triple_point")
        for dy in dy_list:
            if dy < mach_dy:
                result_map[(dx, dy)] = mach_result

    return result_map, mach_curve
