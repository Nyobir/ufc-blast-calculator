"""
Geometry module for UFC 3-340-02 blast wave calculator.

Computes standoff geometry for each grid point relative to the charge:
  - Slant distance R_alpha from charge to surface point
  - Angle of incidence alpha (from the perpendicular line-of-sight)

UFC 3-340-02 applicability limit: the scaled height-of-burst
  Hc / W^(1/3) must exceed 0.397 ft/lb^(1/3) to ensure an air burst
  (not a surface burst). Values at or below this limit fall outside
  the scope of the reflected-pressure charts in the manual.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class GridPoint:
    """Geometry of a single surface point relative to the charge.

    Attributes
    ----------
    dx : float
        Horizontal offset from the perpendicular point on the facade (m).
    dy : float
        Vertical offset from the perpendicular point on the facade (m).
    R_alpha : float
        Slant (true) distance from the charge to this surface point (m).
    alpha_deg : float
        Angle of incidence — angle between the charge-to-point vector and
        the perpendicular line-of-sight, in degrees.  Ranges 0–90°.
    """

    dx: float
    dy: float
    R_alpha: float
    alpha_deg: float


def compute_point_geometry(dx: float, dy: float, R: float) -> GridPoint:
    """Compute the geometry for a single surface point.

    Parameters
    ----------
    dx : float
        Horizontal offset from the perpendicular point (m).
    dy : float
        Vertical offset from the perpendicular point (m).
    R : float
        Perpendicular standoff distance from charge to facade (m).

    Returns
    -------
    GridPoint
        Populated geometry for this surface location.

    Notes
    -----
    The slant distance is::

        R_alpha = sqrt(R^2 + dx^2 + dy^2)

    The angle of incidence is::

        alpha = arccos(R / R_alpha)   [degrees]

    At the perpendicular point (dx=dy=0), R_alpha = R and alpha = 0°.
    """
    R_alpha = math.sqrt(R**2 + dx**2 + dy**2)
    # Guard against floating-point values slightly outside [-1, 1]
    cos_alpha = max(-1.0, min(1.0, R / R_alpha))
    alpha_deg = math.degrees(math.acos(cos_alpha))
    return GridPoint(dx=dx, dy=dy, R_alpha=R_alpha, alpha_deg=alpha_deg)


def generate_grid(
    R: float,
    W: float,
    Hc: float,
    width: float | None = None,
    height: float | None = None,
    step: float = 1.0,
) -> list[GridPoint]:
    """Generate a grid of surface points for blast pressure analysis.

    Two modes are supported:

    **Building mode** (``width`` and ``height`` provided):
        A rectangular grid centered on the perpendicular point, spanning
        ``[-width/2, +width/2]`` horizontally and ``[-height/2, +height/2]``
        vertically, with the given step size.  The number of points along
        each axis is ``floor(dimension / step) + 1``.

    **Open field mode** (``width`` and ``height`` are ``None``):
        The grid expands outward from the perpendicular point in all four
        directions until the angle of incidence exceeds 85°.  Beyond 85°
        the UFC reflected-pressure charts are unreliable (grazing incidence).

    Parameters
    ----------
    R : float
        Perpendicular standoff distance from charge to facade (m).
    W : float
        Charge mass in kg (TNT equivalent) — used only for the applicability
        check.
    Hc : float
        Height of burst above the ground (m) — used for applicability check.
    width : float or None
        Total facade width (m) for building mode.
    height : float or None
        Total facade height (m) for building mode.
    step : float
        Grid spacing (m).  Default 1.0 m.

    Returns
    -------
    list[GridPoint]
        Ordered list of GridPoint objects (row-major: dy changes slowest).

    Raises
    ------
    ValueError
        If ``Hc / W**(1/3)`` does not exceed 0.397, indicating a surface
        or near-surface burst outside UFC 3-340-02 applicability.

    Notes
    -----
    The scaled height-of-burst limit of 0.397 ft/lb^(1/3) is given in
    UFC 3-340-02 Figure 2-15 as the minimum value for air-burst conditions.
    The value is used here in consistent SI-equivalent non-dimensional form:
    since both Hc and W^(1/3) are in SI units the ratio is dimensionally
    identical for the purpose of the threshold check.
    """
    # --- Applicability check -------------------------------------------------
    scaled_hob = Hc / (W ** (1.0 / 3.0))
    if scaled_hob <= 0.397:
        raise ValueError(
            f"Scaled height-of-burst Hc/W^(1/3) = {scaled_hob:.4f} does not "
            f"exceed 0.397 — charge is at or below surface-burst threshold. "
            f"UFC 3-340-02 air-burst charts are not applicable."
        )

    points: list[GridPoint] = []

    if width is not None and height is not None:
        # Building mode: rectangular grid centered on perpendicular point
        half_w = width / 2.0
        half_h = height / 2.0

        # Number of steps along each half-axis (inclusive of centre)
        n_x = round(half_w / step)
        n_y = round(half_h / step)

        for j in range(-n_y, n_y + 1):
            dy = j * step
            for i in range(-n_x, n_x + 1):
                dx = i * step
                points.append(compute_point_geometry(dx, dy, R))
    else:
        # Open field mode: expand until alpha > 85°
        ALPHA_LIMIT = 85.0
        i = 0
        while True:
            dx = i * step
            gp = compute_point_geometry(dx, 0.0, R)
            if i > 0 and gp.alpha_deg > ALPHA_LIMIT:
                break
            # Expand symmetrically in both horizontal and vertical directions
            j = 0
            while True:
                dy = j * step
                gp2 = compute_point_geometry(dx, dy, R)
                if j > 0 and gp2.alpha_deg > ALPHA_LIMIT:
                    break
                # Four-quadrant symmetry (avoid duplicating axes)
                if dx == 0.0 and dy == 0.0:
                    points.append(compute_point_geometry(0.0, 0.0, R))
                elif dx == 0.0:
                    points.append(compute_point_geometry(0.0, dy, R))
                    points.append(compute_point_geometry(0.0, -dy, R))
                elif dy == 0.0:
                    points.append(compute_point_geometry(dx, 0.0, R))
                    points.append(compute_point_geometry(-dx, 0.0, R))
                else:
                    points.append(compute_point_geometry(dx, dy, R))
                    points.append(compute_point_geometry(-dx, dy, R))
                    points.append(compute_point_geometry(dx, -dy, R))
                    points.append(compute_point_geometry(-dx, -dy, R))
                j += 1
            i += 1

    return points
