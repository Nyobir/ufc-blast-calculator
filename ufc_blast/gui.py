"""
PyVista 3-D GUI for UFC 3-340-02 blast facade visualisation.

Facade convention
-----------------
- Facade surface lives in the z = 0 plane.
- Explosion source is at [0, 0, -R] (in front of the facade).
- Grid points have (dx, dy) offsets from the perpendicular centre, mapped to
  PyVista coordinates (dx, dy, 0).

Usage
-----
    from ufc_blast.gui import launch_gui
    launch_gui(W=100, R=30, Hc=5)
"""

from __future__ import annotations

import math

import numpy as np
import pyvista as pv

from ufc_blast.core.blast_params import BlastPointResult, compute_point, friedlander
from ufc_blast.core.geometry import GridPoint, generate_grid


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_result_map(
    grid_points: list[GridPoint],
    W: float,
) -> dict[tuple[float, float], BlastPointResult]:
    """Compute blast parameters for every grid point.

    Returns a mapping  (dx, dy) → BlastPointResult.
    """
    result_map: dict[tuple[float, float], BlastPointResult] = {}
    for gp in grid_points:
        key = (gp.dx, gp.dy)
        if key not in result_map:
            result_map[key] = compute_point(gp.R_alpha, gp.alpha_deg, W)
    return result_map


def _build_structured_grid(
    grid_points: list[GridPoint],
    result_map: dict[tuple[float, float], BlastPointResult],
) -> tuple[pv.StructuredGrid, np.ndarray]:
    """Build a PyVista StructuredGrid from an ordered rectangular grid.

    ``grid_points`` must be in row-major order (dy changes slowest), as
    produced by ``generate_grid`` in building mode.

    Returns
    -------
    mesh : pv.StructuredGrid
    pr_values : np.ndarray  (flat, same order as mesh points)
    """
    dx_vals = sorted({gp.dx for gp in grid_points})
    dy_vals = sorted({gp.dy for gp in grid_points})
    nx = len(dx_vals)
    ny = len(dy_vals)

    # Build coordinate arrays (z = 0 for entire facade)
    dx_arr = np.array(dx_vals)
    dy_arr = np.array(dy_vals)
    DX, DY = np.meshgrid(dx_arr, dy_arr)   # shape (ny, nx)
    DZ = np.zeros_like(DX)

    mesh = pv.StructuredGrid(DX, DY, DZ)

    # Scalar field: Pr_alpha at each mesh node (row-major dy outer, dx inner)
    pr_values = np.array(
        [result_map[(dx, dy)].Pr_alpha for dy in dy_vals for dx in dx_vals],
        dtype=float,
    )
    return mesh, pr_values


def _build_polydata_surface(
    grid_points: list[GridPoint],
    result_map: dict[tuple[float, float], BlastPointResult],
) -> tuple[pv.PolyData, np.ndarray]:
    """Build a PyVista PolyData cloud for open-field (irregular) grids.

    Returns
    -------
    cloud : pv.PolyData
    pr_values : np.ndarray  (same order as cloud.points)
    """
    pts = np.array([[gp.dx, gp.dy, 0.0] for gp in grid_points], dtype=float)
    cloud = pv.PolyData(pts)
    pr_values = np.array(
        [result_map[(gp.dx, gp.dy)].Pr_alpha for gp in grid_points],
        dtype=float,
    )
    return cloud, pr_values


def _nearest_grid_point(
    clicked_xyz: np.ndarray,
    grid_points: list[GridPoint],
) -> GridPoint:
    """Return the GridPoint closest to a clicked PyVista coordinate."""
    clicked_dx = float(clicked_xyz[0])
    clicked_dy = float(clicked_xyz[1])
    best = min(
        grid_points,
        key=lambda gp: (gp.dx - clicked_dx) ** 2 + (gp.dy - clicked_dy) ** 2,
    )
    return best


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def launch_gui(
    W: float,
    R: float,
    Hc: float,
    width: float | None = None,
    height: float | None = None,
    step: float = 1.0,
) -> None:
    """Launch the interactive PyVista blast visualisation.

    Parameters
    ----------
    W : float
        Charge mass (kg TNT equivalent).
    R : float
        Perpendicular standoff distance from charge to facade (m).
    Hc : float
        Height of burst above the ground (m).
    width : float or None
        Facade width (m) for building mode.  ``None`` → open-field mode.
    height : float or None
        Facade height (m) for building mode.  ``None`` → open-field mode.
    step : float
        Grid spacing (m).  Default 1.0 m.
    """
    # ------------------------------------------------------------------
    # 1. Compute blast parameters for the entire grid
    # ------------------------------------------------------------------
    grid_points = generate_grid(R, W, Hc, width=width, height=height, step=step)
    result_map = _build_result_map(grid_points, W)

    # ------------------------------------------------------------------
    # 2. Build PyVista surface mesh
    # ------------------------------------------------------------------
    building_mode = (width is not None) and (height is not None)

    if building_mode:
        mesh, pr_values = _build_structured_grid(grid_points, result_map)
    else:
        mesh, pr_values = _build_polydata_surface(grid_points, result_map)

    mesh["Pr_alpha_kPa"] = pr_values

    # ------------------------------------------------------------------
    # 3. Assemble PyVista plotter
    # ------------------------------------------------------------------
    plotter = pv.Plotter(window_size=[1200, 800] if (width is None and height is None) else None)
    plotter.set_background("white")

    # Facade surface coloured by Pr_alpha
    plotter.add_mesh(
        mesh,
        scalars="Pr_alpha_kPa",
        cmap="YlOrRd",
        show_edges=False,
        opacity=1.0,
        scalar_bar_args={
            "title": "Pr_α (kPa)",
            "vertical": True,
            "position_x": 0.85,
            "position_y": 0.05,
        },
    )

    # ------------------------------------------------------------------
    # 4. Pressure contour rings (~10 levels)
    # ------------------------------------------------------------------
    pr_min = float(pr_values.min())
    pr_max = float(pr_values.max())

    if pr_max - pr_min > 1e-6:
        try:
            contour_mesh = mesh.contour(
                isosurfaces=10,
                scalars="Pr_alpha_kPa",
                rng=(pr_min, pr_max),
            )
            if contour_mesh.n_points > 0:
                plotter.add_mesh(
                    contour_mesh,
                    color="black",
                    line_width=1.5,
                    opacity=0.6,
                    label="Pressure contours",
                )
        except Exception:
            # Contour extraction can fail on degenerate meshes — skip silently
            pass

    # ------------------------------------------------------------------
    # 5. Explosion source (red sphere at [0, 0, -R])
    # ------------------------------------------------------------------
    sphere_radius = max(step * 0.5, R * 0.02)
    source_sphere = pv.Sphere(radius=sphere_radius, center=[0.0, 0.0, -R])
    plotter.add_mesh(source_sphere, color="red", label="Explosion source")

    # ------------------------------------------------------------------
    # 6. Ray lines from explosion to facade corners (building mode only)
    # ------------------------------------------------------------------
    if building_mode and width is not None and height is not None:
        half_w = width / 2.0
        half_h = height / 2.0
        corners = [
            [ half_w,  half_h, 0.0],
            [-half_w,  half_h, 0.0],
            [ half_w, -half_h, 0.0],
            [-half_w, -half_h, 0.0],
        ]
        source_pt = np.array([0.0, 0.0, -R])
        for corner in corners:
            line = pv.Line(source_pt, corner)
            plotter.add_mesh(line, color="orange", line_width=1.5, opacity=0.7)

    # ------------------------------------------------------------------
    # 7. Perpendicular line (blue) from explosion to facade centre
    # ------------------------------------------------------------------
    perp_line = pv.Line([0.0, 0.0, -R], [0.0, 0.0, 0.0])
    plotter.add_mesh(
        perp_line,
        color="blue",
        line_width=2.5,
        label="Perpendicular standoff",
    )

    # ------------------------------------------------------------------
    # 8. Clickable point picking — show blast params + Friedlander plot
    # ------------------------------------------------------------------
    # Pre-build lookup array for fast nearest-neighbour search
    gp_dx = np.array([gp.dx for gp in grid_points], dtype=float)
    gp_dy = np.array([gp.dy for gp in grid_points], dtype=float)

    def _on_click(point: np.ndarray) -> None:  # noqa: ANN001
        """Callback invoked when user clicks on the facade."""
        import matplotlib.pyplot as plt

        clicked_dx = float(point[0])
        clicked_dy = float(point[1])

        # Find nearest grid point index
        dist2 = (gp_dx - clicked_dx) ** 2 + (gp_dy - clicked_dy) ** 2
        idx = int(np.argmin(dist2))
        gp = grid_points[idx]
        res = result_map[(gp.dx, gp.dy)]

        # Print blast parameters
        print("\n--- Blast parameters at clicked point ---")
        print(f"  dx        = {gp.dx:.2f} m")
        print(f"  dy        = {gp.dy:.2f} m")
        print(f"  R_alpha   = {gp.R_alpha:.3f} m")
        print(f"  alpha     = {gp.alpha_deg:.2f} deg")
        print(f"  Ps0       = {res.Ps0:.3f} kPa")
        print(f"  C_alpha   = {res.C_alpha:.4f}")
        print(f"  Pr_alpha  = {res.Pr_alpha:.3f} kPa")
        print(f"  ir_alpha  = {res.ir_alpha:.6f} kPa·s")
        print(f"  tA        = {res.tA:.6f} s")
        print(f"  t0        = {res.t0:.6f} s")
        print(f"  b         = {res.b:.4f}")
        print("------------------------------------------")

        # Friedlander waveform plot
        t_end = res.tA + 1.5 * res.t0
        t = np.linspace(0.0, t_end, 2000)
        P = friedlander(t, res.Pr_alpha, res.tA, res.t0, res.b)

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(t * 1e3, P, color="steelblue", linewidth=2)
        ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Reflected pressure (kPa)")
        ax.set_title(
            f"Friedlander waveform — dx={gp.dx:.1f} m, dy={gp.dy:.1f} m\n"
            f"Pr_α = {res.Pr_alpha:.1f} kPa, t₀ = {res.t0*1e3:.2f} ms, b = {res.b:.2f}"
        )
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        plt.show(block=False)

    plotter.enable_point_picking(
        callback=_on_click,
        use_mesh=True,
        show_message=True,
        tolerance=0.025,
    )

    # ------------------------------------------------------------------
    # 9. Camera, annotations, and show
    # ------------------------------------------------------------------
    plotter.add_axes()
    plotter.add_title(
        f"UFC 3-340-02  |  W = {W} kg TNT  |  R = {R} m  |  Hc = {Hc} m",
        font_size=10,
    )

    # Position camera to see both facade (z=0) and explosion source (z=-R)
    plotter.camera_position = [
        (0.0, -R * 2.5, R * 0.8),   # camera location
        (0.0, 0.0, -R * 0.2),       # focal point
        (0.0, 0.0, 1.0),            # up vector
    ]

    plotter.show()
