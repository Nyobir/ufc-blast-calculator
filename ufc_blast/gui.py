"""
PySide6 + matplotlib GUI for UFC 3-340-02 blast facade visualisation.

Layout
------
- Top bar   : six QDoubleSpinBox inputs (W, R, Hc, width, height, step)
              plus a "Recalculate" button.
- Left  65% : matplotlib FigureCanvas — filled contour plot (contourf) of
              Pr_alpha on the facade with labelled contour lines and a
              building outline rectangle.
- Right 35% : matplotlib FigureCanvas for the Friedlander P(t) waveform
              at a clicked point, plus QLabel rows showing all blast
              parameters.
- Bottom    : QStatusBar — reports calculation time and point count.

Facade coordinate convention
-----------------------------
dx : horizontal offset from the perpendicular centre (m)  — X axis
dy : vertical offset from the perpendicular centre (m)    — Y axis
"""

from __future__ import annotations

import time
from typing import Optional

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QApplication,
    QDoubleSpinBox,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from ufc_blast.core.blast_params import (
    BlastPointResult,
    apply_mach_stem,
    compute_point,
    compute_points_batch,
    friedlander,
    load_ufc_tables,
)
from ufc_blast.core.geometry import GridPoint, determine_burst_type, generate_grid


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def _compute_results(
    grid_points: list[GridPoint],
    W: float,
    burst_type: str = "air",
) -> dict[tuple[float, float], BlastPointResult]:
    """Compute blast parameters for every (dx, dy) grid point using vectorized batch.

    Returns a mapping  (dx, dy) → BlastPointResult.
    """
    R_alphas = np.array([gp.R_alpha for gp in grid_points])
    alpha_degs = np.array([gp.alpha_deg for gp in grid_points])
    results = compute_points_batch(R_alphas, alpha_degs, W, burst_type=burst_type)

    result_map: dict[tuple[float, float], BlastPointResult] = {}
    for gp, res in zip(grid_points, results):
        key = (gp.dx, gp.dy)
        if key not in result_map:
            result_map[key] = res
    return result_map


class _ComputeWorker(QThread):
    """Background thread for blast computation.

    Results are stored as attributes rather than passed through Signal args,
    because PySide6 cannot serialise arbitrary Python dicts across threads.
    """

    finished = Signal()

    def __init__(self, R, W, Hc, width, height, step):
        super().__init__()
        self.R = R
        self.W = W
        self.Hc = Hc
        self.width = width
        self.height = height
        self.step = step
        self.error: str | None = None
        self.grid_points: list | None = None
        self.result_map: dict | None = None
        self.mach_curve: list | None = None
        self.burst_type: str = "air"
        self.scaled_hob: float = 0.0
        self.elapsed: float = 0.0

    def run(self):
        t0 = time.perf_counter()
        try:
            self.burst_type, self.scaled_hob = determine_burst_type(self.Hc, self.W)
            self.grid_points = generate_grid(
                self.R, self.W, self.Hc,
                width=self.width, height=self.height, step=self.step,
            )
            self.result_map = _compute_results(self.grid_points, self.W, burst_type=self.burst_type)

            # Apply Mach stem correction for air bursts only
            if self.burst_type == "air":
                self.result_map, self.mach_curve = apply_mach_stem(
                    self.grid_points, self.result_map,
                    W=self.W, Hc=self.Hc, R=self.R,
                )
            else:
                self.mach_curve = []

            self.elapsed = time.perf_counter() - t0
        except Exception as exc:
            self.error = str(exc)
        self.finished.emit()


def _build_meshgrid(
    grid_points: list[GridPoint],
    result_map: dict[tuple[float, float], BlastPointResult],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reshape the flat grid list into 2-D arrays for matplotlib contourf.

    ``generate_grid`` in building mode produces a row-major list where dy
    varies slowest (outer loop) and dx varies fastest (inner loop), which
    matches numpy's default C order.

    Returns
    -------
    X : np.ndarray, shape (ny, nx)
    Y : np.ndarray, shape (ny, nx)
    Pr : np.ndarray, shape (ny, nx)   — Pr_alpha in kPa
    """
    dx_vals = sorted({gp.dx for gp in grid_points})
    dy_vals = sorted({gp.dy for gp in grid_points})
    nx, ny = len(dx_vals), len(dy_vals)

    dx_arr = np.array(dx_vals)
    dy_arr = np.array(dy_vals)
    X, Y = np.meshgrid(dx_arr, dy_arr)   # shape (ny, nx)

    Pr = np.array(
        [result_map[(dx, dy)].Pr_alpha for dy in dy_vals for dx in dx_vals],
        dtype=float,
    ).reshape(ny, nx)

    return X, Y, Pr


def _nearest_point(
    click_dx: float,
    click_dy: float,
    grid_points: list[GridPoint],
) -> GridPoint:
    """Return the GridPoint closest to the clicked matplotlib coordinate."""
    gp_dx = np.array([gp.dx for gp in grid_points], dtype=float)
    gp_dy = np.array([gp.dy for gp in grid_points], dtype=float)
    dist2 = (gp_dx - click_dx) ** 2 + (gp_dy - click_dy) ** 2
    return grid_points[int(np.argmin(dist2))]


# ---------------------------------------------------------------------------
# Contour canvas
# ---------------------------------------------------------------------------

class ContourCanvas(FigureCanvas):
    """Matplotlib canvas that draws the 2-D Pr_alpha contour map."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        fig = Figure(tight_layout=True)
        self.ax = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # State set after the first draw
        self._click_marker = None
        self._colorbar = None
        self._grid_points: list[GridPoint] = []
        self._result_map: dict[tuple[float, float], BlastPointResult] = {}
        # Callback invoked with the nearest GridPoint when user clicks
        self.on_point_selected = None  # type: ignore[assignment]

        self.mpl_connect("button_press_event", self._on_click)

    # ------------------------------------------------------------------
    # Public
    # ------------------------------------------------------------------

    def draw_contours(
        self,
        grid_points: list[GridPoint],
        result_map: dict[tuple[float, float], BlastPointResult],
        width: float,
        height: float,
        mach_curve: list[tuple[float, float]] | None = None,
    ) -> None:
        """Redraw the contour map with new data, with optional Mach stem zones."""
        self._grid_points = grid_points
        self._result_map = result_map
        self._click_marker = None

        import matplotlib.patches as mpatches

        X, Y, Pr = _build_meshgrid(grid_points, result_map)

        fig = self.figure
        fig.clear()
        self.ax = fig.add_subplot(111)
        ax = self.ax
        self._colorbar = None

        vmin, vmax = float(np.nanmin(Pr)), float(np.nanmax(Pr))

        has_mach = mach_curve is not None and len(mach_curve) > 0

        if has_mach:
            # Build a per-cell Mach height lookup
            mach_dict = {dx: mach_dy for dx, mach_dy in mach_curve}

            # Create masked arrays: above and below the Mach line
            Pr_above = Pr.copy()
            Pr_below = Pr.copy()

            dx_vals = sorted({gp.dx for gp in grid_points})
            dy_vals = sorted({gp.dy for gp in grid_points})

            for j, dy in enumerate(dy_vals):
                for i, dx in enumerate(dx_vals):
                    if dx in mach_dict:
                        if dy < mach_dict[dx]:
                            Pr_above[j, i] = np.nan  # mask from contourf
                        else:
                            Pr_below[j, i] = np.nan  # mask from pcolormesh

            # Upper zone: contourf (rings)
            cf = ax.contourf(X, Y, Pr_above, levels=15, cmap="YlOrRd",
                             vmin=vmin, vmax=vmax, zorder=1)
            cs = ax.contour(X, Y, Pr_above, levels=10, colors="black",
                            linewidths=0.8, zorder=2)
            ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f kPa")

            # Lower zone: pcolormesh (flat bands)
            ax.pcolormesh(X, Y, Pr_below, cmap="YlOrRd",
                          vmin=vmin, vmax=vmax, zorder=1, shading="nearest")

            # Mach curve
            mach_xs = [dx for dx, _ in mach_curve]
            mach_ys = [dy for _, dy in mach_curve]
            ax.plot(mach_xs, mach_ys, "k-", linewidth=2.5, zorder=3, label="Mach stem")

            # Colorbar from contourf
            self._colorbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
            self._colorbar.set_label("Pr_α (kPa)", fontsize=9)
        else:
            # No Mach stem: pure contourf (surface burst or Mach below facade)
            cf = ax.contourf(X, Y, Pr, levels=15, cmap="YlOrRd", zorder=1)
            cs = ax.contour(X, Y, Pr, levels=10, colors="black", linewidths=0.8, zorder=2)
            ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f kPa")
            self._colorbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
            self._colorbar.set_label("Pr_α (kPa)", fontsize=9)

        # Building outline
        half_w = width / 2.0
        half_h = height / 2.0
        ax.add_patch(mpatches.Rectangle(
            (-half_w, -half_h), width, height,
            linewidth=2, edgecolor="black", facecolor="none", zorder=4,
        ))

        # Perpendicular centre marker
        ax.plot(0.0, 0.0, "ko", markersize=6, zorder=5)

        # Clip view to building bounds with small margin
        margin = max(width, height) * 0.1
        ax.set_xlim(-half_w - margin, half_w + margin)
        ax.set_ylim(-half_h - margin, half_h + margin)

        ax.set_xlabel("Horizontal offset (m)", fontsize=9)
        ax.set_ylabel("Vertical offset (m)", fontsize=9)
        ax.set_title("Reflected peak overpressure Pr_α (kPa)", fontsize=10)
        ax.set_aspect("equal", adjustable="box")

        fig.tight_layout()
        self.draw()

    def highlight_point(self, dx: float, dy: float) -> None:
        """Mark the selected grid point with a cyan circle."""
        ax = self.ax
        if self._click_marker is not None:
            try:
                self._click_marker.remove()
            except ValueError:
                pass
        (self._click_marker,) = ax.plot(
            dx, dy, "o",
            color="cyan",
            markersize=10,
            markeredgecolor="black",
            markeredgewidth=1.5,
            zorder=5,
        )
        self.draw()

    # ------------------------------------------------------------------
    # Private
    # ------------------------------------------------------------------

    def _on_click(self, event) -> None:  # noqa: ANN001
        """Forward the nearest grid point to the registered callback."""
        if event.inaxes is not self.ax:
            return
        if not self._grid_points:
            return
        gp = _nearest_point(event.xdata, event.ydata, self._grid_points)
        self.highlight_point(gp.dx, gp.dy)
        if self.on_point_selected is not None:
            res = self._result_map[(gp.dx, gp.dy)]
            self.on_point_selected(gp, res)


# ---------------------------------------------------------------------------
# Friedlander panel
# ---------------------------------------------------------------------------

class FriedlanderPanel(QWidget):
    """Right-side panel: Friedlander P(t) plot with hover crosshair + blast parameter grid."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        # --- Waveform canvas ---
        fig = Figure(tight_layout=True)
        self._ax = fig.add_subplot(111)
        self._canvas = FigureCanvas(fig)
        self._canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(self._canvas, stretch=1)

        # Hover crosshair state
        self._vline = None
        self._hline = None
        self._hover_annotation = None
        self._t_data_ms = None  # time array in ms (for interpolation)
        self._p_data = None     # pressure array in kPa
        self._canvas.mpl_connect("motion_notify_event", self._on_hover)

        # --- Placeholder label ---
        self._label_placeholder = QLabel(
            "Click on the facade to see blast parameters.",
            alignment=Qt.AlignCenter,
        )
        self._label_placeholder.setWordWrap(True)
        layout.addWidget(self._label_placeholder)

        # --- Parameter grid (2 columns) ---
        self._param_labels: dict[str, QLabel] = {}
        # (key, display_label) — arranged in 2 columns
        param_defs = [
            # Column 0          Column 1
            ("coord",   None),  # spans full width
            ("R_alpha", "R_α"), ("alpha",   "α"),
            ("Ps0",     "Ps0"), ("C_alpha", "C_α"),
            ("Pr",      "Pr_α"),("ir",      "ir_α"),
            ("tA",      "tA"),  ("t0",      "t0"),
            ("b",       "b"),
        ]
        self._param_widget = QWidget()
        grid = QGridLayout(self._param_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(2)

        # First row: coord spans both columns
        lbl = QLabel()
        lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self._param_labels["coord"] = lbl
        grid.addWidget(lbl, 0, 0, 1, 2)

        # Remaining params in 2-column grid
        pair_keys = [
            ("R_alpha", "alpha"),
            ("Ps0",     "C_alpha"),
            ("Pr",      "ir"),
            ("tA",      "t0"),
            ("b",       None),
        ]
        for row_idx, (left_key, right_key) in enumerate(pair_keys, start=1):
            lbl_l = QLabel()
            lbl_l.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self._param_labels[left_key] = lbl_l
            grid.addWidget(lbl_l, row_idx, 0)
            if right_key:
                lbl_r = QLabel()
                lbl_r.setTextInteractionFlags(Qt.TextSelectableByMouse)
                self._param_labels[right_key] = lbl_r
                grid.addWidget(lbl_r, row_idx, 1)

        self._param_widget.hide()
        layout.addWidget(self._param_widget)

    # ------------------------------------------------------------------
    # Public
    # ------------------------------------------------------------------

    def update_point(self, gp: GridPoint, res: BlastPointResult, mach_zone: bool = False) -> None:
        """Redraw the waveform and update all parameter labels."""
        self._draw_waveform(gp, res)
        self._update_labels(gp, res, mach_zone=mach_zone)
        self._label_placeholder.hide()
        self._param_widget.show()

    # ------------------------------------------------------------------
    # Private
    # ------------------------------------------------------------------

    def _draw_waveform(self, gp: GridPoint, res: BlastPointResult) -> None:
        ax = self._ax
        ax.cla()

        t_start = max(0.0, res.tA - 5e-3)  # start 5 ms before arrival
        t_end = res.tA + 1.5 * res.t0
        t = np.linspace(t_start, t_end, 2000)
        P = friedlander(t, res.Pr_alpha, res.tA, res.t0, res.b)

        # Store for hover interpolation
        self._t_data_ms = t * 1e3
        self._p_data = P

        ax.plot(self._t_data_ms, P, color="steelblue", linewidth=1.8)
        ax.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
        ax.fill_between(self._t_data_ms, P, 0.0, where=(P > 0.0), alpha=0.2, color="steelblue")
        ax.set_xlabel("Time (ms)", fontsize=9)
        ax.set_ylabel("Reflected pressure (kPa)", fontsize=9)
        ax.set_title(
            f"Friedlander — dx={gp.dx:.1f} m, dy={gp.dy:.1f} m",
            fontsize=9,
        )
        ax.grid(True, alpha=0.25)

        # Reset hover elements
        self._vline = None
        self._hline = None
        self._hover_annotation = None

        self._canvas.draw()

    def _on_hover(self, event) -> None:
        """Show crosshair + values on hover over the Friedlander plot."""
        ax = self._ax
        if event.inaxes is not ax or self._t_data_ms is None:
            # Remove crosshair when outside axes
            if self._vline is not None:
                self._vline.set_visible(False)
                self._hline.set_visible(False)
                self._hover_annotation.set_visible(False)
                self._canvas.draw_idle()
            return

        t_ms = event.xdata
        # Interpolate pressure at cursor time
        p_val = float(np.interp(t_ms, self._t_data_ms, self._p_data))

        if self._vline is None:
            self._vline = ax.axvline(t_ms, color="gray", linewidth=0.8, linestyle="--")
            self._hline = ax.axhline(p_val, color="gray", linewidth=0.8, linestyle="--")
            # Fixed position in axes coordinates (top-left) to avoid plot resizing
            self._hover_annotation = ax.text(
                0.02, 0.97, "",
                transform=ax.transAxes,
                fontsize=8, verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray", alpha=0.9),
            )
        else:
            self._vline.set_xdata([t_ms])
            self._hline.set_ydata([p_val])
            self._vline.set_visible(True)
            self._hline.set_visible(True)
            self._hover_annotation.set_visible(True)

        self._hover_annotation.set_text(f"t = {t_ms:.2f} ms\nP = {p_val:.2f} kPa")

        self._canvas.draw_idle()

    def _update_labels(self, gp: GridPoint, res: BlastPointResult, mach_zone: bool = False) -> None:
        coord_text = f"Point:  dx = {gp.dx:.2f} m,  dy = {gp.dy:.2f} m"
        if mach_zone:
            coord_text += "  (Mach zone)"
        self._param_labels["coord"].setText(coord_text)
        self._param_labels["R_alpha"].setText(f"R_α = {gp.R_alpha:.3f} m")
        self._param_labels["alpha"].setText(f"α = {gp.alpha_deg:.2f}°")
        self._param_labels["Ps0"].setText(f"Ps0 = {res.Ps0:.2f} kPa")
        self._param_labels["C_alpha"].setText(f"C_α = {res.C_alpha:.4f}")
        self._param_labels["Pr"].setText(f"Pr_α = {res.Pr_alpha:.2f} kPa")
        self._param_labels["ir"].setText(f"ir_α = {res.ir_alpha*1e3:.2f} kPa·ms")
        self._param_labels["tA"].setText(f"tA = {res.tA*1e3:.3f} ms")
        self._param_labels["t0"].setText(f"t0 = {res.t0*1e3:.3f} ms")
        self._param_labels["b"].setText(f"b = {res.b:.4f}")


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class BlastWindow(QMainWindow):
    """Main application window."""

    def __init__(
        self,
        W: float,
        R: float,
        Hc: float,
        width: float,
        height: float,
        step: float,
    ) -> None:
        super().__init__()
        self.setWindowTitle("UFC 3-340-02  Blast Facade Calculator")
        self.resize(1300, 750)

        # Ensure UFC tables are loaded once
        load_ufc_tables()

        # ---- Central widget ----
        central = QWidget()
        self.setCentralWidget(central)
        root_layout = QVBoxLayout(central)
        root_layout.setContentsMargins(6, 6, 6, 6)
        root_layout.setSpacing(4)

        # ---- Top input bar ----
        root_layout.addWidget(self._build_input_bar(W, R, Hc, width, height, step))

        # ---- Splitter: contour | Friedlander panel ----
        splitter = QSplitter(Qt.Horizontal)
        splitter.setHandleWidth(4)

        self._contour_canvas = ContourCanvas()
        self._contour_canvas.on_point_selected = self._on_point_selected
        splitter.addWidget(self._contour_canvas)

        self._friedlander_panel = FriedlanderPanel()
        splitter.addWidget(self._friedlander_panel)

        # 65 / 35 split
        splitter.setStretchFactor(0, 65)
        splitter.setStretchFactor(1, 35)
        root_layout.addWidget(splitter, stretch=1)

        # ---- Status bar ----
        self._status = QStatusBar()
        self.setStatusBar(self._status)

        # ---- Initial computation ----
        self._recalculate()

    # ------------------------------------------------------------------
    # Input bar
    # ------------------------------------------------------------------

    def _build_input_bar(
        self,
        W: float,
        R: float,
        Hc: float,
        width: float,
        height: float,
        step: float,
    ) -> QWidget:
        bar = QWidget()
        layout = QHBoxLayout(bar)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        def _spinbox(
            label: str,
            lo: float,
            hi: float,
            val: float,
            decimals: int = 1,
        ) -> QDoubleSpinBox:
            lbl = QLabel(label)
            sb = QDoubleSpinBox()
            sb.setRange(lo, hi)
            sb.setDecimals(decimals)
            sb.setValue(val)
            sb.setFixedWidth(90)
            # Enter key triggers recalculation
            sb.editingFinished.connect(self._recalculate)
            layout.addWidget(lbl)
            layout.addWidget(sb)
            return sb

        self._spin_W      = _spinbox("W (kg):",    1.0,    100_000.0,  W,     1)
        self._spin_R      = _spinbox("R (m):",     1.0,    1_000.0,    R,     1)
        self._spin_Hc     = _spinbox("Hc (m):",    0.0,    100.0,      Hc,    2)
        self._spin_width  = _spinbox("Width (m):", 1.0,    200.0,      width, 1)
        self._spin_height = _spinbox("Height (m):",1.0,    200.0,      height,1)
        self._spin_step   = _spinbox("Step (m):",  0.1,    10.0,       step,  2)

        self._recalc_btn = QPushButton("Recalculate")
        self._recalc_btn.setFixedWidth(110)
        self._recalc_btn.clicked.connect(self._recalculate)
        layout.addWidget(self._recalc_btn)

        # Burst type indicator
        self._burst_type_label = QLabel("")
        self._burst_type_label.setStyleSheet("font-weight: bold; padding-left: 10px;")
        layout.addWidget(self._burst_type_label)

        layout.addStretch()
        return bar

    # ------------------------------------------------------------------
    # Calculation
    # ------------------------------------------------------------------

    def _recalculate(self) -> None:
        """Launch blast computation in a background thread."""
        W      = self._spin_W.value()
        R      = self._spin_R.value()
        Hc     = self._spin_Hc.value()
        width  = self._spin_width.value()
        height = self._spin_height.value()
        step   = self._spin_step.value()

        # Estimate point count for status feedback
        if width and height:
            nx = int(width / step) + 1
            ny = int(height / step) + 1
            est = nx * ny
        else:
            est = "?"

        self._status.showMessage(f"Computing ~{est} points...")
        # Disable recalculate button while working
        self._recalc_btn.setEnabled(False)

        self._worker = _ComputeWorker(R, W, Hc, width, height, step)
        self._worker.finished.connect(self._on_compute_done)
        self._worker_width = width
        self._worker_height = height
        self._worker.start()

    def _on_compute_done(self) -> None:
        """Called when background computation finishes."""
        self._recalc_btn.setEnabled(True)

        if self._worker.error:
            self._status.showMessage(f"Error: {self._worker.error}")
            return

        self._grid_points = self._worker.grid_points
        self._result_map  = self._worker.result_map
        self._mach_curve  = self._worker.mach_curve
        self._width       = self._worker_width
        self._height      = self._worker_height

        # Update burst type label
        bt = self._worker.burst_type
        hob = self._worker.scaled_hob
        if bt == "air":
            self._burst_type_label.setText(f"Air burst (Hc/W^\u215b = {hob:.2f})")
            self._burst_type_label.setStyleSheet("font-weight: bold; padding-left: 10px; color: #2060c0;")
        else:
            self._burst_type_label.setText(f"Surface burst (Hc/W^\u215b = {hob:.2f})")
            self._burst_type_label.setStyleSheet("font-weight: bold; padding-left: 10px; color: #c06020;")

        self._contour_canvas.draw_contours(
            self._grid_points, self._result_map, self._width, self._height,
            mach_curve=self._mach_curve,
        )
        self._status.showMessage(
            f"Calculated {len(self._grid_points)} points in {self._worker.elapsed:.2f} s"
        )

    # ------------------------------------------------------------------
    # Click callback
    # ------------------------------------------------------------------

    def _on_point_selected(
        self,
        gp: GridPoint,
        res: BlastPointResult,
    ) -> None:
        # Check if point is in Mach zone
        mach_zone = False
        if hasattr(self, '_mach_curve') and self._mach_curve:
            mach_dict = {dx: mach_dy for dx, mach_dy in self._mach_curve}
            if gp.dx in mach_dict and gp.dy < mach_dict[gp.dx]:
                mach_zone = True
        self._friedlander_panel.update_point(gp, res, mach_zone=mach_zone)
        self._status.showMessage(
            f"Selected point  dx={gp.dx:.1f} m  dy={gp.dy:.1f} m  "
            f"Pr_α={res.Pr_alpha:.1f} kPa"
            + ("  [Mach zone]" if mach_zone else "")
        )


# ---------------------------------------------------------------------------
# Public API — matches the CLI call signature
# ---------------------------------------------------------------------------

def launch_gui(
    W: float = 200.0,
    R: float = 30.0,
    Hc: float = 5.0,
    width: float = 10.0,
    height: float = 8.0,
    step: float = 1.0,
) -> None:
    """Launch the interactive PySide6 blast facade visualisation.

    Parameters
    ----------
    W : float
        Charge mass (kg TNT equivalent).
    R : float
        Perpendicular standoff distance from charge to facade (m).
    Hc : float
        Height of burst above ground (m).
    width : float
        Facade width (m).
    height : float
        Facade height (m).
    step : float
        Grid spacing (m).
    """
    app = QApplication.instance() or QApplication([])
    window = BlastWindow(
        W=W, R=R, Hc=Hc,
        width=width, height=height,
        step=step,
    )
    window.show()
    app.exec()
