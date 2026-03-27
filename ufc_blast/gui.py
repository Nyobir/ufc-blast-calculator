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

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QDoubleSpinBox,
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
    compute_point,
    friedlander,
    load_ufc_tables,
)
from ufc_blast.core.geometry import GridPoint, generate_grid


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def _compute_results(
    grid_points: list[GridPoint],
    W: float,
) -> dict[tuple[float, float], BlastPointResult]:
    """Compute blast parameters for every (dx, dy) grid point.

    Returns a mapping  (dx, dy) → BlastPointResult.
    """
    result_map: dict[tuple[float, float], BlastPointResult] = {}
    for gp in grid_points:
        key = (gp.dx, gp.dy)
        if key not in result_map:
            result_map[key] = compute_point(gp.R_alpha, gp.alpha_deg, W)
    return result_map


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
    ) -> None:
        """Redraw the contour map with new data."""
        self._grid_points = grid_points
        self._result_map = result_map
        self._click_marker = None

        X, Y, Pr = _build_meshgrid(grid_points, result_map)

        ax = self.ax
        ax.cla()

        # Filled contours
        cf = ax.contourf(X, Y, Pr, levels=15, cmap="YlOrRd", zorder=1)

        # Contour lines
        cs = ax.contour(X, Y, Pr, levels=10, colors="black", linewidths=0.8, zorder=2)
        ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f kPa")

        # Colorbar
        cbar = ax.get_figure().colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Pr_α (kPa)", fontsize=9)

        # Building outline
        half_w = width / 2.0
        half_h = height / 2.0
        rect = ax.add_patch(
            __import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle(
                (-half_w, -half_h),
                width,
                height,
                linewidth=2,
                edgecolor="black",
                facecolor="none",
                zorder=3,
            )
        )

        # Perpendicular centre marker
        ax.plot(0.0, 0.0, "ko", markersize=6, zorder=4, label="Perpendicular pt.")

        ax.set_xlabel("Horizontal offset (m)", fontsize=9)
        ax.set_ylabel("Vertical offset (m)", fontsize=9)
        ax.set_title("Reflected peak overpressure Pr_α (kPa)", fontsize=10)
        ax.set_aspect("equal", adjustable="box")

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
    """Right-side panel: Friedlander P(t) plot + blast parameter labels."""

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
        layout.addWidget(self._canvas)

        # --- Parameter labels ---
        self._label_placeholder = QLabel(
            "Click on the facade to see blast parameters.",
            alignment=Qt.AlignCenter,
        )
        self._label_placeholder.setWordWrap(True)
        layout.addWidget(self._label_placeholder)

        # Individual parameter label widgets (hidden until first click)
        self._param_labels: dict[str, QLabel] = {}
        param_keys = [
            ("coord",   "Point"),
            ("R_alpha", "R_α"),
            ("alpha",   "α"),
            ("Ps0",     "Ps0"),
            ("C_alpha", "C_α"),
            ("Pr",      "Pr_α"),
            ("ir",      "ir_α"),
            ("tA",      "tA"),
            ("t0",      "t0"),
            ("b",       "b"),
        ]
        self._param_widget = QWidget()
        param_layout = QVBoxLayout(self._param_widget)
        param_layout.setContentsMargins(0, 0, 0, 0)
        param_layout.setSpacing(2)
        for key, _ in param_keys:
            lbl = QLabel()
            lbl.setWordWrap(False)
            lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self._param_labels[key] = lbl
            param_layout.addWidget(lbl)
        self._param_widget.hide()
        layout.addWidget(self._param_widget)

    # ------------------------------------------------------------------
    # Public
    # ------------------------------------------------------------------

    def update_point(self, gp: GridPoint, res: BlastPointResult) -> None:
        """Redraw the waveform and update all parameter labels."""
        self._draw_waveform(gp, res)
        self._update_labels(gp, res)

        self._label_placeholder.hide()
        self._param_widget.show()

    # ------------------------------------------------------------------
    # Private
    # ------------------------------------------------------------------

    def _draw_waveform(self, gp: GridPoint, res: BlastPointResult) -> None:
        ax = self._ax
        ax.cla()

        t_end = res.tA + 1.5 * res.t0
        t = np.linspace(0.0, t_end, 2000)
        P = friedlander(t, res.Pr_alpha, res.tA, res.t0, res.b)

        ax.plot(t * 1e3, P, color="steelblue", linewidth=1.8)
        ax.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
        ax.fill_between(t * 1e3, P, 0.0, where=(P > 0.0), alpha=0.2, color="steelblue")
        ax.set_xlabel("Time (ms)", fontsize=9)
        ax.set_ylabel("Reflected pressure (kPa)", fontsize=9)
        ax.set_title(
            f"Friedlander — dx={gp.dx:.1f} m, dy={gp.dy:.1f} m",
            fontsize=9,
        )
        ax.grid(True, alpha=0.25)
        self._canvas.draw()

    def _update_labels(self, gp: GridPoint, res: BlastPointResult) -> None:
        self._param_labels["coord"].setText(
            f"Point:   dx = {gp.dx:.2f} m,  dy = {gp.dy:.2f} m"
        )
        self._param_labels["R_alpha"].setText(f"R_α  =  {gp.R_alpha:.3f} m")
        self._param_labels["alpha"].setText(f"α     =  {gp.alpha_deg:.2f} °")
        self._param_labels["Ps0"].setText(f"Ps0  =  {res.Ps0:.3f} kPa")
        self._param_labels["C_alpha"].setText(f"C_α  =  {res.C_alpha:.4f}")
        self._param_labels["Pr"].setText(f"Pr_α =  {res.Pr_alpha:.3f} kPa")
        self._param_labels["ir"].setText(f"ir_α =  {res.ir_alpha:.6f} kPa·s")
        self._param_labels["tA"].setText(f"tA   =  {res.tA*1e3:.4f} ms")
        self._param_labels["t0"].setText(f"t0   =  {res.t0*1e3:.4f} ms")
        self._param_labels["b"].setText(f"b    =  {res.b:.4f}")


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
        self._spin_Hc     = _spinbox("Hc (m):",    0.1,    100.0,      Hc,    2)
        self._spin_width  = _spinbox("Width (m):", 1.0,    200.0,      width, 1)
        self._spin_height = _spinbox("Height (m):",1.0,    200.0,      height,1)
        self._spin_step   = _spinbox("Step (m):",  0.1,    10.0,       step,  2)

        recalc_btn = QPushButton("Recalculate")
        recalc_btn.setFixedWidth(110)
        recalc_btn.clicked.connect(self._recalculate)
        layout.addWidget(recalc_btn)

        layout.addStretch()
        return bar

    # ------------------------------------------------------------------
    # Calculation
    # ------------------------------------------------------------------

    def _recalculate(self) -> None:
        """Re-run the blast grid computation and refresh the contour plot."""
        W      = self._spin_W.value()
        R      = self._spin_R.value()
        Hc     = self._spin_Hc.value()
        width  = self._spin_width.value()
        height = self._spin_height.value()
        step   = self._spin_step.value()

        t0 = time.perf_counter()
        try:
            grid_points = generate_grid(
                R, W, Hc,
                width=width,
                height=height,
                step=step,
            )
        except ValueError as exc:
            self._status.showMessage(f"Error: {exc}")
            return

        result_map = _compute_results(grid_points, W)
        elapsed = time.perf_counter() - t0

        self._grid_points = grid_points
        self._result_map  = result_map
        self._width       = width
        self._height      = height

        self._contour_canvas.draw_contours(
            grid_points, result_map, width, height
        )
        self._status.showMessage(
            f"Calculated {len(grid_points)} points in {elapsed:.2f} s"
        )

    # ------------------------------------------------------------------
    # Click callback
    # ------------------------------------------------------------------

    def _on_point_selected(
        self,
        gp: GridPoint,
        res: BlastPointResult,
    ) -> None:
        self._friedlander_panel.update_point(gp, res)
        self._status.showMessage(
            f"Selected point  dx={gp.dx:.1f} m  dy={gp.dy:.1f} m  "
            f"Pr_α={res.Pr_alpha:.1f} kPa"
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
