"""
Headless tests for the graphical interface.

Covers:
  - the ``gui-check`` headless smoke test (imports the GUI dependencies,
    loads the UFC tables, runs one facade) exits cleanly
  - the GUI background worker returns exactly what ``compute_facade``
    returns — the GUI holds no calculation logic of its own
  - the GUI worker and the ``compute`` CLI subcommand produce identical
    numbers for identical inputs (single shared computation path)

PySide6 has no wheel on every platform this package is analysed on, so the
whole module is skipped when it cannot be imported.
"""

from __future__ import annotations

import argparse
import csv
import os

import pytest

# Must be set before any Qt import: no display is available under pytest.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6", reason="PySide6 unavailable on this platform")

from ufc_blast.cli import _cmd_compute, _cmd_gui_check  # noqa: E402
from ufc_blast.core.blast_params import compute_facade  # noqa: E402
from ufc_blast.gui import _ComputeWorker  # noqa: E402


# Air burst, small enough to stay fast, large enough to exercise the batch
# lookup path and the per-column Mach correction.
_CASE = {"W": 200.0, "R": 30.0, "Hc": 5.0, "width": 6.0, "height": 6.0, "step": 1.0}


@pytest.fixture(scope="module")
def qapp():
    """A single offscreen QApplication for the module."""
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture(scope="module")
def worker_facade(qapp):
    """Facade result produced through the GUI's own worker object."""
    worker = _ComputeWorker(
        R=_CASE["R"], W=_CASE["W"], Hc=_CASE["Hc"],
        width=_CASE["width"], height=_CASE["height"], step=_CASE["step"],
    )
    worker.run()
    assert worker.error is None, f"GUI worker failed: {worker.error}"
    assert worker.facade is not None
    return worker.facade


class TestGuiCheck:
    """The ``ufc-blast gui-check`` headless smoke test."""

    def test_gui_check_exits_clean(self, capsys):
        rc = _cmd_gui_check(argparse.Namespace())
        out = capsys.readouterr().out
        assert rc == 0
        assert "GUI check OK" in out


class TestGuiSharedComputationPath:
    """The GUI must not compute anything the core does not compute."""

    def test_worker_matches_compute_facade(self, worker_facade):
        core = compute_facade(
            R=_CASE["R"], W=_CASE["W"], Hc=_CASE["Hc"],
            width=_CASE["width"], height=_CASE["height"], step=_CASE["step"],
        )

        assert worker_facade.burst_type == core.burst_type == "air"
        assert worker_facade.scaled_hob == pytest.approx(core.scaled_hob)
        assert len(worker_facade.grid_points) == len(core.grid_points)
        assert worker_facade.mach_curve == core.mach_curve
        assert worker_facade.mach_curve, "air burst must produce a Mach curve"

        for gp in core.grid_points:
            key = (gp.dx, gp.dy)
            gui_res = worker_facade.result_map[key]
            core_res = core.result_map[key]
            assert gui_res.Ps0 == pytest.approx(core_res.Ps0)
            assert gui_res.C_alpha == pytest.approx(core_res.C_alpha)
            assert gui_res.Pr_alpha == pytest.approx(core_res.Pr_alpha)
            assert gui_res.ir_alpha == pytest.approx(core_res.ir_alpha)
            assert gui_res.tA == pytest.approx(core_res.tA)
            assert gui_res.t0 == pytest.approx(core_res.t0)
            assert gui_res.b == pytest.approx(core_res.b)
            assert gui_res.extrapolated == core_res.extrapolated

    def test_worker_matches_cli_compute_csv(self, worker_facade, tmp_path, capsys):
        """The ``compute`` CLI subcommand and the GUI agree point by point."""
        out_csv = tmp_path / "facade.csv"
        args = argparse.Namespace(
            W=_CASE["W"], R=_CASE["R"], Hc=_CASE["Hc"],
            width=_CASE["width"], height=_CASE["height"], step=_CASE["step"],
            output=str(out_csv), format="csv",
        )
        rc = _cmd_compute(args)
        capsys.readouterr()
        assert rc == 0

        with open(out_csv, newline="") as fh:
            rows = list(csv.DictReader(fh))

        assert len(rows) == len(worker_facade.grid_points)

        for row in rows:
            key = (float(row["dx_m"]), float(row["dy_m"]))
            res = worker_facade.result_map[key]
            # The CLI rounds to 4 decimals on write; compare at that precision.
            assert float(row["Pr_alpha_kPa"]) == pytest.approx(res.Pr_alpha, abs=1e-4)
            assert float(row["Ps0_kPa"]) == pytest.approx(res.Ps0, abs=1e-4)
            assert float(row["ir_alpha_kPa_ms"]) == pytest.approx(
                res.ir_alpha * 1000.0, abs=1e-4)
            assert float(row["tA_ms"]) == pytest.approx(res.tA * 1000.0, abs=1e-4)
            assert float(row["t0_ms"]) == pytest.approx(res.t0 * 1000.0, abs=1e-4)
            assert row["extrapolated"] == ",".join(sorted(res.extrapolated))
