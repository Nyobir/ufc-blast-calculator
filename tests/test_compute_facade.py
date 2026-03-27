"""Tests for the shared compute_facade() pipeline.

Both the CLI and GUI call compute_facade() — these tests verify the
single computation path that both interfaces rely on.
"""

from __future__ import annotations

import pytest

from ufc_blast.core.blast_params import (
    BlastPointResult,
    FacadeResult,
    compute_facade,
    load_ufc_tables,
)


@pytest.fixture(autouse=True)
def _load_tables():
    load_ufc_tables()


# ---- Air burst scenario ----
AIR_PARAMS = dict(R=30.0, W=200.0, Hc=5.0, width=10.0, height=8.0, step=2.0)

# ---- Surface burst scenario ----
SURFACE_PARAMS = dict(R=30.0, W=1000.0, Hc=3.0, width=10.0, height=8.0, step=2.0)


class TestComputeFacadeAirBurst:
    """Air burst pipeline: grid → batch compute → Mach stem."""

    def test_returns_facade_result(self):
        facade = compute_facade(**AIR_PARAMS)
        assert isinstance(facade, FacadeResult)
        assert facade.burst_type == "air"
        assert facade.scaled_hob > 0.397

    def test_grid_points_populated(self):
        facade = compute_facade(**AIR_PARAMS)
        assert len(facade.grid_points) > 0

    def test_result_map_covers_all_grid_points(self):
        facade = compute_facade(**AIR_PARAMS)
        for gp in facade.grid_points:
            key = (gp.dx, gp.dy)
            assert key in facade.result_map, f"Missing result for {key}"
            assert isinstance(facade.result_map[key], BlastPointResult)

    def test_mach_curve_present(self):
        """Air burst should produce a Mach curve."""
        facade = compute_facade(**AIR_PARAMS)
        assert isinstance(facade.mach_curve, list)
        # Mach curve entries are (dx, mach_dy) tuples
        for entry in facade.mach_curve:
            assert len(entry) == 2

    def test_mach_zone_uniform_pressure(self):
        """Below the Mach line, all points in a column must have identical Pr."""
        facade = compute_facade(**AIR_PARAMS)
        if not facade.mach_curve:
            pytest.skip("No Mach stem visible for these parameters")

        from collections import defaultdict

        columns = defaultdict(list)
        for gp in facade.grid_points:
            columns[gp.dx].append(gp.dy)

        mach_dict = {dx: mach_dy for dx, mach_dy in facade.mach_curve}
        for dx, mach_dy in mach_dict.items():
            below_prs = [
                facade.result_map[(dx, dy)].Pr_alpha
                for dy in columns[dx]
                if dy < mach_dy
            ]
            if len(below_prs) > 1:
                assert all(
                    abs(pr - below_prs[0]) < 1e-10 for pr in below_prs
                ), f"Non-uniform Pr in Mach zone at dx={dx}"

    def test_deterministic(self):
        """Two calls with same params must produce bit-identical results."""
        f1 = compute_facade(**AIR_PARAMS)
        f2 = compute_facade(**AIR_PARAMS)

        assert f1.burst_type == f2.burst_type
        assert f1.scaled_hob == f2.scaled_hob
        assert len(f1.grid_points) == len(f2.grid_points)
        assert len(f1.mach_curve) == len(f2.mach_curve)

        for gp in f1.grid_points:
            key = (gp.dx, gp.dy)
            r1 = f1.result_map[key]
            r2 = f2.result_map[key]
            assert r1.Ps0 == r2.Ps0, f"Ps0 mismatch at {key}"
            assert r1.Pr_alpha == r2.Pr_alpha, f"Pr_alpha mismatch at {key}"
            assert r1.ir_alpha == r2.ir_alpha, f"ir_alpha mismatch at {key}"
            assert r1.tA == r2.tA, f"tA mismatch at {key}"
            assert r1.t0 == r2.t0, f"t0 mismatch at {key}"
            assert r1.b == r2.b, f"b mismatch at {key}"


class TestComputeFacadeSurfaceBurst:
    """Surface burst pipeline: grid → batch compute (no Mach stem)."""

    def test_returns_surface_burst(self):
        facade = compute_facade(**SURFACE_PARAMS)
        assert facade.burst_type == "surface"
        assert facade.scaled_hob <= 0.397

    def test_no_mach_curve(self):
        """Surface bursts should not have a Mach curve."""
        facade = compute_facade(**SURFACE_PARAMS)
        assert facade.mach_curve == []

    def test_result_map_covers_all_grid_points(self):
        facade = compute_facade(**SURFACE_PARAMS)
        for gp in facade.grid_points:
            key = (gp.dx, gp.dy)
            assert key in facade.result_map
            res = facade.result_map[key]
            assert res.Ps0 > 0.0
            assert res.Pr_alpha > 0.0

    def test_surface_differs_from_air(self):
        """Surface and air bursts at same geometry give different Ps0."""
        # Use same R and W but different Hc to get different burst types
        facade_air = compute_facade(R=30.0, W=200.0, Hc=5.0, width=4.0, height=4.0, step=2.0)
        facade_surface = compute_facade(R=30.0, W=1000.0, Hc=3.0, width=4.0, height=4.0, step=2.0)

        # Centre point (0, 0) — both grids have it
        air_pr = facade_air.result_map[(0.0, 0.0)].Ps0
        surf_pr = facade_surface.result_map[(0.0, 0.0)].Ps0
        assert air_pr != surf_pr, "Air and surface should use different tables"

    def test_deterministic(self):
        """Two calls with same params must produce bit-identical results."""
        f1 = compute_facade(**SURFACE_PARAMS)
        f2 = compute_facade(**SURFACE_PARAMS)

        for gp in f1.grid_points:
            key = (gp.dx, gp.dy)
            r1 = f1.result_map[key]
            r2 = f2.result_map[key]
            assert r1.Ps0 == r2.Ps0
            assert r1.Pr_alpha == r2.Pr_alpha
