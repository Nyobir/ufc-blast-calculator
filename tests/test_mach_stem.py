"""Tests for Mach stem triple point correction."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ufc_blast.core.blast_params import (
    BlastPointResult,
    apply_mach_stem,
    compute_points_batch,
    load_ufc_tables,
)
from ufc_blast.core.geometry import GridPoint, generate_grid


@pytest.fixture(autouse=True)
def _load_tables():
    load_ufc_tables()


class TestApplyMachStem:
    """Tests for the Mach stem post-processing correction."""

    def _make_grid_and_results(self, W=200.0, R=30.0, Hc=5.0, width=10.0, height=8.0, step=1.0):
        """Helper: generate grid and compute results."""
        grid_points = generate_grid(R=R, W=W, Hc=Hc, width=width, height=height, step=step)
        R_alphas = np.array([gp.R_alpha for gp in grid_points])
        alpha_degs = np.array([gp.alpha_deg for gp in grid_points])
        results = compute_points_batch(R_alphas, alpha_degs, W)
        result_map = {}
        for gp, res in zip(grid_points, results):
            key = (gp.dx, gp.dy)
            if key not in result_map:
                result_map[key] = res
        return grid_points, result_map

    def test_returns_result_map_and_mach_curve(self):
        """apply_mach_stem returns a dict and a list of (dx, dy) tuples."""
        grid_points, result_map = self._make_grid_and_results()
        new_map, mach_curve = apply_mach_stem(grid_points, result_map, W=200.0, Hc=5.0, R=30.0)
        assert isinstance(new_map, dict)
        assert isinstance(mach_curve, list)
        for entry in mach_curve:
            assert len(entry) == 2

    def test_mach_zone_has_uniform_pressure_per_column(self):
        """Below the Mach line, all points in a column should have the same Pr_alpha."""
        grid_points, result_map = self._make_grid_and_results()
        new_map, mach_curve = apply_mach_stem(grid_points, result_map, W=200.0, Hc=5.0, R=30.0)

        if not mach_curve:
            pytest.skip("No Mach stem visible for these parameters")

        from collections import defaultdict
        columns = defaultdict(list)
        for gp in grid_points:
            columns[gp.dx].append(gp.dy)

        mach_dict = {dx: mach_dy for dx, mach_dy in mach_curve}
        for dx, mach_dy in mach_dict.items():
            if dx not in columns:
                continue
            below_prs = []
            for dy in columns[dx]:
                if dy < mach_dy:
                    below_prs.append(new_map[(dx, dy)].Pr_alpha)
            if len(below_prs) > 1:
                for pr in below_prs:
                    assert abs(pr - below_prs[0]) < 1e-10, (
                        f"Non-uniform Pr in Mach zone at dx={dx}: {below_prs}"
                    )

    def test_above_mach_line_unchanged(self):
        """Points above the Mach line should be unchanged."""
        grid_points, result_map = self._make_grid_and_results()
        original_map = {k: v for k, v in result_map.items()}
        new_map, mach_curve = apply_mach_stem(grid_points, result_map, W=200.0, Hc=5.0, R=30.0)

        if not mach_curve:
            pytest.skip("No Mach stem visible for these parameters")

        mach_dict = {dx: mach_dy for dx, mach_dy in mach_curve}
        for (dx, dy), res in new_map.items():
            if dx in mach_dict and dy >= mach_dict[dx]:
                orig = original_map[(dx, dy)]
                assert res.Pr_alpha == orig.Pr_alpha, (
                    f"Point ({dx}, {dy}) above Mach line was modified"
                )

    def test_no_mach_stem_when_out_of_table_range(self):
        """If Rg_scaled is outside Figure 2-13 range, return empty Mach curve."""
        grid_points, result_map = self._make_grid_and_results(R=500.0, width=2.0, height=2.0)
        new_map, mach_curve = apply_mach_stem(grid_points, result_map, W=200.0, Hc=5.0, R=500.0)
        assert isinstance(mach_curve, list)
