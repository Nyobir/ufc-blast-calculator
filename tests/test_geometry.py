"""
Tests for ufc_blast.core.geometry module.

Covers:
- Perpendicular point geometry
- 3-4-5 triangle (classic Pythagorean)
- Diagonal offset (R=10, dx=6, dy=8 → 45° angle)
- Symmetry of positive/negative offsets
- Building grid dimensions and centre point
- Surface-burst applicability check
- Open field mode basic sanity
"""

import math

import pytest

from ufc_blast.core.geometry import GridPoint, compute_point_geometry, generate_grid

# ---------------------------------------------------------------------------
# Tolerance for floating-point comparisons
# ---------------------------------------------------------------------------
ABS_TOL = 1e-9
REL_TOL = 1e-9

# Charge parameters that satisfy the applicability check:
#   Hc / W^(1/3) = 10 / 1^(1/3) = 10 > 0.397
VALID_W = 1.0   # kg
VALID_Hc = 10.0  # m


class TestComputePointGeometry:
    """Unit tests for compute_point_geometry."""

    def test_perpendicular_point(self):
        """dx=0, dy=0 should give R_alpha=R and alpha=0°."""
        R = 15.0
        gp = compute_point_geometry(dx=0.0, dy=0.0, R=R)
        assert isinstance(gp, GridPoint)
        assert math.isclose(gp.R_alpha, R, rel_tol=REL_TOL)
        assert math.isclose(gp.alpha_deg, 0.0, abs_tol=ABS_TOL)
        assert gp.dx == 0.0
        assert gp.dy == 0.0

    def test_345_triangle(self):
        """R=4, dx=3 → R_alpha=5, alpha=arccos(4/5)."""
        R = 4.0
        dx = 3.0
        gp = compute_point_geometry(dx=dx, dy=0.0, R=R)
        expected_R_alpha = 5.0  # sqrt(16 + 9)
        expected_alpha = math.degrees(math.acos(4.0 / 5.0))
        assert math.isclose(gp.R_alpha, expected_R_alpha, rel_tol=REL_TOL)
        assert math.isclose(gp.alpha_deg, expected_alpha, rel_tol=REL_TOL)

    def test_diagonal_offset_45_degrees(self):
        """R=10, dx=6, dy=8 → R_alpha=sqrt(200), alpha=45°."""
        R = 10.0
        dx = 6.0
        dy = 8.0
        gp = compute_point_geometry(dx=dx, dy=dy, R=R)
        expected_R_alpha = math.sqrt(200.0)  # sqrt(100 + 36 + 64)
        # cos(alpha) = R / R_alpha = 10 / sqrt(200) = 1/sqrt(2)  →  alpha = 45°
        expected_alpha = 45.0
        assert math.isclose(gp.R_alpha, expected_R_alpha, rel_tol=REL_TOL)
        assert math.isclose(gp.alpha_deg, expected_alpha, abs_tol=1e-10)

    def test_symmetry_horizontal(self):
        """Positive and negative dx give the same R_alpha and alpha."""
        R = 20.0
        dx = 7.5
        gp_pos = compute_point_geometry(dx=+dx, dy=0.0, R=R)
        gp_neg = compute_point_geometry(dx=-dx, dy=0.0, R=R)
        assert math.isclose(gp_pos.R_alpha, gp_neg.R_alpha, rel_tol=REL_TOL)
        assert math.isclose(gp_pos.alpha_deg, gp_neg.alpha_deg, rel_tol=REL_TOL)

    def test_symmetry_vertical(self):
        """Positive and negative dy give the same R_alpha and alpha."""
        R = 20.0
        dy = 5.0
        gp_pos = compute_point_geometry(dx=0.0, dy=+dy, R=R)
        gp_neg = compute_point_geometry(dx=0.0, dy=-dy, R=R)
        assert math.isclose(gp_pos.R_alpha, gp_neg.R_alpha, rel_tol=REL_TOL)
        assert math.isclose(gp_pos.alpha_deg, gp_neg.alpha_deg, rel_tol=REL_TOL)

    def test_symmetry_diagonal(self):
        """Diagonal offsets with mirrored signs give same R_alpha and alpha."""
        R = 30.0
        dx, dy = 4.0, 9.0
        gp_pp = compute_point_geometry(dx=+dx, dy=+dy, R=R)
        gp_pn = compute_point_geometry(dx=+dx, dy=-dy, R=R)
        gp_np = compute_point_geometry(dx=-dx, dy=+dy, R=R)
        gp_nn = compute_point_geometry(dx=-dx, dy=-dy, R=R)
        for gp in (gp_pn, gp_np, gp_nn):
            assert math.isclose(gp.R_alpha, gp_pp.R_alpha, rel_tol=REL_TOL)
            assert math.isclose(gp.alpha_deg, gp_pp.alpha_deg, rel_tol=REL_TOL)


class TestGenerateGridBuildingMode:
    """Tests for generate_grid in building mode (width + height specified)."""

    def test_building_grid_point_count(self):
        """width=10, height=8, step=1 → 11×9 = 99 points."""
        points = generate_grid(
            R=30.0, W=VALID_W, Hc=VALID_Hc,
            width=10.0, height=8.0, step=1.0,
        )
        assert len(points) == 99

    def test_building_grid_centre_point(self):
        """Centre point of the building grid must have dx=0, dy=0."""
        points = generate_grid(
            R=30.0, W=VALID_W, Hc=VALID_Hc,
            width=10.0, height=8.0, step=1.0,
        )
        centre_points = [p for p in points if p.dx == 0.0 and p.dy == 0.0]
        assert len(centre_points) == 1, "Exactly one centre point expected"
        cp = centre_points[0]
        assert math.isclose(cp.alpha_deg, 0.0, abs_tol=ABS_TOL)

    def test_building_grid_geometry_correctness(self):
        """Every point in the building grid satisfies the geometry formulas."""
        R = 25.0
        points = generate_grid(
            R=R, W=VALID_W, Hc=VALID_Hc,
            width=6.0, height=4.0, step=1.0,
        )
        for gp in points:
            expected_R_alpha = math.sqrt(R**2 + gp.dx**2 + gp.dy**2)
            assert math.isclose(gp.R_alpha, expected_R_alpha, rel_tol=1e-12)
            expected_alpha = math.degrees(math.acos(R / expected_R_alpha))
            assert math.isclose(gp.alpha_deg, expected_alpha, abs_tol=1e-10)


class TestGenerateGridApplicability:
    """Tests for the UFC 3-340-02 surface-burst applicability guard."""

    def test_surface_burst_raises(self):
        """Hc/W^(1/3) <= 0.397 must raise ValueError."""
        # W=1000 kg → W^(1/3)≈10, Hc=3 → ratio=0.3 < 0.397
        with pytest.raises(ValueError, match="0.397"):
            generate_grid(R=30.0, W=1000.0, Hc=3.0)

    def test_exact_threshold_raises(self):
        """Ratio exactly equal to 0.397 must also raise (strict inequality)."""
        # Hc = 0.397 * W^(1/3)
        W = 8.0  # W^(1/3) = 2
        Hc = 0.397 * (W ** (1.0 / 3.0))  # exactly 0.397 * 2 = 0.794
        with pytest.raises(ValueError):
            generate_grid(R=30.0, W=W, Hc=Hc)

    def test_valid_hob_does_not_raise(self):
        """Ratio > 0.397 must not raise."""
        generate_grid(R=30.0, W=VALID_W, Hc=VALID_Hc, width=4.0, height=4.0)


class TestGenerateGridOpenField:
    """Tests for generate_grid in open field mode (no width/height)."""

    def test_open_field_generates_points(self):
        """Open field mode should return a non-empty list."""
        points = generate_grid(R=50.0, W=VALID_W, Hc=VALID_Hc)
        assert len(points) > 0

    def test_open_field_max_offset_positive(self):
        """At least one point should have a non-zero offset (grid expands)."""
        points = generate_grid(R=50.0, W=VALID_W, Hc=VALID_Hc)
        max_offset = max(math.sqrt(p.dx**2 + p.dy**2) for p in points)
        assert max_offset > 0.0

    def test_open_field_all_angles_within_85(self):
        """No grid point should exceed 85° angle of incidence."""
        points = generate_grid(R=50.0, W=VALID_W, Hc=VALID_Hc)
        for p in points:
            assert p.alpha_deg <= 85.0 + 1e-9, (
                f"Point ({p.dx}, {p.dy}) has alpha={p.alpha_deg:.2f}° > 85°"
            )

    def test_open_field_includes_perpendicular(self):
        """Perpendicular point (dx=0, dy=0) must be in the open field grid."""
        points = generate_grid(R=50.0, W=VALID_W, Hc=VALID_Hc)
        perp = [p for p in points if p.dx == 0.0 and p.dy == 0.0]
        assert len(perp) == 1
