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

from ufc_blast.core.geometry import GridPoint, compute_point_geometry, determine_burst_type, generate_grid

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

    def test_building_grid_normal_point(self):
        """The normal-incidence point (dx=0, dy=0) lies on the facade when
        the facade reaches the burst height (Hc <= height)."""
        points = generate_grid(
            R=30.0, W=VALID_W, Hc=5.0,
            width=10.0, height=8.0, step=1.0,
        )
        centre_points = [p for p in points if p.dx == 0.0 and p.dy == 0.0]
        assert len(centre_points) == 1, "Exactly one normal point expected"
        cp = centre_points[0]
        assert math.isclose(cp.alpha_deg, 0.0, abs_tol=ABS_TOL)

    def test_building_grid_short_facade_below_burst_height(self):
        """A facade shorter than Hc lies entirely below the burst height."""
        points = generate_grid(
            R=30.0, W=VALID_W, Hc=VALID_Hc,  # Hc = 10 m
            width=10.0, height=8.0, step=1.0,
        )
        assert all(p.dy < 0 for p in points)
        assert not any(p.dx == 0.0 and p.dy == 0.0 for p in points)

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


class TestDetermineBurstType:
    """Tests for burst type auto-detection."""

    def test_air_burst(self):
        """Hc/W^(1/3) > 0.397 → 'air'."""
        burst_type, scaled_hob = determine_burst_type(Hc=10.0, W=1.0)
        assert burst_type == "air"
        assert scaled_hob > 0.397

    def test_surface_burst(self):
        """Hc/W^(1/3) <= 0.397 → 'surface'."""
        burst_type, scaled_hob = determine_burst_type(Hc=3.0, W=1000.0)
        assert burst_type == "surface"
        assert scaled_hob <= 0.397

    def test_exact_threshold_is_surface(self):
        """Ratio exactly equal to 0.397 → 'surface' (strict inequality)."""
        W = 8.0
        Hc = 0.397 * (W ** (1.0 / 3.0))
        burst_type, _ = determine_burst_type(Hc=Hc, W=W)
        assert burst_type == "surface"

    def test_zero_hc_is_surface(self):
        """Hc=0 → surface burst."""
        burst_type, scaled_hob = determine_burst_type(Hc=0.0, W=100.0)
        assert burst_type == "surface"
        assert scaled_hob == 0.0


class TestGenerateGridApplicability:
    """generate_grid now works for both air and surface bursts."""

    def test_surface_burst_no_longer_raises(self):
        """Surface burst parameters should not raise ValueError."""
        points = generate_grid(R=30.0, W=1000.0, Hc=3.0, width=10.0, height=8.0, step=2.0)
        assert len(points) > 0

    def test_valid_hob_still_works(self):
        """Air burst parameters continue to work."""
        points = generate_grid(R=30.0, W=1.0, Hc=10.0, width=4.0, height=4.0)
        assert len(points) > 0


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


class TestGroundPlane:
    """The facade grid must never extend below the ground plane (dy = -Hc)."""

    def test_building_grid_base_rests_on_ground(self):
        """Building mode: lowest row is at the ground plane, top at height-Hc."""
        points = generate_grid(
            R=30.0, W=200.0, Hc=5.0, width=10.0, height=8.0, step=1.0,
        )
        dys = sorted({p.dy for p in points})
        assert math.isclose(dys[0], -5.0, abs_tol=1e-12)
        assert math.isclose(dys[-1], 3.0, abs_tol=1e-12)
        assert all(p.dy >= -5.0 - 1e-12 for p in points)

    def test_building_grid_point_count_preserved(self):
        """Ground-based grid keeps (width/step+1) x (height/step+1) points."""
        points = generate_grid(
            R=30.0, W=200.0, Hc=5.0, width=30.0, height=24.0, step=0.25,
        )
        assert len(points) == 121 * 97

    def test_building_grid_contains_normal_point_when_facade_tall_enough(self):
        """dy=0 (burst height) lies on the facade iff Hc <= height."""
        points = generate_grid(
            R=30.0, W=200.0, Hc=5.0, width=10.0, height=8.0, step=1.0,
        )
        assert any(p.dx == 0.0 and p.dy == 0.0 for p in points)

    def test_open_field_clipped_at_ground(self):
        """Open-field expansion must not descend below the ground plane."""
        points = generate_grid(R=20.0, W=200.0, Hc=5.0, step=5.0)
        assert all(p.dy >= -5.0 - 1e-12 for p in points)
        # ...but still expands upward and sideways as before
        assert any(p.dy > 0 for p in points)
        assert any(p.dx > 0 for p in points)
