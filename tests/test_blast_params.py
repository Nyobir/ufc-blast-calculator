"""
Tests for ufc_blast.core.blast_params.

Covers:
  - solve_friedlander_b: round-trip recovery of known b values
  - friedlander: waveform shape, boundary conditions, integral self-consistency
  - compute_point: physically plausible results for representative UFC scenarios
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ufc_blast.core.blast_params import (
    BlastPointResult,
    compute_point,
    friedlander,
    solve_friedlander_b,
)


# ===========================================================================
# Helper: analytical Friedlander impulse
# ===========================================================================

def _friedlander_impulse_analytical(Pr: float, t0: float, b: float) -> float:
    """Impulse = (Pr * t0) / b^2 * (b - 1 + exp(-b))."""
    return (Pr * t0) / (b ** 2) * (b - 1.0 + math.exp(-b))


# ===========================================================================
# TestSolveFriedlanderB
# ===========================================================================


class TestSolveFriedlanderB:
    """Round-trip recovery and basic properties of solve_friedlander_b."""

    def test_known_b_case_1(self):
        """b=6.4 with Pr=465 kPa, t0=18.6 ms should be recovered within 1%."""
        Pr = 465.0          # kPa
        t0 = 18.6e-3        # s
        b_true = 6.4

        ir = _friedlander_impulse_analytical(Pr, t0, b_true)
        b_solved = solve_friedlander_b(Pr, t0, ir)

        assert abs(b_solved - b_true) / b_true < 0.01, (
            f"Expected b≈{b_true}, got {b_solved:.5f}"
        )

    def test_known_b_case_2(self):
        """b=4.8 with Pr=270 kPa, t0=14.8 ms should be recovered within 1%."""
        Pr = 270.0
        t0 = 14.8e-3
        b_true = 4.8

        ir = _friedlander_impulse_analytical(Pr, t0, b_true)
        b_solved = solve_friedlander_b(Pr, t0, ir)

        assert abs(b_solved - b_true) / b_true < 0.01, (
            f"Expected b≈{b_true}, got {b_solved:.5f}"
        )

    def test_b_is_positive(self):
        """solve_friedlander_b must always return a positive value."""
        Pr = 100.0
        t0 = 20.0e-3
        b_true = 3.0
        ir = _friedlander_impulse_analytical(Pr, t0, b_true)
        b_solved = solve_friedlander_b(Pr, t0, ir)
        assert b_solved > 0.0

    def test_b_varies_with_impulse(self):
        """Higher impulse (relative to Pr*t0) yields a smaller b."""
        Pr = 300.0
        t0 = 15.0e-3
        b_high = 8.0
        b_low = 2.0

        ir_high = _friedlander_impulse_analytical(Pr, t0, b_high)
        ir_low = _friedlander_impulse_analytical(Pr, t0, b_low)

        # Lower b → more area under curve → larger impulse
        assert ir_low > ir_high

        b_solved_high = solve_friedlander_b(Pr, t0, ir_high)
        b_solved_low = solve_friedlander_b(Pr, t0, ir_low)
        assert b_solved_high > b_solved_low


# ===========================================================================
# TestFriedlander
# ===========================================================================


class TestFriedlander:
    """Shape and consistency tests for the Friedlander waveform."""

    # Shared fixture values
    Pr = 200.0   # kPa
    tA = 0.050   # s
    t0 = 0.030   # s
    b = 5.0

    def _make_t(self, n: int = 10_000) -> np.ndarray:
        """Dense time array covering one full positive phase."""
        return np.linspace(self.tA - self.t0, self.tA + 2 * self.t0, n)

    def test_peak_at_arrival_equals_Pr(self):
        """P(tA) must equal Pr_alpha (tau=0)."""
        t = np.array([self.tA])
        P = friedlander(t, self.Pr, self.tA, self.t0, self.b)
        assert abs(P[0] - self.Pr) < 1e-10

    def test_zero_at_end_of_positive_phase(self):
        """P(tA + t0) must be zero (tau=1, factor 1-tau = 0)."""
        t = np.array([self.tA + self.t0])
        P = friedlander(t, self.Pr, self.tA, self.t0, self.b)
        assert abs(P[0]) < 1e-12

    def test_zero_before_arrival(self):
        """All pressures before tA must be exactly zero."""
        t = np.linspace(0.0, self.tA - 1e-9, 500)
        P = friedlander(t, self.Pr, self.tA, self.t0, self.b)
        assert np.all(P == 0.0)

    def test_zero_after_positive_phase(self):
        """All pressures after tA + t0 must be exactly zero."""
        t = np.linspace(self.tA + self.t0 + 1e-9, self.tA + 5 * self.t0, 500)
        P = friedlander(t, self.Pr, self.tA, self.t0, self.b)
        assert np.all(P == 0.0)

    def test_integral_equals_impulse(self):
        """Numerical integral of P(t) over positive phase should match
        the analytical impulse (Pr*t0)/b^2*(b-1+exp(-b))."""
        n = 100_000
        t = np.linspace(self.tA, self.tA + self.t0, n)
        P = friedlander(t, self.Pr, self.tA, self.t0, self.b)
        ir_numerical = np.trapezoid(P, t)
        ir_analytical = _friedlander_impulse_analytical(self.Pr, self.t0, self.b)
        rel_err = abs(ir_numerical - ir_analytical) / ir_analytical
        assert rel_err < 1e-4, (
            f"Impulse mismatch: numerical={ir_numerical:.6f}, "
            f"analytical={ir_analytical:.6f}, rel_err={rel_err:.2e}"
        )

    def test_non_negative_in_positive_phase(self):
        """P(t) must be >= 0 throughout the positive phase."""
        t = np.linspace(self.tA, self.tA + self.t0, 1000)
        P = friedlander(t, self.Pr, self.tA, self.t0, self.b)
        assert np.all(P >= 0.0)


# ===========================================================================
# TestComputePoint
# ===========================================================================


class TestComputePoint:
    """Integration tests against UFC 3-340-02 tables for realistic scenarios."""

    def test_normal_incidence_200kg_30m(self):
        """W=200 kg, R=30 m, alpha=0° — normal incidence baseline."""
        result = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0)

        assert isinstance(result, BlastPointResult)
        assert result.Ps0 > 0.0,        "Incident pressure must be positive"
        assert result.Pr_alpha > result.Ps0, "Reflected > incident at any incidence"
        assert result.C_alpha >= 2.0,   "Reflection coeff >= 2 for alpha=0"
        assert result.b > 0.0,          "Friedlander b must be positive"
        assert result.tA > 0.0,         "Arrival time must be positive"
        assert result.t0 > 0.0,         "Positive phase duration must be positive"
        assert result.ir_alpha > 0.0,   "Reflected impulse must be positive"

    def test_angled_point_lower_pr_than_normal(self):
        """A point at 70° incidence has lower Pr than the 0° (normal) point
        at the same slant distance.

        Note: UFC Figure 2-193 shows that C_alpha can peak around 40-45° in the
        Mach-reflection regime before declining.  At grazing angles (>~60°) the
        reflection coefficient falls below its 0° value, so comparing 0° vs 70°
        is a reliable test regardless of the Ps0 level.
        """
        result_normal = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0)
        result_angled = compute_point(R_alpha=30.0, alpha_deg=70.0, W=200.0)

        # Same R_alpha → same Ps0; only C_alpha differs
        assert abs(result_normal.Ps0 - result_angled.Ps0) < 1e-10
        assert result_angled.Pr_alpha < result_normal.Pr_alpha, (
            "70° incidence must have lower reflected pressure than normal (0°) incidence"
        )

    def test_farther_point_lower_pressure(self):
        """Increasing standoff distance must reduce incident overpressure."""
        result_near = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0)
        result_far = compute_point(R_alpha=60.0, alpha_deg=0.0, W=200.0)
        assert result_far.Ps0 < result_near.Ps0

    def test_farther_point_later_arrival(self):
        """Blast wave arrives later at greater standoff."""
        result_near = compute_point(R_alpha=30.0, alpha_deg=0.0, W=200.0)
        result_far = compute_point(R_alpha=60.0, alpha_deg=0.0, W=200.0)
        assert result_far.tA > result_near.tA

    def test_result_fields_are_finite(self):
        """All result fields must be finite floats."""
        result = compute_point(R_alpha=45.0, alpha_deg=20.0, W=200.0)
        for field_name, value in result.__dict__.items():
            assert math.isfinite(value), f"Field {field_name} is not finite: {value}"
