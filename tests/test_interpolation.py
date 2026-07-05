"""
Tests for the UFC 3-340-02 interpolation engine (Table1D and Table2D).

All tests use the actual digitized data files from data/free_air/.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from ufc_blast.core.interpolation import Table1D, Table2D

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).parent.parent / "data" / "free_air"

INCIDENT_CSV = DATA_DIR / "figure_2_7_incident.csv"
T0_CSV = DATA_DIR / "figure_2_7_t0.csv"
TA_CSV = DATA_DIR / "figure_2_7_tA.csv"
WPD_PS0_CSV = DATA_DIR / "figure_2_7_wpd_ps0.csv"
CALPHA_CSV = DATA_DIR / "figure_2_193_calpha.csv"
IRALPHA_CSV = DATA_DIR / "figure_2_194_iralpha.csv"


# ===========================================================================
# Table1D tests
# ===========================================================================


class TestTable1DFromCSV:
    """Loading and structural sanity checks."""

    def test_loads_incident_csv(self):
        tbl = Table1D.from_csv(INCIDENT_CSV, "z_m_kg13", "ps0_kpa")
        assert len(tbl.x) == 53
        assert len(tbl.y) == 53

    def test_x_is_sorted_ascending(self):
        tbl = Table1D.from_csv(INCIDENT_CSV, "z_m_kg13", "ps0_kpa")
        assert np.all(np.diff(tbl.x) > 0)

    def test_ps0_decreases_as_z_increases(self):
        """Overpressure must fall off with scaled distance."""
        tbl = Table1D.from_csv(INCIDENT_CSV, "z_m_kg13", "ps0_kpa")
        assert np.all(np.diff(tbl.y) < 0), "Ps0 should be monotonically decreasing with Z"


class TestTable1DLookup:
    """Interpolation behaviour."""

    @pytest.fixture(scope="class")
    def incident_table(self):
        return Table1D.from_csv(INCIDENT_CSV, "z_m_kg13", "ps0_kpa")

    def test_exact_point_returns_table_value(self, incident_table):
        """Lookup at an exact grid point must reproduce the stored value."""
        idx = 10  # arbitrary interior point
        x_exact = incident_table.x[idx]
        y_expected = incident_table.y[idx]
        result = incident_table.lookup(x_exact, method="log")
        assert abs(result - y_expected) / y_expected < 1e-6, (
            f"Exact-point lookup: expected {y_expected}, got {result}"
        )

    def test_interpolated_value_between_neighbors(self, incident_table):
        """Midpoint lookup must fall strictly between its two neighbors."""
        idx = 20
        x_lo = incident_table.x[idx]
        x_hi = incident_table.x[idx + 1]
        y_lo = incident_table.y[idx]
        y_hi = incident_table.y[idx + 1]
        x_mid = np.sqrt(x_lo * x_hi)  # geometric midpoint for log space

        result = incident_table.lookup(x_mid, method="log")
        lo, hi = min(y_lo, y_hi), max(y_lo, y_hi)
        assert lo < result < hi, (
            f"Interpolated value {result} not between {lo} and {hi}"
        )

    def test_log_and_linear_give_different_results(self, incident_table):
        """Log-log and linear interpolation must diverge for non-linear data."""
        idx = 15
        x_lo = incident_table.x[idx]
        x_hi = incident_table.x[idx + 1]
        x_mid = (x_lo + x_hi) / 2.0

        result_log = incident_table.lookup(x_mid, method="log")
        result_lin = incident_table.lookup(x_mid, method="linear")
        assert result_log != pytest.approx(result_lin, rel=1e-4), (
            "Log and linear interpolation should differ for curved blast data"
        )

    def test_out_of_range_below_raises(self, incident_table):
        with pytest.raises(ValueError, match="outside table range"):
            incident_table.lookup(incident_table.x[0] - 1e-6, method="log")

    def test_out_of_range_above_raises(self, incident_table):
        with pytest.raises(ValueError, match="outside table range"):
            incident_table.lookup(incident_table.x[-1] + 1e-6, method="log")

    def test_out_of_range_linear_raises(self, incident_table):
        with pytest.raises(ValueError, match="outside table range"):
            incident_table.lookup(0.0, method="linear")

    def test_linear_exact_point_matches_table(self, incident_table):
        idx = 5
        x_exact = incident_table.x[idx]
        y_expected = incident_table.y[idx]
        result = incident_table.lookup(x_exact, method="linear")
        assert abs(result - y_expected) / y_expected < 1e-9

    def test_t0_table_loads_and_interpolates(self):
        """t0 (positive phase duration) table should load without error."""
        tbl = Table1D.from_csv(T0_CSV, "z_m_kg13", "t0_scaled_ms_kg13")
        assert len(tbl.x) == 167
        # Internal point lookup should not raise
        mid_idx = len(tbl.x) // 2
        val = tbl.lookup(tbl.x[mid_idx], method="log")
        assert val > 0

    def test_ta_table_loads_and_interpolates(self):
        """Arrival time table should load without error."""
        tbl = Table1D.from_csv(TA_CSV, "z_m_kg13", "tA_scaled_ms_kg13")
        assert len(tbl.x) == 123
        mid_idx = len(tbl.x) // 2
        val = tbl.lookup(tbl.x[mid_idx], method="log")
        assert val > 0


# ===========================================================================
# Table2D tests
# ===========================================================================


class TestTable2DFromCSV:
    """Loading and structural sanity checks for calpha."""

    def test_loads_calpha_csv(self):
        tbl = Table2D.from_csv(CALPHA_CSV, family_col="ps0_kpa", x_col="angle_deg", value_col="calpha")
        assert len(tbl.family_levels) == 20

    def test_family_levels_sorted(self):
        tbl = Table2D.from_csv(CALPHA_CSV, family_col="ps0_kpa", x_col="angle_deg", value_col="calpha")
        assert np.all(np.diff(tbl.family_levels) > 0)

    def test_loads_iralpha_csv(self):
        tbl = Table2D.from_csv(IRALPHA_CSV, family_col="ps0_kpa", x_col="angle_deg", value_col="iralpha_kpa_ms_kg13")
        assert len(tbl.family_levels) > 0


class TestTable2DLookup:
    """Bivariate interpolation behaviour for calpha."""

    @pytest.fixture(scope="class")
    def calpha(self):
        return Table2D.from_csv(CALPHA_CSV, family_col="ps0_kpa", x_col="angle_deg", value_col="calpha")

    # --- Known corner values from figure 2-193 ---

    def test_calpha_angle0_low_ps0_approx_2(self, calpha):
        """At angle=0 and low pressure, calpha ≈ 2.0 (near-normal reflection)."""
        result = calpha.lookup(x=0.0, family=1.378952)
        assert result == pytest.approx(1.9995, abs=0.01), (
            f"Expected ~2.0, got {result}"
        )

    def test_calpha_angle0_high_ps0_approx_12_24(self, calpha):
        """At angle=0 and very high pressure, calpha ≈ 12.24."""
        result = calpha.lookup(x=0.0, family=34473.8)
        assert result == pytest.approx(12.24, abs=0.05), (
            f"Expected ~12.24, got {result}"
        )

    def test_calpha_decreases_with_angle(self, calpha):
        """Reflection coefficient must decrease as angle of incidence increases."""
        ps0_mid = 344.74  # interior pressure level
        angles = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
        values = [calpha.lookup(x=a, family=ps0_mid) for a in angles]
        # Allow for small numerical noise but expect a clear downward trend
        assert values[0] > values[-1], (
            "calpha at angle=0 should exceed calpha at angle=90"
        )
        # Values should be generally non-increasing
        diffs = np.diff(values)
        assert np.sum(diffs > 0.1) == 0, (
            f"calpha should not increase significantly with angle: {values}"
        )

    def test_calpha_interpolation_between_ps0_levels(self, calpha):
        """Interpolated value at midpoint ps0 must lie between bracket values."""
        ps0_lo = 1.378952
        ps0_hi = 3.44738
        # Geometric midpoint in log space
        ps0_mid = np.exp((np.log(ps0_lo) + np.log(ps0_hi)) / 2.0)

        val_lo = calpha.lookup(x=0.0, family=ps0_lo)
        val_hi = calpha.lookup(x=0.0, family=ps0_hi)
        val_mid = calpha.lookup(x=0.0, family=ps0_mid)

        lo, hi = min(val_lo, val_hi), max(val_lo, val_hi)
        assert lo <= val_mid <= hi, (
            f"Interpolated value {val_mid} not between {lo} and {hi}"
        )

    def test_calpha_ps0_out_of_range_below_raises(self, calpha):
        with pytest.raises(ValueError, match="outside table range"):
            calpha.lookup(x=0.0, family=0.001)

    def test_calpha_ps0_out_of_range_above_raises(self, calpha):
        with pytest.raises(ValueError, match="outside table range"):
            calpha.lookup(x=0.0, family=1e9)

    def test_calpha_angle_clamped_with_warning(self, calpha):
        """Angle slightly out of range should warn and clamp, not raise."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calpha.lookup(x=95.0, family=1.378952)
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)
            assert "clamping" in str(w[0].message).lower()
        assert result > 0

    def test_calpha_exact_ps0_level_returns_table_value(self, calpha):
        """Lookup at an exact ps0 grid point should reproduce the table value."""
        ps0_exact = 137.8952  # exact level present in data (index 6)
        # Find the actual level
        ps0_actual = calpha.family_levels[6]
        x_tbl = calpha.x_tables[float(ps0_actual)]
        y_expected = x_tbl.lookup(0.0, method="linear")
        result = calpha.lookup(x=0.0, family=float(ps0_actual))
        assert result == pytest.approx(y_expected, rel=1e-6)

    def test_log_vs_linear_ps0_interpolation_differ(self, calpha):
        """Log and linear ps0 interpolation should yield different values."""
        ps0_lo = 1.378952
        ps0_hi = 3.44738
        ps0_mid = (ps0_lo + ps0_hi) / 2.0  # arithmetic midpoint favours linear

        result_log = calpha.lookup(x=0.0, family=ps0_mid, method="log")
        result_lin = calpha.lookup(x=0.0, family=ps0_mid, method="linear")
        assert result_log != pytest.approx(result_lin, rel=1e-4)


class TestTable2DFlagged:
    """Tests for extrapolation-aware lookup methods."""

    @pytest.fixture(scope="class")
    def calpha(self):
        return Table2D.from_csv(CALPHA_CSV, family_col="ps0_kpa", x_col="angle_deg", value_col="calpha")

    # -- lookup_flagged --

    def test_lookup_flagged_in_range_returns_false(self, calpha):
        """In-range family should return (value, False)."""
        val, flag = calpha.lookup_flagged(x=0.0, family=6.89476)
        assert flag is False
        assert val > 0.0

    def test_lookup_flagged_below_range_returns_true(self, calpha):
        """Family below table min should return (value, True), not raise."""
        fam_min = calpha.family_levels[0]
        val, flag = calpha.lookup_flagged(x=0.0, family=fam_min * 0.5)
        assert flag is True
        assert val > 0.0

    def test_lookup_flagged_above_range_returns_true(self, calpha):
        """Family above table max should return (value, True), not raise."""
        fam_max = calpha.family_levels[-1]
        val, flag = calpha.lookup_flagged(x=0.0, family=fam_max * 2.0)
        assert flag is True
        assert val > 0.0

    def test_lookup_flagged_matches_lookup_for_in_range(self, calpha):
        """Flagged and unflagged must return identical values for in-range queries."""
        val_flagged, flag = calpha.lookup_flagged(x=10.0, family=6.89476)
        val_plain = calpha.lookup(x=10.0, family=6.89476)
        assert val_flagged == pytest.approx(val_plain, rel=1e-12)
        assert flag is False

    # -- lookup_batch_flagged --

    def test_lookup_batch_flagged_mixed(self, calpha):
        """Batch with mixed in-range and out-of-range families."""
        fam_min = calpha.family_levels[0]
        fam_max = calpha.family_levels[-1]
        xs = np.array([0.0, 0.0, 0.0])
        families = np.array([fam_min * 0.5, 6.89476, fam_max * 2.0])
        vals, flags = calpha.lookup_batch_flagged(xs, families)
        assert flags[0] is True or flags[0] == True  # below range
        assert flags[1] is False or flags[1] == False  # in range
        assert flags[2] is True or flags[2] == True  # above range
        assert all(v > 0.0 for v in vals)

    def test_lookup_batch_flagged_all_in_range(self, calpha):
        """All in-range families should produce all-False mask."""
        xs = np.array([0.0, 10.0, 20.0])
        families = np.array([6.89476, 68.9476, 689.476])
        vals, flags = calpha.lookup_batch_flagged(xs, families)
        assert not np.any(flags)

    def test_lookup_batch_flagged_values_match_unflagged(self, calpha):
        """Flagged batch values must match unflagged batch for same inputs."""
        xs = np.array([0.0, 10.0, 20.0])
        families = np.array([6.89476, 68.9476, 689.476])
        vals_flagged, _ = calpha.lookup_batch_flagged(xs, families)
        vals_plain = calpha.lookup_batch(xs, families)
        np.testing.assert_allclose(vals_flagged, vals_plain, rtol=1e-12)


class TestTable1DFlagged:
    """Flagged 1D lookups: clamp to the table boundary and report it."""

    @pytest.fixture(scope="class")
    def incident_table(self):
        return Table1D.from_csv(INCIDENT_CSV, "z_m_kg13", "ps0_kpa")

    def test_in_range_matches_lookup_and_no_flag(self, incident_table):
        x = float(incident_table.x[5])
        val, flag = incident_table.lookup_flagged(x)
        assert val == pytest.approx(incident_table.lookup(x))
        assert flag is False

    def test_above_range_clamps_and_flags(self, incident_table):
        x_max = float(incident_table.x[-1])
        val, flag = incident_table.lookup_flagged(x_max * 10.0)
        assert val == pytest.approx(float(incident_table.y[-1]))
        assert flag is True

    def test_below_range_clamps_and_flags(self, incident_table):
        x_min = float(incident_table.x[0])
        val, flag = incident_table.lookup_flagged(x_min / 10.0)
        assert val == pytest.approx(float(incident_table.y[0]))
        assert flag is True

    def test_batch_flagged_matches_scalar(self, incident_table):
        xs = np.array([
            float(incident_table.x[0]) / 2.0,
            float(incident_table.x[7]),
            float(incident_table.x[-1]) * 3.0,
        ])
        vals, flags = incident_table.lookup_batch_flagged(xs)
        for i, x in enumerate(xs):
            v, f = incident_table.lookup_flagged(float(x))
            assert vals[i] == pytest.approx(v)
            assert bool(flags[i]) == f
        assert list(flags) == [True, False, True]
