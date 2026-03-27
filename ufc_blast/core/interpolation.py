"""
Interpolation engine for UFC 3-340-02 blast parameter tables.

Supports 1D log-log or linear interpolation (Table1D) and 2D bivariate
interpolation over a (Ps0, angle) grid (Table2D).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import csv

import numpy as np


@dataclass
class Table1D:
    """1D lookup table with configurable interpolation method.

    Parameters
    ----------
    x : np.ndarray
        Independent variable values, sorted ascending.
    y : np.ndarray
        Dependent variable values corresponding to *x*.
    """

    x: np.ndarray
    y: np.ndarray

    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        if self.x.shape != self.y.shape:
            raise ValueError("x and y must have the same length")
        if self.x.ndim != 1:
            raise ValueError("x and y must be 1-D arrays")
        if not np.all(np.diff(self.x) > 0):
            raise ValueError("x values must be strictly ascending")

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_csv(cls, path: str, x_col: str, y_col: str) -> "Table1D":
        """Load a Table1D from a CSV file.

        Parameters
        ----------
        path : str
            Path to the CSV file.
        x_col : str
            Column name for the independent variable.
        y_col : str
            Column name for the dependent variable.
        """
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)

        x_vals = np.array([float(r[x_col]) for r in rows])
        y_vals = np.array([float(r[y_col]) for r in rows])

        # Sort by x ascending
        order = np.argsort(x_vals)
        return cls(x=x_vals[order], y=y_vals[order])

    # ------------------------------------------------------------------
    # Lookup
    # ------------------------------------------------------------------

    def lookup(self, x_val: float, method: str = "log") -> float:
        """Interpolate the table at *x_val*.

        Parameters
        ----------
        x_val : float
            Query point.
        method : str
            ``'log'`` performs log-log interpolation (linear in log space);
            ``'linear'`` performs standard linear interpolation.

        Returns
        -------
        float
            Interpolated value.

        Raises
        ------
        ValueError
            If *x_val* is outside the table range.
        """
        x_min, x_max = self.x[0], self.x[-1]
        if x_val < x_min or x_val > x_max:
            raise ValueError(
                f"x_val={x_val} is outside table range [{x_min}, {x_max}]"
            )

        if method == "log":
            log_x = np.log(self.x)
            log_y = np.log(self.y)
            log_result = np.interp(np.log(x_val), log_x, log_y)
            return float(np.exp(log_result))
        elif method == "linear":
            return float(np.interp(x_val, self.x, self.y))
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'log' or 'linear'.")

    def lookup_batch(self, x_vals: np.ndarray, method: str = "log") -> np.ndarray:
        """Vectorized interpolation for an array of query points.

        Parameters
        ----------
        x_vals : np.ndarray
            Array of query points.
        method : str
            ``'log'`` or ``'linear'``.

        Returns
        -------
        np.ndarray
            Interpolated values, same shape as *x_vals*.
        """
        x_vals = np.asarray(x_vals, dtype=float)
        if method == "log":
            log_x = np.log(self.x)
            log_y = np.log(self.y)
            return np.exp(np.interp(np.log(x_vals), log_x, log_y))
        else:
            return np.interp(x_vals, self.x, self.y)


@dataclass
class Table2D:
    """2D lookup table over a (Ps0, angle) grid.

    The angle axis always uses linear interpolation.  The Ps0 axis uses
    log or linear interpolation depending on *method* passed to
    :meth:`lookup`.

    Parameters
    ----------
    ps0_levels : np.ndarray
        Sorted unique Ps0 values present in the table (kPa).
    angle_tables : dict[float, Table1D]
        Mapping from each Ps0 level to a :class:`Table1D` keyed on angle
        (degrees) → value.
    """

    ps0_levels: np.ndarray
    angle_tables: dict = field(default_factory=dict)  # float -> Table1D

    def __post_init__(self) -> None:
        self.ps0_levels = np.asarray(self.ps0_levels, dtype=float)
        if not np.all(np.diff(self.ps0_levels) > 0):
            raise ValueError("ps0_levels must be strictly ascending")

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_csv(
        cls,
        path: str,
        ps0_col: str,
        angle_col: str,
        value_col: str,
    ) -> "Table2D":
        """Load a Table2D from a long-format CSV file.

        Parameters
        ----------
        path : str
            Path to the CSV file.
        ps0_col : str
            Column name for Ps0 values.
        angle_col : str
            Column name for angle values (degrees).
        value_col : str
            Column name for the interpolated quantity.
        """
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)

        # Group rows by ps0 level
        groups: dict[float, list[tuple[float, float]]] = {}
        for r in rows:
            ps0 = float(r[ps0_col])
            angle = float(r[angle_col])
            value = float(r[value_col])
            groups.setdefault(ps0, []).append((angle, value))

        ps0_levels = np.sort(np.array(list(groups.keys())))

        angle_tables: dict[float, Table1D] = {}
        for ps0 in ps0_levels:
            pairs = sorted(groups[float(ps0)], key=lambda t: t[0])
            angles = np.array([p[0] for p in pairs])
            values = np.array([p[1] for p in pairs])
            angle_tables[float(ps0)] = Table1D(x=angles, y=values)

        return cls(ps0_levels=ps0_levels, angle_tables=angle_tables)

    # ------------------------------------------------------------------
    # Lookup
    # ------------------------------------------------------------------

    def lookup(self, angle: float, ps0: float, method: str = "log") -> float:
        """Bivariate interpolation.

        The angle axis is always linearly interpolated.  The Ps0 axis is
        interpolated in log space when ``method='log'`` and in linear
        space when ``method='linear'``.

        Angle is clamped to the available range with a warning if it falls
        outside.  Ps0 outside the table range raises :class:`ValueError`.

        Parameters
        ----------
        angle : float
            Angle of incidence (degrees).
        ps0 : float
            Peak incident overpressure (kPa).
        method : str
            Interpolation method for the Ps0 axis (``'log'`` or ``'linear'``).

        Returns
        -------
        float
            Interpolated value.

        Raises
        ------
        ValueError
            If *ps0* is outside the table range.
        """
        ps0_min, ps0_max = self.ps0_levels[0], self.ps0_levels[-1]
        if ps0 < ps0_min or ps0 > ps0_max:
            raise ValueError(
                f"ps0={ps0} is outside table range [{ps0_min}, {ps0_max}]"
            )

        # Determine angle bounds across all ps0 slices for clamping
        all_angle_mins = [tbl.x[0] for tbl in self.angle_tables.values()]
        all_angle_maxs = [tbl.x[-1] for tbl in self.angle_tables.values()]
        angle_min = max(all_angle_mins)  # conservative: highest lower bound
        angle_max = min(all_angle_maxs)  # conservative: lowest upper bound

        if angle < angle_min:
            warnings.warn(
                f"angle={angle} below table minimum {angle_min}; clamping.",
                UserWarning,
                stacklevel=2,
            )
            angle = angle_min
        elif angle > angle_max:
            warnings.warn(
                f"angle={angle} above table maximum {angle_max}; clamping.",
                UserWarning,
                stacklevel=2,
            )
            angle = angle_max

        # Find bracketing ps0 levels
        idx = np.searchsorted(self.ps0_levels, ps0)

        if idx == 0:
            # Exact match at lower boundary
            return self.angle_tables[float(self.ps0_levels[0])].lookup(angle, method="linear")

        if idx == len(self.ps0_levels):
            # Exact match at upper boundary
            return self.angle_tables[float(self.ps0_levels[-1])].lookup(angle, method="linear")

        ps0_lo = float(self.ps0_levels[idx - 1])
        ps0_hi = float(self.ps0_levels[idx])

        val_lo = self.angle_tables[ps0_lo].lookup(angle, method="linear")
        val_hi = self.angle_tables[ps0_hi].lookup(angle, method="linear")

        if method == "log":
            # Log-linear interpolation on the Ps0 axis
            t = (np.log(ps0) - np.log(ps0_lo)) / (np.log(ps0_hi) - np.log(ps0_lo))
            return float(np.exp((1 - t) * np.log(val_lo) + t * np.log(val_hi)))
        elif method == "linear":
            t = (ps0 - ps0_lo) / (ps0_hi - ps0_lo)
            return float((1 - t) * val_lo + t * val_hi)
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'log' or 'linear'.")

    def lookup_batch(
        self, angles: np.ndarray, ps0s: np.ndarray, method: str = "log"
    ) -> np.ndarray:
        """Vectorized bivariate interpolation for arrays of (angle, ps0).

        Parameters
        ----------
        angles : np.ndarray
            Angles of incidence (degrees).
        ps0s : np.ndarray
            Peak incident overpressures (kPa).
        method : str
            ``'log'`` or ``'linear'`` for the Ps0 axis.

        Returns
        -------
        np.ndarray
            Interpolated values, same shape as inputs.
        """
        angles = np.asarray(angles, dtype=float)
        ps0s = np.asarray(ps0s, dtype=float)
        result = np.empty_like(angles)

        # Clamp angles to safe range
        all_angle_maxs = [tbl.x[-1] for tbl in self.angle_tables.values()]
        angle_max = min(all_angle_maxs)
        angles_clamped = np.clip(angles, 0.0, angle_max)

        # Find bracketing ps0 indices for all points at once
        idxs = np.searchsorted(self.ps0_levels, ps0s)
        idxs = np.clip(idxs, 1, len(self.ps0_levels) - 1)

        ps0_lo_arr = self.ps0_levels[idxs - 1]
        ps0_hi_arr = self.ps0_levels[idxs]

        # Look up values at lo and hi ps0 for each point
        # Group by (ps0_lo, ps0_hi) pair for efficiency
        for lo_idx in range(len(self.ps0_levels) - 1):
            hi_idx = lo_idx + 1
            mask = idxs == hi_idx
            if not np.any(mask):
                continue

            ps0_lo = float(self.ps0_levels[lo_idx])
            ps0_hi = float(self.ps0_levels[hi_idx])
            tbl_lo = self.angle_tables[ps0_lo]
            tbl_hi = self.angle_tables[ps0_hi]

            a_sub = angles_clamped[mask]
            val_lo = tbl_lo.lookup_batch(a_sub, method="linear")
            val_hi = tbl_hi.lookup_batch(a_sub, method="linear")

            ps0_sub = ps0s[mask]
            if method == "log":
                t = (np.log(ps0_sub) - np.log(ps0_lo)) / (np.log(ps0_hi) - np.log(ps0_lo))
                result[mask] = np.exp((1 - t) * np.log(val_lo) + t * np.log(val_hi))
            else:
                t = (ps0_sub - ps0_lo) / (ps0_hi - ps0_lo)
                result[mask] = (1 - t) * val_lo + t * val_hi

        return result
