"""
Interpolation engine for UFC 3-340-02 blast parameter tables.

Supports 1D log-log or linear interpolation (Table1D) and 2D bivariate
interpolation over a (family, x) grid (Table2D).
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
    """Generic 2D lookup table over a (family, x) grid.

    The x axis always uses linear interpolation.  The family axis uses
    log or linear interpolation depending on *method* passed to
    :meth:`lookup`.

    Parameters
    ----------
    family_levels : np.ndarray
        Sorted unique family-parameter values (e.g. Ps0 levels, Hc/W^⅓).
    x_tables : dict[float, Table1D]
        Mapping from each family level to a :class:`Table1D` keyed on the
        x variable → value.
    """

    family_levels: np.ndarray
    x_tables: dict = field(default_factory=dict)  # float -> Table1D

    def __post_init__(self) -> None:
        self.family_levels = np.asarray(self.family_levels, dtype=float)
        if not np.all(np.diff(self.family_levels) > 0):
            raise ValueError("family_levels must be strictly ascending")

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_csv(
        cls,
        path: str,
        family_col: str,
        x_col: str,
        value_col: str,
    ) -> "Table2D":
        """Load a Table2D from a long-format CSV file.

        Parameters
        ----------
        path : str
            Path to the CSV file.
        family_col : str
            Column name for the family parameter (e.g. ``ps0_kpa``,
            ``hc_scaled_m_kg13``).
        x_col : str
            Column name for the x (lookup) variable (e.g. ``angle_deg``,
            ``rg_scaled_m_kg13``).
        value_col : str
            Column name for the interpolated quantity.
        """
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)

        # Group rows by family level
        groups: dict[float, list[tuple[float, float]]] = {}
        for r in rows:
            family = float(r[family_col])
            x = float(r[x_col])
            value = float(r[value_col])
            groups.setdefault(family, []).append((x, value))

        family_levels = np.sort(np.array(list(groups.keys())))

        x_tables: dict[float, Table1D] = {}
        for fam in family_levels:
            pairs = sorted(groups[float(fam)], key=lambda t: t[0])
            xs = np.array([p[0] for p in pairs])
            values = np.array([p[1] for p in pairs])
            x_tables[float(fam)] = Table1D(x=xs, y=values)

        return cls(family_levels=family_levels, x_tables=x_tables)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _clamp_x(
        x: float,
        tbl: "Table1D | None" = None,
        *,
        x_min: float | None = None,
        x_max: float | None = None,
    ) -> float:
        """Clamp *x* to [x_min, x_max] with a warning if out of range."""
        if tbl is not None:
            x_min = tbl.x[0]
            x_max = tbl.x[-1]
        if x < x_min:
            warnings.warn(
                f"x={x} below table minimum {x_min}; clamping.",
                UserWarning,
                stacklevel=3,
            )
            return x_min
        if x > x_max:
            warnings.warn(
                f"x={x} above table maximum {x_max}; clamping.",
                UserWarning,
                stacklevel=3,
            )
            return x_max
        return x

    # ------------------------------------------------------------------
    # Lookup
    # ------------------------------------------------------------------

    def lookup(self, x: float, family: float, method: str = "log") -> float:
        """Bivariate interpolation.

        The *x* axis is always linearly interpolated.  The *family* axis is
        interpolated in log space when ``method='log'`` and in linear
        space when ``method='linear'``.

        *x* is clamped to the available range with a warning if it falls
        outside.  *family* outside the table range raises :class:`ValueError`.

        Parameters
        ----------
        x : float
            Lookup variable (e.g. angle in degrees, scaled range).
        family : float
            Family parameter (e.g. Ps0 in kPa, scaled Hc).
        method : str
            Interpolation method for the family axis (``'log'`` or ``'linear'``).

        Returns
        -------
        float
            Interpolated value.

        Raises
        ------
        ValueError
            If *family* is outside the table range.
        """
        fam_min, fam_max = self.family_levels[0], self.family_levels[-1]
        if family < fam_min or family > fam_max:
            raise ValueError(
                f"family={family} is outside table range [{fam_min}, {fam_max}]"
            )
        value, _ = self.lookup_flagged(x=x, family=family, method=method)
        return value

    def lookup_flagged(
        self, x: float, family: float, method: str = "log"
    ) -> tuple[float, bool]:
        """Bivariate interpolation with out-of-range flag.

        Identical to :meth:`lookup` except that when *family* falls outside
        the table range the value is clamped to the nearest boundary instead
        of raising :class:`ValueError`.

        Parameters
        ----------
        x : float
            Lookup variable (e.g. angle in degrees, scaled range).
        family : float
            Family parameter (e.g. Ps0 in kPa, scaled Hc).
        method : str
            Interpolation method for the family axis (``'log'`` or ``'linear'``).

        Returns
        -------
        tuple[float, bool]
            ``(value, out_of_range)`` where *out_of_range* is ``True`` when
            *family* was clamped to a table boundary.
        """
        fam_min, fam_max = float(self.family_levels[0]), float(self.family_levels[-1])
        out_of_range = bool(family < fam_min or family > fam_max)
        if out_of_range:
            family = float(np.clip(family, fam_min, fam_max))

        # Find bracketing family levels — same logic as lookup()
        idx = np.searchsorted(self.family_levels, family)

        if idx == 0:
            tbl = self.x_tables[float(self.family_levels[0])]
            x = self._clamp_x(x, tbl)
            return tbl.lookup(x, method="linear"), out_of_range

        if idx == len(self.family_levels):
            tbl = self.x_tables[float(self.family_levels[-1])]
            x = self._clamp_x(x, tbl)
            return tbl.lookup(x, method="linear"), out_of_range

        fam_lo = float(self.family_levels[idx - 1])
        fam_hi = float(self.family_levels[idx])

        # Clamp x to the overlapping range of the two bracketing families
        tbl_lo = self.x_tables[fam_lo]
        tbl_hi = self.x_tables[fam_hi]
        x_min = max(tbl_lo.x[0], tbl_hi.x[0])
        x_max = min(tbl_lo.x[-1], tbl_hi.x[-1])
        x = self._clamp_x(x, x_min=x_min, x_max=x_max)

        val_lo = tbl_lo.lookup(x, method="linear")
        val_hi = tbl_hi.lookup(x, method="linear")

        if method == "log":
            t = (np.log(family) - np.log(fam_lo)) / (np.log(fam_hi) - np.log(fam_lo))
            value = float(np.exp((1 - t) * np.log(val_lo) + t * np.log(val_hi)))
        elif method == "linear":
            t = (family - fam_lo) / (fam_hi - fam_lo)
            value = float((1 - t) * val_lo + t * val_hi)
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'log' or 'linear'.")

        return value, out_of_range

    def lookup_batch_flagged(
        self, xs: np.ndarray, families: np.ndarray, method: str = "log"
    ) -> tuple[np.ndarray, np.ndarray]:
        """Vectorized bivariate interpolation with out-of-range flags.

        Identical to :meth:`lookup_batch` but additionally returns a boolean
        mask indicating which queries had their *family* value clamped.

        Parameters
        ----------
        xs : np.ndarray
            X-axis query points.
        families : np.ndarray
            Family-axis query points.
        method : str
            ``'log'`` or ``'linear'`` for the family axis.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            ``(values, out_of_range)`` where *out_of_range* is a boolean
            array of the same shape as *xs*.
        """
        xs = np.asarray(xs, dtype=float)
        families = np.asarray(families, dtype=float)

        fam_min, fam_max = self.family_levels[0], self.family_levels[-1]
        out_of_range = (families < fam_min) | (families > fam_max)
        families_clamped = np.clip(families, fam_min, fam_max)

        result = np.empty_like(xs)

        # Find bracketing family indices for all points at once
        idxs = np.searchsorted(self.family_levels, families_clamped)
        idxs = np.clip(idxs, 1, len(self.family_levels) - 1)

        # Group by (fam_lo, fam_hi) pair for efficiency
        for lo_idx in range(len(self.family_levels) - 1):
            hi_idx = lo_idx + 1
            mask = idxs == hi_idx
            if not np.any(mask):
                continue

            fam_lo = float(self.family_levels[lo_idx])
            fam_hi = float(self.family_levels[hi_idx])
            tbl_lo = self.x_tables[fam_lo]
            tbl_hi = self.x_tables[fam_hi]

            # Clamp x to the overlapping range of the two bracketing families
            x_min = max(tbl_lo.x[0], tbl_hi.x[0])
            x_max = min(tbl_lo.x[-1], tbl_hi.x[-1])
            x_sub = np.clip(xs[mask], x_min, x_max)

            val_lo = tbl_lo.lookup_batch(x_sub, method="linear")
            val_hi = tbl_hi.lookup_batch(x_sub, method="linear")

            fam_sub = families_clamped[mask]
            if method == "log":
                t = (np.log(fam_sub) - np.log(fam_lo)) / (np.log(fam_hi) - np.log(fam_lo))
                result[mask] = np.exp((1 - t) * np.log(val_lo) + t * np.log(val_hi))
            else:
                t = (fam_sub - fam_lo) / (fam_hi - fam_lo)
                result[mask] = (1 - t) * val_lo + t * val_hi

        return result, out_of_range

    def lookup_batch(
        self, xs: np.ndarray, families: np.ndarray, method: str = "log"
    ) -> np.ndarray:
        """Vectorized bivariate interpolation for arrays of (x, family).

        Parameters
        ----------
        xs : np.ndarray
            X-axis query points.
        families : np.ndarray
            Family-axis query points.
        method : str
            ``'log'`` or ``'linear'`` for the family axis.

        Returns
        -------
        np.ndarray
            Interpolated values, same shape as inputs.
        """
        values, _ = self.lookup_batch_flagged(xs, families, method=method)
        return values
