"""
Tests for the ``batch`` CLI subcommand.

Covers:
  - JSON output has correct keys
  - batch at alpha=0 matches point subcommand output
  - CSV output format (header + data rows)
  - Pressure decreases from alpha=0 to alpha=60
  - CSV file output with -o flag
"""

from __future__ import annotations

import csv
import io
import json
import os
import tempfile

import pytest

from ufc_blast.cli import _build_parser, _cmd_batch, _cmd_point


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_EXPECTED_KEYS = {
    "alpha_deg",
    "R_alpha_m",
    "Ps0_kPa",
    "C_alpha",
    "Pr_alpha_kPa",
    "ir_alpha_kPa_ms",
    "tA_ms",
    "t0_ms",
    "b",
}


def _parse_batch(extra_args: list[str] | None = None) -> object:
    """Parse batch CLI args with sensible defaults."""
    parser = _build_parser()
    base = [
        "batch",
        "--W", "200",
        "--R", "30",
        "--Hc", "5",
        "--alphas", "0,15,30,45,60",
    ]
    if extra_args:
        base.extend(extra_args)
    return parser.parse_args(base)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestBatchJsonOutput:
    """JSON output structure and correctness."""

    def test_json_has_correct_keys(self, capsys):
        args = _parse_batch(["--format", "json"])
        rc = _cmd_batch(args)
        assert rc == 0

        captured = capsys.readouterr()
        rows = json.loads(captured.out)
        assert len(rows) == 5
        for row in rows:
            assert _EXPECTED_KEYS <= set(row.keys())

    def test_alpha_zero_matches_point(self, capsys):
        """batch at alpha=0 should produce the same Pr as the point subcommand."""
        # Run batch with alpha=0 only
        args_batch = _parse_batch(["--alphas", "0", "--format", "json"])
        rc = _cmd_batch(args_batch)
        assert rc == 0
        batch_out = json.loads(capsys.readouterr().out)
        batch_pr = batch_out[0]["Pr_alpha_kPa"]

        # Run point subcommand
        parser = _build_parser()
        args_point = parser.parse_args([
            "point", "--W", "200", "--R", "30", "--Hc", "5", "--alpha", "0",
        ])
        rc = _cmd_point(args_point)
        assert rc == 0
        point_out = capsys.readouterr().out

        # Extract Pr from point textual output
        for line in point_out.splitlines():
            if "Peak reflected Pr" in line:
                # e.g.  "  Peak reflected Pr_α    : 123.456 kPa"
                pr_str = line.split(":")[1].strip().split()[0]
                point_pr = float(pr_str)
                break
        else:
            pytest.fail("Could not find Pr in point output")

        assert abs(batch_pr - point_pr) < 0.01

    def test_pressure_decreases_with_angle(self, capsys):
        """Reflected pressure at alpha=0 should exceed alpha=60."""
        args = _parse_batch(["--format", "json"])
        rc = _cmd_batch(args)
        assert rc == 0

        rows = json.loads(capsys.readouterr().out)
        pr_first = rows[0]["Pr_alpha_kPa"]   # alpha=0
        pr_last = rows[-1]["Pr_alpha_kPa"]    # alpha=60
        assert pr_first > pr_last


class TestBatchCsvOutput:
    """CSV output format tests."""

    def test_csv_stdout_format(self, capsys):
        args = _parse_batch(["--format", "csv"])
        rc = _cmd_batch(args)
        assert rc == 0

        captured = capsys.readouterr()
        reader = csv.DictReader(io.StringIO(captured.out))
        rows = list(reader)
        assert len(rows) == 5
        # Check header names
        assert _EXPECTED_KEYS <= set(reader.fieldnames)

    def test_csv_file_output(self, capsys, tmp_path):
        out_file = str(tmp_path / "batch_output.csv")
        args = _parse_batch(["--format", "csv", "-o", out_file])
        rc = _cmd_batch(args)
        assert rc == 0

        # Status message should go to stderr, not stdout
        captured = capsys.readouterr()
        assert out_file in captured.err

        # File should exist with correct content
        assert os.path.exists(out_file)
        with open(out_file) as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
        assert len(rows) == 5
        assert _EXPECTED_KEYS <= set(reader.fieldnames)
