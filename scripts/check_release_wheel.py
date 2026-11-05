"""Release gate: verify wheel CSV payload and run installed-code smoke checks.

Usage: python scripts/check_release_wheel.py dist/ufc_blast-*.whl
Run with a Python environment containing the package's core dependencies.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python scripts/check_release_wheel.py WHEEL.whl", file=sys.stderr)
        return 2
    wheel = Path(sys.argv[1]).resolve()
    source = Path(__file__).resolve().parent.parent / "data"
    expected = {
        f"ufc_blast/data/{folder}/{p.name}": p
        for folder in ("free_air", "surface")
        for p in (source / folder).glob("*.csv")
    }
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        missing = sorted(expected.keys() - names)
        if missing:
            print("Missing CSV payload in wheel: " + ", ".join(missing), file=sys.stderr)
            return 1
        for member, path in expected.items():
            if hashlib.sha256(archive.read(member)).digest() != hashlib.sha256(path.read_bytes()).digest():
                print(f"CSV payload differs from source: {member}", file=sys.stderr)
                return 1
        with tempfile.TemporaryDirectory(prefix="ufc-wheel-smoke-") as tmp:
            archive.extractall(tmp)
            env = os.environ.copy()
            env["PYTHONPATH"] = str(tmp)
            env["PYTHONNOUSERSITE"] = "1"
            code = (
                "from ufc_blast.core.blast_params import compute_facade; "
                "from ufc_blast.core import blast_params as bp; "
                "assert bp.__file__.startswith(" + repr(str(tmp)) + "); "
                "air=compute_facade(R=30,W=200,Hc=5,width=4,height=4,step=2); "
                "surface=compute_facade(R=30,W=200,Hc=0,width=4,height=4,step=2); "
                "assert air.burst_type=='air' and surface.burst_type=='surface'; "
                "assert len(air.grid_points)==len(surface.grid_points)==9; "
                "print('wheel smoke: air and surface, 9 points each')"
            )
            subprocess.run([sys.executable, "-c", code], cwd=tmp, env=env, check=True)
            subprocess.run(
                [sys.executable, "-m", "ufc_blast", "point", "--W", "200", "--R", "30", "--Hc", "5", "--alpha", "0"],
                cwd=tmp, env=env, check=True, stdout=subprocess.DEVNULL,
            )
    print(f"wheel CSV hashes match {len(expected)} normalized source tables; CLI and core smoke passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
