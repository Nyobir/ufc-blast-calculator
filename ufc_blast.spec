# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for UFC Blast Calculator."""

import sys
from pathlib import Path

block_cipher = None

ROOT = Path(SPECPATH)
DATA_DIR = ROOT / "data" / "free_air"

# Collect all CSV data files
data_files = [(str(f), str(Path("data") / "free_air")) for f in DATA_DIR.glob("*.csv")]

a = Analysis(
    [str(ROOT / "ufc_blast" / "__main__.py")],
    pathex=[str(ROOT)],
    binaries=[],
    datas=data_files,
    hiddenimports=[
        "ufc_blast",
        "ufc_blast.cli",
        "ufc_blast.gui",
        "ufc_blast.core",
        "ufc_blast.core.blast_params",
        "ufc_blast.core.geometry",
        "ufc_blast.core.interpolation",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["pyvista", "vtk", "tkinter"],
    noarchive=False,
    cipher=block_cipher,
)

pyz = PYZ(a.pure, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="ufc-blast",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,  # GUI app — no console window
    disable_windowed_traceback=False,
    argv_emulation=False,
    icon=None,
)
