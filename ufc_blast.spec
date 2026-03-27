# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for UFC Blast Calculator."""

import sys
from pathlib import Path
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

block_cipher = None

ROOT = Path(SPECPATH)
DATA_DIR = ROOT / "data" / "free_air"

# Collect all CSV data files
data_files = [(str(f), str(Path("data") / "free_air")) for f in DATA_DIR.glob("*.csv")]

# Include surface burst data if directory exists
SURFACE_DIR = ROOT / "data" / "surface"
if SURFACE_DIR.exists():
    data_files += [(str(f), str(Path("data") / "surface")) for f in SURFACE_DIR.glob("*.csv")]

# PySide6 Qt plugins (platforms, styles, imageformats) are needed at runtime
data_files += collect_data_files("PySide6", includes=["plugins/**/*"])

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
        # matplotlib Qt backend — imported at runtime, invisible to PyInstaller
        "matplotlib.backends.backend_qtagg",
        "matplotlib.backends.backend_agg",
    ] + collect_submodules("PySide6"),
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
    console=True,  # Show console for debugging — switch to False once stable
    disable_windowed_traceback=False,
    argv_emulation=False,
    icon=None,
)
