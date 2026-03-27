# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for UFC Blast Calculator — directory bundle."""

import sys
from pathlib import Path
from PyInstaller.utils.hooks import collect_data_files

block_cipher = None

ROOT = Path(SPECPATH)
DATA_DIR = ROOT / "data" / "free_air"

# Collect all CSV data files
data_files = [(str(f), str(Path("data") / "free_air")) for f in DATA_DIR.glob("*.csv")]

# Include surface burst data if directory exists
SURFACE_DIR = ROOT / "data" / "surface"
if SURFACE_DIR.exists():
    data_files += [(str(f), str(Path("data") / "surface")) for f in SURFACE_DIR.glob("*.csv")]

# Only bundle the Qt plugins we actually need (platforms, styles, imageformats)
# — skip 3D asset importers, SQL drivers, multimedia, etc.
data_files += collect_data_files("PySide6", includes=[
    "plugins/platforms/*",
    "plugins/styles/*",
    "plugins/imageformats/*",
])

# Unused Qt modules — exclude to cut ~500 MB of DLLs
_qt_excludes = [
    "PySide6.Qt3DAnimation", "PySide6.Qt3DCore", "PySide6.Qt3DExtras",
    "PySide6.Qt3DInput", "PySide6.Qt3DLogic", "PySide6.Qt3DRender",
    "PySide6.QtBluetooth",
    "PySide6.QtCharts",
    "PySide6.QtConcurrent",
    "PySide6.QtDataVisualization",
    "PySide6.QtDBus",
    "PySide6.QtDesigner",
    "PySide6.QtGraphs", "PySide6.QtGraphsWidgets",
    "PySide6.QtHttpServer",
    "PySide6.QtLocation",
    "PySide6.QtMultimedia", "PySide6.QtMultimediaWidgets",
    "PySide6.QtNfc",
    "PySide6.QtNetworkAuth",
    "PySide6.QtPdf", "PySide6.QtPdfWidgets",
    "PySide6.QtPositioning",
    "PySide6.QtQml", "PySide6.QtQmlModels", "PySide6.QtQmlCore",
    "PySide6.QtQuick", "PySide6.QtQuick3D", "PySide6.QtQuickWidgets",
    "PySide6.QtQuickControls2",
    "PySide6.QtRemoteObjects",
    "PySide6.QtScxml",
    "PySide6.QtSensors",
    "PySide6.QtSerialBus", "PySide6.QtSerialPort",
    "PySide6.QtSpatialAudio",
    "PySide6.QtSql",
    "PySide6.QtStateMachine",
    "PySide6.QtSvgWidgets",
    "PySide6.QtTest",
    "PySide6.QtTextToSpeech",
    "PySide6.QtVirtualKeyboard",
    "PySide6.QtWebChannel",
    "PySide6.QtWebEngine", "PySide6.QtWebEngineCore",
    "PySide6.QtWebEngineWidgets",
    "PySide6.QtWebSockets",
    "PySide6.QtXml",
]

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
        # Only the PySide6 modules we actually use
        "PySide6.QtCore",
        "PySide6.QtGui",
        "PySide6.QtWidgets",
        "PySide6.QtOpenGL",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["pyvista", "vtk", "tkinter"] + _qt_excludes,
    noarchive=False,
    cipher=block_cipher,
)

# Strip unused Qt DLLs/resources from the binary list
_strip_prefixes = (
    "Qt6Quick", "Qt6Qml", "Qt6Multimedia", "Qt63D", "Qt6Designer",
    "Qt6Pdf", "Qt6WebEngine", "Qt6Charts", "Qt6DataVis", "Qt6Graphs",
    "Qt6Bluetooth", "Qt6Nfc", "Qt6Sensors", "Qt6Serial", "Qt6Sql",
    "Qt6Test", "Qt6VirtualKeyboard", "Qt6StateMachine", "Qt6HttpServer",
    "Qt6Spatial", "Qt6Scxml", "Qt6RemoteObjects", "Qt6TextToSpeech",
    "Qt6Location", "Qt6Positioning", "Qt6Concurrent",
    "avcodec", "avformat", "avutil", "swresample", "swscale",
    "opengl32sw",
)
a.binaries = [
    b for b in a.binaries
    if not any(Path(b[0]).name.startswith(p) for p in _strip_prefixes)
]

# Also strip qml/ resources/ translations/ directories
a.datas = [
    d for d in a.datas
    if not any(seg in d[0] for seg in ("qml/", "qml\\", "translations/", "translations\\", "resources/", "resources\\"))
]

pyz = PYZ(a.pure, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="ufc-blast",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    console=False,  # GUI app — no console window
    disable_windowed_traceback=False,
    argv_emulation=False,
    icon=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="ufc-blast",
)
