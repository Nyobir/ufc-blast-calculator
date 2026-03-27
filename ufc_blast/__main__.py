"""Allow running as `python -m ufc_blast`."""

import sys


def _run() -> None:
    try:
        from ufc_blast.cli import main
        main()
    except Exception as exc:
        # When frozen as a GUI exe (console=False), stderr is invisible.
        # Show a native error dialog so the user sees what went wrong.
        if getattr(sys, "frozen", False):
            import traceback
            msg = traceback.format_exc()
            try:
                # Try Qt message box if PySide6 loaded
                from PySide6.QtWidgets import QApplication, QMessageBox
                app = QApplication.instance() or QApplication(sys.argv)
                QMessageBox.critical(None, "UFC Blast — Error", msg)
            except Exception:
                # Fallback to ctypes Win32 MessageBox
                try:
                    import ctypes
                    ctypes.windll.user32.MessageBoxW(0, msg, "UFC Blast — Error", 0x10)
                except Exception:
                    pass
            sys.exit(1)
        else:
            raise


_run()
