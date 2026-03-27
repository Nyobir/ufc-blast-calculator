"""Allow running as `python -m ufc_blast`."""

import sys


def _run() -> None:
    try:
        from ufc_blast.cli import main
        main()
    except Exception:
        # When frozen as a GUI exe, stderr may be invisible.
        # Always write a crash log next to the exe first.
        if getattr(sys, "frozen", False):
            import traceback
            msg = traceback.format_exc()

            # 1. Write crash log next to the exe
            try:
                import os
                log_path = os.path.join(os.path.dirname(sys.executable), "ufc-blast-crash.log")
                with open(log_path, "w") as f:
                    f.write(msg)
            except Exception:
                pass

            # 2. Try native Win32 MessageBox (no Qt dependency)
            try:
                import ctypes
                ctypes.windll.user32.MessageBoxW(0, msg, "UFC Blast — Error", 0x10)
            except Exception:
                pass

            sys.exit(1)
        else:
            raise


_run()
