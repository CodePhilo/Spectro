"""Entry point used by PyInstaller and for running from a source checkout."""

from spectro.ui.main_window import main

if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()   # the optimizer's worker processes (frozen .exe)
    raise SystemExit(main())
