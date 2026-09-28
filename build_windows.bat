@echo off
REM Builds dist\Spectro\Spectro.exe (run from the repository root on Windows).
python -m pip install --upgrade pip
python -m pip install -r requirements.txt pyinstaller
python -m PyInstaller --noconfirm --windowed --name Spectro ^
  --collect-submodules sklearn --collect-data sklearn ^
  --hidden-import scipy.special._cdflib ^
  run_spectro.py
echo.
echo Built: dist\Spectro\Spectro.exe
