@echo off
setlocal
cd /d "%~dp0"

python -c "import numpy, pandas, scipy, streamlit" >nul 2>&1
if errorlevel 1 (
  echo iDGRM dependencies are missing. Install them first with:
  echo.
  echo   python -m pip install -e ".[web]"
  echo.
  pause
  exit /b 1
)

python -m streamlit run app.py
endlocal
