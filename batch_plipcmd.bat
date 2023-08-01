@echo off
setlocal enabledelayedexpansion

set "python_script=C:\Users\pschuman\OneDrive - Environmental Protection Agency (EPA)\Profile\Documents\CRADA-Molecular Modeling\Docking Analysis Functions\plipcmd.py"
set /p "DIR_PATH=Enter the path to the folder with the PDB files to analyze: "

cd /d "%DIR_PATH%"

for %%f in ("%DIR_PATH%\*.pdb") do (
  python "!python_script!" -f "%%f" -xv --name "%%~nf"
)

for %%f in ("%DIR_PATH%\*protonated.pdb") do (
    obabel "%%f" -o pdb -O "%DIR_PATH%\%%~nf.pdb" -h
)

pause
