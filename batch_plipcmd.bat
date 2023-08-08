@echo off
setlocal enabledelayedexpansion
REM plip must be installed and added to the PATH environment variable
set /p "DIR_PATH=Enter the path to the folder with the PDB files to analyze: "

cd /d "%DIR_PATH%"

for %%f in ("%DIR_PATH%\*.pdb") do (
  plip -f "%%f" -xv --name "%%~nf"
)

for %%f in ("%DIR_PATH%\*protonated.pdb") do (
    obabel "%%f" -o pdb -O "%DIR_PATH%\%%~nf.pdb" -h
)

REM Move plip results to a folder
if not exist "%DIR_PATH%\plip_results" mkdir "%DIR_PATH%\plip_results"
echo Moving files to the plip_results folder...
move "%DIR_PATH%\*.xml" "%DIR_PATH%\plip_results"
move "%DIR_PATH%\*_protonated.pdb" "%DIR_PATH%\plip_results"

pause
