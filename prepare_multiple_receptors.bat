@echo off
setlocal enabledelayedexpansion

set "python_script=C:\Users\pschuman\OneDrive - Environmental Protection Agency (EPA)\Profile\Documents\CRADA-Molecular Modeling\Docking Analysis Functions\prep_receptor4.py"
set /p "DIR_PATH=Enter the path to the folder with the PDB files to prep: "
set output_folder=%DIR_PATH%\PDBQT_files
mkdir "%output_folder%"

cd /d "%DIR_PATH%"

REM Loop through all PDB files in the directory and convert them to PDBQT files
for %%f in ("%DIR_PATH%\*.pdb") do (
  python "!python_script!" -r "%%f" -v -o "%%~nf.pdbqt" -A "checkhydrogens"
  move "%%~nf.pdbqt" "%output_folder%" > nul 2>&1
)

REM Determine if there are matching PDBQT file names as the PDB, if not flag the files that were not converted
for %%f in ("%DIR_PATH%\*.pdb") do (
  set "pdb_name=%%~nf"
  set "pdbqt_name=!pdb_name!.pdbqt"
  if not exist "%output_folder%\!pdbqt_name!" (
    @echo on 
    echo ^^WARNING: !pdb_name! was not converted to PDBQT format. Consider a manual prep.
    @echo off
  )
)

pause