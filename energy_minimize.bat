@echo off

set /p folder_path=Enter PDB file directory:
cd "%folder_path%"

@echo on

for %%f in ("%folder_path%\*.pdb") do (
    echo Converting "%%f" to mol2 format...
    obabel "%%f" -o mol2 -O "%folder_path%\%%~nf.mol2"
)

for %%f in ("%folder_path%\*.mol2") do (
    echo Energy minimizing "%%f"...
    obminimize -ff GAFF -n 250 -sd -c 1e-6 "%%f"
    obabel "%%f" -o pdb -O "%folder_path%\%%~nf.pdb"
)

pause
