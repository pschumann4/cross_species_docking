@echo off

set /p folder_path=Enter Vina output file directory:

set /p chem=Chemical used for binding simulation:
set output_folder=%folder_path%\vina_output
mkdir "%output_folder%"
cd "%folder_path%"

for %%f in ("%folder_path%"\*_%chem%.pdbqt) do (
    echo Processing file: %%~nxf
    vina_split --input "%%f"
    setlocal enabledelayedexpansion
    for /l %%i in (1,1,100) do (
        if exist "%%~nf_flex_%%i.pdbqt" move "%%~nf_flex_%%i.pdbqt" "%output_folder%" > nul 2>&1
        if exist "%%~nf_rigid_%%i.pdbqt" move "%%~nf_rigid_%%i.pdbqt" "%output_folder%" > nul 2>&1
        if exist "%%~nf_ligand_%%i.pdbqt" move "%%~nf_ligand_%%i.pdbqt" "%output_folder%" > nul 2>&1
    )
    endlocal
)

echo Done.
pause
