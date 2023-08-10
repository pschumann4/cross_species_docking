@echo off

set vina_split_exe="C:\Users\pschuman\OneDrive - Environmental Protection Agency (EPA)\Profile\Documents\CRADA-Molecular Modeling\vina_split.exe"
set /p folder_path=Enter Vina output file directory:

set /p chem=Chemical used for binding simulation: 

for %%f in ("%folder_path%"\*_%chem%.pdbqt) do (
    cd "%folder_path%"
    echo Processing file: %%~nxf
    %vina_split_exe% --input "%%f"
)

echo Done.
pause