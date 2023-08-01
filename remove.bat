@echo off

REM Create a function that takes a directory path as an argument and removes "m1_" and "+2ama+" from all the *.pdb files in that directory

REM Usage: remove.bat <directory path>

set pdb_dir=%1

for %%f in (%pdb_dir%\*.pdb) do (
    set pdb_file=%%f
    set pdb_file=%pdb_file:m1_=%
    set pdb_file=%pdb_file:+2ama+=%
    move %%f %pdb_file%
)

