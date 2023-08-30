@echo off

REM Set the path to the directory containing the configuration files
set /p CONFIG_DIR=Enter the path to the configuration file directory:

REM Loop through all the configuration files in the config directory
for %%f in ("%CONFIG_DIR%"\*_conf.txt) do (
  cd "%CONFIG_DIR%"
  set base_name=%%~nf
  set base_name=%base_name:_conf=%
  @echo on
  echo Running AutoDock Vina on "%%f"...
  @echo off
  vina --config "%%f"
)

pause
