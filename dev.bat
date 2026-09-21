@echo off
setlocal
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0dev.ps1" %*
set "clinreport_exit_code=%ERRORLEVEL%"
if not "%clinreport_exit_code%"=="0" pause
exit /b %clinreport_exit_code%
