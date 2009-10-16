set PROG=%1
:extend_path
shift
if NOT "%1" == "" (
    set PATH=%1;%PATH%
    goto extend_path
)
%PROG%
echo off
if ERRORLEVEL 1 (
    del %PROG%
    exit 1 
)
