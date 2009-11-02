set PROG=%1
:extend_path
shift
if [%1] NEQ [] (
    set PATH=%1;%PATH%
    goto extend_path
)
%PROG%
echo off
if %ERRORLEVEL%==1 (
    del %PROG%
    exit 1 
)
