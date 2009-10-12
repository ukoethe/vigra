set PATH=..\..\src\impex;%PATH%
IF EXIST "%2" del "%2"
%1
echo off
if ERRORLEVEL 1 (
    del %1
    exit 1 
) else ( 
    echo success > "%2"
    exit 0 
)
