%1
echo off
if ERRORLEVEL 1 (
    del %1
    exit 1 
) else ( 
    exit 0 
)
