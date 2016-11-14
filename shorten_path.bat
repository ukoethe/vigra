@echo off

SET MyPath=%PATH%
echo %MyPath%
echo --

setlocal EnableDelayedExpansion

SET TempPath="%MyPath:;=";"%"
SET var=
FOR %%a IN (%TempPath%) DO (
    IF exist %%~sa (
        SET "var=!var!;%%~sa"
    ) ELSE (
        echo %%a does not exist
    )
)

echo --
echo !var:~1!
