@echo off

SET MyPath=%PATH%
rem echo %MyPath%
rem echo --

setlocal EnableDelayedExpansion

SET TempPath="%MyPath:;=";"%"
SET var=
FOR %%a IN (%TempPath%) DO (
    IF exist %%~sa (
        SET "var=!var!;%%~sa"
    ) ELSE (
        rem echo %%a does not exist
    )
)

rem echo --
echo !var:~1!
