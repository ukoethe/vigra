@echo on
call conda info
if errorlevel 1 exit 1

call conda config --add channels conda-forge
if errorlevel 1 exit 1
call conda config --remove channels defaults || true
if errorlevel 1 exit 1
call conda config --show
if errorlevel 1 exit 1

rem cxx compiler version 1.5.0 -> vs2017
rem build currently doesn't work on vs2019
rem ref: https://github.com/ukoethe/vigra/issues/525
call conda create ^
    --quiet --yes ^
    --name vigra ^
    python=%PYTHON_VERSION% pytest c-compiler=1.5.0 cxx-compiler=1.5.0 ^
    zlib jpeg libpng libtiff hdf5 fftw cmake ninja ^
    boost=1.78 boost-cpp=1.78 "numpy<2" h5py sphinx ^
    openexr lemon

if errorlevel 1 exit 1

call activate vigra
if errorlevel 1 exit 1

mkdir build
if errorlevel 1 exit 1

cd build

cmake .. ^
    -G "NMake Makefiles" ^
    -DCMAKE_INSTALL_PREFIX=%CONDA_PREFIX% ^
    -DCMAKE_PREFIX_PATH=%CONDA_PREFIX% ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DPython_ROOT_DIR="%CONDA_PREFIX%" ^
    -DPython_FIND_VIRTUALENV=ONLY ^
    -DBUILD_SHARED_LIBS=ON ^
    -DTEST_VIGRANUMPY=ON ^
    -DWITH_OPENEXR=ON ^
    -DWITH_LEMON=ON

if errorlevel 1 exit 1

nmake
if errorlevel 1 exit 1

:REM nmake check
:REM if errorlevel 1 exit 1
:REM
:REM ctest -V
:REM if errorlevel 1 exit 1
