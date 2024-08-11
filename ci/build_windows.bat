@echo on
call conda info
if errorlevel 1 exit 1

call conda config --add channels conda-forge
if errorlevel 1 exit 1
call conda config --remove channels defaults || true
if errorlevel 1 exit 1
call conda config --show
if errorlevel 1 exit 1

call conda create ^
    --quiet --yes ^
    --name vigra ^
    python=%PYTHON_VERSION% pytest c-compiler cxx-compiler ^
    zlib libjpeg-turbo libpng libtiff hdf5 fftw cmake ninja ^
    boost=1.78 boost-cpp=1.78 numpy h5py sphinx ^
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
