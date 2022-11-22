@echo on
conda info
if errorlevel 1 exit 1

conda config --add channels conda-forge
if errorlevel 1 exit 1
conda config --remove channels defaults || true
if errorlevel 1 exit 1
conda config --show
if errorlevel 1 exit 1

conda create ^
    --quiet --yes ^
    --name vigra ^
    python=%PYTHON_VERSION% c-compiler cxx-compiler ^
    zlib jpeg libpng libtiff hdf5 fftw ^
    boost boost-cpp numpy h5py nose sphinx ^
    openexr lemon
if errorlevel 1 exit 1

call activate vigra
if errorlevel 1 exit 1

mkdir build
if errorlevel 1 exit 1

cd build

cmake .. ^
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
