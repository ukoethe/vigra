set -ex

PYTHON_VERSION="${PYTHON_VERSION:-3.10}"


conda config --add channels conda-forge
conda config --remove channels defaults || true
conda config --show
conda create \
    --quiet --yes \
    --name vigra \
    python=${PYTHON_VERSION} c-compiler cxx-compiler \
    zlib jpeg libpng libtiff hdf5 fftw \
    boost boost-cpp numpy h5py nose sphinx \
    openexr lemon

source $CONDA/bin/activate vigra


if [ ! -d "${CONDA_PREFIX}/include/libpng" ]; then
    # CMake FindPNG seems to look in libpng not libpng16
    # https://gitlab.kitware.com/cmake/cmake/blob/master/Modules/FindPNG.cmake#L55
    ln -s ${CONDA_PREFIX}/include/libpng16 ${CONDA_PREFIX}/include/libpng
fi

mkdir build
cd build

# Use CMAKE_FIND_FRAMEWORK=LAST to ensure that
# we don't find outdated OSX libraries on OSX
cmake .. \
    -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
    -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_ROOT_DIR="${CONDA_PREFIX}" \
    -DPython_FIND_VIRTUALENV=ONLY \
    -DTEST_VIGRANUMPY=ON \
    -DWITH_OPENEXR=ON \
    -DWITH_LEMON=ON \
    -DAUTOEXEC_TESTS=OFF \
    -DCMAKE_FIND_FRAMEWORK=LAST \
    -DPNG_PNG_INCLUDE_DIR=${CONDA_PREFIX}/include

make -j2
make check -j2
ctest -V
