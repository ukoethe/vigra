set -ex

PYTHON_VERSION="${PYTHON_VERSION:-3.10}"

export LDFLAGS="-undefined dynamic_lookup ${LDFLAGS}"

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

mkdir build
cd build

# Use CMAKE_FIND_FRAMEWORK=LAST to ensure that
# we don't find outdated OSX libraries on OSX
cmake .. \
    -DCMAKE_FIND_FRAMEWORK=LAST \
    -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
    -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_ROOT_DIR="${CONDA_PREFIX}" \
    -DPython_FIND_VIRTUALENV=ONLY \
    -DTEST_VIGRANUMPY=ON \
    -DWITH_OPENEXR=ON \
    -DWITH_LEMON=ON \
    -DAUTOEXEC_TESTS=OFF

make -j2
export DYLD_FALLBACK_LIBRARY_PATH="$PREFIX/lib":"${DYLD_FALLBACK_LIBRARY_PATH}"

make check -j2
ctest -V --output-on-failure
