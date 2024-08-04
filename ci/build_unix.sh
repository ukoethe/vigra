set -ex

PYTHON_VERSION="${PYTHON_VERSION:-3.10}"

conda config --add channels conda-forge
conda config --remove channels defaults || true
conda config --show
conda create \
    --quiet --yes \
    --name vigra \
    python=${PYTHON_VERSION} pytest c-compiler cxx-compiler \
    zlib jpeg libpng libtiff hdf5 fftw \
    boost boost-cpp "numpy<2" h5py nose sphinx \
    openexr lemon cmake make

if [[ `uname` == 'Darwin' ]];
then
    export SHLIB_EXT=".dylib"
    export LDFLAGS="-undefined dynamic_lookup ${LDFLAGS}"
else
    export SHLIB_EXT=".so"
fi

source $CONDA/bin/activate vigra

mkdir build
cd build

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
    -DCMAKE_FIND_APPBUNDLE=LAST \
    -DZLIB_INCLUDE_DIR=${CONDA_PREFIX}/include \
    -DZLIB_LIBRARY=${CONDA_PREFIX}/lib/libz${SHLIB_EXT} \
\
    -DPNG_LIBRARY=${CONDA_PREFIX}/lib/libpng${SHLIB_EXT} \
    -DPNG_PNG_INCLUDE_DIR=${CONDA_PREFIX}/include


make -j2
make check -j2
ctest -V
