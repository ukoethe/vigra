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
    openexr lemon cmake make

export LDFLAGS="-undefined dynamic_lookup ${LDFLAGS}"

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

make -j2
make check -j2
ctest -V
