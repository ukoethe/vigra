set -ex

PYTHON_VERSION="${PYTHON_VERSION:-3.10}"

# As of Aug 31, 2025, Python 3.14 is only available in the conda-forge python_rc channel
# We just choose to unconditionally add it here for simplicity of testing

conda config --show
conda create \
    --quiet --yes \
    --override-channels \
    -c conda-forge/label/python_rc \
    -c conda-forge \
    -c nodefaults \
    --name vigra \
    python=${PYTHON_VERSION} pytest c-compiler cxx-compiler \
    zlib libjpeg-turbo libpng libtiff hdf5 fftw \
    libboost-python libboost-python-devel numpy h5py sphinx \
    openexr lemon cmake make ruff

if [[ `uname` == 'Darwin' ]]; then
    export SHLIB_EXT=".dylib"
    # I - hmaarrrfk - forget why the definition of LDFLAGS here is necessary
    export LDFLAGS="-undefined dynamic_lookup ${LDFLAGS}"
else
    export SHLIB_EXT=".so"
fi

# Set the variable CPU_COUNT to take on the default value of 2
export CPU_COUNT=${CPU_COUNT:-2}

source $CONDA/bin/activate vigra

# lint quickly to help find obvious mistakes
( cd vigranumpy && ruff check . )

mkdir -p build
pushd build

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


make -j${CPU_COUNT}
make check -j${CPU_COUNT}
ctest -V
popd
