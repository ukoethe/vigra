set -ex

PYTHON_VERSION="${PYTHON_VERSION:-3.10}"

conda config --add channels conda-forge
conda config --remove channels defaults || true
conda config --show
conda create \
    --quiet --yes \
    --name vigra \
    python=${PYTHON_VERSION} pytest c-compiler cxx-compiler \
    zlib libjpeg-turbo libpng libtiff hdf5 fftw \
    boost boost-cpp numpy h5py sphinx \
    openexr lemon cmake make ruff

if [[ `uname` == 'Darwin' ]]; then
    export SHLIB_EXT=".dylib"
    # I - hmaarrrfk - forget why the definition of LDFLAGS here is necessary
    export LDFLAGS="-undefined dynamic_lookup ${LDFLAGS}"

    # clang14 introduced a new default that turns this flag on on osx arm64 leading
    # to failing tests. Setting it off restores clang13 behavior and seems in line
    # with x86 compilation.
    # ref: https://github.com/llvm/llvm-project/issues/91824
    export CXXFLAGS="-ffp-contract=off ${CXXFLAGS}"
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
