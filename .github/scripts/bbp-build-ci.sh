#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

download_untar () {
    URL=$1

    # Download
    wget $URL 2>/dev/null || curl -O $URL 2>/dev/null

    URL_NO_PROTO=${URL:7}
    URL_REL=${URL_NO_PROTO#*/}
    FILE=`basename "/${URL_REL%%\?*}"`

    # Untar
    tar -xzf ${FILE}

    # Clean up as we go
    rm ${FILE}
}

echo

OLD_DIR=`pwd`
mkdir ${RUNNER_WORKSPACE}/bin
cd ${RUNNER_WORKSPACE}/bin
ln -s /usr/bin/gcc-8 gcc
ln -s /usr/bin/gfortran-8 gfortran
cd ${OLD_DIR}
export PATH=${RUNNER_WORKSPACE}/bin:$PATH

echo "======================GCC===================="
gcc --version

echo "===================GFORTRAN=================="
gfortran --version

echo "===================Python 3=================="
python3 --version
python3 -c "import numpy; print('Numpy: ', numpy.__version__)"
python3 -c "import scipy; print('SciPy: ', scipy.__version__)"
python3 -c "import matplotlib; print('Matplotlib: ', matplotlib.__version__)"
python3 -c "import numba; print('Numba: ', numba.__version__)"
python3 -c "import pyproj; print('PyProj: ', pyproj.__version__)"

# Set basic parameters
VERSION=22.4.0
BASEDIR="${RUNNER_WORKSPACE}"
BBPDIR="${BASEDIR}/bbp/bbp"
SRCDIR="${BBPDIR}/src"
FFTW_BUILD_DIR="${BASEDIR}/fftbuild"
FFTW_INSTALL_DIR="${BASEDIR}/fftw-3.3.8"
export FFTW_INCDIR=${FFTW_INSTALL_DIR}/include
export FFTW_LIBDIR=${FFTW_INSTALL_DIR}/lib

echo "===================FFTW=================="
# Compile FFTW library
mkdir -p ${FFTW_BUILD_DIR}
OLD_DIR=`pwd`
cd ${FFTW_BUILD_DIR}
download_untar https://g-c662a6.a78b8.36fe.data.globus.org/bbp/releases/${VERSION}/fftw-3.3.8.tar.gz
cd fftw-3.3.8
./configure --prefix=${FFTW_INSTALL_DIR} --enable-float
make
make install
cd ${OLD_DIR}

# Create installation directories
echo "=> Creating directory tree..."
mkdir -p ${BASEDIR}/bbp_val
mkdir -p ${BASEDIR}/bbp_gf
mkdir -p ${BASEDIR}/bbp_data
echo

# Compile source distribution
echo "=> Main Broadband Platform Source Distribution"
echo "==> Compiling... (it may take a while)"
OLD_DIR=`pwd`
cd ${SRCDIR}
make
cd ${OLD_DIR}
# Done with main source distribution
echo "==> Source code compiled!"
echo

# Install LABasin500 (CI) region and NR validation packages
echo "==> Installing LA Basin CI region..."
cd ${BASEDIR}/bbp_gf
download_untar https://g-c662a6.a78b8.36fe.data.globus.org/bbp/releases/${VERSION}/labasin500ci-velocity-model-${VERSION}.tar.gz
echo
ls ${BASEDIR}/bbp_gf
echo

echo "==> Installing NR validation package..."
cd ${BASEDIR}/bbp_val
download_untar https://g-c662a6.a78b8.36fe.data.globus.org/bbp/releases/${VERSION}/nr-validation-${VERSION}.tar.gz
echo
ls ${BASEDIR}/bbp_val
echo

echo "==> Build steps completed!"
