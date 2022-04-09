#!/bin/bash

echo "==> Installing Python packages needed by BBP..."

sudo apt install python3-numpy python3-scipy python3-matplotlib python3-numba python3-pyproj libfftw3-dev libfftw3-doc

echo

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
echo "==> Installed!"

# Install LABasin500 (CI) region and NR validation packages
echo "==> LA Basin"
cd ${BASEDIR}/bbp_gf

# download_untar http://hypocenter.usc.edu/research/bbp/versions/${VERSION}/labasin500-velocity-model-${VERSION}.tar.gz ${MD5FILE}

echo "==> NR"
cd ${BASEDIR}/bbp_val

# download_untar http://hypocenter.usc.edu/research/bbp/versions/${VERSION}/nr-validation-${VERSION}.tar.gz ${MD5FILE}

echo "==> Completed!"
