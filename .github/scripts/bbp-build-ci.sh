#!/bin/bash

echo "==> Installing Python packages needed by BBP..."

sudo apt install python3-numpy python3-scipy python3-matplotlib python3-numba python3-pyproj

echo $BBP_BASE_DIR

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

echo $RUNNER_WORKSPACE

ls $RUNNER_WORKSPACE

