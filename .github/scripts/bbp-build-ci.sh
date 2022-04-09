#!/bin/bash

echo "==> Installing Python packages needed by BBP..."

sudo apt install python3-numpy python3-scipy

echo $BBP_DIR


gcc --version

gfortran --version

python --version

python3 --version

python3 -c "import numpy; print(numpy.__version__)"

python3 -c "import scipy; print(scipy.__version__)"

python3 -c "import matplotlib; print(matplotlib.__version__)"

python3 -c "import numba; print(numba.__version__)"

python3 -c "import pyproj; print(pyproj.__version__)"

ls

