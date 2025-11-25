#!/bin/sh

#
# Configuration script to run Broadband validations on USC's HPCC cluster
#

# BBP_DIR should point to the top-level Broadband distribution directory
BBP_DIR=/home/rcf-104/$USER/bbp/15.3.0/bbp;export BBP_DIR
# BBP_GF_DIR should point to the Green's Functions top-level directory
BBP_GF_DIR=/home/rcf-104/$USER/bbp/bbp_gf;export BBP_GF_DIR
# BBP_VAL_DIR should point to the Validation packages' top-level directory
BBP_VAL_DIR=/home/rcf-104/$USER/bbp/bbp_val;export BBP_VAL_DIR

# Do not change anything below this line
python_base_dir=/usr/usc/python/2.7.8; export python_base_dir
PYTHONPATH=$BBP_DIR/comps:/home/scec-00/opt/python/2.7.8/lib/python2.7/site-packages;export PYTHONPATH
BBP_DATA_DIR=$TMPDIR/bbpruns;export BBP_DATA_DIR

# Create BBP_DATA_DIR
mkdir -p $BBP_DATA_DIR

# Use right Python
_bindir=$python_base_dir/bin
_libdir=$python_base_dir/lib

if [ -n "$PATH" ]; then
   PATH=$_bindir:$PATH;export PATH
else
   PATH=$_bindir;export PATH
fi

if [ -n "$LD_LIBRARY_PATH" ]; then
   LD_LIBRARY_PATH=$_libdir:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
else
   LD_LIBRARY_PATH=$_libdir; export LD_LIBRARY_PATH
fi

# Set stack size
ulimit -s unlimited
