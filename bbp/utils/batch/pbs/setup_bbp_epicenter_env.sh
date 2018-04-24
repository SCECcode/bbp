#!/bin/sh

#
# Configuration script to run Broadband validations on the epicenter cluster
#

# BBP_DIR should point to the top-level Broadband distribution directory
BBP_DIR=/home/epicenter/$USER/bbp/14.3.0;export BBP_DIR
# BBP_GF_DIR should point to the Green's Functions top-level directory
BBP_GF_DIR=/home/epicenter/$USER/bbp/bbp_gf;export BBP_GF_DIR
# BBP_VAL_DIR should point to the Validation packages' top-level directory
BBP_VAL_DIR=/home/epicenter/$USER/bbp/bbp_val;export BBP_VAL_DIR

# Do not change anything below this line
PYTHONPATH=$BBP_DIR/comps;export PYTHONPATH
BBP_DATA_DIR=/scratch/$PBS_JOBID/bbpruns;export BBP_DATA_DIR

# Create BBP_DATA_DIR
mkdir -p $BBP_DATA_DIR

if [ -n "$LD_LIBRARY_PATH" ]; then
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/epi/intel/composerxe-2011.4.191/compiler/lib/intel64; export LD_LIBRARY_PATH
else
    LD_LIBRARY_PATH=/usr/epi/intel/composerxe-2011.4.191/compiler/lib/intel64; export LD_LIBRARY_PATH
fi

# Set stack size
ulimit -s unlimited
