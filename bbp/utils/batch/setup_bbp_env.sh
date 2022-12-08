#!/bin/sh

#
# Configuration script to run Broadband validations on USC's Discovery cluster
#

# BBP_DIR should point to the top-level Broadband distribution directory
BBP_DIR=/project/scec_608/$USER/22_4_DEV/22.4.0/bbp;export BBP_DIR
# BBP_GF_DIR should point to the Green's Functions top-level directory
BBP_GF_DIR=/scratch1/$USER/bbp_gf;export BBP_GF_DIR
# BBP_VAL_DIR should point to the Validation packages' top-level directory
BBP_VAL_DIR=/project/scec_608/$USER/22_4_DEV/bbp_val;export BBP_VAL_DIR
BBP_DATA_DIR=/scratch1/$USER/bbp_sims/$SLURM_JOB_ID;export BBP_DATA_DIR

# Do not change anything below this line
PYTHONPATH=$BBP_DIR/comps:$BBP_DIR/comps/PySeismoSoil;export PYTHONPATH

# Create BBP_DATA_DIR
mkdir -p $BBP_DATA_DIR

# Use right Python
module load gcc/8.3.0
module load python/3.7.6
#module load openmpi/4.0.2
#module load fftw/3.3.8-sp

# Set stack size
ulimit -s unlimited
