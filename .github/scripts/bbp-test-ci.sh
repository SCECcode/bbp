#!/bin/bash

# Set basic parameters
VERSION=22.4.0
BASEDIR="${RUNNER_WORKSPACE}"
BBPDIR="${BASEDIR}/bbp/bbp"
SRCDIR="${BBPDIR}/src"

export BBP_DIR=${BBPDIR}
export BBP_GF_DIR=${BASEDIR}/bbp_gf
export BBP_VAL_DIR=${BASEDIR}/bbp_val
export PYTHONPATH=${BBPDIR}/comps:${BBPDIR}/comps/PySeismoSoil:${PYTHONPATH}
export BBP_DATA_DIR=${BASEDIR}/bbp_data
export PATH=${BBPDIR}/comps:${BBPDIR}/utils/batch:$PATH
ulimit -s unlimited

echo
echo "===> Running Unit Tests..."

cd $BBP_DIR/tests
#./UnitTestsCI.py

./test_genslip.py

echo

./test_hfsims.py

echo

./test_rmg.py

echo

./test_irikura_hf.py

echo

./test_bbtoolbox.py

echo



