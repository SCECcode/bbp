#!/bin/bash

# Set basic parameters
VERSION=22.4.0
BASEDIR="${RUNNER_WORKSPACE}"
BBPDIR="${BASEDIR}/bbp/bbp"
SRCDIR="${BBPDIR}/src"

BBP_DIR=${BBPDIR}
BBP_GF_DIR=${BASEDIR}/bbp_gf
BBP_VAL_DIR=${BASEDIR}/bbp_val
PYTHONPATH=${BBPDIR}/comps:${PYTHONPATH}
BBP_DATA_DIR=${BASEDIR}/bbp_data
PATH=${BBPDIR}/comps:${BBPDIR}/utils/batch:$PATH

echo
echo "===> Running Unit Tests..."

cd $BBP_DIR/tests
./UnitTestsCI.py


