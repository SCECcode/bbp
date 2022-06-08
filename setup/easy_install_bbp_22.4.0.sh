#!/bin/bash

#
#
# BSD 3-Clause License
#
# Copyright (c) 2021, University of Southern California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# SCEC Broadband Platform Installation Script
#
#

# Set version number
VERSION=22.4.0

# Set base URL for downloads
BASEURL=https://g-c662a6.a78b8.36fe.data.globus.org/bbp/releases/${VERSION}

die () {
    echo >&2 "$@"
    exit 1
}

download_untar () {
    URL=$1
    MD5=$2

    # Download
    echo "===> Downloading..."
    #wget $URL 2>&1 | grep '%' || curl -O $URL
    curl -O ${URL}

    URL_NO_PROTO=${URL:7}
    URL_REL=${URL_NO_PROTO#*/}
    FILE=`basename "/${URL_REL%%\?*}"`

    # Untar
    echo "===> Uncompressing file..."
    tar -xzf ${FILE}

    # Check MD5SUM
    echo "===> Checking md5sum..."
    PREV_SUM=`cat ${MD5} | grep ${FILE} | cut -d\  -f1`
    SYSTEM="linux"
    md5sum --help > /dev/null 2>&1 || SYSTEM="mac"
    if [ ${SYSTEM} == "linux" ]; then
	CUR_SUM=`md5sum ${FILE} | cut -d\  -f1`
    else
	CUR_SUM=`md5 ${FILE} | cut -d\  -f4`
    fi

    if [ "${CUR_SUM}" != "${PREV_SUM}" ]; then
	echo
	echo "=> ERROR! MD5SUM check failed for ${FILE}"
	echo
	exit 1
    fi

    # Clean up as we go
    rm ${FILE}
}

# Get base directory to install all packages
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR=`echo ${DIR} | rev | cut -d'/' -f3- | rev`
BASEBBP=`echo ${DIR} | rev | cut -d'/' -f2- | rev`
BBPDIR="`echo ${DIR} | rev | cut -d'/' -f2- | rev`/bbp"
SRCDIR="${BBPDIR}/src"
MD5FILE="${BASEBBP}/setup/bbp-${VERSION}-md5.txt"

echo
echo " ====== Welcome to Broadband Platform ${VERSION} installation script ======"
echo
echo " Using destination directory: ${BASEDIR}"
echo

# Parse --force argument, if provided
if [ "$#" -eq 1 ]; then
    if [ $1 == "--force" ]; then
	while true; do
	    echo "=> WARNING!"
	    echo "==> Deleting existing bbp_val / bbp_gf / bbp_data directores!"
	    read -p "==> Do you wish to delete them and all their contents? " yn
	    case $yn in
		[Yy]* ) rm -rf ${BASEDIR}/bbp_val;
			rm -rf ${BASEDIR}/bbp_gf;
			rm -rf ${BASEDIR}/bbp_data;
			break;;
		[Nn]* ) echo "Aborting..."; exit;;
		* ) echo "Please answer yes or no.";;
	    esac
	done
    fi
fi

# Check if directories exist already
if [ -d "${BASEDIR}/bbp_val" ]; then
    echo "=> ${BASEDIR}/bbp_val directory exists!"
    echo "==> Use --force flag to delete bbp_val/bbp_gf/bbp_data or move them manually"
    exit 1
fi
if [ -d "${BASEDIR}/bbp_gf" ]; then
    echo "=> ${BASEDIR}/bbp_gf directory exists!"
    echo "==> Use --force flag to delete bbp_val/bbp_gf/bbp_data or move them manually"
    exit 1
fi
if [ -d "${BASEDIR}/bbp_data" ]; then
    echo "=> ${BASEDIR}/bbp_data directory exists!"
    echo "==> Use --force flag to delete bbp_val/bbp_gf/bbp_data or move them manually"
    exit 1
fi

# Create installation directories
echo "=> Creating directory tree..."
mkdir -p ${BASEDIR}/bbp_val
mkdir -p ${BASEDIR}/bbp_gf
mkdir -p ${BASEDIR}/bbp_data

# Compile source distribution
echo
echo " ====== Setting up Broadband Platform ${VERSION} Source Distribution ======"
echo
echo "==> Compiling... (it may take a few minutes)"
OLD_DIR=`pwd`
cd ${SRCDIR}
make > bbp-build.log 2>&1
make_ret_val=$?
# Check if build was successful
if [ $make_ret_val -ne 0 ]; then
    echo
    echo "****** ERROR: BBP build failed, for more details please check:"
    echo "****** ERROR: ${SRCDIR}/bbp-build.log"
    echo
    cd ${OLD_DIR}
    exit 1
fi

cd ${OLD_DIR}
# Done with main source distribution
echo "==> Build completed!"

# Install velocity model packages
echo
echo " ====== Setting up Broadband Platform ${VERSION} Simulation Regions  ======"
echo
echo " Please select what velocity models (regions) you would like to install,"
echo " using '1' for Yes, or '2' for No:"
echo

# Ask questions first
echo "==> Would you like to install the LA Basin region (25GB)"
echo "    REQUIRED for unit and acceptance tests to run?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) LABASIN=y; break;;
	No ) LABASIN=n; break;;
    esac
done

echo "==> Would you like to install the Northern California region (27GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) NOCAL=y; break;;
	No ) NOCAL=n; break;;
    esac
done

echo "==> Would you like to install the Central California region (14GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) CCAL=y; break;;
	No ) CCAL=n; break;;
    esac
done

echo "==> Would you like to install the Mojave region (27GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) MOJAVE=y; break;;
	No ) MOJAVE=n; break;;
    esac
done

echo "==> Would you like to install the Southern Sierra Nevada 2 region (22GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) SSN2=y; break;;
	No ) SSN2=n; break;;
    esac
done

echo "==> Would you like to install the Central Japan region (15GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) CJAPAN=y; break;;
	No ) CJAPAN=n; break;;
    esac
done

echo "==> Would you like to install the Western Japan region (15GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) WJAPAN=y; break;;
	No ) WJAPAN=n; break;;
    esac
done

echo
echo " ======    Installing Broadband Platform Velocity Model Packages     ======"
echo

if [ "${LABASIN}" == "y" ]; then
    cd ${BASEDIR}/bbp_gf
    echo "==> LA Basin"
    download_untar ${BASEURL}/labasin500-velocity-model-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${NOCAL}" == "y" ]; then
    cd ${BASEDIR}/bbp_gf
    echo "==> Northern California"
    download_untar ${BASEURL}/nocal500-velocity-model-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${CCAL}" == "y" ]; then
    cd ${BASEDIR}/bbp_gf
    echo "==> Central California"
    download_untar ${BASEURL}/centralcal500-velocity-model-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${MOJAVE}" == "y" ]; then
    cd ${BASEDIR}/bbp_gf
    echo "==> Mojave"
    download_untar ${BASEURL}/mojave500-velocity-model-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${SSN2}" == "y" ]; then
    cd ${BASEDIR}/bbp_gf
    echo "==> Southern Sierra Nevada 2"
    download_untar ${BASEURL}/ssn2500-velocity-model-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${CJAPAN}" == "y" ]; then
    cd ${BASEDIR}/bbp_gf
    echo "==> Central Japan"
    download_untar ${BASEURL}/centraljapan500-velocity-model-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${WJAPAN}" == "y" ]; then
    cd ${BASEDIR}/bbp_gf
    echo "==> Western Japan"
    download_untar ${BASEURL}/westernjapan500-velocity-model-${VERSION}.tar.gz ${MD5FILE}
fi

echo
echo "==> Completed!"

echo
echo " ======   Installing Broadband Platform Validation Event Packages    ======"
echo

if [ "${NOCAL}" == "y" ]; then
    cd ${BASEDIR}/bbp_val
    echo "==> LOMAP"
    download_untar ${BASEURL}/lomap-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Alum Rock"
    download_untar ${BASEURL}/alum-rock-validation-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${LABASIN}" == "y" ]; then
    cd ${BASEDIR}/bbp_val
    echo "==> NR"
    download_untar ${BASEURL}/nr-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Whittier Narrows"
    download_untar ${BASEURL}/whittier-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Chino Hills"
    download_untar ${BASEURL}/chino-hills-validation-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${CCAL}" == "y" ]; then
    cd ${BASEDIR}/bbp_val
    echo "==> Parkfield"
    download_untar ${BASEURL}/parkfield-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> San Simeon"
    download_untar ${BASEURL}/sansimeon-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> San Simeon Multi-Segment"
    download_untar ${BASEURL}/sansimeon-ms-validation-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${MOJAVE}" == "y" ]; then
    cd ${BASEDIR}/bbp_val
    echo "==> Landers"
    download_untar ${BASEURL}/landers-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Landers Multi-Segment"
    download_untar ${BASEURL}/landers-ms-draping-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> North Palm Springs"
    download_untar ${BASEURL}/nps-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Hector Mine"
    download_untar ${BASEURL}/hector-mine-validation-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${SSN2}" == "y" ]; then
    cd ${BASEDIR}/bbp_val
    echo "==> Ridgecrest 19A"
    download_untar ${BASEURL}/ridgecrest19a-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Ridgecrest 19B"
    download_untar ${BASEURL}/ridgecrest19b-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Ridgecrest 19C"
    download_untar ${BASEURL}/ridgecrest19c-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Ridgecrest 19C Multi-Segment"
    download_untar ${BASEURL}/ridgecrest19c-ms-validation-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${CJAPAN}" == "y" ]; then
    cd ${BASEDIR}/bbp_val
    echo "==> Chuetsu"
    download_untar ${BASEURL}/chuetsu-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Iwate"
    download_untar ${BASEURL}/iwate-validation-${VERSION}.tar.gz ${MD5FILE}
    echo "==> Niigata"
    download_untar ${BASEURL}/niigata-validation-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${WJAPAN}" == "y" ]; then
    cd ${BASEDIR}/bbp_val
    echo "==> Tottori"
    download_untar ${BASEURL}/tottori-validation-${VERSION}.tar.gz ${MD5FILE}
fi

if [ "${LABASIN}" == "y" ] || [ "${NOCAL}" == "y" ]; then
    cd ${BASEDIR}/bbp_val
    echo "==> GMPEs"
    download_untar ${BASEURL}/gmpe-verification-${VERSION}.tar.gz ${MD5FILE}
fi

echo "==> Completed!"

# All done!

echo
echo "=> All Done!"
echo
echo "Please add the following lines to your bash_profile:"
echo
echo "export BBP_DIR=${BBPDIR}"
echo "export BBP_GF_DIR=${BASEDIR}/bbp_gf"
echo "export BBP_VAL_DIR=${BASEDIR}/bbp_val"
echo "export PYTHONPATH=${BBPDIR}/comps:${BBPDIR}/comps/PySeismoSoil:${PYTHONPATH}"
echo "export BBP_DATA_DIR=${BASEDIR}/bbp_data"
echo "export PATH=${BBPDIR}/comps:${BBPDIR}/utils/batch:\$PATH"
if [ "$(uname)" != "Darwin" ]; then
    echo "ulimit -s unlimited"
fi
