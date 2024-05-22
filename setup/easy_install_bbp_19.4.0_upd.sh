#!/bin/sh

# Set version number
VERSION=19.4.0

BASEURL=https://g-c662a6.a78b8.36fe.data.globus.org/bbp/releases/${VERSION}

die () {
    echo >&2 "$@"
    exit 1
}

download_untar () {
    URL=$1
    MD5=$2

    # Download
    wget $URL 2>/dev/null || curl -O $URL 2>/dev/null

    URL_NO_PROTO=${URL:7}
    URL_REL=${URL_NO_PROTO#*/}
    FILE=`basename "/${URL_REL%%\?*}"`

    # Untar
    tar -xzf $FILE

    # Check MD5SUM
    PREV_SUM=`cat $MD5 | grep $FILE | cut -d\  -f1`
    SYSTEM="linux"
    md5sum --help > /dev/null 2>&1 || SYSTEM="mac"
    if [ $SYSTEM == "linux" ]; then
	CUR_SUM=`md5sum $FILE | cut -d\  -f1`
    else
	CUR_SUM=`md5 $FILE | cut -d\  -f4`
    fi

    if [ "$CUR_SUM" != "$PREV_SUM" ]; then
	echo
	echo "=> ERROR! MD5SUM check failed for $FILE"
	echo
	exit 1
    fi

    # Clean up as we go
    rm $FILE
}

# Get base directory to install all packages
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR=`echo $DIR | rev | cut -d'/' -f3- | rev`
BASEBBP=`echo $DIR | rev | cut -d'/' -f2- | rev`
BBPDIR="`echo $DIR | rev | cut -d'/' -f2- | rev`/bbp"
SRCDIR="$BBPDIR/src"
MD5FILE="$BASEBBP/setup/bbp-$VERSION-md5.txt"

echo
echo " ====== Welcome to Broadband Platform $VERSION installation script ======"
echo
echo " Using destination directory: $BASEDIR"
echo

# Parse --force argument, if provided
if [ "$#" -eq 1 ]; then
    if [ $1 == "--force" ]; then
	while true; do
	    echo "=> WARNING!"
	    echo "==> Deleting existing bbp_val / bbp_gf / bbp_data directores!"
	    read -p "==> Do you wish to delete them and all their contents? " yn
	    case $yn in
		[Yy]* ) rm -rf $BASEDIR/bbp_val;
			rm -rf $BASEDIR/bbp_gf;
			rm -rf $BASEDIR/bbp_data;
			break;;
		[Nn]* ) echo "Aborting..."; exit;;
		* ) echo "Please answer yes or no.";;
	    esac
	done
    fi
fi

# Check if directories exist already
if [ -d "$BASEDIR/bbp_val" ]; then
    echo "=> $BASEDIR/bbp_val directory exists!"
    echo "==> Use --force flag to delete bbp_val/bbp_gf/bbp_data or move them manually"
    exit 1
fi
if [ -d "$BASEDIR/bbp_gf" ]; then
    echo "=> $BASEDIR/bbp_gf directory exists!"
    echo "==> Use --force flag to delete bbp_val/bbp_gf/bbp_data or move them manually"
    exit 1
fi
if [ -d "$BASEDIR/bbp_data" ]; then
    echo "=> $BASEDIR/bbp_data directory exists!"
    echo "==> Use --force flag to delete bbp_val/bbp_gf/bbp_data or move them manually"
    exit 1
fi

# Create installation directories
echo "=> Creating directory tree..."
mkdir -p $BASEDIR/bbp_val
mkdir -p $BASEDIR/bbp_gf
mkdir -p $BASEDIR/bbp_data
echo

# Compile source distribution
echo "=> Main Broadband Platform Source Distribution"
echo "==> Compiling... (it may take a while)"
OLD_DIR=`pwd`
cd $SRCDIR
make > /dev/null 2>&1
cd $OLD_DIR
# Done with main source distribution
echo "==> Installed!"

# Install velocity model packages
echo
echo " Please select what velocity models (regions) you would like to install,"
echo " using '1' for Yes, or '2' for No:"
echo

# Ask questions first
echo "==> Would you like to install the Northern California region (27GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) NOCAL=y; break;;
	No ) NOCAL=n; break;;
    esac
done

echo "==> Would you like to install the LA Basin region (25GB)"
echo "    needed for unit and acceptance tests to run?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) LABASIN=y; break;;
	No ) LABASIN=n; break;;
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

echo "==> Would you like to install the Central Japan region (28GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) CJAPAN=y; break;;
	No ) CJAPAN=n; break;;
    esac
done

echo "==> Would you like to install the Western Japan region (28GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) WJAPAN=y; break;;
	No ) WJAPAN=n; break;;
    esac
done

echo "=> Installing Broadband Platform Velocity Model Packages"

if [ "$NOCAL" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Northern California"
    download_untar ${BASEURL}/nocal500-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$LABASIN" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> LA Basin"
    download_untar ${BASEURL}/labasin500-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$CCAL" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Central California"
    download_untar ${BASEURL}/centralcal500-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$MOJAVE" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Mojave"
    download_untar ${BASEURL}/mojave500-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$CJAPAN" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Central Japan"
    download_untar ${BASEURL}/centraljapan500-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$WJAPAN" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Western Japan"
    download_untar ${BASEURL}/westernjapan500-velocity-model-$VERSION.tar.gz $MD5FILE
fi

echo "==> Completed!"

echo "=> Installing Broadband Platform Validation Packages"

if [ "$NOCAL" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> LOMAP"
    download_untar ${BASEURL}/lomap-validation-$VERSION.tar.gz $MD5FILE
    echo "==> Alum Rock"
    download_untar ${BASEURL}/alum-rock-validation-$VERSION.tar.gz $MD5FILE
fi

if [ "$LABASIN" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> NR"
    download_untar ${BASEURL}/nr-validation-$VERSION.tar.gz $MD5FILE
    echo "==> Whittier Narrows"
    download_untar ${BASEURL}/whittier-validation-$VERSION.tar.gz $MD5FILE
    echo "==> Chino Hills"
    download_untar ${BASEURL}/chino-hills-validation-$VERSION.tar.gz $MD5FILE
fi

if [ "$CCAL" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> Parkfield"
    download_untar ${BASEURL}/parkfield-validation-$VERSION.tar.gz $MD5FILE
    echo "==> San Simeon"
    download_untar ${BASEURL}/sansimeon-validation-$VERSION.tar.gz $MD5FILE
fi

if [ "$MOJAVE" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> Landers"
    download_untar ${BASEURL}/landers-validation-$VERSION.tar.gz $MD5FILE
    echo "==> North Palm Springs"
    download_untar ${BASEURL}/nps-validation-$VERSION.tar.gz $MD5FILE
fi

if [ "$CJAPAN" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> Chuetsu"
    download_untar ${BASEURL}/chuetsu-validation-$VERSION.tar.gz $MD5FILE
    echo "==> Iwate"
    download_untar ${BASEURL}/iwate-validation-$VERSION.tar.gz $MD5FILE
    echo "==> Niigata"
    download_untar ${BASEURL}/niigata-validation-$VERSION.tar.gz $MD5FILE
fi

if [ "$WJAPAN" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> Tottori"
    download_untar ${BASEURL}/tottori-validation-$VERSION.tar.gz $MD5FILE
fi

if [ "$LABASIN" == "y" ] || [ "$NOCAL" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> GMPEs"
    download_untar ${BASEURL}/gmpe-verification-$VERSION.tar.gz $MD5FILE
fi

echo "==> Completed!"

# All done!

echo
echo "=> All Done!"
echo
echo "Please add the following lines to your bash_profile:"
echo
echo "export BBP_DIR=$BBPDIR"
echo "export BBP_GF_DIR=$BASEDIR/bbp_gf"
echo "export BBP_VAL_DIR=$BASEDIR/bbp_val"
echo "export PYTHONPATH=$BBPDIR/comps"
echo "export BBP_DATA_DIR=$BASEDIR/bbp_data"
echo "export PATH=$BBPDIR/comps:$BBPDIR/utils/batch:\$PATH"
if [ "$(uname)" != "Darwin" ]; then
    echo "ulimit -s unlimited"
fi
