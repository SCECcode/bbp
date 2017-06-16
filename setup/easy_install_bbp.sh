#!/bin/sh

# Set version number
VERSION=17.3.0

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
    md5sum 2> /dev/null || SYSTEM="mac"
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
echo " ====== Welcome to Broadband Platform $VERSION installation script ====="
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
echo " Please select what velocity models (regions) you would like to install:"
echo

# Ask questions first
echo "==> Would you like to install the Northern California region (6.5GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) NOCAL=y; break;;
	No ) NOCAL=n; break;;
    esac
done

echo "==> Would you like to install the LA Basin region (6.5GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) LABASIN=y; break;;
	No ) LABASIN=n; break;;
    esac
done

echo "==> Would you like to install the Mojave region (6.5GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) MOJAVE=y; break;;
	No ) MOJAVE=n; break;;
    esac
done

echo "==> Would you like to install the Central Japan region (2.0GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) CJAPAN=y; break;;
	No ) CJAPAN=n; break;;
    esac
done

echo "==> Would you like to install the Western Japan region (2.0GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) WJAPAN=y; break;;
	No ) WJAPAN=n; break;;
    esac
done

echo "==> Would you like to install the Eastern Canada region (3.5GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) CANADA=y; break;;
	No ) CANADA=n; break;;
    esac
done

echo "==> Would you like to install the Eastern United States excluding Mississippi embayment and gulf coast region (3.8GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) CEUS=y; break;;
	No ) CEUS=n; break;;
    esac
done

echo "==> Would you like to install the Central United States including Mississippi embayment region (1.8GB)?"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) CENTRALUS=y; break;;
	No ) CENTRALUS=n; break;;
    esac
done

echo "=> Installing Broadband Platform Velocity Model Packages"

if [ "$NOCAL" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Northern California"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/nocal-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$LABASIN" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> LA Basin"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/labasin-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$MOJAVE" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Mojave"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/mojave-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$CJAPAN" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Central Japan"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/centraljapan-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$WJAPAN" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Western Japan"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/westernjapan-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$CANADA" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Canada"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/canada1000-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$CEUS" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Central Eastern United States"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/ceus1000-velocity-model-$VERSION.tar.gz $MD5FILE
fi

if [ "$CENTRALUS" == "y" ]; then
    cd $BASEDIR/bbp_gf
    echo "==> Central United States"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/centralus-velocity-model-$VERSION.tar.gz $MD5FILE
fi

echo "==> Completed!"

echo "=> Installing Broadband Platform Validation Packages"

if [ "$NOCAL" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> Loma Prieta"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/lomaprieta-validation-$VERSION.tar.gz $MD5FILE
fi

if [ "$LABASIN" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> Northridge"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/northridge-validation-$VERSION.tar.gz $MD5FILE
fi

if [ "$LABASIN" == "y" ] || [ "$NOCAL" == "y" ]; then
    cd $BASEDIR/bbp_val
    echo "==> GMPEs"
    download_untar http://hypocenter.usc.edu/research/bbp/versions/$VERSION/gmpe-verification-$VERSION.tar.gz $MD5FILE
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
echo "ulimit -s unlimited"
