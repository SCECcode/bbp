#!/bin/bash

# For all other worker nodes, the environment variables are set by
# the setup_env.sh script, but for the main node running
# this script, we need to have these environment variables set up here
module load cpu/0.17.3b anaconda3/2021.05
module load gcc/10.2.0 
export BBP_DIR=/expanse/lustre/projects/usc143/qwxdev/apps/expanse/rocky8.8/bbp/bbp/bbp
export BBP_GF_DIR=/expanse/lustre/projects/usc143/qwxdev/apps/expanse/rocky8.8/bbp/bbp_gf
export BBP_VAL_DIR=/expanse/lustre/projects/usc143/qwxdev/apps/expanse/rocky8.8/bbp/bbp_val
export PYTHONPATH=/expanse/lustre/projects/usc143/qwxdev/apps/expanse/rocky8.8/bbp/bbp/bbp/comps:/expanse/lustre/projects/usc143/qwxdev/apps/expanse/rocky8.8/bbp/bbp/bbp/comps/PySeismoSoil:/expanse/lustre/projects/usc143/qwxdev/apps/expanse/rocky8.8/bbp/python3.8/site-packages:$PYTHONPATH
export PATH=/expanse/lustre/projects/usc143/qwxdev/apps/expanse/rocky8.8/bbp/bbp/bbp/comps:/expanse/lustre/projects/usc143/qwxdev/apps/expanse/rocky8.8/bbp/bbp/bbp/utils/batch:$PATH
ulimit -s unlimited

# BBP DATA setup
mkdir -p bbp_data
export BBP_DATA_DIR=`pwd`/bbp_data

# Quakeworx setup
output_dir=outputs/
mkdir -p $output_dir
HOME=`pwd`/bbp_sims

source ./tapisjob.env

######################
## INPUT PARAMETERS ##
######################

if [ "$BBP_HYPOCENTER_RANDOMIZATION" == "YES" ]
then
    HYPOCENTER_OPTION="--hypo-rand"
elif [ "$BBP_HYPOCENTER_RANDOMIZATION" == "NO" ]
then
    HYPOCENTER_OPTION="--no-hypo-rand"
else
    HYPOCENTER_OPTION=""
fi

if [ "$BBP_SITE_RESPONSE" == "GP" ]
then
    SITE_OPTION="--site gp"
elif [ "$BBP_SITE_RESPONSE" == "PySeismoSoil" ]
then
    SITE_OPTION="--site seismosoil"
elif [ "$BBP_SITE_RESPONSE" == "NO" ]
then
    SITE_OPTION=""
else
    echo "[ERROR]: Cannot determine BBP_SITE_RESPONSE, exiting..."
    exit
fi

if [ "$BBP_GENERATE_GMPE" == "1" ]
then
    if [ "$BBP_GMPE_MODEL" == "NGA-WEST1" ]
    then
	GMPE_OPTION="--gmpe-group NGA-West1"
    elif [ "$BBP_GMPE_MODEL" == "NGA-WEST2" ]
    then
	GMPE_OPTION="--gmpe-group NGA-West2"
    elif [ "$BBP_GMPE_MODEL" == "CENA" ]
    then
	GMPE_OPTION="--gmpe-group CENA\ GROUP\ 1"
    else
	echo "[ERROR]: Cannot determine BBP_GMPE_MODEL, exiting..."
	exit
    fi
else
    GMPE_OPTION=""
fi

## BBP PARAMETERS ##
if [ "$BBP_SIMULATION_MODE" == "VALIDATION" ]
then
    echo "[INFO]: Running Broadband Platform validation..."
    # Create the source files and xml workflows for each realization
    bbp_quakeworx_validation.py -d $HOME -e "${BBP_VALIDATION_EVENT}" -c ${BBP_SIMULATION_METHOD} -n ${BBP_NUM_REALIZATIONS} ${HYPOCENTER_OPTION} ${SITE_OPTION} ${GMPE_OPTION}
elif [ "$BBP_SIMULATION_MODE" == "SCENARIO" ]
then
    echo "[INFO]: Running Broadbabd Platform scenario..."
    # Validate input files
    mv ${BBP_SOURCE_DESCRIPTION} qwx-${BBP_SOURCE_DESCRIPTION}
    bbp_quakeworx_validate_src.py -i qwx-${BBP_SOURCE_DESCRIPTION} -o ${BBP_SOURCE_DESCRIPTION}
    if [[ $? -eq 0 ]]; then
	echo "SRC file parsing completed!"
    else
	echo "Failed to parse SRC file... exiting..."
	exit
    fi
    mv ${BBP_STATION_LIST} qwx-${BBP_STATION_LIST}
    bbp_quakeworx_validate_stl.py -i qwx-${BBP_STATION_LIST} -o ${BBP_STATION_LIST}
    if [[ $? -eq 0 ]]; then
	echo "Station file parsing completed!"
    else
	echo "Failed to parse station file... exiting..."
	exit
    fi
    # Create the source files and xml workflows for each realization
    bbp_quakeworx_scenario.py -d $HOME -c ${BBP_SIMULATION_METHOD} -n ${BBP_NUM_REALIZATIONS} ${SITE_OPTION} --station-list ${BBP_STATION_LIST} --source ${BBP_SOURCE_DESCRIPTION} -v ${BBP_VELOCITY_MODEL} ${HYPOCENTER_OPTION}
else
    echo "[ERROR]: Cannot determine BBP_SIMULATION mode, exiting..."
    exit
fi

echo "[INFO]: Broadband Platform simulation is ready to start..."

BBP_DATA_DIR=/expanse/lustre/scratch/$USER/temp_project/temp/bbp_sims/$SLURM_JOB_ID
mkdir -p $BBP_DATA_DIR
SLURM_NODES=`scontrol show hostname $SLURM_JOB_NODELIST | paste -d, -s`
mkdir -p $HOME/Sims/indata
mkdir -p $HOME/Sims/logs
mkdir -p $HOME/Sims/outdata
mkdir -p $HOME/Sims/tmpdata

echo "Processing start"
date
echo "BBP_DATA_DIR = $BBP_DATA_DIR"

cd $HOME

python3 $BBP_DIR/utils/batch/run_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh $HOME/Xml/batch_run_bbp_sims.log $SLURM_NODES 16
echo "Processing end"
date

# Copy BBP output from $BBP_DATA_DIR to $HOME/Sims
cp -frp $BBP_DATA_DIR/outdata/* $HOME/Sims/outdata/.
cp -frp $BBP_DATA_DIR/indata/* $HOME/Sims/indata/.
cp -frp $BBP_DATA_DIR/logs/* $HOME/Sims/logs/.

# Clean up temp results
rm -rf $BBP_DATA_DIR

echo "[INFO]: Broadband Platform run completed!"
exit
