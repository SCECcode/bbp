#!/bin/csh

if( $#argv < 2 ) then

echo "*** ERROR ***"
echo "  Need to specify input & output files: 'srf2vmod.csh infile=<srffile_in> outfile=<srffile_out>'"
echo "exiting ..."
exit

endif

set INFILE = `echo $argv | gawk -v par=infile '{x=-1;for(i=1;i<=NF;i++){n=split($i,a,"=");if(a[1]==par)x=a[2];}}END{print x;}'`

if( -e $INFILE == 0 ) then

echo "*** ERROR: Input file= '$INFILE' not found, exiting ..."
exit

endif

set OUTFILE = `echo $argv | gawk -v par=outfile '{x=-1;for(i=1;i<=NF;i++){n=split($i,a,"=");if(a[1]==par)x=a[2];}}END{print x;}'`

set RWG_BIN_DIR = /home/rgraves/Bin
set VELFILE = tmp_velfile.txt

set MY_UCVM_HOME = /home/rgraves/ucvm-22.7.0
set UCVM_MODEL = cvmsi

srf2xyz lonlatdep=1 infile=$INFILE | \
gawk -v conv=1000.0 'BEGIN{mindep=1.0;}{\
dep=$3*conv; \
if(dep<mindep)dep=mindep; \
printf "%.6f %.6f %.4f\n",$1,$2,dep;}' | \
${MY_UCVM_HOME}/bin/ucvm_query -f ${MY_UCVM_HOME}/conf/ucvm.conf -m ${UCVM_MODEL} | \
gawk -v vsmin=180 'BEGIN{np=-1;}{ \
vp=$7*100.0;vs=$8*100.0;den=$9/1000.0;if(vs<vsmin)vs=vsmin; \
printf "%.5e\t%.5e\t%.5e\n",vp,vs,den;}' \
- > ${VELFILE}

${RWG_BIN_DIR}/srf_adj-vmod velfile=${VELFILE} < $INFILE > $OUTFILE
