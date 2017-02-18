#! /bin/csh

set NT = 1000
set DT = 0.01

set STYPE = 2tri-p10-h50
set STYPE = 2tri-p20-h20
set STYPE = 2tri-p10-h20
set AMP = 1.0
set TZERO = 1.0
set TSHIFT = 0.5

set STYPE = miyatake
set TZERO = 4.0
set TZERO = 1.0
set TZERO = 4.01939
set TSHIFT = 0.5

set STYPE = 2tri-p10-h20
set TZERO = 1.0
set TSHIFT = 0.1

set OUTFILE = $STYPE$TZERO

wcc_genstf stype=$STYPE nt=$NT dt=$DT outfile=$OUTFILE \
	   npeak=1 t0=$TZERO ts=0.0 a0=$AMP tshift=$TSHIFT
