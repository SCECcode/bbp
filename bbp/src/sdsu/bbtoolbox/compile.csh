#!/bin/csh

echo "Starting compiling..."

#compiles without linking the raytracer C-code
icc -w -c ray3DJHfor.c

#compiles and links the whole Fortran95/C codes

set CODE = 'main_bbtoolbox.f90 coda.f90 composition.f90 convolution.f90 fourier.f90 error.f90 geometry.f90 interpolation.f90 io.f90 random.f90 scattering.f90  source.f90 ray3DJHfor.o'
set MODULE = 'module_bbtoolbox.f90 module_interface.f90'

/opt/intel/Compiler/11.1/046/bin/intel64/ifort $MODULE $CODE -o BBtoolboxV15.exe

rm -f *.mod *.o

echo ''
