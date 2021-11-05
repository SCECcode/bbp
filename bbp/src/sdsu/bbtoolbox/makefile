UFLAGS = -w -Wall
MODULES = module_bbtoolbox_0328_lenupdate0421_corr_rupspeed.f90 module_interface_corr.f90
CODES = main_bbtoolbox_0315_corr.f90 coda.f90 composition_lenupdate0421_corr_rupspeed.f90 convolution.f90 fourier.f90 error.f90 geometry.f90 interpolation.f90 io.f90 random.f90 scattering.f90 source.f90 correlation_lenupdate_corr.f90 integ_diff.f90 filtering.f90

all : BBtoolbox.exe

CC=gcc
FC=gfortran

ray3DJHfor.o : ray3DJHfor.c
	${CC} ${UFLAGS} -O0 -c ray3DJHfor.c

BBtoolbox.exe : ray3DJHfor.o ${CODES} ${MODULES}
	${FC} ${MODULES} ${CODES} -O0 ray3DJHfor.o -o BBtoolbox.exe
	cp BBtoolbox.exe ../bin/

clean :
	rm -f *.mod *.o BBtoolbox.exe

