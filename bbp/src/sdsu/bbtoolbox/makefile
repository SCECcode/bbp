UFLAGS = -w -Wall
MODULES = module_bbtoolbox.f90 module_interface.f90
CODES = main_bbtoolbox.f90 coda.f90 composition.f90 convolution.f90 fourier.f90 error.f90 geometry.f90 interpolation.f90 io.f90 random.f90 scattering.f90 source.f90 correlation_readnk.f90 integ_diff.f90 filtering.f90

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

