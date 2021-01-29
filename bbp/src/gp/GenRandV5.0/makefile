HEADS = include.h structure.h function.h defs.h
OBJS = fourg.o iofunc.o misc.o slip.o ruptime.o gslip_srf_subs.o gslip_sliprate_subs.o geoproj_subs.o
SRF_OBJS = ../StandRupFormat/srf_subs.o
FDRT_OBJS = wafront2d-rwg.o

GETPAR = ../getpar/lib
INCPAR = -I ../getpar/include

ifdef FFTW_INCDIR
FFTW_INCFLAGS = -I ${FFTW_INCDIR}
endif
ifdef FFTW_LIBDIR
FFTW_LIBFLAGS = -L ${FFTW_LIBDIR}
endif

LIBS = -lm ${GETPAR}/libget.a
LDLIBS = ${OBJS} ${LIBS} ${FFTW_INCFLAGS} ${FFTW_LIBFLAGS} -lfftw3f
#LDLIBS = ${OBJS} ${LIBS} ${FFTW_INCFLAGS} ${FFTW_LIBFLAGS} -lfftw3 -lfftw3f

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
#LF_FLAGS = 
UFLAGS = -O3

CC = gcc
FC = gfortran

CFLAGS = ${UFLAGS} ${LF_FLAGS} ${FFTW_INCFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-none

##### make options

all : genslip-v5.4.1 srf_gethypo get_rvfrac

genslip-v5.4.1 : genslip-v5.4.1.c ${OBJS} ${SRF_OBJS} ${FDRT_OBJS}
	${CC} ${CFLAGS} -o genslip-v5.4.1 genslip-v5.4.1.c ${LDLIBS} ${INCPAR} ${SRF_OBJS} ${FDRT_OBJS}
	cp genslip-v5.4.1 ../bin/.

genslip-v5.2.2 : genslip-v5.2.2.c ${OBJS} ${SRF_OBJS} ${FDRT_OBJS}
	${CC} ${CFLAGS} -o genslip-v5.2.2 genslip-v5.2.2.c ${LDLIBS} ${INCPAR} ${SRF_OBJS} ${FDRT_OBJS}
	cp genslip-v5.2.2 ../bin/.

srf_gethypo : srf_gethypo.c
	${CC} ${CFLAGS} -o srf_gethypo srf_gethypo.c ${LDLIBS} ${INCPAR} ${SRF_OBJS}
	cp srf_gethypo ../bin/.

get_rvfrac : get_rvfrac.c ${OBJS}
	${CC} ${CFLAGS} -o get_rvfrac get_rvfrac.c ${LDLIBS} ${INCPAR} ${SRF_OBJS}

${OBJS} : ${HEADS}

clean :
	-rm -f *.o genslip-v5.4.1 srf_gethypo get_rvfrac
