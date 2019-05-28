INCLUDES = 
HEADS = params_2.1.f90 srf2grns_2.1.f90 readsrf_2.1.f90 sdrop_2.1.f90 resample_2.1.f90 sort3.f90
OBJS = params_2.1.o srf2grns_2.1.o readsrf_2.1.o sdrop_2.1.o resample_2.1.o sort3.o

CFLAGS = ${INCLUDES}
G++ = gfortran

all: srf2grns

srf2grns: srf2grns_2.1.f90 ${OBJS}
	${G++} ${CFLAGS} -o srf2grns ${OBJS}
	cp srf2grns ../bin/.

${OBJS}: ${HEADS}
	${G++} ${CFLAGS} -c -o $@ $*.f90

clean:
	rm -f ${OBJS} srf2grns *.mod *~
