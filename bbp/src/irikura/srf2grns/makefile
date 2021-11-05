INCLUDES =
HEADERS = params_3.1.f90 srf2grns_3.1.f90 readsrf_3.1.f90 resample_3.1.f90 sort3.f90
OBJS = params_3.1.o srf2grns_3.1.o readsrf_3.1.o resample_3.1.o sort3.o

CFLAGS = ${INCLUDES}
G++ = gfortran

all: srf2grns dtrandom

dtrandom: dtrandom.f90
	${G++} ${CFLAGS} -o dtrandom dtrandom.f90
	cp dtrandom ../bin/.

srf2grns: srf2grns_3.1.f90 ${OBJS}
	${G++} ${CFLAGS} -o srf2grns ${OBJS}
	cp srf2grns ../bin/.

${OBJS}: ${HEADERS}
	${G++} ${CFLAGS} -c -o $@ $*.f90

clean:
	rm -f ${OBJS} srf2grns dtrandom *.mod *~
