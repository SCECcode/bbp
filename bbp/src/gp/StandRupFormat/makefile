HEADS = include.h structure.h function.h defs.h
OBJS = iofunc.o misc.o srf_subs.o sliprate_subs.o
GENRAND_OBJS = ../GenRand/ruptime.o
BAILEY_OBJS = ../JordanBailey/rob_rupm.o \
              ../JordanBailey/ruptime.o \
              ../JordanBailey/stf_subs.o
GEOPROJ_OBJS = ../ModelCords/geoproj_subs.o ../ModelCords/geo_utm.o

GETPAR = ../getpar/lib
INCPAR = -I ../getpar/include

LIBS = -lm ${GETPAR}/libget.a
LDLIBS = ${OBJS} ${LIBS}

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

FC = gfortran
CC = gcc

CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

all: srf2stoch generic_slip2srf fault_seg2gsf srf2moment srf2xyz

test_header : test_header.c ${OBJS} ${GEOPROJ_OBJS}
	$(CC) $(CFLAGS) -o test_header test_header.c ${LDLIBS} ${GEOPROJ_OBJS}

srf2stoch : srf2stoch.c ${OBJS} ${GEOPROJ_OBJS}
	${CC} ${CFLAGS} -o srf2stoch srf2stoch.c ${LDLIBS} ${INCPAR} ${GEOPROJ_OBJS}
	cp srf2stoch ../bin/

generic_slip2srf : generic_slip2srf.c ${OBJS} ${BAILEY_OBJS} ${GEOPROJ_OBJS}
	${CC} -o generic_slip2srf generic_slip2srf.c ${LDLIBS} ${INCPAR} ${BAILEY_OBJS} ${GEOPROJ_OBJS}
	cp generic_slip2srf ../bin/

fault_seg2gsf : fault_seg2gsf.c ${OBJS} ${GEOPROJ_OBJS}
	${CC} ${CFLAGS} -o fault_seg2gsf fault_seg2gsf.c ${LDLIBS} ${INCPAR} ${GEOPROJ_OBJS}
	cp fault_seg2gsf ../bin/

srf2moment : srf2moment.c ${OBJS} ${BAILEY_OBJS} ${GEOPROJ_OBJS}
	${CC} -o srf2moment srf2moment.c ${LDLIBS} ${BAILEY_OBJS} ${INCPAR} ${GEOPROJ_OBJS}
	cp srf2moment ../bin/

srf2xyz : srf2xyz.c ${OBJS} ${GEOPROJ_OBJS}
	${CC} -o srf2xyz srf2xyz.c ${LDLIBS} ${INCPAR} ${GEOPROJ_OBJS}
	cp srf2xyz ../bin/

${OBJS} : ${HEADS}

clean :
	rm -f *.o ${GENRAND_OBJS} ${BAILEY_OBJS} ${GEOPROJ_OBJS} *.o srf2xyz fault_seg2gsf srf2moment generic_slip2srf srf2stoch test_header
