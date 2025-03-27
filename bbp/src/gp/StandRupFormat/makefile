HEADS = include.h structure.h function.h defs.h
OBJS = iofunc.o misc.o srf_subs.o sliprate_subs.o
GENRAND_OBJS = ../GenRand/ruptime.o
BAILEY_OBJS = ../JordanBailey/rob_rupm.o \
              ../JordanBailey/ruptime.o \
              ../JordanBailey/stf_subs.o
GEOPROJ_OBJS = ../ModelCords/geoproj_subs.o ../ModelCords/geo_utm.o

GETPAR = ../getpar/lib
INCPAR = -I ../getpar/include

# if Proj4 libraries ARE NOT available, use the following:
#
PROJ4_FLAGS = -D_NO_PROJ4
PROJ4_LIB_PATHS =
# if Proj4 libraries ARE available, use the following:
#
### PROJ4_FLAGS =
### PROJ4_LIB_PATHS = -lproj -L${PROJ4_LIBDIR} -I ${PROJ4_INCDIR}

LIBS = -lm ${GETPAR}/libget.a
LDLIBS = ${OBJS} ${LIBS} ${PROJ4_LIB_PATHS}

#LF_FLAGS = -D_FILE_OFFSET_BITS=32
#
# use following for large file capability
LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

FC = gfortran
CC = gcc

CFLAGS = ${UFLAGS} ${LF_FLAGS} ${PROJ4_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

all: srf2stoch generic_slip2srf fault_seg2gsf srf2moment srf2xyz srf2mrf srf_downsample

test_header : test_header.c ${OBJS} ${GEOPROJ_OBJS}
	$(CC) $(CFLAGS) -o test_header test_header.c ${LDLIBS} ${GEOPROJ_OBJS}

#srf2stoch : srf2stoch.c ${OBJS} ${GEOPROJ_OBJS}
#	${CC} ${CFLAGS} -o srf2stoch srf2stoch.c ${LDLIBS} ${INCPAR} ${GEOPROJ_OBJS}
#	cp srf2stoch ../bin/

srf2stoch : srf2stoch_sub.c srf2stoch_main.c ${OBJS} ${GEOPROJ_OBJS}
	${CC} ${CFLAGS} -c -o srf2stoch_sub.o srf2stoch_sub.c ${INCPAR}
	${CC} ${CFLAGS} -o srf2stoch srf2stoch_sub.o srf2stoch_main.c ${LDLIBS} ${INCPAR} ${GEOPROJ_OBJS}
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

srf2mrf : srf2mrf.c ${OBJS} ${GEOPROJ_OBJS}
	${CC} -o srf2mrf srf2mrf.c ${LDLIBS} ${INCPAR} ${GEOPROJ_OBJS}
	cp srf2mrf ../bin/

srf_downsample : srf_downsample.c ${OBJS} ${GEOPROJ_OBJS}
	${CC} -o srf_downsample srf_downsample.c ${LDLIBS} ${INCPAR} ${GEOPROJ_OBJS}
	cp srf_downsample ../bin/

${OBJS} : ${HEADS}

clean :
	rm -f *.o ${GENRAND_OBJS} ${BAILEY_OBJS} ${GEOPROJ_OBJS} *.o srf2xyz fault_seg2gsf srf2moment generic_slip2srf srf2stoch test_header srf2mrf srf_downsample
