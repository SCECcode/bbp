OBJS = geoproj_subs.o geo_utm.o

GETPAR = ../getpar/lib
INCPAR = -I ../getpar/include

LIBS = -lm ${GETPAR}/libget.a

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

UFLAGS = -O3

CC = gcc
GFORTRAN = gfortran

##### compile options

CFLAGS = ${UFLAGS} ${LF_FLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

all: xy2ll ll2xy

xy2ll : xy2ll.c ${OBJS}
	${CC} ${CFLAGS} -o xy2ll xy2ll.c ${INCPAR} ${OBJS} ${LIBS}
	cp xy2ll ../bin/

ll2xy: ll2xy.c ${OBJS}
	${CC} ${CFLAGS} -o ll2xy ll2xy.c ${INCPAR} ${OBJS} ${LIBS}
	cp ll2xy ../bin/

geoproj_subs.o: geoproj_subs.c
	${CC} -o geoproj_subs.o ${UFLAGS} -c geoproj_subs.c

geo_utm.o: geo_utm.f
	${GFORTRAN} -o geo_utm.o ${FFLAGS} -c geo_utm.f

clean:
	rm -f *.o xy2ll ll2xy
