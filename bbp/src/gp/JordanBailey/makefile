HEADS = include.h structure.h function.h
OBJS = iofunc.o ruptime.o misc_subs.o greenfunc.o \
       beroza_rupm.o okumura_rupm.o gene_rupm.o rob_rupm.o srf_rupm.o \
       fourg.o stf_subs.o sgt3d_subs.o geoproj_subs.o

SRF_OBJS = ../StandRupFormat/srf_subs.o

LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

GETPAR = ../getpar/lib
LDLIBS = -lm ${GETPAR}/libget.a
INCPAR = -I ../getpar/include

UFLAGS = -O
LDFLAGS = -p
LDFLAGS =

CFLAGS = -O ${LF_FLAGS} ${LDFLAGS}
FFLAGS = ${UFLAGS}

CC = gcc
FC = gfortran

##### make options

all : ray_stimes jbsim-v2.0.0

jbsim032906 : jbsim032906.o ${OBJS} ${SRF_OBJS}
	${CC} ${CFLAGS} -o jbsim032906 jbsim032906.o ${OBJS} ${INCPAR} ${LDLIBS} ${SRF_OBJS}

ray_stimes : ray_stimes.c ${OBJS}
	${CC} ${CFLAGS} -o ray_stimes ray_stimes.c ${OBJS} ${INCPAR} ${LDLIBS}

jbsim3d : jbsim3d.c ${OBJS} ${SRF_OBJS}
	${CC} ${CFLAGS} -o jbsim3d jbsim3d.c ${OBJS} ${INCPAR} ${LDLIBS} ${SRF_OBJS}
	cp jbsim3d ../bin/

jbsim : jbsim.c ${OBJS} ${SRF_OBJS}
	${CC} ${CFLAGS} -o jbsim jbsim.c ${OBJS} ${INCPAR} ${LDLIBS} ${SRF_OBJS}
	cp jbsim ../bin/

jbsim-v2.0.0 : jbsim-v2.0.0.c ${OBJS} ${SRF_OBJS}
	${CC} ${CFLAGS} -o jbsim-v2.0.0 jbsim-v2.0.0.c ${OBJS} ${INCPAR} ${LDLIBS} ${SRF_OBJS}
	cp jbsim-v2.0.0 ../bin/

gen_gflist : gen_gflist.o ${OBJS} ${SRF_OBJS}
	${CC} ${CFLAGS} -o gen_gflist gen_gflist.o ${OBJS} ${INCPAR} ${LDLIBS} ${SRF_OBJS}

bailey2srf : bailey2srf.o ${OBJS}
	${CC} ${CFLAGS} -o bailey2srf bailey2srf.o ${OBJS} ${INCPAR} ${LDLIBS}

${OBJS} : ${HEADS}

clean :
	-rm -f *.o jbsim jbsim3d ray_stimes jbsim-v2.0.0
