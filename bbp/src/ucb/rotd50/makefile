FC=gfortran
FFLAGS = -O2 -ffixed-line-length-none
HEADS = baseline.h
COMMON_OBJS = calcrsp.o fftsub.o ft_th.o sort.o spline.o splint.o
ROTD50_OBJS = ${COMMON_OBJS} rotd50.o
ROTD100_OBJS = ${COMMON_OBJS} rotd100.o

all: rotd50 rotd100

rotd50: ${ROTD50_OBJS}
	${FC} ${FFLAGS} -o rotd50 ${ROTD50_OBJS}

rotd100: ${ROTD100_OBJS}
	${FC} ${FFLAGS} -o rotd100 ${ROTD100_OBJS}

${OBJS}: ${HEADS}

clean:
	rm -f ${ROTD50_OBJS} ${ROTD100_OBJS} rotd50 rotd100 *~
