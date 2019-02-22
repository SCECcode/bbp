HEADS = BasenSpectrum.h calcdistance.h complex.h csvread.h FastFourierTransform.h format.h knetascii.h knetout.h statgreen.h VerticalSHResponse.h
OBJS = BasenSpectrum.o calcdistance.o csvread.o FastFourierTransform.o knetascii.o knetout.o statgreen.o VerticalSHResponse.o

FLAGS = -O0 -g -Wall
G++ = g++

all: statgreen

statgreen: statgreen.cpp ${OBJS}
	${G++} ${FLAGS} -o statgreen ${OBJS}
	cp statgreen ../bin/.

${OBJS}: ${HEADS}

clean:
	rm -f ${OBJS} statgreen *~
