HEADS = BasenSpectrum.h calcdistance.h complex.h csvread.h FastFourierTransform.h format.h knetascii.h knetout.h statgreen.h VerticalSHResponse.h
OBJS = BasenSpectrum.o calcdistance.o csvread.o FastFourierTransform.o knetascii.o knetout.o statgreen.o VerticalSHResponse.o

FLAGS = -O0 -g -Wall
G++ = g++

all: statgreen2

statgreen2: statgreen.cpp ${OBJS}
	${G++} ${FLAGS} -o statgreen2 ${OBJS}
	cp statgreen2 ../bin/.

${OBJS}: ${HEADS}

clean:
	rm -f ${OBJS} statgreen2 *~
