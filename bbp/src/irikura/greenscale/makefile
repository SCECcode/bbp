HEADS = calcdistance.h filter.h greenscale.h knetout.h sekibun.h FastFourierTransform.h format.h knetascii.h saitofilt.h
OBJS = calcdistance.o FastFourierTransform.o filter.o gosei.o greenintensity.o greenscale.o knetascii.o knetout.o saitofilt.o scaleparam.o sekibun.o

FLAGS =
G++ = g++

all: greenscale

greenscale: greenscale.cpp ${OBJS}
	${G++} ${FLAGS} -o greenscale ${OBJS}
	cp greenscale ../bin/.

${OBJS}: ${HEADS}

clean:
	rm -f ${OBJS} greenscale *~
