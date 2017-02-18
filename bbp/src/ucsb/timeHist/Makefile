FC=gfortran -O3
CC=gcc -O3

all:	toSAC BBPtoSAC copy

clean:
	rm -f *.o toSAC BBPtoSAC

toSAC:	toSAC.f sacio.o
	$(FC) -o toSAC toSAC.f sacio.o

BBPtoSAC:	BBPtoSAC.f sacio.o
	$(FC) -o BBPtoSAC BBPtoSAC.f sacio.o

sacio.o: sacio.c
	$(CC) -c sacio.c

copy:
	cp toSAC BBPtoSAC ../bin
