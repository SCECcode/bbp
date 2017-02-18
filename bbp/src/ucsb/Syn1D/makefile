CC = gcc
F90 = gfortran 

all: syn_1d install

syn_1d: syn_1d.f90 syn1d_comm.o sacio.o
	${F90} -O3 -o $@ $@.f90 sacio.o syn1d_comm.o
syn1d_comm.o: syn1d_comm.f90
	${F90} -O3 -c syn1d_comm.f90
sacio.o: sacio.c
	$(CC) -c sacio.c

install:
	cp syn_1d ../bin/.

clean:
	rm -f *.o *.mod syn_1d


