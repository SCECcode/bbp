FC=gfortran -O3
CC=gcc -O3

all: noah_w copy

noah_w: noah_w.f sacio.o
	$(FC) -O3  -o $@ $@.f sacio.o

sacio.o: sacio.c
	$(CC) -c sacio.c

copy:
	mkdir -p ../bin;
	cp noah_w ../bin

clean:
	rm -f *.o noah_w
