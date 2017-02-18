##  If you want to change the size of GF bank, please just modify the size.h
FFLAGS=-O
LINK.f = ifort -O3
FC=ifort -O3 -132
CC=icc -O3
GCC=gcc -O3
SUBS = fft.o Complex.o 
FKSUBS = green_bank.o cmodel.o fk.o kernel.o compound.o haskell.o source.o bessel.o $(SUBS)
FKSUBSF = green_bankFar.o cmodel.o fkFar.o kernel.o compound.o haskell.o source.o bessel.o $(SUBS)
SACLIB = /local/lib
all: green_bank green_bankFar getDuration copy

copy:
	cp green_bank green_bankFar ../bin/

green_bank: \
	$(FKSUBS)
	$(FC) -o green_bank $(FKSUBS)

green_bankFar: \
	$(FKSUBSF)
	$(FC) -o green_bankFar $(FKSUBSF)

fft.o:\
	fft.c 
	$(CC) -c fft.c

Complex.o:\
	Complex.c Complex.h
	$(CC) -c Complex.c

cmodel.o:\
	cmodel.f
	$(FC) -c  cmodel.f

fk.o:\
	fk.f
	$(FC) -c fk.f

fkFar.o:\
	fkFar.f
	$(FC) -c fkFar.f

kernel.o:\
	kernel.f
	$(FC) -c kernel.f

compound.o:\
	compound.f
	$(FC) -c compound.f

haskell.o:\
	haskell.f
	$(FC) -c haskell.f

source.o:\
	source.f
	$(FC) -c source.f
bessel.o:\
	bessel.f
	$(FC) -c bessel.f

green_bank.o: \
	green_bank.f
	$(FC) -c green_bank.f

green_bankFar.o: \
	green_bankFar.f
	$(FC) -c green_bankFar.f

getDuration:	getDuration.c utlDuration.o
		$(GCC) -o getDuration getDuration.c utlDuration.o -lm

utlDuration.o:	utlDuration.c utlDuration.h
		$(GCC) -c utlDuration.c

clean:
	rm -f *.o green_bankFar green_bank getDuration
