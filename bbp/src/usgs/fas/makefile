all: asc2smc smc2fs2

asc2smc: asc2smc.for
	gfortran -o asc2smc asc2smc.for
	cp asc2smc ../bin/

smc2fs2: smc2fs2.for
	gfortran -o smc2fs2 smc2fs2.for
	cp smc2fs2 ../bin

clean:
	rm -rf *.o *~ asc2smc smc2fs2
