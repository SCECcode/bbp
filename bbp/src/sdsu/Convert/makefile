all : bin2bbp bbp2bin

FC=gfortran

clean :
	rm -f bin2bbp bbp2bin *.o

bin2bbp:
	${FC} -o bin2bbp bin2bbp.f90
	cp bin2bbp ../bin/

bbp2bin:
	${FC} -o bbp2bin bbp2bin.f90
	cp bbp2bin ../bin/
