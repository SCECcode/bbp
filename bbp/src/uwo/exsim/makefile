all: exsim14

exsim14: EXSIM14_SCEC.for
	gfortran -O0 -o exsim14 EXSIM14_SCEC.for
	cp exsim14 ../bin/

clean:
	rm -rf *.o exsim14

