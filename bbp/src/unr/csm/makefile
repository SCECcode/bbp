FC=gfortran
FFLAGS = -O3

all: green_99v8 simula

green_99v8:
	${FC} ${FFLAGS} -o green_99v8 B_green_99v8.f
	cp green_99v8 ../bin/.

simula:
	${FC} ${FFLAGS} -o simula C_simula_v12.f
	cp simula ../bin/.

clean:
	rm -f green_99v8 simula
