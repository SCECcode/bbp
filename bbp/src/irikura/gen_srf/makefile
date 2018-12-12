all: gen_srf

gen_srf: gen_srf_4.8.f
	gfortran -o gen_srf gen_srf_4.8.f
	cp gen_srf ../bin/.

clean:
	rm -rf *.o *~ gen_srf
