all: gen_srf

gen_srf: gen_srf_4.7.f
	gfortran -o gen_srf gen_srf_4.7.f
	cp gen_srf ../bin/.

clean:
	rm -rf *.o *~ gen_srf
