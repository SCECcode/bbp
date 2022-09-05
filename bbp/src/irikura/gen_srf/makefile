all: gen_srf gen_srf_segment sum_seg

gen_srf: gen_srf_4.9.3.f
	gfortran -o gen_srf gen_srf_4.9.3.f
	cp gen_srf ../bin/.

gen_srf_segment: gen_srf_segment4.f
	gfortran -o gen_srf_segment gen_srf_segment4.f
	cp gen_srf_segment ../bin/.

sum_seg: sum_seg.f
	gfortran -o sum_seg sum_seg.f
	cp sum_seg ../bin/.

clean:
	rm -rf *.o *~ gen_srf gen_srf_segment sum_seg
