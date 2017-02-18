FC=gfortran -O3

all: stitch stitchBBP deconv deconvBBP copy

stitch: stitch.f90 proc_tool.o
	$(FC) -o $@ $@.f90 proc_tool.o

stitchBBP: stitchBBP.f90 proc_tool.o
	$(FC) -o $@ $@.f90 proc_tool.o

deconv: deconv.f90 proc_tool.o
	$(FC) -o $@ $@.f90 proc_tool.o

deconvBBP: deconvBBP.f90 proc_tool.o
	$(FC) -o $@ $@.f90 proc_tool.o

proc_tool.o: proc_tool.f90
	$(FC) -c proc_tool.f90  

copy:
	mkdir -p ../bin;
	cp deconvBBP stitchBBP stitch ../bin 

clean:
	rm -f *.o *.mod stitchBBP deconvBBP stitch deconv
