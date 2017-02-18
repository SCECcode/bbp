F90     = gfortran -O3
F77     = gfortran -O3

OBJ	      = ffsp_tool.o \
		ffsp_comm.o \
		slip_rate.o \
		spfield.o \

all:		ffsp_v2 install

ffsp_v2:     $(OBJ) ffsp_v2.f90
		$(F90) $(OBJ) -o $@ $@.f90 

ffsp_tool.o:  ffsp_tool.f
		$(F77) -c ffsp_tool.f
ffsp_comm.o:  ffsp_comm.f90
		$(F90) -c ffsp_comm.f90
slip_rate.o:  slip_rate.f90
		$(F90) -c slip_rate.f90
spfield.o:  spfield.f90
		$(F90) -c spfield.f90

install:
	cp ffsp_v2 ../bin/.

clean:		
	rm -f *.o *.mod ffsp_v2

