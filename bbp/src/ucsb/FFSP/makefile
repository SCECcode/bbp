F90     = gfortran -O3
F77     = gfortran -O3

OBJ	      = ffsp_tool.o \
		ffsp_comm.o \
		slip_rate.o \
		spfield_n.o \
		dcf_subs_1.o\

all:		ffsp_dcf_v2 s2m install

ffsp_dcf_v2:     $(OBJ) ffsp_dcf_v2.f90
		$(F90) $(OBJ) -o $@ $@.f90 

s2m:		s2m.f90 s2m_comm.o
		$(F90) s2m.f90 s2m_comm.o -o $@

ffsp_tool.o:  ffsp_tool.f
		$(F77) -c ffsp_tool.f
s2m_comm.o:	s2m_comm.f90
		$(F90) -c s2m_comm.f90
ffsp_comm.o:  ffsp_comm.f90
		$(F90) -c ffsp_comm.f90
slip_rate.o:  slip_rate.f90
		$(F90) -c slip_rate.f90
spfield_n.o:  spfield_n.f90
		$(F90) -c spfield_n.f90
dcf_subs_1.o:  dcf_subs_1.f90
		$(F90) -c dcf_subs_1.f90

install:
	cp ffsp_dcf_v2 ../bin/.
	cp s2m ../bin/.

clean:		
	rm -f *.o *.mod s2m ffsp_dcf_v2

