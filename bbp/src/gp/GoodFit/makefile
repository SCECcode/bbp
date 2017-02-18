GETPAR = ../getpar/lib

LDLIBS = -lm ${GETPAR}/libget.a
UFLAGS = -O3
INCPAR = -I ../getpar/include

CC = gcc
FC = gfortran

CFLAGS = ${UFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

##### make options

all: resid2uncer_varN respect gen_resid_tbl gen_resid_tbl_3comp

resid2uncer_varN:
	$(CC) -o resid2uncer_varN resid2uncer_varN.c ${INCPAR} ${LDLIBS}
	cp resid2uncer_varN ../bin/ 

gen_resid_tbl:
	$(CC) $(UFLAGS) gen_resid_tbl.c ${LDLIBS} ${INCPAR} -o gen_resid_tbl
	cp gen_resid_tbl ../bin/ 

gen_resid_tbl_3comp:
	$(CC) $(UFLAGS) gen_resid_tbl_3comp.c ${LDLIBS} ${INCPAR} -o gen_resid_tbl_3comp
	cp gen_resid_tbl_3comp ../bin/

respect: respect.o pseudo.o
	$(FC) -o respect respect.o pseudo.o ${LDLIBS} ${INCPAR}
	cp respect ../bin/

clean:
	rm -f *.o respect resid2uncer_varN gen_resid_tbl gen_resid_tbl_3comp
