LDLIBS = -lm
UFLAGS = -O2

##### C options

# optimized
CFLAGS = ${UFLAGS} -Xt -fsingle

# debug
#CFLAGS = -g -fsingle

# optimized and profiled
#CFLAGS = -p -xO2 -Xt -fsingle

#CCOMPILE = cc ${CFLAGS} -o $@ $@.c ${LIBS}

FFLAGS = ${UFLAGS} -ffixed-line-length-132

FC = gfortran
#FC = g77

##### make options

all:  hb_high_v6.0.5

#hb_high_v6.0.5_orig : hb_high_v6.0.5.f
#	${FC} ${FFLAGS} -o hb_high_v6.0.5_orig hb_high_v6.0.5.f
#	cp -rp hb_high_v6.0.5_orig ../../bin/

hb_high_v6.0.5 : hb_high_v6.0.5_sub.o hb_high_v6.0.5_main.f
	${FC} ${FFLAGS} -o hb_high_v6.0.5 hb_high_v6.0.5_sub.o hb_high_v6.0.5_main.f
	cp -rp hb_high_v6.0.5 ../../bin/

# Generate object file for use in linking
hb_high_v6.0.5_sub.o : hb_high_v6.0.5_sub.f
	${FC} ${FFLAGS} -c -o hb_high_v6.0.5_sub.o hb_high_v6.0.5_sub.f

clean:
	rm -f *.o hb_high_v6.0.5
