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
 
hb_high_v5.4.3: hb_high_v5.4.3.f
	${FC} ${FFLAGS} -o hb_high_v5.4.3 hb_high_v5.4.3.f
	cp hb_high_v5.4.3 ../../bin/

clean:
	rm -f *.o hb_high_v5.4.3

%.c : ${MODLIB}
	touch $@
