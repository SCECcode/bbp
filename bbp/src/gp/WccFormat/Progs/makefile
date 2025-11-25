COBJS = sacio.o iofunc.o fft1d.o
FOBJS = fourg.o mccamy.o zpass.o

ifdef FFTW_INCDIR
FFTW_INCFLAGS = -I ${FFTW_INCDIR}
endif
ifdef FFTW_LIBDIR
FFTW_LIBFLAGS = -L ${FFTW_LIBDIR}
endif

GETPAR = ../../getpar/lib
INCPAR = -I ../../getpar/include ${FFTW_INCFLAGS}

LIBS = -lm ${GETPAR}/libget.a
LDLIBS = ${COBJS} ${FOBJS} ${LIBS} ${FFTW_INCFLAGS} ${FFTW_LIBFLAGS} -lfftw3 -lfftw3f

#LF_FLAGS = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
LF_FLAGS = -D_GNU_SOURCE
UFLAGS = -O3

CC = gcc
FC = gfortran

##### compile options

CFLAGS = ${UFLAGS} ${LF_FLAGS} ${FFTW_INCFLAGS}
FFLAGS = ${UFLAGS} -ffixed-line-length-132

ALLOBJS = ${COBJS} ${FOBJS}

##### make options

all:wcc2bbp integ_diff wcc_getpeak wcc_tfilter wcc_siteamp09 wcc_siteamp14 wcc_add wcc_resamp_arbdt respect2bbp wcc_rotate

respect2bbp: respect2bbp.c ${COBJS} ${FOBJS}
	${CC} -o $@ $@.c ${INCPAR} ${LDLIBS}
	cp respect2bbp ../../bin/

wcc2bbp: wcc2bbp.c ${COBJS} ${FOBJS}
	${CC} -o $@ $@.c ${INCPAR} ${LDLIBS}
	cp wcc2bbp ../../bin/

#integ_diff: integ_diff.c ${COBJS} ${FOBJS}
#	${CC} -o $@ $@.c ${INCPAR} ${LDLIBS}
#	cp integ_diff ../../bin/

#wcc_getpeak: wcc_getpeak.c ${COBJS} ${FOBJS}
#	${CC} -o $@ $@.c ${INCPAR} ${LDLIBS}
#	cp wcc_getpeak ../../bin/

#wcc_siteamp14: wcc_siteamp14.c ${COBJS} ${FOBJS}
#	${CC} -o $@ $@.c ${INCPAR} ${LDLIBS}
#	cp wcc_siteamp14 ../../bin/

integ_diff: integ_diff_main.c integ_diff_sub.c ${COBJS} ${FOBJS}
	${CC} -c -o integ_diff_sub.o integ_diff_sub.c ${INCPAR}
	${CC} -o integ_diff integ_diff_main.c integ_diff_sub.o ${INCPAR} ${LDLIBS}
	cp integ_diff ../../bin/

wcc_getpeak: wcc_getpeak_sub.c wcc_getpeak_main.c ${COBJS} ${FOBJS}
	${CC} -c -o wcc_getpeak_sub.o wcc_getpeak_sub.c ${INCPAR}
	${CC} -o wcc_getpeak wcc_getpeak_sub.o wcc_getpeak_main.c ${INCPAR} ${LDLIBS}
	cp wcc_getpeak ../../bin/

wcc_siteamp14: wcc_siteamp14_sub.c wcc_siteamp14_main.c ${COBJS} ${FOBJS}
	${CC} -g -c -o wcc_siteamp14_sub.o wcc_siteamp14_sub.c ${INCPAR}
	${CC} -g -o wcc_siteamp14 wcc_siteamp14_sub.o wcc_siteamp14_main.c ${INCPAR} ${LDLIBS}
	cp wcc_siteamp14 ../../bin/

wcc_siteamp09: wcc_siteamp09.c ${COBJS} ${FOBJS}
	${CC} -o wcc_siteamp09 wcc_siteamp09.c ${INCPAR} ${LDLIBS}
	cp wcc_siteamp09 ../../bin/

#wcc_add: wcc_add.c ${COBJS} ${FOBJS}
#	${CC} -o $@ $@.c ${INCPAR} ${LDLIBS}
#	cp wcc_add ../../bin/

wcc_add: wcc_add_sub.c wcc_add_main.c ${COBJS} ${FOBJS}
	${CC} ${CFLAGS} -c -o wcc_add_sub.o wcc_add_sub.c ${INCPAR}
	${CC} ${CFLAGS} -o wcc_add wcc_add_sub.o wcc_add_main.c ${INCPAR} ${LDLIBS}
	cp wcc_add ../../bin/

wcc_siteamp: wcc_siteamp.c ${COBJS} ${FOBJS}
	${CC} -o $@ $@.c ${LDLIBS}
	cp wcc_siteamp ../../bin/

wcc_Xcor: wcc_Xcor.c ${COBJS} ${FOBJS}
	${CC} -o $@ $@.c ${LDLIBS}
	cp wcc_Xcor ../../bin/

wcc_resamp_arbdt: wcc_resamp_arbdt_sub.c wcc_resamp_arbdt_main.c ${COBJS} ${FOBJS}
	${CC} -c -o wcc_resamp_arbdt_sub.o wcc_resamp_arbdt_sub.c ${INCPAR}
	${CC} -o wcc_resamp_arbdt wcc_resamp_arbdt_sub.o wcc_resamp_arbdt_main.c ${INCPAR} ${LDLIBS}
	cp wcc_resamp_arbdt ../../bin/

#wcc_resamp_arbdt: wcc_resamp_arbdt.c ${COBJS} ${FOBJS}
#	${CC} -o $@ $@.c ${INCPAR} ${LDLIBS}
#	cp wcc_resamp_arbdt ../../bin/

#wcc_tfilter: wcc_tfilter.c ${COBJS} ${FOBJS}
#	${CC} -o $@ $@.c ${INCPAR} ${LDLIBS}
#	cp wcc_tfilter ../../bin/

wcc_tfilter: wcc_tfilter_sub.c wcc_tfilter_main.c ${COBJS} ${FOBJS}
	${CC} -c -o wcc_tfilter_sub.o wcc_tfilter_sub.c ${INCPAR}
	${CC} -o wcc_tfilter wcc_tfilter_sub.o wcc_tfilter_main.c ${INCPAR} ${LDLIBS}
	cp wcc_tfilter ../../bin/

wcc_rotate: wcc_rotate.c ${COBJS} ${FOBJS}
	${CC} -o wcc_rotate wcc_rotate.c ${INCPAR} ${LDLIBS}
	cp wcc_rotate ../../bin/

clean:
	rm -f *.o wcc2bbp integ_diff wcc_getpeak wcc_tfilter wcc_siteamp09 wcc_siteamp14 wcc_add wcc_resamp_arbdt respect2bbp wcc_rotate
