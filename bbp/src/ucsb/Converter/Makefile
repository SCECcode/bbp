CC=gcc -O3

all:	srfLL2XYKM srfXY2LLKM statLL2XY srfLL2XYnf convertFaultInput conv3Comp statLLXY conv3Comp2 conv3CompBB convertFaultInputKM copy

clean:
	rm -f srfLL2XYKM srfXY2LLKM statLL2XY srfLL2XYnf convertFaultInput conv3Comp statLLXY conv3Comp2 conv3CompBB convertFaultInputKM *.o

srfLL2XYKM:	srfLL2XYKM.c geoproj_subs.o utlFault.o
	$(CC) -o srfLL2XYKM srfLL2XYKM.c geoproj_subs.o utlFault.o -lm

srfLL2XYnf:	srfLL2XYnf.c geoproj_subs.o utlFault.o
	$(CC) -o srfLL2XYnf srfLL2XYnf.c geoproj_subs.o utlFault.o -lm

srfXY2LLKM:	srfXY2LLKM.c geoproj_subs.o utlFault.o
	$(CC) -o srfXY2LLKM srfXY2LLKM.c geoproj_subs.o utlFault.o -lm

statLL2XY:	statLL2XY.c geoproj_subs.o utlFault.o
	$(CC) -o statLL2XY statLL2XY.c geoproj_subs.o utlFault.o -lm

statLLXY:	statLLXY.c geoproj_subs.o utlFault.o
	$(CC) -o statLLXY statLLXY.c geoproj_subs.o utlFault.o -lm

conv3Comp:	conv3Comp.c utlFault.o utl3Comp.o 
		$(CC) -o conv3Comp conv3Comp.c utlFault.o utl3Comp.o -lm

conv3CompBB:	conv3CompBB.c utlFault.o utl3Comp.o 
		$(CC) -o conv3CompBB conv3CompBB.c utlFault.o utl3Comp.o -lm

conv3Comp2:	conv3Comp2.c utlFault.o utl3Comp.o 
		$(CC) -o conv3Comp2 conv3Comp2.c utlFault.o utl3Comp.o -lm

convertFaultInput:	convertFaultInput.c utlFault.o 
			$(CC) -o convertFaultInput utlFault.o convertFaultInput.c -lm

convertFaultInputKM:	convertFaultInputKM.c utlFault.o 
			$(CC) -o convertFaultInputKM utlFault.o convertFaultInputKM.c -lm

utlFault.o:	utlFault.c utlFault.h include.h
		$(CC) -c utlFault.c

geoproj_subs.o:	geoproj_subs.c include.h
		$(CC) -c geoproj_subs.c

utl3Comp.o:	utl3Comp.c utl3Comp.h
		$(CC) -c utl3Comp.c

statDist.o:	statDist.c statDist.h
		$(CC) -c statDist.c

copy:
	cp statLL2XY conv3Comp conv3CompBB convertFaultInput convertFaultInputKM srfLL2XYKM srfXY2LLKM ../bin
