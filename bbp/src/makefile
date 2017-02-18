all: version
	cd gp; make -f makefile;
	cd sdsu; make -f makefile;
	cd ucsb; make -f makefile;
	cd ucb; make -f makefile;
	cd uwo; make -f makefile;
	cd unr; make -f makefile;
	cd irikura; make -f makefile;
	cd usgs; make -f makefile;

version:
	gcc --version > gcc.version
	gfortran --version > gfortran.version

clean:
	rm -f gcc.version gfortran.version
	cd gp; make -f makefile clean;
	cd sdsu; make -f makefile clean;
	cd ucsb; make -f makefile clean;
	cd ucb; make -f makefile clean;
	cd uwo; make -f makefile clean;
	cd unr; make -f makefile clean;
	cd irikura; make -f makefile clean;
	cd usgs; make -f makefile clean;
