all:
	mkdir -p bin
	cd Converter; make -f Makefile;
	cd getInFile; make -f makefile;
	cd DeconvStitch; make -f makefile;
	cd Syn1D; make -f makefile;
	cd FFSP; make -f makefile;
	cd timeHist; make -f Makefile;
	cd noah_code; make -f makefile;

clean:
	rm -rf bin;
	cd Converter; make -f Makefile clean;
	cd getInFile; make -f makefile clean;
	cd DeconvStitch; make -f makefile clean;
	cd Syn1D; make -f makefile clean;
	cd FFSP; make -f makefile clean;
	cd timeHist; make -f Makefile clean;
	cd noah_code; make -f makefile clean;
