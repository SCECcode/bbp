all:
	mkdir -p bin;
	cd Convert; make -f makefile;
	cd bbtoolbox; make -f makefile;
	cd MOGof; make -f makefile;

clean:
	rm -rf bin
	cd Convert; make -f makefile clean;
	cd bbtoolbox; make -f makefile clean;
	cd MOGof; make -f makefile clean;
