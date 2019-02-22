all:
	mkdir -p bin;
	cd gen_srf; make -f makefile
	cd srf2grns; make -f makefile
	cd statgreen; make -f makefile
	cd statgreen2; make -f makefile
	cd greenscale; make -f makefile
clean:
	rm -rf bin;
	cd gen_srf; make -f makefile clean;
	cd srf2grns; make -f makefile clean;
	cd statgreen; make -f makefile clean;
	cd statgreen2; make -f makefile clean;
	cd greenscale; make -f makefile clean;
