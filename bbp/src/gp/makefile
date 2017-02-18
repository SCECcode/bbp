all:
	mkdir -p bin;
	cd getpar/src; make -f Makefile;
	cd StandRupFormat; make -f makefile;
	cd StochSim/V5.4; make -f makefile;
	cd WccFormat/Progs; make -f makefile;
	cd ModelCords; make -f makefile;
	cd GenRandV5.0; make -f makefile;
	cd GoodFit; make -f makefile;
	cd JordanBailey; make -f makefile;

clean:
	rm -rf bin;
	cd getpar/src; make -f Makefile clean;
	cd StandRupFormat; make -f makefile clean;
	cd StochSim/V5.4; make -f makefile clean;
	cd WccFormat/Progs; make -f makefile clean;
	cd ModelCords; make -f makefile clean;
	cd GenRandV5.0; make -f makefile clean;
	cd GoodFit; make -f makefile clean;
	cd JordanBailey; make -f makefile clean;
