FC=gfortran
#FC=ifort

PCOMP = ${FC} -O3 -o
OCOMP = ${FC} -O3 -c -o

OBJECTS = FileRead.o IntDiffFFT.o SpecResp.o InElastic.o BWfilt.o

all: GOF_MO GOF_MO_NGA GOF_PGX GOF_SpFit GOF_DCumEn GOF_InElFit GOF_SAFit16 GOF_CCFit GOF_FSComp GOF_SpecDurFit

FileRead_NGA.o: FileRead_NGA.f
	$(OCOMP) FileRead_NGA.o	FileRead_NGA.f

CBnga08.o: CBnga08.f
	$(OCOMP) CBnga08.o CBnga08.f

BAnga08.o: BAnga08.f
	$(OCOMP) BAnga08.o BAnga08.f

InElastic.o: InElastic.f
	$(OCOMP) InElastic.o InElastic.f

FileRead.o: FileRead.f
	$(OCOMP) FileRead.o FileRead.f

IntDiffFFT.o: IntDiffFFT.f
	$(OCOMP) IntDiffFFT.o IntDiffFFT.f

SpecResp.o: SpecResp.f
	$(OCOMP) SpecResp.o SpecResp.f

BWfilt.o: BWfilt.f
	$(OCOMP) BWfilt.o BWfilt.f

GOF_MO.o: GOF_MO.f
	$(OCOMP) GOF_MO.o GOF_MO.f

GOF_PGX.o: GOF_PGX.f
	$(OCOMP) GOF_PGX.o GOF_PGX.f

GOF_SpFit.o: GOF_SpFit.f
	$(OCOMP) GOF_SpFit.o GOF_SpFit.f

GOF_DCumEn.o: GOF_DCumEn.f
	$(OCOMP) GOF_DCumEn.o GOF_DCumEn.f

GOF_SAFit16.o: GOF_SAFit16.f
	$(OCOMP) GOF_SAFit16.o GOF_SAFit16.f

GOF_CCFit.o: GOF_CCFit.f
	$(OCOMP) GOF_CCFit.o GOF_CCFit.f

GOF_FSComp.o: GOF_FSComp.f
	$(OCOMP) GOF_FSComp.o GOF_FSComp.f

GOF_SpecDurFit.o: GOF_SpecDurFit.f
	$(OCOMP) GOF_SpecDurFit.o GOF_SpecDurFit.f

GOF_MO_NGA.o: GOF_MO_NGA.f
	$(OCOMP) GOF_MO_NGA.o GOF_MO_NGA.f

GOF_InElFit.o: GOF_InElFit.f
	$(OCOMP) GOF_InElFit.o GOF_InElFit.f

GOF_MO: $(OBJECTS) GOF_MO.o
	$(PCOMP) GOF_MO GOF_MO.o $(OBJECTS)
	cp GOF_MO ../bin/.

GOF_PGX: $(OBJECTS) GOF_PGX.o
	$(PCOMP) GOF_PGX GOF_PGX.o $(OBJECTS)
	cp GOF_PGX ../bin/.

GOF_SpFit: $(OBJECTS) GOF_SpFit.o
	$(PCOMP) GOF_SpFit GOF_SpFit.o $(OBJECTS)
	cp GOF_SpFit ../bin/.

GOF_DCumEn: $(OBJECTS) GOF_DCumEn.o
	$(PCOMP) GOF_DCumEn GOF_DCumEn.o $(OBJECTS)
	cp GOF_DCumEn ../bin/.

GOF_InElFit: $(OBJECTS) GOF_InElFit.o
	$(PCOMP) GOF_InElFit GOF_InElFit.o $(OBJECTS)
	cp GOF_InElFit ../bin/.

GOF_SAFit16: $(OBJECTS) GOF_SAFit16.o
	$(PCOMP) GOF_SAFit16 GOF_SAFit16.o $(OBJECTS)
	cp GOF_SAFit16 ../bin/.

GOF_CCFit: $(OBJECTS) GOF_CCFit.o
	$(PCOMP) GOF_CCFit GOF_CCFit.o $(OBJECTS)
	cp GOF_CCFit ../bin/.

GOF_FSComp: $(OBJECTS) GOF_FSComp.o
	$(PCOMP) GOF_FSComp GOF_FSComp.o $(OBJECTS)
	cp GOF_FSComp ../bin/.

GOF_SpecDurFit: $(OBJECTS) GOF_SpecDurFit.o
	$(PCOMP) GOF_SpecDurFit GOF_SpecDurFit.o $(OBJECTS)
	cp GOF_SpecDurFit ../bin/.

GOF_MO_NGA: $(OBJECTS) CBnga08.o BAnga08.o GOF_MO_NGA.o 
	$(PCOMP) GOF_MO_NGA CBnga08.o BAnga08.o GOF_MO_NGA.o $(OBJECTS)
	cp GOF_MO_NGA ../bin/.

clean:
	rm -f *.o GOF_MO_NGA GOF_CCFit GOF_DCumEn GOF_GSComp GOF_InElFit GOF_MO GOF_PGX GOF_SAFit16 GOF_SpFit GOF_FSComp GOF_SpecDurFit
