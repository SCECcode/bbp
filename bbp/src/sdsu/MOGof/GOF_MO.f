C If you receive a segmentation fault try to set the stack size to a larger value
C i.e. ulimit -s (stack size)
C A program created by John Mayhew starting on September 22, 2008
C Read in the inital parameters into the core portion of the the code
C then send it to the main_comp subroutine to make the comparision.

      program main
      implicit none

      integer            :: ntS, ntR, nt2, ntB, ntT, numSta, numHEAD
      real               :: dtS, dtR, dt2, leng, base, CutL, CutH
      real               :: freqB, dtB, ntBB, dtT
      character*80          :: inputS, inputR
      character*80          :: TAG
      real, dimension(20)   :: VAL

      freqB = 15.0
      dtB = 1.0/(2.0*freqB)

C Get the parameters from PARAM.dat
      call GetParam(inputS, inputR, ntS, ntR, leng, 
     +     TAG, numSta, VAL, numHEAD, CutL, CutH)

C Detirmine the appropriate dt and nt to use for the analysis 
      dtS = leng/real(ntS-1)
      dtR = leng/real(ntR-1)
      dtT = 0.551*CutL*0.9
      ntB = int(leng/dtB) + 1
      ntT = int(leng/dtT) + 1
      ntBB = max(ntS, ntR, ntB, ntT)

      base = LOG10(ntBB)/LOG10(2.0)
	nt2 = 2**int(base+1.0)

      dt2 = leng/real(nt2-1)

C Export values to the MAIN program
      call MAIN_COMP(inputS, inputR, dtS, dtR, dt2, ntS, ntR, 
     +    nt2, leng, TAG, numSta, VAL, numHEAD, CutL, CutH)

      end program main
C --------------------------------------------------------------------
C --------------------------MAIN_COMP---------------------------------
C Subroutine to run MAIN program
C --------------------------------------------------------------------
C --------------------------------------------------------------------
      subroutine MAIN_COMP(inputS, inputR, dtS, dtR, dt2, ntS, 
     +ntR, nt2, leng, TAG, numSta, VAL, numHEAD, CutL, CutH)
      implicit none
C---------------------------------------------------------------------
      integer     :: ntS, ntR, i, j, k, nt2, ii, jj, kk, numHEAD
      integer     :: numSta, staNum, nrec, tempI, ntS2, ntR2
      integer     :: DLAX, DLAY, DLAZ
C---------------------------------------------------------------------
C-------------------------Filter Indicies-----------------------------
      integer     :: L16, H16, L991, H991
C---------------------------------------------------------------------
      real  :: NGA_PGA, NGA_PGV, ttTemp
      real  :: dtS, dtR, dt2, leng
      real  :: tmpxMax, tmpyMax, tmpzMax 
      real  :: siteV, siteVX, siteVY, siteVZ
      real  :: T, tt, r1, r2, r3
      real  :: PGVxR, PGVyR, PGVzR, PGAxR, PGAyR, PGAzR
      real  :: PGVxS, PGVyS, PGVzS, PGAxS, PGAyS, PGAzS
      real  :: PGDxR, PGDyR, PGDzR, PGDxS, PGDyS, PGDzS
      real  :: temp, temp1, temp2, temp3, temp4, temp5
      real  :: temp6, temp7, temp8, temp9, CosT, tempX, tempY
      real  :: Rjb, Rrup, ZvsISO, Vs30, dy, CumDiff, SpecFitT
      real  :: CCX, CCY, CCZ, CVX, CVY, CVZ, CCfit, CVfit
      real  :: EnDur, SpecDurFit, SAfit16, InElFit
      real  :: PSAfit, PGAfit, PGVfit, PGDfit
      real  :: SpecDurFitX, SAfit16X, InElFitX
      real  :: PSAfitX, PGAfitX, PGVfitX, PGDfitX
      real  :: SpecDurFitY, SAfit16Y, InElFitY
      real  :: PSAfitY, PGAfitY, PGVfitY, PGDfitY
      real  :: SpecDurFitZ, SAfit16Z, InElFitZ, InElFitH
      real  :: PSAfitZ, PGAfitZ, PGVfitZ, PGDfitZ
      real  :: PPMCCfit, DCumEnX, DCumEnY, DCumEnZ
      real  :: LSspec, LSaccX, LSvelX, LSdispX
      real  :: LSaccY, LSvelY, LSdispY, TR, TS
      real  :: LSaccZ, LSvelZ, LSdispZ, PSAR, PSAS
      real  :: DLS, DLR, DHS, DHR, dDur
      real  :: SDisCmpX, SDisCmpY, SDisCmpZ 
      real  :: SDurCmpX, SDurCmpY, SDurCmpZ
      real  :: SRCmpX, SRCmpY, SRCmpZ, SpFitVal
      real  :: TO3X, TO3Y, TO3Z, DFX, DFY, DFZ
      real  :: DCumEn, CutL, CutH
      real  :: FScomp, FScompX, FScompY, FScompZ
C---------------------------------------------------------------------
      real, dimension(20)   :: VAL
      real, dimension(ntS)  :: SSx, SSy, SSz
      real, dimension(ntR)  :: SRx, SRy, SRz
      real, dimension(nt2)  :: LFSx, LFSy, LFSz
      real, dimension(nt2)  :: LFRx, LFRy, LFRz
C---------------------------------------------------------------------
      real, dimension(nt2)  :: SSx2, SSy2, SSz2
      real, dimension(nt2)  :: SRx2, SRy2, SRz2
      real, dimension(nt2)  :: ASx, ASy, ASz, ARx, ARy, ARz
      real, dimension(nt2)  :: VSx, VSy, VSz, VRx, VRy, VRz
      real, dimension(nt2)  :: DSx, DSy, DSz, DRx, DRy, DRz
      real, dimension(nt2)  :: outR, outI, zeroNT2, tempFILT
      integer, dimension(16):: OUTPUT2
C---------------------------------------------------------------------
      real, dimension(numSta,991):: rsqSRR, rsqSRS, rsqSDR, rsqSDS
      real, dimension(991):: SRmBIAS, SRmERROR, SDmBIAS, SDmERROR
      real, dimension(1):: PGAmBIAS, PGAmERROR, PGVmBIAS, PGVmERROR
      real, dimension(991):: perList, SRSX, SRSY, SRRX, SRRY
      real, dimension(991):: HSRR, HSRS, CumSRR, CumSRS, VSRR, VSRS
      real, dimension(991):: SDSX, SDSY, SDRX, SDRY, VSDR, VSDS
      real, dimension(991):: SDispLRX, SDispLRY, SDispLRZ
      real, dimension(991):: SDispLSX, SDispLSY, SDispLSZ
      real, dimension(19):: yListS, sigListS
      real, dimension(16):: period, SpecHS
      real, dimension(16):: yListR, sigListR, SpecHR, SpecSZ
      real, dimension(16):: SpecSX, SpecSY, SpecRX, SpecRY, DurRZ
      real, dimension(16):: DurSX, DurRX, DurSY, DurRY, DurSZ
      real, dimension(16):: SpDisRX, SpDisRY, SpDisRZ, SpecRZ
      real, dimension(16)  :: SpDisSX, SpDisSY, SpDisSZ
      real, dimension(16,19):: maxDXR, maxDYR, maxDZR
      real, dimension(16,19):: maxDXS, maxDYS, maxDZS
      real, dimension(16,19):: RatioRX, RatioRY, RatioRZ
      real, dimension(16,19):: RatioSX, RatioSY, RatioSZ
      real, dimension(12)  :: SpCmpMDev
      real, dimension(13)  :: indVec
      real, dimension(19)  :: R
      real, dimension(16,19):: GOFIEX, GOFIEY, GOFIEZ, GOFIEH
C---------------------------------------------------------------------
      character*80          :: inputS, inputR
      character*80          :: TAG, header
C---------------------------------------------------------------------
      real, dimension(nt2):: CumSVin, CumSVout
      real, dimension(nt2):: CumRVin, CumRVout, TO1, TO2
      real, dimension(16):: RatCmpX, RatCmpY, RatCmpZ, RatCmpH
      real, dimension(16):: SpFX, SpFY, SpFZ
      real, dimension(2) :: BSFx, BSFy, BSFz, WSFx, WSFy, WSFz
      real, dimension(numSta,1):: PGDS, PGVS, PGDR, PGVR, PGAR, PGAS
C---------------------------------------------------------------------
      data period/.1,.15,.2,.25,.3,.4,.5,.75,
     +1,1.5,2,3,4,5,7.5,10/
C Deterministic R values set by Baker
C      data R/.1,.2,.3,.4,.5,.6,.7,.8,.9,1,2,3,4,5,6,7,8,9,10/
      data R/1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,
     +7.0,7.5,8.0,8.5,9.0,9.5,10.0/

      data indVec/1,31,51,76,118,201,451,951,961,991,1041,1091,1141/
C---------------------------------------------------------------------
      nrec = 0
      ntS2 = ntS
      ntR2 = ntR
C      if(dtS*ntS .LT. dtR*ntR) ntR2 = nint((dtS*ntS)/dtR)
C      if(dtR*ntR .LT. dtS*ntS) ntS2 = nint((dtR*ntR)/dtS)
C---------------------------------------------------------------------
C Open the output files
         open(852, file='out/SpecDur16GOF_XYZ.out', status='unknown')
         open(86, file='out/SpecDurX.out', status='unknown')
         open(87, file='out/SpecDurY.out', status='unknown')
         open(88, file='out/SpecOutH.out', status='unknown')

         open(93, file='out/SpecDurSX.out', status='unknown')
         open(94, file='out/SpecDurSY.out', status='unknown')
         open(95, file='out/SpecDurRX.out', status='unknown')
         open(96, file='out/SpecDurRY.out', status='unknown')

         open(971, file='out/SpecOutSX.out', status='unknown')
         open(981, file='out/SpecOutSY.out', status='unknown')
         open(972, file='out/SpecOutRX.out', status='unknown')
         open(982, file='out/SpecOutRY.out', status='unknown')

         open(97, file='out/SpecOutHS.out', status='unknown')
         open(98, file='out/SpecOutHR.out', status='unknown')

         open(3030, file='out/ElD_R.out', status='unknown')
         open(3131, file='out/ElD_S.out', status='unknown')

         open(991, file='out/InElS.19.16.NumSta.out', status='unknown')
         open(992, file='out/InElR.19.16.NumSta.out', status='unknown')

         open(98277, file='out/ErrorBias.out', status='unknown')

      open(9875, file='out/FileS.out', status='unknown')
      open(9876, file='out/FileR.out', status='unknown')
      open(98752, file='out/PGA.out', status='unknown')
      open(98762, file='out/PGV.out', status='unknown')
      open(357, file='GOF.log', status='unknown')
      open(38, file='out/GOF.list', status='unknown')

      open(1111, file='out/GOF_PGA.list', status='unknown')
      open(2222, file='out/GOF_PGV.list', status='unknown')
      open(3333, file='out/GOF_PGD.list', status='unknown')
      open(4444, file='out/GOF_PSA.list', status='unknown')
      open(5555, file='out/DELAY.list', status='unknown')
      open(6666, file='out/GOF_SPECFIT.list', status='unknown')
      open(7777, file='out/GOF_ENERGYFIT.list', status='unknown')
      open(8888, file='out/GOF_InElEl.list', status='unknown')
      open(8111, file='out/GOF_InElEl_All.list', status='unknown')
      open(9999, file='out/GOF_SAFIT.list', status='unknown')
      open(111111, file='out/GOF_SPECDUR.list', status='unknown')
      open(121212, file='out/GOF_CROSSCOR.list', status='unknown')
      open(131313, file='out/GOF_COVAR.list', status='unknown')
      open(141414, file='out/GOF_DUR.list', status='unknown')
      open(151515, file='out/GOF_FS.list', status='unknown')

      open(232323, file='out/GOF_SA16X_DATA.dat', status='unknown')
      write(232323, *) '%Period - GOF - Set 1 - Set 2'
      open(242424, file='out/GOF_SA16Y_DATA.dat', status='unknown')
      write(242424, *) '%Period - GOF - Set 1 - Set 2'
      open(252525, file='out/GOF_SA16Z_DATA.dat', status='unknown')
      write(252525, *) '%Period - GOF - Set 1 - Set 2'
C---------------------------------------------------------------------
      write(4444,*) '% Period, Set 1, Set 2'
      write(357,*) '% Set 1 = ', inputR(1:20)
      write(357,*) '% Set 2 = ', inputS(1:20)
      write(38,*) '% Set 1 = ', inputR(1:20)
      write(38,*) '% Set 2 = ', inputS(1:20)
      write(*,*) '% Set 1 = ', inputR(1:20)
      write(*,*) '% Set 2 = ', inputS(1:20)
C----------------Loop Through the Stations----------------------------
      open(1, file=inputS)
      open(2, file=inputR)
C---------------------------------------------------------------------
      do staNum = 1,numSta
         write(*,*) '-------------------------------'
         write(*,*) 'Calculating Metrics for Station  ', staNum

C Read in and resize the Synthetic Seismogram
      if(numHEAD .EQ. 0) then
      else
      do i=1,numHEAD
         read(1,*) header
      enddo
      endif
      do i = 1,ntS2
         read(1,*)  ttTemp, SSx(i), SSy(i), SSz(i)
      enddo
      call ReSize(SSx, SSx2, dtS, dt2, ntS2, nt2)
      call ReSize(SSy, SSy2, dtS, dt2, ntS2, nt2)
      call ReSize(SSz, SSz2, dtS, dt2, ntS2, nt2)

C Read in and resize the raw data seismogram
      if(numHEAD .EQ. 0) then
      else
      do i=1,numHEAD
         read(2,*) header
      enddo
      endif
      do i=1,ntR2
         read(2,*)  ttTemp, SRx(i), SRy(i), SRz(i)
      enddo
      call ReSize(SRx, SRx2, dtR, dt2, ntR2, nt2)
      call ReSize(SRy, SRy2, dtR, dt2, ntR2, nt2)
      call ReSize(SRz, SRz2, dtR, dt2, ntR2, nt2)
C---------------------------------------------------------------------
C Data Filter
      call FILT(SSx2, nt2, 'BP', 1.0/CutH, 1.0/CutL, dt2)
      call FILT(SSy2, nt2, 'BP', 1.0/CutH, 1.0/CutL, dt2)
      call FILT(SSz2, nt2, 'BP', 1.0/CutH, 1.0/CutL, dt2)
      call FILT(SRx2, nt2, 'BP', 1.0/CutH, 1.0/CutL, dt2)
      call FILT(SRy2, nt2, 'BP', 1.0/CutH, 1.0/CutL, dt2)
      call FILT(SRz2, nt2, 'BP', 1.0/CutH, 1.0/CutL, dt2)

      do i = 1,nt2
         tt = real(i-1)*dt2
         write(9876,*) tt, SRx2(i), SRy2(i), SRz2(i)
      enddo 
      do i = 1,nt2
         tt = real(i-1)*dt2
         write(9875,*) tt, SSx2(i), SSy2(i), SSz2(i)
      enddo 
C---------------------------------------------------------------------
C Generate time series of Acceleration, velocity and displacement
      if(TAG .EQ. 'A') then
C If input is Acceleration
C Set 2 data
            ASx(:) = SSx2(:)
            ASy(:) = SSy2(:)
            ASz(:) = SSz2(:)
            call DINTEG(SSx2, dt2, nt2, VSx)
            call DINTEG(SSy2, dt2, nt2, VSy)
            call DINTEG(SSz2, dt2, nt2, VSz)
            call DINTEG(VSx, dt2, nt2, DSx)
            call DINTEG(VSy, dt2, nt2, DSy)
            call DINTEG(VSz, dt2, nt2, DSz)
C Set 1 data
            ARx(:) = SRx2(:)
            ARy(:) = SRy2(:)
            ARz(:) = SRz2(:)
            call DINTEG(SRx2, dt2, nt2, VRx)
            call DINTEG(SRy2, dt2, nt2, VRy)
            call DINTEG(SRz2, dt2, nt2, VRz)
            call DINTEG(VRx, dt2, nt2, DRx)
            call DINTEG(VRy, dt2, nt2, DRy)
            call DINTEG(VRz, dt2, nt2, DRz)
      elseif(TAG .EQ. 'V') then
C If input is Velocity
C Set 2 data
            VSx(:) = SSx2(:)
            VSy(:) = SSy2(:)
            VSz(:) = SSz2(:)
            call DINTEG(SSx2, dt2, nt2, DSx)
            call DINTEG(SSy2, dt2, nt2, DSy)
            call DINTEG(SSz2, dt2, nt2, DSz)
            call DIFF(SSx2, dt2, nt2, ASx)
            call DIFF(SSy2, dt2, nt2, ASy)
            call DIFF(SSz2, dt2, nt2, ASz)
C Set 1 data
            VRx(:) = SRx2(:)
            VRy(:) = SRy2(:)
            VRz(:) = SRz2(:)
            call DINTEG(SRx2, dt2, nt2, DRx)
            call DINTEG(SRy2, dt2, nt2, DRy)
            call DINTEG(SRz2, dt2, nt2, DRz)
            call DIFF(SRx2, dt2, nt2, ARx)
            call DIFF(SRy2, dt2, nt2, ARy)
            call DIFF(SRz2, dt2, nt2, ARz)
      elseif(TAG .EQ. 'D') then
C If input is Displacement
C Set 2 data
            DSx(:) = SSx2(:)
            DSy(:) = SSy2(:)
            DSz(:) = SSz2(:)
            call DIFF(SSx2, dt2, nt2, VSx)
            call DIFF(SSy2, dt2, nt2, VSy)
            call DIFF(SSz2, dt2, nt2, VSz)
            call DIFF(VSx, dt2, nt2, ASx)
            call DIFF(VSy, dt2, nt2, ASy)
            call DIFF(VSz, dt2, nt2, ASz)
C Set 1 data
            DRx(:) = SRx2(:)
            DRy(:) = SRy2(:)
            DRz(:) = SRz2(:)
            call DIFF(SRx2, dt2, nt2, VRx)
            call DIFF(SRy2, dt2, nt2, VRy)
            call DIFF(SRz2, dt2, nt2, VRz)
            call DIFF(VRx, dt2, nt2, ARx)
            call DIFF(VRy, dt2, nt2, ARy)
            call DIFF(VRz, dt2, nt2, ARz)
      endif

C PGA
      PGAxR = MAXVAL(ABS(ARx))
      PGAyR = MAXVAL(ABS(ARy))
      PGAxS = MAXVAL(ABS(ASx))
      PGAyS = MAXVAL(ABS(ASy))
      PGAzR = MAXVAL(ABS(ARz))
      PGAzS = MAXVAL(ABS(ASz))
      write(98752,*) PGAxR, PGAyR, PGAxS, PGAyS

C PGV
      PGVxR = MAXVAL(ABS(VRx))
      PGVyR = MAXVAL(ABS(VRy))
      PGVxS = MAXVAL(ABS(VSx))
      PGVyS = MAXVAL(ABS(VSy))
      PGVzR = MAXVAL(ABS(VRz))
      PGVzS = MAXVAL(ABS(VSz))
      write(98762,*) PGVxR, PGVyR, PGVxS, PGVyS

C PGD
      PGDxR = MAXVAL(ABS(DRx))
      PGDyR = MAXVAL(ABS(DRy))
      PGDxS = MAXVAL(ABS(DSx))
      PGDyS = MAXVAL(ABS(DSy))
      PGDzR = MAXVAL(ABS(DRz))
      PGDzS = MAXVAL(ABS(DSz))
C---------------------------------------------------------------------
C Calculate PGA, PGV, and PGD for both seismograms
      PGAS(staNum,1) = 0.0
      PGAR(staNum,1) = 0.0
      PGVS(staNum,1) = 0.0
      PGVR(staNum,1) = 0.0
      PGDS(staNum,1) = 0.0
      PGDR(staNum,1) = 0.0
      do i = 1,nt2
           if(sqrt(abs(ASx(i)*ASy(i))) .gt. PGAS(staNum,1)) then 
              PGAS(staNum,1)=sqrt(abs(ASx(i)*ASy(i)))
           endif
           if(sqrt(abs(ARx(i)*ARy(i))) .gt. PGAR(staNum,1)) then 
              PGAR(staNum,1)=sqrt(abs(ARx(i)*ARy(i)))
           endif
           if(sqrt(abs(VSx(i)*VSy(i))) .gt. PGVS(staNum,1)) then 
              PGVS(staNum,1)=sqrt(abs(VSx(i)*VSy(i)))
           endif
           if(sqrt(abs(VRx(i)*VRy(i))) .gt. PGVR(staNum,1)) then 
              PGVR(staNum,1)=sqrt(abs(VRx(i)*VRy(i)))
           endif
           if(sqrt(abs(DSx(i)*DSy(i))) .gt. PGDS(staNum,1)) then 
              PGDS(staNum,1)=sqrt(abs(DSx(i)*DSy(i)))
           endif
           if(sqrt(abs(DRx(i)*DRy(i))) .gt. PGDR(staNum,1)) then 
              PGDR(staNum,1)=sqrt(abs(DRx(i)*DRy(i)))
           endif
      enddo
C---------------------------------------------------------------------
C---------------------------------------------------------------------
      L16 = 0
      H16 = 0
      L991 = 0
      H991 = 0
      do ii = 1,16
C Set 2 calculation
          call RESPONSESPEC(ASx, dt2, nt2, period(ii), SpecSX(ii), 
     +         DurSX(ii), SpDisSX(ii))
          call RESPONSESPEC(ASy, dt2, nt2, period(ii), SpecSY(ii), 
     +         DurSY(ii), SpDisSY(ii))
          call RESPONSESPEC(ASz, dt2, nt2, period(ii), SpecSZ(ii), 
     +         DurSZ(ii), SpDisSZ(ii))
          SpecHS(ii) = SQRT(SpecSX(ii)*SpecSY(ii))
C---------------------------------------------------------------------
C Set 1 data calculation
          call RESPONSESPEC(ARx, dt2, nt2, period(ii), SpecRX(ii), 
     +         DurRX(ii), SpDisRX(ii))
          call RESPONSESPEC(ARy, dt2, nt2, period(ii), SpecRY(ii), 
     +         DurRY(ii), SpDisRY(ii))
          call RESPONSESPEC(ARz, dt2, nt2, period(ii), SpecRZ(ii), 
     +         DurRZ(ii), SpDisRZ(ii))
C          write(*,*) SpecSX(ii), SpecRX(ii)
          SpecHR(ii) = SQRT(SpecRX(ii)*SpecRY(ii))
      write(3030,*) period(ii),SpDisRX(ii), SpDisRY(ii),SpDisRZ(ii)
      write(3131,*) period(ii),SpDisSX(ii), SpDisSY(ii),SpDisSZ(ii)
C---------------------------------------------------------------------
C InElastic calculation Set 2
C Calculating dy according to the relationship set forth by Tothong and Cornell 2006
C "An empirical Ground-Motion Attenuation Relation for Inelastic Spectral Displacement"
          do jj = 1,19
            dy = SpDisSX(ii)/R(jj)
            call InElastic(period(ii),dy,ASx,dt2,nt2,maxDXS(ii,jj))
            RatioSX(ii,jj) = maxDXS(ii,jj)/SpDisSX(ii)
            dy = SpDisSY(ii)/R(jj)
            call InElastic(period(ii),dy,ASy,dt2,nt2,maxDYS(ii,jj))
            RatioSY(ii,jj) = maxDYS(ii,jj)/SpDisSY(ii)
            dy = SpDisSZ(ii)/R(jj)
            call InElastic(period(ii),dy,ASz,dt2,nt2,maxDZS(ii,jj))
            RatioSZ(ii,jj) = maxDZS(ii,jj)/SpDisSZ(ii)

      write(991,*) RatioSX(ii,jj), RatioSY(ii,jj), RatioSZ(ii,jj)
C---------------------------------------------------------------------
C InElastic calculation Set 1
                 dy = SpDisRX(ii)/R(jj)
            call InElastic(period(ii),dy,ARx,dt2,nt2,maxDXR(ii,jj))
                 RatioRX(ii,jj) = maxDXR(ii,jj)/SpDisRX(ii)
                 dy = SpDisRY(ii)/R(jj)
            call InElastic(period(ii),dy,ARy,dt2,nt2,maxDYR(ii,jj))
                 RatioRY(ii,jj) = maxDYR(ii,jj)/SpDisRY(ii)
                 dy = SpDisRZ(ii)/R(jj)
            call InElastic(period(ii),dy,ARz,dt2,nt2,maxDZR(ii,jj))
                 RatioRZ(ii,jj) = maxDZR(ii,jj)/SpDisRZ(ii)

      write(992,*) RatioRX(ii,jj), RatioRY(ii,jj), RatioRZ(ii,jj)
          enddo
          if(period(ii) .LT. CutL) L16 = ii 
          if(period(ii) .LE. CutH) H16 = ii 
      enddo
          L16 = L16+1
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C Calculate the Spectral Response Vector 0.05 - 1 sec periods
C HSR(R,S) --> Horizonal Spectral Response (Raw, Synthetic)
C VSR(R,S) --> Vertical Spectral Response (Raw, Synthetic)
      do j = 1,901
           perList(j) = 0.001*real(j-1)+0.1
C Set 2 calculation
          call RESPONSESPEC(ASx, dt2, nt2, perList(j), SRSX(j), SDSX(j), 
     +           SDispLSX(j))
          call RESPONSESPEC(ASy, dt2, nt2, perList(j), SRSY(j), SDSY(j), 
     +           SDispLSY(j))
          HSRS(j) = SQRT(SRSX(j)*SRSY(j))
          call RESPONSESPEC(ASz, dt2, nt2, perList(j), VSRS(j), VSDS(j), 
     +           SDispLSZ(j))
C Set 1 data calculation
          call RESPONSESPEC(ARx, dt2, nt2, perList(j), SRRX(j), SDRX(j), 
     +           SDispLRX(j))
          call RESPONSESPEC(ARy, dt2, nt2, perList(j), SRRY(j), SDRY(j), 
     +           SDispLRY(j))
          HSRR(j) = SQRT(SRRX(j)*SRRY(j))
          call RESPONSESPEC(ARz, dt2, nt2, perList(j), VSRR(j), VSDR(j), 
     +           SDispLRZ(j))
      enddo

C Calculate the Spectral Response Vector 1.1 - 20 sec periods
      do j = 902,991
           perList(j) = 0.1*real(j-901)+1.0
C Set 2 calculation
          call RESPONSESPEC(ASx, dt2, nt2, perList(j), SRSX(j), SDSX(j), 
     +           SDispLSX(j))
          call RESPONSESPEC(ASy, dt2, nt2, perList(j), SRSY(j), SDSY(j), 
     +           SDispLSY(j))
          HSRS(j) = SQRT(SRSX(j)*SRSY(j))
          call RESPONSESPEC(ASz, dt2, nt2, perList(j), VSRS(j), VSDS(j), 
     +           SDispLSZ(j))
C Set 1 data calculation
          call RESPONSESPEC(ARx, dt2, nt2, perList(j), SRRX(j), SDRX(j), 
     +           SDispLRX(j))
          call RESPONSESPEC(ARy, dt2, nt2, perList(j), SRRY(j), SDRY(j), 
     +           SDispLRY(j))
          HSRR(j) = SQRT(SRRX(j)*SRRY(j))
          call RESPONSESPEC(ARz, dt2, nt2, perList(j), VSRR(j), VSDR(j), 
     +           SDispLRZ(j))
      enddo
      do j = 1,991
          if(perList(j) .LT. CutL) L991 = j 
          if(perList(j) .LE. CutH) H991 = j 
      enddo
          L991 = L991+1
          write(*,*) 'LowCut = ', CutL, period(L16), perList(L991)
          write(*,*) 'HighCut = ', CutH, period(H16), perList(H991)
C---------------------------------------------------------------------
C--------------------END OF METRICS CALCULATIONS----------------------
C---------------------------------------------------------------------
      write(*,*) 'Calculating GOF'
C---------------------------------------------------------------------
C  Compute the GOF measure

      call GOF(ABS(ASx), ABS(ARx), nt2, 1, nt2, LSaccX)
      call GOF(ABS(VSx), ABS(VRx), nt2, 1, nt2, LSvelX)
      call GOF(ABS(DSx), ABS(DRx), nt2, 1, nt2, LSdispX)

      call GOF(ABS(ASy), ABS(ARy), nt2, 1, nt2, LSaccY)
      call GOF(ABS(VSy), ABS(VRy), nt2, 1, nt2, LSvelY)
      call GOF(ABS(DSy), ABS(DRy), nt2, 1, nt2, LSdispY)

      call GOF(ABS(ASz), ABS(ARz), nt2, 1, nt2, LSaccZ)
      call GOF(ABS(VSz), ABS(VRz), nt2, 1, nt2, LSvelZ)
      call GOF(ABS(DSz), ABS(DRz), nt2, 1, nt2, LSdispZ)
C---------------------------------------------------------------------
C PGA, PGV, PGD, PSA Fits

	  call GOF(PGAxR, PGAxS, 1, 1, 1, PGAfitX)
	  call GOF(PGAyR, PGAyS, 1, 1, 1, PGAfitY)
	  call GOF(PGAzR, PGAzS, 1, 1, 1, PGAfitZ)
      PGAfit = (PGAfitX+PGAfitY+PGAfitZ)/3.0

	  call GOF(PGVxR, PGVxS, 1, 1, 1, PGVfitX)
	  call GOF(PGVyR, PGVyS, 1, 1, 1, PGVfitY)
	  call GOF(PGVzR, PGVzS, 1, 1, 1, PGVfitZ)
      PGVfit = (PGVfitX+PGVfitY+PGVfitZ)/3.0

	  call GOF(PGDxR, PGDxS, 1, 1, 1, PGDfitX)
	  call GOF(PGDyR, PGDyS, 1, 1, 1, PGDfitY)
	  call GOF(PGDzR, PGDzS, 1, 1, 1, PGDfitZ)
      PGDfit = (PGDfitX+PGDfitY+PGDfitZ)/3.0


      PSAR = 0.0
      PSAS = 0.0
      TR = 0.0
      TS = 0.0
      do i = 1,991
         if(SRRX(i) .GT. PSAR) then 
            PSAR = SRRX(i)
            TR = perList(i)
         endif
         if(SRRY(i) .GT. PSAR) then 
            PSAR = SRRY(i)
            TR = perList(i)
         endif
         if(VSRR(i) .GT. PSAR) then 
            PSAR = VSRR(i)
            TR = perList(i)
         endif

         if(SRSX(i) .GT. PSAS) then 
            PSAS = SRSX(i) 
            TS = perList(i)
         endif
         if(SRSY(i) .GT. PSAS) then 
            PSAS = SRSY(i) 
            TS = perList(i)
         endif
         if(VSRS(i) .GT. PSAS) then 
            PSAS = VSRS(i)
            TS = perList(i)
         endif
      enddo
      
	  call GOF(MAXVAL(SRSX), MAXVAL(SRRX), 1, 1, 1, PSAfitX)
	  call GOF(MAXVAL(SRSY), MAXVAL(SRRY), 1, 1, 1, PSAfitY)
	  call GOF(MAXVAL(VSRR), MAXVAL(VSRS), 1, 1, 1, PSAfitZ)    
      PSAfit = (PSAfitX+PSAfitY+PSAfitZ)/3.0
C---------------------------------------------------------------------
C Calculate Spectral Fit
      SRCmpX = 0.0
      SRCmpY = 0.0
      SRCmpZ = 0.0
      
      call GOF(SRSX, SRRX, 991, L991, H991, SRCmpX)
      call GOF(SRSY, SRRY, 991, L991, H991, SRCmpY)
      call GOF(VSRS, VSRR, 991, L991, H991, SRCmpZ)

      call GOF(SDSX, SDRX, 991, L991, H991, SDurCmpX)
      call GOF(SDSY, SDRY, 991, L991, H991, SDurCmpY)
      call GOF(VSDS, VSDR, 991, L991, H991, SDurCmpZ)

      call GOF(SDispLSX, SDispLRX, 991, L991, H991, SDisCmpX)
      call GOF(SDispLSY, SDispLRY, 991, L991, H991, SDisCmpY)
      call GOF(SDispLSZ, SDispLRZ, 991, L991, H991, SDisCmpZ)

C      call GOF(HSRS, HSRR, 991, L991, 991, LSspec)

C---------------------------------------------------------------------
C Spectral Fit value
      SpFitVal = (SRCmpX+SRCmpY+SRCmpZ)/3.0
C---------------------------------------------------------------------
C Duration Calculation and Comparison
C      do j = 1,nt2
C           CumSVin(j) = (VSx(j)**2.0+VSy(j)**2.0+VSz(j)**2.0)
C           CumRVin(j) = (VRx(j)**2.0+VRy(j)**2.0+VRz(j)**2.0)
C      enddo

C      call DINTEG(CumSVin, dt2, nt2, CumSVout)
C      call DINTEG(CumRVin, dt2, nt2, CumRVout)

      DCumEn = 0.0
      call SDOTP(CumSVout, CumRVout, nt2, CosT)
C      call GOF(CumSVout, CumRVout, nt2, nt2, nt2, DCumEn)

C      do j = 1,nt2
C          CumSVout(j) = CumSVout(j)/CumSVout(nt2)
C           CumRVout(j) = CumRVout(j)/CumRVout(nt2)
C           if(CumSVout(j).le.0.05) DLS = REAL(j-1)*dt2
C           if(CumRVout(j).le.0.05) DLR = REAL(j-1)*dt2
C           if(CumSVout(j).le.0.75) DHS = REAL(j-1)*dt2
C           if(CumRVout(j).le.0.75) DHR = REAL(j-1)*dt2
C      enddo
C      DHS  = DHS-DLS
C      DHR  = DHR-DLR
C      call GOF(DHR, DHS, 1, 1, 1, EnDur)

C Individual vectors comparisons
C X component comparison
      call DINTEG(ABS(VSx)**2.0, dt2, nt2, TO1)
      call DINTEG(ABS(VRx)**2.0, dt2, nt2, TO2)
      call GOF(TO1, TO2, nt2, nt2, nt2, DCumEnX)

      do j = 1,nt2
           TO1(j) = TO1(j)/TO1(nt2)
           TO2(j) = TO2(j)/TO2(nt2)
           if(TO1(j).le.0.05) DLS = REAL(j-1)*dt2
           if(TO2(j).le.0.05) DLR = REAL(j-1)*dt2
           if(TO1(j).le.0.75) DHS = REAL(j-1)*dt2
           if(TO2(j).le.0.75) DHR = REAL(j-1)*dt2
      enddo
      DHS  = DHS-DLS
      DHR  = DHR-DLR
      call GOF(DHR, DHS, 1, 1, 1, DFX)

C Y component comparison
      call DINTEG(ABS(VSy)**2.0, dt2, nt2, TO1)
      call DINTEG(ABS(VRy)**2.0, dt2, nt2, TO2)
      call GOF(TO1, TO2, nt2, nt2, nt2, DCumEnY)

      do j = 1,nt2
           TO1(j) = TO1(j)/TO1(nt2)
           TO2(j) = TO2(j)/TO2(nt2)
           if(TO1(j).le.0.05) DLS = REAL(j-1)*dt2
           if(TO2(j).le.0.05) DLR = REAL(j-1)*dt2
           if(TO1(j).le.0.75) DHS = REAL(j-1)*dt2
           if(TO2(j).le.0.75) DHR = REAL(j-1)*dt2
      enddo
      DHS  = DHS-DLS
      DHR  = DHR-DLR
      call GOF(DHR, DHS, 1, 1, 1, DFY)

C Z component comparison
      call DINTEG(ABS(VSz)**2.0, dt2, nt2, TO1)
      call DINTEG(ABS(VRz)**2.0, dt2, nt2, TO2)
      call GOF(TO1, TO2, nt2, nt2, nt2, DCumEnZ)

      do j = 1,nt2
           TO1(j) = TO1(j)/TO1(nt2)
           TO2(j) = TO2(j)/TO2(nt2)
           if(TO1(j).le.0.05) DLS = REAL(j-1)*dt2
           if(TO2(j).le.0.05) DLR = REAL(j-1)*dt2
           if(TO1(j).le.0.75) DHS = REAL(j-1)*dt2
           if(TO2(j).le.0.75) DHR = REAL(j-1)*dt2
      enddo
      DHS  = DHS-DLS
      DHR  = DHR-DLR
      call GOF(DHR, DHS, 1, 1, 1, DFZ)
      DCumEn = (DCumEnX+DCumEnY+DCumEnZ)/3.0
      EnDur = (DFX+DFY+DFZ)/3.0
C---------------------------------------------------------------------
C Inelastic/Elastic Fit
      InElFitX = 0.0
      InElFitY = 0.0
      InElFitZ = 0.0
      InElFitH = 0.0
      InElFit  = 0.0
      write(8111,*) '%%%%%%STATION = ',staNum, '%%%%%%'
C     real, dimension(16,19):: GOFIEX, GOFIEY, GOFIEZ, GOFIEH 
      do i = 1,16
           RatCmpX(i) = 0.0
           RatCmpY(i) = 0.0
           RatCmpZ(i) = 0.0
           RatCmpH(i) = 0.0
         do j = 1,19
            call GOF(RatioSX(i,j), RatioRX(i,j), 1, 1, 1, GOFIEX(i,j))
            call GOF(RatioSY(i,j), RatioRY(i,j), 1, 1, 1, GOFIEY(i,j))
            call GOF(RatioSZ(i,j), RatioRZ(i,j), 1, 1, 1, GOFIEZ(i,j))
            GOFIEH(i,j) = (GOFIEX(i,j)+GOFIEY(i,j))/2.0
         enddo
         RatCmpX(i) = MINVAL(GOFIEX(i,:))
         RatCmpY(i) = MINVAL(GOFIEY(i,:))
         RatCmpZ(i) = MINVAL(GOFIEZ(i,:))
         RatCmpH(i) = MINVAL(GOFIEH(i,:))
      enddo

      do i = L16,H16

C      tmpxMax = MAXVAL(ABS(RatioSX(i,:) - RatioRX(i,:)))
C      tmpyMax = MAXVAL(ABS(RatioSY(i,:) - RatioRY(i,:)))
C      tmpzMax = MAXVAL(ABS(RatioSZ(i,:) - RatioRZ(i,:)))
C
C
C     
C      do j = 1,19
C         if(ABS(RatioSX(i,j) - RatioRX(i,j)) .EQ. tmpxMax) then
C            call GOF(RatioSX(i,j), RatioRX(i,j), 1, 1, 1, RatCmpX(i))
C         endif
C         if(ABS(RatioSY(i,j) - RatioRY(i,j)) .EQ. tmpyMax) then
C            call GOF(RatioSY(i,j), RatioRY(i,j), 1, 1, 1, RatCmpY(i))
C         endif
C         if(ABS(RatioSZ(i,j) - RatioRZ(i,j)) .EQ. tmpzMax) then
C            call GOF(RatioSZ(i,j), RatioRZ(i,j), 1, 1, 1, RatCmpZ(i))
C         endif
C      enddo
C
C         call GOF(RatioSX(i,:), RatioRX(i,:), 19, 1, 
C     +        19, RatCmpX(i))
C         call GOF(RatioSY(i,:), RatioRY(i,:), 19, 1, 
C     +        19, RatCmpY(i))
C         call GOF(RatioSZ(i,:), RatioRZ(i,:), 19, 1, 
C     +        19, RatCmpZ(i))

         InElFitX = InElFitX + RatCmpX(i)
         InElFitY = InElFitY + RatCmpY(i)
         InElFitZ = InElFitZ + RatCmpZ(i)
         InElFitH = InElFitH + RatCmpH(i)

      write(8111,*) period(i), 100*RatCmpX(i), 100*RatCmpY(i), 
     +100*RatCmpZ(i)
      enddo
      tempI = (H16-L16)+1
      InElFitX = InElFitX/real(tempI)
      InElFitY = InElFitY/real(tempI)
      InElFitZ = InElFitZ/real(tempI)
      InElFit  = InElFitH/real(tempI)
C---------------------------------------------------------------------
C SAfit16
      temp  = 0.0
      temp1 = 0.0
      temp2 = 0.0
      SAfit16X = 0.0
      SAfit16Y = 0.0
      SAfit16Z = 0.0

      do j = L16,H16
         call GOF(SpecRX(j),SpecSX(j),1,1,1,temp)
         call GOF(SpecRY(j),SpecSY(j),1,1,1,temp1)
         call GOF(SpecRZ(j),SpecSZ(j),1,1,1,temp2)

         SAfit16X = SAfit16X + temp
         SAfit16Y = SAfit16Y + temp1
         SAfit16Z = SAfit16Z + temp2
      enddo

      SAfit16X = SAfit16X/real(tempI)
      SAfit16Y = SAfit16Y/real(tempI)
      SAfit16Z = SAfit16Z/real(tempI)
      SAfit16 = (SAfit16X+SAfit16Y+SAfit16Z)/3.0
      
C Set 1 = Reference
C Set 2 = Synthetics
      do j = 1,16
         call GOF(SpecRX(j),SpecSX(j),1,1,1,temp)
         call GOF(SpecRY(j),SpecSY(j),1,1,1,temp1)
         call GOF(SpecRZ(j),SpecSZ(j),1,1,1,temp2)
	  write(232323,*) period(j), temp, SpecRX(j), SpecSX(j)
	  write(242424,*) period(j), temp1, SpecRY(j), SpecSY(j)
	  write(252525,*) period(j), temp2, SpecRZ(j), SpecSZ(j)
      enddo
C---------------------------------------------------------------------
C SpecDurFit
      temp  = 0.0
      temp1 = 0.0
      temp2 = 0.0
      SpecDurFitX = 0.0
      SpecDurFitY = 0.0
      SpecDurFitZ = 0.0

      do j = L16,H16
         call GOF(DurRX(j),DurSX(j),1,1,1,temp)
         call GOF(DurRY(j),DurSY(j),1,1,1,temp1)
         call GOF(DurRZ(j),DurSZ(j),1,1,1,temp2)

         write(852,*) period(j), 100*temp, 100*temp1, 100*temp2

         SpecDurFitX = SpecDurFitX + temp
         SpecDurFitY = SpecDurFitY + temp1
         SpecDurFitZ = SpecDurFitZ + temp2
      enddo

      SpecDurFitX = SpecDurFitX/real(tempI)
      SpecDurFitY = SpecDurFitY/real(tempI)
      SpecDurFitZ = SpecDurFitZ/real(tempI)
      SpecDurFit  = (SpecDurFitX+SpecDurFitY+SpecDurFitZ)/3.0
C---------------------------------------------------------------------
C Cross-Correlation and Covariance Fit
      temp3 = 0.0
      do k = 1,nt2
         LFSx(k) = SSx2(k)
         LFSy(k) = SSy2(k)
         LFSz(k) = SSz2(k)

         LFRx(k) = SRx2(k)
         LFRy(k) = SRy2(k)
         LFRz(k) = SRz2(k)
      enddo

C X component Cross-Correlation, Delay, Covariance
      call CROSSC(SSx2, SRx2, nt2, dt2, CCX, DLAX)
C      call COVAR(SSx2, SRx2, nt2, CVX)

C Y component Cross-Correlation, Delay, Covariance
      call CROSSC(SSy2, SRy2, nt2, dt2, CCY, DLAY)
C      call COVAR(SSy2, SRy2, nt2, CVY)

C Z component Cross-Correlation, Delay, Covariance
      call CROSSC(SSz2, SRz2, nt2, dt2, CCZ, DLAZ)
C      call COVAR(SSz2, SRz2, nt2, CVZ)

C      call ERF((ABS(CVX)+ABS(CVY)+ABS(CVZ))/3.0, CVfit)
      CCfit = (CCX+CCY+CCZ)/3.0
C---------------------------------------------------------------------
C PPMCC fit
      write(5555,*) DLAX, DLAY, DLAZ
      r1 = 0.0
      r2 = 0.0
      r3 = 0.0

C      call PPMCC(LFRx, LFSx, nt2, r1)
C      call PPMCC(LFRy, LFSy, nt2, r2)
C      call PPMCC(LfRz, LFSz, nt2, r3)

C      if(r1 .LT. 0.0) r1 = 0.0
C      if(r2 .LT. 0.0) r2 = 0.0
C      if(r3 .LT. 0.0) r3 = 0.0 
C      PPMCCfit = (r1+r2+r3)/3.0
C---------------------------------------------------------------------
C Fourier Spectrum Comparison
      call GOFFS(SSx2, SRx2, nt2, dt2, 1.0/CutH, 1.0/CutL, FScompX)
      call GOFFS(SSy2, SRy2, nt2, dt2, 1.0/CutH, 1.0/CutL, FScompY)
      call GOFFS(SSz2, SRz2, nt2, dt2, 1.0/CutH, 1.0/CutL, FScompZ)

      FScomp = (FScompX+FScompY+FScompZ)/3.0
C --------------------------------------------------------------------
C Weighting on PGA
C Weighting on PGV
C Weighting on PGD
C Weighting on PSA
C Weighting on Seismogram Fit (weights must be integers)
C Weighting on Spectral Fit
C Weighting on Cumulative Energy Fit
C Weighting on Inelastic/Elastic Fit (16)
C Weighting on Spec Acc (16)
C Weighting on Spec Dur (16)
C Weighting on Pearson's Product-Moment Correlation Coefficient
C Weighting on Data Energy Release Duration (5%-75%)
C Weighting on Cross-Correlation
C Weighting on Covariance
C --------------------------------------------------------------------
      siteV = 100*(
     +   VAL(1) *PGAfit      +
     +   VAL(2) *PGVfit      +
     +   VAL(3) *PGDfit      +
     +   VAL(4) *PSAfit      +
     +   VAL(5) *SpFitVal    +
     +   VAL(6) *DCumEn      +
     +   VAL(7) *InElFit     +
     +   VAL(8) *SAfit16     +
     +   VAL(9) *SpecDurFit  +
     +   VAL(10)*EnDur       +
     +   VAL(11)*CCfit       +
     +   VAL(12)*FScomp      +
     +   VAL(13)*0.0         +
     +   VAL(14)*0.0         +
     +   VAL(15)*0.0         +
     +   VAL(16)*0.0         +
     +   VAL(17)*0.0         +
     +   VAL(18)*0.0         +
     +   VAL(19)*0.0)        /
     +   VAL(20)
C --------------------------------------------------------------------
      siteVX = 100*(
     +   VAL(1) *PGAfitX      +
     +   VAL(2) *PGVfitX      +
     +   VAL(3) *PGDfitX      +
     +   VAL(4) *PSAfitX      +
     +   VAL(5) *SRCmpX       +
     +   VAL(6) *DCumEnX      +
     +   VAL(7) *InElFitX     +
     +   VAL(8) *SAfit16X     +
     +   VAL(9) *SpecDurFitX  +
     +   VAL(10)*DFX          +
     +   VAL(11)*CCX          +
     +   VAL(12)*FScompX      +
     +   VAL(13)*0.0          +
     +   VAL(14)*0.0          +
     +   VAL(15)*0.0          +
     +   VAL(16)*0.0          +
     +   VAL(17)*0.0          +
     +   VAL(18)*0.0          +
     +   VAL(19)*0.0)         /
     +   VAL(20)
C --------------------------------------------------------------------
      siteVY = 100*(
     +   VAL(1) *PGAfitY      +
     +   VAL(2) *PGVfitY      +
     +   VAL(3) *PGDfitY      +
     +   VAL(4) *PSAfitY      +
     +   VAL(5) *SRCmpY       +
     +   VAL(6) *DCumEnY      +
     +   VAL(7) *InElFitY     +
     +   VAL(8) *SAfit16Y     +
     +   VAL(9)*SpecDurFitY  +
     +   VAL(10)*DFY          +
     +   VAL(11)*CCY          +
     +   VAL(12)*FScompY      +
     +   VAL(13)*0.0          +
     +   VAL(14)*0.0          +
     +   VAL(15)*0.0          +
     +   VAL(16)*0.0          +
     +   VAL(17)*0.0          +
     +   VAL(18)*0.0          +
     +   VAL(19)*0.0)         /
     +   VAL(20)
C---------------------------------------------------------------------
      siteVZ = 100*(
     +   VAL(1) *PGAfitZ      +
     +   VAL(2) *PGVfitZ      +
     +   VAL(3) *PGDfitZ      +
     +   VAL(4) *PSAfitZ      +
     +   VAL(5) *SRCmpZ       +
     +   VAL(6) *DCumEnZ      +
     +   VAL(7) *InElFitZ     +
     +   VAL(8) *SAfit16Z     +
     +   VAL(9)*SpecDurFitZ  +
     +   VAL(10)*DFZ          +
     +   VAL(11)*CCZ          +
     +   VAL(12)*FScompZ      +
     +   VAL(13)*0.0          +
     +   VAL(14)*0.0          +
     +   VAL(15)*0.0          +
     +   VAL(16)*0.0          +
     +   VAL(17)*0.0          +
     +   VAL(18)*0.0          +
     +   VAL(19)*0.0)         /
     +   VAL(20)
C --------------------------------------------------------------------
      write(*,*) '-------------------------------'

      SpFitVal = (100*SpFitVal)
      SRCmpX = (100*SRCmpX)
      SRCmpY = (100*SRCmpY)
      SRCmpZ = (100*SRCmpZ)
      write(*,*) 'Spectral Fit Value             = ', NINT(SpFitVal)

      DCumEn = (100*DCumEn)
      DCumEnX = (100*DCumEnX)
      DCumEnY = (100*DCumEnY)
      DCumEnZ = (100*DCumEnZ)
      write(*,*) 'Integated Velocity Vector      = ', NINT(DCumEn)

      PGAfit = (100*PGAfit)
      PGAfitX = (100*PGAfitX)
      PGAfitY = (100*PGAfitY)
      PGAfitZ = (100*PGAfitZ)
      write(*,*) 'PGA Fit                        = ', NINT(PGAfit)

      PGVfit = (100*PGVfit)
      PGVfitX = (100*PGVfitX)
      PGVfitY = (100*PGVfitY)
      PGVfitZ = (100*PGVfitZ)
      write(*,*) 'PGV Fit                        = ', NINT(PGVfit)

      InElFit = (100*InElFit)
      InElFitX = (100*InElFitX)
      InElFitY = (100*InElFitY)
      InElFitZ = (100*InElFitZ)
      write(*,*) 'Inelastic/Elastic Fit (16)     = ', NINT(InElFit)

      SAfit16 = (100*SAfit16)
      SAfit16X = (100*SAfit16X)
      SAfit16Y = (100*SAfit16Y)
      SAfit16Z = (100*SAfit16Z)
      write(*,*) 'Spectral Acceleration Fit (16) = ', NINT(SAfit16)

      SpecDurFit = (100*SpecDurFit)
      SpecDurFitX = (100*SpecDurFitX)
      SpecDurFitY = (100*SpecDurFitY)
      SpecDurFitZ = (100*SpecDurFitZ)
      write(*,*) 'Spectral Duration Fit (16)     = ', NINT(SpecDurFit)

      CCfit = (100*CCfit)
      CCX = (100*CCX)
      CCY = (100*CCY)
      CCZ = (100*CCZ)
      write(*,*) 'Cross-Correlation              = ', NINT(CCfit)

C      CVfit = CVfit
C      CVX = CVX
C      CVY = CVY
C      CVZ = CVZ

      EnDur = (100*EnDur)
      DFX = (100*DFX)
      DFY = (100*DFY)
      DFZ = (100*DFZ)
      write(*,*) 'Duration                       = ', NINT(EnDur)

      FScomp  = (100*FScomp)
      FScompX = (100*FScompX)
      FScompY = (100*FScompY)
      FScompZ = (100*FScompZ)
      write(*,*) 'Fourier Spectrum Fit           = ', NINT(FScomp)

      write(*,*) 'Site Fit Value                 = ', NINT(siteV)
      write(*,*) 'Site Fit Value (X)             = ', NINT(siteVX)
      write(*,*) 'Site Fit Value (Y)             = ', NINT(siteVY)
      write(*,*) 'Site Fit Value (Z)             = ', NINT(siteVZ)
      write(*,*) '-------------------------------'
C---------------------------------------------------------------------
C      (38,GOF.list)
C      (1111,GOF_PGA.list)
C      (2222,GOF_PGV.list)
C      (3333,GOF_PGD.list)
C      (4444,GOF_PSA.list)
C      (6666,GOF_SPECFIT.list)
C      (7777,GOF_ENERGYFIT.list)
C      (8888,GOF_InElEl.list)
C      (9999,GOF_SAFIT.list)
C      (101010,GOF_PPMCC.list)
C      (111111,GOF_SPECDUR.list)
C      (121212,GOF_CROSSCOR.list)
C      (131313,GOF_COVAR.list)
C      (141414,GOF_DUR.list)
C      (151515,GOF_FS.list)
C      (171717,GOF_NGA.list)
C---------------------------------------------------------------------
      write(38,*) siteV, siteVX, siteVY, siteVZ
      write(1111,*) PGAfit, PGAfitX, PGAfitY, PGAfitZ
      write(2222,*) PGVfit, PGVfitX, PGVfitY, PGVfitZ

      PGDfit = (100*PGDfit)
      PGDfitX = (100*PGDfitX)
      PGDfitY = (100*PGDfitY)
      PGDfitZ = (100*PGDfitZ)
      write(3333,*) PGDfit, PGDfitX, PGDfitY, PGDfitZ

      PSAfit = (100*PSAfit)
      PSAfitX = (100*PSAfitX)
      PSAfitY = (100*PSAfitY)
      PSAfitZ = (100*PSAfitZ)
      write(4444,*) PSAfit, PSAfitX, PSAfitY, PSAfitZ

      write(6666,*) SpFitVal, SRCmpX, SRCmpY, SRCmpZ
      write(7777,*) DCumEn, DCumEnX, DCumEnY, DCumEnZ
      write(8888,*) InElFit, InElFitX, InElFitY, InElFitZ
      write(9999,*) SAfit16, SAfit16X, SAfit16Y, SAfit16Z
      write(111111,*) SpecDurFit,SpecDurFitX,SpecDurFitY,SpecDurFitZ
      write(121212,*) CCfit, CCX, CCY, CCZ
C      write(131313,*) CVfit, CVX, CVY, CVZ
      write(141414,*) EnDur, DFX, DFY, DFZ
      write(151515,*) FScomp, FScompX, FScompY, FScompZ
C---------------------------------------------------------------------

      write(357,*) '%------------------------------------'
      write(357,*) '%Station ', staNum 
      write(357,*) '%----------------'
      write(357,*) '%-------------------------------'
      write(357,*) NINT(SRCmpX), '   %Spectral Fit Value X' 
      write(357,*) NINT(DCumEnX), 
     +'%GOF on Integated Velocity Vector X' 
      write(357,*) NINT(PGAfitX),  '   %PGA Fit X' 
      write(357,*) NINT(PGVfitX),  '   %PGV Fit X' 
      write(357,*) NINT(PGDfitX),  '   %PGD Fit X' 
      write(357,*) NINT(InElFitX), '   %Inelastic/Elastic Fit (16) X' 
      write(357,*) NINT(SAfit16X), 
     +'   %Spectral Acceleration Fit (16) X' 
      write(357,*) NINT(SpecDurFitX), 
     +'   %Spectral Duration Fit (16) X' 
      write(357,*) NINT(CCX),    '   %Cross-Correlation X' 
      write(357,*) NINT(DFX),      '   %Duration X' 
      write(357,*) '%-------------------------------' 
      write(357,*) NINT(SRCmpY), '   %Spectral Fit Value Y' 
      write(357,*) NINT(DCumEnY), 
     +'%GOF on Integated Velocity Vector Y' 
      write(357,*) NINT(PGAfitY),  '   %PGA Fit Y' 
      write(357,*) NINT(PGVfitY),  '   %PGV Fit Y' 
      write(357,*) NINT(PGDfitY),  '   %PGD Fit Y' 
      write(357,*) NINT(InElFitY), '   %Inelastic/Elastic Fit (16) Y' 
      write(357,*) NINT(SAfit16Y), 
     +'   %Spectral Acceleration Fit (16) Y' 
      write(357,*) NINT(SpecDurFitY), 
     +'   %Spectral Duration Fit (16) Y' 
      write(357,*) NINT(CCY),    '   %Cross-Correlation Y' 
      write(357,*) NINT(DFY),      '   %Duration Y' 
      write(357,*) '%-------------------------------' 
      write(357,*) NINT(SRCmpZ), '   %Spectral Fit Value Z' 
      write(357,*) NINT(DCumEnZ), 
     +'%GOF on Integated Velocity Vector Z' 
      write(357,*) NINT(PGAfitZ),  '   %PGA Fit Z' 
      write(357,*) NINT(PGVfitZ),  '   %PGV Fit Z' 
      write(357,*) NINT(PGDfitZ),  '   %PGD Fit Z' 
      write(357,*) NINT(InElFitZ), '   %Inelastic/Elastic Fit (16) Z' 
      write(357,*) NINT(SAfit16Z), 
     +'   %Spectral Acceleration Fit (16) Z' 
      write(357,*) NINT(SpecDurFitZ), 
     +'   %Spectral Duration Fit (16) Z'  
      write(357,*) NINT(CCZ),    '   %Cross-Correlation Z' 
      write(357,*) NINT(DFZ),      '   %Duration Z' 
      write(357,*) '%-------------------------------'
      write(357,*) NINT(siteV),      '   %Site Fit Value' 
      write(357,*) NINT(siteVX),     '   %Site Fit Value (X)' 
      write(357,*) NINT(siteVY),     '   %Site Fit Value (Y)' 
      write(357,*) NINT(siteVZ),     '   %Site Fit Value (Z)' 
      write(357,*) '%-------------------------------'
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C calculate Output values
C SpecSX(i), DurSX(i), SpDisSX(i)
C---------------------------------------------------------------------
      write(357,*) '% PGA -- Set 1 -- Set 2 (horizontal)' 
      write(357,*) '%', PGAR(staNum,1), PGAS(staNum,1)
      write(357,*) '% PGV -- Set 1 -- Set 2 (horizontal)' 
      write(357,*) '%', PGVR(staNum,1), PGVS(staNum,1) 
      write(357,*) '% PGD -- Set 1 -- Set 2 (horizontal)' 
      write(357,*) '%', PGDR(staNum,1), PGDS(staNum,1)
      write(357,*) '% PSA -- Period -- Set 1' 
      write(357,*) '%', TR, PSAR
      write(357,*) '% PSA -- Period -- Set 2' 
      write(357,*) '%', TS, PSAS

C Spectrum Fit (16) values
      do i = 1,16
          call GOF(SpecRX(i), SpecSX(i), 1, 1, 1, temp)
          call GOF(DurRX(i), DurSX(i), 1, 1, 1, temp1)
          SpFX(i) = (RatCmpX(i)+temp+(temp1))/3

          call GOF(SpecRY(i), SpecSY(i), 1, 1, 1, temp)
          call GOF(DurRY(i), DurSY(i), 1, 1, 1, temp1)
          SpFY(i) = (RatCmpY(i)+temp+(temp1))/3

          call GOF(SpecRZ(i), SpecSZ(i), 1, 1, 1, temp)
          call GOF(DurRZ(i), DurSZ(i), 1, 1, 1, temp1)
          SpFZ(i) = (RatCmpZ(i)+temp+(temp1))/3

      enddo

C Best and Worst Spectral Fit
      call organize(SpFX, OUTPUT2, 16)
      do i = 1,16
         if(OUTPUT2(i) .EQ. 16) then
            BSFx(1) = period(i)
            BSFx(2) = SpFX(i)
         elseif(OUTPUT2(i) .EQ. 1) then
            WSFx(1) = period(i)
            WSFx(2) = SpFX(i)
         endif
      enddo

      call organize(SpFY, OUTPUT2, 16)
      do i = 1,16
         if(OUTPUT2(i) .EQ. 16) then
            BSFy(1) = period(i)
            BSFy(2) = SpFY(i)
         elseif(OUTPUT2(i) .EQ. 1) then
            WSFy(1) = period(i)
            WSFy(2) = SpFY(i)
         endif
      enddo

      call organize(SpFZ, OUTPUT2, 16)
      do i = 1,16
         if(OUTPUT2(i) .EQ. 16) then
            BSFz(1) = period(i)
            BSFz(2) = SpFZ(i)
         elseif(OUTPUT2(i) .EQ. 1) then
            WSFz(1) = period(i)
            WSFz(2) = SpFZ(i)
         endif
      enddo


      write(357,*) '%Spectrum Fit X (period, BSF, period, WSF)'
      write(357,*) '%', BSFx(1), BSFx(2), WSFx(1), WSFx(2)
      write(357,*) '%Spectrum Fit Y (period, BSF, period, WSF)'
      write(357,*) '%', BSFy(1), BSFy(2), WSFy(1), WSFy(2)
      write(357,*) '%Spectrum Fit Z (period, BSF, period, WSF)'
      write(357,*) '%', BSFz(1), BSFz(2), WSFz(1), WSFz(2)
      write(357,*) ''
      write(357,*) ''
C---------------------Output Declarations-----------------------------
C Generate the output files
         do i = 1,991
         write(86, *) SDSX(i)
         write(87, *) SDSY(i)
         write(88, *) HSRS(i)

         write(93, *) SDSX(i)
         write(94, *) SDSY(i)
         write(95, *) SDRX(i)
         write(96, *) SDRY(i)

         write(971, *) SRSX(i)
         write(981, *) SRSY(i)
         write(972, *) SRRX(i)
         write(982, *) SRRY(i)

         write(97, *) HSRS(i)
         write(98, *) HSRR(i)
         enddo
C---------------------------------------------------------------------
C BB_GOF_Hartzell
C Model Bias (Abramson) "Uncertainty in Numerical Strong Ground motions .."
C---------------------------------------------------------------------
      do k = 1,991
          rsqSRR(staNum,k) = (SQRT(SRRX(k)*SRRY(k)))
          rsqSRS(staNum,k) = (SQRT(SRSX(k)*SRSY(k)))
          rsqSDR(staNum,k) = (SQRT(SDRX(k)*SDRY(k)))
          rsqSDS(staNum,k) = (SQRT(SDSX(k)*SDSY(k)))
      enddo
C---------------------------------------------------------------------
C--------------------END STATION CALCULATION--------------------------
C---------------------------------------------------------------------
      enddo
C---------------------------------------------------------------------
C---------------------------------------------------------------------
      call MFile(numSta, dt2, nt2, CutL, CutH)
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C BB_GOF_Hartzell
C Model Bias (Abramson) "Uncertainty in Numerical Strong Ground motions .."
      call MODBIAS(rsqSRR,rsqSRS,numSta,991,SRmBIAS)
      call MODBIAS(rsqSDR,rsqSDS,numSta,991,SDmBIAS)
      call MODERROR(rsqSRR,rsqSRS,numSta,991,SRmBIAS,SRmERROR)
      call MODERROR(rsqSDR,rsqSDS,numSta,991,SDmBIAS,SDmERROR)

      call MODBIAS(PGAR,PGAS,numSta,1,PGAmBIAS)
      call MODERROR(PGAR,PGAS,numSta,1,PGAmBIAS,PGAmERROR)    

      call MODBIAS(PGVR,PGVS,numSta,1,PGVmBIAS)
      call MODERROR(PGVR,PGVS,numSta,1,PGVmBIAS,PGVmERROR)     

      write(98277,*) '% PGA - Error - Bias'
      write(98277,*) PGAmERROR, PGAmBIAS
      write(98277,*) '% PGV - Error - Bias'
      write(98277,*) PGVmERROR, PGVmBIAS
      write(98277,*) '% Spec Resp - Error - Bias'
      do i = 1,991
          write(98277,*) SRmERROR(i), SRmBIAS(i)
      enddo
      write(98277,*) '% Spec Dur - Error - Bias'
      do i = 1,991
          write(98277,*) SDmERROR(i), SDmBIAS(i)
      enddo
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C---------------------------------------------------------------------

C Close the output files
         close(86)
         close(87)
         close(88)
         close(89)
         close(93)
         close(94)
         close(95)
         close(96)

         close(971)
         close(981)
         close(972)
         close(982)

         close(97)
         close(98)
         close(991)
         close(992)
         close(9875)
         close(9876)

      close(357)
      close(38)
      close(1111)
      close(2222)
      close(3333)
      close(4444)
      close(6666)
      close(7777)
      close(8888)
      close(9999)
      close(101010)
      close(111111)
      close(121212)
      close(131313)
      close(141414)
      close(151515)
      close(232323)
      close(242424)
      close(252525)
C---------------------------------------------------------------------
      end subroutine MAIN_COMP
C --------------------------------------------------------------------
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                    C
C     EEEEEEEEEEE  NNNN         NN  DDDDDD           C
C     EE           NNNNN        NN  DD   DDD         C
C     EE           NN  NN       NN  DD    DDD        C
C     EE           NN   NN      NN  DD      DD       C
C     EEEEEE       NN    NN     NN  DD      DD       C
C     EEEEEE       NN     NN    NN  DD      DD       C
C     EE           NN      NN   NN  DD      DD       C
C     EE           NN       NN  NN  DD    DDD        C
C     EE           NN        NNNNN  DD   DDD         C
C     EEEEEEEEEEE  NN         NNNN  DDDDDD           C
C                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  NOTES
C---------------------------------------------------------------------
C  MAIN_COMP(inputS, inputR, dtS, dtR, dt2, ntS, ntR, nt2, leng, M, del, lam, Ztor, TAG, numSta, VAL)
C  inputS        = Set 2 seismogram file name
C  inputR        = Set 1 seismogram file name
C  dtS           = width of timesteps Set 2
C  dtR           = width of timesteps Set 1
C  dt2           = width of interpolated timesteps
C  ntS           = number of timesteps in Set 2 seismogram
C  ntR           = number of timesteps in Set 1 seismogram
C  nt2           = number of timesteps in interpolated seismogram
C  leng          = length of the seismograms
C  M             = Magnitude
C  del           = dip of the fault
C  lam           = rake of fault displacement 
C  Ztor          = distance to the surface from the top of the fault rupture
C  TAG           = A for acceleration, V for velocity, D for displacement
C  numSta        = number of stations to be analyzed
C  VAL           = 
C---------------------------------------------------------------------
C  CB_2007_NGA(M, T, Rjb, Rrup, Ztor, Zvs, Vs, lambda, delta, arb, y, sigma)
C  M             = Magnitude
C  T             = Period (sec); Use -1 for PGV computation and -10 for PGD computation
C  Rjb           = Joyner-Boore distance (km)
C  Rrup          = Closest distance coseismic rupture (km)
C  Ztor          = Depth to the top of coseismic rupture (km)
C  Zvs           = Depth to the 2.5 km/s shear-wave velocity horizon (km)
C  Vs            = shear wave velocity averaged over top 30 m (m/s)---(Vs30)
C  lambda        = rake angle (degree)
C  delta         = average dip of the rupture place (degree)
C  arb           = 1 for arbitrary component sigma; 0 for average component sigma
C  y             = Median spectral acceleration prediction
C  sigma         = logarithmic standard deviation of spectral acceleration prediction
C---------------------------------------------------------------------
C  FileRead(input, Vx, Vy, Vz, dt, nt, xx, yy, Rjb, Rrup, ZvsISO, Vs30, staNum)
C  input         = seismogram file name
C  Vx            = East component seismogram
C  Vy            = North compontent seismogram
C  Vz            = vertical component seismogram
C  dt            = width of time step interval
C  nt            = number of time steps
C  xx            = x location of the station
C  yy            = y location of the station
C  Rjb           = Joyner-Boore distance (km)
C  Rrup          = Closest distance coseismic rupture (km) 
C  ZvsISO        = Depth to the 2.5 km/s iso-surface
C  Vs30          = Shear Wave velocity in top 30 meters
C  staNum        = Station Number
C---------------------------------------------------------------------
C  GetParam(inputS, inputR, ntS, ntR, leng, M, del, lam, Ztor, TAG, numSta)
C  inputS        = seismogram file name (Set 2)
C  inputR        = seismogram file name (Set 1)
C  ntS           = width of time step interval
C  ntR           = number of time steps
C  leng          = length of seismograms (seconds)
C  M             = Magnitude
C  del           = dip of the fault
C  lam           = rake of fault displacement 
C  Ztor          = distance to the surface from the top of the fault rupture
C  TAG           = A for acceleration, V for velocity, D for displacement
C  numSta        = number of stations to be analyzed
C  VAL           = list of weighting values
C---------------------------------------------------------------------
C  MFile(numSta, dt2, nt2)
C  a subroutine to generate a matlab script for easy plotting
C  numSta        = number of stations in the analysis
C  dt2           = width of time step (secs)
C  nt2           = number of timesteps
C---------------------------------------------------------------------
C  FILT(seism, dt2, nt2, output, f_low, f_high)
C  seism         = seismogram to be filtered
C  dt2           = width of timesteps
C  nt2           = number of timesteps in seismogram
C  output        = filtered seismogram
C  f_low         = lowest frequency to pass
C  f_high        = highest frequency to pass 
C---------------------------------------------------------------------
C  DINTEG(seism, dt2, nt2, output)
C  seism         = seismogram to be integrated
C  dt2           = width of timesteps
C  nt2           = number of timesteps in seismogram
C  output        = integrated seismogram
C---------------------------------------------------------------------
C  DIFF(seism, dt2, nt2, output)
C  seism         = seismogram to be differentiated
C  dt2           = width of timesteps
C  nt2           = number of timesteps in seismogram
C  output        = differentiated seismogram
C---------------------------------------------------------------------
C  FFT(seismR, seismI, dt2, nt2, outR, outI, flag)
C  seismR        = real part of seismogram to be transformed
C  seismI        = imaginary part of seismogram to be transformed
C  dt2           = width of timesteps
C  nt2           = number of timesteps in seismogram
C  outR          = real portion of transformed seismogram
C  outI          = imaginary portion of transformed seismogram
C  flag          = (1 = FFT) (-1 = IFFT)
C---------------------------------------------------------------------
C  RESPONSESPEC(acc, dt, nt, period, spec, specDis)
C  acc           = accelerogram to be analyzed
C  dt            = width of timesteps
C  nt            = number of timesteps in seismogram
C  period        = period to be evaluated
C  spec          = value of spectral acceleration at the given period
C  specDur       = spectral duration
C---------------------------------------------------------------------
C  InElastic( T, dy, alpha, ag, dt, nt, maxA)
C  ag            = accelerogram to be analyzed
C  dt            = width of timesteps
C  nt            = number of timesteps in seismogram
C  T             = period to be evaluated
C  dy             = yield displacements
C  alpha          = strain hardening ratios 
C  maxD          = inelastic spectral displacement
C---------------------------------------------------------------------
C  GOF(IN1, IN2, nt, nts, ntf, SSQ)
C  IN1           = seismogram1 to be compared
C  IN2           = seismogram2 to be compared
C  nt            = number of timesteps in IN1
C  nts           = starting timestep
C  ntf           = finishing timestep
C  SSQ           = GOF measure
C---------------------------------------------------------------------
C  GOFSA(IN1, IN2, nt, nts, ntf, SSQ)
C  IN1           = SA seismogram1 to be compared
C  IN2           = SA seismogram2 to be compared
C  nt            = number of timesteps in IN1
C  nts           = starting timestep
C  ntf           = finishing timestep
C  SSQ           = GOF measure
C---------------------------------------------------------------------
C  GOFFS(INS, INR, nt, dt, freq1, freq2, SSQ)
C  IN1           = SA seismogram1 to be compared
C  IN2           = SA seismogram2 to be compared
C  nt            = number of timesteps
C  nt            = width of the timesteps
C  freq1         = starting frequency
C  freq2         = finishing frequency
C  SSQ           = GOF measure
C---------------------------------------------------------------------
C  NORM(IN1, nt, nts, ntf, nm, pwr)
C  IN1           = seismogram to be interpolated
C  nt            = number of timesteps in IN1
C  nts           = starting timestep
C  ntf           = finishing timestep
C  nm            = norm value
C  pwr           = power of the norm computation -- L(pwr) norm
C---------------------------------------------------------------------
C  SpecComp(SPECS, SPECR, COMPMDEV)
C  SPECS         = Spectral response of the synthetic seismogram (length = 991)
C  SPECR         = Spectral response of the raw seismogram (length = 991)
C  COMPMDEV      = vector with the 12 difference means from the spectral 
C                  simplification comparison
C---------------------------------------------------------------------
C  NGAComp(NGA, SPEC, COMPMDEV)
C  NGA           = Spectral response of the NGA
C  SPEC          = Spectral response of the seismogram (length = 991)
C  COMPMDEV      = vector with the 16 difference means from the spectral 
C                  simplification comparison
C---------------------------------------------------------------------
C  organize(INPUT, OUTPUT2, num)
C  INPUT         = A vector of real values
C  OUTPUT2       = A vector of integer values that give the placeholding value
C                  (sequencially, largest value = 1) for the input vector
C  num           = the total number of placeholdings in the INPUT vector
C---------------------------------------------------------------------
C  MODBIAS(INPUTR, INPUTS, numSta, numPts, B)
C  INPUTS        = Set 2 values
C  INPUTR        = Set 1 values
C  numSta        = The number of independent vectors in the model
C  numPts        = The number of points in each of the model vectors
C  B             = The value of the model bias
C---------------------------------------------------------------------
C  MODERROR(INPUTR, INPUTS, numSta, numPts, B, E)
C  INPUTS        = Set 2 values
C  INPUTR        = Set 1 values
C  numSta        = The number of independent vectors in the model
C  numPts        = The number of points in each of the model vectors
C  B             = The value of the model bias
C  E             = The value of the model error
C---------------------------------------------------------------------
C  ERF(INPUT, OUTPUT)
C  INPUT         = normalized residual value (0 = perfect fit)
C  OUTPUT        = output value generated by the error function (1 = perfect fit)
C---------------------------------------------------------------------
C  CROSSC(IN1, IN2, nt, dt, R, D)
C  IN1           = 1st input vector
C  IN2           = 2nd input vector
C  nt            = number of places in the input vectors
C  R             = Cross-Correlation value
C  D             = The delay value for the best fit (absolute value fit)
C---------------------------------------------------------------------
C  COVAR(IN1, IN2, nt, CV)
C  IN1           = 1st input vector
C  IN2           = 2nd input vector
C  nt            = number of places in the input vectors
C  CV            = Covariance value
C---------------------------------------------------------------------
C  ReSize(seism, output, dt, dt2, nt, nt2)
C  seism         = seismogram to be interpolated
C  output        = interpolated seismogram
C  dt            = width of timesteps
C  dt2           = width of interpolated timesteps
C  nt            = number of timesteps in seismogram
C  nt2           = number of timesteps in interpolated seismogram
C---------------------------------------------------------------------
C  PPMCC(XX, YY, nt, rr)
C  XX            = seismogram 1
C  YY            = seismogram 2
C  nt            = number of timesteps
C  rr            = Pearson Product-Moment Correlation Coefficiant
C---------------------------------------------------------------------
C  SRCC(XX, YY, nt, rho)
C  XX            = seismogram 1
C  YY            = seismogram 2
C  nt            = number of timesteps
C  rho           = Spearman's Rank Correlation Coefficiant
C---------------------------------------------------------------------
C  SDOTP(XX, YY, nt, CosT)
C  XX            = seismogram 1
C  YY            = seismogram 2
C  nt            = number of timesteps
C  CosT          = Cos of the angle between the two vectors
C---------------------------------------------------------------------
C---------------------------------------------------------------------
