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

         open(3030, file='out/ElD_R.out', status='unknown')
         open(3131, file='out/ElD_S.out', status='unknown')

         open(991, file='out/InElS.19.16.NumSta.out', status='unknown')
         open(992, file='out/InElR.19.16.NumSta.out', status='unknown')

      open(9875, file='out/FileS.out', status='unknown')
      open(9876, file='out/FileR.out', status='unknown')
      open(357, file='GOF.log', status='unknown')


      open(8888, file='out/GOF_InElEl.list', status='unknown')
      open(8111, file='out/GOF_InElEl_All.list', status='unknown')

C---------------------------------------------------------------------
      write(357,*) '% Set 1 = ', inputR(1:20)
      write(357,*) '% Set 2 = ', inputS(1:20)
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
C--------------------END OF METRICS CALCULATIONS----------------------
C---------------------------------------------------------------------
      write(*,*) 'Calculating GOF'
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

      write(*,*) '-------------------------------'

      InElFit = (100*InElFit)
      InElFitX = (100*InElFitX)
      InElFitY = (100*InElFitY)
      InElFitZ = (100*InElFitZ)
      write(*,*) 'Inelastic/Elastic Fit (16)     = ', NINT(InElFit)

C      write(*,*) 'Site Fit Value                 = ', NINT(siteV)
C      write(*,*) 'Site Fit Value (X)             = ', NINT(siteVX)
C      write(*,*) 'Site Fit Value (Y)             = ', NINT(siteVY)
C      write(*,*) 'Site Fit Value (Z)             = ', NINT(siteVZ)
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
      write(8888,*) InElFit, InElFitX, InElFitY, InElFitZ

C---------------------------------------------------------------------

      write(357,*) '%------------------------------------'
      write(357,*) '%Station ', staNum 
      write(357,*) '%----------------'
      write(357,*) '%-------------------------------'
      write(357,*) NINT(InElFitX), '   %Inelastic/Elastic Fit (16) X' 
      write(357,*) '%-------------------------------' 
      write(357,*) NINT(InElFitY), '   %Inelastic/Elastic Fit (16) Y' 
      write(357,*) '%-------------------------------' 
      write(357,*) NINT(InElFitZ), '   %Inelastic/Elastic Fit (16) Z' 
      write(357,*) '%-------------------------------'
C---------------------------------------------------------------------
C---------------------------------------------------------------------

C---------------------------------------------------------------------
C--------------------END STATION CALCULATION--------------------------
C---------------------------------------------------------------------
      enddo
C---------------------------------------------------------------------
C---------------------------------------------------------------------

C Close the output files

         close(991)
         close(992)
         close(3030)
         close(3131)
         close(9875)
         close(9876)

      close(357)
      close(8888)
      close(8111)
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
