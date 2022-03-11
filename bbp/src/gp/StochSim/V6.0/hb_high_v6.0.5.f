c hb_high_v5.0
c Needs to be compiled with 'gfortran' (won't work with g77)
c
c 09/23/2010 RWG:
c
c  1) Modifications made to be compatible with 'gfortran'
c
c     - Changed random number generator from F77 rand() to F90 random_number()
c       Calls to rand() have been replaced by function rand_numb(). Seed initialization
c       is done via init_random_seed().
c
c     - Defined variable 'dgamm' in stoc_f() as real*8
c
c  2) Added option to modify velocity profile with normal random perturbations, see grand.
c
c  3) Added option to modify Fourier amplitude spectrum with normal random perturbations, see fgrand.
c

c 09/08/2010 RWG:
c
c  1) Added 'dipfac' parameter to adjust corner frequency, but it is tapered
c     to zero between dmax and dmin. Thus it only applies in full to deeper
c     portion (> dmax) of fault.

c 03/18/2009 RWG:
c
c  1) Changed normalization of random time sequence in stoc_f() such that the average
c     of the POWER spectrum is unity (not amplitude spectrum), this is consistent
c     with Boore (1983).
c  2) Above change reduced average level of motions about 10%.  This is countered
c     by increasing the default corner frequency by 5%.  Thus, zz_default = 2.1.

c hb_high_v3.0
c
c Fall 2008 RWG:
c
c  1) Removed Fc adjustment paramters 'depfac' and 'dipfac'.
c  2) Depth scaling of Fc is now accomodated through scaling of rupture velocity
c     so that it is completely consistent with the low frequency model (and rupture).
c     New parameters are:
c
c       rvfac      = rupture velocity fraction of the local shear velocity; default=0.8
c       shal_rvfac = additional scaling of rupture velocity fraction in the
c                    shallow part of rupture; above dmin, rvf=rvfac*shal_rvfac;
c                    below dmax, rvf=rvfac; linear scaling between dmax and dmin; default=0.7
c
c  3) Set default values for dmin and dmax to
c        dmin_default = 5.0
c        dmax_default = 8.0
c     which is consistent with low frequency specification of rupture velocity and rise time
c     depth scaling.

c hb_high_v2.3
c
c  12/12/05 RWG:
c
c  1) Changed default for 'dipfac' to dipfac=0.0.
c  2) Added additional factor 'fcfac' to specify arbitrary adjustment
c     to corner frequency.  Default is fcfac=0.0

c hb_high_v2.2
c
c  ??/??/05 RWG:
c
c  Added 'depfac' and 'dipfac' as input parameters to specify the
c  amount of adjustment to subfault corner frequency (fc).  By default
c  depfac=0.4 and dipfac=0.0

c hb_high_v2.11-rp
c
c  03/19/04 RWG:
c
c  Discovered that version 2.02 of HF code (hb_high_v2.02) computed
c  average radiation coefficient incorrectly.  Problem is fixed in
c  version 2.03 (hb_high_v2.03).  See "RWG RADPAT FIX" in RADFRQ_lin().

c  12/21/2004 preserve sign of radiation pattern coefficient;
c  see RADFRQ_lin() and highcor_f().

c  multiple-segment fault

c Seismic moment for the small event SMOE
c is calculated based on the average stress on the fault

c       smoe=stress_average*dlm*dx*dw*1.e+21

c  its corner frequency FCE depends on the subfault dimensions and Z factor
c  where z=1.68 for w2 model and 2.67 w3 model

c      fce=zz*vr/dlm/pai

      include 'params.h'
CC      parameter (mm=16348,mmv=30000,mst=6,lv=1000)
CC      parameter (nq=600,np=50)
      DIMENSION RNA(mmv),RNB(mmv),RDNA(mm),DS(3,mmv)                       
      dimension sddp(lv,nq,np),rist(lv,nq,np),rupt(lv,nq,np)
      dimension rlsu(nq,np),phsu(nq,np),thsu(nq,np),dst(nq,np),zet(nq,np),twin(nq,np)
      DIMENSION stdd(mmv,3),dfr(mm),fgrand(mmv),fgrands(mm)
      complex cs(mm,3),sgc

      dimension fn(nlaymax),an(nlaymax)
      dimension rdur(50),dpth(50), dpdr(50)

CFFF Add variable astopq-> This is the along strike distance (km) on
CFFF the top edge of the fault of the reference point (elonq,elatq).
CFFF This replaces ntopq in even_dist1
      dimension strq(lv),dipq(lv),rakeq(lv),elonq(lv),elatq(lv)
      dimension astop(lv),dtop(lv),shyp(lv),dhyp(lv),nx(lv),nw(lv),dx(lv),dw(lv)

      dimension ic(2),c(4),irtype(100)
      real*8 depth,thic,vp,vsh,rho,vp0,vsh0,rho0,vsmoho
      real *4 stime,rpath,p0,qbar,qfexp
      real*4 qp,qs,qs0,grand
                                                      
      CHARACTER*256 asite,slip_model,outname,outname1,outname2,outname3,velfile
      CHARACTER*256 velname,velname1
      CHARACTER*256 dummy
      integer head_lines, j

      character*64 cap,cname(3) 
      
      COMMON/RANDO/ifu 
      common/vmod_in/depth0(nlaymax),thic0(nlaymax),vp0(nlaymax),vsh0(nlaymax),rho0(nlaymax),qp0(nlaymax),qs0(nlaymax),grand(3000)
      common/vmod/depth(nlaymax),thic(nlaymax),vp(nlaymax),vsh(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
                                                              
      pu =3.1415926/180 
      pai=3.1415926 
      IRD=2
      IMDL=2 
      krc = 4

c set time window envelope peak to be at 20% of total duration (e.g. Boore)
      tw_eps = 0.4
      tw_eta = 0.05
      tw_eps = 0.2
      tw_eta = 0.05

c set default factors for scaling of fc
c zz_default=2.1 works when random time sequence spectra is normalized according to Boore (1983).
      zz_default = 2.1

c 20131120 RWG: changed c0 to 2.0 based on SWUS validations
      zz_default = 2.0

c 20171219 RWG: changed c0 to 2.1 and added Czero input parameter replacing fcfac
c    zz=Czero, Czero is consistent with terminology used in papers.  fcfac is now
c    defunct, just hardwired to be =0.0
      Czero_default = 2.1

c 20180822 RWG: changed default Czero to 2.0, being tested ...
      Czero_default = 2.0

      rvfac_default = 0.8

      shal_rvfac_default = 0.6
      shal_dmin_default = 5.0
      shal_dmax_default = 8.0

cc20140422 RWG
      deep_rvfac_default = 0.6
      deep_dmin_default = 15.0
      deep_dmax_default = 20.0

      Calpha_default = 0.1
      fcfac_default = 0.0

      nr = 1000    

      delay=0.
      tvmin=1000.  

      nsfac = 20

      fn( 1) =  0.01
      fn( 2) =  0.02
      fn( 3) =  0.03
      fn( 4) =  0.05
      fn( 5) =  0.07
      fn( 6) =  0.10
      fn( 7) =  0.20
      fn( 8) =  0.30
      fn( 9) =  0.50
      fn(10) =  0.70
      fn(11) =  1.00
      fn(12) =  2.00
      fn(13) =  3.00
      fn(14) =  5.00
      fn(15) =  7.00
      fn(16) = 10.00
      fn(17) = 20.00
      fn(18) = 30.00
      fn(19) = 50.00
      fn(20) = 70.00

      do 4390 i=1,nsfac
         fn(i) = alog(fn(i))
4390  continue

c
c CMP =-90   for E-W comp.
c CMP =0     for N-S

      read(5,*)stress_average
      print*,'stress_average= ',stress_average

      read(5,'(a256)') asite 
      read(5,'(a256)') outname                          
      read(5,*) nrtyp,(irtype(i),i=1,nrtyp)
      read(5,*) isite_amp
      read(5,*) nbu,iftt,flol,fhil   
      read(5,*) irand                                       
      read(5,*) nsite
      read(5,*) duration,dt,fmx,akapp,qfexp

      read(5,*) rvfac,shal_rvfac,deep_rvfac,Czero,Calpha

cRWG  apply default values if input is negative
      if(rvfac.lt.-1.0) rvfac = rvfac_default
      if(shal_rvfac.lt.-1.0) shal_rvfac = shal_rvfac_default
      if(deep_rvfac.lt.-1.0) deep_rvfac = deep_rvfac_default
      if(Czero.lt.-1.0) Czero = Czero_default
      if(Calpha.lt.-1.0) Calpha = Calpha_default

      fcfac = fcfac_default
      fcfac = 0.0

      read(5,*) sm,vr

      read(5,'(a256)') slip_model  
      open(10,file=slip_model)

      read(10,*) nevnt

      zhyp_max = 0.0;
      nstot = 0
      farea_in = 0.0
      do 133 iv=1,nevnt

         read(10,*) elonq(iv),elatq(iv),nx(iv),nw(iv),dx(iv),dw(iv)
         read(10,*) strq(iv),dipq(iv),rakeq(iv),dtop(iv),shyp(iv),dhyp(iv)

	print*,'nx,nw= ',nx(iv),nw(iv)

	 nstot = nx(iv)*nw(iv) + nstot
	 farea_in = nx(iv)*dx(iv)*nw(iv)*dw(iv) + farea_in
	 astop(iv) = 0.5*nx(iv)*dx(iv)

c RWG 201241219: '*' should have been '/' in below, fixed 21041219
	 zhyp = dtop(iv) + dhyp(iv)/sin(dipq(iv)*pu)
	 if(zhyp.gt.zhyp_max) zhyp_max = zhyp

         do j=1,nw(iv)
            read(10,*) (sddp(iv,i,j),i=1,nx(iv))
         end do

         do j=1,nw(iv)
            read(10,*) (rist(iv,i,j),i=1,nx(iv))
         end do

         do j=1,nw(iv)
            read(10,*) (rupt(iv,i,j),i=1,nx(iv))
         end do

133   continue 

      close(10)

	print*,'slip model done'

cc20140422 RWG
      if(zhyp_max.gt.deep_dmin_default) then
         deep_dmin_default = zhyp_max
         deep_dmax_default = deep_dmin_default + 5.0
      endif
            print*,'deep_dmin,deep_dmax,deep_rvfac= ',deep_dmin_default,deep_dmax_default,deep_rvfac

      do iv=1,nevnt
         if(dx(iv).ne.dx(1)) then
            print*,'dx(',iv,')=',dx(iv),' not equal to dx(1)=',dx(1),', exiting...'
	    stop
	 endif
         if(dw(iv).ne.dw(1)) then
            print*,'dw(',iv,')=',dw(iv),' not equal to dw(1)=',dw(1),', exiting...'
	    stop
	 endif
      enddo
      
CCC New way read velocity model from file

      read(5,'(a256)') velfile

C RWG 20140313
c  Option to specify Vs of layer just beneath Moho in case input velocity model
c  extends into upper mantle.  Downgoing rays are reflected from top of last layer
c  in model, so all layers beneath the layer with Vs>vsmoho are truncated.

      vsmoho = -1.0
      read(5,*) vsmoho
      if(vsmoho.le.0.0) vsmoho=999.9

      open(10,file=velfile)
      read(10,*) j0

      jmoho = j0
      do 558 i=1,j0
         read(10,*) thic0(i),vp0(i),vsh0(i),rho0(i),qp0(i),qs0(i)

         depth0(i) = thic0(i)
         if(i.gt.1) depth0(i) = depth0(i) + depth0(i-1)

         if(vsh0(i).ge.vsmoho) then
	    jmoho = i
            thic0(i) = 0.0
            depth0(i) = depth0(i-1)
	    go to 559
	 endif

558   continue
559   continue
      close(10)

      j0 = jmoho
C END 20140313

C RWG 20160803
C   make sure bottom layer has zero thickness so reflected rays are computed properly

      thic0(j0) = 0.0

      nlskip = -99
      vpsig = -1.0
      vshsig = -1.0
      rhosig = -1.0
      qssig = -1.0
      read(5,*) nlskip,vpsig,vshsig,rhosig,qssig,icflag
      read(5,'(a256)') velname                          

      if(icflag.ne.0.and.icflag.ne.1) icflag=1
      if(vpsig.le.0.0.and.vshsig.le.0.0.and.rhosig.le.0.0.and.qssig.le.0.0) nlskip=-99

      print*,'nlskip,vpsig,vshsig,rhosig,qssig,icflag= ',nlskip,vpsig,vshsig,rhosig,qssig,icflag
      print*,'velname= ',velfile
      print*,'j0= ',j0

      rvsig1 = -1.0
      rvfmax = 1.4

      fasig1 = -1.0
      fasig2 = -1.0
      read(5,*) fasig1,fasig2,rvsig1

CCC RWG modified 2014-09-19
      read(5,*) ipdur_model

c set total time duration to be source duration (Ds) + path duration (Dp)
c see Boore and Thompson (BSSA 2014, 2015).

c     path duration options:
c
c     ipdur_model = 0, r0 = 0.0, dpdr = 0.063:	GP2010 formulation
c     ipdur_model = 1, r0 = 0.0, dpdr = 0.070:	WUS modification (trial/error)
c     ipdur_model = 2, r0 = 0.0, dpdr = 0.100:	ENA modification (trial/error)
c     ipdur_model = 11, WUS formulation of BT2014, overpredicts for multiple rays
c     ipdur_model = 12, ENA formulation of BT2015, overpredicts for multiple rays

      if(ipdur_model.le.0) then

         ndur = 1

         rdur(1) = 0.0
         dpth(1) = 0.0
         dpdr(1) = 0.063

      else if(ipdur_model.eq.1) then

         ndur = 1

         rdur(1) = 0.0
         dpth(1) = 0.0
         dpdr(1) = 0.07

      else if(ipdur_model.eq.2) then

         ndur = 1

         rdur(1) = 0.0
         dpth(1) = 0.0
         dpdr(1) = 0.1

      else if(ipdur_model.eq.11) then

         ndur = 6

         rdur(1) = 0.0
         rdur(2) = 7.0
         rdur(3) = 45.0
         rdur(4) = 125.0
         rdur(5) = 175.0
         rdur(6) = 270.0

         dpth(1) = 0.0
         dpth(2) = 2.4
         dpth(3) = 8.4
         dpth(4) = 10.9
         dpth(5) = 17.4
         dpth(6) = 34.2

      else if(ipdur_model.eq.12) then

         ndur = 8

         rdur(1) = 0.0
         rdur(2) = 15.0
         rdur(3) = 35.0
         rdur(4) = 50.0
         rdur(5) = 125.0
         rdur(6) = 200.0
         rdur(7) = 392.0
         rdur(8) = 600.0

         dpth(1) = 0.0
         dpth(2) = 2.6
         dpth(3) = 17.5
         dpth(4) = 25.1
         dpth(5) = 25.1
         dpth(6) = 28.5
         dpth(7) = 46.0
         dpth(8) = 69.1

      endif

      if(ndur.ne.1) then

         do 4389 i=1,ndur-1
            dpdr(i) = (dpth(i+1)-dpth(i))/(rdur(i+1)-rdur(i))
4389     continue
         dpdr(ndur) = dpdr(ndur-1)

      endif

            print*,ipdur_model,ndur
      do kk=1,ndur
            print*,'rdur= ',rdur(kk),'dpth= ',dpth(kk),'dpdr= ',dpdr(kk)
      enddo

CCC END 2014-09-19

CCC RWG modified 2018-08-30
c adjust stress parameter based on magnitude and area relative to
c standard scaling relation
c     stress parameter adjustment options:
c
c     ispar_adjust = 0, no adjustment
c     ispar_adjust = 1, adjust to Leonard (2010) for active tectonic
c     ispar_adjust = 2, adjust to Leonard (2010) for stable continent

      read(5,*) ispar_adjust,targ_mag,fault_area,spar_exp

c Check to see if top layer of vel model has very thin "air" layer
c this is necessary to get the correct FS reflection coefficient if
c using prduct() and sgc with surface reflected rays

      if(depth0(1).gt.0.001.and.vp0(1).gt.0.01) then
         j0 = j0 + 1
         nlskip = nlskip + 1
         do 4433 i=j0,2,-1
            depth0(i) = depth0(i-1)
            thic0(i) = thic0(i-1)
            vp0(i) = vp0(i-1)
            vsh0(i) = vsh0(i-1)
            rho0(i) = rho0(i-1)
            qp0(i) = qp0(i-1)
            qs0(i) = qs0(i-1)
4433     continue

         depth0(1) = 0.0001
         thic0(1) = 0.0001
         vp0(1) = 0.001
         vsh0(1) = 0.0005
         rho0(1) = 0.001

      endif

       print*,'j0= ',j0

      tdur=1.0

       if(fhil.gt.1./2./dt) then
       print*,'fhigh must be < ',  1./2./dt ,fhil,dt
C       pause 
       end if

c avg Subfault dimensions (dlm) and max slip
      dlm = 0.0
      amx2=0.
      do iv=1,nevnt

C         dlm = 0.5*(dx(iv)+dw(iv)) + dlm
         dlm = sqrt(dx(iv)*dw(iv)) + dlm
Cccc         dlm = sqrt(dx(iv)*dw(iv)/pai) + dlm

         do j=1,nw(iv)
            do i=1,nx(iv)
               amx2=amax1(amx2,abs(sddp(iv,i,j)))
            enddo
         enddo
      enddo

      slip_max=amx2
      dlm = dlm/(float(nevnt))

CFFF convert relative slip into relative moment
CFFF also determine total moment if absolute slips are input (sm=-1)

c compute average corner frequency and average risetime here

      islip_weight_avg = 1
      islip_weight_avg = 0

      nstot = 0
      xsum = 0.0
      fce_avg = 0.0
      trise_avg = 0.0
      do 5000 iv=1,nevnt
         dwdj = dw(iv)*sin(dipq(iv)*pu)
         do 4799 j=1,nw(iv)
	 
            zdep = dtop(iv) + (j-0.5)*dwdj
            do 543 k=1,j0
               if(zdep.le.depth0(k)) then
                  bet = vsh0(k)
                  goto 542
               endif
543         continue
542         continue
         
	    xmu = vsh0(k)*vsh0(k)*rho0(k)*dx(iv)*dw(iv)
	    print *,'j= ',j,'  mu= ',xmu

            do 4798 i=1,nx(iv)

               sddp(iv,i,j) = xmu*sddp(iv,i,j)
	       xsum = xsum + sddp(iv,i,j)

c adjust for fault depth-> scale rupture velocity factor
               dmin = shal_dmin_default
               dmax = shal_dmax_default

               rvf = rvfac*shal_rvfac
               if(zdep.ge.dmin.and.zdep.lt.dmax) then
                  rvf = rvfac*(shal_rvfac + (1.0 - shal_rvfac)*(zdep-dmin)/(dmax-dmin))
               else if(zdep.ge.dmax) then
                  rvf = rvfac
               endif

               dmin = deep_dmin_default
               dmax = deep_dmax_default
               if(zdep.ge.dmin.and.zdep.lt.dmax) then
                  rvf = rvfac*(1.0 + (deep_rvfac - 1.0)*(zdep-dmin)/(dmax-dmin))
               else if(zdep.ge.dmax) then
                  rvf = rvfac*deep_rvfac
               endif
c End rvfac depth adjustment

c alphaT adjustment
               avgdip = dipq(iv);

               fD = 0.0
               if(avgdip.le.90.0.and.avgdip.gt.45.0) then
                  fD = 1.0 - (avgdip - 45.0)/45.0
               else if(avgdip.le.45.0.and.avgdip.ge.0.0) then
                  fD = 1.0;
               endif

               avgrak = rakeq(iv)

3256           continue
               if(avgrak.lt.-180.0) then
                  avgrak = avgrak + 360.0
                  go to 3256
               endif

3257           continue
               if(avgrak.gt.180.0) then
                  avgrak = avgrak - 360.0
                  go to 3257
               endif

               fR = 0.0
               if(avgrak.le.180.0.and.avgrak.ge.0.0) then
                  fR = 1.0 - sqrt((avgrak - 90.0)*(avgrak - 90.0))/90.0
               endif

               alphaT = 1.0/(1.0 + fD*fR*Calpha)
c End alphaT

               zz = Czero*(1.0 + fcfac)/alphaT
               fce = zz*rvf*bet/(dlm*pai)

               trise = rist(iv,i,j)

               if(islip_weight_avg.eq.1) then
                  fce = sddp(iv,i,j)*fce
                  trise = sddp(iv,i,j)*trise
               endif

               fce_avg = fce_avg + fce
               trise_avg = trise_avg + trise

               nstot = nstot + 1

c       print*,'rvf1= ',rvf
c      print *,'zz1= ',zz,rvf,bet,dlm,pai
c      print *,'alphaT1= ',alphaT,fD,fR,Calpha
c      print*,'fce1= ',fce

4798        continue
4799     continue
5000  continue

C End corner frequency and risetime

      if(sm.lt.0.0) then
	 sm = 1.0e+20*xsum
      endif

      if(islip_weight_avg.eq.1) then
         xnorm = 1/xsum
      else
         xnorm = 1.0/(float(nstot))
      endif

      fce_avg = fce_avg*xnorm
      trise_avg = trise_avg*xnorm

C trise_avg same as genslip-vX.X.X, based on somerville
ccc      trise_avg = 1.6e-09*exp(alog(sm)/3.0)

C specify fcmain based on trise_avg, fcmain = Fcoef/trise_avg
C
C scaling factor Fcoef is subjective, could be:
C   Fcoef = 1.0/(2*pi), intuitive
C   Fcoef = 5.0/(2*pi), t=0.2*trise is when Mliu reaches 62% of final
C           slip, 62% of final slip corresponds to Czero=2.1 for brune
C   Fcoef = Czero/(2*pi), beresnev & atkinson, 1997

ccc      Fcoef = 5.0/(2*pai)
ccc      Fcoef = 1.0/(2*pai)

      Fcoef = Czero/(2*pai)

      fcmain = Fcoef/trise_avg

      print*,'Maximum slip ',slip_max
      print*,'Average fce ',fce_avg
      print*,'Target fcmain ',fcmain

CFFF end

c normalize relative moments to have average weight of unity

      amx2=0.
      nstot = 0
      do iv=1,nevnt

         do j=1,nw(iv)
         do i=1,nx(iv)

CCC            amx2 = amx2 + sddp(iv,i,j)/(float(nx(iv)*nw(iv)))

	    if(sddp(iv,i,j).gt.0.001) then
               amx2 = amx2 + sddp(iv,i,j)
               nstot = nstot + 1
	    endif

         enddo
         enddo

      enddo

      amx2 = (float(nstot))/amx2
      do iv=1,nevnt

         do j=1,nw(iv)
         do i=1,nx(iv)
            sddp(iv,i,j)=sddp(iv,i,j)*amx2
         enddo
         enddo

      enddo

CCC RWG modified 2018-08-30

c adjust stress parameter based on magnitude and area relative to
c standard scaling relation

c     stress parameter adjustment options:
c
c     ispar_adjust = 0, no adjustment
c     ispar_adjust = 1, adjust to average of Leonard (2010) and Hanks-Bakun (2014) for active tectonic
c     ispar_adjust = 2, adjust to Leonard (2010) for stable continent

      if(targ_mag.le.0.0) then
         targ_mag = 2.0*(alog(sm)/alog(10.0))/3.0 - 10.73
      endif

      if(fault_area.le.0.0) then
         fault_area = farea_in
      endif

      if(ispar_adjust.eq.1) then

         targ_area = exp((targ_mag-3.99)*alog(10.0))

         if(targ_mag.ge.6.71) then
            targ_area = 0.5*(targ_area + exp(3.0*(targ_mag-3.07)*alog(10.0)/4.0))
         endif

         spar_fac = targ_area/fault_area

      else if(ispar_adjust.eq.2) then

         spar_fac = exp((targ_mag-4.19)*alog(10.0))/fault_area

      else

         spar_fac = 1.0

      endif

      spar_fac = exp((spar_exp)*alog(spar_fac))

      stress_average = spar_fac*stress_average

CCC END 2018-08-30

c Seismic moment for the small event SMOE
c is calculated based on the average stress on the fault

      smoe = stress_average*dlm*dlm*dlm*1.e+21
      ratio = sm/(smoe*float(nstot))
      nsum = int(ratio+0.5)

C RWG 04/20/04
C force nsum=1 and use Frankel convolution operator in stoc_f()
      nsum = 1
      bigC1 = sm/(smoe*float(nstot*nsum))

C  RWG 2018-08-29
C scaling by sqrt(ntot) removes dependence on subfault dimension
C by canceling sqrt(N) scaling due to incoherent summation
      bigC1b = sm/(smoe*sqrt(1.0*nstot))

C  RWG 2018-05-27, Method 1:
C
C Changed computation of bigC to remove dependence of result on subfault
C discretization.  Not sure why, but previous implementation had
C dependence of ~dlm scaling.  Method outlined below removes this.
C
C Following Frankel (1995) and Boore (2009), bigC is the square of the ratio
C of the corner frequencies of the small event (fce) to the main event (fcm):
C
C     bigC = (fce/fcm)^2  following Frankel (1995) and Boore (2009)
C
C Approximate these by using formula for 'fc' from Boore (1983):
C
C   fc = 4.9e+06 * Vs * (sdrop/Mo)^1/3
C
C for subevent:   fce = 4.9e+0.6 * Vs * (sdrop/smoe)^1/3
C for main event: fcm = 4.9e+0.6 * Vs * (sdrop/sm)^1/3
C
C bigC = [ 4.9e+06 * Vs * (sdrop/smoe)^1/3 ]^2 / [ 4.9e+06 * Vs * (sdrop/sm)^1/3 ]^2
C      = (sm/smoe)^2/3
C
      bigC2 = exp((2.0/3.0)*alog(sm/smoe))

C  RWG 2018-05-29, Method 2:
C
C Direct computation of bigC using fce_avg and fcmain: bigC = (fce_avg/fcmain)^2

      bigC3 = (fce_avg/fcmain)*(fce_avg/fcmain)

      bigC = 0.5*(bigC2+bigC3)
      bigC = bigC3
      bigC = bigC1
      bigC = bigC2
      bigC = bigC1b

      print*,'sp= ',stress_average,' Me= ',smoe,' nstot= ',nstot
      print*,'targ_mag= ',targ_mag,' fault_area= ',fault_area,' spar_fac= ',spar_fac
      print*,'   target Mo= ',sm
      print*,'   summed Mo= ',smoe*nstot*nsum
      print*,'  bigC1= ',bigC1,'  bigC1b= ',bigC1b,'  bigC2= ',bigC2,' bigC3= ',bigC3
      print*,'  bigC= ',bigC

c      ratio = 1.0
      
c  Stations location
           
      open(1,file=asite)

      ndata=int(duration/dt)
      if(ndata.gt.mmv) ndata=mmv

      ifu=irand
      call init_random_seed(irand)
   
      call RANU2(NR,RNA)
      call RANU2(NR,RNB) 

      loc1= index(outname,' ')-1

      head_lines = 0
      do
         read(1,*) dummy
         !This line checks to see if input starts with a '#' or '%'
         !If so, ignore
         if ((dummy(1:1) .ne. '#') .AND. (dummy(1:1) .ne. '%')) exit
         head_lines=head_lines+1
      enddo

      rewind(1)
      do j=1,head_lines
         read(1,*)
      enddo

      do 555 msite=1,nsite
      read(1,*,end=1) stlon,stlat,cap

      loc2=index(cap,' ')-1
      outname1=outname(1:loc1)//'_'//cap(1:loc2)//'.090'
      outname2=outname(1:loc1)//'_'//cap(1:loc2)//'.000'
      outname3=outname(1:loc1)//'_'//cap(1:loc2)//'.ver'

      cname(1) = '090'
      cname(2) = '000'
      cname(3) = 'ver'

      open(22,file=outname1)
      open(23,file=outname2)
      open(24,file=outname3)

      d10=1000.

      do 20 i=1,mmv 
      do 20 ja=1,3      
      ds(ja,i)=0.      
20    continue
                                                           
      print*,' '                                       
      print*,'Site : ', msite,' ifu= ',ifu,outname

      if(nlskip.ge.0) then

         call grandvel(j0,nlskip,vpsig,vshsig,rhosig,qssig,icflag)

	 if(velname(1:2).ne."-1") then

            loc1= index(velname,' ')-1
            velname1=velname(1:loc1)//'_'//cap(1:loc2)//'.1d'

            open(26,file=velname1)

	    kst = 0
	    if(thic(1).le.0.001.and.vp(1).le.0.01) kst = 1

            write(26,'(i6)') (j0-kst)

	    do 9190 k=1+kst,j0

               write(26,9192)thic(k),vp(k),vsh(k),rho(k),qp(k),qs(k)

9190        continue

            close(26)
	 endif

      else

         do 9213 k=1,j0

            depth(k) = depth0(k)
            thic(k) = thic0(k)
            vp(k) = vp0(k)
            vsh(k) = vsh0(k)
            rho(k) = rho0(k)
            qp(k) = qp0(k)
            qs(k) = qs0(k)

9213     continue

      endif

      if(fasig1.gt.0.0.or.fasig2.gt.0.0.or.rvsig1.gt.0.0) then

         call normal_random_number(mmv,fgrand)
         print *,' fasig1= ',fasig1*fgrand(1)

cCCCC       do ir=1,mmv
cCCCC       print*,'fgrand= ',fgrand(ir)
cCCCC       enddo
      endif

c loop over fault segments
      do 13 iv=1,nevnt

         stra=strq(iv)*pu
         dipa=dipq(iv)*pu
         raka=rakeq(iv)*pu 

c Calculate distance, take-off angle azimuth from the center of the
c  subfault to the station

      call even_dist2(elonq(iv),elatq(iv),stlon,stlat,strq(iv),dipq(iv),
     +   dtop(iv),astop(iv),dx(iv),dw(iv),nx(iv),nw(iv),rlsu,phsu,thsu,dst,zet)

c Time Duration Model

      tmax = 0.0
      d10=10000.
      do 893 j=1,nw(iv)
         do 894 i=1,nx(iv)

            do ksrc=1,j0
               if (depth(ksrc).ge.zet(i,j)) then
                  bet=vsh(ksrc)
                  goto 8966
               end if
            enddo
8966        continue

c get fce, tw0 = 1.0/fce + dlm/vr

            zz = Czero

c adjust for fault depth-> scale rupture velocity factor
            dmin = shal_dmin_default
            dmax = shal_dmax_default

	    rvf = rvfac*shal_rvfac
            if(zet(i,j).ge.dmin.and.zet(i,j).lt.dmax) then
               rvf = rvfac*(shal_rvfac + (1.0 - shal_rvfac)*(zet(i,j)-dmin)/(dmax-dmin))
            else if(zet(i,j).ge.dmax) then
               rvf = rvfac
            endif

            dmin = deep_dmin_default
            dmax = deep_dmax_default
            if(zet(i,j).ge.dmin.and.zet(i,j).lt.dmax) then
               rvf = rvfac*(1.0 + (deep_rvfac - 1.0)*(zet(i,j)-dmin)/(dmax-dmin))
            else if(zet(i,j).ge.dmax) then
               rvf = rvfac*deep_rvfac
            endif

c 20131120 RWG added alphaT adjustment
            avgdip = dipq(iv);

            fD = 0.0
            if(avgdip.le.90.0.and.avgdip.gt.45.0) then
               fD = 1.0 - (avgdip - 45.0)/45.0
            else if(avgdip.le.45.0.and.avgdip.ge.0.0) then
               fD = 1.0;
            endif

            avgrak = rakeq(iv)
3526        continue
            if(avgrak.lt.-180.0) then
               avgrak = avgrak + 360.0
               go to 3526
            endif
3527        continue
            if(avgrak.gt.180.0) then
               avgrak = avgrak - 360.0
               go to 3527
            endif

            fR = 0.0
            if(avgrak.le.180.0.and.avgrak.ge.0.0) then
               fR = 1.0 - sqrt((avgrak - 90.0)*(avgrak - 90.0))/90.0
            endif

            alphaT = 1.0/(1.0 + fD*fR*Calpha)
c End alphaT

            zz = Czero*(1.0 + fcfac)/alphaT

CCC RWG modified 2014-09-19
c
c set total time duration to be source duration (Ds) + path duration (Dp)
c see Boore and Thompson (BSSA 2014, 2015).

c     tw0 = 1.0/fce = Ds
c
c     path duration options:
c
c     ipdur_model = 0, r0 = 0.0, dpdr = 0.063:	GP2010 formulation
c     ipdur_model = 1, r0 = 0.0, dpdr = 0.070:	WUS modification (trial/error)
c     ipdur_model = 2, r0 = 0.0, dpdr = 0.100:	ENA modification (trial/error)
c     ipdur_model = 11, WUS formulation of BT2014, overpredicts for multiple rays
c     ipdur_model = 12, ENA formulation of BT2015, overpredicts for multiple rays

            do 3461 kk=1,ndur
               if(rlsu(i,j).gt.rdur(kk)) then
	          r0 = rdur(kk)
		  d0 = dpth(kk)
		  slp = dpdr(kk)
               endif
3461        continue

            fce = zz*rvf*bet/(dlm*pai)
            tw0 = 1.0/fce
            tw0 = sqrt(bigC)*tw0

            dpath = d0 + slp*(rlsu(i,j)-r0)
            twin(i,j) = tw0 + dpath

CCC END 2014-09-19

	    if(twin(i,j).gt.tmax) tmax = twin(i,j)
            d10=amin1(d10,rlsu(i,j))

894      continue
893   continue

      ntmax = int(2.0*tmax/dt)

      np2 = 2
2776  if(np2.ge.ntmax) go to 2777
      np2=np2*2
      go to 2776
2777  continue
CNNNN      np2 = 8192

      if(np2.gt.mm) then
      
         print *,'np2= ',np2,' > array dimension mm= ',mm
         print *,'need to recompile with larger array size, exiting...'
         go to 9555

      endif

      print *,'np2= ',np2

      if(fasig1.gt.0.0.or.fasig2.gt.0.0) then

         i0=1
         ismth=np2/2
         do 8461 i=1,np2

            i1 = (i-1)/ismth
            i2 = i1 + 1
            fgrands(i) = fgrand(i0+i1) + ((i-1)-i1*ismth)*(fgrand(i0+i2)-fgrand(i0+i1))/(float(ismth))

8461     continue
         fgrands(1) = fgrand(1)

      endif

      nfold = np2/2+1 
      mfold = np2/2-1
      df = 1.0/(np2*dt)
      do 777 i=1,nfold
         dfr(i)=df*float(i-1)
 777  continue

      irandcnt = 1
      do 4 i=1,nx(iv)      
      do 4 j=1,nw(iv)

      print *,'sddp= ',sddp(iv,i,j),iv,i,j
 
         if(sddp(iv,i,j).lt.0.001) goto4  

      do il=1,np2
       stdd(il,1)= 0.
       stdd(il,2)= 0. 
       stdd(il,3)= 0.
      end do

       bet=vsh(1)
       ro=rho(1)

       do ksrc=1,j0
       if (depth(ksrc).ge.zet(i,j)) then
       bet=vsh(ksrc)
       ro=rho(ksrc)
       goto 66
       end if
       enddo 
       print*,'wrong!'
     
66    continue

c Begin subfault corner frequency

c 20171219, changed to use Czero formulation
c
      zz = Czero

c adjust for fault depth-> scale rupture velocity factor
      dmin = shal_dmin_default
      dmax = shal_dmax_default

      rvf = rvfac*shal_rvfac
      if(zet(i,j).ge.dmin.and.zet(i,j).lt.dmax) then
         rvf = rvfac*(shal_rvfac + (1.0 - shal_rvfac)*(zet(i,j)-dmin)/(dmax-dmin))
      else if(zet(i,j).ge.dmax) then
         rvf = rvfac
      endif

            dmin = deep_dmin_default
            dmax = deep_dmax_default
            if(zet(i,j).ge.dmin.and.zet(i,j).lt.dmax) then
               rvf = rvfac*(1.0 + (deep_rvfac - 1.0)*(zet(i,j)-dmin)/(dmax-dmin))
            else if(zet(i,j).ge.dmax) then
               rvf = rvfac*deep_rvfac
            endif

      rvf0 = rvf
      if(rvsig1.gt.0.0) then

         irandcnt = irandcnt + 1
         rvf = rvf0*exp(fgrand(irandcnt)*rvsig1)

         if(rvf.gt.rvfmax) rvf = rvfmax

      endif

       print*,'rvf= ',rvf,rvf0

c 20131120 RWG added alphaT adjustment
            avgdip = dipq(iv);

            fD = 0.0
            if(avgdip.le.90.0.and.avgdip.gt.45.0) then
               fD = 1.0 - (avgdip - 45.0)/45.0
            else if(avgdip.le.45.0.and.avgdip.ge.0.0) then
               fD = 1.0;
            endif

            avgrak = rakeq(iv)
2256        continue
            if(avgrak.lt.-180.0) then
               avgrak = avgrak + 360.0
               go to 2256
            endif
2257        continue
            if(avgrak.gt.180.0) then
               avgrak = avgrak - 360.0
               go to 2257
            endif

            fR = 0.0
            if(avgrak.le.180.0.and.avgrak.ge.0.0) then
               fR = 1.0 - sqrt((avgrak - 90.0)*(avgrak - 90.0))/90.0
            endif

            alphaT = 1.0/(1.0 + fD*fR*Calpha)
c End alphaT

      zz = Czero*(1.0 + fcfac)/alphaT

      fce=zz*rvf*bet/dlm/pai

      print *,'zz= ',zz,rvf,bet,dlm,pai
      print *,'alphaT= ',alphaT,fD,fR,Calpha
      print*,'fce= ',fce
c End subfault corner frequency

c  Rise Time  : Rise
c no modification for the rise time
             rise=rist(iv,i,j)

cNNNN start of ray loop

c calculate S wave arrival time and ray path length for this ray
c geometric decay term is 1./(ray path length)
c set sub_tstart = stime - tw_eps*twin => envelope max conincide with S time
c
c # layers in vel model is j0
c horizontal range is dst(i,j)
c source depth is zet(i,j)
c
c mode = 3 for SV
c mode = 4 for SH
c mode = 5 for P
c
c irtype = 1 for direct ray
c irtype = 2 for moho reflection

c hardwire for SH since only geometric term is considered
      mode = 4
      do 324 ir=1,nrtyp

      irtyp = irtype(ir)
      if(irtyp.eq.0) irtyp = 1

      call gf_amp_tt(j0,zet(i,j),dst(i,j),sgc,irtyp,mode,p0,stime,rpath,qbar)
      sub_tstart = stime - tw_eps*twin(i,j)

c RWG 20201015
c  Added conditional on start time to prevent negative value.
c  Only occurs rarely, but Scott C. ran into issues with some CyberShake
c  runs for sites very close to rupture.
      if(sub_tstart.lt.0.5*stime) sub_tstart = 0.5*stime

      if(irtype(ir).eq.0) then
	 rpath = rlsu(i,j)
	 qbar = rpath/(bet*150.0)
	 stime = rpath/3.7
         sub_tstart = 0.7*stime
      endif

      print *,'tt= ',sub_tstart,' rp= ',rpath,' bet= ',bet,' twin= ',twin(i,j)
      print *,'hs= ',zet(i,j),' range= ',dst(i,j),' rho= ',ro

      tvmin=amin1(tvmin,sub_tstart)  

cXXX      print*,'range= ',dst(i,j),' rpath= ',rpath,' t0=',stime,' p0=',p0
      print*,'irtype= ',irtype(ir),' t0= ',stime,' rpath= ',rpath,' az= ',phsu(i,j)/pu

CFFF the following 9 lines are new
      tw = twin(i,j)
      do 9876 kf=1,3

         fmx1=fmx
         if(fmx1.gt.15.and.kf.eq.3) fmx1=15.

         call stoc_f(np2,rpath,tw,tw_eps,tw_eta,bet,ro,dt,smoe,dlm,fce,fmx1,akapp,cs(1,kf),dfr,qbar,qfexp,bigC)

9876  continue

cNNN need to do something about siteamp if we consider Ref/Trans effects
      if(isite_amp.ne.0) then

	 call get_sitefacs(ksrc,nsfac,fn,an)

cc         do 8997 k=1,nsfac
cc8997     print *,k,ksrc,' z= ',zet(i,j),' fn= ',exp(fn(k)),' an= ',exp(an(k))

         do 8727 k=1,3
            call siteamp(np2,cs(1,k),dfr,nsfac,fn,an)
8727     continue

      endif

      if(fasig1.gt.0.0.or.fasig2.gt.0.0) then

         do 6727 k=1,3
            call famprand(np2,cs(1,k),fgrands,fasig1,fasig2)
6727     continue

      endif

CFFF end

c   average rad coeff freq. dependent 

CFFF the following 14 lines are new (moved from before)

c need to get incidence angle from p0 => sin(i)/vs=p0
c th=i for downgoing ray, th=180-i for upgoing ray
c
       if(bet*p0.gt.1.0) then
	  th = 0.5*pai
       else
	  th = asin(bet*p0)
       endif

      print *,'  th(deg)= ',th/pu,'  th(rad)= ',th,'  bet*p0= ',bet*p0
      if(bet*p0.lt.1.0) print *,'                asin(bet*p0)= ',asin(bet*p0)

CCC       if(irtype(ir).eq.1) th = pai - th
       if(irtype(ir).eq.0) then
          th = thsu(i,j)
       else if(mod(irtype(ir),2).eq.1) then
          th = pai - th
       endif

       pa= phsu(i,j) 

cNNN  cmp = pa gives radial component
c     cmp = pa + 90*pu gives tangential component

      cmp=-90.*pu
      CALL RADFRQ_lin(STRA,DIPA,RAKA,pa,th,DFR,NFOLD,CMP,flol,RNA,RNB,nr,rdna) 
      call highcor_f(nfold,mfold,np2,cs(1,1),stdd(1,1),rdna)

      print *,'   rdna-090= ',rdna(30),' th= ',th/pu,' ir= ',ir

      cmp=0
      CALL RADFRQ_lin(STRA,DIPA,RAKA,pa,th,DFR,NFOLD,CMP,flol,RNA,RNB,nr,rdna)
      call highcor_f(nfold,mfold,np2,cs(1,2),stdd(1,2),rdna)

      print *,'   rdna-000= ',rdna(30),' th= ',th/pu,' ir= ',ir

       call RADV_lin(STRA,DIPA,RAKA,pa,th,DFR,NFOLD,flol,RNA,RNB,nr,rdna) 
       call highcor_f(nfold,mfold,np2,cs(1,3),stdd(1,3),rdna)
CFFF end

CFFF the following 9 lines are new (moved from before)
       if(vr.le.0.0) then
          ratim = rupt(iv,i,j)
       else
          xra = shyp(iv) - (i - 0.5*(nx(iv)+1))*dx(iv)
          yra = dhyp(iv) - (j - 0.5)*dw(iv)
          ratim=SQRT(xra**2+yra**2)/VR

          if(irand.gt.0) then
             ratim=ratim+(rand_numb(0)-0.5)*0.1*ratim
          endif
       endif

       kst=int(ratim/dt) + int(sub_tstart/dt)
CFFF end

       do 1722 k=1,nsum
            si=rand_numb(0)
c            dris=(float(k)-1.0+si)*rise/dt
            dris=si*rise/dt

            k2=nint(dris)
            if(nsum.eq.1) k2=0

CFFF the following 10 lines are modified

       k2 = k2 + kst
       kend = k2 + np2
       if(kend.gt.ndata) kend = ndata

       print*,'k2= ',k2,sub_tstart,ratim

C  RWG 04/20/04
C  remove scaling by 'ratio', now done in stoc_f() using bigC

       DO 1733 li=k2,kend
            DS(1,li) =DS(1,li)+sddp(iv,i,j)*stdd(li-k2,1)
            DS(2,li) =DS(2,li)+sddp(iv,i,j)*stdd(li-k2,2)
            DS(3,li) =DS(3,li)+sddp(iv,i,j)*stdd(li-k2,3)
1733     CONTINUE 

1722     CONTINUE 

CFFF end

324    continue

cNNNN end of ray loop

4      continue
13     continue

      if(iftt.gt.0) then 
         call filter3d(nbu,iftt,flol,fhil,ds,ndata,dt)
      endif

      amx1=0.0
      amx2=0.0
      amx3=0.0
      DO 134 I=1,ndata
      amx1=amax1(amx1,abs(ds(1,i)))
      amx2=amax1(amx2,abs(ds(2,i)))
      amx3=amax1(amx3,abs(ds(3,i)))  
  134 CONTINUE                                                                  
      WRITE(6,*) 'ACC.MAX=',amx1,amx2,amx3
      print*,'Closest Distance ', d10,'  fc ',fce

	ic(1) = 0
	ic(2) = 0

	c(1) = delay
	c(2) = d10
	c(3) = 0.0
	c(4) = 0.0

        do l=1,3
        io=21+l
        write(io,100)cap,cname(l),outname
        write(io,102)ndata,dt,(ic(k),k=1,2),(c(k),k=1,4)

        write(io,'(6e13.5)')(ds(l,i),i=1,ndata)
        close(io)
        end do      
100     format(a8,2x,a4,1x,a56)
102     format(1x,i6,1x,e11.5,2i3,4(1x,f10.4))   

9192    format(f8.4,f8.4,f8.4,f8.4,f8.2,f8.2)
555   continue

9555  continue
                                                                 
1       END                                           
      
        subroutine distazi(evlar,evlor,stlar,stlor,delta,azi,bazi)
  
        pai=3.141592654
  
        evla=evlar/180.*pai
        evlo=evlor/180.*pai
        stla=stlar/180.*pai
        stlo=stlor/180.*pai
 
        
        ae=cos(evla)*cos(evlo)
        as=cos(stla)*cos(stlo)
        be=cos(evla)*sin(evlo)
        bs=cos(stla)*sin(stlo)
        ce=sin(evla)
        cs=sin(stla)

        a=sqrt((ae-as)**2+(be-bs)**2+(ce-cs)**2)/2.0
        at=atan(a/sqrt(-a*a+1.0))*2.0

        delta=at*6371.0
        ae=sin(stla)*cos(evla)-cos(stla)*sin(evla)*cos(stlo-evlo)
        be=sin(evla)*cos(stla)-cos(evla)*sin(stla)*cos(evlo-stlo)

        ae=ae/sin(at)

c        if(abs(ae).ge.1.) ae=ae*0.99999999/abs(ae)
        if(abs(ae).ge.0.9999) ae=ae*1.0/abs(ae)


        be=be/sin(at)
        aa=sqrt(1-ae**2)/ae

        aaa=sqrt(1-be**2)/be
        azi=atan(aa)/pai*180.0
        bazi=atan(aaa)/pai*180.0
        if ((evlo.ge.stlo).and.(evla.gt.stla)) then
        azi=180.+abs(azi)
        else
	if ((evlo.gt.stlo).and.(evla.le.stla)) then
        azi=360.-azi
        else
	if ((evlo.lt.stlo).and.(evla.ge.stla)) then
        azi=180.-abs(azi)
        endif
        endif
        endif
        return
        end 

      subroutine stoc(nfd,r,tdur,betvs,row,dt,smt,dlm,fc,fmx,akapp,nsum,
     +                cw)                         
      parameter (mmv=30000)
      DIMENSION A(16348),W(16348),AS(16348),DFR(16348),GA(16348)
      DIMENSION CW(mmv),CWS(16348),SPF(16348)

      real*8 gsa,gm,dgamm
      COMPLEX AC(16348),ACM(16348) 
      PAI=3.1415926                          
      rp=0.63 
      fs=2.0  
      prtitn=0.71 
      eps=0.2 
      tw=0.                                                                                                                                                                                 
c      fmx=15.0  
      seism=smt*nsum                                                          

      nf=nn/2+1
      nn=2
 776  if(nn.ge.nfd) go to 777
      nn=nn*2
      go to 776
 777  nf=nn/2+1 

      DO 2 I=1,NN 
      AS(I)=0.0 
    2 CONTINUE                                                                  
      ACCMA=0.0 
      AVLMA=0.0                                                                 
      IX=0                                                                      
      Rxx=r*100000.    

c Rise Time model
      AMD=(ALOG10(seism)-17.0)/1.33 
      tw=10.0**(0.31*AMD-0.774)
cXXX
cXXX      print*,'tw: ',tw
cXXX      tw=2.0

      if(r.gt.10.) tw=tw+0.063*(r-10.)
      tw=tw*tdur
      NT=TW/DT                            
      DF=1.0/FLOAT(NN)/DT  
      DO 441 I=1,NF  
      DFR(I)=FLOAT(I-1)*DF 
  441 CONTINUE                                                                  
      NP=NN/2  
c      B=-0.2*ALOG(0.05)/(1.0+0.2*(ALOG(0.2)-1.0))  
c      C=B/0.2/TW  
      B=-eps*ALOG(0.05)/(1.0+eps*(ALOG(eps)-1.0))
      C=B/eps/TW
      gsa=2*b+1.0
      GM=DGAMM(gsa)                          
      AA=SQRT((2.0*C)**(2.0*B+1.0)/GM)
c AA=(2.7182/0.2/TW)**BB  
      DO 10 I=1,NT
      T=FLOAT(I-1)*DT    
      W(I)=AA*T**B*EXP(-C*T)  
   10 CONTINUE                                                                
      BETA=BETvs*100000.   
      CC=RP*FS*PRTITN/(4.0*PAI*ROW*BETA**3) 
      DO 20 I=2,NF     
      FR=dfr(i)  

c Q model
                                                        
c      QV=100+10.0 *FR**1.70     
c      QV=270.0*FR**0.5
       QV=150.0*FR**0.5
 
c Beresnev Northridge     QV=150.0*FR**0.5

      OMG=2.0*PAI*FR   
      OMGC=2.0*PAI*FC 
      OMGM=2.0*PAI*FMX    
      A1=CC*SMT*(OMG**2/(1.0+(OMG/OMGC)**2))  
C     A2=(1.0+(OMG/OMGM)**8)**(-0.5) 

      if (akapp.le.0.) then
      A2=(1.0+(OMG/OMGM)**1)**(-1.0)
      else
      a2=exp(-pai*fr*akapp)
      end if
      A3=EXP(-OMG*Rxx/(2.0*QV*BETA))/Rxx                                       
      AS(I)=A1*A2*A3                                                            
   20 CONTINUE                                                                     
      AM=0.0   
      SD=1.0                                                                    
      CALL normal_random_number(NT,A)                                            
      CALL FLZERO(NT,DT,A) 
c      call aver_remov(nt,A)                                                  
                                                                   
CVVVVV      DO 44 I=1,NN  
CVVVVV      AC(I)=(0.0,0.0)
CVVVVV   44 CONTINUE                                                                  
CVVVVV      DO 45 I=1,NT  
CVVVVV      AC(I)=CMPLX(A(I),0.0) 
CVVVVV   45 CONTINUE                                                                  
CVVVVV      CALL FAST(NN,AC,-1)   
CVVVVV      DO 445 I=1,NF    
CVVVVV      CWS(I)=CABS(AC(I))*DT    
CVVVVV  445 CONTINUE                                                                
CVVVVV      FSA=0.0   
CVVVVV      DO 446 I=1,NF  
CVVVVV      FSA=FSA+CWS(I)  
CVVVVV  446 CONTINUE                                                                  
CVVVVV      FSA=FSA/FLOAT(NF)    

      DO 30 I=1,NT   
      GA(I)=A(I)*W(I)
   30 CONTINUE  
                                                                
      DO 75 I=1,NN  
      AC(I)=CMPLX(0.0,0.0) 
   75 CONTINUE                                                                  
      DO 40 I=1,NT 
      LL=I+kstt  
      AC(LL)=CMPLX(GA(I),0.0) 
   40 CONTINUE  
      CALL FAST(NN,AC,-1)     
      DO 71 I=1,NF   
      SPF(I)=CABS(AC(I))*DT    
   71 CONTINUE                                                                  

      FSA=0.0    
      DO 448 I=1,NF      
      FSA=FSA+SPF(I)       
  448 CONTINUE                                                                  
      FSA=FSA/FLOAT(NF)  
      AMP=1.0/FSA   
      DO 90 I=1,NN                              
      AC(I)=AC(I)/FLOAT(NN)            
   90 CONTINUE                                                                  
      DO 100 I=1,NP      
      ACM(I)=AC(I)*AS(I)*AMP     
      ACM(NN-I+1)=CONJG(AC(I+1)*AS(I+1))*AMP   
  100 CONTINUE                                                                  
      ACM(NF)=AC(NF)*AS(NF)*AMP        
      CALL FAST(NN,ACM,1)                                                       
      DO 110 I=1,NN                                                             
      CW(I)=REAL(ACM(I)) 
110   CONTINUE

      n0=nn/10
      do i=1,n0
      t=float(i-1)/n0*1.5707965
      cw(nn-n0+i)=cw(nn-n0+i)*cos(t)
      end do
      do i=nn+1,nfd
      cw(i)=cw(nn)
      enddo                                                                                             
c      CALL INAC(DT,nn,CW,CWI) 
      return                                                                    
      END

      subroutine stoc_f(np2,r,tw,eps,eta,betvs,row,dt,smt,dlm,fc,fmx,akapp,cw,dfr,qb,qfe,bigC)
      DIMENSION A(np2),W(np2),DFR(1)
      DIMENSION SPF(np2)
      real*8 gsa,a1,a2,a3,as(np2),dgamm,gm
      COMPLEX*8 AC(np2),cw(1) 

      PAI=3.1415926                          
      rp=0.63 

      fc2 = fc*fc

cNNN fs is the free surface factor
c    prtitn is vector partition factor for 2 orthogonal comps [1/sqrt(2)]

      fs=2.0
      prtitn=0.71 

      nf = np2/2+1 
      Rxx=r*100000.    

c      B=-0.2*ALOG(eta)/(1.0+0.2*(ALOG(0.2)-1.0))  
c      C=B/0.2/TW  
      B=-eps*ALOG(eta)/(1.0+eps*(ALOG(eps)-1.0))
      C=B/eps/TW
      gsa=2*b+1.0
      GM=DGAMM(gsa)                          
      AA=SQRT((2.0*C)**(2.0*B+1.0)/GM)
c AA=(2.7182/0.2/TW)**BB  

      DO 10 I=1,np2
      T=FLOAT(I-1)*DT    
      W(I)=AA*T**B*EXP(-C*T)  
   10 CONTINUE                                                                

      BETA=BETvs*100000.   
      CC=RP*FS*PRTITN/(4.0*PAI*ROW*BETA**3) 
      OMGC=2.0*PAI*FC 
      OMGM=2.0*PAI*FMX    
      AS(1)=0.0 
      DO 20 I=2,NF     
      FR=dfr(i)  
      fr2 = fr*fr

c Q model
                                                        
c      QV=100+10.0 *FR**1.70     
c      QV=270.0*FR**0.5
       QV=150.0*FR**0.5

cRWGTEST
CC       qv = 1.0e+15
 
c Beresnev Northridge     QV=150.0*FR**0.5

      OMG=2.0*PAI*FR   
      A1=CC*SMT*(OMG**2/(1.0+(OMG/OMGC)**2))  
C     A2=(1.0+(OMG/OMGM)**8)**(-0.5) 

      if (akapp.le.0.) then
      A2=(1.0+(OMG/OMGM)**1)**(-1.0)
      else
      a2=exp(-pai*fr*akapp)
      end if

cNNN division by Rxx is 1/R geometric spreading factor

      A3=EXP(-OMG*Rxx/(2.0*QV*BETA))/Rxx                                       

c   use input q values via qbar
      A3=EXP(-0.5*OMG*qb*(fr**(-0.5)))/Rxx
      A3=EXP(-0.5*OMG*qb*(fr**(-qfe)))/Rxx

C   Frankel convolution operator
      frank = 1.0
      frank = bigC*(fc2 + fr2)/(fc2 + bigC*fr2)

      AS(I)=A1*A2*A3*frank
   20 CONTINUE                                                                     
      CALL normal_random_number(np2,A)                                            

Cxxxx
Cxxxx        write(6,'(e13.5)')(a(i),i=1,np2)
Cxxxx

      CALL FLZERO(np2,DT,A) 

      DO 40 I=1,np2 
      AC(I)=CMPLX(a(i)*w(i),0.0) 
   40 CONTINUE  
                                                                
      CALL FAST(np2,AC,-1)     

C 03/18/2009 RWG changed normalization of random time sequence such that the average
C of the POWER spectrum is unity (not amplitude spectrum), this is consistent
C with Boore (1983)
C
C RWG      FSA=0.0    
C RWG      DO 71 I=1,nf
C RWG      SPF(I)=CABS(AC(I))*DT    
C RWG      FSA=FSA+SPF(I)       
C RWG   71 CONTINUE
C RWG
C RWG      FSA=FSA/FLOAT(nf)  
C RWG      AMP=1.0/FSA

      fsa = 0.0    
      do 71 i=1,nf
         fsa = fsa + cabs(ac(i))*cabs(ac(i))
   71 continue
      amp = 1.0/(dt*sqrt(fsa/float(nf)))

C RWG

      NP=np2/2
      DO 100 I=1,NP      
      cw(I)=AC(I)*AS(I)*AMP     
      cw(np2-I+1)=CONJG(AC(I+1)*AS(I+1))*AMP   
  100 CONTINUE                                                                  
      cw(1)=CMPLX(0.0,0.0) 
      cw(NF)=AC(NF)*AS(NF)*AMP        

      return                                                                    
      END       
                                 
      SUBROUTINE RADFRQ(RDNA,STRA,DIPA,RAKA,PA,THAA,DFR,NFOLD,                  
     & CMP,RNA,RNB,NR,fr1)   
     
c   Frequency Dependent Aver Rad Coeff.     
                                                          
      DIMENSION RNA(nr),RNB(nr),RDNA(nfold),DFR(nfold)

      pu=3.1415926/180.
c      fr1=1.0
      fr2=3.0
      range=20

C   Theoretical radiation coefficient
                        
                                                                   
      CALL RDATN(STRA,DIPA,RAKA,PA,THAA,RDSHA,RDSVA)                            
      RADVL=RDSVA*COS(THAA)*COS(CMP-PA)+RDSHA*SIN(CMP-PA)
      SIGN=RADVL/ABS(RADVL)
      DO 50 I=1,NFOLD                                                           
      RDNA(I)=RADVL                                                              
  50  CONTINUE   
      
c      print*,'Theory: ',radvl                                                                                  
c Average radiation coefficient
      
      tha1=thaa-range*pu
      tha2=thaa+range*pu
     
      if(tha1.lt.90.*pu) tha1=90.*pu
      if(tha2.gt.180.*pu) tha2=180.*pu 
                         
      RADS=0.0                                                                  
      RADV=0.0                                                                  
                                               
      DO 150 K=1,NR                                                             
      TH=ACOS((1.0-RNA(K))*COS(THA1)+RNA(K)*COS(THA2))  
c      fa=pa+360*pu*rnb(k)
      fa=pa+4*range*pu*(0.5-rnb(k))                       
                                              
      CALL RDATN(STRA,DIPA,RAKA,FA,TH,RDSHA,RDSVA)                              
      RADS=RDSVA*COS(TH)*COS(CMP-FA)+RDSHA*SIN(CMP-FA)                          
      RADV=RADV+ABS(RADS)  
c      print*,k,nr,thaa/pu,tha1/pu,tha2/pu,pa/pu,fa/pu,rads                                                   
  150 CONTINUE 
                   
      RADVH=RADV/FLOAT(NR)*SIGN 
c      print*,radvh,nr
c      stop

c       print*,'Theory, Aver: ',radvl,RADVH      
                                         
      DO 200 I=2,NFOLD                                                          
      IF(DFR(I).LE.fr1) GO TO 210                                               
      IF(DFR(I).GT.fr2) GO TO 220                                               
      RADS=0.0                                                                  
      RADV=0.0   
                                                               
      tha1=thaa-range*pu*(dfr(i)-fr1)/(fr2-fr1)
      tha2=thaa+range*pu*(dfr(i)-fr1)/(fr2-fr1)

      if(tha1.lt.90.*pu) tha1=90.*pu
      if(tha2.gt.180.*pu) tha2=180.*pu 
      
                                        
      DO 300 K=1,NR                                                             
      TH=ACOS((1.0-RNA(K))*COS(THA1)+RNA(K)*COS(THA2)) 
c      fa=pa+360.*pu*(dfr(i)-fr1)/(fr2-fr1)*rnb(k)  
      fa=pa+4*range*pu*(dfr(i)-fr1)/(fr2-fr1)*(0.5-rnb(k))                             
                                                    
      CALL RDATN(STRA,DIPA,RAKA,FA,TH,RDSHA,RDSVA)                              
      RADS=RDSVA*COS(TH)*COS(CMP-FA)+RDSHA*SIN(CMP-FA)                          
      RADV=RADV+ABS(RADS)                                                       
  300 CONTINUE                                                                  
      RDNA(I)=RADV/FLOAT(NR)*SIGN   
      GO TO 200                                                                 
  210 RDNA(I)=RADVL
      GO TO 200                                                                 
  220 RDNA(I)=RADVH                                                         
  200 CONTINUE
                                                                      
      RETURN   
                                                                       
      END 

      SUBROUTINE RADUNI(RDNA,STRA,DIPA,RAKA,PA,THAA,DFR,NFOLD,                  
     & CMP,RNA,RNB,NR)
     
c  Aver Rad. Coeff.
                                                          
      DIMENSION RNA(nr),RNB(nr),RDNA(nfold),DFR(nfold) 
      PU=3.1415926/180
      range=20.   
                                                             
c      CALL RDATN(STRA,DIPA,RAKA,PA,THAA,RDSHA,RDSVA)                            c      RDX=RDSVA*COS(THAA)*COS(CMP-PA)+RDSHA*SIN(CMP-PA) 
c      Rdna(1)=RDX
                                           
       tha1=thaa-range*pu
       tha2=thaa+range*pu
      if(tha1.lt.90.*pu) tha1=90.*pu
      if(tha2.gt.180.*pu) tha2=180.*pu 

                               
c      RADS=0.0 
      RADV=0.0                                                                  
                                               
      DO 150 K=1,NR                                                             
      TH=ACOS((1.0-RNA(K))*COS(THA1)+RNA(K)*COS(THA2))                          
c      fa=pa+120*pu*rnb(k) 
      fa=pa+4*range*pu*(0.5-rnb(k))
                                             
      CALL RDATN(STRA,DIPA,RAKA,FA,TH,RDSHA,RDSVA)                              
      RADS=RDSVA*COS(TH)*COS(CMP-FA)+RDSHA*SIN(CMP-FA)                          
      RADV=RADV+ABS(RADS)                                                       
  150 CONTINUE                             
                                                  
      RADVH=RADV/FLOAT(NR)                                                      
                                               
      DO 100 I=1,NFOLD                                                          
      RDNA(I)=RADVH                                                             
  100 CONTINUE                                                                  
      RETURN                                                                    
      END     
               
                                                             
      SUBROUTINE RADTRL(STRA,DIPA,RAKA,PA,THAA,DFR,NFOLD,                 
     & CMP,rdna) 
     
c Theoretical rad. coeff
      parameter (mm=16348)                               
                  
      DIMENSION Rdna(mm),DFR(mm)
      PU=3.1415926/180   
                                                             
      CALL RDATN(STRA,DIPA,RAKA,PA,THAA,RDSHA,RDSVA)                            
      RDX=RDSVA*COS(THAA)*COS(CMP-PA)+RDSHA*SIN(CMP-PA) 
                                                                              
      DO 101 I=1,NFOLD                                                       
      Rdna(I)=RDX                                                              
  101 CONTINUE                                                                  
      RETURN                                                                    
      END                      

                                                                             
      SUBROUTINE RADFRQ_lin(STRA,DIPA,RAKA,PA,THAA,DFR,NFOLD,                  
     & CMP,fr1,RNA,RNB,nr,rdna)   
      parameter (mm=16348)                               
      DIMENSION Rdna(mm),DFR(mm),RNA(nr),RNB(nr)
      PU=3.1415926/180   

c RADFRQ_lin() calculates "averaged" radiation pattern as a function of
c frequency for a given subfault and receiver.
c
c Inputs are:
c		stra = strike of subfault
c		dipa = dip of subfault
c		raka = rake of subfault
c		pa = azimuth from subfault to receiver
c		thaa = take-off angle of ray from subfault to receiver
c
c Radiation pattern value is stored in variable: rdna
c Theoretical radiation pattern is given by variable: rdx
c Conically averaged radiation pattern is given by variable: radvh
c
c Average radiation pattern is calculated over a range of strike,dip,rake,
c takeoff, and azimuth that bracket the theoretical values.  If the range
c for each of these was allowed to vary between 0 - 360 degrees, then the average
c would represent a full spherical average over the entire focal mechanism.
c We restrict the range so that the average is only performed over a subset
c of values that lie around the theoretical values.
c The idea is that the parameters are more likely to be near their theoretical values
c than in some arbitrary orientation.  We refer to this as a conical average,
c that is, the radiation pattern is averaged over a cone of rays that surround
c the theoretical one.
c
c The range is set using the variable: range
c We set the range to be +/- 45 degrees about the theoretical value for each of the
c above parameters.
c
c rad pattern is all average (conical average) above frequency=fr2
c rad pattern is (partly) theoretical below freq=fr1
c
c Below fr1, rdna will always have at least radmin percentage of average rad
c pattern => radmin=0 all theoretical; radmin=1 all average
c
c  rdna = (rdx + (radvh-rdx)*radmin)
c
c Linear taper applies between fr1 and fr2

      fr1=0.5
      fr2=2.0
      radmin = 0.5

C 02/10/2009 Since using a conical average around theoretical ray, don't allow much
C purely theoretical rad pattern
      radmin = 1.0

c get theoretical radtion pattern for this subfault/station
      CALL RDATN(STRA,DIPA,RAKA,PA,THAA,RDSHA,RDSVA) 
                           
      RDX=RDSVA*COS(THAA)*COS(CMP-PA)+RDSHA*SIN(CMP-PA) 
      print*,'cos(th)= ',cos(thaa),' cos(cmp-pa)= ',cos(cmp-pa),' sin(cmp-pa)= ',sin(cmp-pa)

c  RWG RADPAT FIX 03/19/04
c  version 2.02 average radiadtion pattern=>
c
c      rdx = abs(rdsva)*cos(cmp-pa) + abs(rdsha)*sin(cmp-pa)
c
c  Above can create asymmetry if cos or sin is < 0.0
c  Fix by taking abs() after summing SV and SH

      rdx = rdsva*cos(cmp-pa) + rdsha*sin(cmp-pa)
c
c 12/21/2004
c
c OLD: take abs()
c      rdx = abs(rdx)
c
c NEW: preserve sign
      polarity = 1
      if(rdx.lt.0.0) then
	 polarity = -1
         rdx = -rdx
      endif

       range=10.
      RADV=0.0                                                                  
      DO 150 K=1,NR                                                             
         th = thaa + 9*range*pu*(0.5-rand_numb(0))
         fa = pa + 9*range*pu*(0.5-rand_numb(0))

         strX = stra + 9*range*pu*(0.5-rand_numb(0))
         dipX = dipa + 9*range*pu*(0.5-rand_numb(0))
         rakX = raka + 9*range*pu*(0.5-rand_numb(0))

         CALL RDATN(strX,dipX,rakX,fa,th,RDSHA,RDSVA)

c  RWG RADPAT FIX 03/19/04
c  version 2.02 average radiadtion pattern=>
c
c      rads = abs(rdsva)*cos(cmp-fa) + abs(rdsha)*sin(cmp-fa)
c
c  Above can create asymmetry if cos or sin is < 0.0
c  Fix by taking abs() after summing SV and SH

         rads = rdsva*cos(cmp-fa) + rdsha*sin(cmp-fa)

c square the value to remove the sign, polarity is applied later,
c don't forget to take sqareroot after summation
         RADV=RADV+RADS*RADS

  150 CONTINUE                             
                                                  
      RADVH=sqrt(RADV/FLOAT(NR))

      DO 101 I=1,NFOLD                                                       

         if(dfr(i).le.fr1) then
	    del = radmin

         else if(dfr(i).gt.fr1.and.dfr(i).le.fr2) then
            del = alog(dfr(i)/fr1)/alog(fr2/fr1)
	    if(del.lt.radmin) del = radmin

         else
            del = 1.0

         end if

c 12/21/2004 Use polarity to preserve sign
         rdna(i) = polarity*(rdx + (radvh-rdx)*del)

  101 continue

      RETURN                                                                    
      END

      SUBROUTINE RADFRQ_lin5(STRA,DIPA,RAKA,PA,THAA,DFR,NFOLD,                  
     & CMP,fr1,RNA,RNB,nr,rdna)   

c Theoretical rad. coeff
      parameter (mm=16348)                               
                  
      DIMENSION Rdna(mm),DFR(mm),RNA(nr),RNB(nr)
      PU=3.1415926/180   

      fr1=0.5                  
      fr2=2.7

c     fr2=1.7
c     fr2=1.7  used in the calculations for ORI

      radvh=0.5 
       
       CALL RDATN(STRA,DIPA,RAKA,PA,THAA,RDSHA,RDSVA) 
                           
      RDX=RDSVA*COS(THAA)*COS(CMP-PA)+RDSHA*SIN(CMP-PA) 

       range=10.
       tha1=thaa-range*pu
       tha2=thaa+range*pu
      if(tha1.lt.90.*pu) tha1=90.*pu
      if(tha2.gt.180.*pu) tha2=180.*pu 
                              
      RADV=0.0 
      goto 160                                                                 
                                               
      DO 150 K=1,NR                                                             
      TH=ACOS((1.0-RNA(K))*COS(THA1)+RNA(K)*COS(THA2))    
      fa=360.*pu*rnb(k)
                                             
      CALL RDATN(STRA,DIPA,RAKA,FA,TH,RDSHA,RDSVA)                              
      RADS=RDSVA*COS(TH)*COS(CMP-FA)+RDSHA*SIN(CMP-FA)                       
      RADV=RADV+ABS(RADS)                                                    
  150 CONTINUE                             
                      
                            
      RADVH=RADV/FLOAT(NR)
       
160      radvh=0.6 
         rdx=0.6

      DO 101 I=1,NFOLD                                                       
      Rdna(I)=RDX  
      IF(DFR(I).LE.fr1) GO TO 101
      IF(DFR(I).GT.fr1.and.dfr(i).le.fr2) then
      rdna(i)=rdna(i)+(radvh-rdx)*(dfr(i)-fr1)/(fr2-fr1)
      goto 101
      else
      rdna(i)=radvh
      end if

  101 CONTINUE                                               
              
      RETURN                                                                    
      END                         

      SUBROUTINE RADV_lin(STRA,DIPA,RAKA,PA,THAA,DFR,NFOLD,                  
     & fr1,RNA,RNB,nr,rdna)   

c Theoretical rad. coeff
      parameter (mm=16348)                               
                  
      DIMENSION Rdna(mm),DFR(mm),RNA(nr),RNB(nr)
      PU=3.1415926/180   
                     
      fr2=1.5

c     fr2=1.7
c     fr2=1.7  used in the calculations for ORI

      fr1=0.001                  
      fr2=0.01

      radvh=0.7
                                   
      CALL RDATN(STRA,DIPA,RAKA,PA,THAA,RDSHA,RDSVA)                            
      RDX=RDSVA*sin(THAA) 


       range=40.
       tha1=thaa-range*pu
       tha2=thaa+range*pu
      if(tha1.lt.90.*pu) tha1=90.*pu
      if(tha2.gt.180.*pu) tha2=180.*pu 

                               
      RADV=0.0                                                                  
                                               
      DO 150 K=1,NR                                                             
      TH=ACOS((1.0-RNA(K))*COS(THA1)+RNA(K)*COS(THA2))    
      fa=360*pu*rnb(k)
                                             
      CALL RDATN(STRA,DIPA,RAKA,FA,TH,RDSHA,RDSVA)                              
      RADS=RDSVA*sin(TH)                         
      RADV=RADV+ABS(RADS)                                                       
  150 CONTINUE                             
                                                  
      RADVH=RADV/FLOAT(NR)/2.                                                                                     
      DO 101 I=1,NFOLD                                                       
      Rdna(I)=RDX  
      IF(DFR(I).LE.fr1) GO TO 101
      IF(DFR(I).GT.fr1.and.dfr(i).le.fr2) then
      rdna(i)=rdna(i)+(radvh-rdx)*(dfr(i)-fr1)/(fr2-fr1)
      goto 101
      else
      rdna(i)=radvh
      end if

  101 CONTINUE                                               
              
      RETURN                                                                    
      END                                                                                            
                                  
      subroutine highcor(nf,mf,kla,nnn,stdd,rdna)
      parameter (mmv=30000)
      complex*8 cw1(16348)
      dimension rdna(1),stdd(1)
      
      rp=0.63
      prtitn=0.71
      do 110 i=1,nnn
      cw1(i)=(0.,0.)
110   continue
      do 100 i=1,kla
      cw1(i)=cmplx(stdd(i)/rp/prtitn,0.)
100   continue
      call fast(nnn,cw1,-1)
      do 120 i=1,nnn
      cw1(i)=cw1(i)/float(nnn)
120   continue
      do 140 i=1,nf
      cw1(i)=cw1(i)*abs(rdna(i))
 
140   continue
      do 150 i=1,mf
      m=nf+i
      mm=nf-i
      cw1(m)=cw1(m)*abs(rdna(mm))
150   continue
      call fast(nnn,cw1,1)
      do 160 i=1,nnn
      stdd(i)=real(cw1(i))
160   continue
      return
      end          

      subroutine highcor_f(nf,mf,np2,cw1,stdd,rdna)
      complex*8 cw1(1)
      dimension rdna(1),stdd(1)
      
      rp=0.63
      prtitn=0.71

c 12/21/2004 preserve sign of rdna, OLD: use abs(rdna)

      do 140 i=1,nf
      cw1(i)=cw1(i)*rdna(i)
140   continue

      do 150 i=nf+1,nf+mf
      mm=2*nf-i
      cw1(i)=cw1(i)*rdna(mm)
150   continue

c      do 150 i=1,mf
c      m=nf+i
c      mm=nf-i
c      cw1(m)=cw1(m)*abs(rdna(mm))
c150   continue

      call fast(np2,cw1,1)

      fac = 1.0/(rp*prtitn*float(np2))
      do 160 i=1,np2
      stdd(i)=fac*real(cw1(i))
160   continue

      n0=np2/10
      dd = 3.14159625/(float(n0))
      do 852 i=1,n0
         arg = 0.5*(1.0 + cos(i*dd))
         stdd(np2-n0+i)=stdd(np2-n0+i)*arg
852   continue

      return
      end
                                                                   
      SUBROUTINE RDATN(STR,DIP,RAK,AZ,TH,RDSH,RDSV)                             
C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C          
C                                                                               
C     CALCULATION OF THE RADIATION TERM INDUCED BY THE DIRECT BODY WAVE         
C                                                                               
C     STR : STRIKE  (RADIAN)                                                    
C     DIP : DIP ANGLE                                                           
C     RAK : RAKE ANGLE                                                          
C     AZ  : AZIMTH FROM SORCE TO RECEIVER (CLOCKWISE FROM NORTH)                
C     TH  : INCIDENT ANGLE (FROM DOWN)                                          
C     RDSH: RADIATION FOR SH-WAVE                                               
C     RDSV: RADIATION FOR SV-WAVE                                               
C                                                                               
C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C          
      SR=SIN(RAK)                                                               
      CR=COS(RAK)                                                               
      SD=SIN(DIP)                                                               
      CD=COS(DIP)                                                               
      ST=SIN(TH)                                                                
      CT=COS(TH)                                                                
      SS=SIN(AZ-STR)                                                            
      CS=COS(AZ-STR)      
                                                                  
      RDP=CR*SD*ST**2*2*SS*CS-CR*CD*2*ST*CT*CS                                  
     1    +SR*2*SD*CD*(CT**2-ST**2*SS**2)                                       
     2    +SR*(CD**2-SD**2)*2*ST*CT*SS                                          
      RDSV=SR*(CD**2-SD**2)*(CT**2-ST**2)*SS-CR*CD*(CT**2-ST**2)*CS             
     1    +CR*SD*ST*CT*2*SS*CS                                                  
     2    -SR*SD*CD*2*ST*CT*(1+SS**2)  

                                         
       RDSH=CR*CD*CT*SS+CR*SD*ST*(CS**2-SS**2)
     1    +SR*(CD**2-SD**2)*CT*CS                                               
     2    -SR*SD*CD*ST*2*SS*CS  

     
  
     
c      rdsv=sin(rak)*cos(2*dip)*cos(2*th)*sin(az-str)-cos(rak)*cos(dip)*
c     1 cos(2*th)*cos(az-str)+0.5*cos(rak)*sin(dip)*sin(2*th)*
c     2 sin(2*(az-str))-0.5*sin(rak)*sin(2*dip)*sin(2*th)*
c     3   (1+sin(az-str)*sin(az-str))
     
c      rdsh= cos(rak)*cos(dip)*cos(th)*sin(az-str)+cos(rak)*sin(dip)*
c     1 sin(th)*cos(2*(az-str))+sin(rak)*cos(2*dip)*cos(th)*cos(az-str)-
c     2 0.5*sin(rak)*sin(2*dip)*sin(th)*sin(2*(az-str)) 
     
                                                   
C     RDP=RDP/GM**3                                                             
      RETURN                                                                    
      END   
                               
                                                        
      SUBROUTINE INAC(DT,NN,DDY,DY)                                
      dimension DY(1),DDY(1)    
c integrate
                                                    
      DY(1)=DDY(1)*DT                                                           
      DYMAX=0.                                                                  
      DO  110  M=2,NN                                                           
      DY(M)=DY(M-1)+(DDY(M-1)+DDY(M))/2.*DT   
  110 CONTINUE                                                                  
      RETURN                                                                    
      END                   
                                                             
                                                                
      SUBROUTINE FAST(NNN,ACE,IND)                                              

C RWG 20180815
C changed array dimension declaration to ACE(1) to suppress compiler
C warnings.  I think this is OK.
C     COMPLEX*8 ACE(30000),TEMP,THETA                                            

      COMPLEX*8 ACE(1),TEMP,THETA                                            

      J=1                                                                       
      DO 100 I=1,NNN                                                            
      IF(I.GE.J) GO TO 110                                                      
      TEMP=ACE(J)                                                               
      ACE(J)=ACE(I)                                                             
      ACE(I)=TEMP                                                               
  110 M=NNN/2                                                                   
  120 IF(J.LE.M) GO TO 130                                                      
      J=J-M                                                                     
      M=M/2                                                                     
      IF(M.GE.2) GO TO 120                                                      
  130 J=J+M                                                                     
  100 CONTINUE                                                                  
      KMAX=1                                                                    
  140 IF(KMAX.GE.NNN) RETURN                                                    
      ISTEP=KMAX*2                                                              
      DO 150 K=1,KMAX                                                           
      THETA=CMPLX(0.0,3.141593*FLOAT(IND*(K-1))/FLOAT(KMAX))                    
      DO 160 I=K,NNN,ISTEP                                                      
      J=I+KMAX                                                                  
      TEMP=ACE(J)*CEXP(THETA)                                                   
      ACE(J)=ACE(I)-TEMP                                                        
      ACE(I)=ACE(I)+TEMP                                                        
  160 CONTINUE                                                                  
  150 CONTINUE                                                                  
      KMAX=ISTEP                                                                
      GO TO 140       
      END          

        subroutine RANU2(NRR,RN)
        dimension rn(1)        
        do i=1,nrr
        rn(i)=rand_numb(0)
        end do
        return
        end
       
       	SUBROUTINE RANN2(NTOT,ACC)
	dimension ACC(1)
	COMMON/RANDO/ifu1
	
	J=1
	ACC(1)=0.
	X=RAND(ifu1)
	DO 8 N=2,NTOT
	GOTO (1,2),J
1       X1=RAND(0)

964     continue
        if(x1.eq.0.0) then
           X1=RAND(0)
           go to 964
        endif

        X2=RAND(0)
963     continue
        if(x2.eq.0.0) then
           X2=RAND(0)
           go to 963
        endif

	X2=6.2831853*X2
	X1=-ALOG(X1)
	X1=SQRT(X1+X1)
	W=X1*COS(X2)
	J=2
	GO TO 3
2	W=X1*SIN(X2)
	J=1
3	CONTINUE
	ACC(N)=W
8	CONTINUE
	S=0.
	DO 10 I=1,NTOT
	S=S+ACC(I)*ACC(I)
10	CONTINUE
	S=S/NTOT
	DO 11 I=1,NTOT
	ACC(I)=ACC(I)/SQRT(S)
11	CONTINUE
	RETURN
	END
C

	SUBROUTINE FLZERO(N,DT,A)
	DIMENSION A(1)
	VE=.0

	DE=.0
	A1=DT/2.
	A2=A1*DT/3.
	NSTPS=N-1
	DO 1 I=1,NSTPS
	DE=DE+VE*DT+A2*(2.*A(I)+A(I+1))
1	VE=VE+A1*(A(I)+A(I+1))
	RNSTP=NSTPS
	T=RNSTP*DT
	C1=2./T*(VE-3./T*DE)
	C2=6./T*(2./T*DE-VE)/T
	DO 2 I=3,N
	A3=I-1
	A(I)=A(I)+C1+C2*A3*DT
2	CONTINUE
        
	RETURN
	END
     
       subroutine even_dist1(nevent,xlonq,ylatq,slon,slat,azmq,dipangq,
     +        zm,astopq,dx,dy,nx,nw,nnq,rl,ph,th,dst,zet)


c    Calculates:
c        rl(i,j) : Distance form the center of the subfault to the station
c        th(i,j) : Take-off angle 
c        ph(i,j) : Azimuth

	 include 'params.h'
CCC         parameter (nq=600,np=50)

       dimension rl(nq,np),ph(nq,np),th(nq,np),dst(nq,np),zet(nq,np) 
        dimension xlonq(1),ylatq(1),azmq(1),dipangq(1),astopq(1),
     +nnq(1)

      
        pi=3.14159265
        gh=10000.
        gl=0.      
         alei=0.
         alsi=0.

c  calculates the scaling factors for converting degrees into km

         thei=ylatq(1)
    
         do ii=1,2
         if(ii.eq.1 ) then 
         alsi=alei+1
         thsi=thei
         else
         thsi=thei+1
         alsi=alei
         end if
         
         i=0

        call DELAZ5( THEI, ALEI, THSI, ALSI, DELT, DELTDG,
     +dis, AZES, AZESDG, AZSE, AZSEDG, I)
        az=azesdg
        
        x=dis*sin(pi*az/180.)
        y=dis*cos(pi*az/180.)
        if(ii.eq.1) ddx=x
        if(ii.eq.2) ddy=y
        end do

        i=0
	do 20 k=1,nevent

        az=azmq(k)*pi/180.
        dip=dipangq(k)*pi/180

        astop=astopq(k)
        ylat=ylatq(k)
        xlon=xlonq(k)  
        do 20 kk=1,nnq(k)      
        i=i+1
        do 20 j=1,nw
        
        a1=((j-1)*dy+dy/2)*cos(dip)
        b1=((j-1)*dy+dy/2)*sin(dip)
     
        dlon = ((i-0.5)*dx - astop)*sin(az) + a1*cos(az) 
        dlat = ((i-0.5)*dx - astop)*cos(az) - a1*sin(az)

        stlon = xlon + dlon/ddx
        stlat = ylat + dlat/ddy
      
        zm1=zm+b1 

         call DELAZ5( stlat, stlon, slat, slon, DELT, DELTDG,
     +dis, AZES, AZESDG, AZSE, AZSEDG, 0)

        dst(i,j)=dis
        rl(i,j)=sqrt(dis*dis+zm1*zm1)
        th(i,j)=pi-atan2(dis,zm1)
        ph(i,j)=azes
        zet(i,j)=zm1

c        print*,i,j,stlat,stlon,slat,slon,dis

20      continue
        return
        end

       subroutine even_dist2(xlonq,ylatq,slon,slat,azmq,dipangq,
     +        zm,astop,dx,dy,nx,nw,rl,ph,th,dst,zet)


c    Calculates:
c        rl(i,j) : Distance form the center of the subfault to the station
c        th(i,j) : Take-off angle 
c        ph(i,j) : Azimuth

	 include 'params.h'
CCC         parameter (nq=600,np=50)

       dimension rl(nq,np),ph(nq,np),th(nq,np),dst(nq,np),zet(nq,np) 
      
        pi=3.14159265
        gh=10000.
        gl=0.      
         alei=0.
         alsi=0.

c  calculates the scaling factors for converting degrees into km

         thei=ylatq
    
         do ii=1,2
         if(ii.eq.1 ) then 
         alsi=alei+1
         thsi=thei
         else
         thsi=thei+1
         alsi=alei
         end if
         
         i=0

        call DELAZ5( THEI, ALEI, THSI, ALSI, DELT, DELTDG,
     +dis, AZES, AZESDG, AZSE, AZSEDG, I)
        az=azesdg
        
        x=dis*sin(pi*az/180.)
        y=dis*cos(pi*az/180.)
        if(ii.eq.1) ddx=x
        if(ii.eq.2) ddy=y
        end do

        az=azmq*pi/180.
        dip=dipangq*pi/180

        ylat=ylatq
        xlon=xlonq
        do 20 i=1,nx
        do 20 j=1,nw
        
        a1=((j-1)*dy+dy/2)*cos(dip)
        b1=((j-1)*dy+dy/2)*sin(dip)
     
        dlon = ((i-0.5)*dx - astop)*sin(az) + a1*cos(az) 
        dlat = ((i-0.5)*dx - astop)*cos(az) - a1*sin(az)

        stlon = xlon + dlon/ddx
        stlat = ylat + dlat/ddy
      
        zm1=zm+b1 

         call DELAZ5( stlat, stlon, slat, slon, DELT, DELTDG,
     +dis, AZES, AZESDG, AZSE, AZSEDG, 0)

        dst(i,j)=dis
        rl(i,j)=sqrt(dis*dis+zm1*zm1)
        th(i,j)=pi-atan2(dis,zm1)
        ph(i,j)=azes
        zet(i,j)=zm1

c        print*,i,j,stlat,stlon,slat,slon,dis

20      continue
        return
        end


      SUBROUTINE  DELAZ5( THEI, ALEI, THSI, ALSI, DELT, DELTDG,
     2DELTKM, AZES, AZESDG, AZSE, AZSEDG, I )
       DOUBLE  PRECISION C, AK, D, E, CP, AKP, DP, EP,
     2A, B, G, H, AP, BP, GP, HP
      IF(I) 50, 50, 51
C     IF  COORDINATES ARE GEOGRAPH DEG I=0
C     IF COORDINATES ARE GEOCENT RADIAN  I=1
   50 THE=1.745329252E-2*THEI
      ALE=1.745329252E-2*ALEI
      THS=1.745329252E-2*THSI
      ALS=1.745329252E-2*ALSI
      AAA=0.9931177*TAN(THE)
      THE=ATAN(AAA)
      AAA=0.9931177*TAN(THS)
      THS=ATAN(AAA)
      GO TO 32
   51 THE=THEI
      ALE=ALEI
      THS=THSI
      ALS = ALSI
   32 CONTINUE
      C= SIN(THE)
      AK=-COS(THE)
      D=SIN(ALE)
      E= -COS(ALE)
      A= AK*E
      B= -AK*D
      G=-C*E
      H=C*D
      CP=SIN(THS)
      AKP=-COS(THS)
      DP=SIN(ALS)
      EP = -COS(ALS)
      AP = AKP*EP
      BP=-AKP*DP
      GP=-CP*EP
      HP=CP*DP
      C1=A*AP+B*BP+C*CP
      IF( C1-0.94 )  30, 31, 31
   30 IF(C1+0.94) 28, 28, 29
   29 DELT=ACOS(C1)
   33 DELTKM=6371.0*DELT
      C3 = (AP-D)**2+(BP-E)**2+CP**2-2.0
      C4 = (AP-G)**2+(BP-H)**2+(CP-AK)**2-2.0
      C5 = (A-DP)**2+ (B-EP)**2+C**2-2.0
      C6 = (A-GP)**2+(B-HP)**2+(C-AKP)**2-2.0
      DELTDG = 57.29577951*DELT
      AZES = ATAN2(C3, C4 )
      IF ( AZES ) 80, 81, 81
   80 AZES = 6.283185308+ AZES
   81 AZSE = ATAN2( C5, C6 )
      IF ( AZSE ) 70, 71 , 71
   70 AZSE=6.283185308+AZSE
   71 AZESDG=57.29577951*AZES
      AZSEDG=57.29577951*AZSE
      RETURN
   31 C1=(A-AP)**2+(B-BP)**2+(C-CP)**2
      C1= SQRT(C1)
      C1=C1/2.0
      DELT = ASIN(C1)
      DELT= 2.0*DELT
      GO TO 33
   28 C1=(A+AP)**2+(B+BP)**2+(C+CP)**2
      C1 = SQRT(C1 )
      C1= C1/2.0
      DELT = ACOS(C1)
      DELT = 2.0*DELT
      GO TO 33

      END
      

      FUNCTION DGAMM(X)                                                         
C *****************************************************                           
C *                   GAMMA FUNCTION                  *                           
C *                IN DOUBLE PRECISION                *                           
C *      COPYRIGHT : M.MORI  JUNE 30 1989  V.1        *                           
C *****************************************************                           
      IMPLICIT REAL*8 (A-H,O-Z)                                                 

      DIMENSION C(0:19)                                                         
      DATA IN / 19 /                                                            

C     ---- FOR SINGLE PRECISION ----                                            
C     DATA IN / 10 /                                                            

      DATA C /  1.0                   D0,                                       
     1         -0.42278 43350 98467 1 D0,                                       
     2         -0.23309 37364 21786 7 D0,                                       
     3          0.19109 11013 87691 5 D0,                                       
     4         -0.24552 49000 54000 2 D-1,                                      
     5         -0.17645 24455 01443 2 D-1,                                      
     6          0.80232 73022 26734 7 D-2,                                      
     7         -0.80432 97756 04247 0 D-3,                                      
     8         -0.36083 78162 548     D-3,                                      
     9          0.14559 61421 399     D-3,                                      
     1         -0.17545 85975 17      D-4,                                      
     1         -0.25889 95022 4       D-5,                                      
     2          0.13385 01546 6       D-5,                                      
     3         -0.20547 43152         D-6,                                      
     4         -0.15952 68            D-9,                                      
     5          0.62756 218           D-8,                                      
     6         -0.12736 143           D-8,                                      
     7          0.92339 7             D-10,                                     
     8          0.12002 8             D-10,                                     
     9         -0.42202               D-11 /                                    

      IF (X .GT. 57.0D0) GO TO 901                                              

      XX = X                                                                    
      IF (XX .LE. 1.5D0) THEN                                                   
        IF (XX .GE. 0.5D0) THEN                                                 
          A = XX - 1.0D0                                                        
          FCTR = 1.0D0                                                          
        ELSE                                                                    
          M = INT(XX)                                                           
          A = XX - M                                                            
          IF (A .EQ. 0.0D0) THEN                                                
            GO TO 902                                                           
          ELSE IF (A .GE. -0.5D0) THEN                                          
            MG = IABS(M) + 1                                                    
          ELSE                                                                  
            MG = IABS(M) + 2                                                    
            A = A + 1.0D0                                                       
          END IF                                                                
          Z = 1.0D0                                                             
          DO 20 I = 1, MG                                                       
            Z = Z * XX                                                          
            XX = XX + 1.0D0                                                     
   20     CONTINUE                                                              
          FCTR = 1.0D0 / Z                                                      
        END IF                                                                  

      ELSE                                                                      
        M = INT (XX)                                                            
        A = XX - M                                                              
        IF (A .LE. 0.5D0) THEN                                                  
          MG = M - 1                                                            
        ELSE                                                                    
          MG = M                                                                
          A = A - 1.0D0                                                         
        END IF                                                                  
        Z = 1.0D0                                                               
        DO 30 I = 1, MG                                                         
          Z = Z * (XX - 1.0D0)                                                  
          XX = XX - 1.0D0                                                       
   30   CONTINUE                                                                
        FCTR = Z                                                                
      END IF                                                                    

      Y = C(IN)                                                                 
      DO 10 I = IN - 1, 0, -1                                                   
        Y = C(I) + A * Y                                                        
   10 CONTINUE                                                                  

      DGAMM = FCTR / ((1.0D0 + A) * Y)                                          
      RETURN                                                                    

  901 CONTINUE                                                                  
      WRITE (6,2001) X                                                          
 2001 FORMAT (' (FUNC.DGAMM) X(=',D23.16,')',                                   
     \        ' MUST BE SMALLER THAN 57.0')                                     
      DGAMM = 1.0D75                                                            
      RETURN                                                                    

  902 CONTINUE                                                                  
      WRITE (6,2002) X                                                          
 2002 FORMAT (' (FUNC.DGAMM) INVALID ARGUMENT',                                 
     \        ' X =',D23.16)                                                    
      DGAMM = 1.0D75                                                            
      RETURN                                                                    
      END        

       subroutine aver_remov(nms,s)
       dimension s(1)     
    
      SMA=0    
      DO 73 I=1,nms 
      SMA=SMA+s(I)                                                           
 73   CONTINUE            
      SMA=SMA/nms   
      DO 74 I=1,nms   
      s(I)=s(I)-SMA 
 74   CONTINUE
      return
      end


c
        subroutine filter3d(nbut,ift,flo,fhi,a,n,dt)
        parameter (nphas=0)
        dimension a(3,30000)
        common/butt/yy(30000),xx(30000)

        do k=1,3
        do i=1,n
        xx(i)=a(k,i)
        end do
        call zpass(n,dt,ift,flo,fhi,nphas,nbut)
        do i=1,n
        a(k,i)=xx(i)
        end do
        end do
c
      return
      end
      
      subroutine zpass(n,dt,ift,flo,fhi,nphas,nbut)
c
      common/butt/x(30000),z(30000)
      do 1000 i=1,n
 1000 x(i)=z(i)
      tupi=2.0*3.141592654
      iflag=0
      znyq=1.0/(2.0*dt)
      go to (20,30,40), ift+1
  20  a=1.0
      f1=flo
      w0=tupi*f1
      go to 50
  30  a=-1.0
      f1=fhi
      w0=tupi*(znyq-f1)
      go to 50
  40  a=1.0
      f1=flo
      f2=fhi
      w0=tupi*f2
      w1=tupi*(znyq-f1)
      iflag=1
  50  continue
      iq=(1+(-1)**nbut)/2
      n2b=nbut/2
      do 70 i=1,n2b
      ct=cs2(w0,nbut,tupi,iq,i)
      ict=i
      call butter (a,ct,w0,ict,dt,n)
  70  continue
      if (iq.eq.0) call bttr2 (a,w0,dt,n)
      if(nphas.eq.1) go to 88
      call revers (n)
      do 80 i=1,n2b
      ct=cs2(w0,nbut,tupi,iq,i)
      ict=i
      call butter (a,ct,w0,ict,dt,n)
  80  continue
      if (iq.eq.0) call bttr2 (a,w0,dt,n)
      call revers (n)
  88  if (iflag.eq.0) go to 130
      a=-1.0
      w0=w1
      iflag=0
      go to 50
 130  continue
      do 131 i=1,n
 131  z(i)=x(i)
      return
      end
c
      function cs2(w0,nbut,tupi,iq,i)
c
      cs2=2.0*w0*cos(float(2*i-iq)*tupi/float(4*nbut))
      return
      end
c
      subroutine butter (a,b,w,ict,dt,n)
c
      common/butt/x(30000),z(30000)
      za=2.0*a
      call set (w,dt,n,wp,ak)
      wp2=wp*wp
      b0=wp2
      b1=1.0+b*ak+wp2
      b2=za*(wp2-1.0)
      b3=1.0-b*ak+wp2
      if (a.gt.0.0) go to 20
      yt1=0.0
      yt2=0.0
      if (ict.gt.1) go to 10
      xt1=x(1)
      xt2=x(1)
      go to 30
  10  xt1=0.0
      xt2=0.0
      go to 30
  20  yt1=x(1)
      yt2=x(1)
      xt1=x(1)
      xt2=x(1)
  30  continue
      xt1=0.0
      xt2=0.0
      yt1=0.0
      yt2=0.0
      do 40 i=1,n
      yt3=(b0*(x(i)+za*xt2+xt1)-b2*yt2-b3*yt1)/b1
      xt1=xt2
      xt2=x(i)
      x(i)=yt3
      yt1=yt2
      yt2=yt3
  40  continue
      return
      end
c
      subroutine bttr2 (a,w0,dt,n)
c
      common/butt/x(30000),z(30000)
      call set (w0,dt,n,wp,ak)
      b0=wp+1.0
      b1=a*(wp-1.0)
      if (a.gt.0.0) go to 10
      yt=0.0
      xt=0.0
      go to 20
  10  yt=x(1)
      xt=x(1)
  20  continue
      xt=0.0
      yt=0.0
      do 30 i=1,n
      yt=(wp*(x(i)+a*xt)-b1*yt)/b0
      xt=x(i)
      x(i)=yt
  30  continue
      return
      end
c
      subroutine set (w,dt,n,wp,ak)
c
      wp=tan(w*dt/2.0)
      ak=wp/w
      return
      end
c
      subroutine revers (n)
c
      common/butt/x(30000),z(30000)
      common/larry/y(30000)
      do 1 i=1,n
    1 y(i)=x(n-i+1)
      do 2 i=1,n
    2 x(i)=y(i)
      return
      end
c

      subroutine get_sitefacs(j0,nfreq,fn,an)
      include 'params.h'
      common/vmod/dep(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),dn(nlaymax),qp(nlaymax),qs(nlaymax)
      real*8 dep,th,vp,vs,dn
      real*4 qp,qs
      dimension fn(1),an(1)

      vdsrc = vs(j0)*dn(j0)

      do 6140 kf=1,nfreq

	 stt = 0.25/exp(fn(kf))

	 i = 2
	 zdep = 0.0
	 pz = 0.0
	 tt = 0.0
	 ttp = th(2)/vs(2)

6145     if(ttp.ge.stt.or.i.eq.j0) go to 6146
            
	 zdep = zdep + th(i)
	 pz = pz + dn(i)*th(i)/vs(i)
	 tt = ttp

	 i = i + 1
	 ttp = th(i)/vs(i) + tt

         go to 6145
6146     continue

         bz = (zdep + (stt - tt)*vs(i))/stt
         pz = (pz + dn(i)*(stt - tt))/stt

cc	 print*,'bz=',bz,' pz=',pz,' zd=',zd,' stt=',stt,' tt=',tt

	 an(kf) = 0.5*alog((vdsrc/(bz*pz)))

6140  continue

cPPP      nn = nfreq

      return
      end

      subroutine get_sitefacsOLD(j0,nn,fn,an)
      include 'params.h'
      common/vmod/dep(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),dn(nlaymax),qp(nlaymax),qs(nlaymax)
      real*8 dep,th,vp,vs,dn
      real*4 qp,qs
      dimension fn(1),an(1)
      dimension fz(nlaymax),az(nlaymax)

      vdsrc = vs(j0)*dn(j0)
      nz = 5

      kf = 1
      ttm = 0.0
      dpm = 0.0
      rtm = 0.0

      do 6140 i=2,j0

         dz = (dep(i)-dpm)/nz

         do 6141 j=1,nz

	    zth = j*dz

	    zz = dpm + zth
	    tti = zth/vs(i)
	    tt = ttm + tti

	    bz = zz/tt
	    rz = (rtm + dn(i)*tti)/tt

	    fz(kf) = alog(0.25/tt)
	    az(kf) = 0.5*alog((vdsrc/(bz*rz)))

	    kf = kf + 1

6141     continue

	 ttm = tt
	 dpm = zz
	 rtm = rz*tt

6140  continue

      nn = kf-1
      do 6142 i=1,nn

	 fn(i) = fz(nn+1 - i)
	 an(i) = az(nn+1 - i)

6142  continue

      return
      end

      subroutine siteamp(np2,cw,dfr,nn,fn,an)
      dimension fn(1),an(1),dfr(1)
      COMPLEX*8 cw(1) 

      np = np2/2
      nf = np+1

      fac = 1.0
      ap = 1.0

      kn = 1

      fm = 0.0
      am = an(kn)

      fp = fn(kn)
      ap = an(kn)

      cw(1) = cw(1)*an(1)

      do 2100 i=2,np
	 freq = alog(dfr(i))

         if(freq.gt.fp.and.kn.le.nn) then
9123        fm = fp
	    am = ap

	    kn = kn + 1

	    if(kn.gt.nn) then
	       fp = 1.0e+15
               ap = an(nn)
            else
               fp = fn(kn)
               ap = an(kn)
            endif

         if(freq.gt.fp.and.kn.le.nn) go to 9123

         endif

	 fac = exp(am + (freq-fm)*(ap-am)/(fp-fm))
         cw(i) = cw(i)*fac
2100  continue

      do 2101 i=1,np-1
         cw(np2-i+1) = conjg(cw(i+1))
2101  continue

      cw(nf) = cw(nf)*an(nn)

      return
      end

      subroutine gf_amp_tt(j0,src_depth,range,sgc,itype,md,rp0,stime,rpath,qbar)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      common/rays/nh(1,nlaymax),nm(1,nlaymax),ndeg(nlaymax),nd(nlaymax)
      real*8 dpt,th,vp,vs,rho
      real*4 qp,qs
      complex *16 rcv,cp0,gc
      complex sgc
      real *4 src_depth,range,stime,rpath,rp0,qbar
    
      krec = 2
      ir = 1
      ndeg(ir)=1
      hr=th(1)
      hs=src_depth
      rr=range
      gc = (1.0d+00,0.0d+00)

c find source layer, set to 'ksrc'
      hs_tol = 0.02
      dep = 0.0
      do 1543 ksrc=1,j0
         dep = dep + th(ksrc)

         if(hs.ge.dep.and.(hs-dep).lt.hs_tol) hs = dep + hs_tol
         if(hs.lt.dep) then
            if((dep-hs).lt.hs_tol) hs = dep - hs_tol
            goto 1542
         endif

cc         if(hs.eq.dep) hs = hs + 0.001
cc         if(hs.lt.dep) goto 1542

1543  continue
1542  continue

c Set up ray path description:
c For simplicity, we only track 1 ray at a time, therefore ir=1
c
c    nh() is array containing the layer indices for each ray-path segment
c    nm() is array specifying the mode (p,sv,sh) for each segment
c    nd() contains the total number of ray-path segments
c
c if itype=1,3,5... (odd): upgoing ray, > 1 means Moho multiples
c if itype=2,4,6,...(even): down-going, then Moho reflected ray
c
c Assume Moho is above layer with 0.00 thickness or the deepest layer in the input model

c upgoing ray
      if(mod(itype,2).eq.1) then
      
         l = 0
         do 1922 j=ksrc,krec,-1
	    l = l+1
	    nh(ir,l) = j
	    nm(ir,l) = md

1922     continue
      
c now loop over Moho multiples, if any
	 ktn = (itype-1)/2
         do 5923 kt=1,ktn

            do 3923 j=krec,j0-1
	       l = l+1
	       nh(ir,l) = j
	       nm(ir,l) = md
	       if(th(j+1).eq.0.0) goto 3925

3923        continue
3925        continue

	    kbot = j

            do 3924 j=kbot,krec,-1
      	       l = l+1
   	       nh(ir,l) = j
   	       nm(ir,l) = md

3924        continue

5923     continue
         nd(ir) = l
c end of upgoing ray

c down-going ray
      else
      
         l = 0
         do 2923 j=ksrc,j0-1
	    l = l+1
	    nh(ir,l) = j
	    nm(ir,l) = md
	    if(th(j+1).eq.0.0) goto 2925

2923     continue
2925     continue

	 kbot = j

      print*,'hs,ksrc,j,j0,kbot= ',hs,ksrc,j,j0,kbot

         do 2924 j=kbot,krec,-1
	    l = l+1
	    nh(ir,l) = j
	    nm(ir,l) = md

2924     continue
      
c now loop over Moho multiples, if any
	 ktn = (itype-2)/2
         do 5926 kt=1,ktn

            do 4926 j=krec,j0-1
	       l = l+1
	       nh(ir,l) = j
	       nm(ir,l) = md
	       if(th(j+1).eq.0.0) goto 4927

4926        continue
4927        continue

	    kbot = j

            do 4928 j=kbot,krec,-1
	       l = l+1
	       nh(ir,l) = j
	       nm(ir,l) = md

4928        continue

5926     continue
         nd(ir) = l

      endif

      call trav(ir,hs,hr)
      call pnot(ir,p0,t0,rr)
      call ttime(ir,p0,t0,p1,t1,rr)
      rp0 = p0
      stime = t0

c calculate geometric spreading by 1./(ray path length)
c
      call geom_terms(hs,p0,itype,rpd,qbar)
      rpath = rpd

      print*,'p0= ', p0

      return
      end
c
      function cagcon(p,ir,r)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      real*4 alp,als
c          computes complex time as a function of complex
c          ray parameter.
      complex*16  cagcon,cr,p,a,ea,eb
      common/travel/alp(nlaymax),als(nlaymax),nd,nup
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      real*8 dpt,th,vp,vs,rho,r
      real*4 qp,qs
      a=(0.0d+00,0.0d+00)
      do 10 i=1,nd
      ea=(0.0d+00,0.0d+00)
      eb=(0.0d+00,0.0d+00)
      if(alp(i).gt.0.) ea=cr(p,vp(i))
      if(als(i).gt.0.) eb=cr(p,vs(i))
10    a=a + ea*dble(alp(i))*th(i) + eb*dble(als(i))*th(i)
      cagcon=p*r + a
      return
      end
c
      function cr(p,v)
      implicit real*8 (a-h,o-z)
      complex*16  cr,p
      t1=1.0d-08
      rsq=1.0d+00/(v*v)
      pr=p
      pi=dimag(p)
      a=rsq -pr*pr + pi*pi
      b=-2.0d+00*pi*pr
      d=dsqrt(dsqrt(a*a + b*b))
      if(dabs(pi).lt.t1) go to 10
      phi=datan2(b,a)
      go to 11
10    continue
      phi=0.0d+00
      if(a.lt.0.) phi=3.141592654d+00
11    continue
      e=dcos(phi/2.0d+00)
      f=dsin(phi/2.0d+00)
      if(f.gt.t1) go to 13
      if(e.gt.0.0d+00)  go to 12
13    continue
      e=-e
      f=-f
12    continue
      a=d*e
      b=d*f
      cr=dcmplx(a,b)
      return
      end
c
      function dstdps(p0,ir)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      real*4 alp,als
      complex*16  p,cr,ea,eb,a,b,c
      common/travel/alp(nlaymax),als(nlaymax),nd,nnup
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      real*4 qp,qs
c        calculates second derivative of t with respect to
c        p for first motion approximation.
      a=(0.0d+00,0.0d+00)
      p=dcmplx(p0,0.0d+00)
      do 10 i=1,nd
      b=(0.0d+00,0.0d+00)
      c=(0.0d+00,0.0d+00)
      ea=(0.0d+00,0.0d+00)
      eb=(0.0d+00,0.0d+00)
      if(alp(i).eq.0.) go to 11
      ea=cr(p,vp(i))
      b=-th(i)*dble(alp(i))/(ea*ea*ea*vp(i)*vp(i))
11    continue
      if(als(i).eq.0.) go to 12
      eb=cr(p,vs(i))
      c=-th(i)*dble(als(i))/(eb*eb*eb*vs(i)*vs(i))
12    continue
      a=a+b+c
10    continue
      dstdps=dreal(a)
c     write(13,100) dreal(a)
15    continue
      return
      end
c
      function dtdp(p,ir,r)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      complex*16 dtdp,p,a,b,c,cr,ea,eb
      real*4 alp,als
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      common/travel/alp(nlaymax),als(nlaymax),nd,nup
      real*4 qp,qs
      a=(0.0d+00,0.0d+00)
      do 10 i=1,nd
      b=(0.0d+00,0.0d+00)
      c=(0.0d+00,0.0d+00)
      ea=(0.0d+00,0.0d+00)
      eb=(0.0d+00,0.0d+00)
      if(alp(i).eq.0.) go to 11
      ea=cr(p,vp(i))
      b=th(i)*dble(alp(i))/ea
11    continue
      if(als(i).eq.0.) go to 12
      eb=cr(p,vs(i))
      c=th(i)*dble(als(i))/eb
12    continue
      a=a+b+c
10    continue
      dtdp=r-p*a
      return
      end
c
      subroutine pnot(ir,p0,t0,r)
      include 'params.h'
      implicit real*8 (a-h,o-z)
c     real*8 dpt,th,vp,vs,rho,v,pn,pp,a,r,p0,t0
      real*4 alp, als
      common/travel/alp(nlaymax),als(nlaymax),ndp,nnup
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      complex*16  t,cagcon,dtdp,p
      real*4 qp,qs
c          finding the closest branch cut(i.e. highest velocity)
      v=0.
      do 10 i=1,ndp
      if(alp(i).gt.0.) v=dmax1(v,vp(i))
      if(als(i).gt.0.) v=dmax1(v,vs(i))
10    continue

c   ***RWG
c   fix for case when p0 is very close to 1/v
c
c   OLD
c      eps=1.0d-08
c      p=1.0/v - eps
c
c   NEW

      eps = 1.0d-20
      eps = 1.0d-10
      ptest = 1.0/v

222   jtest = 0
      rp = (ptest - eps)*v
      if(rp.ge.1.0d+00) then
	 eps = 10.0*eps
	 jtest = 1
      endif
      if(jtest.eq.1) go to 222

c add another factor of 10 just to be sure-> problems on Linux
      p = ptest - 10.0*eps

        print*,'eps= ',eps,' p= ',p,' ptest= ',ptest,' rp= ',rp

c   end NEW

      a=dtdp(p,ir,r)
      if(a.lt.0.) go to 11
12    continue
      p0=p
      t=cagcon(p,ir,r)
      t0=t
      return
11    continue
      k=0
      pn=p
      pp=0.
13    k=k+1
      p=(pn+pp)/2.0
      a=dtdp(p,ir,r)
      if((dabs(a).le.0.01).or.(k.ge.40)) go to 12
      if(a.gt.0.) go to 14
      pn=p
      go to 13
14    pp=p
      go to 13
      end
c
      subroutine trav(ir,hs,hr)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      real*4 alp,als
      common/rays/nh(1,nlaymax),nm(1,nlaymax),ndeg(nlaymax),nd(nlaymax)
      common/travel/alp(nlaymax),als(nlaymax),ndeep,nup
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      real*8 dpt,th,vp,vs,rho
      real*4 qp,qs
      common/coff/it(nlaymax),nup1(nlaymax)
      common/rmode/love
      love=1
      if(nm(ir,1).eq.4) love=2
      n=nd(ir)
      do 10 i=1,100
      alp(i)=0.
10    als(i)=0.
      do 11 i=1,n
      if(nm(ir,i).eq.5) alp(nh(ir,i))=alp(nh(ir,i))+1
      if((nm(ir,i).eq.3).or.(nm(ir,i).eq.4)) als(nh(ir,i))=als(nh(ir,i))
     * + 1
11    continue
c          determining ray direction from source.  nup=1, up going;
c          =-1, down going.   lis=source layer, lir=receiver  layer
c          if ndeg(ir) is lt 0, then the ray is up going. this is to
c          remove the ambiguity when source and receiver are in the
c          same layer.
      lis=nh(ir,1)
      lir=nh(ir,n)
      nl=1
      do 12 i=1,n
12    if(nh(ir,i).eq.lis) nl=nl+1
      nup=(-1)**nl
      if(lir.gt.lis) nup=-nup
      if(ndeg(ir).lt.0) nup=+1
      if((n.eq.1).and.(hr.ge.hs)) nup=-1
100   format(1x,'lis=',i5,5x,'nup=',i5)
c          finding interaction types at each interface and direction
c          of each ray segment. nup1 = 1, up ; -1, down.
c          it=0, transmission; =1, reflection; =2, direct ray.
      n1=n-1
      nup1(1)=nup
      if(n.eq.1) go to 20
      do 19 i=1,n1
      k=nh(ir,i)
      m=nh(ir,i+1)
      if(m.eq.k) it(i)=1
      if(m.ne.k) it(i)=0
      if((nup1(i).eq.+1).and.(it(i).eq.1)) nup1(i+1)=-1
      if((nup1(i).eq.-1).and.(it(i).eq.1)) nup1(i+1)=+1
      if((nup1(i).eq.+1).and.(it(i).eq.0)) nup1(i+1)=+1
      if((nup1(i).eq.-1).and.(it(i).eq.0)) nup1(i+1)=-1
19    continue
20    continue
      if(n.eq.1) it(1)=2
c     finding receiver position
      lir1=lir-1
      thtot=0.
      do 22 i=1,lir1
22    thtot=th(i) + thtot
      hrl=hr - thtot
c     adjusting multiplier for receiver layer
c     to include depth
      a1=hrl/th(lir)
      a2=(th(lir)-hrl)/th(lir)
      nupa=nup1(n)
      if(nm(ir,n).eq.5) go to 23
      if((nm(ir,n).eq.3).or.(nm(ir,n).eq.4)) go to 24
23    continue
      if(nupa.eq.1) alp(lir)=alp(lir) - a1
      if(nupa.eq.-1) alp(lir)=alp(lir) - a2
      go to 25
24    continue
      if(nupa.eq.1) als(lir)=als(lir) - a1
      if(nupa.eq.-1) als(lir)=als(lir) - a2
25    continue
c     finding source position
      lis1=lis-1
      thtot=0.
      do 13 i=1,lis1
13    thtot=th(i)+thtot
      hsl=hs-thtot
c     adjusting multiplier for source layer
c     to include depth.
      a1=hsl/th(lis)
      a2=(th(lis)-hsl)/th(lis)
      if(nm(ir,1).eq.5) go to 14
      if((nm(ir,1).eq.3).or.(nm(ir,1).eq.4)) go to 15
14    continue
      if(nup.eq.1) alp(lis)=alp(lis)-a2
      if(nup.eq.-1) alp(lis)=alp(lis)-a1
      go to 16
 15    continue
      if(nup.eq.1) als(lis)=als(lis)-a2
      if(nup.eq.-1) als(lis)=als(lis)-a1
16    continue
c     finding deepest layer that the ray penetrates
      ndeep=0
      do 17 i=1,n
17    ndeep=max0(ndeep,nh(ir,i))
      return
      end
c
      subroutine ttime(ir,p0,t0,p1,t1,r)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      common/coff/it(nlaymax),nup1(nlaymax)
      common/rays/nh(1,nlaymax),nm(1,nlaymax),ndeg(nlaymax),nd(nlaymax)
      real*8 p0,t0,p1,t1,dpt,th,vp,vs,rho,r,va,vb
      real*4 qp,qs
      complex*16 p,t,cagcon
      n=nd(ir)
      p1=p0
      do 99 i=1,n
      nup=nup1(i)
      vb=vs(nh(ir,i))
      va=vb
      if(nm(ir,1).ne.4) va=vp(nh(ir,i))
      p1=dmin1(p1,1.0/va,1.0/vb)
      if(i.eq.n) go to 99
      if(it(i).eq.0) go to 99
      if(nup.eq.1) go to 10
      vb=vs(nh(ir,i) + 1)
      va=vb
      if(nm(ir,1).ne.4) va=vp(nh(ir,i) + 1)
      go to 11
10    vb=vs(nh(ir,i) - 1)
      va=vb
      if(nm(ir,1).ne.4) va=vp(nh(ir,i) - 1)
11    p1=dmin1(p1,1.0/va,1.0/vb)
99    continue
      p=p1
      t=cagcon(p,ir,r)
c      if(prntc) write(13,100) p,t
100   format(1x,'p=',2d15.8,5x,'t=',2d15.8)
      t1=t
c      if(prntc) write(13,101) p1,t1
101   format(1x,'p1=',d15.8,10x,'t1=',d15.8)
      return
      end

      subroutine prduct(ir,p,gc)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      complex*16 p,gencof,gc
      common/coff/it(nlaymax),nup1(nlaymax)
      common/rays/nh(1,nlaymax),nm(1,nlaymax),ndeg(nlaymax),nd(nlaymax)
      common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth
      logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth
c        computes product of generalized reflection and
c        transmission coefficients, then multiples in
c        appropriate receiver functions.
      n=nd(ir)
      if(n.eq.1) go to 20
      n1=n-1
      gc=(1.,0.)
      do 10 i=1,n1
      nup=nup1(i)
      ik=it(i)
      ml=nh(ir,i)
      mr=nm(ir,i)
      kl=nh(ir,i+1)
      kr=nm(ir,i+1)
      gc=gc*gencof(mr,kr,ml,kl,nup,ik,p)
10    continue
      if(prntf) write(13,100) gc
100   format(1x,'gc=',2d15.8)
20    continue
      if(n.eq.1) gc=(1.,0.)
      c=ndeg(ir)
      gc=gc*abs(c)
      return
      end
c
      function gencof(mr,kr,ml,kl,nup,it,p)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      complex*16  gencof,p
      complex*16 tpp,tps,tsp,tss,rpp,rps,rsp,rss
      common/vmod/dpt(nlaymax),th(nlaymax),c(nlaymax),s(nlaymax),d(nlaymax),qp(nlaymax),qs(nlaymax)
      common/rmode/love
      common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth
      logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth
      real*4 qp,qs
      m=ml
      if(it.eq.1) go to 10
c     transmission
      call tranm(p,c(ml),s(ml),d(ml),c(kl),s(kl),d(kl),tpp,tps,tsp,tss)
      k=kl
      if(nup.eq.-1) go to 11
      tps=-tps
      tsp=-tsp
   11 continue
      if(kr.eq.mr) go to 12
      gencof=tps
      if(kr.ne.5) go to 20
      gencof=tsp
      go to 20
   12 continue
      gencof=tss
      if(kr.ne.5) go to 20
      gencof=tpp
      go to 20
   10 continue
c     reflection
      k=kl-1
      if(nup.eq.-1) k=kl+1
      call refft(p,c(ml),s(ml),d(ml),c(k),s(k),d(k),rpp,rps,rsp,rss)
      if(k.gt.ml) go to 13
      rps=-rps
      rsp=-rsp
   13 continue
      if(kr.eq.mr) go to 14
      gencof=rps
      if(kr.ne.5) go to 20
      gencof=rsp
      go to 20
   14 continue
      gencof=rss
      if(kr.ne.5) go to 20
      gencof=rpp
   20 continue
      if(.not.prntf) return
      write(13,100) c(m),s(m),d(m),c(k),s(k),d(k)
100   format(1x,'cb=',d15.8,2x,'sb=',d15.8,2x,'db=',d15.8,2x,'ca=',d15.8
     *,2x,'sa=',d15.8,2x,'da=',d15.8)
      write(13,101) tpp,tps,tsp,tss
      write(13,102) rpp,rps,rsp,rss
101   format(1x,'tpp=',2d10.3,2x,'tps=',2d10.3,2x,'tsp=',2d10.3,2x,'tss=
     *',2d10.3)
102   format(1x,'rpp=',2d10.3,2x,'rps=',2d10.3,2x,'rsp=',2d10.3,2x,'rss=
     *',2d10.3)
      return
      end

c Free surface coefficients for Displacements from p. 140-141 Aki & Richards
c SH is always 2 for horizontal
c
      subroutine recfnct(mr,kr,p0,rcv,ntype)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      complex*16  rcv,ca,cb,cfac1,cfac2,cs2p,cs2s
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      real *8 p0,gb,ga,rfac1,rfac2,cosj,sinj,cosi,sini
      real*4 qp,qs

      if(mr.eq.3) then
	 ga = vp(kr)*p0
	 gb = vs(kr)*p0

	 ga = ga*ga
	 gb = gb*gb

	 rfac1 = (1.0 - 2.0*gb)/(vs(kr)*vs(kr))
	 rfac2 = rfac1*rfac1

	 if(ga.gt.1.0) then
	    ca = cmplx(0.0,sqrt(ga-1.0))
	 else
	    ca = cmplx(sqrt(1.0-ga),0.0)
         endif

	 cb = cmplx(sqrt(1.0-gb),0.0)

	 cfac1 = (4.0*p0*p0/(vp(kr)*vs(kr)))*ca*cb
	 cfac2 = (4.0*p0*rfac1/vp(kr))*cb

	 cs2p = cfac2/(rfac2 + cfac1)
	 cs2s = (rfac2 - cfac1)/(rfac2 + cfac1)

	 sinj = vs(kr)*p0
	 cosj = sqrt(1.0-sinj*sinj)

	 sini = vp(kr)*p0
	 if(sini.gt.1.0) then
	    sini = 0.0
	    cosi = 1.0
	 else
	    cosi = sqrt(1.0-sini*sini)
	 endif

CC      print *,'cs2p= ',cs2p,' cs2s= ',cs2s

c These are the full theoretical ver and rad partitioning formulas (AR p. 141)
c
	 if(ntype.eq.1) rcv = sinj*(1.0 - cs2s) + cosi*cs2p
	 if(ntype.eq.2) rcv = cosj*(1.0 + cs2s) + sini*cs2p
c
c These put all s2s on horizontal and all s2p on vertical
c
CC	 if(ntype.eq.1) rcv = cs2p
CC	 if(ntype.eq.2) rcv = 1.0 + cs2s
      else if(mr.eq.4) then
	 if(ntype.eq.1) rcv = 0.0
	 if(ntype.eq.2) rcv = 2.0
      endif

CC      print *,'rcv= ',rcv

      return
      end

      subroutine refft(p,v3,s1,d1,v4,s2,d2,rpp,rps,rsp,rss)
      implicit real*8 (a-h,o-z)
      common/rmode/love
      real*8 k1,k2,k3,k4,v1,v2,v3,v4,s1,s2,d1,d2
      complex*16  e1,e2,e1p,e2p,c1,c2,c3,c4,c5,c6,ap,bp,cr,bt
      complex*16  rpp,aps,bps,rps,a,b,rss,asp,bsp,rsp,p
      v1=v3
      v2=v4
      if(love.ne.2) go to 11
      v1=s1
      v2=s2
   11 continue
      k4=s2**2*d2/(s1**2*d1)
      b1=.5/(1-k4)
      b2=.5*k4/(k4-1)
      k1=b1/s1**2
      k2=b2/s2**2
      k3=k1+k2
      e1=cr(p,v1)
      e2=cr(p,v2)
      e1p=cr(p,s1)
      e2p=cr(p,s2)
      if(love.eq.2) go to 10
      c1=(p**2)*(k3-p**2)**2
      c2=p**2*e1*e1p*e2p
      c3=(e1*e1p)*(k2-p**2)**2
      c4=e2p*(k1-p*p)**2
      c5=k1*k2*e1*e2p
      c6=k1*k2*e1p
      ap=c1+c3-c5
      bp=c2+c4-c6
      a=-c1+c3-c5
      b=-c2+c4-c6
      bt=ap+bp*e2
      rpp=(a-b*e2)/bt
      aps=2.*p*e1 *(k2-p*p)*(k3-p*p)
      bps=2.*p*e1*(k1-p*p)*e2p
      rps=(aps-bps*e2)/bt
      a=-c1 +c3 +c5
      b=-c2 +c4 +c6
      rss=(a-b*e2)/bt
      asp=2.*p*e1p*(k2-p*p)*(k3-p*p)
      bsp=2.*p*e1p*(k1-p*p)*e2p
      rsp=-(asp-bsp*e2)/bt
      return
 10   r1=d1*s1**2
      r2=d2*s2**2
      rss=(r1*e1p-r2*e2p) / (r1*e1p+r2*e2p)
      rpp=0.
      rps=0.
      rsp=0.
      return
      end
c
      subroutine tranm(p,v3,s1,d1,v4,s2,d2,tpp,tps,tsp,tss)
      implicit real*8 (a-h,o-z)
      common/rmode/love
      real*8 k1,k2,k3,k4,v1,v2,v3,v4,s1,s2,d1,d2,rho1,rho2
      complex*16 cr,e1,e2,e1p,e2p,c1,c2,c3,c4,c5,c6,ap,bp,bs,t
      complex*16 tpp,tps,ts,tsp,tss,p
      v1=v3
      v2=v4
      if(love.ne.2) go to 11
      v1=s1
      v2=s2
   11 continue
      rho1=d1
      rho2=d2
      k4=rho2*s2**2/(rho1*s1**2)
      b1=.5/(1-k4)
      b2=.5*k4/(k4-1)
      k1=b1/s1**2
      k2=b2/s2**2
      k3=k1+k2
      e1=cr(p,v1)
      e2=cr(p,v2)
      e2p=cr(p,s2)
      e1p=cr(p,s1)
      if(love.eq.2) go to 10
      c1=(p**2)*(k3-p**2)**2
      c2=p**2*e1*e1p*e2p
      c3=(e1*e1p)*(k2-p**2)**2
      c4=e2p*(k1-p*p)**2
      c5=k1*k2*e1*e2p
      c6=k1*k2*e1p
      ap=c1+c3-c5
      bp=c2+c4-c6
      bs=ap+e2*bp
      t=2.*k1*e1*(e2p*(k1-p**2)-e1p*(k2-p**2))
      tpp=t/bs
      t=2.*k1*p*e1*(e1p*e2-(k3-p**2))
      tps=t/bs
      ts=2.*k1*p*e1p*((k3-p**2)-e1*e2p)
      tsp=ts/bs
      ts=-2.*k1*e1p*(e1*(k2-p**2)-e2*(k1-p*p))
      tss=ts/bs
      return
 10   r1=d1*s1**2
      r2=d2*s2**2
      tss=(2.*r1*e1p)/ (r1*e1p+r2*e2p)
      tpp=0.
      tps=0.
      tsp=0.
      return
      end

c Geometric spreading factor given by 1.0/(sum of ray segments)
c
      subroutine geom_terms(hs,p0,itype,rp,qb)
      include 'params.h'
      implicit real*8 (a-h,o-z)
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rho(nlaymax),qp(nlaymax),qs(nlaymax)
      common/rays/nh(1,nlaymax),nm(1,nlaymax),ndeg(nlaymax),nd(nlaymax)
      real*8 dpt,th,vp,vs,rho,p0,hs,rp
      real*4 qp,qs,qb

      dep = 0.0
      do 1876 j=2,nh(1,1)-1
         dep = dep + th(j)
1876  continue

cXXX old way, hardwired to direct and 1 down-going Moho
cXXX      if(itype.eq.1) th1 = hs - dep
cXXX      if(itype.eq.2) th1 = dep + th(nh(1,1)) - hs

cXXX new way, hardwired to direct and 1 down-going Moho
      if(mod(itype,2).eq.1) th1 = hs - dep
      if(mod(itype,2).eq.0) th1 = dep + th(nh(1,1)) - hs

      sini = p0*vs(nh(1,1))
      if(sini.ge.1.0) sini = 0.999999
      denom = 1.0/(sqrt(1.0 - sini*sini))

      ri = th1*denom
      ti = ri/vs(nh(1,1))

      rsum = ri
      qb = ti/qs(nh(1,1))
      do 965 j=2,nd(1)

         sini = p0*vs(nh(1,j))
         if(sini.ge.1.0) sini = 0.999999
         denom = 1.0/(sqrt(1.0 - sini*sini))

	 ri = th(nh(1,j))*denom
	 ti = ri/vs(nh(1,j))

	 rsum = rsum + ri
	 qb = qb + ti/qs(nh(1,j))

965   continue

      if(rsum.eq.0.0) rsum = 0.001
      rp = rsum

      return
      end

      subroutine grandvel(j0,nlskip,vpsig,vssig,rhosig,qssig,icf)
      include 'params.h'
      common/vmod_in/dpt0(nlaymax),th0(nlaymax),vp0(nlaymax),vs0(nlaymax),rh0(nlaymax),qp0(nlaymax),qs0(nlaymax),gr(3000)
      common/vmod/dpt(nlaymax),th(nlaymax),vp(nlaymax),vs(nlaymax),rh(nlaymax),qp(nlaymax),qs(nlaymax)
      real*8 dpt,th,vp,vs,rh,vp0,vs0,rh0
      real*4 qp,qs,qs0,gr,strunc

C truncate at +/-strunc
      strunc = 3.0

c icf=correlation flag
c      if icf=1, all vel, den & q variations are exactly correlated
c      if icf=0, vel, den & q variations are completely un-correlated

      np6=(j0)*6
      call normal_random_number(np6,gr)

      do 6753 k=1,nlskip

         dpt(k) = dpt0(k)
         th(k) = th0(k)
         vp(k) = vp0(k)
         vs(k) = vs0(k)
         rh(k) = rh0(k)
         qp(k) = qp0(k)
         qs(k) = qs0(k)

6753  continue

      kt = (1-icf)*5*j0
      ks = (1-icf)*j0
      do 6754 k=nlskip+1,j0

         rr = gr(k+kt)
         if(rr.gt.strunc)rr=strunc
         if(rr.lt.-strunc)rr=-strunc
         th(k) = th0(k)*exp(rr*vpsig)

         rr = gr(k+ks)
         if(rr.gt.strunc)rr=strunc
         if(rr.lt.-strunc)rr=-strunc
         vp(k) = vp0(k)*exp(rr*vpsig)

         rr = gr(k+2*ks)
         if(rr.gt.strunc)rr=strunc
         if(rr.lt.-strunc)rr=-strunc
         vs(k) = vs0(k)*exp(rr*vssig)

         rr = gr(k+3*ks)
         if(rr.gt.strunc)rr=strunc
         if(rr.lt.-strunc)rr=-strunc
         rh(k) = rh0(k)*exp(rr*rhosig)

         rr = gr(k+4*ks)
         if(rr.gt.strunc)rr=strunc
         if(rr.lt.-strunc)rr=-strunc
         qp(k) = qp0(k)*exp(rr*qssig)

         rr = gr(k+4*ks)
         if(rr.gt.strunc)rr=strunc
         if(rr.lt.-strunc)rr=-strunc
         qs(k) = qs0(k)*exp(rr*qssig)

         if(k.gt.1) dpt(k) = th(k) + dpt(k-1)

6754  continue

      return
      end

      subroutine famprand(np2,cw,fgr,sg1,sg2)
      dimension fgr(1)
      COMPLEX*8 cw(1) 

      np = np2/2

      fac1 = exp(fgr(1)*sg1)

      do 2180 i=1,np+1

	 fac2 = exp(fgr(1+i)*sg2)
         cw(i) = cw(i)*fac1*fac2

2180  continue

      do 2181 i=1,np-1
         cw(np2-i+1) = conjg(cw(i+1))
2181  continue

      return
      end



	subroutine init_random_seed(irand)

	integer :: i,n,irand
	integer, dimension(:), allocatable :: seed

	call random_seed(size=n)
	allocate(seed(n))

c	print*,'n= ',n
	do i=1,n
	   seed(i)=irand
	   irand=irand+1
	enddo

	call random_seed(put=seed)

	deallocate(seed)

	return
	end

       
       	subroutine normal_random_number(nr,acc)
	integer nr,j,n
	dimension acc(nr)
	real*4 x1,x2,rand_numb

	j=1
	DO 8 n=1,nr

	   goto (1,2),j
1          x1=rand_numb(0)

964        continue
           if(x1.eq.0.0) then
              x1=rand_numb(0)
              go to 964
           endif

           x2=rand_numb(0)
963        continue
           if(x2.eq.0.0) then
              x2=rand_numb(0)
              go to 963
           endif

	   x2=6.2831853*x2
	   x1=-alog(x1)
	   x1=sqrt(x1+x1)
	   w=x1*cos(x2)
	   j=2
	   go to 3

2	   w=x1*sin(x2)
	   j=1

3	   continue
	   acc(n)=w

8	   continue

	   s=0.0
	   do 10 i=1,nr
	      s=s+acc(i)*acc(i)
10	   continue

	   s=sqrt(nr/s)
	   do 11 i=1,nr
	      acc(i)=acc(i)*s
11	   continue

	return
	end

	function rand_numb(iflag)
	integer iflag
	real*4 rand_numb

	if(iflag.gt.0) call init_random_seed(iflag)
	call random_number(rand_numb)

	return
	end
