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

c     Main

      include 'params.h'

      CHARACTER*256 asite,slip_model,outname,outname1,outname2,outname3,velfile
      CHARACTER*256 velname,velname1
      CHARACTER*256 dummy
      integer head_lines, j
      character*64 cap,cname(3) 

      dimension ic(2),c(4),irtype(100)
      integer dosetup
      real*4 rvfac,shal_rvfac,deep_rvfac,Czero,Calpha
      real*4 vpsig,vshsig,rhosig,qssig,spar_exp
      integer nlskip,icflag
      real*8 depth,thic,vp,vsh,rho,vp0,vsh0,rho0,vsmoho
      real*4 stime,rpath,p0,qbar,qfexp
      real*4 qp,qs,qs0,grand

      integer isite_amp,nbu,iftt,nrtyp
      integer irand,nsite,ipdur_model,ispar_adjust,ndata
      real*4 targ_mag,fault_area,vr
      real*4 flol,fhil,duration,fmx,akapp,fasig1,fasig2,rvsig1
      real*4 sm,stress_average,rvfmax,bigC,orig_stress_average

      DIMENSION RNA(mmv),RNB(mmv),RDNA(mm),DS(3,mmv)                       
      dimension sddp(lv,nq,np),rist(lv,nq,np),rupt(lv,nq,np),sddp_orig(lv,nq,np)
CFFF Add variable astopq-> This is the along strike distance (km) on
CFFF the top edge of the fault of the reference point (elonq,elatq).
CFFF This replaces ntopq in even_dist1
      dimension strq(lv),dipq(lv),rakeq(lv),elonq(lv),elatq(lv)
      dimension astop(lv),dtop(lv),shyp(lv),dhyp(lv),nx(lv),nw(lv),dx(lv),dw(lv)

      
      COMMON/RANDO/ifu 

      common /consts/ irtype,isite_amp,nbu,iftt,flol,
     + fhil,nsite,fmx,akapp,Czero,Calpha,shal_rvfac,deep_rvfac,
     + sm,vr,vsmoho,vpsig,vshsig,rhosig,qssig,icflag,fasig1,fasig2,rvsig1,
     + ipdur_model,ispar_adjust,fault_area,nrtyp,rvfmax,dosetup,bigC,ndata,
     + spar_exp
      save /consts/

      pu =3.1415926/180 
      pai=3.1415926 
      IRD=2
      IMDL=2 
      krc = 4

      bigC = -1.0

      read(5,*)orig_stress_average
      print*,'stress_average= ',orig_stress_average

      read(5,'(a256)') asite 
      read(5,'(a256)') outname                          
      read(5,*) nrtyp,(irtype(i),i=1,nrtyp)
      read(5,*) isite_amp
      read(5,*) nbu,iftt,flol,fhil   
      read(5,*) irand                                       
      read(5,*) nsite
      read(5,*) duration,dt,fmx,akapp,qfexp

      read(5,*) rvfac,shal_rvfac,deep_rvfac,Czero,Calpha

      fcfac = fcfac_default
      fcfac = 0.0

      read(5,*) sm,vr

      read(5,'(a256)') slip_model  
      open(10,file=slip_model)

      read(10,*) nevnt

      if(nevnt.gt.lv) then
         write(6,*) 'nevnt=',nevnt,' and is greater than lv=',lv,'.Change LV in params.h'
         stop
      endif

      zhyp_max = 0.0;
      nstot = 0
      farea_in = 0.0
      do 133 iv=1,nevnt
         read(10,*) elonq(iv),elatq(iv),nx(iv),nw(iv),dx(iv),dw(iv)
         read(10,*) strq(iv),dipq(iv),rakeq(iv),dtop(iv),shyp(iv),dhyp(iv)

	print*,'nx,nw= ',nx(iv),nw(iv)

         do j=1,nw(iv)
            read(10,*) (sddp_orig(iv,i,j),i=1,nx(iv))
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

CCC New way read velocity model from file

      read(5,'(a256)') velfile

C RWG 20140313
c  Option to specify Vs of layer just beneath Moho in case input velocity model
c  extends into upper mantle.  Downgoing rays are reflected from top of last layer
c  in model, so all layers beneath the layer with Vs>vsmoho are truncated.
      vsmoho = -1.0
      read(5,*) vsmoho
      if(vsmoho.le.0.0) vsmoho=999.9

      nlskip = -99
      vpsig = -1.0
      vshsig = -1.0
      rhosig = -1.0
      qssig = -1.0

      read(5,*) nlskip,vpsig,vshsig,rhosig,qssig,icflag
      read(5,'(a256)') velname                          

      rvsig1 = -1.0
      rvfmax = 1.4

      fasig1 = -1.0
      fasig2 = -1.0

      read(5,*) fasig1,fasig2,rvsig1

CCC RWG modified 2014-09-19
      read(5,*) ipdur_model

CCC RWG modified 2018-08-30
c adjust stress parameter based on magnitude and area relative to
c standard scaling relation
c     stress parameter adjustment options:
c
c     ispar_adjust = 0, no adjustment
c     ispar_adjust = 1, adjust to Leonard (2010) for active tectonic
c     ispar_adjust = 2, adjust to Leonard (2010) for stable continent

      read(5,*) ispar_adjust,targ_mag,fault_area,spar_exp

c  Stations location
           
      open(1,file=asite)

      loc1= index(outname,' ')-1

c      ifu=irand
c      call init_random_seed(irand)
c      call RANU2(NR,RNA)
c      call RANU2(NR,RNB) 
      dosetup=1

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

      sddp = sddp_orig
      stress_average = orig_stress_average 
      call hb_high(cap, stlon, stlat, velfile, outname, duration, dt, ds, nevnt,
     + elonq, elatq, targ_mag, nlskip, nx, nw, dx, dw, strq, dipq, rakeq, dtop, shyp,
     + dhyp, sddp, rist, rupt, qfexp, velname, stress_average, rvfac, 
     + irand, msite)

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

555   continue

100     format(a8,2x,a4,1x,a56)
102     format(1x,i6,1x,e11.5,2i3,4(1x,f10.4))   

c     end main
1       end
