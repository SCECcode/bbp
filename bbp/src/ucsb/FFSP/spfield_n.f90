! spfield.f90
! The subroutines in this file are for generating the spatial
! distribution of source parameters on main-fault
!
! Written by Pengcheng Liu
! Copyright (c) 2005 by Pengcheng Liu
!
! Modifications
!      1. frequency band used to constrain the model
!      2. using dcf to constrain the spectral level
!======================================================================
subroutine mainfault(dip,freq_min,freq_max,rv_avg,ratio_rise,nb_taper_TRBL)
! Set the initial parameters for generating the source paramters
! flx,flw-- fault length (in strike) and width (in dip)
! the spectrum between freq_min and freq_max will be used to constrain the model
!
 use sp_sub_f
 use fdtim_2d
 use time_freq
 implicit NONE
 integer,dimension(4), intent(IN):: nb_taper_TRBL
 real, intent(IN):: dip,freq_max,freq_min,rv_avg,ratio_rise
 integer:: i,j,k,lyr,nb,ll,nx1,ny1
 real:: rdip,yij,hk,sumh,sv,sm,rat,ttmp,vs_avg
!
 cft1=ratio_rise
 nsum=nsubx*nsuby
 dx=flx/nsubx
 dy=fwy/nsuby
!
! why 1e-15?  dx dy in km -->10^-6, beta and density in km and ton ->10^-9
!
 smoment= Moment_o*1.e-15
!
! grid_2d <0.2 km, in contrast, 1 km in GP, 0.33km in Frankel
!
 grid_2d=amin1(0.2,  amin1(dx,dy))
! nx_2d: >= nsubx
! nz_2d: >= nsuby

! redefine the hypocenter location

 nx_2d=int((nsubx-1)*dx/grid_2d)+2
 nz_2d=int((nsuby-1)*dy/grid_2d)+2
 ncxp=int(x_hypc/dx) + 1
 ncyp=int(y_hypc/dy) + 1
 cxp=(ncxp-0.5)*dx
 cyp=(ncyp-0.5)*dy
 x_hypc=cxp
 y_hypc=cyp
!
 if(allocated(slip)) deallocate(slip,rstm,rptm,pktm,rpvel,beta,amu, &
                     taper,rtx,rtz,rake_prt,dip_prt,amz_prt)
 allocate(slip(nsum),rstm(nsum),rptm(nsum),pktm(nsum),lrtp(nsum),rpvel(nsum), &
          beta(nsum),amu(nsum),taper(nsum),rtx(nsum),rtz(nsum))
 allocate(rake_prt(nsum),dip_prt(nsum),amz_prt(nsum))
!
 rdip=dip*4.*atan(1.0)/180.0
 do j=1,nsuby
   yij=(j-0.5)*dy-cyp
   hk=yij*sin(rdip)+depth_hypc
   thk(layer)=hk+0.5
   lyr=0
   sumh=0
   do while(hk > sumh)
     lyr=lyr+1
     sumh=sumh+thk(lyr)
   enddo
!
   sv=vvs(lyr)
   sm=roh(lyr)*sv*sv
   do i=1,nsubx
     k=(i-1)*nsuby+j
     beta(k)=sv
     amu(k)=sm
     rtx(k)=(i-0.5)*dx
     rtz(k)=(j-0.5)*dy
   enddo
 enddo
! smoothing
 do j=3,nsuby-2
   ll=j+nsuby
   amu(j)=0.2*sum(amu(ll-2:ll+2))
 enddo
!
 vs_avg=0.0
 do j=1,ncyp
   sv=0.0
   do k=j, ncyp
     sv=sv+1./beta(k)
   enddo
   beta(j)=float(ncyp-j+1)/sv
   vs_avg=vs_avg+beta(j)
 enddo
 do j=nsuby, ncyp+1, -1
   sv=0.0
   do k=ncyp+1, j
     sv=sv+1./beta(k)
   enddo
   beta(j)=float(j-ncyp)/sv
   vs_avg=vs_avg+beta(j)
 enddo
 do j=1,nsuby
   do i=2,nsubx
     k=(i-1)*nsuby+j
     ll=k-nsuby
     beta(k)=beta(ll)
     amu(k) =amu(ll)
   enddo
 enddo
 vs_avg=vs_avg/float(nsuby)
!
!  Setting the parameters for computing random Slip-Amplitude field
!
 select case(id_cl)
 case(1)
!The correlation length vs. magnitude (Mai and Beroza, 2002,JGR)
   slip_clx =10.**(-2.5+Mw/2.)
   slip_cly =10.**(-1.5+Mw/3.)
 case(2)
! The correlation length vs. magnitude (Grave,SCEC 2005)
   slip_clx =10.**(-2.0+Mw/2.)
   slip_cly =10.**(-2.0+Mw/2.)
 case(3)
!The correlation length vs. magnitude (Somerville, 1999,SRL)
   slip_clx =10.**(-1.72+Mw/2.)
   slip_cly =10.**(-1.93+Mw/2.)
 case(4)
! The correlation length vs. magnitude (Graves & Pitaka, 2010)
    slip_clx =10.**(-1.7+Mw/2.)
    slip_cly =10.**(-0.7+Mw/3.)
 end select
!
!  correlation length less than 0.7* fault width/length
!
 slip_cly=amin1(slip_cly, 0.7*fwy)
 slip_clx=amin1(slip_clx, 0.7*flx)
! correlation length along strike shall be less than 2*c_width
 slip_clx=amin1(slip_clx, amax1(2.0*slip_cly, flx/4.0))
 if (flx < fwy) then
   ttmp=slip_cly
   slip_cly=slip_clx
   slip_clx=ttmp
endif
!!! ? check the value for M>7.3 event.
if(Mw >7.3) print *, Mw, slip_clx, slip_cly
 sp1 = 2.0
 sq1 = 5.0
 slip_pw  = 2.0
! hfrd=smoment*fc_main*fc_main
 hfrd=smoment*fc_main_1*fc_main_2
!
! Setting the parameters for computing random Risetime filed
!
 rs_clx = amin1(slip_clx, 10.0)
 rs_cly = amin1(slip_cly, 10.0)
 rs_max=6.0
!
  rs_max= amax1(slip_clx,slip_cly)/2.5 ! asperity duration for 2.5 km/s rupture velocity
!
!   Somerville et al (1999), Chen Ji, 2020
!   Note that the final risetime=rstm+pktm, so rstm = 0.8 risetime is assigned
!
 rstm_mean = 0.8*10**(0.5*Mw-3.34)
 pktm_mean = 0.2 *rstm_mean
!
 rp1= 3.0
 rq1= 2.0
 rs_pw  = 2.0
!
rs_pw = 0.54*slip_pw+0.675 ! using SAL12 relationship, Chen Ji, 2020
!
! Setting the parameters for computing random rupture-velocity filed
!
 rv_clx = amin1(slip_clx, 10.0)
 rv_cly = amin1(slip_cly, 10.0)
!
 rv_pw = 2.0
!
 rv_pw = 0.23*slip_pw + 0.370   ! using SAL12 relationship, Chen Ji, 2020
!
! Setting the parameters for computing random peak time filter_field
! Chen Ji, 2020
!
pk_pw = 0.29*slip_pw + 0.415 ! using SAL12 relationship, Chen Ji, 2020
pk_max= 0.2*rs_max
!
! Setting the parameters for computing rupture velocity
!
 vp1=1.0
 vq1=1.0
 vmin_ratio = 0.6!0.5
 vmax_ratio = 0.85!0.9
 vs_avg=vs_avg*(vq1*vmin_ratio+vp1*vmax_ratio)/(vp1+vq1)
 vs_avg=rv_avg/vs_avg
 vmin_ratio = vmin_ratio*vs_avg
 vmax_ratio = vmax_ratio*vs_avg
 write(*,*) 'ratio= ',vs_avg
!
! Setting the cross-correlation coefficients
! between slip-amplitude and rise-time, and
! between slip-amplitude and rupture-velocity. (default: no correlation?)
!
 cr_sl_rs = 0.75 !0.62
 cr_sl_vc = 0.0  !0.3
! cr_sr_vc = 0.4 ! SAL10 relationship, Chen Ji, 2020
 cr_pk_vc = -0.31
!
! Tapering the slip-amplitudes nearby the edage of fault
 taper= dx*dy
! Top Side
 nb=nb_taper_TRBL(1)
 do j=1,nb
   rat=(j-0.5)/nb
   do i=1,nsubx
     k=(i-1)*nsuby+j
     taper(k)=taper(k)*rat
   enddo
 enddo
! Bottom Side
 nb=nb_taper_TRBL(3)
 do j=1,nb
   rat=(j-0.5)/nb
   do i=1,nsubx
     k=(i-1)*nsuby+(nsuby-j+1)
     taper(k)=taper(k)*rat
   enddo
 enddo
! Right
 nb=nb_taper_TRBL(2)
 do i=1,nb
   rat=(i-0.5)/nb
   do j=1,nsuby
     k=(nsubx-i)*nsuby+j
     taper(k)=taper(k)*rat
   enddo
 enddo
! Left
 nb=nb_taper_TRBL(4)
 do i=1,nb
   rat=(i-0.5)/nb
   do j=1,nsuby
     k=(i-1)*nsuby+j
     taper(k)=taper(k)*rat
   enddo
 enddo
!
! i, j, ttmp: maximum on-fault hypocenter distance
 i=max0(ncxp, nsubx-ncxp)
 j=max0(ncyp, nsuby-ncyp)
 ttmp=sqrt((i*dx)**2+(j*dy)**2)
 vs_avg=rv_avg*2.0*vmin_ratio/(vmin_ratio+vmax_ratio)
 dt=0.001
 i=int(ttmp/vs_avg/dt+1.0)*10
 ntime=1
 do while(ntime<i)
    ntime=ntime*2
enddo
! ntime: length of rise time function
 df=1./(ntime*dt)
 nphf=ntime/2
 if(allocated(freq)) deallocate(freq)
 allocate(freq(nphf))
!
 lnpt=int(log(real(ntime))/log(2.0)+0.1)
 !nfe1=int(0.5*freq_max/df)+1
 !nfe2=int(0.9*freq_max/df)
 nfe1=int(1.1*freq_min/df)+1
 nfe2=int(0.9*freq_max/df)+1
 nbb=int(log(real(nfe1))/log(2.))
 nee=int(log(real(nfe2))/log(2.))+1
 ! In this study, we use 95% seismic moment as source duration criterion
 ! A empirical correction of 0.78 is added (1-sqrt(0.05)~0.78),
 t_precent=0.95
 td_target=0.78/fc_main_1/3.14

!
! Using the acceleration spectrum as the
 io_dva=3
!
 write(*,*) "Finishing the setup of global variable"
 write(*,*) "Magnitude = ",Mw
 write(*,*) "expected double corner frequency = ", fc_main_1,fc_main_2
 write(*,*) "expected Td_95= ",td_target
 write(*,*) "expected average rise time =", rstm_mean
 write(*,*) "ntime = ", ntime, " lnpt = ",lnpt
 write(*,*) "nphf = ", nphf
 write(*,*) "nbb = ", nbb, " nee = ",nee
 write(*,*) "io_dva = ",io_dva
!
end subroutine mainfault
!
!======================================================================
! procedure
!     1. slip distribution: A k^-2 random field that has given comulative moment
!     2. rise time:   A k^-1.75 random field with given rstm_mean
!     A Loop of peaktime_mean
!        peak time:  A k^-1 random field with peaktime_mean
!        estimate peak-slip rate
!        rupture velocity field
!        generate the source spectrum
!        change the peaktime_mean
!     end-loop
!     Note that here we assume that the peaktime << risetime
!     Chen Ji, 2020
!
subroutine random_field(idum1,idum2,idum3,ave_tr,ave_tp,ave_vr,err_spectra)
 use sp_sub_f
 use fdtim_2d
 implicit NONE
 integer, intent(INOUT):: idum1,idum2,idum3
 integer, parameter:: itrmax=10
 integer idum4
 logical:: lloop,hloop
 integer:: i,j,iset,niter,is_stress
 real::ave_tr,ave_tp,ave_vr,err_spectra
 real:: sum1,cmti,correct,cr0,cr1,cr2,csqt,cv0,cv1,cv2, &
        vrup,er_syn,er_target,cr_sr,cr_sv,gasdev

 real,dimension(nsum):: gasv1,gasv2,gasv3,gasv4,stdrop
!
! new parameters using GP approach for local rupture peturbation time (lrpt)
!
 real:: slip_mean, slip_max, slip_white,dt_scale
! new parameters using Gusev approach
 real:: rup_time_max,lrtp_clx,lrtp_cly,lrtp_pw,coef_lrtp
 real:: rayl_sigma
 real:: precent,td,tb,te
!
! Generating 4 Normal distributed white noise
!
 iset=0
 do i=1,nsum
   gasv1(i)=gasdev(idum1,iset)
 enddo
 iset=0
 do i=1,nsum
   gasv2(i)=gasdev(idum2,iset)
 enddo
 iset=0
 do i=1,nsum
   gasv3(i)=gasdev(idum3,iset)
 enddo
 iset=0
 idum4=idum1+nsum
 do i=1,nsum
   gasv4(i)=gasdev(idum4,iset)
 enddo
! 1. slip distribution
!
 slip=gasv1
 call filter_field(slip,slip_clx,slip_cly,slip_pw,dx,dy,nsubx,nsuby)
! mapping to the desired distributions using NORTA technique
! by keeping cumulative probability same, Liu et al. (2006)
! slip amplitude following Cauchy distribution
 call random_slip_am_C(slip,nsum)
!
 stdrop=slip*taper    ! unit is meter
!
 slip=slip*amu*taper
 correct=smoment/sum(slip(1:nsum))
 slip=slip*correct   ! Note: now the unit of slip is not meter.
!
! Stress
! stdrop=slip
!
!call stress_drop(stdrop,dx,dy,nsubx,nsuby)
!call coherent(cr_sr,nsum,slip,stdrop)
!write(*,*) 'correlation between the fault slip and static stress drop = ',cr_sr
!
! 2. rise time distribution
!
 cr0=cr_sl_rs
 cr1=0.0
 cr2=1.0
 if(cr_sl_rs < 0.0) then
   cr1=-1.0
   cr2=0.0
 endif
 niter=0
 lloop=.true.
 do while(lloop)
   niter=niter+1
   csqt=sqrt(1.-cr0*cr0)
   do i=1,nsum
     rstm(i)=abs(cr0)*gasv1(i)+csqt*gasv2(i)
   enddo
   call filter_field(rstm,rs_clx,rs_cly,rs_pw,dx,dy,nsubx,nsuby)
   call random_rise_tm(rstm,nsum,rs_max,rp1,rq1)
   correct=rstm_mean/(sum(rstm(1:nsum))/nsum)
   rstm=rstm*correct
   call coherent(cr_sr,nsum,slip,rstm)
   write(*,*) 'cross rise, slip-amp',cr_sr,cr0,cr_sl_rs
!
   if(abs(cr_sr-cr_sl_rs) < 0.01 .or. niter > itrmax) lloop=.false.
   if(cr_sr > cr_sl_rs) then
     cr2=cr0
   else
     cr1=cr0
   endif
   cr0=0.5*(cr2+cr1)
 enddo
 call coherent(cr_sr,nsum,slip,rstm)
 write(*,*) 'correlation between slip and rise time = ',cr_sr
 call coherent(cr_sr,nsum,rstm,stdrop)
 write(*,*) 'correlation between stress drop and rise time = ',cr_sr
!
! 3. distribution of peak time
! Using SAJ13, peak time has k^~-1 spatial distribution
!
!pktm_mean=0.2*rstm_mean
write(*,*)"creating k^-1 peak time distribution with mean of ", pktm_mean,pk_pw
! SAL10: negligible correlation between fault slip and peak time
! Using a new gaussian vector
pktm=gasv4
call filter_field(pktm,rs_clx,rs_cly,pk_pw,dx,dy,nsubx,nsuby)
!call random_rise_tm(pktm,nsum,pk_max,rp1,rq1)
call random_slip_am_C(pktm,nsum)
! normalization
correct=pktm_mean/(sum(pktm,1)/nsum)
pktm=pktm*correct
call coherent(cr_sv,nsum,pktm,slip)
!
write(*,*)"correct = ",correct
write(*,*)"correlation between pktm and slip= ", cr_sv
call coherent(cr_sv,nsum,pktm,rstm)
write(*,*)"correlation between pktm and rise time= ", cr_sv
!
! rupture initiation time
!
ave_vr=0.7
niter=0
lloop=.true.
do while(lloop)
  niter=niter+1
  do i=1,nsum
    vrup=ave_vr*beta(i)
    rpvel(i)=vrup
    rptm(i)=1.0/vrup
  enddo
  call rup_tim_2d(nsum,nsuby,nsubx,dx,dy,rptm)
  rup_time_max=maxval(rptm)
  write(*,*)"maximum rupture time =",rup_time_max
!  id_lrtp=1
!
! rupture time perturbation
!
! Hisada&Gusev approach 1. GP approach = 2
! Using GP approach, will lead to about 10% increase in rise time.
! Chen Ji, 2020
!
  select case(id_lrtp)
  case(0)
  ! spatially uncorrelated but
    coef_lrtp=0.5
    lrtp=gasv1
    rayl_sigma=rstm_mean/2.0
    call random_rayleigh(lrtp,nsum,rayl_sigma)
    call coherent(cr_sv,nsum,lrtp,slip)
    write(*,*)" case 0 correlation between lrtp and slip= ", cr_sv
    do i=1,nsum
      lrtp(i)=lrtp(i)*((rptm(i)/rup_time_max)**coef_lrtp)
    enddo
    correct=rstm_mean/(sum(lrtp)/nsum)
    lrtp=lrtp*correct
  case(1)
  ! Guesv didn't mention the correlation length, here we use a half of
  ! slip pulse as correlation length, following Bernard (1994). We also
  ! assume the rupture velocity is 2.5 km/sec.
  !
    lrtp_clx=rstm_mean*2.5
    lrtp_cly=lrtp_clx
    lrtp_pw=2.0
    coef_lrtp=0.25
    rayl_sigma=rstm_mean/2.0
  !
    lrtp=gasv3
    call filter_field(lrtp,lrtp_clx,lrtp_cly,lrtp_pw,dx,dy,nsubx,nsuby)
    call random_rayleigh(lrtp,nsum,rayl_sigma)
    write(*,*)"maximum and minimum ",maxval(lrtp),minval(lrtp),rayl_sigma
    call coherent(cr_sv,nsum,lrtp,slip)
    write(*,*)"Gusev correlation between lrtp and slip= ", cr_sv
    write(*,*)"maximum lrtp = ",maxval(lrtp)
    do i=1,nsum
      lrtp(i)=lrtp(i)*((rptm(i)/rup_time_max)**coef_lrtp)
    enddo
    correct=0.5*rstm_mean/(sum(lrtp)/nsum)
    lrtp=-lrtp*correct
    write(*,*)"maximum lrtp = ",maxval(lrtp), correct, rstm_mean
  case(2)
  !
  ! GP approach
  ! GP's approach is deterministic
  !
    slip_mean=sum(slip)/nsum
    slip_max=maxval(slip)
    slip_white=0.05*slip_mean
    dt_scale=1.8e-9*((Moment_o*1.0e7)**(1.0/3.0))
    dt_scale=log(slip_max/slip_mean)*dt_scale*0.5
    coef_lrtp=0.25
    do i=1,nsum
      lrtp(i)=-dt_scale*(log(slip(i)/slip_mean))*((rptm(i)/rup_time_max)**coef_lrtp)
    end do
    call coherent(cr_sv,nsum,lrtp,slip)
    write(*,*)"GP correlation between lrtp and slip= ", cr_sv
  end select

  do i=1,nsum
    rptm(i)=rptm(i)+lrtp(i)
    if(rptm(i).lt.0.0)rptm(i)=0.0
  enddo

  precent=t_precent
  call stf_duration(precent,td,tb,te)
  write(*,*)"duration =",precent,td,tb,te
  if(abs(td/td_target-1.0)<0.025 .or. niter > itrmax) lloop=.false.
  write(*,*)"Before average rupture velocity = ", td, td_target, ave_vr
  ave_vr=ave_vr*td/td_target
  write(*,*)"After average rupture velocity = ", td, td_target, ave_vr
enddo

! As high frequency radiation is more sensitive to rise time,
call misfit_ratio(smoment,fc_main_1,fc_main_2,correct,err_spectra)
write(*,*) 'the misfit in ratio ', correct, (sum(rstm,1)/nsum),err_spectra
if(correct.lt.0.8)then
  do i=1,nsum
    rstm(i)=rstm(i)*((correct/0.8)**1.4)
    pktm(i)=pktm(i)*((correct/0.8)**1.4)
  enddo
  write(*,*) 'the misfit in ratio step-1', correct, (sum(rstm,1)/nsum)
endif
if(correct.gt.1.25)then
  do i=1,nsum
    rstm(i)=rstm(i)*((correct/1.25)**1.4)
    pktm(i)=pktm(i)*((correct/1.25)**1.4)
  enddo
  write(*,*) 'the misfit in ratio step-1', correct, (sum(rstm,1)/nsum)
endif
niter=0
hloop=.true.
do while(hloop)
  call misfit_ratio(smoment,fc_main_1,fc_main_2,correct,err_spectra)
  write(*,*) 'the misfit in ratio step-2', correct, (sum(pktm,1)/nsum),err_spectra
  if(abs(correct-1.) < 0.03 .or. niter > itrmax) then
    hloop=.false.
  else
    pktm=pktm*(correct**2.5)
  endif
  niter=niter+1
enddo

! check the quality of solution
write(*,*)"re-evaluate Quality of selected source model"
call misfit_ratio(smoment,fc_main_1,fc_main_2,correct,err_spectra)
write(*,*)"Misfit_ratio = ", correct,err_spectra
call coherent(cr_sv,nsum,rstm,pktm)
write(*,*) 'cross-correlation between risetime(y) and peak time =',cr_sv
do i=1,nsum
  rstm(i)=rstm(i)+pktm(i)
end do
ave_tr=sum(rstm)/nsum
ave_tp=sum(pktm)/nsum
call coherent(cr_sv,nsum,rstm,slip)
write(*,*) 'cross-correlation between risetime and slip =',cr_sv,'target =',cr_sl_rs
call coherent(cr_sv,nsum,rstm,pktm)
write(*,*) 'cross-correlation between risetime and peak time =',cr_sv
do i=1,nsum
  pktm(i)=pktm(i)/rstm(i)
end do
call stf_synth_output(smoment,fc_main_1,fc_main_2)
!
! scale the unit of seismic moment to N-m
  slip=slip*1.0e+15
!
! Perturbate the rake dip, and azimuth angle
!
  iset=0
  do i=1,nsum
    rake_prt(i)=gasdev(idum2,iset)
  enddo
  call filter_field(rake_prt,slip_clx,slip_cly,2.0,dx,dy,nsubx,nsuby)
  call random_rake(rake_prt,nsum)

  iset=0
  do i=1,nsum
    dip_prt(i)=gasdev(idum2,iset)
  enddo
  call filter_field(dip_prt,slip_clx,slip_cly,2.0,dx,dy,nsubx,nsuby)
  call random_rake(dip_prt,nsum)

  iset=0
  do i=1,nsum
    amz_prt(i)=gasdev(idum2,iset)
  enddo
  call filter_field(amz_prt,slip_clx,slip_cly,1.5,dx,dy,nsubx,nsuby)
  call random_rake(amz_prt,nsum)

end subroutine random_field
!
!======================================================================
! origin name: random_rup_tim
!
subroutine rup_tim_2d(nsum,nsuby,nsubx,dx,dy,rptm)
 use fdtim_2d
 implicit NONE
 integer, intent(IN):: nsum,nsuby,nsubx
 real:: dx,dy,rptm(nsum)
 integer:: i,j,k,jx,jz,jxz,jxzx
 real:: rx,rz,rz1
 real, dimension(nx_2d, nz_2d):: rtm_2d,slow
!
! as input, rptm is local slowness of rupture velocity
! Chen Ji, 2020
! do k=1,nsum
!   rptm(k)=sqrt((rtx(k)-cxp)**2+(rtz(k)-cyp)**2)*rptm(k)
! enddo
! return
!
 do k=1,nz_2d
   rz=float(k-1)*grid_2d/dy
   jz=min0(int(rz)+1,nsuby-1)
   rz=amax1(jz-rz, 0.0)
   rz1=1.-rz
   do i=1,nx_2d
     rx=float(i-1)*grid_2d/dx
     jx=min0(int(rx)+1, nsubx-1)
     rx=amax1(jx-rx, 0.0)
     jxz=(jx-1)*nsuby+jz
     jxzx=jxz+nsuby
     slow(i,k)=rx*(rz*rptm(jxz) +rz1*rptm(jxz+1))+ &
          (1.-rx)*(rz*rptm(jxzx)+rz1*rptm(jxzx+1))
   enddo
 enddo
 call XZPTIM(rtm_2d,slow,nx_2d,nx_2d,nz_2d,grid_2d, &
             cxp,cyp,rptm,rtx,rtz,nsum)

end subroutine rup_tim_2d
!
!
!======================================================================
!
subroutine random_slip_am_C(xp,n)
!  a=a0*averg
!  b=d1*averg
!  xmax=a1*averg
!  af0=atan(-a/b)
!  af1=atan((xmax-a)/b)
!  pdf=1./((af1-af0)*b)/(1.+((x-a)/b)**2)
!
 implicit NONE
 real, parameter:: a0 = 0.30, a1 = 3.5, averg=1.0
 integer:: n,i
 real, dimension(n):: xp
 real:: af0,af1,d1,daf,r

 call gasprb(xp,n)
 call coef_slip_am(a0,a1,af0,af1,d1)
 daf=af1-af0
 do i=1,n
   r=tan(af0+xp(i)*daf)
   xp(i)=(a0+d1*r)*averg
 enddo

end subroutine random_slip_am_C
!
!======================================================================
!
subroutine coef_slip_am(a0,a1,af0,af1,d1)
 implicit NONE
 real, parameter:: err=1.e-5
 real:: a0,a1,af0,af1,d1,x0,x1,ff
 logical:: lloop
 x0=0.01
 x1=a1*4
 lloop=.true.
 do while(lloop)
   d1=0.5*(x0+x1)
   af0=-atan(a0/d1)
   af1= atan((a1-a0)/d1)
! Fitting Average
   ff=-1.0+a0+0.5*d1/(af1-af0)*alog(((a1-a0)**2+d1**2)/(a0**2+d1**2))
!
!   Fitting Varation
!   ff=-1.0+a1*d1/(af1-af0)+2.*a0*(a0+0.5*d1/(af1-af0)* &
!       alog(((a1-a0)**2+d1**2)/(a0**2+d1**2)))-(a0*a0+d1*d1)
!
   if(abs(ff) < err) lloop=.false.
   if(ff > 0.0) then
     x1=d1
   else
     x0=d1
   endif
 enddo

end subroutine coef_slip_am
!
!======================================================================
!
subroutine random_slip_am_B(xp,n,p,q)
!  beta distribution
!  pdf(x)=coef*(x-x0)**c1*(x1-x)**c2, x0 <= x <= x1 <=1
!  pdf(x0)=0.0
!  pdf(x1)=1.0
 implicit NONE
 integer, parameter:: nseg =10000
 integer n,i,j0,j1,jj
 real:: xp(n),p,q
 real:: pacul(nseg),dx,xi,x0,x1,c1,c2,pr
 x0=0.0
 x1=1.0
 call gasprb(xp,n)
 x0=0.0
 x1=1.0
 c1=p-1
 c2=q-1
 dx=(x1-x0)/(nseg-1.)
 pacul(1)=0.0
 do i=2,nseg
   xi=x0+(i-1.5)*dx
   pacul(i)=pacul(i-1)+(xi-x0)**c1*(x1-xi)**c2
 enddo
 do i=2,nseg
   pacul(i)=pacul(i)/pacul(nseg)
 enddo

 do i=1,n
   pr=xp(i)
   j0=1
   j1=nseg
   jj=ifix(nseg*pr-pr)+1
   do while(j1-j0 > 1)
     if(pr > pacul(jj) ) then
       j0=jj
     else
       j1=jj
     endif
     jj=(j0+j1)/2
   enddo
! linearly interpolating between j0 and j0+1
   xp(i)=x0+dx*(j0-1.0+(pr-pacul(j0))/(pacul(j0+1)-pacul(j0)))
 enddo

end subroutine random_slip_am_B
!
!==========================  Rise Time  =============================
!
subroutine random_rise_tm(xp,n,rs_max,p,q)
!  P(x0)=0.0
!  P(x1)=1.0
!  pdf(rs)=1./(x1-x0)**(p+q-1)/beta(p,q)*(x-x0)**(p-1)*(x1-x)**(q-1)
!
 implicit NONE
 integer, parameter:: nseg =10000
 integer:: n,i,j0,j1,jj
 real:: p,q,rs_max,xp(n)
 real:: pacul(nseg),dx,xi,x0,x1,c1,c2,pr,cff

 call gasprb(xp,n)
 x0=1.0
 x1=rs_max
 c1=p-1
 c2=q-1
 cff=1.-rs_max**p

 dx=(x1-x0)/(nseg-1.)
 pacul(1)=0.0
 do i=2,nseg
   xi=x0+(i-1.5)*dx
   pacul(i)=pacul(i-1)+(xi-x0)**c1*(x1-xi)**c2
 enddo
 do i=2,nseg
   pacul(i)=pacul(i)/pacul(nseg)
 enddo
!
 do i=1,n
!   pr=1.-xp(i)
   pr=xp(i)
   j0=1
   j1=nseg
   jj=ifix(nseg*pr-pr)+1
   do while(j1-j0 .gt.1)
     if(pr > pacul(jj) ) then
       j0=jj
     else
       j1=jj
     endif
     jj=(j0+j1)/2
   enddo
! linearly interpolating between j0 and j1
   xp(i)=x0+dx*(j0-1.0+(pr-pacul(j0))/(pacul(j0+1)-pacul(j0)))
 enddo
!
! do i=1,n
!   xp(i)=1./xp(i)
! enddo

end subroutine random_rise_tm
!
!========================== random_rayleigh ====================
! Note: for rayleigh distribution, mode is sigma, 95% of values are less than 2.5*sigma.
! Chen Ji, 2020
subroutine random_rayleigh(xp,n,sigma)
  implicit NONE
  integer:: n,i
  real:: sigma,xp(n),pr
  real:: irayl

  call gasprb(xp,n)
  do i=1,n
    pr=xp(i)
    xp(i)=irayl(pr,sigma)
  end do
end subroutine random_rayleigh
!
!==========================  Rake  =============================
!
subroutine random_rake(xp,n)
 implicit NONE
 integer, parameter:: nseg =10000
 integer:: n,i
 real:: xp(n), pacul(nseg)

 call gasprb(xp,n)
 do i=1,n
   xp(i)=2.*xp(i)-1.0
 enddo

end subroutine random_rake
!
!======================================================================
!
subroutine random_rup_vel(xp,n,x0,x1,p,q)
!  beta : local S-wave velocity
!  rup_vel=beta*x
!  vmin=x0*beta
!  vmax=x1*beta
!  pdf(x)=coef*(x-x0)**c1*(x1-x)**c2, x0 <= x <= x1 <=1
!  pdf(x0)=0.0
!  pdf(x1)=1.0
 implicit NONE
 integer, parameter:: nseg =10000
 integer:: n,i,j0,j1,jj
 real:: xp(n),x0,x1,p,q
 real:: pacul(nseg),c1,c2,dx,xi,pr

 call gasprb(xp,n)
 c1=p-1
 c2=q-1
 dx=(x1-x0)/(nseg-1.)
 pacul(1)=0.0
 do i=2,nseg
   xi=x0+(i-1.5)*dx
   pacul(i)=pacul(i-1)+(xi-x0)**c1*(x1-xi)**c2
 enddo
 do i=2,nseg
   pacul(i)=pacul(i)/pacul(nseg)
 enddo

 do i=1,n
   pr=xp(i)
   j0=1
   j1=nseg
   jj=ifix(nseg*pr-pr)+1
   do while(j1-j0 > 1)
     if(pr .gt. pacul(jj) ) then
       j0=jj
     else
       j1=jj
     endif
     jj=(j0+j1)/2
   enddo
! linearly interpolating between j0 and j0+1
   xp(i)=x0+dx*(j0-1.0+(pr-pacul(j0))/(pacul(j0+1)-pacul(j0)))
 enddo

end subroutine random_rup_vel
!
!======================================================================
!
subroutine filter_field(source,clx,cly,pw2,dx,dy,nx,ny)
 implicit NONE
 integer,intent(IN):: nx,ny
 real, intent(IN):: clx,cly,pw2,dx,dy
 real source(nx*ny),param(nx*ny),speq(nx*2)
 integer:: nxy,nhfx,nhfy,ie,ke1,ke2,kr,k1,k2,j1,j2,i2,i,j,k,ix0,jy0
 real:: dwkx,dwky,cff,sq2,pw,wkx,wky,filter,ru1,ru2,rv1,rv2,cf

 nxy=nx*ny
 nhfx=nx/2
 nhfy=ny/2
! df
 dwkx=clx/(nx*dx)
 dwky=cly/(ny*dy)
 cff=0.0
 sq2=sqrt(1.0+cff)
 pw=0.5*pw2
 kr=0
!-----zero wavenumber kx
 i=1
 wkx=((i-1)*dwkx)**2
 do j=2,nhfy
   wky=((j-1)*dwky)**2
   filter=1./sqrt(1.+(wkx+wky)**pw2)
!   filter=1./(1.+(wkx+wky))**pw
   k1=2*((i-1)*nhfy+j)-1
   k2=k1+1
   param(k1)=-abs(source(kr+1)*filter/sq2)
   param(k2)=cff*source(kr+2)*filter/sq2
   kr=kr+2
 enddo
!-----maximum wavenumber kx
 i=1+nhfx
 wkx=((i-1)*dwkx)**2
 do j=2,nhfy
   wky=((j-1)*dwky)**2
   filter=1./sqrt(1.+(wkx+wky)**pw2)
!   filter=1./(1.+(wkx+wky))**pw
   k1=2*((i-1)*nhfy+j)-1
   k2=k1+1
   param(k1)=-abs(source(kr+1)*filter/sq2)
   param(k2)=cff*source(kr+2)*filter/sq2
   kr=kr+2
 enddo
!----zero wavenumber ky
 j=1
 wky=((j-1)*dwky)**2
 do i=2,nhfx
   wkx=((i-1)*dwkx)**2
   filter=1./sqrt(1.+(wkx+wky)**pw2)
!   filter=1./(1.+(wkx+wky))**pw
   k1=2*((i-1)*nhfy+j)-1
   k2=k1+1
   param(k1)=-abs(source(kr+1)*filter/sq2)
   param(k2)=cff*source(kr+2)*filter/sq2
   kr=kr+2
   ie=nx-i+2
   ke1=2*((ie-1)*nhfy+j)-1
   ke2=ke1+1
   param(ke1)= param(k1)
   param(ke2)=-param(k2)
 enddo
!----maximum wavenumber ky
 j=nhfy+1
 wky=((j-1)*dwky)**2
 do i=2,nhfx
   wkx=((i-1)*dwkx)**2
   filter=1./sqrt(1.+(wkx+wky)**pw2)
!   filter=1./(1.+(wkx+wky))**pw
   k1=2*i-1
   k2=k1+1
   speq(k1)=-abs(source(kr+1)*filter/sq2)
   speq(k2)=cff*source(kr+2)*filter/sq2
   kr=kr+2
   ie=nx-i+2
   ke1=2*ie-1
   ke2=ke1+1
   speq(ke1)= speq(k1)
   speq(ke2)=-speq(k2)
 enddo
!----
 i=1
 j=1
 param(1)=0.0
 param(2)=0.0
 kr=kr+1
!----
 i=nhfx+1
 j=1
 wkx=((i-1)*dwkx)**2
 wky=((j-1)*dwky)**2
 filter=1./sqrt(1.+(wkx+wky)**pw2)
! filter=1./(1.+(wkx+wky))**pw
 k1=2*((i-1)*nhfy+j)-1
 k2=k1+1
 param(k1)=-abs(source(kr+1)*filter)
 param(k2)=0.0
 kr=kr+1
!----
 i=1
 j=nhfy+1
 wkx=((i-1)*dwkx)**2
 wky=((j-1)*dwky)**2
 filter=1./sqrt(1.+(wkx+wky)**pw2)
! filter=1./(1.+(wkx+wky))**pw
 k1=2*((i-1)*nhfy+j)-1
 k2=k1+1
 speq(1)=-abs(source(kr+1)*filter)
 speq(2)=0.0
 kr=kr+1
!----
 i=nhfx+1
 j=nhfy+1
 wkx=((i-1)*dwkx)**2
 wky=((j-1)*dwky)**2
 filter=1./sqrt(1.+(wkx+wky)**pw2)
! filter=1./(1.+(wkx+wky))**pw
 speq(2*i-1)=source(kr+1)*filter
 speq(2*i)  =0.0
 kr=kr+1
!----center
 do i=2,nhfx
 do j=2,nhfy
   wkx=((i-1)*dwkx)**2
   wky=((j-1)*dwky)**2
   filter=1./sqrt(1.+(wkx+wky)**pw2)
!   filter=1./(1.+(wkx+wky))**pw
   ru1=source(kr+1)
   rv1=source(kr+2)
   ru2=source(kr+3)
   rv2=source(kr+4)
   k1=2*((i-1)*nhfy+j)-1
   k2=k1+1
   param(k1)=(ru1-ru2)*filter/2.0
   param(k2)=(rv1+rv2)*filter/2.0
   kr=kr+4
   ie=nx-i+2
   ke1=2*((ie-1)*nhfy+j)-1
   ke2=ke1+1
   param(ke1)=(ru1+ru2)*filter/2.0
   param(ke2)=(rv1-rv2)*filter/2.0
 enddo
 enddo
 call rlft3(param,speq,ny,nx,1,-1)

 ix0=1
 jy0=1
 cf=2./float(nxy)
 k=0
 i2=nx+ix0-1
 j2=ny+jy0-1
 do i=ix0,i2
 do j=jy0,j2
   k=k+1
   source(k)=param((i-1)*ny+j)*cf
 enddo
 enddo

end subroutine filter_field
!
!======================================================================
!
subroutine stress_drop(source,dx,dy,nx,ny)
 implicit NONE
 integer, intent(IN)::nx,ny
 real, intent(IN):: dx,dy
 real:: source(nx*ny),param(nx*ny),speq(nx*2)
 real:: dwkx,dwky,wkx,wky,wkxy,cf
 integer:: nxy,nhfx,nhfy,ie,ke1,ke2,k1,k2,i,j

 nxy=nx*ny
 nhfx=nx/2
 nhfy=ny/2
 dwkx=1.0/(nx*dx)
 dwky=1.0/(ny*dy)
 call rlft3(source,speq,ny,nx,1,+1)
!----
 j=nhfy+1
 wky=((j-1)*dwky)**2
 do i=1,nhfx+1
   wkx=((i-1)*dwkx)**2
   wkxy=sqrt(wkx+wky)
   k1=2*i-1
   k2=k1+1
   speq(k1)=speq(k1)*wkxy
   speq(k2)=speq(k2)*wkxy
   if(i > 1 .and. i <= nhfx) then
     ie=nx-i+2
     ke1=2*ie-1
     ke2=ke1+1
     speq(ke1)= speq(k1)
     speq(ke2)=-speq(k2)
   endif
 enddo
!----
 do i=1,nhfx+1
 do j=1,nhfy
   wkx=((i-1)*dwkx)**2
   wky=((j-1)*dwky)**2
   wkxy=sqrt(wkx+wky)
   k1=2*((i-1)*nhfy+j)-1
   k2=k1+1
   param(k1)=source(k1)*wkxy
   param(k2)=source(k2)*wkxy
   if(i.gt.1 .and. i.le.nhfx) then
     ie=nx-i+2
     ke1=2*((ie-1)*nhfy+j)-1
     ke2=ke1+1
     param(ke1)=source(ke1)*wkxy
     param(ke2)=source(ke2)*wkxy
   endif
 enddo
 enddo
 call rlft3(param,speq,ny,nx,1,-1)

 cf=2./float(nxy)
 do i=1,nxy
   source(i)=param(i)*cf
 enddo

end subroutine stress_drop
!
!=====================================================================
! The name of subroutine is misleading.
! What it calculates is the intergration of acceleration level from 0.5 to 0.9 freq_max
! in other word, keep a_hf between source model and stress parameter
!
! Chen Ji, 2020

subroutine Rd_energy_brune(energy,fc)
 use time_freq
 implicit none
 integer::  i,j
 real:: fc,energy
 real:: moment_rate(nphf),density,vs,coef,pi

 do i=1,nphf
   freq(i)=i*df
   moment_rate(i)=1./(1.+(freq(i)/fc)**2)
 enddo

 density=1.0
 vs=1.0
 pi=3.1415926
 coef=8.*pi/(10.*pi*density*vs**5)*df
 energy=0.0
 do j=nfe1,nfe2
!   energy=energy+(freq(j)*moment_rate(j))**2
!   energy=energy+alog(moment_rate(j)*freq(j)**2)
   energy=energy+(moment_rate(j)*freq(j)**2)
 enddo
! energy=exp(energy/float(nfe2-nfe1+1))*coef
 energy=energy*coef

end subroutine Rd_energy_brune
!
!=====================================================================
!
subroutine Rd_energy_synth(energy,rmt,fc)
 use time_freq
 implicit none
 integer:: i,j
 real:: rmt,fc,energy
 real:: svf(ntime),moment_rate(nphf)
 real:: density,brune,vs,coef,pi,tim,dtmt

 !svf=0.0
 call sum_point_svf(svf)
 open(19,file='calsvf_tim.dat',status='replace')
 write(19,*) ntime
 do i=1,ntime
    tim=(i-1)*dt
    write(19,*) tim,svf(i)/rmt
 enddo
 close(19)
!
 call realft(svf,ntime,-1)
 dtmt=dt/rmt
 do i=1,ntime
   svf(i)=svf(i)*dtmt
 enddo
 do i=1,nphf-1
   moment_rate(i)=sqrt(svf(2*i+1)**2+svf(2*i+2)**2)
 enddo
 moment_rate(nphf)=svf(2)

 open(19,file='calsvf.dat')
 write(19,*) nphf
 do i=1,nphf
   brune=1./(1.+(freq(i)/fc)**2)
   write(19,*) freq(i),moment_rate(i),brune
 end do
 close(19)
!
 density=1.0
 vs=1.0
 pi=3.1415926
 coef=8.*pi/(10.*pi*density*vs**5)*df
 energy=0.0
 do j=nfe1,nfe2
!   energy=energy+(freq(j)*moment_rate(j))**2
!   energy=energy+alog(moment_rate(j)*freq(j)**2)
   energy=energy+(moment_rate(j)*freq(j)**2)
 enddo
 energy=energy*coef

end subroutine Rd_energy_synth
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! The name of subroutine is misleading.
! What it calculates is the intergration of acceleration level from 0.5 to 0.9 freq_max
! in other word, keep a_hf between source model and stress parameter
!
! Chen Ji, 2020

subroutine Rd_energy_dcf(energy,fc1,fc2)
 use time_freq
 implicit none
 integer::  i,j
 real:: fc1,fc2,energy
 real:: moment_rate(nphf),density,vs,coef,pi

 do i=1,nphf
   freq(i)=i*df
   moment_rate(i)=1./((1.+(freq(i)/fc1)**4.0)**0.25)/((1.+(freq(i)/fc2)**4.0)**0.25)
 enddo

 density=1.0
 vs=1.0
 pi=3.1415926
 coef=8.*pi/(10.*pi*density*vs**5)*df
 energy=0.0
 do j=nfe1,nfe2
!   energy=energy+(freq(j)*moment_rate(j))**2
!   energy=energy+alog(moment_rate(j)*freq(j)**2)
   energy=energy+(moment_rate(j)*freq(j)**2)
 enddo
! energy=exp(energy/float(nfe2-nfe1+1))*coef
 energy=energy*coef

end subroutine Rd_energy_dcf
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!=====================================================================
!
subroutine Rd_energy_synth_dcf(energy,rmt,fc1,fc2)
 use time_freq
 implicit none
 integer:: i,j
 real:: rmt,fc1,fc2,energy
 real:: svf(ntime),moment_rate(nphf)
 real:: density,dcf,vs,coef,pi,tim,dtmt

 !svf=0.0
 call sum_point_svf(svf)
 open(19,file='calsvf_tim.dat',status='replace')
 write(19,*) ntime
 do i=1,ntime
    tim=(i-1)*dt
    write(19,*) tim,svf(i)/rmt
 enddo
 close(19)
!
 call realft(svf,ntime,-1)
 dtmt=dt/rmt
 do i=1,ntime
   svf(i)=svf(i)*dtmt
 enddo
 do i=1,nphf-1
   moment_rate(i)=sqrt(svf(2*i+1)**2+svf(2*i+2)**2)
 enddo
 moment_rate(nphf)=svf(2)

 open(19,file='calsvf.dat')
 write(19,*) nphf
 do i=1,nphf
!    brune=1./(1.+(freq(i)/fc)**2)
    dcf=1./((1.+(freq(i)/fc1)**4.0)**0.25)/((1.+(freq(i)/fc2)**4.0)**0.25)
   write(19,*) freq(i),moment_rate(i),dcf
 end do
 close(19)
!
 density=1.0
 vs=1.0
 pi=3.1415926
 coef=8.*pi/(10.*pi*density*vs**5)*df
 energy=0.0
 do j=nfe1,nfe2
!   energy=energy+(freq(j)*moment_rate(j))**2
!   energy=energy+alog(moment_rate(j)*freq(j)**2)
   energy=energy+(moment_rate(j)*freq(j)**2)
 enddo
 energy=energy*coef

end subroutine Rd_energy_synth_dcf
!
!======================================================================
!
subroutine coherent(cr_srv,nsub,slip,rstm)
 implicit NONE
 integer, intent(IN):: nsub
 real, dimension(nsub), intent(IN):: slip, rstm
 real, intent(OUT):: cr_srv
 integer:: i
 real:: avg_s,avg_r,s2,r2,sr,ds,dr
!
 avg_s=0.0
 avg_r=0.0
 do i=1,nsub
   avg_s=avg_s+slip(i)
   avg_r=avg_r+rstm(i)
 enddo
 avg_s=avg_s/float(nsub)
 avg_r=avg_r/float(nsub)
!
 s2=0.0
 r2=0.0
 sr=0.0
 do i=1,nsub
   ds=slip(i)-avg_s
   dr=rstm(i)-avg_r
   s2=s2+ds*ds
   r2=r2+dr*dr
   sr=sr+ds*dr
 enddo
 cr_srv=sr/sqrt(s2)/sqrt(r2)

end subroutine coherent
