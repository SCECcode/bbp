! spfield.f90
! The subroutines in this file are for generating the spatial 
! distribution of source parameters on main-fault
!
! Written by Pengcheng Liu
! Copyright (c) 2005 by Pengcheng Liu
!
!======================================================================
subroutine mainfault(dip,freq_max,rv_avg,ratio_rise,nb_taper_TRBL)
! Set the initial parameters for generating the source paramters
! flx,flw-- fault length (in strike) and width (in dip)
 use sp_sub_f
 use fdtim_2d
 use time_freq
 implicit NONE
 integer,dimension(4), intent(IN):: nb_taper_TRBL
 real, intent(IN):: dip,freq_max,rv_avg,ratio_rise
 integer:: i,j,k,lyr,nb,ll,nx1,ny1
 real:: rdip,yij,hk,sumh,sv,sm,rat,ttmp,vs_avg
!
 cft1=ratio_rise
 nsum=nsubx*nsuby
 dx=flx/nsubx
 dy=fwy/nsuby
 smoment= Moment_o*1.e-15
!
 grid_2d=amin1(0.2,  amin1(dx,dy))
 nx_2d=int((nsubx-1)*dx/grid_2d)+2
 nz_2d=int((nsuby-1)*dy/grid_2d)+2
 ncxp=int(x_hypc/dx) + 1
 ncyp=int(y_hypc/dy) + 1
 cxp=(ncxp-0.5)*dx
 cyp=(ncyp-0.5)*dy
 x_hypc=cxp
 y_hypc=cyp
!        
 if(allocated(slip)) deallocate(slip,rstm,rptm,rpvel,beta,amu, &
                     taper,rtx,rtz,rake_prt,dip_prt,amz_prt)
 allocate(slip(nsum),rstm(nsum),rptm(nsum),rpvel(nsum), &
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
 end select
!
 slip_cly=amin1(slip_cly, 0.7*fwy)
 slip_clx=amin1(slip_clx, 0.7*flx)
 slip_clx=amin1(slip_clx, amax1(2.0*slip_cly, flx/4.0))
 if (flx < fwy) then
   ttmp=slip_cly
   slip_cly=slip_clx
   slip_clx=ttmp
 endif
if(Mw >7.3) print *, Mw, slip_clx, slip_cly
 sp1 = 2.0
 sq1 = 5.0
 slip_pw  = 2.0
 hfrd=smoment*fc_main*fc_main
!
! Setting the parameters for computing random Risetime filed
! 
 rs_clx = amin1(slip_clx, 10.0)
 rs_cly = amin1(slip_cly, 10.0)
 rs_max=6.0
 rp1= 3.0
 rq1= 2.0
 rs_pw  = 2.0
!
! Setting the parameters for computing random rupture-velocity filed
!
 rv_clx = amin1(slip_clx, 10.0)
 rv_cly = amin1(slip_cly, 10.0)
!
 rv_pw = 2.0
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
! between slip-amplitude and rupture-velocity.
!
 cr_sl_rs = 0.6!0.6
 cr_sl_vc = 0.0!0.3
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
 df=1./(ntime*dt)
 nphf=ntime/2
 nfe1=int(0.5*freq_max/df)+1
 nfe2=int(0.9*freq_max/df)
 if(allocated(freq)) deallocate(freq)
 allocate(freq(nphf))
 
end subroutine mainfault
!
!======================================================================
!
subroutine random_field(idum1,idum2,idum3)
 use sp_sub_f
 use fdtim_2d
 implicit NONE
 integer, intent(INOUT):: idum1,idum2,idum3
 integer, parameter:: itrmax=1
 logical:: lloop
 integer:: i,j,iset,niter,is_stress
 real:: sum1,cmti,correct,cr0,cr1,cr2,csqt,cv0,cv1,cv2, &
        vrup,er_syn,er_target,cr_sr,cr_sv,gasdev
 real,dimension(nsum):: gasv1,gasv2,gasv3,stdrop
!
! Generating Normal distributed white noise
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
!
! mapping to the desired distributions using NORTA technique 
!
! slip amplitude
 slip=gasv1
 call filter_field(slip,slip_clx,slip_cly,slip_pw,dx,dy,nsubx,nsuby)
 call random_slip_am_C(slip,nsum)
! call random_slip_am_B(slip,nsum,sp1,sq1)
! 
 stdrop=slip
 slip=slip*amu*taper
 correct=smoment/sum(slip(1:nsum))
 slip=slip*correct

!
! Stress
! stdrop=slip
! call stress_drop(stdrop,dx,dy,nsubx,nsuby)
!
! rise time
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
   sum1=0.0
   do i=1,nsum
     sum1=sum1+(slip(i)/rstm(i)**2)**2
   enddo
   correct=sqrt(sqrt(sum1)/hfrd)
   rstm=rstm*correct
   call coherent(cr_sr,nsum,stdrop,rstm)
!   call coherent(cr_sr,nsum,slip,rstm)
!
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
!
! Estimating peak slip rate
 call peak_slip_rate(stdrop)
!
! rupture velocity
 cv0=cr_sl_vc
 cv1=0.0
 cv2=1.0
 niter=0
 lloop=.true.
 do while(lloop)
   niter=niter+1
   csqt=sqrt(1.-cv0*cv0)
   do i=1,nsum
     rptm(i)=cv0*gasv1(i)+csqt*gasv3(i)
   enddo
   call filter_field(rptm,rv_clx,rv_cly,rv_pw,dx,dy,nsubx,nsuby)
   call random_rup_vel(rptm,nsum,vmin_ratio,vmax_ratio,vp1,vq1)
   do i=1,nsum
     vrup=rptm(i)*beta(i)
     rpvel(i)=vrup
     rptm(i)=1.0/vrup
   enddo
   call coherent(cr_sv,nsum,stdrop,rpvel)
!
   write(*,*) 'cross rvel, peak-slip-rate',cr_sv,cv0,cr_sl_vc
!
   if(abs(cr_sv-cr_sl_vc) < 0.01 .or. niter > itrmax) then
     call random_rup_tim(nsum,nsuby,nsubx,dx,dy,rptm)
     lloop=.false.
   endif
   if(cr_sv > cr_sl_vc) then
     cv2=cv0
   else
     cv1=cv0
   endif
   cv0=0.5*(cv2+cv1)
 enddo
!
! Scale the rise time
 niter=0
 lloop=.true.
 call Rd_energy_brune(er_target,fc_main)
 do while(lloop)
   niter=niter+1
   call Rd_energy_synth(er_syn,smoment,fc_main)
   correct=er_syn/er_target
!
   write(*,*) 'energy ',er_target,er_syn,correct,niter
!
   if(abs(correct-1.) < 0.03 .or. niter > itrmax) then
     lloop=.false.
   else
     rstm=rstm*sqrt(correct)
   endif
 enddo 
!
! scale the unit of seismic moment to N-m
  slip=slip*1.0e+15
!
! Perturbate the rake dip, and azimuth angle
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
!
subroutine random_rup_tim(nsum,nsuby,nsubx,dx,dy,rptm)
 use fdtim_2d
 implicit NONE
 integer, intent(IN):: nsum,nsuby,nsubx
 real:: dx,dy,rptm(nsum)
 integer:: i,j,k,jx,jz,jxz,jxzx
 real:: rx,rz,rz1
 real, dimension(nx_2d, nz_2d):: rtm_2d,slow
!
 do k=1,nsum
   rptm(k)=sqrt((rtx(k)-cxp)**2+(rtz(k)-cyp)**2)*rptm(k)
 enddo
 return
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

end subroutine random_rup_tim
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

 dwkx=clx/(nx*dx)
 dwky=cly/(ny*dy)
 cff=0.0
 sq2=sqrt(1.0+cff)
 pw=0.5*pw2
 kr=0
!-----
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
!-----
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
!----  
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
!----  
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
!----  
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
!
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
