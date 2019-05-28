SUBROUTINE simcoda(station)
!---------------------------------------------------------------------
!
! Description:
!
!   Simulate S-wave coda for synthetic Green's function using Zeng's 
!   multiple S-to-S scattering theory
!
! Dependencies:
!
!   Subroutines fork, scoda, coda, ener2, evwi,four1d
!
! Fundamental input parameters:
!
!   nfcoda    - fft points for the coda envelope            m=mx
!   merr      - tolerance of the coda envelope              err=errx
!   abscoeff  - absorption coefficient                      gi=gix
!   scatcoeff - scattering coefficient  	                gs=gsx
!   fmax      - maximum frequency for coda waves   
!   Q         - Q0 attenuation factor {Q(f)=Q0*f**fdec}
!   fdec      - frequency decay of Q(f) factor	            cn=cnx
!   iseed     - seed number                                
!   nscat     - number of scattering wavelets	
!   hpass	  - highpass corner of the cosine filter        fl=flx
!   trans     - transition bandwidth of the filter          flb=flbx
!   t_len	  - length of time series                       twin=twinx
!   npts	  - number of points computed for coda envelope nt=ntx
!   t0	      - time of first arrival 
!   kappa     - kappa factor akp=akpx
!   sr_hypo   - station-hypocenter distance          	    dist=distx
!   aveVs	  - average S-speed between S-R                 vs=vsx
!   siteVs    - S-wave velocity at the site	                vs0=vs0x
!   srcVs     - S-wave velocity at the source               vs1=vs1x
!   siteR     - density at the site                         rh0=rh0x
!   srcR      - density at the source                       rh1=rh1x
!
! References:
!
!   Zeng, Anderson and Su (1995).  Subevent rake and random 
!   scattering effects in realistic strong ground motion 
!   simulation, Geophy. Res. Lett., 22 17-20
!
! Note:
!
!   Allmost all needed parameters are passed by modules
!
! Authors: W. Imperatori, M. Mai, Y. Zeng 
!
! Modified: July 2009 (v1.4)
!
! Updated: December 2012 (v1.4.1)
!    In scoda subroutine, the last do loop has changed
!    from srcV1.4_OLSEN
!
! Updated: March 2013 (v1.4.2)
!    Add chack if Fmax < Fnyq after calculating nw (Fnyq) and nwm (Fmax).
!    When fmax is a certain number, nwm will be larger than npts.
!    This leads u2, u1 and u having NaN or Infinity, generating NaN
!    in the horizontal components in the BroadBand-synthetics.
!    This has happened when fmax is 20Hz.
!
! Updated: April 2014 (v1.5.5)
!    tmp_npts and tmp_lf_len are in the module tmp_para.
!
! Updated: December 2014 (v1.5.5.2)
!    The dimension of p and t must be dimenstion(ncoda)
!    in coda and scoda
!
! Updated: February 2019 (v2.0)
!    Change dt and df computations using v_npts.

use constants; use def_kind; use scattering; use source_receiver
use waveform, only: lf_len, scattgram, v_npts
use interfaces, only: four1d
use tmp_para

implicit none

! local variables
integer(kind=i_single),intent(in)       :: station 
!complex(kind=r_single),dimension(npts)  :: u2
!complex(kind=r_single),dimension(npts,3):: u1
!complex(kind=r_single),dimension(2*npts):: u
real(kind=r_single)                     :: kappa_local,dt,df,dw,w,x
real(kind=r_single)                     :: c1,c2,c3,tmp,band,db,Qt
integer(kind=i_single)                  :: k,nw,nwm,j,l,ml,mlb
integer(kind=i_single)                  :: mmb,i

! add
complex(kind=r_single),allocatable,dimension(:)   :: u2
complex(kind=r_single),allocatable,dimension(:,:) :: u1
complex(kind=r_single),allocatable,dimension(:)   :: u

!---------------------------------------------------------------------

! compute some parameters
kappa_local=kappa(station)*0.5
!dt = lf_len/(npts-1)
dt = lf_len/(v_npts-1)

! df is set to avoid time aliasing
!df = 1 / (2 * npts * dt)                                        
df = 1 / (2 * v_npts * dt)                                        
dw=pi_double*df

! index of Fnyq 
nw = (npts) + 1

! index of Fmax   (Fmax < Fnyq, guaranteed by subroutine sampling_check)
nwm = nint(fmax/df) + 1

! check Fmax < Fnyq (nwm < nw)
if (npts .lt. nwm) then
   !nwm=npts
   print*,'change Fmax'
   print*,'fmax,npts,nw,df,nwm',fmax,npts,nw,df,nwm
endif

! allocate u2, u1, and u
if (.not. allocated(u2)) then
   allocate(u2(npts), u1(npts,3), u(npts*2))
endif

u2 = zero

do j=2,nwm
   w=(j-1)*dw
   Qt=1.0/(Q*(w/pi_double)**fdec)
   u2(j)=((w-w*log(w/pi_double)*Qt/pi)/aveVs)*cmplx(0.5*Qt,1.0)
enddo

u1 = zero

do l=1,nscat
   x=aveVs*(t0+random_array(1,l)*lf_len)
   c1=random_array(2,l)*2.0-1.0
   c2=random_array(3,l)*2.0-1.0
   c3=random_array(4,l)*2.0-1.0

   do j=1,nwm
      tmp=exp(-x*u2(j))
      u1(j,1)=u1(j,1)+c1*tmp
      u1(j,2)=u1(j,2)+c2*tmp
      u1(j,3)=u1(j,3)+c3*tmp
   enddo

enddo  
      
!  filter the low frequency
ml = nint(hpass/df) + 1          
mlb = nint((hpass-trans)/df) + 1 

do i=1,3

   u = zero

   do j=mlb,ml-1
      w=(j-1)*dw
      c1 = cos(pi_half * ( (ml-j) / (ml-(mlb-1)) )) *df   
      u(j)=u1(j,i)*exp(cmplx(-w*kappa_local,w*t0))*c1
   enddo
 
   do j=ml,nwm
      w=(j-1)*dw
      u(j)=u1(j,i)*exp(cmplx(-w*kappa_local,w*t0)) *df
   enddo

   do j=2,nw-1                        
      u(2*npts+2 - j) = conjg(u(j))
   enddo  

   call four1d(u,-1)

   do j=1,npts
      scattgram(j,i)=real(u(j))
   enddo

enddo  

dw=dw*0.5

! simulate strong motion coda waves for each subfaults
call scoda(sr_hypo(station),dt,station)

END SUBROUTINE simcoda

!=============================================================================

SUBROUTINE fork(lx,cx,signi)

use def_kind

implicit none

integer(kind=i_single),intent(in)                 :: lx
real(kind=r_single),intent(in)                    :: signi
complex(kind=r_single),dimension(lx),intent(inout):: cx

!local
complex(kind=r_single):: carg,cw,ctemp
integer(kind=i_single):: j,i,m,l,istep

j=1

do i=1,lx

   if(i.le.j) then
      ctemp=cx(j)
      cx(j)=cx(i)
      cx(i)=ctemp
   end if

   m=lx/2 
20 if(j.le.m) goto 30
   j=j-m
   m=m/2 
   if(m.ge.1) goto 20
30 j=j+m

enddo

l=1 
40 istep=2*l

do m=1,l
   carg=(0.,1.)*(3.14159265*signi*float(m-1))/float(l)
   cw=cexp(carg)

   do i=m,lx,istep
      ctemp=cw*cx(i+l) 
      cx(i+l)=cx(i)-ctemp
      cx(i)=cx(i)+ctemp 
   enddo

enddo

l=istep
if(l.lt.lx) goto 40

END SUBROUTINE fork

!=============================================================================

SUBROUTINE scoda(dist,dt,station)

!  Program : scoda.f
!  Purpose : simulate strong motion S-wave coda
!                                by Yuehua Zeng at UNR, Oct 1993
!
! Updated: December 2014 (v1.5.5.2)
!    Add if (n1 .gt. ncoda) n1=ncoda
!    to check if n1 becomes greater than ncoda for p(n1)
!
! Updated: December 2016 (v1.6.2)
!    Change dt1 computation. 
!    Add shifting scattgram by minimum initiation time,
!    (tinit, tmp_scatt, sft).
!
! Updated: February 2019 (v2.0)
!    Back to use the original calling coda routine, and no use tmp_para.
!
use constants; use def_kind; use scattering
use waveform, only: lf_len, scattgram
! use tmp_para
use vel_model, only: tinit

implicit none

real(kind=r_single),intent(in):: dist,dt
integer(kind=i_single),intent(in):: station

!local
real(kind=r_single),dimension(ncoda):: p,t
real(kind=r_single),dimension(npts) :: tmp_scatt ! for scattgram shifting, v162
real(kind=r_single):: c0,c1,dt1,t1,pp
integer(kind=i_single):: i1,i2,ii,i,i0,n1
integer(kind=i_single):: sft ! number of scattgram shifting, v162

!	parameter (mm1=8192,mm2=800)
!	dimension pcoda(mm1,3),t(mm2),p(mm2)
!	common/codaw/n,m,err,gi,gs,fl,flb,rh0,rh1,vs0,vs1,vs,twin,dt
!	pi=3.1415926

!  generate the coda waves (1.05 is to avoid noise at time-series end) 

call coda(1.05*lf_len,dist,t,p,dt)

c0=0.85/(pi*4.*dist*aveVs)*sqrt(srcVs*srcR/(siteVs(station)*siteR(station)))

c1=sqrt(3.0*lf_len/float(nscat))*c0

! dt1=t(2)-t(1) ! change v162

i1=ifix(t0/dt+0.000001)
i2=ifix(t(1)/dt+0.000001)

do ii=1,3
   do i=1,ncoda
      if(t0.le.t(i))goto 65
   enddo
 
65 n1=i
   i0=0
   t1=float(i1-1)*dt
 
   do i=1,i2-i1+1
      i0=i0+1
      t1=t1+dt
      scattgram(i,ii)=0.0
   enddo

! ----- changed from srcV1.4_OLSEN ---------------------------------------------

   do i=i0+1,npts
       t1=t1+dt
       if(t1.gt.t(n1)) then
          n1=n1+1
          ! check if n1 <= ncoda (v1552)
          if (n1 .gt. ncoda) then
             ! n1=ncoda
             print*,'should be n1<=ncoda in scoda,n1,ncoda, STOP'
          endif
       endif
       if(n1>1) then
          dt1=t(n1)-t(n1-1) ! add v162
          pp=p(n1-1)+(p(n1)-p(n1-1))*(t1-t(n1-1))/dt1
       else
          pp=p(n1)
       endif
       scattgram(i,ii)=pp*(c1*scattgram(i,ii))
   enddo

!   do i=i0+1,npts
!      t1=t1+dt
!      if(t1.gt.t(n1))n1=n1+1
!      pp=p(n1-1)+(p(n1)-p(n1-1))*(t1-t(n1-1))/dt1
!      scattgram(i,ii)=pp*(c1*scattgram(i,ii))
!   enddo

! ----- end of change ----------------------------------------------------------

! scattgram shift by minimum initiation time (for multiple-segment-rupture)
   if(tinit > 0) then
      sft = nint(tinit/dt)
      tmp_scatt(1:npts)=0.
      tmp_scatt(sft+1:npts)=scattgram(1:npts-sft,ii)
      scattgram(:,ii)=tmp_scatt
   endif

enddo   

END SUBROUTINE scoda

!=============================================================================

SUBROUTINE coda(tw,r,t,p,ddt)

use constants; use def_kind; use scattering

implicit none

real(kind=r_single),intent(in)                 :: tw,r
real(kind=r_single),dimension(ncoda),intent(out):: p,t
real(kind=r_single),intent(in)                 :: ddt

!  Program for calculating the scattered wave energy using equation
!  (24) of Zeng et al. (1991).  Written by Yuehua Zeng in May 1990
!  and updated by Yuehua Zeng on Feb. 9, 1993.
!
! Updated: December 2014 (v1.5.5.2)
!    change t(i) computation that "i" in t(i) must be ncoda 
!    or less than ncoda.
!
! Updated: December 2016 (v1.6.2)
!    Use aveVp for dt and t(1) computations, instead of using aveVs.
!    Add c4 and ratio for a new coda energy computation.
!------------------------------------------------

! local variables
complex(kind=r_single):: atinv,atinv1,atinv3,akw,wg
complex(kind=r_single):: ei2,aker,Es,Es1,Es2
complex(kind=r_single),dimension(1024):: au

real(kind=r_single):: ep,g,dw,dk,dt,w,ak,sinkr,wi
real(kind=r_single):: c0,c1,c2,c3,gamma,r2,vt,tau
real(kind=r_single):: c4,ratio ! add v162
real(kind=r_single):: pd,di,coef

integer(kind=i_single):: nd,ni,m2,i,nk,ik,j,n1

!	subroutine coda(n,m,gs,gi,vs,tw,r,err,t,p,ddt)
!	dimension t(n),p(n)
!	complex au(1024),atinv,atinv1,atinv3,akw,wg,one,zero,ei,
!     +		ei2,aker,Es,Es1,Es2
!
!  compute parameters
!	pi=3.14159265
!	one=cmplx(1.0,0.0)
!	zero=cmplx(0.0,0.0)

ei2=0.5/zeta
ep=1.0e-20
g=scatcoeff+abscoeff
dw=2.0*pi/tw
nd=2*ifix(0.25*(r+aveVs*tw)/r+1)
dk=pi/(r*float(nd))
ni=ifix(20.0/nd+1)*nd+nd/2
m2=nfcoda/2
dt=(tw-r/aveVp)/float(ncoda)
t(1)=r/aveVp
!dt=(tw-r/(aveVs*sqrt(3.0)))/float(ncoda) ! change v162
!t(1)=r/(aveVs*sqrt(3.0)) ! change v162

do i=2,ncoda
   t(i)=t(i-1)+dt
enddo

! check i in t(i) must be ncoda or less than ncoda (v1552)
i=ifix((r/aveVs-t(1))/dt)+1
!t(i+1)=r/aveVs+ddt
!t(i)=r/aveVs-ddt
if (i .lt. ncoda) then
   t(i+1)=r/aveVs+ddt
   t(i)=r/aveVs-ddt
elseif (i .ge. ncoda) then
   t(ncoda)=r/aveVs+ddt
   t(ncoda-1)=r/aveVs-ddt
endif

!  calculate the power spectrum for scattering order higher than two
call evwi

w=-dw

do i=1,m2
   w=w+dw
   wg=cmplx(g+wi/aveVs,w/aveVs)
   ak=0.0
   Es=zero
   Es1=zero
   Es2=zero
   nk=ni
   ik=0

   do j=1,10000
      ak=ak+dk
      akw=cmplx(0.0,ak)/wg
      atinv=ei2*clog((one+akw)/(one-akw))
      atinv1=scatcoeff*atinv/ak
      atinv3=atinv1*atinv1*atinv1
      aker=atinv3*atinv/(one-atinv1)
      sinkr=sin(ak*r)
      Es=Es+sinkr*aker
      if(cabs(aker/(Es+ep)).lt.merr.and.j.gt.20) goto 25
      if(j.eq.nk)then
         ik=ik+1
         Es2=Es2+Es-sinkr*aker*0.5
         if(ik.eq.6)then
            if(cabs((Es2-Es1)/(Es2+ep)).lt.merr)then
               Es=(Es1+Es2)/12.0
               goto 25
            endif
            ik=0
            Es1=Es2
            Es2=zero
         endif
         nk=nk+nd
      endif
   enddo

25 au(i)=Es*dk

enddo

n1=nfcoda+2

do i=1,m2
   au(n1-i)=conjg(au(i))
enddo

au(m2+1)=zero

call fork(nfcoda,au,1.0)

w=0.0
dw=wi*tw/float(nfcoda-1)

do i=2,nfcoda
   w=w+dw
   au(i)=au(i)*exp(w)
enddo

!  loop for calculating time domain scattering energy decay
dt=tw/float(nfcoda)
 c0=0.25/(aveVs*pi*r*r)
 c1=0.25*scatcoeff/(aveVs*r*pi*c0)
 c2=scatcoeff*scatcoeff/(16.0*r*pi*c0)
 c3=0.5/(aveVs*r*pi*pi*tw*c0)
gamma=2.0/(9.0*sqrt(3.0))
r2=r*r
! add c4 as (Vs/Vp)**4, eta_sp = 0.11*eta_s by Zeng et al., 1995, v162
c4=(aveVs/aveVp)
c4=c4**4
! 19.88 km is a reference distance for a new tau, v162
ratio=7.5*r/19.88

do i=1,ncoda
   vt=aveVs*t(i)
   if(vt.lt.r)then
      !tau=exp(-g*vt) ! original
      !p(i)=gamma*c1*alog((vt+r)/(r-vt))/t(i)*tau ! original
      tau=exp(g*vt)/ratio ! v162
      p(i)=c4*gamma*c1*alog((vt+r)/(r-vt))/t(i)*tau ! v162
      goto 55
   endif
   tau=r/vt
 
   call ener2(0.01)

   n1=ifix(t(i)/dt)+1
   pd=(real(au(n1))+(real(au(n1+1))-real(au(n1)))  &
      *(t(i)-dt*float(n1-1))/dt)*exp(abscoeff*vt)
   tau=exp(-scatcoeff*vt)
   p(i)=c1*alog((vt+r)/(vt-r))/t(i)*tau+c2*coef*tau+c3*pd
   p(i)=p(i)*exp(-(0.03-abscoeff)*vt)
55 p(i)=sqrt(abs(p(i))) ! avoid negative p(i) in sqrt()

enddo


!..........................................................................

CONTAINS

   SUBROUTINE ener2(dd)
   !  calculate coefficient of the second-order scattered energy term

   implicit none

   real(kind=r_single),intent(in):: dd
   !local
   integer(kind=i_single)        :: n,i
   real(kind=r_single)           :: dt,t,aux


   n=ifix(tau/dd)+1
   dt=tau/float(n)
   coef=0.0
   t=0.0

   do i=1,n
      t=t+dt
      aux=alog((1.0+t)/(1.0-t))
      aux=aux*aux
      coef=coef+aux
   enddo

   coef=3.0*(0.5*aux-coef)*dt
   coef=coef+tau*pi*pi

   END SUBROUTINE ener2
!.....................................................

   SUBROUTINE evwi
   !  estimating the value of imaginery frequency wi

   implicit none

   real(kind=r_single)   :: dd,vt,ud,du,coef,u
   real(kind=r_single)   :: aux,di,dimax
   integer(kind=i_single):: i,m,j 


   dd=0.01
   vt=r

   do i=1,2
      vt=vt+(aveVs*tw-r)*(i-1)
      ud=sqrt(0.75*scatcoeff*vt)
      m=ifix(ud/dd)+1
      du=ud/float(m)
      coef=0.5
      u=0.0

      do j=1,m
         u=u+du
         aux=exp(-u*u)
         coef=coef+aux
      enddo

      coef=(coef-0.5*aux)*du-ud*exp(-ud*ud)
      u=0.75*scatcoeff/vt
      aux=(1.-exp(-scatcoeff*vt)*(1.+scatcoeff*vt))/coef
      di=aux*exp(1.5*alog(u/real(pi))-u*r*r-abscoeff*vt)
 
      if(i.eq.1)then
         dimax=di
      end if

   enddo

   wi=alog(500*di/dimax)/tw

   END SUBROUTINE evwi 

END SUBROUTINE coda
