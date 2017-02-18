SUBROUTINE set_stf(station)
!----------------------------------------------------------------------------
!
! Description:
!
!   Prepares the source-time-function (or more precisely, the moment-rate 
!   function) to be used during convolution. The following functions are 
!   available: box-car, triangular, Yoffe, Liu, Dreger and user-defined.
!   The code calculates amplitude and rise time of the source-time-function 
!   as follows:
!   a) Computation of rupture area (Ar) from the relation of moment magnitude
!      (Mw) vs rupture area --> internal subroutine CalcAfrM
!   b) Computation of seismic moment (Mo) from the relation of Mo vs Mw --> 
!      internal subroutine fMw2MoN
!   c) Once Mo and Ar are estimated, it calculates slip (D) from the relation
!      Mo=mu*D*Ar
!   d) Computation of Rise-Time (Tr) from the relation of Tr vs Mo --> 
!      internal subroutine CalcTrfromM 
!   e) Computation of source-time-function amplitude (Amp), once rise-time 
!      and slip have been estimated
!
! Dependencies:
!
!   Subroutines CalcAfrM, fMw2MoN, CalcTrfromM, stf_box, stf_tri, stf_yof, 
!   stf_dreg, stf_liu.
!
! References:
!
!   Numerical Recipes
!
! Notes:
!
!   This section is the same for both point-sources and extended-sources.
!
! Authors: W. Imperatori, B. Mena
!
! Modified: January 2009 (v1.3)
!
! Updated: December 2012 (v1.4.1)
!   change taper_len.
!
! Updated: March 2013 (v1.4.2)
!   Add srf_read subroutine to read srf to obtain subfault informations
!   used in acc_spec in composition.f90.
!
!   Add a new source time function in stf_new subroutine. 
!
! Updated: October 2013 (v1.5.2)
!   Add tmp_lf_len and tmp_npts for STF, for output decimation.
!
! Updated: April 2014 (v1.5.5)
!    tmp_npts and tmp_lf_len are in the module tmp_para.
!
! Updated: September 2014 (v1.5.5.1)
!   Change SUBROUTINE set_stf to SUBROUTINE set_stf(station)
!
use constants; use def_kind; use earthquake; use flags
use stf_data; use scattering, only: npts,fmax
use waveform; use fault_area; use tmp_para

implicit none

integer(kind=i_single),intent(in)           :: station

! rigidity, average slip, stf's amplitude, time-step   
real(kind=r_single)                         :: mu,D,Amp,dt
! indexes, taper length 
integer(kind=i_single)                      :: i,j,taper_len
! dummy array
real(kind=r_single),allocatable,dimension(:):: h 
! test dummy array
real(kind=r_single),allocatable,dimension(:):: stf_dreg_t,stf_liu_t 
! rupture area, seismic moment and rise-time (for internal subroutines)
real(kind=r_single)                         :: Mo,Tr

!----------------------------------------------------------------------------

! returns Ar (rupture Area [km^2]) from Mw and mechanism
call CalcAfrM

! returns Mo [Nm] (Seismic Moment) from Mw (moment magnitude)	
call fMw2MoN

! Tr: rise time [sec] from Mo
call CalcTrfromM

! mu: rigidity. Split in two, otherwise can cause run-time error with some compilers
mu=3.3*(10**5); mu=mu*(10**5)   

! average D: slip [m]  
D=Mo/(mu*Ar*(10**6))

! calculating amplitude of velocity source-time-function (moment rate function) of 
! box-car type.
! If triangular source-time-function is chosen, amplitude correction is made inside
! the triangular-function subroutine, so take the same amplitude as for box-car type.
! Amp [m/sec] (i.e. velocity)
Amp=D/Tr


! now create the moment-rate function; for simplicity we fix it here to be 15 sec. long

! total length of source time function (total) = 15 sec (may change for user-defined stf)
total = 15.0
! stf starts at 0 seconds (ton = 0.0)
ton = 0.0

! time-step for STF, assuming maximum frequency for coda waves as Nyquist frequency
if (npts == tmp_npts) dt = lf_len/(npts-1)
if (npts .gt. tmp_npts) dt = tmp_lf_len/(tmp_npts-1)   
!dt = 1 / (2 * fmax)   !it's the same as for coda waves time-series

! number of points of the stf. 
! NOTE: it can be changed if user-defined stf is selected
npts_stf = nint(total/dt) + 1     ! +1 is due to include 0

! read srf file
call srf_read

! allocate arrays for stf amplitudes and its time vector 
if(.not.allocated(t_stf)) allocate(t_stf(npts_stf))
if(.not.allocated(stf)) allocate(stf(npts_stf)) 

! compute selected source-time-function
select case (trim(stf_name))
case('box')
   call stf_box           !BOXCAR
case('tri')
   call stf_tri           !TRIANGULAR
case('yof')
   call stf_yof           !YOFFE
case('dreg')
   call stf_dreg          !DREGER
case('liu')
   call stf_liu           !LIU
case('new')
   call stf_new           !new stf
!case default
!   call MomRateFunct(dt)  !USER-DEFINED
end select

! scale obtained source-time-function by rupture area (only for point-source case)
if(ext_flag == 0) stf=stf*Ar  

! create stf's time-vector

t_stf(1)=0.            
do i=1,(npts_stf-1)
   t_stf(i+1)=t_stf(i)+dt
enddo

! setting up the taper (Hanning window): only the end (last 20%) of convolved time series will
! be tapered
!taper_len=nint(npts*0.5)

! setting up the taper (Hanning window): only the end (10 sec) of convolved time series will
! be tapered
taper_len=nint(10.0/dt)

! allocate local array for building taper
if(.not.allocated(han_win)) allocate(han_win(taper_len))

! equation of Hanning window: h(i)=0.5*(1-cos(2*pi*i)/N) ---> see Numerical Recipes, page: 547
! taking only the second (decaying) half 
do i=taper_len,1,-1
   han_win(taper_len-i+1) = 0.5 * (1 - cos(2 * pi * (i-1) / (2*taper_len)) )
enddo

! deallocate memory
if (allocated(h)) deallocate(h)



!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

CONTAINS


   SUBROUTINE CalcAfrM
   !-----------------------------------------------------------------
   !
   ! Description:
   ! 
   !   Calculates Rupture Area Ar given the moment magnitude Mw 
   !   from empirical relations of Wells & Coppersmith, 1994.
   !   Resultant Rupture Area (Ar) is [km^2]
   !   
   ! Dependencies:
   ! 
   !   None
   !
   ! Notes:
   !
   !   Four mechanism are accepted: strike-slip (ss), reverse-slip
   !   (rs), normal-slip (ns), other (al)
   !
   ! References:
   !
   !
   ! 
   ! Authors: W. Imperatori, B. Mena
   !
   ! Modified: January 2009 (v1.3)
   !

use fault_area
   
   implicit none

   ! local variables
   real(kind=r_single):: a,b

   !-----------------------------------------------------------------

   ! define parameters after mechanism type
   select case(trim(mech))
   case('ss')
      a=-3.42
      b=0.90
   case('rs') 
      a=-3.99
      b=0.98
   case('ns') 
      a=-2.87
      b=0.82
   case('al')
      a=-3.49
      b=0.91
   end select

   Ar=10**(a+b*Mw)    !compute rupture area

   END SUBROUTINE CalcAfrM

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

   SUBROUTINE fMw2MoN
   !----------------------------------------------------------------------------
   !
   ! Description:
   !
   !   Calculates the seismic moment Mo [in Nm] for a given moment magnitude Mw.
   !
   ! Dependencies:
   !
   !   None
   !
   ! Reference: Thorne & Lay, page 384
   !
   ! Authors: W. Imperatori, B. Mena
   !
   ! Modified: January 2009 (v1.3)
   !

   implicit none

   !----------------------------------------------------------------------------
   
   Mo=10**(1.5*Mw+9.05)

   END SUBROUTINE fMw2MoN

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

   SUBROUTINE CalcTrfromM 
   !-----------------------------------------------------------------------
   ! 
   ! Description:
   !
   !   Calculates Rise-Time [second] given the seismic moment (in [N*m], will
   !   be converted into [dyn*cm]).
   !
   ! Dependencies: 
   !
   !   None
   !
   ! References:
   !
   ! Somerville et al., 1999
   ! 
   ! Authors: W. Imperatori, B. Mena
   !
   ! Modified: January 2009 (v1.3)
   !

   implicit none

   ! dummy (local)
   real(kind=r_single) :: Mo2

   !-----------------------------------------------------------------------

   ! converting Mo [Nm] to [dyne-cm]; 1Nm = 10^7 dyne-cm
   Mo2=Mo*(10**7)
   
   ! estimate rise-time
   Tr=2.03*(Mo2**(1.0/3.0))/(10**9)

   END SUBROUTINE CalcTrfromM 

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

   SUBROUTINE stf_box
   !----------------------------------------------------------------
   !
   ! Description:
   !
   !   Calculates box-car type source-time-function in time domain,
   !   given amplitude Amp and rise-time Tr.
   !
   ! Dependencies:
   !
   !   None
   !
   ! Authors: W. Imperatori, B. Mena
   !
   ! Modified: January 2009 (v1.3)
   !
   
   implicit none

   ! index, dummy (locals)
   real(kind=r_single)   :: tof
   integer(kind=i_single):: i
      
   !-----------------------------------------------------------------   
      
   tof=ton+Tr

   do i=1,nint(ton/dt)+1
      stf(i)=0
   enddo

   do i=nint(ton/dt)+2,nint(tof/dt)+1
      stf(i)=Amp
   enddo

   do i=nint(tof/dt)+2,npts_stf
      stf(i)=0
   enddo
   
   END SUBROUTINE stf_box 

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

   SUBROUTINE stf_tri
   !-------------------------------------------------------------------------
   !
   ! Description:
   !
   !   Calculates triangular type source-time function in time domain given
   !   amplitude Amp and rise-time Tr.
   !
   ! Dependencies:
   ! 
   !   None
   !
   ! Authors: W. Imperatori, B. Mena
   !
   ! Modified: January 2009 (v1.3)
   !
 
   implicit none

   ! indexes, dummies (locals)
   integer(kind=i_single):: i
   real(kind=r_single)   :: damp,tmid,tof

   !-------------------------------------------------------------------------	

   tof=ton+Tr
   tmid=ton+(Tr/2)
   damp=(Amp*2)/(nint((Tr/2)/dt))

   do i=1,nint(ton/dt)+1
      stf(i)=0
   enddo

   stf(nint(ton/dt)+2)=damp

   do i=nint(ton/dt)+2,nint(tmid/dt)
      stf(i+1)=stf(i)+damp
   enddo

   do i=nint(tmid/dt)+1,nint(tof/dt)
      stf(i+1)=stf(i)-damp
   enddo

   do i=nint(tof/dt)+2,npts_stf
      stf(i)=0
   enddo

   END SUBROUTINE stf_tri

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

   SUBROUTINE stf_yof 
   !--------------------------------------------------------------
   !
   ! Description:
   !
   !   This subroutine computes Regularized Yoffe function.
   !
   ! Dependencies:
   !
   !   None
   !
   ! References:
   ! 
   !   Adapted from Tinti et al., 2005
   !
   ! Authors: W. Imperatori, B. Mena
   !
   ! Modified: January 2009 (v1.3)
   !

   implicit none
 
   ! indexes, dummies, coefficients (locals)
   integer(kind=i_single)                 :: i,j
   real(kind=r_single)                    :: c2,k,ts
   real(kind=r_single),dimension(npts_stf):: t,c1,c3,c4,c5,c6

   !--------------------------------------------------------------

   ! Rupture time (ton) assumed to be zero
   t(1)=0.
   do i=1,npts_stf-1
      t(i+1)=t(i)+dt
   enddo

   ! compute coefficients
   ts=0.3
   c2=(3.0/8.0)*pi*(Tr**2.0)
   k=2.0/(pi*Tr*ts**2.0)

   c1=((0.5*t+0.25*Tr)*sqrt(t*(Tr-t))+(t*Tr-(Tr**2.0))*asin(sqrt(t/Tr))     &
      -0.75*(Tr**2.0)*atan(sqrt((Tr-t)/t)))

   c3=((ts-t-0.5*Tr)*sqrt((t-ts)*(Tr-t+ts))+Tr*(2.0*Tr-2.0*t+2.0*ts)        &
      *(asin(sqrt((t-ts)/Tr)))+(3.0/2.0)*Tr**2.0*(atan(sqrt((Tr-t+ts)       &
      /(t-ts))))) 

   c4=((-1.0*ts+0.5*t+0.25*Tr)*sqrt((t-2.0*ts)*(Tr-t+2.0*ts))+Tr*((-1.0*Tr) &
      +t-2.0*ts)*asin(sqrt((t-2.0*ts)/Tr))-0.75*Tr**2.0*atan(sqrt((Tr-t     &
      +2.0*ts)/(t-2.0*ts)))) 

   c5=((pi/2.0)*Tr*(t-Tr))

   c6=((pi/2.0)*Tr*(2.0*ts-t+Tr))

   ! here below create Yoffe's stf
   if (Tr .GT. ts*2.0) then
      do j=1,npts_stf
         if (t(j) .LE. 0.0) then
            stf(j)=0.0
         else if (0.0 .LT. t(j) .AND. t(j) .LT. ts) then
            stf(j)=D*k*(c1(j)+c2)
            else if (ts .LE. t(j) .AND. t(j) .LT. 2.0*ts) then
                stf(j)=D*k*(c1(j)-c2+c3(j))
            else if (2.0*ts .LE. t(j) .AND. t(j) .LT. Tr) then
                stf(j)=D*k*(c1(j)+c3(j)+c4(j))
            else if (Tr .LE. t(j) .AND. t(j) .LT. Tr+ts) then
                stf(j)=D*k*(c5(j)+c3(j)+c4(j))
            else if (Tr+ts .LE. t(j) .AND. t(j) .LT. Tr+2.0*ts) then
                stf(j)=D*k*(c4(j)+c6(j))
            else if (Tr+2.0*ts .LE. t(j)) then
                stf(j)=0.0
         end if
      enddo
   endif


   if ((ts .LT. Tr) .AND. (Tr .LT. 2.0*ts)) then
      do j=1,npts_stf
         if (t(j) .LE. 0.0) then
                stf(j)=0.0
             else if (0.0 .LT. t(j) .AND. t(j) .LT. ts) then
                stf(j)=D*k*(c1(j)+c2)
             else if (ts .LE. t(j) .AND. t(j) .LT. Tr) then
                stf(j)=D*k*(c1(j)-c2+c3(j))
             else if (Tr .LE. t(j) .AND. t(j) .LT. 2.0*ts) then
                stf(j)=D*k*(c5(j)+c3(j)-c2)
             else if (2.0*ts .LE. t(j) .AND. t(j) .LT. Tr*ts) then
                stf(j)=D*k*(c5(j)+c3(j)+c4(j))
             else if (Tr*ts .LE. t(j) .AND. t(j) .LT. Tr+2.0*ts) then
                stf(j)=D*k*(c4(j)+c6(j))
             else if (Tr+2.0*ts .LE. t(j)) then
                stf(j)=0.0
         endif
      enddo
   endif

   END SUBROUTINE stf_yof 
   
!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-   

   SUBROUTINE stf_dreg
   !-------------------------------------------------------------------------
   ! 
   ! Description:
   !
   !   Computes Dreger et al., 2007 (SSA poster) type source-time-function. 
   !   Banu Mena, APRIL 2008, SDSU. Healing front + white noise -> Sept.08.
   !
   ! Dependencies:
   !
   !   None
   !
   ! Authors: W. Imperatori, B. Mena
   ! 
   ! Modified: January 2009 (v1.3)
   !

   implicit none

   ! dummies (locals)
   real(kind=r_single),dimension(npts_stf)        :: t,s,ds,randm
   real(kind=r_single),allocatable,dimension(:)   :: heal
   real(kind=r_single),parameter                  :: ksi=0.1
   real(kind=r_single)                            :: max_ss,dss,t1_heal,t2_heal,max_ran
   integer(kind=i_single)                         :: i,n_t1,n_t2
   integer(kind=i_single)                         :: seed,seed_size
   integer(kind=i_single),allocatable,dimension(:):: tmp_seed

   !-------------------------------------------------------------------------

   Tr=Tr*0.3        ! reviewed

   t(1)=0.0
   do i=1,npts_stf-1
      t(i+1)=t(i)+dt
   enddo

   s=(t**ksi)*exp(-t/Tr)

   !linear healing front between t1_heal and t2_heal
   t1_heal=3.0; t2_heal=3.2;
  
   n_t1=nint(t1_heal/dt)+1; n_t2=nint(t2_heal/dt)+1 

   dss=(s(n_t1))/(n_t2-n_t1)

   if(.not.allocated(heal)) allocate(heal(n_t2-n_t1))

   !linear healing function
   do i=1,(n_t2-n_t1)
      heal(i)=s(n_t1)-(dss*i)  
   enddo

   !function until healing
   do i=1,n_t1
      stf(i)=s(i)
   enddo

   !function during healing
   do i=1,(n_t2-n_t1)
      stf(n_t1+i)=heal(i)
   enddo

   !function after healing
   do i=1,(npts_stf-n_t2)
      stf(i+n_t2)=0.0
   enddo

   !adding white noise
   seed=321
   
   ! initialize the generator with a specified seed number
   call random_seed(size=seed_size)
   if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
   tmp_seed=seed
   call random_seed(put=tmp_seed)
  
   ! generate random numbers
   call random_number(randm)    
  
   max_ran=maxval(randm,dim=1)

   stf=stf+(randm/max_ran)*0.05

   ds(1)=stf(1)*dt
   do i=2,npts_stf 
      ds(i)=ds(i-1)+(stf(i-1)+stf(i))/2.*dt
   enddo

   max_ss=maxval(ds,dim=1)

   stf=(D/max_ss)*stf
   
   ! deallocate memory
   if(allocated(heal)) deallocate(heal)

   END SUBROUTINE stf_dreg

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

   SUBROUTINE stf_liu
   !-------------------------------------------------------------------
   !
   ! Description:
   !
   !   This subroutine computes Liu source-time function.
   !
   ! Dependencies:
   !
   !   None
   !
   ! Notes: 
   !
   !   this STF hasn't been widely tested!
   !
   ! Authors: W. Imperatori, B. Mena
   !
   ! Modified: January 2009 (v1.3)
   !

   implicit none

   ! locals
   real(kind=r_single),dimension(npts_stf)        :: t,s,ds,randm
   real(kind=r_single)                            :: t1,t2,cn,max_ss,max_ran
   integer(kind=i_single)                         :: seed,seed_size
   integer(kind=i_single),allocatable,dimension(:):: tmp_seed

   !-------------------------------------------------------------------


   t1=0.13*Tr; t2=Tr-t1
   cn=pi/(1.4*pi*t1+1.2*t1+0.3*pi*t2)

   t(1)=0.0
   do i=1,npts_stf-1
      t(i+1)=t(i)+dt
   enddo

   do i=1,npts_stf
      if (0.0 .LE. t(i) .AND. t(i) .LT. t1) then
         s(i)=cn*(0.7-0.7*cos(pi*t(i)/t1)+0.6*sin(0.5*pi*t(i)/t1))
      else if (t1 .LE. t(i) .AND. t(i) .LT. (2.0*t1)) then
         s(i)=cn*(1.0-0.7*cos(pi*t(i)/t1)+0.3*cos(pi*(t(i)-t1)/t2))
      else if ((2.0*t1) .LE. t(i) .AND. t(i) .LT. Tr) then
         s(i)=cn*(0.3+0.3*cos(pi*(t(i)-t1)/t2))
      else if (t(i) .GE. Tr) then
         s(i)=0.0
      end if
   enddo

   seed=123
   
   ! initialize the generator with a specified seed number
   call random_seed(size=seed_size)
   if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
   tmp_seed=seed
   call random_seed(put=tmp_seed)
   
   ! generate random numbers
   call random_number(randm)
   

   max_ran=maxval(randm,dim=1)

   stf=s+(randm/max_ran)*0.05


   ds(1)=s(1)*dt

   do i=2,npts_stf 
      ds(i)=ds(i-1)+(stf(i-1)+stf(i))/2.*dt
   enddo

   max_ss=maxval(ds,dim=1)

   stf=(D/max_ss)*stf

   END SUBROUTINE stf_liu

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

   SUBROUTINE stf_new
   !-------------------------------------------------------------------
   !
   ! Description:
   !
   !   This subroutine is for a new stf
   !
   ! Dependencies:
   !
   !   None
   !
   ! Notes: 
   !
   !   
   !
   ! Authors: K. Olsen 
   !
   ! Modified: March 2013 (v1.4.2)
   !
   ! Updated: May 2013 (v1.5.1)
   !   Add Tr_sca calculation
   !
   ! Updated: May 2013 (v1.5.2)
   !   Change Tr_sca calculation.
   !   Do not use Tr_sca from scattering.dat, just change the name to Tr_factor.
   !
   ! Updated: November 2013 (v1.5.3)
   !   Change Tr_factor calculation, depending on magnitude and dip.
   !
   ! Updated: September 2014 (v1.5.5.1)
   !   Add Tr_fac calculation for CENA events.
   !
   ! Updated: October 2014 (v1.5.5.1)
   !   Change Tr_fac calculation.
   !
   ! Updated: February 2015 (v1.5.5.3)
   !   Change Tr_fac calculation.
   !   Correct bc calculation for CENA events.
   !   Add tapering on stf.
   !
   ! Updated: May 2015 (v1.5.5.4)
   !   The update is applied for NGA-west2, only WNA and Japanese events,
   !   not for NGA-west1.
   !   Change alt and Tr_fac calculations for Tr.
   !   Add n and n1 parameters for stf computation.
   !
   ! Updated: June 2015 (v1.6)
   !   Merge NGA-west1 and NGA-west2 (v1.5.5.3 and v1.5.5.4).
   !   Move alt computation into srf_read subroutine.
   !  
   use vel_model, only: dip, alt
   use source_receiver, only: sr_rrup
   use flags, only: gs_flag, ngaw_flag

   implicit none

   ! locals
   real(kind=r_single),dimension(npts_stf)        :: t
   real(kind=r_single),dimension(6)               :: P
   real(kind=r_single)                            :: t1
   real(kind=r_single)                            :: Tr_fac
   real(kind=r_single)                            :: bc,rc
   real(kind=r_single)                            :: n, n1
   real(kind=r_single)                            :: distfac
   real(kind=r_single)                            :: mindist
 
   integer(kind=i_single)                         :: taper_sec
   integer(kind=i_single)                         :: taper_len
   real(kind=r_single),allocatable,dimension(:)   :: taper
   !-------------------------------------------------------------------

   ! NGA-west1 western NA and Japan
   if (ngaw_flag == 1) then
      if (gs_flag == 1 .or. gs_flag == 3) then
         rc=0.141-0.090*atan(Mw*1.60-9.55) ! change (02/11/2014)
         Tr_fac=rc*alt
         n=2.0
         n1=0.5
         !print*,'NGAW1 WNA-JPN: rc,dip,alt=',rc,dip,alt
      endif
   endif

   ! NGA-west2 western NA
   if (ngaw_flag == 2 .and. gs_flag == 1) then
      ! Distance dependence for Mw<6.2
      if (Mw < 6.2) then
         distfac=0.1
         mindist=40 ! [km]
         rc=( 0.1449-0.08715*atan(Mw*1.685- 9.45) )*alt
         bc=( 0.1226-0.05900*atan(Mw*1.750- 9.85) )*alt
         if (sr_rrup(station) < mindist) then
            Tr_fac=rc
         elseif (mindist<=sr_rrup(station).and.sr_rrup(station)<50.) then
            Tr_fac=rc + (sr_rrup(station)-mindist)*(bc-rc)*distfac
         elseif (sr_rrup(station) >= 50.) then
            Tr_fac=bc
         endif
      ! No distance dependence for M>=6.2
      elseif (Mw >= 6.2 .and. Mw < 7.5) then ! change (06/03/2015)
         Tr_fac=( 0.1449-0.08715*atan(Mw*1.685- 9.45) )*alt
      elseif (Mw >= 7.5) then
         Tr_fac=( 0.141-0.090*atan(Mw*1.60-9.55) )*alt
      endif

      ! NGA-west2 Rule for n and n1 (04/27/2015)
      if (Mw < 6.2) then
         n=2.75
         n1=0.6
      elseif (Mw >= 6.2 .and. Mw < 6.6) then
         n=-1.875*Mw+14.375
         n1=-0.25*Mw+2.15
      elseif (Mw >= 6.6) then
         n=2.0
         n1=0.5
      endif
      !print*,'NGAW2 WNA: dip, alt, dist= ',dip,alt,sr_rrup(station)
   endif

   ! NGA-west2 Japan (04/27/2015)
   if (ngaw_flag == 2 .and. gs_flag == 3) then
      rc=0.141-0.090*atan(Mw*1.60-9.55)
      Tr_fac=rc*alt
      n=2.75
      n1=0.5
      !print*,'Japanese-NGA-west2: Tr_fac*alt,n,n1= ',Tr_fac,n,n1
   endif

   ! Common NGA-west1 and NGA-west2 CENA
   if (gs_flag == 2) then ! add(09/26/2014)
      ! change rc computation (02/06/2015)
      if (Mw <= 5.95) then
         rc=( 0.1484-0.085*atan(Mw*1.40-7.250) )*alt
      elseif (Mw > 5.95 .and. Mw <= 7.5) then
         rc=( 0.2061-0.155*atan(Mw*1.40-7.250) )*alt
      elseif (Mw > 7.5) then
         rc=( 0.0870-0.060*atan(Mw*1.45-7.250) )*alt
      endif
      bc=( 0.141-0.090*atan(Mw*1.60-9.55) )*alt ! corrected bc (02/13/2015)

      ! change (10/14/2014)
      if (sr_rrup(station) <= 100.0) then
         Tr_fac=rc
      elseif (sr_rrup(station)>100. .and. sr_rrup(station)<=120.) then
         Tr_fac=rc+(bc-rc)*(sr_rrup(station)-100.)*0.05
      elseif (sr_rrup(station) > 120.) then
         Tr_fac=bc
      endif
      n=2.0
      n1=0.5
      !print*,'CENA: dip, alt, dist= ',dip,alt,sr_rrup(station)
   endif

   t(1)=0.0
   do i=1,npts_stf-1
      t(i+1)=t(i)+dt
   enddo

   t1=Tr_fac*Tr

   ! stf=(t**0.5) / (1 + (t/(t1))**2.)
   stf=(t**n1) / (1 + (t/(t1))**n) ! change to use n and n1 (04/27/2015)

   ! Add tapering on stf (02/15/2015)
   ! compute tapering length as a function of magnitude (Tr)
   taper_sec = nint(Tr/dt)
   taper_len = nint(1/dt)

   ! compute taper
   if (.not.allocated(taper)) allocate(taper(taper_len))
   do i=1,taper_len
      taper(i) = 0.5 * (1 + (cos(pi*i/taper_len)) )
   enddo

   ! apply tapering from Tr to Tr+1sec
   stf(taper_sec+1 : taper_sec+taper_len)=                                &
        stf(taper_sec+1 : taper_sec+taper_len) * taper

   ! zeroing later part
   stf(taper_sec+taper_len+1:npts_stf) = 0

   if (allocated(taper)) deallocate(taper)

   END SUBROUTINE stf_new

END SUBROUTINE set_stf

!===================================================================================================   

SUBROUTINE MomRateFunct(dtr)
!--------------------------------------------------------------------------------
! 
! Description:
!
!   Reads a user-defines moment-rate function in the format of Rob Graves.
!   The resulting stf is resampled in stf_npts (i.e. len=15 sec, npts=stf_npts
!   defined in the calling subroutine)
!
! Dependencies:
!
!   Subroutines inac, spline_interp
!
! WARNING: 
!
!   This subroutine has to be still verified!
!
! Authors: R. Graves, W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use interfaces, only: spline_interp; use stf_data

implicit none

! indexes, dummies (locals)
real(kind=r_single),intent(in)                 :: dtr
real(kind=r_single)                            :: r_dum,dtt,max_ton,t1,t2,sum_sdd
integer(kind=i_single)                         :: i_dum,Npts,nt_loc,NT,Tmax,Ntm,ind,i
integer(kind=i_single)                         :: pp,nn,npo,j,indic
character(len=80)                              :: s_dum
real(kind=r_single),allocatable,dimension(:)   :: dt,sv,sd,max_sd,tony,svp,ar 
real(kind=r_single),allocatable,dimension(:)   :: mtime,sdd,mom,svf,svf_new,z
integer(kind=i_single),allocatable,dimension(:):: b,Np
!-------------------------------------------------------------------------------

!open file with stf and read values
open(1,file=trim(stf_name),status='old')
read(1,*); read(1,*); read(1,*); read(1,*)
read(1,*) s_dum,Npts

if(.not.allocated(dt)) allocate(dt(Npts),tony(Npts),max_sd(Npts),ar(Npts))

!loop over stf' # of points
do i = 1, Npts
   read(1,*) r_dum,r_dum,r_dum,i_dum,i_dum,ar(i),tony(i),dt(i)
   read(1,*) i_dum,r_dum,nt_loc,r_dum,i_dum,r_dum,i_dum

   if(.not.allocated(sv)) allocate(sv(nt_loc),sd(nt_loc))

   if (nt_loc /= 0) then
      read(1,*) (sv(j), j=1,nt_loc)
      call inac   !computes displacement
      max_sd(i)=maxval(sd)
   endif

   if(allocated(sv)) deallocate(sv,sd)

enddo

close(1)

! Number of samples in local slip function
Nt=100
dtt=0.1

max_ton=maxval(tony)

Tmax=nint(max_ton+(Nt+1)*dt(1))
Ntm=nint(Tmax/dtt)

if(.not.allocated(mtime)) allocate(mtime(Ntm),Np(Ntm),svp(Ntm),mom(Ntm),svf(Ntm))
if(.not.allocated(b)) allocate(b(Npts),sdd(Npts))

mtime(1)=0.0

!new time vector
do i=1, Ntm-1
   mtime(i+1)=mtime(i)+dtt
enddo

Np(1)=0

do pp=2,Ntm
   t1=mtime(pp-1)
   t2=mtime(pp)
   ind=0

   do i=1,Npts
      if (tony(i) .ge. t1 .and. tony(i) .lt. t2) then
         ind=ind+1
         b(ind)=i
      endif
   enddo
   
   Np(pp)=ind

   if (Np(pp) .eq. 0) then
      
      svf(pp)=0.0
      mom(pp)=0.0
   
   else 
      do j=1,ind
         indic=b(j)
         sdd(j)=max_sd(indic)
      enddo
       
      sum_sdd=sum(sdd)
      svf(1)=0.0
      mom(1)=0.0
      svf(pp)=sum_sdd/Np(pp)
      
!!	mom(pp)=3.3*1e10*sum_sdd*(1/dtt)*1e-2*ar(1)*1e-4
!! ----------------------------------------------------------
!!  the scaling below, w/o rigidity and with 1e-12 simply
!!  creates mom with about the same amplitudes as svf;
!!  since the spectral amplitudes will be scaled later anyway
!!  this does not make any difference for the final results
!!  to get the correctly scaled moment-rate function, simply
!!  scale by a factor 3.3*1e16
! ----------------------------------------------------------

      mom(pp)=sum_sdd*(1/dtt)*ar(1)*1e-12
   endif

enddo
 
! Add 1 sec zeros to the beginning of the source time function 
!nn=NINT(1/dtt)

!npo=pp+nn-1    !npts with initial zero padding

if(.not.allocated(svf_new)) allocate(svf_new(Ntm))

!do i=1,nn
!   z(i)=0.0
!enddo   
      
!do i=1,nn
!   svf_new(i)=z(i)
!enddo  

!do i=nn+1,npo
do i=1,Ntm
!! ------------------------------------------------------------------------
!! two options to construct the "moment-rate function":
!!
!! (1) using svf means the summed slip-velocity at each instant in time
!!     hence resulting in a rather jagged function
!! (2) using mom means in fact moment (i.e. scaled by area and rise time)
!!     hence a somewhat smoother convolution operator
!! ------------------------------------------------------------------------
!!	svf_new(I) = svf(I-nn)
   svf_new(i) = mom(i)
enddo

!mmtime(1)=0.0

!do i=1,npo-1
!   mmtime(i+1)=mmtime(i)+dtt      !old time vector
!enddo

! resampling the moment rate function svf_new
npts_stf=nint((Ntm-1)*dtt/dtr)+1

total=dtr*(npts_stf-1)

if(allocated(t_stf)) deallocate(t_stf,stf)
if(.not.allocated(stf)) allocate(stf(npts_stf),t_stf(npts_stf))

call spline_interp(svf_new,total,npts_stf,Ntm,stf)

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

CONTAINS


   SUBROUTINE inac              
   !------------------------------------------------------------
   ! 
   ! Description:
   !
   !   Inac means 'integration of acceleration', but here is 
   !   used with velocity to obtain displacement
   !
   ! Dependencies:
   !
   !   None
   !
   ! Authors: R. Graves, W. Imperatori
   !
   ! Modified: January 2009 (v1.3)
   !
   
   implicit none

   ! counter (local)
   integer(kind=i_single):: j

   !-------------------------------------------------------------

   sd(1)=sv(1)*dt(i)

   do j=2,nt_loc 
      sd(j)=sd(j-1)+(sv(j-1)+sv(j))/2.*dt(i)
   enddo

   END SUBROUTINE inac


END SUBROUTINE MomRateFunct

!===================================================================================================   

SUBROUTINE srf_read
!--------------------------------------------------------------------------------
! 
! Description:
!
!   Reads a srf file
!
! Dependencies:
!
!   none
!
!
! Author: R. Takedatsu 
!
! Modified: March 2013 (v1.4.2) 
!
! Updated: February 2015 (v1.5.5.3)
!    Add str_fac computation for eastern CA events.
!
! Updated: June 2015 (v1.6)
!    Move alt computation from source.f90 and composition.f90
!    into this subroutine.
!
use stf_data
use geometry, only: n_lay; use vel_model
use scattering, only: str_fac
use flags, only: gs_flag, ngaw_flag

implicit none

real(kind=r_single)                            :: r_dum
integer(kind=i_single)                         :: i_dum,nt_loc,i,j
integer(kind=i_single)                         :: nstk,ndip
character(len=80)                              :: s_dum
real(kind=r_single),allocatable,dimension(:)   :: dt,sv,tony
real(kind=r_single)                            :: ave_slip,mw_test
real(kind=r_single)                            :: arcu
real(kind=r_single)                            :: dtop
real(kind=r_single),allocatable,dimension(:)   :: slip1
real(kind=r_single),allocatable,dimension(:)   :: slip1cu,mu_tmp,mo_tmp
!-------------------------------------------------------------------------------

!open file with srf and read values

open(1,file=trim(srf_name),status='old')
read(1,*) ! read version e.g. 1.0
read(1,*) ! read header block line1 e.g. PLANE 1
read(1,*) ! read header block line2, elon in x,elat in y,nstk,ndp,len,wid
!read(1,*) r_dum,dip,r_dum,r_dum,r_dum ! read header block line3, stk,dip,dtop,shyp,dhyp
read(1,*) r_dum,dip,dtop,r_dum,r_dum ! read header block line3, stk,dip,dtop,shyp,dhyp
read(1,*) s_dum,nsub ! read data block e.g. POINTS, np

!compute str_fac for eastern CA events
if (str_fac == 0.0) then
   str_fac = (225.0 + dtop*16.0)
   str_fac = str_fac * 10**6
endif
if (str_fac < 10**6) then
   print*,'check str_fac should be in bars'
else
   print*,'str_fac = ',str_fac
endif

!compute alt
if (ngaw_flag == 1 .or. gs_flag == 2) then
   if (dip <= 55.) then
      alt=0.82
   elseif (dip > 55. .and. dip < 70.) then
      alt=0.16+0.012*dip
   elseif (dip >= 70.) then
      alt=1.
   endif
elseif (ngaw_flag == 2) then
   if (gs_flag == 1 .or. gs_flag == 3) then
      if (dip <= 45.) then
         alt=0.785
      elseif (dip > 45. .and. dip < 70.) then
         alt=0.0086*dip+0.3980
      elseif (dip >= 70.) then
         alt=1.
      endif
   endif
endif

if(.not.allocated(dt)) then
   allocate(dt(nsub),tony(nsub),slip1(nsub))
endif
if(.not.allocated(fx)) allocate(fx(nsub),fy(nsub),fz(nsub),areas(nsub))

!loop over srf' # of points
do i = 1, nsub
   read(1,*) fx(i),fy(i),fz(i),i_dum,i_dum,areas(i),tony(i),dt(i)
   read(1,*) i_dum,slip1(i),nt_loc,r_dum,i_dum,r_dum,i_dum

   if(.not.allocated(sv)) allocate(sv(nt_loc))

   if (nt_loc /= 0) then
      read(1,*) (sv(j), j=1,nt_loc)
   endif

   if(allocated(sv)) deallocate(sv)

enddo

sum_area=sum(areas)
ave_slip=sum(slip1)/nsub

close(1)


! compute mu (shear modulus) and moment
! allocate arrays
if(.not.allocated(mo_tmp)) allocate(mu_tmp(nsub),slip1cu(nsub),mo_tmp(nsub))
if(.not.allocated(mo_ratio)) allocate(mo_ratio(nsub),vscu(nsub),rhcu(nsub))

do i = 1, nsub
   mu_tmp(i)=0.0
   slip1cu(i)=0.0
   mo_tmp(i)=0.0

   do j = 2, n_lay
      if ( (fz(i) >= depth(j-1)) .and. (fz(i) < depth(j)) ) then
         vscu(i)=vs(j-1)*1000      ! km/s -> m/s
         rhcu(i)=rh(j-1)*1000      ! g/cm3 -> kg/m3
      endif

      if (slip1(i) > 0.0) then
         mu_tmp(i)=rhcu(i)*vscu(i)**2        ! [(kg*m2)/(m3*s2) = kg/(m*s2)]
         slip1cu(i)=slip1(i)/100             ! cm -> m
         arcu=areas(i)/10000                 ! cm2 -> m2
         mo_tmp(i)=mu_tmp(i)*slip1cu(i)*arcu ! [kg*m2/s2 = Nm]
      endif
   enddo

enddo

! calculation of Mw from equation of Mo[Nm]=10**(1.5*Mw+9.05)
total_Mo=sum(mo_tmp)
mw_test=(log10(total_Mo) - 9.05)/1.5

! for acc_spec
sum_area=sum_area/10000  ! cm2 -> m2
do i = 1, nsub
   if (slip1(i) == 0.0) then
      mo_ratio(i)=0.0
   else
      mo_ratio(i)=slip1cu(i)*mu_tmp(i)*sum_area/total_Mo
   endif
enddo

! change units for acc_spec
sum_area=sum_area*10000. ! m2 -> cm2
vscu=vscu*100.           ! m/s -> c/m
rhcu=rhcu/1000.          ! kg/m3 -> g/cm3

! deallocate memory
if(allocated(mo_tmp)) deallocate(mu_tmp,slip1cu,mo_tmp)
!

END SUBROUTINE srf_read
