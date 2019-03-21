SUBROUTINE poly_interp(x1a,x2a,ya,x1,x2,y)
!------------------------------------------------------------------------
!
! Description:
!
!   Returns an interpolated value using polynomial interpolation (2D). If
!   input function is [2x2], bilinear interpolation is performed
!   
!
! Dependencies:
!
!   Subroutine polint
!
! References:
!
!   Numerical Recipes
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind; use interfaces, only: polint

implicit none

! arrays for x,y coordinates	
real(kind=r_single),dimension(:),intent(in)  :: x1a,x2a
! values of function at (x,y)
real(kind=r_single),dimension(:,:),intent(in):: ya
! coordinates of interpolated point
real(kind=r_single),intent(in)               :: x1,x2
! interpolated point
real(kind=r_single),intent(out)              :: y
! counter, vector length 
integer(kind=i_single)                       :: i,m
! temporary arrays 
real(kind=r_single),dimension(size(x1a))     :: ymtmp
real(kind=r_single),dimension(size(x2a))     :: yntmp

!-----------------------------------------------------------------------

! check for data 
if ( size(x1a) /= size(ya,1) ) call error_handling(7,'null','POLY_INTERP (interpolation.f90)')
if ( size(x2a) /= size(ya,2) ) call error_handling(7,'null','POLY_INTERP (interpolation.f90)')

m=size(x1a)

do i=1,m
   yntmp = ya(i,:)
   call polint(x2a,yntmp,x2,ymtmp(i))
enddo

call polint(x1a,ymtmp,x1,y)   


END SUBROUTINE poly_interp

!===================================================================================================

SUBROUTINE polint(xa,ya,x,y)
!----------------------------------------------------------------
!
! Description:
!
!   Returns value through polynomial interpolation (1D), that means
!   the returned value is y=P(x), where P is polynomial of degree
!   N-1 (N is the number of points in xa)
!
! Dependencies:
!
!   None
!
! References:
!
!   Numerical recipes
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

! arrays for x coordinate and y values
real(kind=r_single),dimension(:),intent(in):: xa,ya
! coordinate of interpolated value
real(kind=r_single),intent(in)             :: x
! interpolated value and error estimation
real(kind=r_single),intent(out)            :: y
! counter and dummies (locals)
real(kind=r_single)                        :: dy
integer(kind=i_single)                     :: m,n,ns
real(kind=r_single), dimension(size(xa))   :: c,d,den,ho
integer(kind=i_single),dimension(1)        :: imin

!----------------------------------------------------------------

! check for data consistency 
if ( size(xa) /= size(ya) ) call error_handling(7,'null','POLINT (interpolation.f90)')

n=size(xa)   

! initialize tableau
c=ya
d=ya

ho=xa-x

! find index closest to table entry
imin=minloc(abs(x-xa))
ns=imin(1)   

! initial approximation
y=ya(ns)
ns=ns-1

do m=1,n-1
   
   den(1:n-m)=ho(1:n-m)-ho(1+m:n)      
   
   ! abort if two xa are identical
   if (any(den(1:n-m) == 0.0))  call error_handling(4,'null','POLINT (interpolation.f90)')

   den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
   d(1:n-m)=ho(1+m:n)*den(1:n-m)
   c(1:n-m)=ho(1:n-m)*den(1:n-m)

   if (2*ns < n-m) then
      dy=c(ns+1)
   else
          dy=d(ns)
          ns=ns-1
   endif

   y=y+dy

enddo

END SUBROUTINE polint

!===================================================================================================

SUBROUTINE sampling_check
!----------------------------------------------------------------------
!
! Description:
!
!   Compute number of points (npts) for the interpolated time-series.
!   This value is to guarantee a Nyquist frequency >= the maximum 
!   frequency used during coda wave computation AND a npts that is a
!   2**N value (necessary for FFT/IFFT subroutines)
!
! Dependencies:
!
!   None
!
! Notes:
!
!   For the "LF-HF merging" modality: to assure same time-series length
!   (in time), interpolation is performed here.
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
   
use def_kind; use flags; use scattering, only: fmax,npts
use source_receiver, only: n_stat; use waveform

implicit none   

! local varaibles
real(kind=r_single):: dt,exponent
 
   
!---------------------------------------------------------------------

if (modality_flag == 0) then

   ! compute HF time-sampling
   dt = hf_len / (hf_npts - 1)
   
   ! recompute number of points to closest 2**N number 
   npts = hf_npts   !related number of points

   ! check if this is a 2**N value (if not, adjust npts)
   exponent = log(real(npts)) / log(2.0)

   if (exponent /= nint(exponent))  npts = 2**(ceiling(exponent)) 
   
   if (lf_len /= hf_len .or. hf_npts /= lf_npts)  then
      call error_handling(4,'null','SAMPLING_CHECK (interpolation.f90)')
   endif
   
else if (modality_flag /= 0) then

   ! assume maximum frequency for coda waves as Nyquist frequency
   dt = 1 / (2 * fmax)   !new time-step 
     
   npts = nint(lf_len / dt) + 1   !related number of points

   ! check if this is a 2**N value (if not, adjust npts)
   exponent = log(real(npts)) / log(2.0)

   if (exponent /= nint(exponent))  npts = 2**(ceiling(exponent))
endif

END SUBROUTINE sampling_check

!===================================================================================================

SUBROUTINE spline_interp(in_seis,time_len,interp_npts,orig_npts,out_seis)
!-------------------------------------------------------------------------------------
!
! Description:
!
! Performs spline interpolation on given time-series using a different sampling rate.
! Time series is assumed to be given in an 1 x npts array. 
!
! Dependencies:
!
!   Subroutines spline, splint
!
! References:
!
!   Numerical Recipes
!
! Authors: W. Imperatori, M. Mai
!
! Modified: January 2009 (v1.3)
!
! Updated: February 2019 (v2.0)
!   Change dt_interp and dt_orig computations.
!
use def_kind
use tmp_para;
use waveform, only: lf_dt,v_npts

implicit none 

! number of points of the new (interpolated) and original time-series
integer(kind=i_single),intent(in)           :: interp_npts,orig_npts
! time-length of the time-series
real(kind=r_single),intent(in)              :: time_len
! input time-series
real(kind=r_single),dimension(:),intent(in) :: in_seis       !assumed-shape vec
! output time-series
real(kind=r_single),dimension(:),intent(out):: out_seis     !assumed-shape vec
! local arrays for time vector and second derivative
real(kind=r_single),allocatable,dimension(:):: t_vec,sec_der
! time-vector, time-steps of interpolated and original time-series
real(kind=r_single)                         :: t_interp,dt_interp,dt_orig
! counter
integer(kind=i_single)                      :: j,n
! parameters used for spline interpolation
real(kind=r_single),parameter               :: yp1=2.0e20, ypn=2.0e20

!------------------------------------------------------------------------------------

! time-step of interpolated time-series
!dt_interp = time_len / (interp_npts-1)
n=ceiling(time_len/tmp_lf_len)
dt_interp = tmp_lf_len*n / (v_npts - 1)

! time-step of original time-series
!dt_orig = time_len / (orig_npts-1)
! dt_orig is from the actual dt from LF = lf_dt
dt_orig = lf_dt

if(.not.allocated(t_vec)) allocate(t_vec(orig_npts),sec_der(orig_npts))                 

!create original time vector
do j = 1,orig_npts
   t_vec(j) = dt_orig*(j-1)
enddo

! returns second derivative terms
call spline(t_vec,in_seis,orig_npts,yp1,ypn,sec_der)

! loop over number of points for interpolated time-series
do j = 1,interp_npts

   t_interp = dt_interp*(j-1)     !new time vector
 
   ! spline interpolation 
   call splint(t_vec,in_seis,sec_der,orig_npts,t_interp,out_seis(j))
    
enddo

! deallocate temporary arrays   
if(allocated(t_vec)) deallocate(t_vec,sec_der)

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

CONTAINS


   SUBROUTINE spline(x,y,n,yp1,ypn,y2)
   !---------------------------------------------------------------
   !
   ! Description:
   !
   !   Computes second derivative terms
   !
   ! Dependencies:
   !
   !   None
   !
   ! References:
   !
   !   Numerical Receipes
   ! 
   ! Authors: W. Imperatori, M. Mai
   !
   ! Modified: January 2009 (v1.3)
   !

   implicit none

   ! x-values and y-values of input time-series
   real(kind=r_single),dimension(:),intent(in):: x,y
   ! second derivative terms
   real(kind=r_single),dimension(:),intent(out):: y2
   ! original time-series number or points
   integer(kind=i_single),intent(in):: n
   ! parameters for interpolation
   real(kind=r_single),intent(in):: yp1,ypn
   ! local counters, dummies and indexes
   real(kind=r_single),dimension(30000):: u
   real(kind=r_single):: p,qn,sig,un
   integer(kind=i_single):: i,k
   
   !--------------------------------------------------------------

   if (yp1.gt.1e20) then
      y2(1)=0.
      u(1)=0.
   else
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
   endif

   do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))    &
           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
   enddo
      
   if (ypn.gt.1e20) then
      qn=0.
      un=0.
   else
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
   endif

   y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

   do k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
   enddo
     
   END SUBROUTINE spline

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

   SUBROUTINE splint(xa,ya,y2a,n,x,y)
   !------------------------------------------------------------
   !
   ! Description:
   !
   !   Performs interpolation using spline
   !
   ! Dependencies:
   !
   !   None
   !
   ! References:
   !
   !   Numerical Receipes
   ! 
   ! Authors: W. Imperatori, M. Mai
   !
   ! Modified: January 2009 (v1.3)
   !
  
   implicit none

   ! x-values and y-values of input time-series, second derivative terms
   real(kind=r_single),dimension(:),intent(in):: xa,ya,y2a
   ! original time-series number of points
   integer(kind=i_single),intent(in):: n
   ! x-value of new (interpolated) time-series
   real(kind=r_single),intent(in):: x
   ! y-value of new (interpolated) time-series
   real(kind=r_single),intent(out):: y
   ! local counter and dummies
   integer(kind=i_single):: k,khi,klo
   real(kind=r_single):: a,b,h

   !------------------------------------------------------------

   klo=1
   khi=n

1  if (khi-klo.gt.1) then
      k=(khi+klo)/2
   
      if (xa(k).gt.x) then
         khi=k
      else
         klo=k
      endif
      goto 1
   endif

   h=xa(khi)-xa(klo)

   ! aborting execution
   if (h.eq.0.) call error_handling(4,'null','SPLINT (interpolation.f90)')
   
   
   a=(xa(khi)-x)/h
   b=(x-xa(klo))/h
   y=a*ya(klo)+b*ya(khi)+((a**3-a) * y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      
   END SUBROUTINE splint

     
END SUBROUTINE spline_interp 

!===================================================================================================
