SUBROUTINE conv_stf_scat(station)
!--------------------------------------------------------------------------------
!
! Description:
!
!   Convolves scattered wavelets with selected moment-rate function for each 
!   component of motion (x, y and z). Convolution is performed in time domain.
!
! Dependencies:
!
!   Subroutine convolv.
!
! Notes:
!
!   For extended-source case, Irikura's Empirical Green's Functions (EGF) approach
!   is adopted.
!
! References:
!
!
!
! Authors: W. Imperatori, B. Mena
!
! Modified: January 2009 (v1.3)
!
! Updated: December 2012 (v1.4.1)
!   change taper_len.
!
! Updated: April 2014 (v1.5.5)
!    tmp_npts and tmp_lf_len are in the module tmp_para.

use constants; use def_kind; use flags
use interfaces, only: convlv; use source_receiver
use scattering, only: npts; use stf_data; use waveform
use tmp_para

implicit none
 
! actual station number 
integer(kind=i_single),intent(in)             :: station
! time-step of stf
real(kind=r_single)                           :: dt
! counters, dummies
integer(kind=i_single)                        :: i,j,npts_conv,npts_out,taper_len
! temporary array for extended-source case  
real(kind=r_single),allocatable,dimension(:,:):: stf_ext
real(kind=r_single),allocatable,dimension(:)  :: stf_sum
! adjusted dt
real(kind=r_single)                           :: tmp_dt
!-------------------------------------------------------------------------------

! allocate array for temporary stf (effectively used only for extended-fault case)
allocate(stf_sum(npts_stf))

stf_sum=stf   !assign stf

! EGF approach if extended-source is selected
if (ext_flag == 1) then

   ! allocate array for stf for each subfault
   allocate(stf_ext(n_cell,npts_stf))

   ! weighting due to distance of the subfault to the fault with respect to the hypocenter
   forall(i=1:n_cell, j=1:npts_stf)
      stf_ext(i,j) = ( sr_hypo(station) / sr_cell(i,station) ) * stf(j)
   end forall

   ! summation of the contribution from each subfault
   stf_sum=sum(stf_ext,dim=1)     !column-wise summation (i.e. 1 number for each column)

   deallocate(stf_ext)   !deallocate memory

endif

! compute time-step of both time-series (scattered and stf)
dt = lf_len/(npts-1)
tmp_dt = tmp_lf_len/(tmp_npts-1)

! compute taper length
!!!taper_len=nint(npts*0.5)
if (npts == tmp_npts) taper_len=nint(10.0/dt)
if (npts .gt. tmp_npts) taper_len=nint(10.0/tmp_dt)

! tapering only last part (~20%) of scattered time-series
if (npts == tmp_npts) then
   forall (i=1:taper_len,j=1:3)
      scattgram(npts-taper_len+i,j)=scattgram(npts-taper_len+i,j)*han_win(i)
   end forall
endif
if (npts .gt. tmp_npts) then
   forall (i=1:taper_len,j=1:3)
      scattgram(tmp_npts-taper_len+i,j)                                &
                =scattgram(tmp_npts-taper_len+i,j)*han_win(i)
   end forall
endif

! perform convolution for each component
do i=1,3
   if (npts == tmp_npts) conv_seis(:,i)=convlv(scattgram(:,i),stf_sum)*dt
   if (npts .gt. tmp_npts) conv_seis(:,i)=convlv(scattgram(:,i),stf_sum)*tmp_dt
enddo

! deallocate temporary arrays
deallocate(stf_sum)


END SUBROUTINE conv_stf_scat
 
