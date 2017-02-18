FUNCTION delay(station)
!----------------------------------------------------------------------------
!
! Description:
!
!   Compute delays in S-wave arrival times for each station between deterministic 
!   LF component and stochastic HF component (EXSIM)
!
! Dependencies:
!
!   None
!
! Notes:
!
!   HF first arrival is computed through a running-window algorithm. If delay
!   value is positive, deterministic S-wave onset is located later in time.
!
! Author: W. Imperatori 
!
! Modified: January 2009 (v1.3)
!

use def_kind; use interfaces, only: peaker; use source_receiver; use waveform 

implicit none

! station number
integer(kind=i_single),intent(in):: station
! function
real(kind=r_single)              :: delay

!---------------------------------------------------------------------------

! evaluate arrival time difference
delay = time_s(station) - peaker(hf_seis(:,1),hf_len/(hf_npts-1))


END FUNCTION delay

!===================================================================================================

FUNCTION derive(fun,dt)
!----------------------------------------------------------------------
! 
! Description:
!
!   Compute derivative of a function using lower terms of Taylor series.
!
! Dependencies:
!
!   None.
!
! Note:
!
!   Output time-series will have one point less than the input time-series.
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

! input function
real(kind=r_single),dimension(:),intent(in):: fun
! time-step of input function
real(kind=r_single),intent(in)             :: dt
! derivative (output)
real(kind=r_single),dimension(size(fun)-1) :: derive
! number of points of input time-series
integer(kind=i_single)                     :: npts,i

!----------------------------------------------------------------------

! number of points
npts = size(fun)

do i=1,npts-1
  
   derive(i) = (fun(i+1) - fun(i)) / dt
   
enddo   

END FUNCTION derive

!===================================================================================================

FUNCTION peaker(fun,dt)
!-----------------------------------------------------------------------------
!
! Description:
!
!   Evaluate first break arrival time using a running window algorithm
!
! Dependencies:
!
!   Functions trapz and derive
!
! Note:
!
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind; use interfaces, only: derive, trapz

implicit none

! input time-series
real(kind=r_single),dimension(:),intent(in) :: fun
! time-step of input time-series
real(kind=r_single),intent(in)              :: dt
! output value (representing first break arrival time)
real(kind=r_single)                         :: peaker
! widowing functions
real(kind=r_single)                         :: W2,W1
! windowing functions length, index
integer(kind=i_single)                      :: W2_npts,W1_npts,i
! local time-series
real(kind=r_single),allocatable,dimension(:):: derivative,f



! running window length (sec) and its number of points
W2 = 0.05; W2_npts = W2 / dt +1;

! total window length (sec) over which search is performed and its number of points
W1 = 3.0; W1_npts = W1 / dt +1;

! allocate memory
allocate(f(W1_npts))

! compute energy-ratio function
do i=1,W1_npts
  f(i) = trapz(fun(i:i-1+W2_npts)**2,dt) / trapz(fun(1:i-1+W2_npts)**2,dt) 
enddo

!open(1,file='tau.txt')
!do i=1,W1_npts
!   write(1,*) (i-1)*dt,f(i)
!enddo
!close(1)

! allocate memory
allocate(derivative(W1_npts-1))

! compute derivative
derivative = derive(f,dt)

!open(1,file='deriv.txt')
!do i=1,W1_npts-1
!   write(1,*) (i-1)*dt,derivative(i)
!enddo
!close(1)

! extract index at which peak value occurs
peaker = (maxloc(derivative,dim=1) - 1)*dt

!do i=1,W1_npts-1
!   if (maxval(derivative) == derivative(i)) then
!      !print*, (i-1)*dt
!      peaker = (i-1)*dt
!      exit
!   endif
!enddo   
      
! free memory
deallocate(f,derivative)

END FUNCTION peaker

!===================================================================================================

SUBROUTINE time_distance
!-----------------------------------------------------------------------------
!
! Description:
!
!   Reads P and S travel-times at a given number of locations, computed 
!   using J.Hole's 3D-raytracing program. The input files are binary. It 
!   computes also S-R distances, eventually for each subfault (Empirical
!   Green's Function approach).
!
! Dependencies:
!
!   None. All variables passed by modules.
!
! Notes: 
!
!   Very simple way out to allow for input coordinates given in km has been
!   taken; if you enter positions that DO NOT agree with a sample point in 
!   your grid, the routine will take the P/S-arrival times from the nearest 
!   grid-point (since no high precision is required).
!
! Authors: W. Imperatori, M. Mai
!
! Modified: January 2009 (v1.3)
!
! Updated: September 2014 (v1.5.5.1)
!   Add sr_rrup calculation, finding Rrup distance for each station.

use def_kind; use flags; use geometry; use interfaces, only: poly_interp
use source_receiver

implicit none

! indexes, num. of records
integer(kind=i_single)                        :: j,k,i,nrec
! arrival times from raytracing
real(kind=r_single),allocatable,dimension(:,:):: tp,ts
! grid-indexes for neighbouring points
integer(kind=i_single),dimension(2)           :: x_pos,y_pos

! ----------------------------------------------------------------------------

! allocate arrays for arrival times and S-R distances
allocate(tp(nx,ny),ts(nx,ny),time_p(n_stat),time_s(n_stat))

! allocate arrays for hypocentral distances
allocate(sr_hypo(n_stat))
! allocate array for Rrup distances
allocate(sr_rrup(n_stat))

! compute distances from receivers to source
do k=1,n_stat
   ! hypocenter-station distance
   sr_hypo(k) = sqrt( (xp(k) - hyp_x)**2 + (yp(k) - hyp_y)**2 + (hyp_z)**2 )
enddo

! in case of extended source, calculate distances between each source-receiver couple
if (ext_flag == 1) then
   ! allocate array for cells-receivers distances
   allocate(sr_cell(n_cell,n_stat))
   ! compute distances   
   forall (i=1:n_cell,k=1:n_stat)
      sr_cell(i,k) = sqrt( (xp(k) - x_cell(i))**2 + (yp(k) - y_cell(i))**2 + (z_cell(i))**2 ) 
   end forall

   ! compute Rrup
   do k=1,n_stat
      sr_rrup(k)=999999.
      do i=1,n_cell
         if (sr_rrup(k)>sr_cell(i,k)) sr_rrup(k)=sr_cell(i,k)
      enddo
   enddo

endif   

! open binary files that contain the travel-times
!open(1,file='time3d_P.out',access='direct',recl=nx)
!open(2,file='time3d_S.out',access='direct',recl=nx)
open(1,file='time3d_P.out',access='direct',recl=4)
open(2,file='time3d_S.out',access='direct',recl=4)

!open(3,file='tp.out',status='unknown')
!open(4,file='ts.out',status='unknown')


! read travel times into a potentially large ny x nx array
! using only the information for nz=1, i.e. the surface layer
nrec=0
!do j=1,ny
!   nrec=nrec+1
!   read(1,rec=nrec) (tp(i,j),i=1,nx)
!   read(2,rec=nrec) (ts(i,j),i=1,nx)
!enddo
do j=1,ny
   do i=1,nx
      nrec=nrec+1
      read(1,rec=nrec) tp(i,j)
      read(2,rec=nrec) ts(i,j)
!      print*,'nrec= ',nrec
!      write(3,*) tp(i,j)
!      write(4,*) ts(i,j)
   enddo
enddo

! extract travel-time values
do k=1,n_stat
   
   ! find station's neighbouring points ([2x2] for bilinear interpolation)   
   x_pos = (/ floor( xp(k) / h_step ) + 1, floor( xp(k) / h_step ) + 2 /)
   y_pos = (/ floor( yp(k) / h_step ) + 1, floor( yp(k) / h_step ) + 2 /)
                                                                           
   ! assign travel-times through bilinear interpolation
   call poly_interp(h_step*(x_pos-1),h_step*(y_pos-1),ts(x_pos,y_pos),xp(k),yp(k),time_s(k))
   call poly_interp(h_step*(x_pos-1),h_step*(y_pos-1),tp(x_pos,y_pos),xp(k),yp(k),time_p(k))
                                                                         
enddo

! close input files
close(1);close(2)
!close(3);close(4)

! deallocate memory
deallocate(tp,ts)
 

END SUBROUTINE time_distance 

!===================================================================================================

FUNCTION trapz(fun,dt)
!------------------------------------------------------------------------
! 
! Description:
! 
!    Perform intergration using trapezoidal rule over a given function.
!
! Dependencies:
!
!   None.
!
! Note:
!
!   Here step is equal to function's time-step. No checks for accuracy
!   are implemented. 
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

! function to integrate
real(kind=r_single),dimension(:),intent(in):: fun
! function's time-step
real(kind=r_single),intent(in)             :: dt
! output result
real(kind=r_single)                        :: trapz
! number of points of input time-series
integer(kind=i_single)                     :: npts

!------------------------------------------------------------------------

! find number of points to integrate
npts = size(fun)

trapz = (fun(1) + fun(npts)) / 2.

trapz = trapz + sum(fun(2:npts-1))

trapz = trapz * dt

END FUNCTION trapz
