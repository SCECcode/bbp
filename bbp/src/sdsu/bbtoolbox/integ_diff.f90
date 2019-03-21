SUBROUTINE vel2acc(in_data,dt,nt,out_data)    
!---------------------------------------------------------------------------
!
! Description:
!   from vel2acc.m
!   
!   
!
! Dependencies:
!
!   None.
!
! Author: R. Takedatsu
!
! Modified: February 2019 (v2.0)
!   
   
use def_kind;

implicit none

! local variables
integer(kind=i_single),intent(in)             :: nt
integer(kind=i_single)                        :: k
real(kind=r_single),intent(in)                :: dt              
real(kind=r_single),dimension(nt),intent(in)  :: in_data              
real(kind=r_single),dimension(nt),intent(out) :: out_data              

! Originally whitten by D. Roten in Matlab code
!function [acc]=vel2acc(vel, dt)
!   nt=length(vel);
!   acc=zeros(nt, 1);
!   for k=2:nt
!      acc(k)=(vel(k) - vel(k-1))/dt;
!   end
!   return

!---------------------------------------------------------------------------

out_data(1:nt)=0

do k=2,nt
   out_data(k)=(in_data(k) - in_data(k-1)) / dt
enddo

END SUBROUTINE vel2acc

!===================================================================================================

SUBROUTINE acc2velc(ival,in_data,dt,nt,out_data)    
!---------------------------------------------------------------------------
!
! Description:
!   from acc2velc.m
!   
!   
!
! Dependencies:
!
!   None.
!
! Author: R. Takedatsu 
!
! Modified: February 2019 (v2.0)
! 
   
use def_kind;

implicit none

! local variables
integer(kind=i_single),intent(in)             :: nt
integer(kind=i_single)                        :: k
real(kind=r_single),intent(in)                :: ival ! initial velocity value              
real(kind=r_single),intent(in)                :: dt              
real(kind=r_single),dimension(nt),intent(in)  :: in_data              
real(kind=r_single),dimension(nt),intent(out) :: out_data              

! Originally whitten by D. Roten in Matlab code
!function [velc]=acc2velc(vel, acc, dt)
!   nt=length(acc);
!   velc=zeros(nt, 1);
!   velc(1)= vel(1);
!   for k=2:nt
!       velc(k)=velc(k-1)+acc(k)*dt;
!   end

!---------------------------------------------------------------------------

out_data(1:nt)=0
out_data(1)=ival

do k=2,nt
   out_data(k)=out_data(k-1) + in_data(k)* dt
enddo

END SUBROUTINE acc2velc

!===================================================================================================
