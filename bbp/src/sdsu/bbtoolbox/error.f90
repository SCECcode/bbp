SUBROUTINE error_handling(err_id,object,topo)
!-------------------------------------------------------------------
!
! Description:
!
!   Manages code failures, giving messages and useful informations 
!   when code fails.
!
! Dependencies:
!
!   None.
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
 
use def_kind 
 
implicit none

! error identity number
integer(kind=i_single),intent(in):: err_id
! identify where error occurred
character(len=*),intent(in)      :: topo,object

!-------------------------------------------------------------------

! errors block
select case(err_id)
case(1)
   print*, 'Error while opening file: ', trim(object)
   print*, 'Error in subroutine: ', trim(topo)
   stop
case(2)
   print*, 'Error while writing file: ', trim(object)
   print*, 'Error in subroutine: ', trim(topo)
   stop
case(3)
   print*, 'Error while erasing file: ',trim(object)
   print*, 'Error in subroutine: ', trim(topo)
   stop
case(4)
   print*, 'Error while interpolating time-series: two equal x-values '
   print*, 'Interpolation failed - please check your interpolation coordinates'
   print*, 'Error in subroutine: ', trim(topo)
   stop   
case(5)
   print*, 'Error: not 2**n number of points'
   print*, 'Fast Fourier subroutine will fail - please change npts'
   print*, 'Error in subroutine: ', trim(topo)
   stop
case(6)
   print*, 'Error while merging LF and HF: too large search-window'
   print*, 'Please change target frequency or bandwidth for search'
   print*, 'Error in subroutine: ', trim(topo)
   stop   
case(7)   
   print*, 'Error while interpolating time-series: x and y series have different length'
   print*, 'Interpolation failed - please check your input time-series'
   print*, 'Error in subroutine: ', trim(topo)
   stop 
case(8)
   print*, 'Error: ISOCHRONE approach not yet implemented!'
   print*, 'Please select another operative modality (i.e. 0 or 1)'
   print*, 'Error in subroutine: ', trim(topo)
   stop 
case(9)
   print*, 'Error while convolving time-series'
   print*, 'Input b must have less points than input a'
   print*, 'Error in subroutine: ', trim(topo)
   stop   
case(10)
   print*, 'Error while reading input file' 
   print*, 'Please insert a valid parameter: ',trim(object)
   print*, 'Error in subroutine: ', trim(topo)
   stop 
case(11)
   print*, 'Error: LF and HF time-series have no equal length or time sampling' 
   print*, 'Please check your input time-series'
   print*, 'Error in subroutine: ', trim(topo)
   stop      
end select


END SUBROUTINE error_handling

!===================================================================================================

SUBROUTINE warning_handling(warn_id,object,topo)
!-------------------------------------------------------------------
!
! Description:
!
!   Manages code warnings, giving messages and useful informations 
!   when code encounters difficulties
!
! Dependencies:
!
!   None
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
 
use def_kind 
 
implicit none

! error identity number
integer(kind=i_single),intent(in):: warn_id
! identify where error occurred
character(len=*),intent(in)      :: topo,object

!-------------------------------------------------------------------

select case(warn_id)
case(1)
   print*, 'Warning for tapering broadband: too long taper, roughly shortened'
   print*, 'Warning in subroutine: ', trim(topo)
case(2)
   print*, 'Warning for scaling window: too wide window, roughly adjusted'
   print*, 'Warning in subroutine: ', trim(topo)
case(3)
   print*, 'Warning for reading input time-series: incorrect Rob Graves format'
   print*, 'Number of points not a multiple of 6: resulting time-series will be shorter'
   print*, 'Warning in subroutine: ', trim(topo)
end select

END SUBROUTINE warning_handling