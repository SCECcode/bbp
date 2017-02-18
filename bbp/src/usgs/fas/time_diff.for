
! --------------- BEGIN TIME_DIFF ---------------------------------
      subroutine time_diff(time_start, time_stop, time_elapsed)

! Dates: 02/18/09 - Written by D.M. Boore
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)
! To be used with
! Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
!      character datx*8, timx*10
!      call DATE_AND_TIME( datx, timx )
! Date is returned as 'CCYYMMDD'
! Time is returned as 'hhmmss.sss'

      implicit none
      character, intent(in) :: time_start*(*), time_stop*(*)
      real, intent(out) :: time_elapsed
      real ::  secb, sece 
      integer :: ihb, imb, ihe, ime

      read(time_start(1:10),'(i2,i2,f6.3)') 
     :                       ihb, imb, secb
      read(time_stop(1:10),'(i2,i2,f6.3)') 
     :                       ihe, ime, sece
      time_elapsed = 
     :  3600.0*float(ihe-ihb) + 60.0*float(ime-imb) + sece-secb 

      end subroutine time_diff
! --------------- END TIME_DIFF ---------------------------------


