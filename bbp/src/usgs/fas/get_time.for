!--------------------- BEGIN GET_TIME  ---------------------
      subroutine get_time(time_c)

! Dates: 04/29/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      implicit none


      character time_c*(*)
!      integer*2 ihr, imin, isec, i100th
      character datx*8, timx*10


! Microsoft compiler version:
!      call GETTIM (ihr, imin, isec, i100th)


!      time_c = ' '             ! character time_c*11 for both compilers


!      time_c(3:3) = ':'
!      time_c(6:6) = ':'
!      time_c(9:9) = '.'


!      write(time_c(1:2),'(i2.2)') ihr      
!      write(time_c(4:5),'(i2.2)') imin      
!      write(time_c(7:8),'(i2.2)') isec      
!      write(time_c(10:11),'(i2.2)') i100th


! Lahey compiler version:
!      call TIME(time_c)


! Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
      call DATE_AND_TIME( datx, timx )
! Date is returned as 'CCYYMMDD'
! Time is returned as 'hhmmss.sss'
! time_c should be declared as
! character time_c*11
      time_c = ' '
      time_c = timx(1:2)//':'//timx(3:4) // ':' // timx(5:9)


      return
      end
!--------------------- END GET_TIME  --------------------- 