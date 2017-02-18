!--------------------- BEGIN GET_DATE  ---------------------
      subroutine get_date(date_c)

! Dates: 04/29/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)


      implicit none


      character date_c*(*)
!      integer*2 iyr, imon, iday
      character datx*8, timx*10


! Microsoft compiler version:
!      call GETDAT (iyr, imon, iday)


!      date_c = ' '         ! character date_c*10


!      date_c(3:3) = '/'
!      date_c(6:6) = '/'


!      write(date_c(1:2),'(i2.2)') imon      
!      write(date_c(4:5),'(i2.2)') iday      
!      write(date_c(7:10),'(i4.4)') iyr      


! Lahey compiler version:
!      call DATE(date_c)   ! character date_c*8, but *10 is OK


! Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
      call DATE_AND_TIME( datx, timx )
! Date is returned as 'CCYYMMDD'
! Time is returned as 'hhmmss.sss'
! date_c should be declared as
! character date_c*10
      date_c = ' '
      date_c = datx(5:6) // '/' // datx(7:8) // '/' // datx(1:4)


      return
      end
!--------------------- END GET_DATE  ---------------------