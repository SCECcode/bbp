
! --------------------------- BEGIN GET_LUN ----------------
      subroutine get_lun(lun)

! Finds a logical unit number not in use; returns
! -1 if it cannot find one.

! Dates -- 05/19/98 - Written by D. Boore, following
!                     Larry Baker's suggestion
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      logical isopen
      do i = 99,10,-1
        inquire (unit=i, opened=isopen)
        if(.not.isopen) then
          lun = i
          return
        end if
      end do
      lun = -1

      return
      end
! --------------------------- END GET_LUN ----------------
     
