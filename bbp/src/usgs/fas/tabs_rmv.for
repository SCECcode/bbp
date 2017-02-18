
! ---------------------------------------------------------------- Tabs_Rmv
      subroutine tabs_rmv(string, nchar)

! Replaces tabs in the character string with a blank

! Dates: 11/04/00 - Written by D. Boore
!        10/29/02 - Renamed from rmv_tabs.for
!        05/01/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      character string*(*)


      do i = 1, nchar
        if(ichar(string(i:i)) .eq. 9) then
           string(i:i) = ' '
        end if
      end do

      return
      end
! ---------------------------------------------------------------- Tabs_Rmv

