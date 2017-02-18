
! ------------------------------------------------------------- Get_NPW2
      subroutine get_npw2(npts,signnpw2,npw2)

! Find npw2 (less than npts if signnpw2 < 0)

! Dates: 12/12/00 - Written by D. Boore
!        04/04/10 - Correct error that it does not return
!                   npw2 = npts, if npts is a power of 2.
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      npw2_exp = int( alog(float(npts))/alog(2.0) )
      if (signnpw2 < 0.0) then
        npw2 = 2.0**npw2_exp
      else 
        npw2_temp = 2.0**npw2_exp
        if (npw2_temp == npts) then 
          npw2 = npts
        else
          npw2 = 2.0**(npw2_exp+1)
        end if
      end if

      return
      end
! ------------------------------------------------------------- Get_NPW2

