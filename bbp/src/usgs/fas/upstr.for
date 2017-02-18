! --------------------- BEGIN UPSTR ----------------------------------
      Subroutine UPSTR ( text )
! Converts character string in TEXT to uppercase
! Dates: 03/12/96 - Written by Larry Baker
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

!
      Implicit   None
!
      Character  text*(*)
!
      Integer    j
      Character  ch
!
      Do 1000 j = 1,LEN(text)
         ch = text(j:j)
         If ( LGE(ch,'a') .and. LLE(ch,'z') ) Then
            text(j:j) = CHAR ( ICHAR(ch) - ICHAR('a') + ICHAR('A') )
         End If
 1000    Continue
!
      Return
      End
! --------------------- END UPSTR ----------------------------------
