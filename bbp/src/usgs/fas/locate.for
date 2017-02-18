
!* --------------------- BEGIN LOCATE -----------------
      SUBROUTINE locate(xx,n,x,j)
      
! Comments added by D. Boore on 26feb2010:
!  finds j such that xx(j) < x <= xx(j+1)
!  EXCEPT if x = xx(1), then j = 1 (logically it would be 0 from
!  the above relation, but this is the same returned value of j
!  for a point out of range).
!  Also, if x = xx(n), j = n-1, which is OK
!  Note that j = 0 or j = n indicates that x is out of range.
!
! See the program test_locate.for to test this routine.

! Dates: 04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
!* --------------------- END LOCATE -----------------

