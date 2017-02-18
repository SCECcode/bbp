! --------------------- BEGIN YINTRF ------------------------------------
      function yintrf( x, xin, yin, n)
!
! returns an interpolated value (yintrf) based on straight line
! interpolation of the data in xin and yin.

! Needs Numerical recipe routine locate

!
! dates:  3/14/85 - written
!        11/30/95 - substituted LOCATE instead of starting from beginning
!                   each time
!        03/13/96 - added code to deal with xin increasing or decreasing
!        12/12/00 - Stripped off "locate.for"
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      dimension xin(1), yin(1)
      logical incrs

! Is xin increasing or decreasing?
      incrs = .true.
      if (xin(n) .lt. xin(1)) incrs = .false.

! Set value if x is outside the range of xin:
      if (incrs) then
        if ( x .le. xin(1) ) then
            yintrf = yin(1)
            return
        end if
        if ( x .ge. xin(n) ) then
            yintrf = yin(n)
            return
        end if  
      else
        if ( x .ge. xin(1) ) then
            yintrf = yin(1)
            return
        end if
        if ( x .le. xin(n) ) then
            yintrf = yin(n)
            return
        end if  
      end if

! Locate the proper cell and interpolate:
      call locate(xin, n, x, j)
      yintrf = yin(j) + (x-xin(j))*(yin(j+1) - yin(j))/
     * (xin(j+1)-xin(j))

      return
      end
! --------------------- END YINTRF ------------------------------------

