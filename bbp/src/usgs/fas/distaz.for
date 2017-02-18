
! ------------------------ BEGIN DISTAZ ----------------------
      subroutine distaz( wlongsign, alat, along, blat, blong,
     * rdeg, rkm, az, baz)
! 
! compute distances, azimuths using formulas from
! Bruce Julian.
!
! latest modification: 1/27/84
!
!        11/24/07 - Trap for colocated pair in subroutine distaz
!        05/01/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      if (alat .eq. blat .and. along .eq. blong) then
        rdeg = 0.0
        rkm = 0.0
        az = 0.0
        baz = 0.0
        return
      end if

      pi = 4.0 * atan( 1. )
      dtor = pi/ 180.
!
! convert from degrees to radians and correct sign of
! longitude so that east longitude is positive.
!
      alatr = dtor * alat
      alongr = -dtor * along * wlongsign
      blatr = dtor * blat
      blongr = -dtor * blong * wlongsign
!
! compute geocentric latitudes.
!
      alatr = atan( 0.993305 * tan( alatr ) )
      blatr = atan( 0.993305 * tan( blatr ) )
!
! compute latitude dependent quantities
!
      ca = cos( alatr )
      cb = cos( blatr )
      sa = sin( alatr )
      sb = sin( blatr )
!
! now compute other quantities
!
      a = cb * sin( blongr - alongr )
      b = ca * sb - sa * cb * cos( blongr - alongr )
      cd = ca * cb * cos( blongr - alongr ) + sa * sb
      sd = sqrt( a*a + b*b )
!
! compute distances
!
      rdeg = atan2( sd, cd )/ dtor
      rkm = 111.19 * rdeg
!
! compute azimuth (from a to b) and make it positive.
!
      az = atan2( a, b )/ dtor
      if ( az .lt. 0.0 ) az = az + 360.0
!
! compute back azimuth (from b to a) and make it positive.
!
      a = ca * sin( alongr - blongr )
      b = cb * sa - sb * ca * cos( alongr - blongr )
      baz = atan2( a, b)/ dtor
      if ( baz .lt. 0.0 ) baz = baz + 360.0
!
      return
      end
! ------------------------ END DISTAZ ----------------------

