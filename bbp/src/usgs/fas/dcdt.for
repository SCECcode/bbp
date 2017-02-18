
! ------------------------ begin dcdt -------------------
      subroutine dcdt (y,dt,npts,indx1,indx2,ldc,ldt)
!+
!  dcdt - fits dc or trend between indices indx1 and indx2.
!         then removes dc or detrends whole trace.
!         y is real, dt = delta t.
!         if remove dc, ldc = .true.
!         if detrend, ldt = .true.
!-

! Dates: 12/14/00 - Cleaned up formatting of original program
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      real y(*)
      logical ldc,ldt

      if (.not. ldc .and. .not. ldt) then
        return
      end if

!
!...fit dc and trend between indices indx1 and indx2.
      nsum = indx2-indx1+1
      sumx = 0.0
      sumx2 = 0.0
      sumy = 0.0
      sumxy = 0.0
      do i=indx1,indx2
         xsubi = (i-1)*dt
         sumxy = sumxy+xsubi*y(i)
         sumx = sumx+xsubi
         sumx2 = sumx2+xsubi*xsubi
         sumy = sumy+y(i)
      end do
!
!... remove dc.
      if (ldc) then
        avy = sumy/nsum
        do i=1,npts
          y(i) = y(i)-avy
        end do
! Debug
        write(*,'(a)') ' indx1, indx2, avy'
        write(*, *)      indx1, indx2, avy
! Debug



        return
      endif
!
!... detrend. see draper and smith, p. 10.
      if (ldt) then
        bxy = (sumxy-sumx*sumy/nsum)/(sumx2-sumx*sumx/nsum)
        axy = (sumy-bxy*sumx)/nsum
        qxy = dt*bxy
        do i=1,npts
          y(i) = y(i)-(axy+(i-1)*qxy)
        end do
        return
      endif
!
      return
      end
! ------------------------ end dcdt -------------------

