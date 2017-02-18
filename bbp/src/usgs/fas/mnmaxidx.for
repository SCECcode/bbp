
! ---------------------- BEGIN MNMAXIDX ----------------------
      subroutine mnmaxidx(a,nstrt,nstop,ninc,
     :                    amin,amax,indx_min,indx_max) 

! A rewrite of mnmax, returning the indices of the min and max values
!
! Dates: 10/05/99 - written by D. M. Boore
!        05/16/02 - Replaced "dimension a(1)" with "real a(*)"
!        11/03/12 - Replace if-then with if--then-else

      real a(*)
      amax = a( nstrt)
      indx_max = nstrt
      amin=amax         
      indx_min = nstrt
      do i=nstrt,nstop,ninc
      
        if (a(i) > amax) then
          amax=a(i)
          indx_max = i
        else if (a(i) < amin) then
          amin=a(i)
          indx_min = i
        end if
        
      end do
      
      return                  
      end                     
! ---------------------- END MNMAXIDX ----------------------

