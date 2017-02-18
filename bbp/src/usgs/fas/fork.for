
! ----------------------------- BEGIN FORK --------------------------
      SUBROUTINE FORK(LX,CX,SIGNI)
! FAST FOURIER                                  2/15/69
!                          LX
!    CX(K) = SQRT(1.0/LX)* SUM (CX(J)*EXP(2*PI*SIGNI*I*(J-1)*(K-1)/LX))
!                          J=1                        FOR K=1,2,...,LX
!
!  THE SCALING BETWEEN FFT AND EQUIVALENT CONTINUUM OUTPUTS
!  IS AS FOLLOWS.
!
!
!     GOING FROM TIME TO FREQUENCY:
!             F(W)=DT*SQRT(LX)*CX(K)
!
!                  WHERE W(K)=2.0*PI*(K-1)*DF

!                  and    DF = 1/(LX*DT)
!
!
!     GOING FROM FREQUENCY TO TIME, WHERE THE FREQUENCY
!     SPECTRUM IS GIVEN BY THE DIGITIZED CONTINUUM SPECTRUM:
!
!             F(T)=DF*SQRT(LX)*CX(K)
!
!                  WHERE T(K)=(K-1)*DT
!
!
!  THE RESULT OF THE SEQUENCE...TIME TO FREQUENCY,POSSIBLE MODIFICATIONS
!  OF THE SPECTRUM (FOR FILTERING,ETC.), BACK TO TIME...
!  REQUIRES NO SCALING.
!
!
!  THIS VERSION HAS A SLIGHT MODIFICATION TO SAVE SOME TIME...
!  IT TAKES THE FACTOR 3.1415926*SIGNI/L OUTSIDE A DO LOOP (D.BOORE 12/8
!  FOLLOWING A SUGGESTION BY HENRY SWANGER).
!

! Some brief notes on usage:

! "signi" is a real variable and should be called either with the value "+1.0"
! of "-1.0".  The particular value used depends on the conventions being used
! in the application (e.g., see Aki and Richards, 1980, Box 5.2, pp. 129--130).

! Time to frequency:
! In calling routine,
! 
!       do i = 1, lx
!         cx(i) = CMPLX(y(i), 0.0)
!       end do
!  where y(i) is the time series and lx is a power of 2
! 
!  After calling Fork with the complex array specified above, the following 
! symmetries exist:
! 
!        cx(1)        = dc value (f = 0 * df, where df = 1.0/(lx*dt))
!        cx(lx/2 + 1) = value at Nyquist (f = (lx/2+1-1)*df = 1.0/(2*dt))
!        cx(lx)       = CONJG(cx(2))
!        cx(lx-1)     = CONJG(cx(3))
!         |           =      |
!        cx(lx-i+2)   = CONJG(cx(i))
!         |           =      |
!        cx(lx/2+2)   = CONJG(cx(lx/2))
! 
! where "CONJG" is the Fortran complex conjugate intrinsic function
! 
! This symmetry MUST be preserved if modifications are made in the frequency 
! domain and another call to Fork (with a different sign for signi) is used
! to go back to the time domain.  If the symmetry is not preserved, then the
! time domain array will have nonzero imaginary components.  There is one case
! where advantage can be taken of this, and that is to find the Hilbert 
! transform and the window of a time series with only two calls to Fork (there 
! is a short note in BSSA {GET REFERENCE} discussing this trick, which amounts 
! to zeroing out the last half of the array and multiplying all but the dc and 
! Nyquist values by 2.0; in the time domain, REAL(cx(i)) and AIMAG(cx(i)) 
! contain the filtered (if a filter was applied) and Hilbert transform of the 
! filtered time series, respectively, while CABS(cx(i)) and ATAN2(AIMAG(cx(i)), 
! REAL(cx(i))) are the window and instantaneous phase of the filtered time 
! series, respectively.

! Some references:

! Farnbach, J.S. (1975). The complex envelope in seismic signal analysis, 
! BSSA 65, 951--962. 
! He states that the factor of 2 is applied for i = 2...npw2/2 (his indices 
! start at 0, I've added 1), which is different than the next reference:

! Mitra, S.K. (2001). Digital Signal Processing, McGraw-Hill, New York.
! He gives an algorithm on p. 794 (eq. 11.81), in which the factor of 2 is 
! applied from 0 frequency to just less than Nyquist.

! 
! The easiest way to ensure the proper symmetry is to zero out the
! last half of the array (as discussed above), but the following is what
! I usually use:  
! modify (filter) only half
! of the cx array:
! 
!       do i = 1, lx/2
!         cx(i) = filter(i)*cx(i)
!       end do
! 
! where "filter(i)" is a possibly complex filter function (and recall that 
! the frequency corresponding to i is f = float(i-1)*df).  After this, fill out
! the last half of the array using
!       
!       do i = lx/2+2, lx
!         cx(i) = CONJG(cx(lx+2-j))
!       end do
! 
! Note that nothing is done with the Nyquist value.  I assume (but am not sure!)
! that this value should be 0.0
! 
! Dates: xx/xx/xx - Written by Norm Brenner(?), Jon Claerbout(?)
!        12/21/00 - Replaced hardwired value of pi with pi evaluated here,
!                     and added comments regarding usage.  Also deleted
!                     dimension specification of cx(lx) and replace it with
!                     cx(*) in the type specification statement.  I also
!                     cleaned up the formatting of the subroutine.
!        08/28/01 - Added comment about variable "signi" being real, and 
!                   added "float" in equations for "sc" and "temp", although 
!                   not strictly required.
!        06/19/02 - Added some comments about computing envelopes and
!                   instantaneous frequencies
!        01/19/15 - Modernize code (get rid of go to statements)
             
      complex cx(*),carg,cexp,cw,ctemp

      pi = 4.0*atan(1.0)

      j=1
      sc=sqrt(1./real(lx))

      do i=1,lx
      
        if (i <= j) then
          ctemp=cx(j)*sc
          cx(j)=cx(i)*sc
          cx(i)=ctemp
        end if
        
        m=lx/2        
        
        DO
          if (j <= m) EXIT
          j=j-m
          m=m/2
          if (m < 1) EXIT
        END DO
        
        j = j + m
        
      end do

      l=1
      DO WHILE (l < lx)
        istep=2*l
        temp= pi * signi/real(l)

        do m=1,l
          carg=(0.,1.)*temp*(m-1)
          cw=cexp(carg)
          do i=m,lx,istep
            ctemp=cw*cx(i+l)
            cx(i+l)=cx(i)-ctemp
            cx(i)=cx(i)+ctemp
          end do
        end do

        l=istep
      END DO
 
! Previous code (as of 18 July 2010):
!      pi = 4.0*atan(1.0)
!
!      j=1
!      sc=sqrt(1./float(lx))
!
!      do i=1,lx
!        if(i.gt.j) go to 2
!        ctemp=cx(j)*sc
!        cx(j)=cx(i)*sc
!        cx(i)=ctemp
!2       m=lx/2
!3       if(j.le.m) go to 5
!        j=j-m
!        m=m/2
!        if(m.ge.1) go to 3
!5       j=j+m
!      end do
!
!      l=1
!6     istep=2*l
!      temp= pi * signi/float(l)
!
!      do m=1,l
!        carg=(0.,1.)*temp*(m-1)
!        cw=cexp(carg)
!        do i=m,lx,istep
!          ctemp=cw*cx(i+l)
!          cx(i+l)=cx(i)-ctemp
!          cx(i)=cx(i)+ctemp
!        end do
!      end do
!
!      l=istep
!      if(l.lt.lx) go to 6

      return
      end
! ----------------------------- END FORK --------------------------

