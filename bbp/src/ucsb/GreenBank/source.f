c***************************************************************
c  F-K: @(#) source.f			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c Setup the displacement-stress vector jump across the source
c The coefficients can be all derived from the potential jump in
c Haskell (1964), or directly from Takiuchi and Saito (1972)
c
c  Input:
c	src	--- source type, 0=ex; 1=sf; 2=dc
c	k	--- wavenumber
c	xi	--- mu/(lambda+2*mu) at the source
c	mu	--- shear module
c  Output:
c	s(3,6)	--- source coef. for n=0, 1, 2
c
c  called by:	kernel()	in kernel.f
c
c  modified history:
c	May 12, 2000	Lupei Zhu	initial coding
c	July 17, 2000	Lupei Zhu	change xi, mu to be complex
c***************************************************************
      subroutine source(src, k, xi, mu, s)
      integer i, j, src
      real k
      complex s(3,6), xi, mu

      do i = 1, 3
         do j = 1, 6
	    s(i,j) = 0.
	 enddo
      enddo

      if (src .EQ. 2) then 		! DOUBLE_COUPLE
         s(1,2) = 2.*xi/mu
         s(1,4) = 4.*xi-3.
         s(2,1) = 1./mu
         s(2,5) =-s(2,1)
         s(3,4) = 1.
         s(3,6) =-1.
      else if (src .EQ. 0) then		! EXPLOSION
         s(1,2) = xi/mu
         s(1,4) = 2.*xi
      else if (src .EQ. 1) then		! SINGLE_FORCE
         s(1,3) = -1./k
         s(2,4) = -s(1,3)
         s(2,6) =  s(1,3)
      else				! unknow source type
	 write(0,*)'unknown source type'
	 stop
      endif

      return
      end
