c***************************************************************
c  F-K: @(#) kernel.f			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c subroutine compute free-surface displacement kernels from a point
c source in multi-layer medium using the Haskell propagator matrix.
c
c	reference:
c		Haskell (1964), BSSA
c		Wang and Herrmann (1980), BSSA
c
c	with some modifications:
c	(1) modify the definition of B and K and put back w^2 so all matrix
c	    elements are either dimensionless or in in shear modular unit.
c	(2) suppress overflow by factoring out e^(Re(ra+rb)*d) in matrix
c	    elements.
c	(3) the cordinate is Z up, R outward and T counterclockwise
c Input:
c	k	wavenumber
c	/model/	velocity model (see include file)
c Output:
c	u(i,3)	kernels for azimuthal mode i = 0, 1, 2 in the order of Z, R, T.
c		Note that 1/2*pi*mu is omitted.
c called by
c	main 				in fk.f
c
c subroutines called:
c	layerParameter(), propagateZ()	in haskell.f
c	initialG(), propagateG()	in compound.f
c	setSource()			in compound.f
c***************************************************************
      subroutine kernel(k,u)
      include	'constants.h'
      include	'model.h'
      integer	i, lay
      real	k, k2
      complex	u(3,3), sj(3,6)
      complex*16 g(7), z(3,6)

      k2 = k*k
c propagation - from bottom to surface 
      do lay = mb, 1, -1

	call layerParameter(k, lay)
         
        if (lay.eq.mb) then
	   call initialG(g)
        else 
	   call propagateG(g)
        endif

	if (lay.eq.srclayer) then
	   call source(src,k,ka(lay)/kb(lay),rho(lay)*k2/kb(lay),sj)
	   call setSource(sj, g, z)
	endif

	if (lay.lt.srclayer) call propagateZ(z)

      enddo
      
c displacement kernel, 1 k from k-integration, k2 from k2/w2
      k2 = k*k2
      g(1) = k2/g(1)
      g(6) = k2/g(6)
      do i = 1, 3
         u(i,1) =-z(i,1)*g(1)
         u(i,2) = z(i,2)*g(1)
         u(i,3) = z(i,5)*g(6)
      enddo

      return
      end
