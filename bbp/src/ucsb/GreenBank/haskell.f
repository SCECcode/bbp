c***************************************************************
c  F-K: @(#) haskell.f			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c  this file contains subroutines for computing Haskell propagator
c  matrix and other related subroutines
c	haskell()	compute the Haskell matrix
c	propagateZ()	apply the Haskell matrix to the Z matrix
c	layerParameter() compute some parameters for this layer
c       sh_ch()		compute Cosh(x), Sinh(x)/x, and x*Sinh(x)
c***************************************************************


      subroutine haskellMatrix(a)
c***************************************************************
c compute 4x4 Haskell a for the layer, the layer parameter
c is passed in by common /layer/.
c a(5,5) is for SH because the SH compound matrix is scaled by exb
c***************************************************************
      complex*16 a(5,5)
      include 'layer.h'
c due to the scaling in compound matrix, we have to scale a with exa*exb
      Ca = Ca*exb
      Xa = Xa*exb
      Ya = Ya*exb
      Cb = Cb*exa
      Yb = Yb*exa
      Xb = Xb*exa
c for p-sv, form the 4x4 haskell matrix (p381, Haskell, 1964)
      a(1,1) = r*(Ca-r1*Cb)
      a(1,2) = r*(r1*Ya-Xb)
      a(1,3) = (Cb-Ca)*r/mu2
      a(1,4) = (Xb-Ya)*r/mu2
      a(2,1) = r*(r1*Yb-Xa)
      a(2,2) = r*(Cb-r1*Ca)
      a(2,3) = (Xa-Yb)*r/mu2
      a(2,4) =-a(1,3)
      a(3,1) = mu2*r*r1*(Ca-Cb)
      a(3,2) = mu2*r*(r1*r1*Ya-Xb)
      a(3,3) = a(2,2)
      a(3,4) =-a(1,2)
      a(4,1) = mu2*r*(r1*r1*Yb-Xa)
      a(4,2) =-a(3,1)
      a(4,3) =-a(2,1)
      a(4,4) = a(1,1)
c for sh, the a is scaled by exb (in compound.f too)
      a(5,5) = exb
      return
      end


      subroutine propagateZ(z)
c***************************************************************
c  apply the Haskell matrix to Z matrix
c	Z = Z*a
c***************************************************************
      include 'constants.h'
      integer i, j, l
      complex*16 z(3,6), a(5,5), temp(4)

      call haskellMatrix(a)

      do i = 1, 3
c p-sv
         do j = 1, 4
            temp(j) = zero
            do l = 1, 4
               temp(j) = temp(j) + z(i,l)*a(l,j)
            enddo
         enddo
         do j = 1, 4
            z(i,j) = temp(j)
         enddo
c sh, could be used only once (for S*X) if not for the exb factor
         z(i,5) = z(i,5)*a(5,5)
      enddo
      return
      end


      subroutine layerParameter(k, lay)
c***************************************************************
c compute some parameters for this layer
c	IN:
c		k   --- wave-number square
c		lay --- layer number
c		volocity model passed in common/model/
c	Out:
c		common/layer/
c called by: kernel()		in kernel.f
c***************************************************************
      integer	lay
      real	k
      real*8	k2
      complex*16	kka, kkb
      include 'layer.h'
      include 'model.h'
      k2 = k*k
      kka = ka(lay)/k2
      kkb = kb(lay)/k2
      r = two/kkb
      kd = k*d(lay)
      mu2 = rho(lay)*r
c in above, because we don't have w, we use mu*(w/k)^2 instead of mu,
c the (w/k)^2 was corrected later in kernel.f and fk.f
      ra = zsqrt(one - kka)
      rb = zsqrt(one - kkb)
      r1 = one - one/r
      return
      end


      subroutine sh_ch(c,y,x,ex,a,kd)
c***************************************************************
c compute c=cosh(a*kd); y=sinh(a*kd)/a; x=sinh(a*kd)*a
c and multiply them by ex=exp(-Real(a*kd)) to supress overflow
c
c called by: compoundMatrix()		in compound.f
c***************************************************************
      complex*16 c,y,x,a
      real*8 kd,r,i,ex
      y = kd*a
      r = dreal(y)
      i = dimag(y)
      ex = dexp(-r)
      y = 0.5d0*dcmplx(dcos(i),dsin(i))
      x = ex*ex*dconjg(y)
      c = y + x
      x = y - x
      y = x/a
      x = x*a
      return
      end
