c***************************************************************
c  F-K: @(#) compound.f			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c  This file contains subroutines for computing the compound matrix
c  and related routines.
c
c  For reference, see Wang and Herrmann, BSSA, 1980
c
c       compoundMatrix(a)	compute the compound matrix
c       propagateG()		apply the Haskell matrix to the Z matrix
c       initialG()		initialize g matrix in the half space
c       setSource()		apply g to the source terms
c***************************************************************


      subroutine compoundMatrix(a)
c***************************************************************
c compute 7x7 compound matrix a for the layer (passed in by /layer/
c The upper-left 5x5 is the compound matrix of the Haskell matrix
c with the 3rd row and colummn dropped and the 4th row being replaced
c by (2A41, 2A42, 2A44-1,2A45,2A46) (see W&H, P1035). The lower-right
c 2x2 is the SH part of the Haskell matrix.
c***************************************************************
      real*8 ex
      complex*16 a(7,7)
      complex*16 CaCb, CaXb, CaYb, XaCb, XaXb, YaCb, YaYb, r2, r3
      include 'layer.h'
      call sh_ch(Ca,Ya,Xa,exa,ra,kd)
      call sh_ch(Cb,Yb,Xb,exb,rb,kd)
      CaCb=Ca*Cb
      CaYb=Ca*Yb
      CaXb=Ca*Xb
      XaCb=Xa*Cb
      XaXb=Xa*Xb
      YaCb=Ya*Cb
      YaYb=Ya*Yb
      ex = exa*exb
      r2 = r*r
      r3 = r1*r1
c p-sv, scaled by exa*exb to supress overflow
      a(1,1) = ((one+r3)*CaCb-XaXb-r3*YaYb-two*r1*ex)*r2
      a(1,2) = (XaCb-CaYb)*r/mu2
      a(1,3) = ((one+r1)*(CaCb-ex)-XaXb-r1*YaYb)*r2/mu2
      a(1,4) = (YaCb-CaXb)*r/mu2
      a(1,5) = (two*(CaCb-ex)-XaXb-YaYb)*r2/(mu2*mu2)
      a(2,1) = (r3*YaCb-CaXb)*r*mu2
      a(2,2) = CaCb
      a(2,3) = (r1*YaCb-CaXb)*r
      a(2,4) =-Ya*Xb
      a(2,5) = a(1,4)
      a(3,1) = two*mu2*r2*(r1*r3*YaYb-(CaCb-ex)*(r3+r1)+XaXb)
      a(3,2) = two*r*(r1*CaYb-XaCb)
      a(3,3) = two*(CaCb - a(1,1)) + ex
      a(3,4) =-two*a(2,3)
      a(3,5) =-two*a(1,3)
      a(4,1) = mu2*r*(XaCb-r3*CaYb)
      a(4,2) =-Xa*Yb
      a(4,3) =-a(3,2)/two
      a(4,4) = a(2,2)
      a(4,5) = a(1,2)
      a(5,1) = mu2*mu2*r2*(two*(CaCb-ex)*r3-XaXb-r3*r3*YaYb)
      a(5,2) = a(4,1)
      a(5,3) =-a(3,1)/two
      a(5,4) = a(2,1)
      a(5,5) = a(1,1)
c sh, scaled by exb
      a(6,6) = Cb
      a(6,7) =-two*Yb/mu2
      a(7,6) =-mu2*Xb/two
      a(7,7) = Cb
      return
      end


      subroutine propagateG(g)
c***************************************************************
c propagate 1x7 g matrix upward using the compound matrix
c	g = g*a
c***************************************************************
      include 'constants.h'
      integer i, j
      complex*16 g(7), a(7,7), temp(7)

      call compoundMatrix(a)

c p-sv
      do i = 1, 5
         temp(i) = zero
         do j = 1, 5
            temp(i) = temp(i) + g(j)*a(j,i)
         enddo
      enddo
      do i = 1, 5
         g(i) = temp(i)
      enddo
c sh
      do i = 6, 7
         temp(i) = zero
         do j = 6, 7
            temp(i) = temp(i) + g(j)*a(j,i)
         enddo
      enddo
      do i = 6, 7
         g(i) = temp(i)
      enddo
      return
      end


      subroutine initialG(g)
c***************************************************************
c Initialize the 1x7 g-matrix. The first 5 elements are the first
c row of the compound matrix of inverse(E), with the 3rd column dropped.
c The SH part (placed in g(6) and g(7)) is 5th row of E^-1
c***************************************************************
      complex*16 g(7), delta
      include 'layer.h'
c p-sv, see WH p1034
      delta  = r*(one-ra*rb)-one
      g(1) = mu2*(delta-r1)
      g(2) = ra
      g(3) = delta
      g(4) =-rb
      g(5) = (one+delta)/mu2
c sh, use the 5th row of E^-1
      g(6) =-one
      g(7) = two/(rb*mu2)
      return
      end


      subroutine setSource(s, g, z)
c***************************************************************
c put source term and apply X matrix to it
c  input:
c	s(3,6)	---- source coef.
c	g(7)	---- g vector used to construct matrix X
c		     |	0 -g1 -g2 +g3 |
c		X =  |	g1 0  -g3 -g4 | for P-SV
c		     |	g2 g3  0  -g5 |
c		     | -g3 g4  g5  0  |
c		X =	| g6 g7 |	for SH
c  output:
c	z(3,6)  ---- s*X where s(3,6) is the source coef. n=0,1,2
c***************************************************************
      include 'constants.h'
      integer	i
      complex	s(3,6)
      complex*16 g(7), z(3,6)

      do i=1, 3
c for p-sv, see WH p1018
	 z(i,1)=-s(i,2)*g(1)-s(i,3)*g(2)+s(i,4)*g(3)
	 z(i,2)= s(i,1)*g(1)-s(i,3)*g(3)-s(i,4)*g(4)
	 z(i,3)= s(i,1)*g(2)+s(i,2)*g(3)-s(i,4)*g(5)
	 z(i,4)=-s(i,1)*g(3)+s(i,2)*g(4)+s(i,3)*g(5)
c for sh
	 z(i,5)= s(i,5)*g(6)+s(i,6)*g(7)
      enddo
      return
      end
