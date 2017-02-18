c***************************************************************
c  F-K: @(#) model.h			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c Inlucde header file pertaining the velocity model info
c       velocity model in common block /model/:
c               mb---bottom layer
c           srclayer---source layer (on top of it)
c		src-- source type (0=ex;1=sf;2=dc)
c               ka--- (w/vp)^2
c               kb--- (w/vs)^2
c               d---- thickness of the layer
c               rho--- density of the layer
c		ch--- mu/bulk_modular, this one is for static computation
      integer mb,src,srclayer,nlay,ndis,nt
c max. # of layers and receivers and time length
      parameter(nlay = 500, ndis=630, nt=4096)
      real d(nlay),rho(nlay),ch(nlay)
      complex ka(nlay),kb(nlay)
      common /model/mb,src,srclayer,ka,kb,d,rho,ch
