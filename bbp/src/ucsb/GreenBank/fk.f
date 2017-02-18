c**********************************************************************
c  F-K:	@(#) fk.f			1.0 4/29/2000
c
c  Copyright (c) 1999 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c Compute surface displacements from a point source in layered
c media using the Thompson-Haskell propagator matrix.
c The source can have azimuthal mode number up to n=2, including explosion,
c single force, and double-couple. All source time functions are Dirac delta.
c
c The outputs are the surface displacements (in SAC format) in the order of
c vertical (UP), radial, and tangential (couterclockwise) for n=0, 1, 2
c (i.e. Z0, R0, T0, Z1, ...).
c Their units are (assume v in km/s, rho in g/cm^3, thickness in km):
c	10^-20 cm/(dyne cm)	for double couple source and explosion
c	10^-15 cm/dyne		for single force
c
c The code uses summation for wave-number (k) integration and then FFT for
c omega integration.
c
c subroutines called:
c	kernel() in kernel.f	computing Kernels U(k,w)
c
c Written by Lupei Zhu, Seismo Lab, Caltech (lupei@gps.caltech.edu)
c
c Modified history:
c	03/05/1996  Lupei Zhu	Initial coding
c	03/07/1999  Lupei Zhu	use conditional compiling option
c				to consolidate the options of
c				exact or approx. bessel functions
c	05/06/1999  Lupei Zhu	remove the division of (i*w) which
c				was put for getting the step-response
c				of the double-couple source
c	04/29/2000  Lupei Zhu   determine kmax using source depth hs
c	05/06/2000  Lupei Zhu	use total and count to estimate progress
c	05/06/2000  Lupei Zhu	use array instead of writing temp to disk
c	07/17/2000  Lupei Zhu	treat mu as complex and move source()
c				back into the w-k loops
c	07/26/2000  Lupei Zhu	modify to do sigma correction on real
c				array. The previous on does this on
c				complex array and causes long-period noise
c**********************************************************************
c
c      program fk
      Subroutine sub_bs_dc(nx,x,r_max,lnpt,dt,t0,sum)
      include 'constants.h'
      include 'model.h'
      integer i,j,l,iblank,nfft,nfft2,n,ix,nx,tenPerc,count,total
      INTEGER nCom, wc, tb, smth ,new_old(9)
      real k,omega,dt,dk,dw,sigma,const,phi,hs,xmax,vs
      real dfac,pmin,pmax,kc,lowpass,taper
      real qa(nlay),qb(nlay),a(nlay),b(nlay),x(ndis),t0(ndis)
      complex w,at0,at,grt,u(3,6),cmu_freq
      real aj0,aj1,aj2,z,ddr,ddi
      complex sum(ndis,9,nt), data(nt),tdata(9,nt)
      integer jo,nb,ncom_out,id_src
      real den(nlay),th(nlay),sig(9)
      COMMON/MODEL_new/jo,nb,a,b,den,th,qa,qb,ncom_out
      nCom = 9
c
c	sequence 1   2   3  4    5   6   7   8
c               tss tds rss rds r45 zss zds z45
      id_src=2
      if(id_src.eq.2) ncom_out=8

      new_old(1)=9
      new_old(2)=6
      new_old(3)=8
      new_old(4)=5
      new_old(5)=2
      new_old(6)=7
      new_old(7)=4
      new_old(8)=1
      sig(1)=1.0
      sig(2)=-1.0
      sig(3)=-1.0
      sig(4)=1.0
      sig(5)=1.0
      sig(6)=-1.0
      sig(7)=1.0
      sig(8)=1.0
c***************************************************************
c input velocity model
c      read(*,*)mb,srclayer,src

      mb=jo
      srclayer=nb+1
      src=id_src	
      do i=1,mb
         d(i)=th(i)
         rho(i)=den(i)
      enddo

      if (mb.gt.nlay) then
	 write(0,'(a)') 'Number of layers exceeds maxium #'
	 stop 
      endif
      hs = 0.
c      write(0,'(a)') ' layer thick   Vp     Vs     rho   Qa    Qb'
      do i = 1, mb
c         read(*,*)d(i),a(i),b(i),rho(i),qa(i),qb(i)
c	 write(0,'(i4,4f7.2,2e9.2)')i,d(i),a(i),b(i),rho(i),qa(i),qb(i)
	 if (i .LT. srclayer) hs = hs + d(i)
      enddo
      vs = b(srclayer)
      write(0,'(a15,f8.1,a20,i3)') 'source depth =',hs,
     +				   ' km at top of layer',srclayer

c input sigma, numbe of points, dt, taper, and samples before first-arr
c sigma is the small imaginary part of the complex frequency, input in
c 1/time_length. Its role is to damp the trace (at rate of exp(-sigma*t))
c to reduce the wrape-arround.
c      read(*,*) sigma,nfft,dt,taper,tb,smth
        sigma=2
      nfft=2**lnpt
       taper=0.2
      smth=1
      nfft2 = nfft/2
      if (nfft2*smth .GT. nt) then
	 write(0,'(a)') 'Number of points exceeds maximum allowed'
	 stop 
      endif
      dw = pi2/(nfft*dt)
      sigma = sigma*dw/pi2
      wc = nfft2*(1.-taper)
      if (wc .LT. 1) wc=1
      taper = pi/(nfft2-wc+1)

c Input phase velocity window, dk, and kmax at w=0.
c pmin and pmax are in 1/vs.
c dk is in pi/max(x,hs). Since J(kx) oscillates with period 2pi/x at
c large distance, we require dk < 0.5 to guarentee at least 4 samples per
c period.
c kmax is in 1/hs. Because the kernels decay with k at rate of exp(-k*hs) at
c w=0, we require kmax > 10 to make sure we have have summed enough.
c      read(*,*) pmin,pmax,dk,kc
      pmin=0.0
      pmax=1.11
      dk=0.1
      kc=20
      kc = kc/hs
      pmin = pmin/vs
      pmax = pmax/vs

c input distance ranges
c      read(*,*) nx
      if (nx.gt.ndis) then
	 write(0,'(a)') 'Number of ranges exceeds max.'
	 stop 
      endif
      xmax = r_max
      dk = dk*pi/xmax
      if(dk.gt.0.03)dk=0.03
      const = dk/pi2

c***************************************************************
c*************** do wavenumber integration for each frequency
      z = pmax*nfft2*dw/kc
      k = sqrt(z*z+1)
      total = nfft2*(kc/dk)*0.5*(k+log(k+z)/z)
c     write(0,'(a3,f9.5,a8,f9.2,a8,i9)')'dk',dk,'kmax',kc,'N',total
      write(0,1001)'dk =',dk,'kmax =',kc,'pmax =',pmax,'N =',total
1001  format(a6,f9.5,a8,f6.2,a8,f6.4,a8,i9)
      tenPerc = 0.1*total
      count = 0.
      omega = 0.
      kc = kc*kc
      do j=1,nfft2
	 do ix = 1,nx
            do l =1,nCom
               sum(ix,l,j) = 0.
            enddo
	 enddo
      enddo
      write(0,'(a)') ' start F-K computation ...'
      do j=1, nfft2		! start frequency loop
         w = cmplx(omega,-sigma)	! complex frequency
         do i = 1, mb
c            at = clog(w/pi2)/pi + cmplx(0.,0.5)		! A&R, p182
c            ka(i) = w/(a(i)*(1.+at/qa(i)))
c            kb(i) = w/(b(i)*(1.+at/qb(i)))
c changed by pcl, 10/29/2004
            at = cmu_freq(omega,1.0,qa(i))
            ka(i) = w/(a(i)*csqrt(at))
            at = cmu_freq(omega,1.0,qb(i))
            kb(i) = w/(b(i)*csqrt(at))
            ka(i) = ka(i)*ka(i)
            kb(i) = kb(i)*kb(i)
         enddo
	 k = omega*pmin + 0.5*dk
	 n=(sqrt(kc+(pmax*omega)**2)-k)/dk
	 do i=1,n		! start k-loop
            call kernel(k,u)
c           write(0,'(f6.4,2e12.4)') k,real(u(1,1)),real(u(1,2))
	    do ix=1,nx
               z = k*x(ix)
	       call besselFn(z,aj0,aj1,aj2)
c n=0
               sum(ix,1,j) = sum(ix,1,j) + u(1,1)*aj0
               sum(ix,2,j) = sum(ix,2,j) - u(1,2)*aj1
               sum(ix,3,j) = sum(ix,3,j) - u(1,3)*aj1
c n=1
               grt =    (u(2,2)+u(2,3))*aj1/z
               sum(ix,4,j) = sum(ix,4,j) + u(2,1)*aj1
               sum(ix,5,j) = sum(ix,5,j) + u(2,2)*aj0 - grt
               sum(ix,6,j) = sum(ix,6,j) + u(2,3)*aj0 - grt
c n=2
               grt = 2.*(u(3,2)+u(3,3))*aj2/z
               sum(ix,7,j) = sum(ix,7,j) + u(3,1)*aj2
               sum(ix,8,j) = sum(ix,8,j) + u(3,2)*aj1 - grt
               sum(ix,9,j) = sum(ix,9,j) + u(3,3)*aj1 - grt
	    enddo
	    k = k+dk
	    count=count+1
	    if ( mod(count,tenPerc) .EQ. 0) then
	       write(0,'(i4,a6)') int(100.*count/total)+1, '% done'
	    endif
	 enddo			! end of k-loop
	 lowpass = 1.
	 if (j .GT. wc) lowpass = 0.5*(1.+cos((j-wc)*taper))
	 at0 = const*lowpass/(w*w)	! w^2 from mu*(w/k)^2 instead of mu
	 do ix=1,nx
            phi = omega*t0(ix)
            at = at0*cmplx(cos(phi),sin(phi))
            do l=1,nCom
	       sum(ix,l,j) = sum(ix,l,j)*at
            enddo
	 enddo
	 omega = omega + dw
      enddo			! end of freqency loop
      write(0,'(i9,a40)') count,' 100% done, writing files ... '
      
c***************************************************************
c*************** do inverse fourier transform
      dt = dt/smth
      nfft = smth*nfft
      dfac = exp(sigma*dt)
      do ix=1,nx
	 if (nfft2 .EQ. 1) then
	    write(*,'(f5.1,9e11.3)')x(ix),(real(sum(ix,l,1)),l=1,nCom)
	 else
c            iblank = index(fout(ix),' ')
c	    fout(ix)(iblank+1:iblank+1) = char(0)
c            do l=1,nCom
            do l=1,ncom_out
	       do j=1,nfft2
c		  data(j) = sum(ix,l,j)
		  data(j) = sum(ix,new_old(l),j)
	       enddo
	       do j=nfft2+1,nfft/2
		  data(j) = 0.
	       enddo
               call fftr(data,%val(nfft),%val(-dt))
               z = exp(sigma*t0(ix))		! removing damping due to
	       k = 1				! use of sigma. Damping is
               do j=1,nfft/2			! is w.r.t t=0
                  ddr = real(data(j))*z
		  k = k+1
                  z = z*dfac
		  ddi = aimag(data(j))*z
		  k = k+1
		  z = z*dfac
                  tdata(l,j) = cmplx(ddr,ddi)
              enddo
c               fout(ix)(iblank:iblank) = char(47+l)
c               call wrtsac0(fout(ix),%val(dt),%val(nfft),%val(t0(ix)),
c     &			 %val(x(ix)),tdata)
           enddo
            do l=1,ncom_out
	       do j=1,nfft2
		  sum(ix,l,j) = tdata(l,j)*sig(l)
	       enddo
            enddo
	 endif
      enddo
      
      return
      end
