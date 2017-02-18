c.......................................................................
c
c  Program : green.f
c  Purpose : for calculating the synthetic seismograms in a layered
c            medium over a half space with a double-couple source
c
c  Input :
c            file input from 'stat.dat'
c	     epla, eplo - event location in latitude and longitude
c            ite1     - number of stations
c            station  - station name
c	     sla, slo - station location in latitude and longitude
c	     filename - station name with proper extension
c
c            file input from 'filename.mod'
c            nlay - number of layers
c            h - layer thickness
c            vp, vs - P and S wave velocity
c            Qp, Qs - P and S wave attenuation factor
c            dens - medium density
c            nt - number of time points
c	     twin - time window for the synthetic seismogram
c	     t00  - zero trace time ahead of the first arrival
c	     fmax - maximum frequency for complete synthetics
c	     fmax1- asymptotic synthetics between fmax and fmax1
c      	     ang1, ang2 - orietations of the two horizontal components
c                   measured from north
c	     akp  - kappa at the site
c	     wid  - fault length
c	     flen - fault width
c	     tdep - depth to the top of the fault
c	     sdept- hypocenter depth (reset later by simula.f)
c	     sdep - receiver depth
c	     mx,my- subfault division in length and width, respectively
c	     vss  - average shear wave velocity in the fault area
c	     vr   - fault rupture velocity
c	     ntc  - number of top corners of the fault
c	     tlat(i) - latitude of these top corners, i=1,ntc
c	     tlon(i) - longitude of these top corners, i=1,ntc
c                WARNING: the order of these latitudes and longitudes
c                         must be listed to the direction of strike
c	     dip(i) - corresponding dips for each segment between 
c                     these corners, i=1,ntc-1
c	     rak(i) - same as above but for rakes, i=1,ntc-1
c
c  Adjustable parameters
c
c	     slip - source slip (set to unit for Green's function)
c	     area - source area (set to unit for Green's function)
c	     tr   - source rise time
c	     eps  - tolerance of the wavenumber summation
c	     al   - repeated source interval due to discretization
c                   (if al=0, the program will take defaul value)
c	     smin - tolerance of source time function
c	     damp - imagenary frequency control factor (if damp=0,
c                   the program will take defaul value)
c            kmax - maxmum number of wavenumber summation
c
c  Output :
c            output filename 'filename.grn'
c	     wid,flen,tdep,sdep,mx,my,ntc,vss,vr,akp - same as input
c	     dip(i)  - same as input, i=1,ntc-1
c	     rak(i)  - same as input, i=1,ntc-1
c	     str(i)  - strikes of each fault segment defined by their
c                      top corners, i=1,ntc-1
c	     xtobs(i)- station locations along strikes of the segments
c                      in relative to their top left corners, i=1,ntc-1
c	     ytobs(i)- station locations perpendicular to strikes of 
c                      the segments in relative to their top left 
c                      corners, i=1,ntc-1
c	     swid(i) - distances to the leftmost top corner along the 
c                      top fault trace, i=1,ntc
c	     rdep(i) - depthes of the subfault centers to the surface, 
c                      i=1,msub
c	     ts(i)   - direct S-wave travel times from the center of 
c                      the subfault to the receiver, i=1,msub)
c	     tt(i)   - S-wave time delays for computing the composite 
c                      source time functions of the subfaults
c            fi(i)   - azmuthes of the station in relative to the 
c                      centers of subfaults, i=1,msub
c	     nt,twin - same as the input
c	     tmin    - time shift to all the synthetics
c	     xobs, yobs - station locations along and perpendicular
c                      to the strike in relative to the epicenter
c    	     wso     - distance of the hypocenter along the top fault
c                      trace to the leftmost edge of the fault
c	     fso     - distance of the hypocenter to top of the fault 
c	     sdept   - hypocenter depth 
c	     iep     - fault segment where the hypocenter lays on
c            Uw(j,i,n) - 3-component displacement spectra from each
c                      subfault, j=1,.,nt, i=1,2,3, n=1,.,msub=mx*my,
c                      i=1, 2, 3 correspond to the radial, transverse
c                      (clockwise) and vertical (downward) component
c                      of the displacement synthetics
c
c                                           -- Yuehua Zeng --
c                      Sept., 1992, at UNR, last modified March, 1999
c
c	Modified by Aasha Pancha to read in green.in which is generated
c       by Anderson csim_new.f script (modified for inputs required for
c       this version of green.f
c
c  B_green_99v5.f - created by John Anderson, Feb. 1, 2010.  Input
c   files are incompatible with input from earlier versions, because
c   in this version the velocity model is read in by layers rather
c   than by parameter.  Intent is for greater flexibility - potential
c   to build models with potentially several more layers while keeping
c   the input readable.  Max number of layers is given by parameter
c   ml, which was increased to 30.
c.....................................................................
c
	parameter (ml=30,mm=50,mm1=8192,mm2=mm1*2)
	dimension h(ml),vp(ml),Qp(ml),vs(ml),Qs(ml),dens(ml),t0(mm),
     +		  rdep(mm),r(mm),fi(mm),a1(mm),a2(mm),a3(mm),a4(mm),
     +	          a5(mm),a6(mm),ir(mm),i1(mm),i2(mm),sin1(mm),sin2(mm),
     +            cos1(mm),cos2(mm),tt(mm),tlat(10),tlon(10),swid(10),
     +	          xtobs(10),ytobs(10),dip(10),str(10),rak(10),ths(mm),
     +	          thp(mm),rdep1(mm),na(mm)
        dimension av(mm2)
	complex*16 ka(ml),kb(ml),nu(ml),ga(ml),kbeta(ml),z(ml),mu(ml),
     +	       I11(2,2),I12(2,2),I21(2,2),I22(2,2),I31(2,2),I32(2,2),
     +	       J11(2,2),J12(2,2),J21(2,2),J22(2,2),J31(2,2),J32(2,2),
     +	       Ish11,Ish12,Ish21,Ish22,Ish31,Ish32,cvp(ml),cvs(ml),
     +	       Jsh11,Jsh12,Jsh21,Jsh22,Jsh31,Jsh32,exd(ml,2),d3,d4,d5,
     +	       HTd(ml,2,2),HTu(ml,2,2),HRd(ml,2,2),HRu(0:ml,2,2),d6,
     +	       HTdsh(ml),HTush(ml),HRdsh(ml),HRush(0:ml),exu(0:ml,2),
     +	       Sd(2),Su(2),Sdsh,Sush,Es(2),Eu(2,2),Ed(2,2),Eush,Edsh,
     +	       aj0,aj1,sdepn,ak,ak2,dkn,ci,cw,Uw(mm1,3,mm),src(mm1),
     +	       dUr0(mm),Ur0(mm),dUf1(mm),Uf1(mm),dUr1(mm),Ur1(mm),
     +         aux(2,2),ainv(2,2),auxi(2),tmp,r0(mm),one,zero,rdepn(mm),
     +         dS0(4,mm),S0(4,mm),dS1(6,mm),S1(6,mm),c1,c2,c3,c4,c5,c6,
     +	       c7,c8
	complex au(mm2),ts(mm),tp(mm)
	character filename*20, station*20
        real al0
	data slip,area,tr,eps,al,smin,damp,kmax
     +	     /1.0,1.0,0.01,0.001,0.0,0.00001,1.0,50000/
	common/para/one,zero,ga,nu,mu,kbeta,ka,kb,dens,nlay,ns
	common/matrix/HTd,HTu,HRd,HRu,HTdsh,HTush,HRdsh,HRush,exd,exu
	common/asymp/cw,sdep,vs,vp,nr
c
c  read in the source and station name, latitudes and longitudes
	open(40,file='stat.dat',form='formatted',status='unknown')
c	read(40,'(5x,f9.4,1x,f9.4)')epla,eplo
c	read(40,'(i5)')ite1
	read(40,*)epla,eplo
        write(*,*) epla, eplo
	read(40,*)ite1
        write(*,*) ite1
	do ite=1,ite1
c	read(40,'(a6,1x,f12.5,1x,f12.5)')station,sla,slo
	read(40,*)station,sla,slo
	write(*,*)station,sla,slo
	filename=trim(station)//'.mod'
	open(10,file=filename,form='formatted',status='unknown')
	read(10,*)nlay
	write(*,*)nlay
        do 4 i=1,nlay
        read(10,*) h(i),vp(i),Qp(i),vs(i),Qs(i),dens(i)
        write(*,*) h(i),vp(i),Qp(i),vs(i),Qs(i),dens(i)
    4   continue

c old input for velocity model, replaced in _99v5
c	read(10,*)(h(i),i=1,nlay)
c	read(10,*)(vp(i),i=1,nlay)
c	read(10,*)(Qp(i),i=1,nlay)
c	read(10,*)(vs(i),i=1,nlay)
c	read(10,*)(Qs(i),i=1,nlay)
c	read(10,*)(dens(i),i=1,nlay)

	read(10,*)nt,twin,t00,fmax,fmax1,ang1,ang2,akp
	write(*,*)nt,twin,t00,fmax,fmax1,ang1,ang2,akp
	read(10,*)wid,flen,tdep,sdept,sdep
	write(*,*)wid,flen,tdep,sdept,sdep
	read(10,*)mx,my,vss,vr
	write(*,*)mx,my,vss,vr
	read(10,*)ntc
	write(*,*)ntc
	read(10,*)(tlat(i),i=1,ntc)
	write(*,*)(tlat(i),i=1,ntc)
	read(10,*)(tlon(i),i=1,ntc)
	write(*,*)(tlon(i),i=1,ntc)
	read(10,*)(dip(i),i=1,ntc-1)
	write(*,*)(dip(i),i=1,ntc-1)
	read(10,*)(rak(i),i=1,ntc-1)
	write(*,*)(rak(i),i=1,ntc-1)
	close(10)
c
	filename=trim(station)//'.grn'
	open(20,file=filename,form='unformatted',status='unknown')
	filename=trim(station)//'.gra'
	open(21,file=filename,form='formatted',status='unknown')
	filename=trim(station)//'.gts'
	open(22,file=filename,form='formatted',status='unknown')
	write(20)wid,flen,tdep,sdep,mx,my,ntc,vss,vr,akp
	write(20)(dip(i),i=1,ntc-1)
	write(20)(rak(i),i=1,ntc-1)
	write(21,*)wid,flen,tdep,sdep,mx,my,ntc,vss,vr,akp
	write(21,*)(dip(i),i=1,ntc-1)
	write(21,*)(rak(i),i=1,ntc-1)
c
c  set the parameters
        icpt=1
        write(*,*)icpt
	pi=3.141592653589793
	pi2=pi*2.0
	ci=cmplx(0.0,1.0)
	one=cmplx(1.0,0.0)
	zero=cmplx(0.0,0.0)
	er=0.1e-10
	nw=nt/2
	dt=twin/float(nt)
	df=1./twin
	if(damp.le.0.0)damp=1.0
	wi=-damp*pi/twin
	wm=pi2
	nwm=ifix(fmax/df)+1
	nwm1=ifix(fmax1/df)+1
	d2r=pi/180.
	do 6 i=1,ntc-1
	  dip(i)=dip(i)*d2r
	  rak(i)=rak(i)*d2r
  6	continue
        icpt=2
        write(*,*)icpt
c
c  compute the station x-y coordinate values
	swid(1)=0.0
	j=0
	do 8 i=1,ntc-1
	  call distaz(tlat(i),tlon(i),tlat(i+1),tlon(i+1),swid(i+1),
     +		      str(i))
	  swid(i+1)=swid(i+1)+swid(i)
	  call distaz(tlat(i),tlon(i),sla,slo,dist,azm)
	  xtobs(i)=dist*cos(d2r*(azm-str(i)))
	  ytobs(i)=dist*sin(d2r*(azm-str(i)))
	  if(j.eq.1)goto 8
	  if((tlat(i)-epla)*(epla-tlat(i+1)).ge.0.0.and.
     +	        (tlon(i)-eplo)*(eplo-tlon(i+1)).ge.0.0)then
	    iep=i
	    j=1
	  elseif((tlat(i)-epla)*(epla-tlat(i+1)).ge.0.0.or.
     +	        (tlon(i)-eplo)*(eplo-tlon(i+1)).ge.0.0)then
	    iep=i
	  endif
  8	continue
c	print*,' Fault length = ',swid(ntc),wid
	swid(ntc)=wid
	call distaz(epla,eplo,sla,slo,dist,azm)
	xobs=dist*cos(d2r*(azm-str(iep)))
	yobs=dist*sin(d2r*(azm-str(iep)))
	wso=xtobs(iep)-xobs+swid(iep)
	fso=sdept-tdep
	if(abs(dip(iep)-90.0*d2r).gt.0.1)then
	  fso=(ytobs(iep)-yobs)/cos(dip(iep))
	  sdept=fso*sin(dip(iep))+tdep
	endif
        icpt=3
        write(*,*)icpt
c
c  compute the subfaults' center global coordinates, rupture time and 
c  and travel times
	msub=mx*my
	subw=wid/float(mx)
	subl=flen/float(my)
	x0=0.5*subw
	y0=0.5*subl
	k=2
	do 11 i=1,mx
	  x=x0+float(i-1)*subw
	  do 9 j=k,ntc
	    if(x.le.swid(j))then
	      ir(i)=j-1
	      goto 10
	    endif
  9	  continue
	  k=j
 10	  continue
	  ia=i
	  do 11 j=2,my
	    ia=ia+mx
	    ir(ia)=ir(i)
 11	continue
	do 12 i=1,msub
	  j=mod(i-1,mx)
	  k=ir(i)
	  xsub=x0+float(j)*subw-swid(k)
          ajga=float(i-1)/mx
          iyy=floor(ajga)
	  ysub=y0+float(iyy)*subl
	  rdep(i)=ysub*sin(dip(k))+tdep
	  ysub=ysub*cos(dip(k))
	  x=xtobs(k)-xsub
	  y=ytobs(k)-ysub
	  r(i)=sqrt(x*x+y*y)
	  fi(i)=acos(x/r(i))
	  if(y.gt.0.0)fi(i)=-fi(i)
	  fi(i)=str(k)-fi(i)/d2r
	  a1(i)=(fi(i)-str(k))*d2r
	  i1(i)=i
 12	continue
        icpt=4
        write(*,*)icpt
	rmax=0.0
	do 20 i=1,msub
c     do 15 loop commented out by jga July 4, 2012, reinstated May 1, 2013
	  do 15 j=i+1,msub
	    if(rdep(i).gt.rdep(j))then
	      w=rdep(i)
	      rdep(i)=rdep(j)
	      rdep(j)=w
	      w=r(i)
	      r(i)=r(j)
	      r(j)=w
	      w=fi(i)
	      fi(i)=fi(j)
	      fi(j)=w
	      w=a1(i)
	      a1(i)=a1(j)
	      a1(j)=w
	      k=i1(i)
	      i1(i)=i1(j)
	      i1(j)=k
	      k=ir(i)
	      ir(i)=ir(j)
	      ir(j)=k
	    endif
 15	  continue
c
	  na(i1(i))=i
	  rdep1(i)=rdep(i)
	  tt(i)=sqrt(r(i)*r(i)+rdep(i)*rdep(i))/vss
	  if(rmax.lt.r(i))rmax=r(i)
	  t0(i)=tt(i)*vss
          icpt=5
          write(*,*)icpt, i
	  call ray2(nlay,h,vs,Qs,rdep1(i),sdep,r(i),ts(i),ths(i))
	  call ray2(nlay,h,vp,Qp,rdep1(i),sdep,r(i),tp(i),thp(i))
          icpt=6
          write(*,*)icpt, i
c
c  a1 and a2 are changed sign due to Green's function computation at 
c  the receiver
	  k=ir(i)
	  fi1=a1(i)
	  a1(i)= cos(rak(k))*cos(dip(k))*cos(fi1)
     +		-sin(rak(k))*cos(2.0*dip(k))*sin(fi1)
	  a2(i)=-cos(rak(k))*cos(dip(k))*sin(fi1)
     +		-sin(rak(k))*cos(2.0*dip(k))*cos(fi1)
	  a3(i)= 0.5*sin(rak(k))*sin(2.0*dip(k))
	  a4(i)= 0.5*(cos(rak(k))*sin(dip(k))*sin(2.0*fi1)
     +	        -sin(rak(k))*sin(2.0*dip(k))*(sin(fi1))**2)
	  a5(i)=-a3(i)-a4(i)
	  a6(i)= cos(rak(k))*sin(dip(k))*cos(2.0*fi1)
     +	        -0.5*sin(rak(k))*sin(2.0*dip(k))*sin(2.0*fi1)
	  fi3=(ang1-fi(i))*d2r
	  fi2=(ang2-fi(i))*d2r
	  cos1(i)=cos(fi3)
	  sin1(i)=sin(fi3)
	  cos2(i)=cos(fi2)
	  sin2(i)=sin(fi2)
c          icpt=7
c          write(*,*)icpt,i
20	continue
	write(20)(str(i),i=1,ntc-1)
	write(20)(xtobs(i),i=1,ntc-1)
	write(20)(ytobs(i),i=1,ntc-1)
	write(20)(swid(i),i=1,ntc)
	write(20)(rdep(na(i)),i=1,msub)
	write(20)(ts(na(i)),i=1,msub)
	write(21,*)(str(i),i=1,ntc-1)
	write(21,*)(xtobs(i),i=1,ntc-1)
	write(21,*)(ytobs(i),i=1,ntc-1)
	write(21,*)(swid(i),i=1,ntc)
	write(21,*)(rdep(na(i)),i=1,msub)
	write(21,*)(ts(na(i)),i=1,msub)
c	print*,wso,fso
c
c  set repeated source interval
        al0=al
	if(al0.le.0.0)al0=vp(nlay)*twin+rmax
	dk=pi2/al0
c
	vs1=0.0
	dens1=0.0
	do 22 i=1,nlay-1
	   vs1=vs1+vs(i)
	   dens1=dens1+dens(i)
 22	continue
	if(nlay.gt.1)then
	  vs1=vs1/float(nlay-1)
	  dens1=dens1/float(nlay-1)
	else
	  vs1=vs(1)
	  dens1=dens(1)
	endif
	slip4pi=slip*area/(pi*4.0)
	do 25 i=1,nlay
	   dens(i)=dens(i)/dens1
 25	continue
c
c  we assume the station is at the surface or within the first layer
c  above the source.
	ns=1
c
	ia=1
	i1(1)=1
	depth=0.0
	do 35 i=1,nlay-1
	   depth=depth+h(i)
	   if(rdep(i1(ia)).le.depth)then 
	     ir(ia)=i
	     do 30 j=i1(ia),msub
	       if(rdep(j).le.depth)then
	         rdep(j)=rdep(j)-depth+h(i)
		 if(t0(j).lt.100.0)then
		   t0(j)=real(tp(j))-t00
		 else
		   t0(j)=t0(j)/vp(nlay)-t00
		 endif
	       else
		 i2(ia)=j-1
		 ia=ia+1
		 i1(ia)=j
		 goto 35
	       endif
 30	     continue
	     i2(ia)=msub
	     goto 40
	   endif
 35	continue
	ir(ia)=nlay
	do 38 j=i1(ia),msub
	  rdep(j)=rdep(j)-depth
	  t0(j)=t0(j)/vp(nlay)-t00
 38	continue
	i2(ia)=msub
 40	continue
c
c  An experiment in v11: undo some of the prior work and set t0(j)=0.
        do 41 j=1,msub
        t0(j)=0.0
 41     continue
c  End of experiment. Delete last three lines to undo this test.
        icpt=8
        write(*,*)icpt
c
c  loop for frequency response
	k1=100
	f=-df
	smax=0.0
	do 180 iw=1,nwm1
	f=f+df
	cw=cmplx(f*pi2,wi)
c
c  ramp time function
	src(iw)=(exp(-ci*cw*tr)-1)/(cw*cw*tr)
	if(abs(src(iw)).gt.smax)smax=abs(src(iw))
	if(abs(src(iw))/smax.le.smin)then
	  do 50 j=1,msub
	  do 50 i=1,3
	    uw(iw,i,j)=zero
 50	  continue
	  goto 180
	endif
c
c  compute the complex velocity for P and S waves
	do 55 i=1,nlay
	  cvp(i)=vp(i)/((one-log(cw/wm)/(pi*Qp(i)))
     +	              *cmplx(1.0,-0.5/Qp(i)))
	  cvs(i)=vs(i)/((one-log(cw/wm)/(pi*Qs(i)))
     +	              *cmplx(1.0,-0.5/Qs(i)))
c
	  ka(i)=(vs1/cvp(i))**2
	  kb(i)=(vs1/cvs(i))**2
	  mu(i)=dens(i)/kb(i)
	  z(i)=h(i)*cw/vs1
 55	continue
	if(iw.gt.nwm)goto 170
	do 58 i=1,msub
	  rdepn(i)=rdep(i)*cw/vs1
	  r0(i)=r(i)*cw/vs1
 58	continue
	sdepn=sdep*cw/vs1
	dkn=dk*vs1/cw
c
	do 65 j=1,msub
	  Ur0(j)=zero
	  Ur1(j)=zero
	  Uf1(j)=zero
	  do 60 i=1,4
	    S0(i,j)=zero
	    S1(i,j)=zero
 60	  continue
	  S1(5,j)=zero
	  S1(6,j)=zero
 65	continue
c
c  loop for the wavenumber sum
	ak=-dkn
	do 165 k=1,kmax
	ak=ak+dkn
	ak2=ak*ak
	do 70 j=1,nlay	
	   nu(j)=cdsqrt(ak2-ka(j))
	   ga(j)=cdsqrt(ak2-kb(j))
	   if(real(nu(j)).lt.0.)nu(j)=-nu(j)
	   if(real(ga(j)).lt.0.)ga(j)=-ga(j)
	   kbeta(j)=2.0*ak2-kb(j)
	   if(j.ne.ns)then
	     exd(j,1)=exp(-nu(j)*z(j))
	     exd(j,2)=exp(-ga(j)*z(j))
	     exu(j-1,1)=exd(j,1)
	     exu(j-1,2)=exd(j,2)
	   else
	     exd(j,1)=exp(-nu(j)*(z(j)-sdepn))
	     exd(j,2)=exp(-ga(j)*(z(j)-sdepn))
	     exu(j-1,1)=exp(-nu(j)*sdepn)
	     exu(j-1,2)=exp(-ga(j)*sdepn)
	   endif
 70	continue
	call genrefl(ak)
c
c  calculate the source field and propagate them to the receiver
	do 75 i=1,2
	  do 75 j=1,2
	    aux(i,j)=zero
	    do 75 l=1,2
	       aux(i,j)=aux(i,j)-exu(ns-1,i)*HRu(ns-1,i,l)
     +	               *exd(ns,l)*HRd(ns,l,j)
 75	continue
	call cinv2(aux,ainv)
c
	do 90 n=1,2
	  if(n.eq.1)then
	    Sd(1)=1.0/dens(ns)
  	    Sd(2)=ak/(dens(ns)*ga(ns))
	    Su(1)=-1.0/dens(ns)
  	    Su(2)=-ak/(dens(ns)*ga(ns))
	  else
	    Sd(1)=ak/(dens(ns)*nu(ns))
  	    Sd(2)=1.0/dens(ns)
	    Su(1)=ak/(dens(ns)*nu(ns))
  	    Su(2)=1.0/dens(ns)
	  endif
	  do 80 i=1,2
	    auxi(i)=Sd(i)
	    do 80 j=1,2
	      auxi(i)=auxi(i)+exu(ns-1,i)*HRu(ns-1,i,j)*Su(j)
 80	  continue
	  do 85 i=1,2
	    Es(i)=zero
	    do 85 j=1,2
	      Es(i)=Es(i)+ainv(i,j)*auxi(j)
 85	  continue
	  Ed(n,1)=Es(1)
	  Ed(n,2)=Es(2)
c
 90	continue
c
  	Sdsh=one/(mu(ns)*ga(ns))
  	Sush=Sdsh
	Edsh=(Sdsh+Sush*exu(ns-1,2)*HRush(ns-1))
     +	    /(one-exd(ns,2)*exu(ns-1,2)*HRdsh(ns)*HRush(ns-1))
c
c  loop for propagating the Green's function source field to receivers 
	i=0
	ns1=ns
	do 158 n=1,ia
c
	  nr=ir(n)
	  do 120 m=1,2
	    Es(1)=Ed(m,1)
	    Es(2)=Ed(m,2)
	    do 100 ii=ns1,nr-1
	      do 95 j=1,2
	        Ed(m,j)=zero
	        do 95 l=1,2
	           Ed(m,j)=Ed(m,j)+HTd(ii,j,l)*Es(l)
 95	      continue
	      Es(1)=Ed(m,1)
	      Es(2)=Ed(m,2)
100	    continue
	    Ed(m,1)=Es(1)
	    Ed(m,2)=Es(2)
	    do 110 ii=1,2
	      Eu(m,ii)=zero
	      do 110 j=1,2
	         Eu(m,ii)=Eu(m,ii)+HRd(nr,ii,j)*Ed(m,j)
110	    continue
120	  continue
c
	  do 145 ii=ns1,nr-1
	    Edsh=HTdsh(ii)*Edsh
145	  continue
	  Eush=HRdsh(nr)*Edsh
c
c  p-sv wave 
	  c1=ak
	  c2=ga(nr)
	  c3=nu(nr)
	  c4=2.0*ak*nu(nr)*mu(nr)
	  c5=kbeta(nr)*mu(nr)
	  c6=2.0*ak*ga(nr)*mu(nr)
	  c7=mu(nr)*(2.0*ka(nr)-kb(nr)-ak2)
	  c8=c7-ak2*mu(nr)
	  c7=c7+c7
c
	  I11(1,1)=-c1
	  I11(1,2)=c2
	  I11(2,1)=-c3
	  I11(2,2)=c1
	  I12(1,1)=-c1
	  I12(1,2)=c2
	  I12(2,1)=c3
	  I12(2,2)=-c1
	  I21(1,1)=c4
	  I21(1,2)=-c5
	  I21(2,1)=c5
	  I21(2,2)=-c6
	  I22(1,1)=-c4
	  I22(1,2)=c5
	  I22(2,1)=c5
	  I22(2,2)=-c6
	  I31(1,1)=c7
	  I31(1,2)=c6
	  I31(2,1)=c8
	  I31(2,2)=c6
	  I32(1,1)=c7
	  I32(1,2)=c6
	  I32(2,1)=c8
	  I32(2,2)=c6
c
c  sh wave 
	  c6=ga(nr)*mu(nr)
	  Ish11=one
	  Ish12=one
	  Ish21=-c6
	  Ish22=c6
	  c6=ak*mu(nr)
	  Ish31=c6
	  Ish32=c6
c
	  do 155 m=i1(n),i2(n)
	    i=i+1
	    if(nr.eq.ns1)then
	      d3=exp(-nu(nr)*(rdepn(i)-sdepn))
	      d4=exp(-ga(nr)*(rdepn(i)-sdepn))
	    else
	      d3=exp(-nu(nr)*rdepn(i))
	      d4=exp(-ga(nr)*rdepn(i))
	    endif
	    d5=exp(-nu(nr)*(z(nr)-rdepn(i)))
	    d6=exp(-ga(nr)*(z(nr)-rdepn(i)))
c
	    J11(1,1)=I11(1,1)*d3
	    J11(1,2)=I11(1,2)*d4
	    J11(2,1)=I11(2,1)*d3
	    J11(2,2)=I11(2,2)*d4
	    J12(1,1)=I12(1,1)*d5
	    J12(1,2)=I12(1,2)*d6
	    J12(2,1)=I12(2,1)*d5
	    J12(2,2)=I12(2,2)*d6
	    J21(1,1)=I21(1,1)*d3
	    J21(1,2)=I21(1,2)*d4
	    J21(2,1)=I21(2,1)*d3
	    J21(2,2)=I21(2,2)*d4
	    J22(1,1)=I22(1,1)*d5
	    J22(1,2)=I22(1,2)*d6
	    J22(2,1)=I22(2,1)*d5
	    J22(2,2)=I22(2,2)*d6
	    J31(1,1)=I31(1,1)*d3
	    J31(1,2)=I31(1,2)*d4
	    J31(2,1)=I31(2,1)*d3
	    J31(2,2)=I31(2,2)*d4
	    J32(1,1)=I32(1,1)*d5
	    J32(1,2)=I32(1,2)*d6
	    J32(2,1)=I32(2,1)*d5
	    J32(2,2)=I32(2,2)*d6
c
	    Jsh11=Ish11*d4
	    Jsh12=Ish12*d6
	    Jsh21=Ish21*d4
	    Jsh22=Ish22*d6
	    Jsh31=Ish31*d4
	    Jsh32=Ish32*d6
c
	    dUr0(i)=zero
            dS0(1,i)=zero
            dS0(2,i)=zero
            dS0(4,i)=zero
            dS0(3,i)=zero
	    dUr1(i)=zero
            dS1(1,i)=zero
            dS1(3,i)=zero
            dS1(5,i)=zero
            dS1(4,i)=zero
	    do 150 j=1,2
	      dUr0(i)=dUr0(i)+J11(1,j)*Ed(1,j)+J12(1,j)*Eu(1,j)
              dS0(1,i)=dS0(1,i)+J21(1,j)*Ed(1,j)+J22(1,j)*Eu(1,j)
              dS0(2,i)=dS0(2,i)+J21(2,j)*Ed(1,j)+J22(2,j)*Eu(1,j)
              dS0(4,i)=dS0(4,i)+J31(1,j)*Ed(1,j)+J32(1,j)*Eu(1,j)
              dS0(3,i)=dS0(3,i)+J31(2,j)*Ed(1,j)+J32(2,j)*Eu(1,j)
	      dUr1(i)=dUr1(i)+J11(1,j)*Ed(2,j)+J12(1,j)*Eu(2,j)
              dS1(1,i)=dS1(1,i)+J21(1,j)*Ed(2,j)+J22(1,j)*Eu(2,j)
              dS1(3,i)=dS1(3,i)+J21(2,j)*Ed(2,j)+J22(2,j)*Eu(2,j)
              dS1(5,i)=dS1(5,i)+J31(1,j)*Ed(2,j)+J32(1,j)*Eu(2,j)
              dS1(4,i)=dS1(4,i)+J31(2,j)*Ed(2,j)+J32(2,j)*Eu(2,j)
150	    continue
c
	    dUf1(i)=Jsh11*Edsh+Jsh12*Eush
	    dS1(2,i)=Jsh21*Edsh+Jsh22*Eush
	    dS1(6,i)=Jsh31*Edsh+Jsh32*Eush
c
c  compute the discrete wave-number sum for the single force
c  displacement and stress
	    rk=real(r0(i)*ak)
	    call bessel(rk,20.0,20,besj0,besj1)
	    aj0=besj0*ak
	    aj1=besj1*ak
	    dUr0(i)=dUr0(i)*aj1
	    dS0(1,i)=dS0(1,i)*aj1
	    dS0(2,i)=dS0(2,i)*aj0
	    dS0(3,i)=dS0(3,i)*aj0
	    dS0(4,i)=dS0(4,i)*aj0
c
	    tmp=(dUr1(i)+dUf1(i))/r0(i)*besj1
	    dUr1(i)=-dUr1(i)*aj0+tmp
	    dUf1(i)=-dUf1(i)*aj0+tmp
	    tmp=(dS1(1,i)+dS1(2,i))/r0(i)*besj1
	    dS1(1,i)=-dS1(1,i)*aj0+tmp
	    dS1(2,i)=-dS1(2,i)*aj0+tmp
	    dS1(3,i)=dS1(3,i)*aj1
	    dS1(4,i)=dS1(4,i)*aj1
	    dS1(5,i)=dS1(5,i)*aj1
	    dS1(6,i)=dS1(6,i)*aj1
c
	    Ur0(i)=Ur0(i)+dUr0(i)
	    Ur1(i)=Ur1(i)+dUr1(i)
	    Uf1(i)=Uf1(i)+dUf1(i)
	    do 152 j=1,4
	      S0(j,i)=S0(j,i)+dS0(j,i)
	      S1(j,i)=S1(j,i)+dS1(j,i)
152	    continue
	    S1(5,i)=S1(5,i)+dS1(5,i)
	    S1(6,i)=S1(6,i)+dS1(6,i)
155	  continue
	  ns1=nr
158	continue
	if(k.lt.k1)goto 165
	amax=0.0
	do 160 i=1,msub
	  amax=amax1(amax,abs(dUr0(i))/(abs(Ur0(i))+er),
     +	    abs(dUr1(i))/(abs(Ur1(i))+er),abs(dUf1(i))/(abs(Uf1(i))+er))
160	continue
c
	if(amax.lt.eps)goto 170
165 	continue
c	print*,' Please reset the maximum kmax !!!'
	stop
c
c  generate the displacement at the receiver due to a buried slip
c  source using reciprocity relations
170	i=0
	do 175 n=1,ia
	  nr=ir(n)
	  do 175 m=i1(n),i2(n)
	    i=i+1
	    if(iw.le.nwm)then
	      S0(3,i)=S0(3,i)-2.0*Ur0(i)*mu(nr)/r0(i)
	      S0(4,i)=S0(4,i)-S0(3,i)
	      tmp=2.0*(Ur1(i)+Uf1(i))*mu(nr)/r0(i)
	      S1(4,i)=S1(4,i)-tmp
	      S1(5,i)=S1(5,i)-S1(4,i)
	      S1(6,i)=S1(6,i)-tmp
	      tmp=cw/vs1*slip4pi*dk*exp(ci*cw*t0(i))
	    else
	      call asympt(cvs(nr),cvp(nr),dens,ts(i),tp(i),ths(i),
     +		          thp(i),r(i),rdep1(i),S0(1,i),S1(1,i))
	      tmp=slip4pi*exp(ci*cw*t0(i))
	    endif
c
	    aj0=-(a1(i)*S1(1,i)+a3(i)*S1(3,i)+a4(i)*S1(4,i)
     +		 +a5(i)*S1(5,i))*tmp
	    aj1= (a2(i)*S1(2,i)+a6(i)*S1(6,i))*tmp
	    Uw(iw,1,i)=aj0*cos1(i)+aj1*sin1(i)
	    Uw(iw,2,i)=aj0*cos2(i)+aj1*sin2(i)
	    Uw(iw,3,i)= (a1(i)*S0(1,i)+a3(i)*S0(2,i)+a4(i)*S0(3,i)
     +		       +a5(i)*S0(4,i))*tmp
175	continue
c
	k1=k
c	print*,'nw,iw,k',nw,iw,k
180 	continue
c
c  compute the inverse FFT for the displacement
c  modified calculation of tmin in 99v7
	tmin=t0(1)
	do 200 n=2,msub
	  if(tmin.gt.t0(n))tmin=t0(n)
 200	continue
c
	print*,(t0(n),n=1,msub),tmin
	dw=pi*df
	nt2=nt+nt
	k=nt+2
	write(20)(tt(na(i)),fi(na(i)),i=1,msub),nt,tmin,twin,xobs,yobs,
     +		  wso,fso,sdept,iep
	write(21,*)(tt(na(i)),fi(na(i)),i=1,msub),nt,tmin,twin,xobs,yobs,
     +		  wso,fso,sdept,iep
	do 240 n=1,msub
	do 230 i=1,3
	do 210 j=1,nwm1
	  au(j)=Uw(j,i,na(n))*src(j)*df
	  au(k-j)=conjg(au(j))
210	continue
	do 212 j=nwm1+1,nw
	  au(j)=zero
	  au(k-j)=zero
212	continue
	au(nw+1)=zero
	call fork(nt,au,1.)
	t=0.0
	do 215 j=1,nt
	  t=t+dt
	  au(j)=cmplx(real(au(j))*exp(-wi*t),0.0)
          av(j)=real(au(j))
215	continue
	do 220 j=nt+1,nt2
	  au(j)=zero
          av(j)=zero
220	continue
        write(22,*)(av(j),j=1,nt)
	call fork(nt2,au,-1.)
	w=-dw
	do 225 j=1,nt
	  w=w+dw
	  au(j)=au(j)*exp(-ci*w*(t0(n)-tmin))
225	continue
	write(20)(au(j),j=1,nt)
	write(21,*)(au(j),j=1,nt)
230	continue
240	continue
c
	close(20)
        close(21)
        close(22)
	end do
	close(40)
	stop
	end
c
c.........................................................
c
	subroutine genrefl(ak)
c
c  Purpose : compute the generalized reflection and transmission
c            coefficients
c
c           Yuehua Zeng, Sept. 1992 at UNR, modified May, 1994
c
	parameter (ml=30)
	real*4 dens(ml)
	complex*16 ka(ml),kb(ml),nu(ml),ga(ml),kbeta(ml),exu(0:ml,2),
     +	       Td(ml,2,2),Tu(ml,2,2),Rd(ml,2,2),Ru(ml,2,2),exd(ml,2),
     +	       HTd(ml,2,2),HTu(ml,2,2),HRd(ml,2,2),HRu(0:ml,2,2),
     +	       Tdsh(ml),Tush(ml),Rdsh(ml),Rush(ml),ainv(2,2),det,
     +	       HTdsh(ml),HTush(ml),HRdsh(ml),HRush(0:ml),aux(2,2),
     +	       a1,a2,a3,a4,a5,a6,a7,a8,a,b,c,d,E1,E2,F1,F2,G1,G2,H1,
     +	       H2,ak,ak2,one,zero,mu(ml)
        common/para/one,zero,ga,nu,mu,kbeta,ka,kb,dens,nlay,ns
	common/matrix/HTd,HTu,HRd,HRu,HTdsh,HTush,HRdsh,HRush,exd,exu
c	
c  motion-stress matrices
	ak2=ak*ak
	a=kbeta(1)*kbeta(1)
	b=4.0*ak2*nu(1)*ga(1)
	d=a-b
	c=4.0*ak*kbeta(1)/d
c
c reflection up at free surface
c
c  p-sv wave
        HRu(0,1,1)=-(a+b)/d*exu(0,1)
        HRu(0,1,2)= c*ga(1)*exu(0,2)
        HRu(0,2,1)=-c*nu(1)*exu(0,1)
        HRu(0,2,2)= (a+b)/d*exu(0,2)
c
c  sh wave
        HRush(0)=exu(0,2)
c
c  modified reflection and transmission coefficients
	a1=dens(1)*nu(1)
	a2=dens(1)*ga(1)
	a6=mu(1)*ga(1)
	do 10 j=1,nlay-1
c
c  p-sv wave
	  a=mu(j+1)*kbeta(j+1)-mu(j)*kbeta(j)
	  b=mu(j+1)*kbeta(j+1)-2.*mu(j)*ak2
	  c=mu(j)*kbeta(j)-2.*mu(j+1)*ak2
	  d=2.*(mu(j+1)-mu(j))
	  E1=b*nu(j)+c*nu(j+1)
	  E2=b*nu(j)-c*nu(j+1)
	  F1=b*ga(j)+c*ga(j+1)
	  F2=b*ga(j)-c*ga(j+1)
	  G1=a+d*nu(j)*ga(j+1)
	  G2=a-d*nu(j)*ga(j+1)
	  H1=a+d*nu(j+1)*ga(j)
	  H2=a-d*nu(j+1)*ga(j)
	  det=-E1*F1+G2*H2*ak2
	  a3=-ak*(a*c+b*d*nu(j)*ga(j))
	  a4=ak*(a*b+c*d*nu(j+1)*ga(j+1))
c
	  Td(j,1,1)=2.*a1*F1/det*exd(j,1)
	  Td(j,1,2)=-2.*ak*a2*G2/det*exd(j,2)
	  Td(j,2,1)=-2.*ak*a1*H2/det*exd(j,1)
	  Td(j,2,2)=2.*a2*E1/det*exd(j,2)
	  Ru(j,1,1)=(E2*F1-G2*H1*ak2)/det*exu(j,1)
	  Ru(j,1,2)=2.*ga(j+1)*a3/det*exu(j,2)
	  Ru(j,2,1)=-2.*nu(j+1)*a3/det*exu(j,1)
	  Ru(j,2,2)=(-F2*E1+H2*G1*ak2)/det*exu(j,2)
	  Rd(j,1,1)=(-E2*F1-G1*H2*ak2)/det*exd(j,1)
	  Rd(j,1,2)=2.*ga(j)*a4/det*exd(j,2)
	  Rd(j,2,1)=-2.*nu(j)*a4/det*exd(j,1)
	  Rd(j,2,2)=(F2*E1+H1*G2*ak2)/det*exd(j,2)
	  a1=dens(j+1)*nu(j+1)
	  a2=dens(j+1)*ga(j+1)
	  Tu(j,1,1)=2.*a1*F1/det*exu(j,1)
	  Tu(j,1,2)=2.*ak*a2*H2/det*exu(j,2)
	  Tu(j,2,1)=2.*ak*a1*G2/det*exu(j,1)
	  Tu(j,2,2)=2.*a2*E1/det*exu(j,2)
c
c  sh wave
	  a5=a6
	  a6=mu(j+1)*ga(j+1)
	  det=a5+a6
	  a7=a5/det
	  a8=a6/det
	  Tdsh(j)=2.*a7*exd(j,2)
	  Rush(j)=(a8-a7)*exu(j,2)
	  Rdsh(j)=(a7-a8)*exd(j,2)
	  Tush(j)=2.*a8*exu(j,2)
 10	continue
c
c  generalized reflection and transmission coefficients for
c  layers above the source.
c
c  p-sv wave
	do 60 n=1,ns-1
	  do 20 i=1,2
	  do 20 j=1,2
	    aux(i,j)=zero
	    do 20 k=1,2
	       aux(i,j)=aux(i,j)-Rd(n,i,k)*HRu(n-1,k,j)
 20	  continue
	  call cinv2(aux,ainv)
	  do 30 i=1,2
	    do 30 j=1,2
	      HTu(n,i,j)=zero
	      do 30 k=1,2
	        HTu(n,i,j)=HTu(n,i,j)+ainv(i,k)*Tu(n,k,j)
 30	  continue
	  do 50 i=1,2
	  do 50 j=1,2
	    HRu(n,i,j)=Ru(n,i,j) 
	    do 40 k=1,2
	    do 40 l=1,2
	       HRu(n,i,j)=HRu(n,i,j)+Td(n,i,k)*HRu(n-1,k,l)*HTu(n,l,j)
 40	    continue
 50	  continue
c
c  sh wave
 	  HTush(n)=Tush(n)/(1-Rdsh(n)*HRush(n-1))
	  HRush(n)=Rush(n)+Tdsh(n)*HRush(n-1)*HTush(n)
c
 60	continue
c
c  generalized reflection and transmission coefficients for
c  layers below the source.
	HRd(nlay,1,1)=zero
	HRd(nlay,1,2)=zero
	HRd(nlay,2,1)=zero
	HRd(nlay,2,2)=zero
	HRdsh(nlay)=zero
c
	do 100 n=nlay-1,ns,-1
c
c  p-sv wave
	do 70 i=1,2
	do 70 j=1,2
	  aux(i,j)=zero
	  do 70 k=1,2
	    aux(i,j)=aux(i,j)-Ru(n,i,k)*HRd(n+1,k,j)
 70	continue
	call cinv2(aux,ainv)
	do 80 i=1,2
	do 80 j=1,2
	  HTd(n,i,j)=zero
	  do 80 k=1,2
	    HTd(n,i,j)=HTd(n,i,j)+ainv(i,k)*Td(n,k,j)
 80	continue
	do 90 i=1,2
	do 90 j=1,2
	  HRd(n,i,j)=Rd(n,i,j) 
	  do 90 k=1,2
	  do 90 l=1,2
	    HRd(n,i,j)=HRd(n,i,j)+Tu(n,i,k)*HRd(n+1,k,l)*HTd(n,l,j)
 90	continue
c
c  sh wave
 	HTdsh(n)=Tdsh(n)/(1-Rush(n)*HRdsh(n+1))
	HRdsh(n)=Rdsh(n)+Tush(n)*HRdsh(n+1)*HTdsh(n)
c
100	continue
c	
	return
	end
c
c.................................................................
c
        SUBROUTINE BESASY(N,X,NTERMS,PNX,QNX)
C 
C BESASY GENERATES P AND Q FUNCTIONS IN HANKEL ASYMPTOTIC EXP.
C 
        ANN      = 4*N*N
        X8       = .125/X
        PNX      = 1.
        QNX      = X8*(ANN-1.)
        TERM     = QNX
        AI       = 2.
        AJ       = 3.
        DO 100 K = 1,NTERMS
        TERM     = -TERM*X8*(ANN-AJ**2)/AI
        PNX      = PNX + TERM
        AI       = AI + 1.
        AJ       = AJ + 2.
        TERM     =  TERM*X8*(ANN-AJ**2)/AI
        QNX      = QNX + TERM
        AI       = AI + 1.
        AJ       = AJ + 2.
 100    CONTINUE
        RETURN
        END
c
c................................................................
c
	SUBROUTINE BESSEL(X,XBSY,NTERMS,BESJ0,BESJ1)
C 
C        BESSEL FUNCTIONS OF THE FIRST KIND OF ORDER ZERO AND ONE
C 
C FOR ABS(X).LT.8. : CHEBYSHEV POLYNOMIAL EXPANSIONS ARE USED
C     AS GIVEN IN NATIONAL PHYS. LAB. MATH. TABLES, VOL. 5, 1962.
C 
      DATA BR0/.2827844947E8/,BR1/-.6852659891E7/,BR2/.38831312263E6/, 
     1   BR3/-.90578674277E4/,BR4/.108306963E3/,BR5/-.73485335935/, 
     2   BR6/.29212672487E-2/,BR7/-.65050170571E-5/,BR8/.64538018051E-8/
      DATA BS0/.2827844947E8/,BS1/.21695247743E6/,BS2/.70046825147E3/
      DATA BRR0/.98087274959E7/,BRR1/-.11425325721E7/, 
     1   BRR2/.40946213625E5/,BRR3/-.66660119856E3/, 
     2   BRR4/.57575414035E1/,BRR5/-.27904475519E-1/,
     3   BRR6/.73493132111E-4/,BRR7/-.84306821641E-7/
      DATA BSS0/.19617454991E8/,BSS1/.16711673184E6/,
     1   BSS2/.60777258247E3/
C 
C FOR 8.0.LE.ABS(X).LE.XBSY : CHEBYSHEV RATIONAL POLYS ARE USED STILL
C 
      DATA BA0/2.5323420902E2/,BA1/4.2217704118E1/,BA2/5.2443314672E-1/
      DATA BB0/.44884594896E3/,BB1/.75322048579E2/ 
      DATA BC0/-1.2339445551E1/,BC1/-2.7788921059/,BC2/-4.9517399126E-2/
      DATA BD0/.17496878239E3/,BD1/.4100554523E2/
      DATA BAA0/3.5451899975E2/,BAA1/5.5544843021E1/,
     1   BAA2/6.5223084285E-1/
      DATA BBB0/.62836856631E3/,BBB1/.97300094628E2/
      DATA BCC0/4.4822348228E1/,BCC1/9.7348068764/,BCC2/1.7725579145E-1/
      DATA BDD0/.21185478331E3/,BDD1/.46917127629E2/
C 
C FOR ABS(X).GY.XBSY : HANKEL'S ASYMPTOTIC EXPANSIONS ARE USED
C     AS GIVEN IN ABROMOVITZ AND STEGUN.
C 
      DATA PI2/.63661977236/,CPI4/.70710681/
C 
        IF (X .NE. 0.) GO TO 10
C 
C X.EQ.0.
C 
        BESJ0=1.
        BESJ1=0.
        GO TO 100
C 
 10     AX=ABS(X)
        IF (AX.GE.8.) GO TO 20
C 
C AX.LT.8. (USE CHEBYSHEV EXPANSIONS FOR THIS RANGE)
C 
        B=X*X
        BP0=((((BR8*B+BR7)*B+BR6)*B+BR5)*B+BR4)*B
        BP0=(((BP0+BR3)*B+BR2)*B+BR1)*B+BR0
        BP1=((((BRR7*B+BRR6)*B+BRR5)*B+BRR4)*B+BRR3)*B
        BP1=((BP1+BRR2)*B+BRR1)*B+BRR0
        BESJ0=BP0/(((B+BS2)*B+BS1)*B+BS0)
        BESJ1=X*BP1/(((B+BSS2)*B+BSS1)*B+BSS0)
        GO TO 100
C 
 20     IF (AX.GT.XBSY) GO TO 30
C 
C 8.0.LE.AX.LE.XBSY   (USE CHEYSHEV EXPANSIONS FOR THIS RANGE)
C 
        B=(8./X)**2
        BP0=((BA2*B+BA1)*B+BA0)/((B+BB1)*B+BB0)
        BQ0=((BC2*B+BC1)*B+BC0)/(AX*((B+BD1)*B+BD0))
        BP1=((BAA2*B+BAA1)*B+BAA0)/((B+BBB1)*B+BBB0)
        BQ1=((BCC2*B+BCC1)*B+BCC0)/(AX*((B+BDD1)*B+BDD0))
        BCOS=COS(AX)
        BSIN=SIN(AX)
        BSQRT=SQRT(AX)
        BESJ0=(BCOS*(BP0+BQ0)+BSIN*(BP0-BQ0))/BSQRT
        BESJ1=(BCOS*(BQ1-BP1)+BSIN*(BQ1+BP1))/BSQRT
        IF (X.LT.0.) BESJ1=-BESJ1
        GO TO 100
C 
C AX.GT.XBSY   (USE HANKEL'S ASYMP EXPANSIONS)
C 
 30     B=SQRT(PI2/AX)
        BCOS=COS(AX)
        BSIN=SIN(AX)
        BC=CPI4*( BCOS+BSIN)
        BS=CPI4*(-BCOS+BSIN)
        CALL BESASY(0,AX,NTERMS,BP0,BP1)
        BESJ0=B*(BP0*BC-BP1*BS)
        CALL BESASY(1,AX,NTERMS,BP0,BP1)
        BESJ1=B*(BP0*BS+BP1*BC)
        IF (X.LT.0.) BESJ1=-BESJ1
C 
 100    RETURN
	END
c
c..............................................
c
      SUBROUTINE FORK(LX,CX,SIGNI)                                      
      COMPLEX CX(LX),CARG,CEXP,CW,CTEMP                                 
      J=1                                                               
      DO I=1,LX               
        IF(I.LE.J) THEN
          CTEMP=CX(J)            
          CX(J)=CX(I)            
          CX(I)=CTEMP               
        END IF
        M=LX/2                    
 20     IF(J.LE.M)GO TO 30        
        J=J-M                     
        M=M/2                     
        IF(M.GE.1)GO TO 20        
 30     J=J+M                     
      END DO
      L=1                       
 40   ISTEP=2*L                 
      DO M=1,L               
        CARG=(0.,1.)*(3.14159265*SIGNI*FLOAT(M-1))/FLOAT(L)
        CW=CEXP(CARG)             
        DO I=M,LX,ISTEP        
          CTEMP=CW*CX(I+L)          
          CX(I+L)=CX(I)-CTEMP       
          CX(I)=CX(I)+CTEMP         
        END DO
      END DO
      L=ISTEP
      IF(L.LT.LX)GO TO 40
      RETURN
      END
c
c..........................................................
c
        subroutine cinv2(a,ainv)
        complex*16 a(2,2),ainv(2,2),det,one
	one=cmplx(1.0,0.0)
	a(1,1)=one+a(1,1)
	a(2,2)=one+a(2,2)
        det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        ainv(1,1)=a(2,2)/det
        ainv(1,2)=-a(1,2)/det
        ainv(2,1)=-a(2,1)/det
        ainv(2,2)=a(1,1)/det
        return
        end
c
c...............................................
c
      subroutine distaz(epla,eplo,sla,slo,dist,azi)
c
c     subroutine distaz(epla,eplo,sla,slo,del,dist,azi,baz)
c     a subroutine to calculate great circle distances and
c     azimuth between two points on the earth'h surface. the
c     earth is assummed to be an ellipsoid of revolution.
c     this routine is from a program by lynn peseckis, 1979.
c
c     input parameters :
c
c     epla, eplo ........latitude and longitude of first point
c                        on earth's surface. north latitude
c                        and east longitude is positive, south
c                        latitude and west longitude is negative.
c     sla, slo ..........latitude and longitude of second point
c                        on earth's surface.
c
c     returned parameters :
c
c     del ...............distance in degrees between two points.
c     dist ..............distance in kilometers between two points.
c     az ................azimuth from first point to second.
c     baz ...............back azimuth from second point to first.
c
c                                  Compiled by Yuehua Zeng at USC
c
      real*8 c,bp,ap,gamma,cbp,sbp,cap,sap,abtem,altem,baztem
c
c  basic constants :
      dr=57.295780
      re=6371.0
      geocen = 0.993305458
c
c  begin calculation of distance
      dlon=abs(eplo-slo)
      if(abs(epla)-90.0.eq.0.0)then
        bp=90.0-epla
        ap=90.0-(atan(geocen*tan(sla/dr)))*dr
        if(epla) 170,170,180
      else
        bp=90.0-(atan(geocen*tan(epla/dr)))*dr
      endif
      if(abs(sla)-90.0.eq.0.0)then
        ap=90.0-sla
        if(sla) 180,170,170
      else
        ap=90.0-(atan(geocen*tan(sla/dr)))*dr
      endif
      if(dlon-0.001.gt.0.0)goto 200
      if(epla-sla)170,170,180
170   azi=0.0
      baz=180.0
      goto 190
180   azi=180.0
      baz=0.0
190   del=abs(ap-bp)
      goto 210
200   gamma=(slo-eplo)/dr
      cbp=cos(bp/dr)
      sbp=sin(bp/dr)
      cap=cos(ap/dr)
      sap=sin(ap/dr)
      abtem=cbp*cap+sbp*sap*cos(gamma)
      c=atan(sqrt(1.0-abtem*abtem)/abtem)
      if(abtem.lt.0.0)c=180.0/dr+c
      del=c*dr
      altem=(cap-abtem*cbp)/(sin(c)*sbp)
      azi=atan(sqrt(1.0-altem*altem)/altem)*dr
      if(altem.lt.0.0)azi=180.0+azi
      baztem=(cbp-abtem*cap)/(sin(c)*sap)
      baz=atan(sqrt(1.0-baztem*baztem)/baztem)*dr
      if(baztem.lt.0.0)baz=180.+baz
      if(sin(gamma).lt.0.0)then
        azi=360.-azi
      else
        baz=360.-baz
      endif
210   dist = del*111.18
c
      return
      end
c
c.............................................................
c
c program: ray.f
c purpose: computing two-point ray tracing in a layered medium
c                                Yuehua Zeng, May, 1995 at UNR
c.............................................................
c
	subroutine ray(nlay,h1,vs,Qs,rdep,sdep,r,ct,theta)
	dimension h(30),h1(nlay),vs(nlay),Qs(nlay)
	double precision p,p1,p2,vsp,aux,dp
	complex ct,ct1
c
	pi=3.1415926
	err=0.001
	dp=0.1
        icpt=101
        write(*,*)icpt,rdep,sdep,r
c
c  initial setting for the ray tracing
	h(1)=h1(1)-sdep
	do 5 i=2,nlay
	  h(i)=h1(i)
  5	continue
	depth=rdep-sdep
	do 10 i=1,nlay
	  depth=depth-h(i)
	  if(depth.le.0.0) goto 15
 10	continue
 15	i1=i
	h(i1)=depth+h1(i1)
	p=r/(sqrt(r*r+(rdep-sdep)**2)*vs(i1))
	perr=err/(sqrt(r*r+(rdep-sdep)**2)*vs(i1))
c
        icpt=102
        write(*,*)icpt,p,perr,i1
c
	if(p.lt.perr)goto 60
c
c  iteration for the two points ray tracing
	ite=-1
	r1=500.0
	r2=0.0
c
c  Newton's method
 20	rr=0.0
 	rdp=0.0
	do 25 i=1,i1
	  vsp=vs(i)*p
	  if(vsp.gt.1.0)then
	    p=1.0/vs(i1)-dp
	    dp=dp*0.1
	    goto 20
	  end if
	  aux=1.0/dsqrt(1.0d0-vsp*vsp)
	  rdp=rdp+vs(i)*h(i)*aux*aux*aux
	  rr=rr+vsp*h(i)*aux
 25	continue
c
        icpt=103
        write(*,*)icpt,rdp,rr
c
	ite=ite+1
	if(rr.gt.r)then
	  p2=p
	  r2=rr
	else
	  p1=p
	  r1=rr
	end if
c
        icpt=104
        write(*,*)icpt,ite,r,r1,r2
c
	if(ite.gt.10)then
	  if(r2.gt.r.and.r1.le.r)goto 30
	  p=p+sign(0.0001,r-r2)
	  goto 20
	end if
	if(abs(r-rr).lt.err)goto 60
	p=p+(r-rr)/rdp
	goto 20

c
c  bisection method
 30	p=0.5d0*(p1+p2)
c
        icpt=105
        write(*,*)icpt,p
c
 40	rr=0.0
	do 50 i=1,i1
	  vsp=vs(i)*p
	  rr=rr+vsp*h(i)/dsqrt(1.0d0-vsp*vsp)
 50	continue
	ite=ite+1
	delta=r-rr
	if(abs(delta).lt.err)then
	  goto 60
	elseif(ite.gt.1000)then
	  ct1=0.0
	  ct=0.0
	  do 55 i=1,i1
	    vsp=vs(i)*p1
	    ct1=ct1+h(i)/(dsqrt(1.d0-vsp*vsp)*vs(i))*cmplx(1.,-.5/Qs(i))
	    vsp=vs(i)*p2
	    ct=ct+h(i)/(dsqrt(1.d0-vsp*vsp)*vs(i))*cmplx(1.,-.5/Qs(i))
 55	  continue
	  ct=(ct1*r2+ct*r1)/(r1+r2)
	  theta=asin(vs(1)*0.5*(p1+p2))
	  return
	endif
	if(rr.gt.r)then
	  p2=p
	  r2=rr
	else
	  p1=p
	  r1=rr
	end if
	p=0.5d0*(p1+p2)
	goto 40
c
c  compute the final travel time of the ray
 60	ct=0.0
	do 70 i=1,i1
	  vsp=vs(i)*p
	  ct=ct+h(i)/(dsqrt(1.0d0-vsp*vsp)*vs(i))*cmplx(1.,-.5/Qs(i))
 70	continue
	theta=asin(vs(1)*p)
	return
	end
c
c.............................................................
c
c program: ray2.f
c purpose: computing two-point ray tracing in a layered medium
c                                Yuehua Zeng, May, 1995 at UNR
c Extensively modified, John Anderson, May 1, 2013
c On rare occasions the original program failed to converge. This
c version uses an algorithm of successive division of the range of
c ray parameter p that cannot get stuck in a loop.
c This subroutine reproduces the original functionality of 'ray'.
c It however has a fundamental problem in that it traces only the
c direct upgoing ray from the source to the station. Thus for any
c station where the downgoing ray has the first arrival, this does
c not give the arrival time of the first arrival.
c.............................................................
c
        subroutine ray2(nlay,h1,vs,Qs,rdep,sdep,r,ct,theta)
        dimension h(30),h1(nlay),vs(nlay),Qs(nlay)
        double precision p,p1,p2,vsp,aux,dp
        complex ct,ct1
c
        pi=3.1415926
        err=0.001*r
c         Here err is modified from original Zeng value.
c         Accepted ray must be within 0.1% of correct distance, rather
c         than within 1 m of the correct distance.
        dp=0.1
        icpt=101
        write(*,*)icpt,rdep,sdep,r
c
c  initial setting for the ray tracing
        h(1)=h1(1)-sdep
c            the above statement requires the station is in the top layer.
        do 5 i=2,nlay
          h(i)=h1(i)
  5     continue
        depth=rdep-sdep
        do 10 i=1,nlay
          depth=depth-h(i)
          if(depth.le.0.0) goto 15
 10     continue
 15     i1=i
c          i1 is the layer containing the source.
        h(i1)=depth+h1(i1)
c          h(i1) is the thickness of the source layer above the source.
        p=r/(sqrt(r*r+(rdep-sdep)**2)*vs(i1))
c          p here is an initial guess. So long as the source is not in
c            a low velocity layer, p will be smaller than the ray parameter
c            for the upgoing ray from the source to the station.
        pcrit=1.0/vs(i1)
c
        icpt=102
        write(*,*)icpt,p,pcrit,i1
c
c  test if the estimate above for p predicts a too small distance
 16     rr=0.0
        do 17 i=1,i1
          vsp=vs(i)*p
          rr=rr+vsp*h(i)/dsqrt(1.0d0-vsp*vsp)
 17     continue
c
        icpt=103
        write(*,*)icpt,r,rr
c
        if (rr.gt.r) then
           p=p/2
           go to 16
        end if
c
c  iteration for the right value of p
        ite=-1
        p1=p
        p2=pcrit
c
        icpt=104
        write(*,*)icpt,p1,p2
c
c  Test the midpoint of the current range until convergence is reached.
 20     rr=0.0
        p=(p1+p2)/2
        do 25 i=1,i1
          vsp=vs(i)*p
          aux=1.0/dsqrt(1.0d0-vsp*vsp)
          rr=rr+vsp*h(i)*aux
 25     continue
c
        icpt=105
        write(*,*)icpt,r,p1,p2,p,rr
c
        if(abs(r-rr).lt.err)goto 60
c
        ite=ite+1
        eps=(rr-r)/r
        icpt=106
        write(*,*)icpt,ite,eps
        if(rr.gt.r)then
          p2=p
        else
          p1=p
        end if
        if(ite.lt.25)then
          go to 20
        end if
c
        icpt=107
        p=(p1+p2)/2
        write(*,*)icpt,r,p1,p2,p,rr
c
c  compute the final travel time of the ray
 60     ct=0.0
        do 70 i=1,i1
          vsp=vs(i)*p
          ct=ct+h(i)/(dsqrt(1.0d0-vsp*vsp)*vs(i))*cmplx(1.,-.5/Qs(i))
 70     continue
        theta=asin(vs(1)*p)
        icpt=108
        write(*,*)icpt,r,p,ct,theta
        return
        end
c
c
c.................................................................
c
c  program: asymp.f
c  purpose: computing asymptotic solution for the direct waves
c                                 Yuehua Zeng, May, 1995 at UNR
c.................................................................
c
	subroutine asympt(cvs,cvp,dens,ts,tp,fs,fp,r,rdep,S0,S1)
c
c  Purpose : compute the direct wave amplitudes in frequency domain
c
	parameter (ml=30)
	dimension vp(ml),vs(ml),dens(ml)
	complex ts,tp
	complex*16 cvp,cvs,ci,cw,S0(4),S1(6),exa,exb,c1,c2,c3,c4,c5,
     +	       one,cwi,wp,ws
c
	common/asymp/cw,sdep,vs,vp,nr
	ci=cmplx(0.0,1.0)
	one=cmplx(1.0,0.0)
	z=sdep-rdep
	D2=z*z+r*r
	D=sqrt(D2)
	cwi=ci*cw
	c4=cwi*D/cvp
	c5=cwi*D/cvs
	exa=exp(-cwi*tp)
	exb=exp(-cwi*ts)
	c1=cvs*cvs/(D2*D2*cw*cw)*((one+c5)*exb-(one+c4)*exa)
	c2=cvs*cvs/(cvp*cvp)*exa/D2
	c3=exb/D2
	c4=c2*c4
	c5=c3*c5
c
	wp=9.0*c1+4.0*c2-3.0*c3+c4
	ws=-6.0*c1-2.0*c2+3.0*c3+c5
	sites=sqrt(vs(nr)*dens(nr)/(vs(1)*dens(1)))
	sitep=sqrt(vp(nr)*dens(nr)/(vp(1)*dens(1)))
c
c  plane wave half space response
	zd=z/D
	rd=r/D
	vsp=vs(1)*vs(1)/(vp(1)*vp(1))
	rds=0.8*sin(fs)+0.2*rd
	zds=0.8*cos(fs)+0.2*zd
	sin2f=rds*rds
	cos2f=1.0-2.0*sin2f
	c5=4.0*zds*csqrt(cmplx(vsp,0.0)-sin2f)
	c1=c5/(cos2f*cos2f+sin2f*c5)*sites
	c2=2.0*cos2f/(cos2f*cos2f+sin2f*c5)*sites
	vsp=1.0/vsp
	rdp=0.8*sin(fp)+0.2*rd
	zdp=0.8*cos(fp)+0.2*zd
	sin2f=rdp*rdp
	cos2f=vsp-2.0*sin2f
	c5=4.0*zdp*csqrt(cmplx(vsp,0.0)-sin2f)
	c3=2.0*vsp*cos2f/(cos2f*cos2f+sin2f*c5)*sitep
	c4=vsp*c5/(cos2f*cos2f+sin2f*c5)*sitep
c
	zd2=2.0*zd*zd
	rd2=2.0*rd*rd
	rzd=2.0*rd*zd
c
	S0(1)=(rds*(zd2-1.0))*c1*ws-(zdp*rzd)*c3*wp
	S0(2)=(rds*rzd)*c1*ws+(zdp*zd2)*c3*wp
	S0(3)=-(rds*rzd)*c1*ws+(zdp*rd2)*c3*wp
	S0(4)=cmplx(0.0,0.0)
	S1(2)=-2.0*zd*ws
	S1(6)=2.0*rd*ws
	S1(1)=(zds*(zd2-1.0))*c2*ws+(rdp*rzd)*c4*wp
	S1(3)=(zds*rzd)*c2*ws-(rdp*zd2)*c4*wp
	S1(4)=-(zds*rzd)*c2*ws-(rdp*rd2)*c4*wp
	S1(5)=cmplx(0.0,0.0)
c
	return
	end
