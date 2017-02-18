c.....................................................................
c
c  Program : simula.f
c  Purpose : simulate strong motion by convolving the source pulses
c            with the Green's functions from a layered elastic earth 
c            model.
c
c  Reference : Zeng, Anderson and Yu (1994).  A composite source model 
c              for computing realistic synthetic strong ground motions, 
c              Geophy. Res. Lett., 21 725-728
c              Zeng, Anderson and Su (1995).  Subevent rake and random 
c              scattering effects in realistic strong ground motion 
c              simulation, Geophy. Res. Lett., 22 17-20
c  Input :
c            terminal input :
c            stname - station name
c            
c            file input from 'nuclear.in'
c            ite1 - number of realization of the composite sources
c            perw - hypocenter location on the fault in terms of the 
c                   ratio of its length from the left corner
c            perf - hypocenter location on the fault in terms of the 
c                   ratio of its width from the top
c            vr   - fault rupture velocity
c            ibp  - 1 for the Brune pulse shape, 2 for the Sato & Hirasawa pulse shape
c
c            file input from 'simula.in'
c            mfil - order of the butterworth filter
c            flw  - highpass corner
c            fhi  - lowpass corner
c	     ncoda- number of points computed for the coda envelope
c	     m    - the corrsponding fft points
c	     err  - tolerance of the coda envelope
c	     gi   - obsorption coefficient
c	     gs   - scattering coefficient
c            fmax - maximum frequency computed for coda waves
c	     Qf   - frequency dependent Q for the coda waves and it 
c                   replaces gi to allow for constant Q if cn=0
c        			    Q=Qf*f**cn
c	     cn   - see above equation 
c	     seed - random seed for scattering coda waves
c	     nscat- number of scattering wavelets
c	     fl   - highpass corner of the cosine filter
c	     flb  - transition bandwidth of the filter, the transition
c                   band begins at fl-flb and end at fl
c	     akp1 - kappa of the site.  It only applies to coda waves
c            akp2 - kappa of the site.  If akp1=0.0, it will use the
c                   default kappa passed from the Green's function
c            yes  - if yes='y', it computes the scattering waves
c
c            file input from 'compom.in'
c	     wid  - fault length
c	     flen - fault width
c	     depth- depth of the hypocenter (not used in this version)
c	     smo  - earthquake moment in dyn-cm
c	     sdrp1- low limit dynamic stress drop
c	     sdrp2- high limit dynamic stress drop
c	     amu  - average shear rigidity at the fault plane
c	     rmax - maximum subevent radius
c	     rmin - minimum subevent radius
c	     D    - fractal dimansion
c	     seed - random seeds for composite source generation
c	     plim - limit in percentage of the subevent's radius for 
c                   the periphery of the subevent to exit of the fault
c                   plane
c
c            file input from 'stname.mod', see green.f for explanation
c
c            file input from 'stname.grn', see green.f for explanation
c
c            file input from 'rn.dat', random numbers from 0 to 1
c
c  Output :
c            output filename 'stname.sim(realization number)'
c	     C nt     - number of time points
c	     C tmin   - initial time of the synthetics
c	     C dt     - time increment
c	     C stname - station name
c	      T = dt(j-1) + tmin   
c	      di(j,i)- output displacement, j=1,nt,i=1,3
c                       i=1 oriented at ang1 (see green.f) from north
c                       i=2 oriented at ang2 (see green.f) from north
c                       i=3 vertical downward
c	      ve(j,i) - output velocity,     j=1,nt,i=1,3
c	      ac(j,i) - output acceleration, j=1,nt,i=1,3
c
c                    by Yuehua Zeng at UNR, last modified in June 1995
c	Modified by Aasha Pancha Feb 2005 so output is ascii
c       Nine columns: d1 d2 d3 v1 v2 v3 a1 a2 a3
c       reads in multple stn files
c
c       Modified by John Anderson, December 2009, to read a file of
c       random numbers rather than depending on a random number generator.
c
c       v11 modified by John Anderson, September, 2012, to read a composite
c       source model from a file that is generated external to this code. 
c       
c       v13 modified by John Anderson, April 2013, to read a separate line
c       from 'nuclear.in' for each realization.
c
c       v14 modified by John Anderson, February, 2014, to read a flag
c       from 'nuclear.in' to choose between the Brune source pulse or
c       the Sato and Hirasawa source pulse. Thus there is no need to
c       compile separate versions for these two models as has been done
c       previously.
c
c       v15 has a major change from v14, implemented by John Anderson,
c       February, 2014. In v14, the coda is calculated using the
c       subroutine scoda.f, written by Y. Zeng. In v15, the coda is
c       applied using subroutine scoda2.f, written by John Anderson.
c       'scoda.f' uses Zeng's multiple scattering model, whereas
c       'scoda2.f' uses Aki's single scattering model to generate the
c       coda waves. 'scoda.f' generates the coda "on the fly", and the
c       coda is added to the Green's function. In contrast, 'scoda2.f'
c       applies the coda by a convolution with the Green's function,
c       done in the frequency domain, and the coda transfer function is
c       read from a file of pre-calculated models. The pre-calculation
c       is done using a Matlab code, after the model has been set up.
c
c......................................................................
c
	parameter (mm=25,mm1=8200,mm2=8200*2,mm3=40000)
	dimension pulse(mm2,mm),di(mm1*2,3),ve(mm1*2,3),ac(mm1*2,3),
     +		  sr(mm3),sx(mm3),sy(mm3),sm(mm3),tt(mm),fi(mm),sz(mm),
     +	          dip(10),rak(10),str(10),xtobs(10),ytobs(10),p(mm1,3),
     +	          swid(10),sindip(10),cosdip(10)
        dimension h(30),pv(30),qp(30),sv(30),qs(30),den(30)
	common/codaw/ncoda,m,err,gi,gs,fl,flb,rh0,rh1,vs0,vs1,vs,twin,dt
        common/random/iseed,rnt(100000)
	complex u(mm1*4),u1(mm1,3),u2(mm1),sim(mm1,3),zero,ts(mm)
c	character stname*6(300),filename*20,yes*1
        character*20 stname(300)
        character*20 filename
	character*5 tmpext
        character*1 yes
	pi=3.141592653589793
	pi2=pi+pi
        dtr=pi/180.0
	zero=cmplx(0.0,0.0)
c
	open(70,file='c1_rn.dat',form='formatted',status='old')
	read(70,'(5e16.7)')(rnt(i),i=1,100000)
        close(70)        
c
	open(40,file='c2_input_stn',form='formatted',status='old')
c
        icp=1
        print*,'control point: ',icp
c
c Next statement starts do loop over stations.
c
        do 300 nst=1,300
c  read in station recording information
	write(6,'(a)')' input the station name :'
	read(40,*)stname(nst) 
        print*,stname(nst)
	if (stname(nst).EQ.'ENDX00') stop
	open(50,file='nuclear.in',form='formatted',status='old')
	read(50,*)ite1
c
c Next statement starts do loop over composite models
c
	do 280 ite=1,ite1
	print*,ite
c	open(80,file='seed.dat',form='unformatted',status='scratch')
c	open(80,file='seed.dat',form='unformatted',status='unknown')
c        icp=2
c        print*,'control point: ',icp
c	print*,ite
c
c  input the coda parameters
        open(10,file='simula.in',form='formatted',status='old')
	read(10,*)mfil,flw,fhi
	read(10,*)ncoda,m,err,gi,gs,fmax,Qf,cn,seed,nscat,fl,flb,akp1,
     +		  akp2
	read(10,'(a1)')yes
	close(10)
        iseed=mod(int(seed),99991)
c	print*,iseed
c
c  compute the subfaults' center global coordinates, rupture time and 
c  and travel times
	filename=trim(stname(nst))//'.mod'
        print*,filename
	open(10,file=filename,form='formatted',status='old')
c
	read(10,*)nlay,i,j
c        print*,nlay,i,j

	do 2 i=1,nlay
  2	read(10,*) h(i),pv(i),qp(i),sv(i),qs(i),den(i)

c the following statements could have caused errors for runs started
c with csim_v7 or later.  
c	read(10,*)(h(i),i=1,nlay)
c	do 2 j=1,3
c  2	read(10,*)(sv(i),i=1,nlay)
c	read(10,*)(den(i),i=1,nlay)
c	read(10,*)(den(i),i=1,nlay)

	close(10)
c        icp=3
c        print*,'control point: ',icp
	vs0=sv(1)
	rh0=den(1)
c
c  read in the Green's function
	filename=trim(stname(nst))//'.grn'
        print*,filename
	open(10,file=filename,form='unformatted',status='old')
	read(10)wid,flen,tdep,sdep,mx,my,ntc,vs,vr,akp
	read(10)(dip(i),i=1,ntc-1)
	read(10)(rak(i),i=1,ntc-1)
	read(10)(str(i),i=1,ntc-1)
	read(10)(xtobs(i),i=1,ntc-1)
	read(10)(ytobs(i),i=1,ntc-1)
	read(10)(swid(i),i=1,ntc)
c        icp=4
c        print*,'control point: ',icp
	akp=akp*0.5
	if(akp2.ne.0.0)akp=0.5*akp2
	akp1=akp1*0.5-akp
	subw=wid/float(mx)
	subl=flen/float(my)
	t0=sqrt(subw*subw+subl*subl)*0.5/vr
	do 5 i=1,ntc-1
	  sindip(i)=sin(dtr*dip(i))
	  cosdip(i)=cos(dtr*dip(i))
  5	continue
	msub=mx*my
	read(10)(sz(i),i=1,msub)
	read(10)(ts(i),i=1,msub)
	read(10)(tt(i),fi(i),i=1,msub),nt,tmin,twin,xobs,yobs,
     +	         wso,fso,depth,iep
c        if (ite .eq. 1) then
	read(50,*)perw,perf,vr,ibp
c        endif
c        icp=5
c        print*,'control point: ',icp
	wso=wid*perw
	fso=flen*perf
	dt=twin/float(nt)
	nt2=nt*2
	sc=dt/float(nt2)
	dw=pi/twin
	df=1/twin
c
	nk=nt2+2
	do 10 i=1,3
	  do 10 j=1,nt
	    sim(j,i)=zero
 10	continue
c        icp=6
c        print*,'control point: ',icp
c
c  generate the random numbers
	if(yes.eq.'y'.or.yes.eq.'Y')then
	  k=nt+2
	  nw=nt/2
	  nwm=fmax/df+1
	  if(nwm.gt.nw)nwm=nw
	  dw=dw+dw
	  w=-dw+0.00001
	  do 12 j=1,nwm
	    w=w+dw
	    Qt=1.0/(Qf*(w/pi2)**cn)
	    u2(j)=((w-w*log(w/pi2)*Qt/pi)/vs)*cmplx(0.5*Qt,1.0)
 12	  continue
c        icp=6
c        print*,'control point: ',icp
	  do 20 i=1,3
	    do 13 j=1,nt
	      u(j)=zero
 13	    continue
	    do 15 l=1,nscat
	      x=vs*(tmin+rnt(iseed)*twin)
              iseed=iseed+1
              iseed=mod(iseed,99991)
	      c1=rnt(iseed)*2.0-1.0
              iseed=iseed+1
              iseed=mod(iseed,99991)
	      do 14 j=1,nwm
	        u(j)=u(j)+c1*exp(-x*u2(j))
 14	      continue
 15	    continue
c        icp=7
c        print*,'control point: ',icp
	    w=-dw
	    do 17 j=1,nwm
	      w=w+dw
	      u(j)=u(j)*exp(cmplx(-w*akp1,w*tmin))*df
	      u(k-j)=conjg(u(j))
 17	    continue
	    call fork(nt,u,1.)
	    do 18 j=1,nt
	      p(j,i)=real(u(j))
 18	    continue
 20	  continue
	  dw=dw*0.5
	endif
c        icp=8
c        print*,'control point: ',icp
c
c  read in the composite source model
c

        filename='csevents01.dat'
        if (ite.le.9)then
          write(filename(10:10),'(i1)')ite
        else
          write(filename(9:10),'(i2)')ite
        endif
c        print*,filename



	open(20,file=filename,form='formatted',status='old')
	read(20,*)sdrp,amu,realn
	read(20,*)seed,seed1,seed2
        n=ifix(realn)
        print*,sdrp,amu,n
c        print*,seed,seed1,seed2
        iseed=mod(int(seed),99991)
        do 24 i=1,n
        read(20,*) sx(i),sy(i),sr(i)
c        print*,sx(i),sy(i),sr(i)
   24   continue
	close(20)
c        icp=9
c        print*,'control point: ',icp
c
c	sum=0.0
	do 25 i=1,n
	  sm(i)=(16/7)*(sdrp*1.0e6/amu)*1.0e5*sr(i)**3
c	  sum=sum+sm(i)
 25	continue
c	sum=3.0*smo/sum   %as programmed by Zeng for Sato & Hirasawa pulse
c	sum=smo/sum
c        icp=10
c        print*,'control point: ',icp,n
c
c  loop for subfault source pulse functions computation
	do 30 i=1,msub
	  do 30 j=1,nt2
	    pulse(j,i)=0.0
 30	continue
c        icp=11
c        print*,'control point: ',icp
c
	dt=0.5*dt
	do 50 i=1,n
	  nx=ifix((sx(i))/subw)+1
	  ny=ifix((sy(i))/subl)+1
c        icp=12
c        print*,'control point: ',icp,nx,ny
	  do 35 j=2,ntc
	    if(sx(i).le.swid(j))goto 36
 35	  continue
 36	  k=j-1
	  nsub=(ny-1)*mx+nx
c
c  compute the rupture time for the subevent
	  t=sqrt((sx(i)-wso)**2+(sy(i)-fso)**2)/vr
	  sx(i)=sx(i)-swid(k)
c          print*,t
c
c  generate the global coordinates centered at the epicenter and 
c  compute sub-event moment and corner frequency
	  z=sdep-(sy(i)*sindip(k)+tdep)
	  y=ytobs(k)-sy(i)*cosdip(k)
	  x=xtobs(k)-sx(i)
	  r=sqrt(x*x+y*y+z*z)
c        print*,r
c
c  we should use sm(i)=16/7*sr(i)**3*sdrp.  But the following equation 
c  adjusts the result for probabilistic numerical discretization error
c	  sm(i)=sum*sm(i)
c
c  find the arrival time for the subevent and its brune's pulse
	  t=t+r/vs
c
	  j1=ifix((t-tt(nsub)+t0)/dt)
c  following added by John Anderson, Jan 2010.
c  give user the choice of either the Brune pulse, as in original Zeng et al (1994)
c  or the Sato & Hirasawa pulse, as implemented by Zeng later.
c  ibp=1 for Brune pulse
c  ibp=2 for Sato & Hirasawa pulse
c
c  eventually give user an input flag for this choice.  For now, hard-code.
c
c      ibp=1

      if(ibp.eq.2) go to 43
c        icp=13
c        print*,'control point: ',icp,nx,ny,t,j1

   40 continue
      fc=2.34*vs/(2.0*pi*sr(i))
      wcb=2.0*pi*fc
      nn1=j1+ifix(2.0/(dt*fc))+1
      tj=0.0
      do 41 j=j1,nn1
      tj=tj+dt
      pulse(j,nsub)=pulse(j,nsub)+wcb*wcb*sm(i)*tj*exp(-wcb*tj)
   41 continue
      go to 49
c
c  Sato & Hirasawa pulse.  I could not demonstrate that Zeng programmed this correctly.
c  Referring to the Sato and Hirasawa (1973) paper, I have rewritten this part of 
c  the code, using clearer variable names that more closely follow the equations in 
c  the paper.  The mysterious factor of 3.0 used above was appropriate for this pulse. 
c  It has nothing to do with discretization, but comes from the way the Sato & Hirasawa
c  pulse is normalized.  It does not belong to the Brune pulse.  
c    -John Anderson, Jan 26, 2009.
c
   43 continue
        vrs=vr/vs
        ak=vrs*0.707
        cak=1-ak*ak
        pak=1+ak
        rc=sr(i)/vr
        coef=3*vr/sr(i)
        t1=(1-ak)*rc
        t2=(1+ak)*rc
        nn1=j1+ifix(t1/dt)
        nn2=j1+ifix(t2/dt)
        if(nn1.gt.nt2)nn1=nt2
        if(nn2.gt.nt2)nn2=nt2

c          av=(vrs)**2
c	  ak=sqrt(av-ak*ak)
c	  wc=vr/sr(i)
c	  dtwc=dt*wc
c	  wcm=wc*sm(i)

c	  nn1=j1+ifix((1.0-ak)/dtwc)
c	  nn2=j1+ifix((1.0+ak)/dtwc)
c	  if(nn1.gt.nt2)nn1=nt2
c	  if(nn2.gt.nt2)nn2=nt2
c	  ak1=(1.0-ak*ak)**2
c	  ak2=(1.0+ak)**2

c        icp=13
c        print*,'control point: ',icp
      tj=0.0
      do 44 j=j1,nn1
      tj=tj+dt
      xa=vr*tj/sr(i)
      pulse(j,nsub)=pulse(j,nsub)+sm(i)*coef*xa**2/(cak**2)
 44   continue
      do 45 j=nn1+1,nn2
      tj=tj+dt
      xb=vr*tj/sr(i)
      pulse(j,nsub)=pulse(j,nsub)+
     1    sm(i)*coef*0.25*(1/ak - xb**2/(ak*pak**2))
 45   continue
 49   continue
 50   continue
c
c  Note following mysterious statement from original program
	dt=dt+dt
c        icp=14
c        print*,'control point: ',icp,nst
c
	do 110 i=1,msub
c 
c  compute the FFT of the subfault source pulse function with zeros
c  padded from nt2+1 to 2*nt2 for linear convolution
	  do 60 j=1,nt2
	    u(j)=cmplx(pulse(j,i)*0.5,0.0)
 60	  continue
	  do 70 j=nt2+1,nt2+nt2
	    u(j)=zero
 70	  continue
	  call fork(2*nt2,u,-1.)
	  do 80 j=1,nt
	    u2(j)=u(j)
 80	  continue
c
c  loop for summation of the wave field from each subfault generated 
c  by convolving its source pulses with its Green's function of a 
c  layered elastic solid.
	  do 82 k=1,3
	    read(10)(u1(j,k),j=1,nt)
 82	  continue
c
c  simulate strong motion coda waves for each subfaults
	  if(yes.eq.'y'.or.yes.eq.'Y')then
	  depth=0.0
	  do 85 j=1,nlay
	    depth=depth+h(j)
	    if(sz(i).le.depth)then
	      vs1=sv(j)
	      rh1=den(j)
	      goto 86
	    endif
 85	  continue
 86	  dist=tt(i)*vs
	  call scoda2(nst,i,u1)
c	  call scoda(dist,nt,nscat,ts(i),tmin,u1,p)
	  endif
c
c  sum up the wave field
	  do 100 k=1,3
	    do 90 j=1,nt
	      sim(j,k)=sim(j,k)+u1(j,k)*u2(j)
 90	    continue
100	  continue
c
110	continue
	close(10)
	wcm=twin+twin
	call scat1d(nt2,wcm,fmax,u1(1,1))
c
	do 140 k=1,3
c
c  displacement
	  w=-dw
	  do 115 j=1,nt
	    w=w+dw
	    u2(j)=sim(j,k)*u1(j,1)*sc*exp(w*cmplx(-akp,t0))
	    u(j)=u2(j)
	    u(nk-j)=conjg(u(j))
115	  continue
	  u(nt+1)=zero
	  call fork(nt2,u,1.)
	  do 118 j=1,nt2
	    di(j,k)=real(u(j))
118	  continue
c
c  velocity
	  w=-dw
	  do 120 j=1,nt
	    w=w+dw
	    u(j)=u2(j)*cmplx(0.0,w)
	    u(nk-j)=conjg(u(j))
120	  continue
	  u(nt+1)=zero
	  call fork(nt2,u,1.)
	  do 122 j=1,nt2
	    ve(j,k)=real(u(j))
122	  continue
c
c  acceleration
	  w=-dw
	  do 125 j=1,nt
	    w=w+dw
	    u(j)=-u2(j)*w*w
	    u(nk-j)=conjg(u(j))
125	  continue
	  u(nt+1)=zero
	  call fork(nt2,u,1.)
	  do 130 j=1,nt2
	    ac(j,k)=real(u(j))
130	  continue
140	continue
c
c  filtering the seismogram
	do 160 i=1,3
 	  do 150 j=1,nt2
 	    di(j,i)=di(j,i)-di(1,i)
 	    ve(j,i)=ve(j,i)-ve(1,i)
 	    ac(j,i)=ac(j,i)-ac(1,i)
150	  continue
c	  call pbutter(di(1,i),nt2,mfil,dt,flw,fhi)
c	  call pbutter(ve(1,i),nt2,mfil,dt,flw,fhi)
c	  call pbutter(ac(1,i),nt2,mfil,dt,flw,fhi)
160	continue
c
c  dump out the simulation results
	write(tmpext, '(i2.2)')ite
c	filename=stname(nst)//'.sim00'
	filename=trim(stname(nst))//'.sim'//trim(tmpext)
c	if(ite.le.9)then
c	  write(filename(12:12),'(i1)')ite
c	else
c	  write(filename(11:12),'(i2)')ite
c	endif
c	open(10,file=filename,form='unformatted',status='unknown')
 	open(10,file=filename,form='formatted',status='unknown')
CC	write(10, '(i5)')nt
c 	write(10,'(f10.5)')tmin
CC 	write(10,'(f10.5)')dt
c	write(10, '(a4)')stname
c 	write(10, '(1e12.4)')((di(j,i),j=1,nt),i=1,3)
c 	write(10, '(1e12.4)')((ve(j,i),j=1,nt),i=1,3)
c 	write(10, '(1e12.4)')((ac(j,i),j=1,nt),i=1,3)
	do 245 j=1, nt
	T = dt*(j-1) +tmin
	write(10, '(10e12.4)') T, di(j,1), di(j,2), di(j,3), 
     1    ve(j,1), ve(j,2), ve(j,3), ac(j,1), ac(j,2), ac(j,3)	
245	continue

	close(10)
280     continue
c	end do
	close(50)
300     continue
        close(40)
	stop
	end
c
c..............................................
c
      subroutine fork(lx,cx,signi)
      complex cx(lx),carg,cexp,cw,ctemp
      j=1                             
      do i=1,lx                     
        if(i.le.j) then
          ctemp=cx(j)   
          cx(j)=cx(i)   
          cx(i)=ctemp            
        end if
        m=lx/2                  
 20     if(j.le.m)go to 30     
        j=j-m                 
        m=m/2                
        if(m.ge.1)go to 20  
 30     j=j+m              
      end do
      l=1                 
 40   istep=2*l          
      do m=1,l          
        carg=(0.,1.)*(3.14159265*signi*float(m-1))/float(l)
        cw=cexp(carg)  
        do i=m,lx,istep
          ctemp=cw*cx(i+l) 
          cx(i+l)=cx(i)-ctemp
          cx(i)=cx(i)+ctemp 
        end do
      end do
      l=istep
      if(l.lt.lx)go to 40
      return
      end
c
c.............................................................
c
      subroutine pbutter(xr,lx,mfil,dt,fhi,flo)
c
c  a two pass high and low filter
c
      real*4 xr(lx)
      call bwthhi(xr,xr,lx,mfil,dt,fhi)
      call rev(xr,lx)
      call bwthhi(xr,xr,lx,mfil,dt,fhi)
      call rev(xr,lx)
c
c     call bwthlo(xr,xr,lx,mfil,dt,flo)
c     call rev(xr,lx)
c     call bwthlo(xr,xr,lx,mfil,dt,flo)
c     call rev(xr,lx)
      return
      end
      subroutine rev(x,n)
      dimension x(n)
      m=n/2
      do 1 i=1,m
      j=i-1
      x1=x(i)
      x(i)=x(n-j)
    1 x(n-j)=x1
      return
      end
c
      subroutine bwthhi(x,y,n,iord,dt,w3db)
c
c   11 apr 79
c
c   butterworth high-pass filter.    see gold & rader 1969, p.72.
c
c   unity gain
c   input array     x(n)
c   output array    y(n)    (may be the same as x)
c   iord     order of filter   may be 1, 2 or 3.
c   dt       timestep of series
c   w3db     3db point of filter
c
      dimension x(n),y(n)
      data pi/3.1415926/
c
      is=1
      wc=tan(pi*dt*w3db)
c
c   1st order filter
      if(iord.eq.1)then
      a=1./(1.+wc)
      b=-(1.-wc)*a
c
      x1=x(1)
      y(1)=a*x(1)
      y1=y(1)
      is=2
      do 11 i=is,n
      x0=x(i)
      y(i)=a*(x0-x1)-b*y1
      y1=y(i)
      x1=x0
11    continue
c
c   2nd order filter
      elseif(iord.eq.2)then
      wcr2=wc*sqrt(2.)
      wc2=wc*wc
      a=1./(1.+wcr2+wc2)
      b=-2.*(1.-wc2)/(1.+wcr2+wc2)
      c=(1.-wcr2+wc2)/(1.+wcr2+wc2)
c
      x0=x(1)
      y(1)=a*x0
      y1=y(1)
      x1=x0
      x0=x(2)
      y(2)=a*(x0-x1-x1)-b*y1
      y2=y1
      y1=y(2)
      x2=x1
      x1=x0
      is=3
      do 21 i=is,n
      x0=x(i)
      y(i)=a*(x0-x1-x1+x2)-b*y1-c*y2
      y2=y1
      y1=y(i)
      x2=x1
      x1=x0
21    continue
c
c   3rd order filter
      elseif(iord.eq.3)then
      wc2=wc*wc
      wca=1.+wc+wc2
      wcb=1.-wc+wc2
      wcc=1.-wc2
      wcd=wca*(1.+wc)
      a=1./wcd
      b=(wca*(-1.+wc)-2.*(1.+wc)*wcc)/wcd
      c=(wcb*(1.+wc) +2.*(1.-wc)*wcc)/wcd
      d=wcb*(+wc-1.)/wcd
c
      x0=x(1)
      y(1)=a*x0
      y1=y(1)
      x1=x0
      x0=x(2)
      y(2)=a*(x0-x1-x1-x1)-b*y1
      y2=y1
      y1=y(2)
      x2=x1
      x1=x0
      x0=x(3)
      y(3)=a*(x0-x1-x1-x1+x2+x2+x2)-b*y1-c*y2
      y3=y2
      y2=y1
      y1=y(3)
      x3=x2
      x2=x1
      x1=x0
      is=4
      do 31 i=is,n
      x0=x(i)
      y(i)=a*(x0-x1-x1-x1+x2+x2+x2-x3)-b*y1-c*y2-d*y3
      y3=y2
      y2=y1
      y1=y(i)
      x3=x2
      x2=x1
      x1=x0
31    continue
      else
      write(6,'(a)')' Wrong order for the filter!'
      stop
      endif
c
      return
      end
c
      subroutine bwthlo(x,y,n,iord,dt,w3db)
c
c   11 apr 79
c
c   butterworth low-pass filter.    see gold & rader 1969, p.72.
c
c   unity gain
c   input array     x(n)
c   output array    y(n)    (may be the same as x)
c   iord     order of filter   may be 1, 2 or 3.
c   dt       timestep of series
c   w3db     3db point of filter
c
      dimension x(n),y(n)
      data pi/3.1415926/
c
      is=1
      wc=tan(pi*dt*w3db)
c
c   1st order filter
      if(iord.eq.1)then
      a=wc/(wc+1.)
      b=(wc-1.)/(wc+1.)
c
      x1=x(1)
      y(1)=a*x(1)
      y1=y(1)
      is=2
      do 11 i=is,n
      x0=x(i)
      y(i)=a*(x0+x1)-b*y1
      y1=y(i)
      x1=x0
11    continue
c
c   2nd order filter
      elseif(iord.eq.2)then
      wcr2=wc*sqrt(2.)
      wc2=wc*wc
      a=wc2/(1.+wcr2+wc2)
      b=-2.*(1.-wc2)/(1.+wcr2+wc2)
      c=(1.-wcr2+wc2)/(1.+wcr2+wc2)
c
      x0=x(1)
      y(1)=a*x0
      y1=y(1)
      x1=x0
      x0=x(2)
      y(2)=a*(x0+x1+x1)-b*y1
      y2=y1
      y1=y(2)
      x2=x1
      x1=x0
      is=3
      do 21 i=is,n
      x0=x(i)
      y(i)=a*(x0+x1+x1+x2)-b*y1-c*y2
      y2=y1
      y1=y(i)
      x2=x1
      x1=x0
21    continue
c
c   3rd order filter
      elseif(iord.eq.3)then
      wc2=wc*wc
      wc3=wc*wc2
      wca=1.+wc+wc2
      wcb=1.-wc+wc2
      wcc=1.-wc2
      wcd=wca*(wc+1.)
      a=wc3/wcd
      c=(wcb*(1.+wc)+2.*(1.-wc)*wcc)/wcd
      b=(wca*(wc-1.)-2.*(1.+wc)*wcc)/wcd
      d=wcb*(wc-1.)/wcd
c
      x0=x(1)
      y(1)=a*x0
      y1=y(1)
      x1=x0
      x0=x(2)
      y(2)=a*(x0+x1+x1+x1)-b*y1
      y2=y1
      y1=y(2)
      x2=x1
      x1=x0
      x0=x(3)
      y(3)=a*(x0+x1+x1+x1+x2+x2+x2)-b*y1-c*y2
      y3=y2
      y2=y1
      y1=y(3)
      x3=x2
      x2=x1
      x1=x0
      is=4
      do 31 i=is,n
      x0=x(i)
      y(i)=a*(x0+x1+x1+x1+x2+x2+x2+x3)-b*y1-c*y2-d*y3
      y3=y2
      y2=y1
      y1=y(i)
      x3=x2
      x2=x1
      x1=x0
31    continue
      else
      write(6,'(a)')' Wrong order for the filter!'
      stop
      endif
c
      return
      end
c
c..................................................................
c
	subroutine scoda(dist,nt,nscat,ts,t0,dis,pran)
c...................................................................
c  Program : scoda.f
c  Purpose : simulate strong motion S-wave coda
c
c                                by Yuehua Zeng at UNR, Oct 1993
c...................................................................
c
	parameter (mm1=8200,mm2=400,mm3=mm1*2)
	dimension pran(mm1,3),t(mm2),p(mm2)
	complex u1(mm3),dis(mm1,3),p1(mm1),zero,ts
	common/codaw/n,m,err,gi,gs,fl,flb,rh0,rh1,vs0,vs1,vs,twin,dt
	pi=3.1415926
	pi2=3.1415926*2.0
	pi12=3.1415926/2.0
	zero=cmplx(0.0,0.0)
c
c  filter the low frequency
	nt2=nt*2
	df=1.0/(nt2*dt)
	do 10 j=1,nt
	  p1(j)=1.0
 10	continue
	ml=ifix(fl/df)+1
	mlb=ifix(flb/df)
	if(mlb.eq.0)mlb=1
	do 60 j=1,ml-mlb
	  p1(j)=0.0
 60     continue
	p1(1)=0.0
	band=-1.0
	db=1./float(mlb)
	mmb=ml-mlb+1
	if(mmb.le.1)mmb=2
	do 62 j=mmb,ml+mlb-1
	  band=band+db
	  p1(j)=p1(j)*(0.5+0.5*sin(pi12*band))
 62     continue
c
c  generate the coda waves
	call coda(n,m,gs,gi,vs,2.0*twin,dist,err,t,p,dt)
c
	c0=0.85/(pi*4.*dist*vs)*sqrt(twin*vs1*rh1/(vs0*rh0*nscat))
	c1=sqrt(3.0)*c0
	dt1=t(2)-t(1)
	t00=dist/vs-real(ts)+t0
	i1=ifix(t00/dt)
	i2=ifix(t(1)/dt)
c
	nk=nt2+2
	do 110 ii=1,3
	  do 64 i=1,n
	    if(t00.le.t(i))goto 65
 64	  continue
 65	  n1=i
	  i0=0
	  t1=float(i1-1)*dt
	  do 70 i=1,i2-i1+1
	    i0=i0+1
	    t1=t1+dt
	    u1(i)=zero
 70	  continue
	  do 75 i=i0+1,nt
	    t1=t1+dt
	    if(t1.gt.t(n1))n1=n1+1
	    pp=p(n1-1)+(p(n1)-p(n1-1))*(t1-t(n1-1))/dt1
	    u1(i)=cmplx(pp*(c1*pran(i,ii)),0.0)
 75	  continue
	  do 80 i=nt-10,nt2
	    u1(i)=zero
 80	  continue
	  call fork(nt2,u1,-1.0)
	  do 90 i=1,nt
	    u1(i)=u1(i)*p1(i)
 90	  continue
	  do 100 i=1,nt
	    dis(i,ii)=dis(i,ii)+u1(i)
100	  continue
110	continue
c
	return
	end
C
C.......................................................................
C
C  Program for calculating the scattered wave energy using equation
C  (24) of Zeng et al. (1991).  Written by Yuehua Zeng in May 1990
C  and updated by Yuehua Zeng on Feb. 9, 1993.
C.....................................................................
C
	subroutine coda(n,m,gs,gi,vs,tw,r,err,t,p,ddt)
	dimension t(n),p(n)
	complex au(1024),atinv,atinv1,atinv3,akw,wg,one,zero,ei,
     +		ei2,aker,Es,Es1,Es2
c
c  compute parameters
	pi=3.14159265
	one=cmplx(1.0,0.0)
	zero=cmplx(0.0,0.0)
	ei=cmplx(0.0,1.0)
	ei2=0.5/ei
	ep=1.0e-20
	g=gs+gi
	dw=2.0*pi/tw
	nd=2*ifix(0.25*(r+vs*tw)/r+1)
	dk=pi/(r*float(nd))
	ni=ifix(20.0/nd+1)*nd+nd/2
	m2=m/2
	dt=(tw-r/(vs*sqrt(3.0)))/float(n)
	t(1)=r/(vs*sqrt(3.0))
	do 10 i=2,n
	  t(i)=t(i-1)+dt
 10	continue
	i=ifix((r/vs-t(1))/dt)+1
	t(i+1)=r/vs+ddt
	t(i)=r/vs-ddt
c
c  calculate the power spectrum for scattering order higher than two
	call evwi(gi,gs,vs,r,tw,wi)
	print*,'wi =',wi
	w=-dw
	do 30 i=1,m2
	  w=w+dw
	  wg=cmplx(g+wi/vs,w/vs)
	  ak=0.0
	  Es=zero
	  Es1=zero
	  Es2=zero
	  nk=ni
	  ik=0
	  do 20 j=1,10000
	    ak=ak+dk
	    akw=cmplx(0.0,ak)/wg
	    atinv=ei2*clog((one+akw)/(one-akw))
	    atinv1=gs*atinv/ak
	    atinv3=atinv1*atinv1*atinv1
	    aker=atinv3*atinv/(one-atinv1)
	    sinkr=sin(ak*r)
	    Es=Es+sinkr*aker
	    if(cabs(aker/(Es+ep)).lt.err.and.j.gt.20)goto 25
	    if(j.eq.nk)then
	      ik=ik+1
	      Es2=Es2+Es-sinkr*aker*0.5
	      if(ik.eq.6)then
	        if(cabs((Es2-Es1)/(Es2+ep)).lt.err)then
	          Es=(Es1+Es2)/12.0
		  goto 25
		endif
		ik=0
	        Es1=Es2
	        Es2=zero
	      endif
	      nk=nk+nd
	    endif
 20	  continue
 25	  au(i)=Es*dk
 30	continue
c
	n1=m+2
	do 40 i=1,m2
	  au(n1-i)=conjg(au(i))
 40	continue
	au(m2+1)=zero
	call fork(m,au,1.0)
	w=0.0
	dw=wi*tw/float(m-1)
	do 50 i=2,m
	  w=w+dw
	  au(i)=au(i)*exp(w)
 50	continue
c
c  loop for calculating time domain scattering energy decay
	dt=tw/float(m)
	c0=0.25/(vs*pi*r*r)
	c1=0.25*gs/(vs*r*pi*c0)
	c2=gs*gs/(16.0*r*pi*c0)
	c3=0.5/(vs*r*pi*pi*tw*c0)
	gamma=2.0/(9.0*sqrt(3.0))
	r2=r*r
	do 60 i=1,n
	  vt=vs*t(i)
	  if(vt.lt.r)then
	    tau=exp(-g*vt)
	    p(i)=gamma*c1*alog((vt+r)/(r-vt))/t(i)*tau
	    goto 55
	  endif
	  tau=r/vt
	  call ener2(tau,0.01,coef)
	  n1=ifix(t(i)/dt)+1
	  pd=real(au(n1))+(real(au(n1+1))-real(au(n1)))
     +	     *(t(i)-dt*float(n1-1))/dt*exp(gi*vt)
	  tau=exp(-gs*vt)
	  p(i)=c1*alog((vt+r)/(vt-r))/t(i)*tau+c2*coef*tau+c3*pd
 55	  p(i)=sqrt(p(i))
 60	continue
	return
	end
c
c.................................................................
c
	subroutine ener2(tau,dd,f)
c
c  calculate coefficient of the second-order scattered energy term
	pi=3.1415926
	n=ifix(tau/dd)+1
	dt=tau/float(n)
	f=0.0
	t=0.0
	do 10 i=1,n
	  t=t+dt
	  aux=alog((1.0+t)/(1.0-t))
	  aux=aux*aux
	  f=f+aux
 10	continue
	f=3.0*(0.5*aux-f)*dt
c
	f=f+tau*pi*pi
	return
	end
c
c.....................................................
c
	subroutine evwi(gi,gs,vs,r,tw,wi)
c
c  estimating the value of imaginery frequency wi
	pi=3.14159265
	dd=0.01
	vt=r
	do 20 i=1,2
	  vt=vt+(vs*tw-r)*(i-1)
	  ud=sqrt(0.75*gs*vt)
	  m=ifix(ud/dd)+1
	  du=ud/float(m)
	  coef=0.5
	  u=0.0
	  do 10 j=1,m
	    u=u+du
	    aux=exp(-u*u)
	    coef=coef+aux
 10	  continue
	  coef=(coef-0.5*aux)*du-ud*exp(-ud*ud)
	  u=0.75*gs/vt
	  aux=(1.-exp(-gs*vt)*(1.+gs*vt))/coef
	  di=aux*exp(1.5*alog(u/pi)-u*r*r-gi*vt)
	  if(i.eq.1)then
	    dimax=di
	  end if
 20	continue
	wi=alog(500*di/dimax)/tw
c
	return
	end
c
c.....................................................................
c
c  Purpose : for calculating the synthetic seismograms in a layered
c            1-d medium.
c  Input :
c            nt - number of time points
c            twin, fmax - time window and maximum frequency
c            nlay - number of layer
c            dh - layer thickness
c            v - wave velocity
c            den - medium density
c
c                                           -- Yuehua Zeng --
c                                      Updated in Feb., 1993, at UNR
c
c.....................................................................
c
	subroutine scat1d(nt,twin,fmax,au)
	parameter (m1=500)
	dimension v(m1),den(m1),dh(m1)
	complex HTu,HRu,exd(m1),au(nt/2)
        common/random/iseed,rnt(100000)
c
        open(20,file='scat1d.in',form='formatted',status='old')
	read(20,*)nlay,damp,dha,ddh,va,dv,dena,dden,Q,seed
        iseed=mod(int(seed),99991)
	close(20)
c
c  set the parameters
	pi2=2.0*3.141592653589793
	df=1./twin
	nw=ifix(fmax/df)+1
	if(damp.le.0.0)damp=1.0
	wi=-damp*pi2*0.5/twin
	v15=dv*2.0
	d15=dden*2.0
c
c  generate layer thickness, velocity and density
	v(1)=va
	den(1)=dena
	dh(1)=dha
	do 10 i=2,nlay-1
  5	  v1=2.*rnt(iseed)-1.
          iseed=iseed+1
          iseed=mod(iseed,99991)
	  v(i)=v(i-1)+dv*v1
	  den(i)=den(i-1)+dden*v1
	  if(abs(v(i)-va).gt.v15.or.abs(den(i)-dena).gt.d15)goto 5
	  dh(i)=dha+ddh*(2.*rnt(iseed)-1.)
          iseed=iseed+1
          iseed=mod(iseed,99991)
10	continue
	v(nlay)=va
	den(nlay)=dena
	dh(nlay)=dha
c
	df=df*pi2
	f=-df
	do 40 iw=1,nw
	  f=f+df
	  do 20 i=1,nlay-1
	    hru=dh(i)/v(i)*cmplx(1.0,-0.5/Q)
	    exd(i)=exp(cmplx(wi,-f)*hru)
 20	  continue
	  exd(nlay)=1.0
c
c  generalized reflection and transmission coefficients for
c  layers above the source.
	  au(iw)=1.0
          HRu=exd(1)
	  do 30 j=1,nlay-1
	    vd1=v(j)*den(j)
	    vd2=v(j+1)*den(j+1)
	    r=(vd1-vd2)/(vd1+vd2)
 	    HTu=(1.0-r)/(1.0-r*exd(j)*exd(j)*HRu)
	    HRu=-r+(1.0+r)*exd(j)*exd(j)*HRu*HTu
	    au(iw)=au(iw)*HTu
 30	  continue
 40	continue
	do 50 iw=nw+1,nt/2
	  au(iw)=0.0
 50	continue
c
	return
	end

c
c
ccccccccccccccccccccccccccc
      subroutine scoda2(nst,ngf,u1)
c
ccccccccccccccccccccccccccc
c
c   Program: scoda2.f
c   Purpose: apply a realistic coda to any general seismogram
c
c   Approach:
c     Read three tables of information.
c     Each table has one row for each station, in the sequence
c        listed in stat.dat. The station is indexed by 'nst'.
c     Each column is for a different subfault where the Green's function
c        was calculated.
c     The coda term is found in the frequency domain at a grid of
c        distances.
c
c
c Explanation of variables
c   Input
c     nst: index number of the station where this GF applies
c     ngf: index number of the subfault where this GF applies
c     jcm: index for coda model, usually scatterer density to apply
c     u1: Green's function seismogram, 3 columns, for the 3 directions
c          This is in the frequency domain, complex numbers
c   Obtained from a file.
c     Files tabcoda1, tabcoda2, tabcoda3 (first row)
c       nsmstn: the number of strong motion stations to be modeled
c       nsubf: the number of subfaults (& Green's functions) on the fault
c       nrcoda: the number of distances where coda model is found
c       ncodamod: the number of alternative coda models
c     File tabcoda1
c       res: distance from the subfault center to the station (not used yet)
c     File tabcoda2
c
c     File codafilexxx.dat, where xxx is a number between 001 and 999
c       af: Fourier amplitude spectrum of the coda, complex
c
c
ccccccccccccccccccccccc
c
      parameter(mm1=8200)
      dimension xistn(300),res(300,15),xkres(300,15),xjcm(300,15)
      complex af(mm1),u1(mm1,3)
      character*20 fn
c
c
c      print*,nst
c
c      read in the tables
      open(80,file='tabcoda1.dat',form='formatted',status='old')
      read(80,'(4e16.7)') xnsmstn,xnsubf,xnrcoda,xncodamod
      nsmstn=ifix(xnsmstn)
      nsubf=ifix(xnsubf)
      nrcoda=ifix(xnrcoda)
      ncodamod=ifix(xncodamod)
      do 10 i=1,nsmstn
	 read(80,'(16e16.7)') xistn(i),(res(i,j),j=1,nsubf)
 10   continue
      i=4
c      print*, xistn(i),(res(i,j),j=1,nsubf)
      close(80)

      open(80,file='tabcoda2.dat',form='formatted',status='old')
      read(80,'(2e16.7)') xnsmstn,xnsubf
      do 20 i=1,nsmstn
      read(80,'(16e16.7)') xistn(i),(xkres(i,j),j=1,nsubf)
 20   continue
      i=4
c      print*, xistn(i),(xkres(i,j),j=1,nsubf)
      close(80)

      open(80,file='tabcoda3.dat',form='formatted',status='old')
      read(80,'(2e16.7)') xnsmstn,xnsubf
      do 30 i=1,nsmstn
      read(80,'(16e16.7)') xistn(i),(xjcm(i,j),j=1,nsubf)
 30   continue
      i=4
c      print*, xistn(i),(xjcm(i,j),j=1,nsubf)
      close(80)

      xres=res(nst,ngf)
      kres=ifix(xkres(nst,ngf))
      jcm=ifix(xjcm(nst,ngf))
c      print *, xres,kres, jcm

      kseq=jcm+(kres-1)*ncodamod
c      print *, nst,ngf,xres,kres,jcm,kseq

      fn='codafile000a.dat'
      if (kseq.le.9) then
        write(fn(11:11),'(i1)')kseq
      else if (kseq.le.99) then
        write(fn(10:11),'(i2)')kseq
      else
        write(fn(9:11),'(i3)')kseq
      endif

c      print *, fn

      open(80,file=fn,form='formatted',status='old')
      read(80,'(3e16.7)') rcoda,sdens,xnf
      read(80,'(3e16.7)') sc2,wx1,wx2
      read(80,'(3e16.7)') qd,vs,zmoho
      nf=ifix(xnf)
      do 40 k=1,nf
      read(80,'(3e16.7)') f,a,b
      af(k)=cmplx(a,b)
 40   continue
      close(80)

c      print*,f,a,b
c      print*,af(5)

      do 60 j=1,3
      do 60 i=1,nf
      u1(i,j)=u1(i,j)*af(i)
 60   continue

      return
      end
