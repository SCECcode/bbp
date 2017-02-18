!x
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Program NOAH, version 1.0 (FORTRAN 90)                    c
!                                                                c
!      Program that computes the 1D nonlinear wave propagation   c
!      using a second order staggered grid scheme.  This code    c
!      does not allocates memory dynamically.                    c
!                                                                c
!      This code uses the Towhata-Ishihara-Iai-Iwan-GMR Model    c
!                                                                c
!      Author: Fabian Bonilla, IRSN, 2002                        c
!                                                                c
!      If you have any problems please contact me at:            c
!      Fabian Bonilla                                            c
!      fabian@crustal.ucsb.edu                                   c
!      fabian.bonilla@irsn.fr                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program noah
      implicit none

      integer nx, nt
      integer ix, it, iskip, itt, ir, itt1
      real    dx, dt, Kf, rvelb, dt00

      integer imec, k, itj, ixj, iskix, IMmax, Nspr
      real    S, Sb, S0, ws

      real    rr1, rr2, rr3, sd
      real    minuto0, minuto1, suma

      integer year, mes, dia, hora, minute, seg, mseg

      include 'params.dec'

!      If the user wants to increase the number of springs IMmax
!      (multishear model), or the number of springs in parallel
!      Nspr (Iwan model), then change the two following lines

      parameter (Nspr=50, IMmax=60)

!      Material properties

      real mu(nx), Ka(nx), rho(nx), so(nx), eo(nx)
      real m1(nx), m2(nx), m3(nx), p1(nx), p2(nx), S1(nx)
      real w1(nx), c1(nx), d0(nx), sigma0(nx), tau0(nx)
      real poro(nx), w0(nx), ee0(nx), B(nx)
      real sigx(nx), sigy(nx), sigxy(nx)

      logical lpow(nx), mtype(nx)

!      G/Gmax variables: maximum number of soil types is 100
!      maximum number of data points per soil model is 100

      real    straing(100,100), GoverG(100,100)
      integer soil(nx), Npg(100)

!      Iwan model variables

      real alf(Nspr,IMmax,nx), pbar(2,IMmax,nx)
      real HH(nx,Nspr), RR(nx,Nspr)

!      Q memory variables

      real q(nx,4), chi(nx,4), relaxt(4)

!      Input time history

      real rs(nt)

!      Character variables

      character filout*80
      character ibound*80

!      Constants

      real pi

!      Wave propagation variables

      real u(nx), udot(nx,0:1)

!      Local Variables

      real    pval(7,IMmax,nx), pvalp(11,nx)
      real    dtheta(nx), epsx(nx), epsy(nx), epsxy(nx), epsxy0(nx)
      integer IM(nx)
      logical lkh(nx)

!      Output time histories

      real velth(nt), strainth(nt), stressth(nt), poreth(nt)

!      Timing variables

      character date*8
      character time*10
      character zone*5
      integer   values(8)

!      Allocating memory

      open(14,file='intermed.par',status='old')
         read(14,*) 
         read(14,*) 
         read(14,*) dt
         read(14,*) dx
         read(14,*) iskix
         read(14,*) iskip
      close(14)

      if (iskix .gt. nx) then
         write(*,*) 'Error: iskix > nx (Check intermed.par)'
         stop
      end if

!      Beginning of the process, checking time

      call date_and_time(date,time,zone,values)
      year    = values(1)
      mes     = values(2)
      dia     = values(3)
      hora    = values(5)
      minute  = values(6)
      seg     = values(7)
      mseg    = values(8)
      minuto0 = hora * 60.0 + minute + (seg + mseg/1000.0)/60.0

      call Init(nx,nt,mu,Ka,rho,so,eo,m1,m2,m3,p1,p2,S1,      
     &          w1,c1,d0,sigma0,poro,w0,ee0,sigx,sigy,        
     &          sigxy,tau0,relaxt,q,dx,dt,ir,rvelb,           
     &          Kf,pi,rs,B,ibound,lpow,mtype,filout,
     &          Npg,soil,straing,GoverG,dt00)

      write(*,*)
      write(*,*) 'Space and time node location results (multiplexed data
     &):'
      write(*,*) '   Space shift = ', iskix, ' dx_out = ', nx/iskix, ' (
     &nodes)'
      write(*,*) '   Time shift  = ', iskip, ' dt_out = ', dt*iskip

      write(*,*)
      write(*,*) 'Processing Initial Conditions ...'

      do ix = 1,nx

         u(ix)      = 0.0
         udot(ix,0) = 0.0
         lkh(ix)    = .true.
         S          = 1.0
         IM(ix)     = 12
         dtheta(ix) = pi/IM(ix)

!         call curve(Nspr,eo(ix)*atan(1.0),so(ix)/2.0,HH(ix,1),RR(ix,1))
         call curve_data(Nspr,mu(ix),Npg(soil(ix)),straing(soil(ix),:),
     &                   GoverG(soil(ix),:),HH(ix,:),RR(ix,:))

         call InitPp(tau0(ix)/sigma0(ix),m1(ix),m2(ix),m3(ix),  
     &               p1(ix),p2(ix),S1(ix),w1(ix),w0(ix),        
     &               poro(ix),Sb,S0,ws)

         call strain(IM(ix),ix,so(ix),eo(ix),ee0(ix),Sb,S0,S,   
     &               dtheta(ix),m1(ix),m2(ix),sigma0(ix),       
     &               sigx(ix),sigy(ix),sigxy(ix),               
     &               epsx(ix),epsy(ix),epsxy0(ix))

         do imec = 1,IM(ix)
            pval(1,imec,ix) = 0.0           !xtp     previous strain value
            pval(2,imec,ix) = 0.0           !eb      reversal on backbone
            pval(3,imec,ix) = 0.0           !er      strain reversal
            pval(4,imec,ix) = 0.0           !sr      stress reversal
            pval(5,imec,ix) = 1.0           !bvalue  strain scale factor
            pval(6,imec,ix) = 1.0           !rc      hysteresis shape factor
            pval(7,imec,ix) = int(1)        !j       reversal number 
         end do

         pvalp(1,ix)  = epsx(ix)            !exx
         pvalp(2,ix)  = epsy(ix)            !eyy
         pvalp(3,ix)  = epsxy0(ix)          !exy
         pvalp(4,ix)  = S                   !S
         pvalp(5,ix)  = S0                  !S0
         pvalp(6,ix)  = Sb                  !Sb
         pvalp(7,ix)  = ws                  !ws
         pvalp(8,ix)  = tau0(ix)            !tau
         pvalp(9,ix)  = sigx(ix)            !sxx
         pvalp(10,ix) = sigy(ix)            !syy
         pvalp(11,ix) = sigxy(ix)           !sxy

      end do

      write(*,*)
      write(*,*) 'Processing Wave Propagation ...'

!      Time loop

      itj   = 0
      itt1  = 1
      it    = 1
      do itt=1,nt


         do ix = 1,nx
            u(ix) = u(ix) + dt * udot(ix,it-1)
         end do

         do ix = 1,nx-1

            epsxy(ix) = ( u(ix+1) - u(ix) ) / dx

            if (.not. mtype(ix)) then

             suma = 0.0
             do k=1,4
                chi(ix,k) = chi(ix,k)*exp(-dt/relaxt(k)) +   
     &                      q(ix,k)*epsxy(ix)*               
     &                      (1.0-exp(-dt/relaxt(k)))
                suma      = suma + chi(ix,k)
             end do

             epsxy(ix) = epsxy(ix) - suma
             sigxy(ix) = mu(ix) * epsxy(ix)

            else

               epsxy(ix) = epsxy(ix) + epsxy0(ix)

               call emulti(Nspr,ix,pval(:,:,ix),pvalp(:,ix),           
     &                     pbar(:,:,ix),alf(:,:,ix),HH(ix,:),RR(ix,:), 
     &                     d0(ix),lkh(ix),IM(ix),dtheta(ix),  
     &                     lpow(ix),B(ix),dt,epsx(ix),epsy(ix),     
     &                     epsxy(ix),so(ix),eo(ix),p1(ix),p2(ix),w1(ix),
     &                     S1(ix),c1(ix),m1(ix),m2(ix),m3(ix),      
     &                     ee0(ix),poro(ix),w0(ix),Kf,sigma0(ix),   
     &                     sigx(ix),sigy(ix),sigxy(ix))

            end if

         end do

         sigxy(1) = 0.0

         if (ibound(1:7) .eq. 'elastic') then
            rr1      = dt / (rho(nx-1) * dx)
            rr2      = 1.0 / (1.0 + rr1 * rho(nx) * rvelb/2.0)
            rr3      = 1.0 - rr1 * rho(nx) * rvelb/2.0
            udot(nx,it) = rr2 * (rr3*udot(nx,it-1) + rr1 *     
     &                    ( 2.0 * rho(nx) * rvelb * rs(itt) -  
     &                    sigxy(nx-1) ))
            sigxy(nx)   = rho(nx) * rvelb * ( 2.0 * rs(itt) -  
     &                    0.5*(udot(nx,it)+udot(nx,it-1)) )
         else
            udot(nx,it) = rs(itt)
            sigxy(nx)   = rho(nx) * rvelb * udot(nx,it)
         end if

         do ix = 2,nx-1
            sd          =  ( sigxy(ix) - sigxy(ix-1) ) / dx
            udot(ix,it) =  udot(ix,it-1) + sd * dt / rho(ix)
         end do

         if (mod(itt,iskip) .eq. 0) then
c            write(2,*) udot(ir,it), epsxy(ir), sigxy(ir),      
c     &                 sigma0(ir) * pvalp(4,ir)/1000.0
c            call flush(2)
            velth(itt1)    = udot(ir,it)*100.0
            strainth(itt1) = epsxy(ir)*100.0
            stressth(itt1) = sigxy(ir)/1000.0
            poreth(itt1)   = sigma0(ir) * (1.0-pvalp(4,ir))/1000.0
            itt1 = itt1 + 1

            ixj = 0
            do ix=2,nx,iskix
c               write(3) ix, udot(ix,it), epsxy(ix), sigxy(ix), 
c     &                  sigma0(ix) * pvalp(4,ix)
             itj = itj+1
             ixj = ixj+1
c             call flush(3)
            end do
         end if

         do ix=1,nx
            udot(ix,it-1) = udot(ix,it)
         end do

         it = 1

      end do

!      Writing the time and space windows for the multiplexed data

!      write(1,*) 'number of time (nf) and space (nk) windows'
!      write(1,2) itj/ixj, ixj
!2     format(2(i5,2x))
2     format('        parameter(nf = ',i5,', nk = ',i5,')')
!      close(1)

!      Output of time histories in SAC format at receiver level

      call salida(filout,velth,strainth,stressth,poreth,itt1-1,
     &            dt*iskip,rs,dt00)

!      End of the process

      call date_and_time(date,time,zone,values)
      year    = values(1)
      mes     = values(2)
      dia     = values(3)
      hora    = values(5)
      minute  = values(6)
      seg     = values(7)
      mseg    = values(8)
      minuto1 = hora * 60.0 + minute + (seg + mseg/1000.0)/60.0

      write(*,*)
      write(*,*) 'Elapsed Time (min) = ', minuto1-minuto0
      write(*,*)

      stop ' **** CIAO **** '
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine salida(filin,v,e,s,pp,n,dt,rs,dt00)

!      This subroutine directs the output of NOAH to SAC files
!      located at the receiver position

      implicit none

      integer n, nerr, nchar, i, np, npp, ni,len_trim
      real    dt, fhigh, flowp, dt00
      real    v(n), e(n), s(n), pp(n), rs(n)

      character filin*80, filout*80
!
      nchar = len_trim(filin)
!
!      Lowpassing the velocity to prevent aliasing
      np    = 4
      npp   = 2
      flowp = 0.45/dt00
      fhigh = 0.01

      call xapiir(v,n,'BU',0.0,0.0,np,'LP', 
     &            fhigh,flowp,dt,npp,n)
!
      call interp(dt,v,n,dt00,n,rs,ni)


!      Removing the mean and tapering the output signal.
!      If applying this, you will not see the permanent
!      displacement.  It is up to the user

	

      call rmean(rs,ni)
      call taper(rs,ni)
      do i=1,ni
         v(i) = rs(i)
      end do

!      Output the results in SAC format

!      filout = filin(1:nchar)//'.sac'
!      call wsac1(filout,v,ni,0.0,dt00,nerr)
!      filout = filin(1:nchar) // '.strain.sac'
!      call wsac1(filout,e,n,0.0,dt,nerr)
!      filout = filin(1:nchar) // '.stress.sac'
!      call wsac1(filout,s,n,0.0,dt,nerr)
!      filout = filin(1:nchar) // '.water.sac'
!      call wsac1(filout,pp,n,0.0,dt,nerr)
	filout = filin(1:nchar)
	 open(19, file=filout,status='unknown')
	  write(19,*) ni,dt
	  write(19,925) (v(i), i=1,ni)
	  close(19)
925	  format(5e14.5)	

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp(dt0,y0,n0,dti,nt,yi,ni)

!      Interpolation subroutine

      implicit none

      integer Ndmax, i, j, n0, ni, nt, j0, j1
      real    dt0, dti
      real    difu1, difu2, slop, abdu, ddd
      real    rdt, tt, aa, bb, ab

      parameter(Ndmax=500000)

      real y0(*), yi(*), ydif(Ndmax)

      if (nt .gt. Ndmax) then
         write(*,*) '*** Error in interpolation subroutine ***'
         stop
      end if

      do 10 j=1,n0
        if (j .eq. 1) then
          difu1 = y0(1)
          difu2 = y0(2)-y0(1)
        elseif (j .eq. n0) then
          difu1 = y0(n0)-y0(n0-1)
          difu2 = -y0(n0)
        else
          difu1 = y0(j)-y0(j-1)
          difu2 = y0(j+1)-y0(j)
        endif
        slop = 0.0
        abdu = abs(difu1)
        if ((difu1*difu2) .gt. 0.0) then
          if (abs(difu2) .gt. abdu) then
            ddd = difu1/difu2
            if (ddd .gt. 0.0) slop = 2.0*difu1/(1.0+ddd)
          elseif (abdu .gt. 0.0) then
            ddd = difu2/difu1
            if (ddd .gt. 0.0) slop = 2.0*difu2/(1.0+ddd)
          endif
        endif
        ydif(j) = slop
10    continue

c      write(*,*) n0, dt0, dti	
      ni  = int((n0-1)*dt0/dti)
      if (ni .gt. nt) ni = nt
!      if (ni .lt. 8192) ni = 8192
      rdt = dti/dt0
      do 20 i=1,ni
        tt    = float(i-1)*rdt
        j0    = int(tt)+1
        if (j0 .ge. n0) j0=n0-1
        j1    = j0+1
        aa    = float(j0)-tt
        bb    = 1.0-aa
        ab    = aa*bb
        yi(i) = aa*aa*(3.-2.*aa)*y0(j0)+     
     &          bb*bb*(3.-2.*bb)*y0(j1)+     
     &          ab*(aa*ydif(j0)-bb*ydif(j1))
20    continue

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Init(nx,nt,mu,Ka,rho,so,eo,m1,m2,m3,p1,p2,S1,  
     &                w1,c1,d0,sigma0,poro,w0,ee0,sigx,sigy,   
     &                sigxy,tau0,relaxt,q,dx,dt,ir,rvelb,     
     &                Kf,pi,rs,B,ibound,lpow,mtype,filout,
     &                Npg,soil,straing,GoverG,dt00)

!      This subroutine initializes the data arrays, reads the
!      source acceleration, interpolates, and applies a high pass
!      filter to the input source.

      implicit none

!      Declaring all variables

      integer nx, nt, Nsoil
      real    dx, dt
      real    freqmax, rw, dt0, flowp, fhigh, rir
      integer nlayers, iwat, imin, rsoil, kg, kd, Npd
      real    rtop, rvp, rvs, rrho, rnu, rsigma0
      real    rporo, rko, rphi, rpha, rp1, rp2
      real    rw1, rS1, rc1, rd0
      integer i, n, k
      real    poisson, cM, rvelb, rhob, deltm
      integer ndata, nerr, np, npp, ni, ir, narg, iargc
      real    beg, dt00

!      Variables that will need memory to be allocated

!      Material properties

      real  mu(nx), Ka(nx), top(nx), rho(nx), so(nx), eo(nx)
      real  m1(nx), m2(nx), m3(nx), p1(nx), p2(nx), S1(nx)
      real  w1(nx), c1(nx), d0(nx), sigma0(nx), tau0(nx)
      real  poro(nx), w0(nx), ee0(nx), B(nx), Ko(nx)
      real  sigx(nx), sigy(nx), sigxy(nx)

      logical lpow(nx), mtype(nx)

!      G/Gmax variables: maximum number of soil types is 100
!      maximum number of data points per soil model is 100

      real    straing(100,100), GoverG(100,100)
      integer soil(nx), Npg(100), soil_type(100)

!      Q memory variables

      real q(nx,4), relaxt(4), qwt1(4), qwt2(4)

!      Input time history

      real ai(nt), vi(nt), rs(nt)

!      Character variables

      character fsoil*80
      character fsource*80
      character filout*80
      character ibound*80
      character matcond*80
      character baratin*1
      character ftype*2

!      Constants

      real pi, grav, Kf, rhow

      pi   = 4.0 * atan(1.0)
      grav = 9.8
      Kf   = 2.0E9
      rhow = 1000.0
      rd0  = 0.3

      write(*,*)
      write(*,*) '   Noah: Nonlinear Site Response Analysis     '
      write(*,*) '   --------------------------------------     '
      write(*,*)
      write(*,*) '  Depth (m)    sigma0 (kPa)      eo ( % )     mu (MPa)
     &'

!      First reading of 'input.iwan' to compute sigx, sigy, sigma0

      open (15,file='input.iwan',status='old')
      read (15,*) freqmax
      read (15,'(a)') baratin
      read (15,'(a)') baratin
      read (15,'(a)') baratin
      read (15,'(a)') baratin
      read (15,*) rw
      read (15,*) dt0
      read (15,*) flowp
      read (15,*) fhigh
      read (15,*) ftype
      read (15,*) rir
      read (15,*) nlayers
      do n=1,nlayers
         read (15,'(a)') baratin
         read (15,*)  rtop
         read (15,*)  rvp
         read (15,*)  rvs
         read (15,*)  rrho
         read (15,*)  rnu
         read (15,'(a)') baratin
         read (15,*)  rsigma0
         read (15,*)  rporo
         read (15,*)  rko
         read (15,*)  rphi
         read (15,*)  rpha
         read (15,*)  rp1
         read (15,*)  rp2
         read (15,*)  rw1
         read (15,*)  rS1
         read (15,*)  rc1
         read (15,*)  rsoil
         iwat = int(rw/dx) + 1
         imin = int(rtop/dx) + 1
         do i = imin,nx
            rho(i) = rrho
            top(i) = imin
            Ko(i)  = rko
         end do
      end do

      sigy(1) = 0.0
      do i = 2,nx
         if (top(i) .ge. iwat) rho(i) = rho(i)-rhow
         sigy(i)   = sigy(i-1) + rho(i)*grav*dx
         sigx(i)   = sigy(i) * Ko(i)
         sigma0(i) = (sigx(i) + sigy(i))/2.0
         tau0(i)   = abs(sigx(i) - sigy(i))/2.0
      end do

      sigy(1)   = sigy(2)
      sigx(1)   = sigx(2)
      sigma0(1) = sigma0(2)
      tau0(1)   = tau0(2)

      rewind(15)

!      Second reading of 'input.iwan' to initialize the other parameters

      open (15,file='input.iwan',status='old')
      read (15,*) freqmax
      read (15,'(a)') baratin
      read (15,'(a)') baratin
      read (15,'(a)') baratin
      read (15,'(a)') baratin
      read (15,*) rw
      read (15,*) dt0
      read (15,*) flowp
      read (15,*) fhigh
      read (15,*) ftype
      read (15,*) rir
      read (15,*) nlayers
      do n=1,nlayers
         read (15,'(a)') baratin
         read (15,*)  rtop
         read (15,*)  rvp
         read (15,*)  rvs
         read (15,*)  rrho
         read (15,*)  rnu
         read (15,'(a)') baratin
         read (15,*)  rsigma0
         read (15,*)  rporo
         read (15,*)  rko
         read (15,*)  rphi
         read (15,*)  rpha
         read (15,*)  rp1
         read (15,*)  rp2
         read (15,*)  rw1
         read (15,*)  rS1
         read (15,*)  rc1
         read (15,*)  rsoil
         imin = int(rtop/dx) + 1
         do i = imin,nx
            soil(i)    =  rsoil
            d0(i)      =  rd0
            rho(i)     =  rrho
            poro(i)    =  rporo
            p1(i)      =  rp1
            p2(i)      =  rp2
            w1(i)      =  rw1
            S1(i)      =  rS1
            c1(i)      =  rc1
            so(i)      =  sigma0(i) * sin(rphi*pi/180.0)

!      Computing shear and bulk modulus

            mu(i)      =  rho(i)*rvs**2
            if (rvp .eq. 0.0) then
                 poisson =  0.33
                 cM      =  2.0*mu(i)*(1.0-poisson)/(1.0-2.0*poisson)
            else
                 cM      =  rho(i) * rvp**2
            end if
            Ka(i)      =  cM - 4.0*mu(i)/3.0

            if (rsigma0 .ne. 0.0) then
               mu(i)   = mu(i) * sqrt(sigma0(i)/rsigma0)
               Ka(i)   = Ka(i) * sqrt(sigma0(i)/rsigma0)
               B(i)    = (0.5*Ka(i)/sqrt(rsigma0))**2
               ee0(i)  = sqrt(sigma0(i)/B(i))
               lpow(i) = .true.
            else
               B(i)    = Ka(i)
               ee0(i)  = sigma0(i)/B(i)
               lpow(i) = .false.
            end if

!      Computing m1, m2, m3, eo, w0, B, ee0, Go, To

            m1(i)       =  sin(rphi*pi/180.0)
            m2(i)       =  sin(rpha*pi/180.0)
            m3(i)       =  0.67 * m2(i)
            eo(i)       =  so(i) / mu(i)
            w0(i)       =  so(i) * eo(i)/2.0

!      Linear or nonlinear condition (rnu = 0.0 linear, nonlinear
!       otherwise)

            if (rnu .eq. 0.0) then
               mtype(i) = .false.
               matcond  = 'linear'
            else
               if (poro(i) .le. 0.0) then
                  mtype(i) = .true.
                  matcond  = 'dry nonlinear'
               else
                  mtype(i) = .true.
                  matcond  = 'sat. nonlinear'
               end if
            end if

!      Printing basic variables at the top of each layer

            if (i.eq.imin) then
               write(*,2) dx*(i-1), sigma0(i)/1000.0, eo(i)*100.0, 
     &                    mu(i)/1.0E6, matcond
            end if

         end do
      end do
      read (15,*)
      read (15,*) rvelb
      read (15,*) rhob
      read (15,'(1a)') ibound
      close (15)

      rho(nx) = rhob

      write(*,*)
      if (ibound(1:7) .eq. 'elastic') then
         write(*,*) 'ELASTIC boundary conditions'
      else
         write(*,*) 'RIGID boundary conditions'
      end if

!      Computing variables for Q

      open(16,file='Q.par',status='old')
      do i=1,4
         read(16,*) relaxt(i), qwt1(i), qwt2(i)
      end do
      do n=1,nlayers
         read(16,*) rtop, deltm
         imin = int(rtop/dx) + 1
            do i=imin,nx
             do k=1,4
                  q(i,k) = deltm*(deltm*qwt1(k) + qwt2(k))
             end do
            end do
      end do
      close(16)

      narg = iargc()
      if (narg .eq. 3) then
        call getarg(1,fsoil)
        call getarg(2,fsource)
        call getarg(3,filout)
      elseif (narg .eq. 2) then
        call getarg(1,fsoil)
        call getarg(2,fsource)
        filout = 'pp'
        write(*,*)
        write(*,*) 'OUTPUT will be called pp'
      elseif (narg .eq. 0 .or. narg .eq. 1) then
        write(*,*)
        write(*,*) 'ERROR in input files:'
        write(*,*) 'NOAH USAGE: noah soil_parms input output'
        stop
      endif

!      Reading the soil G/Gmax data. Note that I am using TREMORKA
!      format

      open(16,file=fsoil,status='old')

      read(16,*) Nsoil

      do k=1,Nsoil
         read(16,*) soil_type(k)
         read(16,*) Npg(k)
         do kg=1,Npg(k)
            read(16,*) straing(k,kg), GoverG(k,kg)
            straing(k,kg) = straing(k,kg)/100.0
         end do
         read(16,*) Npd
         do kd=1,Npd
            read(16,*)
         end do
      end do

      close(16)

c      call file_delete('filter.par')
c      open(1,file='filter.par',status='new')
c      call file_delete('output.dat')
c      open(2,file='output.dat',status='new')
c      call file_delete('Nwavefield')
c      open(3,file='Nwavefield',status='new',form='unformatted')

      write(*,*)
      write(*,*) 'Reading the Velocity (m/s)'

      call rsac1(fsource,vi,ndata,beg,dt00,nt,nerr)

      if (flowp .gt. freqmax) then
         write(*,*)
         write(*,*) 'WARNING: flowp > freqmax'
      end if

!      Filtering the data with 4 poles and 2 passes
      write(*,*) ndata, np, ftype, fhigh, flowp, dt00

      !do i=1,ndata
        !write(*,*) vi(i)
      !end do

      np    = 4
      npp   = 2
      call xapiir(vi,ndata,'BU',0.0,0.0,np,ftype, 
     &            fhigh,flowp,dt00,npp,nt)

! added by J. Schmedes: I use input in cm, so go to m
      do i=1,ndata
         vi(i) = vi(i)*0.01 
         !write(*,*) vi(i)
      end do
!      Removing the mean and tapering the input signal

      call rmean(vi,ndata)
      call taper(vi,ndata)
      call rmean(vi,ndata)
      call taper(vi,ndata)

!      Calculating the acceleration and interpolating
      rs(1)=0.0
      do i=2,ndata
         rs(i) =(vi(i) - vi(i-1))/dt00
      end do
      call interp(dt00,rs,ndata,dt,nt,ai,ni)

!      Integrating the acceleration to velocity in m/s
      rs(1)=0.
      do i=2,ni
        rs(i) = rs(i-1) + ai(i-1)*dt
      end do
      call rmean(rs,ni)
      call taper(rs,ni)
      do i=ni,nt
        rs(i)=0.0
      enddo

      ir = int(rir/dx) + 1
      if (ir .gt. nx) stop 'receiver position exceeded maximum depth'

      write(*,*)
      write(*,*) 'receiver (m)   sigma0 (kPa)      eo ( % )     mu (MPa)
     &'
      write(*,2) dx*(ir-1), sigma0(ir)/1000.0, eo(ir)*100.0,  
     &           mu(ir)/1.0E6
2     format(4(f10.3,5x),1x,a14)
      if (ir .eq. 1) ir = 2

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine strain(IM,ix,To,eo,ee0,Sb,S0,S,  
     &                  dtheta,m1,m2,sigma0,      
     &                  sigx,sigy,sigxy,          
     &                  epsx,epsy,epsxy)

!      Computes the initial strains given the initial
!      stress values.  It uses the same rheology (multishear)

      implicit none

      integer Iter
      real    tol

      parameter(Iter=20, tol=1.0E-2)

      integer IM, ix, j, imec
      real    To, eo, ee0, Sb, S0, S, dtheta, m1, m2
      real    sigma0, sigx, sigy, sigxy
      real    epsx, epsy, epsxy
      real    T0, G0, e0, Qv, gv, sxsy, x_new, x_old
      real    difference, value_fun, deriv_fun
      real    thi, sigma, red
      real    xlambda, diff2

      real eps(3), sig(3)

      sig(1) = sigx
      sig(2) = sigy

!      The initial shear strain is zero

      eps(3) = 0.0

!      initial shear stress and shear modulus

        if (S0 .gt. Sb) then
           T0 = To*S
           G0 = (To/eo)*S
      else
           T0 = To*S + (m1-m2)*(Sb-S0)*(0.4/Sb)*sigma0
           e0 = eo/(S0/Sb)
           G0 = T0/e0
        end if

!      Multi-spring model

      Qv = T0/2.0
      gv = (T0/G0)*atan(1.0)

!      We actually solve for the stress difference

      sxsy = (sig(1) - sig(2))/2.0

      if (sxsy .eq. 0.0) then
         x_new = 0.0
         goto 21
      end if

      x_old   = -1.0E-6

      xlambda = 1.0
      diff2   = 1.0E10

      do 2 j=1,Iter
         call shearonly(IM,Qv,gv,dtheta,sxsy,x_old,   
     &                  value_fun,deriv_fun)
         x_new = x_old - value_fun/deriv_fun
         difference = ABS(x_new-x_old)/ABS(x_new)
         if (difference .lt. tol) then
            goto 21
         elseif (difference .lt. diff2) then
            diff2 = difference
            goto 633
         elseif (difference .gt. tol .and. difference .gt. diff2) then
            if (xlambda .gt. 0.0001) then
               xlambda = xlambda/10.0
            else
               xlambda = 0.0001
            end if
            x_new = x_old - xlambda*(value_fun/deriv_fun)
            goto 633
         elseif (j .ge. Iter) then
            write(*,*) 'NoCNR after ', j, ' Iter at ', ix, difference
         end if
633         x_old = x_new
2      continue

21      continue

      eps(1) = (x_new + ee0)/2.0
      eps(2) = ee0 - eps(1)

      sig(3) = 0.0
      sigma  = 0.0
      do imec=1,IM
         thi = dtheta*(imec-1)
         red = (eps(1)-eps(2))*cos(thi)
         call Fbb(Qv,gv,1.0,1.0,1.0,0.0,0.0,red,sigma)
         sig(3) = sig(3) + sigma*dtheta*sin(thi)
      end do

      epsx  = eps(1)
      epsy  = eps(2)
      epsxy = eps(3)

      sigx  = sig(1)
      sigy  = sig(2)
      sigxy = sig(3)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shearonly(IM,Qv,gv,dtheta,sxsy,x_old, 
     &                     fun,deriv)

!      This subroutine computes the function and its
!      derivative for the newton-raphson iteration.
!      In this case the function is the multishear
!      mechanism.

      implicit none

      integer imec, IM
      real    Qv, gv, dtheta, sxsy, x_old, fun, deriv
      real    sigxy, dsigxy, sigma, dsigma, eps, thi

      sigxy  = 0.0
      dsigxy = 0.0
      sigma  = 0.0
      dsigma = 0.0
      do imec=1,IM
         thi = dtheta*(imec-1)
         eps = x_old*cos(thi)
         call Fbb(Qv,gv,1.0,1.0,1.0,0.0,0.0,eps,sigma)
         call Fbb1(Qv,gv,eps,dsigma)
         sigxy  = sigxy + sigma*dtheta*cos(thi)
         dsigxy = dsigxy + dsigma*dtheta*cos(thi)*cos(thi)
      end do

      fun   = sigxy - sxsy
      deriv = dsigxy

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Fbb1(Qv,gv,eps,sigma)

!      This subroutine computes the derivative of the
!      hyperbolic stress-strain model

      implicit none

      real Qv, gv, eps, sigma

      sigma = (Qv/gv) / (1.0 + abs(eps/gv))**2

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine InitPp(r,m1,m2,m3,p1,p2,S1,w1,w0,poro,Sb,S0,ws)

!      Subroutine that computes the initial values of
!      liquefaction front S0, and shear work ws.  It also
!      computes the correction value for the dilatancy
!      threshold Sb

      implicit none

      real m1, m2, m3, m4
      real r, p1, p2, S1, w1, w0, poro, Sb, S0, ws
      real a, b, c, w

      if (r .le. m3 .or. poro .le. 0.0) then
         ws = 0.0
         S0 = 1.0
         Sb = 0.4
      else
         m4 = 1.0 - (m2-m3)/m1
         a  = m4*m4*m1*m1 - m2*m2 - 2.0*m3*m3 + 2.0*m2*m3
         b  = 2.0*r*m3 - 2.0*m1*m1*m4
         c  = m1*m1 - r*r
         S0 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a)
         if (S0 .le. S1) S0 = S1
         if (S0 .ge. 0.4) then
            Sb = 0.4
         else
            Sb = S0
         end if
         if (S0 .ge. 0.4) then
            w = w1 * ((1.0-S0)/0.6)**(1.0/p1)
         else
            if (S0 .le. S1) then
               w = 1.0E10
            else
             w = w1 * ((S0-S1)/(0.4-S1))**(-1.0/p2)
            end if
         end if
         ws = w*w0
      end if

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine Fbb(T0,e0,av,bv,rc,nu,er,eps,sigma)

!      This subroutine computes the hyperbolic stress-strain
!      model.  The different conditions are for the value
!      of the hysteresis shape factor, rc.

      implicit none

      real nu, T0, e0, av, bv, rc, er, eps, sigma
      real rnum, rden

      if (rc .eq. 0.0) then
         sigma = T0 * (eps/e0) / (1.0 + abs(eps/e0)) + nu*er
      elseif (rc .eq. -999.0) then
         sigma = T0 * (eps/e0)*(av/bv) + nu*er
      else
         rnum  = T0 * (eps/e0)*(av/bv)
         rden  = 1.0 + abs(eps/(rc*e0*bv))
         sigma = rnum/rden + nu*er
      end if

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine taper(xt,n)

!       Taper using a Hanning window of 5% for each side

      implicit none

      integer n, j, taperwxt
      real    f0, f1, omega
      real    xt(n), wxt(n)

      f0=0.50
      f1=0.50
      taperwxt=int(0.05*n)
      omega=4.0*atan(1.0)/taperwxt

       do 2 j=1,taperwxt
          wxt(j)=f0-f1*cos(omega*(j-1))
   2   continue

       do 21 j=taperwxt+1,n-taperwxt
          wxt(j)=1.0
   21  continue

       do 22 j=n-taperwxt+1,n
          wxt(j)=1.0-wxt(j-n+taperwxt)
   22  continue

       do 6 j=1,n
          xt(j)=xt(j)*wxt(j)
   6   continue

       return
       end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rmean(yarray,nlen)

!      This subroutine removes the mean of a given time history

      implicit none

      integer nlen, i
      real    yarray(nlen), ave

      call moment(yarray,nlen,ave)
      do i=1,nlen
         yarray(i) = yarray(i) - ave
      end do
      
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine moment(data,n,ave)

!      This subroutine computes the mean of a given data array

      implicit none

      integer n
      real    ave, data(n)
      integer j
      real    s

      if (n.le.1) stop 'n must be at least 2 in moment'
      s = 0.0
      do j=1,n
         s=s+data(j)
      enddo
      ave=s/n

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine file_delete ( file_name )
!
!***********************************************************************
!
!! FILE_DELETE deletes a named file if it exists.
!
!
!  Discussion:
!
!    You might want to call this routine to get rid of any old copy
!    of a file, before trying to open a new copy with the OPEN argument:
!      status = 'new'.
!
!    It's not always safe to open a file with " STATUS = 'UNKNOWN' ".
!    For instance, on the SGI, the most recent version of the FORTRAN
!    compiler seems to go crazy when I open an unformatted direct
!    access file this way.  It creates an enormous file (of somewhat
!    random size).  The problem goes away if I delete any old copy
!    using this routine, and then open a fresh copy with
!    " STATUS = 'NEW' ".  It's a scary world.
!
!  Modified:
!
!    26 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*(*) FILE_NAME, the name of the file to be deleted.
!
      character ctemp*80
      character (len=*) file_name
      integer iunit
      integer lens
      logical lfile
!
!  Does the file exist?
!
      inquire (           
     &  file = file_name, 
     &  exist = lfile )

      if ( .not. lfile ) then
        return
      end if

!
!  Get a free unit number.
!
      call get_unit ( iunit )

      if ( iunit .eq. 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'FILE_DELETE: Warning!'
        write ( *, * ) '  A free FORTRAN unit could not be found.'
        return
      end if

      write ( *, * ) ' '
      ctemp = 'FILE_DELETE: deleting old version of ' // file_name
      lens = len_trim(ctemp)
      write ( *, '(a)' ) ctemp(1:lens)

      open (                
     &   unit = iunit,      
     &   file = file_name,  
     &   status = 'old',    
     &   err = 10 )

      close (               
     &  unit = iunit,       
     &  status = 'delete' )

      return

10    continue

      write ( *, * ) ' '
      write ( *, * ) 'FILE_DELETE: Warning!'
      write ( *, * ) '  Could not open the file.'

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_unit ( iunit )
!
!***********************************************************************
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    If IUNIT < 0, then an I/O error occurred while trying to inquire
!    on the status of unit abs(IUNIT).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
      integer i
      integer iunit
      logical lopen

      iunit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 ) then

          iunit = -i

          inquire (           
     &      unit = i,         
     &      opened = lopen,   
     &      err = 10 )

          if ( .not. lopen ) then
            iunit = i
            return
          end if

        end if

      end do
!
!  No free unit was found.
!
      iunit = 0

      return
!
!  An I/O error occurred during an INQUIRE.
!
10    continue

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine emulti(Nspr,ix,pval,pvalp,pbar,alf,H,RR,   
     &                  d0,lkh,IM,dtheta,lpow,B,dt,           
     &                  epsx,epsy,epsxy,To,eo,p1,p2,w1,S1,c1, 
     &                  m1,m2,m3,ee0,poro,w0,Kf,sigma0,       
     &                  sigmax,sigmay,sigmaxy)

!      Computation of the multishear mechanism for each strain
!      INPUT: variables and strains
!      OUPUT: stresses and memory variables

!      IMPORTANT NOTE: part of the output is also the updated vector
!      'pvalp'.  The value 'pvalp(4)' is the state variable of
!      the liquefaction front 'S'.  Thus, at every time step
!      the excess of pore pressure has to be computed as:
!
!      pore pressure excess = pvalp(4) * sigma0
!
!      This, of course, has to be done either inside of this
!      subroutine or as soon as it is called from outside.  We
!      did the latter in the finite difference code.

      implicit none

      integer ix, IM, imec, ipow, Nspr
      real    d0, dtheta, B, dt, epsx, epsy, epsxy
      real    To, eo, p1, p2, w1, S1, c1, nu
      real    m1, m2, m3, Kf, ee0, poro, W0
      real    sigma0, sigmax, sigmay, sigmaxy
      real    Qv, gv, shear, ep, S0, S, xt, ws
      real    exx, eyy, exy, sxx, syy, sxy
      real    thi, r, r2, r3, taup, ssl, Rc
      real    dexey, dexy, dwst, taunew, sold, snew
      real    dGop, dGo, dws, dwse, w, S2
      real    sigx, sigy, sigxy, Sb, T0, G0, e0
      logical lkh, lpow

      real    eps(3), sig(3), pval(7,IM), pvalp(11)
      real    H(Nspr), RR(Nspr), alf(Nspr,IM), pbar(2,IM)

      nu     = 0.0

      eps(1) = epsx
      eps(2) = epsy
      eps(3) = epsxy

!      Case when no pore pressure is generated
!      This version uses either the GMR or the 1D Iwan model.

      if (poro .le. 0.0) then

         Qv    = To /2.0
         gv    = eo * atan(1.0)
         lkh   = .false.

         sigx  = 0.0
         sigy  = 0.0
         sigxy = 0.0
         shear = 0.0
         ep    = 0.0
         S     = 1.0
         S0    = 1.0

         do imec=1,1 !IM
!            thi = dtheta*(imec-1)
!            xt  = (eps(1)-eps(2))*cos(thi) + eps(3)*sin(thi)
            xt  = eps(3)
!            if (d0 .le. 0.0) then
               call Iwan(Nspr,H,RR,alf(:,imec),pbar(:,imec),xt,shear)
!            else
!               call smulti(ix,d0,lkh,xt,dt,Qv,gv,         
!     &                     nu,pval(:,imec),S0,S,shear)
!            end if
!            sigx  = sigx  + shear*dtheta*cos(thi)
!            sigxy = sigxy + shear*dtheta*sin(thi)
             sigx  = 0.0
             sigxy = shear
         end do

         sigy = sigx

         if (lpow) then
            ipow = 2
         else
            ipow = 1
         end if

         sig(1) = B*((eps(1)+eps(2))-ep)**ipow + sigx
         sig(2) = B*((eps(1)+eps(2))-ep)**ipow - sigy
         sig(3) = sigxy

         sigmax  = sig(1)
         sigmay  = sig(2)
         sigmaxy = sig(3)

         return
      end if

!      Previous values for pore pressure estimation

      exx    = pvalp(1)
      eyy    = pvalp(2)
      exy    = pvalp(3)
      S      = pvalp(4)
      S0     = pvalp(5)
      Sb     = pvalp(6)
      ws     = pvalp(7)
      taup   = pvalp(8)
      sxx    = pvalp(9)
      syy    = pvalp(10)
      sxy    = pvalp(11)

!      plastic deformation

      if (lpow) then
         ipow = 2
         ep  = (poro/Kf)*(1.0-S)*sigma0 - sqrt(sigma0*S/B) + ee0
      else
         ipow = 1
         ep  = (poro/Kf)*(1.0-S)*sigma0 - sigma0*S/B + ee0
      end if

!      initial shear stress and shear modulus

      if (S0 .gt. Sb) then
         T0 = To*S
         G0 = (To/eo)*S
      else
         T0 = To*S + (m1-m2)*(Sb-S0)*(0.4/Sb)*sigma0
         e0 = eo/(S0/Sb)
         G0 = T0/e0
      end if

!      Multi-spring model

      Qv    = T0/2.0
      gv    = (T0/G0)*atan(1.0)

      sigx  = 0.0
      sigy  = 0.0
      sigxy = 0.0
      shear = 0.0

      if (d0 .le. 0.0) call curve(Nspr,gv,Qv,H,RR)
      do imec=1,IM
         thi = dtheta*(imec-1)
         xt  = (eps(1)-eps(2))*cos(thi) + eps(3)*sin(thi)
         if (d0 .le. 0.0) then
            call Iwan(Nspr,H,RR,alf(:,imec),pbar(:,imec),xt,shear)
         else
            call smulti(ix,d0,lkh,xt,dt,Qv,gv,         
     &                  nu,pval(:,imec),S0,S,shear)
         end if
         sigx  = sigx  + shear*dtheta*cos(thi)
         sigxy = sigxy + shear*dtheta*sin(thi)
      end do

      sigy = sigx

      sig(1) = B*((eps(1)+eps(2))-ep)**ipow + sigx
      sig(2) = B*((eps(1)+eps(2))-ep)**ipow - sigy
      sig(3) = sigxy

!      Updating ws, S0, and ep

      r = taup/sigma0

      ssl = 0.4 + (Sb-0.4)*S0/Sb
      if (S .ge. ssl) then
         if (r/S0 .le. m3) then
            Rc = 1.0
         else
            Rc = (m1-r/S)/(m1-m3)
         end if
      else
         if (r .le. ssl*m3) then
            Rc = 1.0
         else
            Rc = (ssl*m1 - r)/(ssl*(m1-m3))
         end if
      end if

      dexey  = (eps(1)-eps(2))-(exx-eyy)
      dexy   = eps(3) - exy
      dwst   = abs(0.5*(sxx-syy)*dexey + sxy*dexy)
      taunew = sqrt(sig(3)**2 + (sig(1)-sig(2))**2 / 4.0)
      sold   = (sxx+syy)/2.0
      snew   = (sig(1)+sig(2))/2.0
      if (lpow) then
         dGop   = (To/eo)*sqrt(sold/sigma0)
         dGo    = (To/eo)*sqrt(snew/sigma0)
      else
         dGop   = (To/eo)*S
         dGo    = (To/eo)*S
      end if
      if (dGop .lt. 1.0E-15 .or. dGo .lt. 1.0E-15) then
         dwse = 0.0
      else
         dwse  = abs(taup*(taunew/dGo-taup/dGop))
      end if
      dws   = dwst - c1*dwse

      if (dws .gt. 0.0) then
         ws = ws + Rc*dws
      end if
      w = ws/w0

      if (w .le. 0.0) then
         S0 = 1.0
      elseif (w .lt. w1) then
         S0 = 1.0-0.6*(w/w1)**p1
      elseif (w .gt. w1) then
         S0 = (0.4-S1)*(w1/w)**p2+S1
      end if

      r  = taunew/sigma0
      r2 = m2*S0
      r3 = m3*S0
      S2 = S0 - (r2-r3)/m1

      if (r .le. r3) S = S0
      if (r .gt. r3) S = S2 + sqrt((S0-S2)**2 + ((r-r3)/m1)**2)

!      Updating the new values for pore pressure

      pvalp(1)  = eps(1)            !exx
      pvalp(2)  = eps(2)            !eyy
      pvalp(3)  = eps(3)            !exy
      pvalp(4)  = S                 !S
      pvalp(5)  = S0                !S0
      pvalp(6)  = Sb                !Sb
      pvalp(7)  = ws                !ws
      pvalp(8)  = taunew            !tau
      pvalp(9)  = sig(1)            !sxx
      pvalp(10) = sig(2)            !syy
      pvalp(11) = sig(3)            !sxy

      sigmax  = sig(1)
      sigmay  = sig(2)
      sigmaxy = sig(3)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine curve(Nspr,eo,To,H,R)

      implicit none

      real    Gmax, eo, To
      real    x0, xu, dx
      integer i, Nspr

      real H(Nspr), R(Nspr), gamma(Nspr), G(Nspr)

      Gmax = To / eo

      x0 = -6.0
      xu = 1.0
      dx = (xu-x0)/Nspr

      do i=1,Nspr
         gamma(i) = 10.0**(x0 + dx*(i-1))
         G(i)     = 1.0 / (1.0 + gamma(i)/eo)
         R(i)     = Gmax * G(i) * gamma(i)
      end do

      H(1)    = Gmax
      H(Nspr) = 0.0
      do i=2,Nspr-1
         H(i) = (G(i+1)*gamma(i+1) - G(i)*gamma(i)) /   
     &            (gamma(i+1) - gamma(i)) * Gmax
      end do

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine curve_data(Nspr,Gmax,Npg,strain,GoverG,H,R)

      implicit none

      real    Gmax
      real    x0, xu, dx
      integer i, Nspr, Npg

      real H(Nspr), R(Nspr), gamma(Nspr), G(Nspr)
      real strain(Npg), GoverG(Npg)

      x0 = -6.0
      xu = 1.0
      dx = (xu-x0)/Nspr

      do i=1,Nspr
         gamma(i) = 10.0**(x0 + dx*(i-1))
         call iahsmbtp(strain,GoverG,Npg,gamma(i),G(i))
         R(i)     = Gmax * G(i) * gamma(i)
      end do

      H(1)    = Gmax
      H(Nspr) = 0.0
      do i=2,Nspr-1
         H(i) = (G(i+1)*gamma(i+1) - G(i)*gamma(i)) /   
     &            (gamma(i+1) - gamma(i)) * Gmax
      end do

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Iwan(Nspr,H,R,alf,pbar,exy,sxy)

      implicit none

      real    dx, dsxy, exyp, sxyp, sxy, x, exy
      integer i, j, Nspr

      real H(Nspr), R(Nspr), alf(Nspr), pbar(2)

      exyp = pbar(1)
      sxyp = pbar(2)

      dx   = exy - exyp
      if (dx .gt. 0.0) then
         x = 1.0
      else
         x = -1.0
      end if

      do i=1,Nspr
         dsxy = H(i)*dx
         if (abs(sxyp+dsxy-alf(i)) .le. R(i)) then
            sxy = sxyp + dsxy
            goto 3
         end if
         dsxy = alf(i) + x*R(i) - sxyp
         sxyp = sxyp + dsxy
         dx   = dx - dsxy/H(i)
      end do

3      continue

      if (i .gt. Nspr) i = Nspr
      if (abs(sxy-alf(i)) .lt. R(i) .or. i .eq. Nspr) i = i-1

      do j=1,i
         alf(j) = sxy - x*R(j)
      end do

      pbar(1) = exy
      pbar(2) = sxy

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine iahsmbtp(x,y,n,xi,yi)
!	A PROGRAM FOR INTERPOLATION BY ARC HYPERBOLIC
!	SINE METHOD BETWEEN TWO POINTS
!	x--ABSCISSA VALUE OF INPUT POINTS;  y--ORDINATE
! 	VALUE OF INPUT POINTS;  xi--INPUT ABSCISSA VALUE;
!	yi--OUTPUT ORDINATE VALUE

        dimension x(n),y(n)
        if(x(1)-1.e-10.lt.1.e-10) x(1)=1.e-10
        if(xi.le.x(1)) then
        yi=y(1)
        return
        else
        j=2
        endif
   30   if(xi.lt.x(j)) then
        ars=alog(xi*1.e+4+((xi*1.e+4)**2+1.)**0.5)
        ars1=alog(x(j-1)*1.e+4+((x(j-1)*1.e+4)**2+1.)**0.5)
        ars2=alog(x(j)*1.e+4+((x(j)*1.e+4)**2+1.)**0.5)
        a=(y(j)*x(j)*1.e+4-y(j-1)*x(j-1)*1.e+4)/(ars2-ars1)
        b=(y(j)*x(j)*1.e+4*ars1-y(j-1)*x(j-1)*1.e+4*ars2)/
     1  (ars1-ars2)
        yi=a/xi/1.e+4*ars+b/xi/1.e+4
        return
        elseif(xi.eq.x(j)) then
        yi=y(j)
        return
        else
        j=j+1
        endif
        if(j.le.n) then
        goto 30
        else
        yi=y(n)
        endif
        return
        end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smulti(ix,d0,lkh,xt,dt,Qv,gv, 
     &                  nu,pval,S0,S,shear)

!      Implementation of the multishear mechanism for each i-th spring

      implicit none

      real    tol

      parameter(tol=1.0E-6)

      real pval(7)

      integer j, ix
      real    nu, d0, xt, dt, Qv, gv, S0, S, shear
      real    xtp, etjp, sttjp, rc, avalue, bvalue
      real    er, etj, sttj, ef, red, sig, x1
      logical lkh

!      Previous values for the GMR

      xtp      = pval(1)
      x1       = pval(2)
      etj      = pval(3)
      sttj     = pval(4)
      bvalue   = pval(5)
      rc       = pval(6)
      j        = int(pval(7)+0.1)

      ef     = abs(x1)
      etjp   = etj
      sttjp  = sttj
      avalue = (bvalue + ef/gv) / (1.0 + ef/gv)
      er     = (xt-xtp)/dt

      if (j .eq. 1) then
         if (tol .ge. abs(xtp/gv)-abs(xt/gv)) then
            goto 70
         else
            goto 15
         end if
      else
         if (tol .ge. abs(x1/gv)-abs(xt/gv)) goto 70
         if (mod(j,2) .eq. 0 .and.   
     &       abs(xtp/gv-(-x1/gv))+tol .ge. abs(xt/gv-(-x1/gv)) ) goto 70
         if (mod(j,2) .eq. 1 .and.  
     &       abs(xtp/gv-(x1/gv))+tol .ge. abs(xt/gv-(x1/gv)) ) goto 70
         goto 15
      end if

!      Reversal point

15    j      = j+1
      etj    = xtp
      if (j .eq. 2) then
         ef = abs(etj)
         x1 = etj
      end if   
      red = etj - etjp
      call Fbb(Qv,gv,avalue,bvalue,rc,nu,er,red,sig)
      sttj = sttjp + sig

70    call CH(ix,j,Qv,gv,ef,etj,nu,er,sttj,avalue,bvalue,rc)
      call Parnew(j,d0,gv,ef,rc,lkh,S0,S,avalue,bvalue)

      red = xt - etj
      call Fbb(Qv,gv,avalue,bvalue,rc,nu,er,red,sig)
      shear = sttj + sig

!      Weak Masing rule

      if (abs(xt/gv) .ge. abs(ef/gv)) then
         j      = 1
         etj    = 0.0
         sttj   = 0.0
         rc     = 1.0
         bvalue = 1.0
         x1     = etj
      end if

!      Updating values of pval(7)

      pval(1) = xt
      pval(2) = x1
      pval(3) = etj
      pval(4) = sttj
      pval(5) = bvalue
      pval(6) = rc
      pval(7) = int(j)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CH(ix,j,Qv,gv,ef,et,nu,er1,sr,avalue,bvalue,rc)

!      This subroutine computes the hysteresis scale factor
!      for the generalized masing rules

      implicit none

      integer j, ix
      real    er, nu, Qv, gv, ef, et, er1, sr, avalue, bvalue, rc
      real    eps, sig, sf, dummy1, dummy2

      avalue = 1.0
      bvalue = 1.0

      er  = sign(1.0,er1)
      eps = er*ef
      call Fbb(Qv,gv,avalue,bvalue,1.0,nu,er1,eps,sig)
      sf  = sig

      if (j .eq. 1 .or. (eps-et) .eq. 0.0) then
         rc = 1.0
         return
      elseif (j .ge. 2) then
           dummy1 = (sf - sr) * abs(eps - et)
           dummy2 = Qv*(eps - et) - gv*(sf - sr)
         if (abs(dummy2) .eq. 0.0) then
            rc = -999.0
         else
            rc = dummy1 / dummy2
         end if
      end if

      if (rc .lt. 0.0 .and. rc .ne. -999.0) rc = 0.0

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Parnew(j,d0,e0,ef,rc,lkh,S0,S,avalue,bvalue)

!      This subroutine computes the values of a and b to
!      scale the stress and strain, respectively, and
!      control the area of the hysteresis loop when a
!      maximum damping ratio, d0, is specified.

      implicit none

      integer j, iflag
      real    d0, e0, ef, rc, S0, S, avalue, bvalue
      logical lkh

      iflag = 0

      if (j .eq. 1 .or. rc .le. 0.0) then
         bvalue = 1.0
         goto 1
      end if

      if (lkh) then
         call dgmr(j,0.1,e0,ef,rc,iflag)
      else
         call dgmr(j,d0,e0,ef,rc,iflag)
      end if

      if (iflag .eq. 1) then
         bvalue = 1.0
         goto 1
      end if

      if (lkh) then
         if (rc .ge. 2.0) then
            bvalue = sqrt(d0/rc)
         else
            bvalue = 1.0
         end if
      else
         call nr(rc,d0,ef,e0,1.0,bvalue)
      end if

1     avalue = (bvalue + ef/e0)/(1.0 + ef/e0)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Damp(rc,d0,ef,e0,b,chi,chip)

!      This subroutine computes the function and its
!      derivative for the newton-raphson method.  The
!      function in this case is the difference between
!      the hyperbolic damping function and the damping
!      for the generalized masing rules.

      implicit none

      real             rc, d0, ef, e0, b, chi, chip
      double precision pi,rcc,d00,x,bb,chi1,chip1

      pi = 4.0D0*datan(1.0D0)

      rcc = dble(rc)
      d00 = dble(d0)
      x   = dble(ef/e0)
      bb  = dble(b)

      chi1 = -((d00*x)/(1.0D0 + x)) +                            
     &         ((1.0D0 + bb/x)*                                    
     &         (-2.0D0/(bb/x + 2.0D0/rcc) + 2.0D0*rcc -            
     &         (bb*rcc**2 * dlog(1.0D0 + 2.0D0*x/(bb*rcc)))/x))/pi

      chi = sngl(chi1)

      chip1 = -((rcc*(-2.0D0*x*(2.0D0*bb**2 * rcc**2 +           
     &          bb*rcc*(6.0D0 + rcc)*x +                           
     &          (2.0D0 + 3.0D0*rcc)*x**2) +                        
     &          rcc*(2.0D0*bb + x)*(bb*rcc + 2.0D0*x)**2*          
     &          dlog(1.0D0 + 2.0D0*x/(bb*rcc))))/                  
     &          (pi*x**2*(bb*rcc + 2.0D0*x)**2))

      chip = sngl(chip1)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nr(rc,d0,ef,e0,x,x_new)

!      This subroutine performs the newton-raphson method
!      to obtain the value of b (strain scaling factor).

      implicit none

      integer Iter
      real    eps

      parameter(Iter=20, eps=1.0E-2)

      integer j
      real    rc, d0, ef, e0, x, x_new, x_old
      real    value_fun, deriv_fun, difference
      real    xlambda, diff2

      x_old = x

      xlambda = 1.0
      diff2   = 1.0E10

      do 2 j=1,Iter
         call Damp(rc,d0,ef,e0,x_old,value_fun,deriv_fun)
         if (deriv_fun .eq. 0.) then
            write(*,*) ' warning ', value_fun, deriv_fun
         end if
         x_new = x_old - xlambda*(value_fun/deriv_fun)
         difference = ABS(x_new-x_old)/ABS(x_new)
         if (x_old .lt. 0.0) then
            write(*,4) j,ef/e0,rc,x_old,x_new,             
     &                 value_fun,deriv_fun,difference
4           format(1x,i2,1x,7(e8.2,1x))
         end if
         if (difference .lt. eps) then
            goto 21
         elseif (difference .lt. diff2) then
            diff2 = difference
            goto 633
         elseif (difference .gt. eps .and. difference .gt. diff2) then
            if (xlambda .gt. 0.0001) then
               xlambda = xlambda/10.0
            else
               xlambda = 0.0001
            end if
            x_new = x_old - xlambda*(value_fun/deriv_fun)
            goto 633
         elseif (j .ge. Iter) then
            write(*,*) 'NoC',j,ef/e0,rc,x_old,x_new,difference
            stop
         end if
633      x_old = x_new
2     continue

21    continue

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dgmr(j,d0,e0,ef,rc,iflag)

!      This subroutine computes the conditions when the scaling
!      factors for the stress and the strain should be computed.

      implicit none

      integer iflag, j
      real    rc, d0, ef, e0, fmax, deriv_fun

      call Damp(rc,d0,ef,e0,1.0,fmax,deriv_fun)

      if (-d0*fmax .lt. 0.0 .and. rc .gt. 1.0) then
         iflag = 0
      else
         iflag = 1
      end if

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C XAPIIR -- SUBROUTINE:   IIR FILTER DESIGN AND IMPLEMENTATION                  
C                                                                               
C  AUTHOR:  Dave Harris                                                         
C                                                                               
C  LAST MODIFIED:  September 12, 1990                                           
C                                                                               
C  ARGUMENTS:                                                                   
C  ----------                                                                   
C                                                                               
C    DATA           REAL ARRAY CONTAINING SEQUENCE TO BE FILTERED               
C                     ORIGINAL DATA DESTROYED, REPLACED BY FILTERED DATA        
C                                                                               
C    NSAMPS         NUMBER OF SAMPLES IN DATA                                   
C                                                                               
C                                                                               
C    APROTO         CHARACTER*2 VARIABLE, CONTAINS TYPE OF ANALOG               
C                     PROTOTYPE FILTER                                          
C                     '(BU)TTER  ' -- BUTTERWORTH FILTER                        
C                     '(BE)SSEL  ' -- BESSEL FILTER                             
C                     'C1      ' -- CHEBYSHEV TYPE I                            
C                     'C2      ' -- CHEBYSHEV TYPE II                           
C                                                                               
C    TRBNDW         TRANSITION BANDWIDTH AS FRACTION OF LOWPASS                 
C                   PROTOTYPE FILTER CUTOFF FREQUENCY.  USED                    
C                   ONLY BY CHEBYSHEV FILTERS.                                  
C                                                                               
C    A              ATTENUATION FACTOR.  EQUALS AMPLITUDE                       
C                   REACHED AT STOPBAND EDGE.  USED ONLY BY                     
C                   CHEBYSHEV FILTERS.                                          
C                                                                               
C    IORD           ORDER (#POLES) OF ANALOG PROTOTYPE                          
C                   NOT TO EXCEED 10 IN THIS CONFIGURATION.  4 - 5              
C                   SHOULD BE AMPLE.                                            
C                                                                               
C    TYPE           CHARACTER*2 VARIABLE CONTAINING FILTER TYPE                 
C                     'LP' -- LOW PASS                                          
C                     'HP' -- HIGH PASS                                         
C                     'BP' -- BAND PASS                                         
C                     'BR' -- BAND REJECT                                       
C                                                                               
C    FLO            LOW FREQUENCY CUTOFF OF FILTER (HERTZ)                      
C                   IGNORED IF TYPE = 'LP'                                      
C                                                                               
C    FHI            HIGH FREQUENCY CUTOFF OF FILTER (HERTZ)                     
C                   IGNORED IF TYPE = 'HP'                                      
C                                                                               
C    TS             SAMPLING INTERVAL (SECONDS)                                 
C                                                                               
C    PASSES           INTEGER VARIABLE CONTAINING THE NUMBER OF PASSES          
C                   1 -- FORWARD FILTERING ONLY                                 
C                   2 -- FORWARD AND REVERSE (I.E. ZERO PHASE) FILTERING        
C                                                                               
C                                                                               
C  SUBPROGRAMS REFERENCED:  BILIN2, BUROOTS, WARP, CUTOFFS, LPTHP, LPTBP,       
C    LP, LPTBR, BEROOTS, C1ROOTS, C2ROOTS, CHEBPARM, DESIGN, APPLY              
C                                                                               
      SUBROUTINE XAPIIR( DATA, NSAMPS, APROTO, TRBNDW, A, IORD,         
     +                   TYPE, FLO, FHI, TS, PASSES, MAX_NT)                                                      
C    
        DIMENSION DATA(MAX_NT)                                               
        CHARACTER*2 TYPE, APROTO                                        
        INTEGER NSAMPS, PASSES, IORD                                    
        REAL*4 TRBNDW, A, FLO, FHI, TS, SN(30), SD(30)                  
        LOGICAL ZP                                                      
C                                                                               
C  Filter designed                                                              
C                                                                               
        CALL DESIGN( IORD, TYPE(1:2), APROTO(1:2), A, TRBNDW,           
     &               FLO, FHI, TS, SN, SD, NSECTS )                     
C                                                                               
C  Filter data                                                                  
C                                                                               
        IF (   PASSES .EQ. 1 ) THEN                                     
          ZP = .FALSE.                                                  
        ELSE                                                            
          ZP = .TRUE.                                                   
        END IF                                                          
        CALL APPLY( DATA, NSAMPS, ZP, SN, SD, NSECTS, MAX_NT)                  
C                                                                               
      RETURN                                                            
      END                                                                       
C                                                               APPLY           
C  Subroutine to apply an iir filter to a data sequence.                        
C    The filter is assumed to be stored as second order sections.               
C    Filtering is in-place.                                                     
C    Zero-phase (forward and reverse) is an option.                             
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    DATA                           Array containing data                       
C                                                                               
C    NSAMPS                         Number of data samples                      
C                                                                               
C    ZP                             Logical variable, true for                  
C                                     zero phase filtering, false               
C                                     for single pass filtering                 
C                                                                               
C    SN                             Numerator polynomials for second            
C                                     order sections.                           
C                                                                               
C    SD                             Denominator polynomials for second          
C                                     order sections.                           
C                                                                               
C    NSECTS                         Number of second-order sections             
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    DATA                          Data array (same as input)                   
C                                                                               
C                                                                               
      SUBROUTINE APPLY( DATA, NSAMPS, ZP, SN, SD, NSECTS, MAX_NT)              
C                                                                               
        REAL*4 SN(1), SD(1), DATA(MAX_NT)                                    
        REAL*4 OUTPUT                                                   
        LOGICAL ZP                                                      
C                                                                               
        JPTR = 1                                                        
        DO    1 J = 1, NSECTS                                           
C                                                                               
          X1 = 0.0                                                      
          X2 = 0.0                                                      
          Y1 = 0.0                                                      
          Y2 = 0.0                                                      
          B0 = SN(JPTR)                                                 
          B1 = SN(JPTR+1)                                               
          B2 = SN(JPTR+2)                                               
          A1 = SD(JPTR+1)                                               
          A2 = SD(JPTR+2)                                               
C                                                                               
          DO    2 I = 1, NSAMPS                                         
C                                                                               
            OUTPUT = B0*DATA(I) + B1*X1 + B2*X2                         
            OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )                         
            Y2 = Y1                                                     
            Y1 = OUTPUT                                                 
            X2 = X1                                                     
            X1 = DATA(I)                                                
            DATA(I) = OUTPUT                                            
C                                                                               
    2     CONTINUE                                                      
C                                                                               
          JPTR = JPTR + 3                                               
C                                                                               
    1   CONTINUE                                                        
C                                                                               
        IF (   ZP ) THEN                                                
C                                                                               
          JPTR = 1                                                      
          DO    3 J = 1, NSECTS                                         
C                                                                               
            X1 = 0.0                                                    
            X2 = 0.0                                                    
            Y1 = 0.0                                                    
            Y2 = 0.0                                                    
            B0 = SN(JPTR)                                               
            B1 = SN(JPTR+1)                                             
            B2 = SN(JPTR+2)                                             
            A1 = SD(JPTR+1)                                             
            A2 = SD(JPTR+2)                                             
C                                                                               
            DO    4 I = NSAMPS, 1, -1                                   
C                                                                               
              OUTPUT = B0*DATA(I) + B1*X1 + B2*X2                       
              OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )                       
              Y2 = Y1                                                   
              Y1 = OUTPUT                                               
              X2 = X1                                                   
              X1 = DATA(I)                                              
              DATA(I) = OUTPUT                                          
C                                                                               
    4       CONTINUE                                                    
C                                                                               
            JPTR = JPTR + 3                                             
C                                                                               
    3     CONTINUE                                                      
C                                                                               
        END IF                                                          
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C  DESIGN -- Subroutine to design IIR digital filters from analog              
C    prototypes.
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    IORD                Filter order (10 MAXIMUM)                              
C                                                                               
C    TYPE                Character*2 variable containing filter type            
C                          LOWPASS (LP)                                         
C                          HIGHPASS (HP)                                        
C                          BANDPASS (BP)                                        
C                          BANDREJECT (BR)                                      
C                                                                               
C    APROTO              Character*2 variable designating analog prototype      
C                          Butterworth (BU)                                     
C                          Bessel (BE)                                          
C                          Chebyshev Type I (C1)                                
C                          Chebyshev Type II (C2)                               
C                                                                               
C    A                   Chebyshev stopband attenuation factor                  
C                                                                               
C    TRBNDW              Chebyshev transition bandwidth (fraction of            
C                          lowpass prototype passband width)                    
C                                                                               
C    FL                  Low-frequency cutoff                                   
C                                                                               
C    FH                  High-frequency cutoff                                  
C                                                                               
C    TS                  Sampling interval (in seconds)                         
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                  Array containing numerator coefficients of             
C                        second-order sections packed head-to-tail.             
C                                                                               
C    SD                  Array containing denominator coefficients              
C                        of second-order sections packed head-to-tail.          
C                                                                               
C    NSECTS              Number of second-order sections.                       
C                                                                               
C                                                                               
      SUBROUTINE DESIGN( IORD, TYPE, APROTO, A, TRBNDW,                 
     &                   FL, FH, TS, SN, SD, NSECTS )                   
C                                                                               
        COMPLEX P(10), Z(10)                                            
        CHARACTER*2 TYPE, APROTO                                        
        CHARACTER*3 STYPE(10)                                           
        REAL*4 SN(1), SD(1)                                             
C                                                                               
C  Analog prototype selection                                                   
C                                                                               
        IF (     APROTO .EQ. 'BU' ) THEN                                
C                                                                               
          CALL BUROOTS( P, STYPE, DCVALUE, NSECTS, IORD )               
C                                                                               
        ELSE IF (    APROTO .EQ. 'BE' ) THEN                            
C                                                                               
          CALL BEROOTS( P, STYPE, DCVALUE, NSECTS, IORD )               
C                                                                               
        ELSE IF (    APROTO .EQ. 'C1' ) THEN                            
C                                                                               
          CALL CHEBPARM( A, TRBNDW, IORD, EPS, RIPPLE )                 
          CALL C1ROOTS( P, STYPE, DCVALUE, NSECTS, IORD, EPS )          
C                                                                               
        ELSE IF (    APROTO .EQ. 'C2' ) THEN                            
C                                                                               
          OMEGAR = 1. + TRBNDW                                          
          CALL C2ROOTS( P, Z, STYPE, DCVALUE, NSECTS, IORD, A, OMEGAR ) 
C                                                                               
        END IF                                                          
C                                                                               
C  Analog mapping selection                                                     
C                                                                               
        IF (     TYPE .EQ. 'BP' ) THEN                                  
C                                                                               
          FLW = WARP( FL*TS/2., 2. )                                    
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LPTBP( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )  
C                                                                               
        ELSE IF (   TYPE .EQ. 'BR' ) THEN                               
C                                                                               
          FLW = WARP( FL*TS/2., 2. )                                    
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LPTBR( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )  
C                                                                               
        ELSE IF (    TYPE .EQ. 'LP' ) THEN                              
C                                                                               
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )               
          CALL CUTOFFS( SN, SD, NSECTS, FHW )                           
C                                                                               
        ELSE IF (    TYPE .EQ. 'HP' ) THEN                              
C                                                                               
          FLW = WARP( FL*TS/2., 2. )                                    
          CALL LPTHP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )            
          CALL CUTOFFS( SN, SD, NSECTS, FLW )                           
C                                                                               
        END IF                                                          
C                                                                               
C  Bilinear analog to digital transformation                                    
C                                                                               
        CALL BILIN2( SN, SD, NSECTS )                                   
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C BUROOTS -- SUBROUTINE TO COMPUTE BUTTERWORTH POLES FOR                        
C   NORMALIZED LOWPASS FILTER                                                   
C                                                                               
C LAST MODIFIED:  SEPTEMBER 7, 1990                                             
C                                                                               
C  OUTPUT ARGUMENTS:                                                            
C  -----------------                                                            
C      P              COMPLEX ARRAY CONTAINING POLES                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL POLES                                          
C                                                                               
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
C                       TYPE:                                                   
C                         (SP)  SINGLE REAL POLE                                
C                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
C                                                                               
C      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
C                                                                               
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
C                                                                               
C  INPUT ARGUMENTS:                                                             
C  ----------------                                                             
C                                                                               
C      IORD           DESIRED FILTER ORDER                                      
C                                                                               
C                                                                               
      SUBROUTINE BUROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )             
C                                                                               
        COMPLEX P(1)                                                    
        INTEGER HALF                                                    
        CHARACTER*3 RTYPE(1)                                            
C                                                                               
        PI=3.14159265                                                   
C                                                                               
        HALF = IORD/2                                                   
C                                                                               
C TEST FOR ODD ORDER, AND ADD POLE AT -1                                        
C                                                                               
        NSECTS = 0                                                      
        IF (    2*HALF .LT. IORD ) THEN                                 
          P(1) = CMPLX( -1., 0. )                                       
          RTYPE(1) = 'SP'                                               
          NSECTS = 1                                                    
        END IF                                                          
C                                                                               
        DO    1  K = 1, HALF                                            
          ANGLE = PI * ( .5 + FLOAT(2*K-1) / FLOAT(2*IORD) )            
          NSECTS = NSECTS + 1                                           
          P(NSECTS) = CMPLX( COS(ANGLE), SIN(ANGLE) )                   
          RTYPE(NSECTS) = 'CP'                                          
    1   CONTINUE                                                        
C                                                                               
        DCVALUE = 1.0                                                   
C                                                                               
      RETURN                                                            
      END  
C                                                                               
C BEROOTS -- SUBROUTINE TO RETURN BESSEL POLES FOR                              
C   NORMALIZED LOWPASS FILTER                                                   
C                                                                               
C LAST MODIFIED:  April 15, 1992. Changed P and RTYPE to adjustable 
C                 array by using an "*" rather than a "1".     
C                                                                               
C  OUTPUT ARGUMENTS:                                                            
C  -----------------                                                            
C      P              COMPLEX ARRAY CONTAINING POLES                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL POLES                                          
C                                                                               
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
C                       TYPE:                                                   
C                         (SP)  SINGLE REAL POLE                                
C                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
C                                                                               
C      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
C                                                                               
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
C                                                                               
C  INPUT ARGUMENTS:                                                             
C  ----------------                                                             
C                                                                               
C      IORD           DESIRED FILTER ORDER                                      
C                                                                               
C                                                                               
      SUBROUTINE BEROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )             
C                                                                               
        COMPLEX P(*)                                                    
        INTEGER NSECTS, IORD                                            
        CHARACTER*3 RTYPE(*)                                            
C                                                                               
        IF (   IORD .EQ. 1 ) THEN                                       
C                                                                               
          P(1) = CMPLX( -1.0, 0.0 )                                     
          RTYPE(1) = 'SP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 2 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -1.1016013,  0.6360098 )                        
          RTYPE(1) = 'CP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 3 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -1.0474091, 0.9992645 )                         
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3226758, 0.0 )                               
          RTYPE(2) = 'SP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 4 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.9952088,  1.2571058 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3700679, 0.4102497 )                         
          RTYPE(2) = 'CP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 5 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.9576766,  1.4711244 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3808774,  0.7179096 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.5023160, 0.0 )                               
          RTYPE(3) = 'SP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 6 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.9306565,  1.6618633 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3818581,  0.9714719 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.5714904,  0.3208964 )                        
          RTYPE(3) = 'CP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 7 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.9098678,  1.8364514 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3789032,  1.1915667 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.6120388,  0.5892445 )                        
          RTYPE(3) = 'CP'                                               
          P(4) = CMPLX( -1.6843682, 0.0 )                               
          RTYPE(4) = 'SP'                                               
C                                                                               
        ELSE IF (  IORD .EQ. 8 ) THEN                                   
C                                                                               
          P(1) = CMPLX( -0.8928710,  1.9983286 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3738431,  1.3883585 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.6369417,  0.8227968 )                        
          RTYPE(3) = 'CP'                                               
          P(4) = CMPLX( -1.7574108,  0.2728679 )                        
          RTYPE(4) = 'CP'                                               
C                                                                               
        END IF                                                          
C                                                                               
        NSECTS = IORD - IORD/2                                          
C                                                                               
        DCVALUE = 1.0                                                   
C                                                                               
C  DONE                                                                         
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C  CHEBPARM - Calculates Chebyshev type I and II design parameters              
C                                                                               
C                                                                               
C  INPUT ARGUMENTS                                                              
C  ---------------                                                              
C                                                                               
C       A                Desired stopband attenuation                           
C                          i.e. max stopband amplitude is 1/ATTEN               
C                                                                               
C       TRBNDW           Transition bandwidth between stop and passbands        
C                          as a fraction of the passband width                  
C                                                                               
C       IORD             Filter order (number of poles)                         
C                                                                               
C                                                                               
C  OUTPUT ARGUMENTS                                                             
C  ----------------                                                             
C                                                                               
C       EPS              Chebyshev passband parameter                           
C                                                                               
C       RIPPLE           Passband ripple                                        
C                                                                               
      SUBROUTINE CHEBPARM( A, TRBNDW, IORD, EPS, RIPPLE )               
                                                                        
          OMEGAR  =  1. + TRBNDW                                        
          ALPHA = ( OMEGAR + SQRT( OMEGAR**2 - 1. ) ) ** IORD           
          G = ( ALPHA**2 + 1. ) / (2.*ALPHA)                            
          EPS = SQRT( A**2 - 1. ) / G                                   
          RIPPLE = 1. / SQRT( 1. + EPS**2 )                             
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C C1ROOTS -- SUBROUTINE TO COMPUTE CHEBYSHEV TYPE I POLES FOR                   
C   NORMALIZED LOWPASS FILTER                                                   
C                                                                               
C LAST MODIFIED:  SEPTEMBER 7, 1990                                             
C                                                                               
C  OUTPUT ARGUMENTS:                                                            
C  -----------------                                                            
C      P              COMPLEX ARRAY CONTAINING POLES                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL POLES                                          
C                                                                               
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
C                       TYPE:                                                   
C                         (SP)  SINGLE REAL POLE                                
C                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
C                                                                               
C      DCVALUE        RESPONSE OF FILTER AT ZERO FREQUENCY                      
C                                                                               
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
C                                                                               
C  INPUT ARGUMENTS:                                                             
C  ----------------                                                             
C                                                                               
C      IORD           DESIRED FILTER ORDER                                      
C                                                                               
C      EPS            CHEBYSHEV PARAMETER RELATED TO PASSBAND RIPPLE            
C                                                                               
      SUBROUTINE C1ROOTS( P, RTYPE, DCVALUE, NSECTS, IORD, EPS )        
C                                                                               
        COMPLEX P(1)                                                    
        INTEGER HALF                                                    
        CHARACTER*3 RTYPE(1)                                            
C                                                                               
        PI = 3.14159265                                                 
        HALF = IORD/2                                                   
C                                                                               
C  INTERMEDIATE DESIGN PARAMETERS                                               
C                                                                               
        GAMMA = ( 1. + SQRT( 1. + EPS*EPS ) ) / EPS                     
        GAMMA = ALOG(GAMMA) / FLOAT(IORD)                               
        GAMMA = EXP(GAMMA)                                              
        S = .5 * ( GAMMA - 1./GAMMA )                                   
        C = .5 * ( GAMMA + 1./GAMMA )                                   
C                                                                               
C  CALCULATE POLES                                                              
C                                                                               
        NSECTS = 0                                                      
        DO    1  I = 1 ,  HALF                                          
          RTYPE(I) = 'CP'                                               
          ANGLE = FLOAT(2*I-1) * PI/FLOAT(2*IORD)                       
          SIGMA = -S * SIN(ANGLE)                                       
          OMEGA =  C * COS(ANGLE)                                       
          P(I) = CMPLX( SIGMA, OMEGA )                                  
          NSECTS = NSECTS + 1                                           
    1   CONTINUE                                                        
        IF (   2*HALF .LT. IORD ) THEN                                  
          RTYPE( HALF + 1 ) = 'SP'                                      
          P(HALF+1) = CMPLX( -S, 0.0 )                                  
          NSECTS = NSECTS + 1                                           
          DCVALUE = 1.0                                                 
        ELSE                                                            
          DCVALUE = 1./SQRT( 1 + EPS**2 )                               
        END IF                                                          
C                                                                               
C  DONE                                                                         
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C C2ROOTS -- SUBROUTINE TO COMPUTE ROOTS FOR NORMALIZED LOWPASS                 
C   CHEBYSHEV TYPE 2 FILTER                                                     
C                                                                               
C LAST MODIFIED:  SEPTEMBER 7, 1990                                             
C                                                                               
C  OUTPUT ARGUMENTS:                                                            
C  -----------------                                                            
C      P              COMPLEX ARRAY CONTAINING POLES                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL POLES                                          
C                                                                               
C      Z              COMPLEX ARRAY CONTAINING ZEROS                            
C                       CONTAINS ONLY ONE FROM EACH                             
C                       COMPLEX CONJUGATE PAIR, AND                             
C                       ALL REAL ZEROS                                          
C                                                                               
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
C                       TYPE:                                                   
C                         (SP)  SINGLE REAL POLE                                
C                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
C                                                                               
C      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
C                                                                               
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
C                                                                               
C  INPUT ARGUMENTS:                                                             
C  ----------------                                                             
C                                                                               
C                                                                               
C      IORD           DESIRED FILTER ORDER                                      
C                                                                               
C      A              STOPBAND ATTENUATION FACTOR                               
C                                                                               
C      OMEGAR         CUTOFF FREQUENCY OF STOPBAND                              
C                     PASSBAND CUTOFF IS AT 1.0 HERTZ                           
C                                                                               
C                                                                               
      SUBROUTINE C2ROOTS( P, Z, RTYPE, DCVALUE, NSECTS, IORD, A, OMEGAR 
     1)                                                                 
C                                                                               
        COMPLEX P(1), Z(1)                                              
        INTEGER HALF                                                    
        CHARACTER*3 RTYPE(1)                                            
C                                                                               
        PI = 3.14159265                                                 
        HALF = IORD/2                                                   
C                                                                               
C  INTERMEDIATE DESIGN PARAMETERS                                               
C                                                                               
        GAMMA = (A+SQRT(A*A-1.))                                        
        GAMMA = ALOG(GAMMA)/FLOAT(IORD)                                 
        GAMMA = EXP(GAMMA)                                              
        S = .5*(GAMMA-1./GAMMA)                                         
        C = .5*(GAMMA+1./GAMMA)                                         
C                                                                               
        NSECTS = 0                                                      
        DO    1 I = 1, HALF                                             
C                                                                               
C  CALCULATE POLES                                                              
C                                                                               
          RTYPE(I) = 'CPZ'                                              
C                                                                               
          ANGLE = FLOAT(2*I-1) * PI/FLOAT(2*IORD)                       
          ALPHA = -S*SIN(ANGLE)                                         
          BETA = C*COS(ANGLE)                                           
          DENOM = ALPHA*ALPHA + BETA*BETA                               
          SIGMA = OMEGAR*ALPHA/DENOM                                    
          OMEGA = -OMEGAR*BETA/DENOM                                    
          P(I) = CMPLX( SIGMA, OMEGA )                                  
C                                                                               
C  CALCULATE ZEROS                                                              
C                                                                               
          OMEGA = OMEGAR/COS(ANGLE)                                     
          Z(I) = CMPLX( 0.0, OMEGA )                                    
C                                                                               
          NSECTS = NSECTS + 1                                           
C                                                                               
    1   CONTINUE                                                        
C                                                                               
C  ODD-ORDER FILTERS                                                            
C                                                                               
        IF (  2*HALF .LT. IORD ) THEN                                   
          RTYPE(HALF+1) = 'SP'                                          
          P(HALF+1) = CMPLX( -OMEGAR/S, 0.0 )                           
          NSECTS = NSECTS + 1                                           
        END IF                                                          
C                                                                               
C  DC VALUE                                                                     
C                                                                               
        DCVALUE = 1.0                                                   
C                                                                               
C  DONE                                                                         
C                                                                               
      RETURN                                                            
      END                                                               
C                                                                               
C WARP -- FUNCTION, APPLIES TANGENT FREQUENCY WARPING TO COMPENSATE             
C         FOR BILINEAR ANALOG -> DIGITAL TRANSFORMATION                         
C                                                                               
C ARGUMENTS:                                                                    
C ----------                                                                    
C                                                                               
C      F       ORIGINAL DESIGN FREQUENCY SPECIFICATION (HERTZ)                  
C      TS      SAMPLING INTERVAL (SECONDS)                                      
C                                                                               
C  LAST MODIFIED:  SEPTEMBER 20, 1990                                           
C                                                                               
      REAL FUNCTION WARP( F , TS )                                      
C                                                                               
        TWOPI = 6.2831853                                               
        ANGLE = TWOPI*F*TS/2.                                           
        WARP = 2.*TAN(ANGLE)/TS                                         
        WARP = WARP/TWOPI                                               
C                                                                               
      RETURN                                                            
      END      
C                                                                               
C  Subroutine to generate second order section parameterization                 
C    from an pole-zero description for lowpass filters.                         
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeros                             
C                                                                               
C    RTYPE                   Character array containing root type information   
C                              (SP)  Single real pole or                        
C                              (CP)  Complex conjugate pole pair                
C                              (CPZ) Complex conjugate pole and zero pairs      
C                                                                               
C    DCVALUE                 Zero-frequency value of prototype filter           
C                                                                               
C    NSECTS                  Number of second-order sections                    
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C                                                                               
      SUBROUTINE LP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )             
C                                                                               
        COMPLEX P(*), Z(*)                                              
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    
C                                                                               
        IPTR = 1                                                        
        DO    1 I = 1, NSECTS                                           
C                                                                               
          IF (   RTYPE(I) .EQ. 'CPZ' ) THEN                             
C                                                                               
            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
     &            / REAL( Z(I) * CONJG( Z(I) ) )                        
            SN( IPTR )     = REAL( Z(I) * CONJG( Z(I) ) ) * SCALE       
            SN( IPTR + 1 ) = -2. * REAL( Z(I) ) * SCALE                 
            SN( IPTR + 2 ) = 1. * SCALE                                 
            SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )               
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
            SN( IPTR )     = SCALE                                      
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )               
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SCALE = -REAL( P(I) )                                       
            SN( IPTR )     = SCALE                                      
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = -REAL( P(I) )                              
            SD( IPTR + 1 ) = 1.                                         
            SD( IPTR + 2 ) = 0.                                         
            IPTR = IPTR + 3                                             
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
        SN(1) = DCVALUE * SN(1)                                         
        SN(2) = DCVALUE * SN(2)                                         
        SN(3) = DCVALUE * SN(3)                                         
                                                                        
C                                                                               
      RETURN                                                            
      END
C                                                    LPTBP                      
C                                                                               
C  Subroutine to convert an prototype lowpass filter to a bandpass filter via   
C    the analog polynomial transformation.  The lowpass filter is               
C    described in terms of its poles and zeros (as input to this routine).      
C    The output consists of the parameters for second order sections.           
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeros                             
C                                                                               
C    RTYPE                   Character array containing type information        
C                              (SP) single real pole  or                        
C                              (CP) complex conjugate pole pair  or             
C                              (CPZ) complex conjugate pole/zero pairs          
C                                                                               
C    DCVALUE                 Zero frequency value of filter                     
C                                                                               
C    NSECTS                  Number of second-order sections upon input         
C                                                                               
C    FL                      Low-frequency cutoff                               
C                                                                               
C    FH                      High-frequency cutoff                              
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C    NSECTS                  Number of second order sections upon output        
C                              This subroutine doubles the number of            
C                              sections.                                        
C                                                                               
C                                                                               
      SUBROUTINE LPTBP( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )  
C                                                                               
        COMPLEX P(*), Z(*), CTEMP, P1, P2, Z1, Z2, S, H                 
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    
C                                                                               
        PI = 3.14159265                                                 
        TWOPI = 2.*PI                                                   
        A = TWOPI*TWOPI*FL*FH                                           
        B = TWOPI*( FH - FL )                                           
C                                                                               
        N = NSECTS                                                      
        NSECTS = 0                                                      
        IPTR = 1                                                        
        DO    1 I = 1, N                                                
C                                                                               
          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN                            
C                                                                               
            CTEMP = ( B*Z(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            Z1 = 0.5*( B*Z(I) + CTEMP )                                 
            Z2 = 0.5*( B*Z(I) - CTEMP )                                 
            CTEMP = ( B*P(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*P(I) + CTEMP )                                 
            P2 = 0.5*( B*P(I) - CTEMP )                                 
            SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z1 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z2 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            CTEMP = ( B*P(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*P(I) + CTEMP )                                 
            P2 = 0.5*( B*P(I) - CTEMP )                                 
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = A                                          
            SD( IPTR + 1 ) = -B*REAL( P(I) )                            
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 1                                         
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
C  Scaling - use the fact that the bandpass filter amplitude at sqrt( omega_l * 
C            equals the amplitude of the lowpass prototype at d.c.              
C                                                                               
        S = CMPLX( 0., SQRT(A) )                                        
        H = CMPLX( 1., 0. )                                             
C                                                                               
        IPTR = 1                                                        
        DO    2 I = 1, NSECTS                                           
          H = H * ( ( SN(IPTR+2)*S + SN(IPTR+1) )*S + SN(IPTR) )        
     &      / ( ( SD(IPTR+2)*S + SD(IPTR+1) )*S + SD(IPTR) )            
          IPTR = IPTR + 3                                               
    2   CONTINUE                                                        
        SCALE = DCVALUE / SQRT( REAL( H ) * CONJG( H ) )                
                                                                        
        SN(1) = SN(1) * SCALE                                           
        SN(2) = SN(2) * SCALE                                           
        SN(3) = SN(3) * SCALE                                           
C                                                                               
      RETURN                                                            
      END                                                               
C                                                    LPTBR                      
C                                                                               
C  Subroutine to convert a lowpass filter to a band reject filter               
C    via an analog polynomial transformation.  The lowpass filter is            
C    described in terms of its poles and zeros (as input to this routine).      
C    The output consists of the parameters for second order sections.           
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeros                             
C                                                                               
C    RTYPE                   Character array containing type information        
C                              (SP)  single real pole or                        
C                              (CP)  complex conjugate pole pair                
C                              (CPZ) complex conjugate pole/zero pairs          
C                                                                               
C    DCVALUE                 Zero-frequency value of prototype filter           
C                                                                               
C    NSECTS                  Number of second-order sections                    
C                              prior to transformation                          
C                                                                               
C    FL                      Low-frequency cutoff                               
C                                                                               
C    FH                      High-frequency cutoff                              
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C    NSECTS                  Number of second order sections following          
C                              transformation.  The number is doubled.          
C                                                                               
C                                                                               
      SUBROUTINE LPTBR( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )  
C                                                                               
        COMPLEX P(*), Z(*), CINV, CTEMP, P1, P2, Z1, Z2                 
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*)                                             
C                                                                               
        PI = 3.14159265                                                 
        TWOPI = 2.*PI                                                   
        A = TWOPI*TWOPI*FL*FH                                           
        B = TWOPI*( FH - FL )                                           
C                                                                               
        N = NSECTS                                                      
        NSECTS = 0                                                      
        IPTR = 1                                                        
        DO    1 I = 1, N                                                
C                                                                               
          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN                            
C                                                                               
            CINV = 1./Z(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            Z1 = 0.5*( B*CINV + CTEMP )                                 
            Z2 = 0.5*( B*CINV - CTEMP )                                 
            CINV = 1./P(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*CINV + CTEMP )                                 
            P2 = 0.5*( B*CINV - CTEMP )                                 
            SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z1 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z2 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            CINV = 1./P(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*CINV + CTEMP )                                 
            P2 = 0.5*( B*CINV - CTEMP )                                 
            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = -A*REAL( P(I) )                            
            SD( IPTR + 1 ) = B                                          
            SD( IPTR + 2 ) = -REAL( P(I) )                              
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 1                                         
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
C  Scaling - use the fact that the bandreject filter amplitude  at d.c.         
C            equals the lowpass prototype amplitude at d.c.                     
C                                                                               
        H = 1.0                                                         
C                                                                               
        IPTR = 1                                                        
        DO    2 I = 1, NSECTS                                           
          H = H * SN(IPTR) / SD(IPTR)                                   
          IPTR = IPTR + 3                                               
    2   CONTINUE                                                        
        SCALE = DCVALUE / ABS(H)                                        
        SN(1) = SN(1) * SCALE                                           
        SN(2) = SN(2) * SCALE                                           
        SN(3) = SN(3) * SCALE                                           
C                                                                               
      RETURN                                                            
      END                                                               
C                                                    LPTHP                      
C                                                                               
C  Subroutine to convert a lowpass filter to a highpass filter via              
C    an analog polynomial transformation.  The lowpass filter is                
C    described in terms of its poles and zeroes (as input to this routine).     
C    The output consists of the parameters for second order sections.           
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeroes                            
C                                                                               
C    RTYPE                   Character array containing root type information   
C                              (SP) single real pole or                         
C                              (CP)  complex conjugate pair                     
C                              (CPZ) complex pole/zero pairs                    
C                                                                               
C    DCVALUE                 Zero-frequency value of prototype filter           
C                                                                               
C    NSECTS                  Number of second-order sections                    
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C                                                                               
      SUBROUTINE LPTHP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )          
C                                                                               
        COMPLEX P(*), Z(*)                                              
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    
C                                                                               
        IPTR = 1                                                        
        DO    1 I = 1, NSECTS                                           
C                                                                               
          IF (     RTYPE(I) .EQ. 'CPZ' ) THEN                           
C                                                                               
            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
     &            / REAL( Z(I) * CONJG( Z(I) ) )                        
            SN( IPTR )     = 1.  *  SCALE                               
            SN( IPTR + 1 ) = -2. * REAL( Z(I) )  *  SCALE               
            SN( IPTR + 2 ) = REAL( Z(I) * CONJG( Z(I) ) )  *  SCALE     
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )               
            IPTR = IPTR + 3                                             
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = SCALE                                      
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )               
            IPTR = IPTR + 3                                             
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SCALE = -REAL( P(I) )                                       
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = SCALE                                      
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -REAL( P(I) )                              
            SD( IPTR + 2 ) = 0.                                         
            IPTR = IPTR + 3                                             
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
        SN(1) = SN(1) * DCVALUE                                         
        SN(2) = SN(2) * DCVALUE                                         
        SN(3) = SN(3) * DCVALUE                                         
C                                                                               
      RETURN                                                            
      END             
C                                                    CUTOFFS                    
C                                                                               
C  Subroutine to alter the cutoff of a filter.  Assumes that the                
C    filter is structured as second order sections.  Changes                    
C    the cutoffs of normalized lowpass and highpass filters through             
C    a simple polynomial transformation.                                        
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    F                       New cutoff frequency                               
C                                                                               
C  Input/Output Arguments:                                                      
C  -----------------------                                                      
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C    NSECTS                  Number of second order sectionsects                
C                                                                               
C                                                                               
      SUBROUTINE CUTOFFS( SN, SD, NSECTS, F )                           
C                                                                               
        REAL*4 SN(1), SD(1)                                             
C                                                                               
        SCALE = 2.*3.14159265*F                                         
C                                                                               
        IPTR = 1                                                        
        DO    1 I = 1, NSECTS                                           
C                                                                               
          SN( IPTR + 1 ) = SN( IPTR + 1 ) / SCALE                       
          SN( IPTR + 2 ) = SN( IPTR + 2 ) / (SCALE*SCALE)               
          SD( IPTR + 1 ) = SD( IPTR + 1 ) / SCALE                       
          SD( IPTR + 2 ) = SD( IPTR + 2 ) / (SCALE*SCALE)               
          IPTR = IPTR + 3                                               
C                                                                               
    1   CONTINUE                                                        
C                                                                               
      RETURN                                                            
      END                                                               

C                                                                               
C  Transforms an analog filter to a digital filter via the bilinear transformati
C    Assumes both are stored as second order sections.  The transformation is   
C    done in-place.                                                             
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    SN                   Array containing numerator polynomial coefficients for
C                           second order sections.  Packed head-to-tail.        
C                                                                               
C    SD                   Array containing denominator polynomial coefficients f
C                           second order sections.  Packed head-to-tail.        
C                                                                               
C    NSECTS               Number of second order sections.                      
C                                                                               
C                                                                               
      SUBROUTINE BILIN2( SN, SD, NSECTS )                               
C                                                                               
        REAL*4 SN(1), SD(1)                                             
C                                                                               
        IPTR = 1                                                        
        DO    1 I = 1, NSECTS                                           
C                                                                               
          A0 = SD(IPTR)                                                 
          A1 = SD(IPTR+1)                                               
          A2 = SD(IPTR+2)                                               
C                                                                               
          SCALE = A2 + A1 + A0                                          
          SD(IPTR)   = 1.                                               
          SD(IPTR+1) = (2.*(A0 - A2)) / SCALE                           
          SD(IPTR+2) = (A2 - A1 + A0) / SCALE                           
C                                                                               
          A0 = SN(IPTR)                                                 
          A1 = SN(IPTR+1)                                               
          A2 = SN(IPTR+2)                                               
C                                                                               
          SN(IPTR)   = (A2 + A1 + A0) / SCALE                           
          SN(IPTR+1) = (2.*(A0 - A2)) / SCALE                           
          SN(IPTR+2) = (A2 - A1 + A0) / SCALE                           
C                                                                               
          IPTR = IPTR + 3                                               
C                                                                               
    1   CONTINUE                                                        
C                                                                               
      RETURN                                                            
      END                           
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer function len_trim(string)
	character *(*) string
	character blank*1, chatmp*80
	parameter (blank=' ')

	do i=1,len(string)
	   if(string(i:i) .ne. blank) goto 11
	enddo
11      n1=i

	do i=n1,len(string)
	   if(string(i:i) .eq. blank) goto 12
	enddo
12	n2=i-1
	nl=n2-n1+1

	chatmp=string
	string=''
	if(nl.gt.0) then
	  string(1:nl)=chatmp(n1:n2)
	else
	  nl=0
	endif
	len_trim=nl

	return
	end
