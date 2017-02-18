      subroutine pseudo (icode)
c

c*******************************************************************************
c
c     subroutine pseudo function:
c         a spectrum computing subroutine which will compute the spectra for
c         one or more components of s strong motion record at several selected
c         dampings.
c
c     history:
c         1968      - the subroutine cmpmax was coded by i.m.idriss
c         1974      - the original main program, spectacular, and all subrou-
c                     tines except cmpmax were coded by a lamont, summer 1974.
c         mar 76    - changes were made to the original main program to change
c                     the resolution option by a.lamont.  this program version
c                     was called spctlm.
c         fall 1980 - the subroutine ucmpmx was coded by r.youngs
c         feb 1982  - modifications were made by r.youngs to remove plotting
c                     option and to introduce computation of spectra for unequal
c                     time step records (subroutine ucmpmx).  also, equal time
c                     step records are interpolated to small time steps when
c                     the period is less than 10*dt and the units were changed
c                     to cm/sec**2, cm/sec and cm on the output.  this version
c                     of the program was called spctln.
c         feb 1985  - modified by roy burger, wcc pasadena, to compute pseudo
c                     velocities for the epri project and renamed pseudo.
c         apr 1986  - modified by nancy smith, wcc pasadena, to accept data in
c                     either g's or cm/sec/sec.
c         may 1986  - modified by nancy smith to be a subroutine to the main
c                     program bigp with all control information determined by
c                     the user answered queries in the main program and passed
c                     to bigp in the named common blocks com1 and com2.
c                     some variable names were changed to provide consistency
c                     with bigp and its other subroutines.
c          jun 1986 - all further modifications are documented in the main
c                     routine, bigp.f.
c
c*******************************************************************************
c
c     * * * *    limitations    * * * *
c
c     1.  the number of points in the accelerogram must be <= 270000
c     2.  the number of periods must be <= 151
c
c*******************************************************************************
c
c      * * * *    commons    * * * *
c
      character*12 cname
      character*65 title1,title2
      integer      ncomp,npts,nptuse,ndamp,nper,ntypper
      real         period(151),damp(50),dt
      logical      flagt1,flagt2
      character*1  unit
      character*15 fmt
      integer      nhead
      real         anorm
      common /com1/ ncomp,cname,npts,nptuse,tws,twe,period,ndamp,
     &    damp,dt,nper,ntypper,flagt1,flagt2,title1,title2
      common /com2/ nhead,anorm,unit,fmt
c
c*******************************************************************************
c
c      * * * *    other variables    * * * *
c
      character*3  ccomp
      character*6  head
      character*10 stanam
      character*80 dummy,head1,head2
c       /* conversion factor to convert cm/sec/sec to g's
      real*8      conv
      logical     untflg
      real        resnew(112) 
      dimension   a(270000),z(3),time(270000),b(151),head(14)
      dimension   rd(151),rv(151),prv(151),aa(151),paa(151),spmax(6)

c
c*******************************************************************************
c
c     * * * *    data initialization    * * * *
c
      parameter (conv=0.001019368)
 
c...input the wcc-pasadena standard periods
 
      data (resnew(m), m=1,112) /
     *     0.010, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017,
     *     0.018, 0.019, 0.020, 0.022, 0.024, 0.026, 0.028, 0.030,
     *     0.032, 0.034, 0.036, 0.038,
     *     0.040, 0.042, 0.044, 0.046, 0.048, 0.050, 0.055, 0.060,
     *     0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100,
     *     0.110, 0.120, 0.130, 0.140, 0.150, 0.160, 0.170, 0.180,
     *     0.190, 0.200, 0.220, 0.240, 0.260, 0.280, 0.300, 0.320,
     *     0.340, 0.360, 0.380, 0.400, 0.420, 0.440, 0.460, 0.480,
     *     0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850,
     *     0.900, 0.950, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500,
     *     1.600, 1.700, 1.800, 1.900, 2.000, 2.200, 2.400, 2.600,
     *     2.800, 3.000, 3.200, 3.400, 3.600, 3.800, 4.000, 4.200,
     *     4.400, 4.600, 4.800, 5.000, 5.500, 6.000, 6.500, 7.000,
     *     7.500, 8.000, 8.500, 9.000, 9.500, 10.00, 11.00, 12.00,
     *     13.00, 14.00, 15.00, 20.00 /
c
c*******************************************************************************
c
c                   p r o g r a m
c
c---------------------------------------------------------------------------
 
      untflg = .false.
      if (unit.eq.' ' .or. unit.eq.'G' .or. unit.eq.'g') untflg=.true.

c...set up array of periods if was not done in calling program

c	   Use only the first 112 elements in resnew:
	if (ntypper.eq.2) then
	   nper = 112
	   do i=1,nper
	      period(i) = resnew(i)
	   enddo
	endif
 
c...print all the options specified
 
   62 print 325, unit,nhead
  325 format (t5,'program spctlr version pseudo 7 apr 86' /
     * t5, 'Listing of options input to the run-- '/
     * t10,'units  = ',a1 /
     * t10,'nhead = ', i2 )

c---------------------------------------------------------------------------
c   read accelerogram and compute earthquake parameters
c---------------------------------------------------------------------------
 
c...   zero the 'a' array
         do 82 ia=1,270000
82       a(ia) = 0.0
c...   read earthquake data, determine input format,normalization fator
c...   and read and print header cards in the accelerogram file
         if (anorm .eq. 0.0) anorm = 1.0
         print 309,cname
309      format ('Component- ',a12)
         if (fmt .eq. '               ') fmt = '(8f9.6)        '

c...   bypass first standard from norma, then get kg (# pts) & dt

      if (flagt1 .and. flagt2) then
           if (fmt.ne.'binary         ' .and. fmt.ne.'mclaren        ') then
              read (3,14) stanam,title1
14            format (a10,5x,a65)
           else
              if (fmt.eq.'binary         ') then
                 read (3) stanam,ccomp,title1
              else
                 read (3) head1,head2
                 read (head1,14) stanam,title1
              endif
           end if
         title2 = stanam //
     &   '                                                       '
      else
         if (.not.flagt1 .and. flagt2) then
           if (fmt.ne.'binary         ' .and. fmt.ne.'mclaren        ') then
              read (3,140) stanam,(head(ia),ia=1,11)
140           format (a10,5x,10a6,a5)
           else
              if (fmt.eq.'binary         ') then
                 read (3) stanam,ccomp,head
              else
                 read (3) head1,head2
                 read (head1,140) stanam,(head(ia),ia=1,11)
              endif
           end if
           title2 = stanam //
     &      '                                                       '
         else
           if (.not.flagt1 .and. .not.flagt2) then
             if (fmt.ne.'binary         ' .and. fmt.ne.'mclaren        ') then
                read (3,140) stanam,(head(ia),ia=1,11)
             else
                if (fmt.eq.'binary         ') then
                   read (3) stanam,ccomp,head
                else
                   read (3) head1,head2
                   read (head1,140) stanam,(head(ia),ia=1,11)
                endif
             end if
           else
             if (fmt.ne.'binary         ' .and. fmt.ne.'mclaren        ') then
                read (3,14) stanam,title1
             else
                if (fmt.eq.'binary         ') then 
                   read (3) stanam,ccomp,title1
                else 
                   read (3) head1,head2 
                   read (head1,14) stanam,title1
                endif
             end if
           endif
         end if
      end if

      write (6,14) stanam,title1
      if (fmt.ne.'binary         ' .and. fmt.ne.'mclaren        ') then
         read (3,*) kg,dt
      else
         if (fmt.eq.'binary         ') then
            read (3) kg,dt
         else
            read (head2,'(i10,f10.7)') kg,dt
         endif
      end if
      write (6,*) kg,dt
 
c --- bypass extra headers user told you about

      keqdt=1
      if (dt .le. 0.0) keqdt=0
      if (nhead .eq. 0) go to 90
      print 17
17    format (t5, 'Header cards from accelerogram deck--')
      do 89 i11 = 1,nhead
         if (fmt.ne.'binary         ' .and. fmt.ne.'mclaren        ') then
            read (3, 8881) dummy
         else 
           if (fmt.eq.'binary         ') read (3) dummy
         end if
         write (6,8881) dummy
8881     format (a80)
   89 continue

c --- read the data

90    continue
      if (keqdt .eq. 1) then
         if (fmt.ne.'binary         ' .and. fmt.ne.'mclaren        ') then
            read (3,fmt) (a(l), l=1,kg)
         else
            read (3) (a(l),l=1,kg)
         end if
      end if
      if (keqdt .eq. 0) read(3,fmt) (time(l),a(l),l=1,kg)

      if (nptuse .eq. 0) then
         npts =  kg
         do 911 i=1,kg
911      time(i)=float(i-1)*dt
         tws = 0.0
         twe = time(kg)
      end if

      if (nptuse.gt.0) then
         if (keqdt .eq. 0) goto 92
         npts = nptuse
         do 91 i=1,nots
91       time(i)=float(i-1)*dt
         tws = 0.0
         twe = time(npts)
      end if

      if (nptuse.lt.0) then
         if (keqdt.eq.0) then
            print *,'error from subroutine pseudo (spectacular:)'
            print *,'   pseudo subroutine not set up to window data'
            print *,'   when time is read in and not computed.  exe-'
            print *,'   cution terminated.'
            icode = 1
            goto 5500
         end if
         if (tws.eq.0) then
            nstart = 1
         else
            nstart = tws/dt
         end if
         nend   = twe/dt
         if (nend .gt. kg) then
            nend = kg
            twe  = nend*dt
         end if
         npts = nend-nstart+1
         do 915  i=1,npts
            jj = i+nstart-1
            time(i) = a(jj)
915      continue
         do 916  i=1,npts
            a(i) = time(i)
            jj = i+nstart-1
            time(i) = float(jj) * dt
916      continue
      end if

c...find max velocity, displacement, and acceleration
c...using the first npts values.
 
92    amaxg = 0.0
      v1=0.
      d1=0.
      vm=0.
      dm=0.
      kgmax = npts - 1
      do 300 k=1, kgmax
         if (keqdt .eq. 0) dt=time(k+1)-time(k)
         if (untflg) then
            vt = v1 + 0.5*dt*(a(k) + a(k + 1)) * anorm
            td = d1+ dt*(v1 + dt * anorm * a(k)/3.+ dt * anorm * a(k+1)/6.)
         else
            vt = v1 + 0.5*dt*(a(k) + a(k + 1)) * (anorm*conv)
            td = d1 + dt*(v1 + dt * (anorm*conv) * a(k)/3.+
     &           dt * (anorm*conv) * a(k+1)/6.)
         end if
         if (abs(vt) .lt. vm) go to 100
         vm=abs(vt)
         tmv = dt*k
         if (keqdt .eq. 0) tmv=time(k+1)
  100    if (abs(td) .lt. dm) go to 200
         dm = abs(td)
         tmd = dt*k
         if (keqdt .eq. 0) tmd=time(k+1)
  200    v1 = vt
         d1 = td
         if (untflg) then
            q = anorm * abs(a(k))
         else
            q = (anorm*conv) * abs(a(k))
         end if
         if (q .lt. amaxg) go to 300
         amaxg = q
         tma = dt * (k-1)
            if (keqdt .eq. 0) tma=time(k+1)
  300 continue

      vm = 981. * vm
      dm = 981. * dm
      amax = amaxg * 981.0
      print 310,dm, tmd, vm, tmv, amax, tma, amaxg, tma, dt, kg, npts
310   format(/, t5, 'data for earthquake-', /,
     1 t21, '********************', t21, 'max. displacement = ',f10.5,
     2 ' cm.', t63, 'at time = ', f10.5, ' secs', /, t21,
     3 'max. velocity     = ', f10.5, ' cm/sec'   , t63, 'at time = ',
     3 f10.5, ' secs',
     4 /, t21, 'max. acceleration = ', f10.5, ' cm/sec**2', t63,
     5 'at time = ', f10.5, ' secs', /, t39, '=', f10.5, ' g''s', t63,
     6 'at time = ', f10.5, ' secs', /,  t21, 'time step = ', f10.5,
     7 ' secs',/ , t21, 'number of points in the accelerogram = ', i5
     8 /, t21, 'number of points used in response computation = ',i5)
      if (keqdt .eq.0) print 332
332   format (t5,'** unequal dt record **')
      if (anorm .eq. 1.0) print 327
327   format (t21, 'accelerogram was not normalized.')
      if (anorm .ne. 1.0) print 328,anorm
328   format (t21, 'accelerogram was normalized by a ',
     *       'factor of ', f6.3, /, t21, 'the above values were ',
     *       'computed from the normalized accelerogram.')
 
c---------------------------------------------------------------------------
c   make checks and print and normalize accelerogram
c---------------------------------------------------------------------------
 
c...check to see that the minimum period is not less than 2 * delta t
 
      tdelt = 2.0 * dt
      if (tdelt.gt.period(1) .or. tdelt.gt.period(nper ))
     &   print 311,tdelt, dt
311   format(t21, 'warning **** spectral values computed at ',
     *      'periods less than ', f6.4, ' secs. are not reliable', /,
     *      t34, 'since the time step is ', f6.4, ' secs.')
 
c... convert acc's to cm/sec**2, and normalize
 
      do 320 i4 = 1,npts
      if (untflg) then
         a(i4) = a(i4) * 981.0 * anorm
      else
         a(i4) = a(i4) * anorm
      end if
  320 continue
 
c...set kg to npts, if necessary
 
      if (npts .gt. 0) kg = npts
 
c----------------------------------------------------------------------------
c   compute response
c---------------------------------------------------------------------------
 
c...for each damping do the following:
 
      do 800 i=1,ndamp
         kug = kg - 1
         dur1 = 0.0
         d = damp(i)
         yy = sqrt(1.-d*d)
         do 600 n= 1, nper 
            w = 4.*asin(1.0)/period(n)
            wd = yy*w
            w2 = w*w
            w3 = w2*w
c
c...compute response
c
       temp = 10*dt 
       print *,keqdt,period(n),temp
            if(keqdt.eq.0 .or. period(n).lt. 10.*dt) go to 350
            call cmpmax (dur1,kug,a,period(n),w,w2,w3,wd,d,dt,z)
            go to 360
350         call ucmpmx(dur1,kug,a,time,period(n),w,w2,w3,wd,d,z)
360         rd(n) = z(1)
            rv(n) = z(2)
            aa(n) = z(3)/981.0
            prv(n)= w*z(1)
            paa(n) = w2*z(1)/ 981.0
600      continue
 
c-------------------------------------------------------------------------
c   print values
c---------------------------------------------------------------------------
 
         do 699 iii=1,6
699      spmax(iii) = 0.

         print 312, d
312      format ( t7,'Damp =',f5.3,//, t45, 15hSPECTRAL VALUES//
     1    ' No    Period     Freq   Rel Disp      Rel Vel     Psu Rel Vel   Abs Accel   Psu Abs Acc    Ratio'/
     2    '        (sec)     (hz)     (cm)        (cm/sec)     (cm/sec)       (g''s)        (g''s)' )
 
         do 700 n=1, nper 
            b(n)=aa(n)/amaxg
            f = 1/period(n)
            if (rd(n) .gt. spmax(1)) spmax(1) = rd(n)
            if (rv(n) .gt. spmax(2)) spmax(2) = rv(n)
            if (prv(n) .gt. spmax(3)) spmax(3) = prv(n)
            if (aa(n) .gt. spmax(4)) spmax(4) = aa(n)
            if (paa(n) .gt. spmax(5)) spmax(5) = paa(n)
            if (b(n) .gt. spmax(6)) spmax(6) = b(n)
            print 322,n,period(n),f, rd(n), rv(n), prv(n), aa(n), paa(n), b(n)
322         format (i3,f10.3,f10.2,6e13.5)
700      continue
 
c       print maximum values
         print 308,(spmax(iii),iii=1,6)
308      format (102('-'),/,'Maximum Values --->    ',6e13.5/)
         idmp=ifix(d*1000.+0.01)
         idmp1 = idmp/100
         idmp2 = (idmp-idmp1*100)/10
         idmp3 = idmp-idmp1*100-idmp2*10
         write (4,301) cname,title1,title2,nper,d
301      format(a12/a65/a65/i4,f6.3,
     *     ' Period  Rel Disp      Rel Vel     Pseudo RV   AbsoluteAcc  Pseudo AA       Ratio' )
         write (4,302) (i9,period(i9),rd(i9),rv(i9),prv(i9),aa(i9),paa(i9),
     *     b(i9), i9 = 1,nper )
302      format (i3, 7e13.5)
 
800   continue
 
5500  continue
      close(3)
      close(4)
      return
      end

      subroutine ucmpmx(dur1,kug,ug,time,pr,w,w2,w3,wd,d,z)
      dimension ug(270000),time(270000),z(3),t(3),c(3),x(2,3)

	print *,dur1,kug,ug(1),time(1),time(kug),pr
      do 10 i=1,3
      x(1,i)=0.
   10 z(i)=0.

      f2=1./w2
      f3=d*w
      f4=1./wd
      f5=f3*f4
      f6=2.*f3

      do 100 k=1,kug
         dt=time(k+1)-time(k)
         ns=ifix(10.*dt/pr-0.01)+1
         dt=dt/float(ns)
   20    f1=2.*d/w3/dt
         e=exp(-f3*dt)
         g1=e*sin(wd*dt)
         g2=e*cos(wd*dt)
         h1=wd*g2-f3*g1
         h2=wd*g1+f3*g2
         dug=(ug(k+1)-ug(k))/float(ns)
         g=ug(k)
         z1=f2*dug
         z3=f1*dug
         z4=z1/dt
         do 100 is=1,ns
            z2=f2*g
            b=x(1,1)+z2-z3
            a=f4*x(1,2)+f5*b+f4*z4
            x(2,1)=a*g1+b*g2+z3-z2-z1
            x(2,2)=a*h1-b*h2-z4
            x(2,3)=-f6*x(2,2)-w2*x(2,1)
            do 80 l=1,3
               c(l)=abs(x(2,l))
               if (c(l) .lt. z(l)) go to 80
               z(l)=c(l)
               t(l)=time(k)+is*dt+dur1
   80	        x(1,l)=x(2,l)
            g=g+dug
  100 continue
c     print 112,pr,(t(l),l=1,3)
c 112 format (' ucmpmx per =',f6.3,5x,'times for maxima --',3x,
c    ,       'td =',f8.4,'   tv =',f8.4,'   ta =',f8.4)
      return
      end
      subroutine cmpmax (dur1,kug,ug,pr,w,w2,w3,wd,d,dt,z)
      dimension ug(270000),x(2,3),t(3),z(3),c(3)
c
      do 10 i=1,3
      x(1,i)=0.
   10 z(i)=0.
      f1 = 2.*d/(w3*dt)
      f2 = 1./w2
      f3 = d*w
      f4 = 1./wd
      f5 = f3*f4
      f6 = 2.*f3
       e = exp(-f3*dt)
      g1 = e*sin(wd*dt)
      g2 = e*cos(wd*dt)
      h1 = wd*g2 - f3*g1
      h2 = wd*g1 + f3*g2
      do 100 k = 1, kug
      dug = ug(k+1) - ug(k)
      z1 = f2*dug
      z2 = f2*ug(k)
      z3 = f1*dug
      z4 = z1/dt
       b = x(1,1) + z2 -z3
       a = f4*x(1,2) + f5*b + f4*z4
      x(2,1) = a*g1 + b*g2 + z3 - z2 - z1
      x(2,2) = a*h1 - b*h2 - z4
      x(2,3)=-f6*x(2,2)-w2*x(2,1)
      do 80 l=1,3
      c(l)=abs(x(2,l))
      if (c(l) .lt. z(l)) go to 80
      z(l)=c(l)
      t(l)=dt*float(k)+dur1
   80 x(1,l)=x(2,l)
  100 continue
c     print 112, pr, (t(l),l=1,3)
c 112 format(' cmpmax per =',f6.3,5x,19htimes for maxima -- ,3x,
c    1 4htd = ,f8.4,3x,4htv = ,f8.4,3x,4hta = ,f8.4)
      return
      end
