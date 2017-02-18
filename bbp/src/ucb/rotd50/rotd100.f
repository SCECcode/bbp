!     ------------------------------------------------------------------
!
!      rotD100.f
!      Computes RotD50 and RotD100 using pairs of orthogonal horizontal components
!      Code developed by Norm Abrahamson, optimized in September 2012 to run faster
!      QA/QC and some clean-up for distribution by Christine Goulet
!      Modified in 9/14 by Scott Callaghan to calculate RotD100 also
!     ------------------------------------------------------------------
!      Input provided in filein:
!             -  Interp: mode of interpolation for small dt 
!                 1. linear interpolation
!                 2. sine wave interpolation (e.g. frequency domain interpolation)
!                 3. cubic spline interpolation 
!             - NPairs: number of pairs to read
!             - NHead: number of header lines in ASCII time series files
!             - Series of file names:
!                  file name component 1 (input)
!                  file name component 2 (input)
!                  file name RotD50 output
!                  file name as-recorded PSa output for both components
!            ** Make sure the input ASCII format is adequate - see how it is read below**
!     ------------------------------------------------------------------

      program Calc_RotD50
      parameter (MAXPTS = 2000000)

      character*80 fileacc1, fileacc2, filein, fileout_rd100
      integer npts1, npts2, npts, npair, nhead, iFlag
      real dt1, dt2, dt, acc1(MAXPTS), acc2(MAXPTS), acc1_0(MAXPTS), acc2_0(MAXPTS)
      real dt10
      real x0(MAXPTS), y0(MAXPTS), u(MAXPTS), y2(MAXPTS)
      real rspTH1(MAXPTS), rspTH2(MAXPTS), rsp1(MAXPTS), rsp2(MAXPTS)
      real x(MAXPTS), y(MAXPTS), rsp_Period(63), famp15(3)
      real rotangle, w(200), sa(1000), workArray(1000), rotD50(3,200), rotD100(3,200), psa5E(200), psa5N(200)
      integer rD100ang(3,200), rD50ang(3,200)
      real saUnsort(1000)
      real damping
      complex cu1(MAXPTS)

      data RSP_Period / 0.010, 0.011, 0.012, 0.013, 0.015, 0.017, 0.020, 0.022, 0.025, 0.029, 
     1               0.032, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.075, 0.085, 
     2               0.100, 0.110, 0.120, 0.130, 0.150, 0.170, 0.200, 0.220, 0.240, 0.260, 
     3               0.280, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.750, 
     4               0.850, 1.000, 1.100, 1.200, 1.300, 1.500, 1.700, 2.000, 2.200, 2.400, 
     5               2.600, 2.800, 3.000, 3.500, 4.000, 4.400, 5.000, 5.500, 6.000, 6.500, 
     6               7.500, 8.500, 10.000 /	

      nFreq = 63
      damping = 0.05
      dt_max = 0.001

!     Convert periods to freq in Radians
      do iFreq=1,nFreq
        w(iFreq) = 2.0*3.14159 / rsp_period(iFreq)
      enddo

!     Read in the input filename filein 
      filein = "rotd100_inp.cfg"
!     Uncomment below to make interactive instead
!     write (*,*) 'Enter the input filename.'
!     read (*,*) filein
      open (30,file=filein, status='old')


!      write (*,'( 2x,''write out interpolations? (0=no, 1=yes)'')')
!     read (*,*) iFlag
      iFlag = 0
      if ( iFlag .eq. 1 ) then
        write (*,'( 2x,'' enter times (sec) of first and last points to write out'')')
        read (*,*) time1, time2
      endif

!     Read interpolation method to be used
      read (30,*) jInterp
!     Read number of pairs and number of header lines
      read (30,*) nPair
      read (30,*) nHead

!     Loop over each pair
      do iPair=1,nPair

!    Uncomment below for screen output
!        write (*,'( 2x,'' set '',i5)') ipair

!       Read Horiz 1 (x) component
        read (30,'(a80)') fileacc1
!    Uncomment below for screen output
!        write (*,'( a70)') fileacc1
        open (32,file=fileacc1,status='old')
        do i=1,nhead-1
          read (32,*)
        enddo
        read (32,*) npts1, dt1
        if (npts1 .gt. MAXPTS) then
          write (*,'( 2x,''NPTS is too large: '',i10)') npts
          stop 99
        endif
        read (32,*) (acc1_0(i),i=1,npts1)
        close (32)

!       Read Horiz 2 (y) component
        read (30,'(a80)') fileacc2
!    Uncomment below for screen output
!        write (*,'( a70)') fileacc2
        open (33,file=fileacc2,status='old')
        do i=1,nhead-1
          read (33,*)
        enddo
        read (33,*) npts2, dt2
        if (npts2 .gt. MAXPTS) then
          write (*,'( 2x,''NPTS is too large: '',i10)') npts
          stop 99
        endif
        read (33,*) (acc2_0(i),i=1,npts2)
        close (33)

!       Check that the two time series have the same number of points.  If not, reset to smaller value
        if (npts1 .lt. npts2) then
          npts0 = npts1
        elseif (npts2 .lt. npts1) then
          npts0 = npts2
        elseif (npts1 .eq. npts2) then
          npts0 = npts1
        endif
      
!       Check that the two time series have the same dt.
        if (dt1 .ne. dt2) then
          write (*,*) 'DT values are not equal!!!'
          write (*,*) 'DT1 = ', dt1
          write (*,*) 'DT2 = ', dt2
        else
          dt = dt1
        endif
        dt0 = dt 

!        Copy to new array for interpolating ( from original acc )
         do i=1,npts0
           acc1(i) = acc1_0(i)
           acc2(i) = acc2_0(i)
         enddo    
         npts = npts0
         dt = dt0

!        Interpolate to finer time step for calculating the Spectral acceleration
          if ( jInterp .ne. 0 ) then
           NN = 2**(int(alog(dt/dt_max)/alog(2.))+1)
            if ( NN*npts .gt. MAXPTS ) then
             write (*,'( 2x,''increase maxpts to '',i10)') nn*npts 
             read (30,'( a80)') fileout
             goto 100
            endif

!    Uncomment below for screen output
!            write (*,'( 2x,''jInterp, dt, interpolation factor '',i5,f10.5,i5)') jinterp, dt, NN

!         Time domain linear interpolation
          if ( jINterp .eq. 1 ) then
            call InterpTime (acc1, dt, npts, dt10, npts10, NN, MAXPTS )
            call InterpTime (acc2, dt, npts, dt10, npts10, NN, MAXPTS )

!         Freq domain interpolation (sine wave)
          elseif (jInterp .eq. 2 ) then
            call InterpFreq (acc1, npts, cu1, NN, npts10, dt )
            call InterpFreq (acc2, npts, cu1, NN, npts10, dt )
            dt10 = dt / NN

!         Time domain cubic spline interpolation
          elseif (jInterp .eq. 3 ) then
            call InterpSpline (acc1, dt, npts, dt10, npts10, NN, y2, x0, y0, u, MAXPTS )
            call InterpSpline (acc2, dt, npts, dt10, npts10, NN, y2, x0, y0, u, MAXPTS )
            dt10 = dt / NN
!            write (*,'(2x,''out of spline'')')
!             write (*,'( 5e15.6)') (acc1(k) ,k=33600,33620)
!           pause

          endif
          npts = npts10    
          dt = dt10
         endif
!         write (*,'( i8, f10.6,2x,''  npts, dt'')') npts, dt
!        if ( iFlag .eq. 1 ) then
!          k1 = int(time1 / dt)
!           k2 = int (time2 /dt ) + 1
!           write (10+jInterp,'( f12.6,e15.6 )') (dt*(k-1), acc1(k),k=33600,34600)
!         endif


!        Loop over each oscilator frequency
         do iFreq=1,nFreq 
!    Uncomment below for screen output
!          write (*,'( i5, f10.3)') iFreq, rsp_period(iFreq)

!         Compute the oscillator time histoires for the two components. 
          call CalcRspTH ( acc1, npts, dt, w(iFreq), damping, rspTH1 )
          call CalcRspTH ( acc2, npts, dt, w(iFreq), damping, rspTH2 )

!         Fill new array with points with amplitude on one component at least SaMin/1.5
!         This sets the points for the rotation to speed up the calculation
          call Calc_Sa ( rspTH1, sa1, npts)
          call Calc_Sa ( rspTH2, sa2, npts )
          test = amin1(sa1, sa2) / 1.5
          j = 1
          do i=1,npts
            amp1 = abs(rspTH1(i))
            amp2 = abs(rspTH2(i))
            if ( amp2 .gt. amp1 ) amp1 = amp2
            if ( amp1 .gt. test .and. iFlag .eq. 0 ) then
              rsp1(j) = rspTH1(i)
              rsp2(j) = rspTH2(i)
              j = j + 1
            endif 
          enddo
          npts1 = j -1
                 
!         Loop over different rotation angles and compute response spectra by rotating the Oscillator TH
          do j=1,90
            rotangle = real(((j-1)*3.14159)/180.0)
            cos1 = cos(rotangle)
            sin1 = sin(rotangle)
            do i=1,npts1
              x(i)=cos1*rsp1(i) - sin1*rsp2(i)
              y(i)=sin1*rsp1(i) + cos1*rsp2(i)
            enddo

!          Find the maximum response for X and Y and load into a single Sa array
            call Calc_Sa ( x, saX, npts1 )
            call Calc_Sa ( y, saY, npts1 )
            sa(j) = saX
            sa(j+90) = SaY
            saUnsort(j) = sa(j)
            saUnsort(j+90) = sa(j+90)
          enddo      

!         Get the as-recorded PSa
          psa5E(iFreq) = sa(1)
          psa5N(iFreq) = sa(91)

!         Sort the Sa array to find rotD100
          n1 = 180
          call SORT(Sa,WorkArray,N1)
          rotD100(jInterp,iFreq) = Sa(180)
          rotD50(jInterp,iFreq) = ( Sa(90) + Sa(91) ) /2.

!         Find the corresponding angle
          do i=1,180
           if ( rotD100(jInterp,iFreq) .eq. saUnsort(i) ) then
            rD100ang(jInterp,iFreq) = i
           endif
           if ( rotD50(jInterp,iFreq) .eq. saUnsort(i) ) then
            rD50ang(jInterp,iFreq) = i
           endif
          enddo
         enddo

c        Find the Famp1.5 (assumes order of freq are high to low)
         do iFreq=2,nFreq 
           shape1 = rotD100(jInterp,iFreq)/rotD100(jInterp,1)

           if ( shape1 .ge. 1.5 ) then
             famp15(jInterp) = 1./rsp_period(iFreq)
             goto 105
           endif
         enddo
  105    continue

!       Open output files for writing
        read (30,'( a80)') fileout_rd100
        open (40,file=fileout_rd100,status='replace')

!       Write RotD100, RotD50 and Psa5 file
        write (40,'(''#'', 2x, ''Psa5_N'', x, ''Psa5_E'', x, ''RotD50'', x, ''RotD100'')')
        write (40,'(''#'', 2x, a80)') fileacc1
        write (40,'(''#'', 2x, a80)') fileacc2
        write (40,'(''#'', 2x, i5, f10.4)') nFreq, damping
        do iFreq=1,nFreq
           write (40,'(f10.4, 1x, e10.5, 1x, e10.5, 1x, e10.5, 1x, e10.5)') rsp_period(iFreq),psa5N(iFreq),psa5E(iFreq),rotD50(jInterp,iFreq),rotD100(jInterp,iFreq)
!          write (40,'(f10.4,1x,e10.5,1x,e10.5,1x,e10.5,1x,i3)') rsp_period(iFreq),psa5N(iFreq),psa5E(iFreq),rotD100(jInterp,iFreq),rD100ang(jInterp,iFreq)
        enddo
!        write (40,'(''#'', 2x, ''Psa5_N'', x, ''Psa5_E'', x,''RotD50'')')
!        write (40,'(''#'', 2x, a80)') fileacc1
!        write (40,'(''#'', 2x, a80)') fileacc2
!        write (40,'(''#'', 2x, i5, f10.4)') nFreq, damping
!        do iFreq=1,nFreq
!         write (40,'(f10.4,1x,e10.5,1x,e10.5,1x,e10.5,1x,i3)') rsp_period(iFreq),psa5N(iFreq),psa5E(iFreq),rotD50(jInterp,iFreq),rD50ang(jInterp,iFreq)
!        enddo
        close (40)

 100    continue
      enddo

      close (30)

!      stop
      end

! ---------------------------------------------------------------------

      subroutine Calc_Sa ( x, Sa, npts )
      real x(1), Sa

      sa = -1E30
      do i=1,npts
        x1 = abs(x(i))
        if ( x1 .gt. Sa ) Sa = x1
      enddo
      return
      end

! ---------------------------------------------------------------------
      Subroutine InterpTime (acc1, dt, npts, dt10, npts10, NN, MAXPTS )

      real acc1(1)
      real acc2(MAXPTS)

      k = 1
      do i=1,npts-1
        do j=1,NN
          dy = (acc1(i+1)-acc1(i))/NN
          acc2(k) = acc1(i) + dy*(j-1)
          k = k + 1
        enddo
      enddo
      acc2(k) = acc1(npts)
      npts10 = k
      dt10 = dt / NN

      do i=1,npts10
        acc1(i) = acc2(i)
      enddo
      return
      end

      

! ---------------------------------------------------------------------

      subroutine InterpSpline (acc1, dt, npts, dt10, npts10, NN, y2, x0, y0, u, MAXPTS )

      real acc1(MAXPTS)
      real y2(MAXPTS), x0(MAXPTS), y0(MAXPTS), u(MAXPTS)

c     Set x array
      do i=1,npts
        x0(i) = i*dt
        y0(i) = acc1(i)
        y2(i) = 0.
      enddo
      yp1 = 0.
      ypn = 0.

      call spline( x0, y0, npts, yp1, ypn, y2, u,MAXPTS)

      k = 1
      do i=1,npts-1
        do j=1,NN
          x_new = x0(i) + dt*float(j-1)/NN
          call splint(x0,y0,y2,npts,x_new,y_new)
c          if ( k .gt. 33600 .and. k .lt. 33620 ) then
c             write (*,'( 3i8, 3f10.4,3e15.6)') k, i,j, x0(i), x0(i+1), x_new, y0(i), y0(i+1), y_new
c          endif 
          acc1(k) = y_new
          k = k + 1
          if (kk .gt. MAXPTS) then
            write (*,'( 4i8)') i, j, k, maxpts 
            stop 99
          endif
        enddo
      enddo
      acc1(k) = y0(npts) 
      dt10 = dt / NN
      npts10 = k

      return
      end

