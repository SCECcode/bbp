c     program model_uncer_plt
c     Reads output from rspresid and creates input files for mplot and
c     files ready for plotting with mplot and files containing
c        1)  the residuals of the validations,
c        2)  the mean residual (model bias),
c        3)  the standard error of a single observation (random and modelling
c            uncertainty).

      character*1  spec_type
      character*80 file(6)
      integer ieq(90000),nresid(90000)
      real per(90000),sigma(900),sigmaZero(900),misfit(90000),bias(90000), 
     1     sigmaBias(90000), CL90(2,900), freq(90000)
      character*80 file1
	include 'inserts'

c  First period will be made max accelerations; 86 periods are in the
c  input data files.
      print *,'Enter number of periods in the input data files:'
      read (*,*) nper
      nper = nper + 1

      print *,'Do you want mplot input files output with periods or frequencies?'
      print *,'Enter "p" or "f"'
      read (*,'(a1)') spec_type
      write (*,'(a1)') spec_type
      specflag=0
      if (spec_type.eq.'f') specflag=1

c  Open mplot input files
      print *,' Enter 6 file names, bias, bias-,+ conf interval,sigma, sigma0,'
      print *,'and residuals for mplot:'
      read (*,'(a80)') (file(i),i=1,6)
      write (*,'(a80)') (file(i),i=1,6)
      open (31,file=file(1))
      open (32,file=file(2))
      open (33,file=file(3))
      open (34,file=file(4))
      open (35,file=file(5))
      open (36,file=file(6))

c     Open residual data files
      write (*,'( 2x,''Enter name of file with residuals that was written by rspresid:'')')
      read (*,'( a60)') file1
      write (*,'( a60)') file1
      open (10,file=file1,status='old')

c     Read the residual data file
      call rdresid ( nsta, nper, npts, per, freq, misfit, ieq )

c     Prepare mplot input files of the residuals
      call biasplt ( misfit, bias, npts, per, freq, sigmabias, ieq, 1,nresid)

c     Compute the mean residual
      call meanResid (nsta,nper,misfit,bias,sigma,sigmaZero,sigmaBias,nresid)

c     Save the modelling plus random standard error
      call writeSigma ( nPer, per, freq, sigma, sigmaZero, bias )

c     Compute the 90% confidence interval of the mean using the T distribution
      call confInter ( CL90, bias, sigma, nsta, nper,nresid )
      
c     Prepare mplot input files of the model bias
      call biasplt ( misfit, bias, nper, per, freq, CL90, ieq, 0,nresid)

c     Plot the se of a single observation
      call varplt ( nper, per, freq, sigma, sigmaZero ,nresid)

      stop
      end

c ----------------------------------------------------------------

      subroutine writeSigma ( nPer, per, freq, sigma, sigmaZero, bias )
      integer nPer
      real per(1), sigma(1), sigmaZero(1), bias(1), freq(1)
      character*80 fileout
      include 'inserts'

      iunit = 15
      write (*,'( 2x,''Enter model_sigma ouput file name'')')
      read (*,'( a60)') fileout
      open (iunit,file=fileout)
      
c     Write the header
      write (iunit,'(''Component =         '',i1)') icomp
      write (iunit,'(i1,''    =   Number of events'')') nEvent
      write (iunit,'(''Events and number of stations:'')')
      do i=1,nEvent
         write (iunit,'(a20,i5)') evName(i),nSim(i)
      enddo
      write (iunit,'( 2x,''Modelling plus random standard error'')')
      write (iunit,'( 2x,''period, sigma, sigmaZero (mean=0), bias'')')
      write (iunit,'( i5,''  nperiod'')') nPer

c     Write the sigmas
      do i=1,nPer
        write (iunit,'( 4f10.4)') per(i), sigma(i), sigmaZero(i), bias(i)
      enddo
      close (iunit)
      return
      end

c ----------------------------------------------------------------

      subroutine rdresid (nsta,nper,npts,per,freq,misfit,ieq )
      real per(1), misfit(1), freq(1)
      integer ieq(1), npts, nper, nsta
      character*80 dummy
      include 'inserts'

c     Read header information from residual file (output from rspresid.f)
	iunit = 10
       read (iunit,'(20x,i1)') icomp
       read (iunit,*) nEvent
       read (iunit,'(a80)') dummy
       do i=1,nEvent
          read (iunit,'(a20,i5)') evName(i),nSim(i)
       enddo

c     Set up eqk number array for each period and station
      nsta = 0
      jj = 1
      do i=1,nEvent
        do j=1,nSim(i)
          do k=1,nper
            ieq(jj) = i
            jj = jj + 1
          enddo
        enddo
        nsta = nsta + nSim(i)
      enddo
      write (*,'( 2x,''nsta ='',i5)') nsta

c     Read residuals
      n1 = nsta * nper
      do i = 1,n1 
        read (iunit,*,end=10) ista, per(i), misfit(i)
        freq(i) = 1./per(i)
        npts = i
      enddo
  10  close (iunit)
      write (*,'( 2x,''npts ='',i5)') npts
      return
      end

c -------------------------------------------------------------

      subroutine meanResid ( nsta, nper, misfit, bias, sigma, sigmaZero, 
     1                       sigmaMean,nresid )
      real misfit(1), bias(1), sigma(1), sigmaZero(1), sigmaMean(1)
      real sum(90000), sum2(90000)
      integer nresid(90000)

      do j=1,nper
        sum(j) = 0.0
        sum2(j) = 0.0
	nresid(j) = 0
      enddo
      jj = 1
      do i=1,nsta
        do j=1,nper

	  if(misfit(jj).gt.-999.0) then
             sum(j)  = sum(j) + misfit(jj)
             sum2(j) = sum2(j) + misfit(jj)**2
	     nresid(j) = nresid(j) + 1
	  endif

          jj = jj + 1

        enddo
      enddo

	print *,'nper,bias,sigma,sigmaZero,sigmaMean:'
      do j=1,nper

	if(nresid(j).gt.1) then
           fn = float ( nresid(j) )
           bias(j) = sum(j) / nresid(j)
           sigma(j) = sqrt( ( sum2(j) - nresid(j)*bias(j)**2 ) / ( nresid(j)-1 ) ) 
           sigmaZero(j) = sqrt( sum2(j) / nresid(j) )
           sigmaMean(j) = sigma(j) / sqrt(fn)
	endif

      enddo
      return
      end

c ---------------------------------------------------------

      subroutine confInter ( CL90, bias, sigma, n, nper,nresid )

      integer n, nper, nresid(90000)
      real bias(1), sigma(1), CL90(2,1), t95(55)
      data t95 / 6.3138, 2.9200, 2.3534, 2.1318, 2.0150, 1.9432, 1.8946,
     1           1.8595, 1.8331, 1.8125, 1.7959, 1.7823, 1.7709, 1.7613,
     2           1.7531, 1.7459, 1.7396, 1.7341, 1.7291, 1.7247, 1.7207,
     3           1.7171, 1.7139, 1.7109, 1.7081, 1.7056, 1.7033, 1.7011,
     4           1.6991, 1.6973, 1.6955, 1.6939, 1.6924, 1.6909, 1.6896,
     5           1.6883, 1.6871, 1.6860, 1.6849, 1.6839, 1.6829, 1.6820,
     6           1.6811, 1.6802, 1.6794, 1.6787, 1.6779, 1.6772, 1.6766,
     7           1.6759, 1.6753, 1.6747, 1.6741, 1.6736, 1.6730/

      do i=1,nper

	 if ( nresid(i) .gt. 1 ) then
	    if ( nresid(i) .gt. 56 ) then
	      tt95 = 1.64
	    else
	      tt95 = t95(nresid(i)-1)
	    endif

            fn = float (nresid(i))
            CL90(1,i) = bias(i) - sigma(i) * tt95 / sqrt(fn) 
            CL90(2,i) = bias(i) + sigma(i) * tt95 / sqrt(fn) 
	 endif

      enddo
      return
      end
        
c ---------------------------------------------------------

      subroutine biasplt ( misfit, bias, npts, per, freq, CL90, ieq,flag,nresid)
      real misfit(1), per(1), bias(1), freq(1)
      real CL90(2,1)
      integer flag, ieq(1),nresid(90000)
      character*80 sym(5)
	include 'inserts'
      data sym /'1','3','2','4','5'/

c     Plot Misfit (skip periods > 2 seconds)
      write (*,'( 2x,''flag ='',i5)') flag
      if (flag .eq. 0 ) goto 55
      write (*,'( 2x,''npts ='',i5)') npts
      do 50 i=1,npts
        if (misfit(i).gt.-999.0) then
           if (specflag) then
              write (36,*) freq(i),misfit(i)
           else
              write (36,*) per(i),misfit(i)
           endif
         endif
  50  continue
      return

c     Plot the bias and 90% C.I. of the bias 
  55  continue
      write (*,'( 2x,''npts ='',i5)') npts
      do 60 i=1,npts
        if (nresid(i).gt.1) then
           if (specflag) then
              write (31,*) freq(i),bias(i)
           else
              write (31,*) per(i),bias(i)
           endif
         endif
  60  continue
      do i=1,npts
        if (nresid(i).gt.1) then
            if (specflag) then
               write (32,*) freq(i),CL90(1,i)
               write (33,*) freq(i),CL90(2,i)
            else
               write (32,*) per(i),CL90(1,i)
               write (33,*) per(i),CL90(2,i)
            endif
         endif
      enddo

      return
      end

c --------------------------------------------------------------------

      subroutine varplt ( nper, per, freq, sigma, sigmaZero,nresid )
      real per(1), freq(1)
      real sigma(1), sigmaZero(1)
      integer nresid(90000)
	include 'inserts'

c     Plot the sigma of a single observation 
      do 60 m=1,nper
        if (nresid(m).gt.1) then
           if (specflag) then
              write (34,*) freq(m),sigma(m)
           else
              write (34,*) per(m),sigma(m)
           endif
        endif 
  60  continue

c     Plot the sigma of a single observation (under hypthesis that mean=0)
      do 70 m=1,nper
        if (nresid(m).gt.1) then
           if (specflag) then 
              write (35,*) freq(m),sigmaZero(m)
           else
              write (35,*) per(m),sigmaZero(m)
           endif 
        endif 
  70  continue

      return
      end
