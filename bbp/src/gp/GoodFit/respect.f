c                              respect.f - revision 1.0
c*******************************************************************************
c
c     * * * *    history    * * * *
c
c   6/86 - nfs - added capability to read data in binary format by
c                entering the word "binary" or "binary" when asked for
c                the data format.  subroutine pseudo checks for this
c                in the variable format and reads accordingly.
c              - added check so that if ncomp .ne. 3, averaging will
c                not be done.
c   7/86 - nfs - modified so second title line defaults to station
c                name and component from first data file header line.
c              - added capability of windowing data
c   3/87 - nfs - revision 2
c              - changed all routines so that averaging in subroutine
c                averag is done using the absolute acceleration rather than
c                the pseudo absolute acceleration
c              - corrected test in main routine after question re
c                scaling by the average of psaa so that program did not
c                think there were multiple dampings when there weren't.
c   9/87 - nfs - revision 3
c              - modified bigp and spcplt so that multiple files can be
c                plotted on the same set of axes.  program now asks additional
c                questions.  see documentation for detail.
c   2/95 - nfs - Major overhaul and streamlining.- renamed from bigp to respect
c
c*******************************************************************************
c
c     * * * *    limitations    * * * *
c
c     1.  the number of points in the accelerogram must be <= 120000
c     2.  the number of periods must be <= 151
c     3.  the number of components must be = 1
c     4.  the number of dampings applied and/or plotted must be <=5
c     5.  the number of x,y coordinates for the target shape must be <=25
c

c*******************************************************************************
c
c     * * * *    Declarations   * * * *

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

      character*1  ysno
      character*256 filnam
      logical      exlog

c*******************************************************************************
c
c     * * * *    processing    * * * *

c --- Determine titles for the run to be used as plot labels

101   print *,'Do you want the main title of this run taken from the'
      print *,'first line of the input data file?  If yes, the title is'
      print *,'assumed to be in columns 16-80 of line 1 of the data file.'
      read (5,825) ysno
825   format (a1)
      write (6,825) ysno
      if (ysno.eq.'y' .or. ysno.eq.'Y') then
         flagt1 = .true.
      else
         flagt1 = .false.
         print *,'Main title for this run:'
         read (5,800) title1
         write (6,800) title1
      end if
800   format (a65)

      print *,'Do you want the 2nd title line taken from line 1 of the'
      print *,'file header (station name)?'
      read (5,825) ysno
      write (6,825) ysno
      if (ysno.eq.'y' .or. ysno.eq.'Y') then
         flagt2 = .true.
      else
         flagt2 = .false.
         print *,'Secondary title for this run:'
         read (5,800) title2
         write (6,800) title2
      end if

c --- Determine the components to be considered

      print *,'Enter name of the component:'
      read (5,810) cname
      write (6,810) cname
810   format (a12)

c --- Input info particular to spectacular (pseudo)

112   print *,'Input data file name:'
      read (5,820) filnam
      write (6,820) filnam
820   format (a256)
      inquire (file=filnam,exist=exlog)
      if (.not.exlog) then
         print *,'***** INPUT DATA FILE DOES NOT EXIST *****'
         stop '***** EXECUTION TERMINATED *****'
      end if

      print *,'Format for the input data.'
      print *,'Enter "binary" if data are in wcc binary format and "mclaren"'
      print *,'if data are in Jim McLaren''s binary format:'
      read (5,815) fmt
815   format (a15)
      if (fmt.ne.'binary         '.and.fmt.ne.'mclaren        ') then
         open (3,file=filnam)
      else
         open (3,file=filnam,form='unformatted')
      end if

      tws = 0.
      twe = 0.
      print *,'Data points to use:'
      print *,'   0 = all'
      print *,'   >1= the number of points to use'
      print *,'   -1= window the data'
      read (5,*) nptuse
      write (6,*) nptuse
      if (nptuse.lt.0) then
         print *, 'start and end times for the window:'
         read (5,*) tws,twe
         write (6,*) tws,twe
      end if

      print *,'Enter response spectral ouput file name:'
      read (5,820) filnam
      write (6,820) filnam
      open (4,file=filnam)
      close (4,status='delete')
      open (4,file=filnam)

114   print *,'what periods do you want to consider in spectacular?'
      print *,'     1 = unused'
      print *,'     2 = wcc-pasadena 106 standard periods'
      print *,'     3 = user inputted periods'
      print *,'     4 = user inputted frequencies'
      read (5,*) ntypper
      write (6,*) ntypper
      if (ntypper.lt.1 .or. ntypper.gt.4) stop 'Number of periods entered is invalid.'
      if (ntypper.eq.1) stop 991
      if (ntypper.eq.2) nper=106
      if (ntypper.eq.3) then
         print *,'How many periods?'
         read (5,*) nper
         write (6,*) nper
         print *,'Enter the periods desired:'
         read (5,*) (period(i),i=1,nper)
         write (6,*) (period(i),i=1,nper)
      end if
      if (ntypper.eq.4) then
         print *,'How many frequencies?'
         read (5,*) nper
         write (6,*) nper
         print *,'Enter the frequencies desired:'
         read (5,*) (period(i),i=1,nper)
         write (6,*) (period(i),i=1,nper)
         do i=1,nper
            period(i) = 1/period(i)
         enddo
         write (6,*) (period(i),i=1,nper)
      end if

      print *,'In what units are the data?'
      print *,'     1 = cm/sec/sec'
      print *,'     2 = g''s'
      read (5,*) nans
      write(6,*) nans
      if (nans.eq.2) then
         unit = 'g'
      else
         unit = 'c'
      end if

      print *,'Number of header cards to skip before the data begin:'
      print *,'It is assumed that the first card is a title and the'
      print *,'second contains the number of points and delta t.'
      print *,'These cards will be automatically read and need not be'
      print *,'skipped.'
      read (5,*) nhead
      write(6,*) nhead

      print *,'How many damping values are to be applied to the data?'
      read (5,*) ndamp
      write(6,*) ndamp
      print *,'Enter the damping values (max=5):'
      read (5,*) (damp(i),i=1,ndamp)
      write(6,*) (damp(i),i=1,ndamp)

      print *,'Enter the normalization factor (0. is converted to 1.)'
      read (5,*) anorm
      write(6,*) anorm

c --- start working

      icode = 0
      call pseudo (icode)
      if (icode.ne.0) then
         close (11)
         close (12)
         stop '***** ERROR IN SPECTACULAR (PSEUDO); EXECUTION TERMINATED *****'
      end if

      print *,'Do you want to process another set of data?'
      read (5,825) ysno
      write(6,825) ysno
      if (ysno.eq.'y' .or. ysno.eq.'Y') goto 101

      stop
      end
