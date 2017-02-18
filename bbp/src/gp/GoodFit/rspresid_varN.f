      program rspresid
c     This program computes the residuals of the simulation.

      integer ieq(250), ista(250), iunit
      real rd, rv, prv, Sa(2,4,200), pSa
      real nSa, per(200), resid(250,200), residBar(200), phcut, plcut
      real sigma(200), sigmaMean(200)
      character*80 fileout, filein, header(5),dummy
	include 'inserts'
        
	bignegres = -1000.0

c  The pga will be determined from the sa/ratio so the number of periods
c	read in is increased by 1.
	print *,'Enter the number of periods in the data files:'
	read (*,*) nPer
	write (*,*) nPer
	nPer = nPer+1

c     Select component and source function and program version used in the 
c     simulations.
	write (*,'( 2x,''Enter comp (1-4) to be considered:'')')
	read (*,*) icomp
       write (*,*) icomp 

	print *,'Enter number of events to be considered in this run:'
	read (*,*) nEvent
	write (*,*) nEvent
	print *,'Enter number of stations in each event:'
	read (*,*) (nSim(iev), iev=1,nEvent)

c     Open file with names of response spectra files
      write (*,'( 2x,''Enter file with spectra file names'')')
      read (*,'( a80)') filein
      open (10,file=filein,status='old')

c     Open file with corner frequencies for each records pair
      write (*,'( 2x,''Enter file with corner frequencies (fmin fmax)'')')
      read (*,'( a80)') filein
      open (17,file=filein,status='old')

c     Initialize ResidBar array
      do 10 i=1,200
        residBar(i) = 0.0
  10  continue

c     Loop over number of events
      n = 1
      do 100 iev=1,nEvent

c       Loop over number of simulations per event
        do 98 iSim = 1,nSim(iev)
          write (*,'( 3i5)') iev, iSim, n

c         Open Observed Spectrum file
          read (10,'( a80)') filein
          open ( 15, file=filein,status='old')

c         Open Simulated Spectrum file
          read (10,'( a80)') filein
          open ( 16, file=filein,status='old')

          read (17,*) fmin,fmax

	  plcut = 0.0
	  if(fmax.gt.0.0) plcut = 1.0/(0.75*fmax)

	  phcut = 1.0e+15
	  if(fmin.gt.0.0) phcut = 1.0/(1.25*fmin)

c         Read Observed and simulated spectra and convert to log acc
          do 90 j=1,2
            do 70 ic=1,3
              ifile = j + 14
	        permx = -1
              read (ifile,'( a80)') (header(k),k=1,4)
              do m=2,nPer
                read (ifile,*) im, per(m), rd, rv, prv,
     1                     Sa(j,ic,m), pSa, nSa
                if (m.eq.2) then
c	           -Calculate the true pga
                  Sa(j,ic,1)=Sa(j,ic,2)/nSa
                  per(1) = 0.005
                endif
	         if (per(m).gt.permx) then
	               permx = per(m)
	               mmax = m
	         endif
              enddo
c	      -If the order of periods in the file is from largest to smallest,
c	       the data points must be renumbered so the pga is put at the end of
c	       the array, not the beginning.
	       if (mmax.eq.2) then
c	         -Save pga
	          tempSa = Sa(j,ic,1)
	          tempPer = per(1)
c	         -Move everything down 1
	          do m=1,nPer-1
	             Sa(j,ic,m) = Sa(j,ic,m+1)
	             per(m) = per(m+1)
	          enddo
c	         -Put pga back in at end of the array
	          Sa(j,ic,nPer) = tempSa
                 per(nPer) = tempPer
c	       The biggest period value must be at one end of the file or the other!
c	            This test is not correct!
c	       else if (mmax.ne.1) then
c	          stop 991
	       endif
              do 60 m=1,nPer
                if ( Sa(j,ic,m) .gt. 0. ) then
                  Sa(j,ic,m) = alog( Sa(j,ic,m) )
                else
                  write (*,'( 2x,''Zero or neg Sa in comp'',i2)') ic
                  goto 70
                endif
 60           continue
 70         continue
            close (ifile)

c           Compute average of two horizontal accelerations
            do 80 m=1,nPer
              Sa(j,4,m) = (Sa(j,2,m) + Sa(j,3,m)) / 2.
  80        continue
  90      continue

c         Compute the residuals
          do 95 m=1,nPer
            resid(n,m) = bignegres
	    if(per(m).le.phcut.and.per(m).ge.plcut) then
               resid(n,m) = Sa(1,icomp,m) - Sa(2,icomp,m)
	    endif
  95      continue
          ieq(n) = iev
          ista(n) = isim
          n = n + 1
  98    continue
 100  continue
      close (17)
      close (10)
      n = n - 1


c     Compute the  mean residual (model bias)
      do 120 m=1,nPer
        residBar(m) = 0.0
        do 115 i=1,n
          residBar(m) = residBar(m) + resid(i,m)/n
  115   continue
  120 continue

c     Compute variance of residuals  (random and model uncertainty)
      do 200 m=1,nPer
        sum = 0.0
        do 150 i=1,n
          sum = sum + ( resid(i,m) - residBar(m) )**2
c         print *,m,i,resid(i,m),residBar(m),sum
 150    continue
        sigma(m) = sqrt(sum/(n-1))
        sigmaMean(m) = sigma(m) / sqrt(float(n))
 200  continue

c     Write results
	write (*,'( 2x,''Enter output file for residuals '')')
	read (*,'( a80)') fileout
	iunit = 20
	open (iunit,file=fileout)
	call wrhdr (iunit)
	do 400 i=1,n
	   do 390 m=1,nPer
	      write (iunit,'( i5,2f12.4)') ista(i), per(m), resid(i,m)
390	   continue
400	continue
 
c     Write the mean residuals and the standard error
	write (*,'( 2x,''Enter output file for mean residuals'')')
	read (*,'( a80)') fileout
	iunit = 19
	open (iunit, file=fileout)
	call wrhdr (iunit)
	do 500 m=1,nPer
	   write (iunit,'( 4f10.4)') per(m), residBar(m), sigma(m),
     1       sigmaMean(m)
500	continue

      close (20)
      close (19)
      stop
      end

c ---------------------------------------------------------
	subroutine wrhdr(iunit)

	integer      iunit
	include 'inserts'

	write (iunit,'(''Component =         '',i1)') icomp
	write (iunit,'(i1,''    =   Number of events'')') nEvent
	write (iunit,'(''Events and number of stations:'')')
	do i=1,nEvent
	   write (iunit,'(a20,i5)') evName(i),nSim(i)
	enddo

	return
	end
