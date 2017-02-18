      program EXSIM14  !! Released in January, 2014, Modified version of EXSIM13 by G. Atkinson and K. Assatourians
                       !! the main modification is about implementing frequency dependent geometrical spreading
                       !! to the subsources. It is done following G.Atkinson's recommendation by multiplying a filter
                       !! to the spectra of subsources. Each subsource is considered as separate source with its
                       !! distance from site.      
c     program EXSIM13  !! Released in January 31, 2013, Modified version of EXSIM12 by G. Atkinson and K. Assatourians
c     program EXSIM12  !! Released in May 3, 2012, Modified version of EXSIM_dmb by G. Atkinson and K. Assatourians
c     EXSIM is a stochastic finite-fault program, following the methodology of:
c     Motazedian,Dariush and Gail M. Atkinson, (2005)
c     "Stochastic Finite-Fault Modeling Based on a Dynamic Corner Frequency",
c     Bulletin of the Seismological Society of America,95,995-1010.
c     EXSIM has been modified several times since its original development.
c     This version was modified from "EXSIM_dmb", where
c     EXSIM_dmb  is a version of EXSIM that was modified by D. Boore as described in Boore (2009 BSSA)
c     Boore's modifications addressed some consistency of terminology issues, and improved low-frequencies.
c
c     This EXSIM is designed to work for a list of observation sites for a
c     user-defined fault geometry and location. Reads the input parameters and
c     outputs acceleration time series, and the average Fourier
c     Spectra and PSA over trials.  It also generates averages.
c
c     An example input file of program parameters is given is Example_EXSIM11.par
c     Other input files for crustal amplifications, site amplifications and empirical filter may also be needed.

 
* Dates: Original EXSIM 2003-2005:  based on modifications to FINSIM
*                   FINSIM written by Beresnev and Atkinson (1997 SRL; also BSSA)
*                   FINSIM modified by EXSIM by Motazedian and Atkinson (2005 BSSA)
*                   EXSIM modified by DM Boore as per log file below
*        08/17/08 - Modifications by D. M. Boore
*        08/17/08 - Corrected typo "ISlipWeigth" and "falg", and allow all
*                   unit slip weights if flag -1.
*        08/18/08 - Read in y (=vrup/beta), write R in PSA_FA output  
*                   Revise column headers of PSA_FA output.
*                   Lots of changes to output, including writing
*                   separate PSA_FA files for each site (to make it
*                   easier to import into a graphics program).
*                   Allow the site coordinates to be entered as 
*                        isitecoordflag   coords
*                                     1  lat,long 
*                                     2  R, Az
*                                     3  N, E
*        08/19/08 - Reverse order of fault input, and allow use of Wells
*                   and Coppersmith if FL=0 or FW=0, put magnitude input before
*                   FL, FW (because the later needs amag).  Also get fault
*                   type (needed for Wells and Coppermsith).
*        08/20/08 - Change averaging algorithm (log avg for PSA, E(FAS^2)^1/2 for FAS)
*                   Return E(FAS^2)^1/2, not log of the value, out of the 
*                   binning routine SAMPLE, called by compute FACCN.
*
*                   Add input variable iflagscalefactor to determine whether to use
*                   scale factor based on the following:

*                   iflagscalefactor   integrand of integral
*                           1            velocity spectrum squared
*                           2            acceleration spectrum squared
*                           3            D. Boore's factor based on acc. sq. high-frequency asymptotes
 
*                   Add input variable to determine if use arithmetic or 
*                   geometric average for PSA
*        08/21/08 - Added scalefactor (h) to computeStochasticWave output and print out 
*                   the scalefactor and corner frequency for each subfault for 
*                   the first simulation.  Also include an option for tapering the 
*                   scaling factor.
*        08/22/08 - Write out the time series for the first realization only.
*        08/25/08 - Read in qmin
*                   In checkFFTpoints, increased upper index so that 2**n=32768.  
*                   Also, doubled dimensions of wave, totalWave, and subwave.
*                   In many subroutines replaced dimensions of arrays passed 
*                   in through the parameter list with "(*)".
*        08/27/08 - Added dur_sub to common /sub/, in order to print it out at some point.
*        11/06/08 - Add option to compute RMS of PSA, and also use same options for computing the average of FAS 
*                   (for study purposes).
*        11/07/08 - Changed way that averages are computed.  Now the accumulated sum is computed
*                   in computeAverageResponse and computeAverageFA and this is used in the main program
*                   to compute an average.  But a modification was needed to allow
*                   for negative values of log10(y) for geometric means.  I set the average to -9.99
*                   if accum_fas = 0.0 or accum_psa = 0.0.
*        11/08/08 - I may be having problems with the rms average of FAS, so compute this differently now.
*        11/15/08 - Increase maximum number of nl, nw to 200.
*        11/25/08 - I am concerned about the calls to ran1.   Note that the distributed version of
*                   exsim uses urand, not ran1.  Apparently I susbtituted ran1 for urand, but I think
*                   I made a mistake.  If used to
*                   generate a series of "random" numbers, ran1 should be called the
*                   first time with a negative integer (idum) as a seed.  
*                   Ran1 performs some initializations, including assigning
*                   values to Saved variables iv and iy, as well as changing
*                   idum.  In exsim ran1 is called for a number of things,
*                   possibly including locateHypocentreRandomly, 
*                   createRandomWeight, computeStochasticWave, and
*                   computing a random delay for the time series from 
*                   each subsource.   I would think that the idum, iv, 
*                   and iy values would need to be in sync with one another,
*                   but this will not be the case if ran1 is called with
*                   a different value for idum that is not the same as in
*                   previous calls.   In the previous version of exsim, iseed
*                   is used as idum in the calls in all calls except for
*                   the call for computing the random delay.  But iseed was 
*                   not constrained to be less than 0.0 initially, so the ran1
*                   initialization was probably not done.  And also there was
*                   one call to ran1 that used islipweight for idum, which
*                   could have values less than, equal, to, or grater than 0.0
*                   (in my runs it was lsess than 0.0).  When less than 0.0, as
*                   it was in all of my runs, the iv and iy would have been
*                   initialzed and would have been used in the next call to
*                   ran1, but the idum would not have been the reset value
*                   of islipweight.  I do not know the consequences of all of
*                   this.  I have changed exsim to contrain iseed to be less 
*                   than 0 initially and have used iseed in all calls to ran1.
*        11/27/08 - Added duration calculation
*        11/28/08 - Use 1/f0 for risetime
*        11/28/08 - Change input file to have a comment line before the name of the scale
*                   factor file.
*        11/30/08 - Move computation of random slipweights 
*                   if Islipweight .eq. 1.0 from main program to getInputParameters
*                   (a more logical place because the slipweights array is filled
*                   in the same place for all values of ISlipWeight).
*        11/30/08 - Move "findSubFaultMoment" out of loop over isite
*        12/01/08 - Use "if then" in loop over subsources in 
*                   finding NumberofActiveSubs
*        12/01/08 - Remove "NumberofActivesSubs" from the arguments passed
*                   to subroutines computeStochasticWave and sourceSpectra, as
*                   the variable is not used.
*        12/01/08 - Added writing of istart, istop, nptsTotalWave, maxPointsOfWave
*                   in WritePar (I'm trying to figure out if findindex is needed and if
*                   code setting istop is correct).
*        11/28/08 - Use 1/firstElementF0 for risetime
*        12/01/08 - Use original exsim risetime
*        12/01/08 - Add option to choose one of three ways of obtaining risetime
*        12/06/08 - Allow iteration of random hypocenters for each station
*        12/08/08 - Read initial seed.
*        12/09/08 - Add scaling_coefficient, following Dariush's suggestion.
*        12/09/08 - Increase length of output file names.
*        12/17/08 - Move computations that do not depend on the loop over
*                   nsims out of the loop.  This required defining f0, etc as
*                   arrays (one value per subsource).  This is computationally efficient.
*                   Another advantage is that
*                   I can calculate the maximum duration for the total wave before doing
*                   the subsource wave calculation, and therefore I can use dynamic
*                   allocation for the time series array dimension.  This required
*                   changing common sub not to include arrays of variables that 
*                   depend on the subsource, passing these things through argument lists.
*                   I renamed common /sub/ to common /params/.  I also removed dur_sub
*                   from this common block.
*                   I also compute dur_sub outside of computeStochasticWave and pass the
*                   value into that subroutine through the argument list.
*        12/17/08 - Major revision involving the time series (use tpads and front and back, do not truncate
*                   the motions from each subsource; use dynamic allocations). NOTE: 
*                   Using Mavroeidis and Papageorgiou has not been tested, and it will probably not work 
*                   because I have not taken into account the added npadl.  I comment out a section of code
*                   before the call to M&P that may have to do with determining the index for which the pulse
*                   starts.  It would be easy to fix this later using t_arrive(i,j).
*        12/22/08 - If nl = nw = 1, over-ride the computation of a random hypocenter and 
*                   set n_hypocenters, i0, j0 = 1, 1, 1
*        12/31/08 - Write sqrt(nl*nw)*fas to help in deciding if that factor can correct for
*                   the incoherent summation.
*        12/31/08 - replace "np" for the taper correction at low frequencies by "ftaper_f0", 
*                   the ratio of ftaper and f0 read from the input file.  Also, there is no
*                   need for "taper_scalefactor", because ftaper_f0 .eq. 0.0 implies
*                   no tapering toward low frequencies.
*        01/04/09 - Apply a new taper that attempts to correct for the incoherent summation,
*                   the decay of the subsource spectra for frequencies less than the subsource
*                   corner frequency (following Frankel, 1995), and the improper long-period spectral
*                   level.
*        01/06/09 - Remove "scaling_coefficient" from input and program.  
*                 - Put computation of c4taper, fc4taper outside of the sourcespectra function.
*                 - Revise format of input, using a header line before each line containing input.
*                 - Include a low-cut filter.
*                 - Compute average pga and pgv.
*        01/10/09 - Increase max number of frequencies to 500, without changing anything else.
*        01/13/09 - "apply_taper" was read in from the parameter file, but it has not been used since
*                   I added the Frankel-like scale factor (on 01/04/09).  Thus I have formally removed it from this
*                   version of the program.
*                 - Add SD, PSV to FAS, PSA output file and change order of FAS and PSA.
*        02/06/09 - Compute and write out husid time series and 95%-5% duration
*                 - Change input to be closer to smsim input
*                 - Read only output file stem name
*        02/14/09 - Correct error in setting i0, j0 if n_hypocenters = 1 and nsites > 1
*        02/16/09 - Change input, use stress scaling to adjust Wells and Coppersmith FL, FW
*                 - Compute RJB, RCD, print all distances for each site.
*                 - Include an option to move the sites to positions relative to
*                   the midpoint of the fault or to the tip of the fault
*                   if move_site .eq. .true. and isitecoordflag .eq. 2 (r and epi)
*                   (reset isitecoordflag to 3 in this case) 
*                   or if move_site .eq. .true.  and isitecoordflag .eq. 3 and one of the 
*                   site locations = 0.0
*                 - Moved writePar before loops, and just after call to getInputParameters.
*        02/17/09 - Modify findDistanceAndAzimuth subroutine
*                 - Modify findSubfaultDistance subroutine
*                 - Change FaultStrikeDeg to FaultStrike
*                 - Delete convertDegreeToRadian subroutine
*                 - Delete findDistance subroutine
*                 - Delete shortestSubfaultDistance subroutine
*                 - Read in dl, dw, and compute nl, nw
*        02/18/09 - Replace calls to get_time with a call to the system routine date_time
*        02/23/09 - Hardwire in choice of 75% rather than 95% in determining the duration.
*                   75% is recommended by Ou and Herrmann (SRL) and is consistent with a few
*                   Husid plots I made.
*        05/28/09 - Specify a hypocenter location by distances along the fault and down the dip.
*                 - Add period column to psa_fs output.
*        12/05/09 - Allow for up to a six-segment path duration (Atkinson and Boore 
*                   use a four-segment function).
*                   The form of the input matches that in the smsim programs.
*        12/06/09 - Pass path duration parameters through calling arguments rather than common.
*                 - Obtain number of frequencies in the crustal and site amp files from the 
*                   first line of the files (not from the params file).
*                 - Move the slip weight params before the site coord params and obtain slip weights
*                   from a text file (if slip weights are to be specified).
*                 - Modify order of site coord params so that it is easy to edit if do not want
*                   to do sims at all sites in the list (just change nsites; the list of sites need
*                   not be edited if the first nsites are those for which simulations are desired).
*                 - changed "y" to "vrup_beta" and "v" to "vrup" (the former names makes it very
*                   difficult to search for those variables).
*
*        May/2012  - D.M. Boore version modified by Assatourians to clean up, with some changes in input conventions
*                    These were mostly changes to simplify inputs and outputs, provide additional flexibility for outputs,
*                    allow use of an empirical filter with variability to modify motions;  changes noted throughout
*                    in a commented version of this code that we have (EXSIM11) but not all such comments are retained below.


!      use dflib                                                                 !Modified/added by KA, Aug.2010

      real wave(:), totalWave(:), vel_total(:), husid(:)
      allocatable :: wave, totalWave, vel_total, husid 
      dimension slipWeights(200,200),
     :          slat(100),slon(100),
     :          dist(100),
     :          accum_fas(500), accum_psa(500),
     :          averagelogPSA(500),averagePSA(500),
     :          averagesqFA(500),averageFA(500),nnFA(500)
      dimension weightedMoment(200,200)
      dimension freq(500),psa(500),fa(500),fh(500),hd1(500),hd2(500)
      
      real subfaultMoment(200,200), f0(200,200), risetime(200,200)
      integer NumberOfActiveSubs(200,200) 
      real subfaultDistance(200,200), actualSlip(200,200), 
     :     dur_sub(200,200), dur_path(200,200),
     :     t_arrive(200,200), delay(200,200), freq_out(50), 
     :     freq4intrp(500), psa4intrp(500)
     
      real dur_75_05, arias
      
      real alogFASsum(500), alogPSAsum(500)
      integer nFASsum(500), nPSAsum(500)
      real avgavgFAS(500), avgavgPSA(500), avgavgPSA_out(50)
      character f4head(50)*5
      character *3                         chr3, chr3_
      
      real amag, r_ref, rlow(10), a_s(10), b_s(10), m_s(10)

      real stutter
      
      real afreq1(500), amp1(500), afreq2(500), amp2(500)
      real afreq3(500), amp3(500), sigma(500)
      dimension inormal(10), xrandom1(10), xrandom2(10),
     *          xrmin(10), xrmax(10)
      dimension siteLocation(3000,2)
      character f_crustal_amps*30,f_site_amps*30,f_empirical_amps*30            !Modified/added by KA, Aug.2010
	character InputFileName*30
      character f_stem*120
      character fpar*120, facc*120, fpsa*120, f_h_f0*120, fhusid*120,
     :          f_times*120, f_dur_75_05*120, f_dist_psa*120
      character ftmp*120                                                        !Modified/added by KA, APR.2012
     
      logical f_exist, rmv_trend, specify_length, specify_width,
     :        move_site, write_site_files, AB14_flag                            !Added by KA, Jan.2014 to handle frequency dependent geom_sprd

      character fault_type*10

      character datx*8, time_start*10, time_stop*10

      character*30 fparrandom

      real rpathdur(6), pathdur(6)
      
      integer(2) n_ofargs                                                       !Modified/added by KA, Aug.2010
      integer(4) i_system                                                       !Modified/added by KA, Feb.2012

      common /params/
     :   rho,beta,iKapa,fmax,
     :   Qmin,Q0,Qeta,                                                          !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
!     :   fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q,                             !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :   iwind,eps,eta,                                                         !Modified/added by KA, Mar.2012, eps/eta added
     :   n_crustal_amps,n_site_amps,n_empirical_amps,totalMoment,               !Modified/added by KA, Aug.2010
     :   c4taper, fc4taper,
     :   nsubs, flocut, nslope, pi, twopi
     
      common /par/  iFFTsub, dt, npadl, npadt
      common /seed/ iseed
      common /fnames/ f_crustal_amps,f_site_amps,f_empirical_amps               !Modified/added by KA, Aug.2010
      common /fltinf/FaultStrike,FaultDip,h,FaultLength,FaultWidth,dl,dw        !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
     :       ,AB14_flag                                                          !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
     
      integer soFarNumberOfSubs
      integer nptsTotalWave
c     nptsTotalWave is the total number of points in the final timeseries.
      integer :: io
      character chtest*80                                                       !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files

 ! Initialize some arrays
      pathdur(:) = 0.0
      rpathdur(:) = 0.0

      nu_ctl = 99
      
      ioPar = nu_ctl - 1
      ioPSA = nu_ctl - 2
      ioHUS = nu_ctl - 3
      ioACC = nu_ctl - 4
      ioresp1 = nu_ctl - 5
      ioresp2 = nu_ctl - 6
      ioa = nu_ctl - 7
      nu_h = nu_ctl - 8
      nu_i0j0 = nu_ctl - 9
      nu_times = nu_ctl - 10
      nu_dur_75_05 = nu_ctl - 11
       
      nu_dist_psa = nu_ctl - 20

      iotmp=nu_ctl-25                                                           !Added by Karen for getting tmp output in one file (20104030)
      ioaccall=nu_ctl-30                                                        !Added by Karen for getting acc output in one file (20101015)
      open(unit=ioaccall,file='ACC.out',status='unknown')                       !Added by Karen for getting acc output in one file (20101015)

    
* Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
      call DATE_AND_TIME( datx, time_start )
* Date is returned as 'CCYYMMDD'
* Time is returned as 'hhmmss.sss'
       

      pi=4.0*atan(1.0)
      twopi = 2.0 * pi

!      n_of_args=nargs()                                                         !Modified/added by KA, Aug.2010
      n_of_args=iargc()+1                                                        !Modified/added by KA, Aug.2010

      f_exist = .false.
      do while (.not. f_exist)
        InputFileName = ' '
        if(n_of_args.gt.1)then                                                  !Modified/added by KA, Aug.2010
            call getarg(1,InputFileName)                                        !Modified/added by KA, Aug.2010
        else                                                                    !Modified/added by KA, Aug.2010
           n_of_args=0                                                          !Modified/added by KA, Aug.2010
           write(*, '(a)', advance='no') 
     :       ' Enter name of input parameter file '//
     :       '(cr = EXSIM12.params): '                                          !Modified/added by KA, Apr.2012
           read(*, '(a)') InputFileName
        endif                                                                   !Modified/added by KA, Aug.2010
        if (InputFileName(1:4) .eq. '    ')                                     !Modified/added by KA, Apr.2012
     :                 InputFileName = 'EXSIM12.params'                         !Modified/added by KA, Apr.2012
        call trim_c(InputFileName, nc_f_in)
        inquire(file=InputFileName(1:nc_f_in), exist=f_exist)
        if (.not. f_exist) then
          write(*,'(a)') ' ******* FILE DOES NOT EXIST ******* '
        end if
      end do
      
      call check_directories()                                                  !Modified/added by KA, Mar.2012: Check/creat directories for outputs
      i_system=system('rm ./psa.out')                                          !Modified/added by KA, Jan.2013: to remove previously created PSA.out file
      open(unit=51,file='psa.out')                                              !Modified/added by KA, Jan.2013: to create PSA.out file
      close(unit=51)                                                            !Modified/added by KA, Jan.2013: to create PSA.out file
      

      
      open(unit=nu_ctl, file=InputFileName(1:nc_f_in), status='unknown')

      print *, ' Control file: '//InputFileName(1:nc_f_in)
c     Read (and echo) the input parameters that are common to all grid points.
      call getInputParameters(nu_ctl,
     :     FaultStrike,FaultDip,h,
     :     FaultLat,FaultLon, move_site, write_site_files,
     :     siteLocation,numberOfSites,isitecoordflag,
     :     fault_type,
     :     FaultLength, FaultWidth, 
     :     dl, dw, nl,nw, nsubs, 
     :     specify_length, specify_width, stress_ref, 
     :     vrup_beta,
     :     hyp_loc_fl, hyp_loc_fw, i0_in,j0_in,n_hypocenters,
     :     tpadl, tpadt, dt, beta, rho, amag,
     :     stress,iKapa,fmax,
!     :     fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q,                           !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :     Qmin,Q0,Qeta,                                                        !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :     r_ref, nsprd_segs, rlow, a_s, b_s, m_s,
     :     rpathdur, pathdur, durslope, ndur_hinges, 
     :     iwind,eps,eta,                                                       !Modified/added by KA, Mar.2012, eps/eta added
     :     flocut, nslope,
     :     nfreq,freq1,freq2,
     :     nfout, freq_out,
     :     f_stem,
     :     f_crustal_amps,  
     :     f_site_amps,
     :     f_empirical_amps,                                                    !Modified/added by KA, Aug.2010
     :     iseed, nsims,damp,
     :     islipweight, slipWeights,
     :     iRow,jColumn,pulsingPercent,
     :     iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak,
     :     iDynamicFlag, 
     :     iflagscalefactor, 
     :     iflagfas_avg, iflagpsa_avg,
     :     i_rise_time,
     :     AB14_flag)                                                           !Added by KA, Jan.2014 to handle frequency dependent geom_sprd

      close(nu_ctl)
      
!     Parsing other input files, at startup

      nu_crustal_amps = 2
      call trim_c(f_crustal_amps, nc_f_crustal_amps)
      open (nu_crustal_amps,
     :     file=f_crustal_amps(1:nc_f_crustal_amps),status='unknown')
!     do while(not(eof(nu_crustal_amps)))                                     !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      do                        !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         read(nu_crustal_amps,'(a80)', iostat=io)chtest !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         if (io < 0) then       ! End of file reached
            exit                ! Leave the loop
         endif
         if(chtest(1:1)/='!')then !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
            backspace(nu_crustal_amps) !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
            exit                !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         endif                  !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      enddo                     !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      read(nu_crustal_amps,*) n_crustal_amps
      do i=1,n_crustal_amps
         read (nu_crustal_amps,*) afreq1(i),amp1(i)
      enddo
      close (nu_crustal_amps)
        
      nu_site_amps = 2
      call trim_c(f_site_amps, nc_f_site_amps)
      open (nu_site_amps,
     :     file=f_site_amps(1:nc_f_site_amps),status='unknown')
!     do while(not(eof(nu_site_amps)))                                        !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      do                        !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         read(nu_site_amps,'(a80)', iostat=io)chtest !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         if (io < 0) then       ! End of file reached
            exit                ! Leave the loop
         endif
         if(chtest(1:1)/='!')then !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
            backspace(nu_site_amps) !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
            exit                !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         endif                  !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      enddo                     !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      read(nu_site_amps,*) n_site_amps
      do i=1,n_site_amps
         read (nu_site_amps,*) afreq2(i),amp2(i)
      enddo
      close (nu_site_amps)
        
      nu_empirical_amps = 2     !Modified/added by KA, Aug.2010
      call trim_c(f_empirical_amps, nc_f_empirical_amps) !Modified/added by KA, Aug.2010
      open (nu_empirical_amps,  !Modified/added by KA, Aug.2010
     :     file=f_empirical_amps(1:nc_f_empirical_amps),
     :     status='unknown')    !Modified/added by KA, Aug.2010
!     do while(not(eof(nu_empirical_amps)))                                   !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      do                        !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         read(nu_empirical_amps,'(a80)', iostat=io)chtest !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         if (io < 0) then       ! End of file reached
            exit                ! Leave this loop
         endif
         if(chtest(1:1)/='!')then !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
            backspace(nu_empirical_amps) !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
            exit                !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
         endif                  !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      enddo                     !Modified/added by KA, Aug.2010, and Apr. 2012 for skipping comment lines in the beginning of filter files
      read(nu_empirical_amps,*) n_empirical_amps !Modified/added by KA, Aug.2010
      do i=1,n_empirical_amps   !Modified/added by KA, Aug.2010
         read (nu_empirical_amps,*) afreq3(i),amp3(i),sigma(i) !Modified/added by KA, Aug.2010, and Feb. 2012
      enddo                     !Modified/added by KA, Aug.2010
      close (nu_empirical_amps) !Modified/added by KA, Aug.2010

!     End modifications

      vrup=vrup_beta * beta

      subFaultRadius=sqrt((dl*dw)/pi)

!%%%%%%%%%%
      riseTime_original=subFaultRadius/vrup
!%%%%%%%%%%

      call trim_c(f_stem, nc_f_stem)
      
      fpar = ' '
      fpar = f_stem(1:nc_f_stem)//'_parameters.out'
      call trim_c(fpar, nc_fpar)
      
* Initialize q_f    (freq-dependent Q model)
      dummy =  q_f_setup(Qmin,Q0,Qeta)                                          !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
!      dummy = q_f_setup(fr1, qr1, s1, ft1, ft2,                                 !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
!     :                  fr2, qr2, s2, c_q)                                      !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition

* Initialize gsprd  (geometric sprading)

      numsource = 1
      dummy = gsprd_f_setup(r_ref,nsprd_segs,rlow,
     :                  a_s,b_s,m_s,
     :                  numsource,amag)
     
      
      npadl = tpadl/dt
      npadt = tpadt/dt


      if(nfreq.eq.-1)then
          call SCEC_freqs(nfreq,freq)
      else
c     Generate some logarithmically-spaced frequencies
          freq(1)=freq1
          if (nfreq.ge.2) then
             frinc=alog(freq2/freq1)/(nfreq-1)
             do k=2,nfreq
                freq(k)=freq1*exp((k-1)*frinc)
             enddo
          endif
c     End of Generate some logarithmically-spaced frequencies
      endif

!      iseed = -iabs(iseed)

      if(i0_in.eq.0.or.j0_in.eq.0)
     :    call locateHypocentreRandomly(i0,j0,nl,nw)

      NoOfEffectiveSubfaults=float(nl)*(pulsingPercent/100.0)
      NoOfEffectiveSubfaults=NoOfEffectiveSubfaults/2.                          !!double side propogation
      if(NoOfEffectiveSubfaults.le.1) NoOfEffectiveSubfaults=1

      totalMoment=10.**(1.5*amag+16.05)
      avgSubFaultMoment=totalMoment/float(nw*nl)

      firstElementF0=  
     :                4.9e+6*beta*(stress/avgSubFaultMoment)**(1.0/3.0)

      F0main=4.9e+6*beta*(stress/totalMoment)**(1.0/3.0)

      call findSubFaultMoment(nl,nw,slipWeights
     :                 ,totalMoment,weightedMoment)

   
      f4head = '00000'
      
      do i = 1, nfout
        if (freq_out(i) .lt. 0.0 .or. freq_out(i) .ge. 10.0) then
           write(f4head(i)(1:5),'(f5.2)') freq_out(i)
        else  
           write(f4head(i)(2:5),'(f4.2)') freq_out(i)
        end if
      end do
        
      f_dur_75_05 = ' '
      f_dur_75_05 = f_stem(1:nc_f_stem)//'_dur_75_05_dist_psa.out'
      call trim_c(f_dur_75_05, nc_f_dur_75_05)
      open(nu_dur_75_05, file='./other/'//f_dur_75_05(1:nc_f_dur_75_05),        !Modified/added by KA, Feb.2012
     :             status='unknown')
      write(nu_dur_75_05,'(2x,a,1x,a,2x,a,2x,a,1x,a,1x,a,
     :             1x,a, 1x,a,1x,a,1x,a,1x,a,5x,a, 3x,a,4x,a,   
     :             50(7x,a4, 3x,a8))') 
     :             'isite','ihyp','i0','j0','dur_75_05(sec)',
     :             'arias(cm/s)','sitecoord(1)', 'sitecoord(2)', 
     :             'isitecoordflag','site_lat_degrees', 
     :             'site_lon_degrees', 'd_jb', 'd_cd2f','d_hyp',
     :             ('freq', 'psa'//f4head(i), i = 1, nfout)
 

      DO isite=1,numberOfSites   ! Loop on Sites around the Fault
      
        SiteLat=siteLocation(isite,1)
        SiteLon= siteLocation(isite,2)
        print *, "Working on Site #",isite,SiteLat, SiteLon

* Compute some distances: 

        call site_coords_in_degrees(
     :         FaultLat, FaultLon, 
     :         SiteLat, SiteLon, isitecoordflag, 
     :         site_lat_degrees, site_lon_degrees ) 

        h_min_c = 3.0 ! Campbell depth to seismogenic region
        w1 = 0.0
        w2 = faultwidth
        s1 = 0.0
        s2 = faultlength
     
        call dist_3df(
     :   site_lat_degrees, site_lon_degrees, 
     :   FaultLat, FaultLon, h, FaultStrike, FaultDip,
     :   w1, w2, s1, s2, 
     :   h_min_c, d_jb, az_jb, d_cd2f, az_cd2f, d_c, az_c,
     :   d_sta_n, d_sta_e, irgn_cd2f, irgn_c, irgn_jb)
 
        nnFA=0
        totalDuration=0
        fi2=0
        call findDistanceAndAzimuth(FaultLat,FaultLon,SiteLat,
     :                           SiteLon,R,fi2,isitecoordflag)
c     Note that R, azm fi2 are the distance, azm. w.r.t. origin (not epi),
c     Note also I must input faultstrike in degrees and its converted to rad.
 
        alogFASsum = 0.0
        alogPSAsum = 0.0
           nFASsum = 0
           nPSAsum = 0
        alogPGAsum = 0.0
        alogPGVsum = 0.0
        
        HypoDistance_mean=0.                                                    !Added by KA in Feb. 2013 to find the average hypocentral distance and report it in *psa.out files
        
        DO ihypo = 1, n_hypocenters

          if (n_hypocenters .gt. 1) then
            call locateHypocentreRandomly(i0,j0,nl,nw)
          else if (i0_in. eq. 0 .or. j0_in .eq. 0) then
            call locateHypocentreRandomly(i0,j0,nl,nw)
          else
            i0 = i0_in
            j0 = j0_in
          end if

          nsubfaults = 0
          
          t_arrive_min=10000.
          t_end_max=0.
          risetime_max = 0.0
            
          DO i=1,nl
            DO j=1,nw
            
              nsubfaults = nsubfaults + 1
              
              subFaultMoment(i,j)=weightedMoment(i,j)
                
              if(iDynamicFlag.eq.0) then
                NumberOfActiveSubs(i,j)=1
              else
                NumberOfActiveSubs(i,j)=NumberOfPulsingSubs(i0,j0,i,j
     :                 ,nl,nw,NoOfEffectiveSubfaults)
                if(NumberOfActiveSubs(i,j).eq.0) then
                  NumberOfActiveSubs(i,j)=1
                end if
              end if
                
              f0(i,j)=
     :           firstElementF0*(NumberOfActiveSubs(i,j)**(-1.0/3.0))
              if (i_rise_time .eq. 1) then
                risetime(i,j) = riseTime_original                
              else if (i_rise_time .eq. 2) then
                riseTime(i,j) = 1.0/f0(i,j)
              else
                print *,' ERROR: invalid value of i_rise_time, = ',
     :                    i_rise_time
                stop
              end if
              
              if(risetime(i,j) .gt. risetime_max) then
                risetime_max = risetime(i,j)
              end if
 
              actualSlip(i,j)=1.e-22*subFaultMoment(i,j)/
     :                            (rho*beta**2.*dl*dw)
              subfaultDistance(i,j)=
     :                 findSubfaultDistance(R,h,FaultStrike,
     :                           fi2,FaultDip,dl,dw,i,j)
     
c     calculate duration of subsource time history at distance
c     subfaultDistance
              call dur_path_cmp(subfaultDistance(i,j),
     :                        rpathdur, pathdur, durslope, ndur_hinges,
     :                  dur_path(i,j))  ! computes path part of duration

              dur_sub(i,j) = dur_path(i,j) + risetime(i,j)
              
              if ((i-i0) .ne. 0 .or. (j-j0) .ne. 0) then
                delay(i,j)=sqrt((dl*(i-i0))**2.+(dw*(j-j0))**2.)/vrup
              else
                delay(i,j)=0.  ! delay from hypocenter to itself
              end if
              
              t_arrive(i,j) = delay(i,j) + subfaultDistance(i,j)/beta
              if (t_arrive(i,j) .lt. t_arrive_min) 
     :                    t_arrive_min = t_arrive(i,j)
               
              t_end= t_arrive(i,j) + dur_sub(i,j)            
              if (t_end .gt. t_end_max) then
                imax = i
                jmax = j
                t_end_max = t_end
              end if
            
            END DO                                                              ! loop over nw to calculate quantities not dependent on nsims
          END DO                                                                ! loop over nl to calculate quantities not dependent on nsims

          nptsTotalWave=(t_end_max - t_arrive_min + tpadl + tpadt + 
     :                        risetime_max)/dt                                  ! need risetime_max because the stutter
                                                                                ! delay can shift the subsource time series
                                                                                ! corresponding to t_end_max by an amount up to
                                                                                ! risetime of that subsource; to be safe I
                                                                                ! use risetimemax rather than risetime(imax,jmax),
                                                                                ! where imax, jmax are the indices of the subsource
                                                                                ! corresponding to t_end_max.
                                                 
          signnpw2 = +1.0
          call get_npw2(nptsTotalWave,signnpw2,npw2TotalWave)
          
          if (write_site_files) then
          
            if (ihypo .eq. 1) then
 
              f_times = ' '
              f_times = f_stem(1:nc_f_stem)//'_times_s000.out'
              write(f_times(nc_f_stem+9:nc_f_stem+11),'(i3.3)') isite
              call trim_c(f_times, nc_f_times)         
       
              end if
              
          end if
          
          HypoDistance=findSubfaultDistance(R,h,FaultStrike,
     :                           fi2,FaultDip,dl,dw,i0,j0)
          HypoDistance_mean=HypoDistance_mean +                                 !Added by KA in Feb. 2013 to find the average hypocentral distance and report it in *psa.out files
     :                      HypoDistance / real(n_hypocenters)                  !Added by KA in Feb. 2013 to find the average hypocentral distance and report it in *psa.out files

c     now do actual generation and summation of subfault time series ***

! Initialize arrays (index over frequency) containing cumulative sums for averages
          averagelogPSA=0.0
          averagePSA=0.0
          averageFA=0.0
          accum_fas = 0.0
          accum_psa = 0.0
          nnFA = 0
          accum_logpga = 0.0
          accum_logpgv = 0.0

          accum_dur_75_05 = 0.0
          accum_arias = 0.0

          DO isim=1,nsims
            print *, ' Site isite, Hypocenter iteration ihypo, '//
     :               'Simulation number isim = ', isite, ihypo, isim

            allocate (totalWave(npw2TotalWave))
            allocate (husid(npw2TotalWave))

            totalWave = 0.0

            if(islipweight.eq.2)then                                            !Added by KA in Feb. 2013 to consider option 2, for variable slip distribution on fault from iteration to iteration
                call createRandomWeight(nl,nw,slipWeights)                      !Added by KA in Feb. 2013 to consider option 2, for variable slip distribution on fault from iteration to iteration
                call findSubFaultMoment(nl,nw,slipWeights                       !Added by KA in Feb. 2013 to consider option 2, for variable slip distribution on fault from iteration to iteration
     :                                    ,totalMoment,weightedMoment)          !Added by KA in Feb. 2013 to consider option 2, for variable slip distribution on fault from iteration to iteration
                subFaultMoment=weightedMoment                                   !Added by KA in Feb. 2013 to consider option 2, for variable slip distribution on fault from iteration to iteration
                actualSlip=1.e-22*subFaultMoment/                               !Added by KA in Feb. 2013 to consider option 2, for variable slip distribution on fault from iteration to iteration
     :                            (rho*beta**2.*dl*dw)                          !Added by KA in Feb. 2013 to consider option 2, for variable slip distribution on fault from iteration to iteration
            endif                                                               !Added by KA in Feb. 2013 to consider option 2, for variable slip distribution on fault from iteration to iteration
            

            DO i=1,nl
              DO j=1,nw
            
                nsubsource = npadl + dur_sub(i,j)/dt + npadt
                signnpw2 = +1.0
                call get_npw2(nsubsource,signnpw2,iFFTsub) 

                allocate (wave(iFFTsub))
                
                wave=0
                
                call computeStochasticWave(
     :               wave,
     :               subFaultMoment(i,j), subfaultDistance(i,j),
     :               F0Main, f0(i,j), risetime(i,j), dur_sub(i,j),
     :               nl,nw,j, 
     :               iflagscalefactor, hij,
     :               afreq1,afreq2,afreq3,amp1,amp2,amp3,sigma)          ! dmb
      
                
c               subsource accelerogram is randomly delayed to
c               simulate complexity in slip process
c               subtract minimum delay to eliminate static time shift

                 
                stutter = Ran1(iseed) * riseTime(i,j)
                ishift = ((t_arrive(i,j) - t_arrive_min) + stutter)/dt
                
                if (npadl+ishift+dur_sub(i,j)+npadt .gt. 
     :                                             npw2TotalWave) then
                   print *,'npadl+ishift+dur_sub(i,j)+'//
     :                     'npadt.gt.npw2TotalWave'
                   print *,  npadl+ishift+dur_sub(i,j)+npadt, 
     :                       npw2TotalWave
                end if   
                 
! No need to call shift---just figure out the shifted index and use it  
! to fill the totalWave vector.
                 
c               add current accelerogram to total wave field

                nhigh = min(iFFTsub, npw2TotalWave - ishift)

                DO k = 1, nhigh
                  totalWave(k+ishift)=totalWave(k+ishift)+wave(k)
                END DO
                
                deallocate (wave)
                
              END DO ! loop over nw
            END DO ! loop over nl
          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     include Mavroeidis and Papageorgiou, 2003 approach
            if( int(iPapaFlag).eq.1)then
              call MavroPapa(totalWave,npw2TotalWave,
     *          iPapaFlag,PapaGama,PapaNu,PapaT0,amag,PapaPeak)
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     ************************** end of summation **********************

c     accelerogram and model parameters are saved only once

            pga=findPeakValue(npw2TotalWave,totalWave)
            accum_logpga = accum_logpga + alog10(pga)

            allocate (vel_total(npw2TotalWave))
            vel_total = 0.0

          
            call acc2v(totalWave, npw2TotalWave, dt, .false., vel_total)

            pgv=findPeakValue(npw2TotalWave,vel_total)
            accum_logpgv = accum_logpgv + alog10(pgv)
            
            deallocate(vel_total)

            rmv_trend = .false.
            husid = 0.0 
            call accsqint(totalWave, npw2TotalWave, dt, rmv_trend, 
     :                    husid)
            husid_end = husid(npw2TotalWave)
            do i = 1, npw2TotalWave
              husid(i) = husid(i)/husid_end
            end do
      
            call locate(husid,npw2TotalWave,0.05,j05) 
            call locate(husid,npw2TotalWave,0.75,j75) 
        
            accum_dur_75_05 = accum_dur_75_05 + 
     :                        alog10(float(j75-j05)*dt)                         ! for geometric mean
            
            accum_arias = accum_arias + 
     :                    alog10((pi/(2.0*981.0))*husid_end)                    ! for geometric mean 
              
              
          if (write_site_files) then                                            !Added by Karen for getting acc output of individual simulations (20120415)

              facc = ' '
              write(chr3,'(i3.3)')isite                                         !Added by Karen for getting acc output of individual simulations (20120415)
              write(chr3_,'(i3.3)')isim                                         !Added by Karen for getting acc output of individual simulations (20120415)
              facc = trim(f_stem)//'S'//chr3//'iter'//chr3_//'.acc'             !Added by Karen for getting acc output of individual simulations (20120415)
              call trim_c(facc,nc_facc)
              open (ioAcc,file='./ACC/'//facc(1:nc_facc),                       !Added by Karen for getting acc output of individual simulations (20120415)
     :              status='unknown')
              call writeAcc(ioAcc, totalWave, npw2TotalWave, pga,               !Added by Karen for getting acc output of individual simulations (20120415)
     :                 isim,isite,HypoDistance,d_cd2f, fpar,                    !Added by Karen for getting acc output of individual simulations (20120415)
     :                 amag,siteLat,siteLon,h,r,fi2*180./pi)                    !Added by Karen for getting acc output of individual simulations (20120415)
              close(ioAcc)                                                      !Added by Karen for getting acc output of individual simulations (20120415)
          endif

!            if (ihypo.eq.1 .and. isim .eq. 1) then                              !Modified/added by KA, Mar.2012 to include all husid time series
            if (ihypo.eq.1) then                                                !Modified/added by KA, Mar.2012 to include all husid time series
            
              if(isim.eq.1)then
                  call writeAcc(ioaccall, totalWave, npw2TotalWave, pga,        !Added by Karen for getting acc output in one file (20101015)
     :                 isim,isite,HypoDistance,d_cd2f, fpar,                    !Added by Karen for getting acc output in one file (20101015)
     :                 amag,siteLat,siteLon,h,r,fi2*180./pi)                    !Added by Karen for getting acc output in one file (20101015)
              endif
              
          if (write_site_files) then

              fhusid = ' '
              fhusid = f_stem(1:nc_f_stem)//'_husid_s000iter000.out'
              write(fhusid(nc_f_stem+9:nc_f_stem+11),'(i3.3)') isite
              write(fhusid(nc_f_stem+16:nc_f_stem+18),'(i3.3)') isim
              call trim_c(fhusid, nc_fhusid)         
              open (ioHus,file='./other/'//fhusid(1:nc_fhusid),                 !Modified/added by KA to direct output in folder "other", Feb.2012
     1                  status='unknown')
              
              dur_75_05_sim1 =  float(j75-j05)*dt 
            
              arias_sim1 = (pi/(2.0*981.0))*husid_end
              
              call writeHUS(ioHus, husid, npw2TotalWave, 
     :                 arias_sim1, dur_75_05_sim1, 
     :                 isim,isite,HypoDistance,d_cd2f,fpar,
     :                 amag,siteLat,siteLon,h,r,fi2*180./pi)                    !Added by Karen for reporting earthquake parameters in husid files, May.2012
              close(ioHUS)
              
            end if
            
          end if

            dur=float(npw2TotalWave-1)*dt-5.0
            psa=0.0
            call computeResponse(totalWave,dt,dur,damp,nfreq,freq,psa)
            call computeAverageResponse(accum_psa,psa,nfreq,nsims,
     :                                 iflagpsa_avg)  ! dmb
            FA=0.0
            call computeFACCN(totalWave,npw2TotalWave,dt,nfreq,freq,FA)
            call computeAverageFA(accum_fas,FA,nfreq,nsims,nnFA,
     :                                 iflagfas_avg)
     
            write (*,"(' *** finished with trial',i4)") isim
            
            deallocate(totalWave)
            deallocate(husid)
         
          if (write_site_files) then
            ftmp = ' '
            ftmp = trim(f_stem)//'S'//chr3//'iter'//chr3_//'.tmp'               !Added by Karen for getting spectral output of every simulation (20120430)
            call trim_c(ftmp,nc_ftmp)                                           !Added by Karen for getting spectral output of every simulation (20120430)
            open (iotmp,file='./PSA/'//ftmp(1:nc_ftmp),status='unknown')        !Added by Karen for getting spectral output of every simulation (20120430)
                                                                                !Added by Karen for getting spectral output of every simulation (20120430)
            call writePSA_FA(iotmp, freq,PGA,PGV,                               !Added by Karen for getting spectral output of every simulation (20120430)
     :           psa,nfreq,isim,0,damp,amag,                                    !Added by Karen for getting spectral output of every simulation (20120430)
     :           siteLat,siteLon,h,r,fi2*180./pi,d_cd2f,FA,                     !Added by Karen for getting spectral output of every simulation (20120430)
     :           isite, d_jb, HypoDistance,fPar, nl, nw)                        !Added by Karen for getting spectral output of every simulation (20120430)!Modified/Added by KA, Jan. 2013 to add Joyner-Boore distance
            call write_slip(iotmp,actualSlip,nl,nw)
            close(iotmp)                                                        !Added by Karen for getting spectral output of every simulation (20120430)
          endif

          END DO  ! End loop over nsims

          averagePGA = 10.0**(accum_logpga/float(nsims))
          alogPGAsum = alogPGAsum + alog10(averagePGA)
          averagePGV = 10.0**(accum_logpgv/float(nsims))
          alogPGVsum = alogPGVsum + alog10(averagePGV)


          DO i = 1, nfreq
            if (iflagpsa_avg .eq. 2 .and. accum_psa(i) .eq. 0.0) then
              averagePSA(i) = -9.99
            else
              if (iflagpsa_avg .eq. 1) then  ! arithmetic
                averagePSA(i) = accum_psa(i)/float(nsims)
              else if (iflagpsa_avg .eq. 2) then  ! geometric
                averagePSA(i) = 10.0**(accum_psa(i)/float(nsims))
              else if (iflagpsa_avg .eq. 3) then  ! rms
                averagePSA(i) = sqrt(accum_psa(i)/float(nsims))
              else
                print *,
     :         ' ERROR, iflagpsa_avg not 1, 2, or 3; = ', iflagpsa_avg
                stop
              end if
            end if
             
            if (iflagfas_avg .eq. 2 .and. accum_fas(i) .eq. 0.0) then   
              averageFA(i) = -9.99
            else
              if (iflagfas_avg .eq. 1) then  ! arithmetic
                averageFA(i) = accum_fas(i)/float(nnFA(i))
              else if (iflagfas_avg .eq. 2) then  ! geometric
                averageFA(i) = 10.0**(accum_fas(i)/float(nnFA(i)))
              else if (iflagfas_avg .eq. 3) then  ! rms
                averageFA(i) = sqrt(accum_fas(i))
              else
                print *,
     :         ' ERROR, iflagfas_avg not 1, 2, or 3; = ', iflagfas_avg
                stop
              end if
            end if
         
            if (averagePSA(i) .gt. 0.0) then
              alogPSAsum(i) = alogPSAsum(i) + alog10(averagePSA(i))
              nPSAsum(i) = nPSAsum(i) + 1
            end if
             
            if (averageFA(i) .gt. 0.0) then
              alogFASsum(i) = alogFASsum(i) + alog10(averageFA(i))
              nFASsum(i) = nFASsum(i) + 1
            end if
             
          END DO  ! loop over nfreq

          dur_75_05 = 10.0**(accum_dur_75_05/float(nsims))
          arias     = 10.0**(accum_arias/float(nsims))

          if(ihypo==n_hypocenters)then
            write(nu_dur_75_05,'(4x,i3,
     :                           2x,i3,
     :                           1x,i3, 1x,i3, 
     :                           6x,f9.2,                                       !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :                           1p, 1x,e10.3)',advance='no')                   !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :                           isite,
     :                           ihypo,
     :                           i0, j0,
     :                           dur_75_05,
     :                           arias
          else
            write(nu_dur_75_05,'(4x,i3,
     :                           2x,i3,
     :                           1x,i3, 1x,i3, 
     :                           6x,f9.2,                                       !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :                           1p, 1x,e10.3)')                                !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :                           isite,
     :                           ihypo,
     :                           i0, j0,
     :                           dur_75_05,
     :                           arias
         endif

        END DO !loop over random hypocenters
       
* Compute average of average FAS and PSA


        avgavgPGA = 10.0**(alogPGAsum/float(n_hypocenters))
        avgavgPGV = 10.0**(alogPGVsum/float(n_hypocenters))

        DO i = 1, nfreq
       
         avgavgPSA(i) = 10.0**(alogPSAsum(i)/float(nPSAsum(i)))
         avgavgFAS(i) = 10.0**(alogFASsum(i)/float(nFASsum(i)))
         
        END DO ! end loop over nfreq
       
        nfreq4intrp = nfreq
        do i = 1, nfreq4intrp
          freq4intrp(i) = freq(i)
          psa4intrp(i) = avgavgPSA(i)
        end do
        
        do i = 1, nfout
        
          if (freq_out(i) .eq. -1.0) then ! pgv
            avgavgPSA_out(i) = avgavgPGV
          else if (freq_out(i) .ge. 99.0) then ! pga
            avgavgPSA_out(i) = avgavgPGA
          else
            avgavgPSA_out(i) = 
     :       yintrf( freq_out(i), freq4intrp, psa4intrp, nfreq4intrp)
          end if
          
        end do
 

        fpsa = ' '
        fpsa = f_stem(1:nc_f_stem)//'_psa_fs_s000.out'
        write(fpsa(nc_f_stem+10:nc_f_stem+12),'(i3.3)') isite
        call trim_c(fpsa, nc_fpsa)         
c        open (ioPsa,file=fpsa(1:nc_fpsa),status='unknown')                      !Modified/added by KA, Mar.2012 to include average spectra in PSA folder
        open (ioPsa,file='./PSA/'//fpsa(1:nc_fpsa),status='unknown')            !Modified/added by KA, Mar.2012 to include average spectra in PSA folder
        call writePSA_FA(ioPSA, freq,avgavgPGA,avgavgPGV,
     :       avgavgPSA,nfreq,nsims,1,damp,amag,                                 !Modified/added by KA, Mar.2012 to include average spectra in PSA folder
     :       siteLat,siteLon,h,r,fi2*180./pi,d_cd2f,avgavgFAS,      
     :       isite, d_jb, HypoDistance_mean, fPar, nl, nw)                      !Modified/Added by KA, Jan. 2013 to add Joyner-Boore distance and report mean hypocentral distance
        close(ioPSA)
!        i_system=system('rm ./psa.out')
        i_system=system('cat ./PSA/'//trim(fpsa(1:nc_fpsa))//        !Modified/added by KA, Jan.2013: to create PSA.out file and append data to it
     1                  ' >> ./psa.out')                                           !Modified/added by KA, Jan.2013: to create PSA.out file and append data to it
        write(*,*)'cat ./PSA/'//trim(fpsa(1:nc_fpsa))//
     1            ' >> ./psa.out'


        write(nu_dur_75_05, '(                                                  !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :    5x,f9.3, 4x,f9.3, 14x,i1,                                             !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :    8x,f9.3, 8x,f9.3,                                                     !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :    1x,f8.2, 1x,f8.2, 1x, f8.2,                                           !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected!Modified/Added by KA, Jan. 13, to add hypocentral distance to the summary file
     :    1p,50(1x,e10.3, 1x,e10.3))')                                          !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :    SiteLat,  SiteLon, isitecoordflag,                                    !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :    site_lat_degrees,  site_lon_degrees,                                  !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
     :     d_jb, d_cd2f, HypoDistance,                                          !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected!Modified/Added by KA, Jan. 13, to add hypocentral distance to the summary file
     :     (freq_out(i), avgavgPSA_out(i), i= 1, nfout)                         !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected
      END DO  ! End loop over sites                                             !Modified/added by KA, Apr.2012 to join summary PSA with duration/arias intensity this feature will work correctly if non-random hypocenter is selected

c
c
c
      open (ioPar,file='./other/'//fpar(1:nc_fpar),status='unknown')            !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly

      call writePar(ioPar,                                                      !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               FaultStrike,FaultDip,h,                                    !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               FaultLat, FaultLon,                                        !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               siteLocation,numberOfSites, move_site,                     !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               fault_type,                                                !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               FaultLength, FaultWidth,nl,nw,dl,dw,                       !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               specify_length, specify_width, stress_ref,                 !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               vrup_beta,                                                 !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :            hyp_loc_fl, hyp_loc_fw, i0_in,j0_in,n_hypocenters,            !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               iseed, nsims,                                              !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               amag,                                                      !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               stress,                                                    !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               pulsingPercent,iDynamicFlag,                               !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               i_rise_time,                                               !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak,                 !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               iflagscalefactor,                                          !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               iflagfas_avg, iflagpsa_avg,                                !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               tpadl, tpadt,                                              !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               r_ref, nsprd_segs, rlow, a_s, b_s, m_s,                    !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :            rpathdur, pathdur, durslope, ndur_hinges,                     !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :               actualSlip                                                 !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
     :                             )                                            !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
                                                                                !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
      close(ioPar)                                                              !Karen&Katsu2010/11/10-Moved to this section to allow all parameters to be reported properly
c
c
c
      print *          
* Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
      call DATE_AND_TIME( datx, time_stop )
* Date is returned as 'CCYYMMDD'
* Time is returned as 'hhmmss.sss'
      call time_diff(time_start, time_stop, time_elapsed)
      print *,  
     :       ' Elapsed time (sec): ', time_elapsed

!      close(nu_h)                                                               !Modified/added by KA, Mar.2012, remove reporting values in a file
      close(nu_dur_75_05)
      close(nu_dist_psa)
      close(ioaccall)                                                           !Added by Karen for getting acc output in one file (20101015)
      STOP
      END
cccccccccccccccccccccccc   End of Main Program    cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------- SUBPROGRAMS -----------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real function amplf (f,afreq,amp,namp)
c       computes amplification value at frequency f by linear
c       interpolation between its adjoining values
        dimension amp(*),afreq(*)

        if (f.le.afreq(1)) then
            amplf=amp(1)
            return
        endif

        if (f.gt.afreq(namp)) then
            amplf=amp(namp)
            return
        endif

        do i=2,namp
            if (f.le.afreq(i)) then
                d=(amp(i)-amp(i-1))/(afreq(i)-afreq(i-1))*(f-afreq(i-1))
                amplf=amp(i-1)+d
                return
            endif
        enddo

        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine avg(n,s,avamp2)
c     calculates average squared amplitude spectrum
      complex s
      dimension s(*)
      sum=0.
      do j=1,n
          amp=cabs(s(j))
          sum=sum+amp*amp
      enddo
      avamp2=sum/float(n)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeAverageFA(accum_fas,FA,nfreq,nsims,nn,
     :                                  iflagfas_avg)
c     calculate average Fourier spectrum    ! revised by dmb
      dimension accum_fas(*),FA(*), nn(500)
      do i=1,nfreq
        if(FA(i).gt.0)then
          nn(i)=nn(i)+1
          if (iflagfas_avg .eq. 1) then ! arithmetic
                accum_fas(i)=
     :              accum_fas(i) + FA(i)
          else if (iflagfas_avg .eq. 2) then ! geometric
                accum_fas(i)=
     :              accum_fas(i) + alog10(FA(i))
          else if (iflagfas_avg .eq. 3) then ! rms
            accum_fas(i)=
     :         (accum_fas(i)*float(nn(i)-1) + FA(i)**2.0 )/float(nn(i))
!          else if (iflagfas_avg .eq. 3) then ! rms
!                accum_fas(i)=
!     :              accum_fas(i) + FA(i)**2.0
          else
            print *,
     :         ' ERROR, iflagfas_avg not 1, 2, or 3; = ', iflagfas_avg
            stop
          end if     
        endif
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeAverageResponse(accum_psa,psa,nfreq,nsims,
     :                                  iflagpsa_avg)
c     calculate average response spectrum   ! revised by dmb
      dimension accum_psa(*),psa(*)
      do i=1,nfreq
          if (iflagpsa_avg .eq. 1) then ! arithmetic
                accum_psa(i)=
     :              accum_psa(i) + psa(i)
          else if (iflagpsa_avg .eq. 2) then ! geometric
                accum_psa(i)=
     :              accum_psa(i) + alog10(psa(i))
          else if (iflagpsa_avg .eq. 3) then ! rms
                accum_psa(i)=
     :              accum_psa(i) + psa(i)**2.0
          else
            print *,
     :         ' ERROR, iflagpsa_avg not 1, 2, or 3; = ', iflagpsa_avg
            stop
          end if     
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeFACCN(totalWave,npts,dt,nfreq,freq, FA)
c
c       Compute Fourier acceleration spectrum for final timeseries.
      dimension totalWave(*), freq(500), FA(500)
      complex fft(:)
      real spectrum(:)
      allocatable :: fft, spectrum
      
      allocate(fft(npts), spectrum(npts))
      
      fft = 0.0

c     Copy totalWave into new array because we need to save original for later.
      do j=1, npts
          spectrum(j)=totalWave(j)
      enddo

      call makeItComplex(spectrum,fft,npts)

c     Take padded time series into freq. domain
      call fork(npts,fft,-1.)
      ncall = npts/2 +1
      df = 1./(dt * float(npts))

c     Scale spectrum and put value back into spectrum array.
!      fft(ncall) = cmplx(real(fft(ncall)),0.)
      fft(ncall) = 0.0  ! ncall = Nyquist, set value to 0.0
      
      do jf = 1, ncall
        spectrum(jf) = dt* sqrt(float(npts)) * cabs(fft(jf))
      enddo

c     Now sample the spectrum into the same freq. bins used for PSA output.
      call sample(spectrum, ncall, df, nfreq, freq, FA)

c     Array FA now contains the nfreq values of sampled faccn.
c     TotalWave input array is unchanged.

      deallocate(fft, spectrum)
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeResponse(a,dt,dur,damp, Nfreq,freq,psa)

c     calculates response spectrum from acceleration time history
* Dates: 11/23/08 - Some changes by DMB (compute pi, use rdcalcdp)
*        12/06/09 - Changed "beta" to "damp_fraction" ("beta" is used elsewhere for
*                   shear-wave velocity)

      dimension Freq(*),psa(*),a(*)
      real pi,freq1,freq2,damp,Maxtim
c     damp is % critical damping. freq2 is maximum freq. to consider.

      pi = 4.0*atan(1.0)
c     Convert % damping to fraction
      damp_fraction=damp*0.01

c     Add 5 sec. to dur for shake-down
      maxtim=dur+5.0
      n=ifix(dur/dt)
      ntime=ifix(maxtim/dt)
      do j=n+1,ntime
          a(j)=0.
      enddo

c      Call new response spectrum routine for each desired freq
      do k=1,nfreq
          omega=2.*pi*freq(k)
          call rdcalcdp(a, ntime, omega, damp_fraction, dt, sd)
          PSA(k)=sd*omega*omega
      enddo
      return
      end
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeStochasticWave(
     :             seism,
     :             subFaultMoment, subfaultDistance,
     :             F0main, f0, risetime, dur_sub,
     :             nl,nw,nw_curr, 
     :             iflagscalefactor, scalingFactor,
     :             afreq1,afreq2,afreq3,amp1,amp2,amp3,sigma)
c     calculates accelerogram from individual subfault by using
c     band-limited white noise method (Boore, BSSA, no.6/83).
c     An w-squared spectrum is assumed

* Dates: 08/20/08 - Added scale factor flag so can choose which factor
*                   to use without recompiling the program
*                   (iflagscalefactor = 1, 2, 3, 4 = vel^2, acc^2, acc, DMB respectively).
*        08/21/08 - Add scalingfactor (h) to output.
*                   Do/do not apply the scalefactor taper, depending on the
*                   value of taper_scalefactor.
*                   Doubled dimension of s
*        12/01/08 - Removed "NumberOfActiveSubs" from argument lists
*        12/17/08 - Added f0 and risetime to argument list and removed iii, jjj, z.
*                   Also, the program used dur rather than dur_sub in ndur, etc.
*                   I changed this.  I also pass dur_sub through the argument list
*                   rather than compute it here.

      dimension seism(*)
      dimension afreq1(*),afreq2(*),afreq3(*),amp1(*),amp2(*),amp3(*)
      dimension sigma(*)
      complex s(:)
      allocatable :: s
      common /params/
     :   rho,beta,iKapa,fmax,
     :   Qmin,Q0,Qeta,                                                          !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
!     :   fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q,                             !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :   iwind,eps,eta,                                                         !Modified/added by KA, Mar.2012, eps/eta added
     :   n_crustal_amps,n_site_amps,n_empirical_amps,totalMoment,               !Modified/added by KA, Aug.2010
     :   c4taper, fc4taper,
     :   nsubs, flocut, nslope, pi, twopi
      common /par/ iFFTsub, dt, npadl, npadt
      common /seed/ iseed  ! Dave added this statement
      common /fltinf/FaultStrike,FaultDip,h,FaultLength,FaultWidth,dl,dw        !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
     :       ,AB14_flag                                                          !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      logical AB14_flag                                                         !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      
c     Specify subfault depth.                                                   !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      SubfaultDep=h+real(nw_curr-0.5)*sin(FaultDip/57.2957795)                  !Added by KA, Jan.2014 to handle frequency dependent geom_sprd

c     specify point-source constants
        tmax=float(iFFTsub)*dt
        df=1.0/tmax
        taper=0.02
        prtitn=1.0/sqrt(2.0)
        rtp=0.55
        fs=2.0
 
c     define spectral constant, for subfaultDistance=1 km
        const=prtitn*rtp*fs*(1.e-20)/(4.*pi*rho*beta**3.)
        nnyq=iFFTsub/2+1
 
      if (dur_sub .lt. 0.0) then
        print *,'STOPPING!!! Negative duration. Check duration model'
        stop
      end if
      
      ndur=dur_sub/dt
      ntaper=int(taper*ndur)
      nstop=ndur+ntaper+ntaper
      if(nstop.gt.iFFTsub) then
          write (*,*) 'PROGRAM STOPS: subsource duration too long.'
          write (*,*) 'Options to fix the problem:'
          write (*,*) '1) Increase iFFTsub or dt, or'
          write (*,*) '2) Reduce distance to observation point, or'
          write (*,*) '3) Change parameters of duration model'
          stop
      endif

      allocate (s(iFFTsub))                                                     ! iFFTsub is passed through the common block "par"

      do i=1,iFFTsub
          seism(i)=0.0
          s(i)=cmplx(0.0,0.0)
          seism(i)=0.0
      enddo

c generate Gaussian white noise with zero mean, unit variance. Noise is
c tapered.

      do i=1,nstop
          if (iwind .eq. 0) then
            call window(i,1,nstop,ntaper,wind)
          else if (iwind .eq. 1) then
            call wind2 (i,1,nstop,ntaper,dur_sub,eps,eta,wind)                  !Modified/added by KA, Mar.2012, eps/eta added
!            call wind2 (i,1,nstop,ntaper,dur_sub,0.2,0.2,wind)                  !Modified/added by KA, Mar.2012, eps/eta added
          else
            print *,'STOPPING: invalid value of iwind (= ', iwind, ')'
            stop
          end if
          s(i+npadl)=wind*cmplx(ggnqf(iseed),0.)
      enddo

c transform to frequency domain
      call fork (iFFTsub,s,-1.)

c find scaling factor, H
      averageMoment=totalMoment/float(nsubs)
      if (iflagscalefactor .eq. 3) then
        ScalingFactor=sqrt(float(nsubs))*(F0main/f0)**2 !! David Boore
      else      
        sum1=0
        sum2=0
        do i=1,nnyq
          f=(i-1)*df
          if (iKapa.eq.0) highc=1./sqrt(1.+(f/fmax)**8.)
          if (iKapa.eq.1) highc=exp(-pi*fmax*f)
          if (iflagscalefactor .eq. 1) then  ! vel spectrum
            s1=const*totalMoment*
     :           (2.*pi*f)**1.0*(1/(1.+(f/F0main)**2.))*
     :           highc
            s2=const*averageMoment*
     :           (2.*pi*f)**1.0*(1/(1.+(f/f0)**2.))*
     :           highc
          else if (iflagscalefactor .eq. 2) then  ! acc spectrum
            s1=const*totalMoment*
     :           (2.*pi*f)**2.0*(1/(1.+(f/F0main)**2.))*
     :           highc
            s2=const*averageMoment*
     :           (2.*pi*f)**2.0*(1/(1.+(f/f0)**2.))*
     :           highc
          else
            print *,' ERROR, iflagscalefactor = ', iflagscalefactor
            print *, '; not a legal value; QUITTING!!!'
            stop
          end if
          sum1=sum1+s1**2/float(nsubs)
          sum2=sum2+s2**2
        enddo      
        scalingFactor=sqrt(sum1/sum2)
      end if

* Compute parameters used to taper spectrum toward low frequencies
      c4taper = sqrt(float(nsubs))/scalingFactor
      fc4taper = f0/sqrt(c4taper)
      
c multiply spectrum by an w-squared spectral shape after normalizing to
c square of unit spectral amplitude
      call avg (nnyq,s,avamp2)
      avamp=sqrt(avamp2)
!      write(*,*)'Dist.=',subfaultDistance,',Depth=',SubfaultDep,
!     :AB14_flag
      do i=1,nnyq
        f=(i-1)*df
        if(AB14_flag)then
        s(i)=sourceSpectra(f,subfaultDistance,const,                            !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
     :          scalingFactor,subFaultMoment,f0,
     :          afreq1,afreq2,afreq3,amp1,amp2,amp3,sigma) *             !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
     :          s(i)*AB14_fac(f,subfaultDistance,SubfaultDep)/avamp !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
c     :                               s(i)/avamp                                !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
        else
        s(i)=sourceSpectra(f,subfaultDistance,const,                            !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
     :          scalingFactor,subFaultMoment,f0,
     :          afreq1,afreq2,afreq3,amp1,amp2,amp3,sigma) *             !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
c     :           s(i)*AB14_fac(f,subfaultDistance,SubfaultDep)/avamp            !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
     :                               s(i)/avamp                                !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
        endif
!        if (i.lt.2) print *, "pl: ", i, s(i)

        if (i.ne.1) s(iFFTsub+2-i)=conjg(s(i))
      enddo
      s(nnyq)=cmplx(real(s(nnyq)),0.)

!      print *, "p0: ", s(1)

c     transform back to time domain
      call fork(iFFTsub,s,1.)
      afact=sqrt(float(iFFTsub))/tmax
      do i=1,iFFTsub
          seism(i)=afact*real(s(i))
      enddo

      deallocate (s)


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine createRandomWeight(nl,nw,slipWeights)
      common /seed/ iseed

      dimension slipWeights(200,200)
           do i=1,nl
              do j=1,nw
c**              draw random number 0 to 1.  This was changed from
c                using =ggnqf(iseed)+1. Also impose nonzero weight.
                 slipWeights(i,j)=Ran1(iseed)
                 if (slipWeights(i,j).le.0.) slipWeights(i,j)=0.001
              enddo
           enddo
       return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dur_path_cmp(subfaultDistance, 
     :                        rpathdur, pathdur, durslope, ndur_hinges,
     :                        dur_path)
c     computes duration dur due to path for a given subsource
* Dates: 12/05/09 - duration parameters contained in arrays and 
* duration computation done as in smsim (subroutine durpath in rv_td_subs)
      real rpathdur(*), pathdur(*)  
      
      common /params/
     :   rho,beta,iKapa,fmax,
     :   Qmin,Q0,Qeta,                                                          !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
!     :   fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q,                             !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :   iwind,eps,eta,                                                         !Modified/added by KA, Mar.2012, eps/eta added
     :   n_crustal_amps,n_site_amps,n_empirical_amps,totalMoment,               !Modified/added by KA, Aug.2010
     :   c4taper, fc4taper,
     :   nsubs, flocut, nslope, pi, twopi
     
      if ( subfaultDistance .le. rpathdur(1) ) then
        dur_path = pathdur(1)
      else if ( subfaultDistance .ge. rpathdur(ndur_hinges) ) then
        dur_path = pathdur(ndur_hinges) + 
     :     (subfaultDistance - rpathdur(ndur_hinges))*durslope
     :              
      else
        call locate(rpathdur, ndur_hinges, subfaultDistance, j)
        dur_path = pathdur(j) + 
     :     (subfaultDistance - rpathdur(j))*( pathdur(j+1)-pathdur(j))
     :              /(rpathdur(j+1)-rpathdur(j))
      end if

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine findDistanceAndAzimuth (FaultLat,FaultLon,SiteLat,
     *                                 SiteLon,epi,azi,isitecoordflag)
c       calculates distance and azimuth between two
c       points using their latitudes and longitudes

* Dates: 02/17/09 - renamed some variables to avoid possible confusion
*                   between angles in degrees and radians.

* The input and output angles are in degrees

        pi = 4.0*atan(1.0)
        d2r = pi/180.0
        
        if (isitecoordflag .eq. 1) then  ! lat,long

          re=6371.
          alat=(FaultLat+SiteLat)/2.
          alat_radians=d2r*alat
          dlat=SiteLat-FaultLat  ! degrees
          dlon=SiteLon-FaultLon  ! degrees
          r1 = d2r*re*(SiteLat-FaultLat)
          epi= d2r*re*sqrt(dlat**2.+cos(alat_radians)**2.*dlon**2.)
          if(r1.ge.epi)then
             azi_radians=pi
             if(dlat.ge.0)azi_radians=0.
          else
             azi_radians=acos(r1/epi)
          endif

          if (dlon.le.0.) azi_radians=2.0*pi-azi_radians
          azi= azi_radians/d2r

        else if (isitecoordflag .eq. 2) then  ! R,Az
        
          epi = sitelat
          azi = sitelon
          
        else  ! assume cartesian coordinates
        
          xn = sitelat
          xe = sitelon
          
          epi = sqrt(xn**2 + xe**2)
          
          azi = 180.0*atan2(xe, xn)/pi
          if (azi .lt. 0.0) then
            azi = 360.0 + azi
          end if
          
        end if
          
          

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function findPeakValue(nptsTotalWave,totalWave)
c     find absolute maximum acceleration in simulated time history
      dimension totalWave(*)
        peakValue=0.
        do i=1,nptsTotalWave
          if(abs(totalWave(i)).gt.peakValue) peakValue=abs(totalWave(i))
        enddo
        findPeakValue=peakValue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real function findSubfaultDistance (R,h,FaultStrike,
     *                  fi2,FaultDip,dl,dw,i,j)
c       computes distance subfaultDistance from center of subfault (i,j)
c       to observation point using formula (1) (Figure 1)

* Note confusion in what dip means in the figure 1 of Beresnev and Atkinson.
* What they label as "delta" is what was used in the original version of
* this subroutine as "dip".  Their figure has delta_1 as the true 
* fault dip, probably added as a result of a review. I have
* replaced "dip" by 90_faultdip here to eliminate confusion (D. Boore, 17Feb09).

        pi = 4.0*atan(1.0)
        d2r = pi/180.0
        
        a90_faultdip_radians = d2r*(90.0-FaultDip)
        phi2_strike_radians = d2r*(fi2-FaultStrike)
        t1=R*cos(phi2_strike_radians)-(2.*i-1)*dl/2.
        t2=R*sin(phi2_strike_radians)-(2.*j-1)*dw/2.*
     :       sin(a90_faultdip_radians)
        t3=-h-(2.*j-1)*dw/2.*cos(a90_faultdip_radians)
        findSubfaultDistance=sqrt(t1**2.+t2**2.+t3**2.)
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine findSubFaultMoment(nl,nw,slipWeights,totalMoment,
     *                             weightedMoment)
      dimension slipWeights(200,200)
      dimension weightedMoment(200,200)

c       calculate moments of subfaults. First calculate the total of all weights
c       (totalSlipWeights), then moment per unit weight (elem), then individual moments and
c       put them back into array slipWeights(i,j)
c       Seismic Moment = stress para* (L**3)

      totalSlipWeights=0.
      do i=1,nl
        do j=1,nw
          totalSlipWeights=totalSlipWeights+slipWeights(i,j)
        enddo
      enddo

      elem=totalMoment/totalSlipWeights
      do i=1,nl
        do j=1,nw
          weightedMoment(i,j)=slipWeights(i,j)*elem
        enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getInputParameters(nu_ctl,
     :    FaultStrike,FaultDip,h,
     :    FaultLat,FaultLon, move_site, write_site_files,
     :    siteLocation,numberOfSites,isitecoordflag,
     :    fault_type,
     :    FaultLength, FaultWidth, 
     :    dl, dw, nl, nw, nsubs, 
     :    specify_length, specify_width, stress_ref, 
     :    vrup_beta,
     :    hyp_loc_fl, hyp_loc_fw, i0_in,j0_in,n_hypocenters,
     :    tpadl, tpadt, dt, beta, rho, amag,
     :    stress,iKapa,fmax,
!     :    fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q,                            !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :    Qmin,Q0,Qeta,                                                         !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :    r_ref, nsprd_segs, rlow, a_s, b_s, m_s,
     :    rdur, dur, durslope, ndur_hinges, 
     :    iwind,eps,eta,                                                        !Modified/added by KA, Mar.2012, eps/eta added
     :    flocut, nslope,
     :    nfreq,freq1,freq2,
     :    nfout, freq_out,
     :    f_stem,
     :    f_crustal_amps, 
     :    f_site_amps,
     :    f_empirical_amps,                                                     !Modified/added by KA, Aug.2010
     :    iseed, nsims,damp,
     :    islipweight, slipWeights,
     :    iRow,jColumn,pulsingPercent,
     :    iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak,
     :    iDynamicFlag, 
     :    iflagscalefactor, 
     :    iflagfas_avg, iflagpsa_avg,
     :    i_rise_time,
     :    AB14_flag)                                                            !Added by KA, Jan.2014 to handle frequency dependent geom_sprd

* Dates: 11/30/08 - Move computation of random slipweights 
*                   if Islipweight .eq. 1.0 from main program to here
*                   (a more logical place).
*        12/02/08 - Get i_rise_time, to determine what type of risetime is used 
     
      real rlow(*), a_s(*), b_s(*), m_s(*), freq_out(*)
      real rdur(*), dur(*)
      dimension slipWeights(200,200), siteLocation(3000,2)
      character f_stem*(*), 
     :  f_crustal_amps*30,f_site_amps*30,f_empirical_amps*30,                   !Modified/added by KA, Aug.2010
     :  aline*60, fault_type*(*),version_ctl*8, version_in*30
      character cmnts2skip(50)*80, buf_in*10, f_slip_weights*80
      
      logical f_exist, specify_length, specify_width,
     :        move_site, write_site_files, AB14_flag
      
      pi = 4.0*atan(1.0)
      d2r = pi/180.0

     
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      aline = ' '
      read (nu_ctl,'(a)')aline
      call trim_c(aline, nc_aline)
!      print *, ' Title: '//aline(1:nc_aline)
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      buf_in = ' '
      read (nu_ctl,'(a)') buf_in
      call upstr(buf_in)
      call trim_c(buf_in, nc_buf_in)
      if (buf_in(1:1) .eq. 'Y') then
         write_site_files = .true.
      else
         write_site_files = .false.
      end if
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) amag,stress,iKapa,fmax
!      if( iKapa.eq.1)then
!        print *,' mag,stress,kappa = ', amag,stress,fmax
!      else
!        print *,' mag,stress,fmax = ', amag,stress,fmax
!      endif
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) FaultLat, FaultLon
!      print *,' FaultLat,FaultLon = ', FaultLat, FaultLon
     
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) FaultStrike, FaultDip, h
!      print *,' FaultStrike,FaultDip(degree),Depth(km) = ',
!     :          FaultStrike,FaultDip,h
     
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      fault_type = ' '
      read(nu_ctl,'(a)') fault_type
      call trim_c(fault_type, nc_fault_type)
      call upstr(fault_type(1:1))
     
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) FaultLength, FaultWidth, dl, dw, stress_ref
      
      stress_factor = (stress_ref/stress)**(1.0/3.0)
 
      specify_length = .true.
      specify_width  = .true.


      if (faultlength .eq. 0.0) then
        specify_length = .false.
        if(fault_type(1:1) .eq. 'S') then
          faultlength = 10.0**(-2.57+0.62*amag)
        else if(fault_type(1:1) .eq. 'R') then
          faultlength = 10.0**(-2.42+0.58*amag)
        else if(fault_type(1:1) .eq. 'N') then
          faultlength = 10.0**(-1.88+0.50*amag)
        else
          faultlength = 10.0**(-2.44+0.59*amag)
        end if
        
!        print *, ' WC FL = ', faultlength
        faultlength = stress_factor * faultlength
!        print *, ' WC FL, after scaling = ', faultlength
       
      end if
      if (faultwidth .eq. 0.0) then
        specify_width  = .false.
        if(fault_type(1:1) .eq. 'S') then
          faultwidth = 10.0**(-0.76+0.27*amag)
        else if(fault_type(1:1) .eq. 'R') then
          faultwidth = 10.0**(-1.61+0.41*amag)
        else if(fault_type(1:1) .eq. 'N') then
          faultwidth = 10.0**(-1.14+0.35*amag)
        else
          faultwidth = 10.0**(-1.01+0.32*amag)
        end if
        
!        print *, ' WC FW = ', faultwidth
        faultwidth = stress_factor * faultwidth
!        print *, ' WC FW, after scaling = ', faultwidth
       
      end if
      
      nl=FaultLength/dl
      if (nl .lt. 1) nl = 1
      nw=FaultWidth/dw
      if (nw .lt. 1) nw = 1
 
      nsubs = nl * nw
        
      dl = FaultLength/float(nl)                                                ! Need to reset dl and dw because nl and nw are integers      
      dw = FaultWidth/float(nw)        

!      write(*,'(a,4f9.3)')'FaultLength,FaultWidth (km)           :'
!     :                          , FaultLength,FaultWidth
!      write(*,'(a,2i8)') 
!     :   'No. of subfaults along length and width of fault plane:', 
!     :    nl,nw
      if(nl.gt.200.or.nw.gt.200) then
        print *,"Burp!, you exceeded the
     :         the array dimentions, (200,200)"
        stop
      end if
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) vrup_beta
      

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) hyp_loc_fl, hyp_loc_fw, n_hypocenters
      if (nl .eq. 1 .and. nw .eq. 1) then   ! a single subsource
        i0_in = 1
        j0_in = 1
        n_hypocenters = 1
      else if (hyp_loc_fl .gt. 0.0 .and. hyp_loc_fw .ge. 0.0 ) then             ! nonrandom hypocenter 
                                                                                ! convert hyp_loc to subfault index
        i0_in = int(hyp_loc_fl/dl)+1
        j0_in = int(hyp_loc_fw/dw)+1
        n_hypocenters = 1                                                       ! not a random hypocenter, force this
      else                                                                      ! not a single subsource and random
        i0_in = 0
        j0_in = 0
      end if
      
        print *,' In GET_PARAMS: n_hypocenters, '//
     :          'hyp_loc_fl, hyp_loc_fw, dl, '//
     :           'dw, i0_in, j0_in ='
        print *, n_hypocenters, 
     :       hyp_loc_fl, hyp_loc_fw, dl, dw, i0_in, j0_in
     
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) i_rise_time

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) tpadl, tpadt, dt
!      write(*,'(a,3(1x,f9.3))')
!     :   'tpadl, tpadt, dt(sec)                     :',
!     :    tpadl, tpadt, dt
     
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) beta,rho
!      write(*,'(a,2f9.3)')"beta,rho                              :"
!     *  , beta,rho
     
* gsprd: 
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl, *) r_ref
      read(nu_ctl, *) nsprd_segs
      if(nsprd_segs.lt.0)then                                                   !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
          AB14_flag= .true.                                                     !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
          nsprd_segs=iabs(nsprd_segs)                                           !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      else                                                                      !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
          AB14_flag= .false.                                                    !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      endif                                                                     !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      do i = 1, nsprd_segs
        b_s(i)=0.                                                               !Modified/added by KA, Mar. 2012. Fix the value to zero
        m_s(i)=6.5                                                              !Modified/added by KA, Mar. 2012. Fix the value to 6.5
        read(nu_ctl, *) rlow(i), a_s(i)                                         !Modified/added by KA, Mar. 2010. remove parameter value setting by user
c       read(nu_ctl, *) rlow(i), a_s(i), b_s(i), m_s(i)                         !Original DMB input, two last terms' values hardwired in above lines
      end do
 
* Q:
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
!      read(nu_ctl, *) fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q                 !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
      read(nu_ctl, *) Qmin,Q0,Qeta                                              !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
      
* Duration:
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) ndur_hinges
      do i = 1, ndur_hinges
        read(nu_ctl,*) rdur(i), dur(i)
      end do
      read(nu_ctl,*) durslope
 
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) iwind, eps, eta
      if (iwind.eq.0) write (*,*)"Windowing info                    :
     * tapered boxcar"
      if (iwind.eq.1) write (*,*)"Windowing info                    :
     * Saragoni-Hart"

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) flocut, nslope

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) damp
      write(*,'(a, 1x,f8.3)')
!     : ' damping: ', 
!     :   damp

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) nfreq,freq1,freq2
!      write(*,'(a,i3,2f8.2)')"No.of freq,lowest and highest freq (Hz)for
!     * PSA :  ",  nfreq,freq1,freq2

      nfreq_max = 500

      if (iabs(nfreq).gt.nfreq_max) then
           print *,' '
           print *, ' nfreq = ', nfreq, 
     :              ' but it cannot exceed ', nfreq_max
           print *,' STOPPING!!!'
           stop
      endif

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)      
      read(nu_ctl,*) nfout
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) (freq_out(i),i=1,nfout)
	
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      f_stem = ' '
      read (nu_ctl,'(a)') f_stem
      call trim_c(f_stem, nc_f_stem)
      
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      f_exist = .false.
      f_crustal_amps = ' '
      read (nu_ctl,'(a)') f_crustal_amps
      call trim_c(f_crustal_amps,nc_f_crustal_amps)
      inquire(file=f_crustal_amps(1:nc_f_crustal_amps), exist=f_exist)
      if (.not. f_exist) then
        print *,' file '//f_crustal_amps(1:nc_f_crustal_amps)//
     :          ' does not exist; QUITTING!!!'
        stop
      end if

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      f_exist = .false.
      f_site_amps = ' '
      read (nu_ctl,'(a)') f_site_amps
      call trim_c(f_site_amps,nc_f_site_amps)
      inquire(file=f_site_amps(1:nc_f_site_amps), exist=f_exist)
      if (.not. f_exist) then
        print *,' file '//f_site_amps(1:nc_f_site_amps)//
     :          ' does not exist; QUITTING!!!'
        stop
      end if

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)                            !Modified/added by KA, Aug.2010
      f_exist = .false.                                                         !Modified/added by KA, Aug.2010
      f_empirical_amps = ' '                                                    !Modified/added by KA, Aug.2010
      read (nu_ctl,'(a)') f_empirical_amps                                      !Modified/added by KA, Aug.2010
      call trim_c(f_empirical_amps,nc_f_empirical_amps)                         !Modified/added by KA, Aug.2010
      inquire(file=f_empirical_amps(1:nc_f_empirical_amps)                      !Modified/added by KA, Aug.2010
     :          ,exist=f_exist)                                                 !Modified/added by KA, Aug.2010
      if (.not. f_exist) then                                                   !Modified/added by KA, Aug.2010
        print *,' file '//f_empirical_amps(1:nc_f_empirical_amps)//             !Modified/added by KA, Aug.2010
     :          ' does not exist; QUITTING!!!'                                  !Modified/added by KA, Aug.2010
        stop                                                                    !Modified/added by KA, Aug.2010
      end if                                                                    !Modified/added by KA, Aug.2010

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*)iDynamicFlag,pulsingPercent
!      write(*,'(a,i8,f8.0)')'Dynamic Flag,pulsing Percent: ',
!     :                           iDynamicFlag,pulsingPercent
     
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) iflagscalefactor                                           ! (1=vel^2; 2=acc^2; 3=asymptotic acc^2 (dmb))
      
* Get types of averages for FAS and PSA
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) iflagfas_avg
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) iflagpsa_avg

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak
      write(*,'(a,i2,4f8.3)')"Analytical flag and its Gama,Nu,T0,peak
     *             :", iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak
     
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) iseed, nsims 
!      write(*,'(a,1x,i6, 1x, i4)')
!     : 'Iseed, No. of trials for PSA calculation: ', 
!     :  iseed, nsims 
      iseed = -iabs(iseed)



      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*) islipweight
      write(*,'(a,i8)')"Slip flag 1=Random, 0=yours
     *       :",islipweight
     
      f_slip_weights = ' '
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,'(a)') f_slip_weights                                        ! need to read a dummy string
                                                                                ! even if slip weights are not
                                                                                ! specified
      call trim_c(f_slip_weights, nc_f_slip_weights)
      

      if (islipweight.eq.-1) then                                               ! assign unity to all cells
        do i = 1, nl
          do j = 1, nw
            slipWeights(i,j) = 1.0
          end do
        end do
      else if (islipweight.eq.0) then
        call get_lun(nu_slip_weights)
        open(nu_slip_weights,file=f_slip_weights(1:nc_f_slip_weights),
     :       status='unknown')
        read(nu_slip_weights,*)((slipWeights(i,j),i=1,nl),j=1,nw)
        close(nu_slip_weights)
        write(*,'(/a)')"       ***  read slip distribution   ***"
      else if (islipweight.eq.1) then
        call createRandomWeight(nl,nw,slipWeights)
      else if (islipweight.eq.2) then                                           ! Lines added by Karen to allow random slip generation in each trial 
        print *,' For every trial a new slip will be generated '                ! Lines added by Karen to allow random slip generation in each trial
      else
        print *,' islipweight = ', islipweight
        print *, 'NOT A LEGAL VALUE; STOPPING'
        stop
      endif

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read (nu_ctl,*)numberOfSites, isitecoordflag
!      write(*,*)"Number of Sites Around the Fault  :",
!     *                   numberOfSites
      print *, ' Type of coordinates (1=lat,long; 2=R,Az; 3=N,E):', 
     :           isitecoordflag
      
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      buf_in = ' '
      read (nu_ctl,'(a)') buf_in
      call upstr(buf_in)
      call trim_c(buf_in, nc_buf_in)
      if (buf_in(1:1) .eq. 'Y') then
         move_site = .true.
      else
         move_site = .false.
      end if
      
      if (move_site .and. FaultStrike .ne. 0.0) then
        print *,' ERROR: Cannot request move site if'//
     :          ' the fault strike .ne. 0.0; QUITTING'
        stop
      end if

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      do i=1,numberOfSites
            read (nu_ctl,*)siteLocation(i,1),siteLocation(i,2)
            print *,' SiteLat,SiteLon for Site #',
     :               i,' = ', siteLocation(i,1),siteLocation(i,2)
      enddo


      isitecoordflag_in = isitecoordflag
      Do i=1,numberOfSites
        if (move_site .and. isitecoordflag_in .eq. 1 ) then
          dum = 0.0
          ! eventually might do something with this case
        else if (move_site .and. isitecoordflag_in .eq. 2 ) then 
          isitecoordflag = 3
          r = siteLocation(i,1)
          az = siteLocation(i,2)
          xn = r*cos(d2r*az)
          xe = r*sin(d2r*az)
          siteLocation(i,1) = faultLength/2.0 + xn
          siteLocation(i,2) = xe
        else if (move_site .and. isitecoordflag_in .eq. 3) then
          if (sitelocation(i,1) .eq. 0.0) then  ! move to midpoint
            siteLocation(i,1) = faultLength/2.0
          end if
          if (siteLocation(i,2) .eq. 0.0) then  ! move to end
            siteLocation(i,1) = faultLength + siteLocation(i,1)
          end if
        end if
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ggnqf(iseed)
c     generates Gaussian white noise by approx. method of Ross (1985)
      usum=0.
      do j=1,12
          usum=usum+Ran1(iseed)
      enddo
      ggnqf=usum-6.
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine locateHypocentreRandomly(i0,j0,nl,nw)
      common /seed/ iseed

c**    draw random hypocentre
      integer i0,j0,nl,nw

      i0= nint(Ran1(iseed)*nl)
      if (i0 .lt. 1) i0=1
      if (i0 .gt. nl) i0=nl

      j0=nint(Ran1(iseed)*nw)
      if (j0 .lt. 1) j0=1
      if (j0 .gt. nw) j0=nw

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makeItComplex(x,cx,noOfSegPoints)
      dimension x(*)
      complex cx(*)
         do jjj = 1, noOfSegPoints
             cx(jjj) = cmplx(x(jjj),0.)
         enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine MavroPapa(totalWave,nptsTotalWave,
     *            iPapaFlag,PapaGama,PapaNu,PapaT0,amag,PapaPeak)
c     include Mavroeidis and Papageorgiou, 2003 approch
      common /par/ iFFTsub,dt, npadl, npadt
      dimension totalWave(*),PapaPulseW(70000),papaW(70000)
      complex   cTotalWave(70000),cPapaPulseW(70000),cPapaW(70000)
      complex   cTemp(70000)
      integer nptsTotalWave
      pi=4.0*atan(1.0)
      peakValue=findPeakValue(nptsTotalWave,totalWave)
      papaA=PapaPeak
      papaTp=10**(-2.9+0.5*aMag)
      papaFp=1.0/papaTp
      Fp=papaFp
      samplingRate=1./dt
      initialnptsTotalWave=nptsTotalWave

      T0=papaT0

cccccc step1 calculate pusle wave in the time domain
      PapaPulseW=0.0
      do i=1,nptsTotalWave
        t=(i-1)*dt
          if(t.ge.(T0-papaGama/2./Fp).and.t.le.(T0+papaGama/2./Fp))then
            c1=papaA*pi*Fp/papaGama
            c2=sin(2*pi*Fp*(t-T0)/papaGama)*cos(2*pi*Fp*(t-T0)+papaNu)
            c3=papaGama*sin(2*pi*Fp*(t-T0)+papaNu)
            c4=1+cos(2*pi*Fp*(t-T0)/papaGama)
            PapaPulseW(i)=c1*(c2+c3*c4)
        else
            PapaPulseW(i)=0.0
        endif
      enddo

cccccc step3-a take totalWave (stochastic) to the frequency domain
      call padding(TotalWave,nptsTotalWave)
      call makeItComplex(totalWave,cTotalWave,nptsTotalWave)
      call fork(nptsTotalWave,cTotalWave,-1.)


cccccc step3-b take PapaPulseW to the frequency domain
      call padding(PapaPulseW,nptsTotalWave)
      call makeItComplex(PapaPulseW,cPapaPulseW,nptsTotalWave)
      call fork(nptsTotalWave,cPapaPulseW,-1.)

cccccc step 4, keep the phase of stocastic,ccalculate abs(cTotalWave)-abs(cPapaPulseW)
cccccc calculate real and imaginary part  of the resulus
      do i=1,nptsTotalWave
        if(real(cTotalWave(i)).ge.0)then !!very importat when using atan
             phase=atan(aimag(cTotalWave(i))/real(cTotalWave(i)))
        else !!very important when using atan
             phase=atan(aimag(cTotalWave(i))/real(cTotalWave(i)))+pi
        endif
            amp=abs(abs(cTotalWave(i))-abs(cPapaPulseW(i)))
          x=amp*cos(phase)
          y=amp*sin(phase)
          cTemp(i)=cmplx(x,y)
      enddo
cccccc take care of complex conjugates
      cTemp(nptsTotalWave/2+1)=cmplx(1.,0.)*
     *                               cTemp(nptsTotalWave/2+1)
      do i=1, nptsTotalWave/2+1
                cTemp(nptsTotalWave+2-i) = conjg(cTemp(i))
      enddo
cccccc go back to the time domain
      call fork(nptsTotalWave,cTemp,+1.)
 
cccccc add them up
      do i=1,nptsTotalWave
         totalWave(i)=real(cTemp(i))+PapaPulseW(i)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function NumberOfPulsingSubs(i0,j0,i,j,nl,nw
     *                                ,NoeffectiveSubfaults)
c**   this function determines how many subfaults are active at this time (i,j).

* Dates: 12/01/08 - DMB added integer specification for Rmin, Rmax
*                   and used iabs rather than abs.

      integer Rmin, Rmax, r,NumberOfPulsingSubs
      n=0
      RMax=max(iabs(i-i0)+1,iabs(j-j0)+1)
      RMin=Rmax-NoeffectiveSubfaults
      if(RMin.lt.0)Rmin=0

      do ii=1,nl
         do jj=1,nw
            r=max(iabs(ii-i0)+1,iabs(jj-j0)+1)
            if(r.gt.RMin.and.r.le.RMax)n=n+1
         enddo
      enddo
      NumberOfPulsingSubs=n

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine padding(x,noOfSegPoints)
      dimension x(70000)

      do ii=1,100
        ncheck=2**(ii)
        if (ncheck.ge.noOfSegPoints) go to 110
        if(ncheck .gt. 70000) then
           write(*,*)' Burp !! Exceeds array dimensions.'
           stop
        endif
      enddo
110   continue

      do j = noOfSegPoints+1, ncheck
            x(j)=0.
      enddo
      noOfSegPoints=ncheck

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*  ------------------- BEGIN RAN1 --------------------------
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software $!6)$-"11j.
*  ------------------- END RAN1 --------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sample(x,npts,df,nfreq,freq,FA)
      dimension x(*), freq(500), FA(500)

c**   takes a large FFT array and samples for nfreq points,
c       in array FA.  Uses avg (sq) FA of 5 around each freq.

      FA=0.0
      
      do jf = 1, nfreq
      
        itarg= ifix(freq(jf)/df)+1
        
        f=df*float(itarg-1)        ! dmb replaced itarg with itarg-1
        
        if (itarg .lt. 4) then
        
            FA(jf)=-9.999   ! dmb
            
!            FA(jf)=9.999

        else
        
c             find spectrum around itarg
          sum = 0.0		          ! dmb
          do ii=itarg-2, itarg+2
              sum=0.2*x(ii)*x(ii) + sum   ! dmb
          enddo
          if(sum.eq.0) then
            FA(jf)= -9.999        ! dmb
          else
            FA(jf) = sqrt(sum)  ! dmb
          end if
        endif                              ! dmb
        
!        if(FA(jf).eq.0)FA(jf)=9.999
!            FA(jf) = alog10(sqrt(FA(jf)))
!        endif

      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine site_coords_in_degrees(
     :         FaultLat, FaultLon, 
     :         SiteLat, SiteLon, isitecoordflag, 
     :         site_lat_degrees, site_lon_degrees ) 
     
        pi = 4.0*atan(1.0)
        dtor = pi/180.0
     
! If site coords not in degrees, convert      
        if (isitecoordflag .eq. 1) then  ! degrees
          site_lat_degrees = sitelat
          site_lon_degrees = sitelon         
        else if (isitecoordflag .eq. 2) then  ! R,Az
          r = sitelat
          az = sitelon
          xn = r*cos(dtor*az)
          xe = r*sin(dtor*az)
          call km2deg_f( xn, xe, FaultLat, FaultLon,  
     :            site_lat_degrees, site_lon_degrees )          
        else  ! assume cartesian coordinates
          xn = sitelat
          xe = sitelon
          call km2deg_f( xn, xe, FaultLat, FaultLon,  
     :            site_lat_degrees, site_lon_degrees )          
        end if          
          
          
      
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function sourceSpectra (f,subfaultDistance,
     :     const,scalingFactor,aMoment,cornerF,
     :     afreq1,afreq2,afreq3,amp1,amp2,amp3,sigma)
c     Computes amplitude value of w-squared spectrum at
c     distance subfaultDistance.
c     Source terms correspond to subfaultDistance=1 km

* Dates: 02/07/09 - Determine gspread and q using gsprd_f and q_f (these
* need to be initialized and deallocated in the main program).

      dimension afreq1(*),afreq2(*),afreq3(*)
      dimension amp1(*),amp2(*),amp3(*),sigma(*)                        !Modified/added by KA, Aug.2010, and Feb. 2012
      common /params/
     :   rho,beta,iKapa,fmax,
     :   Qmin,Q0,Qeta,                                                          !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
!     :   fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q,                             !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :   iwind,eps,eta,                                                         !Modified/added by KA, Mar.2012, eps/eta added
     :   n_crustal_amps,n_site_amps,n_empirical_amps,totalMoment,               !Modified/added by KA, Aug.2010
     :   c4taper, fc4taper,
     :   nsubs, flocut, nslope, pi, twopi
      common /seed/ iseed

      w=twopi*f

      scalingFactor_hf_taper= scalingFactor * 
     :     c4taper*(1.0+(f/cornerF)**2.0)/(1.0+(f/fc4taper)**2.0)
 
      sa=const*aMoment*scalingFactor_hf_taper*w**2.0*
     :                   (1.0/(1.0+(f/cornerF)**2.0))
     
      if (f .eq. 0.0) gamma=0.
      if (f .ne. 0.0) then
          Q=q_f(f)
          gamma=w/(2.*Q*beta)
      endif
      anelas=exp(-gamma*subfaultDistance)
      if (iKapa.eq.0) highc=1./sqrt(1.+(f/fmax)**8.)
      if (iKapa.eq.1) highc=exp(-pi*fmax*f)
      
      ga = gsprd_f(subfaultdistance)
      
      am1=1.
      am2=1.
      am3=1.
      sig3=1.
      if (n_crustal_amps.ne.0) am1=amplf(f,afreq1,amp1,n_crustal_amps)
      if (n_site_amps.ne.0) am2=amplf(f,afreq2,amp2,n_site_amps)
      if (n_empirical_amps.ne.0)then                                            !Modified/added by KA, Aug.2010, Feb.2012
                  am3=amplf(f,afreq3,amp3,n_empirical_amps)                     !Modified/added by KA, Aug.2010 Feb.2012
                  sig3=amplf(f,afreq3,sigma,n_empirical_amps)                   !Modified/added by KA, Feb.2012
                  sig3=10.**(sig3*ggnqf(iseed))                                 !Modified/added by KA, Feb.2012
      endif                                                                     !Modified/added by KA, Feb.2012
      sourceSpectra=sa*ga*highc*anelas*am1*am2*am3*sig3*                        !Modified/added by KA, Aug.2010 Feb.2012
     :              buttrlcf( f, flocut, nslope)                                !Modified/added by KA, Aug.2010
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine window (i,nstart,nstop,ntaper,wind)
c     applies cosine tapered window
c     unit amplitude assumed
c     Subroutine from Dave Boore.
      wind=0.
      if (i.lt.nstart.or.i.gt.nstop) return
      wind=1.
      if(i.ge.nstart+ntaper.and.i.le.nstop-ntaper) return
      pi=3.141593
c     Value of wind goes from 0 to 1 from nstart to nstart+ntaper,
c     then from 1 to 0 from nstop-ntaper to nstop
      dum1=(nstop+nstart)/2.
      dum2=(nstop-nstart-ntaper)/2.
      wind=0.5*(1.-sin(pi*(abs(float(i)-dum1)-dum2)/float(ntaper)))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wind2 (i,nstart,nstop,ntaper,tw,eps,eta,wind)
c     applies Saragoni and Hart (1974) window, with parameters
c     tw (dur), eps (fraction of dur to reach peak), and
c     eta (fraction of peak ampl. at tw)
c     The Saragoni and Hart window is applied between
c     nstart + ntaper, and nstop - ntaper.  Outside these
c     bounds we do a quick cosine taper down to 0
c
      pi=3.141592654
      wind=0.
      if(i.lt.nstart.or.i.gt.nstop) return
      wind=1.

c     Apply Sargoni and Hart window

        b=-eps*alog(eta)/(1.+eps*(alog(eps)-1.))
        c=b/(eps*tw)
        a=(2.7182818/(eps*tw))**b
        if (i.ge.(nstart+ntaper).and.i.le.(nstop-ntaper)) then
           t=(tw/float(nstop-nstart-2*ntaper)) *
     *       (float(i-nstart-ntaper+1))
           wind=a*t**b*exp(-c*t)
           return
        else
c     Cos. taper goes from  0 to 1 from nstart to nstart+ntaper, from
c                           1 to 0 from nstop-ntaper to nstop.

        if (i.lt.(nstart+ntaper)) then
           t1=tw/float(nstop-nstart-2*ntaper)
           wf=a*t1**b*exp(-c*t1)
           wind=abs(sin((float(i-nstart)/float(ntaper))*pi/2.))*wf
        else
           wf=a*tw**b*exp(-c*tw)
           wind=abs(sin((float(nstop-i)/float(ntaper))*pi/2.))*wf
        endif
        return
      endif
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writeAcc(ioAcc, totalWave, npts, peakValue,
     :  isim, isite, HypoDistance, faultDistance, fpar,
     :  amag,siteLat,siteLon,depth,r,azimuth)                                  !Added by Karen for reporting more earthquake parameters in acceleration reports
c       writes acceleration time history into specified ascii file
      common /par/ iFFTsub,dt, npadl, npadt
      character fpar*(*)
      dimension totalWave(*)
      integer  npts

      write (ioAcc,"('**********************************************')")
      write (ioAcc,"('*******  Acceleration Time Series ************')")
      write (ioAcc,"('*** Site #',i4)") isite
      write (ioAcc,"('*** trial ',i4)") isim
      write (ioAcc,"('Input Parameters file =  ',a)") fPar

      write  (ioAcc,50) npts
50    format (i6,' samples')
      write  (ioAcc,30) dt
30    format ('dt: ',f6.3,' sec')
      write  (ioAcc,'(a)') 'data format: (1x,f8.3,1p,2x,e11.4)'
      write  (ioAcc,10) peakValue
10    format ('maximum absolute acceleration:',f7.2)
      amag=anint(amag*10)/10.                                                   !Added by Karen for reporting more earthquake parameters in acceleration reports
      write(ioAcc,'("Mag.          =",f8.2)')amag                               !Added by Karen for reporting more earthquake parameters in acceleration reports
      write(ioAcc,'("siteCoord1    =",f8.2)')siteLat                            !Added by Karen for reporting more earthquake parameters in acceleration reports
      write(ioAcc,'("siteCoord2    =",f8.2)')siteLon                            !Added by Karen for reporting more earthquake parameters in acceleration reports
      write(ioAcc,'("depth         =",f8.2)')depth                              !Added by Karen for reporting more earthquake parameters in acceleration reports
!      write(ioAcc,'("R             =",f8.2)')R                                  !Added by Karen for reporting more earthquake parameters in acceleration reports
!      write(ioAcc,'("azimuth       =",f8.2)')azimuth                            !Added by Karen for reporting more earthquake parameters in acceleration reports
      write(ioAcc,'("fault Dis.(km)=",f8.2)')faultDistance
      write(ioAcc,'("Hypo  Dis.(km)=",f8.2)')HypoDistance
      write  (ioAcc,'(2x,a, 1x,a)') 'time(s)', 'acc(cm/s**2)'

      do i=1,npts
         time=(i-1)*dt
         write (ioAcc,20) time,totalWave(i)
      enddo
20    format (1x,f8.3,1p,2x,e11.4)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writeHUS( ioHus, husid, npts, 
     :                 arias, dur_75_05, 
     :                 isim,isite,HypoDistance,faultDistance,fpar,
     :                 amag,siteLat,siteLon,depth,r,azimuth)                    !Added by Karen for reporting earthquake parameters in husid files, May.2012
*       writes husid time history, arias intensity, and 75%-5% duration
*       into specified ascii file
      common /par/ iFFTsub,dt, npadl, npadt
      character fpar*(*)
      dimension husid(*)
      integer  npts

      write (iohus,"('**********************************************')")
      write (iohus,"('*******  Husid Time Series ************')")
      write (iohus,"('*** Site #',i4)") isite
      write (iohus,"('*** trial ',i4)") isim
      write (iohus,"('Input Parameters file =  ',a)") fPar

      write  (iohus,50) npts
50    format (i6,' samples')
      write  (iohus,30) dt
30    format ('dt: ',f6.3,' sec')
      write  (iohus,'(a)') 'data format: (1x,f8.3,1p,2x,e11.4)'
      write  (iohus,'(a,1x,1pe10.3,a)') ' Arias intensity = ', arias            !Modified by Karen for making output more descriptive, May.2012
     1                                 ,' cm/sec'                               !Modified by Karen for making output more descriptive, May.2012
      write  (iohus,'(a,1x,1pe10.3,a)') ' 75%-5% duration = ', dur_75_05        !Modified by Karen for making output more descriptive, May.2012
     1                                 ,' sec'                                  !Modified by Karen for making output more descriptive, May.2012
      amag=anint(amag*10)/10.                                                   !Added by Karen for reporting earthquake parameters in husid files, May.2012
      write(iohus,'("Mag.             =",f8.2)')amag                            !Added by Karen for reporting earthquake parameters in husid files, May.2012
      write(iohus,'("siteCoord1       =",f8.2)')siteLat                         !Added by Karen for reporting earthquake parameters in husid files, May.2012
      write(iohus,'("siteCoord2       =",f8.2)')siteLon                         !Added by Karen for reporting earthquake parameters in husid files, May.2012
      write(iohus,'("depth            =",f8.2)')depth                           !Added by Karen for reporting earthquake parameters in husid files, May.2012
!      write(iohus,'("R                =",f8.2)')R                               !Added by Karen for reporting earthquake parameters in husid files, May.2012
!      write(iohus,'("azimuth          =",f8.2)')azimuth                         !Added by Karen for reporting earthquake parameters in husid files, May.2012
      write(iohus,'("fault Dis.(km)   =",f8.2)')faultDistance
      write(iohus,'("Hypo  Dis.(km)   =",f8.2)')HypoDistance
      write  (iohus,'(5x,a, 7x,a)') 'time', 'husid'

      do i=1,npts
         time=(i-1)*dt
         write (iohus,'(1x,f8.3,1p,1x,e11.4)') time,husid(i)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writePar(ioPar,
     :               FaultStrike,FaultDip,h,
     :               FaultLat, FaultLon, 
     :               siteLocation,numberOfSites, move_site,
     :               fault_type,
     :               FaultLength, FaultWidth,nl,nw,dl,dw, 
     :               specify_length, specify_width, stress_ref,
     :               vrup_beta,
     :            hyp_loc_fl, hyp_loc_fw, i0_in,j0_in,n_hypocenters,
     :               iseed, nsims,
     :               amag,
     :               stress, 
     :               pulsingPercent,iDynamicFlag,
     :               i_rise_time,
     :               iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak,
     :               iflagscalefactor,   
     :               iflagfas_avg, iflagpsa_avg,
     :               tpadl, tpadt,
     :               r_ref, nsprd_segs, rlow, a_s, b_s, m_s,
     :            rpathdur, pathdur, durslope, ndur_hinges,
     :               actualSlip
     :                             )

c     writes modeling parameters to specified ascii file

* Dates: 12/01/08 - DMB added writing of istart, istop, nptsenvelope, maxPointsOfWave
*        01/20/10 - Changed format descriptor of dt, kappa to f8.4; removed iFFTsub
*                   from output, because it has not been defined at the time that
*                   this subprogram is called.

      character f_crustal_amps*30,f_site_amps*30,f_empirical_amps*30            !Modified/added by KA, Aug.2010
      character typ*7
      common /par/ iFFTsub, dt, npadl, npadt
 
      logical specify_length, specify_width, move_site
      
      character fault_type*(*)
      
      dimension siteLocation(3000,2)
      real actualSlip(200,200)                                                  !Added by Karen for getting acc output in one file (20101110)

      real amag, r_ref, rlow(*), a_s(*), b_s(*), m_s(*)
      real rpathdur(*), pathdur(*)
      
      common /params/
     :   rho,beta,iKapa,fmax,
     :   Qmin,Q0,Qeta,                                                          !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
!     :   fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q,                             !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :   iwind,eps,eta,                                                         !Modified/added by KA, Mar.2012, eps/eta added
     :   n_crustal_amps,n_site_amps,n_empirical_amps,totalMoment,               !Modified/added by KA, Aug.2010
     :   c4taper, fc4taper,
     :   nsubs, flocut, nslope, pi, twopi
      common /fnames/ f_crustal_amps,f_site_amps,f_empirical_amps               !Modified/added by KA, Aug.2010
      
      write(ioPar,'("         modeling parameters     ")' )
      write(ioPar,'("Fault Strike              = ",f8.2)')FaultStrike
      write(ioPar,'("Fault dip                 = ",f8.2)')FaultDip
      write(ioPar,'("Fault depth to upper edge = ",f8.2)')h
      
      if(.not. specify_length) then
        write(ioPar,'(a,1x, f6.1)') 
     :    ' Fault length from Wells and Coppersmith for fault type '//
     :      fault_type(1:1)//', using a reference stress of ', 
     :                             stress_ref
      end if     
      write(ioPar,'("Fault Length              = ",f8.2)')FaultLength
      
      if(.not. specify_width) then
        write(ioPar,'(a,1x, f6.1)') 
     :    ' Fault width from Wells and Coppersmith for fault type '//
     :      fault_type(1:1)//', using a reference stress of ', 
     :                             stress_ref
      end if     
      write(ioPar,'("Fault Width               = ",f8.2)')FaultWidth
      
      write(ioPar,'("ratio of rupture to s-wave velocity = ",f8.2)') 
     :     vrup_beta
      
      write(ioPar,'("FaultLat                  = ",f8.2)')FaultLat
      write(ioPar,'("FaultLon                  = ",f8.2)')FaultLon

      write(ioPar,'("No.of subs along strike   = ",i8)')nl
      write(ioPar,'("No.of subs along dip      = ",i8)')nw
      write(ioPar,'("subfault length           = ",f8.2)')dl
      write(ioPar,'("subfault width            = ",f8.2)')dw
      write(ioPar,'("i_rise_time (1=orig,2=1/f0) = ",i2)')
     :                      i_rise_time
      write(ioPar,'("iseed, nsims = ",1x,i11, 1x,i4)')
     :                      iseed, nsims
      write(ioPar,'("-----------------------------------------------")')
      write(ioPar,'("input hypocenter at position    = ",
     :                                                 1p,2(1x,e10.3))')
     :                           hyp_loc_fl, hyp_loc_fw
      write(ioPar,'("input hypocenter at subfault    = ",2i4)')
     :                           i0_in,j0_in
      write(ioPar,'("n_hypocenters    = ",i4)') n_hypocenters
      write(ioPar,'("Mag.                      = ",f8.2)')amag
      write(ioPar,'("-----------------------------------------------")')
      write(ioPar,'("dt (sec)                  = ",f8.4)')dt
      write(ioPar,'("beta (km/s)               = ",f8.2)')beta
      write(ioPar,'("density (rho), gr/cm3     = ",f8.2)')rho
      write(ioPar,'("pulsing Percentage        = ",f8.2)')pulsingPercent
      write(ioPar,'(a,i1)') 'iflagscalefactor (1=vel^2; 2=acc^2;  '//
     :                     '3=asymptotic acc^2 (dmb)) = ',
     :               iflagscalefactor
      write(ioPar,'("flocut, nslope            = ",1x, f6.3, 1x, i2)') 
     :               flocut, nslope
      write(ioPar,'("iflagfas_avg              = ",i1)')iflagfas_avg
      write(ioPar,'("iflagpsa_avg              = ",i1)')iflagpsa_avg
      write(ioPar,'("stress parameter (bars)   = ",f8.2)')stress
        typ='fmax'
      if (iKapa.eq.1)then
        write(ioPar,'("kappa                     = ",f8.4)')fmax
      endif
      if (iKapa.eq.0)then
        write(ioPar,'("fmax                      = ",f8.2)')fmax
      endif
      write(ioPar,'("-----------------------------------------------")')
      write(ioPar,'("         Corner Frequency    ")' )
      if (iDynamicFlag.eq.1)then
          write(ioPar,'("Dynamic Corner Frequency Flag is  ON = ",i4)')
     *                          iDynamicFlag
      else
          write(ioPar,'("Dynamic Corner Frequency Flag is  OFF = ",i4)')
     *                          iDynamicFlag
          write(ioPar,'("Corner Frequency is static")')
      endif
      write(ioPar,'("-----------------------------------------------")')
      write(ioPar,'(a)')                                                        !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
     :             ' Qmin,          Q0,            eta,         = '             !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
      write(ioPar,*) Qmin, Q0, Qeta                                             !Modified/added by KA, Mar.2012. Return to Q0*F**eta definition
      write(ioPar,'(a)') 
     :  ' Path duration: ndur_hinges, rpathdur, pathdur, durslope: '
      write(ioPar,'(1x,i2)') ndur_hinges
      do i = 1, ndur_hinges
        write(ioPar,'(1x,f5.1, 1x,f6.2)') rpathdur(i), pathdur(i)
      end do
      write(ioPar,'(1x,f6.3)') durslope
      
!      write(ioPar,'(a)')                                                        !Modified/added by KA, Mar.2012. Remove magnitude dependant gsprd effect
!     :    ' gspread: i, nsprd_segs, r_ref, rlow, a_s, b_s, m_s'                 !Modified/added by KA, Mar.2012. Remove magnitude dependant gsprd effect
      write(ioPar,'(a)')                                                        !Modified/added by KA, Mar.2012. Remove magnitude dependant gsprd effect
     :    ' gspread: i, nsprd_segs, r_ref, rlow, a_s'                           !Modified/added by KA, Mar.2012. Remove magnitude dependant gsprd effect
      do i = 1, nsprd_segs
        write(ioPar,*) 
!     :    i, nsprd_segs, r_ref, rlow(i), a_s(i), b_s(i), m_s(i)                 !Modified/added by KA, Mar.2012. Remove magnitude dependant gsprd effect
     :    i, nsprd_segs, r_ref, rlow(i), a_s(i)                                 !Modified/added by KA, Mar.2012. Remove magnitude dependant gsprd effect
      end do

      write(ioPar,'("-----------------------------------------------")')
      if (iwind.eq.0) write (ioPar,120)
120   format ('window applied               = tapered boxcar')
      if (iwind.eq.1) then
      write (ioPar,130)
130   format ('window applied               = Saragoni-Hart')
      write (ioPar,131)eps                                                      !Modified/added by KA, Mar.2012, eps/eta added
131   format ('Saragoni-Hart epsilon        =      ',f8.4)                      !Modified/added by KA, Mar.2012, eps/eta added
      write (ioPar,132)eta                                                      !Modified/added by KA, Mar.2012, eps/eta added
132   format ('Saragoni-Hart eta            =      ',f8.4)                      !Modified/added by KA, Mar.2012, eps/eta added
      endif
      write(ioPar,'("-----------------------------------------------")')
      if (n_crustal_amps.eq.0.and.n_site_amps.eq.0.and.                         !Modified/added by KA, Aug.2010
     :    n_empirical_amps.eq.0)                                                !Modified/added by KA, Aug.2010
     :    write (ioPar,'(a)') '                no local amplification'
      if (n_crustal_amps.ne.0) write (ioPar,200) f_crustal_amps
      if (n_site_amps.ne.0) write (ioPar,200) f_site_amps
      if (n_empirical_amps.ne.0) write (ioPar,200) f_empirical_amps             !Modified/added by KA, Aug.2010
200   format ('amplification as in file ',a30)
      write(ioPar,'("-----------------------------------------------")')
  

      write(ioPar,'("-----------------------------------------------")')
      write(ioPar,'("         Analytical Pulse parameters     ")' )
      write(ioPar,'("-----------------------------------------------")')
      if (iPapaFlag.eq.1)then
          write(ioPar,'("Analytical Flag is  ON   = ",i4)')iPapaFlag
          write(ioPar,'("Analytical Gama          = ",f8.3)')PapaGama
          write(ioPar,'("Analytical Nu            = ",f8.3)')PapaNu
          write(ioPar,'("Analytical T0            = ",f8.3)')PapaT0
          write(ioPar,'("Analytical Peak          = ",f8.3)')PapaPeak
      else
          write(ioPar,'("Analytical Flag is  OFF   =  ",i4)')iPapaFlag
      endif

      write(ioPar,*)
      write(ioPar,'(a,3(1x,f9.3))')
     :   'tpadl, tpadt, dt(sec)                     :',
     :    tpadl, tpadt, dt
     
      do i = 1, numberOfSites
        write(ioPar,'(a, 1x,i3, 1x,a, 1x,f8.2, 1x, f8.2)')
     :   'For site', i, 'siteLocation coordinates 1&2 = ', 
     :                siteLocation(i,1),  siteLocation(i,2)
      end do
      if (move_site) then
        write(ioPar,'(a)')
     : 'Site may have been moved to midpoint or end of '//
     : 'surface projection of upper edge of fault' 
      end if

      write(ioPar,*)                                                            !Added by Karen for getting acc output in one file (20101110)
      write(ioPar,*)'The slip distribution is: '                                !Added by Karen for getting acc output in one file (20101110)
      do i=1,nw                                                                 !Added by Karen for getting acc output in one file (20101110)
        write(ioPar,'(200(f8.3,2x))')(actualSlip(j,i),j=1,nl)                   !Added by Karen for getting acc output in one file (20101110)
      enddo                                                                     !Added by Karen for getting acc output in one file (20101110)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writePSA_FA(ioPSA, freq,averagePGA,averagePGV,
     :  averagePSA,nfreq,nsims,nsimsflag,damp,amag,                             !Modified by Karen for getting proper spectral output of every simulation (20120430)
     :  siteLat,siteLon,depth,r,azimuth,faultDistance,averageFA,isite,
     :  d_jb, HypoDistance,fPar, nl, nw)                                        !Modified/Added by KA, Jan. 2013 to add Joyner-Boore distance
!     :  HypoDistance,fPar)
c       writes simulated average PSA spectrum at nfreq frequencies
c       into specified ascii file. First column is frequency in Hz,
c        second column is PSA value in cm/s**2

* Dates: 12/31/08 - Write sqrt(nl*nw)*psa, sqrt(nl*nw)*fas
*        01/06/09 - Write pgv, pga as first two values
*        01/13/09 - Write sqrt(nl*nw)*psa, sqrt(nl*nw)*fas and add SD, PSV.
*                 - Reverse order of FAS and PSA

      dimension freq(*),averagePSA(*),averageFA(*)
      character fname*20, fPar*(*)
      real omega, omega_sq
      
      pi=4.0*atan(1.0)
      twopi = 2.0 * pi
       
      sqrt_n = sqrt(float(nl*nw))

      write(ioPsa,"('**********************************************')")
      write(ioPsa,"('*******     Site #',i4)") isite
      if(nsimsflag==0)then                                                      !Modified by Karen for getting proper spectral output of every simulation (20120430)
          write(ioPsa,30) nsims                                                 !Modified by Karen for getting proper spectral output of every simulation (20120430)
      else                                                                      !Modified by Karen for getting proper spectral output of every simulation (20120430)
          write(ioPsa,40) nsims                                                 !Modified by Karen for getting proper spectral output of every simulation (20120430)
      endif                                                                     !Modified by Karen for getting proper spectral output of every simulation (20120430)
      write(ioPsa,50) int(damp)
      write(ioPsa,20) nfreq
      amag=anint(amag*10)/10.
      write(ioPsa,'("Mag.          =  ",f8.2)')amag
      write(ioPsa,'("siteCoord1    =  ",f8.2)')siteLat
      write(ioPsa,'("siteCoord2    =  ",f8.2)')siteLon
      write(ioPsa,'("depth         =  ",f8.2)')depth
!      write(ioPsa,'("R             =  ",f8.2)')R
!      write(ioPsa,'("azimuth       =  ",f8.2)')azimuth
      write(ioPsa,'("fault Dis.(km)= ",f8.2)')faultDistance
      write(ioPsa,'("Hypo  Dis.(km)= ",f8.2)')HypoDistance
      write(ioPsa,'("Rjb   Dis.(km)= ",f8.2)')d_jb                              !Modified by Karen for reporting Joyner-Boore distance (20130126)
      write(ioPsa,"('Input Parameters file =  ',a)") fPar

      write (ioPsa,'(1x,a)') 'data format: '//
     :         '(1x,f10.3, 1p, 1x,e11.4, 5x,e11.4,'//
     :         ' 2x,e11.4, 1x,e11.4)'
      write(ioPSA, '(3x,a, 5x,a,
     :                     1x,a, 1x,a,
     :                     1x,a, 3x,a)') 
     :       'freq(Hz)', 'per(s)',
     :       'FAS(cm/sec)', 'avgPSA(cm/s**2)', 
     :                    'avgPSV(cm/s)', 'avgSD(cm)'
 
      
      do i=1, nfreq
      
        omega = twopi*freq(i)
        omega_sq = omega**2.0
        per_dum=1.0/freq(i)
        if(averageFA(i).gt.0) then      ! dmb            
          write (ioPsa,'(1x,f10.3, 1x,f10.3,
     :               1p, 1x,e11.4, 5x,e11.4, 
     :               2x,e11.4, 1x,e11.4)') 
     :              freq(i), per_dum,
     :                averageFA(i), averagePSA(i),
     :                averagePSA(i)/omega, averagePSA(i)/omega_sq
     
        else
        per_dum=1.0/freq(i)
          write (ioPsa,'(1x,f10.3,  1x,f10.3,
     :               9x,a3,    1p, 5x,e11.4, 
     :               1x,e11.4, 1x,e11.4)') 
     :              freq(i), per_dum,
!     :              'NaN', averagePSA(i),                                       !Modified by Karen for not "NaN" amplitude reportings in spectral outputs, insteat "-1" is generated
     :              '-1.', averagePSA(i),                                       !Modified by Karen for not "NaN" amplitude reportings in spectral outputs, insteat "-1" is generated
     :              averagePSA(i)/omega, averagePSA(i)/omega_sq
        endif
      enddo

      write (ioPsa,*)
      write (ioPsa,*)
20    format (i4,' samples')
30    format('PSA and Fouries apectra of trial num. ',i5)
40    format('Average PSA and Fourier spectrum over ',i5,' trials')
50    format ('damping: ',i3,'%')
      return
      end


      subroutine check_directories()
      logical acc_exist, psa_exist, husid_exist                                              !Logical variables used for checking existence of directories/Added by KA, Feb.2012
      integer(4) i_system
!      inquire(directory='ACC',exist=acc_exist)
!      inquire(file='ACC/.',exist=acc_exist) ! gfortran version
!      if(.not. acc_exist)i_system=system('mkdir ACC')
!      inquire(directory='PSA',exist=psa_exist)
!      inquire(file='PSA/.',exist=psa_exist) ! gfortran version
!      if(.not. psa_exist)i_system=system('mkdir PSA')
!      inquire(directory='other',exist=husid_exist)
!      inquire(file='other/.',exist=husid_exist) ! gfortran version
!      if(.not. husid_exist)i_system=system('mkdir other')
      i_system=system('mkdir -p ACC')
      i_system=system('mkdir -p PSA')
      i_system=system('mkdir -p other')
      return
      end

* --------------------------- BEGIN TRIM_C -----------------------
      subroutine trim_c(cstr, nchar)

* strips leading and trailing blanks from cstr, returning the
* result in cstr, which is now nchar characters in length

* Strip off tabs also.

* Here is a sample use in constructing a column header, filled out with 
* periods:

** Read idtag:
*        idtag = ' '
*        read(nu_in, '(1x,a)') idtag
*        call trim_c(idtag, nc_id)
** Set up the column headings:
*        colhead = ' '
*        colhead = idtag(1:nc_id)//'......' ! nc_id + 6 > length of colhead

* Dates: 12/23/97 - written by D. Boore
*        12/08/00 - pad with trailing blanks.  Otherwise some comparisons
*                   of the trimmed character string with other strings
*                   can be in error because the original characters are left
*                   behind after shifting.  For example, here is a string
*                   before and after shifting, using the old version:
*                      col:12345
*                           MTWH  before
*                          MTWHH  after (but nc = 4).
*        03/21/01 - Check for a zero length input string
*        11/09/01 - Change check for zero length string to check for all blanks
*        10/19/09 - Strip off tabs

      character cstr*(*)

      if(cstr .eq. ' ') then
        nchar = 0
        return
      end if

      nend = len(cstr)

! Replace tabs with blanks:

      do i = 1, nend
        if(ichar(cstr(i:i)) .eq. 9) then
           cstr(i:i) = ' '
        end if
      end do



*      if(nend .eq. 0) then
*        nchar = 0
*        return
*      end if

      do i = nend, 1, -1
        if (cstr(i:i) .ne. ' ') then
           nchar2 = i
           goto 10
        end if
      end do

10    continue

      do j = 1, nchar2
        if (cstr(j:j) .ne. ' ') then
          nchar1 = j
          goto 20
        end if
      end do

20    continue
   
      nchar = nchar2 - nchar1 + 1
      cstr(1:nchar) = cstr(nchar1: nchar2)
      if (nchar .lt. nend) then
        do i = nchar+1, nend
          cstr(i:i) = ' '
        end do
      end if

      return
      end
* --------------------------- END TRIM_C -----------------------
* --------------------- BEGIN UPSTR ----------------------------------
      Subroutine UPSTR ( text )
* Converts character string in TEXT to uppercase
* Dates: 03/12/96 - Written by Larry Baker

C
      Implicit   None
C
      Character  text*(*)
C
      Integer    j
      Character  ch
C
      Do 1000 j = 1,LEN(text)
         ch = text(j:j)
         If ( LGE(ch,'a') .and. LLE(ch,'z') ) Then
            text(j:j) = CHAR ( ICHAR(ch) - ICHAR('a') + ICHAR('A') )
         End If
 1000    Continue
C
      Return
      End
* --------------------- END UPSTR ----------------------------------
* ---------------------- BEGIN SKIP -------------------
      subroutine SKIP(lunit, nlines)
        if (nlines .lt. 1) then
          return
        else
          do i = 1, nlines
             read(lunit, *)
          end do
          return
        end if
      end
* ---------------------- END SKIP -------------------
* ------------------------------------------------------------------ skipcmnt
      subroutine skipcmnt(nu, comment, ncomments)

* Skip text comments in the file attached to unit nu, but save skipped 
* comments in character array comment.  Skip at least one line, and more as 
* long as the lines are preceded by "|" or "!".

* Dates: 04/16/01 - Written by D. Boore
*        12/07/01 - Added check for eof
*        11/04/03 - Use trim_c to trim possible leading blank
*        02/03/07 - Initialize comments to blank

      character comment(*)*(*), buf*80

      ncomments = 0
100   buf = ' '
      read (nu,'(a)',end=999) buf
      call trim_c(buf,nc_buf)
      if (buf(1:1) .eq.'!' .or. buf(1:1) .eq.'|' .or. 
     :                     ncomments + 1 .eq. 1) then
        ncomments = ncomments + 1
        comment(ncomments) = ' '
        comment(ncomments) = buf(1:nc_buf)
        goto 100
      else 
        backspace nu
      end if

999   continue
 
      return
      end
* ------------------------------------------------------------------ skipcmnt

* --------------------------- BEGIN GET_LUN ----------------
      subroutine get_lun(lun)

* Finds a logical unit number not in use; returns
* -1 if it cannot find one.

* Dates -- 05/19/98 - Written by D. Boore, following
*                     Larry Baker's suggestion

      logical isopen
      do i = 99,10,-1
        inquire (unit=i, opened=isopen)
        if(.not.isopen) then
          lun = i
          return
        end if
      end do
      lun = -1

      return
      end
* --------------------------- END GET_LUN ----------------
     
*----------------- BEGIN AccSqInt -----------------------------
      subroutine accsqint(acc, npts, dt, rmv_trnd, a_sq_int)


* Form integral of acceleration squared, assuming that the acceleration
* is represented by straight lines connecting the digitized values.  This
* routine can be used to compute Arias intensity, defined as

*            Ixx = (pi/2g)*int(acc^2*dt), integrating from 0.0 to the total
*  duration of the record.  The units of Ixx are 
*  velocity [ l^(-1)t^2*(l/t^2)^2*t ] =  l^(-1+2)*t^(2-4+1) = l*t^(-1) = l/t

* Be sure to use consistent units for the acceleration of gravity (g) and acc.
* I am not sure what is conventionally used, but Wilson (USGS OFR 93-556) 
* uses m/sec.

* Dates: 01/13/99 - Written by D.M. Boore

      real a_sq_int(*), acc(*)
      logical rmv_trnd
      double precision cum, a1, a2, ddt_3

      if (rmv_trnd) then      
* remove trend first
        call dcdt(acc, dt, npts, 1, npts, .false., .true.)
      end if

* compute integral of squared acceleration (assume a_sq_int = 0 for first point)

      ddt_3 = dble(dt/3)

      cum = 0.0

      a_sq_int(1) = sngl(cum)
      do j=2,npts
        a1 = acc(j-1)
        a2 = acc(j)
        cum = cum + (a1**2+a1*a2+a2**2)*ddt_3
        a_sq_int(j) = sngl(cum)
      end do

* high pass filter the velocity (may want to allow this in a future version;
* as it is, the acceleration time series can be filtered, so there is no need
* to do it again).

      return
      end
*----------------- END AccSqInt -----------------------------
* ------------------------ begin dcdt -------------------
      subroutine dcdt (y,dt,npts,indx1,indx2,ldc,ldt)
c+
c  dcdt - fits dc or trend between indices indx1 and indx2.
c         then removes dc or detrends whole trace.
c         y is real, dt = delta t.
c         if remove dc, ldc = .true.
c         if detrend, ldt = .true.
c-

* Dates: 12/14/00 - Cleaned up formatting of original program

      real y(*)
      logical ldc,ldt

      if (.not. ldc .and. .not. ldt) then
        return
      end if

c
c...fit dc and trend between indices indx1 and indx2.
      nsum = indx2-indx1+1
      sumx = 0.0
      sumx2 = 0.0
      sumy = 0.0
      sumxy = 0.0
      do i=indx1,indx2
         xsubi = (i-1)*dt
         sumxy = sumxy+xsubi*y(i)
         sumx = sumx+xsubi
         sumx2 = sumx2+xsubi*xsubi
         sumy = sumy+y(i)
      end do
c
c... remove dc.
      if (ldc) then
        avy = sumy/nsum
        do i=1,npts
          y(i) = y(i)-avy
        end do
* Debug
        write(*,'(a)') ' indx1, indx2, avy'
        write(*, *)      indx1, indx2, avy
* Debug



        return
      endif
c
c... detrend. see draper and smith, p. 10.
      if (ldt) then
        bxy = (sumxy-sumx*sumy/nsum)/(sumx2-sumx*sumx/nsum)
        axy = (sumy-bxy*sumx)/nsum
        qxy = dt*bxy
        do i=1,npts
          y(i) = y(i)-(axy+(i-1)*qxy)
        end do
        return
      endif
c
      return
      end
* ------------------------ end dcdt -------------------
* --------------------- BEGIN LOCATE -----------------
      SUBROUTINE locate(xx,n,x,j)
      
* Comments added by D. Boore on 26feb2010:
*  finds j such that xx(j) < x <= xx(j+1)
*  EXCEPT if x = xx(1), then j = 1 (logically it would be 0 from
*  the above relation, but this is the same returned value of j
*  for a point out of range).
*  Also, if x = xx(n), j = n-1, which is OK
*  Note that j = 0 or j = n indicates that x is out of range.
*
* See the program test_locate.for to test this routine.

      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
* --------------------- END LOCATE -----------------
* --------------------- BEGIN LOCATE_D -----------------
      SUBROUTINE locate_d(xx,n,x,j)
      INTEGER j,n
      double precision x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
* --------------------- END LOCATE_D -----------------
* --------------------- BEGIN ZBRENT -----------------
      FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))then
         write(*,*) 'root must be bracketed for zbrent'
      endif
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      write (*,*) 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
* --------------------- END ZBRENT -----------------
* --------------------------- BEGIN DIST_3DF ---------------------
      subroutine dist_3df(alat_sta, along_sta, 
     :   alat_ref, along_ref, h_ref, strike_f, dip_f, w1, w2, s1, s2, 
     :   h_min_c, d_jb, az_jb, d_cd2f, az_cd2f, d_c, az_c,
     :   d_sta_n, d_sta_e, irgn_cd2f, irgn_c, irgn_jb)

* Computes various distance measures from a station to a fault.  The 
* orientation of the fault is that used by Spudich et al., Yucca Mt. project.
* The fault is assumed to be a rectangle whose upper and lower edges are
* horizontal.

* Input:
*      alat_sta, along_sta:  latitude and longitude of station, in degrees,
*                            west longitude is negative
*      alat_ref, along_ref:  as above, for reference point used in defining
*                            fault orientation
*      h_ref:                depth to reference point
*      strike_f, dip_f:      azimuth and dip of fault, in degrees.  strike_f is
*                            measured positive clockwise from north; dip_f is
*                            measured from the horizontal.  When looking in 
*                            the direction strike_f, a positive dip is down to
*                            the right.
*      w1, w2, s1, s2:       distances from the reference point to the edges of
*                            the fault.  s1 and s2 are the distances along
*                            strike to the near and far edges; w1 and w2 are
*                            the distances along the dip direction to the upper
*                            and lower edges of the fault.
*      h_min_c:              minimum depth for computing Campbell's distance
*                            (usually 3.0 km)

* Output:
*      d_jb, d_cd2f, d_c:    Joyner & Boore, closest distance to fault surface,
*                            and Campbell distance, respectively.
*      az_jb, az_cd2f, az_c: as above, for azimuths (NOT YET IMPLEMENTED IN THIS
*                            SUBROUTINE)
*      d_sta_n, d_sta_e:     north and east components of station location
*                            relative to the reference point
*      irgn_cd2f, etc:       region in fault-plane coordinates used to 
*                            compute distances.  I could include a sketch here,
*                            but I will not take the time.  These output 
*                            variables were included mainly to help me check
*                            the subroutine.


* Dates: 12/06/98 - Written by D. Boore
*        09/17/00 - Bring in subroutines via include statement;
*                   renamed faz and az_f to fstrike, strike_f
*        11/12/01 - Changed specification in headers above to indicate that
*                   west longitude is negative (consistent with revision to
*                   subroutine deg2km_f)

      real dist_sta(3)
      real ix(3), iy(3), iz(3)

      pi = 4.0*atan(1.0)
      dtor = pi/ 180.0

* set up unit vectors in fault coordinates in terms of north, east, down 
* unit vectors.

* Convert angles to radians:
      fstrike =  dtor * strike_f
      fdip = dtor * dip_f

* Initialize arrays:

      do i = 1, 3
         ix(i) = 0.0
         iy(i) = 0.0
         iz(i) = 0.0
      end do

* Compute unit vectors:        ! 1, 2, 3 correspond to n, e, d

      ix(1) = cos(fstrike)
      ix(2) = sin(fstrike)
      ix(3) = 0.0

      iy(1) = -sin(fstrike)*sin(fdip)
      iy(2) =  cos(fstrike)*sin(fdip)
      iy(3) =          -cos(fdip)

      iz(1) = -sin(fstrike)*cos(fdip)
      iz(2) =  cos(fstrike)*cos(fdip)
      iz(3) =           sin(fdip)

* Convert station lat, long into distance north and east:
      call deg2km_f(alat_sta, along_sta, alat_ref, along_ref,
     :              dist_sta(1), dist_sta(2))
      dist_sta(3) = -h_ref    ! note minus sign
      
      d_sta_n = dist_sta(1)
      d_sta_e = dist_sta(2)

* Convert coordinates of reference-to-station vector from n,e,d coordinates
* into fault coordinates:
 
      rx = 0.0
      ry = 0.0
      rz = 0.0
  
      do i = 1, 3
        rx = rx + dist_sta(i) * ix(i)
        ry = ry + dist_sta(i) * iy(i)
        rz = rz + dist_sta(i) * iz(i)
      end do

* Find region and closest distance to fault in the fault plane coordinates:

      call find_h(rx, rz, w1, w2, s1, s2, h_cd2f,
     :            irgn_cd2f)                                                    ! cd2f = Closest Distance
                                                                                !        to Fault

* Now do it for Campbell:

* Define w1 for Campbell (I assume that w2 does not need defining; in other
* words, not all of the fault plane is above the Campbell depth)

      d2top_c = h_min_c
      d2top = h_ref + w1 * iz(3)        ! iz(3) = sin(fdip)
      if ( d2top .lt. d2top_c .and. iz(3) .ne. 0.0) then
        w1_c = (d2top_c - h_ref)/ iz(3)
      else
        w1_c = w1
      end if

      call find_h(rx, rz, w1_c, w2, s1, s2, h_c,      ! c = Campbell
     :            irgn_c)


* Now do it for Joyner-Boore:

* Need to find rx, ry, rz, w1, w2, s1, s2 in terms of coordinates
* of the fault plane projected onto the surface:

      s1_jb = s1
      s2_jb = s2
      w1_jb = w1 * cos(fdip)
      w2_jb = w2 * cos(fdip)

      rx_jb = rx
      rz_jb = -sin(fstrike) * dist_sta(1) + cos(fstrike) * dist_sta(2)

* Then find the region and distance in the plane to the fault surface
      call find_h(rx_jb, rz_jb, 
     :            w1_jb, w2_jb, s1_jb, s2_jb, h_jb,
     :            irgn_jb)

* Now compute the distances:

      d_cd2f = sqrt(h_cd2f**2 + ry   **2)
      d_c    = sqrt(h_c   **2 + ry   **2)
      d_jb   = h_jb

* (Work on azimuths later)
      az_cd2f = 999.9
      az_c    = 999.9
      az_jb   = 999.9

      return

      end
* --------------------------- END DIST_3DF ---------------------

*-------------------- BEGIN FIND_H ----------------------------
      subroutine find_h(rx, rz, w1, w2, s1, s2, h, iregion)

* Now it is easy to see where the station lies with respect to the fault;
* there are 9 possibilities:  

      if (   rx .le. s1 .and. rz .le. w1 ) then
* region 1 (see notes)
        iregion = 1
        h = sqrt( (s1-rx)**2 + (w1-rz)**2 )
       
      else if (   rz .le. w1 .and. rx .ge. s1 .and. rx .le. s2 ) then
* region 2 (see notes)
        iregion = 2
        h = w1 - rz

      else if (   rx .ge. s2 .and. rz .le. w1 ) then
* region 3 (see notes)
        iregion = 3
        h = sqrt( (rx-s2)**2 + (w1-rz)**2 )
        
      else if (   rx .ge. s2 .and. rz .ge. w1 .and. rz .le. w2 ) then
* region 4 (see notes)
        iregion = 4
        h = rx - s2

      else if (   rx .ge. s2 .and. rz .ge. w2 ) then
* region 5 (see notes)
        iregion = 5
        h = sqrt( (rx-s2)**2 + (rz-w2)**2 )
        
      else if (   rz .ge. w2 .and. rx .ge. s1 .and. rx .le. s2 ) then
* region 6 (see notes)
        iregion = 6
        h = rz - w2

      else if (   rz .ge. w2 .and. rx .le. s1 ) then
* region 7 (see notes)
        iregion = 7
        h = sqrt( (s1-rx)**2 + (rz-w2)**2 )
        
      else if (   rx .le. s1 .and. rz .ge. w1 .and. rz .le. w2 ) then
* region 8 (see notes)
        iregion = 8
        h = s1 - rx

      else if (      rx .ge. s1 .and. rx .le. s2 
     :         .and. rz .ge. w1 .and. rz .le. w2 ) then
* region 9 (see notes)
        iregion = 9
        h = 0.0

      else
* reaching this is an error
        write(*,'(a)') ' ERROR: Region not found in find_h'

      end if

      return
      end
*-------------------- END FIND_H ----------------------------

!      include '\forprogs\deg2km_f.for'
!      include '\forprogs\locate.for'

*-------------------- BEGIN KM2DEG ----------------------------
      subroutine km2deg_f( vn_in, ve_in, alat_ref_in, along_ref_in, 
     :               vlat_out, vlong_out )
        
* convert km north and east from a reference point into lat, long

* assumes positive latitude between 0 and 70 degrees
* assumes east longitude is positive
* assumes angles in degrees

* WARNING: NEEDS DOUBLE PRECISION VERSION OF LOCATE (ATTACHED HERE)

* Dates:  10/01/95 - written by D. Boore
*         05/27/98 - Name changed to km2deg_f
*         06/01/98 - Changed to double precision
*         02/14/09 - Changed input to single precision
               
      double precision alat_tbl(71), b_tbl(71), adcoslat_tbl(71)
      double precision vn, ve, alat_ref, along_ref, vlat, vlong
      Data alat_tbl /
     : 0.000000, 1.000000, 2.000000, 3.000000, 4.000000, 5.000000,
     : 6.000000, 7.000000, 8.000000, 9.000000,10.000000,11.000000,
     :12.000000,13.000000,14.000000,15.000000,16.000000,17.000000,
     :18.000000,19.000000,20.000000,21.000000,22.000000,23.000000,
     :24.000000,25.000000,26.000000,27.000000,28.000000,29.000000,
     :30.000000,31.000000,32.000000,33.000000,34.000000,35.000000,
     :36.000000,37.000000,38.000000,39.000000,40.000000,41.000000,
     :42.000000,43.000000,44.000000,45.000000,46.000000,47.000000,
     :48.000000,49.000000,50.000000,51.000000,52.000000,53.000000,
     :54.000000,55.000000,56.000000,57.000000,58.000000,59.000000,
     :60.000000,61.000000,62.000000,63.000000,64.000000,65.000000,
     :66.000000,67.000000,68.000000,69.000000,70.000000
     :/
      Data b_tbl /
     : 1.842808, 1.842813, 1.842830, 1.842858, 1.842898, 1.842950,
     : 1.843011, 1.843085, 1.843170, 1.843265, 1.843372, 1.843488,
     : 1.843617, 1.843755, 1.843903, 1.844062, 1.844230, 1.844408,
     : 1.844595, 1.844792, 1.844998, 1.845213, 1.845437, 1.845668,
     : 1.845907, 1.846153, 1.846408, 1.846670, 1.846938, 1.847213,
     : 1.847495, 1.847781, 1.848073, 1.848372, 1.848673, 1.848980,
     : 1.849290, 1.849605, 1.849992, 1.850242, 1.850565, 1.850890,
     : 1.851217, 1.851543, 1.851873, 1.852202, 1.852531, 1.852860,
     : 1.853188, 1.853515, 1.853842, 1.854165, 1.854487, 1.854805,
     : 1.855122, 1.855433, 1.855742, 1.856045, 1.856345, 1.856640,
     : 1.856928, 1.857212, 1.857490, 1.857762, 1.858025, 1.858283,
     : 1.858533, 1.858775, 1.859008, 1.859235, 1.859452
     :/
      Data adcoslat_tbl /
     : 1.855365, 1.855369, 1.855374, 1.855383, 1.855396, 1.855414,
     : 1.855434, 1.855458, 1.855487, 1.855520, 1.855555, 1.855595,
     : 1.855638, 1.855683, 1.855733, 1.855786, 1.855842, 1.855902,
     : 1.855966, 1.856031, 1.856100, 1.856173, 1.856248, 1.856325,
     : 1.856404, 1.856488, 1.856573, 1.856661, 1.856750, 1.856843,
     : 1.856937, 1.857033, 1.857132, 1.857231, 1.857331, 1.857435,
     : 1.857538, 1.857643, 1.857750, 1.857858, 1.857964, 1.858074,
     : 1.858184, 1.858294, 1.858403, 1.858512, 1.858623, 1.858734,
     : 1.858842, 1.858951, 1.859061, 1.859170, 1.859276, 1.859384,
     : 1.859488, 1.859592, 1.859695, 1.859798, 1.859896, 1.859995,
     : 1.860094, 1.860187, 1.860279, 1.860369, 1.860459, 1.860544,
     : 1.860627, 1.860709, 1.860787, 1.860861, 1.860934
     :/

      pi = 4.0*atan(1.0)
      d2r = pi/ 180.0

      vn = dble(vn_in)
      ve = dble(ve_in)
      alat_ref =  dble(alat_ref_in)
      along_ref =  dble(along_ref_in) 

* interpolate to find proper arc distance:

      call locate_d( alat_tbl, 71, alat_ref, j)
      b = b_tbl(j) + (alat_ref-alat_tbl(j))*
     :  (b_tbl(j+1)-b_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      adcoslat = adcoslat_tbl(j) + (alat_ref-alat_tbl(j))*
     :  (adcoslat_tbl(j+1)-adcoslat_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      a = adcoslat * cos(d2r*alat_ref)

      dlambda = +ve/a ! version with minus used if assume west long is +
*      dlambda = -ve/a ! minus; positve ve corresponds to decrease in long
      dphi    =  vn/b

* convert from minutes of arc to degrees:
      dlambda = dlambda / 60.0
      dphi    = dphi    / 60.0

      vlat  = alat_ref  + dphi
      vlong = along_ref + dlambda

* Consider using the simpler sphere approximation:
*      vlat = alat_ref + vn/(6371.0 * d2r)
*      vlong = along_ref + ve/(6371.0 * d2r * 
*     :        cos(0.5 * (alat_ref + vlat) * d2r))

      vlat_out = sngl(vlat)
      vlong_out = sngl(vlong)
      
      return
      end
*-------------------- END KM2DEG ----------------------------

*-------------------- BEGIN DEG2KM_F ----------------------------
      subroutine deg2km_f( alat_sta, along_sta, alat_ref, along_ref, 
     :                       d_sta_n, d_sta_e   )
        
* convert lat, long into km north and east from a reference point

* assumes latitude between 0 and 70 degrees
* assumes west longitude is negative
* assumes angles in degrees

* Dates:  12/06/98 - written by D. Boore, based on km2deg_f
*         12/18/98 - modified to allow for negative latitudes
*         09/16/00 - Removed double precision
      
      real alat_tbl(71), b_tbl(71), adcoslat_tbl(71)
      real a, b, dphi, dlambda
      Data alat_tbl /
     : 0.000000, 1.000000, 2.000000, 3.000000, 4.000000, 5.000000,
     : 6.000000, 7.000000, 8.000000, 9.000000,10.000000,11.000000,
     :12.000000,13.000000,14.000000,15.000000,16.000000,17.000000,
     :18.000000,19.000000,20.000000,21.000000,22.000000,23.000000,
     :24.000000,25.000000,26.000000,27.000000,28.000000,29.000000,
     :30.000000,31.000000,32.000000,33.000000,34.000000,35.000000,
     :36.000000,37.000000,38.000000,39.000000,40.000000,41.000000,
     :42.000000,43.000000,44.000000,45.000000,46.000000,47.000000,
     :48.000000,49.000000,50.000000,51.000000,52.000000,53.000000,
     :54.000000,55.000000,56.000000,57.000000,58.000000,59.000000,
     :60.000000,61.000000,62.000000,63.000000,64.000000,65.000000,
     :66.000000,67.000000,68.000000,69.000000,70.000000
     :/
      Data b_tbl /
     : 1.842808, 1.842813, 1.842830, 1.842858, 1.842898, 1.842950,
     : 1.843011, 1.843085, 1.843170, 1.843265, 1.843372, 1.843488,
     : 1.843617, 1.843755, 1.843903, 1.844062, 1.844230, 1.844408,
     : 1.844595, 1.844792, 1.844998, 1.845213, 1.845437, 1.845668,
     : 1.845907, 1.846153, 1.846408, 1.846670, 1.846938, 1.847213,
     : 1.847495, 1.847781, 1.848073, 1.848372, 1.848673, 1.848980,
     : 1.849290, 1.849605, 1.849992, 1.850242, 1.850565, 1.850890,
     : 1.851217, 1.851543, 1.851873, 1.852202, 1.852531, 1.852860,
     : 1.853188, 1.853515, 1.853842, 1.854165, 1.854487, 1.854805,
     : 1.855122, 1.855433, 1.855742, 1.856045, 1.856345, 1.856640,
     : 1.856928, 1.857212, 1.857490, 1.857762, 1.858025, 1.858283,
     : 1.858533, 1.858775, 1.859008, 1.859235, 1.859452
     :/
      Data adcoslat_tbl /
     : 1.855365, 1.855369, 1.855374, 1.855383, 1.855396, 1.855414,
     : 1.855434, 1.855458, 1.855487, 1.855520, 1.855555, 1.855595,
     : 1.855638, 1.855683, 1.855733, 1.855786, 1.855842, 1.855902,
     : 1.855966, 1.856031, 1.856100, 1.856173, 1.856248, 1.856325,
     : 1.856404, 1.856488, 1.856573, 1.856661, 1.856750, 1.856843,
     : 1.856937, 1.857033, 1.857132, 1.857231, 1.857331, 1.857435,
     : 1.857538, 1.857643, 1.857750, 1.857858, 1.857964, 1.858074,
     : 1.858184, 1.858294, 1.858403, 1.858512, 1.858623, 1.858734,
     : 1.858842, 1.858951, 1.859061, 1.859170, 1.859276, 1.859384,
     : 1.859488, 1.859592, 1.859695, 1.859798, 1.859896, 1.859995,
     : 1.860094, 1.860187, 1.860279, 1.860369, 1.860459, 1.860544,
     : 1.860627, 1.860709, 1.860787, 1.860861, 1.860934
     :/

      pi = 4.0*atan(1.0)
      d2r = pi/ 180.0

* interpolate to find proper arc distance:

      call locate( alat_tbl, 71, abs(alat_ref), j)

      b = b_tbl(j) + (abs(alat_ref)-alat_tbl(j))*
     :  (b_tbl(j+1)-b_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      adcoslat = adcoslat_tbl(j) + (abs(alat_ref)-alat_tbl(j))*
     :  (adcoslat_tbl(j+1)-adcoslat_tbl(j))/
     :  (alat_tbl(j+1)-alat_tbl(j))

      a = adcoslat * cos(d2r*abs(alat_ref))

* compute lat,long relative to reference:
      dphi = alat_sta - alat_ref
      dlambda = along_sta - along_ref

* convert from degrees to minutes of arc:
      dlambda = dlambda * 60.0
      dphi    = dphi    * 60.0

* compute distances (positive ve corresponds to increase in longitude:
*                    vn positive to the north, ve positive to the east)
      d_sta_e =  a * dlambda 
      d_sta_n =  b * dphi

* Consider replacing the above computation with the following simple
* computation based on assuming that the Earth is a perfect sphere:
*      vn = (alat_sta - alat_ref)*d2r*6371.0
*      ve = (along_sta - along_ref)*d2r*
*     :      cos(0.5*(alat_sta+alat_ref)*d2r)*6371.0

      return
      end
*-------------------- END DEG2KM_F ----------------------------

* --------------------- BEGIN YINTRF ------------------------------------
      function yintrf( x, xin, yin, n)
c
c returns an interpolated value (yintrf) based on straight line
c interpolation of the data in xin and yin.

* Needs Numerical recipe routine locate

c
c dates:  3/14/85 - written
*        11/30/95 - substituted LOCATE instead of starting from beginning
*                   each time
*        03/13/96 - added code to deal with xin increasing or decreasing
*        12/12/00 - Stripped off "locate.for"

      dimension xin(*), yin(*)
      logical incrs

* Is xin increasing or decreasing?
      incrs = .true.
      if (xin(n) .lt. xin(1)) incrs = .false.

* Set value if x is outside the range of xin:
      if (incrs) then
        if ( x .le. xin(1) ) then
            yintrf = yin(1)
            return
        end if
        if ( x .ge. xin(n) ) then
            yintrf = yin(n)
            return
        end if  
      else
        if ( x .ge. xin(1) ) then
            yintrf = yin(1)
            return
        end if
        if ( x .le. xin(n) ) then
            yintrf = yin(n)
            return
        end if  
      end if

* Locate the proper cell and interpolate:
      call locate(xin, n, x, j)
      yintrf = yin(j) + (x-xin(j))*(yin(j+1) - yin(j))/
     * (xin(j+1)-xin(j))

      return
      end
* --------------------- END YINTRF ------------------------------------

* ----------------------------- BEGIN FORK --------------------------
      SUBROUTINE FORK(LX,CX,SIGNI)
C FAST FOURIER                                  2/15/69
C                          LX
C    CX(K) = SQRT(1.0/LX)* SUM (CX(J)*EXP(2*PI*SIGNI*I*(J-1)*(K-1)/LX))
C                          J=1                        FOR K=1,2,...,LX
C
C  THE SCALING BETWEEN FFT AND EQUIVALENT CONTINUUM OUTPUTS
C  IS AS FOLLOWS.
C
C
C     GOING FROM TIME TO FREQUENCY:
C             F(W)=DT*SQRT(LX)*CX(K)
C
C                  WHERE W(K)=2.0*PI*(K-1)*DF

*                  and    DF = 1/(LX*DT)
C
C
C     GOING FROM FREQUENCY TO TIME, WHERE THE FREQUENCY
C     SPECTRUM IS GIVEN BY THE DIGITIZED CONTINUUM SPECTRUM:
C
C             F(T)=DF*SQRT(LX)*CX(K)
*
C                  WHERE T(K)=(K-1)*DT
C
C
C  THE RESULT OF THE SEQUENCE...TIME TO FREQUENCY,POSSIBLE MODIFICATIONS
C  OF THE SPECTRUM (FOR FILTERING,ETC.), BACK TO TIME...
C  REQUIRES NO SCALING.
C
C
C  THIS VERSION HAS A SLIGHT MODIFICATION TO SAVE SOME TIME...
C  IT TAKES THE FACTOR 3.1415926*SIGNI/L OUTSIDE A DO LOOP (D.BOORE 12/8
C  FOLLOWING A SUGGESTION BY HENRY SWANGER).
C

* Some brief notes on usage:

* "signi" is a real variable and should be called either with the value "+1.0"
* of "-1.0".  The particular value used depends on the conventions being used
* in the application (e.g., see Aki and Richards, 1980, Box 5.2, pp. 129--130).

* Time to frequency:
* In calling routine,
* 
*       do i = 1, lx
*         cx(i) = CMPLX(y(i), 0.0)
*       end do
*  where y(i) is the time series and lx is a power of 2
* 
*  After calling Fork with the complex array specified above, the following 
* symmetries exist:
* 
*        cx(1)        = dc value (f = 0 * df, where df = 1.0/(lx*dt))
*        cx(lx/2 + 1) = value at Nyquist (f = (lx/2+1-1)*df = 1.0/(2*dt))
*        cx(lx)       = CONJG(cx(2))
*        cx(lx-1)     = CONJG(cx(3))
*         |           =      |
*        cx(lx-i+2)   = CONJG(cx(i))
*         |           =      |
*        cx(lx/2+2)   = CONJG(cx(lx/2))
* 
* where "CONJG" is the Fortran complex conjugate intrinsic function
* 
* This symmetry MUST be preserved if modifications are made in the frequency 
* domain and another call to Fork (with a different sign for signi) is used
* to go back to the time domain.  If the symmetry is not preserved, then the
* time domain array will have nonzero imaginary components.  There is one case
* where advantage can be taken of this, and that is to find the Hilbert 
* transform and the window of a time series with only two calls to Fork (there 
* is a short note in BSSA {GET REFERENCE} discussing this trick, which amounts 
* to zeroing out the last half of the array and multiplying all but the dc and 
* Nyquist values by 2.0; in the time domain, REAL(cx(i)) and AIMAG(cx(i)) 
* contain the filtered (if a filter was applied) and Hilbert transform of the 
* filtered time series, respectively, while CABS(cx(i)) and ATAN2(AIMAG(cx(i)), 
* REAL(cx(i))) are the window and instantaneous phase of the filtered time 
* series, respectively.

* Some references:

* Farnbach, J.S. (1975). The complex envelope in seismic signal analysis, 
* BSSA 65, 951--962. 
* He states that the factor of 2 is applied for i = 2...npw2/2 (his indices 
* start at 0, I've added 1), which is different than the next reference:

* Mitra, S.K. (2001). Digital Signal Processing, McGraw-Hill, New York.
* He gives an algorithm on p. 794 (eq. 11.81), in which the factor of 2 is 
* applied from 0 frequency to just less than Nyquist.

* 
* The easiest way to ensure the proper symmetry is to zero out the
* last half of the array (as discussed above), but the following is what
* I usually use:  
* modify (filter) only half
* of the cx array:
* 
*       do i = 1, lx/2
*         cx(i) = filter(i)*cx(i)
*       end do
* 
* where "filter(i)" is a possibly complex filter function (and recall that 
* the frequency corresponding to i is f = float(i-1)*df).  After this, fill out
* the last half of the array using
*       
*       do i = lx/2+2, lx
*         cx(i) = CONJG(cx(lx+2-j))
*       end do
* 
* Note that nothing is done with the Nyquist value.  I assume (but am not sure!)
* that this value should be 0.0
* 
* Dates: xx/xx/xx - Written by Norm Brenner(?), Jon Claerbout(?)
*        12/21/00 - Replaced hardwired value of pi with pi evaluated here,
*                     and added comments regarding usage.  Also deleted
*                     dimension specification of cx(lx) and replace it with
*                     cx(*) in the type specification statement.  I also
*                     cleaned up the formatting of the subroutine.
*        08/28/01 - Added comment about variable "signi" being real, and 
*                   added "float" in equations for "sc" and "temp", although 
*                   not strictly required.
*        06/19/02 - Added some comments about computing envelopes and
*                   instantaneous frequencies
             
      complex cx(*),carg,cexp,cw,ctemp

      pi = 4.0*atan(1.0)

      j=1
      sc=sqrt(1./float(lx))

      do i=1,lx
        if(i.gt.j) go to 2
        ctemp=cx(j)*sc
        cx(j)=cx(i)*sc
        cx(i)=ctemp
2       m=lx/2
3       if(j.le.m) go to 5
        j=j-m
        m=m/2
        if(m.ge.1) go to 3
5       j=j+m
      end do

      l=1
6     istep=2*l
      temp= pi * signi/float(l)

      do m=1,l
        carg=(0.,1.)*temp*(m-1)
        cw=cexp(carg)
        do i=m,lx,istep
          ctemp=cw*cx(i+l)
          cx(i+l)=cx(i)-ctemp
          cx(i)=cx(i)+ctemp
        end do
      end do

      l=istep
      if(l.lt.lx) go to 6

      return
      end
* ----------------------------- END FORK --------------------------
! ------------------------------------------------------------- Get_NPW2
      subroutine get_npw2(npts,signnpw2,npw2)

* Find npw2 (less than npts if signnpw2 < 0)

* Dates: 12/12/00 - Written by D. Boore
!        04/04/10 - Correct error that it does not return
!                   npw2 = npts, if npts is a power of 2.

      npw2_exp = int( alog(float(npts))/alog(2.0) )
      if (signnpw2 < 0.0) then
        npw2 = 2.0**npw2_exp
      else 
        npw2_temp = 2.0**npw2_exp
        if (npw2_temp == npts) then 
          npw2 = npts
        else
          npw2 = 2.0**(npw2_exp+1)
        end if
      end if

      return
      end
! ------------------------------------------------------------- Get_NPW2
*----------------- BEGIN RDCALCDP -----------------------------
      subroutine rdcalcdp(acc,na,omega,damp,dt,rd)
* This is a modified version of "Quake.For", originally
* written by J.M. Roesset in 1971 and modified by
* Stavros A. Anagnostopoulos, Oct. 1986.  The formulation is that of
* Nigam and Jennings (BSSA, v. 59, 909-922, 1969).  This modification 
* eliminates the computation of the relative velocity and absolute 
* acceleration; it returns only the relative displacement.  

*   acc = acceleration time series
*    na = length of time series
* omega = 2*pi/per
*  damp = fractional damping (e.g., 0.05)
*    dt = time spacing of input
*    rd = relative displacement of oscillator

* Dates: 05/06/95 - Modified by David M. Boore
*        04/15/96 - Changed name to RD_CALC and added comment lines
*                   indicating changes needed for storing the oscillator 
*                   time series and computing the relative velocity and 
*                   absolute acceleration
!        03/11/01 - Double precision version of Rd_Calc
*        01/31/03 - Moved implicit statement before the type declarations

      implicit real*8 (a - h, o - z)

      real*4 acc(*)

      real*4 omega, damp, dt, rd 

      omt=omega*dt
      d2=1-damp*damp
      d2=dsqrt(d2)
      bom=damp*omega
*      d3 = 2.*bom                 ! for aa
      omd=omega*d2
      om2=omega*omega
      omdt=omd*dt
      c1=1./om2
      c2=2.*damp/(om2*omt)
      c3=c1+c2
      c4=1./(omega*omt)
      ss=dsin(omdt)
      cc=dcos(omdt)
      bomt=damp*omt
      ee=dexp(-bomt)
      ss=ss*ee
      cc=cc*ee
      s1=ss/omd
      s2=s1*bom
      s3=s2+cc
      a11=s3
      a12=s1
      a21=-om2*s1
      a22=cc-s2
      s4=c4*(1.-s3)
      s5=s1*c4+c2
      b11=s3*c3-s5
      b12=-c2*s3+s5-c1
      b21=-s1+s4
      b22=-s4
      rd=0.
*      rv = 0.                                                                   ! for rv
*      aa = 0.                                                                   ! for aa
      n1=na-1
      y=0.
      ydot=0.
      do i=1,n1
        y1=a11*y+a12*ydot+b11*acc(i)+b12*acc(i+1)
        ydot=a21*y+a22*ydot+b21*acc(i)+b22*acc(i+1)
        y=y1                                                                    ! y is the oscillator output at time corresponding to index i
        z=dabs(y)
        if (z.gt.rd) rd=z
*        z1 = dabs(ydot)                                                         ! for rv
*        if (z1.gt.rv) rv = z1                                                   ! for rv
*        ra = -d3*ydot -om2*y1                                                   ! for aa
*        z2 = dabs(ra)                                                           ! for aa
*        if (z2.gt.aa) aa = z2                                                   ! for aa
      end do
      return
      end
*----------------- END RDCALCDP -----------------------------


* --------------- BEGIN TIME_DIFF ---------------------------------
      subroutine time_diff(time_start, time_stop, time_elapsed)

* Dates: 02/18/09 - Written by D.M. Boore
* To be used with
* Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
*      character datx*8, timx*10
*      call DATE_AND_TIME( datx, timx )
* Date is returned as 'CCYYMMDD'
* Time is returned as 'hhmmss.sss'

      implicit none
      character, intent(in) :: time_start*(*), time_stop*(*)
      real, intent(out) :: time_elapsed
      real ::  secb, sece 
      integer :: ihb, imb, ihe, ime

      read(time_start(1:10),'(i2,i2,f6.3)') 
     :                       ihb, imb, secb
      read(time_stop(1:10),'(i2,i2,f6.3)') 
     :                       ihe, ime, sece
      time_elapsed = 
     :  3600.0*float(ihe-ihb) + 60.0*float(ime-imb) + sece-secb 

      end subroutine time_diff
* --------------- END TIME_DIFF ---------------------------------


*  ------------------- BEGIN BUTTRLCF -------------------------------
      function buttrlcf( f, fcut, norder)
c
c Computes the response of an norder, bidirectional
* high-pass Butterworth filter.  This is the filter
* used by the AGRAM processing (the equation was
* provided by April Converse).

* Modification: 3/27/89 - created by modifying HiPassF

      buttrlcf = 1.0
      if ( fcut.eq.0.0 ) return

      buttrlcf = 0.0

      if ( f .eq. 0.0) return

      buttrlcf = 1.0/ (1.0+(fcut/f)**(2.0*norder))

      return
      end
*  ------------------- END BUTTRLCF -------------------------------
*----------------- BEGIN Acc2V -----------------------------
      subroutine acc2v(acc, npts, dt, rmv_trnd, vel)


* Compute velocity time series from acceleration,
* assuming that the acceleration
* is represented by straight lines connecting the digitized values.

* Dates: 01/06/09 - Written by D.M. Boore, patterned after smc2vd

      real acc(*), vel(*) 
      logical rmv_trnd
      double precision cumv, a1, a2,
     : ddt, ddt_2

      if (rmv_trnd) then      
* remove trend first (straight line between first and last points)
* Note: acc is replaced with detrended time series
*        call dcdt(acc, dt, npts, 1, npts, .false., .true.)  ! old routine,
*                                                         ! gives steps at ends

         call rmvtrend(acc, npts)
      end if

* compute velocity and displacement, using analytical formulas based
* on representing the acceleration as a series of straightline segments.

      ddt     = dble(dt)
      ddt_2   = 0.5d0 * dble(dt)

      cumv = 0.0d0

      vel(1) = sngl(cumv)
      do j=2,npts
        a1 = dble(acc(j-1))
        a2 = dble(acc(j))
        cumv = cumv + (a1 + a2)*ddt_2
        vel(j) = sngl(cumv)
      end do

      return
      end
*----------------- END Acc2V -----------------------------
* ----------------------------- BEGIN RMVTREND ----------------
      subroutine rmvtrend(y, n)

* Removes a straightline fit to first and last points, replacing
* the input array with the detrended array

* Dates: 02/09/99 - written by D. Boore


      real y(*)

      y1 = y(1)
      y2 = y(n)
      slope = (y2 - y1)/float(n-1)

      do i = 1, n
        y(i) = y(i) - (y1 + slope*float(i-1))
      end do

      return
      end
* ----------------------------- END RMVTREND ----------------
*----------------- BEGIN GSPRD_F -----------------------------
* Calculate geometric spreading
* DMB added entry points so that the program could be called using a single 
* argument, without passing the other arguments through common.  Using
* a single argument is necessary when called by sme Numerical Recipes programs.

* Use:

* Call the setup entry point:
*      dummy =  gsprd_f_setup(r_ref,nsprd_segs,rlow,
*     :                     a_s,b_s,m_s,
*     :                     numsource,amag)

* Call as a function:
*            gsprd_n = gsprd_f(rn)

* Deallocate arrays:
*      dummy =  gsprd_f_deallocate()





* Dates: 06/07/95 - Written by D.M. Boore
*        07/02/99 - Added magnitude-dependent "depth" from Atkinson
*                   and Silva, which required adding some parameters to
*                   the passed arguments
*        06/05/00 - Added some explanation of r
*        06/08/00 - Make gsprd nondimensional through the use of r_ref, which 
*                   now appears in the definition of variable const
*                   in const_am0_gsprd
*        01/27/02 - Following Silva, parameters added to allow magnitude
*                   dependent slope (to capture finite fault effects)
*        11/27/05 - Remove deff for Atkinson (2005) source
*        04/24/07 - Put "rmod = r" in the if statement
*        11/13/08 - Redo the function so that it can be used in Numerical
*                   Recipes routines, which assume calls to function(x).
*                   Added entry points rather than using common blocks
*                   to do this (using Larry Baker's help).

      function gsprd_f(r)
     
      save
      
      real rlow_init(*), a_s_init(*), 
     :                            b_s_init(*), 
     :                            m_s_init(*)
      real, allocatable :: rlow(:), a_s(:), b_s(:), m_s(:)
      real geff(10)
      
* Note that generally r = hypocentral distance.  For Atkinson and Silva 
* (BSSA 90, 255--274) r is the closest distance to the fault plane ("d" in 
* their paper; below their eq. 4), so that rmod is, in effect, accounting
* source depth twice.  See comments in AS00 section of subroutine
* spect_scale

      
      if (numsource .eq. 9 ) then                                               ! Atkinson and Silva (2000)
        deff = 10.0**(-0.05 + 0.15 * amag)
        rmod = sqrt(r**2 + deff**2)        
      else      
        rmod = r      
      end if
      
      geff(1) = r_ref/rlow(1)                                                   ! usually set r_ref = 1.0 km.  Be careful
                                                                                ! if a different value or different units are
                                                                                ! used.  In particular, using different units
                                                                                ! will require a change in the scaling factor
                                                                                ! of 1.e-20 used in the definition of const in
                                                                                ! const_am0_gsprd

      do i = 2, nsprd_segs
        slope = a_s(i-1) + b_s(i-1)*(amag - m_s(i-1))
        geff(i) = geff(i-1)*(rlow(i)/rlow(i-1))**slope
      end do
      if (rmod .le. rlow(1)) then
        j = 1
      else if (rmod .ge. rlow(nsprd_segs)) then
        j = nsprd_segs
      else
        call locate(rlow, nsprd_segs, rmod, j)
      end if
      slope = a_s(j) + b_s(j)*(amag - m_s(j))

      gsprd_f = (geff(j)) * (rmod/rlow(j))**slope
      
      return

      entry gsprd_f_setup(r_ref_init,nsprd_segs_init,rlow_init,
     :                  a_s_init,b_s_init,m_s_init,
     :                  numsource_init,amag_init)
     
      allocate(rlow(nsprd_segs_init), 
     :                                a_s(nsprd_segs_init), 
     :                                b_s(nsprd_segs_init),  
     :                                m_s(nsprd_segs_init))
      r_ref                    = r_ref_init
      nsprd_segs               = nsprd_segs_init
      rlow       = rlow_init(1:nsprd_segs_init) 
      a_s        = a_s_init(1:nsprd_segs_init) 
      b_s        = b_s_init(1:nsprd_segs_init) 
      m_s        = m_s_init(1:nsprd_segs_init)
      numsource  = numsource_init
      amag       = amag_init
      
      return
      
      entry gsprd_f_deallocate

      deallocate(rlow, a_s, b_s, m_s)
      
      gsprd_f_deallocate = 1.0
      
      return
     
      
      end
*----------------- END GSPRD_F -----------------------------

*----------------- begin q_f -----------------------------
* Calculates anelastic attenuation
* DMB added entry points so that the program could be called using a single 
* argument, without passing the other arguments through common.  using
* a single argument is necessary when called by sme numerical recipes programs.

* use:

* call the setup entry point:
*      dummy = q_f_setup(fr1, qr1, s1, ft1, ft2,
*     :              fr2, qr2, s2, c_q)

* call as a function:
*      q_fq = q_f(fq)



* dates: 06/07/95 - written by d.m. boore
*        12/14/95 - added check for equal transition frequencies
*        05/11/07 - removed "\smsim\" from include statements
*        11/14/08 - redo the function so that it can be used in numerical
*                   recipes routines, which assume calls to function(x).
*                   added entry points rather than using common blocks
*                   to do this (using larry baker's help).

!modified/added by ka, mar.2012. return to q0*f**eta definition
      function q_f(f)                                                           !modified/added by ka, mar.2012. return to q0*f**eta definition
                                                                                !modified/added by ka, mar.2012. return to q0*f**eta definition
      save                                                                      !modified/added by ka, mar.2012. return to q0*f**eta definition
                                                                                !modified/added by ka, mar.2012. return to q0*f**eta definition
      q_f = 9999.0                                                              !modified/added by ka, mar.2012. return to q0*f**eta definition
      if (f .eq. 0.0) return                                                    !modified/added by ka, mar.2012. return to q0*f**eta definition
                                                                                !modified/added by ka, mar.2012. return to q0*f**eta definition
      q_f=q0*f**qeta                                                            !modified/added by ka, mar.2012. return to q0*f**eta definition
      if(q_f.lt.qmin)q_f=qmin                                                   !modified/added by ka, mar.2012. return to q0*f**eta definition
                                                                                !modified/added by ka, mar.2012. return to q0*f**eta definition
      return                                                                    !modified/added by ka, mar.2012. return to q0*f**eta definition
                                                                                !modified/added by ka, mar.2012. return to q0*f**eta definition
      entry q_f_setup(qmin_init,q0_init,qeta_init)                              !modified/added by ka, mar.2012. return to q0*f**eta definition
                                                                                !modified/added by ka, mar.2012. return to q0*f**eta definition
      qmin=     qmin_init                                                       !modified/added by ka, mar.2012. return to q0*f**eta definition
      q0  =     q0_init                                                         !modified/added by ka, mar.2012. return to q0*f**eta definition
      qeta=     qeta_init                                                       !modified/added by ka, mar.2012. return to q0*f**eta definition
                                                                                !modified/added by ka, mar.2012. return to q0*f**eta definition
      return                                                                    !modified/added by ka, mar.2012. return to q0*f**eta definition
                                                                                !modified/added by ka, mar.2012. return to q0*f**eta definition
      end                                                                       !modified/added by ka, mar.2012. return to q0*f**eta definition
*----------------- end q_f -----------------------------

      subroutine SCEC_freqs(nfreq,freq)
      real, dimension(500) :: freq
      integer              :: nfreq
      nfreq=63

      freq(1)=100.00    
      freq(2)=90.909 
      freq(3)=83.333 
      freq(4)=76.923 
      freq(5)=66.667 
      freq(6)=58.824 
      freq(7)=50.000
      freq(8)=45.455 
      freq(9)=40.000     
      freq(10)=34.483
      freq(11)=31.250 
      freq(12)=28.571
      freq(13)=25.000    
      freq(14)=22.222
      freq(15)=20.000    
      freq(16)=18.182
      freq(17)=16.667
      freq(18)=15.385
      freq(19)=13.333
      freq(20)=11.765
      freq(21)=10.000    
      freq(22)=9.091 
      freq(23)=8.333 
      freq(24)=7.692 
      freq(25)=6.667 
      freq(26)=5.882 
      freq(27)=5.000     
      freq(28)=4.545 
      freq(29)=4.167 
      freq(30)=3.846 
      freq(31)=3.571 
      freq(32)=3.333 
      freq(33)=2.857 
      freq(34)=2.500   
      freq(35)=2.222 
      freq(36)=2.000     
      freq(37)=1.818 
      freq(38)=1.667 
      freq(39)=1.538 
      freq(40)=1.333 
      freq(41)=1.176 
      freq(42)=1.000     
      freq(43)=0.909 
      freq(44)=0.833 
      freq(45)=0.769 
      freq(46)=0.6667 
      freq(47)=0.5883 
      freq(48)=0.500   
      freq(49)=0.4546 
      freq(50)=0.4167 
      freq(51)=0.38467 
      freq(52)=0.3572 
      freq(53)=0.3333 
      freq(54)=0.2857 
      freq(55)=0.250  
      freq(56)=0.22725 
      freq(57)=0.200   
      freq(58)=0.18182 
      freq(59)=0.16667 
      freq(60)=0.15385 
      freq(61)=0.13334 
      freq(62)=0.11765
      freq(63)=0.100   
               
      return
      end

      
      
      subroutine write_slip(iunit,actualSlip,nl,nw)
      real actualSlip(200,200)
      write(iunit,*)
      write(iunit,*)'Slip distribution for this run:'
      do j=1,nw
          write(iunit,'(200(f6.3,2x))')(actualSlip(i,j),i=1,nl)
      enddo
      return
      end



      real function AB14_fac(freq,dist,dep)                                     !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      real                                  ::freq, dist, dep                   !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      real                                  ::Clf, Tc, h                        !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
                                                                                !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      pi=4.*atan(1.)                                                            !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      h=dep                                                                     !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
                                                                                !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      if(dist.le.dep)then                                                       !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
          if(h.eq.1.)h=0.99                                                     !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
          Clf=0.2*cos((pi*(dist-dep))/(2.*(1-h)))                               !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      else                                                                      !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
          if(h.eq.50.)h=49.99                                                   !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
          Clf=0.2*cos((pi*(min(dist,50.)-dep))/(2.*(50.-h)))                    !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      endif                                                                     !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
                                                                                !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      Tc=max((1.-1.429*log10(max(freq,1.))),0.)                                 !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
                                                                                !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      AB14_fac=10**(Clf*Tc)                                                     !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
                                                                                !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      return                                                                    !Added by KA, Jan.2014 to handle frequency dependent geom_sprd
      end                                                                       !Added by KA, Jan.2014 to handle frequency dependent geom_sprd

