! ------------------------------------------------------------------- SMC2FS2
      Program SMC2FS2


! Computes the Fourier amplitude spectra for specified time series.
! Writes output in smc or column format.

! The output file name has the same name as the input file, with the specified
! character string appended to the name.

! Input from a control file; here is an example:
!!Control file for program SMC2FS2
!! Revision of program involving a change in the control file on this date:
!   03/10/10
!! As many comment lines as desired, each starting with "!"
!! The string "pp:" indicates a new set of processing parameters
!! to be applied to the following smc files.  The parameters are given on the
!! lines following "pp:", until the next "pp:" line or until "stop" is
!! encountered.
!! NOTE: Use the tapers with caution, choosing them so that important signal
!! is not reduced by the tapering.  This can be particularly a problem with
!! analog data from relatively small earthquakes that triggered near the
!! S-wave arrival.
!!
!! -----------------------------------------------------------------------------
!!
!! Meaning of smoothing input parameters
!!
!! NO SMOOTHING
!! itype = 0
!! SMOOTHING OVER EQUALLY SPACED FREQUENCIES
!! itype = 1: box weighting function
!!   smooth_param = width of box weighting function (Hz)
!! itype = 2: triangular weighting function
!!   smooth_param = width of triangular weighting function (Hz)
!! SMOOTHING OVER LOGARITHMICALLY SPACED FREQUENCIES
!! itype = 3: box weighting function
!!   smooth_param = xi, which is the fraction of a decade for the
!!                  box weighting function
!! itype = 4: triangular weighting function
!!   smooth_param = xi, which is the fraction of a decade for the
!!                  triangular weighting function
!! itype = 5: Konno and Ohmachi weighting function (see BSSA 88, 228-241)
!!   smooth_param = xi, which is the fraction of a decade for which
!!                  the Konno and Ohmachi weighting function is greater
!!                  than 0.043.(it is related to
!!                  their smoothing parameter b by the equation
!!                  b = 4.0/smooth_param, so we have this correspondence between
!!                  b and smooth_param
!!                         b smooth_param
!!                        10         0.40
!!                        20         0.20
!!                        40         0.10
!!
!!                  b = 40 seems to be commonly used, but I do not think that it
!!                  gives enough smoothing; I PREFER SMOOTH_PARAM = 0.2,
!!                  corresponding to b = 20.
!!
!! ipow = power of FAS to be smoothed (2 = smoothing energy spectrum)
!!
!! dx_smooth: Note: need df_smooth for linearly-spaced smoothers,
!! and generally it should be the df from the fft.  For general x data, it is
!! the spacing between x values, assumed to be constant,  The reason for
!! including it as an input parameter is to "fool" the
!! program to do smoothing over a specified number of points by
!! setting df_smooth = 1 and smooth_param = number of points (including
!! points with zero weight at ends; e.g., smooth_param = 5 will
!! give a smoother with weights 0, 1/4, 2/4, 1/4, 0; smooth_param
!! should be odd).
!!
!! -----------------------------------------------------------------------------
!! Meaning of frequency specification parameters:
!!
!!SPECIFY_FREQUENCIES? (y/n):
!! <enter Y or N>
!!FREQUENCY SPECIFICATION:
!!  If specify_frequencies = Y, then enter the
!!    number of frequencies, freq(1), freq(2)..., freq(nfreq)
!!  If specify_frequencies = N, then enter
!!    f_low, f_high, log-spaced (0=N, 1=Y), freq_param
!!         if freq_param = 0.0, there is no interpolation, and the FFT frequencies
!!            are used between f_low and f_high (log-spaced is ignored).
!!         if freq_param /= 0.0 and log-spaced = 0, then freq_param is the spacing of the
!!            interpolated frequencies between f_low and f_high
!!         if freq_param /= 0.0 and log-spaced = 1, then freq_param is the number of
!!            interpolated frequencies between f_low and f_high (NOTE: f_low must be > 0.0)
!! -----------------------------------------------------------------------------
!!
!!Name of summary file:
! smc2fs2.sum
!PP: new set of parameters
!!tskip, tlength
!   0.0 2000.0
!!dc_remove?
!  .true.
!!Length of taper at beginning and end of time series, before adding zeros
!! to make the number of points in the record a power of two.
! 0.0 0.0
!!signnpw2(<0, backup for npw2, no zpad):
! +1.0
!!smoothing: itype, ipow, df_smooth (0 = FFT df), smooth_param
!! (see above for the meaning of these input parameters):
!   0 1 0.0 0.20
!!SPECIFY_FREQUENCIES? (y/n):
!  N
!!FREQUENCY SPECIFICATION
!  0.0 50.0 0 0.0
!!character string to append to filename:
!   .no_smooth.use_fft_frequencies_from_0.0_to_50_hz.fs.col
!!Output in smc format (Y,N)?
!! ***IMPORTANT NOTE: Output cannot be in smc format if use log-spaced
!! frequencies because programs such as smc2asc have not been modified
!! to deal with log-spaced frequency or if specify_frequencies = Y
! n
!!Files to process:
! A11.V01.a5
!PP: new set of parameters
!!tskip, tlength
!   0.0 2000.0
!!dc_remove?
!  .true.
!!Length of taper at beginning and end of time series, before adding zeros
!! to make the number of points in the record a power of two.
! 0.0 0.0
!!signnpw2(<0, backup for npw2, no zpad):
! +1.0
!!smoothing: itype, ipow, df_smooth (0 = FFT df), smooth_param
!! (see above for the meaning of these input parameters):
!   5 1 0.0 0.20
!!SPECIFY_FREQUENCIES? (y/n):
!  Y
!!FREQUENCY SPECIFICATION
!   6
!   0.2 0.5 1.0 2.0 5.0 10.0
!!character string to append to filename:
!   .k-o_smooth_0.20.specify_6_frequencies.fs.col
!!Output in smc format (Y,N)?
!! ***IMPORTANT NOTE: Output cannot be in smc format if use log-spaced
!! frequencies because programs such as smc2asc have not been modified
!! to deal with log-spaced frequency or if specify_frequencies = Y
! n
!!Files to process:
! A11.V01.a5
!PP: new set of parameters
!!tskip, tlength
!   0.0 2000.0
!!dc_remove?
!  .true.
!!Length of taper at beginning and end of time series, before adding zeros
!! to make the number of points in the record a power of two.
! 0.0 0.0
!!signnpw2(<0, backup for npw2, no zpad):
! +1.0
!!smoothing: itype, ipow, df_smooth (0 = FFT df), smooth_param
!! (see above for the meaning of these input parameters):
!   5 1 0.0 0.20
!!SPECIFY_FREQUENCIES? (y/n):
!  N
!!FREQUENCY SPECIFICATION
!  0.0 50.0 0 0.0
!!character string to append to filename:
!   .k-o_smooth_0.20.use_fft_frequencies_from_0.0_to_50_hz.fs.col
!!Output in smc format (Y,N)?
!! ***IMPORTANT NOTE: Output cannot be in smc format if use log-spaced
!! frequencies because programs such as smc2asc have not been modified
!! to deal with log-spaced frequency or if specify_frequencies = Y
! n
!!Files to process:
! A11.V01.a5
!PP: new set of parameters
!!tskip, tlength
!   0.0 2000.0
!!dc_remove?
!  .true.
!!Length of taper at beginning and end of time series, before adding zeros
!! to make the number of points in the record a power of two.
! 0.0 0.0
!!signnpw2(<0, backup for npw2, no zpad):
! +1.0
!!smoothing: itype, ipow, df_smooth (0 = FFT df), smooth_param
!! (see above for the meaning of these input parameters):
!   5 1 0.0 0.20
!!SPECIFY_FREQUENCIES? (y/n):
!  N
!!FREQUENCY SPECIFICATION
!   0.2 10.0 0 0.2
!!character string to append to filename:
!   .k-o_smooth_0.20.lin_spaced_freqs_from_0.2_to_10.0_hz_df_0.2.fs.col
!!Output in smc format (Y,N)?
!! ***IMPORTANT NOTE: Output cannot be in smc format if use log-spaced
!! frequencies because programs such as smc2asc have not been modified
!! to deal with log-spaced frequency or if specify_frequencies = Y
! n
!!Files to process:
! A11.V01.a5
!PP: new set of parameters
!!tskip, tlength
!   0.0 2000.0
!!dc_remove?
!  .true.
!!Length of taper at beginning and end of time series, before adding zeros
!! to make the number of points in the record a power of two.
! 0.0 0.0
!!signnpw2(<0, backup for npw2, no zpad):
! +1.0
!!smoothing: itype, ipow, df_smooth (0 = FFT df), smooth_param
!! (see above for the meaning of these input parameters):
!   5 1 0.0 0.20
!!SPECIFY_FREQUENCIES? (y/n):
!  N
!!FREQUENCY SPECIFICATION
!  0.2 10.0 1 6
!!character string to append to filename:
!   .k-o_smooth_0.20.six_log_spaced_freqs_from_0.2_to_10.0_hz.fs.col
!!Output in smc format (Y,N)?
!! ***IMPORTANT NOTE: Output cannot be in smc format if use log-spaced
!! frequencies because programs such as smc2asc have not been modified
!! to deal with log-spaced frequency or if specify_frequencies = Y
! n
!!Files to process:
! A11.V01.a5
!   stop

! Dates: 06/15/01 - Written by D.M. Boore, patterned after SMC2FAS
!        06/27/01 - Replaced "_r" in f_fs2 with "_f".
!        07/10/01 - Write more stuff regarding processing parameters to
!                   comments in file
!        07/23/01 - Changed filenaming convention.  Now I use the smc file,
!                   with the first two characters of the extension replaced
!                   by "fs".  I did this because otherwise duplicate filenames
!                   were being created.  For example, file.as1 and file.as2
!                   would have had the same output file names before
!                   (file.fs2).
!        10/19/01 - Obtain default name of control file via CR
!        11/11/01 - Write flow, df, fhigh, nfreq to comments
!        12/18/01 - Use double precision version of FORK (FORKDP) (had to
!                   modify absspect.for)
!        12/18/01 - Removed double precision version; it did not seem to
!                   make a difference
!        12/18/01 - Store 0.0 if value is less than what can be written in the
!                   e format.
!        05/31/02 - More changes in file naming conventions.
!        07/17/03 - Add sps to smc_npts parameter list
!        12/08/03 - changed way of naming output (floowing blpadflt): tack on
!                   a specifed character string to the input file name.  This
!                   works because I now compile using LF95.
!        12/13/03 - Allow column file format as an option
!        12/29/03 - Add nfreq to column file output
!        10/29/04 - Modification so will reach desired fhigh
!        07/02/07 - Major change to structure of the control file, and include
!                   new, more generalized smoothing subroutine
!        07/03/07 - I deleted get path stuff
!        10/16/07 - Transparent change to user needed because smcread can
!                   read two input formats
!        03/20/08 - Switch order of smoothing parameters on input.
!        03/20/08 - Add possibility of log-spaced frequencies.
!        03/24/08 - Write elapsed time for processing each record to a summary
!                   file.
!        03/24/08 - Add checks of f_intrp_high > fnyq
!        03/27/08 - Combine interpolation and smoothing to eliminate smoothing
!                   at unnecessary points.  This should speed up the Konno and
!                   Ohmachi smoother.
!        04/21/08 - Combined processing steps into a subroutine.
!        07/11/09 - Redefined the smoothing parameter for Konno & Ohmachi
!                   in subroutine smooth_interpolate, and added version number
!                   to input file.
!        02/09/10 - Replace prctnfrnttaper with length of tapers at front and back.
!                 - This required changes in the control file and FAS subroutines.
!                 - Use this opportunity to "modernize" the code slightly.
!        02/15/10 - Allow specification of individual frequencies.
!                 - Allocated freq_out more intelligently.
!        03/10/10 - Increase length of f_smc and f_fs2 to 150 characters
!        03/13/10 - Use a less confusing (I hope) specification of the frequency
!                   interpolation information.
!        03/16/10 - More debugging.
!        03/17/10 - Trap for smc output but specify_frequencies = Y or log-spaced
!                   frequencies; reset output to a column file in those cases.
!        04/02/10 - Print length of time series for FFT (npw2) in summary file.
!        12/29/10 - Recompile, because now write output freq and fas using more precision.
!        02/07/11 - Delete duplicate date_ctl_correct assignment
!        06/04/11 - Allow for longer f_smc file name
!        05/04/12 - Allow more decimal places in output
!        06/30/13 - Use Fortran 95 extension getcl to get a command line argument that determines if
!                   the default control file name is used (this will be used for calling the main
!                   program within a batch file (c:\forprogs\smc2fs2 -default_ctl)
!        07/08/14 - Use system routine CPU_Time to obtain elapsed time

!
! Dimension and declaration statements:

      real fas(:), fas_out(:), freq(:), freq_out(:), ts(:)
      allocatable :: fas, fas_out, freq, freq_out, ts

      integer npts, nfreq_out, nc_f_smc
      real tskip, tlength, ftaper
      real df_intrp, f_intrp_low, f_intrp_high
      real df_smooth, signnpw2
      real t_smc, t_fft

      character char_head(11)*80, comments(500)*80, cmnts2skip(500)*80
      real real_head(50)
      integer int_head(48)

      character buf*300, buf_upper*4,
     :        f_smc*300,
     :        f_fs2*300, f_sum*300, c2*2, c2up*2,
     :        ext_c*10, string*300, string4name*300,
     :        specify_frequencies*10, cl_arg*50

      logical freq_out_allocated, reset_smc_out

      character f_ctl*60, date*8,
     :        time_begin*10, time_start*10, time_stop*10

      character date_ctl_correct*8, date_ctl_in*30

      logical f_exist, f_smc_exist, dc_remove, datetime_l
      logical smc_format, log_spaced_f

!     Block commented to remove dependency on non standard library
!      f_exist = .false.
!      cl_arg = ' '
!      call getcl(cl_arg)
!      call trim_c(cl_arg, nc_cl_arg)
!      if (nc_cl_arg > 0) then
!        f_ctl = ' '
!        f_ctl = 'smc2fs2.ctl'
!        call trim_c(f_ctl, nc_f_ctl)
!        inquire(file=f_ctl(1:nc_f_ctl), exist=f_exist)
!          if (.not. f_exist) then
!            write(*,'(a)')
!     :             ' ******* CONTROL FILE '//f_ctl(1:nc_f_ctl)//
!     :             ' DOES NOT EXIST; QUITTING! ******* '
!            stop
!          end if
!      else
!        do while (.not. f_exist)
!          f_ctl = ' '
!          write(*, '(a)')
!     :      ' Enter name of control file '//
!     :      '(Enter=smc2fs2.ctl;'//
!     :      ' ctl-brk to quit): '
!          read(*, '(a)') f_ctl
!          if (f_ctl(1:4) .eq. '    ')
!     :      f_ctl = 'smc2fs2.ctl'
!          call trim_c(f_ctl, nc_f_ctl)
!          inquire(file=f_ctl(1:nc_f_ctl), exist=f_exist)
!          if (.not. f_exist) then
!            write(*,'(a)') ' ******* FILE DOES NOT EXIST ******* '
!          end if
!        end do
!      end if

!     Block added to set control file to smc2fs2.ctl
      f_exist = .false.
      f_ctl = 'smc2fs2.ctl'
      call trim_c(f_ctl, nc_f_ctl)
      inquire(file=f_ctl(1:nc_f_ctl), exist=f_exist)
        if (.not. f_exist) then
          write(*,'(a)')
     :             ' ******* CONTROL FILE '//f_ctl(1:nc_f_ctl)//
     :             ' DOES NOT EXIST; QUITTING! ******* '
          stop
        end if

      call get_lun(nu_ctl)
      open(unit=nu_ctl,file=f_ctl(1:nc_f_ctl),status='unknown')

      call skipcmnt(nu_ctl, cmnts2skip, nc_cmnts2skip)
      date_ctl_in = ' '
      read(nu_ctl,'(a)') date_ctl_in
      call trim_c(date_ctl_in,nc_date_ctl_in)

      date_ctl_correct = ' '
      date_ctl_correct = '03/10/10'
      call trim_c(date_ctl_correct,nc_date_ctl_correct)

      if (date_ctl_correct(1:nc_date_ctl_correct) .ne.
     :    date_ctl_in(1:nc_date_ctl_in)) then
        write(*,'(a)')
     :     ' The control file has the wrong date; update your '//
     :       'control file and rerun the program!'
        close(nu_ctl)
        stop
      end if

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)

      f_sum = ' '
      read(nu_ctl,'(a)') f_sum
      call trim_c(f_sum,nc_f_sum)
      call get_lun(nu_sum)
      open(unit=nu_sum,file=f_sum(1:nc_f_sum),status='unknown')

      write(nu_sum,'(a)') ' Output of program SMC2FS2:'
      date = ' '
      time_begin = ' '
! Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
      call DATE_AND_TIME( date, time_begin )
! Date is returned as 'CCYYMMDD'
! Time is returned as 'hhmmss.sss'
      write(nu_sum, '(2a)') ' Date: '//date(1:4)//'/'//
     :                                 date(5:6)//'/'//date(7:8)
      write(nu_sum, '(2a)') ' Time: '//time_begin(1:2)//':'//
     :           time_begin(3:4)//':'//time_begin(5:10)
      write(nu_sum,'(2a)') '   Control file = ', f_ctl(1:nc_f_ctl)

      freq_out_allocated = .false.

      loop_over_files: DO
        buf = ' '
        read(nu_ctl,'(a)',IOSTAT=iostatus) buf
        if (iostatus /= 0) EXIT
        call trim_c(buf,nc_buf)
        buf_upper = ' '
        buf_upper = buf(1:4)
        call upstr(buf_upper)
        if(nc_buf .eq. 0 .or. buf(1:1) .eq. '!') CYCLE
        if(buf_upper(1:4) .eq. 'STOP') EXIT

        processing_parameters: IF(buf_upper(1:3) .eq. 'PP:') THEN  ! a new set of processing parameters

! freq_out previously allocated in a specify-frequency set of processing
! parameters; deallocate for a new set of processing parameters.  Note
! that when the frequencies are not specified, freq_out allocated and
! deallocated for each file.
          IF (freq_out_allocated) THEN
            deallocate(freq_out)
          END IF

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
          read(nu_ctl,*) tskip, tlength_in
          write(nu_sum,'(a)')
     :     ' tskip, tlength_in ='
          write(nu_sum,*)
     :       tskip, tlength_in

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
          read(nu_ctl,*) dc_remove
          write(nu_sum,'(a,1x,l1)')  ' dc_remove = ', dc_remove

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
          read(nu_ctl,*) taper_front, taper_back
          write(nu_sum,'(a)')
     :     ' taper_front, taper_back ='
          write(nu_sum,*)
     :       taper_front, taper_back

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
          read(nu_ctl,*) signnpw2
          write(nu_sum,'(a)')
     :     ' signnpw2 ='
          write(nu_sum,*)
     :       signnpw2

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
          read(nu_ctl,*) itype, ipow, df_in, smooth_param
          write(nu_sum,'(a)')
     :     ' itype, ipow, df_in, smooth_param ='
          write(nu_sum,*)
     :     itype, ipow, df_in, smooth_param

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
          specify_frequencies = ' '
          read(nu_ctl,'(a)',IOSTAT=iostatus) specify_frequencies
          if (iostatus /= 0) EXIT
          call trim_c(specify_frequencies,nc_specify_frequencies)
          call upstr(specify_frequencies)

          IF (specify_frequencies(1:1) == 'Y') THEN

            CALL skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
            READ (nu_ctl,*) nfreq_out

            ALLOCATE ( freq_out(nfreq_out) )
            freq_out_allocated = .true.

            BACKSPACE nu_ctl   ! need to read nfreq_out twice, the first time for
                               ! allocation of freq_out

            READ (nu_ctl,*) nfreq_out, (freq_out(i), i = 1, nfreq_out)

            freq_param = 999.0  ! if 0.0, subroutine smooth_interpolate assumes
                              ! fft-spaced frequencies
            f_intrp_low = 0.0
            f_intrp_high = 0.0
            log_spaced_f = .false.
            WRITE (nu_sum,'(a)')
     :      ' nfreq_out, freq_out = '
            WRITE (nu_sum, '(1x,i3, 8(1x,es10.3))')
     :        nfreq_out, (freq_out(i), i = 1, nfreq_out)

          ELSE

            CALL skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
            READ (nu_ctl,*)
     :            f_intrp_low, f_intrp_high, rlogspacef, freq_param

            IF (freq_param == 0.0) THEN  ! no interpolation
              log_spaced_f = .false.
            ELSE
              IF (rlogspacef == 0.0) THEN
                 log_spaced_f = .false.
              ELSE
                 log_spaced_f = .true.   ! note: df_intrp is set to the number
                                         ! of frequencies later, in the
                                         ! frequency allocation statements, and in
                                         ! get_abs_fas.for
              END IF
            END IF

            WRITE (nu_sum,'(a)')
     :      ' f_intrp_low, f_intrp_high, log_spaced_f, freq_param = '
            WRITE (nu_sum, '(1x, 2(1x,es10.3), 1x,l1, 1x,es10.3)')
     :        f_intrp_low, f_intrp_high, log_spaced_f, freq_param

          END IF

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
          string4name = ' '
          read(nu_ctl,'(a)') string4name
          call trim_c(string4name, nc_string4name)
          write(nu_sum,'(a)')
     :   ' string4name = '//string4name(1:nc_string4name)

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
          buf = ' '
          read(nu_ctl,'(a)',IOSTAT=iostatus) buf
          if (iostatus /= 0) EXIT
          call trim_c(buf,nc_buf)
          buf_upper = ' '
          buf_upper = buf(1:1)
          call upstr(buf_upper(1:1))
          if (buf_upper(1:1) == 'Y') then
            if (log_spaced_f .or.
     :          specify_frequencies(1:1) == 'Y') then
              write(*,*)
              write(*,*) ' Asked for smc output format, but'//
     :          ' this is not allowed for log spaced frequencies'
              write(*,*) '  or specify_frequencies = Y; '//
     :                   'output will be reset to a column file.'
              smc_format = .false.
              reset_smc_out = .true.
            else
              smc_format = .true.
              reset_smc_out = .false.
            end if
          else
            smc_format = .false.
            reset_smc_out = .false.
          end if

          call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)

          CYCLE

        END IF processing_parameters


        f_smc = ' '
        f_smc = buf(1:nc_buf)
        call trim_c(f_smc, nc_f_smc)
!        write(nu_sum,'(a)') ' f_smc: '//f_smc(1:nc_f_smc)
!        write(nu_sum,'(a,1x,i3)') ' nc_f_smc = ', nc_f_smc
        inquire(file=f_smc(1:nc_f_smc), exist=f_smc_exist)
        if (.not. f_smc_exist) then
            write(*,'(a)') ' ****** FILE '//f_smc(1:nc_f_smc)//
     :          ' DOES NOT EXIST; SKIPPING!!! ******* '
          CYCLE
        end if

! Passed tests; proceed

! Extract the time series from the input file:

        print *
        print *,
     :   ' Processing file:', f_smc(1:nc_f_smc)
        write(nu_sum,*)
        write(nu_sum,'(1x,a)')
     :   'Processing file:'//f_smc(1:nc_f_smc)

        call CPU_Time(time_start_processing_record)

! Get number of points on this loop:

      call smc_npts(f_smc(1:nc_f_smc), npts, sps)

      write(nu_sum, '(3x,a, 2(1x,i5))') 'npts = ', npts

! Allocate time series:
      n4_alloc = npts
      ALLOCATE ( ts(n4_alloc) )

! Allocate frequency arrays:
      IF (specify_frequencies(1:1) == 'Y') THEN
        n4_alloc = nfreq_out
        ALLOCATE ( fas_out(n4_alloc) )  ! freq_out was allocated earlier.
      ELSE
        IF (freq_param == 0.0) THEN  ! use fft frequencies
          CALL get_npw2(npts,signnpw2,npw2)
          n4_alloc = npw2
        ELSE                       ! do not use fft frequencies
          IF (log_spaced_f) THEN   ! These statements could be put into the
                                   ! PP: section, as they do not depend on
                                   ! the characterstics of the smc file.
            n4_alloc = int(freq_param)
          ELSE
            df_intrp = freq_param
            n4_alloc = int((f_intrp_high-f_intrp_low)/df_intrp + 1.1)
          END IF
        END IF
        ALLOCATE ( freq_out(n4_alloc), fas_out(n4_alloc)  )
      END IF

! Now read the data from the files:

      tlength = tlength_in
      call SMCRead(f_smc(1:nc_f_smc), tskip, tlength,
     :               ts,
     :               char_head, int_head, real_head, comments)
      sps = real_head(2)
      npts_smc = int_head(17)
      dt = 1.0/sps
      t_smc = float(npts_smc)*dt

      nacc_start = tskip * sps
      nacc_stop  = nacc_start + tlength * sps

      write(nu_sum,'(3x,a,f8.6,a,i5,a,f8.2)')
     :     'signal_in: dt=', dt,
     :     ' npts_smc=', npts_smc,' t_smc=', t_smc
      write(nu_sum,'(3x,a,1x,f8.2)')
     :     'tlength after smcread = ', tlength

! Compute the Fourier spectra:

      call Get_Abs_FAS(
     : ts, sps, npts_smc, dc_remove,
     : taper_front, taper_back, signnpw2, nzpad, npw2, df_fft,
     : itype, ipow, df_in, smooth_param, df_smooth,
     : specify_frequencies,
     : freq_param, f_intrp_low, f_intrp_high, log_spaced_f,
     : freq_out, fas_out, nfreq_out, df_out)

      write(nu_sum,'(3x,a, 1x,i2. 1x,i6, 1x,i6, 1x,es11.4)')
     :      'signnpw2, nzpad, npw2, t4fft = ',
     :      int(signnpw2), nzpad, npw2, real(npw2)*dt

      deallocate( ts )

! Write output:

! Construct output filename:

      f_fs2 = ' '
      f_fs2 = f_smc(1:nc_f_smc)//
     :        string4name(1:nc_string4name)
      call trim_c(f_fs2,nc_f_fs2)

      if (reset_smc_out) then
        f_fs2 = f_fs2(1:nc_f_fs2)//'.reset2col'
        call trim_c(f_fs2,nc_f_fs2)
      end if

      print *,
     :   'Writing output to file:', f_fs2(1:nc_f_fs2)
      write(nu_sum,'(3x,a)') 'Output file = '//f_fs2(1:nc_f_fs2)

      string = ' '
      string = '| Output of program SMC2FS2:'
      call trim_c(string, nc_string)

      datetime_l = .true.

      if (smc_format) then
        call FS2Write(f_fs2, nc_f_fs2, fas_out,
     :              char_head, int_head, real_head, comments,
     :              f_smc, nc_f_smc,
     :              nacc_start, nacc_stop, npw2,
     :              nfreq_out, freq_out(1), freq_out(nfreq_out),
     :              df_out, df_smooth, ftaper,
     :              string,nc_string,datetime_l,
     :              tskip, tlength_in, df_intrp, f_intrp_low,
     :              f_intrp_high, signnpw2, dc_remove, nzpad)
      else
! Open output file:
        call get_lun(nu_fs2)
        open(unit=nu_fs2, file=f_fs2(1:nc_f_fs2), status='unknown')

! Write column heads:
        write(nu_fs2, '(a)') ' Output of SMC2FS2: FAS for file '//
     :                       f_smc(1:nc_f_smc)
        write(nu_fs2, '(a)') ' Nfreq_out = '
        write(nu_fs2, '(1x,i6)') nfreq_out
        write(nu_fs2,'(11x,a, 12x,a)') 'freq', 'fas'

        do j = 1, nfreq_out
          write(nu_fs2,'(1p,2(1x,e14.7))')
     :       freq_out(j), fas_out(j)
        end do

        close(nu_fs2)
      end if

      IF (specify_frequencies == 'Y') THEN
        deallocate (fas_out)    !  freq_out allocated in the PP: section
                                !  use it for all files
      ELSE
        deallocate (fas_out, freq_out)
        freq_out_allocated = .false.
      END IF

      call CPU_Time(time_stop_processing_record)

!DEBUG
      print *,' time_start_processing_record = ',
     :       time_start_processing_record
      print *,' time_stop_processing_record = ',
     :       time_stop_processing_record
!DEBUG

      write(nu_sum, '(3x,a,1x,1p, e10.3)')
     :       'Elapsed time (sec): ',
     :    time_stop_processing_record - time_start_processing_record
      write(nu_sum,*)
      write(nu_sum,*)

      END DO loop_over_files

! Close control and output files:

      close(nu_ctl)
      close(nu_sum)

      stop
      end
! ------------------------------------------------------------------- SMC2FS2


      include 'get_abs_fas.for'
      include 'skip.for'
      include 'skipcmnt.for'
      include 'get_lun.for'
      include 'upstr.for'
      include 'smc_npts.for'
      include 'smcread.for'
      include 'imnmax.for'
      include 'trim_c.for'
      include 'rc_subs.for'
      include 'locate.for'
      include 'yintrf.for'
      include 'fs2write.for'
      include 'get_date.for'
      include 'get_time.for'
      include 'time_diff.for'
