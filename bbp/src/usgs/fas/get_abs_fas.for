
! ---------------------------------------------------------------- Get_Abs_FAS
      subroutine Get_Abs_FAS(
     : ts, sps, npts_smc, dc_remove, 
     : taper_front, taper_back, signnpw2, nzpad, npw2, df_fft,  
     : itype, ipow, df_in, smooth_param, df_smooth,
     : specify_frequencies,
     : freq_param, f_intrp_low, f_intrp_high, log_spaced_f,
     : freq_out, fas_out, nfreq_out, df_out)
     
!INPUT:  ts, sps, npts_smc, dc_remove, 
!        taper_front, taper_back, signnpw2, 
!        itype, ipow, df_in, smooth_param,  
!        specify_frequencies,
!        freq_param, f_intrp_low, f_intrp_high, log_spaced_f,
!        freq_out, nfreq_out (if specify_frequencies = Y)

!OUPUT:  nzpad, npw2, df_fft, df_smooth,
!        freq_out, fas_out, nfreq_out, df_out


! NOTE: if specify_frequencies = Y, then nfreq_out, freq_out contains 
! the specified (input) number of frequencies and the specified frequencies.
 
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
!!   smooth_param = xi, which is the fraction of a decade spanned by the
!!                  first zero crossings for the
!!                  Konno and Ohmachi weighting function 
!!
!! ipow = power of FAS to be smoothed (2 = smoothing energy spectrum)
!!
!! df_smooth: Note: need df_smooth for linearly-spaced smoothers, 
!! and generally it should be the df from the fft.  The reason for
!! including it as an input parameter is to "fool" the
!! program to do smoothing over a specified number of points by
!! setting df = 1 and smooth_param = number of points (including 
!! points with zero weight at ends; e.g., smooth_param = 5 will 
!! give a smoother with weights 0, 1/4, 2/4, 1/4, 0; smooth_param
!! should be odd).
 !
!!smoothing: itype, ipow, df_smooth (0 = FFT df), smooth_param:
!  1  2 1.0 15.0
!!specify frequencies? (y/n):
!  N
!!frequency specification: 
!!  If specify_frequencies = Y, then enter the number 
!!    of frequencies followed by the specific frequency values.
!!  If specify_frequencies = N, then enter interpolation information 
!!    df_intrp(0=FFT freqs), f_intrp_low, f_intrp_high, 
!!    log-spaced(0=F,1=T; integer). NOTE: if log-spaced_f = T, then 
!!    interpret df_intrp as the number of frequencies
!   0.0 0.0 50.0 0  
!!character string to append to filename:
!   .fs
!!Output in smc format (Y,N)?
!! ***IMPORTANT NOTE: Output cannot be in smc format if use log-spaced frequencies
!! because programs such as smc2asc have not been modified to deal with log-spaced
!! frequency.

! Dates: 04/21/08 - Written by D.M. Boore, patterned after SMC2FS2
!        06/10/08 - Skip out of loop if reset nfreq_out because freq_out(j) .gt. fnyq 
!        02/09/10 - Replace prctnfrnttaper with length of tapers at front and back.
!                 - Reorder input slightly
!        02/10/10 - Replaced "go to 8888" with "exit"
!        02/15/10 - Allow specification of individual frequencies. 
!                 - Renamed 'freq' to 'freq_fft' and 'nfreq_fas' to 'nfreq_fft'.
!                 - Deleted code before allocate statement specifying
!                   mspct and n4fas_alloc
!        03/16/10 - Change df_intrp to freq_param in input argument list (because
!                   freq_param is not always the frequency spacing, it may be
!                   less confusing to use freq_param than df_intrp).
!        07/01/13 - Add deallocate statement
 
c
c Dimension and declaration statements:

      real fas_out(*), freq_out(*), ts(*) 
      real fas(:), freq_fft(:)
      allocatable :: fas, freq_fft
      
      character specify_frequencies*(*)

      integer npts, nfreq_out, nc_f_smc
      real tskip, tlength, ftaper
      real df_intrp, f_intrp_low, f_intrp_high
      real df_smooth, signnpw2
      real t_smc, t_fft, freq_param

      logical dc_remove, log_spaced_f

      call get_npw2(npts_smc,signnpw2,npw2)

      allocate( fas(npw2), freq_fft(npw2) )

! Now read the data from the files:

      dt = 1.0/sps
      t_smc = float(npts_smc)*dt
      
      fnyq = 0.5*sps

! Compute the Fourier spectra:

      CALL Abs_Spectra(ts, dt, npts_smc, dc_remove, 
     :  taper_front, taper_back, fas, npw2, signnpw2, 
     :  df_fft, nzpad)

      t_fft = float(npw2)*dt
        
      nfreq_fft = npw2/2 

! Set up the frequency array for output:

! First set up the frequency array for the spectra:

      DO j = 1, nfreq_fft
        freq_fft(j) = float(j-1)*df_fft
      END DO
 
      IF (specify_frequencies(1:1) == 'Y') THEN  ! If 'Y', then nfreq_out, freq_out 
                                            ! have been specified by the calling program
        df_out = -1.0
        m_start = -1
        m_stop  = -1
          
      ELSE  
      
        IF (freq_param == 0.0) THEN  ! use fft frequencies, between f_intrp_low and f_intrp_high

          m_start = 1
          call locate(freq_fft, nfreq_fft, f_intrp_low, m_start)
 
          call locate(freq_fft, nfreq_fft, f_intrp_high, m_stop)
          
          IF (m_stop > nfreq_fft) THEN
            m_stop = nfreq_fft
          END IF
 
          nfreq_out = m_stop - m_start + 1
          DO j = 1, nfreq_out
            freq_out(j) = freq_fft(j + m_start - 1)
          END DO
          df_out = df_fft
        
        ELSE                       ! use interpolated frequencies
      
          IF (log_spaced_f) THEN
            nfreq_out = int(freq_param)
            dlogf = alog10(f_intrp_high/f_intrp_low)/
     :              float(nfreq_out-1)
            DO j = 1, nfreq_out
              freq_out(j) = f_intrp_low*10.0**(float(j-1) * dlogf)
            END DO 
            df_out = -1.0
            m_start = -1
            m_stop = -1
          ELSE
            df_intrp = freq_param
            nfreq_out = int((f_intrp_high-f_intrp_low)/df_intrp + 1.1)
            DO j = 1, nfreq_out
              freq_out(j) = float(j-1) * df_intrp + f_intrp_low
            END DO
            df_out = df_intrp
            m_start = -1
            m_stop = -1
          END IF
        
        END IF
        
      END IF   

      DO j = 1, nfreq_out
        IF(freq_out(j) .gt. fnyq) THEN
          nfreq_out = j -1
          print *, ' ****ATTENTION********* '
          print *, '   nfreq_out has been reset so that'//
     :            ' it is less than the Nyquist frequency'
          print *, ' ****ATTENTION********* '
          EXIT
        END IF
      END DO

! Smooth and interpolate the spectra

      IF (df_in .eq. 0.0) THEN
        df_smooth = df_fft
      ELSE
        df_smooth = df_in
      END IF
           
      CALL Smooth_interpolate (freq_fft, fas, nfreq_fft, 
     :                         freq_out, fas_out, nfreq_out, 
     :                         itype, ipow, df_smooth, smooth_param,
     :                         freq_param, m_start, m_stop)

      deallocate( fas, freq_fft )
 
      RETURN
      END
      
! ---------------------------------------------------------------- Get_Abs_FAS


      include 'get_npw2.for'
      include 'fork.for'
      include 'tprfrctn.for'
      include 'dcdt.for'
      include 'zeropad.for'
      include 'abs_spectra.for'
      include 'smooth_interpolate.for'
 
