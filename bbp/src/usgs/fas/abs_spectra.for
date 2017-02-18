* ----------------------- BEGIN Abs_Spectra -------------------
      subroutine Abs_Spectra( datain, dt, npts, dc_remove, 
     :    taper_front, taper_back, spect, npw2, signnpw2, df_fft, 
     :    nzpad)

c Returns the absolute value of the spectrum.
c The program applies a tapered window to the
c front and back of the time series, pads with zeros,  and computes the
c spectrum.

c Written by D. M. Boore

* Dates:  3/22/89 - Modified from SqrSpectFT, by removing
*                   the smoothing and interpolation.
*        11/08/89 - Pass df_fft out through the argument list.
*         5/02/95 - Changed iprcntfrtaper to iprcnttaper
*         4/07/95 - Replaced prcntfrtaper as iprcntfrtaper (an integer).  
*                   In some previous versions read in prcntfrtaper but 
*                   intended it to be an integer.  I wonder
*                   how many bad results were obtained as a result
*                   of this error?!
*        02/10/99 - Minor cleanup, appended called routines to end
*        12/12/00 - Replaced iprctnfrtaper with frctnfrtaper (and subroutine
*                   fbctpr with tprfrctn).  Also replaced zeropad2 with 
*                   zeropad, which accepts npw2 as in input, rather than an
*                   output argument.  I had to do this in order to use 
*                   allocatable arrays.  I introduced signnpw2, which controls
*                   whether npw2 is less than npts (if signnpw2<0) or greater
*                   than npts, in which case zeros are added.
*        02/01/01 - Base allocation on larger of npts or npw2, and use npw2 
*                   rather than npts in conditioning routines if npw2 < npts.
*        07/10/01 - Add nzpad to output parameters
*        11/10/01 - Correct calculation of nzpad
*        12/18/01 - Uses double precision version of FORK (FORKDP)
*        12/18/01 - Removed double precision version of Fork; it did not seem 
*                   to make a difference.
*        02/09/10 - Replace prctnfrnttaper with length of tapers at front and back.
*                 - Reorder input slightly

      real DATAIN(*), SPECT(*)
      real DATA(:)
      complex CX(:)
      allocatable :: data, cx
      logical dc_remove, detrend

      call get_npw2(npts,signnpw2,npw2)

      n_alloc = max0(npw2, npts)
      allocate( data(n_alloc), cx(n_alloc))

c Fill working array with input data:

      do i = 1, npw2
        data(i) = 0.0
      end do

      n_condition = min0(npw2, npts)

      do i = 1, n_condition
        data(i) = datain(i)
      end do

c
c     remove dc, apply taper, pad with zeros
c
      detrend = .false.
      call dcdt(data, dt, n_condition, 1, n_condition, 
     : dc_remove, detrend)

      data_length = dt * float(n_condition)
      frctnfront = taper_front/data_length
      frctnback =  taper_back/data_length
      call tprfrctn (frctnfront, frctnback, data, n_condition)

      if (npw2 .gt. npts) then
        call zeropad( data, npts, npw2)
        nzpad = npw2 - npts
      else
        nzpad = 0
      end if

c sample rate and frequency spacing:

      sr = 1.0 / dt
      df_fft = sr/ float(npw2)

      npw2d2 = npw2 / 2

c Fill working complex array

      do j = 1, npw2
        cx(j)=cmplx(data(j), 0.0)
      end do

c     FFT to get spectrum

*      call forkdp(npw2,cx,-1.)
      call fork(npw2,cx,-1.)

      do j =1 ,npw2d2
        cx(j)=cx(j)*dt*sqrt(float(npw2))
        spect(j)=cabs( cx(j) )
      end do

      deallocate( data, cx)

      return
      end
* ----------------------- END Abs_Spectra -------------------

