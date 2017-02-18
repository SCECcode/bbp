* -------------------------- BEGIN SMOOTH_interpolate --------------------------
      subroutine Smooth_interpolate (x_in, y_in, npts_in, 
     :                               x_out, y_out, npts_out,  
     :                   itype, ipow, df_smooth, smooth_param,
     :                   freq_param, m_start,m_stop)

! Meaning of smoothing input parameters
!
! NO SMOOTHING
! itype = 0
! SMOOTHING OVER EQUALLY SPACED FREQUENCIES
! itype = 1: box weighting function
!   smooth_param = width of box weighting function (Hz)
! itype = 2: triangular weighting function
!   smooth_param = width of triangular weighting function (Hz)
! SMOOTHING OVER LOGARITHMICALLY SPACED FREQUENCIES
! itype = 3: box weighting function
!   smooth_param = xi, which is the fraction of a decade for the
!                  box weighting function 
! itype = 4: triangular weighting function
!   smooth_param = xi, which is the fraction of a decade for the
!                  triangular weighting function 
! itype = 5: Konno and Ohmachi weighting function (see BSSA 88, 228-241)
!   smooth_param = xi, which is the fraction of a decade for which
!                  the Konno and Ohmachi weighting function is greater
!                  than 0.043.(it is related to
!                  their smoothing parameter b by the equation
!                  b = 4.0/smooth_param, so we have this correspondence between
!                  b and smooth_param
!                         b smooth_param 
!                        10         0.40
!                        20         0.20
!                        40         0.10
!                  
!                  b = 40 seems to be commonly used, but I do not think that it gives
!                  enough smoothing; I prefer smooth_param = 0.2, corresponding to
!                  b = 20. 
!
! ipow = power of FAS to be smoothed (2 = smoothing energy spectrum)
!
! df_smooth: Note: need df_smooth for linearly-spaced smoothers, 
! and generally it should be the df from the fft.  For general x data, it is
! the spacing between x values, assumed to be constant,  The reason for
! including it as an input parameter is to "fool" the
! program to do smoothing over a specified number of points by
! setting df_smooth = 1 and smooth_param = number of points (including 
! points with zero weight at ends; e.g., smooth_param = 5 will 
! give a smoother with weights 0, 1/4, 2/4, 1/4, 0; smooth_param
! should be odd).
!
 
* Dates:  03/27/08 - Written by D. Boore, patterned after smooth.for
*         06/26/09 - Restrict i/ic so that w>=0.1 (rather than go to
*                    first zero crossing of sin X/X).
*         07/11/09 - With this definition (w>=0.1), redefine the meaning of smooth_param
*                    so that it is closer to the fraction of a decade over which the
*                    smoothing takes place.  The width in log frequency corresponding
*                    to the weighting function being greater than 0.1
*                    is related to b by Delta logf = (1/b)*3.5.  Thus
*                    if smooth_param = Delta logf, we have b=3.5/smooth_param
*                    instead of b=2*pi/smooth_param.
*         07/11/09 - I decided to let the smoothing width be determined by 
*                    wmin = 0.043 (as a compromise between 0.01 and 0.1, and because
*                    this corresponds to log(arg)=+-2, a nice round number), 
*                    which means that b=4.0/smooth_param.
!         03/15/10 - Correct error for the case of smoothing and fft frequencies
!                    starting at a frequency greater than 0.0
!        03/16/10 - Change df_intrp to freq_param in input argument list (because
!                   freq_param is not always the frequency spacing, it may be
!                   less confusing to use freq_param than df_intrp).
  
      real x_in(*), y_in(*), x_out(*), y_out(*)
      real work(:)
      allocatable :: work
      
      pi = 4.0*atan(1.0)

      allocate(work(npts_in))
      
      if (itype == 1 .or. itype == 2) then  ! calculate range of indices to use in average
        ic1 = int(aint(0.5*smooth_param/df_smooth))
        ic2 = npts_in - ic1
        b = 0.0
      else if (itype == 5) then  ! convert smoothing parameter to K&O parameters
!        b = 2.0*pi/smooth_param  ! corresponds to zero crossings of K&O weighting function
!        b = 4.6370/smooth_param      ! corresponds to K&O weighting function = 0.01
        b = 4.0/smooth_param      ! corresponds to K&O weighting function = 0.043
!        b = 3.4996/smooth_param      ! corresponds to K&O weighting function = 0.1
        ic1 = 0
        ic2 = 0
      end if
        
      do i = 1, npts_in
        work(i) = y_in(i)**float(ipow)
      end do
      

      if (itype == 0 .and. freq_param == 0.0) then  ! no smoothing, no interpolation
        do j = 1, m_stop - m_start + 1
          y_out(j) = y_in(j+m_start - 1)
        end do
        return
      else if (itype == 0 .and. freq_param /= 0.0) then ! no smoothing, interpolation
        do j = 1, npts_out
          y_out(j) = yintrf(x_out(j), x_in, y_in,
     :                 npts_in)
        end do
        return
      else if (itype /= 0 .and. freq_param == 0.0) then ! smoothing, no interpolation
!        do j = 1, npts_in
        do j = 1, m_stop - m_start + 1
          j4smooth = j + m_start - 1
          y_out(j) = smooth_f(work, npts_in, j4smooth, ic1, ic2, 
     :                        b, itype, smooth_param) 
          y_out(j) = y_out(j)**(1.0/float(ipow))  ! undo power
        end do         
        return
      end if
      
* What's left?  Smoothing and interpolation

      do iy = 1, npts_out
      
* For each x_out, calculate indices of x_in:

        call locate(x_in, npts_in, x_out(iy), j)  ! x_out(ic) will be between x_in(j) and x_in(j+1)
                                                  ! so need to compute smoothed value of y at
                                                  ! these two points and then interpolate.
                                                  
* Compute smoothed values at these bracketing values of x:

        y_smooth_j = 
     :   smooth_f(work, npts_in, j, ic1, ic2, b, itype, smooth_param)                                                 
        y_smooth_jp1 = 
     :   smooth_f(work, npts_in, j+1, ic1, ic2, b, itype, smooth_param)                                                 
        
* And interpolate the smoothed values

        y_out(iy) = y_smooth_j + 
     :             (x_out(iy)-x_in(j))*(y_smooth_jp1 - y_smooth_j)/
     :             (x_in(j+1)-x_in(j))
     
        y_out(iy) = y_out(iy)**(1.0/float(ipow))  ! undo power
     
      end do

      deallocate(work)

      return
      end
* -------------------------- END Smooth_interpolate --------------------------

* -------------------------- BEGIN smooth_f --------------------------        
      real function smooth_f(work, npts_in, ic, ic1, ic2, b, 
     :                  itype, smooth_param)
      
      implicit none

      real, intent(in) :: b, smooth_param, work(*)
      integer, intent(in) :: npts_in, ic, ic1, ic2, itype
      
      real :: sum, wavg, weight, f_fc
      real :: w_triangular, w_logtriangular, w_konno_ohmachi
      integer :: iz1, iz2, j
      
      if (itype == 1 .or. itype == 2) then
          iz1 = ic - ic1
          iz2 = ic + ic1
      else if (itype == 5) then
!          f_fc = 0.00480**(1.0/b)    ! corresponds to K&O weighting function = 0.01
          f_fc = 0.0100**(1.0/b)    ! corresponds to K&O weighting function = 0.043
!          f_fc = 0.01779**(1.0/b)    ! corresponds to K&O weighting function = 0.1
          iz1 = nint(float(ic)*f_fc)      
          iz2 = nint(float(ic)/f_fc)      
          
      else
          iz1 = nint(float(ic)*10.0**(-0.5*smooth_param))
          iz2 = nint(float(ic)*10.0**(0.5*smooth_param))
      end if
        
        if (iz1 < 1) then
          iz1 = 1
        end if
        if (iz2 > npts_in) then
          iz2 = npts_in
        end if
        
        sum = 0.0
        wavg = 0.0
        
        do j = iz1, iz2
          if (itype == 1) then
            weight = 1.0
          else if (itype == 2) then
            weight = w_triangular(j, ic, iz1, iz2, smooth_param)
          else if (itype == 3) then
            weight = 1.0
          else if (itype == 4) then
            weight = w_logtriangular(j, ic, iz1, iz2, smooth_param)
          else if (itype == 5) then          
            weight = w_konno_ohmachi(j, ic, b)
          else
            print *, ' error, invalid itype = ', itype, '; quitting'
            stop
          end if
          
          wavg = wavg + weight          
          sum = sum + weight * work(j)
        end do
        
        smooth_f = (sum/wavg) 
        
        end function smooth_f
* -------------------------- END smooth_f --------------------------        
        
* -------------------------- BEGIN w_triangular --------------------------
      function w_triangular(j, ic, iz1, iz2, smooth_param)
      
      if (j .eq. ic) then
        w_triangular = 1.0
      else if (j .lt. ic) then          
        w_triangular = 
     :      (float(j) - float(iz1))/(float(ic) - float(iz1)) 
      else
        w_triangular = 1 - 
     :      (float(j) - float(ic))/(float(iz2) - float(ic)) 
      end if
        
      return
      
      end
* -------------------------- END w_triangular --------------------------

* -------------------------- BEGIN w_logtriangular --------------------------
      function w_logtriangular(j, ic, iz1, iz2, smooth_param)
      
      if (j .eq. ic) then
        w_logtriangular = 1.0
      else if (j .lt. ic) then          
        w_logtriangular = 
     :      alog10(float(j)/float(iz1))/alog10(float(ic)/float(iz1)) 
      else
        w_logtriangular = 1 - 
     :      alog10(float(j)/float(ic))/alog10(float(iz2)/float(ic)) 
      end if
        
      return
      
      end
* -------------------------- END w_logtriangular --------------------------

* -------------------------- BEGIN w_konno_ohmachi --------------------------
      function w_konno_ohmachi(j,ic,b)
      
      if (j .eq. ic) then
        w_konno_ohmachi = 1.0
      else
        arg = b*alog10(float(j)/float(ic))
        w_konno_ohmachi = (sin(arg)/arg)**4.0    
      end if
      
      return
      
      end
* -------------------------- END w_konno_ohmachi --------------------------
      

      

