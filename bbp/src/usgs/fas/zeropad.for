
! ------------------------------ BEGIN ZEROPAD --------------------
      subroutine ZEROPAD (Y,NIN,NPW2)

c Pads time-series array Y with (NPW2-NIN) zeroes.
*  With this program 
c the window of the data, which determines NIN, can be different 
c for different time series, yet the overall length of 
c the time series used in the FFT can 
c be the same (call get_npw2 with an overall duration NTOTIN), thus 
* guaranteeing 
c that the frequencies for which FFT
c values are computed are the same.
c I assume that the user makes sure that NTOTIN .ge. NIN.

* Dates: 12/12/00 - Written by D. Boore, patterned after zeropad2.for

      real Y(*)
      
      if (npw2 .le. nin) then
        return
      else
        do i=nin+1,npw2
          y(i) = 0.0
        end do
        return
      end if

      end
! ------------------------------ END ZEROPAD --------------------



