*------------------   Begin FS2Write  ----------------------------------
      subroutine FS2Write(f_fs2, nc_f_fs2, fas, 
     :                    char_head, int_head, real_head, comments,
     :                    f_smc, nc_f_smc, 
     :                    nacc_start, nacc_stop, npwr2,
     :                    nfreq, flow, fhigh, df, df_smooth, ftaper,
     :                    string,nc_string,datetime_l,
     :                    tskip, tlength_in, df_intrp, f_intrp_low,
     :                    f_intrp_high, signnpw2, dc_remove, nzpad)

*  Writes Fourier amplitude spectra into SMC format

* If datetime_l = .true., write date/time into the comments
* If nc_string .ne. 0, write "string" into the comments

* The arrays char_head, int_head, real_head must be dimensioned
* by the calling program as follows:
*
*      character char_head(11)*80, comments(N)*80 ! see "N" below
*      integer int_head(48)
*      real real_head(50)

*
* default values for int_head = -32768
*                    real_head = 1.7e+38

* Dates: 06/15/01 - written by D. Boore, patterned after smcwrite
*        07/10/01 - Added more processing parameters to comments in output file
*        11/11/01 - Write flow, df, fhigh, nfreq, npw2, dt to comments
*        12/18/01 - Store 0.0 if value is less than what can be written in the
*                   e format.
*        06/25/09 - Use get_date, get_time rather than datetime
!        12/29/10 - Include more precision in the output freq and fas.
!        08/01/12 - I reset int_head(47) to 8 (high precision) AFTER
!                   writing out int_head!  I corrected this.

      real real_head(*)
      integer int_head(*)
      character*80 char_head(*), comments(*)
      character f_fs2*(*), f_smc*(*), string*(*)
      real fas(*)
      logical datetime_l, dc_remove
      character date_c*10, time_c*11

* Open output file:
      call get_lun(nu_out)
      open(unit=nu_out, file=f_fs2(1:nc_f_fs2), status='unknown')
      
* Write headers:

* First write the 11 lines of text headers:

      char_head(1) = ' '
      char_head(1) = '6 FOURIER AMPLITUDE SPECTRUM'
      
      do i = 1, 11
        write(nu_out,'(a)') char_head(i)
      end do
          
* Fill in the int and real headers:

      int_head(42) = nacc_start
      int_head(43) = nacc_stop
      int_head(44) = npwr2
      int_head(45) = nfreq

      real_head(42) = flow
      real_head(43) = fhigh
      real_head(44) = df
      real_head(45) = df_smooth
      real_head(46) = ftaper

* Write to comments:

      ncomments = int_head(16)

      ncomments = ncomments + 1
      comments(ncomments) = '|'
      ncomments = ncomments + 1
      comments(ncomments) = 
     :  '| Fourier amplitude spectrum computed for file:'
      ncomments = ncomments + 1
      comments(ncomments) = '|   '//f_smc(1:nc_f_smc)

      if(nc_string .ne. 0) then
        ncomments = ncomments + 1
        comments(ncomments) = string(1:nc_string)
      end if

      if(datetime_l) then
        iflag = 0
        call get_date(date_c)
        call get_time(time_c)
        ncomments = ncomments + 1
        comments(ncomments) = 
     :  '|   file written on '//date_c//' at time '//time_c
      end if

      ncomments = ncomments + 1
      comments(ncomments) = 
     :  '|   tskip, tlength_in, ftaper = '
      write(comments(ncomments)(42:80),
     :  '(1x,f7.2,1x,f8.2,1x,f4.1)') 
     :    tskip, tlength_in, ftaper
      ncomments = ncomments + 1
      comments(ncomments) = 
     :  '|   df_intrp(0=FFT freqs), f_intrp_low, f_intrp_high = '
      write(comments(ncomments)(60:80),
     :  '(1x,f3.1,1x,f7.3,1x,f7.3)') 
     :       df_intrp, f_intrp_low, f_intrp_high

      ncomments = ncomments + 1
      comments(ncomments) = 
!        1234567890123456789012345678901234567890
     :  '|   df_smooth(Hz;0=none), '//
     :  'signnpw2(<0= backup for npw2, no zpad) = '
      write(comments(ncomments)(66:80),
     :  '(1x,f5.2,1x,f4.1)') 
     :    df_smooth, signnpw2

      ncomments = ncomments + 1
      comments(ncomments) = 
!        1234567890123456789012345678901234567890
     :  '|   npw2, nzpad, dt = '
      write(comments(ncomments)(23:80),
     :  '(1x,i6,1x,i5, 1x,f8.5)') 
     :       int_head(44), nzpad, 1.0/real_head(2)

      ncomments = ncomments + 1
      comments(ncomments) = 
     :  '|   dc_remove? = '
      write(comments(ncomments)(20:80),
     :  '(1x,l1)') 
     :    dc_remove

      ncomments = ncomments + 1
      comments(ncomments) = 
     :  '|   flow (rhead(42)), df(rhead(44)), fhigh(rhead(43)), '//
     :  'nfreq(ihead(45)) = '
      ncomments = ncomments + 1
      comments(ncomments)(1:1) = '|' 
      write(comments(ncomments)(2:80),
     :  '(6x,f8.4,9x,f8.6,7x,f8.4,13x,i5)') 
     :    real_head(42), real_head(44), real_head(43), int_head(45)

      ncomments = ncomments + 1
      comments(ncomments) = '|'

      int_head(16) = ncomments

      int_head(47) = 8

      write(nu_out, '(8I10)') (int_head(i), i=1, 48)

      write(nu_out, '(5e15.7)') (real_head(i), i = 1, 50)

      if (int_head(16) .gt. 0) then
        do i = 1, int_head(16)
          write(nu_out,'(a)') comments(i)
        end do
      end if
          
* Check for a small number.  There is probably a better way of doing this
* involving writing one line at a time rather than replacing a value of FAS
      do i = 1, int_head(45)
        if (abs(fas(i)) .le. 1.0e-09) then
          write(*,'(a)') '  replacing fas with 0.0; i, fas = '
          write(*,*) i, fas(i)
          fas(i) = 0.0
        end if
      end do
      
      write(nu_out, '(5(1pe14.7e1))') (fas(i), i = 1, int_head(45))

* That should be it; close the output file.

      close(unit=nu_out)

      return
      end
*------------------   End FS2Write   ----------------------------------

 

