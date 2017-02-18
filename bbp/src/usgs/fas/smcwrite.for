
!------------------   Begin SMCWrite  ----------------------------------
      subroutine SMCWrite(f_out, precision,
     :                    y, char_head, itype,
     :                    int_head, real_head, comments)

!  Reformats a time series
!  into SMC format (the format used on the CD-ROM)

! If itype < 0, the first text header is given by char_head(1) and
!               not set to one of the standard headers
!    itype = 0 for 'UNKNOWN'
!    itype = 1 for 'UNCORRECTED ACCELEROGRAM'
!    itype = 2 for 'CORRECTED ACCELEROGRAM'
!    itype = 3 for 'VELOCITY'
!    itype = 4 for 'DISPLACEMENT'
!    itype = 5 for 'RESPONSE SPECTRA'

! The program requires:
!        get_lun.for
!        mnmaxidx.for

! The arrays char_head, int_head, real_head must be dimensioned
! by the calling program as follows:
!
!      character char_head(11)*80, comments(N)*80 ! see "N" below
!      integer int_head(48)
!      real real_head(50)

! Note that
!      int_head(17) = npts
!      real_head(2) = sps
!      int_head(16) = ncomments (dimension array "comments" accordingly
!                                so that N is a number > = ncomments)
!      I have found that if no comments are included, BAP cannot process
!      the resulting file.  For this reason, be sure that at least one comment
!      is included (ncomments = 1, comments(1) = '|').
!
! default values for int_head = -32768
!                    real_head = 1.7e+38

! The minimum required input for char_head is to set each of the 11 variables
! to "*" as in:
!  do i = 1, 11
!    char_head(i) = '*'
!  end do


! Dates: 05/13/99 - revised by D. Boore
!        05/18/99 - continued, adding the comments string array
!        12/15/99 - Renamed from WriteSMC to SMCWrite, and added
!                   itype to argument list and find and write peak
!                   accelerations (itype < 0 to echo input char_head(1)).
!        01/29/00 - Added comments regarding dimensioning of input variables
!        02/18/00 - Substitute 0.0 if abs(y) < 1.0e-09 (otherwise, the field
!                   is stored as ********, because an exponent of -10 is
!                   larger than allowed by the format statement).
!        02/25/00 - Changed output format of peak motion in text headers from
!                   1pe10.3 to 1pe9.3e1
!        02/27/01 - Set char_head(7)(42:60) to blanks before filling with
!                   peak motion
!        08/22/01 - Added comment about minimum requirement for char_head
!        10/15/07 - Use high precision format if precision=HIGH
!        04/24/09 - Changed format of read header output from 5e15.7 to 1p5e15.7
!        04/04/10 - Add call to trim_c, and use nc_f_out.
!        02/02/12 - Compute r_ep if eq and sta coords /= null
!        06/24/13 - Replaced ".eq." with "==", etc.
!                 - Replaced 1pex.x with esx.x
!                 - In writing real_head, revert to the format 5e15.7 rather than 5es15.7,
!                   as the latter might imply more precision than available in the
!                   single precision numbers.
!        10/20/15 - Change high precision output format to include blanks between
!                   numbers (as suggested by Chris Stephens)
!        10/24/00 - For both output precision options substitute 0.0 if abs(y) < 1.0e-09 (otherwise, the field
!                   is stored as ********, because an exponent of -10 is
!                   larger than allowed by the format statement).

      real real_head(*)
      integer int_head(*)
      character*80 char_head(*), comments(*)
      character f_out*(*), chead1(6)*26, precision*(*)
      real y(*)
      data chead1/
     : '0 UNKNOWN                 ', '1 UNCORRECTED ACCELEROGRAM',
     : '2 CORRECTED ACCELEROGRAM  ', '3 VELOCITY                ',
     : '4 DISPLACEMENT            ', '5 RESPONSE SPECTRA        '/

! Open output file:
      call get_lun(nu_out)
      call trim_c(f_out, nc_f_out)
      open(unit=nu_out, file=f_out(1:nc_f_out), status='unknown')

! Write name to screen:
      print *,
     :  ' nu_out, Output file name = ', nu_out, ', ', f_out(1:nc_f_out)

! Get peak motions:

      call mnmaxidx(y, 1, int_head(17), 1,
     :             ymin, ymax, indx_min, indx_max)
      peak_mtn = ymax

      if (abs(ymin) > abs(ymax)) peak_mtn = ymin

! Compute r_ep if eq, sta coords /= null:  NOTE: Need to include distaz in the calling program

      rnull = real_head(1)
      if (real_head(3) /= rnull .and. real_head(3) /= rnull) then  ! eq coords not null
        if (real_head(11) /= rnull .and. real_head(12) /= rnull) then  ! sta coords not null
          eqlat = real_head(3)
          eqlong = real_head(4)
          stalat = real_head(11)
          stalong = real_head(12)
          call distaz( -1.0, eqlat,eqlong,stalat,stalong,
     :    rdeg, rkm, az, baz)
          real_head(17) = rkm
        end if
      end if


! Write headers:

! First write the 11 lines of comments:

      if (itype < 0) then
        write(nu_out,'(a)') char_head(1)
      else
        write(nu_out,'(a)') chead1(itype+1)
      end if

      char_head(7)(34:41) = 'pk mtn ='    ! makes more sense than 'pk acc ='
      char_head(7)(42:60) = ' '
      write(char_head(7)(42:50),'(1pe9.3e1)') peak_mtn

      do i = 2, 11
        write(nu_out,'(a)') char_head(i)
      end do

      if (precision(1:4) == 'HIGH') then
        int_head(47) = 8
      else
        int_head(47) = int_head(1) ! null value
      end if

      write(nu_out, '(8I10)') (int_head(i), i=1, 48)

      dt = 1.0/real_head(2)
      real_head(29) = float(indx_max - 1) * dt
      real_head(30) = ymax
      real_head(31) = float(indx_min - 1) * dt
      real_head(32) = ymin

!      write(nu_out, '(5es15.7)') (real_head(i), i = 1, 50)
!DEBUG
!      write(nu_out, '(a, 1x,e15.7, 1x,f7.2, 1x,es12.5)')
!     :   ' In SMCWRITE: real_head(2), dt = ',
!     :     real_head(2), real_head(2), dt
!DEBUG
      write(nu_out, '(5e15.7)') (real_head(i), i = 1, 50)


      if (int_head(16) > 0) then
        do i = 1, int_head(16)
          write(nu_out,'(a)') comments(i)
        end do
      end if

        do i = 1, int_head(17)
          if (abs(y(i)) <= 1.0e-09) y(i) = 0.0
        end do

      if (precision(1:4) == 'HIGH') then
        write(nu_out, '(5(es14.7e1))') (y(i), i = 1, int_head(17))
      else
!        do i = 1, int_head(17)
!          if (abs(y(i)) < 1.0e-09) y(i) = 0.0
!        end do
        write(nu_out, '(8(es10.4e1))') (y(i), i = 1, int_head(17))
      end if


! That should be it; close the output file.

      close(unit=nu_out)

      return
      end
!------------------   End SMCWrite   ----------------------------------
