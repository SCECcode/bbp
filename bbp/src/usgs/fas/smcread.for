
!------------------   Begin SMCRead   ----------------------------------
      subroutine SMCRead(f_in,  
     :                   tskip, tlength,
     :                   y, char_head, int_head, real_head, comments)

!  If tskip = 0.0, no data are skipped
!  If tlength = 0.0, only headers are read
!  If tlength < 0.0, read all data, less tskip
!  Otherwise tlength of data is read, after skipping tskip, unless 
!  tskip + tlength exceeds the length of the input time series, in which
!  case tlength is adjusted so that all data, less tskip, are read.  Thus
!  setting tlength to a large number (or setting it to a negative number)
!  guarantees that all data (minus tskip) are read.

! The program requires get_lun.for

! The output arrays char_head, int_head, real_head should be dimensioned
! as follows in the calling program:

!      character char_head(11)*80, comments(500)*80 ! "80" is the critical parameter
!      integer int_head(48)
!      real real_head(50)

!  Dates: 07/17/97 - Modified order of input parameters, to put input
!                    first (unit, file, tskip), followed by output.
!         12/18/97 - Read tlength
!         05/19/98 - Determine unit from within program; remove "unit"
!                    from argument list
!         03/09/99 - Change logic to determine amount of record to
!                    read and allow reading of all record, minus tskip,
!                    by setting tlength < 0.0.
!         05/18/99 - Deleted npts, sps from input parameter list.  Added 
!                    character array with comments string (allow for 
!                    50 comments).
!         12/16/99 - Renamed from ReadSMC
!         05/10/00 - Replaced
!                        tlength = float(npts_out/sps)
!                    with
!                        tlength = float(npts_out)/sps
!                    (the error and correction was pointed out by Wang Guoquan
!                    in email on 05/06/2000)
!         01/04/01 - Check for unevenly spaced data (sps = rnull) and
!                    return 0.0 in y
!         02/07/01 - If tlength = 0, only read headers
!         12/14/01 - Fixed logic for case when tlength < 0 and tskip > 0.
!                    Also, program used to reset tlength, but it does not do
!                    that now.  The user can compute tlength as 
!                    float(int_head(17))/real_head(2)
!         01/01/05 - Use "int" in calculation of npts (do not use
!                    implicit conversion from real to integer)
!         10/15/07 - Use higher precision input format if int_head(47) = 8
!         04/04/10 - Add call to trim_c, and use nc_f_in.
!         10/20/15 - Changed format of high precision input.
!                  - This was not necessary, as the old format statement can read
!                    either e14.7 or e14.7e1 (the new format being used in smcwrite for
!                    high precision).
!         12/30/15 - Replaced "int" with "nint" when computing nskip and npts2read
!                  - A few modernizations of syntax.

      real real_head(*)
      integer int_head(*)
      character*80 char_head(*), comments(*)
      character f_in*(*) 
      real y(*)

      
      
      call get_lun(nu)
      call trim_c(f_in, nc_f_in)
      open(unit= nu, file=f_in(1:nc_f_in), status='unknown')

      do i = 1, 11
        read(nu, '(a)') char_head(i)
      end do

      read(nu, '(8I10)') (int_head(i), i=1, 48)

      read(nu, '(5e15.7)') (real_head(i), i = 1, 50)

      sps = real_head(2)

      do i = 1, int_head(16)
        read(nu,'(a)') comments(i)
      end do

      if(real_head(2) == real_head(1)) then   ! unevenly sampled data

         int_head(17) = 2
         do i = 1, int_head(17)
           y(i) = 0.0
         end do

      else

! Skip into the trace:

        nskip = nint(tskip * sps)
!        nskip = int(tskip * sps)

! How much record to read?

        npts_in = int_head(17)
        npts2read = nint(tlength * sps)
!        npts2read = int(tlength * sps)

        if (tlength < 0.0) then  
          npts_out = npts_in - nskip
        else if (nskip + npts2read <= npts_in) then
          npts_out = npts2read
        else
          npts_out = npts_in - nskip
!          tlength = float(npts_out)/sps
        end if

        int_head(17) = npts_out
      
        if (tlength == 0.0) then
          return    ! only read headers
        else
          if (int_head(47) == 8) then
            read(nu, '(5(e14.7))') (y(i), i = 1, npts_out + nskip)
!            read(nu, '(5(e14.7e1))') (y(i), i = 1, npts_out + nskip)
          else
            read(nu, '(8(1pe10.4e1))') (y(i), i = 1, npts_out + nskip)
          end if
          do i = 1, npts_out
            y(i) = y( i + nskip)
          end do
        end if

      end if

      close(unit=nu)

      return
      end
!------------------   End SMCRead   ----------------------------------




