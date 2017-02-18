      Program Asc2SMC

* Reformat into a smc file an ascii file in which there are possible comments 
* preceded by "!" or "|", followed by one header line, followed by the data
* in two columns: time, data

* Sample control file:
* first line below starts the control, with line 1, col 1 = "N" in "Nominal":
*Control file for Asc2SMC      ! first line
*! Revision of program involving a change in the control file on this date:
*   11/23/10
*Name of summary file:
* asc2smc.sum
*!n2skip (-1=headers preceded by !; 0=no headers; otherwise number of headers to skip)
* -1
*!write headers to smc file (even if n2skip > 0)? (Y/N)
* Y
*!sps (0.0 = obtain from input file)
* 0
*!N columns to read, column number for time and data columns 
*!  (for files made using blpadflt, period is in column 1 and sd, pv, pa, rv, 
*!  aa are in columns 2, 3, 4, 5, 6, respectively)
*! Note: if sps .ne. 0.0, then column number for time is ignored (but a placeholder is
*! still needed--e.g., 1 1 1 (read one column, which contains the data; 1 20 1 would be the same)
*! But note: if the data are not in the first column, but only the data column is to be read
*! (because sps will be used to establish the time values),
*! then ncolumns must be the column corresponding to the data.  For example, assume that
*! the data are in column 3 and that columns 1 and 2 contain time and some other variable, but
*! the time column is not to be used (perhaps because accumulated error in creating the column
*! leads to a slight shift in the time values).  Then the input line should be:
*!  3 1 3
* 6 1 3
*Xfactr
* 1.0
*!Read input format (used if the format is such that the values are not separated by spaces,
*!in which case a free format cannot be used for input)?
*  N
*!If yes, specify a format; if not, still need a placeholder
* (3e13.5)
*!For output, use old (standard) smc format or new
*!higher precision format.   Specify "high" for
*!high precision; any other word defaults to standard
*!precision (but some word is needed as a placholder, even if
*!standard precision is desired).
* high
*!String to append to input file name for the output filename.
* .smc8
*Input file name ("stop" in any column to quit):
* \smsim\ts0.col
* \smsim\ts1.col
*STOP

* Dates: 05/01/02 - Written by D. Boore, patterned after col2smc
*        05/16/02 - Add header comments to smc file
*        11/01/03 - Prepare for LF95 (account for longer file names) by using
*                   trim_c.
*                   Add line to ctl indicating if there are any header lines
*        11/03/03 - Add variable "xfactr" to allow multiplicative factor 
*                   (can be used to change the units or the sign of the 
*                   time series values)
*        11/04/03 - Corrected bug (skipped first data value if file has headers)
*        11/13/03 - Allow for long input file names with several periods
*                   by looking for last period and replacing the portion
*                   from there to the end with ".smc"
*        12/01/05 - Closed nu_asc before return for another file.  Without this,
*                   the number of available unit numbers is exceeded after
*                   89 files.
*        12/08/06 - Construct output file name by appending "smc" to end of the
*                   input file name.
*        10/15/07 - Allow for a higher precision file format. If higher precision,
*                   use the extension "smc8", where "8" is the number of digits in
*                   the mantissa.
*        04/30/08 - Specify which of several column are to be used for ghe time and data
*                   axes.  Also, specify string to add to create output filename (no longer
*                   use "smc8" if high precision).
*        05/03/08 - Allow for sps being specified
*        08/17/08 - Replaced header_query in control file with n2skip (to allow
*                   skipping a specified number of headers if they are not preceded by !).
*        12/19/08 - Increased size of filename character strings
!        11/23/10 - Added revision date to control file
!                 - Corrected subtle bug (included "nu_asc" in call to skipcmnt after closing nu_asc)
!        02/02/12 - Add option for writing headers to output file (even if n2skip > 0)
!        05/18/12 - Fix bug found by Tim Huff, in which the sps is not calculated correctly for a 
!                   list if files in which sps varies from file to file.


      real time(200000), ts(200000), dum(10)
      character f_ctl*120, f_sum*120, f_asc*120, f_smc*120, f_asc_4*4
      character char_head(11)*80, comments(100)*80, hdr_cmmnts(100)*79
      character precision*20, cmnts2skip(50)*80
      character extension_c*30, read_format_c*10, format_c*80
      character date_ctl_correct*8, date_ctl_in*30
      character write_headers*10
      integer int_head(48), status
      real real_head(50)
      real sps_in, sps

      logical f_exist

      do i = 2, 11
        char_head(i) = '*'
      end do

      ndimen_max = 200000
      
      int_null = -32768
      do i = 1, 48
        int_head(i) = int_null
      end do
*      int_head(13) = 0
*      int_head(14) = 0

      real_null = 1.7e+38
      do i = 1, 50
        real_head(i) = real_null
      end do

      f_exist = .false.
      do while (.not. f_exist)
        f_ctl = ' '
        write(*, '(a)') 
     :    ' Enter name of control file (cr = asc2smc.ctl): '
        read(*, '(a)') f_ctl
        if (f_ctl(1:4) .eq. '    ') f_ctl = 'asc2smc.ctl'
        call trim_c(f_ctl, nc_f_ctl)
        inquire(file=f_ctl(1:nc_f_ctl), exist=f_exist)
        if (.not. f_exist) then
          write(*,'(a)') ' ******* FILE DOES NOT EXIST ******* '
        end if
      end do
      call get_lun(nu_ctl)
      open(unit=nu_ctl, file=f_ctl(1:nc_f_ctl), status='unknown')

      call skipcmnt(nu_ctl, cmnts2skip,nc_cmnts2skip)
      date_ctl_in = ' '
      read(nu_ctl,'(a)') date_ctl_in
      call trim_c(date_ctl_in,nc_date_ctl_in)
      
      date_ctl_correct = ' '
      date_ctl_correct = '02/02/12'
      call trim_c(date_ctl_correct,nc_date_ctl_correct)

      if (date_ctl_correct(1:nc_date_ctl_correct) /= 
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
      call trim_c(f_sum, nc_f_sum)
      call get_lun(nu_sum)
      open(unit=nu_sum, file=f_sum(1:nc_f_sum), status='unknown')
      write(nu_sum,'(a)') '  Summary program for program ASC2SMC'      

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) n2skip
 
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      write_headers = ' '
      read(nu_ctl,'(a)') write_headers
      call trim_c(write_headers, nc_write_headers)
      call upstr(write_headers)
 
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) sps_in
      write(*,'(a,1x,f7.2)') 
     :  ' sps_in = ', sps_in

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) ncol2read, icolx, icoly
      write(*,'(a,3(1x,i3))') 
     :  ' ncol2read, icolx, icoly = ', ncol2read, icolx, icoly

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read(nu_ctl,*) xfactr

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      read_format_c = ' '
      read(nu_ctl,'(a)') read_format_c
      call trim_c(read_format_c, nc_read_format_c)
      call upstr(read_format_c)
 
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      format_c = ' '
      read(nu_ctl,'(a)') format_c
      call trim_c(format_c, nc_format_c)
  
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      precision = ' '
      read(nu_ctl,'(a)') precision
      call trim_c(precision, nc_precision)
      call upstr(precision(1:nc_precision))

      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)
      extension_c = ' '
      read(nu_ctl,'(a)') extension_c
      call trim_c(extension_c, nc_extension_c)
 
      call skipcmnt(nu_ctl,cmnts2skip,nc_cmnts2skip)

      loop over files: DO

        f_asc = ' '
        read(nu_ctl,'(a)',iostat=status) f_asc
        IF (status /= 0) THEN
          EXIT loop over files
        END IF
        call trim_c(f_asc,nc_f_asc)
        f_asc_4 = ' '
        f_asc_4 = f_asc(1:4)
        call upstr(f_asc_4)
        if(f_asc_4 == 'STOP') EXIT loop over files
      
        f_exist = .false.
        inquire(file=f_asc(1:nc_f_asc), exist=f_exist)
        if (.not. f_exist) then
          write(*,'(a)') ' ******* FILE '//f_asc(1:nc_f_asc)//
     :        ' DOES NOT EXIST *******, going to next file' 
          write(nu_sum,'(a)') ' ******* FILE '//f_asc(1:nc_f_asc)//
     :        ' DOES NOT EXIST *******, going to next file' 
          CYCLE loop over files
        end if


        write(nu_sum,'(2a)') '  Processing file: ', f_asc(1:nc_f_asc)

        print *,' '
        print *,' Processing file: '//f_asc(1:nc_f_asc)

        call get_lun(nu_asc)
        open(unit=nu_asc, file=f_asc(1:nc_f_asc), status='unknown')
      
        write(nu_sum,'(1x,a,1x,i4)') ' nu_asc = ', nu_asc

        hdr_cmmnts = ' '
        if (n2skip == -1) then
          call skipcmnt(nu_asc, hdr_cmmnts, n_hdr_cmmnts)
        else if (n2skip > 0) then 
          do i = 1, n2skip
            read(nu_asc, '(a)') hdr_cmmnts(i)
          end do
          n_hdr_cmmnts = n2skip
        end if

        n_ts = 0
        loop over time points: DO
          if (n_ts == ndimen_max) then
            print *,
     :  ' n_ts exceeds ndimen_max; a shortened version will be saved'
            write(nu_sum,'(a)')
     :  ' n_ts exceeds ndimen_max; a shortened version will be saved'
            EXIT loop over time points
          end if
      
          if (read_format_c(1:nc_read_format_c) == 'Y') then
            read(nu_asc,format_c(1:nc_format_c), iostat=status) 
     :                   (dum(i), i=1, ncol2read)
            IF (status /= 0) THEN
              EXIT loop over time points
            END IF     
          else
            read(nu_asc,*,  iostat=status) (dum(i), i=1, ncol2read)
            IF (status /= 0) THEN
              EXIT loop over time points
            END IF
          end if
      
          n_ts = n_ts + 1
       
          if (sps_in == 0.0) then
            time(n_ts) = dum(icolx)
          else
            time(n_ts) = real(n_ts-1)/sps_in
          end if 
      
          ts(n_ts) = dum(icoly)

        END DO loop over time points
 
      close(nu_asc)
      
* Scale the time series:
      do i = 1, n_ts
        ts(i) = xfactr * ts(i)
      end do

* Get dt:

      dt = time(2) - time(1)    
      sps = 1.0/dt

      write(nu_sum,'(a,1x,i6,1x,f7.2)') '    n_ts, sps = ', n_ts, sps

      int_head(17) = n_ts
      
      real_head(2) = sps

      ncomments = 0
      ncomments = ncomments + 1
      comments(ncomments) = '|'
      if (n2skip /= 0 .and. write_headers(1:1) == 'Y') then  ! The input file has header comments
        do i = 1, n_hdr_cmmnts
          ncomments = ncomments + 1
          comments(ncomments) = '|'//hdr_cmmnts(i)
        end do
        ncomments = ncomments + 1
      end if
      comments(ncomments) = '|'
      ncomments = ncomments + 1
      comments(ncomments) = '| File '//f_asc(1:nc_f_asc)//
     : ' reformatted by program ASC2SMC'
      ncomments = ncomments + 1
      comments(ncomments) = '| Data values scaled by the factor '
      call trim_c(comments(ncomments), nc_c)
      write(comments(ncomments)(nc_c+1:nc_c+11),'(0p,1x,e10.3)')
     :  xfactr
      ncomments = ncomments + 1
      comments(ncomments) = '|'

      int_head(16) = ncomments

      f_smc = ' '
      
!* Find last period:
!      do i = nc_f_asc, 1, -1
!        if (f_asc(i:i) .eq. '.') then
!          indexp = i
!          go to 456
!        end if
!      end do
!      indexp = 1
!456   continue
!      f_smc = f_asc(1:indexp)//'smc'

!      if (precision(1:4) .eq. 'HIGH') then
!        f_smc = f_asc(1:nc_f_asc)//'.smc8'
!      else
!        f_smc = f_asc(1:nc_f_asc)//'.smc'
!      end if

      f_smc = f_asc(1:nc_f_asc)//extension_c(1:nc_extension_c)
      
      call trim_c(f_smc, nc_f_smc)

      write(*,'(a)') ' Writing file '//f_smc(1:nc_f_smc)
      write(nu_sum,'(a)') ' Writing file '//f_smc(1:nc_f_smc)

      call smcwrite(f_smc, precision,
     :              ts, char_head, 0, 
     :              int_head, real_head, comments)
     
      END DO loop over files

      close(nu_ctl)
      close(nu_sum)

      stop
      end

      include 'skip.for'
      include 'smcwrite.for'
      include 'get_lun.for'
      include 'upstr.for'
      include 'trim_c.for'
      include 'mnmaxidx.for'
      include 'skipcmnt.for'
      include 'distaz.for'
