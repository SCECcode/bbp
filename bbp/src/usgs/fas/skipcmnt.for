! ------------------------------------------------------------------ skipcmnt
      subroutine skipcmnt(nu, comment, ncomments)

! Skip text comments in the file attached to unit nu, but save skipped 
! comments in character array comment.  Skip at least one line, and more as 
! long as the lines are preceded by "|" or "!".

! Dates: 04/16/01 - Written by D. Boore
!        12/07/01 - Added check for eof
!        11/04/03 - Use trim_c to trim possible leading blank
!        02/03/07 - Initialize comments to blank
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

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
! ------------------------------------------------------------------ skipcmnt
