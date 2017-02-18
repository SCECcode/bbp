! --------------------------- BEGIN TRIM_C -----------------------
      subroutine trim_c(cstr, nchar)

! strips leading and trailing blanks from cstr, returning the
! result in cstr, which is now nchar characters in length

! Strip off tabs also.

! Here is a sample use in constructing a column header, filled out with 
! periods:

!* Read idtag:
!        idtag = ' '
!        read(nu_in, '(1x,a)') idtag
!        call trim_c(idtag, nc_id)
!* Set up the column headings:
!        colhead = ' '
!        colhead = idtag(1:nc_id)//'......' ! nc_id + 6 > length of colhead

! Dates: 12/23/97 - written by D. Boore
!        12/08/00 - pad with trailing blanks.  Otherwise some comparisons
!                   of the trimmed character string with other strings
!                   can be in error because the original characters are left
!                   behind after shifting.  For example, here is a string
!                   before and after shifting, using the old version:
!                      col:12345
!                           MTWH  before
!                          MTWHH  after (but nc = 4).
!        03/21/01 - Check for a zero length input string
!        11/09/01 - Change check for zero length string to check for all blanks
!        10/19/09 - Strip off tabs
!        04/28/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      character cstr*(*)

      if(cstr .eq. ' ') then
        nchar = 0
        return
      end if

      nend = len(cstr)

! Replace tabs with blanks:

      do i = 1, nend
        if(ichar(cstr(i:i)) .eq. 9) then
           cstr(i:i) = ' '
        end if
      end do



!      if(nend .eq. 0) then
!        nchar = 0
!        return
!      end if

      do i = nend, 1, -1
        if (cstr(i:i) .ne. ' ') then
           nchar2 = i
           goto 10
        end if
      end do

10    continue

      do j = 1, nchar2
        if (cstr(j:j) .ne. ' ') then
          nchar1 = j
          goto 20
        end if
      end do

20    continue
   
      nchar = nchar2 - nchar1 + 1
      cstr(1:nchar) = cstr(nchar1: nchar2)
      if (nchar .lt. nend) then
        do i = nchar+1, nend
          cstr(i:i) = ' '
        end do
      end if

      return
      end
! --------------------------- END TRIM_C -----------------------

