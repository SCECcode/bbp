
! ----------------------- BEGIN RCC -----------------------------------
      subroutine RCC (KPARM,CSTR_IN,NCSTR_IN,C,NC,ibgn,iend)
! RCC - [MUELLER.FS.GEN]CSMGENLB
! Read character string C(1:NC), the KPARMth element in CSTR_IN(1:NCSTR_IN).
!  CSTR_IN contains one-or-more mixed-type fields
!(OLD COMMENT:c  separated by blanks or commas to the right of an optional '='.)
!  KPARMth field contains C to the right of an optional '>'.
! If the KPARMth field is nonexistent or empty, RCC sets NC=0 and returns.
!  This simulates a <carriage-return> in a Q-format read.
!  (An empty field must be separated by commas, not blanks.)

! Dates:  xx/xx/xx - Written by C. Mueller, USGS
!         06/22/98 - Modified by L. Baker
!         11/12/00 - Dave Boore introduced a working array, and a call to 
!                    rmv_tabs (this replaces a tab with a blank, but does not
!                    expand the tabs).  This was done because a file in which 
!                    fields are separated by tabs characters will not be parsed 
!                    correctly; apparently the tab characters (ascii 9) is 
!                    not the same thing as a blank.
!         12/04/00 - added ibgn, iend to rcf, rci, rcc calling arguments
!         06/11/01 - Disable use of "=" to indicate where the parameters
!                    start.  I did this because I sometimes add comments
!                    to the right of the parameters, and if these comments
!                    contain the character "=" then the parameters are not
!                    parsed correctly.  I reasoned that eliminating the check
!                    for "=" is the best solution to the problem (as opposed
!                    to restricting the use of "=" to only indicate 
!                    where the parameters start) because I rarely used
!                    "=" in this way.
!         06/12/01 - Minor modification (i2 = 0 rather than i2 = 1)
!         10/29/02 - Changed "rmv_tabs" to "tabs_rmv"
!         07/25/10 - Increase dimension of cstr
!        05/01/15 - Replaced comment characters * or C with ! (The Fortran 95 standard)

      character CSTR_IN*(*),C*(*),term*1

      character cstr*500

      if (ncstr_in .gt. 500) then
         write(*,'(a, 1x,i5, a)') '  WARNING: IN RCC, NCSTR_IN = ', 
     :          ncstr_in,
     :         ' WHICH IS GREATER THAN THE DIMENSION OF THE'//
     :         ' WORK ARRAY; QUITTING!!!!!'
         return
      end if

      ncstr = ncstr_in

      cstr = ' '
      cstr(1: ncstr) = cstr_in
      call tabs_rmv(cstr, ncstr)

! Eliminate trailing white space.
      N = NCSTR+1
1     N = N-1
      if (CSTR(N:N).eq.' '.or.CSTR(N:N).eq.',') goto 1

! Changes on 6/11/01:
!c Find optional '='.
!      I2 = INDEX(CSTR(1:N),'=')
!      term = '='
      i2 = 0         ! DMB, 6/12/01
      term = ' '     ! DMB, 6/11/01
! Changes on 6/11/01:

! Find start and end of KPARMth field.
      do 4 K=1,KPARM
         I1 = I2
2        I2 = I1
3        I1 = I1+1
         if (CSTR(I1:I1).eq.' ') goto 3
         if (term.eq.' '.and.CSTR(I1:I1).eq.',') then
            term = ','
            goto 2
         end if
         IBLANK = INDEX(CSTR(I1:N),' ')
         ICOMMA = INDEX(CSTR(I1:N),',')
         if (IBLANK.eq.0.and.ICOMMA.eq.0.and.K.lt.KPARM) goto 801
         if (IBLANK.eq.0) IBLANK = N-I1+2
         if (ICOMMA.eq.0) ICOMMA = N-I1+2
         I2 = I1+MIN(IBLANK,ICOMMA)-1
         if (I2.gt.N) then
            term = '.'
         else
            term = CSTR(I2:I2)
         end if
4        continue
! Check for empty field.
      if (I2.eq.I1) goto 801
! Find optional '>'.
      I1 = I1+INDEX(CSTR(I1:I2),'>')
! Get C.
      C = CSTR(I1:I2-1)
      NC = I2-I1
      ibgn = i1
      iend = i2-1
      return
! Errors.
801   NC = 0
      return
      end
! ----------------------- END RCC -----------------------------------


      