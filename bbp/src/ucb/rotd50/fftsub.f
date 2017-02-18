      subroutine rdc (x,npts,iflag)

c     This subroutine removes a dc shift from the data

c     IFLAG = 0    remove the mean
c     IFLAG = 1    remove the mean value of the first 10 points
c     IFLAG = 2    manual set of DC value to be removed

      real x(1),sum,mean
      integer npts,iflag

      if (iflag .eq. 0) then
        sum = 0.0
        do 10 i=1,npts
          sum = x(i) + sum
  10    continue
        mean = sum/float(npts)

      elseif (iflag .eq. 1) then
        sum = 0.0
        do 20 i=1,10
          sum = x(i) + sum
  20    continue
        mean = sum / float(npts)

      else
        write (*,1000)
        read (*,1001) mean

      endif

      do 100 i=1,npts
        x(i) = x(i) - mean
  100 continue

c       write (*,1010) mean
c 1010 format( 2x,'Remove DC of ',f12.8)

      return
 1000 format( 2x,'Enter mean to be removed')
 1001 format( f12.8)
      end

c ---------------------------------------------------------------------------

      subroutine CosTaper (x,npts,tb,te)

c     This subroutine tapers the x array
    
      real x(1), arg
      integer npts,tb,te
      pi = 3.1415926

      if (tb .ne. 0.) then
        n = (npts*tb)/100
        do 10 i=1,n
          arg = pi*float(i-1)/float(n) + pi
          x(i) = x(i)*(1.+cos(arg))/2.
  10    continue
 1000   format( 2x,'Taper beginning ',i2,' percent')

      endif

      if (te .ne. 0.) then
        n = (npts*te)/100
        do 20 i=1,n
          arg = pi*float(i-1)/float(n) + pi
          x(npts-i+1) = x(npts-i+1) * (1.+cos(arg))/2.
  20    continue
 1001   format( 2x,'Taper end ',i2,' percent')

      endif
      return
      end

c ---------------------------------------------------------------------------

      subroutine cool ( signi, n, cx )
c     FFT subroutine.
c     signi = -1.  forward transform
c           =  1.  inverse transform
c     n = log base 2 (npts)

      complex cx(1), carg, temp, cw
      pi = 4. * atan(1.) * signi
      lx = 2**n
      j = 1
      do 30 i=1,lx
        if (i .gt. j) goto 10
        temp = cx(j)
        cx(j) = cx(i)
        cx(i) = temp
  10    m = lx/2
  20    if (j .le. m) goto 25
        j = j-m
        m = m/2
        if (m .ge. 1) goto 20
  25    j = j+m
  30  continue
      l = 1
  40  istep = l+l
      do 50 m=1,l
        carg = cmplx( 0., pi * float(m-1) / float(l) )
        cw = cexp(carg)
        do 45 i=m,lx,istep
          temp = cw * cx(i+l) 
          cx(i+l) = cx(i) - temp
          cx(i) = cx(i) + temp
  45    continue
  50  continue
      l = istep
      if (l .lt. lx) goto 40
      return
      end
