c ----------------------------------------------------------------------
c     This subroutine calculates the response spectrum. 
      subroutine CalcRsp ( acc, npts, dt, w, damping, 
     1              nQ, time, SA, polarity, minTime, maxTime)

      real w(1), damping(1), acc(1), dt, SA(1), time(1), minTime(1), 
     1           maxTime(1)
      integer npts, nQ, i, j , polarity(1), timeIndex

c     LOOP FOR EACH FREQUENCY
      do i=1,nQ
         if ( w(i) .gt. 628. ) then
            call CalcPGA ( acc, npts, sa(i), timeIndex, polarity(i) )
         else
            call CalcOneRsp ( acc, npts, dt, w(i), damping(i), sa(i),
     1                 timeIndex, polarity(i), minTime(i),maxTime(i))
         endif
         time(i) = (timeIndex-1) * dt
      enddo
	  
      return
      end
	  
c ----------------------------------------------------------------------

      subroutine CalcPGA ( acc, npts, sa, timeIndex, polarity  )

      real acc(1), sa
      integer npts, polarity, timeIndex
	  
      accMax = 0.0
      do i=1,npts
         if ( abs(acc(i)) .gt. accMax ) then
            accMax = abs(acc(i) )
            sa = accMax
            if ( acc(i) .ge. 0. ) then
               polarity = 1.
            else 
               polarity = -1
            endif
               timeIndex = i
         endif
      enddo
      write (*,'( 2x,''pga ='',f10.4)') sa

      return
      end
	    
c ----------------------------------------------------------------------

      subroutine CalcRspTH ( acc, npts, dt, w, damping, rsp )
      include 'baseline.h'

      real rsp(1), w, damping, acc(1), dt
      integer npts
	  
c     Compute coeff 
      call coeff ( w, damping, dt )

c     CALCULATE THE RESPONSES
      call brs ( acc, w, damping, npts, rsp )

      return
      end

c ----------------------------------------------------------------------

      subroutine CalcOneRsp ( acc, npts, dt, w, damping, 
     1                 sa, timeIndex, polarity, minTime, maxTime)
      include 'baseline.h'

      real rsp(MAXPTS), w, damping, acc(1), dt, SA, minTime, maxTime
      integer npts, i, j, polarity, timeIndex, iTime, iTime2
	  
c     Compute coeff 
      call coeff ( w, damping, dt )

c     CALCULATE THE RESPONSES
      call brs ( acc, w, damping, npts, rsp )

c     FIND THE MAXIMUM OF THE RESPONSE
      SA = 0.0
      iTime = int( minTime/dt) + 1
      iTime2 = int( maxTime/dt) -1
      do j=iTime,iTime2
         if (abs(rsp(j)) .gt. SA) then
            SA = abs(rsp(j))
            timeIndex = j
            if ( rsp(j) .ge. 0. ) then
               polarity = 1
            else
               polarity = -1
            endif
         endif
      enddo

      return
      end

c ----------------------------------------------------------------------

      subroutine coeff ( w, beta1, dt1 )
      real beta1, dt1, w
      real*8 a11, a12, a21, a22, b11, b12, b21, b22
      real*8 beta, dt, t1, t2, t3, t4, s1, s2
      common /coef/a11,a12,a21,a22,b11,b12,b21,b22

      beta = dble( beta1 )
      dt = dble( dt1 )
c
c     Set up repeated terms
      t1 = sqrt(1.-beta**2)
      t2 = sin (w*t1*dt)
      t3 = cos (w*t1*dt)
      t4 = exp (-beta*w*dt)
      s1 = (2.*beta**2-1.) / (w**2*dt)
      s2 = 2.*beta / (w**3*dt)
c	  write (*,'( 10f10.5)') t1, t2, t3, t4, s1, s2
c
c     calculate the a's
      a11 = t4*(beta*t2/t1+t3)
      a12 = t4*t2 / (w*t1)
      a21 = -t4*w*t2 / t1
      a22 = t4*(t3-beta*t2/t1)
c
c     calculate the b's
      b11 = t4*((s1+beta/w)*t2 / (w*t1) + (s2+1./w**2)*t3) - s2
      b12 = -t4*(s1*t2/(w*t1)+s2*t3) - 1./w**2 + s2
      b21 = (s1+beta/w) * (t3-beta*t2/t1)
      b21 = t4*(b21 - (s2+1./w**2)*(w*t1*t2+beta*w*t3)) + 1./(w**2*dt)
      b22 = s1*(t3-beta*t2/t1)
      b22 = -t4*(b22 - s2*(w*t1*t2+beta*w*t3)) - 1./(w**2*dt)

      return
      end

c ----------------------------------------------------------------------
      subroutine brs (x,w,beta,npts,rsp)
      real x(1), rsp(1), beta
      real w
      real*8 d, v, a, z, ap1, dp1, vp1, t1, t2
      real*8 a11, a12, a21, a22, b11, b12, b21, b22
      common /coef/ a11,a12,a21,a22,b11,b12,b21,b22
c
c     initialize
      t1 = 2.*beta*w
      t2 = w**2
      d = 0.
      v = 0.
      a = 0.
c
c     calculate the response
      do 10 i=1,npts
        ap1 = dble( x(i) )
        dp1 = a11*d + a12*v + b11*a + b12*ap1
        vp1 = a21*d + a22*v + b21*a + b22*ap1
        z = -(t1*vp1 + t2*dp1)

c        absolute acc
c        rsp(i) = sngl( z )

c       Psedo-acceleration
        rsp(i) = sngl( dp1 ) * t2
        a = ap1
        v = vp1
        d = dp1
  10  continue

      return
      end 

c ----------------------------------------------------------------------
