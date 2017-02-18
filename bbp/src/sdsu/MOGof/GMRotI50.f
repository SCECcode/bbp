C -------------------------------------------------------------------------
C John Mayhew
C Thesis 2009
C -------------------------------NGAROT------------------------------------
C Input a seismogram and get out the 17 values coresponding to the 50th 
C percentile
C -------------------------------------------------------------------------
      subroutine NGAROT(AcX, AcY, nt, dt, Per50)
      implicit none
      integer               :: nt, i, j, k, nuM
      real                  :: deg, pi, dt, tmpD, tmpDP, XR, YR
      real, dimension(nt)   :: XROT, YROT, AX, AY, AcX, AcY
      real, dimension(17)   :: Per50
      real, dimension(16)   :: period
      real, dimension(91,17):: values, ValOrg
      integer, dimension(91):: OrgList, penF

      data period/.1,.15,.2,.25,.3,.4,.5,.75,
     +1,1.5,2,3,4,5,7.5,10/

      data pi/3.14159265358979323846/
      do k = 1,16
           call RESPSPECNGA(AcX, dt, nt, period(k), AX)
           call RESPSPECNGA(AcY, dt, nt, period(k), AY)
           do i = 1,91
                deg = real(i-1)*pi/180.0

                do j = 1,nt
                     XROT(j) = COS(deg)*AX(j)+SIN(deg)*AY(j)
                     YROT(j) = COS(deg)*AY(j)-SIN(deg)*AX(j)
                enddo
                values(i,k) = SQRT(MAXVAL(ABS(XROT))*MAXVAL(ABS(YROT)))
           enddo
      enddo
      do i = 1,91
           deg = real(i-1)*pi/180.0

           do j = 1,nt
                XROT(j) = COS(deg)*AcX(j)+SIN(deg)*AcY(j)
                YROT(j) = COS(deg)*AcY(j)-SIN(deg)*AcX(j)
           enddo
           values(i,17) = SQRT(MAXVAL(ABS(XROT))*MAXVAL(ABS(YROT)))
      enddo

      do i = 1,17
           call organize(values(1:91,i), OrgList, 91)
           do j = 1,91
                ValOrg(OrgList(j),i) = values(j,i)
           enddo

           per50(i) = (ValOrg(46,i))
      enddo

C Compute the penalty function
      tmpDP = 10000.0
      nuM = 0

      do j = 1,91
           tmpD = 0.0
           do k = 1,17
                tmpD = tmpD + ((values(j,k)/per50(k))-1)**2
           enddo
           penF(j) = tmpD/17.0
           if(penF(j) .LT. tmpDP) then
              tmpDP = penF(j)
              nuM = j
           endif
      enddo

      do i = 1,17
         per50(i) = values(nuM,i)
      enddo

      end subroutine NGAROT
C -------------------------------------------------------------------------
C -------------------------------------------------------------------------
C ------------------------------organize-----------------------------------
C Compares the two spectral response series
C -------------------------------------------------------------------------
      subroutine organize(INPUT, OUTPUT2, num)
      implicit none
      integer                :: i, j, k, num
      real                   :: pi
      real, dimension(num)   :: INPUT, OUTPUT
      integer, dimension(num):: OUTPUT2
      data pi/3.14159265358979323846/
      
      OUTPUT(1) = minval(INPUT)
      do i = 2,num-1
         OUTPUT(i) = maxval(INPUT)
         do j = 1,num
            if(INPUT(j) .GT. OUTPUT(i-1) .AND. 
     +         INPUT(j) .LT. OUTPUT(i)) then
               OUTPUT(i) = INPUT(j)
            endif
         enddo
      enddo
      OUTPUT(num) = maxval(INPUT)

      do i = 1,num
         do j = 1,num
            if(OUTPUT(i) .EQ. INPUT(j)) then
               OUTPUT2(j) = i
            endif
         enddo
      enddo

      end subroutine organize
C -------------------------------------------------------------------------
C -------------------------------------------------------------------------
C ---------------------------RESPSPECNGA-----------------------------------
C Subroutine to calculate the response spectrum
C -------------------------------------------------------------------------
      subroutine RESPSPECNGA(acc, dt, nt, period, zz)
      implicit none
      integer                 :: tag, nt, ii, k, i, j
      real                    :: beta, omega, dt, period, pga, spec
      real                    :: b2, o2, o3, pi, specDur, DL, DH
      real                    :: specDis, specVel
      real, dimension(nt)     :: acc, x1, x2, a1, a2, zz, v2, v2out
      real, dimension(2,2)    :: A, B
      data pi/3.14159265358979323846/

      omega = 2*pi/period
      beta = 0.05
      b2 = beta**2
      o2 = omega**2
      o3 = omega**3

      A(1,1) = exp(-beta*omega*dt)*((beta/sqrt(1-beta*2))*
     +sin(omega*sqrt(1-b2)*dt)+ cos(omega*sqrt(1-b2)*dt))

      A(1,2) = (exp(-beta*omega*dt)/(omega*sqrt(1-b2)))*
     +sin(omega*sqrt(1-b2)*dt)

      A(2,1) = A(1,2)*(-1)*(o2)

      A(2,2) = exp(-beta*omega*dt)*(-(beta/sqrt(1-beta*2))*
     +sin(omega*sqrt(1-b2)*dt)+ cos(omega*sqrt(1-b2)*dt))

      B(1,1) = exp(-beta*omega*dt)*((((2*b2)-1)/(o2*dt) + 
     +beta/omega)*sin(omega*sqrt(1-b2)*dt)/(omega*sqrt(1-b2)) + 
     +(((2*beta))/(o3*dt) + 1/o2)*cos(omega*sqrt(1-b2)*dt)) -
     +2*beta/((o3)*dt)

      B(1,2) = -exp(-beta*omega*dt)*((((2*b2)-1)/(o2*dt))*
     +sin(omega*sqrt(1-b2)*dt)/(omega*sqrt(1-b2)) + 
     +(((2*beta))/(o3*dt))*cos(omega*sqrt(1-b2)*dt))-
     +(1/o2)+2*beta/((o3)*dt)

      B(2,1) = exp(-beta*omega*dt)*(((((2*b2)-1)/(o2*dt) + 
     +beta/omega)*(cos(omega*sqrt(1-b2)*dt) - (beta/sqrt(1-b2))*
     +sin(omega*sqrt(1-b2)*dt)))- ((2*beta)/(o3*dt) + 1/o2)*
     +((omega*sqrt(1-b2))*sin(omega*sqrt(1-b2)*dt) + (beta*omega)*
     +cos(omega*sqrt(1-b2)*dt))) + 1/(dt*o2)

      B(2,2) = -exp(-beta*omega*dt)*((((2*b2)-1)/(o2*dt))*
     +(cos(omega*sqrt(1-b2)*dt) - 
     +(beta/sqrt(1-b2))*sin(omega*sqrt(1-b2)*dt))- 
     +((2*beta)/(o3*dt))*((omega*sqrt(1-b2))*sin(omega*sqrt(1-b2)*dt)+
     +(beta*omega)*cos(omega*sqrt(1-b2)*dt))) - 1/(dt*o2)

      x1(1)=0
      x2(1)=0

      a1(1)=acc(1)
      a2(1)=acc(2)
      do ii = 2,nt-1
          a1(ii)=acc(ii)
          a2(ii)=acc(ii+1)
      enddo

      do ii = 2,nt
          x1(ii) = A(1,1)*x1(ii-1)+A(1,2)*x2(ii-1)+
     +B(1,1)*a1(ii-1)+B(1,2)*a2(ii-1)
          x2(ii) = A(2,1)*x1(ii-1)+A(2,2)*x2(ii-1)+
     +B(2,1)*a1(ii-1)+B(2,2)*a2(ii-1)
      enddo
     
      do k = 1,nt
          zz(k) = abs(-(2*beta*omega*x2(k)+o2*x1(k)))
      enddo
C--------------------------------------------------------------------------
      end subroutine RESPSPECNGA
C -------------------------------------------------------------------------
C --------------------------RESIZE-----------------------------------------
C Subroutine to resize the seismograms (Spline)
C -------------------------------------------------------------------------
      subroutine ReSize(seism, output, dt, dt2, nt, nt2)
      implicit none
      integer              :: nt, nt2, i, j, k, loc1
      real                 :: a, b, c, d, dt, dt2, xn, xo
      real, dimension(nt)  :: seism, x1, dif1, dif2
      real, dimension(nt2) :: output, x2


      call DIFF(seism, dt, nt, dif1)
      call DIFF(dif1, dt, nt, dif2)

C Generate time vectors
      do i = 1,nt
          x1(i) = real(i-1)*dt
      enddo
      do i = 1,nt2
          x2(i) = real(i-1)*dt2
      enddo

C Evaluate a linear interploation between points
      output(1) = seism(1)
      output(nt2) = seism(nt)
    
      do j = 2,nt2-1
         loc1 = int(dt2*(j-1)/dt)+1
         d = seism(loc1)
         c = dif1(loc1)
         b = dif2(loc1)/2.0
         a = (seism(loc1+1) - b*(dt**2.0) - c*dt - d)/dt**3.0
         output(j) = a*(x2(j)-x1(loc1))**3.0+b*(x2(j)-x1(loc1))**2.0+
     +               c*(x2(j)-x1(loc1))+d
      enddo 

      end subroutine ReSize
C -------------------------------------------------------------------------