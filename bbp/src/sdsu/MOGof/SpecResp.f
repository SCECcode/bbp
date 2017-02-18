C -------------------------------------------------------------------------
C --------------------------RESPONSESPEC-----------------------------------
C Subroutine to calculate the response spectrum
C -------------------------------------------------------------------------
      subroutine RESPONSESPEC(acc, dt, nt, period, spec, specDur, 
     +                        specDis)
      implicit none
      integer                 :: tag, nt, ii, k, i, j
      real                    :: beta, omega, dt, period, pga, spec
      real                    :: b2, o2, o3, pi, specDur, DL, DH
      real                    :: specDurA, DLA, DHA
      real                    :: specDis, specVel
      real, dimension(nt)     :: acc, x1, x2, a1, a2, zz, v2, v2out
      real, dimension(nt)     :: ad, adout
      real, dimension(2,2)    :: A, B
      data pi/3.14159265358979323846/

      omega = 2.0*pi/period
      beta = 0.05
      b2 = beta**2.0
      o2 = omega**2.0
      o3 = omega**3.0

      A(1,1) = exp(-beta*omega*dt)*((beta/sqrt(1.0-b2))*
     +sin(omega*sqrt(1.0-b2)*dt)+ cos(omega*sqrt(1.0-b2)*dt))

      A(1,2) = (exp(-beta*omega*dt)/(omega*sqrt(1.0-b2)))*
     +sin(omega*sqrt(1.0-b2)*dt)

      A(2,1) = A(1,2)*(-1.0)*(o2)

      A(2,2) = exp(-beta*omega*dt)*(-(beta/sqrt(1.0-b2))*
     +sin(omega*sqrt(1.0-b2)*dt)+ cos(omega*sqrt(1.0-b2)*dt))

      B(1,1) = exp(-beta*omega*dt)*((((2.0*b2)-1.0)/(o2*dt) + 
     +beta/omega)*sin(omega*sqrt(1.0-b2)*dt)/(omega*sqrt(1.0-b2)) + 
     +(((2.0*beta))/(o3*dt) + 1.0/o2)*cos(omega*sqrt(1.0-b2)*dt)) -
     +2.0*beta/((o3)*dt)

      B(1,2) = -exp(-beta*omega*dt)*((((2.0*b2)-1.0)/(o2*dt))*
     +sin(omega*sqrt(1.0-b2)*dt)/(omega*sqrt(1.0-b2)) + 
     +(((2.0*beta))/(o3*dt))*cos(omega*sqrt(1.0-b2)*dt))-
     +(1.0/o2)+2.0*beta/((o3)*dt)

      B(2,1) = exp(-beta*omega*dt)*(((((2.0*b2)-1.0)/(o2*dt) + 
     +beta/omega)*(cos(omega*sqrt(1.0-b2)*dt) - (beta/sqrt(1.0-b2))*
     +sin(omega*sqrt(1.0-b2)*dt)))- ((2.0*beta)/(o3*dt) + 1.0/o2)*
     +((omega*sqrt(1.0-b2))*sin(omega*sqrt(1.0-b2)*dt) + (beta*omega)*
     +cos(omega*sqrt(1.0-b2)*dt))) + 1.0/(dt*o2)

      B(2,2) = -exp(-beta*omega*dt)*((((2.0*b2)-1.0)/(o2*dt))*
     +(cos(omega*sqrt(1.0-b2)*dt) - 
     +(beta/sqrt(1.0-b2))*sin(omega*sqrt(1.0-b2)*dt))- 
     +((2.0*beta)/(o3*dt))*((omega*sqrt(1.0-b2))*
     +sin(omega*sqrt(1.0-b2)*dt)+
     +(beta*omega)*cos(omega*sqrt(1.0-b2)*dt))) - 1.0/(dt*o2)

      x1(1)=0.0
      x2(1)=0.0

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
          zz(k) = abs(-(2.0*beta*omega*x2(k)+o2*x1(k)))
          v2(k) = x2(k)**2.0
      enddo
      spec = maxval(zz)
      specDis = maxval(abs(x1))
      specVel = maxval(abs(x2))
C--------------------------------------------------------------------------
C Duration Calculation and Comparison
      call DINTEG(v2, dt, nt, v2out)

      do j = 1,nt
           v2out(j) = v2out(j)/v2out(nt)
           if(v2out(j).le.0.05) DL = REAL(j-1)*dt
           if(v2out(j).le.0.75) DH = REAL(j-1)*dt
      enddo

      specDur = REAL(DH-DL)/period

C--------------------------------------------------------------------------
C      To calculate the Arias Duration
C      do k = 1,nt
C          ad(k) = (-(2*beta*omega*x2(k)+o2*x1(k)))**2
C      enddo
C      call DINTEG(ad, dt, nt, adout)
C      do j = 1,nt
C           adout(j) = adout(j)/adout(nt)
C           if(adout(j).le.0.05) DLA = REAL(j-1)*dt
C           if(adout(j).le.0.95) DHA = REAL(j-1)*dt
C      enddo
C      specDurA = REAL(DHA-DLA)/period
C--------------------------------------------------------------------------
      end subroutine RESPONSESPEC
C -------------------------------------------------------------------------
C -------------------------------------------------------------------------
C --------------------------RESPSPECNGA-----------------------------------
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

      omega = 2.0*pi/period
      beta = 0.05
      b2 = beta**2.0
      o2 = omega**2.0
      o3 = omega**3.0

      A(1,1) = exp(-beta*omega*dt)*((beta/sqrt(1.0-b2))*
     +sin(omega*sqrt(1.0-b2)*dt)+ cos(omega*sqrt(1.0-b2)*dt))

      A(1,2) = (exp(-beta*omega*dt)/(omega*sqrt(1.0-b2)))*
     +sin(omega*sqrt(1.0-b2)*dt)

      A(2,1) = A(1,2)*(-1.0)*(o2)

      A(2,2) = exp(-beta*omega*dt)*(-(beta/sqrt(1.0-b2))*
     +sin(omega*sqrt(1.0-b2)*dt)+ cos(omega*sqrt(1.0-b2)*dt))

      B(1,1) = exp(-beta*omega*dt)*((((2.0*b2)-1.0)/(o2*dt) + 
     +beta/omega)*sin(omega*sqrt(1.0-b2)*dt)/(omega*sqrt(1.0-b2)) + 
     +(((2.0*beta))/(o3*dt) + 1.0/o2)*cos(omega*sqrt(1.0-b2)*dt)) -
     +2.0*beta/((o3)*dt)

      B(1,2) = -exp(-beta*omega*dt)*((((2.0*b2)-1.0)/(o2*dt))*
     +sin(omega*sqrt(1.0-b2)*dt)/(omega*sqrt(1.0-b2)) + 
     +(((2.0*beta))/(o3*dt))*cos(omega*sqrt(1.0-b2)*dt))-
     +(1.0/o2)+2.0*beta/((o3)*dt)

      B(2,1) = exp(-beta*omega*dt)*(((((2.0*b2)-1.0)/(o2*dt) + 
     +beta/omega)*(cos(omega*sqrt(1.0-b2)*dt) - (beta/sqrt(1.0-b2))*
     +sin(omega*sqrt(1.0-b2)*dt)))- ((2.0*beta)/(o3*dt) + 1.0/o2)*
     +((omega*sqrt(1.0-b2))*sin(omega*sqrt(1.0-b2)*dt) + (beta*omega)*
     +cos(omega*sqrt(1.0-b2)*dt))) + 1.0/(dt*o2)

      B(2,2) = -exp(-beta*omega*dt)*((((2.0*b2)-1.0)/(o2*dt))*
     +(cos(omega*sqrt(1.0-b2)*dt) - 
     +(beta/sqrt(1.0-b2))*sin(omega*sqrt(1.0-b2)*dt))- 
     +((2.0*beta)/(o3*dt))*((omega*sqrt(1.0-b2))*
     +sin(omega*sqrt(1.0-b2)*dt)+
     +(beta*omega)*cos(omega*sqrt(1.0-b2)*dt))) - 1.0/(dt*o2)

      x1(1)=0.0
      x2(1)=0.0

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
          zz(k) = abs(-(2.0*beta*omega*x2(k)+o2*x1(k)))
      enddo
C--------------------------------------------------------------------------
      end subroutine RESPSPECNGA
C -------------------------------------------------------------------------

