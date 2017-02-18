C --------------------------FileRead----------------------------------
C Author:  Nicolas Luco
C Last Revised:  11 September 2003
C Reference:  "Dynamics of Structures" (1995) by A.K. Chopra 
C
C Alpha = 2% as used in:
C   "EFFECTS OF EARTHQUAKE RECORD SCALING ON NONLINEAR STRUCTURAL RESPONSE"
C   Report on PEER-LL Program Task 1G00 Addendum (Sub-Task 1 of 3)
C   Nicolas Luco & Paolo Bazzurro
C   AIR Worldwide Corporation (San Francisco, CA)
C
C INPUT:
C ======
C T     = periods                           (scalar )
C z     = damping ratios                    ( " )
C dy    = yield displacements               ( " )
C alpha = strain hardening ratios           ( " )
C ag    = ground acceleration time history  ( length(ag) x 1 )
C Dtg   = time step of ag                   ( scalar )
C Dt    = analyis time step                 ( scalar )
C
C OUTPUT:
C =======
C   d  = relative displacement response time history  ( length(ag) ) 
C   v  = relative velocity response time history      ( " )
C   a  = acceleration response time history           ( " )
C   fs = force response time history                  ( " )
C
C --------------------------------------------------------------------
      subroutine InElastic( T, dy, ag, dt, nt, maxD)
      implicit none
      integer             :: i, j, nt
      real                :: T, z, dy, alpha, dt, AA, BB, fy
      real                :: Sd, Sv, Sa
      real                :: ggg, beta, m, w, pi, c, k, kalpha
      real                :: ki, DPi, Ddi, maxD, maxV, maxA
      real                :: fsmax, fsmin, Df, DR, Dvi, Kki
      real, dimension(nt) :: tt, d, v, a, fs, ag, p
C --------------------------------------------------------------------
      data ggg/0.5/
      data beta/0.16666667/
      data alpha/0.02/
      data z/0.05/
      data pi/3.14159265358979323846/
C --------------------------------------------------------------------
C m*a + c*v + fs(k,fy,kalpha) = p
C --------------------------------------------------------------------
      m = 1.0
      w = 2.0*pi / T
      c = z * (2.0*m*w)
      k = (w**2.0 * m)
      fy = k * dy
      kalpha = k * alpha

      do i = 1,nt
         tt(i) = REAL(i-1)*dt
         p(i) = -ag(i)*m
         d(i) = 0.0
         v(i) = 0.0
         a(i) = 0.0
         fs(i) = 0.0
      enddo

      a(1) = ( p(1) - c*v(1) - fs(1) ) / m
      AA = (1.0/(beta*dt))*m + (ggg/beta)*c
      BB = (1.0/(2.0*beta))*m + dt*((ggg/(2.0*beta))-1.0)*c

      do i = 1,(nt-1)

          DPi = p(i+1)-p(i) + AA*v(i) + BB*a(i)

          ki = k
        if(((DPi .GT. 0.0) .AND. (fs(i) .GE. fy+kalpha*(d(i)-dy))).OR. 
     +   ((DPi .LT. 0.0) .AND. (fs(i) .LE. -fy+kalpha*(d(i)+dy))))then
              ki = kalpha
          endif


C = ki + GAMMA/(BETA*dt)*c + 1/(BETA*dt**2.0)*m
          Kki = ki + AA/dt  
    
          Ddi = DPi / Kki
          fs(i+1) = fs(i) + ki*Ddi
          d(i+1) = d(i) + Ddi

          fsmax =  fy + kalpha*(d(i+1)-dy)
          fsmin = -fy + kalpha*(d(i+1)+dy)

          if( fs(i+1) .GT. fsmax ) then
              fs(i+1) = fsmax
              Df = fs(i+1) - fs(i) + (Kki-ki)*Ddi
              DR = DPi - Df
              Ddi = Ddi + DR / Kki
              d(i+1) = d(i) + Ddi
          endif

          if( fs(i+1) .LT. fsmin ) then
              fs(i+1) = fsmin
              Df = fs(i+1) - fs(i) + (Kki-ki)*Ddi
              DR = DPi - Df
              Ddi = Ddi + DR / Kki
              d(i+1) = d(i) + Ddi
          endif

          Dvi = (ggg/(beta*dt))*Ddi - (ggg/beta)*v(i) + 
     +          dt*(1.0-(ggg/(2.0*beta)))*a(i)
          v(i+1) = v(i) + Dvi
          a(i+1) = (p(i+1) - c*v(i+1) - fs(i+1))/m
    
      enddo

C Spectral values
C ---------------
      Sd = maxval( abs(d) )
      Sv = Sd * (2.0*pi/T)
      Sa = Sd * (2.0*pi/T)**2.0

      maxD = maxval(abs(d))
      maxV = maxval(abs(v))
      maxA = maxval(abs(a))
 
      end subroutine InElastic
C --------------------------------------------------------------------
