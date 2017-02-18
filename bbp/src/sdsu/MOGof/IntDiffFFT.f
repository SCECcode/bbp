C --------------------------------------------------------------------
C --------------------------DINTEG------------------------------------
C Subroutine to perform a definite integral
C --------------------------------------------------------------------
      subroutine DINTEG(seism, dt2, nt2, output)
      implicit none
      integer                 :: nt2, i, j
      real                    :: dt2
      real, dimension(nt2)   :: seism, output
      
      output(1) = 0.0
      do i=2,nt2
         output(i) = output(i-1)+(seism(i-1)*dt2)+
     +              (seism(i)-seism(i-1))*(0.5)*dt2
      enddo

      end subroutine DINTEG
C --------------------------------------------------------------------
C --------------------------DIFF--------------------------------------
C Subroutine to differentiate (4th Order)
C central (-f(i+2)+8*f(i+1)-8*f(i-1)+f(i-2))/(12*hh);
C forward (-f(i+2)+4*f(i+1)-3*f(i))/(2*hh);
C backward (3*f(i)-4*f(i-1)+f(i-2))/(2*hh);
C --------------------------------------------------------------------
      subroutine DIFF(seism, dt2, nt2, output)
      implicit none
      integer                 :: nt2, i
      real                    :: dt2
      real, dimension(nt2)   :: seism, output

C Forward difference
      output(1) = (-seism(3)+4.0*seism(2)-3.0*seism(1))/(2.0*dt2)
      output(2) = (-seism(4)+4.0*seism(3)-3.0*seism(2))/(2.0*dt2)
C Central Difference
      do i=3,nt2-2
        output(i) = (-seism(i+2)+8.0*seism(i+1)-8.0*seism(i-1)+
     +               seism(i-2))/(12.0*dt2)
      enddo
C Backward Difference
      output(nt2-1) = (3.0*seism(nt2-1)-4.0*seism(nt2-2)+
     +                  seism(nt2-3))/(2.0*dt2)
      output(nt2)   = (3.0*seism(nt2)-4.0*seism(nt2-1)+
     +                  seism(nt2-2))/(2.0*dt2)

      end subroutine DIFF

      subroutine MYSIN(val, sum)
      implicit none
      real                  :: myval, val, sum, term, pi
      integer               :: n

      data pi/3.14159265358979323846/

      if (val .ge. 0) then
         myval=mod(val,2*pi)
      endif
      if (val .lt. 0) then
         myval=mod(val,-2*pi)
      endif
      sum=0.0
      term=myval
      do n=1,20
         sum=sum+term
         term=term*(-myval**2/(2*n*(2*n+1)))
      enddo

      end subroutine MYSIN

      subroutine MYCOS(val, sum)
      implicit none
      real                  :: myval, val, sum, term, pi
      integer               :: n

      data pi/3.14159265358979323846/

      myval=val+pi/2
      if (myval .ge. 0) then
         myval=mod(myval,2*pi)
      endif
      if (myval .lt. 0) then
         myval=mod(myval,-2*pi)
      endif
      myval=myval-pi/2

      sum=0.0
      term=1.0
      do n=1,20
         sum=sum+term
         term=term*(-myval**2/(2*n*(2*n-1)))
      enddo

      end subroutine MYCOS


C --------------------------------------------------------------------
C --------------------------FFT---------------------------------------
C Subroutine to FFT a vector
C exp(CMPLX(a,b)) = exp(a)*(cos(b)+i*sin(b))
C http://www.dspguide.com/ch8/6.htm
C---------------------------------------------------------------------
C New routine is modeled off of (and verified) Matlab routine
C --------------------------------------------------------------------
C This subroutine is based off of the psuedo-code in Numerical Recipes
C If sign = 1, routine performs DFT and output is real(odd), complex(even)
C If sign = -1, input is real(odd), complex(even) and output is all real
C --------------------------------------------------------------------
      subroutine FFT(DATA, n, isn)
      implicit none
      integer               :: n, nn, i, j, k, isn, M, Istep, Mmax
      real                  :: WR, WI, WPR, WPI, Wtemp, theta
      real                  :: tempR, tempI, pi, tempcos, tempsin
      real, dimension(n)    :: DATA, XXr, XXi, XXa

      data pi/3.14159265358979323846/

!      write(*,*)'nt =', n

      do k = 1,n
          XXr(k) = 0.0
          XXi(k) = 0.0
          do j = 1,n
             call MYSIN(-2.0*pi*(j-1)*(k-1)/n, tempsin)
             call MYCOS(-2.0*pi*(j-1)*(k-1)/n, tempcos)
!              XXr(k)=XXr(k)+DATA(j)*COS(-2.0*pi*REAL(isn*(j-1)*(k-1))/n)
!              XXi(k)=XXi(k)+DATA(j)*SIN(-2.0*pi*REAL(isn*(j-1)*(k-1))/n)
!              XXr(k)=XXr(k)+DATA(j)*COS(-2.0*pi*(j-1)*(k-1)/n)
!              XXi(k)=XXi(k)+DATA(j)*SIN(-2.0*pi*(j-1)*(k-1)/n)
              XXr(k)=XXr(k)+DATA(j)*tempcos
              XXi(k)=XXi(k)+DATA(j)*tempsin
          enddo
          XXa(k) = sqrt(XXr(k)**2+XXi(k)**2)
!          if(isn .EQ. -1) then
!             XXa(k) = XXa(k)/n
!          endif
      enddo

      DATA(:) = XXa(:)


C      nn = INT(n/2)
C      j = 1
C
C      do i = 1,n,2
C            if(j .GT. i) then
C                  tempR = DATA(j)
C                  tempI = DATA(j+1)
C                  DATA(j) = DATA(i)
C                  DATA(j+1) = DATA(i+1)
C                  DATA(i) = tempR
C                  DATA(i+1) = tempI
C            endif
C            M = nn
C 1000       if((M .GE. 2) .AND. (j .GT. M)) then
C                  j = j-M
C                  M = M/2
C            GOTO 1000
C            endif
C            j = j + M
C      enddo
C      Mmax = 2
C 2000 if(n .GT. Mmax) then
C            Istep = 2*Mmax
C            theta = 2*pi/real(isn*Mmax)
C            WPR = -2.0*SIN(0.50*theta)**2.0
C            WPI = SIN(theta)
C            WR  = 1.0
C            WI  = 0.0
C            do M = 1,Mmax, 2
C                  do i = M, n, Istep
C                        j = i+Mmax
C                        tempR = WR*DATA(j) - WI*DATA(j+1)
C                        tempI = WR*DATA(j+1) + WI*DATA(j)
C                        DATA(j) = DATA(i) - tempR
C                        DATA(j+1) = DATA(i+1) - tempI
C                        DATA(i) = DATA(i) + tempR
C                        DATA(i+1) = DATA(i+1) + tempI
C                  enddo
C                  Wtemp = WR
C                  WR = WR*WPR-WI*WPI+WR
C                  WI = WI*WPR+Wtemp*WPI+WI
C            enddo
C            Mmax = Istep
C      GOTO 2000
C      endif



      end subroutine FFT
C --------------------------------------------------------------------
C --------------------------GOF------------------------------------
C Computes the "goodness of fit" measure
C IN1 = Synthetic Data
C IN2 = Reference Data
C --------------------------------------------------------------------
      subroutine GOF(IN1, IN2, nt, nts, ntf, SSQ)
      implicit none
      integer                :: nts, ntf, i, j, k, nt
      real                   :: SSQ, pi, vv, vv2, STDV
      real                   :: ss, ss2, StDevIN1, z0, tmp
      real                   :: StDevR, StDevIN2, rr, rr2
      real, dimension(nt)    :: IN1, IN2, rVec, IN3, IN4
      real, dimension(20001) :: ErrVec
      data pi/3.14159265358979323846/

      call ERFVEC(ErrVec)

      SSQ      = 0.0
      vv       = 0.0
      vv2      = 0.0
      rr       = 0.0
      rr2      = 0.0
      StDevR   = 0.0
      StDevIN1 = 0.0
      StDevIN2 = 0.0

C Sum the vectors, the residuals and the component wise squares
      do i = nts,ntf
          ss       = ss + ABS(IN1(i))
          ss2      = ss2 + ABS(IN1(i))**2.0
          vv       = vv + ABS(IN2(i))
          vv2      = vv2 + ABS(IN2(i))**2.0
          rr       = rr + ABS(IN1(i)-IN2(i))
          rr2      = rr2 + ABS(IN1(i)-IN2(i))**2.0
          rVec(i)  = (IN1(i)-IN2(i))
      enddo
C Calculate the Means and the Means of the squares
          ss  = ss/(ntf-nts+1.0)
          ss2 = ss2/(ntf-nts+1.0)
          vv  = vv/(ntf-nts+1.0)
          vv2 = vv2/(ntf-nts+1.0)
          rr  = rr/(ntf-nts+1.0)
          rr2 = rr2/(ntf-nts+1.0)
C Calculate the Standard Deviations
          StDevIN1 = SQRT(ss2-ss**2.0)
          StDevIN2 = SQRT(vv2-vv**2.0)
C         StDevR   = SQRT(rr2-rr**2.0)
C Generate z0
C New GOF Measure
      tmp = 0.0
      vv  = 0.0
      rr  = 0.0
      z0  = 0.0
      SSQ = 0.0
      do i = nts,ntf
           rr =  SQRT((IN1(i)-IN2(i))**2.0)
           vv =  (IN1(i)+IN2(i))/2.0
	       if(vv .EQ. 0.0)then
		      tmp = 0.0
           else
		      tmp = rr/vv
           endif
           z0 = z0 + ErrVec(INT(tmp/0.0001)+1)
      enddo
C          z0 = SQRT(rr/vv)
          SSQ = z0/(ntf-nts+1.0)
C          call ERF(z0, SSQ)
C          STDV = StDevR/StDevIN2
      end subroutine GOF
C --------------------------------------------------------------------
C --------------------------GOFFS------------------------------------
C Computes the "goodness of fit" measure for the Fourier Spectrum
C IN1 = Synthetic Data
C IN2 = Reference Data
C --------------------------------------------------------------------
      subroutine GOFFS(INS, INR, nt, dt, freq1, freq2, SSQ)
      implicit none
      integer                :: nts, ntf, i, j, k, nt, wind
      real                   :: SSQ, pi, vv, vv2, STDV, freq1, freq2
      real                   :: ss, ss2, StDevIN1, z0, tmp, dt
      real                   :: StDevR, StDevIN2, rr, rr2
      real                   :: tempR, tempI
      real, dimension(nt)    :: rVec, INR, INS, INTT
      real, dimension(nt/2)  :: IN1, IN2, IN3, IN4, tmp1, tmp2
      real, dimension(20001) :: ErrVec
      data pi/3.14159265358979323846/

      call ERFVEC(ErrVec)

      SSQ      = 0.0
      vv       = 0.0
      vv2      = 0.0
      rr       = 0.0
      rr2      = 0.0
      StDevR   = 0.0
      StDevIN1 = 0.0
      StDevIN2 = 0.0
      
      wind = INT(0.1*dt*nt+1)

      nts = INT(freq1*dt*nt+1)
      ntf = INT(freq2*dt*nt+1)

C calculate the FFT
C      INTT = INR
C      call FFT(INR, nt, 1)
C      do i = 1,nt/2
C            tempR = (i-1)*2+1
C            tempI = (i)*2
C            IN1(i) = SQRT((INR(tempR)**2.0)+(INR(tempI)**2.0))
C      enddo
C      INR = INTT
C      INTT = INS
C
C      call FFT(INS, nt, 1)
C      do i = 1,nt/2
C            tempR = (i-1)*2+1
C            tempI = (i)*2
C            IN2(i) = SQRT((INS(tempR)**2.0)+(INS(tempI)**2.0))
C      enddo
C      INS = INTT

      call FFT(INR, nt, 1)
      call FFT(INS, nt, 1)

      do i = 1,nt/2
            IN1(i) = INR(i)
            IN2(i) = INS(i)

            IN3(i) = IN1(i)
            IN4(i) = IN2(i)
      enddo
C Smooth out the data
      do i = wind+1,(nt/2)-wind
            tmp1(i) = 0.0
            tmp2(i) = 0.0
            do j = -1*wind,wind,1
                  tmp1(i) = tmp1(i) + (IN1(i+j))
                  tmp2(i) = tmp2(i) + (IN2(i+j))
            enddo
            IN3(i) = tmp1(i)/REAL(wind*2.0+1.0)
            IN4(i) = tmp2(i)/REAL(wind*2.0+1.0)
      enddo

C Sum the vectors, the residuals and the component wise squares
C      do i = nts,ntf
C          ss       = ss + ABS(IN1(i))
C          ss2      = ss2 + ABS(IN1(i))**2.0
C          vv       = vv + ABS(IN2(i))
C          vv2      = vv2 + ABS(IN2(i))**2.0
C          rr       = rr + ABS(IN1(i)-IN2(i))
C          rr2      = rr2 + ABS(IN1(i)-IN2(i))**2.0
C          rVec(i)  = (IN1(i)-IN2(i))
C      enddo
C Calculate the Means and the Means of the squares
C          ss  = ss/(ntf-nts+1.0)
C          ss2 = ss2/(ntf-nts+1.0)
C          vv  = vv/(ntf-nts+1.0)
C          vv2 = vv2/(ntf-nts+1.0)
C          rr  = rr/(ntf-nts+1.0)
C          rr2 = rr2/(ntf-nts+1.0)
C Calculate the Standard Deviations
C          StDevIN1 = SQRT(ss2-ss**2.0)
C          StDevIN2 = SQRT(vv2-vv**2.0)
C         StDevR   = SQRT(rr2-rr**2.0)
C Generate z0
C New GOF Measure
      tmp = 0.0
      vv  = 0.0
      rr  = 0.0
      z0  = 0.0
      SSQ = 0.0
      do i = nts,ntf
           rr =  SQRT((IN3(i)-IN4(i))**2.0)
           vv =  (IN3(i)+IN4(i))/2.0
	       if(vv .EQ. 0.0)then
		      tmp = 0.0
           else
		      tmp = rr/vv
           endif
           z0 = z0 + ErrVec(INT(tmp/0.0001)+1)
      enddo
C          z0 = SQRT(rr/vv)
          SSQ = z0/(ntf-nts+1.0)
C          call ERF(z0, SSQ)
C          STDV = StDevR/StDevIN2
      end subroutine GOFFS
C --------------------------------------------------------------------
C --------------------------NORM--------------------------------------
C example ---> x**(2.5) = (x**(2+0.5)) = (x**2)*sqrt(x)
C --------------------------------------------------------------------
      subroutine NORM(IN1, nt, nts, ntf, nm, pwr)
      implicit none
      integer                :: nts, ntf, i, j, k, nt
      real                   :: nm, pwr, pi
      real, dimension(nt)    :: IN1
      data pi/3.14159265358979323846/
     
      nm = 0.0
C Calculate the norm
      do i = nts,ntf
          nm = nm + (IN1(i))**pwr
      enddo
      
      nm = nm**(1.0/pwr)

      end subroutine NORM
C --------------------------------------------------------------------
C --------------------------SpecComp----------------------------------
C Compares the two spectral response series by a binned method
C --------------------------------------------------------------------
      subroutine SpecComp(SPECS, SPECR, COMPMDEV)
      implicit none
      integer                :: i, j, k
      real                   :: pi
      real, dimension(12)    :: COMPMDEV
      integer, dimension(13) :: indVec
      real, dimension(991)  :: SPECS, SPECR, perList
      data pi/3.14159265358979323846/
      data indVec/1,31,51,76,118,201,451,951,961,991,1041,1091,1141/

      do j = 1,901
           perList(j) = 0.001*real(j-1)+0.1
      enddo
      do j = 902,991
           perList(j) = 0.1*real(j-901)+1.0
      enddo
C      write(*,*) 'period 1 ','-period 2-', ' Mean Deviation in Bracket'
      do i = 1,12
         COMPMDEV(i) = 0.0
         do j = indVec(i),indVec(i+1)
              COMPMDEV(i) = COMPMDEV(i) + ABS(SPECS(j)-SPECR(j))
         enddo
         COMPMDEV(i) = COMPMDEV(i)/(indVec(i+1)-indVec(i)+1)
C         write(*,*) perList(indVec(i)),perList(indVec(i+1)), COMPMDEV(i)
C     +              ,SPECS(j)
      enddo

      end subroutine SpecComp
C --------------------------------------------------------------------
C --------------------------NGAComp----------------------------------
C Compares the two spectral response series
C --------------------------------------------------------------------
      subroutine NGAComp(NGA, SPEC, COMPMDEV)
      implicit none
      integer                :: i, j, k
      real                   :: pi
      real, dimension(16)    :: NGA, SPEC, COMPMDEV
      data pi/3.14159265358979323846/
     
      do i = 1,16         
         COMPMDEV(i) = (SPEC(i)-NGA(i))/NGA(i)
      enddo

      end subroutine NGAComp
C --------------------------------------------------------------------
C --------------------------organize----------------------------------
C Compares the two spectral response series
C --------------------------------------------------------------------
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
C --------------------------------------------------------------------
C-------------------------MODEL BIAS----------------------------------
C BB_GOF_Hartzell
C Model Bias Abramson Somerville "Uncertainty in Numerical Strong Ground motions .."
C --------------------------------------------------------------------
      subroutine MODBIAS(INPUTR, INPUTS, numSta, numPts, B)
      implicit none
      integer                       :: i, j, k, numSta, numPts
      real                          :: pi
      real, dimension(numPts)       :: B
      real, dimension(numSta,numPts):: INPUTR, INPUTS
      data pi/3.14159265358979323846/
      
      do i = 1,numPts
         B(i) = 0.0
      enddo
 
      do i = 1,numPts
         do j = 1,numSta
             B(i) = B(i) + LOG(INPUTR(j,i)/INPUTS(j,i))
         enddo 
         B(i) = B(i)/REAL(numSta)
      enddo

      end subroutine MODBIAS
C --------------------------------------------------------------------
C ------------------------MODEL ERROR---------------------------------
C BB_GOF_Hartzell
C Model Bias Abramson Somerville "Uncertainty in Numerical Strong Ground motions .."
C --------------------------------------------------------------------
      subroutine MODERROR(INPUTR, INPUTS, numSta, numPts, B, E)
      implicit none
      integer                       :: i, j, k, numSta, numPts
      real                          :: pi
      real, dimension(numPts)       :: B, E
      real, dimension(numSta,numPts):: INPUTR, INPUTS

      data pi/3.14159265358979323846/
      
      do i = 1,numPts
         E(i) = 0.0
      enddo
 
      do i = 1,numPts
         do j = 1,numSta
             E(i) = E(i) + (LOG(INPUTR(j,i)) - LOG(INPUTS(j,i)) 
     +              - B(i))**2.0  
         enddo 
         E(i) = SQRT(E(i)/REAL(numSta))
      enddo

      end subroutine MODERROR
C --------------------------------------------------------------------
C ----------------------------- ERF ----------------------------------
C Error function described by Scherbaum "On use of Response Spectral-Ref.."
C Returns a value between 0 and 1 ... one = perfect fit, zero error
C --------------------------------------------------------------------
      subroutine ERF(INPUT, OUTPUT)
      implicit none
      integer                       :: i, n
      real                          :: pi, INPUT, OUTPUT, t, z0
      real, dimension(10001)       :: eT, outT

      data pi/3.14159265358979323846/
      z0 = ABS(INPUT)
      OUTPUT = 0.0  
      eT(1) = EXP(-(z0**2.0))

      do i=2,10001
         t     = z0+(REAL(i-1)*0.001)
         eT(i) = EXP(-(t**2.0))
         OUTPUT = OUTPUT+(eT(i-1)*0.001)+
     +              (eT(i)-eT(i-1))*(0.5)*0.001
      enddo

      OUTPUT = OUTPUT*(2.0/SQRT(pi))
      end subroutine ERF
C --------------------------------------------------------------------
C ----------------------------- ERFVEC ----------------------------------
C Error function described by Scherbaum "On use of Response Spectral-Ref.."
C Returns a value between 0 and 1 ... one = perfect fit, zero error
C --------------------------------------------------------------------
      subroutine ERFVEC(OUTPUT)
      implicit none
      integer                      :: i, n
      real                         :: pi, gg, rr, delX
      real, dimension(20001)       :: OUTPUT, qq

      data pi/3.14159265358979323846/
      delX = 0.0001
      OUTPUT(1) = 1.0
      qq(1) = 0.0

      do i = 2,20001
         gg = REAL(i-1)*delX
         rr = EXP(-(gg**2.0))
         qq(i) = qq(i-1)+(2.0/SQRT(pi))*rr*delX
         OUTPUT(i) = 1.0 - qq(i)
      enddo
C     OUTPUT(20001) = 0.0  
      
      end subroutine ERFVEC
C --------------------------------------------------------------------
C ------------------- Cross - Correlation ----------------------------
C Compute the delay and the Cross - Correlation values
C Computes the "normalized" cross-correlation
C --------------------------------------------------------------------
      subroutine CROSSC(IN1, IN2, nt, dt, R, D)
      implicit none
      integer               :: i, nt, D, j, sft
      real                  :: pi, X, Y, XM, YM, XS, YS, X2, Y2, R
      real                  :: P1, P2, Dout, dt
      real, dimension(nt)   :: IN1, IN2, DEN, NUM
      real, dimension(2*nt) :: RR, IN3

      data pi/3.14159265358979323846/
     
      XM = 0.0
      YM = 0.0
      XS = 0.0
      P1 = 0.0
      P2 = 0.0

      do i=1,nt
         if(ABS(IN1(i)) .GT. P1) P1 = ABS(IN1(i))
         if(ABS(IN2(i)) .GT. P2) P2 = ABS(IN2(i))
         XM     = XM + IN1(i)
         YM     = YM + IN2(i)
         IN3(i) = IN2(i)
      enddo

C      XM = XM/nt
C      YM = YM/nt
      XM = 0.0
      YM = 0.0
      X2  = 0.0

      do i=nt+1,nt*2
         X2     = X2 + (IN1(i-nt)-XM)**2.0
         IN3(i) = IN2(i-nt)
      enddo
    
      XS = SQRT(X2)

C Calculate
      do i=0,nt-1
          Y2 = 0.0
          YS = 0.0
          NUM(i+1) = 0.0 
         do j=1+i,nt+i
            NUM(i+1) = NUM(i+1) + (IN1(j-i)-XM)*(IN3(j)-YM)
            Y2 = Y2 + (IN3(j)-YM)**2.0
         enddo
          YS       = SQRT(Y2) 
          DEN(i+1) = XS*YS
          RR(i+1)  = (NUM(i+1)/DEN(i+1))
      enddo

C calculate the cross-correlation without a shift
      X2 = 0.0
      Y2 = 0.0
      XM = 0.0
      do i=1,nt
         X2 = X2 + IN1(i)**2.0
         Y2 = Y2 + IN2(i)**2.0
         XM = XM + IN1(i)*IN2(i)
      enddo
      

      D = 0
      R = 0.0
      do i=1,nt
          if(RR(i) .GT. R) then
             D = (i-1)
          endif
      enddo

      IF(X2 .EQ. 0.0 .OR. Y2 .EQ. 0.0) then
          R = 0.0
      ELSE
          R = XM/(SQRT(X2)*SQRT(Y2))
      ENDIF
      IF(R .LT. 0.0) R = 0.0
      end subroutine CROSSC
C --------------------------------------------------------------------
C ------------------- Covariance ----------------------------
C Compute the covariance
C --------------------------------------------------------------------
      subroutine COVAR(IN1, IN2, nt, CV)
      implicit none
      integer               :: i, nt, D, j
      real                  :: pi, X, Y, XM, YM, CV
      real, dimension(nt)   :: IN1, IN2

      data pi/3.14159265358979323846/
     
      XM = 0.0
      YM = 0.0
      CV = 0.0

      do i=1,nt
         XM     = XM + IN1(i)
         YM     = YM + IN2(i)
      enddo

      XM = XM/nt
      YM = YM/nt

C Calculate
      do i=1,nt
          CV = CV + (IN1(i)-XM)*(IN2(i)-YM)
      enddo 
      CV = CV/nt

      end subroutine COVAR
C --------------------------------------------------------------------
C --------------------------RESIZE------------------------------------
C Subroutine to resize the seismograms (Spline)
C --------------------------------------------------------------------
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
C --------------------------------------------------------------------
C ---------------------------PPMCC------------------------------------
C Pearson Product-Moment Correlation Coefficiant
C Determines how well a linear function can be fit to the variables X and Y
C --------------------------------------------------------------------
      subroutine PPMCC(XX, YY, nt, rr)
      implicit none
      integer              :: nt, i
      real                 :: rr, sx, sy, Xm, Ym
      real, dimension(nt)  :: XX, YY

      Xm = 0.0
      Ym = 0.0
      sx = 0.0
      sy = 0.0
      rr = 0.0

      do i = 1,nt
            Xm = Xm + XX(i)
            Ym = Ym + YY(i)
      enddo

      Xm = Xm/nt
      Ym = Ym/nt

      do i = 1,nt
            sx = sx + (XX(i)-Xm)**2.0
            sy = sy + (YY(i)-Ym)**2.0
      enddo

      sx = SQRT(sx/real(nt-1))
      sy = SQRT(sy/real(nt-1))

      do i = 1,nt
            rr = rr + ((XX(i)-Xm)/sx)*((YY(i)-Ym)/sy)
      enddo

      rr = rr/(nt-1)

C sX = sum of X vec, sXY = sum X*Y vec, sX2 = sum X**2.0 
C      rho = (nt*sXY-sX*sY)/(SQRT(nt*sX2-sX**2)*SQRT(nt*sY2-sY**2))

      end subroutine PPMCC
C --------------------------------------------------------------------
C ---------------------------SRCC------------------------------------
C Spearman's Rank Correlation Coefficiant
C Determines how well a monotonic function can be fit to the variables X and Y
C --------------------------------------------------------------------
      subroutine SRCC(XX, YY, nt, rho)
      implicit none
      integer               :: nt, i, dd
      real               :: rho
      real, dimension(nt):: XX, YY
      integer, dimension(nt):: XI, YI

      call organize(XX, XI, nt)
      call organize(YY, YI, nt)

      dd  = 0.0
      rho = 0.0

      do i = 1,nt
            dd  = dd + (XI(i)-YI(i))**2.0
      enddo

      rho = 1.0 - 6.0*REAL(dd)/real(nt*(nt**2-1))

      end subroutine SRCC
C --------------------------------------------------------------------
C ---------------------------SDOTP------------------------------------
C Scaled Dot Product = Cos(Theta), Theta = angle between XX and YY vectors
C --------------------------------------------------------------------
      subroutine SDOTP(XX, YY, nt, CosT)
      implicit none
      integer               :: nt, i, dd
      real               :: CosT, SSqX, SSqY, SSqXY
      real, dimension(nt):: XX, YY

      SSqX  = 0.0
      SSqY  = 0.0
      SSqXY = 0.0
      CosT  = 0.0

      do i = 1,nt
            SSqX  = SSqX  + XX(i)**2.0
            SSqY  = SSqY  + YY(i)**2.0
            SSqXY = SSqXY + XX(i)*YY(i)
      enddo

      SSqX  = SQRT(SSqX)
      SSqY  = SQRT(SSQY)
      CosT  = SSqXY/(SSqX*SSqY)

      end subroutine SDOTP
C --------------------------------------------------------------------
C ---------------------------NGAROT------------------------------------
C Input a seismogram and get out the 17 values coresponding to the 50th 
C percentile and possibly the standard deviations associated
C --------------------------------------------------------------------
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
C --------------------------------------------------------------------

