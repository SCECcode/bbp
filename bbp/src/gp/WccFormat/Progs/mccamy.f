      SUBROUTINE LOPAS(M,F2,X,NL,IBACK)
C
C
C
C Subroutine LOPAS.
C
C    LOPAS is a Butterworth low pass reCursive filter written by
C Keith McCamy at Lamont-Doherty Geological Observatory, N.Y.  It is
C used in place of the Ormsby filter.
C
C    Call LOPAS (m, f2, x, nl, iback)
C
C        m = a flag;
C            if > 0, a filter is designed and new data used
C            if = 0, the last designed filter is used
C            if < 0, the last filter is used with new data
C        f2 = corner frequency in Nyquist units (fny=2*dt*fcps)
C        x = source data array
C        nl = the number of points in x
C        iback = a flag that indicates direction; 0 = forward,
C            1 = backward
C
C-
C
      DIMENSION  ALR(5), ALI(5), AO(5), A1(5), B1(5), B2(5), X(*)
      DIMENSION  Z(5), Z1(5), Z2(5)
C
      IF (M) 10, 2, 1
   10 N = -M
      GOTO 20
    1 N = M
      FN = N * 2
      RT = 3.1415926 / FN
      FI = -RT / 2.
      WW = X(1)
      DO 12 I = 1,N
         Z1(I) = 0.
         Z2(I) = 0.
         FI = FI + RT
         ALR(I) = -COS(FI)
   12    ALI(I) = SIN(FI)
      DO 122 K = 1,N
         PR = 0.
         PI = 2. * ALI(K)
         PRK = ALR(K)**2 - ALI(K)**2
         DO 1221 I = 1,N
            IF (I .EQ. K) GOTO 1221
            DR = PRK - ALR(K)*ALR(I)*2. + 1.
            DI = 2.*ALR(K)*ALI(K) - 2.*ALI(K)*ALR(I)
            TPR = PR*DR - PI*DI
            PI = PR*DI + PI*DR
            PR = TPR
 1221       CONTINUE
         DR = PR*PR + PI*PI
         AR = PR / DR
         AI = -PI / DR
         A1(K) = 2. * AR
         AO(K) = 2. * (ALR(K)*AR+ALI(K)*AI)
  122    B1(K) = -2. * ALR(K)
C...  HERE ENDS ALL GENERALIZED FORMS
      POO = F2 * 3.1415926/2.
      WU = 2. * SIN(POO)/COS(POO)
      DO 13 K = 1,N
         DR = 4./WU/WU + B1(K)*2./WU + 1.
         TA = (A1(K)*2./WU-AO(K)) / DR
         A1(K) = -(A1(K)*2./WU+AO(K)) / DR
         AO(K) = TA
         B2(K) = (4./WU/WU-B1(K)*2./WU+1.) / DR
   13    B1(K) = (2.-8./WU/WU) / DR
   20 CONTINUE
      WW = X(1)
      IF (IBACK .EQ. 1) WW = X(NL)
      DO 201 I = 1,N
         Z1(I) = 0.
  201    Z2(I) = 0.
    2 DO 22 II = 1,NL
         I = II
         IF (IBACK .EQ. 1) I = NL - II + 1
         ZZ = X(I) + WW
         XX = 0.
         WW = X(I)
         DO 21 J = 1,N
            Z(J) = ZZ - B1(J)*Z1(J) - B2(J)*Z2(J)
            XX = XX + AO(J)*Z(J) + A1(J)*Z1(J)
            Z2(J) = Z1(J)
   21       Z1(J) = Z(J)
   22    X(I) = XX
      PN = N
      FP = F2
      RETURN
      END
      SUBROUTINE HIPAS (M, F2, X, NL, IBACK)
C
C+
C
C Subroutine HIPAS.
C
C    HIPAS uses MCCamy's Butterworth high-pass recursive
C filter.
C
C    Call HIPAS (m, f2, x, nl, iback)
C
C        m = filter flag; if it = 0 then the last
C            designed filter is used; if it > 0 then
C            it designs a filter and starts with new data;
C            if it < 0 then it uses the last filter, but
C            starts on new data
C        f2 = the corner frequency in Nyquist units (fny =
C             2*dt*fcps)
C        x = the source array
C        nl = the number of points in x
C        iback = a flag that indicates direction; 0 = forward,
C            1 = backward
C
C-
C
      DIMENSION  ALR(5), ALI(5), AO(5), A1(5), B1(5), B2(5), X(*)
      DIMENSION  Z(5), Z1(5), Z2(5)
C
      IF(M) 10, 2, 1
C...  N = filter order (roll-off=N*12db/octave)
   10 N = -M
      GOTO 20
    1 N = M
      FN = N * 2
      RT = 3.1415926 / FN
      FI = -RT / 2.
      WW = X(1)
      DO 12 I = 1,N
         Z1(I) = 0.
         Z2(I) = 0.
         FI = FI + RT
         ALR(I) = -COS(FI)
   12    ALI(I) = SIN(FI)
      DO 122 K = 1,N
         PR = 0.
         PI = 2. * ALI(K)
         PRK = ALR(K)**2 - ALI(K)**2
         DO 1221 I = 1,N
            IF(I .EQ. K) GOTO 1221
            DR = PRK - 2.*ALR(K)*ALR(I) + 1.
            DI = 2.*ALR(K)*ALI(K) - 2.*ALI(K)*ALR(I)
            TPR = PR*DR - PI*DI
            PI = PR*DI + PI*DR
            PR = TPR
 1221       CONTINUE
         DR = PR*PR + PI*PI
         AR = PR / DR
         AI = -PI / DR
         A1(K) = 2. * AR
         AO(K) = 2. * (ALR(K)*AR+ALI(K)*AI)
  122    B1(K) = -2. * ALR(K)
C
C   HERE ENDS ALL GENERALIZED FORMS
C   BUILD HIGH PASS
      POO = F2 * 3.1415926 / 2.
      WU = 2. * SIN(POO) / COS(POO)
      DO 13 K = 1,N
         DR = WU*WU/4. + B1(K)*WU/2. + 1.
         TA = (A1(K)*WU/2.-AO(K)) / DR
         A1(K) = (A1(K)*WU/2.+AO(K)) / DR
         AO(K) = TA
         B2(K) = (WU*WU/4.-B1(K)*WU/2.+1.) / DR
   13    B1(K) = (WU*WU/2.-2.) / DR
   20 CONTINUE
      WW = X(1)
      IF (IBACK .EQ. 1) WW = X(NL)
      DO 201 I = 1,N
         Z1(I) = 0.
  201    Z2(I) = 0.
    2 DO 22 II = 1,NL
         I = II
         IF (IBACK .EQ. 1) I = NL - II + 1
         ZZ = X(I) - WW
         XX = 0.
         WW = X(I)
         DO 21 J = 1,N
            Z(J) = ZZ - B1(J)*Z1(J) - B2(J)*Z2(J)
            XX = XX + AO(J)*Z(J) + A1(J)*Z1(J)
            Z2(J) = Z1(J)
   21       Z1(J) = Z(J)
   22    X(I) = XX
      PN = N
      FP = F2
      RETURN
      END
      SUBROUTINE ZPLOP (M, FR2, X, NL)
C
C+
C
C Subroutine ZPLOP.
C
C    ZPLOP Calls McCamys low pass Butterworth routine 2 times,
C once forward and once with the time sequence reversed.
C
C    call ZPLOP (m, fr2, x, nl)
C
C        m  - integer roll-off width for high pass
C        fr2- corner frequency
C        x  - array that will be filtered, equally spaced in time
C        nl - the number of elements in x
C
C
C-
C
      DIMENSION  X(NL)
C
      IBACK=0
      CALL LOPAS(M,FR2,X,NL,IBACK)
      IBACK=1
      CALL LOPAS(M,FR2,X,NL,IBACK)
      RETURN
      END
      SUBROUTINE ZPHIP (M, FR2, X, NL)
C
C+
C
C Subroutine ZPHIP.
C
C    ZPHIP Calls McCamys high pass Butterworth routine 2 times,
C once forward and once with the time sequence reversed.
C
C    call ZPHIP (m, fr2, x, nl)
C
C        m  - integer roll-off width for high pass
C        fr2- corner frequency
C        x  - array that will be filtered, equally spaced in time
C        nl - the number of elements in x
C
C-
C
      DIMENSION  X(NL)
C
      IBACK = 0
      CALL HIPAS(M,FR2,X,NL,IBACK)
      IBACK = 1
      CALL HIPAS(M,FR2,X,NL,IBACK)
      RETURN
      END
