c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE realft(data,n,isign)
c
c  Calculates the Fourier transform of a set of real data points.
c  isign=-1, forward tranform of a real-valued data point.
c     The real-valued first (zero frequency) last (fmax) components
c     of the complex transform are returned as element data(1) and
c     data(2), respectively. n must be a opwer of 2.
c  isign=+1, inverse transform of a complex data array. The result
c     in this case must be multipled by 2/n.
c
      IMPLICIT NONE
      INTEGER isign,n
      REAL data(n)
CU    USES four1
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.-1) then
        c2=-0.5
        theta=-theta
        call four1(data,n/2,isign)
      else
        c2=0.5
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      data(n/2+2)=-data(n/2+2)
      if (isign.eq.-1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,isign)
      endif
      return
      END
C
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE four1(data,nn,isign)
      IMPLICIT NONE
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,5W"1.UY.
c======================================================================
c
      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      implicit NONE
      INTEGER isign,nn1,nn2,nn3
      COMPLEX data(nn1/2,nn2,nn3),speq(nn2,nn3)
c USES fourn
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX c1,c2,h1,h2,w
      c1=cmplx(0.5,0.0)
      c2=cmplx(0.0,-0.5*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do i3=1,nn3
        do i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
        enddo
        enddo
      endif
      do i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
          enddo
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=cmplx(sngl(wr),sngl(wi))
        enddo
      enddo
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END
c
c======================================================================
      SUBROUTINE fourn(data,nn,ndim,isign)
      implicit NONE
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do idim=1,ndim
        ntot=ntot*nn(idim)
      enddo
      nprev=1
      do idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do i1=i2,i2+ip1-2,2
              do i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
              enddo
            enddo
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
            goto 1
          endif
          i2rev=i2rev+ibit
        enddo
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do i3=1,ifp1,ip1
            do i1=i3,i3+ip1-2,2
              do i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
              enddo
            enddo
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
          enddo
          ifp1=ifp2
          goto 2
        endif
        nprev=n*nprev
      enddo
      return
      END
c
c    ===============================================
c    inverse cdf of rayleigh distribution
c    y=sigma*sqrt(-2ln(1-F)), F is cdf of rayleigh distribution
c
      Function irayl(f,sigma)
       implicit NONE
       real f,irayl,sigma
       if(f.ge.0.0.and.f.lt.1.0)then
         irayl=sigma*sqrt(-2.0*log(1.0-f))
       else
         write(*,*)"input error cumlative cdf must within (0,1)"
         write(*,*)"input F= ",f
         stop
       endif
      END
c
c======================================================================
c     rayleigh distribution: (x/sigma^2)exp(-x^2/(2sigma^2))
c     cdf is 1-exp(-x^2/(2sigma^2))
c     mean is 1.253*sigma
c     mode is sigma
c     95% of cdf at x~2.448*sigma
c
      Function rayleigh(idum,sigma)
      implicit NONE
      integer idum
      real rayleigh,u,sigma,ran3
      u=ran3(idum)
      rayleigh=sigma*sqrt(-2.0*log(u))
      return
      END
c
c======================================================================
      FUNCTION gasdev(idum,iset)
      implicit NONE
      INTEGER idum
      REAL gasdev
CU    USES ran3
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran3
      SAVE gset
      if (iset.eq.0) then
1       v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
c
c======================================================================
C
      FUNCTION gammln(xx)
      implicit NONE
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
      return
      END
c
c======================================================================
c
      FUNCTION betafn(z,w)
      implicit NONE
      REAL betafn,z,w
c  USES gammln
      REAL gammln
      betafn=exp(gammln(z)+gammln(w)-gammln(z+w))
      return
      END
c
c======================================================================
c calculate the cumulative probability density related to x
        subroutine gasprb(v2p,n)
c  gasprb=1.0-0.5*erfcc(x/sqrt(2))
        implicit NONE
        integer n,i
        real v2p(n),sum1,sum2,xx,avg,std,prb,x,z,t

        sum1=0.0
        sum2=0.0
        do i=1,n
          xx=v2p(i)
          sum1=sum1+xx
          sum2=sum2+xx*xx
        enddo
        avg=sum1/float(n)
        std=sqrt(sum2/float(n)-avg*avg)
        do i=1,n
          v2p(i)=(v2p(i)-avg)/std
        enddo

        do i=1,n
          x=v2p(i)
          z=abs(x*0.70710678)
          t=1./(1.+0.5*z)
          prb=0.5*t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
     *       t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
     *       t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
          if(x.gt.0.0) prb=1.-prb

          v2p(i)=prb
        enddo
        return
        END

c
c======================================================================
c
        FUNCTION ran3(idum)
c ROUTINE TO GENERATE A UNIFORMLY DISTRIBUTED RANDOM
c NUMBER ON THE INTERVAL [0,1].
c
        implicit NONE
        INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
        REAL ran3,AM,EPS,RNMX
        PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv(NTAB),iy
        SAVE iv,iy
        DATA iv /NTAB*0/, iy /0/
        if (idum.le.0.or.iy.eq.0) then
          idum=max(-idum,1)
          do 11 j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
11        continue
          iy=iv(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        j=1+iy/NDIV
        iy=iv(j)
        iv(j)=idum
        ran3=min(AM*iy,RNMX)
        return
        END

c
c======================================================================
c
c The following subroutine were provided by Dan O'Connell.
c
c======================================================================
c
      SUBROUTINE XZPTIM(time,slow,mxx,nx,nz,H,sx,sz,rtime,rtx,rtz,nreq)
      IMPLICIT NONE
      INTEGER  NX, NZ, ISZ, ISX, I, J, ILEFT, IRIGHT,
     &     ITOP, IBOT, mxx, ifirst, iimax, jimax, ltop,lbot, lleft,
     &     lright, isxm1, iszm1, isxp1,iszp1, irevc, irevr, ictop,
     &     icbot, idum,nreq
      INTEGER nx1,nz1, ifrstc,ifrstr
C
C ... source position and grid boundaries
      REAL SX,SZ,xmin,xmax,zmin,zmax, difmax, dist, testv, htime,
     &     dif,perc
C ... grid spacing and constants to reduce redundant calculations
      REAL H,H2,HTSQR2,HS2I,sval
C ... Variables of interest
      REAL TIME( mxx,nz), SLOW( mxx,nz ), rtime(1:*),rtx(1:*),rtz(1:*)
      H2 = H * H
      HTSQR2 = H * SQRT(2.0)
      HS2I = H / SQRT(2.0)
C
C ... Initialize travel-time array to big value
      DO I = 1, NZ
         DO J = 1, NX
            time(J,I) = 1.E10
         enddo
      enddo
C
C---- calculate grid boundaries
      xmin = 0.0
      zmin = 0.0
      nx1 = nx - 1
      nz1 = nz - 1
      xmax = float(nx1) * H
      zmax = float(nz1) * H
C
C .... Source position cannot be on or outside grid edge
      IF( SX .LT. XMIN .OR. SX .GT. XMAX .OR. SZ .LT. ZMIN
     &     .OR. SZ .GT. ZMAX) THEN
         PRINT *,' SOURCE OUTSIDE GRID.'
         stop
      endif
C
C---- CALCULATE CLOSEST GRID POINT TO SOURCE
      ISX = 1 + NINT((SX - XMIN)/H)
      ISZ = 1 + NINT((SZ - ZMIN)/H)
C
C ... Initialize source point using polar scheme out to 10 x-z rings
      CALL pintim(sx,sz,h,h,nx,nz,slow,mxx,time)
C
C ... Initialize starting positions so that corners of rows and columns
C ... overlap
      iszm1 = ISZ - 1
      ITOP = MAX(1,iszm1-1)
      iszp1 = ISZ + 1
      IBOT = min(nz,iszp1+1)
      isxm1 = ISX - 1
      ILEFT = max(isxm1-1,1)
      isxp1 = ISX + 1
      IRIGHT = min(nx,isxp1+1)
C
C---- set activity flags
      LLEFT = 1
      LRIGHT = 1
      LTOP = 1
      LBOT = 1
C
C---- check for edges that are already inactive due to source location
      IF (ILEFT .eq. 1) LLEFT = 0
      IF (IRIGHT .eq. nx) LRIGHT = 0
      IF (ITOP .eq. 1) LTOP = 0
      IF (IBOT .eq. nz) LBOT = 0
C
C---- calculate until the entire grid is filled with times
      do while (LLEFT .NE. 0 .OR. LRIGHT .NE. 0 .OR. LTOP .NE. 0 .OR.
     &     LBOT .NE. 0)
C
C ...... Initialize refraction flag to false
         ifirst = 0
         ictop = 0
         icbot = 0
C
C ...... Solve for column at left
         IF( LLEFT .NE. 0 ) THEN
            ifirst=1
            do while (ifirst .ne. 0)
               ifirst=0
               CALL COLTIM(time,slow,mxx,nx,nz,ILEFT,ITOP,IBOT,-1,
     &              0,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
C
C---------- check for refraction and the need to reverse the calculation
               if (ifirst .ne. 0) then
C
C---------------- fill out rows to give values to uninitiazed points
                  IF( LTOP .NE. 0 .and. ictop .ne. 0) THEN
                     idum = -1
                     CALL ROWTIM(time,slow,mxx,nx,nz,ITOP,ILEFT,IRIGHT,
     &                    -1,-1,idum,nx1,nz1,h,h2,HTSQR2,HS2I)
                     ictop = 0
                  endif
                  IF( LBOT .NE. 0 .and. icbot .ne. 0) THEN
                     idum = -1
                     CALL ROWTIM(time,slow,mxx,nx,nz,IBOT,ILEFT,IRIGHT,
     &                    1,-1,idum,nx1,nz1,h,h2,HTSQR2,HS2I)
                     icbot = 0
                  endif
C
C------------- reverse, starting column is the one just calculated
                  irevc = ileft + 1
                  do while (ifirst .ne. 0 .and. irevc .lt. isxm1)
                     ifirst = 0
                     call COLTIM(time,slow,mxx,nx,nz,irevc,ITOP,IBOT,1,
     &                    1,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
                     irevc = irevc+1
                  enddo
C
C------------- reverse again to get to starting position
                  do i=irevc-2,ileft,-1
                     call COLTIM(time,slow,mxx,nx,nz,i,ITOP,IBOT,-1,
     &                    -1,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
                  enddo
               endif
            enddo
C
C---------- if we have hit this edge set its flag inactive
            IF (ILEFT .eq. 1) LLEFT = 0
         ENDIF
C
C ...... Solve for column at right
         IF( LRIGHT .NE. 0 ) THEN
            ifirst=1
            do while (ifirst .ne. 0)
               ifirst=0
               CALL COLTIM(time,slow,mxx,nx,nz,iright,ITOP,IBOT,1,
     &              0,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
C
C------------- check for refraction and the need to reverse the calculation
               if (ifirst .ne. 0) then
C
C---------------- fill out rows to give values to uninitiazed points
                  IF( LTOP .NE. 0 .and. ictop .ne. 0) THEN
                     idum = -1
                     CALL ROWTIM(time,slow,mxx,nx,nz,ITOP,ILEFT,IRIGHT,
     &                    -1,-1,idum,nx1,nz1,h,h2,HTSQR2,HS2I)
                     ictop = 0
                  endif
                  IF( LBOT .NE. 0 .and. icbot .ne. 0) THEN
                     idum = -1
                     CALL ROWTIM(time,slow,mxx,nx,nz,IBOT,ILEFT,IRIGHT,
     &                    1,-1,idum,nx1,nz1,h,h2,HTSQR2,HS2I)
                     icbot = 0
                  endif
c                  print *,'reverse right ',iright
C
C---------------- reverse, starting column is the one just calculated
                  irevc = iright-1
                  do while (ifirst .ne. 0 .and. irevc .gt. isxp1)
                     ifirst = 0
                     call COLTIM(time,slow,mxx,nx,nz,irevc,ITOP,IBOT,-1,
     &                    1,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
c                     print *,' reversing right at ',iright,irevc
                     irevc = irevc-1
                  enddo
C
C---------------- reverse again to get to starting position
                  do i=irevc+2,iright
                     call COLTIM(time,slow,mxx,nx,nz,i,ITOP,IBOT,1,
     &                    -1,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
                  enddo
               endif
            enddo
C
C---------- if we have hit this edge set its flag inactive
            IF (IRIGHT .eq. nx) LRIGHT = 0
         ENDIF
C
C ...... Solve row at top
         IF( LTOP .NE. 0 ) THEN
            ifirst=1
            do while (ifirst .ne. 0)
               ifirst=0
               CALL ROWTIM(time,slow,mxx,nx,nz,ITOP,ILEFT,IRIGHT,-1,
     &              0,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
C
C---------- check for refraction and the need to reverse the calculation
               if (ifirst .ne. 0) then
C
C------------- reverse, starting row is the one just calculated
                  irevr = itop+1
                  do while (ifirst .ne. 0 .and. irevr .lt. iszm1)
                     ifirst = 0
                     CALL ROWTIM(time,slow,mxx,nx,nz,irevr,ILEFT,IRIGHT,
     &                    1,1,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
                     irevr = irevr+1
                  enddo
C
C------------- reverse again to get to starting position
                  do i=irevr-2,itop,-1
                     CALL ROWTIM(time,slow,mxx,nx,nz,i,ILEFT,IRIGHT,-1,
     &                    -1,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
                  enddo
               endif
            enddo
C
C---------- if we have hit this edge set its flag inactive
            IF (ITOP .eq. 1) LTOP = 0
         ENDIF
C ...... Solve row at bottom
         IF( LBOT .NE. 0 ) THEN
            ifirst=1
            do while (ifirst .ne. 0)
               ifirst=0
               CALL ROWTIM(time,slow,mxx,nx,nz,IBOT,ILEFT,IRIGHT,1,
     &              0,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
C
C------------- check for refraction and the need to reverse the calculation
               if (ifirst .ne. 0) then
c                  print *,'reverse bot ',ibot
C
C---------------- reverse, starting row is the one just calculated
                  irevr = ibot-1
                  do while (ifirst .ne. 0 .and. irevr .gt. iszp1)
                     ifirst = 0
                     CALL ROWTIM(time,slow,mxx,nx,nz,irevr,ILEFT,IRIGHT,
     &                    -1,1,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
c                     print *,'reversing bot at',ibot,irevr
                     irevr = irevr-1
                  enddo
C
C---------------- reverse again to get to starting position
                  do i=irevr+1,ibot
                     CALL ROWTIM(time,slow,mxx,nx,nz,i,ILEFT,IRIGHT,1,
     &                    -1,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
                  enddo
               endif
            enddo
C
C---------- if we have hit this edge set its flag inactive
            IF (IBOT .eq. nz) LBOT = 0
         ENDIF

C ...... Increment rows and columns that have advanced or set inactive flags
         IF( ILEFT .GT. 1 ) ILEFT = ILEFT - 1
         IF( IRIGHT .lt. nx ) IRIGHT = IRIGHT + 1
         IF( iTOP .gt. 1 ) ITOP = ITOP - 1
         IF( iBOT .lt. nz ) IBOT = IBOT + 1
C
C ... Continue calculation until all grid boundaries are found
      enddo
C
C ... The grid is now filled with travel times
c      print *,' entire grid filled with travel times.'
C
C---- compute requested travel-times
      call inttrt(h,time,slow,mxx,nx,nz,sx,sz,rtx,rtz,rtime,nreq)
C
C ... write out travel-time field
c      open(8,file='ttest.bin',status='unknown',form='unformatted')
c      write(8) nx,nz,1./slow(1,1),h
c      write(8) ((time(j,i),j=1,nx),i=1,nz)
c      close(8)
      return
      END
c
c======================================================================
c
      subroutine pintim(xsrc,zsrc,dx,dz,nx,nz,sloxz,mxx,timexz)
C
C     PINTIM calculates travel times on a polar grid centered on the source
C     position within 2-D cartesian grid of slowesses and interpolates the
C     calculated times onto the portion of the 2-D cartesian grid near the
C     source which provides an accurate initiazation to start the cartesian
C     finite difference propagation of travel times. The parameter SISIZE
C     determines how far into the cartesian grid the polar grid expands from
C     the source point. If SISIZE is increased, MXP (maximum number of angles)
C     and MXR (maximum number of radii) must be increased. As r increases
C     the angle increment must decrease to ensure at least 2 angle increments
C     within the outermost cartesian cell. These calculations are done in
C     in subroutine ANGRGN.
C
C     Only forward propagation is used since the beginning of head waves
C     tends to be picked up properly in a polar geometry. Nearly the
C     entire source region is re-computed in the cartesian finite difference
C     portion anway, and it does do reverersals. PINTIM calculates the
C     spherical wavefronts more accurately near the source, providing more
C     accurate (earlier) times at 45 degrees from the cartesian axes.
C
      IMPLICIT NONE
      INTEGER MXP,MXR
      PARAMETER (MXP = 200, MXR = 32)
      integer nx,nz,mxx
      DOUBLE PRECISION dr,dp,dr2,dr3,dp2,dp2dr2,dp2i,dridpi,drdp22,
     &     drdp,TIME(MXP,MXR),p(1:mxp),cosp,
     &     sang(mxp),angstr,angend,cang(mxp),sring(0:mxp,2),
     &     bigslo,r(0:mxr+1)
      parameter (bigslo = 20.D0)
      real xmin,xmax,zmin,zmax,xsrc,zsrc,dx,dz,sloxz(mxx,nz),
     &     timexz(mxx,nz),SISIZE
      PARAMETER (SISIZE = 10.)
      INTEGER I,J,J1,rftest,nr,np,iangr
C
C---- setup the limits of the cartesian region to intialize
C
      xmin = max(0.0,xsrc - SISIZE*dx)
      zmin = max(0.0,zsrc - SISIZE*dz)
      xmax = min(float(nx-1)*dx,xsrc + SISIZE*dx)
      zmax = min(float(nz-1)*dz,zsrc + SISIZE*dz)
C
C---- setup geometry from cartesion to polar travel-time grid
      call angrgn(xmin,xmax,zmin,zmax,dx,dz,xsrc,zsrc,angstr,angend,
     &     dr,dp,nr,np,r,p,cang,sang)
C
C---- set flag if complete circles are not needed to the endpoints
C---- of each ring see large slowness
      if (abs(angend-angstr) .lt. 4.0D0) THEN
         iangr = 1
      else
         iangr = 0
      endif
C      print *,'nr,np',nr,np
C
C---- call setrng to calculate slownesses for the first two radius rings
      call setrng(xmin,xmax,zmin,zmax,xsrc,zsrc,bigslo,angstr,angend,
     &     dx,dz,dr,dp,nr,np,r,p,cang,sang,sring,mxp,sloxz,mxx,
     &     nx,nz,1,1,iangr)
C
C---- PRE-CALCULATE angles, etc
      cosp=cos(dp)
C
C---- PRE-CALCULATE global grid constants
      call defcon(dr,dp,dr2,dr3,dp2,drdp,dp2dr2,dp2i,dridpi,drdp22)
C
C---- INITIALIZE TRAVEL-TIME FIELD TO INFINITY (almost)
      DO I = 1, nr
         DO J = 1, np
            TIME(J,I) = 1.D10
         enddo
      enddo
C
C---- SETUP SOURCE
C
C---- CALCULATE RADIAL TRAVEL-TIMES
      DO I=1,np
         TIME(I,1) = dr*min(sring(i-1,1),sring(i,1))
      ENDDO
C
C---- calculate possible refraction times
C---- sweep to maximum angle
      rftest = 0
      CALL REFSWP(1,np-1,TIME,MXP,nr,1,drdp,sring(1,1),
     &     sring(1,2),1,rftest)
C---- sweep to minimum angle
      CALL REFSWP(np,2,TIME,MXP,nr,1,drdp,sring(0,1),
     &     sring(0,2),-1,rftest)
C
C---- end source initialization, start ring loop
C
      DO J = 2, nr
         J1 = J - 1
C
C------- call setrng to calculate slownesses for the next radius ring
         call setrng(xmin,xmax,zmin,zmax,xsrc,zsrc,bigslo,angstr,angend,
     &        dx,dz,dr,dp,nr,np,r,p,cang,sang,sring,mxp,sloxz,mxx,
     &        nx,nz,j,0,iangr)
C
C------- calculate forward
         call forwrd(r(j),dr,dp,dr2,dr3,dp2,dp2dr2,dp2i,dridpi,drdp22,
     &        drdp,cosp,sring,time,mxp,nr,j,j1,rftest,np)
      enddo
C
C---- compute times of nodes of cartesian grid that fall within the
C---- initialization polar grid
      call intxzt(dr,dp,r,p,time,mxp,np,nr,xsrc,zsrc,timexz,mxx,nx,nz,
     &     dx,dz,xmin,xmax,zmin,zmax,iangr)
      return
      end
c
c======================================================================
c
      SUBROUTINE COLTIM(time,slow,mxx,nx,nz,NUMCOL,ITOP,IBOT,
     &     IWAY,irev,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
C-----------------------------------------------------------------
C     COLTIM calculates travel times for a column by sweeping up
C     and down the column for all templates until all points
C     in the column are filled with minimum travel times.
C-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER MIDIM
      PARAMETER (MIDIM = 100)
      INTEGER I, J, J1,K,K1, IWAY, NUMCOL, ITOP, IBOT, IREFS, ifirst,
     &     ILASTT, INEXTS, mxx, nx, nz,nx1,nz1, IL, IU,irev,ihit, isave,
     &     nmin,minstr(MIDIM),mnlend(MIDIM), mnuend(MIDIM),indexm(MIDIM)
      REAL h,h2,HTSQR2,HS2I,S1,S2
C .... Variables of interest (memory hogs)
      REAL TIME( mxx, nz), SLOW( mxx, nz ),tmmin(MIDIM)
C
C---- calculate inside row positions for time and slowness
      ILASTT = NUMCOL - IWAY
C
C---- make sure slowness index matches column position and increment direction
      IF (IWAY .LT. 0) THEN
         INEXTS = MAX(1,MIN(NUMCOL,nx1))
      ELSE
         INEXTS = MAX(1,MIN(ILASTT,NX1))
      ENDIF
C
C---- extend base to outer end to ensure a complete base to step forward
      ihit = 0
      IF (IREV .eq. 0) THEN
C
C------- find positions of time minima to order refraction loops below
         call fndcmn(time,mxx,nz,ILASTT,itop,ibot,nmin,minstr,mnlend,
     &        mnuend,tmmin,indexm)
         IREFS = MAX(1,MIN(INEXTS-IWAY,nx1))
         isave = 0
C
C--------loop through number minimum times found
         do i = 1,nmin
C
C---------- loop down
            do j = minstr(indexm(I)),mnuend(indexm(I))
               j1 = J + 1
               IU = MIN(J,nz1)
               ihit = 0
               call reftim(time(ILASTT,j),time(ILASTT,j1),
     &              slow(INEXTS,IU),slow(IREFS,IU),H,ihit)
C
C---------- if a refraction occurs sweep to the end of this column
               if (ihit .ne. 0) then
                  ifirst = 1
                  do k=j+1,ibot-1
                     k1 = k+1
                     IU = MIN(k,nz1)
                     call reftim(time(ILASTT,k),time(ILASTT,k1),
     &                    slow(INEXTS,IU),slow(IREFS,IU),H,ihit)
                  enddo
               endif
            enddo
C
C---------- bottom to top
            do j = minstr(indexm(I)), mnlend(indexm(I)),-1
               j1 = j - 1
               IU = MIN(J1,nz1)
               ihit = 0
               call reftim(time(ILASTT,j),time(ILASTT,j1),
     &              slow(INEXTS,IU),slow(IREFS,IU),H,ihit)
C
C------------- if a refraction occurs sweep to the end of this column
               if (ihit .ne. 0) then
                  ifirst = 1
                  do k=j-1,itop+1,-1
                     k1 = k-1
                     IU = MIN(k1,nz1)
                     call reftim(time(ILASTT,k),time(ILASTT,k1),
     &                    slow(INEXTS,IU),slow(IREFS,IU),H,ihit)
                  enddo
               endif
            enddo
         enddo
C
C------- if a refractor was found, return to do reverse propagation
         if (ifirst .ne. 0) return
      endif
C
      ihit = 0
C
C---- find positions of time minima to order transmission loops below
      call fndcmn(time,mxx,nz,ILASTT,itop,ibot,nmin,minstr,mnlend,
     &     mnuend,tmmin,indexm)
C
C---- 2d refraction outer slowness index
      IREFS = MAX(1,MIN(INEXTS+IWAY,NX1))
C
C---- sweep out from minimum time(s)
      do i = 1,nmin
C
C ...... calculate 1D transmission time at minima
         J = minstr(indexm(I))
         IL = MAX(1,J-1)
         IU = MIN(J,nz1)
C
C------- kludge to make sure a lower time is indicated during a reversal
         IF (slow(INEXTS,IU) .lt. slow(INEXTS,IL)) then
            s1 = slow(INEXTS,IU)
            s2 = slow(INEXTS,IL)
         else
            s1 = slow(INEXTS,IL)
            s2 = slow(INEXTS,IU)
         endif
         call reftim(TIME(ILASTT,J),TIME(NUMCOL,J),s1,s2,H,ihit)
C
C ...... start up from the bottom
         do j = minstr(indexm(I)), mnlend(indexm(I)),-1
            J1 = J - 1
            IL = MAX(1,J1-1)
            IU = max(1,J1)
C
C---------- kludge to make sure a lower time is indicated during a reversal
            IF (slow(INEXTS,IU) .lt. slow(INEXTS,IL)) then
               s1 = slow(INEXTS,IU)
               s2 = slow(INEXTS,IL)
            else
               s1 = slow(INEXTS,IL)
               s2 = slow(INEXTS,IU)
            endif
C
C---------- calculate refraction orthogonal to column
            call reftim(TIME(ILASTT,J1),TIME(NUMCOL,J1),s1,s2,H,ihit)
C
C---------- 2d diffraction
            call diftim(TIME(ILASTT,J),TIME(NUMCOL,J1),slow(INEXTS,IU),
     &           HTSQR2,IHIT)
C
C---------- 2d transmission from base
            call trntim(TIME(ILASTT,J),TIME(ILASTT,J1),TIME(NUMCOL,J),
     &           TIME(NUMCOL,J1),slow(INEXTS,IU),h2,hs2i,ihit,irev)
C
C---------- 2d transmission from side
            call trntim(TIME(ILASTT,J),TIME(NUMCOL,J),TIME(ILASTT,J1),
     &           TIME(NUMCOL,J1),slow(INEXTS,IU),h2,hs2i,ihit,irev)
            isave = 0
C
C---------- edge refraction
            call reftim(time(NUMCOL,j),time(NUMCOL,j1),
     &              slow(IREFS,IU),slow(INEXTS,IU),H,isave)
C
C---------- if a refraction occurs sweep to the end of this column
            if (isave .ne. 0) then
               ifirst = 1
               do k=j-1,itop+1,-1
                  k1 = k-1
                  IU = max(1,k1)
                  call reftim(time(NUMCOL,k),time(NUMCOL,k1),
     &                 slow(IREFS,IU),slow(INEXTS,IU),H,ihit)
               enddo
            endif
         enddo
         do j = minstr(indexm(I)),mnuend(indexm(I))
            J1 = J + 1
            IL = MIN(j1,nz1)
            IU = MIN(J,nz1)
C
C---------- kludge to make sure a lower time is indicated during a reversal
            IF (slow(INEXTS,IU) .lt. slow(INEXTS,IL)) then
               s1 = slow(INEXTS,IU)
               s2 = slow(INEXTS,IL)
            else
               s1 = slow(INEXTS,IL)
               s2 = slow(INEXTS,IU)
            endif
C
C---------- calculate refraction orthogonal to column
            call reftim(TIME(ILASTT,J1),TIME(NUMCOL,J1),s1,s2,H,ihit)
C
C---------- 2d diffraction
            call diftim(TIME(ILASTT,J),TIME(NUMCOL,J1),slow(INEXTS,IU),
     &           HTSQR2,IHIT)
C
C---------- 2d transmission from base
            call trntim(TIME(ILASTT,J),TIME(ILASTT,J1),TIME(NUMCOL,J),
     &           TIME(NUMCOL,J1),slow(INEXTS,IU),h2,hs2i,ihit,irev)
C
C---------- 2d transmission from side
            call trntim(TIME(ILASTT,J),TIME(NUMCOL,J),TIME(ILASTT,J1),
     &           TIME(NUMCOL,J1),slow(INEXTS,IU),h2,hs2i,ihit,irev)
            isave = 0
            call reftim(time(NUMCOL,j),time(NUMCOL,j1),
     &           slow(IREFS,IU),slow(INEXTS,IU),H,isave)
C
C---------- if a refraction occurs sweep to the end of this column
            if (isave .ne. 0) then
               ifirst = 1
               do k=j+1,ibot-1
                  k1 = k+1
                  IU = MIN(k,nz1)
                  call reftim(time(NUMCOL,k),time(NUMCOL,k1),
     &                 slow(IREFS,IU),slow(INEXTS,IU),H,ihit)
               enddo
            endif
         enddo
      enddo
C
C---- if in the process of reversing, indicate if any times were reduced
      if (irev .ne. 0) ifirst = max(ifirst,ihit)
      return
      END
c
c======================================================================
c
      SUBROUTINE rowTIM(time,slow,mxx,nx,nz,NUMrow,ileft,iright,
     &     IWAY,irev,ifirst,nx1,nz1,h,h2,HTSQR2,HS2I)
C-----------------------------------------------------------------
C     ROWTIM calculates travel times for a row by sweeping left
C     and right across the row for all templates until all points
C     in the row are filled with minimum travel times.
C-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER MIDIM
      PARAMETER (MIDIM = 100)
      INTEGER I,J, J1,K,K1, IWAY, numrow, ileft, iright, IREFS, IL, IU,
     &     ILASTT, INEXTS, mxx, nx, nz,nx1,nz1, ifirst, irev,ihit,isave,
     &     nmin,minstr(MIDIM),mnlend(MIDIM), mnuend(MIDIM),indexm(MIDIM)
      REAL h,h2,HTSQR2,HS2I,S1,S2
C ... Variables of interest (memory hogs)
      REAL TIME( mxx, nz), SLOW( mxx, nz ),tmmin(MIDIM)
C
C---- calculate row positions for time and slowness
      ILASTT = NUMROW - IWAY
C
C---- make sure slowness index matches row position and increment direction
      IF (IWAY .LT. 0) THEN
         INEXTS = MAX(1,MIN(numrow,nz1))
      ELSE
         INEXTS = MAX(1,MIN(ILASTT,Nz1))
      ENDIF
C
C---- extend base to bottom of corners to complete ensure a complete base
      ihit = 0
      IF (IREV .LE. 0) THEN
C
C------- find positions of time minima to order refraction loops below
         call fndrmn(time,mxx,nz,ILASTT,ileft,iright,nmin,minstr,mnlend,
     &        mnuend,tmmin,indexm)
         IREFS = MAX(1,MIN(INEXTS-IWAY,NZ1))
         isave = 0
C
C------- loop through number minimum times found
         do i = 1,nmin
C
C---------- left to right
            do j = minstr(indexm(I)),mnuend(indexm(I))
               j1 = J + 1
               IU = MIN(J,nx1)
               ihit = 0
               call reftim(time(j,ILASTT),time(j1,ILASTT),
     &              slow(iu,INEXTS),slow(iu,IREFS),H,ihit)
C
C------------- if a refraction occurs sweep to the end of this row
               if (ihit .ne. 0) then
                  ifirst = 1
                  do k=j+1,iright-1
                     k1 = k+1
                     IU = MIN(k,nx1)
                     call reftim(time(k,ILASTT),time(k1,ILASTT),
     &                    slow(iu,INEXTS),slow(iu,IREFS),H,ihit)
                  enddo
               endif
            enddo
C
C---------- right to left
            do j = minstr(indexm(I)), mnlend(indexm(I)),-1
               j1 = j - 1
               IU = MIN(J1,nx1)
               ihit = 0
               call reftim(time(j,ILASTT),time(j1,ILASTT),
     &              slow(iu,INEXTS),slow(iu,IREFS),H,ihit)
C
C------------- if a refraction occurs sweep to the end of this row
               if (ihit .ne. 0) then
                  ifirst = 1
                  do k=j-1,ileft+1,-1
                     k1 = k-1
                     IU = MIN(k,nx1)
                     call reftim(time(k,ILASTT),time(k1,ILASTT),
     &                    slow(iu,INEXTS),slow(iu,IREFS),H,ihit)
                  enddo
               endif
            enddo
         enddo
C
C------- if a refractor was found, return to do reverse propagation
         if (ifirst .ne. 0) return
      ENDIF
      ihit = 0
C
C---- find positions of time minima to order transmission loops below
      call fndrmn(time,mxx,nz,ILASTT,ileft,iright,nmin,minstr,mnlend,
     &     mnuend,tmmin,indexm)
C
C---- 2d refraction outer slowness index
      IREFS = MAX(1,MIN(INEXTS+IWAY,NZ1))
C
C---- sweep out from minimum time(s)
      do i = 1,nmin
C
C ...... calculate 1D transmission time
         J = minstr(indexm(i))
         IL = MAX(1,J-1)
         IU = MIN(J,nx1)
C
C------- kludge to make sure a lower time is indicated during a reversal
         if (slow(iu,INEXTS) .lt. slow(il,INEXTS)) then
            s1 = slow(iu,INEXTS)
            s2 = slow(il,INEXTS)
         else
            s2 = slow(iu,INEXTS)
            s1 = slow(il,INEXTS)
         endif
         call reftim(TIME(j,ilastt),TIME(j,numrow),s1,s2,H,ihit)
C
C------- right to left
         do j = minstr(indexm(I)), mnlend(indexm(I)),-1
            j1 = j - 1
            IL = max(1,j1-1)
            IU = max(1,J1)
C
C---------- kludge to make sure a lower time is indicated during a reversal
            if (slow(iu,INEXTS) .lt. slow(il,INEXTS)) then
               s1 = slow(iu,INEXTS)
               s2 = slow(il,INEXTS)
            else
               s2 = slow(iu,INEXTS)
               s1 = slow(il,INEXTS)
            endif
C
C---------- calculate refraction orthogonal to row
            call reftim(TIME(j1,ilastt),TIME(j1,numrow),s1,s2,H,ihit)
C
C---------- 2d diffraction
            call diftim(time(j,ilastt),time(j1,numrow),slow(iu,inexts),
     &           HTSQR2,IHIT)
C
C---------- 2d transmission from base
            call trntim(TIME(j,ilastt),TIME(j1,ilastt),TIME(j,numrow),
     &           TIME(j1,numrow),slow(iu,INEXTS),h2,hs2i,ihit,irev)
C
C---------- 2d transmission from side
            call trntim(TIME(j,ilastt),TIME(j,numrow),TIME(j1,ilastt),
     &           TIME(j1,numrow),slow(iu,INEXTS),h2,hs2i,ihit,irev)
            isave = 0
C
C---------- edge refraction
            call reftim(time(j,numrow),time(j1,numrow),
     &           slow(iu,IREFS),slow(iu,INEXTS),H,isave)
C
C---------- if a refraction occurs sweep to the end of this row
            if (isave .ne. 0) then
               ifirst = 1
               do k=j-1,ileft+1,-1
                  k1 = k-1
                  IU = max(1,k1)
                  call reftim(time(k,numrow),time(k1,numrow),
     &                 slow(iu,IREFS),slow(iu,INEXTS),H,ihit)
               enddo
            endif
         enddo
C
C------- sweep from left to right
         do j = minstr(indexm(I)),mnuend(indexm(I))
            J1 = J + 1
            IL = min(j1,nx1)
            IU = MIN(J,nx1)
C
C---------- kludge to make sure a lower time is indicated during a reversal
            if (slow(iu,INEXTS) .lt. slow(il,INEXTS)) then
               s1 = slow(iu,INEXTS)
               s2 = slow(il,INEXTS)
            else
               s2 = slow(iu,INEXTS)
               s1 = slow(il,INEXTS)
            endif
C
C---------- calculate refraction orthogonal to row
            call reftim(TIME(j1,ilastt),TIME(j1,numrow),s1,s2,H,ihit)
C
C---------- 2d diffraction
            call diftim(TIME(j,ilastt),TIME(j1,numrow),slow(IU,INEXTS),
     &           HTSQR2,IHIT)
C
C---------- 2d transmission from base
            call trntim(TIME(j,ilastt),TIME(j1,ilastt),TIME(j,numrow),
     &           TIME(j1,numrow),slow(iu,INEXTS),h2,hs2i,ihit,irev)
C
C---------- 2d transmission from side
            call trntim(TIME(j,ilastt),TIME(j,numrow),TIME(j1,ilastt),
     &           TIME(j1,numrow),slow(iu,INEXTS),h2,hs2i,ihit,irev)
            isave = 0
            call reftim(time(j,numrow),time(j1,numrow),
     &           slow(iu,IREFS),slow(iu,INEXTS),H,isave)
C
C---------- if a refraction occurs sweep to the end of this row
            if (isave .ne. 0) then
               ifirst = 1
               do k=j+1,iright-1
                  k1 = k+1
                  IU = MIN(k,nx1)
                  call reftim(time(k,numrow),time(k1,numrow),
     &                 slow(iu,IREFS),slow(iu,INEXTS),H,ihit)
               enddo
            endif
         enddo
      enddo
C
C---- if we are reversing the calculation assign minimum time hit flag
      if (irev .ne. 0) ifirst = max(ihit,ifirst)
      RETURN
      END
c
c======================================================================
c
      subroutine inttrt(h,time,slow,mxx,nx,nz,srcx,srcz,
     &     xrec,zrec,times,nrec)
C
C     inttim calculates times for requested x,z receiver positions
C     using bilinear interpolation of the cartesian travel-time field
C
      implicit none
      integer nrec,mxx,nx,nz,i,ix,iz,iz1,ix1,ixs,izs
      real time(mxx,nz),slow(mxx,nz),h,hi,za,zb,onemfx,fx,fz,sro
      real xrec(1:*),zrec(1:*),times(1:*),srcx,srcz
C
C---- reduce some divisions
      hi = 1.0/h
C
C---- calculate source position
      ixs = min(1 + int(srcx*hi),nx-1)
      izs = min(1 + int(srcz*hi),nz-1)
C
C---- loop through number of receivers
      do i = 1,nrec
C
C------- calculate index position on cartesian grid
         ix = min(1 + int(xrec(i)*hi),nx-1)
         iz = min(1 + int(zrec(i)*hi),nz-1)
C
C------- do bilinear interpolation
         fx = (xrec(i)-float(ix-1)*h)*hi
         fz = (zrec(i)-float(iz-1)*h)*hi
         ix1 = ix+1
         iz1 = iz+1
         onemfx = 1.0-fx
         za = fx*time(ix1,iz) + time(ix,iz)*onemfx
         zb = fx*time(ix1,iz1) + time(ix,iz1)*onemfx
         times(i) = fz*zb + za*(1.0-fz)
C
C------- check for receiver in source cell and check for minimum time
         if (ixs .eq. ix .and. izs .eq. iz) then
            sro = sqrt((xrec(i)-srcx)**2 + (zrec(i)-srcz)**2)
C
C---------- take min of interpolated time and direct arrival time
            times(i) = min(times(i),sro*slow(ixs,izs))
         endif
      enddo
      return
      end
c
c======================================================================
c
      subroutine angrgn(xmin,xmax,zmin,zmax,dx,dz,xsrc,zsrc,angstr,
     &     angend,dr,dp,nr,np,r,p,cang,sang)
C
C     angrgn determines the minimum and maximum angles required to fit
C     the circular travel-time grid to the rectangular velocity system
C     and resonable values of polar variable incremements
C
      implicit none
      integer nr,np,i
      double precision angstr,angend,dr,dp,dxm,dzm,c,c2,bx,bz,ax,az,
     &     r(0:*),p(1:*),cang(1:*),sang(1:*)
      real xmin,xmax,zmin,zmax,dx,dz,xsrc,zsrc
C
C---- check boundaries
      if (zsrc .eq. zmin) then
         if (xsrc .lt. xmax) then
            angend = 0.D0
            if (xsrc .eq. xmin) then
               angstr = -1.5707963267948966D0
            else
               angstr = -3.1415926535897931D0
            endif
         else
            angend = -1.5707963267948966D0
            angstr = -3.1415926535897931D0
         endif
      else if (zsrc .eq. zmax) then
         if (xsrc .lt. xmax) then
            angstr =0.D0
            if (xsrc .eq. xmin) then
               angend = 1.5707963267949001D0
            else
               angend = 3.1415926535897931D0
            endif
         else
            angstr = 1.5707963267948966D0
            angend = 3.1415926535897931D0
         endif
      else if (xsrc .eq. xmin) then
         angstr = -1.5707963267949001D0
         angend = 1.5707963267949001D0
      else if (xsrc .eq. xmax) then
         angstr = 1.5707963267948966D0
         angend = 4.7123889803846932D0
      else
C
C------- source within model sweep over 2 pi
         angstr = 0.0D0
         angend = 6.2831853071795862D0
c         print *,' Error: source is not located on the grid boundary.'
c         print *,' xmin ',xmin,' xmax ',xmax,' zmin ',zmin,' zmax ',zmax
c         read(*,*)
c         stop
      endif
      dxm = max(xmax - xsrc, xsrc - xmin)
      dzm = max(zmax - zsrc, zsrc - zmin)
      c2 = dxm*dxm + dzm*dzm
      c = sqrt(c2)
      bx = sqrt(c2 - dx*dx)
      bz = sqrt(c2 - dz*dx)
      ax = abs(acos((bx*bx+c2-dx*dx)/(2.D0*bx*c)))
      az = abs(acos((bz*bz+c2-dz*dz)/(2.D0*bz*c)))
C
C---- calculate preliminary delta theta
      dp = 0.5D0*min(ax,az)
C
C---- find delta theta that is <= dp that fits between angstr and angend
      bx = abs(angend-angstr)
      np = 1+ INT(bx/dp)
      dp = bx/dble(np-1)
C      print *,'np ',np
C
C---- make sure the start is exact
      p(1) = angstr
      sang(1) = sin(angstr)
      cang(1) = cos(angstr)
      do i=2,np-1
         p(i) = p(i-1) + dp
         sang(I) = sin(p(i))
         cang(i) = cos(p(i))
      enddo
C
C---- make sure the end is exact
      p(np) = angend
      sang(np) = sin(p(np))
      cang(np) = cos(p(np))
C
C---- find radial increment and number of rings
      dr = 0.5D0*min(dx,dz)
      c = c + 1.5D0*dr
      nr = 1 + int(c/dr)
C      print *,'nr ',nr
      dr = c/dble(nr-1)
      r(0) = 0.D0
      r(1) = dr
      do i=2,nr+1
         r(i) = r(i-1) + dr
      enddo
      end
c
c======================================================================
c
      subroutine crnswp(iostr,ioend,time,mxp,mxr,j,j1,s,
     &     idir,r24,sfac,dr2,rdr4,facnd,drrdpi,tm1cof,tm2cof,tm3cof,dd)
C
C     CRNSWP calculates times to a corner in an outward sweep over a radius
C     ring.
C
      IMPLICIT NONE
      INTEGER iostr,ioend,MXP,MXR,J,I,J1,I1,idir
      DOUBLE PRECISION TIME(MXP,MXR),r24,sfac,dr2,rdr4,facnd,
     &     drrdpi,tm1cof,tm2cof,tm3cof,tmcorn,s(1:*),tdiff,dd,tmp
      EXTERNAL tmcorn
C
C
C------- calculate times within segment
      DO I= iostr,ioend,idir
         I1 = I + idir
         tdiff = time(i1,j1) - time(i,j1)
C
C------- check illumination condition (from Podvin and Lecomte, 1991 JGI)
         if (tdiff .ge. 0.d0 .and. tdiff .le. s(i)*dd) then
            tmp = tmcorn(time(i,j1),time(i1,j1),time(i,j),s(i),r24,
     &           sfac,dr2,rdr4,facnd,drrdpi,tm1cof,tm2cof,tm3cof)
            TIME(I1,J) = min(tmp,time(i1,j))
         endif
      enddo
      return
      end
c
c======================================================================
c
      SUBROUTINE crntim(time1,time2,time3,time4,slocel,slo1,slo2,
     &     ifrstc,ifrstr,h,h2,HTSQR2,HS2I)
C-----------------------------------------------------------------
C     CRNTIM calculates travel times for a corner and returns the
C      minimum travel times in time4.
C
C         slo2
C     4 +      + 3
C
C slo1   slocel
C
C     2 +      + 1 times 1,2,3 are known and 4 is returned
C
C-----------------------------------------------------------------
      IMPLICIT NONE
      REAL time1,time2,time3,time4,slocel,slo1,slo2,h,h2,HTSQR2,HS2I
      INTEGER ifrstc,ifrstr,ihit
C
C---- 2d diffraction
      call diftim(TIME1,TIME4,slocel,HTSQR2,ihit)
C
C---- 2d transmission
      call trntim(TIME1,TIME2,TIME3,TIME4,slocel,h2,hs2i,ihit,0)
      call trntim(TIME1,TIME3,TIME2,TIME4,slocel,h2,hs2i,ihit,0)
C
C---- 2d refraction
      call reftim(time2,time4,slocel,slo1,H,ifrstr)
      call reftim(time3,time4,slocel,slo2,H,ifrstc)
      RETURN
      END
c
c======================================================================
c
      SUBROUTINE defcon(dr,dp,dr2,dr3,dp2,drdp,dp2dr2,dp2i,dridpi,
     &     drdp22)
      IMPLICIT NONE
      DOUBLE PRECISION dr,dp,dr2,dp2,dp2dr2,dp2i,dridpi,drdp22,dr3,
     &     drdp
C
c---- things calculated before starting grid calculations
c
      dr2 = dr*dr
      dr3 = dr2*dr
      dp2 = dp*dp
      dp2dr2 = -dp2/dr2
      dp2i = 1.0D0/dp2
      drdp = dr*dp
      dridpi = 1.0D0/drdp
      drdp22 = dr2*dp2
      RETURN
      END
c
c======================================================================
c
      SUBROUTINE difswp(iostr,ioend,TIME,MXP,MXR,J,J1,dd,s,idir)
C
C     DIFSWP calculates the diagonal diffractions across the cells
C     of a radius ring
C
      IMPLICIT NONE
      INTEGER iostr,ioend,J,idir,I,I1,J1,MXP,MXR
      DOUBLE PRECISION TIME(MXP,MXR),dd,s(1:*),tmp
C
C---- calculate times within segment
      DO I= iostr,ioend,idir
         I1 = I + IDIR
         tmp = dd*s(i) + TIME(I,J1)
C
C------- if time is earlier, save it and set flag
         TIME(I1,J) = min(tmp,TIME(I1,J))
      ENDDO
      RETURN
      END
c
c======================================================================
c
      subroutine DIFTIM(TIMEK,TIMEC,SLOW,HTSQR2,IFIRST)
C
C     DIFTIM calculates the time of a diffraction across a cell
C
      IMPLICIT NONE
      INTEGER IFIRST
      REAL TIMEK,TIMEC,SLOW,TMP,HTSQR2
      TMP = TIMEK + HTSQR2*SLOW
      IF (TMP .LT. TIMEC) THEN
         IFIRST = 1
         TIMEC = TMP
      ENDIF
      RETURN
      END
c
c======================================================================
c
      subroutine fndcmn(time,mxx,nz,nc,itop,ibot,nmin,minstr,mnlend,
     &     mnuend,tmmin,indexm)
      implicit none
      integer mxx,nz,itop,ibot,nmin,minstr(1:*),mnlend(1:*),
     &     mnuend(1:*),i,nc,lookmx,indexm(1:*)
      real time(mxx,nz),tdif,tdifl,sdif,sdifl,tmmin(1:*)
      do i = itop, ibot -1
         tdif = time(nc,i+1)-time(nc,i)
         sdif = sign(1.0,tdif)
         if(i .eq. itop) then
            nmin = 1
            if (tdif .ge. 0.0) then
C
C------------- we are starting at a minimum
               minstr(nmin) = i
               tmmin(nmin) = time(nc,i)
               mnlend(nmin) = i+1
C
C------------- flag other end no yet found
               mnuend(nmin) = -10
               lookmx = 1
            else
C
C------------- we are starting at a maximum
               mnlend(nmin) = i+1
               minstr(nmin) = -10
               mnuend(nmin) = -10
               lookmx = 0
            endif
         else
C
C---------- look for sign change indicating a min or max
            if (sdif .ne. sdifl) then
               if (sdif .gt. 0.0) then
C
C---------------- its a minimum
                  if (lookmx .ne. 0) print *,'trouble, looking max',
     &                 ' found a minimum instead'
                  minstr(nmin) = i
                  tmmin(nmin) = time(nc,i)
                  lookmx = 1
               else
C
C---------------- its a maximum
                  if (lookmx .eq. 0) print *,'trouble, looking min',
     &                 ' found a maximum instead'
                  mnuend(nmin) = i-1
                  nmin = nmin+1
                  mnlend(nmin) = i+1
                  minstr(nmin) = -10
                  mnuend(nmin) = -10
                  lookmx = 0
               endif
            endif
         endif
         tdifl = tdif
         sdifl = sdif
      enddo
C
C---- cleanup loose ends
      if (minstr(nmin) .eq. -10) then
         minstr(nmin) = ibot
         mnuend(nmin) = ibot-1
         tmmin(nmin) = time(nc,ibot)
      else if(mnuend(nmin) .eq. -10) then
         mnuend(nmin) = ibot-1
      endif
C
C---- sort by time and return sort index
      call INDEXX(nmin,tmmin,indexm)
      return
      end
c
c======================================================================
c
      subroutine fndrmn(time,mxx,nz,nc,ileft,iright,nmin,minstr,mnlend,
     &     mnuend,tmmin,indexm)
      implicit none
      integer mxx,nz,ileft,iright,nmin,minstr(1:*),mnlend(1:*),
     &     mnuend(1:*),i,nc,lookmx,indexm(1:*)
      real time(mxx,nz),tdif,tdifl,sdif,sdifl,tmmin(1:*)
      do i = ileft, iright -1
         tdif = time(i+1,nc)-time(i,nc)
         sdif = sign(1.0,tdif)
         if(i .eq. ileft) then
            nmin = 1
            if (tdif .ge. 0.0) then
C
C------------- we are starting at a minimum
               minstr(nmin) = i
               tmmin(nmin) = time(i,nc)
               mnlend(nmin) = i+1
C
C------------- flag other end no yet found
               mnuend(nmin) = -10
               lookmx = 1
            else
C
C------------- we are starting at a maximum
               mnlend(nmin) = i+1
               minstr(nmin) = -10
               mnuend(nmin) = -10
               lookmx = 0
            endif
         else
C
C---------- look for sign change indicating a min or max
            if (sdif .ne. sdifl) then
               if (sdif .gt. 0.0) then
C
C---------------- its a minimum
                  if (lookmx .ne. 0) print *,'trouble, looking max',
     &                 ' found a minimum instead'
                  minstr(nmin) = i
                  tmmin(nmin) = time(i,nc)
                  lookmx = 1
               else
C
C---------------- its a maximum
                  if (lookmx .eq. 0) print *,'trouble, looking min',
     &                 ' found a maximum instead'
                  mnuend(nmin) = i-1
                  nmin = nmin+1
                  mnlend(nmin) = i+1
                  minstr(nmin) = -10
                  mnuend(nmin) = -10
                  lookmx = 0
               endif
            endif
         endif
         tdifl = tdif
         sdifl = sdif
      enddo
C
C---- cleanup loose ends
      if (minstr(nmin) .eq. -10) then
         minstr(nmin) = iright
         mnuend(nmin) = iright-1
         tmmin(nmin) = time(iright,nc)
      else if(mnuend(nmin) .eq. -10) then
         mnuend(nmin) = iright-1
      endif
C
C---- sort by time and return sort index
      call INDEXX(nmin,tmmin,indexm)
      return
      end
c
c======================================================================
c
      subroutine forwrd(r,dr,dp,dr2,dr3,dp2,dp2dr2,dp2i,dridpi,
     &        drdp22,drdp,cosp,sring,time,mxp,mxr,j,j1,rftest,np)
C
C     forwrd calculates outward times for the current ring using
C     all the appropriate time stencils
C
      implicit none
      integer mxp,mxr,j,j1,rftest,np
      double precision r,dr,dp,dr2,dr3,dp2,dp2dr2,dp2i,dridpi,drdp22,
     &     sring(0:mxp,2),time(mxp,mxr),drpr,drpr2,drpr2i,rdr4,r24,
     &     drrdpi,tm1cof,tm2cof,tm3cof,facnd,sfac,dd,drdp,cosp
C
C---- pre-calculate radius dependent parameters
      call RDFCON(r,dr,dp,dr2,dp2,dp2dr2,dp2i,dridpi,drdp22,drpr,
     &     drpr2,drpr2i,rdr4,r24,drrdpi,tm1cof,tm2cof,tm3cof,facnd,
     &     sfac,dd,drdp,cosp)
C
C---- CALCULATE RADIAL TRAVEL-TIMES
      rftest=0
      CALL RADTIM(np,TIME,MXP,MXR,J,J1,dr,sring(0,1))
C
C---- CALCULATE DIFFRACTION TRAVEL-TIMES
C---- sweep to max angle
      call difswp(1,np-1,TIME,MXP,MXR,J,J1,dd,sring(1,1),1)
C---- sweep to min angle
      call difswp(np,2,TIME,MXP,MXR,J,J1,dd,sring(0,1),-1)
C
C---- calculate corner times
C---- sweep to max angle
      call crnswp(1,np-1,TIME,MXP,MXR,J,J1,sring(1,1),1,
     &     r24,sfac,dr2,rdr4,facnd,drrdpi,tm1cof,tm2cof,tm3cof,dd)
C---- sweep to min angle
      call crnswp(np,2,TIME,MXP,MXR,J,J1,sring(0,1),-1,
     &     r24,sfac,dr2,rdr4,facnd,drrdpi,tm1cof,tm2cof,tm3cof,dd)
C
C---- calculate refracted times
C---- sweep to maximum angle
      rftest = 0
      CALL REFSWP(1,np-1,TIME,MXP,MXR,J,drdp,sring(1,1),
     &     sring(1,2),1,rftest)
C---- sweep to minimum angle
      CALL REFSWP(np,2,TIME,MXP,MXR,J,drdp,sring(0,1),
     &     sring(0,2),-1,rftest)
      return
      end
c
c======================================================================
c
      SUBROUTINE INDEXX(N,ARRIN,INDX)
C
C     HEAPSORT indexing routine from numerical recipes
C
      implicit none
      integer n,indx(1:*),I,J,L,IR,INDXT
      real arrin(1:*),Q
      DO  J=1,N
        INDX(J)=J
      enddo
      if(n.lt.2) return
      L=N/2+1
      IR=N
      do while (1 .ne. 0)
         IF(L.GT.1)THEN
            L=L-1
            INDXT=INDX(L)
            Q=ARRIN(INDXT)
         ELSE
            INDXT=INDX(IR)
            Q=ARRIN(INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
               INDX(1)=INDXT
               RETURN
            ENDIF
         ENDIF
         I=L
         J=L+L
         do while (J.LE.IR)
            IF(J.LT.IR)THEN
               IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
            ENDIF
            IF(Q.LT.ARRIN(INDX(J)))THEN
               INDX(I)=INDX(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         enddo
         INDX(I)=INDXT
      enddo
      END
c
c======================================================================
c
      subroutine intxzt(dr,dp,r,p,time,mxp,np,nr,xsrc,zsrc,timexz,mxx,
     &     nx,nz,dx,dz,xmin,xmax,zmin,zmax,iangr)
C
C     inttim calculates times for requested x,z receiver positions
C     using bilinear interpolation of the polar travel-time field
C
      implicit none
      integer mxp,nr,i,ir,ip,ip1,ir1,mxx,nz,np,j,ixs,ixe,
     &     izs,ize,nx,iangr
      double precision dr,dp,p(1:*),time(mxp,nr),rreq,preq,fp,fr,
     &     za,zb,x,z,r(0:*),onemfp,t1,t2
      real xsrc,zsrc,timexz(mxx,nz),dx,dz,xmin,xmax,zmin,zmax,xrec,
     &     zrec
C
C---- compute portion of cartesian grid that needs times from the polar grid
      ixs = max(1,int(xmin/dx))
      ixe = min(nx,1+nint(xmax/dx))
      izs = max(1,int(zmin/dz))
      ize = min(nz,1+nint(zmax/dz))
C      print *,'ixs,ixe,izs,ize'
C      print *,ixs,ixe,izs,ize
C      print *,'nz,zmax,dz,nint(zmax/dz)'
C      print *,nz,zmax,dz,nint(zmax/dz)
C
C---- loop through cartesian subgrid
      do i = izs, ize
         zrec = float(i-1)*dz
         do j = ixs,ixe
            xrec = float(j-1)*dx
            x = xrec - xsrc
C
C---------- make z negative down
            z = zsrc - zrec
C
C---------- calculate polar coordinates
            rreq = sqrt(x*x + z*z)
            if (z .ne. 0.0 .or. x .ne. 0.0) then
               preq = atan2(z,x)
            else
               preq = p(1)
            endif
C
C---------- have to map the angle depending on the signs of the angle range
C---------- or if source position inside model requires 360 degree sweep
            if (p(1) .gt. 0.D0 .or. iangr .eq. 0) then
               if (preq .lt. 0.D0) preq = 6.2831853071795862D0+preq
            else if(p(1) .lt. -3.D0 .and.
     &              preq .eq. 3.1415926535897931D0) then
               preq = -3.1415926535897931D0
            endif
C
C---------- define angle within range if source is at polar origin
            if (rreq .eq. 0.0) preq = p(1)
            if (rreq .ge. 0.0 .and. rreq .le. r(nr) .and.
     &           preq .ge. p(1) .and. preq .le. p(np)) then
C
C------------- calculate index position on polar grid
               ir = int(rreq/dr)
               ip = 1 + int((preq-p(1))/dp)
               if (ip .eq. np) ip = ip - 1
C
C------------- do bilinear interpolation
               fp = (preq-p(ip))/dp
               fr = (rreq-r(ir))/dr
               ip1 = ip+1
               ir1 = ir+1
               onemfp = 1.D0-fp
C
C------------- if time needed is not source point compute interpolation
               if (ir .ne. 0) then
                  za = fp*time(ip1,ir) + time(ip,ir)*onemfp
               else
C
C---------------- polar source time is zero
                  za = 0.d0
               endif
               zb = fp*time(ip1,ir1) + time(ip,ir1)*onemfp
               timexz(j,i) = fr*zb + za*(1.D0-fr)
            endif
         enddo
      enddo
      return
      end
c
c======================================================================
c
      SUBROUTINE radtim(np,TIME,MXP,MXR,J,J1,dr,s)
C
C     RADTIM calculates the simple radial times from one radius ring
C     to the next
C
      IMPLICIT NONE
      INTEGER np,MXP,MXR,J,I,J1
      DOUBLE PRECISION TIME(MXP,MXR),dr,s(0:*),tmp

C
C---- calculate times within ring
      DO I=1,np
         tmp = dr*min(s(i-1),s(i)) + TIME(I,J1)
         TIME(I,J) = min(tmp,TIME(i,j))
      enddo
      RETURN
      END
c
c======================================================================
c
      subroutine raypth(time,slow,mxx,nx,nz,h,xis,zis,xir,zir,xpos,zpos,
     &     dkern,mxa,nrv,nr)
      implicit none
      integer mxx,mxa,nx,nz,nr,nrv(1:*),i,j,isx,isz,nk,it,ix,ix1,
     &     ixm1,iz,iz1,izm1,ixn,ixl,izn,izl,nx1,nz1,k,m,xpos(mxa,nr),
     &     zpos(mxa,nr),nk1
      real time(mxx,nz),slow(mxx,nz),h,xis,zis,xir(1:*),zir(1:*),
     &     dkern(mxa,nr),HI,xi,zi,rinc,
     &     lwzd,upzd1,lwzd1,f,dtdz,dtdx,lfxd,lfxd1,rtxd,rtxd1,xn,zn,fnx,
     &     fnz,xinc,ang,cang,sang,xdif,zdif,angl,angc,xc,xs,zc,zs,dif,
     &     upzd,zinc,xisw,zisw
C
C---- setup source cell index (there is one less cell than time in each dim.)
      hi = 1.0/H
      fnx = float(nx)
      fnz = float(nz)
      nx1 = nx-1
      nz1 = nz-1
      xisw = xis*hi
      zisw = zis*hi
      isx = max(1,min(1+int(xisw),nx1))
      isz = max(1,min(1+int(zisw),nz1))
      xisw = xisw+1.
      zisw = zisw+1.
C
C---- loop through number of receiver positions
      do i=1,nr
         nk=1
         rinc = 0.
         it = 1
         xi = 1.+xir(i)*hi
         zi = 1.+zir(i)*hi
         xpos(1,1) = 0.
         zpos(1,1) = 0.
C
C------- calculate path from receiver to source position until the path
C------- intersects the source cell
         do while (it .ne. 0)
C
C---------- find cell indices about current path position
            ix = min(int(xi),nx)
            ix1 = min(ix+1,nx)
            ixm1 = max(1,ix-1)
            iz = min(int(zi),nz)
            iz1 = min(iz+1,nz)
            izm1 = max(1,iz-1)
C
C---------- compute quantities for z component of the gradient
            upzd = time(ix,izm1)-time(ix,iz)
            lwzd = time(ix,iz)-time(ix,iz1)
            upzd1 = time(ix1,izm1)-time(ix1,iz)
            lwzd1 = time(ix1,iz)-time(ix1,iz1)
C
C---------- handle grid edges gracefully
            if (izm1 .eq. iz) then
               upzd = sign(1.e-6,lwzd)
               upzd1 = sign(1.e-6,lwzd1)
            else if(iz .eq. iz1) then
               lwzd = sign(1.e-6,upzd)
               lwzd1 = sign(1.e-6,upzd1)
            endif
C
C---------- linear interpolation relative to x position
            f = xi - float(ix)
            upzd = f * upzd1 + (1.0-f)* upzd
            lwzd = f * lwzd1 + (1.0-f)* lwzd
C
C---------- calculate z time gradient, checking for time discontinuity
            if (sign(1.0,upzd) .ne. sign(1.0,lwzd)) then
               dtdz = 0.
            else
               dtdz = 0.5*(upzd+lwzd)
            endif
C
C---------- x time gradient
            lfxd = time(ixm1,iz)-time(ix,iz)
            lfxd1 = time(ixm1,iz1)-time(ix,iz1)
            rtxd = time(ix,iz)-time(ix1,iz)
            rtxd1 = time(ix,iz1)-time(ix1,iz1)
            if (ixm1 .eq. ix) then
               lfxd = sign(1.e-6,rtxd)
               lfxd1 = sign(1.e-6,rtxd1)
            else if(ix .eq. ix1) then
               rtxd = sign(1.e-6,lfxd)
               rtxd1 = sign(1.e-6,lfxd1)
            endif
C
C---------- linear interpolation
            f = zi - float(iz)
            lfxd = f * lfxd1 + (1.0-f) * lfxd
            rtxd = f * rtxd1 + (1.0-f) * rtxd
C
C---------- calculate x time gradient, checking for time discontinuity
            if (sign(1.0,lfxd) .ne. sign(1.0,rtxd)) then
               dtdx = 0.
            else
               dtdx = 0.5*(lfxd+rtxd)
            endif
C
C---------- calculate gradient direction
            ang=atan2(dtdz,dtdx)
            cang = cos(ang)
            sang = sin(ang)
C
C---------- calculate position of projected point
            xn = max(1.,min(xi+cang,fnx))
            zn = max(1.,min(zi+sang,fnz))
C
C---------- calculate index location values
            ixn = int(xn)
            ixl = int(xi)
            izn = int(zn)
            izl = int(zi)
C
C---------- if we have left the previous cell calculate its path length
            if (ixn .ne. ixl .or. izn .ne. izl) then
               zdif = zn-zi
               xdif = xn-xi
               if (xi .eq. float(ixl) .and. xdif .lt. 0.0) ixl = ixl-1
               if (zi .eq. float(izl) .and. zdif .lt. 0.0) izl = izl-1
C
C------------- find intersection with the cell wall by looking at angles of
C------------- the ray and the angle to the corner of the sub-cell the ray
C------------- departs through
               angl = atan2(abs(zdif),abs(xdif))
               if (cang .lt. 0.0) then
                  xs = abs(xi - float(ixl))
               else
                  xs = abs(1.+float(ixl)-xi)
               endif
               if (sang .lt. 0.0) then
                  zs = abs(zi - float(izl))
               else
                  zs = abs(1.+float(izl)-zi)
               endif
               angc = atan2(zs,xs)
               if (angl .gt. angc) then
                  zc = zs
                  if (xdif .ne. 0.0) then
                     xc = zc/tan(angl)
                  else
                     xc = 0.
                  endif
               else
                  xc = xs
                  zc = xc*tan(angl)
               endif
C
C------------- check for vertical refractions and use index of fast side
               j = max(1,min(ixl,nx1))
               k = min(ixl,nx1)
               m = min(izl,nz1)
               if (xn .ne. xi) then
                  xpos(nk,i) = j
               else
                  if (slow(k,m) .lt. slow(max(1,ixl-1),m)) then
                     xpos(nk,i) = j
                  else
                     xpos(nk,i) = max(1,ixl-1)
                  endif
               endif
C
C------------- check for horizontal refractions and use index of fast side
               j = max(1,min(izl,nz1))
               k = min(izl,nz1)
               m = min(ixl,nx1)
               if (zn .ne. zi) then
                  zpos(nk,i) = j
               else
                  if (slow(m,k) .lt. slow(m,max(1,izl-1))) then
                     zpos(nk,i) = j
                  else
                     zpos(nk,i) = max(1,izl-1)
                  endif
               endif
               dkern(nk,i) = rinc + h*sqrt(xc*xc+zc*zc)
C
C------------- check for and eliminate duplicates
               nk1 = nk -1
               if (nk1 .ne. 0 .and. xpos(nk1,i) .eq. xpos(nk,i) .and.
     &              zpos(nk1,i) .eq. zpos(nk,i)) then
                  dkern(nk1,i) = dkern(nk1,i) + dkern(nk,i)
                  dkern(nk,i) = 0.0
               endif
               rinc = 0.0
C
C------------- require nonzero contribution
               if (dkern(nk,i) .gt. 0.0) nk = nk+1
C
C------------- if both x and z indices change we probably hit the corner of
C------------- second cell and need to calculate its path length and store it
C------------- before proceeding
               if (ixn .ne. ixl .and. izn .ne. izl .and.
     &              angl .ne. angc) then
C
C---------------- move to boundary and see which of x or z advances
                  if (angl .gt. angc) then
                     xc = abs(xs-xc)
                     zc = xc*tan(angl)
                     if (sang .gt. 0.0) then
                        izl = min(izl+1,nz1)
                     else
                        izl = max(1,izl-1)
                     endif
                     xs = xi + xs*sign(1.0,xdif)
                     zs = zi + (zs + zc)*sign(1.0,zdif)
                  else
                     zc = abs(zs-zc)
                     xc = zc/tan(angl)
                     if (cang .gt. 0.0) then
                        ixl = min(ixl+1,nx1)
                     else
                        ixl = max(1,ixl-1)
                     endif
                     zs = zi + zs*sign(1.0,zdif)
                     xs = xi + (xs + xc)*sign(1.0,xdif)
                  endif
C
C---------------- don't need to check for refraction because the propagation
C---------------- direction must be oblique to the coordinate axes by
C---------------- definition to be within this if section
                  xpos(nk,i) = max(1,min(ixl,nx1))
                  zpos(nk,i) = max(1,min(izl,nz1))
                  dkern(nk,i) = h*sqrt(xc*xc+zc*zc)
C
C---------------- check for and eliminate duplicates
                  nk1 = nk -1
                  if (nk1 .ne. 0 .and. xpos(nk1,i) .eq. xpos(nk,i) .and.
     &                 zpos(nk1,i) .eq. zpos(nk,i)) then
                     dkern(nk1,i) = dkern(nk1,i) + dkern(nk,i)
                     dkern(nk,i) = 0.0
                  endif
                  if (dkern(nk,i) .gt. 0.0) nk = nk+1
               else
C
C---------------- update tail position
                  xs = xi + xc*sign(1.0,xdif)
                  zs = zi + zc*sign(1.0,zdif)
               endif
C
C------------- store leftover length projecting into the next (unfinished) cell
               xinc = abs(xn-xs)
               zinc = abs(zn-zs)
               rinc = rinc + h*sqrt(xinc*xinc+zinc*zinc)
            else
C
C------------- still have not left current cell
               xinc = abs(xn-xi)
               zinc = abs(zn-zi)
               rinc = rinc + h*sqrt(xinc*xinc+zinc*zinc)
            endif
C
C---------- check if we have hit source cell
            xs = xn-xisw
            zs = zn-zisw
            dif = sqrt(xs*xs+zs*zs)
C
C---------- if we are within h of the source position exit the receiver loop
            if (dif .le. 1.0 .or. (xi .eq. xn .and. zi .eq. zn)) it = 0
C            if (dif .le. 1.0) it = 0
C
C---------- next last position becomes current new position
            xi = xn
            zi = zn
         enddo
C
C------- pickup final path incrment to source
C         print *,' number of path segments =',nk,i,xis,zis
         nk = max(1,nk-1)
C
C------- if source cell does not have a contribution yet add it to the list
         if (xpos(nk,i) .ne. isx .or. zpos(nk,i) .ne. isz) then
            nk = nk + 1
            xpos(nk,i) = max(1,min(isx,nx1))
            zpos(nk,i) = max(1,min(isz,nz1))
            dkern(nk,i) = rinc + h*dif
         else
C
C---------- the final segment is added to the source cell total
            dkern(nk,i) = dkern(nk,i)+h*dif+rinc
         endif
         nrv(i) = nk
      enddo
      return
      end
c
c======================================================================
c
      SUBROUTINE rdfcon(r,dr,dp,dr2,dp2,dp2dr2,dp2i,dridpi,drdp22,drpr,
     &     drpr2,drpr2i,rdr4,r24,drrdpi,tm1cof,tm2cof,tm3cof,facnd,
     &     sfac,dd,drdp,cosp)
      IMPLICIT NONE
      DOUBLE PRECISION r,dr,dr2,dp2,dp2dr2,dp2i,dridpi,drdp22,drpr,
     &     drpr2,drpr2i,rdr,r24,ridrpr,drrdpi,tm1cof,tm2cof,tm3cof,
     &     facnd,sfac,rdr4,r2,dd,drdp,cosp,dp
C
c---- things calculated at each radius ring
c
      drpr = dr + r
      drdp = drpr*dp
      drpr2 = drpr*drpr
      drpr2i = -1.0D0/drpr2
      rdr = dr*r
      rdr4 = 4.0D0*rdr
      r2 = r*r
      r24 = 4.0D0*r2
      ridrpr = 1.0D0/(r*drpr)
      drrdpi = -ridrpr*dridpi
      tm1cof = dp2i*(dp2dr2-ridrpr)
      tm2cof = dp2i*(ridrpr+dp2dr2)
      tm3cof = dp2i*(drpr2i-dp2dr2)
      facnd = dp2/(drpr2i+dp2dr2)
      sfac = drdp22 + 2.0D0*dr*dp2*r + dr2 + r2*dp2
C
C---- calculate diffraction distance across cell
      dd = sqrt(dr2 + drpr2 - 2.0D0*drdp*r*cosp)
      return
      end
c
c======================================================================
c
      subroutine refswp(iostr,ioend,time,mxp,mxr,nrn,
     &     drdp,s1,s2,idir,rftest)
      implicit none
C
C     REFSWP calculates the edge refraction at the cell boundary along
C     a sweep of a radius ring
C
      integer mxp,mxr,nrn,i,i1,idir,rftest,j,j1,iostr,ioend
      DOUBLE PRECISION time(mxp,mxr),drdp,s1(1:*),s2(1:*),tmp
C
C---- loop through a ring segment
      DO I = iostr, ioend,idir
         I1 = I+idir
         tmp = drdp*min(s1(i),s2(i))+TIME(I,nrn)
C
C------- check for lower time
         if (tmp .lt. TIME(I1,nrn)) then
            rftest = 1
            TIME(I1,nrn) = tmp
C
C---------- time the rest of the ring (don't need the rest of the ring
C---------- segments since we hit an absorbing boundary first)
            do j = i1, ioend, idir
               j1 = j + idir
               time(j1,nrn) = min(drdp*min(s1(j),s2(j))+TIME(j,nrn),
     &              TIME(j1,nrn))
            enddo
         ENDIF
      enddo
      return
      end
c
c======================================================================
c
      SUBROUTINE REFTIM(TIMEK,TIMEC,SLOW1,SLOW2,H,IFIRST)
C
C     REFTIM calculates the refracted time along cell boundaries
C
      IMPLICIT NONE
      INTEGER IFIRST
      REAL TIMEK,TIMEC,SLOW1,SLOW2,H,TMP
      tmp = TIMEK+MIN(SLOW1,SLOW2)*H
      if (tmp .lt. TIMEC) then
         if (slow1 .lt. slow2) ifirst = 1
         TIMEC = tmp
      endif
      RETURN
      END
c
c======================================================================
c
      subroutine setrng(xmin,xmax,zmin,zmax,xsrc,zsrc,bigslo,angstr,
     &     angend,dx,dz,dr,dp,nr,np,r,p,cang,sang,sring,mxp,sloxz,
     &     mxx,nx,nz,ir,init,iangr)
C
C     setrng determines the indices of the angle limits needed for each
C     ring to cover its portion of the rectangular grid.
C     (for now, however, just use the entire grid defined by nr,np,dr,dp
C      and fill the parts outside the grid with large slownesses)
C
      implicit none
      integer nr,np,mxp,mxx,nx,nz,i,j,nxp1,nzp1,ix,iz,k,j1,
     &     ir,init,inc,m,i1,i2,iangr,nxm1, nzm1
      double precision angstr,angend,dr,dp,r(0:*),p(1:*),cang(1:*),
     &     sang(1:*),sring(0:mxp,2),bigslo,dxi,dzi,sum,sumlst,dzin
      real xmin,xmax,zmin,zmax,xsrc,zsrc,dx,dz,sloxz(mxx,nz)
      dxi = 1.0D0/dble(dx)
      dzi = 1.0D0/dble(dz)
      dzin = -dzi
      nxp1 = nx+1
      nzp1 = nz+1
      nxm1 = nx-1
      nzm1 = nz-1
c      open(11,file='slow.out',status='unknown',form='unformatted')
c      open(11,file='slow.out',status='unknown',form='unformatted')
c      write(11) nr,np,angstr,angend,dr,dp,dxi,dzi
c      write(11) (p(i),i=1,np)
c      write(11) (r(i),i=1,nr)
C
C---- if called for the first time calculate two rings
      if (init .ne. 0) then
         inc = 1
      else
C
C------- Use the result of the previous call to only calculate the new
C------- outer ring by moving the old one inward
         inc = 2
         do i= 0,np
            sring(i,1) = sring(i,2)
         enddo
      endif
C
C---- calculate needed slowness ring(s)
      do i=inc,2
         i1 = i - 1
         i2 = i - 2
C
C------- pre-calculate the first radial edge
         sumlst = 0.D0
         do m=i2,i1
            k = ir+m
C
C---------- calculate indices of x,z slowness cell the polar grid
C---------- point lands on
            ix = int((r(k)*cang(1) + xsrc)*dxi)+1
            iz = int(r(k)*sang(1)*dzin + zsrc*dzi)+1
C
C---------- check if requested cell is within grid or within
C---------- a cell width of the grid (let the grid slowness out a bit)
            if (ix .gt. nxp1 .or. iz .gt. nzp1) then
C
C------------- make points outside, very slow
               sumlst = sumlst + bigslo
            else
               sumlst = sumlst + sloxz(max(1,min(ix,nxm1)),
     &              max(1,min(iz,nzm1)))
            endif
         enddo
C
C------- loop through angles
         do j=0,np
C
C---------- make cells next to edges slow if not a 360 degree sweep
            if (j .eq. 0 .or. j .eq. np .and. iangr .ne. 0) then
               sring(j,i) = bigslo
            else
C
C------------- calculate slowness for polar cell
               j1 = j + 1
C
C------------- wrap angle around for 360 degree sweep
               if (j1 .gt. np) j1 = 2
C
C------------- convert four corners to x and z
               sum = 0.D0
C
C------------- loop through 2 next corners of polar cell
               do m=i2,i1
                  k = ir+m
C
C---------------- calculate indices of x,z slowness cell the polar grid
C---------------- point lands on
                  ix = int((r(k)*cang(j1) + xsrc)*dxi)+1
                  iz = int(r(k)*sang(j1)*dzin + zsrc*dzi)+1
C
C---------------- check if requested cell is within grid or within
C---------------- a cell width of the grid (let the grid slowness out a bit)
                  if (ix .gt. nxp1 .or. iz .gt. nzp1) then
C
C------------------- make points outside, very slow
                     sum = sum + bigslo
                  else
                     sum = sum + sloxz(max(1,min(ix,nxm1)),
     &                                 max(1,min(iz,nzm1)))
                  endif
               enddo
C
C------------- calculate polar cell slowness from the average of the 4 corners
               sring(j,i) = 0.25D0*(sum+sumlst)
C
C------------- store last radial edge value
               sumlst = sum
            endif
         enddo
c         write(11) (sring(j,i),j=0,np)
      enddo
      return
      end
c
c======================================================================
c
      DOUBLE PRECISION FUNCTION tmcorn(tm1,tm2,tm3,s,r24,sfac,dr2,
     &     rdr4,facnd,drrdpi,tm1cof,tm2cof,tm3cof)
      IMPLICIT NONE
      DOUBLE PRECISION tm1,tm2,tm3,s,r24,sfac,dr2,rdr4,facnd,drrdpi,
     &     tm1cof,tm2cof,tm3cof,tm2s,tm3tm2,tm1tm2,tsqr,tsqra
C
c     Solution to (dt/dr)^2 + (dt/dp)^2 = s*s for the time of the far corner
c     of finite difference cell in r and theta where,
c
c     dt/dr = (tm4+tm3 - tm1 -tm2)/(2dr)
c     dt/dp = (tm2-tm1)/(2rdp) + (tm4-tm3)/(2(r+dr)dp)
c
C---- horner expansion solution
c
c---- things calculated at each cell in the routine
      tmcorn = 1.D10
      tm2s = tm2*tm2
      tm3tm2 = tm3*tm2
      tm1tm2 = tm1*tm2
      tsqra = r24*(s*s*sfac - tm3*tm3 + 2.0D0*tm3tm2 - tm2s) +
     &     dr2*(2.0D0*tm1tm2 - tm1*tm1 - tm2s) +
     &     rdr4*(tm3tm2 + tm1tm2 - tm3*tm1 - tm2s)
c      if (tsqra .lt. 0.D0) return
      tsqr = sqrt(abs(tsqra))
      tmcorn = facnd*(tsqr*drrdpi+tm2cof*tm2+tm1cof*tm1+tm3cof*tm3)
      return
      end
c
c======================================================================
c
      subroutine TRNTIM(TIMEM,TIMEN,TIMEB,TIMEC,SLOW,H2,HS2I,IFIRST,
     &     IREV)
C
C     TRNTIM calculates the conditional transmission across a cell
C
      IMPLICIT NONE
      INTEGER IFIRST,IREV
      REAL TIMEM,TIMEN,TIMEB,TIMEC,SLOW,H2,HS2I,TMP,DIFF,DIFF2
      DIFF = TIMEN - TIMEM
      DIFF2 = TIMEB-TIMEN
C
C---- Podvin and Lecomte ilumination condition (with sign test added as 3rd)
      IF (DIFF .LT. 0.0 .OR. DIFF .GT. HS2I*SLOW .OR.
     &     TIMEB-TIMEM .LE. 0.0) RETURN
c      IF (DIFF .LT. 0.0 .OR. DIFF .GT. HS2I*SLOW) RETURN
c      IF (IREV .NE. 0) THEN
C
C------- Podvin and Lecomte stencil
c         TMP = TIMEN + SQRT(H2*SLOW*SLOW - DIFF*DIFF)
c      ELSE
C
C------- Vidale plane-wave stencil (more accurate going forward)
      TMP = TIMEM + SQRT(2.D0*H2*SLOW*SLOW - DIFF2*DIFF2)
c      ENDIF
      IF (TMP .LT. TIMEC) THEN
         IFIRST = 1
         TIMEC = TMP
      ENDIF
      RETURN
      END
