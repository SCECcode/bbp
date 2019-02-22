        parameter (MAXnw=5000, MAXnl=5000) 
	dimension slip(MAXnw,MAXnl)
        dimension ty(2000)
	character*120 slipfile,planar_fault,fltpar,dirin,filein
        character*1 dum1
        character*4, stat
        pi=3.1415926
        pi2=pi/2.

        read(5,'(a100)')dirin
	read(5,'(a100)')planar_fault
        read(5,'(a100)')fltpar
        read(5,*)isegm
        read(5,*)dx,dy

        i1=index(dirin,' ')-1
        i2=index(planar_fault,' ')-1
        filein=dirin(1:i1)//'/'//planar_fault(1:i2)
        open(9,file=fltpar)

        read(9,'(a1)')dum1
        read(9,*)ii
        do k=1,isegm
        read(9,*)stat,alon,alat,zz,xlen1,xlen2,ylen,x0,y0,
     +           stra,dipa,rakea
        enddo
        close(9)

        call convert(alon,alat,xsc,ysc)

        ni=int(xlen1/dx)
        ne=int(xlen2/dx)
        xlen=(ne-ni)*dx
         print*,'Fault Length',xlen
         print*,'Fault Width',ylen

        nx=int(xlen/dx)
        ny=int(ylen/dy)
        print*,'nx,ny = ',nx,ny

c        nx0=int(x0/dx)
c        ny0=int(y0/dy)

        area=dx*dy*1.0e+10
c    fault element area is in cm*cm

        open(34,file=filein)
        read (34,*)dum
        read (34,'(a1)') dum1
        read (34,*) dum,dum,nxf,nyf,dum,dum
        read (34,*) dum,dum,dum,dum,dum
        read (34,'(a1)')dum1

        slipfile=dirin(1:i1)//'/'//stat//'.'//planar_fault(1:i2)
        open(12,file=slipfile)
        version=1.0
        nums=1
        write (12,'(f3.1)')version
        write (12,'("PLANE ",i2)')nums
        write (12,815) alon,alat,nx,ny,xlen,ylen
        write (12,820) stra, dipa, zz, x0,y0
815     format (2f10.4,2i7,2f8.2)
820     format (5f8.2)

        write (12,'(a,i8)') "POINTS",nx*ny

       dl=(nx)*dx
       dw=(ny)*dy
       nx2=nx/2.
       dx2=0.5*dx
       dy2=0.5*dy
       strike=stra*pi/180.
       dip=dipa*pi/180.
       rake=rakea*pi/180.

        nxi=ni
        nxe=ne

c        if (nxe.eq.0) nxe=1
        k=0
        j=0
        do j=1,nyf
           i=0
           k=0
           do ii=1,nxf
        read (34,*)xxf,yyf,zf,strf,dipf,areaf,stimf,dtf,
     +    vrf,densf
        read (34,*) rakef,slip1f,npts1f,slip2f,n0,slip3f,n0
        
        if(npts1f.ne.0)read(34,400) (ty(it),it=1,npts1f)

         y=(j-1)*dy*cos(dip)+dy2
         z=zz+(j-1)*dy*sin(dip)+dy2
           if(ii.gt.nxi.and.ii.le.nxe) then
           i=i+1
            k=k+1
         x=(i)*dx+dy2
         str2=strike
         xx=alon+(x-dl/2)*sin(strike)/xsc+y*cos(str2)/xsc
         yy=alat+(x-dl/2)*cos(strike)/ysc-y*sin(str2)/ysc

        write (12,200)xx,yy,z,stra,dipa,areaf,stimf,dtf,
     +    vrf,densf
        write (12,*) rakef,slip1f,npts1f,slip2f,n0,slip3f,n0
        if(npts1f.ne.0) write (12,400) (ty(it),it=1,npts1f)
        endif
          enddo
        enddo
        if(nx.ne.k) then
        print*,'PROBLEM nx not eq. to k',nx,k
        stop
        endif
200     format(2(1x,f12.6),f9.2,2f8.2,1x,e13.5,1x,f8.4,3(1x,e12.4))
202     format(f10.4,3(1x,f10.2,i4))
300     format(f8.2,3(e13.5,i8))
400     format(6e13.5)
        end

	subroutine latlon2km (arg,latkm,lonkm,rc,g2)
	double precision arg
	real latkm,lonkm,rc,g2
	real cosA,sinA,g2s2,den
	data FONE, FTWO / 1.0, 2.0 /

	cosA = cos(arg)
	sinA = sin(arg)
	g2s2 = g2*sinA*sinA

	den = sqrt(1.0/(1.0+g2s2))
	lonkm = rc*cosA*den
c	latkm = rc * (sqrt(1.0 + g2s2*(2.0+g2))) *den*den*den
	latkm = (rc)*(sqrt((FONE) + g2s2*((FTWO) + (g2))))*den*den*den
	return
	end


         subroutine convert(alon,alat,ddlon,ddlat)
         slon=alon+1
         slat=alat
        i=0
        call DELAZ5( alat, alon, slat, slon, DELT, DELTDG,
     +ddlon, AZES, AZESDG, AZSE, AZSEDG, I )
        slat=slat+1
        slon=alon
        i=0
        call DELAZ5( alat, alon, slat, slon, DELT, DELTDG,
     +ddlat, AZES, AZESDG, AZSE, AZSEDG, I )
1        return
        end

      SUBROUTINE  DELAZ5( THEI, ALEI, THSI, ALSI, DELT, DELTDG,
     2DELTKM, AZES, AZESDG, AZSE, AZSEDG, I )
       DOUBLE  PRECISION C, AK, D, E, CP, AKP, DP, EP,
     2A, B, G, H, AP, BP, GP, HP
      IF(I) 50, 50, 51
C     IF  COORDINATES ARE GEOGRAPH DEG I=0
C     IF COORDINATES ARE GEOCENT RADIAN  I=1
   50 THE=1.745329252E-2*THEI
      ALE=1.745329252E-2*ALEI
      THS=1.745329252E-2*THSI
      ALS=1.745329252E-2*ALSI
      AAA=0.9931177*TAN(THE)
      THE=ATAN(AAA)
      AAA=0.9931177*TAN(THS)
      THS=ATAN(AAA)
      GO TO 32
   51 THE=THEI
      ALE=ALEI
      THS=THSI
      ALS = ALSI
   32 CONTINUE
      C= SIN(THE)
      AK=-COS(THE)
      D=SIN(ALE)
      E= -COS(ALE)
      A= AK*E
      B= -AK*D
      G=-C*E
      H=C*D
      CP=SIN(THS)
      AKP=-COS(THS)
      DP=SIN(ALS)
      EP = -COS(ALS)
      AP = AKP*EP
      BP=-AKP*DP
      GP=-CP*EP
      HP=CP*DP
      C1=A*AP+B*BP+C*CP
      IF( C1-0.94 )  30, 31, 31
   30 IF(C1+0.94) 28, 28, 29
   29 DELT=ACOS(C1)
   33 DELTKM=6371.0*DELT
      C3 = (AP-D)**2+(BP-E)**2+CP**2-2.0
      C4 = (AP-G)**2+(BP-H)**2+(CP-AK)**2-2.0
      C5 = (A-DP)**2+ (B-EP)**2+C**2-2.0
      C6 = (A-GP)**2+(B-HP)**2+(C-AKP)**2-2.0
      DELTDG = 57.29577951*DELT
      AZES = ATAN2(C3, C4 )
      IF ( AZES ) 80, 81, 81
   80 AZES = 6.283185308+ AZES
   81 AZSE = ATAN2( C5, C6 )
      IF ( AZSE ) 70, 71 , 71
   70 AZSE=6.283185308+AZSE
   71 AZESDG=57.29577951*AZES
      AZSEDG=57.29577951*AZSE
      RETURN
   31 C1=(A-AP)**2+(B-BP)**2+(C-CP)**2
      C1= SQRT(C1)
      C1=C1/2.0
      DELT = ASIN(C1)
      DELT= 2.0*DELT
      GO TO 33
   28 C1=(A+AP)**2+(B+BP)**2+(C+CP)**2
      C1 = SQRT(C1 )
      C1= C1/2.0
      DELT = ACOS(C1)
      DELT = 2.0*DELT
      GO TO 33
      END
 
