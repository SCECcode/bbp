c       GEN_SRF4.9
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c Copyright (c) <2016>, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory

c Written by Arben Pitarka, pitarka1@llnl.gov

c LLNL-CODE-687219.
c All rights reserved.

c Redistribution and use in source and binary forms, with or without modification, are
c permitted provided that the following conditions are met:
c • Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
c • Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the documentation 
c   uand/or other materials provided with the distribution.
c • Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
c THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, 
c LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (
c INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
c AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING

c IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c Additional BSD Notice
c 1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore National Laboratory 
c    under Contract No. DE-AC52-07NA27344 with the DOE.
c 2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or implied, 
c    or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, 
c    or represents that its use would not infringe privately-owned rights.
c 3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not necessarily 
c    constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, LLC. The views and 
c     opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence Livermore National Security, LLC, 
c     and shall not be used for advertising or product endorsement purposes.
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c  06/2016  Modified for random asperity location 
c  07/2016  Include three asperities for fault length > 70km
c  07/2016  Improve random asperity location
c  07/2016  Modify parameterization for the three asperity model
c  08/2016  add condition for single asperity fault
c  08/2016  correct bug for equivalent asperity radius when nasp=3
c  08/2016  adjust min distance from top of the fault at 2km
c  08/13/2016 adjust condition for single asperity limited depth  (accomodate 2 apsreities for Northridge)
c  09/01/2016 replace Miyatake's source time function 
c  10/25/2018 corrected fault width condition for single asperity 
c  02/02/2019 make the Vr/Vs ratio an input rupture model parameter
c  03/25/2019 the code outputs the stress drop . output file name: stress_drop.out
c  06/30/2020 fix the min asperity depth measured from the ground surface (instead of top fault edge) 
  
        parameter (MAXnw=6000, MAXnl=6000) 
	dimension slip(MAXnw,MAXnl),stres(MAXnw,MAXnl),alng(MAXnw,MAXnl)
     + ,ty(20000)
        dimension xasp(10),yasp(10),zvel(100),vel(100),den(100)
        dimension beta(MAXnl),rden(MAXnl)
        dimension gam(10),dsma(10),aa(MAXnl)
        dimension awidth(10),n1x(10),n1y(10),aslength(10),aswidth(10)
        double precision sma(10)
        double precision dum, smom, s,sdr,smat,rsmall,rlarge,sdrsa
	character*120 slipfile,file_gmt,velmodel
        pi=3.1415926
        pi2=pi/2.

c	data rperd/0.017453292/

	read(5,'(a256)')slipfile
        read(5,*)emag
        read(5,*)xlen,ylen,stra,dipa,rakea
        read(5,*)alon,alat,zz
        read(5,*)dx,dy
c  hypo coordinates along strike
        read(5,*)x0,y0,dens,vs
        read(5,*)dt
        read(5,*)iseed
        read(5,'(a256)')velmodel
        read(5,*)vel_fract
     
        emag_thre=5.5 
        fc = 6.0
c        if(xlen.le.25.and.ylen.lt.15.) then
        if(xlen.le.25) then
         nasp=1
        else
          if(xlen.lt.70) then
           nasp=2
          else
           nasp=3
          endif
        endif

         print*,'Number of asperities :',nasp

        i5=index(slipfile,' ')-1
        file_gmt=slipfile(1:i5)//'_gmt'

        call convert(alon,alat,xsc,ysc)
c         print*,'xsc',xsc,'ysc',ysc

        dipa10=dipa*pi/180.
        open(15,file=velmodel)
        za=0
        read(15,*)nvel
        do i=1,nvel
        read(15,*)zwdth,vp,vel(i),den(i),qp,qs
        za=za+zwdth
c        zvel(i)=za*sin(dipa10)
        zvel(i)=za  
        enddo
        close(15)

c ======    Irikura Recipe ===========

c  JMA magnitude and Seismic Moment in Nm
         emag_jma = (log10(xlen)+2.9)/0.6
c         smom = 10**(1.17*emag_jma+10.72)
c         write(6,208)smom
c208     format(1x,'Mo= ',e9.3,'Nm')

         width=ylen
c         print*,'Fault Length',xlen
c         print*,'Fault Width',ylen

333      s=xlen*ylen
         write(6,206)s
206      format(1x,'Fault Area Used in the Simulation   =',f8.1,'(km2)')
         if(s.gt.400) then
           if(s.gt.1800) then
           print*,'S>1800m^2'
c Murotani et al. 2015
           smom=s*1.e+17
           write(6,209)smom
209     format(1x,'Mo using Murotani  et al= ',e9.3,'(Nm)')
           else
         print*,'S>400m^2'
c  Irikura and Miyake (2001)
         smom= (s/4.24*1.0e11)**2 * 1.0e-07
         write(6,208)smom
208     format(1x,'Mo using Irikura and Miyake= ',e9.3,'(Nm)')
          endif

         else
         print*,'S<=400m^2'
c  Somerville et al (1999)
         ex=3./2.
         smom= (s/2.23*1.0e15)**ex * 1.0e-07
         write(6,204)smom
204     format(1x,'Mo using Somerville et al. = ',e9.3,'(Nm)')
         endif
         print*,'Fault Area=',s,'km^2'
         print*,'Seismic Moment =',smom,'Nm'

c       Magnitude
          emag0 = 1/1.5*log10(smom*1.e+07)-10.7

c   M<= 5.5
          if(emag.le.emag_thre) then
          emag0=emag
          d34 = 1.5*(emag+10.7)
          smom = 10**(d34)* 1.0e-07
          print*,emag0
          write(6,205)smom
205     format(1x,'Mo using Magnitude = ',e9.3,'(Nm)')

          nasp=1
c asperity area = fault area
          s=xlen*ylen
          smat= s
          rasp=sqrt(smat/pi)

c  average slip
         shearmod=dens*vs*vs
         d=smom/(shearmod*s)*1.e-15
         print*,'Average slip (m)',d

c asperity average slip
         dsmat=d
         print*,'Asperity average slip (m)',dsmat

c  slip in each asperity

        sma(1)=smat
        dsma(1)=dsmat
        print*,'Slip in asperity',dsma(1),'(m)'
        dback=dsmat
         write(6,210)dback

c  Effective stress in the asperity area (MPa)

         rlarge=sqrt(s/pi)
         dum=rasp*rasp*rlarge
         sdrsa=7./16.*smom/dum
         sdrsa=sdrsa*1.0e-15
         sdrback =0 

        print*,'!!!!!!!!!'

         else


        print*,'Magnitude ', emag0
c        print*,'JMA mag',emag_jma
c low frequency level of the acceleration source spectrum
c Nm/s2
         a=2.46*(smom*1.e+07)**(1./3.)
         a=a*1.0e+10

c         print*,'Short Period Level ',a,'Nm/s2'
c    equivalent radius of the fault area  in km
         rlarge=sqrt(s/pi)
         print*,'Equivalent Fault Radius (km)',rlarge

c         print*, 'vs ',vs
c    equivalent radius of the asperity area  in km
         rasp = 7*pi/4.0 *smom /a/rlarge*vs*vs
         print*,'Equivalent Asperity Area Radius (km)',rasp

c   asperity total area
         smat=pi*rasp*rasp
         if(nasp.gt.2) then
          smat=0.22*s
          rasp=sqrt(smat/pi)
         endif
         print*,'Equivalent Asperity Area Radius (km)',rasp
         write(6,207)smat
207      format(1x,'Asperity area (km2)=',f8.1)

c  average slip
         shearmod=dens*vs*vs
         d=smom/(shearmod*s)*1.e-15
         print*,'Average slip (m)',d

c asperity average slip
         dsmat=2*d
         print*,'Asperity average slip (m)',dsmat

c  slip in each asperity
         
        sma(1)=smat
        if(nasp.eq.2) then
           sma(1)=16./22.*smat
           sma(2)=6./22.*smat
c           sma(1)=2./3.*smat
c           sma(2)=1./3.*smat
        endif
        if(nasp.gt.2) then
           sma(1)=2.0/4.*smat
           sma(2)=1.0/4.*smat
           sma(3)=1.0/4.*smat   
        endif

        do i=1,nasp
         r1=(smat/pi)**0.5
         r2=(sma(i)/pi)**0.5
         gam(i)=r2/rasp
c         print*,'gama',i,gam(i)
        enddo
         ss=0
        do i=1,nasp
         ss=ss+gam(i)**3
        enddo
        do i=1,nasp
        dsma(i)=gam(i)/ss*dsmat
        print*,'Slip in asperity',i,dsma(i),'(m)'
        enddo
c average slip in the background area

         smoma=shearmod*dsmat*smat*1.0e+15
         smomback=smom-smoma
         sback=s-smat
         dum1 = smomback*1.e-15
         dum2 = shearmod*sback
         dback=dum1/dum2
         write(6,210)dback
210     format(1x,'Background slip (m)=',f6.2)

c  Effective stress in the asperity area (MPa)

         dum=rasp*rasp*rlarge
         sdrsa=7./16.*smom/dum
         sdrsa=sdrsa*1.0e-15
         if(nasp.gt.2) sdrsa=14.1
         print*,'Effective stress drop in asperity area (MPa):'
     +     ,sdrsa

c  effective stress in the background area

          if(nasp.eq.1) then
c Dan et al. (2002)
         sdrback=dback/dsmat/width*sqrt(smat)
         sdrback=sdrback*sdrsa
         else
         sdrback=(dback/dsmat)/width*sqrt(pi)*rasp*ss
         sdrback=sdrback*sdrsa
         endif

         print*,'Background stress drop (MPa):',sdrback

          dum=rlarge**3
          dum2=7./16.*smom/dum
          if(nasp.gt.2) dum2=3.1*1.0e+15
          print*,'Average stress drop =', dum2*1.0e-15,'MPa'
c  =================================================================
c   dback  :  background slip (cm)
c   dsma(i):  slip in each asperity (cm)
c   sdrback:  background stress drop (Mpa)
c   sdrsa :   asperity stress drop (MPa)
c   convert slip to  cm

c and if on magnitude threshold
      endif

        dback=dback*100.
        sdrsa=sdrsa*10
        sdrback=sdrback*10
        vss =vs

c  output of the slip velocity functions is made at (lx,ly) locations

        nx=int(xlen/dx)
        ny=int(ylen/dy)
c        print*,'nx,ny = ',nx,ny

        nx0=int(x0/dx)
        ny0=int(y0/dy)
        area=dx*dy*1.0e+10
c    fault element area is in cm*cm

       aspect=ylen/xlen
       aspect=0.8
c if nasp<3 aspect = 1
c if nasp=3 aspect =0.8

        if(nasp.le.2) aspect=1.0
        do i=1,nasp
       aleng=sqrt(sma(i)/aspect)
       awidth(i)= aleng*aspect
        n1x(i)=int(aleng/dx)
        n1y(i)=int(awidth(i)/dy)
        aslength(i)=aleng
        aswidth(i)=awidth(i)

        enddo
        do i=1,ny
         do j=1,nx
          slip(i,j)=dback
          stres(i,j)=sdrback
          alng(i,j)=ylen
         enddo
        enddo
c      random location of asperities
       if(emag.gt.emag_thre) then
       call as_loc (iseed,xlen,ylen,zz,dx,dy,nx,ny,nasp,aslength,
     + aswidth, xasp,yasp)
       else
        xasp(1)=0
        yasp(1)=0
       endif

         do i=1,nasp
           i1=int(xasp(i)/dx)
           j1=int(yasp(i)/dy)
           i2=i1+n1x(i)
           j2=j1+n1y(i)

           do j=j1,j2
            do k=i1,i2
              slip(j,k)=dsma(i)*100.
              stres(j,k)=sdrsa
              alng(j,k)=awidth(i)
            enddo
           enddo

          if(i.eq.1) then
            lx1=j1+1
            ly1=i1+1
           endif
           if(i.eq.2) then
            lx2=j1+1
            ly2=i1+1
           endif

          enddo

c  rescale slip to match the seismic moment using the 1D velmod
        z=zz-dy/2.

        do 42 j=1,ny
          z=z+dy
           aa(j)=den(1)*vel(1)*vel(1)
           do k=1,nvel
           if(zvel(k).gt.z*sin(dipa10)) then
           aa(j)=den(k)*vel(k)*vel(k)
           beta(j)=vel(k)
           rden(j)=den(k)
           goto 42
           endif
           enddo
            print*,'vel model problem !'
            stop
42        continue

           sdum=0
         do j=1,ny
           do i=1,nx
           sdum=sdum+aa(j)*slip(j,i)*dx*dy
           enddo
         enddo

          sdum=smom*1.0e07/sdum*1.0e-20
 
          print*,'Seis Mom Scaling Factor:',sdum
        open (12,file=slipfile)
c        open (13,file=file_gmt)
        open (14,file='stf_background.out')
        open (16,file='stf_asperity1.out')
        open (26,file='stf_asperity2.out')
        open (27,file='stress_drop.out')

        version=1.0
        nums=1
        write (12,'(f3.1)')version
        write (12,'("PLANE ",i2)')nums
        write (12,815) alon,alat,nx,ny,xlen,ylen
        mstra =nint(stra)
        mdipa = nint(dipa)

        write (12,820) mstra, mdipa, zz, x0,y0
815     format (f11.6,1x,f11.6,2(1xi6),2(f10.4))
820     format (1x,i3,1x,i3,3(1x,f10.4))
c  CORNER FREQUENCY
        shearmod = dens*vss*vss
c       td = 0.1
c       tb = 1.25*td
c        td=2./fc/pi
c        tb = 1.3*td

        write (12,'(a,i8)') "POINTS",nx*ny

       dl=(nx)*dx
       dw=(ny)*dy
       nx2=nx/2.
       dx2=0.5*dx
       dy2=0.5*dy
       strike=stra*pi/180.
       dip=dipa*pi/180.
       rake=rakea*pi/180.
c       ylat0=alat+dl/2*cos(strike)/ysc
c       xlon0=alon+dl/2*sin(strike)/xsc

        k=0
        do j=1,ny
         y=(j-1+dy2)*dy*cos(dip)
         z=zz+(j-1+dy2)*dy*sin(dip)
         call rupture_vel(z,nvel,zvel,k10)
c         rupvel=vel(k10)*0.8
c         rupvel=vss*0.72
         rupvel=vss*vel_fract

         shearmod=den(k10)*vel(k10)*vel(k10)
           do i=1,nx
            k=k+1
         x=(i-1+dy2)*dx
         str2=strike
         xx=alon+(x-dl/2)*sin(strike)/xsc+y*cos(str2)/xsc
         yy=alat+(x-dl/2)*cos(strike)/ysc-y*sin(str2)/ysc
c         if(j.eq.1)print*,yy,xx,alon,dl/2-x,sin(strike), y*cos(str2)/xsc

         d1=abs(i-(nx2+nx0))*dx
         d2=abs(j-ny0)*dy
         stim=sqrt(d1*d1+d2*d2)/rupvel
              slip1 = slip(j,i)*sdum
              sdrop = stres(j,i)
              slip2 = 0.0
              slip3 = 0.0
              zero  = 0.0
              sdrop = stres(j,i)

c  output stress drop in MPa 
         write(27,*)i,j,sdrop/10.

c  rise time for each asperity
             tr =alng(j,i)/rupvel/2.

             td =1./fc/pi
             tb = 1.3*td
             call srctimfunc(MAXnw,slip1,dt,sdrop,alng(j,i),rupvel,fc,
     &           shearmod,td,tb,tr,ty,npts)

              mstra=nint(stra)
              mdipa=nint(dipa)
              mrak=nint(rakea)

        write (12,200)xx,yy,z,mstra,mdipa,area,stim,dt,
     +    beta(j)*100000.0,rden(j)
        write (12,202) mrak,slip1,npts,slip2,int(zero),
     +    slip3,int(zero)
        write (12,400) (ty(it),it=1,npts)
c        write (12,400) (ty(it),it=1,npts)
c        write(13,*)xx,yy,z,(i-1)*dx, -(j-1)*dy,slip1,stim


        if(j.eq.ny-1.and.i.eq.nx-1) then
           do l=1,npts
           write(14,*)(l-1)*dt,ty(l)
           enddo
        endif
        if(j.eq.lx1.and.i.eq.ly1)then
          do l=1,npts
          write(16,*)(l-1)*dt,ty(l)
          enddo
        endif
        if(j.eq.lx2.and.i.eq.ly2)then
          do l=1,npts
          write(26,*)(l-1)*dt,ty(l)
          enddo
        endif

          enddo
        enddo

200     format(1x,f11.6,1x,f11.6,1x,e12.5,2(1x,i3),1x,
     + e13.5,1x,f8.4,3(1x,e12.5))
202     format(i4,3(1x,f10.4,i5))
300     format(f8.2,3(e13.5,i8))
400     format(6e13.5)
        end


	
	subroutine set_g2 (g2,fc)
	real g2,fc,f
	f = 1.0/fc
c	g2 = (2.0*f - f*f) / ((1.0-f)*(1.0-f))
	g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f))
	return
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

        subroutine srctimfunc(n,slip,dt,sd,aleng,rupvel,fc,shearmod,
     +                            td,tb,tr,st,npts)
c       include 'limits.h'
        dimension st(1)

c  sd:     stress drop (MPa)
c  aleng : asperity characteristic length (km)
c  rupvel: rupture velocity (km/s)
c  fc    : corner frequency
c  shearmod : shear modulus
c  vm    : slip velocity (cm/s)
c  slip  : (cm)
c  npts


        ts=1.5*tr
        vm=sd*sqrt(2*aleng*rupvel*fc)/shearmod
        kk = 1
9876    continue

        eps=(5*tb-6*td)/(1-td/tb)/4.
        b=2.*vm/td*tb*(tb-eps)**0.5*(1-tb/td/2)
c        print *,'Max slip velocity ',vm,' m/s'
c       print *,"td,tb,tr,b,eps:  ",td,tb,tr,b,eps
c       print *," shearmod,aleng,rupvel,fc: ",shearmod,aleng,rupvel,fc

        do 10 i=1,n
         tt=(i-1)*dt
        if(tt.le.tb) then
           st(i)=2.*vm/td*tt*(1-tt/2/td)
        else
           if(tt.le.tr) then
             st(i)=b/(tt-eps)**0.5
             c=st(i)
           else
             if(tt.le.ts) then
               st(i)=c-c/(ts-tr)*(tt-tr)
             else
               st(i)=0.
               npts = i
               goto 20
             endif
           endif
        endif
10      continue
20      continue
        s = 0
        do i=1,npts
           s = s + st(i)*dt
           disp = s
        enddo
        fact = slip/disp
        do i=1,npts
           st(i) = st(i)*fact
        enddo
c       print *,kk,"Factor= ",fact, " Max slip= ",disp,slip
c       print *,"npts=   ",npts
        if(fact.gt.1.01.and.kk.lt.2) then

           tb = tb*0.99
           kk = kk + 1
           go to 9876

        else if(fact.lt.0.99.and.kk.lt.2) then

           tb = tb*1.01
           kk = kk + 1
           go to 9876

        endif

        return
        end



        subroutine rupture_vel(z,nvel,zvel,k)
        dimension zvel(1)
        zm=0
        do i=1,nvel
        if(z.lt.zvel(i)) then
         k=i
         goto 1
        endif
        enddo
        k=nvel
1       return
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

       subroutine as_loc (iseed,fl,fw,zz,dx,dy,nx,ny,n,aslength,
     + aswidth, xa,ya)
       dimension xa(4),ya(4),aslength(4),aswidth(4)
       dimension ik(800,400),asl(4),asw(4)

55       iseed=iseed+20

       xm=rand(iseed)
       xr=rand(0)
       if(xr.le.0.1)xr=xr+0.1
       if(xr.ge.0.9)xr=xr-0.1

       yr=rand(0)
       if(yr.le.0.1)yr=yr+0.1
       if(yr.ge.0.9)yr=yr-0.1

c  asperities length and width
       do i=1,n
       asl(i)=aslength(i)+dx
       asw(i)=aswidth(i)+dy
       enddo

c space remaining to the left and below the first asperity
c location of the first asperity xa(1),ya(1)

c  asperity location should be below 2km from the top edge
c  asperity location should be below 2km from the earth surface

       d9=2.0
       dh=d9-zz
       if(dh.le.0) dh=0.
       x1=(fl-asl(1)-dx)*xr
       y1=(fw-asw(1)-dy)*yr
       if(x1.lt.dx) x1=2*dx
c       if(y1.lt.dy) y1=dy
       if(y1.lt.dh) y1=dh
       xa(1)=x1
       ya(1)=y1

       yr=rand(0)
c       xr=rand(0)

       print*,xr
       print*,yr

       if(n.eq.1) then
c 1 asperity
       return
       else
c 2 asperities

       if(xr.lt.0.5) then
c X start left
        xx=x1-asl(2)
        if (xx.gt.asl(2)) then
        xa(2)=xx*xr
        else
c X right
        xx=fl-x1-asl(1)-asl(2)
        xa(2)=x1+asl(1)+xx*xr
        endif
c Y
c        yy=fw-asw(2)
c        ya(2)=yy*yr
        yy=fw-dh-asw(2)
        ya(2)=yy*yr+dh-dy
c  ===========================
c  xr > 0.5
         else

c X start right
        xx=fl-x1-asl(1)-asl(2)
        if(xx.gt.asl(2)) then
         xa(2)=x1+asl(1)+xx*xr
        else
        xx=x1-asl(2)
        if(xx.le.0) goto 55
         xa(2)=xx*xr
         endif
c Y
c        yy=fw-asw(2)
c        ya(2)=yy*yr+dy
        yy=fw-dh-asw(2)
        ya(2)=yy*yr+dh-dy
        endif
        endif

c==========================

       if(n.eq.3) then

       do i=1,nx
        do j=1,ny
         ik(i,j)=0
        enddo
       enddo

       i1=nint(asl(1)/dx)+1
       k1=nint(asw(1)/dy)+1

       i2=nint(asl(2)/dx)+1
       k2=nint(asw(2)/dy)+1

       i3=nint(asl(3)/dx)+1
       k3=nint(asw(3)/dy)+1

       ii1=nint(xa(1)/dx)-1
       kk1=nint(ya(1)/dy)+1
       ii2=nint(xa(2)/dx)-1
       kk2=nint(ya(2)/dy)+1

       if(ii1.eq.0) ii1=1
       if(ii2.eq.0) ii2=1

       iend=ii1+i1+1
       kend=kk1+k1

c  ik(i,k)=1 in grids within asperity area
c  ik(i,k)=0 in grids ouside asperity area

c  Asperity 1
       do i=ii1,iend
        do k=kk1,kend
          ik(i,k)=1
        enddo
       enddo

c Asperity 2
       iend=ii2+i2+1
       kend=kk2+k2
       do i=ii2,iend
        do k=kk2,kend
          ik(i,k)=1
        enddo
       enddo
c  top of the fault 
         i7=int(dh/dy)+1

         do i=1,nx
           do j=1,i7
            ik(i,j)=1
           enddo
         enddo

         nyend=ny-k3-1
         nxend=nx-i3-1
         n4=nxend

122      xr=rand(0)
         yr=rand(0)

c  closest distance to  the fault edge set to 5.0 km
         g5=int(5.0/dx)*xr
         nxend=n4-int(g5)
         if(xr.lt.0.5) then
          nxend=nxend-3
         endif
         if(nxend.eq.0)nxend=1

          print*,xr

         xd=nxend*xr
         yd=nyend*yr
         ixd=nint(xd)
         iyd=nint(yd)

         do 90 i5=nxend,ixd,-1
         do 80 k5=iyd,nyend

        do i=1,i3
         i6=i5+i
        do k=1,k3
        k6=k5+k
         if(ik(i6,k6).gt.0) goto 80
        enddo
        enddo
       goto 100
80     continue
90     continue
       goto 122
100    xa(3)=float(i5)*dx
       ya(3)=float(k5)*dy
       endif
       return
       end
 
