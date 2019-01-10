subroutine sdrop(len,wid,dens,vs,musum,sdrsa,sdrback,dback)

 implicit none
 real::pi,len,wid,rlarge,rasp,ex,a,mu,vs,dens,musum
 real::smom,smoma,smomback   ! Seismic moment
 real::scale  ! scale to adjust seismic moment 20180822
 real::s,smat,sback,dsmat,dback,sdrsa,sdrback  ! area,slip,sigma
 real::gam1,gam2,gam3,ss,dback0
 integer::nasp

 
  pi=acos(-1.0)
  s=len*wid  ! fault area km2
  rlarge=sqrt(s/pi)
  ! vs (m/s), dens(kg/m3), mu(N/m2)
  mu=dens*vs*vs

! Seismic Moment Nm
  if(s.le.400) then
     ! Somerville et al. 1999
     ex=3./2. ; smom=((s/2.23*1.0e15)**ex)*1.0e-7
  else if(s.le.1800) then
     ! Irikura & Miyake 2001
     smom=((s/4.24*1.0e11)**2)*1.0e-7
  else
     ! Murotani et al. 2015
     smom=s*1.0e17
  endif

  if(len.le.25.and.wid.lt.15) then
     nasp=1
  else
     if(len.lt.70) then
        nasp=2
        ! Sa1:Sa2=16:6
        gam1=sqrt(16./22.) ; gam2=sqrt(6./22.)
        ss=gam1**3+gam2**3
     else
        nasp=3
        ! Sa1:Sa2:Sa3=2:1:1
        gam1=sqrt(2./4.) ; gam2=sqrt(1./4.) ; gam3=sqrt(1./4.)
        ss=gam1**3+gam2**3+gam3**3
     endif
  endif

  write(6,*) ' Number of Asperities : ',nasp

  ex=1./3. ; a=(2.46*(smom*1.0e7)**ex)*1.0e10 ! short-period level
  rasp=(7.*pi/4.0)*smom/a/rlarge*vs*vs*1.0e-6 ! asp radius (km)
  smat=pi*rasp*rasp ! total asperity area (km^2)
  if(nasp.ge.3) then
     smat=0.22*s
     rasp=sqrt(smat/pi)
  endif
  dsmat=2.0*smom/mu/(s*1.0e6) ! average slip (m) in asperities
  smoma=mu*dsmat*smat*1.0e6  ! total moment of asperities (Nm)
  smomback=smom-smoma ! background moment (Nm)
  sback=s-smat  ! background area(km2)
  dback=smomback/mu/(sback*1.0e6)  ! slip(m) in background


  scale=mu*s*1.0e6/musum
  write(6,*) '  scale = ',scale
  dback0=dback
  dback=dback*scale ! 20180822  ! adjust to seismic moment

  ! Effective stress drop in asperities
  sdrsa=(7./16.)*smom*1.0e-15/(rasp*rasp*rlarge) ! MPa
  ! 20180904 adjust difference in Vs and Dens
  if(scale>1.0) then
     sdrsa=sdrsa*scale
  endif
  ! Three asperities : Sa = S x 22%
  if(nasp.ge.3) then
     sdrsa=14.1 ! 3.1MPa/22%
  endif

  if(nasp.eq.1) then
     sdrback=dback0/dsmat/wid*sqrt(smat)*sdrsa  ! MPa
  else
     sdrback=dback0/dsmat/wid*sqrt(pi)*rasp*ss*sdrsa
  endif
  write(6,*) ' Seismic Moment (Nm) : ', smom
  write(6,*) ' Short-period level (Nm/s/s) : ',a
  write(6,*) ' Total asperity area (m^2) : ',smat
  write(6,*) ' Asperity stress drop (MPa) : ',sdrsa
  write(6,*) ' Asperity average slip (m) : ',dsmat
  write(6,*) ' Background stress drop (MPa) : ',sdrback
  write(6,*) ' Background moment (Nm) : ' ,smomback
  write(6,*) ' Background slip (m) : ' ,dback

  return
end subroutine sdrop
