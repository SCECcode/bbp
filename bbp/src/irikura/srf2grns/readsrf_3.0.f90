subroutine readsrf

! 2019.5
! version3.0 = for multi-segment

  use params
  implicit none
  integer::i,nt1,nt2,nt3,j,nseg,k
  real::elon,elat,stk,dip,area,slp1,slp2,slp3,dt,dp,a,b
  real::vsave,densave
  real,dimension(1000)::d1,d2,d3
  character::chara*7

  ! SCEC BBP Standard Rupture Format File
  open(10,file=srffile,action='read')
  ! SRF version = 2.0
  read(10,*)
  ! Header block
  read(10,*) chara,nseg  ! number of segments
  if(nseg >= 2) then
     write(6,*) "ERROR : this program reads single-segment SRF only"
     stop
  endif

  read(10,*) elon,elat,nstk,ndip,len,wid
  read(10,*)
  dx=len*1000./float(nstk) ! subfault length (m)
  dy=wid*1000./float(ndip) ! subfault width (m)
     ! Data block
  read(10,*) chara,np
  write(6,*) 'SRF  number of subfaults : ',np, dx, dy
  
  allocate(lon(np),lat(np),dep(np),tinit(np),slip(np),sdrop(np))
  
  d1=0.0e0 ; d2=0.0e0 ; d3=0.0e0
  vs=0.0e0 ; dens=0.0e0
  musum=0.0

  do i=1,np
     read(10,*) lon(i),lat(i),dp,stk,dip,area,tinit(i),dt,a,b
     read(10,*) rake,slp1,nt1,slp2,nt2,slp3,nt3
     !     write(6,*) ' -- ',i,np
     if(nt1 > 0) then
        read(10,*) (d1(j),j=1,nt1)
     endif
     if(nt2 > 0) then
        read(10,*) (d2(j),j=1,nt2)
     endif
     if(nt3 > 0) then
        read(10,*) (d3(j),j=1,nt3)
     endif
     
     dep(i)=dp*1000.    ! depth (km)
     slip(i)=slp1*0.01 ! slip (m) (in u1 direction)
     !     vs=vs+a*0.01 ; dens=dens+b*1000  ! m/s, kg/m2
     
     ! 20180822
     musum=musum+b*a*a*0.1*dx*dy  ! sum of mu
     
  enddo
  !  vsave=vs/float(np) ; densave=dens/float(np) ! average over the fault
  !  write(6,*) 'VS ave, Dens ave = ', vsave,densave
  
  ! 20180822 temporary values for vs and dens
  vs=3500 ; dens=2700
  write(6,*) 'VS(m/s), Density(kg/m^2) = ', vs, dens

  close(10)


  return

end subroutine readsrf

