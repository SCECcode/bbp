Program stitch
! This code stitches the Low- and High-frequency synthetic 
! ground motions together in wavelet domain.
!
! Here assumed that the inputted ground motions are velocities.
!
! Written by Pengcheng Liu
! Copyright (c) 2007 by Pengcheng Liu
!
 implicit NONE
 character(len=256):: fileinp,fileRec_list,fileVM_list,fileVM, &
                      st_name,path1D,path3D
 integer:: num_sm,id_motion,nst,ks,ism,getlen
 real:: freq_join,freq_max,depth_hyp
 real:: st_x,st_y,st_z,stan_h1,stan_h2,xx00,yy00,afa
 real:: dis_x,dis_y,dist,time0
!
 fileinp='stitch.inp'
 open(unit=8,file=fileinp,status='old')

 write(*,*) 'Enter the name of file listing the stations'
 read(8,'(1a)') fileRec_list

 write(*,*) 'Enter the name of file listng 1D VM for each stations'
 read(8,'(1a)') fileVM_list

 write(*,*) 'Enter the name of folder(and Path) with 1D synthetics'
 read(8,'(1a)') path1D

 write(*,*) 'Enter the name of folder(and Path) with 3D synthetics'
 read(8,'(1a)') path3D

 write(*,*) 'Enter the joint frequency and fmax '
 read(8,*) freq_join, freq_max

 write(*,*) 'Enter the depth of hypocenter (km)'
 read(8,*) depth_hyp

 write(*,*) 'Enter the number of source model'
 read(8,*) num_sm

 write(*,*) 'Outptut Displacement (=1), Vel. (2), or Acc (3)'
 read(8,*) id_motion

 close(8)
!============================================================

 open(10,file=fileVM_list,status='old', position='rewind')
 open(9,file=fileRec_list,status='old', position='rewind')
 read(9,*) nst,xx00,yy00,afa
 do ks=1,nst
   read(10,'(1a)') fileVM
   call cmodel(fileVM,depth_hyp)
!
   read(9,'(1a)') st_name
   read(9,*) st_x,st_y,st_z,stan_h1,stan_h2 
   dis_x=(st_x-xx00)/1000.0
   dis_y=(st_y-yy00)/1000.0
   dist=sqrt(dis_x**2+dis_y**2) 
   call trav(dist,time0)
!   write(*,*) ks,time0
!   write(*,*) st_name
   do ism=1,num_sm
     call ReadVel(ism,stan_h1,stan_h2,st_name,path1D,path3D)
     call align(freq_join,time0)
     call stitch_wt(freq_join)
     call time_history(id_motion,freq_max)
     call output3(ism,stan_h1,stan_h2,st_name)
   enddo  ! end loop of source model (ism)
 enddo  ! end loop of station (ks)
 close(9)
 close(10)
 stop
end Program stitch
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
MODULE ccmodel
 implicit NONE
 save
 integer:: ns,lmax,nlayer,nb
 real:: freq_ref
 real, allocatable, dimension(:):: afa,beta,rho,thick,qp,qs,dep
 real(kind=8), allocatable, dimension(:):: th,v
END MODULE ccmodel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine cmodel(fileVM,depth)
 use ccmodel
 implicit none
 character(len=*), intent(IN):: fileVM
 real, intent(IN):: depth
 integer:: jo,j1,j
!
 open(unit=19,file=fileVM,status='old',position='rewind')
 read(19,*) jo, freq_ref
 if(allocated(beta)) deallocate(afa,beta,rho,thick,qp,qs,dep,th,v)
 j1=jo+1
 allocate(afa(j1),beta(j1),rho(j1),thick(j1),qp(j1),qs(j1))
 allocate(th(j1),v(j1),dep(j1+1))
 do j=1,jo
   read (19,*) afa(j),beta(j),rho(j),thick(j),qp(j),qs(j)
 enddo 
 close(19)
!
 dep(1)=0.0
 do j=1,jo
   dep(j+1)=dep(j)+thick(j)
 enddo 
 dep(jo+1)=amax1(dep(jo+1),depth+10.0)
 nb=2
 do j=1,jo
   if(depth>dep(j)) nb=j+1
 enddo
 if(abs(depth-dep(nb)) < 1.e-5) then
   nlayer=jo
 else if(abs(depth-dep(nb-1)) < 1.e-5) then
   nlayer=jo;   nb=nb-1
 else
   do j=jo,nb-1,-1
     j1=j+1
     afa(j1)  =afa(j)
     beta(j1) =beta(j)
     rho(j1)  =rho(j)
     thick(j1)=thick(j)
     qp(j1)   =qp(j)
     qs(j1)   =qs(j)
   enddo
   thick(nb-1)=depth-dep(nb-1)
   thick(nb)  =dep(nb)-depth
   nlayer=jo+1
 endif
end subroutine cmodel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
MODULE synths_comm 
 implicit NONE
 save
 integer:: nt1,nt3,npt,npt2
 real:: dt1,dt3,dt
 real, allocatable, dimension(:,:):: vel1d,vel3d
END MODULE synths_comm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ReadVel(ism,stan_h1,stan_h2,st_name,path1D,path3D)
 use synths_comm
 implicit NONE
 character(len=*), intent(INOUT):: st_name,path1D,path3D
 integer,intent(IN):: ism
 real,intent(IN):: stan_h1,stan_h2
 character(len=256):: fh1,fh2,fup,chsm,smidx,prefix,suffix
 character(len=3):: ch1,ch2,ch3
 integer:: nc_st,nc_p1,nc_p3,nc_sm,getlen,ntt,i,k
 real, allocatable, dimension(:,:):: syn3d
!
 call num2ch(float(ism)+0.01,smidx)
 call num2ch(stan_h1,ch1)
 call num2ch(stan_h2,ch2)
 ch3='ver'
 nc_p1=getlen(path1D)
 nc_p3=getlen(path3D)
 nc_st=getlen(st_name)
 nc_sm=getlen(smidx)
! Input 1D synthetics
 chsm='.gm1D.'//smidx(1:nc_sm)
 k=getlen(chsm)
 fh1=''; fh2=''; fup=''
 fh1=path1D(1:nc_p1)//st_name(1:nc_st)//'.'//ch1//chsm(1:k)
 fh2=path1D(1:nc_p1)//st_name(1:nc_st)//'.'//ch2//chsm(1:k)
 fup=path1D(1:nc_p1)//st_name(1:nc_st)//'.'//ch3//chsm(1:k)
 open(21,file=fh1,status='old', position='rewind')
 open(22,file=fh2,status='old', position='rewind')
 open(23,file=fup,status='old', position='rewind')
 read(21,*) nt1,dt1
 read(22,*) ntt,dt1
 nt1=min0(nt1,ntt)
 read(23,*) ntt,dt1
 nt1=min0(nt1,ntt)
 if(allocated(vel1d)) deallocate(vel1d)
 allocate(vel1d(nt1,3))
 read(21,*) (vel1d(i,1), i=1,nt1)
 read(22,*) (vel1d(i,2), i=1,nt1)
 read(23,*) (vel1d(i,3), i=1,nt1)
 close(21);  close(22); close(23)
!
! Input 3D synthetics !JORGE changed gm3D to gm1D on line below
 chsm='.gm1D.'//smidx(1:nc_sm)
 k=getlen(chsm)
 fh1=''; fh2=''; fup=''
 fh1=path3D(1:nc_p3)//st_name(1:nc_st)//'.'//ch1//chsm(1:k)
 fh2=path3D(1:nc_p3)//st_name(1:nc_st)//'.'//ch2//chsm(1:k)
 fup=path3D(1:nc_p3)//st_name(1:nc_st)//'.'//ch3//chsm(1:k)
 open(21,file=fh1,status='old', position='rewind')
 open(22,file=fh2,status='old', position='rewind')
 open(23,file=fup,status='old', position='rewind')
 read(21,*) nt3,dt3
 read(22,*) ntt,dt3
 nt3=min0(nt3,ntt)
 read(23,*) ntt,dt3
 nt3=min0(nt3,ntt)
 if(allocated(vel3d)) deallocate(vel3d)
 allocate(vel3d(nt3,3))
 read(21,*) (vel3d(i,1), i=1,nt3)
 read(22,*) (vel3d(i,2), i=1,nt3)
 read(23,*) (vel3d(i,3), i=1,nt3)
! The vertical component from FD is positive down.
! but FK is positive UP
! vel3d(:,3)=-vel3d(:,3)  !JORGE not necessary
 close(21);  close(22); close(23)
!
! Interpolating to the same time step
 npt=ifix((nt3-1)*dt3/dt1)
 ntt=max0(npt,nt3)
 allocate(syn3d(ntt,3))
 do k=1,3
   syn3d(1:nt3,k)=vel3d(1:nt3,k)
   if(abs(dt1-dt3).gt.0.01*dt1) then
     if(dt1.gt.dt3) call xapiir(syn3d(1,k),nt3,'BU',0.0,0.0,4, &       
                             'LP', 0.001, 0.48/dt1, dt3, 2, ntt) 
     call interpol(0.0,dt3,nt3,syn3d(1,k),dt1,npt, ntt)
   endif
 enddo
 npt=min0(npt, nt1)
 deallocate(vel3d);   allocate(vel3d(npt,3))
 vel3d(1:npt, :)=syn3d(1:npt, :)
 nt3=npt;  dt3=dt1
 npt2=2
 do while (npt2<npt*2)
   npt2=npt2*2
 enddo
 deallocate(syn3d)
end subroutine ReadVel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine align(freq_join,ts_start)
 use synths_comm
 implicit NONE
 real, intent(IN):: freq_join,ts_start
 integer:: i,j,k,ntm1,ntm2,nwd1,nwd2,ish
 real:: vmx1,vmx2,crr,cross,sum1,sum2,sum3,shift_half
 real, dimension(npt):: vh1,vh2

 vh1=vel1d(1:npt,1);  vh2=vel1d(1:npt,2)
 shift_half=amin1(3.0,0.9*ts_start,amax1(0.1*ts_start,0.5))
 ntm1=-int(shift_half/dt1+0.5)
 nwd1=max0(int(amax1(0.0,(ts_start-1.0/freq_join))/dt1+1.5),1-ntm1)
 nwd2=min0(int((ts_start+2.0/freq_join)/dt1+1.5),npt)
 ntm2=min0(int(shift_half/dt1+0.5), npt-nwd2)
 call xapiir(vh1, npt, 'BU', 0.0, 0.0, 4, &       
             'BP',freq_join/2.0,freq_join,dt1,2,npt) 
 call xapiir(vh2, npt, 'BU', 0.0, 0.0, 4, &       
             'BP',freq_join/2.0,freq_join,dt1,2,npt) 
 vmx1=0.0;  vmx2=0.0
 do i=nwd1,nwd2
   vmx1=amax1(vmx1,abs(vh1(i)))
   vmx2=amax1(vmx2,abs(vh2(i)))
 enddo
 k=1;  if(vmx2 > vmx1) k=2
 vh1=vel1d(1:npt,k);  vh2=vel3d(1:npt,k)
 call xapiir(vh1, npt, 'BU', 0.0, 0.0, 4, &       
             'BP',freq_join/2.0,freq_join,dt1,2,npt) 
 call xapiir(vh2, npt, 'BU', 0.0, 0.0, 4, &       
             'BP',freq_join/2.0,freq_join,dt1,2,npt) 
 cross=0.0;  ish=0
 do j=ntm1,ntm2
   sum1=0.0;  sum2=0.0;  sum3=0.0
   do i=nwd1,nwd2
     sum1=sum1+vh1(i)*vh2(i+j)
     sum2=sum2+vh1(i)**2
     sum3=sum3+vh2(i+j)**2
   enddo
   crr=sum1/(sqrt(sum2)*sqrt(sum3))
   if(cross < crr) then
     ish=j
     cross = crr
   endif
 enddo  

 do k=1,3
   if(ish < 0) then
     do i=npt,1-ish,-1
       vel3d(i,k)=vel3d(i+ish,k)
     enddo 
     do i=1,-ish
       vel3d(i,k)=0.0
     enddo 
   else
     do i=1,npt-ish
       vel3d(i,k)=vel3d(i+ish,k)
     enddo 
     do i=npt-ish+1,npt
       vel3d(i,k)=0.0
     enddo 
   endif
 enddo
end subroutine align
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine stitch_wt(freq_join)
 use synths_comm
 implicit NONE
 real, intent(IN):: freq_join
 integer:: k,kj,nb,ne
 real:: factor
 real, dimension(npt2):: syn1d,syn3d
! 
 factor=real(npt2)*freq_join
 kj=1;  dt=2.0/factor
 do while(dt <= dt1)
   kj=kj+1
   dt=2.0**kj/factor
 enddo
 kj=kj-1
 dt=2**kj/factor
 ne=2**(kj+1);  nb=ne/2+1
 do k=1,3
   syn1d=0.0;   syn3d=0.0
   syn1d(1:npt)=vel1d(1:npt,k)
   syn3d(1:npt)=vel3d(1:npt,k)
   call interpol(0.0, dt1,npt2,syn1d,dt,npt2,npt2)
   call interpol(0.0, dt3,npt2,syn3d,dt,npt2,npt2)
   call fwd_wt(npt2,syn1d)
   call fwd_wt(npt2,syn3d)
   syn1d(1:ne)=syn3d(1:ne)
   call inv_wt(npt2,syn1d)
!   call xapiir(syn1d,npt2,'BU',0.0,0.0,4, &    !JORGE removing BP       
!               'LP',0.001,0.48/dt1,dt,2,npt2) 
   call interpol(0.0,dt,npt2,syn1d,dt1,npt,npt2)
   vel1d(1:npt, k)=syn1d(1:npt)
 enddo
end subroutine stitch_wt
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine time_history(id_motion)
 use synths_comm
 implicit NONE
 integer, intent(IN):: id_motion
 integer:: k
 do k=1,3
   if(id_motion == 1) then
     call vel2dis(nt1,dt1,vel1d(1,k))
     call vel2dis(nt3,dt3,vel3d(1,k))
   else if(id_motion == 3) then
     call vel2acc(nt1,dt1,vel1d(1,k))
     call vel2acc(nt3,dt3,vel3d(1,k))
   endif
 enddo
end subroutine time_history
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine output3(ism,stan_h1,stan_h2,st_name)
 use synths_comm
 implicit NONE
 character(len=*), intent(IN):: st_name
 integer, intent(IN):: ism
 real, intent(IN):: stan_h1,stan_h2
 character(len=72):: fout,chsm
 character(len=3):: ch1,ch2,ch3
 integer:: i,getlen,nc_st

 call num2ch(float(ism)+0.01,ch3)
 chsm='.gmBB.'//ch3
 ch3='ver'
 call num2ch(stan_h1,ch1)
 call num2ch(stan_h2,ch2)
 nc_st=getlen(st_name)
 fout=''; fout=st_name(1:nc_st)//'.'//ch1//chsm(1:9)
 open(19, file=fout,status='replace',position='rewind')
 write(19,*) npt,dt1
 write(19,925) (vel1d(i,1), i=1,npt)
 close(19)
 fout=''; fout=st_name(1:nc_st)//'.'//ch2//chsm(1:9)
 open(19, file=fout,status='replace',position='rewind')
 write(19,*) npt,dt1
 write(19,925) (vel1d(i,2), i=1,npt)
 close(19)
 fout=''; fout=st_name(1:nc_st)//'.'//ch3//chsm(1:9)
 open(19, file=fout,status='replace',position='rewind')
 write(19,*) npt,dt1
 write(19,925) (vel1d(i,3), i=1,npt)
 close(19)
 925   format(5e14.5) 
end subroutine output3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine vel2acc(nt,dt,acc)
 implicit NONE
 integer, intent(IN) :: nt
 real, intent(IN) :: dt
 real, intent(INOUT), dimension(nt) :: acc
 integer :: i
 do i=1,nt-1
   acc(i)=(acc(i+1)-acc(i))/dt
 enddo 
 acc(nt)=0.0
 do i=nt,2,-1
   acc(i)=0.5*(acc(i)+acc(i-1))
 enddo 
 acc(1)=0.0
end subroutine vel2acc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine vel2dis(nt,dt,acc)
 implicit NONE
 integer, intent(IN) :: nt
 real, intent(IN) :: dt
 real, intent(INOUT), dimension(nt) :: acc
 real :: v1,v2,dt05
 integer :: i
 dt05=0.5*dt;  v1=acc(1);  acc(1)=0.0
 do i=2,nt
   v2=v1;   v1=acc(i)
   acc(i)=acc(i-1)+(v2+v1)*dt05
 enddo 
 return
end subroutine vel2dis
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
function getlen(string)
 implicit NONE
 character(len=*), intent(INOUT):: string
 character(len=256):: chatmp*256
 integer:: nl,getlen

 chatmp=adjustl(string)
 nl=len_trim(chatmp)
 string=''
 string=chatmp(1:nl)
 getlen=nl
end function getlen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine num2ch(r,ch)
 implicit NONE
 real, intent(IN):: r
 character(len=*), intent(OUT):: ch
 integer:: i,n,n2,nl,nc,nzero

 ch='';  n=int(r+0.01)
 nc=3; if(n > 99) nc=3
 nl=10**(nc-1)
 nzero=ichar('0')
 do i=1,nc
   n2=n/nl
   ch(i:i)=char(nzero+n2)
   n=n-n2*nl
   nl=nl/10
 enddo
end subroutine num2ch
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine trav(dis,tmin)
! calculate travel time for horizontal layered model
! it outputs both times for first arrival and direct arrival
!
 use ccmodel
 implicit none
 real, intent(IN):: dis
 real, intent(OUT):: tmin
 real(kind=8):: t, t0, td, x
 complex(kind=8):: pd, p0, p
 integer:: i

 ns=nb-1
 do i=1, nlayer
   th(i)=dble(thick(i))
   v(i)=1.0d0/dble(beta(i))/dble(beta(i))
 enddo

 x=dble(dis)
 lmax = ns
 call find2(x,td,pd)
 t0 = td
 p0 = pd
 do lmax=ns+1, nlayer-1
   call find2(x,t,p)
   if (t < t0) then
     t0=t
     p0=p
   endif
 enddo
 tmin=sngl(t0)
end subroutine trav
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine find2(x,t0,p0)
! solve d(tau)/dp = 0 for p0, tau=p*x+eta*z
! input:  x distance range
! output: p0, t0
!
 use ccmodel
 implicit none
 real(kind=8), intent(IN):: x
 real(kind=8), intent(OUT):: t0 
 complex(kind=8), intent(OUT):: p0
 complex(kind=8):: taup,dtdp,p1,p2
 real(kind=8):: aminm,dtdp0,pm
 integer:: i

 aminm = 1.d-6
 p1 = dcmplx(0.d0,1.0d-20)    ! make p slightly above the real axis
 pm = 9999999.d0
 do i=1,lmax
   pm = dmin1(pm,v(lmax))
 enddo
 p2 = dsqrt(pm) + p1
 p0 = 0.5d0*(p1+p2)          ! the actual p0 lies betw. p1 and p2
 do while( dble(p2-p1) > aminm )
   dtdp0 = dtdp(x,p0)
   if( dabs(dtdp0) < aminm ) exit 
   if( dtdp0 > 0.0d0 ) then
     p1 = p0
   else
     p2 = p0
   endif
   p0 = 0.5d0*(p1+p2)
 enddo
 pm = dmin1(pm,v(lmax+1))
 pm = dsqrt(pm)
 if (lmax>ns .and. pm<dble(p0)) p0 = pm
 t0 = taup(p0,x)
end subroutine find2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function taup(p,x)
! define function tau(p) = p x + eta h
 use ccmodel
 implicit none
 complex(kind=8), intent(IN):: p
 real(kind=8), intent(IN):: x
 complex(kind=8):: pp,taup
 integer:: i

 taup = p*x
 pp = p*p
 do i = 1, ns
   taup=taup+cdsqrt(v(i)-pp)*th(i)
 enddo
 do i = ns+1,lmax
   taup=taup+cdsqrt(v(i)-pp)*th(i)*2.d0
 enddo
end function taup
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function dtdp(x,p)
! define d(tau)/dp
 use ccmodel
 implicit none
 complex(kind=8), intent(IN):: p
 real(kind=8), intent(IN):: x
 complex(kind=8):: pp,dtdp
 integer:: j
 pp = p*p
 dtdp = cmplx(0.0d0,0.0d0)
 do j = 1, ns
   dtdp=dtdp-th(j)/cdsqrt(v(j)-pp)
 enddo
 do j = ns+1, lmax
   dtdp=dtdp-2*th(j)/cdsqrt(v(j)-pp)
 enddo
 dtdp = x + p*dtdp
end function dtdp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module wavelet_base
 implicit NONE
 SAVE
 integer :: num_pt=0, iwt_save=0
 real, allocatable, dimension(:) :: rwt
end module wavelet_base
!
subroutine inv_wt(nlen,u)
! nlen: the size of u. nlen must be a power of 2 (nlen=2**N)
! u: inputting values are wavelet coefficients from fwd_wt
!    outputs are time histories.
! The inputted wavelet coefficients are not multiplied sqrt(dt)
! level(j=0,1, N-1),  input-wavelet-coefficients(u), 
! #-of-wavelet-coefficients 
!     0                     u(2)                               1
!     1                     u(3), u(4)                         2
!     j                     u(2**j+1 : 2**(j+1))               2**j
!
! USES: real_fft.f90
!
 use wavelet_base
 implicit NONE
 integer, intent(IN) :: nlen
 real, intent(INOUT), dimension(nlen):: u 
 real, dimension(2*nlen) :: wcf,fu
 integer :: i1,i2,k1,k2,kmax,kstep
 real :: cff
 real :: cp1,cp2, pi2
 complex :: wave, cc
 if(num_pt /= nlen .or. iwt_save/=1) then
   if(allocated(rwt)) deallocate(rwt)
   allocate(rwt(4*nlen))
   num_pt = nlen  
   iwt_save = 1
   kmax=1
   i2=0
   do while( kmax < nlen)
     do k1=1,2*kmax
       i1=i2+2*k1-1 
       cc=wave(nlen,kmax,k1,1)
       rwt(i1)  =real(cc)
       rwt(i1+1)=aimag(cc)
     enddo
     i2=i2+4*kmax
     kmax=kmax*2
   enddo
 endif
 wcf=0.0
 k1=nlen+1
 k2=k1+1
 i1=nlen/2+1
 i2=i1+1
! j=0
 wcf(1)=u(1)
! j=1, kmax=1
 wcf(3)= u(2)*rwt(3)
 wcf(4)= u(2)*rwt(4)
! j=2, kmax=2
 wcf(3)= wcf(3)+rwt(7)*(u(3)-u(4))
 wcf(4)= wcf(4)+rwt(8)*(u(3)-u(4))
 wcf(5)=        rwt(9)*(u(3)+u(4))
 wcf(6)=        rwt(10)*(u(3)+u(4))
! j=3
 kmax=4
 do 
   do i1=1,kmax
     fu(i1)= u(i1+kmax)
   enddo
   call realft(fu,kmax,-1)
   fu(kmax+1)=fu(2)
   fu(kmax+2)=0.0
   fu(2)=0.0
   k2=2*kmax+2
   do i1=3,kmax,2
     k1=k2-i1
     fu(k1)  = fu(i1)
     fu(k1+1)=-fu(i1+1)
   enddo
   kmax=kmax*2
   kstep=2*kmax-4
   do i2=4,kmax,2
     i1=i2-1
     k1=i1+kstep
     wcf(i1)=wcf(i1)+fu(i1)*rwt(k1)  -fu(i2)*rwt(k1+1)
     wcf(i2)=wcf(i2)+fu(i1)*rwt(k1+1)+fu(i2)*rwt(k1)
   enddo
!
   k2=kmax/2
   kstep=3*kmax-4
   do i2=2,k2,2
     i1=i2-1
     k1=i1+kstep
     wcf(i1+kmax)=fu(i1)*rwt(k1)  -fu(i2)*rwt(k1+1)
     wcf(i2+kmax)=fu(i1)*rwt(k1+1)+fu(i2)*rwt(k1)
   enddo
   if(kmax == nlen) EXIT
 enddo
 call four1(wcf,nlen,1)
 cff=2./float(nlen)
 do i1=1,nlen
   u(i1)=cff*wcf(2*i1-1)
 enddo
 return
end      
!
subroutine fwd_wt(nlen,u)
! nlen: the size of u. nlen must be a power of 2
! u: inputting values are time histories. 
!    outputs are wavelet coefficients at different level
!    (frequencies band).
! The computed wavelet coefficients are not multiplied sqrt(dt)
! level(j=0,1, N-1),  coefficients(u),    #-of-wavelet-coefficients 
!     0                u(2)                               1
!     1                u(3), u(4)                         2
!     j                u(2**j+1 : 2**(j+1))               2**j
!
 use wavelet_base
 implicit NONE
 integer, intent(IN) :: nlen
 real, intent(INOUT), dimension(nlen):: u 
 real, dimension(2*nlen) :: fu
 real, dimension(nlen) :: wcf
 integer :: j,i1,i2,k1,k2,kmax,kstep
 real :: cff
 complex:: wave, cc
 if(num_pt /= nlen .or. iwt_save/=2) then
   if(allocated(rwt)) deallocate(rwt)
   allocate(rwt(4*nlen))
   num_pt = nlen  
   iwt_save = 2
   kmax=1
   i2=0
   do while( kmax < nlen)
     do k1=1,2*kmax
       i1=i2+2*k1-1 
       cc=wave(nlen,kmax,k1,2)
       rwt(i1)  =real(cc)
       rwt(i1+1)=aimag(cc)
     enddo
     i2=i2+4*kmax
     kmax=kmax*2
   enddo
 endif
!
 cff=1./real(nlen)
 do i1=1,nlen
   i2=2*i1
   fu(i2-1) = cff*u(i1)
   fu(i2)   = 0.0
 enddo
 call four1(fu,nlen,-1)
!
! Zero frequency
 u(1)=fu(1)
!j=0, kmax=1
 u(2)=2.*(fu(3)*rwt(3)-fu(4)*rwt(4))
!j=1, kmax=2
 u(3)=2.*(fu(5)*rwt(9)-fu(6)*rwt(10)+fu(3)*rwt(7)-fu(4)*rwt(8))
 u(4)=2.*(fu(5)*rwt(9)-fu(6)*rwt(10)-fu(3)*rwt(7)+fu(4)*rwt(8))
! j=2
 kmax=4
 do while(kmax < nlen)
   kmax=2*kmax
   kstep=2*kmax-4
   do i1=1,kmax,2
     i2=i1+kmax 
     k1=kstep+i1
     k2=kstep+i2
     wcf(i1)  =fu(i1)*rwt(k1)-fu(i1+1)*rwt(k1+1)+ &
               fu(i2)*rwt(k2)-fu(i2+1)*rwt(k2+1)
     wcf(i1+1)=fu(i1)*rwt(k1+1)+fu(i1+1)*rwt(k1)+ &
               fu(i2)*rwt(k2+1)+fu(i2+1)*rwt(k2)
   enddo
   k2=kmax/2
   call four1(wcf,k2,1)
   do k1=1,kmax,2
     k2=k2+1
     u(k2)=2.*wcf(k1)
   enddo 
 enddo
 return
end
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function wave(nlen,kmax,iw,ic)
 implicit NONE
 integer, parameter :: double=kind(0.0d0)
 real(kind=double), parameter :: pi=3.1415926535897d0, &
                   p23=2.0d0*pi/3.0d0,p43=4.0d0*pi/3.0d0
 integer, intent(IN) :: nlen,kmax,iw,ic
 real(kind=double), dimension(2):: g,p
 real(kind=double):: w,hw,ww,wd,wf,wg,wp,fw
 real:: r1,r2
 integer:: i,j,k
 complex:: wave
 w=2.0d0*pi*real(iw-1,double)/real(kmax,double)
 hw=w*0.5d0
 ww=dabs(w)
 wave=cmplx(0.0,0.0)
 if(ww >(2.0d0*p43).or.ww < p23) return
 do i=1,2
   wp=hw*i
   fw=1.0d0
   do j=1,2
     if(j == 1)then
       wf=-wp
     else
       wf=wp
     end if
     do k=1,2
       if(k == 1) wd=p43-wf
       if(k == 2) wd=wf-p23
       if(wd <= 0.0d0)then
         g(k)=0.0d0
       else
        g(k)=exp(-1.0/wd/wd)
       end if
     enddo
     wg=g(1)/(g(1)+g(2))
     fw=fw*wg
   enddo
   p(i)=fw
 enddo
 r1=real(hw)
 r2=real(p(1)-p(2))*float(nlen)/float(kmax)
 wave = cmplx(cos(r1),-sin(r1))*csqrt(cmplx(r2, 0.0))
 if(ic == 1) then
!   w= 2.d0*pi*real(iw-1,double)/real(nlen*dt,double)
!   wave=wave*cmplx(-sngl(w*w),0.0)
 else
   wave = conjg(wave)
 endif
 return
end function wave

