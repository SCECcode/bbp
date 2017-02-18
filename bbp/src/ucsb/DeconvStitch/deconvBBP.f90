MODULE Vmodel_var
 implicit NONE
 save
 integer:: nlayer,nb
 real:: freq_ref
 real, allocatable, dimension(:):: afa,beta,rho,thick,qp,qs,dep
END MODULE Vmodel_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Program deconv
!  
! This code deconvolves  the ground motions at surfac to 
! the given depth
! Written by Pengcheng Liu
! Copyright (c) 2007 by Pengcheng Liu
! modified to J. Schmedes to red 3comp files instead of single components, July 2009
!
 implicit none
 character(len=256):: fileinp,fileRec,fileVM,fileOut,ctmp,junk
 integer:: i,i1,i2,ks,nf,npt,npt2,nfile,id_type_out,comment
 integer:: status,k
 real:: dt,df,factor,v1,v2,depth_out
 real,allocatable, dimension(:):: wave,wave1,wave2,wave3,transf
!
 fileinp='deconvBBP.inp'
 open(unit=8,file=fileinp,status='old')
!
 write(*,*) 'Enter the name of file with 1D velocity Model'
! The deconvolution will be based on this VM
 read(8,'(1a)') fileVM
 write(*,'(1a)') fileVM
!
 write(*,*) 'Enter the depth at where waves to be outputed'
! Note: the unit of depth should be same as that used in 
!       fileVm for thick 
 read(8,*) depth_out
 write(*,*) depth_out
!
 write(*,*)'Enter 1 for outputting Up-going wave' 
 write(*,*)'Enter 2 for outputting Down-going wave' 
 write(*,*)'Enter 3 for outputting Total wave' 
 read(8,*) id_type_out
 write(*,*) id_type_out
! 
 write(*,*) 'Enter the number of data files'
 read(8,*) nfile
 write(*,*) nfile
!
 write(*,*) 'Enter the name of station'
 do ks=1,nfile
   read(8,'(1a)') fileRec
   write(*,*) fileRec(1:len_trim(fileRec))//'.3comp'
   open(unit=13,file=fileRec(1:len_trim(fileRec))//'.3comp',status='old',position='rewind',IOSTAT=status)
   comment=0;	
   read(13,'(1a)') junk
   do while ((junk(1:1).eq.'%').or.(junk(1:1).eq.'#'))
	comment=comment+1;
	 read(13,'(1a)') junk
   end do
   write(*,*) comment	
   npt=1
   status=0
   do while(status.eq.0)
	 read(13,'(1a)',IOSTAT=status) junk
	 npt=npt+1
    end do
    npt=npt-1;
    write(*,*) npt
!
   npt2=1
   do while (npt2 < npt)
     npt2=npt2*2
   enddo
   nf=npt2/2+1
   allocate(wave(npt2),wave1(npt2),wave2(npt2),wave3(npt2),transf(2*nf))
!
   close(13)
   open(unit=13,file=fileRec(1:len_trim(fileRec))//'.3comp',status='old',position='rewind',IOSTAT=status)
   do k=1,comment
	 read(13,'(1a)') junk
   end do
   read(13,*) (wave(i),wave1(i),wave2(i),wave3(i), i=1,npt)
   close(13)
   if (npt2 > npt) then
           do i=npt+1,npt2
                wave(i) = 0.0
                wave1(i) = 0.0
                wave2(i) = 0.0
                wave3(i) = 0.0
           enddo
   endif
   dt=wave(2)-wave(1);
   write(*,*) dt
!
   df=1./(npt2*dt)
   factor=2.0/npt2
   if(ks==1) call cmodel(fileVM,depth_out)
   call linear(nf,df,id_type_out,transf)
!N-S
   call realft(wave1,npt2,-1)
   wave1(1)=wave1(1)*transf(1)
   !write(*,*) wave1(1)
   wave1(2)=wave1(2)*transf(2)
   do i=2,nf-1
     i1=2*i-1
     i2=i1+1
     v1=wave1(i1)
     v2=wave1(i2)
     wave1(i1)=v1*transf(i1)-v2*transf(i2)
     wave1(i2)=v1*transf(i2)+v2*transf(i1)
   enddo
   call realft(wave1,npt2,1)
   wave1(1:npt)=wave1(1:npt)*factor
! Output the result
   fileOut=fileRec(1:len_trim(fileRec))//'.000.de1D.001'; 
   open(unit=14,file=fileOut,status='replace',position='rewind')
   write(14,*) npt,dt
   write(14,'(e14.6)') (wave1(i),i=1,npt)
   close(14)
!E-W
   call realft(wave2,npt2,-1)
   wave2(1)=wave2(1)*transf(1)
   wave2(2)=wave2(2)*transf(2)
   do i=2,nf-1
     i1=2*i-1
     i2=i1+1
     v1=wave2(i1)
     v2=wave2(i2)
     wave2(i1)=v1*transf(i1)-v2*transf(i2)
     wave2(i2)=v1*transf(i2)+v2*transf(i1)
   enddo
   call realft(wave2,npt2,1)
   wave2(1:npt)=wave2(1:npt)*factor
! Output the result
   fileOut=fileRec(1:len_trim(fileRec))//'.090.de1D.001'; 
   open(unit=14,file=fileOut,status='replace',position='rewind')
   write(14,*) npt,dt
   write(14,'(e14.6)') (wave2(i),i=1,npt)
   close(14)
!U-D
   call realft(wave3,npt2,-1)
   wave3(1)=wave3(1)*transf(1)
   wave3(2)=wave3(2)*transf(2)
   do i=2,nf-1
     i1=2*i-1
     i2=i1+1
     v1=wave3(i1)
     v2=wave3(i2)
     wave3(i1)=v1*transf(i1)-v2*transf(i2)
     wave3(i2)=v1*transf(i2)+v2*transf(i1)
   enddo
   call realft(wave3,npt2,1)
   wave3(1:npt)=wave3(1:npt)*factor
! Output the result
   fileOut=fileRec(1:len_trim(fileRec))//'.ver.de1D.001'; 
   open(unit=14,file=fileOut,status='replace',position='rewind')
   write(14,*) npt,dt
   write(14,'(e14.6)') (wave3(i),i=1,npt)
   close(14)
   deallocate(wave,wave1,wave2,wave3, transf)
 enddo
 close(8)
 stop
end program deconv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine cmodel(fileVM,depth)
 use Vmodel_var
 implicit none
 character(len=*), intent(IN):: fileVM
 real, intent(IN):: depth
 integer:: jo,j1,j
!
 open(unit=12,file=fileVM,status='old',position='rewind')
 read(12,*) jo, freq_ref
 if(allocated(beta)) deallocate(afa,beta,rho,thick,qp,qs,dep)
 j1=jo+1
 allocate(afa(j1),beta(j1),rho(j1),thick(j1),qp(j1),qs(j1),dep(j1+1))
 do j=1,jo
   read (12,*) afa(j),beta(j),rho(j),thick(j),qp(j),qs(j)
 enddo 
 close(12)
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
subroutine linear(nf,df,id_type_out,transf)
 use Vmodel_var
 implicit NONE
 integer, intent(IN):: nf,id_type_out
 real, intent(IN):: df
 real, dimension(2*nf), intent(OUT):: transf
!
 complex:: wki,xlnf,w_impact,ctmp,enkh,epkh
 complex, allocatable, dimension(:):: EE,FF,cbeta
 integer:: id_ly_in,id_ly_out,i,j
 real:: e_i,f_i,e_o,f_o,pi,pi_2,dw,wr,ww,gamb

 allocate(EE(nlayer),FF(nlayer),cbeta(nlayer))

 id_ly_in=1
 e_i=1.0
 f_i=1.0
 id_ly_out=nb
 e_o=1.0
 if(id_type_out.eq.2) e_o=0.0
 f_o=1.0
 if(id_type_out.eq.1) f_o=0.0

 pi=4.*atan(1.0)
 pi_2=0.5*pi
 dw=2.*pi*df
 wr=freq_ref*pi*2.0
 do j=1,nf
   ww=(j-1.)*dw
   do i=1,nlayer
! Kjartansson, E. (1979). Constant Q-wave propagation 
! and attenuation, J. Geophys. Res. 84, 4737-4748.
     gamb= atan(1./qs(i))/pi
     xlnf=cmplx(0.0,(dw+ww)/wr)
     cbeta(i) =beta(i)/cos(pi_2*gamb)*xlnf**gamb
   enddo 
!
   EE(1) = (1., 0.)
   FF(1) = (1., 0.)
   do i=1,nlayer-1
     w_impact =0.5*(rho(i)*cbeta(i))/(rho(i+1)*cbeta(i+1))
     wki=cmplx(0.0, ww*thick(i))/cbeta(i)
     EE(i+1)=EE(i)*(0.5+w_impact)*cexp( wki)+ &
             FF(i)*(0.5-w_impact)*cexp(-wki)
     FF(i+1)=EE(i)*(0.5-w_impact)*cexp( wki)+ &
             FF(i)*(0.5+w_impact)*cexp(-wki)
   enddo
   ctmp=(e_o*EE(id_ly_out)+f_o*FF(id_ly_out)) / &
        (e_i*EE(id_ly_in) +f_i*FF(id_ly_in) )
   transf(2*j-1) =real(ctmp)
   transf(2*j  ) =aimag(ctmp)
 end do
 transf(2)=transf(2*nf-1)
 deallocate(EE,FF,cbeta)
end subroutine linear
