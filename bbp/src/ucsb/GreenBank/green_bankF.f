c       Make the Green Function bank
c       Input name
c       Input dep_max dep_min dep_step
c       Input dist_max dist_min d_step
c       Input numdt dt,t_cor
c       Input Green_name
c       output Green

	include "model.h"

c	parameter(block_g=nt*32+100)
c	dimension green(ndis,nt,8),dist(ndis)
        real c(nlay), s(nlay),den(nlay),th(nlay),qa(nlay),qb(nlay)
	real dist(ndis),t0(ndis),tmin(ndis),t_cor,depth,depthMinUsed
	integer block_gg,jo,nb,numdt,ncom_out
	character*60 name, output,junk
        complex sum(ndis,9,nt)
        common/model_new/jo,nb,c,s,den,th,qa,qb,ncom_out

	open(1,file="GreenFar.in")
	
	write(*,*)"Input name of velocity model"
	read(1,1) name

	write(*,*)"input Max., Min., and step of epicenter depth"
	read(1,*) dep_max, dep_min, dep_step

	write(*,*)"input Max., Min., and step of epicenter distance"
	read(1,*) dist_max,dist_min,d_step

	write(*,*)"input number of time-point, time-step, and seconds",
     &            " to save before first-arrival"
	read(1,*) numdt,dt,t_cor   

	write(*,*)"Input the outputfile"
	read(1,1) output
	write(*,*) output

	write(*,*)"Input the minimal depth actually used"
	read(1,1) depthMinUsed
1	format(a)

	lnpt=0
	npt=1
101	lnpt=lnpt+1
	npt=npt*2
	if(npt .lt. numdt) goto 101

	if(npt.gt.nt)then
	write(*,*)  "Please increase <nt> in model.h to ",npt
	stop
	endif
	if(dep_max.lt.dep_min)then
	   pause "The depth region is wrong"
	endif
	if(dist_max.lt.dist_min)then
	   pause "distance region is wrong"
	endif
	if(dist_min.lt.100.)then
	   pause "careful, only for stations further out if shallow sources should be used"
	endif

	nx=int((dist_max-dist_min)/d_step)+1
	nz=int((dep_max-dep_min)/dep_step)+1
	block_gg=4*(npt*8+7)+12

	write(*,*) 'total nx and nz =',nx,nz
	if(nx.gt.ndis) then
	  write(*,*)  "Please increase <ndis> in model.h to ",nx
	  stop
	endif

	open(11,file=output,status='unknown',access='direct',
     1	        recl=block_gg)
	nhfpt=npt/2
	ll=0
	do iz=1,nz
	   depth=dep_min + (iz-1)*dep_step
C inclded this line, use a different depth is too shallow, only for stations far out
	   if(depth.lt.depthMinUsed) then
		depth=depthMinUsed
	   endif
	   do k=1,nx
	      dist(k)=dist_min+(k-1)*d_step
	   enddo
	   io_model=30
	   open(unit=io_model,file=name,status='old')
	   call cmodel(io_model,depth)
	   call trav(dist,nx,tmin)

	   do k=1,nx
	      t0(k)=tmin(k)-t_cor
	   enddo

	   call sub_bs_dc(nx,dist,dist_max,lnpt,dt,t0,sum)
C but write actually the depth that should have been used, otherwise there will be an error in syn1d
	   depth=dep_min + (iz-1)*dep_step
	   do k=1,nx
	      ll=ll+1
	      write(11,rec=ll) iz,k,dist(k),t0(k),depth,dt,npt,
     1	        ((real(sum(k,n_com,ntc)),aimag(sum(k,n_com,ntc)),
     2            ntc=1,nhfpt),n_com=1,8)
	   enddo
	   close(30)

	  write(*,*) 'depth = ',iz,depth
	   if(iz<10) then
		   write(junk,'(I1)') iz
	 endif
	   if(iz>=10) then
		   write(junk,'(I2)') iz
	 endif
	   open(12,file='test'//junk)
	   write(12,*)1
	   write(12,*)"(7e11.4)"
	   do n_com=1,8
	      write(12,*)dist(1),t0(1)
	      write(12,*)npt,dt
	      write(12,*)((real(sum(1,n_com,ntc)),aimag(sum(1,n_com,ntc))),ntc=1,nhfpt)
        
	   enddo	
	close(12)
	enddo
	close(11)
	close(10)

	open(19,file='Green_Bank.inf',status='unknown')
	write(19,*) "Name of velocity model"
	write(19,'(a)') name
	write(19,*) "Maximum, Minimum, and Step of Epicenter Depth"
	write(19,*) dep_max, dep_min, dep_step
	write(19,*)"Maximum, Minimum, and Step of Epicenter Distance"
	write(19,*) dist_max,dist_min,d_step
	write(19,*)"The lnpt, dt, block_gg, t_cor"
	write(19,*) lnpt,dt,block_gg,t_cor
	write(19,*)"The name of file to store Green Bank"
	write(19,'(a)') output
	write(19,*)"The number of GreenFunctions in Distance and Depth"
	write(19,*) nx,nz
	write(19,*) "Velocity Structure model"
	write(19,*) jo,nb
	do jj=1,jo
	  write(19,*) th(jj),c(jj), s(jj),den(jj),qa(jj),qb(jj)
	enddo
	close(19)

	stop
	end	

