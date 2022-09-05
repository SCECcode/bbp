Program single2multiple
!use sp_sub_f
use source_params_all
implicit NONE
!
! By the definition of Liu et al. (2006),
!  xij=(i-0.5)*dx-str_ref
!  yij=(j-0.5)*dy-dip_ref
!
! Then, the single fault plane is defined as
! 1. includes all fault segments
! 2. inlcudes the fault segment with hypocenter
!

character(len=72):: vsfile,spfile,spfile_best,spfile_best_m
integer, dimension(4):: nb_taper_TRBL
integer:: getlen,idum1,idum2,idum3,id_ran1,id_ran2,&
      nsource,i,j,k,nlch
real:: sp2,drp,rdip,rstrike,cosx,sinx,str_ref,dip_ref,area_sub, &
       summ,rmt,rvc,xij,yij,xs,ys,xps,yps,zps,factor,rakei,dipi,&
       pamzi,pdip_max,prake_max
real::xstk,rise_time,x_hypc_n, y_hypc_n,Moment_o
integer::io_model_out,id_min,nfs,id

open(8,file='ffsp.inp',status='old')
!write(7,*) 'Enter the index of slip rate function'
read(8,*)  id_sf_type, freq_min, freq_max

!write(7,*) ' Enter the length and width of main fault'
read(8,*)  flx,fwy

!write(7,*) ' Enter the distance in the strike and down dip'
!write(7,*) ' direction, and the depth of the hypocenter'
read(8,*)  x_hypc, y_hypc, depth_hypc

!write(7,*) 'Enter the coordinates  of hypocenter'
read(8,*)  xref_hypc,yref_hypc

!write(7,*) 'Enter the moment(N-m), corner frequency(Hz) and ', &
!           'Average rupture speed'

read(8,*)  Moment_o,fc_main_1,fc_main_2,rv_avg

!write(7,*) 'Enter the ratio of time between slip-rate increase', &
!           ' and rise time'

read(8,*)  cft1
!write(7,*) cft1

!write(7,*) ' Enter the strike, dip and rake of main fault'
read(8,*)  strike,dip,rake

!write(7,*) 'Enter Max. perturbation on the dip and rake'
read(8,*)  pdip_max,prake_max

!write(7,*) ' Enter the subfaults number along strike adn dip'
read(8,*)  nsubx,nsuby

!write(7,*) ' Enter the numbers of subfault to be tapered at each side'
!           Top-Right-Bottom-Left
read(8,*)  (nb_taper_TRBL(i), i=1,4)

!write(7,*) ' Enter the seed for generating random numbers'
read(8,*)  idum1,idum2,idum3

!write(7,*)'Enter first and last index of random generations'
read(8,*)  id_ran1,id_ran2

!write(7,*) ' Enter the name of file with velocity', &
!           ' structure of source region'
read(8,'(1a)')  vsfile
!write(7,'(1a)') vsfile
!
!  The Parameters for outputing results
!
!write(7,*) 'Enter the angle (degree) from North to X-Axis (clockwise) of SYN. code'
read(8,*)  angle_north_to_x

if(abs(angle_north_to_x).gt.1.0e-6)THEN
  write(*,*)"In the corrent version, this must be fixed to zero"
  stop
endif

!write(7,*) 'Do you output moment(1), slip-area (2), or slip(3) '
read(8,*)  is_moment

!write(7,*) ' Enter the file name of output source parameters'
read(8,'(1a)')  spfile
! read in the information of multiple segments

read(8,*)n_fault_segment
nfs=n_fault_segment
allocate(x_stk_min(nfs),x_stk_max(nfs),y_dip_min(nfs),y_dip_max(nfs), &
        seg_stk(nfs),seg_dip(nfs),seg_rake(nfs),seg_upleft_n(nfs),&
        seg_upleft_e(nfs),seg_upleft_z(nfs),io_reverse(nfs))

do i=1,n_fault_segment
  read(8,*)id_segment
  read(8,*)x_stk_min(i),x_stk_max(i)
  if(x_stk_min(i).gt.x_stk_max(i))THEN
    io_reverse(i)=1
    xstk=x_stk_min(i)
    x_stk_min(i)=x_stk_max(i)
    x_stk_max(i)=xstk
  else
    io_reverse(i)=0
  endif
  read(8,*)y_dip_min(i),y_dip_max(i)
  read(8,*)seg_stk(i),seg_dip(i),seg_rake(i)
! up-left position of fault segment relative to epicenter
  read(8,*)seg_upleft_n(i),seg_upleft_e(i),seg_upleft_z(i)
  write(*,*)seg_upleft_n(i),seg_upleft_e(i),seg_upleft_z(i)
enddo
close(8)
!
open(11,file='source_model.list',status='old')
read(11,118) id, nsubx, nsuby, dx, dy, x_hypc_n, y_hypc_n
read(11,119) xref_hypc, yref_hypc, angle_north_to_x
read(11,'(1a)') spfile_best
close(11)
nlch=len_trim(spfile_best)+1
spfile_best_m=spfile_best
spfile_best_m(nlch:nlch+1)="_m"
open(11,file='source_model.list_m')
write(11,118) id, nsubx, nsuby, dx, dy, x_hypc_n, y_hypc_n
write(11,119) xref_hypc, yref_hypc, angle_north_to_x
write(11,'(1a)') spfile_best_m
close(11)
118     format(I3, 2I5, 4e15.7)
119     format(3e15.7)

call input_source_param(spfile_best,spfile_best_m)

do  i=1,nsubx
do  j=1,nsuby
  k=(i-1)*nsuby+j
  xij=(i-0.5)*dx-str_ref
  yij=(j-0.5)*dy-dip_ref
  xs=xij*cos(rstrike)-yij*sin(rstrike)*cos(rdip)
  ys=xij*sin(rstrike)+yij*cos(rstrike)*cos(rdip)
  xps=xref_hypc  + xs*cosx+ys*sinx
  yps=yref_hypc  - xs*sinx+ys*cosx
  zps=depth_hypc + yij*sin(rdip)
end do
end do

end program single2multiple
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine input_source_param(filesour,fileout)
! Read the source parameters
 use source_params_all
 implicit NONE
 character(len=72):: filesour,fileout
 !integer, intent(IN):: ncx,ncy
 integer:: i,j,nfs,is,ks,nx_ref
 real:: dr,rstrike,rdip,xij,yij,cosx,sinx
 real:: drp,str_ref,dip_ref,xps,zps,yps,xs,ys
!
 write(*,1)filesour
 write(*,1)fileout
 !
 1 format(a)
 open(unit=19, file=filesour,status='old')
 open(unit=29, file=fileout,status='unknown')
 read(19,*) is_moment,npsour,id_sf_type,cft1
 write(29,*) is_moment,npsour,id_sf_type,cft1
!
!nps_sub=nsubx*nsuby

 if(allocated(ps_xn)) then
    deallocate(ps_xn,ps_ye,ps_dz,ps_mt,ps_rp, &
               ps_rs,ps_p2,ps_sk,ps_dp,ps_rk,drtx,drty, &
               xps_sub,yps_sub,zps_sub )
    deallocate(nx_stk_max,nx_stk_min,ny_dip_max,ny_dip_min)
 endif
 allocate(ps_xn(npsour),ps_ye(npsour),ps_dz(npsour),ps_mt(npsour), &
          ps_rp(npsour),ps_rs(npsour),ps_p2(npsour),ps_sk(npsour), &
          ps_dp(npsour),ps_rk(npsour),drtx(npsour),drty(npsour) )
 allocate(xps_sub(nps_sub),yps_sub(nps_sub),zps_sub(nps_sub))
 allocate(nx_stk_max(nps_sub),nx_stk_min(nps_sub),ny_dip_max(nps_sub), &
          ny_dip_min(nps_sub))
!
 write(*,*)npsour

 do ks=1,npsour
   read(19,*) ps_xn(ks),ps_ye(ks),ps_dz(ks),ps_mt(ks),ps_rp(ks), &
              ps_rs(ks),ps_p2(ks),ps_sk(ks),ps_dp(ks),ps_rk(ks)
 enddo
 close(19)
!
! The above do loop over subfaults (or point sources)
!
! 1-3.  are the coornidates in meters (North, East, Down) of each point source.
! 4-7.  are the moment (N-m), rupture time (s), rise time(s), and another source
!       parameters (any value if not used)
!
! 8-10. are the strike, dip, and rake in degree
!

 do ks=1,npsour
   ps_xn(ks)=ps_xn(ks)/1000.0
   ps_ye(ks)=ps_ye(ks)/1000.0
   ps_dz(ks)=ps_dz(ks)/1000.0
 enddo
 drp=acos(-1.0)/180.0
 cosx=cos(angle_north_to_x*drp)
 sinx=sin(angle_north_to_x*drp)
! str_ref=(ncx-0.5)*dx
! dip_ref=(ncy-0.5)*dy
 write(*,*)npsour,cosx,sinx,drp
 dx=dx/1000.0
 dy=dy/1000.0
 do is=1,n_fault_segment
   nx_stk_min(is)=int(x_stk_min(is)/dx-0.5)+1

   nx_stk_max(is)=int(x_stk_max(is)/dx-0.5)+1

   if(io_reverse(is).gt.0)then
     nx_ref=nx_stk_min(is)
   else
     nx_ref=nx_stk_max(is)
   endif
   ny_dip_min(is)=int(y_dip_min(is)/dy-0.5)+1
   ny_dip_max(is)=int(y_dip_max(is)/dy-0.5)+1

   rstrike=seg_stk(is)*drp
   rdip=seg_dip(is)*drp
   cosx=cos(angle_north_to_x*drp)
   sinx=sin(angle_north_to_x*drp)
   write(*,*)x_stk_min(is),x_stk_max(is),is,dx,dy
   write(*,*)"****",nx_stk_max(is),nx_stk_min(is),ny_dip_min(is),ny_dip_max(is)
   write(*,*)seg_upleft_z(is),io_reverse(is)
   write(*,*)seg_rake(is),rake
   write(*,*)seg_dip(is),dip
   write(*,*)seg_stk(is),strike
   do i=1,nx_stk_max(is)-nx_stk_min(is)+1
     do j=1,ny_dip_max(is)-ny_dip_min(is)+1
       if(io_reverse(is).eq.0)then
         ks=(i+nx_stk_min(is)-2)*nsuby+j+ny_dip_min(is)-1
       else
         ks=(nx_stk_max(is)-i)*nsuby+j+ny_dip_min(is)-1
       endif
       xij=(i-0.5)*dx
       yij=(j-0.5)*dy
       xs=xij*cos(rstrike)-yij*sin(rstrike)*cos(rdip)
       ys=xij*sin(rstrike)+yij*cos(rstrike)*cos(rdip)
       xps=seg_upleft_n(is) + xs*cosx+ys*sinx
       yps=seg_upleft_e(is) - xs*sinx+ys*cosx
       zps=seg_upleft_z(is) + yij*sin(rdip)
       if(i.eq.1.and.j.eq.1)write(*,*)i,j,is,ks,zps
       ps_rk(ks)=ps_rk(ks)+(seg_rake(is)-rake)
       ps_dp(ks)=ps_dp(ks)+(seg_dip(is)-dip)
       ps_sk(ks)=ps_sk(ks)+(seg_stk(is)-strike)
       ps_xn(ks)=xps*1000.0
       ps_ye(ks)=yps*1000.0
       ps_dz(ks)=zps*1000.0
!       write(29,109)xps*1000.,yps*1000.,zps*1000.,ps_mt(ks),ps_rp(ks), &
!                  ps_rs(ks),ps_p2(ks),ps_sk(ks),ps_dp(ks),ps_rk(ks)
     end do
   end do
 end do
 do ks=1,npsour
   write(29,109) ps_xn(ks),ps_ye(ks),ps_dz(ks),ps_mt(ks),ps_rp(ks), &
              ps_rs(ks),ps_p2(ks),ps_sk(ks),ps_dp(ks),ps_rk(ks)
 enddo
 close(29)
 109     format(4e15.7,6e12.4)
end subroutine input_source_param
