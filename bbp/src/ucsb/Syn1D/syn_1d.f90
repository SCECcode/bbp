Program syn_1d
!
! This code uses the green functions from green_bank.f and the
! source functions from ffsp_v2.f to calculate the grond motions
! at free surface.
!
! The input parameter from all the files, except 'Green_bank.inf',
! should has the unit of  N-m, m, m/s, kg/m^3,  and/or degree.
!
! X-Y-Z Cartesian coordinate system follows righ-hand rule with Z
! direction positively down
!
! The Output ground motions have positive direction in UP, H1 or H2
! (H1 and H2 are angles (degree) reading in from file 'fileRec').
! The unit of outputing ground motions will be
! m (displacement), m/s (velocity), or m/s/s (acceleration).
!
! Written by Pengcheng Liu
! Copyright (c) 2005 by Pengcheng Liu
!
!======================================================================
 use source_params
 use station_params
 implicit NONE
 character(len=72):: fileinp,fsou_list,filesour,fileRec,file_gf_info
 character(len=72):: st_name,fout1,fout2,fout3,ch1,ch2,ch3
 character(len=5), dimension(3):: chmotion
 character(len=6), dimension(4):: chformat
 integer:: npsx,npsy,getlen,id_motion,id_format,nsrc_model,ism,ks,lc
 real:: tdura,azimu_p,rake_p,dip_p,freq_b,freq_e,rdp,kap
!
 chmotion=(/ 'Displ ','Veloc ','Accel ' /)
 chformat=(/ 'SAC   ','TXT   ','FXDR  ','Binary' /)
!
 fileinp='syn_1d.inp'
 open(unit=8,file=fileinp,status='old', position='rewind')

! write(*,*) 'Enter point source number in strike and dip'
! write(*,*) 'direction for each subfault'
 read(8,*) npsx,npsy

! write(*,*) 'Enter the perturbations for Azimuth, Rake, Dip'
 read(8,*)  azimu_p,rake_p,dip_p

! write(*,*) 'Enter the frequency comtrol points for perturb'
 read(8,*) freq_b, freq_e

! write(*,*) 'Enter the duration of outputing GM'
 read(8,*) kap, tdura

! write(*,*) 'Enter the name of file list source files'
 read(8, '(1a)') fsou_list

! write(*,*) 'Enter the name of file with receiver locations'
 read(8, '(1a)') fileRec

! write(*,*) 'Outptut Displacement (=1), Vel. (2), or Acc (3)'
 read(8,*) id_motion
! write(*,'(1a)') chmotion(id_motion)

! write(*,*) 'Enter:  1 for SAC;  2 TXT; 3 FXDR; 4 Binary'
 read(8,*) id_format
! write(*,'(1a)') chformat(id_format)

 close(8)
!============================================================
!
 open(21,file=fsou_list,status='old', position='rewind')
 read(21,*) nsrc_model,nsubx,nsuby,dxsub,dysub,cxp,cyp
 read(21,*) xref_hypc,yref_hypc,angle_north_x
!
! xref_hypc,yref_hypc are the (x, y) value of hypocenter in the
! coordinate system used for describing  the source parameters.
! angle_north_x is the angle from North to the X-axis
!
!--------------------------------------------------------------
!
!   Input parameters for each station
!
 call input_st_parm(fileRec,xref_hypc,yref_hypc,angle_north_x)
!
!--------------------------------------------------------------
!  Input parameters for Greens function
!
 call Read_Green_Bank_inf(azimu_p,rake_p,dip_p,freq_b,freq_e,tdura)
!
!--------------------------------------------------------------
!
 rdp=atan(1.0)/45.0
 dxsub=dxsub/1000.
 dysub=dysub/1000.
 npsx=2*(npsx/2)+1
 npsy=2*(npsy/2)+1
 angle_north_x= angle_north_x*rdp
!--------------------------------------------------------------
!
! Loop over Source Model
!
 do ism=1,nsrc_model
   read(21,'(1a)') filesour
   call input_source_param(filesour,npsx,npsy)
   do ks=1,nst
     st_x_n=stan_x(ks)
     st_y_e=stan_y(ks)
     st_z_up=stan_z(ks)
     ang_h1=stan_h1(ks)*rdp
     ang_h2=stan_h2(ks)*rdp
     call gf_bound(st_x_n,st_y_e,ss_min,ss_max)
     call synthe(id_motion,st_x_n,st_y_e,ang_h1,ang_h2,ss_min,ss_max)
     st_name =adjustl(stan_name(ks))
     lc=len_trim(st_name)+1
     st_name(lc:lc)='.'
     fout1=''
     fout2=''
     fout3=''
     call num2ch(stan_h1(ks),ch1)
     call num2ch(float(ism),ch3)
     call num2ch(stan_h2(ks),ch2)
     fout1(1:lc+12)=st_name(1:lc)//'ver'//'.gm1D.'//ch3(1:3)
     fout2(1:lc+12)=st_name(1:lc)//ch1(1:3)//'.gm1D.'//ch3(1:3)
     fout3(1:lc+12)=st_name(1:lc)//ch2(1:3)//'.gm1D.'//ch3(1:3)
     if(kap.gt.0.0) then
        call kappa(npt,dt,syn_up,kap) !applying kappa to each seismogram
        call kappa(npt,dt,syn_h1,kap)
        call kappa(npt,dt,syn_h2,kap)
     endif
     call output3(fout1,id_format,npt,dt,syn_up)
     call output3(fout2,id_format,npt,dt,syn_h1)
     call output3(fout3,id_format,npt,dt,syn_h2)
   enddo  ! end loop of stations (ks)
 enddo  ! end loop of source model (ism)
 close(21)
 stop
end program syn_1d
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine synthe(id_motion,st_x,st_y,ang_h1,ang_h2,ss_min,ss_max)
 use source_params
 use gf_comm
 implicit NONE
 integer, intent(IN):: id_motion
 real, intent(IN):: st_x,st_y,ang_h1,ang_h2,ss_min,ss_max
 integer:: nfre,ndur,i_dip,i_str,ips,js,jf,it,k,j,iseed,i0,j0,ii,jj,ij
 integer, save:: idum=-5
 real:: zz_min,zz_max,xp0,yp0,zp0,paz,rake0,dip0,theta0,rptm0,drx0,dry0
 real:: moment,rise,spm2,x_ps,y_ps,depth_ps,rptm,dis,az,baz,cmu,df
 real:: ran1,gasdev

!
 i0=int(ran1(idum)*nsubx)*0 -1
 j0=int(ran1(idum)*nsuby)*0 -1
 iseed=-3
 df=1./(dt*npt)
 nfre=npt/2
 zg_first=dep_max+1.e+5
 zg_last =-1.e+5
 syn_up=0.0;  syn_h1=0.0; syn_h2=0.0
 do i_dip=1,nsuby
   zz_min=depth_min(i_dip)
   zz_max=depth_max(i_dip)
   call read_green(npt,dt,ss_min,ss_max,zz_min,zz_max)
   do i_str=1,nsubx
     ii=mod(i_str+i0, nsubx)+1
     jj=mod(i_dip+j0, nsuby)+1
     ij=(ii-1)*nsuby+jj
     paz=ps_p2(ij)

     js=(i_str-1)*nsuby+i_dip
     xp0=ps_xn(js)
     yp0=ps_ye(js)
     zp0=ps_dz(js)
     rake0 = ps_rk(js)
     dip0  = ps_dp(js)
     theta0= ps_sk(js)
     rptm0 = ps_rp(js)
     drx0  = drtx(js)
     dry0  = drty(js)
     gfun=0.0
     do ips=1,nps_sub
       x_ps=xp0+xps_sub(ips)
       y_ps=yp0+yps_sub(ips)
       depth_ps=zp0+zps_sub(ips)
       rptm=rptm0+drx0*sh_strk(ips)+dry0*sh_dip(ips)
       call distaz(st_x,st_y,x_ps,y_ps,angle_north_x,dis,az,baz)
       call rad_pattern(iseed,az,paz,rake0,dip0,theta0,ang_h1,ang_h2)
       call inter_green(dis,depth_ps,rptm,dt,npt,green_ps)
!  H1
       do k=1,5
       do jf=1,npt
         gfun(jf,1)=gfun(jf,1)+rad_h1(jf,k)*green_ps(jf,k)
       enddo
       enddo
!  H2
       do k=1,5
       do jf=1,npt
         gfun(jf,2)=gfun(jf,2)+rad_h2(jf,k)*green_ps(jf,k)
       enddo
       enddo
! VT
       do k=6,8
       do jf=1,npt
         gfun(jf,3)=gfun(jf,3)+rad_vt(jf,k-5)*green_ps(jf,k)
       enddo
       enddo
     enddo
!
! seiemic moment= shear_modulu * area * slip
! 1.e-15 is for correction from (m/s)^2*(kg/m^3)*m^2
! to (km/s)^2*(g/cm^3)*km^2
     moment=ps_mt(js)/float(nps_sub)*1.0e-15
     if(is_moment /= 1) then
       moment=moment*cmu(depth_ps)
     endif
     rise = ps_rs(js)
     spm2=cft1
     if(id_sf_type == 8 ) spm2 = ps_p2(js)
     if(id_sf_type == 6 ) spm2 = ps_p2(js)
     !write(*,*)'ID is: ',id_sf_type
     call source_fun(id_sf_type,ndur,npt,dt,moment,rise,spm2,fsour)
     call conv_sf_gf(npt,ndur,syn_tmp,fsour,gfun(:,1))
     do it=1,npt
       syn_h1(it)=syn_h1(it)+syn_tmp(it)
     enddo
     call conv_sf_gf(npt,ndur,syn_tmp,fsour,gfun(:,2))
     do it=1,npt
       syn_h2(it)=syn_h2(it)+syn_tmp(it)
     enddo
     call conv_sf_gf(npt,ndur,syn_tmp,fsour,gfun(:,3))
     do it=1,npt
       syn_up(it)=syn_up(it)+syn_tmp(it)
     enddo
   enddo
 enddo
 call vel_2_acc_or_disp(id_motion,npt,dt,syn_up)
 call vel_2_acc_or_disp(id_motion,npt,dt,syn_h1)
 call vel_2_acc_or_disp(id_motion,npt,dt,syn_h2)
end subroutine synthe
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine input_st_parm(fileRec,xref_hypc,yref_hypc,angle_north_x)
 use station_params
 implicit none
 character(len=72):: fileRec
 real, intent(IN):: xref_hypc,yref_hypc,angle_north_x
 integer:: ks
 real:: dra,xs,ys,cs,sn
!
 open(17,file=fileRec,status='old',position='rewind')
 read(17,*) nst,xs0,ys0,afa
!
! xs0,ys0 are the (x, y) value of hypocenter in the
! coordinate system used for describing the location of stations
! afa is the angle from North to the X-axis
!
 if(allocated(stan_name)) then
   deallocate(stan_name,stan_x,stan_y,stan_z,stan_h1,stan_h2)
 endif
 allocate(stan_name(nst),stan_x(nst),stan_y(nst),stan_z(nst), &
          stan_h1(nst),stan_h2(nst))

 dra=atan(1.0)/45.0*(afa-angle_north_x)
 cs=cos(dra)
 sn=sin(dra)
 do ks=1,nst
   read(17,'(1a)') stan_name(ks)
   read(17,*) stan_x(ks),stan_y(ks),stan_z(ks),stan_h1(ks),stan_h2(ks)
   xs=(stan_x(ks)-xs0)/1000.0
   ys=(stan_y(ks)-ys0)/1000.0
   stan_x(ks)=xs*cs-ys*sn+xref_hypc/1000.0
   stan_y(ks)=ys*cs+xs*sn+yref_hypc/1000.0
   stan_z(ks)=stan_z(ks)/1000.0
 enddo
 close(17)
end subroutine input_st_parm
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine Read_Green_Bank_inf(azimu_p,rake_p,dip_p,freq_b,freq_e,tdur)
 use gf_comm
 use source_params
 implicit NONE
 real, intent(IN):: azimu_p,rake_p,dip_p,freq_b,freq_e,tdur
!
 character(len=72):: file_model,temp_char
 integer:: jj,i
 real:: dr,df,ff,varf
!
 open(19,file='Green_Bank.inf',status='old',position='rewind')
 read(19,'(a)') temp_char
 read(19,'(a)') file_model
 read(19,'(a)') temp_char
 read(19,*) dep_max, dep_min, dep_step
 read(19,'(a)') temp_char
 read(19,*) dist_max,dist_min,d_step
 read(19,'(a)') temp_char
 read(19,*) lnpt,dt,block_gg,t_cor
 read(19,'(a)') temp_char
 read(19,'(a)') file_bank
 read(19,'(a)') temp_char
 read(19,*) nx_in,nz_in
 read(19,'(a)') temp_char
 read(19,*) jo,nb
!
 nlay=jo+2
 if(allocated(vp)) deallocate(th,vp,vs,den,qa,qb,depth_min,depth_max)
 allocate(th(nlay),vp(nlay),vs(nlay),den(nlay),qa(nlay),qb(nlay))
 allocate(depth_min(nsuby),depth_max(nsuby))
!
!
!  unit: km, km/s, km, g/cm^3
 do jj=1,jo
   read(19,*) th(jj),vp(jj),vs(jj),den(jj),qa(jj),qb(jj)
 enddo
 close(19)
! Checking
 if(dep_max < dep_min)then
    write(*,*) "The depth region is wrong"
    stop
 endif
 if(dist_max < dist_min)then
    write(*,*)  "distance region is wrong"
    stop
 endif
!
 npt=2**lnpt
 do while((npt-1)*dt < tdur)
   npt=npt*2
 enddo
! Allocating memory
 if(allocated(strk_bud)) then
   deallocate(strk_bud,dip_bud,rake_bud,rd_wl)
   deallocate(syn_up,syn_h1,syn_h2,fsour,syn_tmp)
   deallocate(gfun,green_ps,rad_h1,rad_h2,rad_vt)
 endif
 allocate(strk_bud(npt),dip_bud(npt),rake_bud(npt),rd_wl(npt))
 allocate(syn_up(npt),syn_h1(npt),syn_h2(npt),fsour(npt),syn_tmp(npt))
 allocate(gfun(npt,3),green_ps(npt,8),rad_h1(npt,5),rad_h2(npt,5),rad_vt(npt,3))
!
 dr=atan(1.0)/45.0
 df=1./(dt*npt)
 dw=df*atan(1.0)*8.0
 kfb=ifix(abs(freq_b)/df)-1
 kfe=ifix(abs(freq_e)/df)+2
!
 do i=1,npt/2
  ff=(i-1)*df
  varf=amax1(0.0, amin1(2., 2.*(ff-freq_b)/(freq_e-freq_b)))
  strk_bud(i)=varf*0.5*azimu_p*dr
  dip_bud(i) =varf*dip_p*dr
  rake_bud(i)=varf*rake_p*dr
  rd_wl(i)   =0.25
  if(ff < 1.0) rd_wl(i)=ff*0.25
 enddo
end subroutine Read_Green_Bank_inf
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function cmu(depth)
 use gf_comm
 implicit NONE
 real, intent(IN):: depth
 integer:: ih,ii
 real:: cmu,dzsum
!
 ii=1
 dzsum=th(1)
 do ih=2,jo
   if(depth > dzsum) ii=ih
   dzsum=dzsum+th(ih)
 enddo
 ii=min0(ii,jo)
 cmu=den(ii)*vs(ii)**2
end function cmu
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine gf_bound(st_x_n,st_y_e,ss_min,ss_max)
 use source_params
 use gf_comm
 implicit NONE
 real, intent(IN):: st_x_n,st_y_e
 real, intent(OUT):: ss_min,ss_max
 integer:: iy,ix,js,ks,nfre
 real:: x_ps,y_ps,depth_ps,dis,zz_min,zz_max
!
 ss_min = 100000.0
 ss_max =-100000.0
 do iy=1,nsuby
   zz_min  = 100000.0
   zz_max  =-100000.0
   ngz=3
   do ix=1,nsubx
     js=(ix-1)*nsuby+iy
     do ks=1,nps_sub
       x_ps     = ps_xn(js)+xps_sub(ks)-st_x_n
       y_ps     = ps_ye(js)+yps_sub(ks)-st_y_e
       depth_ps = ps_dz(js)+zps_sub(ks)
       dis=sqrt(x_ps*x_ps+y_ps*y_ps)
       ss_min=amin1(ss_min,dis)
       ss_max=amax1(ss_max,dis)
       zz_min=amin1(zz_min,depth_ps)
       zz_max=amax1(zz_max,depth_ps)
     enddo
   enddo
   depth_min(iy)=zz_min
   depth_max(iy)=zz_max
   ngz=max0(ngz, ifix((zz_max-zz_min)/dep_step+3.05) )
 enddo
 ngr=ifix( (ss_max-ss_min)/d_step+3.05)
 if(allocated(green)) deallocate(csw,trmd, green)
 allocate(csw(2,npt/2), trmd(ngr,ngz), green(npt,8,ngr,ngz))
!
 if((ss_min+0.5*d_step<dist_min).or.(ss_max-0.5*d_step>dist_max)) then
    write(*,*)"ss_min dist_min ss_max dist_max"
    write(*,*) ss_min,dist_min,ss_max,dist_max
    stop
 endif
 zz_min=minval(depth_min)
 zz_max=maxval(depth_max)
 if((zz_min+0.9*dep_step<dep_min).or.(zz_max-0.9*dep_step>dep_max)) then
    write(*,*)"zz_min dep_min zz_max dep_max"
    write(*,*) zz_min,dep_min,zz_max,dep_max
    stop
 endif
end subroutine gf_bound
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine inter_green(rr,zz,rptm,dt,npt,green_ps)
 use gf_comm
 implicit NONE
 integer, intent(IN):: npt
 real, intent(IN):: rr,zz,rptm,dt
 real, dimension(npt,8), intent(OUT):: green_ps
 integer:: ix1,ix2,iz1,iz2,n,n0,n1,n2,k1,k2,ncm,nfre,i1,i2,i
 real:: rx1,rz1,trf,r11,r21,r12,r22,t_ch,wt,cs,sn,gr,gi
!
 ix1=min0(nx_e-nx_b+1, max0(ifix((rr-rg_first)/d_step)+1, 1))
 iz1=min0(nz_e-nz_b+1, max0(ifix((zz-zg_first)/dep_step)+1, 1))
 ix2=min0(ix1+1,nx_e-nx_b+1)
 iz2=min0(iz1+1,nz_e-nz_b+1)
!
 rx1=amin1(1.0, amax1(0.0, ix1-(rr-rg_first)/d_step))
 rz1=amin1(1.0, amax1(0.0, iz1-(zz-zg_first)/dep_step))
 trf=2./float(npt)
 r11=rx1*rz1*trf
 r21=(1.0-rx1)*rz1*trf
 r12=rx1*(1.-rz1)*trf
 r22=(1.0-rx1)*(1.-rz1)*trf
!
 t_ch=rptm+ (r11*trmd(ix1, iz1)+r21*trmd(ix2, iz1) + &
             r12*trmd(ix1, iz2)+r22*trmd(ix2, iz2))/trf
! if(t_ch < 0.0) then
!   n0=-ifix(0.5-t_ch/dt)
!   n1=1-n0
!   n2=npt
! else
!   n0=ifix(0.5+t_ch/dt)
!   n1=1
!   n2=npt-n0
! endif
! k1=n1+n0
! k2=n2+n0
! if(n0>0) green_ps(1:n0, :) =0.0
! green_ps(k1:k2, :)=r11*green(n1:n2, :, ix1, iz1) + &
!                    r21*green(n1:n2, :, ix2, iz1) + &
!                    r12*green(n1:n2, :, ix1, iz2) + &
!                    r22*green(n1:n2, :, ix2, iz2)
! if(n2 < npt) green_ps(n2+1:npt, :)=0.0
!
 nfre=npt/2-1
 do i=1,nfre
   wt=i*dw*t_ch
   csw(1,i)= cos(wt)
   csw(2,i)=-sin(wt)
 enddo
 green_ps=r11*green(:, :, ix1, iz1)+r21*green(:, :, ix2, iz1)+ &
          r12*green(:, :, ix1, iz2)+r22*green(:, :, ix2, iz2)
 do ncm=1,8
   call realft(green_ps(:,ncm),npt,-1)
   i1=1
   do i=1,nfre
     i1=i1+2
     i2=i1+1
     cs=csw(1,i)
     sn=csw(2,i)
     gr=green_ps(i1,ncm)
     gi=green_ps(i2,ncm)
     green_ps(i1,ncm)=gr*cs - gi*sn
     green_ps(i2,ncm)=gr*sn + gi*cs
   enddo
 enddo
end subroutine inter_green
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine input_source_param(filesour,npsx,npsy)
! Read the source parameters
 use source_params
 implicit NONE
 character(len=72):: filesour
 integer, intent(IN):: npsx,npsy
 integer:: i,j,k,ks,ncx,ncy
 real:: dr,dx,dy,rstrike,rdip,xij,yij
!
 open(unit=19, file=filesour,status='old', position='rewind')
 read(19,*) is_moment,npsour,id_sf_type,cft1
!
 nps_sub=npsx*npsy
 if(allocated(ps_xn)) then
    deallocate(ps_xn,ps_ye,ps_dz,ps_mt,ps_rp, &
               ps_rs,ps_p2,ps_sk,ps_dp,ps_rk,drtx,drty, &
               sh_strk,sh_dip,xps_sub,yps_sub,zps_sub )
 endif
 allocate(ps_xn(npsour),ps_ye(npsour),ps_dz(npsour),ps_mt(npsour), &
          ps_rp(npsour),ps_rs(npsour),ps_p2(npsour),ps_sk(npsour), &
          ps_dp(npsour),ps_rk(npsour),drtx(npsour),drty(npsour) )
 allocate(sh_strk(nps_sub),sh_dip(nps_sub),xps_sub(nps_sub), &
          yps_sub(nps_sub),zps_sub(nps_sub))
!
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
 dr=4.0*atan(1.0)/180.0
 do ks=1,npsour
   ps_xn(ks)=ps_xn(ks)/1000.0
   ps_ye(ks)=ps_ye(ks)/1000.0
   ps_dz(ks)=ps_dz(ks)/1000.0
   ps_sk(ks)=ps_sk(ks)*dr
   ps_dp(ks)=ps_dp(ks)*dr
   ps_rk(ks)=ps_rk(ks)*dr
 enddo
!
 drtx=0.0;  drty=0.0
 do i=2,nsubx-1
 do j=1,nsuby
   k=(i-1)*nsuby+j
   drtx(k)=0.5*(ps_rp(k+nsuby)-ps_rp(k-nsuby))/dxsub
 enddo
 enddo
 i=1
 do j=1,nsuby
   k=(i-1)*nsuby+j
   drtx(k)=(ps_rp(k+nsuby)-ps_rp(k))/dxsub
 enddo
 i=nsubx
 do j=1,nsuby
   k=(i-1)*nsuby+j
   drtx(k)=(ps_rp(k)-ps_rp(k-nsuby))/dxsub
 enddo

 do i=1,nsubx
 do j=2,nsuby-1
   k=(i-1)*nsuby+j
   drty(k)=0.5*(ps_rp(k+1)-ps_rp(k-1))/dysub
 enddo
 enddo
 j=1
 do i=1,nsubx
   k=(i-1)*nsuby+j
   drty(k)=(ps_rp(k+1)-ps_rp(k))/dysub
 enddo
 j=nsuby
 do i=1,nsubx
   k=(i-1)*nsuby+j
   drty(k)=(ps_rp(k)-ps_rp(k-1))/dysub
 enddo

 ncx=npsx/2+1
 ncy=npsy/2+1
 dx=dxsub/npsx
 dy=dysub/npsy
 rstrike=ps_sk(1)-angle_north_x
 rdip=ps_dp(1)
 do  i=1,npsx
 do  j=1,npsy
   k=(i-1)*npsy+j
   xij=(i-ncx)*dx
   yij=(j-ncy)*dy
   sh_strk(k)=xij
   sh_dip(k) =yij
   xps_sub(k)=xij*cos(rstrike)-yij*sin(rstrike)*cos(rdip)
   yps_sub(k)=xij*sin(rstrike)+yij*cos(rstrike)*cos(rdip)
   zps_sub(k)=yij*sin(rdip)
 enddo
 enddo
end subroutine input_source_param
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Temporally, case (8) becomes yoffe function for the convenience
! Chen Ji, 2020. It is truncated Kostrov function before
!
subroutine source_fun(id_sf_type,np1,np,dt,moment,rise,sp2,fsour)
 implicit none
 real, parameter:: pi=3.1415926, pi2=2.*pi, pi25=2.5*pi, &
                   c41=0.6*0.5,c42=0.5-c41,c43=2.*c42, &
                   c51=0.5*0.5,c52=0.5-c51
 integer, intent(IN):: id_sf_type,np
 integer, intent(OUT):: np1
 real, intent(IN):: dt,moment,rise,sp2
 real, dimension(np), intent(out):: fsour
 integer:: i,n1,n2
 real:: cnorm,spa,cft4,cft5,tt1,wc,t
 real:: Tp,Te,Tr,rv,psv
 integer:: nn0,nn1,nn2,a_seed
!
  real:: tsin,ty,sn,tsin_min
  integer:: npt_yoffe,nsin,j,nall
  real, dimension(np)::hsin,yoffe
!
 cft4=sp2;  cft5=sp2
 fsour=0.0

 select case (id_sf_type)

!Cases
 case(1)
   wc=pi2/rise
   np1=int(3.0*rise/dt+1.5)
   do i=1,np1
     t=(i-1)*dt
     fsour(i)=wc*wc*t*exp(-wc*t)
   enddo
 case(2)
   wc=4.
   np1=ifix(rise/dt)+1
   do i=1,np1
     t=(i-1)*dt/rise
     fsour(i)=t*(1.-t)**wc
   enddo
!
 case(3)
   np1=ifix(rise/dt)+1
   tt1=0.2*rise
   wc=0.5*pi/tt1
   spa=pi/(rise-tt1)
   do i=1,np1
     t=(i-1)*dt
     if(t < tt1) then
       fsour(i)=sin(wc*t)
     else
       fsour(i)=0.5*(1.+cos(spa*(t-tt1)))
     endif
   enddo
!
 case(4)
   wc=pi/(cft4*rise)
   np1=ifix(rise/dt)+1
   spa=cft4/(1.-cft4)
   do i=1,np1
     t=(i-1)*dt
     tt1=wc*t
     if(tt1 < pi) then
       fsour(i)=c43*sin(0.5*tt1)
     else
       fsour(i)=c42*(1+cos((tt1-pi)*spa))
     endif
     if(tt1 < pi2) then
       fsour(i)=fsour(i)+c41*(1.-cos(tt1))
     endif
   enddo
!
 case(5)
   wc=pi/(cft5*rise)
   np1=ifix(rise/dt)+1
   spa=cft5/(1.-cft5)
   do i=1,np1
     t=(i-1)*dt
     tt1=wc*t
     if(tt1 < pi) then
       fsour(i)=sin(0.5*tt1)
     else
       fsour(i)=c52*(1+cos((tt1-pi)*spa))
!       if(tt1 < pi25) fsour(i)=fsour(i)+c51*(1+cos((tt1-pi)/1.5))
       if(tt1 < pi2) fsour(i)=fsour(i)+c51*(1.0+cos(tt1-pi))
     endif
   enddo
!
 case(6)
   spa=rise+sp2
   np1=ifix(spa/dt)+1
   do i=1,np1
     t=(i-1)*dt
     if(t < rise) then
       fsour(i)=sin(0.5*pi*t/rise)
     else
       fsour(i)=0.5*(1.+cos((t-rise)*pi/sp2))
     endif
   enddo

 case(8)
! modified yoffe function
!   Chen ji, 2020
    fsour=0.0
    if(sp2.ge.1.0.or.sp2.lt.0.0)THEN
      write(*,*)"please check the code, sp2 = ",sp2
      stop
    endif
    tsin_min=2.0*dt
    tsin=rise*sp2
    ty=rise-tsin

    npt_yoffe=int(ty/dt+1.0)+1
    nsin=int(tsin/dt+0.1)+1
    nall=nsin+npt_yoffe
    np1=nall
    do i=1,npt_yoffe
      t=i*dt
      if(t.le.ty)then
          yoffe(i)=sqrt(ty-t)/sqrt(t)
      else
        yoffe(i)=0.0
      endif
    enddo
    if(tsin.ge.tsin_min)then
      do i=1,nsin
        t=(i-1)*dt
        if(t.le.tsin)then
          hsin(i)=sin(t*pi/tsin)
        else
          hsin(i)=0.0
        endif
      enddo
      do i=1,npt_yoffe
        do j=1,nsin
          fsour(i+j-1)=fsour(i+j-1)+yoffe(i)*hsin(j)
        enddo
      enddo
    else
      do i=1,npt_yoffe
        fsour(i)=yoffe(i)
      enddo
    endif

    sn=0.0
    do i=1,nall
      sn=sn+fsour(i)
    enddo
    fsour=fsour/(dt*sn)
!=======================================================
case(7)
  !Tp=sp2*rise
  call random_seed(a_seed)
  call random_number(Tp)
  Tp=0.1*Tp+0.15
  Te=0.8*rise
  Tr=rise
 if(Tr<Tp) Tp=0.25*Tr
  nn0=int(Tp/dt+1.0)
  nn1=int(Te/dt+1.0)
  nn2=int(Tr/dt+1.0)

  psv=sqrt(1.+100./(nn0*dt))
  do i=1,nn0
     t=(i-1)*dt
     fsour(i)=t*psv/Tp*sin(0.5*pi/Tp*t)
  enddo
  do i=nn0+1,nn1
     t=(i-1)*dt
     fsour(i)=sqrt(1.+100./t)
  enddo
  do i=nn1+1,nn2
     t=(i-1)*dt
     fsour(i)=sqrt(1.+100./t)*sin((nn2-i)*dt*pi*0.5/(Tr-Te))
  enddo
  np1=nn2
!
 case default
   np1=0
   fsour=0.0
 end select
!
 if(np1 > 0) then
   cnorm=moment/sum(fsour(1:np1))
   fsour(1:np1)=fsour(1:np1)*cnorm
 endif
end subroutine source_fun
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine conv_sf_gf(n,ndur,syn,fsour,gfun)
 implicit NONE
 integer, intent(IN):: n,ndur
 real, dimension(n):: fsour,gfun,syn
 integer:: ndur1,i,j
 real:: sum1
!
 call realft(gfun,n,1)
 ndur1=ndur-1
 do i=1,ndur1
   sum1=0.0
   do j=1,i
     sum1=sum1+fsour(j)*gfun(i-j+1)
   enddo
   syn(i)=sum1
 enddo
 do i=ndur,n
   sum1=0.0
   do j=1,ndur
     sum1=sum1+fsour(j)*gfun(i-j+1)
   enddo
   syn(i)=sum1
 enddo
end subroutine conv_sf_gf
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine vel_2_acc_or_disp(id_motion,n,dt,gfun)
 implicit NONE
 integer, intent(IN):: id_motion, n
 real, intent(IN):: dt
 real, dimension(n):: gfun
 integer:: i
 real:: dt05,v1,v2
!
 if(id_motion == 1) then
   dt05=0.5*dt
   v1=gfun(1)
   gfun(1)=v1*dt
   do i=2,n
     v2=gfun(i)
     gfun(i) = gfun(i-1)+dt05*(v1+v2)
     v1=v2
   enddo
 endif
 if(id_motion == 3) then
   do i=1,n-1
     gfun(i) = (gfun(i+1)-gfun(i))/dt
   enddo
   gfun(n)=0.0
 endif
end subroutine vel_2_acc_or_disp
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine output3(tmfile,id_format,nump,dt,velc)
 implicit NONE
 character(len=72):: tmfile,file_tmp
 integer, intent(IN):: nump,id_format
 real, intent(IN):: dt
 real, dimension(nump), intent(IN):: velc
 integer:: nerr,myfh,getlen,nlch,j
 integer:: initxdr,ixdrint,ixdrreal,ixdrrmat,ixdrclose
 real:: b0
!
 if(id_format == 1) then
!  SAC
   b0=0.0
   file_tmp=tmfile
   call wsac1(file_tmp,velc,nump,b0,dt,nerr)

 else if(id_format == 2) then
!  TXT
   open(19, file=tmfile, status='unknown', position='rewind')
   write(19,*) nump,dt
   write(19,925) (velc(j), j=1,nump)
   close(19)
   925  format(5e14.5)

 else if(id_format == 3) then
!  FXDR
!   nlch= getlen(tmfile)
!   file_tmp=''
!   file_tmp(1:nlch)=tmfile(1:nlch)
!   myfh = initxdr( file_tmp, 'w', .FALSE. )
!   nerr = ixdrint(  myfh, nump )
!   nerr = ixdrreal( myfh, dt )
!   nerr = ixdrrmat( myfh, nump, velc )
!   nerr = ixdrclose( myfh )

 else if(id_format == 4) then
!  Binary
   open(19, file=tmfile,status='unknown',form='unformatted')
   write(19) nump,dt
   write(19) (velc(j), j=1,nump)
   close(19)
 else
   write(*,*) 'No Format can find for given index'
 endif
end  subroutine output3
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine read_green(npt,dt,d_min,d_max,z_min,z_max)
!  Read Greens function from the Green Function bank
! use source_params
 use gf_comm
 implicit NONE
 integer, intent(IN):: npt
 real, intent(IN):: dt,d_min,d_max,z_min,z_max
 integer:: k,nta,ntb,nnormal,kall_x,kall_z
 integer:: ll,iz,k0,iz0,izz,ixx,npt_c,ntc,nc
 real:: dis0,t0,dep,dt_c,vmove,t_ch
!
! if(z_min>=(zg_first-1.e-7) .and. z_max<=(zg_last+1.e-7)) return
!
 nx_b=max0(ifix((d_min-dist_min)/d_step)+1, 1)
 nx_e=min0(ifix((d_max-dist_min)/d_step)+2, nx_in)
 nz_b=max0(ifix((z_min-dep_min)/dep_step)+1, 1)
 nz_e=min0(ifix((z_max-dep_min)/dep_step)+2, nz_in)
 kall_x=nx_e-nx_b+1
 kall_z=nz_e-nz_b+1
 if(kall_x > ngr .or. kall_z > ngz)then
    ngr=max0(kall_x, ngr)
    ngz=max0(kall_z, ngz)
    deallocate(trmd, green)
    allocate(trmd(ngr,ngz), green(npt,8,ngr,ngz))
    green=0.0
 endif
!
 nta=ifix(0.1*t_cor/dt+0.5)+1
 ntb=max0(ifix(0.8*t_cor/dt+0.5), nta)
 nnormal=ntb-nta+1
!
 open(50,file=file_bank,status='old',access='direct',recl=block_gg)
 do iz=nz_b,nz_e
    izz=iz-nz_b+1
    do k=nx_b,nx_e
      ixx=k-nx_b+1
      ll=(iz-1)*nx_in+k
      read(50,rec=ll)iz0,k0,dis0,t0,dep,dt_c,npt_c, &
         ((green(ntc,nc,ixx,izz),ntc=1,npt_c),nc=1,8)
      if(izz == 1) zg_first=dep
      if(ixx == 1) rg_first=dis0
      trmd(ixx,izz)=t0
      do nc=1,8
        vmove=sum(green(nta:ntb,nc,ixx,izz))/float(nnormal)
        green(1:npt_c,nc,ixx,izz)=green(1:npt_c,nc,ixx,izz)-vmove
      enddo
    enddo
 enddo
 zg_last=dep
 close(50)
end subroutine read_green
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rad_pattern(iseed,az0,paz,rake0,dip0,theta0,ang_h1,ang_h2)
 use source_params
 implicit NONE
 real, intent(IN):: az0,rake0,dip0,theta0,paz,ang_h1,ang_h2
 integer, intent(INOUT):: iseed
 integer, save:: irum=-7
 integer:: nfre,i,j,i1,i2
 real:: az,rake,strike,theta,dip,ran_az,ran_dp,ran_rk,ran1,wl,aj
 real:: cosr,sinr,cosd,sind,cos2d,sin2d,cosa,sina,cos2a,sin2a, &
        sinh1,sinh2,cosh1,cosh2,au(5),a(5)
!
 nfre=npt/2
! ran_dp=ran1(iseed)-0.5
! ran_rk=ran1(iseed)-0.5
! ran_az=ran1(irum)-0.5
 ran_dp=0.0
 ran_rk=0.0
 ran_az=0.0
 do i=1,nfre
   i2=2*i
   i1=i2-1
   if(i==1 .or. (i>kfb .and. i<kfe)) then
     rake=rake0+   ran_rk*rake_bud(i)
     dip=dip0+     ran_dp*dip_bud(i)
     strike=theta0
     cosr=cos(rake)
     sinr=sin(rake)
     cosd=cos(dip)
     sind=sin(dip)
     cos2d=cos(2.*dip)
     sin2d=sin(2.*dip)
!
!     az=az0+ran_az*strk_bud(i)
     az=az0+paz*strk_bud(i)
     theta=az-strike
     cosa=cos(theta)
     sina=sin(theta)
     cos2a=cos(2.*theta)
     sin2a=sin(2.*theta)
!
     sinh1=-sin((az-ang_h1))
     cosh1= cos((az-ang_h1))
     sinh2=-sin((az-ang_h2))
     cosh2= cos((az-ang_h2))
!
     au(3)= 0.5*sinr*sin2d
     au(2)= cosr*cosd*cosa-sinr*cos2d*sina
     au(1)= 0.5*sinr*sin2d*cos2a+cosr*sind*sin2a
     au(5)=-sinr*cos2d*cosa-cosr*cosd*sina
     au(4)= cosr*sind*cos2a-0.5*sinr*sin2d*sin2a
!
     wl=rd_wl(i)
     do j=1,5
       aj=amax1(abs(au(j)), wl)
       a(j)= sign(aj, au(j))
     enddo
!
!   coef_rd(i1,1)=a(4)
!   coef_rd(i2,1)=a(4)
!   coef_rd(i1,2)=a(5)
!   coef_rd(i2,2)=a(5)
!   coef_rd(i1,3)=a(1)
!   coef_rd(i2,3)=a(1)
!   coef_rd(i1,4)=a(2)
!   coef_rd(i2,4)=a(2)
!   coef_rd(i1,5)=a(3)
!   coef_rd(i2,5)=a(3)
!
   endif
! H1
   rad_h1(i1,1)=sinh1*a(4)
   rad_h1(i2,1)=rad_h1(i1,1)
   rad_h1(i1,2)=sinh1*a(5)
   rad_h1(i2,2)=rad_h1(i1,2)
   rad_h1(i1,3)=cosh1*a(1)
   rad_h1(i2,3)=rad_h1(i1,3)
   rad_h1(i1,4)=cosh1*a(2)
   rad_h1(i2,4)=rad_h1(i1,4)
   rad_h1(i1,5)=cosh1*a(3)
   rad_h1(i2,5)=rad_h1(i1,5)
! H2
   rad_h2(i1,1)=sinh2*a(4)
   rad_h2(i2,1)=rad_h2(i1,1)
   rad_h2(i1,2)=sinh2*a(5)
   rad_h2(i2,2)=rad_h2(i1,2)
   rad_h2(i1,3)=cosh2*a(1)
   rad_h2(i2,3)=rad_h2(i1,3)
   rad_h2(i1,4)=cosh2*a(2)
   rad_h2(i2,4)=rad_h2(i1,4)
   rad_h2(i1,5)=cosh2*a(3)
   rad_h2(i2,5)=rad_h2(i1,5)
! vertical
   rad_vt(i1,1)=a(1)
   rad_vt(i2,1)=a(1)
   rad_vt(i1,2)=a(2)
   rad_vt(i2,2)=a(2)
   rad_vt(i1,3)=a(3)
   rad_vt(i2,3)=a(3)
 enddo
end subroutine rad_pattern
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function getlen(string)
 implicit NONE
 character(len=72):: chatmp, string
 integer:: getlen,nl
 chatmp=''
 chatmp=adjustl(string)
 nl=len_trim(chatmp)
 string=''
 string(1:nl)=chatmp(1:nl)
 getlen=nl
end function getlen
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine num2ch(r,ch)
 implicit NONE
 real, intent(IN):: r
 character(len=72),intent(out):: ch
 integer:: i,n,n1,n2,nc,nl,nzero
!
 ch=''
 n=ifix(r+0.001)
 nc=3
 nl=10**(nc-1)
 nzero=ichar('0')
 n2=n
 do i=1,nc
   n1=n2/nl
   ch(i:i)=char(nzero+n1)
   n2=n2-n1*nl
   nl=nl/10
 enddo
end subroutine num2ch
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine distaz(stat_x_n,stat_y_e,x_ps,y_ps,afa,dis,azz,baz)
 implicit NONE
 real, intent(IN):: stat_x_n,stat_y_e,x_ps,y_ps,afa
 real, intent(OUT):: dis,azz,baz
 real:: rdp,dx,dy
!
 rdp=180./3.1415926
 dx=stat_x_n-x_ps
 dy=stat_y_e-y_ps
 dis=sqrt(dx*dx+dy*dy)
 azz=acos(abs(dx)/dis)*rdp
 if(dx < 0.0) azz=180.0-azz
 if(dy < 0.0) azz=360.-azz
 azz=azz+afa*rdp
 if(azz > 360.0) azz=azz-360.0
 baz=azz+180.0
 if(baz > 360.0) baz=baz-360.0
 azz=azz/rdp
 baz=baz/rdp
end subroutine distaz
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE realft(data,n,isign)
!  USES four1
!
!  Calculates the Fourier transform of a set of real data points.
!  n must be a power of 2.
!
!  isign=-1, forward tranform of a real-valued data point.
!     data(2*i+1) and data(2*i+2) return the real and image part of
!     frequency f=i*df. The real-valued first (zero frequency) and
!     last (fmax) components of the complex transform are returned
!     as element data(1) and data(2), respectively.
!  isign=+1, inverse transform of a complex data array. The result
!     in this case must be multipled by 2/n (or 2*df).
!
 implicit NONE
 INTEGER, INTENT(IN) :: isign,n
 REAL, INTENT(INOUT) :: data(n)
 INTEGER :: i,i1,i2,i3,i4,n2p3
 REAL :: c1,c2,h1i,h1r,h2i,h2r,wis,wrs
 REAL (KIND=8) :: theta,wi,wpi,wpr,wr,wtemp
 theta=3.141592653589793d0/dble(n/2)
 c1=0.5
 if (isign==-1) then
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
 do i=2,n/4
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
 enddo
 data(n/2+2)=-data(n/2+2)
 if (isign==-1) then
   h1r=data(1)
   data(1)=h1r+data(2)
   data(2)=h1r-data(2)
 else
   h1r=data(1)
   data(1)=c1*(h1r+data(2))
   data(2)=c1*(h1r-data(2))
   call four1(data,n/2,isign)
 endif
END SUBROUTINE realft
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE four1(data,nn,isign)
 implicit NONE
 INTEGER, INTENT(IN) :: isign,nn
 REAL, DIMENSION(2*NN), INTENT(INOUT):: data
 INTEGER :: i,istep,j,m,mmax,n
 REAL :: tempi,tempr
 REAL (KIND=8) :: theta,wi,wpi,wpr,wr,wtemp
 n=2*nn
 j=1
 do i=1,n,2
   if(j > i)then
     tempr=data(j)
     tempi=data(j+1)
     data(j)=data(i)
     data(j+1)=data(i+1)
     data(i)=tempr
     data(i+1)=tempi
   endif
   m=n/2
   do while ((m>=2).and.(j>m))
     j=j-m
     m=m/2
   enddo
   j=j+m
 enddo
 mmax=2
 do while (n>mmax)
   istep=2*mmax
   theta=6.28318530717959d0/(isign*mmax)
   wpr=-2.d0*sin(0.5d0*theta)**2
   wpi=sin(theta)
   wr=1.d0
   wi=0.d0
   do m=1,mmax,2
     do i=m,n,istep
       j=i+mmax
       tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
       tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
       data(j)=data(i)-tempr
       data(j+1)=data(i+1)-tempi
       data(i)=data(i)+tempr
       data(i+1)=data(i+1)+tempi
     enddo
     wtemp=wr
     wr=wr*wpr-wi*wpi+wr
     wi=wi*wpr+wtemp*wpi+wi
   enddo
   mmax=istep
 enddo
END SUBROUTINE four1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function ran1(idum)
!  Generate a uniformly distributed random number on the interval [0,1].
 implicit NONE
 integer, intent(INOUT):: idum
 integer, parameter:: im=139968, ia=3877, ic=29573
 real:: ran1
!
 idum=mod(idum*ia+ic,im)
 ran1=float(idum)/float(im)
end function ran1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION gasdev(idum,iset)
 implicit NONE
 INTEGER:: idum,iset
 REAL::  gasdev
 REAL:: fac,rsq,v1,v2,ran3
 REAL, SAVE:: gset
 if (iset == 0) then
   do while(1==1)
     v1=2.*ran3(idum)-1.
     v2=2.*ran3(idum)-1.
     rsq=v1**2+v2**2
     if(rsq < 1.0 .and. rsq /= 0.0) exit
   enddo
   fac=sqrt(-2.*log(rsq)/rsq)
   gset=v1*fac
   gasdev=v2*fac
   iset=1
 else
   gasdev=gset
   iset=0
 endif

END FUNCTION gasdev
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION ran3(idum)
! ROUTINE TO GENERATE A UNIFORMLY DISTRIBUTED RANDOM
! NUMBER ON THE INTERVAL [0,1].
!
 implicit NONE
 INTEGER, PARAMETER:: IA=16807,IM=2147483647,IQ=127773,IR=2836, &
                      NTAB=32,NDIV=1+(IM-1)/NTAB
 REAL, PARAMETER:: AM=1./IM, EPS=1.2e-7,RNMX=1.-EPS
 INTEGER:: idum,j,k
 INTEGER, SAVE:: iy=0
 INTEGER, DIMENSION(NTAB), SAVE:: iv=0
 REAL:: ran3

 if (idum <= 0 .or. iy == 0) then
   idum=max(-idum,1)
   do j=NTAB+8,1,-1
     k=idum/IQ
     idum=IA*(idum-k*IQ)-IR*k
     if (idum < 0) idum=idum+IM
     if (j <= NTAB) iv(j)=idum
   enddo
   iy=iv(1)
 endif
 k=idum/IQ
 idum=IA*(idum-k*IQ)-IR*k
 if (idum < 0) idum=idum+IM
 j=1+iy/NDIV
 iy=iv(j)
 iv(j)=idum
 ran3=min(AM*iy,RNMX)

END FUNCTION ran3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine kappa(n,dt,y,kap)
 implicit none
 integer :: i,n
 real :: df, dt, y(n), pi=4.0*ATAN(1.0), kap
 complex :: c(n)
 do i=1,n
    c(i)=CMPLX(y(i),0.0)
 enddo
 call fork(n,c,-1.)
 df=1/(n*dt)
 do i=1,n/2
    c(i)=c(i)*exp(-pi*(i-1)*df*kap)
 enddo
 do i=n/2+2,n
    c(i)=CONJG(c(n+2-i))
 enddo
 call fork(n,c,1.)
 do i=1,n
    y(i)=real(c(i))
 enddo
end
!JORGE calculate Fourier Spectra (from Boore)
! ----------------------------- BEGIN FORK --------------------------
SUBROUTINE FORK(LX,CX,SIGNI)
      complex cx(*),carg,cexp,cw,ctemp

      pi = 4.0*atan(1.0)

      j=1
      sc=sqrt(1./float(lx))

      do i=1,lx
        if(i.gt.j) go to 2
        ctemp=cx(j)*sc
        cx(j)=cx(i)*sc
        cx(i)=ctemp
2       m=lx/2
3       if(j.le.m) go to 5
        j=j-m
        m=m/2
        if(m.ge.1) go to 3
5       j=j+m
      end do

      l=1
6     istep=2*l
      temp= pi * signi/float(l)

      do m=1,l
        carg=(0.,1.)*temp*(m-1)
        cw=cexp(carg)
        do i=m,lx,istep
          ctemp=cw*cx(i+l)
          cx(i+l)=cx(i)-ctemp
          cx(i)=cx(i)+ctemp
        end do
      end do

      l=istep
      if(l.lt.lx) go to 6

      return
      end
! ----------------------------- END FORK --------------------------
