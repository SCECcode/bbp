!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE source_params
 implicit NONE
 save
 integer:: nsubx,nsuby,npsour,nps_sub,is_moment,id_sf_type,npt,kfb,kfe
 integer, allocatable, dimension(:,:):: ind4_sub
 real:: dxsub,dysub,cxp,cyp,xref_hypc,yref_hypc,angle_north_x,dt,cft1
 real, allocatable, dimension(:):: ps_xn,ps_ye,ps_dz,ps_mt,ps_rp, &
                         ps_rs,ps_p2,ps_sk,ps_dp,ps_rk,drtx,drty, &
                         sh_strk,sh_dip,xps_sub,yps_sub,zps_sub
 real, allocatable, dimension(:):: strk_bud,dip_bud,rake_bud,rd_wl, &
                                   fsour,syn_up,syn_h1,syn_h2,syn_tmp
 real, allocatable, dimension(:,:):: gfun,green_ps,rad_h1,rad_h2,rad_vt
 real, allocatable, dimension(:,:):: coef4_ps 
END MODULE source_params
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE station_params
 implicit NONE
 save  
 integer:: nst
 real:: xs0,ys0,afa,st_x_n,st_y_e,st_z_up,ang_h1,ang_h2,ss_min,ss_max
 character(len=72), allocatable, dimension(:):: stan_name
 real, allocatable, dimension(:):: stan_x,stan_y,stan_z,stan_h1,stan_h2
END MODULE station_params
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE gf_comm
 implicit NONE
 save  
 character(len=72):: file_bank
 integer:: lnpt,block_gg,jo,nb,nx_in,nz_in,nlay
 integer:: nz_b,nz_e,nx_b,nx_e,ngz,ngr
 real:: dep_max,dep_min,dep_step,dist_max,dist_min,d_step,t_cor
 real:: rg_first,zg_first,zg_last,dw
 real, allocatable, dimension(:):: vp,vs,den,th,qa,qb,depth_min,depth_max
 real, allocatable, dimension(:,:):: trmd,csw 
 real, allocatable, dimension(:,:,:,:):: green 
END MODULE gf_comm
