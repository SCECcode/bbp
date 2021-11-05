!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE source_params_all
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
!
real::flx,fwy,strike,dip,rake
real::freq_min,freq_max
real::x_hypc, y_hypc, depth_hypc
real::angle_north_to_x,rv_avg
real::dx,dy,fc_main_1,fc_main_2
!
 integer::n_fault_segment,id_segment
 integer, allocatable, dimension(:)::io_reverse,nx_stk_max,nx_stk_min,ny_dip_max,ny_dip_min
 real, allocatable, dimension(:):: x_stk_min,x_stk_max,y_dip_min,y_dip_max
 real, allocatable, dimension(:):: seg_stk,seg_dip,seg_rake
 real, allocatable, dimension(:):: seg_upleft_n,seg_upleft_e,seg_upleft_z

END MODULE source_params_all
