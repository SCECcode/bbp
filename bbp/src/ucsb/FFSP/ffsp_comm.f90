!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE sp_sub_f
 implicit NONE
 save
 integer:: nsubx,nsuby,nsum,ncxp,ncyp,nseg,id_sf_type,layer
 integer:: ir_Mw,ir_Dip,ir_Fl,ir_Fw
!
 integer:: id_cl=1
 integer:: id_lrtp=1
!
 real:: rstm_mean,pktm_mean,td_target,t_precent
 real:: Mw,Moment_o,fc_main_1,fc_main_2,flx,fwy,dx,dy,x_hypc,y_hypc,depth_hypc, &
        lat_hypc,lon_hypc

 real:: smoment,hfrd,sp1,sq1,slip_clx,slip_cly,slip_pw, &
        vmin_ratio,vmax_ratio,vp1,vq1,rv_clx,rv_cly,rv_pw,&
        rs_max,pk_max,rp1,rq1,rs_clx,rs_cly,rs_pw,cr_sl_rs,cr_sl_vc,&
        cr_pk_vc,pk_pw
 real, allocatable, dimension(:):: slip,rstm,rptm,rpvel,pktm,lrtp,&
                                   rake_prt,dip_prt,amz_prt
 real, allocatable, dimension(:):: beta,amu,taper
 real, allocatable, dimension(:):: vvp,vvs,roh,thk,qp,qs
 real, allocatable, dimension(:):: lx_seg,str_seg,dip_seg,rak_seg, &
                                   lon_seg,lat_seg,xfps,yfps,zfps
 integer, allocatable, dimension(:):: nxx_seg
END MODULE sp_sub_f
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE fdtim_2d
 implicit NONE
 save
 integer:: nx_2d,nz_2d
 real cxp,cyp,grid_2d
 real, allocatable, dimension(:):: rtx,rtz
END MODULE fdtim_2d
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE time_freq
 implicit NONE
 save
 integer:: ntime,nphf,nfe1,nfe2
! Add four new parameters for DCF, Chen Ji, 2020
 integer:: lnpt,nbb,nee,io_dva
!
 real:: dt,df,cft1
 real, allocatable, dimension(:):: freq
END MODULE time_freq
