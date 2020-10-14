Program ffsp_dcf_v3
!
! This code output the spatial distribution of source parameters
!
! The relationship between Moment magnitude and moment
!           log10(Mo) = 1.5 Mw+16.1
! The unit of Mo is dyn-cm (Kanamori 1977)
! The unit of input parameters are km (length), km/s (Velocity),
!       g/cm^3 (density), degree (Angle)
!
! This code call the subroutines or functions in source_param2.f
!
! Written by Pengcheng Liu
! Copyright (c) 2005 by Pengcheng Liu
!

! adding the relationships found by Schmedes et al. (2010, 2013)
! Chen Ji, 2020
! A new function is defined to calculate the quality of source model.
! The best model is selected as spfile_best
!======================================================================
 use sp_sub_f
 implicit NONE
!
 character(len=72):: vsfile,spfile,spfile_best,spfile_plot
 character(len=3):: charnum
 integer, dimension(4):: nb_taper_TRBL
 integer:: getlen,idum1,idum2,idum3,id_ran1,id_ran2,&
           is_moment,nlch,nsource,i,j,k
 real:: xref_hypc,yref_hypc,strike,dip,rake,angle_north_to_x, &
        freq_min,freq_max,rv_avg,cft1
 real:: sp2,drp,rdip,rstrike,cosx,sinx,str_ref,dip_ref,area_sub, &
        summ,rmt,rvc,xij,yij,xs,ys,xps,yps,zps,factor,rakei,dipi,&
        pamzi,pdip_max,prake_max
 real::rise_time
 real::ave_tr,ave_tp,ave_vr,err_spectra,pdf,pdf_min
 integer::io_model_out,id_min

 open(8,file='ffsp.inp',status='old')
 open(7,file='ffsp.inf',status='unknown')
!
! Inputing the fault parameters
!
 write(7,*) 'Enter the index of slip rate function'
 read(8,*)  id_sf_type, freq_min, freq_max
 write(7,*) id_sf_type, freq_min, freq_max

 write(7,*) ' Enter the length and width of main fault'
 read(8,*)  flx,fwy
 write(7,*) flx,fwy

 write(7,*) ' Enter the distance in the strike and down dip'
 write(7,*) ' direction, and the depth of the hypocenter'
 read(8,*)  x_hypc, y_hypc, depth_hypc
 write(7,*) x_hypc, y_hypc, depth_hypc

 write(7,*) 'Enter the coordinates  of hypocenter'
 read(8,*)  xref_hypc,yref_hypc
 write(7,*) xref_hypc,yref_hypc

 write(7,*) 'Enter the moment(N-m), corner frequency(Hz) and ', &
            'Average rupture speed'

 read(8,*)  Moment_o,fc_main_1,fc_main_2,rv_avg
 write(7,*) Moment_o,fc_main_1,fc_main_2,rv_avg

 write(7,*) 'Enter the ratio of time between slip-rate increase', &
            ' and rise time'

 read(8,*)  cft1
 write(7,*) cft1

 write(7,*) ' Enter the strike, dip and rake of main fault'
 read(8,*)  strike,dip,rake
 write(7,*) strike,dip,rake

 write(7,*) 'Enter Max. perturbation on the dip and rake'
 read(8,*)  pdip_max,prake_max
 write(7,*) pdip_max,prake_max

 write(7,*) ' Enter the subfaults number along strike adn dip'
 read(8,*)  nsubx,nsuby
 write(7,*) nsubx,nsuby

 write(7,*) ' Enter the numbers of subfault to be tapered at each side'
 !           Top-Right-Bottom-Left
 read(8,*)  (nb_taper_TRBL(i), i=1,4)
 write(7,*) (nb_taper_TRBL(i), i=1,4)

 write(7,*) ' Enter the seed for generating random numbers'
 read(8,*)  idum1,idum2,idum3
 write(7,*) idum1,idum2,idum3

 write(7,*)'Enter first and last index of random generations'
 read(8,*)  id_ran1,id_ran2
 write(7,*) id_ran1,id_ran2
!
!  Enter name of files with 1D velocity structure
!
 write(7,*) ' Enter the name of file with velocity', &
            ' structure of source region'
 read(8,'(1a)')  vsfile
 write(7,'(1a)') vsfile
!
!  The Parameters for outputing results
!
 write(7,*) 'Enter the angle (degree) from North to X-Axis ', &
            '(clockwise) of SYN. code'
 read(8,*)  angle_north_to_x
 write(7,*) angle_north_to_x

 write(7,*) 'Do you output moment(1), slip-area (2), or slip(3) '
 read(8,*)  is_moment
 write(7,*) is_moment

 write(7,*) ' Enter the file name of output source parameters'
 read(8,'(1a)')  spfile
 write(7,'(1a)') spfile
 close(8)

!
 if(Moment_o < 12.0) Moment_o=10**(1.5*Moment_o+9.05)
 Mw=(alog10(Moment_o)-9.05)/1.5

 sp2=1.0
!------------------------------------------------------------------
! Inputting 1D velocity structure
!  assigning values for vvp, vvs, roh, thk, qp qs
 call xmodel(vsfile)
!
!  in spfield: set the initial parameters for generating the source parameters
!  rv_avg: average rupture velocity
 !  freq_max: using to assign nfe1, nfe2 (0.5 - 0.9 freq_max)
 !  cft1: tp/tr
 !  nb_taper_TRBL

 call mainfault(dip,freq_min,freq_max,rv_avg,cft1,nb_taper_TRBL)
!
 nlch=getlen(spfile)+1
 spfile(nlch:nlch)='.'
 drp=4.*atan(1.0)/180.0
 rstrike=strike*drp
 rdip=dip*drp
 cosx=cos(angle_north_to_x*drp)
 sinx=sin(angle_north_to_x*drp)
 str_ref=(ncxp-0.5)*dx
 dip_ref=(ncyp-0.5)*dy
!
 open(31,file='source_model.score')
 write(31,*)id_ran2-id_ran1+1
 open(11,file='source_model.list')
 write(11,118) 1, nsubx, nsuby, &
               dx*1000., dy*1000., x_hypc*1000.0, y_hypc*1000.0
 write(11,119) xref_hypc*1000., yref_hypc*1000., angle_north_to_x
 area_sub=1.e-15
 if(is_moment == 3) area_sub=area_sub/(dx*dy)
!
 write(31,129)rstm_mean+pktm_mean,pktm_mean
 129     format("Target: average Risetime= ",f8.3," average peaktime= ",f8.4)
 130    format(5e15.5)
 pdf_min=1.0e20
 id_min=id_ran1
 io_model_out=0
 spfile_best=spfile
 spfile_best(nlch:nlch+3)='.bst'
 write(11,'(1a)') spfile_best(1:nlch+3)
 close(11)
 do nsource=id_ran1,id_ran2
   call dig2ch(nsource,charnum)
   spfile(nlch+1:nlch+3)=charnum
   write(31,'(1a)') spfile(1:nlch+3)

   call random_field(idum1,idum2,idum3,ave_tr,ave_tp,ave_vr,err_spectra)
   pdf=((log(ave_tr)-log(rstm_mean+pktm_mean))/0.1)**2.0
   pdf=pdf+((log(ave_tp)-log(pktm_mean))/0.1)**2.0
   pdf=pdf+err_spectra

   write(31,130)ave_tr,ave_tp,ave_vr,err_spectra,pdf
   if(pdf.lt.pdf_min)then
     pdf_min=pdf
     id_min=nsource
     io_model_out=1
   endif
   open(12,file=spfile)
   write(12,*) is_moment,nsum,id_sf_type,cft1
   if(io_model_out.eq.1)then
     open(22,file=spfile_best)
     write(22,*) is_moment,nsum,id_sf_type,cft1
   endif
!
   summ=0.0
   rmt=0.0
   rvc=0.0
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
     factor=1.0
     if(is_moment.gt.1) factor=area_sub/amu(k)
     slip(k)=slip(k)*factor
     rakei=rake+rake_prt(k)*35.0

     dipi=dip+dip_prt(k)*pdip_max
     rakei=rake+rake_prt(k)*prake_max
     pamzi=amz_prt(k)
!     rise_time=rstm(k)+pktm(k)
     write(12,109) xps*1000.,yps*1000.,zps*1000., &
          slip(k),rptm(k),rstm(k),pktm(k),strike,dipi,rakei
     if(io_model_out.eq.1)then
        write(22,109) xps*1000.,yps*1000.,zps*1000., &
          slip(k),rptm(k),rstm(k),pktm(k),strike,dipi,rakei
     endif
     summ=summ+rstm(k)
     rmt=rmt+slip(k)
     rvc=rvc+rpvel(k)
   enddo
   enddo
   close(12)
   close(14)
   write(*,*) 'Moment =', rmt
   write(*,*) 'Avg-rst=', summ/nsubx/nsuby
   write(*,*) 'Avg-rpv=', rvc/nsubx/nsuby
   if(io_model_out.eq.1)then
     close(22)
     io_model_out=0
   endif
 enddo
!
 109     format(4e15.7,6e12.4)
 118     format(I3, 2I5, 4e15.7)
 119     format(3e15.7)
 close(7)
 close(31)

end program ffsp_dcf_v3
!
!======================================================================
!
subroutine xmodel(vsfile)
! Inputting 1D velocity structure from file vsfile
!
 use sp_sub_f
 implicit NONE
 character(len=72):: vsfile
 integer:: j
 real:: freq_ref

 open(17,file=vsfile)
! The number of layer
 read(17,*) layer, freq_ref
!
 allocate(vvp(layer+2),vvs(layer+2),roh(layer+2), &
          thk(layer+2),qp(layer+2),qs(layer+2))
!
! P-wave, S-wave velocity, density, thickness, Qp, and Qs of each layer
 do j=1,layer
   read(17,*) vvp(j),vvs(j),roh(j),thk(j),qp(j),qs(j)
 enddo
 close(17)

end subroutine xmodel
!
!======================================================================
!
subroutine dig2ch(n,ch)
 implicit none
 integer:: n,n1,n2,n3,ntmp,nzero
 character(len=3):: ch

 n1=n/100
 ntmp=n-n1*100
 n2=ntmp/10
 n3=ntmp-n2*10

 nzero=ichar('0')
 ch(1:1)=char(nzero+n1)
 ch(2:2)=char(nzero+n2)
 ch(3:3)=char(nzero+n3)

end subroutine dig2ch
!
!=====================================================================
!
function getlen(string)
 implicit NONE
 character(len=*):: string
 integer:: getlen
 string=adjustl(string)
 getlen=len_trim(string)
end function getlen
