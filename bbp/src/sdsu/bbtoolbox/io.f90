SUBROUTINE cleaner
!--------------------------------------------------------------------
!
! Description:
!
!   Erases some files used by the raytracing subroutine
!
! Dependencies:
!
!   None
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

! flag
integer(kind=i_single):: ierr

!--------------------------------------------------------------------

! open and erase files
open(1,file='vel3d_P.bin',status='unknown',iostat=ierr)
if (ierr /= 0) call error_handling(1,'vel3d_P.bin','cleaner (io.f90)')
close(1,status='delete',iostat=ierr)
if (ierr /= 0) call error_handling(3,'vel3d_P.bin','cleaner (io.f90)')

open(1,file='vel3d_S.bin',status='unknown',iostat=ierr)
if (ierr /= 0) call error_handling(1,'vel3d_S.bin','cleaner (io.f90)')
close(1,status='delete',iostat=ierr)
if (ierr /= 0) call error_handling(3,'vel3d_S.bin','cleaner (io.f90)')

!open(1,file='vel3d_P.in',status='unknown',iostat=ierr)
!if (ierr /= 0) call error_handling(1,'vel3d_P.in','cleaner (io.f90)')
!close(1,status='delete',iostat=ierr)
!if (ierr /= 0) call error_handling(3,'vel3d_P.in','cleaner (io.f90)')

!open(1,file='vel3d_S.in',status='unknown',iostat=ierr)
!if (ierr /= 0) call error_handling(1,'vel3d_S.in','cleaner (io.f90)')
!close(1,status='delete',iostat=ierr)
!if (ierr /= 0) call error_handling(3,'vel3d_S.in','cleaner (io.f90)')

open(1,file='time3d_P.out',status='unknown',iostat=ierr)
if (ierr /= 0) call error_handling(1,'time3d_P.out','cleaner (io.f90)')
close(1,status='delete',iostat=ierr)
if (ierr /= 0) call error_handling(3,'time3d_P.out','cleaner (io.f90)')

open(1,file='time3d_S.out',status='unknown',iostat=ierr)
if (ierr /= 0) call error_handling(1,'time3d_S.out','cleaner (io.f90)')
close(1,status='delete',iostat=ierr)
if (ierr /= 0) call error_handling(3,'time3d_S.out','cleaner (io.f90)')


END SUBROUTINE cleaner

!===================================================================================================

SUBROUTINE log_init
!--------------------------------------------------------------------
!
! Description:
!
!   Initializes the run.log file with values shared by all stations
!
! Dependencies:
!
!   None
!
! Notes: 
!
!   Default log-file unit number is 5
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: October  2013 (v1.5.2)
!   Add output decimation factor, time_step.
!   
! Updated: March 2014 (v1.5.4)
!   time_step is set in scattering.dat.
!
! Updated: July 2015 (v1.6.1)
!   Change tmp_dt computation to use fixed lf_len.
!
! Updated: February 2019 (v2.0)
!   Back to use the original tmp_dt computation, no use tmp_lf_len nor tmp_npts.

use def_kind; use earthquake; use flags; use geometry; use matching
use scattering; use source_receiver; use stf_data; use waveform   
use tmp_para
   
implicit none

! stores date 
character(len=8)      :: d
! stores time
character(len=10)     :: t
! flag
integer(kind=i_single):: ierr

! add for output decimation factor
!integer(kind=i_single),parameter              :: time_step = 1
! decimated npts
integer(kind=i_single)                        :: d_npts
! adjusted dt
real(kind=r_single)                           :: tmp_dt
!--------------------------------------------------------------------

d_npts=npts/time_step
tmp_dt = lf_len/(v_npts-1)*time_step

! Formats
80 format(a30,f0.3,T42,a30,f0.3)
81 format(a30,i0,T42,a30,i0)

! open log-file (may overwrite if file already exists) and check it
open(5,file='run.log',status='unknown',iostat=ierr)
if (ierr /= 0) call error_handling(1,'run.log','log_init (io.f90)')

write(5,*) '*** LOG-FILE FOR BROADBAND HYBRID COMPUTATION CODE (P.M.MAI & K.B.OLSEN) ***' 
write(5,*)

! write date&time
call date_and_time(date=d,time=t)
write(5,*) "Current date and time: ",d(7:8),"/",d(5:6),"/",d(1:4),'  ',t(1:2), &
                                  ":",t(3:4),":",t(5:10)
write(5,*)

write(5,*) "Code running in modality: ",modality_flag

write(5,*)

! writes raytracing properties
write(5,*) '------------------------ RAYTRACING PROPERTIES ---------------------------'
write(5,'(a30,a,2(i0,x),i0,a)') 'Domain dimensions [x y z]   : ','[',nint(h_step*(nx-1)), &
                                 nint(h_step*(ny-1)),nint(h_step*(nz-1)),']'
write(5,'(a30,a,2(i0,x),i0,a)') 'Domain discretization (pts) : ','[',nx,ny,nz,']'                                   

write(5,*)

! writes medium properties
write(5,*) '-------------------------- MEDIUM PROPERTIES -----------------------------'

if (vel_flag == 1) then
   write(5,'(a30,i0)') 'Num. of layers for 1D model : ',n_lay
else
   write(5,*) '3D model has been used'
endif

if (trim(verbose_flag) == 'on' .and. vel_flag == 1) then
   write(5,*) 'Resampled 1D model has been written on a separated log-file'
endif

write(5,*)

! writes source properties
write(5,*) '-------------------------- SOURCE PROPERTIES -----------------------------'

if (modality_flag /= 0) then
   if (ext_flag == 1) then
      write(5,'(a28,i0,a24)') 'An extended-source model of ',n_cell,' subfaults has been used'
   else
      write(5,*) 'Point-source model has been used'
   endif
endif

write(5,'(a30,a,2(f0.3,x),f0.3,a)') 'Hypocenter position at      : ', '[',hyp_x,hyp_y,hyp_z,']'

write(5,*)

if (modality_flag /= 0) then
   ! writes source-time function properties
   write(5,*) '---------------------------- STF PROPERTIES ------------------------------'
   write(5,81) 'Npts in STF                 : ',  npts_stf              
   write(5,80) 'Time length                 : ',total, 'Time-step                   : ',          &
                total/(npts_stf-1) 

   write(5,*)
endif

! writes matching parameters
write(5,*) '------------------------- MATCHING PROPERTIES ----------------------------'
write(5,80) 'Target frequency            : ',  targ_fr, 'Bandwidth for search        : ', band_wid 

write(5,*)

if (modality_flag /= 0) then
   ! writes scattering parameters in COMMON for all the stations
   write(5,*) '------------------------ SCATTERING PROPERTIES ---------------------------'

   write(5,80) 'Source density              : ',    srcR, 'Source S-wave speed         : ',  srcVs
   write(5,80) 'Maximum frequency           : ',    fmax, 'Start time for coda waves   : ',     t0
   write(5,81) 'Scattering wavelets         : ',   nscat, 'Number of pts coda envelope : ',  ncoda
   write(5,81) 'Scattering fft points       : ',  nfcoda, 'Seed number for parameters  : ', s_seed  
   write(5,80) 'Coda tolerance              : ',    merr
   write(5,80) 'Absorption coefficient      : ',abscoeff, 'Scattering coefficient      : ',scatcoeff
   write(5,80) 'Q0 value                    : ',       Q, 'Q frequency decay factor    : ',   fdec
   write(5,80) 'High-pass corner frequency  : ',   hpass, 'Bandwidth for filter        : ',  trans

   write(5,*)
endif

! writes time-series properties
write(5,*) '------------------------ TIME-SERIES PROPERTIES --------------------------'
!write(5,81) 'Npts in time-series         : ',  npts
!write(5,80) 'Time length                 : ',lf_len,'Time-step                   : ',lf_len/(npts-1) 
write(5,81) 'Npts in time-series         : ',  d_npts
write(5,80) 'Time length                 : ',d_npts*tmp_dt,'Time-step                   : ',tmp_dt 

write(5,*)

! initialize stations section (continues on "write_log" subroutine)
write(5,*) '------------------------- STATIONS PROPERTIES ----------------------------'

write(5,81) 'Total number of stations    : ', n_stat

write(5,*)

END SUBROUTINE log_init

!===================================================================================================

SUBROUTINE raytracing_input
! ----------------------------------------------------------------------------
!
! Description:
!
!   Prepares the input files needed for J.Hole's 3D-ray-tracing program. 
!   These files, written to disk, are:
!            vel1d_P.in, vel1d_S.in : velocity model for the ray-tracer
!            ray3d_P.in, ray3d_S.in : parameter files for ray-tracer
!            PStimes.in             : paramter file to extract the
!                                     P-S-times for selected stations
!   {i.e. ray-tracing will be performed for both P and S waves}
!   Moreover reads and assignes some stations and source physical properties
!   (siteVs, siteR, srcVs, srcR, kappa).
!
! Dependencies:
!
!   Internal subroutine vel1D and vel3D
!
! Notes:
!
!   Kappa, Vs and rho values may be specified here for each receiver.
!
! Warning:
!
!   Some output files (input for raytracing) will be deleted after their use
!
! Authors: W. Imperatori, M. Mai
!
! Modified: January 2009 (v1.3)
!
! Updated: March 2013 (v1.4.2)
!   Add reading Qp and Qs from a new formatting 1D velocity model file.
!   The depth, vp, vs, rh, Qp and Qs are used by source.f90 and/or composition.f90.
!   The depth, vp, vs and rh became global parameters.
! 
!   Add routines, compute layer thickness, compute average Vs & rho over fault,
!   and compute quarter-wavelength impedance effects (Boore and Joyner, BSSA 1997).
!
! Updated: December 2014 (v1.5.5.2)
!   Add dpt=0 for
!   !compute average Vs, rho over fault

use def_kind; use geometry; use flags; use io_file
use matching, only: match_fr; use scattering
use source_receiver;          use vel_model

implicit none

! counters, indexes, dummies 
integer(kind=i_single)                      :: i,ierr
real(kind=r_scat)                           :: tmp_siteR,tmp_siteVs
real(kind=r_single),allocatable,dimension(:):: fre
real(kind=r_single)                         :: xs,ys,vs_ave,d_ave,dpt
real(kind=r_single)                         :: stt,dl

! -----------------------------------------------------------------------------


tmp_siteR=0.; tmp_siteVs=0.  !initialize variables (useful if 3D model is used)

! open & read formatted 1D VELOCITY MODEL FILE (two header lines,thickness,vp,vs,rho)      
if (vel_flag == 1) then

   open(1,file=trim(vel_file),status='old',iostat=ierr)
   
   ! check-point for file opening
   if (ierr /= 0) call error_handling(1,trim(vel_file),'RAYTRACING_INPUT (io.f90)')

   ! find number of layers   
   read(1,*); read(1,*)
   n_lay=0 
   do
      read(1,*,iostat=ierr)
      if (ierr == -1) exit
      n_lay=n_lay+1
   enddo
                                      
   rewind(1); read(1,*); read(1,*)
 
   ! allocate arrays 
   if(.not.allocated(vp)) then
      allocate(thk(n_lay),fre(0:n_lay),depth(n_lay),vp(n_lay),vs(n_lay),rh(n_lay)) ! fix for fre
      !allocate(thk(n_lay),fre(n_lay),depth(n_lay),vp(n_lay),vs(n_lay),rh(n_lay))
      allocate(Qp(n_lay),Qs(n_lay))
   endif

   ! read physical properties
   do i=1,n_lay
      read(1,*) depth(i),vp(i),vs(i),rh(i),Qp(i),Qs(i)
   enddo   
   close(1)

!compute layer thicknesses
   do i=2,n_lay
      thk(i-1)=depth(i)-depth(i-1)
   enddo   

!compute average Vs, rho over fault
        vs_ave=0.
        d_ave=0.
        dpt=0.    ! add (v1552)
        do i=1,n_lay
            vs_ave=vs_ave+vs(i)*thk(i)
            d_ave=d_ave+rh(i)*thk(i)
            dpt=dpt+thk(i)
      enddo
      d_avef=d_ave/dpt
      vs_avef=vs_ave/dpt

!compute quarter-wavelength impedance effects (Boore and Joyner, BSSA 1997)
        vs_ave=0.
        d_ave=0.
        stt=0.
        dpt=0.
        do i=1,n_lay
          stt=stt+thk(i)/vs(i)
          fre(i)=1./(4.*stt)
          if (fre(i).gt.1.0) then
            vs_ave=vs_ave+vs(i)*thk(i)
            d_ave=d_ave+rh(i)*thk(i)
            dpt=dpt+thk(i)
         else
           dl=thk(i)*(fre(i-1)-1.0)/(fre(i-1)-fre(i))
           dpt=dpt+dl
           d_ave=d_ave+rh(i)*dl
           vs_ave=vs_ave+vs(i)*dl
           goto 10
        endif
      enddo

10    continue
      d_ave=d_ave/dpt
      vs_ave=vs_ave/dpt
      vsd_ave=vs_ave*d_ave

   ! calculate some scattering parameters for the 1D-model case (i.e. source and receiver Vs, 
   ! rho). Receiver Vs and rho are temporary set as constants (may be changed afterwards).
   ! Source values are referred to the hypocenter location. 
   tmp_siteR=rh(1)
   tmp_siteVs=vs(1)    !shallowest layer values
   do i=2,n_lay
      if ( (hyp_z >= depth(i-1)) .and. (hyp_z < depth(i)) ) then
         srcR=rh(i-1)
         srcVs=vs(i-1)
      endif
   enddo

endif
      
! open and read STATIONS FILE
open(1,file=trim(stat_file),status='old',iostat=ierr)

! check-point for file opening
if (ierr /= 0) call error_handling(1,trim(stat_file),'RAYTRACING_INPUT (io.f90)')

! read format and file extensions
read(1,*); read(1,*); read(1,*); read(1,*); read(1,'(a256)') lf_in_dir 
read(1,*); read(1,*)
read(1,*) lf_kind_flag                                                
read(1,*); read(1,*)

if (trim(lf_kind_flag) == 'RGF' .or. trim(lf_kind_flag) == 'CMP') then
   read(1,*) lf_x,lf_y,lf_z                                           
else if (trim(lf_kind_flag) == 'BIN') then
   read(1,*) lf_bin_file 
else if (trim(lf_kind_flag) == '3SF') then
   read(1,*) lf_x    
endif   

read(1,*); read(1,*) 

! find number of stations
n_stat=0 
do
   read(1,*,iostat=ierr)
   if (ierr == -1) exit
   n_stat=n_stat+1
enddo
       
rewind(1)       
                                                                        
do i=1,13
   read(1,*)
enddo   
  
! allocate some arrays for station-related quantities  
if (.not.allocated(xp)) allocate(xp(n_stat),yp(n_stat),stat_name(n_stat)) 
if (.not.allocated(siteR)) allocate(siteR(n_stat),siteVs(n_stat),kappa(n_stat))
! this is allocated here just for practical reasons
if (.not.allocated(match_fr)) allocate(match_fr(3,n_stat))
 
 
! Vs, rho and kappa may be specified here for each station. If [-1] is given,
! default values will be used (default kappa will be assigned later on)
do i=1,n_stat
   read(1,*) xp(i),yp(i),stat_name(i),siteVs(i),siteR(i),kappa(i)
   
   ! For rho and Vs, if no values are specified ([-1]), assign shallowest layer values (default)
   if (siteVs(i) == -1) siteVs(i)=tmp_siteVs
   if (siteR(i) == -1) siteR(i)=tmp_siteR

enddo

close(1)   !close stations-file


! open and read optional 2nd stations file
if (modality_flag == 0) then

   ! deallocate some arrays not needed in this modality
   if (allocated(siteR)) deallocate(siteVs,siteR,kappa)

   ! open and read STATIONS FILE
   open(1,file=trim(opt_stat_file),status='old',iostat=ierr)

   ! check-point for file opening
   if (ierr /= 0) call error_handling(1,trim(opt_stat_file),'RAYTRACING_INPUT (io.f90)')

   ! read format and file extensions
   read(1,*); read(1,*); read(1,*); read(1,*); read(1,'(a256)') opt_dir
   read(1,*); read(1,*)
   read(1,*) hf_kind_flag
   read(1,*); read(1,*)

   if (trim(hf_kind_flag) == 'RGF' .or. trim(hf_kind_flag) == 'CMP') then
      read(1,*) hf_x,hf_y,hf_z
   else if (trim(hf_kind_flag) == 'BIN') then
      read(1,*) hf_bin_file 
   else if (trim(hf_kind_flag) == '3SF') then
      read(1,*) hf_x    
   endif   

   read(1,*); read(1,*) 
                                
   ! allocate array for already computed HF seismograms  
   if (.not.allocated(opt_stat_name)) allocate(opt_stat_name(n_stat))
 
   ! Assign optional stations name
   do i=1,n_stat
      read(1,*) opt_stat_name(i)
   enddo

   close(1)   !close stations-file

endif

! open & read extended source file (if any)
if (ext_flag == 1 .and. modality_flag /=0) then
   open(1,file=trim(ext_file),status='old',iostat=ierr)
   
   ! check-point for file opening
   if (ierr /= 0) call error_handling(1,trim(ext_file),'raytracing_input (io.f90)')
   
   read(1,*) n_cell   !first line specifies number of subfaults
   
   ! allocate extended fault arrays
   if (.not.allocated(x_cell)) allocate(x_cell(n_cell),y_cell(n_cell),z_cell(n_cell))

   do i=1,n_cell         !row-wise ordered file (i.e. along strike)
      read(1,*) x_cell(i),y_cell(i),z_cell(i)
   enddo
   close(1)
endif


! due the formatting issues in the ray-tracing routine, it is most convenient
! to shift the HORIZONTAL coordinates of stations and "shot-point" (hypocenter) such
! that the system starts at the origin, x = 0.0, y = 0.0 {this means one extreme 
! corner having coordinates [0 0 0]}
xs=x_init; ys=y_init
x_init=x_init-xs; y_init=y_init-ys    !i.e. set all *_init equal to zero
x_far=x_far-xs; y_far=y_far-ys        !z-coordinate doesn't change
hyp_x=hyp_x-xs; hyp_y=hyp_y-ys

! shift subfaults (if necessary)
if (ext_flag == 1 .and. modality_flag /=0) then
   do i=1,n_cell
      x_cell(i)=x_cell(i)-xs; y_cell(i)=y_cell(i)-ys
   enddo
endif

! shift stations
do i=1,n_stat
   xp(i)=xp(i)-xs
   yp(i)=yp(i)-ys
enddo

! compute grid size (*_init could be also removed since equal to 0) 
nx=nint((x_far-x_init)/h_step)+1; ny=nint((y_far-y_init)/h_step)+1
nz=nint((z_far-z_init)/h_step)+1 

! pack a couple of arrays to be passed as arguments to the raytracing routine
grid = (/nx,ny,nz/)
hypo = (/hyp_x,hyp_y,hyp_z/)


! create 3D velocity model from 1D one. Again, default names - do not change!
if (vel_flag == 1) then
   call vel1d('vel1d_P.log','vel3d_P.bin',vp)
   call vel1d('vel1d_S.log','vel3d_S.bin',vs)
endif

! read 3D structure file and split it into two parts (Vp and Vs)
if (vel_flag == 0) call vel3d


!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

CONTAINS


   SUBROUTINE vel1d(file1d,file3d,vinp)
   !--------------------------------------------------------------------------
   !
   ! Description:
   !
   !   Converting a set of v(z) datapoints (i.e. 1D model) to an equally sampled
   !   3D velocity file appropriate for use with the ray-tracing code.
   !
   ! Dependencies:
   !
   !   None
   !
   ! Notes:
   !
   !   Log-file with resampled 1D structure available only if verbose_flag 
   !   is  set to 'on'.
   !   Original input structure is resampled LINEARLY along depth according to
   !   ray-tracing space-step. A Compsyn-like input file is assumed.
   !
   ! Authors: W. Imperatori, M. Mai
   !
   ! Modified: January 2009 (v1.3)
   !
   ! Updated: December 2014 (v1.5.5.2)
   !   Clear to open file1d when verbose_flag = 'on'.
   !

   use interfaces, only: polint

   implicit none

   ! names of log 1D file and output 3D file
   character(len=11),intent(in)                :: file1d,file3d
   ! velocity in input
   real(kind=r_single),intent(in),dimension(:) :: vinp
   ! indexes,flag 
   integer(kind=i_single)                      :: i,j,k,ierr,rl      
   ! array for velocity
   real(kind=r_single),dimension(nz)           :: v
   ! variable for depth
   real(kind=r_single)                         :: z_grid
   
   !--------------------------------------------------------------------------

   ! open log-file for structure 
   !if (trim(verbose_flag) == 'on') open(4,file=file1d,status='unknown',iostat=ierr)

   ! check-point for file opening
   !if (ierr /= 0) call error_handling(1,file1d,'vel1d (io.f90)')
   
   ! open log-file for structure (v1552) 
   if (trim(verbose_flag) == 'on') then
      open(4,file=file1d,status='unknown',iostat=ierr)

      ! check-point for file opening
      if (ierr /= 0) call error_handling(1,file1d,'vel1d (io.f90)')
   endif

   ! loop over z-nodes   
   do i=1,nz
   
      z_grid = h_step * (i-1) + z_init  !compute depth of i-th node
      
      if (z_grid == depth(1)) then
         v(1) = vinp(1)
      else if (z_grid >= depth(n_lay)) then
         v(i) = vinp(n_lay)
      else
         do j=1,n_lay-1
            if (z_grid > depth(j) .and. z_grid < depth(j+1)) then
               call polint((/depth(j),depth(j+1)/),(/vinp(j),vinp(j+1)/),z_grid,v(i))
            endif
            if (z_grid == depth(j))  v(i) = vinp(j)
         enddo           
      endif
    
   ! write new structure to log-file   
   if (trim(verbose_flag) == 'on') write(4,*) z_grid,v(i)
           
   enddo        
      
   ! evaluate record length for binary file
   inquire(iolength=rl) (v(1),i=1,nx)              
   print*,'rl= ',rl
   
   ! open output file 
   open(5,file=file3d,form='unformatted',access='direct',recl=rl,status='unknown',iostat=ierr)
   
   ! check-point for file opening
   if (ierr /= 0) call error_handling(1,file3d,'vel1d (io.f90)')
   
   ! write derived 3D model to output binary file
   do k=1,nz
      do j=1,ny
         write (5,rec=(k-1)*ny+j) (v(k),i=1,nx)   
      enddo
   enddo

   ! close output files 
   if (trim(verbose_flag) == 'on') close(4); close(5)


   END SUBROUTINE vel1d 

!<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-
   
   SUBROUTINE vel3d
   !--------------------------------------------------------------------------------
   ! 
   ! Description:
   !
   !   Re-arrange a specified 3D structure model for the ray-tracing code and split
   !   it into 2 output files (P-S waves velocity -> 'vel3d_P/S.out').
   !
   ! Dependencies:
   !
   !   None
   ! 
   ! Notes:
   !   The input binary file MUST have the same dimensions (i.e. same nx,ny,nz) as
   !   those of the input file for the raytracer (there is no check for this inside
   !   the code. The STANDARD input format is Vp (km/s), Vs (km/s), Rho (g/cm^3); 
   !   starting point: llc,surface. 
   !
   !   Vs and rho for receivers may be newly assigned here.
   !
   ! Author: W. Imperatori
   !
   ! Modified: January 2009 (v1.3)
   !
   
   use interfaces, only: poly_interp
   
   implicit none

   ! indexes,flag 
   integer(kind=i_single)                          :: i,j,k,ierr
   ! 3D arrays for parameters
   real(kind=r_single),allocatable,dimension(:,:,:):: vp_3d,vs_3d,rh_3d
   ! stations and source 
   integer(kind=i_single)                          :: src_x,src_y,src_z
   ! grid-indexes for neighbouring points
   integer(kind=i_single),dimension(2)             :: x_pos,y_pos
   ! record length
   integer(kind=i_single)                          :: rl,nrecvs
   !
   real(kind=r_single)                      :: vpmin,vpmax,vsmin,vsmax,rhmin,rhmax

   !--------------------------------------------------------------------------------
   
   ! allocate arrays for 3D structure
   if (.not.allocated(vp_3d)) allocate(vp_3d(nx,ny,nz),vs_3d(nx,ny,nz),rh_3d(nx,ny,nz))

   ! check record length of direct-access file
   inquire(iolength=rl) (vp_3d(1,1,1),i=1,3)                     
   print*,'rl= ',rl

   ! open input 3D velocity model binary file
   open(7,file=trim(vel_file),status='old',access='direct',recl=rl,iostat=ierr)
   
   ! check-point for file opening
   if (ierr /= 0) call error_handling(1,trim(vel_file),'vel3d (io.f90)')

   !STANDARD input format  (3D structure model)
   nrecvs=0
   
   print*,'reading media file'
   !open(177,file='vs3d',access='direct',recl=1)
   open(177,file='vs3d',access='direct',recl=4)
   print*,nx,ny,nz
   do k=1,nz
   vpmax=0.
   vpmin=100000.
   vsmax=0.
   vsmin=100000.
   rhmax=0.
   rhmin=100000.
      do j=1,ny
         do i=1,nx
            read(7,rec=(k-1)*ny*nx+(j-1)*nx+i) vp_3d(i,j,k),vs_3d(i,j,k),rh_3d(i,j,k)            
            if (vp_3d(i,j,k).gt.vpmax) vpmax=vp_3d(i,j,k)
            if (vp_3d(i,j,k).lt.vpmin) vpmin=vp_3d(i,j,k)
            if (vs_3d(i,j,k).gt.vsmax) vsmax=vs_3d(i,j,k)
            if (vs_3d(i,j,k).lt.vsmin) vsmin=vs_3d(i,j,k)
            if (rh_3d(i,j,k).gt.rhmax) rhmax=rh_3d(i,j,k)
            if (rh_3d(i,j,k).lt.rhmin) rhmin=rh_3d(i,j,k)
            nrecvs=nrecvs+1
            if (k.eq.1) write(177,rec=nrecvs) rh_3d(i,j,k)
         enddo
      enddo
   print*,k,vpmin,vpmax,vsmin,vsmax,rhmin,rhmax
   enddo 
   stop

   close(7)
   
   ! evaluate needed record length to write binary file
   inquire(iolength=rl) (vp_3d(i,1,1),i=1,nx)    

   ! open the output files and check 
   open(8,file='vel3d_P.bin',form='unformatted',access='direct',recl=rl,status='unknown', &
        iostat=ierr)
        
   if (ierr /= 0) call error_handling(1,'vel3d_P.bin','vel3d (io.f90)') 
   
   open(9,file='vel3d_S.bin',form='unformatted',access='direct',recl=rl,status='unknown', &
        iostat=ierr)
        
   if (ierr /= 0) call error_handling(1,'vel3d_S.bin','vel3d (io.f90)')      

   !STANDARD raytracer input format (P and S wave-speed)
   do k=1,nz
      do j=1,ny
         write(8,rec=(k-1)*ny+j) (vp_3d(i,j,k),i=1,nx)
         write(9,rec=(k-1)*ny+j) (vs_3d(i,j,k),i=1,nx)
      enddo
   enddo
                                                         
   close(8); close(9)

   ! assign site Vs and rho using bilinear interpolation (if not already assigned via stations file)
   do i=1,n_stat
      
      ! find station's neighbouring points ([2x2] for bilinear interpolation)
      x_pos = (/ floor( xp(i) / h_step ) + 1, floor( xp(i) / h_step ) + 2 /)  
      y_pos = (/ floor( yp(i) / h_step ) + 1, floor( yp(i) / h_step ) + 2 /)
      
      ! assign rho and Vs through bilinear interpolation
      if (siteR(i) <= 0)   call  &
      poly_interp(h_step*(x_pos-1),h_step*(y_pos-1),rh_3d(x_pos,y_pos,1),xp(i),yp(i),siteR(i))
      
      if (siteVs(i) <= 0)  call  &
      poly_interp(h_step*(x_pos-1),h_step*(y_pos-1),vs_3d(x_pos,y_pos,1),xp(i),yp(i),siteVs(i))
      
   enddo
   
   ! source (hypocenter) position in terms of grid points 
   src_x=nint(hyp_x/h_step)+1; src_y=nint(hyp_y/h_step)+1; src_z=nint(hyp_z/h_step)+1
   
   ! extract source parameters
   srcR=rh_3d(src_x,src_y,src_z); srcVs=vs_3d(src_x,src_y,src_z) 
                                                                    
   ! free memory
   if(allocated(vp_3d)) deallocate(vp_3d,vs_3d,rh_3d)

   END SUBROUTINE vel3d

END SUBROUTINE raytracing_input

!===================================================================================================

SUBROUTINE read_inputfile(in_file)
!-----------------------------------------------------------------
!
! Description:
!
!   External subroutine, reads main input file and stores values
!
! Dependencies:
!
!   None
!
! Notes:
!
!   Official units are: km, seconds
!   Allmost all parameters are passed by specific modules
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: March 2013 (v1.4.2)
!   Add reading srf file name from the last line of input file. 
!
! Updated: March 2014 (v1.5.4)
!   Target frequency and the bandwidth are removed from main input file.
!   The targ_fr is computed as a function of Mw.
!
! Updated: December 2016 (v1.6.2)
!   Target frequency is computed in source.f90 with Mw from srf file.
!
! Updated: February 2019 (v2.0)
!   Add reading Kemp_*.bin file name for the correlation.
!
! Updated: March 2019 (v2.1)
!   Add rake for alt computation in source.f90.
!
use def_kind; use earthquake; use flags; use geometry; use io_file
use matching; use stf_data, only: stf_name,srf_name

implicit none

! main input-file name
character(len=*),intent(in):: in_file
! flag 
integer(kind=i_single)     :: ierr
! temporary variable
character(len=10)          :: scratch

!----------------------------------------------------------------

open(1,file=in_file,status='old',iostat=ierr)    !open main input file

! check-point for file opening
if (ierr /= 0) call error_handling(1,in_file,'READ_INPUTFILE (io.f90)')

! read file
read(1,*); read(1,*) modality_flag               !flag for modality status
                                     
! error check
if (modality_flag < 0 .or. modality_flag > 2) then
   call error_handling(10,'modality flag','READ_INPUTFILE (io.f90)')
endif   

! temporary error, to be removed once isochrone approach is implemented
if (modality_flag == 2) call error_handling(1,'null','READ_INPUTFILE (io.f90)')

read(1,*); read(1,*) output_dir                  !output directory

read(1,*); read(1,*) vel_file                    !velocity file and its flag
ierr = index(trim(vel_file),'.bin')
if (ierr /= 0) then
   vel_flag = 0                !3D velocity file (BIN)
else 
   vel_flag = 1                !1D velocity file (ASCII)
endif   
                                                                   
read(1,*); read(1,*) stat_file                   !stations file

if (modality_flag == 0) then
   read(1,*); read(1,*) opt_stat_file            !optional 2nd stations file 
else
   read(1,*); read(1,*)
endif

read(1,*)
read(1,*) scratch
if (trim(scratch) == 'point') then
   ext_flag = 0                                  !fault type: point or extended (and its file)
else
   backspace(1)
   read(1,*) scratch, ext_file         
   ext_flag = 1
endif
                                                      
read(1,*); read(1,*) hyp_x,hyp_y,hyp_z           !hypocenter coordinates

read(1,*); read(1,*) x_init,y_init,z_init,h_step !\
                                                 !  raytracing box corners and grid spacing
read(1,*); read(1,*) x_far,y_far,z_far           !/

read(1,*); read(1,*) scat_file                   !scattering parameters file
!read(1,*); read(1,*) targ_fr,band_wid            !targert matching frequency, bandwidth for search
read(1,*); read(1,*) Mw                          !magnitude of event
read(1,*); read(1,*) mech                        !dominant source mechanism
read(1,*); read(1,*) stf_name                    !STF to be convolved with scatter-coda 

read(1,*); read(1,*) verbose_flag                !flag for verbose output
read(1,*); read(1,*) srf_name                    !SRF file name

read(1,*); read(1,*) rake                        ! read rake for alt factor

read(1,*); read(1,*) corr_file_inf               !correlation matrix file
read(1,*); read(1,*) corr_file_sp1, corr_file_sp2, corr_file_sp3

close(1)

END SUBROUTINE read_inputfile

!===================================================================================================

SUBROUTINE read_seis(station,freq_flag)
!----------------------------------------------------------------------
!
! Description:
!
!   External subroutine, read LF time-series in input. Several formats 
!   are supported (ASCII and binary)
!
! Dependencies:
!
!   None
!
! Formats:
!
!   .RGF -> Rob Graves' format (ASCII):
!   - two-line header: station-name  component
!                         npts          dt
!   - data, written as 6 entries for each line
!
!   .CMP -> Compsyn format (ASCII):
!   - 99 lines header, with dt and npts specified at lines 6 and 7
!
!   - data, one time-step for each line
!
!   .3SF -> Straight-Plain format (ASCII):
!   - no header lines
!
!   - t x y z, one time-step for each line 
!
!   .BIN -> Binary format:
!   - first line: npts, dt, <null>
!
!   -  x y z, one time-step for each line
!
! Notes:
!
!   RGF and CMP FORMATS ASSUME 3-COMPONENTS STORED IN DIFFERENT FILES.
!   3SF and BIN FORMATS ASSUME 3-COMPONENTS STORED IN THE SAME FILE.
!
!   In the binary file, stations have to be stored in the same order as 
!   that specified in the respective station file.
!
!   First time sample is always assumed to be equal to 0. 
!
!   ALL THESE FILES MUST CONTAIN VELOCITY TIME-SERIES IN M/S !
!
! Warning:
!
!   For every ASCII case, file root-name MUST be equal to that specified 
!   inside the stations input file.
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: April 2014 (v1.5.5)
!    Store low friquency files that the lf-length is adjusted by
!    tmp_lf_len.
!    This is only for case('3SF')
!
! Updated: December 2014 (v1.5.5.2)
!    Store low friquency files that the lf-length is adjusted by
!    tmp_lf_len for case('BIN').
!
! Updated: February 2015 (v1.5.5.3)
!    Make clear which parameters are used for each low frequency
!    file format.  
!
! Updated: July 2015 (v1.6.1)
!    Change character length from 4 to 10 for cp_x, cp_y and cp_z
!
! Updated: February 2019 (v2.0)
!   Change to decide lf_length based of tmp_lf_len (102.3750 sec), and dt using v_npts
!   only for 3SF and BIN cases for now,
!   to keep the same dt even though from different/changed lf_length in an event.
!
use def_kind; use flags; use io_file
use source_receiver, only: n_stat,stat_name,opt_stat_name
use waveform; use tmp_para

implicit none

! actual station number
integer(kind=i_single),intent(in)             :: station
! input flag
character(len=2),intent(in)                   :: freq_flag
! local flag
character(len=4)                              :: loc_flag
! local filenames
character(len=256)                            :: name,name_bin
! local time-series
real(kind=r_single),allocatable,dimension(:,:):: timeseries,tmp_timeseries
! local time-series length (seconds)
real(kind=r_single)                           :: ts_len
! local time-series number of points
integer(kind=i_single)                        :: ts_npts
! local file extensions
character(len=10)                             :: cp_x,cp_y,cp_z
! dummies and counters 
character(len=256)                            :: x_name,y_name,z_name,loc_dir,dummy_line
integer(kind=i_single)                        :: k,j,remainder,read_nline
integer(kind=i_single)                        :: pts_count,ierr,i,n,head_lines,wts_npts
real(kind=r_single)                           :: dummy,ts_dt,cm2m,wts_len
!real(kind=r_single),parameter                 :: tmp_lf_len=102.3750

!----------------------------------------------------------------------

! Initial condition test to check wheter LF or optional HF are going to be opened
!if (freq_flag == 'LF') then
!   name = stat_name(station)
!   name_bin = lf_bin_file
!   loc_flag = lf_kind_flag
!   cp_x = lf_x; cp_y = lf_y; cp_z = lf_z
!   loc_dir = lf_in_dir
!else if (freq_flag == 'HF') then
!   name = opt_stat_name(station)
!   name_bin = hf_bin_file
!   loc_flag = hf_kind_flag  
!   cp_x = hf_x; cp_y = hf_y; cp_z = hf_z
!   loc_dir = opt_dir
!endif

! Make clear which parameters are used for the file format
if (freq_flag == 'LF') then
   loc_flag = lf_kind_flag
   loc_dir = lf_in_dir
   if (trim(lf_kind_flag) == 'BIN') then
      name_bin = lf_bin_file
   else
      name = stat_name(station)
      if (trim(lf_kind_flag) == '3SF') then
         cp_x = lf_x
      else
         cp_x = lf_x; cp_y = lf_y; cp_z = lf_z
      endif
   endif
elseif (freq_flag == 'HF') then
   name = opt_stat_name(station)
   name_bin = hf_bin_file
   loc_flag = hf_kind_flag
   cp_x = hf_x; cp_y = hf_y; cp_z = hf_z
   loc_dir = opt_dir
endif

! prepare to read input files according to format)
select case(loc_flag)
case('RGF')               !Rob Graves format
            
   ! combining filenames
   x_name=trim(loc_dir)//trim(name)//cp_x; y_name=trim(loc_dir)//trim(name)//cp_y
   z_name=trim(loc_dir)//trim(name)//cp_z

   ! opening files (one for each component) and check points
   open(1,file=trim(x_name),status='old',iostat=ierr)
   if (ierr /= 0) call error_handling(1,trim(x_name),'READ_SEIS (io.f90)')
   
   open(2,file=trim(y_name),status='old',iostat=ierr)
   if (ierr /= 0) call error_handling(1,trim(y_name),'READ_SEIS (io.f90)')
   
   open(3,file=trim(z_name),status='old',iostat=ierr)
   if (ierr /= 0) call error_handling(1,trim(z_name),'READ_SEIS (io.f90)')
   
   ! read seismograms (npts assumed equal for all components and stations)
   read(1,*); read(1,*) ts_npts, ts_dt
   read(2,*); read(2,*)
   read(3,*); read(3,*) 
   
   ! check for a remainder
   remainder=mod(ts_npts,6)
   read_nline=ts_npts/6
   if (remainder /= 0) then
      read_nline=floor(ts_npts/6.)
      ts_npts = read_nline*6
      call warning_handling(3,'null','READ_SEIS (io.f90)')
   endif
   
   ! allocate array for LF waveforms
   allocate(timeseries(ts_npts,3))

   ts_len=(ts_npts-1)*ts_dt              

   ! reading files...
   do k = 1,read_nline       
      j = 6*k
      read(1,*) timeseries(j-5,1),timeseries(j-4,1),timeseries(j-3,1),  &
                timeseries(j-2,1),timeseries(j-1,1),timeseries(j,1)
      read(2,*) timeseries(j-5,2),timeseries(j-4,2),timeseries(j-3,2),  &
                timeseries(j-2,2),timeseries(j-1,2),timeseries(j,2)
      read(3,*) timeseries(j-5,3),timeseries(j-4,3),timeseries(j-3,3),  &
                timeseries(j-2,3),timeseries(j-1,3),timeseries(j,3)
   enddo

   close(1); close(2); close(3)   !close all files  

case('CMP')               !COMPSYN FORMAT   
         
   ! combining filenames
   x_name=trim(loc_dir)//trim(name)//cp_x; y_name=trim(loc_dir)//trim(name)//cp_y
   z_name=trim(loc_dir)//trim(name)//cp_z

   ! opening files (one for each component) and check points
   open(1,file=trim(x_name),status='old',iostat=ierr)
   if (ierr /= 0) call error_handling(1,trim(x_name),'READ_SEIS (io.f90)')
   
   open(2,file=trim(y_name),status='old',iostat=ierr)
   if (ierr /= 0) call error_handling(1,trim(y_name),'READ_SEIS (io.f90)')
   
   open(3,file=trim(z_name),status='old',iostat=ierr)
   if (ierr /= 0) call error_handling(1,trim(z_name),'READ_SEIS (io.f90)')
   
   ! read seismograms (npts assumed equal for all components and stations)
   do i=1,5
      read(1,*); read(2,*); read(3,*) 
   enddo
   
   read(1,*) ts_dt; read(2,*); read(3,*)
   read(1,*) ts_npts; read(2,*); read(3,*)       
   
   ! allocate array for LF waveforms
   allocate(timeseries(ts_npts,3))

   ts_len=(ts_npts-1)*ts_dt              

   do k=1,(100-7)
      read(1,*); read(2,*); read(3,*)
   enddo
   
   do k=1,ts_npts
      read(1,*) timeseries(k,1); read(2,*) timeseries(k,2); read(3,*) timeseries(k,3)
   enddo

   close(1); close(2); close(3)   !close all files  

case('3SF')           !single ASCII files with 4 columns (time vector and 3 components)
   cm2m = 0.01

   ! merging filenames
   x_name=trim(loc_dir)//trim(name)//cp_x

   ! opening files and check points
   open(1,file=trim(x_name),status='old',iostat=ierr)
   if (ierr /= 0) call error_handling(1,trim(x_name),'READ_SEIS (io.f90)')

   ! Find length of header
   head_lines = 0
   do
      read (1,*) dummy_line
      !This line checks to see if it starts with a '#' or '%'; if so, ignore
      if ((dummy_line(1:1) .ne. '#') .AND. (dummy_line(1:1) .ne. '%')) exit
      head_lines=head_lines+1
   enddo
   rewind(1)

   ! read past header
   do j=1,head_lines
      read(1,*)
   enddo
   
   ! find number of points of time-series
   pts_count=0    
   do
      read(1,*,iostat=ierr)
      if (ierr == -1) exit
      pts_count=pts_count+1
   enddo   
   rewind(1)

   ! assign number or points
   wts_npts=pts_count 

   ! allocate array for LF waveforms
   if (.not.allocated(tmp_timeseries)) allocate (tmp_timeseries(wts_npts,3))

   ! read past header
   do j=1,head_lines
      read(1,*)
   enddo

   ! reading file... compute also time-step 
   do j=1,wts_npts       
      read(1,*) dummy,tmp_timeseries(j,1),tmp_timeseries(j,2),tmp_timeseries(j,3)
      tmp_timeseries(j,1) = cm2m*tmp_timeseries(j,1)
      tmp_timeseries(j,2) = cm2m*tmp_timeseries(j,2)
      tmp_timeseries(j,3) = cm2m*tmp_timeseries(j,3)
      if (j == 2) ts_dt=dummy
   enddo

   wts_len=(wts_npts-1)*ts_dt  !compute time-series length
   real_lf_len = wts_len  ! record real LF length

   ! adjust ts_npts and ts_len
   print*,'wts_npts,ts_dt,wts_len',wts_npts,ts_dt,wts_len
   n = ceiling(wts_len/tmp_lf_len)
   if (mod(wts_len, tmp_lf_len) == 0) then
      ts_npts=wts_npts
      ts_len=wts_len
   elseif (mod(wts_len, tmp_lf_len) .gt. 0) then
      ! n = ceiling(wts_len/tmp_lf_len)
      ts_len=tmp_lf_len*n
      ts_npts=(tmp_lf_len*n/ts_dt) + 1
   endif

   ! compute virtual npts for BBsynthetics used later on to keep the same dt
   ! even LF length changes in an event
   v_npts = (tmp_npts-1)*n+1

   allocate(timeseries(ts_npts,3))
   timeseries(:,:)=0
   timeseries(1:wts_npts,:)=tmp_timeseries(1:wts_npts,:)

   deallocate(tmp_timeseries)

   close(1)   !close file

case('BIN')              !BINARY file format
   
   ! merging filenames
   x_name=trim(loc_dir)//trim(name_bin)
   
   ! open binary file
   !open(1,file=x_name,status='old',access='direct',form='unformatted',recl=3)
   open(1,file=x_name,status='old',access='direct',form='unformatted',recl=12)
   
   ! first line is reserved to: number of points, time-step, dummy
   ! Note: these quantities are assumed to be the same for all the stations
   ! read(1,rec=1) ts_npts, ts_dt, dummy
   read(1,rec=1) wts_npts, ts_dt, dummy
   
   ! allocate array for LF waveforms
   ! allocate(timeseries(ts_npts,3))
   allocate(tmp_timeseries(wts_npts,3))
   
   ! ts_len=(ts_npts-1)*ts_dt  !compute time-series length
   wts_len=(wts_npts-1)*ts_dt  !compute time-series length
   
   ! reading file... (3 components at the same record line -- each record line is a time-step)
   ! do k=2,ts_npts+1  
   !    read(1,rec= (station-1) * ts_npts + k) timeseries(k-1,1),timeseries(k-1,2),timeseries(k-1,3)
   ! enddo
   do k=2,wts_npts+1
      read(1,rec= (station-1) * wts_npts + k) tmp_timeseries(k-1,1), &
                                tmp_timeseries(k-1,2),tmp_timeseries(k-1,3)
   enddo

   ! adjust ts_npts and ts_len
   print*,'wts_npts,ts_dt,wts_len in Bin file',wts_npts,ts_dt,wts_len
   n = ceiling(wts_len/tmp_lf_len)
   if (mod(wts_len, tmp_lf_len) == 0) then
      ts_npts=wts_npts
      ts_len=wts_len
   elseif (mod(wts_len, tmp_lf_len) .gt. 0) then
      ! n = ceiling(wts_len/tmp_lf_len)
      ts_len=tmp_lf_len*n
      ts_npts=(tmp_lf_len*n/ts_dt) + 1
   endif

   ! compute virtual npts for BBsynthetics used later on
   v_npts = (tmp_npts-1)*n+1

   allocate(timeseries(ts_npts,3))
   timeseries(:,:)=0
   timeseries(1:wts_npts,:)=tmp_timeseries(1:wts_npts,:)

   deallocate(tmp_timeseries)
   
   close(1)  !close file

end select

! assign values to correct variables
if (freq_flag == 'LF') then
   lf_len = ts_len
   lf_npts = ts_npts
   lf_dt = ts_dt

   print*,'lf_len,lf_npts,lf_dt,v_npts',lf_len,lf_npts,lf_dt,v_npts

   if (.not.allocated(lf_seis)) allocate (lf_seis(lf_npts,3))

   lf_seis(:,:) = timeseries

else if (freq_flag == 'HF') then
   !hf_len = ts_len
   !hf_npts = ts_npts
   !if (.not.allocated(hf_seis)) allocate (hf_seis(hf_npts,3,n_stat))
   !hf_seis(:,:,station) = timeseries
endif   

deallocate(timeseries)

END SUBROUTINE read_seis

!===================================================================================================

SUBROUTINE write_disk(station,type_flag,in_arr)
!-------------------------------------------------------------------
! Description:
!
!   External subroutine, writes broad-band seismograms to disk in
!   the recently discussed format, with 4 columns and some header 
!   informations (8-lines headers).
!   - two-line header: station-name  component
!                         npts          dt
!   - time, NS, EW, UD
!
! Dependencies:
!
!   None. Many variables passed by modules.
!
! Notes:
!
!   In addition to broad-band, source-time-function, original
!   coda and convolved coda outputs are available. Recently, 
!   binary output (only for broad-band) has been added as well.
!
! Author: W. Imperatori
! 
! Modified: January 2009 (v1.3)
!
! Updated: October 2013 (v1.5.2)
!   Add output decimation factor, time_step.
!
! Updated: March 2014 (v1.5.4)
!   time_step is set in scattering.dat.
!
! Updated: April 2014 (v1.5.5)
!    tmp_npts and tmp_lf_len are in the module tmp_para.
!
! Updated: December 2014 (v1.5.5.2)
!    The d_npts is used for case('BIN').
!    Use d_npts and tmp_dt for the header line in ocd and ccd output files.
!
! Updated: February 2019 (v2.0)
!   Add to use v_npts for dt computation.

!
use def_kind; use flags; use io_file, only: output_dir
use scattering, only: npts,time_step; use source_receiver, only: stat_name
use stf_data, only: npts_stf,total; use waveform, only: lf_len,lf_npts,v_npts,real_lf_len
use tmp_para

implicit none

! actual station number
integer(kind=i_single),intent(in)            :: station
! input time-series
real(kind=r_single),dimension(:,:),intent(in):: in_arr
! kind of output requested
character(len=*),intent(in)                  :: type_flag
! suffix of output file
character(len=4)                             :: suffix
! output file name
character(len=256)                           :: out_name
! time-step and time values
real(kind=r_single)                          :: dt,time,m2cm
! counter
integer(kind=i_single)                       :: i,scratch

! add
!real(kind=r_single),parameter                 :: tmp_lf_len = 102.3750
!integer(kind=i_single),parameter              :: tmp_npts = 32768

! add for output decimation factor
!integer(kind=i_single),parameter              :: time_step = 1
! decimated npts for output
integer(kind=i_single)                       :: d_npts
! adjusted dt for output, by d_npts
real(kind=r_single)                          :: tmp_dt
!-------------------------------------------------------------------

! select kind of output
select case(type_flag)
case('hyb')
   if(lf_kind_flag == 'BIN') then
      suffix='.bin'
   else 
      suffix='.hyb'
   endif     
case('ccd')
   suffix='.ccd'     !scatterogram after convolution
case('ocd')
   suffix='.ocd'     !scatterogram before convolution
case('stf')
   suffix='.stf'
end select

! select time-step accordingly
select case(type_flag)      
case('stf')
   dt = total/(npts_stf-1)     
case default
   !dt = lf_len/(npts-1)
   dt = lf_len/(v_npts-1)
   print*,'lf_len,npts,dt,v_npts,real_lf_len in write_dist=',lf_len,npts,dt,v_npts,real_lf_len
   !!!d_npts=v_npts/time_step  !!d_npts=npts/time_step
   tmp_dt = dt*time_step
   d_npts = ceiling(real_lf_len/tmp_dt +1)
end select

! binary output (only for broad-band)
if (suffix == '.bin') then
   
   out_name=trim(output_dir)//'/BBhyb.bin'
   
   inquire(iolength=scratch) in_arr(1,1),in_arr(1,2),in_arr(1,3)
   if(station == 1) then
      open(1,file=out_name,status='unknown',access='direct',form='unformatted',recl=scratch)
      !write(1,rec=1) npts,dt,1.0
      write(1,rec=1) d_npts,tmp_dt,1.0
   else
      open(1,file=out_name,status='old',access='direct',form='unformatted',recl=scratch)
   endif      
  
   print*,'npts & d_npts in write_disk in .bin file',npts,d_npts
 
   !do i=2,npts+1
   do i=2,d_npts+1
      !write(1,rec=(station-1)*npts +i) in_arr(i-1,1),in_arr(i-1,2),in_arr(i-1,3)
      write(1,rec=(station-1)*d_npts +i) in_arr(i-1,1),in_arr(i-1,2),in_arr(i-1,3)
   enddo
   
! ASCII output   
else
   out_name=trim(output_dir)//'/BB.'//trim(stat_name(station))//suffix
   
   open(1,file=trim(out_name),form='formatted',status='unknown')
 
   ! write HEADERS (8 lines)
   write(1,102) '% -----------------------------------------------'

   select case(type_flag)
   case('hyb')
      write(1,102) '% synthetic broadband seismogram (Mai&Olsen 2008) '
   case('ccd')
      write(1,102) '% scatterogram after convolution (Mai&Olsen 2008) '
   case('ocd')
      write(1,102) '% scatterogram before convolution (Mai&Olsen 2008)'
   case('stf')
      write(1,102) '% source-time function (Mai&Olsen 2008)           '
   end select   

   write(1,103) '% N = 8 header lines'
   write(1,100) '% site: ',trim(stat_name(station))
   select case(type_flag)
   case default
      ! write(1,101) '% NPTS, DT: ',npts,dt
      write(1,101) '% NPTS, DT: ',d_npts,tmp_dt
   case('stf')
      write(1,101) '% NPTS, DT: ',npts_stf,dt
   case('ocd')
      write(1,101) '% NPTS, DT: ',npts,dt
      !!!write(1,101) '% NPTS, DT: ',d_npts,tmp_dt
   case('ccd')
      write(1,101) '% NPTS, DT: ',npts,dt
      !!!write(1,101) '% NPTS, DT: ',d_npts,tmp_dt
   end select

   write(1,104) '%'
   select case(type_flag)
   case default
      !write(1,102) '% time(s)      NS (m/s)      EW(m/s)     UP (m/s)'
      write(1,102) '% time(s)      NS (cm/s)      EW(cm/s)     UP (cm/s)'
   case('stf')
      write(1,105) '% time(s)     Slip-rate(m/s)'
   end select
   write(1,102) '% -----------------------------------------------'
   m2cm = 100.0
           
   ! write all the time-series into a single file
   select case(type_flag)
   case default
     ! do i=1,npts
     !    time = (i-1)*dt
     !    write(1,200) time,in_arr(i,1),in_arr(i,2),in_arr(i,3)
     ! enddo
      do i=1,d_npts
         time = (i-1)*tmp_dt
         write(1,200) time,m2cm*in_arr(i,1),m2cm*in_arr(i,2),m2cm*in_arr(i,3)
      enddo
   case('stf')
      do i=1,npts_stf
         write(1,201) in_arr(i,1),in_arr(i,2)
      enddo
   end select

endif   

! Formats
100   format(A8,A11)
101   format(A12,I6,F15.12)
102   format(A50)
103   format(A20)
104   format(A1)
105   format(A28)
200   format(F9.5,3E15.5)
201   format(F9.5,E15.5)

close(1)  !close output file

END SUBROUTINE write_disk

!===================================================================================================

SUBROUTINE write_log(station)
!-----------------------------------------------------------------------
!
! Description:
!
!   External subroutine, writes several station-specific parameters to a 
!   log-file, already initialized
!
! Dependencies:
!
!   None. Allmost all variables passed via modules
!
! Notes: 
!
!   Default log-file unit number is 5
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: December 2016 (v1.6.2)
!   Change output aveVs to loc_aveVs(station) in run.log
!
use def_kind; use flags; use geometry; use matching; use scattering
use source_receiver

implicit none

! Actual station number
integer(kind=i_single),intent(in):: station

!----------------------------------------------------------------------

! Formats
80 format(a30,f0.3,T42,a30,f0.3)
81 format(a30,i0,T42,a30,i0)
82 format(a30,f0.3,2x,f0.3,2x,f0.3)


write(5,*)  'Station number: ', station,'  [',trim(stat_name(station)),']'

if (modality_flag /= 0) write(5,80) 'Kappa value                 : ', kappa(station)
write(5,82) 'Matching frequencies [P N V]: ',match_fr(1,station),match_fr(2,station),  &
                                             match_fr(3,station)
if (modality_flag /= 0) then
   write(5,80) 'Site density                : ',siteR(station), 'Site S-wave speed           : ',&
                                                siteVs(station)
   write(5,80) 'Average S-wave speed        : ',loc_aveVs(station), 'Hypocenter distance         : ', &
                                                 sr_hypo(station)
   !write(5,80) 'Average S-wave speed        : ', aveVs, 'Hypocenter distance         : ',        &
   !                                              sr_hypo(station)
endif

write(5,80) 'P-wave traveltime           : ',time_p(station), 'S-wave traveltime           : ',  &
                                             time_s(station)

if (modality_flag == 0) then
   write(5,80) 'HF time-series arrival time : ', time_s(station)-tau(station),   &
               'Time-delay value            : ', tau(station)
endif

write(5,*) '--------------------------------------------------------------------------'

END SUBROUTINE write_log


!===================================================================================================
SUBROUTINE read_correlation_matrix_inf
!-----------------------------------------------------------------------
!
! Description:
!
!   Read spatical correlation matrices
!
! Author: N. Wang, December 2018
!
use def_kind
use io_file
use read_correlation_files

implicit none
! indexes
integer(kind=i_single)                      :: i,j,k,ierr
!-----------------------------------------------------------------------

!! Read input K, Cholesky factor of inter-frequency correlation matrix, lower triangular matrix
!print*,'nk,nks in io.f90,before',nk,nks
open(unit=3331,file=trim(corr_file_inf),access='direct',form='unformatted', recl=4,status='old')
read(3331,rec=1,iostat=ierr) nk
read(3331,rec=2,iostat=ierr) nks
print*,'nk nks in io.f90',nk, nks
if (.not.allocated(Kinf)) allocate(Kinf(nk,nk))
k = 2
do j = 1,nk
   do i = 1,nk
      k = k+1
      read(3331,rec=k,iostat=ierr) Kinf(i,j)
   enddo
enddo
close(unit=3331)
print*,'Kinf(1:10,1) in io.f90',Kinf(1:10,1)
print*,'Kinf(1:5,200) in io.f90',Kinf(1:5,200)
print*,'Kinf(1,1:10) in io.f90',Kinf(1,1:10)
END SUBROUTINE read_correlation_matrix_inf



!===================================================================================================
SUBROUTINE read_correlation_matrix_sp
!-----------------------------------------------------------------------
!
! Description:
!
!   Read spatical correlation matrices
!
! Author: N. Wang, December 2018
!
use def_kind
use read_correlation_files
use io_file
!!use source_receiver, only: n_stat

implicit none
! indexes
integer(kind=i_single)                      :: i,j,k,ierr
!-----------------------------------------------------------------------

!! Read input K1 K2 K3, Cholesky factor of frequency correlation matrix, lower triangular matrix
open(unit=3332,file=trim(corr_file_sp1),access='direct',form='unformatted', recl=4,status='old')
read(3332,rec=1,iostat=ierr) nk
read(3332,rec=2,iostat=ierr) nks
if (.not.allocated(Ksp1)) allocate(Ksp1(nk,nk))
k = 2
do j = 1,nk
   do i = 1,nk
      k = k+1
      read(3332,rec=k,iostat=ierr) Ksp1(i,j)
   enddo
enddo
close(unit=3332)

if (.not.allocated(Ksp2)) allocate(Ksp2(nk,nk))
open(unit=3333,file=trim(corr_file_sp2),access='direct',form='unformatted', recl=4*nk*nk,status='old')
read(3333,rec=1,iostat=ierr) ((Ksp2(i,j), i=1,nk), j=1,nk)
close(unit=3333)

if (.not.allocated(Ksp3)) allocate(Ksp3(nk,nk))
open(unit=3334,file=trim(corr_file_sp3),access='direct',form='unformatted', recl=4*nk*nk,status='old')
read(3334,rec=1,iostat=ierr) ((Ksp3(i,j), i=1,nk), j=1,nk)
close(unit=3334)

print*,'Ksp1(200,1:5)', Ksp1(200,1:5)
print*,'Ksp2(200,1:5)', Ksp2(200,1:5)
print*,'Ksp3(200,1:5)', Ksp3(200,1:5)
!!! Read input L1 L2, Cholesky factor of station correlation matrix, upper triangular matrix
!if (.not.allocated(L1)) allocate(L1(n_stat,n_stat))
!open(unit=3334,file='L1.bin',access='direct',form='unformatted', recl=4*n_stat*n_stat,status='old')
!read(3334,rec=1,iostat=ierr) ((L1(i,j), i=1,n_stat), j=1,n_stat)
!close(unit=3334)
!
!if (.not.allocated(L2)) allocate(L2(n_stat,n_stat))
!open(unit=3335,file='L2.bin',access='direct',form='unformatted', recl=4*n_stat*n_stat,status='old')
!read(3335,rec=1,iostat=ierr) ((L2(i,j), i=1,n_stat), j=1,n_stat)
!close(unit=3335)

END SUBROUTINE read_correlation_matrix_sp
