PROGRAM main_broadband
!----------------------------------------------------------------------------------------
!
! Description:
!
!   Code to compute broadband synthetics. It can operate in different modalities to 
!   be specified in the main input file.
!   Modality 0: it merges LF and HF waveforms previously computed. Merging process takes
!               into account possible time-delays, and is based on the approach of Mai
!               and Beroza, 2003.
!   Modality 1: it computes HF scatterograms and then merges them with previously 
!               computed LF waveforms. HF time-series are based on Zeng's multiple 
!               S-to-S scattering theory. Both point and extended source approximations
!               are implemented. Empircal Green functions approach is used for the 
!               extended source case. Some scattering parameters can be specified for
!               each site. 
!   Modality 2: isochrones section -> to be implemented!
!
! Dependencies:
!
!    Several external subroutines and functions. Many variables passed via modules. 
!
! Notes: 
!
!   Modality 0: time-delays are managed through a first-break detecting algorithm.
!               Input time-series must have same numbers of points in time and same
!               time-step.
!   Modality 1: default scattering parameters values can be overridden by the user
!               through a separate input file.
!               Input LF synthetics may be upsampled to avoid aliasing effects (according
!               to fmax). 
!   All modalities: time-series may be re-sampled into 2**N number of points.
!                   Output files are NOT consequently downsampled.
!                   Default log-file unit number for I/O is 5 (open in subroutine log_init
!                   in io.f90)
!
! References:
!
!   Mai and Olsen (2008). ??????
!
!   Mai and Beroza (2003). A hybrid method for calculating near-source, broad-band seismograms: 
!   application to strong motion prediction, PEPI, 137 183-199
!
!   Zeng, Anderson and Su (1995).  Subevent rake and random scattering effects in realistic 
!   strong ground motion simulation, Geophy. Res. Lett., 22 17-20
!
!   Hole and Zelt (1995). 3D finite-difference reflection traveltimes, Geophys. J. Int., 121 
!   427-434
!
!   Numerical Recipes in Fortran90, Cambridge University Press.
!
! Languages: Fortran95 and C
!
! Parallelization: OpenMP (for multicore CPUs)  
!
! Authors: W. Imperatori, B. Mena, M. Mai
!
! Modified: January 2009 (v1.3)
!
! Updated: April 2013 (v1.5)
!
! Updated: September 2014 (v1.5.5.1)
!   Change call set_stf to call set_stf(station).
!
! Updated: May 2015 (v1.5.5.3, gfortran)
!   Add temporary files (velfile and timefile) and PS_flag for raytracing.

use omp_lib

use def_kind; use flags; use geometry; use interfaces, only: spline_interp, write_disk 
use io_file; use scattering; use source_receiver; use stf_data
use waveform; use fault_area

implicit none

! counters
integer(kind=i_single)                        :: station,i,index_delay
! main input-file name 
character(len=256)                            :: input_file
! fixed string variable
character(len=14)                             :: title
! timing barriers, delay for HF
real(kind=r_single)                           :: time_beg,time_end,delay
! temporary array to store time-series
real(kind=r_single),allocatable,dimension(:,:):: in_temp 
! temporary file names
character(len=11)                             :: velfile
character(len=12)                             :: timefile
! Pwave or Swave flag
integer(kind=i_single)                         :: PS_flag
!----------------------------------------------------------------------------------

print*,''
print*,'*** Welcome to the Broad-Band Toolbox (v1.6) ***'
print*,''
print*,'Initialising the code, please type input filename ...'

! reading input file from std-in or shell script
read(*,*) input_file

call cpu_time(time_beg)   !start timing

! read and store input values specified on the input file
call read_inputfile(trim(input_file))         

! prepare input files for raytracing C-code, return some source and station
! parameters for the 1D and the 3D case
call raytracing_input 


!-> start parallel section
! start the raytracer code for P and S waves
!$OMP SECTIONS
!$OMP SECTION           
   PS_flag = 1
   !velfile = 'vel3d_P.bin'
   !timefile = 'time3d_P.out'
   !print*,velfile, timefile
   !call raytracing(hypo,grid,h_step,'vel3d_P.bin','time3d_P.out')
   !call raytracing(hypo,grid,h_step,velfile, timefile)
   call raytracing(hypo,grid,h_step,PS_flag)
   print*,'done Pwave travel time'
!$OMP SECTION        
   PS_flag = 2
   !velfile = 'vel3d_S.bin'
   !timefile = 'time3d_S.out'
   !call raytracing(hypo,grid,h_step,'vel3d_S.bin','time3d_S.out')
   !call raytracing(hypo,grid,h_step,velfile, timefile)
   call raytracing(hypo,grid,h_step,PS_flag)
   print*,'done Swave travel time'
!$OMP END SECTIONS
!<- end parallel section

! read traveltime values for P and S waves at given stations positions
! and calculate source(s)-receiver(s) distance 
call time_distance                             

! delete temporary files from the raytracing code
call cleaner                                  

if (modality_flag /= 0) then
   ! read and/or computes several scattering values  
   call scattering_parameters                    

   ! generate random numbers to be used during coda waves computation
   !call random_sequence                            
endif   

!print*, '---------------------------------------------------'
!print*, 'Reading input low-frequency seismograms...'

! load low-frequency waveforms 
!do station=1,n_stat
!   call read_seis(station,'LF')              
!enddo

! load already computed HF seismograms
if (modality_flag == 0) then
   print*,''
   print*, 'Reading input high-frequency seismograms...'
   do station=1,n_stat
      call read_seis(station,'HF')
   enddo
endif   
            
! define new sampling frequency to avoid aliasing effects
!call sampling_check                           

! set-up or read stf from file
!if (modality_flag /= 0)  call set_stf          

! initialize the LOG file   
!call log_init                                 

! allocating arrays for waveforms
!allocate(bb_seis(npts,3),conv_seis(npts,3))
!if (modality_flag /= 0) allocate(scattgram(npts,3))
!if (modality_flag == 0) allocate(tau(n_stat))

print*, '---------------------------------------------------'
print*, 'Calculating broad-band seismograms...'

if (modality_flag /= 0) then

   !-> start parallel section
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(station,i)
   !loop over stations
   do station=1,n_stat

      print*, '---------------------------------------------------'
      print*, 'Computing broadband for station...'
      print*, station

      ! load low-frequency waveforms 
      call read_seis(station,'LF')              
                              
      ! define new sampling frequency to avoid aliasing effects
      call sampling_check      
      
      ! allocating arrays for waveforms
      if (station == 1) allocate(bb_seis(npts,3),conv_seis(npts,3),scattgram(npts,3))

      ! set-up or read stf from file
      !if (modality_flag /= 0)  call set_stf
      if (modality_flag /= 0)  call set_stf(station)

      ! initialize the LOG file   
      if (station == 1) call log_init  

      ! allocating some arrays before looping
      if (.not.allocated(lf_int)) allocate(lf_int(npts,3))   

      ! calculates average S-wave speed between hypocenter and receiver
      aveVs = sr_hypo(station)/time_s(station)               

      ! generate random numbers to be used during coda waves computation
      call random_sequence 

      ! computes scattered wavelets (i.e. Coda waves) using Zeng's multiple S-to-S theory
      call simcoda(station)                      

      ! convolves scattering seismograms with moment rate function
      call conv_stf_scat(station)                  

      ! interpolates original LF time series to scatt. npts
      do i=1,3	
         call spline_interp(lf_seis(:,i),lf_len,npts,lf_npts,lf_int(:,i))
      enddo 

      ! combines scattering seismograms with LF ones -> broadband 
      call broadband_comp(station)    

      ! write seismograms to disk
      call write_disk(station,'hyb',bb_seis(:,:))

   enddo
   !$OMP END PARALLEL DO
   !<- end parallel section
   
else if (modality_flag == 0) then
  
   !-> start parallel section
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(station,i)
   !loop over stations
   !do station=1,n_stat

      ! allocating some arrays before looping
      !if (.not.allocated(lf_int)) allocate(lf_int(npts,3))   
      
      ! compute lag value for HF seismograms
      !tau(station) = delay(station)             
      !index_delay = nint(abs(tau(station))/(hf_len/(hf_npts-1))) + 1
      
      ! shift HF time-series according to computed delay
      !if (tau(station) > 0) then
         !do i=1,3
            !conv_seis(1:hf_npts,i,station) = hf_seis(:,i,station)
            !hf_seis(1:index_delay-1,i,station) = 0.0
            !hf_seis(index_delay:hf_npts,i,station) = conv_seis(1:hf_npts-index_delay+1,i,station)
         !enddo
      !else if (tau(station) < 0) then
         !do i=1,3
            !conv_seis(1:hf_npts,i,station) = hf_seis(:,i,station)
            !hf_seis(hf_npts-index_delay+2:hf_npts,i,station) = 0.0
            !hf_seis(1:hf_npts-index_delay+1,i,station) = conv_seis(index_delay:hf_npts,i,station)
         !enddo
      !endif   
            
      ! interpolates input LF and HF time series to 2**N npts (here is assumed hf_len == lf_len
      ! and lf_npts == hf_npts)
      !do i=1,3	
         !call spline_interp(lf_seis(:,i,station),lf_len,npts,lf_npts,lf_int(:,i))
         !call spline_interp(hf_seis(:,i,station),hf_len,npts,hf_npts,conv_seis(:,i,station))
      !enddo 
       
      ! combines scattering seismograms with LF ones -> broadband 
      !call broadband_comp(station)              

   !enddo
   !$OMP END PARALLEL DO
   !<- end parallel section

endif

! deallocate memory
deallocate(lf_seis,lf_int)
if (modality_flag == 0) deallocate(hf_seis)
if (verbose_flag == 'off' .and. modality_flag /= 0) deallocate(scattgram,conv_seis)

print*, '---------------------------------------------------'
print*, 'Writing secondary output files on disk...'

do station=1,n_stat
   
   ! write seismograms to disk
   !call write_disk(station,'hyb',bb_seis(:,:,station)) 

   ! verbose stuff
   if (verbose_flag == 'on' .and. modality_flag /= 0) then
   
      if (station == 1) then          !write STF only for one (first) station
         if (.not.allocated(in_temp)) allocate(in_temp(npts_stf,2))
         in_temp(:,1)=t_stf; in_temp(:,2)=stf
         call write_disk(station,'stf',in_temp)
         if (allocated(in_temp)) deallocate(in_temp)
      endif
      call write_disk(station,'ocd',scattgram(:,:))   !coda before convolution
      call write_disk(station,'ccd',conv_seis(:,:))   !coda after convolution
   endif

   ! write other informations on log-file
   call write_log(station)             

enddo


call cpu_time(time_end)    !stop timing

print*, '---------------------------------------------------'
write(5,*)'Computation terminated in: ',(time_end-time_beg),' secs' 

close(5) !Close log-file

! last not-deallocated arrays are automatically deallocated as the program exits
   
END PROGRAM main_broadband
