MODULE def_kind
!
! Definition:
!
!   This is to guarantee portability on different machines. P indicates "precision"
!
! Dependencies:
!
!   All subroutines and modules
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

implicit none

! Real single precision
integer,parameter:: r_single=selected_real_kind(p=6)
! Real double precision
integer,parameter:: r_double=selected_real_kind(p=15)
! Real 'scattering' precision (i.e. only 4 digits after decimal point)
integer,parameter:: r_scat=selected_real_kind(p=4)    
! Integer single representation (i.e. up to 10**9 order)
integer,parameter:: i_single=selected_int_kind(9)  
! Integer double representation (i.e. up to 10**18 order)
integer,parameter:: i_double=selected_int_kind(18)

END MODULE def_kind

!===================================================================================================

MODULE constants
!
! Definition:
!
!   Define several constants used in the code
!
! Dependencies:
!
!   Module def_kind
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

! Pi and 2*Pi
real(kind=r_double),parameter   :: pi=3.141592653589793, pi_double=3.141592653589793*2.
! Pi/2 and degree-to-radiant conversion factor
real(kind=r_double),parameter   :: pi_half=3.141592653589793/2., d2r=3.141592653589793/180.
! Complex triplet
complex(kind=r_single),parameter:: zero=cmplx(0.,0.), one=cmplx(1.0,0.0), zeta=cmplx(0.0,1.0)

END MODULE constants

!===================================================================================================

MODULE earthquake
!
! Description:
!
!   Contains event parameters  
!
! Dependencies:
!
!   Module def_kind
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: March 2019 (v2.1)
!   Add rake for alt computation in source.f90.
!
use def_kind

implicit none

save

! Magnitude of the event (assigned via main input file)
real(kind=r_single):: Mw
! Mechanism of the event (assigned via main input file)
character(len=2)   :: mech
! Rake of the event
real(kind=r_single):: rake


END MODULE earthquake

!===================================================================================================

MODULE flags
!
! Description:
!
!   Various flags 
!
! Dependencies:
!
!   Module def_kind.
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: March 2014 (v1.5.4.1)
!   Add gs_flag (geometric spreading flag).
!
! Updated: June 2015 (v1.6)
!   Add ngaw_flag.
!
! Updated: February 2019 (v2.0)
!   Add merging_flag, merging LFs and HFs in the frequency/time domain.
!   Add infcorr_flag for inter-frequency correlation.
!
use def_kind

implicit none

save

! Flags used in the main input file (velocity model, extended fault, modality flag)
integer(kind=i_single):: vel_flag,ext_flag,modality_flag
! Flag for secondary output (log-file) and LF-HF files format
character(len=4)      :: verbose_flag,lf_kind_flag,hf_kind_flag
! Flag to determine tectonically active/stable region
integer(kind=i_single):: gs_flag
! Flag to determine NGA-west1 or NGA-west2 
integer(kind=i_single):: ngaw_flag
! Flag for merging LFs and HFs in the frequency or time domain
integer(kind=i_single):: merging_flag
! Flag for including inter-frequency correlation, or no
integer(kind=i_single):: infcorr_flag

END MODULE flags

!===================================================================================================

MODULE geometry
!
! Description:
!
!   Store parameters regarding spatial domain and its discretization      
!
! Dependencies:
!
!   Module def_kind
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

save

! Dimensions and grid-spacing of the computational domain (defined in main input file)
real(kind=r_single)                         :: h_step,x_far,y_far,z_far,x_init,y_init,z_init 
! Hypocenter coordinates
real(kind=r_single)                         :: hyp_x,hyp_y,hyp_z
! Number of layers for 1D model (if provided)
integer(kind=i_single)                      :: n_lay    
! Number of grid-points for raytracing code
integer(kind=i_single)                      :: nx,ny,nz
! Arrays for the raytracing code
integer(kind=i_single),dimension(3)         :: grid
real(kind=r_single),dimension(3)            :: hypo
! Array for delay values
real(kind=r_single),allocatable,dimension(:):: tau


END MODULE geometry

!===================================================================================================

MODULE fault_area
!
! Description:
!
!   Store parameters regarding spatial domain and its discretization      
!
! Dependencies:
!
!   Module def_kind
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

save

! Hypocenter coordinates
real(kind=r_single)                         :: Ar

END MODULE fault_area

!===================================================================================================


MODULE io_file
!
! Description:
!
!   Store values related to input/output operations
!
! Dependencies:
!
!   Module def_kind
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: July 2015 (v1.6.1)
!   Change character length from 4 to 10 for lf_ and hf_x, y, and z.
!
! Updated: February 2019 (v2.0)
!   Change character length from 90 to 256.
!   Add k_file to read Kemp_*.bin file name.
!
use def_kind

implicit none

save
 
! Secondary input files (velocity model, stations, extended fault and scattering)
character(len=256):: vel_file,stat_file,ext_file,scat_file
! Optional: 2nd station file and binary input file for HF/LF waveforms
character(len=256):: opt_stat_file,lf_bin_file,hf_bin_file   
! time-series file extensions 
character(len=10) :: lf_x,lf_y,lf_z,hf_x,hf_y,hf_z
! input LF directory and optional (HF) directory for already computed HF files
character(len=256):: lf_in_dir,opt_dir,output_dir
! files for correlation
character(len=256):: k_file

END MODULE io_file

!===================================================================================================

MODULE matching
!
! Description:
!
!   Contains parameters for LF-HF matching operations  
!
! Dependencies:
!
!   Module def_kind
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: March 2014 (v1.5.4)
!   targ_fr and band_wid are removed from SOcomp.par.
!   targ_fr is computed by Mw, and band_wid is fixed.

use def_kind

implicit none

save

! Target frequency and bandwidth (assigned via main input file)
!real(kind=r_single)                           :: targ_fr,band_wid
real(kind=r_single)                           :: targ_fr
real(kind=r_single)                           :: band_wid=0.05
! Actual matching frequency (in the range targ_fr +/- band_wid)
real(kind=r_single),allocatable,dimension(:,:):: match_fr     

END MODULE matching

!===================================================================================================

MODULE scattering
!
! Description:
!
!   Store scattering-related parameters.
!
! Dependencies:
!
!   Module def_kind.
!
! Notes:
!
!   For a complete description of many parameters, see coda.f90
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
!
! Updated: March 2013 (v1.4.2)
!   In scattering.f90, new parameters of fmax, imerg, str_fac, fac, afac and bfac
!   are added to be read from a new form of scattering.dat.
!   In older versions, user-defined fmax was set in this scattering module.
!
! Updated: March 2014 (v1.5.4)
!   fac, ncoda, nfcoda, nscat, merr, abscoeff and trans are fixed parameters,
!   and removed from scattering.dat.
!   time_step is read from scattering.dat.
!
! Updated: December 2016 (v1.6.2)
!   Add aveVp for coda computation.
!   Add loc_aveVs for output in run.log.
!
use def_kind

implicit none

save

! The following parameters MAY change for every receiver
! Average S-wave and P-wave speed
real(kind=r_scat)                         :: aveVs,aveVp
! Site-specific rho, Vs and kappa
real(kind=r_scat),allocatable,dimension(:):: siteR,siteVs,kappa
! Output aveVs in run.log
real(kind=r_scat),allocatable,dimension(:):: loc_aveVs  
! Seed number used to compute scattered wavelets
integer(kind=i_single)                    :: iseed 

!--------add----------------------------------------------------------------
! Flag to determine stle of combining LF and HF, old merging, new with 1 or more subfaults
integer(kind=i_single)                    :: imerg
! Maximum frequency
real(kind=r_scat)                         :: fmax 
! Brune stress parameter (G&P 2010 eqn. 12)
real(kind=r_single)                       :: str_fac
! Scaling parameter for revers fault mechanism
!real(kind=r_single)                       :: fac
! Qk factor 'a' and 'b' (G&P 2010 eqn. 15)
real(kind=r_single)                       :: afac,bfac
! Output decimation factor
integer(kind=i_single)                    :: time_step 
!--------END add-------------------------------------------------------------

! The following parameters DON'T change between receivers and can be changed by the user
! Source rho and Vs 
real(kind=r_scat)                         :: srcR,srcVs   
! Coda envelope and fft points, scattering wavelets
!integer(kind=i_single)                    :: ncoda,nfcoda,nscat
! Seed number used to (optionally) compute some medium scattering properties
integer(kind=i_single)                    :: s_seed
! Highpass corner, transition bandwidth, coda envelope tolerance, absorption coefficient 
!real(kind=r_scat)                         :: hpass,trans,merr,abscoeff               
real(kind=r_scat)                         :: hpass               
! Scattering coefficient, frequency decay, attenuation at f0
real(kind=r_scat)                         :: scatcoeff,fdec,Q                  

! The following parameters are FIXED and should NOT be changed
! Maximum frequency and start-time of scattered wavelets
!real(kind=r_scat),parameter               :: fmax=100.,t0=0.       
real(kind=r_scat),parameter               :: t0=0.       
! Time-domain points of scattered wavelet
integer(kind=i_single)                    :: npts                         
! Scaling parameter for revers fault mechanism
real(kind=r_single)                       :: fac=1.4
! Coda envelope and fft points, scattering wavelets
integer(kind=i_single)                    :: ncoda=500
integer(kind=i_single)                    :: nfcoda=512
integer(kind=i_single)                    :: nscat=1500
! transition bandwidth, coda envelope tolerance, absorption coefficient 
real(kind=r_scat)                         :: merr=0.001               
real(kind=r_scat)                         :: abscoeff=0.01               
real(kind=r_scat)                         :: trans=0.1               

! Random numbers array, necessary to guarantee the same results 
! between serial and parallelized executions
real(kind=i_single),allocatable,dimension(:,:):: random_array

!$OMP THREADPRIVATE(aveVs)

END MODULE scattering

!===================================================================================================

MODULE source_receiver
!
! Description:
!
!   Store different source & stations parameters (X-Y coordinates already shifted, 
!   see io.f90)
!
! Dependencies: 
!
!   Module def_kind
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: September 2014 (v1.5.5.1)
!   Add sr_rrup for CENA Tr_fac calculation in source.f90
!
! Updated: March 2019 (v2.1)
!   Change stat_name and opt_stat_name length from 90 to 128.
!
use def_kind

implicit none

save

! Contains mandatory and optional stations names
character(len=128),allocatable,dimension(:)    :: stat_name, opt_stat_name 
! Number of stations and extended-source cells
integer(kind=i_single)                        :: n_stat,n_cell
! Stations coordinates               
real(kind=r_single),allocatable,dimension(:)  :: xp,yp      
! Extended-source cells coordinates
real(kind=r_single),allocatable,dimension(:)  :: x_cell,y_cell,z_cell  
! S&P-waves traveltimes, hypocentral distances
real(kind=r_single),allocatable,dimension(:)  :: time_s,time_p,sr_hypo   
! Source-receiver distances for each subfault 
real(kind=r_single),allocatable,dimension(:,:):: sr_cell      
! Rrup distances
real(kind=r_single),allocatable,dimension(:)  :: sr_rrup   

END MODULE source_receiver

!===================================================================================================

MODULE stf_data
!
! Description:
!
!   Store Source-Time-Function parameters
!
! Dependencies: 
!
!   Module def_kind
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! Updated: March 2013 (v1.4.2)
!   Change character length for stf_name from len=5 to 25.
!   Add srf_name and Tr_sca.
!
! Updated: February 2019 (v2.0)
!   Change character length from 50 to 256 for srf_name.
!
!
use def_kind

implicit none

save

! Kind of source-time-function selected (from main input file)
character(len=25)                            :: stf_name
! srf file name
character(len=256)                            :: srf_name
! Number of time-points for the source-time-function 
integer(kind=i_single)                      :: npts_stf
! Source-time-function (amplitude values) and its time vector
real(kind=r_single),allocatable,dimension(:):: stf ,t_stf
! Tapering window
real(kind=r_single),allocatable,dimension(:):: han_win
! Source-time-function duration and its starting point (in sec.) 
real(kind=r_single)                         :: total,ton                          
! scaling factor for a new stf to get t1, t1 = Tr_sca * Tr
real(kind=r_single)                         :: Tr_sca                          
                                            
END MODULE stf_data

!===================================================================================================

MODULE waveform
!
! Description:
!
!   Store LF and HF time-series and some related parameters
!
! Dependencies: 
!
!   Module def_kind   
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
!
! Updated: February 2019 (v2.0)
!   Add lf_dt and v_npts.
!
use def_kind

implicit none

save

! LF and (optional) HF waveforms
real(kind=r_single),allocatable,dimension(:,:):: lf_seis,hf_seis
! Scatterograms (Coda waves) 
real(kind=r_single),allocatable,dimension(:,:):: scattgram
! Scatterograms convolved with source-time-function
real(kind=r_single),allocatable,dimension(:,:):: conv_seis 
! Final broad-band seismograms 
real(kind=r_single),allocatable,dimension(:,:):: bb_seis   

! Number of points for LF seismograms 
integer(kind=i_single)                          :: lf_npts                                         
! LF time-series length 
real(kind=r_single)                             :: lf_len 
! LF waveforms interpolated (to avoid spectral aliasing up to fmax)
real(kind=r_single),allocatable,dimension(:,:)  :: lf_int          

! Number of points for already computed HF seismograms 
integer(kind=i_single)                          :: hf_npts                                         
! HF time-series length 
real(kind=r_single)                             :: hf_len 
    
!$OMP THREADPRIVATE(lf_int)

! dt for LFs 
real(kind=r_single)                             :: lf_dt 

! virtual npts for BBsynthetics to keep the same dt
! v_npts = 2^15-1 * n + 1, where n = ceiling(lf_len/tmp_lf_len)
integer(kind=i_single)                          :: v_npts                                         


END MODULE waveform

!===================================================================================================

MODULE vel_model
!
! Description:
!
!   Store parameters in 1d velocity model and srf
!   Originaly stored in io.f90
!   Use the velocity model in source.f90 to calculate mu and Mw
!
! Dependencies: 
!
!   Module def_kind   
!
! Author: R. Takedatsu 
!
! Modified: March 2013 (v1.4.2)
!
! Updated: June 2015 (v1.6)
!   Add alt used in source.f90 and composition.f90
!
! Updated: December 2016 (v1.6.2)
!   Add tinit obtained from source.f90
!
! Updated: February 2019 (v2.0)
!   Add dtop obtained from source.f90

use def_kind

implicit none

save


! number of subfaults in srf file
integer(kind=i_single)                        :: nsub 
! Mo calculated from srf file in [Nm] 
real(kind=r_single)                           :: total_Mo
! fault dip and dtop from srf file
real(kind=r_single)                           :: dip,dtop
! total area of subfaults from srf [cm2]
real(kind=r_single)                           :: sum_area
! ave_vs * ave_d
real(kind=r_single)                           :: vsd_ave,d_avef,vs_avef
! facter alpha_T, scaling the rise time as a function of dip, eqn (9) in G&P(2010)
real(kind=r_single)                           :: alt
! minimum initiation time in srf file
real(kind=r_single)                           :: tinit
! x,y,and depth of each subfault [km] from srf
real(kind=r_single),allocatable,dimension(:)  :: fx,fy,fz
! Vs [cm/s] & density [g/cm3] from velocity model file, for each subfault
real(kind=r_single),allocatable,dimension(:)  :: vscu,rhcu
! each subfault area [cm2] from srf file
real(kind=r_single),allocatable,dimension(:)  :: areas
! depth,vp,vs,density,Qp,Qs from velocity model file. size = n_lay (declared in MODULE geometry)
real(kind=r_single),allocatable,dimension(:)  :: depth,vp,vs,rh,Qp,Qs,thk
! moment ratio; slip(i)*mu(i)*total fault area / total fault moment
real(kind=r_single),allocatable,dimension(:)  :: mo_ratio

END MODULE vel_model

!===================================================================================================

MODULE tmp_para
!
! Description:
!
!   These parameters are used to avoid Segumentation fault
!   in interpolation.f90, at reshape command.
!   Most PlanA and all PlanB events have lf-length=102.375 sec.
!   Lf-lengths for all events are stored based on the lf-length=102.375 sec.
!
! Dependencies: 
!
!   Module def_kind   
!
! Author: R. Takedatsu 
!
! Modified: April 2014 (v1.5.5)
!
!

use def_kind

implicit none

save


! fixed npts for scattergram, 2^15
integer(kind=i_single),parameter              :: tmp_npts=32768
! fixed lf-length for scattergram
real(kind=r_single),parameter                 :: tmp_lf_len=102.3750

END MODULE tmp_para

!===================================================================================================

MODULE read_Kemp
!
! Description:
!
!   These parameters are used to read Kemp.txt, Cholesky factor of frequency
!   correlation matrix, lower triangular matrix.
!
! Dependencies: 
!
!   Module def_kind   
!
! Author: N. Wang 
!
! Modified: February 2019 (v2.0)
!
!

use def_kind

implicit none

save

! length (number of elements in each dimension) of matrix K (lower triangular)
integer(kind=i_single)               :: nk
! starting point of corresponding frequency of Kemp
integer(kind=i_single)               :: nks
! Cholesky factor of frequency correlation matrix, lower triangular matrix
real(kind=r_single),allocatable,dimension(:,:) :: Kemp

END MODULE read_Kemp

!===================================================================================================
