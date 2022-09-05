SUBROUTINE scattering_parameters
!-------------------------------------------------------------------------------
!
! Description:
!
!   Reads user-specified input file for scattering parameters and stores them.
!   Here, default values may be completely or partially overridden
!
! Dependencies:
! 
!   Function rand_pdf. All values passed via modules.
!
! Notes: 
!
!   Absorption coefficient, kappa and Q values follow uniform distribution,
!   scattering coefficient is drawn from a normal distribution instead.
!   For a complete description of scattering parameters, see coda.f90 file.
!   ISEED is a seed used for generating scattering wavelets, S_SEED is used
!   for medium-related scattering properties.
!
! Authors: W. Imperatori, B. Mena
!
! Modified: January 2009 (v1.3)
!
! Updated: March 2013 (v1.4.2)
!   Add new parameters in scattering.dat, fmax, imerg, str_fac, fac, afac, bfac
!   and Tr_sca.
!
! Updated: March 2014 (v1.5.4)
!   Add time_step in scattering.dat.
!   Remove fac, ncoda, nfcoda, nscat, merr, abscoeff and trans from scattering.dat,
!   and become fixed parameters in module_bbtoolbox.f90.
!
! Updated: March 2014 (v1.5.4.1)
!   Add gs_flag in scattering.dat.
!
! Updated: February 2015 (v1.5.5.3)
!    Change str_fac
!    if str_fac = 0.0, if gs_flag = 1, str_fac = 50e6
!                      if gs_flag = 2, str_fac = (225+16*Ztor)e6 in srf_read
!    if str_fac = real, user-defined
!
! Updated: May 2015 (v1.5.5.4)
!    Add gs_flag = 3 for Japanese events
!
! Updated: June 2015 (v1.6)
!    Add ngaw_flag, if ngaw_flag = 1, NGA-west1
!                   if ngaw_flag = 2, NGA-west2
!
! Updated: February 2019 (v2.0)
!   Add merging_flag, for LFs and HFs
!                     if merging_flag = 1, in the frequency domain
!                     if merging_flag = 2, in the time domain
!   Add infcorr_flag, for inter-frequency correlation flag
!                     if infcorr_flag = 0, no correlation
!                     if infcorr_flag = 1, apply correlation
!
! Update: March 2019 (v2.1)
! Author: N. Wang
!   Add cseed, for inter-frequency correlation random seed
!              if cseed = -1, seed defined by system-time
!              else, seed defined by user
!
! Update: April 2020 (v2.1)
! Author: N. Wang
!   Add corr_flag, for correlation flag
!                     if corr_flag = 0, no correlation
!                     if corr_flag = 1, apply only inter-frequency correlation
!                     if corr_flag = 2, apply spatial correlation

use def_kind; use interfaces, only: rand_pdf 
use io_file, only: scat_file; use scattering 
use stf_data, only: Tr_sca; use flags 

implicit none

! scattering coefficient mean value
real(kind=r_scat)                              :: scat_mean
! min, max and std. dev. of scattering coefficient (normal distribution)
real(kind=r_scat),parameter                    :: scat_min=0.005, scat_max=0.050, scat_std=0.005
! min and max of absorption coefficient (uniform distribution)
real(kind=r_scat),parameter                    :: abs_min=0.005, abs_max=0.030
! min and max of Q0 value (uniform distribution)
real(kind=r_scat),parameter                    :: Q_min=100.0, Q_max=500.0
! min and max of kappa value (uniform distribution)
real(kind=r_scat),parameter                    :: kappa_min=0.010, kappa_max=0.1
! dummy variables and indexes
!integer(kind=i_single)                         :: tmp_iseed,tmp_sseed,tmp_ncoda,tmp_imerg
!integer(kind=i_single)                         :: tmp_nfcoda,tmp_nscat,seed_size,i
!real(kind=r_scat)                              :: tmp_Q,tmp_abscoeff,tmp_scatcoeff,tmp_kappa
!real(kind=r_scat)                              :: tmp_hpass,tmp_trans,tmp_merr,tmp_fdec
!real(kind=r_scat)                              :: tmp_fmax,tmp_str_fac,tmp_fac
integer(kind=i_single)                         :: tmp_iseed,tmp_sseed,tmp_imerg,tmp_gsf
integer(kind=i_single)                         :: seed_size,i,tmp_time_step,tmp_ngawf
integer(kind=i_single)                         :: tmp_merging,tmp_corr
real(kind=r_scat)                              :: tmp_Q,tmp_scatcoeff,tmp_kappa
real(kind=r_scat)                              :: tmp_hpass,tmp_fdec
real(kind=r_scat)                              :: tmp_fmax,tmp_str_fac
real(kind=r_scat)                              :: tmp_afac,tmp_bfac,tmp_Tr_sca
! local array necessary if random seed has been selected
integer(kind=i_single),allocatable,dimension(:):: array_seed

!-------------------------------------------------------------------------------


scat_mean=(scat_min+scat_max)/2.    !compute mean scattering coeff. value

! open file with user-specified scattering parameters
open(1,file=scat_file,status='old')
read(1,*);read(1,*);read(1,*)

! here below, seed number for scattering wavelets (see coda.f90).
read(1,*) tmp_iseed   
if (tmp_iseed == 0)  iseed=8888      !default value
if (tmp_iseed == -1) then            !random seed based on system time
   call random_seed(size=seed_size)
   if (.not.allocated(array_seed)) allocate(array_seed(seed_size))
   call random_seed() 
   call random_seed(get=array_seed)
   iseed=array_seed(1)               !to be sure, takes always first element
else if( (tmp_iseed /= 0).and.(tmp_iseed /= -1) ) then
   iseed=tmp_iseed                   !user-defined seed
endif

! here below, seed number for medium scattering parameters
read(1,*) tmp_sseed             
if(tmp_sseed == -1) then            !random seed based on system time
   if (tmp_iseed == -1) then
      s_seed=iseed/37               !this is to avoid very similar values with previous seed (iseed)
   else
      call random_seed(size=seed_size)
      if(.not.allocated(array_seed)) allocate(array_seed(seed_size))
      call random_seed() 
      call random_seed(get=array_seed)
      s_seed=array_seed(1)          !to be sure, takes always first element
   end if
else if( (tmp_sseed /= 0).and.(tmp_sseed /= -1) ) then
   s_seed=tmp_sseed
else
   s_seed=461                       !default value
endif

        
! deallocate and re-allocate local seed array. This array guarantees that the same random number is
! drawn for the same scattering parameter even if the user changes the number of parameters to be 
! computed. Size is 4 because potentially one has to calculate Q, scatcoeff, kappa and abscoeff
if(allocated(array_seed)) deallocate(array_seed)
if(.not.allocated(array_seed)) allocate(array_seed(4))

array_seed(1)=s_seed   !first element is specified seed

! other seeds are higher
do i=2,4
   array_seed(i)=s_seed+array_seed(i-1)
enddo

! Q (i.e. Q0) section
read(1,*) tmp_Q
if (tmp_Q == 0.0) Q=150.0                                        !default value
if (tmp_Q == 1.0) Q=rand_pdf(array_seed(1),Q_min,Q_max,'unif')   !random value   
if ( (tmp_Q /= 0.0).and.(tmp_Q /= 1.0) ) Q=tmp_Q                 !user-defined value

! scattering coefficient section
read(1,*) tmp_scatcoeff
if (tmp_scatcoeff == 0.0) scatcoeff=0.03                                        !default value
if (tmp_scatcoeff == 1.0) then 
   scatcoeff=rand_pdf(array_seed(2),scat_mean,scat_std,'norm')                  !random value
   if (scatcoeff < scat_min) scatcoeff=scat_min
   if (scatcoeff > scat_max) scatcoeff=scat_max
endif
if( (tmp_scatcoeff /= 0.0).and.(tmp_scatcoeff /= 1.0) ) scatcoeff=tmp_scatcoeff !user-defined value

! kappa section. NOTE: assign values only where no values have been previously specified 
! through the stations file 
read(1,*) tmp_kappa

if (tmp_kappa == 0.0) then 
   where (kappa <= 0.0)
      kappa = 0.03                                              !default value 
   end where 
else if(tmp_kappa == 1.0) then
   where (kappa <= 0.0)
      kappa=rand_pdf(array_seed(3),kappa_min,kappa_max,'unif')  !random value
   end where      
else if ( (tmp_kappa /= 0.0) .and. (tmp_kappa /= 1.0) ) then
   where (kappa <= 0.0)
      kappa=tmp_kappa                                           !user-defined value
   end where 
endif    
   
! fdec section   
read(1,*) tmp_fdec
if (tmp_fdec == 0.0) then 
   fdec=0.8          !default value
else
   fdec=tmp_fdec     !user-defined value
endif   

! hpass section
read(1,*) tmp_hpass
if (tmp_hpass == 0.0) then
   hpass=0.2         !default value
else   
   hpass=tmp_hpass   !user-defined value
endif   

! jump to the second part of the scattering file
read(1,*);read(1,*);read(1,*)

!------adding more parameters ----------------------------------------------

! fmax section, previously set in module_bbtoolbox.f90
read(1,*) tmp_fmax
if (tmp_fmax == 0.0) then
   fmax=100.0          ! default value
else
   fmax=tmp_fmax       ! user-defined value
endif

! imerg section
read(1,*) tmp_imerg
if (tmp_imerg > 2) then
   imerg=2             ! default value
   print*,'imerg needs to be 0, 1 or 2'
else
   imerg=tmp_imerg     ! user-defined value
endif

! str_fac section
!read(1,*) tmp_str_fac
!if  (tmp_str_fac == 0.0) then
!   str_fac=50.e6       ! default value
!else
!   str_fac=tmp_str_fac ! user-definded value
!endif

! fac section
!read(1,*) tmp_fac
!if (tmp_fac == 0.0) then
!   fac=1.0             ! defaut value
!else
!   fac=tmp_fac         ! user-defined value
!endif

! afac section
read(1,*) tmp_afac
if (tmp_afac == 0.0) then
   afac=41.0           ! default value
else
   afac=tmp_afac       ! user-defined value
endif

! bfac section
read(1,*) tmp_bfac
if (tmp_bfac == 0.0) then
   bfac=34.0           ! defaul value
else
   bfac=tmp_bfac       ! user-defined value
endif

! Tr_sca section
read(1,*) tmp_Tr_sca
if (tmp_Tr_sca == 0.0) then
   Tr_sca=0.1           ! defaul value
else
   Tr_sca=tmp_Tr_sca    ! user-defined value
endif

! time_step section
read(1,*) tmp_time_step
if (tmp_time_step < 1) then
   print*,'change time_step'
else
   time_step=tmp_time_step    ! user-defined value
endif

! gs_flag section
read(1,*) tmp_gsf
!if (tmp_gsf == 1 .or. tmp_gsf == 2) then
if (tmp_gsf >= 1 .and. tmp_gsf <= 3) then ! add gs_flag=3=Japanese events
   gs_flag=tmp_gsf       ! 1=western CA, 2=eastern North America
else
   print*,'change gs_flag'
endif 

print*,'gs_flag= ',gs_flag

! str_fac section
read(1,*) tmp_str_fac
if  (tmp_str_fac == 0.0) then
   !if (gs_flag == 1) then
   if (gs_flag == 1 .or. gs_flag == 3) then ! add for Japanese events
       str_fac=50.e6
   elseif (gs_flag == 2) then
      str_fac=0.0
   endif
else
   str_fac=tmp_str_fac ! user-definded value
endif

! ngaw_flag section
read(1,*) tmp_ngawf
if (tmp_ngawf == 1 .or. tmp_ngawf == 2) then
   ngaw_flag=tmp_ngawf       ! 1=NGA-west1, 2=NGA-west2
else
   print*,'change ngaw_flag'
endif 

! merging_flag section
read(1,*) tmp_merging
if (tmp_merging == 1 .or. tmp_merging == 2) then
   merging_flag=tmp_merging
else
   print*,'change merging_flag'
endif

! corr_flag section
read(1,*) tmp_corr
if (tmp_corr == 0 .or. tmp_corr == 1 .or. tmp_corr == 2) then
   corr_flag=tmp_corr
else
   print*,'change corr_flag'
endif

! seed number for correlation section
read(1,*) cseed




! jump to the third part of the scattering file
!read(1,*);read(1,*);read(1,*)

!------END adding more parameters ------------------------------------------

! ncoda section
!read(1,*) tmp_ncoda
!if (tmp_ncoda == 0) then
!   ncoda=500         !default value
!else   
!   ncoda=tmp_ncoda   !user-defined value
!endif   

! nfcoda section
!read(1,*) tmp_nfcoda
!if (tmp_nfcoda == 0) then
!   nfcoda=512        !default value
!else
!   nfcoda=tmp_nfcoda !user-defined value
!endif

! nscat section
!read(1,*) tmp_nscat
!if (tmp_nscat == 0) then
!   nscat=1500         !default value
!else   
!   nscat=tmp_nscat    !user-defined value
!endif

! merr section
!read(1,*) tmp_merr
!if (tmp_merr == 0.0) then 
!   merr=0.001          !default value
!else   
!   merr=tmp_merr       !user-defined value
!endif

! abscoeff section
!read(1,*) tmp_abscoeff
!if (tmp_abscoeff == 0.0) abscoeff=0.03                                           !default value
!if (tmp_abscoeff == 1.0) abscoeff=rand_pdf(array_seed(4),abs_min,abs_max,'unif') !random value
!if ( (tmp_abscoeff /= 0.0).and.(tmp_abscoeff /= 1.0) ) abscoeff=tmp_abscoeff     !user-defined value

! trans section
!read(1,*) tmp_trans
!if(tmp_trans == 0.0) then
!   trans=0.1          !default value
!else   
!   trans=tmp_trans    !user-defined value
!endif

! deallocate memory
if(allocated(array_seed)) deallocate(array_seed)

END SUBROUTINE scattering_parameters

!===================================================================================================
