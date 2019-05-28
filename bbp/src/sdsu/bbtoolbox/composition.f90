!contains changes that scales the HF amplitude to the expected value at the merging frequency
!as provided from Graves and Pitarka (2010) and earlier refs. No frequency dependent impedance
!factor included. Arbitrary factor of 1.7 included to match NR, LOMAP, Whittier

SUBROUTINE broadband_comp(station)
!------------------------------------------------------------------------------------
!
! Description:
!
!   External subroutine, calculates the composite broad-band seismograms, given the 
!   low and high frequency time-series. The waveform data are assumed to be velocity
!   seismograms [m/sec] for the three components of motion x, y, and z. It operates
!   in the frequency domain using the approach of Mai and Beroza.
!
! Dependencies:
!
!   Subroutine bb_calc  
!
! References:
!
!   Numerical Recipes
!
! Authors: W. Imperatori, M. Mai, B. Mena
!
! Modified: January 2009 (v1.3)
!
! Updated: March 2013 (v1.4.2)
!   Change taper_len from 10 points to the length between P ans S arrivals.
!
! Updated: May 2013 (v1.5.1)
!   Back to the tapering to the old version.
!   Fix the error,  ratio calculation for imerg=0. 
!
! Updated: October 2013 (v1.5.2)
!   Add output decimation factor, time_step.
!
! Updated: March 2014 (v1.5.4)
!   time_step is read from scattering.dat.
!
! Updated: July 2015 (v1.6.1)
!   Add tmp_para and change dt computation at the first step.
!
! Updated: December 2016 (v1.6.2)
!   Add minimum initiation time on time_p for zero padding.
!
! Updated: February 2019 (v2.0)
!   Add time-domain merging system by R. Graves (Dec. 5 2016), "Theoretical Constraints on
!       the Amplitude Spectra of Matched Filters".
!       The acc_ratio for the time-domain merging is from the ratio in bb_calc, computed from acc_spec.
!   Cahnge taper_len around P-wave arrival time from 10 to 30.
!   Add call infcorr routine.
!
use constants; use def_kind; use flags;
use scattering, only: npts,time_step
use source_receiver; use waveform; use fault_area;  use tmp_para
use vel_model, only: tinit ! add v162
use matching, only: targ_fr

implicit none

! actual station number
integer(kind=i_single),intent(in)              :: station
! indexes, taper length
integer(kind=i_single)                         :: i,ind,taper_len
! for output decimation factor
integer(kind=i_single)                         :: bb_npts,k
! time-step
real(kind=r_single)                            :: dt
! frequency for acc_spec
real(kind=r_single)                            :: acc_spec_freq
! get ratio from bb_calc
real(kind=r_single),allocatable,dimension(:)   :: acc_ratio
! taper
real(kind=r_single),allocatable,dimension(:,:) :: taper
! for tapering and bb_seis
real(kind=r_single),allocatable,dimension(:,:) :: tmp_seis
!------------------------------------------------------------------------------------

! delta-t of both high and low-frequency (already interpolated) time-series   
dt = lf_len/(v_npts-1)

! calculate HFs scaling factor from bb_calc for each component
if (.not.allocated(acc_ratio)) allocate(acc_ratio(3))

acc_spec_freq = 50.0

do i=1,3 
   call bb_calc(i,dt,station,acc_ratio(i),acc_spec_freq)
enddo

!!! ----- time domain merging only ----------

if (merging_flag == 2) then

   print*,'*** start time domain merging ***'

   bb_seis(1:npts,1:3)=0.0

   do i=1,3
      do k=1,npts
         bb_seis(k,i)=lf_int(k,i) + (conv_seis(k,i)*acc_ratio(i))
      enddo
   enddo
endif

!!! ----- end time domain merging --------

if (infcorr_flag == 1) then
   ! calculate corrlated bb seismograms for each component
   call infcorr(npts,dt,station)
endif

if (modality_flag /= 0) then   
   ! find time-index corresponding to P-wave arrival 
   !ind=nint(time_p(station)/dt)
   ind=nint((time_p(station)+tinit)/dt) ! change v162

   ! setting up the taper (rising half of Hanning window) to avoid spikes in acceleration 
   ! at P-wave onset
   ! taper length in points (this is the effective length)
   !taper_len=10
   taper_len=30

   ! reduce taper if too long
   if (floor(taper_len/2.) >= ind) then
      taper_len = 2 * ind - 2
      call warning_handling(1,'null','BROADBAND_COMP (composition.f90)') 
   endif

   ! allocate local array for taper
   if (.not.allocated(taper)) then
       allocate(taper(taper_len,3))
   endif

   ! equation of Hanning window: h(i)=0.5*(1-cos(2*pi*i/N)) ---> see Numerical Recipes, page: 547
   ! computing only the first (rising) half
   do i=1,taper_len
      taper(i,:) = 0.5 * (1 - cos(2 * pi * i / (2*taper_len)) )
   enddo

   ! set to zero all amplitudes values about before the P-wave arrival and apply taper
   !bb_seis(1:ind-floor(taper_len/2.),:)=0.0 
   !bb_seis(ind-floor(taper_len/2.)+1:ind+ceiling(taper_len/2.),:) =     &
   !bb_seis(ind-floor(taper_len/2.)+1:ind+ceiling(taper_len/2.),:)*taper

   if (.not.allocated(tmp_seis)) allocate(tmp_seis(npts,3))
   tmp_seis=bb_seis
   tmp_seis(1:ind-floor(taper_len/2.),:)=0.0
   tmp_seis(ind-floor(taper_len/2.)+1:ind+ceiling(taper_len/2.),:) =     &
   bb_seis(ind-floor(taper_len/2.)+1:ind+ceiling(taper_len/2.),:)*taper

   bb_seis(:,:)=0

   do k=1,3
      bb_npts=0
      do i=1,npts,time_step
         bb_npts=bb_npts+1
         bb_seis(bb_npts,k)=tmp_seis(i,k)
      enddo
   enddo

   if (allocated(tmp_seis)) deallocate(tmp_seis)

   ! deallocate memory
   if (allocated(taper)) then
      deallocate(taper)   
   endif
endif

END SUBROUTINE broadband_comp

!===================================================================================================

!SUBROUTINE bb_calc(comp,dt,station)
SUBROUTINE bb_calc(comp,dt,station,ratio,f_acc_spec)
!-----------------------------------------------------------------------------------
!
! Description:
!
!   External subroutine, merges the HF and LF seismograms into broadband time-series
!   using the approach of Mai and Beroza
!
! Dependencies:
!
!   Subroutine fast_fourier
!
! Notes:
!
!   Input time-series MUST have the same number of points and time-step. The FFT/IFFT
!   subruotine implemented works only with 2**n time-series. 
!   Frequency vector goes from 0 to f_nyq and back to 0+df (npts points)
!
! References:
!
!
!
! Authors: W. Imperatori, M. Mai, B. Mena
!
! Modified: January 2009 (v1.3)
!
! Updated: December 2012 (v1.4.1)
!   Change calculation of averages of hf and lf, and scaling hf by ratio.
!
! Updated: March 2013 (v1.4.2)
!   Change calculation of match_fr, fscale, averages of hf and lf, and ratio.
!   Add function, acc_spec.
!   Add imerg flag.
!
! Updated: April 2014 (v1.5.5)
!    tmp_npts and tmp_lf_len are in the module tmp_para.
!
! Updated: july 2015 (v1.6.1)
!    Do not need tmp_para, tmp_dt for df computation.
!
! Updated: February 2019 (v2.0)
!   Add ratio and f_acc_spec in subroutine, such that
!   bb_calc(comp,dt,station,ratio,f_acc_spec).
!   ratio: for HFs scaling for both frequency-domain merging and time-domain merging.
!   f_acc_spec: specify frequency for acc_spec.
!   Separate ratio computation imerg=0 or imerg=1 or 2.
!   Use v_npts for df computation.
! 
use interfaces, only: four1d
use constants; use def_kind; use matching
use scattering; use waveform; use earthquake
use fault_area 
use flags, only: merging_flag
!use tmp_para

implicit none

! current component and station number
integer(kind=i_single),intent(in)           :: comp,station
! FOR ACC_SPEC TEST
integer(kind=i_single)                      :: sca_index, fscale_11, fscale_22, acc_size1 
! time-step
real(kind=r_single),intent(in)              :: dt,f_acc_spec
! indexes, counters, dummies 
integer(kind=i_single)                      :: index_diff,index_min,fscale_1,fscale_2,acc_size
! flag and counters 
integer(kind=i_single)                      :: trasf_sign,i
! index of real target frequency and Nyquist frequency 
integer(kind=i_single)                      :: targ_index,f_nyq_index
! window for search: half extension, lower and upper bounds, total extension 
integer(kind=i_single)                      :: win_fr,win_low,win_up,win_npts
! array for complex time-series
complex(kind=r_single),dimension(npts)      :: hf_four,lf_four,broad_complex
! arrays for high and low frequency spectra (amplitude and phase) 
real(kind=r_single),dimension(npts)         :: Am_hf,Ph_hf,Am_lf,Ph_lf
! arrays for scaling 
real(kind=r_single),dimension(npts/2 +1)    :: Acc_hf,Acc_lf
! frequency vectors, windowing functions, arrays for phase and amplitude of merged time-series
real(kind=r_single),dimension(npts)         :: f,fn_lf,fn_hf,phase,amp
! arrays for windowed spectra (amplitude and phase) 
real(kind=r_single),dimension(npts)         :: hf_am_win,lf_am_win,hf_ph_win,lf_ph_win
! array for amplitude differences 
real(kind=r_single),allocatable,dimension(:):: diff_amp
! actual target frequency, averaged spectra and their ratio, Nyquist frequency
!real(kind=r_single)                         :: f_targ,Av_hf,Av_lf,ratio,f_nyq,Av_hf1
real(kind=r_single)                         :: f_targ,Av_hf,Av_lf,f_nyq,Av_hf1
! other averaged spectra
real(kind=r_single)                         :: Av_hf11,Av_lf11
! delta-f
real(kind=r_single)                         :: df,A,acc_spec
! ratio
real(kind=r_single),intent(out)              :: ratio

!----------------------------------------------------------------------------------- 

! computing frequency-delta
df = 1 / (v_npts * dt)
!df = 1 / (npts * dt)    !according to NR

! compute nyquist frequency (its index is npts/2 +1) 
! (since npts is multiple of 2, f_nyq = {(npts+2)/2 -1)*df}
f_nyq = 0.5 * (1 / dt)         
f_nyq_index = npts/2 + 1

! create frequency vector
f(1)=0.
do i=1,f_nyq_index-1
   f(i+1)=f(i)+df   
enddo
do i=f_nyq_index,npts-1
   f(i+1)=f(i)-df
enddo   
   
! flag for FFT (i.e. from time to frequency)
trasf_sign=1    

hf_four=cmplx(conv_seis(:,comp),0.0)   !prepare complex array for HF

! calling function for FFT
call four1d(hf_four,trasf_sign)

! extracts Amplitude and Phase spectra (HF)
Am_hf=cabs(hf_four)
Ph_hf=atan2(aimag(hf_four),real(hf_four))

lf_four=cmplx(lf_int(:,comp),0.0)   !prepare complex array for LF	

! calling function for FFT
call four1d(lf_four,trasf_sign)

! extracts Amplitude and Phase spectra (LF)
Am_lf=cabs(lf_four)
Ph_lf=atan2(aimag(lf_four),real(lf_four))

! index of target frequency in frequency vector 
do i=1,f_nyq_index
   if (f(i) >= targ_fr-df .and. f(i) <= targ_fr+df) then
      targ_index=i
      exit
   endif   
enddo

! index of target frequency for scaling in frequency vector 
do i=1,f_nyq_index
   !if (f(i) >= 50.-df .and. f(i) <= 50.+df) then
   if (f(i) >= f_acc_spec-df .and. f(i) <= f_acc_spec+df) then
      sca_index=i
      exit
   endif
enddo

! actual target frequency (user target-fr. "projected" onto frequency vector)
f_targ=f(targ_index)   

! search window half-length (number of points in frequency)	
win_fr = nint(band_wid/df)

! define search region on frequency vector (in points)
win_low = targ_index - win_fr
win_up = targ_index + win_fr
win_npts = (win_up-win_low) + 1      !total window-length for search (i.e. search area)

! check for search window 
if ((win_low < 1) .or. (win_up > f_nyq_index)) then 
   call error_handling(6,'0','BB_CALC (composition.f90)')
endif

! allocating vector where to store amplitude differences
if (.not.allocated(diff_amp)) allocate(diff_amp(win_npts))

! computes spectral amplitude differences inside the search area
diff_amp(1:win_npts)=abs(Am_hf(win_low:win_up)-Am_lf(win_low:win_up))

! find index of minimum inside amplitude-difference vector	
index_diff=minloc(diff_amp,dim=1) 

! deallocate memory 
if (allocated(diff_amp)) deallocate(diff_amp)

! find index of minimum inside positive frequency vector  
index_min=index_diff-1+win_low    

! determine matching frequency, which corresponds to the minimum amplitude difference
match_fr(comp,station)=f(index_min)

! ----> section for AMPLITUDE SCALING (to avoid spectral jumps at the matching frequency) <----

! scaling is based on computing the ratio of LF/HF acceleration amplitude spectra 
! averaged over a wide window

! derive acceleration amplitude spectra from velocity amplitude spectra 
Acc_hf=Am_hf(1:f_nyq_index)*2.*pi*f(1:f_nyq_index) 
Acc_lf=Am_lf(1:f_nyq_index)*2.*pi*f(1:f_nyq_index)
! note: here above the problem WAS that f went from -f_nyq to f_nyq, while Am from 0 to 2*fnyq-1 

if (imerg .eq. 0) then
   ! compute spectral averages over a large window (3 times search window)
   fscale_1=index_min-3*win_fr; if (fscale_1 <= 0) fscale_1=1
   fscale_2=index_min+3*win_fr; if (fscale_2 > f_nyq_index) fscale_2=f_nyq_index
   ! fscale_1=targ_index-3*win_fr; if (fscale_1 <= 0) fscale_1=1
   ! fscale_2=targ_index+3*win_fr; if (fscale_2 > f_nyq_index) fscale_2=f_nyq_index
   acc_size=(fscale_2-fscale_1)+1

   ! check for the above secondary search window
   ! if ( (fscale_1 == 1) .or. (fscale_2 == f_nyq_index) ) then
   !    call warning_handling(2,'0','BB_CALC (composition.f90)')
   ! endif

   ! compute averages
   ! Av_hf = sum(Acc_hf(fscale_1:fscale_2))/(acc_size-1)   
   ! Av_lf = sum(Acc_lf(fscale_1:fscale_2))/(acc_size-1)
   Av_hf = sum(Acc_hf(fscale_1:fscale_2))/(acc_size)   
   Av_lf = sum(Acc_lf(fscale_1:fscale_2))/(acc_size) 
   ratio = Av_lf/Av_hf   !averages ratio
   print*,'Av_lf, Av_hf, ratio for imerg=0: ',Av_lf,Av_hf,ratio

else
   ! compute spectral averages over a large window (3 times search window)
   fscale_11=sca_index-3*win_fr; if (fscale_11 <= 0) fscale_11=1
   fscale_22=sca_index+3*win_fr; if (fscale_22 > f_nyq_index) fscale_22=f_nyq_index
   acc_size1=(fscale_22-fscale_11)+1

   ! compute averages
   Av_hf11 = sum(Acc_hf(fscale_11:fscale_22))/(acc_size1)   
   !Av_lf11 = sum(Acc_lf(fscale_11:fscale_22))/(acc_size1)   

   ! new merging, imerg=1: one big subfault, =2: more subfaults
   Av_hf1=acc_spec(f_acc_spec,station) ! use 50Hz for fmax=100Hz 
   ratio = Av_hf1/Av_hf11
   print*,'imerg,Av_hf1,Av_hf11,ratio for imerg>0= ', &
           imerg,Av_hf1,Av_hf11,ratio
endif

!!! ----- frequency domain merging -----

if (merging_flag == 1) then

   print*,'*** start frequency domain merging ***'

   Am_hf = Am_hf*ratio

! ----------------------> end section for amplitude scaling <---------------------- 

   ! windowing function for low frequency
   where ( f <= match_fr(comp,station) )
      fn_lf=1.0
   elsewhere
      fn_lf=0.
   end where

   ! windowing function for high frequency 
   where ( f > match_fr(comp,station) )
      fn_hf=1.0
   elsewhere
      fn_hf=0.
   end where

   ! multiply phase and amplitude spectra with the windowing function
   hf_ph_win=Ph_hf*fn_hf; lf_ph_win=Ph_lf*fn_lf
   hf_am_win=Am_hf*fn_hf; lf_am_win=Am_lf*fn_lf      

   ! compose low and high frequency phase spectra	
   phase=hf_ph_win+lf_ph_win

   ! compose low and high frequency amplitude spectra
   amp=hf_am_win+lf_am_win

   ! combine amplitude and phase spectra (whole broadband spectrum) 
   broad_complex=amp*exp(zeta*phase)

   ! perform IFFT to find time history from broadband spectrum
   trasf_sign=-1
   call four1d(broad_complex,trasf_sign)

   ! output broadband time-series (after scaling)
   bb_seis(:,comp)=real(broad_complex/npts)

endif

!!! ----- end frequency domain merging -----

END SUBROUTINE bb_calc

FUNCTION acc_spec(f,station)
!------------------------------------------------------------------------------------
!
! Description:
!
!   Function, calculates the spectral acceleration at frequency f.
!   To be used to scale the HFs, rather than scaling to LFs.
!
! Dependencies:
!
!   Subroutine bb_calc  
!
! References:
!
!   Graves and Pitarka (2010) in BSSA
!
! Authors: K. Olsen
!
! Modified: December 2012 (v1.4.1)
!
! Updated: March 2013 (v1.4.2)
!   Use information from velocity model file and srf file.
!
! Updated: November 2013 (v1.5.3)
!   Change alt calculation.
!
! Updated: March 2014 (v1.5.4.1)
!   Add sdec factor related to cell-receiver distances (R_cell).
!
! Updated: April 2014 (v1.5.5)
!   Change sdec factor for gs_flag=2 (CEUS events).
!
! Updated: September 2014 (v1.5.5.1)
!   Change sdec factor for gs_flag=2.
! 
! Updated: February 2015 (v1.5.5.3)
!   Updated sdec computation for gs_flag=2.
!
! Updated: May 2015 (V1.5.5.3, gfortran)
!   Use DtoR for trigonometric functions.
!
! Updated: May 2015 (V1.5.5.4)
!   Add alt computation for NGA-west2.
!
! Updated: June 2015 (V1.6)
!   Move alt computation into srf_read subroutine in source.f90.
!
! Updated: February 2019 (v2.0)
!   Avoid theta = 0 for t_star computation for imerg=1 and imerg=2.
!
use constants; use def_kind; use flags; use scattering
use source_receiver; use waveform; use earthquake
use fault_area
use vel_model; use geometry

implicit none

! function                                                                                                                   
real(kind=r_single)                             :: acc_spec
! number of station
integer(kind=i_single)                          :: station
! integers
integer(kind=i_single)                          :: i,j,hyp_n,k
! frequency
real(kind=r_single)                             :: f
! corner frequency
real(kind=r_single)                             :: f0
! alt factor G&P 2010
!real(kind=r_single)                             :: alt
! source dimension
real(kind=r_single)                             :: w
! Radiation factor
real(kind=r_single)                             :: rad
! Moments
real(kind=r_single)                             :: Mo
! Frankel (1995) factor
real(kind=r_single)                             :: F1
! source-station distance
real(kind=r_single)                             :: R
! radiation scale factor
real(kind=r_single)                             :: C
! source radiation spectrum
real(kind=r_single)                             :: S
! high-frequency spectral decay
real(kind=r_single)                             :: P

! ray-distance in a layer, travel time in a layer
real(kind=r_single)                             :: rayin,tin
! angle between the surface and raypath for each subfault in degree
real(kind=r_single)                             :: theta
! impedance effects, path term, t* effects 
real(kind=r_single)                             :: imp,G,attn
! travel time from hypo to a station [sec]
real(kind=r_single)                             :: t_star1
! source radiation spectrum, rupture velocity, 
real(kind=r_single)                             :: Sref,vr_sub1,Si_sum
! scaling factor
real(kind=r_single)                             :: sca
! parameters in radiation factor
real(kind=r_single)                             :: toa,dp,lam,str,rdm
! parameters in radiation factor
real(kind=r_single)                             :: sv1,sv2,sv3,sh1,sh2
! temp parameters
real(kind=r_single)                             :: t_vsd_ave
! temp parameter for nsub to be real
real(kind=r_single)                             :: rnsub
! scaling factor related to number of subfaults
real(kind=r_single)                             :: subf_sca
! scaling factor related to R_cell
real(kind=r_single)                             :: sdec
! qk factor from G&P (2010) eq 15
real(kind=r_single),allocatable,dimension(:)    :: qk
! each subfault dimention, sqrt(area of each subfault) [cm] 
real(kind=r_single),allocatable,dimension(:)    :: wsub
! distance between each subfault to a station [cm]
real(kind=r_single),allocatable,dimension(:)    :: R_cell
! t* is sum(travel time in a layer/Qs in a layer) [sec] 
real(kind=r_single),allocatable,dimension(:)    :: t_star
! rupture speed at each subfault
real(kind=r_single),allocatable,dimension(:)    :: vr_sub
! degree to radian
real(kind=r_single),parameter                   :: DtoR=pi/180.0

!------------------------------------------------------------------------------------


!---------- common imerg=1,2 ---------------------------------------------------

Mo=10**(1.5*Mw+9.05)*1.e7
t_vsd_ave=vsd_ave*1.e5 ! [km/s] -> [cm/s]

!alpha_T, eqn (9) in G&P(2010)
!if (dip > 45..and.dip < 60.) then
!   alt=0.28+0.012*dip
!elseif (dip <= 45.) then
!   alt=0.82
!elseif (dip >= 60.) then
!   alt=1.
!endif

! new alt calculation (v1.5.3)
!if (dip > 55..and. dip < 70.) then
!   alt=0.16+0.012*dip
!elseif (dip <= 55.) then
!   alt=0.82
!elseif (dip >= 70.) then
!   alt=1.
!endif
!print*,'alt in composition.f90= ',alt

!compute qk factor from G&P (2010) eq 15
if(.not.allocated(qk)) allocate(qk(n_lay))
do i=1,n_lay
   qk(i)=afac+bfac*vs(i)
enddo

!rupture velocity, eqn (4) in G&P(2010)
if (hyp_z > 5..and.hyp_z < 8.) then
   vr_sub1=0.16+0.08*hyp_z*vs_avef
elseif (hyp_z <= 5.) then
   vr_sub1=0.56*vs_avef
elseif (hyp_z >= 8.) then
   vr_sub1=0.8*vs_avef
endif

! start computing t_star1
t_star1=0
   
! angle between the surface and raypath from each subfault
theta=asin(hyp_z/sr_hypo(station))

do j=2, n_lay
   if (hyp_z > depth(j)) then
      rayin=(depth(j)-depth(j-1))/sin(theta) ! travel distance[km]=thickness/sin(theta)
      tin=rayin/vs(j-1) ! Swave travel time in a layer [sec]
      ! using qk factor from  G&P (2010) eqn 15 instead of Qs
      t_star1=t_star1+(tin/qk(j-1)) ! t*=sum(time/qk) [sec]

   elseif (hyp_z < depth(j) .and. hyp_z >= depth(j-1)) then
      if (theta .eq. 0.0) then
         rayin=sr_hypo(station) ! avoid using 0-theta when hyp_z = 0
      else
         rayin=(hyp_z-depth(j-1))/sin(theta)
      endif
      tin=rayin/vs(j-1)
      ! using qk factor from  G&P (2010) eqn 15 instead of Qs
      t_star1=t_star1+(tin/qk(j-1)) ! t*=sum(time/qk) [sec]
      exit
   endif
enddo
! end computing t_star1

!---------- 1 subfault part, but also computed for both imerg=1,2 --------------

!radiation factor from G&P (2010) - should be ~0.7 or so.
rad=0.
do i=1,1000
   call random_number(rdm) !should be 0-1
   toa=rdm*90.*DtoR
   lam=360.*rdm*DtoR 
   dp=90.*rdm*DtoR 
   str=(rdm-0.5)*90.*DtoR
   sh1=cos(lam)*cos(dp)*cos(toa)*sin(str)+cos(lam)*sin(dp)*sin(toa)*cos(2.*str)
   sh2=sin(lam)*cos(2.*dp)*cos(toa)*cos(str)-0.5*sin(lam)*sin(2.*dp)*sin(toa)*sin(2.*str)
   sv1=sin(lam)*cos(2.*dp)*cos(2.*toa)*sin(str)-cos(lam)*cos(dp)*cos(2.*toa)*cos(str)
   sv2=0.5*cos(lam)*sin(dp)*sin(2.*toa)*sin(2.*str)
   sv3=-0.5*sin(lam)*sin(2.*dp)*sin(2.*toa)*(1.+(sin(str))**2)
   rad=rad+abs(sh1+sh2)+abs(sv1+sv2+sv3)
enddo
rad=rad/1000.

w=sqrt(sum_area)
!if (gs_flag .eq. 1) then
if (gs_flag == 1 .or. gs_flag == 3) then 
   f0=(2.1*vr_sub1*1e5)/(alt*pi*w)
elseif (gs_flag .eq. 2) then
   f0=(1.1*vr_sub1*1e5)/(alt*pi*w)
endif
F1=Mo/(str_fac*(w**3))
R=sr_hypo(station)*100000.
C=rad*2./(4.*pi*d_avef*(vs_avef*1e5)**3)
Sref=Mo*(2.*pi*f)**2/(1.+F1*(f/f0)**2)
P=exp(-pi*kappa(1)*f)
imp=sqrt((srcVs*srcR)/(vs_avef*d_avef)) !gross impedance effects, G&P (2010) eqn (14)
attn=exp(-pi*t_star1) !t* effects, G&P (2010), eqn (14)
G=imp*attn/R

acc_spec=fac*Sref*C*P*G
if (imerg.eq.1) return !done in acc_spec when subfault=1
!---------- END OF SUBFAULT=1 --------------------------------------------------


!---------- following for imerg=2 only -----------------------------------------

! subdividing into as many as subfaults as given in srf

! allocation
if(.not.allocated(wsub)) allocate(wsub(nsub),R_cell(nsub),t_star(nsub),vr_sub(nsub))

! sqrt(area of each subfault)
wsub=sqrt(areas) ! [cm]

do i=1,nsub
   ! compute station-subfault distances [km]
   ! xp(:),yp(:) are station coordinates [km]
   R_cell(i) &
      =sqrt((xp(station)-fx(i))**2 + (yp(station)-fy(i))**2 + (fz(i))**2)*100000. ! [cm]

   ! start computing t_star
   t_star(i)=0
   
   ! angle between the surface and raypath from each subfault
   theta=asin(fz(i)/R_cell(i)*100000.)

   do j=2, n_lay
      if (fz(i) > depth(j)) then
         rayin=(depth(j)-depth(j-1))/sin(theta) ! travel distance[km]=thickness/sin(theta)
         tin=rayin/vs(j-1) ! Swave travel time in a layer [sec]
         ! using qk factor from  G&P (2010) eqn 15 instead of Qs
         t_star(i)=t_star(i)+(tin/qk(j-1)) ! t*=sum(time/qk) [sec]

      elseif (fz(i) < depth(j) .and. fz(i) >= depth(j-1)) then
         if (theta .eq. 0.0) then   ! avoid 0-theta when fz(i)=0km
            rayin=R_cell(i)/100000. ! [km]
         else
            rayin=(fz(i)-depth(j-1))/sin(theta)
         endif
         tin=rayin/vs(j-1)
         ! using qk factor from  G&P (2010) eqn 15 instead of Qs
         t_star(i)=t_star(i)+(tin/qk(j-1))
         exit
      endif
   enddo
   ! end computing t_star
 
enddo

!Si_sum=0.
acc_spec=0.
do k=1,nsub

!rupture velocity, eqn (4) in G&P(2010)
   if (fz(k).gt.5..and.fz(k).lt.8.) then
      vr_sub(k)=(0.16+0.08*fz(k))*vscu(k)
   elseif (fz(k).le.5.) then
      vr_sub(k)=0.56*vscu(k)
   elseif (fz(k).ge.8.) then
      vr_sub(k)=0.8*vscu(k)
   endif

!radiation factor from G&P (2010) - should be ~0.7 or so.
   rad=0.
   do i=1,1000
      call random_number(rdm) !should be 0-1
      toa=rdm*90.*DtoR
      lam=360.*rdm*DtoR 
      dp=90.*rdm*DtoR 
      str=(rdm-0.5)*90.*DtoR
      sh1=cos(lam)*cos(dp)*cos(toa)*sin(str)+cos(lam)*sin(dp)*sin(toa)*cos(2.*str)
      sh2=sin(lam)*cos(2.*dp)*cos(toa)*cos(str)-0.5*sin(lam)*sin(2.*dp)*sin(toa)*sin(2.*str)
      sv1=sin(lam)*cos(2.*dp)*cos(2.*toa)*sin(str)-cos(lam)*cos(dp)*cos(2.*toa)*cos(str)
      sv2=0.5*cos(lam)*sin(dp)*sin(2.*toa)*sin(2.*str)
      sv3=-0.5*sin(lam)*sin(2.*dp)*sin(2.*toa)*(1.+(sin(str))**2)
      rad=rad+abs(sh1+sh2)+abs(sv1+sv2+sv3)
   enddo
   rad=rad/1000.

   imp=sqrt(vscu(k)*rhcu(k)/t_vsd_ave) !gross impedance effects, G&P (2010) eqn (14)
   attn=exp(-pi*t_star(k)) !t* effects, G&P (2010), eqn (14)

    
   ! if (gs_flag .eq. 1) then
   if (gs_flag == 1 .or. gs_flag == 3) then
      f0=(2.1*vr_sub(k))/(alt*pi*wsub(k)) !subfault corner frequency, G&P (2010) eqn (13)
   elseif (gs_flag .eq. 2) then
      f0=(1.1*vr_sub(k))/(alt*pi*wsub(k)) !subfault corner frequency, G&P (2010) eqn (13)
   endif
   F1=Mo/(nsub*str_fac*((wsub(k))**3)) !Frankel (1995) factor, G&P (2010) eqn (12)
   C=rad*2./(4.*pi*(rhcu(k))*(vscu(k))**3) !radiation scale factor, G&P (2010) eqn (11)
   S=mo_ratio(k)*(Mo/nsub)*(2.*pi*f)**2/(1.+F1*(f/f0)**2) !radiation spectrum, G&P (2010) eqn (12)

   ! if (gs_flag .eq. 1) sdec = 1.0
   if (gs_flag == 1 .or. gs_flag == 3) sdec = 1.0
   ! updated (01/14/2015)
   if (gs_flag .eq. 2) then
      if (R_cell(k) .le. 16000000.) then
         sdec=1.00
      elseif (R_cell(k).gt.16000000.) then
         sdec=1.10
      endif
   endif
!   if (gs_flag .eq. 2) then
!      if (R_cell(k) .le. 20000000.) then
!         sdec=1.00
!      !elseif (R_cell(k).gt.20000000.) then
!      elseif (R_cell(k).gt.16000000.) then
!         sdec=1.10
!      endif
!   endif

!   if (gs_flag .eq. 2) then
!      if (R_cell(k) .lt. 7000000.) then ! if (R_cell(k) < 70km)
!         sdec=1.02
!      elseif (R_cell(k).ge.7000000..and.R_cell(k).lt.8000000.) then
!         sdec=(-0.004*R_cell(k)/100000.)+1.3
!      elseif (R_cell(k).ge.8000000..and.R_cell(k).lt.20000000.) then
!         sdec=0.98
!      elseif (R_cell(k).ge.20000000..and.R_cell(k).lt.21000000.) then
!         sdec=(0.012*R_cell(k)/100000.)-1.42
!      elseif (R_cell(k).ge.21000000.) then
!         sdec=1.10
!      endif
!   endif
   

   G=imp*attn/(R_cell(k))**sdec !path term, G&P (2010) eqn 14
   P=exp(-pi*kappa(1)*f) ! high-freq spectral decay, G&P (2010) eqn (16)
   acc_spec=acc_spec+S*C*P*G
   !Si_sum=Si_sum+S
enddo

!compute scaling factor from many subfaults (in ADDITION to factor 'fac')
!sca=Sref/Si_sum
rnsub=nsub
subf_sca=sqrt(rnsub)/nsub

acc_spec=acc_spec*subf_sca*fac

! deallocation
if(allocated(wsub)) deallocate(wsub,R_cell) 

END FUNCTION acc_spec
