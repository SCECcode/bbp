SUBROUTINE infcorr(cnpts,dt,station)
!-----------------------------------------------------------------------------------
!
! Description:
!
!   Correlation using Kemp_*.bin
!
! Authors: N. Wang
!
! Modified: February 2019 (v2.0)
!
use interfaces, only: four1d
use constants 
use def_kind 
use waveform 
use read_Kemp
use flags

implicit none

! current component and station number
integer(kind=i_single),intent(in)                     :: cnpts,station
! time-step
real(kind=r_single),intent(in)                        :: dt
! flag and counters
integer(kind=i_single)                                :: trasf_sign,i,kk,jj,m,n,rr
! index of real target frequency and Nyquist frequency
integer(kind=i_single)                                :: targ_index,f_nyq_index
! array for complex time-series
complex(kind=r_single),dimension(cnpts)               :: hf_four_ns,hf_four_ew,hf_four_ud
complex(kind=r_single),dimension(cnpts)               :: broad_complex_ns,broad_complex_ew,broad_complex_ud
! arrays for high and low frequency spectra (amplitude and phase)
real(kind=r_single),dimension(cnpts)                  :: Am_hf_ns,Am_hf_ew,Am_hf_ud,Ph_hf_ns,Ph_hf_ew,Ph_hf_ud
!!real(kind=r_single),dimension(cnpts)                  :: EAS
real(kind=r_single),dimension(cnpts,3)                :: R
real(kind=r_single),dimension(3,3)                  :: C
real(kind=r_single)                                   :: randunifa,randunifb,wrsq
real(kind=r_single),allocatable,dimension(:)          :: corr_cut_ns,corr_cut_ew,corr_cut_ud,corr_ns,corr_ew,corr_ud
! actual target frequency, averaged spectra and their ratio, Nyquist frequency
!!real(kind=r_single)                                   :: f_nyq
! delta-f
real(kind=r_single)                                   :: df,A,sigma,sigma_low

integer(kind=i_single)                                :: seed_size,rseed
integer(kind=i_single),allocatable,dimension(:)       :: tmp_seed
integer, dimension(1:8)                               :: dt_seed

!-----------------------------------------------------------------------------------
print*,'start infcorr'

! computing frequency-delta
!df = 1 / (cnpts * dt)
df = 1 / (v_npts * dt)

! compute nyquist frequency (its index is npts/2 +1)
! (since npts is multiple of 2, f_nyq = {(npts+2)/2 -1)*df}
!! f_nyq = 0.5 * (1 / dt)
f_nyq_index = cnpts/2 + 1

! flag for FFT (i.e. from time to frequency)
trasf_sign=1

hf_four_ns=cmplx(bb_seis(:,1),0.0)
hf_four_ew=cmplx(bb_seis(:,2),0.0)
hf_four_ud=cmplx(bb_seis(:,3),0.0)

! calling function for FFT
call four1d(hf_four_ns,trasf_sign)
call four1d(hf_four_ew,trasf_sign)
call four1d(hf_four_ud,trasf_sign)

! extracts Amplitude and Phase spectra (HF)
Am_hf_ns=cabs(hf_four_ns)
Am_hf_ew=cabs(hf_four_ew)
Am_hf_ud=cabs(hf_four_ud)
!!EAS=sqrt((Am_hf_ns**2+Am_hf_ew**2)/2)
Ph_hf_ns=atan2(aimag(hf_four_ns),real(hf_four_ns))
Ph_hf_ew=atan2(aimag(hf_four_ew),real(hf_four_ew))
Ph_hf_ew=atan2(aimag(hf_four_ud),real(hf_four_ud))

!! log normal random number----------------------
! generate a different random number every time
do i = 1,3

   call random_seed(size=seed_size)
   if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
   !! set fixed random seed for each realization/component
   !!!!   tmp_seed = cseed + (station-1)*3 + comp
   !! set random seed according to system time (different at each run)
   call random_seed(get=tmp_seed)
   call DATE_AND_TIME(values=dt_seed)
   tmp_seed(seed_size)=dt_seed(8); tmp_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)

   call random_seed(put=tmp_seed)
   deallocate(tmp_seed)

   do rr=1,cnpts

      do
        do
           call random_number(randunifa)
           call random_number(randunifb)
           randunifa=2.0*randunifa-1.0
           randunifb=2.0*randunifb-1.0
           wrsq=randunifa**2+randunifb**2
           if (wrsq > 0.0 .and. wrsq < 1.0) exit
        enddo

        wrsq=sqrt(-2.0*log(wrsq)/wrsq)   !Box-Muller transformation
        R(rr,i)=randunifa*wrsq
        if ( abs(R(rr,i)) < 2) exit
      enddo

      if ( abs(R(rr,i)) > 2) then
         print*, 'large randnoise', R(rr,i)
      endif

   enddo

enddo

print*,'randnoise mean',sum(R)/size(R)


!! K*R get correlated random number ------------------------------------
!print*,'nk,nks in correlation.f90',nk,nks

!A=exp(-(sigma^2/2))*exp(K*randn(1,Ncorr)');
sigma=0.5 ! standard deviation of random number
sigma_low=0.1 ! standard deviation of random number for very low frequency part
!K=sigma*K

m=f_nyq_index
allocate(corr_cut_ns(m),corr_cut_ew(m),corr_cut_ud(m))

C(1,1)=1.0
C(1,2)=0.7
C(1,3)=0.7
C(2,1)=0.0
C(2,2)=0.7141
C(2,3)=0.2941
C(3,1)=0.0
C(3,2)=0.0
C(3,3)=0.6508
R(nks+1:nk+nks,:) = matmul ( matmul(Kemp,R(nks+1:nk+nks,:)), C )

!!R(1:m,:) = matmul(R(1:m,:),C)

!! get the exp of correlated random number
! ns
corr_cut_ns(1:nks) = exp(-(sigma_low**2/2))*exp(sigma_low*R(1:nks,1))           ! begining
corr_cut_ns(nks+1:m) = exp(-(sigma**2/2))*exp(sigma*R(nks+1:m,1))
! ew
corr_cut_ew(1:nks) = exp(-(sigma_low**2/2))*exp(sigma_low*2*R(1:nks,2))           ! begining
corr_cut_ew(nks+1:m) = exp(-(sigma**2/2))*exp(sigma*R(nks+1:m,2))
! ud
corr_cut_ud(1:nks) = exp(-(sigma_low**2/2))*exp(sigma_low*2*R(1:nks,3))           ! begining
corr_cut_ud(nks+1:m) = exp(-(sigma**2/2))*exp(sigma*R(nks+1:m,3))

!print*,'ns mean corr_cut(1:nks)',sum(corr_cut_ns(1:nks))/size(corr_cut_ns(1:nks))
!print*,'ew mean corr_cut(1:nks)',sum(corr_cut_ew(1:nks))/size(corr_cut_ew(1:nks))
!print*,'ns mean corr_cut(nk+nks+1:m)',sum(corr_cut_ns(nk+nks+1:m))/size(corr_cut_ns(nk+nks+1:m))
!print*,'ew mean corr_cut(nk+nks+1:m)',sum(corr_cut_ew(nk+nks+1:m))/size(corr_cut_ew(nk+nks+1:m))
!print*,'ns mean corr_cut(nks+1:nk+nks)',sum(corr_cut_ns(nks+1:nk+nks))/size(corr_cut_ns(nks+1:nk+nks))
!print*,'ns mean corr_cut(nks+1:nk+nks)',sum(corr_cut_ew(nks+1:nk+nks))/size(corr_cut_ew(nks+1:nk+nks))

print*,'mean corr_cut',sum(corr_cut_ns)/size(corr_cut_ns),sum(corr_cut_ew)/size(corr_cut_ew),sum(corr_cut_ud)/size(corr_cut_ud)

!!!! only apply random noise to high-frequency part, use 1 for low-frequency part
!!corr_cut(1:101) = 1    ! point 101 corresponding to 1Hz

! generate symmetric corr
allocate(corr_ns(cnpts),corr_ew(cnpts),corr_ud(cnpts))
do i=1,f_nyq_index
corr_ns(i)=corr_cut_ns(i)
corr_ew(i)=corr_cut_ew(i)
corr_ud(i)=corr_cut_ud(i)
enddo
do i=f_nyq_index+1,cnpts
corr_ns(i)=corr_cut_ns(cnpts-i+2)
corr_ew(i)=corr_cut_ew(cnpts-i+2)
corr_ud(i)=corr_cut_ud(cnpts-i+2)
enddo

print*,'be corr mean',sum(corr_ns)/size(corr_ns),sum(corr_ew)/size(corr_ew),sum(corr_ud)/size(corr_ud)

!! end K*R---------------------------------------

!!! end log normal random number----------------------

!! Multiply the exp of correlated random number with Fourier amplitude
Am_hf_ns=corr_ns*Am_hf_ns
Am_hf_ew=corr_ew*Am_hf_ew
Am_hf_ud=corr_ud*Am_hf_ud

! combine amplitude and phase spectra (whole broadband spectrum)
broad_complex_ns=Am_hf_ns*exp(zeta*Ph_hf_ns)
broad_complex_ew=Am_hf_ew*exp(zeta*Ph_hf_ew)
broad_complex_ud=Am_hf_ud*exp(zeta*Ph_hf_ud)

! perform IFFT to find time history from broadband spectrum
trasf_sign=-1
call four1d(broad_complex_ns,trasf_sign)
call four1d(broad_complex_ew,trasf_sign)
call four1d(broad_complex_ud,trasf_sign)

!! output broadband time-series
bb_seis(:,1)=real(broad_complex_ns/cnpts)
bb_seis(:,2)=real(broad_complex_ew/cnpts)
bb_seis(:,3)=real(broad_complex_ud/cnpts)

END SUBROUTINE infcorr
