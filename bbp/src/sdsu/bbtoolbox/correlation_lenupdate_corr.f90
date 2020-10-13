!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!! Inter-frequency Correlation

SUBROUTINE infcorr(cnpts,station)
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
! Update: March 2019 (v2.1)
!   Add cseed, for inter-frequency correlation random seed
!              if cseed = -1, seed defined by system-time
!              else, cseed defined by user and seed value increases 107 for each station

use interfaces, only: four1d
use scattering, only: cseed
use constants
use def_kind
use waveform
use read_correlation_files
use flags

implicit none

! current component and station number
integer(kind=i_single),intent(in)                     :: cnpts,station
!! time-step
!real(kind=r_single),intent(in)                        :: dt
! flag and counters
integer(kind=i_single)                                :: trasf_sign,i,kk,jj,n,rr
! index of real target frequency and Nyquist frequency
integer(kind=i_single)                                :: targ_index, m
! array for complex time-series
complex(kind=r_single),dimension(cnpts)               :: hf_four_ns,hf_four_ew,hf_four_ud
complex(kind=r_single),dimension(cnpts)               :: broad_complex_ns,broad_complex_ew,broad_complex_ud
! arrays for high and low frequency spectra (amplitude and phase)
real(kind=r_single),dimension(cnpts)                  :: Am_hf_ns,Am_hf_ew,Am_hf_ud,Ph_hf_ns,Ph_hf_ew,Ph_hf_ud,f
!!real(kind=r_single),dimension(cnpts)                  :: EAS
real(kind=r_single),dimension(cnpts/2+1,3)            :: R
real(kind=r_single),dimension(3,3)                    :: C
real(kind=r_single)                                   :: randunifa,randunifb,wrsq
real(kind=r_single),allocatable,dimension(:)          :: corr_cut_ns,corr_cut_ew,corr_cut_ud,corr_ns,corr_ew,corr_ud
! actual target frequency, averaged spectra and their ratio, Nyquist frequency
!!real(kind=r_single)                                   :: f_nyq
! delta-f
real(kind=r_single)                                   :: A,sigma,sigma_low

integer(kind=i_single)                                :: seed_size
integer(kind=i_single),allocatable,dimension(:)       :: tmp_seed
integer, dimension(1:8)                               :: dt_seed

!-----------------------------------------------------------------------------------
print*,'start infcorr'

!! computing frequency-delta
!!df = 1 / (cnpts * dt)

! compute nyquist frequency (its index is npts/2 +1)
! (since npts is multiple of 2, f_nyq = {(npts+2)/2 -1)*df}
!! f_nyq = 0.5 * (1 / dt)
m = cnpts/2 + 1

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

! set random seed (cseed) for compute random numbers
if (cseed == -1) then            !random seed based on system time
    call random_seed(size=seed_size)
    if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
    call random_seed(get=tmp_seed)
    call DATE_AND_TIME(values=dt_seed)
    tmp_seed(seed_size)=dt_seed(8); tmp_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
    print*, 'seed value for correlation random number tmp_seed(1) & cseed', tmp_seed(1), cseed
    !print*, 'seed size', seed_size
    !print*, 'tmp_seed', tmp_seed
    call random_seed(put=tmp_seed)
    deallocate(tmp_seed)
else
    call random_seed(size=seed_size)
    if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
    tmp_seed = cseed + 107*station  !set fixed random seed for each realization/station
    print*, 'seed value for correlation random number & cseed', tmp_seed(1), cseed
    !print*, 'seed size', seed_size
    !print*, 'tmp_seed', tmp_seed
    call random_seed(put=tmp_seed)
    deallocate(tmp_seed)
endif


do i = 1,3

   do rr=1,m !!cnpts

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

!A=exp(-(sigma^2/2))*exp(K*randn(1,Ncorr)')
sigma=0.5 ! standard deviation of random number
sigma_low=0.1 ! standard deviation of random number for very low frequency part
!K=sigma*K

C(1,1)=1.0
C(1,2)=0.7
C(1,3)=0.7
C(2,1)=0.0
C(2,2)=0.7141
C(2,3)=0.2941
C(3,1)=0.0
C(3,2)=0.0
C(3,3)=0.6508


!!!!R(1:m,:) = matmul(R(1:m,:),C)
R(nks+1:nk+nks,:) = matmul ( matmul(Kinf,R(nks+1:nk+nks,:)), C )

print*,'randnoise mean after',sum(R)/size(R)
!! get the exp of correlated random number
allocate(corr_cut_ns(m),corr_cut_ew(m),corr_cut_ud(m))

! ns
corr_cut_ns(1:nks) = exp(-(sigma_low**2/2))*exp(sigma_low*R(1:nks,1))           ! begining
corr_cut_ns(nks+1:nk+nks) = exp(-(sigma**2/2))*exp(sigma*R(nks+1:nk+nks,1))
corr_cut_ns(nks+nk+1:m) = exp(-(sigma_low**2/2))*exp(sigma_low*R(nks+nk+1:m,1))
! ew
corr_cut_ew(1:nks) = exp(-(sigma_low**2/2))*exp(sigma_low*R(1:nks,2))           ! begining
corr_cut_ew(nks+1:nk+nks) = exp(-(sigma**2/2))*exp(sigma*R(nks+1:nk+nks,2))
corr_cut_ew(nks+nk+1:m) = exp(-(sigma_low**2/2))*exp(sigma_low*R(nks+nk+1:m,2))
! ud
corr_cut_ud(1:nks) = exp(-(sigma_low**2/2))*exp(sigma_low*R(1:nks,3))           ! begining
corr_cut_ud(nks+1:nk+nks) = exp(-(sigma**2/2))*exp(sigma*R(nks+1:nk+nks,3))
corr_cut_ud(nks+nk+1:m) = exp(-(sigma_low**2/2))*exp(sigma_low*R(nks+nk+1:m,3))

!print*,'ns mean corr_cut(1:nks)',sum(corr_cut_ns(1:nks))/size(corr_cut_ns(1:nks))
!print*,'ew mean corr_cut(1:nks)',sum(corr_cut_ew(1:nks))/size(corr_cut_ew(1:nks))
!print*,'ns mean corr_cut(nk+nks+1:m)',sum(corr_cut_ns(nk+nks+1:m))/size(corr_cut_ns(nk+nks+1:m))
!print*,'ew mean corr_cut(nk+nks+1:m)',sum(corr_cut_ew(nk+nks+1:m))/size(corr_cut_ew(nk+nks+1:m))
!print*,'ns mean corr_cut(nks+1:nk+nks)',sum(corr_cut_ns(nks+1:nk+nks))/size(corr_cut_ns(nks+1:nk+nks))
!print*,'ew mean corr_cut(nks+1:nk+nks)',sum(corr_cut_ew(nks+1:nk+nks))/size(corr_cut_ew(nks+1:nk+nks))

print*,'mean corr_cut',sum(corr_cut_ns)/size(corr_cut_ns),sum(corr_cut_ew)/size(corr_cut_ew),sum(corr_cut_ud)/size(corr_cut_ud)

!!!! only apply random noise to high-frequency part, use 1 for low-frequency part
!!corr_cut(1:101) = 1    ! point 101 corresponding to 1Hz

! generate symmetric corr
allocate(corr_ns(cnpts),corr_ew(cnpts),corr_ud(cnpts))
do i = 1,m
corr_ns(i)=corr_cut_ns(i)
corr_ew(i)=corr_cut_ew(i)
corr_ud(i)=corr_cut_ud(i)
enddo
do i = m+1,cnpts
corr_ns(i)=corr_cut_ns(cnpts-i+2)
corr_ew(i)=corr_cut_ew(cnpts-i+2)
corr_ud(i)=corr_cut_ud(cnpts-i+2)
enddo

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






!-----------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
!! Spatial Correlation

SUBROUTINE sp_corr_random(cnpts)
!-----------------------------------------------------------------------------------
!!Corr using correlation matrices

use interfaces, only: cholesky
use constants 
use def_kind
use scattering, only: cseed
use waveform 
use read_correlation_files
use source_receiver, only: n_stat, xp, yp

implicit none

integer(kind=i_single),intent(in)                     :: cnpts
integer(kind=i_single)                                :: n,i,j,m,rr

! arrays for high and low frequency spectra (amplitude and phase)
real(kind=r_single),dimension(cnpts/2+1,n_stat,9)          :: randnoise
real(kind=r_single),dimension(cnpts/2+1,n_stat,3)          :: R1corr
real(kind=r_single),dimension(cnpts/2+1,n_stat,3)          :: R2corr
real(kind=r_single),dimension(cnpts/2+1,n_stat,3)          :: R3corr
real(kind=r_single),dimension(cnpts/2+1,n_stat)            :: R1ns, R2ns, R3ns, Rns
real(kind=r_single),dimension(cnpts/2+1,n_stat)            :: R1ew, R2ew, R3ew, Rew
real(kind=r_single),dimension(cnpts/2+1,n_stat)            :: R1ud, R2ud, R3ud, Rud
real(kind=r_double),dimension(n_stat,n_stat)              :: L1, L2, S1, S2
real(kind=r_double)                                       :: h
real(kind=r_single)                                   :: randunifa,randunifb,wrsq
real(kind=r_single),allocatable,dimension(:,:)        :: corr_cut_ns, corr_cut_ew, corr_cut_ud
real(kind=r_single)                                   :: sigma,sigma_low
real(kind=r_single),dimension(3,3)                    :: C

integer(kind=i_single)                                :: seed_size,rseed
integer(kind=i_single),allocatable,dimension(:)       :: tmp_seed
integer, dimension(1:8)                               :: dt_seed

!-----------------------------------------------------------------------------------
print*,'start sp_corr_random'
!! computing frequency-delta
!!df = 1 / (cnpts * dt)

! compute nyquist frequency (its index is cnpts/2 +1)
! (since cnpts is multiple of 2, f_nyq = {(cnpts+2)/2 -1)*df}
!f_nyq = 0.5 * (1 / dt)
m = cnpts/2 + 1 !!f_nyq_index
print*,'cnpts',cnpts
print*,'nks,nk',nks,nk

   !! log normal random number----------------------
   ! generate a different random number every time
! set random seed (cseed) for compute random numbers
if (cseed == -1) then            !random seed based on system time
    call random_seed(size=seed_size)
    if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
    call random_seed(get=tmp_seed)
    call DATE_AND_TIME(values=dt_seed)
    tmp_seed(seed_size)=dt_seed(8); tmp_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
    print*, 'seed value for correlation random number tmp_seed(1) & cseed', tmp_seed(1), cseed
    !print*, 'seed size', seed_size
    !print*, 'tmp_seed', tmp_seed
    call random_seed(put=tmp_seed)
    deallocate(tmp_seed)
else
    call random_seed(size=seed_size)
    if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
    tmp_seed = cseed   !set fixed random seed for each realization
    print*, 'seed value for correlation random number & cseed', tmp_seed(1), cseed
    !print*, 'seed size', seed_size
    !print*, 'tmp_seed', tmp_seed
    call random_seed(put=tmp_seed)
    deallocate(tmp_seed)
endif

   do i=1,9
     !print*,'i',i
     do n=1,n_stat
       do rr=1,m !!cnpts

          do
            call random_number(randunifa)
            call random_number(randunifb)
            randunifa=2.0*randunifa-1.0
            randunifb=2.0*randunifb-1.0
            wrsq=randunifa**2+randunifb**2
            if (wrsq > 0.0 .and. wrsq < 1.0) exit
          enddo

          wrsq=sqrt(-2.0*log(wrsq)/wrsq)   !Box-Muller transformation
          randnoise(rr,n,i)=randunifa*wrsq
       enddo
     enddo
   enddo

   !! 3 independent random matrix
   R1corr = randnoise(:,:,1:3)
   R2corr = randnoise(:,:,4:6)
   R3corr = randnoise(:,:,7:9)
   
   !! 0.7 correlation between FAS components
   C(1,1)=1.0
   C(1,2)=0.7
   C(1,3)=0.7
   C(2,1)=0.0
   C(2,2)=0.7141
   C(2,3)=0.2941
   C(3,1)=0.0
   C(3,2)=0.0
   C(3,3)=0.6508

   R1ns = R1corr(:,:,1)*C(1,1) + R1corr(:,:,2)*C(2,1) + R1corr(:,:,3)*C(3,1)
   R1ew = R1corr(:,:,1)*C(1,2) + R1corr(:,:,2)*C(2,2) + R1corr(:,:,3)*C(3,2)
   R1ud = R1corr(:,:,1)*C(1,3) + R1corr(:,:,2)*C(2,3) + R1corr(:,:,3)*C(3,3)
   
   R2ns = R2corr(:,:,1)*C(1,1) + R2corr(:,:,2)*C(2,1) + R2corr(:,:,3)*C(3,1)
   R2ew = R2corr(:,:,1)*C(1,2) + R2corr(:,:,2)*C(2,2) + R2corr(:,:,3)*C(3,2)
   R2ud = R2corr(:,:,1)*C(1,3) + R2corr(:,:,2)*C(2,3) + R2corr(:,:,3)*C(3,3)

   R3ns = R3corr(:,:,1)*C(1,1) + R3corr(:,:,2)*C(2,1) + R3corr(:,:,3)*C(3,1)
   R3ew = R3corr(:,:,1)*C(1,2) + R3corr(:,:,2)*C(2,2) + R3corr(:,:,3)*C(3,2)
   R3ud = R3corr(:,:,1)*C(1,3) + R3corr(:,:,2)*C(2,3) + R3corr(:,:,3)*C(3,3)
   
   !print*,'R1corr(10:12,5:6,:)', R1corr(10:12,5:6,:)
   !print*,'R1ns(10:12,5:6)', R1ns(10:12,5:6)
   !print*,'R1ew(10:12,5:6)', R1ew(10:12,5:6)
   !print*,'R1ud(10:12,5:6)', R1ud(10:12,5:6)
   !print*,'R2corr(10:12,5:6,:)', R2corr(10:12,5:6,:)
   !print*,'R2ns(10:12,5:6)', R2ns(10:12,5:6)
   !print*,'R2ew(10:12,5:6)', R2ew(10:12,5:6)
   !print*,'R2ud(10:12,5:6)', R2ud(10:12,5:6)
   !print*,'R3corr(10:12,5:6,:)', R3corr(10:12,5:6,:)
   !print*,'R3ns(10:12,5:6)', R3ns(10:12,5:6)
   !print*,'R3ew(10:12,5:6)', R3ew(10:12,5:6)
   !print*,'R3ud(10:12,5:6)', R3ud(10:12,5:6)
   !! K1*R1*L1 + K2*R2*L2 + K3*R3---------------------------

   !A=exp(-(sigma^2/2))*exp(K*randn(1,Ncorr)');
   sigma=0.5 ! standard deviation of random number
   sigma_low=0.2 ! standard deviation of random number for very low frequency part
   !K=sigma*K
   
   !! Generate distance correlation matrix
   do i = 1,(n_stat-1)
       S1(i,i) = 1.0
       S2(i,i) = 1.0
       do j = i+1, n_stat
           h = sqrt((xp(i)-xp(j))**2 + (yp(i)-yp(j))**2)
           S1(i,j) = exp(-3.0/10.0*h)
           S1(j,i) = S1(i,j)
           S2(i,j) = exp(-3.0/100.0*h)
           S2(j,i) = S2(i,j)
       enddo
   enddo
   S1(n_stat,n_stat) = 1.0
   S2(n_stat,n_stat) = 1.0
   

   ! Cholesky decomposition
   call cholesky(n_stat,S1,L1)
   call cholesky(n_stat,S2,L2)
   
   L1 = transpose(L1)
   L2 = transpose(L2)


   !!! ns
      allocate(corr_cut_ns(m,n_stat))

      R1ns(nks+1:nk+nks,:) = matmul( matmul(Ksp1,R1ns(nks+1:nk+nks,:)), L1)
      R2ns(nks+1:nk+nks,:) = matmul( matmul(Ksp2,R2ns(nks+1:nk+nks,:)), L2)
      R3ns(nks+1:nk+nks,:) = matmul(Ksp3,R3ns(nks+1:nk+nks,:))

      Rns = R1ns + R2ns + R3ns
      print*,'after KL R1ns(10:12,5:6)',R1ns(10:12,5:6)
      print*,'Rns mean',sum(Rns)/size(Rns)

      !! get the exp of correlated random number
      corr_cut_ns(1:nks,:) = exp(-(sigma_low**2/2))*exp(sigma_low*Rns(1:nks,:)/3)          ! begining
      corr_cut_ns(nks+1:nk+nks,:) = exp(-(sigma**2/2))*exp(sigma*Rns(nks+1:nk+nks,:))
      corr_cut_ns(nk+nks+1:m,:) = exp(-(sigma_low**2/2))*exp(sigma_low*Rns(nk+nks+1:m,:)/3)
      
      !!! ew
      allocate(corr_cut_ew(m,n_stat))

      R1ew(nks+1:nk+nks,:) = matmul( matmul(Ksp1,R1ew(nks+1:nk+nks,:)), L1)
      R2ew(nks+1:nk+nks,:) = matmul( matmul(Ksp2,R2ew(nks+1:nk+nks,:)), L2)
      R3ew(nks+1:nk+nks,:) = matmul(Ksp3,R3ew(nks+1:nk+nks,:))

      Rew = R1ew+ R2ew + R3ew
      print*,'after KL R1ew(10:12,5:6)',R1ew(10:12,5:6)
      print*,'Rew mean',sum(Rew)/size(Rew)

      !! get the exp of correlated random number
      corr_cut_ew(1:nks,:) = exp(-(sigma_low**2/2))*exp(sigma_low*Rew(1:nks,:)/3) ! begining
      corr_cut_ew(nks+1:nk+nks,:) = exp(-(sigma**2/2))*exp(sigma*Rew(nks+1:nk+nks,:))
      corr_cut_ew(nk+nks+1:m,:) = exp(-(sigma_low**2/2))*exp(sigma_low*Rew(nk+nks+1:m,:)/3)

      !!! ud
      allocate(corr_cut_ud(m,n_stat))

      R1ud(nks+1:nk+nks,:) = matmul( matmul(Ksp1,R1ud(nks+1:nk+nks,:)), L1)
      R2ud(nks+1:nk+nks,:) = matmul( matmul(Ksp2,R2ud(nks+1:nk+nks,:)), L2)
      R3ud(nks+1:nk+nks,:) = matmul(Ksp3,R3ud(nks+1:nk+nks,:))

      Rud = R1ud+ R2ud + R3ud
      print*,'after KL R1ud(10:12,5:6)',R1ud(10:12,5:6)
      print*,'Rud mean',sum(Rud)/size(Rud)

      !! get the exp of correlated random number
      corr_cut_ud(1:nks,:) = exp(-(sigma_low**2/2))*exp(sigma_low*Rud(1:nks,:)/3) ! begining
      corr_cut_ud(nks+1:nk+nks,:) = exp(-(sigma**2/2))*exp(sigma*Rud(nks+1:nk+nks,:))
      corr_cut_ud(nk+nks+1:m,:) = exp(-(sigma_low**2/2))*exp(sigma_low*Rud(nk+nks+1:m,:)/3)


      !!!! only apply random noise to high-frequency part, use 1 for low-frequency part
      !!corr_cut(1:101) = 1    ! point 101 corresponding to 1Hz

      ! generate symmetric corr_res
      do i=1,m
         corr_res(i,:,1)=corr_cut_ns(i,:)
         corr_res(i,:,2)=corr_cut_ew(i,:)
         corr_res(i,:,3)=corr_cut_ud(i,:)
      enddo

      do i=m+1,cnpts
         corr_res(i,:,1)=corr_cut_ns(cnpts-i+2,:)
         corr_res(i,:,2)=corr_cut_ew(cnpts-i+2,:)
         corr_res(i,:,3)=corr_cut_ud(cnpts-i+2,:)
      enddo

   !print*,'1, corr_res(10:12,1,1)',corr_res(10:12,1,1)
   !print*,'2, corr_res(10:12,1,2)',corr_res(10:12,1,2)
   !print*,'3, corr_res(10:12,1,3)',corr_res(10:12,1,3)


   print*,'corr 1 mean',sum(corr_res(:,1,1))/size(corr_res(:,1,1))
   print*,'corr 2 mean',sum(corr_res(:,1,2))/size(corr_res(:,1,2))
   print*,'corr 3 mean',sum(corr_res(:,1,3))/size(corr_res(:,1,3))

      !! end K*R---------------------------------------

      !!! end log normal random number----------------------


END SUBROUTINE sp_corr_random




SUBROUTINE add_fourier_residual(station,comp,cnpts)
!-----------------------------------------------------------------------------------

use interfaces, only: four1d
use constants
use def_kind
use waveform
!use source_receiver, only: n_stat
!use scattering, only: cnpts

implicit none

! current component and station number
integer(kind=i_single),intent(in)                     :: station,comp,cnpts

! flag and counters
integer(kind=i_single)                                :: trasf_sign
! array for complex time-series
complex(kind=r_single),dimension(cnpts)                :: hf_four,broad_complex
! arrays for high and low frequency spectra (amplitude and phase)
real(kind=r_single),dimension(cnpts)                   :: Am_hf,Ph_hf

!-----------------------------------------------------------------------------------
print*,'start add_fourier_residual for sp_corr'
! flag for FFT (i.e. from time to frequency)
trasf_sign=1

hf_four=cmplx(bb_seis(:,comp),0.0)

! calling function for FFT
call four1d(hf_four,trasf_sign)

! extracts Amplitude and Phase spectra (HF)
Am_hf=cabs(hf_four)
Ph_hf=atan2(aimag(hf_four),real(hf_four))

!! Multiply the exp of correlated random number with Fourier amplitude
Am_hf=corr_res(:,station,comp)*Am_hf

! combine amplitude and phase spectra (whole broadband spectrum)
broad_complex=Am_hf*exp(zeta*Ph_hf)

! perform IFFT to find time history from broadband spectrum
trasf_sign=-1
call four1d(broad_complex,trasf_sign)

!! output broadband time-series
bb_seis(:,comp)=real(broad_complex/cnpts)


END SUBROUTINE add_fourier_residual

!===================================================================================================

SUBROUTINE cholesky(n,A,G)

use def_kind

implicit none

! formal vars
integer(kind=i_single),intent(in)    :: n
real(kind=r_double),intent(in)    :: A(n,n)
real(kind=r_double),intent(out)   :: G(n,n)
! local vars
integer(kind=i_single) :: i,j      ! iteration counter

! begin loop
G(:,:) = 0.0
  do j = 1, n
     G(j,j) = sqrt( A(j,j) - dot_product(G(j,1:j-1),G(j,1:j-1)) )
     do i = j+1, n
        G(i,j)  = ( A(i,j) - dot_product(G(i,1:j-1),G(j,1:j-1)) ) / G(j,j)
     end do
  end do

END SUBROUTINE cholesky


