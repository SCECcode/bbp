SUBROUTINE random_sequence    
!---------------------------------------------------------------------------
!
! Description:
!
!   Fulfill an array with random numbers to be utilized while computing
!   coda waves (see coda.f90). Necessary for correct OpenMP parallelization.
!
! Dependencies:
!
!   None.
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
! 
   
use def_kind; use scattering, only: iseed, random_array, nscat
use source_receiver, only: n_stat

implicit none

! local variables
integer(kind=i_single)                         :: seed_size
integer(kind=i_single),allocatable,dimension(:):: tmp_seed

!---------------------------------------------------------------------------

! allocate arrays for seeds and random numbers
if (.not.allocated(random_array)) allocate(random_array(4,nscat))

! initialize the generator with a specified seed number
call random_seed(size=seed_size)                       
if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
tmp_seed=iseed
iseed = iseed + 1
call random_seed(put=tmp_seed)

! generate random numbers using previously specified seed. Seed is the same for each receiver
call random_number(random_array)


END SUBROUTINE random_sequence

!===================================================================================================

FUNCTION rand_pdf(seed,val_min,val_max,pdf)
!-----------------------------------------------------------------
!
! Description:
! 
!   External function, returns random numbers drawn from normal
!   or uniform distribution
!
! Dependencies:
!
!   Functions normal_pdf and uniform_pdf
! 
! Note: 
!
!   'pdf' controls the pdf used (normal or uniform)   
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind; use interfaces, only: normal_pdf, uniform_pdf

implicit none

! function
real(kind=r_scat)                :: rand_pdf
! seed number (array shape required by intrinsic function)
integer(kind=i_single),intent(in):: seed
! PDF to be used
character(len=*),intent(in)      :: pdf 
! min-max or mean-std.dev. for selected distribution
real(kind=r_scat),intent(in)     :: val_min,val_max               

!----------------------------------------------------------------


! Pick-up a random number
if (pdf == 'norm') rand_pdf=normal_pdf(val_min,val_max,seed)
if (pdf == 'unif') rand_pdf=uniform_pdf(val_min,val_max,seed)


END FUNCTION rand_pdf


!===================================================================================================


FUNCTION uniform_pdf(a,b,seed)
!--------------------------------------------------------------------
!
! Description:
!
!   Returns a scaled pseudorandom real number drawn from a uniform
!   distribution with limits A and B.
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

! limits of the uniform distribution
real(kind=r_scat),intent(in)                   :: a,b
! seed used to initialize the random numbers generator (array shape required by intrinsic function)
integer(kind=i_single),intent(in)              :: seed
! variable where random number is stored (local)
real(kind=r_single)                            :: harv
! local variables
integer(kind=i_single)                         :: seed_size
integer(kind=r_single),allocatable,dimension(:):: tmp_seed
! function 
real(kind=r_scat)                              :: uniform_pdf

!--------------------------------------------------------------------

if (seed <= 0)  then
   write (*,*) ' '
   write (*,*) ' UNIFORM_PDF - Fatal error!'
   write (*,*) '  Input value of SEED <= 0.'
   stop
endif

! initialize the generator with a specified seed number
call random_seed(size=seed_size)
if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
tmp_seed=seed
call random_seed(put=tmp_seed)

! generate random number
call random_number(harv)

! compute number between A and B
uniform_pdf = a + (b - a) * harv

      
END FUNCTION uniform_pdf


!===================================================================================================


FUNCTION normal_pdf(a,b,seed)
!----------------------------------------------------------------------------------
! 
! Description:
!
!   Returns a scaled pseudonormal R4. The normal probability distribution function
!   (PDF) is sampled, with mean A and standard deviation B.
!
! Dependencies:
!
!   None
!
! References:
!
!   Numerical Recipes
!
! Author: W. Imperatori 
!
! Modified: January 2009 (v1.3)
!

use def_kind; use constants

implicit none

! mean A and standard deviation B
real(kind=r_scat),intent(in)                   :: a,b
! seed number (array shape required by intrinsic function)
integer(kind=i_single),intent(in)              :: seed
! dummies 
real(kind=r_single)                            :: r1,r2,rsq,x
! local variables
integer(kind=i_single)                         :: seed_size
integer(kind=r_single),allocatable,dimension(:):: tmp_seed
! functions
real(kind=r_scat)                              :: normal_pdf

!----------------------------------------------------------------------------------

if (seed <= 0)  then
   write (*,*) ' '
   write (*,*) ' UNIFORM_PDF - Fatal error!'
   write (*,*) '  Input value of SEED <= 0.'
   stop
endif

! initialize the generator with a specified seed number
call random_seed(size=seed_size)
if (.not.allocated(tmp_seed)) allocate(tmp_seed(seed_size))
tmp_seed=seed
call random_seed(put=tmp_seed)

! generate random number
call random_number(r1)

do
   call random_number(r1)
   call random_number(r2)
   r1=2.0*r1-1.0
   r2=2.0*r2-1.0
   rsq=r1**2+r2**2
   if (rsq > 0.0 .and. rsq < 1.0) exit
enddo

rsq=sqrt(-2.0*log(rsq)/rsq)   !Box-Muller transformation
x=r1*rsq

! compute normal with mean a and std.dev. b 
normal_pdf = a + b * x


END FUNCTION normal_pdf
