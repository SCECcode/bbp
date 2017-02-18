FUNCTION convlv(data,respns)
!-----------------------------------------------------------------------------
! 
! Description:
!
!   External function convolving two time-series. Convolution is perfomed in 
!   the frequency domain. Unit time-step is assumed.
!
! Dependencies:
!
!   Subroutine realft
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1253
!
! Notes:
!
!   Respns is assumed to have LESS number of points than Data.
!   The Respns time-series is assumed to NOT have values at negative times.
!   Input time-series are both zero-padded to avoid side-effects. Convolution
!   result is then truncated up to the original number of points before being
!   passed.
!   
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind; use interfaces, only: realft

implicit none

! first input time-series (actually it is NOT modified)
real(kind=r_single),dimension(:),intent(inout):: data
! second input time-series
real(kind=r_single),dimension(:),intent(in)   :: respns
! output function
real(kind=r_single),dimension(size(data))     :: convlv
! local variables
integer(kind=i_single)                        :: n,m
complex(kind=r_single),dimension(size(data))  :: tmpd,tmpr
real(kind=r_single),dimension(2*size(data))   :: pad1,pad2

!-----------------------------------------------------------------------------

! extract sizes	
n=size(data)
m=size(respns)           

! check m > n
if (m > n) call error_handling(9,'null','CONVLV (fourier.f90)')

! check for 2**N number of points
if ( iand(n,n-1) /=0 ) call error_handling(5,'null','CONVLV (fourier.f90)')

! zero-pad input data to avoid side-effects (conservative assumption: zero-padding
! equal to first input length - in theory it could be equal to non-zero length
! of second input vector)
pad1(1:n) = data; pad1(n+1:) = 0. 

! zero-pad also respns to have same number of points in frequency
pad2(1:m)=respns(:); pad2(m+1:)=0. 

! FFT of input time-series
call realft(pad1,1,tmpd)
call realft(pad2,1,tmpr)

! multiply FFTs
tmpr(1)=cmplx(real(tmpd(1))*real(tmpr(1))/n,aimag(tmpd(1))*aimag(tmpr(1))/n,kind=r_single)
tmpr(2:)=tmpd(2:)*tmpr(2:)/n

! IFFT of multiplied Fourier spectra	
call realft(pad2,-1,tmpr)

! return only first n points (exclude padding values)
convlv = pad2(1:n)

END FUNCTION convlv

!===================================================================================================

FUNCTION correl(data1,data2)
!-------------------------------------------------------------------------------
!
! Description:
!
!   Compute correlation between two time series in the frequency domain
!
! Dependencies:
!
!   Subroutine realft
!
! Notes:
!
!   Time-series are not padded here. They should be passed already padded to 
!   avoid side-effects
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1254
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3) 
!
! *** NOT USED IN THIS VERSION ***
!

use def_kind; use interfaces, only: realft

implicit none

! input arrays
real(kind=r_single),dimension(:),intent(in)   :: data1,data2
! function
real(kind=r_single),dimension(size(data1))    :: correl
! local variables
complex(kind=r_single),dimension(size(data1)) :: cdat1,cdat2
integer(kind=i_single)                        :: n
real(kind=r_single),dimension(2*size(data1)) :: pad1,pad2

!------------------------------------------------------------------------------	

n = size(data1)

! checks for equal series length and 2**N npts 
if (n /= size(data2)) call error_handling(9,'null','CORREL (fourier.f90)')   
if (iand(n,n-1) /= 0) call error_handling(5,'null','CORREL (fourier.f90)')

! zero-pad input data to avoid side-effects
pad1(1:n) = data1; pad1(n+1:) = 0. 
pad2(1:n) = data2; pad2(n+1:) = 0. 

! FFT
call realft(pad1,1,cdat1)
call realft(pad2,1,cdat2)

! perform multiplication with complex conjugate
cdat1(1)=cmplx(real(cdat1(1))*real(cdat2(1))/n,aimag(cdat1(1))*aimag(cdat2(1))/n,kind=r_single)
cdat1(2:)=cdat1(2:)*conjg(cdat2(2:))/n

! IFFT
call realft(pad1,-1,cdat1)

! return only N points (exclude padding values)
correl(1:n/2) = pad1(1:n/2); correl(n/2 +1:) = pad1(3*n/2 +1:)

END FUNCTION correl

!===================================================================================================

FUNCTION lag_correl(seq1,seq2)
!--------------------------------------------------------------------
!
! Description:
!
!   Evaluate at which lag value the correlation of two sequences has 
!   its maximum
!
! Dependencies:
!
!   Function correl
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!
! *** NOT USED IN THIS VERSION ***
!

use def_kind; use interfaces, only: correl

implicit none

! input time-series
real(kind=r_single),dimension(:),intent(in):: seq1,seq2
! function
real(kind=r_single)                        :: lag_correl
! local variables
integer(kind=i_single)                     :: max_index,n

!--------------------------------------------------------------------

n = size(seq1)

max_index = maxloc(correl(seq1,seq2),dim=1)

if (max_index > n/2) then
   lag_correl = -(max_index - n - 1) 
else  
   lag_correl = (max_index - 1) 
endif
 
END FUNCTION lag_correl

!===================================================================================================

SUBROUTINE four1d(data,isign)
!---------------------------------------------------------------------
!
! Description:
!
!   External subroutine from Numerical Recipes, prepares input for 
!   FAST FOURIER TRANSFORM (FFT) or inverse (IFFT) computation
!
! Dependencies:
!
!   Subroutine fourrow, function arth
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1239-1240
!
! Notes:
!
!   i_sign=1 -> FFT, i_sign=-1 -> IFFT
!
! Warning:
!
!   This subroutine works ONLY with 2**n number of points
!
! Authors: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use constants; use def_kind; use interfaces, only: fourrow, arth

implicit none

! input-output vector where complex transform will be stored
complex(kind=r_single),dimension(:),intent(inout) :: data
! flag for FFT/IFFT
integer(kind=i_single),intent(in)                 :: isign
! locals
complex(kind=r_single),dimension(:,:),allocatable :: dat,temp
complex(kind=r_double),dimension(:),allocatable   :: w,wp
real(kind=r_double),dimension(:),allocatable      :: theta
integer(kind=i_single)                            :: n,m1,m2,j

!-----------------------------------------------------------------------

n=size(data)                            

if ( iand(n,n-1) /=0 ) call error_handling(5,'null','four1d (fft.f90)')

m1=2**ceiling(0.5*log(real(n,kind=r_single))/0.693147)  
m2=n/m1

allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))

dat=reshape(data,shape(dat))

call fourrow(dat,isign)    !first transform   

theta=arth(0,isign,m1)*pi_double/n
wp=cmplx(-2.0*sin(0.5*theta)**2,sin(theta),kind=r_double)
w=cmplx(1.0,0.0,kind=r_double)

do j=2,m2
   w=w*wp+w
   dat(:,j)=dat(:,j)*w
enddo

temp=transpose(dat)

call fourrow(temp,isign)   !second transform
                                           
data=reshape(temp,shape(data))
                                           
deallocate(dat,w,wp,theta,temp)

END SUBROUTINE four1d

!===================================================================================================

SUBROUTINE fourrow(data,isign)
!---------------------------------------------------------------------
!
! Description:
!
!   External subroutine from Numerical Recipes, calculates FAST FOURIER
!   TRANSFORM (FFT) and its inverse (IFFT) for single-precision complex
!   data
!
! Dependencies:
!
!   Subroutine swap
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1235
!
! Notes:
!
!   i_sign=1 -> FFT, i_sign=-1 -> IFFT
!
! Warning:
!
!   This subroutine works ONLY with 2**n number of points
!
! Authors: W. Imperatori, M. Mai, B. Mena
!
! Modified: January 2009 (v1.3)
!

use constants; use def_kind; use interfaces, only: swap

implicit none

! input/output array
complex(kind=r_single),dimension(:,:),intent(inout) :: data
! flag for FFT/IFFT
integer(kind=i_single),intent(in)                   :: isign
! locals
integer(kind=i_single)                              :: n,i,istep,j,m,mmax,n2
real(kind=r_double)                                 :: theta
complex(kind=r_single), dimension(size(data,1))     :: temp
complex(kind=r_double)                              :: w,wp
complex(kind=r_single)                              :: ws

!---------------------------------------------------------------------

n=size(data,2)                      
                 
if ( iand(n,n-1) /= 0 ) call error_handling(5,'null','fourrow (fft.f90)')

n2=n/2

j=n2

do i=1,n-2
   if (j > i) call swap(data(:,j+1),data(:,i+1))
   m=n2

   do
      if (m < 2 .or. j < m) exit
         j=j-m
         m=m/2
   enddo
 
   j=j+m

enddo

mmax=1

do
   if (n <= mmax) exit
   istep=2*mmax
   theta=pi/(isign*mmax)
   wp=cmplx(-2.0*sin(0.5*theta)**2,sin(theta),kind=r_double)
   w=cmplx(1.0,0.0,kind=r_double)

   do m=1,mmax
      ws=w

      do i=m,n,istep
         j=i+mmax
         temp=ws*data(:,j)
         data(:,j)=data(:,i)-temp
         data(:,i)=data(:,i)+temp
      enddo
      w=w*wp+w
   enddo

   mmax=istep

enddo

END SUBROUTINE fourrow

!===================================================================================================

SUBROUTINE realft(data,isign,zdata)
!-------------------------------------------------------------------
!
! Description:
!
!   External subroutine performing FFT and IFFT of real time-series 
!
! Dependencies:
!
!   subroutine four1d, function zroots_unity
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1243
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind; use interfaces, only: four1d, zroots_unity

implicit none

real(kind=r_single),dimension(:),intent(inout)     :: data
integer(kind=i_single),intent(in)                  :: isign
complex(kind=r_single),dimension(:),optional,target:: zdata
! local variables
integer(kind=i_single)                             :: n,ndum,nh,nq
complex(kind=r_single),dimension(size(data)/4)     :: w
complex(kind=r_single),dimension(size(data)/4-1)   :: h1,h2
complex(kind=r_single),dimension(:),pointer        :: cdata
complex(kind=r_single)                             :: z
real(kind=r_single)                                :: c1=0.5,c2

!------------------------------------------------------------------

n=size(data)                              

if ( iand(n,n-1) /= 0 ) call error_handling(5,'null','realft (fft.f90)')

nh=n/2
nq=n/4

if (present(zdata)) then

   if ( n/2 /= size(zdata) ) call error_handling(7,'null','realft (fft.f90)')
   
   ndum=size(zdata)
   cdata=>zdata
   if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=r_single)
else
   allocate(cdata(n/2))
   cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=r_single)
endif

if (isign == 1) then
   c2=-0.5                             
   call four1d(cdata,+1)
else
   c2=0.5
endif
                                                  
w=zroots_unity(sign(n,isign),n/4)              
w=cmplx(-aimag(w),real(w),kind=r_single)
h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
cdata(2:nq)=h1+w(2:nq)*h2
cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
z=cdata(1)

if (isign == 1) then
   cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=r_single)
else
   cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=r_single)
   call four1d(cdata,-1)                     
endif

if (present(zdata)) then
   if (isign /= 1) then
      data(1:n-1:2)=real(cdata)
      data(2:n:2)=aimag(cdata)
   endif
else
   data(1:n-1:2)=real(cdata)
   data(2:n:2)=aimag(cdata)
   deallocate(cdata)
endif

END SUBROUTINE realft

!===================================================================================================

FUNCTION arth(first,increment,n)
!-------------------------------------------------------
!
! Description:
!
!   Return an arithmetic progression
!
! Dependencies:
!
!   None
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1371
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

integer(kind=i_single),intent(in)  :: first,increment,n
integer(kind=i_single),dimension(n):: arth
! local variables
integer(kind=i_single)             :: k,k2,temp
integer(kind=i_single),parameter   :: NPAR_ARTH=16,NPAR2_ARTH=8

!------------------------------------------------------------
                                  
if (n > 0) arth(1)=first

if (n <= NPAR_ARTH) then

   do k=2,n
      arth(k)=arth(k-1)+increment
   enddo

else

   do k=2,NPAR2_ARTH
      arth(k)=arth(k-1)+increment
   enddo

   temp=increment*NPAR2_ARTH
   k=NPAR2_ARTH

   do
      if (k >= n) exit
         k2=k+k
         arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
         temp=temp+temp
         k=k2
   enddo

endif

END FUNCTION arth

!===================================================================================================

SUBROUTINE swap(a,b)
!-------------------------------------------------------
!
! Description:
!
!   Swap the content of a and b
!
! Dependencies:
!
!   None
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1366-1367
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use def_kind

implicit none

! arrays to be swapped
complex(kind=r_single),dimension(:),intent(inout):: a,b
! local variable
complex(kind=r_single),dimension(size(a))        :: dum

!-------------------------------------------------------
                                  
dum=a
a=b
b=dum

END SUBROUTINE swap

!===================================================================================================

FUNCTION zroots_unity(n,nn)
!-------------------------------------------------------
!
! Description:
!
!   Return nn powers of the n-th root of unity
!
! Dependencies:
!
!   None
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1379
!
! Author: W. Imperatori
!
! Modified: January 2009 (v1.3)
!

use constants; use def_kind

implicit none

! input variables
integer(kind=i_single),intent(in)   :: n,nn
! function
complex(kind=r_single),dimension(nn):: zroots_unity
! local variables
integer(kind=i_single)              :: k
real(kind=r_single)                 :: theta

!-------------------------------------------------------
                               
zroots_unity(1)=1.0
theta=pi_double/n
k=1

do
   if (k >= nn) exit
   zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),kind=r_single)
   zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*zroots_unity(2:min(k,nn-k))
   k=2*k
enddo

END FUNCTION zroots_unity
