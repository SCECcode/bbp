MODULE interfaces
 
INTERFACE 
   FUNCTION arth(first,increment,n)
      use def_kind
      integer(kind=i_single),intent(in)  :: first,increment,n
      integer(kind=i_single),dimension(n):: arth
   END FUNCTION arth
END INTERFACE 

INTERFACE
   FUNCTION convlv(data,respns)
      use def_kind
      implicit none
      real(kind=r_single),dimension(:),intent(inout):: data
      real(kind=r_single),dimension(:),intent(in)   :: respns
      real(kind=r_single),dimension(size(data))     :: convlv
   END FUNCTION convlv   
END INTERFACE

INTERFACE
   FUNCTION correl(data1,data2)
      use def_kind
      implicit none
      real(kind=r_single),dimension(:),intent(in) :: data1,data2
      real(kind=r_single),dimension(size(data1))  :: correl
   END FUNCTION correl
END INTERFACE

INTERFACE
   FUNCTION derive(fun,dt)
      use def_kind
      implicit none
      real(kind=r_single),dimension(:),intent(in):: fun
      real(kind=r_single),intent(in)             :: dt
      real(kind=r_single),dimension(size(fun)-1) :: derive
   END FUNCTION derive
END INTERFACE  

INTERFACE four1d
   SUBROUTINE four1d(data,isign)
      use def_kind
      implicit none
      complex(kind=r_single),dimension(:),intent(inout):: data
      integer(kind=i_single),intent(in)                :: isign
   END SUBROUTINE four1d
END INTERFACE four1d

INTERFACE fourrow
   SUBROUTINE fourrow(data,isign)
      use def_kind
      implicit none
      complex(kind=r_single),dimension(:,:),intent(inout):: data
      integer(kind=i_single),intent(in)                  :: isign
   END SUBROUTINE fourrow
END INTERFACE fourrow   

INTERFACE
   FUNCTION lag_correl(seq1,seq2)
      use def_kind
      implicit none
      real(kind=r_single),dimension(:),intent(in):: seq1,seq2
      real(kind=r_single)                        :: lag_correl
   END FUNCTION lag_correl
END INTERFACE

INTERFACE 
   FUNCTION normal_pdf(a,b,seed)
      use def_kind; use constants
      implicit none
      real(kind=r_scat),intent(in)     :: a,b
      integer(kind=i_single),intent(in):: seed
      real(kind=r_scat)                :: normal_pdf
   END FUNCTION normal_pdf
END INTERFACE

INTERFACE
   FUNCTION peaker(fun,dt)
      use def_kind
      real(kind=r_single),dimension(:),intent(in):: fun
      real(kind=r_single),intent(in)             :: dt
      real(kind=r_single)                        :: peaker
   END FUNCTION peaker
END INTERFACE

INTERFACE polint
   SUBROUTINE polint(xa,ya,x,y)
      use def_kind
      implicit none
      real(kind=r_single),dimension(:),intent(in):: xa,ya
      real(kind=r_single),intent(in)             :: x
      real(kind=r_single),intent(out)            :: y
   END SUBROUTINE polint
END INTERFACE polint

INTERFACE poly_interp
   SUBROUTINE poly_interp(x1a,x2a,ya,x1,x2,y)
      use def_kind
      implicit none
      real(kind=r_single),dimension(:),intent(in)  :: x1a,x2a  
      real(kind=r_single),dimension(:,:),intent(in):: ya
      real(kind=r_single),intent(in)               :: x1,x2
      real(kind=r_single),intent(out)              :: y
   END SUBROUTINE poly_interp
END INTERFACE poly_interp

INTERFACE 
   FUNCTION rand_pdf(seed,val_min,val_max,pdf)
      use def_kind
      implicit none
      real(kind=r_scat)                :: rand_pdf
      integer(kind=i_single),intent(in):: seed
      character(len=*),intent(in)      :: pdf 
      real(kind=r_scat),intent(in)     :: val_min,val_max
   END FUNCTION rand_pdf
END INTERFACE 

INTERFACE realft
   SUBROUTINE realft(data,isign,zdata)
      use def_kind
      implicit none
      real(kind=r_single),dimension(:),intent(inout)     :: data
      integer(kind=i_single),intent(in)                  :: isign
      complex(kind=r_single),dimension(:),optional,target:: zdata
   END SUBROUTINE realft
END INTERFACE realft

INTERFACE spline_interp
   SUBROUTINE spline_interp(in_seis,time_len,interp_npts,orig_npts,out_seis)
      use def_kind
      implicit none 
      integer(kind=i_single),intent(in)           :: interp_npts,orig_npts
      real(kind=r_single),intent(in)              :: time_len
      real(kind=r_single),dimension(:),intent(in) :: in_seis
      real(kind=r_single),dimension(:),intent(out):: out_seis
   END SUBROUTINE spline_interp
END INTERFACE spline_interp

INTERFACE swap
   SUBROUTINE swap(a,b)
      use def_kind
      complex(kind=r_single),dimension(:),intent(inout):: a,b
   END SUBROUTINE swap
END INTERFACE swap

INTERFACE
   FUNCTION trapz(fun,dt)
      use def_kind
      implicit none
      real(kind=r_single),dimension(:),intent(in):: fun
      real(kind=r_single),intent(in)             :: dt
      real(kind=r_single)                        :: trapz
   END FUNCTION trapz
END INTERFACE

INTERFACE 
   FUNCTION uniform_pdf(a,b,seed)
      use def_kind
      implicit none
      real(kind=r_scat),intent(in)     :: a,b
      integer(kind=i_single),intent(in):: seed
      real(kind=r_scat)                :: uniform_pdf
   END FUNCTION uniform_pdf   
END INTERFACE 

INTERFACE write_disk
   SUBROUTINE write_disk(station,type_flag,in_arr)
      use def_kind; use flags; use scattering, only: npts
      use source_receiver, only: stat_name
      use stf_data, only: npts_stf,total
      use waveform, only: lf_len,lf_npts
      implicit none
      integer(kind=i_single),intent(in)            :: station
      real(kind=r_single),dimension(:,:),intent(in):: in_arr
      character(len=*),intent(in)                  :: type_flag
   END SUBROUTINE write_disk
END INTERFACE write_disk

INTERFACE write_disk_mpi
   SUBROUTINE write_disk_mpi(station,type_flag,in_arr)
      use def_kind; use flags; use scattering, only: npts
      use source_receiver, only: stat_name
      use stf_data, only: npts_stf,total
      use waveform, only: lf_len,lf_npts
      implicit none
      integer(kind=i_single),intent(in)            :: station
      real(kind=r_single),dimension(:,:),intent(in):: in_arr
      character(len=*),intent(in)                  :: type_flag
   END SUBROUTINE write_disk_mpi
END INTERFACE write_disk_mpi

INTERFACE 
   FUNCTION zroots_unity(n,nn)
      use def_kind
      integer(kind=i_single),intent(in)   :: n,nn
      complex(kind=r_single),dimension(nn):: zroots_unity
   END FUNCTION zroots_unity    
END INTERFACE 

END MODULE interfaces
