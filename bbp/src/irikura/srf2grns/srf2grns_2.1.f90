! Last updated 2018.7.24 by A. Iwaki

!  This program reads BBP SRF file
!  and generates elem_param.dat for "greenscale"
!  in SGF method used in Irikura Recipe Method 2.
!  Works for single segment only.

program srf2grns

 use params

 implicit none
 real::sdr,sdrsa,sdrback,dback  !stress drop and background slip
 integer::i
 character::outfile*120

  read(5,'(a120)') srffile
  read(5,'(a120)') outfile

! read from BBP .SRF file  ! single segment
  call readsrf

! stress drop by Irikura recipe
  !call sdrop(len,wid,dens,vs,sdrsa,sdrback,dback)
  call sdrop(len,wid,dens,vs,musum,sdrsa,sdrback,dback)

! resampling to 2km interval
  call resample

! elem_param.dat for greenscale
  open(20,file=outfile,form='formatted')
  do i=1,np2
!     write(6,*) lon2(i),lat2(i),dep2(i),tinit2(i),slip2(i)
     if((slip2(i)>0.9*dback).and.(slip2(i)<1.1*dback)) then
        sdr=sdrback
     else
        sdr=sdrsa
     endif
     write(20,22) i,lat2(i),lon2(i),dep2(i),slip2(i),sdr,tinit2(i)
  enddo
  close(20)
22 format(i4,2f12.6,f10.2,2f12.6,f10.4)

end program srf2grns

