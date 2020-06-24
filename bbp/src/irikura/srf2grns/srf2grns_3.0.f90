! Last updated 2019.5 by A. Iwaki

!  This program reads BBP SRF file
!  and generates elem_param.dat for "greenscale"
!  in SGF method used in Irikura Recipe Method 2.

!! VER3.0 !!  !  Works for single- and multi- segments.
program srf2grns

 use params

 implicit none
 real::sdr,sdrsa,sdrback,dback  !stress drop and background slip
 real::a,xlen1,xlen2
 integer::i,iseg,ni,ne,nxf,nyf,npp,j,ii,k,nseg
 character::srffile0*120,outfile*120,fltpar*120,stat*4,dum*1,sdropout*120

  read(5,*) iseg  ! segment number
  read(5,'(a120)') srffile0 ! single-planar fault
  read(5,'(a120)') srffile ! segment
  read(5,'(a120)') fltpar  ! segment.midpoint.4
  read(5,'(a120)') sdropout  ! stress_drop.out
  read(5,'(a120)') outfile ! segment


  ! read from BBP .SRF file  ! multi-segment OK 
  call readsrf

  ! Single Planar Fault
  ! nxf, nyf = # of subfaults along strike & dip for the total fault
  open(10,file=srffile0,form='formatted',status='old')
  read(10,*) a
  read(10,'(a1)') dum
  read(10,*) a,a,nxf,nyf,a,a
  close(10)

  ! Segment Information
  open(11,file=fltpar,form='formatted',status='old')
  read(11,*)
  read(11,*) nseg
  do i=1,iseg
     read(11,*) stat,a,a,a,xlen1,xlen2,a,a,a,a,a,a
  enddo
  close(11)
  ni=int(xlen1*1000./dx) ; ne=int(xlen2*1000./dx)
  npp=(ne-ni)*nyf
  if(np.ne.npp) then
     write(6,*) 'ERROR np != npp ',np,npp,xlen1,xlen2,nyf ; stop
  endif

! Read from "stress_drop.out" 
  sdrop=0.0e0 ; ii=0
  open(12,file=sdropout,form='formatted',status='old')
  write(6,*) 'read from stress_drop.out',np,ni,ne,nyf
  do j=1,nyf
     do i=1,nxf
        read(12,*) k,k,a
        if((i>ni).and.(i<=ne)) then  ! reads only within this segment
           ii=ii+1
           sdrop(ii)=a
        endif
     enddo
  enddo
  write(6,*) ' ii = ',ii
  close(12)

! resampling to 2km interval
  call resample

! elem_param.dat for greenscale
  open(20,file=outfile,form='formatted')
  do i=1,np2
     write(20,22) i,lat2(i),lon2(i),dep2(i),slip2(i),sdrop2(i),tinit2(i)
  enddo
  close(20)
22 format(i4,2f12.6,f10.2,2f12.6,f10.4)

end program srf2grns

