        parameter (MAXnw=5000,MAXnl=5000) 
	dimension slip(MAXnw,MAXnl)
        dimension ty(2000)
	character*120 slipfile,planar_fault,fltpar,dirin,filein
        character*1 dum1
        character*5, stat

	data rperd/0.017453292/

        read(5,'(a100)')dirin
        read(5,'(a100)')planar_fault
        read(5,'(a100)')fltpar
        read(5,*)isegm
        read(5,*)dx,dy

        i1=index(dirin,' ')-1
        i2=index(planar_fault,' ')-1
        slipfile=dirin(1:i1)//'/'//'all_seg.'//planar_fault(1:i2)
        open(12,file=slipfile)
        version=2.0
        write (12,'(f3.1)')version
        write (12,'("PLANE ",i2)')isegm

        open(9,file=fltpar)
        read(9,'(a1)')dum1
        read(9,*)isegm
        do k=1,isegm
        read(9,*)stat,alon,alat,zz,xlen1,xlen2,ylen,x0,y0,
     +           stra,dipa,rakea
        ni=int(xlen1/dx)
        ne=int(xlen2/dx)
        xlen=(ne-ni)*dx
        nx=int(xlen/dx)
        ny=int(ylen/dy)
        write(12,*)alon,alat,nx,ny,xlen,ylen
        if(k.eq.1)then
         write(12,*)stra,dipa,zz,x0,y0
        else
         write(12,*)stra,dipa,zz,-999.9,-999.9
        endif
        enddo
        close(9)

        open(9,file=fltpar)
        read(9,'(a1)')dum1
        read(9,*)isegm
        do k=1,isegm
        read(9,*)stat,alon,alat,zz,xlen1,xlen2,ylen,x0,y0,
     +           stra,dipa,rakea
        i3=index(stat,' ')-1
        slipfile=dirin(1:i1)//'/'//stat(1:i3)//'.'//planar_fault(1:i2)
        print*,k,slipfile
        ni=int(xlen1/dx)
        ne=int(xlen2/dx)
        xlen=(ne-ni)*dx
        nx=int(xlen/dx)
        ny=int(ylen/dy)
        print*,'nx,ny = ',nx,ny
c        nt=nx*ny

        area=dx*dy*1.0e+10
c    fault element area is in cm*cm

        open(34,file=slipfile)
        read (34,*)dum
        read (34,'(a1)') dum1
        read (34,*) dum,dum,nxf,nyf,dum,dum
        read (34,*) dum,dum,dum,dum,dum
        read (34,'(a1)')dum1
        write (12,'(a,i8)') "POINTS",nx*ny
        do j=1,nyf
           do ii=1,nxf
        read (34,*)xx,yy,z,stra,dipa,areaf,stimf,dtf,
     +    vrf,densf
        write (12,200)xx,yy,z,stra,dipa,areaf,stimf,dtf,
     +    vrf,densf
        read (34,*) rakef,slip1f,npts1f,slip2f,n0,slip3f,n0
        write (12,500) rakef,slip1f,npts1f,slip2f,n0,slip3f,n0
        if(npts1f.ne.0) read (34,*) (ty(it),it=1,npts1f)
        if(npts1f.ne.0) write (12,400) (ty(it),it=1,npts1f)
          enddo
        enddo
        close(34)
        enddo
200     format(2(1x,f12.6),f9.2,2f8.2,1x,e13.5,1x,f8.4,3(1x,e12.4))
202     format(f10.4,3(1x,f10.2,i4))
300     format(f8.2,3(e13.5,i8))
400     format(6e13.5)
500     format(f8.2,1x,f8.2,1x,i5,2(1x,f6.2,i3))
        end
 
