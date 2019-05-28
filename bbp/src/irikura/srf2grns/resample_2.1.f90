subroutine resample
  use params

  implicit none
  integer::nskipx,nskipy,nst2,ndip2,i,ii,l,m
  integer,dimension(np)::iloc

  nskipx=dxx/dx ; nskipy=dyy/dy
!  nst2=int(len*1000./dxx) ; ndip2=int(wid*1000./dyy)
! 2018.7.30
!  nst2=int((nstk-1)/(dxx/dx))+1
!  ndip2=int((ndip-1)/(dyy/dy))+1
! 2018.8.6
  nst2=int((nstk-1)/nskipx)+1
  ndip2=int((ndip-1)/nskipy)+1

  np2=nst2*ndip2
  write(6,*) ' # of subfaults NP2',np2
  allocate(lon2(np2),lat2(np2),dep2(np2),tinit2(np2),slip2(np2))

  ! sort by depth,lat,lon
  call sort3(lon,lat,dep,np,iloc)
  i=0 ; ii=0
  do l=1,ndip
     do m=1,nstk
        i=i+1
        if((mod(l,nskipy)==1).and.(mod(m,nskipx)==1)) then
!           write(6,*) '   l, m, = ',l,m
           ii=ii+1
           lon2(ii)=lon(i) ; lat2(ii)=lat(i) ; dep2(ii)=dep(i)
           tinit2(ii)=tinit(iloc(i)) ; slip2(ii)=slip(iloc(i))
        endif
     enddo
  enddo
  if(np2.ne.ii) then
     write(6,*) 'arere ',np2,ii ; stop
  endif
  return
end subroutine resample

     
