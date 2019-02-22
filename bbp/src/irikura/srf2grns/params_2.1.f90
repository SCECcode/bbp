 module params
   integer,save::np,nstk,ndip,np2
   real,parameter::dxx=2000.,dyy=2000.
   real,allocatable,save::lon(:),lat(:),dep(:),tinit(:),slip(:)
   real,allocatable,save::lon2(:),lat2(:),dep2(:),tinit2(:),slip2(:)
   real,save::len,wid,dx,dy,vs,dens,rake,musum
   character::srffile*120
 end module params
