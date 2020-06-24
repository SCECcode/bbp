program dtrandom

 implicit none
 integer::imax,i,ii,seedsize,iseed
 integer,allocatable::seed(:)
 real,parameter::pmax=0.1  ! plus or minus 10 % of original t
 real::t,t0,p,rnd,rlat,rlon,dep,slip,sigma
 character::inpf*60,outf*60

  read(5,*) iseed
  read(5,*) imax
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  seed=iseed
  call random_seed(put=seed)

  open(11,file='rnd.tmp',form='formatted',status='unknown')
  do i=1,imax
     read(5,*) ii,rlat,rlon,dep,slip,sigma,t0
     call random_number(rnd)
     p=rnd*2*pmax-pmax
     t=t0+t0*p
     write(6,22) ii,rlat,rlon,dep,slip,sigma,t
     write(11,*) ii,rnd,p*100,t-t0
  enddo

22 format(i3,1x,2(f14.8,1x),f8.1,1x,2(f6.3,1x),f9.5)
end program dtrandom

