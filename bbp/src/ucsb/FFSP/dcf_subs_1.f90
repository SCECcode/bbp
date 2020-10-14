!
!================================================================
!
subroutine dcf_logmean(fc1,fc2,logmean)
 use time_freq
 implicit none
 integer::  i,j,i1,i2
 real,parameter :: pi=3.14159265,pi2=2.0*pi,pi4=pi2*pi2
 real:: fc1,fc2,energy
 real:: moment_rate(nphf),sum
 real:: logmean(lnpt)
 if(io_dva.ne.1.and.io_dva.ne.2.and.io_dva.ne.3)then
    write(*,*)"io_dva has only three options 1 2 3 "
    write(*,*)"io_dva= ",io_dva
    stop
 endif
 if(nbb.lt.1.or.nee.gt.(lnpt-1))then
    write(*,*)"nb and ne exceed the possible range"
    write(*,*)"nb ne = ",nbb,nee
    stop
 endif
 do i=1,nphf
   freq(i)=(i-1+1.0e-5)*df
   moment_rate(i)=1.0/((1.+(freq(i)/fc1)**4.0)**0.25)/((1.+(freq(i)/fc2)**4.0)**0.25)
   if(io_dva.eq.2)then
      moment_rate(i)=pi2*freq(i)*moment_rate(i)
   elseif(io_dva.eq.3)then
      moment_rate(i)=pi4*freq(i)*freq(i)*moment_rate(i)
   endif
 enddo
 do i=nbb,nee
    i1=2**(i-1)
    i2=2**i-1
    sum=0.0
    do j=i1,i2
       sum=sum+log(moment_rate(j))
    enddo
    logmean(i)=sum/(i2-i1+1)
 enddo
end subroutine dcf_logmean
!
!================================================================
! td is the shortest duration for percent*m0
subroutine stf_duration(percent,td,tb,te)
 use time_freq
 implicit none
 integer:: i,j,ntd,n0,n1,n2,ib,ie,n
 integer:: i_t0,i_t1,i_t2
 real:: percent,td,tb,te,left
 real:: total,target,epslon,x
 real:: svf(ntime),m0(ntime)
 svf=0.0
 m0=0.0
 total=0.0
 call sum_point_svf(svf)
 total=sum(svf)
 epslon=(1-1.0e-6)*total
 target=total*percent
 left=total-target
 n0=1
 n1=1
 n2=1
 i_t0=0
 i_t1=0
 i_t2=0
 do i=2,ntime
   m0(i)=m0(i-1)+svf(i)
   if(m0(i).ge.left.and.i_t0.eq.0)then
     n0=i
     i_t0=1
   endif
   if(m0(i).ge.target.and.i_t1.eq.0)then
     n1=i
     i_t1=1
   endif
   if(m0(i).ge.epslon.and.i_t2.eq.0)then
     n2=i
     i_t2=1
   endif
 end do
 write(*,*)"n0 n1 n2",n0,n1,n2,dt,total
 ntd=n1
 ib=1
 ie=n1
 do i=1,n0
   do j=n1,n2
     n=j-i+1
     x=m0(j)-m0(i)
     if((n.lt.ntd).and.x.ge.target)then
       ntd=n
       ib=i
       ie=j
     endif
   end do
 end do
 td=ntd*dt
 tb=(ib-1)*dt
 te=ie*dt
 return
end subroutine stf_duration
!
!
!================================================================
!
subroutine stf_synth_logmean(rmt,logmean)
 use time_freq
 implicit none
 integer:: i,j,i1,i2
 real,parameter :: pi=3.14159265,pi2=2.0*pi,pi4=pi2*pi2
 real:: rmt,dtmt,fc1,fc2,energy
 real:: svf(ntime),moment_rate(nphf)
 real:: density,vs,coef,tim,dtm,sum
 real:: logmean(lnpt)
 svf=0.0
 call sum_point_svf(svf)
 call realft(svf,ntime,-1)
 dtmt=dt/rmt
 write(*,*)"dtmt=",dtmt
 do i=1,ntime
    svf(i)=svf(i)*dtmt
 enddo

 do i=1,nphf-1
    moment_rate(i)=sqrt(svf(2*(i-1)+1)**2+svf(2*(i-1)+2)**2)
    freq(i)=(i-1+1.0e-5)*df
    if(io_dva.eq.2)then ! velocity
       moment_rate(i)=pi2*freq(i)*moment_rate(i)
    elseif(io_dva.eq.3)then ! acceleration
       moment_rate(i)=pi4*freq(i)*freq(i)*moment_rate(i)
    endif
 enddo
 moment_rate(nphf)=svf(2)
 do i=nbb,nee
    i1=2**(i-1)
    i2=2**i-1
    sum=0.0
    do j=i1,i2
       sum=sum+log(moment_rate(j))
    enddo

    logmean(i)=sum/(i2-i1+1)
 enddo
end subroutine  stf_synth_logmean
!
!================================================================
!
subroutine misfit_ratio(rmt,fc1,fc2,misfit,misfit_2)
  use time_freq
  implicit none
  integer:: i
  real:: logmean_m(lnpt),logmean_s(lnpt)
  real:: misfit,rmt,fc1,fc2,misfit_2
  misfit=0.0
  misfit_2=0.0
  call stf_synth_logmean(rmt,logmean_s)
  call dcf_logmean(fc1,fc2,logmean_m)
  do i=nbb,nee
!     write(*,*)i, logmean_s(i), logmean_m(i)
     misfit=misfit+logmean_s(i)-logmean_m(i)
     misfit_2=misfit_2+(logmean_s(i)-logmean_m(i))**2.0
  enddo
  misfit=exp(misfit/(nee-nbb+1.))
  misfit_2=sqrt(misfit_2/(nee-nbb+1.))
end subroutine misfit_ratio
!
!================================================================
!
subroutine misfit_lsq(rmt,fc1,fc2,misfit)
  use time_freq
  implicit none
  integer:: i
  real:: logmean_m(lnpt),logmean_s(lnpt)
  real:: misfit,rmt,fc1,fc2
  misfit=0.0
  call stf_synth_logmean(rmt,logmean_s)
  call dcf_logmean(fc1,fc2,logmean_m)
  do i=nbb,nee
!     write(*,*)i, logmean_s(i), logmean_m(i)
     misfit=misfit+(logmean_s(i)-logmean_m(i))**2.0
  enddo
  misfit=sqrt(misfit/(nee-nbb+1.))
end subroutine misfit_lsq
!
!================================================================
!
subroutine stf_synth_output(rmt,fc1,fc2)
 use time_freq
 implicit none
 integer:: i,j,i1,i2
 real:: rmt,fc1,fc2,energy
 real:: svf(ntime),moment_rate(nphf),dcf(nphf)
 real:: density,vs,coef,pi,tim,dtmt,sum
 real:: logmean_m(lnpt),logmean_s(lnpt),f2(lnpt)
!svf=0.0
 call sum_point_svf(svf)
 open(19,file='calsvf_tim.dat',status='replace')
 write(19,*) ntime
 do i=1,ntime
    tim=(i-1)*dt
    write(19,*) tim,svf(i)/rmt
 enddo
 close(19)

 call realft(svf,ntime,-1)
 dtmt=dt/rmt
 do i=1,ntime
    svf(i)=svf(i)*dtmt
 enddo
 do i=1,nphf-1
    moment_rate(i)=sqrt(svf(2*(i-1)+1)**2+svf(2*(i-1)+2)**2)
 enddo
 moment_rate(nphf)=svf(nphf-1)
 do i=1,lnpt-1
    i1=2**(i-1)
    i2=2**i-1
    sum=0.0
    do j=i1,i2
       sum=sum+log(moment_rate(j))
    enddo
    logmean_s(i)=sum/(i2-i1+1)
    f2(i)=0.5*(freq(i1)+freq(i2))
 enddo

 open(19,file='calsvf.dat')
 write(19,*) nphf
 do i=1,nphf
    dcf(i)=1./((1.+(freq(i)/fc1)**4.0)**0.25)/((1.+(freq(i)/fc2)**4.0)**0.25)
    write(19,*) freq(i),moment_rate(i),dcf(i)
 enddo
 close(19)
 do i=1,lnpt-1
    i1=2**(i-1)
    i2=2**i-1
    sum=0.0
    do j=i1,i2
       sum=sum+log(dcf(j))
    enddo
    logmean_m(i)=sum/(i2-i1+1)
 enddo
 open(19,file='logsvf.dat')
 write(19,*)lnpt-1
 do i=1,lnpt-1
    write(19,*)f2(i),logmean_s(i),logmean_m(i)
 enddo
 close(19)
end subroutine stf_synth_output
