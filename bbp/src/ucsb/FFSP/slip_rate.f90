!======================================================================
! The subroutines in this file are for generating the slip rate
! functions for each point source, and sum them together to get
! spectra of larget event
!
! Written by Pengcheng Liu
! Copyright (c) 2005 by Pengcheng Liu
!
!  Simply using the modified yoffe function as default
!         Chen Ji, 2020
!======================================================================
subroutine sum_point_svf(svf)
 use time_freq
 use sp_sub_f
 implicit NONE
 integer:: i,j,j0,j1,j2,jmm
 real:: rise_time
 real, dimension(ntime):: svf,svf_s
 real, dimension(2*ntime):: mr
  svf=0.0
  mr=0.0
 do i=1,nsum
   rise_time=rstm(i)+pktm(i)
   cft1=pktm(i)/rise_time
   call svf_yoffe(ntime,dt,rise_time,cft1,svf_s)
   j0=int(rptm(i)/dt+0.5)
   j1=j0+1
   j2=int(rise_time/dt)+j1
   do j=j1,j2
     mr(j)=mr(j)+slip(i)*svf_s(j-j0)
   enddo
 enddo
 do i=1,ntime
    svf(i)=mr(i)
 enddo
end subroutine sum_point_svf
!
!======================================================================
!
subroutine peak_slip_rate(pksr)
 use time_freq
 use sp_sub_f
 implicit NONE
 real, dimension(nsum):: pksr
 integer:: i
!
! just for modified yoffe function
!         Chen Ji, 2020
 do i=1,nsum
   pksr(i)=slip(i)/sqrt(pktm(i))/sqrt(rstm(i))
 end do

end subroutine peak_slip_rate
!
!======================================================================
! modified normalized yoffe function. Chen Ji, 2020
!
subroutine svf_yoffe(nt,dt,rise_time,cft1,svf)
 implicit NONE
 real, parameter:: pi=3.14159265
 integer:: nt,npt_yoffe,nsin,i,j,nall
 real:: cft1,dt,rise_time
 real:: ty,tsin,sn,t,tsin_min
 real,dimension(nt):: svf,yoffe,hsin

 svf=0.0
 if(cft1.ge.1.0)then
    write(*,*)"unphysical rise-time function"
    stop
 endif
 tsin_min=2.0*dt
 tsin=cft1*rise_time
 ty=rise_time-Tsin
 npt_yoffe=int(ty/dt+1)+1
 nsin=int(tsin/dt+0.1)+1
 nall=npt_yoffe+nsin
 do i=1,npt_yoffe
    t=i*dt
    if(t.le.ty)then
       yoffe(i)=sqrt(ty-t)/sqrt(t)
    else
       yoffe(i)=0.0
    endif
 enddo
 if(tsin.ge.tsin_min)then
   do i=1,nsin
     t=(i-1)*dt
     if(t.le.tsin)then
       hsin(i)=sin(t*pi/tsin)
     else
       hsin(i)=0.0
     endif
   enddo
   do i=1,npt_yoffe
     do j=1,nsin
       svf(i+j-1)=svf(i+j-1)+yoffe(i)*hsin(j)
     enddo
   enddo
 else
   do i=1,npt_yoffe
     svf(i)=yoffe(i)
   enddo
 endif
 sn=0.0
 do i=1,nall
    sn=sn+svf(i)
 enddo
 svf=svf/(dt*sn)
end subroutine svf_yoffe
!
!======================================================================
! modified normalized kostrov function at r/v=10
! t0 and t02
subroutine svf_kostrov(nt,dt,rise_time,cft1,svf)
 implicit NONE
 real, parameter:: pi=3.14159265,t0=10.0,t02=t0*t0
 integer:: nt,npt_kostrov,nsin,i,j
 real:: cft1,dt,rise_time
 real:: tk,tsin,sn,t
 real,dimension(nt):: svf,yk,hsin

 svf=0.0
 if(cft1.gt.1.0)then
    write(*,*)"unphysical rise-time function"
    stop
 endif

 tsin=cft1*rise_time
 tk=rise_time-Tsin
 npt_kostrov=int(tk/dt+1.0)+1
 nsin=int(tsin/dt+1.0)+1

 yk(1)=0.0
 do i=2,npt_kostrov
    t=(i-1)*dt
    if(t.le.tk)then
       yk(i)=(t+t0)/sqrt((t+t0)**2.0 -t02)
    else
       yk(i)=0.0
    endif
 enddo

 do i=1,nsin
    t=(i-1)*dt
    if(t.le.tsin)then
       hsin(i)=sin((i-1)*dt*pi/tsin)
    else
       hsin(i)=0.0
    endif
 enddo
 do i=1,npt_kostrov
    do j=1,nsin
        svf(i+j-1)=svf(i+j-1)+yk(i)*hsin(j)
    enddo
 enddo
 do i=1,npt_kostrov+nsin-1
    svf(i)=svf(i+1)
 enddo
 sn=0.0
 do i=1,npt_kostrov+nsin-1
    sn=sn+svf(i)
 enddo
 svf=svf/(dt*sn)
end subroutine svf_kostrov
