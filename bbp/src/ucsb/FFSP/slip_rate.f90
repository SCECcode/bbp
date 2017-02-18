!======================================================================
! The subroutines in this file are for generating the slip rate 
! functions for each point source, and sum them together to get 
! spectra of larget event
!
! Written by Pengcheng Liu
! Copyright (c) 2005 by Pengcheng Liu
!
!======================================================================
subroutine sum_point_svf(svf)
 use time_freq
 use sp_sub_f
 implicit NONE
 integer:: i,j,j0,j1,j2,jmm
 real:: fcor_s
 real, dimension(ntime):: svf,svf_s
 real, dimension(2*ntime):: mr
 do i=1,2*ntime
    mr(i)=0.0
 enddo
 do i=1,ntime
    svf(i)=0.0
 enddo

 do i=1,nsum
   fcor_s=1./rstm(i)
   call slip_rate(id_sf_type,fcor_s,svf_s)
   j0=int(rptm(i)/dt+0.5)
   j1=j0+1
   j2=int(rstm(i)/dt)+j1
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
subroutine slip_rate(id_sf_type,f_cor,svf)
 use time_freq
 implicit NONE
 integer:: id_sf_type
 real:: f_cor,rise_time
 real, dimension(ntime):: svf
       
 if(id_sf_type == 1) then
   rise_time=1.00/f_cor
   call svf_time_1(ntime,dt,rise_time,svf)
   f_cor=3.*rise_time
 else if(id_sf_type == 2) then
   rise_time=1.0/f_cor
   call svf_time_2(ntime,dt,rise_time,svf)
   f_cor=rise_time
 else if(id_sf_type == 3) then
   rise_time=1.0/f_cor
   call svf_time_3(ntime,dt,rise_time,cft1,svf)
   f_cor=rise_time
 else if(id_sf_type == 4) then
   rise_time=1.00/f_cor
   call svf_time_4(ntime,dt,rise_time,cft1,svf)
   f_cor=rise_time
 else if(id_sf_type == 5) then
   rise_time=1.00/f_cor
   call svf_time_5(ntime,dt,rise_time,cft1,svf)
   f_cor=rise_time
 else if(id_sf_type == 6) then
   rise_time=1.0/f_cor
   call svf_time_6(ntime,dt,rise_time,cft1,svf)
   f_cor=rise_time
 else if(id_sf_type == 7) then
   rise_time=1.7/f_cor
   call svf_time_7(ntime,dt,rise_time,svf)
   f_cor=rise_time
 else if(id_sf_type == 8) then
   rise_time=1.0/f_cor
   call svf_time_8(ntime,dt,rise_time,cft1,svf)
   f_cor=rise_time
 endif

end subroutine slip_rate
!
!======================================================================
!
subroutine peak_slip_rate(pksr)
 use time_freq
 use sp_sub_f
 implicit NONE
 real, dimension(nsum):: pksr
 integer:: i, ni
 real:: f_cor,rt,pi2,tt1,tt2,tm,amp
 pi2=8.*atan(1.0)       
 if(id_sf_type == 1) then
   do i=1, nsum
     pksr(i)=slip(i)*exp(-1.0)*pi2/rstm(i)
   enddo
 else if(id_sf_type == 2) then
   do i=1, nsum
     tm=1.0/5.0
     pksr(i)=slip(i)*30.0*tm*(1.-tm)**4/rstm(i)
   enddo
 else if(id_sf_type == 3) then
   do i=1, nsum
     tt1=cft1*rstm(i)
     tt2=rstm(i)-tt1
     pksr(i)=slip(i)*pi2/(4.*tt1+0.5*pi2*tt2)
   enddo
 else if(id_sf_type == 4) then
   do i=1, nsum
     tt1=cft1*rstm(i)
     tt2=rstm(i)-tt1
     pksr(i)=slip(i)*pi2/(0.5*pi2*tt1+0.5*(4.*tt1+0.5*pi2*tt2))
   enddo
 else if(id_sf_type == 5) then
   do i=1, nsum
     tt1=cft1*rstm(i)
     tt2=rstm(i)-tt1
     pksr(i)=slip(i)/(4.*tt1/pi2+0.25*rstm(i))
   enddo
 else if(id_sf_type == 6) then
   do i=1, nsum
     tt1=cft1
     tt2=rstm(i)
     pksr(i)=slip(i)/(4.*tt1/pi2+0.5*tt2)
   enddo
 else
   do i=1, nsum
     pksr(i)=2.*slip(i)/rstm(i)
   enddo
 endif

end subroutine peak_slip_rate
!
!======================================================================
!
subroutine svf_time_1(nt,dt,rise_time,svf)
 implicit NONE
 real, parameter:: pi=3.1415926, pi2=2.*pi
 integer:: nt,np1,i
 real:: wc,ts,sum1,svi,dt,rise_time
 real,dimension(nt):: svf

 svf=0.0
 wc=pi2/rise_time
 np1=3*int(rise_time/dt+1.0)
 sum1=0.0
 do i=1,np1
   ts=(i-1)*dt
   svi=wc*wc*ts*exp(-wc*ts)
   svf(i)=svi
   sum1=sum1+svi
 enddo
 sum1=sum1*dt
 do i=1,np1
   svf(i)=svf(i)/sum1
 enddo

end subroutine svf_time_1
!
!======================================================================
!
subroutine svf_time_2(nt,dt,rise_time,svf)
 implicit none
 integer:: nd,nt,np,i
 real:: dt,rise_time,ts,cof
 real,dimension(nt):: svf

 nd=4
 np=int(rise_time/dt)+1
 cof=(nd+1.0)*(nd+2.0)
 do i=1,np
   ts=(i-1.)*dt/rise_time
   svf(i)=cof*ts*(1.-ts)**nd/rise_time
 enddo
 svf(np+1)=0.0
 return
end subroutine svf_time_2
!
!======================================================================
!
subroutine svf_time_3(nt,dt,rise_time,cft1,svf)
 implicit NONE
 real, parameter:: pi=3.1415926,pi2=2.*pi
 integer:: nt,n1,i
 real:: cft1,dt,rise_time,t21,t22,coef_norm,cf1,cf2,ts
 real,dimension(nt):: svf

 t21=cft1*rise_time
 t22=rise_time-t21
 coef_norm=pi2/(4.*t21+pi*t22)
 n1=int(rise_time/dt)+1
 do i=1,n1
   ts=(i-1.)*dt
   if(ts < t21) then
     svf(i)=coef_norm*sin(0.5*pi*ts/t21)
   else
     svf(i)=coef_norm*0.5*(1.+cos(pi*(ts-t21)/t22))
   endif
 enddo
 svf(n1+1)=0.0

end subroutine svf_time_3
!
!======================================================================
!
subroutine svf_time_4(nt,dt,rise_time,cft1,svf)
 implicit NONE
 real, parameter:: cs1=0.50,cs2=1.-cs1
 real, parameter:: pi=3.1415926, pi2=2.*pi
 integer:: nt,n1,i
 real:: cft1,dt,rise_time,tt1,tt2,coef_norm,cf1,cf2,ts
 real,dimension(nt):: svf
 tt1=cft1*rise_time
 tt2=rise_time-tt1
 coef_norm=pi2/(cs1*2.*pi*tt1+cs2*4.*tt1+cs2*pi*tt2)
 cf2=cs2*coef_norm
 cf1=cs1*coef_norm
 n1=int(rise_time/dt)+1
 do i=1,n1
   ts=(i-1)*dt
   if(ts < tt1) then
     svf(i)=cf2*sin(0.5*pi*ts/tt1)
   else
     svf(i)=cf2*0.5*(1+cos(pi*(ts-tt1)/tt2))
   endif
   if(ts < 2*tt1) then
     svf(i)=svf(i)+cf1*0.5*(1.-cos(pi*ts/tt1))
   endif
 enddo
 svf(n1+1)=0.0

end subroutine svf_time_4
!
!======================================================================
!
subroutine svf_time_5(nt,dt,rise_time,cft1,svf)
 implicit NONE
 real, parameter:: cs1=0.50,cs2=1.-cs1
 real, parameter:: pi=3.1415926,pi2=2.*pi,pid2=0.5*pi
 integer:: nt,n1,i
 real:: cft1,dt,rise_time,t11,t12,t21,t22,coef_norm,cf1,cf2,ts
 real,dimension(nt):: svf

 t11=cft1*rise_time
 t12=1.0*t11
 t21=1.0*t11
 t22=rise_time-t21
 coef_norm=pi2/(cs1*(4.*t11+pi*t12)+cs2*(4.*t21+pi*t22))
 cf1=cs1*coef_norm
 cf2=cs2*coef_norm
 n1=int(rise_time/dt)+1
 do i=1,n1
   ts=(i-1.)*dt
   if(ts < t11) then
     svf(i)=cf1*sin(pid2*ts/t11)
   elseif(ts < t11+t12) then
     svf(i)=cf1*0.5*(1.+cos(pi*(ts-t11)/t12))
   else
     svf(i)=0.0
   endif
   if(ts < t21) then
     svf(i)=svf(i)+cf2*sin(0.5*pi*ts/t21)
   else
     svf(i)=svf(i)+cf2*0.5*(1.+cos(pi*(ts-t21)/t22))
   endif
 enddo
 svf(n1+1)=0.0

end subroutine svf_time_5
!
!======================================================================
!
subroutine svf_time_6(nt,dt,rise_time,cft1,svf)
 implicit NONE
 real, parameter:: pi=3.1415926,pi2=2.*pi
 integer:: nt,n1,i
 real:: cft1,dt,rise_time,cs1,cs2,pid2,t21,t22,coef_norm,cf1,cf2,ts
 real,dimension(nt):: svf

 t21=cft1*rise_time
 t22=rise_time-t21
 coef_norm=pi2/(4.*t21+pi*t22)
 n1=int(rise_time/dt)+1
 do i=1,n1
   ts=(i-1.)*dt
   if(ts < t21) then
     svf(i)=coef_norm*sin(0.5*pi*ts/t21)
   else
     svf(i)=coef_norm*0.5*(1.+cos(pi*(ts-t21)/t22))
   endif
 enddo
 svf(n1+1)=0.0

end subroutine svf_time_6
!
!======================================================================
!
subroutine svf_time_7(nt,dt,rise_time,svf)
 implicit NONE
 integer, parameter:: Nv=4
 real, parameter::  ra=1.0
 integer:: nt,np,n2,i,j,k
 real:: dt,rise_time,aa,ai,tk,t0,t1,tt,cof
 real, dimension(nt):: svf

 np=int(rise_time/dt)+2
 do i=1,np
   svf(i)=0.0
 enddo
 aa=0.0
 do i=1,Nv
   aa=aa+ra**(i-1)
 enddo

 t1= rise_time/2**(Nv-1)
 do k=1,Nv
   ai=ra**(k-1)
   tk=t1*2**(k-1)
   t0=0.5*tk
   n2=int(tk/dt)+1
   cof=4./tk/tk/aa*ai
   do j=1,n2
     tt=(j-1.)*dt
     if(tt < t0) then
       svf(j)=svf(j)+cof*tt
     else
       svf(j)=svf(j)+cof*(tk-tt)
     endif
   enddo
 enddo

end subroutine svf_time_7
!
!======================================================================
!
subroutine svf_time_8(nt,dt,rise_time,cft1,svf)
 implicit NONE
 real, parameter:: pi=3.1415926, pi2=2.*pi
 integer:: nt,np0,np1,np2,i
 real:: cft1,ts,sum1,svi,dt,rise_time,Tp,Te,Tr,psv
 real,dimension(nt):: svf
 integer a_seed
 !Tp=cft1*rise_time
 call random_seed(a_seed)
 CALL random_number(Tp)
 Tp=0.1*Tp+0.15
 Te=0.8*rise_time
 Tr=rise_time
 if(Tr<Tp) Tp=cft1*Tr
 !svf=0.0
 np0=int(Tp/dt+1.0)
 np1=int(Te/dt+1.0)
 np2=int(Tr/dt+1.0)

 psv=sqrt(1.+100./(np0*dt))

 sum1=0.0
 do i=1,np0
    ts=(i-1)*dt
    svi=ts*psv/Tp*sin(0.5*pi/Tp*ts)
    svf(i)=svi
    sum1=sum1+svi
 enddo
 do i=np0+1,np1
    ts=(i-1)*dt
    svi=sqrt(1.+100./ts)
    svf(i)=svi
    sum1=sum1+svi
 enddo
 do i=np1+1,np2
    ts=(i-1)*dt
    svi=sqrt(1.+100./ts)*sin((np2-i)*dt*pi*0.5/(Tr-Te))
    svf(i)=svi
    sum1=sum1+svi
 enddo
 sum1=sum1*dt
 do i=1,np2
   svf(i)=svf(i)/sum1
 enddo
end subroutine svf_time_8
