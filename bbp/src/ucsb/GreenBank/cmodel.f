c        lupei's format. the source at the top of source layer
c
       subroutine cmodel(io_model,depth)
c
       implicit none
       include 'model.h'
       integer jo,nb,j,io_model,jj,j1,ncom_out
       real depth,c(nlay), s(nlay),den(nlay),th(nlay),qa(nlay),qb(nlay)
       real x,cc(nlay),ss(nlay),dd(nlay),tth(nlay),dep(nlay)
       real qqa(nlay),qqb(nlay),unrelaxed_mu,freq_ref
       common/model_new/jo,nb,c,s,den,th,qa,qb,ncom_out
       rewind(io_model)
       read(io_model,*) jo, freq_ref
       do j=1,jo
       read (io_model,*) cc(j),ss(j),dd(j),tth(j),qqa(j),qqb(j)
       enddo 
       rewind(io_model)
       x=0.0
       dep(1)=0.0
       do 84 j=1,jo
          x=x+tth(j)
          dep(j+1)=x
          c(j)=cc(j)
          s(j)=ss(j)
          den(j)=dd(j)
          th(j)=tth(j)
          qa(j)=qqa(j)
          qb(j)=qqb(j)
84     continue
       nb=2
       do 85 j=1,jo
         if(depth.gt.dep(j)) nb=j+1
85     continue

86    th(nb-1)=depth-dep(nb-1)
      do 87 j=nb,jo
      jj=jo+nb-j
      j1=jj+1
      th(j1)=th(jj)
      c(j1)=c(jj)
      s(j1)=s(jj)
      qa(j1)=qqa(jj)
      qb(j1)=qqb(jj)
87    den(j1)=den(jj)
      th(nb)=dep(nb)-depth
      c(nb)=c(nb-1)
      s(nb)=s(nb-1)
      den(nb)=den(nb-1) 
      qa(nb)=qqa(nb-1)
      qb(nb)=qqb(nb-1)
      if(th(nb).le.0.0) th(nb)=0.001
      jo=jo+1
        do jj=1,jo
          c(jj)=c(jj)*sqrt(unrelaxed_mu(freq_ref,1.0,qa(jj)))
          s(jj)=s(jj)*sqrt(unrelaxed_mu(freq_ref,1.0,qb(jj)))
        enddo

        nb=nb-1
        open(25,file='model.out')
        write(25,*)jo,nb,depth
        do j=1,jo
          write(25,*)c(j),s(j),den(j),th(j),qa(j),qb(j)
        enddo
        close(25)
       return
       end
c
       subroutine trav(dis,nx,tmin)
c      program travel
c============================================================
c calculate travel time for horizontal layered model
c it outputs both times for first arrival and direct arrival
c============================================================
      implicit none
      include 'model.h'
      real tmin(*),dis(*)
      real*8 th(50), v(50), t, t0, td, x
      complex*16 pd, p0, p
      integer i, ns, nx, lmax, jo, nb, ncom_out
      real c(nlay),s(nlay),den(nlay),thh(nlay),qa(nlay),qb(nlay)
      common/model_new/jo,nb,c,s,den,thh,qa,qb,ncom_out
      common/ccmodel/th,v,ns,lmax
      ns=nb
      do i=1,jo
        th(i)=thh(i)
        v(i)=1./c(i)**2
      enddo

      do i=1,nx
        x=dis(i)
        lmax = ns
        call find2(x,td,pd)
        t0 = td
        p0 = pd
        do lmax=ns+1,jo-1
          call find2(x,t,p)
          if (t.lt.t0) then
             t0=t
             p0=p
          endif
        enddo
c        write(*,'(f8.2,f8.2,f8.2,f8.4,f8.4)')x,t0,td,real(p0),real(pd)
        tmin(i)=t0
      enddo
      return
      end
c
      subroutine find2(x,t0,p0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     solve d(tau)/dp = 0 for p0, tau=p*x+eta*z
c     input:  x distance range
c     output: p0, t0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer i, ns, lmax
      real*8 aminm, dtdp0, t0, x, pm, th(50), v(50)
      complex*16 taup, dtdp, p0, p1, p2
      common/ccmodel/th,v,ns,lmax
      aminm = 1.e-6
      p1 = cmplx(0.,1.0e-20)    ! make p slightly above the real axis
      pm = 9999999.d0
      do i=1,lmax
        pm = dmin1(pm,v(lmax))
      enddo
      p2 = dsqrt(pm) + p1
      p0 = 0.5d0*(p1+p2)          ! the actual p0 lies betw. p1 and p2
      do while ( dble(p2-p1).gt.aminm )
         dtdp0 = dtdp(x,p0)
         if( dabs(dtdp0).lt.aminm )goto 100
         if( dtdp0.gt.0.d0 )then
            p1 = p0
         else
            p2 = p0
         endif
         p0 = 0.5d0*(p1+p2)
      enddo
 100  pm  = dmin1(pm,v(lmax+1))
      pm = dsqrt(pm)
      if (lmax.gt.ns.and.pm.lt.dble(p0)) p0 = pm
      t0 = taup(p0,x)
      return
      end

      function taup(p,x)
c define function tau(p) = p x + eta h
      implicit none
      integer i, ns, lmax
      real*8 x, th(50), v(50)
      complex*16 p, pp, taup
      common/ccmodel/th,v,ns,lmax
      taup = p*x
      pp = p*p
      do i = 1, ns
         taup=taup+cdsqrt(v(i)-pp)*th(i)
      enddo
      do i = ns+1,lmax
         taup=taup+cdsqrt(v(i)-pp)*th(i)*2.d0
      enddo
      return
      end

      function dtdp(x,p)
c define d(tau)/dp
      implicit none
      integer ns, j, lmax
      real*8 x, th(50), v(50)
      complex*16 p, pp, dtdp
      common/ccmodel/th,v,ns,lmax
      pp = p*p
      dtdp = 0.0d0
      do j = 1, ns
         dtdp=dtdp-th(j)/cdsqrt(v(j)-pp)
      enddo
      do j = ns+1, lmax
         dtdp=dtdp-2.d0*th(j)/cdsqrt(v(j)-pp)
      enddo
      dtdp = x + p*dtdp
      return
      end

C-----------------------------------------------------------------------
        real function unrelaxed_mu(fr,rmu_vm,qq)
        parameter (np=8,Qmin=5, Qmax=5000, ar_ha=1.0, pi2=6.283185)
        real relaxt(np),qwt1(np),qwt2(8),qwt(np)
        complex cxlnf
        data relaxt/0.172333E-02,  0.180701E-02,  0.538887E-02,
     &              0.199322E-01,  0.849833E-01,  0.409335E+00,
     &              0.205951E+01,  0.132629E+02/
        data qwt1/ 0.166958E-01,  0.381644E-01,  0.984666E-02,
     &            -0.136803E-01, -0.285125E-01, -0.537309E-01,
     &            -0.665035E-01, -0.133696E+00/
        data qwt2/ 0.898758E-01,  0.684635E-01,  0.967052E-01,
     &             0.120172E+00,  0.130728E+00,  0.138746E+00,
     &             0.140705E+00,  0.214647E+00/

        qr=alog(qq/Qmin)
        a1= 3.071
        a2= 1.433
        a3=-1.158
        a4= 0.415
        deltm=(a1+a2*qr*qq**a3)/(1.0+a4*qq)
        do i=1,np
          qwt(i)=deltm*(deltm*qwt1(i)+qwt2(i))
        enddo

        ww=fr*pi2
        cxlnf=(1.0, 0.0)
        do i=1,np
          cxlnf=cxlnf-qwt(i)/cmplx(1.,relaxt(i)*ww)
        enddo
        rtmp=cabs(cxlnf)
        unrelaxed_mu=rmu_vm/rtmp*(1.+real(cxlnf)/rtmp)/2.
c unrelaxed_mu=rmu_vm/rtmp
        return
        end
C-----------------------------------------------------------------------
        complex function cmu_freq(ww,cmu_unrelx,qq)
        parameter (np=8,Qmin=5, Qmax=5000, ar_ha=1.0, pi2=6.283185)
        real relaxt(np),qwt1(np),qwt2(8),qwt(np)
        complex cxlnf
        data relaxt/0.172333E-02,  0.180701E-02,  0.538887E-02,
     &              0.199322E-01,  0.849833E-01,  0.409335E+00,
     &              0.205951E+01,  0.132629E+02/
        data qwt1/ 0.166958E-01,  0.381644E-01,  0.984666E-02,
     &            -0.136803E-01, -0.285125E-01, -0.537309E-01,
     &            -0.665035E-01, -0.133696E+00/
        data qwt2/ 0.898758E-01,  0.684635E-01,  0.967052E-01,
     &             0.120172E+00,  0.130728E+00,  0.138746E+00,
     &             0.140705E+00,  0.214647E+00/

        qr=alog(qq/Qmin)
        a1= 3.071
        a2= 1.433
        a3=-1.158
        a4= 0.415
        deltm=(a1+a2*qr*qq**a3)/(1.0+a4*qq)
        do i=1,np
          qwt(i)=deltm*(deltm*qwt1(i)+qwt2(i))
        enddo

        cxlnf=(1.0, 0.0)
        do i=1,np
          cxlnf=cxlnf-qwt(i)/cmplx(1.,relaxt(i)*ww)
        enddo
        cmu_freq=cmu_unrelx*cxlnf
        return
        end
