c
        subroutine zpass(z,n,dt,ift,flo,fhi,nphas,nbut,iff,iq)
c
        dimension z(1)
        dimension x(30000)

        do 1000 i=1,n
 1000   x(i)=z(i)

        tupi=2.0*3.141592654
        iflag=0
        znyq=1.0/(2.0*dt)

        go to (20,30,40), ift+1
  20    a=1.0
        f1=flo
        w0=tupi*f1
        go to 50

  30    a=-1.0
        f1=fhi
        w0=tupi*(znyq-f1)
        go to 50

  40    a=1.0
        f1=flo
        f2=fhi
        w0=tupi*f2
        w1=tupi*(znyq-f1)
        iflag=1

  50    continue

c        iq=(1+(-1)**nbut)/2

	iff = 1
	do 999 j=1,nbut
	   iff = -1*iff
 999    continue
	
        iq=(1+iff)/2

        n2b=nbut/2
        do 70 i=1,n2b
        ct=cs2(w0,nbut,tupi,iq,i)

        ict=i
        call butter (x,a,ct,w0,ict,dt,n)
  70    continue

        if (iq.eq.0) call bttr2 (x,a,w0,dt,n)
        if(nphas.eq.1) go to 88
        call revers (x,n)
        do 80 i=1,n2b
        ct=cs2(w0,nbut,tupi,iq,i)

        ict=i
        call butter (x,a,ct,w0,ict,dt,n)
  80    continue
        if (iq.eq.0) call bttr2 (x,a,w0,dt,n)
        call revers (x,n)
  88    if (iflag.eq.0) go to 130
        a=-1.0
        w0=w1
        iflag=0
        go to 50
 130    continue
        do 131 i=1,n
 131    z(i)=x(i)
        return
        end
c
        function cs2(w0,nbut,tupi,iq,i)
c
        arg = (float(2*i-iq)*tupi/float(4*nbut))
        cs2=2.0*w0*cos(arg)
        return
        end
c
        subroutine cs2x(ct,w0,nbut,tupi,iq,i)
c
        arg = (float(2*i-iq)*tupi/float(4*nbut))
        ct=2.0*w0*cos(arg)
        return
        end
c
        subroutine butter (x,a,b,w,ict,dt,n)
c
        dimension x(1)
        
        za=2.0*a
        call set (w,dt,n,wp,ak)
        wp2=wp*wp
        b0=wp2
        b1=1.0+b*ak+wp2
        b2=za*(wp2-1.0)
        b3=1.0-b*ak+wp2
        if (a.gt.0.0) go to 20
        yt1=0.0
        yt2=0.0
        if (ict.gt.1) go to 10
        xt1=x(1)
        xt2=x(1)
        go to 30
  10    xt1=0.0
        xt2=0.0
        go to 30
  20    yt1=x(1)
        yt2=x(1)
        xt1=x(1)
        xt2=x(1)
  30    do 40 i=1,n
        yt3=(b0*(x(i)+za*xt2+xt1)-b2*yt2-b3*yt1)/b1
        xt1=xt2
        xt2=x(i)
        x(i)=yt3
        yt1=yt2
        yt2=yt3
  40    continue
        return
        end
c
        subroutine bttr2 (x,a,w0,dt,n)
c
        dimension x(1)

        call set (w0,dt,n,wp,ak)
        b0=wp+1.0
        b1=a*(wp-1.0)
        if (a.gt.0.0) go to 10
        yt=0.0
        xt=0.0
        go to 20
  10    yt=x(1)
        xt=x(1)
  20    do 30 i=1,n
        yt=(wp*(x(i)+a*xt)-b1*yt)/b0
        xt=x(i)
        x(i)=yt
  30    continue
        return
        end
c
        subroutine set (w,dt,n,wp,ak)
c
        arg = (w*dt/2.0)
        wp=tan(arg)
        ak=wp/w
        return
        end
c
        subroutine revers (x,n)
c
        dimension x(1)
        dimension y(30000)

        do 1 i=1,n
    1   y(i)=x(n-i+1)
        do 2 i=1,n
    2   x(i)=y(i)
        return
        end
