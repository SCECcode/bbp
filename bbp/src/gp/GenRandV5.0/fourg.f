      subroutine fourg (data,n,isign,work)
c     cooley-tukey fast fourier transform in usasi basic fortran.
c     one-dimensional transform of complex data, arbitrary number of
c     points.  n points can be transformed in time proportional to
c     n*log(n) (for n non-prime), whereas other methods take n**2 time.
c     furthermore, because fewer arithmetic operations are performed,
c     less error is built up.  the transform done is--
c     dimension data(n),transform(n),work(n)
c     complex data,transform,work
c     transform(k) = sum(data(j)*exp(isign*2*pi*i*(j-1)*(k-1)/n)),
c     summed from j = 1 to n for all k from 1 to n.  the transform
c     values are returned to data, replacing the input.  n may be any
c     positive number, but it should be non-prime for speed.  isign =
c     +1 or -1.  a -1 transform followed by a +1 one (or vice versa)
c     returns n times the original data.  work is a one-dimensional
c     complex array of length n used for working storage.
c     running time is proportional to n * (sum of the prime factors of
c     n).  for example, n = 1960, time is t0 * 1960 * (2+2+2+5+7+7).
c     naive methods directly implementing the summation run in time
c     proportional to n**2.  an upper bound for the rms relative error
c     is 3 * 2**(-b) * sum(f**1.5), where b is the number of bits in
c     the floating point fraction and the sum is over the prime
c     factors of n.  written by norman brenner, mit lincoln laboratory,
c     august 1968.  see--ieee transactions on audio and e
c     (june 1967), special issue on the fast fourier transform.
      dimension data(1), work(1), ifact(32)
      twopi=6.283185307*float(isign)
c     factor n into its prime factors, nfact in number.  for example,
c     for n = 1960, nfact = 6 and ifact(if) = 2, 2, 2, 5, 7 and 7.
      if=0
      npart=n
      do 50 id=1,n,2
      idiv=id
      if (id-1) 10,10,20
 10   idiv=2
 20   iquot=npart/idiv
      if (npart-idiv*iquot) 40,30,40
 30   if=if+1
      ifact(if)=idiv
      npart=iquot

      go to 20
 40   if (iquot-idiv) 60,60,50
 50   continue
 60   if (npart-1) 80,80,70
 70   if=if+1
      ifact(if)=npart
 80   nfact=if
c     shuffle the data array by reversing the digits of the index.
c     replace data(i) by data(irev) for all i from 1 to n.  irev-1 is
c     the integer whose digit representation in the multi-radix
c     notation of factors ifact(if) is the reverse of the r
c     of i-1.  for example, if all ifact(if) = 2, then for i-1 = 11001,
c     irev-1 = 10011.  a work array of length n is needed.
      ip0=2
      ip3=ip0*n
      iwork=1
      i3rev=1
      do 110 i3=1,ip3,ip0
      work(iwork)=data(i3rev)
      work(iwork+1)=data(i3rev+1)
      ip2=ip3
      do 100 if=1,nfact
      ip1=ip2/ifact(if)
      i3rev=i3rev+ip1
      if (i3rev-ip2) 110,110,90
 90   i3rev=i3rev-ip2
 100  ip2=ip1
 110  iwork=iwork+ip0
      iwork=1
      do 120 i3=1,ip3,ip0
      data(i3)=work(iwork)
      data(i3+1)=work(iwork+1)
 120  iwork=iwork+ip0
c     phase-shifted fourier transform of length ifact(if).
c     iprod=ip1/ip0
c     irem=n/(ifact(if)*iprod)
c     dimension data(iprod,ifact(if),irem),work(ifact(if))
c     complex data,work
c     data(i1,j2,i3) = sum(data(i1,i2,i3) * w**(i2-1)), summed over
c     i2 = 1 to ifact(if) for all i1 from 1 to iprod, j2 from 1 to

c     ifact(if) and i3 from 1 to irem.
c     w = exp(isign*2*pi*i*(i1-1+iprod*(j2-1))/(iprod*ifact(if))).
      if=0
      ip1=ip0
 130  if (ip1-ip3) 140,240,240
 140  if=if+1
      ifcur=ifact(if)
      ip2=ip1*ifcur
      theta=twopi/float(ifcur)
      sinth=sin(theta/2.)
      rootr=-2.*sinth*sinth
c     cos(theta)-1, for accuracy
      rooti=sin(theta)
      theta=twopi/float(ip2/ip0)
      sinth=sin(theta/2.)
      wstpr=-2.*sinth*sinth
      wstpi=sin(theta)
      wminr=1.
      wmini=0.
      do 230 i1=1,ip1,ip0
      if (ifcur-2) 150,150,170
 150  do 160 i3=i1,ip3,ip2
      j0=i3
      j1=i3+ip1
      tempr=wminr*data(j1)-wmini*data(j1+1)
      tempi=wminr*data(j1+1)+wmini*data(j1)
      data(j1)=data(j0)-tempr
      data(j1+1)=data(j0+1)-tempi
      data(j0)=data(j0)+tempr
 160  data(j0+1)=data(j0+1)+tempi
      go to 220
 170  iwmax=ip0*ifcur
      do 210 i3=i1,ip3,ip2
      i2max=i3+ip2-ip1
      wr=wminr
      wi=wmini
      do 200 iwork=1,iwmax,ip0
      i2=i2max
      sumr=data(i2)
      sumi=data(i2+1)

 180  i2=i2-ip1
      tempr=sumr
      sumr=wr*sumr-wi*sumi+data(i2)
      sumi=wr*sumi+wi*tempr+data(i2+1)
      if (i2-i3) 190,190,180
 190  work(iwork)=sumr
      work(iwork+1)=sumi
      tempr=wr
      wr=wr*rootr-wi*rooti+wr
 200  wi=tempr*rooti+wi*rootr+wi
      iwork=1
      do 210 i2=i3,i2max,ip1
      data(i2)=work(iwork)
      data(i2+1)=work(iwork+1)
 210  iwork=iwork+ip0
 220  tempr=wminr
      wminr=wminr*wstpr-wmini*wstpi+wminr
 230  wmini=tempr*wstpi+wmini*wstpr+wmini
      ip1=ip2
      go to 130
 240  return
      end
