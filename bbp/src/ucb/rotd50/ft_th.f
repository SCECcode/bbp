
      subroutine InterpFreq ( u, npts0, cu1, NN, nTotal, dt )
      real u(1)
      integer npts, m2, npts2, NN
      complex cu1(1)


      m2 = int( alog(float(npts0))/alog(2.) + 0.9999 )
      npts2 = 2**m2  

c     PAD TO POWER OF 2 in the time domain
      do i=npts0+1,npts2
        u(i) = 0.
      enddo  
      npts = npts2

c     FILL COMPLEX ARRAY
      do 10 i=1,npts
        cu1(i) = cmplx(u(i),0.0)
  10  continue

c     CALCULATE FORWARD FFT
      call cool ( -1., m2, cu1 )
      df = 1./(npts*dt)
c      do i=1,npts
c        write (60,'( f10.4, 2f12.4)') (i-1)*df, cu1(i)
c      enddo
      
c     Pad in the frequency domain
      iNyq = npts/2 + 1
      iNyq2 = NN * (npts/2) + 1
c      write (*,'( 2i10)') iNyq, iNyq2
      do i=iNyq+1,iNyq2
        cu1(i) = cmplx(0., 0.)
      enddo

c     Reset nyquist ot half its value (the other half will be added to the neg freq)
      cu1(iNyq) = cu1(iNyq)/2.

c     Load the negative frequencies
      nTotal = 2 * (iNyq2-1)
      do i=iNyq2+1,nTotal
        cu1(i) = conjg (cu1(nTotal+2-i))
      enddo

c     CALCULATE INVERSE FFT
      m3 = m2 + int( alog(float(NN))/alog(2.) + 0.5 )
      call cool ( 1., m3, cu1 )

c     Scale
      do i=1,nTotal
        u(i) = real(cu1(i)) /npts
      enddo
      
      return
      end

