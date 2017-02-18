
! ---------------------- BEGIN IMNMAX ----------------------
      subroutine imnmax(ia,nstrt,nstop,ninc,imin,imax) 
c
* Compute min, max for an integer array

c author: D. M. Boore
c last change: 7/17/97 - Written, based on mnmax.for
c
      integer ia(*)
      imax = ia( nstrt)
      imin=imax         
      do 10 i=nstrt,nstop,ninc
      if(ia(i)-imax) 15,15,20  
20    imax=ia(i)               
      go to 10                
15    if(ia(i)-imin) 25,10,10  
25    imin=ia(i)               
10    continue                
      return                  
      end                     
! ---------------------- END IMNMAX ----------------------
