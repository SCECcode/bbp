c RWG 20141212
c    m (= nx)	number of points in the model in the horizontal direction
c    n (= nz)	number of points in the model in the vertical discretion
c    is		grid index of the source location in the horizontal direction
c    js		grid index of the source location in the vertical direction
c    h		grid step size (in km?)
c    ns		radius to compute tt analytically around source, up to 5 seems good
c    t		2D array for travel times
c
c Changed to pass arrays slwns, ttime, ti and jm from calling routine, so memory can
c be allocated/freed outside of this routine
c
c..............................................................................
      subroutine wfront2d(m,n,is,js,h,ns,ttime,slwns,ntot,ti,jm)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension slwns(ntot),ttime(ntot),ti(ntot),jm(ntot)

c      print *,m,n,is,js
c      print *,h,ns,ntot,ttime(1),slwns(1)

      eps = 1e-5
c     
      ttime(IS + m*(JS-1)) = 0.0
      IF(IS.LT.M) ttime(IS+1 + m*(JS-1)) = ttime(IS + m*(JS-1)) + 0.5*h*(slwns(IS+1 + m*(JS-1))+slwns(IS + m*(JS-1)))
      IF(JS.LT.N) ttime(IS + m*((JS-1)+1)) = ttime(IS + m*(JS-1)) + 0.5*h*(slwns(IS + m*((JS-1)+1))+slwns(IS + m*(JS-1)))
      IF(IS.GT.1) ttime(IS-1 + m*(JS-1)) = ttime(IS + m*(JS-1)) + 0.5*h*(slwns(IS-1 + m*(JS-1))+slwns(IS + m*(JS-1)))
      IF(JS.GT.1) ttime(IS + m*((JS-1)-1)) = ttime(IS + m*(JS-1)) + 0.5*h*(slwns(IS + m*((JS-1)-1))+slwns(IS + m*(JS-1)))
C
C TRAVEL TIME FOR R = 1
C
      nc = 0
      T0 = ttime(IS + m*(JS-1))
      S0 = slwns(IS + m*(JS-1))
      IF(JS.LT.N.AND.IS.GT.1) THEN
        T1 = ttime(IS + m*((JS-1)+1)) 
        T2 = ttime(IS-1 + m*(JS-1))
        slo = savg(m,n,slwns,is-1,js)
        slo1= slo
        if((js+2).le.n) slo1 = savg(m,n,slwns,is-1,js+1)
        slo2= slo
        if((is-2).ge.1) slo2 = savg(m,n,slwns,is-2,js)
        CALL ROOTS(ns,nc,1,T0,T1,T2,SLO,slo1,slo2,H,T3,index)
        ttime(IS-1 + m*(JS)) = T3
      END IF
      IF(JS.GT.1.AND.IS.GT.1) THEN
      S0 = slwns(IS + m*(JS-1))
        T1 = ttime(IS + m*(JS-2))
        T2 = ttime(IS-1 + m*(JS-1))
        slo = savg(m,n,slwns,is-1,js-1)
        slo1= slo
        if((js-2).ge.1) slo1 = savg(m,n,slwns,is-1,js-2)
        slo2= slo
        if((is-2).ge.1) slo2 = savg(m,n,slwns,is-2,js-1)
        CALL ROOTS(ns,nc,1,T0,T1,T2,SLO,slo1,slo2,H,T3,index)
        ttime(IS-1 + m*(JS-2)) = T3
      END IF
      IF(JS.GT.1.AND.IS.LT.M) THEN
      S0 = slwns(IS + m*(JS-1))
        T1 = ttime(IS + m*(JS-2))
        T2 = ttime(IS+1 + m*(JS-1))
        slo = savg(m,n,slwns,is,js-1)
        slo1= slo
        if((js-2).ge.1) slo1 = savg(m,n,slwns,is,js-2)
        slo2= slo
        if((is+2).le.m) slo2 = savg(m,n,slwns,is+1,js-1)
        CALL ROOTS(ns,nc,1,T0,T1,T2,SLO,slo1,slo2,H,T3,index)
        ttime(IS+1 + m*(JS-2)) = T3
      END IF
      IF(JS.LT.N.AND.IS.LT.M) THEN
      S0 = slwns(IS + m*(JS-1))
        T1 = ttime(IS + m*(JS))
        T2 = ttime(IS+1 + m*(JS-1))
        slo = savg(m,n,slwns,is,js)
        slo1= slo
        if((js+2).le.n) slo1 = savg(m,n,slwns,is,js+1)
        slo2= slo
        if((is+2).le.m) slo2 = savg(m,n,slwns,is+1,js)
        CALL ROOTS(ns,nc,1,T0,T1,T2,SLO,slo1,slo2,H,T3,index)
        ttime(IS+1 + m*(JS)) = T3
      END IF
      nc = 0

c      print *,'source done'
C 
C  TRAVEL TIME FOR R > 1
C
      IR = IS-1
      IF((M-IS).GT.IR) IR = M-IS
      IF((N-JS).GT.IR) IR = N-JS
      IF((JS-1).GT.IR) IR = JS-1
C

c      print *,'starting, IR= ',IR

      DO 90 I=2,IR
C
C                    R I G H T   P A R T
C

        JPOS = JS+I
        IF(JPOS.GT.N) GO TO 2
        JB = IS+I-1
        IF(JB.GT.M) JB = M
        JU = IS-I+1
        IF(JU.LT.1) JU = 1

c      print *,'starting right, JPOS= ',JPOS,'JB= ',JB,'JU= ',JU

        DO J=JU,JB
          TI(J) = ttime(J + m*(JPOS-2))
          JM(J) = J
        enddo
        call sort(ju,jb,ti,jm)
        do 30 j=ju,jb        
c
         index = 0
         if(jm(j).eq.m) then
          slo = savg(m,n,slwns,jm(j)-1,jpos-1)
          if(ttime(jm(j)-1 + m*(jpos-1)).eq.0.0) then
            call fd1(ttime(jm(j) + m*(jpos-2)),ttime(jm(j)-1 + m*(jpos-2)),ttime(jm(j)-1 + m*(jpos-2)),slo,slo,h,t3)
          else
            slo1= slo
            if((jpos+1).le.n) slo1 = savg(m,n,slwns,jm(j)-1,jpos)
            call roots(ns,ncode,i,ttime(jm(j)-1 + m*(jpos-2)),ttime(jm(j)-1 + m*(jpos-1)),ttime(jm(j) + m*(jpos-2)),slo,slo1,slo,h,t3,index)
            if(index.eq.1) id =  1
            if(index.eq.2) id = -1
          endif
          ttime(jm(j) + m*(jpos-1)) = t3
          goto 35
         endif
c
         if(jm(j).eq.1) then
          slo = savg(m,n,slwns,jm(j),jpos-1)
          if(ttime(jm(j)+1 + m*(jpos-1)).eq.0.0) then
            call fd1(ttime(jm(j) + m*(jpos-2)),ttime(jm(j)+1 + m*(jpos-2)),ttime(jm(j)+1 + m*(jpos-2)),slo,slo,h,t3)
          else
            slo1 = slo
            if((jpos+1).le.n) slo1 = savg(m,n,slwns,jm(j),jpos)
            call roots(ns,ncode,i,ttime(jm(j)+1 + m*(jpos-2)),ttime(jm(j)+1 + m*(jpos-1)),ttime(jm(j) + m*(jpos-2)),slo,slo1,slo,h,t3,index)
            if(index.eq.1) id = -1
            if(index.eq.2) id =  1
          endif
          ttime(jm(j) + m*(jpos-1)) = t3
          goto 35
         endif
c
         if(ttime(jm(j)+1 + m*(jpos-1)).eq.0.0.and.ttime(jm(j)-1 + m*(jpos-1)).eq.0.0)then
          slo1 = savg(m,n,slwns,jm(j)-1,jpos-1)
          slo2 = savg(m,n,slwns,jm(j)  ,jpos-1)
          call fd1(ttime(jm(j) + m*(jpos-2)),ttime(jm(j)-1 + m*(jpos-2)),ttime(jm(j)+1 + m*(jpos-2)),slo1,slo2,h,t3)
          ttime(jm(j) + m*(jpos-1)) = t3
          goto 30
         endif         
c        
         t3 = 1e+10         
         if(ttime(jm(j)+1 + m*(jpos-1)).gt.0.0)then
          slo = savg(m,n,slwns,jm(j),jpos-1)
          slo1 = slo
          if((jpos+1).le.n) slo1 = savg(m,n,slwns,jm(j),jpos)
          slo2 = savg(m,n,slwns,jm(j)-1,jpos-1)
          call roots(ns,ncode,i,ttime(jm(j)+1 + m*(jpos-2)),ttime(jm(j)+1 + m*(jpos-1)),ttime(jm(j) + m*(jpos-2)),slo,slo1,slo2,h,tsem,inde)
          if(tsem.lt.t3) then
            t3 = tsem
            index = inde
            if(index.eq.1) id= -1
            if(index.eq.2) id=  1
          endif
         endif
         if(ttime(jm(j)-1 + m*(jpos-1)).gt.0.0)then
          slo = savg(m,n,slwns,jm(j)-1,jpos-1)
          slo1 = slo
          if((jpos+1).le.n) slo1 = savg(m,n,slwns,jm(j)-1,jpos)
          slo2 = savg(m,n,slwns,jm(j),jpos-1)
          call roots(ns,ncode,i,ttime(jm(j)-1 + m*(jpos-2)),ttime(jm(j)-1 + m*(jpos-1)),ttime(jm(j) + m*(jpos-2)),slo,slo1,slo2,h,tsem,inde)
          if(tsem.lt.t3) then
            t3 = tsem
            index = inde
            if(index.eq.1) id =  1
            if(index.eq.2) id = -1
          endif
         endif
         ttime(jm(j) + m*(jpos-1)) = t3
35         if(index.eq.2) then
           if(id.eq. 1) nb=jb
           if(id.eq.-1) nb=ju
           do k=jm(j)+id,nb,id
             slo = 0.25*(slwns(k + m*(jpos-1))+slwns(k-id + m*(jpos-1))+slwns(k + m*(jpos-2))+slwns(k-id + m*(jpos-2)))
             call fd2(inde,ttime(k-id + m*(jpos-2)),ttime(k + m*(jpos-2)),ttime(k-id + m*(jpos-1)),slo,slo,slo,h,t3)
             if(t3.lt.ttime(k + m*(jpos-1))) then
               ttime(k + m*(jpos-1)) = t3
             else
               index = 0
               goto 30
             endif
           enddo
         endif
         if(index.eq.1) then
           do k=jpos-1,1,-1
             slo = 0.25*(slwns(jm(j) + m*(k-1))+slwns(jm(j)-id + m*(k-1))+ slwns(jm(j) + m*(k))+slwns(jm(j)-id + m*(k)))
             call fd2(inde,ttime(jm(j)-id + m*(k)),ttime(jm(j) + m*(k)),ttime(jm(j)-id + m*(k-1)),slo,slo,slo,h,t3)
             if(t3.lt.ttime(jm(j) + m*(k-1))) then
               ttime(jm(j) + m*(k-1))=t3
             else
               index = 0
               goto 30
             endif
           enddo   
         endif
c
   30    continue   

c      print *,'right done'

C
C                 A B O V E    P A R T
C 
    2   IPOS = IS-I
        IF(IPOS.LT.1) GO TO 3

c      print *,'starting top'

        JR = JS+I
        IF(JR.GT.N) JR = N
        JL = JS-I+1
        IF(JL.LT.1) JL = 1
        DO 41 J=JL,JR
          TI(J) = ttime(IPOS+1 + m*(J-1))
          JM(J) = J
   41   CONTINUE
        call sort(jl,jr,ti,jm)
        do 40 j=jl,jr
c
         index = 0
         if(jm(j).eq.n) then
          slo = savg(m,n,slwns,ipos,jm(j)-1)
          if(ttime(ipos + m*(jm(j)-2)).eq.0.0) then
            call fd1(ttime(ipos+1 + m*(jm(j)-1)),ttime(ipos+1 + m*(jm(j)-2)),ttime(ipos+1 + m*(jm(j)-2)),slo,slo,h,t3)
          else
            slo2 = slo
            if((ipos-1).ge.1) slo2 = savg(m,n,slwns,ipos-1,jm(j)-1)
            call roots(ns,ncode,i,ttime(ipos+1 + m*(jm(j)-2)),ttime(ipos+1 + m*(jm(j)-1)),ttime(ipos + m*(jm(j)-2)),slo,slo,slo2,h,t3,index)
            if(index.eq.1) id = -1
            if(index.eq.2) id =  1
          endif
          ttime(ipos + m*(jm(j)-1)) = t3
          goto 43
         endif
c
         if(jm(j).eq.1) then
          slo = savg(m,n,slwns,ipos,jm(j))
          if(ttime(ipos + m*(jm(j))).eq.0.0) then
            call fd1(ttime(ipos+1 + m*(jm(j)-1)),ttime(ipos+1 + m*(jm(j))),ttime(ipos+1 + m*(jm(j))),slo,slo,h,t3)
          else
            slo2 = slo
            if((ipos-1).ge.1) slo2 = savg(m,n,slwns,ipos-1,jm(j))
            call roots(ns,ncode,i,ttime(ipos+1 + m*(jm(j))),ttime(ipos+1 + m*(jm(j)-1)),ttime(ipos + m*(jm(j))),slo,slo,slo2,h,t3,index)
            if(index.eq.1) id =  1
            if(index.eq.2) id = -1
          endif
          ttime(ipos + m*(jm(j)-1)) = t3
          goto 43
         endif
c
         if(ttime(ipos + m*(jm(j))).eq.0.0.and.ttime(ipos + m*(jm(j)-2)).eq.0.0)then
          slo1 = savg(m,n,slwns,ipos,jm(j)-1)
          slo2 = savg(m,n,slwns,ipos,jm(j))
          call fd1(ttime(ipos+1 + m*(jm(j)-1)),ttime(ipos+1 + m*(jm(j)-2)),ttime(ipos+1 + m*(jm(j))),slo1,slo2,h,t3)
          ttime(ipos + m*(jm(j)-1)) = t3
          goto 40
         endif         
c        
         t3 = 1e+10
         if(ttime(ipos + m*(jm(j))).gt.0.0)then
          slo  = savg(m,n,slwns,ipos,jm(j))
          slo1 = savg(m,n,slwns,ipos,jm(j)-1)
          slo2 = slo
          if((ipos-1).ge.1) slo2 = savg(m,n,slwns,ipos-1,jm(j))
          call roots(ns,ncode,i,ttime(ipos+1 + m*(jm(j))),ttime(ipos+1 + m*(jm(j)-1)),ttime(ipos + m*(jm(j))),slo,slo1,slo2,h,tsem,inde)
          if(tsem.lt.t3) then
            t3 = tsem
            index = inde
            if(index.eq.1) id =  1
            if(index.eq.2) id = -1
          endif
         endif
         if(ttime(ipos + m*(jm(j)-2)).gt.0.0)then
          slo  = savg(m,n,slwns,ipos,jm(j)-1)
          slo1 = savg(m,n,slwns,ipos,jm(j))
          slo2 = slo
          if((ipos-1).ge.1) slo2 = savg(m,n,slwns,ipos-1,jm(j)-1)
          call roots(ns,ncode,i,ttime(ipos+1 + m*(jm(j)-2)),ttime(ipos+1 + m*(jm(j)-1)),ttime(ipos + m*(jm(j)-2)),slo,slo1,slo2,h,tsem,inde)
          if(tsem.lt.t3) then
            t3=tsem
            index = inde
            if(index.eq.1) id = -1
            if(index.eq.2) id =  1
          endif
         endif
         ttime(ipos + m*(jm(j)-1)) = t3
   43    if(index.eq.1) then
           if(id.eq. 1) nb=jr
           if(id.eq.-1) nb=jl
           do k=jm(j)+id,nb,id
             slo = 0.25*(slwns(ipos   + m*(k-1))+slwns(ipos   + m*(k-id-1))+slwns(ipos+1 + m*(k-1))+slwns(ipos+1 + m*(k-id-1)))
             call fd2(inde,ttime(ipos+1 + m*(k-id-1)),ttime(ipos+1 + m*(k-1)),ttime(ipos + m*(k-id-1)),slo,slo,slo,h,t3)
             if(t3.lt.ttime(ipos + m*(k-1))) then
               ttime(ipos + m*(k-1)) = t3
             else
               goto 40
             endif
           enddo
         endif
         if(index.eq.2) then
           do k=ipos+1,m
             slo = 0.25*(slwns(k   + m*(jm(j)-1))+slwns(k   + m*(jm(j)-id-1))+slwns(k-1 + m*(jm(j)-1))+slwns(k-1 + m*(jm(j)-id-1)))
             call fd2(inde,ttime(k-1 + m*(jm(j)-id-1)),ttime(k-1 + m*(jm(j)-1)),ttime(k + m*(jm(j)-id-1)),slo,slo,slo,h,t3)
             if(t3.lt.ttime(k + m*(jm(j)-1))) then
               ttime(k + m*(jm(j)-1))=t3
             else
               goto 40
             endif
           enddo   
         endif
c
   40   continue   

c      print *,'top done'

C
C                   L E F T   P A R T
C
   3    JPOS = JS-I
        IF(JPOS.LT.1) GO TO 4

c      print *,'starting left'

        JB = IS+I-1
        IF(JB.GT.M) JB = M
        JU = IS-I
        IF(JU.LT.1) JU = 1
        DO J=JU,JB
          TI(J) = ttime(J + m*(JPOS))
          JM(J) = J
        enddo
        call sort(ju,jb,ti,jm)
        do 50 j=ju,jb
c        
         index =0
         if(jm(j).eq.m) then
          slo = savg(m,n,slwns,jm(j)-1,jpos)
          if(ttime(jm(j)-1 + m*(jpos-1)).eq.0.0) then
            call fd1(ttime(jm(j) + m*(jpos)),ttime(jm(j)-1 + m*(jpos)),ttime(jm(j)-1 + m*(jpos)),slo,slo,h,t3)
          else
            slo1= slo
            if((jpos-1).ge.1) slo1 = savg(m,n,slwns,jm(j)-1,jpos-1)
            call roots(ns,ncode,i,ttime(jm(j)-1 + m*(jpos)),ttime(jm(j)-1 + m*(jpos-1)),ttime(jm(j) + m*(jpos)),slo,slo1,slo,h,t3,index)
            if(index.eq.1) id =  1
            if(index.eq.2) id = -1
          endif
          ttime(jm(j) + m*(jpos-1)) = t3
          goto 55
         endif
c
         if(jm(j).eq.1) then
          slo = savg(m,n,slwns,jm(j),jpos)
          if(ttime(jm(j)+1 + m*(jpos-1)).eq.0.0) then
            call fd1(ttime(jm(j) + m*(jpos)),ttime(jm(j)+1 + m*(jpos)),ttime(jm(j)+1 + m*(jpos)),slo,slo,h,t3)
          else
            slo1 = slo
            if((jpos-1).ge.1) slo1 = savg(m,n,slwns,jm(j),jpos-1)
            call roots(ns,ncode,i,ttime(jm(j)+1 + m*(jpos)),ttime(jm(j)+1 + m*(jpos-1)),ttime(jm(j) + m*(jpos)),slo,slo1,slo,h,t3,index)
            if(index.eq.1) id = -1
            if(index.eq.2) id =  1
          endif
          ttime(jm(j) + m*(jpos-1)) = t3
          goto 55
         endif
c
         if(ttime(jm(j)+1 + m*(jpos-1)).eq.0.0.and.ttime(jm(j)-1 + m*(jpos-1)).eq.0.0)then
          slo1 = savg(m,n,slwns,jm(j)-1,jpos)
          slo2 = savg(m,n,slwns,jm(j),jpos)
          call fd1(ttime(jm(j) + m*(jpos)),ttime(jm(j)-1 + m*(jpos)),ttime(jm(j)+1 + m*(jpos)),slo1,slo2,h,t3)
          ttime(jm(j) + m*(jpos-1)) = t3
          goto 50
         endif  
c        
         t3 = 1e+10
         if(ttime(jm(j)+1 + m*(jpos-1)).gt.0.0)then
          slo  = savg(m,n,slwns,jm(j),jpos)
          slo1 = slo
          if((jpos-1).ge.1) slo1 = savg(m,n,slwns,jm(j),jpos-1)
          slo2 = savg(m,n,slwns,jm(j)-1,jpos)
          call roots(ns,ncode,i,ttime(jm(j)+1 + m*(jpos)),ttime(jm(j)+1 + m*(jpos-1)),ttime(jm(j) + m*(jpos)),slo,slo1,slo2,h,tsem,inde)
          if(tsem.lt.t3) then
            t3 = tsem
            index = inde
            if(index.eq.1) id = -1
            if(index.eq.2) id =  1
          endif
         endif
         if(ttime(jm(j)-1 + m*(jpos-1)).gt.0.0) then
          slo  = savg(m,n,slwns,jm(j)-1,jpos)
          slo1 = slo
          if((jpos-1).ge.1) slo1 = savg(m,n,slwns,jm(j)-1,jpos-1)
          slo2 = savg(m,n,slwns,jm(j),jpos)
          call roots(ns,ncode,i,ttime(jm(j)-1 + m*(jpos)),ttime(jm(j)-1 + m*(jpos-1)),ttime(jm(j) + m*(jpos)),slo,slo1,slo2,h,tsem,inde)
          if(tsem.lt.t3) then
            t3 = tsem
            index = inde
            if(index.eq.1) id =  1
            if(index.eq.2) id = -1
          endif
         endif
         ttime(jm(j) + m*(jpos-1)) = t3
   55    if(index.eq.2) then
           if(id.eq. 1) nb=jb
           if(id.eq.-1) nb=ju
           do k=jm(j)+id,nb,id
             slo = 0.25*(slwns(k + m*(jpos-1))+slwns(k-id + m*(jpos-1))+slwns(k + m*(jpos))+slwns(k-id + m*(jpos)))
             call fd2(inde,ttime(k-id + m*(jpos)),ttime(k + m*(jpos)),ttime(k-id + m*(jpos-1)),slo,slo,slo,h,t3)
             if(t3.lt.ttime(k + m*(jpos-1))) then
               ttime(k + m*(jpos-1)) = t3
             else
               index = 0
               goto 50
             endif
           enddo
         endif
         if(index.eq.1) then
           do k=jpos+1,n
             slo = 0.25*(slwns(jm(j) + m*(k-1))+slwns(jm(j)-id + m*(k-1))+slwns(jm(j) + m*(k-2))+slwns(jm(j)-id + m*(k-2)))
             call fd2(inde,ttime(jm(j)-id + m*(k-2)),ttime(jm(j) + m*(k-2)),ttime(jm(j)-id + m*(k-1)),slo,slo,slo,h,t3)
             if(t3.lt.ttime(jm(j) + m*(k-1))) then
               ttime(jm(j) + m*(k-1))=t3
             else
               index = 0
               goto 50
             endif
           enddo   
         endif
c
   50   continue 

c      print *,'left done'

C
C                 B E L O W   P A R T
C
   4    IPOS = IS+I
        IF(IPOS.GT.M) GO TO 90
        JL = JS-I
        IF(JL.LT.1) JL = 1
        JR = JS+I
        IF(JR.GT.N) JR = N

c      print *,'starting bottom, IPOS= ',IPOS,'JL= ',JL,'JR= ',JR

        DO 61 J=JL,JR
          TI(J) = ttime(IPOS-1 + m*(J-1))
          JM(J) = J
   61   CONTINUE
        call sort(jl,jr,ti,jm)
        do 60 j=jl,jr
         index = 0
         if(jm(j).eq.n) then
          slo = savg(m,n,slwns,ipos-1,jm(j)-1)
            if(ttime(ipos + m*(jm(j)-2)).eq.0.0) then
            call fd1(ttime(ipos-1 + m*(jm(j)-1)),ttime(ipos-1 + m*(jm(j)-2)),ttime(ipos-1 + m*(jm(j)-2)),slo,slo,h,t3)
          else
            slo2 = slo
            if((ipos+1).le.m) slo2 = savg(m,n,slwns,ipos,jm(j)-1)
            call roots(ns,ncode,i,ttime(ipos-1 + m*(jm(j)-2)),ttime(ipos-1 + m*(jm(j)-1)),ttime(ipos + m*(jm(j)-2)),slo,slo,slo2,h,t3,index)
            if(index.eq.1) id = -1
            if(index.eq.2) id =  1
          endif
          ttime(ipos + m*(jm(j)-1)) = t3
          goto 63
         endif
c
         if(jm(j).eq.1) then
          slo = savg(m,n,slwns,ipos-1,jm(j))
          if(ttime(ipos + m*(jm(j))).eq.0.0) then
            call fd1(ttime(ipos-1 + m*(jm(j)-1)),ttime(ipos-1 + m*(jm(j))),ttime(ipos-1 + m*(jm(j))),slo,slo,h,t3)
          else
            slo2 = slo
            if((ipos+1).le.m) slo2 = savg(m,n,slwns,ipos,jm(j))
            call roots(ns,ncode,i,ttime(ipos-1 + m*(jm(j))),ttime(ipos-1 + m*(jm(j)-1)),ttime(ipos + m*(jm(j))),slo,slo,slo2,h,t3,index)
            if(index.eq.1) id =  1
            if(index.eq.2) id = -1
          endif
          ttime(ipos + m*(jm(j)-1)) = t3
          goto 63
         endif
c
         if(ttime(ipos + m*(jm(j))).eq.0.0.and.ttime(ipos + m*(jm(j)-2)).eq.0.0)then
          slo1 = savg(m,n,slwns,ipos-1,jm(j)-1)
          slo2 = savg(m,n,slwns,ipos-1,jm(j))
          call fd1(ttime(ipos-1 + m*(jm(j)-1)),ttime(ipos-1 + m*(jm(j)-2)),ttime(ipos-1 + m*(jm(j))),slo1,slo2,h,t3)
          ttime(ipos + m*(jm(j)-1)) = t3
          goto 60
         endif         
c        
         t3 = 1e+10
         if(ttime(ipos + m*(jm(j))).gt.0.0)then
          slo  = savg(m,n,slwns,ipos-1,jm(j))
          slo1 = savg(m,n,slwns,ipos-1,jm(j)-1)
          slo2 = slo
          if((ipos+1).le.m) slo2 = savg(m,n,slwns,ipos,jm(j))
          call roots(ns,ncode,i,ttime(ipos-1 + m*(jm(j))),ttime(ipos-1 + m*(jm(j)-1)),ttime(ipos + m*(jm(j))),slo,slo1,slo2,h,tsem,inde)
          if(tsem.lt.t3) then
            t3 = tsem
            index = inde
            if(index.eq.1) id =  1
            if(index.eq.2) id = -1
          endif
         endif
c
         if(ttime(ipos + m*(jm(j)-2)).gt.0.0)then
          slo  = savg(m,n,slwns,ipos-1,jm(j)-1)
          slo1 = savg(m,n,slwns,ipos-1,jm(j))
          slo2 = slo
          if((ipos+1).le.m) slo2 = savg(m,n,slwns,ipos,jm(j)-1)
          call roots(ns,ncode,i,ttime(ipos-1 + m*(jm(j)-2)),ttime(ipos-1 + m*(jm(j)-1)),ttime(ipos + m*(jm(j)-2)),slo,slo1,slo2,h,tsem,inde)
          
          if(tsem.lt.t3) then
            t3 = tsem
            index = inde
            if(index.eq.1) id = -1
            if(index.eq.2) id =  1
          endif
         endif
         ttime(ipos + m*(jm(j)-1)) = t3
c         goto 60
   63    if(index.eq.1) then
           if(id.eq. 1) nb=jr
           if(id.eq.-1) nb=jl
           do k=jm(j)+id,nb,id
             slo = 0.25*(slwns(ipos   + m*(k-1))+slwns(ipos   + m*(k-id-1))+slwns(ipos-1 + m*(k-1))+slwns(ipos-1 + m*(k-id-1)))
             call fd2(inde,ttime(ipos-1 + m*(k-id-1)),ttime(ipos-1 + m*(k-1)),ttime(ipos + m*(k-id-1)),slo,slo,slo,h,t3)
             if(t3.lt.ttime(ipos + m*(k-1))) then
               ttime(ipos + m*(k-1)) = t3
             else
               index = 0
               goto 60
             endif
           enddo
         endif
         if(index.eq.2) then
           do k=ipos-1,1,-1
             slo = 0.25*(slwns(k   + m*(jm(j)-1))+slwns(k   + m*(jm(j)-id-1))+slwns(k+1 + m*(jm(j)-1))+slwns(k+1 + m*(jm(j)-id-1)))
             call fd2(inde,ttime(k+1 + m*(jm(j)-id-1)),ttime(k+1 + m*(jm(j)-1)),ttime(k + m*(jm(j)-id-1)),slo,slo,slo,h,t3)
             if(t3.lt.ttime(k + m*(jm(j)-1))) then
               ttime(k + m*(jm(j)-1))=t3
             else
               index = 0
               goto 60
             endif
           enddo   
         endif
c
   60  continue   

c      print *,'bottom done'

C        
   90 CONTINUE
      return
      end
c
      function savg(m,n,slwns,i,j)
      implicit real*8 (a-h,o-z)
      dimension slwns(1000)
      savg = 0.25*(slwns(i+m*(j-1))+slwns(i+1+m*(j-1))+slwns(i+m*(j))+slwns(i+1+m*(j)))
      return
      end
c
      SUBROUTINE SORT(MIN,MAX,TI,JM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TI(400), JM(400)
C
      do 10 i=min,max-1
      do j = i+1,max
       if(ti(j).lt.ti(i)) then 
        se    = ti(i)
        ti(i) = ti(j)
        ti(j) = se
        ise   = jm(i)
        jm(i) = jm(j)
        jm(j) = ise
       endif
      enddo
   10 continue
      return
      end
c
      SUBROUTINE ROOTS(ns,ncode,ii,T0,T1,T2,SLO,slo1,slo2,H,T3,index)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX Z1,Z2
C 
      if(ii.gt.ns) then
        call fd2(index,t0,t1,t2,slo,slo1,slo2,h,t3)
        return
      endif
      BC1 = 2.0*(T0/slo**2.0 - T1/slo**2.0)
      BC2 = 2.0*(T0/slo**2.0 - T2/slo**2.0)
      CC1 = (T1/slo)**2.0 - (T0/slo)**2.0 - H**2.0
      CC2 = (T2/slo)**2.0 - (T0/slo)**2.0 - H**2.0
      C = BC1*BC1+BC2*BC2-(2.0*H/slo)**2
      D = 2.0*(BC1*CC1+BC2*CC2)+8.0*T0*(H/slo)**2
      E = CC1*CC1+CC2*CC2 - (2.0*H*T0/slo)**2 
      if(abs(c).lt.1e-05) then
        call fd2(index,t0,t1,t2,slo,slo1,slo2,h,t3)
        return
      endif
      CALL QUADRATIC(E,D,C,Z1,Z2)
      ts1 = real(z1)
      if(abs(ts1).lt.1e-04) ts1=0.0
      ts2 = real(z2)
      if(abs(ts2).lt.1e-04) ts2=0.0
      if(abs(aimag(z1)).gt.1e-4.or.abs(aimag(z2)).gt.1e-4)then
c        write(*,*)'complex roots'
        call fd2(index,t0,t1,t2,slo,slo1,slo2,h,t3)
        return
      endif
      if(abs(ts2).lt.abs(ts1)) then
        ts  = ts1
        ts1 = ts2
        ts2 = ts
      endif
      T3 = (T2-TS1)*(T2-TS1)/(slo**2)+(T1-TS1)*(T1-TS1)/(slo**2)
      T3 = TS1+slo*sqrt(T3 -(T0-TS1)*(T0-TS1)/(Slo**2))
      xs =(t1-ts1)*(t1-ts1)/(slo**2)-(t0-ts1)*(t0-ts1)/(slo**2)
      xs = (xs-h*h)/(2.0*h)
      if(abs(xs).lt.1e-05) xs = 0.0
      zs =(t2-ts1)*(t2-ts1)/(slo**2)-(t0-ts1)*(t0-ts1)/(slo**2)
      zs = (zs-h*h)/(2.0*h)
      if(abs(zs).lt.1e-05) zs = 0.0
      if(xs.lt.0.0.or.zs.lt.0.0.or.ts1.lt.0.0) then
        call fd2(index,t0,t1,t2,slo,slo1,slo2,h,t3)       
c        ncode = 1
      end if
      RETURN
      END
C
C THE ROOTS OF 2TH ORDER POLYNOMIAL
C
      SUBROUTINE QUADRATIC(C,B,A,Z1,Z2)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX Z1,Z2,D1
c
      D2  = B*B-4.0*A*C
      IF(dABS(D2).LT.1E-05) D2 = 0.0
      IF(D2.LT.0.0) D1 = CMPLX(0.0,(-D2)**0.5)
      IF(D2.GE.0.0) D1 = CMPLX(D2**0.5,0.0)
      Z1 = (-B+D1)/(2.0*A)
      Z2 = (-B-D1)/(2.0*A)   
      RETURN
      END      
c
c  FINITE DIFFERENCE EQUATION FOR CALCULATING TRAVELTIME 
c  AFTER RELATIVE MINIMUM POINT (EQ. 6 )
c
      subroutine fd1(t0,t1,t2,slo1,slo2,h,t3)
      IMPLICIT REAL*8 (A-H,O-Z)
      slo = 0.5*(slo1+slo2)
      t3 = h*h*slo*slo-0.25*(t2-t1)*(t2-t1) 
      if(t3.lt.0.0) then
c        write(*,*) 'fd1',t3
        t3 = 1e+10
      else
      t3 = t0 + sqrt(t3)
      endif
c      t3 = 1e+10
      delt = t0-t1
      if(delt.ge.0.0.and.delt.le.(h*slo1/sqrt(2.))) then
        tsem = h*h*slo1*slo1-delt*delt
c        if(tsem.lt.0.0) write(*,*)tsem
        tsem=t0+sqrt(tsem)
        if(tsem.lt.t3) t3 = tsem
      endif
      delt = t0-t2
      if(delt.ge.0.0.and.delt.le.(h*slo2/sqrt(2.))) then
        tsem = h*h*slo2*slo2-delt*delt
c        if(tsem.lt.0.0) write(*,*)tsem
        tsem=t0+sqrt(tsem)
        if(tsem.lt.t3) t3 = tsem
      endif
      delt = t3-t0
      if(delt.gt.h*slo1) then
        tsem = t0+h*slo1
        if(tsem.lt.t3) then
          t3 = tsem
c          write(*,*)'halo11'
        endif
      endif
      if(delt.gt.h*slo2) then
        tsem = t0+h*slo2
        if(tsem.lt.t3) then
          t3 = tsem
c          write(*,*)'halo12'
        endif
      endif
      return
      end
c
c  FINITE DIFFERENCE EQUATION FOR CALCULATING TRAVELTIME
c  IF THERE ARE 3 CALCULATED POINTS
c
      subroutine fd2(index,t0,t1,t2,slo,slo1,slo2,h,t3)
      IMPLICIT REAL*8 (A-H,O-Z)
      index = 0
      eps = 1e-5
      t3 = 1e+10
      t3 = 2.0*h*h*slo*slo-(t2-t1)*(t2-t1)
c      write(*,*)h,slo,t0,t1,t2,t3
c      if(abs(t3).lt.1e-10) t3 = 0.0
      if(t2.lt.t0.or.t1.lt.t0) t3=1e+10
      if(t3.lt.0.0) then
c        write(*,*) 'fd2',t3
        t3 = 1e+10
      else
      t3 = t0 + sqrt(t3)
      index=4
      endif
c PODVIN & LECOMTE'S STENSIL
c      ts = t0 + h*slo*sqrt(2.0)
c      if(ts.lt.t3) t3 = ts
      delt = (t1-t0)
      if(delt.ge.0.0.and.delt.le.(h*slo/sqrt(2.))) then
        ts = h*h*slo*slo-delt*delt
c        if(abs(ts).lt.1e-10) ts=0.0
        if(ts.lt.0.0) ts=1e+10
        ts=t1+sqrt(ts)
        if(ts.lt.t3) then
          t3 = ts
        endif
      endif
      delt = (t2-t0)
      if(delt.ge.0.0.and.delt.le.(h*slo/sqrt(2.))) then
        ts = h*h*slo*slo-delt*delt
c        if(abs(ts).lt.1e-10) ts=0.0
        if(ts.lt.0.0) ts=1e+10
        ts=t2+sqrt(ts)
        if(ts.lt.t3) then
          t3 = ts
        endif
      endif
      delt = t3-t1
      if(delt.gt.h*slo1) then
        ts = t1+h*slo1
        if(ts.lt.t3) then
          index = 1
          t3 = ts
        endif
      endif
      if(delt.gt.h*slo) then
        ts = t1+h*slo
        if(ts.lt.t3) then
          index = 0
c          write(*,*)'1th'
          t3 = ts
        endif
      endif
      delt = t3-t2
      if(delt.gt.h*slo2) then
        ts = t2+h*slo2
        if(ts.lt.t3) then
          index = 2
          t3 = ts
        endif
      endif
      if(delt.gt.h*slo) then
        ts = t2+h*slo
        if(ts.lt.t3) then
          index = 0
c          write(*,*)'2th'
          t3 = ts
        endif
      endif 
      ts = t0 + sqrt(2.)*h*slo
      if(ts.lt.t3) then
        t3 = ts
        index = 0
      endif
      ts = t0 + sqrt(2.)*h*slo1
      if(ts.lt.t3) then
        t3 = ts
        index = 3
c        write(*,*)index,slo1
      endif
      ts = t0 + sqrt(2.)*h*slo2
      if(ts.lt.t3) then
        t3 = ts
        index = 3
      endif
c      write(*,*)t3
      return
      end
c
