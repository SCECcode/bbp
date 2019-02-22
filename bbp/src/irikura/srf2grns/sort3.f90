subroutine sort3(x,y,z,imax,iloc)

! sort -k 3, 2, 1

 implicit none
 integer::imax,i,iloc(imax),k,kk,itmp
 real,dimension(imax)::x,y,z,x0,y0,z0
 real::xtmp,ytmp,ztmp

  z0=z ; y0=y ; x0=x

  do i=1,imax
     iloc(i)=i
  enddo

  do i=1,imax-1
     ztmp=z(i) ; ytmp=y(i) ; xtmp=x(i) ; kk=i
     itmp=iloc(i)
     do k=i+1,imax
        if(ztmp>z(k)) then
           ztmp=z(k) ; ytmp=y(k) ; xtmp=x(k) ; kk=k
           itmp=iloc(k)
        else if(ztmp==z(k)) then
           if(ytmp>y(k)) then
              ztmp=z(k) ; ytmp=y(k) ; xtmp=x(k) ; kk=k
              itmp=iloc(k)
           else if(ytmp==y(k)) then
              if(xtmp>x(k)) then
                 ztmp=z(k) ; ytmp=y(k) ; xtmp=x(k) ; kk=k
                 itmp=iloc(k)
              endif
           endif
        endif
     enddo
     z(kk)=z(i) ; y(kk)=y(i) ; x(kk)=x(i) ; iloc(kk)=iloc(i)
     z(i)=ztmp ; y(i)=ytmp ; x(i)=xtmp ; iloc(i)=itmp
  enddo
  return
end subroutine sort3


                 






