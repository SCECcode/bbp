C converts TXT to SAC 

	character*128 fileinp
	character*128 sta,comp,name
	real syn(16000)
	integer nt
	real dt
	integer nerr,myfh,getlen
        real b0

        name = ''
	write(*,*) 'input file'
	read(*,'(1a)') fileinp
        write(*,*) fileinp
	

	open(17,file=fileinp,status='old')
	read(17,*) nt,dt
	read(17,*) (syn(k),k=1,nt)
	close(17)
	kk=1
1	if(fileinp(kk:kk).ne.'.') then
	kk=kk+1
	goto 1
	endif
	sta=fileinp(1:kk-1)
	comp=fileinp(kk+1:kk+4)
	name(1:kk+16)=fileinp(1:kk+12)//'.sac'
	write(*,*) name
c SAC
          b0=0.0
          call wsac1(name,syn,nt,b0,dt,nerr)
	

	stop
	end
