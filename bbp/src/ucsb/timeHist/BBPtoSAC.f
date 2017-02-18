C converts TXT to SAC 

	character*72 fileinp
	character*72 name
	real syn(100000)
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
	name='output.sac'
	write(*,*) name
c SAC
          b0=0.0
          call wsac1(name,syn,nt,b0,dt,nerr)
	

	stop
	end
