PROGRAM bin2bbp

implicit none

real                           :: dt,dummy,m2cm
real,allocatable,dimension(:,:):: val
integer                        :: n_stat,scratch,npts,k,n,n_args
character(len=10)              :: nstat_char,n_char
character(len=256)             :: input_file,output_file

!------------------------------------------------------------------------------------------

! Get argument first
n_args = iargc()    
if (n_args .ne. 2) then
   print *, 'The program needs 2 arguments'
   stop
endif

m2cm = 100.0

call getarg(1,input_file)    !Name of input BIN file containing time-series 
call getarg(2,nstat_char)    !Number of stations

read(nstat_char,'(I10)') n_stat           

! detect record length
inquire(iolength=scratch) (/1.0,1.0,1.0/)

! open binary file
open(1,file=trim(input_file),status='old',access='direct',form='unformatted',recl=scratch)
   
! first line is reserved to: number of points, time-step, dummy
! Note: these quantities are assumed to be the same for all the stations
read(1,rec=1) npts, dt, dummy         

! allocate array for waveforms
allocate(val(npts,3)) 

do n=1,n_stat
                     
   ! reading file... (3 components at the same record line -- each record line is a time-step)
   do k=2,npts+1  
      read(1,rec= (n-1) * npts + k) val(k-1,1),val(k-1,2),val(k-1,3)
   enddo
    
   if (n < 10) then
      write(n_char,'(I1)') n
   elseif (n >= 10 .and. n < 100) then
      write(n_char,'(I2)') n
   else
      write(n_char,'(I3)') n
   endif 
        
   output_file = trim(n_char)//'.bbp'   
        
   ! open BBP file
   open(2,file=trim(output_file),status='unknown',form='formatted')
 
   ! write header
   do k=1,16
      write(2,'(A)') '#     header'
   enddo   
   write(2,'(A)') '# Column 1: Time (s)'
   write(2,'(A)') '# Column 2: North-south velocity (cm/s)'
   write(2,'(A)') '# Column 3: East-west velocity (cm/s)'
   write(2,'(A)') '# Column 4: Up-down velocity (cm/s)'
   write(2,'(A)') '#'
 
   ! writing output BBP file
   do k=1,npts  
      write(2,"(E12.6,TR4,E12.6,TR4,E12.6,TR4,E12.6)") dt*(k-1),m2cm*val(k,1),m2cm*val(k,2),m2cm*val(k,3) 
      !write(2,"(E12.6,TR4,E12.6,TR4,E12.6,TR4,E12.6)") dt*(k-1),val(k,1),val(k,2),val(k,3)
   enddo
   
   close(2)   
   
enddo

deallocate(val)

close(1)  !close file

   
END PROGRAM bin2bbp  
