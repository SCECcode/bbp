PROGRAM bbp2bin

implicit none

real                           :: dt,dummy,cm2m
real,allocatable,dimension(:,:):: val
integer                        :: n_stat,scratch,npts,k,n,n_args,ierr,head_lines
character(len=10)              :: nstat_char,n_char
character(len=256)             :: input_file,output_file,stat_name,dummy_line

!------------------------------------------------------------------------------------------

! Get argument first
n_args = iargc()    
if (n_args .ne. 3) then
   print *, 'The program needs 3 arguments'
   stop
endif

cm2m = 0.01

call getarg(1,input_file)    !Name of input file (ascii) specifying BBP files to read
call getarg(2,output_file)   !Name of output BIN file 
call getarg(3,nstat_char)    !Number of stations

read(nstat_char,'(I10)') n_stat           

! detect record length
inquire(iolength=scratch) (/1.0,1.0,1.0/)

! open binary file
open(1,file=trim(output_file),status='unknown',access='direct',form='unformatted',recl=scratch)
   
! open ascii file
open(2,file=trim(input_file),status='old',form='formatted') 

! find npts and dt for BBP files
read(2,*) stat_name; rewind(2)

open(3,file=trim(stat_name),status='old',form='formatted')

! find number of points and dt
head_lines=0
do
        read(3,*) dummy_line
        !This line checks to see if it starts with a '#' or '%';  if so, ignore
        if ((dummy_line(1:1) .ne. '#') .AND. (dummy_line(1:1) .ne. '%')) exit
        head_lines=head_lines+1
end do
!write(*,*) head_lines
!do n=1,22
!   read(3,*) dummy_line
!   write(*,*) dummy_line
!enddo
read(3,*) dt,dummy,dummy,dummy; rewind(3); dummy=1.0
do n=1,head_lines
   read(3,*) dummy_line
enddo
npts=0
do
   read(3,*,iostat=ierr)
   if (ierr == -1) exit
   npts=npts+1
enddo

!write(*,*) npts
rewind(3)                                          

! first line is reserved to: number of points, time-step, dummy
! Note: these quantities are assumed to be the same for all the stations
write(1,rec=1) npts, dt, dummy         

! allocate array for waveforms
allocate(val(npts,3)) 

do n=1,n_stat
           
   ! read file name        
   read(2,*) stat_name
   
   ! open file and read its content
   open(3,file=trim(stat_name),status='old',form='formatted')
   do k=1,head_lines
      read(3,*)
   enddo
   do k=1,npts
      read(3,*) dt,val(k,1),val(k,2),val(k,3)
   enddo   
   close(3)
           
   ! writing file... (3 components at the same record line -- each record line is a time-step)
   do k=2,npts+1  
      write(1,rec= (n-1) * npts + k) cm2m*val(k-1,1),cm2m*val(k-1,2),cm2m*val(k-1,3)
      !write(1,rec= (n-1) * npts + k) val(k-1,1),val(k-1,2),val(k-1,3)
   enddo
       
enddo

deallocate(val)

close(1); close(2)  !close file

   
END PROGRAM bbp2bin  
