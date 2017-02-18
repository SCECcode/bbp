cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  complex.f: Dealing with complex numbers on a comput                 c
c                                                                      c
c  taken from: "Projects in Computational Physics" by Landau and Paez  c 
c	       copyrighted by John Wiley and Sons, New York            c      
c                                                                      c
c  written by: students in PH465/565, Computational Physics,           c
c	       at Oregon State University                              c
c              code copyrighted by RH Landau                           c
c  supported by: US National Science Foundation, Northwest Alliance    c
c                for Computational Science and Engineering (NACSE),    c
c                US Department of Energy 	                       c
c								       c
c  UNIX (DEC OSF, IBM AIX): f77 complex.f  			       c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       	Program complex
c
c complex numbers and functions
       	Implicit none
c declarations
        Complex*16 z, zsqrt, zlog,compOne
       	Real*8 i, pi, phi, x, y, zatan, zatan2
       	pi=3.1415926535897932385E0
	compOne=cmplx(-1.0,0.0)
	write(*,*) sqrt(compOne)
c write header for table
       	Write (*,10) 'phi','x','y','sqrt','log','atan','atan2'
       	Write (*,*) ' '
c loop for angle   
 	Do 100 i=0, 2.6, 0.1
           phi    = i * pi
c calculate carthesian representation 
           x      = cos(phi)
           y      = sin(phi)
           z      = cmplx(x,y)
c call functions
           zsqrt  = sqrt(z)
           zlog   = log(z)
           zatan  = atan(y/x)
           zatan2 = atan2(y,x)
c write results
           Write (*,20) i,'*pi',x,y,zsqrt,'i',zlog,'i',zatan,zatan2
 100   	Continue
 10    	Format (a4, 2a9, a14, a18, a14, a10)
 20    	Format (f3.1, a3, 3f9.4, f8.4, a1, f9.4, f8.4,a1, 2f9.4)
	Stop 'complex'
       	End

