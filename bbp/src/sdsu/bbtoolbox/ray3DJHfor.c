/*    punch   */

/*  from Vidale's slug3d
    edited by j.hole 02-08-90    add ability to accept as
       source some input values of travel times on the faces
       of the data volume  
    edited by j.hole 06-09-90    fix bugs in Vidale's code
       remove edges and corners, including
       them in the sides calculation        add ability to choose
       rate of speed of box expansion in each of 6 directions
       **** RENAMED slugjah.c
    edited by j.hole 17-09-90    replace source cube of constant
       velocity with a source cube of linear velocity gradient
    edited by j.hole 11-91       new source option: input a tt field
       that has already been calculated, start calculations on any
       wall of volume, change tt only if a smaller value is found
       (allows rays to travel back into expanding cube)
    MAJOR edits by j.hole 11/91  new fd stencils and stencil logic
       to attempt to allow head waves to travel parallel the faces 
       of the expanding cube (and at some other angles too)
       **** RENAMED bull3d.c
    MAJOR edits by j.hole 01/93  extended bull3d.c logic to 
       systematically test all possible operators
       **** RENAMED bully.c
    MAJOR edits by j.hole 12/93  removed bug in bully.c;  only use
       "new face" and "new edge" stencils on new faces and edges!
       **** RENAMED punch.c
    edited by j.hole 01/95  added ability to automatically detect 
       head waves travelling on the faces of the expanding cube 
       of already-calculated values; then restarts calculations 
       at appropriate wall of the model to allow the waves to 
       travel back into the (previous) expanding cube
****Hole, J.A., and B.C. Zelt, 1995.  3-D finite-difference reflection 
       traveltimes.  Geophys. J. Int., 121, 427-434.
       This reference describes the modifications for head waves 
       made 11-91, 01/93, 12/93, and 01/95.  
    edited by j.hole spring/98   allowed an x0,y0,z0 origin to be input
    edited by j.hole 05/99   double precision on source cube calculation
*/


/*  RECEIVED IN E-MAIL 12-02-90  */
/* 
Message inbox:3 -  Read
From:    <vid@flatte.ucsc.edu>
To:      <hole@geop.ubc.ca>, <vid@flatte.ucsc.edu>
Cc:      <vid@rupture.ucsc.edu>
*** NOTE THAT THESE ADDRESSES ARE OUTDATED ***
Subject: Re:  3-D Travel Times

I can send the code, but this week I am too busy to send sample input and 
output files.
The code is included below.
No special compiling options are necessary.
Two input files are necessary.
One, in ascii, has lines that look like:

nx=100
ny=100
nz=100
xs=10
ys=10
zs=10
h=1
timefile=run0.times
velfile=input.vel

(Notice that there are no spaces.)
The other file is the velocity field (an nx by ny by nz array) that is
written in C binary, and with the parameter file above would have the
name "input.vel".
When the program is finished, the traveltimes to all the grid points will
appear in a file called run0.times.
If the C code was compiled with the command "cc time3d.c -o time3d", and
the parameter file is named "time.par", the
program would run with the command "time3d par=time.par".
If the program gives the correct answer for a uniform wholespace, 
it is probably working.
Send me some mail if something doesn't work, and give me a call (408 459-4585)
if I don't answer the e-mail (our computers are in a state of flux).   */

/* PROGRAM TO CALCULATE TRAVEL TIMES IN 3D MEDIA */
/* author: John E. Vidale, started Sept. 11, 1986*/
/*  restarted, 8-11-88, with more accurate scheme */
/*  last revised, 2-89, by Glenn Nelson for 5x5x5 source */
/* UC Santa Cruz, all rights reserved */

#include    <stdio.h>
#include    <math.h>
#include    <fcntl.h>
#include    <string.h>

#define PI 3.141592654
#define SQR2 1.414213562
#define SQR3 1.732050808
#define SQR6 2.449489743
/* #define	NCUBE	2 */	/* 1 for 3x3x3, 2 for 5x5x5	*/
#define t0(x,y,z)   time0[nxy*(z) + nx*(y) + (x)]
#define s0(x,y,z)   slow0[nxy*(z) + nx*(y) + (x)]
#define	SQR(x)	((x) * (x))
#define	DIST(x,y,z,x1,y1,z1)	sqrt(SQR(x-(x1))+SQR(y-(y1)) + SQR(z-(z1)))

struct sorted
	{ float time; int i1, i2;};

/* FUNCTION DECLARATIONS	*/
int
	compar();
float
	ex0(), fd5(), fd6(), fd7(), fd8(),         /* VIDALE'S STENCILS */
        fdh3d(),fdhne(),fdh2d(),fdhnf();           /* HOLE'S STENCILS */

/*void raytracing_(hypo,grid,step,velfile,timefile)*/
void raytracing_(hypo,grid,step,PS_flag)
	
	float *hypo;
	int   *grid;
	float *step;
	/*char  *velfile;
	char  *timefile;
	float *timing;*/
	int   *PS_flag;
	
	
{
	/* NOTE THAT SEVERAL VARIABLES MUST BE SPECIFIED IN par=xxx FILE, 
	   WHILE OTHERS ARE OPTIONAL:  IF A mstpar STATEMENT READS THE 
	   VARIABLE BELOW, THEN THE VARIABLE IS REQUIRED;  IF A getpar 
	   STATEMENT IS USED BELOW, THEN THE VARIABLE IS OPTIONAL */
	int ac;
        int
		nx,		/* x-dimension of mesh (LEFT-TO-RIGHT) */
		ny,		/* y-dimension of mesh (FRONT-TO-BACK) */
		nz,		/* z-dimension of mesh  (TOP-TO-BOTTOM) */
		iplus=1,	/* rate of expansion of "cube" in the */
		iminus=1,	/*    plus/minus x/y/z direction */
		jplus=1,
		jminus=1,
		kplus=1,
		kminus=1,
		igrow,		/* counter for "cube" growth */
                srctype=1,      /* if 1, source is a point;
                                      2, source is on the walls of the data volume;
				      3, source on wall, time field known; */
		floatsrc=1,	/* if 0, source must be on a grid node 
				      1, source can lie between grid nodes */
                srcwall,        /* if 1, source on x=0 wall, if 2, on x=nx-1 wall
                                   if 3, source on y=0 wall, if 4, on y=ny-1 wall
                                   if 5, source on z=0 wall, if 6, on z=nz-1 wall */
		xs,		/* shot x position (in grid points) */
		ys,		/* shot y position */
		zs,		/* shot depth */
		xx, yy, zz,		/* Used to loop around xs, ys, zs coordinates	*/
		X1, X2, lasti, index, ii, i, j, k, radius, 
	        /*tfint, vfint, wfint, ofint,*/
	        wfint, ofint,
		nxy, nyz, nxz, nxyz, nwall,
		/* counters for the position of the sides of current cube */
		x1, x2, y1, y2, z1, z2,
		/* flags set to 1 until a side has been reached */
		dx1=1, dx2=1, dy1=1, dy2=1, dz1=1, dz2=1, rad0=1,
                /* maximum radius to compute */
	        maxrad,
		/* used in linear velocity gradient cube source */
		xxx,yyy,zzz,
                /* size of source cube     1 for 3x3x3, 2 for 5x5x5...	*/
                NCUBE=2,
		reverse=1,	/* will automatically do up to this number of
			reverse propagation steps to fix waves that travel 
			back into expanding cube */
		headpref=6,	/* if headpref starts > 0, will determine 
			model wall closest to source and will prefer to start
			reverse calculations on opposite wall */
		/* counters for detecting head waves on sides of current cube */
		head,headw[7];
	float
		h,		/* spatial mesh interval (units consistant with vel) */
		fxs,	/* shot position in X (in real units)*/
		fys,	/* shot position in Y (in real units)*/
		fzs,	/* shot position in Z (in real units)*/
		*slow0, *time0, *wall, *times,
		s000, guess, try, slo,
		/* real-unit coordinates of first grid point (for floatsrc=1) */
		x0=0., y0=0., z0=0.,
                /* maximum offset (real units) to compute */
	        maxoff = -1.,
		/* used to detect head waves:  if headwave operator decreases 
		   the previously-computed traveltime by at least 
		   headtest*<~time_across_cell> then the headwave counter is 
		   triggered */
		fhead,headtest=1.e-3,
		/* pointers to times and slownesses */
		*r0, *r1, *r2, *r3,
		*p0, *p1, *p2, *p3, *p4, *p5;
	double
		/* used in linear velocity gradient cube source */
		rx, ry, rz, dvx, dvy, dvz, dv, v0,
		rzc, rxyc, rz1, rxy1, rho, theta1, theta2;
	char
		/*velfile[80],	 file though which velocity structure is input */
		oldtfile[80],	/* file through which old travel times are input */
		/*timefile[80],	 file in which travel times appear at the end */
                wallfile[80];   /* file containing input wall values of traveltimes */

	char  *velfile=malloc(11);
	char  *timefile=malloc(12);
	FILE  *vfint, *tfint;

	/* ARRAY TO ORDER SIDE FOR SOLUTION IN THE RIGHT ORDER */
	struct sorted *sort;

	fprintf(stderr,"Starting slug3d: by J. Vidale, 1988, UCSC\n");
	fprintf(stderr,"Starting punch: by J. Hole, 1993, UBC-Stanford\n");
 
    
    fxs = hypo[0];
    fys = hypo[1];
    fzs = hypo[2];
    nx = grid[0];
    ny = grid[1];
    nz = grid[2];
    h =  *step;
    
    x0 = 0.;
    y0 = 0.;
    z0 = 0.;
    reverse = 0;
    maxoff = 1000.;
    
   
    fxs = (fxs-x0)/h;
	fys = (fys-y0)/h;
	fzs = (fzs-z0)/h;
	xs = (int)(fxs + 0.5);
	ys = (int)(fys + 0.5);
	zs = (int)(fzs + 0.5);

    if (*PS_flag == 1){
       //timefile = "time3d_P.out";
       //velfile =  "vel3d_P.bin";
       strcpy(timefile, "time3d_P.out");
       strcpy(velfile,  "vel3d_P.bin");
       printf("\nTEST1\n %s\n %s\n",timefile,velfile);
    }
    else if (*PS_flag == 2){
       //timefile = "time3d_S.out";
       //velfile =  "vel3d_S.bin";
       strcpy(timefile, "time3d_S.out");
       strcpy(velfile,   "vel3d_S.bin");
       printf("\nTEST2\n %s\n %s\n",timefile,velfile);
    }

    fprintf(stderr,"-- travel-times written to :%s\n",timefile);
    fprintf(stderr,"-- selected velocity file  :%s\n",velfile);


	/* SET MAXIMUM RADIUS TO COMPUTE */
	if (maxoff > 0.) {
	  maxrad = maxoff/h + 1;
	  fprintf(stderr,"WARNING: Computing only to max radius = %d\n",maxrad);
	}
	else   maxrad = 99999999;

	nxy = nx * ny;
	nyz = ny * nz;
	nxz = nx * nz;
	nxyz = nx * ny * nz;
	
	/* FORM AND FILL TT AND SLOWNESS ARRAYS */
	/*if((vfint=open(velfile,O_RDONLY,0664))<=1) {
		fprintf(stderr,"cannot open %s\n",velfile);
		exit(-1);
	}
        if((tfint=open(timefile,O_CREAT|O_WRONLY|O_TRUNC,0664))<=1) {
		fprintf(stderr,"cannot open %s\n",timefile);
		exit(-1);
	}*/
	vfint = fopen(velfile, "rb");
	if(vfint == NULL) {
	         fprintf(stderr,"cannot open %s %d\n",velfile,vfint);
	         exit(-1);
	}

	tfint = fopen(timefile, "wb");
	if(tfint == NULL) {
	        fprintf(stderr,"cannot open %s %d\n",timefile,tfint);
		exit(-1);
	}


	/* ALLOCATE MAIN AND ALTERNATE GRID FOR SLOWNESSES AND TIMES */
	slow0 = (float *) malloc(4*nxyz);
	time0 = (float *) malloc(4*nxyz*1.25);
	times = (float *) malloc(4*nxy);
	
	
	/* MAKE ARRAY SORT LARGE ENOUGH FOR ANY SIDE */
	if(nx <= ny && nx <= nz)  {
		sort = (struct sorted *) malloc(sizeof(struct sorted)*ny*nz);
		nwall = nyz;
	}
	else if(ny <= nx && ny <= nz)  {
		sort = (struct sorted *) malloc(sizeof(struct sorted)*nx*nz);
		nwall = nxz;
	}
	else  {
		sort = (struct sorted *) malloc(sizeof(struct sorted)*nx*ny);
		nwall = nxy;
	}
	wall = (float *) malloc(4*nwall);
	if( slow0 == NULL || time0 == NULL || sort == NULL || wall == NULL) {
		fprintf(stderr,"cannot allocate memory\n");
		exit(-1);
	}
	/* READ IN VELOCITY FILE */
	/*read(vfint,slow0,nxyz*4);*/
	fread(slow0,sizeof(float),nxyz*4,vfint);
	
	
	/* CONVERT TO SLOWNESS TIMES MESH SPACING */
	for(i=0;i<nxyz;i++) slow0[i] = h/slow0[i];

	/* SET TIMES TO DUMMY VALUE ***** JAH ***** BUG IN VIDALE'S CODE */
	for(i=0;i<nxyz*1.25;i++) time0[i] = 1.0e10;

	if (srctype == 1) {			/*  VIDALE'S POINT SOURCE */


		/* HOLE'S NEW LINEAR VELOCITY GRADIENT CUBE (APRIL 1991)*/
		v0 = h/s0(xs,ys,zs);
		for (xx = xs-NCUBE; xx <= xs+NCUBE; xx++) {
			if (xx < 0 || xx >= nx)	continue; 
			for (yy = ys-NCUBE; yy <= ys+NCUBE; yy++) {
				if (yy < 0 || yy >= ny)	continue; 
				for (zz = zs-NCUBE; zz <= zs+NCUBE; zz++) {
					if (zz < 0 || zz >= nz)	continue; 
					if (zz == zs)
					  dvz = 1/s0(xx,yy,zz+1)-1/s0(xs,ys,zs);
					else
					  dvz = (1/s0(xx,yy,zz)-1/s0(xs,ys,zs))/(zz-zs);
					dv = fabs(dvz);
					if (dv == 0.)  {
					  t0(xx,yy,zz) = s0(xs,ys,zs)*DIST(fxs,fys,fzs,xx,yy,zz);
					  continue;
					}
					rzc = -v0/dv;
					rx = h*(xx - fxs);
					ry = h*(yy - fys);
					rz = h*(zz - fzs);
					rz1 = rz*dvz/dv;
					rxy1 = sqrt(rx*rx+ry*ry);
				/*rxy1 = sqrt(rx*rx+ry*ry+rz*rz-rz1*rz1);*/
					if (rxy1<=h/1.e6)
					  t0(xx,yy,zz) = fabs(log((v0+dv*rz1)/v0)/dv);
					else {
					  rxyc = (rz1*rz1+rxy1*rxy1-2*rz1*rzc)/(2*rxy1);
					  rho = sqrt(rzc*rzc+rxyc*rxyc);
					  theta1 = asin(-rzc/rho);
					  /* can't handle asin(1.) ! */
					  if (fabs(rz1-rzc)>=rho)  rho=1.0000001*fabs(rz1-rzc);
					  theta2 = asin((rz1-rzc)/rho);
					  if (rxyc<0) theta1=PI-theta1;
					  if (rxyc<rxy1) theta2=PI-theta2;
					  t0(xx,yy,zz) = log(tan(theta2/2)/tan(theta1/2)) / dv;
				        }
				}
			}
		}

		/* SETS LOCATION OF THE SIDES OF THE CUBE	*/
		radius = NCUBE;
		if(xs > NCUBE) x1 = xs - (NCUBE + 1);
		else{ x1 = -1; dx1 = 0;}
		if(xs < nx-(NCUBE + 1)) x2 = xs + (NCUBE + 1);
		else{ x2 = nx; dx2 = 0;}
		if(ys > NCUBE) y1 = ys - (NCUBE + 1);
		else{ y1 = -1; dy1 = 0;}
		if(ys < ny-(NCUBE + 1)) y2 = ys + (NCUBE + 1);
		else{ y2 = ny; dy2 = 0;}
		if(zs > NCUBE) z1 = zs - (NCUBE + 1);
		else{ z1 = -1; dz1 = 0;}
		if(zs < nz-(NCUBE + 1)) z2 = zs + (NCUBE + 1);
		else{ z2 = nz; dz2 = 0;}
	}
	else if (srctype == 2) {		/*  HOLE'S EXTERNAL SOURCE */

		/* FILL IN WALLS' TIMES FROM EXTERNAL DATAFILE */
		read (wfint,wall,4*nwall);	/* READ X=0 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (k=0; k<nz; k++) {
				for (j=0; j<ny; j++) {
					t0(0,j,k) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ X=NX-1 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (k=0; k<nz; k++) {
				for (j=0; j<ny; j++) {
					t0(nx-1,j,k) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ Y=0 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (k=0; k<nz; k++) {
				for (i=0; i<nx; i++) {
					t0(i,0,k) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ Y=NY-1 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (k=0; k<nz; k++) {
				for (i=0; i<nx; i++) {
					t0(i,ny-1,k) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ Z=0 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (j=0; j<ny; j++) {
				for (i=0; i<nx; i++) {
					t0(i,j,0) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ Z=NZ-1 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (j=0; j<ny; j++) {
				for (i=0; i<nx; i++) {
					t0(i,j,nz-1) = wall[ii];
					ii++;
				}
			}
		}

		/* SET LOCATIONS OF SIDES OF THE CUBE SO THAT CUBE IS A FACE  */
		radius = 1;
		if (srcwall == 1)	x2=1;
		else	{  x2=nx;	dx2=0;  }
		if (srcwall == 2)	x1=nx-2;
		else	{  x1= -1;	dx1=0;  }
		if (srcwall == 3)	y2=1;
		else	{  y2=ny;	dy2=0;  }
		if (srcwall == 4)	y1=ny-2;
		else	{  y1= -1;	dy1=0;  }
		if (srcwall == 5)	z2=1;
		else	{  z2=nz;	dz2=0;  }
		if (srcwall == 6)	z1=nz-2;
		else	{  z1= -1;	dz1=0;  }
	}
	else if (srctype == 3) {                /*  HOLE'S REDO OLD TIMES */
	        /* READ IN OLD TIME FILE */
	        if (srctype == 3)  read(ofint,time0,nxyz*4);
		/* SET LOCATIONS OF SIDES OF THE CUBE SO THAT CUBE IS A FACE */
		radius = 1;
		if (srcwall == 1)	x2=1;
		else	{  x2=nx;	dx2=0;  }
		if (srcwall == 2)	x1=nx-2;
		else	{  x1= -1;	dx1=0;  }
		if (srcwall == 3)	y2=1;
		else	{  y2=ny;	dy2=0;  }
		if (srcwall == 4)	y1=ny-2;
		else	{  y1= -1;	dy1=0;  }
		if (srcwall == 5)	z2=1;
		else	{  z2=nz;	dz2=0;  }
		if (srcwall == 6)	z1=nz-2;
		else	{  z1= -1;	dz1=0;  }
	}
	else  {
		fprintf(stderr,"incorrect value of srctype = %d\n",srctype);
		exit(-1);
	}

	if (headpref>0) {	/* HOLE - PREFERRED REVERSE DIRECTION */
		head = nx*ny*nz;
		if (nx>5 && x2<=head)   {headpref=2;  head=x2;}
		if (nx>5 && (nx-1-x1)<=head)   {headpref=1;  head=nx-1-x1;}
		if (ny>5 && y2<=head)   {headpref=4;  head=y2;}
		if (ny>5 && (ny-1-y1)<=head)   {headpref=3;  head=ny-1-y1;}
		if (nz>5 && z2<=head)   {headpref=6;  head=z2;}
		if (nz>5 && (nz-1-z1)<=head)   {headpref=5;  head=nz-1-z1;}
	}



	/* BIGGER LOOP - HOLE - ALLOWS AUTOMATIC REVERSE PROPAGATION IF 
		HEAD WAVES ARE ENCOUNTERED ON FACES OF EXPANDING CUBE, 
		ALLOWING WAVES TO TRAVEL BACK INTO THE CUBE */
	while ( reverse > -1 )  {

		headw[1]=0; headw[2]=0; headw[3]=0; headw[4]=0;
		headw[5]=0; headw[6]=0;

	/* BIG LOOP */
	while(rad0 && (dx1 || dx2 || dy1 || dy2 || dz1 || dz2))  {

		/* CALCULATE ON PRIMARY (time0) GRID */

		/* TOP SIDE */
      for (igrow=1;igrow<=kminus;igrow++) {  
	if(dz1){
		ii = 0;
		for(j=y1+1; j<=y2-1; j++){
			for(i=x1+1; i<=x2-1; i++){
				sort[ii].time = t0(i,j,z1+1);
				sort[ii].i1 = i;
				sort[ii].i2 = j;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = z1*nxy + X2*nx + X1;
			lasti = (z1+1)*nxy + X2*nx + X1;
			fhead = 0.;
			guess = time0[index];
                        if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,X2,z1+1),
				      t0(X1+1,X2,z1+1),t0(X1+1,X2+1,z1+1),t0(X1,X2+1,z1+1),
				      t0(X1+1,X2,z1  ),t0(X1+1,X2+1,z1  ),t0(X1,X2+1,z1  ),
				      s0(X1,X2,z1), s0(X1,X2,z1+1),
				      s0(X1+1,X2,z1+1),s0(X1+1,X2+1,z1+1),s0(X1,X2+1,z1+1),
				      s0(X1+1,X2,z1  ),s0(X1+1,X2+1,z1  ),s0(X1,X2+1,z1  ));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) {
			  try = fdh3d(              t0(X1,X2,z1+1),
				      t0(X1-1,X2,z1+1),t0(X1-1,X2+1,z1+1),t0(X1,X2+1,z1+1),
				      t0(X1-1,X2,z1  ),t0(X1-1,X2+1,z1  ),t0(X1,X2+1,z1  ),
				      s0(X1,X2,z1), s0(X1,X2,z1+1),
				      s0(X1-1,X2,z1+1),s0(X1-1,X2+1,z1+1),s0(X1,X2+1,z1+1),
				      s0(X1-1,X2,z1  ),s0(X1-1,X2+1,z1  ),s0(X1,X2+1,z1  ));
			  if (try<guess) guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9
			   && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,X2,z1+1),
				      t0(X1+1,X2,z1+1),t0(X1+1,X2-1,z1+1),t0(X1,X2-1,z1+1),
				      t0(X1+1,X2,z1  ),t0(X1+1,X2-1,z1  ),t0(X1,X2-1,z1  ),
				      s0(X1,X2,z1), s0(X1,X2,z1+1),
				      s0(X1+1,X2,z1+1),s0(X1+1,X2-1,z1+1),s0(X1,X2-1,z1+1),
				      s0(X1+1,X2,z1  ),s0(X1+1,X2-1,z1  ),s0(X1,X2-1,z1  ));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9
			   && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(X1,X2,z1+1),
				      t0(X1-1,X2,z1+1),t0(X1-1,X2-1,z1+1),t0(X1,X2-1,z1+1),
				      t0(X1-1,X2,z1  ),t0(X1-1,X2-1,z1  ),t0(X1,X2-1,z1  ),
				      s0(X1,X2,z1), s0(X1,X2,z1+1),
				      s0(X1-1,X2,z1+1),s0(X1-1,X2-1,z1+1),s0(X1,X2-1,z1+1),
				      s0(X1-1,X2,z1  ),s0(X1-1,X2-1,z1  ),s0(X1,X2-1,z1  ));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>y1+1 && X2<y2-1 )  {
			      try = fdhne(t0(X1,X2,z1+1),t0(X1+1,X2,z1+1),t0(X1+1,X2,z1),
					  t0(X1+1,X2-1,z1+1),t0(X1+1,X2+1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1+1,X2,z1+1),s0(X1+1,X2,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 && X2>y1+1 && X2<y2-1 )  {
			      try = fdhne(t0(X1,X2,z1+1),t0(X1-1,X2,z1+1),t0(X1-1,X2,z1),
					  t0(X1-1,X2-1,z1+1),t0(X1-1,X2+1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1-1,X2,z1+1),s0(X1-1,X2,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && X2<ny-1 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,X2,z1+1),t0(X1,X2+1,z1+1),t0(X1,X2+1,z1),
					  t0(X1-1,X2+1,z1+1),t0(X1+1,X2+1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1,X2+1,z1+1),s0(X1,X2+1,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X2>0 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,X2,z1+1),t0(X1,X2-1,z1+1),t0(X1,X2-1,z1),
					  t0(X1-1,X2-1,z1+1),t0(X1+1,X2-1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1,X2-1,z1+1),s0(X1,X2-1,z1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  {
			    try = fdh2d(t0(X1,X2,z1+1),t0(X1+1,X2,z1+1),t0(X1+1,X2,z1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1+1,X2,z1+1),s0(X1+1,X2,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(X1,X2,z1+1),t0(X1-1,X2,z1+1),t0(X1-1,X2,z1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1-1,X2,z1+1),s0(X1-1,X2,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && X2<ny-1 )  {
			    try = fdh2d(t0(X1,X2,z1+1),t0(X1,X2+1,z1+1),t0(X1,X2+1,z1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1,X2+1,z1+1),s0(X1,X2+1,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(X1,X2,z1+1),t0(X1,X2-1,z1+1),t0(X1,X2-1,z1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1,X2-1,z1+1),s0(X1,X2-1,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,X2,z1),t0(X1+1,X2+1,z1),t0(X1,X2+1,z1),
					s0(X1,X2,z1),
					s0(X1+1,X2,z1),s0(X1+1,X2+1,z1),s0(X1,X2+1,z1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9
			     && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,X2,z1),t0(X1+1,X2-1,z1),t0(X1,X2-1,z1),
					s0(X1,X2,z1),
					s0(X1+1,X2,z1),s0(X1+1,X2-1,z1),s0(X1,X2-1,z1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) {
			    try = fdh2d(t0(X1-1,X2,z1),t0(X1-1,X2+1,z1),t0(X1,X2+1,z1),
					s0(X1,X2,z1),
					s0(X1-1,X2,z1),s0(X1-1,X2+1,z1),s0(X1,X2+1,z1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9
			     && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(X1-1,X2,z1),t0(X1-1,X2-1,z1),t0(X1,X2-1,z1),
					s0(X1,X2,z1),
					s0(X1-1,X2,z1),s0(X1-1,X2-1,z1),s0(X1,X2-1,z1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>x1+1 && X1<x2-1 && X2>y1+1 && X2<y2-1 ) {
			    try = fdhnf(t0(X1,X2,z1+1),
					  t0(X1+1,X2,z1+1),t0(X1,X2+1,z1+1),
					  t0(X1-1,X2,z1+1),t0(X1,X2-1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1) );
			    if (try<guess)  guess = try;
			  }
			} 
                          try = t0(X1,X2,z1+1) + .5*(s0(X1,X2,z1)+s0(X1,X2,z1+1));
			  if (try<guess)  guess = try;
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  {
			    try = t0(X1+1,X2,z1) + .5*(s0(X1,X2,z1)+s0(X1+1,X2,z1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-1]<1.e9 && X1>0 )  {
			    try = t0(X1-1,X2,z1) + .5*(s0(X1,X2,z1)+s0(X1-1,X2,z1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nx]<1.e9 && X2<ny-1 )  {
			    try = t0(X1,X2+1,z1) + .5*(s0(X1,X2,z1)+s0(X1,X2+1,z1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nx]<1.e9 && X2>0 )  {
			    try = t0(X1,X2-1,z1) + .5*(s0(X1,X2,z1)+s0(X1,X2-1,z1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index])  {
				time0[index] = guess;
				if (fhead>headtest)  headw[5]++;
			}
		}
		if(z1 == 0) dz1 = 0;
		z1--;
	}
      }
		/* BOTTOM SIDE */
      for (igrow=1;igrow<=kplus;igrow++) {  
	if(dz2){
		ii = 0;
		for(j=y1+1; j<=y2-1; j++){
			for(i=x1+1; i<=x2-1; i++){
				sort[ii].time = t0(i,j,z2-1);
				sort[ii].i1 = i;
				sort[ii].i2 = j;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = z2*nxy + X2*nx + X1;
			lasti = (z2-1)*nxy + X2*nx + X1;
			fhead = 0.;
			guess = time0[index];
                        if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,X2,z2-1),
				      t0(X1+1,X2,z2-1),t0(X1+1,X2+1,z2-1),t0(X1,X2+1,z2-1),
				      t0(X1+1,X2,z2  ),t0(X1+1,X2+1,z2  ),t0(X1,X2+1,z2  ),
				      s0(X1,X2,z2), s0(X1,X2,z2-1),
				      s0(X1+1,X2,z2-1),s0(X1+1,X2+1,z2-1),s0(X1,X2+1,z2-1),
				      s0(X1+1,X2,z2  ),s0(X1+1,X2+1,z2  ),s0(X1,X2+1,z2  ));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) {
			  try = fdh3d(              t0(X1,X2,z2-1),
				      t0(X1-1,X2,z2-1),t0(X1-1,X2+1,z2-1),t0(X1,X2+1,z2-1),
				      t0(X1-1,X2,z2  ),t0(X1-1,X2+1,z2  ),t0(X1,X2+1,z2  ),
				      s0(X1,X2,z2), s0(X1,X2,z2-1),
				      s0(X1-1,X2,z2-1),s0(X1-1,X2+1,z2-1),s0(X1,X2+1,z2-1),
				      s0(X1-1,X2,z2  ),s0(X1-1,X2+1,z2  ),s0(X1,X2+1,z2  ));
			  if (try<guess) guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9
			   && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,X2,z2-1),
				      t0(X1+1,X2,z2-1),t0(X1+1,X2-1,z2-1),t0(X1,X2-1,z2-1),
				      t0(X1+1,X2,z2  ),t0(X1+1,X2-1,z2  ),t0(X1,X2-1,z2  ),
				      s0(X1,X2,z2), s0(X1,X2,z2-1),
				      s0(X1+1,X2,z2-1),s0(X1+1,X2-1,z2-1),s0(X1,X2-1,z2-1),
				      s0(X1+1,X2,z2  ),s0(X1+1,X2-1,z2  ),s0(X1,X2-1,z2  ));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9
			   && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(X1,X2,z2-1),
				      t0(X1-1,X2,z2-1),t0(X1-1,X2-1,z2-1),t0(X1,X2-1,z2-1),
				      t0(X1-1,X2,z2  ),t0(X1-1,X2-1,z2  ),t0(X1,X2-1,z2  ),
				      s0(X1,X2,z2), s0(X1,X2,z2-1),
				      s0(X1-1,X2,z2-1),s0(X1-1,X2-1,z2-1),s0(X1,X2-1,z2-1),
				      s0(X1-1,X2,z2  ),s0(X1-1,X2-1,z2  ),s0(X1,X2-1,z2  ));
			  if (try<guess) guess = try;
			}
                        if(guess > 1.0e9){ 
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>y1+1 && X2<y2-1 )  {
			      try = fdhne(t0(X1,X2,z2-1),t0(X1+1,X2,z2-1),t0(X1+1,X2,z2),
					  t0(X1+1,X2-1,z2-1),t0(X1+1,X2+1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1+1,X2,z2-1),s0(X1+1,X2,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 && X2>y1+1 && X2<y2-1 )  {
			      try = fdhne(t0(X1,X2,z2-1),t0(X1-1,X2,z2-1),t0(X1-1,X2,z2),
					  t0(X1-1,X2-1,z2-1),t0(X1-1,X2+1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1-1,X2,z2-1),s0(X1-1,X2,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && X2<ny-1 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,X2,z2-1),t0(X1,X2+1,z2-1),t0(X1,X2+1,z2),
					  t0(X1-1,X2+1,z2-1),t0(X1+1,X2+1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1,X2+1,z2-1),s0(X1,X2+1,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X2>0 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,X2,z2-1),t0(X1,X2-1,z2-1),t0(X1,X2-1,z2),
					  t0(X1-1,X2-1,z2-1),t0(X1+1,X2-1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1,X2-1,z2-1),s0(X1,X2-1,z2) );
			    if (try<guess)  guess = try;
			  }
		        }
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  {
			    try = fdh2d(t0(X1,X2,z2-1),t0(X1+1,X2,z2-1),t0(X1+1,X2,z2),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1+1,X2,z2-1),s0(X1+1,X2,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(X1,X2,z2-1),t0(X1-1,X2,z2-1),t0(X1-1,X2,z2),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1-1,X2,z2-1),s0(X1-1,X2,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && X2<ny-1 )  {
			    try = fdh2d(t0(X1,X2,z2-1),t0(X1,X2+1,z2-1),t0(X1,X2+1,z2),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1,X2+1,z2-1),s0(X1,X2+1,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(X1,X2,z2-1),t0(X1,X2-1,z2-1),t0(X1,X2-1,z2),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1,X2-1,z2-1),s0(X1,X2-1,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,X2,z2),t0(X1+1,X2+1,z2),t0(X1,X2+1,z2),
					s0(X1,X2,z2),
					s0(X1+1,X2,z2),s0(X1+1,X2+1,z2),s0(X1,X2+1,z2) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9
			     && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,X2,z2),t0(X1+1,X2-1,z2),t0(X1,X2-1,z2),
					s0(X1,X2,z2),
					s0(X1+1,X2,z2),s0(X1+1,X2-1,z2),s0(X1,X2-1,z2) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) {
			    try = fdh2d(t0(X1-1,X2,z2),t0(X1-1,X2+1,z2),t0(X1,X2+1,z2),
					s0(X1,X2,z2),
					s0(X1-1,X2,z2),s0(X1-1,X2+1,z2),s0(X1,X2+1,z2) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9
			     && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(X1-1,X2,z2),t0(X1-1,X2-1,z2),t0(X1,X2-1,z2),
					s0(X1,X2,z2),
					s0(X1-1,X2,z2),s0(X1-1,X2-1,z2),s0(X1,X2-1,z2) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>x1+1 && X1<x2-1 && X2>y1+1 && X2<y2-1 ) {
			    try = fdhnf(t0(X1,X2,z2-1),
					  t0(X1+1,X2,z2-1),t0(X1,X2+1,z2-1),
					  t0(X1-1,X2,z2-1),t0(X1,X2-1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(X1,X2,z2-1) + .5*(s0(X1,X2,z2)+s0(X1,X2,z2-1));
			  if (try<guess)  guess = try;
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  {
			    try = t0(X1+1,X2,z2) + .5*(s0(X1,X2,z2)+s0(X1+1,X2,z2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-1]<1.e9 && X1>0 )  {
			    try = t0(X1-1,X2,z2) + .5*(s0(X1,X2,z2)+s0(X1-1,X2,z2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nx]<1.e9 && X2<ny-1 )  {
			    try = t0(X1,X2+1,z2) + .5*(s0(X1,X2,z2)+s0(X1,X2+1,z2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nx]<1.e9 && X2>0 )  {
			    try = t0(X1,X2-1,z2) + .5*(s0(X1,X2,z2)+s0(X1,X2-1,z2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[6]++;
			}
		}
		if(z2 == nz-1) dz2 = 0;
		z2++;
	}
      }
		/* FRONT SIDE */
      for (igrow=1;igrow<=jminus;igrow++) {  
	if(dy1){
		ii = 0;
		for(k=z1+1; k<=z2-1; k++){
			for(i=x1+1; i<=x2-1; i++){
				sort[ii].time = t0(i,y1+1,k);
				sort[ii].i1 = i;
				sort[ii].i2 = k;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = X2*nxy + y1*nx + X1;
			lasti = X2*nxy + (y1+1)*nx + X1;
			fhead = 0.;
			guess = time0[index];
			if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,y1+1,X2),
				      t0(X1+1,y1+1,X2),t0(X1+1,y1+1,X2+1),t0(X1,y1+1,X2+1),
				      t0(X1+1,y1  ,X2),t0(X1+1,y1  ,X2+1),t0(X1,y1  ,X2+1),
				      s0(X1,y1,X2), s0(X1,y1+1,X2),
				      s0(X1+1,y1+1,X2),s0(X1+1,y1+1,X2+1),s0(X1,y1+1,X2+1),
				      s0(X1+1,y1  ,X2),s0(X1+1,y1  ,X2+1),s0(X1,y1  ,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			  try = fdh3d(              t0(X1,y1+1,X2),
				      t0(X1-1,y1+1,X2),t0(X1-1,y1+1,X2+1),t0(X1,y1+1,X2+1),
				      t0(X1-1,y1  ,X2),t0(X1-1,y1  ,X2+1),t0(X1,y1  ,X2+1),
				      s0(X1,y1,X2), s0(X1,y1+1,X2),
				      s0(X1-1,y1+1,X2),s0(X1-1,y1+1,X2+1),s0(X1,y1+1,X2+1),
				      s0(X1-1,y1  ,X2),s0(X1-1,y1  ,X2+1),s0(X1,y1  ,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,y1+1,X2),
				      t0(X1+1,y1+1,X2),t0(X1+1,y1+1,X2-1),t0(X1,y1+1,X2-1),
				      t0(X1+1,y1  ,X2),t0(X1+1,y1  ,X2-1),t0(X1,y1  ,X2-1),
				      s0(X1,y1,X2), s0(X1,y1+1,X2),
				      s0(X1+1,y1+1,X2),s0(X1+1,y1+1,X2-1),s0(X1,y1+1,X2-1),
				      s0(X1+1,y1  ,X2),s0(X1+1,y1  ,X2-1),s0(X1,y1  ,X2-1));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(X1,y1+1,X2),
				      t0(X1-1,y1+1,X2),t0(X1-1,y1+1,X2-1),t0(X1,y1+1,X2-1),
				      t0(X1-1,y1  ,X2),t0(X1-1,y1  ,X2-1),t0(X1,y1  ,X2-1),
				      s0(X1,y1,X2), s0(X1,y1+1,X2),
				      s0(X1-1,y1+1,X2),s0(X1-1,y1+1,X2-1),s0(X1,y1+1,X2-1),
				      s0(X1-1,y1  ,X2),s0(X1-1,y1  ,X2-1),s0(X1,y1  ,X2-1));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(X1,y1+1,X2),t0(X1+1,y1+1,X2),t0(X1+1,y1,X2),
					  t0(X1+1,y1+1,X2-1),t0(X1+1,y1+1,X2+1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1+1,y1+1,X2),s0(X1+1,y1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(X1,y1+1,X2),t0(X1-1,y1+1,X2),t0(X1-1,y1,X2),
					  t0(X1-1,y1+1,X2-1),t0(X1-1,y1+1,X2+1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1-1,y1+1,X2),s0(X1-1,y1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,y1+1,X2),t0(X1,y1+1,X2+1),t0(X1,y1,X2+1),
					  t0(X1-1,y1+1,X2+1),t0(X1+1,y1+1,X2+1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1,y1+1,X2+1),s0(X1,y1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,y1+1,X2),t0(X1,y1+1,X2-1),t0(X1,y1,X2-1),
					  t0(X1-1,y1+1,X2-1),t0(X1+1,y1+1,X2-1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1,y1+1,X2-1),s0(X1,y1,X2-1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  {
			    try = fdh2d(t0(X1,y1+1,X2),t0(X1+1,y1+1,X2),t0(X1+1,y1,X2),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1+1,y1+1,X2),s0(X1+1,y1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(X1,y1+1,X2),t0(X1-1,y1+1,X2),t0(X1-1,y1,X2),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1-1,y1+1,X2),s0(X1-1,y1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {
			    try = fdh2d(t0(X1,y1+1,X2),t0(X1,y1+1,X2+1),t0(X1,y1,X2+1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1,y1+1,X2+1),s0(X1,y1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(X1,y1+1,X2),t0(X1,y1+1,X2-1),t0(X1,y1,X2-1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1,y1+1,X2-1),s0(X1,y1,X2-1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,y1,X2),t0(X1+1,y1,X2+1),t0(X1,y1,X2+1),
					s0(X1,y1,X2),
					s0(X1+1,y1,X2),s0(X1+1,y1,X2+1),s0(X1,y1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,y1,X2),t0(X1+1,y1,X2-1),t0(X1,y1,X2-1),
					s0(X1,y1,X2),
					s0(X1+1,y1,X2),s0(X1+1,y1,X2-1),s0(X1,y1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			    try = fdh2d(t0(X1-1,y1,X2),t0(X1-1,y1,X2+1),t0(X1,y1,X2+1),
					s0(X1,y1,X2),
					s0(X1-1,y1,X2),s0(X1-1,y1,X2+1),s0(X1,y1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(X1-1,y1,X2),t0(X1-1,y1,X2-1),t0(X1,y1,X2-1),
					s0(X1,y1,X2),
					s0(X1-1,y1,X2),s0(X1-1,y1,X2-1),s0(X1,y1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>x1+1 && X1<x2-1 && X2>z1+1 && X2<z2-1 ) {
			    try = fdhnf(t0(X1,y1+1,X2),
					  t0(X1+1,y1+1,X2),t0(X1,y1+1,X2+1),
					  t0(X1-1,y1+1,X2),t0(X1,y1+1,X2-1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(X1,y1+1,X2) + .5*(s0(X1,y1,X2)+s0(X1,y1+1,X2));
			  if (try<guess)  guess = try;
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  {
			    try = t0(X1+1,y1,X2) + .5*(s0(X1,y1,X2)+s0(X1+1,y1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-1]<1.e9 && X1>0 )  {
			    try = t0(X1-1,y1,X2) + .5*(s0(X1,y1,X2)+s0(X1-1,y1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  {
			    try = t0(X1,y1,X2+1) + .5*(s0(X1,y1,X2)+s0(X1,y1,X2+1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nxy]<1.e9 && X2>0 )  {
			    try = t0(X1,y1,X2-1) + .5*(s0(X1,y1,X2)+s0(X1,y1,X2-1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[3]++;
			}
		}
		if(y1 == 0) dy1 = 0;
		y1--;
	}
      }
		/* BACK SIDE */
      for (igrow=1;igrow<=jplus;igrow++) {  
	if(dy2){
		ii = 0;
		for(k=z1+1; k<=z2-1; k++){
			for(i=x1+1; i<=x2-1; i++){
				sort[ii].time = t0(i,y2-1,k);
				sort[ii].i1 = i;
				sort[ii].i2 = k;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = X2*nxy + y2*nx + X1;
			lasti = X2*nxy + (y2-1)*nx + X1;
			fhead = 0.;
			guess = time0[index];
			if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,y2-1,X2),
				      t0(X1+1,y2-1,X2),t0(X1+1,y2-1,X2+1),t0(X1,y2-1,X2+1),
				      t0(X1+1,y2  ,X2),t0(X1+1,y2  ,X2+1),t0(X1,y2  ,X2+1),
				      s0(X1,y2,X2), s0(X1,y2-1,X2),
				      s0(X1+1,y2-1,X2),s0(X1+1,y2-1,X2+1),s0(X1,y2-1,X2+1),
				      s0(X1+1,y2  ,X2),s0(X1+1,y2  ,X2+1),s0(X1,y2  ,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			  try = fdh3d(              t0(X1,y2-1,X2),
				      t0(X1-1,y2-1,X2),t0(X1-1,y2-1,X2+1),t0(X1,y2-1,X2+1),
				      t0(X1-1,y2  ,X2),t0(X1-1,y2  ,X2+1),t0(X1,y2  ,X2+1),
				      s0(X1,y2,X2), s0(X1,y2-1,X2),
				      s0(X1-1,y2-1,X2),s0(X1-1,y2-1,X2+1),s0(X1,y2-1,X2+1),
				      s0(X1-1,y2  ,X2),s0(X1-1,y2  ,X2+1),s0(X1,y2  ,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,y2-1,X2),
				      t0(X1+1,y2-1,X2),t0(X1+1,y2-1,X2-1),t0(X1,y2-1,X2-1),
				      t0(X1+1,y2  ,X2),t0(X1+1,y2  ,X2-1),t0(X1,y2  ,X2-1),
				      s0(X1,y2,X2), s0(X1,y2-1,X2),
				      s0(X1+1,y2-1,X2),s0(X1+1,y2-1,X2-1),s0(X1,y2-1,X2-1),
				      s0(X1+1,y2  ,X2),s0(X1+1,y2  ,X2-1),s0(X1,y2  ,X2-1));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(X1,y2-1,X2),
				      t0(X1-1,y2-1,X2),t0(X1-1,y2-1,X2-1),t0(X1,y2-1,X2-1),
				      t0(X1-1,y2  ,X2),t0(X1-1,y2  ,X2-1),t0(X1,y2  ,X2-1),
				      s0(X1,y2,X2), s0(X1,y2-1,X2),
				      s0(X1-1,y2-1,X2),s0(X1-1,y2-1,X2-1),s0(X1,y2-1,X2-1),
				      s0(X1-1,y2  ,X2),s0(X1-1,y2  ,X2-1),s0(X1,y2  ,X2-1));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(X1,y2-1,X2),t0(X1+1,y2-1,X2),t0(X1+1,y2,X2),
					  t0(X1+1,y2-1,X2-1),t0(X1+1,y2-1,X2+1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1+1,y2-1,X2),s0(X1+1,y2,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(X1,y2-1,X2),t0(X1-1,y2-1,X2),t0(X1-1,y2,X2),
					  t0(X1-1,y2-1,X2-1),t0(X1-1,y2-1,X2+1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1-1,y2-1,X2),s0(X1-1,y2,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,y2-1,X2),t0(X1,y2-1,X2+1),t0(X1,y2,X2+1),
					  t0(X1-1,y2-1,X2+1),t0(X1+1,y2-1,X2+1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1,y2-1,X2+1),s0(X1,y2,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,y2-1,X2),t0(X1,y2-1,X2-1),t0(X1,y2,X2-1),
					  t0(X1-1,y2-1,X2-1),t0(X1+1,y2-1,X2-1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1,y2-1,X2-1),s0(X1,y2,X2-1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  {
			    try = fdh2d(t0(X1,y2-1,X2),t0(X1+1,y2-1,X2),t0(X1+1,y2,X2),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1+1,y2-1,X2),s0(X1+1,y2,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(X1,y2-1,X2),t0(X1-1,y2-1,X2),t0(X1-1,y2,X2),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1-1,y2-1,X2),s0(X1-1,y2,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {
			    try = fdh2d(t0(X1,y2-1,X2),t0(X1,y2-1,X2+1),t0(X1,y2,X2+1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1,y2-1,X2+1),s0(X1,y2,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(X1,y2-1,X2),t0(X1,y2-1,X2-1),t0(X1,y2,X2-1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1,y2-1,X2-1),s0(X1,y2,X2-1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,y2,X2),t0(X1+1,y2,X2+1),t0(X1,y2,X2+1),
					s0(X1,y2,X2),
					s0(X1+1,y2,X2),s0(X1+1,y2,X2+1),s0(X1,y2,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,y2,X2),t0(X1+1,y2,X2-1),t0(X1,y2,X2-1),
					s0(X1,y2,X2),
					s0(X1+1,y2,X2),s0(X1+1,y2,X2-1),s0(X1,y2,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			    try = fdh2d(t0(X1-1,y2,X2),t0(X1-1,y2,X2+1),t0(X1,y2,X2+1),
					s0(X1,y2,X2),
					s0(X1-1,y2,X2),s0(X1-1,y2,X2+1),s0(X1,y2,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(X1-1,y2,X2),t0(X1-1,y2,X2-1),t0(X1,y2,X2-1),
					s0(X1,y2,X2),
					s0(X1-1,y2,X2),s0(X1-1,y2,X2-1),s0(X1,y2,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>x1+1 && X1<x2-1 && X2>z1+1 && X2<z2-1 ) {
			    try = fdhnf(t0(X1,y2-1,X2),
					  t0(X1+1,y2-1,X2),t0(X1,y2-1,X2+1),
					  t0(X1-1,y2-1,X2),t0(X1,y2-1,X2-1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(X1,y2-1,X2) + .5*(s0(X1,y2,X2)+s0(X1,y2-1,X2));
			  if (try<guess)  guess = try;
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  {
			    try = t0(X1+1,y2,X2) + .5*(s0(X1,y2,X2)+s0(X1+1,y2,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-1]<1.e9 && X1>0 )  {
			    try = t0(X1-1,y2,X2) + .5*(s0(X1,y2,X2)+s0(X1-1,y2,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  {
			    try = t0(X1,y2,X2+1) + .5*(s0(X1,y2,X2)+s0(X1,y2,X2+1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nxy]<1.e9 && X2>0 )  {
			    try = t0(X1,y2,X2-1) + .5*(s0(X1,y2,X2)+s0(X1,y2,X2-1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[4]++;
			}
		}
		if(y2 == ny-1) dy2 = 0;
		y2++;
	}
      }
		/* LEFT SIDE */
      for (igrow=1;igrow<=iminus;igrow++) {  
	if(dx1){
		ii = 0;
		for(k=z1+1; k<=z2-1; k++){
			for(j=y1+1; j<=y2-1; j++){
				sort[ii].time = t0(x1+1,j,k);
				sort[ii].i1 = j;
				sort[ii].i2 = k;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = X2*nxy + X1*nx + x1;
			lasti = X2*nxy + X1*nx + (x1+1);
			fhead = 0.;
			guess = time0[index];
			if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) {
			  try = fdh3d(              t0(x1+1,X1,X2),
				      t0(x1+1,X1+1,X2),t0(x1+1,X1+1,X2+1),t0(x1+1,X1,X2+1),
				      t0(x1  ,X1+1,X2),t0(x1  ,X1+1,X2+1),t0(x1  ,X1,X2+1),
				      s0(x1,X1,X2), s0(x1+1,X1,X2),
				      s0(x1+1,X1+1,X2),s0(x1+1,X1+1,X2+1),s0(x1+1,X1,X2+1),
				      s0(x1  ,X1+1,X2),s0(x1  ,X1+1,X2+1),s0(x1  ,X1,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			  try = fdh3d(              t0(x1+1,X1,X2),
				      t0(x1+1,X1-1,X2),t0(x1+1,X1-1,X2+1),t0(x1+1,X1,X2+1),
				      t0(x1  ,X1-1,X2),t0(x1  ,X1-1,X2+1),t0(x1  ,X1,X2+1),
				      s0(x1,X1,X2), s0(x1+1,X1,X2),
				      s0(x1+1,X1-1,X2),s0(x1+1,X1-1,X2+1),s0(x1+1,X1,X2+1),
				      s0(x1  ,X1-1,X2),s0(x1  ,X1-1,X2+1),s0(x1  ,X1,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) {
			  try = fdh3d(              t0(x1+1,X1,X2),
				      t0(x1+1,X1+1,X2),t0(x1+1,X1+1,X2-1),t0(x1+1,X1,X2-1),
				      t0(x1  ,X1+1,X2),t0(x1  ,X1+1,X2-1),t0(x1  ,X1,X2-1),
				      s0(x1,X1,X2), s0(x1+1,X1,X2),
				      s0(x1+1,X1+1,X2),s0(x1+1,X1+1,X2-1),s0(x1+1,X1,X2-1),
				      s0(x1  ,X1+1,X2),s0(x1  ,X1+1,X2-1),s0(x1  ,X1,X2-1));
			  if (try<guess) guess = try;
			}
			if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(x1+1,X1,X2),
				      t0(x1+1,X1-1,X2),t0(x1+1,X1-1,X2-1),t0(x1+1,X1,X2-1),
				      t0(x1  ,X1-1,X2),t0(x1  ,X1-1,X2-1),t0(x1  ,X1,X2-1),
				      s0(x1,X1,X2), s0(x1+1,X1,X2),
				      s0(x1+1,X1-1,X2),s0(x1+1,X1-1,X2-1),s0(x1+1,X1,X2-1),
				      s0(x1  ,X1-1,X2),s0(x1  ,X1-1,X2-1),s0(x1  ,X1,X2-1));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+nx] < 1.e9 && X1<ny-1 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(x1+1,X1,X2),t0(x1+1,X1+1,X2),t0(x1,X1+1,X2),
					  t0(x1+1,X1+1,X2-1),t0(x1+1,X1+1,X2+1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1+1,X2),s0(x1,X1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X1>0 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(x1+1,X1,X2),t0(x1+1,X1-1,X2),t0(x1,X1-1,X2),
					  t0(x1+1,X1-1,X2-1),t0(x1+1,X1-1,X2+1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1-1,X2),s0(x1,X1-1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>y1+1 && X1<y2-1 )  {
			      try = fdhne(t0(x1+1,X1,X2),t0(x1+1,X1,X2+1),t0(x1,X1,X2+1),
					  t0(x1+1,X1-1,X2+1),t0(x1+1,X1+1,X2+1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1,X2+1),s0(x1,X1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>y1+1 && X1<y2-1 )  {
			      try = fdhne(t0(x1+1,X1,X2),t0(x1+1,X1,X2-1),t0(x1,X1,X2-1),
					  t0(x1+1,X1-1,X2-1),t0(x1+1,X1+1,X2-1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1,X2-1),s0(x1,X1,X2-1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+nx] < 1.e9 && X1<ny-1 )  {
			    try = fdh2d(t0(x1+1,X1,X2),t0(x1+1,X1+1,X2),t0(x1,X1+1,X2),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1+1,X2),s0(x1,X1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(x1+1,X1,X2),t0(x1+1,X1-1,X2),t0(x1,X1-1,X2),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1-1,X2),s0(x1,X1-1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {
			    try = fdh2d(t0(x1+1,X1,X2),t0(x1+1,X1,X2+1),t0(x1,X1,X2+1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1,X2+1),s0(x1,X1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(x1+1,X1,X2),t0(x1+1,X1,X2-1),t0(x1,X1,X2-1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1,X2-1),s0(x1,X1,X2-1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) {
			    try = fdh2d(t0(x1,X1+1,X2),t0(x1,X1+1,X2+1),t0(x1,X1,X2+1),
					s0(x1,X1,X2),
					s0(x1,X1+1,X2),s0(x1,X1+1,X2+1),s0(x1,X1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) {
			    try = fdh2d(t0(x1,X1+1,X2),t0(x1,X1+1,X2-1),t0(x1,X1,X2-1),
					s0(x1,X1,X2),
					s0(x1,X1+1,X2),s0(x1,X1+1,X2-1),s0(x1,X1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			    try = fdh2d(t0(x1,X1-1,X2),t0(x1,X1-1,X2+1),t0(x1,X1,X2+1),
					s0(x1,X1,X2),
					s0(x1,X1-1,X2),s0(x1,X1-1,X2+1),s0(x1,X1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(x1,X1-1,X2),t0(x1,X1-1,X2-1),t0(x1,X1,X2-1),
					s0(x1,X1,X2),
					s0(x1,X1-1,X2),s0(x1,X1-1,X2-1),s0(x1,X1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>y1+1 && X1<y2-1 && X2>z1+1 && X2<z2-1 ) {
			    try = fdhnf(t0(x1+1,X1,X2),
					  t0(x1+1,X1+1,X2),t0(x1+1,X1,X2+1),
					  t0(x1+1,X1-1,X2),t0(x1+1,X1,X2-1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(x1+1,X1,X2) + .5*(s0(x1,X1,X2)+s0(x1+1,X1,X2));
			  if (try<guess)  guess = try;
                          if ( time0[index+nx]<1.e9 && X1<ny-1 )  {
			    try = t0(x1,X1+1,X2) + .5*(s0(x1,X1,X2)+s0(x1,X1+1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nx]<1.e9 && X1>0 )  {
			    try = t0(x1,X1-1,X2) + .5*(s0(x1,X1,X2)+s0(x1,X1-1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  {
			    try = t0(x1,X1,X2+1) + .5*(s0(x1,X1,X2)+s0(x1,X1,X2+1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nxy]<1.e9 && X2>0 )  {
			    try = t0(x1,X1,X2-1) + .5*(s0(x1,X1,X2)+s0(x1,X1,X2-1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[1]++;
			}
		}
		if(x1 == 0) dx1 = 0;
		x1--;
	}
      }
		/* RIGHT SIDE */
      for (igrow=1;igrow<=iplus;igrow++) {  
	if(dx2){
		ii = 0;
		for(k=z1+1; k<=z2-1; k++){
			for(j=y1+1; j<=y2-1; j++){
				sort[ii].time = t0(x2-1,j,k);
				sort[ii].i1 = j;
				sort[ii].i2 = k;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = X2*nxy + X1*nx + x2;
			lasti = X2*nxy + X1*nx + (x2-1);
			fhead = 0.;
			guess = time0[index];
			if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) {
			  try = fdh3d(              t0(x2-1,X1,X2),
				      t0(x2-1,X1+1,X2),t0(x2-1,X1+1,X2+1),t0(x2-1,X1,X2+1),
				      t0(x2  ,X1+1,X2),t0(x2  ,X1+1,X2+1),t0(x2  ,X1,X2+1),
				      s0(x2,X1,X2), s0(x2-1,X1,X2),
				      s0(x2-1,X1+1,X2),s0(x2-1,X1+1,X2+1),s0(x2-1,X1,X2+1),
				      s0(x2  ,X1+1,X2),s0(x2  ,X1+1,X2+1),s0(x2  ,X1,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			  try = fdh3d(              t0(x2-1,X1,X2),
				      t0(x2-1,X1-1,X2),t0(x2-1,X1-1,X2+1),t0(x2-1,X1,X2+1),
				      t0(x2  ,X1-1,X2),t0(x2  ,X1-1,X2+1),t0(x2  ,X1,X2+1),
				      s0(x2,X1,X2), s0(x2-1,X1,X2),
				      s0(x2-1,X1-1,X2),s0(x2-1,X1-1,X2+1),s0(x2-1,X1,X2+1),
				      s0(x2  ,X1-1,X2),s0(x2  ,X1-1,X2+1),s0(x2  ,X1,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) {
			  try = fdh3d(              t0(x2-1,X1,X2),
				      t0(x2-1,X1+1,X2),t0(x2-1,X1+1,X2-1),t0(x2-1,X1,X2-1),
				      t0(x2  ,X1+1,X2),t0(x2  ,X1+1,X2-1),t0(x2  ,X1,X2-1),
				      s0(x2,X1,X2), s0(x2-1,X1,X2),
				      s0(x2-1,X1+1,X2),s0(x2-1,X1+1,X2-1),s0(x2-1,X1,X2-1),
				      s0(x2  ,X1+1,X2),s0(x2  ,X1+1,X2-1),s0(x2  ,X1,X2-1));
			  if (try<guess) guess = try;
			}
			if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(x2-1,X1,X2),
				      t0(x2-1,X1-1,X2),t0(x2-1,X1-1,X2-1),t0(x2-1,X1,X2-1),
				      t0(x2  ,X1-1,X2),t0(x2  ,X1-1,X2-1),t0(x2  ,X1,X2-1),
				      s0(x2,X1,X2), s0(x2-1,X1,X2),
				      s0(x2-1,X1-1,X2),s0(x2-1,X1-1,X2-1),s0(x2-1,X1,X2-1),
				      s0(x2  ,X1-1,X2),s0(x2  ,X1-1,X2-1),s0(x2  ,X1,X2-1));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+nx] < 1.e9 && X1<ny-1 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(x2-1,X1,X2),t0(x2-1,X1+1,X2),t0(x2,X1+1,X2),
					  t0(x2-1,X1+1,X2-1),t0(x2-1,X1+1,X2+1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1+1,X2),s0(x2,X1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X1>0 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(x2-1,X1,X2),t0(x2-1,X1-1,X2),t0(x2,X1-1,X2),
					  t0(x2-1,X1-1,X2-1),t0(x2-1,X1-1,X2+1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1-1,X2),s0(x2,X1-1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>y1+1 && X1<y2-1 )  {
			      try = fdhne(t0(x2-1,X1,X2),t0(x2-1,X1,X2+1),t0(x2,X1,X2+1),
					  t0(x2-1,X1-1,X2+1),t0(x2-1,X1+1,X2+1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1,X2+1),s0(x2,X1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>y1+1 && X1<y2-1 )  {
			      try = fdhne(t0(x2-1,X1,X2),t0(x2-1,X1,X2-1),t0(x2,X1,X2-1),
					  t0(x2-1,X1-1,X2-1),t0(x2-1,X1+1,X2-1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1,X2-1),s0(x2,X1,X2-1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+nx] < 1.e9 && X1<ny-1 )  {
			    try = fdh2d(t0(x2-1,X1,X2),t0(x2-1,X1+1,X2),t0(x2,X1+1,X2),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1+1,X2),s0(x2,X1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(x2-1,X1,X2),t0(x2-1,X1-1,X2),t0(x2,X1-1,X2),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1-1,X2),s0(x2,X1-1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {
			    try = fdh2d(t0(x2-1,X1,X2),t0(x2-1,X1,X2+1),t0(x2,X1,X2+1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1,X2+1),s0(x2,X1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(x2-1,X1,X2),t0(x2-1,X1,X2-1),t0(x2,X1,X2-1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1,X2-1),s0(x2,X1,X2-1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) {
			    try = fdh2d(t0(x2,X1+1,X2),t0(x2,X1+1,X2+1),t0(x2,X1,X2+1),
					s0(x2,X1,X2),
					s0(x2,X1+1,X2),s0(x2,X1+1,X2+1),s0(x2,X1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) {
			    try = fdh2d(t0(x2,X1+1,X2),t0(x2,X1+1,X2-1),t0(x2,X1,X2-1),
					s0(x2,X1,X2),
					s0(x2,X1+1,X2),s0(x2,X1+1,X2-1),s0(x2,X1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			    try = fdh2d(t0(x2,X1-1,X2),t0(x2,X1-1,X2+1),t0(x2,X1,X2+1),
					s0(x2,X1,X2),
					s0(x2,X1-1,X2),s0(x2,X1-1,X2+1),s0(x2,X1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(x2,X1-1,X2),t0(x2,X1-1,X2-1),t0(x2,X1,X2-1),
					s0(x2,X1,X2),
					s0(x2,X1-1,X2),s0(x2,X1-1,X2-1),s0(x2,X1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>y1+1 && X1<y2-1 && X2>z1+1 && X2<z2-1 ) {
			    try = fdhnf(t0(x2-1,X1,X2),
					  t0(x2-1,X1+1,X2),t0(x2-1,X1,X2+1),
					  t0(x2-1,X1-1,X2),t0(x2-1,X1,X2-1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(x2-1,X1,X2) + .5*(s0(x2,X1,X2)+s0(x2-1,X1,X2));
			  if (try<guess)  guess = try;
                          if ( time0[index+nx]<1.e9 && X1<ny-1 )  {
			    try = t0(x2,X1+1,X2) + .5*(s0(x2,X1,X2)+s0(x2,X1+1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nx]<1.e9 && X1>0 )  {
			    try = t0(x2,X1-1,X2) + .5*(s0(x2,X1,X2)+s0(x2,X1-1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  {
			    try = t0(x2,X1,X2+1) + .5*(s0(x2,X1,X2)+s0(x2,X1,X2+1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nxy]<1.e9 && X2>0 )  {
			    try = t0(x2,X1,X2-1) + .5*(s0(x2,X1,X2)+s0(x2,X1,X2-1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[2]++;
			}
		}
		if(x2 == nx-1) dx2 = 0;
		x2++;
	}
      }

		/* UPDATE RADIUS */
		radius++;
		if(radius%10 == 0) fprintf(stderr,"Completed radius = %d\n",radius);
                if(radius == maxrad) rad0 = 0;

	}	/* END BIG LOOP */


	/* TEST IF REVERSE PROPAGATION IS NEEDED */

	if (headw[1]==0 && headw[2]==0 && headw[3]==0 && headw[4]==0 
		     && headw[5]==0 && headw[6]==0)
		reverse=0;
	else {
		head=0;
		if (headw[1]>0) {
			fprintf(stderr,"Head waves found on left: %d\n",headw[1]);
			if (headw[1]>head)  {
				head = headw[1];
				srcwall = 1;
			}
		}
		if (headw[2]>0) {
			fprintf(stderr,"Head waves found on right: %d\n",headw[2]);
			if (headw[2]>head)  {
				head = headw[2];
				srcwall = 2;
			}
		}
		if (headw[3]>0) {
			fprintf(stderr,"Head waves found on front: %d\n",headw[3]);
			if (headw[3]>head)  {
				head = headw[3];
				srcwall = 3;
			}
		}
		if (headw[4]>0) {
			fprintf(stderr,"Head waves found on back: %d\n",headw[4]);
			if (headw[4]>head)  {
				head = headw[4];
				srcwall = 4;
			}
		}
		if (headw[5]>0) {
			fprintf(stderr,"Head waves found on top: %d\n",headw[5]);
			if (headw[5]>head)  {
				head = headw[5];
				srcwall = 5;
			}
		}
		if (headw[6]>0) {
			fprintf(stderr,"Head waves found on bottom: %d\n",headw[6]);
			if (headw[6]>head)  {
				head = headw[6];
				srcwall = 6;
			}
		}
		if (headpref>0 && headw[headpref]>0) {
			fprintf(stderr,"Preference to restart on wall opposite source\n");
			srcwall = headpref;
		}
		/* SET LOCATIONS OF SIDES OF THE CUBE SO THAT CUBE IS A FACE */
		dx1=1; dx2=1; dy1=1; dy2=1; dz1=1; dz2=1; rad0=1;
		radius = 1;
		if (srcwall == 1)	{  x2=1;
			fprintf(stderr,"RESTART at left side of model\n");  }
		else	{  x2=nx;	dx2=0;  }
		if (srcwall == 2)	{ x1=nx-2;
			fprintf(stderr,"RESTART at right side of model\n");  }
		else	{  x1= -1;	dx1=0;  }
		if (srcwall == 3)	{ y2=1;
			fprintf(stderr,"RESTART at front side of model\n");  }
		else	{  y2=ny;	dy2=0;  }
		if (srcwall == 4)	{ y1=ny-2;
			fprintf(stderr,"RESTART at back side of model\n");  }
		else	{  y1= -1;	dy1=0;  }
		if (srcwall == 5)	{ z2=1;
			fprintf(stderr,"RESTART at top side of model\n");  }
		else	{  z2=nz;	dz2=0;  }
		if (srcwall == 6)	{ z1=nz-2;
			fprintf(stderr,"RESTART at bottom side of model\n");  }
		else	{  z1= -1;	dz1=0;  }
		if (reverse == 0)  
			fprintf(stderr,"WARNING:  RESTART CANCELLED by choice of input parameter 'reverse'\n");
	}
	reverse--;

	}	/* END BIGGER LOOP - HOLE */


    /* exctract only traveltime for surficial layer */
    for (j=0;j<ny;j++) {
       for (i=0;i<nx;i++) {
        
           times[nx*j + i] = time0[nxy*(0) + nx*(j) + (i)];
          
       }
    }

    
	/* OUTPUT COMPLETED WAVEFRONT */
	/*write(tfint,time0,nxyz*4); */
	
         //write(tfint,times,nxy*4); 
        //fwrite(times,sizeof(float),nxy*4,tfint);
        fwrite(times,sizeof(float),nxy,tfint);
	
	fclose(vfint);
	fclose(tfint);
	
	
	fprintf(stderr,"wavefront done \n");
}

/* -------------------------------------------------------------------------- */

compar(a,b)
struct sorted *a, *b;
{
	if(a->time > b->time) return(1);
	if(b->time > a->time) return(-1);
	else return(0);
}

/* 3D TRANSMISSION STENCIL
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE
   JAH 11/91 */
float fdh3d(t1,t2,t3,t4,t5,t6,t7,ss0,s1,s2,s3,s4,s5,s6,s7)
     float  t1,t2,t3,t4,t5,t6,t7,ss0,s1,s2,s3,s4,s5,s6,s7;
     /* ss0 at newpoint; s1,t1 adjacent on oldface;
	s2,t2 and s4,t4 on oldface adjacent to s1;
	s3,t3 on oldface diametrically opposite newpoint;
	s5,t5 on newface adjacent to newpoint AND to s2;
	s6,t6 on newface diagonal to newpoint (adjacent to s3);
	s7,t7 on newface adjacent to newpoint AND to s4
	*/
{
  float x,slo;
  double sqrt();
  slo = .125*(ss0+s1+s2+s3+s4+s5+s6+s7);
  x = 6.*slo*slo - (t4-t2)*(t4-t2) - (t2-t6)*(t2-t6) - (t6-t4)*(t6-t4)
                 - (t7-t5)*(t7-t5) - (t5-t1)*(t5-t1) - (t1-t7)*(t1-t7);
  if (x>=0.)  {
    x = t3 + sqrt(.5*x);
    if ( (x<t1) || (x<t2) || (x<t4) || (x<t5) || (x<t6) || (x<t7) )  
      x = 1.e11;   /* ACAUSAL; ABORT */
  }
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */
  return(x);
}

/* 3D STENCIL FOR NEW EDGE
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE
   JAH 11/91 */
float fdhne(t1,t2,t3,t4,t5,ss0,s1,s2,s3)
     float  t1,t2,t3,t4,t5,ss0,s1,s2,s3;
     /* ss0 at newpoint; s1,t1 adjacent on oldface;
	s2,t2 diagonal on oldface; s3,t3 adjacent on newface;
	t4,t5 beside t2 on old face opposite each other */
{
  float x,slo;
  double sqrt();
  slo = .25*(ss0+s1+s2+s3);
  x = 2.*slo*slo - (t3-t1)*(t3-t1) - .5*(t5-t4)*(t5-t4);
  if (x>=0.)  {
    x = t2 + sqrt(x);
    if ( (x<t1) || (x<t3) || (x<t4) || (x<t5) )     /* ACAUSAL; ABORT */
      x = 1.e11;
  }
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */
  return(x);
}

/* 2D TRANSMISSION STENCIL (FOR HEAD WAVES ON FACES OF GRID CELLS)
   STENCIL FROM VIDALE (1988 2D PAPER); CONDITIONS AND OTHER OPTIONS FROM HOLE
   JAH 11/91 */
float fdh2d(t1,t2,t3,ss0,s1,s2,s3)
     float  t1,t2,t3,ss0,s1,s2,s3;
     /* ss0 at newpoint; s1,t1 & s3,t3 adjacent; s2,t2 diagonal
      */
{
  float x,slo;
  double sqrt();
  slo = .25*(ss0+s1+s2+s3);
  x = 2.*slo*slo - (t3-t1)*(t3-t1);
  if (x>=0.)  {
    x = t2 + sqrt(x);
    if ( (x<t1) || (x<t3) )  x = 1.e11;   /* ACAUSAL; ABORT */
  }
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */
  return(x);
}

/* 3D STENCIL FOR NEW FACE
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE
   JAH 11/91 */
float fdhnf(t1,t2,t3,t4,t5,ss0,s1)
     float  t1,t2,t3,t4,t5,ss0,s1;
     /* ss0 at newpoint; s1,t1 adjacent on old face;
	t2,t4 beside t1 on old face and opposite each other;
	t3,t5 beside t1 on old face and opposite each other
	*/
{
  float x,slo;
  double sqrt();
  slo = .5*(ss0+s1);
  x = slo*slo - .25*( (t4-t2)*(t4-t2) + (t5-t3)*(t5-t3) );
  if (x>=0.)  {
    x = t1 + sqrt(x);
    if ( (x<t2) || (x<t3) || (x<t4) || (x<t5) )     /* ACAUSAL; ABORT */
      x = 1.e11;
  }
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */
  return(x);
}
