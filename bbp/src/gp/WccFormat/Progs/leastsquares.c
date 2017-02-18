#include "include.h"
#include "structure.h"
#include "function.h"

#define NPBLOCK 10000
#define NPMAX 10000000

void lsqr(float *x,float *y,float *w,int np,float *c0,int order);

main(int ac,char **av)
{
FILE *fpr, *fpw, *fopfile();
float *xp, *yp, *wp, *c0, x0, y0;
int i, nterms, nx, ix, np, ip;
int getxmin, getxmax;
float dx, xb;
float xmin = 1.0e+16;
float xmax = -1.0e+16;

char infile[128], str[512];
char outfile[128];

nx = 1000;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
mstpar("nterms","d",&nterms);
getpar("xmin","f",&xmin);
getpar("xmax","f",&xmax);
getpar("nx","d",&nx);
endpar();

getxmin = 1;
if(xmin < 1.0e+16)
   getxmin = 0;

getxmax = 1;
if(xmax > -1.0e+16)
   getxmax = 0;

np = NPBLOCK;

c0 = (float *) check_malloc (nterms*sizeof(float));
xp = (float *) check_malloc (np*sizeof(float));
yp = (float *) check_malloc (np*sizeof(float));
wp = (float *) check_malloc (np*sizeof(float));

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

ip = 0;
while(fgets(str,512,fpr) != NULL)
   {
   if(ip+1 == np)
      {
      np = np + NPBLOCK;
      xp = (float *) check_realloc (xp,np*sizeof(float));
      yp = (float *) check_realloc (yp,np*sizeof(float));
      wp = (float *) check_realloc (wp,np*sizeof(float));
      }

   if(sscanf(str,"%f %f %f",&xp[ip],&yp[ip],&wp[ip])==2)
      wp[ip] = 1.0;

   if(xp[ip] < xmin && getxmin)
      xmin = xp[ip];
   if(xp[ip] > xmax && getxmax)
      xmax = xp[ip];

   ip++;
   }
fclose(fpr);

np = ip;
xp = (float *) check_realloc (xp,np*sizeof(float));
yp = (float *) check_realloc (yp,np*sizeof(float));
wp = (float *) check_realloc (wp,np*sizeof(float));

lsqr(xp,yp,wp,np,c0,nterms);

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

dx = (xmax-xmin)/(float)(nx);
for(ix=0;ix<nx;ix++)
   {
   x0 = xmin + (ix+0.5)*dx;

   xb = 1.0;
   y0 = c0[0];
   for(i=1;i<nterms;i++)
      {
      xb = xb*x0;
      y0 = y0 + xb*c0[i];
      }

   fgets(str,512,fpr);
   fprintf(fpw,"%13.5e %13.5e\n",x0,y0);
   }

fclose(fpw);

fprintf(stderr,"nterms= %d:    y = c0",nterms);
for(i=1;i<nterms;i++)
   fprintf(stderr," + c%d*(x^%d)",i,i);

fprintf(stderr,"\n\n");

for(i=0;i<nterms;i++)
   fprintf(stderr,"\tc%d = %13.5e\n",i,c0[i]);
}

void lsqr(float *x,float *y,float *w,int np,float *c0,int order)
{
double *t0, *ts, *b0, *a0, sum, f;
int i, j, k, ip;

t0 = (double *) check_malloc (2*order*sizeof(double));
ts = (double *) check_malloc (order*sizeof(double));
a0 = (double *) check_malloc (order*order*sizeof(double));
b0 = (double *) check_malloc (order*sizeof(double));

for(i=0;i<order*order;i++)
   a0[i] = 0.0;

for(i=0;i<order;i++)
   b0[i] = 0.0;

/*  BRUTE force
for(j=0;j<order;j++)
   {
   for(i=0;i<order;i++)
      {
      sum = 0.0;
      for(ip=0;ip<np;ip++)
         {
	 f = 1.0;
         for(k=0;k<i+j;k++)
            f = f*x[ip];

	 sum = sum + f;
	 }
      a0[i+j*order] = sum;
      }

   sum = 0.0;
   for(ip=0;ip<np;ip++)
      {
      f = 1.0;
      for(k=0;k<j;k++)
	 f = f*x[ip];

      sum = sum + f*y[ip];
      }
   b0[j] = sum;
   }
*/

/*  OLD way without weights 
t0[0] = 1;
for(ip=0;ip<np;ip++)
   {
   t0[1] = x[ip];
   for(i=2;i<2*order;i++)
      t0[i] = t0[i-1]*t0[1];

   ts[0] = y[ip];
   for(i=1;i<order;i++)
      ts[i] = t0[i]*y[ip];

   for(j=0;j<order;j++)
      {
      for(i=0;i<order;i++)
	 a0[i+j*order] = a0[i+j*order] + t0[i+j];

      b0[j] = b0[j] + ts[j];
      }
   }
*/

/*  NEW way w/weights 
*/
t0[0] = 1;
for(ip=0;ip<np;ip++)
   {
   t0[1] = x[ip];
   for(i=2;i<2*order;i++)
      t0[i] = t0[i-1]*t0[1];

   ts[0] = y[ip];
   for(i=1;i<order;i++)
      ts[i] = t0[i]*y[ip];

   for(j=0;j<order;j++)
      {
      for(i=0;i<order;i++)
	 a0[i+j*order] = a0[i+j*order] + t0[i+j]*w[ip];

      b0[j] = b0[j] + ts[j]*w[ip];
      }
   }

gelim_double(a0,order,b0);

for(i=0;i<order;i++)
   c0[i] = b0[i];

free(a0);
free(b0);
free(t0);
free(ts);
}

/*
Gauss elimination without pivoting (no row exchanges):
We solve Ax = b where A is n by n, and x and b have length n
by forming the decomposition A = LU via elimination.
Originally a contains the matrix A and b contains the
vector b.  At the end a contains the lower and upper
triangular matrices L and U and b contains the solution
vector x.  The diagonal of a contains the diagonal of U
(the pivots) since the diagonal elements of U are all 1's.
*/

gelim_double(a,n,b)
double *a, *b;
int n;
   {
	int i, j, jj;
	double *pa, *paj;
	double f, pivot;

         for (j=0; j<n-1; j++)   /* lu decomp of a (no row exchanges) */
	    {
		pa= a + j*n;
		pivot = pa[j];
	 	for (i=j+1; i<n; i++)
		   {
			pa = a + i*n;
			paj= a + j*n;
	     		f = pa[j]/pivot;
			pa[j] = f;
			for (jj=j+1; jj<n; jj++) pa[jj] -= f*paj[jj];
		   }
	   }
         for (i=1; i<n; i++)        /* forward elimination on b */
	    {
		pa = a + i*n;
		for (j=0; j<i; j++) b[i] -= pa[j]*b[j];
	    }
         for (i=n-1; i>-1; i--)        /* back-substitution */
	    {
		pa = a + i*n;
		for (j=n-1; j>i; j--) b[i] -= pa[j]*b[j];
		b[i] = b[i]/pa[i];
	    }
    }
