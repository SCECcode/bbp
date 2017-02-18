#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

void t1t2(float *, int, float *, float *, float *, float *, float *, int);
void gelim_double(double *, int, double *);
void invmat3x3(float *, float *);
void invmat4x4(float *, float *);
void baseline(float *, int, float *, int);
void get_trend(float *, int, float *, float *);
void detrend(float *, int, float *, float *);
void remean(float *, int, float *, float *, float *);
void demean(float *, int, float *, float *);
void tapr(float *, int, int, int, int, int);
void retrend(float *, int, float *);
void differ(float *, int, float *, float *);
void integrate(float *, int, float *, float *);

int size_float = sizeof(float);
int size_int = sizeof(int);

int main (ac, av)
int ac;
char **av;
{
struct statdata shead;
float *seis;
float dt, sec, edist, az, baz;
int nt, hr, min, i, j, nt6;
int integ = 0;
int diff = 0;
int rtrend = 0;
int rbase = 0;
int dmean = 0;
int dtrend = 0;
int rmean = 0;
int taper = 0;
int boorebase = 0;

float tstart = -1.0e+15;
float tlen = -1.0e+15;

float scale = 1.0;

float t1 = -1.0;
float t2 = -1.0;
float tf1 = -1.0;
float tf2 = -1.0;
int v0correct = 0;

float ts0 = 0.0;
float te0 = 0.0;
float ts1 = 9999.0;
float te1 = 9999.0;

float tfront = -1.0;
float tend = -1.0;

float finaldisp = 0.0;
float init_val = 0.0;

char filein[256];
char fileout[256];
char title[128], header[128];

int inbin = 0;
int outbin = 0;

setpar(ac, av);
mstpar("filein","s",filein);
mstpar("fileout","s",fileout);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
getpar("integ","d",&integ);
getpar("diff","d",&diff);
getpar("rtrend","d",&rtrend);
getpar("rbase","d",&rbase);
getpar("dmean","d",&dmean);
getpar("dtrend","d",&dtrend);

getpar("boorebase","d",&boorebase);
if(boorebase)
   {
   getpar("t1","f",&t1);
   getpar("t2","f",&t2);
   getpar("tf1","f",&tf1);
   getpar("tf2","f",&tf2);
   getpar("v0correct","d",&v0correct);
   }

getpar("rmean","d",&rmean);
if(rmean)
   {
   getpar("tstart","f",&tstart);
   getpar("tlen","f",&tlen);
   }

getpar("taper","d",&taper);
if(taper)
   {
   getpar("ts0","f",&ts0);
   getpar("te0","f",&te0);
   getpar("ts1","f",&ts1);
   getpar("te1","f",&te1);

   getpar("tfront","f",&tfront);
   getpar("tend","f",&tend);
   }

getpar("scale","f",&scale);
getpar("finaldisp","f",&finaldisp);
getpar("init_val","f",&init_val);
endpar();

if(finaldisp != 0.0)
   dmean = 1;

seis = NULL;
seis = read_wccseis(filein,&shead,seis,inbin);

if(taper)
   {
   if(tfront > 0.0)
      {
      ts0 = 0.0;
      te0 = tfront;
      }
   if(tend > 0.0)
      {
      ts1 = shead.nt*shead.dt - tend;
      te1 = shead.nt*shead.dt;
      }
   }

if(rtrend) /*  remove linear trend  */
   retrend(seis,shead.nt,&shead.dt);

if(rbase) /*  remove baseline polynomial  */
   baseline(seis,shead.nt,&shead.dt,rbase);

/*
   demean(): adjust 1st point of time history so that when it is integrated,
             the resulting time history will have zero mean as well.
*/
if(dmean)
   demean(seis,shead.nt,&shead.dt,&finaldisp);

/*
   detrend(): adjust 1st point of time history so that when it is integrated,
              the resulting time history will have no static offset.
*/
if(dtrend)
   detrend(seis,shead.nt,&shead.dt,&init_val);

if(rmean) /*  remove mean */
   remean(seis,shead.nt,&shead.dt,&tstart,&tlen);

if(boorebase) /*  remove tri-linear from vel */
   t1t2(seis,shead.nt,&shead.dt,&t1,&t2,&tf1,&tf2,v0correct);

if(integ) /*  integrate  */
   integrate(seis,shead.nt,&shead.dt,&init_val);

if(diff) /*  differentiate  */
   differ(seis,shead.nt,&shead.dt,&init_val);

if(taper) /*  smoothly taper record */
   tapr(seis,shead.nt,(int)(ts0/shead.dt),(int)(te0/shead.dt),(int)(ts1/shead.dt),(int)(te1/shead.dt));

for(i=0;i<shead.nt;i++)
   seis[i] = scale*seis[i];

write_wccseis(fileout,&shead,seis,outbin);
}

void integrate(s,nt,dt,iv)
float *s, *dt, *iv;
int nt;
{
float dt3, dt8, s0, s1, s2, s3, s4;
float s_start, s_end;
int it;

float four = 4.0;
float six = 6.0;

dt3 = (*dt)/3.0;
dt8 = (*dt)/8.0;

/*
   simple rule 0
*/

s0 = *iv;
for(it=0;it<nt;it++)
   {
   s0 = s[it]*(*dt) + s0;
   s[it] = s0;
   }

/*
   trapezoid rule 0

s0 = 0.0;
s0 = 0.5*(s[0])*(*dt);
s0 = *iv;
for(it=1;it<nt;it++)
   {
   s1 = 0.5*(s[it-1] + s[it])*(*dt) + s0;
   s[it-1] = s0;
   s0 = s1;
   }
s[nt-1] = s1;
*/

/*
   trapezoid rule 1

s[0] = s[0]*(*dt);
s[0] = *iv;
for(it=1;it<nt;it++)
   s[it] = s[it-1] + s[it]*(*dt);
*/

/*
   trapezoid rule 2

s0 = dt8*(six*s[0] + s[1]);
s0 = *iv;
for(it=1;it<nt-1;it++)
   {
   s1 = dt8*(s[it-1] + six*s[it] + s[it+1]) + s0;
   s[it-1] = s0;
   s0 = s1;
   }
s[nt-1] = dt8*(s[nt-2] + six*s[nt-1]) + s0;
s[nt-2] = s0;
*/

/*
   Simpson's rule

s0 = dt3*(s[0]);
s1 = dt3*(four*s[0] + s[1]);
for(it=2;it<nt;it++)
   {
   s2 = dt3*(s[it-2] + four*s[it-1] + s[it]) + s0;
   s[it-2] = s0;
   s0 = s1;
   s1 = s2;
   }
s[nt-2] = s1;
s[nt-1] = s2;
*/
}

void differ(s,nt,dt,iv)
float *s, *dt, *iv;
int nt;
{
float c0, c1;
float tot = 0.0;
int it;

c0 = 1.0/(*dt);

/* backward difference -> consistent with simple integration
   assumes s=init_val for t<0 */

for(it=nt-1;it>0;it--)
   s[it] = (s[it] - s[it-1])*c0;

s[0] = (s[0]-(*iv))*c0;

/* forward difference

for(it=0;it<nt-1;it++)
   s[it] = (s[it+1] - s[it])*c0;

s[nt-1] = s[nt-2];
*/
}

void retrend(s,nt,dt)
float *s, *dt;
int nt;
{
float m0, b0, t, y;
float c0 = 0.0;
float c1 = 0.0;
float c2 = 0.0;
float c3 = 0.0;
int it;

for(it=0;it<nt;it++)
   {
   t = it*(*dt);
   c0 = c0 + t;
   c1 = c1 + s[it];
   c2 = c2 + t*t;
   c3 = c3 + t*s[it];
   }

m0 = (c0*c1 - nt*c3)/(c0*c0 - nt*c2);
b0 = (c0*c3 - c1*c2)/(c0*c0 - nt*c2);

for(it=0;it<nt;it++)
   {
   t = it*(*dt);

   y = m0*t + b0;
   s[it] = s[it] - y;
   }
}

void tapr(s,nt,it0,it1,it2,it3)
float *s;
int nt, it0, it1, it2, it3;
{
float fac;
float pi = 3.14159265;
int it;

fac = pi/(float)(it1-it0);
for(it=0;it<it0;it++)
   s[it] = 0.0;
for(it=it0;it<it1;it++)
   s[it] = s[it]*0.5*(1.0 - cos((it-it0)*fac));

if(it2 < nt)
   {
   if(it3 > nt)
      it3 = nt;

   fac = pi/(float)(it3-it2);
   for(it=it2;it<it3;it++)
      s[it] = s[it]*0.5*(1.0 + cos((it-it2)*fac));
   for(it=it3;it<nt;it++)
      s[it] = 0.0;
   }
}

void demean(s,nt,dt,fdisp)
float *s, *dt, *fdisp;
int nt;
{
int it, j, test;
float old, new, oabs, nabs, nsgn;
float difp, tolp, tolm;
float tol = 0.1;
float sum = 0.0;

/*
   adjust start of time history so that when it is integrated,
   the resulting time history will have zero mean as well.
*/

tolp = 1.0 + tol;
tolm = 1.0 - tol;

j = 0;
test = 1;
while(test)
   {
   sum = 0.0;
   for(it=0;it<j;it++)
      sum = sum + (nt-it)*s[it];

   for(it=j+1;it<nt;it++)
      sum = sum + (nt-it)*s[it];

   new = *fdisp - sum/(nt-j);
   nabs = new;
   nsgn = 1.0;
   if(nabs < 0.0)
      {
      nsgn = -1.0;
      nabs = -new;
      }

   old = s[j];
   oabs = old;
   if(oabs < 0.0)
      oabs = -old;

   difp = nabs/oabs;
   if(difp > tolp)
      new = tolp*oabs*nsgn;
   else if(difp < tolm)
      new = tolm*oabs*nsgn;
   else
      test = 0;

   s[j] = new;
   j++;
   }
printf("demean(): j=%d\n",j);
}

void remean(s,nt,dt,ts,tl)
float *s, *dt, *ts, *tl;
int nt;
{
int it, its, ite;
float sum = 0.0;

if(*ts < 0.0)
   its = 0;
else
   its = (*ts)/(*dt);

if(its < 0)
   its = 0;

if(*tl < 0.0)
   ite = nt;
else
   ite = its + (*tl)/(*dt);

if(ite > nt)
   ite = nt;

/*
   remove mean from input time history
*/

sum = 0.0;
for(it=its;it<ite;it++)
   sum = sum + s[it];

sum = sum/(float)(ite-its);
fprintf(stderr,"remean(): its= %d ite= %d mean= %13.5e\n",its,ite,sum);

for(it=0;it<nt;it++)
   s[it] = s[it] - sum;
}

void detrend(s,nt,dt,iv)
float *s, *dt, *iv;
int nt;
{
int it, np;
float *s1, b0;

np = nt/4;

s1 = (float *) check_malloc (nt*sizeof(float));

for(it=0;it<nt;it++)
   s1[it] = s[it];

integrate(s1,nt,dt,iv);
get_trend(s1+(nt-np),np,dt,&b0);

s[0] = s[0] - b0/(*dt);
printf("detrend(): b0=%f\n",b0);
free(s1);
}

void get_trend(s,nt,dt,b0)
float *s, *dt, *b0;
int nt;
{
int it;

*b0 = 0.0;
for(it=0;it<nt;it++)
   *b0 = *b0 + s[it];

*b0 = *b0/((float)(nt));
}

void ck_newline(str,n)
char *str;
int n;
{
int i = 0;

while(i < n && str[i] != '\0')
   i++;
if(str[i-1] != '\n')
   str[i-1] = '\n';
}

void baseline(s,nt,dt,order)
float *s, *dt;
int nt, order;
{
float y;
double t0[11], ts[6];
double x0[6], c0[36];
int i, j, it;

for(i=0;i<36;i++)
   c0[i] = 0.0;

for(i=0;i<6;i++)
   x0[i] = 0.0;

t0[0] = 1;
for(it=0;it<nt;it++)
   {
   t0[1] = it*(*dt);
   for(i=2;i<2*order+1;i++)
      t0[i] = t0[i-1]*t0[1];

   ts[0] = s[it];
   for(i=1;i<order+1;i++)
      ts[i] = t0[i]*s[it];

   for(j=0;j<order+1;j++)
      {
      for(i=0;i<order+1;i++)
	 c0[i+j*(order+1)] = c0[i+j*(order+1)] + t0[i+j];

      x0[j] = x0[j] + ts[j];
      }
   }

/*
for(j=0;j<order+1;j++)
   {
   for(i=0;i<order+1;i++)
      fprintf(stderr,"c[%d,%d]= %13.5e ",j,i,c0[i+j*(order+1)]);

   fprintf(stderr,": x[%d]= %13.5e\n",j,x0[j]);
   }
*/

gelim_double(c0,order+1,x0);

/*
for(i=0;i<order+1;i++)
   fprintf(stderr,"x[%d]= %13.5e\n",i,x0[i]);
   */

for(it=0;it<nt;it++)
   {
   t0[1] = it*(*dt);
   for(i=2;i<order+1;i++)
      t0[i] = t0[i-1]*t0[1];

   y = x0[0];
   for(i=1;i<order+1;i++)
      y = y + x0[i]*t0[i];

   s[it] = s[it] - y;
   }
}

void invmat3x3(inv,mat)
float *inv, *mat;
{
int i;
float ftmp, det;

/* form cofactors */
 
inv[0] = mat[4]*mat[8] - mat[5]*mat[7];
inv[1] = -mat[3]*mat[8] + mat[5]*mat[6];
inv[2] = mat[3]*mat[7] - mat[4]*mat[6];
 
inv[3] = -mat[1]*mat[8] + mat[2]*mat[7];
inv[4] = mat[0]*mat[8] - mat[2]*mat[6];
inv[5] = -mat[0]*mat[7] + mat[1]*mat[6];
 
inv[6] = mat[1]*mat[5] - mat[2]*mat[4];
inv[7] = -mat[0]*mat[5] + mat[2]*mat[3];
inv[8] = mat[0]*mat[4] - mat[1]*mat[3];

/* calculate determinant */
 
det = mat[0]*inv[0] + mat[1]*inv[1] + mat[2]*inv[2];
 
if(det != 0.0)
   {
   det = 1.0/det;
   for(i=0;i<9;i++)
      inv[i] = det*inv[i];
 
   /* transpose */
 
   ftmp = inv[1];
   inv[1] = inv[3];
   inv[3] = ftmp;
 
   ftmp = inv[2];
   inv[2] = inv[6];
   inv[6] = ftmp;

   ftmp = inv[5];
   inv[5] = inv[7];
   inv[7] = ftmp;
   }
else
   {
   fprintf(stderr,"invmat3x3 error: matrix is singular, exiting...\n");
   exit(-9);
   }
}

void invmat4x4(v,m)
float *v, *m;
{
int i;
float ftmp, det;

/* form cofactors */
 
v[0]  =   m[5]*(m[10]*m[15] - m[14]*m[11])
        - m[6]*(m[9]*m[15] - m[13]*m[11])
        + m[7]*(m[9]*m[14] - m[13]*m[10]);
v[1]  = - m[4]*(m[10]*m[15] - m[14]*m[11])
        + m[6]*(m[8]*m[15] - m[12]*m[11])
        - m[7]*(m[8]*m[14] - m[12]*m[10]);
v[2]  =   m[4]*(m[9]*m[15] - m[13]*m[11])
        - m[5]*(m[8]*m[15] - m[12]*m[11])
        + m[7]*(m[8]*m[13] - m[12]*m[9]);
v[3]  = - m[4]*(m[9]*m[14] - m[13]*m[10])
        + m[5]*(m[8]*m[14] - m[12]*m[10])
        - m[6]*(m[8]*m[13] - m[12]*m[9]);
 
v[4]  = - m[1]*(m[10]*m[15] - m[14]*m[11])
        + m[2]*(m[9]*m[15] - m[13]*m[11])
        - m[3]*(m[9]*m[14] - m[13]*m[10]);
v[5]  =   m[0]*(m[10]*m[15] - m[14]*m[11])
        - m[2]*(m[8]*m[15] - m[12]*m[11])
        + m[3]*(m[8]*m[14] - m[12]*m[10]);
v[6]  = - m[0]*(m[9]*m[15] - m[13]*m[11])
        + m[1]*(m[8]*m[15] - m[12]*m[11])
        - m[3]*(m[8]*m[13] - m[12]*m[9]);
v[7]  =   m[0]*(m[9]*m[14] - m[13]*m[10])
        - m[1]*(m[8]*m[14] - m[12]*m[10])
        + m[2]*(m[8]*m[13] - m[12]*m[9]);
  
v[8]  =   m[1]*(m[6]*m[15] - m[14]*m[7])
        - m[2]*(m[5]*m[15] - m[13]*m[7])
        + m[3]*(m[5]*m[14] - m[13]*m[6]);
v[9]  = - m[0]*(m[6]*m[15] - m[14]*m[7])
        + m[2]*(m[4]*m[15] - m[12]*m[7])
        - m[3]*(m[4]*m[14] - m[12]*m[6]);
v[10] =   m[0]*(m[5]*m[15] - m[13]*m[7])
        - m[1]*(m[4]*m[15] - m[12]*m[7])
        + m[3]*(m[4]*m[13] - m[12]*m[5]);
v[11] = - m[0]*(m[5]*m[14] - m[13]*m[6])
        + m[1]*(m[4]*m[14] - m[12]*m[6])
        - m[2]*(m[4]*m[13] - m[12]*m[5]);
  
v[12] = - m[1]*(m[6]*m[11] - m[10]*m[7])
        + m[2]*(m[5]*m[11] - m[9]*m[7])
        - m[3]*(m[5]*m[10] - m[9]*m[6]);
v[13] =   m[0]*(m[6]*m[11] - m[10]*m[7])
        - m[2]*(m[4]*m[11] - m[8]*m[7])
        + m[3]*(m[4]*m[10] - m[8]*m[6]);
v[14] = - m[0]*(m[5]*m[11] - m[9]*m[7])
        + m[1]*(m[4]*m[11] - m[8]*m[7])
        - m[3]*(m[4]*m[9] - m[8]*m[5]);
v[15] =   m[0]*(m[5]*m[10] - m[9]*m[6])
        - m[1]*(m[4]*m[10] - m[8]*m[6])
        + m[2]*(m[4]*m[9] - m[8]*m[5]);
 
/* calculate determinant */
 
det = m[0]*v[0] + m[1]*v[1] + m[2]*v[2] + m[3]*v[3];
 
if(det != 0.0)
   {
   det = 1.0/det;
   for(i=0;i<16;i++)
      v[i] = det*v[i];
 
   /* transpose */
 
   ftmp = v[1];
   v[1] = v[4];
   v[4] = ftmp;
 
   ftmp = v[2];
   v[2] = v[8];
   v[8] = ftmp;
 
   ftmp = v[3];
   v[3] = v[12];
   v[12] = ftmp;
 
   ftmp = v[7];
   v[7] = v[13];
   v[13] = ftmp;
 
   ftmp = v[11];
   v[11] = v[14];
   v[14] = ftmp;
 
   ftmp = v[6];
   v[6] = v[9];
   v[9] = ftmp;
   }
else
   {
   fprintf(stderr,"invmat4x4 error: matrix is singular, exiting...\n");
   exit(-9);
   }
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

void gelim_double(a,n,b)
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

void t1t2(s,nt,dt,t1,t2,tf1,tf2,iv0)
float *s, *dt, *t1, *t2, *tf1, *tf2;
int nt, iv0;
{
int it, it1, it2, itf1, itf2;
float *v;
float m0, af, am;
float b0, t, y;
float c0 = 0.0;
float c1 = 0.0;
float c2 = 0.0;
float c3 = 0.0;

v = (float *)check_malloc(nt*sizeof(float));

if((*t1) < 0.0)
   {
   m0 = 50.0;
   for(it=0;it<nt,(*t1)<0.0;it++)
      {
      if(s[it] > m0 || -s[it] > m0)
	 *t1 = it*(*dt);
      }
   }

if((*t2) < 0.0)
   {
   m0 = 50.0;
   for(it=nt-1;it>=0,(*t2)<0.0;it--)
      {
      if(s[it] > m0 || -s[it] > m0)
	 *t2 = it*(*dt);
      }
   }

if((*t2) > (nt-1)*(*dt))
   *t2 = 0.5*(*t1 + nt*(*dt));

if((*tf2) < 0.0)
   *tf2 = nt*(*dt);
if((*tf1) < 0.0)
   *tf1 = *t2 + 0.1*(*tf2 - *t2);

it1 = (int)((*t1)/(*dt) + 1.5);
it2 = (int)((*t2)/(*dt) + 1.5);
itf1 = (int)((*tf1)/(*dt) + 1.5);
itf2 = (int)((*tf2)/(*dt) + 1.5);

if(it1 < 0 || it1 >= nt)
   return;
if(it2 < 0 || it2 >= nt)
   return;
if(it2 < it1)
   return;

if(itf2 > nt)
   itf2 = nt;

m0 = 0.0;
for(it=0;it<it1/2;it++)
   m0 = m0 + s[it];

m0 = m0/(float)(it1/2);
for(it=0;it<nt;it++)
   s[it] = s[it] - m0;;

m0 = 0.0;
for(it=0;it<nt;it++)
   {
   m0 = s[it]*(*dt) + m0;
   v[it] = m0;
   }

for(it=itf1;it<itf2;it++)
   {
   t = it*(*dt);
   c0 = c0 + t;
   c1 = c1 + v[it];
   c2 = c2 + t*t;
   c3 = c3 + t*v[it];
   }

af = (c0*c1 - (itf2-itf1)*c3)/(c0*c0 - (itf2-itf1)*c2);
b0 = (c0*c3 - c1*c2)/(c0*c0 - (itf2-itf1)*c2);

if(iv0)
   {
   m0 = -b0/af;
   if(m0 > *t1 && m0 < *tf2)
      *t2 = m0;

   it2 = (int)((*t2)/(*dt) + 1.5);
   }

fprintf(stderr,"t1= %f t2= %f tf1= %f tf2= %f\n",*t1,*t2,*tf1,*tf2);

am = (af*(*t2) + b0)/(*t2 - *t1);

for(it=it1;it<it2;it++)
   s[it] = s[it] - am;

for(it=it2;it<nt;it++)
   s[it] = s[it] - af;

free(v);
}
