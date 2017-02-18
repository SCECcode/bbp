#include "include.h"
#include "structure.h"
#include "function.h"

#define TAP_PERC 0.20

struct complex
   {
   float re;
   float im;
   };

int size_float = sizeof(float);
int size_int = sizeof(int);

main(int ac, char **av)
{
struct statdata head1;
float *sx, *space, newdf;
float df, dt, *s;
int nt6, nt_p2, i, j, nt, newnt;
float sec, epi, az, baz, sec1;
int hr, min;
float flo = 0.0;
float fhi = 0.0;

int norm2 = 0;
int norm1 = 0;
int period = 0;
int resamp = 1;

float smoothlen = 0.0;
float tap_per = TAP_PERC;

FILE *fpr, *fpw, *fopfile();
char header[128];
char infile[128];
char outfile[128];
char stitle[128];

int inbin = 0;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("inbin","d",&inbin);
getpar("norm1","d",&norm1);
getpar("norm2","d",&norm2);
getpar("period","d",&period);
getpar("resamp","d",&resamp);
getpar("smoothlen","f",&smoothlen);
getpar("tap_per","f",&tap_per);
endpar();

s = NULL;
s = read_wccseis(infile,&head1,s,inbin);

nt = head1.nt;
dt = head1.dt;

nt_p2 = getnt_p2(nt);
s = (float *) check_realloc (s,nt_p2*size_float);

if(norm1)
   norm_area(s,nt,&dt);
taper_norm(s,&dt,nt,&tap_per);
zero(s+nt,nt_p2-nt);
forfft(s,nt_p2,-1);
ampspec(norm2,(struct complex *)s,nt_p2);

nt = nt_p2/2;
df = 1.0/(nt_p2*dt);

if(smoothlen > 0.0)
   smooth(s,nt,&df,&smoothlen);

if(resamp > 1)
   {
   newnt = nt*resamp;
   newdf = df/resamp;

   space = (float *) check_malloc (4*newnt*size_float);
   sx = (float *) check_malloc (4*newnt*size_float);

   for(j=0;j<nt;j++)
      sx[j] = s[j];
   for(j=0;j<nt;j++)
      sx[nt+j] = s[nt-1-j];

   resample(sx,&df,2*nt,resamp,&newdf,2*newnt,space);
   }
else
   {
   newnt = nt;
   newdf = df;
   sx = (float *) check_malloc (2*newnt*size_float);
   for(j=0;j<newnt;j++)
      sx[j] = s[j];
   }

fprintf(stderr,"nf= %8d df= %13.5e\n",newnt,newdf);

fpw = fopfile(outfile,"w");

if(period)
   {
   for(j=newnt-1;j>=1;j--)
      fprintf(fpw,"%13.5e %13.5e\n",1.0/(j*newdf),sx[j]);
   }
else
   {
   for(j=1;j<newnt;j++)
      fprintf(fpw,"%13.5e %13.5e\n",j*newdf,sx[j]);
   }

fclose(fpw);
}

ampspec(norm,g,n)
struct complex *g;
int n, norm;
{
float *s;
int i;
float amax = 0.0;

s = (float *) g;

for(i=0;i<n/2;i++)
   {
   s[i] = sqrt(g[i].re*g[i].re + g[i].im*g[i].im);
   if(i && s[i] > amax)
      amax = s[i];
   }

if(norm)
   {
   amax = 1.0/amax;
   for(i=0;i<n/2;i++)
      s[i] = s[i]*amax;
   }
/*
*/
}

smooth(s,n,df,len)
float *s, *df, *len;
int n;
{
float *w, *x, fac, sum;
float amp, l2;
int i, j, nw, nw2;
float pi = 3.14159;

nw = (*len)/(*df);
if(nw%2 == 0)
   nw++;

nw2 = nw/2;
l2 = (nw2+1)*(*df);

w = (float *) check_malloc (nw*sizeof(float));
x = (float *) check_malloc (n*sizeof(float));

fac = pi/(float)(nw2+1);
amp = (*df)/(2.0*l2);
sum = 0.0;
for(i=-nw2;i<=nw2;i++)
   {
   w[i+nw2] = amp*(1.0 - cos(i*fac));
   sum = sum + w[i+nw2];
   }

fac = 1.0/sum;
sum = 0.0;
for(j=0;j<nw;j++)
   {
   w[j] = w[j]*fac;
   sum = sum + w[j];
   }

printf("df=%13.5e nw=%d sum=%13.5e\n",*df,nw,sum);

for(i=nw2;i<n-(nw2);i++)
   {
   x[i] = 0.0;
   for(j=0;j<nw;j++)
      x[i] = x[i] + s[j+i-nw2]*w[j];
   }

for(i=nw2;i<n-(nw2);i++)
   s[i] = x[i];

free(w);
free(x);
}

zero(s,n)
float *s;
int n;
{
while(n--)
   {
   s[0] = 0.0;
   s++;
   }
}

taper_norm(g,dt,nt,tap_per)
float *g, *dt, *tap_per;
int nt;
{
float fac, df, arg;
int i;
int ntap;

ntap = nt*(*tap_per);

for(i=0;i<nt-ntap;i++)
   g[i] = g[i]*(*dt);

if(ntap > 0)
   {
   df = 3.14159/(float)(ntap);
   for(i=nt-ntap;i<nt;i++)
      {
      arg = (i-(nt-(ntap+1)))*df;
      fac = (*dt)*0.5*(1.0 + cos(arg));
      g[i] = g[i]*fac;
      }
   }
}

getnt_p2(nt)
int nt;
{
int i = 0;

while(nt > 1)
   {
   if(nt%2)
      nt++;

   nt = nt/2;
   i++;
   }

nt = 1;
while(i--)
   nt = nt*2;

return(nt);
}

norm_area(s,n,dt)
float *s, *dt;
int n;
{
int i, n1;
float area;

n1 = n - 1;

area = 0.5*(*dt)*(s[0] + s[n1]);
for(i=1;i<n1;i++)
   area = area + (*dt)*s[i];

for(i=0;i<n;i++)
   s[i] = s[i]/area;
}

resample(s,olddt,oldnt,isamp,newdt,newnt,p)
float *s, *p, *olddt, *newdt;
int isamp, oldnt, newnt;
{
struct complex *sc;
float fac;
int i, nt_p2;

sc = (struct complex *) s;

nt_p2 = getnt_p2(oldnt);
taper_norm(s,olddt,oldnt);
zero(s+oldnt,(nt_p2)-(oldnt));
forfft(sc,nt_p2,-1);

if(isamp > 0)
   cxzero(sc+(nt_p2/2),(newnt-nt_p2)/2);

invdft(sc,p,newnt,1);

fac = 1.0/((*newdt)*newnt);
taper_norm(s,&fac,newnt);
}
 
cxzero(p,n)
struct complex *p;
int n;
{
float zap = 0.0;
int i;
 
for(i=0;i<n;i++)
   p[i].re = p[i].im = zap;
}
