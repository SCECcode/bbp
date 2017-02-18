#include "include.h"
#include "structure.h"
#include "function.h"

#define AMP 3.0

struct complex
   {
   float re;
   float im;
   };

int size_float = sizeof(float);
int size_int = sizeof(int);

main(ac,av)
int ac;
char **av;
{
struct statdata shead, ghead;
struct complex *sc, *gc;
float *s, *g;
int nt_p2, i, j, nsrc, nt;
int ntshft;
float tstart = 0.0;
float samp;
int isamp;
int npad = 1;
int gaus = 0;
int decon = 0;
int integ = 0;
int norm1 = 1;
int ntap = 25;
int nb;
float tzero;
float t0 = 0.0;
double wlevel = 1.0e-15;

float scale = 1.0;

int ntout = -999;
int inbin1 = 0;
int inbin2 = 0;
int outbin = 0;

int order = 0;
int npass = 2;
float f1 = 0.0;
float f2 = 1.0e+15;

FILE *fpr1, *fpr2, *fpw, *fopfile();
char infile1[128];
char infile2[128];
char outfile[128];

setpar(ac,av);
getpar("integ","d",&integ);
getpar("decon","d",&decon);
getpar("ntout","d",&ntout);
getpar("tstart","f",&tstart);
getpar("npad","d",&npad);
getpar("ntap","d",&ntap);
getpar("norm1","d",&norm1);
getpar("gaus","d",&gaus);
if(gaus)
   {
   mstpar("nb","d",&nb);
   mstpar("tzero","f",&tzero);
   }
else
   {
   mstpar("infile1","s",infile1);
   getpar("inbin1","d",&inbin1);
   }

mstpar("infile2","s",infile2);
getpar("inbin2","d",&inbin2);
mstpar("outfile","s",outfile);
getpar("outbin","d",&outbin);
getpar("scale","f",&scale);
getpar("t0","f",&t0);
getpar("wlevel","F",&wlevel);
getpar("order","d",&order);
if(order)
   {
   mstpar("f1","f",&f1);
   mstpar("f2","f",&f2);
   getpar("npass","d",&npass);
   }
endpar();

g = NULL;
s = NULL;

g = read_wccseis(infile2,&ghead,g,inbin2);

if(ghead.hr < -24 || ghead.hr > 23)
   ghead.hr = 0;
if(ghead.min < -60 || ghead.min > 59)
   ghead.min = 0;

shead.hr = shead.min = 0;
shead.sec = 0.0;
if(!gaus)
   s = read_wccseis(infile1,&shead,s,inbin1);
else
   {
   shead.nt = ghead.nt;
   shead.dt = ghead.dt;
   }

ghead.hr = ghead.hr + shead.hr;
ghead.min = ghead.min + shead.min;
ghead.sec = ghead.sec - t0;

samp = 1;
isamp = 1;
if(ghead.dt != shead.dt)
   {
   samp = ghead.dt/shead.dt;
   if(ghead.dt < shead.dt)
      {
      isamp = (int)(1.0/samp + 0.5);
      if(isamp*ghead.dt != shead.dt)
	 {
         fprintf(stderr,"isamp=%d\n",isamp);
	 /*
	 fprintf(stderr,"*** dt is not an integral of dt2, exiting...\n");
	 exit(-1);
	 */
	 }
      }
   else if(ghead.dt > shead.dt)
      {
      isamp = (int)(samp);
      if(isamp*shead.dt != ghead.dt)
	 {
	 fprintf(stderr,"isamp=%d\n",isamp);
	 /*
	 fprintf(stderr,"*** dt2 is not an integral of dt, exiting...\n");
	 exit(-1);
	 */
	 }
      }
   }

nt = shead.nt + samp*ghead.nt;
nt_p2 = npad*getnt_p2(nt);
s = (float *) check_realloc (s,nt_p2*size_float);
g = (float *) check_realloc (g,nt_p2*size_float);

sc = (struct complex *) s;
gc = (struct complex *) g;

if(gaus)
   makesource(s,&shead.dt,shead.nt,&tzero,nb);

if(norm1)
   norm_area(s,shead.nt,&shead.dt);

if(samp != 1)
   {
   resample(g,&samp,ghead.nt);
   ghead.nt = ghead.nt*samp;
   ghead.dt = shead.dt;
   }

taper_norm(s,&shead.dt,shead.nt,ntap);
zero(s+shead.nt,(nt_p2)-shead.nt);
forfft(sc,nt_p2,-1);

taper_norm(g,&ghead.dt,ghead.nt,ntap);
zero(g+ghead.nt,nt_p2-ghead.nt);
forfft(gc,nt_p2,-1);

if(integ)
   integrate(gc,nt_p2,&ghead.dt);
else
   convolve(sc,gc,nt_p2,decon,(int)(t0/ghead.dt),&wlevel,&ghead.dt);

if(order)
   wfilter(gc,&ghead.dt,nt_p2,&f1,&f2,order,npass);

invfft(gc,nt_p2,1);
norm(g,&ghead.dt,nt_p2);

if(decon)
   nt = nt_p2;
else
   {
   ghead.sec = ghead.sec + shead.sec;
   if(tstart > ghead.sec)
      {
      ntshft = (tstart - ghead.sec)/ghead.dt;

      if(nt+ntshft > nt_p2)
         {
         for(i=0;i<nt_p2-ntshft;i++)
            g[i] = g[i+ntshft];
         for(i=nt_p2-ntshft;i<nt;i++)
            g[i] = 0.0;
         }
      else
         {
         for(i=0;i<nt;i++)
            g[i] = g[i+ntshft];
         }
      ghead.sec = tstart;
      }
   else if(tstart < ghead.sec)
      {
      ntshft = (ghead.sec - tstart)/ghead.dt;
      nt = nt + ntshft;

      if(nt > nt_p2)
         g = (float *) check_realloc (g,nt*size_float);

      for(i=nt-1;i>=ntshft;i--)
         g[i] = g[i-ntshft];
      for(i=ntshft-1;i>=0;i--)
         g[i] = 0.0;

      /*
      for(i=nt-1;i>=ntshft;i--)
         g[i] = g[i-ntshft];
      for(i=ntshft-1;i>=0;i--)
         g[i] = 0.0;
*/

      ghead.sec = tstart;
      }
   }

ghead.nt = nt;
if(ntout > -99)
   {
   if(ntout > nt)
      {
      g = (float *) check_realloc (g,ntout*size_float);

      for(i=nt;i<ntout;i++)
         g[i] = 0.0;
      }

   ghead.nt = ntout;
   }

fprintf(stderr,"**** nt_p2= %d\n",nt_p2);
fprintf(stderr,"**** output nt= %d\n",ghead.nt);

for(i=0;i<ghead.nt;i++)
      g[i] = scale*g[i];

write_wccseis(outfile,&ghead,g,outbin);
}

integrate(g,n,dt)
struct complex *g;
float *dt;
int n;
{
float tmpre, dw, invom;
float twopi = 6.28318531;
float one = 1.0;
int i;

dw = twopi/(n*(*dt));

g[0].re = g[0].im = 0.0;
for(i=1;i<n/2;i++)
   {
   invom = -one/(i*dw);
   tmpre = -g[i].im*invom;
   g[i].im = g[i].re*invom;
   g[i].re = tmpre;
   }
}

convolve(s,g,n,dc,it0,wlev,dt)
struct complex *s, *g;
float *dt;
double *wlev;
int n, dc, it0;
{
double tmpre, tmpim, den;
float arg, cosA, sinA, fac, ff, fmax, df;
float twopi = 6.28318531;
float one = 1.0;
double max = 0.0;
int i;

fmax = 1.0e+15;

if(dc)
   {
   den = s[0].re*s[0].re;
   if(den > max)   
      max = den;
    
   den = s[0].im*s[0].im;
   if(den > max)   
      max = den;
       
   for(i=1;i<n/2;i++) 
      { 
      den = (s[i].re*s[i].re + s[i].im*s[i].im);
      if(den > max) 
         max = den; 
      } 
   max = max*(*wlev);

   den = s[0].re*s[0].re;
   if(den < max)
      den = max;

   den = one/den;
   g[0].re = g[0].re*s[0].re*den;
       
   den = s[0].im*s[0].im;
   if(den < max)
      den = max;
 
   den = one/den; 
   g[0].im = g[0].im*s[0].im*den;

   df = 1.0/((*dt)*n);
   for(i=1;i<n/2;i++)
      {
      ff = i*df;

      if(ff < fmax)
	 {
         den = (s[i].re*s[i].re + s[i].im*s[i].im);
         if(den < max)
	    {
	    den = max;
/*
fprintf(stderr,"%d re=%13.5e m=%13.5e\n",i,s[i].re,s[i].im);
*/
	    }

         den = one/den;
         tmpre   = (g[i].re*s[i].re + g[i].im*s[i].im)*den;
         tmpim = (g[i].im*s[i].re - g[i].re*s[i].im)*den;

         g[i].re = tmpre;
         g[i].im = tmpim;
	 }
      else
	 {
         g[i].re = 0.0;
         g[i].im = 0.0;
	 }
      }


   if(it0)
      {
      fac = -it0*twopi/n;
      arg = 0.5*n*fac;
      cosA = cos(arg);
      g[0].im = g[0].im*cosA;

      for(i=1;i<n/2;i++)
         {
         arg = fac*i;
         cosA = cos(arg);
         sinA = sin(arg);

         tmpre   = g[i].re*cosA - g[i].im*sinA;
         g[i].im = g[i].re*sinA + g[i].im*cosA;
         g[i].re = tmpre;
         }
      }
   }
else
   {
   g[0].re = g[0].re*s[0].re;
   g[0].im = g[0].im*s[0].im;
   for(i=1;i<n/2;i++)
      {
      tmpre   = g[i].re*s[i].re - g[i].im*s[i].im;
      g[i].im = g[i].re*s[i].im + g[i].im*s[i].re;
      g[i].re = tmpre;
      }
   }
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

taper_norm(g,dt,nt,ntap)
float *g, *dt;
int nt, ntap;
{
float fac;
int i;
float pi = 3.141592654;

for(i=0;i<nt-ntap;i++)
   g[i] = g[i]*(*dt);

for(i=nt-ntap;i<nt;i++)
   {
   fac = (*dt)*0.5*(1.0 + cos(pi*((float)(nt-i)/(float)(ntap))));
   g[i] = g[i]*fac;
   }
}

norm(g,dt,nt)
float *g, *dt;
int nt;
{
float fac;

fac = 1.0/((*dt)*nt);
while(nt--)
   {
   g[0] = g[0]*fac;
   g++;
   }
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
 
makesource(s,dt,nt,a,nb)
float *s;
float *dt, *a;
int nt, nb;
{
double exp(), exp_arg;
int it;
float t, t2, tshift, alpha;
float src_amp = 1.0e+05;
 
tshift = AMP*(*a);
alpha = -1.0/((*a)*(*a));
s[0] = 0.0;
for(it=1;it<nt;it++)
   {
   t = (it-1)*(*dt) - tshift;
   t2 = t*t;
      exp_arg = alpha*t2;
      s[it] = src_amp*exp(exp_arg);
      if(nb != 1)
         s[it] = -t*s[it];
   /*
   if(t <= tshift)
      {
      exp_arg = alpha*t2;
      s[it] = src_amp*exp(exp_arg);
      if(nb != 1)
         s[it] = -t*s[it];
      }
   else
      s[it] = 0.0;
   */
   }
}

resample(g,samp,nt)
float *g, *samp;
int nt;
{
float fac;
int i, isamp, it;
float *s;

if(*samp < 1.0)
   {
   isamp = (int)(1.0/(*samp) + 0.5);

   for(it=0;it<nt*(*samp);it++)
      g[it] = g[it*isamp];
   }
else
   {
   s = (float *) check_malloc (((int)(nt*(*samp)))*sizeof(float));
   isamp = (int)(*samp);
   for(it=0;it<nt-1;it++)
      {
      fac = (g[it+1] - g[it])/(*samp);
      for(i=0;i<isamp;i++)
         s[it*isamp + i] = g[it] + i*fac;
      }

   for(it=(nt-1)*isamp;it<(nt*isamp);it++)
      s[it] = g[(nt-1)*isamp-1];

   for(it=0;it<(nt*isamp);it++)
      g[it] = s[it];

   free(s);
   }
}

wfilter(g,dt,n,f1,f2,ord,np)
struct complex *g;
float *dt, *f1, *f2;
int n, ord, np;
{
int i, j;
float df, ff, fac, flo, fhi;
float one = 1.0;

df = 1.0/((*dt)*n);
g[0].im = 0.0;
for(i=1;i<n/2;i++)
   {
   ff = i*df;

   fhi = (*f1)/ff;
   for(j=0;j<ord;j++)
      fhi = fhi*fhi;

   flo = ff/(*f2);
   for(j=0;j<ord;j++)
      flo = flo*flo;

   fac = one/((one + flo)*(one + fhi));
   if(np == 1)
      fac = sqrt(fac);

   g[i].re = fac*g[i].re;
   g[i].im = fac*g[i].im;
   }
}

wfilter2(g,dt,n,f1,f2,f3,f4)
struct complex *g;
float *dt, *f1, *f2, *f3, *f4;
int n;
{
int i;
float df, ff, fac, df1, df2, arg;
float pi = 3.141592654;
float pi2 = 1.570796327;

df1 = *f2 - *f1;
df2 = *f4 - *f3;

df = 1.0/((*dt)*n);
for(i=0;i<n/2;i++)
   {
   ff = i*df;

   if(ff < *f1 || ff > *f4)
      fac = 0.0;
   else if(ff < *f2)
      {
      arg = pi*(ff-(*f1))/df1;
      fac = 0.5*(1.0 - cos(arg));
      }
   else if(ff > *f3)
      {
      arg = pi*(ff-(*f3))/df2;
      fac = 0.5*(1.0 + cos(arg));
      }
   else
      fac = 1.0;

   g[i].re = fac*g[i].re;
   g[i].im = fac*g[i].im;
   }
}
