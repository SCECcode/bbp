#include "include.h"
#include "structure.h"
#include "function.h"

#define TAP_PERC 0.05
#define SPEC_BLOCK 10000

struct complex
   {
   float re;
   float im;
   };

int size_float = sizeof(float);
int size_int = sizeof(int);

main(int ac,char **av)
{
FILE *fpr, *fopfile();
struct statdata head1;
float *s1, vref, vsite, vrat;
int nt_p2;
float *fsp, *asp;
int ns, n;

float tap_per = TAP_PERC;

char str[512];
char infile[128];
char specfile[128];
char outfile[128];

int inbin = 0;
int outbin = 0;

setpar(ac,av);
mstpar("infile","s",infile);
mstpar("specfile","s",specfile);
mstpar("outfile","s",outfile);
getpar("tap_per","f",&tap_per);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
endpar();

fsp = (float *)check_malloc(SPEC_BLOCK*sizeof(float));
asp = (float *)check_malloc(SPEC_BLOCK*sizeof(float));

fsp[0] = 0.0;
asp[0] = 1.0;

fpr = fopfile(specfile,"r");

ns = 1;
n = 1;
while(fgets(str,512,fpr) != NULL)
   {
   if(ns+1 == n*SPEC_BLOCK)
      {
      n++;
      fsp = (float *)check_realloc(fsp,n*SPEC_BLOCK*sizeof(float));
      asp = (float *)check_realloc(asp,n*SPEC_BLOCK*sizeof(float));
      }

   sscanf(str,"%f %f",&fsp[ns],&asp[ns]);
   ns++;
   }

fclose(fpr);

if(ns+1 == n*SPEC_BLOCK)
   {
   n++;
   fsp = (float *)check_realloc(fsp,n*SPEC_BLOCK*sizeof(float));
   asp = (float *)check_realloc(asp,n*SPEC_BLOCK*sizeof(float));
   }

fsp[ns] = 1.0e+15;
asp[ns] = 1.0;
ns++;

for(n=0;n<ns;n++)
   fprintf(stderr,"%13.5e %13.5e\n",fsp[n],asp[n]);

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

nt_p2 = getnt_p2(head1.nt);
s1 = (float *) check_realloc (s1,nt_p2*size_float);

taper_norm(s1,&head1.dt,head1.nt,&tap_per);
zero(s1+head1.nt,(nt_p2)-head1.nt);
forfft(s1,nt_p2,-1);

ampfac(s1,&head1.dt,nt_p2,fsp,asp,ns);

invfft(s1,nt_p2,1);
norm(s1,&head1.dt,nt_p2);

write_wccseis(outfile,&head1,s1,outbin);
}

ampfac(g,dt,n,fsp,asp,nsp)
struct complex *g;
float *dt, *fsp, *asp;
int n, nsp;
{
float fac, df, freq;
int i, j, jb, j1, nsp1;

nsp1 = nsp - 1;

jb = 0;
df = 1.0/(n*(*dt));
for(i=1;i<n/2;i++)
   {
   freq = i*df;

   if(freq < fsp[0])
      fac = 1.0;
   else if(freq >= fsp[nsp1])
      fac = 1.0;
   else
      {
      j = jb;
      while(freq > fsp[j])
         j++;

      jb = j;
      j1 = j - 1;

      fac = asp[j1] + (freq-fsp[j1])*(asp[j]-asp[j1])/(fsp[j]-fsp[j1]);
      }

   /*
   fprintf(stderr,"%13.5e %13.5e\n",freq,fac);
   */

   g[i].re = fac*g[i].re;
   g[i].im = fac*g[i].im;
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

getpeak(s,nt,pga)
float *s, *pga;
int nt;
{
int i;

for(i=0;i<nt;i++)
   {
   if(s[i] > *pga)
      *pga = s[i];
   if(-s[i] > *pga)
      *pga = -s[i];
   }

*pga = *pga/981.0;
}
