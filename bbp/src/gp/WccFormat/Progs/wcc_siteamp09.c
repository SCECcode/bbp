#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

#define TAP_PERC 0.05

int size_float = sizeof(float);
int size_int = sizeof(int);

void borch_ampf(float *,float *,int,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void cb2008_ampf(float *,float *,int,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void ampfac(struct complex *,float *,int);
void zero(float *, int);
void norm(float *, float *, int);
void taper_norm(float *, float *, int, float *);
int getnt_p2(int);
void getpeak(float *, int, float *);

int main(int ac,char **av)
{
struct statdata head1;
float *s1, vref, vsite, vpga;
float *ampf;
int nt_p2;

float tap_per = TAP_PERC;
float pga = -1.0;

float fmin = 0.1;
float fmax = 15.0;

float flowcap = 0.0;   /* ampf for f<flowcap set equal to ampf[f=flowcap], (caps low-freq amplification level) */

char infile[256];
char outfile[256];
char model[128];

int inbin = 0;
int outbin = 0;

float fmidbot = 0.2;     /* bottom-end of middle frequency range */
float fmid = 1.0;        /* center of middle frequency range */
float fhigh = 3.333;     /* center of high frequency range */
float fhightop = 10.0;   /* top-end of high frequency range */

/*
sprintf(model,"borcherdt");
*/
sprintf(model,"cb2008");

setpar(ac,av);

mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
mstpar("vref","f",&vref);
mstpar("vsite","f",&vsite);

getpar("model","s",model);
getpar("pga","f",&pga);
vpga = vref;
getpar("vpga","f",&vpga);
getpar("flowcap","f",&flowcap);

getpar("tap_per","f",&tap_per);
getpar("fmin","f",&fmin);
getpar("fmidbot","f",&fmidbot);
getpar("fmid","f",&fmid);
getpar("fhigh","f",&fhigh);
getpar("fhightop","f",&fhightop);
getpar("fmax","f",&fmax);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);

endpar();

if(strncmp(model,"borcherdt",9) != 0)
   sprintf(model,"cb2008");

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

nt_p2 = getnt_p2(head1.nt);
s1 = (float *) check_realloc (s1,nt_p2*size_float);

ampf = (float *) check_malloc ((nt_p2/2)*size_float);

if(pga < 0.0)
   getpeak(s1,head1.nt,&pga);
else
   fprintf(stderr,"*** External PGA used: ");

fprintf(stderr,"pga= %13.5e\n",pga);

taper_norm(s1,&head1.dt,head1.nt,&tap_per);
zero(s1+head1.nt,(nt_p2)-head1.nt);
forfft((struct complex *)s1,nt_p2,-1);

if(strncmp(model,"cb2008",6) == 0)
   cb2008_ampf(ampf,&head1.dt,nt_p2,&vref,&vsite,&vpga,&pga,&fmin,&fmidbot,&fmid,&fhigh,&fhightop,&fmax,&flowcap);
else
   borch_ampf(ampf,&head1.dt,nt_p2,&vref,&vsite,&pga,&fmin,&fmidbot,&fmid,&fhigh,&fhightop,&fmax);

ampfac((struct complex *)s1,ampf,nt_p2);

invfft((struct complex *)s1,nt_p2,1);
norm(s1,&head1.dt,nt_p2);

write_wccseis(outfile,&head1,s1,outbin);
}

void borch_ampf(float *ampf,float *dt,int n,float *vref,float *vsite,float *pga,float *fmin,float *fmidbot,float *fmid,float *fhigh,float *fhightop,float *fmax)
{
float ma, mv, fa, fv, fac, df, freq;
float vr;
int i;

vr = (*vref)/(*vsite);

if(*pga < 0.1)
   {
   ma = 0.35;
   mv = 0.65;
   }
else if(*pga < 0.2)
   {
   ma = 0.35 - (*pga - 0.1);
   mv = 0.65 - 0.5*(*pga - 0.1);
   }
else if(*pga < 0.3)
   {
   ma = 0.25 - 1.5*(*pga - 0.2);
   mv = 0.60 - 0.7*(*pga - 0.2);
   }
else if(*pga < 0.4)
   {
   ma = 0.10 - 1.5*(*pga - 0.3);
   mv = 0.53 - 0.8*(*pga - 0.3);
   }
else
   {
   ma = -0.05;
   mv = 0.45;
   }

fa = exp(ma*log(vr));
fv = exp(mv*log(vr));

if(*fmin > *fmidbot)
   *fmin = *fmidbot;
if(*fmax < *fhightop)
   *fmin = *fhightop;

df = 1.0/(n*(*dt));
for(i=1;i<n/2;i++)
   {
   freq = i*df;

   if(freq < *fmin)
      fac = 1.0;

   else if(freq < *fmidbot)
      fac = 1.0 + (freq - *fmin)*(fv - 1.0)/(*fmidbot - *fmin);

   else if(freq < *fmid)
      fac = fv;

   else if(freq < *fhigh)
      fac = fv + (freq - *fmid)*(fa - fv)/(*fhigh - *fmid);

   else if(freq < *fhightop)
      fac = fa;

   else if(freq < *fmax)
      fac = fa + (freq - *fhightop)*(1.0 - fa)/(*fmax - *fhightop);

   else
      fac = 1.0;

   ampf[i] = fac;
   }
}

void ampfac(struct complex *g,float *ampf,int n)
{
int i;

for(i=1;i<n/2;i++)
   {
   g[i].re = ampf[i]*g[i].re;
   g[i].im = ampf[i]*g[i].im;
   }
}

void norm(g,dt,nt)
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

void zero(s,n)
float *s;
int n;
{
while(n--)
   {
   s[0] = 0.0;
   s++;
   }
}

void taper_norm(g,dt,nt,tap_per)
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

int getnt_p2(nt)
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

void getpeak(s,nt,pga)
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

void cb2008_ampf(float *ampf,float *dt,int n,float *vref,float *vsite,float *vpga,float *pga,float *fmin,float *fmidbot,float *fmid,float *fhigh,float *fhightop,float *fmax,float *flowcap)
{
float scon_c, scon_n, per[22], c10[22], c11[22], k1[22], k2[22], k3[22];
float ampf0[22];
float a_1100, fs_1100, fs_vpga, fs_vref, fsite;
float df, ampv, afac, freq;
float f0, f1, a0, a1, dadf, ampf_cap;
int i, j;
int nper = 22;

scon_c = 1.88;
scon_n = 1.18;

per[ 0] = 0.0;
per[ 1] = 0.01;
per[ 2] = 0.02;
per[ 3] = 0.03;
per[ 4] = 0.05;
per[ 5] = 0.075;
per[ 6] = 0.10;
per[ 7] = 0.15;
per[ 8] = 0.20;
per[ 9] = 0.25;
per[10] = 0.30;
per[11] = 0.40;
per[12] = 0.50;
per[13] = 0.75;
per[14] = 1.00;
per[15] = 1.50;
per[16] = 2.00;
per[17] = 3.00;
per[18] = 4.00;
per[19] = 5.00;
per[20] = 7.50;
per[21] = 10.0;

c10[ 0] = 1.058;
c10[ 1] = 1.058;
c10[ 2] = 1.102;
c10[ 3] = 1.174;
c10[ 4] = 1.272;
c10[ 5] = 1.438;
c10[ 6] = 1.604;
c10[ 7] = 1.928;
c10[ 8] = 2.194;
c10[ 9] = 2.351;
c10[10] = 2.46;
c10[11] = 2.587;
c10[12] = 2.544;
c10[13] = 2.133;
c10[14] = 1.571;
c10[15] = 0.406;
c10[16] = -0.456;
c10[17] = -0.82;
c10[18] = -0.82;
c10[19] = -0.82;
c10[20] = -0.82;
c10[21] = -0.82;

c11[ 0] = 0.04;
c11[ 1] = 0.04;
c11[ 2] = 0.04;
c11[ 3] = 0.04;
c11[ 4] = 0.04;
c11[ 5] = 0.04;
c11[ 6] = 0.04;
c11[ 7] = 0.04;
c11[ 8] = 0.04;
c11[ 9] = 0.04;
c11[10] = 0.04;
c11[11] = 0.04;
c11[12] = 0.04;
c11[13] = 0.077;
c11[14] = 0.15;
c11[15] = 0.253;
c11[16] = 0.3;
c11[17] = 0.3;
c11[18] = 0.3;
c11[19] = 0.3;
c11[20] = 0.3;
c11[21] = 0.3;

k1[ 0] = 865.0;
k1[ 1] = 865.0;
k1[ 2] = 865.0;
k1[ 3] = 908.0;
k1[ 4] = 1054.0;
k1[ 5] = 1086.0;
k1[ 6] = 1032.0;
k1[ 7] = 878.0;
k1[ 8] = 748.0;
k1[ 9] = 654.0;
k1[10] = 587.0;
k1[11] = 503.0;
k1[12] = 457.0;
k1[13] = 410.0;
k1[14] = 400.0;
k1[15] = 400.0;
k1[16] = 400.0;
k1[17] = 400.0;
k1[18] = 400.0;
k1[19] = 400.0;
k1[20] = 400.0;
k1[21] = 400.0;

k2[ 0] = -1.186;
k2[ 1] = -1.186;
k2[ 2] = -1.219;
k2[ 3] = -1.273;
k2[ 4] = -1.346;
k2[ 5] = -1.471;
k2[ 6] = -1.624;
k2[ 7] = -1.931;
k2[ 8] = -2.188;
k2[ 9] = -2.381;
k2[10] = -2.518;
k2[11] = -2.657;
k2[12] = -2.669;
k2[13] = -2.401;
k2[14] = -1.955;
k2[15] = -1.025;
k2[16] = -0.299;
k2[17] = 0.0;
k2[18] = 0.0;
k2[19] = 0.0;
k2[20] = 0.0;
k2[21] = 0.0;

k3[ 0] = 1.839;
k3[ 1] = 1.839;
k3[ 2] = 1.84;
k3[ 3] = 1.841;
k3[ 4] = 1.843;
k3[ 5] = 1.845;
k3[ 6] = 1.847;
k3[ 7] = 1.852;
k3[ 8] = 1.856;
k3[ 9] = 1.861;
k3[10] = 1.865;
k3[11] = 1.874;
k3[12] = 1.883;
k3[13] = 1.906;
k3[14] = 1.929;
k3[15] = 1.974;
k3[16] = 2.019;
k3[17] = 2.11;
k3[18] = 2.2;
k3[19] = 2.291;
k3[20] = 2.517;
k3[21] = 2.744;

fs_1100 = (c10[0] + k2[0]*scon_n)*log(1100.0/k1[0]);

if((*vpga) < k1[0])   /* 'pga' should really be 'a_1100' below, but this is unknown */
   {
   fs_vpga = c10[0]*log((*vpga)/k1[0]) +
             k2[0]*(log(((*pga) + scon_c*exp(scon_n*log((*vpga)/k1[0])))/((*pga) + scon_c)));
   }
else if((*vpga) < 1100.0)
   fs_vpga = (c10[0] + k2[0]*scon_n)*log((*vpga)/k1[0]);
else
   fs_vpga = (c10[0] + k2[0]*scon_n)*log(1100.0/k1[0]);

a_1100 = (*pga)*exp(fs_1100 - fs_vpga);

ampf_cap = -1.0;
for(i=0;i<nper;i++)
   {
   if((*vsite) < k1[i])
      {
      fsite = c10[i]*log((*vsite)/k1[i]) +
                k2[i]*(log((a_1100 + scon_c*exp(scon_n*log((*vsite)/k1[i])))/(a_1100 + scon_c)));
      }
   else if((*vsite) < 1100.0)
      fsite = (c10[i] + k2[i]*scon_n)*log((*vsite)/k1[i]);
   else
      fsite = (c10[i] + k2[i]*scon_n)*log(1100.0/k1[i]);

   if((*vref) < k1[i])
      {
      fs_vref = c10[i]*log((*vref)/k1[i]) +
                k2[i]*(log((a_1100 + scon_c*exp(scon_n*log((*vref)/k1[i])))/(a_1100 + scon_c)));
      }
   else if((*vref) < 1100.0)
      fs_vref = (c10[i] + k2[i]*scon_n)*log((*vref)/k1[i]);
   else
      fs_vref = (c10[i] + k2[i]*scon_n)*log(1100.0/k1[i]);

   ampf0[i] = exp(fsite - fs_vref);

   if(1.0/per[i] <= (*flowcap))
      {
      if(ampf_cap < 0.0)
         ampf_cap = ampf0[i];
      else
         ampf0[i] = ampf_cap;
      }
   }

   /* go in reverse order so frequencies are increasing */

j = nper - 1;
f0 = 1.0/per[j];
a0 = ampf0[j];
f1 = 1.0/per[j];
a1 = ampf0[j];
dadf = 0.0;

df = 1.0/(n*(*dt));
for(i=1;i<n/2;i++)
   {
   freq = i*df;

   if(freq > f1)
      {
      f0 = f1;
      a0 = a1;

      if(j > 0)
         j--;

      if(per[j] != 0.0)
         f1 = 1.0/per[j];
      else
         f1 = 1000.0;

      a1 = ampf0[j];

      if(f1 != f0)
	 /*
         dadf = (a1-a0)/(f1-f0);
	 */
         dadf = (a1-a0)/log(f1/f0);
      else
         dadf = 0.0;
      }

   /*
   ampv = a0 + dadf*(freq-f0);
   */
   ampv = a0 + dadf*log(freq/f0);

   if(freq < *fmin)
      afac = 1.0;

   else if(freq < *fmidbot)
      /*
      afac = 1.0 + (freq - *fmin)*(ampv - 1.0)/(*fmidbot - *fmin);
      */
      afac = 1.0 + log(freq/(*fmin))*(ampv - 1.0)/log((*fmidbot)/(*fmin));

   else if(freq < *fmid)
      afac = ampv;

   else if(freq < *fhigh)
      afac = ampv;

   else if(freq < *fhightop)
      afac = ampv;

   else if(freq < *fmax)
      /*
      afac = ampv + (freq - *fhightop)*(1.0 - ampv)/(*fmax - *fhightop);
      */
      afac = ampv + log(freq/(*fhightop))*(1.0 - ampv)/log((*fmax)/(*fhightop));

   else
      afac = 1.0;

   ampf[i] = afac;
   }
}
