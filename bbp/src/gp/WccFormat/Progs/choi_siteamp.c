#include "include.h"
#include "structure.h"
#include "function.h"

#define NO_PER_A3 15

struct coeffs
   {
   float per;
   float b;
   float b1;
   float b2;
   float c;
   float vref;
   };

void init_modA3(int,struct coeffs *,float *,float *);

main(int ac,char **av)
{
struct coeffs *modA3;
float df, pga, vref, vsite, vrat;
float fac, freq;
int i;
int nf = 1000;

int nA3 = NO_PER_A3;
float bv = 300.0;

float fmidbot = 0.3;     /* bottom-end of middle frequency range */
float fmid = 1.0;        /* center of middle frequency range */
float fhigh = 3.333;     /* center of high frequency range */
float fhightop = 10.0;   /* top-end of high frequency range */

float fbot = 0.02;
float ftop = 50.0;

float fmin = 0.1;
float fmax = 15.0;

setpar(ac,av);
mstpar("vref","f",&vref);
mstpar("vsite","f",&vsite);
mstpar("pga","f",&pga);
getpar("fmin","f",&fmin);
getpar("fmidbot","f",&fmidbot);
getpar("fmid","f",&fmid);
getpar("fhigh","f",&fhigh);
getpar("fhightop","f",&fhightop);
getpar("fmax","f",&fmax);
getpar("nf","d",&nf);
endpar();

modA3 = (struct coeffs *)check_malloc(nA3*sizeof(struct coeffs));

init_modA3(nA3,modA3,&vsite,&bv);

vrat = vref/vsite;

if(fmin > fmidbot)
   fmin = fmidbot;
if(fmax < fhightop)
   fmin = fhightop;

df = log(ftop/fbot)/(float)(nf-1);
for(i=0;i<nf;i++)
   {
   freq = fbot*exp(i*df);

   /*
   if(freq < fmin)
      fac = 1.0;

   else if(freq < fmidbot)
      fac = 1.0 + (freq - fmin)*(fv - 1.0)/(fmidbot - fmin);

   else if(freq < fmid)
      fac = fv;

   else if(freq < fhigh)
      fac = fv + (freq - fmid)*(fa - fv)/(fhigh - fmid);

   else if(freq < fhightop)
      fac = fa;

   else if(freq < fmax)
      fac = fa + (freq - fhightop)*(1.0 - fa)/(fmax - fhightop);

   else
      fac = 1.0;

   printf("%13.5e %10.5f\n",freq,fac);
   */
   }

vrat = vsite/vref;
for(i=0;i<nA3;i++)
   {
   freq = 1.0/modA3[i].per;

   vrat = vsite/modA3[i].vref;
   fac = modA3[i].c*log(vrat) + modA3[i].b*log(10.0*pga);

   printf("%13.5e %10.5f\n",freq,exp(fac));
   }
}

void init_modA3(int n,struct coeffs *mod,float *vs,float *bv)
{
float f1, f2;
int i;

mod[ 0].per = 0.01;
mod[ 1].per = 0.05;
mod[ 2].per = 0.075;
mod[ 3].per = 0.10;
mod[ 4].per = 0.15;
mod[ 5].per = 0.20;
mod[ 6].per = 0.30;
mod[ 7].per = 0.40;
mod[ 8].per = 0.50;
mod[ 9].per = 0.75;
mod[10].per = 1.00;
mod[11].per = 1.50;
mod[12].per = 2.00;
mod[13].per = 3.00;
mod[14].per = 4.00;

mod[ 0].b1 = -0.55;
mod[ 1].b1 = -0.57;
mod[ 2].b1 = -0.61;
mod[ 3].b1 = -0.57;
mod[ 4].b1 = -0.52;
mod[ 5].b1 = -0.51;
mod[ 6].b1 = -0.51;
mod[ 7].b1 = -0.50;
mod[ 8].b1 = -0.50;
mod[ 9].b1 = -0.49;
mod[10].b1 = -0.49;
mod[11].b1 = -0.48;
mod[12].b1 = -0.46;
mod[13].b1 = -0.42;
mod[14].b1 = -0.40;

mod[ 0].b2 = -0.04;
mod[ 1].b2 = -0.05;
mod[ 2].b2 = -0.11;
mod[ 3].b2 = -0.12;
mod[ 4].b2 = -0.13;
mod[ 5].b2 = -0.07;
mod[ 6].b2 = -0.04;
mod[ 7].b2 = -0.02;
mod[ 8].b2 = -0.02;
mod[ 9].b2 = -0.02;
mod[10].b2 = -0.04;
mod[11].b2 = -0.12;
mod[12].b2 = -0.17;
mod[13].b2 = -0.22;
mod[14].b2 = -0.25;

mod[ 0].c = -0.34;
mod[ 1].c = -0.26;
mod[ 2].c = -0.21;
mod[ 3].c = -0.22;
mod[ 4].c = -0.24;
mod[ 5].c = -0.28;
mod[ 6].c = -0.41;
mod[ 7].c = -0.50;
mod[ 8].c = -0.59;
mod[ 9].c = -0.65;
mod[10].c = -0.68;
mod[11].c = -0.71;
mod[12].c = -0.72;
mod[13].c = -0.72;
mod[14].c = -0.72;

mod[ 0].vref = 501.0;
mod[ 1].vref = 676.0;
mod[ 2].vref = 780.0;
mod[ 3].vref = 643.0;
mod[ 4].vref = 541.0;
mod[ 5].vref = 565.0;
mod[ 6].vref = 610.0;
mod[ 7].vref = 640.0;
mod[ 8].vref = 660.0;
mod[ 9].vref = 703.0;
mod[10].vref = 709.0;
mod[11].vref = 710.0;
mod[12].vref = 710.0;
mod[13].vref = 710.0;
mod[14].vref = 710.0;

f1 = ((*vs) - (*bv))/(180.0 - (*bv));
f1 = f1*f1;
f2 = ((*vs) - 520.0)/240.0;

for(i=0;i<n;i++)
   {
   if((*vs) <= 180.0)
      mod[i].b = mod[i].b1;
   else if((*vs) <= (*bv))
      mod[i].b = mod[i].b2 + f1*(mod[i].b1-mod[i].b2);
   else if((*vs) <= 520.0)
      mod[i].b = mod[i].b2;
   else if((*vs) <= 760.0)
      mod[i].b = mod[i].b2 - f2*mod[i].b2;
   else
      mod[i].b = 0.0;
   }
}
