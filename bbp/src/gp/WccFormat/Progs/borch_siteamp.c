#include "include.h"
#include "structure.h"
#include "function.h"

main(int ac,char **av)
{
float df, pga, vref, vsite, vrat;
float ma, mv, fa, fv, fac, freq;
int i;
int nf = 1000;

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

vrat = vref/vsite;

if(pga < 0.1)
   {
   ma = 0.35;
   mv = 0.65;
   }
else if(pga < 0.2)
   {
   ma = 0.35 - (pga - 0.1);
   mv = 0.65 - 0.5*(pga - 0.1);
   }
else if(pga < 0.3)
   {
   ma = 0.25 - 1.5*(pga - 0.2);
   mv = 0.60 - 0.7*(pga - 0.2);
   }
else if(pga < 0.4)
   {
   ma = 0.10 - 1.5*(pga - 0.3);
   mv = 0.53 - 0.8*(pga - 0.3);
   }
else
   {
   ma = -0.05;
   mv = 0.45;
   }

fa = exp(ma*log(vrat));
fv = exp(mv*log(vrat));

if(fmin > fmidbot)
   fmin = fmidbot;
if(fmax < fhightop)
   fmin = fhightop;

df = log(ftop/fbot)/(float)(nf-1);
for(i=0;i<nf;i++)
   {
   freq = fbot*exp(i*df);

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
   }
}
