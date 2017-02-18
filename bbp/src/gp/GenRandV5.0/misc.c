#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void *check_malloc(size_t len)
{
void *ptr;

ptr = (void *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }

return(ptr);
}

void *check_realloc(void *ptr,size_t len)
{
ptr = (char *) realloc (ptr,len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   exit(-1);
   }

return(ptr);
}

static  long    frandx = 1;

/* frand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double frand(void)
{
frandx = (frandx * 1103515245 + 12345) & 0x7fffffff;
return((double)(frandx)/1073741824.0 - 1.0);
}

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}

double gaus_rand(float *sigma,float *mean,long *seed)
{
double r = 0.0;
double six = 6.0;
double one = 1.0;
double half = 0.5;
int i;

for(i=0;i<12;i++)
   r = r + half*(one + sfrand(seed));

return((double)((r - six)*(*sigma) + *mean));
}

void zapit(s,n)
float *s;
int n;
{
while(n--)
   {
   s[0] = 0.0;
   s++;
   }
}

void set_ne(float *elon,float *elat,float *slon,float *slat,float *sn,float *se)
{
float kperd_n, kperd_e;
double e2, den, g2, lat0;
float cosA, sinA;

double rperd = 0.017453293;
double radius = 6378.139;
double f = 298.256;

f = 1.0/f;
e2 = 2.0*f - f*f;
g2 = e2/((1.0 - f)*(1.0 - f));

lat0 = atan((1.0 - f)*tan((*elat)*rperd));

cosA = cos(lat0);
sinA = sin(lat0);

den = sqrt(1.0/(1.0 + g2*sinA*sinA));
kperd_e = rperd*radius*cosA*den;
kperd_n = rperd*radius*(sqrt(1.0 + g2*sinA*sinA*(2.0 + g2)))*den*den*den;

*sn = ((*slat) - (*elat))*kperd_n;
*se = ((*slon) - (*elon))*kperd_e;
}

void set_ll(float *elon,float *elat,float *slon,float *slat,float *sn,float *se)
{
float kperd_n, kperd_e;
double e2, den, g2, lat0;
float cosA, sinA;

double rperd = 0.017453293;
double radius = 6378.139;
double f = 298.256;

f = 1.0/f;
e2 = 2.0*f - f*f;
g2 = e2/((1.0 - f)*(1.0 - f));

lat0 = atan((1.0 - f)*tan((*elat)*rperd));

cosA = cos(lat0);
sinA = sin(lat0);

den = sqrt(1.0/(1.0 + g2*sinA*sinA));
kperd_e = rperd*radius*cosA*den;
kperd_n = rperd*radius*(sqrt(1.0 + g2*sinA*sinA*(2.0 + g2)))*den*den*den;

*slat = (*sn)/kperd_n + *elat;
*slon = (*se)/kperd_e + *elon;
}

void swap_in_place(int n,char *cbuf)
{
char cv;

while(n--)
   {
   cv = cbuf[0];
   cbuf[0] = cbuf[3];
   cbuf[3] = cv;

   cv = cbuf[1];
   cbuf[1] = cbuf[2];
   cbuf[2] = cv;

   cbuf = cbuf + 4;
   }
}

void gseed(long *seed,struct pointsource *psrc,float *flen,float *fwid,float *dtop,float *mag)
{
*seed = (long)((1000.0*sqrt(psrc[0].lon*psrc[0].lon + psrc[0].lat*psrc[0].lat))
              + (100.0*((*flen) + (*fwid) + (*dtop) + (*mag)))
              + (10.0*sqrt(psrc[0].stk*psrc[0].stk + psrc[0].dip*psrc[0].dip)));
}

void rhypo_uniform(long *seed,float *shyp,float *dhyp,float *flen,float *fwid)
{
/* Simple uniform distribution with no restricted areas */
*shyp = 0.5*(*flen)*sfrand(seed);
*dhyp = 0.5*(*fwid)*(1.0 + sfrand(seed));
}

float rhypo1_lintaper(long *seed,struct hypo_distr_params *hpar)
{
float x0, f0, x1, f1, bigL;
float a0, b0, a1, b1, bigC, bigD, d1;
float xr, y0, y1, yr;

x0 = hpar->x0;
f0 = hpar->f0;
x1 = hpar->x1;
f1 = hpar->f1;
bigL = hpar->xlen;

a0 = (1.0 - f0)/x0;
b0 = f0;
a1 = -(1.0 - f1)/(bigL-x1);
b1 = 1.0 - x1*a1;

bigC = 1.0/(0.5*a0*x0*x0 + b0*x0 + (x1-x0) + 0.5*a1*(bigL*bigL - x1*x1) + b1*(bigL-x1));

bigD = 0.5*a0*x0*x0 + b0*x0 - x0;
d1 = x1 + bigD - (0.5*a1*x1*x1 + b1*x1);

y0 = bigC*(x0+bigD);
y1 = bigC*(x1+bigD);

yr = 0.5*(1.0 + sfrand(seed));

if(yr < y0)
   xr = (sqrt(b0*b0 + 2.0*a0*yr/bigC) - b0)/a0;
else if(yr > y1)
   xr = (sqrt(b1*b1 - 2.0*a1*(d1-yr/bigC)) - b1)/a1;
else
   xr = yr/bigC - bigD;

xr = xr + hpar->xshift;
return(xr);
}
