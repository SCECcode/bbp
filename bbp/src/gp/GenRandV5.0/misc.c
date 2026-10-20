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

void copy_psrc(struct pointsource *ps_orig,struct pointsource *ps_copy,int ns,int nd)
{
int i;

for(i=0;i<ns*nd;i++)
   {
   ps_copy[i].lon = ps_orig[i].lon;
   ps_copy[i].lat = ps_orig[i].lat;
   ps_copy[i].dep = ps_orig[i].dep;
   ps_copy[i].stk = ps_orig[i].stk;
   ps_copy[i].dip = ps_orig[i].dip;
   ps_copy[i].rak = ps_orig[i].rak;
   ps_copy[i].area = ps_orig[i].area;
   ps_copy[i].slip = ps_orig[i].slip;
   ps_copy[i].rupt = ps_orig[i].rupt;
   ps_copy[i].vs = ps_orig[i].vs;
   ps_copy[i].den = ps_orig[i].den;
   ps_copy[i].mu = ps_orig[i].mu;
   ps_copy[i].beta = ps_orig[i].beta;
   ps_copy[i].asp = ps_orig[i].asp;
   ps_copy[i].asp_mask = ps_orig[i].asp_mask;
   ps_copy[i].subevt = ps_orig[i].subevt;
   ps_copy[i].subevt_mask = ps_orig[i].subevt_mask;
   ps_copy[i].aseis = ps_orig[i].aseis;
   ps_copy[i].rvf = ps_orig[i].rvf;
   ps_copy[i].trise_fac = ps_orig[i].trise_fac;
   ps_copy[i].rvf_fac = ps_orig[i].rvf_fac;
   }
}

void get_mean_sigma_c(struct complex *x,int ns,int nd,float *avg,float *sig)
{
int i, j;

*avg = 0.0;
for(i=0;i<ns*nd;i++)
   *avg = *avg + x[i].re;

*avg = (*avg)/(ns*nd);

*sig = 0.0;
for(i=0;i<ns*nd;i++)
   *sig = *sig + (x[i].re - (*avg))*(x[i].re - (*avg));

*sig = sqrt((*sig)/(ns*nd));
}

void hyb_depth_sum(struct complex *x_shal,struct complex *x_deep,int ns2,int nd2,struct pointsource *ps,int ns,int nd,struct hybrid_correlation_length *hc)
{
float dmin, dmax, hc_shal_wt, hc_deep_wt;
int i, j, ip, ip2;

dmin = hc->dep - hc->dep_range;
dmax = hc->dep + hc->dep_range;

for(j=0;j<nd2;j++)
   {
   ip2 = j;
   if(ip2 > nd-1)
      ip2 = nd-1;

   if(ps[ip2*ns].dep <= dmin)
      {
      hc_shal_wt = hc->shal_wt_start;
      hc_deep_wt = hc->deep_wt_start;
      }
   else if(ps[ip2*ns].dep < dmax && ps[ip2*ns].dep > dmin)
      {
      hc_shal_wt = hc->shal_wt_start - (hc->shal_wt_start-hc->shal_wt_end)*((ps[ip2*ns].dep)-dmin)/(dmax-dmin);
      hc_deep_wt = 1.0 - hc_shal_wt;
      }
   else
      {
      hc_shal_wt = hc->shal_wt_end;
      hc_deep_wt = hc->deep_wt_end;
      }

   for(i=0;i<ns2;i++)
      {
      ip = i + j*ns2;

      x_deep[ip].re = hc_deep_wt*x_deep[ip].re + hc_shal_wt*x_shal[ip].re;
      x_deep[ip].im = hc_deep_wt*x_deep[ip].im + hc_shal_wt*x_shal[ip].im;
      }
   }
}

void hyb_depth_sumX(struct complex *x_shal,struct complex *x_deep,int ns2,int nd2,struct pointsource *ps,int ns,int nd,struct hybrid_correlation_length *hc)
{
float dmin1, dmax1, hc_shal_wt, hc_deep_wt;
int i, j, ip, ip2;

dmin1 = hc->dep - hc->dep_range;
dmax1 = hc->dep + hc->dep_range;

for(j=0;j<nd2;j++)
   {
   ip2 = j;
   if(ip2 > nd-1)
      ip2 = nd-1;

   if(ps[ip2*ns].dep <= dmin1)
      {
      hc_shal_wt = 1.0;
      hc_deep_wt = 0.0;
      }
   else if(ps[ip2*ns].dep < dmax1 && ps[ip2*ns].dep > dmin1)
      {
      hc_shal_wt = (dmax1-(ps[ip2*ns].dep))/(dmax1-dmin1);
      hc_deep_wt = 1.0 - hc_shal_wt;
      }
   else
      {
      hc_shal_wt = 0.0;
      hc_deep_wt = 1.0;
      }

   for(i=0;i<ns2;i++)
      {
      ip = i + j*ns2;

      x_deep[ip].re = hc_deep_wt*x_deep[ip].re + hc_shal_wt*x_shal[ip].re;
      x_deep[ip].im = hc_deep_wt*x_deep[ip].im + hc_shal_wt*x_shal[ip].im;
      }
   }
}

void cor_c(struct complex *x1,struct complex *x2,int ns2,int nd2,float *targ_cor)
{
int ip;
float sf;
float x1_avg_r, x1_avg_i, x1_sig_r, x1_sig_i;
float x2_avg_r, x2_avg_i, x2_sig_r, x2_sig_i;
float cor_r, cor_i;

sf = sqrt(1.0 - (*targ_cor)*(*targ_cor));

x1_avg_r = 0.0;
x1_avg_i = 0.0;
x2_avg_r = 0.0;
x2_avg_i = 0.0;
for(ip=0;ip<nd2*ns2;ip++)
   {
   x2[ip].re = (*targ_cor)*x1[ip].re + sf*x2[ip].re;
   x2[ip].im = (*targ_cor)*x1[ip].im + sf*x2[ip].im;

   x1_avg_r = x1_avg_r + x1[ip].re;
   x1_avg_i = x1_avg_i + x1[ip].im;
   x2_avg_r = x2_avg_r + x2[ip].re;
   x2_avg_i = x2_avg_i + x2[ip].im;
   }

x1_avg_r = x1_avg_r/(ns2*nd2);
x1_avg_i = x1_avg_i/(ns2*nd2);
x2_avg_r = x2_avg_r/(ns2*nd2);
x2_avg_i = x2_avg_i/(ns2*nd2);

x1_sig_r = 0.0;
x1_sig_i = 0.0;
x2_sig_r = 0.0;
x2_sig_i = 0.0;
for(ip=0;ip<nd2*ns2;ip++)
   {
   x1_sig_r = x1_sig_r + (x1[ip].re-x1_avg_r)*(x1[ip].re-x1_avg_r);
   x1_sig_i = x1_sig_i + (x1[ip].im-x1_avg_i)*(x1[ip].im-x1_avg_i);
   x2_sig_r = x2_sig_r + (x2[ip].re-x2_avg_r)*(x2[ip].re-x2_avg_r);
   x2_sig_i = x2_sig_i + (x2[ip].im-x2_avg_i)*(x2[ip].im-x2_avg_i);
   }

x1_sig_r = sqrt(x1_sig_r/(ns2*nd2));
x1_sig_i = sqrt(x1_sig_i/(ns2*nd2));
x2_sig_r = sqrt(x2_sig_r/(ns2*nd2));
x2_sig_i = sqrt(x2_sig_i/(ns2*nd2));

cor_r = 0.0;
cor_i = 0.0;
for(ip=0;ip<nd2*ns2;ip++)
   {
   cor_r = cor_r + (x1[ip].re-x1_avg_r)*(x2[ip].re-x2_avg_r);
   cor_i = cor_i + (x1[ip].im-x1_avg_i)*(x2[ip].im-x2_avg_i);
   }

cor_r = cor_r/(ns2*nd2*x1_sig_r*x2_sig_r);
cor_i = cor_i/(ns2*nd2*x1_sig_i*x2_sig_i);

fprintf(stderr,"x1_avg_r= %.5e x1_avg_i= %.5e\n",x1_avg_r,x1_avg_i);
fprintf(stderr,"x1_sig_r= %.5e x1_sig_i= %.5e\n",x1_sig_r,x1_sig_i);
fprintf(stderr,"x2_avg_r= %.5e x2_avg_i= %.5e\n",x2_avg_r,x2_avg_i);
fprintf(stderr,"x2_sig_r= %.5e x2_sig_i= %.5e\n\n",x2_sig_r,x2_sig_i);

fprintf(stderr,"target_cor= %.3f cor_r= %.3f cor_i= %.3f\n\n",*targ_cor,cor_r,cor_i);
}

void check_cor_r(float *x1,char *c1,float *x2,char *c2,int ns2,int nd2,float *targ_cor)
{
int ip, nn;
float x1_avg, x1_sig, x1_min, x1_max;
float x2_avg, x2_sig, x2_min, x2_max;
float cor;

float minv = 0.0;

/*
 * set floor for slip
minv = 0.001;
*/

x1_min = 1.0e+15;
x1_max = -1.0e+15;
x2_min = 1.0e+15;
x2_max = -1.0e+15;

x1_avg = 0.0;
x2_avg = 0.0;
nn = 0;
for(ip=0;ip<nd2*ns2;ip++)
   {
   if(x1[ip] >= minv)
      {
      nn++;

      x1_avg = x1_avg + x1[ip];
      x2_avg = x2_avg + x2[ip];

      if(x1[ip] < x1_min)
         x1_min = x1[ip];
      if(x1[ip] > x1_max)
         x1_max = x1[ip];

      if(x2[ip] < x2_min)
         x2_min = x2[ip];
      if(x2[ip] > x2_max)
         x2_max = x2[ip];
      }
   }

x1_avg = x1_avg/(nn);
x2_avg = x2_avg/(nn);

x1_sig = 0.0;
x2_sig = 0.0;
nn = 0;
for(ip=0;ip<nd2*ns2;ip++)
   {
   if(x1[ip] >= minv)
      {
      nn++;

      x1_sig = x1_sig + (x1[ip]-x1_avg)*(x1[ip]-x1_avg);
      x2_sig = x2_sig + (x2[ip]-x2_avg)*(x2[ip]-x2_avg);
      }
   }

x1_sig = sqrt(x1_sig/(nn));
x2_sig = sqrt(x2_sig/(nn));

cor = 0.0;
nn = 0;
for(ip=0;ip<nd2*ns2;ip++)
   {
   if(x1[ip] >= minv)
      {
      nn++;

      cor = cor + (x1[ip]-x1_avg)*(x2[ip]-x2_avg);
      }
   }

cor = cor/(nn*x1_sig*x2_sig);

fprintf(stderr,"\nChecking correlation between '%s' and '%s'-\n",c1,c2);
fprintf(stderr,"\t%s:\tavg= %.5e sig= %.5e min= %.5e max= %.5e\n",c1,x1_avg,x1_sig,x1_min,x1_max);
fprintf(stderr,"\t%s:\tavg= %.5e sig= %.5e min= %.5e max= %.5e\n",c2,x2_avg,x2_sig,x2_min,x2_max);

fprintf(stderr,"\ttarget cor= %.3f computed cor= %.3f nn/(ns2*nd2)= %.5e\n\n",*targ_cor,cor,(float)(nn)/((float)(ns2*nd2)));
}

void spline_natural(int nt,float *dt,float *y,float *y2,float *u)
{
int it;
float p, invdt;

invdt = 1.0/(*dt);

y2[0] = u[0] = 0.0;
for(it=1;it<nt-1;it++)
   {
   p = 0.5*y2[it-1] + 2.0;
   y2[it] = -0.5/p;

   u[it] = (y[it+1]+y[it-1])*invdt;
   u[it] = (3.0*u[it]*invdt - 0.5*u[it-1])/p;
   }

y2[nt-1] = 0.0;
for(it=nt-2;it>=0;it--)
   y2[it] = y2[it]*y2[it+1] + u[it];
}

void splint_vector(int nt,float *dt,float *y,float *y2,int nt_o,float *dt_o,float *y_o)
{
int it_b, it_f, it_o;
float tb_next, t_o, ht, at, bt, ht26, at3, bt3;

it_b = 0;
it_f = 1;
for(it_o=0;it_o<nt_o;it_o++)
   {
   t_o = it_o*(*dt_o);

   tb_next = (it_b+1)*(*dt);
   while(tb_next < t_o)
      {
      it_b++;
      tb_next = (it_b+1)*(*dt);
      }

   while(it_f*(*dt) <= t_o)
      it_f++;

   if(it_f - it_b <= 0)
      {
      fprintf(stderr,"problem in splint_vector, exiting...\n");
      exit(-1);
      }

   ht = (it_f - it_b)*(*dt);
   at = (it_f*(*dt) - t_o)/ht;
   bt = (t_o - it_b*(*dt))/ht;

   ht26 = ht*ht/6.0;
   at3 = at*at*at;
   bt3 = bt*bt*bt;
   y_o[it_o] = at*y[it_b] + bt*y[it_f] + ((at3-at)*y[it_b] + (bt3-bt)*y[it_f])*ht26;
   }
}
