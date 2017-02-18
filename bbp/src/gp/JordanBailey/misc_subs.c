#include "include.h"
#include "structure.h"
#include "function.h"

extern void fourg_(float *, int *, int*, float *);

void get_radazi(float *az,float *rg,float *de,float *dn,float *x,float *y,float *cs,float *ss,float *e,float *n)
{
float fe, fn;

fe = (*y)*(*cs) + (*x)*(*ss);
fn = -(*y)*(*ss) + (*x)*(*cs);

*de = *e - fe;
*dn = *n - fn;

*az = atan2((*de),(*dn));
*rg = sqrt((*de)*(*de) + (*dn)*(*dn));
}

void get_master_list(struct sgtparams *sgtp,int np,long long *mindx,int *nm)
{
int ip, ig, im, ifnd, sflag, mcnt;
long long ll_int;

mcnt = 0;

for(ip=0;ip<np;ip++)
   {
   for(ig=0;ig<sgtp[ip].nsgt;ig++)
      {
      ifnd = 0;
      for(im=0;im<mcnt;im++)
         {
	 if(sgtp[ip].indx[ig] == mindx[im])
	    {
	    ifnd = 1;
	    break;
	    }
	 }

      if(ifnd == 0)
         {
         mindx[mcnt] = sgtp[ip].indx[ig];
         mcnt++;
         }
      }
   }

sflag = 1;
while(sflag)
   {
   sflag = 0;
   for(im=0;im<mcnt-1;im++)
      {
      if(mindx[im] > mindx[im+1])
         {
         ll_int = mindx[im+1]; mindx[im+1] = mindx[im]; mindx[im] = ll_int;
         sflag = 1;
         }
      }
   }

for(ip=0;ip<np;ip++)
   {
   for(ig=0;ig<sgtp[ip].nsgt;ig++)
      {
      for(im=0;im<mcnt;im++)
         {
	 if(sgtp[ip].indx[ig] == mindx[im])
	    {
	    sgtp[ip].master_ip[ig] = im;
	    break;
	    }
	 }
      }
   }

*nm = mcnt;
}

void get_indx(float *lon,float *lat,float *dep,struct sgtindex *indx,struct geoprojection *gp)
{
float xs, ys, xr, yr;
double invh;

invh = 1.0/(double)(indx->h);

if(gp->geoproj == 0)
   {
   xs = ((*lon) - gp->modellon)*gp->kmlon;
   ys = (gp->modellat - (*lat))*gp->kmlat;

   xr = xs*gp->cosR + ys*gp->sinR - gp->xshift;
   yr = -xs*gp->sinR + ys*gp->cosR - gp->yshift;
   }
else if(gp->geoproj == 1)
   gcproj(&xr,&yr,lon,lat,&gp->erad,&gp->g0,&gp->b0,gp->amat,gp->ainv,1);

indx->xsgt = (int)((double)(xr)*invh + 0.5);
indx->ysgt = (int)((double)(yr)*invh + 0.5);
indx->zsgt = (int)((double)((*dep))*invh + 1.5);
indx->indx = (long long)(indx->xsgt)*(long long)(100000000) + (long long)(indx->ysgt)*(long long)(10000) + (long long)(indx->zsgt);
}

void get_ard_srf(struct standrupformat *srf,int off,int ip,float *az,float *rg,float *z0,float *de,float *dn,float *slon,float *slat,struct geoprojection *gp)
{
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat, xs, ys;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

elon = apval_ptr[ip].lon;
elat = apval_ptr[ip].lat;

if(gp->geoproj == 0)
   set_ne(&elon,&elat,slon,slat,dn,de);
else if(gp->geoproj == 1)
   {
   gcproj(&xs,&ys,&elon,&elat,&gp->erad,&gp->g0,&gp->b0,gp->amat,gp->ainv,1);
   *dn = -(gp->xshift) - xs;
   *de = -(gp->yshift) - ys;
   }

*az = atan2((*de),(*dn));
*rg = sqrt((*de)*(*de) + (*dn)*(*dn));
*z0 = apval_ptr[ip].dep;
}

void resample(s,nt,dt,isamp,ntpad,ntrsmp,newdt,p)
float *s, *p, *dt, *newdt;
int nt, isamp, ntpad, ntrsmp;
{
float df, f, f0, fl, fl2, fac;
int i, j;

int ord = 4;
float one = 1.0;

int minus = -1;
int plus = 1;

taper_norm(s,dt,nt);
zapit(s+nt,(ntpad)-(nt));

for(i=ntpad-1;i>=0;i--)
   {
   s[2*i] = s[i];
   s[2*i + 1] = 0.0;
   }

fourg_(s,&ntpad,&minus,p);

if(isamp > 0)
   zapit(s+ntpad,2*ntrsmp-ntpad);
else if(isamp < 0)
   {
   if(ord)  /* lowpass at 80 % of new Nyquist */
      {
      f0 = 0.8/(2.0*(*newdt));
      df = 1.0/(ntrsmp*(*newdt));
      for(i=1;i<ntrsmp/2;i++)
         {
         f = i*df;

         fl = f/f0;
         fl2 = fl*fl;
         fl = fl2;
         for(j=1;j<ord;j++)
            fl = fl*fl2;

         fac = one/(one + fl);

         s[2*i] = fac*s[2*i];
         s[2*i + 1] = fac*s[2*i + 1];
         }
      }

   s[ntrsmp] = s[ntrsmp+1] = 0.0; /* zero nyquist */
   }

for(i=1;i<ntrsmp/2;i++)  /* replicate with complex-conjugate */
   {
   s[2*(ntrsmp-i)] = s[2*i];
   s[2*(ntrsmp-i) + 1] = -s[2*i + 1];
   }

fourg_(s,&ntrsmp,&plus,p);

for(i=0;i<ntrsmp;i++)
   s[i] = s[2*i];

norm(s,newdt,ntrsmp);
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

void taper_norm(g,dt,nt)
float *g, *dt;
int nt;
{
float fac, df, arg;
int i;
int ntap;

ntap = nt*TAP_PERC;

for(i=0;i<nt-ntap;i++)
   g[i] = g[i]*(*dt);

df = 3.14159/(float)(ntap);
for(i=nt-ntap;i<nt;i++)
   {
   arg = (i-(nt-(ntap+1)))*df;
   fac = (*dt)*0.5*(1.0 + cos(arg));
   g[i] = g[i]*fac;
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

float sin_table[] =
   {
	1.0000000e+00,	/* sin(pi/2) */
	7.0710678e-01,	/* sin(pi/4) */
	3.8268343e-01,	/* sin(pi/8) */
	1.9509032e-01,	/* sin(pi/16) */
	9.8017140e-02,	/* sin(pi/32) */
	4.9067674e-02,	/* sin(pi/64) */
	2.4541228e-02,	/* sin(pi/128) */
	1.2271538e-02,	/* sin(pi/256) */
	6.1358846e-03,	/* sin(pi/512) */
	3.0679568e-03,	/* sin(pi/1024) */
	1.5339802e-03,	/* sin(pi/2048) */
	7.6699032e-04,	/* sin(pi/4096) */
	3.8349519e-04,	/* sin(pi/8192) */
	1.9174760e-04,	/* sin(pi/16384) */
	9.5873799e-05, /* sin(pi/32768) */
	4.7936899e-05, /* sin(pi/65536) */
	2.3968499e-05, /* sin(pi/131072) */
	1.1984224e-05, /* sin(pi/262144) */
	5.9921125e-06, /* sin(pi/524288) */
	2.9960562e-06, /* sin(pi/1048576) */
	1.4980281e-06  /* sin(pi/2097152) */
   };

void forfft(x,n,isign)
register struct complex *x;
int n, isign;
   {
	register struct complex *px, *rx;
	float cn, sn, cd, sd, arg;
	float are, aim, bre, bim, real, imag;
	float *psintab;
	float half = 0.5;
	extern float sin_table[];
	int k;

	cfft_r(x,n/2,isign);

	/* do DC and Nyquist */
	real= x[0].re;
	imag= x[0].im;
	x[0].re= real + imag;
	x[0].im= real - imag;

	/* set up for sine recurrsion */
	psintab= sin_table;
	for(k=4; k<n; k <<= 1) psintab++;
	sd= *psintab++;
	real= *psintab;
	cd= 2.0 * real * real;

	sn= 0.0;
	cn= 1.0;
	if(isign < 0) sd= -sd;
	px= x + 1;
	rx= x + n/2 -1;
	while( px <= rx )
	   {
		real = cd*cn + sd*sn;
		imag = sd*cn - cd*sn;
		cn -= real;
		sn += imag;

		are= half*(px->re + rx->re);
		aim= half*(px->im - rx->im);
		bre= half*(px->im + rx->im);
		bim= half*(rx->re - px->re);

		real= bre*cn - bim*sn;
		imag= bre*sn + bim*cn;

		px->re = are + real;
		px->im = aim + imag;
		rx->re = are - real;
		rx->im = imag - aim;

		px++;
		rx--;
	   }
	if(abs(isign) > 1)
	   {
		x[n/2].re= x[0].im;
		x[0].im= x[n/2].im = 0.0;
	   }
   }

void invfft(x,n,isign)
register struct complex *x;
int n, isign;
   {
	register struct complex *px, *rx;
	float cn, sn, cd, sd;
	float are, aim, bre, bim, real, imag;
	float *psintab;
	extern float sin_table[];
	int k;

	if(abs(isign) > 1) x[0].im= x[n/2].re;

	/* do DC and Nyquist */
	real= x[0].re;
	imag= x[0].im;
	x[0].re= real + imag;
	x[0].im= real - imag;

	/* set up for sine recurrsion */
	psintab= sin_table;
	for(k=4; k<n; k <<= 1) psintab++;
	sd= *psintab++;
	real= *psintab;
	cd= 2.0 * real * real;

	sn= 0.0;
	cn= 1.0;
	if(isign < 0) sd= -sd;
	px= x + 1;
	rx= x + n/2 -1;
	while( px <= rx )
	   {
		real = cd*cn + sd*sn;
		imag = sd*cn - cd*sn;
		cn -= real;
		sn += imag;

		are= (px->re + rx->re);
		aim= (px->im - rx->im);
		bre= (px->re - rx->re);
		bim= (px->im + rx->im);

		real= bre*cn - bim*sn;
		imag= bre*sn + bim*cn;

		px->re = are - imag;
		px->im = aim + real;
		rx->re = are + imag;
		rx->im = real - aim;

		px++;
		rx--;
	   }
	cfft_r(x,n/2,isign);
   }

void cfft_r(x,n,isign)
struct complex *x;
int n,isign;
   {
	register struct complex *px, *qx, *rx;
	struct complex *limit, *qlimit, dtemp;
	float cn, sn, cd, sd, temp, real, imag;
	int m, j, istep;
	float *psintab;
	extern float sin_table[];

	limit= x + n;
	j= 0;
	for(px=x; px<limit; px++)
	   {
		if(px < (qx= x+j))
		   {	dtemp= *qx; *qx= *px; *px= dtemp;	}
		m = n>>1;
		while( m>=1 && j>=m )
		   { j-= m; m>>= 1;    }
		j+= m;
	   }
	rx= x+1;
	for(px=x; px<limit; px+= 2, rx+= 2)
	   {
		temp= rx->re;
		rx->re= px->re -temp;
		px->re += temp;
		temp= rx->im;
		rx->im= px->im -temp;
		px->im += temp;
	   }
	j=2;
	psintab= sin_table;
	while( j < n )
	   {
		istep= j<<1;
		sd= *psintab++;
		temp= *psintab;
		cd= 2.0 * temp * temp;
		cn= 1.0;
		sn= 0.0;
		if( isign < 0 ) sd= -sd;
		qlimit= x+j;
		for(qx=x; qx< qlimit; qx++)
		   {
			for(px=qx; px<limit; px+= istep)
			   {
				rx= px + j;
				real= cn * rx->re - sn * rx->im;
				imag= sn * rx->re + cn * rx->im;
				rx->re = px->re - real;
				rx->im = px->im - imag;
				px->re += real;
				px->im += imag;
			   }
			temp= cd * cn + sd * sn;
			sn += (sd * cn - cd * sn);
			cn -= temp;
		   }
		j= istep;
	   }
	return;
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

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}

void rand_init(float *rt,float *pct,long *seed,int ns,int nd,int nfs,int nfd,int nsmth,int gaus)
{
float *xt, maxr, gmean;
int i, j, k, l, ip;
int ix0, ix1, ix2, ix3, ix4;

xt = (float *) check_malloc (ns*nd*nfs*nfd*sizeof(float));

if(gaus)
   {
   gmean = 0.0;
   maxr = 0.0;
   for(i=0;i<ns*nd*nfs*nfd;i++)
      {
      rt[i] = gaus_rand(pct,&gmean,seed);

      if(rt[i] > maxr)
         maxr = rt[i];
      if(-rt[i] > maxr)
         maxr = -rt[i];
      }

   maxr = *pct/maxr;
   for(i=0;i<ns*nd*nfs*nfd;i++)
      rt[i] = maxr*rt[i];
   }
else
   {
for(i=0;i<ns*nd*nfs*nfd;i++)
   xt[i] = sfrand(seed);

while(nsmth--)
   {
   ix0 = 0;

   ix2 = ix0 + 1;
   ix4 = ix0 + nd*nfd;

   rt[ix0] = 0.2000*(xt[ix0] + xt[ix2] + xt[ix4]);

   for(j=1;j<nd*nfd-1;j++)
      {
      ix0 = j;

      ix1 = ix0 - 1;
      ix2 = ix0 + 1;
      ix4 = ix0 + nd*nfd;

      rt[ix0] = 0.20*(xt[ix0] + xt[ix1] + xt[ix2] + xt[ix4]);
      }

   ix0 = nd*nfd - 1;

   ix1 = ix0 - 1;
   ix4 = ix0 + nd*nfd;

   rt[ix0] = 0.2000*(xt[ix0] + xt[ix1] + xt[ix4]);

   for(i=1;i<ns*nfs-1;i++)
      {
      ix0 = i*nd*nfd;

      ix2 = ix0 + 1;
      ix3 = ix0 - nd*nfd;
      ix4 = ix0 + nd*nfd;

      rt[ix0] = 0.20*(xt[ix0] + xt[ix2] + xt[ix3] + xt[ix4]);

      for(j=1;j<nd*nfd-1;j++)
	 {
	 ix0 = j + i*nd*nfd;

	 ix1 = ix0 - 1;
	 ix2 = ix0 + 1;
	 ix3 = ix0 - nd*nfd;
	 ix4 = ix0 + nd*nfd;

         rt[ix0] = 0.2*(xt[ix0] + xt[ix1] + xt[ix2] + xt[ix3] + xt[ix4]);
	 }

      ix0 = nd*nfd - 1 + i*nd*nfd;

      ix1 = ix0 - 1;
      ix3 = ix0 - nd*nfd;
      ix4 = ix0 + nd*nfd;

      rt[ix0] = 0.20*(xt[ix0] + xt[ix1] + xt[ix3] + xt[ix4]);
      }

   ix0 = (ns*nfs-1)*nd*nfd;

   ix2 = ix0 + 1;
   ix3 = ix0 - nd*nfd;

   rt[ix0] = 0.2000*(xt[ix0] + xt[ix2] + xt[ix3]);

   for(j=1;j<nd*nfd-1;j++)
      {
      ix0 = j + (ns*nfs-1)*nd*nfd;

      ix1 = ix0 - 1;
      ix2 = ix0 + 1;
      ix3 = ix0 - nd*nfd;

      rt[ix0] = 0.20*(xt[ix0] + xt[ix1] + xt[ix2] + xt[ix3]);
      }

   ix0 = ns*nfs*nd*nfd - 1;

   ix1 = ix0 - 1;
   ix3 = ix0 - nd*nfd;

   rt[ix0] = 0.2000*(xt[ix0] + xt[ix1] + xt[ix3]);

   for(i=0;i<ns*nfs;i++)
      {
      for(j=0;j<nd*nfd;j++)
         {
	 ix0 = j + i*nd*nfd;
         xt[ix0] = rt[ix0];
	 }
      }
   }

maxr = 0.0;
for(i=0;i<ns;i++)
   {
   for(k=0;k<nfs;k++)
      {
      for(j=0;j<nd;j++)
         {
         for(l=0;l<nfd;l++)
	    {
	    ip =  l + k*nfd + (j + i*nd)*nfs*nfd;
	    ix0 = l + j*nfd + (k + i*nfs)*nd*nfd;

            rt[ip] = xt[ix0];

	    if(rt[ip] > maxr)
	       maxr = rt[ip];
	    if(-rt[ip] > maxr)
	       maxr = -rt[ip];
            }
         }
      }
   }
fprintf(stderr,"maxr=%13.5f\n",maxr);

maxr = (*pct)/maxr;
for(i=0;i<ns;i++)
   {
   for(j=0;j<nd;j++)
      {
      for(k=0;k<nfs;k++)
	 {
	 for(l=0;l<nfd;l++)
	    {
	    ip = l + k*nfd + (j + i*nd)*nfs*nfd;
            rt[ip] = rt[ip]*maxr;

fprintf(stderr,"rt=%13.5f\n",rt[ip]);
            }
         }
      }
   }
   }

free(xt);
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

double nt_tol(float fnt,int gnt)
{
double diff;

diff = ((double)(fnt) - (double)(gnt));
if(diff < 0.0)
   diff = -diff;

/*
fprintf(stderr,"diff= %15.10e\n",diff);
*/

return(diff);
}
