#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "fftw3.h"

extern void fourg_(struct complex *, int *, int *, float *);
void hermit(struct complex *s0,int nx0,int ny0);

void init_slip_IO(struct complex *sc,int nx2,int ny2,int nx,int ny,float *dx,float *dy,int flip,char *file)
{
FILE *fpr;
char str[1024];
float slip, x0, y0, x1, y1;
int ix0, iy0, ix1, iy1, iystart;
int ix, iy, ip, ip1;
int nr;

for(ip=0;ip<nx2*ny2;ip++)
   {
   sc[ip].re = 0.0;
   sc[ip].im = 0.0;
   }

if(file[0] == '\0')
   {
   for(ip=0;ip<nx2*ny2;ip++)
      {
      sc[ip].re = 0.5;
      sc[ip].im = 0.0;
      }
   }
else
   {
   fpr = fopfile(file,"r");

   while(fgets(str,1024,fpr) != NULL)
      {
      nr = sscanf(str,"%f %f %f %f %f",&slip,&x0,&y0,&x1,&y1);

      if(nr == 1)
         {
         for(ip=0;ip<nx2*ny2;ip++)
            {
            sc[ip].re = slip;
            sc[ip].im = 0.0;
            }
	 }
      else
         {
         iystart = 0;

         ix0 = (int)((0.5*nx*(*dx) + x0)/(*dx) + 0.5);
         if(ix0 < 0)
            ix0 = 0;
         if(ix0 > nx2)
            ix0 = nx2;

         iy0 = iystart + (int)(y0/(*dy) + 0.5);
         if(iy0 < iystart)
            iy0 = iystart;
         if(iy0 > ny2)
            iy0 = ny2;

         ix1 = (int)((0.5*nx*(*dx) + x1)/(*dx) + 0.5);
         if(ix1 < 0)
            ix1 = 0;
         if(ix1 > nx2)
            ix1 = nx2;

         iy1 = iystart + (int)(y1/(*dy) + 0.5) + 1;
         if(iy1 < iystart)
            iy1 = iystart;
         if(iy1 > ny2)
            iy1 = ny2;

         fprintf(stderr,"ix0= %d ix1= %d iy0= %d iy1= %d\n",ix0,ix1,iy0,iy1);

         for(iy=iy0;iy<iy1;iy++)
            {
            for(ix=ix0;ix<ix1;ix++)
               {
	       ip = ix+iy*nx2;
	       sc[ip].re = slip;
	       }
	    }
	 }
      }

   if(flip == 1)
      {
      for(iy=0;iy<ny2/2;iy++)
         {
         for(ix=0;ix<nx2;ix++)
            {
            ip = ix+iy*nx2;
            ip1 = ix+((ny2-1)-iy)*nx2;

            sc[ip1].re = sc[ip].re;
            }
         }
      }
   }
}

void init_slip(struct complex *sc,int nx,int ny,float *st,float *bt)
{
float xdamp, ydamp;
int i;

for(i=0;i<nx*ny;i++)
   {
   sc[i].re = 0.5;
   sc[i].im = 0.0;
   }

/*
taper_slip(sc,nx,ny,st,bt);
*/
}

void taper_slip_all(struct complex *sc,int nx,int ny,float *st,float *bt,float *tt)
{
float xdamp, ydamp;
int ix, iy, xb, yt, yb;

xb = (int)((*st)*nx + 0.5);
if(xb < 0)
   xb = 0;

yt = (int)((*tt)*ny + 0.5);
if(yt < 0)
   yt = 0;

yb = (int)((*bt)*ny + 0.5);
if(yb < 0)
   yb = 0;

for(iy=0;iy<yt;iy++)
   {
   ydamp = (float)(iy+1)/(float)(yt);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < xb)
         xdamp = (float)(ix+1)/(float)(xb);
      if(ix > nx-xb)
         xdamp = (float)(nx-ix)/(float)(xb);

      sc[ix+iy*nx].re = xdamp*ydamp*sc[ix+iy*nx].re;
      }
   }

for(iy=0;iy<yb;iy++)
   {
   ydamp = (float)(iy+1)/(float)(yb);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < xb)
         xdamp = (float)(ix+1)/(float)(xb);
      if(ix > nx-xb)
         xdamp = (float)(nx-ix)/(float)(xb);

      sc[ix+(ny-1-iy)*nx].re = xdamp*ydamp*sc[ix+(ny-1-iy)*nx].re;
      }
   }

for(iy=yt;iy<ny-yb;iy++)
   {
   for(ix=0;ix<xb;ix++)
      {
      xdamp = (float)(ix+1)/(float)(xb);

      sc[ix+iy*nx].re = xdamp*sc[ix+iy*nx].re;
      sc[(nx-1-ix)+iy*nx].re = xdamp*sc[(nx-1-ix)+iy*nx].re;
      }
   }
}

void taper_slip_all_r(float *sr,int nx,int ny,float *st,float *bt,float *tt)
{
float xdamp, ydamp;
int ix, iy, xb, yt, yb;

xb = (int)((*st)*nx + 0.5);
if(xb < 0)
   xb = 0;

yt = (int)((*tt)*ny + 0.5);
if(yt < 0)
   yt = 0;

yb = (int)((*bt)*ny + 0.5);
if(yb < 0)
   yb = 0;

for(iy=0;iy<yt;iy++)
   {
   ydamp = (float)(iy+1)/(float)(yt);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < xb)
         xdamp = (float)(ix+1)/(float)(xb);
      if(ix > nx-xb)
         xdamp = (float)(nx-ix)/(float)(xb);

      sr[ix+iy*nx] = xdamp*ydamp*sr[ix+iy*nx];
      }
   }

for(iy=0;iy<yb;iy++)
   {
   ydamp = (float)(iy+1)/(float)(yb);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < xb)
         xdamp = (float)(ix+1)/(float)(xb);
      if(ix > nx-xb)
         xdamp = (float)(nx-ix)/(float)(xb);

      sr[ix+(ny-1-iy)*nx] = xdamp*ydamp*sr[ix+(ny-1-iy)*nx];
      }
   }

for(iy=yt;iy<ny-yb;iy++)
   {
   for(ix=0;ix<xb;ix++)
      {
      xdamp = (float)(ix+1)/(float)(xb);

      sr[ix+iy*nx] = xdamp*sr[ix+iy*nx];
      sr[(nx-1-ix)+iy*nx] = xdamp*sr[(nx-1-ix)+iy*nx];
      }
   }
}

void taper_slip_all_r_hc(float *sr,int nx,int ny,float *st,float *bt,float *tt,struct hybrid_correlation_length *hc,struct pointsource *ps)
{
float xdamp, ydamp, dmin1, dmax1;
int ix, iy, ixp, ixdeep, ixshal, iyt, iyb;

dmin1 = hc->dep - hc->dep_range;
dmax1 = hc->dep + hc->dep_range;

ixshal = (int)((hc->side_taper)*nx + 0.5);
if(ixshal < 0)
   ixshal = 0;

ixdeep = (int)((*st)*nx + 0.5);
if(ixdeep < 0)
   ixdeep = 0;

iyt = (int)((*tt)*ny + 0.5);
if(iyt < 0)
   iyt = 0;

iyb = (int)((*bt)*ny + 0.5);
if(iyb < 0)
   iyb = 0;

for(iy=0;iy<ny;iy++)
   {
   ydamp = 1.0;
   if(iy < iyt)
      ydamp = (float)(iy+1)/(float)(iyt);
   if(iy > ny-iyb)
      ydamp = (float)(ny-iy)/(float)(iyb);

   if(ps[iy*nx].dep <= dmin1)
      ixp = ixshal;
   else if(ps[iy*nx].dep < dmax1 && ps[iy*nx].dep > dmin1)
      ixp = ixdeep - (int)((ixdeep-ixshal)*(dmax1-(ps[iy*nx].dep))/(dmax1-dmin1) + 0.5);
   else
      ixp = ixdeep;

   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < ixp)
         xdamp = (float)(ix+1)/(float)(ixp);
      if(ix > nx-ixp)
         xdamp = (float)(nx-ix)/(float)(ixp);

      sr[ix+iy*nx] = xdamp*ydamp*sr[ix+iy*nx];
      }
   }
}

void taper_slip(struct complex *sc,int nx,int ny,float *st,float *bt)
{
float xdamp, ydamp;
int ix, iy, xb, yb;

xb = (int)((*st)*nx + 0.5);
if(xb < 0)
   xb = 0;

yb = (int)((*bt)*ny + 0.5);
if(yb < 0)
   yb = 0;

for(iy=0;iy<yb;iy++)
   {
   ydamp = (float)(iy+1)/(float)(yb);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < xb)
         xdamp = (float)(ix+1)/(float)(xb);
      if(ix > nx-xb)
         xdamp = (float)(nx-ix)/(float)(xb);

      sc[ix+iy*nx].re = xdamp*ydamp*sc[ix+iy*nx].re;
      sc[ix+(ny-1-iy)*nx].re = xdamp*ydamp*sc[ix+(ny-1-iy)*nx].re;
      }
   }

for(iy=yb;iy<ny-yb;iy++)
   {
   for(ix=0;ix<xb;ix++)
      {
      xdamp = (float)(ix+1)/(float)(xb);

      sc[ix+iy*nx].re = xdamp*sc[ix+iy*nx].re;
      sc[(nx-1-ix)+iy*nx].re = xdamp*sc[(nx-1-ix)+iy*nx].re;
      }
   }
}

void scale_slip(struct pointsource *ps,struct complex *cs,int nx,int ny,int nys,float *dx,float *dy,float *dtop,float *dip,float *mom,struct velmodel *vm,float *savg,float *smax)
{
float sum, fac, sinD, area;
float zz;
int i, j, k;

float rperd = 0.017453293;

if(*savg < 0.0)
   {
   sinD = sin((*dip)*rperd);
   area = (*dx)*(*dy)*1.0e+10;    /* in CMS units */

   sum = 0.0;
   for(j=0;j<ny;j++)
      {
      zz = (*dtop) + sinD*(j + 0.5)*(*dy);

      k = 0;
      while(zz > vm->dep[k] && k < (vm->nlay)-1)
         k++;

      fac = area*vm->mu[k];
      for(i=0;i<nx;i++)
         sum = sum + fac*cs[i + (j+nys)*nx].re;
      }

   fac = (*mom)/sum;

   *smax = 0.0;
   *savg = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         ps[i + j*nx].slip = fac*cs[i + (j+nys)*nx].re;

         *savg = *savg + ps[i + j*nx].slip;
         if(ps[i + j*nx].slip > *smax)
            *smax = ps[i + j*nx].slip;
         }
      }
   *savg = (*savg)/(float)(nx*ny);
   }
else
   {
   sinD = sin((*dip)*rperd);
   area = (*dx)*(*dy)*1.0e+10;    /* in CMS units */

   sum = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         sum = sum + cs[i + (j+nys)*nx].re;
         }
      }

   fac = (*savg)*(float)(nx*ny)/(sum);

   *smax = 0.0;
   *savg = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         ps[i + j*nx].slip = fac*cs[i + (j+nys)*nx].re;

         *savg = *savg + ps[i + j*nx].slip;
         if(ps[i + j*nx].slip > *smax)
            *smax = ps[i + j*nx].slip;
         }
      }
   *savg = (*savg)/(float)(nx*ny);
   }
}

void scale_slip_r(struct pointsource *ps,float *sr,int nx,int ny,int nys,float *dx,float *dy,float *dtop,float *dip,float *mom,struct velmodel *vm,float *savg,float *smax)
{
float sum, fac, sinD, area;
float zz;
int i, j, k;

float rperd = 0.017453293;

sinD = sin((*dip)*rperd);
area = (*dx)*(*dy)*1.0e+10;    /* in CMS units */

if(*savg < 0.0)
   {
   sum = 0.0;
   for(j=0;j<ny;j++)
      {
      zz = (*dtop) + sinD*(j + 0.5)*(*dy);

      k = 0;
      while(zz > vm->dep[k] && k < (vm->nlay)-1)
         k++;

      fac = area*vm->mu[k];
      for(i=0;i<nx;i++)
         sum = sum + fac*sr[i + (j+nys)*nx];
      }

   fac = (*mom)/sum;

   *smax = 0.0;
   *savg = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         ps[i + j*nx].slip = fac*sr[i + (j+nys)*nx];

         *savg = *savg + ps[i + j*nx].slip;
         if(ps[i + j*nx].slip > *smax)
            *smax = ps[i + j*nx].slip;
         }
      }
   *savg = (*savg)/(float)(nx*ny);
   }
else
   {
   sinD = sin((*dip)*rperd);
   area = (*dx)*(*dy)*1.0e+10;    /* in CMS units */

   sum = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         sum = sum + sr[i + (j+nys)*nx];
         }
      }

   fac = (*savg)*(float)(nx*ny)/(sum);

   *smax = 0.0;
   *savg = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         ps[i + j*nx].slip = fac*sr[i + (j+nys)*nx];

         *savg = *savg + ps[i + j*nx].slip;
         if(ps[i + j*nx].slip > *smax)
            *smax = ps[i + j*nx].slip;
         }
      }
   *savg = (*savg)/(float)(nx*ny);

   sum = 0.0;
   for(j=0;j<ny;j++)
      {
      zz = (*dtop) + sinD*(j + 0.5)*(*dy);

      k = 0;
      while(zz > vm->dep[k] && k < (vm->nlay)-1)
         k++;

      fac = area*vm->mu[k];
      for(i=0;i<nx;i++)
         sum = sum + fac*sr[i + (j+nys)*nx];
      }
   *mom = sum;
   }
}

void scale_slip_r_vsden(struct pointsource *ps,float *sr,int nstk,int ndip,int nys,float *dx,float *dy,float *dtop,float *dip,float *mom,float *savg,float *smax)
{
float sum, fac, sinD, area;
float zz;
int i, j, k;

float rperd = 0.017453293;

sinD = sin((*dip)*rperd);
area = (*dx)*(*dy)*1.0e+10;    /* in CMS units */

if(*savg < 0.0)
   {
   sum = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         sum = sum + area*ps[i + j*nstk].mu*sr[i + (j+nys)*nstk];
      }

   fac = (*mom)/sum;

   *smax = 0.0;
   *savg = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
         ps[i + j*nstk].slip = fac*sr[i + (j+nys)*nstk];

         *savg = *savg + ps[i + j*nstk].slip;
         if(ps[i + j*nstk].slip > *smax)
            *smax = ps[i + j*nstk].slip;
         }
      }
   *savg = (*savg)/(float)(nstk*ndip);
   }
else
   {
   sinD = sin((*dip)*rperd);
   area = (*dx)*(*dy)*1.0e+10;    /* in CMS units */

   sum = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
         sum = sum + sr[i + (j+nys)*nstk];
         }
      }

   fac = (*savg)*(float)(nstk*ndip)/(sum);

   *smax = 0.0;
   *savg = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
         ps[i + j*nstk].slip = fac*sr[i + (j+nys)*nstk];

         *savg = *savg + ps[i + j*nstk].slip;
         if(ps[i + j*nstk].slip > *smax)
            *smax = ps[i + j*nstk].slip;
         }
      }
   *savg = (*savg)/(float)(nstk*ndip);

   sum = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         sum = sum + area*ps[i + j*nstk].mu*sr[i + (j+nys)*nstk];
      }
   *mom = sum;
   }
}

void kfilt(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float phsb;
int ndkc;

float pi = 3.14159265;
float hcoef;

hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

/*

   Transition between deterministic and stochastic parts of spectrum
   are given by

       F = wtS*stoch + wtD*deter

   with
    
       wtD = {1 + k2/Kc2}^-(xp)     (kind of a butterworth filter)
       wtS = 1 - wtD

   and

       k2 = kx*kx + ky*ky          (k-squared)
       Kc2 = (N*dky)*(N*dkx)       (corner wavenumber of transition)

   The parameter N specifies the number of dk's in the corner (somewhat
   like a fraction of the total wavenumber space).  The exponent (xp)
   gives the sharpness of the transition.  Based on very limited
   testing, I came up with

       N = 4
       xp = 2.0

*/

xp = 2.0;
ndkc = 4;

invkc2 = ndkc*(*dky)*ndkc*(*dkx);
invkc2 = 1.0/(invkc2);

xp = 1;
ndkc = 2;

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      phs = pi*sfrand(seed);

      fac1 = sqrt(s0[ip].re*s0[ip].re + s0[ip].im*s0[ip].im);

      phs1 = 0.5*pi;
      if(s0[ip].re != (float)(0.0))
         {
         phs1 = atan(s0[ip].im/s0[ip].re);
         if(s0[ip].re < 0.0)
            phs1 = phs1 + pi;
         }
      else if(s0[ip].im < 0.0)
         phs1 = -0.5*pi;

/* 09/24/2009 I have no idea why I put these conditions here, now they're removed
      while(phs1 > pi)
         phs1 = phs1 - pi;
      while(phs1 < -pi)
         phs1 = phs1 + pi;
*/

      k2 = (kx*kx + ky*ky)*invkc2;
      wtD = exp(-xp*log(1.0 + k2));

      k2 = exp(xp*log(kx*kx/(ndkc*(*dkx)*ndkc*(*dkx)))) + exp(xp*log(ky*ky/(ndkc*(*dky)*ndkc*(*dky))));
      wtD = 1.0/(1.0 + k2);
      k2 = kx*kx/(ndkc*(*dkx)*ndkc*(*dkx)) + ky*ky/(ndkc*(*dky)*ndkc*(*dky));
      wtD = 1.0/(1.0 + exp(xp*log(k2)));

      wtS = 1.0 - wtD;

      s0[ip].re = wtS*fac*cos(phs) + wtD*fac1*cos(phs1);
      s0[ip].im = wtS*fac*sin(phs) + wtD*fac1*sin(phs1);

/*
wtD = 0.0;
wtS = 20.0;
s0[ip].re = fac*gaus_rand(&wtS,&wtD,seed)/sqrt(2.0*wtS*wtS);
s0[ip].im = fac*gaus_rand(&wtS,&wtD,seed)/sqrt(2.0*wtS*wtS);
*/

/*

   OLD STUFF that I tried FOLLOWS from HERE

*/

/*
   Do not alter phase of lowest (k=0) and 2nd lowest (k=dk) wavenumbers
   so average slip and edge taper are not significantly modified
      if(i > 1 && j > 1)
*/

/* 
   Do not alter phase of lowest (k=0) wavenumbers and use average of
   random and determinstic phase for 2nd lowest (k=dk)
   so average slip and edge taper are not significantly modified
*/

/*
      if(i > 1 && j > 1)
         {
         s0[ip].re = fac*cos(phs);
         s0[ip].im = fac*sin(phs);
	 }
      else if((i == 1 && j > 0) || (i > 0 && j == 1))
         {
	 if(i == 1)
	    wtS = sqrt(ky*ky*invkym2);
	 if(j == 1)
	    wtS = sqrt(kx*kx*invkxm2);

	 if(wtS < 0.5)
	    wtS = 0.5;

	 wtD = 1.0 - wtS;

	 fac1 = sqrt(s0[ip].re*s0[ip].re + s0[ip].im*s0[ip].im);

	 phs1 = 0.5*pi;
	 if(s0[ip].re != 0.0)
	    {
	    phs1 = atan(s0[ip].im/s0[ip].re);
	    if(s0[ip].re < 0.0)
	       phs1 = phs1 + pi;
	    }
	 else if(s0[ip].im < 0.0)
	    phs1 = -0.5*pi;

         while(phs1 > pi)
            phs1 = phs1 - pi;
         while(phs1 < -pi)
            phs1 = phs1 + pi;

         s0[ip].re = wtS*fac*cos(phs) + wtD*fac1*cos(phs1);
         s0[ip].im = wtS*fac*sin(phs) + wtD*fac1*sin(phs1);
	 }
*/
      }
   }

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit(s0,nx0,ny0);

/*
for(j=0;j<ny0;j++)
   {
   for(i=0;i<nx0;i++)
      {
      ip = i + j*nx0;

      if(j <= ny0/2)
         ky = j*(*dky);
      else
         ky = (j - ny0)*(*dky);
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      phs = -2.0*pi*kx*10.0;

      fre = s0[ip].re*cos(phs) - s0[ip].im*sin(phs);
      fim = s0[ip].re*sin(phs) + s0[ip].im*cos(phs);

      phs = -2.0*pi*ky*10.0;

      s0[ip].re = fre*cos(phs) - fim*sin(phs);
      s0[ip].im = fre*sin(phs) + fim*cos(phs);
      }
   }
*/
}

void fft2d(struct complex *xc,int n1,int n2,int isgn,float *d1,float *d2)
{
int i, j, ip;
float *space;
struct complex *xtc;
float normf;

normf = (*d1)*(*d2);

space = (float *) check_malloc (2*(n1+n2)*sizeof(float));

for(j=0;j<n2;j++)
   fourg_(xc+j*n1,&n1,&isgn,space);

xtc = (struct complex *) check_malloc (n2*sizeof(struct complex));

for(i=0;i<n1;i++)
   {
   for(j=0;j<n2;j++)
      {
      ip = i + j*n1;

      xtc[j].re = xc[ip].re;
      xtc[j].im = xc[ip].im;
      }

   fourg_(xtc,&n2,&isgn,space);

   for(j=0;j<n2;j++)
      {
      ip = i + j*n1;

      xc[ip].re = normf*xtc[j].re;
      xc[ip].im = normf*xtc[j].im;
      }
   }

free(xtc);
free(space);
}

/* 20160914: fftw from Scott C. */

void fft2d_fftw(struct complex *xc,int n1,int n2,int isgn,float *d1,float *d2)
{
off_t i, j, ip;
float *space;
struct complex *xtc;
float normf;

normf = (*d1)*(*d2);

fftwf_plan plan;
fftwf_complex *arr = check_malloc(sizeof(fftwf_complex)*n1);

if (isgn==-1) {
        plan = fftwf_plan_dft_1d(n1, arr, arr, FFTW_FORWARD, FFTW_ESTIMATE);
} else if (isgn==1){
        plan = fftwf_plan_dft_1d(n1, arr, arr, FFTW_BACKWARD, FFTW_ESTIMATE);
}

for(j=0;j<n2;j++) {
   for (i=0; i<n1; i++) {
        arr[i][0] = (xc+j*n1+i)->re;
        arr[i][1] = (xc+j*n1+i)->im;
   }
   fftwf_execute(plan);
   for (i=0; i<n1; i++) {
        (xc+j*n1+i)->re = arr[i][0];
        (xc+j*n1+i)->im = arr[i][1];
   }
}
arr = check_realloc(arr, n2*sizeof(fftwf_complex));

if (isgn==-1) {
        plan = fftwf_plan_dft_1d(n2, arr, arr, FFTW_FORWARD, FFTW_ESTIMATE);
} else if (isgn==1){
        plan = fftwf_plan_dft_1d(n2, arr, arr, FFTW_BACKWARD, FFTW_ESTIMATE);
}

for(i=0;i<n1;i++)
   {
   for(j=0;j<n2;j++)
      {
      ip = i + j*n1;

      arr[j][0] = xc[ip].re;
      arr[j][1] = xc[ip].im;
      }

   fftwf_execute(plan);

   for(j=0;j<n2;j++)
      {
      ip = i + j*n1;

      xc[ip].re = normf*arr[j][0];
      xc[ip].im = normf*arr[j][1];
      }
   }

fftwf_destroy_plan(plan);
fftwf_free(arr);
}

void kfilt_lw(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag,float *fl,float *fw)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float phsb;
int ndkc;

float pi = 3.14159265;
float hcoef;

hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

/*

   Transition between deterministic and stochastic parts of spectrum
   are given by

       F = wtS*stoch + wtD*deter

   with
    
       wtD = {1 + k2/Kc2}^-(xp)     (kind of a butterworth filter)
       wtS = 1 - wtD

   and

       k2 = kx*kx + ky*ky          (k-squared)
       Kc2 = (N*dky)*(N*dkx)       (corner wavenumber of transition)

   The parameter N specifies the number of dk's in the corner (somewhat
   like a fraction of the total wavenumber space).  The exponent (xp)
   gives the sharpness of the transition.  Based on very limited
   testing, I came up with

       N = 4
       xp = 2.0

*/

xp = 2.0;
ndkc = 4;

invkc2 = ndkc*(*dky)*ndkc*(*dkx);
invkc2 = 1.0/(invkc2);

/* RWG 20130124

updated above to be consistent with GP2010 paper including explicitly putting Flen and Fwid
into filter corners which corrects problem with short faults when flip_at_surface=1

*/

xp = 1;
ndkc = 2;

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      phs = pi*sfrand(seed);

      fac1 = sqrt(s0[ip].re*s0[ip].re + s0[ip].im*s0[ip].im);

      phs1 = 0.5*pi;
      if(s0[ip].re != (float)(0.0))
         {
         phs1 = atan(s0[ip].im/s0[ip].re);
         if(s0[ip].re < 0.0)
            phs1 = phs1 + pi;
         }
      else if(s0[ip].im < 0.0)
         phs1 = -0.5*pi;

      k2 = (*fl)*(*fl)*kx*kx/(ndkc*ndkc) + (*fw)*(*fw)*ky*ky/(ndkc*ndkc);
      wtD = 1.0/(1.0 + exp(xp*log(k2)));

      wtS = 1.0 - wtD;

      s0[ip].re = wtS*fac*cos(phs) + wtD*fac1*cos(phs1);
      s0[ip].im = wtS*fac*sin(phs) + wtD*fac1*sin(phs1);
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit(s0,nx0,ny0);
}

void kfilt_lw_H0(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag,float *fl,float *fw)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float phsb;
int ndkc;

float pi = 3.14159265;
float hcoef;

hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */
hcoef = 1.0; /* H=0.0, hcoef = H + 1 */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

/*

   Transition between deterministic and stochastic parts of spectrum
   are given by

       F = wtS*stoch + wtD*deter

   with
    
       wtD = {1 + k2/Kc2}^-(xp)     (kind of a butterworth filter)
       wtS = 1 - wtD

   and

       k2 = kx*kx + ky*ky          (k-squared)
       Kc2 = (N*dky)*(N*dkx)       (corner wavenumber of transition)

   The parameter N specifies the number of dk's in the corner (somewhat
   like a fraction of the total wavenumber space).  The exponent (xp)
   gives the sharpness of the transition.  Based on very limited
   testing, I came up with

       N = 4
       xp = 2.0

*/

xp = 2.0;
ndkc = 4;

invkc2 = ndkc*(*dky)*ndkc*(*dkx);
invkc2 = 1.0/(invkc2);

/* RWG 20130124

updated above to be consistent with GP2010 paper including explicitly putting Flen and Fwid
into filter corners which corrects problem with short faults when flip_at_surface=1

*/

xp = 1;
ndkc = 2;

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      phs = pi*sfrand(seed);

      fac1 = sqrt(s0[ip].re*s0[ip].re + s0[ip].im*s0[ip].im);

      phs1 = 0.5*pi;
      if(s0[ip].re != (float)(0.0))
         {
         phs1 = atan(s0[ip].im/s0[ip].re);
         if(s0[ip].re < 0.0)
            phs1 = phs1 + pi;
         }
      else if(s0[ip].im < 0.0)
         phs1 = -0.5*pi;

      k2 = (*fl)*(*fl)*kx*kx/(ndkc*ndkc) + (*fw)*(*fw)*ky*ky/(ndkc*ndkc);
      wtD = 1.0/(1.0 + exp(xp*log(k2)));

      wtS = 1.0 - wtD;

      s0[ip].re = wtS*fac*cos(phs) + wtD*fac1*cos(phs1);
      s0[ip].im = wtS*fac*sin(phs) + wtD*fac1*sin(phs1);

      }
   }

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit(s0,nx0,ny0);
}

/* 2015-03-13 RWG
   Following Liu et al (2006), use only stochastic part in kfilt_gaus()
   Appears to work fine as long as factor of 2*pi is backed-out of Somerville wavenumber corners
   and Mai corners are unchanged from published vaules, see comments in main().
*/

void kfilt_gaus(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
int ndkc;

float pi = 3.14159265;
float hcoef, beta2;

float fnorm;
float fzero = 0.0;
float fone = 1.0;

fnorm = 1.0/sqrt(2.0);

hcoef = 2.00; /* H=1.00, hcoef = H + 1 */
hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

beta2 = hcoef;
beta2 = 2.0; /* strictly self-similar */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == FRANKEL_FLAG) /* frankel(2009) scaling */
         {
	 if(amp < 1.0)
	    fac = amp0;
         else
            fac = amp0*exp(-0.5*beta2*log(amp));
	 }

      fre = fnorm*gaus_rand(&fone,&fzero,seed);
      fim = fnorm*gaus_rand(&fone,&fzero,seed);

      s0[ip].re = fac*fre;
      s0[ip].im = fac*fim;
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit(s0,nx0,ny0);
}

void kfilt_rphs(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float phsb;
int ndkc;

float pi = 3.14159265;
float hcoef;

hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      phs = pi*sfrand(seed);

      wtS = 1.0 - wtD;

      s0[ip].re = fac*cos(phs);
      s0[ip].im = fac*sin(phs);
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit(s0,nx0,ny0);
}

/* 2015-03-26 RWG
   Following Shi and Day (2014), kfilt_beta2() produces field with zero mean
   and k^(-beta/2) falloff. Beta/2 = hcoef + 1, hcoef=1 for self-similar.
*/

void kfilt_beta2(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *hcoef,float *lambda_max,float *lambda_min,long *seed)
{
int i, j, ip;
float kx, ky, amp, fre, fim, beta2;
float fac, lmax2, lmin2, k2;

float amp0;

float fnorm;
float fzero = 0.0;
float fone = 1.0;

int ord = 4;

fnorm = 1.0/sqrt(2.0);
beta2 = *hcoef + 1.0;

lmax2 = (*lambda_max)*(*lambda_max);
lmin2 = (*lambda_min)*(*lambda_min);

amp0 = exp(0.5*beta2*log(lmax2));	/* amplitude when k2=1/lmax2 */

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      fre = fnorm*gaus_rand(&fone,&fzero,seed);
      fim = fnorm*gaus_rand(&fone,&fzero,seed);

      if(i == 0 && j == 0)
         amp = 0.0;
      else
         {
         k2 = kx*kx + ky*ky;
         fac = 1.0/((1.0 + exp(ord*log(k2*lmin2)))*(1.0 + exp(-ord*log(k2*lmax2))));
         amp = fac*exp(-0.5*beta2*log(k2));
	 }

      ip = i + j*nx0;

      s0[ip].re = amp*fre;
      s0[ip].im = amp*fim;
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].re = 0.0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit(s0,nx0,ny0);
}

void kfilt_beta2_ll(struct complex *s0,off_t nx0,off_t ny0,float *dkx,float *dky,float *hcoef,float *lambda_max,float *lambda_min,long *seed)
{
off_t i, j, ip;
float kx, ky, amp, fre, fim, beta2;
float fac, lmax2, lmin2, k2;

float amp0;

float fnorm;
float fzero = 0.0;
float fone = 1.0;

int ord = 4;

fnorm = 1.0/sqrt(2.0);
beta2 = *hcoef + 1.0;

lmax2 = (*lambda_max)*(*lambda_max);
lmin2 = (*lambda_min)*(*lambda_min);

amp0 = exp(0.5*beta2*log(lmax2));	/* amplitude when k2=1/lmax2 */

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   ky = j*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      fre = fnorm*gaus_rand(&fone,&fzero,seed);
      fim = fnorm*gaus_rand(&fone,&fzero,seed);

      if(i == 0 && j == 0)
         amp = 0.0;
      else
         {
         k2 = kx*kx + ky*ky;
         fac = 1.0/((1.0 + exp(ord*log(k2*lmin2)))*(1.0 + exp(-ord*log(k2*lmax2))));
         amp = fac*exp(-0.5*beta2*log(k2));
	 }

      ip = i + j*nx0;

      s0[ip].re = amp*fre;
      s0[ip].im = amp*fim;
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].re = 0.0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

hermit_ll(s0,nx0,ny0);
}

/* 2016-10-21 RWG
   added option to high- and low-pass filter random fields
*/

void kfilt_gaus2(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag,float *lambda_max,float *lambda_min)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float lmax2, lmin2;
int ndkc;

int ord = 4;

float hcoef, beta2;

float fnorm;
float fzero = 0.0;
float fone = 1.0;

fnorm = 1.0/sqrt(2.0);

hcoef = 2.00; /* H=1.00, hcoef = H + 1 */
hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

beta2 = hcoef;
beta2 = 2.0; /* strictly self-similar */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

lmax2 = (*lambda_max)*(*lambda_max);
lmin2 = (*lambda_min)*(*lambda_min);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      /* von Karman */
      if(kflag == MAI_FLAG || kflag == SUZUKI_FLAG || kflag == MAI_SOMERVILLE_HYBRID_FLAG)
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == FRANKEL_FLAG) /* frankel(2009) scaling */
         {
	 if(amp < 1.0)
	    fac = amp0;
         else
            fac = amp0*exp(-0.5*beta2*log(amp));
	 }

      k2 = kx*kx + ky*ky;
      if(k2 > 0.0)
         fac = fac/((1.0 + exp(ord*log(k2*lmin2)))*(1.0 + exp(-ord*log(k2*lmax2))));

      fre = fnorm*gaus_rand(&fone,&fzero,seed);
      fim = fnorm*gaus_rand(&fone,&fzero,seed);

      s0[ip].re = fac*fre;
      s0[ip].im = fac*fim;
      }
   }

/* 
   Enforce DC & Nyquists to be real and apply normalization (liu et al., 2006)
*/

/* following lines are what was used for v5.5.2 and earlier
s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;
*/

/* more correct form: following lines are what was used for v5.6.0 and later */
s0[0].re = s0[0].re/fnorm;
s0[0].im = 0.0;

s0[nx0/2].re = s0[nx0/2].re/fnorm;
s0[nx0/2].im = 0.0;

s0[nx0*ny0/2].re = s0[nx0*ny0/2].re/fnorm;
s0[nx0*ny0/2].im = 0.0;

s0[nx0/2+nx0*ny0/2].re = s0[nx0/2+nx0*ny0/2].re/fnorm;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

hermit(s0,nx0,ny0);
}

void kfilt_gaus2_v552(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag,float *lambda_max,float *lambda_min)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float lmax2, lmin2;
int ndkc;

int ord = 4;

float hcoef, beta2;

float fnorm;
float fzero = 0.0;
float fone = 1.0;

fnorm = 1.0/sqrt(2.0);

hcoef = 2.00; /* H=1.00, hcoef = H + 1 */
hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

beta2 = hcoef;
beta2 = 2.0; /* strictly self-similar */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

lmax2 = (*lambda_max)*(*lambda_max);
lmin2 = (*lambda_min)*(*lambda_min);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);
      
      ip = i + j*nx0;
      
      amp = kx*kx*xl2 + ky*ky*yl2;
      
      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);
      
      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
         }
      
      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);
      
      if(kflag == FRANKEL_FLAG) /* frankel(2009) scaling */
         {
         if(amp < 1.0)
            fac = amp0;
         else
            fac = amp0*exp(-0.5*beta2*log(amp));
         }
      
      k2 = kx*kx + ky*ky;
      if(k2 > 0.0)
         fac = fac/((1.0 + exp(ord*log(k2*lmin2)))*(1.0 + exp(-ord*log(k2*lmax2))));
      
      fre = fnorm*gaus_rand(&fone,&fzero,seed);
      fim = fnorm*gaus_rand(&fone,&fzero,seed);
      
      s0[ip].re = fac*fre;
      s0[ip].im = fac*fim;
      }
   }
   
/*    
 *       Enforce DC & Nyquists to be real
 *       */

s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

hermit(s0,nx0,ny0);
}

void kfilt_gaus2_v552x(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag,float *lambda_max,float *lambda_min)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float lmax2, lmin2;
int ndkc;

int ord = 4;

float hcoef, beta2;

float fnorm;
float fzero = 0.0;
float fone = 1.0;

fnorm = 1.0/sqrt(2.0);

hcoef = 2.00; /* H=1.00, hcoef = H + 1 */
hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

beta2 = hcoef;
beta2 = 2.0; /* strictly self-similar */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

lmax2 = (*lambda_max)*(*lambda_max);
lmin2 = (*lambda_min)*(*lambda_min);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      /* von Karman */
      if(kflag == MAI_FLAG || kflag == SUZUKI_FLAG || kflag == MAI_SOMERVILLE_HYBRID_FLAG)
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == FRANKEL_FLAG) /* frankel(2009) scaling */
         {
	 if(amp < 1.0)
	    fac = amp0;
         else
            fac = amp0*exp(-0.5*beta2*log(amp));
	 }

      k2 = kx*kx + ky*ky;
      if(k2 > 0.0)
         fac = fac/((1.0 + exp(ord*log(k2*lmin2)))*(1.0 + exp(-ord*log(k2*lmax2))));

      fre = fnorm*gaus_rand(&fone,&fzero,seed);
      fim = fnorm*gaus_rand(&fone,&fzero,seed);

      s0[ip].re = fac*fre;
      s0[ip].im = fac*fim;
      }
   }

/* 
   Enforce DC & Nyquists to be real and apply normalization (liu et al., 2006)
*/

/* following lines are what was used for v5.5.2 and earlier
*/
s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

hermit(s0,nx0,ny0);
}

/* 2017-xx-xx RWG
   added option to allow detreminstic specification of low wavenumbers/asperities
*/

void kfilter(struct complex *s0,off_t nx0,off_t ny0,float *dkx,float *dky,int ord,float *lambda_max,float *lambda_min)
{
off_t i, j, ip;
float kx, ky, fac, k2;
float lmax2, lmin2;

lmax2 = (*lambda_max)*(*lambda_max);
lmin2 = (*lambda_min)*(*lambda_min);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      k2 = kx*kx + ky*ky;
      fac = 1.0/((1.0 + exp(ord*log(k2*lmin2)))*(1.0 + exp(-ord*log(k2*lmax2))));

      s0[ip].re = fac*s0[ip].re;
      s0[ip].im = fac*s0[ip].im;
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit_ll(s0,nx0,ny0);
}

/* 2017-xx-xx RWG
   added option to allow detreminstic specification of low wavenumbers/asperities
*/

void kfilt_gaus2_asp(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag,float *lambda_max,float *lambda_min)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float lmax2, lmin2;
int ndkc;

int ord = 4;

float pi = 3.14159265;
float hcoef, beta2;

float fnorm;
float fzero = 0.0;
float fone = 1.0;

fnorm = 0.1/sqrt(2.0);
fnorm = 0.3/sqrt(2.0);
fnorm = 1.0/sqrt(2.0);

hcoef = 2.00; /* H=1.00, hcoef = H + 1 */
hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

beta2 = hcoef;
beta2 = 2.0; /* strictly self-similar */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

lmax2 = (*lambda_max)*(*lambda_max);
lmin2 = (*lambda_min)*(*lambda_min);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == FRANKEL_FLAG) /* frankel(2009) scaling */
         {
	 if(amp < 1.0)
	    fac = amp0;
         else
            fac = amp0*exp(-0.5*beta2*log(amp));
	 }

/*   XXXXXX    for asperity preservation
*/
if(amp > 1.0e+20)
	 {
         k2 = kx*kx + ky*ky;
         if(k2 > 0.0)
            fac = fac/((1.0 + exp(ord*log(k2*lmin2)))*(1.0 + exp(-ord*log(k2*lmax2))));

         fre = fnorm*gaus_rand(&fone,&fzero,seed);
         fim = fnorm*gaus_rand(&fone,&fzero,seed);

         s0[ip].re = fac*fre;
         s0[ip].im = fac*fim;
         }
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit(s0,nx0,ny0);
}

/* 2018-06-22 RWG
   added option to apply phase-shift to random field
*/

void shift_phase(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,double *xshift,double *yshift)
{
int i, j, ip;
float kx, ky, fac, amp0;

float xre, yre, tre, xim, yim, tim;
double xarg, yarg, pi;

pi = 4.0*atan(1.0);
xarg = 2.0*pi*(*xshift);
yarg = 2.0*pi*(*yshift);

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      xre = cos(xarg*kx);
      xim = -sin(xarg*kx);
      yre = cos(yarg*ky);
      yim = -sin(yarg*ky);

      tre = xre*s0[ip].re - xim*s0[ip].im;
      tim = xre*s0[ip].im + xim*s0[ip].re;

      s0[ip].re = tre*yre - tim*yim;
      s0[ip].im = tre*yim + tim*yre;
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

/* 
for(j=1;j<=(ny0/2)-1;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0/2)-1;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
*/

hermit(s0,nx0,ny0);
}

void hermit(struct complex *s0,int nx0,int ny0)
{
int i, j;

for(i=1;i<(nx0/2);i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;

   s0[(nx0-i)+nx0*(ny0/2)].re = s0[i+nx0*(ny0/2)].re;
   s0[(nx0-i)+nx0*(ny0/2)].im = -s0[i+nx0*(ny0/2)].im;
   }

for(j=1;j<(ny0/2);j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;

   s0[(nx0/2)+(ny0-j)*nx0].re = s0[(nx0/2)+j*nx0].re;
   s0[(nx0/2)+(ny0-j)*nx0].im = -s0[(nx0/2)+j*nx0].im;
   }

for(j=1;j<ny0/2;j++)
   {
   for(i=1;i<nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
}

void hermit_ll(struct complex *s0,off_t nx0,off_t ny0)
{
off_t i, j;

for(i=1;i<(nx0/2);i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;

   s0[(nx0-i)+nx0*(ny0/2)].re = s0[i+nx0*(ny0/2)].re;
   s0[(nx0-i)+nx0*(ny0/2)].im = -s0[i+nx0*(ny0/2)].im;
   }

for(j=1;j<(ny0/2);j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;

   s0[(nx0/2)+(ny0-j)*nx0].re = s0[(nx0/2)+j*nx0].re;
   s0[(nx0/2)+(ny0-j)*nx0].im = -s0[(nx0/2)+j*nx0].im;
   }

for(j=1;j<ny0/2;j++)
   {
   for(i=1;i<nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
}

void init_asperity_IO(struct pointsource *ps,int nx2,int ny2,int nx,int ny,float *dx,float *dy,char *file,float *savg_old,float *tap_fac)
{
FILE *fpr;
char str[1024], *str_p, word[256], val[256];
float *slip_fac, x0, y0, x1, y1;
float *trise_fac, *rvf_fac;
int *ix0, *iy0, *ix1, *iy1;
int ix, iy, ip, ip1;
int j, k, nr;

float savg_new, area_asp, savg_asp, sum_num, sum_den;
int goflag, ia, nasp, nalloc, ntot;
int nblock = 100;

int ixt, iyt;
float xdamp, ydamp;

nalloc = nblock;
slip_fac = (float *)check_malloc(nalloc*sizeof(float));
trise_fac = (float *)check_malloc(nalloc*sizeof(float));
rvf_fac = (float *)check_malloc(nalloc*sizeof(float));
ix0 = (int *)check_malloc(nalloc*sizeof(int));
ix1 = (int *)check_malloc(nalloc*sizeof(int));
iy0 = (int *)check_malloc(nalloc*sizeof(int));
iy1 = (int *)check_malloc(nalloc*sizeof(int));

fprintf(stderr,"asperity params:\n\n");

fpr = fopfile(file,"r");

sum_num = 1.0;
sum_den = 1.0;
savg_asp = 0.0;
area_asp = 0.0;
ntot = 0;
nasp = 0;
while(fgets(str,1024,fpr) != NULL)
   {
   nr = sscanf(str,"%f %f %f %f %f",&slip_fac[nasp],&x0,&y0,&x1,&y1);
   ntot++;

/* initialize */
   trise_fac[nasp] = 1.0;
   rvf_fac[nasp] = 1.0;

   k = 0;
   str_p = str;
   while(str_p[k] != '\n' && str_p[k] != '\0')
      k++;

   k--;
   while(k > 0)
      {
      j = 0;
      while(strncmp(&str_p[k],"=",1) != 0 && k > 0)
         {
         j++;
         k--;
         }

      if(k == 0)
         break;

      strncpy(val,&str_p[k+1],j);
      val[j] = '\0';

      j = 0;
      k--;
      while(isblank(str_p[k]) == 0 && k > 0)
         {
         j++;
         k--;
         }

      if(k == 0)
         break;

      strncpy(word,&str_p[k+1],j);
      word[j] = '\0';

      if(strncmp(word,"trise_fac",9) == 0)
         trise_fac[nasp] = atof(val);
      if(strncmp(word,"rvf_fac",7) == 0)
         rvf_fac[nasp] = atof(val);
      }

   ix0[nasp] = (int)((0.5*nx*(*dx) + x0)/(*dx) + 0.5);
   if(ix0[nasp] < 0)
      ix0[nasp] = 0;
   if(ix0[nasp] > nx2)
      ix0[nasp] = nx2;

   iy0[nasp] = (int)(y0/(*dy) + 0.5);
   if(iy0[nasp] < 0)
      iy0[nasp] = 0;
   if(iy0[nasp] > ny2)
      iy0[nasp] = ny2;

   ix1[nasp] = (int)((0.5*nx*(*dx) + x1)/(*dx) + 0.5);
   if(ix1[nasp] < 0)
      ix1[nasp] = 0;
   if(ix1[nasp] > nx2)
      ix1[nasp] = nx2;

   iy1[nasp] = (int)(y1/(*dy) + 0.5) + 1;
   if(iy1[nasp] < 0)
      iy1[nasp] = 0;
   if(iy1[nasp] > ny2)
      iy1[nasp] = ny2;

   fprintf(stderr,"\tslip_fac= %.2f trise_fac= %.2f rvf_fac= %.2f\n",slip_fac[nasp],trise_fac[nasp],rvf_fac[nasp]);
   fprintf(stderr,"\tix0= %d ix1= %d iy0= %d iy1= %d a= %.2f km*km\n",ix0[nasp],ix1[nasp],iy0[nasp],iy1[nasp],(ix1[nasp]-ix0[nasp])*(iy1[nasp]-iy0[nasp])*(*dx)*(*dy));

   goflag = 1;
   for(iy=iy0[nasp];iy<iy1[nasp];iy++)
      {
      for(ix=ix0[nasp];ix<ix1[nasp];ix++)
         {
         ip = ix+iy*nx2;

	 if(ps[ip].asp_mask == 1)
	    {
	    goflag = 0;
	    fprintf(stderr,"\tERROR: asperity= %d overlaps with previous asperity, skipping\n",ntot);
	    break;
	    }
         }

      if(goflag == 0)
         break;
      }

   if(goflag)
      {
      for(iy=iy0[nasp];iy<iy1[nasp];iy++)
         {
         for(ix=ix0[nasp];ix<ix1[nasp];ix++)
            {
            ip = ix+iy*nx2;
	    ps[ip].asp_mask = 1;
            }
         }

      sum_num = sum_num - (ix1[nasp]-ix0[nasp])*(iy1[nasp]-iy0[nasp])/(nx2*ny2);
      sum_den = sum_den - slip_fac[nasp]*(ix1[nasp]-ix0[nasp])*(iy1[nasp]-iy0[nasp])/(nx2*ny2);

      savg_asp = savg_asp + slip_fac[nasp]*(ix1[nasp]-ix0[nasp])*(iy1[nasp]-iy0[nasp]);
      area_asp = area_asp + (ix1[nasp]-ix0[nasp])*(iy1[nasp]-iy0[nasp]);

      if(sum_den <= 0.0)
         {
         fprintf(stderr,"\tERROR: asperity scaling too strong, sum_den= %.5e, exiting ...\n",sum_den);
         exit(-1);
         }

      nasp++;
      if(nasp == nalloc)
         {
         nalloc = nalloc + nblock;
         slip_fac = (float *)check_realloc(slip_fac,nalloc*sizeof(float));
         trise_fac = (float *)check_realloc(trise_fac,nalloc*sizeof(float));
         rvf_fac = (float *)check_realloc(rvf_fac,nalloc*sizeof(float));
         ix0 = (int *)check_realloc(ix0,nalloc*sizeof(int));
         ix1 = (int *)check_realloc(ix1,nalloc*sizeof(int));
         iy0 = (int *)check_realloc(iy0,nalloc*sizeof(int));
         iy1 = (int *)check_realloc(iy1,nalloc*sizeof(int));
         }
      }
   }

fclose(fpr);

slip_fac = (float *)check_realloc(slip_fac,nasp*sizeof(float));
trise_fac = (float *)check_realloc(trise_fac,nasp*sizeof(float));
rvf_fac = (float *)check_realloc(rvf_fac,nasp*sizeof(float));
ix0 = (int *)check_realloc(ix0,nasp*sizeof(int));
ix1 = (int *)check_realloc(ix1,nasp*sizeof(int));
iy0 = (int *)check_realloc(iy0,nasp*sizeof(int));
iy1 = (int *)check_realloc(iy1,nasp*sizeof(int));

savg_new = (*savg_old)*sum_num/sum_den;

fprintf(stderr,"\tsavg_old= %13.5e savg_new= %13.5e\n",(*savg_old),savg_new);
fprintf(stderr,"\tnasp= %d\n",nasp);
fprintf(stderr,"\tsavg_asp/savg_new= %.4f\n",savg_asp/area_asp);
fprintf(stderr,"\tarea_asp/area_tot= %.4f\n",area_asp/(nx2*ny2));
fprintf(stderr,"\n");

for(ia=0;ia<nasp;ia++)
   {
   ixt = (int)((*tap_fac)*(ix1[ia]-ix0[ia]) + 0.5);
   if(ixt < 0)
      ixt = 0;

   iyt = (int)((*tap_fac)*(iy1[ia]-iy0[ia]) + 0.5);
   if(iyt < 0)
      iyt = 0;

   /*
   fprintf(stderr,"ixt= %d iyt= %d\n",ixt,iyt);
   */

   for(iy=iy0[ia];iy<iy1[ia];iy++)
      {
      ydamp = 1.0;

      if(iy < iy0[ia]+iyt && iy0[ia] > 1)
         ydamp = (float)(iy+1-iy0[ia])/(float)(iyt);
      else if(iy > iy1[ia]-iyt)
         ydamp = (float)(iy1[ia]-iy)/(float)(iyt);

      if(ydamp < 0.0)
         ydamp = 0.0;
      if(ydamp > 1.0)
         ydamp = 1.0;

      for(ix=ix0[ia];ix<ix1[ia];ix++)
         {
         xdamp = 1.0;

         if(ix < ix0[ia]+ixt)
            xdamp = (float)(ix+1-ix0[ia])/(float)(ixt);
         else if(ix > ix1[ia]-ixt)
            xdamp = (float)(ix1[ia]-ix)/(float)(ixt);

         if(xdamp < 0.0)
            xdamp = 0.0;
         if(xdamp > 1.0)
            xdamp = 1.0;

         ip = ix+iy*nx2;

         ps[ip].asp = xdamp*ydamp*(slip_fac[ia]*savg_new - *savg_old);
         ps[ip].trise_fac = 1.0 + xdamp*ydamp*(trise_fac[ia] - 1.0);
         ps[ip].rvf_fac = 1.0 + xdamp*ydamp*(rvf_fac[ia] - 1.0);
         }
      }
   }
}

void gen_base_c(struct complex *s0,int nx0,int ny0,long *seed)
{
int i, j, ip;
float fzero = 0.0;
float fone = 1.0;

/*
   only do positive half, then use symmetry
*/

for(j=0;j<=ny0/2;j++)
   {
   for(i=0;i<nx0;i++)
      {
      ip = i + j*nx0;

      s0[ip].re = gaus_rand(&fone,&fzero,seed);
      s0[ip].im = gaus_rand(&fone,&fzero,seed);
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

hermit(s0,nx0,ny0);
}

void kfilt_base(struct complex *s0,struct complex *b0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,int kflag,float *lambda_max,float *lambda_min)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
float fre, fim;
float lmax2, lmin2;
int ndkc;
float fnorm, hcoef, beta2;

int ord = 4;

fnorm = 1.0/sqrt(2.0);

hcoef = 2.00; /* H=1.00, hcoef = H + 1 */
hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

beta2 = hcoef;
beta2 = 2.0; /* strictly self-similar */

amp0 = b0[0].re;

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

lmax2 = (*lambda_max)*(*lambda_max);
lmin2 = (*lambda_min)*(*lambda_min);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      /* von Karman */
      if(kflag == MAI_FLAG || kflag == SUZUKI_FLAG || kflag == MAI_SOMERVILLE_HYBRID_FLAG)
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == FRANKEL_FLAG) /* frankel(2009) scaling */
         {
	 if(amp < 1.0)
	    fac = amp0;
         else
            fac = amp0*exp(-0.5*beta2*log(amp));
	 }

      k2 = kx*kx + ky*ky;
      if(k2 > 0.0)
         fac = fac/((1.0 + exp(ord*log(k2*lmin2)))*(1.0 + exp(-ord*log(k2*lmax2))));

      s0[ip].re = fnorm*fac*b0[ip].re;
      s0[ip].im = fnorm*fac*b0[ip].im;
      }
   }

/* 
   Enforce DC & Nyquists to be real
*/

s0[0].re = amp0;
s0[0].im = 0.0;
s0[nx0/2].im = 0.0;
s0[nx0*ny0/2].im = 0.0;
s0[nx0/2+nx0*ny0/2].im = 0.0;

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

hermit(s0,nx0,ny0);
}

void init_rvf_IO(struct pointsource *ps,float *rvfrac0,int nx,int ny,float *dx,float *dy,char *file,float *tap_fac)
{
FILE *fpr;
char str[1024];
float *rvf_fac, x0, y0, x1, y1;
int *ix0, *iy0, *ix1, *iy1;
int ix, iy, ip, ip1;
int nr;

int goflag, ia, nasp, nalloc, ntot;
int nblock = 100;

int ixt, iyt;
float xdamp, ydamp;

nalloc = nblock;
rvf_fac = (float *)check_malloc(nalloc*sizeof(float));
ix0 = (int *)check_malloc(nalloc*sizeof(int));
ix1 = (int *)check_malloc(nalloc*sizeof(int));
iy0 = (int *)check_malloc(nalloc*sizeof(int));
iy1 = (int *)check_malloc(nalloc*sizeof(int));

fprintf(stderr,"rvf params:\n\n");

fpr = fopfile(file,"r");

ntot = 0;
nasp = 0;
while(fgets(str,1024,fpr) != NULL)
   {
   nr = sscanf(str,"%f %f %f %f %f",&rvf_fac[nasp],&x0,&y0,&x1,&y1);
   ntot++;

   ix0[nasp] = (int)((0.5*nx*(*dx) + x0)/(*dx) + 0.5);
   if(ix0[nasp] < 0)
      ix0[nasp] = 0;
   if(ix0[nasp] > nx)
      ix0[nasp] = nx;

   iy0[nasp] = (int)(y0/(*dy) + 0.5);
   if(iy0[nasp] < 0)
      iy0[nasp] = 0;
   if(iy0[nasp] > ny)
      iy0[nasp] = ny;

   ix1[nasp] = (int)((0.5*nx*(*dx) + x1)/(*dx) + 0.5);
   if(ix1[nasp] < 0)
      ix1[nasp] = 0;
   if(ix1[nasp] > nx)
      ix1[nasp] = nx;

   iy1[nasp] = (int)(y1/(*dy) + 0.5) + 1;
   if(iy1[nasp] < 0)
      iy1[nasp] = 0;
   if(iy1[nasp] > ny)
      iy1[nasp] = ny;

   fprintf(stderr,"\tix0= %d ix1= %d iy0= %d iy1= %d rvf= %.6f\n",ix0[nasp],ix1[nasp],iy0[nasp],iy1[nasp],rvf_fac[nasp]);

   nasp++;
   if(nasp == nalloc)
      {
      nalloc = nalloc + nblock;
      rvf_fac = (float *)check_realloc(rvf_fac,nalloc*sizeof(float));
      ix0 = (int *)check_realloc(ix0,nalloc*sizeof(int));
      ix1 = (int *)check_realloc(ix1,nalloc*sizeof(int));
      iy0 = (int *)check_realloc(iy0,nalloc*sizeof(int));
      iy1 = (int *)check_realloc(iy1,nalloc*sizeof(int));
      }
   }

fclose(fpr);

rvf_fac = (float *)check_realloc(rvf_fac,nasp*sizeof(float));
ix0 = (int *)check_realloc(ix0,nasp*sizeof(int));
ix1 = (int *)check_realloc(ix1,nasp*sizeof(int));
iy0 = (int *)check_realloc(iy0,nasp*sizeof(int));
iy1 = (int *)check_realloc(iy1,nasp*sizeof(int));

for(ia=0;ia<nasp;ia++)
   {
   ixt = (int)((*tap_fac)*(ix1[ia]-ix0[ia]) + 0.5);
   if(ixt < 0)
      ixt = 0;

   iyt = (int)((*tap_fac)*(iy1[ia]-iy0[ia]) + 0.5);
   if(iyt < 0)
      iyt = 0;

   /*
   fprintf(stderr,"ixt= %d iyt= %d\n",ixt,iyt);
   */

   for(iy=iy0[ia];iy<iy1[ia];iy++)
      {
      ydamp = 1.0;

      if(iy < iy0[ia]+iyt && iy0[ia] > 1)
         ydamp = (float)(iy+1-iy0[ia])/(float)(iyt);
      else if(iy > iy1[ia]-iyt)
         ydamp = (float)(iy1[ia]-iy)/(float)(iyt);

      if(ydamp < 0.0)
         ydamp = 0.0;
      if(ydamp > 1.0)
         ydamp = 1.0;

      for(ix=ix0[ia];ix<ix1[ia];ix++)
         {
         xdamp = 1.0;

         if(ix < ix0[ia]+ixt)
            xdamp = (float)(ix+1-ix0[ia])/(float)(ixt);
         else if(ix > ix1[ia]-ixt)
            xdamp = (float)(ix1[ia]-ix)/(float)(ixt);

         if(xdamp < 0.0)
            xdamp = 0.0;
         if(xdamp > 1.0)
            xdamp = 1.0;

         ip = ix+iy*nx;

         ps[ip].rvf = (*rvfrac0)*(1.0 - xdamp*ydamp*(1.0-rvf_fac[ia]));
         }
      }
   }
}
