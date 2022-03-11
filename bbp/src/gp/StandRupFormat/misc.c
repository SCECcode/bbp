#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

struct complex
   {
   float re;
   float im;
   };

void get_filtvals(char *cbuf,float *flo,float *fhi,int *ord,int *phs);
void tfilter(float *y, float *dt, int nt, int order,
             float *fhi, float *flo, int phase, float *space);
void cxzero(struct complex *p, int n);
void hp_filter(struct complex *q, struct complex *p,
               int n, float *alpha, float *beta, int sgn);
void lp_filter(struct complex *q, struct complex *p,
               int n, float *alpha, float *beta, int sgn);

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

double gaussian_rand(float *sigma,float *mean,long *seed)
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

int write_xyz(char *file,struct standrupformat *srf,char *type,int ig,int calc_xy,float *xoff,float *yoff,int tsstr,int tsend,int tsinc,float *svmin,float *slipmin,int keepsign,int dump_slip)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat;
float sn, se, arg;
float len, wid, dx, dy, stk, dip, dtop, shypo, dhypo;
float len2, slip, xx, yy, outval, slipval;
float sv, svmax, *stf1, *stf2, *stf3;
float ss, dd, dfac;
int it, maxnt, ts, itdel;
int nx, ny;
int nt1, nt2, nt3;
int i, j, k, npskip;
int igst, ignd, igp;
int nd, place;
char str[512], frmt[32];

float sum1, sum2, sum3;

float flo, fhi, mindt, *space, *sf1, *sf2, *sf3;
int itsh, ord, phs;

float dtdx, dtdy;
int stride, jj;

float a1, a2, a3, a4;
float xr, yr, rr, rhypo_tol;
int astride;

int isgn = 1;

float rperd = 0.017453293;

float *rlon, *rlat, *rdep, zz, x1, y1, x2, y2;
float cosD, sinD;
float erad = 6378.139;
float mrot = -90.0;
double g0 = 0.0;
double b0 = 0.0;
double amat[9], ainv[9];

stride = 3;
stride = 5;
stride = 40;
stride = 4;
stride = 1;

/*

f'(x) = (1/dx)*{(252/3491)*[8*R1 - (1/2)*R2] - (1/6982)*[(64/3)*R3 - (27/4)*R4]} + O(dx^6)

R1 = f(x+1*dx) - f(x-1*dx)
R2 = f(x+2*dx) - f(x-2*dx)
R3 = f(x+3*dx) - f(x-3*dx)
R4 = f(x+4*dx) - f(x-4*dx)

f'(x) = (1/dx)*(a1*R1 + a2*R2 + a3*R3 + a4*R4)

a1 = (252/3491)*8 = 2016/3491;     0.57748496132...
a2 = (252/3491)*(1/2) = -126/3491;     0.03609281008...
a3 = (1/6982)*(64/3) = -32/10473;     0.00305547598587...
a4 = (1/6982)*(27/4) = 27/27928;     0.00096677169865...

*/

a1 = 2016.0/3491.0;
a2 = -126.0/3491.0;
a3 = -32.0/10473.0;
a4 = 27.0/27928.0;
astride = 4;
astride = 1;
astride = 1;
astride = 40;

a1 = 0.5;
a2 = 0.0;
a3 = 0.0;
a4 = 0.0;
/*
*/

/* don't compute secant rupture speed within rhypo_tol of hypocenter */
rhypo_tol = 0.0;
rhypo_tol = 1.0;

if(strcmp(srf[0].type,"PLANE") != 0 || ig >= srf[0].srf_prect.nseg)
   return(-99);

if(ig < 0)
   {
   igst = 0;
   ignd = srf[0].srf_prect.nseg;
   }
else
   {
   igst = ig;
   ignd = igst + 1;
   }

if(strncmp(type,"filter_velocity",15) == 0)
   {
   space = NULL;
   sf1 = NULL;
   sf2 = NULL;
   sf3 = NULL;
   }

if(strcmp(type,"rupture") == 0)
   {
   if(tsstr < 0)
      tsstr = 0;
   if(tsend < 0)
      {
      apval_ptr = srf[0].srf_apnts.apntvals;
      for(k=0;k<srf[0].srf_apnts.np;k++)
         {
	 itdel = (int)(apval_ptr[k].tinit/apval_ptr[k].dt + 0.5);

	 if(tsend < itdel + apval_ptr[k].nt1)
            tsend = itdel + apval_ptr[k].nt1;
	 if(tsend < itdel + apval_ptr[k].nt2)
            tsend = itdel + apval_ptr[k].nt2;
	 if(tsend < itdel + apval_ptr[k].nt3)
            tsend = itdel + apval_ptr[k].nt3;
	 }
      }

   nd = 0;
   place = 1;
   while((tsend-tsstr)/tsinc > (place-1))
      {
      place = place*10;
      nd++;
      }

   sprintf(frmt,"%s%%.%dd",file,nd);
   }
else
   {
   tsstr = 0;
   tsend = 1;
   }

for(ts=tsstr;ts<tsend;ts=ts+tsinc)
   {
   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      {
      if(strcmp(type,"rupture") == 0)
         {
         sprintf(str,frmt,(ts-tsstr)/tsinc);
         fpw = fopfile(str,"w");
         }
      else
         fpw = fopfile(file,"w");
      }

   for(igp=igst;igp<ignd;igp++)
      {
      prseg_ptr = srf[0].srf_prect.prectseg;

      elon = prseg_ptr[igp].elon;
      elat = prseg_ptr[igp].elat;
      nx = prseg_ptr[igp].nstk;
      ny = prseg_ptr[igp].ndip;
      len = prseg_ptr[igp].flen;
      wid = prseg_ptr[igp].fwid;
      stk = prseg_ptr[igp].stk;
      dip = prseg_ptr[igp].dip;
      dtop = prseg_ptr[igp].dtop;
      shypo = prseg_ptr[igp].shyp;
      dhypo = prseg_ptr[igp].dhyp;

      len2 = 0.5*len;
      dx = len/(float)(nx);
      dy = wid/(float)(ny);

      npskip = 0;
      for(i=0;i<igp;i++)
         npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

      apval_ptr = srf[0].srf_apnts.apntvals + npskip;

      if(strncmp(type,"filter_velocity",15) == 0)
         {
	 /* XXXX need to change */
	 flo = 1.0e+15;
	 fhi = 0.0;
	 ord = 4;
	 phs = 0;

	 get_filtvals(type+15,&flo,&fhi,&ord,&phs);

/*
fprintf(stderr,"flo= %.5e fhi= %.5e ord= %d phs= %d\n",flo,fhi,ord,phs);
*/

	 maxnt = 0;
	 mindt = 1.0e+15;
         for(j=0;j<ny;j++)
            {
            for(i=0;i<nx;i++)
               {
	       k = i + j*nx;

	       if(apval_ptr[k].nt1 > maxnt)
	          maxnt = apval_ptr[k].nt1;
	       if(apval_ptr[k].nt2 > maxnt)
	          maxnt = apval_ptr[k].nt2;
	       if(apval_ptr[k].nt3 > maxnt)
	          maxnt = apval_ptr[k].nt3;

	       if(apval_ptr[k].dt < mindt)
	          mindt = apval_ptr[k].dt;
               }
            }

	 itsh = 1.0/(mindt*flo);
	 maxnt = maxnt + 2*itsh;

	 space = (float *) check_realloc (space,4*maxnt*sizeof(float));
	 sf1 = (float *) check_realloc (sf1,maxnt*sizeof(float));
	 sf2 = (float *) check_realloc (sf2,maxnt*sizeof(float));
	 sf3 = (float *) check_realloc (sf3,maxnt*sizeof(float));
         }

      if(strcmp(type,"rough") == 0)
         {
         rlon = (float *)check_malloc(nx*ny*sizeof(float));
         rlat = (float *)check_malloc(nx*ny*sizeof(float));
         rdep = (float *)check_malloc(nx*ny*sizeof(float));

         mrot = stk - 90.0;
         gen_matrices(amat,ainv,&mrot,&elon,&elat);

         cosD = cos(dip*rperd);
         sinD = sin(dip*rperd);

         for(j=0;j<ny;j++)
            {
            dd = (j+0.5)*dy;
            yy = dd*cosD;
            zz = dd*sinD + dtop;
            for(i=0;i<nx;i++)
               {
	       k = i + j*nx;
               xx = (i+0.5)*dx - len2;

               gcproj(&xx,&yy,&(rlon[k]),&(rlat[k]),&erad,&g0,&b0,amat,ainv,0);
               rdep[k] = zz;
               }
            }
         }

      for(j=0;j<ny;j++)
         {
         yy = (j+0.5)*dy;
         for(i=0;i<nx;i++)
            {
            xx = (i+0.5)*dx;
	    k = i + j*nx;

            if(calc_xy)
               {
               set_ne(&(prseg_ptr[igp].elon),&(prseg_ptr[igp].elat),&(apval_ptr[k].lon),&(apval_ptr[k].lat),&sn,&se);

               arg = (apval_ptr[k].stk)*rperd;
               xx = se*sin(arg) + sn*cos(arg) + 0.5*prseg_ptr[igp].flen;

               arg = (apval_ptr[k].dip)*rperd;
               yy = ((apval_ptr[k].dep)-(prseg_ptr[igp].dtop))/sin(arg);
               }

	    xx = xx + *xoff;
	    yy = yy + *yoff;

            if(strcmp(type,"slip") == 0 || dump_slip)
               {
	       slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	                apval_ptr[k].slip2*apval_ptr[k].slip2 +
	                apval_ptr[k].slip3*apval_ptr[k].slip3);

	       if(keepsign)
		  {
		  isgn = 1;

		  if(apval_ptr[k].slip1 < 0.0)
		     isgn = -1;
		  if((apval_ptr[k].slip2*apval_ptr[k].slip2) > (apval_ptr[k].slip1*apval_ptr[k].slip1) && apval_ptr[k].slip2 < 0.0)
		     isgn = -1;
		  if((apval_ptr[k].slip3*apval_ptr[k].slip3) > (apval_ptr[k].slip1*apval_ptr[k].slip1) && (apval_ptr[k].slip3*apval_ptr[k].slip3) > (apval_ptr[k].slip2*apval_ptr[k].slip2) && apval_ptr[k].slip3 < 0.0)
		     isgn = -1;

	          slip = isgn*slip;
		  }

               outval = slip;
               slipval = slip;
	       }

            if(strcmp(type,"secant_speed") == 0)
               {
	       xr = (xx-len2) - shypo;
	       yr = yy - dhypo;
	       xr = ((xx-*xoff)-len2) - shypo;
	       yr = (yy-*yoff) - dhypo;

	       rr = sqrt(xr*xr + yr*yr);

	       if(rr > rhypo_tol && apval_ptr[k].tinit > 0.0)
	          outval = sqrt(xr*xr + yr*yr)/apval_ptr[k].tinit;
	       else
	          outval = 0.0;
	       }

            if(strcmp(type,"secant_speed_frac") == 0)
               {
	       xr = (xx-len2) - shypo;
	       yr = yy - dhypo;
	       xr = ((xx-*xoff)-len2) - shypo;
	       yr = (yy-*yoff) - dhypo;

	       rr = sqrt(xr*xr + yr*yr);

	       if(rr > rhypo_tol && apval_ptr[k].tinit > 0.0)
	          outval = sqrt(xr*xr + yr*yr)/apval_ptr[k].tinit;
	       else
	          outval = 0.0;

	       outval = 1.0e+05*outval/apval_ptr[k].vs;
	       }

            if(strcmp(type,"rupture_speed") == 0)
               {
	       if(i >= stride && i < nx-stride && j >= stride && j < ny-stride)
	          {
		  dtdx = (apval_ptr[k+stride].tinit - apval_ptr[k-stride].tinit)/(2.0*stride*dx);
		  dtdy = (apval_ptr[k+stride*nx].tinit - apval_ptr[k-stride*nx].tinit)/(2.0*stride*dy);

		  outval = sqrt(dtdx*dtdx + dtdy*dtdy);

		  if(outval > 0.0)
		     outval = 1.0/outval;
		  else
		     outval = 0.0;
		  }
	       else
                  outval = 0.0;
	       }

            if(strcmp(type,"rupture_speed_frac") == 0)
               {
	       if(i >= 4*astride && i < nx-4*astride && j >= 4*astride && j < ny-4*astride)
	          {
/*
		  dtdx = 0.0;
		  dtdy = 0.0;
		  for(jj=1;jj<=stride;jj++)
		     {
		     dtdx = dtdx + (apval_ptr[k+jj].tinit - apval_ptr[k-jj].tinit)/(2.0*jj*dx);
		     dtdy = dtdy + (apval_ptr[k+jj*nx].tinit - apval_ptr[k-jj*nx].tinit)/(2.0*jj*dy);
		     }
		  dtdx = dtdx/stride;
		  dtdy = dtdy/stride;

		  dtdx = (apval_ptr[k+stride].tinit - apval_ptr[k-stride].tinit)/(2.0*stride*dx);
		  dtdy = (apval_ptr[k+stride*nx].tinit - apval_ptr[k-stride*nx].tinit)/(2.0*stride*dy);
*/

		  dtdx = (1/dx)*(a1*(apval_ptr[k+1*astride].tinit - apval_ptr[k-1*astride].tinit)
		               + a2*(apval_ptr[k+2*astride].tinit - apval_ptr[k-2*astride].tinit)
		               + a3*(apval_ptr[k+3*astride].tinit - apval_ptr[k-3*astride].tinit)
		               + a4*(apval_ptr[k+4*astride].tinit - apval_ptr[k-4*astride].tinit));

		  dtdy = (1/dy)*(a1*(apval_ptr[k+1*nx*astride].tinit - apval_ptr[k-1*nx*astride].tinit)
		               + a2*(apval_ptr[k+2*nx*astride].tinit - apval_ptr[k-2*nx*astride].tinit)
		               + a3*(apval_ptr[k+3*nx*astride].tinit - apval_ptr[k-3*nx*astride].tinit)
		               + a4*(apval_ptr[k+4*nx*astride].tinit - apval_ptr[k-4*nx*astride].tinit));

		  outval = sqrt(dtdx*dtdx + dtdy*dtdy);

		  if(outval > 0.0)
		     outval = astride/outval;
		  else
		     outval = 9999.0;
		  }
	       else
                  outval = 0.0;

	       outval = 1.0e+05*outval/apval_ptr[k].vs;
	       }

            if(strcmp(type,"integrate_sliprate") == 0)
               {
               stf1 = apval_ptr[k].stf1;
               stf2 = apval_ptr[k].stf2;
               stf3 = apval_ptr[k].stf3;

	       sum1 = 0.0;
               for(it=0;it<apval_ptr[k].nt1;it++)
	          sum1 = sum1 + apval_ptr[k].dt*stf1[it];

	       sum2 = 0.0;
               for(it=0;it<apval_ptr[k].nt2;it++)
	          sum2 = sum2 + apval_ptr[k].dt*stf2[it];

	       sum3 = 0.0;
               for(it=0;it<apval_ptr[k].nt3;it++)
	          sum3 = sum3 + apval_ptr[k].dt*stf3[it];

               outval = sqrt(sum1*sum1 + sum2*sum2 + sum3*sum3);
	       }

            if(strcmp(type,"velocity") == 0)
               {
               stf1 = apval_ptr[k].stf1;
               stf2 = apval_ptr[k].stf2;
               stf3 = apval_ptr[k].stf3;

	       maxnt = apval_ptr[k].nt1;
	       if(apval_ptr[k].nt2 > maxnt)
	          maxnt = apval_ptr[k].nt2;
	       if(apval_ptr[k].nt3 > maxnt)
	          maxnt = apval_ptr[k].nt3;

	       svmax = 0.0;
               for(it=0;it<maxnt;it++)
	          {
	          sv = 0.0;
	          if(it < apval_ptr[k].nt1)
	             sv = sv + stf1[it]*stf1[it];
	          if(it < apval_ptr[k].nt2)
	             sv = sv + stf2[it]*stf2[it];
	          if(it < apval_ptr[k].nt3)
	             sv = sv + stf3[it]*stf3[it];

	          if(sv > svmax)
	             svmax = sv;
	          }

               outval = sqrt(svmax);
	       }

            if(strncmp(type,"filter_velocity",15) == 0)
               {
	       maxnt = apval_ptr[k].nt1;
	       if(apval_ptr[k].nt2 > maxnt)
	          maxnt = apval_ptr[k].nt2;
	       if(apval_ptr[k].nt3 > maxnt)
	          maxnt = apval_ptr[k].nt3;

	       maxnt = maxnt + 2*itsh;
               for(it=0;it<maxnt;it++)
	          {
                  sf1[it] = 0.0;
                  sf2[it] = 0.0;
                  sf3[it] = 0.0;
		  }

               stf1 = apval_ptr[k].stf1;
               for(it=0;it<apval_ptr[k].nt1;it++)
                  sf1[it+itsh] = stf1[it];
               tfilter(sf1,&(apval_ptr[k].dt),maxnt,ord,&fhi,&flo,phs,space);

               stf2 = apval_ptr[k].stf2;
               for(it=0;it<apval_ptr[k].nt2;it++)
                  sf2[it+itsh] = stf2[it];
               tfilter(sf2,&(apval_ptr[k].dt),maxnt,ord,&fhi,&flo,phs,space);

               stf3 = apval_ptr[k].stf3;
               for(it=0;it<apval_ptr[k].nt3;it++)
                  sf3[it+itsh] = stf3[it];
               tfilter(sf3,&(apval_ptr[k].dt),maxnt,ord,&fhi,&flo,phs,space);

	       svmax = 0.0;
               for(it=0;it<maxnt;it++)
	          {
	          sv = sf1[it]*sf1[it] + sf2[it]*sf2[it] + sf3[it]*sf3[it];

	          if(sv > svmax)
	             svmax = sv;
	          }

               outval = sqrt(svmax);
	       }

            if(strcmp(type,"rupture") == 0)
               {
               stf1 = apval_ptr[k].stf1;
               stf2 = apval_ptr[k].stf2;
               stf3 = apval_ptr[k].stf3;

	       nt1 = apval_ptr[k].nt1;
	       while(nt1 > 0 && stf1[nt1-1] < *svmin)
	          nt1--;

	       nt2 = apval_ptr[k].nt2;
	       while(nt2 > 0 && stf2[nt2-1] < *svmin)
	          nt2--;

	       nt3 = apval_ptr[k].nt3;
	       while(nt3 > 0 && stf3[nt3-1] < *svmin)
	          nt3--;

	       itdel = ts - (int)(apval_ptr[k].tinit/apval_ptr[k].dt + 0.5);

	       if(itdel < 1)
	          sv = -1.0e+05;
	       else if(itdel < nt1 || itdel < nt2 || itdel < nt3)
	          {
	          sv = 0.0;
	          if(itdel < nt1)
	             sv = sv + stf1[itdel]*stf1[itdel];
	          if(itdel < nt2)
	             sv = sv + stf2[itdel]*stf2[itdel];
	          if(itdel < nt3)
	             sv = sv + stf3[itdel]*stf3[itdel];

                  sv = sqrt(sv);
		  }
	       else
	          {
		  maxnt = nt1;
		  if(nt2 > maxnt)
		     maxnt = nt2;
		  if(nt3 > maxnt)
		     maxnt = nt3;

		  dfac = 1.0;
		  if((itdel - maxnt) < maxnt)
		     dfac = (float)(itdel-maxnt)/(float)(maxnt);

	          sv = -dfac*sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	                         apval_ptr[k].slip2*apval_ptr[k].slip2 +
	                         apval_ptr[k].slip3*apval_ptr[k].slip3);

		  if(sv > -(*slipmin))
		     sv = -(*slipmin);
		  }

               outval = sv;
	       }

            if(strcmp(type,"tinit") == 0)
               outval = apval_ptr[k].tinit;

            if(strcmp(type,"trise") == 0)
	       {
               outval = apval_ptr[k].nt1*apval_ptr[k].dt;
	       if(apval_ptr[k].nt2*apval_ptr[k].dt > outval)
                  outval = apval_ptr[k].nt2*apval_ptr[k].dt;
	       if(apval_ptr[k].nt3*apval_ptr[k].dt > outval)
                  outval = apval_ptr[k].nt3*apval_ptr[k].dt;
	       }

            if(strcmp(type,"strike") == 0)
               outval = apval_ptr[k].stk;

            if(strcmp(type,"dip") == 0)
               outval = apval_ptr[k].dip;

            if(strcmp(type,"rake") == 0)
	       {
	       ss = cos(rperd*apval_ptr[k].rake)*apval_ptr[k].slip1 - sin(rperd*apval_ptr[k].rake)*apval_ptr[k].slip2;
	       dd = sin(rperd*apval_ptr[k].rake)*apval_ptr[k].slip1 + cos(rperd*apval_ptr[k].rake)*apval_ptr[k].slip2;

	       outval = 90.0;
	       if(ss != 0.0)
	          outval = atan(dd/ss)/rperd;
               if(ss < 0.0)
	          outval = outval + 180.0;

	       if(ss == 0.0 && dd == 0.0)
                  outval = apval_ptr[k].rake;

	       while(outval >= 360.0)
	          outval = outval - 360.0;
	       while(outval < 0.0)
	          outval = outval + 360.0;
	       }

            if(strcmp(type,"rough") == 0)
	       {
               gen_matrices(amat,ainv,&mrot,&(rlon[k]),&(rlat[k]));
               gcproj(&x2,&y2,&(apval_ptr[k].lon),&(apval_ptr[k].lat),&erad,&g0,&b0,amat,ainv,1);

	       zz = apval_ptr[k].dep - rdep[k];

	       outval = -y2*sinD + zz*cosD;
	       }

            if(dump_slip)
               fprintf(fpw,"%12.5e %12.5e %12.5e %12.5e\n",xx,yy,outval,slipval);
	    else
               fprintf(fpw,"%12.5e %12.5e %12.5e\n",xx,yy,outval);
            }
         }

      *xoff = *xoff + len;
      }
   fclose(fpw);
   }

if(strncmp(type,"filter_velocity",15) == 0)
   {
   free(space);
   free(sf1);
   free(sf2);
   free(sf3);
   }

if(strcmp(type,"rough") == 0)
   {
   free(rlon);
   free(rlat);
   free(rdep);
   }
}

void write_maxsvfX(char *file,struct standrupformat *srf,char *type,int ig,float *maxslip)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat;
float sn, se, arg;
float len, wid, dx, dy, stk, dip, dtop, shypo, dhypo;
float len2, slip, xx, yy;
float *stf;
int it, kp, nx, ny;
int i, j, k, npskip, pflag;

float rperd = 0.017453293;

if(strcmp(srf[0].type,"PLANE") == 0 && ig < srf[0].srf_prect.nseg)
   {
   prseg_ptr = srf[0].srf_prect.prectseg;

   elon = prseg_ptr[ig].elon;
   elat = prseg_ptr[ig].elat;
   nx = prseg_ptr[ig].nstk;
   ny = prseg_ptr[ig].ndip;
   len = prseg_ptr[ig].flen;
   wid = prseg_ptr[ig].fwid;
   stk = prseg_ptr[ig].stk;
   dip = prseg_ptr[ig].dip;
   dtop = prseg_ptr[ig].dtop;
   shypo = prseg_ptr[ig].shyp;
   dhypo = prseg_ptr[ig].dhyp;

   len2 = 0.5*len;
   dx = len/(float)(nx);
   dy = wid/(float)(ny);

   npskip = 0;
   for(i=0;i<ig;i++)
      npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

   pflag = 1;
   }
else
   {
   elon = -99;
   elat = -99;
   pflag = 0;

   pflag = 1;
   npskip = 0;

   nx = srf[0].srf_apnts.np;
   ny = 1;
   }

*maxslip = 0.0;
if(pflag)
   {
   apval_ptr = srf[0].srf_apnts.apntvals + npskip;

fprintf(stderr,"nx= %5d ny= %d\n",nx,ny);

   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
	 k = i + j*nx;

         if(strcmp(type,"slip") == 0)
            {
	    slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	             apval_ptr[k].slip2*apval_ptr[k].slip2 +
	             apval_ptr[k].slip3*apval_ptr[k].slip3);

	    if(slip > *maxslip)
	       {
	       *maxslip = slip;
	       kp = k;
	       }
	    }
         }
      }

fprintf(stderr,"kp= %5d\n",kp);

   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(file,"w");

   fprintf(fpw,"svf svf\n");
   fprintf(fpw,"%d %13.5e\n",apval_ptr[kp].nt1,apval_ptr[kp].dt);

   stf = apval_ptr[kp].stf1;
   for(it=0;it<apval_ptr[kp].nt1;it++)
      fprintf(fpw,"%13.5e\n",stf[it]);

   fclose(fpw);
   }
}

void write_maxsvf(char *file,struct standrupformat *srf,char *type,int ig,float *maxslip,int kp_in)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat;
float sn, se, arg;
float len, wid, dx, dy, stk, dip, dtop, shypo, dhypo;
float len2, slip, xx, yy, dels;
float *stfp, *stfout;
int it, kp, nx, ny, ntout;
int i, j, k, npskip, pflag;

float rperd = 0.017453293;

if(strcmp(srf[0].type,"PLANE") == 0 && ig < srf[0].srf_prect.nseg)
   {
   prseg_ptr = srf[0].srf_prect.prectseg;

   elon = prseg_ptr[ig].elon;
   elat = prseg_ptr[ig].elat;
   nx = prseg_ptr[ig].nstk;
   ny = prseg_ptr[ig].ndip;
   len = prseg_ptr[ig].flen;
   wid = prseg_ptr[ig].fwid;
   stk = prseg_ptr[ig].stk;
   dip = prseg_ptr[ig].dip;
   dtop = prseg_ptr[ig].dtop;
   shypo = prseg_ptr[ig].shyp;
   dhypo = prseg_ptr[ig].dhyp;

   len2 = 0.5*len;
   dx = len/(float)(nx);
   dy = wid/(float)(ny);

   npskip = 0;
   for(i=0;i<ig;i++)
      npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

   pflag = 1;
   }
else
   {
   elon = -99;
   elat = -99;
   pflag = 0;

   pflag = 1;
   npskip = 0;

   nx = srf[0].srf_apnts.np;
   ny = 1;
   }

dels = 1.0e+15;
if(pflag)
   {
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
	 k = i + j*nx;

         apval_ptr = &(srf[0].srf_apnts.apntvals[k+npskip]);

         if(strcmp(type,"slip") == 0)
            {
	    slip = sqrt(apval_ptr->slip1*apval_ptr->slip1 +
	             apval_ptr->slip2*apval_ptr->slip2 +
	             apval_ptr->slip3*apval_ptr->slip3);

	    if(*maxslip < 0.0)
	       {
	       if((slip + *maxslip)*(slip + *maxslip) < dels)
	          {
	          kp = k;
		  dels = (slip + *maxslip)*(slip + *maxslip);

fprintf(stderr,"kp= %d dels= %13.5e\n",kp,dels);
		  }
	       }
	    else if(slip > *maxslip)
	       {
	       *maxslip = slip;
	       kp = k;
	       }
	    }
         }
      }

   if(*maxslip < 0.0)
      *maxslip = -(*maxslip);

   if(kp_in >= 0)
      kp = kp_in;

   apval_ptr = &(srf[0].srf_apnts.apntvals[kp+npskip]);

   ntout = 0;
   if(apval_ptr->nt1 > ntout)
      ntout = apval_ptr->nt1;
   if(apval_ptr->nt2 > ntout)
      ntout = apval_ptr->nt2;
   if(apval_ptr->nt3 > ntout)
      ntout = apval_ptr->nt3;

   stfout = (float *)check_malloc(ntout*sizeof(float));

   for(it=0;it<ntout;it++)
      stfout[it] = 0.0;

   if(apval_ptr->nt1)
      {
      stfp = apval_ptr->stf1;
      for(it=0;it<apval_ptr->nt1;it++)
         stfout[it] = stfout[it] + stfp[it]*stfp[it];
      }

   if(apval_ptr->nt2)
      {
      stfp = apval_ptr->stf2;
      for(it=0;it<apval_ptr->nt2;it++)
         stfout[it] = stfout[it] + stfp[it]*stfp[it];
      }

   if(apval_ptr->nt3)
      {
      stfp = apval_ptr->stf3;
      for(it=0;it<apval_ptr->nt3;it++)
         stfout[it] = stfout[it] + stfp[it]*stfp[it];
      }

   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(file,"w");

   fprintf(fpw,"svf svf kp= %d slip= %.2f\n",kp,*maxslip);
   fprintf(fpw,"%d %13.5e\n",ntout,apval_ptr->dt);

   for(it=0;it<ntout;it++)
      fprintf(fpw,"%13.5e\n",sqrt(stfout[it]));

   fclose(fpw);
   free(stfout);
   }
}

void get_vmax2slip(char *file,struct standrupformat *srf,char *type,int ig)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat;
float sn, se, arg;
float len, wid, dx, dy, stk, dip, dtop, shypo, dhypo;
float len2, slip, xx, yy, dels;
float *stfp, *stfout;
int it, kp, nx, ny, ntout;
int i, j, k, npskip, pflag;
float vos, maxv;
int nc;

float rperd = 0.017453293;

if(strcmp(file,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(file,"w");

if(strcmp(srf[0].type,"PLANE") == 0 && ig < srf[0].srf_prect.nseg)
   {
   prseg_ptr = srf[0].srf_prect.prectseg;

   elon = prseg_ptr[ig].elon;
   elat = prseg_ptr[ig].elat;
   nx = prseg_ptr[ig].nstk;
   ny = prseg_ptr[ig].ndip;
   len = prseg_ptr[ig].flen;
   wid = prseg_ptr[ig].fwid;
   stk = prseg_ptr[ig].stk;
   dip = prseg_ptr[ig].dip;
   dtop = prseg_ptr[ig].dtop;
   shypo = prseg_ptr[ig].shyp;
   dhypo = prseg_ptr[ig].dhyp;

   len2 = 0.5*len;
   dx = len/(float)(nx);
   dy = wid/(float)(ny);

   npskip = 0;
   for(i=0;i<ig;i++)
      npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

   pflag = 1;
   }
else
   {
   elon = -99;
   elat = -99;
   pflag = 0;

   pflag = 1;
   npskip = 0;

   nx = srf[0].srf_apnts.np;
   ny = 1;
   }

dels = 1.0e+15;
if(pflag)
   {
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
	 k = i + j*nx;

         apval_ptr = &(srf[0].srf_apnts.apntvals[k+npskip]);

         nc = 0;
         vos = 0.0;
         if(apval_ptr->nt1)
            {
            stfp = apval_ptr->stf1;
            maxv = 0.0;
            for(it=0;it<apval_ptr->nt1;it++)
               {
               if(stfp[it]*stfp[it] > maxv)
                  maxv = stfp[it]*stfp[it];
	       }
            vos = vos + sqrt(maxv/(apval_ptr->slip1*apval_ptr->slip1));
            nc++;
            }

         if(apval_ptr->nt2)
            {
            stfp = apval_ptr->stf2;
            maxv = 0.0;
            for(it=0;it<apval_ptr->nt2;it++)
               {
               if(stfp[it]*stfp[it] > maxv)
                  maxv = stfp[it]*stfp[it];
	       }
            vos = vos + sqrt(maxv/(apval_ptr->slip2*apval_ptr->slip2));
            nc++;
            }

         if(apval_ptr->nt3)
            {
            stfp = apval_ptr->stf3;
            maxv = 0.0;
            for(it=0;it<apval_ptr->nt3;it++)
               {
               if(stfp[it]*stfp[it] > maxv)
                  maxv = stfp[it]*stfp[it];
	       }
            vos = vos + sqrt(maxv/(apval_ptr->slip3*apval_ptr->slip3));
            nc++;
            }

         if(nc)
            fprintf(fpw,"%13.5e\n",vos/(float)(nc));

         }
      }

   fclose(fpw);
   }
}

void read_Fvelmodel(char *vfile,struct velmodel *vm)
{
FILE *fpr, *fopfile();
int i, nr;
char str[512];

int nblock = 50;

if(strcmp(vfile,"NOT_PROVIDED") == 0)
   {
   vm->nlay = 1;
   vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
   vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

   vm->th[0] = 9999.9;
   vm->vp[0] = 6.23;
   vm->vs[0] = 3.60;
   vm->den[0] = 2.80;

   vm->dep[0] = vm->th[0];
   vm->mu[0] = vm->vs[0]*vm->vs[0]*vm->den[0]*1.0e+10;  /* in CMS units */
   }
else
   {
   fpr = fopfile(vfile,"r");

   fgets(str,512,fpr);

   if(str[0] == '#')
      {
      while(str[0] == '#')
         fgets(str,512,fpr);

      vm->nlay = nblock;
      vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
      vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

      i = 0;
      sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

      vm->dep[i] = vm->th[i];
      vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */

      while(fgets(str,512,fpr) != NULL)
         {
         i++;

         if(i == vm->nlay)
            {
	    vm->nlay = vm->nlay + nblock;
	    vm->vp = (float *)check_realloc(vm->vp,vm->nlay*sizeof(float));
	    vm->vs = (double *)check_realloc(vm->vs,vm->nlay*sizeof(double));
	    vm->den = (float *)check_realloc(vm->den,vm->nlay*sizeof(float));
	    vm->th = (float *)check_realloc(vm->th,vm->nlay*sizeof(float));
	    vm->dep = (float *)check_realloc(vm->dep,vm->nlay*sizeof(float));
	    vm->mu = (float *)check_realloc(vm->mu,vm->nlay*sizeof(float));
	    vm->invb2 = (double *)check_realloc(vm->invb2,vm->nlay*sizeof(double));
	    }

         sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

         vm->dep[i] = vm->dep[i-1] + vm->th[i];
         vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
         }

      vm->nlay = i+1;
      vm->vp = (float *)check_realloc(vm->vp,vm->nlay*sizeof(float));
      vm->vs = (double *)check_realloc(vm->vs,vm->nlay*sizeof(double));
      vm->den = (float *)check_realloc(vm->den,vm->nlay*sizeof(float));
      vm->th = (float *)check_realloc(vm->th,vm->nlay*sizeof(float));
      vm->dep = (float *)check_realloc(vm->dep,vm->nlay*sizeof(float));
      vm->mu = (float *)check_realloc(vm->mu,vm->nlay*sizeof(float));
      vm->invb2 = (double *)check_realloc(vm->invb2,vm->nlay*sizeof(double));
      }
   else
      {
      sscanf(str,"%d",&vm->nlay);

      vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
      vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
      vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

      for(i=0;i<vm->nlay;i++)
         {
         fgets(str,512,fpr);
         sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

         if(i==0)
            vm->dep[i] = vm->th[i];
         else
            vm->dep[i] = vm->dep[i-1] + vm->th[i];

         vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
         }
      }

   fclose(fpr);
   }
}

void get_moment(struct standrupformat *srf,struct velmodel *vm,int use_srf_mu)
{
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float sigma, aslip, mslip, tarea, tmom, slip, slip_tol, sw_tarea;
float amu, ssd, c1, c2, c3, len, wid;
float amax, amin, xmu;
float surf_sigma, surf_aslip, surf_mslip;
int np, surf_np;
int j, k;

float conv = 1.0e-11;
float pi = 3.141592654;

amax = -1.0e+15;
amin = 1.0e+15;
aslip = 0.0;
mslip = -1.0;
surf_aslip = 0.0;
surf_mslip = -1.0;
amu = 0.0;
tarea = 0.0;
tmom = 0.0;
np = srf[0].srf_apnts.np;
surf_np = 0;
prseg_ptr = srf[0].srf_prect.prectseg;
apval_ptr = srf[0].srf_apnts.apntvals;

fprintf(stderr,"np= %d\n",np);

for(k=0;k<np;k++)
   {
   slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	             apval_ptr[k].slip2*apval_ptr[k].slip2 +
	             apval_ptr[k].slip3*apval_ptr[k].slip3);

   if(use_srf_mu != 0 && apval_ptr[k].vs > 0.0 && apval_ptr[k].den > 0.0)
      xmu = apval_ptr[k].vs*apval_ptr[k].vs*apval_ptr[k].den;
   else
      {
      j = 0;
      while(vm->dep[j] < apval_ptr[k].dep)
         j++;

      xmu = vm->mu[j];
      }

   aslip = aslip + slip;
   amu = amu + slip*xmu;
   tarea = tarea + apval_ptr[k].area*1.0e-10;
   tmom = tmom + slip*apval_ptr[k].area*xmu;

   if(slip > mslip)
      mslip = slip;

   if(apval_ptr[k].dep < 1.0 && slip > 1)
      {
      surf_aslip = surf_aslip + slip;
      if(slip > surf_mslip)
         surf_mslip = slip;

      surf_np++;
      }

   if(apval_ptr[k].area < amin)
      amin = apval_ptr[k].area;
   if(apval_ptr[k].area > amax)
      amax = apval_ptr[k].area;
   }

amu = amu/aslip;
aslip = aslip/(float)(np);
slip_tol = 0.05*aslip;
sw_tarea = 0.0;
sigma = 0.0;
for(k=0;k<np;k++)
   {
   slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	             apval_ptr[k].slip2*apval_ptr[k].slip2 +
	             apval_ptr[k].slip3*apval_ptr[k].slip3);

   sigma = sigma + (aslip - slip)*(aslip - slip);

   if(slip > slip_tol)
      sw_tarea = sw_tarea + apval_ptr[k].area*1.0e-10;
   }

if(surf_np > 3)
   {
   surf_aslip = surf_aslip/(float)(surf_np);
   surf_sigma = 0.0;
   for(k=0;k<np;k++)
      {
      slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	                   apval_ptr[k].slip2*apval_ptr[k].slip2 +
	                   apval_ptr[k].slip3*apval_ptr[k].slip3);

      if(apval_ptr[k].dep < 1.0 && slip > 1.0)
         surf_sigma = surf_sigma + (surf_aslip - slip)*(surf_aslip - slip);
      }
   }

len = prseg_ptr[0].flen;
wid = prseg_ptr[0].fwid;
/*
ssd1 = conv*2.0*amu*aslip/(pi*wid);

ssd2 = tmom*(7.0/16.0)*exp(1.5*log(1.0e-10*pi/tarea))*1.0e-06;
*/

ssd = tmom*exp(1.5*log(1.0e-10/tarea))*1.0e-06;
c1 = (16.0/7.0)/exp(1.5*log(pi));
c2 = (pi/2.0)*sqrt(wid/len);
c3 = (3.0*pi/4.0)*sqrt(wid/len);

fprintf(stderr,"Total area= %13.5e km*km\n",tarea);
fprintf(stderr,"Slip-weighted rupture area= %13.5e km*km\n",sw_tarea);
fprintf(stderr,"Average area= %13.5e km*km\n",tarea/(float)(np));
fprintf(stderr,"Min. area= %13.5e km*km\n",amin*1.0e-10);
fprintf(stderr,"Max. area= %13.5e km*km\n",amax*1.0e-10);
fprintf(stderr,"Average slip= %10.2f cm (sigma/aslip= %.3f)\n",aslip,sqrt(sigma/(np-1))/aslip);
fprintf(stderr,"Maximum slip= %10.2f cm\n",mslip);
if(surf_np > 3)
   {
   fprintf(stderr,"SURFACE: Average slip= %10.2f cm (sigma/aslip= %.3f)\n",surf_aslip,sqrt(surf_sigma/(surf_np-1))/surf_aslip);
   fprintf(stderr,"SURFACE: Maximum slip= %10.2f cm\n",surf_mslip);
   }
fprintf(stderr,"Total moment= %13.5e dyne-cm (Mw= %.3f)\n",tmom,0.66667*(log(tmom)/log(10))-10.7);
fprintf(stderr,"Static stress drop= %.3f bar (dip-slip: length= %.2f width= %.2f km)\n",ssd/c3,len,wid);
fprintf(stderr,"Static stress drop= %.3f bar (strike-slip: length= %.2f width= %.2f km)\n",ssd/c2,len,wid);
fprintf(stderr,"Static stress drop= %.3f bar (circular crack: area= %13.5e km*km)\n",ssd/c1,tarea);
}

int write_lld(char *file,struct standrupformat *srf,int ig,float *dmin,float *dmax,char *type,int dump_slip)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float outval, svmax, sv, *stf1, *stf2, *stf3;
float slipval;
int maxnt, it;
int nx, ny;
int i, j, k, npskip;
int igst, ignd, igp;
char str[512], frmt[32];

float ss, dd;
float rperd = 0.017453293;

igst = 0;
ignd = 1;
if(strcmp(srf[0].type,"PLANE") == 0)
   {
   if(ig < 0)
      {
      igst = 0;
      ignd = srf[0].srf_prect.nseg;
      }
   else
      {
      igst = ig;
      ignd = igst + 1;
      }
   }

if(strcmp(file,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(file,"w");

for(igp=igst;igp<ignd;igp++)
   {
   npskip = 0;
   nx = srf[0].srf_apnts.np;
   ny = 1;

   if(strcmp(srf[0].type,"PLANE") == 0)
      {
      prseg_ptr = srf[0].srf_prect.prectseg;

      nx = prseg_ptr[igp].nstk;
      ny = prseg_ptr[igp].ndip;

      npskip = 0;
      for(i=0;i<igp;i++)
         npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;
      }

   apval_ptr = srf[0].srf_apnts.apntvals + npskip;

   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         k = i + j*nx;

         if(apval_ptr[k].dep >= *dmin && apval_ptr[k].dep <= *dmax)
	    {
            if(strcmp(type,"slip") == 0 || dump_slip)
               {
	       outval = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	                  apval_ptr[k].slip2*apval_ptr[k].slip2 +
			  apval_ptr[k].slip3*apval_ptr[k].slip3);

	       slipval = outval;
	       }

            else if(strcmp(type,"velocity") == 0)
               {
               stf1 = apval_ptr[k].stf1;
               stf2 = apval_ptr[k].stf2;
               stf3 = apval_ptr[k].stf3;

               maxnt = apval_ptr[k].nt1;
               if(apval_ptr[k].nt2 > maxnt)
                  maxnt = apval_ptr[k].nt2;
               if(apval_ptr[k].nt3 > maxnt)
                  maxnt = apval_ptr[k].nt3;

               svmax = 0.0;
               for(it=0;it<maxnt;it++)
                  {
                  sv = 0.0;
                  if(it < apval_ptr[k].nt1)
                     sv = sv + stf1[it]*stf1[it];
                  if(it < apval_ptr[k].nt2)
                     sv = sv + stf2[it]*stf2[it];
                  if(it < apval_ptr[k].nt3)
                     sv = sv + stf3[it]*stf3[it];

                  sv = sqrt(sv);
                  if(sv > svmax)
                     svmax = sv;
                  }

               outval = svmax;
               }

            if(strcmp(type,"tinit") == 0)
               outval = apval_ptr[k].tinit;

            if(strcmp(type,"trise") == 0)
               {
               outval = apval_ptr[k].nt1*apval_ptr[k].dt;
               if(apval_ptr[k].nt2*apval_ptr[k].dt > outval)
                  outval = apval_ptr[k].nt2*apval_ptr[k].dt;
               if(apval_ptr[k].nt3*apval_ptr[k].dt > outval)
                  outval = apval_ptr[k].nt3*apval_ptr[k].dt;
               }

            else if(strcmp(type,"strike") == 0)
               outval = apval_ptr[k].stk;

            else if(strcmp(type,"dip") == 0)
               outval = apval_ptr[k].dip;

            else if(strcmp(type,"rake") == 0)
	       {
	       ss = cos(rperd*apval_ptr[k].rake)*apval_ptr[k].slip1 - sin(rperd*apval_ptr[k].rake)*apval_ptr[k].slip2;
	       dd = sin(rperd*apval_ptr[k].rake)*apval_ptr[k].slip1 + cos(rperd*apval_ptr[k].rake)*apval_ptr[k].slip2;

	       outval = 90.0;
	       if(ss != 0.0)
	          outval = atan(dd/ss)/rperd;
               if(ss < 0.0)
	          outval = outval + 180.0;

	       if(ss == 0.0 && dd == 0.0)
                  outval = apval_ptr[k].rake;

	       while(outval >= 360.0)
	          outval = outval - 360.0;
	       while(outval < 0.0)
	          outval = outval + 360.0;
	       }

            if(dump_slip)
               fprintf(fpw,"%12.6f %12.6f %12.5e %12.5e %12.5e\n",apval_ptr[k].lon,apval_ptr[k].lat,apval_ptr[k].dep,outval,slipval);
            else
               fprintf(fpw,"%12.6f %12.6f %12.5e %12.5e\n",apval_ptr[k].lon,apval_ptr[k].lat,apval_ptr[k].dep,outval);
	    }
         }
      }
   }

fclose(fpw);
}

void get_moment_rate(struct standrupformat *srf,struct velmodel *vm,int use_srf_mu)
{
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float amu, *stf;
float dt, *mr, *xmr;
int nt, nalloc, nalloc2, np, its, maxnt, it, j, k;

int nblock = 10000;

np = srf[0].srf_apnts.np;
prseg_ptr = srf[0].srf_prect.prectseg;
apval_ptr = srf[0].srf_apnts.apntvals;

mr = NULL;
nalloc = nblock;
mr = (float *) check_realloc(mr,nalloc*sizeof(float));

for(it=0;it<nalloc;it++)
   mr[it] = 0.0;

xmr = NULL;
nalloc2 = nblock;
xmr = (float *) check_realloc(xmr,nalloc2*sizeof(float));

nt = 0;
dt = apval_ptr[0].dt;
for(k=0;k<np;k++)
   {
   if(use_srf_mu != 0 && apval_ptr[k].vs > 0.0 && apval_ptr[k].den > 0.0)
      amu = apval_ptr[k].area*apval_ptr[k].vs*apval_ptr[k].vs*apval_ptr[k].den;
   else
      {
      j = 0;
      while(vm->dep[j] < apval_ptr[k].dep)
         j++;

      amu = apval_ptr[k].area*vm->mu[j];
      }

   maxnt = apval_ptr[k].nt1;
   if(apval_ptr[k].nt2 > maxnt)
      maxnt = apval_ptr[k].nt2;
   if(apval_ptr[k].nt3 > maxnt)
      maxnt = apval_ptr[k].nt3;

   its = (int)(apval_ptr[k].tinit/apval_ptr[k].dt + 0.5);

   if((its + maxnt + 1) > nt)
      nt = its + maxnt + 1;

   if(nt > nalloc)
      {
      nalloc = nalloc + nblock;
      mr = (float *) check_realloc(mr,nalloc*sizeof(float));

      for(it=nalloc-1;it>=nalloc-nblock;it--)
         mr[it] = 0.0;
      }

   if(maxnt > nalloc2)
      {
      nalloc2 = nalloc2 + nblock;
      xmr = (float *) check_realloc(xmr,nalloc2*sizeof(float));
      }

   for(it=0;it<maxnt;it++)
      xmr[it] = 0.0;

   stf = apval_ptr[k].stf1;
   for(it=0;it<apval_ptr[k].nt1;it++)
      xmr[it] = xmr[it] + stf[it]*stf[it];

   stf = apval_ptr[k].stf2;
   for(it=0;it<apval_ptr[k].nt2;it++)
      xmr[it] = xmr[it] + stf[it]*stf[it];

   stf = apval_ptr[k].stf3;
   for(it=0;it<apval_ptr[k].nt3;it++)
      xmr[it] = xmr[it] + stf[it]*stf[it];

   for(it=0;it<maxnt;it++)
      mr[it+its] = mr[it+its] + amu*sqrt(xmr[it]);
   }

fprintf(stdout,"mr mr\n");
fprintf(stdout,"%d %13.5e\n",nt,dt);

for(it=0;it<nt;it++)
   fprintf(stdout,"%13.5e\n",mr[it]);

free(mr);
free(xmr);
}

#define PI 3.141592653589740

void tfilter(y,dt,nt,order,fhi,flo,phase,space)
float *y, *dt, *fhi, *flo, *space;
int nt, order, phase;
{
struct complex *q, *p, *tmpptr;
double trig_arg;
int j, it, i;
float wplo, wphi, cosA, sinA;
float fnyq;
float are, aim;
float one = 1.0;
float two = 2.0;

p = (struct complex *) space;
q = (struct complex *) (space + 2*nt);
cxzero(p,nt);
cxzero(q,nt);

for(it=0;it<nt;it++)
   p[it].re = y[it];

fnyq = one/(two*(*dt));

if((*fhi))              
   {
   if((*fhi) < fnyq)
      wphi = tan(PI*(*fhi)*(*dt));
   else
      wphi = 1.0e+10;
   }   

if((*flo) < fnyq)
   wplo = tan(PI*(*flo)*(*dt));

           /*  forward pass  */ 

for(j=0;j<order;j++)            
   {
   trig_arg = 0.5*(two*j + one)*PI/order;
   cosA = cos(trig_arg);
   sinA = -sin(trig_arg);

   if((*flo) < fnyq)    /* low-pass filter */
      {
      are = wplo*sinA;
      aim = -wplo*cosA;
      lp_filter(q,p,nt,&are,&aim,1);

      tmpptr = p; p = q; q = tmpptr; 
      }

   if((*fhi) > 0.0)    /* high-pass filter */
      {
      are = wphi*sinA;
      aim = -wphi*cosA;
      hp_filter(q,p,nt,&are,&aim,1);

      tmpptr = p; p = q; q = tmpptr; 
      }
   }

              /*  reverse pass to obtain zero-phase response  */

if(phase == 0)                                                   
   {
   for(j=0;j<order;j++)
      {
      trig_arg = 0.5*(two*j + one)*PI/order;
      cosA = cos(trig_arg);
      sinA = -sin(trig_arg);

      if((*flo) < fnyq)    /* low-pass filter */
         {
         are = wplo*sinA;
         aim = -wplo*cosA;
         lp_filter(q,p,nt,&are,&aim,-1);

         tmpptr = p; p = q; q = tmpptr;  
         }

      if((*fhi) > 0.0)    /* high-pass filter */
         {
         are = wphi*sinA;
         aim = -wphi*cosA;
         hp_filter(q,p,nt,&are,&aim,-1);

         tmpptr = p; p = q; q = tmpptr;  
         }
      }
   }

for(it=0;it<nt;it++)
   y[it] = p[it].re;
}

void cxzero(p,n)
struct complex *p;
int n;
{
float zap = 0.0;
int i;

for(i=0;i<n;i++)
   p[i].re = p[i].im = zap;
}

void hp_filter(q,p,n,alpha,beta,sgn)
struct complex *q, *p;
float *alpha, *beta;
int n, sgn;
{
float are, aim, bre, bim;
float tmpre, tmpim, denom;
float one = 1.0;
int i, n1, k;

tmpre = one - (*alpha);
tmpim = -(*beta);
denom = one/(tmpre*tmpre + tmpim*tmpim);
are = tmpre*denom;
aim = -tmpim*denom;

bre = one + (*alpha);
bim = (*beta);

n1 = n - 1;    
if(sgn == 1)
   i = 1;
if(sgn == -1)
   i = n1 - 1;

k = i - sgn;   

q[k].re = (are*p[k].re - aim*p[k].im);
q[k].im = (are*p[k].im + aim*p[k].re);
 
while(n1--)
   {
   tmpre = (p[i].re - p[k].re) + bre*q[k].re - bim*q[k].im;
   tmpim = (p[i].im - p[k].im) + bre*q[k].im + bim*q[k].re;

   q[i].re = are*tmpre - aim*tmpim;                         
   q[i].im = are*tmpim + aim*tmpre;

   i = i + sgn;                     
   k = i - sgn;
   }
}

void lp_filter(q,p,n,alpha,beta,sgn)
struct complex *q, *p;
float *alpha, *beta;
int n, sgn;
{
float are, aim, bre, bim, cre, cim;
float tmpre, tmpim, denom;
float one = 1.0;
int i, n1, k;

tmpre = one - (*alpha);
tmpim = -(*beta);
denom = one/(tmpre*tmpre + tmpim*tmpim);
are = tmpre*denom;
aim = -tmpim*denom;

bre = one + (*alpha);
bim = (*beta);

cre = -(*alpha);
cim = -(*beta);

n1 = n - 1;     
if(sgn == 1)
   i = 1;
if(sgn == -1)
   i = n1 - 1;

k = i - sgn;   

tmpre = cre*p[k].re - cim*p[k].im;
tmpim = cre*p[k].im + cim*p[k].re;

q[k].re = are*tmpre - aim*tmpim;   
q[k].im = are*tmpim + aim*tmpre;

while(n1--)                      
   {
   tmpre = cre*(p[i].re + p[k].re) - cim*(p[i].im + p[k].im)
      + bre*q[k].re - bim*q[k].im;
   tmpim = cre*(p[i].im + p[k].im) + cim*(p[i].re + p[k].re)
      + bre*q[k].im + bim*q[k].re;

   q[i].re = are*tmpre - aim*tmpim;
   q[i].im = are*tmpim + aim*tmpre;

   i = i + sgn;                     
   k = i - sgn;
   }
}

void get_filtvals(char *cbuf,float *flo,float *fhi,int *ord,int *phs)
{
char temp1[1024];
int i, j;

while(cbuf[0] != '\0')
   {
   if(cbuf[0] == '-')
      {
      cbuf++;

      if(strncmp(cbuf,"flo",3) == 0)
         {
	 cbuf = cbuf + 3;
	 j = 0;

         while(cbuf[j] != '-' && cbuf[j] != '\0')
	    j++;

         strncpy(temp1,cbuf,j);
         temp1[j] = '\0';
         *flo = atof(temp1);

         cbuf = cbuf + j;
	 }
      else if(strncmp(cbuf,"fhi",3) == 0)
         {
	 cbuf = cbuf + 3;
	 j = 0;

         while(cbuf[j] != '-' && cbuf[j] != '\0')
	    j++;

         strncpy(temp1,cbuf,j);
         temp1[j] = '\0';
         *fhi = atof(temp1);

         cbuf = cbuf + j;
	 }
      else if(strncmp(cbuf,"ord",3) == 0)
         {
	 cbuf = cbuf + 3;
	 j = 0;

         while(cbuf[j] != '-' && cbuf[j] != '\0')
	    j++;

         strncpy(temp1,cbuf,j);
         temp1[j] = '\0';
         *ord = atoi(temp1);

         cbuf = cbuf + j;
	 }
      else if(strncmp(cbuf,"phs",3) == 0)
         {
	 cbuf = cbuf + 3;
	 j = 0;

         while(cbuf[j] != '-' && cbuf[j] != '\0')
	    j++;

         strncpy(temp1,cbuf,j);
         temp1[j] = '\0';
         *phs = atoi(temp1);

         cbuf = cbuf + j;
	 }
      else
         cbuf[0] = '\0';
      }
   else
      cbuf[0] = '\0';
   }
}

void dump_sliprate(char *fileroot,struct standrupformat *srf,int ig,float *dmin,float *dmax,float *smin)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float outval, svmax, *sv, *stf1, *stf2, *stf3;
float slipval;
int maxnt, it, jt;
int nx, ny;
int i, j, k, npskip;
int igst, ignd, igp;
char outfile[4096], frmt[32];

int nt_default, nt6;

nt_default = 10000;
sv = (float *)check_malloc(nt_default*sizeof(float));

igst = 0;
ignd = 1;
if(strcmp(srf[0].type,"PLANE") == 0)
   {
   if(ig < 0)
      {
      igst = 0;
      ignd = srf[0].srf_prect.nseg;
      }
   else
      {
      igst = ig;
      ignd = igst + 1;
      }
   }

for(igp=igst;igp<ignd;igp++)
   {
   npskip = 0;
   nx = srf[0].srf_apnts.np;
   ny = 1;

   if(strcmp(srf[0].type,"PLANE") == 0)
      {
      prseg_ptr = srf[0].srf_prect.prectseg;

      nx = prseg_ptr[igp].nstk;
      ny = prseg_ptr[igp].ndip;

      npskip = 0;
      for(i=0;i<igp;i++)
         npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;
      }

   apval_ptr = srf[0].srf_apnts.apntvals + npskip;

   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         k = i + j*nx;

	 slipval = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	                  apval_ptr[k].slip2*apval_ptr[k].slip2 +
			  apval_ptr[k].slip3*apval_ptr[k].slip3);

         if(apval_ptr[k].dep >= *dmin && apval_ptr[k].dep <= *dmax && slipval >= *smin)
	    {
            stf1 = apval_ptr[k].stf1;
            stf2 = apval_ptr[k].stf2;
            stf3 = apval_ptr[k].stf3;

            maxnt = apval_ptr[k].nt1;
            if(apval_ptr[k].nt2 > maxnt)
               maxnt = apval_ptr[k].nt2;
            if(apval_ptr[k].nt3 > maxnt)
               maxnt = apval_ptr[k].nt3;

            if(maxnt > nt_default)
               {
               nt_default = maxnt;
               sv = (float *)check_realloc(sv,nt_default*sizeof(float));
               }

            for(it=0;it<maxnt;it++)
               {
               sv[it] = 0.0;
               if(it < apval_ptr[k].nt1)
                  sv[it] = sv[it] + stf1[it]*stf1[it];
               if(it < apval_ptr[k].nt2)
                  sv[it] = sv[it] + stf2[it]*stf2[it];
               if(it < apval_ptr[k].nt3)
                  sv[it] = sv[it] + stf3[it]*stf3[it];

               sv[it] = sqrt(sv[it]);
               }

            sprintf(outfile,"%s-%.6d.stf",fileroot,npskip+k);

            fpw = fopfile(outfile,"w");

            fprintf(fpw,"%-10s %3s %s\n","stf","stf","stitle");
            fprintf(fpw,"%d %12.5e %d %d %12.5e %12.5e %12.5e %12.5e\n",maxnt,
                                           apval_ptr[k].dt,
                                           0,
                                           0,
                                           0.0,
                                           1.0,
                                           0.0,
                                           0.0);

            nt6 = maxnt/6;
            for(it=0;it<nt6;it++)
               {
               for(jt=0;jt<6;jt++)
                  fprintf(fpw,"%13.5e",sv[6*it + jt]);

               fprintf(fpw,"\n");
               }

            if(6*nt6 != maxnt)
               {
               for(it=6*nt6;it<maxnt;it++)
                  fprintf(fpw,"%13.5e",sv[it]);

               fprintf(fpw,"\n");
               }
            fclose(fpw);
	    }
         }
      }
   }
free(sv);
}
