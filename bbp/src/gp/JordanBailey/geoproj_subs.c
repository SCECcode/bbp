#include "include.h"
#include "structure.h"
#include "function.h"

#define         RPERD           0.017453292
#define         ERAD            6378.139
#define         FLAT_CONST      298.256

void gcproj(float *xf,float *yf,float *rlon,float *rlat,float *ref_rad,double *g0,double *b0,double *amat,double *ainv,int gflag)
{
double xp, yp, zp;
double xg, yg, zg;
double arg;
double cosG, sinG;
double cosB, sinB;

double rperd = RPERD;

if(gflag == 0)
   {
   arg = (*xf)/(*ref_rad) - (*b0);
   cosB = cos(arg);
   sinB = sin(arg);

   arg = (*yf)/(*ref_rad) - (*g0);
   cosG = cos(arg);
   sinG = sin(arg);

   arg = sqrt(1.0 + sinB*sinB*sinG*sinG);
   xp = sinG*cosB*arg;
   yp = sinB*cosG*arg;
   zp = sqrt(1.0 - xp*xp - yp*yp);

   xg = xp*amat[0] + yp*amat[1] + zp*amat[2];
   yg = xp*amat[3] + yp*amat[4] + zp*amat[5];
   zg = xp*amat[6] + yp*amat[7] + zp*amat[8];

/*
    RWG 04/10/2012
    Following is incorrect for southern hemisphere.  Replace with the following
    additional conditionals on zg

    arg = sqrt(xg*xg + yg*yg)/zg;
    (*rlat) = 90.0 - atan(arg)/rperd;
*/

   if(zg != (double)(0.0))
      {
      arg = sqrt(xg*xg + yg*yg)/zg;
      (*rlat) = 90.0 - atan(arg)/rperd;
      if(zg < (double)(0.0))
         (*rlat) = (*rlat) - 180.0;
      }
   else
      (*rlat) = 0.0;

   if(xg != (double)(0.0))
      {
      arg = yg/xg;
      (*rlon) = atan(arg)/rperd;
      }
   else
      (*rlon) = 0.0;

   /*
   RWG 05/08/08
   This is incorrect.  Replaced with following conditional on 'xg'.
   if(yg < (double)(0.0))
   */
   if(xg < (double)(0.0))
      (*rlon) = (*rlon) - 180.0;

   while((*rlon) < (double)(-180.0))
      (*rlon) = (*rlon) + 360.0;
   }
else
   {
   arg = (*rlon)*rperd;
   cosG = cos(arg);
   sinG = sin(arg);

   arg = (90.0 - (*rlat))*rperd;
   cosB = cos(arg);
   sinB = sin(arg);

   xg = sinB*cosG;
   yg = sinB*sinG;
   zg = cosB;

   xp = xg*ainv[0] + yg*ainv[1] + zg*ainv[2];
   yp = xg*ainv[3] + yg*ainv[4] + zg*ainv[5];
   zp = xg*ainv[6] + yg*ainv[7] + zg*ainv[8];

   sinG = xp/sqrt(1.0 - yp*yp);
   sinB = yp/sqrt(1.0 - xp*xp);

   *xf = (double)(*ref_rad)*(asin(sinB)+(*b0));
   *yf = (double)(*ref_rad)*(asin(sinG)+(*g0));
   }
}

void gen_matrices(double *amat,double *ainv,float *alpha,float *ref_lon,float *ref_lat)
{
double arg;
double cosA, sinA;
double cosT, sinT;
double cosP, sinP;
double det;

double rperd = RPERD;

arg = (double)(*alpha)*rperd;
cosA = cos(arg);
sinA = sin(arg);

arg = (double)(90.0-*ref_lat)*rperd;
cosT = cos(arg);
sinT = sin(arg);

arg = (double)(*ref_lon)*rperd;
cosP = cos(arg);
sinP = sin(arg);

amat[0] = cosA*cosT*cosP + sinA*sinP;
amat[1] = sinA*cosT*cosP - cosA*sinP;
amat[2] = sinT*cosP;
amat[3] = cosA*cosT*sinP - sinA*cosP;
amat[4] = sinA*cosT*sinP + cosA*cosP;
amat[5] = sinT*sinP;
amat[6] = -cosA*sinT;
amat[7] = -sinA*sinT;
amat[8] = cosT;

det = amat[0]*(amat[4]*amat[8] - amat[7]*amat[5])
    - amat[1]*(amat[3]*amat[8] - amat[6]*amat[5])
    + amat[2]*(amat[3]*amat[7] - amat[6]*amat[4]);

det = 1.0/det;
ainv[0] = det*amat[0];
ainv[1] = det*amat[3];
ainv[2] = det*amat[6];
ainv[3] = det*amat[1];
ainv[4] = det*amat[4];
ainv[5] = det*amat[7];
ainv[6] = det*amat[2];
ainv[7] = det*amat[5];
ainv[8] = det*amat[8];
}

void geocen(float *r,double x)
{
*r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
fprintf(stderr,"%20.10f %20.10f %20.10f\n",*r,x,atan((1.0 - (1.0/FLAT_CONST))*tan(x)));
}

void set_g2(float *g2,float *fc)
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

void latlon2km(float *arg,float *latkm,float *lonkm,float *rc,float *g2)
{
float cosA, sinA, g2s2, den;

float fone = 1.0;
float ftwo = 2.0;

cosA = cos((*arg));
sinA = sin((*arg));
g2s2 = (*g2)*sinA*sinA;

den = sqrt((fone)/((fone) + g2s2));
*lonkm = (*rc)*cosA*den;
*latkm = (*rc)*(sqrt((fone) + g2s2*((ftwo) + (*g2))))*den*den*den;
}

void set_geoproj(struct sgtmaster *sgtmast,struct geoprojection *geop)
{
float latavg, xlen, ylen;

geop->geoproj = sgtmast->geoproj;
geop->modellon = sgtmast->modellon;
geop->modellat = sgtmast->modellat;
geop->modelrot = sgtmast->modelrot;
geop->xshift = sgtmast->xshift;
geop->yshift = sgtmast->yshift;

geop->rperd = RPERD;
geop->erad = ERAD;

xlen = -2.0*geop->xshift;
ylen = -2.0*geop->yshift;

if(geop->geoproj == 0)
   {
   geop->center_origin = 0;  /* assume for now */
   geop->cosR = cos((geop->modelrot)*(geop->rperd));
   geop->sinR = sin((geop->modelrot)*(geop->rperd));

   geop->fc = FLAT_CONST;
   geop->radc = (geop->erad)*(geop->rperd);
   set_g2(&geop->g2,&geop->fc);

   latavg = geop->modellat;
   if(geop->center_origin == 0)  /* backward compatible */
      latavg = geop->modellat - 0.5*(xlen*geop->sinR + ylen*geop->cosR)/111.20;

   geocen(&latavg,(double)(latavg)*(geop->rperd));
   latlon2km(&latavg,&geop->kmlat,&geop->kmlon,&geop->radc,&geop->g2);
   }
else if(geop->geoproj == 1)
   {
   gen_matrices(geop->amat,geop->ainv,&geop->modelrot,&geop->modellon,&geop->modellat);

   geop->g0 = (double)(0.5*ylen)/(double)(geop->erad);
   geop->b0 = (double)(0.5*xlen)/(double)(geop->erad);
   }
}
