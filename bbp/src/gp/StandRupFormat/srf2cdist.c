#include "include.h"
#include "structure.h"
#include "function.h"
#include "../../ModelCords/function.h"
#include "defs.h"

#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292

main(int ac,char **av)
{
FILE *fpr, *fpw;
char sbuf[1024];
float slon, slat;
char srffile[256], outfile[256], infile[256];

struct standrupformat srf;
struct srf_apointvalues *apval_ptr;
int k, ip, np;
int inbin = 0;

float *xf, *yf, *z2, *s2, rjb, minjb, rr, minr, mind, avgdip;

float kperd_n, kperd_e, xs, ys, xr, yr;
float rotate, mlon, mlat, cosR, sinR;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;

double g0, b0;
double amat[9], ainv[9];
int xy2ll = 0;
int ll2xy = 1;

double xr0, yr0, dlon, dlat, dxr, dyr;
int geo2utm = 0;
int utm2geo = 1;
int utm_zone = 11;

/*
   geoproj=0: RWG spherical projection with local kmplat, kmplon
          =1: RWG great circle projection
          =2: UTM coordinate projection
*/
int geoproj = 1;

float xshift;
float yshift;

float xlen = 1000.0;
float ylen = 1000.0;

int cdist = 1;
int dip = 1;
int dtop = 1;
int jbdist = 1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
mstpar("srffile","s",srffile);
getpar("outfile","s",outfile);
getpar("infile","s",infile);
getpar("geoproj","d",&geoproj);
getpar("cdist","d",&cdist);
getpar("dip","d",&dip);
getpar("dtop","d",&dtop);
getpar("jbdist","d",&jbdist);
endpar();

read_srf(&srf,srffile,inbin);
np = srf.srf_apnts.np;
apval_ptr = srf.srf_apnts.apntvals;

mlon = 0.0;
mlat = 0.0;
for(ip=0;ip<np;ip++)
   {
   mlon = mlon + apval_ptr[ip].lon;
   mlat = mlat + apval_ptr[ip].lat;
   }

mlon = mlon/(float)(np);
mlat = mlat/(float)(np);
rotate = -90.0;

fprintf(stderr,"Approximate centroid: %12.5f %12.5f\n",mlon,mlat);

cosR = cos(rotate*rperd);
sinR = sin(rotate*rperd);

xshift = -0.5*xlen;
yshift = -0.5*ylen;

if(geoproj == 0)
   {
   radc = ERAD*RPERD;
   set_g2(&g2,&fc);

   latavg = mlat;
   geocen(&latavg,(double)(latavg*rperd));
   latlon2km(&latavg,&kperd_n,&kperd_e,&radc,&g2);
   }
else if(geoproj == 1)
   {
   gen_matrices(amat,ainv,&rotate,&mlon,&mlat);

   g0 = (double)(0.5*ylen)/(double)(erad);
   b0 = (double)(0.5*xlen)/(double)(erad);
   }
else if(geoproj == 2)
   {
   dlon = mlon;
   dlat = mlat;
   geoutm_(&dlon,&dlat,&xr0,&yr0,&utm_zone,&geo2utm);
   }

xf = (float *)check_malloc(np*sizeof(float));
yf = (float *)check_malloc(np*sizeof(float));
z2 = (float *)check_malloc(np*sizeof(float));
s2 = (float *)check_malloc(np*sizeof(float));

avgdip = 0.0;
for(ip=0;ip<np;ip++)
   {
   if(geoproj == 0)
      {
      xs = (apval_ptr[ip].lon - mlon)*kperd_e;
      ys = (mlat - apval_ptr[ip].lat)*kperd_n;

      xf[ip] = xs*cosR + ys*sinR - xshift;
      yf[ip] = -xs*sinR + ys*cosR - yshift;
      }
   else if(geoproj == 1)
      {
      gcproj(&xf[ip],&yf[ip],&apval_ptr[ip].lon,&apval_ptr[ip].lat,&erad,&g0,&b0,amat,ainv,ll2xy);
      }
   else if(geoproj == 2)
      {
      dlon = apval_ptr[ip].lon;
      dlat = apval_ptr[ip].lat;
      geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&geo2utm);

      xs = 0.001*(dxr - xr0);
      ys = 0.001*(yr0 - dyr);

      xf[ip] = xs*cosR + ys*sinR - xshift;
      yf[ip] = -xs*sinR + ys*cosR - yshift;
      }

   s2[ip] = apval_ptr[ip].slip1*apval_ptr[ip].slip1 + apval_ptr[ip].slip2*apval_ptr[ip].slip2 + apval_ptr[ip].slip3*apval_ptr[ip].slip3;
   z2[ip] = apval_ptr[ip].dep*apval_ptr[ip].dep;
   avgdip = avgdip + apval_ptr[ip].dip;
   }
avgdip = avgdip/(float)(np);

fprintf(stderr,"Processing stations ... ");

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

while(fgets(sbuf,1024,fpr) != NULL)
   {
   sscanf(sbuf,"%f %f",&slon,&slat);

   if(geoproj == 0)
      {
      xs = (slon - mlon)*kperd_e;
      ys = (mlat - slat)*kperd_n;

      xr = xs*cosR + ys*sinR - xshift;
      yr = -xs*sinR + ys*cosR - yshift;
      }
   else if(geoproj == 1)
      {
      gcproj(&xr,&yr,&slon,&slat,&erad,&g0,&b0,amat,ainv,ll2xy);
      }
   else if(geoproj == 2)
      {
      dlon = slon;
      dlat = slat;
      geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&geo2utm);

      xs = 0.001*(dxr - xr0);
      ys = 0.001*(yr0 - dyr);

      xr = xs*cosR + ys*sinR - xshift;
      yr = -xs*sinR + ys*cosR - yshift;
      }

   mind = 1.0e+15;
   minr = 1.0e+15;
   minjb = 1.0e+15;
   for(ip=0;ip<np;ip++)
      {
      if(s2[ip] > 0.01)
         {
         rjb = (xr-xf[ip])*(xr-xf[ip]) + (yr-yf[ip])*(yr-yf[ip]);
         rr = rjb + z2[ip];
         if(rjb < minjb)
            minjb = rjb;
         if(rr < minr)
            minr = rr;
         if(z2[ip] < mind)
            mind = z2[ip];
	 }
      }

   k = 0;
   while(sbuf[k] != '\n' && k < 1023)
      k++;

   sbuf[k] = '\0';

   fprintf(fpw,"%s",sbuf);

   if(cdist)
      fprintf(fpw,"\t%.2f",sqrt(minr));
   if(dtop)
      fprintf(fpw,"\t%.2f",sqrt(mind));
   if(dip)
      fprintf(fpw,"\t%.0f",avgdip);
   if(jbdist)
      fprintf(fpw,"\t%.2f",sqrt(minjb));

   fprintf(fpw,"\n");
   }

fprintf(stderr,"DONE\n");

fclose(fpr);
fclose(fpw);
}
