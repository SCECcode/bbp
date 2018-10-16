#include "include.h"

#define		MAX_STAT 10000
#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)

void *check_malloc(int);
FILE *fopfile(char*, char*);

#include "function.h"
#include "getpar.h"

void geoutm_(double *dlon, double *dlat, double *xr0,
	     double *yr0, int *utm_zone, int *geo2utm);

int main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float mlat, mlon, xazim, mrot;
float kperd_n, kperd_e, ar1, ar2, ar3, ar4;
float xp, yp, plon, plat;
float cosR, sinR;
int nr;

float xlen = 100.0;
float ylen = 100.0;

char infile[512], outfile[512], str[512];
char name[16];

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
int geoproj = 1;  /* default is great-circle way */
int center_origin = 0;

float xshift = -1.0e+15;
float yshift = -1.0e+15;

sprintf(outfile,"stdout");
sprintf(infile,"stdin");

setpar(ac, av);
mstpar("mlat","f",&mlat);
mstpar("mlon","f",&mlon);
mstpar("xazim","f",&xazim);

getpar("xlen","f",&xlen);
getpar("ylen","f",&ylen);
getpar("infile","s",infile);
getpar("outfile","s",outfile);

getpar("geoproj","d",&geoproj);
getpar("center_origin","d",&center_origin);
getpar("xshift","f",&xshift);
getpar("yshift","f",&yshift);

endpar();

mrot = xazim - 90.0;

cosR = cos(mrot*rperd);
sinR = sin(mrot*rperd);

if(center_origin != 0)
   {
   xshift = -0.5*xlen;
   yshift = -0.5*ylen;
   }
else
   {
   if(xshift < -1.0e+14)
      xshift = 0.0;
   if(yshift < -1.0e+14)
      yshift = 0.0;
   }

if(geoproj == 0)
   {
   radc = ERAD*RPERD;
   set_g2(&g2,&fc);

   latavg = mlat;
   if(center_origin == 0)  /* backward compatible */
      latavg = mlat - 0.5*(xlen*sinR + ylen*cosR)/111.20;

   geocen(&latavg,(double)(latavg*rperd));
   latlon2km(&latavg,&kperd_n,&kperd_e,&radc,&g2);

   ar1 = cosR/kperd_e;
   ar2 = sinR/kperd_e;
   ar3 = cosR/kperd_n;
   ar4 = sinR/kperd_n;
   }
else if(geoproj == 1)
   {
   gen_matrices(amat,ainv,&mrot,&mlon,&mlat);

   g0 = (double)(0.5*ylen)/(double)(erad);
   b0 = (double)(0.5*xlen)/(double)(erad);
   g0 = 0.0;
   b0 = 0.0;
   }
else if(geoproj == 2)
   {
   dlon = mlon;
   dlat = mlat;
   geoutm_(&dlon,&dlat,&xr0,&yr0,&utm_zone,&geo2utm);

   fprintf(stderr,"UTM Zone= %d xr0= %13.5e yr0= %13.5e\n",utm_zone,xr0,yr0);
   }

name[0] = '\0';

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

while(fgets(str,512,fpr) != NULL)
   {
   nr = sscanf(str,"%f %f %s",&xp,&yp,name);

      if(geoproj == 0)
         {
         plon = mlon + (xp+xshift)*ar1 - (yp+yshift)*ar2;
         plat = mlat - (xp+xshift)*ar4 - (yp+yshift)*ar3;
         }
      else if(geoproj == 1)
         {
         gcproj(&xp,&yp,&plon,&plat,&erad,&g0,&b0,amat,ainv,xy2ll);
         }
      else if(geoproj == 2)
         {
         dxr = xr0 + 1000.0*((xp+xshift)*cosR - (yp+yshift)*sinR);
         dyr = yr0 - 1000.0*((xp+xshift)*sinR + (yp+yshift)*cosR);
         geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);

	 /*
         dxr = xp;
         dyr = yp;
         geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);
	 */

         plon = dlon;
         plat = dlat;
         }

   if(nr == 3)
      fprintf(fpw,"%12.6f %12.6f %s\n",plon,plat,name);

   if(nr == 2)
      fprintf(fpw,"%12.6f %12.6f\n",plon,plat);
   }

fclose(fpr);
fclose(fpw);
}

double geocenX(x)
double x;
{
double r;
r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
return(r);
}

void set_g2X(g2,fc)
float *g2, *fc;
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

void latlon2kmX(arg,latkm,lonkm,rc,g2)
float *arg, *latkm, *lonkm, *rc, *g2;
{
float cosA, sinA, g2s2, den;

cosA = cos((*arg));
sinA = sin((*arg));
g2s2 = (*g2)*sinA*sinA;

den = sqrt((FONE)/((FONE) + g2s2));
*lonkm = (*rc)*cosA*den;
*latkm = (*rc)*(sqrt((FONE) + g2s2*((FTWO) + (*g2))))*den*den*den;
}

FILE *fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

void *check_malloc(int len)
{
char *ptr;

ptr = (char *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }

return(ptr);
}
