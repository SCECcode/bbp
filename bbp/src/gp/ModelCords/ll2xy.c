#include "include.h"

#define		MAXN 100000
#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)

void *check_malloc(size_t);
void *check_realloc(void *ptr,size_t len);
FILE *fopfile(char*, char*);

#include "function.h"
#include "getpar.h"

void geoutm_(double *dlon, double *dlat, double *xr0,
	     double *yr0, int *utm_zone, int *geo2utm);

int main(int ac,char **av)
{
FILE *fopfile(), *fpw, *fpr;
float mlat, mlon, xazim;
float kperd_n, kperd_e, xs, ys, plon, plat;
int i, mn, ns, nx, ny, test;
float cosR, sinR, xr, yr;

float xlen = 100.0;
float ylen = 100.0;

char infile[512], outfile[512], str[512], name[64];

float mrot = 0.0;
float h = 1.0;

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

float zscale = 1.0;
float zupper = 1.0e+15;
float zlower = -1.0e+15;

sprintf(outfile,"stdout");
sprintf(infile,"stdin");

setpar(ac, av);

mstpar("mlat","f",&mlat);
mstpar("mlon","f",&mlon);
mstpar("xazim","f",&xazim);

getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("xlen","f",&xlen);
getpar("ylen","f",&ylen);
getpar("zscale","f",&zscale);
getpar("zupper","f",&zupper);
getpar("zlower","f",&zlower);

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

   fprintf(stderr,"ke=%12.4f kn=%12.4f latavg=%10.4f\n",kperd_e,kperd_n,latavg/rperd);
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

   fprintf(stderr,"UTM Zone= %d xr0= %15.8e yr0= %15.8e\n",utm_zone,xr0,yr0);
   }

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
   ns = sscanf(str,"%f %f %s",&plon,&plat,name);

   if(geoproj == 0)
      {
      xs = (plon - mlon)*kperd_e;
      ys = (mlat - plat)*kperd_n;

      xr = xs*cosR + ys*sinR - xshift;
      yr = -xs*sinR + ys*cosR - yshift;
      }
   else if(geoproj == 1)
      {
      gcproj(&xr,&yr,&plon,&plat,&erad,&g0,&b0,amat,ainv,ll2xy);
      }
   else if(geoproj == 2)
      {
      dlon = plon;
      dlat = plat;
      geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&geo2utm);

      xs = 0.001*(dxr - xr0);
      ys = 0.001*(yr0 - dyr);

      xr = xs*cosR + ys*sinR - xshift;
      yr = -xs*sinR + ys*cosR - yshift;
      }

/*
   if(ns == 3)
      fprintf(fpw,"%14.6e %14.6e %s\n",xr,yr,name);
   else if(ns == 2)
      fprintf(fpw,"%14.6e %14.6e\n",xr,yr);
*/
   ns = 0;
   while(str[ns] == ' ' || str[ns] == '\t')
      ns++;

   test = 0;
   while(test != 2 && str[ns] != '\n')
      {
      if(str[ns] == ' ' || str[ns] == '\t' || str[ns] == '\n')
         {
         test++;
	 if(test != 2)
	    {
            ns++;
            while(str[ns] == ' ' || str[ns] == '\t')
               ns++;
	    }
	 }
      else
         ns++;
      }

   fprintf(fpw,"%14.6e %14.6e%s",xr,yr,&str[ns]);
   }
fclose(fpr);
fclose(fpw);
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

void *check_malloc(size_t len)
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
