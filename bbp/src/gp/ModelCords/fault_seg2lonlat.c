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

main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float elon, elat, stk, dip, dtop, flen, fwid, l2;
float mlat, mlon, mrot;
float kperd_n, kperd_e, ar1, ar2, ar3, ar4;
float wp, xp, yp, zp, plon, plat;
float cosR, sinR, cosD, sinD;
float ds, dd;
int nstk, ndip, i, j;

char outfile[512];

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

sprintf(outfile,"stdout");

setpar(ac, av);
mstpar("elat","f",&elat);
mstpar("elon","f",&elon);
mstpar("stk","f",&stk);
mstpar("dip","f",&dip);
mstpar("dtop","f",&dtop);
mstpar("flen","f",&flen);
mstpar("fwid","f",&fwid);
mstpar("nstk","d",&nstk);
mstpar("ndip","d",&ndip);

getpar("outfile","s",outfile);
getpar("geoproj","d",&geoproj);

endpar();

mlon = elon;
mlat = elat;
mrot = stk - 90.0;

cosR = cos(mrot*rperd);
sinR = sin(mrot*rperd);

xshift = -0.5*xlen;
yshift = -0.5*ylen;

if(geoproj == 0)
   {
   radc = ERAD*RPERD;
   set_g2(&g2,&fc);

   latavg = mlat;
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
   }
else if(geoproj == 2)
   {
   dlon = mlon;
   dlat = mlat;
   geoutm_(&dlon,&dlat,&xr0,&yr0,&utm_zone,&geo2utm);
   }

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

ds = flen/(float)(nstk);
dd = fwid/(float)(ndip);

l2 = 0.5*flen;

cosD = cos(dip*rperd);
sinD = sin(dip*rperd);

for(j=0;j<ndip;j++)
   {
   wp = (j+0.5)*dd;
   yp = wp*cosD - yshift;
   zp = wp*sinD + dtop;
   for(i=0;i<nstk;i++)
      {
      xp = (i+0.5)*ds - l2 - xshift;

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
         dxr = xr0 + 1000.0*((xp+xshift)*ar1 - (yp+yshift)*ar2);
         dyr = yr0 - 1000.0*((xp+xshift)*ar4 + (yp+yshift)*ar3);
         geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);

         plon = dlon;
         plat = dlat;
         }

      fprintf(fpw,"%11.5f %11.5f %11.5f\n",plon,plat,zp);
      }
   }
fclose(fpw);
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
