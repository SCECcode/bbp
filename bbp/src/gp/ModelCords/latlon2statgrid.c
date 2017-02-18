#include "include.h"

#define		MAX_STAT 100000
#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)

void *check_malloc(size_t);
FILE *fopfile(char*, char*);

#include "function.h"

struct statinfo
   {
   char name[8];
   float lat;
   float lon;
   float xp;
   float yp;
   float zp;
   int ix;
   int iy;
   int iz;
   };

main(int ac,char **av)
{
FILE *fopfile(), *fp;
struct statinfo si[MAX_STAT];
float mlat, mlon;
float kperd_n, kperd_e, xlen, ylen, xs, ys;
int i, ns, nx, ny, test;
float cosR, sinR, xr, yr;

char infile[512], outfile[512], str[512];

float rotate = 0.0;
float h = 1.0;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;

double g0, b0;
double amat[9], ainv[9];
int xy2ll = 0;
int ll2xy = 1;

int read_depth = 0;
int use_all = 0;
int printall = 1;
int xbnd = 0;
int ybnd = 0;

int var_grid = 0;

double xr0, yr0, dlon, dlat, dxr, dyr;
int geo2utm = 0;
int utm2geo = 1;
int utm_zone = 11;

/*
   geoproj=0: RWG spherical projection with local kmplat, kmplon
          =1: RWG great circle projection
          =2: UTM coordinate projection
*/
int geoproj = 0;  /* default is OLD way */
int center_origin = 0;

float xshift = -1.0e+15;
float yshift = -1.0e+15;

setpar(ac, av);
mstpar("mlat","f",&mlat);
mstpar("mlon","f",&mlon);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("rotate","f",&rotate);
getpar("use_all","d",&use_all);
getpar("read_depth","d",&read_depth);

getpar("var_grid","d",&var_grid);

getpar("geoproj","d",&geoproj);
getpar("center_origin","d",&center_origin);
getpar("xshift","f",&xshift);
getpar("yshift","f",&yshift);

if(var_grid == 1)
   {
   mstpar("xlen","f",&xlen);
   mstpar("ylen","f",&ylen);
   }
else
   {
   mstpar("nx","d",&nx);
   mstpar("ny","d",&ny);
   mstpar("h","f",&h);
   getpar("xbnd","d",&xbnd);
   getpar("ybnd","d",&ybnd);
   getpar("printall","d",&printall);

   xlen = nx*h;
   ylen = ny*h;
   }
endpar();

cosR = cos(rotate*rperd);
sinR = sin(rotate*rperd);

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
   gen_matrices(amat,ainv,&rotate,&mlon,&mlat);

   g0 = (double)(0.5*ylen)/(double)(erad);
   b0 = (double)(0.5*xlen)/(double)(erad);
   }
else if(geoproj == 2)
   {
   dlon = mlon;
   dlat = mlat;
   geoutm_(&dlon,&dlat,&xr0,&yr0,&utm_zone,&geo2utm);

   fprintf(stderr,"UTM Zone= %d\n",utm_zone);
   }

fp = fopfile(infile,"r");

ns = 0;
while(fgets(str,512,fp) != NULL)
   {
   if(str[0] != '#')
      {
      if(read_depth)
         sscanf(str,"%f %f %s %f",&si[ns].lon,&si[ns].lat,si[ns].name,&si[ns].zp);
      else
         {
         sscanf(str,"%f %f %s",&si[ns].lon,&si[ns].lat,si[ns].name);
         si[ns].zp = 0.0;
         }

      test = 1;
      if(use_all == 0)
         {
         for(i=0;i<ns;i++)
            {
            if(strcmp(si[i].name,si[ns].name) == 0)
	       {
	       if(si[i].lon == si[ns].lon && si[i].lat == si[ns].lat)
	          {
	          test = 0;
	          fprintf(stderr,"Station '%s' duplicated, only first occurrence used\n",si[ns].name);
	          }
	       else
	          {
	          fprintf(stderr,"Different location for '%s', check output\n",si[ns].name);
	          sprintf(si[ns].name,"%sX",si[ns].name);
	          }
	       }
            }
         }

      if(test)
         {
         if(geoproj == 0)
	    {
            xs = (si[ns].lon - mlon)*kperd_e;
            ys = (mlat - si[ns].lat)*kperd_n;

            xr = xs*cosR + ys*sinR - xshift;
            yr = -xs*sinR + ys*cosR - yshift;
            }
         else if(geoproj == 1)
            {
            gcproj(&xr,&yr,&si[ns].lon,&si[ns].lat,&erad,&g0,&b0,amat,ainv,ll2xy);
	    }
         else if(geoproj == 2)
            {
	    dlon = si[ns].lon;
	    dlat = si[ns].lat;
	    geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&geo2utm);

	    xs = 0.001*(dxr - xr0);
	    ys = 0.001*(yr0 - dyr);

            xr = xs*cosR + ys*sinR - xshift;
            yr = -xs*sinR + ys*cosR - yshift;
	    }

         si[ns].xp = xr;
         si[ns].yp = yr;

         si[ns].ix = (int)(xr/h + 0.5);
         si[ns].iy = (int)(yr/h + 0.5);
         si[ns].iz = (int)(si[ns].zp/h + 1.5);

         if(printall == 0)
	    {
            if(si[ns].ix >= xbnd && si[ns].ix < nx-xbnd && si[ns].iy >= ybnd && si[ns].iy < ny-ybnd)
               ns++;
	    }
         else
	    ns++;
         }
      }
   }
fclose(fp);

fp = fopfile(outfile,"w");

fprintf(fp,"%d\n",ns);
if(var_grid == 1)
   {
   for(i=0;i<ns;i++)
      fprintf(fp,"%10.5f %10.5f %10.5f %s\n",si[i].xp,si[i].yp,si[i].zp,si[i].name);
   }
else
   {
   for(i=0;i<ns;i++)
      fprintf(fp,"%5d %5d %5d %s\n",si[i].ix,si[i].iy,si[i].iz,si[i].name);
   }

fclose(fp);
}

double geocenX(x)
double x;
{
double r;
r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
return(r);
}

set_g2X(g2,fc)
float *g2, *fc;
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

latlon2kmX(arg,latkm,lonkm,rc,g2)
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
