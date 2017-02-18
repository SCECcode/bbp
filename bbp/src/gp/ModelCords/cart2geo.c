#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/file.h>
#include <sys/procfs.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>

#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)
#define         KM2FT           (float)(3280.84)
#define MAX_NLEN 100

main(int ac,char **av)
{
FILE *fpw, *fpr, *fopfile();
float rlon, rlat, east, nort, ref_lon, ref_lat;
float latavg, kmlon, kmlat;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc;
double geocen();

double xr0, yr0, dlon, dlat, dxr, dyr;
int geo2utm = 0;
int utm2geo = 1;
int utm_zone = 11;

int utm = 0;

setpar(ac,av);
getpar("utm","d",&utm);
mstpar("ref_lon","f",&ref_lon);
mstpar("ref_lat","f",&ref_lat);
endpar();

if(utm == 1)
   {
   dlon = ref_lon;
   dlat = ref_lat;
   geoutm_(&dlon,&dlat,&xr0,&yr0,&utm_zone,&geo2utm);
   }
else
   {
   radc = ERAD*RPERD;
   set_g2(&g2,&fc);

   latavg = geocen(ref_lat*rperd);
   latlon2km(&latavg,&kmlat,&kmlon,&radc,&g2);
   }

while(scanf("%f %f",&nort,&east) == 2)
   {
   if(utm == 1)
      {
      dxr = xr0 + 1000.0*east;
      dyr = yr0 + 1000.0*nort;
      geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);
      rlon = dlon;
      rlat = dlat;
      }
   else
      {
      rlon = ref_lon + east/kmlon;
      rlat = ref_lat + nort/kmlat;
      }
   printf("%.5f %.5f\n",rlon,rlat);
   }
}

FILE *fopfile(name,mode)
char *name, *mode;
{
FILE *fp, *fopen();

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

double geocen(x)
double x;
{
double r;
r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
return(r);
}

set_g2(g2,fc)
float *g2, *fc;
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

latlon2km(arg,latkm,lonkm,rc,g2)
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
