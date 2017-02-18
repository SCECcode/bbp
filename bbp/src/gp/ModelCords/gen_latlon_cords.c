#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <sys/fault.h>
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

char *check_malloc();

main(ac,av)
int ac;
char **av;
{
FILE *fpw, *fopfile();
float *flon, *flat, modellon, modellat, modelrot, x, y;
float h, ar1, ar2, ar3, ar4, latavg;
float cosR, sinR, kmlon, kmlat;
float dep, *cc1, *cc2, conv;
int nx, ny, nzout, ix, iy, iz, ip;
char modfile[128], name[128], str[128];

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc;
double geocen();

int latfirst = 1;
int feetout = 1;
int gzip = 1;

float zshift = 0.0;

setpar(ac,av);
mstpar("name","s",name);
mstpar("nx","d",&nx);
mstpar("ny","d",&ny);
mstpar("nzout","d",&nzout);
mstpar("h","f",&h);
mstpar("modellat","f",&modellat);
mstpar("modellon","f",&modellon);
mstpar("modelrot","f",&modelrot);
getpar("latfirst","d",&latfirst);
getpar("feetout","d",&feetout);
getpar("gzip","d",&gzip);
getpar("zshift","f",&zshift);
endpar();

radc = ERAD*RPERD;
set_g2(&g2,&fc);

latavg = geocen(modellat*rperd);
latlon2km(&latavg,&kmlat,&kmlon,&radc,&g2);

cosR = cos(modelrot*RPERD);
sinR = sin(modelrot*RPERD);

ar1 = h*cosR/kmlon;
ar2 = h*sinR/kmlon;
ar3 = h*cosR/kmlat;
ar4 = h*sinR/kmlat;

flon = (float *) check_malloc (nx*ny*sizeof(float));
flat = (float *) check_malloc (nx*ny*sizeof(float));

for(iy=0;iy<ny;iy++)
   {
   fprintf(stderr,"iy=%d\n",iy);

   y = (iy + 0.5);
   for(ix=0;ix<nx;ix++)
      {
      x = (ix + 0.5);
      ip = ix + iy*nx;

      flon[ip] = modellon + x*ar1 - y*ar2;
      flat[ip] = modellat - x*ar4 - y*ar3;
      }
   }

conv = 1.0;
if(feetout)
   {
   conv = KM2FT;
   fprintf(stderr,"**** depth output in FEET\n");
   }
else
   fprintf(stderr,"**** depth output in KM\n");

cc1 = flon;
cc2 = flat;
if(latfirst)
   {
   cc1 = flat;
   cc2 = flon;
   fprintf(stderr,"**** coordinate output format is LAT LON DEP\n");
   }
else
   fprintf(stderr,"**** coordinate output format is LON LAT DEP\n");

if(gzip)
   fprintf(stderr,"**** output files are compressed with 'gzip' (.gz)\n");

for(iz=0;iz<nzout;iz++)
   {
   if(iz < 10)
      sprintf(modfile,"%s00%1d",name,iz);
   else if(iz < 100)
      sprintf(modfile,"%s0%2d",name,iz);
   else
      sprintf(modfile,"%s%3d",name,iz);

   fpw = fopfile(modfile,"w");

   dep = (iz*h + zshift)*conv;

   fprintf(stderr,"iz=%5d dep=%11.4f\n",iz,dep);

   for(iy=0;iy<ny;iy++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 ip = ix + iy*nx;

         fprintf(fpw,"%10.4f %10.4f %10.2f\n",cc1[ip],cc2[ip],dep);
         }
      }
   fclose(fpw);

   if(gzip)
      {
      sprintf(str,"gzip -f %s",modfile);
      system(str);
      }
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

char *check_malloc(len)
int len;
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
