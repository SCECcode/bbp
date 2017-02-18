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
FILE *fpw, *fpr, *fopfile();
float flon, flat, modellon, modellat, modelrot, x, y, z;
float ar1, ar2, ar3, ar4, latavg;
float cosR, sinR, kmlon, kmlat;
int ip, np, n;
char outfile[128], infile[128], str[512], name[64];

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc;
double geocen();

float h = 1.0;

setpar(ac,av);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("h","f",&h);
mstpar("modellat","f",&modellat);
mstpar("modellon","f",&modellon);
mstpar("modelrot","f",&modelrot);
endpar();

radc = ERAD*RPERD;
set_g2(&g2,&fc);

latavg = geocen(modellat*rperd);
latlon2km(&latavg,&kmlat,&kmlon,&radc,&g2);

fprintf(stderr,"lon= %13.5f lat= %13.5f\n",modellon,modellat);
fprintf(stderr,"kmlon= %12.4f kmlat= %12.4f\n",kmlon,kmlat);

cosR = cos(modelrot*RPERD);
sinR = sin(modelrot*RPERD);

ar1 = h*cosR/kmlon;
ar2 = h*sinR/kmlon;
ar3 = h*cosR/kmlat;
ar4 = h*sinR/kmlat;

fpr = fopfile(infile,"r");
fpw = fopfile(outfile,"w");

fgets(str,512,fpr);
sscanf(str,"%d",&np);

for(ip=0;ip<np;ip++)
   {
   fgets(str,512,fpr);
   n = sscanf(str,"%f %f %f %s",&x,&y,&z,&name);

   if(n >= 2)
      {
      if(n < 4)
         name[0] = '\0';

      flon = modellon + x*ar1 - y*ar2;
      flat = modellat - x*ar4 - y*ar3;

      fprintf(fpw,"%10.4f %10.4f %10s\n",flon,flat,name);
      }
   }

fclose(fpr);
fclose(fpw);
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
