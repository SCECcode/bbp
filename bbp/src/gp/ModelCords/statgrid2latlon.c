#include <sys/file.h>
#include <stdio.h>
#include <math.h>

#define		MAX_STAT 10000
#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)

char *check_malloc();

struct statinfo
   {
   char name[8];
   float lat;
   float lon;
   int x;
   int y;
   int z;
   }

main(ac,av)
int ac;
char **av;
{
FILE *fopfile(), *fp;
struct statinfo si[MAX_STAT];
float h, mlat, mlon;
float kperd_n, kperd_e, ar1, ar2, ar3, ar4;
int i, ns;
float cosR, sinR, xr, yr;

char infile[512], outfile[512], str[512];

float rotate = 0.0;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;
double geocen();

setpar(ac, av);
mstpar("mlat","f",&mlat);
mstpar("mlon","f",&mlon);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
mstpar("h","f",&h);
getpar("rotate","f",&rotate);
endpar();

radc = ERAD*RPERD;
set_g2(&g2,&fc);

latavg = geocen(mlat*rperd);
latlon2km(&latavg,&kperd_n,&kperd_e,&radc,&g2);

cosR = cos(rotate*rperd);
sinR = sin(rotate*rperd);

ar1 = h*cosR/kperd_e;
ar2 = h*sinR/kperd_e;
ar3 = h*cosR/kperd_n;
ar4 = h*sinR/kperd_n;

fp = fopfile(infile,"r");
fscanf(fp,"%d",&ns);

for(i=0;i<ns;i++)
   {
   fscanf(fp,"%d %d %d %s",&si[i].x,&si[i].y,&si[i].z,si[i].name);

   xr = (si[i].x + 0.5);
   yr = (si[i].y + 0.5);

   si[i].lon = mlon + xr*ar1 - yr*ar2;
   si[i].lat = mlat - xr*ar4 - yr*ar3;
   }
fclose(fp);

fp = fopfile(outfile,"w");

for(i=0;i<ns;i++)
   fprintf(fp,"%10.4f %10.4f %s\n",si[i].lon,si[i].lat,si[i].name);
fclose(fp);
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
