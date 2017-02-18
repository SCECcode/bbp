#include <sys/file.h>
#include <stdio.h>
#include <math.h>

#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)

void *check_malloc(int);
FILE *fopfile(char*, char*);

main (ac, av)
int ac;
char **av;
{
FILE *fopfile(), *fpw;
float lat, lon, len, wid, strike, dip, kperd_e, kperd_n, latarg;
float flat[4], flon[4];
float cosA, sinA, cosD, sinD, l2, ws, dsn, dse, ddn, dde;
char name[128], fname[512];

float rperd = RPERD;
float fc = FLAT_CONST;
float g2, radc;

setpar(ac, av);
mstpar("lat","f",&lat);
mstpar("lon","f",&lon);
mstpar("len","f",&len);
mstpar("wid","f",&wid);
mstpar("strike","f",&strike);
mstpar("dip","f",&dip);
mstpar("name","s",name);
endpar();

radc = ERAD*RPERD;
set_g2(&g2,&fc);

latarg = lat*rperd;
latlon2km(&latarg,&kperd_n,&kperd_e,&radc,&g2);

cosA = cos(strike*rperd);
sinA = sin(strike*rperd);
cosD = cos((strike+90.0)*rperd);
sinD = sin((strike+90.0)*rperd);

l2 = 0.5*len;
ws = wid*cos(dip*rperd);

dsn = l2*cosA/kperd_n;
dse = l2*sinA/kperd_e;
ddn = ws*cosD/kperd_n;
dde = ws*sinD/kperd_e;

flat[0] = lat + dsn;
flon[0] = lon + dse;

flat[1] = flat[0] + ddn;
flon[1] = flon[0] + dde;

flat[3] = lat - dsn;
flon[3] = lon - dse;

flat[2] = flat[3] + ddn;
flon[2] = flon[3] + dde;

sprintf(fname,"%s.topcords",name);
fpw = fopfile(fname,"w");
fprintf(fpw,"%10.3f %10.3f\n",flon[0],flat[0]);
fprintf(fpw,"%10.3f %10.3f\n",flon[3],flat[3]);
fclose(fpw);

sprintf(fname,"%s.poly",name);
fpw = fopfile(fname,"w");
/*
fprintf(fpw,"5\n");
*/
fprintf(fpw,"%10.3f %10.3f\n",flon[0],flat[0]);
fprintf(fpw,"%10.3f %10.3f\n",flon[1],flat[1]);
fprintf(fpw,"%10.3f %10.3f\n",flon[2],flat[2]);
fprintf(fpw,"%10.3f %10.3f\n",flon[3],flat[3]);
fprintf(fpw,"%10.3f %10.3f\n",flon[0],flat[0]);
fclose(fpw);
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
