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

main(ac,av)
int ac;
char **av;
{
FILE *fopfile(), *fpr, *fpw;
float slat, slon;
float elat, elon;
float kperd_n, kperd_e, xs, ys, rs, azi, dperr;

char infile[512], outfile[512], str[512], sname[8];

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;
double geocen();

int use_all = 0;
int printall = 1;
int xbnd = 0;
int ybnd = 0;

setpar(ac, av);
mstpar("elat","f",&elat);
mstpar("elon","f",&elon);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
endpar();

dperr = 1.0/rperd;
radc = ERAD*RPERD;
set_g2(&g2,&fc);

latavg = elat;
latavg = geocen(latavg*rperd);
latlon2km(&latavg,&kperd_n,&kperd_e,&radc,&g2);

fprintf(stderr,"ke=%12.4f kn=%12.4f latavg=%10.4f\n",kperd_e,kperd_n,latavg/rperd);

fpr = fopfile(infile,"r");
fpw = fopfile(outfile,"w");

while(fgets(str,512,fpr) != NULL)
   {
   sscanf(str,"%f %f %s",&slon,&slat,sname);

   xs = (slon - elon)*kperd_e;
   ys = (slat - elat)*kperd_n;

   rs = sqrt(xs*xs + ys*ys);

   azi = dperr*atan2(xs,ys);

   while(azi < 0.0)
      azi = azi + 360.0;
   while(azi >= 360.0)
      azi = azi - 360.0;

   fprintf(fpw,"%7.2f %6.1f %s\n",rs,azi,sname);
   }
fclose(fpr);
fclose(fpw);

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
