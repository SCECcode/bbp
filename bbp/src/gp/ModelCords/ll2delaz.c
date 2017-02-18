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
double slat, slon;
double elat, elon;
double alpha, beta, cosA, cosB, sinB, sinC, dperr;
float del, azi;

char infile[512], outfile[512], str[512], sname[8];

double rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;
double geocen();

double pi = 3.141592654;

int use_all = 0;
int printall = 1;
int xbnd = 0;
int ybnd = 0;

setpar(ac, av);
mstpar("elat","F",&elat);
mstpar("elon","F",&elon);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
endpar();

dperr = 1.0/rperd;
radc = ERAD*RPERD;
set_g2(&g2,&fc);

if(elon < 0.0)
   elon = 360.0 + elon;

elon = elon*rperd;
elat = geocen(elat*rperd);

fpr = fopfile(infile,"r");
fpw = fopfile(outfile,"w");

while(fgets(str,512,fpr) != NULL)
   {
   sscanf(str,"%F %F %s",&slon,&slat,sname);

   if(slon < 0.0)
      slon = 360.0 + slon;

   alpha = geocen(slat*rperd) - elat;
   beta = slon*rperd - elon;

   if(beta > pi)
      beta = beta - 2.0*pi;
   else if(beta < -pi)
      beta = beta + 2.0*pi;

   if(alpha == 0.0 && beta == 0.0)   /* points are the same */
      {
      del = 0.0;
      azi = 0.0;
      }
   else if(alpha == 0.0)    /* lats are the same */
      {
      fprintf(stderr,"LATS same\n");
      del = beta;
      azi = 90.0;
      if(del < 0.0)
	 {
	 del = -del;
	 azi = 270.0;
	 }
      }
   else if(beta == 0.0)    /* lons are the same */
      {
      fprintf(stderr,"LONS same\n");
      del = alpha;
      azi = 0.0;
      if(del < 0.0)
	 {
	 del = -del;
	 azi = 180.0;
	 }
      }
   else    /* everything else */
      {
      cosA = cos(alpha);
      cosB = cos(beta);

      sinB = sin(beta);
      if(beta < 0.0)
	 sinB = -sinB;

      sinC = sqrt(1.0 - cosA*cosA*cosB*cosB);

      del = dperr*acos(cosA*cosB);
      del = acos(cosA*cosB);
      azi = dperr*acos(sinB*cosA/sinC);

      if(alpha > 0.0 && beta > 0.0)
	 azi = 90.0 - azi;
      else if(alpha < 0.0 && beta > 0.0)
	 azi = 90.0 + azi;
      else if(alpha < 0.0 && beta < 0.0)
	 azi = 270.0 - azi;
      else if(alpha > 0.0 && beta < 0.0)
	 azi = 270.0 + azi;
      }

   fprintf(fpw,"%7.2f %6.1f %s %10.5f %10.5f\n",erad*del,azi,sname,alpha*dperr,beta*dperr);
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
