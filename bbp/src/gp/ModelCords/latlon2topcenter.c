#include <sys/file.h>
#include <stdio.h>
#include <math.h>

main(ac,av)
int ac;
char **av;
{
float elat, elon, mlat, mlon;
float kperd_n, kperd_e, xtop, ytop;
float x0, y0, cosR, sinR;
double cosA, sinA, e2, den, g2, h2, lat0;

float rotate = 0.0;
double rperd = 0.017453292;
double rad = 6378.139;
double f = 298.256;

setpar(ac, av);
mstpar("elat","f",&elat);
mstpar("elon","f",&elon);
mstpar("mlat","f",&mlat);
mstpar("mlon","f",&mlon);
getpar("rotate","f",&rotate);
endpar();

f = 1.0/f;
e2 = 2.0*f - f*f;
h2 = (1.0 - f)*(1.0 - f);
g2 = e2/h2;

/* convert geographical latitude to geocentric latitude */

lat0 = atan((1.0 - f)*tan(0.5*(elat+mlat)*rperd));

cosA = cos(lat0);
sinA = sin(lat0);

den = sqrt(1.0/(1.0 + g2*sinA*sinA));
kperd_e = rperd*rad*cosA*den;
kperd_n = rperd*rad*(sqrt(1.0 + g2*sinA*sinA*(2.0 + g2)))*den*den*den;

cosR = cos(rotate*rperd);
sinR = sin(rotate*rperd);

x0 = (elon - mlon)*kperd_e;
y0 = (mlat - elat)*kperd_n;
xtop = x0*cosR + y0*sinR;
ytop = -x0*sinR + y0*cosR;

printf("%10.4f %10.4f\n",xtop,ytop);
}
