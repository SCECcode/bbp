#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

void read_srf(struct standrupformat *srf,char *file,int bflag);

int main(int ac,char **av)
{
int i, ip, j;
char infile[1024];

struct standrupformat srf1;
struct srf_prectsegments *prseg_ptr1;
struct srf_apointvalues *apval_ptr1;

int inbin = 0;

float hlon, hlat, hdep;
float tmin = 1.0e+15;

sprintf(infile,"stdin");

setpar(ac,av);
getpar("infile","s",infile);
endpar();

read_srf(&srf1,infile,inbin);

apval_ptr1 = srf1.srf_apnts.apntvals;

for(ip=0;ip<srf1.srf_apnts.np;ip++)
   {
   if(apval_ptr1[ip].tinit < tmin)
      {
      hlon = apval_ptr1[ip].lon;
      hlat = apval_ptr1[ip].lat;
      hdep = apval_ptr1[ip].dep;
      tmin = apval_ptr1[ip].tinit;
      }
   }

fprintf(stdout,"%.5f\t%.5f\t%.5f\n",hlon,hlat,hdep);
}
