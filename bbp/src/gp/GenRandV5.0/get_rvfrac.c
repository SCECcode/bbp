#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

int main(int ac,char **av)
{
float rvfrac;
int i, n;

long seed = 0;
float range_rvfrac = 0.05;
float sigma_rvfrac = 0.1;
float mean_rvfrac = 0.75;
float fzero = 0.0;

int gaus_dist = 0;

n = 1;

setpar(ac,av);

getpar("gaus_dist","d",&gaus_dist);

getpar("mean_rvfrac","f",&mean_rvfrac);
getpar("sigma_rvfrac","f",&sigma_rvfrac);
getpar("seed","d",&seed);
getpar("n","d",&n);

endpar();

for(i=0;i<n;i++)
   {
   if(gaus_dist == 1)
      rvfrac = mean_rvfrac*exp(gaus_rand(&sigma_rvfrac,&fzero,&seed));
   else
      rvfrac = mean_rvfrac + range_rvfrac*sfrand(&seed);

   fprintf(stdout,"%.6f %.6f\n",rvfrac,mean_rvfrac);
   }
}
