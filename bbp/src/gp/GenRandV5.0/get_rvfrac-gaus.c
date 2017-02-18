#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float rvfrac;
int i, n;

long seed = 0;
float lnsigma_rvfrac = 0.1;
float median_rvfrac = 0.75;
float fzero = 0.0;

n = 1;

setpar(ac,av);

getpar("median_rvfrac","f",&median_rvfrac);
getpar("lnsigma_rvfrac","f",&lnsigma_rvfrac);
getpar("seed","d",&seed);
getpar("n","d",&n);

endpar();

for(i=0;i<n;i++)
   {
   rvfrac = median_rvfrac*exp(gaussian_rand(&lnsigma_rvfrac,&fzero,&seed));
   fprintf(stdout,"%.6f %.6f\n",median_rvfrac,rvfrac);
   }
}
