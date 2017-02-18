#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float rvfrac;
int i, n;

long seed = 0;
float range_rvfrac = 0.05;
float mean_rvfrac = 0.75;
float fzero = 0.0;

n = 1;

setpar(ac,av);

getpar("mean_rvfrac","f",&mean_rvfrac);
getpar("range_rvfrac","f",&range_rvfrac);
getpar("seed","d",&seed);
getpar("n","d",&n);

endpar();

for(i=0;i<n;i++)
   {
   rvfrac = mean_rvfrac + range_rvfrac*sfrand(&seed);
   fprintf(stdout,"%.6f %.6f\n",mean_rvfrac,rvfrac);
   }
}
