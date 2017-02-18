#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float *val, avg, csigma;
int i, n;

long seed = 0;
float sigma = 1.0;
float fzero = 0.0;
float fone = 1.0;

n = 1;

setpar(ac,av);

getpar("sigma","f",&sigma);
getpar("seed","d",&seed);
getpar("n","d",&n);

endpar();

val = (float *)check_malloc(n*sizeof(float));

avg = 0.0;
for(i=0;i<n;i++)
   {
   val[i] = sigma*gaus_rand(&fone,&fzero,&seed);
   avg = avg + val[i];
   /*
   fprintf(stdout,"%.6f %.6f\n",val[i],avg/i);
   */
   }

avg = avg/n;
csigma = 0.0;
for(i=0;i<n;i++)
   csigma = csigma + (avg-val[i])*(avg-val[i]);

csigma = sqrt(csigma/n);
fprintf(stdout,"%.6f %.6f\n",csigma,sigma);
}
