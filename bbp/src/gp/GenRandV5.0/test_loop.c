#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
double ran;
long k, nr;
long seed = 0;

while(scanf("%d",&nr) == 1)
   {
   printf("nr= %10d: START ... ",nr);
   fflush(stdout);

   for(k=0;k<100*nr;k++)
      ran = sfrand(&seed);

   printf("END\n");
   fflush(stdout);
   }
}
