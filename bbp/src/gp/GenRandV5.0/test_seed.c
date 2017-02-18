#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
double ran;
long seed = 0;

while(scanf("%d",&seed) == 1)
   {
   printf("%20d %20d : ",seed,((seed) * 1103515245 + 12345));
   seed = ((seed) * 1103515245 + 12345) & 0x7fffffff;
   ran = ((double)(seed)/1073741824.0 - 1.0);

   printf("%20d %14.6e\n",seed,ran);
   }
}
