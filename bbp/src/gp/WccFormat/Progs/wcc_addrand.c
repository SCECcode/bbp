#include "include.h"
#include "structure.h"
#include "function.h"

double frand(void);
double sfrand(long *);

int size_float = sizeof(float);
int size_int = sizeof(int);

main(int ac,char **av)
{
struct statdata shead1;
float *s1;

char infile[128];
char outfile[128];

int inbin = 0;
int outbin = 0;

int it;
float add_rand = 0.0;

setpar(ac,av);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("inbin","d",&inbin);
getpar("outbin","d",&outbin);
mstpar("add_rand","f",&add_rand);
endpar();

s1 = NULL;
s1 = read_wccseis(infile,&shead1,s1,inbin);

if(add_rand > (float)(0.0))
   {
   for(it=0;it<shead1.nt;it++)
      s1[it] = s1[it] + add_rand*frand();
   }

write_wccseis(outfile,&shead1,s1,outbin);
}

static  long    frandx = 1;

/* frand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double frand(void)
{
frandx = (frandx * 1103515245 + 12345) & 0x7fffffff;
return((double)(frandx)/1073741824.0 - 1.0);
}

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}
