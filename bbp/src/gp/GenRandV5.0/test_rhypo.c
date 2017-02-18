#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float flen;
struct hypo_distr_params hpar;
int nr;
long seed = 0;

hpar.x0 = 0.2;	/* default tapering starts at 20% of fault length at end end */
hpar.x1 = 0.8;
hpar.f0 = 0.1;	/* default probability at edge is 10% of probability in middle of fault */
hpar.f1 = 0.1;

setpar(ac,av);

mstpar("nr","d",&nr);
mstpar("flen","f",&flen);

getpar("seed","d",&seed);

getpar("x0","f",&hpar.x0);
getpar("f0","f",&hpar.f0);
getpar("x1","f",&hpar.x1);
getpar("f1","f",&hpar.f1);

endpar();

hpar.x0 = hpar.x0*flen;
hpar.x1 = hpar.x1*flen;
hpar.xlen = flen;
hpar.xshift = 0.0;

while(nr--)
   printf("%14.6e\n",rhypo1_lintaper(&seed,&hpar));
}
