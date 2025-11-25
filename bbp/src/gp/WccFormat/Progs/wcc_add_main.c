#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

double frand(void);
double sfrand(long *);

void sum(float *, struct statdata *, float *, float *, float *, struct statdata *, float *, float *, float *, struct statdata *);

int main(int ac,char **av)
{
struct statdata shead1, shead2, shead3;
float* s1, *s2, *p;

int size_float = sizeof(float);

char infile1[256];
char infile2[256];
char outfile[256];

int inbin1 = 0;
int inbin2 = 0;
int outbin = 0;

setpar(ac,av);
mstpar("infile1","s",infile1);
mstpar("infile2","s",infile2);
mstpar("outfile","s",outfile);
getpar("inbin1","d",&inbin1);
getpar("inbin2","d",&inbin2);
getpar("outbin","d",&outbin);
endpar();

s1 = NULL;
s1 = read_wccseis(infile1,&shead1,s1,inbin1);
s2 = NULL;
s2 = read_wccseis(infile2,&shead2,s2,inbin2);

p = (float *) check_malloc ((shead1.nt+shead2.nt)*size_float);

wcc_add(ac, av, s1, &shead1, s2, &shead2, p, &shead3);

strcpy(shead3.stat,shead1.stat);
strcpy(shead3.comp,shead1.comp);
sprintf(shead3.stitle,"summed output");

write_wccseis(outfile,&shead3,p,outbin);
}

