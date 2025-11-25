#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

int main(ac,av)
int ac;
char **av;
{
struct statdata head1;
float *s1;
char infile[128];

int inbin = 0;
int outbin = 0;

sprintf(infile,"stdin");

setpar(ac,av);
getpar("infile","s",infile);
getpar("inbin","d",&inbin);
endpar();

s1 = NULL;
s1 = read_wccseis(infile,&head1,s1,inbin);

float peak = wcc_getpeak(ac, av, s1, &head1);

printf("%10.2f %13.5e %s\n",head1.edist,peak,head1.stat);
}
