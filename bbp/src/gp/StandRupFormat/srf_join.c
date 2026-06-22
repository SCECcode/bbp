#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

int main(int ac,char **av)
{
char infile1[1024], infile2[1024], outfile[1024];
struct standrupformat srf1, srf2, srf3;

int inbin = 0;
int outbin = 0;

int print_command = 1;

sprintf(outfile,"stdout");

setpar(ac,av);
mstpar("infile1","s",infile1);
mstpar("infile2","s",infile2);
getpar("outfile","s",outfile);
getpar("print_command","d",&print_command);
endpar();

read_srf(&srf1,infile1,inbin);
read_srf(&srf2,infile2,inbin);

join_srf(&srf1,&srf2,&srf3,print_command,ac,av);

write_srf(&srf3,outfile,outbin);
}
