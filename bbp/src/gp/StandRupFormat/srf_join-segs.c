#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
char infile1[1024], infile2[1024], outfile[1024];
struct standrupformat srf1, srf2, srf3;

int inbin = 0;
int outbin = 0;

sprintf(outfile,"stdout");

setpar(ac,av);
mstpar("infile1","s",infile1);
mstpar("infile2","s",infile2);
getpar("outfile","s",outfile);
endpar();

read_srf(&srf1,infile1,inbin);
read_srf(&srf2,infile2,inbin);

join_segs(&srf1,&srf2,&srf3);

write_srf(&srf3,outfile,outbin);
}
