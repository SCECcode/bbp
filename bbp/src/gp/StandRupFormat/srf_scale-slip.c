#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
char infile[256], outfile[256];
struct standrupformat srf1, srf2;

int inbin = 0;
int outbin = 0;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

float scale = 1.0;

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("scale","f",&scale);
endpar();

read_srf(&srf1,infile,inbin);

scale_srf(&srf1,&srf2,&scale);

write_srf(&srf2,outfile,outbin);
}
