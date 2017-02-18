#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
char infile[1024], outfile[1024];
struct standrupformat srf0, srf1;
float mindep = -1.0e+15;
float maxdep =  1.0e+15;

int inbin = 0;
int outbin = 0;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("mindep","f",&mindep);
getpar("maxdep","f",&maxdep);
endpar();

read_srf(&srf0,infile,inbin);
select_depths_srf(&srf0,&srf1,&mindep,&maxdep);
write_srf(&srf1,outfile,outbin);
}
