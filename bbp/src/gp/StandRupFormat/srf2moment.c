#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

int main(int ac,char **av)
{
char infile[256], velfile[256];

struct standrupformat srf;
struct velmodel vmod;

int inbin = 0;

int moment_rate = 0;

sprintf(infile,"stdin");
sprintf(velfile,"NOT_PROVIDED");

setpar(ac,av);
getpar("infile","s",infile);
getpar("velfile","s",velfile);
getpar("moment_rate","d",&moment_rate);
endpar();

read_Fvelmodel(velfile,&vmod);
read_srf(&srf,infile,inbin);

if(moment_rate)
   get_moment_rate(&srf,&vmod);
else
   get_moment(&srf,&vmod);
}
