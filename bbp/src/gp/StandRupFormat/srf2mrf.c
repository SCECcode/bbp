#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

int main(int ac,char **av)
{
char infile[1024], outfile[1024];
struct standrupformat srf, mrf;
char velfile[1024];
struct velmodel vmod;

int inbin = 0;
int outbin = 0;

int use_srf_lame = 0;
int print_command = 1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");
sprintf(velfile,"NOT_PROVIDED");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("velfile","s",velfile);
getpar("use_srf_lame","d",&use_srf_lame);
getpar("print_command","d",&print_command);
endpar();

read_Fvelmodel(velfile,&vmod);
read_srf(&srf,infile,inbin);

srf_to_mrf(&srf,&mrf,&vmod,use_srf_lame,print_command,ac,av);

write_srf(&mrf,outfile,outbin);
}
