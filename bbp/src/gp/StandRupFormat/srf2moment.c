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
int use_srf_mu = 0;
int moment_rate = 0;

float dt_faux = -1.0;

/*
   Average slip velocity is determined over time interval 0.0 - T, where
   T is the time when the cummulative slip first exceeds sv_avg_fac*final_slip
*/
float sv_avg_fac = 0.8;

sprintf(infile,"stdin");
sprintf(velfile,"NOT_PROVIDED");

setpar(ac,av);
getpar("infile","s",infile);
getpar("velfile","s",velfile);
getpar("use_srf_mu","d",&use_srf_mu);
getpar("moment_rate","d",&moment_rate);
getpar("dt_faux","f",&dt_faux);
getpar("sv_avg_fac","f",&sv_avg_fac);
endpar();

read_Fvelmodel(velfile,&vmod);
read_srf(&srf,infile,inbin);

if(moment_rate)
   get_moment_rate(&srf,&vmod,use_srf_mu,&dt_faux);
else
   get_moment(&srf,&vmod,use_srf_mu,&sv_avg_fac);
}
