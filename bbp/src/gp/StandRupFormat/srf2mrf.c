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

struct srf_apointvalues *apval_in, *apval_out;
int ip, ig;
size_t msize;

int ncoarsestk = 1;
int ncoarsedip = 1;
int stk_off = 0;
int dip_off = 0;

int inbin = 0;
int outbin = 0;

int use_srf_lame = 0;
int print_command = 1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");
sprintf(velfile,"NOT_PROVIDED");

sprintf(mrf.src_format,"MOMENT");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("velfile","s",velfile);
getpar("mrf_format","s",mrf.src_format); /* valid options: "MOMENT", ("MOMENT-1MECH"), "MOMENT-6MECH" */

getpar("ncoarsestk","d",&ncoarsestk);
getpar("ncoarsedip","d",&ncoarsedip);

stk_off = (int)(0.5*(ncoarsestk-1.0) + 0.5);
dip_off = (int)(0.5*(ncoarsedip-1.0) + 0.5);

getpar("stk_off","d",&stk_off);
getpar("dip_off","d",&dip_off);

getpar("use_srf_lame","d",&use_srf_lame);
getpar("print_command","d",&print_command);
endpar();

read_Fvelmodel(velfile,&vmod);
read_srf(&srf,infile,inbin);

fprintf(stderr,"np= %d\n",srf.srf_apnts.np);
apval_in = srf.srf_apnts.apntvals;
msize = srf.srf_apnts.np*sizeof(struct srf_apointvalues);
for(ip=0;ip<srf.srf_apnts.np;ip++)
   msize = msize + (apval_in[ip].nt1 + apval_in[ip].nt2 + apval_in[ip].nt3)*sizeof(float);

msize = msize + (int)((6.0*msize)/(1.0*ncoarsestk*ncoarsedip));

fprintf(stderr,"memory esimate (Gb)= %.4f\n",1.0e-09*msize);

if(ncoarsestk == 1 && ncoarsedip == 1)
   srf_to_mrf1(&srf,&mrf,&vmod,use_srf_lame,print_command,ac,av);
else
   srf_to_mrf6_dsamp(&srf,&mrf,&vmod,ncoarsestk,ncoarsedip,stk_off,dip_off,use_srf_lame,print_command,ac,av);

write_srf(&mrf,outfile,outbin);
}
