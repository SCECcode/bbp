#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

int main(int ac,char **av)
{
char infile[1024], outfile[1024];
struct standrupformat srf_in, srf_out;

struct srf_apointvalues *apval_in;
int ip;
size_t msize;

int ncoarsestk = 1;
int ncoarsedip = 1;
int stk_off = 0;
int dip_off = 0;

int inbin = 0;
int outbin = 0;

int print_command = 1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

sprintf(srf_out.src_format,"SLIP");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);

getpar("ncoarsestk","d",&ncoarsestk);
getpar("ncoarsedip","d",&ncoarsedip);

stk_off = (int)(0.5*(ncoarsestk-1.0) + 0.5);
dip_off = (int)(0.5*(ncoarsedip-1.0) + 0.5);

getpar("stk_off","d",&stk_off);
getpar("dip_off","d",&dip_off);

getpar("print_command","d",&print_command);
endpar();

read_srf(&srf_in,infile,inbin);

fprintf(stderr,"np= %d\n",srf_in.srf_apnts.np);

apval_in = srf_in.srf_apnts.apntvals;
msize = srf_in.srf_apnts.np*sizeof(struct srf_apointvalues);
for(ip=0;ip<srf_in.srf_apnts.np;ip++)
   msize = msize + (apval_in[ip].nt1 + apval_in[ip].nt2 + apval_in[ip].nt3)*sizeof(float);

msize = msize + (int)((1.0*msize)/(1.0*ncoarsestk*ncoarsedip));

fprintf(stderr,"memory esimate (Gb)= %.4f\n",1.0e-09*msize);

srf_dwnsamp(&srf_in,&srf_out,ncoarsestk,ncoarsedip,stk_off,dip_off,print_command,ac,av);

write_srf(&srf_out,outfile,outbin);
}
