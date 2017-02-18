#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float arg, cosT, sinT, *xt, *yt, *ut, *vt, s1, s2;
int it, nt, i;
char infile[256], outfile[256];

struct standrupformat srf;
struct srf_apointvalues *apval_ptr;

int inbin = 0;
int outbin = 0;

float rake;
float rperd = 0.017453293;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("inbin","d",&inbin);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);
mstpar("rake","f",&rake);
endpar();

xt = NULL;
yt = NULL;

read_srf(&srf,infile,inbin);

apval_ptr = srf.srf_apnts.apntvals;
for(i=0;i<srf.srf_apnts.np;i++)
   {
   arg = (rake - apval_ptr[i].rake)*rperd;
   cosT = cos(arg);
   sinT = sin(arg);

   if(apval_ptr[i].nt1 != 0 || apval_ptr[i].nt2 != 0)
      {
      if(apval_ptr[i].nt1 != 0)
	 {
         nt = apval_ptr[i].nt1;
	 yt = (float *) check_realloc(yt,nt*sizeof(float));
	 ut = apval_ptr[i].stf1;
	 for(it=0;it<nt;it++)
	    yt[it] = ut[it];
	 }
      else
         {
         nt = apval_ptr[i].nt2;
	 yt = (float *) check_realloc(yt,nt*sizeof(float));
	 for(it=0;it<nt;it++)
	    yt[it] = 0.0;
	 }

      if(apval_ptr[i].nt2 != 0)
	 {
	 xt = (float *) check_realloc(xt,nt*sizeof(float));
	 ut = apval_ptr[i].stf2;
	 for(it=0;it<nt;it++)
	    xt[it] = ut[it];
	 }
      else
         {
	 xt = (float *) check_realloc(xt,nt*sizeof(float));
	 for(it=0;it<nt;it++)
	    xt[it] = 0.0;
	 }

      apval_ptr[i].stf1 = (float *) check_realloc(apval_ptr[i].stf1,nt*sizeof(float));
      apval_ptr[i].nt1 = nt;
      apval_ptr[i].stf2 = (float *) check_realloc(apval_ptr[i].stf2,nt*sizeof(float));
      apval_ptr[i].nt2 = nt;
      ut = apval_ptr[i].stf1;
      vt = apval_ptr[i].stf2;

      for(it=0;it<nt;it++)
         {
	 ut[it] = xt[it]*sinT + yt[it]*cosT;
	 vt[it] = xt[it]*cosT - yt[it]*sinT;
	 }
      }

   apval_ptr[i].rake = rake;

   s1 = apval_ptr[i].slip1;
   s2 = apval_ptr[i].slip2;
   apval_ptr[i].slip1 = s2*sinT + s1*cosT;
   apval_ptr[i].slip2 = s2*cosT - s1*sinT;
   }

write_srf(&srf,outfile,outbin);
}
