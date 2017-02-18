#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float *stf1, *stf2;
int it, i, ip, j, k, ig, ik;
char infile[256], outfile[256];

struct standrupformat srf1, srf2;
struct srf_prectsegments *prseg_ptr1, *prseg_ptr2;
struct srf_apointvalues *apval_ptr1, *apval_ptr2;

int inbin = 0;
int outbin = 0;

int npstart = -1;
int npend = -1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("inbin","d",&inbin);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);
getpar("npstart","d",&npstart);
getpar("npend","d",&npend);
endpar();

fprintf(stderr,"READ started ... ");
read_srf(&srf1,infile,inbin);
fprintf(stderr,"DONE\n");

strcpy(srf2.version,srf1.version);

if(strncmp(srf1.type,"PLANE",5) == 0)
   {
   strcpy(srf2.type,srf1.type);
   srf2.srf_prect.nseg = srf1.srf_prect.nseg;

   srf2.srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2.srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_ptr1 = srf1.srf_prect.prectseg;
   prseg_ptr2 = srf2.srf_prect.prectseg;

   for(ig=0;ig<srf2.srf_prect.nseg;ig++)
      {
      prseg_ptr2[ig].elon = prseg_ptr1[ig].elon;
      prseg_ptr2[ig].elat = prseg_ptr1[ig].elat;
      prseg_ptr2[ig].nstk = prseg_ptr1[ig].nstk;
      prseg_ptr2[ig].ndip = prseg_ptr1[ig].ndip;
      prseg_ptr2[ig].flen = prseg_ptr1[ig].flen;
      prseg_ptr2[ig].fwid = prseg_ptr1[ig].fwid;
      prseg_ptr2[ig].stk  = prseg_ptr1[ig].stk;
      prseg_ptr2[ig].dip  = prseg_ptr1[ig].dip;
      prseg_ptr2[ig].dtop = prseg_ptr1[ig].dtop;
      prseg_ptr2[ig].shyp = prseg_ptr1[ig].shyp;
      prseg_ptr2[ig].dhyp = prseg_ptr1[ig].dhyp;
      }
   }

srf2.srf_apnts.np = srf1.srf_apnts.np;
srf2.srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2.srf_apnts.np)*sizeof(struct srf_apointvalues));

apval_ptr1 = srf1.srf_apnts.apntvals;
apval_ptr2 = srf2.srf_apnts.apntvals;

if(npstart < 0)
   npstart = 0;
if(npend < 0 || npend > srf1.srf_apnts.np)
   npend = srf1.srf_apnts.np;

ip = 0;
for(i=0;i<srf1.srf_apnts.np;i++)
   {
   fprintf(stderr,"i= %d ip= %d\n",i,ip);

   if(i >= npstart && i < npend)
      {
	 apval_ptr2[ip].lon = apval_ptr1[i].lon;
	 apval_ptr2[ip].lat = apval_ptr1[i].lat;
         apval_ptr2[ip].dep   = apval_ptr1[i].dep;

         apval_ptr2[ip].stk   = apval_ptr1[i].stk;
         apval_ptr2[ip].dip   = apval_ptr1[i].dip;
         apval_ptr2[ip].area  = apval_ptr1[i].area;

	 apval_ptr2[ip].tinit = apval_ptr1[i].tinit;

         apval_ptr2[ip].dt    = apval_ptr1[i].dt;
         apval_ptr2[ip].rake  = apval_ptr1[i].rake;
         apval_ptr2[ip].slip1 = apval_ptr1[i].slip1;
         apval_ptr2[ip].nt1   = apval_ptr1[i].nt1;
         apval_ptr2[ip].slip2 = apval_ptr1[i].slip2;
         apval_ptr2[ip].nt2   = apval_ptr1[i].nt2;
         apval_ptr2[ip].slip3 = apval_ptr1[i].slip3;
         apval_ptr2[ip].nt3   = apval_ptr1[i].nt3;

         apval_ptr2[ip].stf1 = (float *)check_malloc((apval_ptr2[ip].nt1)*sizeof(float));
         apval_ptr2[ip].stf2 = (float *)check_malloc((apval_ptr2[ip].nt2)*sizeof(float));
         apval_ptr2[ip].stf3 = (float *)check_malloc((apval_ptr2[ip].nt3)*sizeof(float));

	 stf1 = apval_ptr1[i].stf1;
	 stf2 = apval_ptr2[ip].stf1;
	 for(it=0;it<apval_ptr2[ip].nt1;it++)
	    stf2[it] = stf1[it];

	 stf1 = apval_ptr1[i].stf2;
	 stf2 = apval_ptr2[ip].stf2;
	 for(it=0;it<apval_ptr2[ip].nt2;it++)
	    stf2[it] = stf1[it];

	 stf1 = apval_ptr1[i].stf3;
	 stf2 = apval_ptr2[ip].stf3;
	 for(it=0;it<apval_ptr2[ip].nt3;it++)
	    stf2[it] = stf1[it];

      ip++;
      }
   }
srf2.srf_apnts.np = ip;

write_srf(&srf2,outfile,outbin);
}
