#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float *stf1, *stf2;
int it, i, ip, ig, ix, iy, ixp, iyp;
int ncoarsestk, ncoarsedip;
int ntot1, ntot2, *new_nstk, *new_ndip, *old_nstk, *old_ndip;
char infile[256], outfile[256];

struct standrupformat srf1, srf2;
struct srf_prectsegments *prseg_ptr1, *prseg_ptr2;
struct srf_apointvalues *apval_ptr1, *apval_ptr2;

int inbin = 0;
int outbin = 0;

int scale_moment = 0;
double mom1, mom2;
float vslip, momfac;

int stk_off = 0;
int dip_off = 0;

int print_command = 1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("inbin","d",&inbin);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);
mstpar("ncoarsestk","d",&ncoarsestk);
mstpar("ncoarsedip","d",&ncoarsedip);
getpar("stk_off","d",&stk_off);
getpar("dip_off","d",&dip_off);

getpar("scale_moment","d",&scale_moment);

getpar("print_command","d",&print_command);

endpar();

read_srf(&srf1,infile,inbin);

strcpy(srf2.version,srf1.version);
copy_hcmnt(&srf2,&srf1);

if(print_command && atof(srf2.version) >= 2.0)
   load_command_srf(&srf2,ac,av);

if(strncmp(srf1.type,"PLANE",5) == 0)
   {
   strcpy(srf2.type,srf1.type);
   srf2.srf_prect.nseg = srf1.srf_prect.nseg;

   srf2.srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2.srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_ptr1 = srf1.srf_prect.prectseg;
   prseg_ptr2 = srf2.srf_prect.prectseg;

   old_nstk = (int *)check_malloc((srf1.srf_prect.nseg)*sizeof(int));
   old_ndip = (int *)check_malloc((srf1.srf_prect.nseg)*sizeof(int));
   new_nstk = (int *)check_malloc((srf2.srf_prect.nseg)*sizeof(int));
   new_ndip = (int *)check_malloc((srf2.srf_prect.nseg)*sizeof(int));

   for(ig=0;ig<srf2.srf_prect.nseg;ig++)
      {
      prseg_ptr2[ig].elon = prseg_ptr1[ig].elon;
      prseg_ptr2[ig].elat = prseg_ptr1[ig].elat;

      prseg_ptr2[ig].nstk = (int)((1.0*prseg_ptr1[ig].nstk/ncoarsestk + 0.5));
      while(prseg_ptr2[ig].nstk*ncoarsestk > prseg_ptr1[ig].nstk)
         prseg_ptr2[ig].nstk--;

      prseg_ptr2[ig].ndip = (int)((1.0*prseg_ptr1[ig].ndip/ncoarsedip + 0.5));
      while(prseg_ptr2[ig].ndip*ncoarsedip > prseg_ptr1[ig].ndip)
         prseg_ptr2[ig].ndip--;

      prseg_ptr2[ig].flen = prseg_ptr2[ig].nstk*ncoarsestk*(prseg_ptr1[ig].flen/prseg_ptr1[ig].nstk);
      prseg_ptr2[ig].fwid = prseg_ptr2[ig].ndip*ncoarsedip*(prseg_ptr1[ig].fwid/prseg_ptr1[ig].ndip);
      prseg_ptr2[ig].stk  = prseg_ptr1[ig].stk;
      prseg_ptr2[ig].dip  = prseg_ptr1[ig].dip;
      prseg_ptr2[ig].dtop = prseg_ptr1[ig].dtop;
      prseg_ptr2[ig].shyp = prseg_ptr1[ig].shyp;
      prseg_ptr2[ig].dhyp = prseg_ptr1[ig].dhyp;

      old_nstk[ig] = prseg_ptr1[ig].nstk;
      old_ndip[ig] = prseg_ptr1[ig].ndip;
      new_nstk[ig] = prseg_ptr2[ig].nstk;
      new_ndip[ig] = prseg_ptr2[ig].ndip;
      }
   }

srf2.nseg = srf1.nseg;
srf2.np_seg = (int *)check_malloc((srf2.nseg)*sizeof(int));

if(atof(srf2.version) < 2.0)
   {
   for(ig=1;ig<srf2.nseg;ig++)
      {
      old_nstk[0] = old_nstk[0] + old_nstk[ig];
      new_nstk[0] = new_nstk[0] + new_nstk[ig];
      }

   srf2.np_seg[0] = new_nstk[0]*new_ndip[0];
   srf2.srf_apnts.np = new_nstk[0]*new_ndip[0];
   }

else if(atof(srf2.version) >= 2.0)
   {
   srf2.srf_apnts.np = 0;
   for(ig=0;ig<srf2.nseg;ig++)
      {
      srf2.np_seg[ig] = new_nstk[ig]*new_ndip[ig];
      srf2.srf_apnts.np = srf2.srf_apnts.np + srf2.np_seg[ig];
      }
   }

srf2.srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2.srf_apnts.np)*sizeof(struct srf_apointvalues));

apval_ptr1 = srf1.srf_apnts.apntvals;
apval_ptr2 = srf2.srf_apnts.apntvals;

if(scale_moment == 1 && atof(srf2.version) >= 2.0)
   {
   ntot1 = 0;
   mom1 = 0.0;
   for(ig=0;ig<srf1.nseg;ig++)
      {
      for(iy=0;iy<old_ndip[ig];iy++)
         {
         for(ix=0;ix<old_nstk[ig];ix++)
            {
            i = ix + iy*old_nstk[ig] + ntot1;

	    vslip = sqrt(apval_ptr1[i].slip1*apval_ptr1[i].slip1 + apval_ptr1[i].slip2*apval_ptr1[i].slip2 + apval_ptr1[i].slip3*apval_ptr1[i].slip3);
	    mom1 = mom1 + vslip*apval_ptr1[i].vs*apval_ptr1[i].vs*apval_ptr1[i].den*apval_ptr1[i].area;
            }
         }
      ntot1 = ntot1 + srf1.np_seg[ig];
      }
   }

ntot1 = 0;
ntot2 = 0;
mom2 = 0.0;
for(ig=0;ig<srf2.nseg;ig++)
   {
   for(iy=0;iy<new_ndip[ig];iy++)
      {
      /*
      iyp = (int)((float)(iy + 0.5)*ncoarsedip) + dip_off;
      */
      iyp = iy*ncoarsedip + dip_off;
      for(ix=0;ix<new_nstk[ig];ix++)
         {
	 /*
         ixp = (int)((float)(ix + 0.5)*ncoarsestk) + stk_off;
	 */
         ixp = ix*ncoarsestk + stk_off;

         i = ix + iy*new_nstk[ig] + ntot2;
         ip = ixp + iyp*old_nstk[ig] + ntot1;

         apval_ptr2[i].lon   = apval_ptr1[ip].lon;
         apval_ptr2[i].lat   = apval_ptr1[ip].lat;
         apval_ptr2[i].dep   = apval_ptr1[ip].dep;
         apval_ptr2[i].stk   = apval_ptr1[ip].stk;
         apval_ptr2[i].dip   = apval_ptr1[ip].dip;
         apval_ptr2[i].area  = apval_ptr1[ip].area*(ncoarsestk*ncoarsedip);
         apval_ptr2[i].tinit = apval_ptr1[ip].tinit;
         apval_ptr2[i].dt    = apval_ptr1[ip].dt;
         apval_ptr2[i].vs    = apval_ptr1[ip].vs;
         apval_ptr2[i].den    = apval_ptr1[ip].den;

         apval_ptr2[i].rake  = apval_ptr1[ip].rake;
         apval_ptr2[i].slip1 = apval_ptr1[ip].slip1;
         apval_ptr2[i].nt1   = apval_ptr1[ip].nt1;
         apval_ptr2[i].slip2 = apval_ptr1[ip].slip2;
         apval_ptr2[i].nt2   = apval_ptr1[ip].nt2;
         apval_ptr2[i].slip3 = apval_ptr1[ip].slip3;
         apval_ptr2[i].nt3   = apval_ptr1[ip].nt3;

         apval_ptr2[i].stf1 = (float *)check_malloc((apval_ptr2[i].nt1)*sizeof(float));
         apval_ptr2[i].stf2 = (float *)check_malloc((apval_ptr2[i].nt2)*sizeof(float));
         apval_ptr2[i].stf3 = (float *)check_malloc((apval_ptr2[i].nt3)*sizeof(float));

         stf1 = apval_ptr1[ip].stf1;
         stf2 = apval_ptr2[i].stf1;
         for(it=0;it<apval_ptr2[i].nt1;it++)
            stf2[it] = stf1[it];

         stf1 = apval_ptr1[ip].stf2;
         stf2 = apval_ptr2[i].stf2;
         for(it=0;it<apval_ptr2[i].nt2;it++)
            stf2[it] = stf1[it];

         stf1 = apval_ptr1[ip].stf3;
         stf2 = apval_ptr2[i].stf3;
         for(it=0;it<apval_ptr2[i].nt3;it++)
            stf2[it] = stf1[it];

	 vslip = sqrt(apval_ptr2[i].slip1*apval_ptr2[i].slip1 + apval_ptr2[i].slip2*apval_ptr2[i].slip2 + apval_ptr2[i].slip3*apval_ptr2[i].slip3);
	 mom2 = mom2 + vslip*apval_ptr2[i].vs*apval_ptr2[i].vs*apval_ptr2[i].den*apval_ptr2[i].area;
         }
      }
   ntot1 = ntot1 + srf1.np_seg[ig];
   ntot2 = ntot2 + srf2.np_seg[ig];
   }

if(scale_moment == 1 && atof(srf2.version) >= 2.0)
   {
   ntot2 = 0;
   momfac = mom1/mom2;
   for(ig=0;ig<srf2.nseg;ig++)
      {
      for(iy=0;iy<new_ndip[ig];iy++)
         {
         for(ix=0;ix<new_nstk[ig];ix++)
            {
            i = ix + iy*new_nstk[ig] + ntot2;

            apval_ptr2[i].slip1 = momfac*apval_ptr2[i].slip1;
            apval_ptr2[i].slip2 = momfac*apval_ptr2[i].slip2;
            apval_ptr2[i].slip3 = momfac*apval_ptr2[i].slip3;

            stf2 = apval_ptr2[i].stf1;
            for(it=0;it<apval_ptr2[i].nt1;it++)
               stf2[it] = momfac*stf2[it];

            stf2 = apval_ptr2[i].stf2;
            for(it=0;it<apval_ptr2[i].nt2;it++)
               stf2[it] = momfac*stf2[it];

            stf2 = apval_ptr2[i].stf3;
            for(it=0;it<apval_ptr2[i].nt3;it++)
               stf2[it] = momfac*stf2[it];
            }
         }
      ntot2 = ntot2 + srf2.np_seg[ig];
      }

   fprintf(stderr,"momfac= %f\n",momfac);
   }

write_srf(&srf2,outfile,outbin);
}
