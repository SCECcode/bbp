#include "include.h"
#include "structure.h"
#include "function.h"

void get_srfpars(struct standrupformat *srf,int off, int ip,float *rt,float *vs,float *stk,float *dip,float *rak,struct mechparam *mpar)
{
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

mpar->nmech = 0;
mpar->flag[0] = 0;
mpar->flag[1] = 0;
mpar->flag[2] = 0;

if(apval_ptr[ip].nt1 > 0)
   {
   mpar->flag[mpar->nmech] = U1FLAG;
   mpar->nmech = mpar->nmech + 1;
   }
if(apval_ptr[ip].nt2 > 0)
   {
   mpar->flag[mpar->nmech] = U2FLAG;
   mpar->nmech = mpar->nmech + 1;
   }
if(apval_ptr[ip].nt3 > 0)
   {
   mpar->flag[mpar->nmech] = U3FLAG;
   mpar->nmech = mpar->nmech + 1;
   }

*vs = sqrt(apval_ptr[ip].slip1*apval_ptr[ip].slip1
         + apval_ptr[ip].slip2*apval_ptr[ip].slip2
	 + apval_ptr[ip].slip3*apval_ptr[ip].slip3);
*stk = apval_ptr[ip].stk;
*dip = apval_ptr[ip].dip;
*rak = apval_ptr[ip].rake;
*rt = apval_ptr[ip].tinit;
}

void srf_stf(struct standrupformat *srf,int off,int ip,float *s,float *u,float *stf,int nt,float *dt,struct mechparam mp,float *unused)
{
FILE *fpw;
int it, nstf, im;
float sum, *sptr, *uptr;
float *space, *sptr2;

struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float fnt;
int resamp, ntpad, ntrsmp, gnt;

float tol = 1.0e-02;

float pratio_tol, mratio_tol;
float ratio_tol = 0.00001;

pratio_tol = 1.0 + ratio_tol;
mratio_tol = 1.0 - ratio_tol;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

zapit(stf,nt);

space = NULL;
sptr2 = NULL;

for(im=0;im<mp.nmech;im++)
   {
   if(mp.flag[im] == U1FLAG)
      {
      nstf = apval_ptr[ip].nt1;
      sptr = apval_ptr[ip].stf1;
      }
   else if(mp.flag[im] == U2FLAG)
      {
      nstf = apval_ptr[ip].nt2;
      sptr = apval_ptr[ip].stf2;
      }
   else if(mp.flag[im] == U3FLAG)
      {
      nstf = apval_ptr[ip].nt3;
      sptr = apval_ptr[ip].stf3;
      }

/* resample if needed */
   if((*dt)/apval_ptr[ip].dt > pratio_tol || (*dt)/apval_ptr[ip].dt < mratio_tol)
      {
      ntpad = 2*nstf;
      fnt = ntpad*apval_ptr[ip].dt/(*dt);
      gnt = (int)(fnt + 0.5);

      while(nt_tol(fnt,gnt) > tol)
         {
         ntpad++;
         fnt = ntpad*apval_ptr[ip].dt/(*dt);
         gnt = (int)(fnt + 0.5);
         }

      ntrsmp = (int)(fnt);

      if((*dt) < apval_ptr[ip].dt)
	 {
	 space = (float *) check_realloc ((void *)space,2*ntrsmp*sizeof(float));
	 sptr2 = (float *)check_realloc((void *)sptr2,2*ntrsmp*sizeof(float));

         for(it=0;it<2*ntrsmp;it++)
            sptr2[it] = 0.0;

         resamp = 1;
	 }
      else
	 {
	 space = (float *) check_realloc ((void *)space,2*ntpad*sizeof(float));
	 sptr2 = (float *)check_realloc((void *)sptr2,2*ntpad*sizeof(float));

         for(it=0;it<2*ntpad;it++)
            sptr2[it] = 0.0;

         resamp = -1;
	 }

      for(it=0;it<nstf;it++)
         sptr2[it] = sptr[it];

      resample(sptr2,nstf,&apval_ptr[ip].dt,resamp,ntpad,ntrsmp,dt,space);

      sptr = sptr2;
      nstf = ntrsmp;
      }

   if(nstf > nt)
      nstf = nt;

/* add factor of dt to prenormalize convolution */
   for(it=0;it<nstf;it++)
      stf[it] = (*dt)*sptr[it];

   uptr = u + 3*im*nt;
   do_cnvlv(s,uptr,nt,stf,nstf);
   }

free(space);
free(sptr2);
}

void get_srfpars_v2(struct standrupformat *srf,int off, int ip,float *rt,float *vs,struct mechparam *mpar)
{
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

mpar->nmech = 0;
mpar->flag[0] = 0;
mpar->flag[1] = 0;
mpar->flag[2] = 0;

if(apval_ptr[ip].nt1 > 0)
   {
   mpar->flag[mpar->nmech] = U1FLAG;
   mpar->nmech = mpar->nmech + 1;
   }
if(apval_ptr[ip].nt2 > 0)
   {
   mpar->flag[mpar->nmech] = U2FLAG;
   mpar->nmech = mpar->nmech + 1;
   }
if(apval_ptr[ip].nt3 > 0)
   {
   mpar->flag[mpar->nmech] = U3FLAG;
   mpar->nmech = mpar->nmech + 1;
   }

*vs = sqrt(apval_ptr[ip].slip1*apval_ptr[ip].slip1
         + apval_ptr[ip].slip2*apval_ptr[ip].slip2
	 + apval_ptr[ip].slip3*apval_ptr[ip].slip3);

mpar->stk = apval_ptr[ip].stk;
mpar->dip = apval_ptr[ip].dip;
mpar->rak = apval_ptr[ip].rake;

*rt = apval_ptr[ip].tinit;
}
