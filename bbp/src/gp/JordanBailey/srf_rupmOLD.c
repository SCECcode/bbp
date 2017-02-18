#include "include.h"
#include "structure.h"
#include "function.h"

void get_srfpars(struct standrupformat *srf,int i,int j,float *rt,float *vs,float *rk,int nseg,int off,struct mechparam *mpar)
{
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
int ip;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

ip = i + j*(prseg_ptr[nseg].nstk);

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
*rk = apval_ptr[ip].rake;
*rt = apval_ptr[ip].tinit;
}

void srf_stf(struct standrupformat *srf,int i,int j,float *s,float *u,float *stf,int nt,float *dt,int nseg,int off,struct mechparam mp)
{
FILE *fpw;
int ip, it, nstf, im;
float sum, *sptr, *uptr;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

zapit(stf,nt);

ip = i + j*(prseg_ptr[nseg].nstk);

if(apval_ptr[ip].nt1 == 0)
   return;

/* for now, simply copy STF
   should add option to resample to dtout
   */

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

/* add factor of dt to prenormalize convolution */
   for(it=0;it<nstf;it++)
      stf[it] = (*dt)*sptr[it];

   uptr = u + 3*im*nt;
   do_cnvlv(s,uptr,nt,stf,nstf);
   }
}
