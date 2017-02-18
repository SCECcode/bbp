#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float xx, yy, vrup, inv_vrup, dstk, ddip, len2;
int i, ip, j;
char infile[256], outfile[256];

struct standrupformat srf1;
struct srf_prectsegments *prseg_ptr1;
struct srf_apointvalues *apval_ptr1;

int inbin = 0;
int outbin = 0;

int nstk = -1;
int ndip = -1;
float flen = -1.0;
float fwid = -1.0;
float shypo = -1.0e+15;
float dhypo = -1.0e+15;

float savg, smax, rt, sf, svec;
struct velmodel rvmod;
double rayp, rupt_rad;

float dtop = -1.0;
float dip = -1.0;
struct velmodel vmod;
char velfile[128];

float rvfrac = DEFAULT_VR_TO_VS_FRAC;
float shal_vrup = DEFAULT_SHAL_VRUP_FRAC;
float tsfac = DEFAULT_TSFAC;
float htol = 0.1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

velfile[0] = '\0';

setpar(ac,av);
getpar("infile","s",infile);
getpar("inbin","d",&inbin);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);
getpar("nstk","d",&nstk);
getpar("ndip","d",&ndip);
getpar("flen","f",&flen);
getpar("fwid","f",&fwid);
getpar("dtop","f",&dtop);
getpar("dip","f",&dip);
getpar("shypo","f",&shypo);
getpar("dhypo","f",&dhypo);

mstpar("vrup","f",&vrup);
if(vrup < 0.0)
   {
   getpar("rvfrac","f",&rvfrac);
   getpar("shal_vrup","f",&shal_vrup);
   getpar("tsfac","f",&tsfac);
   getpar("velfile","s",velfile);
   }
endpar();

read_srf(&srf1,infile,inbin);

if(strncmp(srf1.type,"PLANE",5) == 0)
   {
   prseg_ptr1 = srf1.srf_prect.prectseg;

   if(nstk < 0.0)
      nstk = prseg_ptr1[0].nstk;
   if(ndip < 0.0)
      ndip = prseg_ptr1[0].ndip;
   if(flen < 0.0)
      flen = prseg_ptr1[0].flen;
   if(fwid < 0.0)
      fwid = prseg_ptr1[0].fwid;
   if(dtop < 0.0)
      dtop = prseg_ptr1[0].dtop;
   if(dip < 0.0)
      dip = prseg_ptr1[0].dip;
   if(shypo < -1.0e+14)
      shypo = prseg_ptr1[0].shyp;
   if(dhypo < -1.0e+14)
      dhypo = prseg_ptr1[0].dhyp;
   }

dstk = flen/nstk;
ddip = fwid/ndip;
len2 = 0.5*flen;
inv_vrup = 1.0/vrup;

if(nstk*ndip != srf1.srf_apnts.np)
   {
   fprintf(stderr,"problem with number of points, exiting...\n");
   exit(-1);
   }

apval_ptr1 = srf1.srf_apnts.apntvals;

if(vrup < 0.0)
   {
   if(velfile[0] != '\0')
      read_velmodel(velfile,&vmod);
   else
      default_velmodel(&vmod);

   conv2vrup(&vmod,&rvmod,&dip,&dtop,&fwid,&rvfrac,&shal_vrup);

   savg = 0.0;
   smax = -1.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
         ip = i + j*nstk;

	 svec = sqrt(apval_ptr1[ip].slip1*apval_ptr1[ip].slip1 +
	             apval_ptr1[ip].slip2*apval_ptr1[ip].slip2 +
	             apval_ptr1[ip].slip3*apval_ptr1[ip].slip3);

         savg = savg + svec;
	 if(svec > smax)
	    smax = svec;
         }
      }

   savg = savg/(float)(nstk*ndip);
   if((smax-savg) != (float)(0.0))
      sf = 1.0/(smax-savg);
   else
      tsfac = sf = 0.0;
   }

for(j=0;j<ndip;j++)
   {
   yy = (j + 0.5)*ddip;
   for(i=0;i<nstk;i++)
      {
      xx = (i + 0.5)*dstk - len2;
      ip = i + j*nstk;

      if(vrup < 0.0)
         {
         get_rupt(&rvmod,&htol,&dhypo,&yy,&shypo,&xx,&rayp,&rupt_rad,&rt);

	 svec = sqrt(apval_ptr1[ip].slip1*apval_ptr1[ip].slip1 +
	             apval_ptr1[ip].slip2*apval_ptr1[ip].slip2 +
	             apval_ptr1[ip].slip3*apval_ptr1[ip].slip3);

         rt = rt + sf*tsfac*(svec - savg);
         }
      else
         rt = inv_vrup*sqrt((xx-shypo)*(xx-shypo) + (yy-dhypo)*(yy-dhypo));

      apval_ptr1[ip].tinit = rt;
      }
   }

write_srf(&srf1,outfile,outbin);
}
