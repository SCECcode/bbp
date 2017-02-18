#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

main(int ac,char **av)
{
float cosS, sinS, cosD, sinD, *stf1, *stf2;
float bigdx, bigdy, bigdx2, bigdy2, dx, dy, clon, clat, cdep, ctinit;
float *tlon0, *tlat0, *tlon1, *tlat1;
float as, dd, dn, de, dz, rad, vr, tmin;
int it, i, ip, j, k, ig, ik;
int nfinestk, nfinedip;
char infile[256], outfile[256];

struct standrupformat srf1, srf2;
struct srf_prectsegments *prseg_ptr1, *prseg_ptr2;
struct srf_apointvalues *apval_ptr1, *apval_ptr2;

int inbin = 0;
int outbin = 0;

float dlen = -1.0;
float dwid = -1.0;

float rperd = 0.017453293;

float hlon = -1.0e+15;
float hlat = -1.0e+15;
float hdep = -1.0e+15;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("inbin","d",&inbin);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);
mstpar("nfinestk","d",&nfinestk);
mstpar("nfinedip","d",&nfinedip);
getpar("dlen","f",&dlen);
getpar("dwid","f",&dwid);
getpar("hlon","f",&hlon);
getpar("hlat","f",&hlat);
getpar("hdep","f",&hdep);
endpar();

read_srf(&srf1,infile,inbin);

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
      prseg_ptr2[ig].nstk = prseg_ptr1[ig].nstk*nfinestk;
      prseg_ptr2[ig].ndip = prseg_ptr1[ig].ndip*nfinedip;
      prseg_ptr2[ig].flen = prseg_ptr1[ig].flen;
      prseg_ptr2[ig].fwid = prseg_ptr1[ig].fwid;
      prseg_ptr2[ig].stk  = prseg_ptr1[ig].stk;
      prseg_ptr2[ig].dip  = prseg_ptr1[ig].dip;
      prseg_ptr2[ig].dtop = prseg_ptr1[ig].dtop;
      prseg_ptr2[ig].shyp = prseg_ptr1[ig].shyp;
      prseg_ptr2[ig].dhyp = prseg_ptr1[ig].dhyp;
      }
   }

srf2.srf_apnts.np = srf1.srf_apnts.np*nfinestk*nfinedip;
srf2.srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2.srf_apnts.np)*sizeof(struct srf_apointvalues));

apval_ptr1 = srf1.srf_apnts.apntvals;
apval_ptr2 = srf2.srf_apnts.apntvals;

tmin = 1.0e+15;
for(i=0;i<srf1.srf_apnts.np;i++)
   {
   if(apval_ptr1[i].tinit < tmin)
      {
      ik = i;
      tmin = apval_ptr1[i].tinit;
      }
   }

if(hlon < -1.0e+14)
   hlon = apval_ptr1[ik].lon;
if(hlat < -1.0e+14)
   hlat = apval_ptr1[ik].lat;
if(hdep < -1.0e+14)
   hdep = apval_ptr1[ik].dep;

for(i=0;i<srf1.srf_apnts.np;i++)
   {
   clon = apval_ptr1[i].lon;
   clat = apval_ptr1[i].lat;
   cdep = apval_ptr1[i].dep;
   ctinit = apval_ptr1[i].tinit;

   cosS = cos(rperd*apval_ptr1[i].stk);
   sinS = sin(rperd*apval_ptr1[i].stk);
   cosD = cos(rperd*apval_ptr1[i].dip);
   sinD = sin(rperd*apval_ptr1[i].dip);

   bigdx = dlen;
   if(bigdx < 0.0)
      bigdx = 1.0e-05*sqrt(apval_ptr1[i].area);
   
   bigdx2 = 0.5*bigdx;
   dx = bigdx/nfinestk;

   bigdy = dwid;
   if(bigdy < 0.0)
      bigdy = 1.0e-05*sqrt(apval_ptr1[i].area);
   
   bigdy2 = 0.5*bigdy;
   dy = bigdy/nfinedip;

   tlon0 = &(apval_ptr1[i].lon);
   tlat0 = &(apval_ptr1[i].lat);
   set_ne(&hlon,&hlat,tlon0,tlat0,&dn,&de);

   dz = apval_ptr1[i].dep - hdep;

   rad = sqrt(dn*dn + de*de + dz*dz);

   vr = 3.0;
   if(rad > 0.0 && ctinit > 0.0)
      vr = rad/ctinit;

   for(j=0;j<nfinedip;j++)
      {
      for(k=0;k<nfinestk;k++)
         {
	 ip = k + j*nfinestk + i*nfinestk*nfinedip;

	 as = (k+0.5)*dx - bigdx2;
	 dd = (j+0.5)*dy - bigdy2;

	 dn = as*cosS - dd*cosD*sinS;
	 de = as*sinS + dd*cosD*cosS;
	 dz = dd*sinD;

	 tlon0 = &(apval_ptr1[i].lon);
	 tlat0 = &(apval_ptr1[i].lat);
	 tlon1 = &(apval_ptr2[ip].lon);
	 tlat1 = &(apval_ptr2[ip].lat);
         set_ll(tlon0,tlat0,tlon1,tlat1,&dn,&de);

         apval_ptr2[ip].dep   = apval_ptr1[i].dep + dz;

         apval_ptr2[ip].stk   = apval_ptr1[i].stk;
         apval_ptr2[ip].dip   = apval_ptr1[i].dip;
         apval_ptr2[ip].area  = apval_ptr1[i].area/(nfinestk*nfinedip);

	 tlon0 = &(apval_ptr2[ip].lon);
	 tlat0 = &(apval_ptr2[ip].lat);
         set_ne(&hlon,&hlat,tlon0,tlat0,&dn,&de);

	 dz = apval_ptr2[ip].dep - hdep;

	 rad = sqrt(dn*dn + de*de + dz*dz);
	 apval_ptr2[ip].tinit = rad/vr;

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
         }
      }
   }

write_srf(&srf2,outfile,outbin);
}
