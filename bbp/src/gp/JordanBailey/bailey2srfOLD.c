#include "include.h"
#include "structure.h"
#include "function.h"

#define NTMAX 10000

main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float *stf;
float x0, y0, z0, dd, vslip, rake;
float len2, ds0, dd0, dsf, ddf;
float area, sn, se, cosA, sinA, slon, slat;
int ig, i, j, k, l, ip, kp, ntot, ntall;
char filelist[256], *infile[20], outfile[256];
char string[1024], filebuf[20*256];

struct rob rrm;
struct standrupformat srf;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float rupvel = -1.0;
float shal_vrup = 1.0;
float xtsfac;
float tsfac = 0.0;
float rvfrac, dt, rt;
float htol = 0.1;
double rayp, rupt_rad;
int it, nt;
struct velmodel vmod;
char modfile[128];

float cosD, sinD, arg;
double rperd = 0.017453293;

int nfinestk = 1;
int nfinedip = 1;

sprintf(srf.version,"1.0-alpha");

nt = NTMAX;

setpar(ac,av);

getpar("version","s",srf.version);

mstpar("filelist","s",filelist);
mstpar("outfile","s",outfile);

mstpar("dt","f",&dt);
getpar("nt","d",&nt);

getpar("nfinestk","d",&nfinestk);
getpar("nfinedip","d",&nfinedip);

getpar("tsfac","f",&tsfac);
getpar("rupvel","f",&rupvel);
if(rupvel < 0.0)
   {
   mstpar("modfile","s",modfile);
   mstpar("rvfrac","f",&rvfrac);
   getpar("shal_vrup","f",&shal_vrup);
   }

endpar();

fpr = fopfile(filelist,"r");

i = 0;
infile[0] = filebuf;
while(fscanf(fpr,"%s",infile[i]) != EOF)
   {
   i++;
   infile[i] = filebuf + i*256;
   }

fclose(fpr);

sprintf(srf.type,"PLANE");

prect_ptr = &srf.srf_prect;
prect_ptr->nseg = i;
prseg_ptr = prect_ptr->prectseg;
prseg_ptr = (struct srf_prectsegments *)check_malloc(prect_ptr->nseg*sizeof(struct srf_prectsegments));

apnts_ptr = &srf.srf_apnts;
apnts_ptr->np = 0;
apval_ptr = apnts_ptr->apntvals;
apval_ptr = (struct srf_apointvalues *)check_malloc(sizeof(struct srf_apointvalues));

for(ig=0;ig<prect_ptr->nseg;ig++)
   {
   fprintf(stderr,"Reading file: %s\n",infile[ig]);
   read_rob(&rrm,infile[ig],&tsfac);

   prseg_ptr[ig].elon = rrm.elon;
   prseg_ptr[ig].elat = rrm.elat;
   prseg_ptr[ig].nstk = rrm.nstk*nfinestk;
   prseg_ptr[ig].ndip = rrm.ndip*nfinedip;
   prseg_ptr[ig].flen = rrm.flen;
   prseg_ptr[ig].fwid = rrm.fwid;
   prseg_ptr[ig].dlen = rrm.flen/prseg_ptr[ig].nstk;
   prseg_ptr[ig].dwid = rrm.fwid/prseg_ptr[ig].ndip;
   prseg_ptr[ig].stk = rrm.stk;
   prseg_ptr[ig].dip = rrm.dip;
   prseg_ptr[ig].dtop = rrm.dtop;
   prseg_ptr[ig].shyp = rrm.shyp;
   prseg_ptr[ig].dhyp = rrm.dhyp;

   ntot = (prseg_ptr[ig].nstk)*(prseg_ptr[ig].ndip);
   area = prseg_ptr[ig].dlen*prseg_ptr[ig].dwid*1.0e+10;

   apval_ptr = (struct srf_apointvalues *)check_realloc(apval_ptr,(apnts_ptr->np+ntot)*sizeof(struct srf_apointvalues));

   if(rupvel < 0.0)
      conv2vrup(modfile,&vmod,&rrm.dip,&rrm.dtop,&rrm.fwid,&rvfrac,&shal_vrup);

   ds0 = rrm.flen/rrm.nstk;
   dd0 = rrm.fwid/rrm.ndip;

   dsf = ds0/nfinestk;
   ddf = dd0/nfinedip;

   len2 = 0.5*rrm.flen;

   arg = rrm.dip*rperd;
   cosD = cos(arg);
   sinD = sin(arg);

   arg = rrm.stk*rperd;
   cosA = cos(arg);
   sinA = sin(arg);

   for(j=0;j<rrm.ndip;j++)
      {
      for(i=0;i<rrm.nstk;i++)
         {
	 kp = i + j*rrm.nstk;

         for(l=0;l<nfinedip;l++)
            {
            for(k=0;k<nfinestk;k++)
               {
	       ip = k + i*nfinestk + l*prseg_ptr[ig].nstk + j*prseg_ptr[ig].nstk*nfinedip + apnts_ptr->np;

               dd = j*dd0 + (l+0.5)*ddf;

               x0 = i*ds0 + (k+0.5)*dsf - len2;
               y0 = dd*cosD;
               z0 = rrm.dtop + dd*sinD;

	       apval_ptr[ip].stf = (float *)check_malloc(nt*sizeof(float));
	       stf = apval_ptr[ip].stf;

	       apval_ptr[ip].dt = dt;
               apval_ptr[ip].nt = gen_rob_stf(&rrm,i,j,stf,nt,&dt,&z0);

	       if(apval_ptr[ip].nt)
	          apval_ptr[ip].stf = (float *)check_realloc(apval_ptr[ip].stf,(apval_ptr[ip].nt)*sizeof(float));
	       else
	          free(apval_ptr[ip].stf);

               get_rrmpars(&rrm,i,j,&x0,&dd,&rt,&vslip,&rake,&xtsfac);

               if(rt < 0.0)
                  {
                  if(rupvel < 0.0)
                     get_rupt(&vmod,&htol,&rrm.dhyp,&dd,&rrm.shyp,&x0,&rayp,&rupt_rad,&rt);
                  else
                     rt = sqrt((rrm.shyp-x0)*(rrm.shyp-x0)+(rrm.dhyp-dd)*(rrm.dhyp-dd))/rupvel;
                  rt = rt + xtsfac;
                  }

               if(rt < 0.0)
                  rt = 0.0;

               se = x0*sinA + y0*cosA;
               sn = x0*cosA - y0*sinA;
	       set_ll(&rrm.elon,&rrm.elat,&slon,&slat,&sn,&se);

	       apval_ptr[ip].as = x0;
	       apval_ptr[ip].dd = dd;
	       apval_ptr[ip].lon = slon;
	       apval_ptr[ip].lat = slat;
	       apval_ptr[ip].dep = z0;
	       apval_ptr[ip].mu = -99;
	       apval_ptr[ip].stk = prseg_ptr[ig].stk;
	       apval_ptr[ip].dip = prseg_ptr[ig].dip;
	       apval_ptr[ip].rake = rrm.rake[kp];
	       apval_ptr[ip].area = area;
	       apval_ptr[ip].slip = rrm.slip[kp];
	       apval_ptr[ip].tinit = rt;
               }
            }
         }
      }

   apnts_ptr->np = apnts_ptr->np + ntot;
   }

fpw = fopfile(outfile,"w");

fprintf(fpw,"%s\n",srf.version);

fprintf(fpw,"%s %d\n",srf.type,prect_ptr->nseg);
for(ig=0;ig<prect_ptr->nseg;ig++)
   {
   fprintf(fpw,"%10.4f %10.4f %5d %5d %8.2f %8.2f\n",prseg_ptr[ig].elon,
                                                     prseg_ptr[ig].elat,
						     prseg_ptr[ig].nstk,
						     prseg_ptr[ig].ndip,
						     prseg_ptr[ig].flen,
						     prseg_ptr[ig].fwid);
   fprintf(fpw,"%4.0f %4.0f %8.2f %8.2f %8.2f\n",prseg_ptr[ig].stk,
                                                 prseg_ptr[ig].dip,
						 prseg_ptr[ig].dtop,
						 prseg_ptr[ig].shyp,
						 prseg_ptr[ig].dhyp);
   }

fprintf(fpw,"POINTS %d\n",apnts_ptr->np);
for(i=0;i<apnts_ptr->np;i++)
   {
   fprintf(fpw,"%10.4f %10.4f %10.4f %10.4f %10.4f %13.5e\n",
					      apval_ptr[i].as,
					      apval_ptr[i].dd,
					      apval_ptr[i].lon,
					      apval_ptr[i].lat,
					      apval_ptr[i].dep,
					      apval_ptr[i].mu);
   fprintf(fpw,"%4.0f %4.0f %4.0f %13.5e %8.2f %10.4f %6d %13.5e\n",
					      apval_ptr[i].stk,
					      apval_ptr[i].dip,
					      apval_ptr[i].rake,
					      apval_ptr[i].area,
					      apval_ptr[i].slip,
					      apval_ptr[i].tinit,
					      apval_ptr[i].nt,
					      apval_ptr[i].dt);
	       
   stf = apval_ptr[i].stf;
   for(it=0;it<apval_ptr[i].nt;it++)
      fprintf(fpw,"%13.5e\n",stf[it]);
   }

fclose(fpw);
}
