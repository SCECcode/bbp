#include "include.h"
#include "structure.h"
#include "../../JordanBailey/structure.h"
#include "function.h"
#include "../../JordanBailey/function.h"
#include "defs.h"

main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float *stf;
float x0, y0, z0, dd, vslip, rake;
float len2, ds0, dd0, dsf, ddf;
float area, sn, se, cosA, sinA, slon, slat;
int ig, i, j, k, l, ip, kp, ntot, ntall, nt6;
char filelist[256], *infile[20], outfile[256];
char string[1024], filebuf[20*256];
char stype[32];

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
struct velmodel vmod, rvmod;
char modfile[128];

float cosD, sinD, arg;
double rperd = 0.017453293;

int nfinestk = 1;
int nfinedip = 1;

int outbin = 0;

sprintf(srf.version,"1.0");
sprintf(stype,"rob");

nt = NTMAX;

setpar(ac,av);

getpar("version","s",srf.version);

mstpar("filelist","s",filelist);
mstpar("outfile","s",outfile);

getpar("outbin","d",&outbin);

mstpar("dt","f",&dt);
getpar("nt","d",&nt);

getpar("nfinestk","d",&nfinestk);
getpar("nfinedip","d",&nfinedip);

getpar("stype","s",stype);

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
prect_ptr->prectseg = (struct srf_prectsegments *)check_malloc(prect_ptr->nseg*sizeof(struct srf_prectsegments));
prseg_ptr = prect_ptr->prectseg;

srf.srf_apnts.np = 0;
srf.srf_apnts.apntvals = NULL;
apnts_ptr = &srf.srf_apnts;

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

   srf.srf_apnts.apntvals = (struct srf_apointvalues *)check_realloc(srf.srf_apnts.apntvals,(apnts_ptr->np+ntot)*sizeof(struct srf_apointvalues));
   apval_ptr = srf.srf_apnts.apntvals + apnts_ptr->np;

   if(rupvel < 0.0)
      {
      read_velmodel(modfile,&vmod);
      conv2vrup(&vmod,&rvmod,&rrm.dip,&rrm.dtop,&rrm.fwid,&rvfrac,&shal_vrup);
      }

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
	       ip = k + i*nfinestk + l*prseg_ptr[ig].nstk + j*prseg_ptr[ig].nstk*nfinedip;

               dd = j*dd0 + (l+0.5)*ddf;

               x0 = i*ds0 + (k+0.5)*dsf - len2;
               y0 = dd*cosD;
               z0 = rrm.dtop + dd*sinD;

	       apval_ptr[ip].stf1 = (float *)check_malloc(nt*sizeof(float));
	       stf = apval_ptr[ip].stf1;

	       apval_ptr[ip].dt = dt;

               if(rrm.slip[kp] > MINSLIP)
                  {
                  if(strcmp(stype,"rob") == 0)
                     apval_ptr[ip].nt1 = gen_rob_stf(&rrm,i,j,stf,nt,&dt,&z0);

                  if(strcmp(stype,"brune") == 0)
                     apval_ptr[ip].nt1 = gen_brune_stf(&rrm.slip[kp],&rrm.trise[kp],stf,nt,&dt,&z0);
                  }
               else
                  apval_ptr[ip].nt1 = 0;

	       if(apval_ptr[ip].nt1)
	          apval_ptr[ip].stf1 = (float *)check_realloc(apval_ptr[ip].stf1,(apval_ptr[ip].nt1)*sizeof(float));
	       else
		  {
	          free(apval_ptr[ip].stf1);
		  apval_ptr[ip].stf1 = NULL;
		  }

               get_rrmpars(&rrm,i,j,&x0,&dd,&rt,&vslip,&rake,&xtsfac);

               if(rt < 0.0)
                  {
                  if(rupvel < 0.0)
                     get_rupt(&rvmod,&htol,&rrm.dhyp,&dd,&rrm.shyp,&x0,&rayp,&rupt_rad,&rt);
                  else
                     rt = sqrt((rrm.shyp-x0)*(rrm.shyp-x0)+(rrm.dhyp-dd)*(rrm.dhyp-dd))/rupvel;
                  rt = rt + xtsfac;
                  }

               if(rt < 0.0)
                  rt = 0.0;

               se = x0*sinA + y0*cosA;
               sn = x0*cosA - y0*sinA;
	       set_ll(&rrm.elon,&rrm.elat,&slon,&slat,&sn,&se);

	       apval_ptr[ip].lon = slon;
	       apval_ptr[ip].lat = slat;
	       apval_ptr[ip].dep = z0;
	       apval_ptr[ip].stk = prseg_ptr[ig].stk;
	       apval_ptr[ip].dip = prseg_ptr[ig].dip;
	       apval_ptr[ip].area = area;
	       apval_ptr[ip].tinit = rt;
	       apval_ptr[ip].rake = rrm.rake[kp];
	       apval_ptr[ip].slip1 = rrm.slip[kp];

	       apval_ptr[ip].slip2 = 0.0;
	       apval_ptr[ip].nt2 = 0;
               apval_ptr[ip].stf2 = NULL;
	       apval_ptr[ip].slip3 = 0.0;
	       apval_ptr[ip].nt3 = 0;
               apval_ptr[ip].stf3 = NULL;
               }
            }
         }
      }

   apnts_ptr->np = apnts_ptr->np + ntot;
   }

write_srf(&srf,outfile,outbin);
}
