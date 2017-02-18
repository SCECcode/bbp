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
float area, tzero;
int i, ip;
char infile[256], outfile[256];
char str[1024];
char stype[32];

struct generic_slip gslip;
struct slippars *spar;
struct standrupformat srf;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float dt;
int it, nt;

int outbin = 0;
int plane_header = 0;

float risetimedep = 0.0;
float risetimefac = 1.0;

sprintf(srf.version,"1.0");
sprintf(stype,"brune");

nt = NTMAX;

prect_ptr = &srf.srf_prect;
prect_ptr->nseg = 1;
prect_ptr->prectseg = (struct srf_prectsegments *)check_malloc(prect_ptr->nseg*sizeof(struct srf_prectsegments));

prseg_ptr = prect_ptr->prectseg;

setpar(ac,av);

getpar("version","s",srf.version);

mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("outbin","d",&outbin);

mstpar("dt","f",&dt);
getpar("nt","d",&nt);
getpar("stype","s",stype);

getpar("risetimefac","f",&risetimefac);
getpar("risetimedep","f",&risetimedep);

getpar("plane_header","d",&plane_header);
if(plane_header)
   {
   sprintf(srf.type,"PLANE");

   mstpar("elon","f",&prseg_ptr[0].elon);
   mstpar("elat","f",&prseg_ptr[0].elat);
   mstpar("nstk","d",&prseg_ptr[0].nstk);
   mstpar("ndip","d",&prseg_ptr[0].ndip);
   mstpar("flen","f",&prseg_ptr[0].flen);
   mstpar("fwid","f",&prseg_ptr[0].fwid);

   prseg_ptr[0].dlen = prseg_ptr[0].flen/prseg_ptr[0].nstk;
   prseg_ptr[0].dwid = prseg_ptr[0].fwid/prseg_ptr[0].ndip;

   mstpar("stk","f",&prseg_ptr[0].stk);
   mstpar("dip","f",&prseg_ptr[0].dip);
   mstpar("dtop","f",&prseg_ptr[0].dtop);
   mstpar("shyp","f",&prseg_ptr[0].shyp);
   mstpar("dhyp","f",&prseg_ptr[0].dhyp);

   area = prseg_ptr[0].dlen*prseg_ptr[0].dwid*1.0e+10;
   }
else
   {
   mstpar("subfault_area","f",&area);
   area = area*1.0e+10;
   }

endpar();

fpr = fopfile(infile,"r");

fgets(str,1024,fpr);
while(strncmp(str,"#",1) == 0)
   fgets(str,1024,fpr);

sscanf(str,"%d",&gslip.np);

gslip.spar = (struct slippars *)check_malloc(gslip.np*sizeof(struct slippars));
spar = gslip.spar;

i = 0;
while(fgets(str,1024,fpr) != NULL)
   {
   sscanf(str,"%f %f %f %f %f %f %f %f",&spar[i].lon,
                                     &spar[i].lat,
                                     &spar[i].dep,
                                     &spar[i].stk,
                                     &spar[i].dip,
                                     &spar[i].rake,
                                     &spar[i].slip,
                                     &spar[i].tinit);
   i++;
   }

fclose(fpr);

srf.srf_apnts.np = gslip.np;
srf.srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((gslip.np)*sizeof(struct srf_apointvalues));
apval_ptr = srf.srf_apnts.apntvals;

for(ip=0;ip<gslip.np;ip++)
   {
   apval_ptr[ip].stf1 = (float *)check_malloc(nt*sizeof(float));
   stf = apval_ptr[ip].stf1;

   apval_ptr[ip].dt = dt;

   if(spar[ip].slip > MINSLIP)
      {
      if(strcmp(stype,"brune") == 0)
	 {
         tzero = 0.025*sqrt(spar[ip].slip);  /* assumes slip in cm */
	 if(risetimefac > 0.0 && spar[ip].dep < risetimedep)
	    tzero = tzero*risetimefac;

         apval_ptr[ip].nt1 = gen_brune_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&spar[ip].dep);
	 }
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

   apval_ptr[ip].lon = spar[ip].lon;
   apval_ptr[ip].lat = spar[ip].lat;
   apval_ptr[ip].dep = spar[ip].dep;
   apval_ptr[ip].stk = spar[ip].stk;
   apval_ptr[ip].dip = spar[ip].dip;
   apval_ptr[ip].area = area;
   apval_ptr[ip].tinit = spar[ip].tinit;
   apval_ptr[ip].rake = spar[ip].rake;
   apval_ptr[ip].slip1 = spar[ip].slip;

   apval_ptr[ip].slip2 = 0.0;
   apval_ptr[ip].nt2 = 0;
   apval_ptr[ip].stf2 = NULL;
   apval_ptr[ip].slip3 = 0.0;
   apval_ptr[ip].nt3 = 0;
   apval_ptr[ip].stf3 = NULL;
   }

write_srf(&srf,outfile,outbin);
}
