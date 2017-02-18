#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

int main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float *stf;
float tzero;
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

struct slippars *spar2;
float depmin, dold, avglon, avglat, se, sn;
float hyplon, hyplat, hypdep, tmin;
int iseg, nseg, nc;

float dmin, dmax, rtfac;
float sabs, snum, sden;
float rt_scalefac = -1.0;
float rtfac_min = 0.05;

long seed = 0;
float rt_rand = -1.0;

float rperd = 0.017453293;

float dt;
int it, nt, nr;

int outbin = 0;
int plane_header = 0;

float risetimedep_range = 1.5;
float risetimedep = 6.5;
float risetimefac = 2.0;
float risetime_dipfac = 0.0;
float risetime, tau1r;
float brune_corner_freq = -1.0;
float fzero = 0.0;

sprintf(srf.version,"1.0");
sprintf(stype,"brune");

nt = NTMAX;

setpar(ac,av);

getpar("version","s",srf.version);

mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
getpar("outbin","d",&outbin);

mstpar("dt","f",&dt);
getpar("nt","d",&nt);
getpar("stype","s",stype);
if(strcmp(stype,"esg2006") == 0 || strcmp(stype,"urs") == 0 || strcmp(stype,"ucsb") == 0 || strcmp(stype,"ucsb2") == 0 || strncmp(stype,"ucsb-T",6) == 0 || strcmp(stype,"ucsb-varT1") == 0 || strncmp(stype,"cos",3) == 0 || strncmp(stype,"seki",4) == 0)
   {
   mstpar("risetime","f",&risetime);
   }

if(strcmp(stype,"brune") == 0)
   {
   getpar("brune_corner_freq","f",&brune_corner_freq);
   }

getpar("risetimefac","f",&risetimefac);
getpar("risetimedep","f",&risetimedep);
getpar("risetimedep_range","f",&risetimedep_range);
getpar("risetime_dipfac","f",&risetime_dipfac);

getpar("rt_scalefac","f",&rt_scalefac);
getpar("rtfac_min","f",&rtfac_min);

getpar("rt_rand","f",&rt_rand);
getpar("seed","d",&seed);

getpar("plane_header","d",&plane_header);

endpar();

dmin = risetimedep - risetimedep_range;
dmax = risetimedep + risetimedep_range;

fpr = fopfile(infile,"r");

fgets(str,1024,fpr);
while(strncmp(str,"#",1) == 0)
   fgets(str,1024,fpr);

sscanf(str,"%d",&gslip.np);

gslip.spar = (struct slippars *)check_malloc(gslip.np*sizeof(struct slippars));
spar = gslip.spar;

snum = 0.0;
sden = 0.0;
tmin = 1.0e+15;
i = 0;
while(fgets(str,1024,fpr) != NULL)
   {
   spar[i].trise = -1.0;
   spar[i].tau1_ratio = -1.0;

   if(risetime > 0.0)
      {
      sscanf(str,"%f %f %f %f %f %f %f %f %f %f %d",&spar[i].lon,
                                     &spar[i].lat,
                                     &spar[i].dep,
                                     &spar[i].ds,
                                     &spar[i].dw,
                                     &spar[i].stk,
                                     &spar[i].dip,
                                     &spar[i].rake,
                                     &spar[i].slip,
                                     &spar[i].tinit,
                                     &spar[i].segno);
      }
   else
      {
      nr = sscanf(str,"%f %f %f %f %f %f %f %f %f %f %d %f %f",&spar[i].lon,
                                     &spar[i].lat,
                                     &spar[i].dep,
                                     &spar[i].ds,
                                     &spar[i].dw,
                                     &spar[i].stk,
                                     &spar[i].dip,
                                     &spar[i].rake,
                                     &spar[i].slip,
                                     &spar[i].tinit,
                                     &spar[i].segno,
                                     &spar[i].trise,
                                     &spar[i].tau1_ratio);

      if(nr <= 12)
         spar[i].tau1_ratio = -1.0;
      if(nr <= 11)
         spar[i].trise = -1.0;
      }

   if(spar[i].tinit < tmin)
      {
      tmin = spar[i].tinit;
      hyplon = spar[i].lon;
      hyplat = spar[i].lat;
      hypdep = spar[i].dep;
      }

   if(spar[i].slip < 0.0)
      {
      spar[i].slip = -spar[i].slip;
      spar[i].rake = 180.0 + spar[i].rake;

      while(spar[i].rake >= 360.0)
         spar[i].rake = spar[i].rake - 360.0;

      while(spar[i].rake < 0.0)
         spar[i].rake = spar[i].rake + 360.0;
      }

   if(spar[i].dep >= dmax)
      rtfac = 1.0;
   else if(spar[i].dep < dmax && spar[i].dep > dmin)
      rtfac = 1.0 + (risetimefac - 1.0)*(dmax-(spar[i].dep))/(dmax-dmin);
   else
      rtfac = risetimefac;

   sabs = sqrt(spar[i].slip*spar[i].slip);

/* this uses slip-weighted averaging (does NOT include depth variation)
*/

   /*
   snum = snum + sabs*sqrt(sabs);
   */
   snum = snum + sabs*exp(0.5*log(sabs));
   sden = sden + sabs;

/* this uses slip-weighted averaging (includes depth variation)

   snum = snum + rtfac*sabs*sqrt(sabs);
   sden = sden + sabs;
*/

/* this uses equal weighting (whole fault average, includes depth variation)

   snum = snum + rtfac*sqrt(sabs);
   sden = sden + 1.0;
*/

   i++;
   }

fclose(fpr);

if(rt_scalefac > 0)
   {
   rt_scalefac = snum/sden;
   fprintf(stderr,"rt_scalefac= %f\n",rt_scalefac);
   }

if(plane_header)
   {
   nseg = 0;
   for(i=0;i<gslip.np;i++)
      {
      if(spar[i].segno > nseg)
         nseg = spar[i].segno;
      }
   nseg++;

   spar2 = (struct slippars *)check_malloc(gslip.np*sizeof(struct slippars));

   sprintf(srf.type,"PLANE");

   prect_ptr = &srf.srf_prect;
   prect_ptr->nseg = nseg;
   prect_ptr->prectseg = (struct srf_prectsegments *)check_malloc(prect_ptr->nseg*sizeof(struct srf_prectsegments));

   prseg_ptr = prect_ptr->prectseg;

   nc = 0;
   for(i=0;i<nseg;i++)
      {
      avglon = 0.0;
      avglat = 0.0;
      prseg_ptr[i].stk = 0.0;
      prseg_ptr[i].dip = 0.0;
      prseg_ptr[i].dlen = 0.0;
      prseg_ptr[i].dwid = 0.0;
      iseg = 0;
      dold = -1;
      prseg_ptr[i].nstk = 0;
      prseg_ptr[i].ndip = 0;

      depmin = 1.0e+15;
      for(ip=0;ip<gslip.np;ip++)
         {
	 if(spar[ip].segno == i && spar[ip].dep < depmin)
	    depmin = spar[ip].dep;
	 }

      for(ip=0;ip<gslip.np;ip++)
         {
	 if(spar[ip].segno == i)
	    {
	    spar2[nc].lon = spar[ip].lon;
	    spar2[nc].lat = spar[ip].lat;
	    spar2[nc].dep = spar[ip].dep;
	    spar2[nc].ds = spar[ip].ds;
	    spar2[nc].dw = spar[ip].dw;
	    spar2[nc].stk = spar[ip].stk;
	    spar2[nc].dip = spar[ip].dip;
	    spar2[nc].rake = spar[ip].rake;
	    spar2[nc].slip = spar[ip].slip;
	    spar2[nc].tinit = spar[ip].tinit;
	    spar2[nc].segno = spar[ip].segno;

            if(spar[ip].dep <= 1.001*depmin && spar[ip].dep >= 0.999*depmin)
	       {
	       avglon = avglon + spar[ip].lon;
	       avglat = avglat + spar[ip].lat;
	       prseg_ptr[i].nstk++;
	       }

            if(spar[ip].dep != dold)
	       {
	       dold = spar[ip].dep;
	       prseg_ptr[i].ndip++;
	       }

	    prseg_ptr[i].dlen = prseg_ptr[i].dlen + spar[ip].ds;
	    prseg_ptr[i].dwid = prseg_ptr[i].dwid + spar[ip].dw;

	    nc++;
	    iseg++;

	    if(iseg == 1)
	       {
	       prseg_ptr[i].stk = spar[ip].stk;
	       prseg_ptr[i].dip = spar[ip].dip;
	       }
	    else
	       {
	       if(spar[ip].stk > prseg_ptr[i].stk/(float)(iseg-1) + 90.0)
	          {
		  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk - 180.0;
		  prseg_ptr[i].dip = prseg_ptr[i].dip + 180.0 - spar[ip].dip;
		  }
	       else if(spar[ip].stk < prseg_ptr[i].stk/(float)(iseg-1) - 90.0)
	          {
		  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk + 180.0;
		  prseg_ptr[i].dip = prseg_ptr[i].dip + 180.0 - spar[ip].dip;
		  }
	       else
	          {
		  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk;
		  prseg_ptr[i].dip = prseg_ptr[i].dip + spar[ip].dip;
		  }
	       }
	    }
         }

      prseg_ptr[i].stk = prseg_ptr[i].stk/(float)(iseg);
      prseg_ptr[i].dip = prseg_ptr[i].dip/(float)(iseg);

      if(prseg_ptr[i].dip > 90.0)
         {
	 prseg_ptr[i].dip = 180.0 - prseg_ptr[i].dip;
	 prseg_ptr[i].stk = prseg_ptr[i].stk - 180.0;
	 }

      while(prseg_ptr[i].stk < 0.0)
         prseg_ptr[i].stk = prseg_ptr[i].stk + 360.0;
      while(prseg_ptr[i].stk >= 360.0)
         prseg_ptr[i].stk = prseg_ptr[i].stk - 360.0;

      prseg_ptr[i].dlen = prseg_ptr[i].dlen/(float)(iseg);
      prseg_ptr[i].dwid = prseg_ptr[i].dwid/(float)(iseg);

      prseg_ptr[i].flen = prseg_ptr[i].nstk*prseg_ptr[i].dlen;
      prseg_ptr[i].fwid = prseg_ptr[i].ndip*prseg_ptr[i].dwid;

      prseg_ptr[i].dtop = depmin - 0.5*prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip);
      if(prseg_ptr[i].dtop < 0.0)
         prseg_ptr[i].dtop = 0.0;

      avglon = avglon/(float)(prseg_ptr[i].nstk);
      avglat = avglat/(float)(prseg_ptr[i].nstk);
      se = -0.5*prseg_ptr[i].dwid*cos(rperd*prseg_ptr[i].dip)*cos(rperd*prseg_ptr[i].stk);
      sn = 0.5*prseg_ptr[i].dwid*cos(rperd*prseg_ptr[i].dip)*sin(rperd*prseg_ptr[i].stk);

      set_ll(&avglon,&avglat,&prseg_ptr[i].elon,&prseg_ptr[i].elat,&sn,&se);

      set_ne(&prseg_ptr[i].elon,&prseg_ptr[i].elat,&hyplon,&hyplat,&sn,&se);

      prseg_ptr[i].shyp = -999.9;
      prseg_ptr[i].dhyp = -999.9;

      prseg_ptr[i].shyp = se*sin(rperd*prseg_ptr[i].stk) + sn*cos(rperd*prseg_ptr[i].stk);
      prseg_ptr[i].dhyp = (hypdep-prseg_ptr[i].dtop)/sin(rperd*prseg_ptr[i].dip);

      fprintf(stderr,"%3d: %11.5f %11.5f %2d %2d %10.4f %10.4f %.1f %.1f %10.4f\n",i,prseg_ptr[i].elon,prseg_ptr[i].elat,prseg_ptr[i].nstk,prseg_ptr[i].ndip,prseg_ptr[i].flen,prseg_ptr[i].fwid,prseg_ptr[i].stk,prseg_ptr[i].dip,prseg_ptr[i].dtop);
      }

   for(ip=0;ip<gslip.np;ip++)
      {
      spar[ip].lon = spar2[ip].lon;
      spar[ip].lat = spar2[ip].lat;
      spar[ip].dep = spar2[ip].dep;
      spar[ip].ds = spar2[ip].ds;
      spar[ip].dw = spar2[ip].dw;
      spar[ip].stk = spar2[ip].stk;
      spar[ip].dip = spar2[ip].dip;
      spar[ip].rake = spar2[ip].rake;
      spar[ip].slip = spar2[ip].slip;
      spar[ip].tinit = spar2[ip].tinit;
      spar[ip].segno = spar2[ip].segno;
      }

   free(spar2);
   }

srf.srf_apnts.np = gslip.np;
srf.srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((gslip.np)*sizeof(struct srf_apointvalues));
apval_ptr = srf.srf_apnts.apntvals;

for(ip=0;ip<gslip.np;ip++)
   {
   apval_ptr[ip].stf1 = (float *)check_malloc(nt*sizeof(float));
   stf = apval_ptr[ip].stf1;

   apval_ptr[ip].dt = dt;

   sabs = sqrt(spar[ip].slip*spar[ip].slip);
   if(sabs > MINSLIP)
      {
      if(spar[ip].dep >= dmax)
         rtfac = 1.0*(1.0 + risetime_dipfac);
      else if(spar[ip].dep < dmax && spar[ip].dep > dmin)
	 {
         rtfac = 1.0 + (risetimefac - 1.0)*(dmax-(spar[ip].dep))/(dmax-dmin);
         rtfac = rtfac + rtfac*risetime_dipfac*((spar[ip].dep)-dmin)/(dmax-dmin);
         }
      else
         rtfac = risetimefac;

      if(rt_rand > 0.0)
         rtfac = rtfac*(1.0 + gaussian_rand(&rt_rand,&fzero,&seed));

      if(rt_scalefac > 0.0)
         {
         rtfac = rtfac*exp(0.5*log(sabs))/rt_scalefac;
	 if(rtfac < rtfac_min)
	    rtfac = rtfac_min;
         }

      if(spar[ip].trise >= 0.0)
         {
         rtfac = 1.0;
	 risetime = spar[ip].trise;
         }

      if(strcmp(stype,"brune") == 0)
	 {
	 if(brune_corner_freq > 0.0)
	    tzero = 1.0/(2.0*3.14159*brune_corner_freq);
	 else
	    {
            tzero = 0.1*exp(-1.0)*sqrt(spar[ip].slip)/(1.5); /* assume slip in cm */
            tzero = 0.1*exp(-1.0)*sqrt(spar[ip].slip)/(1.2); /* assume slip in cm */
	    }

	 tzero = tzero*rtfac;

         apval_ptr[ip].nt1 = gen_brune_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&spar[ip].dep);
	 }
      else if(strcmp(stype,"urs") == 0)
         {
         tzero = rtfac*risetime;

         apval_ptr[ip].nt1 = gen_2tri_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&spar[ip].dep);
         }
      else if(strcmp(stype,"ucsb") == 0)
         {
         tzero = rtfac*risetime;

         apval_ptr[ip].nt1 = gen_ucsb_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&spar[ip].dep);
         }
      else if(strncmp(stype,"ucsb-T",6) == 0)
         {
         tzero = rtfac*risetime;
	 tau1r = 1.0;
         if(strncmp("-T",stype+4,2) == 0)    /* rise time scaling factor given after '-T' */
            tau1r = atof(stype+6);

         apval_ptr[ip].nt1 = gen_ucsbT_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&tau1r);
         }
      else if(strcmp(stype,"ucsb2") == 0)
         {
         tzero = rtfac*risetime;

         apval_ptr[ip].nt1 = gen_ucsb2_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&spar[ip].dep);
         }
      else if(strcmp(stype,"ucsb-varT1") == 0)
         {
         tzero = rtfac*risetime;
	 tau1r = spar[ip].tau1_ratio;
	 if(spar[ip].tau1_ratio < 0.0)
	    tau1r = 0.13;

         apval_ptr[ip].nt1 = gen_ucsbvT_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&tau1r);
         }
      else if(strcmp(stype,"esg2006") == 0)
         {
         tzero = rtfac*risetime;

         apval_ptr[ip].nt1 = gen_esg2006_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&spar[ip].dep);
         }
      else if(strncmp(stype,"cos",3) == 0)
         {
         tzero = rtfac*risetime;

         apval_ptr[ip].nt1 = gen_cos_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&spar[ip].dep);
         }
      else if(strncmp(stype,"seki",4) == 0)
         {
         tzero = rtfac*risetime;

         apval_ptr[ip].nt1 = gen_seki_stf(&spar[ip].slip,&tzero,stf,nt,&dt,&spar[ip].dep);

	 spar[ip].tinit = spar[ip].tinit - 0.25*tzero; /* adjust for non-causal tanh() */
	 if(spar[ip].tinit < 0.0)
	    spar[ip].tinit = 0.0;
         }
      else if(strcmp(stype,"delta") == 0)
         {
         apval_ptr[ip].nt1 = 3;
	 stf[0] = 0.0;
	 stf[1] = spar[ip].slip/dt;
	 stf[2] = 0.0;
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
   apval_ptr[ip].area = spar[ip].ds*spar[ip].dw*1.0e+10;
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
