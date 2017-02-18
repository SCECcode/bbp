#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

int main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpwsv, *fpwrt, *fpwtr;
struct gfheader gfhead[4];

struct sgtmaster sgtmast;
struct geoprojection geop;

float maxgft;
struct gfparam gfpar;
float *gf, *gfmech;
int kg, ig;

float kperd_n, kperd_e;
float elon = -118.0;
float elat = 34.0;
float slat, slon, snorth, seast, xs, ys;
float blon, blat, bdep, dn_avg, de_avg;
double e2, den, g2, lat0;

float strike = 0.0;
float dip = 90.0;
float rake = 180.0;
float dtop = 0.1;

float len, wid;
int i, j, k, l, ip, ip0;

int tshift_timedomain = 0;

int slip_go, jj;

struct beroza brm;
struct okumura orm;
struct gene grm;
struct rob rrm;
struct standrupformat srf;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

struct mechparam mechpar;
int maxmech;

int nstf;
float vslip;

int apv_off = 0;
int nseg = -1;
int inbin = 0;

float tsfac = 0.0;
float tmom = 0.0;
float moment = -1.0;

float rupvel = -1.0;
float shal_vrup = 1.0;
float htol = 0.1;
double rayp, rupt_rad;
float rvfrac, rt, *randt;
struct velmodel vmod, rvmod;

long seed = 1;
int randtime = 0;
float perc_randtime = 0.0;
float delt = 0.0;
int smooth_randt = 2;
int gaus_randt = 0;
int randmech = 0;
float deg_randstk = 0.0;
float deg_randdip = 0.0;
float deg_randrak = 0.0;
float zap = 0.0;
int nn, nrake;

int kp;
float *rwt, sum;
float randslip = 0.0;

float len2, ds0, dd0, dsf, ddf, s2;
int ntsum, maxnt, it, ntp2;
float mindt;
float x0, y0, z0, dd;
float x0c, ddc, avgvrup;
float shypo, dhypo;
int nsubstk, nsubdip;
int nfinestk = 1;
int nfinedip = 1;
int ntout = -99;

float *stf, *seis, *subseis, *se, *sn, *sv;
float cosS, sinS, cosA, sinA;
float scale, arg, cosD, sinD;
float xstr, xdip, xrak;
float area, sfac;
float trise;

float azi, rng, deast, dnorth;
int ncomp = 3;

float *space;
float dtout = -1.0;

int fdw, fdw_n, fdw_e, fdw_v;
char gfpath[2048], gfname[256];
char rtimesfile[2048], modfile[2048], outfile[2048];
char slipfile[2048], rupmodfile[2048], outdir[2048], stat[64], sname[64];
char rupmodtype[2048], trisefile[2048];
char string[2048];

int use_closest_gf = 1;
float min_taper_range = 250.0;
float max_taper_range = 300.0;

int write_ruptimes = 0;
int write_slipvals = 0;
int write_risetime = 0;

int hkt_format = 0;

double rperd = 0.017453293;
float normf = 1.0e+10;  /* km^2 -> cm^2 */

float targetslip = 1.0;  /* slip in cm on each subfault */
float slip_conv = 1.0;  /* input slip in cm on each subfault */

float half = 0.5;
float two = 2.0;

int latloncoords = 0;
int geoproj = 1;

float tstart = 0.0;

rtimesfile[0] = '\0';
slipfile[0] = '\0';
trisefile[0] = '\0';
sname[0] = '\0';
sprintf(rupmodtype,"NULL");

sprintf(gfpar.gftype,"fk");

setpar(ac, av);
getpar("latloncoords","d",&latloncoords);
getpar("geoproj","d",&geoproj);

if(latloncoords == 1)
   {
   getpar("elat","f",&elat);
   getpar("elon","f",&elon);
   mstpar("slat","f",&slat);
   mstpar("slon","f",&slon);
   }
else
   {
   mstpar("snorth","f",&snorth);
   mstpar("seast","f",&seast);
   }

getpar("dtop","f",&dtop);
getpar("strike","f",&strike);
getpar("dip","f",&dip);
getpar("rake","f",&rake);

getpar("rupmodtype","s",rupmodtype);

if(strcmp(rupmodtype,"BEROZA") == 0)
   {
   brm.inc_stk = 1;
   brm.inc_dip = 1;
   brm.generic_risetime = -1.0;
   brm.robstf = 0;

   mstpar("rupmodfile","s",rupmodfile);
   mstpar("npstk","d",&brm.npstk);
   mstpar("npdip","d",&brm.npdip);
   mstpar("inc_stk","d",&brm.inc_stk);
   mstpar("inc_dip","d",&brm.inc_dip);

   mstpar("len","f",&len);
   mstpar("wid","f",&wid);

   getpar("robstf","d",&brm.robstf);
   getpar("generic_risetime","f",&brm.generic_risetime);
   if(brm.robstf == 0 && brm.generic_risetime > 0.0)
      {
      mstpar("generic_pulsedur","f",&brm.generic_pulsedur);
      mstpar("generic_t2","f",&brm.generic_t2);
      }

   getpar("slip_conv","f",&slip_conv);

   mstpar("outdir","s",outdir);
   mstpar("stat","s",stat);
   }
else if(strcmp(rupmodtype,"OKUMURA") == 0)
   {
   mstpar("rupmodfile","s",rupmodfile);

   getpar("slip_conv","f",&slip_conv);

   mstpar("outdir","s",outdir);
   mstpar("stat","s",stat);
   }
else if(strcmp(rupmodtype,"GENE") == 0)
   {
   mstpar("rupmodfile","s",rupmodfile);

   getpar("slip_conv","f",&slip_conv);

   mstpar("outdir","s",outdir);
   mstpar("stat","s",stat);
   }
else if(strcmp(rupmodtype,"ROB") == 0)
   {
   mstpar("rupmodfile","s",rupmodfile);

   getpar("slip_conv","f",&slip_conv);

   mstpar("outdir","s",outdir);
   mstpar("stat","s",stat);

   mstpar("shypo","f",&shypo);
   mstpar("dhypo","f",&dhypo);

   getpar("tsfac","f",&tsfac);

   getpar("rupvel","f",&rupvel);
   if(rupvel < 0.0)
      {
      mstpar("modfile","s",modfile);
      mstpar("rvfrac","f",&rvfrac);
      getpar("shal_vrup","f",&shal_vrup);
      }
   }
else if(strcmp(rupmodtype,"SRF") == 0)
   {
   mstpar("rupmodfile","s",rupmodfile);

   getpar("slip_conv","f",&slip_conv);
   getpar("nseg","d",&nseg);
   getpar("inbin","d",&inbin);

   mstpar("outdir","s",outdir);
   mstpar("stat","s",stat);
   }
else
   {
   mstpar("shypo","f",&shypo);
   mstpar("dhypo","f",&dhypo);

   mstpar("nsubstk","d",&nsubstk);
   mstpar("nsubdip","d",&nsubdip);

   mstpar("len","f",&len);
   mstpar("wid","f",&wid);

   getpar("rupvel","f",&rupvel);
   if(rupvel < 0.0)
      {
      mstpar("modfile","s",modfile);
      mstpar("rvfrac","f",&rvfrac);
      getpar("shal_vrup","f",&shal_vrup);
      }

   getpar("targetslip","f",&targetslip);

   mstpar("outfile","s",outfile);
   getpar("hkt_format","d",&hkt_format);
   }

getpar("nfinestk","d",&nfinestk);
getpar("nfinedip","d",&nfinedip);

mstpar("gftype","s",gfpar.gftype);

if((strncmp(gfpar.gftype,"fk",2) == 0) || (strncmp(gfpar.gftype,"FK",2) == 0))
   {
   gfpar.flag3d = 0;
   gfpar.nc = 8;
   mstpar("gflocs","s",gfpar.gflocs);
   mstpar("gftimes","s",gfpar.gftimes);

   gfpar.swap_flag = 0;
   getpar("gf_swap_bytes","d",&gfpar.swap_flag);
   }
else if((strncmp(gfpar.gftype,"3d",2) == 0) || (strncmp(gfpar.gftype,"3D",2) == 0))
   {
   gfpar.flag3d = 1;
   gfpar.nc = 18;
   mstpar("gflocs","s",gfpar.gflocs);
   mstpar("gfrange_tolerance","f",&gfpar.rtol);

   gfpar.use_depdir = 1;
   getpar("use_gf3d_depdir","d",&gfpar.use_depdir);
   }
else
   {
   fprintf(stderr,"gftype= %s invalid option, exiting...\n",gfpar.gftype);
   exit(-1);
   }

mstpar("gfpath","s",gfpath);
mstpar("gfname","s",gfname);

getpar("use_closest_gf","d",&use_closest_gf);
getpar("min_taper_range","f",&min_taper_range);
getpar("max_taper_range","f",&max_taper_range);

mstpar("maxnt","d",&maxnt);
mstpar("mindt","f",&mindt);
getpar("ntout","d",&ntout);
getpar("dtout","f",&dtout);

getpar("moment","f",&moment);
getpar("tstart","f",&tstart);
getpar("rtimesfile","s",rtimesfile);
getpar("slipfile","s",slipfile);
getpar("trisefile","s",trisefile);

getpar("seed","d",&seed);
getpar("randtime","d",&randtime);
if(randtime >= 1)
   mstpar("perc_randtime","f",&perc_randtime);
if(randtime >= 2)
   getpar("delt","f",&delt);
getpar("smooth_randt","d",&smooth_randt);
getpar("gaus_randt","d",&gaus_randt);

getpar("randslip","f",&randslip);

getpar("randmech","d",&randmech);
if(randmech)
   {
   mstpar("deg_randstk","f",&deg_randstk);
   mstpar("deg_randdip","f",&deg_randdip);
   mstpar("deg_randrak","f",&deg_randrak);
   }

getpar("tshift_timedomain","d",&tshift_timedomain);

getpar("sname","s",sname);

endpar();

fprintf(stderr,"type= %s\n",rupmodtype);

maxmech = 1;
mechpar.nmech = 1;
mechpar.flag[0] = U1FLAG;
mechpar.flag[1] = 0;
mechpar.flag[2] = 0;

if(strcmp(rupmodtype,"BEROZA") == 0)
   {
   len2 = 0.5*len;

   read_beroza(&brm,rupmodfile,&len2);

   nsubstk = (brm.npstk) - 1;
   nsubdip = (brm.npdip) - 1;

   targetslip = slip_conv;
   }
else if(strcmp(rupmodtype,"OKUMURA") == 0)
   {
   read_okumura(&orm,rupmodfile,&len2);

   nsubstk = orm.nstk;
   nsubdip = orm.ndip;

   len = orm.flen;
   wid = orm.fwid;

   targetslip = slip_conv;
   }
else if(strcmp(rupmodtype,"GENE") == 0)
   {
   read_gene(&grm,rupmodfile,&len2);

   nsubstk = grm.nstk;
   nsubdip = grm.ndip;

   len = grm.flen;
   wid = grm.fwid;

   targetslip = slip_conv;
   }
else if(strcmp(rupmodtype,"ROB") == 0)
   {
   read_rob(&rrm,rupmodfile,&tsfac);

/*   07/15/04
     For now, just use the getpar values, eventually we should modify
     in order to use the values read in from the slipmodel
*/

   rrm.elon = elon;
   rrm.elat = elat;
   rrm.stk = strike;
   rrm.dip = dip;
   rrm.dtop = dtop;
   rrm.shyp = shypo;
   rrm.dhyp = dhypo;

   nsubstk = rrm.nstk;
   nsubdip = rrm.ndip;

   len = rrm.flen;
   wid = rrm.fwid;

   len2 = 0.5*len;

   if(rupvel < 0.0)
      {
      read_velmodel(modfile,&vmod);
      conv2vrup(&vmod,&rvmod,&dip,&dtop,&wid,&rvfrac,&shal_vrup);
      }

   targetslip = slip_conv;
   }
else if(strcmp(rupmodtype,"SRF") == 0)
   {
   maxmech = 3;
   nfinestk = 1;
   nfinedip = 1;

   read_srf(&srf,rupmodfile,inbin);
   prect_ptr = &srf.srf_prect;
   prseg_ptr = prect_ptr->prectseg;
   apnts_ptr = &srf.srf_apnts;
   apval_ptr = apnts_ptr->apntvals;

/*
03/03/06

     Default is to use all POINTS specified in SRF file.  This is
     done with nseg=-1.  Values for 'shypo', 'dhypo', 'len', 'wid'
     are not important.

     -or-

     Only use one segment from standard rupture model format;
     specified with 'nseg'.
*/

   if(nseg < 0)  /* do all POINTS */
      {
      shypo = 0.0;
      dhypo = 0.0;

      nsubstk = srf.srf_apnts.np;
      nsubdip = 1;

      len = 10.0;
      wid = 10.0;

      apv_off = 0;
      }
   else
      {
      elon = prseg_ptr[nseg].elon;
      elat = prseg_ptr[nseg].elat;
      strike = prseg_ptr[nseg].stk;
      dip = prseg_ptr[nseg].dip;
      dtop = prseg_ptr[nseg].dtop;
      shypo = prseg_ptr[nseg].shyp;
      dhypo = prseg_ptr[nseg].dhyp;

      nsubstk = prseg_ptr[nseg].nstk;
      nsubdip = prseg_ptr[nseg].ndip;

      len = prseg_ptr[nseg].flen;
      wid = prseg_ptr[nseg].fwid;

   /* reset POINTS pointer to appropriate segment */

      apv_off = 0;
      for(i=0;i<nseg;i++)
         apv_off = apv_off + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

      apval_ptr = apval_ptr + apv_off;
      }

   len2 = 0.5*len;
   targetslip = slip_conv;
   }
else
   {
   len2 = 0.5*len;

   if(rupvel < 0.0)
      {
      read_velmodel(modfile,&vmod);
      conv2vrup(&vmod,&rvmod,&dip,&dtop,&wid,&rvfrac,&shal_vrup);
      }
   }

if(randtime)
   {
   fprintf(stderr,"**** Initiation time randomized\n");
   fprintf(stderr,"          slow variation= +/-%.0f percent\n",100*perc_randtime);
   fprintf(stderr,"          fast variation= +/-%g sec\n",delt);
   }
else
   {
   perc_randtime = 0.0;
   delt = 0.0;
   }

if(randmech)
   {
   fprintf(stderr,"**** strike randomized by +/-%.0f degrees\n",deg_randstk);
   fprintf(stderr,"        dip randomized by +/-%.0f degrees\n",deg_randdip);
   fprintf(stderr,"       rake randomized by +/-%.0f degrees\n",deg_randrak);
   }
else
   {
   deg_randstk = 0.0;
   deg_randdip = 0.0;
   deg_randrak = 0.0;
   }

arg = strike*rperd;
cosS = cos(arg);
sinS = sin(arg);

arg = dip*rperd;
cosD = cos(arg);
sinD = sin(arg);

get_gfpars(&gfpar);

/* calculate lat,lon to km conversions */
/*   XXX
if(latloncoords)
   set_ne(&elon,&elat,&slon,&slat,&snorth,&seast);
*/

if(latloncoords == 1 || strcmp(rupmodtype,"SRF") == 0)
   {
   sgtmast.geoproj = geoproj;
   sgtmast.modelrot = -90;
   sgtmast.xshift = -100.0;
   sgtmast.yshift = -100.0;

   sgtmast.modellon = slon;
   sgtmast.modellat = slat;

   set_geoproj(&sgtmast,&geop);

   if(strcmp(rupmodtype,"SRF") != 0)   /* everything else except SRF */
      {
      if(geoproj == 0)
         set_ne(&elon,&elat,&slon,&slat,&snorth,&seast);
      else if(geoproj == 1)
	 {
         gcproj(&xs,&ys,&elon,&elat,&geop.erad,&geop.g0,&geop.b0,geop.amat,geop.ainv,1);
	 snorth = -(geop.xshift) - xs;
	 seast = -(geop.yshift) - ys;
	 }
      }
   }

if(dtout < 0.0)
   dtout = mindt;

if(dtout < mindt)
   maxnt = (maxnt*mindt/dtout);

ntsum = 2;
while(ntsum < 4*maxnt)
   ntsum = ntsum*2;

if(ntout < 0)
   ntout = ntsum;

gf = (float *) check_malloc (4*gfpar.nc*ntsum*sizeof(float));
gfmech = (float *) check_malloc (maxmech*12*ntsum*sizeof(float));
space = (float *) check_malloc (2*ntsum*sizeof(float));

seis = (float *) check_malloc (3*ntout*sizeof(float));
subseis = (float *) check_malloc (maxmech*3*ntout*sizeof(float));
stf = (float *) check_malloc (ntout*sizeof(float));

/* Calculate subfault responses */

ds0 = len/nsubstk;
dd0 = wid/nsubdip;

dsf = ds0/nfinestk;
ddf = dd0/nfinedip;

area = (len*wid)/(nsubstk*nsubdip);

sfac = targetslip/(nfinestk*nfinedip); /* conversion factor for scaling to moment */
if(gfpar.flag3d == 0)  /* add addtnl factors to convert mu & area for 1d GFs */
   {
   if(strcmp(rupmodtype,"SRF") == 0)
      sfac = sfac*normf;
   else
      sfac = sfac*normf*normf*area;
   }

rwt = (float *) check_malloc (nfinestk*nfinedip*sizeof(float));

if(randtime)
   {
   nn = nsubstk*nsubdip*nfinestk*nfinedip;
   randt = (float *) check_malloc (nn*sizeof(float));

   rand_init(randt,&perc_randtime,&seed,nsubstk,nsubdip,nfinestk,nfinedip,smooth_randt,gaus_randt);
   }

/* open output file */
if(strcmp(rupmodtype,"NULL") == 0)
   {
   if(hkt_format == 1)
      {
      nseg = 1;
      nrake = 1;
      nn = nsubstk*nsubdip;

      sprintf(string,"%s.n",outfile);
      fdw_n = croptrfile(string);
      rite(fdw_n,&nseg,sizeof(int));
      rite(fdw_n,&nrake,sizeof(int));
      rite(fdw_n,&nn,sizeof(int));

      sprintf(string,"%s.e",outfile);
      fdw_e = croptrfile(string);
      rite(fdw_e,&nseg,sizeof(int));
      rite(fdw_e,&nrake,sizeof(int));
      rite(fdw_e,&nn,sizeof(int));

      sprintf(string,"%s.v",outfile);
      fdw_v = croptrfile(string);
      rite(fdw_v,&nseg,sizeof(int));
      rite(fdw_v,&nrake,sizeof(int));
      rite(fdw_v,&nn,sizeof(int));
      }
   else
      fdw = croptrfile(outfile);
   }

if(rtimesfile[0] != '\0')
   {
   write_ruptimes = 1;
   fpwrt = fopfile(rtimesfile,"w");
   }

if(slipfile[0] != '\0')
   {
   write_slipvals = 1;
   fpwsv = fopfile(slipfile,"w");
   }

if(trisefile[0] != '\0')
   {
   write_risetime = 1;
   fpwtr = fopfile(trisefile,"w");
   }

zapit(seis,3*ntout);

for(i=0;i<4;i++)
   {
   gfhead[i].id = -1;  /* initialize: -1 means none read yet */
   gfhead[i].ir = -1;  /* initialize: -1 means none read yet */
   }

tmom = 0.0;
for(i=0;i<nsubstk;i++)
   {
   for(j=0;j<nsubdip;j++)
      {
      sum = 0.0;
      for(l=0;l<nfinedip*nfinestk;l++)
	 {
	 rwt[l] = randslip*sfrand(&seed);
	 sum = sum + rwt[l];
         }

      sum = sum/(float)(nfinedip*nfinestk);
      for(l=0;l<nfinedip*nfinestk;l++)
	 rwt[l] = rwt[l] - sum;

      zapit(subseis,maxmech*3*ntout);

      ip0 = i + j*nsubstk;

      dn_avg = 0.0;
      de_avg = 0.0;
      bdep = 0.0;
      for(k=0;k<nfinestk;k++)
	 {
	 x0 = i*ds0 + (k+0.5)*dsf - len2;

	 for(l=0;l<nfinedip;l++)
	    {
	    dd = j*dd0 + (l+0.5)*ddf;
	    y0 = dd*cosD;
	    z0 = dtop + dd*sinD;

	    kp = l + k*nfinedip;
	    ip = kp + (j + i*nsubdip)*nfinestk*nfinedip;

            if(strcmp(rupmodtype,"BEROZA") == 0)
	       {
               get_brmpars(&brm,i,j,&x0,&dd,&rt,&vslip);
	       trise = brm.tdur[ip0];
	       }
            else if(strcmp(rupmodtype,"OKUMURA") == 0)
	       {
               get_ormpars(&orm,i,j,&x0,&dd,&rt,&vslip);
	       trise = orm.rist[ip0];
	       }
            else if(strcmp(rupmodtype,"GENE") == 0)
	       {
               get_grmpars(&grm,i,j,&x0,&dd,&rt,&vslip,&rake);
	       trise = (grm.nt[ip0]-1)*grm.tdel + grm.trise;
	       }
            else if(strcmp(rupmodtype,"ROB") == 0)
	       {
               get_rrmpars(&rrm,i,j,&x0,&dd,&rt,&vslip,&rake,&tsfac);
	       trise = rrm.trise[ip0];

	       if(rt < 0.0)
		  {
	          if(rupvel < 0.0)
	             get_rupt(&rvmod,&htol,&dhypo,&dd,&shypo,&x0,&rayp,&rupt_rad,&rt);
	          else
	             rt = sqrt((shypo-x0)*(shypo-x0)+(dhypo-dd)*(dhypo-dd))/rupvel;
		  rt = rt + tsfac;
		  }

               if(rt < 0.0)
                  rt = 0.0;
	       }
            else if(strcmp(rupmodtype,"SRF") == 0)
	       {
               get_srfpars(&srf,apv_off,ip0,&rt,&vslip,&strike,&dip,&rake,&mechpar);
	       trise = apval_ptr[ip0].dt*apval_ptr[ip0].nt1;

/*
   For case when nfinestk,nfinedip > 1 =>
   calculate avg. Vr based on subfault center, then re-estimate Tinit 
   when nfinestk = nfinedip = 1, x0c=x0, ddc=dd.
*/
	       if(nfinestk != 1 || nfinedip != 1)
	          {
	          x0c = (i+0.5)*ds0 - len2;
	          ddc = (j+0.5)*dd0;
	          avgvrup = sqrt((shypo-x0c)*(shypo-x0c)+(dhypo-ddc)*(dhypo-ddc))/rt;
	          rt = sqrt((shypo-x0)*(shypo-x0)+(dhypo-dd)*(dhypo-dd))/avgvrup;
		  }
	       }
            else
               {
               vslip = 1.0;

	       if(rupvel < 0.0)
	          get_rupt(&rvmod,&htol,&dhypo,&dd,&shypo,&x0,&rayp,&rupt_rad,&rt);
	       else
	          rt = sqrt((shypo-x0)*(shypo-x0)+(dhypo-dd)*(dhypo-dd))/rupvel;
               }

	    if(randtime)
	       rt = rt*(1.0 + randt[ip]);
	    if(randtime == 2)
	       {
	       rt = rt + delt*sfrand(&seed);
	       if(rt < 0.0)
		  rt = 0.0;
               }

	    if(write_ruptimes == 1)
	       fprintf(fpwrt,"%13.5e %13.5e %13.5e\n",x0+len2,dd,rt);

            vslip = (1.0 + rwt[kp])*vslip;

            if(write_slipvals == 1)
               fprintf(fpwsv,"%13.5e %13.5e %13.5e\n",x0+len2,dd,slip_conv*vslip);

	    if(write_risetime == 1)
	       fprintf(fpwtr,"%13.5e %13.5e %13.5e\n",x0+len2,dd,trise);

            if(strcmp(rupmodtype,"SRF") == 0)
	       get_ard_srf(&srf,apv_off,ip0,&azi,&rng,&z0,&deast,&dnorth,&slon,&slat,&geop);
            else
	       get_radazi(&azi,&rng,&deast,&dnorth,&x0,&y0,&cosS,&sinS,&seast,&snorth);

/* RWG 20140314: Comment out following lines, never look at this output anymore
	    fprintf(stderr,"i=%3d j=%3d k=%3d l=%3d ",i,j,k,l);
	    fprintf(stderr," s=%7.2f d=%7.2f",x0,dd);
	    fprintf(stderr," dn=%10.5f de=%10.5f",dnorth,deast);
	    fprintf(stderr," a=%7.2f r=%7.2f\n",azi,rng);
*/

            dn_avg = dn_avg + dnorth;
            de_avg = de_avg + deast;
            bdep = bdep + z0;

	    if(vslip > (float)(0.0))
	       {
	       find_4gf(gfpar,gfhead,&rng,&z0,&deast,&dnorth);

	       if(use_closest_gf == 1 && rng > min_taper_range)
	          reset_wt_4gf(gfhead,&rng,&min_taper_range,&max_taper_range);

               read_4gf(gfpath,gfname,gf,ntsum,gfhead,gfpar,&maxgft,&maxnt,&dtout,space);

	       if(randmech)
	          {
	          mechpar.stk = strike + deg_randstk*sfrand(&seed);
	          mechpar.dip = dip + deg_randdip*sfrand(&seed);
	          mechpar.rak = rake + deg_randrak*sfrand(&seed);
	          }
	       else
	          {
	          mechpar.stk = strike;
	          mechpar.dip = dip;
	          mechpar.rak = rake;
	          }

               if(strcmp(rupmodtype,"SRF") == 0)
	          scale = sfac*apval_ptr[ip0].area;
	       else
	          scale = sfac;

	       mech_4gf(gfmech,gf,gfhead,gfpar,ntsum,mechpar,&azi,&scale);

/* scale now contains the moment released by this point source */
	       tmom = tmom + vslip*scale;

	       sum_4gf(subseis,ntout,gfmech,gfhead,ntsum,maxnt,&rt,&maxgft,&tstart,tshift_timedomain,mechpar);
	       }
	    }
         }

      z0 = dtop + (j+0.5)*dd0*sinD;

      if(strcmp(rupmodtype,"BEROZA") == 0)
	 beroza_stf(&brm,i,j,seis,subseis,stf,ntout,&dtout,&z0);
      else if(strcmp(rupmodtype,"OKUMURA") == 0)
	 okumura_stf(&orm,i,j,seis,subseis,stf,ntout,&dtout);
      else if(strcmp(rupmodtype,"GENE") == 0)
	 gene_stf(&grm,i,j,seis,subseis,stf,ntout,&dtout);
      else if(strcmp(rupmodtype,"ROB") == 0)
	 rob_stf(&rrm,i,j,seis,subseis,stf,ntout,&dtout,&z0);
      else if(strcmp(rupmodtype,"SRF") == 0)
	 srf_stf(&srf,apv_off,ip0,seis,subseis,stf,ntout,&dtout,mechpar,space);
      else
         {
         sv = subseis;
         sn = subseis + ntout;
         se = subseis + 2*ntout;

	 dn_avg = dn_avg/(nfinestk*nfinedip);
	 de_avg = de_avg/(nfinestk*nfinedip);
	 bdep = bdep/(nfinestk*nfinedip);
         if(geoproj == 0)
            {
            xs = -dn_avg;
            ys = -de_avg;
            set_ll(&slon,&slat,&blon,&blat,&xs,&ys);
            }
         else if(geoproj == 1)
	    {
	    xs = -(geop.xshift) - dn_avg;
	    ys = -(geop.yshift) - de_avg;
            gcproj(&xs,&ys,&blon,&blat,&geop.erad,&geop.g0,&geop.b0,geop.amat,geop.ainv,0);
	    }

         fprintf(stderr,"*** blon=%12.5f blat=%12.5f bdep=%12.5f\n",blon,blat,bdep);

	 if(hkt_format == 1)
	    {
            rite(fdw_n,&ntout,sizeof(int));
            rite(fdw_n,&dtout,sizeof(float));
            rite(fdw_n,&tstart,sizeof(float));
            rite(fdw_n,&blat,sizeof(float));
            rite(fdw_n,&blon,sizeof(float));
            rite(fdw_n,&bdep,sizeof(float));
            rite(fdw_n,&strike,sizeof(float));
            rite(fdw_n,&dip,sizeof(float));
            rite(fdw_n,&rake,sizeof(float));
            rite(fdw_n,sn,ntout*sizeof(float));

            rite(fdw_e,&ntout,sizeof(int));
            rite(fdw_e,&dtout,sizeof(float));
            rite(fdw_e,&tstart,sizeof(float));
            rite(fdw_e,&blat,sizeof(float));
            rite(fdw_e,&blon,sizeof(float));
            rite(fdw_e,&bdep,sizeof(float));
            rite(fdw_e,&strike,sizeof(float));
            rite(fdw_e,&dip,sizeof(float));
            rite(fdw_e,&rake,sizeof(float));
            rite(fdw_e,se,ntout*sizeof(float));

            rite(fdw_v,&ntout,sizeof(int));
            rite(fdw_v,&dtout,sizeof(float));
            rite(fdw_v,&tstart,sizeof(float));
            rite(fdw_v,&blat,sizeof(float));
            rite(fdw_v,&blon,sizeof(float));
            rite(fdw_v,&bdep,sizeof(float));
            rite(fdw_v,&strike,sizeof(float));
            rite(fdw_v,&dip,sizeof(float));
            rite(fdw_v,&rake,sizeof(float));
            rite(fdw_v,sv,ntout*sizeof(float));
	    }
	 else
	    {
            fortran_rite(fdw,1,&ncomp,sizeof(int));

            fortran_rite(fdw,2,&rng,sizeof(float),&tstart,sizeof(float));
            fortran_rite(fdw,2,&ntout,sizeof(int),&dtout,sizeof(float));
            fortran_rite(fdw,1,sn,ntout*sizeof(float));

            fortran_rite(fdw,2,&rng,sizeof(float),&tstart,sizeof(float));
            fortran_rite(fdw,2,&ntout,sizeof(int),&dtout,sizeof(float));
            fortran_rite(fdw,1,se,ntout*sizeof(float));

            fortran_rite(fdw,2,&rng,sizeof(float),&tstart,sizeof(float));
            fortran_rite(fdw,2,&ntout,sizeof(int),&dtout,sizeof(float));
            fortran_rite(fdw,1,sv,ntout*sizeof(float));
	    }
         }
      }
   }

if(strcmp(rupmodtype,"NULL") == 0)
   {
   if(hkt_format == 1)
      {
      close(fdw_n);
      close(fdw_e);
      close(fdw_v);
      }
   else
      close(fdw);
   }
else
   {
   sv = seis;
   sn = seis + ntout;
   se = seis + 2*ntout;

   if(sname[0] == '\0')
      {
      strncpy(sname,stat,7);
      sname[7] = '\0';
      }

   fprintf(stderr,"Total summed moment= %13.5e\n",tmom);

   if(moment > 0.0)
      {
      sfac = moment/tmom;
      for(it=0;it<ntout;it++)
         {
	 sn[it] = sfac*sn[it];
	 se[it] = sfac*se[it];
	 sv[it] = sfac*sv[it];
	 }

      fprintf(stderr,"      Scaled moment= %13.5e\n",moment);
      }

   write_seis(outdir,stat,sname,"000",sn,&dtout,ntout,&tstart);
   write_seis(outdir,stat,sname,"090",se,&dtout,ntout,&tstart);
   write_seis(outdir,stat,sname,"ver",sv,&dtout,ntout,&tstart);
   }

if(write_ruptimes == 1)
   {
   fflush(fpwrt);
   fclose(fpwrt);
   }

if(write_slipvals == 1)
   {
   fflush(fpwsv);
   fclose(fpwsv);
   }

if(write_risetime == 1)
   {
   fflush(fpwtr);
   fclose(fpwtr);
   }
}
