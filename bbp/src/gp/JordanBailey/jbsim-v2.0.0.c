#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

#define NCHAR_STAT 16

int main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpwsv, *fpwrt, *fpwtr;
struct gfheader gfhead[4];

struct sgtmaster sgtmast;
struct geoprojection geop;

float maxgft;
struct gfparam gfpar;
float *gf, *gfmech;

float slat, slon, xs, ys;
int i, ip, np;

int tshift_timedomain = 1;

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

float tmom = 0.0;
float moment = -1.0;
float sum, rt;

float zap = 0.0;

int ntsum, maxnt, it, ntp2;
float mindt;
float z0;
int ntout = -99;

float *stf, *seis, *subseis, *se, *sn, *sv;
float scale, sfac;
float azi, rng, deast, dnorth;
int ncomp = 3;

float *space;
float dtout = -1.0;

int fdw, fdw_n, fdw_e, fdw_v;
char gfpath[2048], gfname[256];
char outfile[2048];
char rupmodfile[2048], outdir[2048], stat[NCHAR_STAT], sname[NCHAR_STAT];
char rupmodtype[2048];
char string[2048];

int use_closest_gf = 1;
float min_taper_range = 0.0;
float max_taper_range = 0.0;
float tadjust = 0.0;

float normf = 1.0e+10;  /* km^2 -> cm^2 */

float half = 0.5;
float two = 2.0;

int latloncoords = 0;
int geoproj = 1;

char statfile[2048];
int nstat, istat;
int nbuf = 1000;
char *stat_buf, *stat_ptr;
float *slon_buf, *slat_buf;

float tstart = 0.0;

sprintf(rupmodtype,"SRF");
statfile[0] = '\0';

sprintf(gfpar.gftype,"fk");

setpar(ac, av);

getpar("statfile","s",statfile);
getpar("latloncoords","d",&latloncoords);
getpar("geoproj","d",&geoproj);

if(statfile[0] != '\0')
   latloncoords = 0;
else
   latloncoords = 1;

if(latloncoords == 1)
   {
   mstpar("slat","f",&slat);
   mstpar("slon","f",&slon);
   }

if(strcmp(rupmodtype,"SRF") == 0)
   {
   mstpar("rupmodfile","s",rupmodfile);

   getpar("nseg","d",&nseg);
   getpar("inbin","d",&inbin);

   mstpar("outdir","s",outdir);
   if(latloncoords == 1)
      mstpar("stat","s",stat);
   }

if((strncmp(gfpar.gftype,"fk",2) == 0) || (strncmp(gfpar.gftype,"FK",2) == 0))
   {
   gfpar.flag3d = 0;
   gfpar.nc = 8;
   mstpar("gflocs","s",gfpar.gflocs);
   mstpar("gftimes","s",gfpar.gftimes);

   gfpar.swap_flag = 0;
   getpar("gf_swap_bytes","d",&gfpar.swap_flag);
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

getpar("tshift_timedomain","d",&tshift_timedomain);

endpar();

maxmech = 1;
mechpar.nmech = 1;
mechpar.flag[0] = U1FLAG;
mechpar.flag[1] = 0;
mechpar.flag[2] = 0;

if(strcmp(rupmodtype,"SRF") == 0)
   {
   maxmech = 3;

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
      np = srf.srf_apnts.np;

      apv_off = 0;
      }
   else
      {
      np = prseg_ptr[nseg].nstk*prseg_ptr[nseg].ndip;

   /* reset POINTS pointer to appropriate segment */

      apv_off = 0;
      for(i=0;i<nseg;i++)
         apv_off = apv_off + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

      apval_ptr = apval_ptr + apv_off;
      }
   }

get_gfpars(&gfpar);

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

if(statfile[0] != '\0')
   {
   stat_buf = (char *) check_malloc (NCHAR_STAT*nbuf*sizeof(char));
   slon_buf = (float *) check_malloc (nbuf*sizeof(float));
   slat_buf = (float *) check_malloc (nbuf*sizeof(float));

   fpr = fopfile(statfile,"r");

   nstat = 0;
   while(fgets(string,2048,fpr) != NULL)
      {
      if(strncmp(string,"#",1) != 0)
         {
	 stat_ptr = stat_buf + NCHAR_STAT*nstat;
	 sscanf(string,"%f %f %s",&slon_buf[nstat],&slat_buf[nstat],stat_ptr);
	 nstat++;

	 if(nstat % nbuf == 0)
	    {
            stat_buf = (char *) check_realloc (stat_buf,NCHAR_STAT*(nstat+nbuf)*sizeof(char));
            slon_buf = (float *) check_realloc (slon_buf,(nstat+nbuf)*sizeof(float));
            slat_buf = (float *) check_realloc (slat_buf,(nstat+nbuf)*sizeof(float));
	    }
	 }
      }

   stat_buf = (char *) check_realloc (stat_buf,NCHAR_STAT*nstat*sizeof(char));
   slon_buf = (float *) check_realloc (slon_buf,nstat*sizeof(float));
   slat_buf = (float *) check_realloc (slat_buf,nstat*sizeof(float));
   }
else if(latloncoords == 1)
   {
   nstat = 1;
   stat_buf = (char *) check_malloc (NCHAR_STAT*sizeof(char));
   slon_buf = (float *) check_malloc (sizeof(float));
   slat_buf = (float *) check_malloc (sizeof(float));

   strcpy(stat_buf,stat);
   slon_buf[0] = slon;
   slat_buf[0] = slat;
   }

for(istat=0;istat<nstat;istat++)
   {
   stat_ptr = stat_buf + NCHAR_STAT*istat;
   fprintf(stderr,"stat= %s\n",stat_ptr);

   sgtmast.geoproj = geoproj;
   sgtmast.modelrot = -90;
   sgtmast.xshift = -100.0;
   sgtmast.yshift = -100.0;

   sgtmast.modellon = slon_buf[istat];
   sgtmast.modellat = slat_buf[istat];

   set_geoproj(&sgtmast,&geop);

   zapit(seis,3*ntout);

   for(i=0;i<4;i++)
      {
      gfhead[i].id = -1;  /* initialize: -1 means none read yet */
      gfhead[i].ir = -1;  /* initialize: -1 means none read yet */
      gfhead[i].read_flag = 0;
      }

   tmom = 0.0;
   for(ip=0;ip<np;ip++)
      {
      zapit(subseis,maxmech*3*ntout);

      get_srfpars_v2(&srf,apv_off,ip,&rt,&vslip,&mechpar);

      if(vslip > (float)(0.0))
         {
         get_ard_srf(&srf,apv_off,ip,&azi,&rng,&z0,&deast,&dnorth,&slon_buf[istat],&slat_buf[istat],&geop);

         find_4gf_v2(gfpar,gfhead,&rng,&z0,&deast,&dnorth,&maxgft);

         if(use_closest_gf == 1 && rng > min_taper_range)
            reset_wt_4gf_v2(gfhead,&rng,&min_taper_range,&max_taper_range,&tadjust);

         read_4gf_v2(gfpath,gfname,gf,ntsum,gfhead,gfpar,&maxnt,&dtout,space);

         scale = normf*apval_ptr[ip].area;
         mech_4gf_v2(gfmech,gf,gfhead,gfpar,ntsum,mechpar,&azi,&scale);
         tmom = tmom + vslip*scale;

         sum_4gf_v2(subseis,ntout,gfmech,gfhead,ntsum,maxnt,&rt,&maxgft,&tstart,tshift_timedomain,mechpar,&tadjust);

         srf_stf(&srf,apv_off,ip,seis,subseis,stf,ntout,&dtout,mechpar,space);
         }
      }

   sv = seis;
   sn = seis + ntout;
   se = seis + 2*ntout;

   strncpy(sname,stat_ptr,7);
   sname[7] = '\0';

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

   write_seis(outdir,stat_ptr,sname,"000",sn,&dtout,ntout,&tstart);
   write_seis(outdir,stat_ptr,sname,"090",se,&dtout,ntout,&tstart);
   write_seis(outdir,stat_ptr,sname,"ver",sv,&dtout,ntout,&tstart);
   }
}
