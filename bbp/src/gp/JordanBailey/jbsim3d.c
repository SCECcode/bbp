#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

struct sgtindex *get_sgtpars(struct sgtfileparams *sgtfpar,struct sgtmaster *sgtmast,struct sgtindex *sgtindx);

int main(int ac,char **av)
{
FILE *fopfile(), *fpr;
struct sgtfileparams sgtfilepar, sgtextract;
struct sgtparams *sgtparms;
struct sgtmaster sgtmast;
struct sgtindex *sgtindx;
struct sgtindex eqindx, statindx;
struct geoprojection geop;
float *gfmech;
float *stf, *seis, *subseis, *se, *sn, *sv;
float rt, scale, slon, slat;
float elon, elat, edep;
float vslip, *space;
float z0, strike, dip, rake;
int fdw, ip, maxmech, nstf, ntsum, maxnt, ig;
float mindt, maxdelta, fweight;

struct sgtheader *sgthead;
float *sgtbuf;

char string[2048], outfile[2048];
char rupmodfile[2048], outdir[2048], sgtdir[2048], stat[64], sname[8];

struct standrupformat srf;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

struct mechparam mechpar;

long long *indx_master;
int nm, non_exact;

int extract_sgt = 0;

float tmom = 0.0;
float dtout = -1.0;
float slip_conv = 1.0;  /* input slip in cm on each subfault */
float tstart = 0.0;

float sdep = 0.0;

int ntout = -99;
int apv_off = 0;
int inbin = 0;

int intmem = 0;
float memlen;
float maxmem = 10000;  /* max RAM in Mb */

int ptol;
int print_tol = 25;

sname[0] = '\0';

sgtfilepar.xfile[0] = '\0';
sgtfilepar.yfile[0] = '\0';
sgtfilepar.zfile[0] = '\0';

sgtfilepar.xfdr = -1;
sgtfilepar.yfdr = -1;
sgtfilepar.zfdr = -1;

sgtextract.xfile[0] = '\0';
sgtextract.yfile[0] = '\0';
sgtextract.zfile[0] = '\0';

sgtextract.xfdr = -1;
sgtextract.yfdr = -1;
sgtextract.zfdr = -1;

sprintf(outdir,".");
sprintf(sgtdir,".");

setpar(ac, av);

mstpar("slat","f",&slat);
mstpar("slon","f",&slon);
getpar("outdir","s",outdir);
mstpar("stat","s",stat);

mstpar("rupmodfile","s",rupmodfile);
getpar("slip_conv","f",&slip_conv);
getpar("inbin","d",&inbin);

getpar("sgt_xfile","s",sgtfilepar.xfile);
getpar("sgt_yfile","s",sgtfilepar.yfile);
getpar("sgt_zfile","s",sgtfilepar.zfile);

if(sgtfilepar.xfile[0] == '\0' && sgtfilepar.yfile[0] == '\0' && sgtfilepar.zfile[0] == '\0')
   {
   fprintf(stderr,"*** need to specify at least one of sgt_xfile, sgt_yfile, or sgt_zfile; exiting ...\n");
   exit(-1);
   }

getpar("extract_sgt","d",&extract_sgt);
if(extract_sgt > 0)
   {
   getpar("sgtdir","s",sgtdir);

   if(sgtfilepar.xfile[0] != '\0')
      mstpar("extract_sgt_xfile","s",sgtextract.xfile);
   if(sgtfilepar.yfile[0] != '\0')
      mstpar("extract_sgt_yfile","s",sgtextract.yfile);
   if(sgtfilepar.zfile[0] != '\0')
      mstpar("extract_sgt_zfile","s",sgtextract.zfile);
   }

getpar("ntout","d",&ntout);
getpar("dtout","f",&dtout);
getpar("tstart","f",&tstart);
getpar("sname","s",sname);

getpar("maxmem","f",&maxmem);

endpar();

read_srf(&srf,rupmodfile,inbin);
prect_ptr = &srf.srf_prect;
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &srf.srf_apnts;
apval_ptr = apnts_ptr->apntvals;

apv_off = 0;

sgtindx = get_sgtpars(&sgtfilepar,&sgtmast,sgtindx);
set_geoproj(&sgtmast,&geop);

eqindx.h = sgtindx[0].h;
statindx.h = sgtindx[0].h;

get_indx(&slon,&slat,&sdep,&statindx,&geop);

/* find sgt locations (indx) for all fault points */

fprintf(stderr,"Find SGTs for this rupture\n");

sgtparms = (struct sgtparams *) check_malloc ((srf.srf_apnts.np)*sizeof(struct sgtparams));

ptol = print_tol;
maxdelta = 0.0;
fweight = 1.0;
non_exact = 0;
for(ip=0;ip<srf.srf_apnts.np;ip++)
   {
   elon = apval_ptr[ip].lon;
   elat = apval_ptr[ip].lat;
   edep = apval_ptr[ip].dep;
   get_indx(&elon,&elat,&edep,&eqindx,&geop);

   find_sgt(&sgtparms[ip],&sgtmast,sgtindx,&eqindx,&statindx,&maxdelta,&fweight);

   if(sgtparms[ip].nsgt != 1)
      non_exact++;

   if((float)(100.0*(float)(ip+1)/(float)(srf.srf_apnts.np)) >= ptol)
      {
      fprintf(stderr," %3d percent done (%d of %d)\n",ptol,ip,srf.srf_apnts.np);
      ptol = ptol + print_tol;
      }
   }

/* find unique and sorted master list of locations (indx) */

indx_master = (long long *) check_malloc (4*(srf.srf_apnts.np)*sizeof(long long));
get_master_list(sgtparms,srf.srf_apnts.np,indx_master,&nm);
indx_master = (long long *) check_realloc (indx_master,nm*sizeof(long long));

fprintf(stderr,"nm= %d non_exact= %d\n",nm,non_exact);
fprintf(stderr,"max_mindelta= %f weight= %f\n",sgtindx[0].h*sqrt(maxdelta),fweight);

if(extract_sgt > 0)
   {
   fprintf(stderr,"Extracting SGTs for this rupture\n");

   sgt_subset(&sgtfilepar,&sgtextract,&sgtmast,sgtindx,nm,indx_master,sgtdir);
   exit(0);
   }

fprintf(stderr,"Constructing synthetic for this rupture\n");

/* try to read all SGTs into memory */

memlen = (sizeof(struct sgtmaster) + nm*(sizeof(struct sgtindex) + sizeof(struct sgtheader) + 18*(sgtmast.nt)*sizeof(float)))*1.0e-06;

fprintf(stderr,"Total memory for SGTs= %.2f Mb\n",memlen);
if(memlen > maxmem)
   {
   fprintf(stderr,"*** Error: SGT memory (%f Mb) > maxmem (%f Mb), check/reset maxmem; exiting ...\n",memlen,maxmem);
   exit(-1);
   }

sgthead = (struct sgtheader *) check_malloc (nm*sizeof(struct sgtheader));
sgtbuf = (float *) check_malloc (18*nm*(sgtmast.nt)*sizeof(float));

get_sgt_subset(&sgtfilepar,&sgtmast,sgtindx,sgthead,nm,indx_master,sgtbuf);

maxnt = sgthead[0].nt;
mindt = sgthead[0].dt;

if(dtout < 0.0)
   dtout = mindt;

if(dtout < mindt)
   maxnt = (maxnt*mindt/dtout);

ntsum = 2;
while(ntsum < 4*maxnt)
   ntsum = ntsum*2;

if(ntout < 0)
   ntout = ntsum;

maxmech = 3;
mechpar.nmech = 1;
mechpar.flag[0] = U1FLAG;
mechpar.flag[1] = 0;
mechpar.flag[2] = 0;

gfmech = (float *) check_malloc (maxmech*12*ntsum*sizeof(float));
space = (float *) check_malloc (2*ntsum*sizeof(float));

seis = (float *) check_malloc (3*ntout*sizeof(float));
subseis = (float *) check_malloc (maxmech*3*ntout*sizeof(float));
stf = (float *) check_malloc (ntout*sizeof(float));

zapit(seis,3*ntout);

ptol = print_tol;
tmom = 0.0;
for(ip=0;ip<srf.srf_apnts.np;ip++)
   {
   zapit(subseis,maxmech*3*ntout);

   get_srfpars(&srf,apv_off,ip,&rt,&vslip,&mechpar.stk,&mechpar.dip,&mechpar.rak,&mechpar);
   scale = slip_conv*apval_ptr[ip].area;

   mech_sgt(gfmech,sgtbuf,sgthead,&sgtparms[ip],ntsum,mechpar,&scale);
   tmom = tmom + vslip*scale;

   sum_sgt(subseis,ntout,gfmech,&sgtparms[ip],sgthead,ntsum,&rt,&tstart,mechpar);
   srf_stf(&srf,apv_off,ip,seis,subseis,stf,ntout,&dtout,mechpar,space);

   if((float)(100.0*(float)(ip+1)/(float)(srf.srf_apnts.np)) >= ptol)
      {
      fprintf(stderr," %3d percent done (%d of %d)\n",ptol,ip,srf.srf_apnts.np);
      ptol = ptol + print_tol;
      }
   }

sv = seis;
sn = seis + ntout;
se = seis + 2*ntout;

if(sname[0] == '\0')
   {
   strncpy(sname,stat,7);
   sname[7] = '\0';
   }

makedir(outdir);
write_seis(outdir,stat,sname,"000",sn,&dtout,ntout,&tstart);
write_seis(outdir,stat,sname,"090",se,&dtout,ntout,&tstart);
write_seis(outdir,stat,sname,"ver",sv,&dtout,ntout,&tstart);

fprintf(stderr,"Total moment= %13.5e\n",tmom);
}
