#include "include.h"
#include "structure.h"
#include "function.h"

#define MAXLINE 256
#define MAXFILE 1000000
#define INCFILE 1000

main(int ac,char **av)
{
FILE *fopfile(), *fplist;
struct gfheader gfhead[4];
float maxgft;
struct gfparam gfpar;
int kg, ig, kfalloc, nfile;

float kperd_n, kperd_e;
float elon = -118.0;
float elat = 34.0;
float slat, slon, snorth, seast;
double e2, den, g2, lat0;

float rt = 0.0;
float strike = 0.0;
float dip = 90.0;
float rake = 180.0;
float dtop = 0.1;

float len, wid;
int i, j, k, l, ip, ip0;

struct standrupformat srf;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

struct mechparam mechpar;
int maxmech;

float vslip;

int apv_off = 0;
int nseg = -1;
int inbin = 0;

int kp;
float len2, ds0, dd0, dsf, ddf, s2;
int ntsum, maxnt, it, ntp2;
float mindt;
float x0, y0, z0, dd;
int nsubstk, nsubdip;
int nfinestk = 1;
int nfinedip = 1;

float scale, arg, cosD, sinD;
float azi, rng, deast, dnorth;
int ncomp = 3;

float *space;
float dtout = -1.0;

int fdw;
char gfpath[128], gfname[64];
char rupmodfile[128], filelist[256];
char string[MAXLINE], *listbuf;

double rperd = 0.017453293;
int latloncoords = 0;
int list_flag = 0;

sprintf(gfpar.gftype,"3d");

setpar(ac, av);
getpar("latloncoords","d",&latloncoords);

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

mstpar("rupmodfile","s",rupmodfile);
getpar("nseg","d",&nseg);
getpar("inbin","d",&inbin);

mstpar("filelist","s",filelist);

mstpar("gftype","s",gfpar.gftype);

if((strncmp(gfpar.gftype,"fk",2) == 0) || (strncmp(gfpar.gftype,"FK",2) == 0))
   {
   gfpar.flag3d = 0;
   gfpar.nc = 8;
   mstpar("gflocs","s",gfpar.gflocs);
   mstpar("gftimes","s",gfpar.gftimes);
   }
else if((strncmp(gfpar.gftype,"3d",2) == 0) || (strncmp(gfpar.gftype,"3D",2) == 0))
   {
   gfpar.flag3d = 1;
   gfpar.nc = 18;
   mstpar("gflocs","s",gfpar.gflocs);
   mstpar("gfrange_tolerance","f",&gfpar.rtol);
   }
else
   {
   fprintf(stderr,"gftype= %s invalid option, exiting...\n",gfpar.gftype);
   exit(-1);
   }

mstpar("gfpath","s",gfpath);
mstpar("gfname","s",gfname);

endpar();

kfalloc = INCFILE;
listbuf = (char *) check_malloc (kfalloc*MAXFILE*sizeof(char));

maxmech = 1;
mechpar.nmech = 1;
mechpar.flag[0] = U1FLAG;
mechpar.flag[1] = 0;
mechpar.flag[2] = 0;

maxmech = 3;
nfinestk = 1;
nfinedip = 1;

read_srf(&srf,rupmodfile,inbin);
prect_ptr = &srf.srf_prect;
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &srf.srf_apnts;
apval_ptr = apnts_ptr->apntvals;

if(nseg < 0)  /* do all POINTS */
   {
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

arg = dip*rperd;
cosD = cos(arg);
sinD = sin(arg);

get_gfpars(&gfpar);

if(latloncoords) /* calculate lat,lon to km conversions */
   set_ne(&elon,&elat,&slon,&slat,&snorth,&seast);

/* Calculate subfault responses */

ds0 = len/nsubstk;
dd0 = wid/nsubdip;

dsf = ds0/nfinestk;
ddf = dd0/nfinedip;

for(i=0;i<4;i++)
   {
   gfhead[i].id = -1;  /* initialize: -1 means none read yet */
   gfhead[i].ir = -1;  /* initialize: -1 means none read yet */
   }

nfile = 0;
for(i=0;i<nsubstk;i++)
   {
   for(j=0;j<nsubdip;j++)
      {
      ip0 = i + j*nsubstk;

      fprintf(stderr,"ip0= %8d/%d nfile= %8d\n",ip0,nsubdip*nsubstk,nfile);

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

            get_srfpars(&srf,apv_off,ip0,&rt,&vslip,&strike,&dip,&rake,&mechpar);
	    get_ard_srf(&srf,apv_off,ip0,&azi,&rng,&z0,&deast,&dnorth,&slon,&slat);
	    find_4gf(gfpar,gfhead,&rng,&z0,&deast,&dnorth);

	    for(ig=0;ig<4;ig++)
	       {
               getname_gf(string,gfname,&gfhead[ig],gfpar);

	       list_flag = check_name(string,listbuf,nfile,MAXLINE);

               if(list_flag == 0)  /* new file */
	          {
		  if(nfile+1 > kfalloc)
		     {
		     kfalloc = kfalloc + INCFILE;
                     listbuf = (char *) check_realloc (listbuf,kfalloc*MAXLINE*sizeof(char));
		     }

                  strcpy(listbuf+nfile*MAXLINE,string);
		  nfile++;
		  }
	       }
	    }
         }
      }
   }

fplist = fopfile(filelist,"w");
for(ig=0;ig<nfile;ig++)
   fprintf(fplist,"%s/%s\n",gfpath,listbuf+ig*MAXLINE);
fclose(fplist);
}
