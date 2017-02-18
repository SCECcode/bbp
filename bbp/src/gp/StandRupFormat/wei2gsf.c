#include "include.h"
#include "structure.h"
#include "function.h"
#include "../../ModelCords/function.h"
#include "defs.h"

#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292

main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float flon, flat, fstk, fdip, ftop, flen, fwid, dlen, dwid, l2;
float *slon, *slat, *sstk, *sdip, *srak, *sslip, *sdep, *sdx, *sdy, *stinit, *strise, *stratio;
float shypo, dhypo, trise_min, trise_interval;
float seg_dip, seg_stk, rup_delay;
int seg_num, seg_nx, seg_ny; 
float *slip0, *rake0, *twin0, *trat0, *vrup0;
float mlat, mlon, mrot;
float kperd_n, kperd_e, ar1, ar2, ar3, ar4;
float sp, wp, xp, yp, zp, plon, plat;
float cosR, sinR, cosD, sinD;
float ds, dd;
float dd0, ds0, sp0, wp0;
int nlen0, nwid0, *segno, nsublen, nsubwid;
int nstk, ndip, i, j, k, mp, ip;

float trise_shift = 0.0;

int input_rupture_time = 0;

char *rchar, infile[512], outfile[512], str[1024];

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;

double g0, b0;
double amat[9], ainv[9];
int xy2ll = 0;
int ll2xy = 1;

double xr0, yr0, dlon, dlat, dxr, dyr;
int geo2utm = 0;
int utm2geo = 1;
int utm_zone = 11;

/*
   geoproj=0: RWG spherical projection with local kmplat, kmplon
          =1: RWG great circle projection
	  =2: UTM coordinate projection
*/
int geoproj = 1;

float xshift;
float yshift;

float xlen = 1000.0;
float ylen = 1000.0;

float slip_scale = 1.0;
float rake_range = 0.0;
long seed = 1;

int get_tau1_ratio = 0;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac, av);

getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("geoproj","d",&geoproj);

mstpar("flon","f",&flon);
mstpar("flat","f",&flat);
mstpar("flen","f",&flen);
mstpar("fwid","f",&fwid);
mstpar("ftop","f",&ftop);
mstpar("fstk","f",&fstk);
mstpar("fdip","f",&fdip);
mstpar("dlen","f",&dlen);
mstpar("nsublen","d",&nsublen);
mstpar("dwid","f",&dwid);
mstpar("nsubwid","d",&nsubwid);
mstpar("shypo","f",&shypo);
mstpar("dhypo","f",&dhypo);

mstpar("trise_interval","f",&trise_interval);
if(trise_interval > 0.0)
   {
   mstpar("trise_min","f",&trise_min);
   getpar("trise_shift","f",&trise_shift);
   }

getpar("get_tau1_ratio","d",&get_tau1_ratio);

getpar("input_rupture_time","d",&input_rupture_time);

endpar();

nlen0 = (int)(flen/dlen + 0.5);
nwid0 = (int)(fwid/dwid + 0.5);

l2 = 0.5*flen;

ds0 = flen/(float)(nlen0);
dd0 = fwid/(float)(nwid0);

slip0 = (float *)check_malloc(nlen0*nwid0*sizeof(float));
rake0 = (float *)check_malloc(nlen0*nwid0*sizeof(float));
twin0 = (float *)check_malloc(nlen0*nwid0*sizeof(float));
trat0 = (float *)check_malloc(nlen0*nwid0*sizeof(float));
vrup0 = (float *)check_malloc(nlen0*nwid0*sizeof(float));

if(strcmp(infile,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(infile,"r");

fgets(str,1024,fpr);
while(strncmp(str,"#",1) == 0)
   {
   rchar = fgets(str,1024,fpr);
   if(rchar == NULL)
      {
      fprintf(stderr,"Unexpected EOF in %s, exiting...\n",infile);
      exit(-99);
      }
   }

sscanf(str,"%d %f %f",&seg_num,&seg_dip,&seg_stk);
fgets(str,1024,fpr);
sscanf(str,"%d %d %f",&seg_nx,&seg_ny,&rup_delay);

for(k=0;k<nwid0;k++)
   {
   wp0 = (k+0.5)*dd0;
   for(j=0;j<nlen0;j++)
      {
      i = j + k*nlen0;
      sp0 = (j+0.5)*ds0 - l2;

      fgets(str,1024,fpr);

      trat0[i] = -1.0;
      if(get_tau1_ratio == 1)
         sscanf(str,"%f %f %f %f %f",&slip0[i],&rake0[i],&twin0[i],&trat0[i],&vrup0[i]);
      else
         sscanf(str,"%f %f %f %f",&slip0[i],&rake0[i],&twin0[i],&vrup0[i]);

      if(trise_interval > 0.0)   /* else twin0 is actual rise_rime in seconds 09/17/2011 */
         twin0[i] = trise_min + trise_interval*(twin0[i] - 1.0);

      if(input_rupture_time != 0) /* rupture time is input, convert to rupture velocity */
	 {
	 if(vrup0[i] > 0.0)
            vrup0[i] = sqrt((sp0-shypo)*(sp0-shypo)+(wp0-dhypo)*(wp0-dhypo))/vrup0[i];
	 else
            vrup0[i] = 2.8;
	 }

      if(twin0[i] < 0.0)
         twin0[i] = 0.0;
      }
   }
fclose(fpr);

nstk = nsublen*nlen0;
ndip = nsubwid*nwid0;

slon = (float *)check_malloc(nstk*ndip*sizeof(float));
slat = (float *)check_malloc(nstk*ndip*sizeof(float));
sdep = (float *)check_malloc(nstk*ndip*sizeof(float));
sdx = (float *)check_malloc(nstk*ndip*sizeof(float));
sdy = (float *)check_malloc(nstk*ndip*sizeof(float));
sstk = (float *)check_malloc(nstk*ndip*sizeof(float));
sdip = (float *)check_malloc(nstk*ndip*sizeof(float));
srak = (float *)check_malloc(nstk*ndip*sizeof(float));
sslip = (float *)check_malloc(nstk*ndip*sizeof(float));
stinit = (float *)check_malloc(nstk*ndip*sizeof(float));
segno = (int *)check_malloc(nstk*ndip*sizeof(int));
strise = (float *)check_malloc(nstk*ndip*sizeof(float));
stratio = (float *)check_malloc(nstk*ndip*sizeof(float));

mlon = flon;
mlat = flat;
mrot = fstk - 90.0;

cosR = cos(mrot*rperd);
sinR = sin(mrot*rperd);

xshift = -0.5*xlen;
yshift = -0.5*ylen;

if(geoproj == 0)
   {
   radc = ERAD*RPERD;
   set_g2(&g2,&fc);

   latavg = mlat;
   geocen(&latavg,(double)(latavg*rperd));
   latlon2km(&latavg,&kperd_n,&kperd_e,&radc,&g2);

   ar1 = cosR/kperd_e;
   ar2 = sinR/kperd_e;
   ar3 = cosR/kperd_n;
   ar4 = sinR/kperd_n;
   }
else if(geoproj == 1)
   {
   gen_matrices(amat,ainv,&mrot,&mlon,&mlat);

   g0 = (double)(0.5*ylen)/(double)(erad);
   b0 = (double)(0.5*xlen)/(double)(erad);
   }
else if(geoproj == 2)
   {
   dlon = mlon;
   dlat = mlat;
   geoutm_(&dlon,&dlat,&xr0,&yr0,&utm_zone,&geo2utm);
   }

ds = flen/(float)(nstk);
dd = fwid/(float)(ndip);

cosD = cos(fdip*rperd);
sinD = sin(fdip*rperd);

for(k=0;k<ndip;k++)
   {
   wp = (k+0.5)*dd;
   yp = wp*cosD - yshift;
   zp = wp*sinD + ftop;
   for(j=0;j<nstk;j++)
      {
      ip = j + k*nstk;
      mp = j/nsublen + (k/nsubwid)*nlen0;

      sp = (j+0.5)*ds - l2;
      xp = sp - xshift;

      if(geoproj == 0)
         {
         plon = mlon + (xp+xshift)*ar1 - (yp+yshift)*ar2;
         plat = mlat - (xp+xshift)*ar4 - (yp+yshift)*ar3;
         }
      else if(geoproj == 1)
         {
         gcproj(&xp,&yp,&plon,&plat,&erad,&g0,&b0,amat,ainv,xy2ll);
         }
      else if(geoproj == 2)
         {
         dxr = xr0 + 1000.0*((xp+xshift)*ar1 - (yp+yshift)*ar2);
         dyr = yr0 - 1000.0*((xp+xshift)*ar4 + (yp+yshift)*ar3);
         geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);

         plon = dlon;
         plat = dlat;
         }

      slon[ip] = plon;
      slat[ip] = plat;
      sdep[ip] = zp;
      sdx[ip] = ds;
      sdy[ip] = dd;
      sstk[ip] = fstk;
      sdip[ip] = fdip;
      srak[ip] = rake0[mp];
      sslip[ip] = slip0[mp];
      stinit[ip] = sqrt((sp-shypo)*(sp-shypo)+(wp-dhypo)*(wp-dhypo))/vrup0[mp] + rup_delay;
      segno[ip] = 0;
      strise[ip] = twin0[mp];
      stratio[ip] = trat0[mp];
      }
   }

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

fprintf(fpw,"# nstk= %d ndip= %d\n",nstk,ndip);
fprintf(fpw,"# flen= %10.4f fwid= %10.4f\n",flen,fwid);
fprintf(fpw,"# LON  LAT  DEP(km)  SUB_DX  SUB_DY  LOC_STK  LOC_DIP  LOC_RAKE  SLIP(cm)  INIT_TIME  SEG_NO  RISE_TIME TAU1_RATIO\n");

fprintf(fpw,"%d\n",nstk*ndip);

for(ip=0;ip<nstk*ndip;ip++)
   fprintf(fpw,"%11.5f %11.5f %8.4f %8.4f %8.4f %6.1f %6.1f %6.1f %8.2f %8.3f %3d %8.3f %10.5f\n",slon[ip],
                                                                                     slat[ip],
                                                                                     sdep[ip],
                                                                                     sdx[ip],
                                                                                     sdy[ip],
                                                                                     sstk[ip],
                                                                                     sdip[ip],
                                                                                     srak[ip],
                                                                                     sslip[ip],
                                                                                     stinit[ip],
                                                                                     segno[ip],
                                                                                     strise[ip],
                                                                                     stratio[ip]);

fclose(fpw);
}
