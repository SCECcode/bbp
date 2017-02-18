#include "include.h"
#include "structure.h"
#include "function.h"
#include "../ModelCords/function.h"
#include "defs.h"
#include "getpar.h"

#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292

extern void geoutm_(double *, double *, double *, double *, int *, int *);

int main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float *seg_lon, *seg_lat, *seg_stk, *seg_dip, *seg_rak, *seg_top, *seg_len, *seg_wid, l2;
float *flon, *flat, *fstk, *fdip, *frak, *fslip, *fdep, *fdx, *fdy, flen, fwid;
float **seg_slip, *seg_sptr;
int nseg;
float mlat, mlon, mrot;
float kperd_n, kperd_e, ar1, ar2, ar3, ar4;
float wp, xp, yp, zp, plon, plat;
float cosR, sinR, cosD, sinD;
float ds, dd;
int *seg_nstk, *seg_ndip, *segno, off_nstk;
int nstk, ndip, i, j, k, mp, ip;

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

int read_slip_vals = 1;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac, av);

getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("geoproj","d",&geoproj);
getpar("slip_scale","f",&slip_scale);
getpar("rake_range","f",&rake_range);

getpar("read_slip_vals","d",&read_slip_vals);

endpar();

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

sscanf(str,"%d",&nseg);

seg_lon = (float *)check_malloc(nseg*sizeof(float));
seg_lat = (float *)check_malloc(nseg*sizeof(float));
seg_top = (float *)check_malloc(nseg*sizeof(float));
seg_stk = (float *)check_malloc(nseg*sizeof(float));
seg_dip = (float *)check_malloc(nseg*sizeof(float));
seg_rak = (float *)check_malloc(nseg*sizeof(float));
seg_len = (float *)check_malloc(nseg*sizeof(float));
seg_wid = (float *)check_malloc(nseg*sizeof(float));
seg_nstk = (int *)check_malloc(nseg*sizeof(int));
seg_ndip = (int *)check_malloc(nseg*sizeof(int));

seg_slip = (float **)check_malloc(nseg*sizeof(float *));

nstk = 0;
for(i=0;i<nseg;i++)
   {
   fgets(str,1024,fpr);
   sscanf(str,"%f %f %f %f %f %f %f %f %d %d",&seg_lon[i],&seg_lat[i],&seg_top[i],&seg_stk[i],&seg_dip[i],&seg_rak[i],&seg_len[i],&seg_wid[i],&seg_nstk[i],&seg_ndip[i]);

   nstk = nstk + seg_nstk[i];

   if(i == 0)
      ndip = seg_ndip[i];
   else if(ndip != seg_ndip[i])
      {
      fprintf(stderr,"problem with unequal seg_ndip, exiting...\n");
      exit(-99);
      }

   seg_slip[i] = (float *)check_malloc(seg_nstk[i]*seg_ndip[i]*sizeof(float *));
   seg_sptr = seg_slip[i];

   if(read_slip_vals)
      {
      for(k=0;k<seg_ndip[i];k++)
         {
         for(j=0;j<seg_nstk[i];j++)
            {
	    ip = j + k*seg_nstk[i];
	    fscanf(fpr,"%f",&seg_sptr[ip]);
            }
         }
      fgets(str,1024,fpr);  /* get rogue newline character */
      }
   else
      {
      for(k=0;k<seg_ndip[i];k++)
         {
         for(j=0;j<seg_nstk[i];j++)
            {
	    ip = j + k*seg_nstk[i];
	    seg_sptr[ip] = -1.0;
            }
         }
      }
   }
fclose(fpr);

flon = (float *)check_malloc(nstk*ndip*sizeof(float));
flat = (float *)check_malloc(nstk*ndip*sizeof(float));
fdep = (float *)check_malloc(nstk*ndip*sizeof(float));
fdx = (float *)check_malloc(nstk*ndip*sizeof(float));
fdy = (float *)check_malloc(nstk*ndip*sizeof(float));
fstk = (float *)check_malloc(nstk*ndip*sizeof(float));
fdip = (float *)check_malloc(nstk*ndip*sizeof(float));
frak = (float *)check_malloc(nstk*ndip*sizeof(float));
fslip = (float *)check_malloc(nstk*ndip*sizeof(float));
segno = (int *)check_malloc(nstk*ndip*sizeof(int));

xlen = 0.0;
ylen = 0.0;

off_nstk = 0;
flen = 0.0;
for(i=0;i<nseg;i++)
   {
   mlon = seg_lon[i];
   mlat = seg_lat[i];
   mrot = seg_stk[i] - 90.0;

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

   seg_sptr = seg_slip[i];

   ds = seg_len[i]/(float)(seg_nstk[i]);
   dd = seg_wid[i]/(float)(seg_ndip[i]);

   l2 = 0.5*seg_len[i];

   cosD = cos(seg_dip[i]*rperd);
   sinD = sin(seg_dip[i]*rperd);

   for(k=0;k<seg_ndip[i];k++)
      {
      wp = (k+0.5)*dd;
      yp = wp*cosD - yshift;
      zp = wp*sinD + seg_top[i];
      for(j=0;j<seg_nstk[i];j++)
         {
	 ip = j + off_nstk + k*nstk;
	 mp = j + k*seg_nstk[i];

         xp = (j+0.5)*ds - l2 - xshift;

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

	 flon[ip] = plon;
	 flat[ip] = plat;
	 fdep[ip] = zp;
	 fdx[ip] = ds;
	 fdy[ip] = dd;
	 fstk[ip] = seg_stk[i];
	 fdip[ip] = seg_dip[i];
	 /*
	 */
	 frak[ip] = seg_rak[i];
	 frak[ip] = seg_rak[i] + rake_range*sfrand(&seed);
	 fslip[ip] = seg_sptr[mp];
	 segno[ip] = i;
         }
      }
   off_nstk = off_nstk + seg_nstk[i];
   flen = flen + seg_len[i];
   fwid = seg_wid[i];
   }

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

fprintf(fpw,"# nstk= %d ndip= %d\n",nstk,ndip);
fprintf(fpw,"# flen= %10.4f fwid= %10.4f\n",flen,fwid);
fprintf(fpw,"# LON  LAT  DEP(km)  SUB_DX  SUB_DY  LOC_STK  LOC_DIP  LOC_RAKE  SLIP(cm)  INIT_TIME  SEG_NO\n");

fprintf(fpw,"%d\n",nstk*ndip);

for(ip=0;ip<nstk*ndip;ip++)
   fprintf(fpw,"%11.5f %11.5f %11.5e %11.5e %11.5e %6.1f %6.1f %6.1f %8.2f %8.3f %3d\n",flon[ip],
                                                                                     flat[ip],
                                                                                     fdep[ip],
                                                                                     fdx[ip],
                                                                                     fdy[ip],
                                                                                     fstk[ip],
                                                                                     fdip[ip],
                                                                                     frak[ip],
                                                                                     slip_scale*fslip[ip],
                                                                                     -1.0,
                                                                                     segno[ip]);

fclose(fpw);
}
