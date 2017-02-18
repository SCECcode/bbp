#include "include.h"

#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)
#define         KM2FT           (float)(3280.84)
#define MAX_NLEN 100

struct fdcf
   {
   float c4[4];
   float c2[2];
   };

struct gridparam
   {
   struct fdcf *cfx0;
   struct fdcf *cfx1;
   struct fdcf *cfy0;
   struct fdcf *cfy1;
   struct fdcf *cfz0;
   struct fdcf *cfz1;
   float *hx;
   float *hy;
   float *hz;
   float *xp;
   float *yp;
   float *zp;
   float hmin;
   float hmax;
   float dt;
   float dthx0;
   float dthxn;
   float dthy0;
   float dthyn;
   float dthz0;
   float dthzn;
   int order;
   int nx;
   int ny;
   int nz;
   };

#include "function.h"

void *check_malloc(size_t);
void *check_realloc(void *,size_t);
void set_gridparams(char *,struct gridparam *, char *,int);
void setcoef(struct fdcf *,struct fdcf *,float *,int);
void gelim_double(double *,int,double *);
void getlens(FILE *,char *,float *,float *,float *,float *,int *);

main(int ac,char **av)
{
struct gridparam gp;
char gridfile[256], gridout[256];
FILE *fpw, *fpr, *fopfile();
float *flon, *flat, modellon, modellat, modelrot, x, y;
float ar1, ar2, ar3, ar4, latavg;
float cosR, sinR, kmlon, kmlat;
float xlen, ylen, zlen;
float xmax, ymax, zmax;
float dep, *cc1, *cc2, conv;
int ix, iy, iz, ip;
char outfile[128], str[512], name[64];

double rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc;
double geocenX();

double g0, b0;
double amat[9], ainv[9];
int xy2ll = 0;
int ll2xy = 1;

int do_coords = 0;
int nzout = 1;

int latfirst = 1;
int feetout = 1;
int gzip = 1;

double xr0, yr0, dlon, dlat, dxr, dyr;
int geo2utm = 0;
int utm2geo = 1;
int utm_zone = 11;

/*
   geoproj=0: RWG spherical projection with local kmplat, kmplon
          =1: RWG great circle projection
          =2: UTM coordinate projection
*/
int geoproj = 0;  /* default is OLD way */
int center_origin = 0;

float xshift = -1.0e+15;
float yshift = -1.0e+15;

sprintf(gridout,"gridout");

setpar(ac,av);
mstpar("gridfile","s",gridfile);
getpar("gridout","s",gridout);
mstpar("modellat","f",&modellat);
mstpar("modellon","f",&modellon);
mstpar("modelrot","f",&modelrot);
getpar("do_coords","d",&do_coords);
if(do_coords)
   {
   mstpar("name","s",name);
   getpar("nzout","d",&nzout);
   getpar("latfirst","d",&latfirst);
   getpar("feetout","d",&feetout);
   getpar("gzip","d",&gzip);
   }

getpar("geoproj","d",&geoproj);
getpar("center_origin","d",&center_origin);
getpar("xshift","f",&xshift);
getpar("yshift","f",&yshift);
endpar();

set_gridparams(gridfile,&gp,gridout,0);  /* free surface = 0 */

xmax = gp.xp[gp.nx-1];
ymax = gp.yp[gp.ny-1];
zmax = gp.zp[gp.nz-1];

xlen = gp.xp[gp.nx-1] + 0.5*(gp.hx[gp.nx-1] + gp.hx[0]);
ylen = gp.yp[gp.ny-1] + 0.5*(gp.hy[gp.ny-1] + gp.hy[0]);
zlen = gp.zp[gp.nz-1] + 0.5*(gp.hz[gp.nz-1] + gp.hz[0]);

flon = (float *) check_malloc (4*sizeof(float));
flat = (float *) check_malloc (4*sizeof(float));

cosR = cos(modelrot*RPERD);
sinR = sin(modelrot*RPERD);

if(center_origin != 0)
   {
   xshift = -0.5*xmax;
   yshift = -0.5*ymax;
   }
else
   {
   if(xshift < -1.0e+14)
      xshift = 0.0;
   if(yshift < -1.0e+14)
      yshift = 0.0;
   }

if(geoproj == 0)
   {
   radc = ERAD*RPERD;
   set_g2(&g2,&fc);

   latavg = modellat;
   if(center_origin == 0)  /* backward compatible */
      latavg = modellat - 0.5*(xlen*sinR + ylen*cosR)/111.20;

   geocen(&latavg,(double)(latavg)*rperd);
   latlon2km(&latavg,&kmlat,&kmlon,&radc,&g2);

   flon[0] = modellon + ((xshift)*cosR - (yshift)*sinR)/kmlon;
   flat[0] = modellat - ((xshift)*sinR + (yshift)*cosR)/kmlat;
   flon[1] = modellon + ((xmax+xshift)*cosR - (yshift)*sinR)/kmlon;
   flat[1] = modellat - ((xmax+xshift)*sinR + (yshift)*cosR)/kmlat;
   flon[2] = modellon + ((xmax+xshift)*cosR - (ymax+yshift)*sinR)/kmlon;
   flat[2] = modellat - ((xmax+xshift)*sinR + (ymax+yshift)*cosR)/kmlat;
   flon[3] = modellon + ((xshift)*cosR - (ymax+yshift)*sinR)/kmlon;
   flat[3] = modellat - ((xshift)*sinR + (ymax+yshift)*cosR)/kmlat;
   }
else if(geoproj == 1)
   {
   gen_matrices(amat,ainv,&modelrot,&modellon,&modellat);

   g0 = (double)(0.5*ymax)/(double)(erad);
   b0 = (double)(0.5*xmax)/(double)(erad);

   x = 0.0;
   y = 0.0;
   gcproj(&x,&y,&flon[0],&flat[0],&erad,&g0,&b0,amat,ainv,xy2ll);

   x = xmax;
   y = 0.0;
   gcproj(&x,&y,&flon[1],&flat[1],&erad,&g0,&b0,amat,ainv,xy2ll);

   x = xmax;
   y = ymax;
   gcproj(&x,&y,&flon[2],&flat[2],&erad,&g0,&b0,amat,ainv,xy2ll);

   x = 0.0;
   y = ymax;
   gcproj(&x,&y,&flon[3],&flat[3],&erad,&g0,&b0,amat,ainv,xy2ll);
   }
else if(geoproj == 2)
   {
   dlon = modellon;
   dlat = modellat;
   geoutm_(&dlon,&dlat,&xr0,&yr0,&utm_zone,&geo2utm);
   
   dxr = xr0 + 1000.0*((xshift)*cosR - (yshift)*sinR);
   dyr = yr0 - 1000.0*((xshift)*sinR + (yshift)*cosR);
   geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);
   flon[0] = dlon;
   flat[0] = dlat;

   dxr = xr0 + 1000.0*((xmax+xshift)*cosR - (yshift)*sinR);
   dyr = yr0 - 1000.0*((xmax+xshift)*sinR + (yshift)*cosR);
   geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);
   flon[1] = dlon;
   flat[1] = dlat;

   dxr = xr0 + 1000.0*((xmax+xshift)*cosR - (ymax+yshift)*sinR);
   dyr = yr0 - 1000.0*((xmax+xshift)*sinR + (ymax+yshift)*cosR);
   geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);
   flon[2] = dlon;
   flat[2] = dlat;

   dxr = xr0 + 1000.0*((xshift)*cosR - (ymax+yshift)*sinR);
   dyr = yr0 - 1000.0*((xshift)*sinR + (ymax+yshift)*cosR);
   geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);
   flon[3] = dlon;
   flat[3] = dlat;

   kmlon = 1.0;
   kmlat = 1.0;
   }

printf("Model origin coordinates:\n");
printf(" lon= %10.5f lat= %10.5f rotate= %7.2f\n\n",modellon,modellat,modelrot);

printf("Model origin shift (cartesian vs. geographic):\n");
printf(" xshift(km)= %12.5f yshift(km)= %12.5f\n\n",xshift,yshift);

printf("Model corners:\n");
printf(" c1= %10.5f %10.5f\n",flon[0],flat[0]);
printf(" c2= %10.5f %10.5f\n",flon[1],flat[1]);
printf(" c3= %10.5f %10.5f\n",flon[2],flat[2]);
printf(" c4= %10.5f %10.5f\n\n",flon[3],flat[3]);

printf("Model Dimensions:\n");
printf(" xlen= %10.4f km\n",xlen);
printf(" ylen= %10.4f km\n",ylen);
printf(" zlen= %10.4f km\n\n",zlen);

if(geoproj == 0)
   {
   printf("Unit Conversions:\n");
   printf(" km/lon= %10.4f\n",kmlon);
   printf(" km/lat= %10.4f\n",kmlat);
   printf(" latavg= %10.4f\n",latavg/rperd);
   }
else if(geoproj == 2)
   printf("UTM Zone= %d\n",utm_zone);

if(do_coords)
   {
   ar1 = cosR/kmlon;
   ar2 = sinR/kmlon;
   ar3 = cosR/kmlat;
   ar4 = sinR/kmlat;

   flon = (float *) check_realloc (flon,gp.nx*gp.ny*sizeof(float));
   flat = (float *) check_realloc (flat,gp.nx*gp.ny*sizeof(float));

   for(iy=0;iy<gp.ny;iy++)
      {
      fprintf(stderr,"iy=%d\n",iy);

      if(geoproj == 0)
         {
         for(ix=0;ix<gp.nx;ix++)
            {
            ip = ix + iy*gp.nx;

            flon[ip] = modellon + (gp.xp[ix]+xshift)*ar1 - (gp.yp[iy]+yshift)*ar2;
            flat[ip] = modellat - (gp.xp[ix]+xshift)*ar4 - (gp.yp[iy]+yshift)*ar3;
            }
         }
      else if(geoproj == 1)
         {
         for(ix=0;ix<gp.nx;ix++)
            {
            ip = ix + iy*gp.nx;

            x = gp.xp[ix];
            y = gp.yp[iy];
            gcproj(&x,&y,&flon[ip],&flat[ip],&erad,&g0,&b0,amat,ainv,xy2ll);
            }
         }
      else if(geoproj == 2)
         {
         for(ix=0;ix<gp.nx;ix++)
            {
            ip = ix + iy*gp.nx;

            dxr = xr0 + 1000.0*((gp.xp[ix]+xshift)*ar1 - (gp.yp[iy]+yshift)*ar2);
            dyr = yr0 - 1000.0*((gp.xp[ix]+xshift)*ar4 + (gp.yp[iy]+yshift)*ar3);
            geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&utm2geo);

            flon[ip] = dlon;
            flat[ip] = dlat;
            }
         }
      }

   conv = 1.0;
   if(feetout)
      {
      conv = KM2FT;
      fprintf(stderr,"**** depth output in FEET\n");
      }
   else
      fprintf(stderr,"**** depth output in KM\n");

   cc1 = flon;
   cc2 = flat;
   if(latfirst)
      {
      cc1 = flat;
      cc2 = flon;
      fprintf(stderr,"**** coordinate output format is LAT LON DEP\n");
      }
   else
      fprintf(stderr,"**** coordinate output format is LON LAT DEP\n");

   if(gzip)
      fprintf(stderr,"**** output files are compressed with 'gzip' (.gz)\n");

   for(iz=0;iz<nzout;iz++)
      {
      fprintf(stderr,"iz=%d\n",iz);

      if(iz == 0)
         sprintf(outfile,"%s",name);
      else
         sprintf(outfile,"%s%.4d",name,iz);

      fpw = fopfile(outfile,"w");

      dep = (gp.zp[iz])*conv;
      for(iy=0;iy<gp.ny;iy++)
         {
         for(ix=0;ix<gp.nx;ix++)
            {
            ip = ix + iy*gp.nx;

            fprintf(fpw,"%10.4f %10.4f %5d %5d\n",cc1[ip],cc2[ip],ix,iy);
            }
         }
      fclose(fpw);

      if(gzip)
         {
         sprintf(str,"gzip -f %s",outfile);
         system(str);
         }
      }
   }
}

FILE *fopfile(name,mode)
char *name, *mode;
{
FILE *fp, *fopen();

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

double geocenX(x)
double x;
{
double r;
r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
return(r);
}

set_g2X(g2,fc)
float *g2, *fc;
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

latlon2kmX(arg,latkm,lonkm,rc,g2)
float *arg, *latkm, *lonkm, *rc, *g2;
{
float cosA, sinA, g2s2, den;

cosA = cos((*arg));
sinA = sin((*arg));
g2s2 = (*g2)*sinA*sinA;

den = sqrt((FONE)/((FONE) + g2s2));
*lonkm = (*rc)*cosA*den;
*latkm = (*rc)*(sqrt((FONE) + g2s2*((FTWO) + (*g2))))*den*den*den;
}

void *check_realloc(void *ptr,size_t len)
{
ptr = (char *) realloc (ptr,len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   exit(-1);
   }

return(ptr);
}

void *check_malloc(size_t len)
{
char *ptr;

ptr = (char *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }

return(ptr);
}

void set_gridparams(char *gfile,struct gridparam *gp, char *gout,int fs)
{
float x0[MAX_NLEN], x1[MAX_NLEN], dxlen[MAX_NLEN];
float y0[MAX_NLEN], y1[MAX_NLEN], dylen[MAX_NLEN];
float z0[MAX_NLEN], z1[MAX_NLEN], dzlen[MAX_NLEN];
float xlen, ylen, zlen;
float xx, yy, zz;
double dd;
int nxlen, nylen, nzlen;
int nx, ny, nz;
int i, j;

FILE *fpr, *fpw, *fopfile();
char str[512];

float hmin = 1.0e+15;
float hmax = -1.0e+15;
float szero = 0.0;

fpr = fopfile(gfile,"r");

while(fgets(str,512,fpr) != NULL)
   {
   if(strncmp(str,"xlen=",5) == 0)
      {
      sscanf(&str[5],"%f",&xlen);
      getlens(fpr,str,&xlen,x0,x1,dxlen,&nxlen);
      }
   else if(strncmp(str,"ylen=",5) == 0)
      {
      sscanf(&str[5],"%f",&ylen);
      getlens(fpr,str,&ylen,y0,y1,dylen,&nylen);
      }
   else if(strncmp(str,"zlen=",5) == 0)
      {
      sscanf(&str[5],"%f",&zlen);
      getlens(fpr,str,&zlen,z0,z1,dzlen,&nzlen);
      }
   }
fclose(fpr);

nx = 0;
for(i=0;i<nxlen;i++)
   nx = nx + (int)((x1[i] - x0[i])/dxlen[i] + 0.5);

gp->hx = (float *) check_malloc (nx*sizeof(float));
gp->xp = (float *) check_malloc (nx*sizeof(float));
gp->cfx0 = (struct fdcf *) check_malloc (nx*sizeof(struct fdcf));
gp->cfx1 = (struct fdcf *) check_malloc (nx*sizeof(struct fdcf));
gp->nx = nx;

j = 0;
dd = x0[0];
for(i=0;i<nx;i++)
   {
   gp->hx[i] = dxlen[j];
   dd = dd + gp->hx[i];
   gp->xp[i] = dd - 0.5*gp->hx[i];

   if(dxlen[j] < hmin)
      hmin = dxlen[j];
   if(dxlen[j] > hmax)
      hmax = dxlen[j];

   if(dd >= x1[j] && j < nxlen-1)
      j++;
   }

xx = gp->xp[0];     /* reset origin */

for(i=0;i<nx;i++)
   gp->xp[i] = gp->xp[i] - xx;

setcoef(gp->cfx0,gp->cfx1,gp->hx,gp->nx);

ny = 0;
for(i=0;i<nylen;i++)
   ny = ny + (int)((y1[i] - y0[i])/dylen[i] + 0.5);

gp->hy = (float *) check_malloc (ny*sizeof(float));
gp->yp = (float *) check_malloc (ny*sizeof(float));
gp->cfy0 = (struct fdcf *) check_malloc (ny*sizeof(struct fdcf));
gp->cfy1 = (struct fdcf *) check_malloc (ny*sizeof(struct fdcf));
gp->ny = ny;

j = 0;
dd = y0[0];
for(i=0;i<ny;i++)
   {
   gp->hy[i] = dylen[j];
   dd = dd + gp->hy[i];
   gp->yp[i] = dd - 0.5*gp->hy[i];

   if(dylen[j] < hmin)
      hmin = dylen[j];
   if(dylen[j] > hmax)
      hmax = dylen[j];

   if(dd >= y1[j] && j < nylen-1)
      j++;
   }

yy = gp->yp[0];     /* reset origin */

for(i=0;i<ny;i++)
   gp->yp[i] = gp->yp[i] - yy;

setcoef(gp->cfy0,gp->cfy1,gp->hy,gp->ny);

nz = 0;
for(i=0;i<nzlen;i++)
   nz = nz + (int)((z1[i] - z0[i])/dzlen[i] + 0.5);

if(fs)
   {
   z0[0] = -dzlen[0];
   nz++;
   }

gp->hz = (float *) check_malloc (nz*sizeof(float));
gp->zp = (float *) check_malloc (nz*sizeof(float));
gp->cfz0 = (struct fdcf *) check_malloc (nz*sizeof(struct fdcf));
gp->cfz1 = (struct fdcf *) check_malloc (nz*sizeof(struct fdcf));
gp->nz = nz;

j = 0;
dd = z0[0];
for(i=0;i<nz;i++)
   {
   gp->hz[i] = dzlen[j];
   dd = dd + gp->hz[i];
   gp->zp[i] = dd - 0.5*gp->hz[i];

   if(dzlen[j] < hmin)
      hmin = dzlen[j];
   if(dzlen[j] > hmax)
      hmax = dzlen[j];

   if(dd >= z1[j] && j < nzlen-1)
      j++;
   }

zz = gp->zp[0];     /* reset origin */
if(fs)
   zz = gp->zp[1];

for(i=0;i<nz;i++)
   gp->zp[i] = gp->zp[i] - zz;

setcoef(gp->cfz0,gp->cfz1,gp->hz,gp->nz);

gp->hmin = hmin;
gp->hmax = hmax;

if(gout[0] != '\0')
   {
   fpw = fopfile(gout,"w");

   fprintf(fpw,"xlen=%.5f\n",xlen);
   fprintf(fpw,"nx=%d\n",nx);
   for(i=0;i<nx;i++)
      {
      fprintf(fpw,"%6d %13.5e %13.5e\n",i,gp->xp[i],gp->hx[i]);
      if(gp->hx[i] <= (float)(0.0))
         {
         fprintf(stderr,"**** Problem with zero hx, check '%s'\n",gout);
         fprintf(stderr,"     exiting...\n");
         exit(-1);
         }
      }

   fprintf(fpw,"ylen=%.5f\n",ylen);
   fprintf(fpw,"ny=%d\n",ny);
   for(i=0;i<ny;i++)
      {
      fprintf(fpw,"%6d %13.5e %13.5e\n",i,gp->yp[i],gp->hy[i]);
      if(gp->hy[i] <= szero)
         {
         fprintf(stderr,"**** Problem with zero hy, check '%s'\n",gout);
         fprintf(stderr,"     exiting...\n");
         exit(-1);
         }
      }

   fprintf(fpw,"zlen=%.5f\n",zlen);
   fprintf(fpw,"nz=%d\n",nz);
   for(i=0;i<nz;i++)
      {
      fprintf(fpw,"%6d %13.5e %13.5e\n",i,gp->zp[i],gp->hz[i]);
      if(gp->hz[i] <= (float)(0.0))
         {
         fprintf(stderr,"**** Problem with zero hz, check '%s'\n",gout);
         fprintf(stderr,"     exiting...\n");
         exit(-1);
         }
      }

   fclose(fpw);
   }
}

void setcoef(struct fdcf *cg0,struct fdcf *cg1,float *hg,int ng)
{
float dd[4];
double am[16], bv[4];
int i, ig;
int sgn[4];

for(ig=0;ig<ng;ig++)
   {
   for(i=0;i<4;i++)
      {
      cg0[ig].c4[i] = 0.0;
      cg1[ig].c4[i] = 0.0;
      }

   for(i=0;i<2;i++)
      {
      cg0[ig].c2[i] = 0.0;
      cg1[ig].c2[i] = 0.0;
      }
   }

sgn[0] = 1;
sgn[1] = 1;
sgn[2] = -1;
sgn[3] = -1;

for(ig=1;ig<ng-1;ig++)
   {
   dd[0] = hg[ig+1] + 0.5*hg[ig];
   dd[1] = 0.5*hg[ig];
   dd[2] = 0.5*hg[ig];
   dd[3] = hg[ig-1] + 0.5*hg[ig];

   for(i=0;i<4;i++)
      {
      am[i]    = 1.0;
      am[i+4]  = sgn[i]*dd[i]*am[i];
      am[i+8]  = sgn[i]*dd[i]*am[i+4];
      am[i+12] = sgn[i]*dd[i]*am[i+8];
      }

   bv[0] = 0.0;
   bv[1] = 1.0;
   bv[2] = 0.0;
   bv[3] = 0.0;

   gelim_double(am,4,bv);

   cg0[ig].c4[0] = bv[0];
   cg0[ig].c4[1] = bv[1];
   cg0[ig].c4[2] = bv[2];
   cg0[ig].c4[3] = bv[3];
   }

for(ig=1;ig<ng-2;ig++)
   {
   dd[0] = hg[ig+1] + 0.5*hg[ig+2];
   dd[1] = 0.5*hg[ig+1];
   dd[2] = 0.5*hg[ig];
   dd[3] = hg[ig] + 0.5*hg[ig-1];

   for(i=0;i<4;i++)
      {
      am[i] = 1.0;
      am[i+4] = sgn[i]*dd[i]*am[i];
      am[i+8] = sgn[i]*dd[i]*am[i+4];
      am[i+12] = sgn[i]*dd[i]*am[i+8];
      }

   bv[0] = 0.0;
   bv[1] = 1.0;
   bv[2] = 0.0;
   bv[3] = 0.0;

   gelim_double(am,4,bv);

   cg1[ig].c4[0] = bv[0];
   cg1[ig].c4[1] = bv[1];
   cg1[ig].c4[2] = bv[2];
   cg1[ig].c4[3] = bv[3];
   }

for(ig=0;ig<ng-1;ig++)
   {
   cg0[ig].c2[0] = 1.0/hg[ig];
   cg0[ig].c2[1] = -cg0[ig].c2[0];

   cg1[ig].c2[0] = 2.0/(hg[ig+1] + hg[ig]);
   cg1[ig].c2[1] = -cg1[ig].c2[0];
   }
}

void getlens(FILE *fp,char *s,float *len,float *g0,float *g1,float *dg,int *n)
{
int i = 0;

fgets(s,512,fp);
sscanf(s,"%f %f %f",&g0[0],&g1[0],&dg[0]);
g0[0] = 0.0;  /* force to be origin */

while(g1[i] < (*len))
   {
   i++;
   fgets(s,512,fp);
   sscanf(s,"%f %f %f",&g0[i],&g1[i],&dg[i]);
   }

*n = i+1;
g1[i] = (*len);  /* force to be total length */
}

/*
Gauss elimination without pivoting (no row exchanges):
We solve Ax = b where A is n by n, and x and b have length n
by forming the decomposition A = LU via elimination.
Originally a contains the matrix A and b contains the
vector b.  At the end a contains the lower and upper
triangular matrices L and U and b contains the solution
vector x.  The diagonal of a contains the diagonal of U
(the pivots) since the diagonal elements of U are all 1's.
*/

void gelim_double(double *a,int n,double *b)
   {
        int i, j, jj;
        double *pa, *paj;
        double f, pivot;

         for (j=0; j<n-1; j++)   /* lu decomp of a (no row exchanges) */
            {
                pa= a + j*n;
                pivot = pa[j];
                for (i=j+1; i<n; i++)
                   {
                        pa = a + i*n;
                        paj= a + j*n;
                        f = pa[j]/pivot;
                        pa[j] = f;
                        for (jj=j+1; jj<n; jj++) pa[jj] -= f*paj[jj];
                   }
           }
         for (i=1; i<n; i++)        /* forward elimination on b */
            {
                pa = a + i*n;
                for (j=0; j<i; j++) b[i] -= pa[j]*b[j];
            }
         for (i=n-1; i>-1; i--)        /* back-substitution */
            {
                pa = a + i*n;
                for (j=n-1; j>i; j--) b[i] -= pa[j]*b[j];
                b[i] = b[i]/pa[i];
            }
    }
