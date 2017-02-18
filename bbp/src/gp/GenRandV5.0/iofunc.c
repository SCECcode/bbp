#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void write_field(char *file,struct pointsource *ps,char *type,int nx,int ny,float *dx,float *dy)
{
FILE *fpw, *fopfile();
float xx, yy;
int i, j;

fpw = fopfile(file,"w");
for(j=0;j<ny;j++)
   {
   yy = (j + 0.5)*(*dy);
   for(i=0;i<nx;i++)
      {
      xx = (i+0.5)*(*dx);

      if(strcmp(type,"slip") == 0)
         fprintf(fpw,"%12.5e %12.5e %12.5e\n",xx,yy,ps[i+j*nx].slip);
      else if(strcmp(type,"rupt") == 0)
         fprintf(fpw,"%12.5e %12.5e %12.5e\n",xx,yy,ps[i+j*nx].rupt);
      }
   }
fclose(fpw);
}

void write_spec(char *file,float *as,struct complex *slip,int nx,int ny,float *dx,float *dy,float *dkx,float *dky,float *xl,float *yl,int kflag)
{
FILE *fpw;
float kx, ky, amp, amp0, lamp;
float xl2, yl2, fac;
int i, j, ip;

float hcoef = 1.8;  /* H=0.8, hcoef = H + 1 */

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

fft2d(slip,nx,ny,-1,dx,dy);

fpw = fopfile(file,"w");

/* this will normalize max spectrum to one */
amp0 = sqrt(slip[0].re*slip[0].re + slip[0].im*slip[0].im);

for(j=0;j<=ny/2;j++)
   {
   if(j<=ny/2)
      ky = j*(*dky);
   else
      ky = (j-ny)*(*dky);

   for(i=0;i<=nx/2;i++)
      {
      if(i<=nx/2)
         kx = i*(*dkx);
      else
         kx = (i-nx)*(*dkx);

      ip = i + j*nx;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = 1.0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = 1.0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = 1.0/sqrt(1.0 + amp*amp);

      amp = sqrt(slip[ip].re*slip[ip].re + slip[ip].im*slip[ip].im);
      amp = amp/amp0;

      lamp = -1.e+20;
      if(amp > 0.0)
	 lamp = log10(amp/fac);

      fprintf(fpw,"%13.5e %13.5e %12.5e %12.5e\n",kx,ky,lamp,180*atan(slip[ip].im/slip[ip].re)/3.14159);

      as[ip] = as[ip] + lamp;
      }
   }
fclose(fpw);

fft2d(slip,nx,ny,1,dkx,dky);
}

void write_avgspec(char *file,float *as,int ns,int nx,int ny,float *dkx,float *dky)
{
FILE *fpw;
float kx, ky, fac;
int i, j, ip;

fac = 1.0/(float)(ns);

fpw = fopfile(file,"w");

for(j=0;j<=ny/2;j++)
   {
   ky = j*(*dky);

   for(i=0;i<=nx/2;i++)
      {
      kx = i*(*dkx);
      ip = i + j*nx;
      fprintf(fpw,"%13.5e %13.5e %12.5e\n",kx,ky,fac*as[ip]);
      }
   }
fclose(fpw);
}

FILE *fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

int opfile_ro(char *name)
{
int fd;
if ((fd = open (name, RDONLY_FLAGS, 0444)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int opfile(char *name)
{
int fd;
if ((fd = open (name, RDWR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int croptrfile(char *name)
{
int fd;
if ((fd = open (name, CROPTR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int reed(int fd, void *pntr, int length)
{
int temp;
if ((temp = read(fd, pntr, length)) < length)
   {
   fprintf (stderr, "READ ERROR\n");
   fprintf (stderr, "%d attempted  %d read\n", length, temp);
   exit(-1);
   }
return(temp);
}

int rite(int fd, void *pntr, int length)
{
int temp;
if ((temp = write(fd, pntr, length)) < length)
   {
   fprintf (stderr, "WRITE ERROR\n");
   fprintf (stderr, "%d attempted  %d written\n", length, temp);
   exit(-1);
   }
return(temp);
}

struct pointsource *read_ruppars(char *file,struct pointsource *psrc,float *mag,int *nx,int *ny,float *dx,float *dy,float *dtop,float *stk,float *dip,float *elon,float *elat)
{
FILE *fpr, *fopfile();
float area;
int i, nn;
char str[1024];

double rperd = 0.017453293;

*dtop = 1.0e+15;
*stk = 0.0;
*dip = 0.0;
*elon = 0.0;
*elat = 0.0;

if(strcmp(file,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(file,"r");

fgets(str,1024,fpr);   /* Probability = <float> */

fgets(str,1024,fpr);   /* Magnitude = <float> */
sscanf(str,"%*s %*s %f",mag);

fgets(str,1024,fpr);   /* GridSpacing = <float> */
sscanf(str,"%*s %*s %f",dx);

if(*dx == (float)(0.0))
   {
   fprintf(stderr,"***** input error\n");
   fprintf(stderr,"      GridSpacing = 0.0, exiting...\n");
   exit(-1);
   }

*dy = *dx;

fgets(str,1024,fpr);   /* NumRows = <int> */
sscanf(str,"%*s %*s %d",ny);

fgets(str,1024,fpr);   /* NumCols = <int> */
sscanf(str,"%*s %*s %d",nx);

fgets(str,1024,fpr);   /* header comment */

psrc = (struct pointsource *)check_realloc(psrc,(*nx)*(*ny)*sizeof(struct pointsource));

area = (*dx)*(*dy)*1.0e+10;  /* km -> cm */

for(i=0;i<(*nx)*(*ny);i++)
   {
   fgets(str,1024,fpr);   /* Lat , Lon , Depth , Rake , Dip , Strike */
   sscanf(str,"%f %f %f %f %f %f",&psrc[i].lat,
                                  &psrc[i].lon,
                                  &psrc[i].dep,
                                  &psrc[i].rak,
                                  &psrc[i].dip,
                                  &psrc[i].stk);

   psrc[i].area = area;

   if(psrc[i].dep < *dtop)
      *dtop = psrc[i].dep;

   *stk = *stk + psrc[i].stk;
   *dip = *dip + psrc[i].dip;
   }
fclose(fpr);

*stk = *stk/((*nx)*(*ny));
*dip = *dip/((*nx)*(*ny));

nn = 0;
for(i=0;i<(*nx)*(*ny);i++)
   {
   if(psrc[i].dep < (*dtop + 0.01))
      {
      *elon = *elon + psrc[i].lon;
      *elat = *elat + psrc[i].lat;
      nn++;
      }
   }

/* adjust for half subfault width */
*dtop = (*dtop) - 0.5*(*dx)*sin((*dip)*rperd);

if(nn == 0)
   nn++;

*elon = *elon/nn;
*elat = *elat/nn;

return(psrc);
}

struct pointsource *set_ruppars(struct pointsource *psrc,float *mag,int *nx,int *ny,float *dx,float *dy,float *dtop,float *stk,float *dip,float *rak,float *elon,float *elat)
{
float area;
float cosA, sinA, cosD, sinD, fwid, flen;
float xx, yy, zz, dd, sn, se;
int i, j, ip;

double rperd = 0.017453293;

psrc = (struct pointsource *)check_realloc(psrc,(*nx)*(*ny)*sizeof(struct pointsource));

area = (*dx)*(*dy)*1.0e+10;  /* km -> cm */
flen = (*nx)*(*dx);
fwid = (*ny)*(*dy);

cosA = cos((*stk)*rperd);
sinA = sin((*stk)*rperd);
cosD = cos((*dip)*rperd);
sinD = sin((*dip)*rperd);
for(j=0;j<(*ny);j++)
   {
   dd = (j + 0.5)*(*dy);
   yy = dd*cosD;
   zz = (*dtop) + dd*sinD;

   for(i=0;i<(*nx);i++)
      {
      ip = i + j*(*nx);
      xx = (i+0.5)*(*dx) - 0.5*flen;

      se = xx*sinA + yy*cosA;
      sn = xx*cosA - yy*sinA;
      set_ll(elon,elat,&psrc[ip].lon,&psrc[ip].lat,&sn,&se);

      psrc[ip].dep = zz;
      psrc[ip].stk = (*stk);
      psrc[ip].dip = (*dip);
      psrc[ip].rak = (*rak);
      psrc[ip].area = area;
      }
   }

return(psrc);
}

struct pointsource *read_gsfpars(char *file,struct pointsource *psrc,struct generic_slip *gslip,float *dx,float *dy,float *dtop,float *dip)
{
FILE *fpr, *fopfile();
int i, nn;
char str[1024];
struct slippars *spar;

double dbldx = 0.0;
double dbldy = 0.0;

double rperd = 0.017453293;

*dtop = 1.0e+15;
*dx = 0.0;
*dy = 0.0;
*dip = 0.0;

if(strcmp(file,"stdin") == 0)
   fpr = stdin;
else
   fpr = fopfile(file,"r");

fgets(str,1024,fpr);
while(strncmp(str,"#",1) == 0)
   fgets(str,1024,fpr);

sscanf(str,"%d",&gslip->np);

psrc = (struct pointsource *)check_realloc(psrc,gslip->np*sizeof(struct pointsource));
gslip->spar = (struct slippars *)check_realloc(gslip->spar,gslip->np*sizeof(struct slippars));
spar = gslip->spar;

i = 0;
while(fgets(str,1024,fpr) != NULL)
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

   psrc[i].lon = spar[i].lon;
   psrc[i].lat = spar[i].lat;
   psrc[i].dep = spar[i].dep;
   psrc[i].stk = spar[i].stk;
   psrc[i].dip = spar[i].dip;
   psrc[i].rak = spar[i].rake;
   psrc[i].area = spar[i].ds*spar[i].dw*1.0e+10;

   if(psrc[i].dep < *dtop)
      *dtop = psrc[i].dep;

   *dip = *dip + psrc[i].dip;
/*
   *dx = *dx + spar[i].ds;
   *dy = *dy + spar[i].dw;
*/

   dbldx = dbldx + spar[i].ds;
   dbldy = dbldy + spar[i].dw;

   i++;
   }
fclose(fpr);

*dip = *dip/(gslip->np);
/*
*dx = *dx/(gslip->np);
*dy = *dy/(gslip->np);
*/

*dx = dbldx/(gslip->np);
*dy = dbldy/(gslip->np);

/* adjust for half subfault width */
*dtop = (*dtop) - 0.5*(*dy)*sin((*dip)*rperd);
if(*dtop < 0.0)
   *dtop = 0.0;

return(psrc);
}
