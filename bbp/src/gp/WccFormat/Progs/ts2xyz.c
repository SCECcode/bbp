#include        "include.h"
#include        "structure.h"
#include        "function.h"

#define MAXL 512
#define         MAX_STAT 10000
#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)

struct xyz
   {
   float x;
   float y;
   float z;
   };

main(int ac,char **av)
{
FILE *fpr, *fopfile();
struct xyz *c0, *c1, *c2;
struct tsheader tshead;
float *val, *val2, *xs, *ys, *zs, *p1, *p2;
float xlen, ylen, zlen, xp, yp, vec;
int *itlist;
int dxts, dyts, dzts;
int fdr, i, i1, i2, ip, ic;
int ix, iy, iz, nx, ny, nz;
int ntts, ts, nxts, nyts, nzts, n1, n2;
char gridfile[512], str[512], infile[512], outfile[512];

float scale = 1.0;
int ncomp = 3;

int getpeak = 0;
int vectorpeak = 0;
int read_header = 0;
int swap_bytes = 0;

float c0mx, c0my, c1mx, c1my, c2mx, c2my;
int c0mt, c1mt, c2mt;
float c0max = 0.0;
float c1max = 0.0;
float c2max = 0.0;

int lonlat = 1;
float mlon, mlat, mrot;
float kperd_n, kperd_e, ar1, ar2, ar3, ar4;
float cosR, sinR;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;
double geocen();

int xyts = 0;
int xzts = 0;
int yzts = 0;

int outbin = 1;
int trv = 0;
float elon, elat, az, x1, y1;
float pi = 3.14159;

int fault_cords = 0;
float cosP, sinP;
float strike = 0.0;
float dip = 90;

float xsh = 0.0;
float ysh = 0.0;
float zsh = 0.0;

setpar(ac, av);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
mstpar("gridfile","s",gridfile);
mstpar("dxts","d",&dxts);
mstpar("dyts","d",&dyts);
mstpar("dzts","d",&dzts);

getpar("read_header","d",&read_header);

getpar("getpeak","d",&getpeak);
getpar("vectorpeak","d",&vectorpeak);
if(getpeak || vectorpeak)
   mstpar("ntts","d",&ntts);
else
   mstpar("ts","d",&ts);

getpar("xyts","d",&xyts);
getpar("xzts","d",&xzts);
getpar("yzts","d",&yzts);

getpar("xsh","f",&xsh);
getpar("ysh","f",&ysh);
getpar("zsh","f",&zsh);

getpar("ncomp","d",&ncomp);
getpar("scale","f",&scale);

getpar("outbin","d",&outbin);
getpar("swap_bytes","d",&swap_bytes);

getpar("lonlat","d",&lonlat);
if(lonlat && read_header == 0)
   {
   mstpar("mlon","f",&mlon);
   mstpar("mlat","f",&mlat);
   mstpar("mrot","f",&mrot);
   }

getpar("trv","d",&trv);
if(trv)
   {
   mstpar("elon","f",&elon);
   mstpar("elat","f",&elat);
   }

getpar("fault_cords","d",&fault_cords);
if(fault_cords)
   {
   mstpar("strike","f",&strike);
   mstpar("dip","f",&dip);
   }

endpar();

fdr = opfile_ro(infile);

if(read_header)
   {
   reed(fdr,&tshead,sizeof(struct tsheader));

   if(swap_bytes)
      {
      swap_in_place(1,(char *)(&tshead.ix0));
      swap_in_place(1,(char *)(&tshead.iy0));
      swap_in_place(1,(char *)(&tshead.iz0));
      swap_in_place(1,(char *)(&tshead.it0));
      swap_in_place(1,(char *)(&tshead.nx));
      swap_in_place(1,(char *)(&tshead.ny));
      swap_in_place(1,(char *)(&tshead.nz));
      swap_in_place(1,(char *)(&tshead.nt));
      swap_in_place(1,(char *)(&tshead.dx));
      swap_in_place(1,(char *)(&tshead.dy));
      swap_in_place(1,(char *)(&tshead.dz));
      swap_in_place(1,(char *)(&tshead.dt));
      swap_in_place(1,(char *)(&tshead.modelrot));
      swap_in_place(1,(char *)(&tshead.modellat));
      swap_in_place(1,(char *)(&tshead.modellon));
      }

   mrot = tshead.modelrot;
   mlon = tshead.modellon;
   mlat = tshead.modellat;
   }

fpr = fopfile(gridfile,"r");

fgets(str,MAXL,fpr);
getf(str,&xlen);

fgets(str,MAXL,fpr);
getd(str,&nx);

nxts = (int)(((nx-1)/(dxts))+1);
xs = (float *) check_malloc (nxts*sizeof(float));

i = 0;
for(ix=0;ix<nx;ix++)
   {
   fgets(str,MAXL,fpr);

   if(ix%dxts == 0 && (xyts || xzts))
      {
      sscanf(str,"%*d %f %*f",&xs[i]);
      xs[i] = xs[i] + xsh;
      i++;
      }
   }

if(nxts == tshead.nx-1)
   {
   nxts = tshead.nx;
   xs = (float *) check_realloc (xs,nxts*sizeof(float));
   xs[nxts-1] = 2.0*xs[nxts-2] - xs[nxts-3];
   }

fgets(str,MAXL,fpr);
getf(str,&ylen);

fgets(str,MAXL,fpr);
getd(str,&ny);

nyts = (int)(((ny-1)/(dyts))+1);
ys = (float *) check_malloc (nyts*sizeof(float));

i = 0;
for(iy=0;iy<ny;iy++)
   {
   fgets(str,MAXL,fpr);

   if(iy%dyts == 0 && (xyts || yzts))
      {
      sscanf(str,"%*d %f %*f",&ys[i]);
      ys[i] = ys[i] + ysh;
      i++;
      }
   }

if(nyts == tshead.ny-1)
   {
   nyts = tshead.ny;
   ys = (float *) check_realloc (ys,nyts*sizeof(float));
   ys[nyts-1] = 2.0*ys[nyts-2] - ys[nyts-3];
   }

fgets(str,MAXL,fpr);
getf(str,&zlen);

fgets(str,MAXL,fpr);
getd(str,&nz);

nzts = (int)(((nz-1)/(dzts))+1);
zs = (float *) check_malloc (nzts*sizeof(float));

i = 0;
for(iz=0;iz<nz;iz++)
   {
   fgets(str,MAXL,fpr);

   if(iz%dzts == 0 && (xzts || yzts))
      {
      sscanf(str,"%*d %f %*f",&zs[i]);
      zs[i] = zs[i] + zsh;
      i++;
      }
   }

if(nzts == tshead.nz-1)
   {
   nzts = tshead.nz;
   zs = (float *) check_realloc (zs,nzts*sizeof(float));
   zs[nzts-1] = 2.0*zs[nzts-2] - zs[nzts-3];
   }

fclose(fpr);

fprintf(stderr,"nxts= %d nyts= %d nzts= %d\n",nxts,nyts,nzts);
fprintf(stderr,"hdnx= %d hdny= %d hdnz= %d\n",tshead.nx,tshead.ny,tshead.nz);

if(lonlat)
   {
   cosR = cos(mrot*rperd);
   sinR = sin(mrot*rperd);

   radc = ERAD*RPERD;
   set_g2(&g2,&fc);

   latavg = mlat - 0.5*(xlen*sinR + ylen*cosR)/111.20;
   latavg = geocen(latavg*rperd);
   latlon2km(&latavg,&kperd_n,&kperd_e,&radc,&g2);

   ar1 = cosR/kperd_e;
   ar2 = sinR/kperd_e;
   ar3 = cosR/kperd_n;
   ar4 = sinR/kperd_n;
   }
else
   {
   mlon = 0.0;
   mlat = 0.0;

   ar1 = 1.0;
   ar2 = 0.0;
   ar3 = -1.0;
   ar4 = 0.0;
   }

if(xyts)
   {
   n1 = nxts;
   p1 = xs;
   n2 = nyts;
   p2 = ys;
   }

if(xzts)
   {
   n1 = nxts;
   p1 = xs;
   n2 = nzts;
   p2 = zs;
   }

if(yzts)
   {
   n1 = nzts;
   p1 = zs;
   n2 = nyts;
   p2 = ys;
   }

fprintf(stderr,"n1= %d n2= %d\n",n1,n2);

val = (float *) check_malloc (3*n1*n2*sizeof(float));
itlist = (int *) check_malloc (3*n1*n2*sizeof(int));
c0 = (struct xyz *) check_malloc (n1*n2*sizeof(struct xyz));
c1 = (struct xyz *) check_malloc (n1*n2*sizeof(struct xyz));
c2 = (struct xyz *) check_malloc (n1*n2*sizeof(struct xyz));

if(getpeak || vectorpeak)
   {
   val2 = (float *) check_malloc (3*n1*n2*sizeof(float));

   for(i=0;i<3*n1*n2;i++)
      val[i] = -1.0;

   for(i=0;i<ntts;i++)
      {
      fprintf(stderr,"%5d of %5d\n",i+1,ntts);

      reed(fdr,val2,ncomp*n1*n2*sizeof(float));
      if(swap_bytes)
         swap_in_place(ncomp*n1*n2,(char *)(val2));

      for(i2=0;i2<n2;i2++)
         {
         for(i1=0;i1<n1;i1++)
            {
            ip = i1 + i2*n1;

	    if(vectorpeak)
	       {
	       vec = (val2[ip]*val2[ip] + val2[ip+n1*n2]*val2[ip+n1*n2] + val2[ip+2*n1*n2]*val2[ip+2*n1*n2]);

	       if(vec > val[ip])
                  {
		  val[ip] = vec;
		  itlist[ip] = i;
		  }
	       }
	    else
	       {
               for(ic=0;ic<3;ic++)
                  {
	          if(val2[ip + ic*n1*n2] > val[ip + ic*n1*n2])
		     {
	             val[ip + ic*n1*n2] = val2[ip + ic*n1*n2];
		     itlist[ip+ic*n1*n2] = i;
		     }

	          else if(-val2[ip + ic*n1*n2] > val[ip + ic*n1*n2])
		     {
	             val[ip + ic*n1*n2] = -val2[ip + ic*n1*n2];
		     itlist[ip+ic*n1*n2] = i;
		     }
	          }
	       }
            }
         }
      }

   free(val2);
   }
else
   {
   for(i=0;i<3*n1*n2;i++)
      itlist[i] = ts;

   lseek(fdr,ncomp*ts*n1*n2*sizeof(float),SEEK_CUR);
   reed(fdr,val,ncomp*n1*n2*sizeof(float));
   if(swap_bytes)
      swap_in_place(ncomp*n1*n2,(char *)(val));
   }

close(fdr);

for(i2=0;i2<n2;i2++)
   {
   for(i1=0;i1<n1;i1++)
      {
      ip = i1 + i2*n1;

      if(vectorpeak)
	 val[ip] = sqrt(val[ip]);

      xp = mlon + p1[i1]*ar1 - p2[i2]*ar2;
      yp = mlat - p1[i1]*ar4 - p2[i2]*ar3;

      cosR = 1.0;
      sinR = 0.0;
      cosP = 1.0;
      sinP = 0.0;

      if(trv)
	 {
	 x1 = (elon - xp)*kperd_e;
	 y1 = (elat - yp)*kperd_n;

	 if(y1 == 0.0)
	    az = 0.5*pi;
	 else
            az = atan(x1/y1);

	 if(y1 < 0.0)
	    az = az + pi;

	 az = az - pi*mrot/180.0;

	 cosR = cos(az);
	 sinR = sin(az);
	 }

      if(fault_cords)
	 {
	 az = pi*(90 + mrot - strike)/180.0;

	 cosR = cos(az);
	 sinR = sin(az);

	 az = pi*(90 - dip)/180.0;

	 cosP = cos(az);
	 sinP = sin(az);
	 }

      c0[ip].x = xp;
      c0[ip].y = yp;
      c1[ip].x = xp;
      c1[ip].y = yp;
      c2[ip].x = xp;
      c2[ip].y = yp;

      c0[ip].z = 0.0;
      c1[ip].z = 0.0;
      c2[ip].z = 0.0;
      if(xp > -1.0e+14 && yp > -1.0e+14)
	 {
	 if(fault_cords == 1)
	    {
            c0[ip].z = (val[ip]*cosR - val[ip + n1*n2]*sinR)*scale;
            c1[ip].z = (val[ip]*sinR + val[ip + n1*n2]*cosR)*cosP*scale
	                - val[ip + 2*n1*n2]*sinP*scale;
            c2[ip].z = (val[ip]*sinR + val[ip + n1*n2]*cosR)*sinP*scale
	                + val[ip + 2*n1*n2]*cosP*scale;
	    }
	 else
	    {
            c0[ip].z = (val[ip]*cosR + val[ip + n1*n2]*sinR)*scale;
            c1[ip].z = -(val[ip]*sinR - val[ip + n1*n2]*cosR)*scale;
            c2[ip].z = val[ip + 2*n1*n2]*scale;
	    }
	 }

      if(c0[ip].z > c0max)
	 {
	 c0max = c0[ip].z;
	 c0mx = xp;
	 c0my = yp;
	 c0mt = itlist[ip];
	 }
      if(-c0[ip].z > c0max)
	 {
	 c0max = -c0[ip].z;
	 c0mx = xp;
	 c0my = yp;
	 c0mt = itlist[ip];
	 }

      if(c1[ip].z > c1max)
	 {
	 c1max = c1[ip].z;
	 c1mx = xp;
	 c1my = yp;
	 c1mt = itlist[ip+n1*n2];
	 }
      if(-c1[ip].z > c1max)
	 {
	 c1max = -c1[ip].z;
	 c1mx = xp;
	 c1my = yp;
	 c1mt = itlist[ip+n1*n2];
	 }

      if(c2[ip].z > c2max)
	 {
	 c2max = c2[ip].z;
	 c2mx = xp;
	 c2my = yp;
	 c2mt = itlist[ip+2*n1*n2];
	 }
      if(-c2[ip].z > c2max)
	 {
	 c2max = -c2[ip].z;
	 c2mx = xp;
	 c2my = yp;
	 c2mt = itlist[ip+2*n1*n2];
	 }
      }
   }

fprintf(stderr,"c0max=%13.5e xp=%10.5f yp=%10.5f it=%d\n",c0max,c0mx,c0my,c0mt);
fprintf(stderr,"c1max=%13.5e xp=%10.5f yp=%10.5f it=%d\n",c1max,c1mx,c1my,c1mt);
fprintf(stderr,"c2max=%13.5e xp=%10.5f yp=%10.5f it=%d\n",c2max,c2mx,c2my,c2mt);

write_out(outfile,0,c0,n1*n2,outbin);
write_out(outfile,1,c1,n1*n2,outbin);
write_out(outfile,2,c2,n1*n2,outbin);
}

getf(s,f)
char *s;
float *f;
{
while(strncmp(s,"=",1) != 0 && s[0] != '\0')
   s++;

if(s[0] == '\0')
   *f = -1;
else
   *f = atof(s+1);
}

getd(s,d)
char *s;
int *d;
{
while(strncmp(s,"=",1) != 0 && s[0] != '\0')
   s++;

if(s[0] == '\0')
   *d = -1;
else
   *d = atoi(s+1);
}

double geocen(x)
double x;
{
double r;
r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
return(r);
}

set_g2(g2,fc)
float *g2, *fc;
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

latlon2km(arg,latkm,lonkm,rc,g2)
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

write_out(outf,code,cc,nn,ob)
char *outf;
struct xyz *cc;
int ob, nn, code;
{
FILE *fpw, *fopfile();
int fdw, i;
char str[512];

sprintf(str,"%s.%d",outf,code);

if(ob)
   {
   fdw = croptrfile(str);
   rite(fdw,cc,nn*sizeof(struct xyz));
   close(fdw);
   }
else
   {
   fpw = fopfile(str,"w");
   for(i=0;i<nn;i++)
      fprintf(fpw,"%10.4f %10.4f %10.4f\n",cc[i].x,cc[i].y,cc[i].z);
   fclose(fpw);
   }
}

void swap_in_place(int n,char *cbuf)
{
char cv;

while(n--)
   {
   cv = cbuf[0];
   cbuf[0] = cbuf[3];
   cbuf[3] = cv;

   cv = cbuf[1];
   cbuf[1] = cbuf[2];
   cbuf[2] = cv;

   cbuf = cbuf + 4;
   }
}
